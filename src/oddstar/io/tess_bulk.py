"""TESS light curves for every star of a live sample (survey mode of the optical-dip search).

Vetting fetches a QLP light curve and an independent TESScut extraction for a star that some channel has already
flagged (`optical_timeseries.fetch_tess_pipeline_pair`). That is two minutes or more per star, and it means the
optical-dip channel never sees a star that no other channel flagged. This module fetches only what the dip search
itself reads, for every star: the QLP `SAP_FLUX` light curve of each sector, from the public MAST bucket.

What is stored. A QLP file is 0.2 to 1.3 MB and a 7,800-star sample has tens of thousands of them, more than the raw
cache holds. Each file is therefore reduced on arrival to what the search uses, good-quality cadences averaged into
30-minute bins, and that reduced series is what the cache keeps (gzip JSON, a few kB). The record carries the URL,
the size and the SHA-256 of the file it was made from, so the reduction can be checked against the archive. A sector
without a QLP file (HTTP 403/404 from the bucket) is stored as "absent", so a re-run does not ask again.

Statuses of one star, as `tess_fetch.status` and `vetting_archive.query_status["tess_qlp"]`:
- "ok": at least one sector gave a light curve;
- "no_coverage": no TIC entry, no TESS sector covers the position, or no sector has a QLP light curve;
- "failed: ...": a lookup or every existing sector failed after the configured attempts (not cached, retried next run);
- "not_fetched": beyond `--tess-max-stars`.
Only "ok" stars are searched; the rest are reported as not searched and excluded from the null result.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.io.physics import propagate_position
from oddstar.schemas import ManifestRecord

BINNED_RELEASE = "qlp_binned_v1"
# Fewer good cadences than this in a sector is not a light curve
MIN_GOOD_CADENCES = 50

ProgressCallback = Callable[[int, int, str, str], None]


def _bulk_settings(cfg: PipelineConfig) -> tuple[int, int, float, int, float]:
    ch = cfg.channels["optical_dips"]
    return (int(ch["tess_bulk_max_sectors"]), int(ch["tess_bulk_max_attempts"]), float(ch["tess_bulk_backoff_s"]),
            int(ch["tess_bulk_workers"]), float(ch["tess_bulk_rate_limit_rps"]))


def reduce_qlp_fits(fits_bytes: bytes) -> dict[str, Any]:
    """Good-quality `SAP_FLUX` cadences of one QLP file averaged into 30-minute bins (flux in file units, not normalised)."""
    t, f = ot.read_qlp_raw(fits_bytes)
    if len(t) < MIN_GOOD_CADENCES:
        return {"n_good": int(len(t)), "cadence_s": None, "time": [], "flux": [], "n_per_bin": []}
    order = np.argsort(t)
    t, f = t[order], f[order]
    cadence_s = float(np.median(np.diff(t))) * 86400.0
    idx = np.floor((t - t[0]) / ot.TESS_BIN_DAYS).astype(int)
    _, inv = np.unique(idx, return_inverse=True)
    cnt = np.bincount(inv)
    tb = np.bincount(inv, weights=t) / cnt
    fb = np.bincount(inv, weights=f) / cnt
    return {"n_good": int(len(t)), "cadence_s": round(cadence_s, 1), "time": [round(float(v), 5) for v in tb],
            "flux": [float(f"{float(v):.7g}") for v in fb], "n_per_bin": [int(v) for v in cnt]}


def _pack(record: dict[str, Any]) -> bytes:
    raw = json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8")
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:  # mtime=0: identical input gives identical bytes
        gz.write(raw)
    return buf.getvalue()


def _unpack(payload: bytes) -> dict[str, Any]:
    return json.loads(gzip.decompress(payload).decode("utf-8"))


def fetch_qlp_sector_binned(
    tic_id: str, sector: int, cache: RawCache, *, sleep: ot.Sleeper = time.sleep,
) -> tuple[dict[str, Any], ManifestRecord, bool]:
    """The reduced light curve of one (TIC, sector), or its "absent" record. Returns (record, manifest, from_cache).

    Raises ArchiveAccessError when the bucket cannot be reached and ValueError when a file that arrived is not a
    readable QLP light curve; neither is stored.
    """
    url = ot.qlp_url(tic_id, sector)
    sig = hashlib.sha256(f"GET {url} -> binned".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("tess", BINNED_RELEASE, sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return _unpack(cache.read_bytes(cached)), cached, True
        except (FileNotFoundError, OSError, ValueError):
            pass  # evicted or unreadable blob: fetch again
    resp = ot._get_with_retries(url, None, timeout_s=60.0, attempts=3, sleep=sleep, ok_statuses=(200, 403, 404))
    if resp.status_code != 200:
        record: dict[str, Any] = {"status": "absent", "http_status": int(resp.status_code), "url": url, "sector": int(sector), "tic_id": str(tic_id)}
    else:
        raw = resp.content
        try:
            reduced = reduce_qlp_fits(raw)
        except (KeyError, OSError, ValueError, TypeError, IndexError) as exc:
            raise ValueError(f"unreadable QLP file {url}: {type(exc).__name__}: {str(exc)[:120]}") from exc
        record = {"status": "ok", "url": url, "sector": int(sector), "tic_id": str(tic_id), "raw_bytes": len(raw),
                  "raw_sha256": hashlib.sha256(raw).hexdigest(), "bin_days": ot.TESS_BIN_DAYS, **reduced}
    man = cache.put(archive="tess", release=BINNED_RELEASE, request_signature=sig,
                    query_text=f"GET {url} -> good-quality SAP_FLUX in 30-min bins", payload=_pack(record),
                    row_count=len(record.get("time", [])), offline_fixture=False)
    return record, man, False


def fetch_star_tess(
    star: dict[str, Any], cache: RawCache, *, max_sectors: int, sleep: ot.Sleeper = time.sleep,
) -> tuple[dict[str, Any], dict[str, Any], int, int]:
    """Light-curve arrays and the fetch record of one star. Returns (lightcurve, record, network requests, cache hits).

    The light curve holds parallel lists `tess_time_btjd`, `tess_flux_qlp` (each sector divided by its median) and
    `tess_sector`. Raises ArchiveAccessError when a lookup fails or every existing sector fails.
    """
    ra, dec = float(star["ra_deg"]), float(star["dec_deg"])
    pmra = float(star.get("pmra_mas_yr", 0.0) or 0.0)
    pmdec = float(star.get("pmdec_mas_yr", 0.0) or 0.0)
    c_ra, c_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2021.0)
    manifests: list[ManifestRecord] = []
    tic, m = ot.lookup_tic_id(ra=c_ra, dec=c_dec, gaia_id=str(star["source_id"]), cache=cache, sleep=sleep)
    manifests.extend(m)
    empty = {"tess_time_btjd": [], "tess_flux_qlp": [], "tess_sector": []}
    if tic is None:
        return empty, {"status": "no_coverage", "note": "no TIC entry at the position", "tic_id": None,
                       "manifest_ids": [x.manifest_id for x in manifests]}, 1, 0
    sectors, m = ot.list_tess_sectors(ra=c_ra, dec=c_dec, cache=cache, sleep=sleep)
    manifests.extend(m)
    used: list[int] = []
    absent: list[int] = []
    short: list[int] = []
    failed: list[str] = []
    skipped: list[int] = []
    t_all: list[float] = []
    f_all: list[float] = []
    s_all: list[int] = []
    n_net = n_hit = 0
    for sector in sectors:
        if len(used) >= max_sectors:
            skipped.append(sector)
            continue
        try:
            rec, man, from_cache = fetch_qlp_sector_binned(tic, sector, cache, sleep=sleep)
        except (ArchiveAccessError, ValueError) as exc:
            failed.append(f"s{sector}: {type(exc).__name__}: {str(exc)[:120]}")
            n_net += 1
            continue
        manifests.append(man)
        if from_cache:
            n_hit += 1
        else:
            n_net += 1
        if rec["status"] != "ok":
            absent.append(sector)
            continue
        if len(rec["time"]) < 20:
            short.append(sector)
            continue
        flux = np.asarray(rec["flux"], dtype=float)
        med = float(np.median(flux))
        if not np.isfinite(med) or med <= 0:
            short.append(sector)
            continue
        used.append(sector)
        t_all.extend(rec["time"])
        f_all.extend((flux / med).round(6).tolist())
        s_all.extend([int(sector)] * len(rec["time"]))
    record = {"tic_id": tic, "sectors_observed": sectors, "sectors_used": used, "sectors_without_qlp": absent,
              "sectors_too_short": short, "sectors_failed": failed, "sectors_beyond_max": skipped, "n_points": len(t_all),
              "manifest_ids": [x.manifest_id for x in manifests]}
    if used:
        record["status"] = "ok"
    elif failed:
        raise ArchiveAccessError(f"every existing TESS sector failed: {'; '.join(failed)[:300]}", retryable=True)
    else:
        record["status"] = "no_coverage"
        record["note"] = "no TESS sector covers the position" if not sectors else "no sector has a usable QLP light curve"
    return {"tess_time_btjd": t_all, "tess_flux_qlp": f_all, "tess_sector": s_all}, record, n_net + 2, n_hit


def enrich_sample_tess(
    stars: list[dict[str, Any]],
    cfg: PipelineConfig,
    cache: RawCache,
    *,
    max_stars: int | None = None,
    progress: ProgressCallback | None = None,
    sleep: Callable[[float], Any] = time.sleep,
    workers: int | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Attach QLP light curves to every sample star. Returns (stars, summary); input order is kept, stars are copies.

    Each processed star gets `lightcurves.tess_time_btjd` / `tess_flux_qlp` / `tess_sector`, a `tess_fetch` record and
    `vetting_archive.query_status["tess_qlp"]` (statuses in the module docstring). Up to
    `channels.optical_dips.tess_bulk_max_sectors` sectors are used per star, the earliest first; the rest are listed
    in `tess_fetch.sectors_beyond_max`. Every answer, including "this sector has no light curve", is cached, so an
    interrupted run resumes without repeating finished work. Failures are not cached.

    `workers` stars are processed at once (default `tess_bulk_workers`), and star fetches from all workers together
    start at least 1 / `tess_bulk_rate_limit_rps` seconds apart. One star is at most two MAST lookups and
    `tess_bulk_max_sectors` file downloads from the public bucket. Results do not depend on the worker count.
    """
    from oddstar.io.neowise_bulk import _RequestPacer

    max_sectors, max_attempts, backoff_s, cfg_workers, rps = _bulk_settings(cfg)
    n_workers = max(1, int(cfg_workers if workers is None else workers))
    pacer = _RequestPacer(1.0 / rps, sleep)
    t_start = time.monotonic()
    n_limit = len(stars) if max_stars is None else max(0, min(int(max_stars), len(stars)))

    def paced_sleep(seconds: float) -> None:
        sleep(seconds)

    def process(star: dict[str, Any]) -> tuple[dict[str, Any], str]:
        s = json.loads(json.dumps(star))
        qs: dict[str, str] = s.setdefault("vetting_archive", {}).setdefault("query_status", {})
        lc = s.setdefault("lightcurves", {})
        status = "failed: not attempted"
        rec: dict[str, Any] = {}
        attempts = 0
        for attempt in range(1, max_attempts + 1):
            attempts = attempt
            t_req = time.monotonic()
            try:
                arrays, rec, n_net, _ = fetch_star_tess(s, cache, max_sectors=max_sectors, sleep=paced_sleep)
                status = str(rec["status"])
                lc.update(arrays)
                if n_net:
                    pacer.after_network_request(t_req)
                break
            except (ArchiveAccessError, ValueError, KeyError, TypeError) as exc:
                status = f"failed: {type(exc).__name__}: {str(exc)[:200]}"
                pacer.after_network_request(t_req)
                if attempt < max_attempts:
                    sleep(backoff_s * (2 ** (attempt - 1)))
        if status.startswith("failed"):
            lc.update({"tess_time_btjd": [], "tess_flux_qlp": [], "tess_sector": []})
            rec = {"status": status}
        rec["attempts"] = attempts
        s["tess_fetch"] = rec
        qs["tess_qlp"] = status
        return s, status

    head = stars[:n_limit]
    if n_workers == 1 or len(head) <= 1:
        results = [process(s) for s in head]
    else:
        with ThreadPoolExecutor(max_workers=n_workers, thread_name_prefix="tess") as pool:
            results = list(pool.map(process, head))

    out: list[dict[str, Any]] = []
    counts = {"ok": 0, "no_coverage": 0, "failed": 0}
    failed_ids: list[str] = []
    n_sectors = 0
    for i, (s, status) in enumerate(results):
        key = "ok" if status == "ok" else ("no_coverage" if status == "no_coverage" else "failed")
        counts[key] += 1
        if key == "failed":
            failed_ids.append(str(s["source_id"]))
        n_sectors += len(s["tess_fetch"].get("sectors_used") or [])
        out.append(s)
        if progress is not None:
            progress(i + 1, n_limit, str(s["source_id"]), status)
    for star in stars[n_limit:]:
        s = json.loads(json.dumps(star))
        s["tess_fetch"] = {"status": "not_fetched", "attempts": 0}
        s.setdefault("lightcurves", {}).update({"tess_time_btjd": [], "tess_flux_qlp": [], "tess_sector": []})
        out.append(s)
    summary = {
        "n_input": len(stars), "n_requested": n_limit, "n_ok": counts["ok"], "n_no_coverage": counts["no_coverage"],
        "n_failed": counts["failed"], "n_not_fetched": len(stars) - n_limit, "failed_ids": failed_ids,
        "n_sectors_used": n_sectors, "max_sectors_per_star": max_sectors, "workers": n_workers,
        "elapsed_s": round(time.monotonic() - t_start, 1),
    }
    return out, summary
