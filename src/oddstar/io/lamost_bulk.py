"""LAMOST low-resolution spectra for every star of a live sample (survey mode of the spectral-line search).

The spectral-line channel compares two independent exposures of one star around H-alpha. No live ingestion provided
spectra, so a live run measured nothing. This module supplies them from the public LAMOST DR12 v2.0 low-resolution
survey (R about 1800, 3700 to 9100 A), which observed many stars on more than one night:

1. the general catalogue is asked which spectra exist for each Gaia source id (batched, cached; every batch also
   asks for a star known to be in the catalogue, so an empty answer is believed only when that star came back);
2. for the best epochs of a star the spectrum file is fetched and reduced on arrival to what the channel reads:
   the window around H-alpha, moved to the star's rest frame in air wavelengths and divided by its continuum. The
   cache keeps the reduced window with the URL, size and SHA-256 of the file it came from.

The two usable epochs with the highest red signal-to-noise, taken on different nights, become `spectra.flux_exp1`
and `spectra.flux_exp2` on a common linear grid.

Statuses of one star, as `spectra_fetch.status` and `vetting_archive.query_status["lamost_lrs"]`:
- "ok": two usable epochs from different nights were loaded;
- "single_epoch": spectra exist but fewer than two usable epochs on different nights, so the repeat test cannot run;
- "no_coverage": the catalogue has no spectrum of this Gaia source (LAMOST observes Dec > -10 deg);
- "failed: ...": the catalogue or a needed spectrum could not be read (not cached, retried next run);
- "not_fetched": beyond `--spectra-max-stars`.
Only "ok" stars are searched; the rest are reported as not searched and excluded from the null result.

Limits that the run must carry: LAMOST has no sky spectrum per fibre in the public file (`sky_flux` stays absent),
and the channel's validation certificate was measured on high-resolution fixtures, not at this resolution.
"""

from __future__ import annotations

import gzip
import hashlib
import io
import json
import re
import time
import xml.etree.ElementTree as ET
import zlib
from collections.abc import Callable
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.schemas import ManifestRecord

LAMOST_RELEASE = "dr12_v2.0"
LAMOST_TAP_URL = "https://www.lamost.org/dr12/v2.0/voservice/tap/sync"
LAMOST_FITS_URL = "https://www.lamost.org/dr12/v2.0/spectrum/fits/{obsid}"
WINDOW_RELEASE = "lrs_halpha_window_v1"  # change when `reduce_lamost_fits` changes what it stores
# A star with spectra in the DR12 v2.0 low-resolution catalogue (three epochs, checked 2026-11-18). It is added to
# every catalogue request; if its rows do not come back, the answer is not trusted.
CANARY_GAIA_ID = "648501102048118528"
CATALOGUE_BATCH = 50
CATALOGUE_COLUMNS = ("obsid", "gaia_source_id", "ra", "dec", "snrr", "lmjd", "obsdate", "class", "subclass", "z", "fibermask")
# LAMOST low-resolution resolving power; the line-spread FWHM handed to the channel is H-alpha / this
LRS_RESOLVING_POWER = 1800.0
HALPHA_AIR_A = 6562.8
# Grid handed to the channel (stellar rest frame, air): the fixture window at close to the native 1.51 A sampling
GRID_START_A, GRID_STOP_A, GRID_STEP_A = 6500.0, 6650.0, 1.5
# The continuum is a straight line fitted to the window outside this half-width around H-alpha (the flux of an
# F to K star falls by about 6 per cent across the window, which a constant would leave in as a tilt)
CONTINUUM_EXCLUDE_HALPHA_A = 15.0
# Points this many robust sigma below / above the line are left out of the next continuum fit (absorption lines pull
# a plain fit down; an emission feature must not pull it up)
CONTINUUM_CLIP_LOW_SIGMA, CONTINUUM_CLIP_HIGH_SIGMA, CONTINUUM_FIT_ITERATIONS = 1.5, 3.0, 4
# An epoch whose window has more than this fraction of masked or zero-weight pixels is not used
MAX_BAD_PIXEL_FRACTION = 0.1
SPEED_OF_LIGHT_KMS = 299_792.458
# Night-sky emission in the window (air, A): the vetting list plus the other OH 6-1 band lines (Osterbrock et al.
# 1996). Positions are good to a few tenths of an angstrom, far inside the exclusion width at this resolution. The
# public file carries no sky spectrum, so every listed line is excluded whether or not it left a residual.
LAMOST_EXTRA_SKY_LINES_A = [6498.73, 6533.04, 6544.03, 6553.62, 6568.78]
_GAIA_ID = re.compile(r"^\d{1,20}$")

ProgressCallback = Callable[[int, int, str, str], None]


def _settings(cfg: PipelineConfig) -> tuple[float, int, float, float]:
    ch = cfg.channels["spectral_lines"]
    return (float(ch["lamost_min_snr_r"]), int(ch["lamost_max_attempts"]), float(ch["lamost_backoff_s"]),
            float(ch["lamost_rate_limit_rps"]))


def vacuum_to_air(wave_vac_a: np.ndarray) -> np.ndarray:
    """IAU standard vacuum-to-air conversion (Morton 2000), valid above 2000 A."""
    s2 = (1e4 / np.asarray(wave_vac_a, dtype=float)) ** 2
    n = 1.0 + 0.0000834254 + 0.02406147 / (130.0 - s2) + 0.00015998 / (38.9 - s2)
    return np.asarray(wave_vac_a, dtype=float) / n


def parse_lamost_votable(payload: bytes, columns: tuple[str, ...] = CATALOGUE_COLUMNS) -> list[dict[str, Any]]:
    """Rows of a LAMOST TAP answer (VOTable TABLEDATA). Raises on anything that is not a readable OK answer."""
    try:
        root = ET.fromstring(payload)
    except ET.ParseError as exc:
        raise ArchiveAccessError(f"LAMOST TAP answer is not XML: {str(exc)[:120]}", retryable=True) from exc
    strip = lambda tag: tag.rsplit("}", 1)[-1]  # noqa: E731 - drop the VOTable namespace
    status = [el for el in root.iter() if strip(el.tag) == "INFO" and el.get("name") == "QUERY_STATUS"]
    if not status or str(status[0].get("value", "")).upper() != "OK":
        detail = (status[0].text or "").strip()[:200] if status else payload[:200].decode("utf-8", errors="replace")
        raise ArchiveAccessError(f"LAMOST TAP query did not answer OK: {detail}", retryable=True)
    fields = [el.get("name") for el in root.iter() if strip(el.tag) == "FIELD"]
    missing = [c for c in columns if c not in fields]
    if missing:
        raise ArchiveAccessError(f"LAMOST TAP answer lacks column(s) {missing}; has {fields}", retryable=False)
    rows: list[dict[str, Any]] = []
    for tr in (el for el in root.iter() if strip(el.tag) == "TR"):
        cells = [(td.text or "").strip() for td in tr if strip(td.tag) == "TD"]
        if len(cells) != len(fields):
            raise ArchiveAccessError(f"LAMOST TAP row has {len(cells)} cells for {len(fields)} fields", retryable=True)
        rows.append({c: cells[fields.index(c)] for c in columns})
    return rows


def _catalogue_batch(ids: list[str], cache: RawCache, *, sleep: ot.Sleeper) -> tuple[list[dict[str, Any]], ManifestRecord | None, bool]:
    asked = sorted(set(ids) | {CANARY_GAIA_ID})
    in_list = ",".join(f"'{i}'" for i in asked)
    query = f"SELECT {','.join(CATALOGUE_COLUMNS)} FROM catalogue WHERE gaia_source_id IN ({in_list})"
    payload, man, hit = ot._cached_get(url=LAMOST_TAP_URL, params={"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": query},
                                       archive="lamost", release=LAMOST_RELEASE, cache=cache, timeout_s=120.0, sleep=sleep)
    try:
        rows = parse_lamost_votable(payload)
        if not any(r["gaia_source_id"] == CANARY_GAIA_ID for r in rows):
            raise ArchiveAccessError("LAMOST catalogue answer lacks the known star asked alongside; an empty answer is not trusted",
                                     retryable=True)
    except ArchiveAccessError:
        if man is not None:
            cache.evict_raw_blob(man.sha256)  # never serve an unreadable or untrusted answer from the cache
        raise
    return rows, man, hit


def lookup_lamost_epochs(
    gaia_ids: list[str], cache: RawCache, *, sleep: ot.Sleeper = time.sleep, min_interval_s: float = 1.0,
) -> tuple[dict[str, list[dict[str, Any]]], dict[str, str], list[str]]:
    """Catalogue epochs of each Gaia source id. Returns (epochs by id, error text by id, manifest ids).

    An id present in `epochs` with an empty list has no spectrum. An id in `errors` could not be looked up (its batch
    failed) and must not be read as "no spectrum".
    """
    for gid in gaia_ids:
        if not _GAIA_ID.match(str(gid)):
            raise ValueError(f"not a Gaia source id: {gid!r}")
    epochs: dict[str, list[dict[str, Any]]] = {}
    errors: dict[str, str] = {}
    manifests: list[str] = []
    todo = [str(g) for g in dict.fromkeys(gaia_ids)]
    for i in range(0, len(todo), CATALOGUE_BATCH):
        batch = todo[i : i + CATALOGUE_BATCH]
        t_req = time.monotonic()
        try:
            rows, man, hit = _catalogue_batch(batch, cache, sleep=sleep)
        except ArchiveAccessError as exc:
            for gid in batch:
                errors[gid] = f"{type(exc).__name__}: {str(exc)[:200]}"
            continue
        if man is not None:
            manifests.append(man.manifest_id)
        for gid in batch:
            epochs[gid] = []
        for r in rows:
            gid = r["gaia_source_id"]
            if gid not in epochs:
                continue  # the canary, unless it is itself a sample star
            try:
                epochs[gid].append({"obsid": int(r["obsid"]), "lmjd": int(r["lmjd"]), "obsdate": r["obsdate"], "snrr": float(r["snrr"]),
                                    "z": float(r["z"]), "class": r["class"], "subclass": r["subclass"], "fibermask": int(r["fibermask"]),
                                    "ra_deg": float(r["ra"]), "dec_deg": float(r["dec"])})
            except ValueError:
                continue  # a row without a usable obsid, date, S/N or redshift is not an epoch
        if not hit:
            wait = min_interval_s - (time.monotonic() - t_req)
            if wait > 0:
                sleep(wait)
    for gid in epochs:
        epochs[gid].sort(key=lambda e: (-e["snrr"], e["obsid"]))
    return epochs, errors, manifests


def reduce_lamost_fits(fits_bytes: bytes) -> dict[str, Any]:
    """The H-alpha window of one LAMOST low-resolution file, in the stellar rest frame (air), continuum-normalised.

    Returns {"usable", "reason", "wavelength_a", "flux", "n_bad_pixels", "z", "helio_rv_kms", "snrr", "lmjd", "date_obs",
    "obsid", "n_exposures"}. Masked or zero-weight pixels are replaced by linear interpolation and counted. The epoch
    is not usable when the file has no pipeline redshift, the window is not covered, or too many pixels are bad.
    """
    from astropy.io import fits

    raw = gzip.decompress(fits_bytes) if fits_bytes[:2] == b"\x1f\x8b" else fits_bytes
    with fits.open(io.BytesIO(raw)) as hdul:
        hdr = hdul[0].header
        tab = hdul[1].data
        wave_vac = np.asarray(tab["WAVELENGTH"][0], dtype=float)
        flux = np.asarray(tab["FLUX"][0], dtype=float)
        ivar = np.asarray(tab["IVAR"][0], dtype=float)
        ormask = np.asarray(tab["ORMASK"][0], dtype=float)
        meta = {"obsid": int(hdr.get("OBSID", 0)), "z": float(hdr.get("Z", -9999.0)), "helio_rv_kms": float(hdr.get("HELIO_RV", 0.0) or 0.0),
                "snrr": float(hdr.get("SNRR", 0.0) or 0.0), "lmjd": int(hdr.get("LMJD", 0) or 0), "date_obs": str(hdr.get("DATE-OBS", "")),
                "n_exposures": int(hdr.get("NEXP", 0) or 0), "vacuum": bool(hdr.get("VACUUM", True))}
    out: dict[str, Any] = dict(meta, usable=False, reason="", wavelength_a=[], flux=[], n_bad_pixels=0)
    if not np.isfinite(meta["z"]) or abs(meta["z"]) > 0.01:
        out["reason"] = f"no usable pipeline redshift (z={meta['z']})"  # stars: |v| < 3000 km/s; -9999 marks a failed fit
        return out
    wave = (vacuum_to_air(wave_vac) if meta["vacuum"] else wave_vac) / (1.0 + meta["z"])
    order = np.argsort(wave)
    wave, flux, ivar, ormask = wave[order], flux[order], ivar[order], ormask[order]
    sel = (wave >= GRID_START_A - 2 * GRID_STEP_A) & (wave <= GRID_STOP_A + 2 * GRID_STEP_A)
    if sel.sum() < 20 or wave[sel][0] > GRID_START_A or wave[sel][-1] < GRID_STOP_A:
        out["reason"] = "the H-alpha window is not covered"
        return out
    w, f = wave[sel], flux[sel]
    bad = ~np.isfinite(f) | (ivar[sel] <= 0) | (ormask[sel] > 0)
    out["n_bad_pixels"] = int(bad.sum())
    if bad.mean() > MAX_BAD_PIXEL_FRACTION:
        out["reason"] = f"{int(bad.sum())} of {len(w)} window pixels are masked"
        return out
    if bad.any():
        f = np.interp(w, w[~bad], f[~bad])
    cont = fit_linear_continuum(w, f)
    if cont is None:
        out["reason"] = "continuum is not positive"
        return out
    out.update(usable=True, wavelength_a=[round(float(v), 4) for v in w], flux=[round(float(v), 5) for v in f / cont])
    return out


def fit_linear_continuum(wave: np.ndarray, flux: np.ndarray) -> np.ndarray | None:
    """Straight-line continuum of the window, fitted away from H-alpha with low and high points clipped; None if not positive."""
    keep = np.abs(wave - HALPHA_AIR_A) > CONTINUUM_EXCLUDE_HALPHA_A
    x = wave - HALPHA_AIR_A
    line = np.full_like(flux, float(np.median(flux[keep])))
    for _ in range(CONTINUUM_FIT_ITERATIONS):
        if keep.sum() < 10:
            break
        slope, icpt = np.polyfit(x[keep], flux[keep], 1)
        line = slope * x + icpt
        resid = flux - line
        sig = 1.4826 * float(np.median(np.abs(resid[keep] - np.median(resid[keep]))))
        if sig <= 0:
            break
        keep = (np.abs(x) > CONTINUUM_EXCLUDE_HALPHA_A) & (resid > -CONTINUUM_CLIP_LOW_SIGMA * sig) & (resid < CONTINUUM_CLIP_HIGH_SIGMA * sig)
    if not np.all(np.isfinite(line)) or float(np.min(line)) <= 0:
        return None
    return line


def fetch_lamost_window(obsid: int, cache: RawCache, *, sleep: ot.Sleeper = time.sleep) -> tuple[dict[str, Any], ManifestRecord, bool]:
    """The reduced H-alpha window of one spectrum. Returns (record, manifest, from_cache); failures are not stored."""
    url = LAMOST_FITS_URL.format(obsid=int(obsid))
    sig = hashlib.sha256(f"GET {url} -> halpha window".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("lamost", WINDOW_RELEASE, sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(gzip.decompress(cache.read_bytes(cached)).decode("utf-8")), cached, True
        except (FileNotFoundError, OSError, ValueError):
            pass  # evicted or unreadable blob: fetch again
    resp = ot._get_with_retries(url, None, timeout_s=60.0, attempts=3, sleep=sleep)
    raw = resp.content
    try:
        record = reduce_lamost_fits(raw)
    except (KeyError, OSError, ValueError, TypeError, IndexError, EOFError, zlib.error) as exc:
        raise ValueError(f"unreadable LAMOST file {url}: {type(exc).__name__}: {str(exc)[:120]}") from exc
    record.update(url=url, raw_bytes=len(raw), raw_sha256=hashlib.sha256(raw).hexdigest())
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb", mtime=0) as gz:  # mtime=0: identical input gives identical bytes
        gz.write(json.dumps(record, sort_keys=True, separators=(",", ":")).encode("utf-8"))
    man = cache.put(archive="lamost", release=WINDOW_RELEASE, request_signature=sig, query_text=f"GET {url} -> H-alpha window, rest frame",
                    payload=buf.getvalue(), row_count=len(record["wavelength_a"]), offline_fixture=False)
    return record, man, False


def sky_lines_in_rest_frame(windows: list[dict[str, Any]]) -> list[float]:
    """Where the night-sky lines fall in each epoch's stellar rest frame (air wavelengths).

    The spectra are heliocentric and were moved to the star's rest frame, so a sky line observed at L sits at
    L x (1 + v_helio / c) / (1 + z). For a fast halo star that is several angstroms from the laboratory wavelength.
    """
    from oddstar.vet.spectral_checks import KNOWN_SKY_LINES_A

    out = []
    for w in windows:
        factor = (1.0 + float(w.get("helio_rv_kms", 0.0)) / SPEED_OF_LIGHT_KMS) / (1.0 + float(w["z"]))
        out += [round(line * factor, 3) for line in [*KNOWN_SKY_LINES_A, *LAMOST_EXTRA_SKY_LINES_A]]
    return sorted(set(out))


def fetch_star_spectra(
    star: dict[str, Any], epochs: list[dict[str, Any]], cache: RawCache, *, min_snr_r: float, sleep: ot.Sleeper = time.sleep,
) -> tuple[dict[str, Any], dict[str, Any], int]:
    """The channel's `spectra` record and the fetch record of one star. Returns (spectra, record, network requests).

    Epochs are tried from the highest red S/N down, one per night, until two usable ones are loaded. Raises
    ArchiveAccessError when a needed file could not be fetched, so a star is never called single-epoch because of
    a failed download.
    """
    record: dict[str, Any] = {"n_epochs_catalogue": len(epochs), "epochs": [], "obsids_used": [], "manifest_ids": []}
    if not epochs:
        return {}, dict(record, status="no_coverage", note="no spectrum of this Gaia source in the LAMOST low-resolution catalogue"), 0
    used: list[dict[str, Any]] = []
    nights: set[int] = set()
    failures: list[str] = []
    n_net = 0
    for ep in sorted(epochs, key=lambda e: (-e["snrr"], e["obsid"])):
        row = {"obsid": ep["obsid"], "lmjd": ep["lmjd"], "obsdate": ep["obsdate"], "snrr": ep["snrr"], "used": False, "note": ""}
        record["epochs"].append(row)
        if len(used) >= 2:
            row["note"] = "not needed"
            continue
        if ep["snrr"] < min_snr_r:
            row["note"] = f"red S/N below {min_snr_r:g}"
            continue
        if ep["lmjd"] in nights:
            row["note"] = "same night as a used epoch"
            continue
        try:
            win, man, hit = fetch_lamost_window(ep["obsid"], cache, sleep=sleep)
        except (ArchiveAccessError, ValueError) as exc:
            failures.append(f"obsid {ep['obsid']}: {type(exc).__name__}: {str(exc)[:120]}")
            row["note"] = "fetch failed"
            n_net += 1
            continue
        n_net += 0 if hit else 1
        record["manifest_ids"].append(man.manifest_id)
        row.update(n_bad_pixels=win["n_bad_pixels"], z=win["z"], helio_rv_kms=win["helio_rv_kms"], raw_sha256=win["raw_sha256"], url=win["url"])
        if not win["usable"]:
            row["note"] = win["reason"]
            continue
        row["used"] = True
        used.append(win)
        nights.add(ep["lmjd"])
    record["n_epochs_usable"] = len(used)
    if len(used) < 2:
        if failures:
            raise ArchiveAccessError(f"LAMOST spectra could not be fetched: {'; '.join(failures)[:300]}", retryable=True)
        return {}, dict(record, status="single_epoch", note=f"{len(used)} usable epoch(s) on different nights; the repeat test needs two"), n_net
    grid = np.arange(GRID_START_A, GRID_STOP_A + 0.5 * GRID_STEP_A, GRID_STEP_A)
    fluxes = [np.interp(grid, np.asarray(w["wavelength_a"]), np.asarray(w["flux"])) for w in used]
    record.update(status="ok", obsids_used=[w["obsid"] for w in used])
    spectra = {
        "wavelength_a": [round(float(v), 3) for v in grid],
        "flux_exp1": [round(float(v), 5) for v in fluxes[0]],
        "flux_exp2": [round(float(v), 5) for v in fluxes[1]],
        "lsf_fwhm_a": round(HALPHA_AIR_A / LRS_RESOLVING_POWER, 3),
        "sky_line_wavelengths_a": sky_lines_in_rest_frame(used),
        # rest wavelength x factor = heliocentric observed wavelength of that epoch
        "rest_to_observed_factors": [round(1.0 + float(w["z"]), 8) for w in used],
        "source": f"LAMOST {LAMOST_RELEASE} low-resolution, obsid {used[0]['obsid']} and {used[1]['obsid']}",
        "frame": "stellar rest frame, air wavelengths, continuum-normalised",
    }
    return spectra, record, n_net


def enrich_sample_lamost(
    stars: list[dict[str, Any]],
    cfg: PipelineConfig,
    cache: RawCache,
    *,
    max_stars: int | None = None,
    progress: ProgressCallback | None = None,
    sleep: Callable[[float], Any] = time.sleep,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Attach two LAMOST epochs around H-alpha to every sample star that has them. Returns (stars, summary).

    Input order is kept and stars are copies. Each processed star gets `spectra` (empty unless "ok"), a
    `spectra_fetch` record, `vetting_archive.query_status["lamost_lrs"]` and, when "ok",
    `vetting_archive.query_status["spectra"] = "ok"` so the spectral vetting checks read the arrays. Requests are
    sequential and at least 1 / `lamost_rate_limit_rps` seconds apart; every answer is cached, so an interrupted
    run resumes. Failures are not cached.
    """
    min_snr_r, max_attempts, backoff_s, rps = _settings(cfg)
    interval = 1.0 / rps
    t_start = time.monotonic()
    n_limit = len(stars) if max_stars is None else max(0, min(int(max_stars), len(stars)))
    head = [json.loads(json.dumps(s)) for s in stars[:n_limit]]
    ids = [str(s["source_id"]) for s in head]
    epochs: dict[str, list[dict[str, Any]]] = {}
    errors: dict[str, str] = {}
    for attempt in range(1, max_attempts + 1):
        todo = [i for i in ids if i not in epochs]
        if not todo:
            break
        got, errors, _ = lookup_lamost_epochs(todo, cache, sleep=sleep, min_interval_s=interval)
        epochs.update(got)
        if errors and attempt < max_attempts:
            sleep(backoff_s * (2 ** (attempt - 1)))

    out: list[dict[str, Any]] = []
    counts = {"ok": 0, "single_epoch": 0, "no_coverage": 0, "failed": 0}
    failed_ids: list[str] = []
    for i, s in enumerate(head):
        sid = str(s["source_id"])
        qs: dict[str, str] = s.setdefault("vetting_archive", {}).setdefault("query_status", {})
        spectra: dict[str, Any] = {}
        rec: dict[str, Any] = {}
        if sid not in epochs:
            status = f"failed: catalogue lookup: {errors.get(sid, 'not attempted')}"
        else:
            status = "failed: not attempted"
            for attempt in range(1, max_attempts + 1):
                t_req = time.monotonic()
                try:
                    spectra, rec, n_net = fetch_star_spectra(s, epochs[sid], cache, min_snr_r=min_snr_r, sleep=sleep)
                    status = str(rec["status"])
                except (ArchiveAccessError, ValueError, KeyError, TypeError) as exc:
                    status, n_net = f"failed: {type(exc).__name__}: {str(exc)[:200]}", 1
                if n_net:
                    wait = interval - (time.monotonic() - t_req)
                    if wait > 0:
                        sleep(wait)
                if not status.startswith("failed"):
                    break
                if attempt < max_attempts:
                    sleep(backoff_s * (2 ** (attempt - 1)))
        if status.startswith("failed"):
            spectra, rec = {}, {"status": status}
            failed_ids.append(sid)
        s["spectra"] = spectra
        s["spectra_fetch"] = rec
        qs["lamost_lrs"] = status
        if status == "ok":
            qs["spectra"] = "ok"
        counts["failed" if status.startswith("failed") else status] += 1
        out.append(s)
        if progress is not None:
            progress(i + 1, n_limit, sid, status)
    for star in stars[n_limit:]:
        s = json.loads(json.dumps(star))
        s["spectra"] = {}
        s["spectra_fetch"] = {"status": "not_fetched"}
        out.append(s)
    summary = {
        "n_input": len(stars), "n_requested": n_limit, "n_ok": counts["ok"], "n_single_epoch": counts["single_epoch"],
        "n_no_coverage": counts["no_coverage"], "n_failed": counts["failed"], "n_not_fetched": len(stars) - n_limit,
        "failed_ids": failed_ids, "min_snr_r": min_snr_r, "survey": f"LAMOST {LAMOST_RELEASE} low-resolution",
        "elapsed_s": round(time.monotonic() - t_start, 1),
    }
    return out, summary
