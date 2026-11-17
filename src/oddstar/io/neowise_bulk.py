"""NEOWISE-R-only bulk light-curve ingest for survey-scale `ir_timedomain` searches.

`enrich_star_live_for_vetting` runs a dozen archive queries per star and is meant for the few flagged targets. A
time-domain search of a whole sample needs only the NEOWISE-R single-exposure rows of each star, binned into visits.
This module provides:

- `fetch_neowise_single_exposures`: the IRSA Gator cone (fast path) -> IRSA TAP (backup) query that vetting section
  "3c" also uses, with the same status strings ("ok", "no_coverage", "failed: ...").
- `bin_neowise_visits`: the visit binning shared with section "3c" (per-visit median fluxes, robust errors, Gaia
  proper-motion track residual vectors, `visit_n_frames`).
- `enrich_sample_neowise`: the survey loop, with the content-addressed `RawCache` (a re-run is free), a few concurrent
  workers under one shared rate limit from `cfg.archives["neowise"]`, bounded retries, a per-star status, resumability through the cache, and a progress
  callback. A star whose fetch failed is left with an empty light curve and status `failed: ...`, so the channel
  reports it as not searched rather than scoring it as quiet.

Archive helpers are looked up on `oddstar.io.live_survey` at call time (not bound at import), so tests that replace
them or `requests` there apply here too.
"""

from __future__ import annotations

import json
import math
import threading
import time
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec, mag_to_flux_mjy, propagate_position
from oddstar.schemas import ManifestRecord

NEOWISE_CATALOG = "neowiser_p1bs_psd"
# Target cone radius used by vetting section 3c (3" Gator cone; 0.00083 deg TAP circle)
NEOWISE_CONE_RADIUS_ARCSEC = 3.0
_TAP_CIRCLE_DEG = 0.00083
# Columns of the section 3c TAP backup query (kept verbatim so its cache signature is unchanged)
VETTING_TAP_COLUMNS = "mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, qual_frame, qi_fact"
# Survey mode keeps only the columns the visit binning reads (plus sigmas/qi_fact for provenance), which cuts the cached
# payload by roughly an order of magnitude relative to Gator's default column set.
BULK_GATOR_SELCOLS = "ra,dec,mjd,w1mpro,w1sigmpro,w2mpro,w2sigmpro,qual_frame,qi_fact"
BULK_TAP_COLUMNS = "ra, dec, mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, qual_frame, qi_fact"
# Frames separated by more than this many days start a new visit (NEOWISE visits are ~1 week every ~6 months)
VISIT_GAP_DAYS = 120.0

EMPTY_NEOWISE: dict[str, list[float]] = {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []}

ProgressCallback = Callable[[int, int, str, str], None]


def fetch_neowise_single_exposures(
    ra: float,
    dec: float,
    cache: RawCache,
    *,
    gator_selcols: str | None = None,
    tap_columns: str = VETTING_TAP_COLUMNS,
) -> tuple[list[dict[str, Any]], list[ManifestRecord], str, bool]:
    """NEOWISE-R single-exposure rows within 3" of (ra, dec): IRSA Gator first, IRSA TAP as backup.

    Returns (rows, manifests, status, from_cache). `status` is "ok", "no_coverage" or "failed: <type>: <msg>"
    (with the Gator error appended when both paths failed), exactly as vetting section 3c records it.
    """
    from oddstar.io import live_survey as ls

    manifests: list[ManifestRecord] = []
    rows: list[dict[str, Any]] = []
    primary_err: str | None = None
    from_cache = False
    try:
        try:
            rows, m_neo, from_cache = ls._gator_cone_cached(
                catalog=NEOWISE_CATALOG, ra=ra, dec=dec, radius_arcsec=NEOWISE_CONE_RADIUS_ARCSEC, archive="neowise",
                release="neowise_r", cache=cache, timeout_s=25, selcols=gator_selcols,
            )
            manifests.append(m_neo)
        except Exception as exc_primary:
            primary_err = ls._failed(exc_primary)  # kept so a backup failure reports both
            q_neo = (
                f"SELECT {tap_columns} "
                f"FROM {NEOWISE_CATALOG} WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {_TAP_CIRCLE_DEG})) "
                f"AND qual_frame > 0"
            )
            rows, m_neo, from_cache = ls._tap_query_cached(
                url=ls.IRSA_TAP_URL, query=q_neo, archive="neowise", release="neowise_r", cache=cache, timeout_s=30
            )
            manifests.append(m_neo)
        status = "ok" if rows else "no_coverage"
    except Exception as exc:
        status = ls._failed(exc) + (f" (gator primary {primary_err})" if primary_err else "")
        rows = []
    return rows, manifests, status, from_cache


def bin_neowise_visits(
    rows: list[dict[str, Any]],
    *,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
) -> tuple[dict[str, Any], list[tuple[float, float]], float] | None:
    """Bin NEOWISE-R single-exposure rows into visits.

    Returns (neowise_dict, visit_windows, target_w1_mag), or None when no row passes the frame-quality and photometry
    filters. `neowise_dict` holds per-visit `mjd`, median `w1_mjy`/`w2_mjy`, robust errors, the Gaia proper-motion
    track residuals (median and per-visit RA*cos(Dec)/Dec vectors, when the rows carry positions) and `visit_n_frames`.
    `visit_windows` are the first/last frame MJD of each visit (used for the vetting comparison-star ensemble) and
    `target_w1_mag` the median single-frame W1 magnitude.
    """
    from oddstar.io.live_survey import _safe_float

    pts: list[tuple[float, float, float, float | None, float | None]] = []
    for row in rows:
        qf = _safe_float(row.get("qual_frame"))
        mjd_val = _safe_float(row.get("mjd"))
        w1_val = _safe_float(row.get("w1mpro"))
        w2_val = _safe_float(row.get("w2mpro"))
        r_ra = _safe_float(row.get("ra"))
        r_dec = _safe_float(row.get("dec"))
        if (qf is None or qf > 0) and mjd_val is not None and w1_val is not None and w2_val is not None:
            pts.append((mjd_val, w1_val, w2_val, r_ra, r_dec))
    pts.sort(key=lambda x: x[0])
    if not pts:
        return None

    # Bin into ~180-day semi-annual visits and check PM-corrected astrometric trajectory
    bins: list[tuple[float, list[tuple[float, float, float, float | None, float | None]]]] = []
    for item in pts:
        if not bins or item[0] - bins[-1][0] > VISIT_GAP_DAYS:
            bins.append((item[0], [item]))
        else:
            bins[-1][1].append(item)
    v_mjd, v_w1, v_w2, v_e1, v_e2, pm_resids = [], [], [], [], [], []
    v_dra: list[float | None] = []
    v_ddec: list[float | None] = []
    v_win: list[tuple[float, float]] = []
    v_nfr: list[int] = []
    for _, grp in bins:
        mjd_med = float(np.median([x[0] for x in grp]))
        v_mjd.append(mjd_med)
        v_win.append((float(grp[0][0]), float(grp[-1][0])))
        v_nfr.append(len(grp))
        f1s = [mag_to_flux_mjy("W1", x[1]) for x in grp]
        f2s = [mag_to_flux_mjy("W2", x[2]) for x in grp]
        m1 = float(np.median(f1s))
        m2 = float(np.median(f2s))
        v_w1.append(m1)
        v_w2.append(m2)
        v_e1.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f1s) - m1))) / math.sqrt(len(grp)), m1 * 0.02))
        v_e2.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f2s) - m2))) / math.sqrt(len(grp)), m2 * 0.02))
        ras_g = [x[3] for x in grp if x[3] is not None]
        decs_g = [x[4] for x in grp if x[4] is not None]
        if ras_g and decs_g:
            jyr = 2000.0 + (mjd_med - 51544.5) / 365.25
            pra, pdec = propagate_position(ra, dec, pmra, pmdec, 2016.0, jyr)
            r_med, d_med = float(np.median(ras_g)), float(np.median(decs_g))
            pm_resids.append(angular_separation_arcsec(r_med, d_med, pra, pdec))
            # Per-visit residual vector from the Gaia PM track (RA*cos(Dec), Dec; arcsec) for neighbour_signal
            v_dra.append(round(((r_med - pra + 540.0) % 360.0 - 180.0) * math.cos(math.radians(pdec)) * 3600.0, 5))
            v_ddec.append(round((d_med - pdec) * 3600.0, 5))
        else:
            v_dra.append(None)
            v_ddec.append(None)
    neo_dict: dict[str, Any] = {"mjd": v_mjd, "w1_mjy": v_w1, "w2_mjy": v_w2, "w1_err_mjy": v_e1, "w2_err_mjy": v_e2}
    if pm_resids:
        neo_dict["median_pm_residual_arcsec"] = float(np.median(pm_resids))
        neo_dict["pm_residual_ra_arcsec"] = v_dra
        neo_dict["pm_residual_dec_arcsec"] = v_ddec
    neo_dict["visit_n_frames"] = v_nfr
    return neo_dict, v_win, float(np.median([x[1] for x in pts]))


def _bulk_settings(cfg: PipelineConfig) -> tuple[float, int, float, int]:
    arch = cfg.archives.get("neowise")
    if arch is None or not arch.enabled:
        raise ArchiveAccessError(
            "NEOWISE archive is not enabled in cfg.archives; cannot run the survey NEOWISE fetch.",
            hint="Enable `archives.neowise` in the pipeline config.",
            retryable=False,
        )
    itd = cfg.channels["ir_timedomain"]
    return (float(arch.rate_limit_rps), int(itd["neowise_bulk_max_attempts"]), float(itd["neowise_bulk_retry_backoff_s"]),
            int(itd["neowise_bulk_workers"]))


class _RequestPacer:
    """Spaces network requests at least `min_interval` apart across all workers.

    After a request that reached the network, the worker reserves the next free start slot and waits for it. With one
    worker this is the plain "wait until `min_interval` has passed since the request started".
    """

    def __init__(self, min_interval: float, sleep: Callable[[float], Any]) -> None:
        self._min_interval = min_interval
        self._sleep = sleep
        self._lock = threading.Lock()
        self._next_start = 0.0

    def after_network_request(self, t_req: float) -> None:
        with self._lock:
            slot = max(self._next_start, t_req + self._min_interval)
            self._next_start = slot + self._min_interval
        wait = slot - time.monotonic()
        if wait > 0:
            self._sleep(wait)


def enrich_sample_neowise(
    stars: list[dict[str, Any]],
    cfg: PipelineConfig,
    cache: RawCache,
    *,
    max_stars: int | None = None,
    progress: ProgressCallback | None = None,
    sleep: Callable[[float], Any] = time.sleep,
    workers: int | None = None,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Attach NEOWISE-R visit light curves to every sample star (NEOWISE queries only).

    Returns (stars, summary). All input stars are returned in their original order, as copies. Each processed star
    gets `neowise` (binned visits; empty lists when nothing usable), `vetting_archive.query_status["neowise"]` and a
    `neowise_fetch` record {status, attempts, n_rows, n_visits, from_cache, manifest_ids}. Status is "ok",
    "no_coverage" (no rows, or no row passing the frame-quality filter) or "failed: ..." after
    `channels.ir_timedomain.neowise_bulk_max_attempts` attempts with exponential backoff. Stars past `max_stars`
    are returned untouched with `neowise_fetch.status == "not_fetched"`; stars already through full vetting
    enrichment (`live_enriched`) keep their data and status.

    Concurrency: `workers` requests run at once (default `channels.ir_timedomain.neowise_bulk_workers`). One IRSA
    cone takes about three seconds to answer, so a single worker needs about seven hours for 7,800 stars; the
    answers do not depend on the order they are fetched in, and results are always returned in input order.

    Rate limit: network requests from all workers together start at least 1 / `archives.neowise.rate_limit_rps`
    seconds apart; cache hits do not wait. Resumability: every completed query (including empty results) is in the
    content-addressed cache, so re-running after an interruption replays finished stars from the cache without
    network traffic or waiting and continues with the unfinished ones. Failed queries are not cached and are retried
    on the next run.

    `progress(done, n_requested, source_id, status)` is called in input order after each processed star.
    """
    rps, max_attempts, backoff_s, cfg_workers = _bulk_settings(cfg)
    n_workers = max(1, int(cfg_workers if workers is None else workers))
    pacer = _RequestPacer(1.0 / rps, sleep)
    t_start = time.monotonic()
    n_limit = len(stars) if max_stars is None else max(0, min(int(max_stars), len(stars)))

    def process(idx: int, star: dict[str, Any]) -> tuple[dict[str, Any], str, str, int, int, bool]:
        """One star -> (copy with results, count key, status, network requests, cache hits, from vetting enrichment)."""
        s = json.loads(json.dumps(star))
        qs: dict[str, str] = s.setdefault("vetting_archive", {}).setdefault("query_status", {})
        if s.get("live_enriched") and "neowise" in qs:
            status = qs["neowise"]
            key = "ok" if status == "ok" else ("no_coverage" if status == "no_coverage" else "failed")
            n_vis = len((s.get("neowise") or {}).get("mjd") or [])
            if key == "ok" and n_vis == 0:
                key, status = "no_coverage", "no_coverage"
            s["neowise_fetch"] = {"status": status, "attempts": 0, "n_visits": n_vis, "source": "vetting_enrichment"}
            return s, key, status, 0, 0, True

        ra = float(s["ra_deg"])
        dec = float(s["dec_deg"])
        pmra = float(s.get("pmra_mas_yr", 0.0) or 0.0)
        pmdec = float(s.get("pmdec_mas_yr", 0.0) or 0.0)

        rows: list[dict[str, Any]] = []
        manifests: list[ManifestRecord] = []
        status = "failed: not attempted"
        from_cache = False
        attempts = 0
        n_net = 0
        n_hit = 0
        for attempt in range(1, max_attempts + 1):
            attempts = attempt
            t_req = time.monotonic()
            rows, manifests, status, from_cache = fetch_neowise_single_exposures(
                ra, dec, cache, gator_selcols=BULK_GATOR_SELCOLS, tap_columns=BULK_TAP_COLUMNS
            )
            if from_cache:
                n_hit += 1
            else:
                n_net += 1
                pacer.after_network_request(t_req)
            if not status.startswith("failed"):
                break
            if attempt < max_attempts:
                sleep(backoff_s * (2 ** (attempt - 1)))

        binned = bin_neowise_visits(rows, ra=ra, dec=dec, pmra=pmra, pmdec=pmdec) if rows else None
        if status == "ok" and binned is None:
            status = "no_coverage"
            note = f"{len(rows)} rows, none passing qual_frame > 0 with W1 and W2 photometry"
        else:
            note = ""
        s["neowise"] = binned[0] if (binned is not None and status == "ok") else json.loads(json.dumps(EMPTY_NEOWISE))
        qs["neowise"] = status
        rec: dict[str, Any] = {
            "status": status,
            "attempts": attempts,
            "n_rows": len(rows),
            "n_visits": len(s["neowise"]["mjd"]),
            "from_cache": bool(from_cache),
            "manifest_ids": [m.manifest_id for m in manifests],
        }
        if note:
            rec["note"] = note
        s["neowise_fetch"] = rec
        key = "ok" if status == "ok" else ("no_coverage" if status == "no_coverage" else "failed")
        return s, key, status, n_net, n_hit, False

    out: list[dict[str, Any]] = []
    counts = {"ok": 0, "no_coverage": 0, "failed": 0}
    failed_ids: list[str] = []
    n_network = 0
    n_cache_hits = 0
    n_already = 0

    def collect(idx: int, res: tuple[dict[str, Any], str, str, int, int, bool]) -> None:
        nonlocal n_network, n_cache_hits, n_already
        s, key, status, n_net, n_hit, already = res
        sid = str(s["source_id"])
        counts[key] += 1
        if key == "failed":
            failed_ids.append(sid)
        n_network += n_net
        n_cache_hits += n_hit
        n_already += 1 if already else 0
        out.append(s)
        if progress is not None:
            progress(idx + 1, n_limit, sid, status)

    todo = stars[:n_limit]
    if n_workers == 1 or len(todo) <= 1:
        for idx, star in enumerate(todo):
            collect(idx, process(idx, star))
    else:
        pool = ThreadPoolExecutor(max_workers=n_workers, thread_name_prefix="neowise")
        try:
            futures = [pool.submit(process, idx, star) for idx, star in enumerate(todo)]
            for idx, fut in enumerate(futures):
                collect(idx, fut.result())
        except BaseException:
            # Interrupted (or a programming error): drop the queued stars; finished ones are already in the cache
            pool.shutdown(wait=False, cancel_futures=True)
            raise
        pool.shutdown(wait=True)
    for star in stars[n_limit:]:
        s = json.loads(json.dumps(star))
        s["neowise_fetch"] = {"status": "not_fetched", "reason": f"beyond max_stars={max_stars}"}
        out.append(s)

    summary: dict[str, Any] = {
        "n_input": len(stars),
        "n_requested": n_limit,
        "n_ok": counts["ok"],
        "n_no_coverage": counts["no_coverage"],
        "n_failed": counts["failed"],
        "n_not_fetched": len(stars) - n_limit,
        "failed_ids": failed_ids,
        "n_network_requests": n_network,
        "n_cache_hits": n_cache_hits,
        "n_from_vetting_enrichment": n_already,
        "rate_limit_rps": rps,
        "max_attempts": max_attempts,
        "workers": n_workers,
        "elapsed_s": round(time.monotonic() - t_start, 3),
    }
    return out, summary
