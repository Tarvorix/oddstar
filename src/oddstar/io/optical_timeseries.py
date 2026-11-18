"""Live optical time-series ingestion for light-curve vetting (TESS QLP/TESScut, ZTF, ASAS-SN).

These products feed the light-curve vetting checks only (`pipeline_agreement`, `systematics_window`,
`neighbour_signal`, `sso_crossing`, `known_variable`, `second_survey`). They are stored under
`star["lightcurves"]["tess_sectors"]`, `star["lightcurves"]["optical_surveys"]`, `star["lightcurves"]["tess_neighbours"]`,
`star["lightcurves"]["sso_queries"]` and `star["neowise"]["comparison"]`, never under the channel inputs
(`tess_flux_qlp`, `ztf_flux`, NEOWISE visit fluxes), so that adding vetting data never changes a certified channel score.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import time
import zipfile
from collections.abc import Callable
from typing import Any

import numpy as np
import requests

from oddstar.errors import ArchiveAccessError
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec, propagate_position
from oddstar.schemas import ManifestRecord

IRSA_ZTF_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
MAST_INVOKE_URL = "https://mast.stsci.edu/api/v0/invoke"
TESSCUT_SECTOR_URL = "https://mast.stsci.edu/tesscut/api/v0.1/sector"
TESSCUT_ASTROCUT_URL = "https://mast.stsci.edu/tesscut/api/v0.1/astrocut"
QLP_S3_BASE = "https://stpubdata.s3.amazonaws.com/mast/hlsp/qlp"

BTJD_OFFSET = 2457000.0
ZTF_FIRST_MJD = 58194.0  # 2018-03-17, start of ZTF public survey operations
ZTF_CONE_LADDER_ARCSEC = (3.0, 2.0)
ZTF_MATCH_RADIUS_ARCSEC = 1.5
ASASSN_CONE_ARCSEC = 15.0
ASASSN_MATCH_RADIUS_ARCSEC = 5.0
TESSCUT_SIZE_PX = 13
TESS_APERTURE_RADIUS_PX = 1.6
TESS_ANNULUS_PX = (4.0, 6.0)
TESS_CENTROID_RADIUS_PX = 2.5  # flux-weighted centroid region (extraction geometry, like the aperture radius)
TESS_BIN_DAYS = 1.0 / 48.0  # 30-minute bins common to every QLP/FFI cadence (1800 s, 600 s, 200 s)

Sleeper = Callable[[float], None]


def _mjd_to_jyear(mjd: float) -> float:
    return 2000.0 + (mjd - 51544.5) / 365.25


def _get_with_retries(
    url: str,
    params: dict[str, Any] | None,
    *,
    timeout_s: float,
    attempts: int = 3,
    backoff_s: float = 3.0,
    sleep: Sleeper = time.sleep,
    ok_statuses: tuple[int, ...] = (200,),
) -> requests.Response:
    """GET with bounded retries and exponential backoff on timeouts, connection errors and 5xx/429."""
    errors: list[str] = []
    for i in range(attempts):
        try:
            resp = requests.get(url, params=params, timeout=timeout_s)
            if resp.status_code in ok_statuses:
                return resp
            errors.append(f"HTTP {resp.status_code}")
            if resp.status_code < 500 and resp.status_code != 429:
                break
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
            errors.append(f"{type(exc).__name__}")
        if i < attempts - 1:
            sleep(backoff_s * (2**i))
    raise ArchiveAccessError(f"GET {url} failed after {len(errors)} attempt(s): {'; '.join(errors)}", retryable=True)


def _cached_get(
    *,
    url: str,
    params: dict[str, Any],
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: float,
    attempts: int = 3,
    sleep: Sleeper = time.sleep,
    ok_statuses: tuple[int, ...] = (200,),
) -> tuple[bytes, ManifestRecord | None, bool]:
    """Fetch raw bytes once into the content-addressed cache; a non-200 allowed status returns b"" uncached."""
    desc = f"GET {url} {json.dumps(params, sort_keys=True)}"
    sig = hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return cache.read_bytes(cached), cached, True
        except FileNotFoundError:
            pass  # evicted blob: refetch
    resp = _get_with_retries(url, params, timeout_s=timeout_s, attempts=attempts, sleep=sleep, ok_statuses=ok_statuses)
    if resp.status_code != 200:
        return b"", None, False
    payload = resp.content
    man = cache.put(
        archive=archive,
        release=release,
        request_signature=sig,
        query_text=desc,
        payload=payload,
        row_count=0,
        offline_fixture=False,
    )
    return payload, man, False


# --------------------------------------------------------------------------------------------- ZTF


def _parse_ztf_csv(text: str) -> list[dict[str, Any]]:
    if not text.strip() or text.lstrip().startswith("<"):
        return []
    return list(csv.DictReader(io.StringIO(text)))


def _ztf_time_pages(now_mjd: float) -> list[tuple[float, float]]:
    pages: list[tuple[float, float]] = []
    lo = ZTF_FIRST_MJD
    while lo < now_mjd:
        hi = min(now_mjd, lo + 365.25)
        pages.append((lo, hi))
        lo = hi
    return pages


def fetch_ztf_lightcurve(
    *,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
    cache: RawCache,
    now_mjd: float | None = None,
    timeout_s: float = 60.0,
    sleep: Sleeper = time.sleep,
) -> tuple[dict[str, Any] | None, list[ManifestRecord]]:
    """Fetch IRSA ZTF light curves robustly: per band, retried, shrinking cone, yearly TIME pages on timeout.

    Rows are kept only when `catflags == 0` and the epoch position lies within 1.5" of the Gaia position
    propagated to that epoch, so a moving star is matched across ZTF object IDs.
    """
    if now_mjd is None:
        now_mjd = 40587.0 + time.time() / 86400.0
    c_ra, c_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2021.5)
    manifests: list[ManifestRecord] = []
    log: list[str] = []
    rows: list[dict[str, Any]] = []
    any_success = False

    for band in ("g", "r"):
        band_rows: list[dict[str, Any]] | None = None
        for radius in ZTF_CONE_LADDER_ARCSEC:
            base = {"POS": f"CIRCLE {c_ra:.6f} {c_dec:+.6f} {radius / 3600.0:.6f}", "BANDNAME": band, "FORMAT": "CSV"}
            try:
                payload, man, _ = _cached_get(
                    url=IRSA_ZTF_URL, params=base, archive="ztf", release="live", cache=cache,
                    timeout_s=timeout_s, attempts=3, sleep=sleep,
                )
                if man is not None:
                    manifests.append(man)
                band_rows = _parse_ztf_csv(payload.decode("utf-8", errors="replace"))
                log.append(f"{band}: full-range {radius:.1f}\" ok")
                any_success = True
                break
            except ArchiveAccessError as exc:
                log.append(f"{band}: full-range {radius:.1f}\" failed ({exc})")
            # Full-range request failed: page by year at this radius
            paged: list[dict[str, Any]] = []
            n_failed = 0
            pages = _ztf_time_pages(now_mjd)
            for lo, hi in pages:
                params = dict(base, TIME=f"{lo:.1f} {hi:.1f}")
                try:
                    payload, man, _ = _cached_get(
                        url=IRSA_ZTF_URL, params=params, archive="ztf", release="live", cache=cache,
                        timeout_s=timeout_s, attempts=3, sleep=sleep,
                    )
                    if man is not None:
                        manifests.append(man)
                    paged.extend(_parse_ztf_csv(payload.decode("utf-8", errors="replace")))
                except ArchiveAccessError:
                    n_failed += 1
            log.append(f"{band}: paged {radius:.1f}\" {len(pages) - n_failed}/{len(pages)} pages ok")
            if n_failed < len(pages):
                band_rows = paged
                any_success = True
                break
        if band_rows:
            rows.extend(band_rows)

    kept: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for r in rows:
        try:
            mjd = float(r["mjd"])
            mag = float(r["mag"])
            magerr = float(r["magerr"])
            catflags = int(float(r.get("catflags") or 0))
            r_ra = float(r["ra"])
            r_dec = float(r["dec"])
        except (KeyError, TypeError, ValueError):
            continue
        if catflags != 0 or not (0.0 < magerr < 0.2):
            continue
        key = (str(r.get("oid")), f"{mjd:.6f}")
        if key in seen:
            continue
        p_ra, p_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, _mjd_to_jyear(mjd))
        if angular_separation_arcsec(r_ra, r_dec, p_ra, p_dec) > ZTF_MATCH_RADIUS_ARCSEC:
            continue
        seen.add(key)
        kept.append({"oid": str(r.get("oid")), "band": str(r.get("filtercode") or ""), "mjd": mjd, "mag": mag, "magerr": magerr})

    if not any_success:
        return {"survey": "ZTF", "n_epochs": 0, "query_log": log, "status": "query_failed"}, manifests
    kept.sort(key=lambda x: x["mjd"])
    summary = summarize_survey(kept, survey="ZTF", season="calendar", group_key="oid")
    summary["query_log"] = log
    return summary, manifests


# ------------------------------------------------------------------------------------------ ASAS-SN


# The pyasassn client sets no socket timeout. On 2026-11-17 a connection the server had already closed kept a live
# search waiting for 27 minutes, until it was stopped by hand. This is an operational limit, not a science threshold.
ASASSN_CLIENT_TIMEOUT_S = 180.0


def _call_with_deadline(fn: Callable[[], Any], timeout_s: float, label: str) -> Any:
    """Run `fn` and return its result, or raise TimeoutError when it has not finished after `timeout_s` seconds.

    For third-party clients that offer no timeout. The call runs in a daemon thread: a call that never returns is
    abandoned (it cannot be interrupted) and does not keep the process alive at exit.
    """
    import threading

    box: dict[str, Any] = {}

    def _run() -> None:
        try:
            box["value"] = fn()
        except BaseException as exc:  # handed to the caller, which decides what a failure means
            box["error"] = exc

    worker = threading.Thread(target=_run, name=f"deadline:{label}", daemon=True)
    worker.start()
    worker.join(timeout_s)
    if worker.is_alive():
        raise TimeoutError(f"{label} gave no answer within {timeout_s:.0f} s")
    if "error" in box:
        raise box["error"]
    return box["value"]


def fetch_asassn_lightcurve(
    *,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
    cache: RawCache,
    client_factory: Callable[[], Any] | None = None,
    timeout_s: float | None = None,
) -> tuple[dict[str, Any] | None, list[ManifestRecord]]:
    """Fetch the ASAS-SN Sky Patrol v2 light curve of the nearest `stellar_main` source.

    The client has no timeout of its own, so the call runs under `_call_with_deadline`; no answer within
    `timeout_s` (default `ASASSN_CLIENT_TIMEOUT_S`) is a recorded `query_failed`, not a stalled run.

    Uses the public `pyasassn` client (package `skypatrol`). Returns None when the client is not installed.
    Keeps `quality == 'G'` with 0 < mag_err < 0.05, matching the adversarial review's selection.
    """
    c_ra, c_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2018.0)
    desc = f"ASASSN_SKYPATROL stellar_main cone({c_ra:.6f}, {c_dec:+.6f}, {ASASSN_CONE_ARCSEC:.1f} arcsec)"
    sig = hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]
    manifests: list[ManifestRecord] = []
    cached = cache.lookup("asassn", "skypatrol_v2", sig)
    rows: list[dict[str, Any]] | None = None
    if cached is not None and not cached.offline_fixture:
        try:
            rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
            manifests.append(cached)
        except FileNotFoundError:
            rows = None
    if rows is None:
        if client_factory is None:
            try:
                from pyasassn.client import SkyPatrolClient  # type: ignore[import-not-found]
            except ImportError:
                return None, manifests
            client_factory = SkyPatrolClient
        def _query() -> Any:
            client = client_factory()
            return client.cone_search(
                ra_deg=c_ra, dec_deg=c_dec, radius=ASASSN_CONE_ARCSEC / 3600.0, units="deg",
                catalog="stellar_main", download=True, threads=1,
            )

        try:
            res = _call_with_deadline(_query, timeout_s if timeout_s is not None else ASASSN_CLIENT_TIMEOUT_S, "ASAS-SN Sky Patrol")
        except Exception as exc:  # network/client failure, or no answer within the deadline
            return {"survey": "ASAS-SN", "n_epochs": 0, "status": "query_failed", "query_log": [f"{type(exc).__name__}: {exc}"]}, manifests
        cat = getattr(res, "catalog_info", None)
        data = getattr(res, "data", None)
        if cat is None or data is None or len(cat) == 0:
            rows = []
        else:
            cat_rows = cat.to_dict("records")
            best = min(cat_rows, key=lambda c: angular_separation_arcsec(c_ra, c_dec, float(c["ra_deg"]), float(c["dec_deg"])))
            sep = angular_separation_arcsec(c_ra, c_dec, float(best["ra_deg"]), float(best["dec_deg"]))
            rows = []
            if sep <= ASASSN_MATCH_RADIUS_ARCSEC:
                sub = data[data["asas_sn_id"] == best["asas_sn_id"]]
                cam_col = "camera" if "camera" in sub.columns else ("cam" if "cam" in sub.columns else None)
                for rec in sub.to_dict("records"):
                    rows.append(
                        {
                            "jd": float(rec["jd"]),
                            "mag": None if rec.get("mag") is None else float(rec["mag"]),
                            "mag_err": None if rec.get("mag_err") is None else float(rec["mag_err"]),
                            "quality": str(rec.get("quality")),
                            "camera": str(rec.get(cam_col)) if cam_col else "",
                            "band": str(rec.get("phot_filter")),
                            "asas_sn_id": str(best["asas_sn_id"]),
                            "sep_arcsec": round(sep, 2),
                        }
                    )
        payload = json.dumps(rows, sort_keys=True).encode("utf-8")
        manifests.append(
            cache.put(archive="asassn", release="skypatrol_v2", request_signature=sig, query_text=desc,
                      payload=payload, row_count=len(rows), offline_fixture=False)
        )

    kept = []
    for r in rows:
        mag, err = r.get("mag"), r.get("mag_err")
        if r.get("quality") != "G" or mag is None or err is None or not (0.0 < err < 0.05) or math.isnan(mag):
            continue
        kept.append({"oid": r.get("asas_sn_id", ""), "camera": r.get("camera", ""), "band": r.get("band", ""),
                     "mjd": float(r["jd"]) - 2400000.5, "mag": float(mag), "magerr": float(err)})
    kept.sort(key=lambda x: x["mjd"])
    return summarize_survey(kept, survey="ASAS-SN", season="july_june", group_key="camera"), manifests


# ---------------------------------------------------------------------------------- survey summary


def summarize_survey(rows: list[dict[str, Any]], *, survey: str, season: str, group_key: str) -> dict[str, Any]:
    """Seasonal nightly-median photometry per band.

    `group_key='oid'` keeps the dominant ZTF object ID per band (field/CCD zero points differ);
    `group_key='camera'` subtracts each ASAS-SN camera's median before combining cameras.
    Seasons are calendar years (`season='calendar'`) or July-June (`season='july_june'`).

    Every season carries the robust error of its median, estimated from the scatter of the nightly medians
    (`1.2533 * 1.4826 * MAD / sqrt(n_nights)`), plus its first/median/last night MJD so vetting can test
    coverage of an event window. Each band also keeps its compact nightly-median series (`nights`, [mjd, mag])
    for simultaneous-coverage tests. A successful query that keeps no rows is `status='no_coverage'`, which
    the vetting checks distinguish from `query_failed`.
    """
    out: dict[str, Any] = {
        "survey": survey,
        "n_epochs": len(rows),
        "status": "ok" if rows else "no_coverage",
        "season_kind": season,
        "bands": {},
    }
    if rows:
        out["mjd_range"] = [rows[0]["mjd"], rows[-1]["mjd"]]
    for band in sorted({r["band"] for r in rows}):
        b_rows = [r for r in rows if r["band"] == band]
        if group_key == "oid":
            counts: dict[str, int] = {}
            for r in b_rows:
                counts[r["oid"]] = counts.get(r["oid"], 0) + 1
            main = max(counts, key=lambda k: counts[k])
            sel = [r for r in b_rows if r["oid"] == main]
            mags = np.asarray([r["mag"] for r in sel])
            offsets = {main: 0.0}
        else:
            sel = b_rows
            offsets = {}
            for cam in {r[group_key] for r in sel}:
                offsets[cam] = float(np.median([r["mag"] for r in sel if r[group_key] == cam]))
            mags = np.asarray([r["mag"] - offsets[r[group_key]] for r in sel])
        nights: dict[int, list[float]] = {}
        night_mjds: dict[int, list[float]] = {}
        night_season: dict[int, int] = {}
        for r, m in zip(sel, mags):
            n = int(math.floor(r["mjd"]))
            nights.setdefault(n, []).append(float(m))
            night_mjds.setdefault(n, []).append(float(r["mjd"]))
            yr = _mjd_to_jyear(r["mjd"])
            year = int(math.floor(yr))
            if season == "july_june" and (yr - year) < 0.5:
                year -= 1
            night_season[n] = year
        seasons: dict[int, list[tuple[float, float]]] = {}
        for n, vals in nights.items():
            seasons.setdefault(night_season[n], []).append((float(np.median(night_mjds[n])), float(np.median(vals))))
        table = []
        for yr, v in sorted(seasons.items()):
            s_mjd = np.asarray([x[0] for x in v])
            s_mag = np.asarray([x[1] for x in v])
            med = float(np.median(s_mag))
            err = None
            if len(s_mag) >= 2:
                err = 1.2533 * 1.4826 * float(np.median(np.abs(s_mag - med))) / math.sqrt(len(s_mag))
            table.append(
                {
                    "season": yr,
                    "n_nights": len(v),
                    "median_mag": round(med, 4),
                    "median_err_mag": None if err is None else round(err, 5),
                    "mjd_first": round(float(s_mjd.min()), 3),
                    "mjd_mid": round(float(np.median(s_mjd)), 3),
                    "mjd_last": round(float(s_mjd.max()), 3),
                }
            )
        nightly = sorted((round(float(np.median(night_mjds[n])), 4), round(float(np.median(vals)), 4)) for n, vals in nights.items())
        good = [s["median_mag"] for s in table if s["n_nights"] >= 5]
        rng = float(max(good) - min(good)) if len(good) >= 2 else None
        band_summary: dict[str, Any] = {
            "n_epochs": len(sel),
            "n_nights": len(nights),
            "seasons": table,
            "seasonal_range_mag": None if rng is None else round(rng, 4),
            "robust_sigma_mag": round(float(1.4826 * np.median(np.abs(mags - np.median(mags)))), 4) if len(mags) else None,
            "nightly_robust_sigma_mag": (
                round(float(1.4826 * np.median(np.abs(np.asarray([x[1] for x in nightly]) - np.median([x[1] for x in nightly])))), 4)
                if nightly else None
            ),
            "nights": [list(x) for x in nightly],
        }
        if group_key == "oid":
            band_summary["object_id"] = main
            band_summary["median_mag"] = round(float(np.median(mags)), 4)
        out["bands"][band] = band_summary
    return out


# --------------------------------------------------------------------------------------------- TESS


def _bin_series(t: np.ndarray, f: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray]:
    if len(t) == 0:
        return t, f
    idx = np.floor((t - t[0]) / width).astype(int)
    uniq, inv = np.unique(idx, return_inverse=True)
    sums = np.bincount(inv, weights=f)
    tsum = np.bincount(inv, weights=t)
    cnt = np.bincount(inv)
    return tsum / cnt, sums / cnt


def lookup_tic_id(*, ra: float, dec: float, gaia_id: str, cache: RawCache, sleep: Sleeper = time.sleep) -> tuple[str | None, list[ManifestRecord]]:
    """Resolve the TIC ID via MAST TIC position search, preferring the row whose GAIA id matches."""
    req = {
        "service": "Mast.Catalogs.Filtered.Tic.Position",
        "params": {"ra": round(ra, 6), "dec": round(dec, 6), "radius": 0.003, "columns": "ID,ra,dec,Tmag,GAIA", "filters": []},
        "format": "json",
    }
    payload, man, _ = _cached_get(
        url=MAST_INVOKE_URL, params={"request": json.dumps(req, sort_keys=True)}, archive="tess", release="tic",
        cache=cache, timeout_s=30.0, sleep=sleep,
    )
    mans = [man] if man is not None else []
    rows = json.loads(payload.decode("utf-8") or "{}").get("data", [])
    if not rows:
        return None, mans
    for r in rows:
        if str(r.get("GAIA") or "") == str(gaia_id):
            return str(r["ID"]), mans
    best = min(rows, key=lambda r: angular_separation_arcsec(ra, dec, float(r["ra"]), float(r["dec"])))
    return str(best["ID"]), mans


def list_tess_sectors(*, ra: float, dec: float, cache: RawCache, sleep: Sleeper = time.sleep) -> tuple[list[int], list[ManifestRecord]]:
    payload, man, _ = _cached_get(
        url=TESSCUT_SECTOR_URL, params={"ra": f"{ra:.6f}", "dec": f"{dec:.6f}"}, archive="tess", release="tesscut",
        cache=cache, timeout_s=30.0, sleep=sleep,
    )
    js = json.loads(payload.decode("utf-8") or "{}")
    sectors = sorted({int(r["sector"]) for r in js.get("results", [])})
    return sectors, ([man] if man is not None else [])


def qlp_url(tic_id: str, sector: int) -> str:
    t16 = f"{int(tic_id):016d}"
    return (
        f"{QLP_S3_BASE}/s{sector:04d}/{t16[0:4]}/{t16[4:8]}/{t16[8:12]}/{t16[12:16]}/"
        f"hlsp_qlp_tess_ffi_s{sector:04d}-{t16}_tess_v01_llc.fits"
    )


def read_qlp_raw(fits_bytes: bytes) -> tuple[np.ndarray, np.ndarray]:
    """Return (BTJD, SAP_FLUX) of good-quality cadences from a QLP light-curve FITS file (undetrended)."""
    from astropy.io import fits

    with fits.open(io.BytesIO(fits_bytes)) as hdul:
        d = hdul[1].data
        t = np.asarray(d["TIME"], dtype=float)
        f = np.asarray(d["SAP_FLUX"], dtype=float)
        q = np.asarray(d["QUALITY"], dtype=int)
    ok = np.isfinite(t) & np.isfinite(f) & (q == 0) & (f > 0)
    return t[ok], f[ok]


def extract_tesscut_photometry(
    zip_bytes: bytes,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
    neighbours: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """Aperture photometry plus the per-cadence systematics series of one TESScut cube.

    Flux: 1.6 px circular aperture at the PM-propagated position minus the 4-6 px annulus median times the
    aperture area. Alongside it the per-cadence annulus background (per pixel), the background-subtracted
    flux-weighted centroid inside 2.5 px (offsets from the expected target pixel, in pixels), the expected target
    pixel itself, and the pixel offsets of any supplied neighbours (Gaia positions) are returned so vetting can
    test scattered light, pointing jumps and difference-image centroid shifts on the same cadences.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    empty = {"time": np.asarray([]), "flux": np.asarray([]), "bkg_per_px": np.asarray([]), "centroid_dcol_px": np.asarray([]),
             "centroid_drow_px": np.asarray([]), "n_ap_px": 0, "target_px": None, "neighbour_px": []}
    with zipfile.ZipFile(io.BytesIO(zip_bytes)) as zf:
        name = next(n for n in zf.namelist() if n.endswith(".fits"))
        raw = zf.read(name)
    with fits.open(io.BytesIO(raw)) as hdul:
        d = hdul[1].data
        t = np.asarray(d["TIME"], dtype=float)
        cube = np.asarray(d["FLUX"], dtype=float)
        q = np.asarray(d["QUALITY"], dtype=int)
        wcs = WCS(hdul[2].header)
    t_ok = t[np.isfinite(t)]
    mid_jyear = 2000.0 + ((float(np.median(t_ok)) if len(t_ok) else 1500.0) + BTJD_OFFSET - 2451545.0) / 365.25
    p_ra, p_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, mid_jyear)
    x0, y0 = (float(v) for v in wcs.world_to_pixel_values(p_ra, p_dec))
    ny, nx = cube.shape[1], cube.shape[2]
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = np.hypot(xx - x0, yy - y0)
    ap = r <= TESS_APERTURE_RADIUS_PX
    ann = (r >= TESS_ANNULUS_PX[0]) & (r <= TESS_ANNULUS_PX[1])
    cen = r <= TESS_CENTROID_RADIUS_PX
    if not ap.any() or not ann.any():
        return empty
    bkg = np.nanmedian(cube[:, ann], axis=1)
    flux = np.nansum(cube[:, ap], axis=1) - bkg * float(ap.sum())
    sub = np.clip(np.nan_to_num(cube[:, cen], nan=0.0) - bkg[:, None], 0.0, None)
    w_sum = np.sum(sub, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        cx = np.sum(sub * xx[cen][None, :], axis=1) / w_sum - x0
        cy = np.sum(sub * yy[cen][None, :], axis=1) / w_sum - y0
    ok = np.isfinite(t) & np.isfinite(flux) & (q == 0) & (flux > 0) & np.isfinite(bkg) & np.isfinite(cx) & np.isfinite(cy)
    nb_px: list[dict[str, Any]] = []
    for nb in neighbours or []:
        try:
            n_ra, n_dec = propagate_position(float(nb["ra_deg"]), float(nb["dec_deg"]), float(nb.get("pmra") or 0.0),
                                             float(nb.get("pmdec") or 0.0), 2016.0, mid_jyear)
        except (KeyError, TypeError, ValueError):
            continue
        nxp, nyp = (float(v) for v in wcs.world_to_pixel_values(n_ra, n_dec))
        nb_px.append({"id": nb.get("id"), "dcol_px": round(nxp - x0, 3), "drow_px": round(nyp - y0, 3), "flux_ratio": nb.get("flux_ratio")})
    return {
        "time": t[ok],
        "flux": flux[ok],
        "bkg_per_px": bkg[ok],
        "centroid_dcol_px": cx[ok],
        "centroid_drow_px": cy[ok],
        "n_ap_px": int(ap.sum()),
        "target_px": [round(x0, 3), round(y0, 3)],
        "neighbour_px": nb_px,
    }


def extract_tesscut_aperture(zip_bytes: bytes, ra: float, dec: float, pmra: float, pmdec: float) -> tuple[np.ndarray, np.ndarray]:
    """Simple-aperture photometry on a TESScut cube: 1.6 px circular aperture minus a 4-6 px annulus median."""
    ph = extract_tesscut_photometry(zip_bytes, ra, dec, pmra, pmdec)
    return ph["time"], ph["flux"]


def match_cadence_index(t1: np.ndarray, t2: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (mask over t1, index into t2) of nearest-neighbour cadence pairs within half the coarser cadence."""
    if len(t1) < 2 or len(t2) < 2:
        return np.zeros(len(t1), dtype=bool), np.asarray([], dtype=int)
    tol = 0.5 * max(float(np.median(np.diff(t1))), float(np.median(np.diff(t2))))
    order = np.argsort(t2)
    t2s = t2[order]
    j = np.clip(np.searchsorted(t2s, t1), 1, len(t2s) - 1)
    left, right = t2s[j - 1], t2s[j]
    jj = np.where(np.abs(t1 - left) <= np.abs(right - t1), j - 1, j)
    keep = np.abs(t2s[jj] - t1) <= tol
    return keep, order[jj[keep]]


def match_cadences(t1: np.ndarray, f1: np.ndarray, t2: np.ndarray, f2: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Keep only cadences present in both series (nearest neighbour within half the coarser cadence)."""
    if len(t1) < 2 or len(t2) < 2:
        return np.asarray([]), np.asarray([]), np.asarray([])
    keep, idx2 = match_cadence_index(t1, t2)
    return t1[keep], f1[keep], f2[idx2]


def _bin_with_counts(t: np.ndarray, width: float) -> tuple[np.ndarray, np.ndarray]:
    idx = np.floor((t - t[0]) / width).astype(int)
    _, inv = np.unique(idx, return_inverse=True)
    return inv, np.bincount(inv)


def fetch_tess_pipeline_pair(
    *,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
    gaia_id: str,
    cache: RawCache,
    max_sectors: int = 8,
    neighbours: list[dict[str, Any]] | None = None,
    sleep: Sleeper = time.sleep,
) -> tuple[dict[str, Any] | None, list[ManifestRecord]]:
    """Fetch QLP raw SAP flux and an independent TESScut aperture extraction for every sector with both.

    Returns 30-minute binned, cadence-matched, per-sector median-normalised series. Detrending is left to the
    `pipeline_agreement` check so both series receive an identical filter. Each sector also carries the number of
    native cadences per 30-min bin, the native cadence, the TESScut annulus background summed over the aperture as
    a fraction of the median target flux (`bkg_frac`), the flux-weighted centroid offsets (pixels) and the pixel
    offsets of the supplied Gaia neighbours. `status` is `ok`, `no_coverage` (no TIC, no sector, or no QLP light
    curve) or `query_failed` (every sector that exists failed to download or parse).
    """
    manifests: list[ManifestRecord] = []
    log: list[str] = []
    c_ra, c_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2021.0)
    tic, m = lookup_tic_id(ra=c_ra, dec=c_dec, gaia_id=gaia_id, cache=cache, sleep=sleep)
    manifests.extend(m)
    if tic is None:
        return {"tic_id": None, "sectors": [], "status": "no_coverage", "query_log": ["TIC lookup returned no rows"]}, manifests
    sectors, m = list_tess_sectors(ra=c_ra, dec=c_dec, cache=cache, sleep=sleep)
    manifests.extend(m)
    out_sectors: list[dict[str, Any]] = []
    n_failed = 0
    n_attempted = 0
    for sector in sectors:
        if len(out_sectors) >= max_sectors:
            log.append(f"s{sector}: skipped (max_sectors={max_sectors})")
            continue
        n_attempted += 1
        try:
            qb, man, _ = _cached_get(url=qlp_url(tic, sector), params={}, archive="tess", release="qlp", cache=cache,
                                     timeout_s=60.0, sleep=sleep, ok_statuses=(200, 403, 404))
            if man is not None:
                manifests.append(man)
            if not qb:
                log.append(f"s{sector}: no QLP light curve")
                continue
            t_q, f_q = read_qlp_raw(qb)
            zb, man, _ = _cached_get(
                url=TESSCUT_ASTROCUT_URL,
                params={"ra": f"{c_ra:.6f}", "dec": f"{c_dec:.6f}", "x": TESSCUT_SIZE_PX, "y": TESSCUT_SIZE_PX, "units": "px", "sector": sector},
                archive="tess", release="tesscut", cache=cache, timeout_s=180.0, attempts=2, sleep=sleep,
            )
            if man is not None:
                manifests.append(man)
            ph = extract_tesscut_photometry(zb, ra, dec, pmra, pmdec, neighbours=neighbours)
        except (ArchiveAccessError, KeyError, OSError, ValueError, StopIteration, zipfile.BadZipFile) as exc:
            n_failed += 1
            log.append(f"s{sector}: query_failed {type(exc).__name__}: {str(exc)[:160]}")
            continue
        t_c = ph["time"]
        if len(t_q) < 2 or len(t_c) < 2:
            log.append(f"s{sector}: only {min(len(t_q), len(t_c))} good cadences")
            continue
        keep, idx_c = match_cadence_index(t_q, t_c)
        t_m, fq_m = t_q[keep], f_q[keep]
        if len(t_m) < 50:
            log.append(f"s{sector}: only {len(t_m)} matched cadences")
            continue
        fc_m = ph["flux"][idx_c]
        inv, cnt = _bin_with_counts(t_m, TESS_BIN_DAYS)

        def _b(v: np.ndarray) -> np.ndarray:
            return np.bincount(inv, weights=v) / cnt

        tb, fqb, fcb = _b(t_m), _b(fq_m), _b(fc_m)
        med_c = float(np.median(fc_m))
        bkg_frac = _b(ph["bkg_per_px"][idx_c]) * float(ph["n_ap_px"]) / max(1e-9, med_c)
        out_sectors.append(
            {
                "sector": sector,
                "time_btjd": [round(float(x), 5) for x in tb],
                "flux_qlp_raw": [round(float(x), 6) for x in fqb / np.median(fqb)],
                "flux_tesscut": [round(float(x), 6) for x in fcb / np.median(fcb)],
                "n_per_bin": [int(x) for x in cnt],
                "cadence_days": round(float(np.median(np.diff(t_m))), 6),
                "bkg_frac": [round(float(x), 6) for x in bkg_frac],
                "centroid_dcol_px": [round(float(x), 5) for x in _b(ph["centroid_dcol_px"][idx_c])],
                "centroid_drow_px": [round(float(x), 5) for x in _b(ph["centroid_drow_px"][idx_c])],
                "target_px": ph["target_px"],
                "neighbour_px": ph["neighbour_px"],
                "n_matched_cadences": int(len(t_m)),
            }
        )
        log.append(f"s{sector}: ok ({len(t_m)} matched cadences, {len(tb)} 30-min bins)")
    if out_sectors:
        status = "ok"
    elif n_attempted and n_failed == n_attempted:
        status = "query_failed"
    else:
        status = "no_coverage"
    if not sectors:
        log.append("TESScut sector list empty")
    return {"tic_id": tic, "sectors": out_sectors, "status": status, "query_log": log}, manifests


# ------------------------------------------------------------------------------ NEOWISE comparison stars


def build_neowise_comparison(
    rows: list[dict[str, Any]],
    *,
    ra: float,
    dec: float,
    visit_windows: list[tuple[float, float]],
    target_w1_mag: float,
    exclude_arcsec: float,
    mag_window: float,
    min_stars: int,
    max_stars: int,
    link_arcsec: float = 1.5,
) -> dict[str, Any]:
    """Field-star ensemble on the target's NEOWISE visit grid from one `neowiser_p1bs_psd` cone.

    Detections (`qual_frame > 0`) are grouped into stars by friends-of-friends linking within `link_arcsec` on the
    tangent plane. Stars within `exclude_arcsec` of the target, outside `|W1 - W1_target| <= mag_window`, or
    detected in fewer than half of the target's visits are dropped; the `max_stars` closest in W1 are kept. For
    each kept star and visit window the median W1/W2 is taken and the star's own median over visits removed. The
    ensemble is the per-visit median of those offsets (None when fewer than `min_stars` stars cover a visit).
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    from scipy.spatial import cKDTree

    def _f(v: Any) -> float | None:
        try:
            x = float(v)
        except (TypeError, ValueError):
            return None
        return None if (math.isnan(x) or x <= -900.0) else x

    det = []
    cos_d = math.cos(math.radians(dec))
    for r in rows:
        r_ra, r_dec, mjd, w1, w2 = (_f(r.get(k)) for k in ("ra", "dec", "mjd", "w1mpro", "w2mpro"))
        qf = _f(r.get("qual_frame"))
        if None in (r_ra, r_dec, mjd, w1, w2) or (qf is not None and qf <= 0):
            continue
        dra = ((r_ra - ra + 540.0) % 360.0 - 180.0) * cos_d * 3600.0  # type: ignore[operator]
        det.append((dra, (r_dec - dec) * 3600.0, mjd, w1, w2))  # type: ignore[operator]
    out: dict[str, Any] = {"n_detections": len(det), "n_stars": 0, "status": "too_few_comparison_stars", "stars": []}
    if len(det) < 2:
        return out
    arr = np.asarray(det, dtype=float)
    tree = cKDTree(arr[:, :2])
    pairs = tree.query_pairs(link_arcsec, output_type="ndarray")
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(arr), len(arr))) if len(pairs) else coo_matrix((len(arr), len(arr)))
    _, labels = connected_components(graph, directed=False)
    n_vis = len(visit_windows)
    cands = []
    for lab in np.unique(labels):
        g = arr[labels == lab]
        x, y = float(np.median(g[:, 0])), float(np.median(g[:, 1]))
        if math.hypot(x, y) <= exclude_arcsec:
            continue
        w1_med = float(np.median(g[:, 3]))
        if abs(w1_med - target_w1_mag) > mag_window:
            continue
        v1, v2 = [], []
        for lo, hi in visit_windows:
            sel = (g[:, 2] >= lo - 2.0) & (g[:, 2] <= hi + 2.0)
            v1.append(float(np.median(g[sel, 3])) if sel.any() else None)
            v2.append(float(np.median(g[sel, 4])) if sel.any() else None)
        n_cov = sum(v is not None for v in v2)
        if n_cov < max(2, 0.5 * n_vis):
            continue
        cands.append((abs(w1_med - target_w1_mag), x, y, w1_med, v1, v2))
    cands.sort(key=lambda c: c[0])
    kept = cands[:max_stars]
    out["n_stars"] = len(kept)
    if not kept:
        return out
    d1 = np.full((len(kept), n_vis), np.nan)
    d2 = np.full((len(kept), n_vis), np.nan)
    for i, (_, x, y, w1_med, v1, v2) in enumerate(kept):
        a1 = np.asarray([np.nan if v is None else v for v in v1])
        a2 = np.asarray([np.nan if v is None else v for v in v2])
        d1[i] = a1 - np.nanmedian(a1)
        d2[i] = a2 - np.nanmedian(a2)
        out["stars"].append({"dx_arcsec": round(x, 1), "dy_arcsec": round(y, 1), "w1_mag": round(w1_med, 3),
                             "w2_amp_mag": round(float(np.nanmax(d2[i]) - np.nanmin(d2[i])), 4)})
    n_cover = np.sum(np.isfinite(d2), axis=0)
    with np.errstate(all="ignore"):
        e1 = np.nanmedian(d1, axis=0)
        e2 = np.nanmedian(d2, axis=0)
    ok = n_cover >= min_stars
    out["ensemble_dw1_mag"] = [round(float(v), 4) if k else None for v, k in zip(e1, ok)]
    out["ensemble_dw2_mag"] = [round(float(v), 4) if k else None for v, k in zip(e2, ok)]
    errs = []
    for j in range(n_vis):
        col = d2[:, j][np.isfinite(d2[:, j])]
        if ok[j] and len(col) >= 2:
            errs.append(round(max(1e-4, 1.2533 * 1.4826 * float(np.median(np.abs(col - np.median(col)))) / math.sqrt(len(col))), 5))
        else:
            errs.append(None)
    out["ensemble_dw2_err_mag"] = errs
    out["n_stars_per_visit"] = [int(n) for n in n_cover]
    out["status"] = "ok" if len(kept) >= min_stars else "too_few_comparison_stars"
    return out


# ------------------------------------------------------------------------------------------ SkyBoT


def skybot_cone_search(**kwargs: Any) -> tuple[list[dict[str, Any]], ManifestRecord]:
    """Thin wrapper around `oddstar.io.skybot.skybot_cone_search` (imported lazily; tests monkeypatch this)."""
    from oddstar.io.skybot import skybot_cone_search as _impl  # type: ignore[import-not-found]

    return _impl(**kwargs)


def query_sso_at_events(
    events: list[dict[str, Any]],
    *,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
    radius_arcsec: float,
    cache: RawCache,
    sleep: Sleeper = time.sleep,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Query SkyBoT (observer C57 = TESS) at the start, middle and end of every flagged TESS event bin.

    Three epochs per bin keep a body moving up to ~28"/h inside `radius_arcsec` (2 TESS pixels) at one of the
    sampled epochs. Returns one record per event with the bodies found (closest approach kept per name) and its
    query status; any failed epoch marks that event `query_failed`.
    """
    manifests: list[ManifestRecord] = []
    records: list[dict[str, Any]] = []
    for ev in events:
        rec: dict[str, Any] = {"sector": ev.get("sector"), "btjd": ev.get("btjd"), "radius_arcsec": radius_arcsec,
                               "observer": "C57", "status": "ok", "bodies": [], "query_log": []}
        bodies: dict[str, dict[str, Any]] = {}
        for t_b in (float(ev["btjd_start"]), float(ev["btjd"]), float(ev["btjd_end"])):
            epoch_jd = t_b + BTJD_OFFSET
            jyr = 2000.0 + (epoch_jd - 2451545.0) / 365.25
            p_ra, p_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, jyr)
            try:
                found, man = skybot_cone_search(ra_deg=p_ra, dec_deg=p_dec, radius_arcsec=radius_arcsec, epoch_jd=epoch_jd,
                                                observer="C57", cache=cache, sleep=sleep)
            except ImportError as exc:
                rec["status"] = "client_unavailable"
                rec["query_log"].append(f"JD {epoch_jd:.4f}: {type(exc).__name__}: {str(exc)[:160]}")
                break
            except Exception as exc:  # noqa: BLE001 - the SkyBoT client contract is "raises on failure"; record any failure
                rec["status"] = "query_failed"
                rec["query_log"].append(f"JD {epoch_jd:.4f}: {type(exc).__name__}: {str(exc)[:160]}")
                continue
            manifests.append(man)
            rec["query_log"].append(f"JD {epoch_jd:.4f}: ok ({len(found)} bodies)")
            for b in found:
                name = str(b.get("name"))
                prev = bodies.get(name)
                if prev is None or float(b.get("sep_arcsec") or 1e9) < float(prev.get("sep_arcsec") or 1e9):
                    bodies[name] = dict(b, epoch_jd=round(epoch_jd, 5))
        rec["bodies"] = sorted(bodies.values(), key=lambda b: float(b.get("sep_arcsec") or 1e9))
        records.append(rec)
    return {"events": records}, manifests
