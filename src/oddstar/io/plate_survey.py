"""Historical plate sources around every star of a live sample (survey mode of the vanishing search).

The vanishing channel asks whether a source recorded on mid-century photographic plates is still in the sky. A
sample star is in Gaia by construction, so it cannot be the unit of that search. The unit is a plate source: every
USNO-B1.0 object within `plate_field_radius_arcsec` of a sample position that was detected on both first-epoch
plates (blue and red) and is not flagged as lying on a diffraction spike. For each field three catalogues are read:

1. USNO-B1.0 (VizieR I/284), the plate sources. A position is given at epoch 2000 when the catalogue measured a
   proper motion, and at the mean epoch of the plates when it did not (proper motion listed as zero).
2. Gaia DR3, moved to that epoch with its proper motions. A Gaia source within `match_radius_arcsec` of a plate
   source is its modern counterpart. The catalogue's own proper motion can be wrong (mismatched detections), which
   puts its epoch-2000 position far from a star that is still there. So a Gaia source also counts when its path
   over the plate era (`PLATE_ERA`), drawn with the Gaia proper motion, passes within `match_radius_arcsec` of the
   plate source at its mean plate epoch.
3. Pan-STARRS1 DR2 mean objects with at least `plate_ps1_min_detections` detections. An object within
   `match_radius_arcsec` is a counterpart. The local depth is the `plate_depth_percentile` of the r-band PSF
   magnitudes of the field's objects, and needs `plate_depth_min_objects` of them.

A plate source is **searched** only when all three answers were read, the field has a Pan-STARRS depth, and the
source is bright enough that its absence would be significant (depth minus first red magnitude of at least
`min_depth_margin_mag`). Everything else is counted, with the reason, and excluded from the null result:

- "no_coverage": the field has no Pan-STARRS depth (south of Dec -30 deg, or too few objects);
- "failed": the Gaia or Pan-STARRS query of the field failed (not cached, retried next run);
- "too_faint": the source is closer than `min_depth_margin_mag` to the local depth;
- "near_bright_star": the source lies inside the halo radius of a bright Gaia star (`plate_bright_star_halo`),
  where plates show halos, spikes and ghosts that both first-epoch plates repeat, and is not that star itself.
A field whose USNO-B query failed yields no sources at all; it is listed in the summary (`failed_field_ids`).

Limits the run must carry: USNO-B contains false sources (plate flaws, halos of bright stars, mismatched
detections); the two first-epoch plates were usually taken the same night, so they are independent emulsions but
not independent epochs; no plate-epoch Solar System ephemeris, high-resolution imaging or variable-star lookup is
made here, so those vetting checks stay `not_evaluated`.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Callable
from typing import Any

import numpy as np
import requests

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io import live_survey as ls
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec, propagate_position

USNO_TABLE = "I/284/out"
USNO_COLUMNS = "USNO-B1.0,RAJ2000,DEJ2000,Epoch,pmRA,pmDE,Ndet,Flags,B1mag,R1mag,B2mag,R2mag,Imag"
# A position where USNO-B1.0 certainly has sources (COSMOS field centre). A VizieR mirror that does not hold a table
# answers with zero rows for every cone, so an empty field is accepted from a mirror only after it returned rows here.
USNO_CANARY_CONE = (150.119000, 2.205800, 40.0)
USNO_ROW_LIMIT = 2000
USNO_POSITION_EPOCH = 2000.0  # RAJ2000/DEJ2000 of I/284 are at epoch 2000 when a proper motion was applied (note 2)
GAIA_EPOCH = 2016.0
PS1_COLUMNS = "[objID,raMean,decMean,nDetections,rMeanPSFMag]"
PS1_ROW_LIMIT = 20000
PS1_DEC_LIMIT_DEG = -30.0  # the 3-pi survey footprint
# Gaia sources are read at least this far beyond the field edge, so a counterpart that moved across the edge is still
# found (the largest halo radius of `plate_bright_star_halo` when that is larger)
GAIA_EDGE_MARGIN_ARCSEC = 30.0
# First and last plate epochs of USNO-B1.0 (POSS-I began in 1949, the last southern plates were taken in 2002)
PLATE_ERA = (1949.0, 2002.0)

ProgressCallback = Callable[[int, int, str, str], None]
_usno_mirror_serves: dict[str, bool] = {}


def _settings(cfg: PipelineConfig) -> dict[str, Any]:
    ch = cfg.channels["vanishing"]
    return {
        "radius": float(ch["plate_field_radius_arcsec"]), "match": float(ch["match_radius_arcsec"]),
        "min_margin": float(ch["min_depth_margin_mag"]), "ps1_min_det": int(ch["plate_ps1_min_detections"]),
        "depth_pct": float(ch["plate_depth_percentile"]), "depth_min_n": int(ch["plate_depth_min_objects"]),
        "attempts": int(ch["plate_max_attempts"]), "backoff": float(ch["plate_backoff_s"]), "rps": float(ch["plate_rate_limit_rps"]),
        "halo": sorted((float(g), float(r)) for g, r in ch["plate_bright_star_halo"]),
    }


def gaia_margin_arcsec(st: dict[str, Any]) -> float:
    """How far beyond the field edge Gaia is read: far enough for every halo radius."""
    return max([GAIA_EDGE_MARGIN_ARCSEC, *[r for _, r in st["halo"]]])


def halo_radius_arcsec(g_mag: float | None, halo: list[tuple[float, float]]) -> float | None:
    """Halo radius of a star of magnitude `g_mag`: that of the first limit it is brighter than; None for a faint star."""
    if g_mag is None:
        return None
    return next((r for g_lim, r in halo if g_mag < g_lim), None)


def _float(val: Any) -> float | None:
    try:
        out = float(val)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


# ---------------------------------------------------------------------------------------------------- USNO-B1.0


def _usno_get(base: str, ra: float, dec: float, radius_arcsec: float, timeout_s: float) -> list[dict[str, Any]]:
    resp = requests.get(base, params={"-source": USNO_TABLE, "-c": f"{ra:.6f} {dec:+.6f}", "-c.r": f"{radius_arcsec:.2f}", "-c.u": "arcsec",
                                      "-out": USNO_COLUMNS, "-out.max": str(USNO_ROW_LIMIT)}, timeout=timeout_s)
    body = ls._strip_leaked_cgi_headers(resp.text) if resp.status_code == 200 else ""
    if resp.status_code != 200 or not body.strip().startswith("#"):
        raise ArchiveAccessError(f"HTTP {resp.status_code}" + ("" if resp.status_code != 200 else " with a non-VizieR body"), retryable=True)
    rows = ls._parse_vizier_tsv(body)
    if rows and "USNO-B1.0" not in rows[0]:
        raise ArchiveAccessError(f"answer lacks the USNO-B1.0 column (has {sorted(rows[0])[:6]})", retryable=True)
    if len(rows) >= USNO_ROW_LIMIT:
        raise ArchiveAccessError(f"{len(rows)} rows reach the row limit; the field answer may be cut short", retryable=False)
    return rows


def usno_cone(ra: float, dec: float, radius_arcsec: float, cache: RawCache, *, timeout_s: float = 30.0) -> tuple[list[dict[str, Any]], str, bool]:
    """USNO-B1.0 rows of one field. Returns (rows, manifest id, from_cache); raises ArchiveAccessError (not cached)."""
    desc = f"VIZIER_TSV {USNO_TABLE} [{USNO_COLUMNS}] cone({ra:.6f}, {dec:+.6f}, {radius_arcsec:.2f} arcsec)"
    sig = hashlib.sha256(f"{'|'.join(ls.VIZIER_TSV_URLS)}|{desc}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("usnob1", "I_284", sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(cache.read_bytes(cached).decode("utf-8")), cached.manifest_id, True
        except (FileNotFoundError, OSError, ValueError):
            pass
    errors: list[str] = []
    for base in ls.VIZIER_TSV_URLS:
        try:
            rows = _usno_get(base, ra, dec, radius_arcsec, timeout_s)
            if not rows:
                if base not in _usno_mirror_serves:
                    _usno_mirror_serves[base] = bool(_usno_get(base, *USNO_CANARY_CONE, timeout_s))
                if not _usno_mirror_serves[base]:
                    errors.append(f"{base}: no rows even at the canary position (table not held)")
                    continue
            man = cache.put(archive="usnob1", release="I_284", request_signature=sig, query_text=desc,
                            payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
            return rows, man.manifest_id, False
        except (ArchiveAccessError, requests.exceptions.RequestException) as exc:
            errors.append(f"{base}: {type(exc).__name__}: {str(exc)[:120]}")
    raise ArchiveAccessError(f"USNO-B1.0 cone failed on every VizieR mirror: {'; '.join(errors)}", retryable=True)


def parse_usno_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Typed plate sources. A row without an id or a position is dropped."""
    out = []
    for r in rows:
        ra, dec = _float(r.get("RAJ2000")), _float(r.get("DEJ2000"))
        uid = str(r.get("USNO-B1.0") or "").strip()
        if not uid or ra is None or dec is None:
            continue
        first = [m for m in (_float(r.get("B1mag")), _float(r.get("R1mag"))) if m is not None]
        out.append({
            "usno_id": uid, "ra_deg": ra, "dec_deg": dec, "mean_epoch": _float(r.get("Epoch")),
            "pmra_mas_yr": _float(r.get("pmRA")) or 0.0, "pmdec_mas_yr": _float(r.get("pmDE")) or 0.0,
            "n_detections": int(_float(r.get("Ndet")) or 0), "flags": str(r.get("Flags") or "").strip(),
            "b1_mag": _float(r.get("B1mag")), "r1_mag": _float(r.get("R1mag")), "b2_mag": _float(r.get("B2mag")),
            "r2_mag": _float(r.get("R2mag")), "i_mag": _float(r.get("Imag")), "n_first_epoch_plates": len(first),
            "on_spike": "s" in str(r.get("Flags") or ""),
        })
        u = out[-1]
        # Without a catalogue proper motion nothing was moved to 2000: the position is that of the mean plate epoch
        u["position_epoch"] = USNO_POSITION_EPOCH if (u["pmra_mas_yr"] or u["pmdec_mas_yr"] or u["mean_epoch"] is None) else u["mean_epoch"]
        # Where the plates saw it: the catalogue position taken back to the mean plate epoch with the catalogue's own motion
        u["plate_epoch"] = u["position_epoch"] if u["mean_epoch"] is None else u["mean_epoch"]
        u["plate_ra_deg"], u["plate_dec_deg"] = propagate_position(ra, dec, u["pmra_mas_yr"], u["pmdec_mas_yr"], u["position_epoch"], u["plate_epoch"])
    return out


# ------------------------------------------------------------------------------------------------ Gaia and PS1


def gaia_cone(ra: float, dec: float, radius_arcsec: float, cfg: PipelineConfig, cache: RawCache, *, sleep: Any = None) -> tuple[list[dict[str, Any]], str, bool]:
    q = ("SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag FROM gaiadr3.gaia_source "
         f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_arcsec / 3600.0:.6f}))")
    rows, man, hit = ls._tap_query_cached(urls=ls.GAIA_TAP_URLS, query=q, archive="gaia", release=cfg.gaia.release, cache=cache, timeout_s=30, sleep=sleep)
    out = []
    for r in rows:
        g_ra, g_dec = _float(r.get("ra")), _float(r.get("dec"))
        if g_ra is None or g_dec is None:
            continue
        pmra, pmdec = _float(r.get("pmra")) or 0.0, _float(r.get("pmdec")) or 0.0
        out.append({"source_id": str(r.get("source_id")), "ra_deg": g_ra, "dec_deg": g_dec, "pmra_mas_yr": pmra, "pmdec_mas_yr": pmdec,
                    "g_mag": _float(r.get("phot_g_mean_mag"))})
    return out, man.manifest_id, hit


def ps1_cone(ra: float, dec: float, radius_arcsec: float, min_detections: int, cache: RawCache, *, sleep: ot.Sleeper = time.sleep,
             timeout_s: float = 30.0) -> tuple[list[dict[str, Any]], str, bool]:
    """Pan-STARRS1 DR2 mean objects of one field with at least `min_detections` detections. Failures are not cached."""
    params = {"ra": f"{ra:.6f}", "dec": f"{dec:.6f}", "radius": f"{radius_arcsec / 3600.0:.6f}", "nDetections.gte": str(int(min_detections)),
              "columns": PS1_COLUMNS, "pagesize": str(PS1_ROW_LIMIT)}
    desc = f"GET {ls.MAST_PS1_DR2_URL} {json.dumps(params, sort_keys=True)}"
    sig = hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("deep_optical", "ps1_dr2_field", sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(cache.read_bytes(cached).decode("utf-8")), cached.manifest_id, True
        except (FileNotFoundError, OSError, ValueError):
            pass
    resp = ot._get_with_retries(ls.MAST_PS1_DR2_URL, params, timeout_s=timeout_s, attempts=3, sleep=sleep)
    try:
        js = json.loads(resp.content.decode("utf-8"))
        cols = [c["name"] for c in js["info"]]
        rows = [dict(zip(cols, raw)) for raw in js["data"]]
    except (ValueError, KeyError, TypeError, UnicodeDecodeError) as exc:
        raise ArchiveAccessError(f"Pan-STARRS answer is not the expected JSON: {type(exc).__name__}: {str(exc)[:120]}", retryable=True) from exc
    if rows and not {"raMean", "decMean", "nDetections", "rMeanPSFMag"} <= set(cols):
        raise ArchiveAccessError(f"Pan-STARRS answer lacks requested columns (has {cols})", retryable=True)
    if len(rows) >= PS1_ROW_LIMIT:
        raise ArchiveAccessError(f"{len(rows)} Pan-STARRS rows reach the page size; the field answer may be cut short", retryable=False)
    man = cache.put(archive="deep_optical", release="ps1_dr2_field", request_signature=sig, query_text=desc,
                    payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
    return rows, man.manifest_id, False


def local_depth_mag(ps1_rows: list[dict[str, Any]], percentile: float, min_objects: int) -> tuple[float | None, int]:
    """Local Pan-STARRS depth: a high percentile of the r-band PSF magnitudes of the field's objects. (depth, n used)."""
    mags = [m for m in (_float(r.get("rMeanPSFMag")) for r in ps1_rows) if m is not None and 5.0 < m < 30.0]  # -999 marks no r detection
    if len(mags) < min_objects:
        return None, len(mags)
    return round(float(np.percentile(mags, percentile)), 3), len(mags)


# ------------------------------------------------------------------------------------------------- one field


def _nearest(ra: float, dec: float, cands: list[dict[str, Any]], ra_key: str, dec_key: str) -> tuple[dict[str, Any] | None, float | None]:
    best, best_sep = None, None
    for c in cands:
        if abs(float(c[dec_key]) - dec) * 3600.0 > 120.0:
            continue
        sep = angular_separation_arcsec(ra, dec, float(c[ra_key]), float(c[dec_key]))
        if best_sep is None or sep < best_sep:
            best, best_sep = c, sep
    return best, best_sep


def track_separation(c: dict[str, Any], ra: float, dec: float) -> tuple[float, float]:
    """Smallest separation (arcsec) between a position and a Gaia star's path over `PLATE_ERA`, and the epoch of it."""
    cos_dec = max(1e-6, math.cos(math.radians(dec)))
    x0 = ((c["ra_deg"] - ra + 180.0) % 360.0 - 180.0) * cos_dec * 3600.0  # offset of the star in 2016, arcsec
    y0 = (c["dec_deg"] - dec) * 3600.0
    vx, vy = c["pmra_mas_yr"] / 1000.0, c["pmdec_mas_yr"] / 1000.0  # arcsec per year
    v2 = vx * vx + vy * vy
    epoch = GAIA_EPOCH if v2 == 0.0 else GAIA_EPOCH - (x0 * vx + y0 * vy) / v2
    epoch = min(PLATE_ERA[1], max(PLATE_ERA[0], epoch))
    e_ra, e_dec = propagate_position(c["ra_deg"], c["dec_deg"], c["pmra_mas_yr"], c["pmdec_mas_yr"], GAIA_EPOCH, epoch)
    return angular_separation_arcsec(ra, dec, e_ra, e_dec), epoch


def build_field_units(star: dict[str, Any], usno: list[dict[str, Any]], gaia: list[dict[str, Any]] | None, ps1: list[dict[str, Any]] | None,
                      st: dict[str, Any], *, modern_status: str = "ok") -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Plate-source units of one field, each with its status, and the counts of what was left out.

    `gaia` / `ps1` are None when that catalogue was not read (`modern_status` then says why: "failed: ..." or
    "no_coverage"). A unit's `plate_fetch.status` is "ok", "too_faint", "near_bright_star", "no_coverage" or
    "failed: ...".
    """
    counts = {"n_usno_sources": len(usno), "n_on_spike": 0, "n_one_first_epoch_plate": 0, "n_no_first_epoch_plate": 0}
    sid = str(star["source_id"])
    f_ra, f_dec = float(star["ra_deg"]), float(star["dec_deg"])
    depth, n_depth = (None, 0) if ps1 is None else local_depth_mag(ps1, st["depth_pct"], int(st["depth_min_n"]))
    ps1_real = [] if ps1 is None else [p for p in ps1 if _float(p.get("raMean")) is not None and _float(p.get("decMean")) is not None]
    units: list[dict[str, Any]] = []
    for u in usno:
        if u["on_spike"]:
            counts["n_on_spike"] += 1
            continue
        if u["n_first_epoch_plates"] < 2:
            counts["n_one_first_epoch_plate" if u["n_first_epoch_plates"] == 1 else "n_no_first_epoch_plate"] += 1
            continue
        unit: dict[str, Any] = {
            "source_id": f"USNO-B1.0 {u['usno_id']}", "plate_source": True, "field_source_id": sid,
            "field_sep_arcsec": round(angular_separation_arcsec(f_ra, f_dec, u["ra_deg"], u["dec_deg"]), 2),
            "field_star_g_mag": (star.get("mags") or {}).get("G"),
            "ra_deg": u["ra_deg"], "dec_deg": u["dec_deg"], "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0,
            "mags": {"G": u["r1_mag"]}, "distance_pc": None, "teff_k": None, "live_enriched": False,
            "usno": u, "vanishing_data": {}, "vetting_archive": {"query_status": {"usno_b1": "ok"}},
        }
        status = modern_status
        if status == "ok" and depth is None:
            status = "no_coverage"
        if status == "ok":
            at_epoch = []
            for c in gaia or []:
                e_ra, e_dec = propagate_position(c["ra_deg"], c["dec_deg"], c["pmra_mas_yr"], c["pmdec_mas_yr"], GAIA_EPOCH, u["position_epoch"])
                at_epoch.append(dict(c, ra_epoch_deg=e_ra, dec_epoch_deg=e_dec))
            g, g_sep = _nearest(u["ra_deg"], u["dec_deg"], at_epoch, "ra_epoch_deg", "dec_epoch_deg")
            p, p_sep = _nearest(u["ra_deg"], u["dec_deg"], ps1_real, "raMean", "decMean")
            # Each Gaia star's path over the plate era against where the plates saw the source
            tracks = sorted(((*track_separation(c, u["plate_ra_deg"], u["plate_dec_deg"]), c) for c in gaia or []
                             if abs(c["dec_deg"] - u["plate_dec_deg"]) * 3600.0 <= 600.0), key=lambda t: (t[0], t[2]["source_id"]))
            t_sep, t_epoch, t_star = tracks[0] if tracks else (None, None, None)
            at_catalogue_epoch = g is not None and g_sep is not None and g_sep <= st["match"]
            on_track = t_star is not None and t_sep is not None and t_sep <= st["match"]
            in_gaia = at_catalogue_epoch or on_track
            in_ps1 = p is not None and p_sep is not None and p_sep <= st["match"]
            epoch, e_ra, e_dec = u["position_epoch"], u["ra_deg"], u["dec_deg"]
            matched = g if at_catalogue_epoch else t_star if on_track else None
            if on_track and not at_catalogue_epoch:
                # compared where and when the star's path meets the plate position
                epoch, e_ra, e_dec = round(float(t_epoch), 2), u["plate_ra_deg"], u["plate_dec_deg"]
            if matched is not None:
                # The channel moves this position back to the plate epoch with these proper motions
                unit.update(ra_deg=matched["ra_deg"], dec_deg=matched["dec_deg"], pmra_mas_yr=matched["pmra_mas_yr"], pmdec_mas_yr=matched["pmdec_mas_yr"])
                if matched["g_mag"] is not None:
                    unit["mags"] = {"G": matched["g_mag"]}
            unit["vanishing_data"] = {
                "usno_epoch": epoch, "usno_mean_epoch": u["mean_epoch"], "usno_ra_deg": e_ra, "usno_dec_deg": e_dec,
                "usno_mag_r": u["r1_mag"], "early_plates_count": u["n_first_epoch_plates"], "present_in_gaia": bool(in_gaia),
                "gaia_match": "catalogue_epoch" if at_catalogue_epoch else "track" if on_track else None,
                "present_in_panstarrs": bool(in_ps1), "panstarrs_depth_mag": depth, "panstarrs_depth_n_objects": n_depth,
                "nearest_gaia": None if g is None else {"source_id": g["source_id"], "sep_arcsec": round(float(g_sep), 2), "g_mag": g["g_mag"]},
                "nearest_gaia_track": {
                    "n_gaia_stars": len(tracks), "plate_epoch": u["plate_epoch"], "plate_era": list(PLATE_ERA),
                    "source_id": None if t_star is None else t_star["source_id"], "sep_arcsec": None if t_sep is None else round(float(t_sep), 2),
                    "epoch": None if t_epoch is None else round(float(t_epoch), 2), "g_mag": None if t_star is None else t_star["g_mag"]},
                "nearest_panstarrs": None if p is None else {"objID": str(p.get("objID")), "sep_arcsec": round(float(p_sep), 2),
                                                             "n_detections": p.get("nDetections"), "r_mag": _float(p.get("rMeanPSFMag"))},
            }
            unit["vetting_archive"]["query_status"].update({"modern_counterpart": "ok", "panstarrs_depth": "ok"})
            # Inside the halo of a bright star and not that star itself: the plates cannot be trusted there
            halos = [(sep, c, r) for sep, _, c in tracks for r in [halo_radius_arcsec(c["g_mag"], st["halo"])]
                     if r is not None and st["match"] < sep <= r]
            if halos:
                h_sep, h_star, h_radius = halos[0]
                status = "near_bright_star"
                unit["plate_fetch"] = {"bright_star": {"source_id": h_star["source_id"], "g_mag": h_star["g_mag"],
                                                       "sep_arcsec": round(float(h_sep), 2), "halo_radius_arcsec": h_radius}}
            elif depth - float(u["r1_mag"]) < st["min_margin"]:
                status = "too_faint"
        unit["plate_fetch"] = {**unit.get("plate_fetch", {}), "status": status}
        units.append(unit)
    return units, counts


def plate_search_status(unit: dict[str, Any]) -> str:
    """ "searched" or the reason a plate source was not: "too_faint", "near_bright_star", "no_coverage", "failed"."""
    status = str((unit.get("plate_fetch") or {}).get("status") or "")
    if status == "ok" and (unit.get("vanishing_data") or {}).get("panstarrs_depth_mag") is not None:
        return "searched"
    if status.startswith("failed"):
        return "failed"
    return status if status in {"too_faint", "near_bright_star", "no_coverage"} else "no_plate_record"


def search_coverage_summary(units: list[dict[str, Any]], cfg: PipelineConfig) -> dict[str, Any]:
    """Counts of searched and not-searched plate sources, in the layout the other survey channels use."""
    ids: dict[str, list[str]] = {"too_faint": [], "near_bright_star": [], "no_coverage": [], "failed": [], "no_plate_record": []}
    n_searched = 0
    for u in units:
        cat = plate_search_status(u)
        if cat == "searched":
            n_searched += 1
        else:
            ids[cat].append(str(u["source_id"]))
    out: dict[str, Any] = {"n_input": len(units), "n_searched": n_searched, "n_not_searched": len(units) - n_searched}
    for cat, lst in ids.items():
        out[f"n_{cat}"] = len(lst)
    out["not_searched_ids"] = {cat: lst for cat, lst in ids.items() if lst}
    return out


# ---------------------------------------------------------------------------------------------------- the survey


def fetch_plate_units(
    stars: list[dict[str, Any]],
    cfg: PipelineConfig,
    cache: RawCache,
    *,
    max_fields: int | None = None,
    progress: ProgressCallback | None = None,
    sleep: Callable[[float], Any] = time.sleep,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Plate-source units of every field. Returns (units, summary).

    Fields are processed in input order, sequentially and at least 1 / `plate_rate_limit_rps` seconds apart per
    uncached request; every answer is cached, so an interrupted run resumes. A failed query is retried
    `plate_max_attempts` times and is never cached.
    """
    st = _settings(cfg)
    interval = 1.0 / st["rps"]
    t_start = time.monotonic()
    n_limit = len(stars) if max_fields is None else max(0, min(int(max_fields), len(stars)))
    units: list[dict[str, Any]] = []
    totals = {"n_usno_sources": 0, "n_on_spike": 0, "n_one_first_epoch_plate": 0, "n_no_first_epoch_plate": 0}
    fields = {"ok": 0, "no_coverage": 0, "failed_modern": 0, "failed_usno": 0}
    failed_fields: list[str] = []
    failed_errors: dict[str, str] = {}
    seen: set[str] = set()
    n_duplicates = 0

    def _retry(fn: Callable[[], tuple[list[dict[str, Any]], str, bool]]) -> list[dict[str, Any]]:
        """Rows of one catalogue query, retried; answers from the cache are not rate-limited."""
        last: Exception | None = None
        for attempt in range(1, int(st["attempts"]) + 1):
            t_req = time.monotonic()
            from_cache = False
            try:
                rows, _, from_cache = fn()
                return rows
            except (ArchiveAccessError, requests.exceptions.RequestException, ValueError, KeyError, TypeError) as exc:
                last = exc
                if isinstance(exc, ArchiveAccessError) and not exc.retryable:
                    break
            finally:
                wait = interval - (time.monotonic() - t_req)
                if wait > 0 and not from_cache:
                    sleep(wait)
            if attempt < st["attempts"]:
                sleep(st["backoff"] * (2 ** (attempt - 1)))
        raise ArchiveAccessError(f"{type(last).__name__}: {str(last)[:200]}", retryable=True)

    for i, star in enumerate(stars[:n_limit]):
        sid = str(star["source_id"])
        ra, dec = float(star["ra_deg"]), float(star["dec_deg"])
        try:
            usno = parse_usno_rows(_retry(lambda: usno_cone(ra, dec, st["radius"], cache)))
        except ArchiveAccessError as exc:
            fields["failed_usno"] += 1
            failed_fields.append(sid)
            failed_errors[sid] = f"usno_b1: {str(exc)[:300]}"
            if progress is not None:
                progress(i + 1, n_limit, sid, f"failed: usno_b1: {str(exc)[:120]}")
            continue
        gaia: list[dict[str, Any]] | None = None
        ps1: list[dict[str, Any]] | None = None
        modern = "ok"
        if dec < PS1_DEC_LIMIT_DEG:
            modern = "no_coverage"
        elif any(u["n_first_epoch_plates"] >= 2 and not u["on_spike"] for u in usno):
            try:
                gaia = _retry(lambda: gaia_cone(ra, dec, st["radius"] + gaia_margin_arcsec(st), cfg, cache, sleep=sleep))
                ps1 = _retry(lambda: ps1_cone(ra, dec, st["radius"] + st["match"], int(st["ps1_min_det"]), cache, sleep=sleep))
            except ArchiveAccessError as exc:
                gaia, ps1, modern = None, None, f"failed: modern catalogues: {str(exc)[:160]}"
        field_units, counts = build_field_units(star, usno, gaia, ps1, st, modern_status=modern)
        for k in totals:
            totals[k] += counts[k]
        statuses = {plate_search_status(u) for u in field_units}
        if modern.startswith("failed"):
            fields["failed_modern"] += 1
            failed_fields.append(sid)
            failed_errors[sid] = modern[len("failed: "):]
        elif "no_coverage" in statuses or modern == "no_coverage":
            fields["no_coverage"] += 1
        else:
            fields["ok"] += 1
        for u in field_units:
            if u["source_id"] in seen:  # two sample stars closer than a field diameter share plate sources
                n_duplicates += 1
                continue
            seen.add(u["source_id"])
            units.append(u)
        if progress is not None:
            progress(i + 1, n_limit, sid, modern.split(":")[0])
    cov = search_coverage_summary(units, cfg)
    summary = {
        "n_fields_input": len(stars), "n_fields_requested": n_limit, "n_fields_ok": fields["ok"], "n_fields_no_coverage": fields["no_coverage"],
        "n_fields_failed_modern": fields["failed_modern"], "n_fields_failed_usno": fields["failed_usno"],
        "n_fields_not_fetched": len(stars) - n_limit, "failed_field_ids": failed_fields, "failed_field_errors": failed_errors, **totals,
        "n_units": len(units), "n_units_searched": cov["n_searched"], "n_duplicate_units_dropped": n_duplicates,
        "field_radius_arcsec": st["radius"], "elapsed_s": round(time.monotonic() - t_start, 1),
    }
    return units, summary
