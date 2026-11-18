"""Breakthrough Listen Open Data Archive: which sample stars lie inside a radio pointing (file index only).

The narrowband radio channel needs observed cadences. The only public narrowband data are Breakthrough Listen
pointings (Green Bank and Parkes), a few thousand positions on the sky. This module asks the archive's file index,
for each star, which files were recorded within a small radius, and keeps those whose half-power beam contains
the star. It downloads no radio data: one fine-resolution scan is 15 to 110 GB.

Statuses of one star (`radio_coverage.status`):
- "in_beam": at least one radio file has the star inside its half-power beam;
- "unverified_pointing": no file is known to hold the star in its beam, but a file whose southern declination the
  index may list wrongly (see `true_declination`) could;
- "near_pointing": files exist within the query radius, but the star is outside every half-power beam;
- "no_coverage": the index holds no radio file within the query radius;
- "failed: ...": the index could not be read (not cached, retried next run).
A star is covered only when "in_beam". Everything else must be reported as not searched by the radio channel.

`plan_star_scans` goes one step further for a covered star, still from the index alone: which fine-resolution scans
contain it (ON), and which scans can serve as OFF: scans the telescope took just before and after on other targets
with the same receiver, or other beams of the Parkes multibeam recorded at the same instant. A narrowband signal
from the star's direction must be in the ON scan and in no OFF scan. The plan gives the files and their sizes, so
the download can be decided before any of it is made.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

from oddstar.errors import ArchiveAccessError
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec
from oddstar.schemas import ManifestRecord

BL_QUERY_URL = "http://seti.berkeley.edu/opendata/api/query-files"  # the API is HTTP only (https redirects to http)
# Larger than the half-power radius of the widest beam in the archive (Parkes at 1.2 GHz: 8 arcmin)
BL_QUERY_RADIUS_DEG = 0.25
# Row limit sent with every request. The index sorts its candidate rows by right ascension, applies the limit, and
# only then cuts to the requested circle (measured 2026-11-18: at a position with files inside 0.25 deg, limit=3
# answered zero rows and limit=500 answered the files). A short answer therefore does not prove that nothing was
# dropped, so the limit is set far above any real answer and an answer anywhere near it is refused.
BL_ROW_LIMIT = 20000
BL_TRUSTED_ROWS = BL_ROW_LIMIT // 4
# Single-dish diameters (m). A file from any other telescope (e.g. the optical APF) is not a radio pointing.
DISH_DIAMETER_M = {"GBT": 100.0, "PARKES": 64.0}
# Half-power beam width of a tapered single dish: about 1.2 wavelength / diameter (GBT 1.4 GHz: 8.8', Parkes: 13.8')
BEAM_FACTOR = 1.2
SPEED_OF_LIGHT_M_S = 299_792_458.0
MIN_REQUEST_INTERVAL_S = 0.5
# Finest spectral product of a scan (about 3 Hz channels): the only one a narrowband drift search can use.
# GBT and the Parkes multibeam name it product 0000 (.0000.h5, .0000.fil, .0000.hires.hdf). The Parkes ultra-wideband
# receiver writes one .fil per 128 MHz sub-band named by start time; its size (18 or 19 spectra of 2**25 channels of
# 4 bytes plus a header) shows 3.8 Hz channels.
# The Parkes single-beam receivers of 2017 name it "<...>_fine.h5".
FINE_RESOLUTION_URL = re.compile(r"\.0000\.(h5|fil)$|\.0000\.hires\.hdf$|/PKSUWL/.*T\d\d:\d\d:\d\d\.\d{3}\.fil$|/Parkes_[^/]*_fine\.h5$")
# Scans of one observing block follow each other within minutes; an OFF scan must start within this of the ON scan
OFF_SCAN_WINDOW_DAYS = 0.02
# An OFF pointing must be at least this many half-power radii from the star, so the star is far outside its beam
OFF_MIN_SEPARATION_BEAM_RADII = 6.0
# Files of one scan (one per recorder node, band or beam) share a start time to within this
SAME_SCAN_MJD_TOL = 5.0 / 86400.0
# The Parkes multibeam records 13 beams at once; this many other beams of the same instant make a comparison
MIN_SIMULTANEOUS_OFF_BEAMS = 2

# --- Southern declinations of the index (measured 2026-11-18).
# Right ascensions and northern declinations are right. Many southern declinations are not: for a true declination
# -(D + f), with D whole degrees and 0 <= f < 1, the index lists -(D + c - f). The fraction of a degree is mirrored,
# and c depends on how the file family was entered:
# - Green Bank "guppi_<day>_<second>_<target>_<scan>.gpuspec" products and Parkes ultra-wideband files: c = 40' 40"
#   (HIP 3026, HIP 3909 and HIP 17147 at Green Bank, HIP 1837 at Parkes, each against its Gaia DR3 position);
# - Parkes single-beam "Parkes_<day>_<second>_<target>_fine.h5" files: c = 0, the degrees were taken as negative
#   and the minutes and seconds added as positive (HIP 8102, HIP 113368);
# - Parkes multibeam files (collate_mb): right as listed (G238.23-51.24 and G306.02-45.34: the listed beam is
#   29 arcmin from the position its Galactic name gives, the spacing of the inner ring of beams).
# A southern file of any other family is "unverified": its declination is not trusted in either direction.
# True declinations between -1 and 0 degrees cannot be recovered (the degree field carries no sign there), so a
# pointing in that strip is not found.
DEC_MIRROR_40_40 = 40.0 / 60.0 + 40.0 / 3600.0
DEC_MIRRORS = {"mirror_40_40": DEC_MIRROR_40_40, "mirror_0": 0.0}
# (telescope, URL pattern, rule); the first match decides
_DEC_RULES: tuple[tuple[str, re.Pattern[str], str], ...] = (
    ("PARKES", re.compile(r"/collate_mb"), "as_listed"),
    ("PARKES", re.compile(r"/PKSUWL/"), "mirror_40_40"),
    ("PARKES", re.compile(r"/Parkes_\d+_\d+_[^/]*_fine\.h5$"), "mirror_0"),
    ("GBT", re.compile(r"guppi_\d+_\d+_[^/]*\.gpuspec\."), "mirror_40_40"),
)

ProgressCallback = Callable[[int, int, str, str], None]


def unmirror_declination(index_decl: float, c: float) -> float:
    """True declination of a southern index value listed as -(D + c - f)."""
    a = abs(index_decl)
    whole = math.ceil(round(a - c, 6))
    return -(2.0 * whole + c - a)


def mirror_declination(dec: float, c: float) -> float | None:
    """Where the index lists a true southern declination under mirror constant `c`; None in the strip above -1."""
    a = abs(dec)
    whole = math.floor(a)
    return None if dec >= 0.0 or whole == 0 else -(whole + c - (a - whole))


def declination_rule(telescope: str, url: str, index_decl: float) -> str:
    """ "as_listed", "mirror_40_40", "mirror_0" or "unverified" for one file record."""
    if index_decl >= 0.0:
        return "as_listed"
    for tel, pattern, rule in _DEC_RULES:
        if tel == str(telescope).upper() and pattern.search(url):
            return rule
    return "unverified"


def true_declination(telescope: str, url: str, index_decl: float) -> tuple[float, str]:
    """(declination, rule) of a file record. An "unverified" record keeps the listed value, which must not be trusted."""
    rule = declination_rule(telescope, url, index_decl)
    if rule in DEC_MIRRORS:
        return unmirror_declination(index_decl, DEC_MIRRORS[rule]), rule
    return index_decl, rule


def candidate_declinations(index_decl: float) -> list[float]:
    """Every declination an unverified southern record could truly have: as listed, or mirrored either way."""
    return [index_decl, *[unmirror_declination(index_decl, c) for c in DEC_MIRRORS.values()]]


def index_query_declinations(dec: float, radius_deg: float) -> list[float]:
    """Declinations at which the index must be asked to find every file within `radius_deg` of a true declination.

    North of the equator that is the declination itself. South of it, a file can also be listed at the mirrored
    value of either family, and the part of the circle that lies in a neighbouring whole degree is mirrored into
    a different place, which is asked separately.
    """
    out = [dec]
    if dec >= 0.0:
        return out
    for c in DEC_MIRRORS.values():
        for d in (dec, dec - radius_deg, dec + radius_deg):
            m = mirror_declination(d, c)
            if m is not None and all(abs(m - o) > 1e-6 for o in out) and (d == dec or math.floor(abs(d)) != math.floor(abs(dec))):
                out.append(m)
    return out


def half_power_radius_arcsec(telescope: str, center_freq_mhz: float) -> float | None:
    """Half of the half-power beam width at the file's centre frequency, or None for a non-radio or unknown telescope."""
    dish = DISH_DIAMETER_M.get(str(telescope).upper())
    if dish is None or not (center_freq_mhz and center_freq_mhz > 0):
        return None
    hpbw_rad = BEAM_FACTOR * SPEED_OF_LIGHT_M_S / (center_freq_mhz * 1e6) / dish
    return 0.5 * math.degrees(hpbw_rad) * 3600.0


def query_bl_files(
    ra: float, dec: float, cache: RawCache, *, radius_deg: float = BL_QUERY_RADIUS_DEG, sleep: ot.Sleeper = time.sleep,
) -> tuple[list[dict[str, Any]], ManifestRecord | None, bool]:
    """File records within `radius_deg` of a position. Returns (rows, manifest, from_cache).

    Raises ArchiveAccessError when the index cannot be read or does not answer "success"; such an answer is never
    stored as an empty result.
    """
    params = {"target": "", "pos-ra": f"{ra:.5f}", "pos-dec": f"{dec:.5f}", "pos-rad": f"{radius_deg:.3f}", "limit": str(BL_ROW_LIMIT)}
    desc = f"GET {BL_QUERY_URL} {json.dumps(params, sort_keys=True)}"
    sig = hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("breakthrough_listen", "opendata_api", sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(cache.read_bytes(cached).decode("utf-8")), cached, True
        except (FileNotFoundError, OSError, ValueError):
            pass  # evicted or unreadable blob: ask again
    resp = ot._get_with_retries(BL_QUERY_URL, params, timeout_s=30.0, attempts=3, sleep=sleep)
    try:
        body = json.loads(resp.content.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise ArchiveAccessError(f"Breakthrough Listen index answer is not JSON: {str(exc)[:120]}", retryable=True) from exc
    if not isinstance(body, dict) or body.get("result") != "success" or not isinstance(body.get("data"), list):
        raise ArchiveAccessError(f"Breakthrough Listen index did not answer 'success': {str(body)[:160]}", retryable=True)
    rows = body["data"]
    if len(rows) >= BL_TRUSTED_ROWS:
        raise ArchiveAccessError(f"{len(rows)} rows are too close to the index row limit {BL_ROW_LIMIT}; the answer may be cut short",
                                 retryable=False)
    man = cache.put(archive="breakthrough_listen", release="opendata_api", request_signature=sig, query_text=desc,
                    payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
    return rows, man, False


def query_star_rows(
    ra: float, dec: float, cache: RawCache, *, radius_deg: float = BL_QUERY_RADIUS_DEG, sleep: ot.Sleeper = time.sleep,
    min_interval_s: float = 0.0,
) -> tuple[list[dict[str, Any]], list[ManifestRecord], bool]:
    """Every file record that could lie within `radius_deg` of a true position, wherever the index lists it.

    Asks the index at each declination of `index_query_declinations` and merges the answers (one record per file).
    Every request that went to the network is followed by a pause that fills `min_interval_s`.
    Returns (rows, manifests, all_from_cache).
    """
    merged: dict[Any, dict[str, Any]] = {}
    manifests: list[ManifestRecord] = []
    all_cached = True
    for q_dec in index_query_declinations(dec, radius_deg):
        t_req = time.monotonic()
        hit = False
        try:
            rows, man, hit = query_bl_files(ra, q_dec, cache, radius_deg=radius_deg, sleep=sleep)
        finally:
            wait = min_interval_s - (time.monotonic() - t_req)
            if not hit and wait > 0:
                sleep(wait)
        all_cached = all_cached and hit
        if man is not None:
            manifests.append(man)
        for r in rows:
            merged.setdefault(r.get("id", r.get("url")), r)
    return list(merged.values()), manifests, all_cached


def query_bl_time_window(
    telescope: str, mjd_start: float, mjd_end: float, cache: RawCache, *, sleep: ot.Sleeper = time.sleep,
) -> tuple[list[dict[str, Any]], ManifestRecord | None, bool]:
    """Every file record of one telescope that starts inside a time window. Returns (rows, manifest, from_cache)."""
    params = {"target": "", "telescopes": str(telescope).lower(), "time-start": f"{mjd_start:.5f}", "time-end": f"{mjd_end:.5f}",
              "limit": str(BL_ROW_LIMIT)}
    desc = f"GET {BL_QUERY_URL} {json.dumps(params, sort_keys=True)}"
    sig = hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("breakthrough_listen", "opendata_api", sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(cache.read_bytes(cached).decode("utf-8")), cached, True
        except (FileNotFoundError, OSError, ValueError):
            pass
    resp = ot._get_with_retries(BL_QUERY_URL, params, timeout_s=30.0, attempts=3, sleep=sleep)
    try:
        body = json.loads(resp.content.decode("utf-8"))
    except (ValueError, UnicodeDecodeError) as exc:
        raise ArchiveAccessError(f"Breakthrough Listen index answer is not JSON: {str(exc)[:120]}", retryable=True) from exc
    if not isinstance(body, dict) or body.get("result") != "success" or not isinstance(body.get("data"), list):
        raise ArchiveAccessError(f"Breakthrough Listen index did not answer 'success': {str(body)[:160]}", retryable=True)
    rows = body["data"]
    if len(rows) >= BL_TRUSTED_ROWS:
        raise ArchiveAccessError(f"{len(rows)} rows are too close to the index row limit {BL_ROW_LIMIT}; the time window answer may be cut short",
                                 retryable=False)
    man = cache.put(archive="breakthrough_listen", release="opendata_api", request_signature=sig, query_text=desc,
                    payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
    return rows, man, False


def _fine_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Fine-resolution radio files with a usable position, time and frequency.

    The declination is the true one (`true_declination`); a file whose declination is unverified is left out, so it
    is neither an ON scan nor a comparison scan.
    """
    out = []
    for r in rows:
        try:
            rec = {"url": str(r["url"]), "target": str(r.get("target", "")), "telescope": str(r.get("telescope", "")).upper(),
                   "center_freq_mhz": float(r["center_freq"]), "ra_deg": float(r["ra"]), "mjd": float(r["mjd"]),
                   "size_bytes": int(r.get("size") or 0), "md5sum": str(r.get("md5sum") or "")}
            rec["dec_deg"], rec["dec_rule"] = true_declination(rec["telescope"], rec["url"], float(r["decl"]))
        except (KeyError, TypeError, ValueError):
            continue
        if (rec["dec_rule"] != "unverified" and FINE_RESOLUTION_URL.search(rec["url"])
                and half_power_radius_arcsec(rec["telescope"], rec["center_freq_mhz"]) is not None):
            out.append(rec)
    return out


def _file_entry(r: dict[str, Any]) -> dict[str, Any]:
    return {"url": r["url"], "size_bytes": r["size_bytes"], "md5sum": r["md5sum"], "center_freq_mhz": r["center_freq_mhz"]}


def _one_file_per_frequency(files: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """The archive often holds a scan's fine product twice (HDF5 and filterbank); keep one per centre frequency, HDF5 first."""
    best: dict[float, dict[str, Any]] = {}
    for f in sorted(files, key=lambda x: (0 if x["url"].endswith((".h5", ".hdf")) else 1, x["url"])):
        best.setdefault(round(f["center_freq_mhz"], 1), f)
    return [best[k] for k in sorted(best)]


def _group_scans(files: list[dict[str, Any]]) -> list[list[dict[str, Any]]]:
    """Files of one target that start together (within `SAME_SCAN_MJD_TOL` of the first) are one scan."""
    groups: list[list[dict[str, Any]]] = []
    for f in sorted(files, key=lambda x: (x["target"], x["mjd"])):
        if groups and groups[-1][0]["target"] == f["target"] and f["mjd"] - groups[-1][0]["mjd"] <= SAME_SCAN_MJD_TOL:
            groups[-1].append(f)
        else:
            groups.append([f])
    return groups


def _off_entry(target: str, files: list[dict[str, Any]], mjd_on: float, sep_arcsec: float) -> dict[str, Any]:
    o_mjd = min(f["mjd"] for f in files)
    return {"target": target, "mjd": o_mjd, "minutes_from_on": round((o_mjd - mjd_on) * 1440.0, 1), "sep_arcsec": round(sep_arcsec, 1),
            "files": [_file_entry(f) for f in files], "bytes": sum(f["size_bytes"] for f in files)}


def plan_star_scans(
    ra: float, dec: float, rows: list[dict[str, Any]], cache: RawCache, *, sleep: ot.Sleeper = time.sleep,
) -> list[dict[str, Any]]:
    """ON scans that contain the position, each with the OFF scans a comparison needs.

    `rows` is the index answer around the position (`query_star_rows`). A scan is the set of fine-resolution files of
    one target that start together; only files whose half-power beam contains the position are ON files. An OFF
    scan is another target observed by the same telescope, at least `OFF_MIN_SEPARATION_BEAM_RADII` half-power radii
    from the position, holding a file at the frequency of every ON file (to 1 MHz). Two kinds are listed: the nearest
    scan before and after within `OFF_SCAN_WINDOW_DAYS`, and beams recorded at the same instant (Parkes multibeam).
    `comparison` names the set a search should use: "simultaneous_beams" when at least
    `MIN_SIMULTANEOUS_OFF_BEAMS` exist, else "before_after" when both exist, else "incomplete".
    `bytes_on_plus_offs` is what that comparison downloads.
    """
    in_beam: list[dict[str, Any]] = []
    for r in _fine_rows(rows):
        radius = half_power_radius_arcsec(r["telescope"], r["center_freq_mhz"])
        sep = angular_separation_arcsec(ra, dec, r["ra_deg"], r["dec_deg"])
        if radius is not None and sep <= radius:
            in_beam.append(dict(r, sep_arcsec=sep))
    plans: list[dict[str, Any]] = []
    for group in sorted(_group_scans(in_beam), key=lambda g: g[0]["mjd"]):
        target, tel, mjd = group[0]["target"], group[0]["telescope"], group[0]["mjd"]
        on_files = _one_file_per_frequency([f for f in group if f["telescope"] == tel])
        freqs = [round(f["center_freq_mhz"], 1) for f in on_files]
        max_radius = max(half_power_radius_arcsec(tel, f["center_freq_mhz"]) or 0.0 for f in on_files)
        window, _, _ = query_bl_time_window(tel, mjd - OFF_SCAN_WINDOW_DAYS, mjd + OFF_SCAN_WINDOW_DAYS, cache, sleep=sleep)
        far = [r for r in _fine_rows(window) if r["telescope"] == tel and r["target"] != target
               and angular_separation_arcsec(ra, dec, r["ra_deg"], r["dec_deg"]) >= OFF_MIN_SEPARATION_BEAM_RADII * max_radius]
        before: dict[str, Any] | None = None
        after: dict[str, Any] | None = None
        simultaneous: list[dict[str, Any]] = []
        for o_group in _group_scans(far):
            per_freq = {round(o["center_freq_mhz"], 1): o for o in _one_file_per_frequency(o_group)}
            matched = [next((o for k, o in per_freq.items() if abs(k - f) <= 1.0), None) for f in freqs]
            if any(m is None for m in matched):
                continue  # this scan does not cover every ON frequency
            sep = angular_separation_arcsec(ra, dec, o_group[0]["ra_deg"], o_group[0]["dec_deg"])
            entry = _off_entry(o_group[0]["target"], [m for m in matched if m is not None], mjd, sep)
            if abs(entry["mjd"] - mjd) <= SAME_SCAN_MJD_TOL:
                simultaneous.append(entry)
            elif entry["mjd"] < mjd:
                before = entry if before is None or entry["mjd"] > before["mjd"] else before
            else:
                after = entry if after is None or entry["mjd"] < after["mjd"] else after
        simultaneous = sorted(simultaneous, key=lambda e: (e["sep_arcsec"], e["target"]))[:MIN_SIMULTANEOUS_OFF_BEAMS]
        if len(simultaneous) >= MIN_SIMULTANEOUS_OFF_BEAMS:
            comparison, used = "simultaneous_beams", simultaneous
        elif before is not None and after is not None:
            comparison, used = "before_after", [before, after]
        else:
            comparison, used = "incomplete", [e for e in (before, after) if e is not None] + simultaneous
        bytes_on = sum(f["size_bytes"] for f in on_files)
        plans.append({
            "target": target, "telescope": tel, "mjd": mjd, "sep_arcsec": round(min(f["sep_arcsec"] for f in on_files), 1),
            "center_freqs_mhz": freqs, "on_files": [_file_entry(f) for f in on_files], "bytes_on": bytes_on,
            "off_before": before, "off_after": after, "off_simultaneous": simultaneous,
            "comparison": comparison, "complete": comparison != "incomplete",
            "bytes_on_plus_offs": bytes_on + sum(e["bytes"] for e in used),
        })
    return plans


def beam_coverage(ra: float, dec: float, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Group file records into pointings and say which contain the position inside their half-power beam.

    A pointing is (target, telescope, receiver band to the nearest 100 MHz of the centre frequency). Its separation
    is the smallest separation of any of its files; it is in beam when that is within the half-power radius.
    Separations use the true declination of each file (`true_declination`). A pointing with a file of unverified
    declination is never in beam; its separation is the smallest over the declinations the file could have, and
    `possibly_in_beam` says whether one of them would hold the position.
    """
    pointings: dict[tuple[str, str, int], dict[str, Any]] = {}
    n_radio = 0
    for r in rows:
        tel = str(r.get("telescope", "")).upper()
        try:
            f_mhz = float(r.get("center_freq"))
            p_ra, listed = float(r.get("ra")), float(r.get("decl"))
        except (TypeError, ValueError):
            continue
        radius = half_power_radius_arcsec(tel, f_mhz)
        if radius is None:
            continue
        n_radio += 1
        p_dec, rule = true_declination(tel, str(r.get("url", "")), listed)
        if rule == "unverified":
            sep = min(angular_separation_arcsec(ra, dec, p_ra, d) for d in candidate_declinations(listed))
        else:
            sep = angular_separation_arcsec(ra, dec, p_ra, p_dec)
        key = (str(r.get("target", "")), tel, int(round(f_mhz / 100.0)) * 100)
        p = pointings.setdefault(key, {"target": key[0], "telescope": tel, "band_mhz": key[2], "n_files": 0, "total_bytes": 0,
                                       "sep_arcsec": sep, "half_power_radius_arcsec": radius, "mjd_min": None, "mjd_max": None,
                                       "file_types": [], "dec_rules": []})
        if rule not in p["dec_rules"]:
            p["dec_rules"].append(rule)
        p["n_files"] += 1
        p["total_bytes"] += int(r.get("size") or 0)
        if sep < p["sep_arcsec"]:
            p["sep_arcsec"], p["half_power_radius_arcsec"] = sep, radius
        mjd = r.get("mjd")
        if isinstance(mjd, (int, float)):
            p["mjd_min"] = mjd if p["mjd_min"] is None else min(p["mjd_min"], mjd)
            p["mjd_max"] = mjd if p["mjd_max"] is None else max(p["mjd_max"], mjd)
        ft = str(r.get("file_type", ""))
        if ft and ft not in p["file_types"]:
            p["file_types"].append(ft)
    out = []
    for p in sorted(pointings.values(), key=lambda x: (x["sep_arcsec"], x["target"], x["band_mhz"])):
        inside = bool(p["sep_arcsec"] <= p["half_power_radius_arcsec"])
        p["dec_rules"] = sorted(p["dec_rules"])
        p["dec_unverified"] = "unverified" in p["dec_rules"]
        p["in_beam"] = inside and not p["dec_unverified"]
        p["possibly_in_beam"] = inside and p["dec_unverified"]
        p["sep_arcsec"] = round(p["sep_arcsec"], 1)
        p["half_power_radius_arcsec"] = round(p["half_power_radius_arcsec"], 1)
        p["file_types"] = sorted(p["file_types"])
        out.append(p)
    in_beam = [p for p in out if p["in_beam"]]
    possible = [p for p in out if p["possibly_in_beam"]]
    status = "in_beam" if in_beam else "unverified_pointing" if possible else "near_pointing" if out else "no_coverage"
    return {"status": status, "n_index_rows": len(rows), "n_radio_files": n_radio, "index_answer_full": len(rows) >= BL_TRUSTED_ROWS,
            "n_pointings": len(out), "n_pointings_in_beam": len(in_beam), "n_pointings_unverified": sum(p["dec_unverified"] for p in out),
            "pointings": out}


def sample_radio_coverage(
    stars: list[dict[str, Any]], cache: RawCache, *, sleep: Callable[[float], Any] = time.sleep,
    progress: ProgressCallback | None = None, min_interval_s: float = MIN_REQUEST_INTERVAL_S,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Radio-pointing coverage of every star. Returns (per-star records in input order, summary).

    A star south of the equator is asked for at several declinations (`query_star_rows`). Cached answers are not
    rate-limited, so an interrupted census resumes quickly. A failed star is recorded as
    "failed: ..." and counted separately; it is neither covered nor uncovered.
    """
    records: list[dict[str, Any]] = []
    counts = {"in_beam": 0, "unverified_pointing": 0, "near_pointing": 0, "no_coverage": 0, "failed": 0}
    t_start = time.monotonic()
    for i, s in enumerate(stars):
        sid = str(s["source_id"])
        ra, dec = float(s["ra_deg"]), float(s["dec_deg"])
        try:
            rows, mans, _ = query_star_rows(ra, dec, cache, sleep=sleep, min_interval_s=min_interval_s)
            rec = beam_coverage(ra, dec, rows)
            rec["manifest_ids"] = [m.manifest_id for m in mans]
        except ArchiveAccessError as exc:
            rec = {"status": f"failed: {str(exc)[:200]}", "pointings": []}
        rec.update({"source_id": sid, "ra_deg": ra, "dec_deg": dec})
        key = "failed" if rec["status"].startswith("failed") else rec["status"]
        counts[key] += 1
        records.append(rec)
        if progress is not None:
            progress(i + 1, len(stars), sid, rec["status"])
    summary = {
        "n_input": len(stars), "n_in_beam": counts["in_beam"], "n_unverified_pointing": counts["unverified_pointing"],
        "n_near_pointing": counts["near_pointing"], "n_no_coverage": counts["no_coverage"], "n_failed": counts["failed"],
        "in_beam_ids": [r["source_id"] for r in records if r["status"] == "in_beam"],
        "unverified_pointing_ids": [r["source_id"] for r in records if r["status"] == "unverified_pointing"],
        "near_pointing_ids": [r["source_id"] for r in records if r["status"] == "near_pointing"],
        "failed_ids": [r["source_id"] for r in records if r["status"].startswith("failed")],
        "query_radius_deg": BL_QUERY_RADIUS_DEG, "elapsed_s": round(time.monotonic() - t_start, 1),
    }
    return records, summary


def main(argv: list[str] | None = None) -> int:
    """`python -m oddstar.io.bl_opendata --region halo --limit 1000`: census of a live sample, written under runs/."""
    from oddstar.config import load_config, repo_root
    from oddstar.io.live_survey import fetch_live_halo_bulge_sample

    ap = argparse.ArgumentParser(description="Which stars of a live sample lie inside a Breakthrough Listen radio pointing")
    ap.add_argument("--region", default="halo", choices=["halo", "bulge", "local"])
    ap.add_argument("--limit", type=int, default=1000)
    ap.add_argument("--after-random-index", type=int, default=None)
    ap.add_argument("--out", default=None, help="output JSON (default runs/radio-coverage-<region>-<limit>/bl_coverage.json)")
    ap.add_argument("--plan", action="store_true",
                    help="read the census written earlier and list, for each star in a beam, its ON scans with their OFF scans (index only)")
    args = ap.parse_args(argv)
    root = repo_root()
    cfg, short_hash, _ = load_config(root=root)
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    if args.plan:
        census_path = Path(args.out) if args.out else root / cfg.storage.runs_dir / f"radio-coverage-{args.region}-{args.limit}" / "bl_coverage.json"
        census = json.loads(census_path.read_text(encoding="utf-8"))
        stars_out = []
        for rec in census["stars"]:
            if rec.get("status") != "in_beam":
                continue
            ra, dec = float(rec["ra_deg"]), float(rec["dec_deg"])
            try:
                rows, _, _ = query_star_rows(ra, dec, cache)
                scans = plan_star_scans(ra, dec, rows, cache)
                entry: dict[str, Any] = {"status": "ok", "scans": scans}
            except ArchiveAccessError as exc:
                entry = {"status": f"failed: {str(exc)[:200]}", "scans": []}
            complete = [s for s in entry["scans"] if s["complete"]]
            entry.update({"source_id": rec["source_id"], "ra_deg": ra, "dec_deg": dec, "n_on_scans": len(entry["scans"]),
                          "n_complete_scans": len(complete),
                          "smallest_complete_bytes": min((s["bytes_on_plus_offs"] for s in complete), default=None)})
            stars_out.append(entry)
            print(f"  plan {rec['source_id']}: {entry['n_on_scans']} ON scan(s), {entry['n_complete_scans']} with a full comparison", file=sys.stderr)
        plan_path = census_path.with_name("bl_scan_plan.json")
        plan_path.write_text(json.dumps({"census": str(census_path), "stars": stars_out}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps({"out": str(plan_path), "stars": [{k: e[k] for k in ("source_id", "status", "n_on_scans", "n_complete_scans",
                                                                             "smallest_complete_bytes")} for e in stars_out]}, indent=2))
        return 0
    stars, _ = fetch_live_halo_bulge_sample(cfg, cache, limit=args.limit, region=args.region, after_random_index=args.after_random_index)

    def _progress(done: int, total: int, sid: str, status: str) -> None:
        if done == total or done % 100 == 0:
            print(f"  radio index {done}/{total} (last {sid}: {status.split(':')[0]})", file=sys.stderr)

    records, summary = sample_radio_coverage(stars, cache, progress=_progress)
    out = Path(args.out) if args.out else root / cfg.storage.runs_dir / f"radio-coverage-{args.region}-{args.limit}" / "bl_coverage.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {"region": args.region, "limit": args.limit, "after_random_index": args.after_random_index, "config_hash": short_hash,
               "summary": summary, "stars": records}
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(out), **summary}, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
