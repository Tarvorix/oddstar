"""Breakthrough Listen Open Data Archive: which sample stars lie inside a radio pointing (file index only).

The narrowband radio channel needs observed cadences. The only public narrowband data are Breakthrough Listen
pointings (Green Bank and Parkes), a few thousand positions on the sky. This module asks the archive's file index,
for each star, which files were recorded within a small radius, and keeps those whose half-power beam contains
the star. It downloads no radio data: one fine-resolution scan is 15 to 110 GB.

Statuses of one star (`radio_coverage.status`):
- "in_beam": at least one radio file has the star inside its half-power beam;
- "near_pointing": files exist within the query radius, but the star is outside every half-power beam;
- "no_coverage": the index holds no radio file within the query radius;
- "failed: ...": the index could not be read (not cached, retried next run).
A star is covered only when "in_beam". Everything else must be reported as not searched by the radio channel.

`plan_star_scans` goes one step further for a covered star, still from the index alone: which fine-resolution scans
contain it (ON), and which scans the telescope took just before and after on other targets with the same receiver
(OFF). A narrowband signal from the star's direction must be in the ON scan and in neither OFF scan. The plan gives
the files and their sizes, so the download can be decided before any of it is made.
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
# The index returns at most this many rows per request; an answer of this size may be cut short
BL_ROW_LIMIT = 500
# Single-dish diameters (m). A file from any other telescope (e.g. the optical APF) is not a radio pointing.
DISH_DIAMETER_M = {"GBT": 100.0, "PARKES": 64.0}
# Half-power beam width of a tapered single dish: about 1.2 wavelength / diameter (GBT 1.4 GHz: 8.8', Parkes: 13.8')
BEAM_FACTOR = 1.2
SPEED_OF_LIGHT_M_S = 299_792_458.0
MIN_REQUEST_INTERVAL_S = 0.5
# Finest spectral product of a scan (about 3 Hz channels): the only one a narrowband drift search can use
FINE_RESOLUTION_URL = re.compile(r"\.0000\.(h5|fil)$")
# Scans of one observing block follow each other within minutes; an OFF scan must start within this of the ON scan
OFF_SCAN_WINDOW_DAYS = 0.02
# An OFF pointing must be at least this many half-power radii from the star, so the star is far outside its beam
OFF_MIN_SEPARATION_BEAM_RADII = 6.0
# Files of one scan (one per recorder node or product) share a start time to within this
SAME_SCAN_MJD_TOL = 2.0 / 86400.0

ProgressCallback = Callable[[int, int, str, str], None]


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
    man = cache.put(archive="breakthrough_listen", release="opendata_api", request_signature=sig, query_text=desc,
                    payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
    return rows, man, False


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
    if len(rows) >= BL_ROW_LIMIT:
        raise ArchiveAccessError(f"{len(rows)} rows reach the index row limit; the time window answer may be cut short", retryable=False)
    man = cache.put(archive="breakthrough_listen", release="opendata_api", request_signature=sig, query_text=desc,
                    payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
    return rows, man, False


def _fine_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Fine-resolution radio files with a usable position, time and frequency."""
    out = []
    for r in rows:
        try:
            rec = {"url": str(r["url"]), "target": str(r.get("target", "")), "telescope": str(r.get("telescope", "")).upper(),
                   "center_freq_mhz": float(r["center_freq"]), "ra_deg": float(r["ra"]), "dec_deg": float(r["decl"]), "mjd": float(r["mjd"]),
                   "size_bytes": int(r.get("size") or 0), "md5sum": str(r.get("md5sum") or "")}
        except (KeyError, TypeError, ValueError):
            continue
        if FINE_RESOLUTION_URL.search(rec["url"]) and half_power_radius_arcsec(rec["telescope"], rec["center_freq_mhz"]) is not None:
            out.append(rec)
    return out


def _file_entry(r: dict[str, Any]) -> dict[str, Any]:
    return {"url": r["url"], "size_bytes": r["size_bytes"], "md5sum": r["md5sum"], "center_freq_mhz": r["center_freq_mhz"]}


def plan_star_scans(
    ra: float, dec: float, rows: list[dict[str, Any]], cache: RawCache, *, sleep: ot.Sleeper = time.sleep,
) -> list[dict[str, Any]]:
    """ON scans that contain the position, each with the nearest earlier and later OFF scan of the same receiver.

    `rows` is the index answer around the position (`query_bl_files`). A scan is the set of fine-resolution files of
    one target that start together; only files whose half-power beam contains the position are ON files. An OFF
    scan is another target observed by the same telescope within `OFF_SCAN_WINDOW_DAYS`, at least
    `OFF_MIN_SEPARATION_BEAM_RADII` half-power radii from the position, holding a file at the frequency of every ON
    file (to 1 MHz). `complete` says both OFF scans exist; `bytes_on_plus_offs` is what a full comparison downloads.
    """
    scans: dict[tuple[str, str, int], list[dict[str, Any]]] = {}
    for r in _fine_rows(rows):
        radius = half_power_radius_arcsec(r["telescope"], r["center_freq_mhz"])
        sep = angular_separation_arcsec(ra, dec, r["ra_deg"], r["dec_deg"])
        if radius is not None and sep <= radius:
            scans.setdefault((r["target"], r["telescope"], int(round(r["mjd"] / SAME_SCAN_MJD_TOL))), []).append(dict(r, sep_arcsec=sep))
    plans: list[dict[str, Any]] = []
    for (target, tel, _), on_files in sorted(scans.items(), key=lambda kv: min(f["mjd"] for f in kv[1])):
        mjd = min(f["mjd"] for f in on_files)
        freqs = sorted({round(f["center_freq_mhz"], 1) for f in on_files})
        max_radius = max(half_power_radius_arcsec(tel, f["center_freq_mhz"]) or 0.0 for f in on_files)
        window, _, _ = query_bl_time_window(tel, mjd - OFF_SCAN_WINDOW_DAYS, mjd + OFF_SCAN_WINDOW_DAYS, cache, sleep=sleep)
        others: dict[tuple[str, int], list[dict[str, Any]]] = {}
        for r in _fine_rows(window):
            if r["telescope"] != tel or r["target"] == target or abs(r["mjd"] - mjd) <= SAME_SCAN_MJD_TOL:
                continue
            if angular_separation_arcsec(ra, dec, r["ra_deg"], r["dec_deg"]) < OFF_MIN_SEPARATION_BEAM_RADII * max_radius:
                continue
            others.setdefault((r["target"], int(round(r["mjd"] / SAME_SCAN_MJD_TOL))), []).append(r)
        offs = {"before": None, "after": None}
        for (o_target, _), o_files in others.items():
            matched = [next((o for o in o_files if abs(o["center_freq_mhz"] - f) <= 1.0), None) for f in freqs]
            if any(m is None for m in matched):
                continue  # this scan does not cover every ON frequency
            o_mjd = min(o["mjd"] for o in o_files)
            side = "before" if o_mjd < mjd else "after"
            best = offs[side]
            if best is None or abs(o_mjd - mjd) < abs(best["mjd"] - mjd):
                offs[side] = {"target": o_target, "mjd": o_mjd, "minutes_from_on": round((o_mjd - mjd) * 1440.0, 1),
                              "files": [_file_entry(m) for m in matched if m is not None],
                              "bytes": sum(m["size_bytes"] for m in matched if m is not None)}
        bytes_on = sum(f["size_bytes"] for f in on_files)
        complete = offs["before"] is not None and offs["after"] is not None
        plans.append({
            "target": target, "telescope": tel, "mjd": mjd, "sep_arcsec": round(min(f["sep_arcsec"] for f in on_files), 1),
            "center_freqs_mhz": freqs, "on_files": [_file_entry(f) for f in sorted(on_files, key=lambda x: x["center_freq_mhz"])],
            "bytes_on": bytes_on, "off_before": offs["before"], "off_after": offs["after"], "complete": complete,
            "bytes_on_plus_offs": bytes_on + sum(o["bytes"] for o in offs.values() if o is not None),
        })
    return plans


def beam_coverage(ra: float, dec: float, rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Group file records into pointings and say which contain the position inside their half-power beam.

    A pointing is (target, telescope, receiver band to the nearest 100 MHz of the centre frequency). Its separation
    is the smallest separation of any of its files; it is in beam when that is within the half-power radius.
    """
    pointings: dict[tuple[str, str, int], dict[str, Any]] = {}
    n_radio = 0
    for r in rows:
        tel = str(r.get("telescope", "")).upper()
        try:
            f_mhz = float(r.get("center_freq"))
            p_ra, p_dec = float(r.get("ra")), float(r.get("decl"))
        except (TypeError, ValueError):
            continue
        radius = half_power_radius_arcsec(tel, f_mhz)
        if radius is None:
            continue
        n_radio += 1
        sep = angular_separation_arcsec(ra, dec, p_ra, p_dec)
        key = (str(r.get("target", "")), tel, int(round(f_mhz / 100.0)) * 100)
        p = pointings.setdefault(key, {"target": key[0], "telescope": tel, "band_mhz": key[2], "n_files": 0, "total_bytes": 0,
                                       "sep_arcsec": sep, "half_power_radius_arcsec": radius, "mjd_min": None, "mjd_max": None,
                                       "file_types": []})
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
        p["in_beam"] = bool(p["sep_arcsec"] <= p["half_power_radius_arcsec"])
        p["sep_arcsec"] = round(p["sep_arcsec"], 1)
        p["half_power_radius_arcsec"] = round(p["half_power_radius_arcsec"], 1)
        p["file_types"] = sorted(p["file_types"])
        out.append(p)
    in_beam = [p for p in out if p["in_beam"]]
    status = "in_beam" if in_beam else ("near_pointing" if out else "no_coverage")
    return {"status": status, "n_index_rows": len(rows), "n_radio_files": n_radio, "index_answer_full": len(rows) >= BL_ROW_LIMIT,
            "n_pointings": len(out), "n_pointings_in_beam": len(in_beam), "pointings": out}


def sample_radio_coverage(
    stars: list[dict[str, Any]], cache: RawCache, *, sleep: Callable[[float], Any] = time.sleep,
    progress: ProgressCallback | None = None, min_interval_s: float = MIN_REQUEST_INTERVAL_S,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Radio-pointing coverage of every star. Returns (per-star records in input order, summary).

    Cached answers are not rate-limited, so an interrupted census resumes quickly. A failed star is recorded as
    "failed: ..." and counted separately; it is neither covered nor uncovered.
    """
    records: list[dict[str, Any]] = []
    counts = {"in_beam": 0, "near_pointing": 0, "no_coverage": 0, "failed": 0}
    t_start = time.monotonic()
    for i, s in enumerate(stars):
        sid = str(s["source_id"])
        ra, dec = float(s["ra_deg"]), float(s["dec_deg"])
        t_req = time.monotonic()
        try:
            rows, man, hit = query_bl_files(ra, dec, cache, sleep=sleep)
            rec = beam_coverage(ra, dec, rows)
            rec["manifest_id"] = man.manifest_id if man is not None else None
        except ArchiveAccessError as exc:
            rec, hit = {"status": f"failed: {str(exc)[:200]}", "pointings": []}, False
        rec.update({"source_id": sid, "ra_deg": ra, "dec_deg": dec})
        key = "failed" if rec["status"].startswith("failed") else rec["status"]
        counts[key] += 1
        records.append(rec)
        if progress is not None:
            progress(i + 1, len(stars), sid, rec["status"])
        if not hit:
            wait = min_interval_s - (time.monotonic() - t_req)
            if wait > 0:
                sleep(wait)
    summary = {
        "n_input": len(stars), "n_in_beam": counts["in_beam"], "n_near_pointing": counts["near_pointing"],
        "n_no_coverage": counts["no_coverage"], "n_failed": counts["failed"],
        "in_beam_ids": [r["source_id"] for r in records if r["status"] == "in_beam"],
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
                rows, _, _ = query_bl_files(ra, dec, cache)
                scans = plan_star_scans(ra, dec, rows, cache)
                entry: dict[str, Any] = {"status": "ok", "scans": scans}
            except ArchiveAccessError as exc:
                entry = {"status": f"failed: {str(exc)[:200]}", "scans": []}
            complete = [s for s in entry["scans"] if s["complete"]]
            entry.update({"source_id": rec["source_id"], "ra_deg": ra, "dec_deg": dec, "n_on_scans": len(entry["scans"]),
                          "n_complete_scans": len(complete),
                          "smallest_complete_bytes": min((s["bytes_on_plus_offs"] for s in complete), default=None)})
            stars_out.append(entry)
            print(f"  plan {rec['source_id']}: {entry['n_on_scans']} ON scan(s), {entry['n_complete_scans']} with both OFF scans", file=sys.stderr)
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
