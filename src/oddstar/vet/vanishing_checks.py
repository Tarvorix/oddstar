"""Vanishing-source vetting checks (Section 6).

Input rule (R9, lessons L007/L008): a check returns `pass` or `fail` only from inputs that were actually obtained.

- A live object (`live_enriched` key) counts an input only when `vetting_archive.query_status[<query>] == "ok"`.
  The plate-source survey of the vanishing search (`oddstar.io.plate_survey`) writes `usno_b1`,
  `modern_counterpart` and `panstarrs_depth` for a USNO-B1.0 source it compared with Gaia and Pan-STARRS;
  `vet --live` adds `plate_images` (what the digitized plates show at the position), which `plate_defect` needs. For a
  plate source without a Gaia counterpart, `proper_motion` passes only on the survey's record of how close any Gaia
  star's path over the plate era comes to the plate position (`nearest_gaia_track`). No code
  path fetches a plate-epoch Solar System ephemeris, modern high-resolution imaging or a variable-star catalogue
  cross-match, so `sso_plate`, `blend_split` and `known_variable` return `not_evaluated` for live objects and name
  the missing data product. A live sample star (not a plate source) has none of these inputs.
- Synthetic fixtures and hand-built dicts count an input only when its key is explicitly present and not None.
- The SIMBAD part of `known_variable` follows `oddstar.vet.common.simbad_lookup` (failed lookup ->
  `not_evaluated`, "no entry" -> a real negative answer).

Every evaluated result stores the input values it used in its metrics (copied into the committed evidence).
"""

from __future__ import annotations

import math

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.io.physics import angular_separation_arcsec, propagate_position
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result, simbad_lookup, simbad_note

GAIA_EPOCH_JYEAR = 2016.0

# vanishing_data key -> (query name in vetting_archive.query_status, data product named in not_evaluated notes)
_PRODUCTS: dict[str, tuple[str, str]] = {
    "usno_epoch": ("usno_b1", "USNO-B1.0 / historical plate record (epoch)"),
    "usno_ra_deg": ("usno_b1", "USNO-B1.0 / historical plate record (position)"),
    "usno_dec_deg": ("usno_b1", "USNO-B1.0 / historical plate record (position)"),
    "usno_mag_r": ("usno_b1", "USNO-B1.0 / historical plate record (magnitude)"),
    "early_plates_count": ("usno_b1", "USNO-B1.0 / historical plate detections (independent plate count)"),
    "present_in_gaia": ("modern_counterpart", "modern Gaia counterpart match at the plate position"),
    "nearest_gaia_track": ("modern_counterpart", "paths of the Gaia stars of the field over the plate era"),
    "panstarrs_depth_mag": ("panstarrs_depth", "Pan-STARRS local limiting depth"),
    "plate_images": ("plate_images", "digitized first-epoch plate images at the catalogued position"),
    "plate_neighbours": ("plate_neighbours", "first-epoch plate images at the nearby Pan-STARRS objects"),
    "plate_first_place": ("plate_first_place", "Gaia and Pan-STARRS objects at the place where the first-epoch plates show the source"),
    "sso_on_plate": ("plate_epoch_ephemeris", "Solar System ephemeris at the plate epoch"),
    "blend_split_sources": ("highres_imaging", "modern high-resolution imaging (blend split)"),
    "known_lpv": ("variable_star_catalogue", "variable-star catalogue cross-match (LPV/Mira)"),
}
NO_LIVE_INGESTION = "the vanishing channel has no live ingestion for this product"


def _is_live(star: dict[str, Any]) -> bool:
    return "live_enriched" in star


def _query_status(star: dict[str, Any], query: str) -> str:
    return ((star.get("vetting_archive") or {}).get("query_status") or {}).get(query) or "not run"


def _vd_input(star: dict[str, Any], key: str) -> Any:
    """Explicit `vanishing_data[key]`; for live stars only when its query is recorded as "ok"."""
    if _is_live(star) and _query_status(star, _PRODUCTS[key][0]) != "ok":
        return None
    return (star.get("vanishing_data") or {}).get(key)


def _star_input(star: dict[str, Any], key: str) -> float | None:
    val = star.get(key)
    return None if val is None else float(val)


def _missing_note(star: dict[str, Any], keys: list[str]) -> str:
    parts = []
    for k in keys:
        if k in _PRODUCTS:
            query, product = _PRODUCTS[k]
            txt = f"{product} [vanishing_data.{k}"
            if _is_live(star):
                txt += f"; query_status {query}={_query_status(star, query)}; {NO_LIVE_INGESTION}"
            parts.append(txt + "]")
        else:
            parts.append(f"Gaia {k}")
    return "Missing input: " + "; ".join(parts)


def _not_evaluated(out_dir: Path, check: str, star: dict[str, Any], missing: list[str], metrics: dict[str, Any]) -> VettingResult:
    m = dict(metrics, missing_inputs=missing, live_star=_is_live(star))
    return make_vet_result(out_dir, check=check, status="not_evaluated", note=_missing_note(star, missing), metrics=m)


def check_proper_motion(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Fail when the plate-to-Gaia offset is a high-proper-motion star that moved (Gaia counterpart present).

    A plate source of the survey without a Gaia counterpart has no star to move. It passes only when the survey
    recorded that no Gaia star's path over the plate era comes within the match radius of the plate position. A live
    plate source also fails when a present-day object is at the place where the first-epoch plates show it
    (`first_place_counterpart`: a false catalogue proper motion) or when a nearby Pan-STARRS object is the source,
    moved (`moved_star_verdict`), and passes only when both tests could be made and found nothing.
    """
    inputs = {k: _vd_input(star, k) for k in ("usno_epoch", "usno_ra_deg", "usno_dec_deg", "present_in_gaia")}
    gaia = {k: _star_input(star, k) for k in ("ra_deg", "dec_deg", "pmra_mas_yr", "pmdec_mas_yr")}
    missing = [k for k, v in {**inputs, **gaia}.items() if v is None]
    metrics: dict[str, Any] = {"input_values": {**inputs, **gaia}}
    if missing:
        return _not_evaluated(out_dir, "proper_motion", star, missing, metrics)
    if star.get("plate_source") and not bool(inputs["present_in_gaia"]):
        track = _vd_input(star, "nearest_gaia_track")
        metrics["input_values"]["nearest_gaia_track"] = track
        if not isinstance(track, dict):
            return _not_evaluated(out_dir, "proper_motion", star, ["nearest_gaia_track"], metrics)
        match_rad = float(cfg.channels["vanishing"]["match_radius_arcsec"])
        sep = track.get("sep_arcsec")
        moved_star = sep is not None and float(sep) <= match_rad
        metrics.update({"nearest_track_sep_arcsec": sep, "match_radius_arcsec": match_rad, "n_gaia_stars": track.get("n_gaia_stars"),
                        "fast_pm_resolved": moved_star})
        era = track.get("plate_era") or ["?", "?"]
        where = (f"nearest path {float(sep):.2f}\" (Gaia DR3 {track.get('source_id')}, at {track.get('epoch')})" if sep is not None
                 else "no Gaia star in the field")
        gaia_note = (f"No Gaia counterpart; paths of {track.get('n_gaia_stars')} Gaia star(s) over {era[0]}-{era[1]} against the plate position: "
                     f"{where} (match radius {match_rad:.1f}\")")
        if moved_star or not _is_live(star):
            return make_vet_result(out_dir, check="proper_motion", status="fail" if moved_star else "pass", note=gaia_note, metrics=metrics)
        # The catalogue position can be wrong (a false catalogue proper motion): the survey's comparison again, at the
        # place where the first-epoch plates show the source
        still = first_place_counterpart(star, cfg)
        metrics["first_place_test"] = still
        if still["verdict"]:
            return make_vet_result(out_dir, check="proper_motion", status="fail", metrics=metrics, note=f"{gaia_note}. {still['note']}")
        # A star too faint for Gaia can still have moved: the Pan-STARRS objects nearby against the first-epoch plates
        moved = moved_star_verdict(star, cfg)
        metrics["moved_star_test"] = moved
        metrics["input_values"]["plate_neighbours"] = moved.get("neighbours")
        if moved["verdict"]:
            return make_vet_result(out_dir, check="proper_motion", status="fail", metrics=metrics, note=f"{gaia_note}. {moved['note']}")
        if moved["verdict"] is None:
            metrics.update(missing_inputs=moved["missing"], live_star=True)
            return make_vet_result(out_dir, check="proper_motion", status="not_evaluated", metrics=metrics,
                                   note=f"{gaia_note}. Not settled for a star too faint for Gaia: {moved['note']}")
        if still["verdict"] is None:
            metrics.update(missing_inputs=still["missing"], live_star=True)
            return make_vet_result(out_dir, check="proper_motion", status="not_evaluated", metrics=metrics,
                                   note=f"{gaia_note}. Not settled at the place where the plates show the source: {still['note']}")
        return make_vet_result(out_dir, check="proper_motion", status="pass", metrics=metrics,
                               note=f"{gaia_note}. {still['note']}. {moved['note']}")
    usno_epoch = float(inputs["usno_epoch"])
    usno_ra, usno_dec = float(inputs["usno_ra_deg"]), float(inputs["usno_dec_deg"])
    pred_ra, pred_dec = propagate_position(
        gaia["ra_deg"], gaia["dec_deg"], gaia["pmra_mas_yr"], gaia["pmdec_mas_yr"], GAIA_EPOCH_JYEAR, usno_epoch
    )
    pm_resid = angular_separation_arcsec(pred_ra, pred_dec, usno_ra, usno_dec)
    raw_sep = angular_separation_arcsec(gaia["ra_deg"], gaia["dec_deg"], usno_ra, usno_dec)
    min_raw = float(cfg.vetting["vanishing_pm_min_raw_sep_arcsec"])
    max_resid = float(cfg.channels["vanishing"]["max_pm_residual_arcsec"])
    present = bool(inputs["present_in_gaia"])
    is_fast_pm = raw_sep > min_raw and pm_resid <= max_resid and present
    metrics.update({
        "raw_sep_arcsec": raw_sep, "pm_residual_arcsec": pm_resid, "fast_pm_resolved": is_fast_pm,
        "predicted_plate_epoch_ra_deg": pred_ra, "predicted_plate_epoch_dec_deg": pred_dec,
        "min_raw_sep_arcsec": min_raw, "max_pm_residual_arcsec": max_resid,
    })
    return make_vet_result(
        out_dir,
        check="proper_motion",
        status="fail" if is_fast_pm else "pass",
        note=(f"Plate epoch {usno_epoch:.1f}: raw plate-Gaia sep={raw_sep:.2f}\", PM-propagated residual={pm_resid:.2f}\" "
              f"(PM {gaia['pmra_mas_yr']:.1f}, {gaia['pmdec_mas_yr']:.1f} mas/yr; Gaia counterpart present={present})"),
        metrics=metrics,
    )


def plate_image_verdict(row: dict[str, Any], cfg: PipelineConfig) -> bool | None:
    """Does one plate image show a point of light at the catalogued position? None when the row cannot say.

    A detection is a local maximum within `plate_image_max_offset_arcsec` of the position whose aperture signal is at
    least `plate_image_min_snr` times the measured scatter of the same aperture around it. An image that could not be
    registered can show the source but cannot show its absence.
    """
    vc = cfg.channels["vanishing"]
    if row.get("status") != "ok":
        return None
    peak = row.get("peak")
    seen = bool(isinstance(peak, dict) and peak.get("snr") is not None and float(peak["snr"]) >= float(vc["plate_image_min_snr"])
                and float(peak["offset_arcsec"]) <= float(vc["plate_image_max_offset_arcsec"]))
    if not seen and row.get("registered") is False:
        # Where the position falls on an image that could not be lined up on Gaia stars is uncertain by more than the
        # offset limit (first-epoch red plates are off by up to 4 arcsec), so an empty aperture there says nothing
        return None
    return seen


def first_epoch_seen(rows: list[dict[str, Any]], cfg: PipelineConfig) -> tuple[int, int]:
    """(first-epoch images that show a point of light at the measured position, first-epoch images that could say)."""
    before = float(cfg.channels["vanishing"]["plate_image_first_epoch_before"])
    said = [v for r in rows if isinstance(r, dict) and r.get("plate_epoch") is not None and float(r["plate_epoch"]) < before
            for v in [plate_image_verdict(r, cfg)] if v is not None]
    return sum(1 for v in said if v), len(said)


PANSTARRS_EPOCH = 2012.0  # mean epoch of the Pan-STARRS1 3-pi survey positions, to about two years


def _peak_offset_from(row: dict[str, Any], ra: float, dec: float) -> tuple[float, float] | None:
    """Where a row's peak is on the sky, as arcsec east and north of (ra, dec)."""
    peak = row.get("peak")
    if not isinstance(peak, dict) or row.get("ra_measured_deg") is None or peak.get("dx_arcsec") is None:
        return None
    cos_dec = max(1e-6, math.cos(math.radians(dec)))
    east = ((float(row["ra_measured_deg"]) - ra + 180.0) % 360.0 - 180.0) * cos_dec * 3600.0 - float(peak["dx_arcsec"])  # east is towards smaller x
    north = (float(row["dec_measured_deg"]) - dec) * 3600.0 + float(peak["dy_arcsec"])
    return east, north


def first_place_counterpart(star: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    """Is a present-day object at the place where the first-epoch plates show the source?

    The survey looks for a Gaia or Pan-STARRS counterpart at the catalogue position. USNO-B1.0 gives some sources a
    false proper motion, by joining detections of different objects into one track; the catalogue position is then
    several arcsec from the object the first plates recorded, and the survey finds nothing there although the object
    never left. This repeats the survey's own comparison, with the survey's `match_radius_arcsec`, at the mean sky
    place of the first-epoch plate images that show the source: a Gaia star that was there at the plate's date
    (moved back with its Gaia proper motion), or a Pan-STARRS object there now.
    `verdict` is True (a counterpart: the source did not vanish), False (none, or no first-epoch image shows the
    source, so there is no place to compare at) or None (open; `missing` names the input).
    """
    match = float(cfg.channels["vanishing"]["match_radius_arcsec"])
    before = float(cfg.channels["vanishing"]["plate_image_first_epoch_before"])
    images = _vd_input(star, "plate_images")
    out: dict[str, Any] = {"match_radius_arcsec": match, "place": None, "counterparts": []}
    if not isinstance(images, dict):
        return dict(out, verdict=None, missing=["plate_images"], note=_missing_note(star, ["plate_images"]))
    shown = [r for r in images.get("rows") or [] if isinstance(r, dict) and r.get("plate_epoch") is not None
             and float(r["plate_epoch"]) < before and plate_image_verdict(r, cfg)]
    if not shown:
        return dict(out, verdict=False, missing=[], note="No first-epoch plate image shows the source, so there is no measured place to compare at")
    place = _vd_input(star, "plate_first_place")
    if not isinstance(place, dict) or not isinstance(place.get("position"), list):
        return dict(out, verdict=None, missing=["plate_first_place"], note=_missing_note(star, ["plate_first_place"]))
    ra0, dec0 = float(place["position"][0]), float(place["position"][1])
    offs = [o for o in (_peak_offset_from(r, ra0, dec0) for r in shown) if o is not None]
    if not offs:
        return dict(out, verdict=None, missing=["plate_images"], note="The first-epoch rows that show the source carry no peak position")
    e1, n1 = sum(o[0] for o in offs) / len(offs), sum(o[1] for o in offs) / len(offs)
    t1 = sum(float(r["plate_epoch"]) for r in shown) / len(shown)
    cos_dec = max(1e-6, math.cos(math.radians(dec0)))

    def _dist(ra: float, dec: float) -> float:
        east = ((ra - ra0 + 180.0) % 360.0 - 180.0) * cos_dec * 3600.0
        return math.hypot(east - e1, (dec - dec0) * 3600.0 - n1)

    found: list[dict[str, Any]] = []
    for g in place.get("gaia") or []:
        d = _dist(float(g["ra_then_deg"]), float(g["dec_then_deg"]))
        if d <= match:
            found.append({"catalogue": "Gaia DR3", "id": str(g.get("source_id")), "sep_arcsec": round(d, 2), "mag": g.get("g_mag"), "band": "G",
                          "pm_mas_yr": round(math.hypot(float(g.get("pmra_mas_yr") or 0.0), float(g.get("pmdec_mas_yr") or 0.0)), 1)})
    for p in place.get("ps1") or []:
        d = _dist(float(p["ra_deg"]), float(p["dec_deg"]))
        if d <= match:
            found.append({"catalogue": "Pan-STARRS", "id": str(p.get("objID")), "sep_arcsec": round(d, 2), "mag": p.get("r_mag"), "band": "r",
                          "n_detections": p.get("n_detections")})
    found.sort(key=lambda f: (f["sep_arcsec"], f["id"]))
    out.update(place={"east_arcsec": round(e1, 2), "north_arcsec": round(n1, 2), "of_position": [ra0, dec0], "epoch": round(t1, 1),
                      "n_images": len(offs)},
               counterparts=found, n_gaia_in_cone=len(place.get("gaia") or []), n_ps1_in_cone=len(place.get("ps1") or []))
    where = f"the place where {len(offs)} first-epoch plate image(s) ({t1:.0f}) show the source"
    if not found:
        return dict(out, verdict=False, missing=[], note=f"No Gaia star was, and no Pan-STARRS object is, within {match:.1f}\" of {where}")
    f = found[0]
    mag = "" if f.get("mag") is None else f", {f['band']} {float(f['mag']):.2f}"
    when = "was" if f["catalogue"] == "Gaia DR3" else "is"
    return dict(out, verdict=True, missing=[], note=(
        f"{f['catalogue']} {f['id']}{mag} {when} {f['sep_arcsec']:.2f}\" from {where}: the object the plates recorded is still there, "
        "and it is the catalogue position that is off"))


def plate_track(rows: list[dict[str, Any]], ra: float, dec: float, cfg: PipelineConfig) -> dict[str, Any]:
    """Three epochs on a line: did the plate source travel towards a present-day position (ra, dec)?

    The first-epoch images that show the source give where it was (P1, t1); the object is where it would be now (P3,
    2012). On each later plate that shows a point of light near the source (registered, `plate_image_min_snr`), that
    light is either where a star on its way from P1 to P3 would be at the plate's date (`d_line`), or still at P1
    (`d_stay`). `kind` is "moved" when at least `plate_track_min_images` later images are within
    `plate_track_max_residual_arcsec` of the line and fewer are at P1, "stayed" for the reverse, "too_close" when
    P1 and P3 are nearer than the match radius (the two cases then look alike), else "open".
    On controls (2026-11-18): 43 of 57 known moving stars "moved" and none "stayed"; 73 of 79 still stars with a
    neighbour "stayed" and none "moved".
    """
    vc = cfg.channels["vanishing"]
    before, tol, need = float(vc["plate_image_first_epoch_before"]), float(vc["plate_track_max_residual_arcsec"]), int(vc["plate_track_min_images"])
    first = [(o, float(r["plate_epoch"])) for r in rows if isinstance(r, dict) and r.get("plate_epoch") is not None
             and float(r["plate_epoch"]) < before and plate_image_verdict(r, cfg) for o in [_peak_offset_from(r, ra, dec)] if o is not None]
    out: dict[str, Any] = {"n_first_epoch": len(first), "max_residual_arcsec": tol, "min_images": need, "later": []}
    if not first:
        return dict(out, kind="open", why="no first-epoch image shows the source")
    e1, n1 = sum(o[0] for o, _ in first) / len(first), sum(o[1] for o, _ in first) / len(first)
    t1 = sum(t for _, t in first) / len(first)
    length = math.hypot(e1, n1)
    out.update(first_epoch=round(t1, 1), first_to_now_arcsec=round(length, 2))
    n_line = n_stay = 0
    for r in rows:
        if not isinstance(r, dict) or r.get("status") != "ok" or r.get("plate_epoch") is None or float(r["plate_epoch"]) < before:
            continue
        peak = r.get("peak")
        if r.get("registered") is not True or not isinstance(peak, dict) or peak.get("snr") is None \
                or float(peak["snr"]) < float(vc["plate_image_min_snr"]):
            continue
        o = _peak_offset_from(r, ra, dec)
        if o is None:
            continue
        frac = (float(r["plate_epoch"]) - t1) / (PANSTARRS_EPOCH - t1)
        d_line = math.hypot(o[0] - e1 * (1.0 - frac), o[1] - n1 * (1.0 - frac))  # the object itself is at (0, 0)
        d_stay = math.hypot(o[0] - e1, o[1] - n1)
        n_line += d_line <= tol
        n_stay += d_stay <= tol
        out["later"].append({"survey": r.get("survey"), "plate_epoch": r.get("plate_epoch"), "d_line_arcsec": round(d_line, 2),
                             "d_stay_arcsec": round(d_stay, 2)})
    out.update(n_on_line=int(n_line), n_stayed=int(n_stay))
    if length < float(vc["match_radius_arcsec"]):
        return dict(out, kind="too_close", why=f"the source on the plates is only {length:.2f}\" from the object")
    if n_line >= need and n_stay < need:
        return dict(out, kind="moved")
    if n_stay >= need and n_line < need:
        return dict(out, kind="stayed")
    return dict(out, kind="open", why=f"{len(out['later'])} later image(s) usable: {int(n_line)} on the line, {int(n_stay)} at the first-epoch place")


def moved_star_verdict(star: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    """Is a nearby Pan-STARRS object the plate source itself, having moved since the plates were taken?

    A star too faint for Gaia can move several arcsec in sixty years. Pan-STARRS then holds an object of about the
    plate magnitude (within `plate_moved_max_dmag`) between the match radius and `plate_moved_radius_arcsec` of the
    plate position. Such an object is another star when a first-epoch plate shows it where it is now (none of 60
    known moving stars was; control run of 2026-11-18). Otherwise the later plates decide (`plate_track`): the
    source on its way to the object means the same star, moved; the source still in place means another star.
    Empty first-epoch plates at the object's position are not taken as proof of motion: that rule called 10 of 100
    still stars with a neighbour moved.
    `verdict` is True (moved), False (no object is the source, moved) or None (open; `missing` names the input).
    """
    vc = cfg.channels["vanishing"]
    min_plates, max_dmag = int(vc["min_early_plates"]), float(vc["plate_moved_max_dmag"])
    images, nb = _vd_input(star, "plate_images"), _vd_input(star, "plate_neighbours")
    r1 = _vd_input(star, "usno_mag_r")
    if not isinstance(images, dict):
        return {"verdict": None, "missing": ["plate_images"], "note": _missing_note(star, ["plate_images"]), "neighbours": None}
    rows = [r for r in images.get("rows") or [] if isinstance(r, dict)]
    n_seen, n_said = first_epoch_seen(rows, cfg)
    out: dict[str, Any] = {"source_first_epoch_seen": n_seen, "source_first_epoch_images": n_said, "plate_moved_max_dmag": max_dmag, "neighbours": None}
    if n_said >= min_plates and n_seen == 0:
        # Nothing on the plates to have moved; what the catalogue entry is, is plate_defect's question
        return dict(out, verdict=False, missing=[], note=f"None of the {n_said} first-epoch plate images shows the source, so there is no "
                                                         "plate image of a star that could have moved")
    if not isinstance(nb, dict) or r1 is None:
        missing = ["plate_neighbours"] if not isinstance(nb, dict) else ["usno_mag_r"]
        return dict(out, verdict=None, missing=missing, note=_missing_note(star, missing))
    table, open_, moved = [], [], []
    for n in nb.get("neighbours") or []:
        o_seen, o_said = first_epoch_seen(n.get("rows") or [], cfg)
        dmag = None if n.get("r_mag") is None else round(float(n["r_mag"]) - float(r1), 2)
        track = None
        if dmag is not None and abs(dmag) > max_dmag:
            kind = "other_brightness"
        elif o_seen >= 1:
            kind = "on_the_plates_where_it_is_now"
        else:
            kind = "could_be_the_source_moved"
            if n.get("ra_deg") is not None and n.get("dec_deg") is not None:
                track = plate_track(rows, float(n["ra_deg"]), float(n["dec_deg"]), cfg)
                kind = {"moved": "the_source_moved", "stayed": "source_stayed_in_place"}.get(track["kind"], kind)
        table.append({"objID": n.get("objID"), "sep_arcsec": n.get("sep_arcsec"), "r_mag": n.get("r_mag"), "dmag": dmag,
                      "n_detections": n.get("n_detections"), "first_epoch_seen": o_seen, "first_epoch_images": o_said, "kind": kind,
                      "track": track})
        if kind == "could_be_the_source_moved":
            open_.append(table[-1])
        elif kind == "the_source_moved":
            moved.append(table[-1])
    out["neighbours"] = table
    out["n_objects_in_radius"] = nb.get("n_objects_in_radius")
    if moved:
        t = moved[0]
        tr = t["track"]
        return dict(out, verdict=True, missing=[], note=(
            f"Pan-STARRS object {t['objID']} is {tr['first_to_now_arcsec']}\" from where the first-epoch plates ({tr['first_epoch']:.0f}) show the "
            f"source (about {1000.0 * tr['first_to_now_arcsec'] / (PANSTARRS_EPOCH - tr['first_epoch']):.0f} mas/yr), is {t['dmag']:+.2f} mag from "
            f"the plate magnitude, and {tr['n_on_line']} later plate image(s) show the source on its way there, within "
            f"{tr['max_residual_arcsec']:g}\" of the straight line: the same star, moved"))
    if open_:
        t = open_[0]
        mag = "has no r magnitude" if t["dmag"] is None else f"is {t['dmag']:+.2f} mag from the plate magnitude"
        why = f"; later plates: {t['track']['why']}" if t.get("track") and t["track"].get("why") else ""
        return dict(out, verdict=None, missing=["plate_neighbours"], note=(
            f"Pan-STARRS object {t['objID']} lies {t['sep_arcsec']}\" from the plate position, {mag}, and no first-epoch plate shows it "
            f"where it is now ({t['first_epoch_seen']} of {t['first_epoch_images']} images): it may be the plate source, moved{why}"))
    if int(nb.get("n_objects_in_radius") or 0) > len(table):
        return dict(out, verdict=None, missing=["plate_neighbours"],
                    note=f"{nb.get('n_objects_in_radius')} Pan-STARRS objects in the radius, only {len(table)} measured")
    return dict(out, verdict=False, missing=[], note=(
        f"Of {len(table)} Pan-STARRS object(s) between {nb.get('inner_radius_arcsec')}\" and {nb.get('radius_arcsec')}\" none is the "
        "source, moved: " + ("; ".join(f"{t['objID']} at {t['sep_arcsec']}\" {t['kind']}" for t in table) or "none in the radius")))


def check_plate_defect(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Is the plate source real? Fail on a source the first-epoch plates do not hold.

    The catalogue must list it on `min_early_plates` first-epoch plates. For a plate source of the live survey that
    is not enough: USNO-B1.0 lists sources that are pieces of bright-star halos, shared plate flaws and mismatched
    detections of moving stars. Such a source passes only when the digitized first-epoch plates themselves show a
    point of light at its position on `min_early_plates` plates (`oddstar.io.plate_images`). It fails when
    `min_early_plates` first-epoch images could say and none shows it, provided the catalogue magnitude is no
    fainter than `plate_image_judge_max_r_mag`. It is not evaluated when the images disagree, when too few could
    say (not read, or not registered and empty), or when the source is too faint for a missing image to count.
    """
    vc = cfg.channels["vanishing"]
    n_raw = _vd_input(star, "early_plates_count")
    min_plates = int(vc["min_early_plates"])
    metrics: dict[str, Any] = {"input_values": {"early_plates_count": n_raw}, "min_required": min_plates}
    if n_raw is None:
        return _not_evaluated(out_dir, "plate_defect", star, ["early_plates_count"], metrics)
    n_plates = int(n_raw)
    metrics["early_plates_count"] = n_plates
    if n_plates < min_plates or not star.get("plate_source"):
        return make_vet_result(
            out_dir,
            check="plate_defect",
            status="pass" if n_plates >= min_plates else "fail",
            note=f"Early-epoch detections on {n_plates} independent plates/bands (min {min_plates})",
            metrics=metrics,
        )
    images = _vd_input(star, "plate_images")
    if not isinstance(images, dict):
        return _not_evaluated(out_dir, "plate_defect", star, ["plate_images"], metrics)
    first_before = float(vc["plate_image_first_epoch_before"])
    rows = [r for r in images.get("rows") or [] if isinstance(r, dict)]
    table = []
    for r in rows:
        peak = r.get("peak") or {}
        table.append({"survey": r.get("survey"), "survey_returned": r.get("survey_returned"), "plate_id": r.get("plate_id"),
                      "plate_epoch": r.get("plate_epoch"), "status": r.get("status"), "snr_at_position": r.get("snr_at_position"),
                      "peak_snr": peak.get("snr"), "peak_offset_arcsec": peak.get("offset_arcsec"), "peak_fwhm_arcsec": peak.get("fwhm_arcsec"),
                      "field_fwhm_median_arcsec": r.get("field_fwhm_median_arcsec"), "detected": plate_image_verdict(r, cfg)})
    first = [t for t in table if t["plate_epoch"] is not None and float(t["plate_epoch"]) < first_before and t["detected"] is not None]
    later = [t for t in table if t["plate_epoch"] is not None and float(t["plate_epoch"]) >= first_before and t["detected"] is not None]
    n_first_seen = sum(1 for t in first if t["detected"])
    metrics.update({"input_values": {"early_plates_count": n_raw, "plate_images": table,
                                     "measurement_version": images.get("measurement_version"), "errors": images.get("errors") or {}},
                    "plate_image_min_snr": float(vc["plate_image_min_snr"]), "plate_image_max_offset_arcsec": float(vc["plate_image_max_offset_arcsec"]),
                    "plate_image_first_epoch_before": first_before, "n_first_epoch_images": len(first), "n_first_epoch_seen": n_first_seen,
                    "n_later_images": len(later), "n_later_seen": sum(1 for t in later if t["detected"])})
    seen_txt = "; ".join(f"{t['survey_returned'] or t['survey']} {t['plate_epoch']:.1f}: " + (
        f"S/N {t['peak_snr']} at {t['peak_offset_arcsec']}\"" if t["peak_snr"] is not None else "no peak") for t in first) or "none readable"
    later_txt = f"; later plates: seen on {metrics['n_later_seen']} of {len(later)}" if later else ""
    limits = f"(limits: S/N {metrics['plate_image_min_snr']:g} within {metrics['plate_image_max_offset_arcsec']:g}\")"
    if n_first_seen >= min_plates:
        return make_vet_result(out_dir, check="plate_defect", status="pass", metrics=metrics,
                               note=f"Catalogued on {n_plates} first-epoch plates and the plate images show a point of light at the position on "
                                    f"{n_first_seen} of {len(first)} {limits}: {seen_txt}{later_txt}")
    if len(first) < min_plates:
        metrics.update(missing_inputs=["plate_images"], live_star=_is_live(star))
        return make_vet_result(out_dir, check="plate_defect", status="not_evaluated", metrics=metrics,
                               note=f"Only {len(first)} first-epoch plate image(s) could be measured at the position (need {min_plates}): {seen_txt}")
    r1, judge = _vd_input(star, "usno_mag_r"), float(vc["plate_image_judge_max_r_mag"])
    metrics.update({"usno_mag_r": r1, "plate_image_judge_max_r_mag": judge})
    if n_first_seen > 0:
        # Ordinary stars show on one first-epoch image and not the other about one time in eight; blank sky never did
        metrics.update(missing_inputs=["plate_images"], live_star=_is_live(star))
        return make_vet_result(out_dir, check="plate_defect", status="not_evaluated", metrics=metrics,
                               note=f"The first-epoch plate images disagree: a point of light at the position on {n_first_seen} of {len(first)} "
                                    f"{limits}: {seen_txt}{later_txt}")
    if r1 is None or float(r1) > judge:
        # Ordinary stars this faint are themselves often missed on these plates, so a missing image says nothing
        metrics.update(missing_inputs=["plate_images"], live_star=_is_live(star))
        return make_vet_result(out_dir, check="plate_defect", status="not_evaluated", metrics=metrics,
                               note=f"The plate images show a point of light at the position on {n_first_seen} of {len(first)} first-epoch "
                                    f"image(s) {limits}, but at plate magnitude {r1} a missing image is not evidence (judged only to {judge:g}): "
                                    f"{seen_txt}{later_txt}")
    return make_vet_result(out_dir, check="plate_defect", status="fail", metrics=metrics,
                           note=f"Catalogued on {n_plates} first-epoch plates, but none of the {len(first)} first-epoch plate images shows a "
                                f"point of light at the position {limits}: {seen_txt}{later_txt}. A catalogue entry the plates do not hold")


def check_depth(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    inputs = {k: _vd_input(star, k) for k in ("usno_mag_r", "panstarrs_depth_mag")}
    metrics: dict[str, Any] = {"input_values": inputs}
    missing = [k for k, v in inputs.items() if v is None]
    if missing:
        return _not_evaluated(out_dir, "depth", star, missing, metrics)
    usno_mag, ps1_depth = float(inputs["usno_mag_r"]), float(inputs["panstarrs_depth_mag"])
    margin = ps1_depth - usno_mag
    min_margin = float(cfg.channels["vanishing"]["min_depth_margin_mag"])
    metrics.update({"usno_mag": usno_mag, "later_depth_mag": ps1_depth, "depth_margin_mag": margin, "min_depth_margin_mag": min_margin})
    return make_vet_result(
        out_dir,
        check="depth",
        status="pass" if margin >= min_margin else "fail",
        note=f"Later survey depth ({ps1_depth:.1f} mag) exceeds early mag ({usno_mag:.1f}) by {margin:.1f} mag (min {min_margin:.1f})",
        metrics=metrics,
    )


def check_known_variable_vanishing(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Fail on a catalogued large-amplitude variable; pass needs both the SIMBAD lookup and the LPV catalogue flag."""
    var_types = {str(t) for t in cfg.channels["vanishing"]["vetting_variable_simbad_types"]}
    sim = simbad_lookup(star)
    lpv = _vd_input(star, "known_lpv")
    simbad_var = sim["evaluated"] and sim["simbad_type"] in var_types
    metrics: dict[str, Any] = {
        "input_values": {"known_lpv": lpv, "simbad_type": sim["simbad_type"], "simbad_query_status": sim["query_status"], "simbad_entry": sim["entry"]},
        "known_lpv": None if lpv is None else bool(lpv), "simbad_type": sim["simbad_type"], "variable_simbad_types": sorted(var_types),
    }
    if (lpv is not None and bool(lpv)) or simbad_var:
        why = []
        if lpv is not None and bool(lpv):
            why.append("variable-star catalogue lists an LPV/Mira")
        if simbad_var:
            why.append(f"SIMBAD type {sim['simbad_type']}")
        return make_vet_result(out_dir, check="known_variable", status="fail",
                               note="Explained by large-amplitude variable: " + "; ".join(why), metrics=metrics)
    if lpv is None:
        return _not_evaluated(out_dir, "known_variable", star, ["known_lpv"], metrics)
    if not sim["evaluated"]:
        return make_vet_result(out_dir, check="known_variable", status="not_evaluated",
                               note=f"Variable-star catalogue flag clear, but {simbad_note(sim)}", metrics=metrics)
    return make_vet_result(
        out_dir,
        check="known_variable",
        status="pass",
        note=f"Not a catalogued Mira/LPV/flare variable ({simbad_note(sim)}; variable-star catalogue flag clear)",
        metrics=metrics,
    )


def check_sso_plate(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    sso = _vd_input(star, "sso_on_plate")
    metrics: dict[str, Any] = {"input_values": {"sso_on_plate": sso}}
    if sso is None:
        return _not_evaluated(out_dir, "sso_plate", star, ["sso_on_plate"], metrics)
    hit = bool(sso)
    metrics["sso_on_plate"] = hit
    return make_vet_result(
        out_dir,
        check="sso_plate",
        status="fail" if hit else "pass",
        note="Asteroid matched at historical plate epoch" if hit else "No Solar System body on historical plate epoch",
        metrics=metrics,
    )


def check_blend_split(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    n_raw = _vd_input(star, "blend_split_sources")
    metrics: dict[str, Any] = {"input_values": {"blend_split_sources": n_raw}}
    if n_raw is None:
        return _not_evaluated(out_dir, "blend_split", star, ["blend_split_sources"], metrics)
    n_split = int(n_raw)
    metrics["blend_split_sources"] = n_split
    return make_vet_result(
        out_dir,
        check="blend_split",
        status="fail" if n_split > 1 else "pass",
        note=f"Old plate blob resolved into {n_split} modern sources" if n_split > 1 else "Single unresolved component in modern high-res imaging",
        metrics=metrics,
    )


def run_vanishing_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [
        check_proper_motion,
        check_plate_defect,
        check_depth,
        check_known_variable_vanishing,
        check_sso_plate,
        check_blend_split,
    ]
    return [fn(star, cfg, out_dir) for fn in checks]
