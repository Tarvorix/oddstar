"""Vanishing-source vetting checks (Section 6).

Input rule (R9, lessons L007/L008): a check returns `pass` or `fail` only from inputs that were actually obtained.

- The vanishing channel has NO live ingestion: no code path fetches a USNO-B1.0 / historical plate record,
  a Pan-STARRS depth, a plate-epoch Solar System ephemeris, modern high-resolution imaging or a variable-star
  catalogue cross-match for a live star. A live star (`live_enriched` key) counts such an input only when
  `vetting_archive.query_status[<query>] == "ok"`, which nothing writes today, so every check here returns
  `not_evaluated` for live stars and names the missing data product.
- Synthetic fixtures and hand-built dicts count an input only when its key is explicitly present and not None.
- The SIMBAD part of `known_variable` follows `oddstar.vet.common.simbad_lookup` (failed lookup ->
  `not_evaluated`, "no entry" -> a real negative answer).

Every evaluated result stores the input values it used in its metrics (copied into the committed evidence).
"""

from __future__ import annotations

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
    "panstarrs_depth_mag": ("panstarrs_depth", "Pan-STARRS local limiting depth"),
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
    """Fail when the plate-to-Gaia offset is a high-proper-motion star that moved (Gaia counterpart present)."""
    inputs = {k: _vd_input(star, k) for k in ("usno_epoch", "usno_ra_deg", "usno_dec_deg", "present_in_gaia")}
    gaia = {k: _star_input(star, k) for k in ("ra_deg", "dec_deg", "pmra_mas_yr", "pmdec_mas_yr")}
    missing = [k for k, v in {**inputs, **gaia}.items() if v is None]
    metrics: dict[str, Any] = {"input_values": {**inputs, **gaia}}
    if missing:
        return _not_evaluated(out_dir, "proper_motion", star, missing, metrics)
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


def check_plate_defect(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    n_raw = _vd_input(star, "early_plates_count")
    min_plates = int(cfg.channels["vanishing"]["min_early_plates"])
    metrics: dict[str, Any] = {"input_values": {"early_plates_count": n_raw}, "min_required": min_plates}
    if n_raw is None:
        return _not_evaluated(out_dir, "plate_defect", star, ["early_plates_count"], metrics)
    n_plates = int(n_raw)
    metrics["early_plates_count"] = n_plates
    return make_vet_result(
        out_dir,
        check="plate_defect",
        status="pass" if n_plates >= min_plates else "fail",
        note=f"Early-epoch detections on {n_plates} independent plates/bands (min {min_plates})",
        metrics=metrics,
    )


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
