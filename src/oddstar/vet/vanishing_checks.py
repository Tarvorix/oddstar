"""Vanishing-source vetting checks (Section 6)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.io.physics import angular_separation_arcsec, propagate_position
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result


def check_proper_motion(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vd = star.get("vanishing_data", {})
    usno_epoch = float(vd.get("usno_epoch", 1953.4))
    usno_ra = float(vd.get("usno_ra_deg", star["ra_deg"]))
    usno_dec = float(vd.get("usno_dec_deg", star["dec_deg"]))
    pred_ra, pred_dec = propagate_position(
        float(star["ra_deg"]),
        float(star["dec_deg"]),
        float(star.get("pmra_mas_yr", 0.0)),
        float(star.get("pmdec_mas_yr", 0.0)),
        2016.0,
        usno_epoch,
    )
    pm_resid = angular_separation_arcsec(pred_ra, pred_dec, usno_ra, usno_dec)
    raw_sep = angular_separation_arcsec(float(star["ra_deg"]), float(star["dec_deg"]), usno_ra, usno_dec)
    is_fast_pm = (raw_sep > 3.0) and (pm_resid <= float(cfg.channels["vanishing"].get("max_pm_residual_arcsec", 1.5))) and bool(vd.get("present_in_gaia", True))
    return make_vet_result(
        out_dir,
        check="proper_motion",
        status="fail" if is_fast_pm else "pass",
        note=f"Raw epoch sep={raw_sep:.2f}\", PM-propagated residual={pm_resid:.2f}\"",
        metrics={"raw_sep_arcsec": raw_sep, "pm_residual_arcsec": pm_resid, "fast_pm_resolved": is_fast_pm},
    )


def check_plate_defect(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vd = star.get("vanishing_data", {})
    n_plates = int(vd.get("early_plates_count", 2))
    min_plates = int(cfg.channels["vanishing"].get("min_early_plates", 2))
    ok = n_plates >= min_plates
    return make_vet_result(
        out_dir,
        check="plate_defect",
        status="pass" if ok else "fail",
        note=f"Early-epoch detections on {n_plates} independent plates/bands (min {min_plates})",
        metrics={"early_plates_count": n_plates, "min_required": min_plates},
    )


def check_depth(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vd = star.get("vanishing_data", {})
    usno_mag = float(vd.get("usno_mag_r", 17.0))
    ps1_depth = float(vd.get("panstarrs_depth_mag", 21.8))
    margin = ps1_depth - usno_mag
    min_margin = float(cfg.channels["vanishing"].get("min_depth_margin_mag", 2.0))
    ok = margin >= min_margin
    return make_vet_result(
        out_dir,
        check="depth",
        status="pass" if ok else "fail",
        note=f"Later survey depth ({ps1_depth:.1f} mag) exceeds early mag ({usno_mag:.1f}) by {margin:.1f} mag",
        metrics={"usno_mag": usno_mag, "later_depth_mag": ps1_depth, "depth_margin_mag": margin},
    )


def check_known_variable_vanishing(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vd = star.get("vanishing_data", {})
    is_var = bool(vd.get("known_lpv", False)) or str(star.get("simbad_type", "")) in {"Mira", "LPV*", "CV*", "Flare*"}
    return make_vet_result(
        out_dir,
        check="known_variable",
        status="fail" if is_var else "pass",
        note="Not a catalogued Mira/LPV/flare variable" if not is_var else f"Explained by large-amplitude variable ({star.get('simbad_type')})",
        metrics={"known_lpv": bool(vd.get("known_lpv", False)), "simbad_type": star.get("simbad_type")},
    )


def check_sso_plate(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vd = star.get("vanishing_data", {})
    sso = bool(vd.get("sso_on_plate", False))
    return make_vet_result(
        out_dir,
        check="sso_plate",
        status="fail" if sso else "pass",
        note="No Solar System body on historical plate epoch" if not sso else "Asteroid matched at historical plate epoch",
        metrics={"sso_on_plate": sso},
    )


def check_blend_split(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vd = star.get("vanishing_data", {})
    n_split = int(vd.get("blend_split_sources", 1))
    bad = n_split > 1
    return make_vet_result(
        out_dir,
        check="blend_split",
        status="fail" if bad else "pass",
        note="Single unresolved component in modern high-res imaging" if not bad else f"Old plate blob resolved into {n_split} modern sources",
        metrics={"blend_split_sources": n_split},
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
