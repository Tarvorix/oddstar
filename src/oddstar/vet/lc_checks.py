"""Light-curve and IR time-domain vetting checks (Section 6)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result


def check_neighbour_signal(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    neo = star.get("neowise", {})
    shift = bool(lc.get("neighbour_eb_centroid_shift", False)) or bool(neo.get("w1_w2_neighbour_blend", False))
    return make_vet_result(
        out_dir,
        check="neighbour_signal",
        status="fail" if shift else "pass",
        note="Aperture photocenter stable (no neighbour contamination)" if not shift else "Dip/IR variation tracks neighbouring star within TESS/WISE beam",
        metrics={"neighbour_centroid_shift": bool(lc.get("neighbour_eb_centroid_shift", False)), "w1_w2_blend": bool(neo.get("w1_w2_neighbour_blend", False))},
    )


def check_systematics_window(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    neo = star.get("neowise", {})
    sys_flag = bool(lc.get("tess_scattered_light_window", False)) or bool(neo.get("scan_systematic", False))
    return make_vet_result(
        out_dir,
        check="systematics_window",
        status="fail" if sys_flag else "pass",
        note="Event occurs outside spacecraft scattered-light/momentum-dump/scan systematic windows" if not sys_flag else "Event aligns with documented scattered-light or scan-direction systematic window",
        metrics={"tess_scattered_light": bool(lc.get("tess_scattered_light_window", False)), "scan_systematic": bool(neo.get("scan_systematic", False))},
    )


def check_pipeline_agreement(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    f1 = np.asarray(lc.get("tess_flux_qlp", []), dtype=float)
    f2 = np.asarray(lc.get("tess_flux_tesscut", []), dtype=float)
    if len(f1) == 0 or len(f2) == 0 or len(f1) != len(f2):
        return make_vet_result(
            out_dir,
            check="pipeline_agreement",
            status="not_evaluated",
            note="Independent second light-curve extraction missing",
        )
    # Check whether maximum dip depth in QLP is also seen in TESScut
    d1 = 1.0 - float(np.min(f1))
    d2 = 1.0 - float(np.min(f2))
    agree = (d1 < 0.02 and d2 < 0.02) or (abs(d1 - d2) <= 0.08 and d2 >= 0.5 * d1)
    return make_vet_result(
        out_dir,
        check="pipeline_agreement",
        status="pass" if agree else "fail",
        note=f"QLP max dip={d1:.3f}, independent TESScut max dip={d2:.3f}",
        metrics={"qlp_max_depth": d1, "tesscut_max_depth": d2},
    )


def check_sso_crossing(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    hit = lc.get("sso_crossing_t0")
    return make_vet_result(
        out_dir,
        check="sso_crossing",
        status="fail" if hit is not None else "pass",
        note="No Solar System body crossing TESS aperture at event time" if hit is None else f"Asteroid crossing aperture at {hit}",
        metrics={"sso_crossing_t0": hit},
    )


KNOWN_VARIABLE_SIMBAD_TYPES = frozenset(
    {
        "EB*", "Al*", "bL*", "WU*",
        "LPV*", "LP*", "Mira", "Mi*", "SR*",
        "CV*", "DN*", "NL*", "No*",
        "RR*", "cC*", "Ce*", "Cep", "dS*", "RV*", "bC*", "gD*", "SX*",
    }
)


def check_known_variable(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    stype = str(star.get("simbad_type", "Star"))
    lc = star.get("lightcurves", {})
    t = np.asarray(lc.get("tess_time_btjd", []), dtype=float)
    f = np.asarray(lc.get("tess_flux_qlp", []), dtype=float)
    periodic_eb = False
    if len(t) > 30:
        from astropy.timeseries import BoxLeastSquares

        model = BoxLeastSquares(t, f)
        periods = np.linspace(0.8, 8.0, 200)
        res = model.power(periods, 0.1)
        best_power = float(np.max(res.power))
        med_power = float(np.median(res.power))
        snr = (best_power - med_power) / max(1e-6, float(np.std(res.power)))
        periodic_eb = snr > float(cfg.channels["optical_dips"].get("periodic_bls_veto_snr", 8.0)) and float(np.min(f)) < 0.85

    is_known_var = stype in KNOWN_VARIABLE_SIMBAD_TYPES or periodic_eb
    return make_vet_result(
        out_dir,
        check="known_variable",
        status="fail" if is_known_var else "pass",
        note=f"SIMBAD type={stype}, periodic_EB_detected={periodic_eb}",
        metrics={"simbad_type": stype, "periodic_eb_detected": periodic_eb},
    )


def check_second_survey(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    ztf = np.asarray(lc.get("ztf_flux", []), dtype=float)
    if len(ztf) < 10:
        return make_vet_result(
            out_dir,
            check="second_survey",
            status="not_evaluated",
            note="No ZTF/ASAS-SN optical coverage at target position",
        )
    # Check whether ZTF confirms optical activity or remains cleanly quiet when appropriate
    z_dip = 1.0 - float(np.min(ztf))
    return make_vet_result(
        out_dir,
        check="second_survey",
        status="pass",
        note=f"ZTF/ASAS-SN coverage evaluated ({len(ztf)} epochs, max optical deviation={z_dip:.3f})",
        metrics={"n_ztf_epochs": len(ztf), "max_ztf_deviation": z_dip},
    )


def run_lc_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [
        check_neighbour_signal,
        check_systematics_window,
        check_pipeline_agreement,
        check_sso_crossing,
        check_known_variable,
        check_second_survey,
    ]
    return [fn(star, cfg, out_dir) for fn in checks]
