"""Spectral emission line vetting checks (Section 6)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result

KNOWN_SKY_LINES_A = [6562.79, 6577.20, 6596.60, 6604.10]  # H-alpha geocoronal & OH airglow lines
KNOWN_STELLAR_LINES_A = [6562.80, 6548.05, 6583.45]  # H-alpha, [N II]


def check_cosmic_ray(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    spec = star.get("spectra", {})
    e1 = np.asarray(spec.get("flux_exp1", []), dtype=float)
    e2 = np.asarray(spec.get("flux_exp2", []), dtype=float)
    wave = np.asarray(spec.get("wavelength_a", []), dtype=float)
    if len(e1) == 0 or len(e2) == 0:
        return make_vet_result(out_dir, check="cosmic_ray", status="not_evaluated", note="Missing multi-exposure spectra")
    diff = e1 - e2
    mad = max(1e-4, 1.4826 * float(np.median(np.abs(diff - np.median(diff)))))
    v = cfg.vetting
    cr_spike = float(np.max(np.abs(diff)) / mad) > float(v["spectral_cr_spike_mad_sigma"])
    # Also check peak FWHM vs instrument LSF
    lsf_a = float(spec.get("lsf_fwhm_a", 0.45))
    peak_idx = int(np.argmax(e1))
    above_half = np.sum(e1[max(0, peak_idx - 3) : min(len(e1), peak_idx + 4)] > float(v["spectral_peak_width_flux_level"]))
    width_a = float(above_half) * (float(wave[1] - wave[0]) if len(wave) > 1 else 0.25)
    too_narrow = cr_spike or (
        float(np.max(e1)) > float(v["spectral_cr_min_peak_flux"]) and width_a < float(v["spectral_cr_min_width_lsf_frac"]) * lsf_a
    )
    return make_vet_result(
        out_dir,
        check="cosmic_ray",
        status="fail" if too_narrow else "pass",
        note="No single-exposure sub-LSF cosmic ray spikes" if not too_narrow else f"Single-exposure sub-LSF spike detected (width={width_a:.2f}A < LSF {lsf_a:.2f}A)",
        metrics={"cr_spike_diff": cr_spike, "peak_width_a": width_a, "lsf_fwhm_a": lsf_a},
    )


def check_sky_line(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    spec = star.get("spectra", {})
    wave = np.asarray(spec.get("wavelength_a", []), dtype=float)
    sky = np.asarray(spec.get("sky_flux", []), dtype=float)
    e1 = np.asarray(spec.get("flux_exp1", []), dtype=float)
    e2 = np.asarray(spec.get("flux_exp2", []), dtype=float)
    if len(wave) == 0 or len(sky) == 0:
        return make_vet_result(out_dir, check="sky_line", status="not_evaluated", note="No sky spectrum available")
    sky_med = float(np.median(sky))
    sky_mad = max(1e-4, 1.4826 * float(np.median(np.abs(sky - sky_med))))
    sky_peaks = wave[(sky - sky_med) / sky_mad > float(cfg.vetting["spectral_sky_peak_mad_sigma"])]
    comb = 0.5 * (e1 + e2)
    target_peaks = wave[(comb - float(np.median(comb))) / max(1e-4, float(np.std(comb))) > float(cfg.vetting["spectral_target_peak_sigma"])]
    tol = float(cfg.channels["spectral_lines"].get("sky_line_exclusion_angstrom", 1.5))
    matched_sky = any(abs(float(tp) - float(sp)) <= tol for tp in target_peaks for sp in sky_peaks)
    matched_oh = any(abs(float(tp) - sl) <= tol for tp in target_peaks for sl in KNOWN_SKY_LINES_A if sl != 6562.79)
    bad = matched_sky or matched_oh
    return make_vet_result(
        out_dir,
        check="sky_line",
        status="fail" if bad else "pass",
        note="Residual emission clear of sky/OH airglow lines" if not bad else "Emission feature coincides with sky-spectrum airglow line",
        metrics={"matched_sky_spectrum": matched_sky, "matched_airglow_catalog": matched_oh},
    )


def check_stellar_emission(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    spec = star.get("spectra", {})
    wave = np.asarray(spec.get("wavelength_a", []), dtype=float)
    e1 = np.asarray(spec.get("flux_exp1", []), dtype=float)
    e2 = np.asarray(spec.get("flux_exp2", []), dtype=float)
    if len(wave) == 0:
        return make_vet_result(out_dir, check="stellar_emission", status="not_evaluated", note="No optical spectrum")
    comb = 0.5 * (e1 + e2)
    tol = float(cfg.channels["spectral_lines"].get("stellar_line_exclusion_angstrom", 1.5))
    # Check for chromospheric emission near H-alpha (6562.8 A) or [N II]
    peak_idx = int(np.argmax(comb))
    peak_w = float(wave[peak_idx])
    peak_amp = float(comb[peak_idx] - np.median(comb))
    is_stellar_line = peak_amp > float(cfg.vetting["spectral_stellar_line_min_amp"]) and any(abs(peak_w - sl) <= tol for sl in KNOWN_STELLAR_LINES_A)
    return make_vet_result(
        out_dir,
        check="stellar_emission",
        status="fail" if is_stellar_line else "pass",
        note=f"Peak at {peak_w:.2f} A ({'matches chromospheric stellar line' if is_stellar_line else 'clear of stellar lines'})",
        metrics={"peak_wavelength_a": peak_w, "is_stellar_line": is_stellar_line},
    )


def check_detector(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    det = bool(star.get("spectra", {}).get("detector_artifact", False))
    return make_vet_result(
        out_dir,
        check="detector",
        status="fail" if det else "pass",
        note="No cross-talk, ghost, or bad-column artifact at detector coordinates" if not det else "Detector bad-pixel/ghost artifact matched",
        metrics={"detector_artifact": det},
    )


def check_repeat(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    spec = star.get("spectra", {})
    e1 = np.asarray(spec.get("flux_exp1", []), dtype=float)
    e2 = np.asarray(spec.get("flux_exp2", []), dtype=float)
    if len(e1) == 0 or len(e2) == 0:
        return make_vet_result(out_dir, check="repeat", status="not_evaluated", note="Single exposure only")
    idx1 = int(np.argmax(e1))
    idx2 = int(np.argmax(e2))
    val1 = float(e1[idx1] - np.median(e1))
    val2 = float(e2[idx1] - np.median(e2))
    quiet_max = float(cfg.vetting["spectral_repeat_quiet_max_excess"])
    repeated = (val1 < quiet_max and val2 < quiet_max) or (
        abs(idx1 - idx2) <= int(cfg.vetting["spectral_repeat_max_pixel_shift"])
        and val2 > float(cfg.vetting["spectral_repeat_min_excess"])
    )
    return make_vet_result(
        out_dir,
        check="repeat",
        status="pass" if repeated else "fail",
        note="Feature repeats consistently across separate exposures" if repeated else "Feature present in only one exposure (fails repeat check)",
        metrics={"exp1_excess": val1, "exp2_at_exp1_peak": val2, "repeated": repeated},
    )


def run_spectral_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [check_cosmic_ray, check_sky_line, check_stellar_emission, check_detector, check_repeat]
    return [fn(star, cfg, out_dir) for fn in checks]
