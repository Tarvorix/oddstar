"""Spectral emission line vetting checks (Section 6).

Input rule (R9, lessons L007/L008): a check returns `pass` or `fail` only from inputs that were actually obtained.

- The spectral_lines channel has NO live ingestion: no code path fetches per-exposure optical spectra, a sky
  spectrum, the instrument line-spread function or a detector artifact map for a live star (`live_enriched` key
  sets `spectra: {}`). A live star counts a spectral input only when `vetting_archive.query_status["spectra"] ==
  "ok"`, which nothing writes today, so every check here returns `not_evaluated` for live stars and names the
  missing data product.
- Synthetic fixtures and hand-built dicts count an input only when its key is explicitly present, not None and
  (for arrays) non-empty with matching lengths. There is no default LSF width, pixel size or artifact flag.

Every evaluated result stores the input values it used (array sizes, wavelength range, LSF, flags) in its metrics.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result

KNOWN_SKY_LINES_A = [6562.79, 6577.20, 6596.60, 6604.10]  # H-alpha geocoronal & OH airglow lines
KNOWN_STELLAR_LINES_A = [6562.80, 6548.05, 6583.45]  # H-alpha, [N II]

SPECTRA_QUERY = "spectra"
_PRODUCTS = {
    "wavelength_a": "wavelength grid of the optical spectrum",
    "flux_exp1": "first-exposure optical spectrum",
    "flux_exp2": "second independent exposure",
    "sky_flux": "sky spectrum of the same exposure",
    "lsf_fwhm_a": "instrument line-spread function FWHM",
    "detector_artifact": "detector bad-pixel / ghost / cross-talk map lookup",
}
NO_LIVE_INGESTION = "the spectral_lines channel has no live ingestion of survey spectra"


def _is_live(star: dict[str, Any]) -> bool:
    return "live_enriched" in star


def _spec(star: dict[str, Any]) -> dict[str, Any]:
    if _is_live(star):
        qs = ((star.get("vetting_archive") or {}).get("query_status") or {}).get(SPECTRA_QUERY)
        if qs != "ok":
            return {}
    return star.get("spectra") or {}


def _array(spec: dict[str, Any], key: str) -> np.ndarray | None:
    val = spec.get(key)
    if val is None:
        return None
    arr = np.asarray(val, dtype=float)
    return arr if arr.ndim == 1 and len(arr) > 0 else None


def _arrays(star: dict[str, Any], keys: list[str]) -> tuple[dict[str, np.ndarray], list[str]]:
    """Requested arrays, or the list of keys that are missing / empty / of mismatched length."""
    spec = _spec(star)
    got = {k: _array(spec, k) for k in keys}
    missing = [k for k, v in got.items() if v is None]
    if not missing:
        n = {len(v) for v in got.values() if v is not None}
        if len(n) > 1:
            missing = [f"{k} (length mismatch: " + ", ".join(f"{kk}={len(vv)}" for kk, vv in got.items() if vv is not None) + ")"
                       for k in keys[:1]]
    return {k: v for k, v in got.items() if v is not None}, missing


def _scalar(star: dict[str, Any], key: str) -> Any:
    return _spec(star).get(key)


def _not_evaluated(out_dir: Path, check: str, star: dict[str, Any], missing: list[str], metrics: dict[str, Any] | None = None) -> VettingResult:
    parts = [f"{_PRODUCTS.get(k.split(' ')[0], k)} [spectra.{k}]" for k in missing]
    note = "Missing input: " + "; ".join(parts)
    if _is_live(star):
        qs = ((star.get("vetting_archive") or {}).get("query_status") or {}).get(SPECTRA_QUERY) or "not run"
        note += f" (query_status {SPECTRA_QUERY}={qs}; {NO_LIVE_INGESTION})"
    m = dict(metrics or {}, missing_inputs=missing, live_star=_is_live(star))
    return make_vet_result(out_dir, check=check, status="not_evaluated", note=note, metrics=m)


def _grid(wave: np.ndarray) -> dict[str, float]:
    return {"n_pixels": int(len(wave)), "wavelength_min_a": float(wave[0]), "wavelength_max_a": float(wave[-1]),
            "pixel_a": float(np.median(np.diff(wave))) if len(wave) > 1 else float("nan")}


def check_cosmic_ray(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    arrs, missing = _arrays(star, ["wavelength_a", "flux_exp1", "flux_exp2"])
    if not missing and len(arrs["wavelength_a"]) < 2:
        missing = ["wavelength_a (fewer than 2 pixels)"]
    if missing:
        return _not_evaluated(out_dir, "cosmic_ray", star, missing)
    wave, e1, e2 = arrs["wavelength_a"], arrs["flux_exp1"], arrs["flux_exp2"]
    v = cfg.vetting
    diff = e1 - e2
    mad = max(1e-4, 1.4826 * float(np.median(np.abs(diff - np.median(diff)))))
    spike_sigma = float(np.max(np.abs(diff)) / mad)
    cr_spike = spike_sigma > float(v["spectral_cr_spike_mad_sigma"])
    peak_idx = int(np.argmax(e1))
    pixel_a = float(wave[1] - wave[0])
    above_half = np.sum(e1[max(0, peak_idx - 3) : min(len(e1), peak_idx + 4)] > float(v["spectral_peak_width_flux_level"]))
    width_a = float(above_half) * pixel_a
    peak_flux = float(np.max(e1))
    lsf_raw = _scalar(star, "lsf_fwhm_a")
    metrics: dict[str, Any] = {
        "input_values": dict(_grid(wave), lsf_fwhm_a=lsf_raw),
        "cr_spike_diff": cr_spike, "max_abs_diff_over_mad": spike_sigma, "diff_mad": mad, "peak_flux_exp1": peak_flux,
        "peak_wavelength_a": float(wave[peak_idx]), "peak_width_a": width_a, "lsf_fwhm_a": lsf_raw,
        "thresholds": {k: float(v[k]) for k in ("spectral_cr_spike_mad_sigma", "spectral_peak_width_flux_level", "spectral_cr_min_peak_flux",
                                                 "spectral_cr_min_width_lsf_frac")},
    }
    if cr_spike:
        return make_vet_result(out_dir, check="cosmic_ray", status="fail", metrics=metrics,
                               note=f"Single-exposure spike: |exp1 - exp2| reaches {spike_sigma:.1f} x MAD (> {float(v['spectral_cr_spike_mad_sigma']):.1f})")
    if lsf_raw is None:
        return _not_evaluated(out_dir, "cosmic_ray", star, ["lsf_fwhm_a"], metrics)
    lsf_a = float(lsf_raw)
    too_narrow = peak_flux > float(v["spectral_cr_min_peak_flux"]) and width_a < float(v["spectral_cr_min_width_lsf_frac"]) * lsf_a
    return make_vet_result(
        out_dir,
        check="cosmic_ray",
        status="fail" if too_narrow else "pass",
        note=(f"Single-exposure sub-LSF spike detected (width={width_a:.2f}A < {float(v['spectral_cr_min_width_lsf_frac']):.2f} x LSF {lsf_a:.2f}A)"
              if too_narrow else
              f"No single-exposure sub-LSF cosmic ray spike (max |exp1-exp2| = {spike_sigma:.1f} x MAD, peak width {width_a:.2f}A, LSF {lsf_a:.2f}A)"),
        metrics=metrics,
    )


def check_sky_line(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    arrs, missing = _arrays(star, ["wavelength_a", "flux_exp1", "flux_exp2"])
    if missing:
        return _not_evaluated(out_dir, "sky_line", star, missing)
    wave, e1, e2 = arrs["wavelength_a"], arrs["flux_exp1"], arrs["flux_exp2"]
    comb = 0.5 * (e1 + e2)
    target_sigma = float(cfg.vetting["spectral_target_peak_sigma"])
    target_peaks = wave[(comb - float(np.median(comb))) / max(1e-4, float(np.std(comb))) > target_sigma]
    tol = float(cfg.channels["spectral_lines"]["sky_line_exclusion_angstrom"])
    matched_oh = any(abs(float(tp) - sl) <= tol for tp in target_peaks for sl in KNOWN_SKY_LINES_A if sl != 6562.79)
    metrics: dict[str, Any] = {
        "input_values": _grid(wave), "n_target_peaks": int(len(target_peaks)),
        "target_peak_wavelengths_a": [round(float(x), 2) for x in target_peaks[:20]], "matched_airglow_catalog": matched_oh,
        "exclusion_angstrom": tol, "target_peak_sigma": target_sigma,
    }
    sky_arrs, sky_missing = _arrays(star, ["wavelength_a", "sky_flux"])
    if sky_missing:
        if matched_oh:
            return make_vet_result(out_dir, check="sky_line", status="fail", metrics=metrics,
                                   note="Emission feature coincides with a catalogued OH airglow line (no sky spectrum available)")
        return _not_evaluated(out_dir, "sky_line", star, ["sky_flux"], metrics)
    sky = sky_arrs["sky_flux"]
    sky_med = float(np.median(sky))
    sky_mad = max(1e-4, 1.4826 * float(np.median(np.abs(sky - sky_med))))
    sky_sigma = float(cfg.vetting["spectral_sky_peak_mad_sigma"])
    sky_peaks = wave[(sky - sky_med) / sky_mad > sky_sigma]
    matched_sky = any(abs(float(tp) - float(sp)) <= tol for tp in target_peaks for sp in sky_peaks)
    bad = matched_sky or matched_oh
    metrics.update({"matched_sky_spectrum": matched_sky, "n_sky_peaks": int(len(sky_peaks)),
                    "sky_peak_wavelengths_a": [round(float(x), 2) for x in sky_peaks[:20]], "sky_peak_mad_sigma": sky_sigma})
    return make_vet_result(
        out_dir,
        check="sky_line",
        status="fail" if bad else "pass",
        note="Emission feature coincides with sky-spectrum airglow line" if bad else "Residual emission clear of sky/OH airglow lines",
        metrics=metrics,
    )


def check_stellar_emission(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    arrs, missing = _arrays(star, ["wavelength_a", "flux_exp1", "flux_exp2"])
    if missing:
        return _not_evaluated(out_dir, "stellar_emission", star, missing)
    wave, e1, e2 = arrs["wavelength_a"], arrs["flux_exp1"], arrs["flux_exp2"]
    comb = 0.5 * (e1 + e2)
    tol = float(cfg.channels["spectral_lines"]["stellar_line_exclusion_angstrom"])
    min_amp = float(cfg.vetting["spectral_stellar_line_min_amp"])
    # Chromospheric emission near H-alpha (6562.8 A) or [N II]
    peak_idx = int(np.argmax(comb))
    peak_w = float(wave[peak_idx])
    peak_amp = float(comb[peak_idx] - np.median(comb))
    is_stellar_line = peak_amp > min_amp and any(abs(peak_w - sl) <= tol for sl in KNOWN_STELLAR_LINES_A)
    return make_vet_result(
        out_dir,
        check="stellar_emission",
        status="fail" if is_stellar_line else "pass",
        note=f"Peak at {peak_w:.2f} A, amplitude {peak_amp:.3f} ({'matches chromospheric stellar line' if is_stellar_line else 'clear of stellar lines'})",
        metrics={"input_values": _grid(wave), "peak_wavelength_a": peak_w, "peak_amplitude": peak_amp, "is_stellar_line": is_stellar_line,
                 "min_amplitude": min_amp, "exclusion_angstrom": tol},
    )


def check_detector(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    raw = _scalar(star, "detector_artifact")
    if raw is None:
        return _not_evaluated(out_dir, "detector", star, ["detector_artifact"], {"input_values": {"detector_artifact": None}})
    det = bool(raw)
    return make_vet_result(
        out_dir,
        check="detector",
        status="fail" if det else "pass",
        note="Detector bad-pixel/ghost artifact matched" if det else "No cross-talk, ghost, or bad-column artifact at detector coordinates",
        metrics={"input_values": {"detector_artifact": det}, "detector_artifact": det},
    )


def check_repeat(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    arrs, missing = _arrays(star, ["flux_exp1", "flux_exp2"])
    if missing:
        return _not_evaluated(out_dir, "repeat", star, missing)
    e1, e2 = arrs["flux_exp1"], arrs["flux_exp2"]
    idx1 = int(np.argmax(e1))
    idx2 = int(np.argmax(e2))
    val1 = float(e1[idx1] - np.median(e1))
    val2 = float(e2[idx1] - np.median(e2))
    v = cfg.vetting
    quiet_max = float(v["spectral_repeat_quiet_max_excess"])
    max_shift = int(v["spectral_repeat_max_pixel_shift"])
    min_excess = float(v["spectral_repeat_min_excess"])
    repeated = (val1 < quiet_max and val2 < quiet_max) or (abs(idx1 - idx2) <= max_shift and val2 > min_excess)
    return make_vet_result(
        out_dir,
        check="repeat",
        status="pass" if repeated else "fail",
        note="Feature repeats consistently across separate exposures" if repeated else "Feature present in only one exposure (fails repeat check)",
        metrics={"input_values": {"n_pixels": int(len(e1))}, "exp1_excess": val1, "exp2_at_exp1_peak": val2, "exp1_peak_index": idx1,
                 "exp2_peak_index": idx2, "repeated": repeated,
                 "thresholds": {"quiet_max_excess": quiet_max, "max_pixel_shift": max_shift, "min_excess": min_excess}},
    )


def run_spectral_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [check_cosmic_ray, check_sky_line, check_stellar_emission, check_detector, check_repeat]
    return [fn(star, cfg, out_dir) for fn in checks]
