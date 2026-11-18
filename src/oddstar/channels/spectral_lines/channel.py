"""Channel 7: Unresolved non-stellar emission lines in optical spectra (Zuckerman et al. 2023)."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import nnls

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec
from oddstar.vet.spectral_checks import KNOWN_SKY_LINES_A, KNOWN_STELLAR_LINES_A, LINE_EXCLUSION_LSF_FWHM

# Fewer wavelength points than this is not a spectrum the search can use
MIN_SPECTRUM_POINTS = 20
# Line-spread FWHM (A) at which the photospheric H-alpha templates below are defined; a wider instrument profile
# broadens them in quadrature by the excess
REFERENCE_LSF_FWHM_A = 0.45
COVERAGE_CATEGORIES = ("searched", "single_epoch", "no_coverage", "failed", "not_fetched", "no_spectrum")


def spectra_search_status(star: dict[str, Any]) -> str:
    """Whether the line search had data for this star: "searched" or the reason it did not.

    "searched" needs a wavelength grid of at least `MIN_SPECTRUM_POINTS` points and two exposures of the same
    length, the condition `measure` applies. Otherwise the survey fetch record says why (`single_epoch`,
    `no_coverage`, `failed`, `not_fetched`); a star with no fetch record and no spectra is `no_spectrum`. A star
    that was not searched must never be counted in a null result.
    """
    spec = star.get("spectra") or {}
    n = len(spec.get("wavelength_a") or [])
    if n >= MIN_SPECTRUM_POINTS and len(spec.get("flux_exp1") or []) == n and len(spec.get("flux_exp2") or []) == n:
        return "searched"
    status = str((star.get("spectra_fetch") or {}).get("status") or "")
    if status.startswith("failed"):
        return "failed"
    if status in {"single_epoch", "no_coverage", "not_fetched"}:
        return status
    return "no_spectrum"


def search_coverage_summary(stars: list[dict[str, Any]], cfg: PipelineConfig) -> dict[str, Any]:
    """Counts of searched and not-searched stars, in the layout the other survey channels use."""
    ids: dict[str, list[str]] = {c: [] for c in COVERAGE_CATEGORIES if c != "searched"}
    n_searched = 0
    for s in stars:
        cat = spectra_search_status(s)
        if cat == "searched":
            n_searched += 1
        else:
            ids[cat].append(str(s["source_id"]))
    out: dict[str, Any] = {"n_input": len(stars), "n_searched": n_searched, "n_not_searched": len(stars) - n_searched}
    for cat, lst in ids.items():
        out[f"n_{cat}"] = len(lst)
    out["not_searched_ids"] = {cat: lst for cat, lst in ids.items() if lst}
    return out


class SpectralLinesChannel:
    name = "spectral_lines"

    def __init__(
        self,
        cfg: PipelineConfig,
        root: Path,
        *,
        stellar_line_veto_enabled: bool = True,
    ) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["spectral_lines"]
        self.root = root
        self.stellar_line_veto_enabled = stellar_line_veto_enabled

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec:
        chash = compute_config_hash(cfg.model_dump())[:12]
        spec, _ = build_parent_sample(
            self.name,
            cfg,
            chash,
            self.root,
            extra_stars=kwargs.get("extra_stars"),
            limit=kwargs.get("limit"),
        )
        return spec

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        """Largest emission residual that appears in both exposures, after a template fit, away from known lines.

        `spectra.lsf_fwhm_a` sets the matched-filter width; the H-alpha templates are broadened by it and the
        exclusion zones around known lines are at least `LINE_EXCLUSION_LSF_FWHM` times it. `spectra.sky_line_wavelengths_a`, when
        present, is where the night-sky lines fall on this star's wavelength grid (a survey spectrum moved to the
        stellar rest frame has them displaced); every listed position is excluded. Without it the laboratory
        airglow list is used.
        """
        spec = star.get("spectra") or {}
        wave = np.asarray(spec.get("wavelength_a", []), dtype=float)
        e1 = np.asarray(spec.get("flux_exp1", []), dtype=float)
        e2 = np.asarray(spec.get("flux_exp2", []), dtype=float)
        sky = np.asarray(spec.get("sky_flux", []), dtype=float)

        extra: dict[str, Any] = {"search_status": spectra_search_status(star)}
        if len(wave) < MIN_SPECTRUM_POINTS or len(e1) != len(wave) or len(e2) != len(wave):
            vals = {"max_residual_snr": 0.0, "peak_wavelength_a": 0.0, "repeat_ratio": 0.0, "matched_sky": 0.0, "matched_stellar": 0.0, "stellar_veto_enabled": 1.0 if self.stellar_line_veto_enabled else 0.0}
        else:
            # Build library template matrix (continuum + photospheric H-alpha Voigt/Gaussian profiles of widths 1.4-2.2 A)
            lsf_fwhm = float(spec.get("lsf_fwhm_a", REFERENCE_LSF_FWHM_A))
            extra_var = max(0.0, lsf_fwhm**2 - REFERENCE_LSF_FWHM_A**2) / 2.355**2
            t1 = np.ones_like(wave)
            t2 = -np.exp(-0.5 * (wave - 6562.8) ** 2 / (1.8**2 + extra_var))
            t3 = -np.exp(-0.5 * (wave - 6562.8) ** 2 / (1.4**2 + extra_var))
            t4 = (wave - np.mean(wave)) / 100.0
            A = np.column_stack([t1, t2, t3, t4])

            # Fit non-negative template combination on minimum (robust against cosmic rays in 1 exposure)
            comb = np.minimum(e1, e2)
            # Fit on shifted positive data
            coeffs, _ = nnls(A, comb)
            model = A @ coeffs

            r1 = e1 - model
            r2 = e2 - model
            # Robust noise per exposure
            sig1 = max(1e-4, 1.4826 * float(np.median(np.abs(r1 - np.median(r1)))))
            sig2 = max(1e-4, 1.4826 * float(np.median(np.abs(r2 - np.median(r2)))))

            # Matched filter at LSF width (~0.3 A sigma)
            lsf_sig = max(0.15, lsf_fwhm / 2.355)
            dw = float(wave[1] - wave[0])
            kern_x = np.arange(-4, 5) * dw
            kern = np.exp(-0.5 * (kern_x / lsf_sig) ** 2)
            kern /= math.sqrt(float(np.sum(kern**2)))

            snr1 = np.convolve(r1 / sig1, kern, mode="same")
            snr2 = np.convolve(r2 / sig2, kern, mode="same")
            joint_snr = np.minimum(snr1, snr2)

            # Mask out sky lines & optionally stellar lines
            sky_tol = max(float(self.ch_cfg["sky_line_exclusion_angstrom"]), LINE_EXCLUSION_LSF_FWHM * lsf_fwhm)
            st_tol = max(float(self.ch_cfg["stellar_line_exclusion_angstrom"]), LINE_EXCLUSION_LSF_FWHM * lsf_fwhm)
            mask = np.ones_like(wave, dtype=bool)
            if len(sky) == len(wave):
                sky_sig = max(1e-4, 1.4826 * float(np.median(np.abs(sky - np.median(sky)))))
                for sw in wave[(sky - np.median(sky)) / sky_sig > 6.0]:
                    mask &= np.abs(wave - sw) > sky_tol
            sky_lines = spec.get("sky_line_wavelengths_a")
            if sky_lines is None:
                # Laboratory airglow positions; geocoronal H-alpha coincides with the stellar line handled below
                sky_lines = [sl for sl in KNOWN_SKY_LINES_A if sl != 6562.79]
            for sl in sky_lines:
                mask &= np.abs(wave - float(sl)) > sky_tol
            if self.stellar_line_veto_enabled:
                for st in KNOWN_STELLAR_LINES_A:
                    mask &= np.abs(wave - st) > st_tol

            masked_snr = np.where(mask, joint_snr, 0.0)
            best_idx = int(np.argmax(masked_snr))
            best_snr = float(masked_snr[best_idx])
            best_wave = float(wave[best_idx])

            matched_sky = 0.0 if mask[best_idx] else 1.0
            # How much of the window was open to the search, and where the peak sits on the detector side: the
            # same observed wavelength in many stars points at the sky or the instrument, not at the stars
            extra["searched_window_fraction"] = round(float(np.mean(mask)), 3)
            factors = [float(x) for x in (spec.get("rest_to_observed_factors") or [])]
            if factors:
                extra["peak_observed_wavelengths_a"] = [round(best_wave * x, 2) for x in factors]
            matched_stellar = 1.0 if any(abs(best_wave - st) <= st_tol for st in KNOWN_STELLAR_LINES_A) else 0.0
            repeat_ratio = min(float(snr1[best_idx]), float(snr2[best_idx])) / max(1e-3, max(float(snr1[best_idx]), float(snr2[best_idx])))

            vals = {
                "max_residual_snr": float(best_snr),
                "peak_wavelength_a": float(best_wave),
                "repeat_ratio": float(repeat_ratio),
                "matched_sky": float(matched_sky),
                "matched_stellar": float(matched_stellar),
                "stellar_veto_enabled": 1.0 if self.stellar_line_veto_enabled else 0.0,
            }

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=vals,
            uncertainties={"max_residual_snr": 1.0},
            epochs=["multi-exposure"],
            caveats=["Requires multi-exposure repeatability and night-sky/stellar line checks."],
            extra=extra,
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        thr = float(self.ch_cfg["residual_snr_threshold"])
        flagged = v["max_residual_snr"] >= thr and v["repeat_ratio"] >= 0.4
        raw_stat = float(max(0.0, v["max_residual_snr"]))
        from scipy.stats import norm

        fap = float(max(1e-12, norm.sf(raw_stat)))
        status = str(m.extra.get("search_status", "searched"))
        summary = f"Spectral residual SNR={v['max_residual_snr']:.2f} at {v['peak_wavelength_a']:.1f}A (stellar_veto={int(v['stellar_veto_enabled'])})"
        if status != "searched" and not flagged:
            summary = f"NOT SEARCHED ({status}): no pair of exposures to compare"
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=bool(flagged),
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=summary,
            feature_contributions={
                "max_residual_snr": float(v["max_residual_snr"]),
                "repeat_ratio": float(v["repeat_ratio"]),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "spectral_lines_measurement.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        path.write_text(txt, encoding="utf-8")
        return [
            EvidenceRef(
                kind="spectral_residual_measurement",
                path=str(path),
                description="NNLS template-subtracted LSF matched-filter residual spectrum summary",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        ]
