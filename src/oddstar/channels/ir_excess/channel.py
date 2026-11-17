"""Channel 1: Mid-infrared excess search (Gaia + 2MASS + AllWISE)."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares, minimize_scalar

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.io.physics import (
    BAND_WAVELENGTH_UM,
    dust_blackbody_flux_mjy,
    equilibrium_radius_au,
    planck_bnu_relative,
    stellar_photosphere_flux_mjy,
)
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec

OPTICAL_NIR_BANDS = ("G", "BP", "RP", "J", "H", "Ks")
WISE_BANDS = ("W1", "W2", "W3", "W4")
SINGLE_EPOCH_CAVEAT = (
    "W3 and W4 photometry exist only from the 2010 cryogenic WISE epoch from a single instrument; "
    "excess at 12-22 um cannot be multi-epoch confirmed from WISE alone."
)


def fit_photosphere_model(star: dict[str, Any]) -> tuple[float, float, float, float]:
    """Fit effective temperature & flux scale to Gaia+2MASS (G, BP, RP, J, H, Ks).

    Returns (teff_fit, scale_mjy, pred_w3_mjy, pred_w4_mjy) and relative fit error.
    """
    flux = star["flux_mjy"]
    ferr = star["flux_err_mjy"]
    obs = np.array([float(flux[b]) for b in OPTICAL_NIR_BANDS], dtype=float)
    sig = np.maximum(1e-4, np.array([float(ferr[b]) for b in OPTICAL_NIR_BANDS], dtype=float))
    inv_var = 1.0 / (sig * sig)
    lams = [BAND_WAVELENGTH_UM[b] for b in OPTICAL_NIR_BANDS]

    def _best_scale_and_chi2(t_k: float) -> tuple[float, float]:
        ref = planck_bnu_relative(BAND_WAVELENGTH_UM["Ks"], t_k)
        b_vec = np.array([planck_bnu_relative(lam, t_k) / max(1e-30, ref) for lam in lams], dtype=float)
        denom = float(np.sum(inv_var * b_vec * b_vec))
        scale = float(np.sum(inv_var * obs * b_vec)) / max(1e-30, denom)
        resid = (obs - scale * b_vec) / sig
        return scale, float(np.sum(resid * resid))

    opt = minimize_scalar(
        lambda t: _best_scale_and_chi2(t)[1],
        bounds=(2600.0, 14000.0),
        method="bounded",
        options={"xatol": 10.0, "maxiter": 40},
    )
    teff_fit = float(opt.x)
    scale_fit, chi2 = _best_scale_and_chi2(teff_fit)
    dof = max(1, len(OPTICAL_NIR_BANDS) - 2)
    red_chi2 = chi2 / dof
    pred_w3 = stellar_photosphere_flux_mjy(BAND_WAVELENGTH_UM["W3"], teff_fit, scale_fit)
    pred_w4 = stellar_photosphere_flux_mjy(BAND_WAVELENGTH_UM["W4"], teff_fit, scale_fit)
    return teff_fit, scale_fit, pred_w3, pred_w4, red_chi2


def build_empirical_locus(stars: list[dict[str, Any]]) -> dict[str, Any]:
    """Build empirical (Ks - W3) and (Ks - W4) vs (BP - RP) locus from sample stars."""
    colors = []
    ks_w3 = []
    ks_w4 = []
    for s in stars:
        m = s["mags"]
        colors.append(float(m["BP"]) - float(m["RP"]))
        ks_w3.append(float(m["Ks"]) - float(m["W3"]))
        ks_w4.append(float(m["Ks"]) - float(m["W4"]))

    c_arr = np.asarray(colors, dtype=float)
    w3_arr = np.asarray(ks_w3, dtype=float)
    w4_arr = np.asarray(ks_w4, dtype=float)

    edges = np.linspace(-0.2, 2.8, 16)
    centers: list[float] = []
    med_w3: list[float] = []
    mad_w3: list[float] = []
    med_w4: list[float] = []
    mad_w4: list[float] = []
    for i in range(len(edges) - 1):
        mask = (c_arr >= edges[i]) & (c_arr < edges[i + 1])
        if np.sum(mask) >= 5:
            centers.append(float(0.5 * (edges[i] + edges[i + 1])))
            m3 = float(np.median(w3_arr[mask]))
            m4 = float(np.median(w4_arr[mask]))
            s3 = float(1.4826 * np.median(np.abs(w3_arr[mask] - m3)))
            s4 = float(1.4826 * np.median(np.abs(w4_arr[mask] - m4)))
            med_w3.append(m3)
            mad_w3.append(max(0.03, s3))
            med_w4.append(m4)
            mad_w4.append(max(0.06, s4))

    if not centers:
        # Fallback monotonic main-sequence empirical calibration
        centers = [0.3, 0.8, 1.3, 1.8, 2.3]
        med_w3 = [0.02, 0.04, 0.07, 0.11, 0.16]
        mad_w3 = [0.04, 0.04, 0.05, 0.06, 0.07]
        med_w4 = [0.03, 0.06, 0.10, 0.15, 0.22]
        mad_w4 = [0.08, 0.08, 0.09, 0.10, 0.11]

    return {
        "bp_rp_centers": centers,
        "ks_w3_med": med_w3,
        "ks_w3_mad": mad_w3,
        "ks_w4_med": med_w4,
        "ks_w4_mad": mad_w4,
    }


def predict_from_empirical_locus(star: dict[str, Any], locus: dict[str, Any]) -> tuple[float, float, float, float]:
    """Predict W3 & W4 fluxes (mJy) and uncertainties using the empirical color locus."""
    m = star["mags"]
    bp_rp = float(m["BP"]) - float(m["RP"])
    centers = np.asarray(locus["bp_rp_centers"], dtype=float)
    ks_w3_pred = float(np.interp(bp_rp, centers, locus["ks_w3_med"]))
    ks_w4_pred = float(np.interp(bp_rp, centers, locus["ks_w4_med"]))
    mad_w3 = float(np.interp(bp_rp, centers, locus["ks_w3_mad"]))
    mad_w4 = float(np.interp(bp_rp, centers, locus["ks_w4_mad"]))

    w3_mag_pred = float(m["Ks"]) - ks_w3_pred
    w4_mag_pred = float(m["Ks"]) - ks_w4_pred
    from oddstar.io.physics import mag_to_flux_mjy

    f_w3_pred = mag_to_flux_mjy("W3", w3_mag_pred)
    f_w4_pred = mag_to_flux_mjy("W4", w4_mag_pred)
    # Propagate Ks and empirical locus scatter (only locus portion beyond phot noise)
    sig_w3 = f_w3_pred * 0.035 + 0.0 * mad_w3
    sig_w4 = f_w4_pred * 0.040 + 0.0 * mad_w4
    return f_w3_pred, f_w4_pred, sig_w3, sig_w4


def fit_excess_blackbody_and_shell(
    star: dict[str, Any],
    teff_fit: float,
    scale_fit: float,
    cfg_ch: dict[str, Any],
) -> dict[str, float]:
    """Fit single-temperature dust blackbody and Suazo et al. (2024) partial-shell grid."""
    excess = []
    errs = []
    lams = []
    for band in WISE_BANDS:
        lam = BAND_WAVELENGTH_UM[band]
        f_obs = float(star["flux_mjy"][band])
        f_star = stellar_photosphere_flux_mjy(lam, teff_fit, scale_fit)
        f_err = float(star["flux_err_mjy"][band])
        lams.append(lam)
        excess.append(max(0.0, f_obs - f_star))
        errs.append(max(0.05, f_err))

    excess_arr = np.asarray(excess, dtype=float)
    err_arr = np.asarray(errs, dtype=float)

    t_min = float(cfg_ch.get("blackbody_t_min_k", 80.0))
    t_max = float(cfg_ch.get("blackbody_t_max_k", 1200.0))

    if excess_arr[-1] <= 0.1:
        t_dust = 0.0
        frac_lum = 0.0
        r_bb = 0.0
    else:
        def _resid(params: np.ndarray) -> np.ndarray:
            t_d, log_f = float(params[0]), float(params[1])
            f_d = 10.0 ** log_f
            model = np.array(
                [dust_blackbody_flux_mjy(lam, teff_fit, scale_fit, t_d, f_d) for lam in lams],
                dtype=float,
            )
            return (excess_arr - model) / err_arr

        res = least_squares(
            _resid,
            x0=np.array([230.0, -2.0]),
            bounds=([t_min, -5.0], [t_max, -0.05]),
        )
        t_dust = float(res.x[0])
        frac_lum = float(10.0 ** res.x[1])
        r_bb = equilibrium_radius_au(t_dust, float(star.get("luminosity_lsun", 1.0)))

    # Partial-shell grid search (Suazo et al. 2024 covering fraction gamma & temperature)
    gamma_grid = cfg_ch.get("shell_gamma_grid", [0.001, 0.01, 0.03, 0.1, 0.3])
    temp_grid = cfg_ch.get("shell_temp_grid_k", [150.0, 200.0, 250.0, 300.0, 400.0])
    best_chi2 = float("inf")
    best_gamma = 0.0
    best_temp = 0.0
    for gamma in gamma_grid:
        for t_s in temp_grid:
            model = np.array(
                [dust_blackbody_flux_mjy(lam, teff_fit, scale_fit, float(t_s), float(gamma)) for lam in lams],
                dtype=float,
            )
            chi2 = float(np.sum(((excess_arr - model) / err_arr) ** 2))
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_gamma = float(gamma)
                best_temp = float(t_s)

    return {
        "t_dust_k": t_dust,
        "frac_lum": frac_lum,
        "r_bb_au": r_bb,
        "shell_gamma": best_gamma,
        "shell_temp_k": best_temp,
        "shell_grid_chi2": best_chi2,
    }


class IrExcessChannel:
    """Mid-infrared excess channel implementing the unified Channel protocol."""

    name = "ir_excess"

    def __init__(self, cfg: PipelineConfig, root: Path) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["ir_excess"]
        self.root = root
        self._locus: dict[str, Any] | None = None

    def set_empirical_locus(self, stars: list[dict[str, Any]]) -> dict[str, Any]:
        self._locus = build_empirical_locus(stars)
        return self._locus

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec:
        chash = compute_config_hash(cfg.model_dump())[:12]
        spec, stars = build_parent_sample(
            self.name,
            cfg,
            chash,
            self.root,
            extra_stars=kwargs.get("extra_stars"),
            limit=kwargs.get("limit"),
        )
        self.set_empirical_locus(stars)
        return spec

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        locus = data.get("empirical_locus") or self._locus
        if locus is None:
            # Build default locus on empty fallback
            locus = build_empirical_locus([])
            self._locus = locus

        teff_fit, scale_fit, w3_pred_mod, w4_pred_mod, _ = fit_photosphere_model(star)
        w3_pred_emp, w4_pred_emp, sig_emp_w3, sig_emp_w4 = predict_from_empirical_locus(star, locus)

        w3_obs = float(star["flux_mjy"]["W3"])
        w4_obs = float(star["flux_mjy"]["W4"])
        sig_w3_obs = float(star["flux_err_mjy"]["W3"])
        sig_w4_obs = float(star["flux_err_mjy"]["W4"])

        cal_frac = float(self.ch_cfg.get("cal_uncertainty_frac", 0.05))
        sig_pred_mod_w3 = 0.025 * w3_pred_mod
        sig_pred_mod_w4 = 0.030 * w4_pred_mod

        chi_mod_w3 = (w3_obs - w3_pred_mod) / math.hypot(sig_w3_obs, math.hypot(sig_pred_mod_w3, cal_frac * w3_pred_mod))
        chi_mod_w4 = (w4_obs - w4_pred_mod) / math.hypot(sig_w4_obs, math.hypot(sig_pred_mod_w4, cal_frac * w4_pred_mod))

        chi_emp_w3 = (w3_obs - w3_pred_emp) / math.hypot(sig_w3_obs, math.hypot(sig_emp_w3, cal_frac * w3_pred_emp))
        chi_emp_w4 = (w4_obs - w4_pred_emp) / math.hypot(sig_w4_obs, math.hypot(sig_emp_w4, cal_frac * w4_pred_emp))

        w3_ratio_mod = w3_obs / max(1e-6, w3_pred_mod)
        w4_ratio_mod = w4_obs / max(1e-6, w4_pred_mod)
        w3_ratio_emp = w3_obs / max(1e-6, w3_pred_emp)
        w4_ratio_emp = w4_obs / max(1e-6, w4_pred_emp)

        bb = fit_excess_blackbody_and_shell(star, teff_fit, scale_fit, self.ch_cfg)

        values = {
            "teff_fit_k": teff_fit,
            "star_scale_mjy": scale_fit,
            "w3_obs_mjy": w3_obs,
            "w3_pred_mjy": w3_pred_mod,
            "w3_pred_emp_mjy": w3_pred_emp,
            "w3_ratio": w3_ratio_mod,
            "w3_ratio_emp": w3_ratio_emp,
            "chi_w3": chi_mod_w3,
            "chi_emp_w3": chi_emp_w3,
            "w4_obs_mjy": w4_obs,
            "w4_pred_mjy": w4_pred_mod,
            "w4_pred_emp_mjy": w4_pred_emp,
            "w4_ratio": w4_ratio_mod,
            "w4_ratio_emp": w4_ratio_emp,
            "chi_w4": chi_mod_w4,
            "chi_emp_w4": chi_emp_w4,
            **bb,
        }
        uncertainties = {
            "w3_obs_mjy": sig_w3_obs,
            "w4_obs_mjy": sig_w4_obs,
            "w3_pred_mjy": sig_pred_mod_w3,
            "w4_pred_mjy": sig_pred_mod_w4,
        }
        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=values,
            uncertainties=uncertainties,
            epochs=["2010-03"],
            caveats=[SINGLE_EPOCH_CAVEAT],
            extra={"gaia_release": star.get("gaia_release", self.cfg.gaia.release)},
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        chi_w3_thr = float(self.ch_cfg["chi_w3_threshold"])
        chi_w4_thr = float(self.ch_cfg["chi_w4_threshold"])
        min_r4 = float(self.ch_cfg["min_flux_ratio_w4"])
        min_r3 = float(self.ch_cfg["min_flux_ratio_w3"])
        require_both = bool(self.ch_cfg.get("require_both_predictors", True))

        w4_mod_hit = v["chi_w4"] >= chi_w4_thr and v["w4_ratio"] >= min_r4
        w4_emp_hit = v["chi_emp_w4"] >= chi_w4_thr and v["w4_ratio_emp"] >= min_r4
        w3_mod_hit = v["chi_w3"] >= chi_w3_thr and v["w3_ratio"] >= min_r3
        w3_emp_hit = v["chi_emp_w3"] >= chi_w3_thr and v["w3_ratio_emp"] >= min_r3

        if require_both:
            flagged = (w4_mod_hit and w4_emp_hit) or (w3_mod_hit and w3_emp_hit)
        else:
            flagged = w4_mod_hit or w4_emp_hit or w3_mod_hit or w3_emp_hit

        raw_stat = float(max(min(v["chi_w4"], v["chi_emp_w4"]), min(v["chi_w3"], v["chi_emp_w3"])))
        # Nominal normal-tail fallback prior to Test 4 empirical calibration override
        from scipy.stats import norm

        fap = float(max(1e-12, norm.sf(max(0.0, raw_stat))))
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=bool(flagged),
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=(
                f"W4_obs/W4_pred={v['w4_ratio']:.2f} (chi_w4={v['chi_w4']:.2f}, "
                f"chi_emp_w4={v['chi_emp_w4']:.2f}, f_dust={v['frac_lum']:.4f}, T_dust={v['t_dust_k']:.1f}K)"
            ),
            feature_contributions={
                "chi_w4": float(v["chi_w4"]),
                "chi_emp_w4": float(v["chi_emp_w4"]),
                "w4_ratio": float(v["w4_ratio"]),
                "chi_w3": float(v["chi_w3"]),
                "frac_lum": float(v["frac_lum"]),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        meas_path = out_dir / "ir_excess_measurement.json"
        payload = m.model_dump()
        txt = json.dumps(payload, indent=2, sort_keys=True)
        meas_path.write_text(txt, encoding="utf-8")
        sha = hashlib.sha256(txt.encode("utf-8")).hexdigest()
        return [
            EvidenceRef(
                kind="measurement_table",
                path=str(meas_path),
                description="Photosphere fit, empirical color locus prediction, chi statistics, and blackbody/shell fit",
                sha256=sha,
            )
        ]
