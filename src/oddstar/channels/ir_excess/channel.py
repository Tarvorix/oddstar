"""Channel 1: Mid-infrared excess search (Gaia + 2MASS + AllWISE).

Two photosphere predictors are evaluated for every star:

* PRIMARY (used for all excess quantification): an empirical, dereddened (Ks - Wx) versus (BP - RP)_0 colour locus
  for W1..W4. It is built from the sample itself with iterative outlier clipping (`sample_locus`), or loaded from
  the most recent locus persisted for the same Gaia release and sample kind (`persisted_locus`), or, only when
  neither exists, taken from the packaged Pecaut & Mamajek dwarf sequence (`literature_fallback`, or
  `literature_fallback_unverified` if those numbers were not retrieved from the source).
* SECONDARY (diagnostic, and the second leg of the two-predictor detection rule): a Planck fit to dereddened
  G, BP, RP, J, H, Ks with the empirical WISE colour correction of `physics.stellar_photosphere_flux_mjy`.
  The uncorrected Planck predictions are also reported (`*_pred_planck_mjy`).

Detection still requires both predictors to agree (`require_both_predictors`), exactly as before.
"""

from __future__ import annotations

import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar, nnls

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.io.physics import (
    BAND_WAVELENGTH_UM,
    blackbody_colour_temperature_k,
    colour_excess,
    dust_blackbody_flux_mjy,
    dwarf_ks_minus_wise_sequence,
    equilibrium_radius_au,
    extinction_mag,
    extinction_two_case_summary,
    literature_photosphere_label,
    load_dwarf_colour_table,
    mag_to_flux_mjy,
    planck_bnu_relative,
    planck_photosphere_flux_mjy,
    star_applied_ebv,
    stellar_photosphere_flux_mjy,
)
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec

OPTICAL_NIR_BANDS = ("G", "BP", "RP", "J", "H", "Ks")
WISE_BANDS = ("W1", "W2", "W3", "W4")
LOCUS_FORMAT_VERSION = 2
LOCUS_DIRNAME = "_photosphere_loci"
LN10_OVER_2P5 = math.log(10.0) / 2.5
SECONDARY_PHOTOSPHERE = "planck_fit_dereddened_plus_dwarf_wise_colour_correction"
SINGLE_EPOCH_CAVEAT = (
    "W3 and W4 photometry exist only from the 2010 cryogenic WISE epoch from a single instrument; "
    "excess at 12-22 um cannot be multi-epoch confirmed from WISE alone."
)


# ------------------------------------------------------------------------------------------------ helpers


def star_sample_kind(star: dict[str, Any]) -> str:
    """`live` for records built from archive rows (they always carry `live_enriched`), otherwise `synthetic`."""
    return "live" if "live_enriched" in star else "synthetic"


def _mag_err(star: dict[str, Any], band: str) -> float:
    me = (star.get("mag_errs") or {}).get(band)
    if me is not None:
        return float(me)
    f = float(star["flux_mjy"][band])
    return 1.0857 * float(star["flux_err_mjy"][band]) / max(1e-9, f)


def _applied_ebv(star: dict[str, Any]) -> tuple[float, float | None, str]:
    ebv, src = star_applied_ebv(star)
    return (0.0 if ebv is None else float(ebv)), ebv, src


def dereddened_bp_rp(star: dict[str, Any]) -> float:
    ebv_used, _, _ = _applied_ebv(star)
    m = star["mags"]
    return float(m["BP"]) - float(m["RP"]) - colour_excess("BP", "RP", ebv_used)


def _locus_params(cfg_ch: dict[str, Any]) -> dict[str, Any]:
    return {
        "bp_rp_min": float(cfg_ch["locus_bp_rp_min"]),
        "bp_rp_max": float(cfg_ch["locus_bp_rp_max"]),
        "n_bins": int(cfg_ch["locus_n_bins"]),
        "min_stars": int(cfg_ch["locus_min_stars_per_bin"]),
        "clip_sigma": float(cfg_ch["locus_clip_sigma"]),
        "clip_max_iter": int(cfg_ch["locus_clip_max_iter"]),
        "scatter_floor": float(cfg_ch["locus_intrinsic_scatter_floor_mag"]),
    }


def _locus_id(locus: dict[str, Any]) -> str:
    core = {k: locus[k] for k in ("format_version", "source", "bands", "n_stars", "gaia_release", "sample_kind") if k in locus}
    return hashlib.sha256(json.dumps(core, sort_keys=True).encode("utf-8")).hexdigest()[:16]


# ------------------------------------------------------------------------------------------------ locus


def _clip_bin(vals: np.ndarray, k: float, max_iter: int) -> np.ndarray:
    """Iterative median/MAD sigma clipping; returns the boolean keep-mask."""
    keep = np.ones(len(vals), dtype=bool)
    for _ in range(max(1, max_iter)):
        med = float(np.median(vals[keep]))
        sig = 1.4826 * float(np.median(np.abs(vals[keep] - med)))
        if sig <= 0.0:
            break
        new_keep = np.abs(vals - med) <= k * sig
        if new_keep.sum() < 3 or np.array_equal(new_keep, keep):
            keep = new_keep if new_keep.sum() >= 3 else keep
            break
        keep = new_keep
    return keep


def build_empirical_locus(
    stars: list[dict[str, Any]],
    cfg_ch: dict[str, Any],
    *,
    meta: dict[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Dereddened (Ks - W1..W4) vs (BP - RP)_0 locus from sample stars, or None if no bin has enough stars.

    Per BP-RP bin and band: iterative `locus_clip_sigma` median/MAD clipping (so stars with excess do not bias the
    median), the median colour, the robust total scatter of the kept stars, the median photometric variance
    sigma_Ks^2 + sigma_W^2, and the intrinsic scatter sqrt(max(floor^2, scatter^2 - median variance)), i.e. the
    locus width beyond measurement noise (it contains real star-to-star and calibration scatter).
    """
    p = _locus_params(cfg_ch)
    rows: list[tuple[float, dict[str, float], dict[str, float]]] = []
    for s in stars:
        try:
            m = s["mags"]
            ebv_used, _, _ = _applied_ebv(s)
            bp_rp0 = float(m["BP"]) - float(m["RP"]) - colour_excess("BP", "RP", ebv_used)
            cols: dict[str, float] = {}
            var: dict[str, float] = {}
            for b in WISE_BANDS:
                cols[b] = (float(m["Ks"]) - float(m[b])) - colour_excess("Ks", b, ebv_used)
                var[b] = _mag_err(s, "Ks") ** 2 + _mag_err(s, b) ** 2
        except (KeyError, TypeError, ValueError):
            continue
        if all(math.isfinite(v) for v in (bp_rp0, *cols.values())):
            rows.append((bp_rp0, cols, var))
    if not rows:
        return None

    c_arr = np.asarray([r[0] for r in rows], dtype=float)
    edges = np.linspace(p["bp_rp_min"], p["bp_rp_max"], p["n_bins"] + 1)
    bands: dict[str, dict[str, list[float]]] = {
        b: {"bp_rp": [], "median": [], "scatter": [], "intrinsic_scatter": [], "median_phot_err": [], "n": [], "n_clipped": []}
        for b in WISE_BANDS
    }
    used_ids: set[int] = set()
    for i in range(len(edges) - 1):
        last = i == len(edges) - 2
        mask = (c_arr >= edges[i]) & ((c_arr <= edges[i + 1]) if last else (c_arr < edges[i + 1]))
        idx = np.nonzero(mask)[0]
        if len(idx) < p["min_stars"]:
            continue
        for b in WISE_BANDS:
            vals = np.asarray([rows[j][1][b] for j in idx], dtype=float)
            var = np.asarray([rows[j][2][b] for j in idx], dtype=float)
            keep = _clip_bin(vals, p["clip_sigma"], p["clip_max_iter"])
            if keep.sum() < 3:
                continue
            med = float(np.median(vals[keep]))
            scat = 1.4826 * float(np.median(np.abs(vals[keep] - med)))
            med_var = float(np.median(var[keep]))
            intrinsic = math.sqrt(max(p["scatter_floor"] ** 2, scat * scat - med_var))
            bk = bands[b]
            bk["bp_rp"].append(float(0.5 * (edges[i] + edges[i + 1])))
            bk["median"].append(med)
            bk["scatter"].append(scat)
            bk["intrinsic_scatter"].append(intrinsic)
            bk["median_phot_err"].append(math.sqrt(med_var))
            bk["n"].append(int(keep.sum()))
            bk["n_clipped"].append(int(len(idx) - keep.sum()))
            used_ids.update(int(j) for j in idx)
    if any(len(bands[b]["bp_rp"]) == 0 for b in WISE_BANDS):
        return None

    locus: dict[str, Any] = {
        "format_version": LOCUS_FORMAT_VERSION,
        "source": "sample_locus",
        "n_stars": len(used_ids),
        "dereddened": True,
        "colour_axis": "(BP-RP)_0 using ebv_star",
        "bp_rp_edges": [float(x) for x in edges],
        "bin_half_width": float(0.5 * (edges[1] - edges[0])),
        "clip_sigma": p["clip_sigma"],
        "clip_max_iter": p["clip_max_iter"],
        "intrinsic_scatter_floor_mag": p["scatter_floor"],
        "bands": bands,
        **(meta or {}),
    }
    locus["locus_id"] = _locus_id(locus)
    return locus


def literature_locus(cfg_ch: dict[str, Any]) -> dict[str, Any]:
    """Locus from the packaged Pecaut & Mamajek dwarf sequence (no sample needed).

    The table gives mean colours only; the per-band scatter is the assumed `literature_locus_scatter_mag`.
    Ks-W2/W3/W4 are tabulated only to K5V, so redder stars are clamped and flagged as extrapolated.
    """
    table = load_dwarf_colour_table()
    scat = [float(x) for x in cfg_ch["literature_locus_scatter_mag"]]
    bands: dict[str, dict[str, list[float]]] = {}
    for b, sc in zip(WISE_BANDS, scat):
        _, bp_rp, vals = dwarf_ks_minus_wise_sequence(b)
        order = np.argsort(bp_rp)
        bands[b] = {
            "bp_rp": [float(bp_rp[i]) for i in order],
            "median": [float(vals[i]) for i in order],
            "scatter": [sc] * len(order),
            "intrinsic_scatter": [sc] * len(order),
            "n": [0] * len(order),
        }
    locus: dict[str, Any] = {
        "format_version": LOCUS_FORMAT_VERSION,
        "source": literature_photosphere_label(),
        "n_stars": 0,
        "dereddened": True,
        "colour_axis": "intrinsic (BP-RP)_0",
        "bin_half_width": 0.0,
        "bands": bands,
        "reference": table.get("reference"),
        "url": table.get("url"),
        "table_version": table.get("table_version"),
        "retrieved_utc_date": table.get("retrieved_utc_date"),
        "provenance": table.get("provenance"),
        "scatter_note": "per-band scatter is an assumed configuration value (literature_locus_scatter_mag), not from the table",
    }
    locus["locus_id"] = _locus_id(locus)
    return locus


def predict_from_empirical_locus(
    star: dict[str, Any],
    locus: dict[str, Any],
    cfg_ch: dict[str, Any],
) -> dict[str, dict[str, Any]]:
    """Primary photosphere: predicted W1..W4 fluxes (observed frame, mJy) with uncertainties from the locus.

    sigma_pred(mag) = hypot(sigma_Ks, intrinsic locus scatter [, extrapolation term]); the prediction is
    Ks_0 - (Ks - W)_0 re-reddened with A_W, so the reddening applied to the star is used consistently.
    """
    ebv_used, _, _ = _applied_ebv(star)
    m = star["mags"]
    bp_rp0 = float(m["BP"]) - float(m["RP"]) - colour_excess("BP", "RP", ebv_used)
    ks0 = float(m["Ks"]) - extinction_mag("Ks", ebv_used)
    ks_err = _mag_err(star, "Ks")
    half = float(locus.get("bin_half_width", 0.0))
    extrap_sig = float(cfg_ch["locus_extrapolation_scatter_mag"])
    out: dict[str, dict[str, Any]] = {}
    for b in WISE_BANDS:
        bk = locus["bands"][b]
        x = np.asarray(bk["bp_rp"], dtype=float)
        ks_w0 = float(np.interp(bp_rp0, x, bk["median"]))
        intr = float(np.interp(bp_rp0, x, bk["intrinsic_scatter"]))
        extrap = bool(bp_rp0 < x[0] - half or bp_rp0 > x[-1] + half)
        sig_mag = math.hypot(ks_err, intr)
        if extrap:
            sig_mag = math.hypot(sig_mag, extrap_sig)
        w_pred_mag = ks0 - ks_w0 + extinction_mag(b, ebv_used)
        flux = mag_to_flux_mjy(b, w_pred_mag)
        out[b] = {
            "flux": flux,
            "sigma": flux * LN10_OVER_2P5 * sig_mag,
            "sigma_mag": sig_mag,
            "ks_w0": ks_w0,
            "intrinsic_scatter_mag": intr,
            "extrapolated": extrap,
        }
    return out


# ------------------------------------------------------------------------------------------------ secondary model


def fit_photosphere_model(star: dict[str, Any]) -> tuple[float, float, float, float, float]:
    """Fit Teff & Ks flux scale to dereddened G, BP, RP, J, H, Ks (Planck shape; `ebv_star` applied).

    Returns (teff_fit, scale_mjy [dereddened Ks scale], pred_w3_mjy, pred_w4_mjy [observed frame, secondary model
    with WISE colour correction and re-reddened], reduced chi2).
    """
    ebv_used, _, _ = _applied_ebv(star)
    flux = star["flux_mjy"]
    ferr = star["flux_err_mjy"]
    dered = np.array([10.0 ** (0.4 * extinction_mag(b, ebv_used)) for b in OPTICAL_NIR_BANDS], dtype=float)
    obs = np.array([float(flux[b]) for b in OPTICAL_NIR_BANDS], dtype=float) * dered
    sig = np.maximum(1e-4, np.array([float(ferr[b]) for b in OPTICAL_NIR_BANDS], dtype=float) * dered)
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
    pred = secondary_model_predictions(teff_fit, scale_fit, ebv_used)
    return teff_fit, scale_fit, pred["W3"], pred["W4"], red_chi2


def secondary_model_predictions(teff_k: float, scale_mjy: float, ebv: float, *, planck_only: bool = False) -> dict[str, float]:
    """Observed-frame W1..W4 fluxes of the secondary model (or of the bare Planck curve)."""
    fn = planck_photosphere_flux_mjy if planck_only else stellar_photosphere_flux_mjy
    return {
        b: fn(BAND_WAVELENGTH_UM[b], teff_k, scale_mjy) * 10.0 ** (-0.4 * extinction_mag(b, ebv))
        for b in WISE_BANDS
    }


# ------------------------------------------------------------------------------------------------ dust


def _dust_basis(teff_fit: float, scale_fit: float, t_d: float) -> np.ndarray:
    return np.array(
        [dust_blackbody_flux_mjy(BAND_WAVELENGTH_UM[b], teff_fit, scale_fit, t_d, 1.0) for b in WISE_BANDS],
        dtype=float,
    )


def _colour_temperature(
    excess: dict[str, tuple[float, float]],
    b1: str,
    b2: str,
    k_sig: float,
    t_min: float,
    t_max: float,
) -> dict[str, float]:
    e1, s1 = excess[b1]
    e2, s2 = excess[b2]
    key = f"t_color_{b1.lower()}{b2.lower()}"
    if not (e1 > 0.0 and e2 > 0.0 and e1 / s1 >= k_sig and e2 / s2 >= k_sig):
        return {f"{key}_k": 0.0, f"{key}_lo_k": 0.0, f"{key}_hi_k": 0.0, f"{key}_valid": 0.0, f"{key}_hi_unbounded": 0.0}
    r = e2 / e1
    sr = r * math.hypot(s1 / e1, s2 / e2)
    t, st = blackbody_colour_temperature_k(b1, b2, r, t_min_k=t_min, t_max_k=t_max)
    t_lo, _ = blackbody_colour_temperature_k(b1, b2, r + sr, t_min_k=t_min, t_max_k=t_max)
    if r - sr > 0.0:
        t_hi, st_hi = blackbody_colour_temperature_k(b1, b2, r - sr, t_min_k=t_min, t_max_k=t_max)
    else:
        t_hi, st_hi = t_max, "above_max"
    return {
        f"{key}_k": float(t),
        f"{key}_lo_k": float(t_lo),
        f"{key}_hi_k": float(t_hi),
        f"{key}_valid": 1.0 if st == "ok" else 0.0,
        f"{key}_hi_unbounded": 1.0 if st_hi != "ok" else 0.0,
    }


def fit_excess_blackbody_and_shell(
    star: dict[str, Any],
    teff_fit: float,
    scale_fit: float,
    cfg_ch: dict[str, Any],
    excess: dict[str, tuple[float, float]] | None = None,
) -> dict[str, float]:
    """Single-temperature dust blackbody (with a Delta-chi2 temperature interval), W1-W2 and W3-W4 excess colour
    temperatures, a two-temperature fractional luminosity, and the Suazo et al. (2024) partial-shell grid.

    `excess` maps band -> (excess flux mJy, 1-sigma); when omitted it is computed against the secondary model
    (kept for callers of the old signature). The blackbody is fitted only to bands whose excess is significant
    at >= `excess_significance_sigma`; at least two such bands are needed to constrain a temperature.
    """
    if excess is None:
        ebv_used, _, _ = _applied_ebv(star)
        pred = secondary_model_predictions(teff_fit, scale_fit, ebv_used)
        excess = {}
        for b in WISE_BANDS:
            e_obs = float(star["flux_err_mjy"][b])
            excess[b] = (float(star["flux_mjy"][b]) - pred[b], math.hypot(e_obs, 0.03 * pred[b]))

    k_sig = float(cfg_ch["excess_significance_sigma"])
    t_min = float(cfg_ch.get("blackbody_t_min_k", 80.0))
    t_max = float(cfg_ch.get("blackbody_t_max_k", 1200.0))
    d_chi2 = float(cfg_ch["dust_temp_interval_delta_chi2"])
    lum = float(star.get("luminosity_lsun", 1.0) or 1.0)

    e_all = np.asarray([excess[b][0] for b in WISE_BANDS], dtype=float)
    s_all = np.asarray([max(1e-6, excess[b][1]) for b in WISE_BANDS], dtype=float)
    sig_mask = (e_all > 0.0) & (e_all / s_all >= k_sig)
    n_fit = int(sig_mask.sum())

    res: dict[str, float] = {
        "t_dust_k": 0.0,
        "t_dust_lo_k": 0.0,
        "t_dust_hi_k": 0.0,
        "t_dust_at_bound": 0.0,
        "frac_lum": 0.0,
        "frac_lum_lo": 0.0,
        "frac_lum_hi": 0.0,
        "r_bb_au": 0.0,
        "dust_fit_n_bands": float(n_fit),
        "dust_fit_chi2": 0.0,
        "dust_fit_red_chi2": 0.0,
        "dust_fit_interval_delta_chi2": 0.0,
    }
    if n_fit >= 2:
        e = e_all[sig_mask]
        w = 1.0 / (s_all[sig_mask] ** 2)

        def _prof(t_d: float) -> tuple[float, float, float]:
            g = _dust_basis(teff_fit, scale_fit, t_d)[sig_mask]
            gg = float(np.sum(w * g * g))
            if gg <= 0.0:
                return 0.0, float("inf"), 0.0
            f = max(0.0, float(np.sum(w * e * g)) / gg)
            chi2 = float(np.sum(w * (e - f * g) ** 2))
            return f, chi2, 1.0 / math.sqrt(gg)

        grid = np.geomspace(t_min, t_max, 241)
        prof = [_prof(float(t)) for t in grid]
        chi = np.asarray([p[1] for p in prof], dtype=float)
        i0 = int(np.argmin(chi))
        lo_t = float(grid[max(0, i0 - 1)])
        hi_t = float(grid[min(len(grid) - 1, i0 + 1)])
        ref = minimize_scalar(lambda t: _prof(float(t))[1], bounds=(lo_t, hi_t), method="bounded", options={"xatol": 0.5})
        t_best = float(ref.x) if ref.fun <= chi[i0] else float(grid[i0])
        f_best, chi_best, sf_best = _prof(t_best)
        # A poor single-temperature fit (reduced chi2 > 1) widens the interval by the reduced chi2 (error scaling),
        # so the quoted interval is not artificially narrow when one blackbody cannot describe all bands.
        dof = n_fit - 2
        red = chi_best / dof if dof > 0 else 0.0
        width = d_chi2 * max(1.0, red)
        ok = chi <= chi_best + width
        res["dust_fit_red_chi2"] = red
        res["dust_fit_interval_delta_chi2"] = width
        # contiguous region around the best grid point
        j_lo = i0
        while j_lo > 0 and ok[j_lo - 1]:
            j_lo -= 1
        j_hi = i0
        while j_hi < len(grid) - 1 and ok[j_hi + 1]:
            j_hi += 1
        f_lo = min([f_best - sf_best] + [prof[j][0] - prof[j][2] for j in range(j_lo, j_hi + 1)])
        f_hi = max([f_best + sf_best] + [prof[j][0] + prof[j][2] for j in range(j_lo, j_hi + 1)])
        res.update(
            {
                "t_dust_k": t_best,
                "t_dust_lo_k": float(grid[j_lo]),
                "t_dust_hi_k": float(grid[j_hi]),
                "t_dust_at_bound": 1.0 if (j_lo == 0 or j_hi == len(grid) - 1) else 0.0,
                "frac_lum": f_best,
                "frac_lum_lo": max(0.0, f_lo),
                "frac_lum_hi": f_hi,
                "r_bb_au": equilibrium_radius_au(t_best, lum) if f_best > 0.0 else 0.0,
                "dust_fit_chi2": chi_best,
            }
        )

    # Colour temperatures of the W1-W2 and W3-W4 excess (a single blackbody need not describe both)
    ct_min = float(cfg_ch["colour_temp_min_k"])
    ct_max = float(cfg_ch["colour_temp_max_k"])
    ct12 = _colour_temperature(excess, "W1", "W2", k_sig, ct_min, ct_max)
    ct34 = _colour_temperature(excess, "W3", "W4", k_sig, ct_min, ct_max)
    res.update(ct12)
    res.update(ct34)

    # Two-temperature (hot W1-W2 + warm W3-W4) fractional luminosity, non-negative linear fit to all four bands
    res.update({"frac_lum_hot": 0.0, "frac_lum_warm": 0.0, "frac_lum_two_temp": 0.0})
    if ct12["t_color_w1w2_valid"] and ct34["t_color_w3w4_valid"]:
        a = np.vstack(
            [
                _dust_basis(teff_fit, scale_fit, ct12["t_color_w1w2_k"]),
                _dust_basis(teff_fit, scale_fit, ct34["t_color_w3w4_k"]),
            ]
        ).T / s_all[:, None]
        coef, _ = nnls(a, np.maximum(0.0, e_all) / s_all)
        res.update({"frac_lum_hot": float(coef[0]), "frac_lum_warm": float(coef[1]), "frac_lum_two_temp": float(coef.sum())})

    # Partial-shell grid search (Suazo et al. 2024 covering fraction gamma & temperature) on all four bands
    gamma_grid = cfg_ch.get("shell_gamma_grid", [0.001, 0.01, 0.03, 0.1, 0.3])
    temp_grid = cfg_ch.get("shell_temp_grid_k", [150.0, 200.0, 250.0, 300.0, 400.0])
    e_pos = np.maximum(0.0, e_all)
    best_chi2 = float("inf")
    best_gamma = 0.0
    best_temp = 0.0
    for t_s in temp_grid:
        g = _dust_basis(teff_fit, scale_fit, float(t_s))
        for gamma in gamma_grid:
            chi2 = float(np.sum(((e_pos - float(gamma) * g) / s_all) ** 2))
            if chi2 < best_chi2:
                best_chi2 = chi2
                best_gamma = float(gamma)
                best_temp = float(t_s)
    res.update({"shell_gamma": best_gamma, "shell_temp_k": best_temp, "shell_grid_chi2": best_chi2})
    return res


# ------------------------------------------------------------------------------------------------ channel


class IrExcessChannel:
    """Mid-infrared excess channel implementing the unified Channel protocol."""

    name = "ir_excess"

    def __init__(self, cfg: PipelineConfig, root: Path) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["ir_excess"]
        self.root = root
        self._locus: dict[str, Any] | None = None
        self._fallback_cache: dict[tuple[str, str], dict[str, Any]] = {}
        self.last_persisted_locus_path: Path | None = None

    # ---- locus persistence -------------------------------------------------------------------------------

    def _loci_dir(self) -> Path:
        return self.root / self.cfg.storage.runs_dir / LOCUS_DIRNAME

    def persist_locus(self, locus: dict[str, Any]) -> Path | None:
        """Write a sample locus to runs/_photosphere_loci/ (idempotent per locus content; refreshes `last_built_utc`)."""
        if locus.get("source") != "sample_locus":
            return None
        d = self._loci_dir()
        name = f"ir_excess_locus-{locus.get('gaia_release', 'na')}-{locus.get('sample_kind', 'na')}-{locus['locus_id']}.json"
        payload = dict(locus)
        payload["source"] = "persisted_locus"
        payload["original_source"] = "sample_locus"
        payload["last_built_utc"] = datetime.now(timezone.utc).isoformat()
        try:
            d.mkdir(parents=True, exist_ok=True)
            path = d / name
            path.write_text(json.dumps(payload, indent=1, sort_keys=True) + "\n", encoding="utf-8")
        except OSError:
            return None
        self.last_persisted_locus_path = path
        return path

    def load_persisted_locus(self, gaia_release: str, sample_kind: str) -> dict[str, Any] | None:
        """Most recent persisted locus (by `last_built_utc`) for this Gaia release and sample kind."""
        d = self._loci_dir()
        if not d.exists():
            return None
        best: tuple[str, str] | None = None
        best_obj: dict[str, Any] | None = None
        for p in d.glob(f"ir_excess_locus-{gaia_release}-{sample_kind}-*.json"):
            try:
                obj = json.loads(p.read_text(encoding="utf-8"))
            except (OSError, ValueError):
                continue
            if obj.get("format_version") != LOCUS_FORMAT_VERSION or obj.get("gaia_release") != gaia_release:
                continue
            key = (str(obj.get("last_built_utc", "")), p.name)
            if best is None or key > best:
                best, best_obj = key, obj
        if best_obj is not None:
            best_obj["source"] = "persisted_locus"
            best_obj["persisted_file"] = f"{self.cfg.storage.runs_dir}/{LOCUS_DIRNAME}/{best[1]}" if best else None
        return best_obj

    def _fallback_locus(self, gaia_release: str, sample_kind: str) -> dict[str, Any]:
        key = (gaia_release, sample_kind)
        if key not in self._fallback_cache:
            self._fallback_cache[key] = self.load_persisted_locus(gaia_release, sample_kind) or literature_locus(self.ch_cfg)
        return self._fallback_cache[key]

    def set_empirical_locus(self, stars: list[dict[str, Any]]) -> dict[str, Any]:
        """Build the sample locus (persisting it); fall back to a persisted, then the literature, locus."""
        release = str(stars[0].get("gaia_release", self.cfg.gaia.release)) if stars else str(self.cfg.gaia.release)
        kind = "live" if any(star_sample_kind(s) == "live" for s in stars) else "synthetic"
        meta = {
            "gaia_release": release,
            "sample_kind": kind,
            "config_hash": compute_config_hash(self.cfg.model_dump())[:12],
            "sample_size": len(stars),
        }
        locus = build_empirical_locus(stars, self.ch_cfg, meta=meta)
        if locus is None:
            locus = self._fallback_locus(release, kind)
        elif locus["n_stars"] >= int(self.ch_cfg["locus_persist_min_stars"]):
            self.persist_locus(locus)
        self._locus = locus
        return locus

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

    # ---- measurement -------------------------------------------------------------------------------------

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        locus = data.get("empirical_locus") or self._locus
        if locus is None:
            locus = self._fallback_locus(str(star.get("gaia_release", self.cfg.gaia.release)), star_sample_kind(star))
        phot_source = str(locus.get("source", "unknown"))
        cal_frac = float(self.ch_cfg.get("cal_uncertainty_frac", 0.05))
        ebv_used, ebv_val, ebv_src = _applied_ebv(star)

        # Secondary model (Planck shape + WISE colour correction), and the bare Planck curve as a diagnostic
        teff_fit, scale_fit, _, _, red_chi2 = fit_photosphere_model(star)
        mod = secondary_model_predictions(teff_fit, scale_fit, ebv_used)
        planck = secondary_model_predictions(teff_fit, scale_fit, ebv_used, planck_only=True)
        # Primary: empirical locus
        emp = predict_from_empirical_locus(star, locus, self.ch_cfg)

        obs = {b: float(star["flux_mjy"][b]) for b in WISE_BANDS}
        sig_obs = {b: float(star["flux_err_mjy"][b]) for b in WISE_BANDS}

        values: dict[str, float] = {
            "teff_fit_k": teff_fit,
            "star_scale_mjy": scale_fit,
            "photosphere_fit_red_chi2": red_chi2,
            "ebv_star_applied": ebv_used,
            "bp_rp_0": dereddened_bp_rp(star),
            "locus_n_stars": float(locus.get("n_stars", 0) or 0),
        }
        uncertainties: dict[str, float] = {}
        excess: dict[str, tuple[float, float]] = {}
        any_extrap = False
        for b in WISE_BANDS:
            lb = b.lower()
            ph, sph = emp[b]["flux"], emp[b]["sigma"]
            exc = obs[b] - ph
            sexc = math.hypot(sig_obs[b], sph)
            excess[b] = (exc, sexc)
            any_extrap = any_extrap or bool(emp[b]["extrapolated"])
            values[f"{lb}_obs_mjy"] = obs[b]
            values[f"{lb}_phot_mjy"] = ph
            values[f"{lb}_excess_mjy"] = exc
            values[f"{lb}_excess_sigma"] = exc / max(1e-9, sexc)
            values[f"{lb}_excess_ratio"] = obs[b] / max(1e-6, ph)
            values[f"{lb}_pred_model_mjy"] = mod[b]
            values[f"{lb}_excess_model_mjy"] = obs[b] - mod[b]
            values[f"{lb}_pred_planck_mjy"] = planck[b]
            values[f"{lb}_locus_extrapolated"] = 1.0 if emp[b]["extrapolated"] else 0.0
            uncertainties[f"{lb}_obs_mjy"] = sig_obs[b]
            uncertainties[f"{lb}_phot_mjy"] = sph
            uncertainties[f"{lb}_excess_mjy"] = sexc

        # Two-predictor detection statistics (structure unchanged): W3/W4 against each predictor
        sig_pred_mod = {"W3": 0.025 * mod["W3"], "W4": 0.030 * mod["W4"]}
        for b in ("W3", "W4"):
            lb = b.lower()
            chi_mod = (obs[b] - mod[b]) / math.hypot(sig_obs[b], math.hypot(sig_pred_mod[b], cal_frac * mod[b]))
            # Empirical leg: the locus intrinsic scatter already contains the star-to-star calibration scatter of
            # the sample, so no separate calibration term is added (it would be counted twice).
            chi_emp = excess[b][0] / max(1e-9, excess[b][1])
            values.update(
                {
                    f"{lb}_pred_mjy": mod[b],
                    f"{lb}_pred_emp_mjy": emp[b]["flux"],
                    f"{lb}_ratio": obs[b] / max(1e-6, mod[b]),
                    f"{lb}_ratio_emp": obs[b] / max(1e-6, emp[b]["flux"]),
                    f"chi_{lb}": chi_mod,
                    f"chi_emp_{lb}": chi_emp,
                }
            )
            uncertainties[f"{lb}_pred_mjy"] = sig_pred_mod[b]
            uncertainties[f"{lb}_pred_emp_mjy"] = emp[b]["sigma"]

        # Dust quantification on the PRIMARY (empirical) excess
        dust = fit_excess_blackbody_and_shell(star, teff_fit, scale_fit, self.ch_cfg, excess=excess)
        values.update(dust)
        uncertainties["frac_lum"] = 0.5 * (dust["frac_lum_hi"] - dust["frac_lum_lo"])
        uncertainties["t_dust_k"] = 0.5 * (dust["t_dust_hi_k"] - dust["t_dust_lo_k"])

        # Secondary-model dust fit, diagnostic only
        excess_mod = {b: (obs[b] - mod[b], math.hypot(sig_obs[b], 0.03 * mod[b])) for b in WISE_BANDS}
        dust_mod = fit_excess_blackbody_and_shell(star, teff_fit, scale_fit, self.ch_cfg, excess=excess_mod)
        values["t_dust_model_k"] = dust_mod["t_dust_k"]
        values["frac_lum_model"] = dust_mod["frac_lum"]

        vcfg = self.cfg.vetting
        ext_cases = extinction_two_case_summary(
            star,
            min_los_excess_ebv=float(vcfg["extinction_los_excess_summary_ebv"]),
            sfd_to_sf11=float(vcfg["sfd_to_sf11_factor"]),
        )

        caveats = [SINGLE_EPOCH_CAVEAT]
        if phot_source.startswith("literature_fallback"):
            caveats.append(
                f"Photosphere from the packaged dwarf colour table ({phot_source}); no sample or persisted locus was "
                "available, so the locus scatter is an assumed value."
            )
        if any_extrap:
            caveats.append("Star's (BP-RP)_0 lies outside the locus coverage in at least one band; prediction clamped and uncertainty inflated.")
        if ebv_val is None:
            caveats.append("No reddening estimate for this star (ebv_star missing): photometry used without dereddening.")
        if ext_cases is not None:
            caveats.append(f"Reddening uncertain: {ext_cases['triggered_by']}.")

        extra = {
            "gaia_release": star.get("gaia_release", self.cfg.gaia.release),
            "photosphere_source": phot_source,
            "photosphere_secondary": SECONDARY_PHOTOSPHERE,
            "locus_n_stars": int(locus.get("n_stars", 0) or 0),
            "locus_id": locus.get("locus_id"),
            "locus_gaia_release": locus.get("gaia_release"),
            "locus_sample_kind": locus.get("sample_kind"),
            "locus_config_hash": locus.get("config_hash"),
            "ebv_star": ebv_val,
            "ebv_star_source": ebv_src,
            "extinction_cases": ext_cases,
        }
        if phot_source.startswith("literature_fallback"):
            extra["literature_reference"] = locus.get("reference")
            extra["literature_url"] = locus.get("url")
            extra["literature_retrieved_utc_date"] = locus.get("retrieved_utc_date")
        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=values,
            uncertainties=uncertainties,
            epochs=["2010-03"],
            caveats=caveats,
            extra=extra,
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
        src = str(m.extra.get("photosphere_source", "unknown"))
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=bool(flagged),
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=(
                f"W4_obs/W4_phot={v['w4_ratio_emp']:.2f} (chi_emp_w4={v['chi_emp_w4']:.2f}, model chi_w4={v['chi_w4']:.2f}, "
                f"f_dust={v['frac_lum']:.4f}, T_dust={v['t_dust_k']:.0f}K [{v['t_dust_lo_k']:.0f}-{v['t_dust_hi_k']:.0f}]; "
                f"photosphere={src}, N_locus={int(v.get('locus_n_stars', 0))})"
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
                description=(
                    f"Primary photosphere ({m.extra.get('photosphere_source', 'unknown')}) and secondary Planck-based model, "
                    "per-band excess with uncertainties, chi statistics, and blackbody/shell/colour-temperature fits"
                ),
                sha256=sha,
            )
        ]
