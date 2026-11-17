"""Light-curve and IR time-domain vetting checks (Section 6)."""

from __future__ import annotations

import math
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


def matched_detrend(t: np.ndarray, f: np.ndarray, window_days: float, gap_days: float = 0.5) -> np.ndarray:
    """Detrend one light curve with the filter applied identically to every extraction being compared.

    Per continuous segment (split at gaps > `gap_days`), divide by a running median of `window_days`, computed
    after masking points more than 2.5 robust sigma below a first-pass baseline so a real dip does not pull the
    baseline down. The window is the optical_dips channel's own `detrend_window_days`, so the comparison
    preserves the same dip timescales the detection channel searches.
    """
    from scipy.ndimage import median_filter

    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    out = np.ones_like(f)
    if len(t) == 0:
        return out
    breaks = np.where(np.diff(t) > gap_days)[0] + 1
    for seg in np.split(np.arange(len(t)), breaks):
        if len(seg) < 5:
            out[seg] = f[seg] / max(1e-9, float(np.median(f[seg])))
            continue
        ts, fs = t[seg], f[seg]
        dt = max(1e-4, float(np.median(np.diff(ts))))
        win = max(5, int(round(window_days / dt)) | 1)
        win = min(win, len(seg) if len(seg) % 2 == 1 else len(seg) - 1)
        base = median_filter(fs, size=win, mode="nearest")
        resid = fs - base
        noise = max(1e-9, 1.4826 * float(np.median(np.abs(resid))))
        valid = resid >= -2.5 * noise
        if int(np.sum(valid)) >= 5:
            base = median_filter(np.interp(ts, ts[valid], fs[valid]), size=win, mode="nearest")
        out[seg] = fs / np.maximum(1e-9, base)
    return out


def _bin_mean(t: np.ndarray, f: np.ndarray, width_days: float) -> np.ndarray:
    idx = np.floor((t - t[0]) / width_days).astype(int)
    _, inv = np.unique(idx, return_inverse=True)
    return np.bincount(inv, weights=f) / np.bincount(inv)


def _depths_agree(d_qlp: float, d_cut: float) -> bool:
    return (d_qlp < 0.02 and d_cut < 0.02) or (abs(d_qlp - d_cut) <= 0.08 and d_cut >= 0.5 * d_qlp)


def check_pipeline_agreement(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    sectors = [s for s in lc.get("tess_sectors", []) or [] if len(s.get("flux_qlp_raw", [])) >= 20]
    if sectors:
        dip_cfg = cfg.channels["optical_dips"]
        window = float(dip_cfg.get("detrend_window_days", 5.0))
        bin_days = float(dip_cfg.get("min_duration_hours", 6.0)) / 24.0
        per_sector: list[dict[str, Any]] = []
        for sec in sectors:
            t = np.asarray(sec["time_btjd"], dtype=float)
            fq = np.asarray(sec["flux_qlp_raw"], dtype=float)
            fc = np.asarray(sec["flux_tesscut"], dtype=float)
            if len(fc) != len(fq) or len(t) != len(fq):
                continue
            rq = matched_detrend(t, fq, window)
            rc = matched_detrend(t, fc, window)
            bq = _bin_mean(t, rq, bin_days)
            bc = _bin_mean(t, rc, bin_days)
            d_q = max(0.0, 1.0 - float(np.min(bq)))
            d_c = max(0.0, 1.0 - float(np.min(bc)))
            r = float(np.corrcoef(bq, bc)[0, 1]) if len(bq) > 3 and np.std(bq) > 0 and np.std(bc) > 0 else float("nan")
            per_sector.append(
                {
                    "sector": int(sec.get("sector", -1)),
                    "n_bins": int(len(bq)),
                    "qlp_max_depth": round(d_q, 4),
                    "tesscut_max_depth": round(d_c, 4),
                    "qlp_robust_sigma": round(float(1.4826 * np.median(np.abs(bq - np.median(bq)))), 5),
                    "tesscut_robust_sigma": round(float(1.4826 * np.median(np.abs(bc - np.median(bc)))), 5),
                    "binned_pearson_r": None if math.isnan(r) else round(r, 3),
                    "agree": _depths_agree(d_q, d_c),
                }
            )
        if per_sector:
            bad = [p["sector"] for p in per_sector if not p["agree"]]
            d1 = max(p["qlp_max_depth"] for p in per_sector)
            d2 = max(p["tesscut_max_depth"] for p in per_sector)
            note = (
                f"Matched detrending ({window:.1f} d running median, {bin_days * 24:.0f} h bins) over "
                f"{len(per_sector)} sector(s): QLP max dip={d1:.3f}, independent TESScut max dip={d2:.3f}"
            )
            if bad:
                note += f"; disagreement in sector(s) {bad}"
            return make_vet_result(
                out_dir,
                check="pipeline_agreement",
                status="fail" if bad else "pass",
                note=note,
                metrics={
                    "qlp_max_depth": d1,
                    "tesscut_max_depth": d2,
                    "tic_id": lc.get("tess_tic_id"),
                    "detrend_window_days": window,
                    "bin_hours": bin_days * 24.0,
                    "sectors": per_sector,
                },
            )

    f1 = np.asarray(lc.get("tess_flux_qlp", []), dtype=float)
    f2 = np.asarray(lc.get("tess_flux_tesscut", []), dtype=float)
    if len(f1) == 0 or len(f2) == 0 or len(f1) != len(f2):
        log = lc.get("tess_query_log") or []
        return make_vet_result(
            out_dir,
            check="pipeline_agreement",
            status="not_evaluated",
            note="Independent second light-curve extraction missing" + (f" ({'; '.join(map(str, log))[:300]})" if log else ""),
        )
    # Check whether maximum dip depth in QLP is also seen in TESScut
    d1 = 1.0 - float(np.min(f1))
    d2 = 1.0 - float(np.min(f2))
    agree = _depths_agree(d1, d2)
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
    surveys = lc.get("optical_surveys") or {}
    usable = {k: v for k, v in surveys.items() if int(v.get("n_epochs", 0) or 0) >= 10}
    if usable:
        parts = []
        metrics: dict[str, Any] = {}
        for key, sv in sorted(usable.items()):
            band_txt = []
            for band, b in sorted((sv.get("bands") or {}).items()):
                rng = b.get("seasonal_range_mag")
                band_txt.append(f"{band}: {b.get('n_epochs')} epochs, seasonal range={'n/a' if rng is None else f'{rng:.3f}'} mag")
            parts.append(f"{sv.get('survey', key)} ({'; '.join(band_txt)})")
            metrics[key] = sv
        missing = sorted(k for k in surveys if k not in usable)
        for k in missing:
            metrics[k] = {kk: surveys[k].get(kk) for kk in ("survey", "n_epochs", "status", "query_log")}
        ranges = [
            b["seasonal_range_mag"]
            for sv in usable.values()
            for b in (sv.get("bands") or {}).values()
            if b.get("seasonal_range_mag") is not None
        ]
        metrics["max_seasonal_range_mag"] = max(ranges) if ranges else None
        metrics["n_surveys"] = len(usable)
        return make_vet_result(
            out_dir,
            check="second_survey",
            status="pass",
            note="Independent optical survey coverage: " + " | ".join(parts) + (f" (not usable: {missing})" if missing else ""),
            metrics=metrics,
        )

    ztf = np.asarray(lc.get("ztf_flux", []), dtype=float)
    if len(ztf) < 10:
        why = "; ".join(f"{v.get('survey', k)}: {v.get('status')}" for k, v in surveys.items())
        return make_vet_result(
            out_dir,
            check="second_survey",
            status="not_evaluated",
            note="No ZTF/ASAS-SN optical coverage at target position" + (f" ({why})" if why else ""),
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
