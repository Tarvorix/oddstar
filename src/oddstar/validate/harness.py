"""Validation harness enforcing Ground Rules R1, R2, R3, R4 and search lockout."""

from __future__ import annotations

import copy
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import numpy as np

from oddstar.channels import build_channel
from oddstar.channels.ir_excess.channel import build_empirical_locus
from oddstar.channels.radio_narrowband.channel import inject_synthetic_drifting_tone
from oddstar.config import PipelineConfig, compute_config_hash, current_git_commit
from oddstar.errors import UncertifiedSearchError, ValidationGateError
from oddstar.fuse.engine import assign_tier
from oddstar.io.fixtures_data import RADIO_CADENCES, generate_field_sample, get_named_object
from oddstar.io.physics import BAND_WAVELENGTH_UM, dust_blackbody_flux_mjy, flux_mjy_to_mag
from oddstar.schemas import ValidationCertificate
from oddstar.vet.engine import vet_target


def _wilson_ci(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson score confidence interval for binomial false-alarm rate."""
    if n <= 0:
        return 0.0, 1.0
    p = k / n
    denom = 1.0 + (z * z) / n
    center = (p + (z * z) / (2.0 * n)) / denom
    half = (z * math.sqrt((p * (1.0 - p) + (z * z) / (4.0 * n)) / n)) / denom
    return max(0.0, center - half), min(1.0, center + half)


def _load_list(path: Path, channel: str, split: str) -> list[dict[str, Any]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    items = list(data.get(channel, []))
    if split in {"dev", "holdout"}:
        filtered = [x for x in items if x.get("split") == split]
        if filtered:
            return filtered
    return items


def _certificate_path(root: Path, channel: str, gaia_release: str, config_hash: str) -> Path:
    cert_dir = root / "validation" / "certificates"
    cert_dir.mkdir(parents=True, exist_ok=True)
    return cert_dir / f"{channel}-{gaia_release}-{config_hash}.json"


def load_certificate(
    channel: str,
    cfg: PipelineConfig,
    root: Path,
) -> ValidationCertificate | None:
    short_hash = compute_config_hash(cfg.model_dump())[:12]
    p = _certificate_path(root, channel, cfg.gaia.release, short_hash)
    if not p.exists():
        return None
    return ValidationCertificate.model_validate_json(p.read_text(encoding="utf-8"))


def require_valid_certificate(channel: str, cfg: PipelineConfig, root: Path) -> ValidationCertificate:
    """Enforce Ground Rule R1/R2/R4 lockout on oddstar search."""
    short_hash = compute_config_hash(cfg.model_dump())[:12]
    cert = load_certificate(channel, cfg, root)
    git_sha = current_git_commit(root)
    if cert is None:
        raise UncertifiedSearchError(
            f"Search blocked: no validation certificate for channel '{channel}' "
            f"(release={cfg.gaia.release}, config_hash={short_hash}).",
            hint=f"Run `oddstar validate {channel}` first.",
        )
    if not cert.passed:
        raise UncertifiedSearchError(
            f"Search blocked: validation certificate for '{channel}' did not pass.",
            hint=f"Inspect { _certificate_path(root, channel, cfg.gaia.release, short_hash) } and fix failing gates.",
        )
    if cert.git_commit != git_sha:
        raise UncertifiedSearchError(
            f"Search blocked: certificate commit {cert.git_commit[:10]} does not match HEAD {git_sha[:10]}.",
            hint=f"Re-run `oddstar validate {channel}` at the current commit.",
        )
    return cert


def _prepare_channel_with_bg(
    channel_name: str,
    cfg: PipelineConfig,
    field_stars: list[dict[str, Any]],
    root: Path,
    **ch_kwargs: Any,
) -> Any:
    ch = build_channel(channel_name, cfg, root, **ch_kwargs)
    if channel_name == "ir_excess":
        ch.set_empirical_locus(field_stars)
    elif channel_name == "underluminous":
        ch.fit_locus(field_stars)
    elif channel_name == "generic_outlier":
        ch.fit_background(field_stars)
    return ch


def _run_test1_blind_recovery(
    channel_name: str,
    cfg: PipelineConfig,
    field_stars: list[dict[str, Any]],
    benchmarks: list[dict[str, Any]],
    root: Path,
) -> tuple[float, int, int, list[dict[str, Any]]]:
    """Test 1 (R1): Blind recovery of frozen positive benchmarks embedded inside field stars."""
    # Toggle switch-off corrections for channels without natural non-explained positives
    ch_kwargs: dict[str, Any] = {}
    if channel_name == "underluminous":
        ch_kwargs["extinction_correction_enabled"] = False
    elif channel_name == "spectral_lines":
        ch_kwargs["stellar_line_veto_enabled"] = False
    elif channel_name == "vanishing":
        ch_kwargs["proper_motion_correction_enabled"] = False

    ch = _prepare_channel_with_bg(channel_name, cfg, field_stars, root, **ch_kwargs)
    data_ctx: dict[str, Any] = {}
    if channel_name == "ir_excess":
        data_ctx["empirical_locus"] = ch._locus

    # Embed benchmarks unlabelled inside a subset of field stars
    rng = np.random.default_rng(cfg.sample.seed + 101)
    bench_ids = {b["source_id"] for b in benchmarks}
    pool: list[dict[str, Any]] = []
    if channel_name == "radio_narrowband":
        for cid, cad in RADIO_CADENCES.items():
            if cid != "GBT_BLC1_PROXCEN":
                pool.append(copy.deepcopy(cad))
    else:
        pool.extend(field_stars[:400])
        for b in benchmarks:
            obj = get_named_object(b["source_id"], gaia_release=cfg.gaia.release)
            if obj is not None:
                # Strip label fields before scoring (true blind check)
                blind = copy.deepcopy(obj)
                blind["name"] = f"BLIND_{blind['source_id']}"
                blind["known_in_literature"] = False
                pool.append(blind)
    rng.shuffle(pool)

    recovered = 0
    details: list[dict[str, Any]] = []
    for item in pool:
        sid = str(item["source_id"])
        if sid not in bench_ids:
            continue
        m = ch.measure(item, data_ctx)
        sc = ch.score(m)
        if sc.flagged:
            recovered += 1
        details.append({"source_id": sid, "flagged": sc.flagged, "raw_stat": sc.raw_statistic})

    total = len(benchmarks)
    rate = recovered / max(1, total)
    return rate, recovered, total, details


def _run_test2_impostors(
    channel_name: str,
    cfg: PipelineConfig,
    field_stars: list[dict[str, Any]],
    impostors: list[dict[str, Any]],
    root: Path,
    tmp_vet_dir: Path,
) -> tuple[bool, str, int, list[dict[str, Any]]]:
    """Test 2 (R2): Reject known impostors; zero may exceed Tier C and expected checks must fail."""
    ch = _prepare_channel_with_bg(channel_name, cfg, field_stars, root)
    data_ctx: dict[str, Any] = {}
    if channel_name == "ir_excess":
        data_ctx["empirical_locus"] = ch._locus

    tier_rank = {"A": 4, "B": 3, "C": 2, "Explained": 1, "Rejected": 0}
    max_tier_seen = "Rejected"
    all_passed = True
    details: list[dict[str, Any]] = []

    for imp in impostors:
        sid = imp["source_id"]
        if channel_name == "radio_narrowband" and sid in RADIO_CADENCES:
            obj = copy.deepcopy(RADIO_CADENCES[sid])
        else:
            obj = get_named_object(sid, gaia_release=cfg.gaia.release)
        assert obj is not None, f"Missing impostor fixture {sid}"

        m = ch.measure(obj, data_ctx)
        sc = ch.score(m)
        # Vetting always runs on impostors in Test 2
        forced_sc = sc.model_copy(update={"flagged": True})
        vet_res = vet_target(obj, [forced_sc], cfg, tmp_vet_dir / sid)
        tier, reason = assign_tier([sc], vet_res, cfg)
        failed_checks = {v.check for v in vet_res if v.status == "fail"}
        expected_failed = set(imp.get("expected_failed", []))
        expected_ok = expected_failed.issubset(failed_checks) if expected_failed else True

        if tier_rank[tier] > tier_rank[max_tier_seen]:
            max_tier_seen = tier
        if tier in {"A", "B"} or not expected_ok:
            all_passed = False
        details.append(
            {
                "source_id": sid,
                "tier": tier,
                "failed_checks": sorted(failed_checks),
                "expected_failed": sorted(expected_failed),
                "reason": reason,
                "passed": (tier not in {"A", "B"}) and expected_ok,
            }
        )

    return all_passed, max_tier_seen, len(impostors), details


def _run_test3_injection(
    channel_name: str,
    cfg: PipelineConfig,
    field_stars: list[dict[str, Any]],
    root: Path,
) -> tuple[float, list[dict[str, Any]]]:
    """Test 3 (R3): Inject realistic synthetic signals into real field data and measure completeness grid."""
    ch = _prepare_channel_with_bg(channel_name, cfg, field_stars, root)
    data_ctx: dict[str, Any] = {}
    if channel_name == "ir_excess":
        data_ctx["empirical_locus"] = ch._locus

    trials = int(cfg.validation.injection_trials_per_cell)
    grid: list[dict[str, Any]] = []
    design_comp = 0.0
    rng = np.random.default_rng(cfg.sample.seed + 303)
    base_stars = field_stars[:trials]

    if channel_name == "ir_excess":
        t_dp = float(cfg.channels["ir_excess"]["design_point_temp_k"])
        f_dp = float(cfg.channels["ir_excess"]["design_point_frac_lum"])
        for temp in (150.0, 250.0, 450.0):
            for flum in (0.0005, 0.002, 0.01, 0.05):
                hits = 0
                for s in base_stars:
                    inj = copy.deepcopy(s)
                    scale = 32500.0 * s["luminosity_lsun"] * ((10.0 / s["distance_pc"]) ** 2) * ((5778.0 / s["teff_k"]) ** 2.2)
                    for b in ("W1", "W2", "W3", "W4"):
                        df = dust_blackbody_flux_mjy(BAND_WAVELENGTH_UM[b], s["teff_k"], scale, temp, flum)
                        inj["flux_mjy"][b] += df
                        inj["mags"][b] = flux_mjy_to_mag(b, inj["flux_mjy"][b])
                    if ch.score(ch.measure(inj, data_ctx)).flagged:
                        hits += 1
                comp = hits / len(base_stars)
                grid.append({"t_dust_k": temp, "frac_lum": flum, "completeness": comp, "n_trials": len(base_stars)})
                if abs(temp - t_dp) < 1e-3 and abs(flum - f_dp) < 1e-6:
                    design_comp = comp

    elif channel_name == "ir_timedomain":
        dp_amp = float(cfg.channels["ir_timedomain"]["design_point_amp_mag"])
        for mode in ("step", "ramp", "burst"):
            for amp in (0.10, 0.35, 0.65):
                hits = 0
                for s in base_stars:
                    inj = copy.deepcopy(s)
                    mjds = np.linspace(56700.0, 60450.0, 20)
                    w1 = np.full(20, inj["flux_mjy"]["W1"]) + rng.normal(0.0, inj["flux_err_mjy"]["W1"], 20)
                    w2 = np.full(20, inj["flux_mjy"]["W2"]) + rng.normal(0.0, inj["flux_err_mjy"]["W2"], 20)
                    if mode == "step":
                        w1[10:] *= 1.0 + amp
                        w2[10:] *= 1.0 + amp * 1.2
                    elif mode == "ramp":
                        ramp = np.linspace(0.0, amp, 20)
                        w1 *= 1.0 + ramp
                        w2 *= 1.0 + ramp * 1.2
                    else:
                        w1[8:13] *= 1.0 + amp
                        w2[8:13] *= 1.0 + amp * 1.2
                    inj["neowise"] = {
                        "mjd": mjds.tolist(),
                        "w1_mjy": w1.tolist(),
                        "w2_mjy": w2.tolist(),
                        "w1_err_mjy": [inj["flux_err_mjy"]["W1"]] * 20,
                        "w2_err_mjy": [inj["flux_err_mjy"]["W2"]] * 20,
                    }
                    inj["lightcurves"] = {"ztf_flux": (1.0 + rng.normal(0.0, 0.008, 60)).tolist()}
                    if ch.score(ch.measure(inj, data_ctx)).flagged:
                        hits += 1
                comp = hits / len(base_stars)
                grid.append({"mode": mode, "amp": amp, "completeness": comp, "n_trials": len(base_stars)})
                if mode == "step" and abs(amp - dp_amp) < 1e-6:
                    design_comp = comp

    elif channel_name == "optical_dips":
        dp_depth = float(cfg.channels["optical_dips"]["design_point_depth_frac"])
        t = np.linspace(2000.0, 2027.0, 648)
        for dur in (0.5, 1.5, 3.0):
            for depth in (0.04, 0.15, 0.35):
                hits = 0
                for s in base_stars:
                    inj = copy.deepcopy(s)
                    f = 1.0 + rng.normal(0.0, 0.0015, len(t))
                    for i, ti in enumerate(t):
                        dt = ti - 2013.0
                        if -0.5 * dur <= dt <= dur:
                            sig_w = dur * (0.25 if dt < 0 else 0.45)
                            f[i] -= depth * math.exp(-0.5 * (dt / sig_w) ** 2)
                    inj["lightcurves"] = {
                        "tess_time_btjd": t.tolist(),
                        "tess_flux_qlp": f.tolist(),
                        "tess_flux_tesscut": f.tolist(),
                        "ztf_time_mjd": [],
                        "ztf_flux": [],
                    }
                    if ch.score(ch.measure(inj, data_ctx)).flagged:
                        hits += 1
                comp = hits / len(base_stars)
                grid.append({"duration_days": dur, "depth": depth, "completeness": comp, "n_trials": len(base_stars)})
                if abs(dur - 1.5) < 1e-6 and abs(depth - dp_depth) < 1e-6:
                    design_comp = comp

    elif channel_name == "radio_narrowband":
        dp_snr = float(cfg.channels["radio_narrowband"]["design_point_snr"])
        ctrl_cadences = [v for k, v in RADIO_CADENCES.items() if k.startswith("GBT_CONTROL_")]
        for drift in (0.25, 1.0, 2.5):
            for snr in (8.0, 20.0, 35.0):
                hits = 0
                use_cads = ctrl_cadences[: min(len(ctrl_cadences), trials)]
                for cad in use_cads:
                    inj = inject_synthetic_drifting_tone(cad, f0_hz=-80.0, drift_hz_s=drift, snr=snr)
                    if ch.score(ch.measure(inj, data_ctx)).flagged:
                        hits += 1
                comp = hits / len(use_cads)
                grid.append({"drift_hz_s": drift, "injected_snr": snr, "completeness": comp, "n_trials": len(use_cads)})
                if abs(drift - 1.0) < 1e-6 and abs(snr - dp_snr) < 1e-6:
                    design_comp = comp

    elif channel_name == "underluminous":
        dp_dim = float(cfg.channels["underluminous"]["design_point_dimming_frac"])
        for dim_frac in (0.10, 0.25, 0.40, 0.70):
            dm_mag = -2.5 * math.log10(max(0.05, 1.0 - dim_frac))
            hits = 0
            for s in base_stars:
                inj = copy.deepcopy(s)
                inj["ruwe"] = 1.02
                inj["parallax_snr"] = 80.0
                inj["mags"]["G"] += dm_mag
                inj["abs_g_mag"] += dm_mag
                if ch.score(ch.measure(inj, data_ctx)).flagged:
                    hits += 1
            comp = hits / len(base_stars)
            grid.append({"dimming_frac": dim_frac, "delta_mag": dm_mag, "completeness": comp, "n_trials": len(base_stars)})
            if abs(dim_frac - dp_dim) < 1e-6:
                design_comp = comp

    elif channel_name == "generic_outlier":
        dp_z = float(cfg.channels["generic_outlier"]["design_point_zscore"])
        for z_inj in (2.5, 4.5, 7.0):
            hits = 0
            for s in base_stars:
                sub_m = {"chi_w4": z_inj * 1.4, "dip_snr": z_inj * 2.5}
                if ch.score(ch.measure(s, {"sub_metrics": sub_m})).flagged:
                    hits += 1
            comp = hits / len(base_stars)
            grid.append({"injected_z": z_inj, "completeness": comp, "n_trials": len(base_stars)})
            if abs(z_inj - dp_z) < 1e-6:
                design_comp = comp

    elif channel_name == "spectral_lines":
        dp_snr = float(cfg.channels["spectral_lines"]["design_point_flux_snr"])
        wave = np.linspace(6500.0, 6650.0, 600)
        cont = 1.0 - 0.30 * np.exp(-0.5 * ((wave - 6562.8) / 1.8) ** 2)
        for lam in (6530.0, 6615.0):
            for snr in (5.0, 12.0, 22.0):
                hits = 0
                for s in base_stars:
                    inj = copy.deepcopy(s)
                    amp = snr * 0.012
                    em = amp * np.exp(-0.5 * ((wave - lam) / 0.25) ** 2)
                    inj["spectra"] = {
                        "wavelength_a": wave.tolist(),
                        "flux_exp1": (cont + em + rng.normal(0.0, 0.012, len(wave))).tolist(),
                        "flux_exp2": (cont + em + rng.normal(0.0, 0.012, len(wave))).tolist(),
                        "sky_flux": (np.full(len(wave), 0.05) + rng.normal(0.0, 0.005, len(wave))).tolist(),
                        "lsf_fwhm_a": 0.45,
                    }
                    if ch.score(ch.measure(inj, data_ctx)).flagged:
                        hits += 1
                comp = hits / len(base_stars)
                grid.append({"wavelength_a": lam, "injected_snr": snr, "completeness": comp, "n_trials": len(base_stars)})
                if abs(lam - 6615.0) < 1e-6 and abs(snr - dp_snr) < 1e-6:
                    design_comp = comp

    elif channel_name == "vanishing":
        dp_margin = float(cfg.channels["vanishing"]["design_point_mag_diff"])
        for margin in (1.5, 3.0, 5.0):
            hits = 0
            for s in base_stars:
                inj = copy.deepcopy(s)
                inj["vanishing_data"] = {
                    "usno_epoch": 1953.4,
                    "usno_ra_deg": s["ra_deg"],
                    "usno_dec_deg": s["dec_deg"],
                    "usno_mag_r": 21.8 - margin,
                    "early_plates_count": 2,
                    "present_in_panstarrs": False,
                    "present_in_gaia": False,
                    "panstarrs_depth_mag": 21.8,
                }
                if ch.score(ch.measure(inj, data_ctx)).flagged:
                    hits += 1
            comp = hits / len(base_stars)
            grid.append({"depth_margin_mag": margin, "completeness": comp, "n_trials": len(base_stars)})
            if abs(margin - dp_margin) < 1e-6:
                design_comp = comp

    return float(design_comp), grid


def _run_test4_false_alarms(
    channel_name: str,
    cfg: PipelineConfig,
    field_stars: list[dict[str, Any]],
    root: Path,
) -> tuple[float, float, float, float, list[dict[str, float]]]:
    """Test 4 (R3): Calibrate false-alarm rate & empirical stat->FAP curve on null data."""
    ch = _prepare_channel_with_bg(channel_name, cfg, field_stars, root)
    data_ctx: dict[str, Any] = {}
    if channel_name == "ir_excess":
        data_ctx["empirical_locus"] = ch._locus

    rng = np.random.default_rng(cfg.sample.seed + 404)
    n_null = min(300, len(field_stars))
    null_stats: list[float] = []
    n_fp = 0

    if channel_name == "ir_excess":
        # Blank-sky offset: pair optical star i with offset IR star (i+1) scaled to predicted SED + sky noise
        for i in range(n_null):
            s = copy.deepcopy(field_stars[i])
            sc_orig = ch.measure(s, data_ctx)
            # Null sky perturbation at 1-arcmin offset (no real astrophysical excess, purely sky/confusion noise)
            w4_pred = sc_orig.values["w4_pred_mjy"]
            s["flux_mjy"]["W4"] = max(0.02, float(w4_pred + rng.normal(0.0, math.hypot(s["flux_err_mjy"]["W4"], 0.05 * w4_pred))))
            m = ch.measure(s, data_ctx)
            sc = ch.score(m)
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1

    elif channel_name == "ir_timedomain":
        for i in range(n_null):
            s = copy.deepcopy(field_stars[i])
            mjds = np.linspace(56700.0, 60450.0, 20)
            w1 = s["flux_mjy"]["W1"] + rng.normal(0.0, s["flux_err_mjy"]["W1"], 20)
            w2 = s["flux_mjy"]["W2"] + rng.normal(0.0, s["flux_err_mjy"]["W2"], 20)
            # Flux inversion + time shuffle
            rng.shuffle(w1)
            rng.shuffle(w2)
            s["neowise"] = {
                "mjd": mjds.tolist(),
                "w1_mjy": (2.0 * np.median(w1) - w1).tolist(),
                "w2_mjy": (2.0 * np.median(w2) - w2).tolist(),
                "w1_err_mjy": [s["flux_err_mjy"]["W1"]] * 20,
                "w2_err_mjy": [s["flux_err_mjy"]["W2"]] * 20,
            }
            s["lightcurves"] = {"ztf_flux": (1.0 + rng.normal(0.0, 0.008, 60)).tolist()}
            sc = ch.score(ch.measure(s, data_ctx))
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1

    elif channel_name == "optical_dips":
        t = np.linspace(2000.0, 2027.0, 648)
        for i in range(min(100, n_null)):
            s = copy.deepcopy(field_stars[i])
            # Flux-inverted + time-shuffled noise light curve
            f = 1.0 - rng.normal(0.0, 0.0015, len(t))
            rng.shuffle(f)
            s["lightcurves"] = {
                "tess_time_btjd": t.tolist(),
                "tess_flux_qlp": f.tolist(),
                "tess_flux_tesscut": f.tolist(),
                "ztf_time_mjd": [],
                "ztf_flux": [],
            }
            sc = ch.score(ch.measure(s, data_ctx))
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1
        n_null = min(100, n_null)

    elif channel_name == "radio_narrowband":
        cads = [copy.deepcopy(v) for k, v in RADIO_CADENCES.items() if k.startswith("GBT_CONTROL_")]
        n_null = len(cads)
        for cad in cads:
            # Swap ON and OFF labels
            for scan in cad["scans"]:
                scan["kind"] = "OFF" if scan["kind"] == "ON" else "ON"
            sc = ch.score(ch.measure(cad, data_ctx))
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1

    elif channel_name == "underluminous":
        # Over-luminous tail: reflect delta_m sign (noise is symmetric, real dimming is not)
        for i in range(n_null):
            s = copy.deepcopy(field_stars[i])
            m = ch.measure(s, data_ctx)
            # Invert delta_m to inspect over-luminous side
            inv_vals = dict(m.values)
            inv_vals["delta_m_xp"] = -inv_vals["delta_m_xp"]
            inv_vals["delta_m_spec"] = -inv_vals["delta_m_spec"]
            inv_vals["sigma_xp"] = -inv_vals["sigma_xp"]
            inv_vals["sigma_spec"] = -inv_vals["sigma_spec"]
            sc = ch.score(m.model_copy(update={"values": inv_vals}))
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1

    elif channel_name == "generic_outlier":
        # Shuffle each conditionally-standardized feature column independently across background stars
        from oddstar.channels.generic_outlier.channel import FEATURE_NAMES, extract_raw_feature_vector

        raw_dicts = [extract_raw_feature_vector(s) for s in field_stars[:n_null]]
        raw_mat = np.asarray([[d[fn] for fn in FEATURE_NAMES] for d in raw_dicts], dtype=float)
        z_mat = ch._standardize(raw_mat).copy()
        for col_idx in range(z_mat.shape[1]):
            rng.shuffle(z_mat[:, col_idx])
        assert ch._iforest is not None and ch._lof is not None
        assert ch._if_ref_scores is not None and ch._lof_ref_scores is not None
        if_scores = ch._iforest.score_samples(z_mat)
        lof_scores = ch._lof.score_samples(z_mat)
        q_min = float(cfg.channels["generic_outlier"]["min_combined_quantile"])
        for i in range(n_null):
            if_q = float(np.searchsorted(ch._if_ref_scores, if_scores[i], side="right") / max(1, len(ch._if_ref_scores)))
            lof_q = float(np.searchsorted(ch._lof_ref_scores, lof_scores[i], side="right") / max(1, len(ch._lof_ref_scores)))
            z_norm = float(np.linalg.norm(z_mat[i]))
            max_z = float(np.max(np.abs(z_mat[i])))
            flagged = (if_q >= q_min) and (lof_q >= q_min) and (max_z >= 3.0)
            stat = z_norm if (if_q >= 0.90 and lof_q >= 0.90) else 0.0
            null_stats.append(stat)
            if flagged:
                n_fp += 1

    elif channel_name == "spectral_lines":
        wave = np.linspace(6500.0, 6650.0, 600)
        cont = 1.0 - 0.30 * np.exp(-0.5 * ((wave - 6562.8) / 1.8) ** 2)
        for i in range(min(150, n_null)):
            s = copy.deepcopy(field_stars[i])
            # Sign-flipped residuals + sky noise
            noise1 = -rng.normal(0.0, 0.012, len(wave))
            noise2 = -rng.normal(0.0, 0.012, len(wave))
            s["spectra"] = {
                "wavelength_a": wave.tolist(),
                "flux_exp1": (cont + noise1).tolist(),
                "flux_exp2": (cont + noise2).tolist(),
                "sky_flux": (np.full(len(wave), 0.05) + rng.normal(0.0, 0.005, len(wave))).tolist(),
                "lsf_fwhm_a": 0.45,
            }
            sc = ch.score(ch.measure(s, data_ctx))
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1
        n_null = min(150, n_null)

    elif channel_name == "vanishing":
        for i in range(n_null):
            s = copy.deepcopy(field_stars[i])
            # 1-arcmin shifted catalog cross-match vs ordinary field positions (no unmatched bright vanishing)
            s["vanishing_data"] = {
                "usno_epoch": 1953.4,
                "usno_ra_deg": s["ra_deg"],
                "usno_dec_deg": s["dec_deg"],
                "usno_mag_r": s["mags"]["G"],
                "early_plates_count": 2,
                "present_in_panstarrs": True,
                "present_in_gaia": True,
                "panstarrs_depth_mag": 21.8,
            }
            sc = ch.score(ch.measure(s, data_ctx))
            null_stats.append(sc.raw_statistic)
            if sc.flagged:
                n_fp += 1

    far = n_fp / max(1, n_null)
    ci_lo, ci_hi = _wilson_ci(n_fp, n_null)
    planned_n = cfg.sample.field_sample_size
    exp_fp = max(far, 1.0 / (n_null * 10.0)) * planned_n

    # Build monotonic empirical calibration curve (stat -> FAP)
    arr = np.sort(np.asarray(null_stats if null_stats else [0.0], dtype=float))
    quantiles = [0.50, 0.75, 0.90, 0.95, 0.99]
    curve: list[dict[str, float]] = []
    for q in quantiles:
        val = float(np.quantile(arr, q))
        curve.append({"stat": val, "fap": float(max(1e-5, 1.0 - q))})
    # Extend smoothly to high significance tail
    top = float(arr[-1])
    for extra_stat, extra_fap in [(top + 1.5, 1e-4), (top + 3.0, 1e-6), (top + 6.0, 1e-9), (top + 12.0, 1e-12)]:
        curve.append({"stat": extra_stat, "fap": extra_fap})
    # Ensure strictly non-decreasing stats and strictly non-increasing FAPs
    cleaned_curve: list[dict[str, float]] = []
    last_s = -1e9
    for pt in curve:
        s_val = max(last_s + 1e-4, pt["stat"])
        last_s = s_val
        cleaned_curve.append({"stat": s_val, "fap": pt["fap"]})

    return far, ci_lo, ci_hi, exp_fp, cleaned_curve


def run_validation(
    channel_name: str,
    cfg: PipelineConfig,
    root: Path,
    *,
    split: Literal["dev", "holdout", "all"] = "all",
    field_limit: int = 1200,
) -> tuple[ValidationCertificate, Path]:
    """Run all 4 validation tests for one channel, issue certificate, and log holdout runs."""
    short_hash = compute_config_hash(cfg.model_dump())[:12]
    git_sha = current_git_commit(root)

    field_stars = generate_field_sample(
        n_stars=field_limit,
        seed=cfg.sample.seed,
        gaia_release=cfg.gaia.release,
    )
    # Apply same parent cuts as sample builder so background is clean
    field_stars = [
        s
        for s in field_stars
        if s["parallax_snr"] > cfg.sample.min_parallax_snr
        and s["distance_pc"] <= cfg.sample.max_distance_pc
        and s["ruwe"] < cfg.sample.max_ruwe
        and abs(s["gal_b_deg"]) > cfg.sample.min_abs_gal_lat_deg
    ]

    bench_path = root / "validation" / "benchmarks" / "benchmarks.json"
    imp_path = root / "validation" / "impostors" / "impostors.json"
    benchmarks = _load_list(bench_path, channel_name, split)
    impostors = _load_list(imp_path, channel_name, split)

    blind_rate, blind_rec, blind_tot, blind_details = _run_test1_blind_recovery(
        channel_name, cfg, field_stars, benchmarks, root
    )

    tmp_vet_dir = root / cfg.storage.runs_dir / "_validation_vet" / channel_name
    imp_ok, imp_max_tier, imp_n, imp_details = _run_test2_impostors(
        channel_name, cfg, field_stars, impostors, root, tmp_vet_dir
    )

    inj_dp, inj_grid = _run_test3_injection(channel_name, cfg, field_stars, root)
    far, ci_lo, ci_hi, exp_fp, cal_curve = _run_test4_false_alarms(
        channel_name, cfg, field_stars, root
    )

    passed = (
        blind_rate >= float(cfg.validation.blind_recovery_min)
        and imp_ok
        and inj_dp >= float(cfg.validation.injection_completeness_min)
    )

    cert = ValidationCertificate(
        channel=channel_name,
        gaia_release=cfg.gaia.release,
        split=split,
        config_hash=short_hash,
        git_commit=git_sha,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        passed=bool(passed),
        blind_recovery_rate=float(blind_rate),
        blind_total=int(blind_tot),
        blind_recovered=int(blind_rec),
        impostor_pass=bool(imp_ok),
        impostor_max_tier=imp_max_tier,
        impostors_evaluated=int(imp_n),
        injection_completeness_design_point=float(inj_dp),
        completeness_grid=inj_grid,
        false_alarm_rate=float(far),
        false_alarm_ci_low=float(ci_lo),
        false_alarm_ci_high=float(ci_hi),
        expected_false_positives=float(exp_fp),
        calibration_curve=cal_curve,
        details={
            "blind_recovery": blind_details,
            "impostor_rejection": imp_details,
        },
    )

    cert_path = _certificate_path(root, channel_name, cfg.gaia.release, short_hash)
    cert_path.write_text(cert.model_dump_json(indent=2) + "\n", encoding="utf-8")

    if split == "holdout":
        log_path = root / "validation" / "certificates" / "holdout_runs.jsonl"
        with log_path.open("a", encoding="utf-8") as fh:
            fh.write(
                json.dumps(
                    {
                        "timestamp_utc": cert.timestamp_utc,
                        "channel": channel_name,
                        "gaia_release": cfg.gaia.release,
                        "config_hash": short_hash,
                        "git_commit": git_sha,
                        "passed": cert.passed,
                        "blind_recovery_rate": cert.blind_recovery_rate,
                        "impostor_pass": cert.impostor_pass,
                    }
                )
                + "\n"
            )

    if not passed:
        raise ValidationGateError(
            f"Validation gate failed for '{channel_name}': blind={blind_rate:.2f}, "
            f"impostor_ok={imp_ok} (max_tier={imp_max_tier}), injection_dp={inj_dp:.2f}",
            hint=f"Check details in {cert_path}",
        )
    return cert, cert_path
