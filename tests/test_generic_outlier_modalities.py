"""generic_outlier: absent modalities carry no evidence (L003) and present ones are judged against their own noise model."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oddstar.channels import build_channel
from oddstar.channels.generic_outlier.channel import (
    MODALITY_FEATURES,
    _IsolationForest,
    extract_modality_presence,
    extract_raw_feature_vector,
    modality_significance,
)
from oddstar.config import load_config, repo_root
from oddstar.io import generate_field_sample, get_named_object
from oddstar.validate import run_validation
from oddstar.vet import vet_target


@pytest.fixture(scope="module")
def cfg() -> Any:
    c, _, _ = load_config()
    return c


@pytest.fixture(scope="module")
def background(cfg: Any) -> list[dict[str, Any]]:
    stars = generate_field_sample(n_stars=900, seed=cfg.sample.seed, gaia_release=cfg.gaia.release)
    return [s for s in stars if s["parallax_snr"] > cfg.sample.min_parallax_snr and s["ruwe"] < cfg.sample.max_ruwe
            and abs(s["gal_b_deg"]) > cfg.sample.min_abs_gal_lat_deg and s["distance_pc"] <= cfg.sample.max_distance_pc]


@pytest.fixture(scope="module")
def channel(cfg: Any, background: list[dict[str, Any]]) -> Any:
    ch = build_channel("generic_outlier", cfg, repo_root())
    ch.fit_background(background)
    return ch


def _with_noise_modalities(star: dict[str, Any], rng: np.random.Generator) -> dict[str, Any]:
    s = copy.deepcopy(star)
    s["lightcurves"] = {"tess_flux_qlp": (1.0 + rng.normal(0.0, 0.002, 1300)).tolist()}
    w2, e2 = s["flux_mjy"]["W2"], s["flux_err_mjy"]["W2"]
    s["neowise"] = {"w2_mjy": (w2 + rng.normal(0.0, e2, 20)).tolist(), "w2_err_mjy": [e2] * 20}
    s["spectra"] = {"flux_exp1": (1.0 + rng.normal(0.0, 0.01, 600)).tolist(), "flux_exp2": (1.0 + rng.normal(0.0, 0.01, 600)).tolist()}
    s["vanishing_data"] = {"usno_mag_r": 15.0, "present_in_panstarrs": True, "present_in_gaia": True}
    return s


def test_absent_modalities_are_exactly_zero_and_marked_absent(cfg: Any, background: list[dict[str, Any]], channel: Any) -> None:
    star = background[0]
    presence = extract_modality_presence(star)
    assert all(presence[fn] is None for fn in MODALITY_FEATURES)
    sig = modality_significance(extract_raw_feature_vector(star), presence, cfg.channels["generic_outlier"])
    assert all(sig[fn] == 0.0 for fn in MODALITY_FEATURES)
    m = channel.measure(star, {})
    assert all(m.values[f"z_{fn}"] == 0.0 and m.values[f"present_{fn}"] == 0.0 for fn in MODALITY_FEATURES)
    assert any("Modalities not observed" in c for c in m.caveats)


def test_ordinary_present_modalities_are_not_outliers(background: list[dict[str, Any]], channel: Any) -> None:
    """Mirror of L003: having a light curve, spectra or NEOWISE visits with only noise must not make a star an outlier
    against a background that lacks those observations (before the fix every such star had |z| ~ 10-15)."""
    rng = np.random.default_rng(11)
    base = background[:200]
    noisy = [_with_noise_modalities(s, rng) for s in base]
    z = np.array([[channel.measure(s, {}).values[f"z_{fn}"] for fn in ("ir_var_chi2", "dip_snr", "spec_line_snr")] for s in noisy])
    assert np.all(np.abs(z) < 2.5)
    n_bare = sum(channel.score(channel.measure(s, {})).flagged for s in base)
    n_noisy = sum(channel.score(channel.measure(s, {})).flagged for s in noisy)
    assert n_noisy <= n_bare + 1


def test_forest_isolates_a_query_on_a_feature_constant_in_the_background() -> None:
    """Liu et al. draw the split attribute from all attributes; a constant one still isolates an out-of-hull query."""
    rng = np.random.default_rng(3)
    x = np.column_stack([rng.normal(0.0, 1.0, 400), rng.normal(0.0, 1.0, 400), np.zeros(400)])
    forest = _IsolationForest(n_trees=200, subsample_size=256, seed=5)
    forest.fit(x)
    ref = np.sort(forest.score_samples(x))
    s_out, s_in = forest.score_samples(np.array([[0.0, 0.0, 40.0], [0.0, 0.0, 0.5]]))
    assert np.searchsorted(ref, s_out) / len(ref) >= 0.99
    assert np.searchsorted(ref, s_in) / len(ref) < 0.9


def test_single_modality_benchmarks_are_flagged_by_both_detectors(cfg: Any, channel: Any) -> None:
    for sid in ("6049824505359595776", "755500000000000003"):  # dip-only and vanishing-only benchmarks
        obj = get_named_object(sid, gaia_release=cfg.gaia.release)
        assert obj is not None
        m = channel.measure(obj, {})
        assert m.values["iforest_quantile"] >= 0.99 and m.values["lof_quantile"] >= 0.99
        assert channel.score(m).flagged


def test_single_exposure_spectral_spike_is_a_lead_that_vetting_rejects(cfg: Any, channel: Any, tmp_path: Path) -> None:
    imp = get_named_object("744400000000000003", gaia_release=cfg.gaia.release)  # SL-CosmicRay-1
    assert imp is not None
    sc = channel.score(channel.measure(imp, {}))
    assert sc.feature_contributions["spec_line_snr"] > 10.0
    vet = vet_target(imp, [sc], cfg, tmp_path / "ev")
    assert {v.check: v.status for v in vet}["cosmic_ray"] == "fail"


def test_generic_outlier_passes_all_validation_gates(cfg: Any, tmp_path: Path) -> None:
    root = repo_root()
    for sub in ("benchmarks", "impostors"):
        (tmp_path / "validation" / sub).mkdir(parents=True)
        (tmp_path / "validation" / sub / f"{sub}.json").write_text(
            (root / "validation" / sub / f"{sub}.json").read_text(encoding="utf-8"), encoding="utf-8")
    cert, _ = run_validation("generic_outlier", cfg, tmp_path, split="all", field_limit=1200)
    assert cert.passed
    assert cert.blind_recovered == cert.blind_total == 6
    assert cert.impostor_pass
    assert cert.injection_completeness_design_point >= cfg.validation.injection_completeness_min
    assert cert.code_fingerprint and len(cert.code_fingerprint) == 64
