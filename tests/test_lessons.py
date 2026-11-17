"""Regression tests linked directly to lessons/ entries (Ground Rule R7)."""

from __future__ import annotations

from pathlib import Path

from oddstar.channels import build_channel
from oddstar.config import load_config
from oddstar.fuse import assign_tier, build_hypothesis_table
from oddstar.io import get_named_object
from oddstar.vet import vet_target


def test_lesson_l001_single_epoch_w4_background_galaxy_cap(tmp_path: Path) -> None:
    """L001: Passing radio_counterpart must NOT clear background_galaxy or promote above Tier C
    when deep_imaging or second_instrument is not_evaluated.
    """
    cfg, _, _ = load_config()
    star = get_named_object("4582145718700356096")
    assert star is not None

    ch = build_channel("ir_excess", cfg, tmp_path)
    meas = ch.measure(star, {})
    score = ch.score(meas)
    assert score.flagged is True

    vet = vet_target(star, [score], cfg, tmp_path / "ev")
    by_check = {v.check: v for v in vet}
    assert by_check["radio_counterpart"].status == "pass"
    assert by_check["deep_imaging"].status == "not_evaluated"
    assert by_check["second_instrument"].status == "not_evaluated"

    tier, _ = assign_tier([score], vet, cfg)
    assert tier == "C"

    hypotheses = {h.name: h for h in build_hypothesis_table([score], vet)}
    assert hypotheses["background_galaxy"].status == "untested"
    assert "radio_counterpart passed" in hypotheses["background_galaxy"].note
    assert hypotheses["unexplained"].status == "untested"


def test_lesson_l002_missing_plate_data_no_spurious_vanishing_on_high_pm_stars(tmp_path: Path) -> None:
    """L002: Missing vanishing_data must NOT default USNO coords to J2016 Gaia coords
    or emit a spurious proper-motion vanishing trigger on high-PM live stars.
    """
    cfg, _, _ = load_config()
    star = {
        "source_id": "5628404318169472",
        "ra_deg": 43.54827286,
        "dec_deg": 5.88271224,
        "pmra_mas_yr": 82.5,
        "pmdec_mas_yr": -63.8,
        "mags": {"G": 13.43},
        "vanishing_data": {},
    }
    ch = build_channel("vanishing", cfg, tmp_path)
    meas = ch.measure(star, {})
    score = ch.score(meas)

    assert meas.values["early_plates_count"] == 0.0
    assert meas.values["sep_arcsec"] == 0.0
    assert meas.values["modern_matched"] == 1.0
    assert score.flagged is False
    assert any("No historical plate catalog" in c for c in meas.caveats)


def test_lesson_l003_missing_modality_routing_and_rv_binary_honesty(tmp_path: Path) -> None:
    """L003: Missing modalities in generic_outlier must stay 0.0 (no random jitter triggering
    phantom light-curve holes), metal-poor outliers must run metallicity/evolved checks -> Explained,
    and faint-star Gaia RVS (rv_faint_unconstrained=True) must mark binary=not_evaluated rather than pass.
    """
    from oddstar.channels.generic_outlier.channel import extract_raw_feature_vector
    from oddstar.io.physics import propagate_position
    from oddstar.schemas import ChannelScore
    from oddstar.vet.ir_checks import check_binary

    cfg, _, _ = load_config()
    wise_ra, wise_dec = propagate_position(115.8556, 72.8083, 150.0, -1200.0, 2016.0, 2010.2)
    subdwarf = {
        "source_id": "1111149103659931264",
        "ra_deg": 115.8556,
        "dec_deg": 72.8083,
        "pmra_mas_yr": 150.0,
        "pmdec_mas_yr": -1200.0,
        "abs_g_mag": 8.47,
        "teff_k": 4100.0,
        "feh": -1.09,
        "ruwe": 1.02,
        "mags": {"G": 12.53, "BP": 13.32, "RP": 11.66, "J": 10.58, "Ks": 9.85, "W3": 9.66, "W4": 8.38},
        "flux_mjy": {"Ks": 70.0, "W4": 3.7},
        "flux_err_mjy": {"W4": 1.0},
        "allwise": {"ra_deg": wise_ra, "dec_deg": wise_dec, "n_neighbours_gaia_xmatch": 1},
        "neowise": {},
        "lightcurves": {},
        "spectra": {},
        "vanishing_data": {},
        "vetting_archive": {"nss_flag": 0, "rv_scatter_kms": 6.12, "rv_faint_unconstrained": True},
    }
    raw = extract_raw_feature_vector(subdwarf)
    assert raw["dip_snr"] == 0.0
    assert raw["ir_var_chi2"] == 0.0
    assert raw["spec_line_snr"] == 0.0
    assert raw["vanished_margin"] == 0.0
    assert raw["radio_snr"] == 0.0

    # Even if an outlier score had high bp_rp / dip contributions, missing lightcurves do not trigger lc_checks,
    # and metallicity check explains metal-poor subdwarfs.
    score = ChannelScore(
        source_id="1111149103659931264",
        channel="generic_outlier",
        raw_statistic=7.35,
        false_alarm_prob=7.1e-5,
        flagged=True,
        rank_score=4.15,
        summary="Outlier",
        feature_contributions={"bp_rp": 3.52, "dip_snr": 3.18},
    )
    vet = vet_target(subdwarf, [score], cfg, tmp_path / "ev_sub")
    by_check = {v.check: v for v in vet}
    assert by_check["metallicity"].status == "fail"
    assert "pipeline_agreement" not in by_check
    assert "second_survey" not in by_check
    tier, reason = assign_tier([score], vet, cfg)
    assert tier == "Explained"
    assert reason and "metallicity" in reason

    bin_res = check_binary(subdwarf, cfg, tmp_path / "ev_bin")
    assert bin_res.status == "not_evaluated"


def test_lesson_l004_known_pulsating_variable_simbad_types(tmp_path: Path) -> None:
    """L004: Classical pulsating and evolved SIMBAD variable types (e.g. RR*, LP*, Cep, dS*)
    must fail check_known_variable and resolve generic_outlier targets as Explained rather than Tier B.
    """
    from oddstar.schemas import ChannelScore
    from oddstar.vet.lc_checks import check_known_variable

    cfg, _, _ = load_config()
    rr_star = {
        "source_id": "6164530246002222208",
        "ra_deg": 204.4323,
        "dec_deg": -35.5543,
        "abs_g_mag": 0.503,
        "teff_k": 6439.0,
        "feh": -0.876,
        "ruwe": 1.389,
        "simbad_type": "RR*",
        "mags": {"G": 12.09, "BP": 12.33, "RP": 11.66, "J": 11.08, "Ks": 10.80, "W3": 10.74, "W4": 8.61},
        "allwise": {"ra_deg": 204.4323, "dec_deg": -35.5543, "n_neighbours_gaia_xmatch": 1},
        "lightcurves": {},
        "vetting_archive": {"nss_flag": 0, "rv_scatter_kms": 53.6},
    }
    kv_res = check_known_variable(rr_star, cfg, tmp_path / "ev_kv")
    assert kv_res.status == "fail"

    score = ChannelScore(
        source_id="6164530246002222208",
        channel="generic_outlier",
        raw_statistic=7.86,
        false_alarm_prob=3.77e-5,
        flagged=True,
        rank_score=1.71,
        summary="Outlier",
        feature_contributions={"delta_m_sigma": -7.49, "ruwe": 1.46},
    )
    vet = vet_target(rr_star, [score], cfg, tmp_path / "ev_rr")
    tier, reason = assign_tier([score], vet, cfg)
    assert tier == "Explained"
    assert reason and "known_variable" in reason
