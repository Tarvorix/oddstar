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


def test_lesson_l005_channel_scoped_metallicity_wise_routing_and_southern_ls_dr10(tmp_path: Path, monkeypatch: Any) -> None:
    """L005: (1) metallicity=fail must NOT explain away an active ir_excess flag in assign_tier;
    (2) WISE-driven generic_outlier flags (ks_w3, w3_w4) must route through run_ir_vetting (rejecting AABU/neighbours);
    (3) check_wise_quality catches hh00 W1/W2 halo contamination and low w4nm/w4m frame fractions;
    (4) southern Dec <= -30 targets query DECam ls_dr10.tractor via CSV TAP without a Dec > -30 gate.
    """
    import requests
    from oddstar.io.cache import RawCache
    from oddstar.io.live_survey import enrich_star_live_for_vetting
    from oddstar.schemas import ChannelScore

    cfg, _, _ = load_config()

    # 1. metallicity=fail must NOT mask an active ir_excess flag
    ir_star = {
        "source_id": "6739572937544816000",
        "ra_deg": 292.0619,
        "dec_deg": -36.5529,
        "abs_g_mag": 4.02,
        "teff_k": 5745.0,
        "feh": -1.44,
        "ruwe": 0.96,
        "mags": {"G": 13.93, "BP": 14.33, "RP": 13.36, "J": 12.66, "Ks": 12.20, "W1": 12.06, "W2": 12.05, "W3": 10.68, "W4": 8.05},
        "flux_mjy": {"W4": 5.0},
        "allwise": {"ra_deg": 292.0619, "dec_deg": -36.5529, "cc_flags": "0000", "ext_flg": 0, "ph_qual": "AAAB", "w4rchi2": 0.92, "w4m": 14, "n_neighbours_gaia_xmatch": 1},
        "wise_single_frames_w4_mjy": [5.0] * 14,
        "vetting_archive": {"deep_imaging_available": False, "second_instrument_ir": None, "neighbours_in_w4_beam": [], "radio_sources_within_10arcsec": [], "moving_objects_2010": [], "nss_flag": 0, "rv_scatter_kms": 0.2},
    }
    sc_ir = ChannelScore(source_id="6739572937544816000", channel="ir_excess", raw_statistic=7.3, false_alarm_prob=3.9e-7, flagged=True, rank_score=12.0, summary="IR")
    sc_out = ChannelScore(source_id="6739572937544816000", channel="generic_outlier", raw_statistic=9.1, false_alarm_prob=1e-12, flagged=True, rank_score=12.0, summary="Out", feature_contributions={"chi_w4": 3.9, "ks_w3": 2.5})
    vet_ir = vet_target(ir_star, [sc_ir, sc_out], cfg, tmp_path / "ev_ir_feh")
    assert {v.check: v.status for v in vet_ir}["metallicity"] == "fail"
    tier_ir, _ = assign_tier([sc_ir, sc_out], vet_ir, cfg)
    assert tier_ir == "C"

    # 2. WISE-driven generic_outlier (ks_w3=-4.83, AABU, 17 neighbours) routes through IR vetting -> Rejected
    lmc_star = {
        "source_id": "4651844963414920576",
        "ra_deg": 82.4412,
        "dec_deg": -70.9409,
        "abs_g_mag": 3.86,
        "teff_k": 5980.0,
        "feh": -0.86,
        "ruwe": 0.89,
        "mags": {"G": 12.07, "BP": 12.37, "RP": 11.58, "J": 11.04, "Ks": 10.64, "W1": 10.62, "W2": 10.66, "W3": 11.87, "W4": 9.69},
        "flux_mjy": {"W4": 1.1},
        "allwise": {"ra_deg": 82.4411, "dec_deg": -70.9407, "cc_flags": "0000", "ext_flg": 0, "ph_qual": "AABU", "w4rchi2": 1.07, "w4m": 46, "n_neighbours_gaia_xmatch": 1},
        "vetting_archive": {"deep_imaging_available": False, "neighbours_in_w4_beam": [{"id": "n1", "sep_arcsec": 5.0, "flux_frac": 0.30}], "radio_sources_within_10arcsec": [], "moving_objects_2010": []},
    }
    sc_lmc = ChannelScore(source_id="4651844963414920576", channel="generic_outlier", raw_statistic=5.29, false_alarm_prob=0.0078, flagged=True, rank_score=2.11, summary="Out", feature_contributions={"ks_w3": -4.83, "w3_w4": -1.78})
    vet_lmc = vet_target(lmc_star, [sc_lmc], cfg, tmp_path / "ev_lmc")
    by_lmc = {v.check: v for v in vet_lmc}
    assert by_lmc["wise_quality"].status == "fail"
    assert by_lmc["neighbours"].status == "fail"
    tier_lmc, _ = assign_tier([sc_lmc], vet_lmc, cfg)
    assert tier_lmc == "Rejected"

    # 3. hh00 cc_flags and low w4nm/w4m frame detection fraction fail wise_quality
    from oddstar.vet.ir_checks import check_wise_quality
    hh_star = dict(ir_star, allwise=dict(ir_star["allwise"], cc_flags="hh00"))
    assert check_wise_quality(hh_star, cfg, tmp_path / "ev_hh").status == "fail"
    low_nm_star = dict(ir_star, allwise=dict(ir_star["allwise"], ph_qual="AABB", w4m=41, w4nm=3))
    assert check_wise_quality(low_nm_star, cfg, tmp_path / "ev_nm").status == "fail"

    # 4. Southern Dec <= -30 queries NOIRLab DECam ls_dr10.tractor via CSV
    class _Resp:
        def __init__(self, status_code: int, text: str = "", json_data: Any = None) -> None:
            self.status_code = status_code
            self.text = text
            self._json = json_data or {"metadata": [], "data": []}
        def json(self) -> Any:
            return self._json

    seen_ls = []
    def fake_post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _Resp:
        q = (data or {}).get("QUERY", "")
        if "ls_dr10.tractor" in q:
            seen_ls.append((data or {}).get("FORMAT"))
            return _Resp(200, text="ls_id,ra,dec,type,flux_g,flux_r,flux_z,flux_w3,flux_w4,flux_ivar_w4\n1,292.061395,-36.551841,EXP,0,0,0,131.8,0,0\n")
        return _Resp(200)

    def fake_get(url: str, params: dict[str, Any] | None = None, timeout: int = 15) -> _Resp:
        return _Resp(200, text="# empty\n")

    monkeypatch.setattr(requests, "post", fake_post)
    monkeypatch.setattr(requests, "get", fake_get)
    cache = RawCache(tmp_path / "raw_ls", max_cache_gb=1.0)
    south_star = {
        "source_id": "6739572937544816000",
        "name": "Gaia DR3 6739572937544816000",
        "ra_deg": 292.061921,
        "dec_deg": -36.552873,
        "mags": {"G": 13.93, "W4": 8.05},
        "flux_mjy": {"W4": 5.0},
        "allwise": {"ra_deg": 292.0619, "dec_deg": -36.5528},
        "neowise": {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []},
        "lightcurves": {"ztf_flux": []},
        "vetting_archive": {"radio_sources_within_10arcsec": [], "deep_imaging_available": False, "deep_imaging_offset_galaxy": None},
        "wise_single_frames_w4_mjy": [],
        "live_enriched": False,
    }
    enriched, _ = enrich_star_live_for_vetting(south_star, cfg, cache)
    assert seen_ls == ["csv"]
    # The seeing-limited LS DR10 catalogue is a neighbour screen only (it can never clear deep_imaging). The
    # extended Tractor source ~4" away is recorded; with no photosphere (no Ks/Teff) the target excess is
    # unknown, so it is not promoted to a contaminant.
    va = enriched["vetting_archive"]
    assert va["query_status"]["deep_optical_ls_dr10"] == "ok"
    assert any("Legacy Survey DR10" in s for s in va["neighbour_screen"]["surveys"])
    ext = va["neighbour_screen"]["extended_neighbours"]
    assert len(ext) == 1 and ext[0]["ls_type"] == "EXP" and 3.0 < ext[0]["sep_arcsec"] < 5.0
    assert "deep_imaging_available" not in va and "deep_imaging_offset_galaxy" not in va
