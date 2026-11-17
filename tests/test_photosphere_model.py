"""Empirical-locus photosphere for ir_excess: sample/persisted/literature sources and the candidate regression."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oddstar.channels import build_channel
from oddstar.channels.ir_excess.channel import build_empirical_locus, literature_locus, predict_from_empirical_locus
from oddstar.config import load_config
from oddstar.dossier.writer import write_dossier
from oddstar.io import generate_field_sample, get_named_object
from oddstar.io.live_survey import _row_to_star_dict
from oddstar.io.physics import (
    BAND_WAVELENGTH_UM,
    dwarf_ks_minus_wise_sequence,
    load_dwarf_colour_table,
    planck_photosphere_flux_mjy,
    stellar_photosphere_flux_mjy,
)

CANDIDATE_ROW: dict[str, Any] = {
    "source_id": 5628404318169472,
    "ra": 43.54827,
    "dec": 5.88271,
    "b": -45.61,
    "parallax": 4.3149,
    "parallax_error": 0.0162,
    "parallax_over_error": 267.14,
    "pmra": 37.914,
    "pmdec": -97.175,
    "ruwe": 1.078,
    "phot_g_mean_mag": 13.430,
    "phot_bp_mean_mag": 14.041,
    "phot_rp_mean_mag": 12.687,
    "teff_gspphot": 4510.0,
    "logg_gspphot": 4.6,
    "mh_gspphot": 0.036,
    "ebpminrp_gspphot": 0.039,
    "j_m": 11.762,
    "j_msigcom": 0.022,
    "h_m": 11.143,
    "h_msigcom": 0.027,
    "ks_m": 10.971,
    "ks_msigcom": 0.019,
    "w1mpro": 10.665,
    "w1mpro_error": 0.023,
    "w2mpro": 10.276,
    "w2mpro_error": 0.021,
    "w3mpro": 8.701,
    "w3mpro_error": 0.026,
    "w4mpro": 7.290,
    "w4mpro_error": 0.102,
    "ph_qual": "AAAA",
    "tmass_ph_qual": "AAA",
}


def candidate_star() -> dict[str, Any]:
    star = _row_to_star_dict(dict(CANDIDATE_ROW))
    star["vetting_archive"]["dust"] = {"ebv_sfd_los": 0.2414, "ebv_sf11_los": None, "i100_mjy_sr": None}
    star["vetting_archive"]["query_status"] = {"irsa_dust": "ok"}
    return star


def test_literature_table_provenance_and_coverage() -> None:
    t = load_dwarf_colour_table()
    assert t["provenance"] in {"retrieved", "memory"}
    assert t["url"].endswith("EEM_dwarf_UBVIJHK_colors_Teff.txt")
    assert t["retrieved_utc_date"]
    teffs, bprp, vals = dwarf_ks_minus_wise_sequence("W1")
    assert len(teffs) > 40 and teffs == sorted(teffs)
    # W2..W4 stop at the first untabulated row (K5V), never interpolating across the gap
    for band in ("W2", "W3", "W4"):
        teffs_b, _, _ = dwarf_ks_minus_wise_sequence(band)
        assert min(teffs_b) == 4440.0


def test_photosphere_reproduces_dwarf_colours_and_planck_is_kept() -> None:
    # At a tabulated Teff the corrected photosphere has the tabulated Ks - Wx colour
    from oddstar.io.physics import BAND_ZERO_POINT_JY

    row = next(r for r in load_dwarf_colour_table()["rows"] if r["spt"] == "K2V")
    teff = float(row["teff_k"])
    f_ks = stellar_photosphere_flux_mjy(BAND_WAVELENGTH_UM["Ks"], teff, 100.0)
    for band, col in (("W1", row["ks_w1"]), ("W2", row["ks_w1"] + row["w1_w2"]), ("W3", row["ks_w1"] + row["w1_w3"]), ("W4", row["ks_w1"] + row["w1_w4"])):
        f_w = stellar_photosphere_flux_mjy(BAND_WAVELENGTH_UM[band], teff, 100.0)
        ks_w = 2.5 * np.log10((f_w / BAND_ZERO_POINT_JY[band]) / (f_ks / BAND_ZERO_POINT_JY["Ks"]))
        assert abs(ks_w - col) < 1e-6
    # Shortward of Ks the model is the Planck curve
    for band in ("G", "BP", "RP", "J", "H", "Ks"):
        lam = BAND_WAVELENGTH_UM[band]
        assert stellar_photosphere_flux_mjy(lam, teff, 100.0) == pytest.approx(planck_photosphere_flux_mjy(lam, teff, 100.0))


def test_sample_locus_clips_excess_stars_and_uses_locus_scatter() -> None:
    cfg, _, _ = load_config()
    ch_cfg = cfg.channels["ir_excess"]
    stars = generate_field_sample(1500, seed=cfg.sample.seed)
    base = build_empirical_locus(stars, ch_cfg)
    assert base is not None and base["source"] == "sample_locus" and base["n_stars"] > 1000
    # Contaminate 10% of stars with a large W3/W4 excess: the clipped median must barely move
    dirty = copy.deepcopy(stars)
    for s in dirty[::10]:
        s["mags"]["W3"] -= 1.5
        s["mags"]["W4"] -= 2.5
    contaminated = build_empirical_locus(dirty, ch_cfg)
    assert contaminated is not None
    for band in ("W3", "W4"):
        a = np.asarray(base["bands"][band]["median"])
        b = np.asarray(contaminated["bands"][band]["median"])
        assert np.max(np.abs(a - b)) < 0.03
        assert sum(contaminated["bands"][band]["n_clipped"]) > 100
    # Locus scatter enters the predicted uncertainty (it used to be multiplied by 0.0)
    pred = predict_from_empirical_locus(stars[0], base, ch_cfg)
    for band in ("W1", "W2", "W3", "W4"):
        assert pred[band]["intrinsic_scatter_mag"] > 0.02
        assert pred[band]["sigma_mag"] > pred[band]["intrinsic_scatter_mag"]


def test_locus_persisted_after_sample_run_and_reused_in_single_star_mode(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    ch = build_channel("ir_excess", cfg, tmp_path)
    ch.sample(cfg, limit=600)
    assert ch._locus["source"] == "sample_locus"
    path = ch.last_persisted_locus_path
    assert path is not None and path.exists()
    saved = json.loads(path.read_text(encoding="utf-8"))
    for key in ("bands", "n_stars", "gaia_release", "config_hash", "sample_size", "sample_kind", "locus_id"):
        assert key in saved
    for band in ("W1", "W2", "W3", "W4"):
        assert {"bp_rp", "median", "scatter", "intrinsic_scatter", "n"} <= set(saved["bands"][band])

    star = get_named_object("4582145718700356096")
    # A fresh channel with no sample loads the persisted locus for the same release and sample kind
    fresh = build_channel("ir_excess", cfg, tmp_path)
    m = fresh.measure(star, {})
    assert m.extra["photosphere_source"] == "persisted_locus"
    assert m.extra["locus_n_stars"] == saved["n_stars"] and m.values["locus_n_stars"] == saved["n_stars"]
    # A live star never reuses a synthetic locus: it falls back to the literature table
    live_like = dict(copy.deepcopy(star), live_enriched=False)
    m_live = fresh.measure(live_like, {})
    assert m_live.extra["photosphere_source"].startswith("literature_fallback")
    # Elsewhere (no persisted locus at all) the literature fallback is used and labelled
    other = build_channel("ir_excess", cfg, tmp_path / "empty")
    m2 = other.measure(star, {})
    assert m2.extra["photosphere_source"] == literature_locus(cfg.channels["ir_excess"])["source"]
    assert m2.extra["locus_n_stars"] == 0
    assert any("dwarf colour table" in c for c in m2.caveats)


def test_fixture_open_case_dust_recovered_on_empirical_photosphere(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    ch = build_channel("ir_excess", cfg, tmp_path)
    ch.set_empirical_locus(generate_field_sample(1500, seed=cfg.sample.seed))
    m = ch.measure(get_named_object("4582145718700356096"), {})
    v = m.values
    assert m.extra["photosphere_source"] == "sample_locus"
    # Fixture injected T_dust = 226 K, L_dust/L_star = 0.0132
    assert v["t_dust_lo_k"] <= 226.0 <= v["t_dust_hi_k"] or abs(v["t_dust_k"] - 226.0) < 15.0
    assert 0.011 <= v["frac_lum"] <= 0.016
    assert v["frac_lum_lo"] <= v["frac_lum"] <= v["frac_lum_hi"]


def test_candidate_regression_w1_w2_excess_and_labels(tmp_path: Path) -> None:
    """Gaia DR3 5628404318169472: K-dwarf W1/W2 excess is ~3 and ~6 mJy on dwarf colours, not 1.3 and 3.8 (Planck)."""
    cfg, _, _ = load_config()
    star = candidate_star()
    assert star["ebv_star_source"] == "gaia_gspphot" and abs(star["ebv_star"] - 0.039 / 1.31) < 1e-9
    assert star["ebv_sfd"] is None and star["ebv_planck"] is None

    ch = build_channel("ir_excess", cfg, tmp_path)  # no sample, nothing persisted -> literature fallback
    m = ch.measure(star, {})
    v, u = m.values, m.uncertainties
    assert m.extra["photosphere_source"] in {"literature_fallback", "literature_fallback_unverified"}
    assert m.extra["locus_n_stars"] == 0
    assert m.extra["ebv_star_source"] == "gaia_gspphot"

    assert 2.5 <= v["w1_excess_mjy"] <= 3.5
    assert 5.5 <= v["w2_excess_mjy"] <= 6.7
    assert v["w1_excess_mjy"] > 2.0 * 1.29 * 0.9 and v["w2_excess_mjy"] > 1.4 * 3.81
    for b in ("w1", "w2", "w3", "w4"):
        assert u[f"{b}_phot_mjy"] > 0.0 and u[f"{b}_excess_mjy"] > u[f"{b}_obs_mjy"]
    # Excess ratios are higher than the biased Planck values (5.54x, 18.31x)
    assert v["w3_excess_ratio"] > 1.15 * 5.54
    assert v["w4_excess_ratio"] > 1.15 * 18.31
    # Bare Planck diagnostic is kept and is brighter than the empirical photosphere at W1/W2
    assert v["w1_pred_planck_mjy"] > v["w1_phot_mjy"] * 0.99
    assert v["w2_pred_planck_mjy"] > 1.15 * v["w2_phot_mjy"]
    # Dust temperatures: single blackbody with an interval, plus distinct W1-W2 and W3-W4 colour temperatures
    assert v["t_dust_lo_k"] < v["t_dust_k"] < v["t_dust_hi_k"]
    assert v["t_color_w1w2_valid"] == 1.0 and v["t_color_w3w4_valid"] == 1.0
    assert v["t_color_w1w2_k"] > 1.5 * v["t_color_w3w4_k"]
    assert v["t_color_w1w2_lo_k"] < v["t_color_w1w2_k"] < v["t_color_w1w2_hi_k"]
    assert v["frac_lum"] > 0.0 and v["frac_lum_two_temp"] > 0.0
    # Still detected under the unchanged two-predictor rule
    assert ch.score(m).flagged
    assert "photosphere=literature_fallback" in ch.score(m).summary

    # Two-case extinction summary on the measurement
    ext = m.extra["extinction_cases"]
    assert ext is not None
    cases = {c["case"]: c for c in ext["cases"]}
    assert set(cases) == {"applied", "full_line_of_sight"}
    assert cases["applied"]["spectral_type"].startswith("K4") or cases["applied"]["spectral_type"].startswith("K3")
    assert cases["full_line_of_sight"]["spectral_type"] in {"K1V", "K2V", "K3V"}
    assert cases["full_line_of_sight"]["bp_rp_0"] < cases["applied"]["bp_rp_0"] - 0.15
    assert cases["full_line_of_sight"]["m_g"] < cases["applied"]["m_g"] - 0.3


def test_candidate_dossier_carries_labels_and_extinction_cases(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    star = candidate_star()
    star["lightcurves"] = {}
    chans = [build_channel(n, cfg, tmp_path) for n in ("ir_excess", "underluminous")]
    ms = [c.measure(star, {}) for c in chans]
    scs = [c.score(m) for c, m in zip(chans, ms)]
    under = ms[1]
    assert under.extra["ebv_star_source"] == "gaia_gspphot"
    assert under.extra["extinction_cases"] is not None
    dossier, json_p, md_p = write_dossier(
        source_id=str(star["source_id"]),
        star=star,
        measurements=ms,
        scores=scs,
        vetting=[],
        sample_stars=generate_field_sample(50, seed=cfg.sample.seed),
        sample_stats_by_channel={},
        certs_by_channel={},
        reviews=[],
        run_id="test",
        manifests=[],
        cfg=cfg,
        root=tmp_path,
    )
    assert dossier.extinction_cases is not None
    assert dossier.star.ebv_star_source == "gaia_gspphot"
    ir = next(c for c in dossier.channels if c.name == "ir_excess")
    assert ir.labels["photosphere_source"].startswith("literature_fallback")
    md = md_p.read_text(encoding="utf-8")
    assert "Extinction: Two Cases" in md and "full_line_of_sight" in md
    assert "photosphere_source=literature_fallback" in md
    assert "W1-W2 excess colour temperature" in md
