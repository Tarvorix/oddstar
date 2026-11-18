"""Offline tests: every light-curve vetting check evaluates fetched data or returns not_evaluated, and can fail.

No network: archive responses are fakes, SkyBoT goes through the monkeypatched wrapper
`oddstar.io.optical_timeseries.skybot_cone_search`, and the regression test reads the published review evidence in
`review/evidence/2026-11-17/` for Gaia DR3 5628404318169472.
"""

from __future__ import annotations

import csv
import io
import json
import math
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import requests

from oddstar.config import load_config, repo_root
from oddstar.errors import ArchiveAccessError
from oddstar.io import live_survey
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.io.physics import mag_to_flux_mjy
from oddstar.vet.lc_checks import (
    BTJD_TO_MJD,
    check_known_variable,
    check_neighbour_signal,
    check_pipeline_agreement,
    check_second_survey,
    check_sso_crossing,
    check_systematics_window,
    run_lc_vetting,
)

RA, DEC, PMRA, PMDEC = 43.548272861, 5.882712239, 37.914, -97.175
EVIDENCE = repo_root() / "review" / "evidence" / "2026-11-17"


@pytest.fixture(scope="module")
def cfg() -> Any:
    c, _, _ = load_config()
    return c


# ------------------------------------------------------------------------------------------ builders


def tess_sector(
    rng: np.random.Generator,
    *,
    sector: int = 31,
    dip: float = 0.0,
    dip_t: float = 2150.5,
    dip_len: float = 0.75,
    dip_on: str = "target",
    cut_dip_t: float | None = None,
    q: float = 0.3,
    nb: tuple[float, float] = (1.0, 0.0),
    bkg_spike_t: float | None = None,
    jump_t: float | None = None,
    extra: list[tuple[float, float]] | None = None,
) -> dict[str, Any]:
    """Two-orbit 30-min sector (600 s native cadence) with background, centroid and one Gaia neighbour.

    `dip` is the fractional depth of the total (target + neighbour) flux, applied to the target or the neighbour.
    `cut_dip_t` puts the TESScut dip at another time than the QLP one (a non-coincident pair).
    """
    t = np.concatenate([np.arange(2144.5, 2156.5, 1 / 48), np.arange(2158.0, 2170.0, 1 / 48)])
    if extra:
        t = np.sort(np.concatenate([t, [e[0] for e in extra]]))

    def fluxes(t0: float) -> tuple[np.ndarray, np.ndarray]:
        ft, fn = np.ones(len(t)), np.full(len(t), q)
        in_dip = np.abs(t - t0) <= 0.5 * dip_len
        loss = dip * (1.0 + q)
        if dip_on == "target":
            ft[in_dip] -= loss
        else:
            fn[in_dip] -= loss
        return ft, fn

    ft, fn = fluxes(dip_t)
    tot_q = ft + fn
    ft_c, fn_c = fluxes(dip_t if cut_dip_t is None else cut_dip_t)
    tot_c = ft_c + fn_c
    fq = tot_q * (1.0 + rng.normal(0, 0.0008, len(t)))
    fc = tot_c * (1.0 + rng.normal(0, 0.0012, len(t)))
    if extra:
        for te, depth in extra:
            fq[np.argmin(np.abs(t - te))] *= 1.0 - depth
    cx = fn_c * nb[0] / tot_c + rng.normal(0, 0.002, len(t))
    cy = fn_c * nb[1] / tot_c + rng.normal(0, 0.002, len(t))
    if jump_t is not None:
        cx[t >= jump_t] += 0.2
    bkg = 0.01 + rng.normal(0, 0.0005, len(t))
    if bkg_spike_t is not None:
        bkg[np.abs(t - bkg_spike_t) <= 0.2] += 0.05
    return {
        "sector": sector,
        "time_btjd": t.tolist(),
        "flux_qlp_raw": (fq / np.median(fq)).tolist(),
        "flux_tesscut": (fc / np.median(fc)).tolist(),
        "n_per_bin": [3] * len(t),
        "cadence_days": 600.0 / 86400.0,
        "bkg_frac": bkg.tolist(),
        "centroid_dcol_px": cx.tolist(),
        "centroid_drow_px": cy.tolist(),
        "target_px": [6.0, 6.0],
        "neighbour_px": [{"id": "Gaia DR3 1", "dcol_px": nb[0], "drow_px": nb[1], "flux_ratio": q}],
    }


def live_star(**lc: Any) -> dict[str, Any]:
    return {"source_id": "1", "ra_deg": RA, "dec_deg": DEC, "pmra_mas_yr": PMRA, "pmdec_mas_yr": PMDEC, "mags": {"G": 13.4},
            "simbad_type": "Star", "lightcurves": dict(lc), "neowise": {}}


def reference_neowise() -> dict[str, Any]:
    rows = list(csv.DictReader(open(EVIDENCE / "5628404318169472-neowise-visits.csv", encoding="utf-8")))
    w1 = [mag_to_flux_mjy("W1", float(r["w1_median_vega_mag"])) for r in rows]
    w2 = [mag_to_flux_mjy("W2", float(r["w2_median_vega_mag"])) for r in rows]
    # Visit flux errors from the per-visit MAD (as the live enrichment does), floored at 0.5%
    e1 = [f * max(0.005, 1.2533 * float(r["w1_mad_mag"]) / math.sqrt(int(r["n_frames"]))) / 1.0857 for f, r in zip(w1, rows)]
    e2 = [f * max(0.005, 1.2533 * float(r["w2_mad_mag"]) / math.sqrt(int(r["n_frames"]))) / 1.0857 for f, r in zip(w2, rows)]
    return {
        "mjd": [float(r["mjd"]) for r in rows], "w1_mjy": w1, "w2_mjy": w2, "w1_err_mjy": e1, "w2_err_mjy": e2,
        "pm_residual_ra_arcsec": [float(r["ra_cos_dec_residual_arcsec"]) for r in rows],
        "pm_residual_dec_arcsec": [float(r["dec_residual_arcsec"]) for r in rows],
    }


def reference_surveys() -> dict[str, Any]:
    """ZTF yearly medians (dominant object ID per band, as `summarize_survey` keeps) and ASAS-SN camera-subtracted seasons."""
    zrows = list(csv.DictReader(open(EVIDENCE / "5628404318169472-ztf-annual.csv", encoding="utf-8")))
    zbands: dict[str, Any] = {}
    for band in ("zg", "zr"):
        br = [r for r in zrows if r["filtercode"] == band]
        main = max({r["oid"] for r in br}, key=lambda o: sum(1 for r in br if r["oid"] == o))
        seasons = [{"season": int(r["year"]), "n_nights": int(r["n"]), "median_mag": float(r["median_mag"]),
                    "median_err_mag": 1.2533 * float(r["sample_std_mag"]) / math.sqrt(int(r["n"])), "mjd_mid": float(r["visit_median_mjd"])}
                   for r in br if r["oid"] == main]
        zbands[band] = {"n_epochs": sum(s["n_nights"] for s in seasons), "seasons": seasons, "object_id": main}
    arows = [r for r in csv.DictReader(open(EVIDENCE / "asassn-seasonal.csv", encoding="utf-8")) if r["target"] == "5628404318169472"]
    abands: dict[str, Any] = {}
    for band in ("V", "g"):
        seasons = [{"season": int(r["season_start_year"]), "n_nights": int(r["n_nights"]),
                    "median_mag": float(r["camera_subtracted_nightly_median_mag"]),
                    "median_err_mag": 1.2533 * float(r["camera_subtracted_nightly_robust_sigma_mag"]) / math.sqrt(int(r["n_nights"]))}
                   for r in arows if r["band"] == band]
        abands[band] = {"n_epochs": sum(int(r["n_exposures"]) for r in arows if r["band"] == band), "seasons": seasons}
    return {
        "ztf": {"survey": "ZTF", "status": "ok", "season_kind": "calendar", "n_epochs": 1, "bands": zbands},
        "asassn": {"survey": "ASAS-SN", "status": "ok", "season_kind": "july_june", "n_epochs": 1, "bands": abands},
    }


def reference_star() -> dict[str, Any]:
    star = live_star(optical_surveys=reference_surveys())
    star["source_id"] = "5628404318169472"
    star["neowise"] = reference_neowise()
    return star


def step_neowise(rng: np.random.Generator, *, amp_mag: float = 0.3, n: int = 21) -> dict[str, Any]:
    mjd = 56700.0 + np.arange(n) * 182.6
    m = np.where(np.arange(n) >= n // 2, -amp_mag, 0.0) + rng.normal(0, 0.01, n)
    w2 = 100.0 * 10 ** (-0.4 * m)
    return {"mjd": mjd.tolist(), "w1_mjy": (w2 * 1.3).tolist(), "w2_mjy": w2.tolist(),
            "w1_err_mjy": (w2 * 0.01).tolist(), "w2_err_mjy": (w2 * 0.01).tolist()}


def seasons_from(mags_by_year: dict[int, float], err: float = 0.004, n: int = 40) -> list[dict[str, Any]]:
    out = []
    for y, m in sorted(mags_by_year.items()):
        lo = 51544.5 + (y - 2000) * 365.25
        out.append({"season": y, "n_nights": n, "median_mag": m, "median_err_mag": err, "mjd_first": lo + 30, "mjd_mid": lo + 180, "mjd_last": lo + 330})
    return out


class _Resp:
    def __init__(self, status_code: int, content: bytes = b"") -> None:
        self.status_code = status_code
        self.content = content
        self.text = content.decode("utf-8", errors="replace")


# ------------------------------------------------------------------------------- pipeline_agreement


def test_pipeline_agreement_requires_time_coincidence_and_filled_bins(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(1)
    # Moved literals keep their values
    od = cfg.channels["optical_dips"]
    assert (od["pipeline_quiet_depth_frac"], od["pipeline_depth_abs_tol_frac"], od["pipeline_depth_min_ratio"]) == (0.02, 0.08, 0.5)

    same = check_pipeline_agreement(live_star(tess_sectors=[tess_sector(rng, dip=0.06)]), cfg, tmp_path / "s")
    assert same.status == "pass" and same.evidence
    dc = same.metrics["deepest_coincident"]
    assert dc["sector"] == 31 and abs(dc["btjd"] - 2150.5) < 0.3 and dc["coincident_exceeds_min_depth"]
    assert "sector 31 at BTJD 2150" in same.note and "same TESS FFI pixels" in same.note
    sec = same.metrics["sectors"][0]
    for key in ("btjd", "depth_qlp", "depth_tesscut", "n_cadences", "edge_distance_days", "coincident_exceeds_min_depth"):
        assert key in sec["deepest_qlp"] and key in sec["deepest_tesscut"]
    assert sec["binned_pearson_r"] is not None and sec["correlation_used"] is True

    # Equal 6% dips at different times: the deepest QLP bin has no TESScut counterpart -> fail
    apart = check_pipeline_agreement(live_star(tess_sectors=[tess_sector(rng, dip=0.06, cut_dip_t=2163.0)]), cfg, tmp_path / "a")
    assert apart.status == "fail" and "sector(s) [31]" in apart.note and "not matched by TESScut" in apart.note

    # A 20% QLP-only value in a bin holding 3 of 36 expected cadences is dropped, not counted as a dip
    lone = tess_sector(rng, extra=[(2156.78, 0.20)])
    res = check_pipeline_agreement(live_star(tess_sectors=[lone]), cfg, tmp_path / "l")
    assert res.status == "pass", res.note
    assert res.metrics["sectors"][0]["n_bins_dropped_underfilled"] >= 1


# ------------------------------------------------------------------------------- systematics_window


def test_systematics_window_tess_and_neowise(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(2)
    # No data at all (live star without TESS and NEOWISE) -> not_evaluated, never pass
    none = check_systematics_window(live_star(), cfg, tmp_path / "n")
    assert none.status == "not_evaluated" and none.evidence is None

    # Query failures are named in the note
    failed = live_star(tess_sectors=[], tess_status="query_failed", tess_query_log=["s4: query_failed ArchiveAccessError: HTTP 503"])
    failed["neowise"] = dict(step_neowise(rng), comparison={"status": "query_failed", "n_stars": 0, "query_log": ["ArchiveAccessError: Gator HTTP 500"]})
    res_f = check_systematics_window(failed, cfg, tmp_path / "f")
    assert res_f.status == "not_evaluated" and "query_failed" in res_f.note and "HTTP 503" in res_f.note

    # Clean mid-orbit dip -> pass with evidence
    clean = check_systematics_window(live_star(tess_sectors=[tess_sector(rng, dip=0.06)]), cfg, tmp_path / "c")
    assert clean.status == "pass" and clean.evidence and Path(clean.evidence).exists()
    assert clean.metrics["tess"]["events"]

    # Dip coinciding with a scattered-light background excursion -> fail
    bkg = check_systematics_window(live_star(tess_sectors=[tess_sector(rng, dip=0.06, bkg_spike_t=2150.5)]), cfg, tmp_path / "b")
    assert bkg.status == "fail" and "background excursion" in bkg.note
    # Dip within the edge margin of the orbit end (not running into it: the matched running median with
    # mode="nearest" padding absorbs a dip that reaches the last cadence) -> fail
    edge = check_systematics_window(live_star(tess_sectors=[tess_sector(rng, dip=0.06, dip_t=2156.05, dip_len=0.5)]), cfg, tmp_path / "e")
    assert edge.status == "fail" and "segment edge" in edge.note
    # Dip coinciding with a pointing jump -> fail
    jump = check_systematics_window(live_star(tess_sectors=[tess_sector(rng, dip=0.06, jump_t=2150.5)]), cfg, tmp_path / "j")
    assert jump.status == "fail" and "centroid jump" in jump.note

    # NEOWISE: target W2 step shared by the field ensemble -> scan/zero-point systematic -> fail
    neo = step_neowise(rng)
    tgt = -2.5 * np.log10(np.asarray(neo["w2_mjy"]) / np.median(neo["w2_mjy"]))
    shared = live_star()
    shared["neowise"] = dict(neo, comparison={"status": "ok", "n_stars": 5, "ensemble_dw2_mag": [round(float(v), 4) for v in tgt],
                                              "ensemble_dw2_err_mag": [0.005] * len(tgt)})
    res_s = check_systematics_window(shared, cfg, tmp_path / "ns")
    assert res_s.status == "fail" and "scan/zero-point systematic" in res_s.note
    assert res_s.metrics["neowise"]["metrics"]["n_comparison_stars"] == 5
    # Flat ensemble -> pass
    flat = live_star()
    flat["neowise"] = dict(neo, comparison={"status": "ok", "n_stars": 5,
                                            "ensemble_dw2_mag": [round(float(v), 4) for v in rng.normal(0, 0.005, len(tgt))],
                                            "ensemble_dw2_err_mag": [0.005] * len(tgt)})
    res_p = check_systematics_window(flat, cfg, tmp_path / "np")
    assert res_p.status == "pass" and "not shared by field stars" in res_p.note


def test_untested_neowise_variation_cannot_pass_on_tess_alone(tmp_path: Path, cfg: Any) -> None:
    """A star flagged for an infrared change, with no comparison stars, must not pass because its TESS data are clean."""
    rng = np.random.default_rng(5)
    too_few = {"status": "too_few_comparison_stars", "n_stars": 1, "query_log": ['Gator neowiser_p1bs_psd 600" cone: 9000 rows']}
    varying = live_star(tess_sectors=[tess_sector(rng)])
    varying["neowise"] = dict(step_neowise(rng), comparison=dict(too_few))
    res = check_systematics_window(varying, cfg, tmp_path / "v")
    assert res.status == "not_evaluated" and res.note.startswith("The star varies in NEOWISE")
    assert res.metrics["neowise_variation_untested"] is True and res.metrics["tess"]["status"] == "pass"
    assert "too_few_comparison_stars (1 stars)" in res.note

    # No comparison ensemble at all (fetch never ran) is the same situation
    missing = live_star(tess_sectors=[tess_sector(rng)])
    missing["neowise"] = step_neowise(rng)
    assert check_systematics_window(missing, cfg, tmp_path / "m").status == "not_evaluated"

    # A star that is quiet in NEOWISE has nothing infrared to attribute: the TESS result stands
    quiet = live_star(tess_sectors=[tess_sector(rng)])
    quiet["neowise"] = dict(step_neowise(rng, amp_mag=0.0), comparison=dict(too_few))
    res_q = check_systematics_window(quiet, cfg, tmp_path / "q")
    assert res_q.status == "pass" and res_q.metrics["neowise_variation_untested"] is False

    # A failing TESS part still fails, and an explicit fixture flag still counts as an evaluation
    bad = live_star(tess_sectors=[tess_sector(rng, dip=0.06, bkg_spike_t=2150.5)])
    bad["neowise"] = dict(step_neowise(rng), comparison=dict(too_few))
    assert check_systematics_window(bad, cfg, tmp_path / "b").status == "fail"
    flagged = live_star(tess_sectors=[tess_sector(rng)])
    flagged["neowise"] = dict(step_neowise(rng), scan_systematic=False)
    assert check_systematics_window(flagged, cfg, tmp_path / "f").status == "pass"


def test_build_neowise_comparison_groups_stars_on_visit_grid() -> None:
    rng = np.random.default_rng(3)
    windows = [(56700.0 + 182.6 * k, 56701.5 + 182.6 * k) for k in range(12)]
    rows: list[dict[str, Any]] = []
    stars = [(40.0, 10.0, 10.6), (-60.0, 30.0, 10.4), (20.0, -80.0, 10.9), (90.0, 90.0, 10.5), (0.5, 0.3, 10.5), (70.0, -20.0, 7.0)]
    for k, (lo, _) in enumerate(windows):
        for sx, sy, mag in stars:
            for _ in range(8):
                rows.append({"ra": RA + (sx + rng.normal(0, 0.2)) / 3600.0 / math.cos(math.radians(DEC)),
                             "dec": DEC + (sy + rng.normal(0, 0.2)) / 3600.0, "mjd": lo + rng.uniform(0, 1.5),
                             "w1mpro": mag + rng.normal(0, 0.02), "w2mpro": mag - 0.05 + (0.1 if k >= 6 else 0.0) + rng.normal(0, 0.02),
                             "qual_frame": 10})
    comp = ot.build_neowise_comparison(rows, ra=RA, dec=DEC, visit_windows=windows, target_w1_mag=10.5, exclude_arcsec=8.0,
                                       mag_window=1.0, min_stars=3, max_stars=10)
    # The target (0.6") and the W1=7 star are excluded; four comparison stars remain
    assert comp["status"] == "ok" and comp["n_stars"] == 4
    ens = np.asarray(comp["ensemble_dw2_mag"], dtype=float)
    assert len(ens) == 12 and abs(float(np.mean(ens[6:]) - np.mean(ens[:6])) - 0.1) < 0.03
    assert all(e is not None for e in comp["ensemble_dw2_err_mag"])


def test_variation_carried_by_saturated_neowise_bands_fails(tmp_path: Path, cfg: Any) -> None:
    """Config 1.7.0: brighter than W1 = 8.0 / W2 = 7.0 the single-exposure fluxes carry a bias the quoted errors do not hold."""
    from oddstar.fuse.engine import assign_tier
    from oddstar.schemas import ChannelScore
    from oddstar.vet import vet_target

    rng = np.random.default_rng(11)
    neo = step_neowise(rng)                                   # median W2 8.09, W1 8.44: not saturated, a 0.3 mag step in both
    bright = {k: ([v * 40.0 for v in vals] if k.endswith("_mjy") else vals) for k, vals in neo.items()}   # 4 mag brighter
    # Saturated in both bands, and no comparison stars exist for a star this bright: the check fails on saturation
    star = live_star()
    star["neowise"] = dict(bright, comparison={"status": "too_few_comparison_stars", "n_stars": 0, "query_log": ["w1mpro between 3.98 and 6.48: 255 rows"]})
    res = check_systematics_window(star, cfg, tmp_path / "sat")
    assert res.status == "fail" and "saturated photometry" in res.note and "W2 3.80 (limit 7.0)" in res.note
    sat = res.metrics["neowise"]["saturation"]
    assert sat["bands_varying"] == ["w1", "w2"] and sat["saturated"] == {"w1": True, "w2": True}
    # The same star without the brightening is `not_evaluated` as before (no comparison stars), never a pass
    faint = live_star()
    faint["neowise"] = dict(neo, comparison={"status": "too_few_comparison_stars", "n_stars": 0})
    assert check_systematics_window(faint, cfg, tmp_path / "unsat").status == "not_evaluated"
    # W1 saturated, W2 not, and the variation is in both: W2 still carries it, so saturation does not decide
    mixed = dict(neo, w1_mjy=[v * 2.0 for v in neo["w1_mjy"]], w1_err_mjy=[v * 2.0 for v in neo["w1_err_mjy"]])   # W1 7.69
    half = live_star()
    half["neowise"] = dict(mixed, comparison={"status": "too_few_comparison_stars", "n_stars": 0})
    res_h = check_systematics_window(half, cfg, tmp_path / "half")
    assert res_h.status == "not_evaluated" and "saturated photometry" not in res_h.note
    # A saturated star that does not vary is not failed
    flat = {k: ([float(np.median(vals))] * len(vals) if k in ("w1_mjy", "w2_mjy") else vals) for k, vals in bright.items()}
    quiet = live_star()
    quiet["neowise"] = flat
    assert "saturated" not in check_systematics_window(quiet, cfg, tmp_path / "quiet").note
    # In the tiers: a failed systematics_window is an artifact, so the flag is Rejected
    score = ChannelScore(source_id="1", channel="ir_timedomain", raw_statistic=30.0, false_alarm_prob=1e-6, flagged=True, rank_score=30.0, summary="x")
    star.update({"abs_g_mag": 3.4, "teff_k": 5009.0, "teff_k_source": "gspphot"})
    tier, why = assign_tier([score], vet_target(star, [score], cfg, tmp_path / "tier"), cfg)
    assert tier == "Rejected" and "systematics_window" in why


def test_infrared_change_of_a_luminous_giant_is_explained(tmp_path: Path, cfg: Any) -> None:
    """Config 1.7.0: an ir_timedomain flag also gets the `evolved` check; before, a variable giant could reach Tier B."""
    from oddstar.fuse.engine import assign_tier
    from oddstar.schemas import ChannelScore
    from oddstar.vet import vet_target
    from oddstar.vet.ir_checks import check_evolved

    rng = np.random.default_rng(12)
    score = ChannelScore(source_id="1", channel="ir_timedomain", raw_statistic=30.0, false_alarm_prob=1e-6, flagged=True, rank_score=30.0, summary="x")

    def tier_of(tag: str, **hr: Any) -> tuple[str, str | None, dict[str, str]]:
        star = live_star()
        star["neowise"] = step_neowise(rng)
        star.update(hr)
        vet = vet_target(star, [score], cfg, tmp_path / tag)
        return (*assign_tier([score], vet, cfg), {v.check: v.status for v in vet})

    # An M giant (M_G -2.9, 3600 K): long-period variability is ordinary
    tier, why, checks = tier_of("m", abs_g_mag=-2.9, teff_k=3612.0, teff_k_source="gspphot")
    assert checks["evolved"] == "fail" and tier == "Explained" and "evolved" in why
    # A luminous star that the catalogue calls 5,600 K behind 4 mag of extinction is no main-sequence star either:
    # at M_G < 0.5 the main sequence is hotter than 9,000 K. The limit was 5,200 K until 1.7.0.
    tier, why, checks = tier_of("g", abs_g_mag=-4.8, teff_k=5597.0, teff_k_source="gspphot")
    assert checks["evolved"] == "fail" and tier == "Explained"
    # A dwarf passes `evolved`; an A star at M_G 0.3 is on the main sequence and passes too
    assert tier_of("k", abs_g_mag=6.5, teff_k=4510.0, teff_k_source="gspphot")[2]["evolved"] == "pass"
    hot = live_star()
    hot.update({"abs_g_mag": 0.3, "teff_k": 9800.0, "teff_k_source": "gspphot"})
    res = check_evolved(hot, cfg, tmp_path / "a")
    assert res.status == "pass" and "not a luminous giant" in res.note
    # No luminosity, or a giant-bright star with a default temperature: `not_evaluated`, so the tier is capped at C
    tier, why, checks = tier_of("none")
    assert checks["evolved"] == "not_evaluated" and tier == "C"
    tier, why, checks = tier_of("def", abs_g_mag=-1.0, teff_k=5500.0, teff_k_source="default_5500K")
    assert checks["evolved"] == "not_evaluated" and tier == "C"


# --------------------------------------------------------------------------------- neighbour_signal


def test_neighbour_signal_tess_centroid(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(4)
    none = check_neighbour_signal(live_star(), cfg, tmp_path / "n")
    assert none.status == "not_evaluated"
    failed = check_neighbour_signal(live_star(tess_sectors=[], tess_status="query_failed", tess_query_log=["TIC lookup HTTP 503"]), cfg, tmp_path / "f")
    assert failed.status == "not_evaluated" and "HTTP 503" in failed.note

    on_target = live_star(tess_sectors=[tess_sector(rng, dip=0.06, dip_on="target")],
                          tess_neighbours={"status": "ok", "radius_arcsec": 42.0, "sources": [{"id": "Gaia DR3 1", "sep_arcsec": 21.0, "flux_ratio": 0.3}]})
    res_t = check_neighbour_signal(on_target, cfg, tmp_path / "t")
    assert res_t.status == "pass" and res_t.evidence, res_t.note
    assert res_t.metrics["tess"]["events"][0]["deficit_offset_px"] < 0.1
    assert "1 source(s)" in res_t.note

    on_nb = live_star(tess_sectors=[tess_sector(rng, dip=0.06, dip_on="neighbour")])
    res_n = check_neighbour_signal(on_nb, cfg, tmp_path / "nb")
    assert res_n.status == "fail" and "Gaia DR3 1" in res_n.note
    ev = res_n.metrics["tess"]["events"][0]
    assert abs(ev["deficit_location_px"][0] - 1.0) < 0.15 and ev["nearest_neighbour"]["dist_px"] < 0.15

    quiet = check_neighbour_signal(live_star(tess_sectors=[tess_sector(rng)]), cfg, tmp_path / "q")
    assert quiet.status == "pass" and "no dip event above" in quiet.note


def test_neighbour_signal_neowise_reference_star_moves_with_its_flux(tmp_path: Path, cfg: Any) -> None:
    star = reference_star()
    res = check_neighbour_signal(star, cfg, tmp_path / "ref")
    assert res.status == "pass", res.note
    m = res.metrics["neowise"]["metrics"]
    assert abs(m["median_residual_arcsec"] - 0.053) < 0.003
    assert m["static_blend_frac"] < 0.1 and m["residual_w2_sigma"] < 3.0
    assert m["stationary_source_pull_predicted_arcsec"] > 0.25 and m["measured_pull_arcsec"] < 0.05
    assert "moves with the star" in res.note

    # Synthetic counter-case: the W2 brightening sits on a stationary source at the star's 2010.5 position
    rng = np.random.default_rng(5)
    mjd = np.asarray(star["neowise"]["mjd"])
    jyr = 2000.0 + (mjd - 51544.5) / 365.25
    blend = np.where(jyr > 2017.0, 0.4, 0.0)
    frac = blend / (1.0 + blend)
    mu = np.asarray([PMRA, PMDEC]) / 1000.0
    resid = -frac[:, None] * (jyr[:, None] - 2010.5) * mu[None, :] + rng.normal(0, 0.015, (len(mjd), 2))
    blended = reference_star()
    blended["neowise"]["w2_mjy"] = (100.0 * (1.0 + blend)).tolist()
    blended["neowise"]["pm_residual_ra_arcsec"] = resid[:, 0].tolist()
    blended["neowise"]["pm_residual_dec_arcsec"] = resid[:, 1].tolist()
    res_b = check_neighbour_signal(blended, cfg, tmp_path / "blend")
    assert res_b.status == "fail" and "blended stationary source" in res_b.note


# ------------------------------------------------------------------------------------- sso_crossing


def test_sso_crossing_statuses(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(6)
    assert check_sso_crossing(live_star(), cfg, tmp_path / "n").status == "not_evaluated"
    # TESS loaded, no event above threshold: legitimate pass (see the comment in check_sso_crossing)
    quiet = check_sso_crossing(live_star(tess_sectors=[tess_sector(rng), tess_sector(rng, sector=32)]), cfg, tmp_path / "q")
    assert quiet.status == "pass" and "No dip event above threshold" in quiet.note and "2 TESS sector(s)" in quiet.note and quiet.evidence

    sec = tess_sector(rng, dip=0.06)
    # Event present but no SkyBoT query recorded -> not_evaluated
    assert check_sso_crossing(live_star(tess_sectors=[sec]), cfg, tmp_path / "m").status == "not_evaluated"

    from oddstar.vet.lc_checks import tess_dip_events

    events, _ = tess_dip_events({"tess_sectors": [sec]}, cfg)
    assert len(events) == 1
    calls: list[dict[str, Any]] = []

    def fake_ok(**kw: Any) -> Any:
        calls.append(kw)
        return [{"name": "2001 AB", "ra_deg": RA, "dec_deg": DEC, "sep_arcsec": 30.0, "v_mag": 21.5}], None

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    import oddstar.io.optical_timeseries as otm

    orig = otm.skybot_cone_search
    try:
        otm.skybot_cone_search = fake_ok  # type: ignore[assignment]
        sso, _ = ot.query_sso_at_events(events, ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, radius_arcsec=42.0, cache=cache)
    finally:
        otm.skybot_cone_search = orig  # type: ignore[assignment]
    assert len(calls) == 3 and all(c["observer"] == "C57" and c["radius_arcsec"] == 42.0 for c in calls)
    assert abs(calls[1]["epoch_jd"] - (events[0]["btjd"] + 2457000.0)) < 1e-6
    faint = check_sso_crossing(live_star(tess_sectors=[sec], sso_queries=sso), cfg, tmp_path / "p")
    assert faint.status == "pass" and faint.evidence  # V=21.5 body is far too faint to matter for a G=13.4 star

    bright = json.loads(json.dumps(sso))
    bright["events"][0]["bodies"][0]["v_mag"] = 13.0
    res_b = check_sso_crossing(live_star(tess_sectors=[sec], sso_queries=bright), cfg, tmp_path / "b")
    assert res_b.status == "fail" and "2001 AB" in res_b.note

    broken = json.loads(json.dumps(sso))
    broken["events"][0].update(status="query_failed", query_log=["JD 2459150.0: ArchiveAccessError: SkyBoT HTTP 502"])
    res_q = check_sso_crossing(live_star(tess_sectors=[sec], sso_queries=broken), cfg, tmp_path / "x")
    assert res_q.status == "not_evaluated" and "HTTP 502" in res_q.note


def test_enrichment_records_failures_and_queries_skybot(tmp_path: Path, cfg: Any, monkeypatch: Any) -> None:
    rng = np.random.default_rng(7)
    sec = tess_sector(rng, dip=0.06)

    def ztf_fail(**kw: Any) -> Any:
        raise ArchiveAccessError("IRSA ZTF HTTP 503", retryable=True)

    def tap_fail(**kw: Any) -> Any:
        raise requests.exceptions.ConnectionError("Gaia TAP unreachable")

    monkeypatch.setattr(live_survey, "fetch_ztf_lightcurve", ztf_fail)
    monkeypatch.setattr(live_survey, "fetch_asassn_lightcurve", lambda **kw: (None, []))
    monkeypatch.setattr(live_survey, "_tap_query_cached", tap_fail)
    monkeypatch.setattr(live_survey, "fetch_tess_pipeline_pair",
                        lambda **kw: ({"tic_id": "1", "sectors": [sec], "status": "ok", "query_log": ["s31: ok"]}, []))

    def skybot_down(**kw: Any) -> Any:
        raise ArchiveAccessError("SkyBoT HTTP 502", retryable=True)

    monkeypatch.setattr(ot, "skybot_cone_search", skybot_down)
    star = live_star()
    star["source_id"] = "1"
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    live_survey._enrich_optical_timeseries(star, cache, cfg)
    lc = star["lightcurves"]
    assert lc["optical_surveys"]["ztf"]["status"] == "query_failed" and "HTTP 503" in lc["optical_surveys"]["ztf"]["query_log"][0]
    assert lc["optical_surveys"]["asassn"]["status"] == "client_unavailable"
    assert lc["tess_neighbours"]["status"] == "query_failed" and lc["tess_neighbours"]["radius_arcsec"] == 42.0
    assert lc["tess_status"] == "ok" and lc["sso_queries"]["status"] == "query_failed"
    res = check_sso_crossing(star, cfg, tmp_path / "ev")
    assert res.status == "not_evaluated" and "SkyBoT HTTP 502" in res.note
    surv = check_second_survey(star, cfg, tmp_path / "ev2")
    assert surv.status == "not_evaluated" and "query_failed" in surv.note


# ------------------------------------------------------------------------------------ known_variable


def _periodic_sector(rng: np.random.Generator, kind: str) -> dict[str, Any]:
    t = np.concatenate([np.arange(2144.5, 2156.5, 1 / 48), np.arange(2158.0, 2170.0, 1 / 48)])
    f = 1.0 + rng.normal(0, 0.001, len(t))
    if kind == "rotation":
        f += 0.01 * np.sin(2 * np.pi * t / 3.2)
    elif kind == "eb":
        ph = ((t - 2145.0) % 2.1) / 2.1
        f[np.abs(ph - 0.5) < 0.02] -= 0.25
    elif kind == "single_dip":
        f[np.abs(t - 2150.0) < 0.5] -= 0.25
    return {"sector": 31, "time_btjd": t.tolist(), "flux_qlp_raw": f.tolist(), "flux_tesscut": f.tolist()}


def test_known_variable_rotation_reported_eb_fails(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(8)
    nothing = {"source_id": "1", "lightcurves": {}}
    assert check_known_variable(nothing, cfg, tmp_path / "n").status == "not_evaluated"

    rot = check_known_variable(live_star(tess_sectors=[_periodic_sector(rng, "rotation")]), cfg, tmp_path / "r")
    assert rot.status == "pass" and rot.evidence
    assert abs(rot.metrics["rotation"]["period_days"] - 3.2) < 0.15 and rot.metrics["rotation"]["best_fap"] < 1e-3
    assert "rotation/spot modulation candidate" in rot.note

    eb = check_known_variable(live_star(tess_sectors=[_periodic_sector(rng, "eb")]), cfg, tmp_path / "e")
    assert eb.status == "fail" and eb.metrics["periodic_eb_sectors"] == [31]

    dip = check_known_variable(live_star(tess_sectors=[_periodic_sector(rng, "single_dip")]), cfg, tmp_path / "d")
    assert dip.status == "pass" and not dip.metrics["periodic_eb_detected"]


# ------------------------------------------------------------------------------------- second_survey


def test_second_survey_reference_star_dims_while_ir_brightens(tmp_path: Path, cfg: Any) -> None:
    res = check_second_survey(reference_star(), cfg, tmp_path / "ref")
    assert res.status == "pass", res.note
    assert res.metrics["optical_response"] == "dimming_with_ir_brightening"
    bands = {(b["survey"], b["band"]): b for b in res.metrics["contexts"]["ir_timedomain"]["bands"]}
    assert abs(bands[("ZTF", "zg")]["pearson_r"] + 0.84) < 0.03 and abs(bands[("ZTF", "zr")]["pearson_r"] + 0.84) < 0.03
    assert bands[("ZTF", "zg")]["slope_dopt_dw2"] < 0 and bands[("ZTF", "zg")]["significant_change"]
    # Coverage of the IR event window is stated per survey (ZTF starts 2018, ASAS-SN V ends 2019)
    assert "ZTF zg covers 2018.0" in res.note and "ASAS-SN V covers 2017.1-2019.5" in res.note
    assert "significant optical change detected" in res.note and "clean" not in res.note.lower()


def test_second_survey_counter_case_and_coverage(tmp_path: Path, cfg: Any) -> None:
    ref = reference_star()
    neo = ref["neowise"]
    # Counter-case: optical brightens by the same magnitude as W2 (photospheric/blended variable) -> fail
    wm = np.asarray(neo["mjd"])
    w2m = -2.5 * np.log10(np.asarray(neo["w2_mjy"]) / np.median(neo["w2_mjy"]))
    years = range(2015, 2024)
    mids = {y: 51544.5 + (y - 2000) * 365.25 + 180 for y in years}
    bright = {y: 14.3 + float(np.interp(mids[y], wm, w2m)) for y in years}
    star = live_star(optical_surveys={"ztf": {"survey": "ZTF", "status": "ok", "season_kind": "calendar", "n_epochs": 900,
                                              "bands": {"zg": {"n_epochs": 900, "seasons": seasons_from(bright)}}}})
    star["neowise"] = neo
    res = check_second_survey(star, cfg, tmp_path / "c")
    assert res.status == "fail" and res.metrics["optical_response"] == "brightening_with_ir_brightening"
    assert "contradicting circumstellar dust" in res.note

    # Flat optical through the IR event -> response none -> pass, never described as clean
    flat = live_star(optical_surveys={"ztf": {"survey": "ZTF", "status": "ok", "season_kind": "calendar", "n_epochs": 900,
                                              "bands": {"zg": {"n_epochs": 900, "seasons": seasons_from({y: 14.3 for y in years})}}}})
    flat["neowise"] = neo
    res_f = check_second_survey(flat, cfg, tmp_path / "f")
    assert res_f.status == "pass" and res_f.metrics["optical_response"] == "none"

    # Seasons only after the IR event window -> not_evaluated naming the window
    late = live_star(optical_surveys={"ztf": {"survey": "ZTF", "status": "ok", "season_kind": "calendar", "n_epochs": 900,
                                              "bands": {"zg": {"n_epochs": 900, "seasons": seasons_from({2023: 14.3, 2024: 14.31})}}}})
    late["neowise"] = neo
    res_l = check_second_survey(late, cfg, tmp_path / "l")
    assert res_l.status == "not_evaluated" and "IR event window 2017.1-2022.0" in res_l.note

    # No survey at all -> not_evaluated
    assert check_second_survey(live_star(), cfg, tmp_path / "n").status == "not_evaluated"


def test_second_survey_simultaneous_tess_dip(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(9)
    sec = tess_sector(rng, dip=0.06)
    ev_mjd = 2150.5 + BTJD_TO_MJD
    base = [[ev_mjd - 40 + k, 14.30 + float(rng.normal(0, 0.005))] for k in range(0, 80, 2)]

    def star_with(night_mag: float | None) -> dict[str, Any]:
        nights = [n for n in base if abs(n[0] - ev_mjd) > 1.0]
        if night_mag is not None:
            nights.append([ev_mjd + 0.05, night_mag])
        nights.sort()
        band = {"n_epochs": len(nights), "seasons": seasons_from({2020: 14.30}), "nights": nights, "nightly_robust_sigma_mag": 0.005}
        return live_star(tess_sectors=[sec], optical_surveys={"ztf": {"survey": "ZTF", "status": "ok", "season_kind": "calendar",
                                                                      "n_epochs": len(nights), "bands": {"zg": band}}})

    res_c = check_second_survey(star_with(14.30), cfg, tmp_path / "c")
    assert res_c.status == "fail" and "not on the sky" in res_c.note
    res_p = check_second_survey(star_with(14.365), cfg, tmp_path / "p")
    assert res_p.status == "pass" and "consistent with TESS" in res_p.note
    res_n = check_second_survey(star_with(None), cfg, tmp_path / "n")
    assert res_n.status == "not_evaluated" and "do not constrain the dip" in res_n.note


# --------------------------------------------------------------------------- fetchers and summaries


def test_summarize_survey_season_errors_and_no_coverage() -> None:
    rng = np.random.default_rng(10)
    rows = [{"oid": "a", "band": "zg", "mjd": 58500.0 + i * 1.3, "mag": 14.3 + rng.normal(0, 0.02), "magerr": 0.02} for i in range(400)]
    s = ot.summarize_survey(rows, survey="ZTF", season="calendar", group_key="oid")
    seas = s["bands"]["zg"]["seasons"]
    assert all(x["median_err_mag"] is not None and 0 < x["median_err_mag"] < 0.01 for x in seas if x["n_nights"] >= 5)
    assert all(x["mjd_first"] <= x["mjd_mid"] <= x["mjd_last"] for x in seas)
    assert len(s["bands"]["zg"]["nights"]) == s["bands"]["zg"]["n_nights"] and s["season_kind"] == "calendar"
    assert ot.summarize_survey([], survey="ZTF", season="calendar", group_key="oid")["status"] == "no_coverage"


def _fake_tess(n: int = 300, bkg_spike: bool = True) -> tuple[bytes, bytes]:
    from astropy.io import fits
    from astropy.wcs import WCS

    t = 1410.0 + np.arange(n) / 48.0
    q = fits.BinTableHDU.from_columns([
        fits.Column(name="TIME", format="D", array=t),
        fits.Column(name="SAP_FLUX", format="E", array=np.full(n, 5000.0)),
        fits.Column(name="QUALITY", format="J", array=np.zeros(n, dtype=int)),
    ])
    qlp = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), q]).writeto(qlp)
    w = WCS(naxis=2)
    w.wcs.crpix = [7.0, 7.0]
    w.wcs.crval = [RA, DEC]
    w.wcs.cdelt = [-21.0 / 3600, 21.0 / 3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    yy, xx = np.mgrid[0:13, 0:13]
    star = 1000.0 * np.exp(-((xx - 6.0) ** 2 + (yy - 6.0) ** 2) / (2 * 0.7**2))
    cube = np.repeat((star + 20.0)[None, :, :], n, axis=0)
    if bkg_spike:
        cube[150:160] += 30.0  # scattered-light excursion: +30 e/s/px everywhere
    tab = fits.BinTableHDU.from_columns([
        fits.Column(name="TIME", format="D", array=t),
        fits.Column(name="FLUX", format="169E", dim="(13,13)", array=cube),
        fits.Column(name="QUALITY", format="J", array=np.zeros(n, dtype=int)),
    ])
    fb = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), tab, fits.ImageHDU(np.ones((13, 13), dtype=np.int32), header=w.to_header())]).writeto(fb)
    zb = io.BytesIO()
    with zipfile.ZipFile(zb, "w") as zf:
        zf.writestr("astrocut.fits", fb.getvalue())
    return qlp.getvalue(), zb.getvalue()


def test_fetch_tess_pair_keeps_background_centroid_and_neighbours(tmp_path: Path, monkeypatch: Any) -> None:
    qlp_bytes, zip_bytes = _fake_tess()
    state = {"cut_status": 200}

    def fake_get(url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        if url == ot.MAST_INVOKE_URL:
            return _Resp(200, json.dumps({"data": [{"ID": "387517877", "ra": RA, "dec": DEC, "GAIA": "1"}]}).encode())
        if url == ot.TESSCUT_SECTOR_URL:
            return _Resp(200, json.dumps({"results": [{"sector": "4"}]}).encode())
        if url.startswith(ot.QLP_S3_BASE):
            return _Resp(200, qlp_bytes)
        if url == ot.TESSCUT_ASTROCUT_URL:
            return _Resp(state["cut_status"], zip_bytes if state["cut_status"] == 200 else b"")
        return _Resp(404)

    monkeypatch.setattr(requests, "get", fake_get)
    nb = [{"id": "Gaia DR3 2", "ra_deg": RA + 21.0 / 3600 / math.cos(math.radians(DEC)), "dec_deg": DEC, "flux_ratio": 0.1}]
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    out, mans = ot.fetch_tess_pipeline_pair(ra=RA, dec=DEC, pmra=0.0, pmdec=0.0, gaia_id="1", cache=cache, neighbours=nb, sleep=lambda s: None)
    assert out is not None and out["status"] == "ok" and mans
    sec = out["sectors"][0]
    n = len(sec["time_btjd"])
    for key in ("bkg_frac", "centroid_dcol_px", "centroid_drow_px", "n_per_bin"):
        assert len(sec[key]) == n
    bkg = np.asarray(sec["bkg_frac"])
    tb = np.asarray(sec["time_btjd"])
    spike = (tb > 1410.0 + 150.5 / 48) & (tb < 1410.0 + 158.5 / 48)
    assert spike.sum() >= 5 and float(bkg[spike].min()) > float(np.median(bkg)) + 0.05  # excursion kept as a background series
    assert abs(float(np.median(sec["centroid_dcol_px"]))) < 0.05 and abs(float(np.median(sec["centroid_drow_px"]))) < 0.05
    # RA increases to the east = decreasing column for this WCS (cdelt1 < 0): one pixel away
    assert abs(sec["neighbour_px"][0]["dcol_px"] + 1.0) < 0.05 and abs(sec["neighbour_px"][0]["drow_px"]) < 0.05

    # Every existing sector failing to download is query_failed, not "no coverage"
    state["cut_status"] = 500
    out2, _ = ot.fetch_tess_pipeline_pair(ra=RA + 1.0, dec=DEC, pmra=0.0, pmdec=0.0, gaia_id="1",
                                          cache=RawCache(tmp_path / "raw2", max_cache_gb=1.0), sleep=lambda s: None)
    assert out2 is not None and out2["status"] == "query_failed" and "query_failed" in out2["query_log"][0]


def test_fixture_flags_still_decide_and_live_absence_is_not_pass(tmp_path: Path, cfg: Any) -> None:
    """Explicit legacy flags are evaluations; a live star with none of the inputs never passes any LC check."""
    fixture_like = {"source_id": "9", "simbad_type": "Star", "mags": {"G": 12.0},
                    "lightcurves": {"neighbour_eb_centroid_shift": False, "tess_scattered_light_window": False, "sso_crossing_t0": None},
                    "neowise": {"scan_systematic": False, "w1_w2_neighbour_blend": True}}
    by = {r.check: r.status for r in run_lc_vetting(fixture_like, cfg, tmp_path / "fx")}
    assert by["neighbour_signal"] == "fail" and by["systematics_window"] == "pass" and by["sso_crossing"] == "pass"

    bare = {"source_id": "8", "ra_deg": RA, "dec_deg": DEC, "lightcurves": {"tess_time_btjd": [], "tess_flux_qlp": [], "tess_flux_tesscut": [],
                                                                           "ztf_time_mjd": [], "ztf_flux": []},
            "neowise": {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []}}
    res = {r.check: r.status for r in run_lc_vetting(bare, cfg, tmp_path / "bare")}
    assert all(s == "not_evaluated" for s in res.values()), res
