"""Remaining default-value passes removed (lessons L007/L008): vanishing, spectral, radio, SIMBAD, generic path, edges.

For every changed check: (a) explicit good inputs pass with an evidence file that records the inputs used,
(b) explicit bad inputs fail, (c) missing inputs give `not_evaluated`, (d) a live star whose query failed (or was
never made) gives `not_evaluated`, even when stale values are present. Offline: archives are fakes, no network.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pytest
import requests

from oddstar.channels import build_channel
from oddstar.config import load_config
from oddstar.io.fixtures_data import RADIO_CADENCES, generate_field_sample, get_named_object
from oddstar.io.physics import propagate_position
from oddstar.schemas import ChannelScore
from oddstar.vet import ir_checks as irc
from oddstar.vet import lc_checks as lcc
from oddstar.vet import radio_checks as rc
from oddstar.vet import spectral_checks as sc
from oddstar.vet import vanishing_checks as vc
from oddstar.vet.common import simbad_lookup
from oddstar.vet.engine import vet_target
from oddstar.vet.underluminous_checks import check_metallicity

RA, DEC = 43.548272861, 5.882712239
FAILED = "failed: ConnectionError: simulated outage"


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


def _evidence(res: Any) -> dict[str, Any]:
    assert res.evidence and Path(res.evidence).exists(), res.note
    return json.loads(Path(res.evidence).read_text(encoding="utf-8"))


def _is_pass(res: Any) -> dict[str, Any]:
    assert res.status == "pass", res.note
    return _evidence(res)


def _is_fail(res: Any) -> None:
    assert res.status == "fail", res.note
    _evidence(res)


def _is_ne(res: Any, *words: str) -> None:
    assert res.status == "not_evaluated", (res.status, res.note)
    assert res.evidence is None
    for w in words:
        assert w in res.note, (w, res.note)


# ------------------------------------------------------------------------------------------------ vanishing


def _van_star(**vd: Any) -> dict[str, Any]:
    """Hand-built fixture with every vanishing input stated explicitly (a clean, not-explained vanishing)."""
    base = {"usno_epoch": 1953.4, "usno_ra_deg": RA, "usno_dec_deg": DEC, "usno_mag_r": 16.5, "early_plates_count": 2,
            "present_in_panstarrs": False, "present_in_gaia": False, "panstarrs_depth_mag": 22.0, "sso_on_plate": False,
            "blend_split_sources": 1, "known_lpv": False}
    base.update(vd)
    return {"source_id": "755500000000000099", "synthetic_fixture": True, "ra_deg": RA, "dec_deg": DEC, "pmra_mas_yr": 2.0,
            "pmdec_mas_yr": -1.0, "simbad_type": "Star", "vanishing_data": base}


def _live(star: dict[str, Any], **qs: str) -> dict[str, Any]:
    out = copy.deepcopy(star)
    out.pop("synthetic_fixture", None)
    out["live_enriched"] = True
    out["vetting_archive"] = {"query_status": dict(qs)}
    return out


def _drop(star: dict[str, Any], key: str) -> dict[str, Any]:
    out = copy.deepcopy(star)
    out["vanishing_data"].pop(key)
    return out


VAN_CASES: list[tuple[Callable[..., Any], dict[str, Any], str]] = [
    (vc.check_plate_defect, {"early_plates_count": 1}, "early_plates_count"),
    (vc.check_depth, {"usno_mag_r": 21.0, "panstarrs_depth_mag": 21.8}, "panstarrs_depth_mag"),
    (vc.check_sso_plate, {"sso_on_plate": True}, "sso_on_plate"),
    (vc.check_blend_split, {"blend_split_sources": 3}, "blend_split_sources"),
    (vc.check_known_variable_vanishing, {"known_lpv": True}, "known_lpv"),
]


@pytest.mark.parametrize("fn,bad,key", VAN_CASES, ids=[c[0].__name__ for c in VAN_CASES])
def test_vanishing_checks_need_explicit_inputs(fn: Callable[..., Any], bad: dict[str, Any], key: str, cfg: Any, tmp_path: Path) -> None:
    ev = _is_pass(fn(_van_star(), cfg, tmp_path / "a"))
    assert ev["metrics"]["input_values"][key] is not None  # the evidence shows the value the pass rests on
    _is_fail(fn(_van_star(**bad), cfg, tmp_path / "b"))
    _is_ne(fn(_drop(_van_star(), key), cfg, tmp_path / "c"), "Missing input", key)
    # Live star: stale vanishing values are ignored while the product's query has not succeeded
    live = _live(_van_star(), simbad="ok", usno_b1=FAILED)
    _is_ne(fn(live, cfg, tmp_path / "d"), "no live ingestion")


def test_vanishing_proper_motion(cfg: Any, tmp_path: Path) -> None:
    ev = _is_pass(vc.check_proper_motion(_van_star(), cfg, tmp_path / "a"))
    assert ev["metrics"]["input_values"]["usno_epoch"] == 1953.4 and ev["metrics"]["input_values"]["pmra_mas_yr"] == 2.0
    # A 230 mas/yr star whose plate position is its PM-propagated position: explained by proper motion
    fast = _van_star(present_in_gaia=True)
    fast.update(pmra_mas_yr=185.0, pmdec_mas_yr=-140.0)
    p_ra, p_dec = propagate_position(RA, DEC, 185.0, -140.0, 2016.0, 1953.4)
    fast["vanishing_data"].update(usno_ra_deg=p_ra, usno_dec_deg=p_dec)
    _is_fail(vc.check_proper_motion(fast, cfg, tmp_path / "b"))
    for key in ("usno_epoch", "usno_ra_deg", "present_in_gaia"):
        _is_ne(vc.check_proper_motion(_drop(_van_star(), key), cfg, tmp_path / f"c_{key}"), key)
    no_pm = _van_star()
    no_pm["pmra_mas_yr"] = None  # was read as 0.0 mas/yr
    _is_ne(vc.check_proper_motion(no_pm, cfg, tmp_path / "c_pm"), "pmra_mas_yr")
    _is_ne(vc.check_proper_motion(_live(fast, usno_b1=FAILED), cfg, tmp_path / "d"), "usno_b1=failed")


def test_vanishing_known_variable_simbad_states(cfg: Any, tmp_path: Path) -> None:
    _is_fail(vc.check_known_variable_vanishing(dict(_van_star(), simbad_type="Mira"), cfg, tmp_path / "mira"))
    no_simbad = _van_star()
    no_simbad.pop("simbad_type")
    _is_ne(vc.check_known_variable_vanishing(no_simbad, cfg, tmp_path / "ns"), "SIMBAD lookup not available")
    no_entry = dict(no_simbad, vetting_archive={"query_status": {"simbad": "ok"}, "simbad_entry": False})
    ev = _is_pass(vc.check_known_variable_vanishing(no_entry, cfg, tmp_path / "ne"))
    assert "SIMBAD: no entry" in ev["note"] and ev["metrics"]["input_values"]["simbad_query_status"] == "ok"
    # Live: the catalogue flag is never fetched, so no pass; a real SIMBAD answer can still explain the source
    live = _live(_van_star(), simbad="ok")
    _is_ne(vc.check_known_variable_vanishing(live, cfg, tmp_path / "l1"), "variable_star_catalogue")
    _is_fail(vc.check_known_variable_vanishing(dict(live, simbad_type="Mira"), cfg, tmp_path / "l2"))
    _is_ne(vc.check_known_variable_vanishing(dict(_live(_van_star(), simbad=FAILED), simbad_type="Mira"), cfg, tmp_path / "l3"))


def test_vanishing_live_star_is_never_evaluated(cfg: Any, tmp_path: Path) -> None:
    # A live star as built by live ingestion: vanishing_data is empty and no plate query exists
    live = {"source_id": "1", "live_enriched": True, "ra_deg": RA, "dec_deg": DEC, "pmra_mas_yr": 37.9, "pmdec_mas_yr": -97.2,
            "simbad_type": None, "vanishing_data": {}, "vetting_archive": {"query_status": {"simbad": "ok"}, "simbad_entry": False}}
    res = vc.run_vanishing_vetting(live, cfg, tmp_path)
    assert {r.status for r in res} == {"not_evaluated"}
    assert all("no live ingestion" in r.note for r in res)


def test_vanishing_named_fixtures_still_decide(cfg: Any, tmp_path: Path) -> None:
    ok = get_named_object("755500000000000003")
    assert {r.status for r in vc.run_vanishing_vetting(ok, cfg, tmp_path / "ok")} == {"pass"}
    for sid, chk in (("755500000000000001", "proper_motion"), ("755500000000000005", "known_variable"), ("755500000000000006", "plate_defect")):
        by = {r.check: r.status for r in vc.run_vanishing_vetting(get_named_object(sid), cfg, tmp_path / sid)}
        assert by[chk] == "fail", (sid, by)


# ------------------------------------------------------------------------------------------------- spectral


def _spec_star(*, cr: bool = False, airglow: bool = False, halpha: bool = False, **over: Any) -> dict[str, Any]:
    rng = np.random.default_rng(11)
    wave = np.linspace(6500.0, 6650.0, 600)
    cont = 1.0 - 0.30 * np.exp(-0.5 * ((wave - 6562.8) / 1.8) ** 2)
    e1 = cont + rng.normal(0.0, 0.012, len(wave))
    e2 = cont + rng.normal(0.0, 0.012, len(wave))
    sky = np.full(len(wave), 0.05) + rng.normal(0.0, 0.005, len(wave))
    if cr:
        e1[260] += 1.2
    if airglow:
        oh = 0.45 * np.exp(-0.5 * ((wave - 6577.2) / 0.25) ** 2)
        e1, e2, sky = e1 + oh, e2 + oh, sky + 1.4 * oh
    if halpha:
        em = 0.65 * np.exp(-0.5 * ((wave - 6562.8) / 0.35) ** 2)
        e1, e2 = e1 + em, e2 + em
    spectra = {"wavelength_a": wave.tolist(), "flux_exp1": e1.tolist(), "flux_exp2": e2.tolist(), "sky_flux": sky.tolist(),
               "lsf_fwhm_a": 0.45, "detector_artifact": False}
    spectra.update(over)
    return {"source_id": "744400000000000099", "synthetic_fixture": True, "spectra": spectra}


def _spec_drop(star: dict[str, Any], key: str) -> dict[str, Any]:
    out = copy.deepcopy(star)
    out["spectra"].pop(key)
    return out


SPEC_CASES: list[tuple[Callable[..., Any], dict[str, Any], list[str]]] = [
    (sc.check_cosmic_ray, {"cr": True}, ["lsf_fwhm_a", "flux_exp2"]),
    (sc.check_sky_line, {"airglow": True}, ["sky_flux", "flux_exp1"]),
    (sc.check_stellar_emission, {"halpha": True}, ["flux_exp1", "wavelength_a"]),
    (sc.check_detector, {"detector_artifact": True}, ["detector_artifact"]),
    (sc.check_repeat, {"cr": True}, ["flux_exp2"]),
]


@pytest.mark.parametrize("fn,bad,missing", SPEC_CASES, ids=[c[0].__name__ for c in SPEC_CASES])
def test_spectral_checks_need_explicit_inputs(fn: Callable[..., Any], bad: dict[str, Any], missing: list[str], cfg: Any, tmp_path: Path) -> None:
    ev = _is_pass(fn(_spec_star(), cfg, tmp_path / "a"))
    assert ev["metrics"]["input_values"]
    _is_fail(fn(_spec_star(**bad), cfg, tmp_path / "b"))
    for key in missing:
        _is_ne(fn(_spec_drop(_spec_star(), key), cfg, tmp_path / f"c_{key}"), "Missing input", key)
    # Empty arrays are missing data, not a clean spectrum (sky_line used to pass on them)
    empty = _spec_star(flux_exp1=[], flux_exp2=[], wavelength_a=[], sky_flux=[], detector_artifact=None)
    _is_ne(fn(empty, cfg, tmp_path / "e"))
    live = _live(_spec_star(), spectra=FAILED)
    _is_ne(fn(live, cfg, tmp_path / "d"), "spectra=failed", "no two-epoch survey spectrum")


def test_spectral_edge_cases(cfg: Any, tmp_path: Path) -> None:
    # A cosmic-ray difference spike fails even without the LSF (a failing available input still fails)
    _is_fail(sc.check_cosmic_ray(_spec_drop(_spec_star(cr=True), "lsf_fwhm_a"), cfg, tmp_path / "cr"))
    ev = _is_pass(sc.check_cosmic_ray(_spec_star(), cfg, tmp_path / "ok"))
    assert ev["metrics"]["lsf_fwhm_a"] == 0.45 and ev["metrics"]["input_values"]["n_pixels"] == 600
    # Mismatched exposure lengths are not evaluated
    short = _spec_star()
    short["spectra"]["flux_exp2"] = short["spectra"]["flux_exp2"][:300]
    _is_ne(sc.check_repeat(short, cfg, tmp_path / "len"), "length mismatch")
    # Catalogued OH line without a sky spectrum still fails; a clean one without sky is not evaluated
    _is_fail(sc.check_sky_line(_spec_drop(_spec_star(airglow=True), "sky_flux"), cfg, tmp_path / "oh"))
    # A live star as built by live ingestion (spectra = {}) is never evaluated
    res = sc.run_spectral_vetting({"source_id": "1", "live_enriched": True, "spectra": {}, "vetting_archive": {"query_status": {}}}, cfg, tmp_path / "l")
    assert {r.status for r in res} == {"not_evaluated"} and all("spectra=not run" in r.note for r in res)


# ---------------------------------------------------------------------------------------------------- radio


def _measured(cid: str, cfg: Any, tmp: Path) -> dict[str, Any]:
    cad = copy.deepcopy(RADIO_CADENCES[cid])
    build_channel("radio_narrowband", cfg, tmp).measure(cad, {})  # writes off_scan_hits / drift_hz_s
    return cad


RADIO_CASES = [
    (rc.check_off_scan, "off_scan_hits"),
    (rc.check_zero_drift, "drift_hz_s"),
    (rc.check_rfi_band, "in_rfi_band"),
    (rc.check_multi_target, "multi_target_match"),
    (rc.check_reobservation, "reobserved_confirmed"),
]


@pytest.mark.parametrize("fn,key", RADIO_CASES, ids=[c[0].__name__ for c in RADIO_CASES])
def test_radio_checks_need_explicit_inputs(fn: Callable[..., Any], key: str, cfg: Any, tmp_path: Path) -> None:
    good = _measured("GBT_VOYAGER1_20151230", cfg, tmp_path)
    ev = _is_pass(fn(good, cfg, tmp_path / "a"))
    assert ev["metrics"]["input_values"][key] is not None
    bad = _measured("GBT_BLC1_PROXCEN", cfg, tmp_path)
    if fn is rc.check_reobservation:
        # No confirmed second epoch is an open hole that caps the tier at C, never a pass and not a fail
        _is_ne(fn(bad, cfg, tmp_path / "b"), "caps tier at C")
    elif fn is rc.check_zero_drift:
        _is_fail(fn(dict(bad, drift_hz_s=0.0), cfg, tmp_path / "b"))
    else:
        _is_fail(fn(bad, cfg, tmp_path / "b"))
    unmeasured = copy.deepcopy(good)
    unmeasured.pop(key)
    _is_ne(fn(unmeasured, cfg, tmp_path / "c"), "Missing input", key)
    live = dict(_live({"source_id": "1"}, breakthrough_listen=FAILED), **{key: good[key]})
    _is_ne(fn(live, cfg, tmp_path / "d"), "breakthrough_listen=failed", "no ON/OFF radio scans were compared")
    # A live star whose scans were compared: its inputs count, and a missing one is named as not produced by the search
    searched = dict(_live({"source_id": "1"}, breakthrough_listen="ok"), **{key: good[key]})
    _is_pass(fn(searched, cfg, tmp_path / "e"))
    _is_ne(fn(_live({"source_id": "1"}, breakthrough_listen="ok"), cfg, tmp_path / "f"), "breakthrough_listen=ok", "does not produce this input")


def test_radio_cadence_without_measurement_is_not_evaluated(cfg: Any, tmp_path: Path) -> None:
    # The raw fixture cadence carries no measured off-scan hits or drift (they used to default to 0 hits / 0.0 Hz/s)
    res = {r.check: r for r in rc.run_radio_vetting(copy.deepcopy(RADIO_CADENCES["GBT_CONTROL_000"]), cfg, tmp_path)}
    _is_ne(res["off_scan"], "channel measurement")
    _is_ne(res["zero_drift"], "channel measurement")


# -------------------------------------------------------------------------------------- SIMBAD lookup states


def test_simbad_lookup_three_states() -> None:
    assert simbad_lookup({"simbad_type": "EB*"}) == {"evaluated": True, "simbad_type": "EB*", "entry": True, "query_status": "explicit"}
    assert simbad_lookup({})["evaluated"] is False
    live_failed = {"live_enriched": True, "simbad_type": "Star", "vetting_archive": {"query_status": {"simbad": FAILED}}}
    assert simbad_lookup(live_failed) == {"evaluated": False, "simbad_type": None, "entry": None, "query_status": FAILED}
    live_none = {"live_enriched": True, "simbad_type": None, "vetting_archive": {"query_status": {"simbad": "ok"}, "simbad_entry": False}}
    assert simbad_lookup(live_none) == {"evaluated": True, "simbad_type": None, "entry": False, "query_status": "ok"}
    live_star = {"live_enriched": True, "simbad_type": "Star", "vetting_archive": {"query_status": {"simbad": "ok"}, "simbad_entry": True}}
    assert simbad_lookup(live_star)["simbad_type"] == "Star" and simbad_lookup(live_star)["entry"] is True
    assert simbad_lookup({"live_enriched": True, "simbad_type": "Star"})["query_status"] == "not run"


def _simbad_enriched(tmp: Path, *, rows: list[list[Any]] | None, fail: bool = False) -> dict[str, Any]:
    from test_ir_vetting_live import FakeArchives, _enrich, _tap

    class _Fake(FakeArchives):
        def post(self, url: str, data: dict[str, Any] | None = None, timeout: float = 0) -> Any:
            q = str((data or {}).get("QUERY", ""))
            if rows is not None and self._service_post(url, q) == "simbad":
                self.calls.append(("simbad", {"url": url, "QUERY": q}))
                return _tap(["main_id", "otype_txt", "sep_arcsec"], rows)
            return super().post(url, data, timeout)

    return _enrich(tmp, _Fake(fail={"simbad"} if fail else None))


def test_live_enrichment_keeps_simbad_states_distinct(tmp_path: Path, cfg: Any) -> None:
    typed = _simbad_enriched(tmp_path / "t", rows=[["V* XX Tst", "EB*", 0.2]])
    empty = _simbad_enriched(tmp_path / "e", rows=[])
    failed = _simbad_enriched(tmp_path / "f", rows=None, fail=True)
    assert typed["simbad_type"] == "EB*" and typed["vetting_archive"]["query_status"]["simbad"] == "ok"
    assert typed["vetting_archive"]["simbad_entry"] is True
    assert empty["simbad_type"] is None and empty["vetting_archive"]["query_status"]["simbad"] == "ok"
    assert empty["vetting_archive"]["simbad_entry"] is False
    assert failed["simbad_type"] is None and failed["vetting_archive"]["query_status"]["simbad"].startswith("failed: ")
    assert "simbad_entry" not in failed["vetting_archive"]
    # Checks reading simbad_type: EB* fails, "no entry" is a real negative answer, a failed lookup is not evaluated
    _is_fail(lcc.check_known_variable(typed, cfg, tmp_path / "kt"))
    ev = _is_pass(lcc.check_known_variable(empty, cfg, tmp_path / "ke"))
    assert ev["metrics"]["simbad_query_status"] == "ok" and ev["metrics"]["simbad_entry"] is False
    _is_ne(lcc.check_known_variable(failed, cfg, tmp_path / "kf"), "simbad=failed")
    _is_pass(check_metallicity(empty, cfg, tmp_path / "me"))
    _is_ne(check_metallicity(failed, cfg, tmp_path / "mf"), "SIMBAD lookup not available")


def test_known_variable_requires_simbad(cfg: Any, tmp_path: Path) -> None:
    rng = np.random.default_rng(4)
    t = np.concatenate([np.arange(2144.5, 2156.5, 1 / 48), np.arange(2158.0, 2170.0, 1 / 48)])
    quiet = 1.0 + rng.normal(0, 0.001, len(t))
    sec = {"sector": 31, "time_btjd": t.tolist(), "flux_qlp_raw": quiet.tolist(), "flux_tesscut": quiet.tolist()}
    eb = quiet.copy()
    eb[np.abs((((t - 2145.0) % 2.1) / 2.1) - 0.5) < 0.02] -= 0.25
    eb_sec = dict(sec, flux_qlp_raw=eb.tolist(), flux_tesscut=eb.tolist())
    hand = {"source_id": "1", "simbad_type": "Star", "lightcurves": {"tess_sectors": [sec]}}
    ev = _is_pass(lcc.check_known_variable(hand, cfg, tmp_path / "a"))
    assert ev["metrics"]["simbad_type"] == "Star" and ev["metrics"]["simbad_query_status"] == "explicit"
    _is_fail(lcc.check_known_variable(dict(hand, simbad_type="RR*"), cfg, tmp_path / "b"))
    # Quiet TESS light curve, no SIMBAD answer: used to pass as SIMBAD type "Star"
    _is_ne(lcc.check_known_variable({"source_id": "1", "lightcurves": {"tess_sectors": [sec]}}, cfg, tmp_path / "c"), "SIMBAD lookup not available")
    live = _live(hand, simbad=FAILED)
    _is_ne(lcc.check_known_variable(live, cfg, tmp_path / "d"), "simbad=failed")
    # The TESS period search still fails a periodic eclipsing binary without SIMBAD
    _is_fail(lcc.check_known_variable(dict(live, lightcurves={"tess_sectors": [eb_sec]}), cfg, tmp_path / "e"))


def test_metallicity_inputs(cfg: Any, tmp_path: Path) -> None:
    ev = _is_pass(check_metallicity({"feh": -0.1, "simbad_type": "Star"}, cfg, tmp_path / "a"))
    assert ev["metrics"]["feh"] == -0.1 and ev["metrics"]["feh_subdwarf_max"] == -1.0
    _is_fail(check_metallicity({"feh": -1.8, "simbad_type": "Star"}, cfg, tmp_path / "b"))
    _is_fail(check_metallicity({"feh": -0.2, "simbad_type": "Subdwarf"}, cfg, tmp_path / "b2"))
    _is_ne(check_metallicity({"feh": -0.1}, cfg, tmp_path / "c"), "SIMBAD")
    _is_ne(check_metallicity({"simbad_type": "Star"}, cfg, tmp_path / "c2"), "feh")
    live = {"live_enriched": True, "feh": -0.1, "simbad_type": "Star", "vetting_archive": {"query_status": {"simbad": FAILED}}}
    _is_ne(check_metallicity(live, cfg, tmp_path / "d"), "simbad=failed")


# --------------------------------------------------------------------------------- generic_outlier path (IR)


def _ir_star(**over: Any) -> dict[str, Any]:
    w_ra, w_dec = propagate_position(RA, DEC, 37.9, -97.2, 2016.0, 2010.2)
    star = {"source_id": "1", "ra_deg": RA, "dec_deg": DEC, "pmra_mas_yr": 37.9, "pmdec_mas_yr": -97.2, "abs_g_mag": 6.1, "teff_k": 4500.0,
            "teff_k_source": "gspphot", "ruwe": 1.02, "feh": -0.1, "simbad_type": "Star",
            "allwise": {"ra_deg": w_ra, "dec_deg": w_dec, "n_neighbours_gaia_xmatch": 1},
            "vetting_archive": {"nss_flag": 0, "rv_scatter_kms": 0.4, "radio_sources_within_10arcsec": [],
                                "second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}}
    star.update(over)
    return star


def test_generic_path_ir_checks_need_inputs(cfg: Any, tmp_path: Path) -> None:
    ev = _is_pass(irc.check_xmatch_integrity(_ir_star(), cfg, tmp_path / "x1"))
    assert ev["metrics"]["pmra_mas_yr"] == 37.9
    _is_fail(irc.check_xmatch_integrity(_ir_star(allwise={**_ir_star()["allwise"], "n_neighbours_gaia_xmatch": 2}), cfg, tmp_path / "x2"))
    _is_ne(irc.check_xmatch_integrity(_ir_star(pmra_mas_yr=None), cfg, tmp_path / "x3"), "proper motion")  # was read as 0.0

    _is_pass(irc.check_evolved(_ir_star(), cfg, tmp_path / "e1"))
    _is_fail(irc.check_evolved(_ir_star(abs_g_mag=0.0, teff_k=4600.0), cfg, tmp_path / "e2"))
    _is_ne(irc.check_evolved(_ir_star(abs_g_mag=None), cfg, tmp_path / "e3"), "Absolute G")  # was read as M_G = 4.8
    _is_ne(irc.check_evolved(_ir_star(abs_g_mag=0.0, teff_k=None), cfg, tmp_path / "e4"), "Teff")

    _is_pass(irc.check_binary(_ir_star(), cfg, tmp_path / "b1"))
    _is_fail(irc.check_binary(_ir_star(ruwe=2.4), cfg, tmp_path / "b2"))
    _is_ne(irc.check_binary(_ir_star(ruwe=None), cfg, tmp_path / "b3"), "RUWE not available")  # was read as 1.0
    live = _live(_ir_star(), gaia_meta=FAILED)
    live["vetting_archive"].update(nss_flag=0, rv_scatter_kms=0.4)
    _is_ne(irc.check_binary(live, cfg, tmp_path / "b4"), "gaia_meta=failed")

    _is_pass(irc.check_second_instrument(_ir_star(), cfg, tmp_path / "s1"))
    bad = _ir_star()
    bad["vetting_archive"]["second_instrument_ir"] = {"instrument": "AKARI", "interpretation": "no_excess_at_second_instrument"}
    _is_fail(irc.check_second_instrument(bad, cfg, tmp_path / "s2"))
    bare = _ir_star()
    bare["vetting_archive"]["second_instrument_ir"] = {"instrument": "AKARI"}  # was read as "confirmed"
    _is_ne(irc.check_second_instrument(bare, cfg, tmp_path / "s3"), "no comparison outcome")

    _is_pass(irc.check_radio_counterpart(_ir_star(), cfg, tmp_path / "r1"))
    rad = _ir_star()
    rad["vetting_archive"]["radio_sources_within_10arcsec"] = [{"survey": "NVSS", "sep_arcsec": 2.0}]
    _is_fail(irc.check_radio_counterpart(rad, cfg, tmp_path / "r2"))
    rad["vetting_archive"]["radio_sources_within_10arcsec"] = [{"survey": "NVSS"}]  # was read as 99"
    _is_ne(irc.check_radio_counterpart(rad, cfg, tmp_path / "r3"), "without a separation")

    fx = generate_field_sample(n_stars=2, seed=3)[0]
    assert irc.check_vision_rubric(fx, cfg, tmp_path / "v1").status == "pass"
    fx_bad = copy.deepcopy(fx)
    fx_bad["allwise"]["centroid_offset_arcsec"] = 6.0
    assert irc.check_vision_rubric(fx_bad, cfg, tmp_path / "v2").status == "fail"
    fx_bad["allwise"].pop("centroid_offset_arcsec")  # was read as 0.0"
    _is_ne(irc.check_vision_rubric(fx_bad, cfg, tmp_path / "v3"), "explicit centroid offset")


def test_generic_outlier_dispatch_does_not_pass_without_inputs(cfg: Any, tmp_path: Path) -> None:
    bare = {"source_id": "1", "ra_deg": RA, "dec_deg": DEC, "lightcurves": {}, "allwise": {"ra_deg": RA, "dec_deg": DEC, "n_neighbours_gaia_xmatch": 1}}
    score = ChannelScore(source_id="1", channel="generic_outlier", raw_statistic=5.0, false_alarm_prob=1e-4, flagged=True, rank_score=4.0,
                         summary="t", feature_contributions={"delta_m_sigma": 3.0})
    res = {r.check: r for r in vet_target(bare, [score], cfg, tmp_path)}
    for chk in ("xmatch_integrity", "known_variable", "metallicity", "evolved"):
        _is_ne(res[chk])
    assert not any(r.status == "pass" for r in res.values()), {k: v.status for k, v in res.items() if v.status == "pass"}
    # A live-style light curve (per-sector TESS data only) now routes a dip-driven outlier to the light-curve checks
    rng = np.random.default_rng(23)
    lc_star = dict(bare, simbad_type="Star", lightcurves={"tess_sectors": [_sector(rng, 2150.5, 0.75)]})
    dip = score.model_copy(update={"feature_contributions": {"dip_snr": 9.0}})
    checks = {r.check for r in vet_target(lc_star, [dip], cfg, tmp_path / "lc")}
    assert {"pipeline_agreement", "systematics_window", "neighbour_signal", "sso_crossing"} <= checks


# ------------------------------------------------------------------------------------------- edge dips (LC)


def _sector(rng: np.random.Generator, dip_t: float | None, dip_len: float = 0.5, depth: float = 0.06) -> dict[str, Any]:
    """Two-orbit 30-min sector (orbits 2144.5-2156.48 and 2158.0-2169.98) with background and centroid series."""
    from test_lc_vetting_live import tess_sector

    if dip_t is None:
        return tess_sector(rng)
    return tess_sector(rng, dip=depth, dip_t=dip_t, dip_len=dip_len)


def test_dip_at_segment_edge_is_not_evaluated(cfg: Any, tmp_path: Path) -> None:
    rng = np.random.default_rng(21)
    zone = cfg.channels["optical_dips"]["detrend_edge_zone_window_frac"] * cfg.channels["optical_dips"]["detrend_window_days"]
    assert zone == 2.5
    # A 6% dip over the first 0.5 d of the sector: the running median erases it ...
    edge = _sector(rng, 2144.75)
    t = np.asarray(edge["time_btjd"])
    det = lcc.matched_detrend(t, np.asarray(edge["flux_qlp_raw"]), 5.0, 0.5, 2.5)
    in_dip = np.abs(t - 2144.75) <= 0.2
    assert 1.0 - float(np.median(det[in_dip])) < 0.02  # 6% dip reduced below the 3% detection threshold
    a = lcc.analyze_tess_sector(edge, cfg)
    assert not a["summary"]["deepest_coincident"]["coincident_exceeds_min_depth"]
    cands = a["summary"]["detrend_edge_zone"]["candidates"]
    assert cands and min(c["edge_distance_days"] for c in cands) < 0.5
    assert max(c["depth_qlp_vs_interior"] for c in cands) > 0.04  # the interior reference still sees it
    # ... so no detrend-based check may pass on "no dip found" (each used to pass)
    star = {"source_id": "1", "ra_deg": RA, "dec_deg": DEC, "mags": {"G": 13.4}, "simbad_type": "Star", "neowise": {},
            "lightcurves": {"tess_sectors": [edge]}}
    for fn in (lcc.check_pipeline_agreement, lcc.check_systematics_window, lcc.check_neighbour_signal, lcc.check_sso_crossing):
        r = fn(star, cfg, tmp_path / fn.__name__)
        _is_ne(r, "detrend edge zone")
        assert json.dumps(r.metrics)
    # A dip 1.1 d before the inter-orbit data gap (outside the 0.5 d systematics margin) is in the zone too
    gap = dict(star, lightcurves={"tess_sectors": [_sector(rng, 2155.35, 0.75)]})
    for fn in (lcc.check_pipeline_agreement, lcc.check_systematics_window, lcc.check_neighbour_signal):
        _is_ne(fn(gap, cfg, tmp_path / f"g_{fn.__name__}"), "segment edge or data gap")
    # The same dip mid-orbit is evaluated and passes
    mid = dict(star, lightcurves={"tess_sectors": [_sector(rng, 2150.5, 0.75)]})
    for fn in (lcc.check_pipeline_agreement, lcc.check_systematics_window, lcc.check_neighbour_signal):
        ev = _is_pass(fn(mid, cfg, tmp_path / f"m_{fn.__name__}"))
        assert "detrend_edge_zone_candidates" in json.dumps(ev["metrics"])
    # A quiet sector has no edge-zone candidate and still passes
    _is_pass(lcc.check_sso_crossing(dict(star, lightcurves={"tess_sectors": [_sector(rng, None)]}), cfg, tmp_path / "q"))
    # The zone width is read from config
    narrow = cfg.model_copy(deep=True)
    narrow.channels["optical_dips"]["detrend_edge_zone_window_frac"] = 0.1
    _is_pass(lcc.check_pipeline_agreement(gap, narrow, tmp_path / "narrow"))


def test_edge_dip_blocks_second_survey_pass(cfg: Any, tmp_path: Path) -> None:
    rng = np.random.default_rng(22)

    def series(dip_t: float) -> list[dict[str, Any]]:
        # Nightly ground-survey medians every 0.25 d across the sector, 6% fainter during the TESS dip
        nights = [[t + lcc.BTJD_TO_MJD, 0.067 if abs(t - dip_t) <= 0.2 else 0.0] for t in np.arange(2144.5, 2170.0, 0.25)]
        return [{"key": "ztf", "survey": "ZTF", "band": "zg", "kind": "calendar", "seasons": [], "nights": nights, "night_sigma": 0.01}]

    mid = lcc._dip_context({"tess_sectors": [_sector(rng, 2150.5, 0.75)]}, series(2150.5), cfg)
    assert mid is not None and mid["status"] == "pass", mid
    edge = lcc._dip_context({"tess_sectors": [_sector(rng, 2144.75)]}, series(2144.75), cfg)
    assert edge is not None and edge["status"] == "not_evaluated" and "detrend edge zone" in edge["note"]  # was None (no dip)
