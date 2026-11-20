"""The dip the search flagged must be in the light curves the checks read, or they cannot pass.

Regression for Gaia DR3 23767631757133952 (nearby 1,300-star run, config 1.12.0): the search flagged a 26% box in
TESS sector 42, where the star lies on two bad CCD columns. The second extraction returned no cadence for that
sector, so vetting loaded sectors 43, 70 and 71 only, found them quiet, and four checks passed. ZTF observed the
star inside the box and saw no dimming. Offline: every light curve here is synthetic.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oddstar.config import load_config
from oddstar.vet.lc_checks import (
    BTJD_TO_MJD,
    SEARCH_DIP_MIN_BOX_COVERAGE,
    check_neighbour_signal,
    check_pipeline_agreement,
    check_second_survey,
    check_sso_crossing,
    check_systematics_window,
    search_dip_from_measurement,
    search_dip_gap,
)
from test_lc_vetting_live import live_star, seasons_from, tess_sector

# The search's box of the regression star: sector 42, 48 h wide
SEARCH_DIP = {"sector": 42, "btjd": 2469.93118, "depth": 0.2615, "duration_hours": 48.0, "snr": 10.09}
FETCH_LOG = ["s42: only 0 good cadences", "s43: ok (2550 matched cadences, 964 30-min bins)"]


@pytest.fixture(scope="module")
def cfg() -> Any:
    c, _, _ = load_config()
    return c


def _quiet_sectors(seed: int = 1) -> list[dict[str, Any]]:
    rng = np.random.default_rng(seed)
    return [tess_sector(rng, sector=43), tess_sector(rng, sector=70)]


def _ztf(nights: list[list[float]]) -> dict[str, Any]:
    band = {"n_epochs": len(nights), "seasons": seasons_from({2021: 13.63}), "nights": sorted(nights), "nightly_robust_sigma_mag": 0.014}
    return {"ztf": {"survey": "ZTF", "status": "ok", "season_kind": "calendar", "n_epochs": len(nights), "bands": {"zg": band}}}


def test_measurement_gives_the_search_dip_only_for_a_flagged_box_in_a_real_sector(cfg: Any) -> None:
    vals = {"dip_snr": 10.09, "dip_depth": 0.2615, "duration_hours": 48.0, "periodic_bls_snr": 3.82,
            "best_time_btjd": 2469.93118, "best_sector": 42.0}
    sd = search_dip_from_measurement(vals, cfg)
    assert sd == SEARCH_DIP
    # below the channel thresholds, periodic, a fixture sector, or an old measurement without the box time: nothing to hand over
    assert search_dip_from_measurement(dict(vals, dip_snr=5.0), cfg) is None
    assert search_dip_from_measurement(dict(vals, dip_depth=0.01), cfg) is None
    assert search_dip_from_measurement(dict(vals, periodic_bls_snr=9.0), cfg) is None
    assert search_dip_from_measurement(dict(vals, best_sector=0.0), cfg) is None
    assert search_dip_from_measurement({k: v for k, v in vals.items() if k != "best_time_btjd"}, cfg) is None


def test_four_tess_checks_cannot_pass_on_sectors_that_do_not_hold_the_flagged_dip(tmp_path: Path, cfg: Any) -> None:
    secs = _quiet_sectors()
    checks = (check_pipeline_agreement, check_systematics_window, check_neighbour_signal, check_sso_crossing)
    # Without the search's dip the quiet sectors pass, as before (a star flagged on another channel)
    for fn in checks:
        assert fn(live_star(tess_sectors=secs), cfg, tmp_path / "plain").status == "pass", fn.__name__
    star = live_star(tess_sectors=secs, search_dip=SEARCH_DIP, tess_query_log=FETCH_LOG)
    gap = search_dip_gap(star["lightcurves"], cfg)
    assert gap is not None and gap["vetted_sectors"] == [43, 70] and gap["box_coverage"] == 0.0
    assert gap["fetch_log"] == ["s42: only 0 good cadences"]
    for fn in checks:
        res = fn(star, cfg, tmp_path / "gap")
        assert res.status == "not_evaluated", fn.__name__
        assert "sector 42" in res.note and "not among the vetted sectors [43, 70]" in res.note, fn.__name__
        assert "only 0 good cadences" in res.note, fn.__name__
        assert res.metrics["search_dip_gap"]["search_dip"] == SEARCH_DIP, fn.__name__


def test_no_vetting_sector_at_all_is_reported_with_the_flagged_dip(tmp_path: Path, cfg: Any) -> None:
    star = live_star(search_dip=SEARCH_DIP, tess_query_log=FETCH_LOG, tess_status="no_sector_usable")
    for fn in (check_pipeline_agreement, check_systematics_window, check_neighbour_signal, check_sso_crossing):
        res = fn(star, cfg, tmp_path / "none")
        assert res.status == "not_evaluated" and "the dip the search flagged (sector 42" in res.note, fn.__name__


def test_a_failing_check_still_fails_when_the_flagged_dip_is_in_another_sector(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(3)
    # sector 43 has a dip in one extraction only: the disagreement is a fail whatever happened in sector 42
    bad = tess_sector(rng, sector=43, dip=0.08, cut_dip_t=2162.0)
    star = live_star(tess_sectors=[bad], search_dip=SEARCH_DIP)
    assert check_pipeline_agreement(star, cfg, tmp_path / "f").status == "fail"


def test_a_dip_inside_the_vetted_sector_leaves_the_checks_as_they_were(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(4)
    sec = tess_sector(rng, sector=31, dip=0.06)
    sd = {"sector": 31, "btjd": 2150.5, "depth": 0.06, "duration_hours": 18.0, "snr": 30.0}
    star = live_star(tess_sectors=[sec], search_dip=sd)
    assert search_dip_gap(star["lightcurves"], cfg) is None
    res = check_pipeline_agreement(star, cfg, tmp_path / "in")
    assert res.status == "pass" and res.metrics["search_dip_gap"] is None
    # the same sector, but the box lies in the orbit gap where no cadence was matched
    in_gap = dict(sd, btjd=2157.25, duration_hours=24.0)
    star_gap = live_star(tess_sectors=[sec], search_dip=in_gap)
    gap = search_dip_gap(star_gap["lightcurves"], cfg)
    assert gap is not None and gap["box_coverage"] < SEARCH_DIP_MIN_BOX_COVERAGE and "of the box lies in vetted bins" in gap["note"]
    assert check_pipeline_agreement(star_gap, cfg, tmp_path / "og").status == "not_evaluated"


def test_second_survey_tests_the_search_box_against_ground_nights(tmp_path: Path, cfg: Any) -> None:
    rng = np.random.default_rng(5)
    ev_mjd = SEARCH_DIP["btjd"] + BTJD_TO_MJD
    base = [[ev_mjd - 60 + k, 13.631 + float(rng.normal(0, 0.014))] for k in range(0, 120, 2) if abs(k - 60) > 3]
    secs = _quiet_sectors(6)

    def star_with(night_mag: float | None) -> dict[str, Any]:
        nights = list(base) + ([[ev_mjd, night_mag]] if night_mag is not None else [])
        return live_star(tess_sectors=secs, search_dip=SEARCH_DIP, optical_surveys=_ztf(nights))

    # ZTF saw the star at the centre of the box at its normal brightness: a 26% dip would be 0.33 mag
    res_c = check_second_survey(star_with(13.633), cfg, tmp_path / "c")
    assert res_c.status == "fail" and "not on the sky" in res_c.note and "the search's box" in res_c.note
    assert res_c.metrics["contexts"]["optical_dips"]["search_dip_gap"]["search_dip"] == SEARCH_DIP
    # ZTF saw it 0.33 mag fainter: the dimming is on the sky
    res_p = check_second_survey(star_with(13.961), cfg, tmp_path / "p")
    assert res_p.status == "pass" and "consistent with TESS" in res_p.note
    # no night inside the box: the surveys cannot say, and the note names the dip instead of "no flagged signal"
    res_n = check_second_survey(star_with(None), cfg, tmp_path / "n")
    assert res_n.status == "not_evaluated" and "the search's box" in res_n.note and "No flagged signal" not in res_n.note
    # a night just outside the box is not used for it (the box is tested without padding)
    nights = list(base) + [[ev_mjd + 1.1, 13.633]]
    res_o = check_second_survey(live_star(tess_sectors=secs, search_dip=SEARCH_DIP, optical_surveys=_ztf(nights)), cfg, tmp_path / "o")
    assert res_o.status == "not_evaluated"
    # the same nights without the search's dip: nothing flagged in the vetted sectors, as before
    res_plain = check_second_survey(live_star(tess_sectors=secs, optical_surveys=_ztf(list(base) + [[ev_mjd, 13.633]])), cfg, tmp_path / "x")
    assert res_plain.status == "not_evaluated" and "No flagged signal to test" in res_plain.note


def test_vet_hands_the_search_dip_to_the_checks(tmp_path: Path, cfg: Any) -> None:
    from oddstar.cli import _attach_search_dip

    star = live_star(tess_sectors=_quiet_sectors(7))
    _attach_search_dip(star, {"dip_snr": 10.09, "dip_depth": 0.2615, "duration_hours": 48.0, "periodic_bls_snr": 3.82,
                              "best_time_btjd": 2469.93118, "best_sector": 42.0}, cfg)
    assert star["lightcurves"]["search_dip"] == SEARCH_DIP
    assert check_pipeline_agreement(star, cfg, tmp_path / "v").status == "not_evaluated"
    quiet = live_star(tess_sectors=_quiet_sectors(8))
    _attach_search_dip(quiet, {"dip_snr": 2.0, "dip_depth": 0.001, "duration_hours": 6.0, "periodic_bls_snr": 1.0,
                               "best_time_btjd": 2469.9, "best_sector": 42.0, "secular_fade_snr": 6.0}, cfg)
    assert "search_dip" not in quiet["lightcurves"]
