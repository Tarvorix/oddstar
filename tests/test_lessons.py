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
