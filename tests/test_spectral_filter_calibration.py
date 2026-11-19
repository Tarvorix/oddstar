"""spectral_lines: the matched-filter output is judged against its own baseline and scatter (config 1.12.0).

A low-resolution spectrum of an F, G or K star has absorption lines the template does not hold. The stretch of
continuum between two of them is a bump as wide as the line spread, the same in every exposure. The stars here are
built that way: no emission line unless the test adds one.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from oddstar.channels import build_channel
from oddstar.channels.spectral_lines import channel as spectral_channel
from oddstar.channels.spectral_lines.channel import MIN_OPEN_PIXELS_FOR_CALIBRATION, _calibrated_filter_snr
from oddstar.config import load_config, repo_root

# Stronger photospheric lines of the window and just outside it (air, A; solar equivalent width, mA)
METAL_LINES = [(6493.78, 130), (6494.98, 165), (6496.47, 65), (6496.90, 98), (6498.94, 46), (6499.65, 86), (6516.08, 54),
               (6518.37, 57), (6546.24, 105), (6569.22, 73), (6572.78, 34), (6574.23, 28), (6575.02, 65), (6586.31, 40),
               (6592.91, 125), (6593.87, 85), (6597.56, 44), (6609.11, 65), (6633.75, 70), (6643.63, 93), (6663.44, 80)]
SKY_REST = [6498.73, 6533.04, 6544.03, 6553.62, 6562.79, 6568.78, 6577.20, 6596.60, 6604.10]
WAVE = np.arange(6500.0, 6650.0 + 0.1, 1.5)
LSF_FWHM = 3.65
LSF_SIGMA = LSF_FWHM / 2.355


@pytest.fixture(scope="module")
def channel() -> Any:
    return build_channel("spectral_lines", load_config()[0], repo_root())


def _star(seed: int, *, snr: tuple[float, float] = (250.0, 175.0), line_scale: float = 1.0, tilt: float = 0.0,
          metal_lines: bool = True, emission: tuple[float, float] | None = None) -> dict[str, Any]:
    rng = np.random.default_rng(seed)
    flux = np.ones_like(WAVE)
    if metal_lines:
        for centre, ew_ma in METAL_LINES:
            flux -= line_scale * (ew_ma / 1000.0) / (LSF_SIGMA * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((WAVE - centre) / LSF_SIGMA) ** 2)
        flux -= 0.25 * np.exp(-0.5 * ((WAVE - 6562.8) / 2.4) ** 2) + 0.03 / (1 + ((WAVE - 6562.8) / 8.0) ** 2)
    flux *= 1 + tilt * (WAVE - 6575.0) / 75.0
    if emission:
        flux += emission[1] * np.exp(-0.5 * ((WAVE - emission[0]) / LSF_SIGMA) ** 2)
    return {"source_id": str(seed), "spectra": {
        "wavelength_a": WAVE.tolist(), "flux_exp1": (flux + rng.normal(0, 1 / snr[0], len(WAVE))).tolist(),
        "flux_exp2": (flux + rng.normal(0, 1 / snr[1], len(WAVE))).tolist(), "lsf_fwhm_a": LSF_FWHM,
        "sky_line_wavelengths_a": SKY_REST}}


@pytest.mark.parametrize("line_scale", [1.0, 1.8, 4.0])
def test_continuum_between_unmodelled_absorption_lines_is_not_a_line(channel: Any, line_scale: float) -> None:
    for seed in range(20):
        m = channel.measure(_star(seed, line_scale=line_scale), {})
        assert m.values["max_residual_snr"] < 3.5, (seed, m.values)
        assert not channel.score(m).flagged
        for cal in m.extra["filter_calibration"]:
            assert cal["pixel_noise"] > 0 and cal["filter_scatter"] > 0
            assert cal["n_open_pixels"] >= MIN_OPEN_PIXELS_FOR_CALIBRATION


def test_without_the_absorption_templates_the_baseline_carries_the_offset(channel: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Continuum and H-alpha alone sit below the continuum between the lines; the baseline of the filter output takes that up."""
    monkeypatch.setattr(spectral_channel, "PHOTOSPHERIC_LINES_A", ())
    for seed in range(20):
        m = channel.measure(_star(seed, line_scale=1.8), {})
        assert m.values["max_residual_snr"] < 3.5, (seed, m.values)
        for cal in m.extra["filter_calibration"]:
            assert cal["baseline"] > 0.8 * cal["pixel_noise"]  # an offset of about one pixel noise under every peak


def test_a_continuum_that_falls_to_the_red_is_fitted(channel: Any) -> None:
    for seed in range(20):
        m = channel.measure(_star(seed, tilt=-0.06, line_scale=0.5), {})
        assert m.values["max_residual_snr"] < 3.5, (seed, m.values)


def test_a_line_well_above_the_structure_still_flags(channel: Any) -> None:
    for where in (6508.0, 6620.0):
        for seed in range(20):
            m = channel.measure(_star(seed, emission=(where, 0.20)), {})
            assert channel.score(m).flagged, (where, seed, m.values)
            assert abs(m.values["peak_wavelength_a"] - where) <= 1.5
            assert m.values["repeat_ratio"] > 0.5


def test_a_line_in_white_noise_still_flags_and_noise_alone_does_not(channel: Any) -> None:
    hits = sum(channel.score(channel.measure(_star(s, snr=(60.0, 60.0), metal_lines=False, emission=(6620.0, 0.16)), {})).flagged for s in range(40))
    assert hits >= 38
    quiet = [channel.measure(_star(s, snr=(60.0, 60.0), metal_lines=False), {}).values["max_residual_snr"] for s in range(40)]
    assert max(quiet) < 3.5
    # the two exposures are compared with a template fitted to their minimum: that offset is in the baseline, not the peak
    assert float(np.median(quiet)) < 1.6


def test_calibration_lowers_a_significance_and_never_raises_it() -> None:
    rng = np.random.default_rng(5)
    open_window = np.ones(60, dtype=bool)
    quiet = rng.normal(0.0, 0.002, 60)  # scatters by less than the pixel noise
    snr, rec = _calibrated_filter_snr(quiet, 0.01, open_window)
    assert np.allclose(snr, (quiet - np.median(quiet)) / 0.01)
    assert rec["filter_scatter"] < rec["pixel_noise"] == 0.01
    rough = rng.normal(0.03, 0.05, 60)  # an offset and more scatter than the pixel noise
    snr, rec = _calibrated_filter_snr(rough, 0.01, open_window)
    assert 0.0 < rec["baseline"] < 0.06 and rec["filter_scatter"] > 0.03
    assert np.all(np.abs(snr) <= np.abs(rough - np.median(rough)) / 0.01)
    # too few open pixels: no baseline, the pixel noise alone
    few = np.zeros(60, dtype=bool)
    few[: MIN_OPEN_PIXELS_FOR_CALIBRATION - 1] = True
    snr, rec = _calibrated_filter_snr(rough, 0.01, few)
    assert rec["baseline"] == 0.0 and rec["filter_scatter"] == 0.0
    assert np.allclose(snr, rough / 0.01)
