"""The point-source FWHM calibration tool: selection, statistics, the pre-registered rule, and failure handling."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from oddstar.config import load_config
from oddstar.io import cutouts
from oddstar.io import live_survey as ls
from oddstar.io.cache import RawCache
from oddstar.validate import atlas_psf as ap


def test_selection_query_is_fixed_and_candidate_free() -> None:
    q = ap.selection_query(150.0, 2.2)
    for cut in ("ext_flg = 0", "cc_flags = '0000'", "nb = 1", "na = 0", "w3sat = 0", "w4sat = 0", "w3snr >= 20", "w4snr >= 10",
                "w3rchi2 < 2", "w4rchi2 < 2", "w3mpro - w4mpro BETWEEN -0.2 AND 0.3"):
        assert cut in q
    assert "CIRCLE('ICRS', 150.0000, 2.2000, 1.0)" in q
    # the reference star of the live runs (RA 43.55, Dec +5.88) is in none of the calibration fields
    assert all(abs(ra - 43.548) > 5.0 or abs(dec - 5.883) > 5.0 for ra, dec in ap.CALIBRATION_FIELDS)
    from astropy import units as u
    from astropy.coordinates import SkyCoord

    assert all(abs(SkyCoord(ra * u.deg, dec * u.deg).galactic.b.deg) > 30.0 for ra, dec in ap.CALIBRATION_FIELDS)


def test_band_stats_and_rule() -> None:
    st = ap.band_stats([10.0, 11.0, 12.0, 13.0, 14.0])
    assert st["n"] == 5 and st["median_arcsec"] == 12.0
    assert st["robust_sigma_arcsec"] == pytest.approx(1.4826)
    assert ap.band_stats([]) == {"n": 0}

    tight = {"n": 20, "median_arcsec": 12.0, "robust_sigma_arcsec": 0.3}
    wide = {"n": 20, "median_arcsec": 20.0, "robust_sigma_arcsec": 2.1}
    assert ap.proposed_ratio({"W3": tight, "W4": tight}, 1.3) == pytest.approx(1.3)  # never below the current value
    assert ap.proposed_ratio({"W3": tight, "W4": wide}, 1.3) == pytest.approx(1.45)  # 1 + 4 * 2.1 / 20 = 1.42, rounded up
    assert ap.proposed_ratio({"W3": tight, "W4": {"n": 14, "median_arcsec": 20.0, "robust_sigma_arcsec": 1.0}}, 1.3) is None


def test_run_calibration_uses_the_pipeline_measurement_and_records_failures(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    queries: list[str] = []

    def fake_tap(*, url: str, query: str, archive: str, release: str, cache: RawCache, timeout_s: int) -> Any:
        queries.append(query)
        if len(queries) == 2:
            raise RuntimeError("IRSA busy")
        k = len(queries)
        rows = [{"designation": f"J{k:02d}{i:02d}", "ra": 10.0 * k + 0.01 * i, "dec": 1.0, "w3mpro": 5.0, "w3snr": 40.0, "w3rchi2": 1.0,
                 "w4mpro": 4.95, "w4snr": 15.0, "w4rchi2": 1.0} for i in (3, 1, 2)]
        return rows, None, False

    fetched: list[tuple[float, float]] = []

    def fake_fetch(*, ra: float, dec: float, size_arcsec: float, cache: RawCache, sleep: Any) -> Any:
        fetched.append((ra, dec))
        if len(fetched) == 1:
            raise RuntimeError("cutout service down")
        return {"W3": {}, "W4": {}}, [], []

    def fake_measure(images: Any, *, ra_pred: float, dec_pred: float, vetting: Any) -> Any:
        n = len(fetched)
        return {"W3": {"status": "ok", "fwhm_arcsec": 12.0 + 0.01 * n, "fwhm_ref_arcsec": 9.19, "fwhm_ratio": 1.3, "flux_snr": 30.0},
                "W4": {"status": "ok", "fwhm_arcsec": 21.0 + 0.01 * n, "fwhm_ref_arcsec": 16.97, "fwhm_ratio": 1.24, "flux_snr": 2.0 if n == 2 else 12.0}}

    monkeypatch.setattr(ls, "_tap_query_cached", fake_tap)
    monkeypatch.setattr(cutouts, "fetch_allwise_cutouts", fake_fetch)
    monkeypatch.setattr(cutouts, "measure_wise_cutouts", fake_measure)

    rec = ap.run_calibration(cfg, cache, max_per_field=2, sleep=lambda s: None)

    assert len(queries) == len(ap.CALIBRATION_FIELDS)
    failed = [f for f in rec["fields"] if f["status"] != "ok"]
    assert len(failed) == 1 and "IRSA busy" in failed[0]["status"]
    n_ok_fields = len(ap.CALIBRATION_FIELDS) - 1
    assert rec["n_stars"] == 2 * n_ok_fields
    # the first two rows by designation are used in each field
    assert [s["designation"] for s in rec["stars"][:2]] == ["J0101", "J0102"]
    assert "cutout service down" in rec["stars"][0]["status"]
    # one star failed to fetch; one more has a W4 fit below the minimum SNR and is left out of W4 only
    assert rec["stats"]["W3"]["n"] == rec["n_stars"] - 1
    assert rec["stats"]["W4"]["n"] == rec["n_stars"] - 2
    assert rec["stats"]["W3"]["theoretical_reference_arcsec"] == pytest.approx(6.5 * 1.414)
    assert rec["stats"]["W3"]["median_over_theoretical_reference"] > 1.3
    assert rec["current_extended_fwhm_ratio"] == pytest.approx(1.3)
    assert rec["proposed_extended_fwhm_ratio"] == pytest.approx(1.3)
