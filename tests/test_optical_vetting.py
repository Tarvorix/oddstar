"""L006: offline tests for live optical time-series vetting (matched TESS detrending, robust ZTF, ASAS-SN)."""

from __future__ import annotations

import io
import zipfile
from pathlib import Path
from typing import Any

import numpy as np
import requests

from oddstar.config import load_config
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.vet.lc_checks import check_pipeline_agreement, check_second_survey, matched_detrend

RA, DEC, PMRA, PMDEC = 43.548272861, 5.882712239, 37.914, -97.175


class _Resp:
    def __init__(self, status_code: int, text: str = "", content: bytes | None = None) -> None:
        self.status_code = status_code
        self.text = text
        self.content = content if content is not None else text.encode("utf-8")


def _ztf_csv(band: str, oid: str, n: int, ra_off_arcsec: float = 0.0, mjd0: float = 58400.0) -> str:
    lines = ["oid,expid,hjd,mjd,mag,magerr,catflags,filtercode,ra,dec"]
    rng = np.random.default_rng(3)
    for i in range(n):
        mjd = mjd0 + i * 9.0
        jyr = 2000.0 + (mjd - 51544.5) / 365.25
        p_ra, p_dec = ot.propagate_position(RA, DEC, PMRA, PMDEC, 2016.0, jyr)
        mag = 14.35 + rng.normal(0, 0.02)
        flags = 32768 if i % 10 == 0 else 0
        lines.append(f"{oid},1,{mjd},{mjd},{mag:.4f},0.015,{flags},z{band},{p_ra + ra_off_arcsec / 3600.0:.7f},{p_dec:.7f}")
    return "\n".join(lines) + "\n"


def test_ztf_retries_shrinks_and_pages_after_timeouts(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    calls: list[dict[str, Any]] = []

    def fake_get(url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        p = dict(params or {})
        calls.append(p)
        # Every full-range request times out (the failure seen in live vetting); yearly pages succeed.
        if "TIME" not in p:
            raise requests.exceptions.ReadTimeout("read timeout")
        lo = float(p["TIME"].split()[0])
        if lo > 58194.0 + 365.25 * 3:  # only the first three pages hold data in this fixture
            return _Resp(200, "oid,expid,hjd,mjd,mag,magerr,catflags,filtercode,ra,dec\n")
        good = _ztf_csv(p["BANDNAME"], f"45210910000172{p['BANDNAME']}", 12, mjd0=lo + 1.0)
        # A second object 4" away along the cone edge must be rejected by the PM-track match
        far = _ztf_csv(p["BANDNAME"], "999", 5, ra_off_arcsec=4.0, mjd0=lo + 2.0).split("\n", 1)[1]
        return _Resp(200, good + far)

    monkeypatch.setattr(requests, "get", fake_get)
    summary, mans = ot.fetch_ztf_lightcurve(ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, cache=cache, now_mjd=61000.0, sleep=lambda s: None)
    assert summary is not None and summary["status"] == "ok"
    full_range = [c for c in calls if "TIME" not in c]
    assert len(full_range) == 6  # 3 retries x first cone radius x 2 bands; paging then succeeds at that radius
    assert any("paged 3.0\"" in e for e in summary["query_log"])
    # 3 pages x 12 rows, minus every 10th flagged row (2 per page) per band, far object rejected
    assert summary["bands"]["zg"]["n_epochs"] == 30
    assert summary["bands"]["zr"]["n_epochs"] == 30
    assert summary["bands"]["zg"]["object_id"] == "45210910000172g"
    assert summary["bands"]["zg"]["seasonal_range_mag"] is not None
    assert mans  # every successful page is in the content-addressed cache

    # Second call: full-range requests still time out, but every yearly page is served from the cache
    def only_full_range(url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        assert "TIME" not in (params or {}), "paged request should come from cache"
        raise requests.exceptions.ReadTimeout("read timeout")

    monkeypatch.setattr(requests, "get", only_full_range)
    summary2, _ = ot.fetch_ztf_lightcurve(ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, cache=cache, now_mjd=61000.0, sleep=lambda s: None)
    assert summary2 is not None and summary2["bands"]["zg"]["n_epochs"] == 30


def test_ztf_total_failure_is_reported_not_silent(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)

    def always_timeout(*a: Any, **k: Any) -> _Resp:
        raise requests.exceptions.ReadTimeout("read timeout")

    monkeypatch.setattr(requests, "get", always_timeout)
    summary, _ = ot.fetch_ztf_lightcurve(ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, cache=cache, now_mjd=58900.0, sleep=lambda s: None)
    assert summary is not None and summary["status"] == "query_failed"
    cfg, _, _ = load_config()
    res = check_second_survey({"lightcurves": {"optical_surveys": {"ztf": summary}}}, cfg, tmp_path / "ev")
    assert res.status == "not_evaluated"
    assert "query_failed" in res.note


def test_asassn_client_selection_and_camera_offsets(tmp_path: Path) -> None:
    import pandas as pd

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    rng = np.random.default_rng(5)
    jd = 2458320.5 + np.arange(400) * 4.0
    cams = np.where(np.arange(400) % 2 == 0, "bd", "bh")
    mag = 14.33 + np.where(cams == "bd", 0.05, -0.04) + rng.normal(0, 0.01, 400)
    data = pd.DataFrame(
        {
            "asas_sn_id": [7] * 400 + [8] * 10,
            "jd": list(jd) + list(jd[:10]),
            "mag": list(mag) + [12.0] * 10,
            "mag_err": [0.02] * 400 + [0.02] * 10,
            "quality": ["G"] * 395 + ["B"] * 5 + ["G"] * 10,
            "camera": list(cams) + ["bd"] * 10,
            "phot_filter": ["g"] * 410,
        }
    )
    cat = pd.DataFrame({"asas_sn_id": [7, 8], "ra_deg": [RA + 0.4 / 3600, RA + 12 / 3600], "dec_deg": [DEC - 0.6 / 3600, DEC]})

    class FakeClient:
        def cone_search(self, **kw: Any) -> Any:
            assert kw["catalog"] == "stellar_main" and kw["download"] is True
            return type("R", (), {"catalog_info": cat, "data": data})()

    summary, mans = ot.fetch_asassn_lightcurve(ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, cache=cache, client_factory=FakeClient)
    assert summary is not None and summary["n_epochs"] == 395
    g = summary["bands"]["g"]
    # Camera zero-point offsets (0.09 mag apart) are removed before the seasonal comparison
    assert g["seasonal_range_mag"] < 0.02
    assert all(s["season"] >= 2018 for s in g["seasons"])
    assert mans

    # Missing client is reported as None (caller records client_unavailable), never as coverage
    import importlib.util

    if importlib.util.find_spec("pyasassn") is None:
        summary_none, _ = ot.fetch_asassn_lightcurve(ra=RA + 1, dec=DEC, pmra=0, pmdec=0, cache=cache)
        assert summary_none is None


def _sector(rng: np.random.Generator, *, qlp_dip: float = 0.0, cut_dip: float = 0.0, cut_ramp: float = 0.03) -> dict[str, Any]:
    t = np.concatenate([np.arange(0.0, 12.0, 1 / 48), np.arange(13.5, 26.0, 1 / 48)]) + 2160.0
    # Slow, different trends in each extraction (QLP raw SAP vs. TESScut scattered-light ramps)
    q = 1.0 + 0.004 * np.sin((t - t[0]) / 3.0) + rng.normal(0, 0.002, len(t))
    c = 1.0 + cut_ramp * np.exp(-((t - t[0]) % 13.5) / 1.5) + rng.normal(0, 0.006, len(t))
    dip = (t > 2166.0) & (t < 2166.75)
    q[dip] *= 1.0 - qlp_dip
    c[dip] *= 1.0 - cut_dip
    return {"sector": 4, "time_btjd": t.tolist(), "flux_qlp_raw": q.tolist(), "flux_tesscut": c.tolist()}


def test_pipeline_agreement_matched_detrending(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    rng = np.random.default_rng(11)

    # Quiet star: different raw systematics in each extraction, identical filter -> both < 2% -> pass
    quiet = {"lightcurves": {"tess_tic_id": "387517877", "tess_sectors": [_sector(rng), _sector(rng)]}}
    res = check_pipeline_agreement(quiet, cfg, tmp_path / "q")
    assert res.status == "pass", res.note
    assert res.metrics["qlp_max_depth"] < 0.02 and res.metrics["tesscut_max_depth"] < 0.02
    assert len(res.metrics["sectors"]) == 2 and res.evidence

    # Without matched detrending the raw TESScut ramp alone would read as a >2% "dip"
    raw_c = np.asarray(quiet["lightcurves"]["tess_sectors"][0]["flux_tesscut"])
    assert 1.0 - raw_c.min() / np.median(raw_c) > 0.02

    # Real 10% dip seen by both extractions -> pass
    both = {"lightcurves": {"tess_sectors": [_sector(rng, qlp_dip=0.10, cut_dip=0.09)]}}
    assert check_pipeline_agreement(both, cfg, tmp_path / "b").status == "pass"

    # 10% dip only in QLP -> pipeline artefact -> fail, naming the sector
    qlp_only = {"lightcurves": {"tess_sectors": [_sector(rng, qlp_dip=0.10)]}}
    res_f = check_pipeline_agreement(qlp_only, cfg, tmp_path / "f")
    assert res_f.status == "fail" and "sector(s) [4]" in res_f.note

    # No sectors and no legacy arrays -> not_evaluated carrying the fetch log
    none = {"lightcurves": {"tess_sectors": [], "tess_query_log": ["s4: no QLP light curve"]}}
    res_n = check_pipeline_agreement(none, cfg, tmp_path / "n")
    assert res_n.status == "not_evaluated" and "no QLP" in res_n.note


def test_matched_detrend_preserves_dip_and_flattens_ramp() -> None:
    t = np.arange(0.0, 12.0, 1 / 48)
    f = 1.0 + 0.05 * np.exp(-t / 1.5)
    f[(t > 6.0) & (t < 6.5)] *= 0.9
    r = matched_detrend(t, f, 5.0)
    assert abs(1.0 - float(np.median(r))) < 0.005
    assert 0.08 < 1.0 - float(np.min(r)) < 0.12


def test_second_survey_reports_both_surveys(tmp_path: Path) -> None:
    """Since config 1.1.0 coverage alone no longer passes: the survey must be compared with a flagged signal."""
    cfg, _, _ = load_config()
    seasons = [
        {"season": y, "n_nights": 40, "median_mag": 14.35 + (0.06 if y == 2021 else 0.0), "median_err_mag": 0.004,
         "mjd_first": 51544.5 + (y - 2000) * 365.25 + 30, "mjd_mid": 51544.5 + (y - 2000) * 365.25 + 180,
         "mjd_last": 51544.5 + (y - 2000) * 365.25 + 330}
        for y in range(2018, 2025)
    ]
    star = {
        "lightcurves": {
            "optical_surveys": {
                "ztf": {"survey": "ZTF", "n_epochs": 1200, "status": "ok", "season_kind": "calendar",
                        "bands": {"zg": {"n_epochs": 600, "seasonal_range_mag": 0.06, "seasons": seasons}}},
                "asassn": {"survey": "ASAS-SN", "n_epochs": 0, "status": "client_unavailable", "query_log": []},
            }
        }
    }
    # No NEOWISE variation and no TESS dip: nothing to confirm or contradict -> not_evaluated, coverage still reported
    res = check_second_survey(star, cfg, tmp_path / "ev")
    assert res.status == "not_evaluated" and "No flagged signal" in res.note
    assert res.metrics["n_surveys"] == 1 and abs(res.metrics["max_seasonal_range_mag"] - 0.06) < 1e-9
    assert "not usable: ['asassn']" in res.note


def _fake_tess_files(tmp_path: Path) -> tuple[bytes, bytes]:
    from astropy.io import fits
    from astropy.wcs import WCS

    n = 300
    t = 1410.0 + np.arange(n) / 48.0
    q = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="TIME", format="D", array=t),
            fits.Column(name="SAP_FLUX", format="E", array=np.full(n, 5000.0)),
            fits.Column(name="QUALITY", format="J", array=np.zeros(n, dtype=int)),
        ]
    )
    qlp = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), q]).writeto(qlp)

    w = WCS(naxis=2)
    w.wcs.crpix = [7.0, 7.0]  # 1-based: target at 0-based pixel (6, 6)
    w.wcs.crval = [RA, DEC]
    w.wcs.cdelt = [-21.0 / 3600, 21.0 / 3600]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    yy, xx = np.mgrid[0:13, 0:13]
    star = 1000.0 * np.exp(-((xx - 6.0) ** 2 + (yy - 6.0) ** 2) / (2 * 0.7**2))
    cube = np.repeat((star + 20.0)[None, :, :], n, axis=0)
    cube[150:160] *= 0.9
    tab = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="TIME", format="D", array=t),
            fits.Column(name="FLUX", format="169E", dim="(13,13)", array=cube),
            fits.Column(name="QUALITY", format="J", array=np.zeros(n, dtype=int)),
        ]
    )
    ap = fits.ImageHDU(np.ones((13, 13), dtype=np.int32), header=w.to_header())
    fb = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), tab, ap]).writeto(fb)
    zb = io.BytesIO()
    with zipfile.ZipFile(zb, "w") as zf:
        zf.writestr("astrocut_43.5_5.9_13x13px.fits", fb.getvalue())
    return qlp.getvalue(), zb.getvalue()


def test_qlp_and_tesscut_readers(tmp_path: Path) -> None:
    qlp_bytes, zip_bytes = _fake_tess_files(tmp_path)
    tq, fq = ot.read_qlp_raw(qlp_bytes)
    assert len(tq) == 300 and float(fq[0]) == 5000.0
    tc, fc = ot.extract_tesscut_aperture(zip_bytes, RA, DEC, 0.0, 0.0)
    assert len(tc) == 300
    # Local annulus removes the 20 e/s pedestal; the injected 10% dimming is recovered
    assert abs(float(fc[150]) / float(fc[0]) - 0.9) < 0.01
    tm, a, b = ot.match_cadences(tq, fq, tc, fc)
    assert len(tm) == 300
    assert ot.qlp_url("387517877", 4).endswith("s0004/0000/0003/8751/7877/hlsp_qlp_tess_ffi_s0004-0000000387517877_tess_v01_llc.fits")
