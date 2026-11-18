"""Survey-scale TESS light curves for optical_dips (offline, fake HTTP) and honest not-searched accounting."""

from __future__ import annotations

import io
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from oddstar.channels import build_channel
from oddstar.channels.optical_dips.channel import search_coverage_summary, tess_search_status
from oddstar.cli import CHANNELS_WITHOUT_LIVE_SOURCE, dispatch
from oddstar.config import load_config, repo_root
from oddstar.errors import ArchiveAccessError
from oddstar.io import optical_timeseries as ot
from oddstar.io import tess_bulk as tb
from oddstar.io.cache import RawCache


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


def _qlp_fits(t0: float = 1500.0, days: float = 12.0, cadence_s: float = 600.0, level: float = 5000.0,
              dip: tuple[float, float, float] | None = None, seed: int = 1, bad_every: int = 0) -> bytes:
    """A QLP-like light-curve file: TIME, SAP_FLUX, QUALITY. `dip` is (centre BTJD, full width in days, depth)."""
    from astropy.io import fits

    rng = np.random.default_rng(seed)
    t = t0 + np.arange(int(days * 86400.0 / cadence_s)) * cadence_s / 86400.0
    f = level * (1.0 + rng.normal(0.0, 0.002, len(t)))
    if dip is not None:
        f = np.where(np.abs(t - dip[0]) <= dip[1] / 2.0, f * (1.0 - dip[2]), f)
    q = np.zeros(len(t), dtype=np.int32)
    if bad_every:
        q[::bad_every] = 4096
    hdu = fits.BinTableHDU.from_columns([
        fits.Column(name="TIME", format="D", array=t), fits.Column(name="SAP_FLUX", format="E", array=f.astype(np.float32)),
        fits.Column(name="QUALITY", format="J", array=q)])
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    return buf.getvalue()


class _Resp:
    def __init__(self, status_code: int, content: bytes = b"") -> None:
        self.status_code = status_code
        self.content = content


class _Bucket:
    """Fake MAST: TIC ids and sector lists per star, QLP files per (TIC, sector). Counts every file request."""

    def __init__(self) -> None:
        self.tic: dict[str, str | None] = {}
        self.sectors: dict[str, list[int]] = {}
        self.files: dict[str, bytes | int | Exception] = {}
        self.file_requests: list[str] = []
        self.lookups = 0

    def add(self, gaia_id: str, tic: str | None, sectors: dict[int, bytes | int | Exception]) -> None:
        self.tic[gaia_id] = tic
        self.sectors[gaia_id] = sorted(sectors)
        if tic is not None:
            for sec, body in sectors.items():
                self.files[ot.qlp_url(tic, sec)] = body

    def install(self, monkeypatch: pytest.MonkeyPatch) -> None:
        by_pos: dict[tuple[float, float], str] = {}

        def lookup_tic_id(*, ra: float, dec: float, gaia_id: str, cache: RawCache, sleep: Any = None) -> tuple[str | None, list]:
            self.lookups += 1
            if isinstance(self.tic[gaia_id], Exception):
                raise self.tic[gaia_id]  # type: ignore[misc]
            by_pos[(round(ra, 6), round(dec, 6))] = gaia_id
            return self.tic[gaia_id], []

        def list_tess_sectors(*, ra: float, dec: float, cache: RawCache, sleep: Any = None) -> tuple[list[int], list]:
            return list(self.sectors[by_pos[(round(ra, 6), round(dec, 6))]]), []

        def get(url: str, params: Any, **kwargs: Any) -> _Resp:
            self.file_requests.append(url)
            body = self.files[url]
            if isinstance(body, Exception):
                raise body
            return _Resp(body) if isinstance(body, int) else _Resp(200, body)

        monkeypatch.setattr(ot, "lookup_tic_id", lookup_tic_id)
        monkeypatch.setattr(ot, "list_tess_sectors", list_tess_sectors)
        monkeypatch.setattr(ot, "_get_with_retries", get)


def _star(sid: str, ra: float = 10.0, dec: float = -20.0) -> dict[str, Any]:
    return {"source_id": sid, "ra_deg": ra + 0.01 * int(sid), "dec_deg": dec, "pmra_mas_yr": 3.0, "pmdec_mas_yr": -2.0,
            "lightcurves": {"tess_time_btjd": [], "tess_flux_qlp": [], "tess_flux_tesscut": [], "ztf_time_mjd": [], "ztf_flux": []}}


def _survey() -> tuple[list[dict[str, Any]], _Bucket]:
    b = _Bucket()
    # 1: two quiet sectors and one with no QLP file; 2: a 1-day, 8% dip in the middle of its second sector;
    # 3: no TIC entry; 4: every sector fails; 5: sectors exist but none has a QLP file; 6: no sector covers it
    b.add("1", "111", {3: _qlp_fits(1385.0, seed=3), 30: _qlp_fits(2115.0, seed=4), 57: 404})
    b.add("2", "222", {5: _qlp_fits(1438.0, seed=5), 32: _qlp_fits(2174.0, seed=6, dip=(2180.0, 1.0, 0.08))})
    b.add("3", None, {})
    b.add("4", "444", {7: ArchiveAccessError("bucket timeout", retryable=True), 8: ArchiveAccessError("bucket timeout", retryable=True)})
    b.add("5", "555", {9: 403, 10: 404})
    b.add("6", "666", {})
    return [_star(str(i)) for i in range(1, 7)], b


def test_reduce_keeps_good_cadences_in_half_hour_bins() -> None:
    red = tb.reduce_qlp_fits(_qlp_fits(days=2.0, cadence_s=600.0, bad_every=3))
    assert red["n_good"] == 192 and red["cadence_s"] == pytest.approx(600.0, abs=1.0)  # 288 cadences, every third flagged
    assert len(red["time"]) == len(red["flux"]) == len(red["n_per_bin"]) == 96
    assert set(red["n_per_bin"]) == {2} and sum(red["n_per_bin"]) == 192
    assert np.allclose(np.diff(red["time"]), ot.TESS_BIN_DAYS, atol=0.005)
    assert np.median(red["flux"]) == pytest.approx(5000.0, rel=0.002)
    short = tb.reduce_qlp_fits(_qlp_fits(days=0.2))
    assert short["time"] == [] and short["n_good"] < tb.MIN_GOOD_CADENCES


def test_sector_answers_are_cached_including_missing_files_but_not_failures(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    b = _Bucket()
    b.add("1", "111", {3: _qlp_fits(), 4: 404, 5: b"not a fits file", 6: ArchiveAccessError("timeout", retryable=True)})
    b.install(monkeypatch)
    cache = RawCache(tmp_path / "raw")
    rec, man, hit = tb.fetch_qlp_sector_binned("111", 3, cache)
    assert rec["status"] == "ok" and not hit and rec["raw_bytes"] > 30_000 and len(rec["raw_sha256"]) == 64
    assert man.byte_size < rec["raw_bytes"] / 3, "the cache keeps the reduced series, not the file"
    rec2, man2, hit2 = tb.fetch_qlp_sector_binned("111", 3, cache)
    assert hit2 and rec2 == rec and man2.manifest_id == man.manifest_id
    absent, _, _ = tb.fetch_qlp_sector_binned("111", 4, cache)
    assert absent["status"] == "absent" and absent["http_status"] == 404
    assert tb.fetch_qlp_sector_binned("111", 4, cache)[2] is True
    n = len(b.file_requests)
    for sector, err in ((5, ValueError), (6, ArchiveAccessError)):
        for _ in range(2):  # failures are not stored: the second call asks the archive again
            with pytest.raises(err):
                tb.fetch_qlp_sector_binned("111", sector, cache)
    assert len(b.file_requests) == n + 4


def test_star_light_curve_is_normalised_per_sector_and_limited_to_max_sectors(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    b = _Bucket()
    b.add("1", "111", {1: _qlp_fits(1325.0, level=5000.0), 2: 404, 3: _qlp_fits(1385.0, level=9000.0, seed=2),
                       4: _qlp_fits(1411.0, seed=3), 5: _qlp_fits(1438.0, days=0.2)})
    b.install(monkeypatch)
    lc, rec, _, _ = tb.fetch_star_tess(_star("1"), RawCache(tmp_path / "raw"), max_sectors=2)
    assert rec["status"] == "ok" and rec["tic_id"] == "111"
    assert rec["sectors_used"] == [1, 3] and rec["sectors_without_qlp"] == [2] and rec["sectors_beyond_max"] == [4, 5]
    assert len(lc["tess_time_btjd"]) == len(lc["tess_flux_qlp"]) == len(lc["tess_sector"]) == rec["n_points"]
    sec = np.asarray(lc["tess_sector"])
    flux = np.asarray(lc["tess_flux_qlp"])
    assert sorted(set(sec.tolist())) == [1, 3]
    assert np.median(flux[sec == 1]) == pytest.approx(1.0, abs=1e-3) and np.median(flux[sec == 3]) == pytest.approx(1.0, abs=1e-3)
    lc_all, rec_all, _, _ = tb.fetch_star_tess(_star("1"), RawCache(tmp_path / "raw2"), max_sectors=8)
    assert rec_all["sectors_used"] == [1, 3, 4] and rec_all["sectors_too_short"] == [5]


def test_survey_fetch_reports_every_status_and_resumes_from_cache(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, b = _survey()
    b.install(monkeypatch)
    cache = RawCache(tmp_path / "raw")
    seen: list[tuple[int, int, str, str]] = []
    out, summary = tb.enrich_sample_tess(stars, cfg, cache, sleep=lambda s: None, progress=lambda *a: seen.append(a))
    assert [s["source_id"] for s in out] == ["1", "2", "3", "4", "5", "6"]
    statuses = [s["tess_fetch"]["status"].split(":")[0] for s in out]
    assert statuses == ["ok", "ok", "no_coverage", "failed", "no_coverage", "no_coverage"]
    assert [s["vetting_archive"]["query_status"]["tess_qlp"].split(":")[0] for s in out] == statuses
    assert summary["n_ok"] == 2 and summary["n_no_coverage"] == 3 and summary["n_failed"] == 1 and summary["failed_ids"] == ["4"]
    assert summary["n_sectors_used"] == 4 and summary["n_not_fetched"] == 0
    assert out[3]["tess_fetch"]["attempts"] == int(cfg.channels["optical_dips"]["tess_bulk_max_attempts"])
    assert out[3]["lightcurves"]["tess_time_btjd"] == []
    assert stars[0]["lightcurves"]["tess_time_btjd"] == [], "input stars are not modified"
    assert [a[0] for a in seen] == [1, 2, 3, 4, 5, 6] and seen[-1][1] == 6

    n = len(b.file_requests)
    out2, summary2 = tb.enrich_sample_tess(stars, cfg, RawCache(tmp_path / "raw"), sleep=lambda s: None, workers=1)
    failing = [u for u in b.file_requests[n:]]
    assert all("444" in u.replace("/", "") or "0444" in u for u in failing), "only the star that failed is asked for again"
    strip = lambda rows: [{k: v for k, v in s.items()} for s in rows]  # noqa: E731
    assert json.dumps(strip(out2), sort_keys=True) == json.dumps(strip(out), sort_keys=True), "worker count and cache do not change results"

    out3, summary3 = tb.enrich_sample_tess(stars, cfg, RawCache(tmp_path / "raw"), max_stars=2, sleep=lambda s: None)
    assert summary3["n_requested"] == 2 and summary3["n_not_fetched"] == 4
    assert [s["tess_fetch"]["status"] for s in out3[2:]] == ["not_fetched"] * 4


def test_failed_lookup_is_a_failed_star_not_a_quiet_one(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    b = _Bucket()
    b.add("1", "111", {3: _qlp_fits()})
    b.tic["1"] = ArchiveAccessError("MAST unreachable", retryable=True)  # type: ignore[assignment]
    b.install(monkeypatch)
    out, summary = tb.enrich_sample_tess([_star("1")], cfg, RawCache(tmp_path / "raw"), sleep=lambda s: None)
    assert out[0]["tess_fetch"]["status"].startswith("failed: ArchiveAccessError") and summary["n_failed"] == 1
    assert tess_search_status(out[0], cfg) == "failed"


def test_dip_search_finds_the_dip_and_counts_only_searched_stars(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, b = _survey()
    b.install(monkeypatch)
    out, _ = tb.enrich_sample_tess(stars, cfg, RawCache(tmp_path / "raw"), sleep=lambda s: None)
    never_fetched = _star("7")
    ch = build_channel("optical_dips", cfg, tmp_path)
    sample = [*out, never_fetched]
    meas = {s["source_id"]: ch.measure(s, {}) for s in sample}
    scores = {sid: ch.score(m) for sid, m in meas.items()}

    assert scores["2"].flagged and not scores["1"].flagged
    v = meas["2"].values
    assert v["best_sector"] == 32.0 and abs(v["best_time_btjd"] - 2180.0) < 0.6
    assert v["dip_depth"] == pytest.approx(0.08, abs=0.015) and v["n_segments_searched"] == 2.0
    assert meas["1"].values["n_segments_searched"] == 2.0 and meas["1"].values["dip_depth"] < 0.01
    for sid, cat in (("3", "no_coverage"), ("4", "failed"), ("5", "no_coverage"), ("6", "no_coverage"), ("7", "no_light_curve")):
        assert meas[sid].extra["search_status"] == cat and meas[sid].values["n_points_searched"] == 0.0
        assert scores[sid].flagged is False and scores[sid].summary.startswith(f"NOT SEARCHED ({cat})")

    cov = search_coverage_summary(sample, cfg)
    assert cov["n_input"] == 7 and cov["n_searched"] == 2 and cov["n_not_searched"] == 5
    assert cov["n_no_coverage"] == 3 and cov["n_failed"] == 1 and cov["n_no_light_curve"] == 1 and cov["n_not_fetched"] == 0
    assert cov["not_searched_ids"] == {"no_coverage": ["3", "5", "6"], "failed": ["4"], "no_light_curve": ["7"]}


def _lc(times: np.ndarray, flux: np.ndarray, sector: np.ndarray) -> dict[str, Any]:
    return {"source_id": "9", "lightcurves": {"tess_time_btjd": times.tolist(), "tess_flux_qlp": flux.tolist(), "tess_sector": sector.tolist()}}


def test_segments_are_searched_alone_and_edges_are_not_searched(tmp_path: Path, cfg: Any) -> None:
    ch = build_channel("optical_dips", cfg, tmp_path)
    rng = np.random.default_rng(11)
    step = ot.TESS_BIN_DAYS
    t1 = 1500.0 + np.arange(600) * step  # 12.5 days
    t2 = 1513.5 + np.arange(600) * step  # after a one-day orbit gap
    t3 = 2200.0 + np.arange(600) * step  # another sector, years later
    t = np.r_[t1, t2, t3]
    sec = np.r_[np.full(1200, 5.0), np.full(600, 31.0)]
    noise = rng.normal(0.0, 0.002, len(t))

    # Offsets between segments (per-sector normalisation never matches exactly) are not dips
    stepped = 1.0 + noise + np.r_[np.zeros(600), np.full(600, -0.06), np.full(600, 0.05)]
    m = ch.measure(_lc(t, stepped, sec), {})
    assert m.values["n_segments_searched"] == 3.0 and not ch.score(m).flagged and m.values["dip_depth"] < 0.01

    # A settling ramp in the first hours of a segment sits in the unsearched edge zone
    ramp = 1.0 + noise.copy()
    ramp[600:612] -= np.linspace(0.12, 0.0, 12)
    ramp[1188:1200] -= 0.10
    m = ch.measure(_lc(t, ramp, sec), {})
    assert not ch.score(m).flagged
    margin = float(cfg.channels["optical_dips"]["systematics_edge_margin_days"])
    assert m.values["n_points_searched"] == 3 * float(((t1 >= t1[0] + margin) & (t1 <= t1[-1] - margin)).sum())

    # The same dip in the middle of a segment is found, with its sector and time
    mid = 1.0 + noise.copy()
    mid[1480:1528] -= 0.10
    m = ch.measure(_lc(t, mid, sec), {})
    assert ch.score(m).flagged and m.values["best_sector"] == 31.0 and abs(m.values["best_time_btjd"] - t[1504]) < 0.3
    assert m.values["dip_depth"] == pytest.approx(0.10, abs=0.01)

    # One bad point is not a dip: depth is the median inside the box
    spike = 1.0 + noise.copy()
    spike[300] = 0.2
    m = ch.measure(_lc(t, spike, sec), {})
    assert not ch.score(m).flagged and m.values["dip_depth"] < 0.02

    # Too little data is not searched
    short = _lc(t1[:15], np.ones(15), np.full(15, 5.0))
    assert tess_search_status(short, cfg) == "no_light_curve" and ch.measure(short, {}).values["n_segments_searched"] == 0.0
    edge_only = _lc(t1[:40], np.ones(40), np.full(40, 5.0))  # 0.8 days: nothing is left inside two 0.5-day margins
    assert tess_search_status(edge_only, cfg) == "no_light_curve"


def test_live_search_refuses_channels_without_a_real_data_source(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    assert set(CHANNELS_WITHOUT_LIVE_SOURCE) == {"spectral_lines", "vanishing", "radio_narrowband"}
    for channel in CHANNELS_WITHOUT_LIVE_SOURCE:
        code, payload = dispatch(["search", channel, "--live", "--limit", "5", "--run-id", f"t-{channel}", "--json"])
        assert code != 0 and payload["ok"] is False and "no real-data source" in json.dumps(payload), channel
        assert not (repo_root() / "runs" / f"t-{channel}").exists(), "a refused search leaves no run directory"


def test_live_dip_search_records_coverage_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """search -> null-report -> run record on a faked live sample: only stars with a light curve count as searched."""
    import shutil

    from oddstar import cli
    from oddstar.run_records import build_run_record

    stars, b = _survey()
    for s in stars:
        s.update({"distance_pc": 150.0, "teff_k": 5500.0})
    b.install(monkeypatch)
    enriched: list[str] = []
    monkeypatch.setattr(cli, "fetch_live_halo_bulge_sample", lambda cfg, cache, **kw: ([json.loads(json.dumps(s)) for s in stars], []))
    monkeypatch.setattr(cli, "enrich_star_live_for_vetting", lambda star, cfg, cache: (enriched.append(star["source_id"]) or dict(star), []))
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    monkeypatch.setattr("time.sleep", lambda s: None)
    run_dir = repo_root() / "runs" / "test-tess-survey"
    try:
        argv = ["search", "optical_dips", "--live", "--region", "halo", "--limit", "6", "--run-id", "test-tess-survey", "--json"]
        code, payload = dispatch([*argv, "--defer-enrichment"])
        assert code == 0, payload
        assert enriched == [] and not (run_dir / "live_stars.json").exists()
        d = payload["data"]
        assert d["n_input"] == 6 and d["n_searched"] == 2 and d["flagged_ids"] == ["2"]
        assert d["coverage"]["n_no_coverage"] == 3 and d["coverage"]["n_failed"] == 1 and d["tess_fetch"]["n_ok"] == 2
        assert any("4 of 6 stars not searched" in w for w in payload["warnings"])

        code, payload = dispatch(argv)
        assert code == 0 and enriched == ["2"] and (run_dir / "live_stars.json").exists()
        saved = json.loads((run_dir / "search_optical_dips.json").read_text(encoding="utf-8"))
        assert saved["n_searched"] == 2 and saved["not_searched_ids"] == {"no_coverage": ["3", "5", "6"], "failed": ["4"]}
        assert saved["tess_fetch"]["failed_ids"] == ["4"] and len(saved["sample_stars"]) == 6

        code, payload = dispatch(["null-report", "test-tess-survey", "--json"])
        assert code == 0, payload
        rep = json.loads((run_dir / "null_report_optical_dips.json").read_text(encoding="utf-8"))
        assert rep["n_searched"] == 2 and rep["n_flagged"] == 1

        rec, samples = build_run_record(run_dir)
        s = rec["searches"][0]
        assert s["coverage"]["n_searched"] == 2 and s["coverage"]["n_not_searched"] == 4 and s["coverage"]["n_no_light_curve"] == 0
        assert s["not_searched_ids"] == ["3", "4", "5", "6"] and s["tess_fetch"]["n_failed"] == 1 and samples[0]["n"] == 6
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)
