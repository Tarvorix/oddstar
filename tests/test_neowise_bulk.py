"""Survey-scale NEOWISE-R ingest for ir_timedomain (offline, fake HTTP) and honest not-searched accounting."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import pytest
import requests

from oddstar.channels import build_channel
from oddstar.channels.ir_timedomain.channel import neowise_search_status, search_coverage_summary
from oddstar.config import load_config
from oddstar.dossier.null_report import write_null_report
from oddstar.errors import ArchiveAccessError
from oddstar.io import live_survey as ls
from oddstar.io.cache import RawCache
from oddstar.io.neowise_bulk import bin_neowise_visits, enrich_sample_neowise
from oddstar.io.physics import propagate_position
from oddstar.schemas import SampleCutStep, ValidationCertificate

RA, DEC, PMRA, PMDEC = 201.123456, -33.654321, -45.2, 12.7

# sha256 of json.dumps(out["neowise"]) from enrich_star_live_for_vetting section 3c on `_fake_neowise_rows()`, recorded
# with the pre-refactor inline implementation (commit f1913c3). The shared bin_neowise_visits must reproduce it.
SECTION_3C_GOLDEN_SHA256 = "5fc15e2b934f264ce5841d64c57a8d70c14db16bbbaa02f1f95f1799a6e845f0"


def _fake_neowise_rows(ra: float = RA, dec: float = DEC, pmra: float = PMRA, pmdec: float = PMDEC, n_visits: int = 21,
                       step_mag: float = 0.0) -> list[dict[str, Any]]:
    rows = []
    for v in range(n_visits):
        t0 = 56650.0 + 182.6 * v
        for k in range(13):
            jyr = 2000.0 + (t0 - 51544.5) / 365.25
            pra, pdec = propagate_position(ra, dec, pmra, pmdec, 2016.0, jyr)
            dra = (0.11 * math.sin(v + k)) / 3600.0 / math.cos(math.radians(dec))
            ddec = (0.09 * math.cos(2 * v - k)) / 3600.0
            jump = -step_mag if v >= n_visits // 2 else 0.0
            rows.append({"ra": f"{pra + dra:.7f}", "dec": f"{pdec + ddec:.7f}", "mjd": f"{t0 + 0.13 * k:.6f}",
                         "w1mpro": f"{10.512 + jump + 0.02 * math.sin(3 * v + k):.3f}", "w1sigmpro": "0.021",
                         "w2mpro": f"{10.488 + jump + 0.03 * math.cos(v * k):.3f}" if (v, k) != (4, 2) else None,
                         "w2sigmpro": "0.025", "qual_frame": "0" if (v, k) == (7, 5) else "10", "qi_fact": "1.0"})
    return rows


class _Resp:
    def __init__(self, status_code: int, text: str = "", json_data: Any = None) -> None:
        self.status_code = status_code
        self.text = text if text or json_data is None else json.dumps(json_data)
        self._json = json_data

    def json(self) -> Any:
        return self._json if self._json is not None else json.loads(self.text)


def _ipac(cols: list[str], rows: list[dict[str, Any]]) -> str:
    w = {c: max([len(c)] + [len(str(r.get(c) if r.get(c) is not None else "null")) for r in rows]) + 2 for c in cols}
    lines = ["\\fixlen = T"] + ["|" + "|".join(h.rjust(w[c]) if h else c.rjust(w[c]) for c in cols) + "|" for h in ("", "char", " ", "null")]
    lines += [" " + " ".join(str(r.get(c) if r.get(c) is not None else "null").rjust(w[c]) for c in cols) + " " for r in rows]
    return "\n".join(lines) + "\n"


class FakeIrsa:
    """IRSA Gator/TAP stand-in: per-position behaviour 'ok' (rows), 'empty', 'fail' or 'flaky' (fails once)."""

    COLS = ["ra", "dec", "mjd", "w1mpro", "w1sigmpro", "w2mpro", "w2sigmpro", "qual_frame", "qi_fact"]

    def __init__(self, behaviour: dict[str, str], rows_by_key: dict[str, list[dict[str, Any]]]) -> None:
        self.behaviour = behaviour
        self.rows_by_key = rows_by_key
        self.calls: list[tuple[str, str]] = []
        self._flaky_seen: set[str] = set()

    @staticmethod
    def key(ra: float, dec: float) -> str:
        return f"{ra:.6f} {dec:+.6f}"

    def _mode(self, key: str) -> str:
        mode = self.behaviour[key]
        if mode == "flaky":
            if key in self._flaky_seen:
                return "ok"
            self._flaky_seen.add(key)
            return "fail"
        return mode

    def get(self, url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        p = dict(params or {})
        assert "Gator" in url and p["catalog"] == "neowiser_p1bs_psd", f"unexpected GET {url} {p}"
        assert float(p["radius"]) == 3.0, "survey mode must not run the comparison-star cone"
        key = p["objstr"]
        self.calls.append(("gator", key))
        mode = self._mode(key)
        if mode == "fail":
            raise requests.exceptions.ConnectionError(f"simulated Gator outage for {key}")
        rows = self.rows_by_key.get(key, []) if mode == "ok" else []
        return _Resp(200, _ipac(self.COLS, rows))

    def post(self, url: str, data: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        q = str((data or {}).get("QUERY", ""))
        assert "neowiser_p1bs_psd" in q, f"unexpected POST {q[:80]}"
        self.calls.append(("tap", q))
        return _Resp(503, "Service Unavailable")


@pytest.fixture
def cfg() -> Any:
    c, _, _ = load_config()
    return c


def _star(sid: str, ra: float, dec: float) -> dict[str, Any]:
    return {"source_id": sid, "name": f"Gaia DR3 {sid}", "ra_deg": ra, "dec_deg": dec, "pmra_mas_yr": PMRA, "pmdec_mas_yr": PMDEC,
            "mags": {"W1": 10.5, "W2": 10.5}, "flux_mjy": {}, "flux_err_mjy": {}, "lightcurves": {"ztf_flux": []},
            "neowise": {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []},
            "vetting_archive": {"query_status": {}}, "live_enriched": False}


def test_section_3c_output_is_byte_identical_after_refactor(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    cache = RawCache(tmp_path / "raw")

    def no_net(*a: Any, **k: Any) -> Any:
        raise requests.exceptions.ConnectionError("offline test")

    monkeypatch.setattr(requests, "get", no_net)
    monkeypatch.setattr(requests, "post", no_net)

    def gator(**kw: Any) -> Any:
        if kw["catalog"] == "neowiser_p1bs_psd" and kw["radius_arcsec"] == 3.0:
            man = cache.put(archive="neowise", release="neowise_r", request_signature="fake", query_text="q", payload=b"[]", row_count=0)
            return _fake_neowise_rows(), man, False
        raise ArchiveAccessError("offline")  # the comparison cone records this text, which is part of the golden hash

    monkeypatch.setattr(ls, "_gator_cone_cached", gator)
    star = {"source_id": "123", "name": "x", "ra_deg": RA, "dec_deg": DEC, "pmra_mas_yr": PMRA, "pmdec_mas_yr": PMDEC,
            "mags": {"G": 11.0, "W1": 10.5, "W2": 10.5}, "flux_mjy": {}, "flux_err_mjy": {}, "vetting_archive": {"query_status": {}}}
    out, _ = ls.enrich_star_live_for_vetting(star, cfg, cache, enrich_optical=False, sleep=lambda s: None)
    assert out["vetting_archive"]["query_status"]["neowise"] == "ok"
    assert out["neowise"]["comparison"]["status"] == "query_failed"
    assert hashlib.sha256(json.dumps(out["neowise"]).encode()).hexdigest() == SECTION_3C_GOLDEN_SHA256


def test_bin_neowise_visits_drops_bad_frames_and_keeps_pm_vectors() -> None:
    binned = bin_neowise_visits(_fake_neowise_rows(), ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC)
    assert binned is not None
    neo, windows, w1_mag = binned
    assert len(neo["mjd"]) == 21 and len(windows) == 21
    # one frame with qual_frame=0 and one with null W2 are dropped
    assert sum(neo["visit_n_frames"]) == 21 * 13 - 2
    assert len(neo["pm_residual_ra_arcsec"]) == len(neo["pm_residual_dec_arcsec"]) == 21
    assert neo["median_pm_residual_arcsec"] < 0.2
    assert 10.4 < w1_mag < 10.6
    assert list(neo)[:5] == ["mjd", "w1_mjy", "w2_mjy", "w1_err_mjy", "w2_err_mjy"]
    assert bin_neowise_visits([{"mjd": "56650.0", "w1mpro": "10.0", "w2mpro": "10.0", "qual_frame": "0"}],
                              ra=RA, dec=DEC, pmra=0.0, pmdec=0.0) is None


def _survey(n_ok: int = 2) -> tuple[list[dict[str, Any]], FakeIrsa]:
    stars = []
    behaviour: dict[str, str] = {}
    rows: dict[str, list[dict[str, Any]]] = {}
    modes = ["ok"] * n_ok + ["empty", "fail", "flaky"]
    for i, mode in enumerate(modes):
        ra, dec = 150.0 + i, 20.0 + 0.5 * i
        s = _star(f"10{i}", ra, dec)
        stars.append(s)
        key = FakeIrsa.key(ra, dec)
        behaviour[key] = mode
        rows[key] = _fake_neowise_rows(ra, dec, PMRA, PMDEC, step_mag=0.3 if i == 0 else 0.0)
    return stars, FakeIrsa(behaviour, rows)


def test_enrich_sample_neowise_statuses_retries_rate_limit_and_summary(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, fake = _survey()
    monkeypatch.setattr(requests, "get", fake.get)
    monkeypatch.setattr(requests, "post", fake.post)
    sleeps: list[float] = []
    progress: list[tuple[int, int, str, str]] = []
    out, summary = enrich_sample_neowise(stars, cfg, RawCache(tmp_path / "raw"), progress=lambda *a: progress.append(a), sleep=sleeps.append)

    by_id = {s["source_id"]: s for s in out}
    assert [s["source_id"] for s in out] == [s["source_id"] for s in stars]
    assert by_id["100"]["neowise_fetch"]["status"] == "ok" and len(by_id["100"]["neowise"]["mjd"]) == 21
    assert "comparison" not in by_id["100"]["neowise"], "the comparison-star cone is vetting-only"
    assert by_id["102"]["neowise_fetch"]["status"] == "no_coverage" and by_id["102"]["neowise"]["mjd"] == []
    failed = by_id["103"]
    assert failed["neowise_fetch"]["status"].startswith("failed: ") and failed["neowise_fetch"]["attempts"] == 3
    assert failed["vetting_archive"]["query_status"]["neowise"] == failed["neowise_fetch"]["status"]
    assert failed["neowise"]["mjd"] == []
    assert by_id["104"]["neowise_fetch"]["status"] == "ok" and by_id["104"]["neowise_fetch"]["attempts"] == 2
    assert summary["n_requested"] == 5 and summary["n_ok"] == 3 and summary["n_no_coverage"] == 1 and summary["n_failed"] == 1
    assert summary["failed_ids"] == ["103"] and summary["n_not_fetched"] == 0 and summary["elapsed_s"] >= 0.0
    assert [p[0] for p in progress] == [1, 2, 3, 4, 5] and all(p[1] == 5 for p in progress)
    # Bounded retries: 3 Gator + 3 TAP attempts for the dead star, nothing more
    assert sum(1 for c in fake.calls if c == ("gator", FakeIrsa.key(153.0, 21.5))) == 3
    assert sum(1 for c in fake.calls if c[0] == "tap") == 3 + 1
    # Rate limit (2 rps -> 0.5 s spacing) and exponential backoff (2 s, 4 s) went through the injected sleep
    rps = cfg.archives["neowise"].rate_limit_rps
    assert any(0.0 < s_ <= 1.0 / rps for s_ in sleeps)
    backoff = float(cfg.channels["ir_timedomain"]["neowise_bulk_retry_backoff_s"])
    assert sleeps.count(backoff) == 2 and sleeps.count(2.0 * backoff) == 1


def test_enrich_sample_neowise_resumes_from_cache_without_refetching(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, fake = _survey()
    monkeypatch.setattr(requests, "get", fake.get)
    monkeypatch.setattr(requests, "post", fake.post)
    cache_dir = tmp_path / "raw"

    def kill_after_two(done: int, total: int, sid: str, status: str) -> None:
        if done == 2:
            raise KeyboardInterrupt

    with pytest.raises(KeyboardInterrupt):
        enrich_sample_neowise(stars, cfg, RawCache(cache_dir), progress=kill_after_two, sleep=lambda s: None, workers=1)
    first_keys = [c[1] for c in fake.calls]
    assert first_keys == [FakeIrsa.key(150.0, 20.0), FakeIrsa.key(151.0, 20.5)]

    fake.calls.clear()
    sleeps: list[float] = []
    out, summary = enrich_sample_neowise(stars, cfg, RawCache(cache_dir), sleep=sleeps.append)  # fresh cache object = new process
    resumed_keys = {c[1] for c in fake.calls if c[0] == "gator"}
    assert not (resumed_keys & set(first_keys)), "finished stars must come from the cache"
    assert summary["n_cache_hits"] == 2 and summary["n_ok"] == 3

    # A third run is free: only the dead star is retried (failures are never cached)
    fake.calls.clear()
    out3, summary3 = enrich_sample_neowise(stars, cfg, RawCache(cache_dir), sleep=lambda s: None)
    assert {c[1] for c in fake.calls if c[0] == "gator"} == {FakeIrsa.key(153.0, 21.5)}
    assert summary3["n_network_requests"] == 3 and summary3["n_cache_hits"] == 4
    assert out3[0]["neowise"] == out[0]["neowise"]


def test_max_stars_marks_the_rest_not_fetched(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, fake = _survey()
    monkeypatch.setattr(requests, "get", fake.get)
    monkeypatch.setattr(requests, "post", fake.post)
    out, summary = enrich_sample_neowise(stars, cfg, RawCache(tmp_path / "raw"), max_stars=2, sleep=lambda s: None)
    assert summary["n_requested"] == 2 and summary["n_not_fetched"] == 3
    assert [s["neowise_fetch"]["status"] for s in out[2:]] == ["not_fetched"] * 3
    cov = search_coverage_summary(out, cfg)
    assert cov["n_searched"] == 2 and cov["n_not_fetched"] == 3


def test_ir_timedomain_reports_failed_and_uncovered_stars_as_not_searched(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, fake = _survey()
    monkeypatch.setattr(requests, "get", fake.get)
    monkeypatch.setattr(requests, "post", fake.post)
    out, _ = enrich_sample_neowise(stars, cfg, RawCache(tmp_path / "raw"), sleep=lambda s: None)
    never_queried = _star("999", 10.0, 10.0)  # live record that never went through a NEOWISE query
    stale_failed = dict(out[1], source_id="998", vetting_archive={"query_status": {"neowise": "failed: Timeout: x"}})  # data but failed
    sample = [*out, never_queried, stale_failed]

    ch = build_channel("ir_timedomain", cfg, tmp_path)
    scores = {s["source_id"]: ch.score(ch.measure(s, {})) for s in sample}
    meas = {s["source_id"]: ch.measure(s, {}) for s in sample}
    assert scores["100"].flagged, "a real 0.3 mag step must still be found"
    for sid, cat in (("102", "no_coverage"), ("103", "failed"), ("999", "not_fetched"), ("998", "failed")):
        assert meas[sid].extra["search_status"] == cat
        assert meas[sid].values["searched"] == 0.0
        assert scores[sid].flagged is False and scores[sid].summary.startswith(f"NOT SEARCHED ({cat})")
    assert neowise_search_status(out[0], 6)[0] == "searched"

    cov = search_coverage_summary(sample, cfg)
    assert cov["n_input"] == 7 and cov["n_searched"] == 3
    assert cov["n_no_coverage"] == 1 and cov["n_failed"] == 2 and cov["n_not_fetched"] == 1
    assert cov["not_searched_ids"]["failed"] == ["103", "998"]

    cert = ValidationCertificate(
        channel="ir_timedomain", gaia_release="dr3", split="all", config_hash="x", git_commit="x", timestamp_utc="x", passed=True,
        blind_recovery_rate=1.0, blind_total=1, blind_recovered=1, impostor_pass=True, impostor_max_tier="Explained",
        impostors_evaluated=1, injection_completeness_design_point=1.0, completeness_grid=[{"completeness": 1.0}],
        false_alarm_rate=0.0, false_alarm_ci_low=0.0, false_alarm_ci_high=0.01, expected_false_positives=0.0,
    )
    for s in sample:
        s.setdefault("distance_pc", 120.0)
        s.setdefault("teff_k", 5600.0)
    rep, jp, mp = write_null_report(
        channel="ir_timedomain", run_id="t", sample_stars=sample,
        sample_cuts=[SampleCutStep(cut_name="live", description="d", count_before=7, count_after=7)], healpix_cells=[1],
        n_flagged=1, tier_breakdown={"C": 0}, cert=cert, manifests=[], cfg=cfg, root=tmp_path,
    )
    assert rep.n_searched == 3
    assert {c.cut_name: c.count_before - c.count_after for c in rep.sample_cuts[1:]} == {
        "not_searched_failed": 2, "not_searched_no_coverage": 1, "not_searched_not_fetched": 1}
    assert all(ul["effective_n_stars"] <= 3 for ul in rep.upper_limits)
    assert sum(rep.distance_distribution.values()) == 3
    assert "Not Searched**: 4" in mp.read_text(encoding="utf-8")
    with pytest.raises(ValueError, match="absent from sample_stars"):
        write_null_report(channel="ir_timedomain", run_id="t2", sample_stars=sample[:2], sample_cuts=[], healpix_cells=[], n_flagged=0,
                          tier_breakdown={}, cert=cert, manifests=[], cfg=cfg, root=tmp_path, not_searched={"failed": ["103"]})


def test_concurrent_workers_give_the_same_stars_in_input_order(tmp_path: Path, cfg: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    """Several workers change how long the fetch takes, not what it returns."""
    import threading

    stars, fake = _survey()
    stars = [*stars, *[_star(str(200 + i), 160.0 + i, 22.0) for i in range(12)]]
    for st in stars[5:]:
        fake.behaviour[FakeIrsa.key(st["ra_deg"], st["dec_deg"])] = "empty"
    in_flight = {"now": 0, "max": 0}
    lock = threading.Lock()
    real_get = fake.get

    def slow_get(url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> Any:
        with lock:
            in_flight["now"] += 1
            in_flight["max"] = max(in_flight["max"], in_flight["now"])
        try:
            import time as _t

            _t.sleep(0.02)
            return real_get(url, params, timeout)
        finally:
            with lock:
                in_flight["now"] -= 1

    monkeypatch.setattr(requests, "get", slow_get)
    monkeypatch.setattr(requests, "post", fake.post)
    progress: list[int] = []
    par, sum_par = enrich_sample_neowise(stars, cfg, RawCache(tmp_path / "par"), sleep=lambda s: None, workers=4,
                                         progress=lambda done, *a: progress.append(done))
    assert sum_par["workers"] == 4 and in_flight["max"] > 1, "requests must overlap"
    assert progress == list(range(1, len(stars) + 1))

    fake2 = _survey()[1]
    fake2.behaviour.update(fake.behaviour)
    monkeypatch.setattr(requests, "get", fake2.get)
    monkeypatch.setattr(requests, "post", fake2.post)
    seq, sum_seq = enrich_sample_neowise(stars, cfg, RawCache(tmp_path / "seq"), sleep=lambda s: None, workers=1)
    assert [s["source_id"] for s in par] == [s["source_id"] for s in stars]

    def strip(star: dict[str, Any]) -> dict[str, Any]:
        d = json.loads(json.dumps(star))
        d["neowise_fetch"].pop("manifest_ids", None)
        st = d["neowise_fetch"]["status"]
        if st.startswith("failed"):  # the failure text carries no order-dependent detail, but keep the comparison to its class
            d["neowise_fetch"]["status"] = d["vetting_archive"]["query_status"]["neowise"] = "failed"
        return d

    assert [strip(a) for a in par] == [strip(b) for b in seq]
    for k in ("n_ok", "n_no_coverage", "n_failed", "failed_ids", "n_network_requests", "n_cache_hits"):
        assert sum_par[k] == sum_seq[k], k
    # Every manifest line written by the workers is whole and the cache replays the run without the network
    cache = RawCache(tmp_path / "par")
    assert not cache.recovered_truncated_manifest_line
    monkeypatch.setattr(requests, "get", lambda *a, **k: (_ for _ in ()).throw(requests.exceptions.ConnectionError("offline")))
    again, sum_again = enrich_sample_neowise(stars, cfg, cache, sleep=lambda s: None, workers=4)
    assert sum_again["n_ok"] == sum_par["n_ok"] and sum_again["n_no_coverage"] == sum_par["n_no_coverage"]


def test_pacer_spaces_request_starts_across_workers() -> None:
    from oddstar.io.neowise_bulk import _RequestPacer

    waits: list[float] = []
    pacer = _RequestPacer(0.5, waits.append)
    import time as _t

    t0 = _t.monotonic()
    for _ in range(3):  # three workers finish requests that all started at t0
        pacer.after_network_request(t0)
    assert len(waits) == 3 and 0.4 < waits[0] <= 0.5 and 0.9 < waits[1] <= 1.0 and 1.4 < waits[2] <= 1.5
