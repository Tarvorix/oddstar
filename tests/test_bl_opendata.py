"""Radio-pointing census against the Breakthrough Listen file index (offline; the index is a double)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest

from oddstar.errors import ArchiveAccessError
from oddstar.io import bl_opendata as bl
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache


class _Resp:
    def __init__(self, body: Any, status_code: int = 200) -> None:
        self.status_code = status_code
        self.content = body if isinstance(body, bytes) else json.dumps(body).encode("utf-8")


def _file(target: str, ra: float, dec: float, freq: float, tel: str = "GBT", size: int = 17_000_000_000, mjd: float = 57668.3,
          ftype: str = "HDF5") -> dict[str, Any]:
    return {"target": target, "ra": ra, "decl": dec, "center_freq": freq, "telescope": tel, "size": size, "mjd": mjd, "file_type": ftype,
            "url": f"https://bldata.berkeley.edu/x/{target}.h5", "md5sum": "0" * 32}


class _Index:
    def __init__(self, files: list[dict[str, Any]]) -> None:
        self.files = files
        self.calls: list[dict[str, Any]] = []
        self.fail_at: set[str] = set()

    def get(self, url: str, params: Any, **kw: Any) -> _Resp:
        self.calls.append(dict(params))
        assert url == bl.BL_QUERY_URL and params["target"] == ""
        ra, dec, rad = float(params["pos-ra"]), float(params["pos-dec"]), float(params["pos-rad"])
        if params["pos-ra"] in self.fail_at:
            raise ArchiveAccessError("index unreachable", retryable=True)
        hits = [f for f in self.files if ot.angular_separation_arcsec(ra, dec, f["ra"], f["decl"]) <= rad * 3600.0]
        return _Resp({"data": hits, "result": "success"})


def test_half_power_radius_follows_dish_and_frequency() -> None:
    assert bl.half_power_radius_arcsec("GBT", 1400.0) == pytest.approx(265.0, abs=5.0)       # 8.8' beam
    assert bl.half_power_radius_arcsec("Parkes", 1400.0) == pytest.approx(414.0, abs=8.0)    # 13.8' beam
    assert bl.half_power_radius_arcsec("GBT", 6000.0) == pytest.approx(61.9, abs=1.5)
    assert bl.half_power_radius_arcsec("APF", 5500.0) is None, "an optical telescope is not a radio pointing"
    assert bl.half_power_radius_arcsec("GBT", 0.0) is None


def test_star_is_covered_only_inside_the_half_power_beam() -> None:
    ra, dec = 55.0, -3.0
    rows = [
        _file("HIP1", ra, dec + 200.0 / 3600.0, 1400.0), _file("HIP1", ra, dec + 200.0 / 3600.0, 1420.0, size=264_000_000),
        _file("HIP1", ra, dec + 200.0 / 3600.0, 6000.0),           # same pointing, C band: 200" is outside a 62" radius
        _file("HIP2", ra, dec - 600.0 / 3600.0, 1400.0, tel="GBT"),  # 10' away: outside
        _file("HIP3", ra, dec + 5.0 / 3600.0, 5500.0, tel="APF"),    # optical spectrum, ignored
        {"target": "broken", "ra": None, "decl": dec, "center_freq": 1400.0, "telescope": "GBT"},
    ]
    cov = bl.beam_coverage(ra, dec, rows)
    assert cov["status"] == "in_beam" and cov["n_radio_files"] == 4 and cov["n_pointings"] == 3 and cov["n_pointings_in_beam"] == 1
    first = cov["pointings"][0]
    assert first["target"] == "HIP1" and first["band_mhz"] == 1400 and first["in_beam"] and first["n_files"] == 2
    assert first["total_bytes"] == 17_264_000_000 and first["sep_arcsec"] == pytest.approx(200.0, abs=0.5)
    assert [p["in_beam"] for p in cov["pointings"]] == [True, False, False]
    assert bl.beam_coverage(ra, dec, rows[2:])["status"] == "near_pointing"
    assert bl.beam_coverage(ra, dec, rows[4:])["status"] == "no_coverage"
    assert bl.beam_coverage(ra, dec, [])["status"] == "no_coverage"


def test_index_answers_are_cached_and_a_bad_answer_is_never_an_empty_result(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    idx = _Index([_file("HIP1", 10.0, 20.0, 1400.0)])
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    cache = RawCache(tmp_path / "raw")
    rows, man, hit = bl.query_bl_files(10.01, 20.0, cache)
    assert len(rows) == 1 and not hit and man is not None and idx.calls[0]["limit"] == str(bl.BL_ROW_LIMIT)
    assert bl.query_bl_files(10.01, 20.0, cache)[2] is True and len(idx.calls) == 1
    empty, _, _ = bl.query_bl_files(100.0, -40.0, cache)
    assert empty == [] and bl.query_bl_files(100.0, -40.0, cache)[2] is True, "a real empty answer is cached"

    for body in ({"result": "error", "data": []}, {"result": "success"}, b"<html>busy</html>"):
        monkeypatch.setattr(ot, "_get_with_retries", lambda url, params, _b=body, **kw: _Resp(_b))
        for _ in range(2):  # not stored: asked again each time
            with pytest.raises(ArchiveAccessError):
                bl.query_bl_files(200.0, 5.0, cache)


def test_census_counts_each_status_and_resumes(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    idx = _Index([_file("HIP1", 10.0, 20.0, 1400.0), _file("HIP2", 50.0, -10.0, 6000.0)])
    idx.fail_at = {"80.00000"}
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    stars = [
        {"source_id": "1", "ra_deg": 10.0 + 100.0 / 3600.0, "dec_deg": 20.0},   # inside the L-band beam
        {"source_id": "2", "ra_deg": 50.0, "dec_deg": -10.0 + 300.0 / 3600.0},  # 5' from a C-band pointing: near, not in beam
        {"source_id": "3", "ra_deg": 200.0, "dec_deg": 45.0},                   # nothing nearby
        {"source_id": "4", "ra_deg": 80.0, "dec_deg": 0.0},                     # index fails
    ]
    slept: list[float] = []
    seen: list[str] = []
    recs, summary = bl.sample_radio_coverage(stars, RawCache(tmp_path / "raw"), sleep=slept.append,
                                             progress=lambda d, t, sid, st: seen.append(f"{d}/{t}:{sid}:{st.split(':')[0]}"))
    assert [r["status"].split(":")[0] for r in recs] == ["in_beam", "near_pointing", "no_coverage", "failed"]
    assert summary["n_in_beam"] == 1 and summary["n_near_pointing"] == 1 and summary["n_no_coverage"] == 1 and summary["n_failed"] == 1
    assert summary["in_beam_ids"] == ["1"] and summary["near_pointing_ids"] == ["2"] and summary["failed_ids"] == ["4"]
    assert seen == ["1/4:1:in_beam", "2/4:2:near_pointing", "3/4:3:no_coverage", "4/4:4:failed"]
    assert len(slept) == 4, "every network request is followed by a pause"

    n = len(idx.calls)
    slept.clear()
    recs2, _ = bl.sample_radio_coverage(stars, RawCache(tmp_path / "raw"), sleep=slept.append)
    assert len(idx.calls) == n + 1 and len(slept) == 1, "only the failed star is asked for again"
    assert [r["status"].split(":")[0] for r in recs2] == ["in_beam", "near_pointing", "no_coverage", "failed"]
