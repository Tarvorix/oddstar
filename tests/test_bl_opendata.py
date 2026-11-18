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
    ra, dec = 55.0, 3.0
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

    crowded = {"result": "success", "data": [_file("HIP1", 10.0, 20.0, 1400.0)] * bl.BL_TRUSTED_ROWS}
    for body in ({"result": "error", "data": []}, {"result": "success"}, b"<html>busy</html>", crowded):
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
    assert len(slept) == len(idx.calls) == 8, "every network request is followed by a pause (the southern star takes five)"

    n = len(idx.calls)
    slept.clear()
    recs2, _ = bl.sample_radio_coverage(stars, RawCache(tmp_path / "raw"), sleep=slept.append)
    assert len(idx.calls) == n + 1 and len(slept) == 1, "only the failed star is asked for again"
    assert [r["status"].split(":")[0] for r in recs2] == ["in_beam", "near_pointing", "no_coverage", "failed"]


def _scan(target: str, ra: float, dec: float, mjd: float, freq: float = 1500.0, tel: str = "GBT", product: str = "0000",
          size: int = 15_000_000_000) -> dict[str, Any]:
    return {"target": target, "ra": ra, "decl": dec, "center_freq": freq, "telescope": tel, "size": size, "mjd": mjd,
            "file_type": "HDF5", "md5sum": "a" * 32,
            "url": f"https://bldata.berkeley.edu/x/blc00_guppi_{int(mjd)}_{int((mjd % 1) * 86400)}_{target}_{int(freq)}.gpuspec.{product}.h5"}


class _TimeIndex(_Index):
    """Answers both position queries and time-window queries."""

    def get(self, url: str, params: Any, **kw: Any) -> _Resp:
        if "time-start" not in params:
            return super().get(url, params, **kw)
        self.calls.append(dict(params))
        t0, t1 = float(params["time-start"]), float(params["time-end"])
        hits = [f for f in self.files if t0 <= f["mjd"] <= t1 and f["telescope"].lower() == params["telescopes"]]
        return _Resp({"data": hits, "result": "success"})


def test_scan_plan_pairs_each_on_scan_with_the_off_scans_around_it(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ra, dec, t = 9.7, 13.4, 57585.508
    step = 5.5 / 1440.0
    files = [
        _scan("HIPON", ra, dec + 200.0 / 3600.0, t), _scan("HIPON", ra, dec + 200.0 / 3600.0, t, product="0002", size=90_000_000),
        _scan("HIPB1", ra + 3.0, dec, t - step), _scan("HIPB2", ra - 4.0, dec + 1.0, t - 2 * step),
        _scan("HIPA1", ra + 5.0, dec - 2.0, t + step), _scan("HIPA2", ra + 6.0, dec, t + 2 * step),
        _scan("HIPCBAND", ra + 2.0, dec, t + step / 2, freq=6000.0),          # nearer in time but another receiver
        _scan("HIPCLOSE", ra, dec + 400.0 / 3600.0, t - step / 2),            # the star is too near this pointing to be an OFF
        _scan("HIPPKS", ra + 3.0, dec, t + step / 3, tel="Parkes"),
        _scan("HIPLATE", ra + 3.0, dec, t + 0.5),                              # hours later: not part of the block
    ]
    idx = _TimeIndex(files)
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    cache = RawCache(tmp_path / "raw")
    rows, _, _ = bl.query_bl_files(ra, dec, cache)
    plans = bl.plan_star_scans(ra, dec, rows, cache)
    assert len(plans) == 1
    p = plans[0]
    assert p["target"] == "HIPON" and p["complete"] and p["center_freqs_mhz"] == [1500.0] and p["sep_arcsec"] == pytest.approx(200.0, abs=0.5)
    assert len(p["on_files"]) == 1 and p["bytes_on"] == 15_000_000_000, "only the fine-resolution product is planned"
    assert p["off_before"]["target"] == "HIPB1" and p["off_before"]["minutes_from_on"] == pytest.approx(-5.5, abs=0.1)
    assert p["off_after"]["target"] == "HIPA1" and p["off_after"]["minutes_from_on"] == pytest.approx(5.5, abs=0.1)
    assert p["bytes_on_plus_offs"] == 45_000_000_000 and p["comparison"] == "before_after" and p["off_simultaneous"] == []
    n = len(idx.calls)
    assert bl.plan_star_scans(ra, dec, rows, cache)[0] == p and len(idx.calls) == n, "the time-window answer is cached"

    # without a later scan on another target the plan says so instead of inventing a comparison
    idx2 = _TimeIndex([f for f in files if f["target"] not in ("HIPA1", "HIPA2")])
    monkeypatch.setattr(ot, "_get_with_retries", idx2.get)
    cache2 = RawCache(tmp_path / "raw2")
    rows2, _, _ = bl.query_bl_files(ra, dec, cache2)
    p2 = bl.plan_star_scans(ra, dec, rows2, cache2)[0]
    assert not p2["complete"] and p2["comparison"] == "incomplete" and p2["off_after"] is None and p2["off_before"]["target"] == "HIPB1"
    assert p2["bytes_on_plus_offs"] == 30_000_000_000
    assert bl.plan_star_scans(200.0, 45.0, [], cache2) == []


def test_on_scans_that_follow_each_other_become_one_cadence(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A, B, A, C, A, D on one target: one plan with the middle ON scan as reference, the other two and all three OFF scans."""
    ra, dec, t = 183.0, 13.26, 59426.8885
    step = 5.3 / 1440.0
    files = []
    for k, (name, d_ra) in enumerate((("NGCX", 0.0), ("HIPB", 1.0), ("NGCX", 0.0), ("HIPC", -0.3), ("NGCX", 0.0), ("HIPD", 0.5))):
        on = name == "NGCX"
        for freq in (845.2, 1032.7):
            files.append(_scan(name, ra + (0.05 if on else d_ra), dec + (0.0 if on else 1.2), t + k * step, freq=freq, size=3_200_000_000))
    files.append(_scan("NGCX", ra + 0.05, dec, t + 0.5, freq=845.2))       # half a day later: its own scan, and nothing to compare it with
    idx = _TimeIndex(files)
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    cache = RawCache(tmp_path / "raw")
    rows, _, _ = bl.query_bl_files(ra, dec, cache)
    plans = bl.plan_star_scans(ra, dec, rows, cache)
    assert [(p["comparison"], p["complete"]) for p in plans] == [("cadence", True), ("incomplete", False)]
    p = plans[0]
    assert p["n_on_scans"] == 3 and p["mjd"] == pytest.approx(t + 2 * step) and p["center_freqs_mhz"] == [845.2, 1032.7] and len(p["on_files"]) == 2
    assert [e["minutes_from_on"] for e in p["on_other"]] == pytest.approx([-10.6, 10.6], abs=0.05) and all(len(e["files"]) == 2 for e in p["on_other"])
    assert [(e["target"], e["minutes_from_on"]) for e in p["off_cadence"]] == [("HIPB", pytest.approx(-5.3, abs=0.05)), ("HIPC", pytest.approx(5.3, abs=0.05)),
                                                                              ("HIPD", pytest.approx(15.9, abs=0.05))]
    assert p["off_before"] is None and p["off_after"] is None and p["bytes_on_plus_offs"] == 12 * 3_200_000_000
    assert "on_other" not in plans[1]


def test_on_and_off_files_are_of_the_same_product(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The archive holds two fine products of different channel width (.0000.h5 and .0000.hires.hdf); they are never paired."""
    ra, dec, t = 150.0, 20.0, 58225.6
    step = 5.0 / 1440.0

    def hires(row: dict[str, Any]) -> dict[str, Any]:
        return dict(row, url=row["url"].replace(".gpuspec.0000.h5", ".0000.hires.hdf"), size=500_000_000)

    assert bl._product("http://x/a_0001.0000.hires.hdf") == "hires" and bl._product("http://x/a_0001.0000.h5") == "fine"
    # Scan 1: the ON scan exists as hires only, its neighbours as h5 only: nothing to compare it with
    files = [hires(_scan("ONA", ra + 0.02, dec, t)), _scan("OFF1", ra + 2.0, dec, t - step), _scan("OFF2", ra + 2.0, dec, t + step)]
    # Scan 2, an hour later: ON as hires only, neighbours in both products: the hires neighbours are used
    t2 = t + 1.0 / 24.0
    files += [hires(_scan("ONA", ra + 0.02, dec, t2))]
    for name, dt in (("OFF3", -step), ("OFF4", step)):
        row = _scan(name, ra + 2.0, dec, t2 + dt)
        files += [row, hires(row)]
    # Scan 3, another hour later: ON and neighbours in both products: searched once, as h5
    t3 = t + 2.0 / 24.0
    for name, d_ra, dt in (("ONA", 0.02, 0.0), ("OFF5", 2.0, -step), ("OFF6", 2.0, step)):
        row = _scan(name, ra + d_ra, dec, t3 + dt)
        files += [row, hires(row)]
    idx = _TimeIndex(files)
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    cache = RawCache(tmp_path / "raw")
    rows, _, _ = bl.query_bl_files(ra, dec, cache)
    plans = bl.plan_star_scans(ra, dec, rows, cache)
    assert [(p["comparison"], p["complete"]) for p in plans] == [("incomplete", False), ("before_after", True), ("before_after", True)]
    assert plans[0]["off_before"] is None and plans[0]["off_after"] is None
    for p, product in ((plans[1], "hires"), (plans[2], "fine")):
        urls = [f["url"] for e in (p, p["off_before"], p["off_after"]) for f in (e["on_files"] if e is p else e["files"])]
        assert len(urls) == 3 and {bl._product(u) for u in urls} == {product}


def test_fine_resolution_files_are_recognised_for_each_recorder() -> None:
    fine = [
        "http://x/blc00_guppi_57585_43906_HIP3026_0022.gpuspec.0000.h5", "http://x/blc00_guppi_57585_43906_HIP3026_0022.gpuspec.0000.fil",
        "http://x/collate_mb_2/PKS_0295/blc23/blc23_guppi_58233_09671_008272_G239.52-51.04_0001.0000.hires.hdf",
        "http://x/collate_mb/PKS_0338/blc04/guppi_58336_70277_045474_G306.02-45.34_0001.0000.h5",
        "http://x/PKSUWL/blc01/PKSUWL/blc01_HIP1837_S_2019-05-19T19:11:11.000.fil",
    ]
    coarse = [
        "http://x/blc00_guppi_57585_43906_HIP3026_0022.gpuspec.0002.h5", "http://x/a_0001.0000.midres.hdf", "http://x/a_0001.0000.midres.stokes.hdf",
        "http://x/a_0001.0001.fil", "http://x/other/blc01_HIP1837_S_2019-05-19T19:11:11.000.fil",
    ]
    assert all(bl.FINE_RESOLUTION_URL.search(u) for u in fine)
    assert not any(bl.FINE_RESOLUTION_URL.search(u) for u in coarse)


def test_scan_plan_uses_beams_of_the_same_instant_and_one_file_per_frequency(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    ra, dec, t = 57.56, -36.46, 58233.1038
    sec = 1.0 / 86400.0

    def beam(name: str, d_ra: float, d_dec: float, mjd: float = t) -> list[dict[str, Any]]:
        out = []
        for i, freq in enumerate((1282.75, 1436.75)):
            f = _scan(name, ra + d_ra, dec + d_dec, mjd + i * sec, freq=freq, tel="Parkes", size=2_768_000_000)  # bands start 1 s apart
            f["url"] = f"http://x/collate_mb_2/blc{i}/{name}_0001.0000.hires.hdf"
            out.append(f)
        dup = dict(out[0], url=out[0]["url"].replace(".0000.hires.hdf", ".0000.fil"), size=3_000_000_000)   # same product as filterbank
        mid = dict(out[0], url=out[0]["url"].replace(".hires.", ".midres."), size=90_000_000)
        return out + [dup, mid]

    files = (beam("G239.40-51.24", 0.0, -0.05) + beam("G238.90-50.90", 1.0, 0.3) + beam("G239.90-51.60", -1.1, -0.2)
             + beam("G240.50-52.00", 2.5, 1.0) + beam("G239.00-51.30", 0.3, 0.0)                  # 14 arcmin away: too near to be an OFF
             + beam("G250.00-50.00", 3.0, 0.0, mjd=t + 6.0 / 1440.0))
    idx = _TimeIndex(files)
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    cache = RawCache(tmp_path / "raw")
    rows, _, _ = bl.query_bl_files(ra, dec, cache)
    plans = bl.plan_star_scans(ra, dec, rows, cache)
    assert len(plans) == 1, "two bands that start a second apart are one scan"
    p = plans[0]
    assert p["target"] == "G239.40-51.24" and p["center_freqs_mhz"] == pytest.approx([1282.75, 1436.75], abs=0.06)
    assert len(p["on_files"]) == 2 and all(f["url"].endswith(".0000.hires.hdf") for f in p["on_files"]) and p["bytes_on"] == 2 * 2_768_000_000
    assert p["comparison"] == "simultaneous_beams" and p["complete"]
    assert [e["target"] for e in p["off_simultaneous"]] == ["G238.90-50.90", "G239.90-51.60"], "the two nearest beams beyond the exclusion radius"
    assert all(e["minutes_from_on"] == 0.0 and len(e["files"]) == 2 for e in p["off_simultaneous"])
    assert p["off_before"] is None and p["off_after"]["target"] == "G250.00-50.00"
    assert p["bytes_on_plus_offs"] == 6 * 2_768_000_000


GBT_URL = "https://bldata.berkeley.edu/p/spliced_blc0001020304050607_guppi_57585_43576_{}_0021.gpuspec.0000.h5"
UWL_URL = "http://blpd2.ssl.berkeley.edu/PKSUWL/blc01/PKSUWL/blc01_{}_2019-05-19T19:11:11.000.fil"
PKS_FINE_URL = "http://blpd8.ssl.berkeley.edu/dl/Parkes_57840_11779_{}_fine.h5"
MULTIBEAM_URL = "http://blpd12.ssl.berkeley.edu/collate_mb_2/PKS_0295/blc10/blc10_guppi_58233_07909_007538_{}_0001.0000.hires.hdf"


def test_southern_index_declinations_are_put_right() -> None:
    # index values and catalogue positions measured 2026-11-18 (see the comment above DEC_MIRROR_40_40)
    for tel, url, listed, true, rule in [
        ("GBT", GBT_URL.format("Hip3909"), -10.0325006, -10.64528, "mirror_40_40"),
        ("GBT", GBT_URL.format("Hip3026"), -8.3698361111, -8.30794, "mirror_40_40"),
        ("GBT", GBT_URL.format("HIP17147"), -3.4600025, -3.21778, "mirror_40_40"),
        ("Parkes", UWL_URL.format("HIP1837_S"), -33.51083, -33.16695, "mirror_40_40"),
        ("Parkes", PKS_FINE_URL.format("HIP113368_S"), -28.3777777778, -29.62222, "mirror_0"),
        ("Parkes", PKS_FINE_URL.format("HIP8102_S"), -14.0625, -15.9375, "mirror_0"),
        ("Parkes", MULTIBEAM_URL.format("G238.23-51.24"), -35.93489, -35.93489, "as_listed"),
        ("Parkes", "http://blpd0.ssl.berkeley.edu/other/a_0001.0000.fil", -40.5, -40.5, "unverified"),
        ("GBT", "https://bldata.berkeley.edu/other/a.0000.h5", -40.5, -40.5, "unverified"),
        ("GBT", "https://bldata.berkeley.edu/other/a.0000.h5", 40.5, 40.5, "as_listed"),
        ("GBT", GBT_URL.format("HIP1"), 12.25, 12.25, "as_listed"),
    ]:
        dec, how = bl.true_declination(tel, url, listed)
        assert how == rule and dec == pytest.approx(true, abs=2e-5), url
    for c in bl.DEC_MIRRORS.values():
        for dec in (-1.0, -1.5, -20.0, -20.999999, -45.677778, -45.7, -89.3):
            listed = bl.mirror_declination(dec, c)
            assert listed is not None and bl.unmirror_declination(listed, c) == pytest.approx(dec, abs=1e-9)
        assert bl.mirror_declination(-0.5, c) is None and bl.mirror_declination(3.0, c) is None
    assert bl.candidate_declinations(-14.0625) == pytest.approx([-14.0625, -14.615278, -15.9375], abs=1e-6)


def test_a_southern_position_is_asked_for_wherever_the_index_may_list_it() -> None:
    assert bl.index_query_declinations(20.0, 0.25) == [20.0]
    assert bl.index_query_declinations(-10.5, 0.25) == pytest.approx([-10.5, -(10 + bl.DEC_MIRROR_40_40 - 0.5), -9.5])
    # the circle reaches into the ninth degree, which is mirrored somewhere else
    assert bl.index_query_declinations(-10.1, 0.25) == pytest.approx(
        [-10.1, -(10 + bl.DEC_MIRROR_40_40 - 0.1), -(9 + bl.DEC_MIRROR_40_40 - 0.85), -9.9, -8.15])
    assert bl.index_query_declinations(-0.5, 0.25) == [-0.5], "the strip above -1 degree cannot be recovered"


def test_southern_stars_are_found_through_the_mirrored_declination(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    def listed(dec: float, rule: str) -> float:
        m = bl.mirror_declination(dec, bl.DEC_MIRRORS[rule])
        assert m is not None
        return m

    def rec(target: str, url: str, ra: float, decl: float, freq: float, tel: str, n: int) -> dict[str, Any]:
        return {"id": n, "target": target, "ra": ra, "decl": decl, "center_freq": freq, "telescope": tel, "size": 10**9, "mjd": 57585.5,
                "file_type": "HDF5", "md5sum": "0" * 32, "url": url}

    idx = _Index([
        rec("Hip3026", GBT_URL.format("Hip3026"), 9.6329, listed(-8.30794, "mirror_40_40"), 1475.0, "GBT", 1),
        rec("HIP1837_S", UWL_URL.format("HIP1837_S"), 5.81625, listed(-33.16695, "mirror_40_40"), 768.0, "Parkes", 2),
        rec("HIP8102_S", PKS_FINE_URL.format("HIP8102_S"), 26.017, listed(-15.9375, "mirror_0"), 3094.0, "Parkes", 3),
        rec("G1", MULTIBEAM_URL.format("G1"), 57.56, -35.93, 1282.75, "Parkes", 4),
        rec("odd", "http://blpd0.ssl.berkeley.edu/other/odd_0001.0000.fil", 100.0, listed(-50.3, "mirror_0"), 1400.0, "Parkes", 5),
        # a Green Bank file listed where a star truly is, but truly 0.6 degree away
        rec("HIP9", GBT_URL.format("HIP9"), 150.0, -20.2, 1475.0, "GBT", 6),
    ])
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    stars = [
        {"source_id": "hip3026", "ra_deg": 9.6332, "dec_deg": -8.3117},
        {"source_id": "hip1837", "ra_deg": 5.8163, "dec_deg": -33.1670},
        {"source_id": "hip8102", "ra_deg": 26.017, "dec_deg": -15.9375},
        {"source_id": "multibeam", "ra_deg": 57.56, "dec_deg": -35.95},
        {"source_id": "odd", "ra_deg": 100.0, "dec_deg": -50.3},
        {"source_id": "coincidence", "ra_deg": 150.0, "dec_deg": -20.2},
    ]
    recs, summary = bl.sample_radio_coverage(stars, RawCache(tmp_path / "raw"), sleep=lambda s: None)
    by = {r["source_id"]: r for r in recs}
    assert [by[k]["status"] for k in ("hip3026", "hip1837", "hip8102", "multibeam")] == ["in_beam"] * 4
    assert by["hip3026"]["pointings"][0]["sep_arcsec"] == pytest.approx(13.6, abs=0.5) and by["hip3026"]["pointings"][0]["dec_rules"] == ["mirror_40_40"]
    assert by["odd"]["status"] == "unverified_pointing" and by["odd"]["pointings"][0]["possibly_in_beam"] and not by["odd"]["pointings"][0]["in_beam"]
    assert by["coincidence"]["status"] == "near_pointing", "the listed declination of a mirrored file is not where the telescope pointed"
    assert summary["n_in_beam"] == 4 and summary["n_unverified_pointing"] == 1 and summary["unverified_pointing_ids"] == ["odd"]
    assert len(by["hip3026"]["manifest_ids"]) == 3


def test_scan_plan_of_a_southern_star_uses_true_declinations(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    c = bl.DEC_MIRROR_40_40
    ra, dec, t = 9.6332, -8.3117, 57585.508
    files = [_scan("Hip3026", 9.6329, bl.mirror_declination(-8.30794, c), t), _scan("Hip3909", 12.53, bl.mirror_declination(-10.64528, c), t - 0.004),
             _scan("Hip3909", 12.53, bl.mirror_declination(-10.64528, c), t + 0.004)]
    idx = _TimeIndex(files)
    monkeypatch.setattr(ot, "_get_with_retries", idx.get)
    cache = RawCache(tmp_path / "raw")
    rows, _, _ = bl.query_star_rows(ra, dec, cache)
    plans = bl.plan_star_scans(ra, dec, rows, cache)
    assert len(plans) == 1 and plans[0]["comparison"] == "before_after" and plans[0]["sep_arcsec"] == pytest.approx(13.6, abs=0.5)
    assert plans[0]["off_before"]["target"] == "Hip3909" and plans[0]["off_before"]["sep_arcsec"] > 10000.0
