"""Offline tests for live archive mirror fallback, Gator/VizieR parsers, and MAST PS1 filtering."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pytest
import requests

from oddstar.config import load_config
from oddstar.errors import ArchiveAccessError
from oddstar.io import live_survey
from oddstar.io.cache import RawCache
from oddstar.io.live_survey import (
    _parse_ipac_table,
    _parse_vizier_tsv,
    _tap_query_cached,
    enrich_star_live_for_vetting,
)


class _DummyResp:
    def __init__(self, status_code: int, text: str = "", json_data: dict[str, Any] | None = None) -> None:
        self.status_code = status_code
        self.text = text
        self._json = json_data

    def json(self) -> dict[str, Any]:
        assert self._json is not None
        return self._json


def test_tap_query_mirror_fallback_and_cache(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    calls: list[str] = []

    def fake_post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        calls.append(url)
        if "primary.example" in url:
            raise requests.exceptions.ConnectTimeout("primary down")
        return _DummyResp(200, json_data={"metadata": [{"name": "id"}], "data": [[12345]]})

    monkeypatch.setattr(requests, "post", fake_post)
    rows, _, hit1 = _tap_query_cached(
        urls=["https://primary.example/sync", "https://backup.example/sync"],
        query="SELECT 1",
        archive="simbad",
        release="live",
        cache=cache,
    )
    assert hit1 is False
    assert rows == [{"id": 12345}]
    assert calls == ["https://primary.example/sync", "https://backup.example/sync"]

    # Second call uses RawCache with zero network
    rows2, _, hit2 = _tap_query_cached(
        urls=["https://primary.example/sync", "https://backup.example/sync"],
        query="SELECT 1",
        archive="simbad",
        release="live",
        cache=cache,
    )
    assert hit2 is True
    assert rows2 == [{"id": 12345}]
    assert len(calls) == 2


def test_parse_ipac_and_vizier_tsv() -> None:
    ipac = (
        "\\fixlen = T\n"
        "|      ra|     dec|w4mpro|cc_flags|\n"
        "|  double|  double|double|    char|\n"
        "|     deg|     deg|   mag|        |\n"
        "|    null|    null|  null|    null|\n"
        " 43.54827  5.88271  7.290     0000\n"
        " 43.54830  5.88280   null     0000\n"
    )
    rows = _parse_ipac_table(ipac)
    assert len(rows) == 2
    assert rows[0]["w4mpro"] == "7.290"
    assert rows[1]["w4mpro"] is None

    tsv = (
        "# VizieR\n"
        "RAJ2000\tDEJ2000\tS1_4\n"
        "deg\tdeg\tmJy\n"
        "-------\t-------\t---\n"
        "43.5482\t5.8827\t12.4\n"
    )
    vrows = _parse_vizier_tsv(tsv)
    assert vrows == [{"RAJ2000": "43.5482", "DEJ2000": "5.8827", "S1_4": "12.4"}]


def test_mast_ps1_filters_single_detection_ghosts(tmp_path: Path, monkeypatch: Any) -> None:
    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    star: dict[str, Any] = {
        "source_id": "5628404318169472",
        "name": "Gaia DR3 5628404318169472",
        "ra_deg": 43.54827,
        "dec_deg": 5.88271,
        "mags": {"G": 13.43, "W4": 7.29},
        "flux_mjy": {"W4": 10.15},
        "allwise": {"ra_deg": 43.5482, "dec_deg": 5.8828},
        "neowise": {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []},
        "lightcurves": {"ztf_flux": []},
        "vetting_archive": {"radio_sources_within_10arcsec": [], "deep_imaging_available": False, "deep_imaging_offset_galaxy": None},
        "wise_single_frames_w4_mjy": [],
        "live_enriched": False,
    }

    def fake_post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        return _DummyResp(200, json_data={"metadata": [], "data": []})

    def fake_get(url: str, params: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        if "panstarrs" in url:
            return _DummyResp(
                200,
                json_data={
                    "info": [{"name": "raMean"}, {"name": "decMean"}, {"name": "nDetections"}, {"name": "rMeanPSFMag"}, {"name": "rMeanKronMag"}],
                    "data": [
                        [43.54827, 5.88271, 67, 13.476, 13.477],
                        [43.54950, 5.88320, 1, 21.1, -999.0],
                    ],
                },
            )
        if "vizier" in url:
            return _DummyResp(200, text="# VizieR empty\n")
        return _DummyResp(404, text="not found")

    monkeypatch.setattr(requests, "post", fake_post)
    monkeypatch.setattr(requests, "get", fake_get)

    enriched, _ = enrich_star_live_for_vetting(star, cfg, cache)
    va = enriched["vetting_archive"]
    # The PS1 catalogue is recorded as a neighbour screen only; it no longer sets a "deep imaging available" pass.
    assert va["query_status"]["deep_optical_ps1"] == "ok"
    screen = va["neighbour_screen"]
    assert screen["surveys"] and screen["surveys"][0].startswith("PS1 DR2")
    assert screen["contaminant"] is None
    # The nDetections=1 ghost ~4.8" away is filtered out, not recorded as a faint point source.
    assert screen["faint_point_sources"] == [] and screen["extended_neighbours"] == []
    assert "deep_imaging_available" not in va


def test_live_region_queries_halo_and_bulge(tmp_path: Path, monkeypatch: Any) -> None:
    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    seen_queries: list[str] = []

    def fake_post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        if data and "QUERY" in data:
            seen_queries.append(data["QUERY"])
        return _DummyResp(200, json_data={"metadata": [], "data": []})

    monkeypatch.setattr(requests, "post", fake_post)
    from oddstar.io.live_survey import fetch_live_halo_bulge_sample

    fetch_live_halo_bulge_sample(cfg, cache, limit=100, region="halo", after_random_index=999)
    fetch_live_halo_bulge_sample(cfg, cache, limit=100, region="bulge")
    # Halo and bulge go straight to the Gaia-only index query (the single-shot join timed out live on every attempt);
    # the cursor is applied to the index rows, not in SQL.
    assert len(seen_queries) == 2
    assert all(q.startswith("SELECT g.source_id, g.random_index FROM gaiadr3.gaia_source") for q in seen_queries)
    assert "ABS(g.b) > 20.0" in seen_queries[0]
    assert "g.mh_gspphot < -0.8" in seen_queries[0]
    assert "g.l <= 12.0 OR g.l >= 348.0" in seen_queries[1]
    assert "ABS(g.b) BETWEEN 3.0 AND 10.0" in seen_queries[1]


def test_live_region_uses_a_cached_single_shot_answer(tmp_path: Path, monkeypatch: Any) -> None:
    """A single-shot sample stored by an earlier run still defines that sample; nothing is re-queried."""
    from oddstar.io import live_survey as ls

    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    queries: list[str] = []

    def record(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        queries.append((data or {}).get("QUERY", ""))
        return _DummyResp(200, json_data={"metadata": [], "data": []})

    monkeypatch.setattr(requests, "post", record)
    monkeypatch.setattr(ls, "_fetch_live_two_stage", lambda **kw: ([], []))
    real_lookup = cache.lookup
    seen_sigs: list[tuple[str, str, str]] = []

    def spy(archive: str, release: str, sig: str) -> Any:
        seen_sigs.append((archive, release, sig))
        return real_lookup(archive, release, sig)

    monkeypatch.setattr(cache, "lookup", spy)
    ls.fetch_live_halo_bulge_sample(cfg, cache, limit=100, region="halo")
    archive, release, sig = next(s for s in seen_sigs if s[1] == f"{cfg.gaia.release}_halo")
    assert queries == []  # no single-shot attempt when nothing is cached
    cache.put(archive=archive, release=release, request_signature=sig, query_text="single-shot halo", payload=b"[]", row_count=0,
              offline_fixture=False)
    monkeypatch.setattr(ls, "_fetch_live_two_stage", lambda **kw: (_ for _ in ()).throw(AssertionError("two-stage must not run")))
    stars, mans = ls.fetch_live_halo_bulge_sample(cfg, cache, limit=100, region="halo")
    assert stars == [] and len(mans) == 1 and queries == []


def test_tap_busy_service_is_retried_after_waits(tmp_path: Path, monkeypatch: Any) -> None:
    """HTTP 429 means "queue full": the round is retried after each configured wait, then the failure says so."""
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    state = {"n": 0}
    slept: list[float] = []

    def busy_then_ok(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        state["n"] += 1
        if state["n"] <= 2:
            return _DummyResp(429)
        return _DummyResp(200, json_data={"metadata": [{"name": "source_id"}], "data": [[7]]})

    monkeypatch.setattr(requests, "post", busy_then_ok)
    rows, _, hit = _tap_query_cached(urls=["https://tap.example/sync"], query="SELECT 1", archive="gaia", release="dr3", cache=cache,
                                     busy_waits_s=(5.0, 9.0, 13.0), sleep=slept.append)
    assert rows == [{"source_id": 7}] and hit is False and slept == [5.0, 9.0]

    monkeypatch.setattr(requests, "post", lambda url, data=None, timeout=15: _DummyResp(429))
    slept.clear()
    with pytest.raises(ArchiveAccessError, match=r"after 2 wait\(s\) on a busy service.*HTTP 429"):
        _tap_query_cached(urls=["https://tap.example/sync"], query="SELECT 2", archive="gaia", release="dr3", cache=cache,
                          busy_waits_s=(5.0, 9.0), sleep=slept.append)
    assert slept == [5.0, 9.0]

    # A plain failure (not "busy") is not waited on
    monkeypatch.setattr(requests, "post", lambda url, data=None, timeout=15: _DummyResp(500))
    slept.clear()
    with pytest.raises(ArchiveAccessError):
        _tap_query_cached(urls=["https://tap.example/sync"], query="SELECT 3", archive="gaia", release="dr3", cache=cache,
                          busy_waits_s=(5.0,), sleep=slept.append)
    assert slept == []


class _MirrorResp:
    def __init__(self, status_code: int, content: bytes, content_type: str) -> None:
        self.status_code, self.content, self.headers = status_code, content, {"Content-Type": content_type}
        self.text = content.decode("utf-8")

    def json(self) -> Any:
        return json.loads(self.content)


_AIP_VOTABLE = b"""<?xml version="1.0" encoding="UTF-8"?>
<VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE type="results">
<INFO name="QUERY_STATUS" value="OK"/>
<TABLE>
<FIELD name="source_id" ucd="meta.id;meta.main;meta.ref" ID="datalinkID" datatype="long"/>
<FIELD name="ra" unit="deg" datatype="double"/>
<FIELD name="phot_g_mean_mag" unit="mag" datatype="float"/>
<FIELD name="pmra" datatype="double"/>
<FIELD name="flag" datatype="boolean"/>
<FIELD name="designation" datatype="char" arraysize="*"/>
<DATA><TABLEDATA>
<TR><TD>1697778723611368320</TD><TD>234.75674234671334</TD><TD>15.739729</TD><TD></TD><TD>true</TD><TD>Gaia DR3 1</TD></TR>
<TR><TD>1697778341358555776</TD><TD>234.73681651867895</TD><TD>19.423586</TD><TD>NaN</TD><TD>f</TD><TD>Gaia DR3 2</TD></TR>
</TABLEDATA></DATA>
</TABLE></RESOURCE></VOTABLE>"""


def test_gaia_mirrors_answer_when_the_esa_archive_is_in_maintenance(tmp_path: Path, monkeypatch: Any) -> None:
    """Seen 2026-11-18: ESA answered HTTP 200 with its downtime page; Heidelberg writes JSON rows as objects; Potsdam
    refuses FORMAT=json and answers FORMAT=votable with TABLEDATA."""
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    esa, ari, aip = live_survey.GAIA_TAP_URLS
    asked: list[tuple[str, str]] = []
    ari_rows = [{"source_id": 1697766384169254400, "ra": 234.95209613797536, "phot_g_mean_mag": 20.489115},
                {"source_id": 1697766388465292800, "ra": 234.93439619221073, "phot_g_mean_mag": None}]
    ari_body = json.dumps({"metadata": [{"name": "source_id", "datatype": "LONG"}, {"name": "ra"}, {"name": "phot_g_mean_mag"}],
                           "data": ari_rows}).encode("utf-8")

    def post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _MirrorResp:
        asked.append((url, (data or {})["FORMAT"]))
        if url == esa:
            return _MirrorResp(200, b"<html><body>Maintenance ongoing!</body></html>", "text/html; charset=UTF-8")
        if url == ari:
            return state["ari"]
        if (data or {})["FORMAT"] != "votable":
            return _MirrorResp(400, b"<VOTABLE/>", "application/xml; charset=utf-8")
        return state["aip"]

    state = {"ari": _MirrorResp(200, ari_body, "application/json;charset=UTF-8"), "aip": _MirrorResp(200, _AIP_VOTABLE, "application/xml")}
    monkeypatch.setattr(requests, "post", post)
    rows, _, hit = _tap_query_cached(urls=live_survey.GAIA_TAP_URLS, query="SELECT a", archive="gaia", release="dr3", cache=cache)
    assert rows == ari_rows and hit is False and asked == [(esa, "json"), (ari, "json")]

    # Heidelberg down too: Potsdam is asked for a VOTable, and its cells are typed by the FIELD list
    state["ari"] = _MirrorResp(500, b"", "text/plain")
    asked.clear()
    rows, _, _ = _tap_query_cached(urls=live_survey.GAIA_TAP_URLS, query="SELECT b", archive="gaia", release="dr3", cache=cache)
    assert asked == [(esa, "json"), (ari, "json"), (aip, "votable")]
    assert rows == [
        {"source_id": 1697778723611368320, "ra": 234.75674234671334, "phot_g_mean_mag": 15.739729, "pmra": None, "flag": True, "designation": "Gaia DR3 1"},
        {"source_id": 1697778341358555776, "ra": 234.73681651867895, "phot_g_mean_mag": 19.423586, "pmra": None, "flag": False, "designation": "Gaia DR3 2"},
    ]

    # Nothing answers: the error names the maintenance page instead of a JSON decoding error
    state["aip"] = _MirrorResp(200, _AIP_VOTABLE.replace(b'value="OK"', b'value="OVERFLOW"'), "application/xml")
    with pytest.raises(ArchiveAccessError) as err:
        _tap_query_cached(urls=live_survey.GAIA_TAP_URLS, query="SELECT c", archive="gaia", release="dr3", cache=cache)
    # (the ESA archive has failed twice by now, so it is asked after the mirrors and named last)
    for part in ("an HTML page instead of a TAP answer", "HTTP 500", "QUERY_STATUS=OVERFLOW"):
        assert part in str(err.value)

    # An unreadable answer is never an empty one
    for bad in (b"not xml", _AIP_VOTABLE.replace(b"<TD>true</TD>", b""), _AIP_VOTABLE.replace(b"TABLEDATA", b"BINARY2"),
                b"<VOTABLE><RESOURCE><INFO name='QUERY_STATUS' value='OK'/></RESOURCE></VOTABLE>"):
        with pytest.raises(ValueError):
            live_survey._rows_from_tap_votable(bad)
    with pytest.raises(ValueError):
        live_survey._rows_from_tap_json({"metadata": [{"name": "a"}, {"name": "b"}], "data": [{"a": 1}]})
    with pytest.raises(ValueError):
        live_survey._rows_from_tap_json({"metadata": [{"name": "a"}, {"name": "b"}], "data": [[1]]})
    empty = _AIP_VOTABLE[:_AIP_VOTABLE.index(b"<TR>")] + _AIP_VOTABLE[_AIP_VOTABLE.index(b"</TABLEDATA>"):]
    assert live_survey._rows_from_tap_votable(empty) == []


def test_live_region_two_stage_fallback_on_timeout(tmp_path: Path, monkeypatch: Any) -> None:
    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    seen_queries: list[str] = []

    def fake_post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        q = (data or {}).get("QUERY", "")
        seen_queries.append(q)
        if "SELECT TOP" in q:
            return _DummyResp(408)
        if "SELECT g.source_id, g.random_index FROM gaiadr3.gaia_source" in q:
            return _DummyResp(
                200,
                json_data={
                    "metadata": [{"name": "source_id"}, {"name": "random_index"}],
                    "data": [[111, 500], [222, 1200], [333, 1500]],
                },
            )
        if "WHERE g.source_id IN (222,333)" in q:
            return _DummyResp(
                200,
                json_data={
                    "metadata": [
                        {"name": "source_id"}, {"name": "ra"}, {"name": "dec"}, {"name": "l"}, {"name": "b"},
                        {"name": "parallax"}, {"name": "parallax_error"}, {"name": "parallax_over_error"},
                        {"name": "pmra"}, {"name": "pmdec"}, {"name": "ruwe"},
                        {"name": "phot_g_mean_mag"}, {"name": "phot_bp_mean_mag"}, {"name": "phot_rp_mean_mag"},
                        {"name": "teff_gspphot"}, {"name": "logg_gspphot"}, {"name": "mh_gspphot"},
                        {"name": "ebpminrp_gspphot"}, {"name": "random_index"},
                        {"name": "wise_designation"}, {"name": "w1mpro"}, {"name": "w1mpro_error"},
                        {"name": "w2mpro"}, {"name": "w2mpro_error"}, {"name": "w3mpro"}, {"name": "w3mpro_error"},
                        {"name": "w4mpro"}, {"name": "w4mpro_error"}, {"name": "cc_flags"}, {"name": "ext_flag"},
                        {"name": "ph_qual"}, {"name": "wise_ra"}, {"name": "wise_dec"},
                        {"name": "j_m"}, {"name": "j_msigcom"}, {"name": "h_m"}, {"name": "h_msigcom"},
                        {"name": "ks_m"}, {"name": "ks_msigcom"}, {"name": "tmass_ph_qual"},
                    ],
                    "data": [
                        [222, 10.0, 20.0, 120.0, 35.0, 2.0, 0.05, 40.0, 30.0, 40.0, 1.0, 11.0, 11.5, 10.4, 5200.0, 4.5, -1.2, 0.02, 1200, "J0040", 9.8, 0.02, 9.8, 0.02, 9.7, 0.04, 8.9, 0.2, "0000", 0, "AAAA", 10.0, 20.0, 10.1, 0.02, 9.9, 0.02, 9.8, 0.02, "AAA"],
                    ],
                },
            )
        return _DummyResp(400)

    monkeypatch.setattr(requests, "post", fake_post)
    from oddstar.io.live_survey import fetch_live_halo_bulge_sample

    stars, _ = fetch_live_halo_bulge_sample(cfg, cache, limit=10, region="halo", after_random_index=999)
    assert [s["source_id"] for s in stars] == ["222"]
    assert any("WHERE g.source_id IN (222,333)" in q for q in seen_queries)
    assert not any("SELECT TOP" in q for q in seen_queries)


# ------------------------------------------------------------------------------- key-range slices

import re  # noqa: E402


class _SlicedArchive:
    """A stand-in Gaia table that answers `lo < g.<key> <= hi` pages the way a TAP service does."""

    def __init__(self, keys: list[int], *, key: str = "random_index", fail_wider_than: int | None = None) -> None:
        self.rows = sorted(keys)
        self.key = key
        self.fail_wider_than = fail_wider_than
        self.queries: list[str] = []
        self.ranges: list[tuple[int, int]] = []

    def post(self, url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        q = (data or {}).get("QUERY", "")
        self.queries.append(q)
        m = re.search(rf"g\.{self.key} > (-?\d+) AND g\.{self.key} <= (\d+)", q)
        top = re.search(r"SELECT TOP (\d+) ", q)
        if "ari.uni-heidelberg" in url or "aip.de" in url or not m or not top:
            return _DummyResp(400)
        lo, hi = int(m.group(1)), int(m.group(2))
        if self.fail_wider_than is not None and hi - lo > self.fail_wider_than:
            return _DummyResp(408)
        assert f"ORDER BY g.{self.key} ASC" in q
        self.ranges.append((lo, hi))
        hit = [k for k in self.rows if lo < k <= hi][: int(top.group(1))]
        return _DummyResp(200, json_data={"metadata": [{"name": "source_id"}, {"name": self.key}] if self.key != "source_id"
                                          else [{"name": "source_id"}], "data": [[k, k] if self.key != "source_id" else [k] for k in hit]})


def _sliced(tmp_path: Path, archive: _SlicedArchive, monkeypatch: Any, *, n: int, after: int | None = None, key: str = "random_index") -> list[int]:
    from oddstar.io import live_survey as ls

    cfg, _, _ = load_config()
    monkeypatch.setattr(requests, "post", archive.post)
    rows, mans = ls._fetch_live_sliced(where_gaia="g.ruwe < 1.4", key=key, reg="bulge", n_lim=n, after_key=after, cfg=cfg,
                                       cache=RawCache(tmp_path / "raw", max_cache_gb=1.0))
    assert len(mans) == len(archive.ranges)
    return [int(r[key]) for r in rows]


def test_slices_return_the_first_n_by_key_in_a_sparse_region(tmp_path: Path, monkeypatch: Any) -> None:
    keys = [7 + 40_000_003 * i for i in range(40)]  # about one star per 40 million random indices
    arc = _SlicedArchive(keys)
    assert _sliced(tmp_path, arc, monkeypatch, n=10) == keys[:10]
    widths = [hi - lo for lo, hi in arc.ranges]
    assert widths[0] == 1_000_000 and widths == sorted(widths) and widths[-1] > widths[0]  # short pages widen the range
    assert all(a[1] == b[0] for a, b in zip(arc.ranges, arc.ranges[1:]))  # ranges are contiguous: nothing is skipped


def test_slices_page_through_a_dense_region_and_honour_the_cursor(tmp_path: Path, monkeypatch: Any) -> None:
    keys = list(range(100, 700_000, 137))  # about 5,100 stars inside the first slice: every page comes back full
    arc = _SlicedArchive(keys)
    assert _sliced(tmp_path, arc, monkeypatch, n=2500) == keys[:2500]
    assert len(arc.ranges) == 3 and arc.ranges[1][0] == keys[999] and arc.ranges[2][0] == keys[1999]

    arc2 = _SlicedArchive(keys)
    cursor = keys[1234]
    got = _sliced(tmp_path / "b", arc2, monkeypatch, n=300, after=cursor)
    assert got == keys[1235:1535] and "g.random_index > " + str(cursor) in arc2.queries[0]


def test_slices_narrow_after_a_timeout_and_give_up_after_too_many(tmp_path: Path, monkeypatch: Any) -> None:
    keys = [5 + 9_001 * i for i in range(3000)]
    arc = _SlicedArchive(keys, fail_wider_than=300_000)  # the 1,000,000-wide first slice times out; a quarter of it answers
    assert _sliced(tmp_path, arc, monkeypatch, n=120) == keys[:120]
    assert arc.ranges[0] == (-1, 249_999)

    dead = _SlicedArchive(keys, fail_wider_than=0)
    with pytest.raises(ArchiveAccessError):
        _sliced(tmp_path / "dead", dead, monkeypatch, n=10)
    assert dead.ranges == []


def test_slices_stop_at_the_end_of_the_key_range(tmp_path: Path, monkeypatch: Any) -> None:
    arc = _SlicedArchive([11, 2_000_000_000])
    assert _sliced(tmp_path, arc, monkeypatch, n=50) == [11, 2_000_000_000]
    assert arc.ranges[-1][1] == 2**33 and len(arc.ranges) < 20


def test_an_index_at_the_row_cap_is_not_used(tmp_path: Path, monkeypatch: Any) -> None:
    """A region index that comes back at the archive's row cap is truncated; the sample is read in slices instead."""
    from oddstar.io import live_survey as ls

    cfg, _, _ = load_config()
    monkeypatch.setattr(ls, "GAIA_SYNC_ROW_CAP", 3)
    queries: list[str] = []

    def post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        q = (data or {}).get("QUERY", "")
        queries.append(q)
        if q.startswith("SELECT g.source_id, g.random_index FROM gaiadr3.gaia_source"):
            return _DummyResp(200, json_data={"metadata": [{"name": "source_id"}, {"name": "random_index"}], "data": [[1, 900], [2, 800], [3, 700]]})
        return _DummyResp(200, json_data={"metadata": [], "data": []})

    monkeypatch.setattr(requests, "post", post)
    stars, mans = ls.fetch_live_halo_bulge_sample(cfg, RawCache(tmp_path / "raw", max_cache_gb=1.0), limit=100, region="bulge")
    assert stars == []
    assert not any("g.source_id IN (" in q for q in queries)  # the truncated index never picks the stars
    sliced = [q for q in queries if "AND g.random_index > " in q]
    assert sliced and all("g.l <= 12.0 OR g.l >= 348.0" in q and "ORDER BY g.random_index ASC" in q for q in sliced)
    assert len(mans) == 1 + len(sliced)


def test_local_region_reads_source_id_slices_when_nothing_is_cached(tmp_path: Path, monkeypatch: Any) -> None:
    from oddstar.io import live_survey as ls

    cfg, _, _ = load_config()
    queries: list[str] = []

    def post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        queries.append((data or {}).get("QUERY", ""))
        return _DummyResp(200, json_data={"metadata": [], "data": []})

    monkeypatch.setattr(requests, "post", post)
    stars, _ = ls.fetch_live_halo_bulge_sample(cfg, RawCache(tmp_path / "raw", max_cache_gb=1.0), limit=300, region="local")
    assert stars == [] and queries
    assert all("AND g.source_id > " in q and "ORDER BY g.source_id ASC" in q and "SELECT TOP 1000 " in q for q in queries)
    assert "g.source_id > -1 AND g.source_id <= " + str(2**52 - 1) in queries[0]
    assert str(2**63 - 1) in queries[-1]


def test_a_mirror_that_keeps_failing_is_tried_last_and_a_fallback_mirror_keeps_the_cache_key(tmp_path: Path, monkeypatch: Any) -> None:
    from oddstar.io import live_survey as ls

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    calls: list[str] = []
    down = {"https://hangs.example/sync", "https://gone.example/sync"}

    def fake_post(url: str, data: dict[str, Any] | None = None, timeout: int = 15) -> _DummyResp:
        calls.append(url)
        if url == "https://hangs.example/sync":
            raise requests.exceptions.ReadTimeout("no answer")
        if url in down:
            return _DummyResp(404, json_data={})
        return _DummyResp(200, json_data={"metadata": [{"name": "n"}], "data": [[len(calls)]]})

    monkeypatch.setattr(requests, "post", fake_post)
    urls = ["https://hangs.example/sync", "https://gone.example/sync"]
    kw = {"urls": urls, "fallback_urls": ["https://third.example/sync"], "archive": "simbad", "release": "live", "cache": cache}

    # the listed order while nothing is known about the endpoints, and again after one failure each
    for n in (1, 2):
        calls.clear()
        _tap_query_cached(query=f"SELECT {n}", **kw)
        assert calls == [*urls, "https://third.example/sync"]
    # two failures in a row: the working mirror comes first and the failing ones are not asked at all
    calls.clear()
    _tap_query_cached(query="SELECT 3", **kw)
    assert calls == ["https://third.example/sync"]

    # the key of the answer does not name the fallback mirror: the same query without it is a cache hit
    calls.clear()
    rows, _, hit = _tap_query_cached(query="SELECT 3", urls=urls, archive="simbad", release="live", cache=cache)
    assert hit is True and calls == [] and rows == [{"n": 1}]

    # a mirror that works again is back in its place after one good answer
    down.clear()
    monkeypatch.setattr(requests, "post", lambda url, data=None, timeout=15: (calls.append(url), _DummyResp(200, json_data={"metadata": [{"name": "n"}], "data": [[7]]}))[1])
    ls._note_tap_endpoint("https://hangs.example/sync", ok=True)
    calls.clear()
    _tap_query_cached(query="SELECT 4", **kw)
    assert calls == ["https://hangs.example/sync"]

    # a busy answer (HTTP 429) is not counted as a failure
    ls.reset_tap_endpoint_health()
    monkeypatch.setattr(requests, "post", lambda url, data=None, timeout=15: (calls.append(url), _DummyResp(429 if "hangs" in url else 200, json_data={"metadata": [{"name": "n"}], "data": [[8]]}))[1])
    for n in (5, 6, 7):
        calls.clear()
        _tap_query_cached(query=f"SELECT {n}", sleep=lambda s: None, **kw)
        assert calls[0] == "https://hangs.example/sync"
