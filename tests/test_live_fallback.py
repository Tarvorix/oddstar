"""Offline tests for live archive mirror fallback, Gator/VizieR parsers, and MAST PS1 filtering."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import requests

from oddstar.config import load_config
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
    assert enriched["vetting_archive"]["deep_imaging_available"] is True
    assert enriched["vetting_archive"]["deep_imaging_offset_galaxy"] is None


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
    assert len(seen_queries) == 2
    assert "ABS(g.b) > 20.0" in seen_queries[0]
    assert "g.mh_gspphot < -0.8" in seen_queries[0]
    assert "g.random_index > 999" in seen_queries[0]
    assert "ORDER BY g.random_index ASC" in seen_queries[0]
    assert "g.l <= 12.0 OR g.l >= 348.0" in seen_queries[1]
    assert "ABS(g.b) BETWEEN 3.0 AND 10.0" in seen_queries[1]


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
