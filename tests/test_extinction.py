"""Extinction bookkeeping: ebv_star vs line-of-sight keys, check_extinction, 3D dust query, honest defaults."""

from __future__ import annotations

import base64
import copy
import json
from pathlib import Path
from typing import Any

import numpy as np
import requests

from oddstar.channels import build_channel
from oddstar.config import load_config
from oddstar.io import dust3d
from oddstar.io import live_survey as ls
from oddstar.io.cache import RawCache
from oddstar.io.fixtures_data import get_named_object
from oddstar.io.physics import extinction_two_case_summary, line_of_sight_ebv, star_applied_ebv
from oddstar.vet.underluminous_checks import (
    check_astrometry,
    check_companion,
    check_extinction,
    check_metallicity,
    check_parameter_agreement,
)

from test_photosphere_model import CANDIDATE_ROW, candidate_star


class _Resp:
    def __init__(self, status_code: int, text: str = "") -> None:
        self.status_code = status_code
        self.text = text
        self.content = text.encode("utf-8")


def _b64_array(arr: np.ndarray) -> dict[str, Any]:
    a = np.ascontiguousarray(arr, dtype="<f4")
    return {
        "_type": "np.ndarray",
        "dtype": {"_type": "np.dtype", "descr": "<f4"},
        "shape": list(a.shape),
        "data": base64.b64encode(a.tobytes()).decode("ascii"),
    }


# ------------------------------------------------------------------------------------------- star dict keys


def test_live_row_keeps_gaia_reddening_out_of_line_of_sight_keys() -> None:
    star = ls._row_to_star_dict(dict(CANDIDATE_ROW))
    assert star["ebv_star"] == 0.039 / 1.31 and star["ebv_star_source"] == "gaia_gspphot"
    assert star["ebv_sfd"] is None and star["ebv_planck"] is None
    # M_G uses the applied value
    assert abs(star["abs_g_mag"] - (13.430 - 5 * np.log10(1000 / 4.3149 / 10) - 3.1 * star["ebv_star"] * 0.789)) < 1e-9

    missing = dict(CANDIDATE_ROW)
    missing["ebpminrp_gspphot"] = None
    s2 = ls._row_to_star_dict(missing)
    assert s2["ebv_star"] is None and s2["ebv_star_source"] == "missing"
    assert s2["abs_g_extinction_applied"] is False
    assert star_applied_ebv(s2) == (None, "missing")


def test_fixtures_apply_ebv_star_equal_to_ebv_sfd() -> None:
    for sid in ("733300000000000001", "115725517896200320"):
        s = get_named_object(sid)
        assert s is not None
        assert s["ebv_star"] == s["ebv_sfd"] and s["ebv_star_source"] == "fixture"


def test_line_of_sight_prefers_vetting_dust_and_converts_sfd() -> None:
    star = candidate_star()
    los = line_of_sight_ebv(star, sfd_to_sf11=0.86)
    assert los is not None and abs(los["ebv_los"] - 0.2414 * 0.86) < 1e-12
    star["vetting_archive"]["dust"]["ebv_sf11_los"] = 0.205
    assert line_of_sight_ebv(star)["ebv_los"] == 0.205


def test_underluminous_channel_uses_ebv_star_and_reports_cases(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    ch = build_channel("underluminous", cfg, tmp_path)
    star = candidate_star()
    m = ch.measure(star, {})
    assert abs(m.values["a_g_applied"] - 3.1 * star["ebv_star"] * 0.789) < 1e-12
    assert m.extra["ebv_star_source"] == "gaia_gspphot"
    ext = m.extra["extinction_cases"]
    assert ext is not None and len(ext["cases"]) == 2
    assert "no 3D dust map" in ext["summary"]
    # Below the trigger no summary is produced
    small = copy.deepcopy(star)
    small["vetting_archive"]["dust"]["ebv_sfd_los"] = 0.04
    assert extinction_two_case_summary(small, min_los_excess_ebv=0.05) is None


# ------------------------------------------------------------------------------------------- check_extinction


def test_check_extinction_states(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    star = candidate_star()
    star["live_enriched"] = True

    # Line-of-sight total far above the applied value and no 3D map: uncertain -> not_evaluated, never pass
    r = check_extinction(star, cfg, tmp_path / "a")
    assert r.status == "not_evaluated" and "intrinsic colour and luminosity are uncertain" in r.note

    # Only one estimate (dust query failed) -> not_evaluated
    failed = copy.deepcopy(star)
    failed["vetting_archive"]["query_status"] = {"irsa_dust": "failed: timeout"}
    r = check_extinction(failed, cfg, tmp_path / "b")
    assert r.status == "not_evaluated" and "Fewer than two" in r.note

    # 3D map agreeing with the applied value -> pass (with evidence file)
    agree = copy.deepcopy(star)
    agree["vetting_archive"]["dust3d"] = {"status": "ok", "ebv": 0.04, "ebv_err": 0.02, "map": "bayestar19"}
    r = check_extinction(agree, cfg, tmp_path / "c")
    assert r.status == "pass" and r.evidence is not None and Path(r.evidence).exists()

    # 3D map putting the whole column in front of the star -> applied value disagrees -> fail
    behind = copy.deepcopy(star)
    behind["vetting_archive"]["dust3d"] = {"status": "ok", "ebv": 0.24, "ebv_err": 0.02, "map": "bayestar19"}
    r = check_extinction(behind, cfg, tmp_path / "d")
    assert r.status == "fail" and "applied vs 3D" in r.note

    # Gaia value present but nothing else at all -> not_evaluated (it used to pass by comparing a number to itself)
    lone = ls._row_to_star_dict(dict(CANDIDATE_ROW))
    assert check_extinction(lone, cfg, tmp_path / "e").status == "not_evaluated"

    # Fixture semantics unchanged: heavy-dust impostor fails, clean benchmark passes
    assert check_extinction(get_named_object("733300000000000001"), cfg, tmp_path / "f").status == "fail"
    assert check_extinction(get_named_object("115725517896200320"), cfg, tmp_path / "g").status == "pass"


# ------------------------------------------------------------------------------------------- missing inputs


def test_underluminous_checks_do_not_pass_on_defaults(tmp_path: Path) -> None:
    cfg, _, _ = load_config()
    live = ls._row_to_star_dict(dict(CANDIDATE_ROW))
    live["vetting_archive"]["nss_flag"] = 0
    live["vetting_archive"]["galex_uv_excess_mag"] = 0.0
    live["vetting_archive"]["query_status"] = {}
    # Values present but their queries were not recorded as ok -> not_evaluated
    r = check_astrometry(live, cfg, tmp_path / "a")
    assert r.status == "not_evaluated" and "gaia_meta" in r.note
    r = check_companion(live, cfg, tmp_path / "b")
    assert r.status == "not_evaluated" and "galex" in r.note
    live["vetting_archive"]["query_status"] = {"gaia_meta": "ok", "galex": "ok"}
    assert check_astrometry(live, cfg, tmp_path / "c").status == "pass"
    assert check_companion(live, cfg, tmp_path / "d").status == "pass"
    # Live teff_spec_k is a copy of the Gaia value: not an independent check
    assert check_parameter_agreement(live, cfg, tmp_path / "e").status == "not_evaluated"
    # A failing available input still fails even when another input is missing
    bad = copy.deepcopy(live)
    bad["ruwe"] = 2.5
    bad["vetting_archive"]["query_status"] = {}
    assert check_companion(bad, cfg, tmp_path / "f").status == "fail"

    # Hand-built dicts: absent/None means missing; explicit values count
    hand = {"ruwe": 1.0, "parallax_snr": 80.0, "teff_k": 5000.0, "vetting_archive": {}}
    assert check_astrometry(hand, cfg, tmp_path / "g").status == "not_evaluated"
    assert check_companion(hand, cfg, tmp_path / "h").status == "not_evaluated"
    assert check_metallicity(hand, cfg, tmp_path / "i").status == "not_evaluated"
    assert check_parameter_agreement(hand, cfg, tmp_path / "j").status == "not_evaluated"
    hand["vetting_archive"] = {"nss_flag": 0, "galex_uv_excess_mag": 0.1}
    hand.update({"feh": 0.0, "teff_spec_k": 5050.0})
    for fn in (check_astrometry, check_companion, check_metallicity, check_parameter_agreement):
        assert fn(hand, cfg, tmp_path / f"k_{fn.__name__}").status == "pass"


# ------------------------------------------------------------------------------------------- 3D dust query


def test_bayestar_query_offline(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    calls: list[dict[str, Any]] = []

    def fake_post(url: str, data: bytes | None = None, headers: dict[str, str] | None = None, timeout: float = 0) -> _Resp:
        body = json.loads(data or b"{}")
        calls.append({"url": url, "body": body, "headers": headers})
        # 16/50/84 percentiles at d - sigma, d, d + sigma (Bayestar19 units)
        arr = np.array([[0.02, 0.03, 0.04], [0.025, 0.035, 0.045], [0.03, 0.04, 0.05]])
        return _Resp(200, json.dumps(_b64_array(arr)))

    monkeypatch.setattr(requests, "post", fake_post)
    rec, mans = dust3d.query_bayestar19_at_distance(
        ra_deg=43.54827, dec_deg=5.88271, distance_pc=231.76, parallax_mas=4.3149, parallax_err_mas=0.0162,
        cache=cache, ebv_per_unit=0.884, sleep=lambda s: None,
    )
    assert rec["status"] == "ok" and mans
    assert calls[0]["url"] == dust3d.BAYESTAR19_QUERY_URL
    body = calls[0]["body"]
    assert body["mode"] == "percentile" and body["pct"] == [16.0, 50.0, 84.0]
    assert body["coords"]["frame"] == "galactic" and body["coords"]["distance"]["unit"] == "kpc"
    assert abs(rec["ebv"] - 0.884 * 0.035) < 1e-6
    assert rec["ebv_err"] > 0.884 * 0.01 - 1e-9
    assert rec["ebv_conversion_systematic"] > 0.0

    # Second identical call is served from the content-addressed cache
    def no_network(*a: Any, **k: Any) -> _Resp:
        raise AssertionError("should come from cache")

    monkeypatch.setattr(requests, "post", no_network)
    rec2, _ = dust3d.query_bayestar19_at_distance(
        ra_deg=43.54827, dec_deg=5.88271, distance_pc=231.76, parallax_mas=4.3149, parallax_err_mas=0.0162,
        cache=cache, ebv_per_unit=0.884, sleep=lambda s: None,
    )
    assert rec2["ebv"] == rec["ebv"]


def test_bayestar_query_coverage_and_failures(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    rec, mans = dust3d.query_bayestar19_at_distance(
        ra_deg=10.0, dec_deg=-45.0, distance_pc=100.0, cache=cache, ebv_per_unit=0.884, sleep=lambda s: None
    )
    assert rec["status"] == "no_coverage" and not mans

    def plain_nan(url: str, data: bytes | None = None, headers: Any = None, timeout: float = 0) -> _Resp:
        return _Resp(200, json.dumps([[None, None, None]] * 3).replace("null", "NaN"))

    monkeypatch.setattr(requests, "post", plain_nan)
    rec, _ = dust3d.query_bayestar19_at_distance(
        ra_deg=10.0, dec_deg=20.0, distance_pc=100.0, cache=cache, ebv_per_unit=0.884, sleep=lambda s: None
    )
    assert rec["status"] == "no_coverage"

    def server_error(url: str, data: bytes | None = None, headers: Any = None, timeout: float = 0) -> _Resp:
        return _Resp(503, "busy")

    monkeypatch.setattr(requests, "post", server_error)
    rec, _ = dust3d.query_bayestar19_at_distance(
        ra_deg=11.0, dec_deg=20.0, distance_pc=100.0, cache=cache, ebv_per_unit=0.884, sleep=lambda s: None
    )
    assert rec["status"].startswith("failed:") and rec["ebv"] is None


def test_enrichment_block_stores_dust3d(tmp_path: Path, monkeypatch: Any) -> None:
    cfg, _, _ = load_config()
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)

    def fake_get(url: str, params: Any = None, timeout: float = 0, **k: Any) -> _Resp:
        return _Resp(500, "")

    def fake_post(url: str, data: Any = None, headers: Any = None, timeout: float = 0, **k: Any) -> _Resp:
        if url == dust3d.BAYESTAR19_QUERY_URL:
            return _Resp(200, json.dumps({"_type": "np.ndarray", "dtype": {"_type": "np.dtype", "descr": "<f8"}, "shape": [3, 3], "value": [[0.1, 0.12, 0.14]] * 3}))
        return _Resp(500, "")

    monkeypatch.setattr(requests, "get", fake_get)
    monkeypatch.setattr(requests, "post", fake_post)
    star = ls._row_to_star_dict(dict(CANDIDATE_ROW))
    out, _ = ls.enrich_star_live_for_vetting(star, cfg, cache, enrich_optical=False)
    d3 = out["vetting_archive"]["dust3d"]
    assert d3["status"] == "ok" and abs(d3["ebv"] - cfg.vetting["bayestar19_ebv_per_unit"] * 0.12) < 1e-9
    assert out["live_enriched"] is True
