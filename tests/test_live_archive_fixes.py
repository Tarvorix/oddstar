"""Regressions for the problems the 2026-11-17 live archive run exposed (all offline; every archive is a double).

Each test reproduces the answer an archive really gave that day and checks the pipeline now reads it correctly:
an empty answer from a service that does not hold the table, a "no bodies" HTTP 204, a lost asteroid whose
ephemeris uncertainty spanned degrees, VOTable columns keyed by opaque IDs, a server that cannot be connected to,
and an unfinished MAST job must never be recorded as a real negative (or a real match).
"""

from __future__ import annotations

import copy
import io
import json
from pathlib import Path
from typing import Any

import pytest
import requests

from oddstar.config import load_config
from oddstar.errors import ArchiveAccessError
from oddstar.io import dust3d, high_energy
from oddstar.io import live_survey as ls
from oddstar.io.cache import RawCache
from oddstar.io.skybot import skybot_cone_search
from oddstar.io.vetting_queries import fetch_subarcsec_imaging
from oddstar.vet import ir_checks as irc

from test_ir_vetting_live import DEC, NO_SLEEP, PMDEC, PMRA, RA, FakeArchives, _enrich, _heasarc_votable, _hms, _dms, _Resp, _sia_votable

SKYBOT_HEADER = ("# Flag: 1\n# Ticket: 1\n"
                 "# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec) | dRA(arcsec/h) | dDEC(arcsec/h) | Dg(ua) | Dh(ua)\n")


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


# ------------------------------------------------------------------------------------------------ SkyBoT


def test_skybot_http_204_is_an_empty_field_and_is_cached(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    calls: list[str] = []

    def no_content(url: str, params: Any = None, timeout: float = 0) -> _Resp:
        calls.append(url)
        return _Resp(204)

    monkeypatch.setattr(requests, "get", no_content)
    kw = dict(ra_deg=43.548, dec_deg=5.8827, radius_arcsec=600.0, epoch_jd=2455216.6, observer="C51", cache=cache, sleep=NO_SLEEP)
    bodies, man = skybot_cone_search(**kw)
    assert bodies == [] and man.row_count == 0 and len(calls) == 1
    bodies2, _ = skybot_cone_search(**kw)
    assert bodies2 == [] and len(calls) == 1  # second call served from the cache


def test_http_204_is_a_failure_for_services_that_do_not_use_it(tmp_path: Path, monkeypatch: Any) -> None:
    from oddstar.io.archive_http import cached_get

    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(204))
    with pytest.raises(ArchiveAccessError):
        cached_get(url="https://archive.example/q", params={"a": 1}, archive="x", release="y", cache=RawCache(tmp_path / "raw", max_cache_gb=1.0),
                   timeout_s=5.0, attempts=1, sleep=NO_SLEEP)


def _sso(tmp_path: Path, monkeypatch: Any, cfg: Any, body_line: str) -> dict[str, Any]:
    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(200, SKYBOT_HEADER + body_line))
    sso, mans = ls._search_sso(ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, frame_mjds=[55216.10, 55216.17, 55216.23],
                               vcfg=cfg.vetting, cache=RawCache(tmp_path / "raw", max_cache_gb=1.0), sleep=NO_SLEEP)
    assert mans
    return sso


def test_lost_body_with_huge_ephemeris_uncertainty_is_not_a_match(tmp_path: Path, monkeypatch: Any, cfg: Any) -> None:
    """Live: 2021 VU98 sat 561" from the star with a 225,012" uncertainty; adding that to the match radius failed the check."""
    line = f"- | 2021 VU98 | {_hms(RA)} | {_dms(DEC + 560.88 / 3600.0)} | MB>Inner | 22.9 | 225012.5 | 560.88 | 0.0 | 0.0 | 1.9 | 2.6\n"
    sso = _sso(tmp_path, monkeypatch, cfg, line)
    assert sso["matches"] == []
    assert sso["max_pos_err_arcsec"] == float(cfg.vetting["sso_max_pos_err_arcsec"])
    unc = sso["unconstrained_bodies"]
    assert [u["name"] for u in unc] == ["2021 VU98"] and unc[0]["pos_err_arcsec"] == 225012.5

    star = {"source_id": "1", "live_enriched": True,
            "vetting_archive": {"sso_search": sso, "query_status": {"wise_l1b": "ok", "skybot": "ok"}}}
    res = irc.check_moving_object(star, cfg, tmp_path / "ev")
    assert res.status == "pass" and "tested at their nominal position only" in res.note
    assert res.metrics["unconstrained_bodies"][0]["name"] == "2021 VU98"


def test_well_located_body_within_its_uncertainty_still_matches(tmp_path: Path, monkeypatch: Any, cfg: Any) -> None:
    r_match = float(cfg.vetting["sso_match_radius_arcsec"])
    off = r_match + 10.0  # outside the bare radius, inside radius + a 20" uncertainty (below the cap)
    line = f"4321 | Nearby | {_hms(RA)} | {_dms(DEC + off / 3600.0)} | MB>Middle | 17.9 | 20.0 | {off:.1f} | 0.0 | 0.0 | 1.6 | 2.5\n"
    sso = _sso(tmp_path, monkeypatch, cfg, line)
    assert sso["matches"] and sso["matches"][0]["name"] == "Nearby" and sso["unconstrained_bodies"] == []


def test_skybot_status_text_has_one_failed_prefix(tmp_path: Path) -> None:
    star = _enrich(tmp_path, FakeArchives(http500={"skybot"}))
    status = star["vetting_archive"]["query_status"]["skybot"]
    assert status.startswith("failed: ") and "failed: failed:" not in status


# --------------------------------------------------------------------------------------- IRSA image search


def test_sia_votable_columns_are_read_by_name_not_id() -> None:
    from astropy.io.votable import parse

    payload = _sia_votable(RA, DEC)
    by_id = parse(io.BytesIO(payload)).get_first_table().to_table()
    assert "access_url" not in by_id.colnames  # the double reproduces IRSA's opaque column IDs
    by_name = parse(io.BytesIO(payload)).get_first_table().to_table(use_names_over_ids=True)
    assert {"access_url", "s_ra", "s_dec"} <= set(by_name.colnames)


def test_cutouts_are_measured_from_an_id_keyed_sia_answer(tmp_path: Path) -> None:
    star = _enrich(tmp_path, FakeArchives())
    assert star["vetting_archive"]["query_status"]["cutouts_wise"] == "ok"
    assert set(star["vetting_archive"]["cutouts"]["bands"]) >= {"W3", "W4"}


# ----------------------------------------------------------------------------------------------- Bayestar19


def test_bayestar_unreachable_opens_the_circuit(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    calls: list[Any] = []

    def refuse(url: str, data: Any = None, headers: Any = None, timeout: Any = 0) -> _Resp:
        calls.append(timeout)
        raise requests.exceptions.ConnectTimeout("connect timed out")

    monkeypatch.setattr(requests, "post", refuse)
    kw = dict(dec_deg=20.0, distance_pc=100.0, cache=cache, ebv_per_unit=0.884, sleep=NO_SLEEP)
    rec, _ = dust3d.query_bayestar19_at_distance(ra_deg=10.0, **kw)
    assert rec["status"].startswith("failed:") and "could not connect" in rec["status"] and rec["ebv"] is None
    assert len(calls) == 2 and calls[0][0] == dust3d.CONNECT_TIMEOUT_S  # short connect timeout, one retry only
    rec2, _ = dust3d.query_bayestar19_at_distance(ra_deg=11.0, **kw)
    assert rec2["status"].startswith("failed:") and "unreachable earlier in this run" in rec2["status"]
    assert len(calls) == 2  # no further connection attempts while the circuit is open
    dust3d.reset_unreachable()
    dust3d.query_bayestar19_at_distance(ra_deg=12.0, **kw)
    assert len(calls) == 4


# ------------------------------------------------------------------------------------------- Legacy Survey


def test_legacy_survey_uses_the_positional_index(tmp_path: Path) -> None:
    fake = FakeArchives()
    _enrich(tmp_path, fake)
    q = next(p["QUERY"] for svc, p in fake.calls if svc == "ls")
    assert "WHERE 't' = q3c_radial_query(ra, dec," in q and " BETWEEN " not in q


def test_datalab_allwise_fallback_uses_the_accepted_predicate() -> None:
    """Data Lab's parser rejects a bare `WHERE q3c_radial_query(...)`; both Data Lab queries must compare with 't'."""
    src = Path(ls.__file__).read_text(encoding="utf-8")
    bare = [ln for ln in src.splitlines() if "q3c_radial_query(" in ln and "WHERE" in ln and "'t' = q3c_radial_query(" not in ln]
    assert bare == []


# ------------------------------------------------------------------------------- GALEX and X-ray catalogues


def test_galex_json_parsing_rules() -> None:
    ok = {"status": "COMPLETE", "data": [{"ra": 1.0, "dec": 2.0, "nuv_mag": 21.5, "fuv_mag": -999.0}], "paging": {"rowsTotal": 1}}
    rows, paging = high_energy.parse_galex_json(json.dumps(ok).encode())
    assert rows[0]["nuv_mag"] == 21.5 and paging["rowsTotal"] == 1
    assert high_energy._mag(-999.0) is None and high_energy._mag("x") is None and high_energy._mag(21.5) == 21.5
    for bad in ({"status": "EXECUTING", "data": []}, {"status": "COMPLETE"}, {"status": "ERROR", "data": [], "msg": "boom"}):
        with pytest.raises(ArchiveAccessError):
            high_energy.parse_galex_json(json.dumps(bad).encode())
    with pytest.raises(ArchiveAccessError):
        high_energy.parse_galex_json(b"<html>maintenance</html>")


def test_galex_truncated_cone_is_an_error(tmp_path: Path, monkeypatch: Any) -> None:
    body = {"status": "COMPLETE", "data": [{"ra": 1.0, "dec": 2.0, "nuv_mag": 21.5}], "paging": {"rowsTotal": 9000}}
    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(200, json_data=body))
    with pytest.raises(ArchiveAccessError, match="truncated"):
        high_energy.fetch_galex_sources(ra_deg=1.0, dec_deg=2.0, radius_arcsec=600.0, cache=RawCache(tmp_path / "raw", max_cache_gb=1.0),
                                        sleep=NO_SLEEP)


def test_galex_sources_nearest_first_with_missing_bands(tmp_path: Path, monkeypatch: Any) -> None:
    body = {"status": "COMPLETE", "paging": {"rowsTotal": 2}, "data": [
        {"ra": 10.01, "dec": 20.0, "nuv_mag": 19.0, "nuv_magerr": 0.1, "fuv_mag": 20.0, "fuv_magerr": 0.2, "survey": "MIS", "nuv_exptime": 1500.0},
        {"ra": 10.0, "dec": 20.0003, "nuv_mag": 22.4, "nuv_magerr": 0.4, "fuv_mag": -999.0, "fuv_magerr": -999.0, "survey": "AIS", "nuv_exptime": 200.0},
    ]}
    seen: list[dict[str, Any]] = []

    def get(url: str, params: Any = None, timeout: float = 0) -> _Resp:
        seen.append(json.loads(params["request"]))
        return _Resp(200, json_data=body)

    monkeypatch.setattr(requests, "get", get)
    rows, man = high_energy.fetch_galex_sources(ra_deg=10.0, dec_deg=20.0, radius_arcsec=600.0,
                                                cache=RawCache(tmp_path / "raw", max_cache_gb=1.0), sleep=NO_SLEEP)
    assert seen[0]["service"] == "Mast.Galex.Catalog" and abs(seen[0]["params"]["radius"] - 600.0 / 3600.0) < 1e-6
    assert [r["survey"] for r in rows] == ["AIS", "MIS"] and rows[0]["fuv_mag"] is None and rows[0]["sep_arcsec"] < 2.0
    assert man.row_count is not None


def test_heasarc_votable_rules() -> None:
    cols = ["name", "ra", "dec", "count_rate"]
    rows = high_energy.parse_heasarc_votable(
        _heasarc_votable(cols, [{"name": "2RXS J1", "ra": "10.0", "dec": "20.0", "count_rate": None}]), cols)
    assert rows == [{"name": "2RXS J1", "ra": 10.0, "dec": 20.0, "count_rate": None}]
    assert high_energy.parse_heasarc_votable(_heasarc_votable(cols, []), cols) == []
    with pytest.raises(ArchiveAccessError, match="lacks column"):
        high_energy.parse_heasarc_votable(_heasarc_votable(cols, []), [*cols, "hardness_ratio_1"])
    err = (b'<?xml version="1.0"?><VOTABLE><RESOURCE type="results"><INFO name="QUERY_STATUS" value="ERROR">'
           b"Table not found: rass2rxs_x</INFO></RESOURCE></VOTABLE>")
    with pytest.raises(ArchiveAccessError, match="Table not found"):
        high_energy.parse_heasarc_votable(err, cols)
    with pytest.raises(ArchiveAccessError, match="no QUERY_STATUS"):
        high_energy.parse_heasarc_votable(b"<html>maintenance</html>", cols)


def test_heasarc_cone_builds_a_safe_query(tmp_path: Path, monkeypatch: Any) -> None:
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    cols = ["name", "ra", "dec", "count_rate"]
    seen: list[dict[str, Any]] = []

    def get(url: str, params: Any = None, timeout: float = 0) -> _Resp:
        seen.append({"url": url, **params})
        return _Resp(200, content=_heasarc_votable(cols, [
            {"name": "far", "ra": "10.0", "dec": "20.01", "count_rate": "0.2"},
            {"name": "near", "ra": "10.0", "dec": "20.001", "count_rate": "0.1"},
        ]))

    monkeypatch.setattr(requests, "get", get)
    rows, _ = high_energy.heasarc_cone(table="rass2rxs", columns=cols, ra_deg=10.0, dec_deg=20.0, radius_arcsec=40.0, cache=cache,
                                       release="t", sleep=NO_SLEEP)
    assert [r["name"] for r in rows] == ["near", "far"] and abs(rows[0]["sep_arcsec"] - 3.6) < 0.05
    assert seen[0]["url"] == high_energy.HEASARC_TAP_URL and "FORMAT" not in seen[0]
    assert "FROM rass2rxs WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 10.000000, 20.000000, 0.011111))" in seen[0]["QUERY"]
    with pytest.raises(ArchiveAccessError, match="Unsafe"):
        high_energy.heasarc_cone(table="rass2rxs; DROP", columns=cols, ra_deg=1, dec_deg=2, radius_arcsec=3, cache=cache, release="t")
    with pytest.raises(ArchiveAccessError, match="ra and dec"):
        high_energy.heasarc_cone(table="rass2rxs", columns=["name"], ra_deg=1, dec_deg=2, radius_arcsec=3, cache=cache, release="t")


def test_enrichment_takes_uv_and_xray_from_the_owning_archives(tmp_path: Path) -> None:
    fake = FakeArchives()
    star = _enrich(tmp_path, fake)
    youth = star["vetting_archive"]["youth"]
    assert youth["galex"]["catalogue"] == "MAST Mast.Galex.Catalog"
    assert youth["galex"]["n_sources_in_coverage_cone"] is None  # a match at the star needs no coverage query
    assert sum(1 for svc, _ in fake.calls if svc == "galex") == 1
    assert youth["galex"]["match"]["nuv_mag"] == 22.44 and youth["galex"]["match"]["survey"] == "AIS"
    assert youth["rosat_2rxs"]["catalogue"].startswith("HEASARC rass2rxs") and youth["rosat_2rxs"]["match"] is None
    services = [svc for svc, _ in fake.calls]
    assert "galex" in services and "rosat" in services
    bad = _enrich(tmp_path / "bad", FakeArchives(bad={"rosat"}))
    assert bad["vetting_archive"]["youth"]["rosat_2rxs"]["match"]["rate_cts_s"] == 0.10


def test_field_never_observed_by_galex_is_no_coverage_not_a_non_detection(tmp_path: Path, cfg: Any) -> None:
    star = _enrich(tmp_path, FakeArchives(not_held={"galex"}))
    assert star["vetting_archive"]["query_status"]["galex"] == "no_coverage"
    assert star["vetting_archive"]["youth"]["galex"]["status"] == "no_coverage"
    star["distance_pc"] = 60.0
    res = irc.check_youth(star, cfg, tmp_path / "ev")
    assert res.status == "not_evaluated"
    assert res.metrics["uv"]["verdict"] == "non_constraining" and res.metrics["uv"]["status"] == "no_coverage"


def test_observed_field_without_a_source_at_the_star_is_ok_with_no_match(tmp_path: Path, cfg: Any) -> None:
    """The coverage cone of a deep field overflows one MAST page; any row there still shows the field was observed."""
    fake = FakeArchives(no_target={"galex"})
    star = _enrich(tmp_path, fake)
    galex = star["vetting_archive"]["youth"]["galex"]
    assert star["vetting_archive"]["query_status"]["galex"] == "ok"
    assert galex["status"] == "ok" and galex["match"] is None and galex["n_sources_in_coverage_cone"] == 1
    radii = sorted(round(json.loads(p["request"])["params"]["radius"] * 3600.0) for svc, p in fake.calls if svc == "galex")
    assert radii == [round(cfg.vetting["youth_galex_match_radius_arcsec"]), round(cfg.vetting["youth_galex_coverage_radius_arcsec"])]
    uv = irc._youth_uv(star, galex, cfg)
    assert uv["verdict"] == "non_constraining" and "observed the field" in uv["reason"]


def test_uv_youth_uses_the_reddening_applied_to_the_star(tmp_path: Path, cfg: Any) -> None:
    star = _enrich(tmp_path, FakeArchives())
    galex = star["vetting_archive"]["youth"]["galex"]
    assert star.get("ebv_sfd") is None  # live stars carry no stellar value under the line-of-sight key
    uv = irc._youth_uv(star, galex, cfg)
    assert uv["ebv_used"] == star["ebv_star"] and uv["ebv_source"] == star["ebv_star_source"]
    redder = copy.deepcopy(star)
    redder["ebv_star"] = float(star["ebv_star"]) + 0.3
    assert irc._youth_uv(redder, galex, cfg)["nuv_excess_mag"] != uv["nuv_excess_mag"]
    unknown = copy.deepcopy(star)
    unknown["ebv_star"] = None
    out = irc._youth_uv(unknown, galex, cfg)
    assert out["verdict"] == "non_constraining" and "no stellar reddening estimate" in out["reason"]


# ------------------------------------------------------------------------------------------- VizieR mirrors


def test_mirror_that_does_not_hold_the_table_is_not_a_negative(tmp_path: Path) -> None:
    """Live: the CfA mirror answered 200 with metadata and zero rows for tables it does not serve."""
    star = _enrich(tmp_path, FakeArchives(not_held={"radio"}))
    status = star["vetting_archive"]["query_status"]["radio_nvss"]
    assert status.startswith("failed: ") and "canary" in status
    ls._vizier_mirror_serves.clear()  # the verdict on a mirror is remembered for the rest of a run
    assert _enrich(tmp_path / "ok", FakeArchives())["vetting_archive"]["query_status"]["radio_nvss"] == "ok"


def test_vizier_table_without_a_canary_cannot_return_an_empty_answer(tmp_path: Path, monkeypatch: Any) -> None:
    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(200, "#\n#RESOURCE=x\nRAJ2000\tDEJ2000\ndeg\tdeg\n---\t---\n"))
    with pytest.raises(ArchiveAccessError, match="no canary position"):
        ls._vizier_tsv_cone_cached(source="J/unknown/table", ra=1.0, dec=2.0, radius_arcsec=5.0, archive="x", release="y",
                                   cache=RawCache(tmp_path / "raw", max_cache_gb=1.0))


def test_leaked_cgi_headers_are_stripped() -> None:
    body = "#\n#RESOURCE=yCat\nRAJ2000\tDEJ2000\ndeg\tdeg\n---\t---\n10.0\t20.0\n"
    leaked = "Content-Type: text/plain\nX-Powered-By: vizier\n\n" + body
    assert ls._strip_leaked_cgi_headers(leaked) == body
    assert ls._strip_leaked_cgi_headers(body) == body
    html = "<html>\nerror\n#not a table\n"
    assert ls._strip_leaked_cgi_headers(html) == html
    assert ls._parse_vizier_tsv(ls._strip_leaked_cgi_headers(leaked)) == [{"RAJ2000": "10.0", "DEJ2000": "20.0"}]


def test_vizier_header_without_rows_is_empty() -> None:
    """A table printed with names, units and dashes but no data must not turn its units line into a row."""
    assert ls._parse_vizier_tsv("#\n#RESOURCE=x\nRAJ2000\tDEJ2000\ndeg\tdeg\n---\t---\n") == []
    assert ls._parse_vizier_tsv("#\n# no table at all\n") == []
    one = ls._parse_vizier_tsv("RAJ2000\tDEJ2000\ndeg\tdeg\n---\t---\n1.5\t\n")
    assert one == [{"RAJ2000": "1.5", "DEJ2000": None}]


# --------------------------------------------------------------------------------------------------- MAST


def test_cached_payload_that_fails_todays_validation_is_fetched_again(tmp_path: Path, monkeypatch: Any) -> None:
    """An unfinished MAST job cached by the older, laxer rule must not be served forever."""
    from oddstar.io.archive_http import cached_get, request_signature

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    url, params = "https://archive.example/q", {"a": 1}
    desc, sig = request_signature(url, params)
    cache.put(archive="x", release="y", request_signature=sig, query_text=desc, payload=b'{"status": "EXECUTING"}', row_count=0)
    calls: list[str] = []

    def get(u: str, params: Any = None, timeout: float = 0) -> _Resp:
        calls.append(u)
        return _Resp(200, '{"status": "COMPLETE"}')

    def validate(b: bytes) -> None:
        if b"COMPLETE" not in b:
            raise ArchiveAccessError("not complete", retryable=True)

    monkeypatch.setattr(requests, "get", get)
    kw = dict(url=url, params=params, archive="x", release="y", cache=cache, timeout_s=5.0, sleep=NO_SLEEP, validate=validate)
    payload, _, hit = cached_get(**kw)
    assert payload == b'{"status": "COMPLETE"}' and hit is False and len(calls) == 1
    payload, _, hit = cached_get(**kw)
    assert payload == b'{"status": "COMPLETE"}' and hit is True and len(calls) == 1
    assert RawCache(tmp_path / "raw", max_cache_gb=1.0).lookup("x", "y", sig).sha256 == _.sha256  # replacement survives a reload


def test_vizier_answers_cached_before_the_mirror_check_are_not_reused(tmp_path: Path) -> None:
    import hashlib

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    desc = f"VIZIER_TSV VIII/65/nvss cone({1.0:.6f}, {2.0:+.6f}, {10.0:.2f} arcsec)"
    old_sig = hashlib.sha256(f"{'|'.join(ls.VIZIER_TSV_URLS)}|{desc}".encode("utf-8")).hexdigest()[:24]
    cache.put(archive="radio_continuum", release="nvss", request_signature=old_sig, query_text=desc, payload=b"[]", row_count=0)
    fake = FakeArchives()
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(requests, "get", fake.get)
        rows, _, hit = ls._vizier_tsv_cone_cached(source="VIII/65/nvss", ra=1.0, dec=2.0, radius_arcsec=10.0, archive="radio_continuum",
                                                  release="nvss", cache=cache)
    assert rows == [] and hit is False
    assert [svc for svc, _ in fake.calls] == ["radio", "radio"]  # the cone, then the canary that shows the mirror holds NVSS


def test_unfinished_mast_job_is_not_an_empty_result(tmp_path: Path, monkeypatch: Any) -> None:
    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(200, json_data={"status": "EXECUTING", "data": []}))
    with pytest.raises(ArchiveAccessError):
        fetch_subarcsec_imaging(ra_deg=RA, dec_deg=DEC, radius_arcsec=10.0, cache=RawCache(tmp_path / "raw", max_cache_gb=1.0), sleep=NO_SLEEP)


# ------------------------------------------------------------------------------- TAP JSON layouts

def _irsa_votable_json(rows: list[list[str]], status: str = "OK") -> dict[str, Any]:
    """The layout IRSA TAP returns for FORMAT=json (copied from a live answer on 2026-11-17): every value a string."""
    fields = [("designation", "char", "*"), ("ra", "double", None), ("w3mpro", "float", None), ("nb", "long", None), ("cc_flags", "char", "*")]
    return {"VOTABLE": {"<xmlattr>": {"version": "1.3"}, "RESOURCE_ARRAY": [{
        "<xmlattr>": {"type": "results"},
        "INFO_ARRAY": [{"<xmlattr>": {"name": "QUERY_STATUS", "value": status}}],
        "TABLE": {"FIELD_ARRAY": [{"<xmlattr>": {"name": n, "datatype": d, **({"arraysize": a} if a else {})}, "DESCRIPTION": n}
                                  for n, d, a in fields],
                  "DATA": {"TABLEDATA": rows}}}]}}


def test_irsa_votable_json_is_read_with_types_and_nulls() -> None:
    rows = ls._rows_from_tap_json(_irsa_votable_json([["J025411.56+055258.2", "43.5481993", "8.701", "1", "0000"],
                                                      ["J025410.03+055309.9", "43.5418064", "", "2", "Dd00"]]))
    assert rows == [{"designation": "J025411.56+055258.2", "ra": 43.5481993, "w3mpro": 8.701, "nb": 1, "cc_flags": "0000"},
                    {"designation": "J025410.03+055309.9", "ra": 43.5418064, "w3mpro": None, "nb": 2, "cc_flags": "Dd00"}]
    assert ls._rows_from_tap_json(_irsa_votable_json([])) == []  # a real empty answer stays an empty answer
    assert ls._rows_from_tap_json({"metadata": [{"name": "a"}], "data": [[1], [2]]}) == [{"a": 1}, {"a": 2}]


@pytest.mark.parametrize("payload", [
    {"VOTABLE": {"RESOURCE_ARRAY": [{"INFO_ARRAY": []}]}},           # no result table
    {"something": "else"},                                           # neither layout
    {"metadata": [{"name": "a"}]},                                   # half of the Gaia layout
    _irsa_votable_json([["only", "three", "values"]]),              # row does not match its fields
    _irsa_votable_json([], status="ERROR"),
])
def test_an_unreadable_tap_answer_is_an_error_not_an_empty_answer(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, payload: dict[str, Any]) -> None:
    with pytest.raises(ValueError):
        ls._rows_from_tap_json(payload)

    class R:
        status_code = 200
        headers: dict[str, str] = {}

        def json(self) -> Any:
            return payload

    monkeypatch.setattr(requests, "post", lambda url, data=None, timeout=15: R())
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    with pytest.raises(ArchiveAccessError):
        ls._tap_query_cached(url=ls.IRSA_TAP_URL, query="SELECT 1", archive="allwise", release="t", cache=cache, sleep=lambda s: None)
    assert cache.lookup("allwise", "t", ls._tap_signature([ls.IRSA_TAP_URL], "SELECT 1")) is None  # nothing stored


def test_irsa_answers_stored_before_the_layout_was_read_are_not_served(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    import hashlib

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    query = "SELECT designation FROM allwise_p3as_psd"
    old_sig = hashlib.sha256(f"{ls.IRSA_TAP_URL}|{query}".encode("utf-8")).hexdigest()[:24]
    cache.put(archive="allwise", release="t", request_signature=old_sig, query_text=query, payload=b"[]", row_count=0, offline_fixture=False)

    class R:
        status_code = 200
        headers: dict[str, str] = {}

        def json(self) -> Any:
            return _irsa_votable_json([["J1", "1.5", "7.0", "1", "0000"]])

    monkeypatch.setattr(requests, "post", lambda url, data=None, timeout=15: R())
    rows, _, hit = ls._tap_query_cached(url=ls.IRSA_TAP_URL, query=query, archive="allwise", release="t", cache=cache)
    assert hit is False and rows[0]["designation"] == "J1"
    # Gaia keys are unchanged, so the samples and indexes already stored still answer
    gaia_sig = hashlib.sha256(f"{'|'.join(ls.GAIA_TAP_URLS)}|{query}".encode("utf-8")).hexdigest()[:24]
    assert ls._tap_signature(ls.GAIA_TAP_URLS, query) == gaia_sig


def test_asassn_client_that_never_answers_is_a_recorded_failure(tmp_path: Path) -> None:
    import threading
    import time

    from oddstar.io import optical_timeseries as ot

    release = threading.Event()

    class Hung:
        def cone_search(self, **kw: Any) -> Any:
            release.wait(30.0)  # stands in for a read on a connection the server already closed
            raise RuntimeError("released")

    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    t0 = time.monotonic()
    summary, mans = ot.fetch_asassn_lightcurve(ra=RA, dec=DEC, pmra=PMRA, pmdec=PMDEC, cache=cache, client_factory=Hung, timeout_s=0.3)
    release.set()
    assert time.monotonic() - t0 < 5.0
    assert summary is not None and summary["status"] == "query_failed" and "no answer within" in summary["query_log"][0]
    assert mans == []  # a timeout is not stored as an answer
    assert ot._call_with_deadline(lambda: 7, 1.0, "x") == 7
    with pytest.raises(KeyError):
        ot._call_with_deadline(lambda: {}["missing"], 1.0, "x")
    assert ot.ASASSN_CLIENT_TIMEOUT_S >= 60.0
