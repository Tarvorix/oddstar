"""Plate sources around live sample positions for the vanishing search (offline, fake HTTP)."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pytest

from oddstar.channels import build_channel
from oddstar.cli import CHANNELS_WITHOUT_LIVE_SOURCE, dispatch
from oddstar.config import load_config, repo_root
from oddstar.errors import ArchiveAccessError
from oddstar.io import live_survey as ls
from oddstar.io import optical_timeseries as ot
from oddstar.io import plate_survey as ps
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec
from oddstar.vet.vanishing_checks import run_vanishing_vetting

COLS = ps.USNO_COLUMNS.split(",")


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


def _offset(ra: float, dec: float, east_arcsec: float, north_arcsec: float) -> tuple[float, float]:
    return ra + east_arcsec / 3600.0 / math.cos(math.radians(dec)), dec + north_arcsec / 3600.0


def _usno_row(uid: str, ra: float, dec: float, *, b1: float | None = 16.0, r1: float | None = 15.0, pm: tuple[int, int] = (0, 0),
              epoch: float = 1975.0, flags: str = "...", ndet: int = 5) -> dict[str, Any]:
    return {"USNO-B1.0": uid, "RAJ2000": f"{ra:.6f}", "DEJ2000": f"{dec:+.6f}", "Epoch": f"{epoch:.1f}", "pmRA": str(pm[0]), "pmDE": str(pm[1]),
            "Ndet": str(ndet), "Flags": flags, "B1mag": "" if b1 is None else f"{b1:.2f}", "R1mag": "" if r1 is None else f"{r1:.2f}",
            "B2mag": "16.10", "R2mag": "15.10", "Imag": "14.50"}


def _tsv(rows: list[dict[str, Any]]) -> str:
    head = "#\n#   VizieR Astronomical Server\n#Name: I/284/out\n"
    lines = ["\t".join(COLS), "\t".join(" " for _ in COLS), "\t".join("-" * 5 for _ in COLS)]
    lines += ["\t".join(str(r.get(c, "")) for c in COLS) for r in rows]
    return head + "\n".join(lines) + "\n"


class _Resp:
    def __init__(self, status_code: int = 200, text: str = "", content: bytes = b"") -> None:
        self.status_code, self.text, self.content = status_code, text, content


class _Sky:
    """Fake VizieR (USNO-B1.0), Gaia TAP and Pan-STARRS: rows per field centre. Counts every request."""

    def __init__(self) -> None:
        self.usno: dict[tuple[float, float], list[dict[str, Any]] | Exception] = {}
        self.gaia: dict[tuple[float, float], list[dict[str, Any]]] = {}
        self.ps1: dict[tuple[float, float], list[list[Any]] | Exception] = {}
        self.requests: list[str] = []

    @staticmethod
    def _key(ra: float, dec: float) -> tuple[float, float]:
        return round(float(ra), 4), round(float(dec), 4)

    def install(self, monkeypatch: pytest.MonkeyPatch) -> None:
        def vizier_get(url: str, params: dict[str, Any], timeout: float) -> _Resp:
            ra, dec = (float(x) for x in params["-c"].split())
            self.requests.append(f"usno {ra:.4f} {dec:.4f}")
            assert params["-source"] == ps.USNO_TABLE and params["-out"] == ps.USNO_COLUMNS
            if self._key(ra, dec) == self._key(*ps.USNO_CANARY_CONE[:2]):
                return _Resp(text=_tsv([_usno_row("0922-0238660", ra, dec)]))
            body = self.usno.get(self._key(ra, dec), [])
            if isinstance(body, Exception):
                raise body
            return _Resp(text=_tsv(body))

        def gaia_tap(*, urls: list[str], query: str, archive: str, release: str, cache: RawCache, timeout_s: int, sleep: Any = None) -> Any:
            ra, dec, rad = (float(x) for x in query.split("CIRCLE('ICRS', ")[1].rstrip(")").split(", "))
            self.requests.append(f"gaia {ra:.4f} {dec:.4f}")
            assert abs(rad * 3600.0 - 150.0) < 0.01  # field radius plus the largest halo radius

            class _Man:
                manifest_id = "gaia-manifest"

            return self.gaia.get(self._key(ra, dec), []), _Man(), False

        def ps1_get(url: str, params: dict[str, Any] | None, **kw: Any) -> _Resp:
            assert url == ls.MAST_PS1_DR2_URL and params is not None
            assert params["columns"] == ps.PS1_COLUMNS and params["nDetections.gte"] == "2"
            self.requests.append(f"ps1 {float(params['ra']):.4f} {float(params['dec']):.4f}")
            body = self.ps1.get(self._key(params["ra"], params["dec"]), [])
            if isinstance(body, Exception):
                raise body
            info = [{"name": c} for c in ("objID", "raMean", "decMean", "nDetections", "rMeanPSFMag")]
            return _Resp(content=json.dumps({"info": info, "data": body}).encode("utf-8"))

        monkeypatch.setattr(ps.requests, "get", vizier_get)
        monkeypatch.setattr(ls, "_tap_query_cached", gaia_tap)
        monkeypatch.setattr(ot, "_get_with_retries", ps1_get)
        ps._usno_mirror_serves.clear()


F1, F2, F3, F4, F5 = (150.0, 20.0), (60.0, -45.0), (200.0, 10.0), (220.0, 30.0), (10.0, 40.0)


def _fields() -> tuple[list[dict[str, Any]], _Sky]:
    """Five fields: F1 full answers, F2 south of Pan-STARRS, F3 Pan-STARRS fails, F4 USNO fails, F5 no plate source at all."""
    sky = _Sky()
    stars = [{"source_id": str(i + 1), "ra_deg": f[0], "dec_deg": f[1], "mags": {"G": 12.0}, "distance_pc": 300.0, "teff_k": 5500.0}
             for i, f in enumerate((F1, F2, F3, F4, F5))]
    a, b, c, d, e, f, g = (_offset(*F1, x, y) for x, y in ((0, 0), (20, 10), (-25, 15), (8, -8), (30, -30), (-40, 5), (10, 40)))
    sky.usno[sky._key(*F1)] = [
        _usno_row("A", *a, pm=(4, -6)),  # the sample star itself: still there
        _usno_row("B", *b, r1=17.0, b1=18.0, ndet=2),  # on both first-epoch plates, nothing there now
        _usno_row("C", *c, pm=(500, 0)),  # high proper motion; catalogue position at 2000
        _usno_row("D", *d, flags=".s."),  # on a diffraction spike of the sample star
        _usno_row("E", *e, r1=None),  # one first-epoch plate only
        _usno_row("F", *f, r1=20.6, b1=21.5),  # still there or not, too close to the modern depth to tell
        _usno_row("G", *g, epoch=1975.0),  # no catalogue proper motion: position is that of 1975
    ]
    c16 = _offset(*c, 0.5 * 16.0, 0.0)  # 500 mas/yr for 16 years
    g16 = _offset(*g, 0.1 * 41.0, 0.0)  # 100 mas/yr since 1975: 4.1 arcsec from the plate position in 2016
    sky.gaia[sky._key(*F1)] = [
        {"source_id": 11, "ra": a[0], "dec": a[1], "pmra": 4.0, "pmdec": -6.0, "phot_g_mean_mag": 12.0},
        {"source_id": 13, "ra": c16[0], "dec": c16[1], "pmra": 500.0, "pmdec": 0.0, "phot_g_mean_mag": 14.8},
        {"source_id": 17, "ra": g16[0], "dec": g16[1], "pmra": 100.0, "pmdec": 0.0, "phot_g_mean_mag": 15.2},
        {"source_id": 19, "ra": _offset(*b, 7.0, 0.0)[0], "dec": b[1], "pmra": None, "pmdec": None, "phot_g_mean_mag": 20.4},  # 7" from B
    ]
    filler = [[f"p{i}", *_offset(*F1, -50 + 8 * i, -45), 12, 19.0 + 0.25 * i] for i in range(12)]  # r = 19.0 .. 21.75
    sky.ps1[sky._key(*F1)] = filler + [["pa", *a, 60, 12.3], ["pf", *f, 9, 20.7], ["junk", *_offset(*F1, 55, 0), 2, -999.0]]
    sky.usno[sky._key(*F2)] = [_usno_row("S1", *F2), _usno_row("S2", *_offset(*F2, 15, 0), r1=16.0)]
    sky.usno[sky._key(*F3)] = [_usno_row("T1", *F3)]
    sky.ps1[sky._key(*F3)] = ArchiveAccessError("GET failed after 3 attempt(s): HTTP 503", retryable=True)
    sky.usno[sky._key(*F4)] = ArchiveAccessError("HTTP 500", retryable=True)
    return stars, sky


def test_usno_rows_are_typed_and_spurious_classes_marked() -> None:
    rows = ps.parse_usno_rows(ls._parse_vizier_tsv(_tsv([_usno_row("X", 10.0, 5.0, pm=(4, -6)), _usno_row("Y", 10.1, 5.0, r1=None, flags="Ms."),
                                                       _usno_row("Z", 10.2, 5.0, epoch=1968.4), {"USNO-B1.0": "", "RAJ2000": "1", "DEJ2000": "1"}])))
    assert [r["usno_id"] for r in rows] == ["X", "Y", "Z"]
    x, y, z = rows
    assert x["n_first_epoch_plates"] == 2 and x["r1_mag"] == 15.0 and x["pmra_mas_yr"] == 4.0 and not x["on_spike"] and x["position_epoch"] == 2000.0
    assert y["n_first_epoch_plates"] == 1 and y["r1_mag"] is None and y["on_spike"]
    assert z["position_epoch"] == 1968.4  # no catalogue proper motion: the position was never moved to 2000


def test_local_depth_needs_enough_objects_and_ignores_missing_magnitudes() -> None:
    rows = [{"rMeanPSFMag": 19.0 + 0.25 * i} for i in range(12)] + [{"rMeanPSFMag": -999.0}, {"rMeanPSFMag": None}]
    depth, n = ps.local_depth_mag(rows, 90.0, 10)
    assert n == 12 and depth == pytest.approx(21.475)
    assert ps.local_depth_mag(rows[:9], 90.0, 10) == (None, 9)


def test_plate_units_statuses_matches_and_resume(cfg: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    stars, sky = _fields()
    sky.install(monkeypatch)
    cache = RawCache(tmp_path / "raw")
    seen: list[tuple[str, str]] = []
    units, summary = ps.fetch_plate_units(stars, cfg, cache, sleep=lambda s: None, progress=lambda d, t, sid, st: seen.append((sid, st.split(":")[0])))
    by = {u["source_id"].split()[-1]: u for u in units}
    # D (spike) and E (one first-epoch plate) are not units; the field whose plate catalogue failed has none
    assert sorted(by) == ["A", "B", "C", "F", "G", "S1", "S2", "T1"]
    assert {k: ps.plate_search_status(u) for k, u in by.items()} == {
        "A": "searched", "B": "searched", "C": "searched", "G": "searched", "F": "too_faint", "S1": "no_coverage", "S2": "no_coverage", "T1": "failed"}
    assert summary["n_fields_ok"] == 2 and summary["n_fields_no_coverage"] == 1 and summary["n_fields_failed_modern"] == 1
    assert summary["n_fields_failed_usno"] == 1 and summary["failed_field_ids"] == ["3", "4"] and summary["n_units_searched"] == 4
    errs = summary["failed_field_errors"]
    assert sorted(errs) == ["3", "4"] and sum(e.startswith("usno_b1: ") for e in errs.values()) == 1, "each failed field keeps its error text"
    assert summary["n_usno_sources"] == 10 and summary["n_on_spike"] == 1 and summary["n_one_first_epoch_plate"] == 1
    assert seen == [("1", "ok"), ("2", "no_coverage"), ("3", "failed"), ("4", "failed"), ("5", "ok")]

    vd = {k: by[k]["vanishing_data"] for k in ("A", "B", "C", "G", "F")}
    assert vd["A"]["present_in_gaia"] and vd["A"]["present_in_panstarrs"] and vd["A"]["panstarrs_depth_mag"] == pytest.approx(21.425)  # 14 objects with an r magnitude
    assert not vd["B"]["present_in_gaia"] and not vd["B"]["present_in_panstarrs"] and vd["B"]["early_plates_count"] == 2
    assert vd["B"]["nearest_gaia"]["sep_arcsec"] == pytest.approx(7.0, abs=0.05)  # kept for review, too far to be the same object
    # C: Gaia sees it 8 arcsec away in 2016; moved back to 2000 it is the plate source
    assert vd["C"]["present_in_gaia"] and vd["C"]["usno_epoch"] == 2000.0 and vd["C"]["nearest_gaia"]["sep_arcsec"] < 0.1
    assert by["C"]["pmra_mas_yr"] == 500.0
    # G: no catalogue proper motion, so the comparison is made at the mean plate epoch, where Gaia's star was
    assert vd["G"]["present_in_gaia"] and vd["G"]["usno_epoch"] == 1975.0 and vd["G"]["nearest_gaia"]["sep_arcsec"] < 0.1
    assert vd["F"]["present_in_panstarrs"] and by["F"]["plate_fetch"]["status"] == "too_faint"
    # A source that was not compared carries no comparison values, so nothing can be read from it
    assert by["S1"]["vanishing_data"] == {} and by["T1"]["vanishing_data"] == {} and by["T1"]["plate_fetch"]["status"].startswith("failed: modern")
    assert by["B"]["vetting_archive"]["query_status"] == {"usno_b1": "ok", "modern_counterpart": "ok", "panstarrs_depth": "ok"}
    assert by["S1"]["vetting_archive"]["query_status"] == {"usno_b1": "ok"}

    cov = ps.search_coverage_summary(units, cfg)
    assert cov["n_input"] == 8 and cov["n_searched"] == 4 and cov["n_too_faint"] == 1 and cov["n_no_coverage"] == 2 and cov["n_failed"] == 1
    assert cov["not_searched_ids"] == {"too_faint": ["USNO-B1.0 F"], "no_coverage": ["USNO-B1.0 S1", "USNO-B1.0 S2"], "failed": ["USNO-B1.0 T1"]}

    # The channel flags only the plate source with no modern counterpart
    ch = build_channel("vanishing", cfg, repo_root())
    flagged = [u["source_id"] for u in units if ps.plate_search_status(u) == "searched" and ch.score(ch.measure(u, {})).flagged]
    assert flagged == ["USNO-B1.0 B"]
    vet = {r.check: r.status for r in run_vanishing_vetting(by["B"], cfg, tmp_path / "vet")}
    # The catalogue row alone cannot show that a plate source is real, or that it did not move: without the plate
    # images, no verdict on either
    assert vet == {"proper_motion": "not_evaluated", "plate_defect": "not_evaluated", "depth": "pass", "known_variable": "not_evaluated",
                   "sso_plate": "not_evaluated", "blend_split": "not_evaluated"}
    assert {r.status for r in run_vanishing_vetting(by["S1"], cfg, tmp_path / "vet2")} == {"not_evaluated"}

    # Second run: answers come from the cache; the two failed queries are asked again (failures are never stored)
    n_before = len(sky.requests)
    sky.ps1[sky._key(*F3)] = [[f"q{i}", *_offset(*F3, 5 * i, 20), 5, 19.0 + 0.2 * i] for i in range(12)] + [["t1", *F3, 40, 14.9]]
    units2, summary2 = ps.fetch_plate_units(stars, cfg, cache, sleep=lambda s: None)
    again = sky.requests[n_before:]
    assert [r.split()[0] for r in again].count("ps1") == 1 and not any(r.startswith("usno 150.") for r in again)
    assert summary2["n_fields_failed_modern"] == 0 and summary2["n_fields_failed_usno"] == 1 and summary2["n_units_searched"] == 5

    units3, summary3 = ps.fetch_plate_units(stars, cfg, cache, max_fields=1, sleep=lambda s: None)
    assert summary3["n_fields_requested"] == 1 and summary3["n_fields_not_fetched"] == 4 and len(units3) == 5


def test_moving_stars_are_matched_along_their_path_and_bright_star_halos_are_not_searched(cfg: Any, tmp_path: Path) -> None:
    st = ps._settings(cfg)
    assert ps.gaia_margin_arcsec(st) == 90.0
    assert [ps.halo_radius_arcsec(g, st["halo"]) for g in (9.1, 10.0, 10.9, 11.99, 12.0, None)] == [90.0, 50.0, 50.0, 20.0, None, None]
    ch = build_channel("vanishing", cfg, repo_root())
    ps1 = [{"objID": f"p{i}", "raMean": 150.0 + 0.05, "decMean": 20.0 + 0.05 + i * 1e-3, "nDetections": 12, "rMeanPSFMag": 19.0 + 0.25 * i} for i in range(12)]

    def units(star_g: float, gaia: list[dict[str, Any]], rows: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
        star = {"source_id": "1", "ra_deg": 150.0, "dec_deg": 20.0, "mags": {"G": star_g}}
        out, _ = ps.build_field_units(star, ps.parse_usno_rows(rows), gaia, ps1, st)
        return {u["source_id"].split()[-1]: u for u in out}

    # A 10th-magnitude star moving half an arcsecond a year (measured case: HIP 59490). The plates saw it where it
    # was in 1960; the catalogue gave it a wrong proper motion, so its epoch-2000 position is 20 arcsec off the star.
    pm = (-218.0, -439.0)
    fast = {"source_id": "900", "ra_deg": 150.0, "dec_deg": 20.0, "pmra_mas_yr": pm[0], "pmdec_mas_yr": pm[1], "g_mag": 10.02}
    in_1960 = _offset(150.0, 20.0, -pm[0] / 1000.0 * 56.0, -pm[1] / 1000.0 * 56.0)
    wrong_pm = (-206, 114)
    at_2000 = _offset(*in_1960, wrong_pm[0] / 1000.0 * 33.4, wrong_pm[1] / 1000.0 * 33.4)
    by = units(10.02, [fast], [
        _usno_row("H", *at_2000, r1=9.77, b1=10.7, pm=wrong_pm, epoch=1966.6, ndet=3),
        _usno_row("I", *_offset(*in_1960, 30.0, 0.0), r1=17.0, b1=18.0, epoch=1954.0, ndet=2),   # 28 arcsec from the star's path
        _usno_row("J", *_offset(150.0, 20.0, -55.0, -20.0), r1=18.0, b1=19.0, epoch=1954.0, ndet=2),  # 64 arcsec from it
    ])
    h = by["H"]["vanishing_data"]
    assert angular_separation_arcsec(150.0, 20.0, *at_2000) > 25.0 and h["nearest_gaia"]["sep_arcsec"] > 10.0
    assert h["present_in_gaia"] and h["gaia_match"] == "track" and h["nearest_gaia_track"]["sep_arcsec"] < 0.1
    assert h["usno_epoch"] == pytest.approx(1960.0, abs=0.1) and by["H"]["pmdec_mas_yr"] == -439.0
    assert ps.plate_search_status(by["H"]) == "searched" and not ch.score(ch.measure(by["H"], {})).flagged, "the star itself is not in its own halo"
    assert ps.plate_search_status(by["I"]) == "near_bright_star"
    assert by["I"]["plate_fetch"]["bright_star"] == {"source_id": "900", "g_mag": 10.02, "sep_arcsec": pytest.approx(28.0, abs=0.5), "halo_radius_arcsec": 50.0}
    assert ps.plate_search_status(by["J"]) == "searched" and ch.score(ch.measure(by["J"], {})).flagged
    # With the plate images read: the source is on both first-epoch plates and no Pan-STARRS object could be it, moved
    seen = {"status": "ok", "registered": True, "plate_epoch": 1954.0, "peak": {"snr": 25.0, "offset_arcsec": 0.4}}
    j = json.loads(json.dumps(by["J"]))
    j["vanishing_data"]["plate_images"] = {"rows": [dict(seen, survey="poss1_red"), dict(seen, survey="poss1_blue")], "errors": {}}
    j["vanishing_data"]["plate_neighbours"] = {"radius_arcsec": 10.0, "inner_radius_arcsec": 3.0, "n_objects_in_radius": 0, "neighbours": []}
    j["vetting_archive"]["query_status"].update({"plate_images": "ok", "plate_neighbours": "ok"})
    vet = {r.check: r for r in run_vanishing_vetting(j, cfg, tmp_path / "vet")}
    assert vet["plate_defect"].status == "pass"
    assert vet["proper_motion"].status == "pass" and vet["proper_motion"].metrics["nearest_track_sep_arcsec"] == pytest.approx(64.0, abs=1.5)
    assert "paths of 1 Gaia star(s) over 1949.0-2002.0" in vet["proper_motion"].note
    assert run_vanishing_vetting(by["J"], cfg, tmp_path / "vet0")[0].status == "not_evaluated", "not without the plate images"
    # Without the record of the paths (a run made before config 1.6.0) the check cannot pass
    old = json.loads(json.dumps(by["J"]))
    del old["vanishing_data"]["nearest_gaia_track"]
    assert {r.check: r.status for r in run_vanishing_vetting(old, cfg, tmp_path / "vet2")}["proper_motion"] == "not_evaluated"
    # A plate source on the path of a faint moving star that the catalogue epoch misses fails the check when not matched
    moved = json.loads(json.dumps(by["J"]))
    moved["vanishing_data"]["nearest_gaia_track"]["sep_arcsec"] = 1.2
    assert {r.check: r.status for r in run_vanishing_vetting(moved, cfg, tmp_path / "vet3")}["proper_motion"] == "fail"

    # An 11.5-magnitude star that does not move: 20 arcsec halo
    still = {"source_id": "901", "ra_deg": 150.0, "dec_deg": 20.0, "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0, "g_mag": 11.5}
    by = units(11.5, [still], [
        _usno_row("K", *_offset(150.0, 20.0, 0.0, 15.0), r1=18.0, b1=19.0, epoch=1954.0, ndet=2),
        _usno_row("L", *_offset(150.0, 20.0, 0.0, 26.0), r1=18.0, b1=19.0, epoch=1954.0, ndet=2),
        _usno_row("M", 150.0, 20.0, r1=11.6, b1=12.2),
    ])
    assert [ps.plate_search_status(by[k]) for k in "KLM"] == ["near_bright_star", "searched", "searched"]
    assert by["M"]["vanishing_data"]["gaia_match"] == "catalogue_epoch"
    cov = ps.search_coverage_summary(list(by.values()), cfg)
    assert cov["n_near_bright_star"] == 1 and cov["not_searched_ids"] == {"near_bright_star": ["USNO-B1.0 K"]}


def test_an_empty_field_is_believed_only_from_a_mirror_that_holds_the_catalogue(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []

    def empty_everywhere(url: str, params: dict[str, Any], timeout: float) -> _Resp:
        calls.append(url)
        return _Resp(text=_tsv([]))

    monkeypatch.setattr(ps.requests, "get", empty_everywhere)
    ps._usno_mirror_serves.clear()
    with pytest.raises(ArchiveAccessError, match="table not held"):
        ps.usno_cone(10.0, 40.0, 60.0, RawCache(tmp_path / "raw"))
    assert len(calls) == 2 * len(ls.VIZIER_TSV_URLS)  # the field and the canary position on each mirror
    ps._usno_mirror_serves.clear()


def test_live_search_refuses_only_the_radio_channel() -> None:
    assert set(CHANNELS_WITHOUT_LIVE_SOURCE) == set()


def test_live_vanishing_search_end_to_end(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """search -> vet -> null-report -> run record: the unit is a plate source and only compared sources count."""
    import shutil

    from oddstar import cli
    from oddstar.run_records import build_run_record

    stars, sky = _fields()
    sky.install(monkeypatch)
    enriched: list[str] = []
    monkeypatch.setattr(cli, "fetch_live_halo_bulge_sample", lambda cfg, cache, **kw: ([json.loads(json.dumps(s)) for s in stars], []))
    monkeypatch.setattr(cli, "enrich_star_live_for_vetting", lambda star, cfg, cache: (enriched.append(star["source_id"]) or dict(star), []))
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    monkeypatch.setattr("time.sleep", lambda s: None)
    run_dir = repo_root() / "runs" / "test-plate-survey"
    try:
        code, payload = dispatch(["search", "vanishing", "--live", "--region", "halo", "--limit", "5", "--run-id", "test-plate-survey", "--json"])
        assert code == 0, payload
        d = payload["data"]
        assert d["n_input"] == 8 and d["n_searched"] == 4 and d["flagged_ids"] == ["USNO-B1.0 B"] and enriched == []
        assert d["coverage"]["n_too_faint"] == 1 and d["coverage"]["n_no_coverage"] == 2 and d["coverage"]["n_failed"] == 1
        assert d["plate_fetch"]["n_fields_failed_usno"] == 1
        assert any("4 of 8 stars not searched" in w for w in payload["warnings"]) and any("1 of 5 fields returned no plate catalogue" in w for w in payload["warnings"])
        saved = json.loads((run_dir / "search_vanishing.json").read_text(encoding="utf-8"))
        assert [c["cut_name"] for c in saved["cuts"]] == ["gaia_halo_live", "usno_b1_plate_sources"] and saved["cuts"][1]["count_after"] == 8
        summaries = {s["source_id"]: s["summary"] for s in saved["scores"]}
        assert summaries["USNO-B1.0 S1"].startswith("NOT SEARCHED (no_coverage)") and summaries["USNO-B1.0 F"].startswith("NOT SEARCHED (too_faint)")
        kept = json.loads((run_dir / "live_stars.json").read_text(encoding="utf-8"))
        assert list(kept) == ["USNO-B1.0 B"] and kept["USNO-B1.0 B"]["field_source_id"] == "1"

        code, payload = dispatch(["vet", "test-plate-survey", "--json"])
        assert code == 0, payload
        vet = json.loads((run_dir / "vetting.json").read_text(encoding="utf-8"))["vetting"]["USNO-B1.0 B"]
        status = {v["check"]: v["status"] for v in vet}
        assert status["depth"] == "pass" and status["plate_defect"] == "not_evaluated" and status["sso_plate"] == "not_evaluated"

        code, payload = dispatch(["null-report", "test-plate-survey", "--json"])
        assert code == 0, payload
        rep = json.loads((run_dir / "null_report_vanishing.json").read_text(encoding="utf-8"))
        assert rep["n_searched"] == 4 and rep["n_flagged"] == 1

        rec, samples = build_run_record(run_dir)
        s = rec["searches"][0]
        assert s["coverage"]["n_searched"] == 4 and s["coverage"]["n_too_faint"] == 1 and s["plate_fetch"]["n_units"] == 8
        assert len(s["not_searched_ids"]) == 4 and samples[0]["n"] == 8
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)
