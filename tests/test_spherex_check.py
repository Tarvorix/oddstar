"""The spherex_colours check: the colour rule of docs/spherex.md, its comparison stars and the stored record."""

from __future__ import annotations

import json

import pytest

from oddstar.config import load_config
from oddstar.io import spherex as sx
from oddstar.io.cache import RawCache
from oddstar.vet.ir_checks import check_spherex_colours, spherex_colour_verdict


def _rec(sid, ks, w1, w2, k, s1, s2, *, complete=True, n=(25, 24, 76), err=0.003):
    """A colour record from catalogue Ks, W1, W2 and SPHEREx window magnitudes k, w1, w2."""
    return {"source_id": str(sid), "catalogue": {"Ks": ks, "W1": w1, "W2": w2},
            "windows": {"k": {"n": n[0], "mag": k}, "w1": {"n": n[1], "mag": s1}, "w2": {"n": n[2], "mag": s2}},
            "k_w2": round(k - s2, 4), "w1_w2": round(s1 - s2, 4), "k_w2_err": err, "w1_w2_err": err, "complete": complete}


def _comparison(n=12, *, complete=None, spread=0.01):
    """Bare photospheres: SPHEREx minus catalogue is +0.07 in K - W2 and +0.11 in W1 - W2, with a little scatter."""
    out = []
    for i in range(n):
        d = spread * ((i % 5) - 2) / 2.0
        ks, w1, w2 = 10.50 + 0.05 * i, 10.44 + 0.05 * i, 10.49 + 0.05 * i
        # the catalogue colours are Ks - W2 = 0.01 and W1 - W2 = -0.05
        out.append(_rec(1000 + i, ks, w1, w2, w2 + 0.01 + 0.07 + d, w2 - 0.05 + 0.11 - d, w2,
                        complete=True if complete is None else complete[i]))
    return out


def _star(d_kw2, d_w1w2, **kw):
    """A star whose catalogue colours are Ks - W2 = 0.70 and W1 - W2 = 0.39, measured by SPHEREx with these differences."""
    return _rec("42", 10.97, 10.66, 10.27, 9.80 + 0.70 + d_kw2, 9.80 + 0.39 + d_w1w2, 9.80, **kw)


def test_the_rule_gives_the_three_verdicts_and_says_why_when_it_gives_none() -> None:
    cfg = load_config()[0]
    v = cfg.vetting
    assert (v["spherex_z_limit"], v["spherex_max_scatter_mag"], v["spherex_min_comparison_stars"], v["spherex_min_points_per_window"]) == (3.0, 0.10, 8, 5)
    comp = _comparison()
    same = spherex_colour_verdict(_star(0.07, 0.11), comp, cfg)
    assert same["verdict"] == "reproduced" and same["n_comparison_usable"] == 12
    assert same["colours"]["k_w2"]["median"] == pytest.approx(0.07, abs=0.003) and abs(same["colours"]["k_w2"]["z"]) < 1.0
    assert same["colours"]["k_w2"]["catalogue"] == pytest.approx(0.70, abs=1e-3)
    # Redder in one colour and unchanged in the other: a change, not a confirmation
    red = spherex_colour_verdict(_star(0.41, 0.13), comp, cfg)
    assert red["verdict"] == "redder_now" and red["colours"]["k_w2"]["z"] > 10 and abs(red["colours"]["w1_w2"]["z"]) < 3
    # Bluer in either colour decides, whatever the other does
    assert spherex_colour_verdict(_star(-0.30, 0.11), comp, cfg)["verdict"] == "not_reproduced"
    assert spherex_colour_verdict(_star(0.41, -0.20), comp, cfg)["verdict"] == "not_reproduced"
    # Comparison stars count only with a complete run and all three windows
    flags = [True] * 7 + [False] * 5
    few = spherex_colour_verdict(_star(0.41, 0.13), _comparison(complete=flags), cfg)
    assert few["verdict"] is None and few["reason"].startswith("7 of 12 comparison star(s)") and few["n_comparison_usable"] == 7
    eight = spherex_colour_verdict(_star(0.41, 0.13), _comparison(complete=[True] * 8 + [False] * 4), cfg)
    assert eight["verdict"] == "redder_now" and eight["comparison_ids"] == [str(1000 + i) for i in range(8)]
    holed = _comparison()
    for c in holed[:5]:
        c["windows"]["w1"] = None
        c["w1_w2"] = None
    assert spherex_colour_verdict(_star(0.41, 0.13), holed, cfg)["n_comparison_usable"] == 7
    # Too much scatter among the comparison stars, and too few measurements of the star
    wide = spherex_colour_verdict(_star(0.41, 0.13), _comparison(spread=0.5), cfg)
    assert wide["verdict"] is None and "scatter by more than 0.10 mag" in wide["reason"]
    thin = spherex_colour_verdict(_star(0.41, 0.13, n=(25, 4, 76)), comp, cfg)
    assert thin["verdict"] is None and "fewer than 5 SPHEREx measurements in window(s) w1" in thin["reason"]
    assert spherex_colour_verdict(dict(_star(0.41, 0.13), catalogue={"Ks": None, "W1": 10.66, "W2": 10.27}), comp, cfg)["verdict"] is None


def test_the_numbers_of_the_first_star_give_the_verdict_written_into_the_report() -> None:
    """Gaia DR3 5628404318169472 against the twelve comparison stars of 2026-11-18 (docs/spherex.md): redder now."""
    cfg = load_config()[0]
    rows = [(10.168, 10.105, 10.144, 10.186, 10.143, 10.090), (10.253, 10.205, 10.266, 10.278, 10.259, 10.209),
            (10.779, 10.718, 10.797, 10.785, 10.834, 10.725), (10.567, 10.459, 10.523, 10.553, 10.530, 10.462),
            (10.701, 10.612, 10.677, 10.703, 10.675, 10.628), (11.250, 11.199, 11.232, 11.272, 11.245, 11.168),
            (11.027, 10.967, 11.047, 11.028, 11.017, 10.974), (10.927, 10.841, 10.924, 10.926, 10.911, 10.851),
            (10.034, 10.005, 10.059, 10.098, 10.054, 10.015), (10.378, 10.316, 10.371, 10.385, 10.364, 10.326),
            (10.418, 10.353, 10.394, 10.431, 10.401, 10.354), (11.155, 11.105, 11.157, 11.165, 11.121, 11.084)]
    comp = [_rec(i, *r) for i, r in enumerate(rows)]
    star = dict(_rec("5628404318169472", 10.971, 10.665, 10.276, 10.914, 10.335, 9.810, err=0.0032), k_w2=1.1039, w1_w2=0.5252)
    res = spherex_colour_verdict(star, comp, cfg)
    assert res["verdict"] == "redder_now"
    assert res["colours"]["k_w2"]["z"] == pytest.approx(19.5, abs=0.5) and res["colours"]["w1_w2"]["z"] == pytest.approx(0.97, abs=0.1)
    assert res["colours"]["k_w2"]["median"] == pytest.approx(0.073, abs=0.002) and res["colours"]["w1_w2"]["scatter"] == pytest.approx(0.026, abs=0.002)


def _live(record, status="ok"):
    star = {"source_id": "42", "live_enriched": True, "vetting_archive": {"query_status": {"spherex_colours": status}}}
    if record is not None:
        star["vetting_archive"]["spherex_colours"] = record
    return star


def test_the_check_passes_fails_or_names_what_is_missing(tmp_path) -> None:
    cfg = load_config()[0]
    comp = _comparison()

    def run(star):
        r = check_spherex_colours(star, cfg, tmp_path / f"v{len(list(tmp_path.iterdir()))}")
        return r.status, r.note, r.metrics

    status, note, m = run(_live({"target": _star(0.07, 0.11), "comparison": comp, "as_of": "2026-11-18"}))
    assert status == "pass" and note.startswith("Reproduced: ") and "12 comparison stars, limit 3" in note and m["as_of"] == "2026-11-18"
    status, note, m = run(_live({"target": _star(0.41, 0.13), "comparison": comp}))
    assert status == "pass" and note.startswith("Redder now: ") and "a change, not a confirmation" in note and m["rule"]["verdict"] == "redder_now"
    status, note, _ = run(_live({"target": _star(-0.30, 0.11), "comparison": comp}))
    assert status == "fail" and note.startswith("Not reproduced: ")
    # No verdict from the rule, no stored record, and a record the query status does not vouch for
    status, note, _ = run(_live({"target": _star(0.41, 0.13), "comparison": comp[:5]}))
    assert status == "not_evaluated" and "5 of 5 comparison star(s)" in note
    status, note, m = run(_live(None, status="not_measured"))
    assert status == "not_evaluated" and "oddstar spherex verdict" in note and m["query_status"] == {"spherex_colours": "not_measured"}
    assert run(_live({"target": _star(0.07, 0.11), "comparison": comp}, status="not_measured"))[0] == "not_evaluated"
    # A fixture star (not live) without the record is not evaluated either, and one with it is judged
    assert run({"source_id": "fixture"})[0] == "not_evaluated"
    assert run({"source_id": "fixture", "vetting_archive": {"spherex_colours": {"target": _star(0.07, 0.11), "comparison": comp}}})[0] == "pass"


def _photometry(k=28.7, w1=22.7, w2=20.5, *, complete=True):
    def avg(flux, n, zp):
        import math

        return {"n": n, "flux_mjy": flux, "flux_err_mjy": 0.05, "mag": -2.5 * math.log10(flux * 1e-3 / zp), "mag_err": 2.5 / math.log(10.0) * 0.05 / flux}

    block = {"n": 125, "mjd_min": 60886.1, "mjd_max": 61287.7, "k": avg(k, 25, 666.7), "w1": avg(w1, 24, 309.54), "w2": avg(w2, 76, 171.787)}
    return {"n_images": 230, "n_measured": 210, "n_failed": 0 if complete else 2, "complete": complete, "summary": {"all": block, "visits": [], "n_misfit": 1}}


def test_colour_record_keeps_what_the_rule_needs() -> None:
    star = {"source_id": 7, "ra_deg": 43.5, "dec_deg": 5.9, "mags": {"BP": 14.041, "RP": 12.687, "Ks": 10.971, "W1": 10.665, "W2": 10.276}}
    rec = sx.colour_record(star, _photometry())
    assert rec["source_id"] == "7" and rec["bp_rp"] == pytest.approx(1.354) and rec["catalogue"] == {"Ks": 10.971, "W1": 10.665, "W2": 10.276}
    assert rec["windows"]["k"]["n"] == 25 and rec["windows"]["k"]["mag"] == pytest.approx(10.914, abs=0.003)
    assert rec["k_w2"] == pytest.approx(rec["windows"]["k"]["mag"] - rec["windows"]["w2"]["mag"], abs=1e-4)
    assert rec["k_w2_err"] == pytest.approx((rec["windows"]["k"]["mag_err"] ** 2 + rec["windows"]["w2"]["mag_err"] ** 2) ** 0.5, abs=2e-4)
    assert rec["complete"] and rec["n_failed"] == 0 and rec["n_misfit"] == 1 and rec["mjd_max"] == 61287.7
    # A window with no measurement leaves the colours that need it empty
    res = _photometry()
    res["summary"]["all"]["w1"] = None
    hole = sx.colour_record(star, res)
    assert hole["windows"]["w1"] is None and hole["w1_w2"] is None and hole["k_w2"] is not None


def test_comparison_stars_are_the_nearest_of_like_colour_and_brightness(monkeypatch, tmp_path) -> None:
    cfg = load_config()[0]
    asked: dict = {"cones": [], "batches": []}

    def cand(sid, east_deg, *, ph="AAAA", cc="0000", ext=0, clean=True):
        return {"source_id": str(sid), "ra_deg": 43.5 + east_deg, "dec_deg": 0.0, "mags": {"BP": 14.0, "RP": 12.7, "W1": 10.6},
                "allwise": {"ph_qual": ph, "cc_flags": cc, "ext_flg": ext}, "2mass_clean": clean}

    # Sources with an AllWISE and 2MASS record in the W1 window: 14 good ones 0.1 deg apart, the star itself, and five near ones
    records = {str(s["source_id"]): s for s in
               [cand(100 + i, 0.1 * (i + 1)) for i in range(14)]
               + [cand(1, 0.0), cand(2, 0.01, ph="ABAA"), cand(3, 0.02, cc="0d00"), cand(4, 0.03, ext=1), cand(5, 0.04, clean=False), cand(6, 0.05, ph="AACC")]}
    # Gaia sources of the right colour: those, and 30 faint ones (no record in the W1 window) spread over the first degree
    gaia = {sid: s["ra_deg"] - 43.5 for sid, s in records.items()} | {str(900 + i): 0.033 * (i + 1) - 0.001 for i in range(30)}

    def fake_tap(*, urls, query, archive, release, cache, timeout_s=45, sleep=None, **kw):
        radius = float(query.split("CIRCLE('ICRS', 43.500000, 0.000000, ")[1].split(")")[0])
        asked["cones"].append((radius, release, query))
        return [{"source_id": int(sid), "ra": 43.5 + east, "dec": 0.0} for sid, east in gaia.items() if east <= radius + 1e-9], None, False

    def fake_catalogue(where_sql, cfg, cache, *, release, sleep=None):
        ids = where_sql.split("IN (")[1].split(")")[0].split(",")
        asked["batches"].append((ids, where_sql, release))
        return [dict(records[i]) for i in ids if i in records]

    from oddstar.io import live_survey as ls
    monkeypatch.setattr(ls, "_tap_query_cached", fake_tap)
    monkeypatch.setattr(sx, "catalogue_stars", fake_catalogue)
    monkeypatch.setattr(sx, "COMPARISON_BATCH", 10)
    target = {"source_id": "1", "ra_deg": 43.5, "dec_deg": 0.0, "mags": {"BP": 14.041, "RP": 12.687, "W1": 10.665}}
    sel = sx.comparison_stars(target, cfg, RawCache(tmp_path / "raw"))
    # The star itself, a poor W2, a contaminated W2, an extended source and poor 2MASS photometry are left out; W3 and W4 do not matter
    assert [s["source_id"] for s in sel["stars"]] == ["6"] + [str(100 + i) for i in range(11)]
    assert sel["stars"][0]["comparison_sep_arcsec"] == pytest.approx(180.0, abs=0.5)
    assert sel["bp_rp_range"] == [1.154, 1.554] and sel["w1_range"] == [10.065, 11.265] and sel["radius_deg"] == 3.0
    # The cone grew from a quarter of the radius to half of it, and the stars were read nearest first, ten at a time
    assert [c[0] for c in asked["cones"]] == [0.75, 1.5] and all(c[1].endswith("_spherex_cone") for c in asked["cones"])
    assert "BETWEEN 1.1540 AND 1.5540" in asked["cones"][0][2] and "allwise" not in asked["cones"][0][2]
    assert all(len(b[0]) <= 10 and "1" not in b[0] for b in asked["batches"]) and asked["batches"][0][0][:5] == ["2", "3", "4", "900", "5"]
    assert all("aw.w1mpro BETWEEN 10.0650 AND 11.2650" in b[1] and b[2].endswith("_spherex_comparison") for b in asked["batches"])
    seen = [i for b in asked["batches"] for i in b[0]]
    assert len(seen) == len(set(seen)) == sel["n_examined"] and sel["cone_deg"] == 1.5
    # It stopped with the batch that brought the twelfth clean star: nothing nearer than the last star chosen was left unread
    last_sep_deg = sel["stars"][-1]["comparison_sep_arcsec"] / 3600.0
    assert sel["examined_out_to_deg"] >= last_sep_deg and all(i in seen for i, east in gaia.items() if east <= last_sep_deg and i != "1")
    assert sel["n_clean"] >= 12 and sel["n_matched"] >= sel["n_clean"] and sel["n_gaia_in_cone"] == sum(1 for e in gaia.values() if e <= 1.5)
    # The answer is the same when everything is read in one batch
    monkeypatch.setattr(sx, "COMPARISON_BATCH", 10_000)
    whole = sx.comparison_stars(target, cfg, RawCache(tmp_path / "raw2"))
    assert [s["source_id"] for s in whole["stars"]] == [s["source_id"] for s in sel["stars"]]
    # A Gaia list at the archive's row limit is not the whole cone and is refused
    monkeypatch.setattr(ls, "GAIA_SYNC_ROW_CAP", 5)
    with pytest.raises(Exception, match="row limit"):
        sx.comparison_stars(target, cfg, RawCache(tmp_path / "raw3"))
    monkeypatch.undo()
    with pytest.raises(Exception, match="lacks BP, RP or W1"):
        sx.comparison_stars({"source_id": "9", "ra_deg": 1.0, "dec_deg": 1.0, "mags": {"BP": 14.0}}, cfg, RawCache(tmp_path / "raw"))


def test_the_stored_record_reaches_vetting(tmp_path) -> None:
    cache = RawCache(tmp_path / "raw")
    star = {"source_id": "42", "vetting_archive": {"query_status": {"simbad": "ok"}}}
    assert sx.attach_colours(json.loads(json.dumps(star)), cache)["vetting_archive"]["query_status"]["spherex_colours"] == "not_measured"
    payload = {"target": _star(0.41, 0.13), "comparison": _comparison(), "as_of": "2026-11-18"}
    sx.store_colour_inputs(cache, "42", payload)
    got = sx.attach_colours(json.loads(json.dumps(star)), cache)
    assert got["vetting_archive"]["query_status"] == {"simbad": "ok", "spherex_colours": "ok"} and got["vetting_archive"]["spherex_colours"] == payload
    # A later measurement replaces the earlier one
    sx.store_colour_inputs(cache, "42", dict(payload, as_of="2027-01-01"))
    assert sx.stored_colour_inputs(cache, "42")["as_of"] == "2027-01-01" and sx.stored_colour_inputs(cache, "43") is None


def test_spherex_verdict_command_measures_stores_and_judges(monkeypatch, tmp_path) -> None:
    import shutil

    from oddstar import cli
    from oddstar.cli import dispatch
    from oddstar.config import repo_root

    target = {"source_id": "42", "ra_deg": 43.5, "dec_deg": 5.9, "pmra_mas_yr": 38.0, "pmdec_mas_yr": -97.0,
              "mags": {"G": 13.43, "BP": 14.041, "RP": 12.687, "Ks": 10.971, "W1": 10.665, "W2": 10.276}}
    comps = [{"source_id": str(500 + i), "ra_deg": 43.6 + 0.1 * i, "dec_deg": 5.9, "comparison_sep_arcsec": 360.0 * (i + 1),
              "mags": {"G": 13.0, "BP": 13.7, "RP": 12.4, "Ks": 10.50, "W1": 10.44, "W2": 10.49}} for i in range(12)]
    measured: list[tuple[str, tuple[int, ...]]] = []

    def fake_photometry(sid, ra, dec, cache, *, as_of, detectors, **kw):
        measured.append((sid, tuple(detectors)))
        if sid == "42":
            return _photometry()
        if sid == "511":
            raise RuntimeError("HTTP 503 from the archive")
        # bare photospheres of the catalogue colours, with the usual window offsets; one run short of two cutouts
        k = 666.7e3 * 10 ** (-0.4 * (10.49 + 0.01 + 0.07))
        w1 = 309.54e3 * 10 ** (-0.4 * (10.49 - 0.05 + 0.11))
        w2 = 171.787e3 * 10 ** (-0.4 * 10.49)
        return _photometry(k, w1, w2, complete=sid != "510")

    monkeypatch.setattr(sx, "target_star", lambda sid, cfg, cache, **kw: dict(target) if sid == "42" else None)
    monkeypatch.setattr(sx, "comparison_stars", lambda t, cfg, cache, **kw: {"stars": [dict(c) for c in comps], "n_examined": 30, "examined_out_to_deg": 0.5, "n_clean": 20,
                                                                           "radius_deg": 3.0, "bp_rp_range": [1.154, 1.554], "w1_range": [10.065, 11.265]})
    monkeypatch.setattr(sx, "gaia_neighbours", lambda *a, **kw: ())
    monkeypatch.setattr(sx, "forced_photometry", fake_photometry)
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    run_dir = repo_root() / "runs" / "test-spherex-verdict"
    try:
        code, payload = dispatch(["spherex", "verdict", "--ids", "42,43", "--as-of", "2026-11-18", "--run-id", "test-spherex-verdict", "--dry-run", "--json"])
        assert code == 0 and "[dry-run]" in payload["summary"] and not measured
        code, payload = dispatch(["spherex", "verdict", "--ids", "42,43", "--as-of", "2026-11-18", "--run-id", "test-spherex-verdict", "--json"])
        assert code == 0, payload
        first, second = payload["data"]["stars"]
        assert first["status"] == "measured" and first["verdict"] == "redder_now" and first["rule"]["n_comparison_usable"] == 10
        assert second == {"source_id": "43", "status": "no_catalogue_record"}
        assert [m[0] for m in measured] == ["42"] + [str(500 + i) for i in range(12)] and all(m[1] == (3, 4, 5, 6) for m in measured)
        assert first["comparison"][11]["error"].startswith("RuntimeError: HTTP 503") and first["comparison"][10]["complete"] is False
        assert first["comparison"][0]["comparison_sep_arcsec"] == 360.0 and first["selection"]["n_clean"] == 20
        assert json.loads((run_dir / "spherex_verdict.json").read_text())["stars"][0]["verdict"] == "redder_now"
        # What vetting reads afterwards
        stored = sx.stored_colour_inputs(RawCache(tmp_path / "raw"), "42")
        assert stored["as_of"] == "2026-11-18" and stored["detectors"] == [3, 4, 5, 6] and len(stored["comparison"]) == 12
        star = sx.attach_colours({"source_id": "42", "live_enriched": True}, RawCache(tmp_path / "raw"))
        cfg = load_config()[0]
        res = check_spherex_colours(star, cfg, tmp_path / "vet")
        assert res.status == "pass" and res.note.startswith("Redder now: ")
        code, payload = dispatch(["spherex", "verdict", "--ids", "abc", "--json"])
        assert code != 0 and payload["code"] == "USAGE_ERROR"
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)
