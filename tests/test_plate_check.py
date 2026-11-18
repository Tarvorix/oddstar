"""The plate_defect check on plate-image rows (waits for the limits of config 1.8.0)."""

from __future__ import annotations

import pytest

from oddstar.io import plate_images as pi

RA0, DEC0 = 150.0, 20.0


def _unit(rows, *, status="ok", errors=None, early=2, r1=17.5):
    return {"source_id": "USNO-B1.0 0900-0000001", "plate_source": True, "live_enriched": False,
            "usno": {"ra_deg": RA0, "dec_deg": DEC0, "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0, "position_epoch": 1975.0},
            "vanishing_data": {"early_plates_count": early, "usno_mag_r": r1,
                               "plate_images": {"measurement_version": pi.MEASUREMENT_VERSION, "rows": rows, "errors": errors or {}}},
            "vetting_archive": {"query_status": {"usno_b1": "ok", "plate_images": status}}}


def _row(survey, epoch, snr, offset, status="ok"):
    peak = None if snr is None else {"snr": snr, "offset_arcsec": offset, "dx_arcsec": -offset, "dy_arcsec": 0.0, "fwhm_arcsec": 4.0}
    return {"survey": survey, "survey_returned": survey.upper(), "plate_id": "X", "plate_epoch": epoch, "status": status,
            "ra_measured_deg": RA0, "dec_measured_deg": DEC0, "snr_at_position": snr, "peak": peak, "field_fwhm_median_arcsec": 4.1}


def test_plate_defect_reads_the_plates_for_a_survey_source(tmp_path) -> None:
    """A catalogue row is not enough for a live plate source: the first-epoch plate images decide."""
    from oddstar.config import load_config
    from oddstar.vet.vanishing_checks import check_plate_defect, plate_image_verdict

    cfg = load_config()[0]
    vc = cfg.channels["vanishing"]
    snr, off = float(vc["plate_image_min_snr"]), float(vc["plate_image_max_offset_arcsec"])
    assert float(vc["plate_image_first_epoch_before"]) == 1960.0
    good, far, faint = (snr + 3.0, off - 0.5), (snr + 30.0, off + 1.0), (snr - 1.0, 0.2)
    assert plate_image_verdict(_row("poss1_red", 1953.0, *good), cfg) is True
    assert plate_image_verdict(_row("poss1_red", 1953.0, *far), cfg) is False, "light beside the position is not the source"
    assert plate_image_verdict(_row("poss1_red", 1953.0, *faint), cfg) is False
    assert plate_image_verdict(_row("poss1_red", 1953.0, None, None), cfg) is False
    assert plate_image_verdict(_row("poss1_red", 1953.0, None, None, status="off_cutout"), cfg) is None
    # An image that could not be lined up on Gaia stars can show the source, but not its absence
    assert plate_image_verdict({**_row("poss1_red", 1953.0, *good), "registered": False}, cfg) is True
    assert plate_image_verdict({**_row("poss1_red", 1953.0, *far), "registered": False}, cfg) is None
    assert plate_image_verdict({**_row("poss1_red", 1953.0, None, None), "registered": True}, cfg) is False

    def run(rows, **kw):
        r = check_plate_defect(_unit(rows, **kw), cfg, tmp_path / f"v{len(list(tmp_path.iterdir()))}")
        return r.status, r.note, r.metrics

    later = [_row("poss2ukstu_red", 1991.0, *good), _row("poss2ukstu_blue", 1989.0, None, None)]
    # Both first-epoch plates show it: a real source of the 1950s sky
    status, note, m = run([_row("poss1_red", 1953.0, *good), _row("poss1_blue", 1953.0, *good), *later])
    assert status == "pass" and m["n_first_epoch_seen"] == 2 and m["n_later_seen"] == 1 and m["n_later_images"] == 2
    assert "point of light at the position on 2 of 2" in note and "later plates: seen on 1 of 2" in note
    # The plates were read and do not hold it: a catalogue entry, not a source
    status, note, m = run([_row("poss1_red", 1953.0, *far), _row("poss1_blue", 1953.0, None, None), *later])
    assert status == "fail" and m["n_first_epoch_seen"] == 0 and "A catalogue entry the plates do not hold" in note
    status, note, m = run([_row("poss1_red", 1953.0, *good), _row("poss1_blue", 1953.0, *faint)])
    assert status == "not_evaluated" and m["n_first_epoch_seen"] == 1 and "disagree" in note, "one plate is not two, and not none"
    # An image that could not be registered and is empty does not count towards "none shows it"
    status, note, m = run([{**_row("poss1_red", 1953.0, None, None), "registered": False}, _row("poss1_blue", 1953.0, None, None)])
    assert status == "not_evaluated" and m["n_first_epoch_images"] == 1 and "Only 1 first-epoch plate image" in note
    # A source so faint that ordinary stars are missed on these plates: a missing image is not evidence
    judge = float(vc["plate_image_judge_max_r_mag"])
    status, note, m = run([_row("poss1_red", 1953.0, *far), _row("poss1_blue", 1953.0, None, None)], r1=judge + 0.3)
    assert status == "not_evaluated" and "a missing image is not evidence" in note and m["plate_image_judge_max_r_mag"] == judge
    status, _, _ = run([_row("poss1_red", 1953.0, *good), _row("poss1_blue", 1953.0, *good)], r1=judge + 0.3)
    assert status == "pass", "a faint source the plates do show is a source"
    # Not read, or not measurable: no verdict
    status, note, _ = run([_row("poss1_red", 1953.0, *good)], status="failed: poss1_blue: HTTP 503", errors={"poss1_blue": "HTTP 503"})
    assert status == "not_evaluated" and "plate_images" in note
    status, note, _ = run([_row("poss1_red", 1953.0, *good), _row("poss1_blue", 1953.0, None, None, status="off_cutout")])
    assert status == "not_evaluated" and "Only 1 first-epoch plate image" in note
    # A plate returned for a first-epoch request that is not a first-epoch plate does not count as one
    status, _, m = run([_row("poss1_red", 1953.0, *good), _row("poss1_blue", 1978.0, *good)])
    assert status == "not_evaluated" and m["n_first_epoch_images"] == 1 and m["n_later_seen"] == 1
    # The catalogue count still fails first, without reading anything
    status, note, _ = run([], early=1)
    assert status == "fail" and "Early-epoch detections on 1" in note


def _live_unit(rows, neighbours, *, nb_status="ok", n_in_radius=None, r1=17.5, place_status="ok", gaia=(), ps1=()):
    u = _unit(rows, r1=r1)
    u["live_enriched"] = True
    u.update(ra_deg=RA0, dec_deg=DEC0, pmra_mas_yr=0.0, pmdec_mas_yr=0.0)
    u["vanishing_data"].update({
        "usno_epoch": 1953.0, "usno_ra_deg": RA0, "usno_dec_deg": DEC0, "present_in_gaia": False, "present_in_panstarrs": False,
        "nearest_gaia_track": {"n_gaia_stars": 4, "sep_arcsec": 41.2, "source_id": "1", "epoch": 1953.0, "plate_era": [1949.0, 2002.0]}})
    u["vetting_archive"]["query_status"].update({"modern_counterpart": "ok", "plate_neighbours": nb_status, "plate_first_place": place_status})
    if place_status == "ok":
        u["vanishing_data"]["plate_first_place"] = {"position": [RA0, DEC0], "epoch": 1953.0, "radius_arcsec": 10.0, "gaia": list(gaia),
                                                    "ps1": list(ps1)}
    if neighbours is not None:
        u["vanishing_data"]["plate_neighbours"] = {"radius_arcsec": 10.0, "inner_radius_arcsec": 3.0, "neighbours": neighbours,
                                                   "n_objects_in_radius": len(neighbours) if n_in_radius is None else n_in_radius}
    return u


def test_a_faint_star_that_may_have_moved_is_not_passed_as_vanished(tmp_path) -> None:
    """No Gaia star to carry a proper motion, but Pan-STARRS may hold the star a few arcsec on."""
    from oddstar.config import load_config
    from oddstar.vet.vanishing_checks import check_proper_motion

    cfg = load_config()[0]
    vc = cfg.channels["vanishing"]
    snr, off, dmag = float(vc["plate_image_min_snr"]), float(vc["plate_image_max_offset_arcsec"]), float(vc["plate_moved_max_dmag"])
    seen, empty = (snr + 5.0, off - 0.5), (None, None)
    source = [_row("poss1_red", 1953.0, *seen), _row("poss1_blue", 1953.0, *seen)]

    def nb(sep, r_mag, a, b, objid="1"):
        return {"objID": objid, "sep_arcsec": sep, "r_mag": r_mag, "n_detections": 40,
                "rows": [_row("poss1_red", 1953.0, *a), _row("poss1_blue", 1953.0, *b)]}

    def run(unit):
        r = check_proper_motion(unit, cfg, tmp_path / f"v{len(list(tmp_path.iterdir()))}")
        return r.status, r.note, r.metrics

    # Same brightness, 5.5 arcsec away now, and not there in 1953: it may be the star, moved (about 92 mas/yr)
    status, note, m = run(_live_unit(source, [nb(5.5, 17.9, empty, empty)]))
    assert status == "not_evaluated" and "it may be the plate source, moved" in note
    assert m["moved_star_test"]["neighbours"][0]["kind"] == "could_be_the_source_moved" and m["missing_inputs"] == ["plate_neighbours"]
    # A neighbour the 1953 plates show where it is now is another star, and so is one of another brightness
    status, note, m = run(_live_unit(source, [nb(5.5, 17.9, seen, empty), nb(7.0, 17.5 + dmag + 0.2, empty, empty, "2")]))
    assert status == "pass" and [t["kind"] for t in m["moved_star_test"]["neighbours"]] == ["on_the_plates_where_it_is_now", "other_brightness"]
    status, note, _ = run(_live_unit(source, []))
    assert status == "pass" and "none in the radius" in note
    # What was not read gives no verdict: the neighbours, a neighbour without a magnitude, more objects than measured
    for unit in (_live_unit(source, None, nb_status="failed: Pan-STARRS: HTTP 503"),
                 _live_unit(source, [nb(5.5, None, empty, empty)]),
                 _live_unit(source, [nb(5.5, 17.5 + dmag + 1.0, empty, empty)], n_in_radius=6)):
        status, note, m = run(unit)
        assert status == "not_evaluated" and m["missing_inputs"] == ["plate_neighbours"], note
    # Nothing on the plates at the source's own position: no star to have moved, and plate_defect says what the entry is
    status, note, _ = run(_live_unit([_row("poss1_red", 1953.0, *empty), _row("poss1_blue", 1953.0, *empty)], [nb(5.5, 17.9, empty, empty)]))
    assert status == "pass" and "no plate image of a star that could have moved" in note
    # A Gaia star whose path crosses the plate position still decides first
    unit = _live_unit(source, [])
    unit["vanishing_data"]["nearest_gaia_track"]["sep_arcsec"] = 1.1
    assert run(unit)[0] == "fail"


def _sky_row(survey, epoch, east, north, *, snr=20.0, registered=True):
    """A row whose peak lies (east, north) arcsec from (RA0, DEC0); the row was measured at (RA0, DEC0)."""
    return {"survey": survey, "survey_returned": survey.upper(), "plate_id": "X", "plate_epoch": epoch, "status": "ok", "registered": registered,
            "ra_measured_deg": RA0, "dec_measured_deg": DEC0, "snr_at_position": snr,
            "peak": {"snr": snr, "offset_arcsec": round((east**2 + north**2) ** 0.5, 2), "dx_arcsec": -east, "dy_arcsec": north, "fwhm_arcsec": 4.0}}


def test_the_later_plates_say_whether_the_source_moved_to_a_nearby_object(tmp_path) -> None:
    """Three epochs on a line: 1953 at the plate position, the 1990s part of the way, Pan-STARRS at the end."""
    import math

    from oddstar.config import load_config
    from oddstar.vet.vanishing_checks import PANSTARRS_EPOCH, check_proper_motion, plate_track

    cfg = load_config()[0]
    vc = cfg.channels["vanishing"]
    tol = float(vc["plate_track_max_residual_arcsec"])
    assert int(vc["plate_track_min_images"]) == 2
    # The object is 6 arcsec east and 2 north of where the 1953 plates show the source
    obj = {"objID": "7", "ra_deg": RA0 + 6.0 / 3600.0 / math.cos(math.radians(DEC0)), "dec_deg": DEC0 + 2.0 / 3600.0, "sep_arcsec": 6.32,
           "r_mag": 17.8, "n_detections": 40, "rows": [_row("poss1_red", 1953.0, None, None), _row("poss1_blue", 1953.0, None, None)]}
    first = [_sky_row("poss1_red", 1953.0, 0.2, -0.1), _sky_row("poss1_blue", 1953.0, -0.2, 0.1)]

    def at(epoch, slip=0.0):
        f = (epoch - 1953.0) / (PANSTARRS_EPOCH - 1953.0)
        return 6.0 * f + slip, 2.0 * f

    on_the_way = [_sky_row("poss2ukstu_red", 1990.0, *at(1990.0, 0.3)), _sky_row("poss2ukstu_blue", 1987.0, *at(1987.0, -0.4)),
                  _sky_row("poss2ukstu_ir", 1996.0, *at(1996.0))]
    in_place = [_sky_row("poss2ukstu_red", 1990.0, 0.3, 0.2), _sky_row("poss2ukstu_blue", 1987.0, -0.1, -0.4), _sky_row("poss2ukstu_ir", 1996.0, 0.0, 0.3)]

    tr = plate_track(first + on_the_way, obj["ra_deg"], obj["dec_deg"], cfg)
    assert tr["kind"] == "moved" and tr["n_on_line"] == 3 and tr["n_stayed"] == 0 and tr["first_to_now_arcsec"] == pytest.approx(6.32, abs=0.05)
    assert all(row["d_line_arcsec"] <= tol < row["d_stay_arcsec"] for row in tr["later"])
    tr = plate_track(first + in_place, obj["ra_deg"], obj["dec_deg"], cfg)
    assert tr["kind"] == "stayed" and tr["n_on_line"] == 0 and tr["n_stayed"] == 3
    # One later image is not enough, an unregistered one does not count, and a short line cannot tell the two apart
    assert plate_track(first + on_the_way[:1], obj["ra_deg"], obj["dec_deg"], cfg)["kind"] == "open"
    unreg = [dict(r, registered=False) for r in on_the_way]
    assert plate_track(first + unreg, obj["ra_deg"], obj["dec_deg"], cfg)["kind"] == "open"
    near_ra = RA0 + 1.5 / 3600.0 / math.cos(math.radians(DEC0))
    assert plate_track(first + in_place, near_ra, DEC0, cfg)["kind"] == "too_close"
    assert plate_track(on_the_way, obj["ra_deg"], obj["dec_deg"], cfg)["kind"] == "open", "no first-epoch position, no line"

    def run(rows, nbs):
        r = check_proper_motion(_live_unit(rows, nbs), cfg, tmp_path / f"v{len(list(tmp_path.iterdir()))}")
        return r.status, r.note, r.metrics

    status, note, m = run(first + on_the_way, [obj])
    assert status == "fail" and "the same star, moved" in note and "about 107 mas/yr" in note
    assert m["moved_star_test"]["neighbours"][0]["kind"] == "the_source_moved"
    status, note, m = run(first + in_place, [obj])
    assert status == "pass" and m["moved_star_test"]["neighbours"][0]["kind"] == "source_stayed_in_place"
    status, note, m = run(first + on_the_way[:1], [obj])
    assert status == "not_evaluated" and "later plates: 1 later image(s) usable" in note


def test_an_object_at_the_place_the_first_plates_show_means_the_source_never_left(tmp_path) -> None:
    """A false catalogue proper motion puts the catalogue position arcsec away from an object that stayed put."""
    import math

    from oddstar.config import load_config
    from oddstar.vet.vanishing_checks import check_proper_motion, first_place_counterpart

    cfg = load_config()[0]
    match = float(cfg.channels["vanishing"]["match_radius_arcsec"])
    cos_dec = math.cos(math.radians(DEC0))
    # The 1953 plates show the source 1.2 arcsec east and 0.8 south of where they were measured
    first = [_sky_row("poss1_red", 1953.0, 1.3, -0.9), _sky_row("poss1_blue", 1953.0, 1.1, -0.7)]

    def at(east, north):
        return RA0 + east / 3600.0 / cos_dec, DEC0 + north / 3600.0

    def gaia(east, north, sid="77"):
        ra, dec = at(east, north)
        return {"source_id": sid, "ra_then_deg": ra, "dec_then_deg": dec, "sep_then_arcsec": 1.0, "pmra_mas_yr": 30.0, "pmdec_mas_yr": -40.0, "g_mag": 17.6}

    def ps1(east, north, objid="88", r_mag=19.7):
        ra, dec = at(east, north)
        return {"objID": objid, "ra_deg": ra, "dec_deg": dec, "sep_arcsec": 1.0, "n_detections": 60, "r_mag": r_mag}

    def run(unit):
        r = check_proper_motion(unit, cfg, tmp_path / f"v{len(list(tmp_path.iterdir()))}")
        return r.status, r.note, r.metrics

    # A Gaia star that was there in 1953
    status, note, m = run(_live_unit(first, [], gaia=[gaia(1.4, -0.6)]))
    assert status == "fail" and "Gaia DR3 77, G 17.60 was 0.28\" from the place" in note and "is still there" in note
    assert m["first_place_test"]["place"] == {"east_arcsec": 1.2, "north_arcsec": -0.8, "of_position": [RA0, DEC0], "epoch": 1953.0, "n_images": 2}
    assert m["first_place_test"]["counterparts"][0]["pm_mas_yr"] == 50.0 and "moved_star_test" not in m
    # A Pan-STARRS object there now, with or without an r magnitude; the nearest is named
    status, note, m = run(_live_unit(first, [], ps1=[ps1(1.2, -0.8 + match - 0.2, "89", None), ps1(1.0, -1.0)]))
    assert status == "fail" and "Pan-STARRS 88, r 19.70 is 0.28\" from the place" in note
    assert [c["id"] for c in m["first_place_test"]["counterparts"]] == ["88", "89"]
    # Objects outside the survey's match radius of that place are not counterparts, also when they are near the measured position
    far = _live_unit(first, [], gaia=[gaia(1.2 - match - 0.3, -0.8)], ps1=[ps1(1.2, -0.8 + match + 0.3)])
    status, note, m = run(far)
    assert status == "pass" and f"within {match:.1f}\" of the place where 2 first-epoch plate image(s) (1953) show the source" in note
    assert m["first_place_test"]["verdict"] is False and m["first_place_test"]["n_gaia_in_cone"] == 1
    # The catalogues at that place not read: no pass, whatever the neighbours say
    for st in ("failed: HTTP 503", "no_coverage", "not run"):
        status, note, m = run(_live_unit(first, [], place_status=st))
        assert status == "not_evaluated" and m["missing_inputs"] == ["plate_first_place"] and "Not settled at the place" in note, st
    # No first-epoch image shows the source: there is no place to compare at, and the test does not stand in the way
    empty = [_row("poss1_red", 1953.0, None, None), _row("poss1_blue", 1953.0, None, None)]
    res = first_place_counterpart(_live_unit(empty, [], place_status="failed: HTTP 503"), cfg)
    assert res["verdict"] is False and "no measured place" in res["note"]
    # One first-epoch image is enough for a place
    one = _live_unit(first[:1], [], ps1=[ps1(1.3, -0.9)])
    assert first_place_counterpart(one, cfg)["place"]["n_images"] == 1 and run(one)[0] == "fail"
