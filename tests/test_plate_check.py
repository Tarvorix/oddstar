"""The plate_defect check on plate-image rows (waits for the limits of config 1.8.0)."""

from __future__ import annotations

from oddstar.io import plate_images as pi

RA0, DEC0 = 150.0, 20.0


def _unit(rows, *, status="ok", errors=None, early=2, r1=17.5):
    return {"source_id": "USNO-B1.0 0900-0000001", "plate_source": True, "live_enriched": False,
            "usno": {"ra_deg": RA0, "dec_deg": DEC0, "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0, "position_epoch": 1975.0},
            "vanishing_data": {"early_plates_count": early, "usno_mag_r": r1,
                               "plate_images": {"measurement_version": pi.MEASUREMENT_VERSION, "rows": rows, "errors": errors or {}}},
            "vetting_archive": {"query_status": {"usno_b1": "ok", "plate_images": status}}}


def _row(survey, epoch, snr, offset, status="ok"):
    peak = None if snr is None else {"snr": snr, "offset_arcsec": offset, "fwhm_arcsec": 4.0}
    return {"survey": survey, "survey_returned": survey.upper(), "plate_id": "X", "plate_epoch": epoch, "status": status,
            "snr_at_position": snr, "peak": peak, "field_fwhm_median_arcsec": 4.1}


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


def _live_unit(rows, neighbours, *, nb_status="ok", n_in_radius=None, r1=17.5):
    u = _unit(rows, r1=r1)
    u["live_enriched"] = True
    u.update(ra_deg=RA0, dec_deg=DEC0, pmra_mas_yr=0.0, pmdec_mas_yr=0.0)
    u["vanishing_data"].update({
        "usno_epoch": 1953.0, "usno_ra_deg": RA0, "usno_dec_deg": DEC0, "present_in_gaia": False, "present_in_panstarrs": False,
        "nearest_gaia_track": {"n_gaia_stars": 4, "sep_arcsec": 41.2, "source_id": "1", "epoch": 1953.0, "plate_era": [1949.0, 2002.0]}})
    u["vetting_archive"]["query_status"].update({"modern_counterpart": "ok", "plate_neighbours": nb_status})
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
    assert status == "not_evaluated" and "it may be the plate source, moved" in note and "about 92 mas/yr" in note
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
