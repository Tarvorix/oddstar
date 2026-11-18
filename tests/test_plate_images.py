"""Plate-image measurement (`oddstar.io.plate_images`) on synthetic Digitized Sky Survey cutouts."""

from __future__ import annotations

import io
import math

import numpy as np
import pytest

from oddstar.errors import ArchiveAccessError
from oddstar.io import plate_images as pi
from oddstar.io.cache import RawCache

RA0, DEC0 = 150.0, 20.0


def _cutout(stars=(), *, pix=1.7, n=71, sky=4000.0, grain=60.0, fwhm=4.0, seed=0, date="1951-09-08T07:98:00", survey="POSSI-E",
            gradient=0.0) -> bytes:
    """A plate cutout as the service delivers it: int16 density, stars brighter than sky, a TAN sky solution.

    `stars` are (east offset arcsec, north offset arcsec, peak density above sky). The grain is smoothed over 3
    pixels, as plate grain is. `gradient` adds a slope in density per pixel along x (the wing of a bright star).
    """
    from astropy.io import fits
    from scipy.ndimage import uniform_filter

    rng = np.random.default_rng(seed)
    img = sky + uniform_filter(rng.normal(0.0, grain * 3.0, (n, n)), size=3, mode="wrap")
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    c = (n - 1) / 2.0
    img += gradient * (xx - c)
    sig = fwhm / 2.3548 / pix
    for east, north, peak in stars:
        # east is towards smaller x (CD1_1 negative)
        img += peak * np.exp(-((xx - (c - east / pix)) ** 2 + (yy - (c + north / pix)) ** 2) / (2.0 * sig**2))
    hdu = fits.PrimaryHDU(np.clip(np.round(img), -32000, 32000).astype(">i2"))
    h = hdu.header
    h["CTYPE1"], h["CTYPE2"] = "RA---TAN", "DEC--TAN"
    h["CRPIX1"], h["CRPIX2"] = c + 1.0, c + 1.0
    h["CRVAL1"], h["CRVAL2"] = RA0, DEC0
    h["CD1_1"], h["CD1_2"], h["CD2_1"], h["CD2_2"] = -pix / 3600.0, 0.0, 0.0, pix / 3600.0
    h["SURVEY"], h["PLATEID"], h["REGION"], h["DATE-OBS"], h["EXPOSURE"] = survey, "07XY", "XE123", date, 50.0
    h["EMULSION"], h["FILTER"], h["BANDPASS"], h["PLTGRADE"] = "103aE", "plexi", 8, "A"
    buf = io.BytesIO()
    hdu.writeto(buf)
    return buf.getvalue()


def test_plate_epoch_reads_the_date_and_ignores_an_invalid_time() -> None:
    assert pi.plate_epoch({"DATE-OBS": "1951-09-08T07:98:00"}) == pytest.approx(1951.685, abs=0.003)
    assert pi.plate_epoch({"DATE-OBS": "", "EPOCH": 1990.75}) == 1990.75
    assert pi.plate_epoch({"DATE-OBS": "garbage", "EPOCH": 2000.0}) == 2000.0
    assert pi.plate_epoch({"EPOCH": 1855.0}) is None and pi.plate_epoch({}) is None


def test_a_star_at_the_position_is_found_on_it() -> None:
    row = pi.measure_plate_cutout(_cutout([(0.8, -0.5, 900.0), (30.0, 25.0, 6000.0), (-35.0, 10.0, 5000.0)]), RA0, DEC0)
    assert row["status"] == "ok" and row["survey_returned"] == "POSSI-E" and row["plate_id"] == "07XY"
    assert row["pixel_arcsec"] == pytest.approx(1.7, abs=0.01) and row["plate_epoch"] == pytest.approx(1951.69, abs=0.01)
    assert row["snr_at_position"] > 8 and row["n_noise_apertures"] >= pi.MIN_NOISE_APERTURES
    p = row["peak"]
    assert p["offset_arcsec"] == pytest.approx(math.hypot(0.8, 0.5), abs=0.6) and p["snr"] > 8
    assert p["dx_arcsec"] < 0 < 1 and p["fwhm_arcsec"] == pytest.approx(4.0, abs=1.5)
    assert row["n_field_stars"] == 2 and row["field_fwhm_median_arcsec"] == pytest.approx(4.0, abs=1.5)


def test_an_empty_position_shows_nothing() -> None:
    row = pi.measure_plate_cutout(_cutout([(30.0, 25.0, 6000.0)], seed=3), RA0, DEC0)
    assert row["status"] == "ok" and abs(row["snr_at_position"]) < 3.5
    assert row["peak"] is None or row["peak"]["snr"] < 4


def test_light_beside_the_position_is_reported_where_it_is() -> None:
    # A star 6 arcsec away: the position itself is empty, and the nearest peak says how far off the light is
    row = pi.measure_plate_cutout(_cutout([(6.0, 0.0, 3000.0)], seed=4, pix=1.01, n=119), RA0, DEC0)
    assert row["status"] == "ok" and row["peak"]["offset_arcsec"] == pytest.approx(6.0, abs=0.6) and row["peak"]["snr"] > 10
    assert row["snr_at_position"] < row["peak"]["snr"] / 3
    # On the wing of a bright star outside the cutout: a slope, no point source, and the noise estimate carries the slope
    wing = pi.measure_plate_cutout(_cutout([], seed=5, gradient=40.0), RA0, DEC0)
    assert wing["status"] == "ok" and abs(wing["snr_at_position"]) < 3.5
    assert wing["peak"] is None or wing["peak"]["snr"] < 4


def test_the_position_is_moved_to_the_plate_date() -> None:
    # 500 mas/yr north over 48.3 years from 2000.0 back to 1951.7: the star was 24.2 arcsec south of its 2000 position
    payload = _cutout([(0.0, -24.16, 1500.0)], seed=6)
    still = pi.measure_plate_cutout(payload, RA0, DEC0)
    assert still["peak"] is None or still["peak"]["offset_arcsec"] > 3
    moved = pi.measure_plate_cutout(payload, RA0, DEC0, pmdec_mas_yr=500.0, position_epoch=2000.0)
    assert moved["status"] == "ok" and moved["peak"]["offset_arcsec"] < 1.0 and moved["snr_at_position"] > 8
    assert moved["dec_measured_deg"] == pytest.approx(DEC0 - 24.16 / 3600.0, abs=2e-5)
    undated = pi.measure_plate_cutout(_cutout([], date=""), RA0, DEC0, pmdec_mas_yr=500.0, position_epoch=2000.0)
    assert undated["status"] == "no_plate_date"


def test_a_position_at_the_edge_is_not_measured() -> None:
    row = pi.measure_plate_cutout(_cutout([]), RA0 + 50.0 / 3600.0 / math.cos(math.radians(DEC0)), DEC0)
    assert row["status"] == "off_cutout" and "peak" not in row


def _ref(east: float, north: float, pmra: float = 0.0, pmdec: float = 0.0) -> tuple[float, float, float, float]:
    """A reference star given by its true offset from the field centre (arcsec)."""
    return (RA0 + east / 3600.0 / math.cos(math.radians(DEC0)), DEC0 + north / 3600.0, pmra, pmdec)


FIELD = [(60.0, 40.0), (-70.0, 25.0), (35.0, -80.0), (-45.0, -60.0), (85.0, -15.0)]
SHIFT = (2.8, 1.5)  # where this plate's delivered sky solution puts every star, relative to where it is (east, north)


def _shifted_cutout(target_peak: float = 700.0, **kw) -> bytes:
    stars = [(e + SHIFT[0], n + SHIFT[1], 5000.0) for e, n in FIELD]
    if target_peak:
        stars.append((SHIFT[0], SHIFT[1], target_peak))
    return _cutout(stars, n=141, **kw)


def test_a_cutout_is_registered_on_its_reference_stars() -> None:
    payload = _shifted_cutout(seed=11)
    refs = tuple(_ref(e, n) for e, n in FIELD)
    as_delivered = pi.measure_plate_cutout(payload, RA0, DEC0)
    assert as_delivered["registered"] is False and as_delivered["registration"]["n_references"] == 0
    assert as_delivered["peak"]["offset_arcsec"] == pytest.approx(math.hypot(*SHIFT), abs=0.6)
    row = pi.measure_plate_cutout(payload, RA0, DEC0, references=refs)
    reg = row["registration"]
    assert row["registered"] is True and reg["n_stars"] == 5 and reg["scatter_arcsec"] < 0.5
    # east is towards smaller x
    assert reg["dx_arcsec"] == pytest.approx(-SHIFT[0], abs=0.4) and reg["dy_arcsec"] == pytest.approx(SHIFT[1], abs=0.4)
    assert row["peak"]["offset_arcsec"] < 0.7 and row["snr_at_position"] > as_delivered["snr_at_position"]
    assert row["x"] != row["x_delivered"]
    # An empty position stays empty: the registration is taken from other stars, never from the position itself
    empty = pi.measure_plate_cutout(_shifted_cutout(target_peak=0.0, seed=12), RA0, DEC0, references=refs)
    assert empty["registered"] is True and abs(empty["snr_at_position"]) < 3.5
    assert empty["peak"] is None or empty["peak"]["snr"] < 4


def test_reference_stars_move_to_the_plate_date() -> None:
    # Each reference star has moved 200 mas/yr east since the plate: 64.3 years from 1951.7 to 2016.0 is 12.9 arcsec
    moved = 0.2 * (2016.0 - 1951.685)
    refs = tuple(_ref(e + moved, n, pmra=200.0) for e, n in FIELD)
    row = pi.measure_plate_cutout(_shifted_cutout(seed=13), RA0, DEC0, references=refs)
    assert row["registered"] is True and row["registration"]["n_stars"] == 5 and row["peak"]["offset_arcsec"] < 0.7


def test_too_few_reference_stars_leave_the_position_as_delivered() -> None:
    payload = _shifted_cutout(seed=14)
    # The star at the position is not a reference for itself, a pair 6 arcsec apart is not used, one star has no image
    refs = (_ref(0.0, 0.0), _ref(*FIELD[0]), _ref(FIELD[0][0] + 6.0, FIELD[0][1]), _ref(*FIELD[1]), _ref(*FIELD[2]), _ref(-90.0, 90.0))
    row = pi.measure_plate_cutout(payload, RA0, DEC0, references=refs[:5])
    assert row["registered"] is False and row["registration"] == {"n_references": 5, "n_usable": 2, "n_stars": 2}
    assert row["x"] == row["x_delivered"] and row["peak"]["offset_arcsec"] == pytest.approx(math.hypot(*SHIFT), abs=0.6)
    three = pi.measure_plate_cutout(payload, RA0, DEC0, references=refs[3:] + (_ref(*FIELD[3]),))
    assert three["registered"] is True and three["registration"]["n_usable"] == 4 and three["registration"]["n_stars"] == 3


def test_reference_stars_are_chosen_by_brightness_and_are_part_of_what_a_stored_row_stands_for() -> None:
    gaia = [{"ra_deg": RA0, "dec_deg": DEC0, "pmra_mas_yr": 1.0, "pmdec_mas_yr": None, "g_mag": g} for g in (11.9, 12.0, 15.5, 19.0, 19.1, None)]
    refs = pi.reference_stars(gaia)
    assert len(refs) == 3 and refs[0] == (RA0, DEC0, 1.0, 0.0)
    assert pi._signature("poss1_red", RA0, DEC0, 0.0, 0.0, None)[0] != pi._signature("poss1_red", RA0, DEC0, 0.0, 0.0, None, refs)[0]


def test_rows_are_stored_and_a_bad_answer_is_not(monkeypatch, tmp_path) -> None:
    calls: list[str] = []

    class Resp:
        def __init__(self, content: bytes) -> None:
            self.content = content

    def fake_get(url, params, **kw):
        calls.append(url)
        if "poss1_blue" in url:
            return Resp(b"<html>No plates found</html>")
        if "poss2ukstu_ir" in url:
            raise ArchiveAccessError("GET failed after 3 attempt(s): HTTP 503", retryable=True)
        return Resp(_cutout([(0.0, 0.0, 900.0)], seed=7))

    monkeypatch.setattr(pi.ot, "_get_with_retries", fake_get)
    cache = RawCache(tmp_path / "raw")
    out = pi.measure_all_surveys(RA0, DEC0, cache, sleep=lambda s: None)
    assert [r["survey"] for r in out["rows"]] == ["poss1_red", "poss2ukstu_red", "poss2ukstu_blue"] and out["complete"] is False
    assert "not a FITS file" in out["errors"]["poss1_blue"] and "HTTP 503" in out["errors"]["poss2ukstu_ir"]
    assert len(calls) == 5 and "v=poss1_red&r=150.000000&d=20.000000&e=J2000&h=4.0&w=4.0&f=fits" in calls[0]
    # A second pass fetches only what has no stored row
    again = pi.measure_all_surveys(RA0, DEC0, cache, sleep=lambda s: None)
    assert again["rows"] == out["rows"] and len(calls) == 7
    with pytest.raises(ValueError, match="unknown plate survey"):
        pi.measure_position("sdss", RA0, DEC0, cache)


def test_enrich_unit_attaches_rows_and_says_what_was_read(monkeypatch, tmp_path) -> None:
    class Resp:
        def __init__(self, content: bytes) -> None:
            self.content = content

    def fake_get(url, params, **kw):
        if "poss1_blue" in url:
            raise ArchiveAccessError("GET failed after 3 attempt(s): HTTP 503", retryable=True)
        return Resp(_cutout([(0.0, 0.0, 900.0)], seed=8))

    monkeypatch.setattr(pi.ot, "_get_with_retries", fake_get)
    cache = RawCache(tmp_path / "raw")
    slept: list[float] = []
    unit = {"source_id": "USNO-B1.0 X", "plate_source": True, "usno": {"ra_deg": RA0, "dec_deg": DEC0, "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0,
                                                                      "position_epoch": 1975.0}, "vanishing_data": {"early_plates_count": 2}}
    out = pi.enrich_unit(unit, cache, min_interval_s=0.5, sleep=slept.append)
    assert out["vetting_archive"]["query_status"]["plate_images"].startswith("failed: poss1_blue: ")
    assert len(out["vanishing_data"]["plate_images"]["rows"]) == 4 and len(slept) == 5
    # Stored rows are not asked again and cost no wait
    slept.clear()
    again = pi.enrich_unit(dict(unit), cache, min_interval_s=0.5, sleep=slept.append)
    assert len(again["vanishing_data"]["plate_images"]["rows"]) == 4 and len(slept) == 1
    bare = pi.enrich_unit({"source_id": "x"}, cache)
    assert bare["vetting_archive"]["query_status"] == {"plate_images": "no_plate_record"}


def test_enrich_unit_registers_on_gaia_stars_when_given_a_config(monkeypatch, tmp_path) -> None:
    from oddstar.io import plate_survey

    class Resp:
        content = _shifted_cutout(seed=15)

    cones: list[float] = []

    down: list[bool] = []

    def fake_cone(ra, dec, radius_arcsec, cfg, cache, *, sleep=None):
        cones.append(radius_arcsec)
        if down:
            raise ArchiveAccessError("Gaia TAP failed on every mirror", retryable=True)
        rows = [{"source_id": str(i), "ra_deg": r[0], "dec_deg": r[1], "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0, "g_mag": 16.0}
                for i, r in enumerate(_ref(e, n) for e, n in FIELD)]
        return rows, "manifest", False

    from oddstar.config import load_config

    cfg = load_config()[0]
    monkeypatch.setattr(pi.ot, "_get_with_retries", lambda url, params, **kw: Resp())
    monkeypatch.setattr(plate_survey, "gaia_cone", fake_cone)
    monkeypatch.setattr(plate_survey, "ps1_cone", lambda *a, **kw: ([], "m", False))
    cache = RawCache(tmp_path / "raw")
    unit = {"source_id": "USNO-B1.0 Y", "plate_source": True, "usno": {"ra_deg": RA0, "dec_deg": DEC0}}
    out = pi.enrich_unit(json_copy(unit), cache, cfg=cfg, min_interval_s=0.0, sleep=lambda s: None)
    assert out["vanishing_data"]["plate_neighbours"]["n_objects_in_radius"] == 0
    res = out["vanishing_data"]["plate_images"]
    assert cones == [120.0, 120.0] and res["n_reference_stars"] == 5 and out["vetting_archive"]["query_status"]["plate_images"] == "ok"
    assert all(r["registered"] and r["peak"]["offset_arcsec"] < 1.2 for r in res["rows"])
    # The catalogues at the place the first-epoch plates were measured at (no reference star is within 10 arcsec of it)
    assert out["vetting_archive"]["query_status"]["plate_first_place"] == "ok"
    assert out["vanishing_data"]["plate_first_place"] == {"position": [RA0, DEC0], "epoch": res["rows"][0]["plate_epoch"],
                                                          "radius_arcsec": pi.FIRST_PLACE_CONE_ARCSEC, "gaia": [], "ps1": []}
    # Without the reference stars nothing is measured, and the status says why
    down.append(True)
    failed = pi.enrich_unit(json_copy(unit), cache, cfg=cfg, min_interval_s=0.0, sleep=lambda s: None)
    assert failed["vetting_archive"]["query_status"]["plate_images"].startswith("failed: reference stars: Gaia TAP failed")
    assert "plate_images" not in failed.get("vanishing_data", {})


def test_neighbours_are_measured_on_the_first_epoch_plates(monkeypatch, tmp_path) -> None:
    from oddstar.config import load_config
    from oddstar.io import plate_survey

    cfg = load_config()[0]
    asked: list[str] = []

    class Resp:
        content = _cutout([(0.0, 0.0, 900.0)], n=141, seed=16)

    def fake_get(url, params, **kw):
        asked.append(url)
        return Resp()

    def obj(objid, east, north, r_mag, n_det=30):
        ra, dec, _, _ = _ref(east, north)
        return {"objID": objid, "raMean": ra, "decMean": dec, "nDetections": n_det, "rMeanPSFMag": r_mag}

    cones: list[float] = []

    def fake_ps1(ra, dec, radius_arcsec, min_detections, cache, **kw):
        cones.append(radius_arcsec)
        # one inside the match radius (the source itself would be there), two in the ring, one without an r magnitude
        return [obj(11, 1.0, 0.5, 18.0), obj(12, 6.0, 0.0, 18.2), obj(13, 0.0, -8.0, -999.0), obj(14, 30.0, 0.0, 17.0)], "m", False

    monkeypatch.setattr(pi.ot, "_get_with_retries", fake_get)
    monkeypatch.setattr(plate_survey, "ps1_cone", fake_ps1)
    cache = RawCache(tmp_path / "raw")
    unit = {"source_id": "USNO-B1.0 Z", "plate_source": True, "usno": {"ra_deg": RA0, "dec_deg": DEC0, "plate_ra_deg": RA0, "plate_dec_deg": DEC0}}
    out = pi.attach_neighbours(unit, cfg, cache, (), min_interval_s=0.0, sleep=lambda s: None)
    nb = out["vanishing_data"]["plate_neighbours"]
    assert cones == [float(cfg.channels["vanishing"]["plate_moved_radius_arcsec"])] and nb["n_objects_in_radius"] == 2
    assert [(n["objID"], n["sep_arcsec"], n["r_mag"]) for n in nb["neighbours"]] == [("12", 6.0, 18.2), ("13", 8.0, None)]
    assert all([r["survey"] for r in n["rows"]] == list(pi.FIRST_EPOCH_SURVEYS) for n in nb["neighbours"]) and len(asked) == 4
    assert out["vetting_archive"]["query_status"]["plate_neighbours"] == "ok"
    south = pi.attach_neighbours({"usno": {"ra_deg": RA0, "dec_deg": -45.0}}, cfg, cache, ())
    assert south["vetting_archive"]["query_status"]["plate_neighbours"] == "no_coverage"

    def down(*a, **kw):
        raise ArchiveAccessError("HTTP 503", retryable=True)

    monkeypatch.setattr(plate_survey, "ps1_cone", down)
    failed = pi.attach_neighbours({"usno": {"ra_deg": RA0, "dec_deg": DEC0}}, cfg, cache, ())
    assert failed["vetting_archive"]["query_status"]["plate_neighbours"].startswith("failed: Pan-STARRS: HTTP 503")


def json_copy(obj):
    import json

    return json.loads(json.dumps(obj))


def test_first_place_reads_gaia_moved_back_and_panstarrs_around_the_measured_position(monkeypatch, tmp_path) -> None:
    from oddstar.config import load_config
    from oddstar.io import plate_survey

    cfg = load_config()[0]
    # The catalogue position (epoch 2000) is 12 arcsec north of where the 1953 plates were measured: a catalogue motion of 255 mas/yr
    pos_1953 = _ref(0.0, -12.0)[:2]
    rows = [{"survey": "poss1_red", "plate_epoch": 1953.0, "ra_measured_deg": pos_1953[0], "dec_measured_deg": pos_1953[1]},
            {"survey": "poss1_blue", "plate_epoch": 1953.0, "ra_measured_deg": pos_1953[0], "dec_measured_deg": pos_1953[1]},
            {"survey": "poss2ukstu_red", "plate_epoch": 1990.0, "ra_measured_deg": RA0, "dec_measured_deg": DEC0}]
    asked: dict[str, tuple] = {}

    def fake_gaia(ra, dec, radius_arcsec, cfg, cache, *, sleep=None):
        asked["gaia"] = (ra, dec, radius_arcsec)

        def star(sid, east, north, pmra, pmdec):
            r = _ref(east, north)
            return {"source_id": sid, "ra_deg": r[0], "dec_deg": r[1], "pmra_mas_yr": pmra, "pmdec_mas_yr": pmdec, "g_mag": 17.5}

        # "1" is 2.5 arcsec east of the place now and was 0.02 arcsec from it in 1953 (40 mas/yr east for 63 years);
        # "2" is at the place now but was 31 arcsec away in 1953; "3" does not move and is 30 arcsec off
        return [star("1", 2.5, -12.0, 40.0, 0.0), star("2", 0.0, -12.0, 500.0, 0.0), star("3", 30.0, 0.0, 0.0, 0.0)], "m", False

    def fake_ps1(ra, dec, radius_arcsec, min_detections, cache, **kw):
        asked["ps1"] = (ra, dec, radius_arcsec, min_detections)

        def obj(objid, east, north, r_mag):
            r = _ref(east, north)
            return {"objID": objid, "raMean": r[0], "decMean": r[1], "nDetections": 20, "rMeanPSFMag": r_mag}

        return [obj(21, 0.5, -12.5, 19.1), obj(22, 4.0, -12.0, -999.0), obj(23, 0.0, 0.0, 18.0)], "m", False

    monkeypatch.setattr(plate_survey, "gaia_cone", fake_gaia)
    monkeypatch.setattr(plate_survey, "ps1_cone", fake_ps1)
    cache = RawCache(tmp_path / "raw")

    def unit(rows=rows, dec=DEC0):
        return {"source_id": "USNO-B1.0 W", "plate_source": True, "usno": {"ra_deg": RA0, "dec_deg": dec},
                "vanishing_data": {"plate_images": {"rows": rows}}}

    out = pi.attach_first_place(unit(), cfg, cache)
    place = out["vanishing_data"]["plate_first_place"]
    assert out["vetting_archive"]["query_status"]["plate_first_place"] == "ok"
    assert place["position"] == [round(pos_1953[0], 7), round(pos_1953[1], 7)] and place["epoch"] == 1953.0
    # Gaia is read from the cone the cutouts were registered on, the Pan-STARRS objects around the 1953 position
    assert asked["gaia"] == (RA0, DEC0, 120.0) and asked["ps1"][:2] == pos_1953 and asked["ps1"][2:] == (pi.FIRST_PLACE_CONE_ARCSEC, 2)
    assert [(g["source_id"], g["sep_then_arcsec"]) for g in place["gaia"]] == [("1", 0.02)]
    assert [(p["objID"], p["sep_arcsec"], p["r_mag"]) for p in place["ps1"]] == [("21", 0.71, 19.1), ("22", 4.0, None)]
    # No first-epoch row, the southern sky, and an archive that is down
    none = pi.attach_first_place(unit(rows=rows[2:]), cfg, cache)
    assert none["vetting_archive"]["query_status"]["plate_first_place"] == "no_first_epoch_image" and "plate_first_place" not in none["vanishing_data"]
    south_rows = [dict(rows[0], dec_measured_deg=-45.0)]
    assert pi.attach_first_place(unit(rows=south_rows, dec=-45.0), cfg, cache)["vetting_archive"]["query_status"]["plate_first_place"] == "no_coverage"

    def down(*a, **kw):
        raise ArchiveAccessError("HTTP 503", retryable=True)

    monkeypatch.setattr(plate_survey, "ps1_cone", down)
    failed = pi.attach_first_place(unit(), cfg, cache)
    assert failed["vetting_archive"]["query_status"]["plate_first_place"] == "failed: HTTP 503" and "plate_first_place" not in failed["vanishing_data"]
    # A measured position far from the catalogue position gets its own Gaia cone
    far = _ref(0.0, -90.0)[:2]
    monkeypatch.setattr(plate_survey, "ps1_cone", fake_ps1)
    pi.attach_first_place(unit(rows=[dict(rows[0], ra_measured_deg=far[0], dec_measured_deg=far[1])]), cfg, cache)
    assert asked["gaia"] == (far[0], far[1], 120.0)
