"""SPHEREx forced photometry (offline; stamps are synthetic)."""

from __future__ import annotations

import numpy as np
import pytest

from oddstar.io import spherex as sx

OV = 5


def _epsf(fwhm_pix: float = 1.1, n: int = 33) -> np.ndarray:
    """An effective PSF: a Gaussian already integrated over a pixel, sampled 5 times per pixel, centre on the middle sample."""
    from scipy.special import erf

    sig = fwhm_pix / 2.3548
    off = (np.arange(n) - (n - 1) / 2.0) / OV          # pixel offsets of the samples
    one_d = 0.5 * (erf((off + 0.5) / (sig * np.sqrt(2.0))) - erf((off - 0.5) / (sig * np.sqrt(2.0))))
    return np.outer(one_d, one_d)


def test_rendered_psf_is_normalised_and_follows_the_pixel_phase() -> None:
    e = _epsf()
    centred = sx.render_epsf(e, (OV, OV), 0.0, 0.0)
    assert centred.shape == (9, 9) and centred.sum() == pytest.approx(1.0)
    assert np.unravel_index(np.argmax(centred), centred.shape) == (4, 4) and np.allclose(centred, centred.T, atol=1e-9)
    # A source 0.4 pixel to the right puts more light in the right-hand neighbour; the centroid follows it
    shifted = sx.render_epsf(e, (OV, OV), 0.4, -0.3)
    yy, xx = np.mgrid[-4:5, -4:5]
    assert shifted[4, 5] > shifted[4, 3] and shifted[3, 4] > shifted[5, 4]
    assert float((shifted * xx).sum()) == pytest.approx(0.4, abs=0.02) and float((shifted * yy).sum()) == pytest.approx(-0.3, abs=0.02)
    # The pixel fractions are those of the integrated Gaussian itself (cubic interpolation of the 5x samples)
    from scipy.special import erf
    sig = 1.1 / 2.3548
    def frac(c: float) -> float:
        return float(0.5 * (erf((c + 0.5) / (sig * np.sqrt(2))) - erf((c - 0.5) / (sig * np.sqrt(2)))))
    assert shifted[4, 4] == pytest.approx(frac(-0.4) * frac(0.3), rel=0.01)
    with pytest.raises(ValueError):
        sx.render_epsf(np.zeros((33, 33)), (OV, OV), 0.0, 0.0)


def test_fit_recovers_a_point_source_on_a_sloped_background() -> None:
    rng = np.random.default_rng(7)
    psf = sx.render_epsf(_epsf(), (OV, OV), 0.23, -0.41)
    yy, xx = np.mgrid[0:9, 0:9]
    sigma = 0.02
    truth = 3.7
    amps = []
    for _ in range(200):
        img = truth * psf + 0.50 + 0.01 * (xx - 4) - 0.02 * (yy - 4) + rng.normal(0.0, sigma, psf.shape)
        fit = sx.fit_point_source(img, np.full(psf.shape, sigma**2), np.ones(psf.shape, bool), psf)
        amps.append(fit["amplitude"])
    assert np.mean(amps) == pytest.approx(truth, abs=3 * fit["amplitude_err"] / np.sqrt(200))
    assert np.std(amps) == pytest.approx(fit["amplitude_err"], rel=0.2), "the formal error is the scatter of repeated fits"
    assert fit["background"] == pytest.approx(0.50, abs=0.02) and 0.5 < fit["reduced_chi2"] < 1.6 and fit["n_pixels"] == 81
    # Masked pixels are left out, and the flux is still recovered from the rest
    good = np.ones(psf.shape, bool)
    good[4, 5] = good[0, :] = False
    img = truth * psf + 0.5
    img[4, 5] = 99.0
    fit = sx.fit_point_source(img, np.full(psf.shape, sigma**2), good, psf)
    assert fit["amplitude"] == pytest.approx(truth, rel=1e-6) and fit["n_pixels"] == 81 - 10 and fit["psf_fraction_used"] < 1.0
    # No usable pixel, or shapes that differ: no number
    with pytest.raises(ValueError):
        sx.fit_point_source(img, np.full(psf.shape, sigma**2), np.zeros(psf.shape, bool), psf)
    with pytest.raises(ValueError):
        sx.fit_point_source(img[:5], np.full(psf.shape, sigma**2), good, psf)


def test_band_average_weights_by_error_and_reports_scatter() -> None:
    pts = [{"wavelength_um": 3.1, "flux_mjy": 10.0, "flux_err_mjy": 0.1, "mjd": 60800.0},
           {"wavelength_um": 3.5, "flux_mjy": 12.0, "flux_err_mjy": 0.2, "mjd": 60990.0},
           {"wavelength_um": 4.5, "flux_mjy": 30.0, "flux_err_mjy": 0.1, "mjd": 60801.0},
           {"wavelength_um": 3.3, "flux_mjy": 99.0, "flux_err_mjy": 0.0, "mjd": 60802.0}]   # no error: not used
    w1 = sx.band_average(pts, 2.8, 3.9)
    assert w1 is not None and w1["n"] == 2 and w1["flux_mjy"] == pytest.approx(10.4) and w1["flux_err_mjy"] == pytest.approx(0.0894, rel=1e-3)
    assert w1["scatter_mjy"] == pytest.approx(np.sqrt(2.0)) and (w1["mjd_min"], w1["mjd_max"]) == (60800.0, 60990.0)
    assert sx.band_average(pts, 1.0, 2.0) is None and sx.band_average(pts, 4.0, 5.0)["scatter_mjy"] is None


# ------------------------------------------------------------------ a synthetic Level 2 cutout, laid out like a real one

RA0, DEC0 = 45.165009, 0.200053
PIX_DEG = 6.15 / 3600.0
FLAG_CARDS = {0: "TRANSIENT", 1: "OVERFLOW", 2: "SUR_ERROR", 6: "NONFUNC", 10: "HOT", 19: "OUTLIER", 21: "SOURCE", 22: "GHOST"}


def _cutout(*, flux_mjy: float = 95.0, star_xy: tuple[float, float] = (8.37, 7.81), n: int = 17, zodi: float = 0.2,
            sky: float = 0.04, sigma: float = 0.004, neighbour: tuple[float, float, float] | None = None,
            bad_pixels: dict[tuple[int, int], int] | None = None, seed: int = 3) -> bytes:
    """IMAGE, FLAGS, VARIANCE, ZODI, EPSF (21 x 21 zones) and WCS-WAVE, with the star at 0-based pixel `star_xy`."""
    import io

    from astropy.io import fits

    rng = np.random.default_rng(seed)
    omega = (PIX_DEG * np.pi / 180.0) ** 2
    hdr = fits.Header()
    hdr["BUNIT"] = "MJy / sr"
    hdr["DETECTOR"] = 4
    hdr["MJD-OBS"], hdr["MJD-AVG"] = 60882.0511, 60882.0518
    # Sky: a plain tangent plane with the star's pixel as the reference (x to the west, y to the north)
    hdr["CTYPE1"], hdr["CTYPE2"] = "RA---TAN", "DEC--TAN"
    hdr["CRPIX1"], hdr["CRPIX2"] = star_xy[0] + 1.0, star_xy[1] + 1.0
    hdr["CRVAL1"], hdr["CRVAL2"] = RA0, DEC0
    hdr["CDELT1"], hdr["CDELT2"] = -PIX_DEG, PIX_DEG
    hdr["CUNIT1"], hdr["CUNIT2"] = "deg", "deg"
    # Detector: the cutout's first pixel is 0-based detector pixel (249, 1885)
    for ax, crpix in ((1, -248.0), (2, -1884.0)):
        hdr[f"CRPIX{ax}A"], hdr[f"CDELT{ax}A"], hdr[f"CRVAL{ax}A"] = crpix, 1.0, 0.0
        hdr[f"CTYPE{ax}A"], hdr[f"CUNIT{ax}A"] = "LINEAR", "pixel"
        hdr[f"CRPIX{ax}W"], hdr[f"CDELT{ax}W"], hdr[f"CRVAL{ax}W"] = crpix, 1.0, 1.0
        hdr[f"CTYPE{ax}W"], hdr[f"CUNIT{ax}W"] = "WAVE-TAB", "um"
        hdr[f"PS{ax}_0W"], hdr[f"PS{ax}_1W"] = "WCS-WAVE", "VALUES"
    hdr["WCSNAMEA"] = "0-based active pixel coordinates"
    hdr["PS1_2W"], hdr["PS2_2W"] = "X", "Y"
    hdr["PV1_3W"], hdr["PV2_3W"] = 1, 2
    # Wavelength rises with detector y from 2.40 to 3.85 micron (a linear variable filter); bandwidth = wavelength / 35
    xs, ys = np.linspace(1.0, 2040.0, 13), np.linspace(1.0, 2040.0, 20)
    values = np.zeros((20, 13, 2))
    values[:, :, 0] = (2.40 + (3.85 - 2.40) * (ys - 1.0) / 2039.0)[:, None]
    values[:, :, 1] = values[:, :, 0] / 35.0
    wave_tab = fits.BinTableHDU.from_columns([
        fits.Column(name="X", format="13E", array=xs[None, :]), fits.Column(name="Y", format="20E", array=ys[None, :]),
        fits.Column(name="VALUES", format="520E", dim="(2,13,20)", array=values[None, ...])], name="WCS-WAVE")
    # EPSF zones: narrower towards low x, so that a wrong zone gives a wrong flux
    centres = (np.arange(21) + 0.5) * 2040.0 / 21.0
    rows_x, rows_y = np.meshgrid(centres, centres)
    epsfs = np.stack([_epsf(0.9 + 0.5 * cx / 2040.0) for cx in rows_x.ravel()])
    epsf_tab = fits.BinTableHDU.from_columns([
        fits.Column(name="XCENTER", format="E", array=rows_x.ravel()), fits.Column(name="YCENTER", format="E", array=rows_y.ravel()),
        fits.Column(name="CWAVE", format="E", array=2.40 + 1.45 * rows_y.ravel() / 2040.0),
        fits.Column(name="EPSF", format="1089D", dim="(33,33)", array=epsfs)], name="EPSF")
    epsf_tab.header["OVSMPX"], epsf_tab.header["OVSMPY"] = OV, OV

    def render(x: float, y: float, f_mjy: float) -> np.ndarray:
        ix, iy = int(round(x)), int(round(y))
        zone = int(np.argmin((rows_x.ravel() - (x + 249.0)) ** 2 + (rows_y.ravel() - (y + 1885.0)) ** 2))
        big = sx.render_epsf(epsfs[zone], (OV, OV), x - ix, y - iy, half=6)
        out = np.zeros((n, n))
        ys_, xs_ = np.mgrid[iy - 6:iy + 7, ix - 6:ix + 7]
        ok = (ys_ >= 0) & (ys_ < n) & (xs_ >= 0) & (xs_ < n)
        out[ys_[ok], xs_[ok]] = big[ok] * f_mjy * 1e-9 / omega
        return out

    yy, xx = np.mgrid[0:n, 0:n]
    image = zodi + sky + 0.002 * (xx - 8) + render(star_xy[0], star_xy[1], flux_mjy) + rng.normal(0.0, sigma, (n, n))
    flags = np.zeros((n, n), dtype=np.int32)
    flags[np.hypot(xx - star_xy[0], yy - star_xy[1]) < 2.0] = 1 << 21
    if neighbour is not None:
        image += render(*neighbour)
    for (px, py), value in (bad_pixels or {}).items():
        flags[py, px] |= value
        image[py, px] += 50.0
    flag_hdu = fits.ImageHDU(flags, name="FLAGS")
    for i, (bit, name) in enumerate(FLAG_CARDS.items()):
        flag_hdu.header[f"MSKN{i:04d}"], flag_hdu.header[f"MSKM{i:04d}"] = name, 1 << bit
    primary = fits.PrimaryHDU()
    primary.header["VERSION"] = "7.0.5"
    hdul = fits.HDUList([primary, fits.ImageHDU(image.astype(np.float32), header=hdr, name="IMAGE"), flag_hdu,
                         fits.ImageHDU(np.full((n, n), sigma**2, dtype=np.float32), name="VARIANCE"),
                         fits.ImageHDU(np.full((n, n), zodi, dtype=np.float32), name="ZODI"), epsf_tab, wave_tab])
    buf = io.BytesIO()
    hdul.writeto(buf)
    return buf.getvalue()


def test_cutout_gives_flux_wavelength_and_detector_position() -> None:
    row = sx.measure_cutout(_cutout(), RA0, DEC0)
    assert row["status"] == "ok" and row["detector"] == 4 and row["mjd"] == pytest.approx(60882.0518) and row["pipeline_version"] == "7.0.5"
    assert (row["x_cutout"], row["y_cutout"]) == (pytest.approx(8.37, abs=2e-3), pytest.approx(7.81, abs=2e-3))
    assert (row["x_detector"], row["y_detector"]) == (pytest.approx(257.37, abs=0.01), pytest.approx(1892.81, abs=0.01))
    # Wavelength of the star's own detector row, read through the file's lookup table
    assert row["wavelength_um"] == pytest.approx(2.40 + 1.45 * (1893.81 - 1.0) / 2039.0, abs=2e-4)
    assert row["bandwidth_um"] == pytest.approx(row["wavelength_um"] / 35.0, rel=1e-3)
    assert row["pixel_solid_angle_sr"] == pytest.approx((PIX_DEG * np.pi / 180.0) ** 2, rel=1e-4)
    assert row["flux_mjy"] == pytest.approx(95.0, abs=4 * row["flux_err_mjy"]) and row["flux_err_mjy"] < 0.2
    assert row["flags_core"] == 1 << 21 and row["n_pixels"] == 81 and 0.5 < row["reduced_chi2"] < 1.6
    assert row["background_mjy_sr"] == pytest.approx(0.04 + 0.002 * (8 - 8), abs=0.01)   # the zodiacal model is taken off first


def test_flagged_pixels_are_left_out_and_a_flagged_core_gives_no_number() -> None:
    bits = {name: 1 << bit for bit, name in FLAG_CARDS.items()}
    # A cosmic-ray hit in the wing: masked, the flux is unchanged. A ghost bit is reported, not masked.
    row = sx.measure_cutout(_cutout(bad_pixels={(11, 10): bits["TRANSIENT"] | bits["SUR_ERROR"]}), RA0, DEC0)
    assert row["status"] == "ok" and row["n_pixels"] == 80 and row["flux_mjy"] == pytest.approx(95.0, abs=4 * row["flux_err_mjy"])
    # The same hit on the star's own pixels: no measurement, and the flags say why
    row = sx.measure_cutout(_cutout(bad_pixels={(8, 8): bits["TRANSIENT"]}), RA0, DEC0)
    assert row["status"] == "core_flagged" and row["flags_core"] & bits["TRANSIENT"] and "flux_mjy" not in row
    # A star within four pixels of the cutout's edge (the detector's edge): no measurement
    assert sx.measure_cutout(_cutout(star_xy=(2.2, 7.8)), RA0, DEC0)["status"] == "off_cutout"
    # A file that is not a Level 2 cutout is an error, not a zero
    from astropy.io import fits
    import io
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.zeros((5, 5)), name="IMAGE")]).writeto(buf)
    with pytest.raises(ValueError, match="not a SPHEREx Level 2 cutout"):
        sx.measure_cutout(buf.getvalue(), RA0, DEC0)


def test_a_known_neighbour_in_the_stamp_is_fitted_with_the_star() -> None:
    # A neighbour 3 pixels (18 arcsec) east with a third of the star's flux
    data = _cutout(neighbour=(8.37 - 3.0, 7.81, 30.0))
    nb = (RA0 + 3.0 * PIX_DEG / np.cos(np.radians(DEC0)), DEC0)
    alone = sx.measure_cutout(data, RA0, DEC0)
    both = sx.measure_cutout(data, RA0, DEC0, neighbours=(nb,))
    assert both["n_neighbours_fitted"] == 1 and both["flux_mjy"] == pytest.approx(95.0, abs=4 * both["flux_err_mjy"])
    assert both["reduced_chi2"] < 1.6 < alone["reduced_chi2"], "ignoring the neighbour leaves a misfit that the fit quality shows"
    # A neighbour closer than a pixel cannot be separated: it is counted as blended, not fitted
    close = sx.measure_cutout(data, RA0, DEC0, neighbours=((RA0 + 0.4 * PIX_DEG, DEC0),))
    assert close["n_blended"] == 1 and close["n_neighbours_fitted"] == 0


# ------------------------------------------------------------------ archive access (no network: the GET is replaced)


class _Resp:
    def __init__(self, content: bytes = b"", payload: object = None) -> None:
        self.content, self._payload, self.status_code = content, payload, 200

    def json(self) -> object:
        return self._payload


def _tap_answer(rows: list[tuple]) -> dict:
    """IRSA's layout: a VOTable written out as JSON, every value a string."""
    names = [("uri", "char"), ("planeid", "char"), ("energy_bandpassname", "char"), ("time_bounds_lower", "double"),
             ("time_bounds_upper", "double"), ("collection", "char"), ("observationid", "char")]
    return {"VOTABLE": {"RESOURCE_ARRAY": [{"INFO_ARRAY": [{"<xmlattr>": {"name": "QUERY_STATUS", "value": "OK"}}], "TABLE": {
        "FIELD_ARRAY": [{"<xmlattr>": {"name": n, "datatype": t}} for n, t in names],
        "DATA": {"TABLEDATA": [[str(v) for v in r] for r in rows]}}}]}}


def _archive(monkeypatch, tmp_path, rows, cutouts, calls, fail=()):
    from oddstar.errors import ArchiveAccessError
    from oddstar.io.cache import RawCache

    def fake_get(url, params, **kwargs):
        calls.append(url)
        if url == sx.IRSA_TAP_URL:
            assert "spherex.plane" in params["QUERY"] and "CONTAINS(POINT('ICRS', 45.165009, 0.200" in params["QUERY"]
            return _Resp(payload=_tap_answer(rows))
        name = url.split("?")[0].rsplit("/", 1)[1]
        if name in fail:
            raise ArchiveAccessError(f"GET {url} failed after 4 attempt(s): HTTP 503", retryable=True)
        assert url.startswith("https://irsa.ipac.caltech.edu/ibe/data/spherex/") and "?center=45.165" in url and url.endswith("d&size=0.03")
        return _Resp(content=cutouts[name])

    monkeypatch.setattr(sx.ot, "_get_with_retries", fake_get)
    return RawCache(tmp_path / "raw")


def _row(name: str, band: int, mjd: float, collection: str = "spherex_dr1") -> tuple:
    return (f"ibe/data/spherex/dr1/level2/2025W30_2B/l2b/{band}/{name}", f"plane-{name}", f"SPHEREx-D{band}", mjd, mjd + 0.0013,
            collection, name[:-5])


def test_every_exposure_is_measured_once_and_kept(monkeypatch, tmp_path) -> None:
    rows = [_row("b.fits", 4, 60890.0), _row("a.fits", 4, 60882.0), _row("c.fits", 6, 61080.0, "spherex_qr3"),
            ("ibe/data/spherex/dr1/level3/cube.fits", "p", "SPHEREx-D4", 60882.0, 60883.0, "spherex_dr1_cube", "cube")]
    cutouts = {"a.fits": _cutout(flux_mjy=95.0, seed=1), "b.fits": _cutout(flux_mjy=95.0, seed=2, star_xy=(2.0, 8.0)),
               "c.fits": _cutout(flux_mjy=60.0, seed=3)}
    calls: list[str] = []
    cache = _archive(monkeypatch, tmp_path, rows, cutouts, calls)
    res = sx.forced_photometry("9281425163264", RA0, DEC0, cache, as_of="2026-11-18", workers=1)
    assert res["n_images"] == 3 and res["n_measured"] == 2 and res["complete"] and res["status_counts"] == {"ok": 2, "off_cutout": 1}
    assert [p["uri"].rsplit("/", 1)[1] for p in res["points"]] == ["a.fits", "c.fits"] and res["points"][1]["collection"] == "spherex_qr3"
    assert res["not_measured"][0]["status"] == "off_cutout" and res["collections"] == ["spherex_dr1", "spherex_qr3"]
    assert len(calls) == 4   # one image list, three cutouts
    # A second call reads the stored rows: nothing is fetched again, and the answer is the same
    again = sx.forced_photometry("9281425163264", RA0, DEC0, cache, as_of="2026-11-18", workers=1)
    assert len(calls) == 4 and again["points"] == res["points"]
    # A new `as_of` asks the archive for the list again, and only for the list
    sx.forced_photometry("9281425163264", RA0, DEC0, cache, as_of="2026-12-01", workers=1)
    assert len(calls) == 5 and calls[-1] == sx.IRSA_TAP_URL
    # `detectors` limits the bands
    only6 = sx.forced_photometry("9281425163264", RA0, DEC0, cache, as_of="2026-11-18", detectors=(6,), workers=1)
    assert only6["n_images"] == 1 and only6["detectors"] == ["SPHEREx-D6"]


def test_a_cutout_that_cannot_be_fetched_is_counted_and_fetched_on_the_next_call(monkeypatch, tmp_path) -> None:
    rows = [_row("a.fits", 4, 60882.0), _row("b.fits", 4, 60890.0)]
    cutouts = {"a.fits": _cutout(seed=1), "b.fits": _cutout(seed=2)}
    calls: list[str] = []
    cache = _archive(monkeypatch, tmp_path, rows, cutouts, calls, fail=("b.fits",))
    res = sx.forced_photometry("1", RA0, DEC0, cache, as_of="x", workers=2)
    assert res["n_failed"] == 1 and not res["complete"] and res["n_measured"] == 1 and "HTTP 503" in res["errors"][0]
    cache = _archive(monkeypatch, tmp_path, rows, cutouts, calls)
    n_before = len(calls)
    res = sx.forced_photometry("1", RA0, DEC0, cache, as_of="x", workers=2)
    assert res["complete"] and res["n_measured"] == 2 and len(calls) == n_before + 1, "only the missing cutout is fetched"
    # A file that is not a cutout (an error page) is a failure of that image, never a measurement
    cache = _archive(monkeypatch, tmp_path, [_row("d.fits", 5, 60900.0)], {"d.fits": b"<html>Server busy</html>"}, calls)
    res = sx.forced_photometry("1", RA0, DEC0, cache, as_of="y", workers=1)
    assert res["n_failed"] == 1 and "not readable" in res["errors"][0]


def test_the_star_is_measured_where_its_proper_motion_puts_it(monkeypatch, tmp_path) -> None:
    # 2 arcsec per year to the north: 9.56 years after 2016.0 the star is 19.1 arcsec (3.1 pixels) north of its Gaia position
    mjd = 60882.0
    years = 2000.0 + (mjd - 51544.5) / 365.25 - 2016.0
    dec_then = DEC0 - 2.0 * years / 3600.0
    calls: list[str] = []
    cache = _archive(monkeypatch, tmp_path, [_row("a.fits", 4, mjd)], {"a.fits": _cutout(seed=5)}, calls)
    moved = sx.forced_photometry("1", RA0, dec_then, cache, as_of="x", pmdec_mas_yr=2000.0, workers=1)
    assert moved["points"][0]["dec_deg"] == pytest.approx(DEC0, abs=1e-6)
    assert moved["points"][0]["flux_mjy"] == pytest.approx(95.0, abs=4 * moved["points"][0]["flux_err_mjy"])


def test_summary_gives_window_means_per_visit() -> None:
    def pt(wave: float, flux: float, mjd: float) -> dict:
        return {"wavelength_um": wave, "flux_mjy": flux, "flux_err_mjy": 0.01 * flux, "mjd": mjd}
    k0, w1_0, w2_0 = 666.7e3 * 10 ** (-0.4 * 10.0), 309.54e3 * 10 ** (-0.4 * 10.0), 171.787e3 * 10 ** (-0.4 * 9.5)
    pts = [pt(2.1, k0, 60880.0), pt(2.3, k0, 60885.0), pt(3.2, w1_0, 60890.0), pt(4.5, w2_0, 60895.0),
           pt(3.4, w1_0, 61080.0), pt(4.6, w2_0 * 10 ** (-0.4 * 0.3), 61085.0), pt(1.2, 5.0, 61081.0)]
    s = sx.summarize(pts)
    assert s["all"]["n"] == 7 and len(s["visits"]) == 2 and [v["n"] for v in s["visits"]] == [4, 3]
    first, second = s["visits"]
    assert first["k"]["n"] == 2 and first["k"]["mag"] == pytest.approx(10.0) and first["w1"]["mag"] == pytest.approx(10.0)
    assert first["w1_w2"] == pytest.approx(0.5) and first["k_w1"] == pytest.approx(0.0) and first["k_w2"] == pytest.approx(0.5)
    assert second["w1_w2"] == pytest.approx(0.2) and second["k"] is None and second["k_w1"] is None
    assert first["w2"]["mag_err"] == pytest.approx(2.5 / np.log(10) * 0.01)
    assert sx.summarize([]) == {"all": {"n": 0, "k": None, "w1": None, "w2": None, "k_w1": None, "w1_w2": None, "k_w2": None}, "visits": []}


def test_command_measures_chosen_stars_and_writes_one_file_per_star(monkeypatch, tmp_path) -> None:
    """`oddstar spherex photometry`: stars at Tier C or better by default, `--ids` for any star of the sample."""
    import json
    import shutil

    from oddstar import cli
    from oddstar.cli import dispatch
    from oddstar.io.cache import RawCache
    from oddstar.config import repo_root

    rows = [_row("a.fits", 4, 60882.0), _row("c.fits", 6, 61080.0, "spherex_qr3")]
    calls: list[str] = []
    _archive(monkeypatch, tmp_path, rows, {"a.fits": _cutout(seed=1), "c.fits": _cutout(flux_mjy=60.0, seed=3)}, calls)
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / "raw"))
    asked: list[str] = []

    def fake_neighbours(sid, ra, dec, g_mag, cache, *, release, **kw):
        asked.append(f"{sid}:{g_mag}")
        return ()

    monkeypatch.setattr(sx, "gaia_neighbours", fake_neighbours)
    stars = [{"source_id": "1", "ra_deg": RA0, "dec_deg": DEC0, "pmra_mas_yr": 0.0, "pmdec_mas_yr": 0.0, "mags": {"G": 10.6}},
             {"source_id": "2", "ra_deg": RA0, "dec_deg": DEC0, "mags": {}}]
    runs = repo_root() / "runs"
    src, out_dir = runs / "test-spherex-src", runs / "spherex-test-spherex-src"
    try:
        src.mkdir(parents=True)
        score = {"source_id": "2", "channel": "ir_excess", "raw_statistic": 9.0, "false_alarm_prob": 1e-4, "flagged": True, "rank_score": 9.0,
                 "summary": "x"}
        (src / "search_ir_excess.json").write_text(json.dumps({
            "channel": "ir_excess", "sample_stars": stars, "flagged_ids": ["2"], "scores": [score]}), encoding="utf-8")
        code, payload = dispatch(["spherex", "photometry", "--from-run", "test-spherex-src", "--dry-run", "--json"])
        assert code == 0 and payload["data"]["source_ids"] == ["2"] and payload["data"]["selection"] == "tier C or better", payload
        assert not out_dir.exists() and calls == [], "a dry run fetches and writes nothing"
        code, payload = dispatch(["spherex", "photometry", "--from-run", "test-spherex-src", "--ids", "1", "--as-of", "2026-11-18",
                                  "--workers", "1", "--json"])
        assert code == 0, payload
        d = payload["data"]
        assert d["status_counts"] == {"measured": 1} and d["photometry_version"] == sx.PHOTOMETRY_VERSION and asked == ["1:10.6"]
        assert d["stars"][0]["n_measured"] == 2 and d["stars"][0]["summary"]["all"]["n"] == 2 and "points" not in d["stars"][0]
        saved = json.loads((out_dir / "spherex" / "1.json").read_text(encoding="utf-8"))
        assert len(saved["points"]) == 2 and saved["points"][0]["flux_mjy"] == pytest.approx(95.0, rel=0.02) and saved["as_of"] == "2026-11-18"
        assert json.loads((out_dir / "spherex_photometry.json").read_text(encoding="utf-8"))["n_selected"] == 1
        code, payload = dispatch(["spherex", "photometry", "--from-run", "test-spherex-src", "--ids", "1", "--detectors", "7", "--json"])
        assert code != 0 and "1 to 6" in json.dumps(payload)
        code, payload = dispatch(["spherex", "photometry", "--from-run", "test-spherex-src", "--ids", "99", "--json"])
        assert code != 0 and "No position" in json.dumps(payload)
    finally:
        shutil.rmtree(src, ignore_errors=True)
        shutil.rmtree(out_dir, ignore_errors=True)


def test_neighbour_list_keeps_stars_bright_enough_to_matter(monkeypatch, tmp_path) -> None:
    from oddstar.io import live_survey
    from oddstar.io.cache import RawCache

    rows = [{"source_id": 1, "ra": RA0, "dec": DEC0, "pmra": 0.0, "pmdec": 0.0, "phot_g_mean_mag": 10.6},      # the star itself
            {"source_id": 2, "ra": RA0 + 0.004, "dec": DEC0, "pmra": None, "pmdec": None, "phot_g_mean_mag": 13.0},
            {"source_id": 3, "ra": RA0, "dec": DEC0 + 0.004, "pmra": 0.0, "pmdec": 3600.0, "phot_g_mean_mag": 15.5},
            {"source_id": 4, "ra": RA0 - 0.004, "dec": DEC0, "pmra": 0.0, "pmdec": 0.0, "phot_g_mean_mag": 19.0},  # more than 5 mag fainter
            {"source_id": 5, "ra": RA0, "dec": DEC0 - 0.004, "pmra": 0.0, "pmdec": 0.0, "phot_g_mean_mag": None}]
    seen: list[str] = []

    def fake_tap(*, urls, query, **kw):
        seen.append(query)
        return rows, None, False

    monkeypatch.setattr(live_survey, "_tap_query_cached", fake_tap)
    cache = RawCache(tmp_path / "raw")
    out = sx.gaia_neighbours("1", RA0, DEC0, 10.6, cache, release="dr3")
    assert "CIRCLE('ICRS', 45.165009, 0.200053, 0.008333)" in seen[0] and len(out) == 2
    assert (RA0 + 0.004, DEC0) in out and any(d == pytest.approx(DEC0 + 0.004 + 0.01) for _, d in out), "moved by its proper motion to 2026"
    assert len(sx.gaia_neighbours("1", RA0, DEC0, None, cache, release="dr3")) == 4, "with no magnitude for the star, every source counts"
