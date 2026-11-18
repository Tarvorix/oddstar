"""SPHEREx spectral images (IRSA): forced point-source photometry at a star's position.

SPHEREx maps the whole sky at 0.75 to 5.0 micron through linear variable filters: every exposure images each sky
position at one wavelength, set by where the position falls on the detector, and the survey returns to a
position about every six months with all wavelengths. For a star this gives a 0.75 to 5 micron spectral energy
distribution from an instrument that is not WISE, and epochs from 2025 on, after NEOWISE ended in 2024.

What is measured, per exposure that contains the position:

1. The image list comes from IRSA's TAP service (`spherex.artifact` joined to `spherex.plane`, footprint `poly`).
2. A cutout of each Level 2 spectral image is fetched (`<uri>?center=<ra>,<dec>d&size=<deg>`). It holds the
   calibrated surface brightness (IMAGE, MJy/sr, zodiacal light included), FLAGS, VARIANCE, the zodiacal model
   (ZODI), the effective point spread function in 21 x 21 detector zones (EPSF, 5x oversampled) and the lookup
   table of the wavelength of each pixel (alternate WCS "W").
3. The zodiacal model is subtracted, and a point source at the star's fixed position plus a sloped plane
   background is fitted to a small stamp, weighted by the variance, with flagged pixels left out. Only the
   amplitude and the plane are free (and one amplitude for each known neighbour in the stamp, at its fixed
   position). The amplitude times the solid angle of a pixel is the flux density.

The star's position is its Gaia position moved by its proper motion to the date of the exposure.

The cutouts are not kept (about 4 MB each whatever their size, because each carries the whole PSF table; a star
has a few hundred): the measured row of each exposure is stored in the raw cache under `PHOTOMETRY_VERSION`.

A SPHEREx pixel is 6.15 arcsec wide, as coarse as WISE: this photometry confirms or contradicts the brightness of
the blended light at the position. It cannot separate a neighbour that WISE could not separate.
"""

from __future__ import annotations

import hashlib
import io
import json
import math
import time
from concurrent.futures import ThreadPoolExecutor
from typing import Any, Callable

import numpy as np

from oddstar.errors import ArchiveAccessError
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.io.physics import propagate_position

ARCHIVE = "spherex"
PHOTOMETRY_VERSION = "spherex_forced_v1"
IRSA_BASE_URL = "https://irsa.ipac.caltech.edu/"
IRSA_TAP_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
# Side of the cutout asked for (degrees): about 18 pixels of 6.15 arcsec, room for the fit stamp at any pixel phase
CUTOUT_SIZE_DEG = 0.03
# The fit stamp is (2 * STAMP_HALF + 1) pixels on a side; the effective PSF spans 33 / 5 = 6.6 pixels
STAMP_HALF = 4
# Fewer good pixels than this in the stamp, or a flagged pixel within this radius of the star: no measurement
MIN_GOOD_PIXELS = 30
CORE_RADIUS_PIX = 1.5
# A neighbour closer to the star than this cannot be told from it at 6.15 arcsec per pixel: its light stays in the
# star's measurement and the row counts it under `n_blended`. Farther neighbours inside the stamp get their own amplitude.
BLEND_RADIUS_PIX = 1.0
# Flag bits (by the names the FLAGS header gives them) whose pixels are left out of the fit. SOURCE marks pixels
# that hold a catalogued source, the star itself included, and is not a defect. Ghost, halo, crosstalk and dichroic
# bits mark possible contamination of pixels that still hold a valid measurement: they are reported on the row
# (`flags_core`), not masked.
BAD_FLAG_NAMES = ("TRANSIENT", "OVERFLOW", "SUR_ERROR", "NONFUNC", "MISSING_DATA", "HOT", "COLD", "PHANMISS",
                  "NONLINEAR", "PERSIST", "OUTLIER", "STREAK", "BLOOM", "SNOWBALL")
# Exposures closer in time than this belong to one visit (SPHEREx returns to a position about every six months
# and takes a few weeks to sweep all wavelengths across it)
VISIT_GAP_DAYS = 60.0
# Wavelength windows of the summary (micron): the K band, and the parts of the WISE W1 and W2 bands SPHEREx covers
# (it ends at 5.0 micron; W2 runs to 5.3). Extraction geometry, not a threshold.
SED_WINDOWS_UM = {"k": (2.0, 2.4), "w1": (3.0, 3.8), "w2": (4.1, 5.0)}
# Vega zero points (Jy) used to quote a window mean as a magnitude: 2MASS Ks, WISE W1, WISE W2. A window mean is not
# an integral through the 2MASS or WISE response, so these magnitudes differ from catalogue ones by a colour term
# of a few hundredths of a magnitude; compare stars with stars measured the same way.
WINDOW_ZERO_POINT_JY = {"k": 666.7, "w1": 309.54, "w2": 171.787}

# Date at which a moving star's position is tested against the image footprints
LIST_EPOCH = 2026.0

Sleeper = Callable[[float], None]


def render_epsf(epsf: np.ndarray, oversampling: tuple[int, int], dx: float, dy: float, half: int = STAMP_HALF) -> np.ndarray:
    """The fraction of a point source's flux in each detector pixel of a (2 half + 1)^2 stamp.

    `epsf` is the effective PSF sampled `oversampling` (y, x) times per pixel with its centre on the middle
    sample; the source sits at (`dx`, `dy`) pixels from the centre of the stamp's middle pixel. An effective PSF
    is already integrated over a pixel, so a pixel's value is the PSF read at that pixel's offset from the source.
    The stamp is normalised to a total of one.
    """
    from scipy.ndimage import map_coordinates

    arr = np.asarray(epsf, dtype=float)
    if arr.ndim != 2 or not np.all(np.isfinite(arr)) or arr.sum() <= 0.0:
        raise ValueError("effective PSF is not a finite 2-D array with positive sum")
    ov_y, ov_x = int(oversampling[0]), int(oversampling[1])
    cy, cx = (arr.shape[0] - 1) / 2.0, (arr.shape[1] - 1) / 2.0
    offs = np.arange(-half, half + 1, dtype=float)
    yy, xx = np.meshgrid(offs - dy, offs - dx, indexing="ij")
    stamp = map_coordinates(arr, [yy * ov_y + cy, xx * ov_x + cx], order=3, mode="constant", cval=0.0)
    stamp = np.clip(stamp, 0.0, None)
    total = float(stamp.sum())
    if total <= 0.0:
        raise ValueError("effective PSF renders to nothing on the stamp")
    return stamp / total


def fit_point_source(image: np.ndarray, variance: np.ndarray, good: np.ndarray, psf: np.ndarray,
                     others: tuple[np.ndarray, ...] | list[np.ndarray] = ()) -> dict[str, float]:
    """Weighted linear fit of `amplitude * psf + plane` to a stamp. Returns amplitude, its error and the fit quality.

    The plane (constant, x slope, y slope) takes up the sky left after the zodiacal model: across a few pixels the
    wavelength changes, and with it the background. `good` marks the pixels used. `others` are the rendered PSFs of
    neighbours at their fixed positions; each gets a free amplitude of its own. Raises ValueError when the fit is
    not determined.
    """
    img, var, ok, model = (np.asarray(a) for a in (image, variance, good, psf))
    extra = [np.asarray(o, dtype=float) for o in others]
    if not (img.shape == var.shape == ok.shape == model.shape) or img.ndim != 2 or any(o.shape != img.shape for o in extra):
        raise ValueError("stamp arrays differ in shape")
    half_y, half_x = (img.shape[0] - 1) / 2.0, (img.shape[1] - 1) / 2.0
    yy, xx = np.mgrid[0:img.shape[0], 0:img.shape[1]]
    use = ok.astype(bool) & np.isfinite(img) & np.isfinite(var) & (var > 0.0)
    n = int(use.sum())
    n_par = 4 + len(extra)
    if n <= n_par:
        raise ValueError("no more usable pixels than free parameters")
    design = np.stack([model[use], np.ones(n), (xx[use] - half_x), (yy[use] - half_y)] + [o[use] for o in extra], axis=1)
    w = 1.0 / var[use]
    lhs = design.T @ (design * w[:, None])
    try:
        cov = np.linalg.inv(lhs)
    except np.linalg.LinAlgError as exc:
        raise ValueError("fit is singular") from exc
    beta = cov @ (design.T @ (img[use] * w))
    resid = img[use] - design @ beta
    chi2 = float(np.sum(resid * resid * w))
    if not np.all(np.isfinite(beta)) or cov[0, 0] <= 0.0:
        raise ValueError("fit is not finite")
    return {"amplitude": float(beta[0]), "amplitude_err": float(math.sqrt(cov[0, 0])), "background": float(beta[1]),
            "chi2": chi2, "n_pixels": n, "reduced_chi2": chi2 / (n - n_par), "psf_fraction_used": float(model[use].sum())}


def band_average(points: list[dict[str, Any]], lo_um: float, hi_um: float) -> dict[str, Any] | None:
    """Weighted mean flux density of the measurements whose central wavelength lies in [lo, hi) micron.

    The scatter of the points about the mean is reported beside the formal error: a stellar spectrum is not flat
    across a WISE band, so the scatter holds the spectral slope as well as the noise.
    """
    sel = [p for p in points if lo_um <= float(p["wavelength_um"]) < hi_um and float(p["flux_err_mjy"]) > 0.0]
    if not sel:
        return None
    f = np.array([float(p["flux_mjy"]) for p in sel])
    w = 1.0 / np.array([float(p["flux_err_mjy"]) for p in sel]) ** 2
    mean = float(np.sum(f * w) / np.sum(w))
    return {"n": len(sel), "flux_mjy": mean, "flux_err_mjy": float(1.0 / math.sqrt(np.sum(w))),
            "scatter_mjy": float(np.std(f, ddof=1)) if len(sel) > 1 else None, "median_mjy": float(np.median(f)),
            "mjd_min": float(min(p["mjd"] for p in sel)), "mjd_max": float(max(p["mjd"] for p in sel))}


# ------------------------------------------------------------------------------------------ one cutout


def flag_bits(header: Any) -> dict[str, int]:
    """Name -> mask value of the bit planes of a FLAGS extension (cards MSKNnnnn name, MSKMnnnn value)."""
    out: dict[str, int] = {}
    for key in header:
        if not str(key).startswith("MSKN"):
            continue
        suffix = str(key)[4:]
        name = str(header[key]).strip().upper()
        if f"MSKM{suffix}" in header:
            out[name] = int(header[f"MSKM{suffix}"])
        elif suffix.isdigit():
            out[name] = 1 << int(suffix)
    return out


def _pixel_solid_angle_sr(wcs: Any, x: float, y: float) -> float:
    """Solid angle of the detector pixel at (x, y), from the image's own sky solution with its distortion."""
    ra, dec = wcs.all_pix2world([x - 0.5, x + 0.5, x, x], [y, y, y - 0.5, y + 0.5], 0)
    cosd = math.cos(math.radians(float(np.mean(dec))))
    ax = ((ra[1] - ra[0] + 180.0) % 360.0 - 180.0) * cosd
    ay = dec[1] - dec[0]
    bx = ((ra[3] - ra[2] + 180.0) % 360.0 - 180.0) * cosd
    by = dec[3] - dec[2]
    return abs(ax * by - ay * bx) * (math.pi / 180.0) ** 2


def measure_cutout(fits_bytes: bytes, ra: float, dec: float, neighbours: tuple[tuple[float, float], ...] = ()) -> dict[str, Any]:
    """Forced photometry of the position in one Level 2 cutout. Returns the measured row (`status` "ok" or the reason
    there is no number: "off_cutout", "core_flagged", "too_few_pixels", "fit_undetermined").

    `neighbours` are (ra, dec) of other catalogued stars; those that fall in the stamp are fitted with it.
    Raises ValueError when the file does not have the layout of a Level 2 cutout.
    """
    import warnings

    from astropy import units as u
    from astropy.io import fits
    from astropy.wcs import WCS, FITSFixedWarning

    with warnings.catch_warnings():
        # The header carries the distortion (SIP) of the sky solution beside the wavelength solution "W", which has
        # no distortion of its own; astropy warns about that on every file
        warnings.simplefilter("ignore", FITSFixedWarning)
        warnings.filterwarnings("ignore", message=".*SIP.*")
        with fits.open(io.BytesIO(fits_bytes), memmap=False) as hdul:
            names = {h.name for h in hdul}
            missing = {"IMAGE", "FLAGS", "VARIANCE", "ZODI", "EPSF", "WCS-WAVE"} - names
            if missing:
                raise ValueError(f"not a SPHEREx Level 2 cutout: no {sorted(missing)} extension (has {sorted(names)})")
            hdr = hdul["IMAGE"].header
            if str(hdr.get("BUNIT", "")).replace(" ", "") != "MJy/sr":
                raise ValueError(f"IMAGE unit is {hdr.get('BUNIT')!r}, not MJy/sr")
            image = np.asarray(hdul["IMAGE"].data, dtype=float)
            flags = np.asarray(hdul["FLAGS"].data).astype(np.int64)
            variance = np.asarray(hdul["VARIANCE"].data, dtype=float)
            zodi = np.asarray(hdul["ZODI"].data, dtype=float)
            sky = WCS(hdr)
            x, y = (float(v) for v in sky.all_world2pix(ra, dec, 0))
            row: dict[str, Any] = {
                "version": PHOTOMETRY_VERSION, "detector": int(hdr["DETECTOR"]), "mjd": float(hdr.get("MJD-AVG", hdr.get("MJD-OBS"))),
                "pipeline_version": str(hdul[0].header.get("VERSION", "")), "x_cutout": round(x, 3), "y_cutout": round(y, 3),
            }
            ix, iy, h = int(round(x)), int(round(y)), STAMP_HALF
            if not (np.isfinite(x) and np.isfinite(y)) or ix - h < 0 or iy - h < 0 or ix + h >= image.shape[1] or iy + h >= image.shape[0]:
                return {**row, "status": "off_cutout"}   # the star lies within a few pixels of the detector's edge
            # Position on the full detector ("0-based active pixel coordinates", alternate system A, linear)
            xp = float(hdr.get("CRVAL1A", 0.0)) + float(hdr.get("CDELT1A", 1.0)) * (x + 1.0 - float(hdr["CRPIX1A"]))
            yp = float(hdr.get("CRVAL2A", 0.0)) + float(hdr.get("CDELT2A", 1.0)) * (y + 1.0 - float(hdr["CRPIX2A"]))
            wave_wcs = WCS(hdr, fobj=hdul, key="W")
            wave_wcs.sip = None
            wave, width = wave_wcs.pixel_to_world(x, y)
            wavelength_um, bandwidth_um = float(u.Quantity(wave).to_value(u.um)), float(u.Quantity(width).to_value(u.um))
            table = hdul["EPSF"].data
            ehdr = hdul["EPSF"].header
            zone = int(np.argmin((np.asarray(table["XCENTER"], dtype=float) - xp) ** 2 + (np.asarray(table["YCENTER"], dtype=float) - yp) ** 2))
            epsf = np.asarray(table["EPSF"][zone], dtype=float)
            oversampling = (int(ehdr["OVSMPY"]), int(ehdr["OVSMPX"]))
            bits = flag_bits(hdul["FLAGS"].header)
            omega_sr = _pixel_solid_angle_sr(sky, x, y)
            others_xy = [tuple(float(v) for v in sky.all_world2pix(nra, ndec, 0)) for nra, ndec in neighbours]
    if not bits or "SOURCE" not in bits:
        raise ValueError("FLAGS extension names no bit planes (MSKN cards)")
    if not (np.isfinite(wavelength_um) and 0.5 < wavelength_um < 6.0 and omega_sr > 0.0):
        raise ValueError(f"wavelength {wavelength_um} micron or pixel solid angle {omega_sr} sr is not usable")
    bad_mask = 0
    for name in BAD_FLAG_NAMES:
        bad_mask |= bits.get(name, 0)
    sl = (slice(iy - h, iy + h + 1), slice(ix - h, ix + h + 1))
    stamp_flags = flags[sl]
    good = ((stamp_flags & bad_mask) == 0) & np.isfinite(image[sl]) & np.isfinite(variance[sl]) & (variance[sl] > 0.0) & np.isfinite(zodi[sl])
    offs = np.arange(-h, h + 1, dtype=float)
    yy, xx = np.meshgrid(offs, offs, indexing="ij")
    core = np.hypot(xx - (x - ix), yy - (y - iy)) <= CORE_RADIUS_PIX
    row.update({
        "wavelength_um": wavelength_um, "bandwidth_um": bandwidth_um, "x_detector": round(xp, 2), "y_detector": round(yp, 2),
        "epsf_zone": zone, "epsf_cwave_um": float(table["CWAVE"][zone]), "pixel_solid_angle_sr": omega_sr,
        "flags_core": int(np.bitwise_or.reduce(stamp_flags[core])), "n_good_pixels": int(good.sum()),
    })
    if not bool(np.all(good[core])):
        return {**row, "status": "core_flagged"}
    if int(good.sum()) < MIN_GOOD_PIXELS:
        return {**row, "status": "too_few_pixels"}
    psf = render_epsf(epsf, oversampling, x - ix, y - iy)
    extra, n_blended = [], 0
    for ox, oy in others_xy:
        if not (np.isfinite(ox) and np.isfinite(oy)):
            continue
        if math.hypot(ox - x, oy - y) < BLEND_RADIUS_PIX:
            n_blended += 1
        elif abs(ox - ix) <= h + 1 and abs(oy - iy) <= h + 1:
            extra.append(render_epsf(epsf, oversampling, ox - ix, oy - iy))
    try:
        fit = fit_point_source(image[sl] - zodi[sl], variance[sl], good, psf, extra)
    except ValueError:
        return {**row, "n_blended": n_blended, "n_neighbours_fitted": len(extra), "status": "fit_undetermined"}
    to_mjy = omega_sr * 1.0e9   # MJy/sr x sr = MJy; x 1e9 = mJy
    # The formal error assumes the PSF model is exact. On a bright star it is not, and the misfit shows as a reduced
    # chi-square above one: the error is scaled up by its square root, never down.
    inflate = math.sqrt(max(1.0, fit["reduced_chi2"]))
    row.update({
        "status": "ok", "flux_mjy": fit["amplitude"] * to_mjy, "flux_err_mjy": fit["amplitude_err"] * inflate * to_mjy,
        "flux_err_formal_mjy": fit["amplitude_err"] * to_mjy, "background_mjy_sr": fit["background"],
        "reduced_chi2": fit["reduced_chi2"], "n_pixels": fit["n_pixels"], "psf_fraction_used": fit["psf_fraction_used"],
        "n_blended": n_blended, "n_neighbours_fitted": len(extra),
    })
    return row


# ------------------------------------------------------------------------------------------ archive access


def image_query(ra: float, dec: float) -> str:
    """ADQL for every SPHEREx spectral image whose footprint holds the position."""
    return ("SELECT a.uri, p.planeid, p.energy_bandpassname, p.time_bounds_lower, p.time_bounds_upper, o.collection, o.observationid "
            "FROM spherex.artifact a JOIN spherex.plane p ON a.planeid = p.planeid JOIN spherex.observation o ON p.obsid = o.obsid "
            f"WHERE 1 = CONTAINS(POINT('ICRS', {ra:.6f}, {dec:.6f}), p.poly)")


def list_images(ra: float, dec: float, cache: RawCache, *, as_of: str, timeout_s: float = 90.0,
                sleep: Sleeper = time.sleep) -> list[dict[str, Any]]:
    """Level 2 spectral images that hold the position, oldest first: uri, bandpass (SPHEREx-D1..D6), mjd, collection.

    The archive grows every week. `as_of` (a date or release label) is part of the cache key: the same label reads
    the same stored list again, a new label asks the archive again.
    """
    from oddstar.io.live_survey import _rows_from_tap_json

    query = image_query(ra, dec)
    sig = hashlib.sha256(f"{IRSA_TAP_URL}|{query}|{as_of}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(ARCHIVE, f"image_list:{as_of}", sig)
    rows: list[dict[str, Any]] | None = None
    if cached is not None and not cached.offline_fixture:
        try:
            rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        except (FileNotFoundError, OSError, ValueError):
            rows = None
    if rows is None:
        resp = ot._get_with_retries(IRSA_TAP_URL, {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": query},
                                    timeout_s=timeout_s, attempts=4, sleep=sleep)
        try:
            rows = _rows_from_tap_json(resp.json())
        except ValueError as exc:
            raise ArchiveAccessError(f"IRSA TAP answer for the SPHEREx image list is not readable: {exc}", retryable=True) from exc
        cache.put(archive=ARCHIVE, release=f"image_list:{as_of}", request_signature=sig, query_text=query,
                  payload=json.dumps(rows, sort_keys=True).encode("utf-8"), row_count=len(rows), offline_fixture=False)
    out = []
    for r in rows:
        uri = str(r.get("uri") or "")
        if "/level2/" not in uri or not uri.endswith(".fits") or r.get("time_bounds_lower") is None:
            continue
        out.append({"uri": uri.lstrip("/"), "bandpass": str(r.get("energy_bandpassname") or ""), "mjd": float(r["time_bounds_lower"]),
                    "collection": str(r.get("collection") or ""), "observationid": str(r.get("observationid") or "")})
    out.sort(key=lambda r: (r["mjd"], r["uri"]))
    return out


def cutout_url(uri: str, ra: float, dec: float, size_deg: float = CUTOUT_SIZE_DEG) -> str:
    return f"{IRSA_BASE_URL}{uri}?center={ra:.7f},{dec:.7f}d&size={size_deg}"


def _row_signature(uri: str, ra: float, dec: float, neighbours: tuple[tuple[float, float], ...]) -> tuple[str, str]:
    desc = json.dumps({"version": PHOTOMETRY_VERSION, "uri": uri, "ra": round(ra, 7), "dec": round(dec, 7), "size_deg": CUTOUT_SIZE_DEG,
                       "neighbours": [[round(a, 6), round(b, 6)] for a, b in neighbours]}, sort_keys=True)
    return hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24], desc


def measure_image(img: dict[str, Any], ra: float, dec: float, cache: RawCache, *, neighbours: tuple[tuple[float, float], ...] = (),
                  timeout_s: float = 60.0, sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """The measured row of one image at the position: from the cache, or from a fresh cutout (which is not kept)."""
    sig, desc = _row_signature(img["uri"], ra, dec, neighbours)
    cached = cache.lookup(ARCHIVE, PHOTOMETRY_VERSION, sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(cache.read_bytes(cached).decode("utf-8"))
        except (FileNotFoundError, OSError, ValueError):
            pass
    # IRSA answers 503 "Server busy, retry later." under load: four attempts, 5 to 20 s apart
    resp = ot._get_with_retries(cutout_url(img["uri"], ra, dec), None, timeout_s=timeout_s, attempts=4, backoff_s=5.0, sleep=sleep)
    try:
        row = measure_cutout(resp.content, ra, dec, neighbours)
    except (ValueError, KeyError, OSError) as exc:
        raise ArchiveAccessError(f"SPHEREx cutout {img['uri']} is not readable: {type(exc).__name__}: {exc}", retryable=False) from exc
    row.update({"uri": img["uri"], "bandpass": img["bandpass"], "collection": img["collection"], "observationid": img["observationid"],
                "ra_deg": round(ra, 7), "dec_deg": round(dec, 7)})
    cache.put(archive=ARCHIVE, release=PHOTOMETRY_VERSION, request_signature=sig, query_text=desc,
              payload=json.dumps(row, sort_keys=True).encode("utf-8"), row_count=1, offline_fixture=False)
    return row


def forced_photometry(source_id: str, ra: float, dec: float, cache: RawCache, *, as_of: str, pmra_mas_yr: float = 0.0,
                      pmdec_mas_yr: float = 0.0, ref_epoch: float = 2016.0, neighbours: tuple[tuple[float, float], ...] = (),
                      detectors: tuple[int, ...] | None = None, workers: int = 2, sleep: Sleeper = time.sleep,
                      log: Callable[[str], None] = lambda msg: None) -> dict[str, Any]:
    """Every SPHEREx measurement of one star: one row per exposure, and the summary of `summarize`.

    `ra`, `dec` are at `ref_epoch` (Gaia DR3: 2016.0); each exposure is measured at the position the proper motion
    gives for its date. `detectors` (1 to 6) limits the bands fetched. A cutout that cannot be fetched or read
    after the retries is counted under `n_failed` with its error, and `complete` is then false: the rows that were
    measured stay valid, and a later call fetches only what is missing.
    """
    # The footprint test uses the position near the middle of the survey so far; a footprint is 3.5 degrees wide
    ra_list, dec_list = propagate_position(ra, dec, pmra_mas_yr, pmdec_mas_yr, ref_epoch, LIST_EPOCH)
    images = list_images(ra_list, dec_list, cache, as_of=as_of, sleep=sleep)
    if detectors is not None:
        wanted = {f"SPHEREx-D{int(d)}" for d in detectors}
        images = [im for im in images if im["bandpass"] in wanted]

    def one(im: dict[str, Any]) -> dict[str, Any]:
        epoch = 2000.0 + (im["mjd"] - 51544.5) / 365.25
        ra_t, dec_t = propagate_position(ra, dec, pmra_mas_yr, pmdec_mas_yr, ref_epoch, epoch)
        try:
            return measure_image(im, ra_t, dec_t, cache, neighbours=neighbours, sleep=sleep)
        except ArchiveAccessError as exc:
            return {"uri": im["uri"], "bandpass": im["bandpass"], "mjd": im["mjd"], "status": "failed", "error": str(exc)[:300]}

    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        for n, row in enumerate(pool.map(one, images), 1):
            rows.append(row)
            if n % 25 == 0 or n == len(images):
                log(f"SPHEREx {source_id}: {n} of {len(images)} images")
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
    points = sorted((r for r in rows if r["status"] == "ok"), key=lambda r: (r["mjd"], r["wavelength_um"]))
    failed = [r for r in rows if r["status"] == "failed"]
    return {
        "source_id": str(source_id), "version": PHOTOMETRY_VERSION, "as_of": as_of, "ra_deg": ra, "dec_deg": dec, "ref_epoch": ref_epoch,
        "n_images": len(images), "n_measured": len(points), "n_failed": len(failed), "complete": not failed,
        "status_counts": counts, "errors": sorted({r["error"] for r in failed})[:5],
        "collections": sorted({im["collection"] for im in images}), "detectors": sorted({im["bandpass"] for im in images}),
        "n_neighbours": len(neighbours), "summary": summarize(points), "points": points,
        "not_measured": [{k: r.get(k) for k in ("uri", "status", "mjd")} for r in rows if r["status"] not in ("ok", "failed")],
    }


# ------------------------------------------------------------------------------------------ summary


def _visits(points: list[dict[str, Any]], gap_days: float = VISIT_GAP_DAYS) -> list[list[dict[str, Any]]]:
    out: list[list[dict[str, Any]]] = []
    for p in sorted(points, key=lambda r: r["mjd"]):
        if out and p["mjd"] - out[-1][-1]["mjd"] <= gap_days:
            out[-1].append(p)
        else:
            out.append([p])
    return out


def summarize(points: list[dict[str, Any]]) -> dict[str, Any]:
    """Window means (K, W1 part, W2 part) over all exposures and per visit, with the colours between them.

    A magnitude here is -2.5 log10 of a window mean over a Vega zero point (`WINDOW_ZERO_POINT_JY`): useful to
    compare one star with another measured the same way, and one visit with the next. It is not a 2MASS or WISE
    catalogue magnitude.
    """
    def block(sel: list[dict[str, Any]]) -> dict[str, Any]:
        out: dict[str, Any] = {"n": len(sel)}
        if sel:
            out.update({"mjd_min": float(min(p["mjd"] for p in sel)), "mjd_max": float(max(p["mjd"] for p in sel)),
                        "wavelength_min_um": float(min(p["wavelength_um"] for p in sel)),
                        "wavelength_max_um": float(max(p["wavelength_um"] for p in sel))})
        mags: dict[str, float] = {}
        for name, (lo, hi) in SED_WINDOWS_UM.items():
            avg = band_average(sel, lo, hi)
            out[name] = avg
            if avg is not None and avg["flux_mjy"] > 0.0:
                avg["mag"] = -2.5 * math.log10(avg["flux_mjy"] * 1.0e-3 / WINDOW_ZERO_POINT_JY[name])
                avg["mag_err"] = 2.5 / math.log(10.0) * avg["flux_err_mjy"] / avg["flux_mjy"]
                mags[name] = avg["mag"]
        out["k_w1"] = mags["k"] - mags["w1"] if "k" in mags and "w1" in mags else None
        out["w1_w2"] = mags["w1"] - mags["w2"] if "w1" in mags and "w2" in mags else None
        out["k_w2"] = mags["k"] - mags["w2"] if "k" in mags and "w2" in mags else None
        return out

    return {"all": block(points), "visits": [block(v) for v in _visits(points)]}


__all__ = ["ARCHIVE", "BAD_FLAG_NAMES", "CUTOUT_SIZE_DEG", "PHOTOMETRY_VERSION", "SED_WINDOWS_UM", "band_average", "cutout_url",
           "fit_point_source", "flag_bits", "forced_photometry", "image_query", "list_images", "measure_cutout", "measure_image",
           "render_epsf", "summarize"]
