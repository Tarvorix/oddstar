"""Digitized plate images at the position of a plate source: is there a point of light where the catalogue says?

The vanishing search takes its plate sources from USNO-B1.0, a catalogue made by matching detections across five
plates. It contains false sources: pieces of the halos and spikes of bright stars, plate flaws that two plates
happen to share, and detections of a moving star matched wrongly. The catalogue row cannot tell these from a star
that was there. The plates can: the Digitized Sky Survey serves the scans the catalogue was made from.

For one position and one survey this module fetches a 2 arcmin cutout from STScI (about 70 kB, not kept), measures
it and stores the measured row. A row says what the plate shows at the position:

- `snr_at_position`: the mean density above sky in a `APERTURE_RADIUS_ARCSEC` aperture at the position, in units of
  the scatter of the same aperture placed at `SKY_ANNULUS_ARCSEC` around it. Plate grain is correlated from pixel
  to pixel, so the noise of an aperture is measured, not computed from the pixel scatter.
- `peak`: the local maximum nearest to the position within `SEARCH_RADIUS_ARCSEC`, with its centroid offset, the
  same signal-to-noise taken at its centroid, and its width beside the median width of the other stars of the
  cutout (`field_fwhm_median_arcsec`). The image's own sky solution is used as delivered; on the first-epoch
  plates it is good to about 1.5 arcsec (measured on 30 fields, 2026-11-18).

The position is the catalogue position moved to the date of each plate with the catalogue's proper motion.

Nothing here decides whether a source is real. The limits on signal-to-noise and offset live in the config and are
applied by the vetting check that reads these rows.

Surveys (`v` parameter of the service) and what they return: `poss1_red` and `poss1_blue` are the first-epoch
Palomar plates (1949 to 1958) north of Dec -30 deg or so. `poss2ukstu_red`, `poss2ukstu_blue` and `poss2ukstu_ir`
return the second-epoch plates: Palomar POSS-II in the north, and UK Schmidt plates in the south, whose blue
survey (SERC-J) dates from 1975 to 1983. The `SURVEY` header key says which plate came back.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import time
from collections.abc import Callable
from typing import Any

import numpy as np

from oddstar.errors import ArchiveAccessError
from oddstar.io import optical_timeseries as ot
from oddstar.io.cache import RawCache
from oddstar.io.cutouts import _clipped_stats, read_fits_image
from oddstar.io.physics import propagate_position

ARCHIVE = "historical_plates"
MEASUREMENT_VERSION = "dss_plate_v1"
DSS_URL = "https://archive.stsci.edu/cgi-bin/dss_search"
FIRST_EPOCH_SURVEYS = ("poss1_red", "poss1_blue")
SECOND_EPOCH_SURVEYS = ("poss2ukstu_red", "poss2ukstu_blue", "poss2ukstu_ir")
SURVEYS = FIRST_EPOCH_SURVEYS + SECOND_EPOCH_SURVEYS
CUTOUT_ARCMIN = 2.0
APERTURE_RADIUS_ARCSEC = 2.5
SKY_ANNULUS_ARCSEC = (8.0, 20.0)
SEARCH_RADIUS_ARCSEC = 8.0
PEAK_MIN_SIGMA = 3.0  # a local maximum counts when it stands this many pixel scatters above the sky
FIELD_STAR_MIN_SNR = 10.0  # stars of the cutout whose width is the reference for a point source
MIN_NOISE_APERTURES = 8
PLATE_YEARS = (1940.0, 2005.0)

Sleeper = Callable[[float], Any]


def cutout_url(survey: str, ra: float, dec: float, size_arcmin: float = CUTOUT_ARCMIN) -> str:
    return (f"{DSS_URL}?v={survey}&r={ra:.6f}&d={dec:.6f}&e=J2000&h={size_arcmin:.1f}&w={size_arcmin:.1f}"
            "&f=fits&c=none&fov=NONE&v3=")


def plate_epoch(header: Any) -> float | None:
    """Decimal year of the plate. Only the date of `DATE-OBS` is read: its time of day is often not a valid time."""
    m = re.match(r"\s*(\d{4})-(\d{2})-(\d{2})", str(header.get("DATE-OBS") or ""))
    if m:
        year, month, day = int(m.group(1)), int(m.group(2)), int(m.group(3))
        if 1 <= month <= 12 and 1 <= day <= 31:
            return round(year + ((month - 1) * 30.4375 + (day - 1)) / 365.25, 3)
    try:
        epoch = float(header.get("EPOCH"))
    except (TypeError, ValueError):
        return None
    return round(epoch, 3) if PLATE_YEARS[0] <= epoch <= PLATE_YEARS[1] else None


def _aperture_mean(img: np.ndarray, xx: np.ndarray, yy: np.ndarray, x0: float, y0: float, r_px: float) -> float | None:
    sel = (xx - x0) ** 2 + (yy - y0) ** 2 <= r_px**2
    return float(img[sel].mean()) if sel.sum() >= 3 else None


def _moments(img: np.ndarray, xx: np.ndarray, yy: np.ndarray, x0: float, y0: float, r_px: float, sky: float) -> tuple[float, float, float] | None:
    """Centroid (x, y) and FWHM in pixels of the light above sky within `r_px` of a pixel."""
    sel = (xx - x0) ** 2 + (yy - y0) ** 2 <= r_px**2
    w = np.clip(img[sel] - sky, 0.0, None)
    total = float(w.sum())
    if total <= 0:
        return None
    cx, cy = float((w * xx[sel]).sum() / total), float((w * yy[sel]).sum() / total)
    var = float((w * ((xx[sel] - cx) ** 2 + (yy[sel] - cy) ** 2)).sum() / total) / 2.0
    return cx, cy, 2.3548 * math.sqrt(max(var, 0.0))


def _local_maxima(img: np.ndarray, floor: float) -> list[tuple[int, int]]:
    """Pixels (y, x) that are the maximum of their 5x5 neighbourhood in the 3x3-smoothed image and exceed `floor`."""
    from scipy.ndimage import maximum_filter, uniform_filter

    smooth = uniform_filter(img, size=3, mode="nearest")
    peaks = (smooth == maximum_filter(smooth, size=5, mode="nearest")) & (smooth > floor)
    return [(int(y), int(x)) for y, x in zip(*np.nonzero(peaks))]


def measure_plate_cutout(payload: bytes, ra: float, dec: float, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
                         position_epoch: float | None = None) -> dict[str, Any]:
    """What one plate cutout shows at a position. `ra`, `dec` are at `position_epoch`; with a proper motion the
    position is moved to the plate's date before it is measured.

    `status` is "ok", "off_cutout" (the position is closer to the edge than the sky annulus needs),
    "no_plate_date" (a proper motion was given and the plate has no date) or "no_noise_estimate".
    """
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    if not payload.startswith(b"SIMPLE"):
        raise ValueError(f"not a FITS file: {payload[:120]!r}")
    img, hdr = read_fits_image(payload)
    epoch = plate_epoch(hdr)
    row: dict[str, Any] = {
        "status": "ok", "survey_returned": str(hdr.get("SURVEY") or "").strip() or None, "plate_id": str(hdr.get("PLATEID") or "").strip() or None,
        "region": str(hdr.get("REGION") or "").strip() or None, "bandpass": hdr.get("BANDPASS"), "emulsion": str(hdr.get("EMULSION") or "").strip() or None,
        "filter": str(hdr.get("FILTER") or "").strip() or None, "date_obs": str(hdr.get("DATE-OBS") or "").strip() or None,
        "plate_epoch": epoch, "exposure_min": hdr.get("EXPOSURE"), "plate_grade": str(hdr.get("PLTGRADE") or "").strip() or None,
        "shape": [int(img.shape[0]), int(img.shape[1])],
    }
    moving = bool(pmra_mas_yr or pmdec_mas_yr)
    if moving and (epoch is None or position_epoch is None):
        row["status"] = "no_plate_date"
        return row
    if moving:
        ra, dec = propagate_position(ra, dec, pmra_mas_yr, pmdec_mas_yr, float(position_epoch), float(epoch))
    row["ra_measured_deg"], row["dec_measured_deg"] = round(float(ra), 7), round(float(dec), 7)
    wcs = WCS(hdr).celestial
    pix = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
    x0, y0 = (float(v) for v in wcs.all_world2pix(ra, dec, 0))
    row.update({"pixel_arcsec": round(pix, 3), "x": round(x0, 2), "y": round(y0, 2)})
    ny, nx = img.shape
    r_out = SKY_ANNULUS_ARCSEC[1] / pix
    if not (r_out <= x0 <= nx - 1 - r_out and r_out <= y0 <= ny - 1 - r_out):
        row["status"] = "off_cutout"
        return row
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    rr = np.hypot(xx - x0, yy - y0) * pix
    sky, sky_sig = _clipped_stats(img[(rr >= SKY_ANNULUS_ARCSEC[0]) & (rr <= SKY_ANNULUS_ARCSEC[1])])
    r_ap = APERTURE_RADIUS_ARCSEC / pix
    # The same aperture on rings around the position: the scatter of its mean is the noise of one measurement
    ring_means: list[float] = []
    rad = SKY_ANNULUS_ARCSEC[0] + APERTURE_RADIUS_ARCSEC
    while rad + APERTURE_RADIUS_ARCSEC <= SKY_ANNULUS_ARCSEC[1]:
        n_around = max(1, int(math.floor(2.0 * math.pi * rad / (2.0 * APERTURE_RADIUS_ARCSEC))))
        for k in range(n_around):
            ang = 2.0 * math.pi * k / n_around
            m = _aperture_mean(img, xx, yy, x0 + rad / pix * math.cos(ang), y0 + rad / pix * math.sin(ang), r_ap)
            if m is not None:
                ring_means.append(m)
        rad += 2.0 * APERTURE_RADIUS_ARCSEC
    ring_med, ring_sig = _clipped_stats(np.asarray(ring_means)) if ring_means else (float("nan"), float("nan"))
    row.update({"sky": round(sky, 1), "sky_sigma": round(sky_sig, 2), "n_noise_apertures": len(ring_means),
                "aperture_noise": None if not math.isfinite(ring_sig) else round(ring_sig, 3)})
    if len(ring_means) < MIN_NOISE_APERTURES or not math.isfinite(ring_sig) or ring_sig <= 0 or not math.isfinite(sky_sig) or sky_sig <= 0:
        row["status"] = "no_noise_estimate"
        return row

    def snr_at(x: float, y: float) -> float | None:
        m = _aperture_mean(img, xx, yy, x, y, r_ap)
        return None if m is None else round((m - ring_med) / ring_sig, 2)

    row["snr_at_position"] = snr_at(x0, y0)
    maxima = _local_maxima(img, sky + PEAK_MIN_SIGMA * sky_sig)
    found: list[dict[str, Any]] = []
    for py, px in maxima:
        mom = _moments(img, xx, yy, float(px), float(py), r_ap + 1.0, sky)
        if mom is None:
            continue
        cx, cy, fwhm_px = mom
        found.append({"offset_arcsec": round(math.hypot(cx - x0, cy - y0) * pix, 2), "dx_arcsec": round((cx - x0) * pix, 2),
                      "dy_arcsec": round((cy - y0) * pix, 2), "snr": snr_at(cx, cy), "fwhm_arcsec": round(fwhm_px * pix, 2),
                      "peak_above_sky": round(float(img[py, px] - sky), 1)})
    near = sorted((p for p in found if p["offset_arcsec"] <= SEARCH_RADIUS_ARCSEC), key=lambda p: p["offset_arcsec"])
    row["peak"] = near[0] if near else None
    row["n_peaks_in_search_radius"] = len(near)
    # Width of a point source on this plate: the other bright stars of the cutout
    others = [p["fwhm_arcsec"] for p in found if p["offset_arcsec"] > SEARCH_RADIUS_ARCSEC and (p["snr"] or 0.0) >= FIELD_STAR_MIN_SNR]
    row["n_field_stars"] = len(others)
    row["field_fwhm_median_arcsec"] = round(float(np.median(others)), 2) if others else None
    return row


def _signature(survey: str, ra: float, dec: float, pmra: float, pmdec: float, position_epoch: float | None) -> tuple[str, str]:
    desc = json.dumps({"version": MEASUREMENT_VERSION, "survey": survey, "ra": round(ra, 6), "dec": round(dec, 6), "pmra": round(pmra, 2),
                       "pmdec": round(pmdec, 2), "position_epoch": position_epoch, "size_arcmin": CUTOUT_ARCMIN}, sort_keys=True)
    return hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24], desc


def measure_position(survey: str, ra: float, dec: float, cache: RawCache, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
                     position_epoch: float | None = None, timeout_s: float = 30.0, sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """The measured row of one survey at a position: from the cache, or from a fresh cutout (which is not kept).

    The cutout is centred on the position as given; it is 2 arcmin wide, so the position at the plate's date lies
    inside it for any proper motion a catalogue holds. A response that is not a FITS image raises
    `ArchiveAccessError` and nothing is stored.
    """
    if survey not in SURVEYS:
        raise ValueError(f"unknown plate survey {survey!r}; one of {SURVEYS}")
    sig, desc = _signature(survey, ra, dec, pmra_mas_yr, pmdec_mas_yr, position_epoch)
    cached = cache.lookup(ARCHIVE, MEASUREMENT_VERSION, sig)
    if cached is not None and not cached.offline_fixture:
        try:
            return json.loads(cache.read_bytes(cached).decode("utf-8"))
        except (FileNotFoundError, OSError, ValueError):
            pass
    resp = ot._get_with_retries(cutout_url(survey, ra, dec), None, timeout_s=timeout_s, attempts=3, backoff_s=3.0, sleep=sleep)
    try:
        row = measure_plate_cutout(resp.content, ra, dec, pmra_mas_yr=pmra_mas_yr, pmdec_mas_yr=pmdec_mas_yr, position_epoch=position_epoch)
    except (ValueError, KeyError, OSError, ArchiveAccessError) as exc:
        raise ArchiveAccessError(f"plate cutout {survey} at {ra:.5f} {dec:.5f} is not readable: {type(exc).__name__}: {str(exc)[:200]}",
                                 retryable=False) from exc
    row.update({"survey": survey, "ra_deg": round(ra, 6), "dec_deg": round(dec, 6), "measurement_version": MEASUREMENT_VERSION})
    cache.put(archive=ARCHIVE, release=MEASUREMENT_VERSION, request_signature=sig, query_text=desc,
              payload=json.dumps(row, sort_keys=True).encode("utf-8"), row_count=1, offline_fixture=False)
    return row


def measure_all_surveys(ra: float, dec: float, cache: RawCache, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
                        position_epoch: float | None = None, surveys: tuple[str, ...] = SURVEYS, min_interval_s: float = 0.0,
                        sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """One row per survey at a position, and whether all of them were read.

    A survey that could not be fetched or read is listed in `errors` and has no row; `complete` is then false.
    """
    rows: list[dict[str, Any]] = []
    errors: dict[str, str] = {}
    for survey in surveys:
        t0 = time.monotonic()
        try:
            rows.append(measure_position(survey, ra, dec, cache, pmra_mas_yr=pmra_mas_yr, pmdec_mas_yr=pmdec_mas_yr,
                                         position_epoch=position_epoch, sleep=sleep))
        except ArchiveAccessError as exc:
            errors[survey] = str(exc)[:300]
        wait = min_interval_s - (time.monotonic() - t0)
        if wait > 0:
            sleep(wait)
    return {"measurement_version": MEASUREMENT_VERSION, "surveys": list(surveys), "rows": rows, "errors": errors, "complete": not errors}
