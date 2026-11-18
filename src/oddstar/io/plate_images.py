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
  cutout (`field_fwhm_median_arcsec`).

The sky solution delivered with a cutout is a local approximation of the plate solution. On 203 ordinary stars
(2026-11-18) it put one in eight first-epoch red images 2.5 to 3.5 arcsec from the star. So each cutout is first
registered on the Gaia stars it holds (`REG_G_RANGE`, moved to the plate's date): the median offset between where
the solution puts them and where their images are is added to the position before it is measured. A row says how
many stars carried the registration (`registration`), and `registered` is false when fewer than `REG_MIN_STARS`
did; the position is then measured where the delivered solution puts it.

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
from oddstar.io.physics import angular_separation_arcsec, propagate_position

ARCHIVE = "historical_plates"
MEASUREMENT_VERSION = "dss_plate_v2"
DSS_URL = "https://archive.stsci.edu/cgi-bin/dss_search"
FIRST_EPOCH_SURVEYS = ("poss1_red", "poss1_blue")
SECOND_EPOCH_SURVEYS = ("poss2ukstu_red", "poss2ukstu_blue", "poss2ukstu_ir")
SURVEYS = FIRST_EPOCH_SURVEYS + SECOND_EPOCH_SURVEYS
CUTOUT_ARCMIN = 4.0  # wide enough to hold a handful of Gaia stars to register the image on
APERTURE_RADIUS_ARCSEC = 2.5
SKY_ANNULUS_ARCSEC = (8.0, 20.0)
SEARCH_RADIUS_ARCSEC = 8.0
PEAK_MIN_SIGMA = 3.0  # a local maximum counts when it stands this many pixel scatters above the sky
FIELD_STAR_MIN_SNR = 10.0  # stars of the cutout whose width is the reference for a point source
MIN_NOISE_APERTURES = 8
# Registration of a cutout on Gaia stars (the delivered sky solution is off by up to 3.5 arcsec on first-epoch plates)
REG_G_RANGE = (12.0, 19.0)  # brighter stars are saturated on the plates, fainter ones are not reliably there
REG_MATCH_ARCSEC = 6.0  # a reference star's image is looked for this far from where the delivered solution puts it
REG_EXCLUDE_ARCSEC = 10.0  # reference stars this close to the measured position, or to each other, are not used
REG_MIN_PEAK_SIGMA = 10.0
REG_MIN_STARS = 3
GAIA_EPOCH = 2016.0
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


def _register(bright: list[tuple[float, float]], refs_xy: list[tuple[float, float]], x0: float, y0: float, pix: float) -> dict[str, Any]:
    """Median offset (pixels) between where the delivered solution puts the reference stars (`refs_xy`) and their
    images (`bright`, centroids of the cutout's bright peaks). (x0, y0) is the position to be measured."""
    excl = REG_EXCLUDE_ARCSEC / pix
    usable = [(xr, yr) for i, (xr, yr) in enumerate(refs_xy)
              if math.hypot(xr - x0, yr - y0) > excl
              and all(math.hypot(xr - xo, yr - yo) > excl for j, (xo, yo) in enumerate(refs_xy) if j != i)]
    offsets: list[tuple[float, float]] = []
    for xr, yr in usable:
        near = [(cx, cy) for cx, cy in bright if math.hypot(cx - xr, cy - yr) * pix <= REG_MATCH_ARCSEC]
        if len(near) == 1:  # two bright peaks near one prediction: which is the star cannot be told
            offsets.append((near[0][0] - xr, near[0][1] - yr))
    out: dict[str, Any] = {"n_references": len(refs_xy), "n_usable": len(usable), "n_stars": len(offsets)}
    if len(offsets) < REG_MIN_STARS:
        return dict(out, registered=False, dx_px=0.0, dy_px=0.0)
    arr = np.asarray(offsets)
    dx, dy = float(np.median(arr[:, 0])), float(np.median(arr[:, 1]))
    resid = np.hypot(arr[:, 0] - dx, arr[:, 1] - dy) * pix
    return dict(out, registered=True, dx_px=dx, dy_px=dy, dx_arcsec=round(dx * pix, 2), dy_arcsec=round(dy * pix, 2),
                scatter_arcsec=round(float(1.4826 * np.median(resid)), 2), max_residual_arcsec=round(float(resid.max()), 2))


def measure_plate_cutout(payload: bytes, ra: float, dec: float, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
                         position_epoch: float | None = None, references: tuple[tuple[float, float, float, float], ...] = ()) -> dict[str, Any]:
    """What one plate cutout shows at a position. `ra`, `dec` are at `position_epoch`; with a proper motion the
    position is moved to the plate's date before it is measured.

    `references` are Gaia stars (ra, dec, pmra, pmdec at epoch 2016.0) to register the image on; see the module text.

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
    row.update({"pixel_arcsec": round(pix, 3), "x_delivered": round(x0, 2), "y_delivered": round(y0, 2)})
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    r_ap = APERTURE_RADIUS_ARCSEC / pix

    # Registration: the bright peaks of the whole cutout, against its own sky level and pixel scatter
    g_sky, g_sig = _clipped_stats(img.ravel())
    bright: list[tuple[float, float]] = []
    if references and math.isfinite(g_sig) and g_sig > 0:
        for py, px in _local_maxima(img, g_sky + REG_MIN_PEAK_SIGMA * g_sig):
            mom = _moments(img, xx, yy, float(px), float(py), r_ap + 1.0, g_sky)
            if mom is not None:
                bright.append((mom[0], mom[1]))
    refs_xy: list[tuple[float, float]] = []
    for r_ra, r_dec, r_pmra, r_pmdec in references:
        if epoch is not None:
            r_ra, r_dec = propagate_position(r_ra, r_dec, r_pmra, r_pmdec, GAIA_EPOCH, float(epoch))
        xr, yr = (float(v) for v in wcs.all_world2pix(r_ra, r_dec, 0))
        if 3.0 <= xr <= nx - 4.0 and 3.0 <= yr <= ny - 4.0:
            refs_xy.append((xr, yr))
    reg = _register(bright, refs_xy, x0, y0, pix)
    x0, y0 = x0 + reg.pop("dx_px"), y0 + reg.pop("dy_px")
    row.update({"registered": reg.pop("registered"), "registration": reg, "x": round(x0, 2), "y": round(y0, 2)})

    r_out = SKY_ANNULUS_ARCSEC[1] / pix
    if not (r_out <= x0 <= nx - 1 - r_out and r_out <= y0 <= ny - 1 - r_out):
        row["status"] = "off_cutout"
        return row
    rr = np.hypot(xx - x0, yy - y0) * pix
    sky, sky_sig = _clipped_stats(img[(rr >= SKY_ANNULUS_ARCSEC[0]) & (rr <= SKY_ANNULUS_ARCSEC[1])])
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


References = tuple[tuple[float, float, float, float], ...]


def _signature(survey: str, ra: float, dec: float, pmra: float, pmdec: float, position_epoch: float | None,
               references: References = ()) -> tuple[str, str]:
    refs = hashlib.sha256(json.dumps([[round(v, 6) for v in r] for r in references]).encode("utf-8")).hexdigest()[:16]
    desc = json.dumps({"version": MEASUREMENT_VERSION, "survey": survey, "ra": round(ra, 6), "dec": round(dec, 6), "pmra": round(pmra, 2),
                       "pmdec": round(pmdec, 2), "position_epoch": position_epoch, "size_arcmin": CUTOUT_ARCMIN,
                       "n_references": len(references), "references": refs}, sort_keys=True)
    return hashlib.sha256(desc.encode("utf-8")).hexdigest()[:24], desc


def reference_stars(gaia_rows: list[dict[str, Any]]) -> References:
    """The Gaia stars a cutout can be registered on, from rows with `ra_deg`, `dec_deg`, proper motions and `g_mag`."""
    out = [(float(g["ra_deg"]), float(g["dec_deg"]), float(g.get("pmra_mas_yr") or 0.0), float(g.get("pmdec_mas_yr") or 0.0))
           for g in gaia_rows if g.get("g_mag") is not None and REG_G_RANGE[0] <= float(g["g_mag"]) <= REG_G_RANGE[1]]
    return tuple(sorted(out))


def stored_row(survey: str, ra: float, dec: float, cache: RawCache, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
               position_epoch: float | None = None, references: References = ()) -> dict[str, Any] | None:
    """The row measured earlier for this survey, position and set of reference stars, or None."""
    sig, _ = _signature(survey, ra, dec, pmra_mas_yr, pmdec_mas_yr, position_epoch, references)
    cached = cache.lookup(ARCHIVE, MEASUREMENT_VERSION, sig)
    if cached is None or cached.offline_fixture:
        return None
    try:
        return json.loads(cache.read_bytes(cached).decode("utf-8"))
    except (FileNotFoundError, OSError, ValueError):
        return None


def measure_position(survey: str, ra: float, dec: float, cache: RawCache, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
                     position_epoch: float | None = None, references: References = (), timeout_s: float = 30.0,
                     sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """The measured row of one survey at a position: from the cache, or from a fresh cutout (which is not kept).

    The cutout is centred on the position as given; it is 4 arcmin wide, so the position at the plate's date lies
    inside it for any proper motion a catalogue holds. A response that is not a FITS image raises
    `ArchiveAccessError` and nothing is stored.
    """
    if survey not in SURVEYS:
        raise ValueError(f"unknown plate survey {survey!r}; one of {SURVEYS}")
    stored = stored_row(survey, ra, dec, cache, pmra_mas_yr=pmra_mas_yr, pmdec_mas_yr=pmdec_mas_yr, position_epoch=position_epoch,
                        references=references)
    if stored is not None:
        return stored
    sig, desc = _signature(survey, ra, dec, pmra_mas_yr, pmdec_mas_yr, position_epoch, references)
    resp = ot._get_with_retries(cutout_url(survey, ra, dec), None, timeout_s=timeout_s, attempts=3, backoff_s=3.0, sleep=sleep)
    try:
        row = measure_plate_cutout(resp.content, ra, dec, pmra_mas_yr=pmra_mas_yr, pmdec_mas_yr=pmdec_mas_yr, position_epoch=position_epoch,
                                   references=references)
    except (ValueError, KeyError, OSError, ArchiveAccessError) as exc:
        raise ArchiveAccessError(f"plate cutout {survey} at {ra:.5f} {dec:.5f} is not readable: {type(exc).__name__}: {str(exc)[:200]}",
                                 retryable=False) from exc
    row.update({"survey": survey, "ra_deg": round(ra, 6), "dec_deg": round(dec, 6), "measurement_version": MEASUREMENT_VERSION})
    cache.put(archive=ARCHIVE, release=MEASUREMENT_VERSION, request_signature=sig, query_text=desc,
              payload=json.dumps(row, sort_keys=True).encode("utf-8"), row_count=1, offline_fixture=False)
    return row


def measure_all_surveys(ra: float, dec: float, cache: RawCache, *, pmra_mas_yr: float = 0.0, pmdec_mas_yr: float = 0.0,
                        position_epoch: float | None = None, references: References = (), surveys: tuple[str, ...] = SURVEYS,
                        min_interval_s: float = 0.0, sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """One row per survey at a position, and whether all of them were read.

    A survey that could not be fetched or read is listed in `errors` and has no row; `complete` is then false.
    """
    rows: list[dict[str, Any]] = []
    errors: dict[str, str] = {}
    for survey in surveys:
        t0 = time.monotonic()
        kw = {"pmra_mas_yr": pmra_mas_yr, "pmdec_mas_yr": pmdec_mas_yr, "position_epoch": position_epoch, "references": references}
        row = stored_row(survey, ra, dec, cache, **kw) if survey in SURVEYS else None
        if row is not None:  # nothing was asked of the archive, so nothing to wait for
            rows.append(row)
            continue
        try:
            rows.append(measure_position(survey, ra, dec, cache, sleep=sleep, **kw))
        except ArchiveAccessError as exc:
            errors[survey] = str(exc)[:300]
        wait = min_interval_s - (time.monotonic() - t0)
        if wait > 0:
            sleep(wait)
    return {"measurement_version": MEASUREMENT_VERSION, "surveys": list(surveys), "rows": rows, "errors": errors, "complete": not errors,
            "n_reference_stars": len(references)}


def field_references(ra: float, dec: float, cfg: Any, cache: RawCache, *, sleep: Sleeper | None = None) -> References:
    """The Gaia stars of the cutout around a position, to register it on (one cached cone query)."""
    from oddstar.io.plate_survey import gaia_cone

    rows, _, _ = gaia_cone(ra, dec, CUTOUT_ARCMIN * 60.0 / 2.0, cfg, cache, sleep=sleep)
    return reference_stars(rows)


def enrich_unit(unit: dict[str, Any], cache: RawCache, *, cfg: Any = None, min_interval_s: float = 0.5,
                sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """Attach what the plates show to a plate-source unit of the vanishing survey, for vetting.

    Writes `vanishing_data.plate_images` and `vetting_archive.query_status.plate_images`: "ok" when both first-epoch
    surveys returned a readable image, "failed: ..." otherwise (the rows that were read are kept), and
    "no_plate_record" for a unit without its catalogue row. The position is the catalogue position with the
    catalogue's proper motion, moved to each plate's date. With `cfg` the cutouts are registered on the Gaia
    stars around the position; when that query fails the status says so and nothing is measured, so that a row
    is never a registered one on one day and an unregistered one on the next.
    """
    u = unit.get("usno") or {}
    qs = unit.setdefault("vetting_archive", {}).setdefault("query_status", {})
    if u.get("ra_deg") is None or u.get("dec_deg") is None:
        qs["plate_images"] = "no_plate_record"
        return unit
    ra, dec = float(u["ra_deg"]), float(u["dec_deg"])
    references: References = ()
    if cfg is not None:
        try:
            references = field_references(ra, dec, cfg, cache, sleep=sleep)
        except ArchiveAccessError as exc:
            qs["plate_images"] = f"failed: reference stars: {str(exc)[:160]}"
            return unit
    res = measure_all_surveys(ra, dec, cache, pmra_mas_yr=float(u.get("pmra_mas_yr") or 0.0),
                              pmdec_mas_yr=float(u.get("pmdec_mas_yr") or 0.0), position_epoch=u.get("position_epoch"),
                              references=references, min_interval_s=min_interval_s, sleep=sleep)
    unit.setdefault("vanishing_data", {})["plate_images"] = res
    missing = [sv for sv in FIRST_EPOCH_SURVEYS if sv in res["errors"]]
    qs["plate_images"] = "ok" if not missing else "failed: " + "; ".join(f"{sv}: {res['errors'][sv][:120]}" for sv in missing)
    if cfg is not None:
        attach_first_place(unit, cfg, cache, sleep=sleep)
        attach_neighbours(unit, cfg, cache, references, min_interval_s=min_interval_s, sleep=sleep)
    return unit


FIRST_PLACE_CONE_ARCSEC = 10.0  # not a limit of a check: wide enough for the offset limit plus the match radius


def attach_first_place(unit: dict[str, Any], cfg: Any, cache: RawCache, *, sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """The present-day catalogue objects around the place where the first-epoch plates were measured.

    The survey compares a plate source with Gaia and Pan-STARRS at the catalogue position. USNO-B1.0 gives some
    sources a false proper motion (detections of different objects joined into one track), and the catalogue
    position is then several arcsec from where the source was on the first plates. This reads the Gaia stars, moved
    back to the first plate's date, and the Pan-STARRS objects within `FIRST_PLACE_CONE_ARCSEC` of the position the
    first-epoch plates were measured at, so that vetting can repeat the survey's comparison at the place where the
    plates show the source. Writes `vanishing_data.plate_first_place` and
    `vetting_archive.query_status.plate_first_place` ("ok", "no_first_epoch_image", "no_coverage" south of the
    Pan-STARRS footprint, or "failed: ...").
    """
    from oddstar.io.plate_survey import GAIA_EPOCH as GAIA_POSITION_EPOCH, PS1_DEC_LIMIT_DEG, _float, gaia_cone, ps1_cone

    u = unit.get("usno") or {}
    qs = unit.setdefault("vetting_archive", {}).setdefault("query_status", {})
    vc = cfg.channels["vanishing"]
    rows = [r for r in ((unit.get("vanishing_data") or {}).get("plate_images") or {}).get("rows") or []
            if isinstance(r, dict) and r.get("survey") in FIRST_EPOCH_SURVEYS and r.get("ra_measured_deg") is not None
            and r.get("plate_epoch") is not None]
    if not rows:
        qs["plate_first_place"] = "no_first_epoch_image"
        return unit
    ra, dec, epoch = float(rows[0]["ra_measured_deg"]), float(rows[0]["dec_measured_deg"]), float(rows[0]["plate_epoch"])
    if dec < PS1_DEC_LIMIT_DEG:
        qs["plate_first_place"] = "no_coverage"
        return unit
    half = CUTOUT_ARCMIN * 60.0 / 2.0
    try:
        # The cone the cutouts were registered on (stored) holds every star that was within reach of the place
        centre = (float(u["ra_deg"]), float(u["dec_deg"]))
        if angular_separation_arcsec(centre[0], centre[1], ra, dec) > half / 2.0:
            centre = (ra, dec)
        gaia, _, _ = gaia_cone(centre[0], centre[1], half, cfg, cache, sleep=sleep)
        ps1, _, _ = ps1_cone(ra, dec, FIRST_PLACE_CONE_ARCSEC, int(vc["plate_ps1_min_detections"]), cache, sleep=sleep)
    except ArchiveAccessError as exc:
        qs["plate_first_place"] = f"failed: {str(exc)[:160]}"
        return unit
    g_out: list[dict[str, Any]] = []
    for g in gaia:
        then_ra, then_dec = propagate_position(g["ra_deg"], g["dec_deg"], g["pmra_mas_yr"], g["pmdec_mas_yr"], GAIA_POSITION_EPOCH, epoch)
        sep = angular_separation_arcsec(ra, dec, then_ra, then_dec)
        if sep <= FIRST_PLACE_CONE_ARCSEC:
            g_out.append({"source_id": g["source_id"], "ra_then_deg": round(then_ra, 7), "dec_then_deg": round(then_dec, 7),
                          "sep_then_arcsec": round(sep, 2), "pmra_mas_yr": g["pmra_mas_yr"], "pmdec_mas_yr": g["pmdec_mas_yr"],
                          "g_mag": g.get("g_mag")})
    p_out: list[dict[str, Any]] = []
    for p in ps1:
        p_ra, p_dec = _float(p.get("raMean")), _float(p.get("decMean"))
        if p_ra is None or p_dec is None:
            continue
        sep = angular_separation_arcsec(ra, dec, p_ra, p_dec)
        if sep <= FIRST_PLACE_CONE_ARCSEC:
            r_mag = _float(p.get("rMeanPSFMag"))
            p_out.append({"objID": str(p.get("objID")), "ra_deg": p_ra, "dec_deg": p_dec, "sep_arcsec": round(sep, 2),
                          "n_detections": p.get("nDetections"), "r_mag": r_mag if r_mag is not None and 5.0 < r_mag < 30.0 else None})
    unit.setdefault("vanishing_data", {})["plate_first_place"] = {
        "position": [round(ra, 7), round(dec, 7)], "epoch": epoch, "radius_arcsec": FIRST_PLACE_CONE_ARCSEC,
        "gaia": sorted(g_out, key=lambda g: (g["sep_then_arcsec"], g["source_id"])),
        "ps1": sorted(p_out, key=lambda p: (p["sep_arcsec"], p["objID"]))}
    qs["plate_first_place"] = "ok"
    return unit


MAX_NEIGHBOURS = 4


def attach_neighbours(unit: dict[str, Any], cfg: Any, cache: RawCache, references: References, *, min_interval_s: float = 0.5,
                      sleep: Sleeper = time.sleep) -> dict[str, Any]:
    """What the first-epoch plates show where the nearby Pan-STARRS objects are now.

    A plate source with no modern counterpart at its position may be a star that has since moved a few arcsec:
    Pan-STARRS then holds an object nearby that the first-epoch plates do not show at its present place. This reads
    the Pan-STARRS objects between the match radius and `plate_moved_radius_arcsec` of the plate position (the
    nearest `MAX_NEIGHBOURS`) and measures the first-epoch plates at each. Writes
    `vanishing_data.plate_neighbours` and `vetting_archive.query_status.plate_neighbours`.
    """
    from oddstar.io.plate_survey import PS1_DEC_LIMIT_DEG, _float, ps1_cone

    u = unit.get("usno") or {}
    qs = unit.setdefault("vetting_archive", {}).setdefault("query_status", {})
    vc = cfg.channels["vanishing"]
    ra = float(u.get("plate_ra_deg", u["ra_deg"]))
    dec = float(u.get("plate_dec_deg", u["dec_deg"]))
    if dec < PS1_DEC_LIMIT_DEG:
        qs["plate_neighbours"] = "no_coverage"
        return unit
    radius, match = float(vc["plate_moved_radius_arcsec"]), float(vc["match_radius_arcsec"])
    try:
        ps1, _, _ = ps1_cone(ra, dec, radius, int(vc["plate_ps1_min_detections"]), cache, sleep=sleep)
    except ArchiveAccessError as exc:
        qs["plate_neighbours"] = f"failed: Pan-STARRS: {str(exc)[:160]}"
        return unit
    near: list[dict[str, Any]] = []
    for p in ps1:
        p_ra, p_dec = _float(p.get("raMean")), _float(p.get("decMean"))
        if p_ra is None or p_dec is None:
            continue
        sep = angular_separation_arcsec(ra, dec, p_ra, p_dec)
        if match < sep <= radius:
            r_mag = _float(p.get("rMeanPSFMag"))
            near.append({"objID": str(p.get("objID")), "ra_deg": p_ra, "dec_deg": p_dec, "sep_arcsec": round(sep, 2),
                         "n_detections": p.get("nDetections"), "r_mag": r_mag if r_mag is not None and 5.0 < r_mag < 30.0 else None})
    near.sort(key=lambda n: (n["sep_arcsec"], n["objID"]))
    failed: list[str] = []
    for n in near[:MAX_NEIGHBOURS]:
        res = measure_all_surveys(n["ra_deg"], n["dec_deg"], cache, references=references, surveys=FIRST_EPOCH_SURVEYS,
                                  min_interval_s=min_interval_s, sleep=sleep)
        n["rows"], n["errors"] = res["rows"], res["errors"]
        failed += [f"{n['objID']} {sv}: {msg[:100]}" for sv, msg in res["errors"].items()]
    unit.setdefault("vanishing_data", {})["plate_neighbours"] = {
        "measurement_version": MEASUREMENT_VERSION, "position": [round(ra, 7), round(dec, 7)], "radius_arcsec": radius,
        "inner_radius_arcsec": match, "n_objects_in_radius": len(near), "neighbours": near[:MAX_NEIGHBOURS]}
    qs["plate_neighbours"] = "ok" if not failed else "failed: " + "; ".join(failed[:3])
    return unit
