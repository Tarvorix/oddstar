"""Real-pixel image cutouts for infrared vetting: AllWISE Atlas W1-W4 (IRSA) and Pan-STARRS1 stacks (STScI).

Fetch:
- AllWISE Atlas tiles via IRSA SIA v2 `https://irsa.ipac.caltech.edu/SIA?COLLECTION=wise_allwise&POS=circle+ra+dec+r`
  (VOTable). For each band the intensity image (`-w<b>-int-3` in `access_url`) of the tile nearest the target is
  cut out through the IBE cutout syntax `<access_url>?center=<ra>,<dec>deg&size=<n>arcsec&gzip=false`.
- Pan-STARRS1 stack: `https://ps1images.stsci.edu/cgi-bin/ps1filenames.py?ra=&dec=&filters=r&type=stack` then
  `https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra=&dec=&size=<px>&format=fits&red=<filename>` (Dec > -30 only).
Every FITS payload goes through the content-addressed RawCache; the star dict only keeps the blob path/sha256
and the derived metrics, never pixels.

Measure (per WISE band, `measure_band`):
- Gaussian + tilted-plane fit around the Gaia position propagated to the WISE epoch; centroid offset from
  that position and a centroid uncertainty sigma = hypot(formal fit error, sigma_psf / SNR_empirical), where
  SNR_empirical uses the measured scatter of beam-sized aperture sums in a blank annulus. The Atlas pixels are
  correlated (the Atlas co-adds use the PSF as interpolation kernel), so the formal per-pixel error alone
  would be far too small; the empirical aperture scatter includes that correlation and any cirrus structure.
- Background: robust per-pixel rms and the clipped rms of aperture sums (aperture radius = FWHM_eff/2, which
  holds 50% of a Gaussian point source) converted with the image MAGZP to a point-source-equivalent 1-sigma
  flux in mJy ("mJy/beam" below).
- Extendedness: fitted FWHM / reference FWHM, with reference FWHM = native band PSF FWHM x Atlas broadening.
- Local maxima above `peak_snr_min` per-pixel sigma within one beam FWHM of the target, split into the peak
  at the fitted centroid and secondary peaks.
"""

from __future__ import annotations

import io
import math
import re
import time
from typing import Any

import numpy as np

from oddstar.errors import ArchiveAccessError
from oddstar.io.archive_http import Sleeper, cached_get
from oddstar.io.cache import RawCache
from oddstar.io.physics import BAND_ZERO_POINT_JY, angular_separation_arcsec
from oddstar.schemas import ManifestRecord

IRSA_SIA_URL = "https://irsa.ipac.caltech.edu/SIA"
PS1_FILENAMES_URL = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
PS1_FITSCUT_URL = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
PS1_PIXEL_ARCSEC = 0.25
WISE_BANDS = ("W1", "W2", "W3", "W4")
GAUSS_FWHM_PER_SIGMA = 2.0 * math.sqrt(2.0 * math.log(2.0))


def _s(val: Any) -> str:
    if isinstance(val, bytes):
        return val.decode("utf-8", errors="replace")
    return str(val)


def read_fits_image(payload: bytes) -> tuple[np.ndarray, Any]:
    """Return (2-D float image, header) of the first HDU holding image data."""
    from astropy.io import fits

    with fits.open(io.BytesIO(payload)) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, "ndim", 0) >= 2:
                data = np.asarray(hdu.data, dtype=float)
                while data.ndim > 2:
                    data = data[0]
                return data, hdu.header.copy()
    raise ArchiveAccessError("FITS payload holds no image HDU", retryable=False)


def _validate_fits(payload: bytes) -> None:
    if not payload.startswith(b"SIMPLE"):
        raise ArchiveAccessError(f"Cutout response is not FITS: {payload[:120]!r}", retryable=True)
    read_fits_image(payload)


def parse_allwise_sia(payload: bytes, ra: float, dec: float) -> dict[str, str]:
    """Map band -> access_url of the intensity image of the tile whose centre is nearest (ra, dec)."""
    from astropy.io.votable import parse

    try:
        table = parse(io.BytesIO(payload)).get_first_table().to_table()
    except Exception as exc:  # astropy raises several types on a non-VOTable body
        raise ArchiveAccessError(f"IRSA SIA response is not a VOTable: {type(exc).__name__}: {exc}", retryable=True) from exc
    cols = set(table.colnames)
    if "access_url" not in cols:
        raise ArchiveAccessError("IRSA SIA VOTable has no access_url column", retryable=False)
    best: dict[str, tuple[float, str]] = {}
    for row in table:
        url = _s(row["access_url"]).strip()
        m = re.search(r"-w([1-4])-int-3", url)
        if m is None:
            continue  # uncertainty / coverage / other products
        band = f"W{m.group(1)}"
        dist = 0.0
        if "s_ra" in cols and "s_dec" in cols:
            try:
                dist = angular_separation_arcsec(ra, dec, float(row["s_ra"]), float(row["s_dec"]))
            except (TypeError, ValueError):
                dist = 0.0
        if band not in best or dist < best[band][0]:
            best[band] = (dist, url)
    return {b: u for b, (_, u) in best.items()}


def fetch_allwise_cutouts(
    *,
    ra: float,
    dec: float,
    size_arcsec: float,
    cache: RawCache,
    timeout_s: float = 60.0,
    sleep: Sleeper = time.sleep,
) -> tuple[dict[str, dict[str, Any]], list[ManifestRecord], list[str]]:
    """Fetch AllWISE Atlas W1-W4 cutouts. Returns ({band: {data, header, blob_path, sha256}}, manifests, log)."""
    manifests: list[ManifestRecord] = []
    log: list[str] = []
    sia_payload, man, _ = cached_get(
        url=IRSA_SIA_URL,
        params={"COLLECTION": "wise_allwise", "POS": f"circle {ra:.6f} {dec:.6f} 0.002"},
        archive="allwise",
        release="atlas_sia",
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
        validate=lambda b: None if parse_allwise_sia(b, ra, dec) is not None else None,
    )
    manifests.append(man)
    urls = parse_allwise_sia(sia_payload, ra, dec)
    images: dict[str, dict[str, Any]] = {}
    for band in WISE_BANDS:
        url = urls.get(band)
        if url is None:
            log.append(f"{band}: no Atlas intensity image in SIA response")
            continue
        try:
            fb, man, _ = cached_get(
                url=url,
                params={"center": f"{ra:.6f},{dec:.6f}deg", "size": f"{size_arcsec:.0f}arcsec", "gzip": "false"},
                archive="allwise",
                release="atlas_cutout",
                cache=cache,
                timeout_s=timeout_s,
                sleep=sleep,
                validate=_validate_fits,
            )
            manifests.append(man)
            data, header = read_fits_image(fb)
            images[band] = {"data": data, "header": header, "blob_path": man.blob_path, "sha256": man.sha256}
            log.append(f"{band}: ok {data.shape[1]}x{data.shape[0]} px")
        except Exception as exc:
            log.append(f"{band}: {type(exc).__name__}: {str(exc)[:160]}")
    return images, manifests, log


def fetch_ps1_stack_cutout(
    *,
    ra: float,
    dec: float,
    size_arcsec: float,
    cache: RawCache,
    band: str = "r",
    timeout_s: float = 60.0,
    sleep: Sleeper = time.sleep,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Fetch one Pan-STARRS1 stack cutout (display and visual inspection only)."""
    manifests: list[ManifestRecord] = []
    listing, man, _ = cached_get(
        url=PS1_FILENAMES_URL,
        params={"ra": f"{ra:.6f}", "dec": f"{dec:.6f}", "filters": band, "type": "stack"},
        archive="deep_optical",
        release="ps1_stack_list",
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
    )
    manifests.append(man)
    lines = [ln.split() for ln in listing.decode("utf-8", errors="replace").splitlines() if ln.strip()]
    if len(lines) < 2 or "filename" not in lines[0]:
        raise ArchiveAccessError(f"PS1 ps1filenames.py returned no stack for filter {band}", retryable=False)
    fname = lines[1][lines[0].index("filename")]
    fb, man, _ = cached_get(
        url=PS1_FITSCUT_URL,
        params={"ra": f"{ra:.6f}", "dec": f"{dec:.6f}", "size": int(round(size_arcsec / PS1_PIXEL_ARCSEC)), "format": "fits", "red": fname},
        archive="deep_optical",
        release="ps1_stack_cutout",
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
        validate=_validate_fits,
    )
    manifests.append(man)
    data, _ = read_fits_image(fb)
    return {"status": "ok", "filter": band, "blob_path": man.blob_path, "sha256": man.sha256, "shape": list(data.shape)}, manifests


# ------------------------------------------------------------------------------------------- measurement


def _clipped_stats(vals: np.ndarray, nsig: float = 3.0, iters: int = 5) -> tuple[float, float]:
    v = vals[np.isfinite(vals)]
    if v.size == 0:
        return float("nan"), float("nan")
    for _ in range(iters):
        med = float(np.median(v))
        sig = float(1.4826 * np.median(np.abs(v - med)))
        if sig <= 0:
            break
        keep = np.abs(v - med) <= nsig * sig
        if keep.all():
            break
        v = v[keep]
    med = float(np.median(v))
    return med, float(1.4826 * np.median(np.abs(v - med)))


def _aperture_sums(img: np.ndarray, xx: np.ndarray, yy: np.ndarray, x0: float, y0: float, *, pix: float,
                   r_in: float, r_out: float, r_ap: float, bg: float) -> np.ndarray:
    """Sums of (img - bg) in non-overlapping apertures of radius r_ap (arcsec) tiling an annulus."""
    sums: list[float] = []
    r_ap_px = r_ap / pix
    rad = r_in + r_ap
    n_ap_px = math.pi * r_ap_px**2
    while rad + r_ap <= r_out:
        rad_px = rad / pix
        n_around = max(1, int(math.floor(2.0 * math.pi * rad / (2.0 * r_ap))))
        for k in range(n_around):
            ang = 2.0 * math.pi * k / n_around
            cx, cy = x0 + rad_px * math.cos(ang), y0 + rad_px * math.sin(ang)
            sel = (xx - cx) ** 2 + (yy - cy) ** 2 <= r_ap_px**2
            vals = img[sel]
            good = np.isfinite(vals)
            if sel.sum() == 0 or good.sum() < 0.8 * n_ap_px:
                continue
            sums.append(float(np.sum(vals[good] - bg) * (sel.sum() / good.sum())))
        rad += 2.0 * r_ap
    return np.asarray(sums, dtype=float)


def measure_band(
    data: np.ndarray,
    header: Any,
    *,
    band: str,
    ra_pred: float,
    dec_pred: float,
    native_fwhm_arcsec: float,
    psf_broadening: float,
    beam_arcsec: float,
    annulus_arcsec: tuple[float, float],
    peak_snr_min: float,
) -> dict[str, Any]:
    """Quantitative image metrics for one band (see module docstring)."""
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    from scipy.ndimage import maximum_filter
    from scipy.optimize import least_squares

    img = np.asarray(data, dtype=float)
    wcs = WCS(header).celestial
    pix = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
    magzp = header.get("MAGZP")
    dn_mjy = (BAND_ZERO_POINT_JY[band] * 1000.0 * 10.0 ** (-0.4 * float(magzp))) if magzp is not None else None
    x0, y0 = (float(v) for v in wcs.world_to_pixel_values(ra_pred, dec_pred))
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(float)
    r_arc = np.hypot(xx - x0, yy - y0) * pix
    out: dict[str, Any] = {"band": band, "pix_scale_arcsec": pix, "magzp": None if magzp is None else float(magzp),
                           "pred_ra_deg": ra_pred, "pred_dec_deg": dec_pred}
    if not (0 <= x0 < nx and 0 <= y0 < ny):
        out["status"] = "target_off_image"
        return out

    r_in, r_out = annulus_arcsec
    ann = np.isfinite(img) & (r_arc >= r_in) & (r_arc <= r_out)
    if ann.sum() < 50:
        out["status"] = "annulus_off_image"
        return out
    # Remove a large-scale tilted plane fitted to the annulus (iteratively clipped) so a smooth gradient across
    # the cutout is not counted as beam-scale background structure.
    ax, ay, az = xx[ann] - x0, yy[ann] - y0, img[ann]
    keep = np.ones(az.size, dtype=bool)
    coef = np.array([float(np.median(az)), 0.0, 0.0])
    for _ in range(5):
        amat = np.column_stack([np.ones(keep.sum()), ax[keep], ay[keep]])
        coef, *_ = np.linalg.lstsq(amat, az[keep], rcond=None)
        resid = az - (coef[0] + coef[1] * ax + coef[2] * ay)
        sig = 1.4826 * float(np.median(np.abs(resid[keep] - np.median(resid[keep]))))
        new_keep = np.abs(resid) <= 3.0 * max(sig, 1e-12)
        if (new_keep == keep).all():
            break
        keep = new_keep
    plane = coef[0] + coef[1] * (xx - x0) + coef[2] * (yy - y0)
    flat = img - plane
    _, bg_sig = _clipped_stats(flat[ann])
    bg_med = float(coef[0])
    fwhm_ref = native_fwhm_arcsec * psf_broadening
    r_ap = fwhm_ref / 2.0
    ap = _aperture_sums(flat, xx, yy, x0, y0, pix=pix, r_in=r_in, r_out=r_out, r_ap=r_ap, bg=0.0)
    _, ap_sig = _clipped_stats(ap) if ap.size >= 8 else (float("nan"), float("nan"))
    # point-source-equivalent 1-sigma noise: aperture of radius FWHM/2 holds 50% of a Gaussian's flux
    ps_sig_dn = ap_sig / 0.5 if np.isfinite(ap_sig) else float("nan")
    out.update(
        {
            "bg_median_dn": bg_med,
            "bg_rms_pix_dn": bg_sig,
            "n_bg_apertures": int(ap.size),
            "bg_rms_mjy_beam": (ps_sig_dn * dn_mjy) if (dn_mjy is not None and np.isfinite(ps_sig_dn)) else None,
            "fwhm_ref_arcsec": fwhm_ref,
        }
    )

    # Gaussian + tilted plane fit in a box of +-1.5 FWHM_ref
    half = max(3, int(math.ceil(1.5 * fwhm_ref / pix)))
    xlo, xhi = max(0, int(round(x0)) - half), min(nx, int(round(x0)) + half + 1)
    ylo, yhi = max(0, int(round(y0)) - half), min(ny, int(round(y0)) + half + 1)
    sub = img[ylo:yhi, xlo:xhi]
    sx, sy = xx[ylo:yhi, xlo:xhi], yy[ylo:yhi, xlo:xhi]
    good = np.isfinite(sub)
    if good.sum() < 20:
        out["status"] = "fit_box_empty"
        return out
    zs, xs, ys = sub[good], sx[good], sy[good]
    sig0 = fwhm_ref / GAUSS_FWHM_PER_SIGMA / pix

    def model(p: np.ndarray) -> np.ndarray:
        a, xc, yc, s, c0, c1, c2 = p
        return a * np.exp(-((xs - xc) ** 2 + (ys - yc) ** 2) / (2.0 * s * s)) + c0 + c1 * (xs - x0) + c2 * (ys - y0)

    p0 = np.array([max(1e-6, float(np.nanmax(sub)) - bg_med), x0, y0, sig0, bg_med, 0.0, 0.0])
    lo = [-np.inf, x0 - half, y0 - half, 0.3 * sig0, -np.inf, -np.inf, -np.inf]
    hi = [np.inf, x0 + half, y0 + half, 4.0 * sig0, np.inf, np.inf, np.inf]
    try:
        res = least_squares(lambda p: model(p) - zs, p0, bounds=(lo, hi))
    except Exception as exc:
        out["status"] = f"fit_failed: {type(exc).__name__}"
        return out
    a, xc, yc, s, _, _, _ = (float(v) for v in res.x)
    dof = max(1, zs.size - 7)
    s2 = float(np.sum(res.fun**2)) / dof
    try:
        cov = np.linalg.inv(res.jac.T @ res.jac) * max(s2, bg_sig**2)
        ex, ey = math.sqrt(max(0.0, cov[1, 1])), math.sqrt(max(0.0, cov[2, 2]))
    except np.linalg.LinAlgError:
        ex = ey = float("nan")
    flux_dn = 2.0 * math.pi * a * s * s
    snr = flux_dn / ps_sig_dn if (np.isfinite(ps_sig_dn) and ps_sig_dn > 0) else float("nan")
    sig_formal = math.sqrt(0.5 * (ex * ex + ey * ey)) * pix if np.isfinite(ex) else float("nan")
    sig_emp = (fwhm_ref / GAUSS_FWHM_PER_SIGMA) / snr if (np.isfinite(snr) and snr > 0) else float("inf")
    sig_c = math.hypot(sig_formal if np.isfinite(sig_formal) else 0.0, sig_emp)
    c_ra, c_dec = (float(v) for v in wcs.pixel_to_world_values(xc, yc))
    offset = angular_separation_arcsec(ra_pred, dec_pred, c_ra, c_dec)
    fwhm_fit = s * GAUSS_FWHM_PER_SIGMA * pix

    # Local maxima within one beam FWHM of the target
    fp = max(3, int(round(0.5 * fwhm_ref / pix)) | 1)
    filled = np.where(np.isfinite(flat), flat, -np.inf)
    peaks_mask = (filled == maximum_filter(filled, size=fp, mode="nearest")) & np.isfinite(flat)
    peaks_mask &= flat >= peak_snr_min * bg_sig
    peaks_mask &= r_arc <= beam_arcsec
    pk_y, pk_x = np.nonzero(peaks_mask)
    peaks = []
    for px_, py_ in zip(pk_x, pk_y):
        d_cen = math.hypot(px_ - xc, py_ - yc) * pix
        peaks.append({"dist_from_centroid_arcsec": d_cen, "dist_from_target_arcsec": float(r_arc[py_, px_]),
                      "value_dn": float(flat[py_, px_]), "snr_pix": float(flat[py_, px_] / bg_sig)})
    primary = [p for p in peaks if p["dist_from_centroid_arcsec"] < 0.5 * fwhm_ref]
    secondary = [p for p in peaks if p["dist_from_centroid_arcsec"] >= 0.5 * fwhm_ref]
    if primary:
        prim_val = max(p["value_dn"] for p in primary)
    else:
        ix, iy = int(round(min(max(xc, 0), nx - 1))), int(round(min(max(yc, 0), ny - 1)))
        prim_val = float(flat[iy, ix])
    secondary.sort(key=lambda p: -p["value_dn"])
    out.update(
        {
            "status": "ok",
            "centroid_ra_deg": c_ra,
            "centroid_dec_deg": c_dec,
            "offset_arcsec": offset,
            "offset_sigma_arcsec": sig_c,
            "offset_sigma_formal_arcsec": sig_formal,
            "fwhm_arcsec": fwhm_fit,
            "fwhm_ratio": fwhm_fit / fwhm_ref,
            "flux_dn": flux_dn,
            "flux_mjy": flux_dn * dn_mjy if dn_mjy is not None else None,
            "flux_snr": snr,
            "n_secondary_peaks": len(secondary),
            "brighter_secondary_peak": any(p["value_dn"] > prim_val for p in secondary),
            "secondary_peaks": [{k: round(v, 3) for k, v in p.items()} for p in secondary[:5]],
            "primary_peak_value_dn": prim_val,
        }
    )
    return out


def measure_wise_cutouts(
    images: dict[str, dict[str, Any]],
    *,
    ra_pred: float,
    dec_pred: float,
    vetting: dict[str, float],
) -> dict[str, dict[str, Any]]:
    """Run `measure_band` on each fetched WISE band with thresholds/PSF values from `cfg.vetting`."""
    native = {
        "W1": float(vetting["w1_beam_arcsec"]),
        "W2": float(vetting["w2_beam_arcsec"]),
        "W3": float(vetting["w3_beam_arcsec"]),
        "W4": float(vetting["w4_beam_arcsec"]),
    }
    out: dict[str, dict[str, Any]] = {}
    for band, im in images.items():
        try:
            m = measure_band(
                im["data"],
                im["header"],
                band=band,
                ra_pred=ra_pred,
                dec_pred=dec_pred,
                native_fwhm_arcsec=native[band],
                psf_broadening=float(vetting["wise_atlas_psf_broadening"]),
                beam_arcsec=native[band],
                annulus_arcsec=(float(vetting["cirrus_annulus_inner_arcsec"]), float(vetting["cirrus_annulus_outer_arcsec"])),
                peak_snr_min=float(vetting["vision_peak_snr_min"]),
            )
        except Exception as exc:
            m = {"band": band, "status": f"measure_failed: {type(exc).__name__}: {str(exc)[:120]}"}
        m["blob_path"] = im.get("blob_path")
        m["sha256"] = im.get("sha256")
        out[band] = {k: _json_safe(v) for k, v in m.items()}
    return out


def _json_safe(v: Any) -> Any:
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, (float, np.floating)):
        f = float(v)
        return f if math.isfinite(f) else None
    if isinstance(v, np.integer):
        return int(v)
    return v
