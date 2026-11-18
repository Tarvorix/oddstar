"""Publication-grade evidence plot and fixed-scale multi-band cutout renderers."""

from __future__ import annotations

import hashlib
import math
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from oddstar.io.physics import (
    BAND_WAVELENGTH_UM,
    dust_blackbody_flux_mjy,
    stellar_photosphere_flux_mjy,
)


def render_sed_plot(
    star: dict[str, Any],
    measurement_values: dict[str, float],
    out_path: Path,
) -> Path:
    """Render observed SED (Gaia + 2MASS + WISE), fitted photosphere, and dust blackbody."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    teff = float(measurement_values.get("teff_fit_k", star.get("teff_k", 5800.0)))
    scale = float(measurement_values.get("star_scale_mjy", star["flux_mjy"]["Ks"]))
    t_dust = float(measurement_values.get("t_dust_k", 0.0))
    frac_lum = float(measurement_values.get("frac_lum", 0.0))

    lams = np.logspace(np.log10(0.42), np.log10(30.0), 200)
    f_star_curve = np.array([stellar_photosphere_flux_mjy(float(l), teff, scale) for l in lams])
    f_dust_curve = (
        np.array([dust_blackbody_flux_mjy(float(l), teff, scale, t_dust, frac_lum) for l in lams])
        if (t_dust > 0.0 and frac_lum > 0.0)
        else np.zeros_like(lams)
    )
    f_tot_curve = f_star_curve + f_dust_curve

    obs_lams = [BAND_WAVELENGTH_UM[b] for b in BAND_WAVELENGTH_UM]
    obs_flux = [float(star["flux_mjy"][b]) for b in BAND_WAVELENGTH_UM]
    obs_err = [float(star["flux_err_mjy"][b]) for b in BAND_WAVELENGTH_UM]

    fig, ax = plt.subplots(figsize=(7.5, 4.8), dpi=130)
    ax.loglog(lams, f_star_curve, "--", color="#1f77b4", lw=1.4, label=f"Photosphere (Teff={teff:.0f} K)")
    if np.any(f_dust_curve > 0):
        ax.loglog(lams, f_dust_curve, ":", color="#d62728", lw=1.4, label=f"Blackbody excess (T={t_dust:.0f} K, f={frac_lum:.4f})")
        ax.loglog(lams, f_tot_curve, "-", color="#222222", lw=1.5, alpha=0.8, label="Total model")
    ax.errorbar(obs_lams, obs_flux, yerr=obs_err, fmt="o", color="#000000", ecolor="#555555", ms=5, capsize=2, label="Observed (Gaia/2MASS/WISE)")
    ax.set_xlabel("Wavelength (um)")
    ax.set_ylabel("Flux density (mJy)")
    ax.set_title(f"SED — {star.get('name', star['source_id'])}")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="best")
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def render_cutouts_panel(
    star: dict[str, Any],
    out_path: Path,
    *,
    fov_arcsec: float = 40.0,
) -> Path:
    """Cutout evidence panel: real fetched pixels when available, a synthetic panel only for synthetic fixtures,
    otherwise a placeholder stating that no pixels were fetched (never invented sky images for a real star)."""
    if star.get("synthetic_fixture"):
        return render_synthetic_fixture_cutouts(star, out_path, fov_arcsec=fov_arcsec)
    real = render_real_cutouts_panel(star, out_path)
    if real is not None:
        return real
    return _render_no_pixels_placeholder(star, out_path)


def _render_no_pixels_placeholder(star: dict[str, Any], out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    qs = (star.get("vetting_archive") or {}).get("query_status") or {}
    fig, ax = plt.subplots(figsize=(6.0, 2.4), dpi=110)
    ax.axis("off")
    txt = (
        f"No cutout pixels fetched for {star.get('name', star.get('source_id'))}\n"
        f"cutouts_wise: {qs.get('cutouts_wise', 'not run')}\ncutouts_ps1: {qs.get('cutouts_ps1', 'not run')}"
    )
    ax.text(0.5, 0.5, txt, ha="center", va="center", fontsize=9, wrap=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def _stamp(blob_path: str, ra: float, dec: float, half_arcsec: float) -> tuple[np.ndarray, float, float, Any] | None:
    """Crop a cached FITS cutout to +/-half_arcsec around (ra, dec). Returns (stamp, x, y, wcs) or None."""
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales

    from oddstar.io.cutouts import read_fits_image

    try:
        data, header = read_fits_image(Path(blob_path).read_bytes())
        wcs = WCS(header).celestial
        pix = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
        x0, y0 = (float(v) for v in wcs.world_to_pixel_values(ra, dec))
    except Exception:
        return None
    h = max(3, int(round(half_arcsec / pix)))
    ix, iy = int(round(x0)), int(round(y0))
    x_lo, y_lo = max(0, ix - h), max(0, iy - h)
    stamp = data[y_lo : iy + h + 1, x_lo : ix + h + 1]
    if stamp.size == 0:
        return None
    return stamp, x0 - x_lo, y0 - y_lo, (wcs, x_lo, y_lo, pix)


def render_real_cutouts_panel(star: dict[str, Any], out_path: Path, *, half_fov_arcsec: float = 40.0) -> Path | None:
    """Render the fetched AllWISE Atlas W1-W4 (and PS1 stack) cutouts with the predicted Gaia position (cyan +),
    the fitted image centroid (lime x) and the band FWHM circle, annotated with the measured metrics.
    Returns None when no cached cutout pixels are readable."""
    cut = (star.get("vetting_archive") or {}).get("cutouts") or {}
    bands = cut.get("bands") or {}
    ra, dec = cut.get("pred_ra_deg"), cut.get("pred_dec_deg")
    if ra is None or dec is None:
        return None
    panels: list[tuple[str, tuple[np.ndarray, float, float, Any], dict[str, Any]]] = []
    ps1 = cut.get("ps1") or {}
    if ps1.get("blob_path"):
        st = _stamp(ps1["blob_path"], float(star["ra_deg"]), float(star["dec_deg"]), min(half_fov_arcsec, 25.0))
        if st is not None:
            panels.append((f"PS1 {ps1.get('filter', 'r')} stack", st, {}))
    for b in ("W1", "W2", "W3", "W4"):
        m = bands.get(b) or {}
        if not m.get("blob_path"):
            continue
        st = _stamp(m["blob_path"], float(ra), float(dec), half_fov_arcsec)
        if st is not None:
            panels.append((f"AllWISE Atlas {b}", st, m))
    if not panels:
        return None
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ncol = 3
    nrow = int(math.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 3.4 * nrow), dpi=120, squeeze=False)
    for ax in axes.flat:
        ax.axis("off")
    for ax, (title, (img, x0, y0, (wcs, x_lo, y_lo, pix)), m) in zip(axes.flat, panels):
        finite = img[np.isfinite(img)]
        vmin, vmax = (np.percentile(finite, [5.0, 99.5]) if finite.size else (0.0, 1.0))
        norm_img = np.arcsinh(np.clip((img - vmin) / max(1e-9, vmax - vmin), 0.0, 1.0) * 5.0)
        ax.imshow(norm_img, origin="lower", cmap="magma")
        ax.plot(x0, y0, "+", color="cyan", ms=12, mew=1.4)
        label = title
        if m.get("status") == "ok" and m.get("centroid_ra_deg") is not None:
            cx, cy = (float(v) for v in wcs.world_to_pixel_values(m["centroid_ra_deg"], m["centroid_dec_deg"]))
            ax.plot(cx - x_lo, cy - y_lo, "x", color="lime", ms=9, mew=1.4)
            if m.get("fwhm_ref_arcsec"):
                ax.add_patch(plt.Circle((x0, y0), 0.5 * float(m["fwhm_ref_arcsec"]) / pix, fill=False, edgecolor="lime", lw=0.9, ls="--"))
            sig = m.get("offset_sigma_arcsec")
            sig_txt = f"{float(sig):.2f}" if sig is not None else "n/a"
            snr = m.get("flux_snr")
            snr_txt = f"{float(snr):.1f}" if snr is not None else "n/a"
            ratio = m.get("fwhm_ratio")
            ratio_txt = f"{float(ratio):.2f}" if ratio is not None else "n/a"
            label += f"\noff {float(m['offset_arcsec']):.2f}\"±{sig_txt}  S/N {snr_txt}  FWHM/star {ratio_txt}"
        elif m:
            label += f"\n{m.get('status')}"
        ax.set_title(label, fontsize=7.5)
        ax.axis("on")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.suptitle(
        f"Real cutouts — {star.get('name', star.get('source_id'))} (+ Gaia at epoch {cut.get('epoch_jyear', float('nan')):.2f}, x fitted centroid)",
        fontsize=9,
    )
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def render_synthetic_fixture_cutouts(
    star: dict[str, Any],
    out_path: Path,
    *,
    fov_arcsec: float = 40.0,
) -> Path:
    """SYNTHETIC multi-band cutouts drawn from a fixture star's explicit inputs. Not real pixels: fixtures only."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    seed_val = int.from_bytes(hashlib.sha256(str(star["source_id"]).encode()).digest()[:4], "big")
    rng = np.random.default_rng(seed_val)
    n_pix = 64
    scale_arcsec = fov_arcsec / n_pix
    yy, xx = np.mgrid[0:n_pix, 0:n_pix]
    cx, cy = (n_pix - 1) / 2.0, (n_pix - 1) / 2.0

    vet_arch = star.get("vetting_archive", {})
    aw = star.get("allwise", {})
    offset_arcsec = float(aw.get("centroid_offset_arcsec", 0.2))
    bg_rms = float(aw.get("cirrus_bg_rms_mjy", 0.4))
    has_galaxy = vet_arch.get("deep_imaging_offset_galaxy") is not None
    gal_sep = float(vet_arch.get("deep_imaging_offset_galaxy", {}).get("sep_arcsec", offset_arcsec)) if has_galaxy else 0.0

    bands = [
        ("Deep Optical", 1.2, 0.0, False),
        ("2MASS Ks", 1.8, 0.0, False),
        ("WISE W1 (3.4 um)", 3.0, offset_arcsec * 0.4, True),
        ("WISE W2 (4.6 um)", 3.2, offset_arcsec * 0.6, True),
        ("WISE W3 (12 um)", 4.0, offset_arcsec * 0.9, True),
        ("WISE W4 (22 um)", 5.5, offset_arcsec, True),
    ]

    fig, axes = plt.subplots(2, 3, figsize=(9.0, 6.0), dpi=120)
    for ax, (title, psf_arcsec, dx_arcsec, show_beam) in zip(axes.flat, bands):
        psf_pix = max(1.0, psf_arcsec / scale_arcsec)
        dx_pix = dx_arcsec / scale_arcsec
        img = rng.normal(0.0, 0.08 + 0.03 * bg_rms, (n_pix, n_pix))
        if bg_rms > 3.0 and "WISE" in title:
            # Structured Galactic cirrus filament gradient
            img += 0.8 * np.sin((xx + yy) / 8.0)
        # Central star
        star_amp = 1.0 if not show_beam else min(2.5, float(star["flux_mjy"].get("W4", 1.0)) / 10.0 + 0.6)
        img += star_amp * np.exp(-0.5 * (((xx - (cx + dx_pix)) ** 2 + (yy - cy) ** 2) / (psf_pix**2)))
        # Optional background galaxy / neighbour
        if has_galaxy:
            gx = cx + gal_sep / scale_arcsec
            gy = cy + 0.3 * gal_sep / scale_arcsec
            img += 0.95 * np.exp(-0.5 * (((xx - gx) ** 2 + (yy - gy) ** 2) / ((psf_pix * 1.35) ** 2)))
        for nb in vet_arch.get("neighbours_in_w4_beam", []):
            sep = float(nb.get("sep_arcsec", 4.0)) / scale_arcsec
            img += 0.8 * np.exp(-0.5 * (((xx - (cx - sep)) ** 2 + (yy - (cy + 0.4 * sep)) ** 2) / (psf_pix**2)))

        # Standard asinh stretch
        vmin, vmax = np.percentile(img, [5.0, 99.5])
        norm_img = np.arcsinh(np.clip((img - vmin) / max(1e-5, vmax - vmin), 0.0, 1.0) * 5.0)
        ax.imshow(norm_img, origin="lower", cmap="magma")
        ax.plot(cx, cy, "+", color="cyan", ms=10, mew=1.4)
        if show_beam:
            beam_r_pix = (12.0 if "W4" in title else 6.5) / scale_arcsec
            circ = plt.Circle((cx, cy), beam_r_pix, fill=False, edgecolor="lime", lw=1.0, ls="--")
            ax.add_patch(circ)
        ax.set_title(title, fontsize=9)
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle(f"SYNTHETIC fixture cutouts (not real pixels, {fov_arcsec:.0f}\" FOV) — {star.get('name', star['source_id'])}", fontsize=11)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def render_lightcurves_plot(star: dict[str, Any], out_path: Path) -> Path:
    """Render all available optical and mid-infrared light curves on aligned panels."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lc = star.get("lightcurves", {})
    neo = star.get("neowise", {})

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.0, 7.0), dpi=120)
    t_tess = lc.get("tess_time_btjd", [])
    f_qlp = lc.get("tess_flux_qlp", [])
    f_cut = lc.get("tess_flux_tesscut", [])
    if t_tess:
        ax1.plot(t_tess, f_qlp, ".", ms=2.5, alpha=0.8, label="TESS QLP")
        ax1.plot(t_tess, f_cut, "-", lw=0.8, alpha=0.65, label="TESScut aperture")
        ax1.legend(fontsize=8)
    ax1.set_ylabel("TESS rel. flux")
    ax1.set_title(f"Time-Domain Photometry — {star.get('name', star['source_id'])}")
    ax1.grid(True, alpha=0.25)

    t_ztf = lc.get("ztf_time_mjd", [])
    f_ztf = lc.get("ztf_flux", [])
    if t_ztf:
        ax2.plot(t_ztf, f_ztf, "o", ms=2.5, color="#2ca02c", alpha=0.8, label="ZTF / ASAS-SN optical")
        ax2.legend(fontsize=8)
    ax2.set_ylabel("Optical rel. flux")
    ax2.grid(True, alpha=0.25)

    mjd_neo = neo.get("mjd", [])
    w1 = neo.get("w1_mjy", [])
    w2 = neo.get("w2_mjy", [])
    if mjd_neo:
        ax3.plot(mjd_neo, w1, "s-", ms=3.5, lw=1.0, label="NEOWISE W1 (3.4 um)")
        ax3.plot(mjd_neo, w2, "^-", ms=3.5, lw=1.0, label="NEOWISE W2 (4.6 um)")
        ax3.legend(fontsize=8)
    ax3.set_ylabel("Flux (mJy)")
    ax3.set_xlabel("Time (BTJD / MJD)")
    ax3.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def render_hr_diagram(
    star: dict[str, Any],
    sample_stars: list[dict[str, Any]],
    out_path: Path,
) -> Path:
    """Render Gaia color-magnitude HR diagram with parent sample and target star."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    colors = [float(s["mags"]["BP"]) - float(s["mags"]["RP"]) for s in sample_stars[:2000]]
    absmags = [float(s["abs_g_mag"]) for s in sample_stars[:2000]]

    t_col = float(star["mags"]["BP"]) - float(star["mags"]["RP"])
    t_mg = float(star["abs_g_mag"])

    fig, ax = plt.subplots(figsize=(6.2, 5.2), dpi=120)
    if colors:
        ax.scatter(colors, absmags, s=6, c="#888888", alpha=0.35, label="Parent sample")
    ax.scatter([t_col], [t_mg], s=95, c="#d62728", marker="*", zorder=5, label=star.get("name", star["source_id"]))
    ax.invert_yaxis()
    ax.set_xlabel("Gaia BP - RP (mag)")
    ax.set_ylabel("Absolute G magnitude (M_G)")
    ax.set_title("HR-Diagram Position")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def render_statistic_distribution(
    channel_name: str,
    target_stat: float,
    sample_stats: list[float],
    threshold: float,
    out_path: Path,
) -> Path:
    """Render parent-sample distribution of channel statistic with target star vertical line."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.asarray(sample_stats if sample_stats else [0.0, 0.5, -0.5], dtype=float)
    fig, ax = plt.subplots(figsize=(6.5, 4.2), dpi=120)
    ax.hist(arr, bins=40, density=True, color="#4c72b0", alpha=0.7, label="Parent sample")
    ax.axvline(threshold, color="#ff7f0e", ls="--", lw=1.4, label=f"Threshold ({threshold:.2f})")
    ax.axvline(target_stat, color="#d62728", lw=2.0, label=f"Target ({target_stat:.2f})")
    ax.set_xlabel(f"{channel_name} statistic")
    ax.set_ylabel("Density")
    ax.set_title(f"{channel_name} Statistic Distribution")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path


def render_waterfall_plot(
    cadence: dict[str, Any],
    out_path: Path,
) -> Path:
    """Render stacked ON/OFF radio filterbank waterfall cutout."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    scans = cadence.get("scans", [])
    fig, axes = plt.subplots(len(scans), 1, figsize=(7.0, 1.3 * max(1, len(scans))), dpi=110, sharex=True)
    if len(scans) == 1:
        axes = [axes]
    for ax, sc in zip(axes, scans):
        arr = np.asarray(sc["power"], dtype=float)
        ax.imshow(arr, aspect="auto", origin="lower", cmap="viridis")
        ax.set_ylabel(sc["kind"], fontsize=8)
        ax.set_yticks([])
    axes[-1].set_xlabel("Frequency channel")
    fig.suptitle(f"Radio Cadence Waterfall — {cadence['source_id']}", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    return out_path
