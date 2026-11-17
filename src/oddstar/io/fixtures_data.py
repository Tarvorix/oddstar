"""Deterministic astrophysical star catalog, named benchmarks, impostors, and open cases.

All coordinates, proper motions, SEDs, light curves, spectra, radio filterbanks, and
vetting archive tables are computed consistently from astrophysical parameters so offline
validation, reproduction, and testing exercise exact pipeline mathematics without shortcuts.
"""

from __future__ import annotations

import hashlib
import math
from typing import Any

import numpy as np

from oddstar.io.physics import (
    BAND_EXTINCTION_AV,
    BAND_WAVELENGTH_UM,
    dust_blackbody_flux_mjy,
    flux_mjy_to_mag,
    mag_err_to_flux_err_mjy,
    propagate_position,
    stellar_photosphere_flux_mjy,
)
from oddstar.io.storage import equatorial_to_galactic, healpix_nest_to_radec


def _make_base_star(
    source_id: str,
    *,
    name: str = "",
    ra_deg: float | None = None,
    dec_deg: float | None = None,
    distance_pc: float = 120.0,
    teff_k: float = 5800.0,
    logg: float = 4.4,
    feh: float = 0.0,
    luminosity_lsun: float | None = None,
    ruwe: float = 1.02,
    pmra: float = 12.0,
    pmdec: float = -8.5,
    ebv_sfd: float = 0.02,
    ebv_planck: float = 0.02,
    t_dust_k: float | None = None,
    frac_lum: float = 0.0,
    w4_extra_mult: float = 1.0,
    gaia_release: str = "dr3",
    simbad_type: str = "Star",
    known_in_literature: bool = False,
    seed_offset: int = 0,
    flags: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Generate a complete, physically consistent multi-archive record for one star."""
    seed_val = int.from_bytes(hashlib.sha256(f"{source_id}:{seed_offset}".encode()).digest()[:4], "big")
    rng = np.random.default_rng(seed_val)

    # Derive coordinates from Gaia source_id HEALPix if numeric and not explicitly given
    if ra_deg is None or dec_deg is None:
        if source_id.isdigit() and int(source_id) > (1 << 35):
            r_dec = healpix_nest_to_radec(int(source_id) >> 35)
            ra_deg = r_dec[0] if ra_deg is None else ra_deg
            dec_deg = r_dec[1] if dec_deg is None else dec_deg
        else:
            ra_deg = float(rng.uniform(0.0, 360.0)) if ra_deg is None else ra_deg
            dec_deg = float(rng.uniform(-75.0, 75.0)) if dec_deg is None else dec_deg

    _, gal_b = equatorial_to_galactic(ra_deg, dec_deg)

    # Mass-luminosity relation on main sequence if not overridden
    if luminosity_lsun is None:
        luminosity_lsun = max(0.01, (teff_k / 5778.0) ** 4.8)

    parallax_mas = 1000.0 / max(1.0, distance_pc)
    parallax_snr = float(rng.uniform(45.0, 220.0))
    if gaia_release == "dr4":
        parallax_snr *= 1.41
    parallax_err = parallax_mas / parallax_snr

    # Apparent Ks flux scale from L / d^2
    # At d=10 pc, Sun (L=1, Teff=5778) has Ks ~ 3.28 => F_Ks ~ 32500 mJy
    star_scale_mjy = 32500.0 * luminosity_lsun * ((10.0 / distance_pc) ** 2) * ((5778.0 / teff_k) ** 2.2)

    av_true = 3.1 * ebv_sfd
    flux_mjy: dict[str, float] = {}
    flux_err_mjy: dict[str, float] = {}
    mags: dict[str, float] = {}
    mag_errs: dict[str, float] = {}

    rel_err = {
        "G": 0.004,
        "BP": 0.008,
        "RP": 0.008,
        "J": 0.020,
        "H": 0.022,
        "Ks": 0.024,
        "W1": 0.025,
        "W2": 0.028,
        "W3": 0.045,
        "W4": 0.085,
    }

    for band, lam in BAND_WAVELENGTH_UM.items():
        f_star = stellar_photosphere_flux_mjy(lam, teff_k, star_scale_mjy)
        ext_factor = 10.0 ** (-0.4 * av_true * BAND_EXTINCTION_AV[band])
        f_total = f_star * ext_factor
        if t_dust_k is not None and frac_lum > 0.0 and band in {"W1", "W2", "W3", "W4"}:
            f_dust = dust_blackbody_flux_mjy(lam, teff_k, star_scale_mjy, t_dust_k, frac_lum)
            if band == "W4":
                f_dust *= w4_extra_mult
            f_total += f_dust
        sig_frac = rel_err[band]
        if band == "W4" and f_total < 5.0:
            sig_frac = min(0.22, max(0.085, 0.45 / max(1.0, f_total)))
        sigma = f_total * sig_frac
        f_obs = max(0.05, float(f_total + rng.normal(0.0, sigma)))
        flux_mjy[band] = f_obs
        flux_err_mjy[band] = max(0.02, sigma)
        mags[band] = flux_mjy_to_mag(band, f_obs)
        mag_errs[band] = min(0.35, 1.0857 * flux_err_mjy[band] / max(1e-6, f_obs))

    # Absolute G magnitude consistent with parallax & extinction
    abs_g = mags["G"] - 5.0 * math.log10(distance_pc / 10.0) - av_true * BAND_EXTINCTION_AV["G"]

    fdict = dict(flags or {})
    # Compute WISE J2010 propagated position
    wise_ra, wise_dec = propagate_position(ra_deg, dec_deg, pmra, pmdec, 2016.0, 2010.2)
    centroid_offset_arcsec = float(fdict.get("centroid_offset_arcsec", abs(rng.normal(0.18, 0.10))))

    # Single-exposure 2010 W4 frames (12 frames)
    n_frames = 12
    w4_frame_base = flux_mjy["W4"]
    if fdict.get("single_frame_glitch"):
        w4_frames = [
            float(stellar_photosphere_flux_mjy(22.088, teff_k, star_scale_mjy) + rng.normal(0, 0.2))
            for _ in range(n_frames - 1)
        ] + [w4_frame_base * 3.0]
    else:
        w4_frames = [
            max(0.05, float(w4_frame_base + rng.normal(0.0, flux_err_mjy["W4"] * 1.4)))
            for _ in range(n_frames)
        ]

    # NEOWISE-R visits (2014-2024, 20 visits)
    visit_mjds = np.linspace(56700.0, 60450.0, 20)
    w1_visits = np.full(20, flux_mjy["W1"]) + rng.normal(0.0, flux_err_mjy["W1"], 20)
    w2_visits = np.full(20, flux_mjy["W2"]) + rng.normal(0.0, flux_err_mjy["W2"], 20)
    if fdict.get("ir_var_amp"):
        amp = float(fdict["ir_var_amp"])
        mode = fdict.get("ir_var_mode", "ramp")
        if mode == "burst":
            w1_visits[10:15] *= 1.0 + amp
            w2_visits[10:15] *= 1.0 + amp * 1.25
        elif mode == "step":
            w1_visits[10:] *= 1.0 + amp
            w2_visits[10:] *= 1.0 + amp * 1.2
        elif mode == "lpv":
            w1_visits *= 1.0 + amp * np.sin(2.0 * np.pi * (visit_mjds - 56700.0) / 420.0)
            w2_visits *= 1.0 + amp * np.sin(2.0 * np.pi * (visit_mjds - 56700.0) / 420.0)
        else:
            ramp = np.linspace(0.0, amp, 20)
            w1_visits *= 1.0 + ramp
            w2_visits *= 1.0 + ramp * 1.25

    # Optical light curves (TESS 27d @ 30min + ZTF/ASAS-SN 6yr)
    tess_time = np.linspace(2000.0, 2027.0, 648)  # 1 hour cadence for compact offline speed
    tess_flux_a = 1.0 + rng.normal(0.0, 0.0015, len(tess_time))
    tess_flux_b = 1.0 + rng.normal(0.0, 0.0016, len(tess_time))
    ztf_time = np.linspace(58200.0, 60400.0, 240)
    ztf_flux = 1.0 + rng.normal(0.0, 0.008, len(ztf_time))

    if fdict.get("optical_dips"):
        dip_spec = fdict["optical_dips"]
        depth = float(dip_spec.get("depth", 0.20))
        asym = float(dip_spec.get("asymmetry", 0.4))
        t0 = float(dip_spec.get("t0", 2012.0))
        dur = float(dip_spec.get("duration_days", 1.8))
        for i, t in enumerate(tess_time):
            dt = t - t0
            if -0.5 * dur <= dt <= dur:
                sigma_w = dur * (0.25 if dt < 0 else 0.25 * (1.0 + asym * 2.0))
                dip = depth * math.exp(-0.5 * (dt / max(0.05, sigma_w)) ** 2)
                tess_flux_a[i] -= dip
                if not fdict.get("detrend_artifact_only"):
                    tess_flux_b[i] -= dip
        # Also inject into ZTF if long-baseline
        if dip_spec.get("in_ztf", True):
            ztf_dip_idx = (ztf_time > 59350.0) & (ztf_time < 59600.0)
            ztf_flux[ztf_dip_idx] -= depth * 0.85

    if fdict.get("periodic_eb"):
        period = float(fdict.get("eb_period_days", 2.4))
        phase = ((tess_time - 2000.5) % period) / period
        eclipse = np.exp(-0.5 * ((phase - 0.5) / 0.03) ** 2) * 0.25
        tess_flux_a -= eclipse
        tess_flux_b -= eclipse
        z_phase = ((ztf_time - 58200.5) % period) / period
        ztf_flux -= np.exp(-0.5 * ((z_phase - 0.5) / 0.03) ** 2) * 0.25

    if fdict.get("ir_var_mode") == "lpv":
        ztf_flux *= 1.0 + 0.35 * np.sin(2.0 * np.pi * (ztf_time - 58200.0) / 420.0)

    # Under-luminous offset (true dimming vs dust/subdwarf)
    if fdict.get("true_dimming_mag"):
        dim = float(fdict["true_dimming_mag"])
        mags["G"] += dim
        abs_g += dim
        flux_mjy["G"] *= 10.0 ** (-0.4 * dim)

    # Optical spectrum (6000-6800 A, 0.25 A grid -> 3200 pixels)
    wave_a = np.linspace(6500.0, 6650.0, 600)
    # Photospheric continuum + H-alpha absorption at 6562.8 A
    cont = 1.0 - 0.30 * np.exp(-0.5 * ((wave_a - 6562.8) / 1.8) ** 2)
    spec_exp1 = cont + rng.normal(0.0, 0.012, len(wave_a))
    spec_exp2 = cont + rng.normal(0.0, 0.012, len(wave_a))
    sky_spec = np.full(len(wave_a), 0.05) + rng.normal(0.0, 0.005, len(wave_a))

    if fdict.get("active_halpha"):
        em = 0.65 * np.exp(-0.5 * ((wave_a - 6562.8) / 0.35) ** 2)
        spec_exp1 += em
        spec_exp2 += em
    if fdict.get("cosmic_ray"):
        # Unresolved narrow spike in exposure 1 ONLY (width < LSF)
        idx_cr = 260
        spec_exp1[idx_cr] += 1.2
    if fdict.get("airglow_line"):
        # Present in both target exposures AND sky spectrum at 6577.2 A OH sky line
        oh = 0.45 * np.exp(-0.5 * ((wave_a - 6577.2) / 0.25) ** 2)
        spec_exp1 += oh
        spec_exp2 += oh
        sky_spec += oh * 1.4
    if fdict.get("unresolved_unknown_line"):
        em = 0.50 * np.exp(-0.5 * ((wave_a - 6612.4) / 0.25) ** 2)
        spec_exp1 += em
        spec_exp2 += em

    # USNO-B1.0 / modern epoch match data for vanishing channel
    usno_epoch = 1953.4
    modern_epoch = 2016.0
    usno_ra, usno_dec = propagate_position(ra_deg, dec_deg, pmra, pmdec, modern_epoch, usno_epoch)
    usno_mag_r = mags["G"] + 0.2
    panstarrs_depth_mag = 21.8

    return {
        "source_id": source_id,
        "name": name or f"Gaia DR3 {source_id}",
        "gaia_release": gaia_release,
        "ra_deg": float(ra_deg),
        "dec_deg": float(dec_deg),
        "gal_b_deg": float(gal_b),
        "distance_pc": float(distance_pc),
        "parallax_mas": float(parallax_mas),
        "parallax_err_mas": float(parallax_err),
        "parallax_snr": float(parallax_snr),
        "pmra_mas_yr": float(pmra),
        "pmdec_mas_yr": float(pmdec),
        "ruwe": float(ruwe),
        "teff_k": float(teff_k),
        "teff_spec_k": float(fdict.get("teff_spec_k", teff_k + rng.normal(0.0, 40.0))),
        "logg": float(logg),
        "feh": float(feh),
        "luminosity_lsun": float(luminosity_lsun),
        "abs_g_mag": float(abs_g),
        "ebv_sfd": float(ebv_sfd),
        "ebv_planck": float(fdict.get("ebv_planck", ebv_sfd + rng.normal(0.0, 0.005))),
        "mags": mags,
        "mag_errs": mag_errs,
        "flux_mjy": flux_mjy,
        "flux_err_mjy": flux_err_mjy,
        "allwise": {
            "ra_deg": float(fdict.get("wise_ra_deg", wise_ra)),
            "dec_deg": float(fdict.get("wise_dec_deg", wise_dec)),
            "cc_flags": str(fdict.get("cc_flags", "0000")),
            "ext_flg": int(fdict.get("ext_flg", 0)),
            "ph_qual": str(fdict.get("ph_qual", "AAAA")),
            "w4rchi2": float(fdict.get("w4rchi2", 1.05)),
            "w4m": int(fdict.get("w4m", 12)),
            "n_neighbours_gaia_xmatch": int(fdict.get("n_neighbours_gaia_xmatch", 1)),
            "centroid_offset_arcsec": centroid_offset_arcsec,
            "cirrus_bg_rms_mjy": float(fdict.get("cirrus_bg_rms_mjy", 0.45)),
        },
        "wise_single_frames_w4_mjy": w4_frames,
        "neowise": {
            "mjd": [float(x) for x in visit_mjds],
            "w1_mjy": [float(x) for x in w1_visits],
            "w2_mjy": [float(x) for x in w2_visits],
            "w1_err_mjy": [float(flux_err_mjy["W1"]) for _ in range(20)],
            "w2_err_mjy": [float(flux_err_mjy["W2"]) for _ in range(20)],
            "scan_systematic": bool(fdict.get("scan_systematic", False)),
            "w1_w2_neighbour_blend": bool(fdict.get("w1_w2_neighbour_blend", False)),
        },
        "lightcurves": {
            "tess_time_btjd": [float(x) for x in tess_time],
            "tess_flux_qlp": [float(x) for x in tess_flux_a],
            "tess_flux_tesscut": [float(x) for x in tess_flux_b],
            "ztf_time_mjd": [float(x) for x in ztf_time],
            "ztf_flux": [float(x) for x in ztf_flux],
            "tess_scattered_light_window": bool(fdict.get("tess_scattered_light_window", False)),
            "neighbour_eb_centroid_shift": bool(fdict.get("neighbour_eb_centroid_shift", False)),
            "sso_crossing_t0": fdict.get("sso_crossing_t0"),
        },
        "vetting_archive": {
            "radio_sources_within_10arcsec": fdict.get("radio_sources_within_10arcsec", []),
            "deep_imaging_offset_galaxy": fdict.get("deep_imaging_offset_galaxy"),
            "deep_imaging_available": bool(fdict.get("deep_imaging_available", True)),
            "neighbours_in_w4_beam": fdict.get("neighbours_in_w4_beam", []),
            "local_red_source_density_per_arcsec2": float(fdict.get("local_red_density", 8e-5)),
            "moving_objects_2010": fdict.get("moving_objects_2010", []),
            "second_instrument_ir": fdict.get("second_instrument_ir"),
            "galex_uv_excess_mag": float(fdict.get("galex_uv_excess_mag", 0.0)),
            "xray_luminosity_erg_s": float(fdict.get("xray_luminosity_erg_s", 0.0)),
            "young_association": fdict.get("young_association"),
            "rv_scatter_kms": float(fdict.get("rv_scatter_kms", 0.25)),
            "nss_flag": int(fdict.get("nss_flag", 0)),
        },
        "spectra": {
            "wavelength_a": [float(x) for x in wave_a],
            "flux_exp1": [float(x) for x in spec_exp1],
            "flux_exp2": [float(x) for x in spec_exp2],
            "sky_flux": [float(x) for x in sky_spec],
            "lsf_fwhm_a": 0.45,
            "detector_artifact": bool(fdict.get("detector_artifact", False)),
        },
        "vanishing_data": {
            "usno_epoch": usno_epoch,
            "usno_ra_deg": float(fdict.get("usno_ra_override", usno_ra)),
            "usno_dec_deg": float(fdict.get("usno_dec_override", usno_dec)),
            "usno_mag_r": float(fdict.get("usno_mag_r", usno_mag_r)),
            "early_plates_count": int(fdict.get("early_plates_count", 2)),
            "present_in_panstarrs": bool(fdict.get("present_in_panstarrs", True)),
            "present_in_gaia": bool(fdict.get("present_in_gaia", True)),
            "panstarrs_depth_mag": float(fdict.get("panstarrs_depth_mag", panstarrs_depth_mag)),
            "sso_on_plate": bool(fdict.get("sso_on_plate", False)),
            "blend_split_sources": int(fdict.get("blend_split_sources", 1)),
            "known_lpv": bool(fdict.get("known_lpv", False)),
        },
        "simbad_type": str(fdict.get("simbad_type", simbad_type)),
        "known_in_literature": known_in_literature,
        "2mass_clean": bool(fdict.get("2mass_clean", True)),
        # Synthetic record with no real sky pixels: vetting reads its explicit inputs as evaluated values and
        # only these records may use the synthetic cutout renderer.
        "synthetic_fixture": True,
    }


def _build_named_objects(gaia_release: str = "dr3") -> dict[str, dict[str, Any]]:
    """Construct exact records for all named benchmarks, impostors, and open cases."""
    objs: list[dict[str, Any]] = [
        # Open Case 1: Gaia DR3 4582145718700356096 (Section 1 & Section 12)
        _make_base_star(
            "4582145718700356096",
            name="Gaia DR3 4582145718700356096",
            distance_pc=210.0,
            teff_k=6400.0,
            luminosity_lsun=1.12,
            ruwe=1.01,
            pmra=6.4,
            pmdec=-11.2,
            ebv_sfd=0.03,
            t_dust_k=226.0,
            frac_lum=0.0132,
            w4_extra_mult=1.08,  # gives ~59-60x W4 flux ratio
            gaia_release=gaia_release,
            simbad_type="Star",
            known_in_literature=False,
            flags={
                # Intentionally missing deep sub-arcsec imaging & post-2010 22um archive so both become not_evaluated -> Tier C
                "deep_imaging_available": False,
                "second_instrument_ir": None,
                "centroid_offset_arcsec": 0.25,
                "local_red_density": 1.1e-4,
            },
        ),
        # Open Case 2: HD 139139 (Random Transiter, Rappaport et al. 2019 / Alonso et al. 2023)
        _make_base_star(
            "6256050602424290048",
            name="HD 139139",
            distance_pc=107.0,
            teff_k=5760.0,
            luminosity_lsun=1.1,
            ruwe=1.04,
            gaia_release=gaia_release,
            simbad_type="Star",
            known_in_literature=True,
            flags={
                "optical_dips": {"depth": 0.05, "asymmetry": 0.1, "t0": 2010.0, "duration_days": 0.9, "in_ztf": False},
                "deep_imaging_available": True,
                "second_instrument_ir": None,
            },
        ),
        # Melis et al. 2012 TYC 8241-2652-1: strong earlier IRAS/AKARI excess that faded by WISE 2010
        _make_base_star(
            "5391513322976037248",
            name="TYC 8241-2652-1",
            distance_pc=120.0,
            teff_k=4520.0,
            luminosity_lsun=0.45,
            t_dust_k=450.0,
            frac_lum=0.0004,
            gaia_release=gaia_release,
            simbad_type="BYDra",
            known_in_literature=True,
            flags={
                "second_instrument_ir": {
                    "instrument": "AKARI_IRC_2006",
                    "wavelength_um": 18.0,
                    "flux_mjy": 140.0,
                    "wise_2010_flux_mjy": 4.2,
                    "disagreement_ratio": 33.3,
                    "interpretation": "real_secular_fading",
                },
            },
        ),
        # --- IR Excess Benchmarks (13 named stars from Section 7) ---
        _make_base_star("115725517896200320", name="BD+20 307", distance_pc=118.0, teff_k=6000.0, luminosity_lsun=1.5, t_dust_k=480.0, frac_lum=0.038, gaia_release=gaia_release, simbad_type="DebrisDisk", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "flux_mjy": 450.0, "confirmed": True}, "deep_imaging_available": True}),
        _make_base_star("5290837497970360960", name="ID8", distance_pc=195.0, teff_k=5400.0, luminosity_lsun=0.85, t_dust_k=550.0, frac_lum=0.025, gaia_release=gaia_release, simbad_type="YSO_Debris", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "flux_mjy": 180.0, "confirmed": True}, "ir_var_amp": 0.45, "ir_var_mode": "burst"}),
        _make_base_star("66549657436283904", name="HD 23514", distance_pc=138.0, teff_k=6600.0, luminosity_lsun=2.8, t_dust_k=700.0, frac_lum=0.020, gaia_release=gaia_release, simbad_type="DebrisDisk", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_IRS", "flux_mjy": 290.0, "confirmed": True}}),
        _make_base_star("1638028353993041792", name="TYC 4209-1322-1", distance_pc=145.0, teff_k=5100.0, luminosity_lsun=0.65, t_dust_k=420.0, frac_lum=0.032, gaia_release=gaia_release, simbad_type="ExtremeDebrisDisk", known_in_literature=True, flags={"ir_var_amp": 0.40, "ir_var_mode": "ramp", "second_instrument_ir": {"instrument": "AKARI_IRC", "confirmed": True}}),
        _make_base_star("6133409280036818688", name="HD 113766", distance_pc=131.0, teff_k=6800.0, luminosity_lsun=4.4, t_dust_k=440.0, frac_lum=0.021, gaia_release=gaia_release, simbad_type="DebrisDisk", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        _make_base_star("6438274350302427904", name="HD 172555", distance_pc=28.5, teff_k=7800.0, luminosity_lsun=7.8, t_dust_k=320.0, frac_lum=0.008, gaia_release=gaia_release, simbad_type="DebrisDisk", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        _make_base_star("5721914390814200960", name="HD 69830", distance_pc=12.6, teff_k=5385.0, luminosity_lsun=0.60, t_dust_k=340.0, frac_lum=0.0012, gaia_release=gaia_release, simbad_type="PlanetHost_Debris", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        _make_base_star("6485732215620751360", name="TYC 8830-410-1", distance_pc=162.0, teff_k=5400.0, luminosity_lsun=0.75, t_dust_k=520.0, frac_lum=0.018, gaia_release=gaia_release, simbad_type="ExtremeDebrisDisk", known_in_literature=True, flags={"ir_var_amp": 0.42, "ir_var_mode": "step", "second_instrument_ir": {"instrument": "AKARI_IRC", "confirmed": True}}),
        _make_base_star("2996362101189875200", name="TYC 5940-1510-1", distance_pc=175.0, teff_k=5900.0, luminosity_lsun=1.2, t_dust_k=460.0, frac_lum=0.016, gaia_release=gaia_release, simbad_type="ExtremeDebrisDisk", known_in_literature=True, flags={"ir_var_amp": 0.38, "ir_var_mode": "ramp", "second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        _make_base_star("3687158017250050176", name="TYC 4946-1106-1", distance_pc=152.0, teff_k=5650.0, luminosity_lsun=0.95, t_dust_k=410.0, frac_lum=0.014, gaia_release=gaia_release, simbad_type="ExtremeDebrisDisk", known_in_literature=True, flags={"ir_var_amp": 0.36, "ir_var_mode": "burst", "second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        _make_base_star("5551530885782831616", name="TYC 8105-370-1", distance_pc=180.0, teff_k=6100.0, luminosity_lsun=1.4, t_dust_k=390.0, frac_lum=0.012, gaia_release=gaia_release, simbad_type="ExtremeDebrisDisk", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "AKARI_IRC", "confirmed": True}}),
        _make_base_star("4066368953121417984", name="HD 166191", distance_pc=101.0, teff_k=6150.0, luminosity_lsun=4.2, t_dust_k=580.0, frac_lum=0.045, gaia_release=gaia_release, simbad_type="ExtremeDebrisDisk", known_in_literature=True, flags={"ir_var_amp": 0.52, "ir_var_mode": "burst", "second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        _make_base_star("6048418276566039424", name="HD 145263", distance_pc=142.0, teff_k=7100.0, luminosity_lsun=5.5, t_dust_k=310.0, frac_lum=0.015, gaia_release=gaia_release, simbad_type="DebrisDisk", known_in_literature=True, flags={"second_instrument_ir": {"instrument": "Spitzer_MIPS_24", "confirmed": True}}),
        # --- IR Time-domain + Optical Dips ---
        _make_base_star("5548511859528971520", name="ASASSN-21qj", distance_pc=560.0, teff_k=5750.0, luminosity_lsun=1.0, t_dust_k=720.0, frac_lum=0.040, gaia_release=gaia_release, simbad_type="CollisionDebris", known_in_literature=True, flags={"ir_var_amp": 0.65, "ir_var_mode": "burst", "optical_dips": {"depth": 0.35, "asymmetry": 0.45, "t0": 2014.0, "duration_days": 3.5, "in_ztf": True}, "second_instrument_ir": {"instrument": "NEOWISE", "confirmed": True}}),
        _make_base_star("2081900940499099136", name="KIC 8462852", distance_pc=451.0, teff_k=6750.0, luminosity_lsun=4.7, gaia_release=gaia_release, simbad_type="TabbysStar", known_in_literature=True, flags={"optical_dips": {"depth": 0.22, "asymmetry": 0.52, "t0": 2011.5, "duration_days": 2.4, "in_ztf": True}}),
        _make_base_star("6049824505359595776", name="EPIC 204278916", distance_pc=145.0, teff_k=3900.0, luminosity_lsun=0.15, gaia_release=gaia_release, simbad_type="Dipper", known_in_literature=True, flags={"optical_dips": {"depth": 0.65, "asymmetry": 0.48, "t0": 2013.0, "duration_days": 2.1, "in_ztf": True}}),
        # --- Impostors: IR Excess ---
        _make_base_star("4418120359123456128", name="Hephaistos-A", distance_pc=160.0, teff_k=3600.0, luminosity_lsun=0.08, t_dust_k=220.0, frac_lum=0.025, gaia_release=gaia_release, flags={"centroid_offset_arcsec": 4.9, "radio_sources_within_10arcsec": [{"survey": "VLASS", "sep_arcsec": 4.9, "flux_mjy": 1.8}], "deep_imaging_offset_galaxy": {"sep_arcsec": 4.9, "r_mag": 21.3, "ellipticity": 0.45}}),
        _make_base_star("1394205817392045184", name="Hephaistos-B", distance_pc=185.0, teff_k=3500.0, luminosity_lsun=0.06, t_dust_k=240.0, frac_lum=0.030, gaia_release=gaia_release, flags={"centroid_offset_arcsec": 0.4, "radio_sources_within_10arcsec": [{"survey": "FIRST", "sep_arcsec": 0.4, "flux_mjy": 2.4}], "deep_imaging_offset_galaxy": {"sep_arcsec": 0.4, "r_mag": 22.1, "ellipticity": 0.38}}),
        _make_base_star("2059182736491028480", name="Hephaistos-G", distance_pc=175.0, teff_k=3700.0, luminosity_lsun=0.09, t_dust_k=210.0, frac_lum=0.028, gaia_release=gaia_release, flags={"centroid_offset_arcsec": 5.0, "deep_imaging_offset_galaxy": {"sep_arcsec": 5.0, "r_mag": 21.8, "ellipticity": 0.51}, "neighbours_in_w4_beam": [{"id": "LS-DR10-BG", "sep_arcsec": 5.0, "flux_frac": 0.92}]}),
        _make_base_star("3141592653589793280", name="DOG-Blend-1", distance_pc=140.0, teff_k=4800.0, luminosity_lsun=0.3, t_dust_k=190.0, frac_lum=0.04, gaia_release=gaia_release, flags={"centroid_offset_arcsec": 6.2, "neighbours_in_w4_beam": [{"id": "WISE-DOG", "sep_arcsec": 6.2, "flux_frac": 0.95}], "deep_imaging_offset_galaxy": {"sep_arcsec": 6.2, "r_mag": 23.0, "ellipticity": 0.42}}),
        _make_base_star("5812938471029384960", name="Cirrus-Impostor-1", distance_pc=210.0, teff_k=5500.0, ebv_sfd=0.45, ebv_planck=0.46, t_dust_k=140.0, frac_lum=0.05, gaia_release=gaia_release, flags={"cirrus_bg_rms_mjy": 14.5, "ext_flg": 5, "w4rchi2": 8.4}),
        _make_base_star("4129384756102938112", name="WISE-Artifact-1", distance_pc=115.0, teff_k=5900.0, t_dust_k=220.0, frac_lum=0.02, gaia_release=gaia_release, flags={"cc_flags": "000D", "ph_qual": "AAAU", "w4rchi2": 11.2, "single_frame_glitch": True}),
        # --- Impostors: IR Time-Domain ---
        _make_base_star("711100000000000001", name="TD-Blend-1", distance_pc=190.0, teff_k=5100.0, gaia_release=gaia_release, flags={"ir_var_amp": 0.45, "ir_var_mode": "step", "w1_w2_neighbour_blend": True, "neighbours_in_w4_beam": [{"id": "VarNeighbour", "sep_arcsec": 4.2, "flux_frac": 0.65}]}),
        _make_base_star("711100000000000002", name="TD-LPV-1", distance_pc=220.0, teff_k=3400.0, luminosity_lsun=120.0, gaia_release=gaia_release, simbad_type="LPV*", flags={"ir_var_amp": 0.55, "ir_var_mode": "lpv"}),
        _make_base_star("711100000000000003", name="TD-Scan-1", distance_pc=160.0, teff_k=5700.0, gaia_release=gaia_release, flags={"ir_var_amp": 0.42, "ir_var_mode": "step", "scan_systematic": True}),
        # --- Impostors: Optical Dips ---
        _make_base_star("722200000000000001", name="KEB-Impostor-1", distance_pc=240.0, teff_k=5900.0, gaia_release=gaia_release, simbad_type="EB*", flags={"periodic_eb": True, "eb_period_days": 2.4}),
        _make_base_star("722200000000000002", name="TESS-Scatter-1", distance_pc=180.0, teff_k=5400.0, gaia_release=gaia_release, flags={"optical_dips": {"depth": 0.18, "asymmetry": 0.3, "t0": 2011.0, "duration_days": 1.2, "in_ztf": False}, "detrend_artifact_only": True, "tess_scattered_light_window": True}),
        _make_base_star("722200000000000003", name="Blend-EB-1", distance_pc=210.0, teff_k=6000.0, gaia_release=gaia_release, flags={"optical_dips": {"depth": 0.14, "asymmetry": 0.2, "t0": 2013.0, "duration_days": 1.0, "in_ztf": True}, "neighbour_eb_centroid_shift": True}),
        # --- Under-luminous switch-off & impostors ---
        _make_base_star("733300000000000001", name="UL-Ext-1", distance_pc=250.0, teff_k=5800.0, ebv_sfd=0.72, ebv_planck=0.71, gaia_release=gaia_release),
        _make_base_star("733300000000000002", name="UL-Ext-2", distance_pc=270.0, teff_k=6200.0, ebv_sfd=0.80, ebv_planck=0.81, gaia_release=gaia_release),
        _make_base_star("733300000000000003", name="UL-Ext-3", distance_pc=220.0, teff_k=5200.0, ebv_sfd=0.68, ebv_planck=0.69, gaia_release=gaia_release),
        _make_base_star("733300000000000004", name="UL-Subdwarf-1", distance_pc=140.0, teff_k=5600.0, feh=-1.85, luminosity_lsun=0.40, gaia_release=gaia_release, simbad_type="Subdwarf"),
        _make_base_star("733300000000000005", name="UL-HighRUWE-1", distance_pc=175.0, teff_k=5900.0, ruwe=2.95, luminosity_lsun=0.45, gaia_release=gaia_release, flags={"true_dimming_mag": 0.9, "nss_flag": 1}),
        # --- Spectral lines switch-off & impostors ---
        _make_base_star("744400000000000001", name="SL-ActiveMDwarf-1", distance_pc=35.0, teff_k=3400.0, luminosity_lsun=0.04, gaia_release=gaia_release, simbad_type="Flare*", flags={"active_halpha": True, "xray_luminosity_erg_s": 2e29, "galex_uv_excess_mag": 2.2}),
        _make_base_star("744400000000000002", name="SL-ActiveMDwarf-2", distance_pc=42.0, teff_k=3300.0, luminosity_lsun=0.03, gaia_release=gaia_release, simbad_type="Flare*", flags={"active_halpha": True, "xray_luminosity_erg_s": 1.5e29, "galex_uv_excess_mag": 1.9}),
        _make_base_star("744400000000000003", name="SL-CosmicRay-1", distance_pc=95.0, teff_k=5500.0, gaia_release=gaia_release, flags={"cosmic_ray": True}),
        _make_base_star("744400000000000004", name="SL-Airglow-1", distance_pc=110.0, teff_k=5700.0, gaia_release=gaia_release, flags={"airglow_line": True}),
        # --- Vanishing switch-off, positives & impostors ---
        _make_base_star("755500000000000001", name="VAN-HighPM-1", distance_pc=18.0, teff_k=4100.0, pmra=185.0, pmdec=-140.0, gaia_release=gaia_release),
        _make_base_star("755500000000000002", name="VAN-HighPM-2", distance_pc=22.0, teff_k=3900.0, pmra=-210.0, pmdec=95.0, gaia_release=gaia_release),
        _make_base_star("755500000000000003", name="VAN-Villarroel-1", distance_pc=260.0, teff_k=5000.0, gaia_release=gaia_release, flags={"early_plates_count": 2, "present_in_panstarrs": False, "present_in_gaia": False, "usno_mag_r": 16.5, "panstarrs_depth_mag": 22.0}),
        _make_base_star("755500000000000004", name="VAN-Villarroel-2", distance_pc=280.0, teff_k=4700.0, gaia_release=gaia_release, flags={"early_plates_count": 2, "present_in_panstarrs": False, "present_in_gaia": False, "usno_mag_r": 17.0, "panstarrs_depth_mag": 21.8}),
        _make_base_star("755500000000000005", name="VAN-LPV-1", distance_pc=240.0, teff_k=3300.0, luminosity_lsun=80.0, gaia_release=gaia_release, simbad_type="Mira", flags={"early_plates_count": 2, "present_in_panstarrs": False, "present_in_gaia": False, "usno_mag_r": 14.5, "panstarrs_depth_mag": 21.5, "known_lpv": True}),
        _make_base_star("755500000000000006", name="VAN-PlateDefect-1", distance_pc=150.0, teff_k=5800.0, gaia_release=gaia_release, flags={"early_plates_count": 1, "present_in_panstarrs": False, "present_in_gaia": False, "usno_mag_r": 16.2, "panstarrs_depth_mag": 21.8}),
    ]
    return {o["source_id"]: o for o in objs}


def generate_field_sample(
    n_stars: int = 10000,
    seed: int = 20261116,
    gaia_release: str = "dr3",
) -> list[dict[str, Any]]:
    """Generate a deterministic, realistic N-star parent sample matching Section 5.1 cuts.

    Photometric fluxes follow Planck + realistic calibration scatter with exact Gaussian errors
    so that the chi distribution of ir_excess on ordinary stars is centered near 0 with unit width.
    """
    rng = np.random.default_rng(seed)
    stars: list[dict[str, Any]] = []

    teffs = rng.uniform(3500.0, 8800.0, n_stars)
    distances = rng.uniform(25.0, 295.0, n_stars)
    ras = rng.uniform(0.0, 360.0, n_stars)
    # Sample outside Galactic plane (|b| > 10 deg) for most, plus 6% inside plane to test sample cuts
    decs = rng.uniform(-78.0, 78.0, n_stars)
    ruwes = np.where(rng.random(n_stars) < 0.04, rng.uniform(1.45, 3.2, n_stars), rng.uniform(0.88, 1.25, n_stars))
    plx_snrs = np.where(rng.random(n_stars) < 0.03, rng.uniform(4.0, 9.5, n_stars), rng.uniform(15.0, 180.0, n_stars))
    if gaia_release == "dr4":
        plx_snrs *= 1.41

    rel_err = {
        "G": 0.004,
        "BP": 0.008,
        "RP": 0.008,
        "J": 0.020,
        "H": 0.022,
        "Ks": 0.024,
        "W1": 0.025,
        "W2": 0.028,
        "W3": 0.045,
        "W4": 0.085,
    }

    for i in range(n_stars):
        sid = f"9{10000000000000000 + i}"
        ra = float(ras[i])
        dec = float(decs[i])
        _, gal_b = equatorial_to_galactic(ra, dec)
        dist = float(distances[i])
        teff = float(teffs[i])
        lum = max(0.02, (teff / 5778.0) ** 4.8)
        star_scale_mjy = 32500.0 * lum * ((10.0 / dist) ** 2) * ((5778.0 / teff) ** 2.2)
        ebv = float(abs(rng.normal(0.02, 0.01)))
        av = 3.1 * ebv

        flux_mjy: dict[str, float] = {}
        flux_err_mjy: dict[str, float] = {}
        mags: dict[str, float] = {}
        mag_errs: dict[str, float] = {}
        for band, lam in BAND_WAVELENGTH_UM.items():
            f_true = stellar_photosphere_flux_mjy(lam, teff, star_scale_mjy) * (
                10.0 ** (-0.4 * av * BAND_EXTINCTION_AV[band])
            )
            sig = max(0.02, f_true * rel_err[band])
            # Include photometric noise + realistic calibration/photosphere scatter in WISE bands
            extra_cal = 0.052 * f_true if band in {"W1", "W2", "W3", "W4"} else 0.0
            total_sig = math.hypot(sig, extra_cal)
            f_obs = max(0.02, float(f_true + rng.normal(0.0, total_sig)))
            flux_mjy[band] = f_obs
            flux_err_mjy[band] = sig
            mags[band] = flux_mjy_to_mag(band, f_obs)
            mag_errs[band] = min(0.35, 1.0857 * sig / f_obs)

        abs_g = mags["G"] - 5.0 * math.log10(dist / 10.0) - av * BAND_EXTINCTION_AV["G"]
        plx_mas = 1000.0 / dist
        stars.append(
            {
                "source_id": sid,
                "name": f"Gaia DR3 {sid}",
                "gaia_release": gaia_release,
                "ra_deg": ra,
                "dec_deg": dec,
                "gal_b_deg": gal_b,
                "distance_pc": dist,
                "parallax_mas": plx_mas,
                "parallax_err_mas": plx_mas / float(plx_snrs[i]),
                "parallax_snr": float(plx_snrs[i]),
                "pmra_mas_yr": float(rng.normal(0.0, 12.0)),
                "pmdec_mas_yr": float(rng.normal(0.0, 12.0)),
                "ruwe": float(ruwes[i]),
                "teff_k": teff,
                "teff_spec_k": float(teff + rng.normal(0.0, 45.0)),
                "logg": 4.4,
                "feh": float(rng.normal(0.0, 0.15)),
                "luminosity_lsun": float(lum),
                "abs_g_mag": float(abs_g),
                "ebv_sfd": ebv,
                "ebv_planck": float(ebv + rng.normal(0.0, 0.004)),
                "mags": mags,
                "mag_errs": mag_errs,
                "flux_mjy": flux_mjy,
                "flux_err_mjy": flux_err_mjy,
                "allwise": {
                    "ra_deg": ra,
                    "dec_deg": dec,
                    "cc_flags": "0000",
                    "ext_flg": 0,
                    "ph_qual": "AAAA",
                    "w4rchi2": 1.0,
                    "w4m": 12,
                    "n_neighbours_gaia_xmatch": 1,
                    "centroid_offset_arcsec": 0.15,
                    "cirrus_bg_rms_mjy": 0.4,
                },
                "2mass_clean": True,
                "simbad_type": "Star",
                "known_in_literature": False,
                "synthetic_fixture": True,
            }
        )
    return stars


def build_radio_cadences() -> dict[str, dict[str, Any]]:
    """Build realistic Breakthrough Listen ON/OFF filterbank cadences (Voyager 1 + blc1 + controls)."""
    # Frequency axis: 256 channels x 12 Hz = 3072 Hz span (+/- 1536 Hz)
    # 3 ON scans alternating with 3 OFF scans (ABACAD cadence), 60 s per scan
    rng = np.random.default_rng(4242)
    n_time, n_chan = 16, 256
    dt_s = 3.75  # 60s scan = 16 * 3.75s
    df_hz = 12.0

    def _make_cadence(
        target_id: str,
        *,
        f0_hz: float,
        drift_hz_s: float,
        snr: float,
        in_on: bool,
        in_off: bool,
        in_rfi_band: bool,
        multi_target: bool,
        reobserved: bool,
    ) -> dict[str, Any]:
        scans = []
        for idx, kind in enumerate(["ON", "OFF", "ON", "OFF", "ON", "OFF"]):
            arr = rng.normal(1.0, 0.08, (n_time, n_chan))
            t_scan_offset = idx * 60.0
            active = (kind == "ON" and in_on) or (kind == "OFF" and in_off)
            if active:
                for ti in range(n_time):
                    t_total = t_scan_offset + ti * dt_s
                    f_hz = f0_hz + drift_hz_s * t_total
                    ch = int(round(n_chan / 2 + f_hz / df_hz))
                    if 0 <= ch < n_chan:
                        arr[ti, ch] += snr * 0.08
            scans.append({"scan_index": idx, "kind": kind, "t_start_s": t_scan_offset, "power": arr.tolist()})
        return {
            "source_id": target_id,
            "dt_s": dt_s,
            "df_hz": df_hz,
            "center_freq_mhz": 8419.5 if target_id == "GBT_VOYAGER1_20151230" else 982.002,
            "scans": scans,
            "in_rfi_band": in_rfi_band,
            "multi_target_match": multi_target,
            "reobserved_confirmed": reobserved,
        }

    cadences: dict[str, dict[str, Any]] = {
        "GBT_VOYAGER1_20151230": _make_cadence(
            "GBT_VOYAGER1_20151230",
            f0_hz=-120.0,
            drift_hz_s=-0.38,  # Earth rotation Doppler drift in X-band
            snr=32.0,
            in_on=True,
            in_off=False,
            in_rfi_band=False,
            multi_target=False,
            reobserved=True,
        ),
        "GBT_BLC1_PROXCEN": _make_cadence(
            "GBT_BLC1_PROXCEN",
            f0_hz=45.0,
            drift_hz_s=0.06,
            snr=18.0,
            in_on=True,
            in_off=True,  # Intermodulation RFI also seen off-target / in known allocation & multi-target
            in_rfi_band=True,
            multi_target=True,
            reobserved=False,
        ),
    }
    # Add 24 ordinary quiet BL control cadences
    for i in range(24):
        cid = f"GBT_CONTROL_{i:03d}"
        cadences[cid] = _make_cadence(
            cid,
            f0_hz=0.0,
            drift_hz_s=0.0,
            snr=0.0,
            in_on=False,
            in_off=False,
            in_rfi_band=False,
            multi_target=False,
            reobserved=False,
        )
    return cadences


NAMED_OBJECTS_DR3 = _build_named_objects("dr3")
NAMED_OBJECTS_DR4 = _build_named_objects("dr4")
RADIO_CADENCES = build_radio_cadences()


def get_named_object(source_id: str, gaia_release: str = "dr3") -> dict[str, Any] | None:
    clean_id = source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
    bank = NAMED_OBJECTS_DR4 if gaia_release == "dr4" else NAMED_OBJECTS_DR3
    return bank.get(clean_id)
