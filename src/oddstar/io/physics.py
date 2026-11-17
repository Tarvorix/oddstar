"""Astrophysical radiative, photometric, and astrometric utilities."""

from __future__ import annotations

import math
from typing import Final

# Fundamental physical constants (CODATA SI)
H_PLANCK: Final[float] = 6.62607015e-34
C_LIGHT: Final[float] = 2.99792458e8
K_BOLTZ: Final[float] = 1.380649e-23
HC_OVER_K_UM_K: Final[float] = (H_PLANCK * C_LIGHT / K_BOLTZ) * 1e6  # ~14387.7688 um*K

# Effective wavelengths (microns) and Vega zero-point flux densities (Jy)
BAND_WAVELENGTH_UM: Final[dict[str, float]] = {
    "G": 0.622,
    "BP": 0.511,
    "RP": 0.777,
    "J": 1.235,
    "H": 1.662,
    "Ks": 2.159,
    "W1": 3.353,
    "W2": 4.603,
    "W3": 11.561,
    "W4": 22.088,
}

BAND_ZERO_POINT_JY: Final[dict[str, float]] = {
    "G": 3229.0,
    "BP": 3552.0,
    "RP": 2555.0,
    "J": 1594.0,
    "H": 1024.0,
    "Ks": 666.7,
    "W1": 309.54,
    "W2": 171.787,
    "W3": 31.674,
    "W4": 8.363,
}

# Extinction coefficients A_lambda / A_V (Cardelli/Wang-Chen IR-aware curve)
BAND_EXTINCTION_AV: Final[dict[str, float]] = {
    "G": 0.789,
    "BP": 1.002,
    "RP": 0.589,
    "J": 0.282,
    "H": 0.175,
    "Ks": 0.112,
    "W1": 0.056,
    "W2": 0.043,
    "W3": 0.032,
    "W4": 0.019,
}


def planck_bnu_relative(wavelength_um: float, temp_k: float) -> float:
    """Relative specific intensity B_nu(lambda, T) normalized in mJy-scale units."""
    x = HC_OVER_K_UM_K / max(1e-6, wavelength_um * max(1.0, temp_k))
    if x > 700.0:
        return 0.0
    # 2 h c / lambda^3 in arbitrary unit proportional to Jy/sr: (1 / um^3) / (exp(x) - 1)
    return (1.0 / (wavelength_um**3)) / math.expm1(x)


def stellar_photosphere_flux_mjy(
    wavelength_um: float,
    teff_k: float,
    scale_mjy: float,
) -> float:
    """Compute photospheric flux density (mJy) normalized at Ks-scale factor."""
    ref = planck_bnu_relative(BAND_WAVELENGTH_UM["Ks"], teff_k)
    if ref <= 0.0:
        return 0.0
    return scale_mjy * (planck_bnu_relative(wavelength_um, teff_k) / ref)


def dust_blackbody_flux_mjy(
    wavelength_um: float,
    teff_k: float,
    star_scale_mjy: float,
    t_dust_k: float,
    frac_lum: float,
) -> float:
    """Compute optically thin / partial-shell blackbody thermal emission flux density (mJy).

    Uses exact bolometric energy balance:
      F_nu,dust(lam) / F_nu,star(lam) = frac_lum * (T_star / T_dust)^4 * B_nu(lam, T_dust) / B_nu(lam, T_star).
    """
    f_star = stellar_photosphere_flux_mjy(wavelength_um, teff_k, star_scale_mjy)
    b_star = planck_bnu_relative(wavelength_um, teff_k)
    b_dust = planck_bnu_relative(wavelength_um, t_dust_k)
    if b_star <= 0.0 or b_dust <= 0.0 or frac_lum <= 0.0:
        return 0.0
    return f_star * frac_lum * ((teff_k / t_dust_k) ** 4) * (b_dust / b_star)


def mag_to_flux_mjy(band: str, mag: float) -> float:
    """Convert Vega magnitude to flux density in mJy."""
    zp = BAND_ZERO_POINT_JY[band]
    return zp * 1000.0 * (10.0 ** (-0.4 * mag))


def flux_mjy_to_mag(band: str, flux_mjy: float) -> float:
    """Convert flux density in mJy to Vega magnitude."""
    zp = BAND_ZERO_POINT_JY[band]
    val = max(1e-12, flux_mjy / (zp * 1000.0))
    return -2.5 * math.log10(val)


def mag_err_to_flux_err_mjy(flux_mjy: float, mag_err: float) -> float:
    """Propagate magnitude uncertainty into flux density uncertainty (mJy)."""
    return flux_mjy * (math.log(10.0) / 2.5) * mag_err


def equilibrium_radius_au(t_dust_k: float, luminosity_lsun: float) -> float:
    """Blackbody equilibrium radius r = (278 K / T)^2 * sqrt(L / L_sun) in au."""
    return ((278.0 / max(1.0, t_dust_k)) ** 2) * math.sqrt(max(1e-6, luminosity_lsun))


def propagate_position(
    ra_deg: float,
    dec_deg: float,
    pmra_mas_yr: float,
    pmdec_mas_yr: float,
    epoch_from: float,
    epoch_to: float,
) -> tuple[float, float]:
    """Propagate ICRS coordinates between epochs using Gaia proper motion (mas/yr)."""
    dt_yr = epoch_to - epoch_from
    cos_dec = max(1e-6, math.cos(math.radians(dec_deg)))
    dra_deg = (pmra_mas_yr * dt_yr) / (3600000.0 * cos_dec)
    ddec_deg = (pmdec_mas_yr * dt_yr) / 3600000.0
    return (ra_deg + dra_deg) % 360.0, max(-90.0, min(90.0, dec_deg + ddec_deg))


def angular_separation_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Exact Vincenty spherical angular separation in arcseconds."""
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    dr = r2 - r1
    num = math.hypot(
        math.cos(d2) * math.sin(dr),
        math.cos(d1) * math.sin(d2) - math.sin(d1) * math.cos(d2) * math.cos(dr),
    )
    den = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(dr)
    return math.degrees(math.atan2(num, den)) * 3600.0
