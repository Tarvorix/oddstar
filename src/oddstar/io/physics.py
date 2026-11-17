"""Astrophysical radiative, photometric, and astrometric utilities."""

from __future__ import annotations

import json
import math
from functools import lru_cache
from importlib import resources
from typing import Any, Final

import numpy as np

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


def planck_photosphere_flux_mjy(
    wavelength_um: float,
    teff_k: float,
    scale_mjy: float,
) -> float:
    """Pure Planck photosphere (mJy) normalised so that the flux at Ks equals `scale_mjy`."""
    ref = planck_bnu_relative(BAND_WAVELENGTH_UM["Ks"], teff_k)
    if ref <= 0.0:
        return 0.0
    return scale_mjy * (planck_bnu_relative(wavelength_um, teff_k) / ref)


# ---------------------------------------------------------------------------------------------------------------
# Empirical dwarf colour sequence (Pecaut & Mamajek 2013 online table; see the provenance header of the data file)
# ---------------------------------------------------------------------------------------------------------------

DWARF_COLOUR_TABLE_FILE: Final[str] = "pecaut_mamajek_dwarf_colours.json"
WISE_BAND_NAMES: Final[tuple[str, ...]] = ("W1", "W2", "W3", "W4")
R_V: Final[float] = 3.1


@lru_cache(maxsize=1)
def load_dwarf_colour_table() -> dict[str, Any]:
    """Load the packaged dwarf colour table (with its provenance header)."""
    txt = resources.files("oddstar.refdata").joinpath(DWARF_COLOUR_TABLE_FILE).read_text(encoding="utf-8")
    return json.loads(txt)


def literature_photosphere_label() -> str:
    """`photosphere_source` label for anything computed from the packaged literature table.

    Values typed from memory (provenance != "retrieved") are always labelled `literature_fallback_unverified`.
    """
    prov = str(load_dwarf_colour_table().get("provenance", ""))
    return "literature_fallback" if prov == "retrieved" else "literature_fallback_unverified"


@lru_cache(maxsize=8)
def dwarf_ks_minus_wise_sequence(band: str) -> tuple[list[float], list[float], list[float]]:
    """Return (teff_k ascending, bp_rp, Ks-band) rows of the dwarf table where Ks-band is tabulated.

    Ks-W1 is tabulated A0V-M5V. Ks-W2/W3/W4 = Ks-W1 + (W1-W2/W3/W4) exist only for the contiguous block of rows
    where the W1-Wx colour is tabulated (A0V-K5V in the 2022.04.16 table); isolated later entries are ignored so
    that no colour is interpolated across an untabulated gap.
    """
    rows = sorted(load_dwarf_colour_table()["rows"], key=lambda r: float(r["teff_k"]), reverse=True)
    key = {"W1": None, "W2": "w1_w2", "W3": "w1_w3", "W4": "w1_w4"}[band]
    teffs: list[float] = []
    bprp: list[float] = []
    vals: list[float] = []
    for r in rows:  # hottest -> coolest; stop at the first gap
        if r.get("ks_w1") is None or r.get("bp_rp") is None:
            break
        if key is None:
            v = float(r["ks_w1"])
        else:
            if r.get(key) is None:
                break
            v = float(r["ks_w1"]) + float(r[key])
        teffs.append(float(r["teff_k"]))
        bprp.append(float(r["bp_rp"]))
        vals.append(v)
    order = np.argsort(teffs)
    return [teffs[i] for i in order], [bprp[i] for i in order], [vals[i] for i in order]


def planck_ks_minus_band(band: str, teff_k: float) -> float:
    """Vega Ks - band colour of a pure Planck photosphere with this code's effective wavelengths and zero points."""
    b_ks = planck_bnu_relative(BAND_WAVELENGTH_UM["Ks"], teff_k)
    b_x = planck_bnu_relative(BAND_WAVELENGTH_UM[band], teff_k)
    ratio = (b_x / max(1e-300, b_ks)) * (BAND_ZERO_POINT_JY["Ks"] / BAND_ZERO_POINT_JY[band])
    return 2.5 * math.log10(max(1e-300, ratio))


@lru_cache(maxsize=8192)
def wise_colour_correction_mag(teff_k: float) -> tuple[float, float, float, float]:
    """Empirical minus Planck (Ks - Wx) colour for W1..W4 at this Teff (mag; positive = W brighter than Planck).

    Teff outside the tabulated range of a band is clamped to the nearest tabulated row.
    """
    out: list[float] = []
    for band in WISE_BAND_NAMES:
        teffs, _, vals = dwarf_ks_minus_wise_sequence(band)
        t = min(max(float(teff_k), teffs[0]), teffs[-1])
        emp = float(np.interp(t, teffs, vals))
        out.append(emp - planck_ks_minus_band(band, float(teff_k)))
    return out[0], out[1], out[2], out[3]


def _wise_correction_at_wavelength(wavelength_um: float, teff_k: float) -> float:
    lam_ks = BAND_WAVELENGTH_UM["Ks"]
    if wavelength_um <= lam_ks:
        return 0.0
    knots = [math.log(lam_ks)] + [math.log(BAND_WAVELENGTH_UM[b]) for b in WISE_BAND_NAMES]
    corr = [0.0, *wise_colour_correction_mag(round(float(teff_k), 1))]
    return float(np.interp(math.log(wavelength_um), knots, corr))


def stellar_photosphere_flux_mjy(
    wavelength_um: float,
    teff_k: float,
    scale_mjy: float,
) -> float:
    """Photospheric flux density (mJy) normalised at Ks: Planck shape plus an empirical WISE colour correction.

    Shortward of Ks this is the Planck curve. Longward of Ks the Planck curve is multiplied by 10^(0.4 dC(lambda)),
    where dC is (Ks - Wx)_empirical(Teff) - (Ks - Wx)_Planck(Teff) at the four WISE bands (Pecaut & Mamajek dwarf
    sequence, packaged in oddstar/refdata) interpolated in log wavelength, and held at its W4 value beyond 22 um.
    At each WISE band the model therefore reproduces the tabulated dwarf Ks - Wx colour for that Teff. Used both by
    the synthetic fixtures and by the secondary photosphere model of `ir_excess`. `planck_photosphere_flux_mjy`
    is the uncorrected Planck curve.
    """
    base = planck_photosphere_flux_mjy(wavelength_um, teff_k, scale_mjy)
    if base <= 0.0:
        return 0.0
    return base * 10.0 ** (0.4 * _wise_correction_at_wavelength(wavelength_um, teff_k))


def dust_blackbody_flux_mjy(
    wavelength_um: float,
    teff_k: float,
    star_scale_mjy: float,
    t_dust_k: float,
    frac_lum: float,
) -> float:
    """Compute optically thin / partial-shell blackbody thermal emission flux density (mJy).

    Uses exact bolometric energy balance against the Planck-shaped stellar SED (normalised at Ks):
      F_nu,dust(lam) = F_nu,Planck star(lam) * frac_lum * (T_star / T_dust)^4 * B_nu(lam, T_dust) / B_nu(lam, T_star).
    The empirical WISE colour correction of the photosphere is deliberately not applied here: it redistributes a
    negligible fraction of the stellar bolometric flux, and keeping the dust model independent of it means
    `frac_lum` keeps one meaning whichever photosphere predictor is used.
    """
    f_star = planck_photosphere_flux_mjy(wavelength_um, teff_k, star_scale_mjy)
    b_star = planck_bnu_relative(wavelength_um, teff_k)
    b_dust = planck_bnu_relative(wavelength_um, t_dust_k)
    if b_star <= 0.0 or b_dust <= 0.0 or frac_lum <= 0.0:
        return 0.0
    return f_star * frac_lum * ((teff_k / t_dust_k) ** 4) * (b_dust / b_star)


def blackbody_colour_temperature_k(
    band_short: str,
    band_long: str,
    flux_ratio_long_over_short: float,
    *,
    t_min_k: float = 30.0,
    t_max_k: float = 3000.0,
) -> tuple[float, str]:
    """Blackbody temperature whose F_nu(long)/F_nu(short) equals the given ratio.

    Returns (T, status) with status "ok", "above_max" (ratio at or below the T_max value, e.g. Rayleigh-Jeans
    limit) or "below_min".
    """
    l1 = BAND_WAVELENGTH_UM[band_short]
    l2 = BAND_WAVELENGTH_UM[band_long]

    def _ratio(t: float) -> float:
        return planck_bnu_relative(l2, t) / max(1e-300, planck_bnu_relative(l1, t))

    r = float(flux_ratio_long_over_short)
    if not math.isfinite(r) or r <= 0.0:
        return t_max_k, "above_max"
    # Ratio decreases monotonically with temperature
    if r <= _ratio(t_max_k):
        return t_max_k, "above_max"
    if r >= _ratio(t_min_k):
        return t_min_k, "below_min"
    lo, hi = math.log(t_min_k), math.log(t_max_k)
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if _ratio(math.exp(mid)) > r:
            lo = mid
        else:
            hi = mid
    return float(math.exp(0.5 * (lo + hi))), "ok"


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


# ---------------------------------------------------------------------------------------------------------------
# Extinction bookkeeping
#   ebv_star / ebv_star_source : reddening applied to the star (gaia_gspphot, map3d_<name>, assumed_zero, fixture,
#                                missing)
#   ebv_sfd / ebv_planck       : true line-of-sight (total column) map values only, or None
#   vetting_archive.dust       : line-of-sight values from vetting enrichment (ebv_sfd_los, ebv_sf11_los, i100)
#   vetting_archive.dust3d     : 3D-map reddening at the star's distance (oddstar.io.dust3d)
# ---------------------------------------------------------------------------------------------------------------


def extinction_mag(band: str, ebv: float | None) -> float:
    """A_band (mag) for E(B-V) with R_V = 3.1 and the IR-aware A_lambda/A_V coefficients above (0 if ebv is None)."""
    if ebv is None:
        return 0.0
    return R_V * float(ebv) * BAND_EXTINCTION_AV[band]


def colour_excess(band_a: str, band_b: str, ebv: float | None) -> float:
    """E(a - b) for E(B-V)."""
    return extinction_mag(band_a, ebv) - extinction_mag(band_b, ebv)


def star_applied_ebv(star: dict[str, Any]) -> tuple[float | None, str]:
    """Reddening applied to the star and its source label.

    Records without the explicit `ebv_star` key (built before this key existed) fall back to `ebv_sfd`,
    labelled `legacy_ebv_sfd` so the substitution is visible downstream.
    """
    if "ebv_star" in star:
        val = star.get("ebv_star")
        src = str(star.get("ebv_star_source") or ("missing" if val is None else "unspecified"))
        return (None if val is None else float(val)), src
    if star.get("ebv_sfd") is not None:
        return float(star["ebv_sfd"]), "legacy_ebv_sfd"
    return None, "missing"


def line_of_sight_ebv(star: dict[str, Any], *, sfd_to_sf11: float = 0.86) -> dict[str, Any] | None:
    """Total line-of-sight E(B-V) through the Galaxy, on the Schlafly & Finkbeiner (2011) scale.

    Priority: vetting enrichment `vetting_archive.dust` (ebv_sf11_los, else ebv_sfd_los x sfd_to_sf11), then the
    top-level `ebv_sfd` map value (x sfd_to_sf11). Returns None when no line-of-sight value exists.
    """
    dust = (star.get("vetting_archive") or {}).get("dust") or {}
    sf11 = dust.get("ebv_sf11_los")
    sfd = dust.get("ebv_sfd_los")
    if sf11 is not None:
        return {"ebv_los": float(sf11), "ebv_sfd_los": None if sfd is None else float(sfd), "source": "vetting_archive.dust.ebv_sf11_los"}
    if sfd is not None:
        return {"ebv_los": float(sfd) * sfd_to_sf11, "ebv_sfd_los": float(sfd), "source": f"vetting_archive.dust.ebv_sfd_los x {sfd_to_sf11}"}
    if star.get("ebv_sfd") is not None:
        return {"ebv_los": float(star["ebv_sfd"]) * sfd_to_sf11, "ebv_sfd_los": float(star["ebv_sfd"]), "source": f"ebv_sfd x {sfd_to_sf11}"}
    return None


def map3d_ebv(star: dict[str, Any]) -> dict[str, Any] | None:
    """3D-map E(B-V) at the star's distance from `vetting_archive.dust3d` when the query succeeded."""
    d3 = (star.get("vetting_archive") or {}).get("dust3d") or {}
    if d3.get("status") != "ok" or d3.get("ebv") is None:
        return None
    return {
        "ebv": float(d3["ebv"]),
        "ebv_err": float(d3.get("ebv_err") or 0.0),
        "map": str(d3.get("map", "map3d")),
    }


def dwarf_sequence_at_bp_rp(bp_rp0: float) -> dict[str, Any]:
    """Nearest dwarf spectral type, interpolated Teff and M_G at an intrinsic BP-RP (packaged dwarf table)."""
    rows = [r for r in load_dwarf_colour_table()["rows"] if r.get("bp_rp") is not None and r.get("m_g") is not None]
    rows.sort(key=lambda r: float(r["bp_rp"]))
    x = [float(r["bp_rp"]) for r in rows]
    nearest = min(rows, key=lambda r: abs(float(r["bp_rp"]) - bp_rp0))
    in_range = x[0] <= bp_rp0 <= x[-1]
    return {
        "spectral_type": str(nearest["spt"]),
        "teff_k": float(np.interp(bp_rp0, x, [float(r["teff_k"]) for r in rows])),
        "m_g_sequence": float(np.interp(bp_rp0, x, [float(r["m_g"]) for r in rows])),
        "in_table_range": bool(in_range),
    }


def _extinction_case(star: dict[str, Any], label: str, ebv: float, ebv_source: str) -> dict[str, Any]:
    m = star["mags"]
    dist = max(1.0, float(star["distance_pc"]))
    e_bprp = colour_excess("BP", "RP", ebv)
    bp_rp0 = float(m["BP"]) - float(m["RP"]) - e_bprp
    a_g = extinction_mag("G", ebv)
    m_g = float(m["G"]) - 5.0 * math.log10(dist / 10.0) - a_g
    seq = dwarf_sequence_at_bp_rp(bp_rp0)
    d_mg = m_g - seq["m_g_sequence"]
    return {
        "case": label,
        "ebv": round(float(ebv), 4),
        "ebv_source": ebv_source,
        "e_bp_rp": round(e_bprp, 4),
        "a_g": round(a_g, 4),
        "bp_rp_0": round(bp_rp0, 4),
        "m_g": round(m_g, 4),
        "spectral_type": seq["spectral_type"],
        "teff_k_from_colour": round(seq["teff_k"], 0),
        "m_g_dwarf_sequence": round(seq["m_g_sequence"], 4),
        "delta_m_g_vs_sequence": round(d_mg, 4),
        "luminosity_vs_sequence": "overluminous" if d_mg < 0.0 else "underluminous",
    }


def extinction_two_case_summary(
    star: dict[str, Any],
    *,
    min_los_excess_ebv: float,
    sfd_to_sf11: float = 0.86,
) -> dict[str, Any] | None:
    """Intrinsic colour, M_G and implied dwarf class under (i) the applied reddening and (ii) the full column.

    Produced only when the line-of-sight E(B-V) exceeds the applied E(B-V) by more than `min_los_excess_ebv`.
    The dwarf sequence used for the class and the reference M_G is the packaged Pecaut & Mamajek table.
    """
    if "mags" not in star or star.get("distance_pc") is None:
        return None
    applied, applied_src = star_applied_ebv(star)
    los = line_of_sight_ebv(star, sfd_to_sf11=sfd_to_sf11)
    if los is None:
        return None
    a_val = 0.0 if applied is None else float(applied)
    excess = float(los["ebv_los"]) - a_val
    if excess <= float(min_los_excess_ebv):
        return None
    cases = [
        _extinction_case(star, "applied", a_val, applied_src if applied is not None else "missing (treated as 0)"),
        _extinction_case(star, "full_line_of_sight", float(los["ebv_los"]), str(los["source"])),
    ]
    d3 = map3d_ebv(star)
    if d3 is not None:
        cases.append(_extinction_case(star, "map3d_at_distance", d3["ebv"], d3["map"]))
    c1, c2 = cases[0], cases[1]

    def _desc(c: dict[str, Any]) -> str:
        return (
            f"E(B-V)={c['ebv']:.3f}: (BP-RP)_0={c['bp_rp_0']:.3f} ({c['spectral_type']}, ~{c['teff_k_from_colour']:.0f} K), "
            f"M_G={c['m_g']:.2f}, {abs(c['delta_m_g_vs_sequence']):.2f} mag {'brighter' if c['delta_m_g_vs_sequence'] < 0 else 'fainter'} "
            f"than the dwarf sequence at that colour"
        )

    where = (
        f"a 3D map gives E(B-V)={d3['ebv']:.3f}+/-{d3['ebv_err']:.3f} at the star's distance"
        if d3 is not None
        else "no 3D dust map constrains how much of the column lies in front of the star, so its intrinsic colour and luminosity are uncertain between these cases"
    )
    return {
        "triggered_by": f"line-of-sight E(B-V) {los['ebv_los']:.3f} exceeds applied {a_val:.3f} by {excess:.3f} (> {min_los_excess_ebv:.3f})",
        "ebv_applied": round(a_val, 4),
        "ebv_applied_source": applied_src,
        "ebv_line_of_sight": round(float(los["ebv_los"]), 4),
        "ebv_line_of_sight_source": los["source"],
        "ebv_sfd_line_of_sight": los.get("ebv_sfd_los"),
        "cases": cases,
        "summary": f"(i) applied {_desc(c1)}; (ii) full line-of-sight column {_desc(c2)}; {where}.",
    }
