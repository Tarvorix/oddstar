"""Live multi-archive TAP ingestion and vetting enrichment for real sky searches."""

from __future__ import annotations

import hashlib
import json
import math
from typing import Any

import numpy as np
import requests

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io.cache import RawCache
from oddstar.io.physics import (
    BAND_EXTINCTION_AV,
    angular_separation_arcsec,
    mag_err_to_flux_err_mjy,
    mag_to_flux_mjy,
)
from oddstar.schemas import ManifestRecord

GAIA_TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
IRSA_TAP_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
VIZIER_TAP_URL = "http://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
SIMBAD_TAP_URL = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"


def _tap_query_cached(
    *,
    url: str,
    query: str,
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: int = 45,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Execute an ADQL sync TAP query with content-addressed RawCache provenance."""
    sig = hashlib.sha256(f"{url}|{query}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    resp = requests.post(
        url,
        data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": query},
        timeout=timeout_s,
    )
    if resp.status_code != 200:
        raise ArchiveAccessError(
            f"TAP query to {archive} ({url}) returned HTTP {resp.status_code}: {resp.text[:300]}",
            retryable=True,
        )
    data = resp.json()
    cols = [m["name"] for m in data.get("metadata", [])]
    raw_rows = data.get("data", [])
    rows = [dict(zip(cols, r)) for r in raw_rows]
    payload = json.dumps(rows, sort_keys=True).encode("utf-8")
    manifest = cache.put(
        archive=archive,
        release=release,
        request_signature=sig,
        query_text=query,
        payload=payload,
        row_count=len(rows),
        offline_fixture=False,
    )
    return rows, manifest, False


def _row_to_star_dict(row: dict[str, Any], gaia_release: str = "dr3") -> dict[str, Any]:
    """Convert joined Gaia+AllWISE+2MASS TAP row into oddstar star dictionary."""
    sid = str(row["source_id"])
    ra = float(row["ra"])
    dec = float(row["dec"])
    gal_b = float(row.get("b", 30.0) if row.get("b") is not None else 30.0)
    plx = max(0.1, float(row["parallax"]))
    plx_err = max(1e-4, float(row.get("parallax_error") or (plx / max(1.0, float(row.get("parallax_over_error") or 20.0)))))
    plx_snr = float(row.get("parallax_over_error") or (plx / plx_err))
    dist_pc = 1000.0 / plx
    pmra = float(row.get("pmra") or 0.0)
    pmdec = float(row.get("pmdec") or 0.0)
    ruwe = float(row.get("ruwe") or 1.0)
    teff = float(row.get("teff_gspphot") or 5500.0)
    logg = float(row.get("logg_gspphot") or 4.4)
    feh = float(row.get("mh_gspphot") or -0.5)
    ebv = max(0.0, float(row.get("ebpminrp_gspphot") or 0.03) / 1.31)
    av = 3.1 * ebv

    mags = {
        "G": float(row["phot_g_mean_mag"]),
        "BP": float(row["phot_bp_mean_mag"]),
        "RP": float(row["phot_rp_mean_mag"]),
        "J": float(row["j_m"]),
        "H": float(row["h_m"]),
        "Ks": float(row["ks_m"]),
        "W1": float(row["w1mpro"]),
        "W2": float(row["w2mpro"]),
        "W3": float(row["w3mpro"]),
        "W4": float(row["w4mpro"]),
    }

    # Handle AllWISE upper limits (where *mpro_error is None -> assign large uncertainty so chi reflects upper limit)
    mag_errs = {
        "G": 0.005,
        "BP": 0.010,
        "RP": 0.010,
        "J": float(row.get("j_msigcom") or 0.025),
        "H": float(row.get("h_msigcom") or 0.025),
        "Ks": float(row.get("ks_msigcom") or 0.025),
        "W1": float(row.get("w1mpro_error") or 0.03),
        "W2": float(row.get("w2mpro_error") or 0.03),
        "W3": float(row.get("w3mpro_error") or 0.50),
        "W4": float(row.get("w4mpro_error") or 1.50),
    }

    flux_mjy: dict[str, float] = {}
    flux_err_mjy: dict[str, float] = {}
    for band, mag in mags.items():
        flx = mag_to_flux_mjy(band, mag)
        ferr = mag_err_to_flux_err_mjy(flx, mag_errs[band])
        flux_mjy[band] = flx
        flux_err_mjy[band] = max(0.02, ferr)

    abs_g = mags["G"] - 5.0 * math.log10(max(1.0, dist_pc) / 10.0) - av * BAND_EXTINCTION_AV["G"]
    lum_lsun = 10.0 ** (-0.4 * (abs_g - 4.68))

    wise_ra = float(row.get("wise_ra") or ra)
    wise_dec = float(row.get("wise_dec") or dec)
    ph_qual = str(row.get("ph_qual") or "AAAA")
    tmass_ph = str(row.get("tmass_ph_qual") or "AAA")

    return {
        "source_id": sid,
        "name": f"Gaia DR3 {sid}",
        "gaia_release": gaia_release,
        "ra_deg": ra,
        "dec_deg": dec,
        "gal_b_deg": gal_b,
        "distance_pc": dist_pc,
        "parallax_mas": plx,
        "parallax_err_mas": plx_err,
        "parallax_snr": plx_snr,
        "pmra_mas_yr": pmra,
        "pmdec_mas_yr": pmdec,
        "ruwe": ruwe,
        "teff_k": teff,
        "teff_spec_k": teff,
        "logg": logg,
        "feh": feh,
        "luminosity_lsun": lum_lsun,
        "abs_g_mag": abs_g,
        "ebv_sfd": ebv,
        "ebv_planck": ebv,
        "mags": mags,
        "mag_errs": mag_errs,
        "flux_mjy": flux_mjy,
        "flux_err_mjy": flux_err_mjy,
        "v_tan_kms": 4.74047 * math.hypot(pmra, pmdec) / plx,
        "allwise_designation": str(row.get("wise_designation") or ""),
        "allwise": {
            "ra_deg": wise_ra,
            "dec_deg": wise_dec,
            "cc_flags": str(row.get("cc_flags") or "0000"),
            "ext_flg": int(row.get("ext_flag") or 0),
            "ph_qual": ph_qual,
            "w4rchi2": float(row.get("w4rchi2") or 1.0),
            "w4m": int(row.get("w4m") or 12),
            "n_neighbours_gaia_xmatch": 1,
            "centroid_offset_arcsec": angular_separation_arcsec(ra, dec, wise_ra, wise_dec),
            "cirrus_bg_rms_mjy": 0.4,
        },
        "wise_single_frames_w4_mjy": [],
        "neowise": {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []},
        "lightcurves": {"tess_time_btjd": [], "tess_flux_qlp": [], "tess_flux_tesscut": [], "ztf_time_mjd": [], "ztf_flux": []},
        "vetting_archive": {
            "radio_sources_within_10arcsec": [],
            "deep_imaging_offset_galaxy": None,
            "deep_imaging_available": False,
            "neighbours_in_w4_beam": [],
            "local_red_source_density_per_arcsec2": 8e-5,
            "moving_objects_2010": [],
            "second_instrument_ir": None,
            "galex_uv_excess_mag": 0.0,
            "xray_luminosity_erg_s": 0.0,
            "young_association": None,
            "rv_scatter_kms": 0.3,
            "nss_flag": 0,
        },
        "spectra": {},
        "vanishing_data": {},
        "simbad_type": "Star",
        "known_in_literature": False,
        "2mass_clean": all(c in "ABCD" for c in tmass_ph[:3]),
        "live_enriched": False,
    }


def fetch_live_halo_bulge_sample(
    cfg: PipelineConfig,
    cache: RawCache,
    *,
    limit: int = 300,
) -> tuple[list[dict[str, Any]], list[ManifestRecord]]:
    """Fetch a real sample of old high-velocity Galactic Halo / Bulge-Halo stars from ESA Gaia DR3."""
    n_lim = max(10, min(2000, int(limit)))
    # Kinematic Halo / Thick-Disk / Bulge-halo cut: v_tan > 65 km/s, clean parallax & RUWE, |b| > 12 deg,
    # unique AllWISE + 2MASS matches, real W1-W3 detections (w3mpro_error NOT NULL)
    query = f"""
    SELECT TOP {n_lim}
      g.source_id, g.ra, g.dec, g.b, g.parallax, g.parallax_error, g.parallax_over_error,
      g.pmra, g.pmdec, g.ruwe, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
      g.teff_gspphot, g.logg_gspphot, g.mh_gspphot, g.ebpminrp_gspphot,
      aw.designation AS wise_designation,
      aw.w1mpro, aw.w1mpro_error, aw.w2mpro, aw.w2mpro_error,
      aw.w3mpro, aw.w3mpro_error, aw.w4mpro, aw.w4mpro_error,
      aw.cc_flags, aw.ext_flag, aw.ph_qual, aw.ra AS wise_ra, aw.dec AS wise_dec,
      tm.j_m, tm.j_msigcom, tm.h_m, tm.h_msigcom, tm.ks_m, tm.ks_msigcom, tm.ph_qual AS tmass_ph_qual
    FROM gaiadr3.gaia_source AS g
    JOIN gaiadr3.allwise_best_neighbour AS awn ON g.source_id = awn.source_id
    JOIN gaiadr1.allwise_original_valid AS aw ON awn.allwise_oid = aw.allwise_oid
    JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS tmn ON g.source_id = tmn.source_id
    JOIN gaiadr3.tmass_psc_xsc_join AS tmj ON tmn.clean_tmass_psc_xsc_oid = tmj.clean_tmass_psc_xsc_oid
    JOIN gaiadr1.tmass_original_valid AS tm ON tmj.original_psc_source_id = tm.designation
    WHERE g.parallax_over_error > 15.0
      AND g.parallax BETWEEN 2.0 AND 25.0
      AND g.ruwe < 1.4
      AND ABS(g.b) > 12.0
      AND g.phot_g_mean_mag BETWEEN 7.5 AND 13.5
      AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
      AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 65.0
      AND awn.number_of_neighbours = 1
      AND awn.number_of_mates = 0
      AND tmn.number_of_neighbours = 1
      AND tmn.number_of_mates = 0
      AND aw.w1mpro_error IS NOT NULL
      AND aw.w2mpro_error IS NOT NULL
      AND aw.w3mpro_error IS NOT NULL
      AND tm.j_msigcom IS NOT NULL
      AND tm.h_msigcom IS NOT NULL
      AND tm.ks_msigcom IS NOT NULL
    ORDER BY g.source_id ASC
    """
    rows, man, _ = _tap_query_cached(
        url=GAIA_TAP_URL,
        query=" ".join(query.split()),
        archive="gaia",
        release=cfg.gaia.release,
        cache=cache,
        timeout_s=60,
    )
    stars = [_row_to_star_dict(r, gaia_release=cfg.gaia.release) for r in rows]
    return stars, [man]


def fetch_single_star_live(
    source_id: str,
    cfg: PipelineConfig,
    cache: RawCache,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Fetch a single real star from ESA Gaia + 2MASS + AllWISE and enrich for vetting."""
    clean = source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
    query = f"""
    SELECT TOP 1
      g.source_id, g.ra, g.dec, g.b, g.parallax, g.parallax_error, g.parallax_over_error,
      g.pmra, g.pmdec, g.ruwe, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
      g.teff_gspphot, g.logg_gspphot, g.mh_gspphot, g.ebpminrp_gspphot,
      aw.designation AS wise_designation,
      aw.w1mpro, aw.w1mpro_error, aw.w2mpro, aw.w2mpro_error,
      aw.w3mpro, aw.w3mpro_error, aw.w4mpro, aw.w4mpro_error,
      aw.cc_flags, aw.ext_flag, aw.ph_qual, aw.ra AS wise_ra, aw.dec AS wise_dec,
      tm.j_m, tm.j_msigcom, tm.h_m, tm.h_msigcom, tm.ks_m, tm.ks_msigcom, tm.ph_qual AS tmass_ph_qual
    FROM gaiadr3.gaia_source AS g
    JOIN gaiadr3.allwise_best_neighbour AS awn ON g.source_id = awn.source_id
    JOIN gaiadr1.allwise_original_valid AS aw ON awn.allwise_oid = aw.allwise_oid
    JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS tmn ON g.source_id = tmn.source_id
    JOIN gaiadr3.tmass_psc_xsc_join AS tmj ON tmn.clean_tmass_psc_xsc_oid = tmj.clean_tmass_psc_xsc_oid
    JOIN gaiadr1.tmass_original_valid AS tm ON tmj.original_psc_source_id = tm.designation
    WHERE g.source_id = {clean}
    """
    rows, man, _ = _tap_query_cached(
        url=GAIA_TAP_URL,
        query=" ".join(query.split()),
        archive="gaia",
        release=cfg.gaia.release,
        cache=cache,
        timeout_s=45,
    )
    if not rows:
        raise ArchiveAccessError(f"Gaia DR3 {clean} not found with unique AllWISE+2MASS match", retryable=False)
    star = _row_to_star_dict(rows[0], gaia_release=cfg.gaia.release)
    enriched, extra_mans = enrich_star_live_for_vetting(star, cfg, cache)
    return enriched, [man, *extra_mans]


def enrich_star_live_for_vetting(
    star: dict[str, Any],
    cfg: PipelineConfig,
    cache: RawCache,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Enrich one real star with live SIMBAD, Gaia neighbours, IRSA AllWISE/NEOWISE, and VizieR vetting data."""
    if star.get("live_enriched"):
        return star, []

    ra = float(star["ra_deg"])
    dec = float(star["dec_deg"])
    sid = str(star["source_id"])
    manifests: list[ManifestRecord] = []
    out = json.loads(json.dumps(star))

    # 1. Gaia DR3 neighbours inside 12" W4 beam
    try:
        q_nbr = (
            f"SELECT source_id, ra, dec, phot_g_mean_mag, phot_rp_mean_mag "
            f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.003333))"
        )
        nbr_rows, m_nbr, _ = _tap_query_cached(url=GAIA_TAP_URL, query=q_nbr, archive="gaia", release=cfg.gaia.release, cache=cache, timeout_s=25)
        manifests.append(m_nbr)
        nbs = []
        for nr in nbr_rows:
            if str(nr["source_id"]) == sid:
                continue
            sep = angular_separation_arcsec(ra, dec, float(nr["ra"]), float(nr["dec"]))
            nbs.append({"id": f"Gaia DR3 {nr['source_id']}", "sep_arcsec": round(sep, 2), "g_mag": nr.get("phot_g_mean_mag"), "flux_frac": 0.10})
        out["vetting_archive"]["neighbours_in_w4_beam"] = nbs
        out["allwise"]["n_neighbours_gaia_xmatch"] = 1
    except Exception:
        pass

    # 2. SIMBAD classification via SIMBAD TAP
    try:
        q_sim = (
            f"SELECT b.main_id, b.otype_txt, DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', {ra:.6f}, {dec:.6f}))*3600.0 AS sep_arcsec "
            f"FROM basic AS b WHERE 1=CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0014)) "
            f"ORDER BY sep_arcsec ASC"
        )
        s_rows, m_sim, _ = _tap_query_cached(url=SIMBAD_TAP_URL, query=q_sim, archive="simbad", release="live", cache=cache, timeout_s=20)
        manifests.append(m_sim)
        if s_rows:
            out["simbad_type"] = str(s_rows[0].get("otype_txt") or "Star")
            out["name"] = str(s_rows[0].get("main_id") or out["name"])
            out["known_in_literature"] = True
    except Exception:
        pass

    # 3. IRSA AllWISE exact w4rchi2 / w4m & NEOWISE-R multi-epoch W1/W2 light curve
    try:
        q_aw = (
            f"SELECT ra, dec, w4rchi2, w4m, w3rchi2, cc_flags, ext_flg, ph_qual "
            f"FROM allwise_p3as_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.001))"
        )
        aw_rows, m_aw, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_aw, archive="allwise", release="allwise", cache=cache, timeout_s=30)
        manifests.append(m_aw)
        if aw_rows:
            best = min(aw_rows, key=lambda r: angular_separation_arcsec(ra, dec, float(r["ra"]), float(r["dec"])))
            out["allwise"]["w4rchi2"] = float(best.get("w4rchi2") or 1.0)
            out["allwise"]["w4m"] = int(best.get("w4m") or 12)
    except Exception:
        pass

    try:
        q_neo = (
            f"SELECT mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, qual_frame, qi_fact "
            f"FROM neowiser_p1bs_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.00083)) "
            f"AND qual_frame > 0"
        )
        neo_rows, m_neo, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_neo, archive="neowise", release="neowise_r", cache=cache, timeout_s=35)
        manifests.append(m_neo)
        if neo_rows:
            pts = sorted((float(r["mjd"]), r) for r in neo_rows if r.get("w1mpro") is not None and r.get("w2mpro") is not None)
            if pts:
                # Bin into ~180-day semi-annual visits
                bins: list[tuple[float, list[dict[str, Any]]]] = []
                for mjd_val, row in pts:
                    if not bins or mjd_val - bins[-1][0] > 120.0:
                        bins.append((mjd_val, [row]))
                    else:
                        bins[-1][1].append(row)
                v_mjd, v_w1, v_w2, v_e1, v_e2 = [], [], [], [], []
                for _, grp in bins:
                    v_mjd.append(float(np.median([float(x["mjd"]) for x in grp])))
                    f1s = [mag_to_flux_mjy("W1", float(x["w1mpro"])) for x in grp]
                    f2s = [mag_to_flux_mjy("W2", float(x["w2mpro"])) for x in grp]
                    m1 = float(np.median(f1s))
                    m2 = float(np.median(f2s))
                    v_w1.append(m1)
                    v_w2.append(m2)
                    v_e1.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f1s) - m1))) / math.sqrt(len(grp)), m1 * 0.02))
                    v_e2.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f2s) - m2))) / math.sqrt(len(grp)), m2 * 0.02))
                out["neowise"] = {"mjd": v_mjd, "w1_mjy": v_w1, "w2_mjy": v_w2, "w1_err_mjy": v_e1, "w2_err_mjy": v_e2}
    except Exception:
        pass

    # 4. VizieR radio continuum (NVSS & FIRST) + Pan-STARRS deep optical check
    try:
        q_rad = (
            f"SELECT RAJ2000, DEJ2000, S1_4 AS flux_mjy FROM \"VIII/65/nvss\" "
            f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0028))"
        )
        r_rows, m_rad, _ = _tap_query_cached(url=VIZIER_TAP_URL, query=q_rad, archive="radio_continuum", release="nvss", cache=cache, timeout_s=20)
        manifests.append(m_rad)
        hits = []
        for rr in r_rows:
            sep = angular_separation_arcsec(ra, dec, float(rr["RAJ2000"]), float(rr["DEJ2000"]))
            hits.append({"survey": "NVSS", "sep_arcsec": round(sep, 2), "flux_mjy": rr.get("flux_mjy")})
        out["vetting_archive"]["radio_sources_within_10arcsec"] = hits
    except Exception:
        pass

    try:
        if dec > -30.0:
            q_ps1 = (
                f"SELECT RAJ2000, DEJ2000, rmag FROM \"II/349/ps1\" "
                f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0028))"
            )
            p_rows, m_ps1, _ = _tap_query_cached(url=VIZIER_TAP_URL, query=q_ps1, archive="deep_optical", release="ps1", cache=cache, timeout_s=20)
            manifests.append(m_ps1)
            out["vetting_archive"]["deep_imaging_available"] = len(p_rows) > 0
            # Check for faint offset red background sources (> 1.8" and < 8" away, r > 17)
            off = [
                (angular_separation_arcsec(ra, dec, float(pr["RAJ2000"]), float(pr["DEJ2000"])), pr)
                for pr in p_rows
            ]
            faint_off = [(sep, pr) for sep, pr in off if 1.8 <= sep <= 8.0 and (pr.get("rmag") is None or float(pr["rmag"]) > 17.5)]
            if faint_off:
                sep0, pr0 = min(faint_off, key=lambda x: x[0])
                out["vetting_archive"]["deep_imaging_offset_galaxy"] = {"sep_arcsec": round(sep0, 2), "r_mag": pr0.get("rmag")}
        else:
            out["vetting_archive"]["deep_imaging_available"] = False
    except Exception:
        pass

    out["live_enriched"] = True
    return out, manifests
