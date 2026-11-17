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
GAIA_TAP_URLS = [
    "https://gea.esac.esa.int/tap-server/tap/sync",
    "https://gaia.ari.uni-heidelberg.de/tap/sync",
    "https://gaia.aip.de/tap/sync",
]
IRSA_GATOR_URL = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
IRSA_TAP_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
IRSA_ZTF_URL = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
NOIRLAB_TAP_URL = "https://datalab.noirlab.edu/tap/sync"
MAST_PS1_DR2_URL = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.json"
VIZIER_TAP_URL = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
VIZIER_TAP_URLS = [
    "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
]
VIZIER_TSV_URLS = [
    "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",
    "https://vizier.cds.unistra.fr/viz-bin/asu-tsv",
]
SIMBAD_TAP_URL = "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync"
SIMBAD_TAP_URLS = [
    "https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
    "https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
]


def _parse_ipac_table(text: str) -> list[dict[str, Any]]:
    """Parse an IRSA Gator IPAC ASCII table (outfmt=1) into a list of row dicts."""
    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("\\")]
    if len(lines) < 4 or not lines[0].lstrip().startswith("|"):
        return []
    header = lines[0]
    spans: list[tuple[str, int, int]] = []
    pos = header.find("|")
    while pos != -1:
        nxt = header.find("|", pos + 1)
        if nxt == -1:
            break
        col = header[pos + 1 : nxt].strip()
        if col:
            spans.append((col, pos + 1, nxt))
        pos = nxt
    if not spans:
        return []

    rows: list[dict[str, Any]] = []
    for ln in lines:
        if ln.lstrip().startswith("|"):
            continue
        row: dict[str, Any] = {}
        for idx, (col, start, end) in enumerate(spans):
            raw = ln[start:] if idx == len(spans) - 1 else ln[start:end]
            val = raw.strip()
            row[col] = None if val in {"", "null", "NULL", "NaN", "nan"} else val
        rows.append(row)
    return rows


def _parse_vizier_tsv(text: str) -> list[dict[str, Any]]:
    """Parse a VizieR asu-tsv response into a list of dicts."""
    lines = [
        ln
        for ln in text.splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    data_lines = [ln for ln in lines if not set(ln.replace("\t", "").strip()).issubset({"-", " "})]
    if len(data_lines) < 2:
        return []
    cols = [c.strip() for c in data_lines[0].split("\t")]
    # Skip unit row immediately under header when it contains standard unit strings or blanks
    start_idx = 2 if len(data_lines) >= 3 else 1
    rows: list[dict[str, Any]] = []
    for ln in data_lines[start_idx:]:
        vals = [v.strip() for v in ln.split("\t")]
        if len(vals) < len(cols):
            vals.extend([""] * (len(cols) - len(vals)))
        row = {c: (None if v == "" else v) for c, v in zip(cols, vals)}
        rows.append(row)
    return rows


def _tap_query_cached(
    *,
    url: str | None = None,
    urls: list[str] | None = None,
    query: str,
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: int = 45,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Execute an ADQL sync TAP query with mirror fallback and content-addressed RawCache provenance."""
    candidates = urls or ([url] if url else [])
    if not candidates:
        raise ArchiveAccessError(f"No TAP URLs configured for {archive}", retryable=False)

    sig = hashlib.sha256(f"{'|'.join(candidates)}|{query}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    errors: list[str] = []
    for endpoint in candidates:
        try:
            resp = requests.post(
                endpoint,
                data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "json", "QUERY": query},
                timeout=timeout_s,
            )
            if resp.status_code != 200:
                errors.append(f"{endpoint}: HTTP {resp.status_code}")
                continue
            data = resp.json()
            cols = [m["name"] for m in data.get("metadata", [])]
            raw_rows = data.get("data", [])
            rows = [dict(zip(cols, r)) for r in raw_rows]
            # Some GAVO DaCHS mirrors echo header row when FORMAT=json isn't honored
            if rows and cols and str(rows[0].get(cols[0])) == cols[0]:
                errors.append(f"{endpoint}: non-JSON TAP echo")
                continue
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
        except Exception as exc:
            errors.append(f"{endpoint}: {type(exc).__name__}: {exc}")

    raise ArchiveAccessError(
        f"TAP query to {archive} failed across {len(candidates)} endpoint(s): {'; '.join(errors)}",
        retryable=True,
    )


def _gator_cone_cached(
    *,
    catalog: str,
    ra: float,
    dec: float,
    radius_arcsec: float,
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: int = 20,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Execute an IRSA Gator cone query with RawCache provenance."""
    query_desc = f"GATOR {catalog} cone({ra:.6f}, {dec:+.6f}, {radius_arcsec:.2f} arcsec)"
    sig = hashlib.sha256(f"{IRSA_GATOR_URL}|{query_desc}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    resp = requests.get(
        IRSA_GATOR_URL,
        params={
            "catalog": catalog,
            "spatial": "cone",
            "radius": f"{radius_arcsec:.2f}",
            "radunits": "arcsec",
            "objstr": f"{ra:.6f} {dec:+.6f}",
            "outfmt": "1",
        },
        timeout=timeout_s,
    )
    if resp.status_code != 200 or '[struct stat="ERROR"' in resp.text[:400]:
        raise ArchiveAccessError(
            f"IRSA Gator query for {catalog} failed (HTTP {resp.status_code}): {resp.text[:250]}",
            retryable=True,
        )
    rows = _parse_ipac_table(resp.text)
    payload = json.dumps(rows, sort_keys=True).encode("utf-8")
    manifest = cache.put(
        archive=archive,
        release=release,
        request_signature=sig,
        query_text=query_desc,
        payload=payload,
        row_count=len(rows),
        offline_fixture=False,
    )
    return rows, manifest, False


def _mast_ps1_cone_cached(
    *,
    ra: float,
    dec: float,
    radius_deg: float,
    cache: RawCache,
    timeout_s: int = 20,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Query MAST Pan-STARRS DR2 mean object catalog REST JSON API with caching."""
    query_desc = f"MAST_PS1_DR2 cone({ra:.6f}, {dec:+.6f}, {radius_deg:.6f} deg)"
    sig = hashlib.sha256(f"{MAST_PS1_DR2_URL}|{query_desc}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup("deep_optical", "ps1_dr2", sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    resp = requests.get(
        MAST_PS1_DR2_URL,
        params={"ra": f"{ra:.6f}", "dec": f"{dec:.6f}", "radius": f"{radius_deg:.6f}"},
        timeout=timeout_s,
    )
    if resp.status_code != 200:
        raise ArchiveAccessError(
            f"MAST PS1 DR2 returned HTTP {resp.status_code}: {resp.text[:250]}",
            retryable=True,
        )
    js = resp.json()
    cols = [c["name"] for c in js.get("info", [])]
    rows = [dict(zip(cols, raw)) for raw in js.get("data", [])]
    payload = json.dumps(rows, sort_keys=True).encode("utf-8")
    manifest = cache.put(
        archive="deep_optical",
        release="ps1_dr2",
        request_signature=sig,
        query_text=query_desc,
        payload=payload,
        row_count=len(rows),
        offline_fixture=False,
    )
    return rows, manifest, False


def _vizier_tsv_cone_cached(
    *,
    source: str,
    ra: float,
    dec: float,
    radius_arcsec: float,
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: int = 15,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Query VizieR asu-tsv across primary and CfA mirrors with RawCache provenance."""
    query_desc = f"VIZIER_TSV {source} cone({ra:.6f}, {dec:+.6f}, {radius_arcsec:.2f} arcsec)"
    sig = hashlib.sha256(f"{'|'.join(VIZIER_TSV_URLS)}|{query_desc}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    errors: list[str] = []
    for base in VIZIER_TSV_URLS:
        try:
            resp = requests.get(
                base,
                params={
                    "-source": source,
                    "-c": f"{ra:.6f} {dec:+.6f}",
                    "-c.r": f"{radius_arcsec:.2f}",
                    "-c.u": "arcsec",
                    "-out.max": "50",
                },
                timeout=timeout_s,
            )
            if resp.status_code != 200 or not resp.text.strip().startswith("#"):
                errors.append(f"{base}: HTTP {resp.status_code}")
                continue
            rows = _parse_vizier_tsv(resp.text)
            payload = json.dumps(rows, sort_keys=True).encode("utf-8")
            manifest = cache.put(
                archive=archive,
                release=release,
                request_signature=sig,
                query_text=query_desc,
                payload=payload,
                row_count=len(rows),
                offline_fixture=False,
            )
            return rows, manifest, False
        except Exception as exc:
            errors.append(f"{base}: {type(exc).__name__}: {exc}")

    raise ArchiveAccessError(
        f"VizieR TSV query for {source} failed across mirrors: {'; '.join(errors)}",
        retryable=True,
    )


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
    region: str = "halo",
    after_random_index: int | None = None,
) -> tuple[list[dict[str, Any]], list[ManifestRecord]]:
    """Fetch a real sample of old Galactic Halo, Galactic Bulge, or local high-PM stars from Gaia DR3."""
    n_lim = max(10, min(2000, int(limit)))
    reg = region.lower().strip()
    extra_cursor = f" AND g.random_index > {int(after_random_index)}" if after_random_index is not None else ""

    if reg == "halo":
        # True all-sky Galactic Halo: |b| > 20 deg, metal-poor [Fe/H] < -0.8, high v_tan > 150 km/s, randomized over whole sky
        where_region = f"""
          g.parallax_over_error > 10.0
          AND g.parallax BETWEEN 0.5 AND 20.0
          AND g.ruwe < 1.4
          AND ABS(g.b) > 20.0
          AND g.phot_g_mean_mag BETWEEN 7.5 AND 14.5
          AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
          AND g.mh_gspphot < -0.8
          AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 150.0
          {extra_cursor}
        """
        order_clause = "ORDER BY g.random_index ASC"
    elif reg == "bulge":
        # Galactic Bulge window around l ~ 0 deg, |b| 3..10 deg, bulge distances (~2.2-12.5 kpc)
        where_region = f"""
          (g.l <= 12.0 OR g.l >= 348.0)
          AND ABS(g.b) BETWEEN 3.0 AND 10.0
          AND g.parallax BETWEEN 0.08 AND 0.45
          AND g.parallax_over_error > 5.0
          AND g.ruwe < 1.4
          AND g.phot_g_mean_mag BETWEEN 8.0 AND 14.5
          AND g.teff_gspphot BETWEEN 3500.0 AND 7500.0
          {extra_cursor}
        """
        order_clause = "ORDER BY g.random_index ASC"
    elif reg == "local":
        # Legacy local high-velocity cut (v_tan > 65 km/s, |b| > 12 deg, 40-500 pc)
        where_region = """
          g.parallax_over_error > 15.0
          AND g.parallax BETWEEN 2.0 AND 25.0
          AND g.ruwe < 1.4
          AND ABS(g.b) > 12.0
          AND g.phot_g_mean_mag BETWEEN 7.5 AND 13.5
          AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
          AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 65.0
        """
        order_clause = "ORDER BY g.source_id ASC"
    else:
        raise ArchiveAccessError(f"Unsupported live region '{region}'; expected 'halo', 'bulge', or 'local'", retryable=False)

    query = f"""
    SELECT TOP {n_lim}
      g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.parallax_over_error,
      g.pmra, g.pmdec, g.ruwe, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
      g.teff_gspphot, g.logg_gspphot, g.mh_gspphot, g.ebpminrp_gspphot, g.random_index,
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
    WHERE {where_region}
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
    {order_clause}
    """
    rows, man, _ = _tap_query_cached(
        urls=GAIA_TAP_URLS,
        query=" ".join(query.split()),
        archive="gaia",
        release=f"{cfg.gaia.release}_{reg}",
        cache=cache,
        timeout_s=90,
    )
    stars = [_row_to_star_dict(r, gaia_release=cfg.gaia.release) for r in rows]
    return stars, [man]


def fetch_single_star_live(
    source_id: str,
    cfg: PipelineConfig,
    cache: RawCache,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Fetch a single real star from ESA Gaia + 2MASS + AllWISE (with mirror fallback) and enrich for vetting."""
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
        urls=GAIA_TAP_URLS,
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


def _safe_float(val: Any) -> float | None:
    if val is None:
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    if math.isnan(f) or f <= -900.0:
        return None
    return f


def enrich_star_live_for_vetting(
    star: dict[str, Any],
    cfg: PipelineConfig,
    cache: RawCache,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Enrich one real star using primary + backup mirrors across Gaia, SIMBAD, IRSA Gator/TAP, MAST PS1, VizieR, and Data Lab."""
    if star.get("live_enriched"):
        return star, []

    ra = float(star["ra_deg"])
    dec = float(star["dec_deg"])
    sid = str(star["source_id"])
    manifests: list[ManifestRecord] = []
    out = json.loads(json.dumps(star))

    # 1. Gaia DR3 neighbours inside 12" W4 beam + Gaia RV/NSS/variability metadata
    try:
        q_nbr = (
            f"SELECT source_id, ra, dec, phot_g_mean_mag, phot_rp_mean_mag "
            f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.003333))"
        )
        nbr_rows, m_nbr, _ = _tap_query_cached(urls=GAIA_TAP_URLS, query=q_nbr, archive="gaia", release=cfg.gaia.release, cache=cache, timeout_s=25)
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

    try:
        q_meta = (
            f"SELECT non_single_star, radial_velocity_error, phot_variable_flag "
            f"FROM gaiadr3.gaia_source WHERE source_id = {sid}"
        )
        m_rows, m_meta, _ = _tap_query_cached(urls=GAIA_TAP_URLS, query=q_meta, archive="gaia", release=cfg.gaia.release, cache=cache, timeout_s=20)
        manifests.append(m_meta)
        if m_rows:
            mr = m_rows[0]
            if mr.get("non_single_star") is not None:
                out["vetting_archive"]["nss_flag"] = int(mr["non_single_star"])
            rv_err = _safe_float(mr.get("radial_velocity_error"))
            if rv_err is not None:
                out["vetting_archive"]["rv_scatter_kms"] = rv_err
            if mr.get("phot_variable_flag"):
                out["vetting_archive"]["gaia_phot_variable_flag"] = str(mr["phot_variable_flag"])
    except Exception:
        pass

    # 2. SIMBAD classification (CfA Harvard primary, Strasbourg CDS backup)
    try:
        q_sim = (
            f"SELECT b.main_id, b.otype_txt, DISTANCE(POINT('ICRS', b.ra, b.dec), POINT('ICRS', {ra:.6f}, {dec:.6f}))*3600.0 AS sep_arcsec "
            f"FROM basic AS b WHERE 1=CONTAINS(POINT('ICRS', b.ra, b.dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0014)) "
            f"ORDER BY sep_arcsec ASC"
        )
        s_rows, m_sim, _ = _tap_query_cached(urls=SIMBAD_TAP_URLS, query=q_sim, archive="simbad", release="live", cache=cache, timeout_s=15)
        manifests.append(m_sim)
        if not s_rows:
            q_sim_id = (
                f"SELECT b.main_id, b.otype_txt FROM basic AS b "
                f"JOIN ident AS i ON b.oid = i.oidref WHERE i.id = 'Gaia DR3 {sid}'"
            )
            s_rows, m_sim_id, _ = _tap_query_cached(urls=SIMBAD_TAP_URLS, query=q_sim_id, archive="simbad", release="live", cache=cache, timeout_s=15)
            manifests.append(m_sim_id)
        if s_rows:
            out["simbad_type"] = str(s_rows[0].get("otype_txt") or "Star")
            out["name"] = str(s_rows[0].get("main_id") or out["name"])
            out["known_in_literature"] = True
    except Exception:
        pass

    # 3a. AllWISE coadd quality (IRSA Gator fast path -> IRSA TAP -> NOIRLab Data Lab backup)
    aw_rows: list[dict[str, Any]] = []
    try:
        aw_rows, m_aw, _ = _gator_cone_cached(
            catalog="allwise_p3as_psd",
            ra=ra,
            dec=dec,
            radius_arcsec=5.0,
            archive="allwise",
            release="allwise",
            cache=cache,
            timeout_s=20,
        )
        manifests.append(m_aw)
    except Exception:
        try:
            q_aw = (
                f"SELECT ra, dec, w4rchi2, w4m, w3rchi2, cc_flags, ext_flg, ph_qual "
                f"FROM allwise_p3as_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.001))"
            )
            aw_rows, m_aw, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_aw, archive="allwise", release="allwise", cache=cache, timeout_s=25)
            manifests.append(m_aw)
        except Exception:
            try:
                q_aw_dl = (
                    f"SELECT ra, dec, w4rchi2, w4m, cc_flags, ext_flg, ph_qual "
                    f"FROM allwise.source WHERE q3c_radial_query(ra, dec, {ra:.6f}, {dec:.6f}, 0.0015)"
                )
                aw_rows, m_aw, _ = _tap_query_cached(url=NOIRLAB_TAP_URL, query=q_aw_dl, archive="allwise", release="datalab", cache=cache, timeout_s=20)
                manifests.append(m_aw)
            except Exception:
                aw_rows = []

    if aw_rows:
        best = min(aw_rows, key=lambda r: angular_separation_arcsec(ra, dec, float(r["ra"]), float(r["dec"])))
        w4r = _safe_float(best.get("w4rchi2"))
        w4m = _safe_float(best.get("w4m"))
        if w4r is not None:
            out["allwise"]["w4rchi2"] = w4r
        if w4m is not None:
            out["allwise"]["w4m"] = int(w4m)
        if best.get("cc_flags"):
            out["allwise"]["cc_flags"] = str(best["cc_flags"])
        if best.get("ph_qual"):
            out["allwise"]["ph_qual"] = str(best["ph_qual"])
        if best.get("ext_flg") is not None:
            out["allwise"]["ext_flg"] = int(float(best["ext_flg"]))

    # 3b. WISE 2010 4-band Single-Exposure (L1b) W4 frames (IRSA Gator fast path -> IRSA TAP)
    l1b_rows: list[dict[str, Any]] = []
    try:
        l1b_rows, m_l1b, _ = _gator_cone_cached(
            catalog="allsky_4band_p1bs_psd",
            ra=ra,
            dec=dec,
            radius_arcsec=5.0,
            archive="allwise",
            release="allsky_4band_l1b",
            cache=cache,
            timeout_s=25,
        )
        manifests.append(m_l1b)
    except Exception:
        try:
            q_l1b = (
                f"SELECT mjd, w4mpro, w4sigmpro, w4rchi2, cc_flags, qual_frame "
                f"FROM allsky_4band_p1bs_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0014))"
            )
            l1b_rows, m_l1b, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_l1b, archive="allwise", release="allsky_4band_l1b", cache=cache, timeout_s=25)
            manifests.append(m_l1b)
        except Exception:
            l1b_rows = []

    if l1b_rows:
        w4_frames_mjy: list[float] = []
        for lr in l1b_rows:
            w4_mag = _safe_float(lr.get("w4mpro"))
            if w4_mag is not None:
                w4_frames_mjy.append(float(mag_to_flux_mjy("W4", w4_mag)))
        if w4_frames_mjy:
            out["wise_single_frames_w4_mjy"] = w4_frames_mjy

    # 3c. NEOWISE-R multi-epoch W1/W2 light curve (IRSA Gator fast path -> IRSA TAP)
    neo_rows: list[dict[str, Any]] = []
    try:
        neo_rows, m_neo, _ = _gator_cone_cached(
            catalog="neowiser_p1bs_psd",
            ra=ra,
            dec=dec,
            radius_arcsec=3.0,
            archive="neowise",
            release="neowise_r",
            cache=cache,
            timeout_s=25,
        )
        manifests.append(m_neo)
    except Exception:
        try:
            q_neo = (
                f"SELECT mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, qual_frame, qi_fact "
                f"FROM neowiser_p1bs_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.00083)) "
                f"AND qual_frame > 0"
            )
            neo_rows, m_neo, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_neo, archive="neowise", release="neowise_r", cache=cache, timeout_s=30)
            manifests.append(m_neo)
        except Exception:
            neo_rows = []

    if neo_rows:
        pts: list[tuple[float, float, float]] = []
        for row in neo_rows:
            qf = _safe_float(row.get("qual_frame"))
            mjd_val = _safe_float(row.get("mjd"))
            w1_val = _safe_float(row.get("w1mpro"))
            w2_val = _safe_float(row.get("w2mpro"))
            if (qf is None or qf > 0) and mjd_val is not None and w1_val is not None and w2_val is not None:
                pts.append((mjd_val, w1_val, w2_val))
        pts.sort(key=lambda x: x[0])
        if pts:
            # Bin into ~180-day semi-annual visits
            bins: list[tuple[float, list[tuple[float, float, float]]]] = []
            for mjd_val, w1_val, w2_val in pts:
                if not bins or mjd_val - bins[-1][0] > 120.0:
                    bins.append((mjd_val, [(mjd_val, w1_val, w2_val)]))
                else:
                    bins[-1][1].append((mjd_val, w1_val, w2_val))
            v_mjd, v_w1, v_w2, v_e1, v_e2 = [], [], [], [], []
            for _, grp in bins:
                v_mjd.append(float(np.median([x[0] for x in grp])))
                f1s = [mag_to_flux_mjy("W1", x[1]) for x in grp]
                f2s = [mag_to_flux_mjy("W2", x[2]) for x in grp]
                m1 = float(np.median(f1s))
                m2 = float(np.median(f2s))
                v_w1.append(m1)
                v_w2.append(m2)
                v_e1.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f1s) - m1))) / math.sqrt(len(grp)), m1 * 0.02))
                v_e2.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f2s) - m2))) / math.sqrt(len(grp)), m2 * 0.02))
            out["neowise"] = {"mjd": v_mjd, "w1_mjy": v_w1, "w2_mjy": v_w2, "w1_err_mjy": v_e1, "w2_err_mjy": v_e2}

    # 4. Second-instrument 18-25 um archival check (AKARI IRC/FIS, Spitzer SEIP)
    for cat, instr_name in (("slphotdr4", "Spitzer_MIPS24"), ("akari_irc", "AKARI_IRC_18um"), ("akari_fis", "AKARI_FIS")):
        try:
            sec_rows, m_sec, _ = _gator_cone_cached(
                catalog=cat,
                ra=ra,
                dec=dec,
                radius_arcsec=10.0,
                archive=cat,
                release="live",
                cache=cache,
                timeout_s=15,
            )
            manifests.append(m_sec)
            if sec_rows and out["vetting_archive"].get("second_instrument_ir") is None:
                out["vetting_archive"]["second_instrument_ir"] = {
                    "instrument": instr_name,
                    "n_matches": len(sec_rows),
                    "interpretation": "confirmed",
                }
                break
        except Exception:
            pass

    # 5. Radio continuum (NVSS via VizieR CfA/CDS TSV fast path -> VizieR TAP backup)
    r_rows: list[dict[str, Any]] = []
    try:
        r_rows, m_rad, _ = _vizier_tsv_cone_cached(
            source="VIII/65/nvss",
            ra=ra,
            dec=dec,
            radius_arcsec=10.0,
            archive="radio_continuum",
            release="nvss",
            cache=cache,
            timeout_s=15,
        )
        manifests.append(m_rad)
    except Exception:
        try:
            q_rad = (
                f"SELECT RAJ2000, DEJ2000, S1_4 AS flux_mjy FROM \"VIII/65/nvss\" "
                f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0028))"
            )
            r_rows, m_rad, _ = _tap_query_cached(urls=VIZIER_TAP_URLS, query=q_rad, archive="radio_continuum", release="nvss", cache=cache, timeout_s=15)
            manifests.append(m_rad)
        except Exception:
            r_rows = []

    if r_rows is not None:
        hits = []
        for rr in r_rows:
            r_ra = _safe_float(rr.get("RAJ2000"))
            r_dec = _safe_float(rr.get("DEJ2000"))
            if r_ra is None or r_dec is None:
                continue
            sep = angular_separation_arcsec(ra, dec, r_ra, r_dec)
            hits.append({"survey": "NVSS", "sep_arcsec": round(sep, 2), "flux_mjy": rr.get("S1_4") or rr.get("flux_mjy")})
        out["vetting_archive"]["radio_sources_within_10arcsec"] = hits

    # 6. Deep sub-arcsec optical imaging (MAST PS1 DR2 primary -> VizieR CfA TSV -> VizieR TAP -> NOIRLab LS DR10)
    if dec > -30.0:
        deep_evaluated = False
        offset_galaxy: dict[str, Any] | None = None
        try:
            p_rows, m_ps1, _ = _mast_ps1_cone_cached(ra=ra, dec=dec, radius_deg=0.0028, cache=cache, timeout_s=20)
            manifests.append(m_ps1)
            # Filter out single-frame diffraction/ghost artifacts (require nDetections >= 2)
            valid_rows = [
                pr
                for pr in p_rows
                if _safe_float(pr.get("raMean")) is not None
                and _safe_float(pr.get("decMean")) is not None
                and int(pr.get("nDetections") or 0) >= 2
            ]
            if valid_rows:
                deep_evaluated = True
                for pr in valid_rows:
                    sep = angular_separation_arcsec(ra, dec, float(pr["raMean"]), float(pr["decMean"]))
                    r_psf = _safe_float(pr.get("rMeanPSFMag"))
                    r_kron = _safe_float(pr.get("rMeanKronMag"))
                    if sep < 1.5 and r_psf is not None and r_kron is not None and (r_psf - r_kron) > 0.35:
                        offset_galaxy = {"sep_arcsec": round(sep, 2), "r_mag": r_kron, "extended_psf_kron_diff": round(r_psf - r_kron, 3)}
                        break
                    if 1.8 <= sep <= 8.0 and (r_psf is None or r_psf > 17.5):
                        offset_galaxy = {"sep_arcsec": round(sep, 2), "r_mag": r_psf}
                        break
        except Exception:
            pass

        if not deep_evaluated:
            try:
                p_rows, m_ps1, _ = _vizier_tsv_cone_cached(
                    source="II/349/ps1",
                    ra=ra,
                    dec=dec,
                    radius_arcsec=10.0,
                    archive="deep_optical",
                    release="ps1_vizier",
                    cache=cache,
                    timeout_s=15,
                )
                manifests.append(m_ps1)
                if p_rows:
                    deep_evaluated = True
                    for pr in p_rows:
                        p_ra = _safe_float(pr.get("RAJ2000"))
                        p_dec = _safe_float(pr.get("DEJ2000"))
                        if p_ra is None or p_dec is None:
                            continue
                        sep = angular_separation_arcsec(ra, dec, p_ra, p_dec)
                        r_mag = _safe_float(pr.get("rmag"))
                        if 1.8 <= sep <= 8.0 and (r_mag is None or r_mag > 17.5):
                            offset_galaxy = {"sep_arcsec": round(sep, 2), "r_mag": r_mag}
                            break
            except Exception:
                pass

        if not deep_evaluated:
            try:
                q_ls = (
                    f"SELECT ls_id, ra, dec, type, mag_r FROM ls_dr10.tractor "
                    f"WHERE q3c_radial_query(ra, dec, {ra:.6f}, {dec:.6f}, 0.0028)"
                )
                ls_rows, m_ls, _ = _tap_query_cached(url=NOIRLAB_TAP_URL, query=q_ls, archive="deep_optical", release="ls_dr10", cache=cache, timeout_s=20)
                manifests.append(m_ls)
                if ls_rows:
                    deep_evaluated = True
                    for lr in ls_rows:
                        sep = angular_separation_arcsec(ra, dec, float(lr["ra"]), float(lr["dec"]))
                        mtype = str(lr.get("type") or "")
                        rmag = _safe_float(lr.get("mag_r"))
                        if (sep < 1.5 and mtype not in {"PSF", "DUP", ""}) or (1.8 <= sep <= 8.0 and (rmag is None or rmag > 17.5)):
                            offset_galaxy = {"sep_arcsec": round(sep, 2), "r_mag": rmag, "ls_type": mtype}
                            break
            except Exception:
                pass

        out["vetting_archive"]["deep_imaging_available"] = deep_evaluated
        out["vetting_archive"]["deep_imaging_offset_galaxy"] = offset_galaxy
    else:
        out["vetting_archive"]["deep_imaging_available"] = False

    out["live_enriched"] = True
    return out, manifests
