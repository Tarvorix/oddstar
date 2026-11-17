"""Live multi-archive TAP ingestion and vetting enrichment for real sky searches."""

from __future__ import annotations

import hashlib
import json
import math
import time
from typing import Any

import numpy as np
import requests

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError
from oddstar.io.cache import RawCache
from oddstar.io.optical_timeseries import (
    fetch_asassn_lightcurve,
    fetch_tess_pipeline_pair,
    fetch_ztf_lightcurve,
)
from oddstar.io.physics import (
    BAND_EXTINCTION_AV,
    angular_separation_arcsec,
    mag_err_to_flux_err_mjy,
    mag_to_flux_mjy,
    propagate_position,
)
from oddstar.io.cutouts import fetch_allwise_cutouts, fetch_ps1_stack_cutout, measure_wise_cutouts
from oddstar.io.skybot import skybot_cone_search
from oddstar.io.storage import equatorial_to_galactic
from oddstar.io.vetting_queries import fetch_irsa_dust, fetch_subarcsec_imaging, target_excess_mjy
from oddstar.schemas import ManifestRecord

GAIA_TAP_URL = "https://gea.esac.esa.int/tap-server/tap/sync"
GAIA_TAP_URLS = [
    "https://gea.esac.esa.int/tap-server/tap/sync",
    "https://gaia.ari.uni-heidelberg.de/tap/sync",
    "https://gaia.aip.de/tap/sync",
]
IRSA_GATOR_URL = "https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query"
IRSA_TAP_URL = "https://irsa.ipac.caltech.edu/TAP/sync"
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


def _tap_csv_query_cached(
    *,
    url: str,
    query: str,
    archive: str,
    release: str,
    cache: RawCache,
    timeout_s: int = 25,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Execute a TAP sync query with FORMAT=csv (used for NOIRLab DataLab) and RawCache provenance."""
    import csv
    import io

    sig = hashlib.sha256(f"{url}|csv|{query}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    resp = requests.post(
        url,
        data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query},
        timeout=timeout_s,
    )
    if resp.status_code != 200 or resp.text.lstrip().startswith("<"):
        raise ArchiveAccessError(
            f"CSV TAP query to {archive} failed (HTTP {resp.status_code}): {resp.text[:250]}",
            retryable=True,
        )
    reader = csv.DictReader(io.StringIO(resp.text))
    rows = [
        {k: (None if v in {"", "null", "NULL", "NaN", "nan"} else v) for k, v in r.items()}
        for r in reader
    ]
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
    selcols: str | None = None,
) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
    """Execute an IRSA Gator cone query with RawCache provenance (optionally restricted to `selcols`)."""
    query_desc = f"GATOR {catalog} cone({ra:.6f}, {dec:+.6f}, {radius_arcsec:.2f} arcsec)"
    if selcols:
        query_desc += f" selcols={selcols}"
    sig = hashlib.sha256(f"{IRSA_GATOR_URL}|{query_desc}".encode("utf-8")).hexdigest()[:24]
    cached = cache.lookup(archive, release, sig)
    if cached is not None and not cached.offline_fixture:
        rows = json.loads(cache.read_bytes(cached).decode("utf-8"))
        return rows, cached, True

    gator_params = {
        "catalog": catalog,
        "spatial": "cone",
        "radius": f"{radius_arcsec:.2f}",
        "radunits": "arcsec",
        "objstr": f"{ra:.6f} {dec:+.6f}",
        "outfmt": "1",
    }
    if selcols:
        gator_params["selcols"] = selcols
    resp = requests.get(IRSA_GATOR_URL, params=gator_params, timeout=timeout_s)
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
    teff_source = "gspphot" if row.get("teff_gspphot") else "default_5500K"
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

    prop_ra, prop_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2010.2)
    ph_qual = str(row.get("ph_qual") or "AAAA")
    tmass_ph = str(row.get("tmass_ph_qual") or "AAA")

    # Only values present in the joined Gaia/AllWISE/2MASS row are recorded here. Anything that needs a further
    # archive query (W4 profile-fit chi2, frame counts, image centroids, background rms, neighbours, radio,
    # imaging, ephemerides, youth indicators, RV/NSS) is left ABSENT so vetting reports `not_evaluated` until
    # `enrich_star_live_for_vetting` records a successful query in `vetting_archive.query_status`.
    allwise_rec: dict[str, Any] = {
        "cc_flags": str(row.get("cc_flags") or "0000"),
        "ext_flg": int(row.get("ext_flag") or 0),
        "ph_qual": ph_qual,
    }
    if row.get("wise_ra") is not None and row.get("wise_dec") is not None:
        wise_ra = float(row["wise_ra"])
        wise_dec = float(row["wise_dec"])
        allwise_rec["ra_deg"] = wise_ra
        allwise_rec["dec_deg"] = wise_dec
        # AllWISE catalogue position (dominated by W1/W2), not a W3/W4 image centroid
        allwise_rec["centroid_offset_arcsec"] = angular_separation_arcsec(prop_ra, prop_dec, wise_ra, wise_dec)
        allwise_rec["centroid_source"] = "allwise_catalogue_position"
    if row.get("allwise_n_neighbours") is not None:
        allwise_rec["n_neighbours_gaia_xmatch"] = int(row["allwise_n_neighbours"])

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
        "teff_k_source": teff_source,
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
        "allwise": allwise_rec,
        "wise_single_frames_w4_mjy": [],
        "neowise": {"mjd": [], "w1_mjy": [], "w2_mjy": [], "w1_err_mjy": [], "w2_err_mjy": []},
        "lightcurves": {"tess_time_btjd": [], "tess_flux_qlp": [], "tess_flux_tesscut": [], "ztf_time_mjd": [], "ztf_flux": []},
        "vetting_archive": {
            "second_instrument_ir": None,
            "query_status": {},
        },
        "spectra": {},
        "vanishing_data": {},
        "simbad_type": "Star",
        "known_in_literature": False,
        "2mass_clean": all(c in "ABCD" for c in tmass_ph[:3]),
        "live_enriched": False,
    }


def _select_columns_sql() -> str:
    return """
      g.source_id, g.ra, g.dec, g.l, g.b, g.parallax, g.parallax_error, g.parallax_over_error,
      g.pmra, g.pmdec, g.ruwe, g.phot_g_mean_mag, g.phot_bp_mean_mag, g.phot_rp_mean_mag,
      g.teff_gspphot, g.logg_gspphot, g.mh_gspphot, g.ebpminrp_gspphot, g.random_index,
      aw.designation AS wise_designation,
      aw.w1mpro, aw.w1mpro_error, aw.w2mpro, aw.w2mpro_error,
      aw.w3mpro, aw.w3mpro_error, aw.w4mpro, aw.w4mpro_error,
      aw.cc_flags, aw.ext_flag, aw.ph_qual, aw.ra AS wise_ra, aw.dec AS wise_dec,
      awn.number_of_neighbours AS allwise_n_neighbours,
      tm.j_m, tm.j_msigcom, tm.h_m, tm.h_msigcom, tm.ks_m, tm.ks_msigcom, tm.ph_qual AS tmass_ph_qual
    """


def _join_and_phot_cuts_sql() -> str:
    return """
    FROM gaiadr3.gaia_source AS g
    JOIN gaiadr3.allwise_best_neighbour AS awn ON g.source_id = awn.source_id
    JOIN gaiadr1.allwise_original_valid AS aw ON awn.allwise_oid = aw.allwise_oid
    JOIN gaiadr3.tmass_psc_xsc_best_neighbour AS tmn ON g.source_id = tmn.source_id
    JOIN gaiadr3.tmass_psc_xsc_join AS tmj ON tmn.clean_tmass_psc_xsc_oid = tmj.clean_tmass_psc_xsc_oid
    JOIN gaiadr1.tmass_original_valid AS tm ON tmj.original_psc_source_id = tm.designation
    """


def _phot_quality_where_sql() -> str:
    return """
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
    """


def _fetch_live_two_stage(
    *,
    base_where_gaia: str,
    reg: str,
    n_lim: int,
    after_random_index: int | None,
    cfg: PipelineConfig,
    cache: RawCache,
) -> tuple[list[dict[str, Any]], list[ManifestRecord]]:
    """Two-stage Gaia-index -> AllWISE/2MASS batch join for paginated live queries."""
    idx_query = " ".join(
        f"SELECT g.source_id, g.random_index FROM gaiadr3.gaia_source AS g WHERE {base_where_gaia}".split()
    )
    idx_rows, idx_man, _ = _tap_query_cached(
        urls=GAIA_TAP_URLS,
        query=idx_query,
        archive="gaia",
        release=f"{cfg.gaia.release}_{reg}_index",
        cache=cache,
        timeout_s=90,
    )
    ordered = sorted(idx_rows, key=lambda r: int(r["random_index"]))
    if after_random_index is not None:
        ordered = [r for r in ordered if int(r["random_index"]) > int(after_random_index)]

    manifests: list[ManifestRecord] = [idx_man]
    matched_rows: list[dict[str, Any]] = []
    pos = 0
    while pos < len(ordered) and len(matched_rows) < n_lim:
        remaining_needed = n_lim - len(matched_rows)
        chunk_size = min(500, max(100, int(remaining_needed * 1.8)))
        chunk = ordered[pos : pos + chunk_size]
        pos += len(chunk)
        id_csv = ",".join(str(int(r["source_id"])) for r in chunk)
        batch_query = " ".join(
            f"""
            SELECT {_select_columns_sql()}
            {_join_and_phot_cuts_sql()}
            WHERE g.source_id IN ({id_csv})
            {_phot_quality_where_sql()}
            ORDER BY g.random_index ASC
            """.split()
        )
        b_rows, b_man, _ = _tap_query_cached(
            urls=GAIA_TAP_URLS,
            query=batch_query,
            archive="gaia",
            release=f"{cfg.gaia.release}_{reg}_batch",
            cache=cache,
            timeout_s=60,
        )
        manifests.append(b_man)
        matched_rows.extend(b_rows)

    matched_rows.sort(key=lambda r: int(r["random_index"]))
    return matched_rows[:n_lim], manifests


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
        base_where_gaia = """
          g.parallax_over_error > 10.0
          AND g.parallax BETWEEN 0.5 AND 20.0
          AND g.ruwe < 1.4
          AND ABS(g.b) > 20.0
          AND g.phot_g_mean_mag BETWEEN 7.5 AND 14.5
          AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
          AND g.mh_gspphot < -0.8
          AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 150.0
        """
        where_region = f"{base_where_gaia} {extra_cursor}"
        order_clause = "ORDER BY g.random_index ASC"
    elif reg == "bulge":
        # Galactic Bulge window around l ~ 0 deg, |b| 3..10 deg, bulge distances (~2.2-12.5 kpc)
        base_where_gaia = """
          (g.l <= 12.0 OR g.l >= 348.0)
          AND ABS(g.b) BETWEEN 3.0 AND 10.0
          AND g.parallax BETWEEN 0.08 AND 0.45
          AND g.parallax_over_error > 5.0
          AND g.ruwe < 1.4
          AND g.phot_g_mean_mag BETWEEN 8.0 AND 14.5
          AND g.teff_gspphot BETWEEN 3500.0 AND 7500.0
        """
        where_region = f"{base_where_gaia} {extra_cursor}"
        order_clause = "ORDER BY g.random_index ASC"
    elif reg == "local":
        # Legacy local high-velocity cut (v_tan > 65 km/s, |b| > 12 deg, 40-500 pc)
        base_where_gaia = """
          g.parallax_over_error > 15.0
          AND g.parallax BETWEEN 2.0 AND 25.0
          AND g.ruwe < 1.4
          AND ABS(g.b) > 12.0
          AND g.phot_g_mean_mag BETWEEN 7.5 AND 13.5
          AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
          AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 65.0
        """
        where_region = base_where_gaia
        order_clause = "ORDER BY g.source_id ASC"
    else:
        raise ArchiveAccessError(f"Unsupported live region '{region}'; expected 'halo', 'bulge', or 'local'", retryable=False)

    # If a two-stage Gaia index is already cached for paginated queries, use it directly without a 60s sync timeout.
    if reg in {"halo", "bulge"} and after_random_index is not None:
        idx_query = " ".join(
            f"SELECT g.source_id, g.random_index FROM gaiadr3.gaia_source AS g WHERE {base_where_gaia}".split()
        )
        idx_sig = hashlib.sha256(f"{'|'.join(GAIA_TAP_URLS)}|{idx_query}".encode("utf-8")).hexdigest()[:24]
        cached_idx = cache.lookup("gaia", f"{cfg.gaia.release}_{reg}_index", idx_sig)
        if cached_idx is not None and not cached_idx.offline_fixture:
            rows, mans = _fetch_live_two_stage(
                base_where_gaia=base_where_gaia,
                reg=reg,
                n_lim=n_lim,
                after_random_index=after_random_index,
                cfg=cfg,
                cache=cache,
            )
            stars = [_row_to_star_dict(r, gaia_release=cfg.gaia.release) for r in rows]
            return stars, mans

    query = f"""
    SELECT TOP {n_lim}
      {_select_columns_sql()}
    {_join_and_phot_cuts_sql()}
    WHERE {where_region}
    {_phot_quality_where_sql()}
    {order_clause}
    """
    try:
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
    except ArchiveAccessError:
        if reg not in {"halo", "bulge"}:
            raise
        rows, mans = _fetch_live_two_stage(
            base_where_gaia=base_where_gaia,
            reg=reg,
            n_lim=n_lim,
            after_random_index=after_random_index,
            cfg=cfg,
            cache=cache,
        )
        stars = [_row_to_star_dict(r, gaia_release=cfg.gaia.release) for r in rows]
        return stars, mans


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
      awn.number_of_neighbours AS allwise_n_neighbours,
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


def _enrich_optical_timeseries(out: dict[str, Any], cache: RawCache, cfg: PipelineConfig | None = None) -> list[ManifestRecord]:
    """Attach the vetting-only optical time series: TESS QLP/TESScut pairs (with background and centroid series),
    Gaia neighbours inside the TESS aperture, SkyBoT queries at flagged TESS event times, and ZTF + ASAS-SN
    seasonal photometry.

    Every fetch records its outcome so checks can tell no coverage from a failed query: survey summaries carry
    `status` in {ok, no_coverage, query_failed, client_unavailable} plus `query_log`; TESS sets
    `lightcurves.tess_status`/`tess_query_log`; the Gaia cone sets `lightcurves.tess_neighbours.status`; SkyBoT
    sets `lightcurves.sso_queries.status`. Only archive/network/parse errors are caught (by type, recorded);
    programming errors propagate.
    """
    import zipfile

    from oddstar.io.optical_timeseries import query_sso_at_events

    # Archive, network and payload-parsing failures. AttributeError/TypeError cover malformed or stubbed HTTP responses
    # and third-party client payloads; they are recorded as query_failed with their type in the note, never silent.
    _OPTICAL_FETCH_ERRORS = (ArchiveAccessError, requests.RequestException, ValueError, KeyError, OSError, zipfile.BadZipFile,
                             StopIteration, AttributeError, TypeError)

    ra = float(out["ra_deg"])
    dec = float(out["dec_deg"])
    pmra = float(out.get("pmra_mas_yr", 0.0) or 0.0)
    pmdec = float(out.get("pmdec_mas_yr", 0.0) or 0.0)
    sid = str(out["source_id"])
    lc = out.setdefault("lightcurves", {})
    surveys: dict[str, Any] = {}
    manifests: list[ManifestRecord] = []

    def _failed(survey: str, exc: BaseException) -> dict[str, Any]:
        return {"survey": survey, "n_epochs": 0, "status": "query_failed", "query_log": [f"{type(exc).__name__}: {str(exc)[:240]}"]}

    try:
        ztf, m = fetch_ztf_lightcurve(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec, cache=cache)
        manifests.extend(m)
        surveys["ztf"] = ztf if ztf is not None else {"survey": "ZTF", "n_epochs": 0, "status": "query_failed", "query_log": ["no summary returned"]}
    except _OPTICAL_FETCH_ERRORS as exc:
        surveys["ztf"] = _failed("ZTF", exc)
    try:
        asn, m = fetch_asassn_lightcurve(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec, cache=cache)
        manifests.extend(m)
        surveys["asassn"] = asn if asn is not None else {
            "survey": "ASAS-SN", "n_epochs": 0, "status": "client_unavailable",
            "query_log": ["pyasassn (pip package 'skypatrol') not installed"],
        }
    except _OPTICAL_FETCH_ERRORS as exc:
        surveys["asassn"] = _failed("ASAS-SN", exc)
    lc["optical_surveys"] = surveys

    # Gaia DR3 sources inside the TESS aperture region (vetting.tess_neighbour_radius_px TESS pixels) for neighbour_signal
    neighbours: list[dict[str, Any]] = []
    if cfg is not None:
        radius_arcsec = float(cfg.vetting["tess_neighbour_radius_px"]) * float(cfg.vetting["tess_pixel_arcsec"])
        nb_info: dict[str, Any] = {"radius_arcsec": radius_arcsec, "sources": [], "status": "ok", "query_log": []}
        try:
            q_tn = (
                f"SELECT source_id, ra, dec, pmra, pmdec, phot_g_mean_mag, phot_rp_mean_mag "
                f"FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, {radius_arcsec / 3600.0:.6f}))"
            )
            tn_rows, m_tn, _ = _tap_query_cached(urls=GAIA_TAP_URLS, query=q_tn, archive="gaia", release=cfg.gaia.release, cache=cache, timeout_s=25)
            manifests.append(m_tn)
            t_rp = next((_safe_float(r.get("phot_rp_mean_mag")) for r in tn_rows if str(r.get("source_id")) == sid), None)
            if t_rp is None:
                t_rp = _safe_float((out.get("mags") or {}).get("RP"))
            for r in tn_rows:
                if str(r.get("source_id")) == sid:
                    continue
                n_ra, n_dec, rp = _safe_float(r.get("ra")), _safe_float(r.get("dec")), _safe_float(r.get("phot_rp_mean_mag"))
                if n_ra is None or n_dec is None:
                    continue
                ratio = None if (rp is None or t_rp is None) else round(10 ** (-0.4 * (rp - t_rp)), 5)
                neighbours.append({"id": f"Gaia DR3 {r['source_id']}", "ra_deg": n_ra, "dec_deg": n_dec,
                                   "pmra": _safe_float(r.get("pmra")) or 0.0, "pmdec": _safe_float(r.get("pmdec")) or 0.0,
                                   "sep_arcsec": round(angular_separation_arcsec(ra, dec, n_ra, n_dec), 2),
                                   "g_mag": _safe_float(r.get("phot_g_mean_mag")), "rp_mag": rp, "flux_ratio": ratio})
            nb_info["sources"] = sorted(neighbours, key=lambda x: x["sep_arcsec"])
            nb_info["query_log"].append(f"ok ({len(neighbours)} sources)")
        except _OPTICAL_FETCH_ERRORS as exc:
            nb_info.update({"status": "query_failed", "query_log": [f"{type(exc).__name__}: {str(exc)[:240]}"]})
        lc["tess_neighbours"] = nb_info

    try:
        tess, m = fetch_tess_pipeline_pair(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec, gaia_id=sid, cache=cache, neighbours=neighbours)
        manifests.extend(m)
        if tess is not None:
            lc["tess_tic_id"] = tess.get("tic_id")
            lc["tess_sectors"] = tess.get("sectors", [])
            lc["tess_query_log"] = tess.get("query_log", [])
            lc["tess_status"] = tess.get("status", "ok")
    except _OPTICAL_FETCH_ERRORS as exc:
        lc["tess_status"] = "query_failed"
        lc["tess_query_log"] = [f"{type(exc).__name__}: {str(exc)[:240]}"]

    # SkyBoT at the flagged TESS event times (deepest coincident bin per sector above optical_dips.min_dip_depth_frac)
    if cfg is not None and lc.get("tess_sectors"):
        from oddstar.vet.lc_checks import tess_dip_events

        events, _ = tess_dip_events(lc, cfg)
        if events:
            radius = float(cfg.vetting["sso_tess_radius_px"]) * float(cfg.vetting["tess_pixel_arcsec"])
            sso, m = query_sso_at_events(events, ra=ra, dec=dec, pmra=pmra, pmdec=pmdec, radius_arcsec=radius, cache=cache)
            manifests.extend(m)
            sts = {e["status"] for e in sso["events"]}
            sso["status"] = "ok" if sts == {"ok"} else ("client_unavailable" if "client_unavailable" in sts else "query_failed")
            lc["sso_queries"] = sso
        else:
            lc["sso_queries"] = {"status": "no_event", "events": []}
    return manifests


def _failed(exc: BaseException) -> str:
    """Query-status string for a failed archive query."""
    msg = " ".join(str(exc).split())
    return f"failed: {type(exc).__name__}: {msg[:300]}"


# SIMBAD object types (basic.otype short codes and otype_txt long forms) that are youth indicators by themselves
YOUNG_STAR_OTYPES = {"Y*O", "YSO", "Y*?", "YSO_Candidate", "TT*", "TTauri*", "TT?", "TTauri*_Candidate", "Or*", "Orion_V*", "Ae*", "Ae*_Candidate"}
# SIMBAD parent (h_link) object types of young stellar groups vs clusters of unknown age
YOUNG_GROUP_OTYPES = {"As*", "Association", "MGr", "MouvGroup", "SFR", "StarFormingRegion"}
CLUSTER_OTYPES = {"OpC", "OpenCluster", "Cl*", "Cluster*", "C?*", "ClG", "GlC", "GlobCluster"}


def _vizier_value(row: dict[str, Any], *keys: str) -> float | None:
    for k in keys:
        if k in row:
            v = _safe_float(row.get(k))
            if v is not None:
                return v
    return None


def _vizier_radec(row: dict[str, Any]) -> tuple[float | None, float | None]:
    ra = _vizier_value(row, "RAJ2000", "RA_ICRS", "_RAJ2000", "RAdeg", "RA")
    de = _vizier_value(row, "DEJ2000", "DE_ICRS", "_DEJ2000", "DEdeg", "DE", "Dec")
    return ra, de


def _neighbour_ir_flux_estimate(
    g_mag: Any, rp_mag: Any, target_w3_mjy: float | None, target_w4_mjy: float | None, vcfg: dict[str, float]
) -> dict[str, Any]:
    """Conservative W3/W4 photospheric flux of a Gaia neighbour relative to the target's observed W3/W4 flux.

    Colour assumption (configs `neighbour_assumed_*`): the neighbour is given the colours of a mid/late-M dwarf,
    RP - W3 = 3.7 and RP - W4 = 4.0 (Vega), with G - RP = 1.3 when RP is missing. These values are approximate
    (not re-derived from a published table in this change) and were chosen to be redder than any normal FGK or
    early-M star, so for a *stellar* neighbour the estimate over-predicts its mid-IR flux. It does not bound a
    galaxy or a dusty neighbour; those are screened by `deep_imaging`/`vision_rubric`. A plain Rayleigh-Jeans
    extrapolation from RP is NOT conservative for cool stars (it under-predicts a 3000 K star's 22 um flux by
    ~70x), which is why a fixed red colour is used instead.
    """
    g_f, rp_f = _safe_float(g_mag), _safe_float(rp_mag)
    if rp_f is None and g_f is None:
        return {"flux_frac": None, "flux_frac_basis": "no Gaia G/RP photometry"}
    basis = "RP" if rp_f is not None else "G"
    rp_eff = rp_f if rp_f is not None else float(g_f) - float(vcfg["neighbour_assumed_g_minus_rp"])
    f3 = mag_to_flux_mjy("W3", rp_eff - float(vcfg["neighbour_assumed_rp_minus_w3"]))
    f4 = mag_to_flux_mjy("W4", rp_eff - float(vcfg["neighbour_assumed_rp_minus_w4"]))
    fr3 = f3 / target_w3_mjy if target_w3_mjy and target_w3_mjy > 0 else None
    fr4 = f4 / target_w4_mjy if target_w4_mjy and target_w4_mjy > 0 else None
    fracs = [f for f in (fr3, fr4) if f is not None]
    return {
        "w3_est_mjy": round(f3, 6),
        "w4_est_mjy": round(f4, 6),
        "flux_frac_w3": None if fr3 is None else round(fr3, 6),
        "flux_frac_w4": None if fr4 is None else round(fr4, 6),
        "flux_frac": max(fracs) if fracs else None,
        "flux_frac_basis": f"{basis} with assumed M-dwarf colours",
    }


def _group_visits(mjds: list[float], gap_days: float) -> list[list[float]]:
    visits: list[list[float]] = []
    for m in sorted(mjds):
        if not visits or m - visits[-1][-1] > gap_days:
            visits.append([m])
        else:
            visits[-1].append(m)
    return visits


def _search_sso(
    *,
    ra: float,
    dec: float,
    pmra: float,
    pmdec: float,
    frame_mjds: list[float],
    vcfg: dict[str, float],
    cache: RawCache,
    sleep: Any,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """SkyBoT (WISE observer C51) at the median epoch of each 2010 visit, propagated to every frame epoch."""
    match_r = float(vcfg["sso_match_radius_arcsec"])
    max_rate = float(vcfg["sso_max_rate_arcsec_per_hr"])
    visits_out: list[dict[str, Any]] = []
    matches: list[dict[str, Any]] = []
    manifests: list[ManifestRecord] = []
    for frames in _group_visits(frame_mjds, float(vcfg["sso_visit_gap_days"])):
        ep_mjd = float(np.median(frames))
        half_span_h = max(abs(f - ep_mjd) for f in frames) * 24.0
        radius = match_r + max_rate * max(0.5, half_span_h)
        c_ra, c_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2000.0 + (ep_mjd - 51544.5) / 365.25)
        visit: dict[str, Any] = {
            "epoch_mjd": ep_mjd,
            "n_frames": len(frames),
            "frame_mjd_range": [min(frames), max(frames)],
            "query_radius_arcsec": round(radius, 1),
        }
        try:
            bodies, man = skybot_cone_search(
                ra_deg=c_ra, dec_deg=c_dec, radius_arcsec=radius, epoch_jd=ep_mjd + 2400000.5,
                observer="C51", cache=cache, sleep=sleep,
            )
            manifests.append(man)
        except Exception as exc:
            visit["status"] = _failed(exc)
            visits_out.append(visit)
            continue
        visit["status"] = "ok"
        visit["n_bodies"] = len(bodies)
        ambiguous: list[str] = []
        closest: float | None = None
        for b in bodies:
            dra, ddec = b.get("dra_arcsec_hr"), b.get("ddec_arcsec_hr")
            err = float(b.get("pos_err_arcsec") or 0.0)
            for f in frames:
                dt_h = (f - ep_mjd) * 24.0
                s_ra, s_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2000.0 + (f - 51544.5) / 365.25)
                if dra is None or ddec is None:
                    # Motion unknown: only decidable at the visit epoch itself
                    sep = angular_separation_arcsec(s_ra, s_dec, b["ra_deg"], b["dec_deg"])
                    if sep - max_rate * abs(dt_h) <= match_r + err and abs(dt_h) > 1e-3:
                        ambiguous.append(b["name"])
                        continue
                    b_ra, b_dec = b["ra_deg"], b["dec_deg"]
                else:
                    cosd = max(1e-6, math.cos(math.radians(b["dec_deg"])))
                    b_ra = b["ra_deg"] + float(dra) * dt_h / 3600.0 / cosd
                    b_dec = b["dec_deg"] + float(ddec) * dt_h / 3600.0
                sep = angular_separation_arcsec(s_ra, s_dec, b_ra, b_dec)
                closest = sep if closest is None else min(closest, sep)
                if sep <= match_r + err:
                    matches.append({"name": b["name"], "class": b.get("class"), "v_mag": b.get("v_mag"),
                                    "frame_mjd": f, "sep_arcsec": round(sep, 2), "pos_err_arcsec": err})
        if ambiguous:
            visit["status"] = f"ambiguous: no rates for {sorted(set(ambiguous))[:5]}"
        visit["closest_approach_arcsec"] = None if closest is None else round(closest, 1)
        visits_out.append(visit)
    return {"frame_mjds": sorted(frame_mjds), "visits": visits_out, "matches": matches,
            "observer": "C51", "match_radius_arcsec": match_r}, manifests


def enrich_star_live_for_vetting(
    star: dict[str, Any],
    cfg: PipelineConfig,
    cache: RawCache,
    *,
    enrich_optical: bool = True,
    sleep: Any = time.sleep,
) -> tuple[dict[str, Any], list[ManifestRecord]]:
    """Enrich one real star using primary + backup mirrors across Gaia, SIMBAD, IRSA, MAST, VizieR, IMCCE and Data Lab.

    Every query records its outcome in `vetting_archive.query_status[name]` as "ok", "no_coverage" or
    "failed: <ExceptionType>: <msg>". Vetting inputs are written only from successful queries, so a failed
    query leaves the input absent and the check reports `not_evaluated` naming the query.
    """
    if star.get("live_enriched"):
        return star, []

    ra = float(star["ra_deg"])
    dec = float(star["dec_deg"])
    pmra = float(star.get("pmra_mas_yr", 0.0) or 0.0)
    pmdec = float(star.get("pmdec_mas_yr", 0.0) or 0.0)
    sid = str(star["source_id"])
    manifests: list[ManifestRecord] = []
    out = json.loads(json.dumps(star))
    va: dict[str, Any] = out.setdefault("vetting_archive", {})
    # Legacy seeing-limited "deep imaging" verdicts are superseded by `neighbour_screen` + `subarcsec_imaging`.
    va.pop("deep_imaging_available", None)
    va.pop("deep_imaging_offset_galaxy", None)
    out.setdefault("allwise", {})
    qs: dict[str, str] = va.setdefault("query_status", {})
    vcfg = cfg.vetting
    flux = out.get("flux_mjy") or {}

    # 1. Gaia DR3 neighbours inside the 12" W4 beam with estimated W3/W4 flux fractions
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
            est = _neighbour_ir_flux_estimate(nr.get("phot_g_mean_mag"), nr.get("phot_rp_mean_mag"), flux.get("W3"), flux.get("W4"), vcfg)
            nbs.append({"id": f"Gaia DR3 {nr['source_id']}", "sep_arcsec": round(sep, 2), "g_mag": nr.get("phot_g_mean_mag"),
                        "rp_mag": nr.get("phot_rp_mean_mag"), **est})
        va["neighbours_in_w4_beam"] = nbs
        qs["gaia_neighbours"] = "ok"
    except Exception as exc:
        qs["gaia_neighbours"] = _failed(exc)

    # 1b. Gaia RV / NSS / variability metadata
    try:
        q_meta = (
            f"SELECT non_single_star, radial_velocity_error, rv_nb_transits, rv_method_used, "
            f"rv_chisq_pvalue, rv_amplitude_robust, phot_variable_flag "
            f"FROM gaiadr3.gaia_source WHERE source_id = {sid}"
        )
        m_rows, m_meta, _ = _tap_query_cached(urls=GAIA_TAP_URLS, query=q_meta, archive="gaia", release=cfg.gaia.release, cache=cache, timeout_s=20)
        manifests.append(m_meta)
        if not m_rows:
            qs["gaia_meta"] = "no_coverage"
        else:
            mr = m_rows[0]
            if mr.get("non_single_star") is not None:
                va["nss_flag"] = int(mr["non_single_star"])
            rv_amp = _safe_float(mr.get("rv_amplitude_robust"))
            rv_err = _safe_float(mr.get("radial_velocity_error"))
            rv_n = int(mr.get("rv_nb_transits") or 0)
            rv_method = int(mr.get("rv_method_used") or 0)
            rv_p = _safe_float(mr.get("rv_chisq_pvalue"))
            if rv_amp is not None:
                va["rv_scatter_kms"] = rv_amp
                va["rv_faint_unconstrained"] = False
            elif rv_err is not None and rv_n >= 2:
                # radial_velocity_error is uncertainty on the median RV; scale to epoch scatter ~ sqrt(2N/pi)*err
                va["rv_scatter_kms"] = rv_err * math.sqrt((2.0 * rv_n) / math.pi)
                va["rv_faint_unconstrained"] = bool(rv_method == 2 and rv_p is None)
            elif rv_err is not None:
                va["rv_scatter_kms"] = rv_err
            if rv_p is not None:
                va["rv_chisq_pvalue"] = rv_p
            if mr.get("phot_variable_flag"):
                va["gaia_phot_variable_flag"] = str(mr["phot_variable_flag"])
            qs["gaia_meta"] = "ok"
    except Exception as exc:
        qs["gaia_meta"] = _failed(exc)

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
        qs["simbad"] = "ok"
    except Exception as exc:
        qs["simbad"] = _failed(exc)

    # 3a. AllWISE coadd quality (IRSA Gator fast path -> IRSA TAP -> NOIRLab Data Lab backup)
    aw_rows: list[dict[str, Any]] = []
    aw_errors: list[str] = []
    aw_ok = False
    try:
        aw_rows, m_aw, _ = _gator_cone_cached(
            catalog="allwise_p3as_psd", ra=ra, dec=dec, radius_arcsec=5.0, archive="allwise", release="allwise", cache=cache, timeout_s=20
        )
        manifests.append(m_aw)
        aw_ok = True
    except Exception as exc:
        aw_errors.append(f"gator {_failed(exc)}")
        try:
            q_aw = (
                f"SELECT ra, dec, w4rchi2, w4m, w3rchi2, cc_flags, ext_flg, ph_qual "
                f"FROM allwise_p3as_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.001))"
            )
            aw_rows, m_aw, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_aw, archive="allwise", release="allwise", cache=cache, timeout_s=25)
            manifests.append(m_aw)
            aw_ok = True
        except Exception as exc2:
            aw_errors.append(f"irsa_tap {_failed(exc2)}")
            try:
                q_aw_dl = (
                    f"SELECT ra, dec, w4rchi2, w4m, w4nm, w4snr, cc_flags, ext_flg, ph_qual "
                    f"FROM allwise.source WHERE q3c_radial_query(ra, dec, {ra:.6f}, {dec:.6f}, 0.0015)"
                )
                aw_rows, m_aw, _ = _tap_query_cached(url=NOIRLAB_TAP_URL, query=q_aw_dl, archive="allwise", release="datalab", cache=cache, timeout_s=20)
                manifests.append(m_aw)
                aw_ok = True
            except Exception as exc3:
                aw_errors.append(f"datalab {_failed(exc3)}")
                aw_rows = []
    aw_rows = [r for r in aw_rows if _safe_float(r.get("ra")) is not None and _safe_float(r.get("dec")) is not None]
    if not aw_ok:
        qs["allwise_coadd"] = "failed: " + " | ".join(aw_errors)
    elif not aw_rows:
        qs["allwise_coadd"] = "no_coverage"
    else:
        qs["allwise_coadd"] = "ok"
        best = min(aw_rows, key=lambda r: angular_separation_arcsec(ra, dec, float(r["ra"]), float(r["dec"])))
        for key, cast in (("w4rchi2", float), ("w4m", int), ("w4nm", int), ("w3nm", int), ("w4snr", float), ("w3snr", float)):
            val = _safe_float(best.get(key))
            if val is not None:
                out["allwise"][key] = cast(val)
        if best.get("cc_flags"):
            out["allwise"]["cc_flags"] = str(best["cc_flags"])
        if best.get("ph_qual"):
            out["allwise"]["ph_qual"] = str(best["ph_qual"])
        if best.get("ext_flg") is not None:
            out["allwise"]["ext_flg"] = int(float(best["ext_flg"]))

    # 3b. WISE 2010 4-band Single-Exposure (L1b) frames (IRSA Gator fast path -> IRSA TAP)
    l1b_rows: list[dict[str, Any]] = []
    l1b_errors: list[str] = []
    l1b_ok = False
    try:
        l1b_rows, m_l1b, _ = _gator_cone_cached(
            catalog="allsky_4band_p1bs_psd", ra=ra, dec=dec, radius_arcsec=5.0, archive="allwise", release="allsky_4band_l1b", cache=cache, timeout_s=25
        )
        manifests.append(m_l1b)
        l1b_ok = True
    except Exception as exc:
        l1b_errors.append(f"gator {_failed(exc)}")
        try:
            q_l1b = (
                f"SELECT mjd, w4mpro, w4sigmpro, w4rchi2, cc_flags, qual_frame "
                f"FROM allsky_4band_p1bs_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0014))"
            )
            l1b_rows, m_l1b, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_l1b, archive="allwise", release="allsky_4band_l1b", cache=cache, timeout_s=25)
            manifests.append(m_l1b)
            l1b_ok = True
        except Exception as exc2:
            l1b_errors.append(f"irsa_tap {_failed(exc2)}")
            l1b_rows = []
    frame_mjds: list[float] = []
    if not l1b_ok:
        qs["wise_l1b"] = "failed: " + " | ".join(l1b_errors)
    else:
        w4_frames_mjy: list[float] = []
        for lr in l1b_rows:
            mjd_v = _safe_float(lr.get("mjd"))
            if mjd_v is not None:
                frame_mjds.append(mjd_v)
            w4_mag = _safe_float(lr.get("w4mpro"))
            w4_sig = _safe_float(lr.get("w4sigmpro"))
            w4_snr = _safe_float(lr.get("w4snr"))
            if w4_mag is not None and w4_sig is not None and (w4_snr is None or w4_snr >= 2.0):
                w4_frames_mjy.append(float(mag_to_flux_mjy("W4", w4_mag)))
            elif w4_mag is not None:
                # Upper-limit / non-detection frame contributes 0.0 detected flux to single-frame consistency
                w4_frames_mjy.append(0.0)
        qs["wise_l1b"] = "ok" if l1b_rows else "no_coverage"
        if w4_frames_mjy:
            out["wise_single_frames_w4_mjy"] = w4_frames_mjy
        va["wise_l1b_frame_mjds"] = sorted(frame_mjds)

    # 3c. NEOWISE-R multi-epoch W1/W2 light curve (IRSA Gator fast path -> IRSA TAP)
    neo_rows: list[dict[str, Any]] = []
    neo_primary_err: str | None = None
    try:
        try:
            neo_rows, m_neo, _ = _gator_cone_cached(
                catalog="neowiser_p1bs_psd", ra=ra, dec=dec, radius_arcsec=3.0, archive="neowise", release="neowise_r", cache=cache, timeout_s=25
            )
            manifests.append(m_neo)
        except Exception as exc_primary:
            neo_primary_err = _failed(exc_primary)  # kept so a backup failure reports both
            q_neo = (
                f"SELECT mjd, w1mpro, w1sigmpro, w2mpro, w2sigmpro, qual_frame, qi_fact "
                f"FROM neowiser_p1bs_psd WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.00083)) "
                f"AND qual_frame > 0"
            )
            neo_rows, m_neo, _ = _tap_query_cached(url=IRSA_TAP_URL, query=q_neo, archive="neowise", release="neowise_r", cache=cache, timeout_s=30)
            manifests.append(m_neo)
        qs["neowise"] = "ok" if neo_rows else "no_coverage"
    except Exception as exc:
        qs["neowise"] = _failed(exc) + (f" (gator primary {neo_primary_err})" if neo_primary_err else "")
        neo_rows = []

    if neo_rows:
        pts: list[tuple[float, float, float, float | None, float | None]] = []
        for row in neo_rows:
            qf = _safe_float(row.get("qual_frame"))
            mjd_val = _safe_float(row.get("mjd"))
            w1_val = _safe_float(row.get("w1mpro"))
            w2_val = _safe_float(row.get("w2mpro"))
            r_ra = _safe_float(row.get("ra"))
            r_dec = _safe_float(row.get("dec"))
            if (qf is None or qf > 0) and mjd_val is not None and w1_val is not None and w2_val is not None:
                pts.append((mjd_val, w1_val, w2_val, r_ra, r_dec))
        pts.sort(key=lambda x: x[0])
        if pts:
            # Bin into ~180-day semi-annual visits and check PM-corrected astrometric trajectory
            bins: list[tuple[float, list[tuple[float, float, float, float | None, float | None]]]] = []
            for item in pts:
                if not bins or item[0] - bins[-1][0] > 120.0:
                    bins.append((item[0], [item]))
                else:
                    bins[-1][1].append(item)
            v_mjd, v_w1, v_w2, v_e1, v_e2, pm_resids = [], [], [], [], [], []
            v_dra: list[float | None] = []
            v_ddec: list[float | None] = []
            v_win: list[tuple[float, float]] = []
            v_nfr: list[int] = []
            for _, grp in bins:
                mjd_med = float(np.median([x[0] for x in grp]))
                v_mjd.append(mjd_med)
                v_win.append((float(grp[0][0]), float(grp[-1][0])))
                v_nfr.append(len(grp))
                f1s = [mag_to_flux_mjy("W1", x[1]) for x in grp]
                f2s = [mag_to_flux_mjy("W2", x[2]) for x in grp]
                m1 = float(np.median(f1s))
                m2 = float(np.median(f2s))
                v_w1.append(m1)
                v_w2.append(m2)
                v_e1.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f1s) - m1))) / math.sqrt(len(grp)), m1 * 0.02))
                v_e2.append(max(0.05, float(1.4826 * np.median(np.abs(np.asarray(f2s) - m2))) / math.sqrt(len(grp)), m2 * 0.02))
                ras_g = [x[3] for x in grp if x[3] is not None]
                decs_g = [x[4] for x in grp if x[4] is not None]
                if ras_g and decs_g:
                    jyr = 2000.0 + (mjd_med - 51544.5) / 365.25
                    pra, pdec = propagate_position(ra, dec, pmra, pmdec, 2016.0, jyr)
                    r_med, d_med = float(np.median(ras_g)), float(np.median(decs_g))
                    pm_resids.append(angular_separation_arcsec(r_med, d_med, pra, pdec))
                    # Per-visit residual vector from the Gaia PM track (RA*cos(Dec), Dec; arcsec) for neighbour_signal
                    v_dra.append(round(((r_med - pra + 540.0) % 360.0 - 180.0) * math.cos(math.radians(pdec)) * 3600.0, 5))
                    v_ddec.append(round((d_med - pdec) * 3600.0, 5))
                else:
                    v_dra.append(None)
                    v_ddec.append(None)
            neo_dict: dict[str, Any] = {"mjd": v_mjd, "w1_mjy": v_w1, "w2_mjy": v_w2, "w1_err_mjy": v_e1, "w2_err_mjy": v_e2}
            if pm_resids:
                neo_dict["median_pm_residual_arcsec"] = float(np.median(pm_resids))
                neo_dict["pm_residual_ra_arcsec"] = v_dra
                neo_dict["pm_residual_dec_arcsec"] = v_ddec
            neo_dict["visit_n_frames"] = v_nfr
            # Field-star ensemble on the same visit grid for systematics_window (vetting-only; channel inputs unchanged)
            from oddstar.io.optical_timeseries import build_neowise_comparison

            itd = cfg.channels["ir_timedomain"]
            comp_r = float(itd["neowise_comparison_radius_arcsec"])
            try:
                c_rows, m_cmp, _ = _gator_cone_cached(
                    catalog="neowiser_p1bs_psd", ra=ra, dec=dec, radius_arcsec=comp_r,
                    archive="neowise", release="neowise_r_field", cache=cache, timeout_s=120,
                )
                manifests.append(m_cmp)
                comp = build_neowise_comparison(
                    c_rows, ra=ra, dec=dec, visit_windows=v_win,
                    target_w1_mag=float(np.median([x[1] for x in pts])),
                    exclude_arcsec=float(itd["neowise_comparison_exclude_arcsec"]),
                    mag_window=float(itd["neowise_comparison_mag_window"]),
                    min_stars=int(itd["neowise_comparison_min_stars"]),
                    max_stars=int(itd["neowise_comparison_max_stars"]),
                )
                comp["query_log"] = [f"Gator neowiser_p1bs_psd {comp_r:.0f}\" cone: {len(c_rows)} rows"]
            except (ArchiveAccessError, requests.RequestException, ValueError, KeyError, OSError) as exc:
                comp = {"status": "query_failed", "n_stars": 0, "query_log": [f"{type(exc).__name__}: {str(exc)[:240]}"]}
            comp["radius_arcsec"] = comp_r
            neo_dict["comparison"] = comp
            out["neowise"] = neo_dict

    # 4. Second-instrument 18-25 um archival check (Spitzer, AKARI IRC/FIS). A catalogue row within 10" is
    #    recorded as an unverified match: the photometry is not compared with W4 here.
    sec_errors: list[str] = []
    for cat, instr_name in (("slphotdr4", "Spitzer_MIPS24"), ("akari_irc", "AKARI_IRC_18um"), ("akari_fis", "AKARI_FIS")):
        try:
            sec_rows, m_sec, _ = _gator_cone_cached(catalog=cat, ra=ra, dec=dec, radius_arcsec=10.0, archive=cat, release="live", cache=cache, timeout_s=15)
            manifests.append(m_sec)
            qs[f"second_instrument_{cat}"] = "ok"
            if sec_rows and va.get("second_instrument_ir") is None:
                va["second_instrument_ir"] = {
                    "instrument": instr_name,
                    "n_matches": len(sec_rows),
                    "interpretation": "catalogue_match_unverified",
                }
                break
        except Exception as exc:
            qs[f"second_instrument_{cat}"] = _failed(exc)
            sec_errors.append(f"{cat}: {_failed(exc)}")
    qs["second_instrument"] = "ok" if (va.get("second_instrument_ir") is not None or not sec_errors) else "failed: " + " | ".join(sec_errors)

    # 5. Radio continuum (NVSS via VizieR CfA/CDS TSV fast path -> VizieR TAP backup)
    rad_primary_err: str | None = None
    try:
        try:
            r_rows, m_rad, _ = _vizier_tsv_cone_cached(
                source="VIII/65/nvss", ra=ra, dec=dec, radius_arcsec=10.0, archive="radio_continuum", release="nvss", cache=cache, timeout_s=15
            )
        except Exception as exc_primary:
            rad_primary_err = _failed(exc_primary)  # kept so a backup failure reports both
            q_rad = (
                f"SELECT RAJ2000, DEJ2000, S1_4 AS flux_mjy FROM \"VIII/65/nvss\" "
                f"WHERE 1=CONTAINS(POINT('ICRS', RAJ2000, DEJ2000), CIRCLE('ICRS', {ra:.6f}, {dec:.6f}, 0.0028))"
            )
            r_rows, m_rad, _ = _tap_query_cached(urls=VIZIER_TAP_URLS, query=q_rad, archive="radio_continuum", release="nvss", cache=cache, timeout_s=15)
        manifests.append(m_rad)
        hits = []
        for rr in r_rows:
            r_ra = _safe_float(rr.get("RAJ2000"))
            r_dec = _safe_float(rr.get("DEJ2000"))
            if r_ra is None or r_dec is None:
                continue
            sep = angular_separation_arcsec(ra, dec, r_ra, r_dec)
            hits.append({"survey": "NVSS", "sep_arcsec": round(sep, 2), "flux_mjy": rr.get("S1_4") or rr.get("flux_mjy")})
        va["radio_sources_within_10arcsec"] = hits
        qs["radio_nvss"] = "ok"
    except Exception as exc:
        qs["radio_nvss"] = _failed(exc) + (f" (vizier_tsv primary {rad_primary_err})" if rad_primary_err else "")

    # 6. Seeing-limited neighbour screen (PS1 DR2 mean catalogue for Dec > -30, Legacy Survey DR10 Tractor).
    #    It can provide positive evidence of a contaminant but never clears `deep_imaging` (lesson L001).
    core_r = float(vcfg["deep_screen_core_radius_arcsec"])
    ext_min = float(vcfg["deep_screen_ps1_psf_minus_kron_min"])
    faint_r = float(vcfg["deep_screen_faint_r_mag"])
    flux_frac_min = float(vcfg["deep_screen_min_flux_frac_of_excess"])
    beam_w4 = float(vcfg["w4_beam_arcsec"])
    exc_w3 = target_excess_mjy(out, "W3")
    exc_w4 = target_excess_mjy(out, "W4")
    screen: dict[str, Any] = {
        "surveys": [],
        "contaminant": None,
        "extended_neighbours": [],
        "faint_point_sources": [],
        "target_excess_mjy": {"W3": exc_w3, "W4": exc_w4},
    }

    def _note_contaminant(item: dict[str, Any]) -> None:
        if screen["contaminant"] is None:
            screen["contaminant"] = item

    if dec > -30.0:
        ps1_done = False
        try:
            p_rows, m_ps1, _ = _mast_ps1_cone_cached(ra=ra, dec=dec, radius_deg=0.0028, cache=cache, timeout_s=20)
            manifests.append(m_ps1)
            # Filter out single-frame diffraction/ghost artifacts (require nDetections >= 2)
            valid_rows = [
                pr for pr in p_rows
                if _safe_float(pr.get("raMean")) is not None and _safe_float(pr.get("decMean")) is not None and int(pr.get("nDetections") or 0) >= 2
            ]
            qs["deep_optical_ps1"] = "ok" if valid_rows else "no_coverage"
            if valid_rows:
                ps1_done = True
                screen["surveys"].append("PS1 DR2 mean catalogue (MAST)")
                for pr in valid_rows:
                    sep = angular_separation_arcsec(ra, dec, float(pr["raMean"]), float(pr["decMean"]))
                    r_psf = _safe_float(pr.get("rMeanPSFMag"))
                    r_kron = _safe_float(pr.get("rMeanKronMag"))
                    diff = (r_psf - r_kron) if (r_psf is not None and r_kron is not None) else None
                    item = {"survey": "PS1", "sep_arcsec": round(sep, 2), "r_psf": r_psf, "r_kron": r_kron,
                            "psf_minus_kron": None if diff is None else round(diff, 3)}
                    if diff is not None and diff > ext_min:
                        if sep < core_r:
                            _note_contaminant({**item, "reason": f"extended (r_PSF - r_Kron > {ext_min}) inside {core_r}\""})
                        elif sep <= beam_w4:
                            screen["extended_neighbours"].append(item)
                    elif 1.8 <= sep <= 8.0 and (r_psf is None or r_psf > faint_r):
                        screen["faint_point_sources"].append(item)
        except Exception as exc:
            qs["deep_optical_ps1"] = _failed(exc)
        if not ps1_done:
            try:
                p_rows, m_ps1v, _ = _vizier_tsv_cone_cached(
                    source="II/349/ps1", ra=ra, dec=dec, radius_arcsec=10.0, archive="deep_optical", release="ps1_vizier", cache=cache, timeout_s=15
                )
                manifests.append(m_ps1v)
                rows_v = [pr for pr in p_rows if None not in _vizier_radec(pr)]
                qs["deep_optical_ps1_vizier"] = "ok" if rows_v else "no_coverage"
                if rows_v:
                    screen["surveys"].append("PS1 DR1 (VizieR II/349, no morphology)")
                for pr in rows_v:
                    p_ra, p_dec = _vizier_radec(pr)
                    sep = angular_separation_arcsec(ra, dec, float(p_ra), float(p_dec))
                    r_mag = _safe_float(pr.get("rmag"))
                    if 1.8 <= sep <= 8.0 and (r_mag is None or r_mag > faint_r):
                        screen["faint_point_sources"].append({"survey": "PS1-VizieR", "sep_arcsec": round(sep, 2), "r_mag": r_mag})
            except Exception as exc:
                qs["deep_optical_ps1_vizier"] = _failed(exc)
    else:
        qs["deep_optical_ps1"] = "no_coverage"

    try:
        cos_dec = max(0.1, math.cos(math.radians(dec)))
        dra = (12.0 / 3600.0) / cos_dec
        ddec = 12.0 / 3600.0
        q_ls = (
            f"SELECT ls_id, ra, dec, type, flux_g, flux_r, flux_z, flux_w3, flux_w4, flux_ivar_w4 "
            f"FROM ls_dr10.tractor "
            f"WHERE ra BETWEEN {ra - dra:.6f} AND {ra + dra:.6f} "
            f"AND dec BETWEEN {dec - ddec:.6f} AND {dec + ddec:.6f}"
        )
        ls_rows, m_ls, _ = _tap_csv_query_cached(url=NOIRLAB_TAP_URL, query=q_ls, archive="deep_optical", release="ls_dr10", cache=cache, timeout_s=25)
        manifests.append(m_ls)
        qs["deep_optical_ls_dr10"] = "ok" if ls_rows else "no_coverage"
        if ls_rows:
            screen["surveys"].append("Legacy Survey DR10 Tractor catalogue (with unWISE forced W3/W4)")
        for lr in ls_rows:
            l_ra = _safe_float(lr.get("ra"))
            l_dec = _safe_float(lr.get("dec"))
            if l_ra is None or l_dec is None:
                continue
            sep = angular_separation_arcsec(ra, dec, l_ra, l_dec)
            mtype = str(lr.get("type") or "")
            fr = _safe_float(lr.get("flux_r"))
            rmag = (22.5 - 2.5 * math.log10(fr)) if (fr is not None and fr > 0.0) else None
            # Tractor fluxes are AB nanomaggies: 1 nmgy = 3.631 uJy
            fw3 = (_safe_float(lr.get("flux_w3")) or 0.0) * 3.631e-3
            fw4 = (_safe_float(lr.get("flux_w4")) or 0.0) * 3.631e-3
            ivar_w4 = _safe_float(lr.get("flux_ivar_w4")) or 0.0
            if sep < core_r and ivar_w4 > 0.0:
                va["unwise_w4_snr"] = round((_safe_float(lr.get("flux_w4")) or 0.0) * math.sqrt(ivar_w4), 2)
            non_psf = mtype not in {"PSF", "DUP", ""}
            item = {"survey": "LS DR10", "sep_arcsec": round(sep, 2), "ls_type": mtype, "r_mag": None if rmag is None else round(rmag, 2),
                    "unwise_w3_mjy": round(fw3, 4), "unwise_w4_mjy": round(fw4, 4)}
            if non_psf and sep < core_r:
                _note_contaminant({**item, "reason": f"non-PSF Tractor source inside {core_r}\""})
            elif non_psf and sep <= beam_w4:
                comparable = (exc_w4 is not None and exc_w4 > 0 and fw4 >= flux_frac_min * exc_w4) or (
                    exc_w3 is not None and exc_w3 > 0 and fw3 >= flux_frac_min * exc_w3
                )
                if comparable:
                    _note_contaminant({**item, "reason": f"non-PSF source inside the W4 beam with unWISE flux >= {flux_frac_min:.0%} of the target excess"})
                else:
                    screen["extended_neighbours"].append({**item, "excess_known": exc_w4 is not None or exc_w3 is not None})
            elif not non_psf and 1.8 <= sep <= 8.0 and (rmag is None or rmag > faint_r):
                screen["faint_point_sources"].append(item)
    except Exception as exc:
        qs["deep_optical_ls_dr10"] = _failed(exc)
    va["neighbour_screen"] = screen

    # 6b. Sub-arcsecond imaging coverage (HST/JWST) from MAST CAOM
    try:
        r_sub = float(vcfg["subarcsec_imaging_search_radius_arcsec"])
        obs, m_caom = fetch_subarcsec_imaging(ra_deg=ra, dec_deg=dec, radius_arcsec=r_sub, cache=cache, sleep=sleep)
        manifests.append(m_caom)
        va["subarcsec_imaging"] = {"observations": obs[:25], "n_observations": len(obs), "search_radius_arcsec": r_sub,
                                   "collections": ["HST", "JWST"], "assessment": None}
        qs["subarcsec_imaging_mast"] = "ok"
    except Exception as exc:
        qs["subarcsec_imaging_mast"] = _failed(exc)

    # 7. Local surface density of W3/W4 sources for chance alignment (AllWISE annulus)
    try:
        r_in = float(vcfg["chance_alignment_annulus_inner_arcsec"])
        r_out = float(vcfg["chance_alignment_annulus_outer_arcsec"])
        min_snr = float(vcfg["chance_alignment_min_snr"])
        d_rows, m_den, _ = _gator_cone_cached(
            catalog="allwise_p3as_psd", ra=ra, dec=dec, radius_arcsec=r_out, archive="allwise", release="allwise_density", cache=cache,
            timeout_s=60, selcols="ra,dec,w3mpro,w3sigmpro,w3snr,w4mpro,w4sigmpro,w4snr",
        )
        manifests.append(m_den)
        srcs = []
        n_rows = 0
        for dr in d_rows:
            d_ra, d_dec = _safe_float(dr.get("ra")), _safe_float(dr.get("dec"))
            if d_ra is None or d_dec is None:
                continue
            n_rows += 1
            sep = angular_separation_arcsec(ra, dec, d_ra, d_dec)
            if sep < r_in or sep > r_out:
                continue
            s3, s4 = _safe_float(dr.get("w3snr")), _safe_float(dr.get("w4snr"))
            if (s3 is None or s3 < min_snr) and (s4 is None or s4 < min_snr):
                continue
            m3, m4 = _safe_float(dr.get("w3mpro")), _safe_float(dr.get("w4mpro"))
            srcs.append({"sep_arcsec": round(sep, 1), "w3snr": s3, "w4snr": s4,
                         "w3_mjy": None if m3 is None else round(mag_to_flux_mjy("W3", m3), 4),
                         "w4_mjy": None if m4 is None else round(mag_to_flux_mjy("W4", m4), 4)})
        if n_rows == 0:
            qs["allwise_density"] = "no_coverage"
        else:
            va["chance_alignment"] = {"catalogue": "AllWISE allwise_p3as_psd", "annulus_arcsec": [r_in, r_out], "min_snr": min_snr,
                                      "n_rows_in_cone": n_rows, "sources": srcs}
            qs["allwise_density"] = "ok"
    except Exception as exc:
        qs["allwise_density"] = _failed(exc)

    # 8. Youth indicators: GALEX NUV, ROSAT 2RXS, eRASS1 and SIMBAD group membership
    youth: dict[str, Any] = {}
    g_ra, g_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2008.0)  # mid GALEX AIS epoch
    try:
        g_r = float(vcfg["youth_galex_match_radius_arcsec"])
        g_rows, m_gal, _ = _vizier_tsv_cone_cached(source="II/335/galex_ais", ra=g_ra, dec=g_dec, radius_arcsec=g_r, archive="high_energy", release="galex_guvcat_ais", cache=cache)
        manifests.append(m_gal)
        cands = []
        for gr in g_rows:
            gra, gde = _vizier_radec(gr)
            nuv = _vizier_value(gr, "NUVmag", "NUV", "nuv_mag")
            if gra is None or gde is None or nuv is None:
                continue
            cands.append({"sep_arcsec": round(angular_separation_arcsec(g_ra, g_dec, gra, gde), 2), "nuv_mag": nuv,
                          "nuv_err": _vizier_value(gr, "e_NUVmag", "e_NUV", "nuv_magerr"),
                          "fuv_mag": _vizier_value(gr, "FUVmag", "FUV", "fuv_mag")})
        youth["galex"] = {"status": "ok", "catalogue": "VizieR II/335/galex_ais (GUVcat AIS)", "match_radius_arcsec": g_r,
                          "match": min(cands, key=lambda c: c["sep_arcsec"]) if cands else None}
        qs["galex"] = "ok"
    except Exception as exc:
        qs["galex"] = _failed(exc)
        youth["galex"] = {"status": qs["galex"]}

    r_ra, r_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 1990.8)  # ROSAT all-sky survey epoch
    try:
        x_r = float(vcfg["youth_rosat_match_radius_arcsec"])
        x_rows, m_x, _ = _vizier_tsv_cone_cached(source="IX/47/2rxs", ra=r_ra, dec=r_dec, radius_arcsec=x_r, archive="high_energy", release="rosat_2rxs", cache=cache)
        manifests.append(m_x)
        cands = []
        for xr in x_rows:
            xra, xde = _vizier_radec(xr)
            rate = _vizier_value(xr, "Rate", "CRate", "CR", "RATE", "Cts")
            if xra is None or xde is None or rate is None:
                continue
            cands.append({"sep_arcsec": round(angular_separation_arcsec(r_ra, r_dec, xra, xde), 1), "rate_cts_s": rate,
                          "hr1": _vizier_value(xr, "HR1", "HR_1")})
        youth["rosat_2rxs"] = {"status": "ok", "catalogue": "VizieR IX/47/2rxs", "match_radius_arcsec": x_r,
                               "match": min(cands, key=lambda c: c["sep_arcsec"]) if cands else None}
        qs["xray_2rxs"] = "ok"
    except Exception as exc:
        qs["xray_2rxs"] = _failed(exc)
        youth["rosat_2rxs"] = {"status": qs["xray_2rxs"]}

    gal_l, _ = equatorial_to_galactic(ra, dec)
    if not (180.0 <= gal_l < 360.0):
        qs["xray_erass1"] = "no_coverage"
        youth["erass1"] = {"status": "no_coverage", "note": f"eRASS1 DR covers Galactic l 180-360 deg; l={gal_l:.1f}"}
    else:
        e_ra, e_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, 2020.1)
        try:
            e_r = float(vcfg["youth_erass1_match_radius_arcsec"])
            e_rows, m_e, _ = _vizier_tsv_cone_cached(source="J/A+A/682/A34/erass1-m", ra=e_ra, dec=e_dec, radius_arcsec=e_r, archive="high_energy", release="erass1_main", cache=cache)
            manifests.append(m_e)
            cands = []
            for er in e_rows:
                era, ede = _vizier_radec(er)
                fl = _vizier_value(er, "Flux", "FluxP1", "ML_FLUX_1", "F0.2-2.3", "Fx")
                if era is None or ede is None or fl is None:
                    continue
                cands.append({"sep_arcsec": round(angular_separation_arcsec(e_ra, e_dec, era, ede), 1), "flux_erg_cm2_s": fl})
            youth["erass1"] = {"status": "ok", "catalogue": "VizieR J/A+A/682/A34/erass1-m", "match_radius_arcsec": e_r,
                               "match": min(cands, key=lambda c: c["sep_arcsec"]) if cands else None}
            qs["xray_erass1"] = "ok"
        except Exception as exc:
            qs["xray_erass1"] = _failed(exc)
            youth["erass1"] = {"status": qs["xray_erass1"]}

    try:
        q_mem = (
            f"SELECT p.main_id, p.otype, h.membership FROM h_link AS h "
            f"JOIN basic AS p ON h.parent = p.oid JOIN ident AS i ON i.oidref = h.child "
            f"WHERE i.id = 'Gaia DR3 {sid}'"
        )
        mem_rows, m_mem, _ = _tap_query_cached(urls=SIMBAD_TAP_URLS, query=q_mem, archive="simbad", release="live_membership", cache=cache, timeout_s=15)
        manifests.append(m_mem)
        parents = [{"main_id": str(r.get("main_id")), "otype": str(r.get("otype") or ""), "membership": r.get("membership")} for r in mem_rows]
        youth["association"] = {
            "status": "ok",
            "source": "SIMBAD h_link parents",
            "parents": parents,
            "young_groups": [p for p in parents if p["otype"] in YOUNG_GROUP_OTYPES],
            "clusters": [p for p in parents if p["otype"] in CLUSTER_OTYPES],
            "simbad_object_type": out.get("simbad_type") if qs.get("simbad") == "ok" else None,
            "young_object_type": (out.get("simbad_type") in YOUNG_STAR_OTYPES) if qs.get("simbad") == "ok" else None,
        }
        qs["young_association"] = "ok"
    except Exception as exc:
        qs["young_association"] = _failed(exc)
        youth["association"] = {"status": qs["young_association"]}
    va["youth"] = youth

    # 9. Known Solar System bodies at the 2010 WISE single-exposure epochs (IMCCE SkyBoT, observer C51)
    if not l1b_ok:
        qs["skybot"] = "failed: no frame epochs because wise_l1b failed"
    elif not frame_mjds:
        qs["skybot"] = "no_coverage"
        va["sso_search"] = {"frame_mjds": [], "visits": [], "matches": [], "observer": "C51"}
    else:
        sso, m_sso = _search_sso(ra=ra, dec=dec, pmra=pmra, pmdec=pmdec, frame_mjds=frame_mjds, vcfg=vcfg, cache=cache, sleep=sleep)
        manifests.extend(m_sso)
        va["sso_search"] = sso
        bad = [v["status"] for v in sso["visits"] if v["status"] != "ok"]
        qs["skybot"] = "ok" if not bad else "failed: " + " | ".join(bad)

    # 10. Line-of-sight Galactic dust (IRSA DUST service). Stored apart from the per-star `ebv_sfd` used for
    #     stellar extinction, which is left untouched here.
    try:
        dust, m_dust = fetch_irsa_dust(ra_deg=ra, dec_deg=dec, cache=cache, sleep=sleep)
        manifests.append(m_dust)
        va["dust"] = dust
        qs["irsa_dust"] = "ok"
    except Exception as exc:
        qs["irsa_dust"] = _failed(exc)

    # 11. Real-pixel cutouts: AllWISE Atlas W1-W4 (measured) and a PS1 r stack (display only, Dec > -30)
    if frame_mjds:
        wise_epoch = 2000.0 + (float(np.mean(frame_mjds)) - 51544.5) / 365.25
        epoch_basis = "mean 2010 L1b frame epoch"
    else:
        wise_epoch = 2010.4
        epoch_basis = "nominal 4-band cryogenic mid-epoch"
    p_ra, p_dec = propagate_position(ra, dec, pmra, pmdec, 2016.0, wise_epoch)
    size = float(vcfg["vision_cutout_size_arcsec"])
    cut: dict[str, Any] = {"epoch_jyear": wise_epoch, "epoch_basis": epoch_basis, "pred_ra_deg": p_ra, "pred_dec_deg": p_dec,
                           "size_arcsec": size, "bands": {}, "log": []}
    try:
        images, m_cut, log = fetch_allwise_cutouts(ra=p_ra, dec=p_dec, size_arcsec=size, cache=cache, sleep=sleep)
        manifests.extend(m_cut)
        cut["log"] = log
        cut["bands"] = measure_wise_cutouts(images, ra_pred=p_ra, dec_pred=p_dec, vetting=vcfg)
        missing = [b for b in ("W3", "W4") if b not in images]
        qs["cutouts_wise"] = "ok" if not missing else f"failed: no Atlas cutout for {missing} ({'; '.join(log)})"
    except Exception as exc:
        qs["cutouts_wise"] = _failed(exc)
    if dec > -30.0:
        try:
            ps1_cut, m_pc = fetch_ps1_stack_cutout(ra=ra, dec=dec, size_arcsec=min(60.0, size), cache=cache, sleep=sleep)
            manifests.extend(m_pc)
            cut["ps1"] = ps1_cut
            qs["cutouts_ps1"] = "ok"
        except Exception as exc:
            qs["cutouts_ps1"] = _failed(exc)
    else:
        qs["cutouts_ps1"] = "no_coverage"
    va["cutouts"] = cut

    # 12. Optical time series for light-curve vetting only (pipeline_agreement, second_survey).
    #     Stored under lightcurves.tess_sectors / lightcurves.optical_surveys so channel inputs stay unchanged.
    if enrich_optical:
        opt_mans = _enrich_optical_timeseries(out, cache, cfg)
        manifests.extend(opt_mans)

    out["live_enriched"] = True
    return out, manifests
