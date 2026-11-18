"""Ultraviolet and X-ray catalogue lookups for the youth check, from the archives that own the catalogues.

These replaced VizieR mirror queries after the 2026-11-17 live test: the CfA VizieR mirror answered HTTP 200 with
table metadata and zero rows for GALEX and eRASS1 everywhere on the sky (it does not serve those tables), which the
pipeline recorded as "queried, no counterpart", and the 2RXS catalogue identifier in use did not exist. An empty
answer is only evidence when the service is known to hold the table, so each lookup here goes to the owning archive
and raises on anything other than a well-formed answer.

- GALEX: MAST `Mast.Galex.Catalog` (GET https://mast.stsci.edu/api/v0/invoke?request=<json>), radius in degrees,
  JSON rows with `ra`, `dec`, `nuv_mag`, `nuv_magerr`, `fuv_mag`, `fuv_magerr`, `survey`, `nuv_exptime`; missing
  magnitudes are -999. Verified live 2026-11-17.
- ROSAT 2RXS: HEASARC TAP table `rass2rxs` (`name`, `ra`, `dec`, `count_rate` [ct/s], `count_rate_error`,
  `hardness_ratio_1`, `exposure`). Verified live 2026-11-17.
- eROSITA eRASS1 main catalogue: HEASARC TAP table `erass1main` (`name`, `ra`, `dec`, `b1_flux` and `b1_flux_error`
  [erg cm^-2 s^-1, band 1 = 0.2-2.3 keV]). Table and column names verified live 2026-11-17.

HEASARC's TAP service returns VOTable only (FORMAT=csv is refused); a query error is an `INFO name="QUERY_STATUS"
value="ERROR"` element and no table.
"""

from __future__ import annotations

import io
import json
import math
import re
import time
from typing import Any

import numpy as np

from oddstar.errors import ArchiveAccessError
from oddstar.io.archive_http import Sleeper, cached_get
from oddstar.io.cache import RawCache
from oddstar.io.physics import angular_separation_arcsec
from oddstar.schemas import ManifestRecord

MAST_INVOKE_URL = "https://mast.stsci.edu/api/v0/invoke"
HEASARC_TAP_URL = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
GALEX_PAGE_SIZE = 5000
GALEX_MISSING_BELOW = -90.0  # MAST marks a missing magnitude as -999 (or -99)
_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _mag(val: Any) -> float | None:
    try:
        v = float(val)
    except (TypeError, ValueError):
        return None
    return None if (math.isnan(v) or v <= GALEX_MISSING_BELOW) else v


def parse_galex_json(payload: bytes) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Rows and paging block of a `Mast.Galex.Catalog` answer; raises unless the answer is complete and well formed."""
    try:
        js = json.loads(payload.decode("utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError) as exc:
        raise ArchiveAccessError(f"MAST GALEX JSON parse failed: {exc}", retryable=True) from exc
    if not isinstance(js, dict) or not isinstance(js.get("data"), list):
        raise ArchiveAccessError(f"MAST GALEX response has no 'data' list: {str(js)[:200]}", retryable=True)
    if str(js.get("status", "")).upper() != "COMPLETE":
        raise ArchiveAccessError(f"MAST GALEX status={js.get('status')!r}: {js.get('msg', '')}", retryable=True)
    return [dict(r) for r in js["data"]], dict(js.get("paging") or {})


def fetch_galex_sources(
    *,
    ra_deg: float,
    dec_deg: float,
    radius_arcsec: float,
    cache: RawCache,
    timeout_s: float = 60.0,
    sleep: Sleeper = time.sleep,
    allow_truncated: bool = False,
) -> tuple[list[dict[str, Any]], ManifestRecord]:
    """GALEX catalogue sources (all surveys) within `radius_arcsec`, nearest first.

    Each row: sep_arcsec, ra_deg, dec_deg, nuv_mag, nuv_err, fuv_mag, fuv_err, survey, nuv_exptime_s. Magnitudes are
    None where GALEX has no measurement in that band. Raises ArchiveAccessError when MAST reports more rows than one
    page holds (the cone would be silently truncated), unless `allow_truncated` (for a caller that only asks whether
    the field holds any source at all and can work from the first page).
    """
    req = {
        "service": "Mast.Galex.Catalog",
        "params": {"ra": round(float(ra_deg), 6), "dec": round(float(dec_deg), 6), "radius": round(radius_arcsec / 3600.0, 6)},
        "format": "json",
        "pagesize": GALEX_PAGE_SIZE,
        "page": 1,
    }
    payload, man, _ = cached_get(
        url=MAST_INVOKE_URL,
        params={"request": json.dumps(req, sort_keys=True)},
        archive="high_energy",
        release="galex_mast",
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
        validate=lambda b: None if parse_galex_json(b) is not None else None,
    )
    rows, paging = parse_galex_json(payload)
    total = paging.get("rowsTotal")
    if total is not None and int(total) > len(rows) and not allow_truncated:
        raise ArchiveAccessError(f"MAST GALEX cone truncated: {len(rows)} of {total} rows returned", retryable=False)
    out: list[dict[str, Any]] = []
    for r in rows:
        try:
            s_ra, s_dec = float(r["ra"]), float(r["dec"])
        except (KeyError, TypeError, ValueError):
            continue
        out.append(
            {
                "sep_arcsec": round(angular_separation_arcsec(ra_deg, dec_deg, s_ra, s_dec), 2),
                "ra_deg": s_ra,
                "dec_deg": s_dec,
                "nuv_mag": _mag(r.get("nuv_mag")),
                "nuv_err": _mag(r.get("nuv_magerr")),
                "fuv_mag": _mag(r.get("fuv_mag")),
                "fuv_err": _mag(r.get("fuv_magerr")),
                "survey": str(r.get("survey") or ""),
                "nuv_exptime_s": _mag(r.get("nuv_exptime")),
            }
        )
    out.sort(key=lambda x: x["sep_arcsec"])
    return out, man


def parse_heasarc_votable(payload: bytes, columns: list[str]) -> list[dict[str, Any]]:
    """Rows of a HEASARC TAP VOTable answer; raises on a QUERY_STATUS error or a missing table or column."""
    text_head = payload[:4000].decode("utf-8", errors="replace")
    m = re.search(r'name="QUERY_STATUS"\s+value="([A-Za-z]+)"[^>]*>([^<]*)<', text_head)
    if m is None:
        raise ArchiveAccessError(f"HEASARC TAP answer has no QUERY_STATUS: {text_head[:200]!r}", retryable=True)
    if m.group(1).upper() != "OK":
        raise ArchiveAccessError(f"HEASARC TAP query error: {m.group(2).strip()[:300]}", retryable=False)
    from astropy.io.votable import parse

    try:
        table = parse(io.BytesIO(payload)).get_first_table().to_table(use_names_over_ids=True)
    except Exception as exc:  # astropy raises several types on a malformed body
        raise ArchiveAccessError(f"HEASARC TAP VOTable not readable: {type(exc).__name__}: {exc}", retryable=True) from exc
    missing = [c for c in columns if c not in table.colnames]
    if missing:
        raise ArchiveAccessError(f"HEASARC TAP answer lacks column(s) {missing}; has {table.colnames}", retryable=False)
    rows: list[dict[str, Any]] = []
    for rec in table:
        row: dict[str, Any] = {}
        for c in columns:
            v = rec[c]
            if v is np.ma.masked:
                row[c] = None
                continue
            if isinstance(v, bytes):
                v = v.decode("utf-8", errors="replace")
            if isinstance(v, str):
                row[c] = v.strip()
                continue
            try:
                f = float(v)
            except (TypeError, ValueError):
                row[c] = None
                continue
            row[c] = None if math.isnan(f) else f
        rows.append(row)
    return rows


def heasarc_cone(
    *,
    table: str,
    columns: list[str],
    ra_deg: float,
    dec_deg: float,
    radius_arcsec: float,
    cache: RawCache,
    release: str,
    timeout_s: float = 60.0,
    sleep: Sleeper = time.sleep,
) -> tuple[list[dict[str, Any]], ManifestRecord]:
    """Cone search of one HEASARC table through its TAP service; rows carry `sep_arcsec`, nearest first.

    `columns` must include `ra` and `dec`. A table or column HEASARC does not have is a query error (raised), so an
    empty list always means "the table was searched and has no row in the cone".
    """
    for name in [table, *columns]:
        if not _IDENT.match(name):
            raise ArchiveAccessError(f"Unsafe HEASARC identifier {name!r}", retryable=False)
    if "ra" not in columns or "dec" not in columns:
        raise ArchiveAccessError("heasarc_cone needs the ra and dec columns", retryable=False)
    adql = (
        f"SELECT {', '.join(columns)} FROM {table} "
        f"WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra_deg:.6f}, {dec_deg:.6f}, {radius_arcsec / 3600.0:.6f}))"
    )
    payload, man, _ = cached_get(
        url=HEASARC_TAP_URL,
        params={"REQUEST": "doQuery", "LANG": "ADQL", "QUERY": adql},
        archive="high_energy",
        release=release,
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
        validate=lambda b: None if parse_heasarc_votable(b, columns) is not None else None,
    )
    rows = parse_heasarc_votable(payload, columns)
    out = []
    for r in rows:
        if r.get("ra") is None or r.get("dec") is None:
            continue
        out.append({**r, "sep_arcsec": round(angular_separation_arcsec(ra_deg, dec_deg, float(r["ra"]), float(r["dec"])), 1)})
    out.sort(key=lambda x: x["sep_arcsec"])
    return out, man
