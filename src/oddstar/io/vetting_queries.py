"""Single-purpose archive queries used only by infrared vetting checks (line-of-sight dust, HST/JWST coverage).

- IRSA Galactic Dust Reddening and Extinction service
  `https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=<ra>+<dec>+equ+j2000` (XML). We read the
  reference-pixel (line-of-sight) SFD (1998) and Schlafly & Finkbeiner (2011) E(B-V) and the IRAS/SFD 100 um
  intensity; the region means (over the service's default `regSize` box) are kept only as context.
- MAST CAOM position search `Mast.Caom.Filtered.Position` through `https://mast.stsci.edu/api/v0/invoke`,
  restricted to HST and JWST image products, to learn whether sub-arcsecond imaging of the position exists.
"""

from __future__ import annotations

import json
import re
import time
import xml.etree.ElementTree as ET
from typing import Any

from oddstar.errors import ArchiveAccessError
from oddstar.io.archive_http import Sleeper, cached_get
from oddstar.io.cache import RawCache
from oddstar.io.physics import BAND_WAVELENGTH_UM, stellar_photosphere_flux_mjy
from oddstar.schemas import ManifestRecord


def photosphere_flux_mjy(star: dict[str, Any], band: str) -> float | None:
    """Photospheric flux in `band` from a Planck photosphere at `teff_k` normalised to the observed Ks flux.

    This is the vetting-side estimate used to size the infrared excess (for chance-alignment and neighbour-screen
    flux cuts); it is not the channel's SED fit. Returns None when Teff or Ks is unavailable.
    """
    teff = star.get("teff_k")
    f_ks = (star.get("flux_mjy") or {}).get("Ks")
    if teff is None or f_ks is None or float(f_ks) <= 0.0:
        return None
    return float(stellar_photosphere_flux_mjy(BAND_WAVELENGTH_UM[band], float(teff), float(f_ks)))


def target_excess_mjy(star: dict[str, Any], band: str) -> float | None:
    """Observed minus photospheric flux (mJy) in `band`, or None when either is unavailable."""
    obs = (star.get("flux_mjy") or {}).get(band)
    phot = photosphere_flux_mjy(star, band)
    if obs is None or phot is None:
        return None
    return float(obs) - phot

IRSA_DUST_URL = "https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust"
MAST_INVOKE_URL = "https://mast.stsci.edu/api/v0/invoke"

_FLOAT_RE = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?")


def _first_float(text: str | None) -> float | None:
    if text is None:
        return None
    m = _FLOAT_RE.search(text)
    return float(m.group(0)) if m else None


def parse_irsa_dust_xml(text: str) -> dict[str, float | None]:
    """Parse the nph-dust XML. Raises ArchiveAccessError on an error status or a missing SFD value."""
    try:
        root = ET.fromstring(text.strip())
    except ET.ParseError as exc:
        raise ArchiveAccessError(f"IRSA dust XML parse failed: {exc}", retryable=True) from exc
    status = (root.attrib.get("status") or "").lower()
    if status and status != "ok":
        msg = root.findtext(".//message") or text[:200]
        raise ArchiveAccessError(f"IRSA dust service status={status}: {msg.strip()}", retryable=False)
    out: dict[str, float | None] = {
        "ebv_sfd_los": None,
        "ebv_sf11_los": None,
        "i100_mjy_sr": None,
        "ebv_sfd_region_mean": None,
        "ebv_sf11_region_mean": None,
        "i100_region_mean_mjy_sr": None,
    }
    for res in root.iter("result"):
        desc = (res.findtext("desc") or "").lower()
        stats = res.find("statistics")
        if stats is None:
            continue
        if "e(b-v)" in desc or "reddening" in desc:
            out["ebv_sfd_los"] = _first_float(stats.findtext("refPixelValueSFD"))
            out["ebv_sf11_los"] = _first_float(stats.findtext("refPixelValueSandF"))
            out["ebv_sfd_region_mean"] = _first_float(stats.findtext("meanValueSFD"))
            out["ebv_sf11_region_mean"] = _first_float(stats.findtext("meanValueSandF"))
        elif "100 micron" in desc or "100um" in desc:
            out["i100_mjy_sr"] = _first_float(stats.findtext("refPixelValue"))
            out["i100_region_mean_mjy_sr"] = _first_float(stats.findtext("meanValue"))
    if out["ebv_sfd_los"] is None:
        raise ArchiveAccessError("IRSA dust response has no reference-pixel SFD E(B-V)", retryable=False)
    return out


def fetch_irsa_dust(
    *, ra_deg: float, dec_deg: float, cache: RawCache, timeout_s: float = 30.0, sleep: Sleeper = time.sleep
) -> tuple[dict[str, float | None], ManifestRecord]:
    """Line-of-sight dust values at (ra, dec) from the IRSA DUST service."""
    params = {"locstr": f"{ra_deg:.6f} {dec_deg:+.6f} equ j2000"}
    payload, man, _ = cached_get(
        url=IRSA_DUST_URL,
        params=params,
        archive="dustmaps",
        release="irsa_dust",
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
        validate=lambda b: None if parse_irsa_dust_xml(b.decode("utf-8", errors="replace")) else None,
    )
    return parse_irsa_dust_xml(payload.decode("utf-8", errors="replace")), man


CAOM_COLUMNS = "obs_collection,instrument_name,filters,obs_id,target_name,t_exptime,t_min,s_ra,s_dec,dataproduct_type,calib_level"


def fetch_subarcsec_imaging(
    *,
    ra_deg: float,
    dec_deg: float,
    radius_arcsec: float,
    cache: RawCache,
    timeout_s: float = 60.0,
    sleep: Sleeper = time.sleep,
) -> tuple[list[dict[str, Any]], ManifestRecord]:
    """HST/JWST image observations whose footprint covers a cone around (ra, dec), from MAST CAOM."""
    req = {
        "service": "Mast.Caom.Filtered.Position",
        "format": "json",
        "params": {
            "columns": CAOM_COLUMNS,
            "filters": [
                {"paramName": "obs_collection", "values": ["HST", "JWST"]},
                {"paramName": "dataproduct_type", "values": ["image"]},
            ],
            "position": f"{ra_deg:.6f}, {dec_deg:.6f}, {radius_arcsec / 3600.0:.6f}",
        },
    }

    def _parse(payload: bytes) -> list[dict[str, Any]]:
        try:
            js = json.loads(payload.decode("utf-8") or "{}")
        except (json.JSONDecodeError, UnicodeDecodeError) as exc:
            raise ArchiveAccessError(f"MAST CAOM JSON parse failed: {exc}", retryable=True) from exc
        if not isinstance(js, dict) or "data" not in js:
            raise ArchiveAccessError(f"MAST CAOM response has no 'data': {str(js)[:200]}", retryable=True)
        # Only a COMPLETE answer is a result: EXECUTING carries no (or partial) rows and would read as "no observations"
        if str(js.get("status", "")).upper() != "COMPLETE":
            raise ArchiveAccessError(f"MAST CAOM status={js.get('status')!r}: {js.get('msg', '')}", retryable=True)
        return [dict(r) for r in js["data"]]

    payload, man, _ = cached_get(
        url=MAST_INVOKE_URL,
        params={"request": json.dumps(req, sort_keys=True)},
        archive="deep_optical",
        release="mast_caom_hst_jwst",
        cache=cache,
        timeout_s=timeout_s,
        sleep=sleep,
        validate=lambda b: None if _parse(b) is not None else None,
    )
    rows = _parse(payload)
    obs = [
        {
            "collection": str(r.get("obs_collection") or ""),
            "instrument": str(r.get("instrument_name") or ""),
            "filters": str(r.get("filters") or ""),
            "obs_id": str(r.get("obs_id") or ""),
            "exptime_s": r.get("t_exptime"),
            "t_min_mjd": r.get("t_min"),
        }
        for r in rows
    ]
    return obs, man
