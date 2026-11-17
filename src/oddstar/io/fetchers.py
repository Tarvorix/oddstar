"""Stage 1 public archive fetchers with rate limiting, retries, release adapters, and cache."""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.errors import ArchiveAccessError, PolicyViolationError
from oddstar.io.cache import RawCache
from oddstar.io.fixtures_data import RADIO_CADENCES, generate_field_sample, get_named_object
from oddstar.schemas import ManifestRecord


@dataclass(frozen=True)
class ArchiveSpec:
    name: str
    description: str
    access_method: str
    table_or_service: str


ARCHIVE_SPECS: dict[str, ArchiveSpec] = {
    "gaia": ArchiveSpec(
        "gaia",
        "Positions, parallaxes, proper motions, G/BP/RP photometry, XP spectra, RUWE, NSS",
        "ESA Gaia Archive TAP via astroquery.gaia",
        "gaiadr3.gaia_source / gaiadr4.gaia_source",
    ),
    "2mass": ArchiveSpec("2mass", "J, H, Ks photometry", "IRSA TAP", "fp_psc"),
    "allwise": ArchiveSpec("allwise", "W1-W4 (3.4, 4.6, 12, 22 um) photometry & flags", "IRSA TAP", "allwise_p3as_psd"),
    "wise_single_exposure": ArchiveSpec(
        "wise_single_exposure",
        "2010 cryogenic single-frame photometry & cutouts",
        "IRSA TAP + SIA",
        "allsky_4band_p1bs_psd",
    ),
    "neowise": ArchiveSpec(
        "neowise",
        "NEOWISE-R W1/W2 2013-2024 time series & unWISE/unTimely coadds",
        "IRSA TAP + unWISE",
        "neowiser_p1bs_psd",
    ),
    "tess": ArchiveSpec("tess", "TESS full-frame light curves & TESScut cutouts", "MAST astroquery.mast / TESScut", "TESS"),
    "kepler": ArchiveSpec("kepler", "Kepler and K2 high-precision light curves", "MAST astroquery.mast", "Kepler/K2"),
    "ztf": ArchiveSpec("ztf", "Northern-sky g and r light curves since 2018", "IRSA ZTF light-curve API", "ztf_objects"),
    "asassn": ArchiveSpec("asassn", "All-sky decade-baseline optical light curves", "ASAS-SN Sky Patrol", "skypatrol_v2"),
    "simbad": ArchiveSpec("simbad", "Known object types & bibliography", "astroquery.simbad", "SIMBAD"),
    "vizier": ArchiveSpec("vizier", "Thousands of published astronomical catalogs", "astroquery.vizier", "VizieR"),
    "radio_continuum": ArchiveSpec(
        "radio_continuum",
        "VLASS, NVSS, FIRST, RACS, LoTSS, TGSS, SUMSS radio catalogs",
        "VizieR + CIRADA cutouts",
        "radio_surveys",
    ),
    "deep_optical": ArchiveSpec(
        "deep_optical",
        "Pan-STARRS, Legacy Surveys / DECam, SDSS imaging & catalogs",
        "MAST + NOIRLab DataLab + SDSS",
        "ps1_ls_sdss",
    ),
    "high_energy": ArchiveSpec(
        "high_energy",
        "GALEX UV and ROSAT, eROSITA DR1, XMM, Chandra X-ray catalogs",
        "MAST + VizieR + eROSITA TAP",
        "galex_xray",
    ),
    "mid_far_ir": ArchiveSpec(
        "mid_far_ir",
        "Spitzer, AKARI, IRAS, Herschel point-source catalogs",
        "IRSA + VizieR",
        "spitzer_akari_iras_herschel",
    ),
    "spherex": ArchiveSpec(
        "spherex",
        "0.75-5 um weekly quick-release spectrophotometric images",
        "IRSA SPHEREx TAP/SIA",
        "spherex_qr",
    ),
    "dustmaps": ArchiveSpec("dustmaps", "SFD, Planck, and 3D Bayestar extinction maps", "dustmaps / IRSA", "sfd_planck"),
    "ephemeris": ArchiveSpec(
        "ephemeris",
        "Solar-system moving-object ephemerides at WISE/plate epochs",
        "IRSA MOST + IMCCE SkyBoT",
        "skybot_most",
    ),
    "breakthrough_listen": ArchiveSpec(
        "breakthrough_listen",
        "Green Bank & Parkes filterbank/HDF5 & APF optical spectra",
        "Breakthrough Listen Open Data Archive",
        "bl_opendata",
    ),
    "survey_spectra": ArchiveSpec(
        "survey_spectra",
        "SDSS, LAMOST, GALAH, APOGEE, Gaia XP/RVS, APF spectra",
        "SDSS/LAMOST/GALAH/Gaia TAP",
        "stellar_spectra",
    ),
    "usnob1": ArchiveSpec("usnob1", "20th-century photographic plate catalog", "VizieR I/284", "I/284/out"),
    "historical_plates": ArchiveSpec("historical_plates", "Harvard DASCH and Digitized Sky Survey plates", "DASCH + DSS SIA", "dasch_dss"),
    "hsc_ssp": ArchiveSpec(
        "hsc_ssp",
        "Deep Subaru HSC-SSP imaging (OFF by default; requires registration)",
        "NAOJ HSC-SSP TAP",
        "pdr3_wide",
    ),
}


class ArchiveFetcher:
    """Rate-limited, retrying fetcher writing through content-addressed RawCache."""

    def __init__(self, archive: str, cfg: PipelineConfig, cache: RawCache) -> None:
        if archive not in ARCHIVE_SPECS:
            raise ArchiveAccessError(
                f"Unknown archive: {archive}",
                hint=f"Registered archives: {', '.join(sorted(ARCHIVE_SPECS))}",
                retryable=False,
            )
        self.archive = archive
        self.spec = ARCHIVE_SPECS[archive]
        self.cfg = cfg
        self.cache = cache
        self.policy = cfg.archives[archive]
        self._last_request_ts = 0.0
        self.backend_calls = 0

    def release(self) -> str:
        if self.archive == "gaia":
            return self.cfg.gaia.release
        return self.policy.release

    def _respect_rate_limit(self) -> None:
        min_interval = 1.0 / max(0.1, self.policy.rate_limit_rps)
        elapsed = time.monotonic() - self._last_request_ts
        if elapsed < min_interval and self._last_request_ts > 0.0:
            time.sleep(min(0.05, min_interval - elapsed))
        self._last_request_ts = time.monotonic()

    @staticmethod
    def build_signature(mode: str, query: dict[str, Any]) -> str:
        canonical = json.dumps({"mode": mode, "query": query}, sort_keys=True, separators=(",", ":"))
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()[:24]

    def fetch(
        self,
        *,
        ids: list[str] | None = None,
        cone: tuple[float, float, float] | None = None,
        adql: str | None = None,
        live: bool = False,
        max_retries: int = 3,
    ) -> tuple[list[dict[str, Any]], ManifestRecord, bool]:
        """Fetch records from cache or archive backend.

        Returns (rows, manifest_record, cache_hit).
        When cache_hit is True, zero network/backend calls are performed.
        """
        if not self.policy.enabled:
            raise PolicyViolationError(
                f"Archive '{self.archive}' is disabled in configuration.",
                hint="Registered-access archives (like hsc_ssp) require explicit human acceptance of terms.",
            )

        query_spec: dict[str, Any] = {}
        if ids is not None:
            query_spec["ids"] = [str(i).strip() for i in ids]
        if cone is not None:
            query_spec["cone"] = {"ra": cone[0], "dec": cone[1], "radius_arcsec": cone[2]}
        if adql is not None:
            query_spec["adql"] = adql

        sig = self.build_signature("query", query_spec)
        rel = self.release()
        cached = self.cache.lookup(self.archive, rel, sig)
        if cached is not None:
            rows = json.loads(self.cache.read_bytes(cached).decode("utf-8"))
            return rows, cached, True

        last_err: Exception | None = None
        for attempt in range(max_retries):
            try:
                self._respect_rate_limit()
                self.backend_calls += 1
                if live:
                    rows = self._query_live_archive(ids=ids, cone=cone, adql=adql, release=rel)
                else:
                    rows = self._query_offline_archive(ids=ids, cone=cone, adql=adql, release=rel)
                payload = json.dumps(rows, sort_keys=True).encode("utf-8")
                query_text = adql or json.dumps(query_spec, sort_keys=True)
                manifest = self.cache.put(
                    archive=self.archive,
                    release=rel,
                    request_signature=sig,
                    query_text=query_text,
                    payload=payload,
                    row_count=len(rows),
                    offline_fixture=not live,
                )
                return rows, manifest, False
            except Exception as exc:
                last_err = exc
                time.sleep(0.05 * (2**attempt))

        raise ArchiveAccessError(
            f"Archive query failed for {self.archive} after {max_retries} attempts: {last_err}",
            hint="Verify network connectivity or rerun offline without --live",
        )

    def _query_offline_archive(
        self,
        *,
        ids: list[str] | None,
        cone: tuple[float, float, float] | None,
        adql: str | None,
        release: str,
    ) -> list[dict[str, Any]]:
        gaia_rel = self.cfg.gaia.release
        if self.archive == "breakthrough_listen" and ids:
            rows = []
            for tid in ids:
                if tid in RADIO_CADENCES:
                    rows.append(RADIO_CADENCES[tid])
            if rows:
                return rows

        if ids:
            field_bank = {
                s["source_id"]: s
                for s in generate_field_sample(
                    n_stars=max(200, len(ids)),
                    seed=self.cfg.sample.seed,
                    gaia_release=gaia_rel,
                )
            }
            out = []
            for raw_id in ids:
                clean = raw_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
                star = get_named_object(clean, gaia_release=gaia_rel) or field_bank.get(clean)
                if star is None:
                    # Synthesize deterministic star indexed to this ID
                    from oddstar.io.fixtures_data import _make_base_star

                    star = _make_base_star(clean, gaia_release=gaia_rel)
                out.append(self._project_for_archive(star))
            return out

        if cone is not None:
            ra0, dec0, rad_arcsec = cone
            from oddstar.io.physics import angular_separation_arcsec

            candidates = list(
                generate_field_sample(
                    n_stars=500, seed=self.cfg.sample.seed, gaia_release=gaia_rel
                )
            )
            results = []
            for star in candidates:
                sep = angular_separation_arcsec(ra0, dec0, star["ra_deg"], star["dec_deg"])
                if sep <= rad_arcsec:
                    proj = self._project_for_archive(star)
                    proj["separation_arcsec"] = sep
                    results.append(proj)
            return results

        # ADQL bulk sample query
        n_rows = min(200, self.cfg.sample.field_sample_size)
        if adql and "LIMIT" in adql.upper():
            try:
                tok = adql.upper().split("LIMIT")[1].strip().split()[0]
                n_rows = min(10000, max(1, int(tok)))
            except ValueError:
                pass
        return [
            self._project_for_archive(s)
            for s in generate_field_sample(n_stars=n_rows, seed=self.cfg.sample.seed, gaia_release=gaia_rel)
        ]

    def _project_for_archive(self, star: dict[str, Any]) -> dict[str, Any]:
        """Return archive-specific view while keeping source_id & coordinates."""
        base = {
            "source_id": star["source_id"],
            "archive": self.archive,
            "release": self.release(),
            "ra_deg": star["ra_deg"],
            "dec_deg": star["dec_deg"],
        }
        if self.archive == "gaia":
            base.update(
                {
                    "gaia_release": star["gaia_release"],
                    "parallax_mas": star["parallax_mas"],
                    "parallax_err_mas": star["parallax_err_mas"],
                    "parallax_snr": star["parallax_snr"],
                    "distance_pc": star["distance_pc"],
                    "pmra_mas_yr": star["pmra_mas_yr"],
                    "pmdec_mas_yr": star["pmdec_mas_yr"],
                    "ruwe": star["ruwe"],
                    "teff_k": star["teff_k"],
                    "logg": star["logg"],
                    "feh": star["feh"],
                    "luminosity_lsun": star["luminosity_lsun"],
                    "mags": {k: star["mags"][k] for k in ("G", "BP", "RP")},
                    "flux_mjy": {k: star["flux_mjy"][k] for k in ("G", "BP", "RP")},
                    "flux_err_mjy": {k: star["flux_err_mjy"][k] for k in ("G", "BP", "RP")},
                }
            )
        elif self.archive == "2mass":
            base.update(
                {
                    "clean_2mass": star.get("2mass_clean", True),
                    "mags": {k: star["mags"][k] for k in ("J", "H", "Ks")},
                    "flux_mjy": {k: star["flux_mjy"][k] for k in ("J", "H", "Ks")},
                    "flux_err_mjy": {k: star["flux_err_mjy"][k] for k in ("J", "H", "Ks")},
                }
            )
        elif self.archive == "allwise":
            base.update(star["allwise"])
            base["mags"] = {k: star["mags"][k] for k in ("W1", "W2", "W3", "W4")}
            base["flux_mjy"] = {k: star["flux_mjy"][k] for k in ("W1", "W2", "W3", "W4")}
            base["flux_err_mjy"] = {k: star["flux_err_mjy"][k] for k in ("W1", "W2", "W3", "W4")}
        elif self.archive == "wise_single_exposure":
            base["w4_frames_mjy"] = star.get("wise_single_frames_w4_mjy", [])
        elif self.archive == "neowise":
            base["neowise"] = star.get("neowise", {})
        elif self.archive in {"tess", "kepler", "ztf", "asassn"}:
            base["lightcurves"] = star.get("lightcurves", {})
        elif self.archive in {"simbad", "vizier"}:
            base["simbad_type"] = star.get("simbad_type", "Star")
            base["known_in_literature"] = star.get("known_in_literature", False)
            base["name"] = star.get("name", "")
        elif self.archive == "dustmaps":
            base["ebv_sfd"] = star.get("ebv_sfd", 0.02)
            base["ebv_planck"] = star.get("ebv_planck", 0.02)
        elif self.archive in {"radio_continuum", "deep_optical", "high_energy", "mid_far_ir", "spherex", "ephemeris"}:
            base["vetting_archive"] = star.get("vetting_archive", {})
        elif self.archive == "survey_spectra":
            base["spectra"] = star.get("spectra", {})
            base["teff_spec_k"] = star.get("teff_spec_k", star.get("teff_k"))
        elif self.archive in {"usnob1", "historical_plates"}:
            base["vanishing_data"] = star.get("vanishing_data", {})
        else:
            base["star"] = star
        return base

    def _query_live_archive(
        self,
        *,
        ids: list[str] | None,
        cone: tuple[float, float, float] | None,
        adql: str | None,
        release: str,
    ) -> list[dict[str, Any]]:
        """Execute live archive query using astroquery / pyvo (used in @pytest.mark.network runs)."""
        if self.archive == "gaia":
            from astroquery.gaia import Gaia

            table_name = "gaiadr4.gaia_source" if release == "dr4" else "gaiadr3.gaia_source"
            if adql:
                job = Gaia.launch_job(adql)
                tbl = job.get_results()
                return [dict(zip(tbl.colnames, row)) for row in tbl]
            if ids:
                id_csv = ",".join(str(x) for x in ids)
                q = f"SELECT TOP 500 source_id, ra, dec, parallax, parallax_error, pmra, pmdec, ruwe, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag FROM {table_name} WHERE source_id IN ({id_csv})"
                tbl = Gaia.launch_job(q).get_results()
                return [dict(zip(tbl.colnames, row)) for row in tbl]
        if self.archive in {"allwise", "2mass", "neowise", "wise_single_exposure"}:
            from astroquery.ipac.irsa import Irsa

            if cone is not None:
                from astropy.coordinates import SkyCoord
                import astropy.units as u

                coord = SkyCoord(ra=cone[0] * u.deg, dec=cone[1] * u.deg, frame="icrs")
                tbl = Irsa.query_region(coord, catalog=self.spec.table_or_service, radius=cone[2] * u.arcsec)
                return [dict(zip(tbl.colnames, row)) for row in tbl]
        if self.archive == "simbad" and ids:
            from astroquery.simbad import Simbad

            out = []
            for sid in ids:
                tbl = Simbad.query_object(f"Gaia DR3 {sid}")
                if tbl is not None:
                    out.extend(dict(zip(tbl.colnames, row)) for row in tbl)
            return out
        raise ArchiveAccessError(
            f"Live query requires explicit cone, ids, or ADQL supported by {self.archive}",
            retryable=False,
        )


def list_archives_info(cfg: PipelineConfig) -> list[dict[str, Any]]:
    """Return structured list of all registered archives and their policy status."""
    rows = []
    for name, spec in ARCHIVE_SPECS.items():
        pol = cfg.archives[name]
        rel = cfg.gaia.release if name == "gaia" else pol.release
        rows.append(
            {
                "archive": name,
                "enabled": pol.enabled,
                "release": rel,
                "rate_limit_rps": pol.rate_limit_rps,
                "requires_registration": pol.requires_registration,
                "access_method": spec.access_method,
                "table_or_service": spec.table_or_service,
                "description": spec.description,
            }
        )
    return rows
