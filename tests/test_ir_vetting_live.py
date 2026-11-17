"""Offline tests for live infrared vetting: every archive is a fake HTTP double (no network).

Covers the new fetchers (SkyBoT, IRSA DUST, IRSA SIA + IBE cutouts, PS1 stack, MAST CAOM), the enrichment's
`query_status` recording, and for each changed check: (a) no query -> not_evaluated, (b) failed query ->
not_evaluated naming the failure, (c) clean data -> pass with an evidence file, (d) bad data -> fail.
The reference star is Gaia DR3 5628404318169472 (K dwarf at 231.8 pc, W4 = 7.29).
"""

from __future__ import annotations

import copy
import io
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pytest
import requests

from oddstar.config import load_config
from oddstar.errors import ArchiveAccessError
from oddstar.io import cutouts as cut_mod
from oddstar.io.cache import RawCache
from oddstar.io.live_survey import _row_to_star_dict, enrich_star_live_for_vetting
from oddstar.io.physics import BAND_ZERO_POINT_JY, mag_to_flux_mjy, propagate_position
from oddstar.io.skybot import parse_skybot_text, skybot_cone_search
from oddstar.io.vetting_queries import fetch_irsa_dust, fetch_subarcsec_imaging, parse_irsa_dust_xml
from oddstar.vet import ir_checks as irc

SID = "5628404318169472"
RA, DEC, PMRA, PMDEC = 43.548272861, 5.882712239, 37.914, -97.175
PLX = 1000.0 / 231.76
NO_SLEEP = lambda s: None  # noqa: E731
MAGZP = {"W1": 20.5, "W2": 19.5, "W3": 18.0, "W4": 13.0}
NOISE_DN = {"W1": 2.0, "W2": 2.0, "W3": 1.0, "W4": 0.1}
WISE_MAGS = {"W1": 11.00, "W2": 11.02, "W3": 10.40, "W4": 7.29}
ALL_SERVICES = {
    "gaia_neighbours", "gaia_meta", "simbad", "membership", "allwise_coadd", "wise_l1b", "neowise", "second_instrument",
    "radio", "ps1_cat", "ps1_vizier", "ls", "caom", "allwise_density", "galex", "rosat", "erass1", "skybot", "dust",
    "sia", "wise_cutout", "ps1_list", "ps1_cut",
}


# ----------------------------------------------------------------------------------------- fake responses


class _Resp:
    def __init__(self, status_code: int, text: str = "", content: bytes | None = None, json_data: Any = None) -> None:
        self.status_code = status_code
        if json_data is not None and not text:
            text = json.dumps(json_data)
        self.text = text if text or content is None else content.decode("latin-1")
        self.content = content if content is not None else text.encode("utf-8")
        self._json = json_data

    def json(self) -> Any:
        if self._json is None:
            return json.loads(self.text)
        return self._json


def _ipac(cols: list[str], rows: list[dict[str, Any]]) -> str:
    w = {c: max([len(c)] + [len(str(r.get(c, "null"))) for r in rows]) + 2 for c in cols}
    lines = ["\\fixlen = T"] + ["|" + "|".join(h.rjust(w[c]) if h else c.rjust(w[c]) for c in cols) + "|" for h in ("", "char", " ", "null")]
    lines += [" " + " ".join(str(r.get(c, "null")).rjust(w[c]) for c in cols) + " " for r in rows]
    return "\n".join(lines) + "\n"


def _tsv(cols: list[str], rows: list[dict[str, Any]]) -> str:
    if not rows:
        return "#\n# VizieR: no rows in the cone\n"
    lines = ["# VizieR fake", "\t".join(cols), "\t".join("u" for _ in cols), "\t".join("-----" for _ in cols)]
    lines += ["\t".join("" if r.get(c) is None else str(r[c]) for c in cols) for r in rows]
    return "\n".join(lines) + "\n"


def _tap(cols: list[str], rows: list[list[Any]]) -> _Resp:
    return _Resp(200, json_data={"metadata": [{"name": c} for c in cols], "data": rows})


def _hms(ra_deg: float) -> str:
    h = ra_deg / 15.0
    hh, rem = int(h), (h - int(h)) * 60.0
    return f"{hh:02d} {int(rem):02d} {(rem - int(rem)) * 60.0:05.2f}"


def _dms(dec_deg: float) -> str:
    s = "-" if dec_deg < 0 else "+"
    d = abs(dec_deg)
    dd, rem = int(d), (d - int(d)) * 60.0
    return f"{s}{dd:02d} {int(rem):02d} {(rem - int(rem)) * 60.0:04.1f}"


def _wise_fits(band: str, ra_c: float, dec_c: float, size_arcsec: float, *, src_dra_arcsec: float = 0.0, src_ddec_arcsec: float = 0.0,
               fwhm_scale: float = 1.0, seed: int = 0) -> bytes:
    from astropy.io import fits
    from astropy.wcs import WCS

    pix = 1.375
    n = int(size_arcsec / pix) | 1
    hdr = fits.Header()
    hdr.update({"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRVAL1": ra_c, "CRVAL2": dec_c, "CRPIX1": (n + 1) / 2.0,
                "CRPIX2": (n + 1) / 2.0, "CDELT1": -pix / 3600.0, "CDELT2": pix / 3600.0, "MAGZP": MAGZP[band], "BAND": int(band[1])})
    wcs = WCS(hdr)
    s_ra = ra_c + src_dra_arcsec / 3600.0 / math.cos(math.radians(dec_c))
    s_dec = dec_c + src_ddec_arcsec / 3600.0
    x0, y0 = (float(v) for v in wcs.world_to_pixel_values(s_ra, s_dec))
    beam = {"W1": 6.1, "W2": 6.4, "W3": 6.5, "W4": 12.0}[band] * 1.414 * fwhm_scale
    sig = beam / 2.3548 / pix
    dn_mjy = BAND_ZERO_POINT_JY[band] * 1000.0 * 10.0 ** (-0.4 * MAGZP[band])
    amp = mag_to_flux_mjy(band, WISE_MAGS[band]) / dn_mjy / (2.0 * math.pi * sig * sig)
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    rng = np.random.default_rng(seed + int(band[1]))
    img = 100.0 + rng.normal(0.0, NOISE_DN[band], (n, n)) + amp * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2.0 * sig * sig))
    buf = io.BytesIO()
    fits.PrimaryHDU(img.astype(np.float32), hdr).writeto(buf)
    return buf.getvalue()


def _ps1_fits(ra_c: float, dec_c: float, n: int) -> bytes:
    from astropy.io import fits

    hdr = fits.Header()
    hdr.update({"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRVAL1": ra_c, "CRVAL2": dec_c, "CRPIX1": (n + 1) / 2.0,
                "CRPIX2": (n + 1) / 2.0, "CDELT1": -0.25 / 3600.0, "CDELT2": 0.25 / 3600.0})
    yy, xx = np.mgrid[0:n, 0:n].astype(float)
    c = (n - 1) / 2.0
    img = 50.0 * np.exp(-((xx - c) ** 2 + (yy - c) ** 2) / (2.0 * 2.5**2)) + np.random.default_rng(1).normal(0, 0.5, (n, n))
    buf = io.BytesIO()
    fits.PrimaryHDU(img.astype(np.float32), hdr).writeto(buf)
    return buf.getvalue()


def _sia_votable(ra: float, dec: float) -> bytes:
    from astropy.io.votable import from_table, writeto
    from astropy.table import Table

    base = "https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/04/0435/0435p060_ac51/0435p060_ac51"
    far = "https://irsa.ipac.caltech.edu/ibe/data/wise/allwise/p3am_cdd/04/0450/0450p060_ac51/0450p060_ac51"
    rows = []
    for b in (1, 2, 3, 4):
        rows.append((f"{base}-w{b}-int-3.fits", ra + 0.05, dec - 0.05))
        rows.append((f"{base}-w{b}-unc-3.fits.gz", ra + 0.05, dec - 0.05))
        rows.append((f"{far}-w{b}-int-3.fits", ra + 1.2, dec))
    t = Table(rows=rows, names=("access_url", "s_ra", "s_dec"), dtype=(str, float, float))
    buf = io.BytesIO()
    writeto(from_table(t), buf)
    return buf.getvalue()


DUST_XML = """<?xml version="1.0" encoding="UTF-8"?>
<results status="ok">
  <input><locstr>{ra} {dec} equ j2000</locstr><regSize>5.0 (deg)</regSize></input>
  <result>
    <desc>E(B-V) Reddening</desc>
    <data><image>https://irsa.ipac.caltech.edu/workspace/x/ebv.fits</image></data>
    <statistics>
      <refPixelValueSFD>{ebv} (mag)</refPixelValueSFD>
      <refPixelValueSandF>{sf} (mag)</refPixelValueSandF>
      <meanValueSFD>0.2387 (mag)</meanValueSFD>
      <meanValueSandF>0.2053 (mag)</meanValueSandF>
    </statistics>
  </result>
  <result>
    <desc>100 Micron Emission</desc>
    <statistics>
      <refPixelValue>{i100} (MJy/sr)</refPixelValue>
      <meanValue>8.12 (MJy/sr)</meanValue>
    </statistics>
  </result>
  <result><desc>Dust Temperature</desc><statistics><refPixelValue>17.9 (K)</refPixelValue></statistics></result>
</results>
"""


class FakeArchives:
    """Routes `requests.get`/`requests.post` to per-service fakes. `fail` names services that raise ConnectionError."""

    def __init__(self, *, fail: set[str] | None = None, bad: set[str] | None = None, caom_rows: list[dict[str, Any]] | None = None,
                 http500: set[str] | None = None) -> None:
        self.fail = set(fail or ())
        self.bad = set(bad or ())
        self.http500 = set(http500 or ())
        self.caom_rows = caom_rows or []
        self.calls: list[tuple[str, dict[str, Any]]] = []

    # --- classification
    def _service_get(self, url: str, p: dict[str, Any]) -> str:
        if "Gator" in url:
            cat = p.get("catalog")
            if cat == "allwise_p3as_psd":
                return "allwise_density" if "selcols" in p else "allwise_coadd"
            return {"allsky_4band_p1bs_psd": "wise_l1b", "neowiser_p1bs_psd": "neowise"}.get(str(cat), "second_instrument")
        if "asu-tsv" in url:
            return {"VIII/65/nvss": "radio", "II/349/ps1": "ps1_vizier", "II/335/galex_ais": "galex", "IX/47/2rxs": "rosat",
                    "J/A+A/682/A34/erass1-m": "erass1"}[p["-source"]]
        if "panstarrs" in url:
            return "ps1_cat"
        if "api/v0/invoke" in url:
            return "caom"
        if "skybot" in url:
            return "skybot"
        if "DUST" in url:
            return "dust"
        if url == cut_mod.IRSA_SIA_URL:
            return "sia"
        if "-int-3" in url:
            return "wise_cutout"
        if "ps1filenames" in url:
            return "ps1_list"
        if "fitscut" in url:
            return "ps1_cut"
        raise AssertionError(f"unrouted GET {url}")

    def _service_post(self, url: str, q: str) -> str:
        if "h_link" in q:
            return "membership"
        if "FROM basic" in q:
            return "simbad"
        if "non_single_star" in q:
            return "gaia_meta"
        if "gaiadr3.gaia_source" in q:
            return "gaia_neighbours"
        if "ls_dr10" in q:
            return "ls"
        if "allwise" in q:
            return "allwise_coadd"
        if "allsky_4band" in q:
            return "wise_l1b"
        if "neowiser" in q:
            return "neowise"
        if "nvss" in q:
            return "radio"
        raise AssertionError(f"unrouted POST {url}: {q[:80]}")

    def _guard(self, svc: str) -> _Resp | None:
        if svc in self.fail:
            raise requests.exceptions.ConnectionError(f"simulated outage of {svc}")
        if svc in self.http500:
            return _Resp(503, text="Service Unavailable")
        return None

    # --- GET
    def get(self, url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        p = dict(params or {})
        svc = self._service_get(url, p)
        self.calls.append((svc, p))
        g = self._guard(svc)
        if g is not None:
            return g
        bad = svc in self.bad
        if svc == "allwise_coadd":
            p_ra, p_dec = propagate_position(RA, DEC, PMRA, PMDEC, 2016.0, 2010.2)
            row = {"ra": f"{p_ra:.7f}", "dec": f"{p_dec:.7f}", "cc_flags": "00D0" if bad else "0000", "ext_flg": 0, "ph_qual": "AAAA",
                   "w4rchi2": 0.93, "w4m": 12, "w4nm": 12, "w3nm": 12, "w4snr": 14.9, "w3snr": 41.0}
            return _Resp(200, _ipac(list(row), [row]))
        if svc == "wise_l1b":
            rows = []
            for i, mjd0 in enumerate((55216.10, 55398.40)):
                for k in range(8):
                    rows.append({"mjd": f"{mjd0 + 0.0662 * k:.5f}", "w4mpro": f"{7.29 + 0.03 * ((k + i) % 3 - 1):.3f}", "w4sigmpro": "0.110", "w4snr": "9.8"})
            return _Resp(200, _ipac(["mjd", "w4mpro", "w4sigmpro", "w4snr"], rows))
        if svc in {"neowise", "second_instrument"}:
            return _Resp(200, _ipac(["ra", "dec", "mjd"], []))
        if svc == "allwise_density":
            rows = [{"ra": f"{RA:.6f}", "dec": f"{DEC:.6f}", "w3mpro": "10.400", "w3snr": "41", "w4mpro": "7.290", "w4snr": "14.9"}]
            n_bright = 900 if bad else 2
            rng = np.random.default_rng(7)
            for k in range(n_bright + 40):
                r = 30.0 + 260.0 * rng.random()
                th = 2.0 * math.pi * rng.random()
                bright = k < n_bright
                rows.append({"ra": f"{RA + r * math.cos(th) / 3600.0 / math.cos(math.radians(DEC)):.6f}", "dec": f"{DEC + r * math.sin(th) / 3600.0:.6f}",
                             "w3mpro": "9.000" if bright else "12.500", "w3snr": "20" if bright else "1.5",
                             "w4mpro": "6.500" if bright else "9.200", "w4snr": "8" if bright else "1.2"})
            return _Resp(200, _ipac(["ra", "dec", "w3mpro", "w3snr", "w4mpro", "w4snr"], rows))
        if svc == "radio":
            rows = [{"RAJ2000": f"{RA + 2.0 / 3600.0:.6f}", "DEJ2000": f"{DEC:.6f}", "S1.4": "4.1"}] if bad else []
            return _Resp(200, _tsv(["RAJ2000", "DEJ2000", "S1.4"], rows))
        if svc == "ps1_vizier":
            return _Resp(200, _tsv(["RAJ2000", "DEJ2000", "rmag"], []))
        if svc == "galex":
            c_ra, c_dec = (float(v) for v in p["-c"].split())
            row = {"RAJ2000": f"{c_ra:.6f}", "DEJ2000": f"{c_dec + 1.9 / 3600.0:.6f}", "NUVmag": "18.00" if bad else "22.44",
                   "e_NUVmag": "0.05" if bad else "0.43", "FUVmag": None}
            return _Resp(200, _tsv(["RAJ2000", "DEJ2000", "NUVmag", "e_NUVmag", "FUVmag"], [row]))
        if svc == "rosat":
            c_ra, c_dec = (float(v) for v in p["-c"].split())
            rows = [{"RAJ2000": f"{c_ra:.5f}", "DEJ2000": f"{c_dec + 8.0 / 3600.0:.5f}", "Rate": "0.10", "HR1": "0.0"}] if bad else []
            return _Resp(200, _tsv(["RAJ2000", "DEJ2000", "Rate", "HR1"], rows))
        if svc == "erass1":
            raise AssertionError("eRASS1 must not be queried at Galactic l < 180")
        if svc == "ps1_cat":
            return _Resp(200, json_data={
                "info": [{"name": n} for n in ("raMean", "decMean", "nDetections", "rMeanPSFMag", "rMeanKronMag")],
                "data": [[RA, DEC, 61, 13.21, 13.20]],
            })
        if svc == "caom":
            req = json.loads(p["request"])
            assert req["service"] == "Mast.Caom.Filtered.Position"
            return _Resp(200, json_data={"status": "COMPLETE", "data": self.caom_rows})
        if svc == "skybot":
            if bad:
                ra_c, dec_c = float(p["-ra"]), float(p["-dec"])
                txt = ("# Flag: 1\n# Ticket: 1\n"
                       "# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec) | dRA(arcsec/h) | dDEC(arcsec/h) | Dg(ua) | Dh(ua)\n"
                       f"99999 | Testbody | {_hms(ra_c)} | {_dms(dec_c)} | MB>Middle | 17.9 | 0.15 | 0.4 | 0.0 | 0.0 | 1.6 | 2.5\n")
                return _Resp(200, txt)
            return _Resp(200, "# Flag: 0\n# Ticket: 1\n")
        if svc == "dust":
            return _Resp(200, DUST_XML.format(ra=RA, dec=DEC, ebv="0.6120" if bad else "0.2414", sf="0.5263" if bad else "0.2076",
                                              i100="21.3" if bad else "8.25"))
        if svc == "sia":
            ra_s, dec_s = (float(v) for v in p["POS"].split()[1:3])
            return _Resp(200, content=_sia_votable(ra_s, dec_s))
        if svc == "wise_cutout":
            c_ra, c_dec = (float(v) for v in p["center"].removesuffix("deg").split(","))
            band = "W" + url.split("-int-3")[0][-1]
            off = 6.0 if (bad and band in {"W3", "W4"}) else 0.0
            return _Resp(200, content=_wise_fits(band, c_ra, c_dec, float(p["size"].removesuffix("arcsec")), src_ddec_arcsec=off))
        if svc == "ps1_list":
            return _Resp(200, "projcell subcell ra dec filter mjd type filename shortname badflag\n"
                              "1234 056 43.548 5.883 r 0 stack /rings.v3.skycell/1234/056/rings.v3.skycell.1234.056.stk.r.unconv.fits "
                              "rings.v3.skycell.1234.056.stk.r.unconv.fits 0\n")
        if svc == "ps1_cut":
            return _Resp(200, content=_ps1_fits(float(p["ra"]), float(p["dec"]), int(p["size"])))
        raise AssertionError(svc)

    # --- POST (TAP)
    def post(self, url: str, data: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        q = str((data or {}).get("QUERY", ""))
        svc = self._service_post(url, q)
        self.calls.append((svc, {"url": url, "QUERY": q}))
        g = self._guard(svc)
        if g is not None:
            return g
        bad = svc in self.bad
        if svc == "gaia_neighbours":
            rows = [[int(SID), RA, DEC, 13.43, 12.69]]
            if bad:
                rows.append([5628404318169999, RA, DEC + 5.0 / 3600.0, 10.1, 9.0])
            return _tap(["source_id", "ra", "dec", "phot_g_mean_mag", "phot_rp_mean_mag"], rows)
        if svc == "gaia_meta":
            return _tap(["non_single_star", "radial_velocity_error", "rv_nb_transits", "rv_method_used", "rv_chisq_pvalue",
                         "rv_amplitude_robust", "phot_variable_flag"], [[1 if bad else 0, 1.1, 14, 1, 0.42, 3.0, "NOT_AVAILABLE"]])
        if svc == "simbad":
            return _tap(["main_id", "otype_txt", "sep_arcsec"], [[f"Gaia DR3 {SID}", "*", 0.1]])
        if svc == "membership":
            return _tap(["main_id", "otype", "membership"], [])
        if svc == "ls":
            if bad:
                return _Resp(200, "ls_id,ra,dec,type,flux_g,flux_r,flux_z,flux_w3,flux_w4,flux_ivar_w4\n"
                                  f"7,{RA + 1.0 / 3600.0:.6f},{DEC:.6f},EXP,2.0,5.0,9.0,900.0,2500.0,0.01\n")
            return _Resp(200, "ls_id,ra,dec,type,flux_g,flux_r,flux_z,flux_w3,flux_w4,flux_ivar_w4\n")
        raise AssertionError(svc)


# -------------------------------------------------------------------------------------------- fixtures


def _gaia_row() -> dict[str, Any]:
    w_ra, w_dec = propagate_position(RA, DEC, PMRA, PMDEC, 2016.0, 2010.2)
    return {
        "source_id": SID, "ra": RA, "dec": DEC, "b": -45.6, "parallax": PLX, "parallax_over_error": 85.0, "pmra": PMRA, "pmdec": PMDEC,
        "ruwe": 1.02, "phot_g_mean_mag": 13.43, "phot_bp_mean_mag": 14.04, "phot_rp_mean_mag": 12.69, "teff_gspphot": 4510.0,
        "logg_gspphot": 4.6, "mh_gspphot": -0.1, "ebpminrp_gspphot": 0.07, "j_m": 11.67, "h_m": 11.18, "ks_m": 11.08,
        "w1mpro": WISE_MAGS["W1"], "w2mpro": WISE_MAGS["W2"], "w3mpro": WISE_MAGS["W3"], "w4mpro": WISE_MAGS["W4"],
        "w1mpro_error": 0.023, "w2mpro_error": 0.021, "w3mpro_error": 0.06, "w4mpro_error": 0.073,
        "cc_flags": "0000", "ext_flag": 0, "ph_qual": "AAAA", "wise_ra": w_ra, "wise_dec": w_dec, "allwise_n_neighbours": 1,
    }


def _enrich(tmp: Path, fake: FakeArchives) -> dict[str, Any]:
    cfg, _, _ = load_config()
    cache = RawCache(tmp / "raw", max_cache_gb=1.0)
    star = _row_to_star_dict(_gaia_row())
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(requests, "get", fake.get)
        mp.setattr(requests, "post", fake.post)
        out, mans = enrich_star_live_for_vetting(star, cfg, cache, enrich_optical=False, sleep=NO_SLEEP)
    if not fake.fail:
        assert mans, "every successful query must leave a manifest"
    return out


@pytest.fixture(scope="module")
def cfg() -> Any:
    return load_config()[0]


@pytest.fixture(scope="module")
def ref_star(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return _enrich(tmp_path_factory.mktemp("ref"), FakeArchives())


@pytest.fixture(scope="module")
def bad_star(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    bad = {"gaia_neighbours", "gaia_meta", "allwise_coadd", "radio", "ls", "allwise_density", "galex", "rosat", "skybot", "dust", "wise_cutout"}
    return _enrich(tmp_path_factory.mktemp("bad"), FakeArchives(bad=bad))


@pytest.fixture(scope="module")
def failed_star(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    return _enrich(tmp_path_factory.mktemp("failed"), FakeArchives(fail=ALL_SERVICES))


def _run(check: str, star: dict[str, Any], cfg: Any, out: Path) -> Any:
    return getattr(irc, f"check_{check}")(star, cfg, out)


def _assert_pass_with_evidence(res: Any) -> None:
    assert res.status == "pass", res.note
    assert res.evidence is not None and Path(res.evidence).exists()


# ------------------------------------------------------------------------------------------- fetchers


def test_skybot_parse_and_search(tmp_path: Path, monkeypatch: Any) -> None:
    txt = ("# Flag: 1\n# Ticket: 1\n"
           "# Num | Name | RA(h) | DE(deg) | Class | Mv | Err(arcsec) | d(arcsec) | dRA(arcsec/h) | dDEC(arcsec/h) | Dg(ua) | Dh(ua)\n"
           "4321 | Testbody | 02 54 11.59 | +05 52 57.8 | MB>Inner | 18.5 | 0.2 | 3.1 | 20.5 | -5.25 | 1.5 | 2.4\n")
    bodies = parse_skybot_text(txt, ra_deg=43.548, dec_deg=5.8827, epoch_jd=2455216.6)
    assert len(bodies) == 1
    b = bodies[0]
    assert {"name", "ra_deg", "dec_deg", "sep_arcsec", "v_mag"} <= set(b)
    assert b["name"] == "Testbody" and b["v_mag"] == 18.5 and b["dra_arcsec_hr"] == 20.5 and b["ddec_arcsec_hr"] == -5.25
    assert abs(b["ra_deg"] - 43.54829) < 1e-4 and abs(b["dec_deg"] - 5.88272) < 1e-4
    assert parse_skybot_text("# Flag: 0\n", ra_deg=0, dec_deg=0, epoch_jd=2455000.0) == []
    with pytest.raises(ArchiveAccessError):
        parse_skybot_text("# Flag: -1\n# error: bad epoch\n", ra_deg=0, dec_deg=0, epoch_jd=2455000.0)
    with pytest.raises(ArchiveAccessError):
        parse_skybot_text("<html>maintenance</html>", ra_deg=0, dec_deg=0, epoch_jd=2455000.0)

    seen: list[dict[str, Any]] = []

    def fake_get(url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        seen.append(dict(params or {}))
        return _Resp(200, txt)

    monkeypatch.setattr(requests, "get", fake_get)
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    out, man = skybot_cone_search(ra_deg=43.548, dec_deg=5.8827, radius_arcsec=60.0, epoch_jd=2455216.6, observer="C51", cache=cache, sleep=NO_SLEEP)
    assert out[0]["name"] == "Testbody" and man.archive == "ephemeris"
    assert seen[0]["-loc"] == "C51" and abs(float(seen[0]["-rd"]) - 60.0 / 3600.0) < 1e-6
    # Cached: second call makes no request
    skybot_cone_search(ra_deg=43.548, dec_deg=5.8827, radius_arcsec=60.0, epoch_jd=2455216.6, observer="C51", cache=cache, sleep=NO_SLEEP)
    assert len(seen) == 1

    def flaky(url: str, params: dict[str, Any] | None = None, timeout: float = 0) -> _Resp:
        return _Resp(503, "busy")

    monkeypatch.setattr(requests, "get", flaky)
    with pytest.raises(ArchiveAccessError):
        skybot_cone_search(ra_deg=10.0, dec_deg=0.0, radius_arcsec=30.0, epoch_jd=2455300.0, cache=cache, sleep=NO_SLEEP)


def test_irsa_dust_parse_and_fetch(tmp_path: Path, monkeypatch: Any) -> None:
    vals = parse_irsa_dust_xml(DUST_XML.format(ra=RA, dec=DEC, ebv="0.2414", sf="0.2076", i100="8.25"))
    assert vals["ebv_sfd_los"] == 0.2414 and vals["ebv_sf11_los"] == 0.2076 and vals["i100_mjy_sr"] == 8.25
    assert vals["ebv_sfd_region_mean"] == 0.2387
    with pytest.raises(ArchiveAccessError):
        parse_irsa_dust_xml('<results status="error"><message>Invalid location</message></results>')
    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(200, DUST_XML.format(ra=RA, dec=DEC, ebv="0.2414", sf="0.2076", i100="8.25")))
    d, man = fetch_irsa_dust(ra_deg=RA, dec_deg=DEC, cache=RawCache(tmp_path / "raw", max_cache_gb=1.0), sleep=NO_SLEEP)
    assert d["ebv_sfd_los"] == 0.2414 and man.archive == "dustmaps"


def test_mast_caom_subarcsec_listing(tmp_path: Path, monkeypatch: Any) -> None:
    rows = [{"obs_collection": "HST", "instrument_name": "ACS/WFC", "filters": "F814W", "obs_id": "jabc01010", "t_exptime": 2200.0, "t_min": 57000.1}]
    fake = FakeArchives(caom_rows=rows)
    monkeypatch.setattr(requests, "get", fake.get)
    obs, _ = fetch_subarcsec_imaging(ra_deg=RA, dec_deg=DEC, radius_arcsec=3.0, cache=RawCache(tmp_path / "raw", max_cache_gb=1.0), sleep=NO_SLEEP)
    assert obs == [{"collection": "HST", "instrument": "ACS/WFC", "filters": "F814W", "obs_id": "jabc01010", "exptime_s": 2200.0, "t_min_mjd": 57000.1}]


def test_allwise_cutouts_fetch_and_measure(tmp_path: Path, monkeypatch: Any, cfg: Any) -> None:
    fake = FakeArchives()
    monkeypatch.setattr(requests, "get", fake.get)
    cache = RawCache(tmp_path / "raw", max_cache_gb=1.0)
    images, mans, log = cut_mod.fetch_allwise_cutouts(ra=RA, dec=DEC, size_arcsec=200.0, cache=cache, sleep=NO_SLEEP)
    assert set(images) == {"W1", "W2", "W3", "W4"}, log
    # The nearest tile's intensity image is used (not the unc image, not the far tile)
    cut_urls = [p for s, p in fake.calls if s == "wise_cutout"]
    assert len(cut_urls) == 4 and all(p["size"] == "200arcsec" for p in cut_urls)
    m = cut_mod.measure_wise_cutouts(images, ra_pred=RA, dec_pred=DEC, vetting=cfg.vetting)
    for b in ("W3", "W4"):
        assert m[b]["status"] == "ok"
        assert m[b]["offset_arcsec"] < 1.0 and m[b]["offset_sigma_arcsec"] < 1.0
        assert 0.85 < m[b]["fwhm_ratio"] < 1.15
        assert m[b]["flux_snr"] > 20 and m[b]["bg_rms_mjy_beam"] < 1.0
        assert not m[b]["brighter_secondary_peak"]
        assert Path(m[b]["blob_path"]).exists()
    assert abs(m["W4"]["flux_mjy"] - mag_to_flux_mjy("W4", 7.29)) / mag_to_flux_mjy("W4", 7.29) < 0.1
    # An offset source is measured as offset
    off = cut_mod.read_fits_image(_wise_fits("W4", RA, DEC, 200.0, src_ddec_arcsec=6.0))
    mo = cut_mod.measure_wise_cutouts({"W4": {"data": off[0], "header": off[1]}}, ra_pred=RA, dec_pred=DEC, vetting=cfg.vetting)
    assert abs(mo["W4"]["offset_arcsec"] - 6.0) < 1.0
    # SIA garbage is a failure, not an empty result
    monkeypatch.setattr(requests, "get", lambda url, params=None, timeout=0: _Resp(200, "<html>oops</html>"))
    with pytest.raises(ArchiveAccessError):
        cut_mod.fetch_allwise_cutouts(ra=10.0, dec=10.0, size_arcsec=200.0, cache=cache, sleep=NO_SLEEP)


# ----------------------------------------------------------------------------------- enrichment status


def test_enrichment_records_query_status(ref_star: dict[str, Any], failed_star: dict[str, Any]) -> None:
    qs = ref_star["vetting_archive"]["query_status"]
    for name in ("gaia_neighbours", "gaia_meta", "simbad", "allwise_coadd", "wise_l1b", "radio_nvss", "deep_optical_ps1",
                 "subarcsec_imaging_mast", "allwise_density", "galex", "xray_2rxs", "young_association", "skybot", "irsa_dust",
                 "cutouts_wise", "cutouts_ps1", "second_instrument"):
        assert qs[name] == "ok", (name, qs[name])
    assert qs["xray_erass1"] == "no_coverage"
    assert qs["deep_optical_ls_dr10"] == "no_coverage" and qs["neowise"] == "no_coverage"
    fq = failed_star["vetting_archive"]["query_status"]
    for name in ("gaia_neighbours", "gaia_meta", "allwise_coadd", "wise_l1b", "radio_nvss", "deep_optical_ls_dr10", "subarcsec_imaging_mast",
                 "allwise_density", "galex", "xray_2rxs", "young_association", "skybot", "irsa_dust", "cutouts_wise", "cutouts_ps1"):
        assert fq[name].startswith("failed: "), (name, fq[name])
    assert "ConnectionError" in fq["irsa_dust"] and "ConnectionError" in fq["cutouts_wise"]
    # A failed query writes no vetting input
    va = failed_star["vetting_archive"]
    for key in ("neighbours_in_w4_beam", "nss_flag", "radio_sources_within_10arcsec", "chance_alignment", "dust", "sso_search"):
        assert key not in va, key
    assert "w4rchi2" not in failed_star["allwise"]


# --------------------------------------------------------------------------------- reference regression


def test_reference_star_regression(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    res = {r.check: r for r in irc.run_ir_vetting(ref_star, cfg, tmp_path / "ev")}
    status = {k: v.status for k, v in res.items()}
    # Chance alignment is measured from the AllWISE annulus, not a default density
    ca = res["chance_alignment"]
    assert ca.status == "pass" and "Measured AllWISE density" in ca.note
    assert ca.metrics["bands"]["W4"]["n_sources"] == 2 and ca.metrics["chance_alignment_prob"] < 0.01
    # Cirrus reads the line-of-sight IRSA value, not the stellar ebv_sfd
    cir = res["cirrus"]
    assert cir.metrics["ebv_sfd_los"] == 0.2414 and cir.metrics["ebv_sf11_los"] == 0.2076 and cir.metrics["i100_mjy_sr"] == 8.25
    assert cir.status == "pass" and "0.2414" in cir.note
    assert ref_star["vetting_archive"]["dust"]["ebv_sfd_los"] == 0.2414
    # The Gaia per-star reddening lives in ebv_star; ebv_sfd is reserved for true line-of-sight map values
    assert ref_star["ebv_sfd"] is None and ref_star["ebv_star"] == pytest.approx(0.07 / 1.31)
    assert ref_star["ebv_star_source"] == "gaia_gspphot"
    # ROSAT non-detection at 232 pc cannot reach 1e29 erg/s; GALEX NUV 22.44+/-0.43 is photospheric
    y = res["youth"]
    assert y.status == "not_evaluated"
    ros = y.metrics["xray"]["surveys"]["rosat_2rxs"]
    assert ros["verdict"] == "non_constraining" and 6.0e29 < ros["lx_upper_limit_erg_s"] < 7.0e29
    assert y.metrics["xray"]["surveys"]["erass1"]["verdict"] == "no_coverage"
    assert y.metrics["uv"]["verdict"] == "not_young" and abs(y.metrics["uv"]["sep_arcsec"] - 1.9) < 0.2
    # PS1 alone cannot clear deep_imaging
    di = res["deep_imaging"]
    assert di.status == "not_evaluated" and di.metrics["screened_surveys"] == ["PS1 DR2 mean catalogue (MAST)"]
    assert "cannot clear the beam" in di.note
    # Real-pixel checks
    assert status["centroid_offset"] == "pass" and res["centroid_offset"].metrics["source"].startswith("W3/W4 Atlas image")
    assert status["vision_rubric"] == "pass"
    png = Path(res["vision_rubric"].metrics["cutout_png"])
    assert png.exists() and png.stat().st_size > 5000
    assert status["moving_object"] == "pass" and len(res["moving_object"].metrics["visits"]) == 2
    assert status["second_instrument"] == "not_evaluated"
    for chk in ("xmatch_integrity", "wise_quality", "neighbours", "radio_counterpart", "extragalactic_colours", "single_frames", "evolved", "binary"):
        assert status[chk] == "pass", (chk, res[chk].note)
    for r in res.values():
        if r.status in {"pass", "fail"}:
            assert r.evidence and Path(r.evidence).exists()
        json.dumps(r.metrics)  # evidence metrics must be JSON-serialisable


def test_all_queries_failed_gives_not_evaluated_with_reason(failed_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    res = {r.check: r for r in irc.run_ir_vetting(failed_star, cfg, tmp_path / "ev")}
    for chk in ("wise_quality", "centroid_offset", "neighbours", "radio_counterpart", "deep_imaging", "chance_alignment", "cirrus",
                "single_frames", "moving_object", "youth", "binary", "vision_rubric"):
        assert res[chk].status == "not_evaluated", (chk, res[chk].status, res[chk].note)
        assert "failed: " in res[chk].note, (chk, res[chk].note)
    assert "ConnectionError" in res["cirrus"].note and "irsa_dust=failed" in res["cirrus"].note
    # The catalogue offset is still reported when no image centroid is available
    assert res["centroid_offset"].metrics["catalogue_offset_arcsec"] is not None
    # Nothing passes on missing data except checks evaluated from the Gaia/AllWISE join row itself
    assert {k for k, v in res.items() if v.status == "pass"} <= {"xmatch_integrity", "extragalactic_colours", "evolved"}


def test_bad_data_fails(bad_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    res = {r.check: r for r in irc.run_ir_vetting(bad_star, cfg, tmp_path / "ev")}
    for chk in ("wise_quality", "centroid_offset", "neighbours", "radio_counterpart", "deep_imaging", "chance_alignment", "cirrus",
                "moving_object", "youth", "binary", "vision_rubric"):
        assert res[chk].status == "fail", (chk, res[chk].status, res[chk].note)
        assert res[chk].evidence and Path(res[chk].evidence).exists()
    assert res["chance_alignment"].metrics["chance_alignment_prob"] > 0.01
    assert "Testbody" in res["moving_object"].note
    assert res["cirrus"].metrics["ebv_sfd_los"] == 0.612
    assert res["youth"].metrics["uv"]["verdict"] == "young" and res["youth"].metrics["xray"]["verdict"] == "young"
    assert "LS DR10" in res["deep_imaging"].note
    assert res["neighbours"].metrics["blend_flux_fraction"] > 0.25


# --------------------------------------------------------------------- per-check (a) no query, edge cases


NO_QUERY = {
    # check: (keys to delete from vetting_archive, query_status names to delete)
    "neighbours": (["neighbours_in_w4_beam"], ["gaia_neighbours"]),
    "radio_counterpart": (["radio_sources_within_10arcsec"], ["radio_nvss"]),
    "chance_alignment": (["chance_alignment"], ["allwise_density"]),
    "cirrus": (["dust", "cutouts"], ["irsa_dust", "cutouts_wise"]),
    "moving_object": (["sso_search"], ["skybot"]),
    "youth": (["youth"], ["galex", "xray_2rxs", "xray_erass1", "young_association"]),
    "vision_rubric": (["cutouts"], ["cutouts_wise", "cutouts_ps1"]),
    "centroid_offset": (["cutouts"], ["cutouts_wise"]),
    "binary": (["nss_flag", "rv_scatter_kms"], ["gaia_meta"]),
    "deep_imaging": (["neighbour_screen", "subarcsec_imaging"], ["deep_optical_ps1", "deep_optical_ls_dr10", "subarcsec_imaging_mast"]),
    "second_instrument": (["second_instrument_ir"], ["second_instrument"]),
}


@pytest.mark.parametrize("check", sorted(NO_QUERY))
def test_no_query_is_not_evaluated(check: str, ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    keys, queries = NO_QUERY[check]
    for k in keys:
        star["vetting_archive"].pop(k, None)
    for q in queries:
        star["vetting_archive"]["query_status"].pop(q, None)
    res = _run(check, star, cfg, tmp_path)
    assert res.status == "not_evaluated", res.note
    assert res.evidence is None
    if check != "second_instrument":
        assert "not run" in res.note


def test_flat_inputs_ignored_when_query_failed(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    """A stale legacy value must not be trusted once its producing query is recorded as failed."""
    star = copy.deepcopy(ref_star)
    star["vetting_archive"]["query_status"]["radio_nvss"] = "failed: ReadTimeout: VizieR timed out"
    assert star["vetting_archive"]["radio_sources_within_10arcsec"] == []
    r = irc.check_radio_counterpart(star, cfg, tmp_path)
    assert r.status == "not_evaluated" and "ReadTimeout" in r.note
    star["vetting_archive"]["query_status"]["gaia_meta"] = "failed: HTTPError: 500"
    assert irc.check_binary(star, cfg, tmp_path).status == "not_evaluated"


def test_wise_quality_needs_profile_fit_but_fails_on_flags(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    for k in ("w4rchi2", "w4m", "w4nm"):
        star["allwise"].pop(k)
    star["vetting_archive"]["query_status"]["allwise_coadd"] = "failed: ConnectionError: x"
    r = irc.check_wise_quality(star, cfg, tmp_path)
    assert r.status == "not_evaluated" and "allwise_coadd=failed" in r.note
    star["allwise"]["cc_flags"] = "000H"
    assert irc.check_wise_quality(star, cfg, tmp_path).status == "fail"


def test_xmatch_needs_neighbour_count(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    _assert_pass_with_evidence(irc.check_xmatch_integrity(star, cfg, tmp_path))
    star["allwise"]["n_neighbours_gaia_xmatch"] = 2
    assert irc.check_xmatch_integrity(star, cfg, tmp_path).status == "fail"
    star["allwise"].pop("n_neighbours_gaia_xmatch")
    assert irc.check_xmatch_integrity(star, cfg, tmp_path).status == "not_evaluated"


def test_centroid_catalogue_only_is_not_evaluated_and_reports_offset(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["vetting_archive"]["cutouts"]["bands"] = {}
    star["allwise"]["centroid_offset_arcsec"] = 9.0  # even a large catalogue offset is not a W4 centroid
    r = irc.check_centroid_offset(star, cfg, tmp_path)
    assert r.status == "not_evaluated" and r.metrics["catalogue_offset_arcsec"] == 9.0 and "catalogue" in r.note


def test_centroid_non_constraining_when_sigma_large(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    for b in ("W3", "W4"):
        star["vetting_archive"]["cutouts"]["bands"][b]["offset_sigma_arcsec"] = 3.0  # 3 sigma = 9" > 0.35 beam
    r = irc.check_centroid_offset(star, cfg, tmp_path)
    assert r.status == "not_evaluated" and r.metrics["bands"]["W4"]["verdict"] == "non_constraining"


def test_chance_alignment_needs_excess(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["flux_mjy"].pop("Ks")
    r = irc.check_chance_alignment(star, cfg, tmp_path)
    assert r.status == "not_evaluated" and "excess" in r.note


def test_youth_passes_only_when_both_constraining(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["distance_pc"] = 60.0  # ROSAT limit 1e-13 -> 4.3e28 erg/s < 1e29
    r = irc.check_youth(star, cfg, tmp_path)
    _assert_pass_with_evidence(r)
    assert r.metrics["xray"]["surveys"]["rosat_2rxs"]["verdict"] == "not_young"
    # A GALEX non-detection leaves the UV side unconstrained
    star["vetting_archive"]["youth"]["galex"]["match"] = None
    assert irc.check_youth(star, cfg, tmp_path).status == "not_evaluated"
    # A young-group membership fails on its own
    star["vetting_archive"]["youth"]["association"]["young_groups"] = [{"main_id": "Test group", "otype": "MGr", "membership": 90}]
    assert irc.check_youth(star, cfg, tmp_path).status == "fail"


def test_deep_imaging_lists_hst_but_does_not_pass(tmp_path: Path, cfg: Any) -> None:
    rows = [{"obs_collection": "HST", "instrument_name": "WFC3/IR", "filters": "F160W", "obs_id": "ib1234010", "t_exptime": 1200.0, "t_min": 56000.0}]
    star = _enrich(tmp_path, FakeArchives(caom_rows=rows))
    r = irc.check_deep_imaging(star, cfg, tmp_path / "ev")
    assert r.status == "not_evaluated" and "ib1234010" in r.note and "must be inspected" in r.note
    # Only an actual assessment of the image can pass ...
    star["vetting_archive"]["subarcsec_imaging"]["assessment"] = {"assessed": True, "observation": "ib1234010", "assessed_by": "reviewer",
                                                                   "contaminant": None}
    _assert_pass_with_evidence(irc.check_deep_imaging(star, cfg, tmp_path / "ev"))
    # ... or fail
    star["vetting_archive"]["subarcsec_imaging"]["assessment"]["contaminant"] = "resolved galaxy 0.8\" NE"
    assert irc.check_deep_imaging(star, cfg, tmp_path / "ev").status == "fail"


def test_moving_object_failed_visit_is_not_evaluated(tmp_path: Path, cfg: Any) -> None:
    star = _enrich(tmp_path, FakeArchives(http500={"skybot"}))
    r = irc.check_moving_object(star, cfg, tmp_path / "ev")
    assert r.status == "not_evaluated" and "HTTP 503" in r.note
    assert star["vetting_archive"]["query_status"]["skybot"].startswith("failed: ")


def test_cirrus_partial_inputs(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["vetting_archive"]["cutouts"]["bands"] = {}
    r = irc.check_cirrus(star, cfg, tmp_path)
    assert r.status == "not_evaluated" and "measured W3/W4 background" in r.note
    star["vetting_archive"]["dust"]["ebv_sfd_los"] = 0.4  # one input over the limit is enough to fail
    assert irc.check_cirrus(star, cfg, tmp_path).status == "fail"


def test_neighbour_flux_estimate_uses_colour_assumption(bad_star: dict[str, Any], cfg: Any) -> None:
    nb = bad_star["vetting_archive"]["neighbours_in_w4_beam"][0]
    assert nb["flux_frac_basis"] == "RP with assumed M-dwarf colours"
    assert nb["w4_est_mjy"] == pytest.approx(mag_to_flux_mjy("W4", 9.0 - cfg.vetting["neighbour_assumed_rp_minus_w4"]), rel=1e-4)
    assert nb["flux_frac_w4"] == pytest.approx(nb["w4_est_mjy"] / mag_to_flux_mjy("W4", 7.29), rel=1e-3)


def test_extragalactic_colours_reports_excess_colours(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    r = irc.check_extragalactic_colours(ref_star, cfg, tmp_path)
    _assert_pass_with_evidence(r)
    assert r.metrics["excess_w3_w4"] is not None and "hot circumstellar dust" in r.note
    star = copy.deepcopy(ref_star)
    star["mags"]["W1"], star["mags"]["W2"], star["mags"]["W3"] = 12.0, 11.0, 8.0
    assert irc.check_extragalactic_colours(star, cfg, tmp_path).status == "fail"


def test_vision_rubric_extended_source_fails(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["vetting_archive"]["cutouts"]["bands"]["W4"]["fwhm_ratio"] = 1.8
    r = irc.check_vision_rubric(star, cfg, tmp_path)
    assert r.status == "fail" and "extended" in r.note
    star["vetting_archive"]["cutouts"]["bands"]["W4"]["fwhm_ratio"] = 1.0
    star["vetting_archive"]["cutouts"]["bands"]["W4"]["flux_snr"] = 2.0
    assert irc.check_vision_rubric(star, cfg, tmp_path).status == "not_evaluated"


def test_second_instrument_states(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["vetting_archive"]["second_instrument_ir"] = {"instrument": "Spitzer_MIPS24", "n_matches": 1, "interpretation": "catalogue_match_unverified"}
    assert irc.check_second_instrument(star, cfg, tmp_path).status == "not_evaluated"
    star["vetting_archive"]["second_instrument_ir"]["interpretation"] = "no_excess_at_second_instrument"
    assert irc.check_second_instrument(star, cfg, tmp_path).status == "fail"


def test_evolved_with_default_teff_is_not_evaluated(ref_star: dict[str, Any], cfg: Any, tmp_path: Path) -> None:
    star = copy.deepcopy(ref_star)
    star["abs_g_mag"], star["teff_k"], star["teff_k_source"] = 0.0, 5500.0, "default_5500K"
    assert irc.check_evolved(star, cfg, tmp_path).status == "not_evaluated"
    star["teff_k"], star["teff_k_source"] = 4600.0, "gspphot"
    assert irc.check_evolved(star, cfg, tmp_path).status == "fail"


def test_fixture_cutouts_are_labelled_synthetic(cfg: Any, tmp_path: Path) -> None:
    from oddstar.dossier.plots import render_cutouts_panel
    from oddstar.io.fixtures_data import generate_field_sample

    fx = generate_field_sample(n_stars=3, seed=1)[0]
    assert fx["synthetic_fixture"] is True
    r = irc.check_vision_rubric(fx, cfg, tmp_path)
    assert r.metrics["synthetic_fixture"] is True and "SYNTHETIC" in r.note
    # A real star without fetched pixels gets a placeholder, never a synthetic sky image
    star = _row_to_star_dict(_gaia_row())
    p = render_cutouts_panel(star, tmp_path / "placeholder.png")
    assert p.exists()
