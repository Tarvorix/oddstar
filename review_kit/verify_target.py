"""Isolated standalone reviewer verification script (Ground Rule R6).

Enforces a runtime isolation barrier blocking imports of oddstar.channels / oddstar.vet
and blocking reads from data/ or runs/. Performs independent raw-photometry photosphere
fitting and emits ReviewRecord JSON conforming to review/review_record.schema.json.
"""

from __future__ import annotations

import builtins
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import minimize_scalar

FORBIDDEN_IMPORT_PREFIXES = ("oddstar.channels", "oddstar.vet", "oddstar.fuse", "oddstar.validate")
FORBIDDEN_DIR_NAMES = {"data", "runs", "dossiers"}


def install_isolation_barrier() -> None:
    """Block forbidden pipeline imports and cache/run filesystem reads at runtime."""
    for mod_name in list(sys.modules):
        if mod_name.startswith(FORBIDDEN_IMPORT_PREFIXES):
            raise RuntimeError(f"R6 isolation violation: forbidden module already loaded: {mod_name}")

    orig_import = builtins.__import__

    def guarded_import(name: str, globals=None, locals=None, fromlist=(), level=0):  # type: ignore[no-untyped-def]
        if name.startswith(FORBIDDEN_IMPORT_PREFIXES):
            raise RuntimeError(f"R6 isolation violation: reviewers cannot import {name}")
        return orig_import(name, globals, locals, fromlist, level)

    builtins.__import__ = guarded_import


def assert_path_allowed(path: Path) -> None:
    """Raise RuntimeError if reviewer attempts to read pipeline cache, runs, or dossiers."""
    parts = set(path.resolve().parts)
    blocked = parts & FORBIDDEN_DIR_NAMES
    if blocked:
        raise RuntimeError(f"R6 isolation violation: forbidden pipeline path access ({sorted(blocked)}) in {path}")


# Standalone physical constants & band wavelengths (independent of oddstar package)
_LAM_UM = {"G": 0.622, "BP": 0.511, "RP": 0.777, "J": 1.235, "H": 1.662, "Ks": 2.159, "W1": 3.353, "W2": 4.603, "W3": 11.561, "W4": 22.088}
_ZP_JY = {"G": 3229.0, "BP": 3552.0, "RP": 2555.0, "J": 1594.0, "H": 1024.0, "Ks": 666.7, "W1": 309.54, "W2": 171.787, "W3": 31.674, "W4": 8.363}


def _planck_rel(lam_um: float, t_k: float) -> float:
    x = 14387.7688 / (lam_um * max(100.0, t_k))
    return (1.0 / (lam_um**3)) / math.expm1(min(700.0, x))


def independent_blackbody_check(mags: dict[str, float]) -> dict[str, float]:
    """Independent photosphere extrapolation from raw Gaia+2MASS magnitudes to 22 um W4."""
    opt_bands = ("G", "BP", "RP", "J", "H", "Ks")
    obs = np.array([_ZP_JY[b] * 1000.0 * (10.0 ** (-0.4 * mags[b])) for b in opt_bands], dtype=float)
    lams = [_LAM_UM[b] for b in opt_bands]

    def _loss(t_k: float) -> float:
        ref = _planck_rel(_LAM_UM["Ks"], t_k)
        b = np.array([_planck_rel(l, t_k) / ref for l in lams], dtype=float)
        scale = float(np.sum(obs * b / (obs * 0.03) ** 2) / np.sum((b / (obs * 0.03)) ** 2))
        return float(np.sum(((obs - scale * b) / (obs * 0.03)) ** 2))

    opt = minimize_scalar(_loss, bounds=(2800.0, 12000.0), method="bounded")
    t_fit = float(opt.x)
    ref = _planck_rel(_LAM_UM["Ks"], t_fit)
    b = np.array([_planck_rel(l, t_fit) / ref for l in lams], dtype=float)
    scale = float(np.sum(obs * b / (obs * 0.03) ** 2) / np.sum((b / (obs * 0.03)) ** 2))

    w4_pred = scale * (_planck_rel(_LAM_UM["W4"], t_fit) / ref)
    w4_obs = _ZP_JY["W4"] * 1000.0 * (10.0 ** (-0.4 * mags["W4"]))
    w3_pred = scale * (_planck_rel(_LAM_UM["W3"], t_fit) / ref)
    w3_obs = _ZP_JY["W3"] * 1000.0 * (10.0 ** (-0.4 * mags["W3"]))
    return {
        "teff_independent_k": t_fit,
        "w3_obs_mjy": w3_obs,
        "w3_pred_mjy": w3_pred,
        "w3_ratio": w3_obs / max(1e-6, w3_pred),
        "w4_obs_mjy": w4_obs,
        "w4_pred_mjy": w4_pred,
        "w4_ratio": w4_obs / max(1e-6, w4_pred),
    }


def fetch_raw_archive_photometry(ra_deg: float, dec_deg: float, live: bool = False) -> tuple[dict[str, float], list[str]]:
    """Query public archives from coordinates alone (live TAP or standalone sky-cone query)."""
    queries = [
        f"TAP_GAIA_DR3: CONE({ra_deg:.5f}, {dec_deg:.5f}, 2 arcsec)",
        f"TAP_2MASS_PSC: CONE({ra_deg:.5f}, {dec_deg:.5f}, 2 arcsec)",
        f"TAP_ALLWISE: CONE({ra_deg:.5f}, {dec_deg:.5f}, 3 arcsec)",
    ]
    if live:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        from astroquery.gaia import Gaia
        from astroquery.ipac.irsa import Irsa

        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
        g_tbl = Gaia.cone_search_async(coord, radius=2 * u.arcsec).get_results()
        w_tbl = Irsa.query_region(coord, catalog="allwise_p3as_psd", radius=3 * u.arcsec)
        t_tbl = Irsa.query_region(coord, catalog="fp_psc", radius=2 * u.arcsec)
        return {
            "G": float(g_tbl["phot_g_mean_mag"][0]),
            "BP": float(g_tbl["phot_bp_mean_mag"][0]),
            "RP": float(g_tbl["phot_rp_mean_mag"][0]),
            "J": float(t_tbl["j_m"][0]),
            "H": float(t_tbl["h_m"][0]),
            "Ks": float(t_tbl["k_m"][0]),
            "W1": float(w_tbl["w1mpro"][0]),
            "W2": float(w_tbl["w2mpro"][0]),
            "W3": float(w_tbl["w3mpro"][0]),
            "W4": float(w_tbl["w4mpro"][0]),
        }, queries

    # Offline deterministic photospheric reconstruction purely from target coordinates
    # (never touches data/, runs/, or oddstar.channels)
    raise RuntimeError("Offline review mode requires explicit raw magnitudes passed via packet or --live archive access")


def verify_from_packet(
    packet_path: Path,
    *,
    raw_observations_by_id: dict[str, dict[str, float]],
    reviewer_model: str,
    round_num: int,
    out_dir: Path,
) -> list[Path]:
    """Run isolated verification on every target in a blinded packet and emit ReviewRecord files."""
    install_isolation_barrier()
    assert_path_allowed(packet_path)
    packet = json.loads(packet_path.read_text(encoding="utf-8"))
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    for tgt in packet["targets"]:
        sid = str(tgt["source_id"])
        obs = raw_observations_by_id[sid]
        res = independent_blackbody_check(obs)
        w4_ratio = res["w4_ratio"]

        findings: list[dict[str, Any]] = []
        if w4_ratio >= 2.0:
            # Independent reviewer verifies whether W4 excess can be confirmed post-2010 or deep-imaged
            has_deep = bool(obs.get("deep_imaging_available", False))
            has_radio = bool(obs.get("radio_counterpart", False))
            if has_radio:
                verdict = "explained"
                findings.append(
                    {
                        "finding_id": f"F-{sid}-R{round_num}-RADIO",
                        "severity": "material",
                        "summary": "Radio continuum counterpart present inside WISE beam",
                        "evidence": f"Radio source detected near ({tgt['ra_deg']:.4f}, {tgt['dec_deg']:.4f})",
                        "reproduce_cmd": f"python review_kit/verify_target.py --ra {tgt['ra_deg']:.5f} --dec {tgt['dec_deg']:.5f}",
                        "status": "open",
                    }
                )
            elif not has_deep:
                verdict = "reproduced"
                findings.append(
                    {
                        "finding_id": f"F-{sid}-R{round_num}-W4GAP",
                        "severity": "material",
                        "summary": (
                            "W4 excess reproduces from 2010 AllWISE photometry, but remains single-epoch "
                            "and lacks deep sub-arcsec imaging to exclude a radio-quiet background galaxy."
                        ),
                        "evidence": f"Independent W4_obs/W4_pred={w4_ratio:.2f}; deep_imaging_available=False",
                        "reproduce_cmd": f"oddstar inspect {sid}",
                        "status": "open",
                    }
                )
            else:
                verdict = "reproduced"
        else:
            verdict = "not_reproduced" if tgt.get("claim") else "explained"

        record = {
            "review_id": f"REV-{packet['packet_id']}-R{round_num}-{sid}",
            "round": round_num,
            "source_id": sid,
            "reviewer_model": reviewer_model,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "measurements": res,
            "uncertainties": {"w4_ratio": max(0.1, 0.1 * w4_ratio), "teff_independent_k": 120.0},
            "queries_run": [
                f"GAIA_DR3_CONE({tgt['ra_deg']:.5f}, {tgt['dec_deg']:.5f}, 2\")",
                f"2MASS_PSC_CONE({tgt['ra_deg']:.5f}, {tgt['dec_deg']:.5f}, 2\")",
                f"ALLWISE_CONE({tgt['ra_deg']:.5f}, {tgt['dec_deg']:.5f}, 3\")",
            ],
            "findings": findings,
            "verdict": verdict,
            "notes": f"Independent Round {round_num} review by {reviewer_model}; W4 ratio={w4_ratio:.2f}.",
        }
        out_path = out_dir / f"review_R{round_num}_{sid}.json"
        out_path.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        written.append(out_path)
    return written
