"""Measure how wide a point source comes out of the pipeline's own fit to AllWISE Atlas cutouts.

`vision_rubric` calls a band "extended" when the fitted FWHM exceeds a reference by `vision_extended_fwhm_ratio`.
The reference in config 1.1.x is a theoretical one (native beam x sqrt(2), i.e. a Gaussian beam convolved with
itself). The first live run measured the reference star at 1.40 x that value in W3 and 1.31 x in W4, which is either
real extended emission or a reference that does not describe what this fitter returns for a point source in an
Atlas image (the real WISE beam is not Gaussian, and the fit is a Gaussian plus plane in a fixed box). Only a
measurement on known point sources can tell the two apart.

This tool selects ordinary bright stars from the AllWISE catalogue in fixed high-latitude fields, fetches their
Atlas cutouts with the pipeline's fetcher, runs the pipeline's `measure_wise_cutouts` on them unchanged, and writes
the distribution of fitted W3 and W4 FWHM to `validation/calibration/atlas_point_source_fwhm.json`.

Selection (fixed here before any calibration result was seen; no candidate is involved):
- fields: `CALIBRATION_FIELDS`, 1.0 deg radius each, all at |b| > 30 deg and none near a pipeline candidate;
- AllWISE: `ext_flg = 0`, `cc_flags = '0000'`, `nb = 1`, `na = 0` (single, unblended, no artifact flags),
  `w3sat = 0`, `w4sat = 0`, `w3snr >= 20`, `w4snr >= 10`, `w3rchi2 < 2`, `w4rchi2 < 2`, and a photospheric
  colour `-0.2 <= W3 - W4 <= 0.3` (excludes galaxies and dusty stars);
- the first `max_per_field` rows by designation in each field.

Rule for using the result (also fixed before the result was seen):
- the extendedness reference of a band becomes the median fitted FWHM of the calibration stars in that band;
- `vision_extended_fwhm_ratio` becomes max(1.3, 1 + 4 x robust sigma / median) over W3 and W4, rounded up to 0.05;
- if fewer than 15 stars per band are measured, nothing is changed and the check stays as it is.
The config change is made by hand in a separate commit with a CHANGELOG entry that cites the output file.
"""

from __future__ import annotations

import argparse
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig, compute_config_hash, load_config, repo_root
from oddstar.io.cache import RawCache

CALIBRATION_FIELDS: list[tuple[float, float]] = [
    (150.0, 2.2), (200.0, 30.0), (230.0, 10.0), (130.0, 45.0), (0.0, 20.0),
    (20.0, -30.0), (340.0, -15.0), (60.0, -40.0), (170.0, -20.0), (250.0, 40.0), (190.0, 60.0), (310.0, -50.0),
]
FIELD_RADIUS_DEG = 1.0
# IRSA took 4 minutes to answer this selection for one 1.5 deg field when it was first tried; a slow answer is not a failure
SELECTION_TIMEOUT_S = 600
MIN_STARS_PER_BAND = 15
OUTPUT_REL = "validation/calibration/atlas_point_source_fwhm.json"
BANDS = ("W3", "W4")
SELECT_COLUMNS = "designation, ra, dec, w3mpro, w3snr, w3rchi2, w4mpro, w4snr, w4rchi2"


def selection_query(ra: float, dec: float) -> str:
    return (
        f"SELECT {SELECT_COLUMNS} FROM allwise_p3as_psd "
        f"WHERE CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra:.4f}, {dec:.4f}, {FIELD_RADIUS_DEG})) = 1 "
        "AND ext_flg = 0 AND cc_flags = '0000' AND nb = 1 AND na = 0 AND w3sat = 0 AND w4sat = 0 "
        "AND w3snr >= 20 AND w4snr >= 10 AND w3rchi2 < 2 AND w4rchi2 < 2 "
        "AND w3mpro - w4mpro BETWEEN -0.2 AND 0.3"
    )


def band_stats(values: list[float]) -> dict[str, Any]:
    if not values:
        return {"n": 0}
    arr = np.asarray(values, dtype=float)
    med = float(np.median(arr))
    sig = float(1.4826 * np.median(np.abs(arr - med)))
    return {"n": int(arr.size), "median_arcsec": med, "robust_sigma_arcsec": sig, "min_arcsec": float(arr.min()),
            "max_arcsec": float(arr.max()), "p16_arcsec": float(np.percentile(arr, 16)), "p84_arcsec": float(np.percentile(arr, 84))}


def proposed_ratio(stats: dict[str, dict[str, Any]], current_ratio: float) -> float | None:
    """max(current, 1 + 4 sigma / median) over the bands, rounded up to 0.05; None when a band has too few stars."""
    worst = current_ratio
    for band in BANDS:
        st = stats.get(band) or {}
        if int(st.get("n", 0)) < MIN_STARS_PER_BAND:
            return None
        worst = max(worst, 1.0 + 4.0 * float(st["robust_sigma_arcsec"]) / float(st["median_arcsec"]))
    return math.ceil(worst * 20.0 - 1e-9) / 20.0


def run_calibration(cfg: PipelineConfig, cache: RawCache, *, max_per_field: int = 6, sleep: Any = None) -> dict[str, Any]:
    """Select the calibration stars, measure them with the pipeline's own cutout code, and return the record."""
    import time

    from oddstar.io import live_survey as ls
    from oddstar.io.cutouts import fetch_allwise_cutouts, measure_wise_cutouts

    do_sleep = sleep if sleep is not None else time.sleep
    vcfg = cfg.vetting
    size = float(vcfg["vision_cutout_size_arcsec"])
    min_snr = float(vcfg["vision_min_fit_snr"])
    fields: list[dict[str, Any]] = []
    stars: list[dict[str, Any]] = []
    for f_ra, f_dec in CALIBRATION_FIELDS:
        query = selection_query(f_ra, f_dec)
        field: dict[str, Any] = {"ra_deg": f_ra, "dec_deg": f_dec, "radius_deg": FIELD_RADIUS_DEG, "query": query}
        try:
            rows, _, _ = ls._tap_query_cached(url=ls.IRSA_TAP_URL, query=query, archive="allwise", release="atlas_psf_calibration",
                                              cache=cache, timeout_s=SELECTION_TIMEOUT_S)
        except Exception as exc:  # one failed field must not hide the others; it is recorded
            field["status"] = ls._failed(exc)
            fields.append(field)
            continue
        rows = sorted(rows, key=lambda r: str(r["designation"]))
        field.update({"status": "ok", "n_selected": len(rows), "n_used": min(len(rows), max_per_field)})
        fields.append(field)
        for row in rows[:max_per_field]:
            ra, dec = float(row["ra"]), float(row["dec"])
            star: dict[str, Any] = {"designation": str(row["designation"]).strip(), "ra_deg": ra, "dec_deg": dec,
                                    "catalogue": {k: (None if row.get(k) is None else float(row[k]))
                                                  for k in ("w3mpro", "w3snr", "w3rchi2", "w4mpro", "w4snr", "w4rchi2")},
                                    "bands": {}}
            try:
                images, _, log = fetch_allwise_cutouts(ra=ra, dec=dec, size_arcsec=size, cache=cache, sleep=do_sleep)
                meas = measure_wise_cutouts(images, ra_pred=ra, dec_pred=dec, vetting=vcfg)
                star["status"] = "ok"
                for band in BANDS:
                    m = meas.get(band) or {"status": "no_cutout"}
                    star["bands"][band] = {k: m.get(k) for k in ("status", "fwhm_arcsec", "fwhm_ref_arcsec", "fwhm_ratio", "flux_snr",
                                                                 "offset_arcsec", "brighter_secondary_peak", "sha256")}
            except Exception as exc:
                star["status"] = ls._failed(exc)
            stars.append(star)
    stats: dict[str, dict[str, Any]] = {}
    for band in BANDS:
        vals = [float(s["bands"][band]["fwhm_arcsec"]) for s in stars
                if s.get("status") == "ok" and s["bands"].get(band, {}).get("status") == "ok"
                and s["bands"][band].get("fwhm_arcsec") is not None and float(s["bands"][band].get("flux_snr") or 0.0) >= min_snr]
        st = band_stats(vals)
        if vals:
            ref = float(vcfg[f"{band.lower()}_beam_arcsec"]) * float(vcfg["wise_atlas_psf_broadening"])
            st["theoretical_reference_arcsec"] = ref
            st["median_over_theoretical_reference"] = st["median_arcsec"] / ref
        stats[band] = st
    current = float(vcfg["vision_extended_fwhm_ratio"])
    return {
        "purpose": "Fitted FWHM of known point sources in AllWISE Atlas cutouts, measured with oddstar.io.cutouts.measure_wise_cutouts",
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "config_hash": compute_config_hash(cfg.model_dump())[:12],
        "max_per_field": max_per_field,
        "min_stars_per_band": MIN_STARS_PER_BAND,
        "min_fit_snr": min_snr,
        "fields": fields,
        "n_stars": len(stars),
        "stats": stats,
        "current_extended_fwhm_ratio": current,
        "proposed_extended_fwhm_ratio": proposed_ratio(stats, current),
        "stars": stars,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--max-per-field", type=int, default=6)
    ap.add_argument("--out", default=None, help=f"Output path (default {OUTPUT_REL})")
    args = ap.parse_args(argv)
    cfg, _, _ = load_config()
    root = repo_root()
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    record = run_calibration(cfg, cache, max_per_field=args.max_per_field)
    out = Path(args.out) if args.out else root / OUTPUT_REL
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"out": str(out), "n_stars": record["n_stars"], "stats": record["stats"],
                      "proposed_extended_fwhm_ratio": record["proposed_extended_fwhm_ratio"],
                      "failed_fields": [f for f in record["fields"] if f["status"] != "ok"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
