"""Light-curve and IR time-domain vetting checks (Section 6).

Every check here either evaluates fetched data or returns `not_evaluated` with a note naming what was missing or
which query failed, and every check can `fail`. Synthetic fixtures carry explicit legacy flags
(`neighbour_eb_centroid_shift`, `tess_scattered_light_window`, `scan_systematic`, `w1_w2_neighbour_blend`,
`sso_crossing_t0`); a flag that is present counts as an evaluation, an absent flag never does.

Live inputs (all vetting-only, never channel inputs):
- `lightcurves.tess_sectors`: per-sector 30-min QLP raw SAP + TESScut aperture flux on matched cadences, with
  native cadences per bin, TESScut annulus background (`bkg_frac`) and flux-weighted centroid offsets.
- `lightcurves.tess_neighbours`: Gaia DR3 sources within `vetting.tess_neighbour_radius_px` TESS pixels.
- `lightcurves.sso_queries`: SkyBoT results at the flagged TESS event times.
- `lightcurves.optical_surveys`: ZTF / ASAS-SN seasonal and nightly medians.
- `neowise`: visit fluxes, per-visit PM-track residual vectors and the field-star `comparison` ensemble.

SIMBAD (`known_variable`): read through `oddstar.vet.common.simbad_lookup`; a failed or missing lookup cannot support a
pass, "no SIMBAD entry" is a real negative answer.

Detrend edge zone: within `detrend_edge_zone_window_frac` x `detrend_window_days` of a segment edge or data gap the
matched running median can erase or distort a dip. `analyze_tess_sector` screens that zone against the raw interior
flux, and the detrend-based checks (pipeline_agreement, systematics_window, neighbour_signal, sso_crossing,
second_survey) report `not_evaluated` instead of `pass` while a dip candidate lies there.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result, simbad_lookup, simbad_note

BTJD_TO_MJD = 2457000.0 - 2400000.5


# ----------------------------------------------------------------------------------------- helpers


def _od(cfg: PipelineConfig, key: str) -> float:
    return float(cfg.channels["optical_dips"][key])


def _it(cfg: PipelineConfig, key: str) -> float:
    return float(cfg.channels["ir_timedomain"][key])


def _vt(cfg: PipelineConfig, key: str) -> float:
    return float(cfg.vetting[key])


def _robust_sigma(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def _pearson(x: np.ndarray, y: np.ndarray) -> float | None:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3 or np.std(x) <= 0 or np.std(y) <= 0:
        return None
    return float(np.corrcoef(x, y)[0, 1])


def _r(v: float | None, nd: int = 4) -> float | None:
    if v is None or not math.isfinite(v):
        return None
    return round(float(v), nd)


def _combine(parts: list[str]) -> str:
    """Tri-state combination: any fail -> fail; else any pass -> pass; else not_evaluated."""
    if "fail" in parts:
        return "fail"
    if "pass" in parts:
        return "pass"
    return "not_evaluated"


def _tess_log(lc: dict[str, Any]) -> str:
    bits = []
    if lc.get("tess_status"):
        bits.append(f"status={lc['tess_status']}")
    log = lc.get("tess_query_log") or []
    if log:
        bits.append("; ".join(map(str, log))[:300])
    return f" ({'; '.join(bits)})" if bits else ""


def matched_detrend(
    t: np.ndarray, f: np.ndarray, window_days: float, gap_days: float = 0.5, mask_sigma: float = 2.5
) -> np.ndarray:
    """Detrend one light curve with the filter applied identically to every extraction being compared.

    Per continuous segment (split at gaps > `gap_days`), divide by a running median of `window_days`, computed
    after masking points more than `mask_sigma` robust sigma below a first-pass baseline so a real dip does not pull
    the baseline down (pipeline callers pass `vetting.lc_detrend_low_mask_sigma`; the keyword defaults equal the
    configured values and only serve direct callers). The window is the optical_dips channel's own `detrend_window_days`, so the comparison
    preserves the same dip timescales the detection channel searches.

    Within half a window of a segment edge the `mode="nearest"` padding repeats the edge cadences, so a dip there can be
    erased or made shallower; `_edge_zone_screen` covers that zone and the checks report it as `not_evaluated`.
    """
    from scipy.ndimage import median_filter

    t = np.asarray(t, dtype=float)
    f = np.asarray(f, dtype=float)
    out = np.ones_like(f)
    if len(t) == 0:
        return out
    breaks = np.where(np.diff(t) > gap_days)[0] + 1
    for seg in np.split(np.arange(len(t)), breaks):
        if len(seg) < 5:
            out[seg] = f[seg] / max(1e-9, float(np.median(f[seg])))
            continue
        ts, fs = t[seg], f[seg]
        dt = max(1e-4, float(np.median(np.diff(ts))))
        win = max(5, int(round(window_days / dt)) | 1)
        win = min(win, len(seg) if len(seg) % 2 == 1 else len(seg) - 1)
        base = median_filter(fs, size=win, mode="nearest")
        resid = fs - base
        noise = max(1e-9, 1.4826 * float(np.median(np.abs(resid))))
        valid = resid >= -mask_sigma * noise
        if int(np.sum(valid)) >= 5:
            base = median_filter(np.interp(ts, ts[valid], fs[valid]), size=win, mode="nearest")
        out[seg] = fs / np.maximum(1e-9, base)
    return out


def _bin_mean(t: np.ndarray, f: np.ndarray, width_days: float) -> np.ndarray:
    idx = np.floor((t - t[0]) / width_days).astype(int)
    _, inv = np.unique(idx, return_inverse=True)
    return np.bincount(inv, weights=f) / np.bincount(inv)


def _segments(t: np.ndarray, gap_days: float) -> list[tuple[float, float]]:
    if len(t) == 0:
        return []
    breaks = np.where(np.diff(t) > gap_days)[0] + 1
    return [(float(t[s[0]]), float(t[s[-1]])) for s in np.split(np.arange(len(t)), breaks) if len(s)]


def _edge_distance(tc: float, segs: list[tuple[float, float]]) -> float:
    """Days from `tc` to the nearest edge of the continuous segment holding it (0 when in a gap)."""
    for lo, hi in segs:
        if lo <= tc <= hi:
            return min(tc - lo, hi - tc)
    return 0.0


def _edge_zone_days(cfg: PipelineConfig) -> float:
    """Width of the detrend edge zone: `detrend_edge_zone_window_frac` x the matched running-median window."""
    return _od(cfg, "detrend_edge_zone_window_frac") * _od(cfg, "detrend_window_days")


def _edge_zone_screen(
    t: np.ndarray, fq: np.ndarray, fc: np.ndarray, w: np.ndarray, inv: np.ndarray, kept: np.ndarray,
    dq: np.ndarray, dc: np.ndarray, b_mid: np.ndarray, segs: list[tuple[float, float]], cfg: PipelineConfig,
) -> dict[str, Any]:
    """Dip candidates inside the detrend edge zone, measured without the running median.

    Within `_edge_zone_days` of a segment edge or data gap the running median (`mode="nearest"` padding) is
    dominated by the edge cadences themselves, so a dip there can be erased or made shallower. Each edge-zone
    cadence is therefore also compared with the median raw flux of the adjacent interior stretch (one detrend window
    just inside the zone; the whole segment outside the zone, or the whole segment, when that stretch is too short).
    A kept bin inside the zone is a candidate when either series, detrended or against the interior reference,
    reaches `min_dip_depth_frac`. Detrend-based checks report `not_evaluated` instead of `pass` when one exists.
    """
    zone = _edge_zone_days(cfg)
    window = _od(cfg, "detrend_window_days")
    gap = _od(cfg, "tess_segment_gap_days")
    thr = _od(cfg, "min_dip_depth_frac")
    eq = np.full(len(t), np.nan)
    ec = np.full(len(t), np.nan)
    breaks = np.where(np.diff(t) > gap)[0] + 1
    for seg in np.split(np.arange(len(t)), breaks):
        if len(seg) == 0:
            continue
        ts = t[seg]
        for dist in (ts - ts[0], ts[-1] - ts):
            in_zone = dist < zone
            if not in_zone.any():
                continue
            ref = (dist >= zone) & (dist < zone + window)
            if int(ref.sum()) < 5:
                ref = ~in_zone if int((~in_zone).sum()) >= 5 else np.ones(len(seg), dtype=bool)
            for f, e in ((fq, eq), (fc, ec)):
                level = max(1e-9, float(np.median(f[seg][ref])))
                idx = seg[in_zone]
                e[idx] = np.fmax(e[idx], 1.0 - f[idx] / level)
    cands: list[dict[str, Any]] = []
    n_zone = 0
    max_ref = 0.0
    for i in kept:
        dist = _edge_distance(float(b_mid[i]), segs)
        if dist >= zone:
            continue
        n_zone += 1
        sel = inv == i
        ref_depths = []
        for e in (eq, ec):
            ok = sel & np.isfinite(e)
            ref_depths.append(float(np.sum(w[ok] * e[ok]) / np.sum(w[ok])) if ok.any() and np.sum(w[ok]) > 0 else float("nan"))
        rq, rc = ref_depths
        finite = [d for d in ref_depths if math.isfinite(d)]
        if finite:
            max_ref = max(max_ref, *finite)
        if max([float(dq[i]), float(dc[i]), *finite]) >= thr:
            cands.append({
                "btjd": round(float(b_mid[i]), 4), "edge_distance_days": round(dist, 3),
                "depth_qlp_detrended": round(float(dq[i]), 4), "depth_tesscut_detrended": round(float(dc[i]), 4),
                "depth_qlp_vs_interior": _r(rq), "depth_tesscut_vs_interior": _r(rc),
            })
    return {"zone_days": round(zone, 3), "n_edge_zone_bins": n_zone, "max_depth_vs_interior": round(max_ref, 4), "candidates": cands}


def _edge_note(cands: list[dict[str, Any]], cfg: PipelineConfig) -> str:
    c0 = max(cands, key=lambda c: max(c["depth_qlp_detrended"], c["depth_tesscut_detrended"], c["depth_qlp_vs_interior"] or 0.0,
                                      c["depth_tesscut_vs_interior"] or 0.0))
    return (f"{len(cands)} dip candidate bin(s) lie within {_edge_zone_days(cfg):.2f} d of a segment edge or data gap "
            f"(detrend edge zone = {_od(cfg, 'detrend_edge_zone_window_frac'):.2f} x the {_od(cfg, 'detrend_window_days'):.1f} d window), e.g. "
            f"{'sector ' + str(c0['sector']) + ' ' if 'sector' in c0 else ''}BTJD {c0['btjd']:.3f}, {c0['edge_distance_days']:.2f} d from the edge "
            f"(QLP/TESScut depth vs interior {c0['depth_qlp_vs_interior']}/{c0['depth_tesscut_vs_interior']}, detrended "
            f"{c0['depth_qlp_detrended']:.3f}/{c0['depth_tesscut_detrended']:.3f}); the matched running median can erase or distort a dip "
            f"there, so the detrend-based test is not evaluated for it")


def _sector_edge_candidates(a: dict[str, Any]) -> list[dict[str, Any]]:
    return [dict(c, sector=a["summary"]["sector"]) for c in a["summary"]["detrend_edge_zone"]["candidates"]]


def tess_edge_zone_candidates(lc: dict[str, Any], cfg: PipelineConfig, *, include_legacy: bool = False) -> list[dict[str, Any]]:
    """Detrend edge-zone dip candidates of every usable TESS sector (see `_edge_zone_screen`)."""
    secs, _ = _analysis_sectors(lc, include_legacy=include_legacy)
    out: list[dict[str, Any]] = []
    for sec in secs:
        a = analyze_tess_sector(sec, cfg)
        if a and a["usable"]:
            out.extend(_sector_edge_candidates(a))
    return out


def _depths_agree(d_primary: float, d_secondary: float, cfg: PipelineConfig) -> bool:
    """Depth agreement of one bin: both below the quiet level, or within the absolute tolerance and ratio."""
    quiet = _od(cfg, "pipeline_quiet_depth_frac")
    tol = _od(cfg, "pipeline_depth_abs_tol_frac")
    ratio = _od(cfg, "pipeline_depth_min_ratio")
    return (d_primary < quiet and d_secondary < quiet) or (abs(d_primary - d_secondary) <= tol and d_secondary >= ratio * d_primary)


def _analysis_sectors(lc: dict[str, Any], *, include_legacy: bool = True) -> tuple[list[dict[str, Any]], str]:
    secs = [s for s in lc.get("tess_sectors", []) or [] if len(s.get("flux_qlp_raw", [])) >= 20]
    if secs:
        return secs, "tess_sectors"
    if include_legacy:
        t = list(lc.get("tess_time_btjd", []) or [])
        f1 = list(lc.get("tess_flux_qlp", []) or [])
        f2 = list(lc.get("tess_flux_tesscut", []) or [])
        if len(t) >= 20 and len(t) == len(f1) == len(f2):
            return [{"sector": 0, "time_btjd": t, "flux_qlp_raw": f1, "flux_tesscut": f2, "legacy": True}], "legacy_arrays"
    return [], "none"


def analyze_tess_sector(sec: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any] | None:
    """Matched-detrended, time-coincident 6 h bin comparison of the QLP and TESScut series of one sector.

    Bins holding fewer than `pipeline_min_bin_fill_frac` of the expected native cadences are dropped. For the
    deepest QLP bin and the deepest TESScut bin, the other series is read in that same bin. The coincident depth
    of a bin is min(QLP depth, TESScut depth); its maximum is the sector's deepest coincident feature. Keys
    starting with `_` hold arrays for the other checks and are never written to evidence.
    """
    t = np.asarray(sec.get("time_btjd", []), dtype=float)
    fq = np.asarray(sec.get("flux_qlp_raw", []), dtype=float)
    fc = np.asarray(sec.get("flux_tesscut", []), dtype=float)
    sector = int(sec.get("sector", -1))
    if len(t) < 20 or len(fq) != len(t) or len(fc) != len(t):
        return None
    window = _od(cfg, "detrend_window_days")
    gap = _od(cfg, "tess_segment_gap_days")
    bin_days = _od(cfg, "min_duration_hours") / 24.0
    thr = _od(cfg, "min_dip_depth_frac")
    mask_sigma = _vt(cfg, "lc_detrend_low_mask_sigma")
    rq = matched_detrend(t, fq, window, gap, mask_sigma)
    rc = matched_detrend(t, fc, window, gap, mask_sigma)
    n_per = np.asarray(sec.get("n_per_bin") or [], dtype=float)
    if len(n_per) != len(t):
        n_per = np.ones(len(t))
        cadence = float(np.median(np.diff(t)))
    else:
        cadence = float(sec.get("cadence_days") or np.median(np.diff(t)))
    idx = np.floor((t - t[0]) / bin_days).astype(int)
    uniq, inv = np.unique(idx, return_inverse=True)
    wsum = np.bincount(inv, weights=n_per)
    cnt = wsum  # native cadences per bin
    bq = np.bincount(inv, weights=n_per * rq) / wsum
    bc = np.bincount(inv, weights=n_per * rc) / wsum
    b_mid = np.bincount(inv, weights=n_per * t) / wsum
    b_lo = t[0] + uniq * bin_days
    b_hi = b_lo + bin_days
    expected = bin_days / max(1e-6, cadence)
    keep = cnt >= _od(cfg, "pipeline_min_bin_fill_frac") * expected
    kept = np.where(keep)[0]
    base = {"sector": sector, "n_bins_total": int(len(uniq)), "n_bins_kept": int(len(kept)),
            "n_bins_dropped_underfilled": int(len(uniq) - len(kept)), "expected_cadences_per_bin": round(expected, 1)}
    if len(kept) < 4:
        return {"usable": False, "summary": dict(base, reason="fewer than 4 sufficiently filled bins")}
    dq = 1.0 - bq
    dc = 1.0 - bc
    co = np.minimum(dq, dc)
    iq = int(kept[np.argmax(dq[kept])])
    ic = int(kept[np.argmax(dc[kept])])
    ico = int(kept[np.argmax(co[kept])])
    segs = _segments(t, gap)

    def _rec(i: int) -> dict[str, Any]:
        return {
            "btjd": round(float(b_mid[i]), 4),
            "btjd_start": round(float(b_lo[i]), 4),
            "btjd_end": round(float(b_hi[i]), 4),
            "depth_qlp": round(float(dq[i]), 4),
            "depth_tesscut": round(float(dc[i]), 4),
            "depth_coincident": round(float(co[i]), 4),
            "n_cadences": int(round(float(cnt[i]))),
            "edge_distance_days": round(_edge_distance(float(b_mid[i]), segs), 3),
            "in_detrend_edge_zone": bool(_edge_distance(float(b_mid[i]), segs) < _edge_zone_days(cfg)),
            "qlp_exceeds_min_depth": bool(dq[i] >= thr),
            "tesscut_exceeds_min_depth": bool(dc[i] >= thr),
            "coincident_exceeds_min_depth": bool(co[i] >= thr),
        }

    r = _pearson(bq[kept], bc[kept])
    quiet = _od(cfg, "pipeline_quiet_depth_frac")
    feature = bool(dq[iq] >= quiet or dc[ic] >= quiet)
    agree_q = _depths_agree(float(dq[iq]), float(dc[iq]), cfg)
    agree_c = _depths_agree(float(dc[ic]), float(dq[ic]), cfg)
    min_r = _od(cfg, "pipeline_min_binned_corr")
    corr_ok = (not feature) or (r is not None and r >= min_r)
    reasons = []
    if not agree_q:
        reasons.append(f"deepest QLP bin BTJD {b_mid[iq]:.3f} (QLP {dq[iq]:.3f}) not matched by TESScut ({dc[iq]:.3f})")
    if not agree_c:
        reasons.append(f"deepest TESScut bin BTJD {b_mid[ic]:.3f} (TESScut {dc[ic]:.3f}) not matched by QLP ({dq[ic]:.3f})")
    if not corr_ok:
        reasons.append(f"binned r={'n/a' if r is None else f'{r:.2f}'} < {min_r:.2f} with a feature above {quiet:.2f}")
    thr_bins = [int(i) for i in kept if dq[i] >= thr or dc[i] >= thr]
    event_bins = [int(i) for i in kept if co[i] >= thr]
    summary = dict(
        base,
        deepest_qlp=_rec(iq),
        deepest_tesscut=_rec(ic),
        deepest_coincident=_rec(ico),
        qlp_robust_sigma=round(_robust_sigma(bq[kept]), 5),
        tesscut_robust_sigma=round(_robust_sigma(bc[kept]), 5),
        binned_pearson_r=_r(r, 3),
        correlation_used=feature,
        agree=bool(agree_q and agree_c and corr_ok),
        disagreement=reasons,
        n_bins_above_min_depth=len(thr_bins),
        legacy_arrays=bool(sec.get("legacy", False)),
        detrend_edge_zone=_edge_zone_screen(t, fq, fc, n_per, inv, kept, dq, dc, b_mid, segs, cfg),
    )
    return {
        "usable": True, "summary": summary, "_t": t, "_rq": rq, "_rc": rc, "_inv": inv, "_b_lo": b_lo, "_b_hi": b_hi,
        "_b_mid": b_mid, "_dq": dq, "_dc": dc, "_co": co, "_ico": ico, "_thr_bins": thr_bins, "_event_bins": event_bins,
        "_segs": segs,
    }


def tess_dip_events(lc: dict[str, Any], cfg: PipelineConfig) -> tuple[list[dict[str, Any]], int]:
    """Flagged TESS events: per real sector, the deepest coincident 6 h bin when it exceeds `min_dip_depth_frac`.

    Returns (events, number of usable sectors). Legacy fixture arrays are not used (their times are synthetic).
    """
    secs, _ = _analysis_sectors(lc, include_legacy=False)
    events: list[dict[str, Any]] = []
    n_usable = 0
    for sec in secs:
        a = analyze_tess_sector(sec, cfg)
        if not a or not a["usable"]:
            continue
        n_usable += 1
        dc = a["summary"]["deepest_coincident"]
        if dc["coincident_exceeds_min_depth"]:
            events.append(dict(dc, sector=a["summary"]["sector"], depth=dc["depth_coincident"]))
    return events, n_usable


# ------------------------------------------------------------------------------- pipeline_agreement


def check_pipeline_agreement(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    secs, source = _analysis_sectors(lc)
    window = _od(cfg, "detrend_window_days")
    bin_days = _od(cfg, "min_duration_hours") / 24.0
    thr = _od(cfg, "min_dip_depth_frac")
    if not secs:
        return make_vet_result(out_dir, check="pipeline_agreement", status="not_evaluated",
                               note="Independent second light-curve extraction missing" + _tess_log(lc))
    analyses = [analyze_tess_sector(s, cfg) for s in secs]
    usable = [a for a in analyses if a and a["usable"]]
    skipped = [a["summary"] for a in analyses if a and not a["usable"]]
    if not usable:
        return make_vet_result(out_dir, check="pipeline_agreement", status="not_evaluated",
                               note=f"No TESS sector has >= 4 sufficiently filled {bin_days * 24:.0f} h bins" + _tess_log(lc),
                               metrics={"sectors_skipped": skipped})
    per_sector = [a["summary"] for a in usable]
    bad = [p["sector"] for p in per_sector if not p["agree"]]
    best = max(per_sector, key=lambda p: p["deepest_coincident"]["depth_coincident"])
    bc = best["deepest_coincident"]
    src_txt = "legacy fixture arrays" if source == "legacy_arrays" else f"sector {best['sector']}"
    d1 = max(p["deepest_qlp"]["depth_qlp"] for p in per_sector)
    d2 = max(p["deepest_tesscut"]["depth_tesscut"] for p in per_sector)
    note = (
        f"Deepest coincident feature: {src_txt} at BTJD {bc['btjd']:.3f} (QLP depth {bc['depth_qlp']:.3f}, TESScut depth "
        f"{bc['depth_tesscut']:.3f}, {bc['n_cadences']} cadences, {bc['edge_distance_days']:.2f} d from segment edge; "
        f"{'exceeds' if bc['coincident_exceeds_min_depth'] else 'below'} min_dip_depth_frac={thr:.3f}). "
        f"Matched detrending ({window:.1f} d running median, {bin_days * 24:.0f} h bins with >= "
        f"{_od(cfg, 'pipeline_min_bin_fill_frac') * 100:.0f}% of expected cadences) over {len(per_sector)} sector(s), deepest bins "
        f"compared at the same time in both series; binned r required >= {_od(cfg, 'pipeline_min_binned_corr'):.2f} only when a "
        f"feature exceeds {_od(cfg, 'pipeline_quiet_depth_frac'):.2f} (sector r: "
        f"{', '.join(str(p['binned_pearson_r']) for p in per_sector)})"
    )
    if bad:
        why = "; ".join(f"s{p['sector']}: {', '.join(p['disagreement'])}" for p in per_sector if not p["agree"])
        note += f"; disagreement in sector(s) {bad}: {why}"
    edge_c = [c for a in usable for c in _sector_edge_candidates(a)]
    if edge_c and not bad:
        note += "; " + _edge_note(edge_c, cfg)
    note += (". Both extractions use the same TESS FFI pixels, so agreement validates the light-curve extraction, "
             "not the astrophysical origin of the signal")
    return make_vet_result(
        out_dir,
        check="pipeline_agreement",
        status="fail" if bad else ("not_evaluated" if edge_c else "pass"),
        note=note,
        metrics={
            "source": source,
            "qlp_max_depth": d1,
            "tesscut_max_depth": d2,
            "deepest_coincident": dict(bc, sector=best["sector"]),
            "tic_id": lc.get("tess_tic_id"),
            "detrend_window_days": window,
            "bin_hours": bin_days * 24.0,
            "sectors": per_sector,
            "sectors_skipped": skipped,
            "shared_pixels_caveat": True,
            "detrend_edge_zone_days": _edge_zone_days(cfg),
            "detrend_edge_zone_candidates": edge_c,
        },
    )


# ------------------------------------------------------------------------------- systematics_window


def _tess_aux(sec: dict[str, Any], n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray] | None:
    keys = ("bkg_frac", "centroid_dcol_px", "centroid_drow_px")
    if not all(len(sec.get(k) or []) == n for k in keys):
        return None
    return tuple(np.asarray(sec[k], dtype=float) for k in keys)  # type: ignore[return-value]


def _tess_systematics(lc: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    secs, _ = _analysis_sectors(lc, include_legacy=False)
    if not secs:
        return {"status": "not_evaluated", "note": "TESS: no sectors loaded" + _tess_log(lc)}
    thr = _od(cfg, "min_dip_depth_frac")
    margin = _od(cfg, "systematics_edge_margin_days")
    bsig = _od(cfg, "systematics_bkg_sigma")
    bfrac = _od(cfg, "systematics_bkg_depth_frac")
    jsig = _od(cfg, "systematics_centroid_jump_sigma")
    jmin = _od(cfg, "systematics_centroid_jump_min_px")
    bin_days = _od(cfg, "min_duration_hours") / 24.0
    table: list[dict[str, Any]] = []
    edge_c: list[dict[str, Any]] = []
    n_aux = 0
    for sec in secs:
        a = analyze_tess_sector(sec, cfg)
        if not a or not a["usable"]:
            continue
        t = a["_t"]
        aux = _tess_aux(sec, len(t))
        if aux is None:
            continue
        n_aux += 1
        edge_c.extend(_sector_edge_candidates(a))
        bkg, cx, cy = aux
        b_med, b_sig = float(np.median(bkg)), max(1e-9, _robust_sigma(bkg))
        dt = np.diff(t)
        steps = np.hypot(np.diff(cx), np.diff(cy))
        contiguous = dt <= _vt(cfg, "lc_tess_contiguous_dt_factor") * float(np.median(dt))
        s_med = float(np.median(steps[contiguous])) if contiguous.any() else 0.0
        s_sig = max(1e-6, _robust_sigma(steps[contiguous])) if contiguous.any() else 1e-6
        tested = sorted(set(a["_thr_bins"]) | {a["_ico"]})
        for i in tested:
            in_bin = a["_inv"] == i
            above = i in a["_thr_bins"]
            depth = float(max(a["_dq"][i], a["_dc"][i])) if above else float(a["_co"][i])
            ex = float(np.max(bkg[in_bin]) - b_med) if in_bin.any() else 0.0
            z = ex / b_sig
            bkg_flag = bool(z >= bsig and ex >= bfrac * max(depth, 0.0) and ex > 0)
            edge = _edge_distance(float(a["_b_mid"][i]), a["_segs"])
            edge_flag = bool(edge < margin)
            lo, hi = float(a["_b_lo"][i]) - bin_days, float(a["_b_hi"][i]) + bin_days
            win = contiguous & (t[:-1] >= lo) & (t[1:] <= hi)
            jump = float(np.max(steps[win])) if win.any() else 0.0
            jz = (jump - s_med) / s_sig
            jump_flag = bool(jz >= jsig and jump >= jmin)
            table.append({
                "sector": a["summary"]["sector"], "btjd": round(float(a["_b_mid"][i]), 4), "depth": round(depth, 4),
                "above_min_depth": bool(above or a["_co"][i] >= thr),
                "deepest_coincident": bool(i == a["_ico"]),
                "bkg_excess_frac": round(ex, 6), "bkg_excess_sigma": round(z, 2), "bkg_flag": bkg_flag,
                "edge_distance_days": round(edge, 3), "edge_flag": edge_flag,
                "centroid_jump_px": round(jump, 4), "centroid_jump_sigma": round(jz, 2), "jump_flag": jump_flag,
            })
    if n_aux == 0:
        return {"status": "not_evaluated", "note": "TESS: sectors lack per-cadence background/centroid series" + _tess_log(lc)}
    flagged = [e for e in table if e["above_min_depth"] and (e["bkg_flag"] or e["edge_flag"] or e["jump_flag"])]
    sub = [e for e in table if not e["above_min_depth"] and (e["bkg_flag"] or e["edge_flag"] or e["jump_flag"])]
    n_ev = sum(1 for e in table if e["above_min_depth"])
    if flagged:
        f0 = flagged[0]
        why = [w for w, k in (("background excursion", "bkg_flag"), (f"within {margin:.2f} d of a segment edge", "edge_flag"),
                              ("centroid jump", "jump_flag")) if f0[k]]
        note = (f"TESS: {len(flagged)} event bin(s) above {thr:.3f} coincide with systematics, e.g. sector {f0['sector']} BTJD "
                f"{f0['btjd']:.3f} (depth {f0['depth']:.3f}): {', '.join(why)}")
        status = "fail"
    else:
        note = (f"TESS: {n_ev} event bin(s) above {thr:.3f} in {n_aux} sector(s) with background/centroid series; none coincides "
                f"with a background excursion (>= {bsig:.0f} sigma and >= {bfrac:.2f} x depth), a segment edge (< {margin:.2f} d) or "
                f"a centroid jump (>= {jsig:.0f} sigma and >= {jmin:.2f} px)")
        if sub:
            note += f"; {len(sub)} sub-threshold deepest bin(s) coincide with a systematic (not a flagged event)"
        status = "pass"
        if edge_c:
            note += "; " + _edge_note(edge_c, cfg)
            status = "not_evaluated"
    return {"status": status, "note": note, "events": table, "n_sectors": n_aux, "detrend_edge_zone_candidates": edge_c,
            "edge_blocked": bool(edge_c) and status != "fail"}


def _ir_variation(neo: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    """W1/W2 visit variability using the ir_timedomain channel's own thresholds (chi2, slope, jump SNR)."""
    mjd = np.asarray(neo.get("mjd", []) or [], dtype=float)
    w1 = np.asarray(neo.get("w1_mjy", []) or [], dtype=float)
    w2 = np.asarray(neo.get("w2_mjy", []) or [], dtype=float)
    n = len(mjd)
    if n < 3 or len(w2) != n or len(w1) != n or np.any(w2 <= 0):
        return {"active": False, "n_visits": n, "usable": False}
    e1 = np.maximum(1e-3, np.asarray(neo.get("w1_err_mjy") or np.ones(n), dtype=float))
    e2 = np.maximum(1e-3, np.asarray(neo.get("w2_err_mjy") or np.ones(n), dtype=float))
    chi1 = float(np.mean(((w1 - np.median(w1)) / e1) ** 2))
    chi2 = float(np.mean(((w2 - np.median(w2)) / e2) ** 2))
    t_yr = (mjd - np.mean(mjd)) / 365.25
    denom = max(1e-6, float(np.sum(t_yr**2 / e2**2)))
    slope = float(np.sum(t_yr * (w2 - np.mean(w2)) / e2**2)) / denom
    slope_snr = abs(slope) * math.sqrt(denom)
    jumps = np.abs(np.diff(w2)) / np.sqrt(e2[1:] ** 2 + e2[:-1] ** 2)
    jump_snr = float(np.max(jumps)) if len(jumps) else 0.0
    active = (chi1 >= _it(cfg, "reduced_chi2_threshold") or chi2 >= _it(cfg, "reduced_chi2_threshold")
              or slope_snr >= _it(cfg, "slope_snr_threshold") or jump_snr >= _it(cfg, "jump_snr_threshold"))
    return {"active": bool(active) and n >= int(_it(cfg, "min_visits")), "usable": True, "n_visits": n,
            "chi2_w1": round(chi1, 2), "chi2_w2": round(chi2, 2), "slope_snr": round(slope_snr, 2), "jump_snr": round(jump_snr, 2)}


def _neowise_systematics(neo: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    comp = neo.get("comparison")
    if not neo.get("mjd"):
        return {"status": "not_evaluated", "note": "NEOWISE: no visit light curve"}
    if not comp:
        return {"status": "not_evaluated", "note": "NEOWISE: no comparison-star ensemble fetched"}
    if comp.get("status") != "ok":
        log = "; ".join(map(str, comp.get("query_log") or []))[:240]
        return {"status": "not_evaluated", "note": f"NEOWISE: comparison ensemble {comp.get('status')} "
                f"({comp.get('n_stars', 0)} stars){(' - ' + log) if log else ''}"}
    w2 = np.asarray(neo.get("w2_mjy", []), dtype=float)
    ens = comp.get("ensemble_dw2_mag") or []
    ens_err = comp.get("ensemble_dw2_err_mag") or [None] * len(ens)
    if len(ens) != len(w2) or np.any(w2 <= 0):
        return {"status": "not_evaluated", "note": "NEOWISE: comparison ensemble is not on the target visit grid"}
    tgt = -2.5 * np.log10(w2 / np.median(w2))
    ok = np.asarray([e is not None for e in ens])
    n_ok = int(ok.sum())
    if n_ok < int(_it(cfg, "min_visits")):
        return {"status": "not_evaluated", "note": f"NEOWISE: ensemble covers only {n_ok} visits"}
    e = np.asarray([v for v in ens if v is not None], dtype=float)
    ee = np.asarray([(x if x is not None else np.nan) for x, k in zip(ens_err, ok) if k], dtype=float)
    tt = tgt[ok]
    r = _pearson(tt, e)
    amp_t = float(np.ptp(tt))
    amp_e = float(np.ptp(e))
    if np.all(np.isfinite(ee)) and np.all(ee > 0):
        ens_chi2 = float(np.mean(((e - np.median(e)) / ee) ** 2))
    else:
        ens_chi2 = float("nan")
    ens_sig = bool(math.isfinite(ens_chi2) and ens_chi2 >= _it(cfg, "reduced_chi2_threshold"))
    var = _ir_variation(neo, cfg)
    metrics = {"n_comparison_stars": comp.get("n_stars"), "target_w2_amp_mag": round(amp_t, 4), "ensemble_w2_amp_mag": round(amp_e, 4),
               "ensemble_reduced_chi2": _r(ens_chi2, 2), "target_ensemble_pearson_r": _r(r, 3), "n_visits": n_ok,
               "target_variation": var}
    corr_thr = _it(cfg, "systematics_ensemble_corr")
    flag_corr = bool(r is not None and r >= corr_thr and amp_e >= _it(cfg, "systematics_ensemble_corr_min_amp_frac") * amp_t)
    flag_amp = bool(ens_sig and amp_e >= _it(cfg, "systematics_ensemble_amp_frac") * amp_t)
    base = (f"NEOWISE: {comp.get('n_stars')} comparison stars on the target visit grid, ensemble W2 amplitude {amp_e:.3f} mag "
            f"(reduced chi2 {metrics['ensemble_reduced_chi2']}) vs target {amp_t:.3f} mag, r={metrics['target_ensemble_pearson_r']}")
    if not var.get("active"):
        return {"status": "pass", "note": base + "; target W2 shows no significant variation to attribute", "metrics": metrics}
    if flag_corr or flag_amp:
        why = []
        if flag_corr:
            why.append(f"target variation correlates with the field ensemble (r >= {corr_thr:.2f})")
        if flag_amp:
            why.append(f"ensemble varies by >= {_it(cfg, 'systematics_ensemble_amp_frac'):.2f} x the target amplitude")
        return {"status": "fail", "note": base + "; scan/zero-point systematic: " + " and ".join(why), "metrics": metrics}
    return {"status": "pass", "note": base + "; target W2 variation is not shared by field stars in the same frames", "metrics": metrics}


def check_systematics_window(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    neo = star.get("neowise", {}) or {}
    statuses: list[str] = []
    notes: list[str] = []
    metrics: dict[str, Any] = {}
    legacy = {k: bool(d[k]) for d, k in ((lc, "tess_scattered_light_window"), (neo, "scan_systematic")) if k in d}
    if legacy:
        metrics["explicit_flags"] = legacy
        if any(legacy.values()):
            statuses.append("fail")
            notes.append("Event aligns with documented scattered-light or scan-direction systematic window (explicit flag)")
        else:
            statuses.append("pass")
            notes.append("Explicit flags: event outside scattered-light/momentum-dump/scan systematic windows")
    tess = _tess_systematics(lc, cfg)
    neow = _neowise_systematics(neo, cfg)
    for part, key in ((tess, "tess"), (neow, "neowise")):
        statuses.append(part["status"])
        notes.append(part["note"])
        metrics[key] = {k: v for k, v in part.items() if k != "note"}
    status = _combine(statuses)
    if status == "pass" and tess.get("edge_blocked"):
        status = "not_evaluated"
        notes.insert(0, "A TESS dip candidate inside the detrend edge zone is not evaluated, so the check cannot pass")
    elif status == "not_evaluated" and not tess.get("edge_blocked"):
        notes.insert(0, "Neither TESS background/centroid series nor a NEOWISE comparison ensemble available")
    return make_vet_result(out_dir, check="systematics_window", status=status, note=" | ".join(notes), metrics=metrics)


# --------------------------------------------------------------------------------- neighbour_signal


def _tess_neighbour(lc: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    secs, _ = _analysis_sectors(lc, include_legacy=False)
    if not secs:
        return {"status": "not_evaluated", "note": "TESS: no sectors loaded" + _tess_log(lc)}
    thr = _od(cfg, "min_dip_depth_frac")
    win_d = _od(cfg, "neighbour_centroid_window_days")
    s_thr = _od(cfg, "neighbour_source_offset_sigma")
    min_px = _od(cfg, "neighbour_source_offset_min_px")
    nb_info = lc.get("tess_neighbours") or {}
    sources = nb_info.get("sources") or []
    rows: list[dict[str, Any]] = []
    edge_c: list[dict[str, Any]] = []
    n_aux = 0
    for sec in secs:
        a = analyze_tess_sector(sec, cfg)
        if not a or not a["usable"]:
            continue
        t = a["_t"]
        aux = _tess_aux(sec, len(t))
        if aux is None:
            continue
        n_aux += 1
        edge_c.extend(_sector_edge_candidates(a))
        _, cx, cy = aux
        excl = set(a["_thr_bins"]) | set(a["_event_bins"])
        out_mask_base = ~np.isin(a["_inv"], list(excl))
        for i in a["_event_bins"]:
            in_m = a["_inv"] == i
            out_m = out_mask_base & (np.abs(t - a["_b_mid"][i]) <= win_d)
            n_in, n_out = int(in_m.sum()), int(out_m.sum())
            depth = float(1.0 - np.mean(a["_rc"][in_m])) if n_in else 0.0
            if n_in < 2 or n_out < 10 or depth <= 0:
                continue
            c_in = np.asarray([np.mean(cx[in_m]), np.mean(cy[in_m])])
            c_out = np.asarray([np.mean(cx[out_m]), np.mean(cy[out_m])])
            err = np.asarray([_robust_sigma(cx[out_m]), _robust_sigma(cy[out_m])]) * math.sqrt(1.0 / n_in + 1.0 / n_out)
            amp = (1.0 - depth) / depth
            src = c_out - amp * (c_in - c_out)  # difference-image (flux deficit) location
            src_err = np.maximum(1e-6, amp * err)
            offset = float(np.hypot(*src))
            off_sig = float(np.hypot(*(src / src_err)))
            nearest = None
            for nb in sec.get("neighbour_px") or []:
                d = math.hypot(src[0] - float(nb["dcol_px"]), src[1] - float(nb["drow_px"]))
                if nearest is None or d < nearest["dist_px"]:
                    nearest = {"id": nb.get("id"), "dist_px": round(d, 3), "flux_ratio": nb.get("flux_ratio")}
            rows.append({
                "sector": a["summary"]["sector"], "btjd": round(float(a["_b_mid"][i]), 4), "depth_tesscut": round(depth, 4),
                "centroid_shift_px": [round(float(v), 5) for v in (c_in - c_out)],
                "deficit_location_px": [round(float(v), 3) for v in src], "deficit_offset_px": round(offset, 3),
                "deficit_offset_sigma": round(off_sig, 2), "n_in": n_in, "n_out": n_out,
                "flag": bool(off_sig >= s_thr and offset >= min_px), "nearest_neighbour": nearest,
            })
    if n_aux == 0:
        return {"status": "not_evaluated", "note": "TESS: sectors lack flux-weighted centroid series" + _tess_log(lc)}
    nb_txt = ""
    if nb_info:
        if nb_info.get("status") == "ok":
            brightest = max((s.get("flux_ratio") or 0.0 for s in sources), default=0.0)
            nb_txt = (f"; Gaia DR3 within {nb_info.get('radius_arcsec')}\": {len(sources)} source(s), brightest flux ratio "
                      f"{brightest:.3f}")
        else:
            nb_txt = f"; Gaia neighbour cone {nb_info.get('status')}"
    edge_txt = ("; " + _edge_note(edge_c, cfg)) if edge_c else ""
    edge_m = {"detrend_edge_zone_candidates": edge_c, "edge_blocked": bool(edge_c)}
    pass_status = "not_evaluated" if edge_c else "pass"
    if not rows:
        return {"status": pass_status, "note": f"TESS: no dip event above {thr:.3f} in {n_aux} sector(s) with centroid series, so no dip to "
                f"attribute to a neighbour{nb_txt}{edge_txt}", "events": [], "n_sectors": n_aux, **edge_m}
    flagged = [r for r in rows if r["flag"]]
    if flagged:
        f0 = flagged[0]
        nn = f0["nearest_neighbour"]
        note = (f"TESS: flux deficit of sector {f0['sector']} BTJD {f0['btjd']:.3f} located {f0['deficit_offset_px']:.2f} px from the "
                f"target ({f0['deficit_offset_sigma']:.1f} sigma)" + (f", nearest Gaia neighbour {nn['id']} at {nn['dist_px']:.2f} px" if nn else "")
                + nb_txt)
        return {"status": "fail", "note": note, "events": rows, "n_sectors": n_aux, **dict(edge_m, edge_blocked=False)}
    worst = max(rows, key=lambda r: r["deficit_offset_sigma"])
    return {"status": pass_status, "note": f"TESS: difference-image flux deficit of {len(rows)} dip event(s) consistent with the target "
            f"(largest offset {worst['deficit_offset_px']:.2f} px, {worst['deficit_offset_sigma']:.1f} sigma; limit {s_thr:.0f} sigma and "
            f"{min_px:.2f} px){nb_txt}{edge_txt}", "events": rows, "n_sectors": n_aux, **edge_m}


def _lsq(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    dof = max(1, len(y) - X.shape[1])
    s2 = float(np.sum((y - X @ beta) ** 2)) / dof
    cov = s2 * np.linalg.pinv(X.T @ X)
    return beta, np.sqrt(np.maximum(0.0, np.diag(cov)))


def _neowise_neighbour(star: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    neo = star.get("neowise", {}) or {}
    ra_r, dec_r = neo.get("pm_residual_ra_arcsec"), neo.get("pm_residual_dec_arcsec")
    if not ra_r or not dec_r:
        return {"status": "not_evaluated", "note": "NEOWISE: no per-visit PM-track residual vectors"}
    mjd = np.asarray(neo.get("mjd", []), dtype=float)
    w2 = np.asarray(neo.get("w2_mjy", []), dtype=float)
    if not (len(mjd) == len(w2) == len(ra_r) == len(dec_r)):
        return {"status": "not_evaluated", "note": "NEOWISE: residual vectors not aligned with visits"}
    ok = np.asarray([a is not None and b is not None for a, b in zip(ra_r, dec_r)]) & (w2 > 0)
    n = int(ok.sum())
    if n < int(_it(cfg, "min_visits")):
        return {"status": "not_evaluated", "note": f"NEOWISE: only {n} visits with residual vectors"}
    x = np.asarray([v for v, k in zip(ra_r, ok) if k], dtype=float)
    y = np.asarray([v for v, k in zip(dec_r, ok) if k], dtype=float)
    m = mjd[ok]
    f = w2[ok]
    jyr = 2000.0 + (m - 51544.5) / 365.25
    t_yr = jyr - float(np.mean(jyr))
    phase = ((m - m[0]) / 365.25) % 1.0
    parity = np.where((phase >= 0.25) & (phase < 0.75), 0.5, -0.5)  # ascending vs descending scan visits
    wfrac = f / float(np.median(f)) - 1.0
    X = np.column_stack([np.ones(n), parity, t_yr, wfrac])
    bx, sx = _lsq(X, x)
    by, sy = _lsq(X, y)
    if star.get("pmra_mas_yr") is None or star.get("pmdec_mas_yr") is None:
        return {"status": "not_evaluated", "note": "NEOWISE: no Gaia proper motion to define the PM track (pmra/pmdec missing)"}
    pmra = float(star["pmra_mas_yr"])
    pmdec = float(star["pmdec_mas_yr"])
    mu = math.hypot(pmra, pmdec) / 1000.0
    sig = _it(cfg, "neighbour_blend_sigma")
    fmax = _vt(cfg, "neighbour_blend_flux_frac")
    metrics: dict[str, Any] = {
        "n_visits": n, "median_residual_arcsec": round(float(np.median(np.hypot(x, y))), 4),
        "pm_arcsec_yr": round(mu, 4), "residual_trend_arcsec_yr": [_r(bx[2]), _r(by[2])],
        "residual_per_w2_frac_arcsec": [_r(bx[3]), _r(by[3])], "scan_parity_offset_arcsec": [_r(bx[1]), _r(by[1])],
    }
    flags = []
    track = mu * float(np.ptp(jyr))
    if track >= _it(cfg, "neighbour_min_track_arcsec") and mu > 0:
        u = np.asarray([pmra, pmdec]) / (mu * 1000.0)
        v_along = float(u[0] * bx[2] + u[1] * by[2])
        v_err = math.sqrt((u[0] * sx[2]) ** 2 + (u[1] * sy[2]) ** 2)
        f_static = -v_along / mu
        f_err = max(1e-9, v_err / mu)
        metrics.update({"static_blend_frac": round(f_static, 4), "static_blend_frac_err": round(f_err, 4), "track_arcsec": round(track, 3)})
        if f_static >= fmax and f_static / f_err >= sig:
            flags.append(f"residual lags the proper-motion track: stationary flux share {f_static:.2f} +/- {f_err:.2f} (>= {fmax:.2f})")
        a_txt = f"along-track residual trend implies a stationary flux share {f_static:.3f} +/- {f_err:.3f} over a {track:.2f}\" track"
    else:
        a_txt = f"along-track test not applicable (track {track:.2f}\" < {_it(cfg, 'neighbour_min_track_arcsec'):.2f}\")"
    k = np.asarray([bx[3], by[3]])
    k_err = np.maximum(1e-9, np.asarray([sx[3], sy[3]]))
    k_sig = float(np.hypot(*(k / k_err)))
    w_rng = float(np.ptp(wfrac))
    pull_meas = float(np.hypot(*k)) * w_rng
    pull_err = float(np.hypot(*(k_err))) * w_rng
    i_pk = int(np.argmax(np.abs(wfrac)))
    d_pk = mu * max(0.0, float(jyr[i_pk]) - 2010.5)
    pull_pred = abs(float(f[i_pk] - np.min(f)) / float(f[i_pk])) * d_pk if f[i_pk] > 0 else 0.0
    metrics.update({"residual_w2_sigma": round(k_sig, 2), "measured_pull_arcsec": round(pull_meas, 4), "measured_pull_err_arcsec": round(pull_err, 4),
                    "stationary_source_pull_predicted_arcsec": round(pull_pred, 4), "star_offset_from_2010_arcsec_at_peak": round(d_pk, 3)})
    if k_sig >= sig:
        flags.append(f"residual correlates with W2 flux ({k_sig:.1f} sigma, {pull_meas:.3f}\" pull over the W2 range)")
    b_txt = (f"residual-W2 dependence {k_sig:.1f} sigma ({pull_meas:.3f} +/- {pull_err:.3f}\" over the W2 range; a stationary source at the "
             f"star's 2010.5 position carrying the variable flux would pull {pull_pred:.3f}\")")
    head = f"NEOWISE: {n} visits, median PM-track residual {metrics['median_residual_arcsec']:.3f}\"; {a_txt}; {b_txt}"
    if flags:
        return {"status": "fail", "note": head + "; blended stationary source: " + "; ".join(flags), "metrics": metrics}
    tail = ""
    if pull_err > 0 and pull_pred - pull_meas >= sig * pull_err:
        tail = "; the variable flux moves with the star"
    return {"status": "pass", "note": head + tail, "metrics": metrics}


def check_neighbour_signal(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    neo = star.get("neowise", {}) or {}
    statuses: list[str] = []
    notes: list[str] = []
    metrics: dict[str, Any] = {}
    legacy = {k: bool(d[k]) for d, k in ((lc, "neighbour_eb_centroid_shift"), (neo, "w1_w2_neighbour_blend")) if k in d}
    if legacy:
        metrics["explicit_flags"] = legacy
        if any(legacy.values()):
            statuses.append("fail")
            notes.append("Dip/IR variation tracks neighbouring star within TESS/WISE beam (explicit flag)")
        else:
            statuses.append("pass")
            notes.append("Explicit flags: aperture photocenter stable (no neighbour contamination)")
    tess = _tess_neighbour(lc, cfg)
    neow = _neowise_neighbour(star, cfg)
    for part, key in ((tess, "tess"), (neow, "neowise")):
        statuses.append(part["status"])
        notes.append(part["note"])
        metrics[key] = {k: v for k, v in part.items() if k != "note"}
    if lc.get("tess_neighbours"):
        metrics["tess_gaia_neighbours"] = lc["tess_neighbours"]
    status = _combine(statuses)
    if status == "pass" and tess.get("edge_blocked"):
        status = "not_evaluated"
        notes.insert(0, "A TESS dip candidate inside the detrend edge zone is not evaluated, so the check cannot pass")
    elif status == "not_evaluated" and not tess.get("edge_blocked"):
        notes.insert(0, "Neither TESS centroid series nor NEOWISE residual vectors available")
    return make_vet_result(out_dir, check="neighbour_signal", status=status, note=" | ".join(notes), metrics=metrics)


# ------------------------------------------------------------------------------------- sso_crossing


def check_sso_crossing(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    if "sso_crossing_t0" in lc:  # explicit (fixture) flag: present means evaluated
        hit = lc.get("sso_crossing_t0")
        return make_vet_result(
            out_dir, check="sso_crossing", status="fail" if hit is not None else "pass",
            note="No Solar System body crossing TESS aperture at event time (explicit flag)" if hit is None else f"Asteroid crossing aperture at {hit}",
            metrics={"sso_crossing_t0": hit, "explicit_flag": True},
        )
    secs, _ = _analysis_sectors(lc, include_legacy=False)
    if not secs:
        return make_vet_result(out_dir, check="sso_crossing", status="not_evaluated",
                               note="No TESS light curve loaded, so no event time to test" + _tess_log(lc))
    events, n_usable = tess_dip_events(lc, cfg)
    thr = _od(cfg, "min_dip_depth_frac")
    if n_usable == 0:
        return make_vet_result(out_dir, check="sso_crossing", status="not_evaluated",
                               note="TESS sectors loaded but none has enough filled bins to locate an event" + _tess_log(lc))
    edge_c = tess_edge_zone_candidates(lc, cfg)
    if not events:
        # The hypothesis tested here is "a Solar System body caused the flagged TESS dip". With TESS data actually
        # loaded and analysed for every fetched sector and no coincident dip above the channel threshold, there is no
        # TESS event an SSO could have caused, so the check legitimately passes. Without TESS data it stays
        # not_evaluated (above), so a pass never stands in for missing data; nor does a dip candidate in the detrend
        # edge zone, where the running median may have erased the event.
        if edge_c:
            return make_vet_result(out_dir, check="sso_crossing", status="not_evaluated",
                                   note=f"No dip event above threshold ({thr:.3f}) found by the detrend in {n_usable} TESS sector(s), but "
                                        + _edge_note(edge_c, cfg),
                                   metrics={"n_sectors": n_usable, "events": [], "detrend_edge_zone_candidates": edge_c})
        return make_vet_result(out_dir, check="sso_crossing", status="pass",
                               note=f"No dip event above threshold ({thr:.3f}) in {n_usable} TESS sector(s); no SSO-caused TESS event",
                               metrics={"n_sectors": n_usable, "events": [], "detrend_edge_zone_candidates": []})
    records = (lc.get("sso_queries") or {}).get("events") or []
    radius = _vt(cfg, "sso_tess_radius_px") * _vt(cfg, "tess_pixel_arcsec")
    g_mag = (star.get("mags") or {}).get("G")
    min_ratio = _vt(cfg, "sso_min_flux_ratio_of_depth")
    table: list[dict[str, Any]] = []
    missing: list[str] = []
    failed: list[str] = []
    hits: list[str] = []
    for ev in events:
        rec = next((r for r in records if r.get("sector") == ev["sector"] and abs(float(r.get("btjd") or -1e9) - ev["btjd"]) < 1e-3), None)
        if rec is None:
            missing.append(f"s{ev['sector']} BTJD {ev['btjd']:.3f}")
            continue
        if rec.get("status") != "ok":
            failed.append(f"s{ev['sector']} BTJD {ev['btjd']:.3f}: {rec.get('status')} ({'; '.join(rec.get('query_log') or [])[:160]})")
            continue
        relevant = []
        for b in rec.get("bodies") or []:
            v = b.get("v_mag")
            ratio = None if (v is None or g_mag is None) else 10 ** (-0.4 * (float(v) - float(g_mag)))
            if ratio is None or ratio >= min_ratio * float(ev["depth"]):
                relevant.append({"name": b.get("name"), "sep_arcsec": b.get("sep_arcsec"), "v_mag": v, "flux_ratio": _r(ratio, 5)})
        table.append({"sector": ev["sector"], "btjd": ev["btjd"], "depth": ev["depth"], "n_bodies": len(rec.get("bodies") or []), "relevant": relevant})
        if relevant:
            hits.append(f"s{ev['sector']} BTJD {ev['btjd']:.3f}: {relevant[0]['name']} at {relevant[0]['sep_arcsec']}\"")
    metrics = {"radius_arcsec": radius, "observer": "C57", "events": table, "missing": missing, "query_failed": failed,
               "detrend_edge_zone_candidates": edge_c}
    if hits:
        return make_vet_result(out_dir, check="sso_crossing", status="fail",
                               note=f"Solar System body within {radius:.0f}\" at TESS event time: " + "; ".join(hits), metrics=metrics)
    if missing or failed:
        why = []
        if missing:
            why.append("no SkyBoT query recorded for " + ", ".join(missing))
        if failed:
            why.append("SkyBoT query failed for " + "; ".join(failed))
        return make_vet_result(out_dir, check="sso_crossing", status="not_evaluated", note="; ".join(why), metrics=metrics)
    if edge_c:
        return make_vet_result(out_dir, check="sso_crossing", status="not_evaluated", metrics=metrics,
                               note=f"SkyBoT (TESS, C57) found no relevant body at {len(table)} event(s), but " + _edge_note(edge_c, cfg))
    return make_vet_result(
        out_dir, check="sso_crossing", status="pass",
        note=(f"SkyBoT (TESS, C57) at {len(table)} event(s): no body within {radius:.0f}\" bright enough (flux >= {min_ratio:.2f} x depth "
              f"relative to G={g_mag}) at the event times"),
        metrics=metrics,
    )


# ------------------------------------------------------------------------------------ known_variable


KNOWN_VARIABLE_SIMBAD_TYPES = frozenset(
    {
        "EB*", "Al*", "bL*", "WU*",
        "LPV*", "LP*", "Mira", "Mi*", "SR*",
        "CV*", "DN*", "NL*", "No*",
        "RR*", "cC*", "Ce*", "Cep", "dS*", "RV*", "bC*", "gD*", "SX*",
    }
)


def _light_detrend(t: np.ndarray, f: np.ndarray, gap: float) -> np.ndarray:
    """Per continuous segment: divide by the median and remove a linear trend only, preserving rotation."""
    out = np.ones_like(f)
    breaks = np.where(np.diff(t) > gap)[0] + 1
    for seg in np.split(np.arange(len(t)), breaks):
        ys = f[seg] / max(1e-9, float(np.median(f[seg])))
        if len(seg) >= 10:
            p = np.polyfit(t[seg] - t[seg][0], ys, 1)
            ys = ys - np.polyval(p, t[seg] - t[seg][0]) + 1.0
        out[seg] = ys
    return out


def _sector_period_search(sec: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any] | None:
    from astropy.timeseries import BoxLeastSquares, LombScargle

    t = np.asarray(sec.get("time_btjd", []), dtype=float)
    f = np.asarray(sec.get("flux_qlp_raw", []), dtype=float)
    if len(t) < 50 or len(f) != len(t):
        return None
    y = _light_detrend(t, f, _od(cfg, "tess_segment_gap_days"))
    span = float(t[-1] - t[0])
    p_min = _od(cfg, "rotation_min_period_days")
    p_max = min(_od(cfg, "rotation_max_period_days"), 0.5 * span)
    out: dict[str, Any] = {"sector": int(sec.get("sector", -1))}
    if p_max > p_min:
        ls = LombScargle(t, y)
        freq, power = ls.autopower(minimum_frequency=1.0 / p_max, maximum_frequency=1.0 / p_min, samples_per_peak=10)
        j = int(np.argmax(power))
        fap = float(ls.false_alarm_probability(power[j], method="baluev", minimum_frequency=1.0 / p_max, maximum_frequency=1.0 / p_min))
        amp = float(np.percentile(y, 95) - np.percentile(y, 5))
        out.update({"ls_period_days": round(1.0 / float(freq[j]), 4), "ls_power": round(float(power[j]), 4), "ls_fap": float(f"{fap:.3e}"),
                    "amplitude_5_95_frac": round(amp, 5)})
    bls_min = _vt(cfg, "lc_bls_min_period_days")
    bls_max = min(_vt(cfg, "lc_bls_max_period_days"), 0.5 * span)
    if bls_max > bls_min:
        model = BoxLeastSquares(t, y)
        periods = np.linspace(bls_min, bls_max, 200)
        res = model.power(periods, _vt(cfg, "lc_bls_duration_days"))
        snr = (float(np.max(res.power)) - float(np.median(res.power))) / max(1e-6, float(np.std(res.power)))
        jb = int(np.argmax(res.power))
        period, t0, dur, depth = float(res.period[jb]), float(res.transit_time[jb]), float(res.duration[jb]), float(res.depth[jb])
        # A periodic eclipse needs >= 2 observed eclipses reaching half the model depth, and at least half of the
        # predicted eclipse epochs that have data must show it (a single aperiodic dip folds into one or two epochs only)
        n_ecl = n_pred = 0
        k0 = int(math.floor((t[0] - t0) / period)) - 1
        for k in range(k0, k0 + int(span / period) + 3):
            tc = t0 + k * period
            sel = np.abs(t - tc) <= 0.5 * dur
            if sel.sum() >= 2:
                n_pred += 1
                if (1.0 - float(np.mean(y[sel]))) >= _vt(cfg, "lc_eb_eclipse_min_depth_frac_of_model") * depth:
                    n_ecl += 1
        eb = bool(snr > _od(cfg, "periodic_bls_veto_snr") and depth >= _od(cfg, "periodic_eb_min_depth_frac")
                  and n_ecl >= int(_vt(cfg, "lc_eb_min_eclipses")) and n_ecl >= _vt(cfg, "lc_eb_min_eclipse_epoch_frac") * n_pred)
        out.update({"bls_period_days": round(period, 4), "bls_snr": round(snr, 2), "bls_depth": round(depth, 4), "bls_n_eclipses": n_ecl,
                    "bls_n_predicted": n_pred, "periodic_eb": eb})
    return out


def check_known_variable(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    lc = star.get("lightcurves", {})
    secs_any, _ = _analysis_sectors(lc)
    sim = simbad_lookup(star)
    if not sim["evaluated"] and not secs_any:
        return make_vet_result(out_dir, check="known_variable", status="not_evaluated",
                               note=f"{simbad_note(sim)} and no TESS light curve to search for periodicity" + _tess_log(lc),
                               metrics={"simbad_query_status": sim["query_status"]})
    stype = sim["simbad_type"]
    t = np.asarray(lc.get("tess_time_btjd", []), dtype=float)
    f = np.asarray(lc.get("tess_flux_qlp", []), dtype=float)
    periodic_eb = False
    if len(t) > 30:
        from astropy.timeseries import BoxLeastSquares

        model = BoxLeastSquares(t, f)
        periods = np.linspace(_vt(cfg, "lc_bls_min_period_days"), _vt(cfg, "lc_bls_max_period_days"), 200)
        res = model.power(periods, _vt(cfg, "lc_bls_duration_days"))
        best_power = float(np.max(res.power))
        med_power = float(np.median(res.power))
        snr = (best_power - med_power) / max(1e-6, float(np.std(res.power)))
        periodic_eb = snr > _od(cfg, "periodic_bls_veto_snr") and float(np.min(f)) < 1.0 - _od(cfg, "periodic_eb_min_depth_frac")

    metrics: dict[str, Any] = {"simbad_type": stype, "simbad_query_status": sim["query_status"], "simbad_entry": sim["entry"],
                               "periodic_eb_detected": periodic_eb}
    rot_txt = ""
    secs, _ = _analysis_sectors(lc, include_legacy=False)
    if secs:
        per = [p for p in (_sector_period_search(s, cfg) for s in secs) if p]
        fap_max = _od(cfg, "rotation_fap_max")
        det = [p for p in per if p.get("ls_fap") is not None and p["ls_fap"] < fap_max]
        eb_secs = [p["sector"] for p in per if p.get("periodic_eb")]
        rot: dict[str, Any] = {"n_sectors_searched": len(per), "n_sectors_detected": len(det), "fap_max": fap_max}
        if det:
            pm = float(np.median([p["ls_period_days"] for p in det]))
            match = sum(1 for p in det if abs(p["ls_period_days"] - pm) <= _od(cfg, "rotation_period_match_frac") * pm)
            rot.update({"period_days": round(pm, 4), "best_fap": min(p["ls_fap"] for p in det), "n_sectors_consistent": match})
            rot_txt = (f"; rotation/spot modulation candidate P={pm:.3f} d (best FAP {rot['best_fap']:.1e}, {match}/{len(per)} sectors "
                       f"consistent; reported, not a veto)")
        else:
            rot_txt = f"; no Lomb-Scargle period with FAP < {fap_max:.0e} in {len(per)} TESS sector(s)"
        metrics["rotation"] = rot
        metrics["tess_period_search"] = per
        if eb_secs:
            periodic_eb = True
            metrics["periodic_eb_sectors"] = eb_secs
        metrics["periodic_eb_detected"] = periodic_eb
    else:
        metrics["rotation"] = {"n_sectors_searched": 0, "note": "no TESS sectors loaded" + _tess_log(lc)}

    gaia_var = (star.get("vetting_archive") or {}).get("gaia_phot_variable_flag")
    if gaia_var:
        metrics["gaia_phot_variable_flag"] = gaia_var
        rot_txt += f"; Gaia phot_variable_flag={gaia_var}"
    is_known_var = (sim["evaluated"] and stype in KNOWN_VARIABLE_SIMBAD_TYPES) or periodic_eb
    # A TESS period search can find an eclipsing binary but cannot show that a star is not a catalogued variable,
    # so a pass needs the SIMBAD lookup to have succeeded ("no entry" is a real negative answer).
    status = "fail" if is_known_var else ("pass" if sim["evaluated"] else "not_evaluated")
    return make_vet_result(
        out_dir,
        check="known_variable",
        status=status,
        note=f"{simbad_note(sim)}, periodic_EB_detected={periodic_eb}{rot_txt}",
        metrics=metrics,
    )


# ------------------------------------------------------------------------------------- second_survey


def _season_span(s: dict[str, Any], kind: str) -> tuple[float, float]:
    if s.get("mjd_first") is not None and s.get("mjd_last") is not None:
        return float(s["mjd_first"]), float(s["mjd_last"])
    yr = float(s["season"]) + (0.5 if kind == "july_june" else 0.0)
    lo = 51544.5 + (yr - 2000.0) * 365.25
    return lo, lo + 365.25


def _season_mid(s: dict[str, Any], kind: str) -> float:
    if s.get("mjd_mid") is not None:
        return float(s["mjd_mid"])
    lo, hi = _season_span(s, kind)
    return 0.5 * (lo + hi)


def _jyear(mjd: float) -> float:
    return 2000.0 + (mjd - 51544.5) / 365.25


def _band_series(lc: dict[str, Any]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Per survey/band season tables and nightly medians, plus the status of every survey that has none."""
    series: list[dict[str, Any]] = []
    statuses: dict[str, Any] = {}
    for key, sv in sorted((lc.get("optical_surveys") or {}).items()):
        kind = str(sv.get("season_kind") or ("july_june" if "ASAS" in str(sv.get("survey", key)).upper() else "calendar"))
        bands = sv.get("bands") or {}
        statuses[key] = {kk: sv.get(kk) for kk in ("survey", "n_epochs", "status", "query_log")}
        for band, b in sorted(bands.items()):
            if b.get("seasons"):
                series.append({"key": key, "survey": str(sv.get("survey", key)), "band": band, "kind": kind, "seasons": b["seasons"],
                               "nights": b.get("nights") or [], "night_sigma": b.get("nightly_robust_sigma_mag") or b.get("robust_sigma_mag"),
                               "n_epochs": b.get("n_epochs")})
    t = np.asarray(lc.get("ztf_time_mjd", []) or [], dtype=float)
    fz = np.asarray(lc.get("ztf_flux", []) or [], dtype=float)
    if not series and len(t) >= 10 and len(t) == len(fz) and np.all(fz > 0):
        mags = -2.5 * np.log10(fz / np.median(fz))
        seasons = []
        years = np.floor(np.asarray([_jyear(x) for x in t])).astype(int)
        for yr in sorted(set(years.tolist())):
            sel = years == yr
            mm, tt = mags[sel], t[sel]
            med = float(np.median(mm))
            err = 1.2533 * _robust_sigma(mm) / math.sqrt(len(mm)) if len(mm) >= 2 else None
            seasons.append({"season": int(yr), "n_nights": int(sel.sum()), "median_mag": med, "median_err_mag": err,
                            "mjd_first": float(tt.min()), "mjd_mid": float(np.median(tt)), "mjd_last": float(tt.max())})
        series.append({"key": "ztf_legacy", "survey": "ZTF (legacy normalised flux)", "band": "z", "kind": "calendar", "seasons": seasons,
                       "nights": [[float(a), float(b)] for a, b in zip(t, mags)], "night_sigma": _robust_sigma(mags), "n_epochs": len(t)})
        statuses["ztf_legacy"] = {"survey": "ZTF", "n_epochs": int(len(t)), "status": "ok"}
    return series, statuses


def _good_seasons(sr: dict[str, Any], cfg: PipelineConfig) -> list[dict[str, Any]]:
    min_n = int(_it(cfg, "second_survey_min_season_nights"))
    floor = _it(cfg, "second_survey_season_floor_mag")
    out = []
    for s in sr["seasons"]:
        if int(s.get("n_nights", 0) or 0) < min_n or s.get("median_mag") is None:
            continue
        err = s.get("median_err_mag")
        if err is None and sr.get("night_sigma"):
            err = 1.2533 * float(sr["night_sigma"]) / math.sqrt(int(s["n_nights"]))
        if err is None:
            continue
        lo, hi = _season_span(s, sr["kind"])
        out.append({"season": s["season"], "mag": float(s["median_mag"]), "err": math.hypot(float(err), floor),
                    "mid": _season_mid(s, sr["kind"]), "lo": lo, "hi": hi, "n_nights": int(s["n_nights"])})
    return out


def _w2_rel_mag(neo: dict[str, Any]) -> tuple[np.ndarray, np.ndarray] | None:
    mjd = np.asarray(neo.get("mjd", []) or [], dtype=float)
    w2 = np.asarray(neo.get("w2_mjy", []) or [], dtype=float)
    if len(mjd) < 3 or len(w2) != len(mjd) or np.any(w2 <= 0):
        return None
    order = np.argsort(mjd)
    return mjd[order], -2.5 * np.log10(w2[order] / np.median(w2))


def _ir_event_window(mjd: np.ndarray, m: np.ndarray, cfg: PipelineConfig) -> tuple[float, float]:
    """Visits deviating from the baseline (median of the first `lc_ir_event_baseline_visits` visits) by at least
    `lc_ir_event_dev_frac` of the maximum deviation span the IR event; a window shorter than
    `lc_ir_event_min_window_days` is padded by `lc_ir_event_pad_days` on both sides."""
    base = float(np.median(m[: min(int(_vt(cfg, "lc_ir_event_baseline_visits")), len(m))]))
    dev = np.abs(m - base)
    sel = np.where(dev >= _vt(cfg, "lc_ir_event_dev_frac") * float(np.max(dev)))[0]
    lo, hi = float(mjd[sel[0]]), float(mjd[sel[-1]])
    if hi - lo < _vt(cfg, "lc_ir_event_min_window_days"):
        pad = _vt(cfg, "lc_ir_event_pad_days")
        lo, hi = lo - pad, hi + pad
    return lo, hi


def _band_response(sr: dict[str, Any], good: list[dict[str, Any]], w2: tuple[np.ndarray, np.ndarray] | None, cfg: PipelineConfig) -> dict[str, Any]:
    mags = np.asarray([g["mag"] for g in good])
    errs = np.asarray([g["err"] for g in good])
    i_max, i_min = int(np.argmax(mags)), int(np.argmin(mags))
    rng = float(mags[i_max] - mags[i_min]) if len(mags) >= 2 else 0.0
    rng_sig = rng / max(1e-9, math.hypot(errs[i_max], errs[i_min])) if len(mags) >= 2 else 0.0
    significant = bool(rng_sig >= _it(cfg, "second_survey_change_sigma"))
    out: dict[str, Any] = {
        "survey": sr["survey"], "band": sr["band"], "n_seasons": len(good),
        "seasons": [[g["season"], round(g["mag"], 4), round(g["err"], 4)] for g in good],
        "range_mag": round(rng, 4), "range_sigma": round(rng_sig, 2), "significant_change": significant,
        "n_overlap_w2": 0, "pearson_r": None, "slope_dopt_dw2": None, "slope_err": None, "correlation_tested": False,
    }
    if w2 is not None:
        wm, wv = w2
        ov = [g for g in good if wm[0] <= g["mid"] <= wm[-1]]
        out["n_overlap_w2"] = len(ov)
        if len(ov) >= int(_it(cfg, "second_survey_min_overlap_seasons")):
            x = np.interp([g["mid"] for g in ov], wm, wv)
            y = np.asarray([g["mag"] for g in ov])
            e = np.asarray([g["err"] for g in ov])
            r = _pearson(x, y)
            if r is not None:
                W = 1.0 / e**2
                X = np.column_stack([np.ones(len(x)), x])
                cov = np.linalg.pinv(X.T @ (X * W[:, None]))
                beta = cov @ (X.T @ (W * y))
                chi2r = float(np.sum(W * (y - X @ beta) ** 2)) / max(1, len(y) - 2)
                b_err = math.sqrt(max(0.0, cov[1, 1]) * max(1.0, chi2r))
                out.update({"pearson_r": round(r, 3), "slope_dopt_dw2": round(float(beta[1]), 4), "slope_err": round(b_err, 4),
                            "correlation_tested": True, "w2_range_mag_at_seasons": round(float(np.ptp(x)), 4)})
    if not significant:
        out["response"] = "none"
    elif out["correlation_tested"] and abs(out["pearson_r"]) >= _it(cfg, "second_survey_min_abs_corr") and \
            abs(out["slope_dopt_dw2"]) >= _it(cfg, "second_survey_slope_sigma") * max(1e-9, out["slope_err"]):
        # Magnitudes: optical fainter (larger) while W2 brighter (smaller) gives a negative slope
        out["response"] = "dimming_with_ir_brightening" if out["slope_dopt_dw2"] < 0 else "brightening_with_ir_brightening"
    else:
        out["response"] = "uncorrelated_variation"
    return out


def _coverage(good: list[dict[str, Any]], lo: float, hi: float) -> tuple[float, float, float]:
    ivs = sorted((max(lo, g["lo"]), min(hi, g["hi"])) for g in good if g["hi"] >= lo and g["lo"] <= hi)
    covered = 0.0
    cur_lo, cur_hi = None, None
    for a, b in ivs:
        if cur_hi is None or a > cur_hi:
            if cur_hi is not None:
                covered += cur_hi - cur_lo  # type: ignore[operator]
            cur_lo, cur_hi = a, b
        else:
            cur_hi = max(cur_hi, b)
    if cur_hi is not None:
        covered += cur_hi - cur_lo  # type: ignore[operator]
    first = ivs[0][0] if ivs else float("nan")
    last = ivs[-1][1] if ivs else float("nan")
    return covered / max(1e-9, hi - lo), first, last


def _ir_context(star: dict[str, Any], series: list[dict[str, Any]], cfg: PipelineConfig) -> dict[str, Any] | None:
    neo = star.get("neowise", {}) or {}
    var = _ir_variation(neo, cfg)
    if not var.get("active"):
        return None
    w2 = _w2_rel_mag(neo)
    if w2 is None:
        return None
    lo, hi = _ir_event_window(w2[0], w2[1], cfg)
    win_txt = f"{_jyear(lo):.1f}-{_jyear(hi):.1f}"
    bands, cov_txt = [], []
    for sr in series:
        good = _good_seasons(sr, cfg)
        frac, first, last = _coverage(good, lo, hi)
        if len(good) < 2 or frac <= 0:
            cov_txt.append(f"{sr['survey']} {sr['band']}: no season overlapping")
            continue
        resp = _band_response(sr, good, w2, cfg)
        resp.update({"event_window_coverage_frac": round(frac, 3), "covered_from_jyear": round(_jyear(first), 2), "covered_to_jyear": round(_jyear(last), 2)})
        bands.append(resp)
        cov_txt.append(f"{sr['survey']} {sr['band']} covers {_jyear(first):.1f}-{_jyear(last):.1f} ({frac * 100:.0f}%)")
    ctx: dict[str, Any] = {"event_window_mjd": [round(lo, 1), round(hi, 1)], "event_window_jyear": win_txt, "ir_variation": var, "bands": bands}
    if not bands:
        ctx.update({"status": "not_evaluated", "note": f"IR event window {win_txt}: no optical season overlaps it ({'; '.join(cov_txt)})"})
        return ctx
    sig = [b for b in bands if b["significant_change"]]
    comp = _it(cfg, "second_survey_comparable_frac")
    contra = [b for b in sig if b["response"] == "brightening_with_ir_brightening" and b["slope_dopt_dw2"] >= comp]
    if not sig:
        overall = "none"
    else:
        labels = [b["response"] for b in sig]
        prio = ["brightening_with_ir_brightening", "uncorrelated_variation", "dimming_with_ir_brightening"]
        overall = max(prio, key=lambda lab: (labels.count(lab), -prio.index(lab)))
    ctx["optical_response"] = overall

    def _btxt(b: dict[str, Any]) -> str:
        s = f"{b['survey']} {b['band']} range {b['range_mag']:.3f} mag ({b['range_sigma']:.1f} sigma)"
        if b["correlation_tested"]:
            s += f", r={b['pearson_r']:+.2f} vs W2, dOpt/dW2={b['slope_dopt_dw2']:+.3f}+/-{b['slope_err']:.3f}"
        return s + f" [{b['response']}]"

    detail = "; ".join(_btxt(b) for b in bands)
    head = f"IR event window {win_txt} ({'; '.join(cov_txt)}): optical_response={overall}"
    if contra:
        ctx.update({"status": "fail", "note": head + f" - optical brightens with W2 by >= {comp:.2f} x the IR change, contradicting circumstellar "
                    f"dust (stellar/accretion variability or a blended variable): " + detail})
    elif overall == "none":
        ctx.update({"status": "pass", "note": head + f" - no seasonal optical change at >= {_it(cfg, 'second_survey_change_sigma'):.0f} sigma: " + detail})
    else:
        ctx.update({"status": "pass", "note": head + " - significant optical change detected that does not contradict a dust interpretation: " + detail})
    return ctx


def _dip_context(lc: dict[str, Any], series: list[dict[str, Any]], cfg: PipelineConfig) -> dict[str, Any] | None:
    secs, _ = _analysis_sectors(lc)
    events = []
    edge_c: list[dict[str, Any]] = []
    for sec in secs:
        a = analyze_tess_sector(sec, cfg)
        if not a or not a["usable"]:
            continue
        edge_c.extend(_sector_edge_candidates(a))
        if a["summary"]["deepest_coincident"]["coincident_exceeds_min_depth"]:
            events.append(dict(a["summary"]["deepest_coincident"], sector=a["summary"]["sector"]))
    if not events:
        if edge_c:
            return {"tess_events": [], "detrend_edge_zone_candidates": edge_c, "status": "not_evaluated",
                    "note": "No TESS dip event found by the detrend, but " + _edge_note(edge_c, cfg)}
        return None
    pad = _od(cfg, "second_survey_simultaneous_pad_days")
    ratio = _od(cfg, "second_survey_simultaneous_min_depth_ratio")
    nsig = _od(cfg, "second_survey_night_dip_sigma")
    min_nights = int(_od(cfg, "second_survey_min_dip_nights"))
    thr = _od(cfg, "min_dip_depth_frac")
    simul, contra, confirm_sim, indep = [], [], [], []
    for sr in series:
        nights = np.asarray(sr["nights"], dtype=float).reshape(-1, 2) if sr["nights"] else np.zeros((0, 2))
        sig = float(sr.get("night_sigma") or 0.0)
        if len(nights) < 5 or sig <= 0:
            continue
        base = float(np.median(nights[:, 1]))
        for ev in events:
            lo = ev["btjd_start"] + BTJD_TO_MJD - pad
            hi = ev["btjd_end"] + BTJD_TO_MJD + pad
            sel = (nights[:, 0] >= lo) & (nights[:, 0] <= hi)
            if not sel.any():
                continue
            dmag = float(np.max(nights[sel, 1]) - base)
            d_frac = 1.0 - 10 ** (-0.4 * dmag)
            sig_frac = sig * math.log(10) / 2.5
            rec = {"survey": sr["survey"], "band": sr["band"], "sector": ev["sector"], "btjd": ev["btjd"], "tess_depth": ev["depth_coincident"],
                   "survey_depth": round(d_frac, 4), "survey_sigma": round(sig_frac, 4), "n_nights": int(sel.sum())}
            simul.append(rec)
            if d_frac + nsig * sig_frac < ratio * ev["depth_coincident"]:
                contra.append(rec)
            else:
                confirm_sim.append(rec)
        faint = (nights[:, 1] - base) >= max(nsig * sig, 2.5 * math.log10(1.0 / (1.0 - thr)))
        run, best = 0, 0
        for fl in faint:
            run = run + 1 if fl else 0
            best = max(best, run)
        if best >= min_nights:
            indep.append(f"{sr['survey']} {sr['band']} ({best} consecutive nights >= {thr:.2f} fainter)")
    ev_txt = ", ".join(f"s{e['sector']} BTJD {e['btjd']:.2f} (depth {e['depth_coincident']:.3f})" for e in events)
    ctx: dict[str, Any] = {"tess_events": events, "simultaneous": simul, "independent_dimming": indep, "detrend_edge_zone_candidates": edge_c}
    if contra:
        c = contra[0]
        ctx.update({"status": "fail", "note": f"TESS events {ev_txt}: simultaneous {c['survey']} {c['band']} shows depth {c['survey_depth']:.3f} "
                    f"(+/- {c['survey_sigma']:.3f}) < {ratio:.2f} x TESS depth {c['tess_depth']:.3f}; the dip is not on the sky"})
    elif confirm_sim:
        c = confirm_sim[0]
        ctx.update({"status": "pass", "note": f"TESS events {ev_txt}: simultaneous {c['survey']} {c['band']} depth {c['survey_depth']:.3f} is "
                    f"consistent with TESS {c['tess_depth']:.3f}"})
    elif indep:
        ctx.update({"status": "pass", "note": f"TESS events {ev_txt}: no simultaneous survey night, but an independent survey confirms dimming "
                    f"episodes: {'; '.join(indep)}"})
    else:
        ctx.update({"status": "not_evaluated", "note": f"TESS events {ev_txt}: no survey night within the event bins and no independent "
                    f"dimming episode, so ground surveys do not constrain the dip"})
    if ctx["status"] == "pass" and edge_c:
        ctx.update({"status": "not_evaluated", "note": ctx["note"] + "; but " + _edge_note(edge_c, cfg)})
    return ctx


def check_second_survey(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Do independent optical surveys confirm, contradict, or fail to constrain the flagged signal?

    Two signal contexts are evaluated from the data present (the check does not know which channel flagged):
    the NEOWISE IR variation (when it passes the ir_timedomain thresholds) and any TESS dip event above
    `min_dip_depth_frac`. Any failing context fails the check; any context that cannot be constrained keeps it
    `not_evaluated`; it passes only when every applicable context passes.
    """
    lc = star.get("lightcurves", {})
    series, statuses = _band_series(lc)
    ranges = []
    for sr in series:
        good = _good_seasons(sr, cfg)
        if len(good) >= 2:
            ranges.append(max(g["mag"] for g in good) - min(g["mag"] for g in good))
    not_usable = sorted(k for k in statuses if not any(sr["key"] == k for sr in series))
    surv_txt = "; ".join(f"{v.get('survey', k)}: {v.get('status')}" + (f" ({'; '.join(map(str, v.get('query_log') or []))[:160]})"
                                                                       if v.get("status") not in ("ok", None) else "")
                         for k, v in statuses.items())
    metrics: dict[str, Any] = {
        "n_surveys": len({sr["survey"] for sr in series}),
        "max_seasonal_range_mag": round(max(ranges), 4) if ranges else None,
        "survey_status": statuses,
    }
    tail = f" (not usable: {not_usable})" if not_usable else ""
    if not series:
        return make_vet_result(out_dir, check="second_survey", status="not_evaluated",
                               note="No ZTF/ASAS-SN optical coverage at target position" + (f" ({surv_txt})" if surv_txt else ""), metrics=metrics)
    ctxs = {"ir_timedomain": _ir_context(star, series, cfg), "optical_dips": _dip_context(lc, series, cfg)}
    applicable = {k: v for k, v in ctxs.items() if v is not None}
    metrics["contexts"] = applicable
    if "ir_timedomain" in applicable and applicable["ir_timedomain"].get("optical_response"):
        metrics["optical_response"] = applicable["ir_timedomain"]["optical_response"]
    if not applicable:
        return make_vet_result(out_dir, check="second_survey", status="not_evaluated",
                               note=("No flagged signal to test: NEOWISE W2 variation below ir_timedomain thresholds (or absent) and no TESS dip "
                                     f"above {_od(cfg, 'min_dip_depth_frac'):.3f}; survey coverage: {surv_txt}{tail}"), metrics=metrics)
    sts = [c["status"] for c in applicable.values()]
    status = "fail" if "fail" in sts else ("not_evaluated" if "not_evaluated" in sts else "pass")
    note = " | ".join(f"{k}: {c['note']}" for k, c in applicable.items()) + tail
    return make_vet_result(out_dir, check="second_survey", status=status, note=note, metrics=metrics)


def run_lc_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [
        check_neighbour_signal,
        check_systematics_window,
        check_pipeline_agreement,
        check_sso_crossing,
        check_known_variable,
        check_second_survey,
    ]
    return [fn(star, cfg, out_dir) for fn in checks]
