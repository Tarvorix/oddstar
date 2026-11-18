"""Channel 3: Aperiodic and asymmetric optical dimming events (TESS, Kepler/K2, ZTF, ASAS-SN)."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from astropy.timeseries import BoxLeastSquares
from scipy.ndimage import median_filter

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec


# Fewer points than this in a segment (after the edge margins are removed) is not searched
MIN_SEARCH_POINTS = 20
COVERAGE_CATEGORIES = ("searched", "no_coverage", "failed", "not_fetched", "no_light_curve")


def tess_search_status(star: dict[str, Any], cfg: PipelineConfig) -> str:
    """Whether the dip search had data for this star: "searched" or the reason it did not.

    "searched" needs at least one light-curve segment with `MIN_SEARCH_POINTS` points left inside the edge margins,
    the same condition `measure` applies. Otherwise the survey fetch record says why (`no_coverage`, `failed`,
    `not_fetched`); a star with no fetch record and no light curve is `no_light_curve`. A star that was not
    searched must never be counted in a null result.
    """
    ch = cfg.channels["optical_dips"]
    lc = star.get("lightcurves", {}) or {}
    t = np.asarray(lc.get("tess_time_btjd", []), dtype=float)
    sec = np.asarray(lc.get("tess_sector", []), dtype=float)
    if len(sec) != len(t):
        sec = np.zeros(len(t))
    if len(t) >= MIN_SEARCH_POINTS and len(lc.get("tess_flux_qlp", [])) == len(t):
        margin = float(ch["systematics_edge_margin_days"])
        cuts = np.nonzero((np.diff(t) > float(ch["tess_segment_gap_days"])) | (np.diff(sec) != 0))[0] + 1
        for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, len(t)]):
            ts = t[lo:hi]
            if len(ts) >= MIN_SEARCH_POINTS and int(((ts >= ts[0] + margin) & (ts <= ts[-1] - margin)).sum()) >= MIN_SEARCH_POINTS:
                return "searched"
    status = str((star.get("tess_fetch") or {}).get("status") or "")
    if status.startswith("failed"):
        return "failed"
    if status in {"no_coverage", "not_fetched"}:
        return status
    return "no_light_curve"


def search_coverage_summary(stars: list[dict[str, Any]], cfg: PipelineConfig) -> dict[str, Any]:
    """Counts of searched and not-searched stars, in the layout the infrared time-domain channel uses."""
    ids: dict[str, list[str]] = {c: [] for c in COVERAGE_CATEGORIES if c != "searched"}
    n_searched = 0
    for s in stars:
        cat = tess_search_status(s, cfg)
        if cat == "searched":
            n_searched += 1
        else:
            ids[cat].append(str(s["source_id"]))
    out: dict[str, Any] = {"n_input": len(stars), "n_searched": n_searched, "n_not_searched": len(stars) - n_searched}
    for cat, lst in ids.items():
        out[f"n_{cat}"] = len(lst)
    out["not_searched_ids"] = {cat: lst for cat, lst in ids.items() if lst}
    return out


class OpticalDipsChannel:
    name = "optical_dips"

    def __init__(self, cfg: PipelineConfig, root: Path) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["optical_dips"]
        self.root = root

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec:
        chash = compute_config_hash(cfg.model_dump())[:12]
        spec, _ = build_parent_sample(
            self.name,
            cfg,
            chash,
            self.root,
            extra_stars=kwargs.get("extra_stars"),
            limit=kwargs.get("limit"),
        )
        return spec

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        """Search one star's TESS light curve for a dip, segment by segment.

        A light curve is cut wherever the sector changes or cadences are more than `tess_segment_gap_days` apart
        (orbit gaps, sector gaps). Each segment is detrended and searched alone, because a running median or a
        box filter laid across a gap mixes unrelated data. The first and last `systematics_edge_margin_days` of a
        segment are not searched: thermal settling and scattered light sit there, and the vetting check
        `systematics_window` cannot clear an event in that zone anyway. `dip_depth` is the median deficit inside
        the best box, so one bad cadence is not a dip.
        """
        lc = star.get("lightcurves", {})
        t = np.asarray(lc.get("tess_time_btjd", []), dtype=float)
        f_qlp = np.asarray(lc.get("tess_flux_qlp", []), dtype=float)
        f_cut = np.asarray(lc.get("tess_flux_tesscut", []), dtype=float)
        sec = np.asarray(lc.get("tess_sector", []), dtype=float)
        ztf_t = np.asarray(lc.get("ztf_time_mjd", []), dtype=float)
        ztf_f = np.asarray(lc.get("ztf_flux", []), dtype=float)
        if len(sec) != len(t):
            sec = np.zeros(len(t))
        if len(f_qlp) != len(t):
            t, f_qlp, sec = t[:0], f_qlp[:0], sec[:0]

        vals = {"dip_snr": 0.0, "dip_depth": 0.0, "asymmetry": 0.0, "duration_hours": 0.0, "periodic_bls_snr": 0.0,
                "secular_fade_snr": 0.0, "pipeline_agree": 0.0, "n_points_searched": 0.0, "n_segments_searched": 0.0,
                "best_time_btjd": 0.0, "best_sector": 0.0}
        gap = float(self.ch_cfg["tess_segment_gap_days"])
        margin = float(self.ch_cfg["systematics_edge_margin_days"])
        best: dict[str, Any] | None = None
        parts: list[tuple[float, np.ndarray, np.ndarray]] = []  # (sector, times, detrended flux) of each searched segment
        if len(t) >= MIN_SEARCH_POINTS:
            cuts = np.nonzero((np.diff(t) > gap) | (np.diff(sec) != 0))[0] + 1
            for lo, hi in zip(np.r_[0, cuts], np.r_[cuts, len(t)]):
                ts, fs = t[lo:hi], f_qlp[lo:hi]
                inner = (ts >= ts[0] + margin) & (ts <= ts[-1] - margin)
                if len(ts) < MIN_SEARCH_POINTS or int(inner.sum()) < MIN_SEARCH_POINTS:
                    continue
                dt = float(np.median(np.diff(ts)))
                win_pts = max(15, int(round(float(self.ch_cfg["detrend_window_days"]) / max(1e-4, dt))) | 1)
                # Robust sliding median baseline preserving multi-day dips
                baseline = median_filter(fs, size=win_pts, mode="nearest")
                # Iterative upper clip so deep wide dips don't pull the baseline down
                resid_init = fs - baseline
                noise_est = max(1e-4, 1.4826 * float(np.median(np.abs(resid_init))))
                masked = np.where(resid_init < -2.5 * noise_est, np.nan, fs)
                valid = ~np.isnan(masked)
                if np.sum(valid) >= 10:
                    interp_f = np.interp(ts, ts[valid], masked[valid])
                    baseline = median_filter(interp_f, size=win_pts, mode="nearest")
                rel = fs / np.maximum(1e-4, baseline)
                sigma = max(1e-4, 1.4826 * float(np.median(np.abs(rel[inner] - np.median(rel[inner])))))
                parts.append((float(sec[lo]), ts[inner], rel[inner]))
                vals["n_points_searched"] += float(inner.sum())
                vals["n_segments_searched"] += 1.0
                # Matched-filter boxcar search over dip durations (6h to 4 days), box centres inside the margins only
                for dur_days in (0.25, 0.5, 1.0, 2.0, 3.5):
                    w = max(2, int(round(dur_days / max(1e-4, dt))))
                    if w > len(ts):
                        continue
                    smoothed_deficit = np.convolve(1.0 - rel, np.ones(w) / w, mode="same")
                    snr_curve = np.where(inner, smoothed_deficit * math.sqrt(w) / sigma, -np.inf)
                    idx = int(np.argmax(snr_curve))
                    if best is None or snr_curve[idx] > best["snr"]:
                        a, b = max(0, idx - w // 2), min(len(ts), idx - w // 2 + w)
                        best = {"snr": float(snr_curve[idx]), "dur_h": dur_days * 24.0, "lo": int(lo), "a": a, "b": b,
                                "rel": rel, "t": ts, "idx": idx, "sector": float(sec[lo])}

        if best is not None:
            rel, a, b = best["rel"], best["a"], best["b"]
            dip_depth = max(0.0, float(1.0 - np.median(rel[a:b])))
            # Ingress vs egress asymmetry around the deepest point of the best box
            deep_idx = a + int(np.argmin(rel[a:b]))
            half_lvl = 1.0 - 0.5 * max(0.0, float(1.0 - rel[deep_idx]))
            left = deep_idx
            while left > 0 and rel[left] < half_lvl:
                left -= 1
            right = deep_idx
            while right < len(rel) - 1 and rel[right] < half_lvl:
                right += 1
            dur_in = max(1, deep_idx - left)
            dur_out = max(1, right - deep_idx)
            asym = abs(dur_out - dur_in) / (dur_in + dur_out)

            # Periodicity veto with BoxLeastSquares over the searched points of the sector that holds the best box.
            # One sector is one continuous month; this period grid cannot keep phase across the years between sectors.
            t_all = np.concatenate([p[1] for p in parts if p[0] == best["sector"]])
            rel_all = np.concatenate([p[2] for p in parts if p[0] == best["sector"]])
            bls = BoxLeastSquares(t_all, rel_all)
            periods = np.linspace(0.7, 9.0, 250)
            bls_res = bls.power(periods, 0.1)
            bls_snr = (float(np.max(bls_res.power)) - float(np.median(bls_res.power))) / max(1e-6, float(np.std(bls_res.power)))

            # Independent extraction check (only when a second extraction on the same cadences is loaded)
            if len(f_cut) == len(t):
                cut_seg = f_cut[best["lo"] : best["lo"] + len(rel)]
                d_cut = max(0.0, float(1.0 - np.median(cut_seg[a:b])))
                pipe_agree = 1.0 if (dip_depth < 0.02 or d_cut >= 0.4 * dip_depth) else 0.0
            else:
                pipe_agree = 1.0

            vals.update({
                "dip_snr": max(0.0, float(best["snr"])),
                "dip_depth": float(dip_depth),
                "asymmetry": float(asym),
                "duration_hours": float(best["dur_h"]),
                "periodic_bls_snr": float(bls_snr),
                "pipeline_agree": float(pipe_agree),
                "best_time_btjd": float(best["t"][best["idx"]]),
                "best_sector": float(best["sector"]),
            })

        # Long-term secular fading in ZTF (independent of TESS coverage)
        if len(ztf_f) >= 20 and len(ztf_f) == len(ztf_t):
            zt_yr = (ztf_t - np.mean(ztf_t)) / 365.25
            slope = float(np.polyfit(zt_yr, ztf_f, 1)[0])
            slope_err = float(np.std(ztf_f)) / max(1e-4, math.sqrt(float(np.sum(zt_yr**2))))
            if slope < 0:
                vals["secular_fade_snr"] = abs(slope) / max(1e-4, slope_err)

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=vals,
            uncertainties={"dip_depth": 0.005},
            epochs=["TESS/Kepler", "ZTF/ASAS-SN"],
            caveats=["TESS 21 arcsec pixels require centroid/neighbour vetting."],
            extra={"search_status": tess_search_status(star, self.cfg)},
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        dip_hit = (
            v["dip_snr"] >= float(self.ch_cfg["dip_snr_threshold"])
            and v["dip_depth"] >= float(self.ch_cfg["min_dip_depth_frac"])
            and v["duration_hours"] >= float(self.ch_cfg["min_duration_hours"])
            and v["periodic_bls_snr"] < float(self.ch_cfg["periodic_bls_veto_snr"])
        )
        fade_hit = v["secular_fade_snr"] >= float(self.ch_cfg["secular_fade_snr_threshold"])
        flagged = bool(dip_hit or fade_hit)

        raw_stat = float(max(v["dip_snr"] if v["periodic_bls_snr"] < float(self.ch_cfg["periodic_bls_veto_snr"]) else 0.0, v["secular_fade_snr"]))
        from scipy.stats import norm

        fap = float(max(1e-12, norm.sf(max(0.0, raw_stat))))
        summary = f"Optical dip depth={v['dip_depth']:.3f}, SNR={v['dip_snr']:.2f}, asym={v['asymmetry']:.2f}, BLS_SNR={v['periodic_bls_snr']:.2f}"
        status = str(m.extra.get("search_status", "searched"))
        if status != "searched" and not flagged:
            # No TESS segment was searched: say so, so a zero is never read as a quiet star
            summary = f"NOT SEARCHED ({status}): no TESS light-curve segment to search"
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=flagged,
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=summary,
            feature_contributions={
                "dip_snr": float(v["dip_snr"]),
                "dip_depth": float(v["dip_depth"]),
                "asymmetry": float(v["asymmetry"]),
                "secular_fade_snr": float(v["secular_fade_snr"]),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "optical_dips_measurement.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        path.write_text(txt, encoding="utf-8")
        return [
            EvidenceRef(
                kind="lightcurve_measurement",
                path=str(path),
                description="Detrended dip matched-filter SNR, BLS periodicity check, and secular fade slope",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        ]
