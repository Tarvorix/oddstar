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
        lc = star.get("lightcurves", {})
        t = np.asarray(lc.get("tess_time_btjd", []), dtype=float)
        f_qlp = np.asarray(lc.get("tess_flux_qlp", []), dtype=float)
        f_cut = np.asarray(lc.get("tess_flux_tesscut", []), dtype=float)
        ztf_t = np.asarray(lc.get("ztf_time_mjd", []), dtype=float)
        ztf_f = np.asarray(lc.get("ztf_flux", []), dtype=float)

        if len(t) < 20:
            vals = {"dip_snr": 0.0, "dip_depth": 0.0, "asymmetry": 0.0, "duration_hours": 0.0, "periodic_bls_snr": 0.0, "secular_fade_snr": 0.0, "pipeline_agree": 0.0}
        else:
            dt = float(np.median(np.diff(t)))
            win_pts = max(15, int(round(float(self.ch_cfg["detrend_window_days"]) / max(1e-4, dt))) | 1)
            # Robust sliding upper-quantile/median baseline preserving multi-day dips
            baseline = median_filter(f_qlp, size=win_pts, mode="nearest")
            # Iterative upper clip so deep wide dips don't pull the baseline down
            resid_init = f_qlp - baseline
            noise_est = max(1e-4, 1.4826 * float(np.median(np.abs(resid_init))))
            masked = np.where(resid_init < -2.5 * noise_est, np.nan, f_qlp)
            valid = ~np.isnan(masked)
            if np.sum(valid) >= 10:
                interp_f = np.interp(t, t[valid], masked[valid])
                baseline = median_filter(interp_f, size=win_pts, mode="nearest")

            rel = f_qlp / np.maximum(1e-4, baseline)
            sigma = max(1e-4, 1.4826 * float(np.median(np.abs(rel - np.median(rel)))))

            # Matched-filter boxcar search over dip durations (6h to 4 days)
            best_snr = 0.0
            best_idx = 0
            best_dur_h = 0.0
            for dur_days in (0.25, 0.5, 1.0, 2.0, 3.5):
                w = max(2, int(round(dur_days / max(1e-4, dt))))
                kernel = np.ones(w) / w
                smoothed_deficit = np.convolve(1.0 - rel, kernel, mode="same")
                snr_curve = smoothed_deficit * math.sqrt(w) / sigma
                idx = int(np.argmax(snr_curve))
                if snr_curve[idx] > best_snr:
                    best_snr = float(snr_curve[idx])
                    best_idx = idx
                    best_dur_h = dur_days * 24.0

            dip_depth = max(0.0, float(1.0 - np.min(rel)))
            # Measure ingress vs egress asymmetry around deepest point
            deep_idx = int(np.argmin(rel))
            half_lvl = 1.0 - 0.5 * dip_depth
            left = deep_idx
            while left > 0 and rel[left] < half_lvl:
                left -= 1
            right = deep_idx
            while right < len(rel) - 1 and rel[right] < half_lvl:
                right += 1
            dur_in = max(1, deep_idx - left)
            dur_out = max(1, right - deep_idx)
            asym = abs(dur_out - dur_in) / (dur_in + dur_out)

            # Check periodicity with BoxLeastSquares
            bls = BoxLeastSquares(t, rel)
            periods = np.linspace(0.7, 9.0, 250)
            bls_res = bls.power(periods, 0.1)
            bls_snr = (float(np.max(bls_res.power)) - float(np.median(bls_res.power))) / max(1e-6, float(np.std(bls_res.power)))

            # Independent extraction check
            if len(f_cut) == len(t):
                d_cut = max(0.0, float(1.0 - np.min(f_cut)))
                pipe_agree = 1.0 if (dip_depth < 0.02 or d_cut >= 0.4 * dip_depth) else 0.0
            else:
                pipe_agree = 1.0

            # Long-term secular fading in ZTF
            secular_fade_snr = 0.0
            if len(ztf_f) >= 20:
                zt_yr = (ztf_t - np.mean(ztf_t)) / 365.25
                slope = float(np.polyfit(zt_yr, ztf_f, 1)[0])
                slope_err = float(np.std(ztf_f)) / max(1e-4, math.sqrt(float(np.sum(zt_yr**2))))
                if slope < 0:
                    secular_fade_snr = abs(slope) / max(1e-4, slope_err)

            vals = {
                "dip_snr": float(best_snr),
                "dip_depth": float(dip_depth),
                "asymmetry": float(asym),
                "duration_hours": float(best_dur_h),
                "periodic_bls_snr": float(bls_snr),
                "secular_fade_snr": float(secular_fade_snr),
                "pipeline_agree": float(pipe_agree),
            }

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=vals,
            uncertainties={"dip_depth": 0.005},
            epochs=["TESS/Kepler", "ZTF/ASAS-SN"],
            caveats=["TESS 21 arcsec pixels require centroid/neighbour vetting."],
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
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=flagged,
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=f"Optical dip depth={v['dip_depth']:.3f}, SNR={v['dip_snr']:.2f}, asym={v['asymmetry']:.2f}, BLS_SNR={v['periodic_bls_snr']:.2f}",
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
