"""Channel 2: Mid-infrared time-domain variability (NEOWISE-R + ZTF/ASAS-SN)."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec


class IrTimedomainChannel:
    name = "ir_timedomain"

    def __init__(self, cfg: PipelineConfig, root: Path) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["ir_timedomain"]
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
        neo = star.get("neowise", {})
        lc = star.get("lightcurves", {})
        mjd = np.asarray(neo.get("mjd", []), dtype=float)
        w1 = np.asarray(neo.get("w1_mjy", []), dtype=float)
        w2 = np.asarray(neo.get("w2_mjy", []), dtype=float)
        e1 = np.maximum(1e-3, np.asarray(neo.get("w1_err_mjy", np.ones_like(w1)), dtype=float))
        e2 = np.maximum(1e-3, np.asarray(neo.get("w2_err_mjy", np.ones_like(w2)), dtype=float))

        if len(mjd) < 3:
            values = {"chi2_w1": 0.0, "chi2_w2": 0.0, "slope_snr": 0.0, "max_jump_snr": 0.0, "color_change_mag": 0.0, "optical_chi2": 1.0, "n_visits": float(len(mjd))}
        else:
            med1, med2 = float(np.median(w1)), float(np.median(w2))
            chi2_w1 = float(np.mean(((w1 - med1) / e1) ** 2))
            chi2_w2 = float(np.mean(((w2 - med2) / e2) ** 2))

            # Linear trend SNR across epochs (W2 most sensitive to warm dust changes)
            t_yr = (mjd - np.mean(mjd)) / 365.25
            denom = max(1e-6, float(np.sum((t_yr**2) / (e2**2))))
            slope = float(np.sum((t_yr * (w2 - np.mean(w2))) / (e2**2))) / denom
            slope_err = 1.0 / math.sqrt(denom)
            slope_snr = abs(slope) / max(1e-6, slope_err)

            # Max single-visit jump SNR
            jumps = np.abs(np.diff(w2)) / np.sqrt(e2[1:] ** 2 + e2[:-1] ** 2)
            max_jump_snr = float(np.max(jumps)) if len(jumps) else 0.0

            # W1 - W2 color excursion (magnitude scale)
            col = -2.5 * np.log10(np.maximum(1e-4, w1) / np.maximum(1e-4, w2))
            color_change_mag = float(np.max(col) - np.min(col))

            # Optical quietness check over overlapping baseline
            ztf = np.asarray(lc.get("ztf_flux", []), dtype=float)
            if len(ztf) > 5:
                z_med = float(np.median(ztf))
                # Exclude late optical eclipses that follow an earlier IR brightening (ASASSN-21qj behavior)
                early_ztf = ztf[: max(5, len(ztf) // 2)]
                opt_chi2 = float(np.mean(((early_ztf - np.median(early_ztf)) / 0.01) ** 2))
                opt_full_chi2 = float(np.mean(((ztf - z_med) / 0.01) ** 2))
                # Use early-epoch optical quietness if an IR burst preceded optical eclipse
                optical_chi2 = min(opt_chi2, opt_full_chi2)
            else:
                optical_chi2 = 1.0

            values = {
                "chi2_w1": chi2_w1,
                "chi2_w2": chi2_w2,
                "slope_snr": slope_snr,
                "max_jump_snr": max_jump_snr,
                "color_change_mag": color_change_mag,
                "optical_chi2": optical_chi2,
                "n_visits": float(len(mjd)),
            }

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values=values,
            uncertainties={"chi2_w2": 0.2},
            epochs=["2014-2024"],
            caveats=["NEOWISE-R covers W1 (3.4um) and W2 (4.6um) only; W3/W4 have no post-2010 baseline."],
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        ir_stat = max(v["chi2_w1"], v["chi2_w2"], v["slope_snr"], v["max_jump_snr"])
        ir_active = (
            v["chi2_w1"] >= float(self.ch_cfg["reduced_chi2_threshold"])
            or v["chi2_w2"] >= float(self.ch_cfg["reduced_chi2_threshold"])
            or v["slope_snr"] >= float(self.ch_cfg["slope_snr_threshold"])
            or v["max_jump_snr"] >= float(self.ch_cfg["jump_snr_threshold"])
        )
        opt_quiet = v["optical_chi2"] <= float(self.ch_cfg["optical_flat_max_chi2"])
        enough_visits = v["n_visits"] >= float(self.ch_cfg["min_visits"])
        flagged = bool(ir_active and opt_quiet and enough_visits)

        fap = float(max(1e-12, math.exp(-0.5 * max(0.0, ir_stat))))
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=float(ir_stat),
            false_alarm_prob=fap,
            flagged=flagged,
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=f"IR time-domain stat={ir_stat:.2f} (chi2_w2={v['chi2_w2']:.2f}, optical_chi2={v['optical_chi2']:.2f})",
            feature_contributions={
                "chi2_w2": float(v["chi2_w2"]),
                "slope_snr": float(v["slope_snr"]),
                "max_jump_snr": float(v["max_jump_snr"]),
                "optical_chi2": float(v["optical_chi2"]),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "ir_timedomain_measurement.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        path.write_text(txt, encoding="utf-8")
        return [
            EvidenceRef(
                kind="timeseries_measurement",
                path=str(path),
                description="NEOWISE W1/W2 visit statistics and overlapping optical variability metrics",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        ]
