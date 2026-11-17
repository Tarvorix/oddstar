"""Channel 8: Vanishing and appearing sources across 20th-century plates & modern surveys."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.io.physics import angular_separation_arcsec, propagate_position
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec


class VanishingChannel:
    name = "vanishing"

    def __init__(
        self,
        cfg: PipelineConfig,
        root: Path,
        *,
        proper_motion_correction_enabled: bool = True,
    ) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["vanishing"]
        self.root = root
        self.proper_motion_correction_enabled = proper_motion_correction_enabled

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
        vd = star.get("vanishing_data", {})
        usno_epoch = float(vd.get("usno_epoch", 1953.4))
        usno_ra = float(vd.get("usno_ra_deg", star["ra_deg"]))
        usno_dec = float(vd.get("usno_dec_deg", star["dec_deg"]))
        usno_mag = float(vd.get("usno_mag_r", star["mags"]["G"] + 0.2))
        n_plates = int(vd.get("early_plates_count", 2))
        ps1_present = bool(vd.get("present_in_panstarrs", True))
        gaia_present = bool(vd.get("present_in_gaia", True))
        ps1_depth = float(vd.get("panstarrs_depth_mag", 21.8))

        if self.proper_motion_correction_enabled:
            pred_ra, pred_dec = propagate_position(
                float(star["ra_deg"]),
                float(star["dec_deg"]),
                float(star.get("pmra_mas_yr", 0.0)),
                float(star.get("pmdec_mas_yr", 0.0)),
                2016.0,
                usno_epoch,
            )
        else:
            pred_ra, pred_dec = float(star["ra_deg"]), float(star["dec_deg"])

        sep_arcsec = angular_separation_arcsec(pred_ra, pred_dec, usno_ra, usno_dec)
        match_rad = float(self.ch_cfg["match_radius_arcsec"])
        modern_matched = (ps1_present or gaia_present) and (sep_arcsec <= match_rad)
        depth_margin = ps1_depth - usno_mag

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values={
                "early_plates_count": float(n_plates),
                "modern_matched": 1.0 if modern_matched else 0.0,
                "sep_arcsec": float(sep_arcsec),
                "depth_margin_mag": float(depth_margin),
                "usno_mag_r": float(usno_mag),
                "panstarrs_depth_mag": float(ps1_depth),
                "pm_correction_enabled": 1.0 if self.proper_motion_correction_enabled else 0.0,
            },
            uncertainties={"sep_arcsec": 0.4},
            epochs=[f"USNO-B1.0 ({usno_epoch:.1f})", "Pan-STARRS/Gaia"],
            caveats=["Historical plate mismatches are overwhelmingly artifacts/high-PM stars without deep multi-plate vetting."],
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        min_plates = float(self.ch_cfg["min_early_plates"])
        min_margin = float(self.ch_cfg["min_depth_margin_mag"])
        vanished = (
            v["early_plates_count"] >= min_plates
            and v["modern_matched"] == 0.0
            and v["depth_margin_mag"] >= min_margin
        )
        raw_stat = float(v["depth_margin_mag"] + (0.5 * v["sep_arcsec"] if v["modern_matched"] == 0.0 else 0.0)) if vanished else 0.0
        fap = float(max(1e-12, math.exp(-1.5 * raw_stat))) if raw_stat > 0 else 0.5
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=bool(vanished),
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=f"Vanishing matched={int(v['modern_matched'])}, sep={v['sep_arcsec']:.2f}\", margin={v['depth_margin_mag']:.1f} mag (PM_corr={int(v['pm_correction_enabled'])})",
            feature_contributions={
                "depth_margin_mag": float(v["depth_margin_mag"]),
                "sep_arcsec": float(v["sep_arcsec"]),
                "vanished_margin": float(raw_stat),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "vanishing_measurement.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        path.write_text(txt, encoding="utf-8")
        return [
            EvidenceRef(
                kind="epoch_astrometry_measurement",
                path=str(path),
                description="Proper-motion propagated epoch astrometry and local survey depth comparison",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        ]
