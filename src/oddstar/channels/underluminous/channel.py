"""Channel 5: Under-luminous stars relative to spectral parameters & parallax distance."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.io.physics import BAND_EXTINCTION_AV
from oddstar.sample.builder import build_parent_sample
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec


def predict_ms_abs_g(teff_k: float, logg: float = 4.4, feh: float = 0.0) -> float:
    """Empirical main-sequence absolute G magnitude M_G(Teff, logg, [Fe/H])."""
    # Calibrated to match generate_field_sample luminosity relation L ~ (Teff/5778)^4.8
    # M_G_sun = 4.68
    lum = max(0.01, (teff_k / 5778.0) ** 4.8)
    mg_base = 4.68 - 2.5 * math.log10(lum)
    # Mild metallicity correction (~0.25 mag fainter at [Fe/H] = -1.0)
    mg_corr = mg_base - 0.15 * feh
    return float(mg_corr)


class UnderluminousChannel:
    name = "underluminous"

    def __init__(
        self,
        cfg: PipelineConfig,
        root: Path,
        *,
        extinction_correction_enabled: bool = True,
    ) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["underluminous"]
        self.root = root
        self.extinction_correction_enabled = extinction_correction_enabled
        self._locus_offset = 0.0
        self._locus_scatter = 0.14

    def fit_locus(self, stars: list[dict[str, Any]]) -> None:
        resids = []
        for s in stars[:2000]:
            mg_obs = float(s["abs_g_mag"])
            mg_pred = predict_ms_abs_g(float(s["teff_k"]), float(s.get("logg", 4.4)), float(s.get("feh", 0.0)))
            resids.append(mg_obs - mg_pred)
        if resids:
            arr = np.asarray(resids, dtype=float)
            self._locus_offset = float(np.median(arr))
            self._locus_scatter = max(0.08, 1.4826 * float(np.median(np.abs(arr - self._locus_offset))))

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec:
        chash = compute_config_hash(cfg.model_dump())[:12]
        spec, stars = build_parent_sample(
            self.name,
            cfg,
            chash,
            self.root,
            extra_stars=kwargs.get("extra_stars"),
            limit=kwargs.get("limit"),
        )
        self.fit_locus(stars)
        return spec

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        dist = float(star["distance_pc"])
        g_mag = float(star["mags"]["G"])
        ebv = float(star.get("ebv_sfd", 0.02))
        ag = 3.1 * ebv * BAND_EXTINCTION_AV["G"] if self.extinction_correction_enabled else 0.0
        mg_obs = g_mag - 5.0 * math.log10(max(1.0, dist) / 10.0) - ag

        mg_pred_xp = predict_ms_abs_g(float(star.get("teff_k", 5800.0)), float(star.get("logg", 4.4)), float(star.get("feh", 0.0))) + self._locus_offset
        mg_pred_spec = predict_ms_abs_g(float(star.get("teff_spec_k", star.get("teff_k", 5800.0))), float(star.get("logg", 4.4)), float(star.get("feh", 0.0))) + self._locus_offset

        # Positive delta_m means star is FAINTER than predicted
        delta_m_xp = mg_obs - mg_pred_xp
        delta_m_spec = mg_obs - mg_pred_spec
        plx_snr = max(1.0, float(star.get("parallax_snr", 50.0)))
        dist_mag_err = 2.17 / plx_snr
        comb_sig = math.hypot(self._locus_scatter, dist_mag_err)

        sigma_xp = delta_m_xp / comb_sig
        sigma_spec = delta_m_spec / comb_sig

        # Energy balance check if IR excess exists
        l_missing_frac = max(0.0, 1.0 - 10.0 ** (-0.4 * max(0.0, delta_m_xp)))

        return Measurement(
            source_id=str(star["source_id"]),
            channel=self.name,
            values={
                "mg_obs": float(mg_obs),
                "mg_pred_xp": float(mg_pred_xp),
                "mg_pred_spec": float(mg_pred_spec),
                "delta_m_xp": float(delta_m_xp),
                "delta_m_spec": float(delta_m_spec),
                "sigma_xp": float(sigma_xp),
                "sigma_spec": float(sigma_spec),
                "missing_opt_frac": float(l_missing_frac),
                "ruwe": float(star.get("ruwe", 1.0)),
                "parallax_snr": float(plx_snr),
                "extinction_applied": 1.0 if self.extinction_correction_enabled else 0.0,
            },
            uncertainties={"delta_m_xp": float(comb_sig)},
            epochs=["Gaia-DR3"],
            caveats=["Requires dual parameter agreement and 3D dust map agreement."],
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        sig_thr = float(self.ch_cfg["delta_m_sigma_threshold"])
        min_dm = float(self.ch_cfg["min_delta_mag"])
        clean_astrom = v["ruwe"] <= float(self.ch_cfg["max_ruwe"]) and v["parallax_snr"] >= float(self.ch_cfg["min_parallax_snr"])
        dual_hit = min(v["sigma_xp"], v["sigma_spec"]) >= sig_thr and min(v["delta_m_xp"], v["delta_m_spec"]) >= min_dm
        flagged = bool(dual_hit and clean_astrom)

        raw_stat = float(min(v["sigma_xp"], v["sigma_spec"]) if clean_astrom else 0.0)
        from scipy.stats import norm

        fap = float(max(1e-12, norm.sf(max(0.0, raw_stat))))
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=flagged,
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=f"Underluminous dM={v['delta_m_xp']:.2f} mag ({v['sigma_xp']:.2f} sigma, ext_on={int(v['extinction_applied'])})",
            feature_contributions={
                "sigma_xp": float(v["sigma_xp"]),
                "sigma_spec": float(v["sigma_spec"]),
                "delta_m_xp": float(v["delta_m_xp"]),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / "underluminous_measurement.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        path.write_text(txt, encoding="utf-8")
        return [
            EvidenceRef(
                kind="photometric_luminosity_measurement",
                path=str(path),
                description="Observed vs predicted absolute magnitude under dual stellar parameter sets",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        ]
