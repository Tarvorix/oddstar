"""Channel 4: Narrowband Doppler-drifting radio search (Breakthrough Listen cadences)."""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig, compute_config_hash
from oddstar.dossier.plots import render_waterfall_plot
from oddstar.io.fixtures_data import RADIO_CADENCES
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleCutStep, SampleSpec


def dedoppler_search_scan(
    power: np.ndarray,
    dt_s: float,
    df_hz: float,
    max_abs_drift: float = 4.0,
) -> tuple[float, float, float]:
    """Run incoherent de-doppler drift search on a 2D (time, freq) filterbank block.

    Returns (max_snr, best_freq_offset_hz, best_drift_hz_s).
    """
    n_t, n_f = power.shape
    # Noise estimate from robust MAD of time-integrated spectrum at zero drift
    zero_int = np.mean(power, axis=0)
    med = float(np.median(zero_int))
    sig = max(1e-6, 1.4826 * float(np.median(np.abs(zero_int - med))))

    drifts = np.linspace(-max_abs_drift, max_abs_drift, 161)
    t_rel = np.arange(n_t, dtype=float) * dt_s
    best_snr = 0.0
    best_ch = n_f // 2
    best_drift = 0.0
    for d in drifts:
        shifts = np.round((d * t_rel) / df_hz).astype(int)
        acc = np.zeros(n_f, dtype=float)
        for ti, sh in enumerate(shifts):
            acc += np.roll(power[ti], -sh)
        acc /= n_t
        snr_spec = (acc - med) / sig
        idx = int(np.argmax(snr_spec))
        if snr_spec[idx] > best_snr:
            best_snr = float(snr_spec[idx])
            best_ch = idx
            best_drift = float(d)

    freq_hz = (best_ch - n_f / 2.0) * df_hz
    return best_snr, freq_hz, best_drift


def inject_synthetic_drifting_tone(
    cadence: dict[str, Any],
    *,
    f0_hz: float,
    drift_hz_s: float,
    snr: float,
) -> dict[str, Any]:
    """Setigen-equivalent synthetic drifting narrowband signal injection into real filterbank noise."""
    dt_s = float(cadence["dt_s"])
    df_hz = float(cadence["df_hz"])
    out = json.loads(json.dumps(cadence))
    for sc in out["scans"]:
        if sc["kind"] != "ON":
            continue
        arr = np.asarray(sc["power"], dtype=float)
        n_t, n_f = arr.shape
        t_start = float(sc["t_start_s"])
        for ti in range(n_t):
            f_hz = f0_hz + drift_hz_s * (t_start + ti * dt_s)
            ch = int(round(n_f / 2.0 + f_hz / df_hz))
            if 0 <= ch < n_f:
                arr[ti, ch] += snr * 0.08
        sc["power"] = arr.tolist()
    return out


class RadioNarrowbandChannel:
    name = "radio_narrowband"

    def __init__(self, cfg: PipelineConfig, root: Path) -> None:
        self.cfg = cfg
        self.ch_cfg = cfg.channels["radio_narrowband"]
        self.root = root

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec:
        chash = compute_config_hash(cfg.model_dump())[:12]
        ids = sorted(RADIO_CADENCES.keys())
        return SampleSpec(
            channel=self.name,
            gaia_release=cfg.gaia.release,
            config_hash=chash,
            seed=cfg.sample.seed,
            cuts=[SampleCutStep(cut_name="bl_cadences", description="ON/OFF cadences available", count_before=len(ids), count_after=len(ids))],
            total_count=len(ids),
            source_ids=ids,
        )

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement:
        cadence = star if "scans" in star else data.get("cadence")
        if cadence is None:
            return Measurement(
                source_id=str(star["source_id"]),
                channel=self.name,
                values={"max_on_snr": 0.0, "min_on_snr": 0.0, "off_scan_hits": 0.0, "drift_hz_s": 0.0, "on_persistence": 0.0},
                epochs=["breakthrough_listen"],
                caveats=["No Breakthrough Listen cadence available at target coordinates."],
            )

        dt_s = float(cadence["dt_s"])
        df_hz = float(cadence["df_hz"])
        max_drift = float(self.ch_cfg["max_abs_drift_hz_s"])
        min_snr = float(self.ch_cfg["min_snr"])
        tol_hz = float(self.ch_cfg["freq_match_tol_hz"])

        on_hits = []
        off_hits = []
        for sc in cadence["scans"]:
            arr = np.asarray(sc["power"], dtype=float)
            snr, f_hz, drift = dedoppler_search_scan(arr, dt_s, df_hz, max_abs_drift=max_drift)
            entry = {"scan_index": sc["scan_index"], "t_start_s": sc["t_start_s"], "snr": snr, "f_hz": f_hz, "drift_hz_s": drift}
            if sc["kind"] == "ON":
                on_hits.append(entry)
            else:
                off_hits.append(entry)

        # Check frequency/drift consistency across ON scans
        ref = max(on_hits, key=lambda x: x["snr"]) if on_hits else {"snr": 0.0, "f_hz": 0.0, "drift_hz_s": 0.0, "t_start_s": 0.0}
        on_persist = 0
        eff_tol = max(tol_hz, abs(ref["drift_hz_s"]) * 30.0)
        for h in on_hits:
            if h["snr"] >= min_snr:
                expected_f = ref["f_hz"] + ref["drift_hz_s"] * (h["t_start_s"] - ref["t_start_s"])
                if abs(h["f_hz"] - expected_f) <= eff_tol or abs(h["drift_hz_s"] - ref["drift_hz_s"]) <= 0.4:
                    on_persist += 1

        off_count = 0
        for h in off_hits:
            if h["snr"] >= min_snr:
                expected_f = ref["f_hz"] + ref["drift_hz_s"] * (h["t_start_s"] - ref["t_start_s"])
                if abs(h["f_hz"] - expected_f) <= eff_tol:
                    off_count += 1

        min_on_snr = min((h["snr"] for h in on_hits), default=0.0)
        max_on_snr = max((h["snr"] for h in on_hits), default=0.0)
        on_frac = on_persist / max(1, len(on_hits))

        # Populate vetting context fields directly on cadence dict if missing
        cadence["off_scan_hits"] = off_count
        cadence["drift_hz_s"] = ref["drift_hz_s"]

        return Measurement(
            source_id=str(cadence["source_id"]),
            channel=self.name,
            values={
                "max_on_snr": float(max_on_snr),
                "min_on_snr": float(min_on_snr),
                "off_scan_hits": float(off_count),
                "drift_hz_s": float(ref["drift_hz_s"]),
                "on_persistence": float(on_frac),
                "center_freq_mhz": float(cadence.get("center_freq_mhz", 1420.0)),
            },
            uncertainties={"drift_hz_s": 0.02},
            epochs=["BL-cadence"],
            caveats=["Requires second-epoch independent reobservation above Tier C."],
        )

    def score(self, m: Measurement) -> ChannelScore:
        v = m.values
        min_snr = float(self.ch_cfg["min_snr"])
        min_drift = float(self.ch_cfg["min_drift_abs_hz_s"])
        flagged = (
            v["min_on_snr"] >= min_snr
            and v["on_persistence"] >= 0.99
            and v["off_scan_hits"] == 0.0
            and abs(v["drift_hz_s"]) >= min_drift
        )
        raw_stat = float(v["min_on_snr"] if (v["off_scan_hits"] == 0.0 and v["on_persistence"] > 0.5) else 0.0)
        from scipy.stats import norm

        fap = float(max(1e-12, norm.sf(max(0.0, raw_stat))))
        return ChannelScore(
            source_id=m.source_id,
            channel=self.name,
            raw_statistic=raw_stat,
            false_alarm_prob=fap,
            flagged=bool(flagged),
            rank_score=float(max(0.0, -math.log10(fap))),
            summary=f"Radio min_ON_SNR={v['min_on_snr']:.1f}, OFF_hits={int(v['off_scan_hits'])}, drift={v['drift_hz_s']:.3f} Hz/s",
            feature_contributions={
                "min_on_snr": float(v["min_on_snr"]),
                "off_scan_hits": float(v["off_scan_hits"]),
                "drift_hz_s": float(v["drift_hz_s"]),
            },
        )

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]:
        out_dir.mkdir(parents=True, exist_ok=True)
        refs = []
        json_path = out_dir / "radio_hits.json"
        txt = json.dumps(m.model_dump(), indent=2, sort_keys=True)
        json_path.write_text(txt, encoding="utf-8")
        refs.append(
            EvidenceRef(
                kind="radio_hit_table",
                path=str(json_path),
                description="De-doppler ON/OFF cadence hit metrics",
                sha256=hashlib.sha256(txt.encode("utf-8")).hexdigest(),
            )
        )
        if "scans" in star:
            wf_path = render_waterfall_plot(star, out_dir / "radio_waterfall.png")
            refs.append(
                EvidenceRef(
                    kind="radio_waterfall",
                    path=str(wf_path),
                    description="Stacked ON/OFF filterbank waterfall cutout",
                    sha256=hashlib.sha256(wf_path.read_bytes()).hexdigest(),
                )
            )
        return refs
