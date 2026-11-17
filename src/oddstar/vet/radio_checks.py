"""Radio narrowband vetting checks (Section 6)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result


def check_off_scan(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    hits_off = int(target.get("off_scan_hits", 0))
    ok = hits_off == 0
    return make_vet_result(
        out_dir,
        check="off_scan",
        status="pass" if ok else "fail",
        note=f"Off-target scan detections at matching frequency/drift: {hits_off}",
        metrics={"off_scan_hits": hits_off},
    )


def check_zero_drift(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    drift = float(target.get("drift_hz_s", 0.0))
    min_abs = float(cfg.channels["radio_narrowband"].get("min_drift_abs_hz_s", 0.05))
    ok = abs(drift) >= min_abs
    return make_vet_result(
        out_dir,
        check="zero_drift",
        status="pass" if ok else "fail",
        note=f"Measured drift rate={drift:.4f} Hz/s (min |drift|={min_abs} Hz/s)",
        metrics={"drift_hz_s": drift, "min_abs_drift_hz_s": min_abs},
    )


def check_rfi_band(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    in_rfi = bool(target.get("in_rfi_band", False))
    return make_vet_result(
        out_dir,
        check="rfi_band",
        status="fail" if in_rfi else "pass",
        note="Center frequency outside known satellite/terrestrial RFI allocations" if not in_rfi else "Signal falls inside known terrestrial/satellite RFI band",
        metrics={"in_rfi_band": in_rfi, "center_freq_mhz": target.get("center_freq_mhz")},
    )


def check_multi_target(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    multi = bool(target.get("multi_target_match", False))
    return make_vet_result(
        out_dir,
        check="multi_target",
        status="fail" if multi else "pass",
        note="Signal absent toward unrelated targets in session" if not multi else "Same frequency/drift signature observed toward unrelated targets in observing session",
        metrics={"multi_target_match": multi},
    )


def check_reobservation(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    reobs = bool(target.get("reobserved_confirmed", False))
    return make_vet_result(
        out_dir,
        check="reobservation",
        status="pass" if reobs else "not_evaluated",
        note="Confirmed in independent second-epoch radio cadence" if reobs else "No independent second-epoch reobservation confirmed (caps tier at C)",
        metrics={"reobserved_confirmed": reobs},
    )


def run_radio_vetting(target: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [
        check_off_scan,
        check_zero_drift,
        check_rfi_band,
        check_multi_target,
        check_reobservation,
    ]
    return [fn(target, cfg, out_dir) for fn in checks]
