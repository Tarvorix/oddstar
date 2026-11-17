"""Unified vetting orchestrator mapping channels/features to check suites."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.schemas import ChannelScore, VettingResult
from oddstar.vet.ir_checks import check_xmatch_integrity, run_ir_vetting
from oddstar.vet.lc_checks import check_known_variable, run_lc_vetting
from oddstar.vet.radio_checks import run_radio_vetting
from oddstar.vet.spectral_checks import run_spectral_vetting
from oddstar.vet.underluminous_checks import run_underluminous_vetting
from oddstar.vet.vanishing_checks import run_vanishing_vetting


def vet_target(
    star_or_cadence: dict[str, Any],
    scores: list[ChannelScore],
    cfg: PipelineConfig,
    out_dir: Path,
) -> list[VettingResult]:
    """Run every applicable vetting check for the channels that flagged or evaluated the target."""
    out_dir.mkdir(parents=True, exist_ok=True)
    flagged_channels = {s.channel for s in scores if s.flagged}
    evaluated_channels = {s.channel for s in scores}
    target_channels = flagged_channels or evaluated_channels

    results_by_name: dict[str, VettingResult] = {}

    def _merge(items: list[VettingResult]) -> None:
        for item in items:
            # Preserve fail > not_evaluated > pass if multiple channels invoke same check
            prev = results_by_name.get(item.check)
            if prev is None:
                results_by_name[item.check] = item
            elif item.status == "fail":
                results_by_name[item.check] = item
            elif item.status == "not_evaluated" and prev.status == "pass":
                results_by_name[item.check] = item

    if "ir_excess" in target_channels:
        _merge(run_ir_vetting(star_or_cadence, cfg, out_dir))
    if "ir_timedomain" in target_channels or "optical_dips" in target_channels:
        _merge(run_lc_vetting(star_or_cadence, cfg, out_dir))
    if "radio_narrowband" in target_channels:
        _merge(run_radio_vetting(star_or_cadence, cfg, out_dir))
    if "underluminous" in target_channels:
        _merge(run_underluminous_vetting(star_or_cadence, cfg, out_dir))
    if "spectral_lines" in target_channels:
        _merge(run_spectral_vetting(star_or_cadence, cfg, out_dir))
    if "vanishing" in target_channels:
        _merge(run_vanishing_vetting(star_or_cadence, cfg, out_dir))
    if "generic_outlier" in target_channels:
        # Section 6: generic_outlier runs xmatch_integrity, known_variable, and feature-linked checks
        _merge([check_xmatch_integrity(star_or_cadence, cfg, out_dir), check_known_variable(star_or_cadence, cfg, out_dir)])
        outlier_scores = [s for s in scores if s.channel == "generic_outlier"]
        drivers = set()
        for os_score in outlier_scores:
            for feat, contrib in os_score.feature_contributions.items():
                if abs(contrib) >= 2.0:
                    drivers.add(feat)
        if any(k in drivers for k in ("chi_w4", "w4_ratio", "chi_w3")):
            _merge(run_ir_vetting(star_or_cadence, cfg, out_dir))
        if any(k in drivers for k in ("dip_snr", "ir_var_chi2")):
            _merge(run_lc_vetting(star_or_cadence, cfg, out_dir))
        if "delta_m_sigma" in drivers:
            _merge(run_underluminous_vetting(star_or_cadence, cfg, out_dir))
        if "spec_line_snr" in drivers:
            _merge(run_spectral_vetting(star_or_cadence, cfg, out_dir))
        if "vanished_margin" in drivers:
            _merge(run_vanishing_vetting(star_or_cadence, cfg, out_dir))

    return list(results_by_name.values())
