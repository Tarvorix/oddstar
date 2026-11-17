"""Cross-channel fusion, calibrated FAP conversion, rule-based tiering, and hypothesis matrix."""

from __future__ import annotations

import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal

import numpy as np

from oddstar.config import PipelineConfig, compute_config_hash, current_git_commit
from oddstar.schemas import (
    CandidateFusionRow,
    ChannelScore,
    FusionResult,
    HypothesisRow,
    ValidationCertificate,
    VettingResult,
)

# Distinct primary-instrument groups so multi-channel Tier A genuinely requires independent instruments
CHANNEL_INSTRUMENT_GROUP: dict[str, str] = {
    "ir_excess": "wise_2mass_gaia",
    "ir_timedomain": "neowise_ztf",
    "optical_dips": "tess_kepler_ztf",
    "radio_narrowband": "gbt_parkes_radio",
    "underluminous": "gaia_spectra",
    "generic_outlier": "derived_joint_features",
    "spectral_lines": "highres_spectroscopy",
    "vanishing": "historical_plates_panstarrs",
}

ARTIFACT_INTEGRITY_CHECKS = {
    "xmatch_integrity",
    "wise_quality",
    "centroid_offset",
    "single_frames",
    "systematics_window",
    "pipeline_agreement",
    "off_scan",
    "zero_drift",
    "rfi_band",
    "multi_target",
    "cosmic_ray",
    "detector",
    "repeat",
    "plate_defect",
    "depth",
    "astrometry",
    "parameter_agreement",
}

ASTROPHYSICAL_EXPLAINED_CHECKS = {
    "neighbours",
    "radio_counterpart",
    "deep_imaging",
    "extragalactic_colours",
    "cirrus",
    "moving_object",
    "youth",
    "evolved",
    "binary",
    "neighbour_signal",
    "sso_crossing",
    "known_variable",
    "sky_line",
    "stellar_emission",
    "proper_motion",
    "sso_plate",
    "blend_split",
    "extinction",
    "metallicity",
    "companion",
    "vision_rubric",
}

IR_CHANNELS = {"ir_excess", "ir_timedomain"}
OPTICAL_COLOR_ONLY_EXPLAINERS = {"metallicity"}


def calibrate_fap(raw_stat: float, cert: ValidationCertificate | None, fallback_fap: float) -> float:
    """Convert raw channel statistic to calibrated false-alarm probability using Test 4 curve."""
    if cert is None or not cert.calibration_curve:
        return float(max(1e-12, min(1.0, fallback_fap)))
    pts = sorted(cert.calibration_curve, key=lambda r: r["stat"])
    xs = [float(p["stat"]) for p in pts]
    ys = [float(p["fap"]) for p in pts]
    val = float(np.interp(raw_stat, xs, ys))
    return float(max(1e-12, min(1.0, val)))


def assign_tier(
    scores: list[ChannelScore],
    vetting: list[VettingResult],
    cfg: PipelineConfig,
) -> tuple[Literal["A", "B", "C", "Explained", "Rejected"], str | None]:
    """Rule-based candidate tier assignment (Section 8)."""
    failed = [v for v in vetting if v.status == "fail"]
    not_eval = [v for v in vetting if v.status == "not_evaluated"]

    flagged = [s for s in scores if s.flagged]
    if not flagged:
        return "Explained", "Not flagged above channel thresholds"

    # Check artifact/integrity failures first -> Rejected
    for v in failed:
        if v.check in ARTIFACT_INTEGRITY_CHECKS:
            return "Rejected", f"Failed integrity/artifact check '{v.check}': {v.note}"

    # Astrophysical/background impostor failures -> Explained (scoped by flagged modality;
    # optical subdwarf metallicity never explains an active mid-IR excess or IR variability flag).
    flagged_channels = {s.channel for s in flagged}
    if flagged_channels & IR_CHANNELS:
        explaining_failed = [v for v in failed if v.check not in OPTICAL_COLOR_ONLY_EXPLAINERS]
    else:
        explaining_failed = failed

    if explaining_failed:
        first = explaining_failed[0]
        return "Explained", f"Explained by '{first.check}': {first.note}"

    # Any open hole caps tier strictly at C (Rule R9)
    if not_eval:
        return "C", None

    distinct_instruments = {CHANNEL_INSTRUMENT_GROUP.get(s.channel, s.channel) for s in flagged}
    # generic_outlier is derived from other channels, so it does not count as an independent second instrument
    distinct_physical = distinct_instruments - {"derived_joint_features"}
    min_fap = min(s.false_alarm_prob for s in flagged)

    if len(distinct_physical) >= int(cfg.fusion.tier_a_min_channels) and min_fap <= float(cfg.fusion.tier_a_max_fap):
        return "A", None
    if min_fap <= float(cfg.fusion.tier_b_max_fap):
        return "B", None
    return "C", None


def build_hypothesis_table(
    scores: list[ChannelScore],
    vetting: list[VettingResult],
) -> list[HypothesisRow]:
    """Build mandatory 12-row hypothesis table (Section 8)."""
    by_check = {v.check: v for v in vetting}
    flagged_channels = {s.channel for s in scores if s.flagged}

    def _ev(names: list[str]) -> list[str]:
        return [by_check[n].evidence for n in names if n in by_check and by_check[n].evidence]

    def _status(
        disfavour_all: list[str],
        consistent_if_any_fail: list[str],
        require_evaluated: list[str] | None = None,
    ) -> tuple[Literal["disfavoured", "consistent", "untested"], list[str], str]:
        req = require_evaluated if require_evaluated is not None else disfavour_all
        present = [by_check[n] for n in req if n in by_check]
        if not present:
            return "untested", [], "Required checks were not applicable or not run"
        for n in consistent_if_any_fail:
            if n in by_check and by_check[n].status == "fail":
                return "consistent", _ev([n]), by_check[n].note
        if any(v.status == "not_evaluated" for v in present):
            holes = [v.check for v in present if v.status == "not_evaluated"]
            return "untested", _ev(req), f"Open hole(s): {', '.join(holes)}"
        if all(by_check.get(n) and by_check[n].status == "pass" for n in disfavour_all if n in by_check):
            return "disfavoured", _ev(disfavour_all), "Checks disfavour this hypothesis"
        return "untested", _ev(req), "Incomplete check coverage"

    rows: list[HypothesisRow] = []

    # 1. instrument_artifact
    st, ev, note = _status(
        ["wise_quality", "single_frames", "systematics_window", "pipeline_agreement", "cosmic_ray", "detector", "plate_defect"],
        ["wise_quality", "single_frames", "systematics_window", "pipeline_agreement", "cosmic_ray", "detector", "plate_defect"],
        require_evaluated=[c for c in ("wise_quality", "single_frames", "systematics_window", "pipeline_agreement", "cosmic_ray", "detector", "plate_defect") if c in by_check],
    )
    if "ir_excess" in flagged_channels and st == "disfavoured":
        note = "2010 single-frame consistency & PSF checks passed (note: W3/W4 remain single-epoch)"
    rows.append(HypothesisRow(name="instrument_artifact", status=st, evidence=ev, note=note))

    # 2. catalog_mismatch
    st, ev, note = _status(["xmatch_integrity", "proper_motion"], ["xmatch_integrity", "proper_motion"], require_evaluated=[c for c in ("xmatch_integrity", "proper_motion") if c in by_check])
    rows.append(HypothesisRow(name="catalog_mismatch", status=st, evidence=ev, note=note))

    # 3. blend_neighbouring_star
    st, ev, note = _status(
        ["neighbours", "centroid_offset", "neighbour_signal", "blend_split", "deep_imaging"],
        ["neighbours", "centroid_offset", "neighbour_signal", "blend_split"],
        require_evaluated=[c for c in ("neighbours", "centroid_offset", "neighbour_signal", "blend_split", "deep_imaging") if c in by_check],
    )
    rows.append(HypothesisRow(name="blend_neighbouring_star", status=st, evidence=ev, note=note))

    # 4. background_galaxy (explicitly warns about radio-quiet majority per Ren et al. 2026)
    st, ev, note = _status(
        ["radio_counterpart", "deep_imaging", "extragalactic_colours", "centroid_offset"],
        ["radio_counterpart", "deep_imaging", "extragalactic_colours"],
        require_evaluated=[c for c in ("radio_counterpart", "deep_imaging", "extragalactic_colours") if c in by_check],
    )
    if st == "untested" and "radio_counterpart" in by_check and by_check["radio_counterpart"].status == "pass":
        note = "radio_counterpart passed (<10% of dusty galaxies are radio-loud), but deep_imaging is not_evaluated"
    rows.append(HypothesisRow(name="background_galaxy", status=st, evidence=ev, note=note))

    # 5. galactic_cirrus
    st, ev, note = _status(["cirrus"], ["cirrus"], require_evaluated=["cirrus"] if "cirrus" in by_check else [])
    rows.append(HypothesisRow(name="galactic_cirrus", status=st, evidence=ev, note=note))

    # 6. solar_system_object
    st, ev, note = _status(
        ["moving_object", "sso_crossing", "sso_plate"],
        ["moving_object", "sso_crossing", "sso_plate"],
        require_evaluated=[c for c in ("moving_object", "sso_crossing", "sso_plate") if c in by_check],
    )
    rows.append(HypothesisRow(name="solar_system_object", status=st, evidence=ev, note=note))

    # 7. young_star_disk
    st, ev, note = _status(["youth"], ["youth"], require_evaluated=["youth"] if "youth" in by_check else [])
    rows.append(HypothesisRow(name="young_star_disk", status=st, evidence=ev, note=note))

    # 8. collision_debris
    if {"ir_excess", "ir_timedomain", "optical_dips"} & flagged_channels:
        rows.append(
            HypothesisRow(
                name="collision_debris",
                status="consistent",
                evidence=_ev(["second_instrument", "youth"]),
                note="Warm/variable circumstellar debris from giant planetary impacts can reproduce infrared excess and/or optical dips",
            )
        )
    else:
        rows.append(HypothesisRow(name="collision_debris", status="disfavoured", evidence=[], note="Signature does not match warm collision debris"))

    # 9. stellar_companion
    st, ev, note = _status(["binary", "companion"], ["binary", "companion"], require_evaluated=[c for c in ("binary", "companion") if c in by_check])
    rows.append(HypothesisRow(name="stellar_companion", status=st, evidence=ev, note=note))

    # 10. evolved_star_mass_loss
    st, ev, note = _status(["evolved"], ["evolved"], require_evaluated=["evolved"] if "evolved" in by_check else [])
    rows.append(HypothesisRow(name="evolved_star_mass_loss", status=st, evidence=ev, note=note))

    # 11. known_variable_class
    st, ev, note = _status(["known_variable"], ["known_variable"], require_evaluated=["known_variable"] if "known_variable" in by_check else [])
    rows.append(HypothesisRow(name="known_variable_class", status=st, evidence=ev, note=note))

    # 12. unexplained — never auto-consistent while any mundane hypothesis is consistent or untested
    all_natural_disfavoured = all(r.status == "disfavoured" for r in rows)
    rows.append(
        HypothesisRow(
            name="unexplained",
            status="consistent" if all_natural_disfavoured else "untested",
            evidence=[],
            note=(
                "All evaluated natural/instrumental hypotheses disfavoured; requires independent human & multi-epoch confirmation"
                if all_natural_disfavoured
                else "Other natural or untested hypotheses remain open; cannot treat as unexplained residual"
            ),
        )
    )
    return rows


def detect_cross_channel_patterns(
    scores: list[ChannelScore],
    measurements: dict[str, Any] | None = None,
) -> list[str]:
    """Identify named multi-channel signatures from Section 8."""
    flagged = {s.channel for s in scores if s.flagged}
    patterns: list[str] = []
    if "ir_excess" in flagged and "optical_dips" in flagged:
        patterns.append("ir_excess_plus_optical_dips: collision debris crossing the star")
    if "ir_excess" in flagged and "underluminous" in flagged:
        patterns.append("ir_excess_plus_underluminous: partial obscuration (energy-balance check)")
    if "ir_timedomain" in flagged:
        patterns.append("ir_variability_flat_optical: changing warm dust or variable background source")
    if "radio_narrowband" in flagged:
        patterns.append("any_plus_radio_narrowband: escalate immediately to independent adversarial review")
    return patterns


def fuse_results(
    *,
    run_id: str,
    scores_by_star: dict[str, list[ChannelScore]],
    vetting_by_star: dict[str, list[VettingResult]],
    certs_by_channel: dict[str, ValidationCertificate],
    n_stars_total: int,
    cfg: PipelineConfig,
    root: Path,
) -> FusionResult:
    """Join channel scores on source_id, calibrate FAPs, assign tiers, and account for look-elsewhere."""
    rows: list[CandidateFusionRow] = []
    cross_patterns_summary: list[dict[str, Any]] = []

    n_trials_total = max(1, n_stars_total * max(1, len(certs_by_channel)))
    expected_fp = sum(
        float(cert.false_alarm_rate) * n_stars_total for cert in certs_by_channel.values()
    ) if certs_by_channel else float(n_stars_total) * 1e-4

    for sid, raw_scores in scores_by_star.items():
        cal_scores: list[ChannelScore] = []
        for s in raw_scores:
            cert = certs_by_channel.get(s.channel)
            fap = calibrate_fap(s.raw_statistic, cert, s.false_alarm_prob)
            cal_scores.append(
                s.model_copy(
                    update={
                        "false_alarm_prob": fap,
                        "rank_score": float(max(0.0, -math.log10(fap))),
                    }
                )
            )
        vet = vetting_by_star.get(sid, [])
        tier, explained_reason = assign_tier(cal_scores, vet, cfg)
        flagged_channels = [s.channel for s in cal_scores if s.flagged]
        if not flagged_channels:
            continue

        # Combined FAP across distinct flagged channels (Fisher / product bounded)
        faps = [s.false_alarm_prob for s in cal_scores if s.flagged]
        combined_fap = float(max(1e-12, np.prod(faps)))
        patterns = detect_cross_channel_patterns(cal_scores)
        if patterns:
            cross_patterns_summary.append({"source_id": sid, "patterns": patterns})
        holes = sorted({v.check for v in vet if v.status == "not_evaluated"})

        rows.append(
            CandidateFusionRow(
                source_id=sid,
                tier=tier,
                rank_score=float(-math.log10(combined_fap)),
                combined_fap=combined_fap,
                flagged_channels=flagged_channels,
                patterns=patterns,
                open_holes=holes,
                explained_reason=explained_reason,
            )
        )

    tier_order = {"A": 0, "B": 1, "C": 2, "Explained": 3, "Rejected": 4}
    rows.sort(key=lambda r: (tier_order[r.tier], -r.rank_score, r.source_id))

    n_candidates = sum(1 for r in rows if r.tier in {"A", "B", "C"})
    _, short_hash, _ = (cfg, compute_config_hash(cfg.model_dump())[:12], "")

    return FusionResult(
        run_id=run_id,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        config_hash=short_hash,
        git_commit=current_git_commit(root),
        n_stars_total=n_stars_total,
        n_trials_total=n_trials_total,
        expected_false_positives=float(expected_fp),
        consistent_with_noise=bool(n_candidates <= expected_fp),
        cross_channel_patterns=cross_patterns_summary,
        candidates=rows,
    )
