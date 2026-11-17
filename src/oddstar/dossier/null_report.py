"""Null-result report generator with completeness contours and Clopper-Pearson upper limits."""

from __future__ import annotations

import math
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scipy.stats import beta

from oddstar.config import PipelineConfig, compute_config_hash, current_git_commit
from oddstar.schemas import NullReport, SampleCutStep, ValidationCertificate


def clopper_pearson_upper_limit(k_detections: int, n_trials: int, confidence: float = 0.95) -> float:
    """Exact one-sided Clopper-Pearson upper limit on occurrence fraction."""
    if n_trials <= 0:
        return 1.0
    if k_detections >= n_trials:
        return 1.0
    return float(beta.ppf(confidence, k_detections + 1, n_trials - k_detections))


def write_null_report(
    *,
    channel: str,
    run_id: str,
    sample_stars: list[dict[str, Any]],
    sample_cuts: list[SampleCutStep],
    healpix_cells: list[int],
    n_flagged: int,
    tier_breakdown: dict[str, int],
    cert: ValidationCertificate,
    manifests: list[str],
    cfg: PipelineConfig,
    root: Path,
    not_searched: dict[str, list[str]] | None = None,
) -> tuple[NullReport, Path, Path]:
    """Emit JSON and Markdown null-result report for a channel run (Section 11).

    `not_searched` maps a reason (e.g. "no_coverage", "failed", "not_fetched", "insufficient_visits") to the source
    IDs of sample stars the channel could not search. They are excluded from `n_searched`, the distributions and the
    upper-limit denominators, and recorded as explicit sample cuts. For `ir_timedomain` it is derived from the star
    records when not given, so a star without NEOWISE coverage or with a failed fetch is never counted as quiet.
    """
    if not_searched is None and channel == "ir_timedomain":
        from oddstar.channels.ir_timedomain.channel import search_coverage_summary

        not_searched = search_coverage_summary(sample_stars, cfg)["not_searched_ids"]
    excluded: set[str] = set()
    extra_cuts: list[SampleCutStep] = []
    n_remaining = len(sample_stars)
    sample_ids = {str(s.get("source_id")) for s in sample_stars}
    unknown = sorted({str(x) for v in (not_searched or {}).values() for x in v} - sample_ids)
    if unknown:
        # Never let a not-searched star silently fall back into the searched denominator (e.g. a regenerated sample)
        raise ValueError(f"not_searched lists {len(unknown)} source_id(s) absent from sample_stars: {', '.join(unknown[:5])}")
    for reason in sorted(not_searched or {}):
        ids = [i for i in dict.fromkeys(str(x) for x in not_searched[reason]) if i not in excluded]
        if not ids:
            continue
        excluded.update(ids)
        extra_cuts.append(
            SampleCutStep(
                cut_name=f"not_searched_{reason}",
                description=f"Excluded from the searched denominator: {reason} ({len(ids)} stars)",
                count_before=n_remaining,
                count_after=n_remaining - len(ids),
            )
        )
        n_remaining -= len(ids)
    searched_stars = [s for s in sample_stars if str(s.get("source_id")) not in excluded]
    sample_cuts = [*sample_cuts, *extra_cuts]
    n_searched = len(searched_stars)
    dist_bins: Counter[str] = Counter()
    spec_bins: Counter[str] = Counter()
    for s in searched_stars:
        d = float(s.get("distance_pc") or 100.0)
        if d < 50:
            dist_bins["0-50pc"] += 1
        elif d < 100:
            dist_bins["50-100pc"] += 1
        elif d < 200:
            dist_bins["100-200pc"] += 1
        else:
            dist_bins["200-300pc"] += 1

        teff = float(s.get("teff_k") or 5800.0)
        if teff >= 7300:
            spec_bins["A_or_hotter"] += 1
        elif teff >= 6000:
            spec_bins["F"] += 1
        elif teff >= 5200:
            spec_bins["G"] += 1
        elif teff >= 3900:
            spec_bins["K"] += 1
        else:
            spec_bins["M"] += 1

    contour_50 = [row for row in cert.completeness_grid if float(row.get("completeness", 0.0)) >= 0.50]
    contour_90 = [row for row in cert.completeness_grid if float(row.get("completeness", 0.0)) >= 0.90]

    # Compute effective upper limits at 50% and 90% sensitivity regimes
    n_survivors = sum(tier_breakdown.get(t, 0) for t in ("A", "B", "C"))
    upper_limits = []
    for level, cells in ((0.50, contour_50), (0.90, contour_90)):
        eff_n = int(round(n_searched * level))
        ul_95 = clopper_pearson_upper_limit(n_survivors, max(1, eff_n), confidence=0.95)
        upper_limits.append(
            {
                "sensitivity_level": level,
                "effective_n_stars": eff_n,
                "surviving_candidates": n_survivors,
                "upper_limit_fraction_95pct": ul_95,
                "n_grid_cells": len(cells),
            }
        )

    short_hash = compute_config_hash(cfg.model_dump())[:12]
    report = NullReport(
        channel=channel,
        run_id=run_id,
        gaia_release=cfg.gaia.release,
        config_hash=short_hash,
        git_commit=current_git_commit(root),
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        sample_cuts=sample_cuts,
        sky_coverage_healpix=healpix_cells,
        distance_distribution=dict(dist_bins),
        spectral_type_distribution=dict(spec_bins),
        completeness_grid=cert.completeness_grid,
        contour_50=contour_50,
        contour_90=contour_90,
        n_searched=n_searched,
        n_flagged=n_flagged,
        expected_false_alarms=float(cert.false_alarm_rate * n_searched),
        tier_breakdown=tier_breakdown,
        upper_limits=upper_limits,
        manifests=manifests,
        reproduce=f"oddstar search {channel} --run-id {run_id}",
    )

    out_dir = root / cfg.storage.runs_dir / run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / f"null_report_{channel}.json"
    md_path = out_dir / f"null_report_{channel}.md"
    json_path.write_text(report.model_dump_json(indent=2) + "\n", encoding="utf-8")

    md_lines = [
        f"# Null-Result & Sensitivity Report: `{channel}` (`{run_id}`)",
        "",
        f"- **Gaia Release**: `{cfg.gaia.release}` | **Config Hash**: `{short_hash}` | **Commit**: `{report.git_commit[:10]}`",
        f"- **Stars Searched**: {n_searched} across {len(healpix_cells)} HEALPix cells (`nside={cfg.storage.healpix_nside}`)",
        f"- **Not Searched**: {len(excluded)}"
        + (" (" + ", ".join(f"{c.cut_name.removeprefix('not_searched_')}={c.count_before - c.count_after}" for c in extra_cuts) + ")" if extra_cuts else "")
        + " (excluded from every denominator below)",
        f"- **Flagged**: {n_flagged} (expected false alarms: {report.expected_false_alarms:.2f})",
        f"- **Tier Breakdown**: {tier_breakdown}",
        "",
        "## Upper Limits (95% Confidence)",
        "",
        "| Sensitivity | Effective N | Surviving Candidates | 95% Upper Limit Fraction |",
        "| --- | --- | --- | --- |",
    ]
    for ul in upper_limits:
        md_lines.append(
            f"| {ul['sensitivity_level']:.0%} | {ul['effective_n_stars']} | {ul['surviving_candidates']} | {ul['upper_limit_fraction_95pct']:.5f} |"
        )
    md_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
    return report, json_path, md_path
