"""Candidate dossier generator (dossier.json, dossier.md, and evidence plots)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Literal

from oddstar.config import PipelineConfig, compute_config_hash, current_git_commit
from oddstar.dossier.plots import (
    render_cutouts_panel,
    render_hr_diagram,
    render_lightcurves_plot,
    render_sed_plot,
    render_statistic_distribution,
    render_waterfall_plot,
)
from oddstar.fuse.engine import (
    assign_tier,
    build_hypothesis_table,
    calibrate_fap,
    detect_cross_channel_patterns,
)
from oddstar.schemas import (
    ChannelDossierEntry,
    ChannelScore,
    Coords,
    Dossier,
    DossierReviewSummary,
    Measurement,
    Provenance,
    ReviewStatus,
    StarMetadata,
    ValidationCertificate,
    VettingResult,
)


def _what_would_settle(
    vetting: list[VettingResult],
    scores: list[ChannelScore],
) -> list[str]:
    actions: list[str] = []
    not_eval = {v.check for v in vetting if v.status == "not_evaluated"}
    flagged = {s.channel for s in scores if s.flagged}

    if "second_instrument" in not_eval or "ir_excess" in flagged:
        actions.append("Second-epoch photometry near 22 um (e.g. JWST MIRI) to test post-2010 stability")
        actions.append("Mid-infrared spectrum (5-25 um) to distinguish silicate debris emission from blackbody dust or PAH lines")
    if "deep_imaging" in not_eval:
        actions.append("Sub-arcsecond deep optical/NIR imaging (HST/JWST/adaptive optics) inside the 12\" WISE beam to rule out radio-quiet background galaxies")
    if "reobservation" in not_eval or "radio_narrowband" in flagged:
        actions.append("Independent second-epoch radio cadence confirming signal persistence and drift")
    if "second_survey" in not_eval or "optical_dips" in flagged:
        actions.append("Multi-band high-cadence optical monitoring to test chromatic dust extinction vs grey occultation")
    if not actions:
        actions.append("No unresolved archive holes remain; independent follow-up observation only if required")
    return actions


def _build_four_lead_lines(
    source_id: str,
    tier: str,
    scores: list[ChannelScore],
    vetting: list[VettingResult],
    what_would_settle: list[str],
) -> list[str]:
    """Generate the mandatory four plain-word opening lines (Section 11)."""
    flagged = [s for s in scores if s.flagged]
    if flagged:
        meas_desc = "; ".join(f"[{s.channel}] {s.summary}" for s in flagged)
    else:
        meas_desc = "No search channel crossed detection threshold."

    passed = [v.check for v in vetting if v.status == "pass"]
    failed = [v.check for v in vetting if v.status == "fail"]
    not_eval = [v.check for v in vetting if v.status == "not_evaluated"]

    ruled_out = (
        f"Checks passed ({len(passed)}): {', '.join(passed[:6])}" + ("..." if len(passed) > 6 else "")
        if passed
        else "No vetting checks passed."
    )
    if failed:
        ruled_out += f" | Explained/rejected by: {', '.join(failed)}."

    not_ruled = (
        f"Open holes ({len(not_eval)} not evaluated): {', '.join(not_eval)}. Tier={tier}."
        if not_eval
        else f"Zero open archive holes (all applicable checks evaluated). Tier={tier}."
    )
    settle = f"Next observation needed: {what_would_settle[0]}."
    return [
        f"What was measured: {meas_desc}",
        f"What has been ruled out: {ruled_out}",
        f"What has not been ruled out: {not_ruled}",
        f"What observation would settle it: {settle}",
    ]


def write_dossier(
    *,
    source_id: str,
    star: dict[str, Any],
    measurements: list[Measurement],
    scores: list[ChannelScore],
    vetting: list[VettingResult],
    sample_stars: list[dict[str, Any]],
    sample_stats_by_channel: dict[str, list[float]],
    certs_by_channel: dict[str, ValidationCertificate],
    reviews: list[DossierReviewSummary],
    run_id: str,
    manifests: list[str],
    cfg: PipelineConfig,
    root: Path,
    review_status: ReviewStatus | None = None,
) -> tuple[Dossier, Path, Path]:
    """Generate dossier.json, dossier.md, and standard evidence plots under dossiers/<source_id>/."""
    clean_id = source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
    dossier_dir = root / cfg.storage.dossiers_dir / clean_id
    ev_dir = dossier_dir / "evidence"
    ev_dir.mkdir(parents=True, exist_ok=True)

    # Calibrate scores
    cal_scores: list[ChannelScore] = []
    for s in scores:
        cert = certs_by_channel.get(s.channel)
        fap = calibrate_fap(s.raw_statistic, cert, s.false_alarm_prob)
        cal_scores.append(s.model_copy(update={"false_alarm_prob": fap}))

    tier: Literal["A", "B", "C", "Explained", "Rejected"]
    tier, _ = assign_tier(cal_scores, vetting, cfg)
    # Apply Section 9 reviewer disagreement downgrade if any open finding exists
    if any(r.verdict in {"not_reproduced", "unresolved"} for r in reviews) and tier in {"A", "B"}:
        tier = "C"

    hypotheses = build_hypothesis_table(cal_scores, vetting)
    patterns = detect_cross_channel_patterns(cal_scores)
    settle = _what_would_settle(vetting, cal_scores)
    leads = _build_four_lead_lines(clean_id, tier, cal_scores, vetting, settle)

    # Render standard evidence plots
    meas_by_ch = {m.channel: m for m in measurements}
    generated_plots: dict[str, str] = {}

    if "scans" in star:
        wf = render_waterfall_plot(star, ev_dir / "radio_waterfall.png")
        generated_plots["radio_narrowband"] = str(wf.relative_to(dossier_dir))
    else:
        ir_vals = meas_by_ch["ir_excess"].values if "ir_excess" in meas_by_ch else {}
        sed_p = render_sed_plot(star, ir_vals, ev_dir / "sed.png")
        cut_p = render_cutouts_panel(star, ev_dir / "wise_cutouts.png")
        lc_p = render_lightcurves_plot(star, ev_dir / "lightcurves.png")
        hr_p = render_hr_diagram(star, sample_stars, ev_dir / "hr_diagram.png")
        generated_plots["sed"] = str(sed_p.relative_to(dossier_dir))
        generated_plots["cutouts"] = str(cut_p.relative_to(dossier_dir))
        generated_plots["lightcurves"] = str(lc_p.relative_to(dossier_dir))
        generated_plots["hr_diagram"] = str(hr_p.relative_to(dossier_dir))

    ch_entries: list[ChannelDossierEntry] = []
    for s in cal_scores:
        m = meas_by_ch.get(s.channel)
        ev_list: list[str] = []
        if s.channel == "ir_excess" and "sed" in generated_plots:
            ev_list.extend([generated_plots["sed"], generated_plots["cutouts"], generated_plots["hr_diagram"]])
        elif s.channel in {"ir_timedomain", "optical_dips"} and "lightcurves" in generated_plots:
            ev_list.append(generated_plots["lightcurves"])
        elif s.channel == "radio_narrowband" and "radio_narrowband" in generated_plots:
            ev_list.append(generated_plots["radio_narrowband"])

        # Channel statistic distribution plot
        stat_vals = sample_stats_by_channel.get(s.channel, [0.0, 0.5, -0.5])
        dist_p = render_statistic_distribution(
            s.channel,
            s.raw_statistic,
            stat_vals,
            5.0 if s.channel == "ir_excess" else 4.5,
            ev_dir / f"stat_{s.channel}.png",
        )
        ev_list.append(str(dist_p.relative_to(dossier_dir)))

        ch_entries.append(
            ChannelDossierEntry(
                name=s.channel,
                measurements=m.values if m else {},
                raw_statistic=s.raw_statistic,
                false_alarm_prob=s.false_alarm_prob,
                flagged=s.flagged,
                epochs=m.epochs if m else [],
                caveats=m.caveats if m else [],
                evidence=ev_list,
            )
        )

    short_hash = compute_config_hash(cfg.model_dump())[:12]
    cert_paths = [
        f"validation/certificates/{c.channel}-{c.gaia_release}-{c.config_hash}.json"
        for c in certs_by_channel.values()
    ]

    rel_prefix = "Gaia DR4" if cfg.gaia.release == "dr4" else "Gaia DR3"
    display_id = f"{rel_prefix} {clean_id}" if clean_id.isdigit() else clean_id

    dossier = Dossier(
        source_id=display_id,
        coords=Coords(
            ra_deg=float(star.get("ra_deg", 0.0)),
            dec_deg=float(star.get("dec_deg", 0.0)),
            epoch="J2016.0",
        ),
        tier=tier,
        provenance=Provenance(
            run_id=run_id,
            git_commit=current_git_commit(root),
            config_hash=short_hash,
            certificates=cert_paths,
            manifests=manifests,
        ),
        star=StarMetadata(
            distance_pc=star.get("distance_pc"),
            teff_k=star.get("teff_k"),
            luminosity_lsun=star.get("luminosity_lsun"),
            ruwe=star.get("ruwe"),
            simbad_type=star.get("simbad_type"),
            known_in_literature=bool(star.get("known_in_literature", False)),
            gaia_release=cfg.gaia.release,
        ),
        channels=ch_entries,
        vetting=vetting,
        hypotheses=hypotheses,
        cross_channel_patterns=patterns,
        reviews=reviews,
        review_status=review_status,
        what_would_settle_it=settle,
        lead_lines=leads,
        reproduce=f"oddstar inspect {clean_id} --config configs/default.yaml",
    )

    json_path = dossier_dir / "dossier.json"
    md_path = dossier_dir / "dossier.md"
    json_path.write_text(dossier.model_dump_json(indent=2) + "\n", encoding="utf-8")

    # Render Markdown dossier leading with the 4 required lines
    lines = [
        f"# Candidate Dossier: {display_id} (Tier {tier})",
        "",
    ]
    for lead in leads:
        lines.append(f"- **{lead.split(':', 1)[0]}**: {lead.split(':', 1)[1].strip()}")
    lines.extend(
        [
            "",
            f"- **Coordinates**: RA={dossier.coords.ra_deg:.5f} deg, Dec={dossier.coords.dec_deg:.5f} deg ({dossier.coords.epoch})",
            f"- **Provenance**: `run_id={run_id}`, `config_hash={short_hash}`, `gaia={cfg.gaia.release}`",
            f"- **Reproduce**: `{dossier.reproduce}`",
            "",
            "## Channel Measurements",
            "",
        ]
    )
    for ch in ch_entries:
        lines.append(f"### `{ch.name}` (flagged={ch.flagged}, FAP={ch.false_alarm_prob:.3e})")
        if ch.caveats:
            lines.append(f"> Caveat: {'; '.join(ch.caveats)}")
        for k, val in list(ch.measurements.items())[:12]:
            if val is not None:
                lines.append(f"- `{k}`: {val:.4g}")
        lines.append("")

    lines.extend(["## Vetting Checks", "", "| Check | Status | Note |", "| --- | --- | --- |"])
    for v in vetting:
        lines.append(f"| `{v.check}` | `{v.status}` | {v.note} |")

    lines.extend(["", "## Hypothesis Matrix", "", "| Hypothesis | Status | Note |", "| --- | --- | --- |"])
    for h in hypotheses:
        lines.append(f"| `{h.name}` | `{h.status}` | {h.note} |")

    if reviews:
        # Review status wording comes from oddstar.review.status (R10: independence is stated only when measured)
        lines.extend(["", "## Review Records", ""])
        if review_status is not None:
            lines.append(f"- **Review status**: {review_status.summary} ({review_status.independence_basis})")
        for r in reviews:
            fam = f", family `{r.model_family}`" if r.model_family else ""
            lines.append(f"- Round {r.round} ({r.reviewer_model}{fam}): verdict=`{r.verdict}` (`{r.record}`)")

    lines.extend(["", "## What Would Settle It", ""])
    for item in settle:
        lines.append(f"- {item}")

    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dossier, json_path, md_path
