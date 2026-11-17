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
from oddstar.io.physics import extinction_two_case_summary, star_applied_ebv
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


_LABEL_KEYS = (
    "photosphere_source",
    "photosphere_secondary",
    "locus_n_stars",
    "locus_id",
    "locus_gaia_release",
    "locus_sample_kind",
    "literature_reference",
    "literature_url",
    "literature_retrieved_utc_date",
    "ebv_star",
    "ebv_star_source",
)


def _channel_labels(m: Measurement | None) -> dict[str, Any]:
    """Provenance labels carried from Measurement.extra (photosphere source, locus size, applied reddening)."""
    if m is None:
        return {}
    return {k: m.extra[k] for k in _LABEL_KEYS if k in m.extra}


def _ir_excess_md(v: dict[str, float | None], unc: dict[str, float]) -> list[str]:
    def _f(key: str) -> float:
        val = v.get(key)
        return float(val) if val is not None else float("nan")

    out = [
        "",
        "| Band | Observed (mJy) | Photosphere (mJy) | Excess (mJy) | Excess / sigma | Obs / phot | Secondary model (mJy) | Bare Planck (mJy) |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for b in ("w1", "w2", "w3", "w4"):
        out.append(
            f"| {b.upper()} | {_f(f'{b}_obs_mjy'):.3f} +/- {unc.get(f'{b}_obs_mjy', float('nan')):.3f} "
            f"| {_f(f'{b}_phot_mjy'):.3f} +/- {unc.get(f'{b}_phot_mjy', float('nan')):.3f} "
            f"| {_f(f'{b}_excess_mjy'):.3f} +/- {unc.get(f'{b}_excess_mjy', float('nan')):.3f} "
            f"| {_f(f'{b}_excess_sigma'):.1f} | {_f(f'{b}_excess_ratio'):.2f} | {_f(f'{b}_pred_model_mjy'):.3f} | {_f(f'{b}_pred_planck_mjy'):.3f} |"
        )
    out.append("")
    if _f("t_dust_k") > 0.0:
        out.append(
            f"- Single-blackbody dust fit ({int(_f('dust_fit_n_bands'))} significant bands): T_dust = {_f('t_dust_k'):.0f} K "
            f"[{_f('t_dust_lo_k'):.0f}, {_f('t_dust_hi_k'):.0f}] (Delta chi2 interval), L_dust/L_star = {_f('frac_lum'):.4g} "
            f"[{_f('frac_lum_lo'):.4g}, {_f('frac_lum_hi'):.4g}], chi2 = {_f('dust_fit_chi2'):.2f} (reduced {_f('dust_fit_red_chi2'):.2f}; "
            f"interval uses Delta chi2 = {_f('dust_fit_interval_delta_chi2'):.2f}), r_bb = {_f('r_bb_au'):.2f} au"
        )
    else:
        out.append("- Single-blackbody dust fit: fewer than two bands with significant excess; temperature not constrained")
    for key, lab in (("t_color_w1w2", "W1-W2"), ("t_color_w3w4", "W3-W4")):
        if _f(f"{key}_valid") > 0.0:
            hi = f"{_f(f'{key}_hi_k'):.0f}" + ("+ (unbounded)" if _f(f"{key}_hi_unbounded") > 0.0 else "")
            out.append(f"- {lab} excess colour temperature: {_f(f'{key}_k'):.0f} K [{_f(f'{key}_lo_k'):.0f}, {hi}]")
        else:
            out.append(f"- {lab} excess colour temperature: not determined (excess not significant in both bands)")
    if _f("frac_lum_two_temp") > 0.0:
        out.append(
            f"- Two-temperature L_dust/L_star: hot {_f('frac_lum_hot'):.4g} + warm {_f('frac_lum_warm'):.4g} = {_f('frac_lum_two_temp'):.4g}"
        )
    out.append(
        f"- Detection statistics (two predictors): chi_emp_w3={_f('chi_emp_w3'):.2f}, chi_w3={_f('chi_w3'):.2f}, "
        f"chi_emp_w4={_f('chi_emp_w4'):.2f}, chi_w4={_f('chi_w4'):.2f}"
    )
    return out


def _extinction_md(ext: dict[str, Any]) -> list[str]:
    out = [
        "## Extinction: Two Cases",
        "",
        f"> {ext['triggered_by']}. {ext['summary']}",
        "",
        "| Case | E(B-V) | Source | (BP-RP)_0 | M_G | Dwarf class | Teff from colour (K) | M_G - dwarf sequence |",
        "| --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for c in ext["cases"]:
        out.append(
            f"| {c['case']} | {c['ebv']:.3f} | {c['ebv_source']} | {c['bp_rp_0']:.3f} | {c['m_g']:.2f} | {c['spectral_type']} "
            f"| {c['teff_k_from_colour']:.0f} | {c['delta_m_g_vs_sequence']:+.2f} ({c['luminosity_vs_sequence']}) |"
        )
    out.append("")
    return out


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
                labels=_channel_labels(m),
            )
        )

    short_hash = compute_config_hash(cfg.model_dump())[:12]
    cert_paths = [
        f"validation/certificates/{c.channel}-{c.gaia_release}-{c.config_hash}.json"
        for c in certs_by_channel.values()
    ]

    ebv_applied, ebv_applied_src = star_applied_ebv(star) if "mags" in star else (None, None)
    ext_cases = (
        extinction_two_case_summary(
            star,
            min_los_excess_ebv=float(cfg.vetting["extinction_los_excess_summary_ebv"]),
            sfd_to_sf11=float(cfg.vetting["sfd_to_sf11_factor"]),
        )
        if "mags" in star
        else None
    )

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
            ebv_star=ebv_applied,
            ebv_star_source=ebv_applied_src,
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
        extinction_cases=ext_cases,
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
        ]
    )
    if ebv_applied_src is not None:
        lines.append(
            f"- **Applied reddening**: E(B-V)={'n/a' if ebv_applied is None else f'{ebv_applied:.3f}'} (`ebv_star_source={ebv_applied_src}`)"
        )
    lines.append("")
    if ext_cases is not None:
        lines.extend(_extinction_md(ext_cases))
    lines.extend(["## Channel Measurements", ""])
    for ch in ch_entries:
        lines.append(f"### `{ch.name}` (flagged={ch.flagged}, FAP={ch.false_alarm_prob:.3e})")
        if ch.caveats:
            lines.append(f"> Caveat: {'; '.join(ch.caveats)}")
        if ch.labels:
            lines.append("- Labels: " + ", ".join(f"`{k}={v}`" for k, v in ch.labels.items() if v is not None))
        if ch.name == "ir_excess" and "w1_phot_mjy" in ch.measurements:
            lines.extend(_ir_excess_md(ch.measurements, meas_by_ch["ir_excess"].uncertainties if "ir_excess" in meas_by_ch else {}))
        else:
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
