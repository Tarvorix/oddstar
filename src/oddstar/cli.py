"""Single JSON-first CLI entry point for the oddstar pipeline (Section 10)."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from oddstar import __version__
from oddstar.channels import ALL_CHANNEL_NAMES, build_channel
from oddstar.config import PipelineConfig, compute_config_hash, current_git_commit, load_config, repo_root
from oddstar.dossier import write_dossier, write_null_report
from oddstar.errors import ArtifactNotFoundError, OddstarError, PolicyViolationError
from oddstar.fuse import fuse_results
from oddstar.io import (
    ArchiveFetcher,
    RADIO_CADENCES,
    RawCache,
    TableStore,
    enrich_star_live_for_vetting,
    fetch_live_halo_bulge_sample,
    fetch_single_star_live,
    generate_field_sample,
    get_named_object,
    list_archives_info,
    ra_dec_to_healpix,
)
from oddstar.review import add_lesson, build_review_packet, ingest_review_record
from oddstar.sample import build_parent_sample, resolve_star_record
from oddstar.schemas import (
    CLIResponse,
    ChannelScore,
    DossierReviewSummary,
    Measurement,
    RegistryEntry,
    ReviewRecord,
    RunState,
    SCHEMA_MODELS,
    SampleSpec,
    ValidationCertificate,
    VettingResult,
    export_all_schemas,
)
from oddstar.validate import load_certificate, require_valid_certificate, run_validation
from oddstar.vet import vet_target


def _emit(resp: CLIResponse, as_json: bool) -> None:
    if as_json:
        print(resp.model_dump_json(indent=2))
    else:
        print(resp.summary)
        for art in resp.artifacts:
            print(f"  artifact: {art}")
        for w in resp.warnings:
            print(f"  warning: {w}")
        if resp.next:
            print(f"  next: {resp.next}")


def _cmd_doctor(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    lock_file = root / "requirements.lock"
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    data = {
        "version": __version__,
        "python_version": sys.version.split()[0],
        "config_hash": short_hash,
        "git_commit": current_git_commit(root),
        "gaia_release": cfg.gaia.release,
        "lock_file_present": lock_file.exists(),
        "archives_registered": len(cfg.archives),
        "channels_registered": len(ALL_CHANNEL_NAMES),
        "cache": cache.stats(),
        "schemas_dir_exists": (root / "schemas").exists(),
    }
    warnings = []
    if not lock_file.exists():
        warnings.append("requirements.lock not found; generate with pip freeze > requirements.lock")
    return CLIResponse(
        ok=True,
        summary=(
            f"oddstar {__version__} healthy (Python {data['python_version']}, "
            f"cache_cap={cfg.storage.max_cache_gb} GB, gaia={cfg.gaia.release}, config={short_hash})"
        ),
        artifacts=[str(lock_file)] if lock_file.exists() else [],
        warnings=warnings,
        next="oddstar status --json",
        data=data,
    )


def _cmd_sources_list(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    rows = list_archives_info(cfg)
    return CLIResponse(
        ok=True,
        summary=f"{len(rows)} public archives registered (gaia release={cfg.gaia.release}).",
        data={"archives": rows},
        next="oddstar fetch gaia --ids 4582145718700356096 --json",
    )


def _cmd_fetch(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    ids = [x.strip() for x in args.ids.split(",") if x.strip()] if args.ids else None
    cone = (float(args.cone[0]), float(args.cone[1]), float(args.cone[2])) if args.cone else None
    if args.dry_run:
        return CLIResponse(
            ok=True,
            summary=f"[dry-run] Would fetch from '{args.archive}' (ids={ids}, cone={cone}, adql={bool(args.adql)})",
            data={"dry_run": True, "archive": args.archive},
        )
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    fetcher = ArchiveFetcher(args.archive, cfg, cache)
    rows, manifest, cache_hit = fetcher.fetch(ids=ids, cone=cone, adql=args.adql, live=bool(args.live))
    if args.limit is not None:
        rows = rows[: args.limit]
    return CLIResponse(
        ok=True,
        summary=(
            f"Fetched {len(rows)} rows from '{args.archive}' ({fetcher.release()}) "
            f"[cache_hit={cache_hit}, backend_calls={fetcher.backend_calls}, sha256={manifest.sha256[:12]}]"
        ),
        artifacts=[manifest.blob_path, str(cache.manifest_path)],
        data={
            "row_count": len(rows),
            "cache_hit": cache_hit,
            "backend_calls": fetcher.backend_calls,
            "manifest": manifest.model_dump(),
        },
        next=f"oddstar sample build ir_excess --json",
    )


def _cmd_sample_build(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(
            ok=True,
            summary=f"[dry-run] Would build parent sample for '{args.channel}' (limit={args.limit}, live={bool(args.live)})",
        )
    if args.live:
        cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
        n_lim = int(args.limit if args.limit is not None else 300)
        stars, _ = fetch_live_halo_bulge_sample(
            cfg,
            cache,
            limit=n_lim,
            region=args.region,
            after_random_index=args.after_random_index,
        )
        nside = cfg.storage.healpix_nside
        store = TableStore(root / cfg.storage.tables_dir, nside=nside)
        store.write_partitioned(f"sample_live_{args.channel}_{cfg.gaia.release}", stars)
        cells = sorted({ra_dec_to_healpix(float(s["ra_deg"]), float(s["dec_deg"]), nside=nside) for s in stars})
        spec = SampleSpec(
            channel=args.channel,
            gaia_release=cfg.gaia.release,
            config_hash=short_hash,
            seed=cfg.sample.seed,
            cuts=[],
            total_count=len(stars),
            healpix_cells=cells,
            source_ids=[str(s["source_id"]) for s in stars],
        )
    else:
        ch = build_channel(args.channel, cfg, root)
        spec = ch.sample(cfg, limit=args.limit)
    out_dir = root / cfg.storage.runs_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"sample_{args.channel}_{cfg.gaia.release}.json"
    out_path.write_text(spec.model_dump_json(indent=2) + "\n", encoding="utf-8")
    return CLIResponse(
        ok=True,
        summary=f"Built parent sample for '{args.channel}' ({spec.total_count} stars across {len(spec.healpix_cells)} HEALPix cells, live={bool(args.live)})",
        artifacts=[str(out_path)],
        data={"sample": spec.model_dump(exclude={"source_ids"})},
        next=f"oddstar validate {args.channel} --json",
    )


def _cmd_validate(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(
            ok=True,
            summary=f"[dry-run] Would run 4 validation tests for '{args.channel}' (split={args.split})",
        )
    field_limit = args.limit if args.limit is not None else min(1200, cfg.validation.field_sample_size)
    cert, cert_path = run_validation(
        args.channel,
        cfg,
        root,
        split=args.split,
        field_limit=field_limit,
    )
    return CLIResponse(
        ok=True,
        summary=(
            f"Certified '{args.channel}' ({cfg.gaia.release}, split={args.split}): "
            f"blind={cert.blind_recovered}/{cert.blind_total} ({cert.blind_recovery_rate:.0%}), "
            f"impostor_max_tier={cert.impostor_max_tier}, "
            f"inj_dp={cert.injection_completeness_design_point:.0%}, "
            f"FAR={cert.false_alarm_rate:.4f}"
        ),
        artifacts=[str(cert_path)],
        data={"certificate": cert.model_dump(exclude={"details"})},
        next=f"oddstar search {args.channel} --json",
    )


def _cmd_search(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    cert = require_valid_certificate(args.channel, cfg, root)
    run_id = args.run_id or f"run-{args.channel}-{cfg.gaia.release}-{short_hash}"
    run_dir = root / cfg.storage.runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    search_out = run_dir / f"search_{args.channel}.json"

    if args.dry_run:
        return CLIResponse(
            ok=True,
            summary=f"[dry-run] Certificate valid; would run '{args.channel}' search into {run_id}",
        )

    seed = int(args.seed if args.seed is not None else cfg.sample.seed)
    n_stars = int(args.limit if args.limit is not None else min(1000, cfg.sample.field_sample_size))
    is_live = bool(args.live)
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)

    extra_stars: list[dict[str, Any]] = []
    if args.extra_ids:
        for raw_id in args.extra_ids.split(","):
            cid = raw_id.strip()
            if cid:
                if is_live:
                    st, _ = fetch_single_star_live(cid, cfg, cache)
                    extra_stars.append(st)
                else:
                    extra_stars.append(resolve_star_record(cid, cfg))

    ch = build_channel(args.channel, cfg, root)
    if args.channel == "radio_narrowband":
        stars = [dict(v) for v in RADIO_CADENCES.values()]
        cuts = []
        healpix_cells = []
    elif is_live:
        stars, _ = fetch_live_halo_bulge_sample(
            cfg,
            cache,
            limit=n_stars,
            region=args.region,
            after_random_index=args.after_random_index,
        )
        if extra_stars:
            by_id = {s["source_id"]: s for s in stars}
            for es in extra_stars:
                by_id[es["source_id"]] = es
            stars = list(by_id.values())
        if args.channel == "ir_timedomain":
            # Enrich up to 25 targets with live NEOWISE-R multi-epoch data per Section 5.2 cost control
            stars = [enrich_star_live_for_vetting(s, cfg, cache)[0] for s in stars[: min(25, len(stars))]]
        nside = cfg.storage.healpix_nside
        cuts = [
            {
                "cut_name": f"gaia_{args.region}_live",
                "description": f"Live Gaia DR3 x AllWISE x 2MASS region='{args.region}' (n={len(stars)}, after_random_index={args.after_random_index})",
                "count_before": len(stars),
                "count_after": len(stars),
            }
        ]
        healpix_cells = sorted({ra_dec_to_healpix(float(s["ra_deg"]), float(s["dec_deg"]), nside=nside) for s in stars})
    else:
        spec, stars = build_parent_sample(
            args.channel,
            cfg,
            short_hash,
            root,
            extra_stars=extra_stars,
            limit=n_stars,
        )
        cuts = [c.model_dump() for c in spec.cuts]
        healpix_cells = spec.healpix_cells

    if args.channel == "ir_excess":
        ch.set_empirical_locus(stars)
    elif args.channel == "underluminous":
        ch.fit_locus(stars)
    elif args.channel == "generic_outlier":
        ch.fit_background(stars)

    measurements: list[dict[str, Any]] = []
    scores: list[dict[str, Any]] = []
    flagged_ids: list[str] = []
    stars_by_id = {str(s["source_id"]): s for s in stars}
    for star in stars:
        m = ch.measure(star, {})
        sc = ch.score(m)
        # Calibrate FAP immediately against the active certificate
        from oddstar.fuse.engine import calibrate_fap

        cal_fap = calibrate_fap(sc.raw_statistic, cert, sc.false_alarm_prob)
        sc = sc.model_copy(update={"false_alarm_prob": cal_fap})
        if is_live and args.channel == "ir_excess" and sc.flagged:
            aw = star.get("allwise", {})
            phq = str(aw.get("ph_qual", "AAAA"))
            ccf = str(aw.get("cc_flags", "0000"))
            w4_upper = len(phq) >= 4 and phq[3] in {"U", "X"}
            ir_artifact = len(ccf) >= 4 and (ccf[2] != "0" or ccf[3] != "0")
            if w4_upper or ir_artifact:
                reason = "W4_upper_limit" if w4_upper else f"cc_flags={ccf}"
                sc = sc.model_copy(
                    update={
                        "flagged": False,
                        "summary": f"{sc.summary} [filtered_live_ingestion: {reason}]",
                    }
                )
        measurements.append(m.model_dump())
        scores.append(sc.model_dump())
        if sc.flagged:
            flagged_ids.append(sc.source_id)

    if is_live and flagged_ids:
        live_flagged: dict[str, Any] = {}
        for fsid in flagged_ids:
            enriched, _ = enrich_star_live_for_vetting(stars_by_id[fsid], cfg, cache)
            live_flagged[fsid] = enriched
        (run_dir / "live_stars.json").write_text(json.dumps(live_flagged, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    payload = {
        "run_id": run_id,
        "channel": args.channel,
        "gaia_release": cfg.gaia.release,
        "config_hash": short_hash,
        "git_commit": current_git_commit(root),
        "seed": seed,
        "certificate": f"validation/certificates/{args.channel}-{cfg.gaia.release}-{short_hash}.json",
        "n_searched": len(stars),
        "n_flagged": len(flagged_ids),
        "flagged_ids": flagged_ids,
        "cuts": cuts,
        "healpix_cells": healpix_cells,
        "measurements": measurements,
        "scores": scores,
    }
    search_out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    state = RunState(
        run_id=run_id,
        stage="search",
        channel=args.channel,
        gaia_release=cfg.gaia.release,
        config_hash=short_hash,
        git_commit=current_git_commit(root),
        seed=seed,
        status="completed",
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        artifacts=[str(search_out)],
        summary={"n_searched": len(stars), "n_flagged": len(flagged_ids)},
    )
    (run_dir / "run_state.json").write_text(state.model_dump_json(indent=2) + "\n", encoding="utf-8")

    return CLIResponse(
        ok=True,
        summary=f"Searched {len(stars)} targets on '{args.channel}' -> {len(flagged_ids)} flagged (run_id={run_id})",
        artifacts=[str(search_out)],
        data={"run_id": run_id, "n_searched": len(stars), "n_flagged": len(flagged_ids), "flagged_ids": flagged_ids},
        next=f"oddstar vet {run_id} --json",
    )


def _cmd_vet(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    run_dir = root / cfg.storage.runs_dir / args.run_id
    search_files = sorted(run_dir.glob("search_*.json"))
    if not search_files:
        raise ArtifactNotFoundError(f"No search outputs found in {run_dir}")

    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would vet flagged targets in {args.run_id}")

    vet_by_star: dict[str, list[dict[str, Any]]] = {}
    ev_base = run_dir / "vet_evidence"
    live_stars_file = run_dir / "live_stars.json"
    live_stars_map = json.loads(live_stars_file.read_text(encoding="utf-8")) if live_stars_file.exists() else {}
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb) if args.live else None
    for sf in search_files:
        sdata = json.loads(sf.read_text(encoding="utf-8"))
        scores_by_id = {row["source_id"]: ChannelScore.model_validate(row) for row in sdata["scores"]}
        meas_by_id = {m["source_id"]: m.get("values", {}) for m in sdata.get("measurements", [])}
        for sid in sdata.get("flagged_ids", []):
            if sid in live_stars_map:
                star = dict(live_stars_map[sid])
            elif args.live and cache is not None and sid not in RADIO_CADENCES:
                star, _ = fetch_single_star_live(sid, cfg, cache)
            elif sid in RADIO_CADENCES:
                star = dict(RADIO_CADENCES[sid])
            else:
                star = resolve_star_record(sid, cfg)
            mvals = meas_by_id.get(sid, {})
            if "drift_hz_s" in mvals:
                star["drift_hz_s"] = mvals["drift_hz_s"]
            if "off_scan_hits" in mvals:
                star["off_scan_hits"] = int(mvals["off_scan_hits"])
            sc = scores_by_id[sid]
            vres = vet_target(star, [sc], cfg, ev_base / sid)
            vet_by_star.setdefault(sid, []).extend([v.model_dump() for v in vres])

    vet_out = run_dir / "vetting.json"
    vet_out.write_text(json.dumps({"run_id": args.run_id, "vetting": vet_by_star}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return CLIResponse(
        ok=True,
        summary=f"Vetted {len(vet_by_star)} flagged targets in {args.run_id}",
        artifacts=[str(vet_out)],
        data={"run_id": args.run_id, "vetted_count": len(vet_by_star)},
        next=f"oddstar fuse {args.run_id} --json",
    )


def _cmd_fuse(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    scores_by_star: dict[str, list[ChannelScore]] = {}
    vet_by_star: dict[str, list[VettingResult]] = {}
    certs_by_ch: dict[str, ValidationCertificate] = {}
    n_stars_max = 0

    for rid in args.run_ids:
        rdir = root / cfg.storage.runs_dir / rid
        if not rdir.exists():
            raise ArtifactNotFoundError(f"Run directory not found: {rdir}")
        for sf in sorted(rdir.glob("search_*.json")):
            sdata = json.loads(sf.read_text(encoding="utf-8"))
            ch = sdata["channel"]
            n_stars_max = max(n_stars_max, int(sdata.get("n_searched", 0)))
            cert = load_certificate(ch, cfg, root)
            if cert is not None:
                certs_by_ch[ch] = cert
            for row in sdata["scores"]:
                sc = ChannelScore.model_validate(row)
                if sc.flagged:
                    scores_by_star.setdefault(sc.source_id, []).append(sc)
        vf = rdir / "vetting.json"
        if vf.exists():
            vdata = json.loads(vf.read_text(encoding="utf-8"))
            for sid, vlist in vdata.get("vetting", {}).items():
                for vrow in vlist:
                    vet_by_star.setdefault(sid, []).append(VettingResult.model_validate(vrow))

    out_run_id = args.out_run_id or f"fuse-{cfg.gaia.release}-{short_hash}"
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would fuse {len(scores_by_star)} flagged targets across {len(args.run_ids)} runs")

    fusion = fuse_results(
        run_id=out_run_id,
        scores_by_star=scores_by_star,
        vetting_by_star=vet_by_star,
        certs_by_channel=certs_by_ch,
        n_stars_total=n_stars_max,
        cfg=cfg,
        root=root,
    )
    out_dir = root / cfg.storage.runs_dir / out_run_id
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "fusion.json"
    out_path.write_text(fusion.model_dump_json(indent=2) + "\n", encoding="utf-8")

    # Append fusion outcomes to registry ledger
    reg_dir = root / cfg.storage.registry_dir
    reg_dir.mkdir(parents=True, exist_ok=True)
    ledger = reg_dir / "ledger.jsonl"
    now = datetime.now(timezone.utc).isoformat()
    with ledger.open("a", encoding="utf-8") as fh:
        for cand in fusion.candidates:
            ev_type = "explained_anomaly" if cand.tier == "Explained" else "candidate_registered"
            entry = RegistryEntry(
                entry_id=f"REG-FUSE-{out_run_id}-{cand.source_id}",
                timestamp_utc=now,
                source_id=cand.source_id,
                event_type=ev_type,
                tier=cand.tier,
                status=cand.explained_reason or f"Tier {cand.tier}",
                provenance={"run_id": out_run_id, "config_hash": short_hash},
                details=cand.model_dump(),
            )
            fh.write(entry.model_dump_json() + "\n")

    tier_counts: dict[str, int] = {}
    for c in fusion.candidates:
        tier_counts[c.tier] = tier_counts.get(c.tier, 0) + 1

    return CLIResponse(
        ok=True,
        summary=f"Fused {len(fusion.candidates)} flagged targets across {len(certs_by_ch)} channels (tiers={tier_counts})",
        artifacts=[str(out_path), str(ledger)],
        data={"fusion": fusion.model_dump(), "tier_counts": tier_counts},
        next=f"oddstar null-report {args.run_ids[0]} --json",
    )


def _inspect_single_target(
    source_id: str,
    cfg: PipelineConfig,
    short_hash: str,
    root: Path,
    *,
    run_id: str | None = None,
    live: bool = False,
) -> tuple[Any, Path, Path]:
    clean_id = source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
    rid = run_id or f"inspect-{clean_id}-{cfg.gaia.release}-{short_hash}"
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)

    if live:
        star, live_mans = fetch_single_star_live(clean_id, cfg, cache)
        sample_stars, bg_mans = fetch_live_halo_bulge_sample(cfg, cache, limit=1000)
        manifests = [m.manifest_id for m in (*live_mans, *bg_mans)]
    else:
        # Record stage-1 raw fetch manifests across primary archives
        manifests = []
        for arc in ("gaia", "2mass", "allwise", "neowise", "tess", "simbad", "radio_continuum", "deep_optical", "mid_far_ir"):
            fetcher = ArchiveFetcher(arc, cfg, cache)
            _, man, _ = fetcher.fetch(ids=[clean_id])
            manifests.append(man.manifest_id)

        star = RADIO_CADENCES.get(clean_id) or resolve_star_record(clean_id, cfg)
        sample_stars = generate_field_sample(n_stars=300, seed=cfg.sample.seed, gaia_release=cfg.gaia.release)

    measurements: list[Measurement] = []
    scores: list[ChannelScore] = []
    stats_by_ch: dict[str, list[float]] = {}
    certs_by_ch: dict[str, ValidationCertificate] = {}

    eval_channels = ("radio_narrowband",) if "scans" in star else (
        "ir_excess",
        "ir_timedomain",
        "optical_dips",
        "underluminous",
        "spectral_lines",
        "vanishing",
        "generic_outlier",
    )

    for cname in eval_channels:
        ch = build_channel(cname, cfg, root)
        if cname == "ir_excess":
            ch.set_empirical_locus(sample_stars)
        elif cname == "underluminous":
            ch.fit_locus(sample_stars)
        elif cname == "generic_outlier":
            ch.fit_background(sample_stars)
        m = ch.measure(star, {})
        sc = ch.score(m)
        measurements.append(m)
        scores.append(sc)
        cert = load_certificate(cname, cfg, root)
        if cert is not None:
            certs_by_ch[cname] = cert
        # Background histogram sample (fast 50 stars)
        if cname != "radio_narrowband":
            stats_by_ch[cname] = [ch.score(ch.measure(s, {})).raw_statistic for s in sample_stars[:40]]
        else:
            stats_by_ch[cname] = [0.0, 1.0, 2.0]

    ev_dir = root / cfg.storage.dossiers_dir / clean_id / "evidence"
    vetting = vet_target(star, scores, cfg, ev_dir)

    # Load any existing ingested reviews for this target
    reviews: list[DossierReviewSummary] = []
    rec_dir = root / "review" / "records"
    if rec_dir.exists():
        for rf in sorted(rec_dir.glob(f"*{clean_id}*.json")):
            rr = ReviewRecord.model_validate_json(rf.read_text(encoding="utf-8"))
            reviews.append(
                DossierReviewSummary(
                    round=rr.round,
                    record=str(rf.relative_to(root)),
                    verdict=rr.verdict,
                    reviewer_model=rr.reviewer_model,
                )
            )

    dossier, json_p, md_p = write_dossier(
        source_id=clean_id,
        star=star,
        measurements=measurements,
        scores=scores,
        vetting=vetting,
        sample_stars=sample_stars,
        sample_stats_by_channel=stats_by_ch,
        certs_by_channel=certs_by_ch,
        reviews=reviews,
        run_id=rid,
        manifests=manifests,
        cfg=cfg,
        root=root,
    )
    return dossier, json_p, md_p


def _cmd_inspect(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would inspect {args.source_id} across all channels and write dossier")
    dossier, json_p, md_p = _inspect_single_target(args.source_id, cfg, short_hash, root, run_id=args.run_id, live=bool(args.live))
    flagged = [c.name for c in dossier.channels if c.flagged]
    holes = [v.check for v in dossier.vetting if v.status == "not_evaluated"]
    return CLIResponse(
        ok=True,
        summary=f"Inspected {dossier.source_id} -> Tier {dossier.tier} (flagged={flagged}, open_holes={holes})",
        artifacts=[str(json_p), str(md_p)],
        data={"tier": dossier.tier, "flagged_channels": flagged, "open_holes": holes, "lead_lines": dossier.lead_lines},
        next=f"oddstar review packet {args.source_id} --round 1 --json",
    )


def _cmd_dossier(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would generate dossier for {args.source_id}")
    dossier, json_p, md_p = _inspect_single_target(args.source_id, cfg, short_hash, root, run_id=args.run_id, live=bool(args.live))
    return CLIResponse(
        ok=True,
        summary=f"Dossier generated for {dossier.source_id} (Tier {dossier.tier})",
        artifacts=[str(json_p), str(md_p)],
        data={"tier": dossier.tier, "dossier_path": str(json_p)},
    )


def _cmd_review_packet(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    ids = [x.strip() for x in args.ids.split(",") if x.strip()]
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would build Round {args.round} review packet for {ids}")
    pkt, pkt_p, key_p = build_review_packet(ids, cfg, root, round_num=args.round, claim=args.claim)
    return CLIResponse(
        ok=True,
        summary=f"Built Round {args.round} review packet {pkt.packet_id} ({len(pkt.targets)} targets; key isolated)",
        artifacts=[str(pkt_p)],
        data={"packet_id": pkt.packet_id, "round": pkt.round, "target_count": len(pkt.targets)},
        next=f"oddstar review ingest <review_record.json> --json",
    )


def _cmd_review_ingest(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would ingest review record {args.file}")
    res = ingest_review_record(Path(args.file), cfg, root)
    return CLIResponse(
        ok=True,
        summary=f"Ingested review {res['review_id']} (discrepancy={res['discrepancy']}, open_findings={len(res['open_findings'])})",
        artifacts=[res["stored_path"]],
        data=res,
        next="oddstar lessons add --help" if res["open_findings"] else "oddstar status --json",
    )


def _cmd_lessons_add(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would create lesson {args.lesson_id} for finding {args.finding_id}")
    rec, md_p, json_p = add_lesson(
        lesson_id=args.lesson_id,
        finding_id=args.finding_id,
        source_id=args.source_id,
        title=args.title,
        what_happened=args.what_happened,
        generic_lesson=args.generic_lesson,
        automated_check=args.automated_check,
        linked_test=args.linked_test,
        added_to_impostors=bool(args.added_to_impostors),
        skill_updated=args.skill,
        cfg=cfg,
        root=root,
    )
    return CLIResponse(
        ok=True,
        summary=f"Recorded lesson {rec.lesson_id} closing finding {rec.finding_id} -> {rec.linked_test}",
        artifacts=[str(md_p), str(json_p)],
        data={"lesson": rec.model_dump()},
        next="oddstar status --json",
    )


def _cmd_null_report(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    run_dir = root / cfg.storage.runs_dir / args.run_id
    search_files = sorted(run_dir.glob("search_*.json"))
    if not search_files:
        raise ArtifactNotFoundError(f"No search run found at {run_dir}")
    sdata = json.loads(search_files[0].read_text(encoding="utf-8"))
    ch_name = args.channel or sdata["channel"]
    cert = load_certificate(ch_name, cfg, root)
    if cert is None:
        raise ArtifactNotFoundError(f"Certificate missing for {ch_name}")

    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would generate null-report for {args.run_id} ({ch_name})")

    stars = generate_field_sample(n_stars=int(sdata["n_searched"]), seed=int(sdata.get("seed", cfg.sample.seed)), gaia_release=cfg.gaia.release)
    tier_breakdown = {"A": 0, "B": 0, "C": 0, "Explained": 0, "Rejected": 0}
    vf = run_dir / "vetting.json"
    if vf.exists():
        vdata = json.loads(vf.read_text(encoding="utf-8")).get("vetting", {})
        from oddstar.fuse.engine import assign_tier

        scores_by_id = {r["source_id"]: ChannelScore.model_validate(r) for r in sdata["scores"]}
        for sid in sdata.get("flagged_ids", []):
            vres = [VettingResult.model_validate(v) for v in vdata.get(sid, [])]
            t, _ = assign_tier([scores_by_id[sid]], vres, cfg)
            tier_breakdown[t] = tier_breakdown.get(t, 0) + 1
    else:
        tier_breakdown["C"] = int(sdata["n_flagged"])

    from oddstar.schemas import SampleCutStep

    cuts = [SampleCutStep.model_validate(c) for c in sdata.get("cuts", [])]
    rep, jp, mp = write_null_report(
        channel=ch_name,
        run_id=args.run_id,
        sample_stars=stars,
        sample_cuts=cuts,
        healpix_cells=sdata.get("healpix_cells", []),
        n_flagged=int(sdata["n_flagged"]),
        tier_breakdown=tier_breakdown,
        cert=cert,
        manifests=[],
        cfg=cfg,
        root=root,
    )
    return CLIResponse(
        ok=True,
        summary=f"Null report emitted for '{ch_name}' ({rep.n_searched} searched, {rep.n_flagged} flagged, {rep.expected_false_alarms:.2f} expected FA)",
        artifacts=[str(jp), str(mp)],
        data={"null_report": rep.model_dump(exclude={"completeness_grid"})},
    )


def _cmd_status(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    certs_dir = root / "validation" / "certificates"
    certified: dict[str, list[str]] = {}
    if certs_dir.exists():
        for cf in sorted(certs_dir.glob("*.json")):
            if cf.name == "dr3_vs_dr4_comparison.json":
                continue
            c = ValidationCertificate.model_validate_json(cf.read_text(encoding="utf-8"))
            if c.passed:
                certified.setdefault(c.gaia_release, []).append(c.channel)

    open_findings: list[str] = []
    closed_findings: set[str] = set()
    ledger = root / cfg.storage.registry_dir / "ledger.jsonl"
    if ledger.exists():
        for line in ledger.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("event_type") == "finding_opened":
                open_findings.append(row["details"]["finding_id"])
            elif row.get("event_type") == "finding_closed":
                closed_findings.add(row["details"]["finding_id"])
    unresolved = [fid for fid in open_findings if fid not in closed_findings]

    runs_dir = root / cfg.storage.runs_dir
    run_ids = sorted(p.name for p in runs_dir.iterdir() if p.is_dir() and not p.name.startswith("_")) if runs_dir.exists() else []
    uncertified_dr3 = [c for c in ALL_CHANNEL_NAMES if c not in certified.get("dr3", [])]
    next_cmd = f"oddstar validate {uncertified_dr3[0]} --json" if uncertified_dr3 else "oddstar inspect 4582145718700356096 --json"
    return CLIResponse(
        ok=True,
        summary=(
            f"Certified channels DR3={len(certified.get('dr3', []))}/8, DR4={len(certified.get('dr4', []))}/8; "
            f"active runs={len(run_ids)}; open findings={len(unresolved)}"
        ),
        data={
            "certified": certified,
            "uncertified_dr3": uncertified_dr3,
            "active_runs": run_ids,
            "open_findings": unresolved,
        },
        next=next_cmd,
    )


def _cmd_show(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    path = Path(args.path)
    if not path.is_absolute():
        path = root / path
    if not path.exists():
        raise ArtifactNotFoundError(f"Artifact not found: {path}")

    limit = max(1, int(args.rows))
    offset = max(0, int(args.offset))
    text = path.read_text(encoding="utf-8")
    if path.suffix == ".jsonl":
        lines = [json.loads(l) for l in text.splitlines() if l.strip()]
        page = lines[offset : offset + limit]
        return CLIResponse(ok=True, summary=f"Showing rows {offset}..{offset+len(page)} of {len(lines)} from {path.name}", data={"total": len(lines), "rows": page})
    if path.suffix == ".json":
        obj = json.loads(text)
        if isinstance(obj, list):
            page = obj[offset : offset + limit]
            return CLIResponse(ok=True, summary=f"Showing items {offset}..{offset+len(page)} of {len(obj)} from {path.name}", data={"total": len(obj), "rows": page})
        return CLIResponse(ok=True, summary=f"Showing JSON object from {path.name}", data={"content": obj})

    lines_txt = text.splitlines()
    page_txt = lines_txt[offset : offset + limit]
    return CLIResponse(ok=True, summary=f"Showing lines {offset}..{offset+len(page_txt)} of {len(lines_txt)} from {path.name}", data={"total": len(lines_txt), "lines": page_txt})


def _cmd_schema(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.export:
        written = export_all_schemas(root)
        return CLIResponse(
            ok=True,
            summary=f"Exported {len(written)} JSON Schemas to schemas/ and review/review_record.schema.json",
            artifacts=[str(p) for p in written],
            data={"schemas": [p.name for p in written]},
        )
    if args.type:
        if args.type not in SCHEMA_MODELS:
            raise ArtifactNotFoundError(f"Unknown schema type: {args.type}", hint=f"Available: {', '.join(sorted(SCHEMA_MODELS))}")
        return CLIResponse(
            ok=True,
            summary=f"JSON Schema for '{args.type}'",
            data={"schema": SCHEMA_MODELS[args.type].model_json_schema()},
        )
    return CLIResponse(
        ok=True,
        summary=f"{len(SCHEMA_MODELS)} schemas registered",
        data={"types": sorted(SCHEMA_MODELS)},
    )


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--json", action="store_true", help="Emit structured JSON output")
    common.add_argument("--dry-run", action="store_true", help="Show what would run without modifying state")
    common.add_argument("--limit", type=int, default=None, help="Bound row/star count")
    common.add_argument("--config", type=str, default=None, help="Path to pipeline YAML config")
    common.add_argument("--release", type=str, default=None, choices=["dr3", "dr4"], help="Override Gaia data release")
    common.add_argument("--live", action="store_true", help="Query live public sky archives instead of offline fixtures")
    common.add_argument("--region", type=str, default="halo", choices=["halo", "bulge", "local"], help="Live sky region preset")
    common.add_argument("--after-random-index", type=int, default=None, help="Keyset pagination cursor for live Gaia random_index")

    parser = argparse.ArgumentParser(prog="oddstar", description="Reproducible multi-channel anomalous star search pipeline", parents=[common])
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("doctor", parents=[common], help="Check environment, lock file, cache cap, and archives")

    p_src = sub.add_parser("sources", parents=[common], help="Public archive registry")
    src_sub = p_src.add_subparsers(dest="sources_cmd", required=True)
    src_sub.add_parser("list", parents=[common], help="List registered archives")

    p_fetch = sub.add_parser("fetch", parents=[common], help="Stage 1 raw archive fetch into cache")
    p_fetch.add_argument("archive")
    p_fetch.add_argument("--ids", type=str, default=None)
    p_fetch.add_argument("--cone", nargs=3, type=float, default=None, metavar=("RA", "DEC", "ARCSEC"))
    p_fetch.add_argument("--adql", type=str, default=None)

    p_sample = sub.add_parser("sample", parents=[common], help="Stage 2 parent sample construction")
    sample_sub = p_sample.add_subparsers(dest="sample_cmd", required=True)
    p_sb = sample_sub.add_parser("build", parents=[common], help="Build parent sample for channel")
    p_sb.add_argument("channel")

    p_val = sub.add_parser("validate", parents=[common], help="Run 4-part validation harness and issue certificate")
    p_val.add_argument("channel")
    p_val.add_argument("--split", default="all", choices=["dev", "holdout", "all"])

    p_search = sub.add_parser("search", parents=[common], help="Stage 3 search (locked unless certified)")
    p_search.add_argument("channel")
    p_search.add_argument("--run-id", default=None)
    p_search.add_argument("--seed", type=int, default=None)
    p_search.add_argument("--extra-ids", default=None)

    p_vet = sub.add_parser("vet", parents=[common], help="Stage 4 vetting checks on flagged targets")
    p_vet.add_argument("run_id")

    p_fuse = sub.add_parser("fuse", parents=[common], help="Stage 5 cross-channel fusion and ranking")
    p_fuse.add_argument("run_ids", nargs="+")
    p_fuse.add_argument("--out-run-id", default=None)

    p_dos = sub.add_parser("dossier", parents=[common], help="Stage 6 dossier writer")
    p_dos.add_argument("source_id")
    p_dos.add_argument("--run-id", default=None)

    p_ins = sub.add_parser("inspect", parents=[common], help="One-star full pipeline mode")
    p_ins.add_argument("source_id")
    p_ins.add_argument("--run-id", default=None)

    p_rev = sub.add_parser("review", parents=[common], help="Stage 7 isolated adversarial review")
    rev_sub = p_rev.add_subparsers(dest="review_cmd", required=True)
    p_rp = rev_sub.add_parser("packet", parents=[common], help="Build blinded reviewer packet")
    p_rp.add_argument("ids")
    p_rp.add_argument("--round", type=int, default=1)
    p_rp.add_argument("--claim", default=None)
    p_ri = rev_sub.add_parser("ingest", parents=[common], help="Validate and record reviewer JSON")
    p_ri.add_argument("file")

    p_les = sub.add_parser("lessons", parents=[common], help="Lessons loop (R7)")
    les_sub = p_les.add_subparsers(dest="lessons_cmd", required=True)
    p_la = les_sub.add_parser("add", parents=[common], help="Add lesson and link test")
    p_la.add_argument("--lesson-id", required=True)
    p_la.add_argument("--finding-id", required=True)
    p_la.add_argument("--source-id", required=True)
    p_la.add_argument("--title", required=True)
    p_la.add_argument("--what-happened", required=True)
    p_la.add_argument("--generic-lesson", required=True)
    p_la.add_argument("--automated-check", required=True)
    p_la.add_argument("--linked-test", required=True)
    p_la.add_argument("--added-to-impostors", action="store_true")
    p_la.add_argument("--skill", default=None)

    p_null = sub.add_parser("null-report", parents=[common], help="Emit sensitivity and null-result report")
    p_null.add_argument("run_id")
    p_null.add_argument("--channel", default=None)

    sub.add_parser("status", parents=[common], help="Report certificates, runs, and open findings")

    p_show = sub.add_parser("show", parents=[common], help="Safely inspect artifacts in bounded pages")
    p_show.add_argument("path")
    p_show.add_argument("--rows", type=int, default=20)
    p_show.add_argument("--offset", type=int, default=0)

    p_sch = sub.add_parser("schema", parents=[common], help="Inspect or export JSON Schemas")
    p_sch.add_argument("type", nargs="?", default=None)
    p_sch.add_argument("--export", action="store_true")

    return parser


def dispatch(argv: list[str] | None = None, *, root: Path | None = None) -> tuple[int, dict[str, Any]]:
    """Execute CLI command programmatically and return (exit_code, response_dict)."""
    parser = build_parser()
    args = parser.parse_args(argv)
    base = root or repo_root()
    as_json = bool(args.json)

    try:
        cfg_path = Path(args.config) if args.config else None
        cfg, short_hash, _ = load_config(cfg_path, root=base, gaia_release_override=args.release)

        cmd = args.command
        if cmd == "doctor":
            resp = _cmd_doctor(args, cfg, short_hash, base)
        elif cmd == "sources" and args.sources_cmd == "list":
            resp = _cmd_sources_list(args, cfg, short_hash, base)
        elif cmd == "fetch":
            resp = _cmd_fetch(args, cfg, short_hash, base)
        elif cmd == "sample" and args.sample_cmd == "build":
            resp = _cmd_sample_build(args, cfg, short_hash, base)
        elif cmd == "validate":
            resp = _cmd_validate(args, cfg, short_hash, base)
        elif cmd == "search":
            resp = _cmd_search(args, cfg, short_hash, base)
        elif cmd == "vet":
            resp = _cmd_vet(args, cfg, short_hash, base)
        elif cmd == "fuse":
            resp = _cmd_fuse(args, cfg, short_hash, base)
        elif cmd == "dossier":
            resp = _cmd_dossier(args, cfg, short_hash, base)
        elif cmd == "inspect":
            resp = _cmd_inspect(args, cfg, short_hash, base)
        elif cmd == "review" and args.review_cmd == "packet":
            resp = _cmd_review_packet(args, cfg, short_hash, base)
        elif cmd == "review" and args.review_cmd == "ingest":
            resp = _cmd_review_ingest(args, cfg, short_hash, base)
        elif cmd == "lessons" and args.lessons_cmd == "add":
            resp = _cmd_lessons_add(args, cfg, short_hash, base)
        elif cmd == "null-report":
            resp = _cmd_null_report(args, cfg, short_hash, base)
        elif cmd == "status":
            resp = _cmd_status(args, cfg, short_hash, base)
        elif cmd == "show":
            resp = _cmd_show(args, cfg, short_hash, base)
        elif cmd == "schema":
            resp = _cmd_schema(args, cfg, short_hash, base)
        else:
            raise PolicyViolationError(f"Unsupported command combination: {cmd}")
        return 0, resp.model_dump()
    except OddstarError as exc:
        return exc.exit_code, exc.to_dict()
    except Exception as exc:
        err = OddstarError(str(exc), hint="Run `oddstar doctor --json` to check environment")
        return err.exit_code, err.to_dict()


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    # Quick check for --json before full parse so errors respect JSON mode cleanly
    raw_argv = argv if argv is not None else sys.argv[1:]
    as_json = "--json" in raw_argv
    code, payload = dispatch(raw_argv)
    if code == 0:
        resp = CLIResponse.model_validate(payload)
        _emit(resp, as_json)
    else:
        if as_json:
            print(json.dumps(payload, indent=2), file=sys.stderr)
        else:
            print(f"ERROR [{payload['code']}]: {payload['message']}", file=sys.stderr)
            if payload.get("hint"):
                print(f"Hint: {payload['hint']}", file=sys.stderr)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
