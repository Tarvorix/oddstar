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
from oddstar.errors import ArchiveAccessError, ArtifactNotFoundError, OddstarError, PolicyViolationError
from oddstar.fuse import fuse_results
from oddstar.io import (
    ArchiveFetcher,
    RADIO_CADENCES,
    RawCache,
    TableStore,
    enrich_star_live_for_vetting,
    fetch_gaia_astrometry,
    fetch_live_halo_bulge_sample,
    fetch_single_star_live,
    generate_field_sample,
    get_named_object,
    list_archives_info,
    ra_dec_to_healpix,
)
from oddstar.evidence import export_candidate_evidence, verify_evidence
from oddstar.registry import REGISTRATION_EVENTS, append_entry, clean_source_id, read_ledger
from oddstar.review import (
    add_lesson,
    amend_lesson,
    build_review_packet,
    ingest_all_stored_records,
    ingest_review_record,
    record_lesson_commit_shas,
    register_lesson_file,
)
from oddstar.review.status import model_family, open_findings_by_source, load_review_records, review_status, reviewed_source_ids
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

# Channels whose inputs no live ingestion provides yet, with the data each one still needs. `search --live` refuses
# them: measuring empty inputs for every star would report a search that did not happen. Empty since config 1.5.0:
# every channel has a real-data source, and each reports the stars it could not search.
CHANNELS_WITHOUT_LIVE_SOURCE: dict[str, str] = {}


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
    if args.live and args.channel in CHANNELS_WITHOUT_LIVE_SOURCE:
        # These channels read data no live ingestion provides. A live run would measure empty inputs for every star
        # and report "searched N, flagged 0" without having looked at anything, which is not a null result.
        raise ArchiveAccessError(
            f"Channel '{args.channel}' has no real-data source yet ({CHANNELS_WITHOUT_LIVE_SOURCE[args.channel]}); "
            "a live search would count stars as searched without any data",
            retryable=False,
        )
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
    neowise_fetch: dict[str, Any] | None = None
    tess_fetch: dict[str, Any] | None = None
    spectra_fetch: dict[str, Any] | None = None
    plate_fetch: dict[str, Any] | None = None
    radio_fetch: dict[str, Any] | None = None
    if args.channel == "radio_narrowband" and not is_live:
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
        if args.channel == "optical_dips":
            # QLP light curves for every sample star (cached, rate-limited, retried, resumable). Stars beyond
            # --tess-max-stars, without TESS coverage, or whose fetch failed stay in the sample as not searched.
            from oddstar.io.tess_bulk import enrich_sample_tess

            def _tess_progress(done: int, total: int, sid: str, status: str) -> None:
                if not args.json and (done == total or done % 100 == 0):
                    print(f"  tess {done}/{total} (last {sid}: {status.split(':')[0]})", file=sys.stderr)

            stars, tess_fetch = enrich_sample_tess(stars, cfg, cache, max_stars=args.tess_max_stars, progress=_tess_progress)
        n_fields = len(stars)
        if args.channel == "radio_narrowband":
            # Breakthrough Listen fine-resolution scans for every sample star inside a telescope beam whose ON and
            # OFF files are in the data directory (index answers cached; radio files are never downloaded here).
            # Every other star stays in the sample as not searched, with the reason.
            from oddstar.io.bl_radio import default_data_dir, enrich_sample_radio

            def _radio_progress(done: int, total: int, sid: str, status: str) -> None:
                if not args.json and (done == total or done % 100 == 0 or status == "ok"):
                    print(f"  radio {done}/{total} (last {sid}: {status.split(':')[0]})", file=sys.stderr)

            def _radio_blocks(sid: str, done: int, total: int) -> None:
                if not args.json and (done == total or done % 16 == 0):
                    print(f"  radio {sid}: block {done}/{total}", file=sys.stderr)

            radio_dir = Path(args.radio_data_dir) if args.radio_data_dir else default_data_dir(root, cfg)
            stars, radio_fetch = enrich_sample_radio(stars, cfg, cache, ch, radio_dir, max_stars=args.radio_max_stars,
                                                     progress=_radio_progress, block_progress=_radio_blocks)
        if args.channel == "vanishing":
            # The unit of this search is a plate source, not a sample star (a sample star is in Gaia by
            # construction): every USNO-B1.0 source with two first-epoch plates around each sample position,
            # compared with Gaia and Pan-STARRS (cached, rate-limited, retried, resumable). From here on `stars`
            # holds plate sources; those that could not be compared stay in the list as not searched.
            from oddstar.io.plate_survey import fetch_plate_units

            def _plate_progress(done: int, total: int, sid: str, status: str) -> None:
                if not args.json and (done == total or done % 100 == 0):
                    print(f"  plate fields {done}/{total} (last {sid}: {status.split(':')[0]})", file=sys.stderr)

            stars, plate_fetch = fetch_plate_units(stars, cfg, cache, max_fields=args.plate_max_fields, progress=_plate_progress)
        if args.channel == "spectral_lines":
            # Two LAMOST low-resolution epochs around H-alpha for every sample star that has them (cached,
            # rate-limited, retried, resumable). Stars beyond --spectra-max-stars, without a spectrum, with one
            # usable epoch, or whose fetch failed stay in the sample as not searched.
            from oddstar.io.lamost_bulk import enrich_sample_lamost

            def _spectra_progress(done: int, total: int, sid: str, status: str) -> None:
                if not args.json and (done == total or done % 100 == 0):
                    print(f"  spectra {done}/{total} (last {sid}: {status.split(':')[0]})", file=sys.stderr)

            stars, spectra_fetch = enrich_sample_lamost(stars, cfg, cache, max_stars=args.spectra_max_stars, progress=_spectra_progress)
        if args.channel == "ir_timedomain":
            # NEOWISE-R single-exposure rows only, for every sample star (cached, rate-limited, retried, resumable).
            # Full vetting enrichment still runs below, for flagged stars only. Stars beyond --neowise-max-stars or
            # whose fetch failed stay in the sample and are reported as not searched.
            from oddstar.io.neowise_bulk import enrich_sample_neowise

            def _neowise_progress(done: int, total: int, sid: str, status: str) -> None:
                if not args.json and (done == total or done % 100 == 0):
                    print(f"  neowise {done}/{total} (last {sid}: {status.split(':')[0]})", file=sys.stderr)

            stars, neowise_fetch = enrich_sample_neowise(
                stars, cfg, cache, max_stars=args.neowise_max_stars, progress=_neowise_progress
            )
        nside = cfg.storage.healpix_nside
        cuts = [
            {
                "cut_name": f"gaia_{args.region}_live",
                "description": f"Live Gaia DR3 x AllWISE x 2MASS region='{args.region}' (n={n_fields}, after_random_index={args.after_random_index})",
                "count_before": n_fields,
                "count_after": n_fields,
            }
        ]
        if plate_fetch is not None:
            cuts.append({
                "cut_name": "usno_b1_plate_sources",
                "description": (f"USNO-B1.0 sources within {plate_fetch['field_radius_arcsec']:g}\" of each sample position with both first-epoch "
                                f"plates and no spike flag ({plate_fetch['n_usno_sources']} catalogue sources in {plate_fetch['n_fields_requested']} fields)"),
                "count_before": int(plate_fetch["n_usno_sources"]),
                "count_after": len(stars),
            })
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
        if plate_fetch is not None:
            from oddstar.io.plate_survey import plate_search_status

            plate_status = plate_search_status(star)
            if plate_status != "searched":
                # No modern comparison was made for this plate source, so it can be neither flagged nor cleared
                sc = sc.model_copy(update={"flagged": False, "summary": f"NOT SEARCHED ({plate_status}): plate source not compared with modern catalogues"})
        if radio_fetch is not None:
            from oddstar.io.bl_radio import radio_search_status

            radio_status = radio_search_status(star)
            if radio_status != "searched":
                sc = sc.model_copy(update={"flagged": False, "summary": f"NOT SEARCHED ({radio_status}): no ON/OFF radio scans were compared"})
        # Calibrate FAP immediately against the active certificate
        from oddstar.fuse.engine import calibrate_fap

        cal_fap = calibrate_fap(sc.raw_statistic, cert, sc.false_alarm_prob, float(cfg.fusion.fap_floor))
        sc = sc.model_copy(update={"false_alarm_prob": cal_fap})
        if is_live and args.channel == "ir_excess" and sc.flagged:
            # filtered_live_ingestion: W4 upper limits and contaminated W3/W4 coadds never count as live flags
            # (ph_qual letters and cc_flags bands from `live_ingestion`; WISE band n is character n-1 of the flag strings)
            aw = star.get("allwise", {})
            phq = str(aw.get("ph_qual", "AAAA"))
            ccf = str(aw.get("cc_flags", "0000"))
            w4_upper = len(phq) >= 4 and phq[3] in set(cfg.live_ingestion.ir_excess_w4_upper_limit_ph_qual)
            ir_artifact = len(ccf) >= 4 and any(ccf[int(band) - 1] != "0" for band in cfg.live_ingestion.ir_excess_artifact_cc_flag_bands)
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

    keep_as_searched = args.channel in ("vanishing", "radio_narrowband") or (args.channel == "spectral_lines" and args.defer_enrichment)
    if is_live and flagged_ids and keep_as_searched:
        # The spectral checks read only the spectra the survey fetch attached, the vanishing checks only the plate
        # record and its modern comparison, and the radio checks only the cadence cut from the scans; `vet --live`
        # cannot fetch any of them again for a single object. A plate source is not a Gaia star, so the per-star
        # archive enrichment does not apply to it at all.
        (run_dir / "live_stars.json").write_text(
            json.dumps({fsid: stars_by_id[fsid] for fsid in flagged_ids}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if is_live and flagged_ids and not args.defer_enrichment and args.channel not in ("vanishing", "radio_narrowband"):
        live_flagged: dict[str, Any] = {}
        for fsid in flagged_ids:
            enriched, _ = enrich_star_live_for_vetting(stars_by_id[fsid], cfg, cache)
            live_flagged[fsid] = enriched
        (run_dir / "live_stars.json").write_text(json.dumps(live_flagged, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    coverage: dict[str, Any] | None = None
    if args.channel == "ir_timedomain":
        from oddstar.channels.ir_timedomain.channel import search_coverage_summary

        coverage = search_coverage_summary(stars, cfg)
    elif args.channel == "optical_dips" and is_live:
        from oddstar.channels.optical_dips.channel import search_coverage_summary as dips_coverage_summary

        coverage = dips_coverage_summary(stars, cfg)
    elif args.channel == "spectral_lines" and is_live:
        from oddstar.channels.spectral_lines.channel import search_coverage_summary as spectra_coverage_summary

        coverage = spectra_coverage_summary(stars, cfg)
    elif plate_fetch is not None:
        from oddstar.io.plate_survey import search_coverage_summary as plate_coverage_summary

        coverage = plate_coverage_summary(stars, cfg)
    elif radio_fetch is not None:
        from oddstar.io.bl_radio import search_coverage_summary as radio_coverage_summary

        coverage = radio_coverage_summary(stars, cfg)
    n_searched = coverage["n_searched"] if coverage is not None else len(stars)

    payload = {
        "run_id": run_id,
        "channel": args.channel,
        "gaia_release": cfg.gaia.release,
        "config_hash": short_hash,
        "git_commit": current_git_commit(root),
        "seed": seed,
        "certificate": f"validation/certificates/{args.channel}-{cfg.gaia.release}-{short_hash}.json",
        "n_searched": n_searched,
        "n_flagged": len(flagged_ids),
        "flagged_ids": flagged_ids,
        "cuts": cuts,
        "healpix_cells": healpix_cells,
        "measurements": measurements,
        "scores": scores,
    }
    if coverage is not None:
        # n_input, n_searched, n_not_searched, n_no_coverage, n_failed, n_not_fetched, n_insufficient_visits, not_searched_ids
        payload.update(coverage)
        if neowise_fetch is not None:
            payload["neowise_fetch"] = neowise_fetch
        if tess_fetch is not None:
            payload["tess_fetch"] = tess_fetch
        if spectra_fetch is not None:
            payload["spectra_fetch"] = spectra_fetch
        if plate_fetch is not None:
            payload["plate_fetch"] = plate_fetch
        if radio_fetch is not None:
            payload["radio_fetch"] = radio_fetch
            # Counts and evaluated windows of every searched star, flagged or not
            payload["radio_search"] = {str(s["source_id"]): s["radio_search"] for s in stars if s.get("radio_search")}
    if is_live or coverage is not None:
        # The real searched sample, so the null report is built from these stars and not from regenerated fixtures
        # With the position, so that a follow-up command (radio recordings, SPHEREx images) can find any star of the sample
        payload["sample_stars"] = [
            {"source_id": str(s["source_id"]), "distance_pc": s.get("distance_pc"), "teff_k": s.get("teff_k"),
             **{k: s[k] for k in ("ra_deg", "dec_deg", "pmra_mas_yr", "pmdec_mas_yr") if s.get(k) is not None},
             **({"g_mag": s["mags"]["G"]} if (s.get("mags") or {}).get("G") is not None else {})}
            for s in stars
        ]
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
        summary={"n_searched": n_searched, "n_flagged": len(flagged_ids)},
    )
    (run_dir / "run_state.json").write_text(state.model_dump_json(indent=2) + "\n", encoding="utf-8")

    warnings: list[str] = []
    if coverage is not None and coverage["n_not_searched"]:
        warnings.append(
            f"{coverage['n_not_searched']} of {coverage['n_input']} stars not searched "
            f"(no_coverage={coverage.get('n_no_coverage', 0)}, failed={coverage.get('n_failed', 0)}, not_fetched={coverage.get('n_not_fetched', 0)}, "
            f"insufficient_visits={coverage.get('n_insufficient_visits', 0)}, no_light_curve={coverage.get('n_no_light_curve', 0)}, "
            f"single_epoch={coverage.get('n_single_epoch', 0)}, no_spectrum={coverage.get('n_no_spectrum', 0)}, "
            f"too_faint={coverage.get('n_too_faint', 0)}, near_bright_star={coverage.get('n_near_bright_star', 0)}, not_downloaded={coverage.get('n_not_downloaded', 0)}, "
            f"no_comparison={coverage.get('n_no_comparison', 0)}); "
            "they are excluded from the null result"
        )
    if plate_fetch is not None and plate_fetch["n_fields_failed_usno"]:
        warnings.append(f"{plate_fetch['n_fields_failed_usno']} of {plate_fetch['n_fields_requested']} fields returned no plate catalogue "
                        "(query failed); their plate sources are unknown and are in no count")
    return CLIResponse(
        ok=True,
        summary=f"Searched {n_searched} of {len(stars)} targets on '{args.channel}' -> {len(flagged_ids)} flagged (run_id={run_id})",
        artifacts=[str(search_out)],
        data={
            "run_id": run_id,
            "n_searched": n_searched,
            "n_input": len(stars),
            "n_flagged": len(flagged_ids),
            "flagged_ids": flagged_ids,
            **(
                {"coverage": {k: v for k, v in coverage.items() if k != "not_searched_ids"}, "neowise_fetch": neowise_fetch,
                 "tess_fetch": tess_fetch, "spectra_fetch": spectra_fetch, "plate_fetch": plate_fetch, "radio_fetch": radio_fetch}
                if coverage is not None
                else {}
            ),
        },
        warnings=warnings,
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
                if args.live and cache is not None and star.get("plate_source"):
                    # What the digitized plates show at a plate source's position: about five small cutouts, stored rows
                    from oddstar.io.plate_images import enrich_unit

                    star = enrich_unit(json.loads(json.dumps(star)), cache, cfg=cfg)
            elif args.live and cache is not None and sid not in RADIO_CADENCES:
                star, _ = fetch_single_star_live(sid, cfg, cache)
            elif sid in RADIO_CADENCES:
                star = dict(RADIO_CADENCES[sid])
            else:
                star = resolve_star_record(sid, cfg)
            if args.live and cache is not None and not star.get("plate_source") and sid not in RADIO_CADENCES:
                from oddstar.io.spherex import attach_colours

                # the stored SPHEREx comparison, when `oddstar spherex verdict` made one for this star
                star = attach_colours(json.loads(json.dumps(star)), cache)
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
    meas_by_star: dict[str, list[dict[str, Any]]] = {}
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
            meas_rows = {str(m.get("source_id")): m for m in sdata.get("measurements", [])}
            for row in sdata["scores"]:
                sc = ChannelScore.model_validate(row)
                if sc.flagged:
                    scores_by_star.setdefault(sc.source_id, []).append(sc)
                    if sc.source_id in meas_rows:
                        meas_by_star.setdefault(sc.source_id, []).append(meas_rows[sc.source_id])
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

    (out_dir / "fuse_inputs.json").write_text(json.dumps({"run_id": out_run_id, "input_run_ids": list(args.run_ids)}, indent=2) + "\n", encoding="utf-8")

    # Append fusion outcomes to the registry ledger, each with its committed evidence directory
    ledger = root / cfg.storage.registry_dir / "ledger.jsonl"
    now = datetime.now(timezone.utc).isoformat()
    evidence_dirs: list[str] = []
    for cand in fusion.candidates:
        ev_type = "explained_anomaly" if cand.tier == "Explained" else "candidate_registered"
        exp = _export_for_registration(
            root, cfg, short_hash, source_id=cand.source_id, run_id=out_run_id, tier=cand.tier,
            vetting=vet_by_star.get(cand.source_id, []), measurements=meas_by_star.get(cand.source_id, []),
            scores=scores_by_star.get(cand.source_id, []),
        )
        evidence_dirs.append(exp.rel_dir)
        entry = RegistryEntry(
            entry_id=f"REG-FUSE-{out_run_id}-{cand.source_id}",
            timestamp_utc=now,
            source_id=cand.source_id,
            event_type=ev_type,
            tier=cand.tier,
            status=cand.explained_reason or f"Tier {cand.tier}",
            provenance={"run_id": out_run_id, "config_hash": short_hash, "git_commit": current_git_commit(root)},
            details={**cand.model_dump(), "vetting": exp.vetting_rows, "evidence_scope": exp.index.scope},
            evidence_dir=exp.rel_dir,
            evidence_index_sha256=exp.index_sha256,
        )
        append_entry(root, entry, cfg.storage.registry_dir)

    tier_counts: dict[str, int] = {}
    for c in fusion.candidates:
        tier_counts[c.tier] = tier_counts.get(c.tier, 0) + 1

    return CLIResponse(
        ok=True,
        summary=f"Fused {len(fusion.candidates)} flagged targets across {len(certs_by_ch)} channels (tiers={tier_counts})",
        artifacts=[str(out_path), str(ledger), *evidence_dirs],
        data={"fusion": fusion.model_dump(), "tier_counts": tier_counts, "evidence_dirs": evidence_dirs},
        next=f"oddstar null-report {args.run_ids[0]} --json",
    )


def _export_for_registration(
    root: Path,
    cfg: PipelineConfig,
    short_hash: str,
    *,
    source_id: str,
    run_id: str,
    tier: str,
    vetting: list[VettingResult],
    measurements: list[dict[str, Any]],
    scores: list[ChannelScore],
    images: list[Path] | None = None,
) -> Any:
    """Committed evidence for a registration (all checks for A/B/C/Explained; failing cause only for Rejected)."""
    return export_candidate_evidence(
        root=root,
        cfg=cfg,
        source_id=source_id,
        run_id=run_id,
        tier=tier,
        vetting=vetting,
        measurements=measurements,
        scores=scores,
        images=images or [],
        config_hash=short_hash,
        config_hash_full=compute_config_hash(cfg.model_dump()),
        git_commit=current_git_commit(root),
        pipeline_version=__version__,
    )


def _inspect_single_target(
    source_id: str,
    cfg: PipelineConfig,
    short_hash: str,
    root: Path,
    *,
    run_id: str | None = None,
    live: bool = False,
    register: bool = False,
) -> tuple[Any, Path, Path]:
    clean_id = source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
    rid = run_id or f"inspect-{clean_id}-{cfg.gaia.release}-{short_hash}"
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)

    if live:
        star, live_mans = fetch_single_star_live(clean_id, cfg, cache)
        from oddstar.io.spherex import attach_colours

        star = attach_colours(star, cache)  # the stored SPHEREx comparison, when `oddstar spherex verdict` made one
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
                    model_family=model_family(rr.reviewer_model),
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
        review_status=review_status(root, clean_id, registry_dir=cfg.storage.registry_dir) if reviews else None,
    )
    if register:
        _register_inspection(root, cfg, short_hash, dossier, json_p, measurements, scores, vetting)
    return dossier, json_p, md_p


def _register_inspection(
    root: Path,
    cfg: PipelineConfig,
    short_hash: str,
    dossier: Any,
    json_p: Path,
    measurements: list[Measurement],
    scores: list[ChannelScore],
    vetting: list[VettingResult],
) -> None:
    """`inspect --register`: append the dossier tier to the ledger with its committed evidence directory."""
    sid = clean_source_id(dossier.source_id)
    rid = dossier.provenance.run_id
    flagged = {c.name for c in dossier.channels if c.flagged}
    images = [json_p.parent / ev for c in dossier.channels if c.name in flagged for ev in c.evidence]
    exp = _export_for_registration(
        root, cfg, short_hash, source_id=sid, run_id=rid, tier=dossier.tier, vetting=vetting,
        measurements=[m.model_dump() for m in measurements], scores=scores, images=images,
    )
    holes = sorted({v.check for v in vetting if v.status == "not_evaluated"})
    append_entry(
        root,
        RegistryEntry(
            entry_id=f"REG-INSPECT-{rid}-{sid}-{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')}",
            timestamp_utc=datetime.now(timezone.utc).isoformat(),
            source_id=sid,
            event_type="explained_anomaly" if dossier.tier == "Explained" else "candidate_registered",
            tier=dossier.tier,
            status=f"Tier {dossier.tier} (inspect)",
            provenance={"run_id": rid, "config_hash": short_hash, "git_commit": current_git_commit(root), "mode": "inspect"},
            details={
                "source_id": sid,
                "tier": dossier.tier,
                "flagged_channels": sorted(flagged),
                "open_holes": holes,
                "vetting": exp.vetting_rows,
                "evidence_scope": exp.index.scope,
            },
            evidence_dir=exp.rel_dir,
            evidence_index_sha256=exp.index_sha256,
        ),
        cfg.storage.registry_dir,
    )


def _cmd_inspect(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would inspect {args.source_id} across all channels and write dossier")
    dossier, json_p, md_p = _inspect_single_target(
        args.source_id, cfg, short_hash, root, run_id=args.run_id, live=bool(args.live), register=bool(args.register)
    )
    flagged = [c.name for c in dossier.channels if c.flagged]
    holes = [v.check for v in dossier.vetting if v.status == "not_evaluated"]
    return CLIResponse(
        ok=True,
        summary=f"Inspected {dossier.source_id} -> Tier {dossier.tier} (flagged={flagged}, open_holes={holes})",
        artifacts=[str(json_p), str(md_p)],
        data={
            "tier": dossier.tier,
            "flagged_channels": flagged,
            "open_holes": holes,
            "lead_lines": dossier.lead_lines,
            "registered": bool(args.register),
            "review_status": dossier.review_status.summary if dossier.review_status else None,
        },
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
    if args.all:
        if args.dry_run:
            return CLIResponse(ok=True, summary="[dry-run] Would ingest every stored record under review/records/")
        results = ingest_all_stored_records(cfg, root)
        appended = [e for r in results for e in r["appended_entries"]]
        return CLIResponse(
            ok=True,
            summary=f"Ingested {len(results)} stored review records ({sum(not r['already_ingested'] for r in results)} new; {len(appended)} ledger events appended)",
            artifacts=[r["stored_path"] for r in results],
            data={"results": results, "appended_entries": appended},
            next="oddstar status --json",
        )
    if not args.file:
        raise PolicyViolationError("review ingest needs a record file or --all")
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


def _cmd_lessons_register(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    path = Path(args.file)
    if not path.is_absolute():
        path = root / path
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would register lesson file {path}")
    rec, appended = register_lesson_file(path, cfg, root)
    return CLIResponse(
        ok=True,
        summary=f"Registered lesson {rec.lesson_id} ({len(appended)} ledger events appended)",
        artifacts=[str(path)],
        data={"lesson": rec.model_dump(), "appended_entries": appended},
        next="oddstar status --json",
    )


def _cmd_lessons_amend(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    updates: dict[str, Any] = {}
    for item in args.set or []:
        key, _, raw = item.partition("=")
        try:
            updates[key.strip()] = json.loads(raw)
        except json.JSONDecodeError:
            updates[key.strip()] = raw
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would amend lesson {args.lesson_id}: {sorted(updates)}")
    rec, appended = amend_lesson(args.lesson_id, updates, reason=args.reason, cfg=cfg, root=root, md_note=args.md_note)
    return CLIResponse(
        ok=True,
        summary=f"Amended lesson {rec.lesson_id} ({len(appended)} ledger events appended)",
        artifacts=[str(root / "lessons" / f"{rec.lesson_id}.json")],
        data={"lesson": rec.model_dump(), "appended_entries": appended},
    )


def _cmd_lessons_record_sha(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    if args.dry_run:
        return CLIResponse(ok=True, summary="[dry-run] Would set each lesson commit_sha to the commit that introduced its files")
    res = record_lesson_commit_shas(cfg, root)
    changed = [r for r in res if r["action"] == "amended"]
    pending = [r["lesson_id"] for r in res if r["action"] == "not committed yet"]
    return CLIResponse(
        ok=True,
        summary=f"Lesson commit_sha: {len(changed)} amended, {len(pending)} not committed yet",
        data={"results": res},
        warnings=[f"{lid}: lesson files not committed yet; run again after the commit that adds them" for lid in pending],
    )


def _cmd_evidence_export(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    sid = clean_source_id(args.source_id)
    vetting, measurements, scores, tier, run_id, images = _collect_candidate_artifacts(root, cfg, sid, args.run_id)
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would export evidence for {sid} from run {run_id} (tier {tier}, {len(vetting)} checks)")
    exp = _export_for_registration(
        root, cfg, short_hash, source_id=sid, run_id=run_id, tier=tier, vetting=vetting, measurements=measurements,
        scores=scores, images=images,
    )
    return CLIResponse(
        ok=True,
        summary=f"Exported evidence for {sid} (run {run_id}, tier {tier}, {len(exp.index.files)} files) to {exp.rel_dir}",
        artifacts=[str(exp.directory / "index.json")],
        data={"evidence_dir": exp.rel_dir, "index_sha256": exp.index_sha256, "index": exp.index.model_dump(mode="json")},
        warnings=["Export only: no ledger entry was written (fuse and inspect --register register candidates)"],
        next=f"oddstar evidence verify {sid} --json",
    )


def _cmd_report(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    from oddstar.report import stale_documents, write_documents

    if args.report_cmd == "check" or args.dry_run:
        stale = stale_documents(root, cfg.storage.registry_dir)
        if stale and args.report_cmd == "check":
            raise PolicyViolationError(
                "Generated documents differ from a fresh build: " + ", ".join(stale),
                hint="Run `oddstar report build` and commit the result; edit commentary only in docs/notes/",
            )
        return CLIResponse(ok=True, summary=f"{len(stale)} generated document(s) would change" if stale else "Generated documents are current",
                           data={"stale": stale})
    changed = write_documents(root, cfg.storage.registry_dir)
    return CLIResponse(ok=True, summary=f"Rebuilt generated documents ({len(changed)} changed)", artifacts=changed, data={"changed": changed})


def _cmd_runs_export(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    from oddstar.run_records import build_run_record, export_run_record

    if args.dry_run:
        for rid in args.run_ids:
            build_run_record(root / cfg.storage.runs_dir / rid)
        return CLIResponse(ok=True, summary=f"[dry-run] Would export {len(args.run_ids)} run record(s)")
    written: list[str] = []
    rows: list[dict[str, Any]] = []
    for rid in args.run_ids:
        out = export_run_record(root, rid, runs_dir=cfg.storage.runs_dir, registry_dir=cfg.storage.registry_dir, replace=args.replace)
        written.extend(out["written"])
        rec = out["record"]
        rows.append({
            "run_id": rid,
            "record_path": out["record_path"],
            "record_sha256": out["record_sha256"],
            "searches": [{"channel": s["channel"], "n_sample": s["sample"]["n"], "n_searched": s["n_searched"], "n_flagged": s["n_flagged"]}
                         for s in rec.get("searches", [])],
            "fusion_tier_counts": (rec.get("fusion") or {}).get("tier_counts"),
        })
    return CLIResponse(
        ok=True,
        summary=f"Exported {len(rows)} run record(s) ({len(written)} file(s) written) under {cfg.storage.registry_dir}/",
        artifacts=sorted(set(written)),
        data={"records": rows},
        next="oddstar runs verify --json",
    )


def _cmd_runs_verify(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    from oddstar.run_records import verify_run_records

    rep = verify_run_records(root, cfg.storage.registry_dir)
    if rep["problems"]:
        raise PolicyViolationError(
            "Run record verification failed: " + "; ".join(rep["problems"][:20]),
            hint="Re-export the run with `oddstar runs export <run_id> --replace` from the machine that holds runs/",
        )
    return CLIResponse(ok=True, summary=f"{rep['n_records']} run record(s) and {rep['n_samples']} sample file(s) verified", data=rep)


def _cmd_evidence_verify(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    rep = verify_evidence(root, cfg.storage, args.source_id)
    scope = clean_source_id(args.source_id) if args.source_id else "all sources"
    if rep["problems"]:
        raise PolicyViolationError(
            f"Evidence verification failed for {scope}: " + "; ".join(rep["problems"][:20]),
            hint="Evidence files are immutable once registered; re-export to a new run instead of editing them",
        )
    return CLIResponse(
        ok=True,
        summary=(
            f"Evidence verified for {scope}: {len(rep['verified_dirs'])} directories, "
            f"{rep['ledger_entries_with_evidence']} ledger entries with evidence, "
            f"{rep['legacy_without_evidence']} legacy_without_evidence (registered before the evidence store; no committed evidence)"
        ),
        data={k: v for k, v in rep.items() if k != "legacy_entry_ids"},
    )


def _collect_candidate_artifacts(
    root: Path, cfg: PipelineConfig, sid: str, run_id: str | None
) -> tuple[list[VettingResult], list[dict[str, Any]], list[ChannelScore], str, str, list[Path]]:
    """Vetting rows, measurements, scores, tier, run id and images of a candidate from run or dossier artefacts.

    `--run-id` may name a vet/search run, a fuse run (its `fuse_inputs.json` lists the input runs) or an inspect run
    (the dossier). Without it the source's dossier is used.
    """
    from oddstar.fuse.engine import assign_tier
    from oddstar.schemas import Dossier

    runs_dir = root / cfg.storage.runs_dir
    run_dirs: list[Path] = []
    if run_id and (runs_dir / run_id).is_dir():
        fi = runs_dir / run_id / "fuse_inputs.json"
        if fi.exists():
            run_dirs = [runs_dir / r for r in json.loads(fi.read_text(encoding="utf-8")).get("input_run_ids", [])]
        else:
            run_dirs = [runs_dir / run_id]
    if run_dirs:
        vetting: list[VettingResult] = []
        meas: list[dict[str, Any]] = []
        scores: list[ChannelScore] = []
        for rd in run_dirs:
            for sf in sorted(rd.glob("search_*.json")):
                sdata = json.loads(sf.read_text(encoding="utf-8"))
                meas += [m for m in sdata.get("measurements", []) if str(m.get("source_id")) == sid]
                scores += [ChannelScore.model_validate(r) for r in sdata.get("scores", []) if str(r.get("source_id")) == sid]
            vf = rd / "vetting.json"
            if vf.exists():
                vetting += [VettingResult.model_validate(v) for v in json.loads(vf.read_text(encoding="utf-8")).get("vetting", {}).get(sid, [])]
        if not scores:
            raise ArtifactNotFoundError(f"No scores for {sid} in run(s) {[d.name for d in run_dirs]}")
        tier, _ = assign_tier(scores, vetting, cfg)
        return vetting, meas, scores, tier, str(run_id), []
    d_json = root / cfg.storage.dossiers_dir / sid / "dossier.json"
    if not d_json.exists():
        raise ArtifactNotFoundError(f"No run '{run_id}' and no dossier for {sid}", hint=f"Run `oddstar inspect {sid}` or pass --run-id")
    dossier = Dossier.model_validate_json(d_json.read_text(encoding="utf-8"))
    if run_id and dossier.provenance.run_id != run_id:
        raise ArtifactNotFoundError(f"Dossier for {sid} belongs to run {dossier.provenance.run_id}, not {run_id}")
    scores = [
        ChannelScore(source_id=sid, channel=c.name, raw_statistic=float(c.raw_statistic or 0.0), false_alarm_prob=float(c.false_alarm_prob or 1.0),
                     flagged=c.flagged, rank_score=0.0, summary="from dossier")
        for c in dossier.channels
    ]
    meas = [{"channel": c.name, "values": c.measurements, "epochs": c.epochs, "caveats": c.caveats} for c in dossier.channels]
    images = [d_json.parent / ev for c in dossier.channels if c.flagged for ev in c.evidence]
    return list(dossier.vetting), meas, scores, dossier.tier, dossier.provenance.run_id, images


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
        closes_finding=not bool(args.no_close),
        related_finding_ids=list(args.related_finding or []),
        finding_source=args.finding_source or "",
        additional_linked_tests=list(args.also_test or []),
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

    if "sample_stars" in sdata:
        stars = sdata["sample_stars"]  # the real searched sample (live or offline), with its not-searched stars
    else:
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
        not_searched=sdata.get("not_searched_ids") if "sample_stars" in sdata else None,
    )
    return CLIResponse(
        ok=True,
        summary=f"Null report emitted for '{ch_name}' ({rep.n_searched} searched, {rep.n_flagged} flagged, {rep.expected_false_alarms:.2f} expected FA)",
        artifacts=[str(jp), str(mp)],
        data={"null_report": rep.model_dump(exclude={"completeness_grid"})},
    )


_TIER_RANK = {"A": 3, "B": 2, "C": 1}


def _stars_and_tiers_of_runs(run_ids: list[str], cfg: PipelineConfig, root: Path) -> tuple[dict[str, dict[str, Any]], dict[str, str], list[str]]:
    """What earlier runs hold about their stars, for a follow-up command.

    Returns the records that carry a position (from `sample_stars` of the search outputs and from `live_stars.json`,
    which holds the flagged stars in full), the best tier of each flagged star (C when the run was not vetted), and
    the ids of every star of the samples, with or without a position.
    """
    from oddstar.fuse.engine import assign_tier

    stars: dict[str, dict[str, Any]] = {}
    tiers: dict[str, str] = {}
    all_ids: dict[str, None] = {}

    def _keep(record: dict[str, Any], sid: str) -> None:
        if record.get("ra_deg") is None or record.get("dec_deg") is None:
            return
        rec = {k: record[k] for k in ("ra_deg", "dec_deg", "pmra_mas_yr", "pmdec_mas_yr") if record.get(k) is not None}
        g = record.get("g_mag") if record.get("g_mag") is not None else (record.get("mags") or {}).get("G")
        stars.setdefault(sid, {"source_id": sid, **rec, "g_mag": g})

    for run_id in run_ids:
        run_dir = root / cfg.storage.runs_dir / run_id
        search_files = sorted(run_dir.glob("search_*.json"))
        if not search_files:
            raise ArtifactNotFoundError(f"No search outputs found in {run_dir}")
        vf = run_dir / "vetting.json"
        vdata = json.loads(vf.read_text(encoding="utf-8")).get("vetting", {}) if vf.exists() else None
        lf = run_dir / "live_stars.json"
        if lf.exists():
            for sid, record in json.loads(lf.read_text(encoding="utf-8")).items():
                _keep(record, str(sid))
        for sf in search_files:
            sdata = json.loads(sf.read_text(encoding="utf-8"))
            for star in sdata.get("sample_stars", []):
                all_ids[str(star["source_id"])] = None
                _keep(star, str(star["source_id"]))
            scores_by_id = {r["source_id"]: ChannelScore.model_validate(r) for r in sdata["scores"]}
            for sid in sdata.get("flagged_ids", []):
                all_ids[str(sid)] = None
                if vdata is None:
                    tier = "C"  # flagged and not vetted yet: nothing has explained it
                else:
                    tier, _ = assign_tier([scores_by_id[sid]], [VettingResult.model_validate(v) for v in vdata.get(sid, [])], cfg)
                if _TIER_RANK.get(tier, 0) > _TIER_RANK.get(tiers.get(str(sid), ""), 0):
                    tiers[str(sid)] = tier
    return stars, tiers, list(all_ids)


def _positions_for(chosen: list[str], stars: dict[str, dict[str, Any]], run_ids: list[str], cfg: PipelineConfig,
                   cache: RawCache) -> list[str]:
    """Fill `stars` for the chosen ids the runs hold no position for, from Gaia by source id. Returns the ids looked up.

    Search outputs written before positions were recorded hold a position for flagged stars only.
    """
    missing = [sid for sid in chosen if sid not in stars]
    looked_up = [sid for sid in missing if sid.isdigit()]
    if looked_up:
        stars.update(fetch_gaia_astrometry(looked_up, cfg, cache))
    unknown = [sid for sid in chosen if sid not in stars]
    if unknown:
        raise ArtifactNotFoundError(f"No position for {unknown[:5]} in the samples of {run_ids} or in Gaia",
                                    hint="Name a Gaia source id, or a run whose search sample contains the star")
    return looked_up


def _cmd_radio_followup(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    """Search the Breakthrough Listen recordings of chosen stars: the expensive half of the radio channel, on request.

    `search radio_narrowband --live` only records which stars lie in a telescope beam and reads results that are
    already stored. This command fetches and searches the recordings, one file at a time, for: the stars of the
    given runs that vetting left at `--min-tier` or better (default), the stars named with `--ids`, or every star
    of the runs' samples (`--all-covered`). A star outside every beam is recorded as "no_coverage".
    """
    from oddstar.io.bl_radio import default_data_dir, followup_star

    stars, tiers, all_ids = _stars_and_tiers_of_runs(args.from_run, cfg, root)
    if args.ids:
        chosen = [i.strip() for i in args.ids.split(",") if i.strip()]
        selection = "ids"
    elif args.all_covered:
        chosen = sorted(all_ids)
        selection = "all_covered"
    else:
        chosen = sorted(sid for sid, tier in tiers.items() if _TIER_RANK[tier] >= _TIER_RANK[args.min_tier])
        selection = f"tier {args.min_tier} or better"

    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    _positions_for(chosen, stars, args.from_run, cfg, cache)
    positions = {sid: (float(stars[sid]["ra_deg"]), float(stars[sid]["dec_deg"])) for sid in chosen}
    channel = build_channel("radio_narrowband", cfg, root)
    data_dir = Path(args.radio_data_dir) if args.radio_data_dir else default_data_dir(root, cfg)
    records = []
    for n, sid in enumerate(chosen, 1):
        ra, dec = positions[sid]
        rec = followup_star(sid, ra, dec, cfg, cache, channel, data_dir, all_scan_sets=args.scan_sets == "all",
                            keep_windows=args.keep_windows, min_free_bytes=int(args.min_free_gb * 1e9), dry_run=args.dry_run,
                            log=lambda msg: print(f"  {msg}", file=sys.stderr, flush=True))
        rec["tier_before"] = tiers.get(sid)
        records.append(rec)
        print(f"  radio follow-up {n}/{len(chosen)}: {sid}: {rec['status'].split(':')[0]}", file=sys.stderr, flush=True)
    counts: dict[str, int] = {}
    for rec in records:
        key = rec["status"].split(":")[0]
        counts[key] = counts.get(key, 0) + 1
    payload = {
        "run_id": args.run_id or f"radio-followup-{args.from_run[0]}", "from_runs": list(args.from_run), "selection": selection,
        "config_hash": short_hash, "git_commit": current_git_commit(root), "dry_run": bool(args.dry_run),
        "n_selected": len(chosen), "status_counts": counts,
        "n_windows_flagged_outside_rfi_bands": sum(r.get("n_windows_flagged_outside_rfi_bands", 0) for r in records),
        "stars": records,
    }
    artifacts: list[str] = []
    if not args.dry_run:
        out_dir = root / cfg.storage.runs_dir / payload["run_id"]
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / "radio_followup.json"
        out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        artifacts.append(str(out))
    return CLIResponse(
        ok=True,
        summary=(("[dry-run] " if args.dry_run else "") + f"Radio follow-up of {len(chosen)} star(s) ({selection}): "
                 + (", ".join(f"{v} {k}" for k, v in sorted(counts.items())) or "nothing selected")),
        artifacts=artifacts,
        data=payload,
        next="oddstar search radio_narrowband --live --json" if counts.get("searched") else "",
    )


def _cmd_radio_beams(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    """Test the flagged radio windows of chosen stars against the other beams recorded at the same instant.

    For a scan of a multibeam receiver (Parkes multibeam: 13 beams) the other beams point at other parts of the sky.
    Every window of a star's stored search that is flagged outside the allocated bands is read from each of them,
    a few kilohertz per beam by HTTP range requests: no scan file is downloaded. The result is stored with the
    search; the next `search radio_narrowband --live` and `vet` then evaluate the check `multi_target`.
    Stars: those the given runs left at `--min-tier` or better (default), or the ones named with `--ids`.
    """
    from oddstar.io.bl_radio import beam_test_star

    stars, tiers, _ = _stars_and_tiers_of_runs(args.from_run, cfg, root)
    if args.ids:
        chosen = [i.strip() for i in args.ids.split(",") if i.strip()]
        selection = "ids"
    else:
        chosen = sorted(sid for sid, tier in tiers.items() if _TIER_RANK[tier] >= _TIER_RANK[args.min_tier])
        selection = f"tier {args.min_tier} or better"
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    _positions_for(chosen, stars, args.from_run, cfg, cache)
    channel = build_channel("radio_narrowband", cfg, root)
    records = []
    for n, sid in enumerate(chosen, 1):
        rec = beam_test_star(sid, float(stars[sid]["ra_deg"]), float(stars[sid]["dec_deg"]), cfg, cache, channel, redo=args.redo,
                             log=lambda msg, _sid=sid: print(f"  {_sid}: {msg}", file=sys.stderr, flush=True))
        rec["tier_before"] = tiers.get(sid)
        records.append(rec)
        print(f"  beam test {n}/{len(chosen)}: {sid}: {rec['status'].split(':')[0]}", file=sys.stderr, flush=True)
    counts: dict[str, int] = {}
    for rec in records:
        key = rec["status"].split(":")[0]
        counts[key] = counts.get(key, 0) + 1
    totals = {k: sum(int(r.get(k, 0)) for r in records)
              for k in ("n_windows_tested", "n_windows_in_far_beams", "n_windows_in_no_far_beam", "n_windows_not_judged")}
    payload = {
        "run_id": args.run_id or f"radio-beams-{args.from_run[0]}", "from_runs": list(args.from_run), "selection": selection,
        "config_hash": short_hash, "git_commit": current_git_commit(root), "n_selected": len(chosen), "status_counts": counts, **totals,
        "stars": records,
    }
    out_dir = root / cfg.storage.runs_dir / payload["run_id"]
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "radio_beams.json"
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return CLIResponse(
        ok=True,
        summary=(f"Beam test of {len(chosen)} star(s) ({selection}): " + (", ".join(f"{v} {k}" for k, v in sorted(counts.items())) or "nothing selected")
                 + f"; {totals['n_windows_tested']} window(s) read: {totals['n_windows_in_far_beams']} also in a far beam, "
                 f"{totals['n_windows_in_no_far_beam']} in no far beam, {totals['n_windows_not_judged']} not judged"),
        artifacts=[str(out)],
        data=payload,
        next="oddstar search radio_narrowband --live --json" if counts.get("tested") else "",
    )


def _cmd_spherex_photometry(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    """Measure chosen stars on every SPHEREx spectral image that holds them (0.75 to 5 micron, 2025 on).

    For the stars of the given runs that vetting left at `--min-tier` or better (default), or the stars named with
    `--ids`. Each exposure is one cutout of about 4 MB (not kept) and one stored row; a star has a few hundred.
    Writes one file per star with every measurement, and a summary of all stars.
    """
    from oddstar.io import spherex as sx

    stars, tiers, _ = _stars_and_tiers_of_runs(args.from_run, cfg, root)
    if args.ids:
        chosen = [i.strip() for i in args.ids.split(",") if i.strip()]
        selection = "ids"
    else:
        chosen = sorted(sid for sid, tier in tiers.items() if _TIER_RANK[tier] >= _TIER_RANK[args.min_tier])
        selection = f"tier {args.min_tier} or better"
    detectors = None
    if args.detectors:
        parts = [d.strip() for d in args.detectors.split(",") if d.strip()]
        if not parts or any(not d.isdigit() or not 1 <= int(d) <= 6 for d in parts):
            raise OddstarError(f"--detectors takes numbers 1 to 6, got {args.detectors!r}", code="USAGE_ERROR",
                               hint="Example: --detectors 3,4,5,6")
        detectors = tuple(sorted({int(d) for d in parts}))
    as_of = args.as_of or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    run_name = args.run_id or f"spherex-{args.from_run[0]}"
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would measure {len(chosen)} star(s) ({selection}) on SPHEREx images as of {as_of}",
                           data={"run_id": run_name, "selection": selection, "source_ids": chosen, "as_of": as_of})

    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)
    from_gaia = _positions_for(chosen, stars, args.from_run, cfg, cache)
    out_dir = root / cfg.storage.runs_dir / run_name
    (out_dir / "spherex").mkdir(parents=True, exist_ok=True)
    records, artifacts = [], []
    for n, sid in enumerate(chosen, 1):
        star = stars[sid]
        ra, dec = float(star["ra_deg"]), float(star["dec_deg"])
        rec: dict[str, Any] = {"source_id": sid, "tier_before": tiers.get(sid), "position_from": "gaia" if sid in from_gaia else "run"}
        try:
            neighbours = () if args.no_neighbours else sx.gaia_neighbours(
                sid, ra, dec, star.get("g_mag"), cache, release=cfg.gaia.release)
            res = sx.forced_photometry(
                sid, ra, dec, cache, as_of=as_of, pmra_mas_yr=float(star.get("pmra_mas_yr") or 0.0),
                pmdec_mas_yr=float(star.get("pmdec_mas_yr") or 0.0), neighbours=neighbours, detectors=detectors, workers=args.workers,
                log=lambda msg: print(f"  {msg}", file=sys.stderr, flush=True))
            path = out_dir / "spherex" / f"{sid}.json"
            path.write_text(json.dumps(res, indent=1, sort_keys=True) + "\n", encoding="utf-8")
            artifacts.append(str(path))
            rec.update({k: res[k] for k in ("n_images", "n_measured", "n_failed", "complete", "status_counts", "errors", "collections",
                                            "n_neighbours", "summary")})
            rec["status"] = "no_coverage" if res["n_images"] == 0 else ("measured" if res["complete"] else "incomplete")
        except Exception as exc:  # one star's archive failure must not stop the others; it is recorded with its cause
            rec["status"] = f"failed: {type(exc).__name__}: {exc}"[:400]
        records.append(rec)
        print(f"  SPHEREx {n}/{len(chosen)}: {sid}: {rec['status'].split(':')[0]}", file=sys.stderr, flush=True)
    counts: dict[str, int] = {}
    for rec in records:
        key = rec["status"].split(":")[0]
        counts[key] = counts.get(key, 0) + 1
    payload = {"run_id": run_name, "from_runs": list(args.from_run), "selection": selection, "as_of": as_of, "detectors": detectors,
               "photometry_version": sx.PHOTOMETRY_VERSION, "config_hash": short_hash, "git_commit": current_git_commit(root),
               "n_selected": len(chosen), "status_counts": counts, "stars": records}
    out = out_dir / "spherex_photometry.json"
    out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return CLIResponse(
        ok=True,
        summary=f"SPHEREx photometry of {len(chosen)} star(s) ({selection}): " + (", ".join(f"{v} {k}" for k, v in sorted(counts.items())) or "nothing selected"),
        artifacts=[str(out)] + artifacts,
        data=payload,
    )


def _cmd_spherex_verdict(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    """Measure stars and their comparison stars on the SPHEREx images and apply the colour rule of `docs/spherex.md`.

    For each id: its catalogue record (Gaia with AllWISE and 2MASS), the nearest comparison stars of about its colour
    and brightness (`vetting.spherex_comparison_*`), forced photometry of all of them on detectors 3 to 6 (a few
    hundred cutouts per star; stored rows are not fetched again), and the verdict. The measured colours are stored,
    and `vet --live` and `inspect --live` read them for the `spherex_colours` check.
    """
    from oddstar.io import spherex as sx
    from oddstar.vet.ir_checks import spherex_colour_verdict

    ids = [i.strip().replace("Gaia DR3 ", "") for i in (args.ids or "").split(",") if i.strip()]
    if not ids or any(not i.isdigit() for i in ids):
        raise OddstarError(f"--ids takes Gaia source ids, got {args.ids!r}", code="USAGE_ERROR", hint="Example: --ids 5628404318169472")
    as_of = args.as_of or datetime.now(timezone.utc).strftime("%Y-%m-%d")
    run_name = args.run_id or f"spherex-verdict-{ids[0]}"
    n_comp = int(cfg.vetting["spherex_n_comparison_stars"])
    if args.dry_run:
        return CLIResponse(ok=True, summary=f"[dry-run] Would measure {len(ids)} star(s) and up to {n_comp} comparison stars each on SPHEREx "
                                            f"images as of {as_of}", data={"run_id": run_name, "source_ids": ids, "as_of": as_of})
    cache = RawCache(root / cfg.storage.raw_dir, max_cache_gb=cfg.storage.max_cache_gb)

    def measure(star: dict[str, Any]) -> dict[str, Any]:
        sid, ra, dec = str(star["source_id"]), float(star["ra_deg"]), float(star["dec_deg"])
        neighbours = sx.gaia_neighbours(sid, ra, dec, (star.get("mags") or {}).get("G"), cache, release=cfg.gaia.release)
        res = sx.forced_photometry(
            sid, ra, dec, cache, as_of=as_of, pmra_mas_yr=float(star.get("pmra_mas_yr") or 0.0),
            pmdec_mas_yr=float(star.get("pmdec_mas_yr") or 0.0), neighbours=neighbours, detectors=sx.COLOUR_DETECTORS, workers=args.workers,
            log=lambda msg: print(f"  {msg}", file=sys.stderr, flush=True))
        return sx.colour_record(star, res)

    def measure_whole(star: dict[str, Any], *, tolerate: bool) -> dict[str, Any]:
        """Measure a star, and again while images are missing: a later pass fetches only what the archive refused before."""
        rec: dict[str, Any] = {}
        for attempt in range(1, passes + 1):
            try:
                rec = measure(star)
            except Exception as exc:  # one comparison star's archive failure is recorded; it then does not count
                if not tolerate:
                    raise
                rec = {"source_id": str(star["source_id"]), "complete": False, "error": f"{type(exc).__name__}: {exc}"[:300]}
            rec["passes"] = attempt
            if rec.get("complete") or rec.get("n_images") == 0 or attempt == passes:
                break
            print(f"    {rec['source_id']}: {rec.get('n_failed', 'all')} image(s) not fetched in pass {attempt}; pass {attempt + 1} asks for those only",
                  file=sys.stderr, flush=True)
        return rec

    passes = max(1, int(args.passes))
    records: list[dict[str, Any]] = []
    for n, sid in enumerate(ids, 1):
        rec: dict[str, Any] = {"source_id": sid}
        try:
            target = sx.target_star(sid, cfg, cache)
            if target is None:
                rec["status"] = "no_catalogue_record"
            else:
                sel = sx.comparison_stars(target, cfg, cache)
                print(f"  SPHEREx verdict {n}/{len(ids)}: {sid}: {len(sel['stars'])} comparison star(s), the nearest of {sel['n_clean']} clean "
                      f"among {sel['n_examined']} Gaia source(s) examined out to {sel['examined_out_to_deg']:g} deg", file=sys.stderr, flush=True)
                t_rec = measure_whole(target, tolerate=False)
                comparison: list[dict[str, Any]] = []
                for k, st in enumerate(sel["stars"], 1):
                    c_rec = measure_whole(st, tolerate=True)
                    c_rec["comparison_sep_arcsec"] = st.get("comparison_sep_arcsec")
                    comparison.append(c_rec)
                    state = "complete" if c_rec.get("complete") else f"incomplete ({c_rec.get('n_failed', 'all')} image(s) not fetched)"
                    print(f"    comparison {k}/{len(sel['stars'])}: {c_rec['source_id']}: {state}", file=sys.stderr, flush=True)
                payload = {"target": t_rec, "comparison": comparison, "selection": {k: v for k, v in sel.items() if k != "stars"},
                           "as_of": as_of, "photometry_version": sx.PHOTOMETRY_VERSION, "detectors": list(sx.COLOUR_DETECTORS),
                           "config_hash": short_hash, "git_commit": current_git_commit(root)}
                rec["manifest_id"] = sx.store_colour_inputs(cache, sid, payload)
                rule = spherex_colour_verdict(t_rec, comparison, cfg)
                rec.update(payload)
                rec.update({"status": "no_coverage" if not t_rec.get("n_images") else "measured", "verdict": rule["verdict"],
                            "reason": rule["reason"], "rule": rule})
        except Exception as exc:  # one star's archive failure must not stop the others; it is recorded with its cause
            rec["status"] = f"failed: {type(exc).__name__}: {exc}"[:400]
        records.append(rec)
        print(f"  SPHEREx verdict {n}/{len(ids)}: {sid}: {rec['status'].split(':')[0]}, verdict {rec.get('verdict')}", file=sys.stderr, flush=True)
    out_dir = root / cfg.storage.runs_dir / run_name
    out_dir.mkdir(parents=True, exist_ok=True)
    payload_out = {"run_id": run_name, "as_of": as_of, "photometry_version": sx.PHOTOMETRY_VERSION, "config_hash": short_hash,
                   "git_commit": current_git_commit(root), "stars": records}
    out = out_dir / "spherex_verdict.json"
    out.write_text(json.dumps(payload_out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    words = ", ".join(f"{r['source_id']}: {r.get('verdict') or r.get('reason') or r['status']}" for r in records)
    return CLIResponse(ok=True, summary=f"SPHEREx colour rule for {len(ids)} star(s): {words}", artifacts=[str(out)], data=payload_out)


def _cmd_status(args: argparse.Namespace, cfg: PipelineConfig, short_hash: str, root: Path) -> CLIResponse:
    # A channel counts as certified only when its certificate passed for the current config hash and the code it
    # validated is the code that would run a search now (the same condition `search` enforces).
    from oddstar.validate import all_certificate_status

    certified: dict[str, list[str]] = {}
    certificate_rows: dict[str, list[dict[str, Any]]] = {}
    for rel in ("dr3", "dr4"):
        rel_cfg = cfg if cfg.gaia.release == rel else cfg.model_copy(update={"gaia": cfg.gaia.model_copy(update={"release": rel})})
        certificate_rows[rel] = all_certificate_status(rel_cfg, root)
        certified[rel] = [r["channel"] for r in certificate_rows[rel] if r["usable"]]

    rows = read_ledger(root, cfg.storage.registry_dir)
    records = load_review_records(root)
    open_by_source = open_findings_by_source(rows, records)
    unresolved = sorted(f for fs in open_by_source.values() for f in fs)
    ev = verify_evidence(root, cfg.storage)
    cache: dict[str, Any] = {}
    review = {sid: review_status(root, sid, registry_dir=cfg.storage.registry_dir, _cache=cache).model_dump() for sid in reviewed_source_ids(root)}
    n_registrations = sum(1 for r in rows if r.get("event_type") in REGISTRATION_EVENTS)

    runs_dir = root / cfg.storage.runs_dir
    run_ids = sorted(p.name for p in runs_dir.iterdir() if p.is_dir() and not p.name.startswith("_")) if runs_dir.exists() else []
    uncertified_dr3 = [c for c in ALL_CHANNEL_NAMES if c not in certified.get("dr3", [])]
    next_cmd = f"oddstar validate {uncertified_dr3[0]} --json" if uncertified_dr3 else "oddstar inspect 4582145718700356096 --json"
    warnings = [f"evidence: {p}" for p in ev["problems"]]
    for rel, cert_rows in certificate_rows.items():
        for r in cert_rows:
            if not r["usable"]:
                warnings.append(f"certificate {r['channel']} {rel}: {'; '.join(r['reasons'])}")
    if ev["legacy_without_evidence"]:
        warnings.append(
            f"{ev['legacy_without_evidence']} of {n_registrations} ledger registrations are legacy_without_evidence: "
            "written before the committed evidence store, so the evidence behind their pass rows is not in the repository"
        )
    for sid, st in review.items():
        warnings.append(f"review {sid}: {st['summary']} ({st['independence_basis']})")
    return CLIResponse(
        ok=True,
        summary=(
            f"Certified channels DR3={len(certified.get('dr3', []))}/8, DR4={len(certified.get('dr4', []))}/8; "
            f"active runs={len(run_ids)}; open findings={len(unresolved)} across {len(open_by_source)} sources; "
            f"registrations with committed evidence={ev['ledger_entries_with_evidence']}, "
            f"legacy_without_evidence={ev['legacy_without_evidence']}"
        ),
        warnings=warnings,
        data={
            "certified": certified,
            "certificates": certificate_rows,
            "uncertified_dr3": uncertified_dr3,
            "active_runs": run_ids,
            "open_findings": unresolved,
            "open_findings_by_source": open_by_source,
            "review_status": review,
            "evidence": {k: v for k, v in ev.items() if k != "legacy_entry_ids"},
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
    p_search.add_argument(
        "--neowise-max-stars",
        type=int,
        default=None,
        help="ir_timedomain --live: fetch NEOWISE for at most N sample stars (default: all); the rest are reported as not searched",
    )
    p_search.add_argument(
        "--defer-enrichment",
        action="store_true",
        help="--live: do not fetch vetting data for flagged stars now (minutes per star); `vet --live` fetches it later",
    )
    p_search.add_argument(
        "--plate-max-fields",
        type=int,
        default=None,
        help="vanishing --live: read plate sources around at most N sample positions (default: all)",
    )
    p_search.add_argument(
        "--radio-data-dir",
        default=None,
        help="radio_narrowband --live: folder holding downloaded Breakthrough Listen scans (default: data/breakthrough_listen_scans)",
    )
    p_search.add_argument(
        "--radio-max-stars",
        type=int,
        default=None,
        help="radio_narrowband --live: look up radio scans for at most N sample stars (default: all); the rest are reported as not searched",
    )
    p_search.add_argument(
        "--spectra-max-stars",
        type=int,
        default=None,
        help="spectral_lines --live: look up LAMOST spectra for at most N sample stars (default: all); the rest are reported as not searched",
    )
    p_search.add_argument(
        "--tess-max-stars",
        type=int,
        default=None,
        help="optical_dips --live: fetch TESS light curves for at most N sample stars (default: all); the rest are reported as not searched",
    )

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
    p_ins.add_argument("--register", action="store_true", help="Register the dossier tier in the ledger with committed evidence")

    p_ev = sub.add_parser("evidence", parents=[common], help="Committed evidence store (registry/evidence)")
    ev_sub = p_ev.add_subparsers(dest="evidence_cmd", required=True)
    p_ee = ev_sub.add_parser("export", parents=[common], help="Export compact evidence for a candidate from run/dossier artefacts")
    p_ee.add_argument("source_id")
    p_ee.add_argument("--run-id", default=None)
    p_evv = ev_sub.add_parser("verify", parents=[common], help="Re-hash evidence and check ledger references")
    p_evv.add_argument("source_id", nargs="?", default=None)

    p_rep = sub.add_parser("report", parents=[common], help="Generated documents (docs/candidates.md, reviewed_candidate.md)")
    rep_sub = p_rep.add_subparsers(dest="report_cmd", required=True)
    rep_sub.add_parser("build", parents=[common], help="Rebuild the candidate documents from committed outputs")
    rep_sub.add_parser("check", parents=[common], help="Fail if a generated document differs from a fresh build")

    p_runs = sub.add_parser("runs", parents=[common], help="Committed run records (registry/runs, registry/samples)")
    runs_sub = p_runs.add_subparsers(dest="runs_cmd", required=True)
    p_rx = runs_sub.add_parser("export", parents=[common], help="Write the compact committed record of search or fusion runs")
    p_rx.add_argument("run_ids", nargs="+")
    p_rx.add_argument("--replace", action="store_true", help="Overwrite an existing record whose content differs")
    runs_sub.add_parser("verify", parents=[common], help="Check every committed run record against its sample files")

    p_rev = sub.add_parser("review", parents=[common], help="Stage 7 isolated adversarial review")
    rev_sub = p_rev.add_subparsers(dest="review_cmd", required=True)
    p_rp = rev_sub.add_parser("packet", parents=[common], help="Build blinded reviewer packet")
    p_rp.add_argument("ids")
    p_rp.add_argument("--round", type=int, default=1)
    p_rp.add_argument("--claim", default=None)
    p_ri = rev_sub.add_parser("ingest", parents=[common], help="Validate and record reviewer JSON")
    p_ri.add_argument("file", nargs="?", default=None)
    p_ri.add_argument("--all", action="store_true", help="Ingest every stored record under review/records/ (idempotent)")

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
    p_la.add_argument("--no-close", action="store_true", help="Link the finding without closing it")
    p_la.add_argument("--related-finding", action="append", default=None, help="Other finding ID addressed (not closed)")
    p_la.add_argument("--finding-source", default=None, help="Where the finding is recorded when not in a review record")
    p_la.add_argument("--also-test", action="append", default=None, help="Additional linked test")
    p_lr = les_sub.add_parser("register", parents=[common], help="Register an existing lesson file in the ledger")
    p_lr.add_argument("file")
    p_lm = les_sub.add_parser("amend", parents=[common], help="Correct lesson metadata and record the change")
    p_lm.add_argument("lesson_id")
    p_lm.add_argument("--set", action="append", default=None, metavar="FIELD=VALUE", help="VALUE is JSON or a plain string")
    p_lm.add_argument("--reason", required=True)
    p_lm.add_argument("--md-note", default=None)
    les_sub.add_parser("record-sha", parents=[common], help="Set commit_sha to the commit that introduced each lesson")

    p_radio = sub.add_parser("radio", parents=[common], help="Radio recordings, searched on request")
    radio_sub = p_radio.add_subparsers(dest="radio_cmd", required=True)
    p_rf = radio_sub.add_parser("followup", parents=[common],
                                help="Fetch and search the Breakthrough Listen recordings of chosen stars, one file at a time")
    p_rf.add_argument("--from-run", action="append", required=True, metavar="RUN_ID",
                      help="search run whose sample holds the stars (repeatable); its flagged stars and vetting give the tiers")
    p_rf.add_argument("--min-tier", choices=["A", "B", "C"], default="C", help="follow up stars at this tier or better (default C)")
    p_rf.add_argument("--ids", type=str, default=None, help="comma-separated source ids to follow up instead of the tier selection")
    p_rf.add_argument("--all-covered", action="store_true", help="follow up every star of the runs' samples that lies in a beam")
    p_rf.add_argument("--scan-sets", choices=["smallest", "all"], default="smallest",
                      help="per star: the complete scan set with the fewest bytes (default) or every complete scan set")
    p_rf.add_argument("--keep-windows", action="store_true", help="keep the windows cut from the scans (for re-scoring without a download)")
    p_rf.add_argument("--min-free-gb", type=float, default=15.0, help="do not fetch a scan set with less free disk than this (default 15)")
    p_rf.add_argument("--radio-data-dir", type=str, default=None, help="where scan files pass through (default: data/breakthrough_listen_scans)")
    p_rf.add_argument("--run-id", type=str, default=None, help="name of the output folder under runs/ (default radio-followup-<first run>)")

    p_rb = radio_sub.add_parser("beams", parents=[common],
                                help="Test flagged radio windows against the other beams recorded at the same instant (range requests, no download)")
    p_rb.add_argument("--from-run", action="append", required=True, metavar="RUN_ID",
                      help="search run whose sample holds the stars (repeatable); its flagged stars and vetting give the tiers")
    p_rb.add_argument("--min-tier", choices=["A", "B", "C"], default="C", help="test stars at this tier or better (default C)")
    p_rb.add_argument("--ids", type=str, default=None, help="comma-separated source ids to test instead of the tier selection")
    p_rb.add_argument("--redo", action="store_true", help="read the beams again even when a beam test is stored")
    p_rb.add_argument("--run-id", type=str, default=None, help="name of the output folder under runs/ (default radio-beams-<first run>)")

    p_sx = sub.add_parser("spherex", parents=[common], help="SPHEREx 0.75 to 5 micron spectral images")
    sx_sub = p_sx.add_subparsers(dest="spherex_cmd", required=True)
    p_sxp = sx_sub.add_parser("photometry", parents=[common],
                              help="Measure chosen stars on every SPHEREx image that holds them (one small cutout per exposure)")
    p_sxp.add_argument("--from-run", action="append", required=True, metavar="RUN_ID",
                       help="search run whose sample holds the stars (repeatable); its flagged stars and vetting give the tiers")
    p_sxp.add_argument("--min-tier", choices=["A", "B", "C"], default="C", help="measure stars at this tier or better (default C)")
    p_sxp.add_argument("--ids", type=str, default=None, help="comma-separated source ids to measure instead of the tier selection")
    p_sxp.add_argument("--detectors", type=str, default=None,
                       help="comma-separated detectors 1 to 6 (0.75-1.1, 1.1-1.6, 1.6-2.4, 2.4-3.8, 3.8-4.4, 4.4-5.0 micron); default all")
    p_sxp.add_argument("--as-of", type=str, default=None,
                       help="label of the image list (default: today's date, UTC); the same label reads the stored list again")
    p_sxp.add_argument("--workers", type=int, default=2, help="cutouts fetched at a time (default 2; the archive answers 503 when pushed)")
    p_sxp.add_argument("--no-neighbours", action="store_true", help="do not fit Gaia neighbours in the stamp with the star")
    p_sxp.add_argument("--run-id", type=str, default=None, help="name of the output folder under runs/ (default spherex-<first run>)")

    p_sxv = sx_sub.add_parser("verdict", parents=[common],
                              help="Measure stars and their comparison stars and apply the colour rule (feeds the spherex_colours check)")
    p_sxv.add_argument("--ids", type=str, required=True, help="comma-separated Gaia source ids")
    p_sxv.add_argument("--as-of", type=str, default=None,
                       help="label of the image list (default: today's date, UTC); the same label reads the stored list again")
    p_sxv.add_argument("--workers", type=int, default=2, help="cutouts fetched at a time (default 2; the archive answers 503 when pushed)")
    p_sxv.add_argument("--passes", type=int, default=2,
                       help="times a star is gone through while images are missing (default 2; a later pass fetches only what was refused)")
    p_sxv.add_argument("--run-id", type=str, default=None, help="name of the output folder under runs/ (default spherex-verdict-<first id>)")

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
        elif cmd == "lessons" and args.lessons_cmd == "register":
            resp = _cmd_lessons_register(args, cfg, short_hash, base)
        elif cmd == "lessons" and args.lessons_cmd == "amend":
            resp = _cmd_lessons_amend(args, cfg, short_hash, base)
        elif cmd == "lessons" and args.lessons_cmd == "record-sha":
            resp = _cmd_lessons_record_sha(args, cfg, short_hash, base)
        elif cmd == "evidence" and args.evidence_cmd == "export":
            resp = _cmd_evidence_export(args, cfg, short_hash, base)
        elif cmd == "evidence" and args.evidence_cmd == "verify":
            resp = _cmd_evidence_verify(args, cfg, short_hash, base)
        elif cmd == "report":
            resp = _cmd_report(args, cfg, short_hash, base)
        elif cmd == "runs" and args.runs_cmd == "export":
            resp = _cmd_runs_export(args, cfg, short_hash, base)
        elif cmd == "runs" and args.runs_cmd == "verify":
            resp = _cmd_runs_verify(args, cfg, short_hash, base)
        elif cmd == "radio" and args.radio_cmd == "followup":
            resp = _cmd_radio_followup(args, cfg, short_hash, base)
        elif cmd == "radio" and args.radio_cmd == "beams":
            resp = _cmd_radio_beams(args, cfg, short_hash, base)
        elif cmd == "spherex" and args.spherex_cmd == "photometry":
            resp = _cmd_spherex_photometry(args, cfg, short_hash, base)
        elif cmd == "spherex" and args.spherex_cmd == "verdict":
            resp = _cmd_spherex_verdict(args, cfg, short_hash, base)
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
