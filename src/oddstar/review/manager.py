"""Adversarial review packet builder, reviewer record ingest, and R7 lessons loop."""

from __future__ import annotations

import hashlib
import json
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig, current_git_commit
from oddstar.errors import ArtifactNotFoundError, PolicyViolationError
from oddstar.io.fixtures_data import generate_field_sample, get_named_object
from oddstar.schemas import (
    Dossier,
    DossierReviewSummary,
    LessonRecord,
    RegistryEntry,
    ReviewPacket,
    ReviewPacketTarget,
    ReviewRecord,
)


def _append_registry(root: Path, cfg: PipelineConfig, entry: RegistryEntry) -> Path:
    reg_dir = root / cfg.storage.registry_dir
    reg_dir.mkdir(parents=True, exist_ok=True)
    ledger = reg_dir / "ledger.jsonl"
    with ledger.open("a", encoding="utf-8") as fh:
        fh.write(entry.model_dump_json() + "\n")
    return ledger


def build_review_packet(
    candidate_ids: list[str],
    cfg: PipelineConfig,
    root: Path,
    *,
    round_num: int = 1,
    claim: str | None = None,
) -> tuple[ReviewPacket, Path, Path]:
    """Build a blinded reviewer packet mixed with controls and known objects (Section 9)."""
    prompt_map = {
        1: "review/prompts/round1_blind.md",
        2: "review/prompts/round2_adversarial.md",
        3: "review/prompts/round3_different_eyes.md",
    }
    if round_num not in prompt_map:
        raise PolicyViolationError(f"Unsupported review round: {round_num}", hint="Use round 1, 2, or 3")

    prompt_rel = prompt_map[round_num]
    prompt_path = root / prompt_rel
    prompt_text = prompt_path.read_text(encoding="utf-8") if prompt_path.exists() else ""

    packet_seed = int.from_bytes(
        hashlib.sha256(f"{','.join(candidate_ids)}:R{round_num}".encode()).digest()[:4], "big"
    )
    rng = np.random.default_rng(packet_seed)
    packet_id = f"PKT-R{round_num}-{hashlib.sha256(f'{candidate_ids}:{round_num}'.encode()).hexdigest()[:10]}"

    items_with_roles: list[tuple[dict[str, Any], str, str | None]] = []
    for cid in candidate_ids:
        clean = cid.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
        star = get_named_object(clean, gaia_release=cfg.gaia.release)
        if star is None:
            from oddstar.io.fixtures_data import _make_base_star

            star = _make_base_star(clean, gaia_release=cfg.gaia.release)
        c_text = claim
        if round_num in (2, 3) and not c_text:
            # Bare numerical measurement claim without pipeline interpretation
            c_text = (
                f"Reported mid-infrared 22 um flux for Gaia DR3 {clean} exceeds photospheric "
                f"prediction at W4 (see archive photometry at RA={star['ra_deg']:.4f}, Dec={star['dec_deg']:.4f})."
            )
        items_with_roles.append((star, "candidate", c_text))

    if round_num == 1:
        # Mix equal number of ordinary control stars + known positives + known impostors
        controls = [
            s
            for s in generate_field_sample(n_stars=50, seed=cfg.sample.seed + 777, gaia_release=cfg.gaia.release)
            if abs(s["gal_b_deg"]) > 10.0
        ][: max(1, len(candidate_ids))]
        for ctrl in controls:
            items_with_roles.append((ctrl, "control", None))
        for kp in ("115725517896200320", "6133409280036818688"):
            obj = get_named_object(kp, gaia_release=cfg.gaia.release)
            if obj is not None:
                items_with_roles.append((obj, "known_positive", None))
        for ki in ("4418120359123456128", "4129384756102938112"):
            obj = get_named_object(ki, gaia_release=cfg.gaia.release)
            if obj is not None:
                items_with_roles.append((obj, "known_impostor", None))

        perm = rng.permutation(len(items_with_roles))
        items_with_roles = [items_with_roles[i] for i in perm]

    targets: list[ReviewPacketTarget] = []
    answer_key: list[dict[str, str]] = []
    for idx, (star, role, c_text) in enumerate(items_with_roles):
        tid = f"{packet_id}-T{idx:03d}"
        targets.append(
            ReviewPacketTarget(
                packet_target_id=tid,
                source_id=str(star["source_id"]),
                ra_deg=float(star["ra_deg"]),
                dec_deg=float(star["dec_deg"]),
                claim=c_text if round_num in (2, 3) else None,
            )
        )
        answer_key.append({"packet_target_id": tid, "source_id": str(star["source_id"]), "role": role})

    packet = ReviewPacket(
        packet_id=packet_id,
        round=round_num,
        created_at_utc=datetime.now(timezone.utc).isoformat(),
        prompt_template_path=prompt_rel,
        prompt_text=prompt_text,
        targets=targets,
    )

    pkt_dir = root / "review" / "packets"
    pkt_dir.mkdir(parents=True, exist_ok=True)
    pkt_path = pkt_dir / f"{packet_id}.json"
    pkt_path.write_text(packet.model_dump_json(indent=2) + "\n", encoding="utf-8")

    # Separate answer key stored inside runs/review_keys (inaccessible to review_kit/)
    key_dir = root / cfg.storage.runs_dir / "review_keys"
    key_dir.mkdir(parents=True, exist_ok=True)
    key_path = key_dir / f"{packet_id}_key.json"
    key_path.write_text(
        json.dumps({"packet_id": packet_id, "round": round_num, "key": answer_key}, indent=2) + "\n",
        encoding="utf-8",
    )

    # Round 3 also produces a one-page human checklist over not_evaluated dossier rows
    if round_num == 3:
        checklist_lines = [f"# Round 3 Human Checklist — {packet_id}", ""]
        for cid in candidate_ids:
            clean = cid.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
            d_path = root / cfg.storage.dossiers_dir / clean / "dossier.json"
            checklist_lines.append(f"## Target `Gaia DR3 {clean}`")
            if d_path.exists():
                d_obj = Dossier.model_validate_json(d_path.read_text(encoding="utf-8"))
                holes = [v for v in d_obj.vetting if v.status == "not_evaluated"]
                if holes:
                    for h in holes:
                        checklist_lines.append(f"- [ ] Verify open hole `{h.check}`: {h.note}")
                else:
                    checklist_lines.append("- [x] All archival vetting checks evaluated.")
            else:
                checklist_lines.append("- [ ] Generate candidate dossier and inspect `not_evaluated` rows.")
            checklist_lines.append("")
        (pkt_dir / f"{packet_id}_human_checklist.md").write_text("\n".join(checklist_lines), encoding="utf-8")

    return packet, pkt_path, key_path


def ingest_review_record(
    record_path: Path,
    cfg: PipelineConfig,
    root: Path,
) -> dict[str, Any]:
    """Validate and record a reviewer JSON record, checking reliability and opening findings."""
    if not record_path.exists():
        raise ArtifactNotFoundError(f"Review record file not found: {record_path}")
    rec = ReviewRecord.model_validate_json(record_path.read_text(encoding="utf-8"))

    rec_dir = root / "review" / "records"
    rec_dir.mkdir(parents=True, exist_ok=True)
    stored_path = rec_dir / f"{rec.review_id}.json"
    stored_path.write_text(rec.model_dump_json(indent=2) + "\n", encoding="utf-8")

    discrepancy = rec.verdict in {"not_reproduced", "unresolved"} or any(
        f.severity == "material" and f.status == "open" for f in rec.findings
    )
    now = datetime.now(timezone.utc).isoformat()
    clean_id = rec.source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()

    # Record in append-only registry
    _append_registry(
        root,
        cfg,
        RegistryEntry(
            entry_id=f"REG-{rec.review_id}",
            timestamp_utc=now,
            source_id=clean_id,
            event_type="review_ingested",
            tier="C" if discrepancy else None,
            status=rec.verdict,
            provenance={"git_commit": current_git_commit(root), "round": rec.round, "reviewer_model": rec.reviewer_model},
            details={"review_path": str(stored_path.relative_to(root)), "discrepancy": discrepancy},
        ),
    )

    for f in rec.findings:
        if f.status == "open":
            _append_registry(
                root,
                cfg,
                RegistryEntry(
                    entry_id=f"REG-{f.finding_id}",
                    timestamp_utc=now,
                    source_id=clean_id,
                    event_type="finding_opened",
                    status="open",
                    provenance={"review_id": rec.review_id, "round": rec.round},
                    details=f.model_dump(),
                ),
            )

    # Update dossier reviews and enforce Tier C downgrade if discrepancy exists
    dossier_json = root / cfg.storage.dossiers_dir / clean_id / "dossier.json"
    if dossier_json.exists():
        dossier = Dossier.model_validate_json(dossier_json.read_text(encoding="utf-8"))
        existing_rounds = {r.round: r for r in dossier.reviews}
        existing_rounds[rec.round] = DossierReviewSummary(
            round=rec.round,
            record=str(stored_path.relative_to(root)),
            verdict=rec.verdict,
            reviewer_model=rec.reviewer_model,
        )
        new_reviews = [existing_rounds[k] for k in sorted(existing_rounds)]
        new_tier = "C" if (discrepancy and dossier.tier in {"A", "B"}) else dossier.tier
        updated = dossier.model_copy(update={"reviews": new_reviews, "tier": new_tier})
        dossier_json.write_text(updated.model_dump_json(indent=2) + "\n", encoding="utf-8")

    return {
        "review_id": rec.review_id,
        "stored_path": str(stored_path),
        "discrepancy": discrepancy,
        "open_findings": [f.finding_id for f in rec.findings if f.status == "open"],
    }


def add_lesson(
    *,
    lesson_id: str,
    finding_id: str,
    source_id: str,
    title: str,
    what_happened: str,
    generic_lesson: str,
    automated_check: str,
    linked_test: str,
    added_to_impostors: bool,
    skill_updated: str | None,
    cfg: PipelineConfig,
    root: Path,
) -> tuple[LessonRecord, Path, Path]:
    """Create lesson Markdown + JSON record, close finding in registry, and validate linked test."""
    test_file = root / linked_test.split("::")[0]
    if not test_file.exists():
        raise ArtifactNotFoundError(
            f"Linked test file does not exist: {test_file}",
            hint="R7 requires every lesson to link an existing automated test or check.",
        )

    lessons_dir = root / "lessons"
    lessons_dir.mkdir(parents=True, exist_ok=True)
    rec = LessonRecord(
        lesson_id=lesson_id,
        finding_id=finding_id,
        source_id=source_id,
        title=title,
        what_happened=what_happened,
        generic_lesson=generic_lesson,
        automated_check=automated_check,
        linked_test=linked_test,
        added_to_impostors=added_to_impostors,
        skill_updated=skill_updated,
        commit_sha=current_git_commit(root),
    )

    json_path = lessons_dir / f"{lesson_id}.json"
    md_path = lessons_dir / f"{lesson_id}.md"
    json_path.write_text(rec.model_dump_json(indent=2) + "\n", encoding="utf-8")
    md_path.write_text(
        "\n".join(
            [
                f"# {lesson_id}: {title}",
                "",
                f"- **Finding ID**: `{finding_id}`",
                f"- **Target**: `{source_id}`",
                f"- **Automated Check / Vetting Rule**: `{automated_check}`",
                f"- **Linked Test**: `{linked_test}`",
                f"- **Added to Impostors**: `{added_to_impostors}`",
                "",
                "## What Happened",
                what_happened,
                "",
                "## Generic Lesson",
                generic_lesson,
                "",
                "## Automation & Regression Test",
                f"Enforced via `{automated_check}` and verified in `{linked_test}`.",
                "",
            ]
        ),
        encoding="utf-8",
    )

    now = datetime.now(timezone.utc).isoformat()
    _append_registry(
        root,
        cfg,
        RegistryEntry(
            entry_id=f"REG-{lesson_id}",
            timestamp_utc=now,
            source_id=source_id,
            event_type="lesson_added",
            status="recorded",
            provenance={"git_commit": rec.commit_sha},
            details=rec.model_dump(),
        ),
    )
    _append_registry(
        root,
        cfg,
        RegistryEntry(
            entry_id=f"REG-CLOSE-{finding_id}",
            timestamp_utc=now,
            source_id=source_id,
            event_type="finding_closed",
            status="closed",
            provenance={"lesson_id": lesson_id},
            details={"finding_id": finding_id, "linked_test": linked_test, "automated_check": automated_check},
        ),
    )
    return rec, md_path, json_path
