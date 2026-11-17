"""Adversarial review packet builder, reviewer record ingest, and R7 lessons loop."""

from __future__ import annotations

import hashlib
import json
import subprocess
import uuid
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig, current_git_commit
from oddstar.errors import ArtifactNotFoundError, PolicyViolationError
from oddstar.registry import append_entry, clean_source_id, entry_ids, read_ledger
from oddstar.review.status import independence, load_review_records, model_family
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
    return append_entry(root, entry, cfg.storage.registry_dir)


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
    """Validate and record a reviewer JSON record, checking reliability and opening findings.

    Idempotent: a record is stored once under `review/records/<review_id>.json` (a record already stored there is
    read in place and never rewritten; different content under the same review_id is refused), and the
    `review_ingested` / `finding_opened` ledger events are appended only when not already present.
    Every event carries the normalised `model_family`; a round-3 event also records `independent_family` (true only
    when its recognised family differs from every recognised family of rounds 1-2 for the same source).
    """
    if not record_path.exists():
        raise ArtifactNotFoundError(f"Review record file not found: {record_path}")
    rec = ReviewRecord.model_validate_json(record_path.read_text(encoding="utf-8"))

    rec_dir = root / "review" / "records"
    rec_dir.mkdir(parents=True, exist_ok=True)
    stored_path = rec_dir / f"{rec.review_id}.json"
    if record_path.resolve() != stored_path.resolve():
        if stored_path.exists():
            existing = ReviewRecord.model_validate_json(stored_path.read_text(encoding="utf-8"))
            if existing != rec:
                raise PolicyViolationError(
                    f"Review record {rec.review_id} is already stored with different content",
                    hint="Review records are immutable; give a new review a new review_id",
                )
        else:
            stored_path.write_text(rec.model_dump_json(indent=2) + "\n", encoding="utf-8")

    discrepancy = rec.verdict in {"not_reproduced", "unresolved"} or any(
        f.severity == "material" and f.status == "open" for f in rec.findings
    )
    now = datetime.now(timezone.utc).isoformat()
    clean_id = clean_source_id(rec.source_id)
    family = model_family(rec.reviewer_model)
    rows = read_ledger(root, cfg.storage.registry_dir)
    known = entry_ids(rows)
    opened = {(r.get("details") or {}).get("finding_id") for r in rows if r.get("event_type") == "finding_opened"}
    appended: list[str] = []

    entry_id = f"REG-{rec.review_id}"
    provenance: dict[str, Any] = {
        "git_commit": current_git_commit(root),
        "round": rec.round,
        "reviewer_model": rec.reviewer_model,
        "model_family": family,
    }
    if rec.round == 3:
        same_source = [r for _, r in load_review_records(root) if clean_source_id(r.source_id) == clean_id and r.round in (1, 2)]
        indep, basis, _ = independence([*same_source, rec])
        provenance["independent_family"] = bool(
            indep and family not in {model_family(r.reviewer_model) for r in same_source}
        )
        provenance["earlier_round_families"] = sorted({model_family(r.reviewer_model) for r in same_source})
        provenance["independence_basis"] = basis
    if entry_id not in known:
        _append_registry(
            root,
            cfg,
            RegistryEntry(
                entry_id=entry_id,
                timestamp_utc=now,
                source_id=clean_id,
                event_type="review_ingested",
                tier="C" if discrepancy else None,
                status=rec.verdict,
                provenance=provenance,
                details={"review_path": str(stored_path.relative_to(root)), "discrepancy": discrepancy},
            ),
        )
        appended.append(entry_id)

    for f in rec.findings:
        fid_entry = f"REG-{f.finding_id}"
        if f.status == "open" and fid_entry not in known and f.finding_id not in opened:
            _append_registry(
                root,
                cfg,
                RegistryEntry(
                    entry_id=fid_entry,
                    timestamp_utc=now,
                    source_id=clean_id,
                    event_type="finding_opened",
                    status="open",
                    provenance={"review_id": rec.review_id, "round": rec.round, "model_family": family},
                    details=f.model_dump(),
                ),
            )
            appended.append(fid_entry)

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
            model_family=family,
        )
        new_reviews = [existing_rounds[k] for k in sorted(existing_rounds)]
        new_tier = "C" if (discrepancy and dossier.tier in {"A", "B"}) else dossier.tier
        from oddstar.review.status import review_status

        updated = dossier.model_copy(
            update={"reviews": new_reviews, "tier": new_tier, "review_status": review_status(root, clean_id, registry_dir=cfg.storage.registry_dir)}
        )
        dossier_json.write_text(updated.model_dump_json(indent=2) + "\n", encoding="utf-8")

    return {
        "review_id": rec.review_id,
        "stored_path": str(stored_path),
        "discrepancy": discrepancy,
        "open_findings": [f.finding_id for f in rec.findings if f.status == "open"],
        "model_family": family,
        "independent_family": provenance.get("independent_family"),
        "appended_entries": appended,
        "already_ingested": not appended,
    }


def ingest_all_stored_records(cfg: PipelineConfig, root: Path) -> list[dict[str, Any]]:
    """Ingest every stored record under review/records/ (rounds in order), skipping events already in the ledger."""
    recs = sorted(load_review_records(root), key=lambda pr: (pr[1].round, pr[1].timestamp_utc, pr[1].review_id))
    return [ingest_review_record(p, cfg, root) for p, _ in recs]


# --------------------------------------------------------------------------------------------- lessons (R7)


def _lesson_path(root: Path, lesson_id: str) -> Path:
    return root / "lessons" / f"{lesson_id}.json"


def _check_tests_exist(root: Path, tests: list[str]) -> None:
    for t in tests:
        test_file = root / t.split("::")[0]
        if not test_file.exists():
            raise ArtifactNotFoundError(
                f"Linked test file does not exist: {test_file}",
                hint="R7 requires every lesson to link an existing automated test or check.",
            )


def _lesson_markdown(rec: LessonRecord, body: str | None = None) -> str:
    head = [
        f"# {rec.lesson_id}: {rec.title}",
        "",
        f"- **Finding ID**: `{rec.finding_id}`" + ("" if rec.closes_finding else " (linked, not closed by this lesson)"),
    ]
    if rec.related_finding_ids:
        head.append("- **Related findings (addressed, not closed)**: " + ", ".join(f"`{f}`" for f in rec.related_finding_ids))
    if rec.finding_source:
        head.append(f"- **Finding source**: {rec.finding_source}")
    head += [
        f"- **Target**: `{rec.source_id}`",
        f"- **Automated Check / Vetting Rule**: `{rec.automated_check}`",
        f"- **Linked Test**: `{rec.linked_test}`",
    ]
    if rec.additional_linked_tests:
        head.append("- **Additional Tests**: " + ", ".join(f"`{t}`" for t in rec.additional_linked_tests))
    head.append(f"- **Added to Impostors**: `{rec.added_to_impostors}`")
    sha = rec.commit_sha or "not recorded yet (filled by `oddstar lessons record-sha` after the commit that adds this lesson)"
    head.append(f"- **Commit (introduced the lesson files)**: `{sha}`" + (f" ({rec.commit_sha_note})" if rec.commit_sha_note else ""))
    if body is None:
        body = "\n".join(
            [
                "## What Happened",
                rec.what_happened,
                "",
                "## Generic Lesson",
                rec.generic_lesson,
                "",
                "## Automation & Regression Test",
                f"Enforced via `{rec.automated_check}` and verified in `{rec.linked_test}`.",
                "",
            ]
        )
    return "\n".join(head) + "\n\n" + body.lstrip("\n")


def _md_body(md_path: Path) -> str | None:
    """Everything after the header bullet list of an existing lesson Markdown file."""
    if not md_path.exists():
        return None
    lines = md_path.read_text(encoding="utf-8").splitlines()
    i = 0
    if lines and lines[0].startswith("# "):
        i = 1
    while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith("- **")):
        i += 1
    return "\n".join(lines[i:]) + "\n"


def _write_lesson_files(root: Path, rec: LessonRecord, body: str | None = None) -> tuple[Path, Path]:
    lessons_dir = root / "lessons"
    lessons_dir.mkdir(parents=True, exist_ok=True)
    json_path = lessons_dir / f"{rec.lesson_id}.json"
    md_path = lessons_dir / f"{rec.lesson_id}.md"
    json_path.write_text(rec.model_dump_json(indent=2) + "\n", encoding="utf-8")
    md_path.write_text(_lesson_markdown(rec, body), encoding="utf-8")
    return md_path, json_path


def register_lesson(rec: LessonRecord, cfg: PipelineConfig, root: Path) -> list[str]:
    """Append `lesson_added` (and `finding_closed` when the lesson closes its finding); idempotent."""
    _check_tests_exist(root, [rec.linked_test, *rec.additional_linked_tests])
    rows = read_ledger(root, cfg.storage.registry_dir)
    known = entry_ids(rows)
    now = datetime.now(timezone.utc).isoformat()
    appended: list[str] = []
    if f"REG-{rec.lesson_id}" not in known:
        _append_registry(
            root,
            cfg,
            RegistryEntry(
                entry_id=f"REG-{rec.lesson_id}",
                timestamp_utc=now,
                source_id=rec.source_id,
                event_type="lesson_added",
                status="recorded",
                provenance={"git_commit": rec.commit_sha or current_git_commit(root), "commit_sha_pending": not rec.commit_sha},
                details=rec.model_dump(),
            ),
        )
        appended.append(f"REG-{rec.lesson_id}")
    close_id = f"REG-CLOSE-{rec.finding_id}"
    if rec.closes_finding and close_id not in known:
        _append_registry(
            root,
            cfg,
            RegistryEntry(
                entry_id=close_id,
                timestamp_utc=now,
                source_id=rec.source_id,
                event_type="finding_closed",
                status="closed",
                provenance={"lesson_id": rec.lesson_id},
                details={"finding_id": rec.finding_id, "linked_test": rec.linked_test, "automated_check": rec.automated_check},
            ),
        )
        appended.append(close_id)
    return appended


def register_lesson_file(json_path: Path, cfg: PipelineConfig, root: Path) -> tuple[LessonRecord, list[str]]:
    """Register an existing lesson file (written earlier without ledger events) without rewriting it."""
    if not json_path.exists():
        raise ArtifactNotFoundError(f"Lesson file not found: {json_path}")
    rec = LessonRecord.model_validate_json(json_path.read_text(encoding="utf-8"))
    return rec, register_lesson(rec, cfg, root)


def amend_lesson(
    lesson_id: str,
    updates: dict[str, Any],
    *,
    reason: str,
    cfg: PipelineConfig,
    root: Path,
    md_note: str | None = None,
) -> tuple[LessonRecord, list[str]]:
    """Correct lesson metadata in place and record the change.

    The lesson JSON and the header of its Markdown are rewritten (the Markdown body is kept; `md_note` is appended as
    a section). When the lesson is already registered, a `lesson_updated` event records old and new values; when the
    amended lesson now closes a finding that has no `finding_closed` event, that event is appended too.
    """
    path = _lesson_path(root, lesson_id)
    if not path.exists():
        raise ArtifactNotFoundError(f"Lesson file not found: {path}")
    old = LessonRecord.model_validate_json(path.read_text(encoding="utf-8"))
    unknown = set(updates) - set(LessonRecord.model_fields)
    if unknown:
        raise PolicyViolationError(f"Unknown lesson fields: {sorted(unknown)}")
    new = LessonRecord.model_validate({**old.model_dump(), **updates})
    _check_tests_exist(root, [new.linked_test, *new.additional_linked_tests])
    changes = {k: [getattr(old, k), getattr(new, k)] for k in updates if getattr(old, k) != getattr(new, k)}
    md_path = root / "lessons" / f"{lesson_id}.md"
    body = _md_body(md_path)
    if md_note:
        body = (body or "").rstrip("\n") + "\n\n" + md_note.strip("\n") + "\n"
    if not changes and not md_note:
        return old, []
    _write_lesson_files(root, new, body)
    rows = read_ledger(root, cfg.storage.registry_dir)
    known = entry_ids(rows)
    appended: list[str] = []
    if f"REG-{lesson_id}" in known and changes:
        n = sum(1 for e in known if e.startswith(f"REG-UPDATE-{lesson_id}-")) + 1
        eid = f"REG-UPDATE-{lesson_id}-{n}"
        _append_registry(
            root,
            cfg,
            RegistryEntry(
                entry_id=eid,
                timestamp_utc=datetime.now(timezone.utc).isoformat(),
                source_id=new.source_id,
                event_type="lesson_updated",
                status="amended",
                provenance={"lesson_id": lesson_id, "git_commit": current_git_commit(root)},
                details={"lesson_id": lesson_id, "changes": changes, "reason": reason},
            ),
        )
        appended.append(eid)
        if new.closes_finding and f"REG-CLOSE-{new.finding_id}" not in known:
            appended += register_lesson(new, cfg, root)
    return new, appended


def lesson_introducing_commit(root: Path, json_path: Path) -> str | None:
    """Commit that added the lesson JSON (`git log --diff-filter=A`), or None when not committed / no history."""
    try:
        res = subprocess.run(
            ["git", "log", "--diff-filter=A", "--format=%H", "--", str(json_path.relative_to(root))],
            cwd=str(root), capture_output=True, text=True, check=False,
        )
    except OSError:
        return None
    shas = [x for x in res.stdout.split() if x]
    return shas[-1] if res.returncode == 0 and shas else None


def record_lesson_commit_shas(cfg: PipelineConfig, root: Path) -> list[dict[str, Any]]:
    """Set every lesson's `commit_sha` to the commit that introduced its files (amending registered lessons)."""
    results = []
    for jp in sorted((root / "lessons").glob("*.json")):
        rec = LessonRecord.model_validate_json(jp.read_text(encoding="utf-8"))
        intro = lesson_introducing_commit(root, jp)
        if intro is None:
            results.append({"lesson_id": rec.lesson_id, "commit_sha": rec.commit_sha, "action": "not committed yet"})
            continue
        if rec.commit_sha == intro:
            results.append({"lesson_id": rec.lesson_id, "commit_sha": intro, "action": "unchanged"})
            continue
        _, appended = amend_lesson(
            rec.lesson_id,
            {"commit_sha": intro, "commit_sha_note": ""},
            reason="commit_sha is defined as the commit that introduced the lesson's files (git log --diff-filter=A)",
            cfg=cfg,
            root=root,
        )
        results.append({"lesson_id": rec.lesson_id, "commit_sha": intro, "previous": rec.commit_sha, "action": "amended", "ledger": appended})
    return results


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
    closes_finding: bool = True,
    related_finding_ids: list[str] | None = None,
    finding_source: str = "",
    additional_linked_tests: list[str] | None = None,
) -> tuple[LessonRecord, Path, Path]:
    """Create lesson Markdown + JSON record, register it (closing its finding), and validate linked tests.

    `commit_sha` stays empty: it is the commit that will introduce these files, which does not exist yet. Run
    `oddstar lessons record-sha` after that commit.
    """
    _check_tests_exist(root, [linked_test, *(additional_linked_tests or [])])
    if _lesson_path(root, lesson_id).exists() and f"REG-{lesson_id}" in entry_ids(read_ledger(root, cfg.storage.registry_dir)):
        raise PolicyViolationError(f"Lesson {lesson_id} already exists", hint="Use `oddstar lessons amend` to correct it")
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
        commit_sha="",
        closes_finding=closes_finding,
        related_finding_ids=list(related_finding_ids or []),
        finding_source=finding_source,
        additional_linked_tests=list(additional_linked_tests or []),
    )
    md_path, json_path = _write_lesson_files(root, rec)
    register_lesson(rec, cfg, root)
    return rec, md_path, json_path
