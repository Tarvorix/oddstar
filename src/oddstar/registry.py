"""Append-only registry ledger access shared by the CLI, review loop, evidence store and CI checks."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from oddstar.schemas import RegistryEntry

LEDGER_NAME = "ledger.jsonl"
REGISTRATION_EVENTS = {"candidate_registered", "explained_anomaly"}


def ledger_path(root: Path, registry_dir: str = "registry") -> Path:
    return root / registry_dir / LEDGER_NAME


def read_ledger(root: Path, registry_dir: str = "registry") -> list[dict[str, Any]]:
    """All ledger rows in file order (blank lines skipped)."""
    path = ledger_path(root, registry_dir)
    if not path.exists():
        return []
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


def append_entry(root: Path, entry: RegistryEntry, registry_dir: str = "registry") -> Path:
    """Append one entry (never rewrites earlier lines).

    Rows are compact `model_dump_json` lines like the existing ledger lines; the evidence reference fields are
    written only when set, so non-registration events keep their earlier shape.
    """
    path = ledger_path(root, registry_dir)
    path.parent.mkdir(parents=True, exist_ok=True)
    unset = {k for k in ("evidence_dir", "evidence_index_sha256") if getattr(entry, k) is None}
    with path.open("a", encoding="utf-8") as fh:
        fh.write(entry.model_dump_json(exclude=unset) + "\n")
    return path


def entry_ids(rows: list[dict[str, Any]]) -> set[str]:
    return {str(r.get("entry_id")) for r in rows}


def clean_source_id(source_id: str) -> str:
    return source_id.replace("Gaia DR3 ", "").replace("Gaia DR4 ", "").strip()
