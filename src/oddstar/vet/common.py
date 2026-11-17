"""Shared vetting helper enforcing evidence-backed tri-state results."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from oddstar.schemas import VettingResult


def make_vet_result(
    out_dir: Path,
    *,
    check: str,
    status: Literal["pass", "fail", "not_evaluated"],
    note: str,
    metrics: dict[str, Any] | None = None,
    evidence_override: Path | None = None,
) -> VettingResult:
    """Write check evidence to disk and return a validated VettingResult.

    Enforces AGENTS.md rule: never return 'pass' or 'fail' without an existing evidence file.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    if status == "not_evaluated" and evidence_override is None:
        return VettingResult(
            check=check,
            status="not_evaluated",
            evidence=None,
            note=note,
            metrics=metrics or {},
        )

    ev_path = evidence_override or (out_dir / f"vet_{check}.json")
    if evidence_override is None:
        payload = {
            "check": check,
            "status": status,
            "note": note,
            "metrics": metrics or {},
        }
        ev_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    elif not ev_path.exists():
        raise FileNotFoundError(f"Expected evidence file does not exist: {ev_path}")

    return VettingResult(
        check=check,
        status=status,
        evidence=str(ev_path),
        note=note,
        metrics=metrics or {},
    )
