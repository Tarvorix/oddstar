"""Shared vetting helpers: evidence-backed tri-state results and the SIMBAD lookup outcome."""

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
    inputs: dict[str, Any] | None = None,
) -> VettingResult:
    """Write check evidence to disk and return a validated VettingResult.

    Enforces AGENTS.md rule: never return 'pass' or 'fail' without an existing evidence file. `inputs` optionally
    records the raw-cache manifest IDs, SHA-256 values and query URLs the check consumed; the committed evidence
    store (`oddstar.evidence`) copies them, together with any it finds in the metrics, when a candidate is registered.
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
        if inputs:
            payload["inputs"] = inputs
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


def simbad_lookup(star: dict[str, Any]) -> dict[str, Any]:
    """Outcome of the SIMBAD classification lookup, keeping "failed", "no entry" and "type X" distinct.

    - A recorded `vetting_archive.query_status["simbad"]` decides for every star: "ok" means the lookup succeeded
      and `simbad_type` is its answer (None = SIMBAD has no entry, a real negative answer); anything else means
      the lookup did not succeed and `simbad_type` is ignored.
    - A live star (`live_enriched` key) without a recorded status was never looked up.
    - A fixture or hand-built dict without a recorded status counts an explicit, non-None `simbad_type`.

    Returns {"evaluated", "simbad_type", "entry", "query_status"}; `entry` is None when it cannot be told (a
    record enriched before `vetting_archive.simbad_entry` was written).
    """
    va = star.get("vetting_archive") or {}
    status = (va.get("query_status") or {}).get("simbad")
    stype = star.get("simbad_type")
    if status is not None:
        if status != "ok":
            return {"evaluated": False, "simbad_type": None, "entry": None, "query_status": str(status)}
        entry = va.get("simbad_entry")
        if entry is None and stype is not None:
            entry = True
        return {"evaluated": True, "simbad_type": None if stype is None else str(stype), "entry": entry, "query_status": "ok"}
    if "live_enriched" in star:
        return {"evaluated": False, "simbad_type": None, "entry": None, "query_status": "not run"}
    if stype is None:
        return {"evaluated": False, "simbad_type": None, "entry": None, "query_status": "absent (no explicit simbad_type)"}
    return {"evaluated": True, "simbad_type": str(stype), "entry": True, "query_status": "explicit"}


def simbad_note(sim: dict[str, Any]) -> str:
    """Short text for a SIMBAD lookup outcome from `simbad_lookup`."""
    if not sim["evaluated"]:
        return f"SIMBAD lookup not available (query_status simbad={sim['query_status']})"
    if sim["simbad_type"] is None:
        return "SIMBAD: no entry" if sim["entry"] is not True else "SIMBAD: entry without object type"
    return f"SIMBAD type={sim['simbad_type']}"
