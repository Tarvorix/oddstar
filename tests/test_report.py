"""The candidate documents are built from committed outputs, so their numbers cannot drift from the pipeline's."""

from __future__ import annotations

import json
from pathlib import Path

from oddstar.report import build_candidates_md, build_reviewed_md, stale_documents, write_documents
from oddstar.run_records import export_run_record

from test_run_records import _fake_run


def _live(root: Path, run_id: str, **kw: object) -> None:
    rd = _fake_run(root, run_id, **kw)  # type: ignore[arg-type]
    sf = next(rd.glob("search_*.json"))
    data = json.loads(sf.read_text())
    data["cuts"][0]["cut_name"] = "gaia_halo_live"
    data["cuts"][0]["description"] = "Live Gaia DR3 x AllWISE x 2MASS region='halo' (n=6, after_random_index=None)"
    sf.write_text(json.dumps(data))
    export_run_record(root, run_id)


def _register(root: Path, sid: str, tier: str, run_id: str) -> None:
    ev = root / "registry" / "evidence" / sid / run_id
    ev.mkdir(parents=True)
    vet = [("moving_object", "pass", "No known Solar System body within 15\" at any of 19 frame epochs"),
           ("binary", "not_evaluated", "Gaia RVS scatter 6.13 km/s is not constraining | pipes are escaped")]
    if tier == "Rejected":
        vet = [("centroid_offset", "fail", "W4 centroid 9.1\" from the star")]
    files = []
    for check, status, note in vet:
        name = f"vet_{check}.json"
        (ev / name).write_text(json.dumps({"check": check, "status": status, "note": note, "metrics": {}}))
        files.append({"path": name, "sha256": "0" * 64, "bytes": 1, "kind": "vetting", "check": check})
    if tier != "Rejected":
        (ev / "measurement_ir_excess.json").write_text(json.dumps({
            "channel": "ir_excess", "source_id": sid, "values": {"w4_excess_mjy": 9.5532, "t_dust_k": 629.0, "curve": [1, 2, 3], "label": "K4V"},
            "uncertainties": {"w4_excess_mjy": 0.95}, "caveats": ["E(B-V) 0.030 applied to the star (Gaia GSP-Phot)"],
            "score": {"raw_statistic": 9.87654, "false_alarm_prob": 1.2e-9, "flagged": True, "summary": "W4 excess 16.9x photosphere"}}))
        files.append({"path": "measurement_ir_excess.json", "sha256": "0" * 64, "bytes": 1, "kind": "measurement", "channel": "ir_excess"})
    (ev / "index.json").write_text(json.dumps({
        "source_id": sid, "run_id": run_id, "tier": tier, "scope": "full", "flagged_channels": ["ir_excess"], "config_hash": "abc",
        "git_commit": "deadbeef0000", "pipeline_version": "x", "timestamp_utc": "2026-11-18T00:00:00+00:00",
        "vetting": [{"check": c, "status": s, "evidence": f"vet_{c}.json"} for c, s, _ in vet], "files": files, "omitted_files": []}))
    row = {"entry_id": f"REG-{sid}-{run_id}", "timestamp_utc": "2026-11-18T00:00:00+00:00", "source_id": sid, "event_type": "candidate_registered",
           "tier": tier, "status": "registered", "provenance": {}, "details": {}, "evidence_dir": f"registry/evidence/{sid}/{run_id}",
           "evidence_index_sha256": "0" * 64}
    with (root / "registry" / "ledger.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(row) + "\n")


def _repo(tmp_path: Path) -> Path:
    _live(tmp_path, "v111-halo-ir_excess", flagged=("3", "4"))
    _live(tmp_path, "v111-halo-ir_timedomain", flagged=("3",), timedomain=True)
    rd = _fake_run(tmp_path, "fixture-run", n=9)
    sf = next(rd.glob("search_*.json"))
    data = json.loads(sf.read_text())
    data["cuts"][0]["cut_name"] = "synthetic_field"
    sf.write_text(json.dumps(data))
    export_run_record(tmp_path, "fixture-run")
    (tmp_path / "registry" / "ledger.jsonl").write_text(json.dumps({
        "entry_id": "OLD-1", "timestamp_utc": "2026-11-16T00:00:00+00:00", "source_id": "9", "event_type": "candidate_registered", "tier": "C",
        "status": "registered", "provenance": {}, "details": {}}) + "\n")
    _register(tmp_path, "3", "C", "v111-fuse")
    _register(tmp_path, "4", "Rejected", "v111-fuse")
    (tmp_path / "docs" / "channels").mkdir(parents=True)
    for ch in ("ir_excess", "ir_timedomain", "vanishing"):
        (tmp_path / "docs" / "channels" / f"{ch}.md").write_text("x")
    return tmp_path


def test_counts_come_from_run_records(tmp_path: Path) -> None:
    md = build_candidates_md(_repo(tmp_path))
    assert "6 distinct stars were in at least one live sample" in md
    assert "`ir_excess` 6, `ir_timedomain` 4" in md  # two stars of the time-domain run were not searched
    assert "2 stars were flagged by at least one channel. Registered outcome of those: 1 C, 1 Rejected." in md
    assert "2 (fetch failed 1, no coverage 1)" in md
    assert "Channels with no committed live search: `vanishing`" in md
    assert "`fixture-run`" in md and "fixture (synthetic) searches and are not counted above" in md
    assert "Of 3 registrations, 2 carry committed evidence and 1 are `legacy_without_evidence`" in md


def test_candidate_section_shows_every_check_and_its_evidence(tmp_path: Path) -> None:
    md = build_candidates_md(_repo(tmp_path))
    assert "### `Gaia DR3 3`: Tier C" in md and "#### Vetting: 1 pass, 0 fail, 1 not evaluated" in md
    assert "| `binary` | not_evaluated | Gaia RVS scatter 6.13 km/s is not constraining \\| pipes are escaped | [vet_binary.json](../registry/evidence/3/v111-fuse/vet_binary.json) |" in md
    assert "Not ruled out (each `not_evaluated` check is an open hole, never a pass): `binary`." in md
    assert "| `w4_excess_mjy` | 9.553 | 0.95 |" in md and "| `label` | K4V |  |" in md and "curve" not in md
    assert "statistic 9.877, false-alarm probability 1.2e-09, flagged. W4 excess 16.9x photosphere" in md
    assert "Caveat: E(B-V) 0.030 applied to the star (Gaia GSP-Phot)" in md
    assert "Review status: no review rounds, no model family recorded, no human checklist." in md
    assert "| `Gaia DR3 4` | `v111-fuse` | `centroid_offset`: W4 centroid 9.1\" from the star |" in md
    assert "Gaia DR3 9" not in md  # a legacy registration without evidence supports no statement


def test_reviewed_document_lists_only_surviving_tiers(tmp_path: Path) -> None:
    md = build_reviewed_md(_repo(tmp_path))
    assert "1 star(s) are registered at Tier A, B or C with committed evidence" in md
    assert "## `Gaia DR3 3`: Tier C" in md and "Gaia DR3 4" not in md
    assert "### What is not ruled out" in md and "- `binary`: Gaia RVS scatter" in md
    assert "### Checks that passed" in md and "`moving_object`" in md


def test_notes_are_included_and_labelled(tmp_path: Path) -> None:
    root = _repo(tmp_path)
    (root / "docs" / "notes").mkdir()
    (root / "docs" / "notes" / "candidates_notes.md").write_text("The W3/W4 slope needs two dust temperatures.\n")
    md = build_candidates_md(root)
    assert md.rstrip().endswith("The W3/W4 slope needs two dust temperatures.")
    assert "## Notes (written by hand, not pipeline output)" in md


def test_build_is_deterministic_and_staleness_is_detected(tmp_path: Path) -> None:
    root = _repo(tmp_path)
    assert sorted(stale_documents(root)) == ["docs/candidates.md", "reviewed_candidate.md"]
    assert sorted(write_documents(root)) == ["docs/candidates.md", "reviewed_candidate.md"]
    assert stale_documents(root) == [] and write_documents(root) == []
    p = root / "reviewed_candidate.md"
    p.write_text(p.read_text().replace("Tier C", "Tier A"))
    assert stale_documents(root) == ["reviewed_candidate.md"]
