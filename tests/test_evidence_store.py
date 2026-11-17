"""Committed evidence store (registry/evidence): export, size guard, tamper detection, CI check, CLI."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path
from typing import Any

import pytest

from oddstar.ci_checks import check_evidence_store
from oddstar.cli import dispatch
from oddstar.config import load_config, repo_root
from oddstar.evidence import (
    canonical_json,
    export_candidate_evidence,
    reduce_large_values,
    sha256_file,
    verify_evidence,
)
from oddstar.registry import read_ledger
from oddstar.schemas import ChannelScore, EvidenceIndex, VettingResult

SID = "5628404318169472"


@pytest.fixture(scope="module")
def cfg() -> Any:
    c, _, _ = load_config()
    return c


def _workspace(tmp: Path) -> Path:
    root = tmp / "ws"
    (root / "configs").mkdir(parents=True)
    src = repo_root()
    shutil.copy(src / "configs" / "default.yaml", root / "configs" / "default.yaml")
    (root / "pyproject.toml").write_text("[project]\nname='ws'\n", encoding="utf-8")
    (root / "registry").mkdir()
    (root / "registry" / "ledger.jsonl").write_text("", encoding="utf-8")
    return root


def _png(path: Path, n_bytes: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + bytes(n_bytes))
    return path


def _score(flagged: bool = True, channel: str = "ir_excess") -> ChannelScore:
    return ChannelScore(source_id=SID, channel=channel, raw_statistic=9.0, false_alarm_prob=1e-6, flagged=flagged,
                        rank_score=6.0, summary="W4 excess")


def _vet(tmp: Path) -> list[VettingResult]:
    ev_dir = tmp / "runs" / "r1" / "vet_evidence" / SID
    ev_dir.mkdir(parents=True)
    series = [0.001 * i for i in range(5000)]
    (ev_dir / "vet_single_frames.json").write_text(json.dumps({
        "check": "single_frames", "status": "pass", "note": "19/19 frames consistent",
        "metrics": {"w4_frames_mjy": series},
        "inputs": {"manifest_ids": ["MAN-abc"], "query_urls": ["https://irsa.ipac.caltech.edu/TAP/sync?QUERY=x"]},
    }), encoding="utf-8")
    png = _png(ev_dir / "wise_cutouts.png", 2000)
    return [
        VettingResult(check="single_frames", status="pass", evidence=str(ev_dir / "vet_single_frames.json"), note="19/19 frames consistent",
                      metrics={"n_frames": 19, "w4_frames_mjy": series, "cutout_sha256": "a" * 64}),
        VettingResult(check="vision_rubric", status="pass", evidence=str(png), note="clean", metrics={"cutout_path": str(png)}),
        VettingResult(check="deep_imaging", status="not_evaluated", evidence=None, note="no sub-arcsec imaging"),
        VettingResult(check="wise_quality", status="fail", evidence=None, note="cc_flags=00D0", metrics={"cc_flags": "00D0"}),
    ]


# ------------------------------------------------------------------------------------------ size guard


def test_size_guard_reduces_long_arrays_and_strings_generally() -> None:
    series = [float(i) for i in range(1000)]
    pixels = [[float(i + j) for j in range(40)] for i in range(40)]  # 1600 scalars in a 40-row outer list
    obj = {"short": [1, 2, 3], "series": series, "nested": {"pixels": pixels}, "blob": "x" * 5000, "label": "ok"}
    red = reduce_large_values(obj, max_elements=64, max_chars=500)
    assert red["short"] == [1, 2, 3] and red["label"] == "ok"
    s = red["series"]
    assert s["_reduced"] == "array" and s["length"] == 1000 and s["n_scalars"] == 1000
    assert s["sha256"] == hashlib.sha256(canonical_json(series)).hexdigest()
    assert s["summary"]["min"] == 0.0 and s["summary"]["max"] == 999.0 and s["summary"]["median"] == 499.5
    p = red["nested"]["pixels"]
    assert p["_reduced"] == "array" and p["shape"] == [40, 40] and p["n_scalars"] == 1600
    assert red["blob"]["_reduced"] == "string" and red["blob"]["length"] == 5000
    # Non-numeric long lists keep a short head instead of a numeric summary
    words = reduce_large_values(["w"] * 100, max_elements=10, max_chars=500)
    assert words["_reduced"] == "array" and "summary" not in words and len(words["head"]) == 3


# ---------------------------------------------------------------------------------------------- export


def test_export_writes_compact_files_index_and_inputs(tmp_path: Path, cfg: Any) -> None:
    root = _workspace(tmp_path)
    vet = _vet(root)
    big = _png(root / "dossiers" / SID / "evidence" / "sed.png", int(cfg.evidence.max_image_bytes) + 10)
    small = _png(root / "dossiers" / SID / "evidence" / "stat_ir_excess.png", 500)
    meas = [{"channel": "ir_excess", "values": {"w4_ratio": 18.3}, "extra": {"sed_flux": list(range(500))}}]
    exp = export_candidate_evidence(
        root=root, cfg=cfg, source_id=f"Gaia DR3 {SID}", run_id="fuse-test", tier="C", vetting=vet, measurements=meas,
        scores=[_score(), _score(False, "generic_outlier")], images=[big, small], config_hash="abc123def456",
        git_commit="0" * 40, pipeline_version="0.1.0",
    )
    assert exp.rel_dir == f"registry/evidence/{SID}/fuse-test"
    idx = EvidenceIndex.model_validate_json((exp.directory / "index.json").read_text(encoding="utf-8"))
    assert sha256_file(exp.directory / "index.json") == exp.index_sha256
    assert idx.tier == "C" and idx.scope == "full" and idx.flagged_channels == ["ir_excess"] and idx.config_hash == "abc123def456"
    names = {f.path for f in idx.files}
    assert {"vet_single_frames.json", "vet_vision_rubric.json", "vet_deep_imaging.json", "vet_wise_quality.json",
            "measurement_ir_excess.json", "measurement_generic_outlier.json", "stat_ir_excess.png",
            "vision_rubric_wise_cutouts.png"} <= names
    assert "sed.png" not in names and idx.omitted_files[0].name == "sed.png" and idx.omitted_files[0].bytes == big.stat().st_size
    for f in idx.files:
        assert sha256_file(exp.directory / f.path) == f.sha256
        assert (exp.directory / f.path).stat().st_size <= max(int(cfg.evidence.max_json_bytes), int(cfg.evidence.max_image_bytes))

    sf = json.loads((exp.directory / "vet_single_frames.json").read_text(encoding="utf-8"))
    assert sf["status"] == "pass" and sf["metrics"]["n_frames"] == 19
    assert sf["metrics"]["w4_frames_mjy"]["_reduced"] == "array" and sf["metrics"]["w4_frames_mjy"]["length"] == 5000
    assert sf["inputs"]["manifest_ids"] == ["MAN-abc"]
    assert sf["inputs"]["query_urls"] == ["https://irsa.ipac.caltech.edu/TAP/sync?QUERY=x"]
    assert sf["inputs"]["sha256"] == ["a" * 64]
    assert sf["source_evidence"]["sha256"] == sha256_file(Path(vet[0].evidence))  # type: ignore[arg-type]
    m = json.loads((exp.directory / "measurement_ir_excess.json").read_text(encoding="utf-8"))
    assert m["values"]["w4_ratio"] == 18.3 and m["extra"]["sed_flux"]["_reduced"] == "array" and m["score"]["flagged"] is True
    rows = {r["check"]: r for r in exp.vetting_rows}
    assert rows["single_frames"]["evidence"] == "vet_single_frames.json" and rows["deep_imaging"]["status"] == "not_evaluated"
    assert verify_evidence(root)["problems"] == []


def test_rejected_exports_only_the_failing_cause(tmp_path: Path, cfg: Any) -> None:
    root = _workspace(tmp_path)
    exp = export_candidate_evidence(
        root=root, cfg=cfg, source_id=SID, run_id="fuse-rej", tier="Rejected", vetting=_vet(root),
        measurements=[{"channel": "ir_excess", "values": {"w4_ratio": 1.0}}], scores=[_score()],
        config_hash="abc", git_commit="0" * 40, pipeline_version="0.1.0",
    )
    assert exp.index.scope == "rejection_only"
    assert [f.path for f in exp.index.files] == ["vet_wise_quality.json"]
    assert exp.vetting_rows == [{"check": "wise_quality", "status": "fail", "evidence": "vet_wise_quality.json"}]
    assert exp.index.flagged_channels == ["ir_excess"]


# ------------------------------------------------------------------------------- fuse, verify, CI check


def _fuse_workspace(tmp_path: Path, *, fail_check: str | None = None) -> Path:
    root = _workspace(tmp_path)
    run = root / "runs" / "vet-run"
    run.mkdir(parents=True)
    ev = run / "vet_evidence" / SID
    ev.mkdir(parents=True)
    (ev / "vet_cirrus.json").write_text(json.dumps({"check": "cirrus", "status": "pass", "note": "E(B-V)=0.06", "metrics": {}}), encoding="utf-8")
    score = _score()
    (run / "search_ir_excess.json").write_text(json.dumps({
        "run_id": "vet-run", "channel": "ir_excess", "n_searched": 300,
        "scores": [score.model_dump()], "measurements": [{"source_id": SID, "channel": "ir_excess", "values": {"w4_ratio": 18.3}}],
        "flagged_ids": [SID],
    }), encoding="utf-8")
    rows = [
        VettingResult(check="cirrus", status="pass", evidence=str(ev / "vet_cirrus.json"), note="E(B-V)=0.06").model_dump(),
        VettingResult(check="deep_imaging", status="not_evaluated", note="no imaging").model_dump(),
    ]
    if fail_check:
        (ev / f"vet_{fail_check}.json").write_text(json.dumps({"check": fail_check, "status": "fail"}), encoding="utf-8")
        rows.append(VettingResult(check=fail_check, status="fail", evidence=str(ev / f"vet_{fail_check}.json"), note="bad").model_dump())
    (run / "vetting.json").write_text(json.dumps({"run_id": "vet-run", "vetting": {SID: rows}}), encoding="utf-8")
    return root


def test_fuse_registers_candidate_with_evidence_reference(tmp_path: Path) -> None:
    root = _fuse_workspace(tmp_path)
    # A legacy registration (no evidence reference) is counted, never failed
    legacy = {"entry_id": "REG-FUSE-old-1", "timestamp_utc": "2026-11-16T00:00:00+00:00", "source_id": "1",
              "event_type": "candidate_registered", "tier": "C", "status": "Tier C", "provenance": {}, "details": {}}
    (root / "registry" / "ledger.jsonl").write_text(json.dumps(legacy) + "\n", encoding="utf-8")
    code, pay = dispatch(["fuse", "vet-run", "--out-run-id", "fuse-t", "--json"], root=root)
    assert code == 0, pay
    entry = read_ledger(root)[-1]
    assert entry["event_type"] == "candidate_registered" and entry["tier"] == "C"
    assert entry["evidence_dir"] == f"registry/evidence/{SID}/fuse-t"
    assert entry["evidence_index_sha256"] == sha256_file(root / entry["evidence_dir"] / "index.json")
    assert {r["check"]: r["evidence"] for r in entry["details"]["vetting"]} == {"cirrus": "vet_cirrus.json", "deep_imaging": "vet_deep_imaging.json"}
    notices: list[str] = []
    assert check_evidence_store(root, notices) == []
    assert "1 ledger entries with committed evidence, 1 legacy_without_evidence" in notices[0]
    code, st = dispatch(["status", "--json"], root=root)
    assert code == 0 and st["data"]["evidence"]["legacy_without_evidence"] == 1
    assert st["data"]["evidence"]["ledger_entries_with_evidence"] == 1
    code, ver = dispatch(["evidence", "verify", SID, "--json"], root=root)
    assert code == 0 and ver["data"]["problems"] == []

    # Re-fusing the same run id never overwrites the referenced directory
    code, _ = dispatch(["fuse", "vet-run", "--out-run-id", "fuse-t", "--json"], root=root)
    assert read_ledger(root)[-1]["evidence_dir"] == f"registry/evidence/{SID}/fuse-t.2"
    assert check_evidence_store(root) == []


def test_fuse_rejected_registration_exports_failing_check_only(tmp_path: Path) -> None:
    root = _fuse_workspace(tmp_path, fail_check="wise_quality")
    code, pay = dispatch(["fuse", "vet-run", "--out-run-id", "fuse-r", "--json"], root=root)
    assert code == 0, pay
    entry = read_ledger(root)[-1]
    assert entry["tier"] == "Rejected" and entry["details"]["evidence_scope"] == "rejection_only"
    assert [r["check"] for r in entry["details"]["vetting"]] == ["wise_quality"]
    assert sorted(p.name for p in (root / entry["evidence_dir"]).iterdir()) == ["index.json", "vet_wise_quality.json"]
    assert check_evidence_store(root) == []


def test_verification_fails_on_tampered_or_missing_evidence(tmp_path: Path) -> None:
    root = _fuse_workspace(tmp_path)
    assert dispatch(["fuse", "vet-run", "--out-run-id", "fuse-x", "--json"], root=root)[0] == 0
    d = root / "registry" / "evidence" / SID / "fuse-x"
    tampered = d / "vet_cirrus.json"
    obj = json.loads(tampered.read_text(encoding="utf-8"))
    obj["status"] = "fail"
    tampered.write_text(json.dumps(obj), encoding="utf-8")
    problems = check_evidence_store(root)
    assert any("vet_cirrus.json: SHA-256 mismatch" in p for p in problems)
    code, err = dispatch(["evidence", "verify", SID, "--json"], root=root)
    assert code != 0 and err["ok"] is False and "SHA-256 mismatch" in err["message"]

    # Editing index.json to match the tampered file breaks the ledger reference instead
    idx = EvidenceIndex.model_validate_json((d / "index.json").read_text(encoding="utf-8"))
    for f in idx.files:
        if f.path == "vet_cirrus.json":
            f.sha256 = sha256_file(tampered)
    (d / "index.json").write_text(idx.model_dump_json(indent=2), encoding="utf-8")
    assert any("index.json SHA-256 differs from the ledger reference" in p for p in check_evidence_store(root))

    # A pass row whose evidence file is gone fails the CI check
    shutil.copy(d / "index.json", tmp_path / "idx.bak")
    tampered.unlink()
    problems = check_evidence_store(root)
    assert any("pass row 'cirrus' has no evidence file" in p for p in problems)
    assert any("listed in index.json but missing" in p for p in problems)


def test_ci_check_requires_vetting_rows_and_files(tmp_path: Path) -> None:
    root = _fuse_workspace(tmp_path)
    assert dispatch(["fuse", "vet-run", "--out-run-id", "fuse-y", "--json"], root=root)[0] == 0
    rows = read_ledger(root)
    forged = dict(rows[-1], entry_id="REG-forged", details={**rows[-1]["details"], "vetting": [
        {"check": "neighbours", "status": "pass", "evidence": None}]})
    with (root / "registry" / "ledger.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(forged) + "\n")
    problems = check_evidence_store(root)
    assert any("pass row 'neighbours' has no evidence file" in p for p in problems)


def test_cli_evidence_export_from_vet_run(tmp_path: Path) -> None:
    root = _fuse_workspace(tmp_path)
    code, pay = dispatch(["evidence", "export", SID, "--run-id", "vet-run", "--json"], root=root)
    assert code == 0, pay
    assert pay["data"]["evidence_dir"] == f"registry/evidence/{SID}/vet-run"
    assert read_ledger(root) == []  # export alone never writes the ledger
    assert verify_evidence(root)["problems"] == []
    code, err = dispatch(["evidence", "export", "123", "--json"], root=root)
    assert code != 0 and err["ok"] is False


def test_inspect_register_exports_dossier_evidence(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    code, pay = dispatch(["inspect", "4582145718700356096", "--json"], root=root)
    assert code == 0 and read_ledger(root) == []  # plain inspect does not register
    code, pay = dispatch(["inspect", "4582145718700356096", "--register", "--json"], root=root)
    assert code == 0, pay
    entry = read_ledger(root)[-1]
    assert entry["tier"] == "C" and entry["provenance"]["mode"] == "inspect"
    d = root / entry["evidence_dir"]
    idx = EvidenceIndex.model_validate_json((d / "index.json").read_text(encoding="utf-8"))
    assert "ir_excess" in idx.flagged_channels
    assert any(f.kind == "image" for f in idx.files) or idx.omitted_files
    assert {r["check"] for r in entry["details"]["vetting"] if r["status"] in ("pass", "fail")} <= {f.check for f in idx.files if f.check}
    assert check_evidence_store(root) == []
