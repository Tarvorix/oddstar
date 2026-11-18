"""Committed run records: a search run's coverage and findings must be checkable from the repository alone."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from oddstar.ci_checks import check_run_records
from oddstar.cli import dispatch
from oddstar.config import repo_root
from oddstar.errors import PolicyViolationError
from oddstar.run_records import build_run_record, export_run_record, source_ids_sha256, verify_run_records


def _fake_run(root: Path, run_id: str, *, n: int = 6, flagged: tuple[str, ...] = ("3",), timedomain: bool = False) -> Path:
    rd = root / "runs" / run_id
    rd.mkdir(parents=True)
    ids = [str(i) for i in range(1, n + 1)]
    channel = "ir_timedomain" if timedomain else "ir_excess"
    payload: dict = {
        "run_id": run_id, "channel": channel, "gaia_release": "dr3", "config_hash": "abc", "git_commit": "deadbeef", "seed": 7,
        "certificate": f"validation/certificates/{channel}-dr3-abc.json",
        "cuts": [{"cut_name": "gaia_halo_live", "description": "test", "count_before": n, "count_after": n}],
        "n_searched": n, "n_flagged": len(flagged), "flagged_ids": list(flagged),
        "scores": [{"source_id": i, "channel": channel, "raw_statistic": 9.0 if i in flagged else 0.5, "false_alarm_prob": 0.001,
                    "flagged": i in flagged, "rank_score": 1.0, "summary": "s", "feature_contributions": {}, "evidence": []} for i in ids],
        "measurements": [{"source_id": i, "big": list(range(50))} for i in ids],
    }
    if timedomain:
        payload.update({"n_input": n, "n_searched": n - 2, "n_not_searched": 2, "n_no_coverage": 1, "n_failed": 1, "n_not_fetched": 0,
                        "n_insufficient_visits": 0, "not_searched_ids": ["1", "2"],
                        "neowise_fetch": {"n_ok": n - 2, "n_failed": 1, "n_no_coverage": 1, "workers": 4},
                        "sample_stars": [{"source_id": i, "ra_deg": 1.0} for i in ids]})
    (rd / f"search_{channel}.json").write_text(json.dumps(payload), encoding="utf-8")
    (rd / "run_state.json").write_text(json.dumps({"run_id": run_id, "timestamp_utc": "2026-11-17T00:00:00+00:00"}), encoding="utf-8")
    (rd / "vetting.json").write_text(json.dumps({"vetting": {sid: [
        {"check": "moving_object", "status": "pass", "note": "x " * 400, "evidence": "e.json"},
        {"check": "binary", "status": "not_evaluated", "note": "no RV"}] for sid in flagged}}), encoding="utf-8")
    (rd / f"null_report_{channel}.json").write_text(json.dumps({
        "n_searched": payload["n_searched"], "n_flagged": len(flagged), "expected_false_alarms": 0.4, "tier_breakdown": {"C": 1},
        "upper_limits": [{"x": 1}], "contour_50": [], "contour_90": [], "completeness_grid": [{"huge": 1}] * 100}), encoding="utf-8")
    return rd


def test_record_flattens_not_searched_ids_grouped_by_reason(tmp_path: Path) -> None:
    """A real search payload groups the not-searched ids by reason; the record must list the ids, not the reasons."""
    rd = _fake_run(tmp_path, "v-test-grouped", timedomain=True)
    sf = rd / "search_ir_timedomain.json"
    payload = json.loads(sf.read_text(encoding="utf-8"))
    payload["not_searched_ids"] = {"no_coverage": ["2"], "failed": ["1"]}
    payload["tess_fetch"] = {"n_ok": 4, "n_failed": 1}
    sf.write_text(json.dumps(payload), encoding="utf-8")
    rec, _ = build_run_record(rd)
    s = rec["searches"][0]
    assert s["not_searched_ids"] == ["1", "2"] and s["tess_fetch"] == {"n_ok": 4, "n_failed": 1}
    assert "n_no_light_curve" in s["coverage"]
    export_run_record(tmp_path, "v-test-grouped")
    assert verify_run_records(tmp_path)["problems"] == []


def test_record_keeps_counts_flagged_stars_and_vetting_without_measurements(tmp_path: Path) -> None:
    rd = _fake_run(tmp_path, "v-test-itd", timedomain=True)
    rec, samples = build_run_record(rd)
    s = rec["searches"][0]
    assert rec["run_id"] == "v-test-itd" and s["channel"] == "ir_timedomain"
    assert s["sample"]["n"] == 6 and s["n_searched"] == 4 and s["coverage"]["n_failed"] == 1 and s["not_searched_ids"] == ["1", "2"]
    assert [f["source_id"] for f in s["flagged"]] == ["3"] and s["flagged"][0]["raw_statistic"] == 9.0
    assert s["vetting"]["3"][1] == {"check": "binary", "status": "not_evaluated", "note": "no RV"}
    assert len(s["vetting"]["3"][0]["note"]) <= 400
    assert s["null_report"]["expected_false_alarms"] == 0.4 and "completeness_grid" not in s["null_report"]
    assert "measurements" not in json.dumps(rec) and "big" not in json.dumps(rec)
    assert samples[0]["source_ids"] == ["1", "2", "3", "4", "5", "6"]
    assert samples[0]["source_ids_sha256"] == source_ids_sha256(["6", "5", "4", "3", "2", "1"]) == s["sample"]["source_ids_sha256"]


def test_export_shares_one_sample_file_and_verifies(tmp_path: Path) -> None:
    _fake_run(tmp_path, "v-a")
    _fake_run(tmp_path, "v-b", flagged=())
    out_a = export_run_record(tmp_path, "v-a")
    out_b = export_run_record(tmp_path, "v-b")
    assert len(list((tmp_path / "registry" / "samples").glob("*.json"))) == 1  # same six stars, stored once
    assert len(out_a["written"]) == 2 and len(out_b["written"]) == 1
    rep = verify_run_records(tmp_path)
    assert rep == {"n_records": 2, "n_samples": 1, "problems": []}
    assert export_run_record(tmp_path, "v-a")["record_sha256"] == out_a["record_sha256"]  # idempotent


def test_changed_run_needs_replace(tmp_path: Path) -> None:
    rd = _fake_run(tmp_path, "v-a")
    export_run_record(tmp_path, "v-a")
    (rd / "vetting.json").unlink()
    with pytest.raises(PolicyViolationError, match="different content"):
        export_run_record(tmp_path, "v-a")
    export_run_record(tmp_path, "v-a", replace=True)
    assert "vetting" not in json.loads((tmp_path / "registry" / "runs" / "v-a.json").read_text())["searches"][0]


def test_verification_catches_edited_records_and_samples(tmp_path: Path) -> None:
    _fake_run(tmp_path, "v-a", timedomain=True)
    export_run_record(tmp_path, "v-a")
    rec_p = tmp_path / "registry" / "runs" / "v-a.json"
    good = rec_p.read_text()
    rec = json.loads(good)

    rec["searches"][0]["n_searched"] = 6  # claim the two unsearched stars were searched
    rec_p.write_text(json.dumps(rec))
    probs = verify_run_records(tmp_path)["problems"]
    assert any("coverage.n_searched differs" in p for p in probs)

    rec_p.write_text(good)
    sample_p = next((tmp_path / "registry" / "samples").glob("*.json"))
    sample = json.loads(sample_p.read_text())
    sample["source_ids"].append("999")
    sample_p.write_text(json.dumps(sample))
    probs = verify_run_records(tmp_path)["problems"]
    assert any("does not hash" in p for p in probs)
    assert check_run_records(tmp_path) == probs

    sample_p.unlink()
    assert any("is missing" in p for p in verify_run_records(tmp_path)["problems"])
    rec_p.unlink()
    shutil.rmtree(tmp_path / "registry" / "samples")
    assert verify_run_records(tmp_path) == {"n_records": 0, "n_samples": 0, "problems": []}


def test_fusion_record_and_orphan_sample(tmp_path: Path) -> None:
    rd = tmp_path / "runs" / "fuse-x"
    rd.mkdir(parents=True)
    (rd / "fusion.json").write_text(json.dumps({
        "run_id": "fuse-x", "timestamp_utc": "t", "config_hash": "abc", "git_commit": "deadbeef", "n_stars_total": 6, "n_trials_total": 12,
        "expected_false_positives": 0.1, "consistent_with_noise": True, "cross_channel_patterns": [],
        "candidates": [{"source_id": "3", "tier": "C", "rank_score": 1.0, "combined_fap": 0.01, "flagged_channels": ["ir_excess"],
                        "patterns": [], "open_holes": ["binary"], "explained_reason": None},
                       {"source_id": "4", "tier": "Rejected", "rank_score": 0.1, "combined_fap": 0.5, "flagged_channels": ["ir_excess"],
                        "patterns": [], "open_holes": [], "explained_reason": "moving_object"}]}), encoding="utf-8")
    (rd / "fuse_inputs.json").write_text(json.dumps({"run_id": "fuse-x", "input_run_ids": ["v-a", "v-b"]}), encoding="utf-8")
    out = export_run_record(tmp_path, "fuse-x")
    fus = out["record"]["fusion"]
    assert fus["tier_counts"] == {"C": 1, "Rejected": 1} and fus["input_run_ids"] == ["v-a", "v-b"] and fus["candidates"][0]["open_holes"] == ["binary"]
    assert verify_run_records(tmp_path)["problems"] == []
    (tmp_path / "registry" / "samples").mkdir()
    (tmp_path / "registry" / "samples" / "0000000000000000.json").write_text("{}")
    assert any("not referenced" in p for p in verify_run_records(tmp_path)["problems"])


def test_cli_exports_a_real_search_run_and_the_repository_records_verify() -> None:
    root = repo_root()
    code, payload = dispatch(["search", "ir_excess", "--run-id", "test-run-record", "--seed", "20261116", "--limit", "60", "--json"])
    assert code == 0, payload
    code, payload = dispatch(["runs", "export", "test-run-record", "--dry-run", "--json"])
    assert code == 0 and "Would export 1" in payload["summary"]
    rec, samples = build_run_record(root / "runs" / "test-run-record")
    assert rec["searches"][0]["n_searched"] == samples[0]["n"] and rec["searches"][0]["search_payload_sha256"]
    code, payload = dispatch(["runs", "verify", "--json"])
    assert code == 0, payload
    assert check_run_records(root) == []
