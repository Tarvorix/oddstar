"""Certificates carry a code fingerprint; search refuses a stale one and CI reports stale committed certificates."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Any

import pytest

from oddstar.ci_certificates import check_certificates_current
from oddstar.cli import dispatch
from oddstar.config import load_config, repo_root
from oddstar.errors import EXIT_UNCERTIFIED, UncertifiedSearchError
from oddstar.schemas import ValidationCertificate
from oddstar.validate import (
    all_certificate_status,
    certificate_status,
    channel_source_files,
    compute_code_fingerprint,
    require_valid_certificate,
    run_validation,
)
from oddstar.validate.fingerprint import default_source_root


def _workspace(tmp_path: Path, *, with_src: bool = False) -> Path:
    root = repo_root()
    ws = tmp_path / "ws"
    (ws / "configs").mkdir(parents=True)
    shutil.copy(root / "configs" / "default.yaml", ws / "configs" / "default.yaml")
    for sub in ("benchmarks", "impostors"):
        (ws / "validation" / sub).mkdir(parents=True)
        shutil.copy(root / "validation" / sub / f"{sub}.json", ws / "validation" / sub / f"{sub}.json")
    if with_src:
        shutil.copytree(default_source_root() / "oddstar", ws / "src" / "oddstar", ignore=shutil.ignore_patterns("__pycache__"))
    return ws


def test_file_list_is_the_channel_closure_plus_validation_machinery() -> None:
    files = channel_source_files("generic_outlier")
    assert files == sorted(files)
    for must in (
        "oddstar/channels/generic_outlier/channel.py",
        "oddstar/channels/underluminous/channel.py",  # imported by generic_outlier itself
        "oddstar/io/physics.py",
        "oddstar/io/fixtures_data.py",
        "oddstar/validate/harness.py",
        "oddstar/vet/engine.py",
        "oddstar/fuse/engine.py",
        "oddstar/channels/__init__.py",
    ):
        assert must in files, must
    for never in ("oddstar/channels/optical_dips/channel.py", "oddstar/channels/ir_excess/channel.py", "oddstar/cli.py",
                  "oddstar/ci_checks.py", "oddstar/io/live_survey.py"):
        assert never not in files, never
    assert "oddstar/channels/ir_excess/channel.py" in channel_source_files("ir_excess")
    assert "oddstar/channels/underluminous/channel.py" not in channel_source_files("vanishing")


def test_fingerprint_is_deterministic_line_ending_neutral_and_channel_local(tmp_path: Path) -> None:
    src = tmp_path / "src"
    shutil.copytree(default_source_root() / "oddstar", src / "oddstar", ignore=shutil.ignore_patterns("__pycache__"))
    fp_go, files_go = compute_code_fingerprint("generic_outlier", src)
    assert fp_go == compute_code_fingerprint("generic_outlier")[0]
    fp_od = compute_code_fingerprint("optical_dips", src)[0]

    physics = src / "oddstar" / "io" / "physics.py"
    physics.write_bytes(physics.read_bytes().replace(b"\n", b"\r\n"))
    assert compute_code_fingerprint("generic_outlier", src)[0] == fp_go

    od = src / "oddstar" / "channels" / "optical_dips" / "channel.py"
    od.write_text(od.read_text(encoding="utf-8") + "\n# behaviour change\n", encoding="utf-8")
    assert compute_code_fingerprint("generic_outlier", src)[0] == fp_go
    assert compute_code_fingerprint("optical_dips", src)[0] != fp_od

    ul = src / "oddstar" / "channels" / "underluminous" / "channel.py"
    ul.write_text(ul.read_text(encoding="utf-8") + "\nX = 1\n", encoding="utf-8")
    fp2, files2 = compute_code_fingerprint("generic_outlier", src)
    assert fp2 != fp_go
    assert [k for k in files2 if files2[k] != files_go.get(k)] == ["oddstar/channels/underluminous/channel.py"]


def test_search_is_blocked_when_the_certificate_code_fingerprint_is_stale(tmp_path: Path) -> None:
    ws = _workspace(tmp_path)
    cfg, short_hash, _ = load_config(ws / "configs" / "default.yaml", root=ws)
    cert, cert_path = run_validation("vanishing", cfg, ws, split="all", field_limit=300)
    assert cert.code_fingerprint == compute_code_fingerprint("vanishing")[0]
    assert set(cert.details["code_fingerprint_files"]) == set(channel_source_files("vanishing"))
    assert require_valid_certificate("vanishing", cfg, ws).code_fingerprint == cert.code_fingerprint
    assert certificate_status("vanishing", cfg, ws)["usable"] is True
    code, _ = dispatch(["search", "vanishing", "--limit", "50", "--json"], root=ws)
    assert code == 0

    stale = json.loads(cert_path.read_text(encoding="utf-8"))
    stale["code_fingerprint"] = "0" * 64
    stale["details"]["code_fingerprint_files"]["oddstar/channels/vanishing/channel.py"] = "f" * 64
    cert_path.write_text(json.dumps(stale), encoding="utf-8")
    with pytest.raises(UncertifiedSearchError, match="stale.*changed oddstar/channels/vanishing/channel.py"):
        require_valid_certificate("vanishing", cfg, ws)
    code, payload = dispatch(["search", "vanishing", "--limit", "50", "--json"], root=ws)
    assert code == EXIT_UNCERTIFIED and payload["code"] == "UNCERTIFIED_CHANNEL"
    assert "stale" in payload["message"] and "validate vanishing" in payload["hint"]
    st = certificate_status("vanishing", cfg, ws)
    assert st["exists"] and st["passed"] and st["stale"] and not st["usable"]

    legacy = dict(stale)
    legacy.pop("code_fingerprint")
    cert_path.write_text(json.dumps(legacy), encoding="utf-8")
    with pytest.raises(UncertifiedSearchError, match="no code_fingerprint"):
        require_valid_certificate("vanishing", cfg, ws)

    rows = {r["channel"]: r for r in all_certificate_status(cfg, ws)}
    assert len(rows) == 8 and rows["ir_excess"]["exists"] is False and rows["vanishing"]["stale"] is True


def test_ci_reports_committed_certificates_that_do_not_match_the_working_tree(tmp_path: Path) -> None:
    ws = _workspace(tmp_path, with_src=True)
    cfg, _, _ = load_config(ws / "configs" / "default.yaml", root=ws)
    _, cert_path = run_validation("vanishing", cfg, ws, split="all", field_limit=300)
    # the comparison report next to the certificates is not a certificate
    (cert_path.parent / "dr3_vs_dr4_comparison.json").write_text("{}", encoding="utf-8")
    assert check_certificates_current(ws) == []

    ch_file = ws / "src" / "oddstar" / "channels" / "vanishing" / "channel.py"
    ch_file.write_text(ch_file.read_text(encoding="utf-8") + "\n# edited\n", encoding="utf-8")
    out = check_certificates_current(ws)
    assert len(out) == 1 and "code_fingerprint" in out[0] and "changed oddstar/channels/vanishing/channel.py" in out[0]

    ch_file.write_text(ch_file.read_text(encoding="utf-8").replace("\n# edited\n", ""), encoding="utf-8")
    cfg_file = ws / "configs" / "default.yaml"
    cfg_file.write_text(cfg_file.read_text(encoding="utf-8").replace("min_early_plates: 2", "min_early_plates: 3"), encoding="utf-8")
    out = check_certificates_current(ws)
    assert len(out) == 1 and "config_hash" in out[0]

    data: dict[str, Any] = json.loads(cert_path.read_text(encoding="utf-8"))
    data.pop("code_fingerprint")
    cert_path.write_text(json.dumps(data), encoding="utf-8")
    ValidationCertificate.model_validate(data)
    assert any("no code_fingerprint" in v for v in check_certificates_current(ws))
