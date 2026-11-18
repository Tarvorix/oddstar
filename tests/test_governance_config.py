"""R4 governance: config version vs CHANGELOG, thresholds moved to config unchanged, commit separation, lesson SHAs."""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from pathlib import Path

import pytest

from oddstar.ci_checks import check_config_version, check_git_commit_separation, check_lesson_commits
from oddstar.config import load_config, repo_root
from oddstar.io.live_survey import _phot_quality_where_sql, _region_where_sql
from oddstar.schemas import LessonRecord

pytestmark = pytest.mark.skipif(shutil.which("git") is None, reason="git not installed")

# Query text before the cuts moved to `live_ingestion` (copied from live_survey.py at 7bf6ed5). The raw-cache key
# includes the query, so the text must stay identical.
OLD_PHOT = """
      AND awn.number_of_neighbours = 1
      AND awn.number_of_mates = 0
      AND tmn.number_of_neighbours = 1
      AND tmn.number_of_mates = 0
      AND aw.w1mpro_error IS NOT NULL
      AND aw.w2mpro_error IS NOT NULL
      AND aw.w3mpro_error IS NOT NULL
      AND tm.j_msigcom IS NOT NULL
      AND tm.h_msigcom IS NOT NULL
      AND tm.ks_msigcom IS NOT NULL
    """
OLD_REGIONS = {
    "halo": """
          g.parallax_over_error > 10.0
          AND g.parallax BETWEEN 0.5 AND 20.0
          AND g.ruwe < 1.4
          AND ABS(g.b) > 20.0
          AND g.phot_g_mean_mag BETWEEN 7.5 AND 14.5
          AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
          AND g.mh_gspphot < -0.8
          AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 150.0
        """,
    "bulge": """
          (g.l <= 12.0 OR g.l >= 348.0)
          AND ABS(g.b) BETWEEN 3.0 AND 10.0
          AND g.parallax BETWEEN 0.08 AND 0.45
          AND g.parallax_over_error > 5.0
          AND g.ruwe < 1.4
          AND g.phot_g_mean_mag BETWEEN 8.0 AND 14.5
          AND g.teff_gspphot BETWEEN 3500.0 AND 7500.0
        """,
    "local": """
          g.parallax_over_error > 15.0
          AND g.parallax BETWEEN 2.0 AND 25.0
          AND g.ruwe < 1.4
          AND ABS(g.b) > 12.0
          AND g.phot_g_mean_mag BETWEEN 7.5 AND 13.5
          AND g.teff_gspphot BETWEEN 3800.0 AND 7500.0
          AND 4.74047 * SQRT(g.pmra * g.pmra + g.pmdec * g.pmdec) / g.parallax > 65.0
        """,
}

ENV = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
       "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid"}


def _git(root: Path, *args: str) -> str:
    return subprocess.run(["git", *args], cwd=root, env=ENV, check=True, capture_output=True, text=True).stdout.strip()


def _commit(root: Path, files: dict[str, str], msg: str) -> str:
    for rel, text in files.items():
        p = root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    _git(root, "add", "-A")
    _git(root, "commit", "-q", "-m", msg)
    return _git(root, "rev-parse", "HEAD")


def _repo(tmp: Path) -> Path:
    root = tmp / "repo"
    root.mkdir()
    _git(root, "init", "-q")
    return root


def test_live_sample_sql_is_unchanged_by_the_move_to_config() -> None:
    cfg = load_config()[0]
    assert _phot_quality_where_sql(cfg) == OLD_PHOT
    for reg, text in OLD_REGIONS.items():
        assert _region_where_sql(reg, cfg) == text, reg


def test_moved_threshold_values_are_unchanged() -> None:
    cfg = load_config()[0]
    v = cfg.vetting
    assert cfg.version == "1.4.0"
    expected = {
        "xmatch_max_sep_arcsec": 2.5, "wise_w4_min_frames": 5, "wise_w4_max_rchi2": 3.5, "dog_w1_w2_min": 0.45,
        "dog_w2_w3_min": 3.0, "lc_bls_min_period_days": 0.8, "lc_bls_max_period_days": 8.0, "lc_bls_duration_days": 0.1,
        "spectral_cr_spike_mad_sigma": 8.0, "spectral_sky_peak_mad_sigma": 6.0, "vanishing_pm_min_raw_sep_arcsec": 3.0,
        "outlier_driver_min_abs_contribution": 2.0, "binary_rv_scatter_max_kms": 5.0, "evolved_teff_max_k": 5200.0,
    }
    assert {k: v[k] for k in expected} == expected
    assert (cfg.fusion.fap_floor, cfg.fusion.fallback_false_alarm_rate_per_star) == (1e-12, 1e-4)
    li = cfg.live_ingestion
    assert li.ir_excess_w4_upper_limit_ph_qual == ["U", "X"] and li.ir_excess_artifact_cc_flag_bands == [3, 4]
    assert not hasattr(li, "ir_timedomain_neowise_enrich_max_stars")  # the 25-star cap is gone: ir_timedomain fetches every sample star


def test_remaining_default_pass_keys_and_no_literal_fallbacks_in_vetting_source() -> None:
    """Vanishing/spectral/radio vetting read decision values from config and inputs only when explicitly present."""
    cfg = load_config()[0]
    assert cfg.channels["vanishing"]["vetting_variable_simbad_types"] == ["Mira", "LPV*", "CV*", "Flare*"]
    assert cfg.channels["optical_dips"]["detrend_edge_zone_window_frac"] == 0.5
    assert cfg.vetting["metallicity_subdwarf_feh_max"] == -1.0
    # Values these modules read with literal fallbacks before 1.1.0 are config keys with unchanged values
    van, rad, spl = cfg.channels["vanishing"], cfg.channels["radio_narrowband"], cfg.channels["spectral_lines"]
    assert (van["min_early_plates"], van["min_depth_margin_mag"], van["max_pm_residual_arcsec"]) == (2, 2.0, 1.5)
    assert rad["min_drift_abs_hz_s"] == 0.05
    assert (spl["sky_line_exclusion_angstrom"], spl["stellar_line_exclusion_angstrom"]) == (1.5, 1.5)
    vet = repo_root() / "src" / "oddstar" / "vet"
    # `.get("key", <number/bool/"Star">)`: a default input value, or a config fallback that hides a missing key
    literal_default = re.compile(r"""\.get\(\s*["'][^"']+["']\s*,\s*(-?[0-9.]+|True|False|["']Star["'])\s*\)""")
    for name in ("vanishing_checks.py", "spectral_checks.py", "radio_checks.py", "common.py"):
        text = (vet / name).read_text(encoding="utf-8")
        assert not literal_default.findall(text), (name, literal_default.findall(text))
        assert "cfg.channels" not in text or not re.search(r"cfg\.channels\[[^\]]+\]\.get\(", text), name
    lc = (vet / "lc_checks.py").read_text(encoding="utf-8")
    assert 'get("simbad_type", "Star")' not in lc and '"simbad_type": "Star"' not in (repo_root() / "src" / "oddstar" / "io" / "live_survey.py").read_text(encoding="utf-8")


def test_config_version_must_match_newest_changelog_heading(tmp_path: Path) -> None:
    assert check_config_version(repo_root()) == []
    (tmp_path / "configs").mkdir()
    (tmp_path / "configs" / "default.yaml").write_text('version: "1.0.0"\n', encoding="utf-8")
    (tmp_path / "configs" / "CHANGELOG.md").write_text("# Log\n\n## [1.1.0] - x\n- a\n\n## [1.0.0] - y\n", encoding="utf-8")
    errs = check_config_version(tmp_path)
    assert len(errs) == 1 and "'1.0.0' differs from newest CHANGELOG entry '1.1.0'" in errs[0]


def test_commit_separation_checks_head_and_base_range_after_cutoff(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GITHUB_BASE_REF", raising=False)
    monkeypatch.delenv("ODDSTAR_BASE_REF", raising=False)
    root = _repo(tmp_path)
    # Before the cutoff: a mixed commit is historical and not reported
    _commit(root, {"configs/default.yaml": "a: 1\n", "registry/ledger.jsonl": "{}\n"}, "old mixed")
    cutoff = _commit(root, {"README.md": "x\n"}, "cutoff")
    notes: list[str] = []
    assert check_git_commit_separation(root, notes) == []
    assert any("no R4 cutoff recorded" in n for n in notes)

    base = _commit(root, {"configs/CHANGELOG.md": f"## [1.0.0]\nR4 commit-separation cutoff: `{cutoff}`\n"}, "record cutoff")
    notes = []
    assert check_git_commit_separation(root, notes) == []
    assert any("1 commit(s) after cutoff" in n for n in notes)

    bad = _commit(root, {"configs/default.yaml": "a: 2\n", "registry/ledger.jsonl": "{}\n{}\n"}, "mixed after cutoff")
    errs = check_git_commit_separation(root)
    assert len(errs) == 1 and bad[:10] in errs[0] and "configs/default.yaml" in errs[0]

    _commit(root, {"README.md": "y\n"}, "clean head")
    assert check_git_commit_separation(root) == []  # HEAD alone is clean
    monkeypatch.setenv("ODDSTAR_BASE_REF", base)
    errs = check_git_commit_separation(root)
    assert len(errs) == 1 and bad[:10] in errs[0]  # PR range catches the earlier commit

    # Working-tree rule unchanged: uncommitted config + registry need a CHANGELOG change too
    monkeypatch.delenv("ODDSTAR_BASE_REF")
    _git(root, "reset", "-q", "--hard", base)
    (root / "configs" / "default.yaml").write_text("a: 3\n", encoding="utf-8")
    (root / "registry" / "ledger.jsonl").write_text("{}\n{}\n{}\n", encoding="utf-8")
    assert any("without configs/CHANGELOG.md" in e for e in check_git_commit_separation(root))


def _lesson(sha: str = "") -> str:
    return LessonRecord(lesson_id="L900", finding_id="F-1-X", source_id="1", title="t", what_happened="w",
                        generic_lesson="g", automated_check="c", linked_test="tests/x.py::t", commit_sha=sha).model_dump_json()


def test_lesson_commit_sha_must_be_the_introducing_commit(tmp_path: Path) -> None:
    root = _repo(tmp_path)
    _commit(root, {"README.md": "x\n"}, "init")
    notes: list[str] = []
    intro = _commit(root, {"lessons/L900.json": _lesson()}, "add lesson")
    assert check_lesson_commits(root, notes) == []
    assert any("commit_sha not recorded yet" in n and intro[:10] in n for n in notes)

    _commit(root, {"lessons/L900.json": _lesson("0" * 40)}, "wrong sha")
    assert any("does not resolve" in e for e in check_lesson_commits(root))
    first = _git(root, "rev-list", "--max-parents=0", "HEAD")
    _commit(root, {"lessons/L900.json": _lesson(first)}, "other sha")
    assert any("not the commit that introduced" in e for e in check_lesson_commits(root))
    _commit(root, {"lessons/L900.json": _lesson(intro)}, "right sha")
    assert check_lesson_commits(root) == []

    shallow = tmp_path / "shallow"
    subprocess.run(["git", "clone", "-q", "--depth", "1", f"file://{root}", str(shallow)], env=ENV, check=True)
    notes = []
    assert check_lesson_commits(shallow, notes) == []
    assert any("shallow clone" in n and "NOT checked" in n for n in notes)


def test_repository_lesson_shas_resolve() -> None:
    root = repo_root()
    notes: list[str] = []
    errs = check_lesson_commits(root, notes)
    if any("NOT checked" in n for n in notes):
        pytest.skip("; ".join(notes))
    assert errs == []
