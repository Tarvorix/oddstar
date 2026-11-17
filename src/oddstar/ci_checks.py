"""CI and Ground Rule governance checks (R4, R6, R7, AGENTS.md)."""

from __future__ import annotations

import ast
import json
import os
import re
import subprocess
from pathlib import Path

import yaml

from oddstar.ci_certificates import check_certificates_current
from oddstar.errors import PolicyViolationError
from oddstar.schemas import LessonRecord, ReviewRecord, SCHEMA_MODELS

FORBIDDEN_TERMS = [
    r"\btechnosignature\s+detected\b",
    r"\bmegastructure\b",
    r"\balien\b",
    r"\bextraterrestrial\s+intelligence\b",
]


def check_forbidden_nomenclature(root: Path) -> list[str]:
    """Scan source, docs, lessons, and dossiers for AGENTS.md forbidden terms."""
    violations: list[str] = []
    scan_dirs = ["src", "review_kit", "docs", "lessons", "dossiers"]
    for d in scan_dirs:
        base = root / d
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if not p.is_file() or p.suffix not in {".py", ".md", ".json", ".yaml", ".yml"}:
                continue
            if p.name == "ci_checks.py":
                continue
            text = p.read_text(encoding="utf-8", errors="ignore")
            for pat in FORBIDDEN_TERMS:
                for m in re.finditer(pat, text, flags=re.IGNORECASE):
                    violations.append(f"{p.relative_to(root)}: matched forbidden term '{m.group(0)}'")
    return violations


def check_review_kit_isolation(root: Path) -> list[str]:
    """Enforce R6: review_kit must never import oddstar.channels or oddstar.vet or reference data/runs."""
    violations: list[str] = []
    rk = root / "review_kit"
    if not rk.exists():
        return ["review_kit directory missing"]
    for py_file in rk.rglob("*.py"):
        src = py_file.read_text(encoding="utf-8")
        tree = ast.parse(src, filename=str(py_file))
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    if alias.name.startswith("oddstar.channels") or alias.name.startswith("oddstar.vet"):
                        violations.append(f"{py_file.name}: forbidden import {alias.name}")
            elif isinstance(node, ast.ImportFrom):
                mod = node.module or ""
                if mod.startswith("oddstar.channels") or mod.startswith("oddstar.vet"):
                    violations.append(f"{py_file.name}: forbidden import from {mod}")
    return violations


def check_lessons_and_findings(root: Path) -> list[str]:
    """Enforce R7: every closed finding in reviews/registry must link an existing test/check or reason."""
    violations: list[str] = []
    registry_file = root / "registry" / "ledger.jsonl"
    lessons_dir = root / "lessons"
    if registry_file.exists():
        for idx, line in enumerate(registry_file.read_text(encoding="utf-8").splitlines(), start=1):
            if not line.strip():
                continue
            row = json.loads(line)
            if row.get("event_type") == "finding_closed":
                details = row.get("details", {})
                linked_test = details.get("linked_test")
                reason = details.get("unautomated_reason")
                if not linked_test and not reason:
                    violations.append(f"registry line {idx}: closed finding without linked_test or reason")
                elif linked_test:
                    test_path = root / linked_test.split("::")[0]
                    if not test_path.exists():
                        violations.append(f"registry line {idx}: linked_test missing at {test_path}")

    if lessons_dir.exists():
        for lf in lessons_dir.glob("*.json"):
            rec = LessonRecord.model_validate_json(lf.read_text(encoding="utf-8"))
            for t in [rec.linked_test, *rec.additional_linked_tests]:
                test_file = root / t.split("::")[0]
                if not test_file.exists():
                    violations.append(f"{lf.name}: linked_test {t} does not exist")

    reviews_dir = root / "review" / "records"
    if reviews_dir.exists():
        for rf in reviews_dir.glob("*.json"):
            rr = ReviewRecord.model_validate_json(rf.read_text(encoding="utf-8"))
            for f in rr.findings:
                if f.status == "closed" and not (f.linked_test or f.unautomated_reason):
                    violations.append(f"{rf.name}: finding {f.finding_id} closed without linked_test or reason")
    return violations


def check_schemas_up_to_date(root: Path) -> list[str]:
    """Ensure generated JSON Schemas on disk match current Pydantic models."""
    violations: list[str] = []
    schema_dir = root / "schemas"
    for name, model_cls in SCHEMA_MODELS.items():
        path = schema_dir / f"{name}.schema.json"
        if not path.exists():
            violations.append(f"Missing schema file: {path.relative_to(root)}")
            continue
        expected = json.dumps(model_cls.model_json_schema(), indent=2, sort_keys=True) + "\n"
        actual = path.read_text(encoding="utf-8")
        if actual != expected:
            violations.append(f"Schema drift detected in {path.relative_to(root)}")
    return violations


R4_THRESHOLD_PREFIXES = ("configs/default.yaml", "validation/benchmarks/", "validation/impostors/")
R4_CUTOFF_RE = re.compile(r"R4 commit-separation cutoff[^`\n]*`([0-9a-f]{7,40})`", flags=re.IGNORECASE)


def _git(root: Path, *args: str, strip: bool = True) -> tuple[int, str]:
    try:
        res = subprocess.run(["git", *args], cwd=str(root), capture_output=True, text=True, check=False)
    except OSError:
        return 127, ""
    return res.returncode, res.stdout.strip() if strip else res.stdout


def _r4_pair(files: list[str]) -> tuple[list[str], list[str]]:
    thr = [f for f in files if f.startswith(R4_THRESHOLD_PREFIXES)]
    reg = [f for f in files if f.startswith("registry/") and f.endswith(".jsonl")]
    return thr, reg


def r4_cutoff_commit(root: Path) -> str | None:
    """The recorded cutoff (configs/CHANGELOG.md) after which committed history must keep R4 separation."""
    cl = root / "configs" / "CHANGELOG.md"
    if not cl.exists():
        return None
    m = R4_CUTOFF_RE.search(cl.read_text(encoding="utf-8"))
    return m.group(1) if m else None


def check_git_commit_separation(root: Path, notices: list[str] | None = None) -> list[str]:
    """Enforce R4: thresholds/benchmarks/impostors must not change in the same commit as registry results.

    Three scopes: (1) the uncommitted working tree (unchanged rule: allowed only together with a
    configs/CHANGELOG.md change); (2) the files changed by HEAD (for a merge HEAD, each non-merge commit it brings in);
    (3) in CI for pull requests (`GITHUB_BASE_REF`, or `ODDSTAR_BASE_REF`), each non-merge commit in base..HEAD.
    Committed history is checked only after the cutoff commit recorded in configs/CHANGELOG.md, since earlier history
    is not rewritten. A committed violation has no CHANGELOG exemption.
    """
    notes = notices if notices is not None else []
    violations: list[str] = []
    rc, out = _git(root, "status", "--porcelain", strip=False)  # leading status column is significant
    if rc != 0:
        notes.append("r4_separation: not a git work tree; commit separation not checked")
        return []
    files = [line[3:].strip() for line in out.splitlines() if len(line) > 3]
    thr, reg = _r4_pair(files)
    if thr and reg and (root / "configs" / "CHANGELOG.md").exists():
        # Uncommitted tree: allowed only when CHANGELOG is modified as well (unchanged rule)
        if "configs/CHANGELOG.md" not in files:
            violations.append("R4 violation: config/benchmark modified alongside registry without configs/CHANGELOG.md")

    rc, head = _git(root, "rev-parse", "--verify", "HEAD^{commit}")
    if rc != 0:
        return violations
    cutoff_ref = r4_cutoff_commit(root)
    rc, cutoff = _git(root, "rev-parse", "--verify", f"{cutoff_ref}^{{commit}}") if cutoff_ref else (1, "")
    if rc != 0:
        shallow = _git(root, "rev-parse", "--is-shallow-repository")[1] == "true"
        notes.append(
            "r4_separation: committed history NOT checked: "
            + (f"cutoff commit {cutoff_ref} not in history" + (" (shallow clone; use fetch-depth: 0)" if shallow else "")
               if cutoff_ref else "no R4 cutoff recorded in configs/CHANGELOG.md")
        )
        return violations

    commits: list[str] = []
    rc, parents = _git(root, "rev-list", "--parents", "-n", "1", "HEAD")
    if len(parents.split()) > 2:
        commits += _git(root, "rev-list", "--no-merges", "HEAD^1..HEAD")[1].split()
    else:
        commits.append(head)
    base = os.environ.get("ODDSTAR_BASE_REF") or (f"origin/{os.environ['GITHUB_BASE_REF']}" if os.environ.get("GITHUB_BASE_REF") else "")
    if base:
        rc, rng = _git(root, "rev-list", "--no-merges", f"{base}..HEAD")
        if rc == 0:
            commits += rng.split()
        else:
            notes.append(f"r4_separation: base ref {base} not available; only HEAD checked")
    seen: set[str] = set()
    checked = 0
    for c in commits:
        if c in seen:
            continue
        seen.add(c)
        if c == cutoff or _git(root, "merge-base", "--is-ancestor", c, cutoff)[0] == 0:
            continue  # at or before the cutoff: history is not rewritten
        checked += 1
        rc, changed = _git(root, "diff-tree", "--no-commit-id", "--name-only", "-r", "--root", c)
        thr, reg = _r4_pair(changed.split())
        if thr and reg:
            violations.append(
                f"R4 violation: commit {c[:10]} changes {', '.join(thr)} together with {', '.join(reg)} "
                "(commit thresholds/benchmarks/impostors and registry results separately)"
            )
    notes.append(f"r4_separation: {checked} commit(s) after cutoff {cutoff[:10]} checked")
    return violations


# ------------------------------------------------------------------------------------- evidence store (A)


def check_evidence_store(root: Path, notices: list[str] | None = None) -> list[str]:
    """Every ledger entry referencing evidence must match its directory and hashes; pass/fail rows need files.

    Registrations written before the evidence store are counted as legacy_without_evidence and never fail.
    """
    from oddstar.evidence import verify_evidence

    rep = verify_evidence(root)
    if notices is not None:
        notices.append(
            f"evidence: {rep['ledger_entries_with_evidence']} ledger entries with committed evidence, "
            f"{rep['legacy_without_evidence']} legacy_without_evidence (registered before the evidence store; "
            "the evidence behind their pass rows is not in the repository)"
        )
    return list(rep["problems"])


def check_run_records(root: Path, notices: list[str] | None = None) -> list[str]:
    """Every committed run record matches its sample file (ids hash and count) and its own counts."""
    from oddstar.run_records import verify_run_records

    rep = verify_run_records(root)
    if notices is not None:
        notices.append(f"run_records: {rep['n_records']} run record(s), {rep['n_samples']} sample file(s)")
    return list(rep["problems"])


# ------------------------------------------------------------------------------------ review ingestion (B)


def check_review_ingestion(root: Path) -> list[str]:
    """Every stored review record is ingested; every open record finding is opened and not silently closed."""
    from oddstar.registry import read_ledger

    violations: list[str] = []
    rec_dir = root / "review" / "records"
    if not rec_dir.exists():
        return violations
    rows = read_ledger(root)
    ingested_ids = {str(r.get("entry_id")) for r in rows if r.get("event_type") == "review_ingested"}
    ingested_paths = {str((r.get("details") or {}).get("review_path")) for r in rows if r.get("event_type") == "review_ingested"}
    lessons: dict[str, LessonRecord] = {}
    for lf in (root / "lessons").glob("*.json") if (root / "lessons").exists() else []:
        rec = LessonRecord.model_validate_json(lf.read_text(encoding="utf-8"))
        lessons[rec.lesson_id] = rec
    events: dict[str, list[tuple[int, dict]]] = {}
    for n, r in enumerate(rows):
        fid = (r.get("details") or {}).get("finding_id")
        if fid and r.get("event_type") in ("finding_opened", "finding_closed"):
            events.setdefault(fid, []).append((n, r))
    for rf in sorted(rec_dir.glob("*.json")):
        rr = ReviewRecord.model_validate_json(rf.read_text(encoding="utf-8"))
        rel = str(rf.relative_to(root))
        if f"REG-{rr.review_id}" not in ingested_ids and rel not in ingested_paths:
            violations.append(f"{rel}: no review_ingested ledger event (run `oddstar review ingest --all`)")
        for f in rr.findings:
            if f.status != "open":
                continue
            ev = events.get(f.finding_id, [])
            opened = [n for n, r in ev if r.get("event_type") == "finding_opened"]
            if not opened:
                violations.append(f"{rel}: open finding {f.finding_id} has no finding_opened ledger event")
                continue
            for n, r in ev:
                if r.get("event_type") != "finding_closed" or n < opened[0]:
                    continue
                lid = (r.get("provenance") or {}).get("lesson_id")
                lesson = lessons.get(str(lid))
                if lesson is None or lesson.finding_id != f.finding_id or not lesson.closes_finding:
                    violations.append(
                        f"{rel}: open finding {f.finding_id} closed in the ledger (line {n + 1}) without a lesson that closes it"
                    )
    return violations


# ---------------------------------------------------------------------------------- lesson metadata (C)


def check_lesson_commits(root: Path, notices: list[str] | None = None) -> list[str]:
    """Each lesson `commit_sha` must resolve in git history and be the commit that introduced the lesson files.

    Skips with a visible notice when history is unavailable or truncated (shallow clone). An empty `commit_sha` is
    allowed while the commit that adds the lesson does not exist yet (a commit cannot contain its own SHA); once it
    exists the notice names it for `oddstar lessons record-sha`.
    """
    notes = notices if notices is not None else []
    lessons_dir = root / "lessons"
    if not lessons_dir.exists():
        return []
    rc, _ = _git(root, "rev-parse", "--verify", "HEAD^{commit}")
    if rc != 0:
        notes.append("lesson_commits: no git history available; lesson commit_sha values NOT checked")
        return []
    if _git(root, "rev-parse", "--is-shallow-repository")[1] == "true":
        notes.append("lesson_commits: shallow clone (history truncated); lesson commit_sha values NOT checked. Use fetch-depth: 0")
        return []
    violations: list[str] = []
    for lf in sorted(lessons_dir.glob("*.json")):
        rec = LessonRecord.model_validate_json(lf.read_text(encoding="utf-8"))
        rc, intro = _git(root, "log", "--diff-filter=A", "--format=%H", "--", str(lf.relative_to(root)))
        intro_sha = intro.split()[-1] if intro.split() else ""
        if not rec.commit_sha:
            notes.append(
                f"lesson_commits: {lf.name} commit_sha not recorded yet"
                + (f"; introduced by {intro_sha[:10]} (run `oddstar lessons record-sha`)" if intro_sha else " (lesson not committed yet)")
            )
            continue
        rc, full = _git(root, "rev-parse", "--verify", "--quiet", f"{rec.commit_sha}^{{commit}}")
        if rc != 0:
            violations.append(f"{lf.name}: commit_sha {rec.commit_sha} does not resolve in git history")
            continue
        if _git(root, "merge-base", "--is-ancestor", full, "HEAD")[0] != 0:
            violations.append(f"{lf.name}: commit_sha {rec.commit_sha[:10]} is not an ancestor of HEAD")
        elif intro_sha and full != intro_sha:
            violations.append(
                f"{lf.name}: commit_sha {rec.commit_sha[:10]} is not the commit that introduced the lesson files ({intro_sha[:10]})"
            )
    return violations


# ------------------------------------------------------------------------------------ config version (E2)


def check_config_version(root: Path) -> list[str]:
    """configs/default.yaml `version` must equal the newest `## [x.y.z]` heading of configs/CHANGELOG.md."""
    cfg_p = root / "configs" / "default.yaml"
    cl_p = root / "configs" / "CHANGELOG.md"
    if not cfg_p.exists() or not cl_p.exists():
        return []
    version = str((yaml.safe_load(cfg_p.read_text(encoding="utf-8")) or {}).get("version", ""))
    m = re.search(r"^##\s+\[([^\]]+)\]", cl_p.read_text(encoding="utf-8"), flags=re.MULTILINE)
    if not m:
        return ["configs/CHANGELOG.md has no '## [version]' heading"]
    if m.group(1).strip() != version:
        return [f"configs/default.yaml version {version!r} differs from newest CHANGELOG entry {m.group(1).strip()!r}"]
    return []


# ---------------------------------------------------------------------------- review independence (D, R10)

INDEPENDENT_REVIEW_RE = re.compile(
    r"\b(?:independent(?:ly)?)\b(?:[\s*_`]+[\w-]+){0,3}?[\s*_`]+review(?:s|ed|ers?)?\b"
    r"|\bindependently[\s*_`]+reviewed\b",
    flags=re.IGNORECASE,
)
_NEGATION_RE = re.compile(r"\b(?:not|no|without|lacks?|lacking|never|non)[\s-]+(?:\S+\s+){0,2}$", flags=re.IGNORECASE)


def _md_units(text: str) -> list[tuple[int, str]]:
    """Paragraphs, with list items and table rows as separate units: (first line number, text)."""
    units: list[tuple[int, str]] = []
    cur: list[str] = []
    start = 1
    for i, line in enumerate(text.splitlines() + [""], start=1):
        item = re.match(r"^\s*(?:[-*+]|\d+[.)]|\|)\s", line)
        if not line.strip() or item:
            if cur:
                units.append((start, " ".join(cur)))
                cur = []
            if item:
                cur, start = [line], i
            continue
        if not cur:
            start = i
        cur.append(line)
    return units


def check_review_independence_wording(root: Path) -> list[str]:
    """Tracked Markdown under docs/ or the repo root must not call a source's review independent unless the
    round-3 record of that source comes from a different, recognised model family (`independent_family`).

    A unit (paragraph, list item or table row) is flagged when it contains "independent ... review" (up to three
    words between, not negated) and names a reviewed source ID whose `independent_family` is false.
    """
    from oddstar.review.status import SOURCE_ID_RE, review_status, reviewed_source_ids

    reviewed = set(reviewed_source_ids(root))
    if not reviewed:
        return []
    rc, out = _git(root, "ls-files", "*.md")
    if rc == 0:
        paths = [root / p for p in out.splitlines() if p.startswith("docs/") or "/" not in p]
    else:
        paths = sorted([*root.glob("*.md"), *(root / "docs").rglob("*.md")])
    cache: dict = {}
    status: dict[str, object] = {}
    violations: list[str] = []
    for p in paths:
        if not p.exists():
            continue
        for line_no, unit in _md_units(p.read_text(encoding="utf-8")):
            hits = [m for m in INDEPENDENT_REVIEW_RE.finditer(unit) if not _NEGATION_RE.search(unit[: m.start()])]
            if not hits:
                continue
            for sid in sorted(set(SOURCE_ID_RE.findall(unit)) & reviewed):
                if sid not in status:
                    status[sid] = review_status(root, sid, _cache=cache)
                st = status[sid]
                if not st.independent_family:  # type: ignore[attr-defined]
                    violations.append(
                        f"{p.relative_to(root)}:{line_no}: '{hits[0].group(0)}' for {sid}, but independent_family is false "
                        f"({st.summary}; {st.independence_basis})"  # type: ignore[attr-defined]
                    )
    return violations


def run_all_ci_checks(root: Path, notices: list[str] | None = None) -> dict[str, list[str]]:
    """Run all governance checks; `notices` (when given) collects visible skip/coverage notices."""
    results = {
        "nomenclature": check_forbidden_nomenclature(root),
        "review_isolation": check_review_kit_isolation(root),
        "lessons_findings": check_lessons_and_findings(root),
        "schemas": check_schemas_up_to_date(root),
        "r4_separation": check_git_commit_separation(root, notices),
        "config_version": check_config_version(root),
        "evidence_store": check_evidence_store(root, notices),
        "run_records": check_run_records(root, notices),
        "review_ingestion": check_review_ingestion(root),
        "lesson_commits": check_lesson_commits(root, notices),
        "review_independence": check_review_independence_wording(root),
        "certificates": check_certificates_current(root),
    }
    return results


def main() -> int:
    from oddstar.config import repo_root

    root = repo_root()
    notices: list[str] = []
    report = run_all_ci_checks(root, notices)
    for n in notices:
        print(f"NOTICE: {n}")
    failures = [f"[{k}] {v}" for k, errs in report.items() for v in errs]
    if failures:
        raise PolicyViolationError(
            "CI governance checks failed:\n" + "\n".join(failures),
            hint="Fix the reported violations and regenerate schemas with `oddstar schema --export`",
        )
    print("All CI governance checks passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
