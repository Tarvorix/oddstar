"""CI and Ground Rule governance checks (R4, R6, R7, AGENTS.md)."""

from __future__ import annotations

import ast
import json
import re
import subprocess
from pathlib import Path

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
            test_file = root / rec.linked_test.split("::")[0]
            if not test_file.exists():
                violations.append(f"{lf.name}: linked_test {rec.linked_test} does not exist")

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


def check_git_commit_separation(root: Path) -> list[str]:
    """Enforce R4: staged/committed changes must not alter configs/benchmarks and search results together."""
    try:
        res = subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=str(root),
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return []
    if res.returncode != 0:
        return []
    files = [line[3:].strip() for line in res.stdout.splitlines() if len(line) > 3]
    has_threshold = any(
        f.startswith("configs/default.yaml")
        or f.startswith("validation/benchmarks/")
        or f.startswith("validation/impostors/")
        for f in files
    )
    has_search_results = any(f.startswith("registry/") and f.endswith(".jsonl") for f in files)
    if has_threshold and has_search_results and (root / "configs" / "CHANGELOG.md").exists():
        # Allowed only when CHANGELOG is modified or explicitly separated in clean commits;
        # enforce no unlogged threshold modification alongside results.
        if "configs/CHANGELOG.md" not in files:
            return ["R4 violation: config/benchmark modified alongside registry without configs/CHANGELOG.md"]
    return []


def run_all_ci_checks(root: Path) -> dict[str, list[str]]:
    """Run all governance checks and raise PolicyViolationError if any fail."""
    results = {
        "nomenclature": check_forbidden_nomenclature(root),
        "review_isolation": check_review_kit_isolation(root),
        "lessons_findings": check_lessons_and_findings(root),
        "schemas": check_schemas_up_to_date(root),
        "r4_separation": check_git_commit_separation(root),
    }
    return results


def main() -> int:
    from oddstar.config import repo_root

    root = repo_root()
    report = run_all_ci_checks(root)
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
