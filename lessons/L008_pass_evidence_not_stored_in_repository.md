# L008_pass_evidence_not_stored_in_repository: Evidence behind vetting passes was not stored in the repository

- **Finding ID**: `F-REPO-20261117-EVIDENCE-NOT-STORED`
- **Finding source**: repository audit 2026-11-17 (no review record)
- **Target**: `5628404318169472`
- **Automated Check / Vetting Rule**: `oddstar.evidence.export_candidate_evidence + oddstar.evidence.verify_evidence + oddstar.ci_checks.check_evidence_store`
- **Linked Test**: `tests/test_evidence_store.py::test_verification_fails_on_tampered_or_missing_evidence`
- **Additional Tests**: `tests/test_evidence_store.py::test_fuse_registers_candidate_with_evidence_reference`, `tests/test_evidence_store.py::test_ci_check_requires_vetting_rows_and_files`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `not recorded yet (filled by `oddstar lessons record-sha` after the commit that adds this lesson)`

## What Happened
Registry ledger entries recorded tier and check status, but the per-check notes, metrics, consumed manifest IDs / SHA-256 values / query URLs and images that supported each pass lived only under the git-ignored runs/ and dossiers/ directories. A pass row in the ledger could therefore not be checked from a clone of the repository. The registrations made before the evidence store (including 5628404318169472) remain without committed evidence and are reported as legacy_without_evidence.

## Generic Lesson
Every registration commits a compact, size-guarded evidence directory (registry/evidence/<source_id>/<run_id>/) with an index of file SHA-256 values, and the ledger entry references the directory and the index SHA-256. CI verifies the hashes and that every pass/fail row has its evidence file; entries without evidence are counted, not hidden.

## Automation & Regression Test
Enforced via `oddstar.evidence.export_candidate_evidence + oddstar.evidence.verify_evidence + oddstar.ci_checks.check_evidence_store` and verified in `tests/test_evidence_store.py::test_verification_fails_on_tampered_or_missing_evidence`.
