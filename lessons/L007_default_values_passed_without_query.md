# L007_default_values_passed_without_query: Vetting checks passed on default input values when no archive query had run

- **Finding ID**: `F-REPO-20261117-DEFAULT-PASS`
- **Finding source**: repository audit 2026-11-17 (no review record)
- **Target**: `5628404318169472`
- **Automated Check / Vetting Rule**: `vetting_archive.query_status rule in oddstar.vet.ir_checks and oddstar.vet.lc_checks (missing or failed query -> not_evaluated)`
- **Linked Test**: `tests/test_ir_vetting_live.py::test_no_query_is_not_evaluated`
- **Additional Tests**: `tests/test_lc_vetting_live.py::test_fixture_flags_still_decide_and_live_absence_is_not_pass`, `tests/test_ir_vetting_live.py::test_flat_inputs_ignored_when_query_failed`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `ba94372ce1d6b547b89dd9da9066ab25ddcdfe57`

## What Happened
Before commits 7d3702f (infrared) and 0e5956a (light curve), live enrichment filled vetting inputs with benign defaults when no archive query supplied them (for example n_neighbours_gaia_xmatch=1, cirrus_bg_rms_mjy=0.4, galex_uv_excess_mag=0.0, xray_luminosity_erg_s=0.0, rv_scatter_kms=0.3, nss_flag=0), and several light-curve checks read fixture-only flags. Checks then reported pass on those defaults, so a pass row did not show that the archive had been queried. This applied to the live candidates, including 5628404318169472.

## Generic Lesson
A vetting check may report pass or fail only from data that a recorded query returned. Live enrichment records the outcome of every query in vetting_archive.query_status; when the query did not run or failed, the check reports not_evaluated and names the missing query, and a stale input value is not used.

## Automation & Regression Test
Enforced via `vetting_archive.query_status rule in oddstar.vet.ir_checks and oddstar.vet.lc_checks (missing or failed query -> not_evaluated)` and verified in `tests/test_ir_vetting_live.py::test_no_query_is_not_evaluated`.
