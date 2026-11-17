# L006_matched_detrending_and_robust_optical_surveys: Fetch TESS QLP/TESScut, ZTF and ASAS-SN live for LC vetting, with matched detrending and retried/paged ZTF queries

- **Finding ID**: `F-5628404318169472-20261117-R2-OPTICAL`
- **Finding source**: review/records/REV-20261117-R2-5628404318169472.json
- **Target**: `5628404318169472`
- **Automated Check / Vetting Rule**: `oddstar.vet.lc_checks.check_second_survey + oddstar.vet.lc_checks.check_pipeline_agreement + oddstar.vet.lc_checks.matched_detrend + oddstar.io.optical_timeseries.fetch_ztf_lightcurve + oddstar.io.optical_timeseries.fetch_asassn_lightcurve + oddstar.io.optical_timeseries.fetch_tess_pipeline_pair`
- **Linked Test**: `tests/test_lc_vetting_live.py::test_second_survey_reference_star_dims_while_ir_brightens`
- **Additional Tests**: `tests/test_optical_vetting.py::test_pipeline_agreement_matched_detrending`, `tests/test_lc_vetting_live.py::test_second_survey_counter_case_and_coverage`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `2b24823d0b8909d526951058b1d972a91c98d3a8`

## What Happened
The sole Tier C survivor 5628404318169472 carried two not_evaluated holes that code could close. pipeline_agreement never received a second TESS extraction because live enrichment fetched no TESS data, and the adversarial review's hand comparison of QLP vs. TESScut (raw r = 0.43-0.81) used different detrending on each series. second_survey never received ZTF or ASAS-SN because live enrichment had no optical survey fetch, and the reviewer's own IRSA ZTF cone requests timed out at 40 s unless narrowed to one band list and a 6 arcsec cone.

## Generic Lesson
A cross-pipeline comparison is only meaningful when both series go through the identical detrending filter at matched cadences, using the detection channel's own window. Archive time-series queries that can time out need bounded retries, a shrinking cone and time-window paging, and every failure mode must surface in the check note instead of silently reading as no coverage. Vetting-only time series are stored apart from channel inputs so new vetting data never changes a certified channel score.

## Method
- **TESS**: TIC resolved by MAST TIC position search (Gaia ID match preferred). For every TESScut sector with a QLP light curve: QLP `SAP_FLUX` (undetrended, `QUALITY == 0`) and an independent TESScut aperture extraction (1.6 px circular aperture at the PM-propagated position, 4-6 px annulus median background, `QUALITY == 0`) are matched cadence-by-cadence and binned to 30 min. `pipeline_agreement` then divides both by the same running median (`optical_dips.detrend_window_days`, per continuous segment, with low outliers masked), bins to `optical_dips.min_duration_hours`, and applies the existing depth-agreement rule per sector. Any sector disagreeing fails the check.
- **ZTF**: IRSA light-curve API per band (`g`, `r`); 3 attempts with exponential backoff at a 3" cone, then yearly `TIME` pages, then a 2" cone. Rows need `catflags == 0` and an epoch position within 1.5" of the Gaia position propagated to that epoch. Seasonal (calendar-year) nightly-median ranges are reported for the dominant object ID per band.
- **ASAS-SN**: Sky Patrol v2 `stellar_main` cone (15") through the public `pyasassn` client (pip package `skypatrol`, optional extra `oddstar[surveys]`), nearest source within 5", `quality == 'G'`, `0 < mag_err < 0.05`, per-camera median removed, July-June seasons.
- No configuration thresholds changed: the detrending window and bin width reuse existing `optical_dips` settings, and the depth-agreement rule is unchanged.

## Result on 5628404318169472 (live, 2026-11-17)
- `pipeline_agreement`: **pass**. TIC 387517877, sectors 4, 31, 42, 43, 44, 70, 71. Max 6 h dip QLP 0.028 vs. TESScut 0.027; per-sector binned Pearson r = 0.76-0.99.
- `second_survey`: **pass**. ZTF full-range queries succeeded on the first try (1,196 clean rows). Yearly nightly-median magnitudes match the adversarial review's hand-fetched ZTF table within 0.003 mag for every year (2018-2024 range 0.060 g / 0.053 r vs. review 0.060 / 0.055). ASAS-SN camera-subtracted seasonal medians match the review's within 0.0013 mag in every V and g season (seasonal range 0.062 g / 0.048 V).
- Tier stays **C**; open holes reduce from four to two (`second_instrument`, `binary`).

## Automation & Regression Test
Enforced via `oddstar.vet.lc_checks.check_pipeline_agreement + oddstar.vet.lc_checks.matched_detrend + oddstar.vet.lc_checks.check_second_survey + oddstar.io.optical_timeseries.fetch_ztf_lightcurve + oddstar.io.optical_timeseries.fetch_asassn_lightcurve + oddstar.io.optical_timeseries.fetch_tess_pipeline_pair` and verified in `tests/test_optical_vetting.py`.

## Finding linkage (amended 2026-11-17)
- Finding: `F-5628404318169472-20261117-R2-OPTICAL` (round-2 record `review/records/REV-20261117-R2-5628404318169472.json`): ASAS-SN finds small optical seasonal changes across the infrared event. The earlier ID `F-5628404318169472-ADV-L006` matched no review record.
- What closes it: `check_second_survey` now measures the optical change instead of passing on its absence. On the reference star it reports ASAS-SN g range 0.060 mag at 6.9 sigma, response `dimming_with_ir_brightening`, note "significant optical change detected" (test `tests/test_lc_vetting_live.py::test_second_survey_reference_star_dims_while_ir_brightens`). That response classification was added in commit 0e5956a, after this lesson; the second_survey **pass** reported above under "Result" predates it and by itself did not address the finding.
- Not addressed: `F-5628404318169472-20261117-R2-PRIOR-IR` (VarWISE prior catalogue entry) and `F-5628404318169472-20261117-R3-UNRESOLVED-BLEND` stay open. The matched detrending of `pipeline_agreement` partly covers unticked human checklist item 4 of `review/adversarial-2026-11-17.md`; that item stays for a human.
