# L005_channel_scoped_vetting_and_southern_decam: Scope metallicity away from IR channels, route all WISE outliers, enforce L1b/cc_flags/w4nm quality, and enable southern DECam ls_dr10

- **Finding ID**: `F-6739572937544816000-ADV-L005`
- **Target**: `6739572937544816000`
- **Automated Check / Vetting Rule**: `oddstar.fuse.engine.assign_tier + oddstar.vet.engine.vet_target + oddstar.vet.ir_checks.check_wise_quality + oddstar.io.live_survey.enrich_star_live_for_vetting`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l005_channel_scoped_metallicity_wise_routing_and_southern_ls_dr10`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `9eb6dc4a97d7319d4241d471571294e69bb4537e`

## What Happened
Adversarial audit of 7,800-star survey exposed four pipeline leaks: (1) metallicity=fail masked ir_excess on Halo Batch 4 star 6739572937544816000 in fusion; (2) generic_outlier targets driven by ks_w3/w3_w4 (LMC star 4651844963414920576) bypassed IR vetting; (3) L1b single-exposure upper-limit frames without w4sigmpro inflated single_frames consistency and low w4nm/w4m or hh00 W1/W2 halos passed wise_quality (1598166882848977024, 6739572937544816000); (4) ls_dr10 deep imaging was gated to Dec > -30 deg and sent JSON format to DataLab.

## Generic Lesson
Scope metallicity explanations strictly away from ir_excess/ir_timedomain, route all WISE-driven generic_outlier leads through run_ir_vetting, require detected L1b sigmas/w4nm fraction and clean 4-band cc_flags, and query DECam ls_dr10.tractor via CSV across all southern sky fields.

## Automation & Regression Test
Enforced via `oddstar.fuse.engine.assign_tier + oddstar.vet.engine.vet_target + oddstar.vet.ir_checks.check_wise_quality + oddstar.io.live_survey.enrich_star_live_for_vetting` and verified in `tests/test_lessons.py::test_lesson_l005_channel_scoped_metallicity_wise_routing_and_southern_ls_dr10`.
