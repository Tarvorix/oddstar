# L003_adversarial_audit_and_missing_modality_routing: Prevent synthetic noise leakage, W4 upper-limit leakage, and optimistic Gaia RVS binary clearance

- **Finding ID**: `F-1111149103659931264-ADV-L003`
- **Related findings (addressed, not closed)**: `F-1111149103659931264-20261117-R2-KNOWN-SUBDWARF`
- **Finding source**: pre-review adversarial audit of 2026-11-17 (no review record; the review records start after it). Item (1) concerns 1111149103659931264, item (4) the RV binary test on 5628404318169472.
- **Target**: `1111149103659931264`
- **Automated Check / Vetting Rule**: `oddstar.channels.generic_outlier.channel.extract_raw_feature_vector + oddstar.vet.engine.vet_target + oddstar.vet.ir_checks.check_binary`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l003_missing_modality_routing_and_rv_binary_honesty`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `315e9b017742ff25318a9399dc31e7d14eeec303`

## What Happened
An adversarial audit identified four pipeline weaknesses: (1) generic_outlier injected random Gaussian noise when lightcurves/spectra/vanishing_data were absent, producing phantom open holes and retaining metal-poor halo subdwarf Gaia DR3 1111149103659931264 at Tier C; (2) live sample ingest allowed W4='U' upper limits (w4mpro_error IS NULL) to flood ir_excess flags; (3) _row_to_star_dict compared unpropagated 2016 Gaia vs 2010 WISE coords (0.58 arcsec drift); and (4) Gaia DR3 radial_velocity_error (median uncertainty) was treated as epoch scatter on faint rv_method_used=2 stars without rv_chisq_pvalue.

## Generic Lesson
Use deterministic 0.0 for absent modalities in generic_outlier, run metallicity/evolved checks on all generic_outlier stars, require w4mpro_error IS NOT NULL in live excess searches, propagate proper motion before centroid offsets, and convert median RV uncertainty via sqrt(2N/pi) (marking faint unconstrained Gaia RVS stars as binary=not_evaluated).

## Automation & Regression Test
Enforced via `oddstar.channels.generic_outlier.channel.extract_raw_feature_vector + oddstar.vet.engine.vet_target + oddstar.vet.ir_checks.check_binary` and verified in `tests/test_lessons.py::test_lesson_l003_missing_modality_routing_and_rv_binary_honesty`.
