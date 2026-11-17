# L004_known_pulsating_variable_simbad_types: Classify RR* and classical pulsating SIMBAD types in check_known_variable

- **Finding ID**: `F-6164530246002222208-HALO-B2`
- **Target**: `6164530246002222208`
- **Automated Check / Vetting Rule**: `oddstar.vet.lc_checks.check_known_variable`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l004_known_pulsating_variable_simbad_types`
- **Added to Impostors**: `False`

## What Happened
In All-Sky Galactic Halo Batch 2 (live-allsky-halo-b2-out), horizontal-branch RR Lyrae variable Gaia DR3 6164530246002222208 (SIMBAD RR*, M_G=0.50, Teff=6439 K, [Fe/H]=-0.88, RV scatter=53.6 km/s) landed at Tier B because check_known_variable only checked {EB*, LPV*, Mira, CV*}.

## Generic Lesson
Include RR*, LP*, Cepheid, delta Scuti, SR*, and eclipsing/cataclysmic SIMBAD codes in check_known_variable so known pulsating variables are resolved as Explained.

## Automation & Regression Test
Enforced via `oddstar.vet.lc_checks.check_known_variable` and verified in `tests/test_lessons.py::test_lesson_l004_known_pulsating_variable_simbad_types`.
