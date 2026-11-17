# L001_single_epoch_w4_background_galaxy_cap: Single-epoch W4 excess and radio-quiet background galaxy ambiguity

- **Finding ID**: `F-4582145718700356096-R1-W4GAP`
- **Target**: `4582145718700356096`
- **Automated Check / Vetting Rule**: `oddstar.fuse.engine.build_hypothesis_table + assign_tier`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l001_single_epoch_w4_background_galaxy_cap`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `f7ec8a65eed728799f6a046395e46b0a1c23a1bb`

## What Happened
Three independent review rounds confirmed the ~56-60x 22um WISE excess for Gaia DR3 4582145718700356096 from raw archive photometry, while identifying that radio_counterpart passing does not exclude the >90% radio-quiet dusty galaxy population without deep sub-arcsec optical/NIR imaging.

## Generic Lesson
Passing radio_counterpart never moves background_galaxy out of untested when deep_imaging is not_evaluated; single-epoch W4 candidates missing deep imaging or second-instrument 22um data must remain capped at Tier C.

## Automation & Regression Test
Enforced via `oddstar.fuse.engine.build_hypothesis_table + assign_tier` and verified in `tests/test_lessons.py::test_lesson_l001_single_epoch_w4_background_galaxy_cap`.
