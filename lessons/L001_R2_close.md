# L001_R2_close: Round 2 confirmation of single-epoch W4 Tier-C cap

- **Finding ID**: `F-4582145718700356096-R2-W4GAP`
- **Target**: `4582145718700356096`
- **Automated Check / Vetting Rule**: `oddstar.fuse.engine.assign_tier`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l001_single_epoch_w4_background_galaxy_cap`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `f7ec8a65eed728799f6a046395e46b0a1c23a1bb`

## What Happened
Round 2 adversarial review confirmed W4 excess and unresolved deep_imaging hole.

## Generic Lesson
Same generic check L001: Tier capped at C while deep_imaging or second_instrument is not_evaluated.

## Automation & Regression Test
Enforced via `oddstar.fuse.engine.assign_tier` and verified in `tests/test_lessons.py::test_lesson_l001_single_epoch_w4_background_galaxy_cap`.
