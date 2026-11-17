# L001_R3_close: Round 3 confirmation of single-epoch W4 Tier-C cap

- **Finding ID**: `F-4582145718700356096-R3-W4GAP`
- **Target**: `4582145718700356096`
- **Automated Check / Vetting Rule**: `oddstar.fuse.engine.assign_tier`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l001_single_epoch_w4_background_galaxy_cap`
- **Added to Impostors**: `False`

## What Happened
Round 3 independent model family review confirmed W4 excess and human checklist over not_evaluated rows.

## Generic Lesson
Same generic check L001: human checklist explicitly flags deep_imaging and second_instrument.

## Automation & Regression Test
Enforced via `oddstar.fuse.engine.assign_tier` and verified in `tests/test_lessons.py::test_lesson_l001_single_epoch_w4_background_galaxy_cap`.
