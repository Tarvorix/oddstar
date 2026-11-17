# L002_missing_plate_data_high_pm_vanishing: Missing historical plate record must not default to J2016 coordinates in vanishing channel

- **Finding ID**: `F-5628404318169472-LIVE-VANISH`
- **Target**: `5628404318169472`
- **Automated Check / Vetting Rule**: `oddstar.channels.vanishing.channel.VanishingChannel.measure`
- **Linked Test**: `tests/test_lessons.py::test_lesson_l002_missing_plate_data_no_spurious_vanishing_on_high_pm_stars`
- **Added to Impostors**: `False`
- **Commit (introduced the lesson files)**: `c678696db380cc1995290c47ea90f547bb8d502e`

## What Happened
During live single-star inspection of high-proper-motion halo dwarf Gaia DR3 5628404318169472 (v_tan=114.6 km/s) with empty vanishing_data, VanishingChannel.measure defaulted usno_ra_deg/usno_dec_deg to J2016 Gaia coordinates and early_plates_count to 2 while back-propagating proper motion to 1953.4, producing a synthetic 6.53 arcsec mismatch and false vanishing trigger.

## Generic Lesson
Channels comparing multi-epoch catalogs must never substitute modern coordinates for absent historical epoch measurements or assume non-zero plate counts when archival plate data are missing.

## Automation & Regression Test
Enforced via `oddstar.channels.vanishing.channel.VanishingChannel.measure` and verified in `tests/test_lessons.py::test_lesson_l002_missing_plate_data_no_spurious_vanishing_on_high_pm_stars`.
