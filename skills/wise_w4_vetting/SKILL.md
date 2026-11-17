---
name: wise_w4_vetting
description: Adversarial vetting rules for single-epoch WISE W4 (22 um) infrared-excess candidates derived from Lesson L001 and Ren et al. (2024, 2026).
---

# WISE W4 Single-Epoch & Background-Galaxy Vetting Skill

1. **Never treat clean radio vetting as proof against a background galaxy**: Fewer than 10% of dusty star-forming galaxies inside a 12" WISE W4 beam are radio-detected (`Ren et al. 2026`).
2. **Enforce `not_evaluated` hole accounting**: If deep sub-arcsec optical/NIR imaging (`HST`, `JWST`, or deep `Legacy/HSC/AO`) or post-2010 overlapping 22 um photometry (`second_instrument`) is absent, `deep_imaging` and `second_instrument` MUST remain `not_evaluated` and the candidate MUST stay capped at Tier `C`.
3. **Always emit the single-epoch caveat**: Cryogenic W3/W4 ended in 2010; no public archive confirms post-2010 22 um stability.
