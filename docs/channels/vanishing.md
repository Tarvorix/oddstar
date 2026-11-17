# Channel: `vanishing`

## Physical Signature
Point sources confidently detected on two or more 20th-century photographic plates/bands (USNO-B1.0) that are absent in deeper modern surveys (Pan-STARRS, Gaia) after proper-motion propagation.

## Inputs
- Historical USNO-B1.0/DASCH/DSS positions and magnitudes paired with Pan-STARRS and Gaia DR3/DR4 astrometry and local depth limits.

## Algorithm
1. Propagate modern Gaia coordinates back to historical plate epochs using `(pmra, pmdec)`.
2. Require `>= 2` independent early plate detections (`early_plates_count >= 2`) and modern survey non-detection with local depth exceeding early magnitude by `>= 2.0 mag`.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.vanishing`:
- `min_early_plates: 2`, `min_depth_margin_mag: 2.0`, `match_radius_arcsec: 3.0`, `max_pm_residual_arcsec: 1.5`

## Known Impostors
- High-proper-motion stars that moved beyond the static cross-match radius, plate emulsion flaws/scratches, large-amplitude Mira/LPV variables, asteroids, and blended sources resolved into multiple modern stars.

## Benchmark Objects
- Bidirectional switch-off test on fast-PM stars + single-epoch sources from Villarroel et al. (2020).

## Failure Modes
- Nearly all raw historical-to-modern catalog mismatches (~99.9%) are high-PM stars or single-plate emulsion flaws.

## Lessons Applied
- Mandatory multi-plate verification and proper-motion track back-propagation.
