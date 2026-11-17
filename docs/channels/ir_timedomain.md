# Channel: `ir_timedomain`

## Physical Signature
Year-scale mid-infrared brightening, fading, steps, or bursts (3.4–4.6 um) while optical photometry stays steady or dims years later (collision-debris cascade signature).

## Inputs
- NEOWISE-R single-exposure `W1` and `W2` photometry binned into ~180-day visits (2013–2024).
- Overlapping ZTF `g`/`r`, ASAS-SN, and TESS optical light curves.

## Algorithm
1. Group NEOWISE exposures into semi-annual visit medians with MAD uncertainties.
2. Compute reduced `chi2_w1`, `chi2_w2`, linear slope SNR, largest single-visit jump SNR, and `W1 - W2` color excursion.
3. Require optical curve quietness (`optical_chi2 <= optical_flat_max_chi2`) during the baseline or early pre-eclipse window.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.ir_timedomain`:
- `reduced_chi2_threshold: 5.0`, `slope_snr_threshold: 4.0`, `jump_snr_threshold: 4.5`
- `optical_flat_max_chi2: 2.5`, `min_visits: 6`

## Known Impostors
- Bright blended W1/W2 neighbours, Mira/LPV long-period variables, and scan-direction detector systematics.

## Benchmark Objects
- Variable disks in Moor et al. (2021), `TYC 4209-1322-1`, and `ASASSN-21qj`.

## Failure Modes
- Single-exposure NEOWISE queries are heavy; search runs preferentially on `ir_excess`-flagged targets plus controls.

## Lessons Applied
- Enforce `neighbour_signal` and `systematics_window` checks against bright W1/W2 blends and scan-angle systematics.

## Light-curve vetting rules for IR variability (config 1.1.0)
Fixtures set explicit `w1_w2_neighbour_blend` / `scan_systematic` flags (present = evaluated). Live stars are vetted on data:

- **`systematics_window`** (NEOWISE part): enrichment fetches `neowiser_p1bs_psd` within `neowise_comparison_radius_arcsec`, links
  detections into stars (1.5" friends-of-friends), drops the target (`neowise_comparison_exclude_arcsec`), keeps up to
  `neowise_comparison_max_stars` stars within `neowise_comparison_mag_window` W1 mag covering >= half the target visits, and forms
  the per-visit ensemble median W1/W2 offset on the target's visit grid (`neowise.comparison`, `None` where fewer than
  `neowise_comparison_min_stars` stars). When the target W2 variation passes this channel's own thresholds, the check fails as a
  scan/zero-point systematic if the target correlates with the ensemble (r >= `systematics_ensemble_corr` with ensemble amplitude
  >= `systematics_ensemble_corr_min_amp_frac` x target) or the ensemble varies significantly (reduced chi2 >= `reduced_chi2_threshold`)
  by >= `systematics_ensemble_amp_frac` x the target amplitude. Fewer than `neowise_comparison_min_stars` stars or a failed query
  gives `not_evaluated` for this part.
- **`neighbour_signal`** (NEOWISE part): per-visit RA*cos(Dec)/Dec residuals from the Gaia PM track (`neowise.pm_residual_*_arcsec`)
  are fitted with offset + scan-direction parity + linear trend + W2 fractional flux. (a) A stationary blend makes the residual lag
  the PM track: the along-track trend gives a stationary flux share, failing when >= `vetting.neighbour_blend_flux_frac` at
  >= `neighbour_blend_sigma` (tested when the track is >= `neighbour_min_track_arcsec`). (b) A variable stationary source makes the
  residual follow the W2 flux: fail when that dependence is >= `neighbour_blend_sigma`. The note also gives the pull a stationary
  source at the star's 2010.5 position carrying the variable flux would produce. For Gaia DR3 5628404318169472 (review evidence of
  2026-11-17) the median residual is 0.053", the stationary share 0.03 +/- 0.01 and the predicted pull 0.35" against a measured
  0.009 +/- 0.017": the variable flux moves with the star (pass).
- **`second_survey`** (IR context): per survey/band seasons with >= `second_survey_min_season_nights` nights carry a robust error of
  the median (from nightly medians) plus `second_survey_season_floor_mag` in quadrature. The max-min seasonal range is significant at
  >= `second_survey_change_sigma`. W2 is interpolated to the season epochs; with >= `second_survey_min_overlap_seasons` overlapping
  seasons a weighted fit gives dOpt/dW2 and Pearson r. A band is `dimming_with_ir_brightening` (negative slope) or
  `brightening_with_ir_brightening` (positive slope) when |r| >= `second_survey_min_abs_corr` and |slope| >=
  `second_survey_slope_sigma` sigma, `uncorrelated_variation` when significant otherwise, `none` when not significant. The IR event
  window is the span of visits deviating from the first-three-visit baseline by >= half the maximum deviation; the note states
  which part each survey covers (ZTF from 2018, ASAS-SN V to ~2019). `fail` when a band brightens with W2 at a slope >=
  `second_survey_comparable_frac` (optical change comparable to the IR one: photospheric/accretion variability or a blended
  variable, contradicting circumstellar dust); `not_evaluated` when no season overlaps the window; otherwise `pass` with the
  classification and numbers, never described as clean when a significant change was found.
