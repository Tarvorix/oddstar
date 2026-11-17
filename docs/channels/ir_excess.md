# Channel: `ir_excess`

## Physical Signature
Mid-infrared excess emission above a star's photospheric spectrum at 12 um (WISE W3) and 22 um (WISE W4), probing circumstellar warm dust (~100–700 K) or re-radiated waste heat at terrestrial-to-asteroid-belt separations.

## Inputs
- Gaia DR3/DR4 astrometry (`ra`, `dec`, `parallax`, `ruwe`) and `G`, `BP`, `RP` photometry.
- 2MASS `J`, `H`, `Ks` photometry.
- AllWISE `W1`–`W4` photometry, profile-fit `chi2`, contamination/extension flags, and 2010 cryogenic single-exposure frames.

## Algorithm
1. Apply parent-sample cuts (`parallax_snr > 10`, `distance_pc <= 300`, `ruwe < 1.4`, `|b| > 10 deg`, unique AllWISE cross-match, clean 2MASS).
2. Predict photospheric W3 and W4 fluxes independently via (a) weighted Planck SED fit to `G, BP, RP, J, H, Ks` and (b) empirical `(Ks - W3)` and `(Ks - W4)` vs `(BP - RP)` median/MAD locus.
3. Compute excess statistics `chi_w3`, `chi_w4` and flux ratios under both predictors.
4. Fit single-temperature dust blackbody `(T_dust, f_dust, r_bb)` and Suazo et al. (2024) partial-shell covering-fraction grid.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.ir_excess`:
- `chi_w3_threshold: 5.0`, `chi_w4_threshold: 5.0`
- `min_flux_ratio_w4: 2.0`, `min_flux_ratio_w3: 1.5`
- `cal_uncertainty_frac: 0.05`, `require_both_predictors: true`

## Known Impostors
- Radio-loud and radio-quiet dusty background galaxies inside the 12" WISE W4 beam (Ren et al. 2024, 2026).
- Dust-obscured galaxies (DOGs) and red neighbours.
- Galactic dust cirrus filaments, AllWISE diffraction spikes (`cc_flags`), and single-frame glitches.

## Benchmark Objects
- `BD+20 307`, `ID8`, `HD 23514`, `TYC 4209-1322-1`, `HD 113766`, `HD 172555`, `HD 69830`, and the six extreme debris disks in Moor et al. (2021).
- Open case workflow target: `Gaia DR3 4582145718700356096`.

## Failure Modes
- W3 and W4 photometry exist solely from the 2010 cryogenic WISE mission (single epoch, single instrument).
- Over 90% of background galaxies contaminating WISE beams are radio-quiet, so passing `radio_counterpart` never clears `background_galaxy` without deep sub-arcsec optical/NIR imaging.

## Lessons Applied
- `L001`: Single-epoch 22 um W4 excesses lacking deep sub-arcsec imaging and second-instrument 22 um follow-up remain capped at Tier C with `background_galaxy` and `deep_imaging` flagged explicitly as open holes.
