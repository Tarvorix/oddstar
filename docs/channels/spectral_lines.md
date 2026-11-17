# Channel: `spectral_lines`

## Physical Signature
Unresolved emission features narrower than or equal to the instrument line-spread function at wavelengths lacking known stellar, telluric, or sky emission lines.

## Inputs
- Multi-exposure optical spectra (`flux_exp1`, `flux_exp2`), sky-fiber spectra, and instrument LSF width from SDSS, LAMOST, GALAH, APOGEE, Gaia RVS, and APF.

## Algorithm
1. Fit non-negative linear combination of stellar continuum/absorption templates (`nnls`) on the per-pixel minimum across exposures and subtract.
2. Run LSF-width Gaussian matched filter across residual spectra, requiring multi-exposure persistence.
3. Veto known night-sky/airglow lines and chromospheric stellar lines (e.g., H-alpha).

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.spectral_lines`:
- `residual_snr_threshold: 6.0`, `min_exposures_repeat: 2`, `sky_line_exclusion_angstrom: 1.5`, `stellar_line_exclusion_angstrom: 1.5`

## Known Impostors
- Cosmic-ray hits, OH airglow / geocoronal emission lines (Lipman et al. 2019), chromospheric M-dwarf H-alpha emission, and detector cross-talk/bad pixels.

## Benchmark Objects
- Bidirectional switch-off test: disabling stellar-line veto flags active M dwarfs; enabling it explains them.

## Failure Modes
- Single-exposure cosmic rays mimic unresolved laser lines unless `repeat` and sub-LSF width checks run.

## Lessons Applied
- Require multi-exposure persistence and joint sky-spectrum matching.
