# Channel: `spectral_lines`

## Physical Signature
Unresolved emission features narrower than or equal to the instrument line-spread function at wavelengths lacking known stellar, telluric, or sky emission lines.

## Inputs
- Multi-exposure optical spectra (`flux_exp1`, `flux_exp2`), sky-fiber spectra, and instrument LSF width from SDSS, LAMOST, GALAH, APOGEE, Gaia RVS, and APF.

## Algorithm
1. Fit non-negative linear combination of stellar continuum/absorption templates (`nnls`) on the per-pixel minimum across exposures and subtract. The templates are a constant, a slope of either sign, two H-alpha profiles and (config 1.12.0) 21 of the stronger photospheric lines of 6493 to 6664 A with free depths (`PHOTOSPHERIC_LINES_A`).
2. Run LSF-width Gaussian matched filter across residual spectra, requiring multi-exposure persistence.
3. Veto known night-sky/airglow lines and chromospheric stellar lines (e.g., H-alpha).
4. Config 1.12.0: the significance of an exposure at a pixel is `(filter output - baseline) / max(per-pixel noise, scatter)`. Baseline and scatter are the median and half the 16th-to-84th percentile range of the filter output over the open pixels, taken without the pixels within 2 line-spread FWHM of the strongest peak. The channel statistic is the smaller of the two exposures' values at the best pixel. Each measurement stores the four numbers per exposure in `extra.filter_calibration`.

## Live survey search (config 1.3.0)
`oddstar search spectral_lines --live` loads two epochs of every sample star from the LAMOST DR12 v2.0 low-resolution
survey before measuring (`oddstar/io/lamost_bulk.py`):

1. The general catalogue is asked, in batches of 50 Gaia source ids, which spectra exist. Every batch also asks for a
   star known to be in the catalogue; an answer without it is not trusted, is not cached, and the batch counts as failed.
2. Epochs are tried from the highest red signal-to-noise down, one per night, skipping epochs below `lamost_min_snr_r`,
   until two usable ones are loaded. A spectrum is reduced on arrival to the window 6500 to 6650 A: vacuum to air, moved
   to the stellar rest frame with the pipeline redshift, masked pixels interpolated, divided by a straight-line continuum
   fitted away from H-alpha. The cache keeps the reduced window with the URL, size and SHA-256 of the file.
3. The two epochs become `flux_exp1` and `flux_exp2` on a 1.5 A grid with `lsf_fwhm_a = 3.65` (R = 1800).
   `sky_line_wavelengths_a` lists where the night-sky lines (the vetting list plus the other OH 6-1 lines of the window)
   fall in the star's rest frame for each epoch. For a halo star at 300 km/s that is 6.6 A from the laboratory position.

At this resolution the channel differs from the high-resolution case in three ways, all driven by `lsf_fwhm_a`: the
H-alpha templates are broadened by the line spread, a known line is excluded out to 1.5 line-spread FWHM (5.5 A), and
every listed sky position is excluded, geocoronal H-alpha included. About half of the window stays open to the search;
each measurement records the fraction (`extra.searched_window_fraction`) and the observed wavelength of its peak in both
epochs (`extra.peak_observed_wavelengths_a`).

A star is **searched** only when two epochs from different nights were loaded. Every other star is listed with its
reason and excluded from the null result:

| Reason | Meaning |
|---|---|
| `single_epoch` | spectra exist, but fewer than two usable epochs on different nights |
| `no_coverage` | no spectrum of this Gaia source in the catalogue (LAMOST observes north of about Dec -10) |
| `failed` | the catalogue batch or a needed file failed after `lamost_max_attempts`; not cached, retried on the next run |
| `not_fetched` | beyond `--spectra-max-stars` |
| `no_spectrum` | no fetch record and no usable spectra |

With `--defer-enrichment` the flagged stars are still written to `live_stars.json` with their spectra, because
`vet --live` cannot fetch the spectra again for one star and the spectral checks read nothing else.

What a live flag means, and what it does not:
- The certificate of this channel was measured on fixtures at a line-spread width of 0.45 A. Nothing was injected at
  R = 1800, so the completeness and the false-alarm rate of a LAMOST search are not measured. A flag is a lead.
- Up to config 1.11.0 the template was a continuum and H-alpha only. The pseudo-continuum between the unmodelled
  absorption lines is then a repeatable bump of 1.6 to 3.3 per cent in nearly every star (1,535 LAMOST spectra,
  strongest at rest 6504 to 6512 A), and it reached the threshold in spectra of S/N above 170: the one star
  flagged in the 1,300 nearby stars (29259241366225024, S/N 6.00 at 6507.5 A) was this bump. Config 1.12.0 fits
  the stronger metal lines and judges the filter output against its own baseline and scatter (see Algorithm and
  the changelog). What is left of a star's spectrum after the fit still limits the search: in a spectrum of high
  signal-to-noise the smallest line that flags is set by the scatter of the filter output, not by photon noise.
- No sky spectrum and no detector map come with the public file: `sky_line` returns `not_evaluated` unless the peak
  sits on a listed sky position (then `fail`), and `detector` is always `not_evaluated`. A sky or instrument feature
  that is not in the list shows as flagged peaks at the same observed wavelength in several stars.
- The two exposures of the fixtures are the same night; here they are different nights, so a line that varies between
  nights can fail `repeat` or `cosmic_ray`.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.spectral_lines`:
- `residual_snr_threshold: 6.0`, `min_exposures_repeat: 2`, `sky_line_exclusion_angstrom: 1.5`, `stellar_line_exclusion_angstrom: 1.5`
- Live fetch: `lamost_min_snr_r: 20.0` (an epoch below this red signal-to-noise is not used), and the operational keys
  `lamost_max_attempts: 3`, `lamost_backoff_s: 2.0`, `lamost_rate_limit_rps: 2.0`
- The two exclusion widths are minimum values: a known line is excluded out to 1.5 x `lsf_fwhm_a` when that is larger

## Known Impostors
- Cosmic-ray hits, OH airglow / geocoronal emission lines (Lipman et al. 2019), chromospheric M-dwarf H-alpha emission, and detector cross-talk/bad pixels.

## Benchmark Objects
- Bidirectional switch-off test: disabling stellar-line veto flags active M dwarfs; enabling it explains them.

## Failure Modes
- Single-exposure cosmic rays mimic unresolved laser lines unless `repeat` and sub-LSF width checks run.

## Lessons Applied
- Require multi-exposure persistence and joint sky-spectrum matching.
