# Channel: `generic_outlier`

## Physical Signature
Hypothesis-free discovery of targets occupying sparse/isolated regions of joint multi-channel photometric, astrometric, spectral, and time-domain feature space.

## Inputs
- Joint feature vector combining colors (`BP-RP`, `J-Ks`, `Ks-W3`, `W3-W4`), channel statistics (`chi_w4`, `ir_var_chi2`, `dip_snr`, `delta_m_sigma`, `spec_line_snr`, `vanished_margin`, `radio_snr`), and `ruwe`.

## Algorithm
1. Conditionally standardize the always-measured features (`BP-RP`, `J-Ks`, `Ks-W3`, `W3-W4`, `chi_w4`, `delta_m_sigma`, `ruwe`) by median and MAD within `BP-RP` color bins.
2. Express each *modality* feature (`ir_var_chi2`, `dip_snr`, `spec_line_snr`, `vanished_margin`, `radio_snr`) as a significance against its own noise model, in units of the noise sigma, when the observation exists, and as exactly `0.0` when it does not (no evidence either way, L003):
   - `dip_snr`, `spec_line_snr` (most extreme of n points in robust sigma): minus the expected maximum of n Gaussian points (Blom approximation; n from the data, else `modality_null_default_n_points`);
   - `ir_var_chi2` (reduced chi^2 of the NEOWISE visits about their median): rms excess `sqrt(chi2_r) - 1`;
   - `vanished_margin` (mag below the modern survey limit): divided by `vanished_margin_plate_sigma_mag`;
   - `radio_snr`: already a signal-to-noise ratio about a zero null.
   An ordinary light curve, spectrum or NEOWISE series therefore sits near 0 instead of being an outlier merely because the background stars lack that observation.
3. Score with seeded Isolation Forest (`n_trees >= 200`, `subsample_size = 256`) and `cKDTree`-backed Local Outlier Factor (`k = 20`). As in Liu et al. (2008) each split attribute is drawn from all features; an attribute that is constant in a node's data gives a pass-through node that still isolates a query lying more than 1.5 standardized units outside the data.
4. Flag targets falling in the top fraction (`>= 0.99` quantile) of both detectors simultaneously with `max |z| >= 3`, and output per-feature standardized contributions (`present_<feature>` records which modalities were observed).

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.generic_outlier`:
- `iforest_n_trees: 100`, `lof_k_neighbors: 20`, `min_combined_quantile: 0.99`
- `modality_null_default_n_points: 1000`, `vanished_margin_plate_sigma_mag: 0.3` (since 1.1.0)

## Known Impostors
- Combined impostors from every channel; flagged stars automatically inherit vetting checks tied to driving features.

## Benchmark Objects
- Positives from the other seven channels (`BD+20 307`, `TYC 4209-1322-1`, `KIC 8462852`, `ASASSN-21qj`, `EPIC 204278916`, `VAN-Villarroel-1`).

## Failure Modes
- Unconditioned colors bias cool/hot stars into outlier tails; conditional `BP-RP` standardization eliminates this artifact.
- Modality presence is not evidence. With absent modalities at `0.0` (L003) the background has no spread in the modality columns. Before 1.1.0 the channel then (a) standardized a present modality against that zero background with the 0.25 MAD floor, so any star with an ordinary light curve or spectrum had `|z| ~ 10-15`, and (b) its Isolation Forest skipped constant attributes, so it could never isolate a star on a modality feature. Benchmarks whose only anomaly is a modality feature (`EPIC 204278916`: dip; `VAN-Villarroel-1`: vanishing) were missed (blind recovery 4/6), while every star with spectra looked spectrally anomalous. The noise-model significances and the Liu et al. attribute draw remove both effects.
- The noise models are analytic (Gaussian extremes, chi^2). Real light curves with red noise or bad cadences have heavier tails, so modality-driven leads are routed to the corresponding domain vetting suite (`|contribution| >= 2`) and never stand on the outlier score alone. No background star in the current fixtures or live samples carries these modalities, so an empirical (population) calibration of the modality columns is not available yet.
- `spec_line_snr` is the strongest single-exposure emission feature (not the minimum of the two exposures): single-exposure cosmic rays become leads that `cosmic_ray` vetting rejects, which is what the frozen impostor `SL-CosmicRay-1` requires.

## Lessons Applied
- Channel produces leads only; every driving feature triggers full domain vetting.
- `L003`: Missing observational modalities (`lightcurves`, `spectra`, `vanishing_data`, `neowise`, `radio`) map deterministically to `0.0` rather than random Gaussian noise; `vet_target` always runs `metallicity` and `evolved` checks on flagged `generic_outlier` stars and only invokes modality-specific vetting when real observations exist.
