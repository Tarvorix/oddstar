# Channel: `generic_outlier`

## Physical Signature
Hypothesis-free discovery of targets occupying sparse/isolated regions of joint multi-channel photometric, astrometric, spectral, and time-domain feature space.

## Inputs
- Joint feature vector combining colors (`BP-RP`, `J-Ks`, `Ks-W3`, `W3-W4`), channel statistics (`chi_w4`, `ir_var_chi2`, `dip_snr`, `delta_m_sigma`, `spec_line_snr`, `vanished_margin`, `radio_snr`), and `ruwe`.

## Algorithm
1. Conditionally standardize features by median and MAD within `BP-RP` color bins.
2. Score with seeded Isolation Forest (`n_trees >= 200`, `subsample_size = 256`) and `cKDTree`-backed Local Outlier Factor (`k = 20`).
3. Flag targets falling in the top fraction (`>= 0.99` quantile) of both detectors simultaneously, and output per-feature standardized contributions.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.generic_outlier`:
- `iforest_n_trees: 100`, `lof_k_neighbors: 20`, `min_combined_quantile: 0.99`

## Known Impostors
- Combined impostors from every channel; flagged stars automatically inherit vetting checks tied to driving features.

## Benchmark Objects
- Positives from the other seven channels (`BD+20 307`, `TYC 4209-1322-1`, `KIC 8462852`, `ASASSN-21qj`, `EPIC 204278916`, `VAN-Villarroel-1`).

## Failure Modes
- Unconditioned colors bias cool/hot stars into outlier tails; conditional `BP-RP` standardization eliminates this artifact.

## Lessons Applied
- Channel produces leads only; every driving feature triggers full domain vetting.
