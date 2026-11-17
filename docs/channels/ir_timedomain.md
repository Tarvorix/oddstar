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
