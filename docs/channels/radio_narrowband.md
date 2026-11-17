# Channel: `radio_narrowband`

## Physical Signature
Narrowband, non-zero Doppler-drifting radio signals persistent across every ON-target scan and absent from every OFF-target scan.

## Inputs
- Breakthrough Listen open-data filterbank/HDF5 ABACAD cadences from Green Bank and Parkes.

## Algorithm
1. Execute incoherent de-doppler drift search (`dedoppler_search_scan`) over `[-4.0, +4.0] Hz/s`.
2. Require hit presence (`SNR >= 10`) in all ON scans following consistent drift trajectory `f(t)` and zero detections in OFF scans.
3. Enforce process-and-evict cache policy on raw filterbank blobs unless flagged.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.radio_narrowband`:
- `max_abs_drift_hz_s: 4.0`, `min_drift_abs_hz_s: 0.05`, `min_snr: 10.0`, `freq_match_tol_hz: 500.0`

## Known Impostors
- Terrestrial/satellite RFI (`blc1`), zero-drift ground transmitters, and multi-target session interference.

## Benchmark Objects
- `Voyager 1` (2015-12-30 Green Bank X-band observation).

## Failure Modes
- One-off detections without independent second-epoch reobservation cannot exceed Tier C.

## Lessons Applied
- Enforced `reobservation` gate capping single-epoch cadences strictly at Tier C.
