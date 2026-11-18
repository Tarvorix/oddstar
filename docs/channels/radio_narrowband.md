# Channel: `radio_narrowband`

## Physical Signature
Narrowband, non-zero Doppler-drifting radio signals persistent across every ON-target scan and absent from every OFF-target scan.

## Inputs
- Breakthrough Listen open-data filterbank/HDF5 ABACAD cadences from Green Bank and Parkes (built-in cadences for
  validation; downloaded fine-resolution scans for the live search, see below).

## Algorithm
1. Execute incoherent de-doppler drift search (`dedoppler_search_scan`) over `[-4.0, +4.0] Hz/s`.
2. Require hit presence (`SNR >= 10`) in all ON scans following consistent drift trajectory `f(t)` and zero detections in OFF scans.
3. Enforce process-and-evict cache policy on raw filterbank blobs unless flagged.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.radio_narrowband`:
- `max_abs_drift_hz_s: 4.0`, `min_drift_abs_hz_s: 0.05`, `min_snr: 10.0`, `freq_match_tol_hz: 500.0`

## Live survey search (config 1.5.0)
`oddstar search radio_narrowband --live` works on real Breakthrough Listen data in three steps. Only the last one
decides anything, and it is the channel code above.

1. **Which stars have data** (`python -m oddstar.io.bl_opendata`). The archive's file index is asked for every
   sample position. A star counts as covered only inside the half-power beam of a pointing (radius 0.6 wavelength
   / dish diameter). `--plan` lists, for each covered star, the fine-resolution files that contain it (ON) and the
   comparison scans: the nearest scan before and after on another target with the same receiver, or other beams of
   the Parkes multibeam recorded at the same instant. Only the index is read.
2. **Download** (`python -m oddstar.io.bl_radio --plan <bl_scan_plan.json> --source-id <id>`). One scan set is tens
   of gigabytes, so nothing is downloaded by a search; this command fetches one star's set into
   `data/breakthrough_listen_scans`, resumes a broken transfer and checks size and MD5.
3. **Search** (`oddstar/io/bl_radio.py`). The ON scan is read in blocks of `live_block_channels`. In each block every
   stretch of `live_window_channels` is tested with the channel's statistic and drift grid; a stretch that reaches
   `live_prefilter_snr` is cut from the ON scan and, at the same sky frequencies, from every OFF scan, and that
   small cadence is measured and scored by the channel. Frequencies above the one where the star leaves the
   half-power beam are not searched. The star's record keeps the cadence of one window (flagged outside the
   allocated transmitter bands first, then flagged, then the strongest ON-only excess) and the list of evaluated
   windows; `in_rfi_band` is set from the allocation table in `bl_radio.RFI_ALLOCATIONS_MHZ`.

Not-searched reasons, each counted separately: `n_no_coverage` (inside no beam), `n_not_downloaded` (inside a beam
with a full scan set in the archive, files not on disk), `n_no_comparison` (inside a beam, no OFF scan in the
archive), `n_failed` (index or file unreadable).

Known limits of the live search:
- A scan set with one ON scan cannot test persistence. Interference that happened to be on during the ON scan only
  is flagged; the `reobservation` check stays `not_evaluated`, which caps the tier at C.
- The trial drift rates are 0.05 Hz/s apart while one channel over one scan is 0.01 Hz/s, so a line between two
  trial rates is smeared over a few channels and loses sensitivity.
- One window gives one hit per scan (its strongest). A weak line within 5.7 kHz of a stronger one is not seen, and
  an OFF-scan counterpart weaker than another line in its window is missed.
- A star off the pointing centre is observed at reduced gain (half at the beam edge), which is not corrected.
- The certificate of this channel is measured on the built-in cadences, not on archive scans. The false-alarm rate
  of the live search has not been measured; a flag is a lead.

## Known Impostors
- Terrestrial/satellite RFI (`blc1`), zero-drift ground transmitters, and multi-target session interference.

## Benchmark Objects
- `Voyager 1` (2015-12-30 Green Bank X-band observation).

## Failure Modes
- One-off detections without independent second-epoch reobservation cannot exceed Tier C.

## Lessons Applied
- Enforced `reobservation` gate capping single-epoch cadences strictly at Tier C.
