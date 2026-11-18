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
   the Parkes multibeam recorded at the same instant. ON and comparison files are always of the same product
   (`.hires.hdf` files differ in channel width from `.h5` and `.fil` files). Only the index is read.
   The index lists many southern declinations wrongly: the fraction of a degree is mirrored, differently for each
   file family (`bl_opendata.true_declination`, measured against catalogue positions of named targets). A southern
   star is therefore asked for at every declination the index may use, and each file is put back at its true
   declination before any separation is computed. A southern file of a family that was not measured is
   "unverified": it never makes a star covered, and a star that only such a file could contain is counted as
   `unverified_pointing`.
2. **Download** (`python -m oddstar.io.bl_radio --plan <bl_scan_plan.json> --source-id <id>`). One scan set is tens
   of gigabytes, so nothing is downloaded by a search; this command fetches one star's set into
   `data/breakthrough_listen_scans`, resumes a broken transfer and checks size and MD5. With `--search` it instead
   downloads, searches and deletes one scan at a time (disk for one file is enough) and stores the result where
   the survey finds it, so the star counts as searched without its files on disk.
3. **Search** (`oddstar/io/bl_radio.py`). The ON scan is read in blocks of `live_block_channels`. In each block every
   stretch of `live_window_channels` is tested with the channel's statistic and drift grid; a stretch that reaches
   `live_prefilter_snr` is cut from the ON scan and, at the same sky frequencies, from every OFF scan, and that
   small cadence is measured and scored by the channel. Frequencies above the one where the star leaves the
   half-power beam are not searched. The star's record keeps the cadence of one window (flagged outside the
   allocated transmitter bands first, then flagged, then the strongest ON-only excess) and the list of evaluated
   windows; `in_rfi_band` is set from the allocation table in `bl_radio.RFI_ALLOCATIONS_MHZ`.
   From 698 to 4200 MHz the table lists every band the United States table of frequency allocations gives to a
   transmitting service, so `rfi_band` can pass there only in the passive bands 1400 to 1427, 1660.5 to 1668.4 and
   2690 to 2700 MHz; above 4200 MHz it lists two satellite downlink bands only and is not complete (rule and cost
   in `configs/CHANGELOG.md`).
   Stored results are keyed by `bl_radio.SEARCH_VERSION`; a change of the search, or of the table inside the
   frequency range of a stored result, changes that version.
   The windows cut by `--search` stay under `data/breakthrough_listen_scans/_work`, keyed by `bl_radio.CUT_VERSION`:
   after a change that only touches the scoring, the same command scores them again without fetching a scan.
   They can be deleted once no such change is expected (about 1.6 GB per scan of a Green Bank L-band set).

   **Cadences.** When the archive holds several ON scans of one target that follow each other (A, B, A, C, A), the
   plan joins them (`comparison: "cadence"`): candidates are found in the middle ON scan, the same sky frequencies
   are cut from the other ON scans and from every OFF scan between and around them, and the channel then requires
   the signal in every ON scan and in no OFF scan (`require_all_on_scans`). A lone ON scan is compared with the
   OFF scans before and after it. The cut around a hit is `live_window_channels` wide, or the next power-of-two
   multiple when a signal at `max_abs_drift_hz_s` would leave that between the start of the reference scan and
   the end of the farthest scan (three 5-minute ON scans 10.6 minutes apart: 4096 channels).

Not-searched reasons, each counted separately: `n_no_coverage` (inside no beam), `n_not_downloaded` (inside a beam
with a full scan set in the archive, files not on disk), `n_no_comparison` (inside a beam, no OFF scan in the
archive), `n_failed` (index or file unreadable).

Known limits of the live search:
- Pointings with a true declination between -1 and 0 degrees cannot be recovered from the index and are not found.
  The declination rules rest on a few named targets per file family (listed in `bl_opendata.py`).
- The index and the Parkes file servers answer over plain HTTP only. A machine that can reach HTTPS only can read
  Green Bank files (`bldata.berkeley.edu`) but neither the index nor Parkes files.
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
