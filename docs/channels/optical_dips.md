# Channel: `optical_dips`

## Physical Signature
Aperiodic, asymmetric, or unusually deep optical dimming events and multi-year secular fading inconsistent with planetary transits or regular eclipsing binaries.

## Inputs
- Two independent TESS extractions (public QLP / SPOC and pipeline TESScut aperture photometry), Kepler/K2 light curves, ZTF, and ASAS-SN.

## Algorithm
1. Cut the light curve into segments where the sector changes or cadences are more than `tess_segment_gap_days` apart.
2. Detrend each segment alone with an iterative upper-clipped sliding median filter preserving multi-day dips.
3. Run a multi-duration boxcar matched filter (6 h to 3.5 days) inside each segment, leaving out the first and last
   `systematics_edge_margin_days`. Depth is the median deficit inside the best box; asymmetry is measured around its deepest point.
4. Compute a Box-Least-Squares (BLS) periodogram on the sector holding the best box to veto strictly periodic eclipsing binaries.
5. Fit the long-baseline ZTF linear fading slope when ZTF photometry is loaded.

## Live survey search (config 1.2.0)
`oddstar search optical_dips --live` fetches a TESS light curve for every sample star before measuring
(`oddstar/io/tess_bulk.py`): the QLP `SAP_FLUX` file of up to `tess_bulk_max_sectors` sectors from the public MAST bucket,
good-quality cadences averaged into 30-minute bins, each sector divided by its median. The cache keeps the reduced series
with the URL, size and SHA-256 of the file it came from, and remembers sectors that have no QLP file, so an interrupted run
resumes. `--tess-max-stars N` limits the fetch; `--defer-enrichment` skips the per-star vetting fetch for flagged stars
(`vet --live` does it later), which is how to look at the flag rate first. `vet --live` fetches three flagged stars at a time (`--fetch-workers`). Measured on 2026-11-19 with one star at a
time: about 4 minutes of archive queries per star that is not cached, and up to 10 more when the ZTF service hangs
(the fetch gives up after 600 s and records `query_failed`; a later run asks again).

A star is **searched** only when at least one segment has 20 points inside the edge margins. Every other star is listed
with its reason and excluded from the null result:

| Reason | Meaning |
|---|---|
| `no_coverage` | no TIC entry, no TESS sector at the position, or no sector with a QLP light curve (QLP is magnitude-limited) |
| `failed` | a lookup or every existing sector failed after `tess_bulk_max_attempts`; not cached, retried on the next run |
| `not_fetched` | beyond `--tess-max-stars` |
| `no_light_curve` | no fetch record and no usable light curve |

The search uses one extraction (QLP). The second extraction (TESScut) and the background and centroid series are fetched
for flagged stars only, by vetting. `radio_narrowband` has no live data source yet and `search --live` refuses it
(`spectral_lines` has one since config 1.3.0, `vanishing` since 1.4.0).

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.optical_dips`:
- `dip_snr_threshold: 6.0`, `min_dip_depth_frac: 0.03`, `min_duration_hours: 6.0`
- `periodic_bls_veto_snr: 8.0`, `secular_fade_snr_threshold: 4.5`

## Known Impostors
- Kepler/TESS eclipsing binaries, blended background EBs in 21" TESS pixels, scattered-light/momentum-dump artifacts, and detrending-only glitches.

## Benchmark Objects
- `KIC 8462852`, `EPIC 204278916`, `ASASSN-21qj`.
- Open case workflow target: `HD 139139`.

## Failure Modes
- Large 21" TESS pixels blend nearby EBs; `pipeline_agreement` and `neighbour_signal` centroid checks are mandatory.

## Lessons Applied
- Require dual-extraction agreement (`QLP` and `TESScut`) before vetting passes. Since `L006`, live vetting fetches QLP `SAP_FLUX` and an independent TESScut aperture extraction per sector, and `pipeline_agreement` applies the identical `detrend_window_days` running median to both at matched cadences before comparing `min_duration_hours`-binned dip depths.
- `second_survey` (since `L006`) uses live ZTF (retried, shrinking cone, yearly paging, PM-track matched) and ASAS-SN Sky Patrol seasonal photometry (`pip install oddstar[surveys]`).

## Light-curve vetting rules (config 1.1.0)
Every check below evaluates fetched data or returns `not_evaluated` naming what was missing or which query failed, and
every check can `fail`. Synthetic fixtures set explicit flags (`neighbour_eb_centroid_shift`, `tess_scattered_light_window`,
`sso_crossing_t0`); a flag that is present counts as evaluated, an absent flag never does. All inputs live under
`lightcurves.tess_sectors` / `lightcurves.optical_surveys` / `lightcurves.tess_neighbours` / `lightcurves.sso_queries`, so
channel scores do not change.

- **`pipeline_agreement`**: both series get the same `detrend_window_days` running median (segments split at
  `tess_segment_gap_days`), then `min_duration_hours` bins; bins with fewer than `pipeline_min_bin_fill_frac` of the expected
  native cadences are dropped. In each sector the deepest QLP bin is compared with TESScut **in the same bin**, and the deepest
  TESScut bin with QLP in the same bin, using `pipeline_quiet_depth_frac`, `pipeline_depth_abs_tol_frac` and
  `pipeline_depth_min_ratio` (the former literals 0.02 / 0.08 / 0.5). When either series has a feature above the quiet level the
  binned Pearson r must be >= `pipeline_min_binned_corr`. Per sector the evidence records the BTJD, both depths, the cadence count,
  the distance from the nearest segment edge and whether `min_dip_depth_frac` is exceeded for the deepest QLP, TESScut and
  coincident (min of both depths) bins. The note names the sector and time of the deepest coincident feature. Both extractions
  come from the same FFI pixels, so a pass validates the extraction, not the astrophysical origin.
- **`systematics_window`** (TESS part): tested bins are the deepest coincident bin and every bin above `min_dip_depth_frac`. A bin
  above threshold fails when it coincides with an annulus-background excursion (>= `systematics_bkg_sigma` robust sigma and
  >= `systematics_bkg_depth_frac` x the depth, in units of the target flux), lies within `systematics_edge_margin_days` of a segment
  edge, or has a 30-min centroid step >= `systematics_centroid_jump_sigma` robust sigma and >= `systematics_centroid_jump_min_px`
  within one bin of it. Sub-threshold deepest bins are reported, never failed. The NEOWISE part is described in
  `ir_timedomain.md`; `not_evaluated` only when neither part has data.
- **`neighbour_signal`** (TESS part): for every bin whose coincident depth exceeds `min_dip_depth_frac`, the flux-weighted centroid
  (background-subtracted, 2.5 px) inside the bin is compared with the out-of-dip centroid within `neighbour_centroid_window_days`.
  The flux-deficit (difference-image) location is `c_out - (1-D)/D (c_in - c_out)`; the raw centroid moves *away* from the source
  of a dip, the deficit location sits *on* it. Fail when the deficit lies >= `neighbour_source_offset_min_px` from the target at
  >= `neighbour_source_offset_sigma`; the nearest Gaia DR3 source within `vetting.tess_neighbour_radius_px` pixels is named. With
  TESS loaded and no dip above threshold the TESS part passes ("no dip to attribute to a neighbour").
- **`sso_crossing`**: SkyBoT (observer `C57`) is queried during enrichment at the start, middle and end of the deepest coincident
  bin of every sector whose coincident depth exceeds `min_dip_depth_frac`, within `vetting.sso_tess_radius_px` x
  `vetting.tess_pixel_arcsec`. Fail when a body with flux >= `vetting.sso_min_flux_ratio_of_depth` x the depth (relative to G; a body
  without V counts) is inside the radius; `not_evaluated` when a query failed or is missing. With TESS sectors loaded and no event
  above threshold the check passes ("no SSO-caused TESS event"); with no TESS data it is `not_evaluated`.
- **`known_variable`**: SIMBAD class veto as before, plus per sector a linear-per-segment detrend (rotation preserved), Lomb-Scargle
  between `rotation_min_period_days` and `rotation_max_period_days` (Baluev FAP; detections below `rotation_fap_max`, consistency
  within `rotation_period_match_frac`) and BLS for eclipses. A periodic EB needs BLS SNR > `periodic_bls_veto_snr`, depth >=
  `periodic_eb_min_depth_frac`, and >= 2 observed eclipses covering at least half of the predicted epochs with data. Spot
  rotation is reported, never failed. No SIMBAD class and no light curve gives `not_evaluated`.
- **`second_survey`** (TESS-dip context): ZTF/ASAS-SN nightly medians within the event bin (+/- `second_survey_simultaneous_pad_days`)
  fail the check when the survey depth plus `second_survey_night_dip_sigma` sigma is below
  `second_survey_simultaneous_min_depth_ratio` x the TESS depth, and confirm it otherwise. Without a simultaneous night, an
  independent dimming episode (>= `second_survey_min_dip_nights` consecutive nights fainter by `min_dip_depth_frac` at
  `second_survey_night_dip_sigma`) confirms the star dims; otherwise `not_evaluated`. The IR context is in `ir_timedomain.md`.
- Limitation: the matched running median uses `mode="nearest"` padding (as the channel does), which absorbs a dip that runs into
  the last cadence of a segment.

## The flagged dip must be in the vetted light curves (config 1.13.0)
The search reads the QLP extraction of every sector. The vetting light curves need the QLP and the TESScut extraction
on matched cadences, so a sector the search used can be missing from them (TESScut returned no usable cadence, a
fetch failed). Until config 1.13.0 the checks then judged the sectors they had, found them quiet and passed.

`vet` and `inspect` now hand the search's dip to the checks as `lightcurves.search_dip` (sector, box centre, box
width, depth, S/N), taken from the search measurement when it holds a dip above the channel thresholds in a real
sector. `search_dip_gap` lays the box over the kept 6 h bins of that sector in the vetting light curves. When less
than half of the box (`SEARCH_DIP_MIN_BOX_COVERAGE`) lies in kept bins:

- `pipeline_agreement`, `systematics_window`, `neighbour_signal` and `sso_crossing` return `not_evaluated` instead of
  `pass`. The note names the sector, the box, the vetted sectors and the fetch log line of the missing sector; the
  metrics hold `search_dip_gap`. A check that fails on the sectors it has still fails.
- `second_survey` tests the box itself: a ZTF or ASAS-SN night inside the box (no padding) is compared with the box
  depth by the rule used for vetted events. A night at normal brightness fails the check; no night inside the box
  leaves it `not_evaluated`.

The case that led to the rule: Gaia DR3 23767631757133952 in the nearby 1,300-star run (config 1.12.0), box depth
0.26 and S/N 10.09 in sector 42, where the star lies on two bad CCD columns (aperture sum negative, QLP scatter 0.21
to 0.38 per day against 0.001 to 0.005 in its other sectors). TESScut returned no cadence with positive flux for that
sector, vetting loaded sectors 43, 70 and 71, four checks passed and the star was Tier C. ZTF observed it inside the
box at normal brightness. Of the 68 flagged stars of that run it was the only one whose flagged sector was not vetted.

Not covered: the search still accepts such a sector. A sector whose flux is noise on a level near zero should not be
searched at all; that rule needs a threshold and a measurement on real sectors.
