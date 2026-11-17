# Channel: `optical_dips`

## Physical Signature
Aperiodic, asymmetric, or unusually deep optical dimming events and multi-year secular fading inconsistent with planetary transits or regular eclipsing binaries.

## Inputs
- Two independent TESS extractions (public QLP / SPOC and pipeline TESScut aperture photometry), Kepler/K2 light curves, ZTF, and ASAS-SN.

## Algorithm
1. Detrend light curves with an iterative upper-clipped sliding median filter preserving multi-day dips.
2. Compute Box-Least-Squares (BLS) periodogram to veto strictly periodic eclipsing binaries.
3. Run multi-duration boxcar matched filter for dip SNR, depth, ingress/egress asymmetry, and long-baseline ZTF linear fading slope.

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
