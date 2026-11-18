# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

## [1.5.0] - 2026-11-18

### The radio search reads real Breakthrough Listen scans
Until now `search radio_narrowband --live` was refused: the channel had only its built-in cadences. A fine-resolution
archive scan is one block of 16 spectra by hundreds of millions of channels, which the channel cannot take whole.

- New search of scans on local disk (`oddstar/io/bl_radio.py`, used by `search radio_narrowband --live`). It walks the
  ON scan in blocks, cuts every window that stands out from the ON and OFF scans, and hands each to the channel as an
  ordinary cadence. See `docs/channels/radio_narrowband.md` for the steps, the not-searched reasons
  (`n_no_coverage`, `n_not_downloaded`, `n_no_comparison`, `n_failed`) and the known limits.
- New keys under `channels.radio_narrowband`, all chosen before any archive scan was searched:
  - `live_block_channels: 1048576`. One coarse channel of the Green Bank fine-resolution product; memory only.
  - `live_window_channels: 2048`. 5.7 kHz at 2.8 Hz per channel: wider than the largest drift (4 Hz/s) carries a
    line over one scan plus the gap to the neighbouring scan, so the OFF window holds where the line would be.
  - `live_prefilter_snr: 8.0`. Below the channel's `min_snr: 10.0`, so a window is never lost because the block
    search and the channel estimate the noise from slightly different stretches. Noise alone reached 6.6 in a
    1,048,576-channel test block.
  - `live_max_listed_windows: 200` limits the list kept per scan; the counts are of all windows.
- The channel code did not change, and no threshold changed. `vet/radio_checks.py` changed in wording only: a live
  star counts its inputs when `query_status["breakthrough_listen"]` is "ok", which the live search now writes; the
  note of a missing input says whether no scans were compared or the search does not produce that input.
  `multi_target` and `reobservation` stay `not_evaluated` for a live star.
- `search --live` no longer refuses any channel.
- The archive index applies its row limit before cutting to the requested circle, so the earlier limit of 500 could
  drop files silently. Requests now use 20000 and an answer of 5000 rows or more is refused.

Known limits: the certificates are measured on the built-in cadences; the false-alarm rate of the live search on
archive scans is not measured; a scan set with a single ON scan cannot test persistence, so interference present
only during that scan is flagged and stays a lead until a second epoch exists.

## [1.4.0] - 2026-11-18

### The vanishing search reads real plate sources
Until now no code fetched a plate record for a live object and `search vanishing --live` was refused (1.2.0). The
channel had only ever run on fixtures. A live sample star is in Gaia by construction and cannot vanish from it, so the
live search has a different unit: the USNO-B1.0 plate sources around each sample position.

- New survey fetch (`oddstar/io/plate_survey.py`, used by `search vanishing --live`): USNO-B1.0, Gaia DR3 and
  Pan-STARRS1 DR2 for a field around every sample position. See `docs/channels/vanishing.md` for the steps and the
  not-searched reasons (`n_too_faint`, `n_no_coverage`, `n_failed`).
- New keys under `channels.vanishing`, all chosen before any field was searched:
  - `plate_field_radius_arcsec: 60.0`. A high-latitude field of this size holds a handful of plate sources; 1,000
    fields are 0.87 square degrees.
  - `plate_ps1_min_detections: 2`. A Pan-STARRS object with fewer detections is not taken as a counterpart or into the
    depth. A low value is the cautious choice for a search that flags absence: a doubtful modern object still counts
    as "something is there".
  - `plate_depth_percentile: 90.0` and `plate_depth_min_objects: 10`. The local depth is the 90th percentile of the
    r-band PSF magnitudes of the field's Pan-STARRS objects (22.6 in the one field looked at to check the query
    format). It is where detections become rare, not where they are complete; the existing `min_depth_margin_mag: 2.0`
    is applied on top, so a searched plate source is at least 2 mag brighter than that.
  - `plate_max_attempts: 3`, `plate_backoff_s: 2.0`, `plate_rate_limit_rps: 2.0` are operational.
- The channel code did not change. Thresholds `min_early_plates`, `min_depth_margin_mag`, `match_radius_arcsec` and
  `max_pm_residual_arcsec` are as in 1.3.0. `vet/vanishing_checks.py` changed in its description only (which inputs a
  live object can have); every certificate value of all eight channels is identical to 1.3.0.
- `search --live` now refuses only `radio_narrowband`.

Known limits:
- The certificate describes a Gaia star paired with a plate record. The rate of false USNO-B1.0 entries (plate flaws,
  halos of bright stars, merged detections) is not part of it, so the false-alarm rate of a live run is not measured
  and its flags are leads. The null report of such a run carries the fixture grid.
- The two first-epoch plates are independent emulsions, usually of the same night, not independent epochs.
- `sso_plate`, `blend_split` and `known_variable` have no live input and stay `not_evaluated`, so a live flag cannot
  clear every check.
- South of Dec -30 there is no Pan-STARRS coverage and nothing is searched.
- `version: "1.4.0"` (was "1.3.0"). The config hash changes; every certificate is regenerated in a separate commit.

## [1.3.0] - 2026-11-18

### The spectral-line search reads real survey spectra
Until now no code fetched a spectrum for a live star, and `search spectral_lines --live` was refused (1.2.0). The
channel had only ever run on fixtures.

- New survey fetch (`oddstar/io/lamost_bulk.py`, used by `search spectral_lines --live`): two epochs per star from the
  LAMOST DR12 v2.0 low-resolution survey, reduced to the window around H-alpha in the stellar rest frame. See
  `docs/channels/spectral_lines.md` for the steps and the not-searched reasons (`n_single_epoch`, `n_no_coverage`,
  `n_failed`, `n_not_fetched`, `n_no_spectrum`). Only a star with two usable epochs from different nights is searched.
- New keys under `channels.spectral_lines`: `lamost_min_snr_r: 20.0` is a selection cut (an epoch with a lower red
  signal-to-noise is not used; the usual quality cut for LAMOST low-resolution work). It was chosen before any star was
  searched. `lamost_max_attempts: 3`, `lamost_backoff_s: 2.0`, `lamost_rate_limit_rps: 2.0` are operational.
- `SpectralLinesChannel.measure` now follows the line-spread width it is given. The H-alpha templates are broadened in
  quadrature by the excess over 0.45 A (the width they were defined at). A known line is excluded out to
  1.5 x `lsf_fwhm_a` when that exceeds the configured width: a matched filter still answers a line one FWHM away with
  a quarter of its peak, so at R = 1800 the configured 1.5 A let a strong sky line flag from just outside its zone.
  When `spectra.sky_line_wavelengths_a` is present, those positions are excluded (all of them, geocoronal H-alpha
  included) in place of the laboratory list; a spectrum in the stellar rest frame has its sky lines displaced by the
  star's velocity. Exposure arrays of unequal length are not measured. `sky_line` and `stellar_emission` in vetting
  use the same list and width. No threshold changed (`residual_snr_threshold`, the 0.4 repeat ratio, both exclusion
  widths and every `vetting.spectral_*` value are as in 1.2.0).
- Effect on validation: none. The fixtures have `lsf_fwhm_a = 0.45` and no `sky_line_wavelengths_a`, so every code path
  above reduces to the 1.2.0 arithmetic. Every certificate value of all eight channels was compared with 1.2.0: blind
  recovery, impostor tiers, completeness grids and false-alarm rates are identical, and one `spectral_lines`
  calibration-curve statistic per release moved by about 1e-14 (rounding of the rewritten template expression).
- `search --live` still refuses `vanishing` and `radio_narrowband`.

Known limits:
- The certificate describes the channel at a line-spread width of 0.45 A. Nothing was injected at R = 1800: the
  completeness and false-alarm rate of a LAMOST search are not measured, and its flags are leads, not calibrated
  detections. The null report of such a run carries the fixture grid and must not be read as a LAMOST sensitivity.
- The public LAMOST file has no sky spectrum per fibre and no detector map: `sky_line` can fail but not pass, and
  `detector` is `not_evaluated`, so a LAMOST flag cannot clear every check.
- About half of the 150 A window is excluded around sky and stellar lines at this resolution.
- `version: "1.3.0"` (was "1.2.0"). The config hash changes; every certificate is regenerated in a separate commit.
  Runs made under 1.2.0 (`v120-*`) keep their hash; to fuse channels they are run again under 1.3.0, from the cache.

## [1.2.0] - 2026-11-18

Two changes. Both were decided from how the pipeline behaved, not from which stars it flagged, and neither can
turn a failed or open check into a pass.

### The optical-dip search reads real TESS light curves for every sample star
Until now a live `search optical_dips` measured empty arrays for every star and reported "searched N, flagged 0":
a live sample star only received a light curve after some other channel had flagged it. The channel was therefore
never run on real data.

- New survey fetch (`oddstar/io/tess_bulk.py`, used by `search optical_dips --live`): the QLP `SAP_FLUX` light curve
  of each sector, good-quality cadences averaged into 30-minute bins, each sector divided by its median. New
  operational keys under `channels.optical_dips` (not thresholds): `tess_bulk_max_sectors: 8` (the sector limit of the
  vetting fetch, so vetting sees the sectors the search saw), `tess_bulk_max_attempts: 3`, `tess_bulk_backoff_s: 2.0`,
  `tess_bulk_workers: 4`, `tess_bulk_rate_limit_rps: 4.0`.
- A star counts as searched only when at least one light-curve segment has 20 points inside the edge margins. A star
  with no TIC entry, no sector, no QLP file, a failed fetch, or one beyond `--tess-max-stars` is reported as not
  searched, with the reason, and is excluded from the null result (`n_no_coverage`, `n_failed`, `n_not_fetched`,
  `n_no_light_curve`).
- `OpticalDipsChannel.measure` now works on segments. The light curve is cut where the sector changes or cadences
  are more than `tess_segment_gap_days` apart, and each segment is detrended and searched alone. The old code ran one
  median filter and one box filter across orbit gaps and sector boundaries, where the per-sector normalisation leaves
  a step. The first and last `systematics_edge_margin_days` of a segment are not searched: thermal settling and
  scattered light sit there, and `systematics_window` cannot clear an event in that zone. `dip_depth` is the median
  deficit inside the best box (was: the single lowest point within a few cadences of the peak), so one bad cadence is
  not a dip. The periodicity veto runs on the sector that holds the best box. No threshold changed
  (`dip_snr_threshold`, `min_dip_depth_frac`, `min_duration_hours`, `periodic_bls_veto_snr` are as in 1.1.2).
- Effect on validation (fixtures, dr3 and dr4 identical): blind recovery 3/3, impostors at most Explained, false-alarm
  rate 0, design point (0.15 depth) 100%, all unchanged. Completeness at the shallowest injected depth fell, because
  the median inside a box wider than the dip is below `min_dip_depth_frac`:

  | duration (days) | depth | 1.1.2 | 1.2.0 |
  |---|---|---|---|
  | 0.5 | 0.04 | 1.00 | 0.70 |
  | 1.5 | 0.04 | 1.00 | 0.85 |

  All other grid cells are unchanged (3.0 days at 0.04 was and is 0.00). The loss is real and is the price of not
  flagging single low points; the null report carries the new grid.
- `search --live` now refuses `spectral_lines`, `vanishing` and `radio_narrowband`. No live ingestion provides their
  inputs, so a live run would again report a search that did not happen (`radio_narrowband --live` crashed instead).

### A star that varies in NEOWISE cannot pass `systematics_window` on its TESS data
The first live `ir_timedomain` runs (1,800 stars, 30 flagged) showed the NEOWISE half of `systematics_window`
almost never ran: the field-star ensemble had 0 or 1 stars for nearly every target, the half returned
`not_evaluated`, and the check passed on the TESS half alone. A check that has not looked at the infrared frames was
passing stars flagged for an infrared change. Eight stars of one run shared one faint W2 trend, which is what an
untested zero-point drift would look like.

- `check_systematics_window`: when the star's NEOWISE light curve varies by the channel's own thresholds and the
  NEOWISE half is `not_evaluated`, the check is `not_evaluated`, whatever the TESS half says. A fixture's explicit
  `scan_systematic` flag still counts as an evaluation. This can only move a result from `pass` to `not_evaluated`.
- `channels.ir_timedomain.neowise_comparison_radius_arcsec: 600.0` (was 150.0). A 150" cone around a W1 = 11.9 halo
  star holds one other star within `neowise_comparison_mag_window` (1.0 mag); 450" holds seven (AllWISE count around
  Gaia DR3 3580358134642783872, made to size the cone, not to judge that star). The zero-point and scan systematics
  the ensemble tests are common to a 47' frame, so stars 10' away see them. `neowise_comparison_mag_window`,
  `neowise_comparison_min_stars`, `neowise_comparison_max_stars`, `neowise_comparison_exclude_arcsec` and the three
  `systematics_ensemble_*` thresholds are unchanged.
- The cone is filtered by the archive to the target's W1 +/- (window + 0.25 mag) and to six columns; unfiltered it
  would be hundreds of thousands of single exposures. The margin keeps a star at the edge of the window from losing
  part of its exposures.
- Stars brighter than about W1 = 8 will usually still have no comparison stars. They now stay `not_evaluated`.

Known limits:
- Results of the 1.1.2 `ir_timedomain` runs (`v113-*`) were vetted with the old rule; their `systematics_window`
  passes for varying stars are not evidence. They must be vetted again under 1.2.0 before any tier is quoted.
- `neighbour_signal` combines its TESS and NEOWISE halves the same way and was not changed here.
- `version: "1.2.0"` (was "1.1.2"). The config hash changes; every certificate is regenerated in a separate commit.

## [1.1.2] - 2026-11-17

### The image-shape check is judged against measured stars
`vision_rubric` calls a band "extended" when its fitted FWHM exceeds a reference by `vision_extended_fwhm_ratio`.
In 1.1.0 and 1.1.1 the reference was theoretical: native beam x `wise_atlas_psf_broadening` (9.19" in W3, 16.97"
in W4). The first live run measured the reference star at 1.40 x that in W3 and 1.31 x in W4 and failed it.
Whether that meant extended emission or a reference that does not describe this fit could only be settled on
known point sources, so the fit was run on ordinary stars.

Calibration (`python -m oddstar.validate.atlas_psf`, output committed in ab488cb before this change):
63 bright, isolated, unsaturated stars with photospheric W3-W4 in 12 fixed fields at |b| > 30 deg, none near a
pipeline candidate, each measured by the unchanged `measure_wise_cutouts`. All 12 selection queries and all 63
cutout sets succeeded.

| Band | n | median fitted FWHM | robust sigma | min - max | median / theoretical reference |
|---|---|---|---|---|---|
| W3 | 63 | 12.971" | 0.099" | 12.85" - 13.32" | 1.41 |
| W4 | 63 | 20.597" | 0.487" | 17.84" - 23.76" | 1.21 |

So an ordinary star comes out of this fit 1.41 x (W3) and 1.21 x (W4) wider than the theoretical reference, and the
1.1.1 check would have called essentially every star extended in W3.

The selection and the rule for using the result were written in `oddstar/validate/atlas_psf.py` and committed
(160a8f9) before the calibration ran. Applied as written:

- `vetting.w3_atlas_point_fwhm_arcsec: 12.971` and `vetting.w4_atlas_point_fwhm_arcsec: 20.597` (new): the median
  fitted FWHM of the calibration stars. `fwhm_ratio` is now fitted FWHM over this value. The geometry of the
  measurement (aperture, fit box, peak search) still uses native beam x `wise_atlas_psf_broadening`, so fitted
  widths are unchanged and stay comparable with the calibration file.
- `vetting.vision_extended_fwhm_ratio: 1.3` (unchanged). The rule was max(1.3, 1 + 4 x robust sigma / median) over
  the two bands, rounded up to 0.05; 1 + 4 sigma / median is 1.03 in W3 and 1.09 in W4, so the rule returns 1.3.
- `version: "1.1.2"` (was "1.1.1"). The config hash changes; every certificate is regenerated in a separate commit.

Known limits, recorded so they are not mistaken for results:
- 1.3 is loose against the measured scatter (about 30 sigma in W3, 12 sigma in W4). The rule could only keep or
  raise the threshold. Tightening it is a separate decision that must be made on calibration data alone.
- Fitted width rises slightly with S/N (Spearman rho +0.43 in W3, +0.32 in W4), and the W4 scatter is larger for
  faint stars: the 19 calibration stars with fit S/N below 25 span 17.8" - 23.8" (robust sigma 0.81").
- The calibration stars have catalogue W4 S/N of 10 or more. A fainter target is outside the calibrated range.

No search result is part of this change.

## [1.1.1] - 2026-11-17

### Fixes from the first live archive run of 1.1.0
The 2026-11-17 live run of the reference star against the real archives exposed answers the 1.1.0 code misread.
No channel threshold changes. Three keys are new; the config hash changes, so every certificate is regenerated
(in a separate commit).

- `version: "1.1.1"` (was "1.1.0").
- `vetting.sso_max_pos_err_arcsec: 30.0` (new). SkyBoT lists lost and single-opposition bodies whose ephemeris
  uncertainty reaches degrees (seen live: 2021 VU98, 225,012" uncertainty, nominal position 561" from the star).
  The match test adds a body's uncertainty to `sso_match_radius_arcsec`, so such a body matched any target in the
  field and failed `moving_object`. A body with an uncertainty above this cap is now tested at its nominal position
  only and listed under `unconstrained_bodies` in the check's metrics and note. 30" is twice the 15" match radius
  and 2.5 times the 12" W4 beam: beyond that the catalogue position no longer says whether the body was in the beam.
- `vetting.youth_galex_coverage_radius_arcsec: 600.0` (new). GALEX did not observe the whole sky. A field with no
  GALEX source of any kind within this radius was not observed (`galex: no_coverage`), which is different from an
  observed field with no NUV source at the star. An observed AIS field holds of order a hundred catalogue sources
  within 10 arcmin.
- `channels.ir_timedomain.neowise_bulk_workers: 4` (new). Concurrent IRSA requests of the survey NEOWISE fetch; all
  workers share the `archives.neowise.rate_limit_rps` limit. Results do not depend on the worker count (tested).

Source changes without a config key (recorded here because they change which archive answers a check):
- GALEX now comes from MAST `Mast.Galex.Catalog` and ROSAT 2RXS / eRASS1 from the HEASARC TAP tables `rass2rxs` /
  `erass1main`. The CfA VizieR mirror answered HTTP 200 with zero rows for GALEX and eRASS1 everywhere (it does not
  serve those tables), which was recorded as "queried, no counterpart", and the 2RXS identifier `IX/47/2rxs` does
  not exist. Match radii and flux limits are unchanged.
- The remaining VizieR `asu-tsv` tables (NVSS, PS1 fallback) accept an empty answer from a mirror only after that
  mirror returned rows at a known position for the table (`VIZIER_CANARY_CONES` in `io/live_survey.py`).
- The Legacy Survey DR10 neighbour query is a 17" `q3c_radial_query` cone (circumscribes the former +/-12" box,
  which timed out); rows are still filtered by separation downstream.
- SkyBoT HTTP 204 is "no body in the field", not a failure. IRSA image-search columns are read by name. A MAST
  answer counts only when its status is `COMPLETE`. `_youth_uv` dereddens with the reddening applied to the star
  (`ebv_star`), not the line-of-sight key, which is empty for live stars since 1.1.0.
- Operational only (no effect on results): Bayestar19 stops being asked for 30 minutes after a failed connection;
  a TAP service answering HTTP 429 is retried after 20, 60, 180 and 480 s; halo and bulge samples skip the
  single-shot join that timed out and use the index + batch route unless a single-shot answer is already cached.

## [1.1.0] - 2026-11-17

### Infrared vetting inputs
Measured infrared vetting inputs (`vetting:` section only; channel thresholds unchanged). Config hash changes from
`ce6e224fd31d` to `9d274caf416b`, so channels need re-validation before search (certificates not regenerated here).
- `chance_alignment_p_max: 0.01` - fail threshold for P(unrelated W3/W4 source within the confusion radius); replaces the hard-coded 0.01.
- `chance_alignment_annulus_inner_arcsec: 15.0`, `chance_alignment_annulus_outer_arcsec: 300.0` - AllWISE annulus for the local density (inner radius excludes the target's own blend).
- `chance_alignment_min_snr: 3.0` - W3/W4 S/N for a source to count towards the density.
- `chance_alignment_min_flux_frac_of_excess: 0.5` - only sources at least this fraction of the target's excess could mimic it.
- `youth_xray_lx_min_erg_s: 1.0e+29` - X-ray youth threshold; was hard-coded in `ir_checks.py`.
- `youth_nuv_locus_sigma_mag: 0.5` - adopted uncertainty of the approximate inactive NUV-G locus.
- `youth_uv_constraint_nsigma: 2.0` - UV test is constraining only if excess + n sigma stays below `youth_uv_excess_mag`.
- `youth_galex_match_radius_arcsec: 3.0`, `youth_rosat_match_radius_arcsec: 30.0`, `youth_erass1_match_radius_arcsec: 15.0` - cross-match radii (GALEX ~5" PSF; ROSAT PSPC survey positions ~10-30"; eRASS1 ~5-15").
- `youth_rosat_flux_limit_erg_cm2_s: 1.0e-13`, `youth_erass1_flux_limit_erg_cm2_s: 5.0e-14` - survey flux limits used to turn a non-detection into an L_X upper limit (approximate all-sky typical values).
- `sso_visit_gap_days: 3.0` - L1b frames closer than this belong to one visit (one SkyBoT query per visit).
- `sso_max_rate_arcsec_per_hr: 120.0` - widens each visit's SkyBoT cone so bodies reaching the target within the visit are returned.
- `w1_beam_arcsec: 6.1`, `w2_beam_arcsec: 6.4` - native W1/W2 PSF FWHM for cutout measurements (W3 6.5", W4 12.0" already present).
- `wise_atlas_psf_broadening: 1.414` - Atlas images are PSF-convolved, so the reference FWHM is ~sqrt(2) x native.
- `vision_cutout_size_arcsec: 200.0` - Atlas cutout size (room for the 30-90" background annulus).
- `vision_centroid_max_sigma: 3.0` - k in the k sigma centroid test and the chance-alignment radius.
- `vision_extended_fwhm_ratio: 1.3` - fitted / reference FWHM above which a band is called extended.
- `vision_min_fit_snr: 5.0` - minimum point-source S/N for the rubric to evaluate a band.
- `vision_peak_snr_min: 3.0` - per-pixel significance for local maxima (secondary peaks).
- `cirrus_annulus_inner_arcsec: 30.0`, `cirrus_annulus_outer_arcsec: 90.0` - blank annulus for the measured background rms.
- `deep_screen_core_radius_arcsec: 1.5` - non-PSF catalogue source inside this radius is a contaminant.
- `deep_screen_ps1_psf_minus_kron_min: 0.35` - PS1 r_PSF - r_Kron above which a source is extended.
- `deep_screen_faint_r_mag: 17.5` - faint point sources 1.8-8" fainter than this are recorded (no longer an automatic failure).
- `deep_screen_min_flux_frac_of_excess: 0.5` - a non-PSF LS source in the beam is a contaminant if its unWISE W3/W4 flux reaches this fraction of the excess.
- `subarcsec_imaging_search_radius_arcsec: 3.0` - MAST CAOM HST/JWST search radius.
- `neighbour_assumed_g_minus_rp: 1.3`, `neighbour_assumed_rp_minus_w3: 3.7`, `neighbour_assumed_rp_minus_w4: 4.0` - conservative (red, M-dwarf-like) colours used to estimate a Gaia neighbour's W3/W4 flux.
- Note: the top-level `version:` field was still "1.0.0" when the entries above were written; it is bumped to "1.1.0" below.

### Governance and thresholds moved from source
Decision cuts that were literals in the source now live in `configs/default.yaml`. **Every value is unchanged**, so
no vetting outcome, sample or channel score changes; the live sample SQL text is byte-identical (tested in
`tests/test_governance_config.py`). The config hash changes from `10bc507f1b29` to `c4ad1ae04d32`, so validation
certificates must be regenerated for the new hash (not done here).

- `version: "1.1.0"` (was "1.0.0"). CI now requires `version` to equal the newest `## [x.y.z]` heading here (`check_config_version`).
- R4 commit-separation cutoff: `7bf6ed5f4be3050e91b99f87503bbb20078b9604`. `check_git_commit_separation` checks HEAD and, in pull-request CI, every non-merge commit of the base range, but only commits after this one; earlier history is not rewritten.
- `vetting` (from `vet/engine.py` `vet_target`): `outlier_driver_min_abs_contribution: 2.0` (generic_outlier driver contribution that routes to a check).
- `vetting` (from `vet/ir_checks.py`):
  - `check_xmatch_integrity`: `xmatch_max_sep_arcsec: 2.5`.
  - `check_wise_quality`: `wise_w4_min_frames: 5`, `wise_w4_min_detected_frame_frac: 0.15`, `wise_w4_max_rchi2: 3.5`, `wise_max_ext_flg: 1`.
  - `check_extragalactic_colours`: `dog_w1_w2_min: 0.45`, `dog_w2_w3_min: 3.0`.
  - `check_single_frames`: `single_frame_consistency_frac: 0.5`, `single_frame_coadd_floor_frac: 0.6`, `single_frame_max_to_median: 2.2`.
  - `_youth_uv`: `youth_nuv_default_err_mag: 0.3` (NUV error adopted when GALEX gives none).
  - `check_evolved`: `evolved_teff_max_k: 5200.0`.
  - `check_binary`: `binary_rv_scatter_max_kms: 5.0`.
- `vetting` (from `vet/lc_checks.py`):
  - `matched_detrend` / `analyze_tess_sector`: `lc_detrend_low_mask_sigma: 2.5` (the function keeps 2.5 as its default argument for direct callers).
  - `_tess_systematics`: `lc_tess_contiguous_dt_factor: 1.5`.
  - `_sector_period_search` and `check_known_variable` (BLS): `lc_bls_min_period_days: 0.8`, `lc_bls_max_period_days: 8.0`, `lc_bls_duration_days: 0.1`.
  - `_sector_period_search` (eclipse confirmation): `lc_eb_eclipse_min_depth_frac_of_model: 0.5`, `lc_eb_min_eclipses: 2`, `lc_eb_min_eclipse_epoch_frac: 0.5`.
  - `_ir_event_window`: `lc_ir_event_baseline_visits: 3`, `lc_ir_event_dev_frac: 0.5`, `lc_ir_event_min_window_days: 180.0`, `lc_ir_event_pad_days: 90.0`.
- `vetting` (from `vet/spectral_checks.py`):
  - `check_cosmic_ray`: `spectral_cr_spike_mad_sigma: 8.0`, `spectral_peak_width_flux_level: 1.3`, `spectral_cr_min_peak_flux: 1.5`, `spectral_cr_min_width_lsf_frac: 0.7`.
  - `check_sky_line`: `spectral_sky_peak_mad_sigma: 6.0`, `spectral_target_peak_sigma: 3.5`.
  - `check_stellar_emission`: `spectral_stellar_line_min_amp: 0.20`.
  - `check_repeat`: `spectral_repeat_quiet_max_excess: 0.25`, `spectral_repeat_max_pixel_shift: 2`, `spectral_repeat_min_excess: 0.20`.
- `vetting` (from `vet/vanishing_checks.py` `check_proper_motion`): `vanishing_pm_min_raw_sep_arcsec: 3.0`.
- `fusion` (from `fuse/engine.py` `calibrate_fap` / `fuse_results`): `fap_floor: 1.0e-12`, `fallback_false_alarm_rate_per_star: 1.0e-4` (`calibrate_fap` keeps `floor=1e-12` as its default argument for callers that pass none).
- `live_ingestion` (new section):
  - from `io/live_survey.py` `_phot_quality_where_sql`: `allwise_number_of_neighbours: 1`, `allwise_number_of_mates: 0`, `tmass_number_of_neighbours: 1`, `tmass_number_of_mates: 0`.
  - from `io/live_survey.py` `fetch_live_halo_bulge_sample` region presets, `regions.halo`: `min_parallax_over_error: 10.0`, `parallax_min_mas: 0.5`, `parallax_max_mas: 20.0`, `max_ruwe: 1.4`, `min_abs_b_deg: 20.0`, `g_mag_min: 7.5`, `g_mag_max: 14.5`, `teff_min_k: 3800.0`, `teff_max_k: 7500.0`, `max_mh_gspphot: -0.8`, `min_vtan_kms: 150.0`.
  - `regions.bulge`: `gal_l_half_width_deg: 12.0` (the query keeps `l <= 12.0 OR l >= 348.0`), `abs_b_min_deg: 3.0`, `abs_b_max_deg: 10.0`, `parallax_min_mas: 0.08`, `parallax_max_mas: 0.45`, `min_parallax_over_error: 5.0`, `max_ruwe: 1.4`, `g_mag_min: 8.0`, `g_mag_max: 14.5`, `teff_min_k: 3500.0`, `teff_max_k: 7500.0`.
  - `regions.local`: `min_parallax_over_error: 15.0`, `parallax_min_mas: 2.0`, `parallax_max_mas: 25.0`, `max_ruwe: 1.4`, `min_abs_b_deg: 12.0`, `g_mag_min: 7.5`, `g_mag_max: 13.5`, `teff_min_k: 3800.0`, `teff_max_k: 7500.0`, `min_vtan_kms: 65.0`.
  - from `cli.py` `filtered_live_ingestion`: `ir_excess_w4_upper_limit_ph_qual: ["U", "X"]` (W4 `ph_qual` codes treated as upper limits), `ir_excess_artifact_cc_flag_bands: [3, 4]` (bands whose `cc_flags` character marks an artifact), `ir_timedomain_neowise_enrich_max_stars: 25` (NEOWISE enrichment cap).
- `storage.evidence_dir: "registry/evidence"` and a new `evidence` section (`max_array_elements: 64`, `max_string_chars: 2000`, `max_json_bytes: 65536`, `max_image_bytes: 262144`): size guards of the committed evidence store. Storage policy, not decision thresholds.
- Deliberately left in source: statistical and unit constants (1.4826, 1.2533, 2.5 log10, 4.74047, 365.25, MJD/BTJD offsets, catalogue epochs 2016.0 / 2010.2, R_V 3.1, physical constants), numerical guards (1e-9, 1e-6, 1e-4 floors), array indices and window offsets, minimum sample sizes for a statistic to be computed (points per sector or segment, bins, seasons, nights), the BLS grid size 200, the `0.5 x span` period bound, data defaults used when an input is missing (spectral `lsf`/pixel size, `abs_g`/`teff` defaults, USNO epoch/magnitude, PS1 depth, plate count), enrichment-fetch query radii, timeouts, chunk sizes and the `n_lim` clamp in `live_survey.py`, `_row_to_star_dict` and the end of `enrich_star_live_for_vetting` (out of scope here), `vet/underluminous_checks.py`, and the dossier writer's threshold text.

### Star model: empirical photosphere and explicit extinction
- Added `channels.ir_excess` empirical-locus keys (moved from hard-coded values where they existed; no existing threshold changed): `locus_bp_rp_min: -0.2`, `locus_bp_rp_max: 2.8`, `locus_n_bins: 15`, `locus_min_stars_per_bin: 5` (same binning as before), `locus_clip_sigma: 3.0`, `locus_clip_max_iter: 5`, `locus_intrinsic_scatter_floor_mag: 0.01`, `locus_extrapolation_scatter_mag: 0.10`, `locus_persist_min_stars: 100`, `literature_locus_scatter_mag: [0.03, 0.04, 0.05, 0.08]` (assumed W1..W4 scatter for the literature fallback; the table has none).
- Added `channels.ir_excess` dust-quantification keys: `excess_significance_sigma: 3.0`, `dust_temp_interval_delta_chi2: 1.0`, `colour_temp_min_k: 30.0`, `colour_temp_max_k: 3000.0`.
- Added `vetting` extinction keys: `extinction_heavy_dust_ebv: 0.30` (previously hard-coded in `check_extinction`, same value), `extinction_los_excess_summary_ebv: 0.05`, `sfd_to_sf11_factor: 0.86`, `bayestar19_ebv_per_unit: 0.884` (Green et al. 2019 Eq. 1; Eq. 2 gives 0.996; to be confirmed in a live run).
- Behaviour: the ir_excess primary photosphere is now the dereddened W1..W4 empirical locus (`photosphere_source` recorded); the empirical-predictor chi now includes the locus intrinsic scatter (was multiplied by 0.0) and no longer adds `cal_uncertainty_frac` (already inside the locus scatter); fixture photospheres carry the dwarf WISE colours; `ebv_star` replaces `ebv_sfd` wherever the star is dereddened.

### Light-curve vetting inputs
Light-curve vetting checks now evaluate fetched data or return `not_evaluated`, and every check can fail. All keys are
vetting-only additions; no channel detection threshold changed, so channel scores are unaffected. The config hash
changes, so validation certificates must be regenerated for the new hash.

- `channels.optical_dips` (pipeline_agreement): `pipeline_quiet_depth_frac: 0.02`, `pipeline_depth_abs_tol_frac: 0.08`,
  `pipeline_depth_min_ratio: 0.5` (moved from literals in `lc_checks._depths_agree`, values unchanged);
  `pipeline_min_bin_fill_frac: 0.5` (6 h bins with fewer than half the expected native cadences are dropped);
  `pipeline_min_binned_corr: 0.5` (binned QLP/TESScut Pearson r required only when a feature exceeds the quiet level);
  `tess_segment_gap_days: 0.5` (segment split used by `matched_detrend` and edge distances; was the default argument).
- `channels.optical_dips` (systematics_window, TESS): `systematics_edge_margin_days: 0.5`, `systematics_bkg_sigma: 5.0`,
  `systematics_bkg_depth_frac: 0.25`, `systematics_centroid_jump_sigma: 5.0`, `systematics_centroid_jump_min_px: 0.05`.
- `channels.optical_dips` (neighbour_signal, TESS): `neighbour_centroid_window_days: 2.0`, `neighbour_source_offset_sigma: 3.0`,
  `neighbour_source_offset_min_px: 0.3`.
- `channels.optical_dips` (known_variable): `periodic_eb_min_depth_frac: 0.15` (was the literal `min(flux) < 0.85`, unchanged);
  `rotation_min_period_days: 0.1`, `rotation_max_period_days: 13.0`, `rotation_fap_max: 0.001`, `rotation_period_match_frac: 0.1`.
- `channels.optical_dips` (second_survey, TESS events): `second_survey_simultaneous_pad_days: 0.25`,
  `second_survey_simultaneous_min_depth_ratio: 0.5`, `second_survey_night_dip_sigma: 3.0`, `second_survey_min_dip_nights: 2`.
- `channels.ir_timedomain` (systematics_window, NEOWISE): `neowise_comparison_radius_arcsec: 150.0`,
  `neowise_comparison_exclude_arcsec: 8.0`, `neowise_comparison_mag_window: 1.0`, `neowise_comparison_min_stars: 3`,
  `neowise_comparison_max_stars: 10`, `systematics_ensemble_corr: 0.7`, `systematics_ensemble_corr_min_amp_frac: 0.2`,
  `systematics_ensemble_amp_frac: 0.5`.
- `channels.ir_timedomain` (neighbour_signal, NEOWISE): `neighbour_blend_sigma: 3.0`, `neighbour_min_track_arcsec: 0.3`
  (the stationary-flux-share limit reuses `vetting.neighbour_blend_flux_frac: 0.25`).
- `channels.ir_timedomain` (second_survey): `second_survey_min_season_nights: 5` (was the literal in `summarize_survey`),
  `second_survey_season_floor_mag: 0.005`, `second_survey_change_sigma: 5.0`, `second_survey_min_overlap_seasons: 4`,
  `second_survey_min_abs_corr: 0.6`, `second_survey_slope_sigma: 3.0`, `second_survey_comparable_frac: 0.5`.
- `vetting`: `tess_neighbour_radius_px: 2.0` (Gaia cone for TESS-aperture neighbours), `sso_tess_radius_px: 2.0`
  (SkyBoT radius at TESS event times = 2 x `tess_pixel_arcsec`), `sso_min_flux_ratio_of_depth: 0.1`.

### generic_outlier re-certification, code fingerprints, survey-scale NEOWISE (wt/survey-scale)
Keys added only; no existing threshold or gate changed. The config hash changes, so every channel must be re-validated;
certificates now also carry `code_fingerprint` and go stale when the validated code changes.
- `channels.generic_outlier.modality_null_default_n_points: 1000` - series length assumed for the noise expectation of a
  `dip_snr` / `spec_line_snr` value supplied without its series (`sub_metrics`), i.e. E[max of 1000 N(0,1)] = 3.23 (Blom approximation).
- `channels.generic_outlier.vanished_margin_plate_sigma_mag: 0.3` - photographic-plate magnitude uncertainty (USNO-B1.0,
  Monet et al. 2003) used to express `vanished_margin` in sigma.
- Behaviour (generic_outlier): a present modality feature is now an excess over its own noise expectation (absent stays 0.0,
  L003); the Isolation Forest draws split attributes from all features (Liu et al. 2008), so a feature constant in the
  background can still isolate a query; `spec_line_snr` is the strongest single-exposure feature (repeat test left to the
  spectral vetting suite, like `dip_snr`). See `docs/channels/generic_outlier.md`.
- `channels.ir_timedomain.neowise_bulk_max_attempts: 3`, `neowise_bulk_retry_backoff_s: 2.0` - bounded retries (exponential
  backoff 2 s, 4 s) of the survey NEOWISE-only fetch (`oddstar.io.neowise_bulk`); its rate limit is `archives.neowise.rate_limit_rps`.
- Behaviour (ir_timedomain): stars without NEOWISE coverage, with a failed or missing NEOWISE fetch, or with fewer than
  `min_visits` visits are reported as not searched and excluded from null-report denominators (detection thresholds unchanged).
- Removed `live_ingestion.ir_timedomain_neowise_enrich_max_stars: 25`. A live `search ir_timedomain` used to fetch NEOWISE
  for the first 25 sample stars only and drop the rest; it now fetches every sample star through the bulk path
  (`--neowise-max-stars N` bounds a run explicitly, and the stars beyond N are reported as not searched).

### Remaining default-value passes removed
Vanishing, spectral and radio vetting, the SIMBAD-based checks and the generic_outlier path no longer pass on
default input values (lessons L007/L008): a check passes or fails only on inputs that are explicitly present (fixtures,
hand-built dicts) or written by a query recorded as `ok` (live stars); otherwise it reports `not_evaluated` and names the
missing data product. No existing threshold value changed. The config hash changes, so validation certificates must
be regenerated for the new hash (not done here).

- `channels.vanishing.vetting_variable_simbad_types: ["Mira", "LPV*", "CV*", "Flare*"]` - SIMBAD types that explain a
  vanished plate source in `known_variable`; was a literal set in `vet/vanishing_checks.py` (unchanged).
- `channels.optical_dips.detrend_edge_zone_window_frac: 0.5` (new) - dip candidates within this fraction of
  `detrend_window_days` (2.5 d) of a TESS segment edge or data gap, where the matched running median can erase or
  distort a dip, make `pipeline_agreement`, `systematics_window`, `neighbour_signal`, `sso_crossing` and the TESS context
  of `second_survey` report `not_evaluated` instead of `pass` (a `fail` stands). Edge-zone bins are also compared with the
  raw interior flux so an erased dip is still seen.
- `vetting.metallicity_subdwarf_feh_max: -1.0` - [Fe/H] below which `check_metallicity` calls a subdwarf; was the
  literal `-1.0` in `vet/underluminous_checks.py` (unchanged).
- Literal fallbacks of config keys removed from vetting source (`.get(key, default)` -> required key, values unchanged):
  `channels.vanishing.min_early_plates` (2), `min_depth_margin_mag` (2.0), `max_pm_residual_arcsec` (1.5),
  `channels.radio_narrowband.min_drift_abs_hz_s` (0.05), `channels.spectral_lines.sky_line_exclusion_angstrom` (1.5),
  `stellar_line_exclusion_angstrom` (1.5).
- Input defaults removed (now `not_evaluated` when absent): vanishing `usno_epoch` 1953.4, `usno_mag_r` 17.0,
  `panstarrs_depth_mag` 21.8, `early_plates_count` (2 with any plate record), `known_lpv` False, `sso_on_plate` False,
  `blend_split_sources` 1, `present_in_gaia` True, proper motion 0.0; spectral `lsf_fwhm_a` 0.45, pixel size 0.25 A,
  `detector_artifact` False, empty flux arrays; radio `off_scan_hits` 0, `drift_hz_s` 0.0, `in_rfi_band` False,
  `multi_target_match` False, `reobserved_confirmed` False; SIMBAD type "Star"; generic path `abs_g_mag` 4.8,
  `teff_k` 5800, `ruwe` 1.0, proper motion 0.0 (xmatch, NEOWISE PM track), radio row separation 99",
  second-instrument interpretation "confirmed", fixture cutout offset 0.0" / background 0.4 mJy.
- SIMBAD: live ingestion records `query_status["simbad"]` and `vetting_archive.simbad_entry`; `simbad_type` is None when
  the lookup failed or SIMBAD has no entry. `known_variable` (light curve and vanishing) and `metallicity` return
  `not_evaluated` when the lookup did not succeed; "no entry" is a successful negative answer.
- `vet/engine.py`: a generic_outlier target with a dip/IR driver is routed to the light-curve checks also when it carries
  only per-sector TESS data (`lightcurves.tess_sectors`, the live format), not just the legacy fixture arrays.
- Fixture-only channels for live stars: vanishing, spectral_lines and radio_narrowband have no live ingestion, so every
  check of those suites returns `not_evaluated` for a live star (a SIMBAD-listed Mira can still fail vanishing
  `known_variable`).

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.
