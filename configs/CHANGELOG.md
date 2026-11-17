# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

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

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.

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
