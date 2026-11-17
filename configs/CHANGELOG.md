# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.

## [1.1.0] - 2026-11-17
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
