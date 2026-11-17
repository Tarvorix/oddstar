# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

## [1.1.0] - 2026-11-17

### Star model: empirical photosphere and explicit extinction (worktree `wt/star-model`)
- Added `channels.ir_excess` empirical-locus keys (moved from hard-coded values where they existed; no existing threshold changed): `locus_bp_rp_min: -0.2`, `locus_bp_rp_max: 2.8`, `locus_n_bins: 15`, `locus_min_stars_per_bin: 5` (same binning as before), `locus_clip_sigma: 3.0`, `locus_clip_max_iter: 5`, `locus_intrinsic_scatter_floor_mag: 0.01`, `locus_extrapolation_scatter_mag: 0.10`, `locus_persist_min_stars: 100`, `literature_locus_scatter_mag: [0.03, 0.04, 0.05, 0.08]` (assumed W1..W4 scatter for the literature fallback; the table has none).
- Added `channels.ir_excess` dust-quantification keys: `excess_significance_sigma: 3.0`, `dust_temp_interval_delta_chi2: 1.0`, `colour_temp_min_k: 30.0`, `colour_temp_max_k: 3000.0`.
- Added `vetting` extinction keys: `extinction_heavy_dust_ebv: 0.30` (previously hard-coded in `check_extinction`, same value), `extinction_los_excess_summary_ebv: 0.05`, `sfd_to_sf11_factor: 0.86`, `bayestar19_ebv_per_unit: 0.884` (Green et al. 2019 Eq. 1; Eq. 2 gives 0.996; to be confirmed in a live run).
- Behaviour: the ir_excess primary photosphere is now the dereddened W1..W4 empirical locus (`photosphere_source` recorded); the empirical-predictor chi now includes the locus intrinsic scatter (was multiplied by 0.0) and no longer adds `cal_uncertainty_frac` (already inside the locus scatter); fixture photospheres carry the dwarf WISE colours; `ebv_star` replaces `ebv_sfd` wherever the star is dereddened.

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.
