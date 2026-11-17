# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

## [1.1.0] - 2026-11-17
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
- Note: the top-level `version:` field is still "1.0.0"; bump it if the project ties it to this changelog.

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.
