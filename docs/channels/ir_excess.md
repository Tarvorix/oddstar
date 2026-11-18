# Channel: `ir_excess`

## Physical Signature
Mid-infrared excess emission above a star's photospheric spectrum at 12 um (WISE W3) and 22 um (WISE W4), probing circumstellar warm dust (~100–700 K) or re-radiated waste heat at terrestrial-to-asteroid-belt separations.

## Inputs
- Gaia DR3/DR4 astrometry (`ra`, `dec`, `parallax`, `ruwe`) and `G`, `BP`, `RP` photometry.
- 2MASS `J`, `H`, `Ks` photometry.
- AllWISE `W1`–`W4` photometry, profile-fit `chi2`, contamination/extension flags, and 2010 cryogenic single-exposure frames.

## Algorithm
1. Apply parent-sample cuts (`parallax_snr > 10`, `distance_pc <= 300`, `ruwe < 1.4`, `|b| > 10 deg`, unique AllWISE cross-match, clean 2MASS).
2. **Primary photosphere (all excess quantification): empirical colour locus.** Dereddened `(Ks - W1)`, `(Ks - W2)`, `(Ks - W3)`, `(Ks - W4)` against `(BP - RP)_0`, all corrected with the star's applied reddening `ebv_star` (R_V = 3.1, `physics.BAND_EXTINCTION_AV`). Per BP-RP bin and band the locus stores the median after iterative median/MAD clipping at `locus_clip_sigma` (so stars with excess do not bias it), the robust scatter, the median photometric error, the intrinsic scatter `sqrt(max(floor^2, scatter^2 - median(sigma_Ks^2 + sigma_W^2)))`, and N. A star's prediction is `W = Ks_0 - (Ks - W)_0 + A_W` with `sigma_pred = hypot(sigma_Ks, intrinsic scatter)` (plus `locus_extrapolation_scatter_mag` outside the locus coverage).
   - Source, recorded on every measurement as `photosphere_source` together with `locus_n_stars` and `locus_id`:
     - `sample_locus`: built from the stars of the current sample. After a sample run the locus (bins, medians, scatter, N, Gaia release, sample kind `synthetic`/`live`, config hash, sample size) is written to `runs/_photosphere_loci/ir_excess_locus-<release>-<kind>-<locus_id>.json` when it uses at least `locus_persist_min_stars` stars.
     - `persisted_locus`: single-star mode with no sample loads the most recent persisted locus for the same Gaia release and sample kind (a live star never uses a synthetic-fixture locus).
     - `literature_fallback`: only if neither exists, the Pecaut & Mamajek (2013) mean dwarf sequence packaged in `src/oddstar/refdata/pecaut_mamajek_dwarf_colours.json` (Ks-W1 A0V-M5V; W1-W2/W3/W4 tabulated only B5V-K5V, so redder stars are clamped and flagged). The per-band scatter is the assumed `literature_locus_scatter_mag`. If the packaged table were ever replaced by values not retrieved from the source (`provenance != "retrieved"`), the label becomes `literature_fallback_unverified`. Retrieval: 2026-11-17 through a web-fetch tool that passes text through a summarising model (three fetches agree); diff against the original file before publication.
3. **Secondary photosphere (diagnostic and second detection leg).** Planck fit to dereddened `G, BP, RP, J, H, Ks`, then `physics.stellar_photosphere_flux_mjy`, which imposes the tabulated dwarf `Ks - Wx` colour at the fitted Teff on the Planck shape. The bare Planck predictions are also reported (`w1..w4_pred_planck_mjy`). A pure Planck curve anchored at Ks overpredicts a K dwarf by ~13% at W1 and ~30% at W2.
4. Per band W1..W4: observed flux, primary photosphere with uncertainty, excess with uncertainty, excess significance and `obs/phot` ratio (`w*_obs_mjy`, `w*_phot_mjy`, `w*_excess_mjy`, `w*_excess_sigma`, `w*_excess_ratio`; uncertainties in `Measurement.uncertainties`).
5. Detection (unchanged structure): `chi_w3`, `chi_w4`, `w3_ratio`, `w4_ratio` against the secondary model and `chi_emp_w3`, `chi_emp_w4`, `w3_ratio_emp`, `w4_ratio_emp` against the primary locus; with `require_both_predictors` a band must pass under both. The empirical chi uses `hypot(sigma_obs, sigma_pred)`; the locus intrinsic scatter already contains the sample's star-to-star calibration scatter, so `cal_uncertainty_frac` is applied to the secondary model only.
6. Dust on the primary excess: single-temperature blackbody fitted only to bands whose excess is significant at `excess_significance_sigma` (>= 2 bands needed), with `T_dust` and `L_dust/L_star` intervals from a chi2 profile (`dust_temp_interval_delta_chi2`, widened by the reduced chi2 when one blackbody fits badly); separately the colour temperatures of the W1-W2 excess and of the W3-W4 excess with intervals (`t_color_w1w2_*`, `t_color_w3w4_*`), a two-temperature fractional luminosity (`frac_lum_hot + frac_lum_warm`), and the Suazo et al. (2024) partial-shell grid. The same fit against the secondary model is reported as `t_dust_model_k`, `frac_lum_model`.
7. If the line-of-sight dust column exceeds the applied reddening by more than `vetting.extinction_los_excess_summary_ebv`, the measurement carries `extinction_cases` (intrinsic BP-RP, M_G and dwarf class under the applied and the full line-of-sight reddening); the dossier prints it.

## Synthetic Fixtures
Fixture photospheres use the same `physics.stellar_photosphere_flux_mjy` as the secondary model (Planck shape plus the tabulated dwarf WISE colours), so both predictors and the literature fallback describe the fixtures consistently, and the fixture sample locus reproduces the dwarf colours. Injected dust (`dust_blackbody_flux_mjy`) is normalised on the Planck-shaped stellar SED and is unaffected by the colour correction.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.ir_excess`:
- `chi_w3_threshold: 5.0`, `chi_w4_threshold: 5.0`
- `min_flux_ratio_w4: 2.0`, `min_flux_ratio_w3: 1.5`
- `cal_uncertainty_frac: 0.05`, `require_both_predictors: true`
- Locus: `locus_bp_rp_min`, `locus_bp_rp_max`, `locus_n_bins`, `locus_min_stars_per_bin`, `locus_clip_sigma`, `locus_clip_max_iter`, `locus_intrinsic_scatter_floor_mag`, `locus_extrapolation_scatter_mag`, `locus_persist_min_stars`, `literature_locus_scatter_mag`
- Dust: `excess_significance_sigma`, `dust_temp_interval_delta_chi2`, `colour_temp_min_k`, `colour_temp_max_k`

## Infrared Vetting Inputs (config 1.1.0)
Every check in `oddstar.vet.ir_checks` is tri-state. For live stars an input is used only when the query that produces
it is recorded as `"ok"` in `vetting_archive.query_status[<query>]` (otherwise `"no_coverage"` or
`"failed: <Type>: <msg>"`); a missing or failed query gives `not_evaluated` and the note names the query and its status.
Synthetic fixtures (`synthetic_fixture: true`) carry explicit inputs and are the only stars that may use the synthetic
cutout renderer, the fixture centroid offset, or the stellar `ebv_sfd` as a cirrus proxy.

| Check | Live input (query name) | Pass / fail rule |
|---|---|---|
| `centroid_offset` | W3/W4 Gaussian + tilted-plane centroids on AllWISE Atlas cutouts (`cutouts_wise`) | fail if offset > `centroid_offset_max_beam_frac` x beam (and constraining); pass needs a constraining W4 centroid. The AllWISE catalogue position (W1/W2-dominated) is reported only and leaves the check `not_evaluated`. |
| `chance_alignment` | AllWISE sources with W3/W4 S/N >= `chance_alignment_min_snr` and flux >= `chance_alignment_min_flux_frac_of_excess` x excess in a 15"-300" annulus (`allwise_density`) | P = 1 - exp(-n pi r^2) per band, r = max(k sigma_centroid, frac x beam) or beam/2 without a constraining centroid; fail if P > `chance_alignment_p_max`. |
| `spherex_colours` | The stored SPHEREx comparison of the star (`vetting_archive.spherex_colours`, status `spherex_colours`), measured by `oddstar spherex verdict`: the star and its nearest comparison stars on SPHEREx detectors 3 to 6 | Since 1.10.0. The rule of `docs/spherex.md`: SPHEREx colour minus catalogue colour (K - W2, W1 - W2), against the median and robust scatter of the same difference over the comparison stars. pass on "reproduced" and on "redder now" (reported as a change); fail on "not reproduced" (bluer in SPHEREx by more than `spherex_z_limit`); `not_evaluated` without a stored measurement or when the rule gives no verdict. Speaks for 2 to 5 micron only. |
| `youth` | GALEX NUV (`galex`, MAST `Mast.Galex.Catalog`; `no_coverage` when no GALEX source lies within `youth_galex_coverage_radius_arcsec`), ROSAT 2RXS (`xray_2rxs`, HEASARC `rass2rxs`), eRASS1 main (`xray_erass1`, HEASARC `erass1main`, Galactic l 180-360 only), SIMBAD `h_link` parents (`young_association`) | fail on any youth indicator; pass only when at least two of three tests say "not young": UV, X-ray (an X-ray non-detection whose survey-limit L_X exceeds `youth_xray_lx_min_erg_s` is not a pass) and, since 1.8.0, kinematics (tangential velocity of `youth_old_vtan_min_kms` or more with a parallax measured to `youth_old_vtan_min_parallax_snr`; a slow star gives no verdict). The NUV-G inactive locus in `ir_checks.NUV_G_LOCUS` is approximate and carries `youth_nuv_locus_sigma_mag`. |
| `moving_object` | IMCCE SkyBoT (observer C51) at the median epoch of each 2010 L1b visit, propagated to every frame (`skybot`, `wise_l1b`) | fail on a body within `sso_match_radius_arcsec` (plus its ephemeris uncertainty, when that is at most `sso_max_pos_err_arcsec`) at any frame; pass only when every visit query succeeded. Bodies with a larger uncertainty are tested at their nominal position and listed as `unconstrained_bodies`. |
| `vision_rubric` | Real AllWISE W1-W4 Atlas cutouts (IRSA SIA + IBE) and a PS1 r stack for Dec > -30 | fail on a centroid offset beyond both k sigma and the beam fraction, fitted FWHM / the measured point-source FWHM of the band (`w3_atlas_point_fwhm_arcsec`, `w4_atlas_point_fwhm_arcsec`; 63 ordinary stars, `validation/calibration/atlas_point_source_fwhm.json`) > `vision_extended_fwhm_ratio`, or a brighter secondary W4 peak; the evidence PNG is the real panel. |
| `cirrus` | IRSA DUST reference-pixel SFD E(B-V) (`irsa_dust`, stored in `vetting_archive.dust`) and measured W3/W4 Atlas background rms (mJy per point source) | fail if either exceeds `cirrus_ebv_max` / `cirrus_bg_rms_max_mjy`; pass needs both. |
| `deep_imaging` | PS1 DR2 / LS DR10 catalogues (neighbour screen) and MAST CAOM HST/JWST image coverage (`subarcsec_imaging_mast`) | catalogues can only fail (non-PSF source in the core or with comparable unWISE flux); pass requires an assessed sub-arcsec image. Existing but unassessed HST/JWST imaging is listed and leaves `not_evaluated`. |
| `neighbours` | Gaia DR3 sources inside 12" with W3/W4 flux estimated from RP (or G) assuming M-dwarf colours `neighbour_assumed_*` (`gaia_neighbours`) | fail if the summed flux fraction > `neighbour_blend_flux_frac`. |
| `extragalactic_colours` | total WISE colours (rule unchanged) | photosphere-subtracted excess colours are reported only; hot circumstellar dust shares the AGN colour region. |

## Known Impostors
- Radio-loud and radio-quiet dusty background galaxies inside the 12" WISE W4 beam (Ren et al. 2024, 2026).
- Dust-obscured galaxies (DOGs) and red neighbours.
- Galactic dust cirrus filaments, AllWISE diffraction spikes (`cc_flags`), and single-frame glitches.

## Benchmark Objects
- `BD+20 307`, `ID8`, `HD 23514`, `TYC 4209-1322-1`, `HD 113766`, `HD 172555`, `HD 69830`, and the six extreme debris disks in Moor et al. (2021).
- Open case workflow target: `Gaia DR3 4582145718700356096`.

## Failure Modes
- W3 and W4 photometry exist solely from the 2010 cryogenic WISE mission (single epoch, single instrument).
- Over 90% of background galaxies contaminating WISE beams are radio-quiet, so passing `radio_counterpart` never clears `background_galaxy` without deep sub-arcsec optical/NIR imaging.

## Lessons Applied
- `L001`: Single-epoch 22 um W4 excesses lacking deep sub-arcsec imaging and second-instrument 22 um follow-up remain capped at Tier C with `background_galaxy` and `deep_imaging` flagged explicitly as open holes.
- `L001` (1.1.0): seeing-limited PS1/Legacy Survey catalogues are a neighbour screen only and can never clear `deep_imaging`.
- `L003`: Live Gaia ingest requires `w4mpro_error IS NOT NULL` to prevent `W4 = 'U'` upper-limit leakage, propagates stellar proper motion (`2016.0 -> 2010.2`) before computing `centroid_offset_arcsec`, checks time-resolved NEOWISE proper-motion tracking across visits, and scales `radial_velocity_error * sqrt(2N/pi)` (marking faint `rv_method_used = 2` stars without p-values as `binary = not_evaluated`).
