# Channel: `ir_excess`

## Physical Signature
Mid-infrared excess emission above a star's photospheric spectrum at 12 um (WISE W3) and 22 um (WISE W4), probing circumstellar warm dust (~100–700 K) or re-radiated waste heat at terrestrial-to-asteroid-belt separations.

## Inputs
- Gaia DR3/DR4 astrometry (`ra`, `dec`, `parallax`, `ruwe`) and `G`, `BP`, `RP` photometry.
- 2MASS `J`, `H`, `Ks` photometry.
- AllWISE `W1`–`W4` photometry, profile-fit `chi2`, contamination/extension flags, and 2010 cryogenic single-exposure frames.

## Algorithm
1. Apply parent-sample cuts (`parallax_snr > 10`, `distance_pc <= 300`, `ruwe < 1.4`, `|b| > 10 deg`, unique AllWISE cross-match, clean 2MASS).
2. Predict photospheric W3 and W4 fluxes independently via (a) weighted Planck SED fit to `G, BP, RP, J, H, Ks` and (b) empirical `(Ks - W3)` and `(Ks - W4)` vs `(BP - RP)` median/MAD locus.
3. Compute excess statistics `chi_w3`, `chi_w4` and flux ratios under both predictors.
4. Fit single-temperature dust blackbody `(T_dust, f_dust, r_bb)` and Suazo et al. (2024) partial-shell covering-fraction grid.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.ir_excess`:
- `chi_w3_threshold: 5.0`, `chi_w4_threshold: 5.0`
- `min_flux_ratio_w4: 2.0`, `min_flux_ratio_w3: 1.5`
- `cal_uncertainty_frac: 0.05`, `require_both_predictors: true`

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
| `youth` | GALEX AIS NUV (`galex`, VizieR II/335/galex_ais), ROSAT 2RXS (`xray_2rxs`, IX/47/2rxs), eRASS1 main (`xray_erass1`, J/A+A/682/A34/erass1-m, Galactic l 180-360 only), SIMBAD `h_link` parents (`young_association`) | fail on any youth indicator; pass only when both the UV and the X-ray tests are constraining (an X-ray non-detection whose survey-limit L_X exceeds `youth_xray_lx_min_erg_s` is not a pass). The NUV-G inactive locus in `ir_checks.NUV_G_LOCUS` is approximate and carries `youth_nuv_locus_sigma_mag`. |
| `moving_object` | IMCCE SkyBoT (observer C51) at the median epoch of each 2010 L1b visit, propagated to every frame (`skybot`, `wise_l1b`) | fail on a body within `sso_match_radius_arcsec` at any frame; pass only when every visit query succeeded. |
| `vision_rubric` | Real AllWISE W1-W4 Atlas cutouts (IRSA SIA + IBE) and a PS1 r stack for Dec > -30 | fail on a centroid offset beyond both k sigma and the beam fraction, FWHM / (native x `wise_atlas_psf_broadening`) > `vision_extended_fwhm_ratio`, or a brighter secondary W4 peak; the evidence PNG is the real panel. |
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
