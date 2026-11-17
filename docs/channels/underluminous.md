# Channel: `underluminous`

## Physical Signature
Main-sequence stars appearing systematically fainter in extinction-corrected absolute magnitude than predicted by two independent stellar-parameter determinations.

## Inputs
- Gaia DR3/DR4 parallax, RUWE, `G` photometry, Gaia XP parameters, survey spectroscopic parameters (`teff_spec_k`, `logg`, `[Fe/H]`).
- Reddening, with explicit names:
  - `ebv_star` / `ebv_star_source`: the E(B-V) applied to the star (`gaia_gspphot` = Gaia GSP-Phot `ebpminrp_gspphot / 1.31` in live mode, `map3d_<name>`, `assumed_zero`, `fixture`, or `missing` with `ebv_star = None`; there is no silent default).
  - `ebv_sfd` / `ebv_planck`: true line-of-sight map values only (None in live rows until enrichment provides them).
  - `vetting_archive.dust` (`ebv_sfd_los`, `ebv_sf11_los`, `i100_mjy_sr`): line-of-sight values from vetting enrichment; compared on the Schlafly & Finkbeiner (2011) scale (`sfd_to_sf11_factor`).
  - `vetting_archive.dust3d`: Bayestar19 E(B-V) at the star's distance with uncertainty (`oddstar.io.dust3d`, Argonaut API, Dec > -30 deg only; `E(B-V) = bayestar19_ebv_per_unit x Bayestar19`, 0.884 per Green et al. 2019 Eq. 1, 0.996 per Eq. 2).

## Algorithm
1. Build empirical `M_G(Teff, logg, [Fe/H])` locus from the parent sample.
2. Compute `M_G,obs` corrected with `A_G` from `ebv_star` (no correction, with a caveat, when it is missing) and residual `Delta M = M_G,obs - M_G,pred` under both XP and spectroscopic parameters.
3. Flag stars faint under both parameter sources with clean astrometry (`RUWE < 1.4`, `parallax_snr > 10`).
4. When the line-of-sight column exceeds the applied reddening by more than `vetting.extinction_los_excess_summary_ebv`, report `extinction_cases`: intrinsic BP-RP, M_G, dwarf class and offset from the Pecaut & Mamajek dwarf sequence under (i) the applied `ebv_star` and (ii) the full line-of-sight column (and (iii) the 3D-map value when available).

## Vetting: `extinction`
Compares `ebv_star` with the 3D-map value at the star's distance and with the line-of-sight total (and SFD with Planck when both maps exist).
- `fail`: any estimate above `vetting.extinction_heavy_dust_ebv` (the correction dominates), applied vs 3D difference above `max(extinction_disagreement_mag, 2 sigma_3D)`, applied value above the whole column by more than the tolerance, or two 2D maps disagreeing.
- `not_evaluated`: fewer than two estimates, or the column exceeds the applied value by more than `extinction_los_excess_summary_ebv` with no 3D map saying where the dust lies (intrinsic colour and luminosity uncertain).
- `pass`: otherwise.
Other checks in `vet/underluminous_checks.py` never pass on defaults: for live stars an input from a vetting query counts only if `vetting_archive.query_status[<galex|gaia_meta|irsa_dust>] == "ok"`; for fixtures and hand-built dicts an explicit value counts; missing inputs give `not_evaluated` naming the missing query. Live `teff_spec_k` is a copy of the Gaia Teff, so `parameter_agreement` is `not_evaluated` for live stars.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.underluminous`:
- `delta_m_sigma_threshold: 4.5`, `min_delta_mag: 0.35`, `max_ruwe: 1.4`
Under `vetting`: `extinction_disagreement_mag: 0.20`, `extinction_heavy_dust_ebv: 0.30` (was hard-coded), `extinction_los_excess_summary_ebv: 0.05`, `sfd_to_sf11_factor: 0.86`, `bayestar19_ebv_per_unit: 0.884`.

## Known Impostors
- Heavily dust-extinguished stars, metal-poor halo/thick-disk subdwarfs, and high-RUWE unresolved binaries.

## Benchmark Objects
- Bidirectional switch-off validation: disabling extinction correction flags reddened stars behind heavy dust; re-enabling explains them.

## Failure Modes
- Disagreeing 3D dust maps or biased `Teff` can mimic faintness.
- A star whose applied reddening (e.g. Gaia GSP-Phot) is much lower than the line-of-sight column has an uncertain intrinsic colour and luminosity until a 3D map places the dust.

## Lessons Applied
- Require two independent parameter sets and agreement between independent extinction estimates.
