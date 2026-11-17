# Channel: `underluminous`

## Physical Signature
Main-sequence stars appearing systematically fainter in extinction-corrected absolute magnitude than predicted by two independent stellar-parameter determinations.

## Inputs
- Gaia DR3/DR4 parallax, RUWE, `G` photometry, Gaia XP parameters, survey spectroscopic parameters (`teff_spec_k`, `logg`, `[Fe/H]`), and SFD/Planck/Bayestar dust maps.

## Algorithm
1. Build empirical `M_G(Teff, logg, [Fe/H])` locus from the parent sample.
2. Compute extinction-corrected `M_G,obs` and residual `Delta M = M_G,obs - M_G,pred` under both XP and spectroscopic parameters.
3. Flag stars faint under both parameter sources with clean astrometry (`RUWE < 1.4`, `parallax_snr > 10`).

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.underluminous`:
- `delta_m_sigma_threshold: 4.5`, `min_delta_mag: 0.35`, `max_ruwe: 1.4`

## Known Impostors
- Heavily dust-extinguished stars, metal-poor halo/thick-disk subdwarfs, and high-RUWE unresolved binaries.

## Benchmark Objects
- Bidirectional switch-off validation: disabling extinction correction flags reddened stars behind heavy dust; re-enabling explains them.

## Failure Modes
- Disagreeing 3D dust maps or biased `Teff` can mimic faintness.

## Lessons Applied
- Require two independent parameter sets and agreement between two extinction maps.
