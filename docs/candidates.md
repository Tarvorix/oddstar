# `oddstar` Living Candidates Registry & Research Log

This document is the living scientific log of all real-sky candidates discovered, vetted, and followed up by `oddstar`, together with an incrementally updated **Research Paper Draft** section. Every measurement, query, and vetting result below is tied to frozen pipeline thresholds (`config_hash=ce6e224fd31d`) and content-addressed raw archive caches (`data/raw/`).

---

## 1. Master Candidates List

| Source ID | Designation / Region | Coords (J2016.0 RA, Dec) | Dist / $v_{\text{tan}}$ | Flagged Channels | Tier | Status / Summary |
| --- | --- | --- | --- | --- | --- | --- |
| [`Gaia DR3 5628404318169472`](../dossiers/5628404318169472/dossier.md) | `AllWISE J025411.56+055258.2` *(Local Thick-Disk/Halo, Cetus)* | `43.54827°, +5.88271°` ($b = -45.61^\circ$) | `231.76 pc`, `114.60 km/s` | `ir_excess`, `ir_timedomain`, `generic_outlier` | **Tier C** | Isolated high-velocity K4V dwarf with 6× W3 / 19.6× W4 excess ($T_{\text{dust}} = 585.8\text{ K}$, $f_{\text{dust}} = 4.4\%$) and a real multi-year NEOWISE-R W1/W2 brightening (+0.40 mag in W2, 2017–2021) followed by fading (2022–2024). **19/22** vetting checks passed; 3 open holes (`second_instrument`, `pipeline_agreement`, `second_survey`). |
| [`Gaia DR3 1111149103659931264`](../dossiers/1111149103659931264/dossier.md) | `SIMBAD LM*` *(All-Sky Halo)* | `115.85563°, +72.80826°` ($b = +29.65^\circ$) | `58.80 pc`, `341.10 km/s` | `generic_outlier` | **Tier C** *(cool subdwarf lead)* | Nearby extreme-velocity metal-poor (`[Fe/H]=-1.09`) K/M subdwarf (`Teff=4100 K`, `M_G=8.47`). No mid-IR excess (`chi_w4=0.37`) and 9 flat NEOWISE visits (`chi2_w2=0.03`); flagged on cool optical color (`z_bp_rp=3.52`) and optical/spectral vetting holes (`pipeline_agreement`, `second_survey`, spectral). |
| [`Gaia DR3 3899891606837752576`](../dossiers/3899891606837752576/dossier.md) | All-Sky Halo ($[\text{Fe/H}] = -2.00$) | `179.57271°, +8.37526°` ($b = +67.42^\circ$) | `1305.90 pc`, `174.57 km/s` | `generic_outlier` | **Explained** | **Background galaxy blend (`deep_imaging = fail`, `3.59"` offset in MAST PS1 DR2)**; 22 flat NEOWISE visits (`chi2_w2=0.51`). |
| [`Gaia DR3 2427069874188580480`](../dossiers/2427069874188580480/dossier.md) | `SIMBAD PM*` *(All-Sky Halo, $[\text{Fe/H}] = -1.49$)* | `9.63321°, -8.31170°` ($b = -70.94^\circ$) | `117.40 pc`, `304.38 km/s` | `generic_outlier` | **Explained** | **Background source at `6.86"` in MAST PS1 DR2 (`deep_imaging = fail`)**; W4 is low-SNR (`ph_qual=AAAB`, no IR excess). |
| [`Gaia DR3 1642833413151302400`](../dossiers/1642833413151302400/dossier.md) | All-Sky Halo ($[\text{Fe/H}] = -0.84$) | `239.90169°, +65.83897°` ($b = +41.63^\circ$) | `904.87 pc`, `172.18 km/s` | `generic_outlier` | **Explained** | **Background source at `2.85"` in MAST PS1 DR2 (`deep_imaging = fail`)**; evolved luminosity ($M_G = 0.37$), `55/55` WISE single frames consistent. |
| `Gaia DR3 5743059538967112576`, `3962473193054146944` | All-Sky Halo | North Galactic Cap / Hydra | `1367–1709 pc`, `287–406 km/s` | `generic_outlier` | **Explained** | Automatically explained by `metallicity` (`[Fe/H] = -1.01, -1.69` metal-poor subdwarf locus). |
| `Gaia DR3 4753637578074800768`, `6517329683012486912` | All-Sky Halo | `RA 44.4° / 339.7°`, `Dec -47.5° / -48.9°` | `1156–1331 pc`, `[Fe/H] -1.13/-1.37` | `ir_excess` | **Rejected** | Failed `wise_quality` (`ph_qual = AAAU` W4 upper limit; `cc_flags = 00P0`, `ph_qual = AABU` diffraction halo + upper limit). |
| [`Gaia DR3 4120051454339127808`](../dossiers/4120051454339127808/dossier.md) | `SIMBAD LP*` *(Galactic Bulge)* | `264.67916°, -19.28223°` ($b = +6.40^\circ$) | `2940.7 pc` | `generic_outlier` | **Explained** | SIMBAD long-period variable (`LP*`) behind heavy Galactic Bulge reddening (`E(B-V) = 0.803`, `extinction = fail`). |
| `Gaia DR3 4050034592591468544` | Galactic Bulge | Galactic Center region | `2.2–10.6 kpc` sample | `ir_excess` | **Explained** | Crowded bulge field: `neighbours = fail` (**22 Gaia DR3 neighbours inside 12" W4 beam**, `220%` blend flux). |
| `4050404539548308992`, `4089516951284506240`, `4036202422728871680`, `4107777331159599488`, `4050084306784858624`, `4059323576110862592` | Galactic Bulge | `RA 258°–271°`, `Dec -43.5°..-14.2°` | `2.2–10.6 kpc` | `ir_excess` *(+`generic_outlier` on 2)* | **Rejected** | Failed `wise_quality` on halo/diffraction contamination flags (`cc_flags = ddHH, 00hD, hhp0, hhHH, dhh0, hhdH`, `w4rchi2` up to `74.74`) and/or `AAAB`. |
| `Gaia DR3 8362885440918400` | `AllWISE J025505.03+074731.4` *(Local Halo/Thick-Disk)* | `43.77112°, +7.79202°` ($b = -44.00^\circ$) | `217.41 pc`, `84.80 km/s` | `ir_excess`, `generic_outlier` | **Rejected** | Rejected on `wise_quality` (`ext_flg = 3`, `ph_qual = AAAB`). |

---

## 2. Survey Search Runs Completed

| Run ID | Region & Sample Selection | Channels Run | $N_{\text{searched}}$ | Flagged | Surviving Candidates & 95% CL Upper Limit (`ir_excess`) |
| --- | --- | --- | ---: | ---: | --- |
| `live-halo-fusion` (`live-halo-ir/ul/out`) | Local Thick-Disk/Halo box (`Cetus/Taurus`, `RA 39.5–53.0°`, `Dec 0.2–9.6°`, `b -49..-38°`, `d 65–498 pc`, `v_tan > 65 km/s`) | `ir_excess`, `underluminous`, `generic_outlier` | 300 | 2 | **1 Tier C (`5628404318169472`)**, 1 Rejected; UL $\le 3.123\%$ (50%) / $\le 1.745\%$ (90% sens.) |
| `live-allsky-halo-fusion` (`live-allsky-halo-ir/ul/out`) | **True All-Sky Galactic Halo** (`RA 0.2–359.5°`, `Dec -89.0..+87.5°`, `759` HEALPix cells, `|b| > 20°`, `[Fe/H] -3.54..-0.80`, `v_tan > 150 km/s`, `d 59–1992 pc`, `ORDER BY random_index`) | `ir_excess`, `underluminous`, `generic_outlier` | 1,000 | 8 (`2` IR, `0` UL, `6` outlier) | **0 `ir_excess` survivors** (`2` Rejected); 5 `generic_outlier` **Explained** after live backup vetting (`deep_imaging` / `metallicity`), **1** subdwarf lead (`1111149103659931264`). **UL $\le 0.332\%$ at 90% sensitivity ($N_{\text{eff}} = 900$)** / $\le 0.597\%$ at 50%. |
| `live-bulge-fusion` (`live-bulge-ir/ul/out`) | **True Galactic Bulge** (`RA 249.3–282.0°`, `Dec -43.5..-14.2°`, `l` within $\pm 12^\circ$, `|b| 3–10°`, `41` HEALPix cells, `d 2.23–10.56 kpc`, `ORDER BY random_index`) | `ir_excess`, `underluminous`, `generic_outlier` | 500 | 8 (`7` IR, `0` UL, `3` outlier) | **0 survivors** after live inspection (`6` Rejected `wise_quality` artifacts, `1` Explained `22` neighbours, `1` Explained `LP*` + `E(B-V)=0.803`). **UL $\le 0.664\%$ at 90% sensitivity ($N_{\text{eff}} = 450$)** / $\le 1.191\%$ at 50%. |

---

## 3. Detailed Candidate Dossier: `Gaia DR3 5628404318169472`

### 3.1 Identification, Kinematics, & Stellar Parameters
- **Identifiers**: `Gaia DR3 5628404318169472`, `AllWISE J025411.56+055258.2`, `PS1 objID 115050435482009905`
- **Coordinates (J2016.0)**: $\text{RA} = 43.54827286^\circ$, $\text{Dec} = +5.88271224^\circ$ ($b = -45.609^\circ$)
- **Distance & Kinematics (ESA Gaia DR3)**:
  - Parallax: $\varpi = 4.3149 \pm 0.0162\text{ mas}$ ($\varpi/\sigma_\varpi = 267.14$, $d = 231.76\text{ pc}$)
  - Proper motion: $\mu_{\alpha*} = +37.914\text{ mas yr}^{-1}$, $\mu_\delta = -97.175\text{ mas yr}^{-1}$, $v_{\text{tan}} = 114.60\text{ km s}^{-1}$
  - Astrometric quality: $\text{RUWE} = 1.0783$, `non_single_star = 0`, `astrometric_excess_noise = 0.028 mas` ($\text{sig} = 0.56$), `ipd_frac_multi_peak = 0`
  - Radial velocity: $3.78 \pm 2.13\text{ km s}^{-1}$ (`rv_nb_transits = 13`)
- **Photospheric Parameters**:
  - $T_{\text{eff}} = 4509.8\text{ K}$, $\log g = 4.55$, $[\text{Fe/H}] = +0.036\text{ dex}$, $L = 0.1817\,L_\odot$, $M_G = 6.532\text{ mag}$, $E(B-V)_{\text{SFD}} = 0.030\text{ mag}$
  - Gaia DSC class probabilities: `star = 0.99983`, `galaxy = 5.18e-13`, `quasar = 1.04e-13`
- **SIMBAD Classification**: Not previously catalogued as a peculiar variable or extragalactic source (`simbad_type = Star`, `known_in_literature = False`; queried via CfA Harvard SIMBAD TAP mirror).

### 3.2 Multi-Band Photometry & SED Measurements
- **Gaia DR3 + 2MASS + AllWISE magnitudes**:
  - $G = 13.430$, $BP = 14.041$, $RP = 12.687$
  - $J = 11.762 \pm 0.022$, $H = 11.143 \pm 0.027$, $K_s = 10.971 \pm 0.019$ (`tmass_ph_qual = AAA`)
  - $W1 = 10.665 \pm 0.023$, $W2 = 10.276 \pm 0.021$, $W3 = 8.701 \pm 0.026$, $W4 = 7.290 \pm 0.102$ (`ph_qual = AAAA`, `cc_flags = 0000`, `ext_flg = 0`)
- **Pan-STARRS DR2 (MAST `nDetections = 67`)**:
  - $g = 14.258$, $r = 13.476$ ($r_{\text{Kron}} = 13.478$, $r_{\text{PSF}} - r_{\text{Kron}} = -0.001\text{ mag}$), $i = 13.125$, $z = 12.939$, $y = 12.810$
- **Channel 1 (`ir_excess`, `FAP = 1.0e-12`)**:
  - Photospheric fit: $T_{\text{eff,fit}} = 4095\text{ K}$, scale $= 23.35\text{ mJy}$
  - $W3_{\text{obs}} = 10.48\text{ mJy}$ vs $W3_{\text{pred}} = 1.752\text{ mJy}$ ($5.98\times$, $\chi_{W3} = 32.40$, empirical $\chi_{\text{emp},W3} = 34.30$)
  - $W4_{\text{obs}} = 10.15\text{ mJy}$ vs $W4_{\text{pred}} = 0.5174\text{ mJy}$ ($19.61\times$, $\chi_{W4} = 10.10$, empirical $\chi_{\text{emp},W4} = 9.05$)
  - Blackbody excess fit: $T_{\text{dust}} = 585.8\text{ K}$, $f_{\text{dust}} = 0.0440$ ($4.4\%$ of stellar bolometric luminosity), $R_{\text{bb}} = 0.096\text{ AU}$
- **Channel 6 (`generic_outlier`, `FAP = 1.0e-12`)**:
  - $\text{iForest quantile} = 1.000$ (`score = 0.738`), $\text{LOF quantile} = 1.000$ (`score = 51.52`), $\|z\| = 165.01$ (`max_abs_z = 160.1`, driven by `z_ir_var_chi2 = 160.1`, `z_chi_w4 = 38.84`, `z_ks_w3 = 8.76`).

### 3.3 Time-Domain Behavior (2010 WISE, 2014–2024 NEOWISE-R, & Gaia DR3 Optical)
1. **2010 WISE Cryogenic Single-Exposure Frames (`allsky_4band_p1bs_psd`)**:
   - Recovered **19/19 individual exposures** spanning January 2010 (`MJD 55219.8–55225.4`) and August 2010 (`MJD 55414.1–55414.3`).
   - All 19 frames show consistent detections in $W3$ ($8.539–8.862\text{ mag}$, $\text{SNR } 18.2–23.8$) and $W4$ ($6.661–7.781\text{ mag}$, $\text{SNR } 2.1–5.8$, `cc_flags = 0000` on 18/19 frames; `max/median = 1.52`, `w4rchi2 = 1.02`, `w4m = 17`).
2. **2014–2024 NEOWISE-R Light Curve (`neowiser_p1bs_psd`, 283 frames in 21 semi-annual visits)**:
   - Channel 2 (`ir_timedomain`, `FAP = 1.0e-12`) metrics: `chi2_w1 = 13.62`, `chi2_w2 = 41.19`, `slope_snr = 10.77`, `max_jump_snr = 8.16`, `color_change_mag = 0.197`, `n_visits = 21`.
   - **Visit-by-visit trajectory**:
     - **2014–2016** (`MJD 56686.9–57605.7`, 6 visits): $W1 \approx 10.54–10.64$, $W2 \approx 10.18–10.28$
     - **2017–2021 outburst/brightening** (`MJD 57772.2–59224.9`, 9 visits): smoothly rose to $W1 = 10.406 \pm 0.010$ (`MJD 58493.7`) and $W2 = 9.871 \pm 0.017$ (`MJD 59224.9`), a **~0.24 mag brightening in $W1$** and **~0.41 mag brightening in $W2$** (redder $W1 - W2 \approx 0.55\text{ mag}$ at peak).
     - **2022–2024 decay** (`MJD 59432.5–60320.4`, 6 visits): steadily faded back toward $W1 = 10.589–10.636$, $W2 = 10.162–10.210$.

| Epoch (MJD) | $N_{\text{frames}}$ | Median $W1 \pm \sigma$ (mag) | Median $W2 \pm \sigma$ (mag) |
| ---: | ---: | --- | --- |
| 56686.9 | 10 | $10.586 \pm 0.013$ | $10.218 \pm 0.037$ |
| 56876.7 | 13 | $10.541 \pm 0.017$ | $10.181 \pm 0.025$ |
| 57048.0 | 17 | $10.631 \pm 0.021$ | $10.262 \pm 0.025$ |
| 57238.9 | 21 | $10.624 \pm 0.018$ | $10.215 \pm 0.084$ |
| 57410.7 | 12 | $10.642 \pm 0.018$ | $10.282 \pm 0.019$ |
| 57605.7 | 14 | $10.630 \pm 0.017$ | $10.264 \pm 0.017$ |
| 57772.2 | 9 | $10.500 \pm 0.017$ | $10.011 \pm 0.019$ |
| 57969.8 | 13 | $10.453 \pm 0.019$ | $9.967 \pm 0.018$ |
| 58132.1 | 11 | $10.454 \pm 0.010$ | $9.960 \pm 0.012$ |
| 58336.8 | 20 | $10.445 \pm 0.014$ | $9.959 \pm 0.078$ |
| 58493.7 | 12 | $10.406 \pm 0.010$ | $9.874 \pm 0.025$ |
| 58701.1 | 11 | $10.444 \pm 0.018$ | $9.945 \pm 0.035$ |
| 58860.7 | 10 | $10.453 \pm 0.026$ | $9.896 \pm 0.011$ |
| 59065.3 | 14 | $10.448 \pm 0.010$ | $9.941 \pm 0.022$ |
| 59224.9 | 11 | $10.421 \pm 0.012$ | $9.871 \pm 0.017$ |
| 59432.5 | 11 | $10.519 \pm 0.009$ | $10.067 \pm 0.026$ |
| 59589.5 | 21 | $10.488 \pm 0.010$ | $9.981 \pm 0.020$ |
| 59796.6 | 12 | $10.552 \pm 0.016$ | $10.085 \pm 0.016$ |
| 59956.2 | 13 | $10.567 \pm 0.023$ | $10.129 \pm 0.013$ |
| 60160.7 | 12 | $10.636 \pm 0.015$ | $10.210 \pm 0.014$ |
| 60320.4 | 10 | $10.589 \pm 0.014$ | $10.162 \pm 0.023$ |

3. **Optical Variability Context (`gaiadr3.vari_summary` & `vari_classifier_result`)**:
   - **26 transits** spanning `913.55 days` (`G = 13.4216–13.4492`, peak-to-peak `0.0276 mag`, `std_dev = 0.0075 mag`, `BP std_dev = 0.0105 mag`, `RP std_dev = 0.0081 mag`), classified as **`SOLAR_LIKE`** (low-amplitude spot modulation, not an eclipsing binary, LPV, or AGN).
   - Thus, while Gaia DR3 optical variability stays within `~2.5%` peak-to-peak during 2014–2017, the mid-infrared $W2$ emission brightened by `~46%` over 2017–2021.

### 3.4 Complete Vetting Ledger (19 Pass, 0 Fail, 3 Open Holes)
- **Passed (19)**: `xmatch_integrity` (`0.064"`), `wise_quality` (`0000`, `ext_flg=0`, `AAAA`, `w4rchi2=1.02`, `w4m=17`), `centroid_offset` (`0.58"`), `neighbours` (`0` Gaia DR3 neighbours in `12"` cone), `radio_counterpart` (`0` NVSS sources within `10"` via VizieR CfA), `deep_imaging` (MAST PS1 DR2 clean isolated PSF, single-detection diffraction ghosts `nDetections=1` excluded), `chance_alignment` (`P=0.0090`), `extragalactic_colours` (`W1-W2=0.39`, `W2-W3=1.57`), `cirrus` (`E(B-V)=0.030`), `single_frames` (`19/19` consistent), `moving_object` (none), `youth` (`pass`), `evolved` (`M_G=6.53` main sequence), `binary` (`RUWE=1.08`, `NSS=0`, `RV scatter=2.13 km/s`), `vision_rubric`, `neighbour_signal`, `systematics_window`, `sso_crossing`, `known_variable`.
- **Open Holes (`not_evaluated` → Tier C cap)**:
  1. `second_instrument` (`akari_irc`, `akari_fis`, and Spitzer `slphotdr4` returned 0 detections at that depth)
  2. `pipeline_agreement` (no independent TESS QLP vs TESScut light-curve pair loaded yet)
  3. `second_survey` (`irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves` timed out at `12–25s`; ZTF/ASAS-SN optical coverage overlapping 2017–2024 is needed to confirm optical quiescence during the 2017–2021 IR flare).

### 3.5 What Was Built & Fixed During This Investigation
1. **Live Multi-Archive Ingest & Backup Mirrors** ([live_survey.py](file:///Users/kylebullock/oddstar/src/oddstar/io/live_survey.py)): Implemented `fetch_live_halo_bulge_sample`, `fetch_single_star_live`, and `enrich_star_live_for_vetting` with `RawCache` SHA-256 provenance and automatic fallback across ESA Gaia / ARI Heidelberg / Gaia@AIP, SIMBAD CfA / Strasbourg, IRSA Gator (`allwise_p3as_psd`, `allsky_4band_p1bs_psd`, `neowiser_p1bs_psd`, `akari_irc`, `akari_fis`, `slphotdr4`) / IRSA TAP / NOIRLab Data Lab, MAST PS1 DR2 (with `nDetections >= 2` ghost filtering), and VizieR CfA `asu-tsv`.
2. **Lesson `L002`** ([L002_missing_plate_data_high_pm_vanishing.md](file:///Users/kylebullock/oddstar/lessons/L002_missing_plate_data_high_pm_vanishing.md), [channel.py](file:///Users/kylebullock/oddstar/src/oddstar/channels/vanishing/channel.py)): Caught and fixed a bug exposed by `5628404318169472` (`v_tan = 114.6 km/s`) where missing archival plate coordinates (`vanishing_data={}`) defaulted to J2016 coordinates and produced a synthetic `6.53"` proper-motion offset.

---

## 4. Research Paper Section (Living Draft)

> **Working Title**: *A Reproducible Multi-Channel Search for Mid-Infrared and Time-Domain Stellar Outliers in High-Velocity Solar-Neighborhood Stars*

### Abstract *(Draft)*
We present initial live-sky results from `oddstar`, a reproducible, provenance-tracked pipeline designed to identify statistically significant multi-band stellar anomalies while enforcing explicit tri-state archival vetting (`pass`, `fail`, `not_evaluated`). Searching a kinematically selected sample of 300 high-tangential-velocity stars ($v_{\text{tan}} > 65\text{ km s}^{-1}$, $|b| > 12^\circ$, $\varpi/\sigma_\varpi > 15$, $\text{RUWE} < 1.4$) drawn from Gaia DR3 cross-matched to 2MASS and AllWISE across three channels (`ir_excess`, `underluminous`, and `generic_outlier`; 900 total trials), we find zero under-luminous candidates and two mid-infrared excess flags. Automated quality gating rejects `Gaia DR3 8362885440918400` (`ext_flg = 3`, `ph_qual = AAAB`). The surviving candidate, **`Gaia DR3 5628404318169472`** ($d = 231.76\text{ pc}$, $v_{\text{tan}} = 114.60\text{ km s}^{-1}$, $T_{\text{eff}} = 4510\text{ K}$, main-sequence K4V), exhibits a $5.98\times$ excess at $12\,\mu\text{m}$ ($\chi_{W3} = 32.4$) and $19.61\times$ excess at $22\,\mu\text{m}$ ($\chi_{W4} = 10.1$), well fit by warm dust at $T_{\text{dust}} = 585.8\text{ K}$ with fractional luminosity $f_{\text{dust}} = 0.044$. Follow-up querying of MAST Pan-STARRS DR2 (`67` detections, isolated PSF) and all 19 cryogenic WISE single-exposure $W4$ frames rules out diffraction spikes, frame glitches, and resolved background galaxies. Crucially, 283 NEOWISE-R $W1/W2$ exposures (2014–2024) reveal a coherent, multi-year mid-infrared brightening ($\Delta W2 \approx -0.41\text{ mag}$ peaking in 2019–2021) followed by a gradual decay toward the 2010–2016 baseline, triggering the `ir_timedomain` channel ($\chi^2_{W2} = 41.19$). Three vetting gates (`second_instrument`, `pipeline_agreement`, `second_survey`) remain open, keeping the target at Tier C pending optical light-curve confirmation over 2017–2024 and independent $5–25\,\mu\text{m}$ spectroscopy.

### 1. Sample Selection & Methodology
Public-archive queries are executed through a content-addressed cache (`RawCache`) enforcing frozen pre-registered selection cuts (`config_hash = ce6e224fd31d`). Targets satisfy $\varpi \in [2, 25]\text{ mas}$ ($40–500\text{ pc}$), $\varpi/\sigma_\varpi > 15$, $\text{RUWE} < 1.4$, $|b| > 12^\circ$, $G \in [7.5, 13.5]\text{ mag}$, $T_{\text{eff}} \in [3800, 7500]\text{ K}$, $v_{\text{tan}} = 4.74047\,\mu/\varpi > 65\text{ km s}^{-1}$, unique `AllWISE` and `2MASS` cross-matches (`number_of_neighbours = 1`, `number_of_mates = 0`), and measured photometric uncertainties in $JHK_s$ and $W1–W3$. Synthetic injection testing over this sample establishes 95% confidence upper limits on unvetted/Tier-C excesses in the surveyed footprint (HEALPix pixels `1320, 1321, 1383, 1384, 1448`) of $\le 1.745\%$ at 90% completeness ($N_{\text{eff}} = 270$) and $\le 3.123\%$ at 50% completeness ($N_{\text{eff}} = 150$).

### 2. Candidate `Gaia DR3 5628404318169472`
Combining photometry from Gaia DR3, 2MASS, AllWISE, Pan-STARRS DR2, and 21 binned visits of NEOWISE-R establishes that `Gaia DR3 5628404318169472` is an unblended, astrometrically single ($\text{RUWE} = 1.078$, `non_single_star = 0`), near-solar-metallicity ($[\text{Fe/H}] = +0.036$) K4V dwarf ($M_G = 6.532$) moving with thick-disk/halo-boundary tangential velocity ($v_{\text{tan}} = 114.60\text{ km s}^{-1}$). Its $3.4–22\,\mu\text{m}$ spectral energy distribution displays an unusually large warm infrared fractional luminosity ($f_{\text{dust}} = 4.4\%$, $T_{\text{dust}} = 586\text{ K}$) for a main-sequence star without youth indicators or giant-branch luminosity. Multi-epoch NEOWISE-R photometry rules out a static background-galaxy blend hypothesis on two independent grounds: (a) Pan-STARRS DR2 resolves no offset galaxy inside $12''$ ($r_{\text{PSF}} - r_{\text{Kron}} = -0.001\text{ mag}$), and (b) the mid-infrared $W1/W2$ flux brightened and faded coherently over a 6-year timescale (`slope_snr = 10.77`, `chi2_w2 = 41.19`).

### 3. Open Hypotheses & Required Observations
Under `oddstar`'s 12-hypothesis matrix, `catalog_mismatch`, `blend_neighbouring_star`, `background_galaxy`, `galactic_cirrus`, `solar_system_object`, `young_star_disk`, `stellar_companion`, `evolved_star_mass_loss`, and `known_variable_class` are all currently `disfavoured`, while **`collision_debris`** (warm transient dust generated by a giant planetary collision or tidally disrupted planetesimal in a ~0.1 AU orbit) is `consistent`. Three required next measurements will resolve the remaining `untested` hypotheses (`instrument_artifact`, `unexplained`):
1. **Co-eval optical light curve (2017–2024)** from ZTF/ASAS-SN/TESS (`second_survey`, `pipeline_agreement`) to determine whether the 2017–2021 $W1/W2$ brightening was accompanied by an optical dipping/dimming event (chromatic dust extinction) or occurred optically quiescence-free (analogous to face-on/non-transiting post-impact debris plumes).
2. **Optical spectroscopy & RV monitoring** to check for chromospheric/accretion H$\alpha$ emission, Li I $\lambda 6708$ absorption, or a hidden short-period companion consistent with the $2.13\text{ km s}^{-1}$ Gaia RV uncertainty.
3. **Second-epoch $5–25\,\mu\text{m}$ spectroscopy/photometry** (e.g., JWST MIRI) to test whether the 2010 $W3/W4$ excess tracks the post-2021 $W1/W2$ decay and to probe warm silicate/silica emission features characteristic of high-velocity rocky vapor-condensate debris.

---

## 5. Changelog
- **2026-11-17 (Run 3 — True All-Sky Halo & Galactic Bulge)**: Added `--region halo|bulge|local` and `--after-random-index` to [live_survey.py](file:///Users/kylebullock/oddstar/src/oddstar/io/live_survey.py) and [cli.py](file:///Users/kylebullock/oddstar/src/oddstar/cli.py). Searched **1,000** true metal-poor Halo stars across `759` all-sky HEALPix pixels (`[Fe/H] <= -0.8`, `v_tan > 150 km/s`) and **500** Galactic Bulge stars (`l` within $\pm 12^\circ$, `2.2–10.6 kpc`). Zero unvetted `ir_excess` or `underluminous` candidates survived in either sample (giving 95% upper limits of $\le 0.332\%$ across the all-sky Halo and $\le 0.664\%$ in the Bulge); live MAST PS1 DR2 + IRSA Gator + SIMBAD CfA inspection explained background blends (`3899891606837752576`, `2427069874188580480`, `1642833413151302400`) and Bulge dust/`LP*` (`4120051454339127808`), leaving only cool halo subdwarf `1111149103659931264` (no IR excess) alongside our primary candidate `Gaia DR3 5628404318169472`.
- **2026-11-17 (Runs 1–2)**: Initialized living candidates log; added 300-star Cetus/Taurus high-velocity live survey results (`live-halo-ir`, `live-halo-ul`, `live-halo-out`, `live-halo-fusion`); added `Gaia DR3 5628404318169472` (Tier C) and rejected impostor `Gaia DR3 8362885440918400`; updated with backup-archive verification (MAST PS1 DR2 `deep_imaging=pass`, WISE L1b `single_frames=pass`, NEOWISE-R 21-visit `ir_timedomain=True`, Gaia DR3 `vari_summary`).
