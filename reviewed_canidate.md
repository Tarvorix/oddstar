# Reviewed Surviving Candidate: `Gaia DR3 5628404318169472`

Out of **1,800 stars** searched across the local thick-disk/halo (`300`), all-sky Galactic Halo (`1,000`), and Galactic Bulge (`500`), **`Gaia DR3 5628404318169472`** is the **only surviving candidate**, retained at **Tier C** after three rounds of independent adversarial review ([adversarial-2026-11-17.md](review/adversarial-2026-11-17.md)) and Lesson [`L003`](lessons/L003_adversarial_audit_and_missing_modality_routing.md) pipeline hardening. *(All other flagged stars—including `Gaia DR3 1111149103659931264` / `G 251-44`, confirmed as an `sdK7` subdwarf—have been resolved as **Explained** or **Rejected** in [candidates.md](docs/candidates.md).)*

---

## 1. Quick Summary

| Field | Value |
| --- | --- |
| **Gaia DR3 Source ID** | `Gaia DR3 5628404318169472` |
| **Cross-Identifications** | `AllWISE J025411.56+055258.2`, `VarWISE 2540580983592406`, `AAVSO VSX ROT` (`oid=2850403`), `TIC 387517877`, `PS1 115050435482009905` |
| **Sky Position (J2016.0)** | $\text{RA} = 43.54827286^\circ$, $\text{Dec} = +5.88271224^\circ$ ($b = -45.609^\circ$) |
| **Distance & Tangential Velocity** | $\varpi = 4.3149 \pm 0.0162\text{ mas}$ ($d = 231.76\text{ pc}$), $\mu = (+37.914, -97.175)\text{ mas yr}^{-1}$, $v_{\text{tan}} = 114.60\text{ km s}^{-1}$ |
| **Stellar Parameters (Gaia DR3)** | $T_{\text{eff}} = 4509.8\text{ K}$ (GSP-Phot, K4V), $\log g = 4.55$, $[\text{M/H}] = +0.036\text{ dex}$, $M_G = 6.532$ (`6.605` unreddened), $L = 0.1817\,L_\odot$ (FLAME `0.223 L_sun`) |
| **Pipeline Tier & Disposition** | **Tier C** (`ir_excess`, `ir_timedomain`, `generic_outlier` flagged; **19 passed**, **0 failed**, **4 `not_evaluated` open holes**) |
| **Full Dossier & Review Evidence** | [dossier.md](dossiers/5628404318169472/dossier.md) · [adversarial-2026-11-17.md](review/adversarial-2026-11-17.md) · [R1](review/records/REV-20261117-R1-5628404318169472.json) / [R2](review/records/REV-20261117-R2-5628404318169472.json) / [R3](review/records/REV-20261117-R3-5628404318169472.json) |

---

## 2. Multi-Archive Measurements

### Cryogenic Mid-IR Excess (AllWISE 2010 & Single-Exposure L1b)
- **Optical & NIR baseline**: Gaia $G = 13.430$, $BP = 14.041$, $RP = 12.687$; Pan-STARRS DR2 (`67` detections) $g = 14.258$, $r = 13.476$ ($r_{\text{PSF}} - r_{\text{Kron}} = -0.001\text{ mag}$, clean point source); 2MASS (`AAA`) $J = 11.762 \pm 0.022$, $H = 11.143 \pm 0.027$, $K_s = 10.971 \pm 0.019$.
- **AllWISE (`AAAA`, `cc_flags=0000`, `ext_flg=0`)**: $W1 = 10.665 \pm 0.023$, $W2 = 10.276 \pm 0.021$, $W3 = 8.701 \pm 0.026$, $W4 = 7.290 \pm 0.102\text{ mag}$.
- **Photosphere-fixed excess ($T_{\text{eff}} = 4510\text{ K}$ normalized at $K_s = 10.971$)**:
  - $W3_{\text{obs}} = 10.48 \pm 0.25\text{ mJy}$ vs. $W3_{\text{pred}} = 1.89\text{ mJy}$ (**$5.54\times$**, $\Delta W3 = +8.59\text{ mJy}$).
  - $W4_{\text{obs}} = 10.15 \pm 0.95\text{ mJy}$ vs. $W4_{\text{pred}} = 0.55\text{ mJy}$ (**$18.31\times$**, $\Delta W4 = +9.59\text{ mJy}$).
  - Smaller baseline excess is already present at $W1$ ($\Delta W1 \approx +1.29\text{ mJy}$) and $W2$ ($\Delta W2 \approx +3.81\text{ mJy}$).
- **Cryogenic repeatability & localization**: **19/19** WISE single-exposure frames (`14` in Jan 2010, `5` in Aug 2010) consistently recover $W3$ ($\text{SNR } 18.2–23.8$) and $W4$ ($\text{SNR } 2.1–5.8$, `w4rchi2 = 1.02`, `w4m = 17`). AllWISE J2010.57 position agrees to `0.059"–0.064"` with Gaia proper-motion propagation, and direct atlas pixel Gaussian fits yield $W3/W4$ offsets of `0.34"` / `0.44"` (within `~0.6"` W4 blank-patch error).

### Multi-Year Mid-IR Brightening & Astrometric Shear Test (NEOWISE-R 2014–2024)
- **Light curve (`21` semi-annual visits, `283` frames)**: From `2014–2016` quiescence ($W1 \approx 10.630$, $W2 \approx 10.264$) to the `2017–2021` peak ($W1 = 10.406 \pm 0.010$ at MJD `58493.7`, $W2 = 9.871 \pm 0.017$ at MJD `59224.9`), the star brightened by `~0.24 mag` in $W1$ and `~0.41 mag` in $W2$ (`chi2_w2 = 41.19`, `slope_snr = 10.77`, isolated flare color $(W1 - W2)_{\text{flare}} \approx +1.06\text{ mag}$) before decaying back toward baseline across `2022–2024` ($W2 = 10.162–10.210$). Meanwhile field star `J025409.21+055101.5` (`121"` away) stayed within `0.045 mag` full range ([visit CSV](review/evidence/2026-11-17/5628404318169472-neowise-visits.csv), [control CSV](review/evidence/2026-11-17/J025409.21+055101.5-neowise-control.csv)).
- **Proper-motion tracking during flare**: Fitting linear proper motion to all `21` NEOWISE-R visit centroids gives $(+39.10, -93.38)\text{ mas yr}^{-1}$ (matching Gaia's `+37.91, -97.17`), with median residual **`0.053"`** across all visits (`0.051"` during `2018–2021` flare vs. `0.054"` at baseline; `0.026"` differential shift in the 277-frame track check), while stationary 2010 separation reaches `1.003"` at peak and `1.407"` by `2024.03`. Displaced (`> 0.2"–0.3"`) stationary background AGN blends are disfavored.

### Optical Coverage (Gaia DR3, ZTF, ASAS-SN, & TESS)
- **Gaia DR3 `vari_summary`**: `26` transits (`913.55 d`) span `0.0276 mag` peak-to-peak (`std = 0.0075 mag`, `SOLAR_LIKE` score `0.185`; AAVSO VSX lists `ROT`).
- **ZTF & ASAS-SN** ([ZTF CSV](review/evidence/2026-11-17/5628404318169472-ztf-annual.csv), [ASAS-SN CSV](review/evidence/2026-11-17/asassn-seasonal.csv)): IRSA ZTF yearly medians (`2018–2024`, `1,325` rows) span `0.060 mag` (`g`) and `0.055 mag` (`r`), and differential ASAS-SN photometry (`1,772` coadded epochs) spans `0.065 mag` in `g`, ruling out a sustained `~0.4 mag` optical counterpart to the `2017–2021` $W2$ brightening.
- **TESS (`TIC 387517877`)**: Covered in Sectors `4, 31, 42–44, 70, 71` (`0.83–1.94%` central 90% QLP range; preliminary TESScut raw correlations `r = 0.43–0.81`).

---

## 3. Vetting Status & Critical Caveats

- **Passed checks (19)**: `xmatch_integrity` (`0.064"` propagated), `wise_quality` (`0000`, `ext_flg=0`, `AAAA`), `centroid_offset` (`0.06"` propagated), `neighbours` (`0` inside `12"`), `radio_counterpart` (`0` NVSS within `10"`), `deep_imaging` (PS1 DR2 PSF), `chance_alignment` (`P=0.0090`), `extragalactic_colours` (`W1-W2=0.39`, `W2-W3=1.57`), `cirrus` (`E(B-V)=0.030`), `single_frames` (`19/19`), `moving_object`, `youth` *(archival UV/X-ray/association only; age remains open)*, `evolved`, `vision_rubric`, `neighbour_signal`, `systematics_window`, `sso_crossing`, `known_variable`, `metallicity`.
- **Open holes (`not_evaluated` → Tier C cap; 4 checks)**:
  1. `second_instrument` — No archival Spitzer, AKARI, or IRAS `12–22 um` counterpart (`SPHEREx` stops at `5 um`).
  2. `binary` — Gaia DR3 faint-star RVS (`rv_method_used = 2`, `13` transits, `SNR = 11.16`, `rv_chisq_pvalue = None`) gives `3.78 +/- 2.13 km/s` uncertainty on the combined median (implied scatter `~6.13 km/s`); despite `RUWE = 1.078` and `NSS = 0`, close SB1 companions remain unconstrained.
  3. `pipeline_agreement` — QLP vs. TESScut light curves exist across 7 TESS sectors, needing matched detrending before ingestion.
  4. `second_survey` — IRSA ZTF CGI timed out in programmatic live vetting (mitigated offline via [adversarial-2026-11-17.md](review/adversarial-2026-11-17.md), leaving onset/short dips less constrained).
- **Physical modeling caveats from review**:
  - The automated `585.8 K` / `f_dust = 4.4%` single-blackbody fit underpredicts $W4$ by ~`2x` when normalized to $W3$; the flat $12–22\,\mu\text{m}$ $F_\nu$ spectrum ($\Delta W3/\Delta W4 \approx 0.90$) corresponds to `~305 K` (`~1.7%` conditional luminosity) across $W3–W4$ alone, requiring multi-temperature dust, silicate emission, or a subarcsecond blend to link `2010 W3/W4` and `2017–2021 W1/W2`.
  - `generic_outlier` shares WISE features with `ir_excess` and `ir_timedomain`, so its `1e-12` analytic floor is not an independent third detection.

---

## 4. What Would Settle It

1. **High-resolution `5–25 um` imaging & spectroscopy (e.g., JWST MIRI)** to verify `< 0.3"` spatial association with the K4V star, test whether the `2010` $W3/W4$ flux decayed post-`2021`, and check for silicate emission vs. an obscured subarcsecond background source.
2. **Optical spectroscopy & multi-epoch RV monitoring (`R >= 20,000`)** to measure H$\alpha$, Li I $\lambda 6708$, spectral subtype, and orbital binarity across phases.
3. **Subarcsecond AO optical/NIR imaging + matched TESS detrending** to bound close companions/faint red blends and close `pipeline_agreement`.
