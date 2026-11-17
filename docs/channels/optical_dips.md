# Channel: `optical_dips`

## Physical Signature
Aperiodic, asymmetric, or unusually deep optical dimming events and multi-year secular fading inconsistent with planetary transits or regular eclipsing binaries.

## Inputs
- Two independent TESS extractions (public QLP / SPOC and pipeline TESScut aperture photometry), Kepler/K2 light curves, ZTF, and ASAS-SN.

## Algorithm
1. Detrend light curves with an iterative upper-clipped sliding median filter preserving multi-day dips.
2. Compute Box-Least-Squares (BLS) periodogram to veto strictly periodic eclipsing binaries.
3. Run multi-duration boxcar matched filter for dip SNR, depth, ingress/egress asymmetry, and long-baseline ZTF linear fading slope.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.optical_dips`:
- `dip_snr_threshold: 6.0`, `min_dip_depth_frac: 0.03`, `min_duration_hours: 6.0`
- `periodic_bls_veto_snr: 8.0`, `secular_fade_snr_threshold: 4.5`

## Known Impostors
- Kepler/TESS eclipsing binaries, blended background EBs in 21" TESS pixels, scattered-light/momentum-dump artifacts, and detrending-only glitches.

## Benchmark Objects
- `KIC 8462852`, `EPIC 204278916`, `ASASSN-21qj`.
- Open case workflow target: `HD 139139`.

## Failure Modes
- Large 21" TESS pixels blend nearby EBs; `pipeline_agreement` and `neighbour_signal` centroid checks are mandatory.

## Lessons Applied
- Require dual-extraction agreement (`QLP` and `TESScut`) before vetting passes. Since `L006`, live vetting fetches QLP `SAP_FLUX` and an independent TESScut aperture extraction per sector, and `pipeline_agreement` applies the identical `detrend_window_days` running median to both at matched cadences before comparing `min_duration_hours`-binned dip depths.
- `second_survey` (since `L006`) uses live ZTF (retried, shrinking cone, yearly paging, PM-track matched) and ASAS-SN Sky Patrol seasonal photometry (`pip install oddstar[surveys]`).
