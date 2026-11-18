# Configuration Changelog

All threshold, sample-cut, and prompt changes must be logged here under Ground Rule R4 before running validation or search.

## [1.9.0] - 2026-11-18

### Vanishing vetting: the survey's comparison is repeated at the place where the first plates show the source
Run `v180-halo-vanishing` (config 1.8.0) left five plate sources that both first-epoch plate images show and that
have no counterpart at the catalogue position. All five were looked at on 2026-11-18, and each is an ordinary
object that never left: three stars and two extended objects (galaxies by their Legacy Survey profiles), on every
plate at one place (two of the stars drift at their Gaia proper motions of 22 and 42 mas/yr). USNO-B1.0 gives each
a proper motion of 92 to 553 mas/yr that the plates do not show: detections of different objects joined into one
track. The catalogue position for epoch 2000, where the survey looks for a Pan-STARRS counterpart, is then 4 to 25
arcsec from the object, and the mean-epoch plate position, where it looks for the paths of Gaia stars, 2 to 7
arcsec. The survey finds nothing within its 3 arcsec, and the later plates, measured where the false motion puts
the source, are empty. The rule below was prompted by those five sources, which is said here plainly. It adds no
limit: it is the survey's own comparison, with the survey's own `match_radius_arcsec`, made at another place.

- **No new keys.** The version changes because `proper_motion` gives different verdicts.
- **`proper_motion` for a plate source of the live survey without a Gaia counterpart.** `vet --live` now also reads
  the Gaia stars, moved back to the first plate's date with their Gaia proper motions, and the Pan-STARRS objects
  (at least `plate_ps1_min_detections` detections) around the position the first-epoch plates were measured at
  (`vanishing_data.plate_first_place`, status `plate_first_place`). The check fails ("the object the plates
  recorded is still there") when one of them is within `match_radius_arcsec` of the mean sky place of the
  first-epoch images that show the source. It passes only when that comparison and the moved-star test of 1.8.0
  could both be made and found nothing; when the catalogues at that place could not be read it is `not_evaluated`.
  A source that no first-epoch image shows has no measured place, and the comparison does not apply to it.
- **What the rule costs.** A source that did vanish is wrongly explained when an unrelated object lies within the
  match radius of its place. The survey's comparison at the catalogue position has the same cost. Both rates are
  measured on control positions before this version is applied to a flagged source: see the control numbers below.
- **Not changed.** The search, its flags and `plate_defect`. The later-plate measurement still uses the catalogue
  motion; it is reported, and no check reads it except through the track test of 1.8.0, which works from the peaks'
  sky places.

## [1.8.0] - 2026-11-18

### Vanishing vetting: the plate images decide whether a catalogue entry is a source
The halo run `v160-halo-vanishing` flagged 127 of 2,859 plate sources. Their catalogue rows were looked at on
2026-11-18: 83 of the 127 are listed on two plates only, and three in four of all two-plate entries of the fields
are flagged, so most flags are entries of the USNO-B1.0 catalogue that are not sources. `plate_defect` could not
say so: it counted the plates the catalogue lists. This version makes it read the digitized plates
(`oddstar.io.plate_images`, measurement `dss_plate_v2`). The rule was prompted by those 127 flags, which is said
here plainly. Its limits were fixed on control positions only, never on a flagged source: 203 ordinary plate
sources with a Pan-STARRS counterpart, spread over plate magnitude, and 69 blank positions in the same fields
(1,360 cutouts, 2026-11-18).

- **New keys** under `channels.vanishing`: `plate_image_min_snr: 5.0`, `plate_image_max_offset_arcsec: 2.0`,
  `plate_image_first_epoch_before: 1960.0`, `plate_image_judge_max_r_mag: 19.5`, `plate_moved_radius_arcsec: 10.0`,
  `plate_moved_max_dmag: 1.5`, `plate_track_max_residual_arcsec: 1.0`, `plate_track_min_images: 2`.
- **`plate_defect` for a plate source of the live survey.** An image shows the source when a local maximum with
  aperture signal-to-noise of at least 5.0 lies within 2.0 arcsec of the position, after the image has been
  registered on the Gaia stars it holds. An image that could not be registered can show the source but not its
  absence. Pass: both first-epoch images show it. Fail (tier Rejected): two first-epoch images could say, neither
  shows it, and the catalogue's red magnitude is 19.5 or brighter. Not evaluated: the images disagree, fewer than
  two could say, or the source is fainter than 19.5.
- **Measured on the controls with these limits.** Blank positions (69): 57 called absent, 12 without a verdict, none
  shown on either image; no blank position gave a detection at any limit tried (signal-to-noise 4 to 6, offset up
  to 3.0 arcsec). Ordinary sources brighter than 19.5 (190): 148 pass, 22 seen on one image only, 18 with fewer than
  two images that could say, 2 called absent (1.1%). By magnitude, pass / one image / absent / too few images:
  under 16: 36/3/0/1; 16 to 17: 36/1/0/3; 17 to 18: 30/3/1/6; 18 to 19: 29/8/0/3; 19 to 19.5: 17/7/1/5. From 19.5 to
  20 (10 sources) 3 pass and 1 is called absent, and of 3 sources fainter than 20, 2 are called absent: hence the
  magnitude limit for a fail.
- **Why registration.** The sky solution delivered with a cutout put 38 of 200 first-epoch red images 2 to 4 arcsec
  from the stars, consistently across a plate. Without registration one bright ordinary star in ten was missed.
- **`proper_motion` for a plate source without a Gaia counterpart.** A star too faint for Gaia has no Gaia path to
  test, but it can have moved several arcsec since the 1950s. The check now also reads the Pan-STARRS objects
  between the match radius and 10 arcsec of the plate position whose r magnitude is within 1.5 mag of the plate
  magnitude (60 known moving stars: r minus plate magnitude has median +0.19, 5th and 95th percentile -0.65 and
  +0.72). Such an object is another star when a first-epoch plate shows it where it is now (true of none of the 60
  moving stars). Otherwise the later plates decide, with new keys `plate_track_max_residual_arcsec: 1.0` and
  `plate_track_min_images: 2`: where the first-epoch images show the source, where Pan-STARRS has the object
  (epoch 2012), and the date of each later plate give the place a star on its way would have on that plate. Fail
  ("the same star, moved", tier Explained) when two later images show the source within 1.0 arcsec of that place
  and fewer than two show it still at its first-epoch place; the object is ruled out when two later images show
  the source still in place; `not_evaluated` otherwise, and when the two places are less than the match radius
  apart.
- **Measured on controls for the track rule** (60 Gaia stars moving 57 to 165 mas/yr with a USNO-B entry, 100 still
  USNO-B sources with an unrelated Pan-STARRS neighbour 3 to 10 arcsec away; none a flagged source). Moving stars
  with a first-epoch position at least 3 arcsec from the object (57): 43 "moved", 0 "stayed", 14 open (8 without a
  usable later image). Still stars with a first-epoch position at least 3 arcsec from the neighbour (79): 73
  "stayed", 0 "moved", 6 open. The same counts hold for limits of 0.75 and 1.0 arcsec.
- **A rule that was tried and dropped.** "The object is absent from the first-epoch plates at its present place, so
  it moved" called 10 of 100 still stars moved and found 17 of the 60 moving stars. It is not used.

### Youth: a star moving at thick-disc or halo speed is not young
The one infrared candidate left after 1.7.0, Gaia DR3 5628404318169472, stood at `youth: not_evaluated` because the
ROSAT limit at its distance cannot reach the X-ray threshold. Its proper motion and parallax give 115 km/s across
the sky. The check had no rule for motion; this version adds one, prompted by that star and said so here. The
limit is not fitted to it: it is the usual kinematic boundary of the thin disc.

- **New keys** `vetting.youth_old_vtan_min_kms: 85.0` and `vetting.youth_old_vtan_min_parallax_snr: 10.0`. Stars
  younger than about a billion years share the motion of the gas they formed from to a few tens of km/s; a
  tangential velocity of 85 km/s or more is thick-disc or halo motion. Runaway stars are the exception and are
  hot, massive stars, which this channel's sample cuts exclude.
- **Rule.** Kinematics is a third subtest beside UV and X-ray, and can only say "not young". `youth` passes when at
  least two of the three say "not young" and none says "young" (before: both UV and X-ray had to). A slow star
  gives no verdict from its motion.

### Not changed
- Every channel's search thresholds, so no search result changes. Vetting results of vanishing flags and of
  5628404318169472 made under 1.7.0 are redone under this version before any tier is quoted.
- Not covered yet: `sso_plate`, `blend_split`, `known_variable` for plate sources, which still have no live data.

## [1.7.0] - 2026-11-18

### Infrared-change vetting: saturated NEOWISE photometry and luminous giants
On 2026-11-18 every infrared-change (`ir_timedomain`) flag that vetting had not explained was looked at again
under the corrected NEOWISE systematics test (runs `v160-local-td`, `v160-bulge-td`, and the halo run). Two bulge
stars stood at Tier B and three stars at Tier C. Gaia DR3 gives four of the five a giant's surface gravity, and
four are far brighter than the NEOWISE saturation limits. The vetting of this channel had no rule for either
case: it ran the six light-curve checks only. This version adds both rules. They were written after looking at
those five stars, which is said here plainly; their limits are not fitted to them but taken from the NEOWISE
documentation and from the position of the main sequence.

- **Saturation** (`systematics_window`, NEOWISE part). New keys `channels.ir_timedomain.neowise_saturation_w1_mag: 8.0`
  and `neowise_saturation_w2_mag: 7.0`: "The characteristic transition from non-saturated to saturated detections
  occurs near W1=8.0 mag and W2=7.0 mags" (NEOWISE Explanatory Supplement, II.1.c.iv). Brighter than that the
  profile-fit fluxes are overestimated, by more than 1 mag near 4th magnitude, and the bias changed with the
  temperature of the telescope and detectors (II.1.c.iv.a). The quoted errors do not hold it, so a saturated star
  passes the channel's chi-square, slope and jump thresholds without varying. Rule: when the star's median
  magnitude is brighter than the limit in every band that carries its variation (W1: its chi-square; W2: its
  chi-square, the slope, the jump), the check fails as an artifact (tier Rejected). A star with an unsaturated
  band that varies is tested as before. Such stars also have no comparison stars of their brightness, so until now
  they stayed `not_evaluated` at Tier C for good.
- **Luminous giants.** An `ir_timedomain` flag now also gets the existing `evolved` check. A giant above the
  horizontal branch varies in the infrared by tenths of a magnitude over months to years (long-period and
  semi-regular variables); `evolved` failing gives tier Explained.
- **`vetting.evolved_teff_max_k`: 5200 -> 8000.** `evolved` fails for M_G < `evolved_mg_bright_limit` (0.5) and
  Teff below this limit. At M_G < 0.5 the main sequence is hotter than about 9,500 K, so a star that bright and
  cooler than 8,000 K has left it whatever its exact temperature. The old limit let a star pass that Gaia GSP-Phot
  puts at 5,600 K behind 4.3 mag of extinction with log g = 1.5 (Gaia DR3 4144479235457826816): the photometric
  temperature of a reddened giant is not reliable to 400 K. The check also applies to `ir_excess` and
  `generic_outlier` flags: a star of M_G < 0.5 between 5,200 and 8,000 K that was passed as "main sequence" there
  now fails as a luminous giant. The note of a passing star says "not a luminous giant" instead of "main sequence"
  (a clump giant or subgiant is neither).
- Not changed: every channel's search thresholds, so no search result changes. Vetting results of infrared-change
  flags made under 1.6.0 are redone under this version before any tier is quoted.
- Not covered yet: giants below the limit (red clump, subgiants) and crowding in the bulge fields have no rule of
  their own.

## Radio planner: ON and OFF files are of the same product - 2026-11-18, configuration unchanged (1.6.0)

The first Parkes scan set that was fetched (the calibrator field 1934-638, beam 3, of 2018-04-17) could not be
searched: its ON file was a `.0000.hires.hdf` product and its OFF files `.0000.h5` products, which differ in
channel width. The planner (`bl_opendata.plan_star_scans`) matched OFF files by frequency only. It now pairs files
of the same product only (`.hires.hdf` against `.hires.hdf`; `.h5` and `.fil` against each other) and gives one
plan for a scan held as two products: the product with a full comparison, then the one with more frequencies,
then `.h5`/`.fil`. A scan whose neighbours exist only as the other product is listed as incomplete. No search
result existed for such a set, and no stored result changes.

## Radio search version 4 (`bl_local_search_v4`): the allocation table is complete up to 4200 MHz - 2026-11-18, configuration unchanged (1.6.0)

The owner chose the rule of version 3 on 2026-11-18 ("List every allocated band"). The Parkes ultra-wideband
receiver, with which one star of the halo sample was observed, covers 704 to 4032 MHz, and the table was complete
only up to 2900 MHz. It is completed up to 4200 MHz by the same rule, before any scan above 2900 MHz was searched
(no result exists there, so no flag was looked at). Four bands are new: 2900 to 3100 MHz (maritime radionavigation
and radiolocation), 3100 to 3300 MHz (radiolocation), 3300 to 3550 MHz (radiolocation, amateur, mobile broadband)
and 3550 to 3700 MHz (citizens broadband radio service, radiolocation); 3700 to 4200 MHz was listed already. The
United States table gives no band between 2900 and 4200 MHz to passive services only, so in that range every
narrowband signal is marked as interference. The three passive bands below 2900 MHz are unchanged. Above 4200 MHz
the table is not complete.

The rule is taken from the United States table also for Parkes scans. Between 698 and 4200 MHz this changes
nothing: outside the three passive bands, which are protected worldwide, every frequency is listed whatever the
national plan gives it to.

`SEARCH_VERSION` goes to `bl_local_search_v4`. No result below 2900 MHz changes; the stored version 3 results are
scored again from the saved windows without a download.

## Radio search version 3 (`bl_local_search_v3`): the allocation table is complete from 698 to 2900 MHz - 2026-11-18, configuration unchanged (1.6.0)

### Why
The table `oddstar.io.bl_radio.RFI_ALLOCATIONS_MHZ` was extended twice in one night, each time after a search
showed flags in a band it did not list (entries below). A third search then did the same: the 2017 S-band scan of
HIP 11729 gave 64 flagged windows of 4,080, 49 of them outside the table: 35 at 1920 to 1930 MHz (cordless
telephones), 1 at 1820 MHz, 4 in the PCS downlink band, 6 in the satellite telemetry band 2200 to 2290 MHz and 3
in the 2400 MHz band of wireless networks. A table that grows where flags appear is fitted to the flags. It is
replaced by a rule that does not look at them.

### The rule
From 698 to 2900 MHz every band that the United States table of frequency allocations (47 CFR 2.106, federal and
non-federal columns) gives to a transmitting service is listed. In that range only the passive bands are not
listed: 1400 to 1427 MHz, 1660.5 to 1668.4 MHz and 2690 to 2700 MHz. Seventeen bands are new: 1390 to 1400,
1427 to 1435, 1668.4 to 1670, 1780 to 1850, 1915 to 1930, 1930 to 2000, 2000 to 2025, 2025 to 2110, 2110 to 2180,
2200 to 2290, 2290 to 2300, 2300 to 2320, 2345 to 2360, 2360 to 2400, 2400 to 2483.5, 2500 to 2690 and 2700 to
2900 MHz. A test checks that the listed bands and the passive bands tile 698 to 2900 MHz without gap or overlap.
Above 2900 MHz only the satellite downlink bands are listed; the table has to be completed there before scans of
those receivers are searched.

### What it costs
Between 698 and 2900 MHz a narrowband signal can pass `rfi_band` only in the three passive bands (44.9 MHz in
all). Everywhere else it is marked as interference whatever the ON and OFF scans show, also when a full cadence
shows it in every ON scan and in no OFF scan (the six lines at 825 MHz toward the NGC4168 field are such a case).
Only a re-observation could tell such a signal from a transmitter, and the archive holds none for these stars.
This is a choice between missing a real signal in a transmitter band and listing transmitters as leads; it is the
owner's to make and was put to him on 2026-11-18.

`SEARCH_VERSION` goes to `bl_local_search_v3`: no version 2 result is read again. The windows saved under
`_work` are scored again without a download; a scan set whose saved windows were deleted is fetched again.

## Radio: cadences and the bands from 698 to 960 MHz - 2026-11-18, configuration unchanged (1.6.0)

### Radio: ON scans that follow each other are searched as one cadence
The planner made one comparison per ON scan (that scan against the OFF scan before and after it), also where the
archive holds several ON scans of the same target in a row. The channel's own rule (`require_all_on_scans: true`)
could then not act: with one ON scan, "in every ON scan" is always true. Now ON scans of one target, telescope and
frequency set whose starts follow each other within 0.02 day become one plan (`comparison: "cadence"`); the signal
has to be in every ON scan and in no OFF scan. No threshold changes. A lone ON scan is handled as before and its
stored result keeps its key.

Basis: the two complete ON scans of the NGC4168 pointing that holds Gaia DR3 3920615854833203328 (HIP 59490;
Green Bank, 2021-07-31, 751 to 1126 MHz), searched one at a time, gave 255 and 1,327 flagged windows of 1,840 and
6,459 evaluated. The archive holds three ON scans of that pointing that start within 22 minutes.

The cut around a hit stays `live_window_channels` (2048) for a lone ON scan with OFF scans 5 to 6 minutes away.
For a cadence it is doubled until a signal at `max_abs_drift_hz_s` stays inside from the start of the reference
scan to the end of the farthest scan (4096 channels for three 5-minute ON scans 10.6 minutes apart); otherwise a
fast-drifting signal would leave the cut in the outer scans and be lost.

### Radio: four more transmitter bands, 698 to 960 MHz
Same table and same reasoning as in the entry below. All 906 flagged windows of those two searches that fell
outside the table lie between 751 and 960 MHz (173 and 654 in 751 to 939 MHz, 79 in 939 to 960 MHz); the other 676
were inside it (960 to 1126 MHz, aeronautical radionavigation). 698 to 960 MHz is allocated without a gap to
transmitters (47 CFR 2.106): 698 to 806 MHz mobile telephone and public safety, 806 to 902 MHz land mobile, mobile
telephone and air-ground, 902 to 928 MHz industrial, scientific and medical devices and radiolocation, 928 to
960 MHz paging, fixed and land mobile. These four bands are added. As before, the bands were chosen after looking
at the flags, so this field is not a test of them.

`SEARCH_VERSION` stays `bl_local_search_v2`: the added bands lie below 960 MHz, and the only stored results that
reach below 960 MHz are the two one-at-a-time searches named above, which the cadence plan replaces (it has another
key). No stored result that is still read changes a label.

## Radio search version 2 (`bl_local_search_v2`) - 2026-11-18, configuration unchanged (1.6.0)

### Radio: six transmitter bands added to the allocation table
The table `oddstar.io.bl_radio.RFI_ALLOCATIONS_MHZ` sets `in_rfi_band`, which decides the vetting check `rfi_band`.
It is code, not a configuration key, so the configuration version does not change; the stored search results are
keyed by `SEARCH_VERSION`, which goes from `bl_local_search_v1` to `bl_local_search_v2`, so no version 1 result is
read again.

Basis: the first search of real data (Gaia DR3 2427069874188580480 = HIP 3026, Green Bank L band, one ON scan of
5 minutes with the scans of HIP 3909 before and after; 322,961,408 channels, 12,333 windows evaluated) flagged 614
windows, 525 inside the version 1 table and 89 outside it. All 89 lie in bands where transmitters work that the
version 1 table did not list (it held satellite downlinks and navigation bands only):

| band (MHz) | service (47 CFR 2.106) | flags of the 89 |
|---|---|---|
| 1300 to 1350 | aeronautical radionavigation (air route surveillance radar) | 7 |
| 1350 to 1390 | radiolocation; GPS L3 at 1381.05 | 34 |
| 1435 to 1525 | aeronautical mobile telemetry | 5 |
| 1626.5 to 1660.5 | mobile satellite uplinks (Inmarsat) | 27 |
| 1710 to 1780 | mobile telephone uplinks (AWS) | 3 |
| 1850 to 1915 | mobile telephone uplinks (PCS) | 13 |

The bands were chosen after looking at that one star's flags, so that star is not a test of the table; the next
stars searched are. The protected radio astronomy bands (1400 to 1427 and 1660.5 to 1670 MHz) and everything else
outside the table still pass `rfi_band`. Cost: a real narrowband signal inside one of these bands is marked as
interference, as it already was in the satellite bands.

Also in this search version (no effect on any number): the stored result no longer depends on whether the index
lists a file as http or https; files on `bldata.berkeley.edu` are fetched over https; an index MD5 of "-" is treated
as no MD5 (the size is still checked); the scoring phase reports progress.

## [1.6.0] - 2026-11-18

### Vanishing: moving stars and bright-star halos in the plate-source search
The first full run of `search vanishing --live` on the 1,000-star halo sample (run `v140-halo-vanishing`, config
1.4.0) searched 3,012 plate sources and flagged 207. Read-only analysis of that run's cached answers showed two
classes of false flags, and this version removes them. The halo radii below were fitted to that same run, so the
next run on this sample is not an independent test of them.

- **Moving stars.** A Gaia star now also counts as the counterpart of a plate source when its path from 1949 to
  2002 (Gaia proper motion) passes within `match_radius_arcsec` of the plate-epoch position of the source. Before,
  the comparison was made only at the catalogue's position epoch, which is wrong when USNO-B1.0 lists a wrong proper
  motion. In `v140-halo-vanishing`, 18 of the 207 flagged sources have a Gaia path within 3.0", among them the
  sample star HIP 59490 (0.49"/yr) and seven more of the ten brightest flags. No threshold changed
  (`match_radius_arcsec: 3.0` is reused), so the other two of those ten, the sample stars Gaia DR3
  2665857750766417536 and 1744837439030364416, whose paths pass 3.45" and 4.02" from their plate positions, stay
  flagged; the `proper_motion` note of such a flag names the star and the distance. Cost: a plate source that truly vanished on the path of a moving star is
  matched to that star; a star at 0.5"/yr sweeps about 1.4% of a field.
- **Bright-star halos.** New key `channels.vanishing.plate_bright_star_halo: [[10.0, 90.0], [11.0, 50.0], [12.0, 20.0]]`
  ([G limit, radius in arcsec]). A plate source within the radius of a Gaia star brighter than the limit, and not
  that star itself, is not searched (`near_bright_star`, counted in `n_near_bright_star`). Basis: flagged fraction
  of searched sources by magnitude of the sample star in `v140-halo-vanishing`: 18 of 30 for G < 10 (at every
  distance out to the 60" field edge), 22 of 60 for G 10 to 11 (0 of 9 beyond 50"), 14 of 118 for G 11 to 12
  (6 of 14 at 10" to 20", 6 of 73 beyond), 3 to 5% for fainter stars. Gaia is now read 90" beyond the field edge
  (was 30"), so the field answers of Gaia are fetched again.
- `vet/vanishing_checks.py`: for a plate source without a Gaia counterpart, `proper_motion` used to pass by
  comparing the plate position with itself. It now passes only on the recorded distance of the nearest Gaia path
  (`vanishing_data.nearest_gaia_track`), fails when that is within the match radius, and is `not_evaluated` when
  the record is missing (units of runs before this version).
- The channel code did not change. Flags of this channel remain leads: the remaining rate (3 to 5% of searched
  sources in the 1.4.0 run) is dominated by sources seen on the first-epoch plates only, and the false-alarm rate
  of the live plate search is still not measured.

## [1.5.0] - 2026-11-18

### The radio search reads real Breakthrough Listen scans
Until now `search radio_narrowband --live` was refused: the channel had only its built-in cadences. A fine-resolution
archive scan is one block of 16 spectra by hundreds of millions of channels, which the channel cannot take whole.

- New search of scans on local disk (`oddstar/io/bl_radio.py`, used by `search radio_narrowband --live`). It walks the
  ON scan in blocks, cuts every window that stands out from the ON and OFF scans, and hands each to the channel as an
  ordinary cadence. See `docs/channels/radio_narrowband.md` for the steps, the not-searched reasons
  (`n_no_coverage`, `n_not_downloaded`, `n_no_comparison`, `n_failed`) and the known limits.
- New keys under `channels.radio_narrowband`, all chosen before any archive scan was searched:
  - `live_block_channels: 1048576`. One coarse channel of the Green Bank fine-resolution product; memory only.
  - `live_window_channels: 2048`. 5.7 kHz at 2.8 Hz per channel: wider than the largest drift (4 Hz/s) carries a
    line over one scan plus the gap to the neighbouring scan, so the OFF window holds where the line would be.
  - `live_prefilter_snr: 8.0`. Below the channel's `min_snr: 10.0`, so a window is never lost because the block
    search and the channel estimate the noise from slightly different stretches. Noise alone reached 6.6 in a
    1,048,576-channel test block.
  - `live_max_listed_windows: 200` limits the list kept per scan; the counts are of all windows.
- The channel code did not change, and no threshold changed. `vet/radio_checks.py` changed in wording only: a live
  star counts its inputs when `query_status["breakthrough_listen"]` is "ok", which the live search now writes; the
  note of a missing input says whether no scans were compared or the search does not produce that input.
  `multi_target` and `reobservation` stay `not_evaluated` for a live star.
- `search --live` no longer refuses any channel.
- The archive index applies its row limit before cutting to the requested circle, so the earlier limit of 500 could
  drop files silently. Requests now use 20000 and an answer of 5000 rows or more is refused.

Known limits: the certificates are measured on the built-in cadences; the false-alarm rate of the live search on
archive scans is not measured; a scan set with a single ON scan cannot test persistence, so interference present
only during that scan is flagged and stays a lead until a second epoch exists.

## [1.4.0] - 2026-11-18

### The vanishing search reads real plate sources
Until now no code fetched a plate record for a live object and `search vanishing --live` was refused (1.2.0). The
channel had only ever run on fixtures. A live sample star is in Gaia by construction and cannot vanish from it, so the
live search has a different unit: the USNO-B1.0 plate sources around each sample position.

- New survey fetch (`oddstar/io/plate_survey.py`, used by `search vanishing --live`): USNO-B1.0, Gaia DR3 and
  Pan-STARRS1 DR2 for a field around every sample position. See `docs/channels/vanishing.md` for the steps and the
  not-searched reasons (`n_too_faint`, `n_no_coverage`, `n_failed`).
- New keys under `channels.vanishing`, all chosen before any field was searched:
  - `plate_field_radius_arcsec: 60.0`. A high-latitude field of this size holds a handful of plate sources; 1,000
    fields are 0.87 square degrees.
  - `plate_ps1_min_detections: 2`. A Pan-STARRS object with fewer detections is not taken as a counterpart or into the
    depth. A low value is the cautious choice for a search that flags absence: a doubtful modern object still counts
    as "something is there".
  - `plate_depth_percentile: 90.0` and `plate_depth_min_objects: 10`. The local depth is the 90th percentile of the
    r-band PSF magnitudes of the field's Pan-STARRS objects (22.6 in the one field looked at to check the query
    format). It is where detections become rare, not where they are complete; the existing `min_depth_margin_mag: 2.0`
    is applied on top, so a searched plate source is at least 2 mag brighter than that.
  - `plate_max_attempts: 3`, `plate_backoff_s: 2.0`, `plate_rate_limit_rps: 2.0` are operational.
- The channel code did not change. Thresholds `min_early_plates`, `min_depth_margin_mag`, `match_radius_arcsec` and
  `max_pm_residual_arcsec` are as in 1.3.0. `vet/vanishing_checks.py` changed in its description only (which inputs a
  live object can have); every certificate value of all eight channels is identical to 1.3.0.
- `search --live` now refuses only `radio_narrowband`.

Known limits:
- The certificate describes a Gaia star paired with a plate record. The rate of false USNO-B1.0 entries (plate flaws,
  halos of bright stars, merged detections) is not part of it, so the false-alarm rate of a live run is not measured
  and its flags are leads. The null report of such a run carries the fixture grid.
- The two first-epoch plates are independent emulsions, usually of the same night, not independent epochs.
- `sso_plate`, `blend_split` and `known_variable` have no live input and stay `not_evaluated`, so a live flag cannot
  clear every check.
- South of Dec -30 there is no Pan-STARRS coverage and nothing is searched.
- `version: "1.4.0"` (was "1.3.0"). The config hash changes; every certificate is regenerated in a separate commit.

## [1.3.0] - 2026-11-18

### The spectral-line search reads real survey spectra
Until now no code fetched a spectrum for a live star, and `search spectral_lines --live` was refused (1.2.0). The
channel had only ever run on fixtures.

- New survey fetch (`oddstar/io/lamost_bulk.py`, used by `search spectral_lines --live`): two epochs per star from the
  LAMOST DR12 v2.0 low-resolution survey, reduced to the window around H-alpha in the stellar rest frame. See
  `docs/channels/spectral_lines.md` for the steps and the not-searched reasons (`n_single_epoch`, `n_no_coverage`,
  `n_failed`, `n_not_fetched`, `n_no_spectrum`). Only a star with two usable epochs from different nights is searched.
- New keys under `channels.spectral_lines`: `lamost_min_snr_r: 20.0` is a selection cut (an epoch with a lower red
  signal-to-noise is not used; the usual quality cut for LAMOST low-resolution work). It was chosen before any star was
  searched. `lamost_max_attempts: 3`, `lamost_backoff_s: 2.0`, `lamost_rate_limit_rps: 2.0` are operational.
- `SpectralLinesChannel.measure` now follows the line-spread width it is given. The H-alpha templates are broadened in
  quadrature by the excess over 0.45 A (the width they were defined at). A known line is excluded out to
  1.5 x `lsf_fwhm_a` when that exceeds the configured width: a matched filter still answers a line one FWHM away with
  a quarter of its peak, so at R = 1800 the configured 1.5 A let a strong sky line flag from just outside its zone.
  When `spectra.sky_line_wavelengths_a` is present, those positions are excluded (all of them, geocoronal H-alpha
  included) in place of the laboratory list; a spectrum in the stellar rest frame has its sky lines displaced by the
  star's velocity. Exposure arrays of unequal length are not measured. `sky_line` and `stellar_emission` in vetting
  use the same list and width. No threshold changed (`residual_snr_threshold`, the 0.4 repeat ratio, both exclusion
  widths and every `vetting.spectral_*` value are as in 1.2.0).
- Effect on validation: none. The fixtures have `lsf_fwhm_a = 0.45` and no `sky_line_wavelengths_a`, so every code path
  above reduces to the 1.2.0 arithmetic. Every certificate value of all eight channels was compared with 1.2.0: blind
  recovery, impostor tiers, completeness grids and false-alarm rates are identical, and one `spectral_lines`
  calibration-curve statistic per release moved by about 1e-14 (rounding of the rewritten template expression).
- `search --live` still refuses `vanishing` and `radio_narrowband`.

Known limits:
- The certificate describes the channel at a line-spread width of 0.45 A. Nothing was injected at R = 1800: the
  completeness and false-alarm rate of a LAMOST search are not measured, and its flags are leads, not calibrated
  detections. The null report of such a run carries the fixture grid and must not be read as a LAMOST sensitivity.
- The public LAMOST file has no sky spectrum per fibre and no detector map: `sky_line` can fail but not pass, and
  `detector` is `not_evaluated`, so a LAMOST flag cannot clear every check.
- About half of the 150 A window is excluded around sky and stellar lines at this resolution.
- `version: "1.3.0"` (was "1.2.0"). The config hash changes; every certificate is regenerated in a separate commit.
  Runs made under 1.2.0 (`v120-*`) keep their hash; to fuse channels they are run again under 1.3.0, from the cache.

## [1.2.0] - 2026-11-18

Two changes. Both were decided from how the pipeline behaved, not from which stars it flagged, and neither can
turn a failed or open check into a pass.

### The optical-dip search reads real TESS light curves for every sample star
Until now a live `search optical_dips` measured empty arrays for every star and reported "searched N, flagged 0":
a live sample star only received a light curve after some other channel had flagged it. The channel was therefore
never run on real data.

- New survey fetch (`oddstar/io/tess_bulk.py`, used by `search optical_dips --live`): the QLP `SAP_FLUX` light curve
  of each sector, good-quality cadences averaged into 30-minute bins, each sector divided by its median. New
  operational keys under `channels.optical_dips` (not thresholds): `tess_bulk_max_sectors: 8` (the sector limit of the
  vetting fetch, so vetting sees the sectors the search saw), `tess_bulk_max_attempts: 3`, `tess_bulk_backoff_s: 2.0`,
  `tess_bulk_workers: 4`, `tess_bulk_rate_limit_rps: 4.0`.
- A star counts as searched only when at least one light-curve segment has 20 points inside the edge margins. A star
  with no TIC entry, no sector, no QLP file, a failed fetch, or one beyond `--tess-max-stars` is reported as not
  searched, with the reason, and is excluded from the null result (`n_no_coverage`, `n_failed`, `n_not_fetched`,
  `n_no_light_curve`).
- `OpticalDipsChannel.measure` now works on segments. The light curve is cut where the sector changes or cadences
  are more than `tess_segment_gap_days` apart, and each segment is detrended and searched alone. The old code ran one
  median filter and one box filter across orbit gaps and sector boundaries, where the per-sector normalisation leaves
  a step. The first and last `systematics_edge_margin_days` of a segment are not searched: thermal settling and
  scattered light sit there, and `systematics_window` cannot clear an event in that zone. `dip_depth` is the median
  deficit inside the best box (was: the single lowest point within a few cadences of the peak), so one bad cadence is
  not a dip. The periodicity veto runs on the sector that holds the best box. No threshold changed
  (`dip_snr_threshold`, `min_dip_depth_frac`, `min_duration_hours`, `periodic_bls_veto_snr` are as in 1.1.2).
- Effect on validation (fixtures, dr3 and dr4 identical): blind recovery 3/3, impostors at most Explained, false-alarm
  rate 0, design point (0.15 depth) 100%, all unchanged. Completeness at the shallowest injected depth fell, because
  the median inside a box wider than the dip is below `min_dip_depth_frac`:

  | duration (days) | depth | 1.1.2 | 1.2.0 |
  |---|---|---|---|
  | 0.5 | 0.04 | 1.00 | 0.70 |
  | 1.5 | 0.04 | 1.00 | 0.85 |

  All other grid cells are unchanged (3.0 days at 0.04 was and is 0.00). The loss is real and is the price of not
  flagging single low points; the null report carries the new grid.
- `search --live` now refuses `spectral_lines`, `vanishing` and `radio_narrowband`. No live ingestion provides their
  inputs, so a live run would again report a search that did not happen (`radio_narrowband --live` crashed instead).

### A star that varies in NEOWISE cannot pass `systematics_window` on its TESS data
The first live `ir_timedomain` runs (1,800 stars, 30 flagged) showed the NEOWISE half of `systematics_window`
almost never ran: the field-star ensemble had 0 or 1 stars for nearly every target, the half returned
`not_evaluated`, and the check passed on the TESS half alone. A check that has not looked at the infrared frames was
passing stars flagged for an infrared change. Eight stars of one run shared one faint W2 trend, which is what an
untested zero-point drift would look like.

- `check_systematics_window`: when the star's NEOWISE light curve varies by the channel's own thresholds and the
  NEOWISE half is `not_evaluated`, the check is `not_evaluated`, whatever the TESS half says. A fixture's explicit
  `scan_systematic` flag still counts as an evaluation. This can only move a result from `pass` to `not_evaluated`.
- `channels.ir_timedomain.neowise_comparison_radius_arcsec: 600.0` (was 150.0). A 150" cone around a W1 = 11.9 halo
  star holds one other star within `neowise_comparison_mag_window` (1.0 mag); 450" holds seven (AllWISE count around
  Gaia DR3 3580358134642783872, made to size the cone, not to judge that star). The zero-point and scan systematics
  the ensemble tests are common to a 47' frame, so stars 10' away see them. `neowise_comparison_mag_window`,
  `neowise_comparison_min_stars`, `neowise_comparison_max_stars`, `neowise_comparison_exclude_arcsec` and the three
  `systematics_ensemble_*` thresholds are unchanged.
- The cone is filtered by the archive to the target's W1 +/- (window + 0.25 mag) and to six columns; unfiltered it
  would be hundreds of thousands of single exposures. The margin keeps a star at the edge of the window from losing
  part of its exposures.
- Stars brighter than about W1 = 8 will usually still have no comparison stars. They now stay `not_evaluated`.

Known limits:
- Results of the 1.1.2 `ir_timedomain` runs (`v113-*`) were vetted with the old rule; their `systematics_window`
  passes for varying stars are not evidence. They must be vetted again under 1.2.0 before any tier is quoted.
- `neighbour_signal` combines its TESS and NEOWISE halves the same way and was not changed here.
- `version: "1.2.0"` (was "1.1.2"). The config hash changes; every certificate is regenerated in a separate commit.

## [1.1.2] - 2026-11-17

### The image-shape check is judged against measured stars
`vision_rubric` calls a band "extended" when its fitted FWHM exceeds a reference by `vision_extended_fwhm_ratio`.
In 1.1.0 and 1.1.1 the reference was theoretical: native beam x `wise_atlas_psf_broadening` (9.19" in W3, 16.97"
in W4). The first live run measured the reference star at 1.40 x that in W3 and 1.31 x in W4 and failed it.
Whether that meant extended emission or a reference that does not describe this fit could only be settled on
known point sources, so the fit was run on ordinary stars.

Calibration (`python -m oddstar.validate.atlas_psf`, output committed in ab488cb before this change):
63 bright, isolated, unsaturated stars with photospheric W3-W4 in 12 fixed fields at |b| > 30 deg, none near a
pipeline candidate, each measured by the unchanged `measure_wise_cutouts`. All 12 selection queries and all 63
cutout sets succeeded.

| Band | n | median fitted FWHM | robust sigma | min - max | median / theoretical reference |
|---|---|---|---|---|---|
| W3 | 63 | 12.971" | 0.099" | 12.85" - 13.32" | 1.41 |
| W4 | 63 | 20.597" | 0.487" | 17.84" - 23.76" | 1.21 |

So an ordinary star comes out of this fit 1.41 x (W3) and 1.21 x (W4) wider than the theoretical reference, and the
1.1.1 check would have called essentially every star extended in W3.

The selection and the rule for using the result were written in `oddstar/validate/atlas_psf.py` and committed
(160a8f9) before the calibration ran. Applied as written:

- `vetting.w3_atlas_point_fwhm_arcsec: 12.971` and `vetting.w4_atlas_point_fwhm_arcsec: 20.597` (new): the median
  fitted FWHM of the calibration stars. `fwhm_ratio` is now fitted FWHM over this value. The geometry of the
  measurement (aperture, fit box, peak search) still uses native beam x `wise_atlas_psf_broadening`, so fitted
  widths are unchanged and stay comparable with the calibration file.
- `vetting.vision_extended_fwhm_ratio: 1.3` (unchanged). The rule was max(1.3, 1 + 4 x robust sigma / median) over
  the two bands, rounded up to 0.05; 1 + 4 sigma / median is 1.03 in W3 and 1.09 in W4, so the rule returns 1.3.
- `version: "1.1.2"` (was "1.1.1"). The config hash changes; every certificate is regenerated in a separate commit.

Known limits, recorded so they are not mistaken for results:
- 1.3 is loose against the measured scatter (about 30 sigma in W3, 12 sigma in W4). The rule could only keep or
  raise the threshold. Tightening it is a separate decision that must be made on calibration data alone.
- Fitted width rises slightly with S/N (Spearman rho +0.43 in W3, +0.32 in W4), and the W4 scatter is larger for
  faint stars: the 19 calibration stars with fit S/N below 25 span 17.8" - 23.8" (robust sigma 0.81").
- The calibration stars have catalogue W4 S/N of 10 or more. A fainter target is outside the calibrated range.

No search result is part of this change.

## [1.1.1] - 2026-11-17

### Fixes from the first live archive run of 1.1.0
The 2026-11-17 live run of the reference star against the real archives exposed answers the 1.1.0 code misread.
No channel threshold changes. Three keys are new; the config hash changes, so every certificate is regenerated
(in a separate commit).

- `version: "1.1.1"` (was "1.1.0").
- `vetting.sso_max_pos_err_arcsec: 30.0` (new). SkyBoT lists lost and single-opposition bodies whose ephemeris
  uncertainty reaches degrees (seen live: 2021 VU98, 225,012" uncertainty, nominal position 561" from the star).
  The match test adds a body's uncertainty to `sso_match_radius_arcsec`, so such a body matched any target in the
  field and failed `moving_object`. A body with an uncertainty above this cap is now tested at its nominal position
  only and listed under `unconstrained_bodies` in the check's metrics and note. 30" is twice the 15" match radius
  and 2.5 times the 12" W4 beam: beyond that the catalogue position no longer says whether the body was in the beam.
- `vetting.youth_galex_coverage_radius_arcsec: 600.0` (new). GALEX did not observe the whole sky. A field with no
  GALEX source of any kind within this radius was not observed (`galex: no_coverage`), which is different from an
  observed field with no NUV source at the star. An observed AIS field holds of order a hundred catalogue sources
  within 10 arcmin.
- `channels.ir_timedomain.neowise_bulk_workers: 4` (new). Concurrent IRSA requests of the survey NEOWISE fetch; all
  workers share the `archives.neowise.rate_limit_rps` limit. Results do not depend on the worker count (tested).

Source changes without a config key (recorded here because they change which archive answers a check):
- GALEX now comes from MAST `Mast.Galex.Catalog` and ROSAT 2RXS / eRASS1 from the HEASARC TAP tables `rass2rxs` /
  `erass1main`. The CfA VizieR mirror answered HTTP 200 with zero rows for GALEX and eRASS1 everywhere (it does not
  serve those tables), which was recorded as "queried, no counterpart", and the 2RXS identifier `IX/47/2rxs` does
  not exist. Match radii and flux limits are unchanged.
- The remaining VizieR `asu-tsv` tables (NVSS, PS1 fallback) accept an empty answer from a mirror only after that
  mirror returned rows at a known position for the table (`VIZIER_CANARY_CONES` in `io/live_survey.py`).
- The Legacy Survey DR10 neighbour query is a 17" `q3c_radial_query` cone (circumscribes the former +/-12" box,
  which timed out); rows are still filtered by separation downstream.
- SkyBoT HTTP 204 is "no body in the field", not a failure. IRSA image-search columns are read by name. A MAST
  answer counts only when its status is `COMPLETE`. `_youth_uv` dereddens with the reddening applied to the star
  (`ebv_star`), not the line-of-sight key, which is empty for live stars since 1.1.0.
- Operational only (no effect on results): Bayestar19 stops being asked for 30 minutes after a failed connection;
  a TAP service answering HTTP 429 is retried after 20, 60, 180 and 480 s; halo and bulge samples skip the
  single-shot join that timed out and use the index + batch route unless a single-shot answer is already cached.

## [1.1.0] - 2026-11-17

### Infrared vetting inputs
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
- Note: the top-level `version:` field was still "1.0.0" when the entries above were written; it is bumped to "1.1.0" below.

### Governance and thresholds moved from source
Decision cuts that were literals in the source now live in `configs/default.yaml`. **Every value is unchanged**, so
no vetting outcome, sample or channel score changes; the live sample SQL text is byte-identical (tested in
`tests/test_governance_config.py`). The config hash changes from `10bc507f1b29` to `c4ad1ae04d32`, so validation
certificates must be regenerated for the new hash (not done here).

- `version: "1.1.0"` (was "1.0.0"). CI now requires `version` to equal the newest `## [x.y.z]` heading here (`check_config_version`).
- R4 commit-separation cutoff: `7bf6ed5f4be3050e91b99f87503bbb20078b9604`. `check_git_commit_separation` checks HEAD and, in pull-request CI, every non-merge commit of the base range, but only commits after this one; earlier history is not rewritten.
- `vetting` (from `vet/engine.py` `vet_target`): `outlier_driver_min_abs_contribution: 2.0` (generic_outlier driver contribution that routes to a check).
- `vetting` (from `vet/ir_checks.py`):
  - `check_xmatch_integrity`: `xmatch_max_sep_arcsec: 2.5`.
  - `check_wise_quality`: `wise_w4_min_frames: 5`, `wise_w4_min_detected_frame_frac: 0.15`, `wise_w4_max_rchi2: 3.5`, `wise_max_ext_flg: 1`.
  - `check_extragalactic_colours`: `dog_w1_w2_min: 0.45`, `dog_w2_w3_min: 3.0`.
  - `check_single_frames`: `single_frame_consistency_frac: 0.5`, `single_frame_coadd_floor_frac: 0.6`, `single_frame_max_to_median: 2.2`.
  - `_youth_uv`: `youth_nuv_default_err_mag: 0.3` (NUV error adopted when GALEX gives none).
  - `check_evolved`: `evolved_teff_max_k: 5200.0`.
  - `check_binary`: `binary_rv_scatter_max_kms: 5.0`.
- `vetting` (from `vet/lc_checks.py`):
  - `matched_detrend` / `analyze_tess_sector`: `lc_detrend_low_mask_sigma: 2.5` (the function keeps 2.5 as its default argument for direct callers).
  - `_tess_systematics`: `lc_tess_contiguous_dt_factor: 1.5`.
  - `_sector_period_search` and `check_known_variable` (BLS): `lc_bls_min_period_days: 0.8`, `lc_bls_max_period_days: 8.0`, `lc_bls_duration_days: 0.1`.
  - `_sector_period_search` (eclipse confirmation): `lc_eb_eclipse_min_depth_frac_of_model: 0.5`, `lc_eb_min_eclipses: 2`, `lc_eb_min_eclipse_epoch_frac: 0.5`.
  - `_ir_event_window`: `lc_ir_event_baseline_visits: 3`, `lc_ir_event_dev_frac: 0.5`, `lc_ir_event_min_window_days: 180.0`, `lc_ir_event_pad_days: 90.0`.
- `vetting` (from `vet/spectral_checks.py`):
  - `check_cosmic_ray`: `spectral_cr_spike_mad_sigma: 8.0`, `spectral_peak_width_flux_level: 1.3`, `spectral_cr_min_peak_flux: 1.5`, `spectral_cr_min_width_lsf_frac: 0.7`.
  - `check_sky_line`: `spectral_sky_peak_mad_sigma: 6.0`, `spectral_target_peak_sigma: 3.5`.
  - `check_stellar_emission`: `spectral_stellar_line_min_amp: 0.20`.
  - `check_repeat`: `spectral_repeat_quiet_max_excess: 0.25`, `spectral_repeat_max_pixel_shift: 2`, `spectral_repeat_min_excess: 0.20`.
- `vetting` (from `vet/vanishing_checks.py` `check_proper_motion`): `vanishing_pm_min_raw_sep_arcsec: 3.0`.
- `fusion` (from `fuse/engine.py` `calibrate_fap` / `fuse_results`): `fap_floor: 1.0e-12`, `fallback_false_alarm_rate_per_star: 1.0e-4` (`calibrate_fap` keeps `floor=1e-12` as its default argument for callers that pass none).
- `live_ingestion` (new section):
  - from `io/live_survey.py` `_phot_quality_where_sql`: `allwise_number_of_neighbours: 1`, `allwise_number_of_mates: 0`, `tmass_number_of_neighbours: 1`, `tmass_number_of_mates: 0`.
  - from `io/live_survey.py` `fetch_live_halo_bulge_sample` region presets, `regions.halo`: `min_parallax_over_error: 10.0`, `parallax_min_mas: 0.5`, `parallax_max_mas: 20.0`, `max_ruwe: 1.4`, `min_abs_b_deg: 20.0`, `g_mag_min: 7.5`, `g_mag_max: 14.5`, `teff_min_k: 3800.0`, `teff_max_k: 7500.0`, `max_mh_gspphot: -0.8`, `min_vtan_kms: 150.0`.
  - `regions.bulge`: `gal_l_half_width_deg: 12.0` (the query keeps `l <= 12.0 OR l >= 348.0`), `abs_b_min_deg: 3.0`, `abs_b_max_deg: 10.0`, `parallax_min_mas: 0.08`, `parallax_max_mas: 0.45`, `min_parallax_over_error: 5.0`, `max_ruwe: 1.4`, `g_mag_min: 8.0`, `g_mag_max: 14.5`, `teff_min_k: 3500.0`, `teff_max_k: 7500.0`.
  - `regions.local`: `min_parallax_over_error: 15.0`, `parallax_min_mas: 2.0`, `parallax_max_mas: 25.0`, `max_ruwe: 1.4`, `min_abs_b_deg: 12.0`, `g_mag_min: 7.5`, `g_mag_max: 13.5`, `teff_min_k: 3800.0`, `teff_max_k: 7500.0`, `min_vtan_kms: 65.0`.
  - from `cli.py` `filtered_live_ingestion`: `ir_excess_w4_upper_limit_ph_qual: ["U", "X"]` (W4 `ph_qual` codes treated as upper limits), `ir_excess_artifact_cc_flag_bands: [3, 4]` (bands whose `cc_flags` character marks an artifact), `ir_timedomain_neowise_enrich_max_stars: 25` (NEOWISE enrichment cap).
- `storage.evidence_dir: "registry/evidence"` and a new `evidence` section (`max_array_elements: 64`, `max_string_chars: 2000`, `max_json_bytes: 65536`, `max_image_bytes: 262144`): size guards of the committed evidence store. Storage policy, not decision thresholds.
- Deliberately left in source: statistical and unit constants (1.4826, 1.2533, 2.5 log10, 4.74047, 365.25, MJD/BTJD offsets, catalogue epochs 2016.0 / 2010.2, R_V 3.1, physical constants), numerical guards (1e-9, 1e-6, 1e-4 floors), array indices and window offsets, minimum sample sizes for a statistic to be computed (points per sector or segment, bins, seasons, nights), the BLS grid size 200, the `0.5 x span` period bound, data defaults used when an input is missing (spectral `lsf`/pixel size, `abs_g`/`teff` defaults, USNO epoch/magnitude, PS1 depth, plate count), enrichment-fetch query radii, timeouts, chunk sizes and the `n_lim` clamp in `live_survey.py`, `_row_to_star_dict` and the end of `enrich_star_live_for_vetting` (out of scope here), `vet/underluminous_checks.py`, and the dossier writer's threshold text.

### Star model: empirical photosphere and explicit extinction
- Added `channels.ir_excess` empirical-locus keys (moved from hard-coded values where they existed; no existing threshold changed): `locus_bp_rp_min: -0.2`, `locus_bp_rp_max: 2.8`, `locus_n_bins: 15`, `locus_min_stars_per_bin: 5` (same binning as before), `locus_clip_sigma: 3.0`, `locus_clip_max_iter: 5`, `locus_intrinsic_scatter_floor_mag: 0.01`, `locus_extrapolation_scatter_mag: 0.10`, `locus_persist_min_stars: 100`, `literature_locus_scatter_mag: [0.03, 0.04, 0.05, 0.08]` (assumed W1..W4 scatter for the literature fallback; the table has none).
- Added `channels.ir_excess` dust-quantification keys: `excess_significance_sigma: 3.0`, `dust_temp_interval_delta_chi2: 1.0`, `colour_temp_min_k: 30.0`, `colour_temp_max_k: 3000.0`.
- Added `vetting` extinction keys: `extinction_heavy_dust_ebv: 0.30` (previously hard-coded in `check_extinction`, same value), `extinction_los_excess_summary_ebv: 0.05`, `sfd_to_sf11_factor: 0.86`, `bayestar19_ebv_per_unit: 0.884` (Green et al. 2019 Eq. 1; Eq. 2 gives 0.996; to be confirmed in a live run).
- Behaviour: the ir_excess primary photosphere is now the dereddened W1..W4 empirical locus (`photosphere_source` recorded); the empirical-predictor chi now includes the locus intrinsic scatter (was multiplied by 0.0) and no longer adds `cal_uncertainty_frac` (already inside the locus scatter); fixture photospheres carry the dwarf WISE colours; `ebv_star` replaces `ebv_sfd` wherever the star is dereddened.

### Light-curve vetting inputs
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

### generic_outlier re-certification, code fingerprints, survey-scale NEOWISE (wt/survey-scale)
Keys added only; no existing threshold or gate changed. The config hash changes, so every channel must be re-validated;
certificates now also carry `code_fingerprint` and go stale when the validated code changes.
- `channels.generic_outlier.modality_null_default_n_points: 1000` - series length assumed for the noise expectation of a
  `dip_snr` / `spec_line_snr` value supplied without its series (`sub_metrics`), i.e. E[max of 1000 N(0,1)] = 3.23 (Blom approximation).
- `channels.generic_outlier.vanished_margin_plate_sigma_mag: 0.3` - photographic-plate magnitude uncertainty (USNO-B1.0,
  Monet et al. 2003) used to express `vanished_margin` in sigma.
- Behaviour (generic_outlier): a present modality feature is now an excess over its own noise expectation (absent stays 0.0,
  L003); the Isolation Forest draws split attributes from all features (Liu et al. 2008), so a feature constant in the
  background can still isolate a query; `spec_line_snr` is the strongest single-exposure feature (repeat test left to the
  spectral vetting suite, like `dip_snr`). See `docs/channels/generic_outlier.md`.
- `channels.ir_timedomain.neowise_bulk_max_attempts: 3`, `neowise_bulk_retry_backoff_s: 2.0` - bounded retries (exponential
  backoff 2 s, 4 s) of the survey NEOWISE-only fetch (`oddstar.io.neowise_bulk`); its rate limit is `archives.neowise.rate_limit_rps`.
- Behaviour (ir_timedomain): stars without NEOWISE coverage, with a failed or missing NEOWISE fetch, or with fewer than
  `min_visits` visits are reported as not searched and excluded from null-report denominators (detection thresholds unchanged).
- Removed `live_ingestion.ir_timedomain_neowise_enrich_max_stars: 25`. A live `search ir_timedomain` used to fetch NEOWISE
  for the first 25 sample stars only and drop the rest; it now fetches every sample star through the bulk path
  (`--neowise-max-stars N` bounds a run explicitly, and the stars beyond N are reported as not searched).

### Remaining default-value passes removed
Vanishing, spectral and radio vetting, the SIMBAD-based checks and the generic_outlier path no longer pass on
default input values (lessons L007/L008): a check passes or fails only on inputs that are explicitly present (fixtures,
hand-built dicts) or written by a query recorded as `ok` (live stars); otherwise it reports `not_evaluated` and names the
missing data product. No existing threshold value changed. The config hash changes, so validation certificates must
be regenerated for the new hash (not done here).

- `channels.vanishing.vetting_variable_simbad_types: ["Mira", "LPV*", "CV*", "Flare*"]` - SIMBAD types that explain a
  vanished plate source in `known_variable`; was a literal set in `vet/vanishing_checks.py` (unchanged).
- `channels.optical_dips.detrend_edge_zone_window_frac: 0.5` (new) - dip candidates within this fraction of
  `detrend_window_days` (2.5 d) of a TESS segment edge or data gap, where the matched running median can erase or
  distort a dip, make `pipeline_agreement`, `systematics_window`, `neighbour_signal`, `sso_crossing` and the TESS context
  of `second_survey` report `not_evaluated` instead of `pass` (a `fail` stands). Edge-zone bins are also compared with the
  raw interior flux so an erased dip is still seen.
- `vetting.metallicity_subdwarf_feh_max: -1.0` - [Fe/H] below which `check_metallicity` calls a subdwarf; was the
  literal `-1.0` in `vet/underluminous_checks.py` (unchanged).
- Literal fallbacks of config keys removed from vetting source (`.get(key, default)` -> required key, values unchanged):
  `channels.vanishing.min_early_plates` (2), `min_depth_margin_mag` (2.0), `max_pm_residual_arcsec` (1.5),
  `channels.radio_narrowband.min_drift_abs_hz_s` (0.05), `channels.spectral_lines.sky_line_exclusion_angstrom` (1.5),
  `stellar_line_exclusion_angstrom` (1.5).
- Input defaults removed (now `not_evaluated` when absent): vanishing `usno_epoch` 1953.4, `usno_mag_r` 17.0,
  `panstarrs_depth_mag` 21.8, `early_plates_count` (2 with any plate record), `known_lpv` False, `sso_on_plate` False,
  `blend_split_sources` 1, `present_in_gaia` True, proper motion 0.0; spectral `lsf_fwhm_a` 0.45, pixel size 0.25 A,
  `detector_artifact` False, empty flux arrays; radio `off_scan_hits` 0, `drift_hz_s` 0.0, `in_rfi_band` False,
  `multi_target_match` False, `reobserved_confirmed` False; SIMBAD type "Star"; generic path `abs_g_mag` 4.8,
  `teff_k` 5800, `ruwe` 1.0, proper motion 0.0 (xmatch, NEOWISE PM track), radio row separation 99",
  second-instrument interpretation "confirmed", fixture cutout offset 0.0" / background 0.4 mJy.
- SIMBAD: live ingestion records `query_status["simbad"]` and `vetting_archive.simbad_entry`; `simbad_type` is None when
  the lookup failed or SIMBAD has no entry. `known_variable` (light curve and vanishing) and `metallicity` return
  `not_evaluated` when the lookup did not succeed; "no entry" is a successful negative answer.
- `vet/engine.py`: a generic_outlier target with a dip/IR driver is routed to the light-curve checks also when it carries
  only per-sector TESS data (`lightcurves.tess_sectors`, the live format), not just the legacy fixture arrays.
- Fixture-only channels for live stars: vanishing, spectral_lines and radio_narrowband have no live ingestion, so every
  check of those suites returns `not_evaluated` for a live star (a SIMBAD-listed Mira can still fail vanishing
  `known_variable`).

## [1.0.0] - 2026-11-16
- Initial frozen baseline configuration for all 8 search channels, 23 archive adapters, and 42 vetting checks.
- Set `storage.max_cache_gb: 10.0` per workstation disk budget decision (with process-and-evict retention of manifest SHA-256 hashes for heavy raw radio/pixel payloads).
- Set `archives.hsc_ssp.enabled: false` pending explicit human license/registration acceptance.
