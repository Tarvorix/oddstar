# Channel: `vanishing`

## Physical Signature
Point sources confidently detected on two or more 20th-century photographic plates/bands (USNO-B1.0) that are absent in deeper modern surveys (Pan-STARRS, Gaia) after proper-motion propagation.

## Inputs
- Historical USNO-B1.0/DASCH/DSS positions and magnitudes paired with Pan-STARRS and Gaia DR3/DR4 astrometry and local depth limits.

## Algorithm
1. Propagate modern Gaia coordinates back to historical plate epochs using `(pmra, pmdec)`.
2. Require `>= 2` independent early plate detections (`early_plates_count >= 2`) and modern survey non-detection with local depth exceeding early magnitude by `>= 2.0 mag`.

## Live survey search (config 1.4.0, matching changed in 1.6.0)
A sample star is in Gaia by construction, so it cannot be the unit of this search. `oddstar search vanishing --live`
searches **plate sources**: the USNO-B1.0 objects around each sample position (`oddstar/io/plate_survey.py`).
`--limit` is the number of sample positions; `n_searched` counts plate sources.

For each field of radius `plate_field_radius_arcsec` three catalogues are read (cached, rate-limited, retried):

1. USNO-B1.0 (VizieR I/284). A source is a unit when it was detected on both first-epoch plates (blue and red) and is
   not flagged as lying on a diffraction spike. Its position is at epoch 2000 when the catalogue lists a proper motion,
   and at the mean plate epoch when it lists none.
2. Gaia DR3 out to 90" beyond the field edge (the largest halo radius below), moved to the epoch of the plate
   position with its proper motions. A source within `match_radius_arcsec` is the modern counterpart; the unit then
   carries the Gaia position and proper motion, which is what the channel propagates.
   Since 1.6.0 a Gaia star is also the counterpart when its path from 1949 to 2002, drawn with the Gaia proper
   motion, passes within `match_radius_arcsec` of where the plates saw the source (the catalogue position taken
   back to the mean plate epoch with the catalogue's own motion). The catalogue's proper motion can be wrong, and
   its epoch-2000 position is then far from a star that never left: in the first full halo run the brightest
   "vanished" sources were the fast-moving sample stars themselves. The comparison epoch of such a match is the
   epoch of closest approach (`vanishing_data.gaia_match` is `"track"`).
3. Pan-STARRS1 DR2 mean objects with at least `plate_ps1_min_detections` detections. An object within
   `match_radius_arcsec` is a counterpart. The local depth is the `plate_depth_percentile` of the r-band PSF magnitudes
   of the field's objects and needs `plate_depth_min_objects` of them.

A plate source is **searched** only when all three answers were read, the field has a depth, and the depth exceeds its
first red magnitude by `min_depth_margin_mag`. Every other one is listed with its reason and excluded from the null
result:

| Reason | Meaning |
|---|---|
| `too_faint` | closer than `min_depth_margin_mag` to the local Pan-STARRS depth: its absence would not be significant |
| `near_bright_star` | inside the halo radius of a bright Gaia star (`plate_bright_star_halo`: 90" for G < 10, 50" for G < 11, 20" for G < 12) and not that star itself; halos, spikes and ghosts repeat on both first-epoch plates |
| `no_coverage` | the field has no Pan-STARRS depth (south of Dec -30, or too few objects with an r magnitude) |
| `failed` | the Gaia or Pan-STARRS query of the field failed after `plate_max_attempts`; not cached, retried next run |

A field whose USNO-B query failed yields no plate sources at all. It is counted in `plate_fetch.n_fields_failed_usno`,
listed in `failed_field_ids`, and the search prints a warning, because its sources are in no count.

Flagged plate sources are written to `live_stars.json` as searched; the per-star archive enrichment does not apply to
an object that is not a Gaia star. Vetting evaluates `proper_motion`, `plate_defect` and `depth` from the survey record
and, since 1.8.0, from the digitized plates (`oddstar/io/plate_images.py`): `vet --live` fetches a 4 arcmin cutout of
each of the five Digitized Sky Survey plates at the source, registers it on the Gaia stars it holds, and stores what it
shows at the position (about five cutouts and four seconds per source; only the measured rows are kept).

- `plate_defect` passes when both first-epoch images show a point of light within `plate_image_max_offset_arcsec` at
  `plate_image_min_snr` or more. It fails (tier Rejected) when two first-epoch images could say, neither shows it, and
  the source is no fainter than `plate_image_judge_max_r_mag`. It is `not_evaluated` when the images disagree, when
  fewer than two could say (not delivered, or not registered and empty), or when the source is fainter than that.
  On controls: 0 of 69 blank positions pass and 57 fail; of 190 ordinary sources brighter than 19.5, 148 pass, 2 fail.
- `proper_motion`, for a plate source without a Gaia counterpart, fails when a Gaia star's path over the plate era
  passes within the match radius (`vanishing_data.nearest_gaia_track`; a unit of a run made before 1.6.0 has no such
  record and is `not_evaluated`). Since 1.8.0 it also reads the Pan-STARRS objects between the match radius and
  `plate_moved_radius_arcsec` whose r magnitude is within `plate_moved_max_dmag` of the plate magnitude. One that a
  first-epoch plate shows where it is now is another star. For any other, the later plates decide: the check fails
  ("the same star, moved") when `plate_track_min_images` later images show the source within
  `plate_track_max_residual_arcsec` of where a star travelling from the first-epoch place to the object would be at
  the plate's date; the object is ruled out when they show the source still in place; otherwise the check is
  `not_evaluated`. On controls: 43 of 57 known moving stars fail and none is ruled out; 73 of 79 still stars with a
  neighbour are ruled out and none fails (see `configs/CHANGELOG.md`, 1.8.0).
  Since 1.9.0 the check first repeats the survey's comparison at the place where the first-epoch plates show the
  source (`vanishing_data.plate_first_place`): USNO-B1.0 gives some sources a false proper motion, and the catalogue
  position is then several arcsec from an object that never left. The check fails ("the object the plates recorded
  is still there") when a Gaia star, moved back to the plate's date, or a Pan-STARRS object is within
  `match_radius_arcsec` of that place, and it passes only when this comparison and the moved-star test could both be
  made and found nothing. On controls: the counterpart is found for 328 of 328 sources known to be there, and 9.4%
  of 5,640 places without a source have an unrelated object within the match radius (see `configs/CHANGELOG.md`,
  1.9.0).
`sso_plate`, `blend_split` and `known_variable` stay `not_evaluated`: no plate-epoch ephemeris, high-resolution imaging
or variable-star lookup is made.

What a live flag means, and what it does not:
- USNO-B1.0 contains false sources: plate flaws, halos and spikes of bright stars the flag misses, and detections of
  different objects merged into one entry. The sample positions are bright stars, so their own halos are in every
  field; each unit records its distance from the sample star and that star's magnitude (`field_sep_arcsec`,
  `field_star_g_mag`).
- The two first-epoch plates of a field were usually exposed the same night. They are independent emulsions, which
  rejects most plate flaws, but not independent epochs: an asteroid or a satellite glint can be on both.
- A flag says that no catalogued modern source lies within `match_radius_arcsec`. It does not say the sky is empty
  there; only the images can.
- The halo radii were read off the first full halo run (config 1.4.0, 3,012 searched plate sources): the flagged
  fraction was 60% around sample stars with G < 10, 37% for G 10 to 11 (none beyond 50"), 12% for G 11 to 12 (the
  excess inside 20"), and 3 to 5% around fainter stars. They are a mask fitted to that sample, not a model of the
  plates, and the 3 to 5% that remains is the rate of unexplained flags per searched source, mostly sources seen
  on the first-epoch plates only.
- The certificate of this channel was measured on fixtures (a Gaia star paired with a plate record). The rate of false
  USNO-B entries per field is not part of it, so the expected number of false alarms of a live run is not measured.

## Thresholds in Config
Defined in `configs/default.yaml` under `channels.vanishing`:
- `min_early_plates: 2`, `min_depth_margin_mag: 2.0`, `match_radius_arcsec: 3.0`, `max_pm_residual_arcsec: 1.5`
- Live plate-source search: `plate_field_radius_arcsec: 60.0`, `plate_ps1_min_detections: 2`, `plate_depth_percentile: 90.0`,
  `plate_depth_min_objects: 10`, `plate_bright_star_halo: [[10.0, 90.0], [11.0, 50.0], [12.0, 20.0]]`, and the operational keys `plate_max_attempts: 3`, `plate_backoff_s: 2.0`,
  `plate_rate_limit_rps: 2.0`
- Plate images in vetting (since 1.8.0): `plate_image_min_snr: 5.0`, `plate_image_max_offset_arcsec: 2.0`,
  `plate_image_first_epoch_before: 1960.0`, `plate_image_judge_max_r_mag: 19.5`, `plate_moved_radius_arcsec: 10.0`,
  `plate_moved_max_dmag: 1.5`, `plate_track_max_residual_arcsec: 1.0`, `plate_track_min_images: 2`

## Known Impostors
- High-proper-motion stars that moved beyond the static cross-match radius, plate emulsion flaws/scratches, large-amplitude Mira/LPV variables, asteroids, and blended sources resolved into multiple modern stars.

## Benchmark Objects
- Bidirectional switch-off test on fast-PM stars + single-epoch sources from Villarroel et al. (2020).

## Failure Modes
- Nearly all raw historical-to-modern catalog mismatches (~99.9%) are high-PM stars or single-plate emulsion flaws.

## Lessons Applied
- Mandatory multi-plate verification and proper-motion track back-propagation.
- `L002_missing_plate_data_high_pm_vanishing`: When archival `vanishing_data` is absent (e.g., live multi-archive queries without USNO-B1.0 plates), do not default early-epoch coordinates to J2016 Gaia coordinates or assume `early_plates_count = 2`, preventing synthetic PM-propagation mismatches on high-PM stars.
