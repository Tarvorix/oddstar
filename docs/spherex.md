# SPHEREx photometry

SPHEREx (launched 2025) images the whole sky at 0.75 to 5.0 micron through linear variable filters and returns to
each position about every six months. For a star this gives two things no other source in the pipeline has: a
measurement at 3 to 5 micron from an instrument that is not WISE, and epochs after NEOWISE ended in 2024.

Code: `src/oddstar/io/spherex.py`. Command: `oddstar spherex photometry`.

## What is measured

- **Image list.** IRSA TAP (`spherex.artifact`, `spherex.plane`, `spherex.observation`): every Level 2 spectral
  image whose footprint holds the star. A star has a few hundred (369 for a star at RA 45, Dec 0 in November
  2026: 242 in data release 1 and 127 in quick release 3, six detectors, July 2025 to September 2026).
- **One cutout per exposure.** `<image>?center=<ra>,<dec>d&size=0.03`. A cutout is about 4 MB whatever its size,
  because it carries the detector's whole table of point spread functions. Cutouts are not kept.
- **Wavelength.** Each exposure sees the star at one wavelength, set by where the star falls on the detector: a
  detector covers a range (detector 4: 2.4 to 3.8 micron), and the wavelength comes from the file's own lookup
  table at the star's pixel.
- **Flux.** The zodiacal model in the file is subtracted. The file's effective point spread function for the
  star's part of the detector, placed at the star's position, and a sloped background are fitted to a 9 x 9
  pixel stamp with the file's variance as weights. The position is the Gaia position moved by the proper motion
  to the date of the exposure and is not fitted. Gaia neighbours within 30 arcsec that are no more than 5
  magnitudes fainter get an amplitude of their own; a neighbour closer than one pixel cannot be separated and is
  counted as blended.
- **Flags.** Pixels with a defect bit (cosmic-ray hit, saturation, non-functional, hot, cold, persistence,
  outlier, streak and the like) are left out. A defect within 1.5 pixels of the star gives no measurement
  (`core_flagged`). So does the dichroic bit on the star's core: it marks the edge of a detector's band, where
  two detectors measured the same star 20 percent apart at the same wavelength. Ghost, halo and crosstalk bits
  are reported on the row and not masked.
- **Saturation.** Stars brighter than about G = 11 saturate at the short wavelengths (the overflow bit on the
  core): the first star measured, G = 10.6, was usable in 3 of 30 exposures on detectors 1 and 2, 11 of 30 on
  detector 3, and nearly all on detectors 4 to 6.
- **Error.** The formal error of the fit, scaled up by the square root of the reduced chi-square when that is
  above one: on a bright star the model of the point spread function, not the noise, limits the fit.

Each measured row is stored in the raw cache (`spherex`, release `spherex_forced_v1`), keyed by image, position
and neighbour list. A repeat fetches nothing. A cutout that cannot be fetched is counted (`n_failed`,
`complete: false`) and is fetched by the next call.

## The summary

Per star and per visit (exposures less than 60 days apart): the weighted mean flux density in three windows,
K (2.0 to 2.4 micron), the W1 part (3.0 to 3.8) and the W2 part (4.1 to 5.0; SPHEREx ends at 5.0, W2 at 5.3), each
also as a magnitude over the 2MASS or WISE zero point, and the colours between them.

A measurement whose reduced chi-square is more than 5 times the median of the star's measurements on the same
detector is left out of the summary and counted (`n_misfit`). The first star measured had one such exposure in
22: the star lay 0.36 pixel from where that image's sky solution put it, and the forced fit read 25 percent low.

A window mean is not an integral through the 2MASS or WISE response. These magnitudes differ from catalogue
magnitudes by a colour term of a few hundredths of a magnitude: compare a star with stars measured the same way,
and one visit with the next.

## Limits

- A SPHEREx pixel is 6.15 arcsec, as coarse as WISE W1 and W2. This photometry confirms or contradicts the
  brightness of the blended light at the position. It cannot separate a neighbour that WISE could not separate.
- SPHEREx ends at 5 micron. It says nothing about an excess seen only at 12 or 22 micron.
- A star within four pixels of a detector edge is not measured in that exposure (`off_cutout`).

## Status

The reader and the fit are tested on synthetic files laid out like a real cutout (the layout was read from a real
file on 2026-11-18). On 22 real exposures of one star (2026-11-18) the point spread function as delivered fits
better than its mirror images on five of six detectors (on the sixth, three exposures, one mirror image is
slightly better), and the best position shift has a median of zero on every detector. Reduced chi-square is 1.3
to 2.4 at 3.8 to 5 micron and 6 to 70 below that, where this star is bright: there the model of the point spread
function limits the fit, and the error is scaled accordingly. The photometry has not yet been compared with 2MASS
and WISE on a set of ordinary stars, and no vetting check reads it yet. Both come before it is applied to a flagged star's tier: thresholds are fixed on ordinary
stars first.

## Decision rule for the first flagged star

Written on 2026-11-18 at 15:55 UTC, before any SPHEREx measurement of Gaia DR3 5628404318169472 or of the
comparison stars had been reported. It is not changed after the numbers arrive; if it turns out to be unusable,
that is said, and a new rule is tested on stars that played no part in writing it.

The question: WISE and 2MASS give this star Ks - W2 and W1 - W2 colours about 0.4 mag redder than its spectral
type explains. Does a second instrument, in 2025 and 2026, measure the same colours?

1. Comparison set: the first twelve stars by source id of the local sample (run `v170-local-td`) that no channel
   flagged, with W1 between 10.0 and 11.2 and BP - RP between 1.15 and 1.55. The flagged star has BP - RP 1.354.
   Gaia DR3 9281425163264 (BP - RP 0.77, the star the method was first tried on) is measured again but is not
   part of the set.
2. For every star and each colour c in (k_w2, w1_w2): D_c = the SPHEREx colour of the summary over all visits
   minus the same colour from the catalogue magnitudes (2MASS Ks, AllWISE W1 and W2). D_c holds the difference
   between window means and catalogue bands, which stars of one spectral type share, and the noise.
3. Over the comparison set: m_c = median of D_c, s_c = 1.4826 x the median absolute deviation of D_c.
4. For the flagged star: z_c = (D_c - m_c) / sqrt(s_c^2 + e_c^2), where e_c is the error of its SPHEREx colour.
5. The rule gives a verdict only when: at least 8 comparison stars have all three windows measured and a complete
   run; s_c is at most 0.10 mag for both colours; the flagged star has at least 5 measurements in each window.
   Otherwise: no verdict, and the reason.
6. Verdict:
   - both |z_c| <= 3: **reproduced**. SPHEREx measures the colours WISE and 2MASS measured, so the red colours
     are not a fault of one instrument's photometry of this star.
   - z_c < -3 for either colour: **not reproduced**. The star is bluer in SPHEREx than in the catalogues by more
     than the comparison stars allow. The catalogue colour then needs an explanation (a change since 2010 to
     2024, or contamination that one instrument sees and the other does not) before it counts as an excess.
   - z_c > 3 for either colour and none below -3: **redder now**. Reported as a change, not as a confirmation.
7. Whatever the verdict, it speaks only for wavelengths up to 5 micron and for blended light within one 6 arcsec
   pixel. It does not test the 12 and 22 micron excess, and it cannot separate a neighbour WISE could not.

The limit of 3 and the scatter ceiling of 0.10 mag are conventions chosen here, not fitted to anything. When this
rule becomes a vetting check, both go into the config with a version entry.

## The verdict for the first flagged star (2026-11-18)

Run `spherex-validation-local` measured the twelve comparison stars (8,061 s; 8 with a complete run, 4 short of one or
two cutouts the archive refused). Over all twelve: D(k_w2) has median +0.073 and scatter 0.017, D(w1_w2) has median
+0.110 and scatter 0.026. The flagged star has D(k_w2) = +0.409 and D(w1_w2) = +0.136, so z(k_w2) = +19.6 and
z(w1_w2) = +1.0 (over the 8 complete stars only: +18.2 and +0.8). The conditions of step 5 hold. Verdict:
**redder now**. Against the comparison stars the star is 0.38 mag brighter in the W1 part and 0.41 mag brighter in the
W2 part than AllWISE measured in 2010, and 0.07 mag brighter in K than 2MASS. The numbers per star are in
`reports/5628404318169472_report.md`.

## The vetting check `spherex_colours` (since config 1.10.0)

The rule above is now a check of the infrared channels (`oddstar.vet.ir_checks.check_spherex_colours`), with two
changes that make it usable for any star. Both were made after the first verdict and do not change it.

- **Comparison stars** are no longer taken from a run's sample. They are the nearest stars on the sky, within
  `spherex_comparison_radius_deg` (3.0), whose BP - RP is within `spherex_comparison_bp_rp_halfwidth_mag` (0.2) and
  whose W1 is within `spherex_comparison_w1_halfwidth_mag` (0.6) of the star's, with clean AllWISE W1 and W2 (quality
  A, no contamination flag, not extended) and clean 2MASS photometry: the nearest `spherex_n_comparison_stars` (12).
  Nearby stars were observed by SPHEREx in the same weeks through the same zodiacal light. Nothing about a star's
  infrared colour enters the choice. The archive is asked in two steps, because one query joining Gaia, AllWISE and
  2MASS over a 3 degree cone timed out on the first live run (HTTP 408): first the Gaia sources of the right colour
  in a cone (a quarter of the radius, then half, then all of it, only as far as needed), then the AllWISE and 2MASS
  records of those sources, nearest first, 500 at a time, until twelve clean stars are in hand. Every source nearer
  than the last one read has then been looked at, so the stars are the same as one query of the whole cone would
  give. A Gaia list that comes back at the archive's 50,000 row limit is refused, since it is not the whole cone.
- **Which comparison stars count.** Only those with a complete run and all three windows enter the median and the
  scatter, and at least `spherex_min_comparison_stars` (8) of them are needed.

The limits are `spherex_z_limit` (3.0), `spherex_max_scatter_mag` (0.10) and `spherex_min_points_per_window` (5).

| Verdict | Check status | Meaning |
| --- | --- | --- |
| reproduced | pass | SPHEREx measures the 2 to 5 micron colours the catalogues give |
| redder now | pass | SPHEREx measures the star redder: the catalogue colour is not a fault of one instrument, and the note says this is a change, not a confirmation |
| not reproduced | fail | SPHEREx measures the star bluer than the catalogues by more than the comparison stars allow |
| none | not_evaluated | no stored measurement, too few usable comparison stars, too much scatter, or too few measurements of the star |

The check answers one question: is the catalogue's 2 to 5 micron colour a fault of one instrument's photometry of
this star? It does not test a 12 or 22 micron excess (`second_instrument` still asks for that), and it speaks for
the blended light within one 6 arcsec pixel.

Vetting does not fetch SPHEREx images itself: a star and twelve comparison stars are about 3,000 cutouts and two
hours. `oddstar spherex verdict` measures them and stores the colours; `vet --live` and `inspect --live` read the
stored record (`vetting_archive.spherex_colours`, query status `spherex_colours`: "ok" or "not_measured").

## Commands

```
oddstar spherex verdict --ids <gaia_id>[,<gaia_id>] [--as-of 2026-11-18] [--passes 2] [--run-id <name>]
```

Measures each star and its comparison stars on detectors 3 to 6, applies the rule, stores the colours for vetting
and writes `runs/<run-id>/spherex_verdict.json` (default run id `spherex-verdict-<first id>`). Stored rows are not
fetched again, so a repeat after a refused cutout fetches only what is missing. The command does that itself: a
star with images missing is gone through again (`--passes`, 2 by default), because the rule counts a comparison
star only when its run is complete and the archive refuses cutouts when it is busy. On the first live run
(2026-11-18, evening) the archive took about 9 s per cutout and 30 minutes per star, against 8 minutes in the
morning, and every one of the first five comparison stars was short of some images after one pass. Give a later
run of the same star the same `--as-of`, so that it works from the same image list.


```
oddstar spherex photometry --from-run <run> [--ids a,b] [--min-tier C] [--detectors 3,4,5,6] [--as-of 2026-11-18]
```

Writes `runs/spherex-<run>/spherex/<source_id>.json` (every measurement) and `spherex_photometry.json` (the
summaries). About 370 cutouts and 1.5 GB of transfer per star, two at a time.

A star's position comes from the run: from `sample_stars` of its search output, or from `live_stars.json` for a
flagged star. Search outputs written before 2026-11-18 hold a position for flagged stars only; for any other star
named with `--ids` the position, proper motion and G magnitude are looked up in Gaia by source id (each star's
record says which, in `position_from`). `oddstar radio followup` finds positions the same way.
