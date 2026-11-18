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

## Command

```
oddstar spherex photometry --from-run <run> [--ids a,b] [--min-tier C] [--detectors 3,4,5,6] [--as-of 2026-11-18]
```

Writes `runs/spherex-<run>/spherex/<source_id>.json` (every measurement) and `spherex_photometry.json` (the
summaries). About 370 cutouts and 1.5 GB of transfer per star, two at a time.

A star's position comes from the run: from `sample_stars` of its search output, or from `live_stars.json` for a
flagged star. Search outputs written before 2026-11-18 hold a position for flagged stars only; for any other star
named with `--ids` the position, proper motion and G magnitude are looked up in Gaia by source id (each star's
record says which, in `position_from`). `oddstar radio followup` finds positions the same way.
