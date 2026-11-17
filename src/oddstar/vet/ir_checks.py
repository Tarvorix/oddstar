"""All 15 infrared vetting checks plus vision cutout rubric (Section 6).

Input semantics (R9: `not_evaluated` never counts as `pass`):

- Live stars (built by `oddstar.io.live_survey._row_to_star_dict`, recognised by the `live_enriched` key) use only
  inputs written by an archive query whose outcome `vetting_archive.query_status[<name>]` is "ok". A missing
  input or a failed query gives `not_evaluated`, and the note names the query and its recorded status.
- Synthetic fixture stars (`synthetic_fixture: True`) and hand-built dicts carry explicit inputs, which are read as
  evaluated values. An absent or None input still means "not queried" and gives `not_evaluated`.
- Fixture-only semantics (no real sky position or pixels): `allwise.centroid_offset_arcsec` stands for the W4
  image centroid, top-level `ebv_sfd` for the line-of-sight dust, and only fixtures may use the synthetic cutout
  renderer.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.dossier.plots import render_real_cutouts_panel, render_synthetic_fixture_cutouts
from oddstar.io.physics import (
    BAND_EXTINCTION_AV,
    angular_separation_arcsec,
    flux_mjy_to_mag,
    propagate_position,
)
from oddstar.io.vetting_queries import target_excess_mjy
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result

PC_CM = 3.0856775814913673e18


# ----------------------------------------------------------------------------------------------- helpers


def _is_live(star: dict[str, Any]) -> bool:
    return "live_enriched" in star


def _is_fixture(star: dict[str, Any]) -> bool:
    return bool(star.get("synthetic_fixture"))


def _va(star: dict[str, Any]) -> dict[str, Any]:
    return star.get("vetting_archive") or {}


def _qs(star: dict[str, Any]) -> dict[str, str]:
    return _va(star).get("query_status") or {}


def _status_of(star: dict[str, Any], name: str) -> str:
    return _qs(star).get(name) or "not run"


def _query_ok(star: dict[str, Any], name: str) -> bool:
    return _qs(star).get(name) == "ok"


def _flat_input(star: dict[str, Any], key: str, query: str) -> Any:
    """Explicit `vetting_archive[key]`; for live stars only when `query` succeeded."""
    if _is_live(star) and not _query_ok(star, query):
        return None
    return _va(star).get(key)


def _not_evaluated(
    out_dir: Path,
    check: str,
    reason: str,
    star: dict[str, Any],
    queries: list[str],
    metrics: dict[str, Any] | None = None,
) -> VettingResult:
    m = dict(metrics or {})
    if queries:
        m["query_status"] = {q: _status_of(star, q) for q in queries}
        reason = f"{reason} [queries: " + "; ".join(f"{q}={_status_of(star, q)}" for q in queries) + "]"
    return make_vet_result(out_dir, check=check, status="not_evaluated", note=reason, metrics=m)


def _cutout_bands(star: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return ((_va(star).get("cutouts") or {}).get("bands")) or {}


def _band_ok(m: dict[str, Any] | None) -> bool:
    return bool(m) and m.get("status") == "ok" and m.get("offset_arcsec") is not None


# ---------------------------------------------------------------------------------------------- checks


def check_xmatch_integrity(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    aw = star.get("allwise") or {}
    n_neigh = aw.get("n_neighbours_gaia_xmatch")
    if aw.get("ra_deg") is None or aw.get("dec_deg") is None or n_neigh is None:
        return _not_evaluated(
            out_dir,
            "xmatch_integrity",
            "AllWISE position or Gaia-AllWISE best-neighbour count not available",
            star,
            [],
            {"allwise_ra_deg": aw.get("ra_deg"), "gaia_xmatch_count": n_neigh},
        )
    prop_ra, prop_dec = propagate_position(
        float(star["ra_deg"]),
        float(star["dec_deg"]),
        float(star.get("pmra_mas_yr", 0.0)),
        float(star.get("pmdec_mas_yr", 0.0)),
        2016.0,
        2010.2,
    )
    sep = angular_separation_arcsec(prop_ra, prop_dec, float(aw["ra_deg"]), float(aw["dec_deg"]))
    ok = (int(n_neigh) == 1) and (sep <= float(cfg.vetting["xmatch_max_sep_arcsec"]))
    return make_vet_result(
        out_dir,
        check="xmatch_integrity",
        status="pass" if ok else "fail",
        note=f"Gaia-to-WISE propagated sep={sep:.3f} arcsec, neighbours={int(n_neigh)}",
        metrics={"propagated_sep_arcsec": sep, "gaia_xmatch_count": int(n_neigh)},
    )


def check_wise_quality(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    aw = star.get("allwise") or {}
    if aw.get("cc_flags") is None or aw.get("ph_qual") is None:
        return _not_evaluated(out_dir, "wise_quality", "AllWISE cc_flags/ph_qual not available", star, ["allwise_coadd"])
    cc = str(aw["cc_flags"])
    ph_qual = str(aw["ph_qual"])
    ext_flg = int(aw.get("ext_flg") or 0)
    w4rchi2 = aw.get("w4rchi2")
    w4m = aw.get("w4m")
    w4nm = aw.get("w4nm")
    # Fail on artifact/halo flags in any WISE band (W1/W2 halos corrupt local background),
    # extended flag > 1, poor quality U/X, bad PSF chi2, or spurious W4 coadds with <15% frame detections.
    min_frames = int(cfg.vetting["wise_w4_min_frames"])
    bad_cc = any(c in "DHOPdhop" for c in cc)
    bad_qual = len(ph_qual) >= 4 and ph_qual[3] in {"U", "X"}
    low_frame_det = (
        w4nm is not None
        and w4m is not None
        and int(w4m) >= min_frames
        and len(ph_qual) >= 4
        and ph_qual[3] in {"A", "B", "C"}
        and (int(w4nm) / max(1, int(w4m))) < float(cfg.vetting["wise_w4_min_detected_frame_frac"])
    )
    bad_chi2 = w4rchi2 is not None and float(w4rchi2) >= float(cfg.vetting["wise_w4_max_rchi2"])
    few_frames = w4m is not None and int(w4m) < min_frames
    metrics = {
        "cc_flags": cc,
        "ext_flg": ext_flg,
        "ph_qual": ph_qual,
        "w4rchi2": None if w4rchi2 is None else float(w4rchi2),
        "w4m": None if w4m is None else int(w4m),
        "w4nm": None if w4nm is None else int(w4nm),
        "low_w4_frame_detection_fraction": low_frame_det,
    }
    note = f"cc_flags={cc}, ext_flg={ext_flg}, ph_qual={ph_qual}"
    note += f", w4rchi2={float(w4rchi2):.2f}" if w4rchi2 is not None else ", w4rchi2=n/a"
    note += f", w4m={int(w4m)}" if w4m is not None else ", w4m=n/a"
    if w4nm is not None and w4m is not None:
        note += f", w4nm={int(w4nm)}/{int(w4m)}"
    if bad_cc or bad_qual or ext_flg > int(cfg.vetting["wise_max_ext_flg"]) or bad_chi2 or few_frames or low_frame_det:
        return make_vet_result(out_dir, check="wise_quality", status="fail", note=note, metrics=metrics)
    if w4rchi2 is None or w4m is None:
        return _not_evaluated(out_dir, "wise_quality", f"{note}; W4 profile-fit chi2/frame count not available", star, ["allwise_coadd"], metrics)
    return make_vet_result(out_dir, check="wise_quality", status="pass", note=note, metrics=metrics)


def _centroid_verdicts(star: dict[str, Any], cfg: PipelineConfig) -> dict[str, dict[str, Any]]:
    """Per-band (W3, W4) image-centroid verdicts: pass / fail / non_constraining."""
    k = float(cfg.vetting["vision_centroid_max_sigma"])
    frac = float(cfg.vetting["centroid_offset_max_beam_frac"])
    out: dict[str, dict[str, Any]] = {}
    for band, beam_key in (("W3", "w3_beam_arcsec"), ("W4", "w4_beam_arcsec")):
        m = _cutout_bands(star).get(band)
        if not _band_ok(m):
            continue
        beam = float(cfg.vetting[beam_key])
        off = float(m["offset_arcsec"])
        sig = m.get("offset_sigma_arcsec")
        limit = frac * beam
        if sig is None or k * float(sig) > limit:
            verdict = "non_constraining"
        elif off > limit:
            verdict = "fail"
        else:
            verdict = "pass"
        out[band] = {
            "offset_arcsec": off,
            "offset_sigma_arcsec": None if sig is None else float(sig),
            "beam_fraction": off / beam,
            "limit_arcsec": limit,
            "k_sigma_arcsec": None if sig is None else k * float(sig),
            "verdict": verdict,
        }
    return out


def check_centroid_offset(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """W3/W4 image centroid offset from the Gaia position propagated to the WISE epoch.

    The AllWISE catalogue position is dominated by W1/W2 and cannot answer whether the 12-22 um emission is
    centred on the star, so a catalogue-only offset is reported but leaves this check `not_evaluated`.
    """
    aw = star.get("allwise") or {}
    beam_w4 = float(cfg.vetting["w4_beam_arcsec"])
    max_frac = float(cfg.vetting["centroid_offset_max_beam_frac"])
    verdicts = _centroid_verdicts(star, cfg)
    if verdicts:
        metrics = {"source": "W3/W4 Atlas image centroids", "bands": verdicts, "limit_fraction": max_frac,
                   "catalogue_offset_arcsec": aw.get("centroid_offset_arcsec")}
        parts = [
            f"{b} offset={v['offset_arcsec']:.2f}\"+/-{(v['offset_sigma_arcsec'] or float('nan')):.2f}\" "
            f"(limit {v['limit_arcsec']:.2f}\", {v['verdict']})"
            for b, v in verdicts.items()
        ]
        failed = [b for b, v in verdicts.items() if v["verdict"] == "fail"]
        if failed:
            return make_vet_result(out_dir, check="centroid_offset", status="fail",
                                   note=f"Image centroid offset beyond {max_frac:.2f} beam in {failed}: " + "; ".join(parts), metrics=metrics)
        if verdicts.get("W4", {}).get("verdict") == "pass":
            return make_vet_result(out_dir, check="centroid_offset", status="pass",
                                   note="W4 image centroid consistent with the propagated Gaia position: " + "; ".join(parts), metrics=metrics)
        return _not_evaluated(out_dir, "centroid_offset",
                              "W4 image centroid not measured or too uncertain to bound the offset: " + "; ".join(parts),
                              star, ["cutouts_wise"], metrics)
    if _is_fixture(star) and aw.get("centroid_offset_arcsec") is not None:
        offset = float(aw["centroid_offset_arcsec"])
        frac = offset / max(1e-6, beam_w4)
        ok = frac <= max_frac
        return make_vet_result(
            out_dir,
            check="centroid_offset",
            status="pass" if ok else "fail",
            note=f"WISE W4 centroid offset={offset:.2f} arcsec ({frac:.3f} of {beam_w4:.1f}\" beam) [synthetic fixture value]",
            metrics={"centroid_offset_arcsec": offset, "beam_fraction": frac, "limit_fraction": max_frac, "source": "synthetic_fixture"},
        )
    cat_off = aw.get("centroid_offset_arcsec")
    metrics = {"source": "allwise_catalogue_position", "catalogue_offset_arcsec": cat_off,
               "catalogue_beam_fraction": None if cat_off is None else float(cat_off) / beam_w4, "limit_fraction": max_frac}
    cat_txt = f"AllWISE catalogue position offset={float(cat_off):.2f}\" ({float(cat_off) / beam_w4:.3f} beam)" if cat_off is not None else "no catalogue offset"
    return _not_evaluated(
        out_dir,
        "centroid_offset",
        f"No W3/W4 image centroid; {cat_txt} is a W1/W2-dominated catalogue-position test and does not locate the 22 um emission",
        star,
        ["cutouts_wise"],
        metrics,
    )


def check_neighbours(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    nbs = _flat_input(star, "neighbours_in_w4_beam", "gaia_neighbours")
    if nbs is None:
        return _not_evaluated(out_dir, "neighbours", "Gaia neighbours inside the W4 beam not queried", star, ["gaia_neighbours"])
    max_flux_frac = float(cfg.vetting["neighbour_blend_flux_frac"])
    known = [float(n["flux_frac"]) for n in nbs if n.get("flux_frac") is not None]
    n_unknown = sum(1 for n in nbs if n.get("flux_frac") is None)
    blend_frac = float(sum(known))
    metrics = {"neighbour_count": len(nbs), "blend_flux_fraction": blend_frac, "n_without_flux_estimate": n_unknown,
               "limit_fraction": max_flux_frac, "neighbours": nbs}
    basis = {str(n.get("flux_frac_basis")) for n in nbs if n.get("flux_frac_basis")}
    basis_txt = f" (flux basis: {', '.join(sorted(basis))})" if basis else ""
    note = f"{len(nbs)} catalogued neighbours inside W4 beam contributing {blend_frac:.2%} of the target W3/W4 flux{basis_txt}"
    if blend_frac > max_flux_frac:
        return make_vet_result(out_dir, check="neighbours", status="fail", note=note, metrics=metrics)
    if n_unknown:
        return _not_evaluated(out_dir, "neighbours", f"{note}; {n_unknown} neighbour(s) without a flux estimate", star, ["gaia_neighbours"], metrics)
    return make_vet_result(out_dir, check="neighbours", status="pass", note=note, metrics=metrics)


def check_radio_counterpart(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    rad = float(cfg.vetting["radio_search_radius_arcsec"])
    srcs = _flat_input(star, "radio_sources_within_10arcsec", "radio_nvss")
    if srcs is None:
        return _not_evaluated(out_dir, "radio_counterpart", "Radio continuum catalogue not queried", star, ["radio_nvss"])
    hits = [r for r in srcs if float(r.get("sep_arcsec", 99.0)) <= rad]
    ok = len(hits) == 0
    return make_vet_result(
        out_dir,
        check="radio_counterpart",
        status="pass" if ok else "fail",
        note=(
            f"No radio continuum counterparts within {rad:.1f} arcsec (remember <10% of dusty galaxies are radio-loud)"
            if ok
            else f"Radio counterpart detected at {hits[0]['sep_arcsec']:.2f} arcsec ({hits[0].get('survey', 'radio')})"
        ),
        metrics={"search_radius_arcsec": rad, "radio_hits": hits},
    )


def check_deep_imaging(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Sub-arcsecond imaging test for a background galaxy inside the W4 beam (lesson L001, skill wise_w4_vetting).

    Seeing-limited catalogues (PS1, Legacy Survey) are a neighbour screen only: they can show a contaminant
    (`fail`) but never clear the position. `pass` needs an actual assessment of HST/JWST/AO/deep imaging.
    """
    va = _va(star)
    if not _is_live(star):
        if "deep_imaging_available" not in va and "deep_imaging_offset_galaxy" not in va:
            return _not_evaluated(out_dir, "deep_imaging", "No deep imaging input", star, [])
        gal = va.get("deep_imaging_offset_galaxy")
        if gal is not None:
            return make_vet_result(
                out_dir, check="deep_imaging", status="fail",
                note=f"Resolved background galaxy in deep imaging at sep={gal.get('sep_arcsec')} arcsec (r={gal.get('r_mag')})",
                metrics={"deep_imaging_available": True, "offset_galaxy": gal},
            )
        if not bool(va.get("deep_imaging_available", False)):
            return make_vet_result(
                out_dir, check="deep_imaging", status="not_evaluated",
                note="No deep sub-arcsec optical/NIR imaging available at this sky position",
                metrics={"deep_imaging_available": False},
            )
        return make_vet_result(
            out_dir, check="deep_imaging", status="pass",
            note="Deep optical/NIR imaging shows clean isolated point source" + (" [synthetic fixture value]" if _is_fixture(star) else ""),
            metrics={"deep_imaging_available": True, "offset_galaxy": None},
        )

    queries = ["deep_optical_ps1", "deep_optical_ls_dr10", "subarcsec_imaging_mast"]
    screen = va.get("neighbour_screen") or {}
    sub = va.get("subarcsec_imaging") if _query_ok(star, "subarcsec_imaging_mast") else None
    contam = screen.get("contaminant")
    metrics = {
        "screened_surveys": screen.get("surveys", []),
        "contaminant": contam,
        "extended_neighbours": screen.get("extended_neighbours", []),
        "faint_point_sources_1p8_8arcsec": screen.get("faint_point_sources", []),
        "subarcsec_observations": (sub or {}).get("observations", []),
        "subarcsec_assessment": (sub or {}).get("assessment"),
    }
    if contam is not None:
        return make_vet_result(
            out_dir, check="deep_imaging", status="fail",
            note=f"Seeing-limited screen found a likely contaminant: {contam.get('survey')} source at {contam.get('sep_arcsec')}\" ({contam.get('reason')})",
            metrics=metrics,
        )
    assessment = (sub or {}).get("assessment")
    if isinstance(assessment, dict) and assessment.get("assessed"):
        if assessment.get("contaminant"):
            return make_vet_result(
                out_dir, check="deep_imaging", status="fail",
                note=f"Sub-arcsec imaging ({assessment.get('observation')}) shows a contaminant: {assessment.get('contaminant')}",
                metrics=metrics,
            )
        return make_vet_result(
            out_dir, check="deep_imaging", status="pass",
            note=f"Sub-arcsec imaging ({assessment.get('observation')}) assessed by {assessment.get('assessed_by')}: isolated point source inside the W4 beam",
            metrics=metrics,
        )
    n_faint = len(screen.get("faint_point_sources", []))
    n_ext = len(screen.get("extended_neighbours", []))
    screened = ", ".join(screen.get("surveys", [])) or "nothing"
    if sub and sub.get("n_observations", 0) > 0:
        obs = sub.get("observations", [])
        listing = "; ".join(f"{o.get('collection')}/{o.get('instrument')} {o.get('filters')} {o.get('obs_id')}" for o in obs[:6])
        return _not_evaluated(
            out_dir, "deep_imaging",
            f"{sub['n_observations']} HST/JWST image observation(s) cover the position and must be inspected ({listing}); "
            f"seeing-limited screen ({screened}) found no contaminant but cannot clear the beam",
            star, queries, metrics,
        )
    return _not_evaluated(
        out_dir, "deep_imaging",
        f"No sub-arcsec imaging assessed; seeing-limited screen ({screened}) found no contaminant "
        f"({n_ext} extended neighbour(s) recorded without comparable mid-IR flux, {n_faint} faint point source(s) 1.8-8\" recorded) "
        f"but cannot clear the beam (L001)",
        star, queries, metrics,
    )


def _poisson_upper_95(n: int) -> float:
    from scipy.stats import chi2

    return float(chi2.ppf(0.95, 2 * (n + 1)) / 2.0)


def check_chance_alignment(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Probability that an unrelated W3/W4 source bright enough to mimic the excess lies within the confusion radius.

    P = 1 - exp(-n * pi * r^2), n = local surface density of AllWISE sources with S/N >= min_snr and flux >=
    `chance_alignment_min_flux_frac_of_excess` x the target's excess (15"-300" annulus by default). The radius r is
    what the centroid test cannot exclude: max(k sigma_centroid, frac x beam) when the band's image centroid was
    measured and constraining, otherwise half the beam FWHM.
    """
    v = cfg.vetting
    p_max = float(v["chance_alignment_p_max"])
    frac_beam = float(v["centroid_offset_max_beam_frac"])
    k = float(v["vision_centroid_max_sigma"])
    va = _va(star)
    if _is_live(star):
        ca = va.get("chance_alignment") if _query_ok(star, "allwise_density") else None
        if ca is None:
            return _not_evaluated(out_dir, "chance_alignment", "Local W3/W4 source density not measured", star, ["allwise_density"])
        r_in, r_out = (float(x) for x in ca["annulus_arcsec"])
        area = math.pi * (r_out**2 - r_in**2)
        min_snr = float(v["chance_alignment_min_snr"])
        f_frac = float(v["chance_alignment_min_flux_frac_of_excess"])
        verdicts = _centroid_verdicts(star, cfg)
        bands: dict[str, Any] = {}
        for band, beam_key in (("W3", "w3_beam_arcsec"), ("W4", "w4_beam_arcsec")):
            exc = target_excess_mjy(star, band)
            if exc is None or exc <= 0.0:
                bands[band] = {"evaluated": False, "reason": "no positive excess estimate"}
                continue
            f_min = f_frac * exc
            key_f, key_s = f"{band.lower()}_mjy", f"{band.lower()}snr"
            n = sum(
                1 for s in ca.get("sources", [])
                if s.get(key_s) is not None and float(s[key_s]) >= min_snr and s.get(key_f) is not None and float(s[key_f]) >= f_min
            )
            beam = float(v[beam_key])
            cv = verdicts.get(band)
            if cv is not None and cv["verdict"] != "non_constraining" and cv["offset_sigma_arcsec"] is not None:
                r = max(k * cv["offset_sigma_arcsec"], frac_beam * beam)
                r_basis = "max(k*sigma_centroid, frac*beam) from image centroid"
            else:
                r = beam / 2.0
                r_basis = "half beam FWHM (no constraining image centroid)"
            dens = n / area
            p = 1.0 - math.exp(-dens * math.pi * r * r)
            p_up = 1.0 - math.exp(-(_poisson_upper_95(n) / area) * math.pi * r * r)
            bands[band] = {"evaluated": True, "excess_mjy": exc, "flux_cut_mjy": f_min, "n_sources": n, "area_arcsec2": area,
                           "density_per_arcsec2": dens, "radius_arcsec": r, "radius_basis": r_basis,
                           "chance_alignment_prob": p, "chance_alignment_prob_95up": p_up}
        evaluated = {b: d for b, d in bands.items() if d.get("evaluated")}
        metrics = {"bands": bands, "p_max": p_max, "annulus_arcsec": [r_in, r_out], "min_snr": min_snr,
                   "n_rows_in_cone": ca.get("n_rows_in_cone")}
        if not evaluated:
            return _not_evaluated(out_dir, "chance_alignment", "No positive W3/W4 excess estimate to compare against", star, ["allwise_density"], metrics)
        p_worst = max(d["chance_alignment_prob"] for d in evaluated.values())
        metrics["chance_alignment_prob"] = p_worst
        txt = "; ".join(
            f"{b}: {d['n_sources']} src >= {d['flux_cut_mjy']:.2f} mJy in {d['area_arcsec2']:.0f} arcsec2 "
            f"(n={d['density_per_arcsec2']:.2e}/arcsec2), r={d['radius_arcsec']:.2f}\" -> P={d['chance_alignment_prob']:.4f}"
            for b, d in evaluated.items()
        )
        return make_vet_result(
            out_dir, check="chance_alignment", status="fail" if p_worst > p_max else "pass",
            note=f"Measured AllWISE density: {txt} (limit P <= {p_max})", metrics=metrics,
        )

    density = va.get("local_red_source_density_per_arcsec2")
    if density is None:
        return _not_evaluated(out_dir, "chance_alignment", "No local red-source density input", star, [])
    beam_w4 = float(v["w4_beam_arcsec"])
    r = frac_beam * beam_w4 if _is_fixture(star) else beam_w4 / 2.0
    p_align = 1.0 - math.exp(-float(density) * math.pi * r * r)
    return make_vet_result(
        out_dir, check="chance_alignment", status="fail" if p_align > p_max else "pass",
        note=f"Background red-source chance alignment P={p_align:.4f} within {r:.1f}\" (density {float(density):.1e}/arcsec2, limit {p_max})"
        + (" [synthetic fixture value]" if _is_fixture(star) else ""),
        metrics={"chance_alignment_prob": p_align, "radius_arcsec": r, "local_density_arcsec2": float(density), "p_max": p_max},
    )


def _excess_colours(star: dict[str, Any]) -> dict[str, float | None]:
    mags: dict[str, float | None] = {}
    for band in ("W1", "W2", "W3", "W4"):
        exc = target_excess_mjy(star, band) if band in (star.get("flux_mjy") or {}) else None
        mags[band] = flux_mjy_to_mag(band, exc) if (exc is not None and exc > 0.0) else None
    def _d(a: str, b: str) -> float | None:
        return (mags[a] - mags[b]) if (mags[a] is not None and mags[b] is not None) else None  # type: ignore[operator]
    return {"excess_w1_w2": _d("W1", "W2"), "excess_w2_w3": _d("W2", "W3"), "excess_w3_w4": _d("W3", "W4")}


def check_extragalactic_colours(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    m = star["mags"]
    w1_w2 = float(m["W1"]) - float(m["W2"])
    w2_w3 = float(m["W2"]) - float(m["W3"])
    w1_w2_lim = float(cfg.vetting["extragalactic_w1_w2_min"])
    w2_w3_lim = float(cfg.vetting["extragalactic_w2_w3_min"])
    # Stern/Mateos AGN/ULIRG locus AND high-excess galaxy color check (total-source colours)
    agn_like = w1_w2 >= w1_w2_lim and w2_w3 >= w2_w3_lim
    # If blended with an offset galaxy or DOG, flag extragalactic colors if star is faint & W1-W2 red
    nbs = _flat_input(star, "neighbours_in_w4_beam", "gaia_neighbours")
    dog_colours = w1_w2 >= float(cfg.vetting["dog_w1_w2_min"]) and w2_w3 >= float(cfg.vetting["dog_w2_w3_min"])
    dog_like = dog_colours and nbs is not None and len(nbs) > 0
    exc = _excess_colours(star)
    exc_txt = ", ".join(f"{k.replace('excess_', '').replace('_', '-').upper()}={v:+.2f}" for k, v in exc.items() if v is not None) or "n/a"
    note = (
        f"WISE colors W1-W2={w1_w2:.2f}, W2-W3={w2_w3:.2f}; photosphere-subtracted excess colours {exc_txt} "
        f"(reported only: hot circumstellar dust occupies the same colour region as AGN, so excess colours alone do not discriminate)"
    )
    metrics = {"w1_w2": w1_w2, "w2_w3": w2_w3, "agn_like": agn_like, "dog_like": dog_like, **exc}
    if agn_like or dog_like:
        return make_vet_result(out_dir, check="extragalactic_colours", status="fail", note=note, metrics=metrics)
    if dog_colours and nbs is None:
        return _not_evaluated(out_dir, "extragalactic_colours", f"{note}; DOG-like colours but beam neighbours not queried", star, ["gaia_neighbours"], metrics)
    return make_vet_result(out_dir, check="extragalactic_colours", status="pass", note=note, metrics=metrics)


def check_cirrus(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Line-of-sight dust (IRSA DUST reference pixel, SFD 1998) and measured W3/W4 Atlas background rms.

    Background rms is the point-source-equivalent 1-sigma noise of beam-sized apertures in a blank annulus of the
    Atlas cutout (mJy, see `oddstar.io.cutouts`); AllWISE `w3sigsk/w4sigsk` are not used.
    """
    ebv_max = float(cfg.vetting["cirrus_ebv_max"])
    bg_max = float(cfg.vetting["cirrus_bg_rms_max_mjy"])
    aw = star.get("allwise") or {}
    va = _va(star)
    if _is_fixture(star) and va.get("dust") is None and not _cutout_bands(star):
        ebv = star.get("ebv_sfd")
        bg = aw.get("cirrus_bg_rms_mjy")
        if ebv is None or bg is None:
            return _not_evaluated(out_dir, "cirrus", "Fixture lacks ebv_sfd or background rms", star, [])
        ok = (float(ebv) <= ebv_max) and (float(bg) <= bg_max)
        return make_vet_result(
            out_dir, check="cirrus", status="pass" if ok else "fail",
            note=f"SFD E(B-V)={float(ebv):.3f}, local W3/W4 background RMS={float(bg):.2f} mJy [synthetic fixture values]",
            metrics={"ebv_sfd": float(ebv), "cirrus_bg_rms_mjy": float(bg), "source": "synthetic_fixture"},
        )
    dust = va.get("dust") if (not _is_live(star) or _query_ok(star, "irsa_dust")) else None
    ebv = (dust or {}).get("ebv_sfd_los")
    bg_by_band = {b: m.get("bg_rms_mjy_beam") for b, m in _cutout_bands(star).items() if b in ("W3", "W4") and m.get("bg_rms_mjy_beam") is not None}
    if not bg_by_band and not _is_live(star) and aw.get("cirrus_bg_rms_mjy") is not None:
        bg_by_band = {"W3/W4": aw["cirrus_bg_rms_mjy"]}
    bg = max(float(x) for x in bg_by_band.values()) if bg_by_band else None
    metrics = {
        "ebv_sfd_los": ebv,
        "ebv_sf11_los": (dust or {}).get("ebv_sf11_los"),
        "i100_mjy_sr": (dust or {}).get("i100_mjy_sr"),
        "bg_rms_mjy_beam": bg_by_band,
        "ebv_max": ebv_max,
        "bg_rms_max_mjy": bg_max,
        "stellar_ebv_not_used": star.get("ebv_sfd"),
    }
    parts = []
    if ebv is not None:
        parts.append(f"line-of-sight SFD E(B-V)={float(ebv):.4f} (S&F {metrics['ebv_sf11_los']}, I100={metrics['i100_mjy_sr']} MJy/sr)")
    if bg is not None:
        parts.append("Atlas background rms " + ", ".join(f"{b}={float(x):.2f}" for b, x in bg_by_band.items()) + " mJy/beam")
    note = "; ".join(parts)
    queries = ["irsa_dust", "cutouts_wise"]
    if ebv is None and bg is None:
        return _not_evaluated(out_dir, "cirrus", "Neither line-of-sight dust nor a measured W3/W4 background is available", star, queries, metrics)
    if (ebv is not None and float(ebv) > ebv_max) or (bg is not None and bg > bg_max):
        return make_vet_result(out_dir, check="cirrus", status="fail", note=f"{note} exceeds limits (E(B-V) <= {ebv_max}, rms <= {bg_max} mJy)", metrics=metrics)
    if ebv is None or bg is None:
        missing = "line-of-sight dust" if ebv is None else "measured W3/W4 background"
        return _not_evaluated(out_dir, "cirrus", f"{note} within limits, but {missing} not available", star, queries, metrics)
    return make_vet_result(out_dir, check="cirrus", status="pass", note=f"{note} within limits (E(B-V) <= {ebv_max}, rms <= {bg_max} mJy)", metrics=metrics)


def check_single_frames(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    frames = list(star.get("wise_single_frames_w4_mjy") or [])
    if not frames:
        return _not_evaluated(out_dir, "single_frames", "No 2010 single-exposure WISE frames available", star, ["wise_l1b"] if _is_live(star) else [])
    arr = np.asarray(frames, dtype=float)
    med = float(np.median(arr))
    consistent = int(np.sum(arr >= float(cfg.vetting["single_frame_consistency_frac"])
                            * max(med, float(star["flux_mjy"]["W4"]) * float(cfg.vetting["single_frame_coadd_floor_frac"]))))
    frac = consistent / len(arr)
    min_frac = float(cfg.vetting["single_frame_min_fraction"])
    # Also check whether one single outlier spike dominates the mean over median
    max_to_med = float(np.max(arr) / max(0.05, med))
    ok = (frac >= min_frac) and (max_to_med < float(cfg.vetting["single_frame_max_to_median"]))
    return make_vet_result(
        out_dir,
        check="single_frames",
        status="pass" if ok else "fail",
        note=f"{consistent}/{len(arr)} single frames consistent ({frac:.2%}), max/median={max_to_med:.2f}",
        metrics={"n_frames": len(arr), "consistent_fraction": frac, "max_to_median": max_to_med},
    )


def check_moving_object(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    va = _va(star)
    r_match = float(cfg.vetting["sso_match_radius_arcsec"])
    if _is_live(star):
        sso = va.get("sso_search")
        queries = ["wise_l1b", "skybot"]
        if sso is None:
            return _not_evaluated(out_dir, "moving_object", "SkyBoT ephemeris search not run", star, queries)
        metrics = {"frame_mjds": sso.get("frame_mjds", []), "visits": sso.get("visits", []), "matches": sso.get("matches", []),
                   "match_radius_arcsec": r_match, "observer": sso.get("observer")}
        if sso.get("matches"):
            mt = sso["matches"][0]
            return make_vet_result(
                out_dir, check="moving_object", status="fail",
                note=f"Known Solar System body {mt['name']} ({mt.get('class')}, V={mt.get('v_mag')}) at {mt['sep_arcsec']}\" "
                f"at frame MJD {mt['frame_mjd']:.4f} ({len(sso['matches'])} frame match(es))",
                metrics=metrics,
            )
        if not sso.get("frame_mjds"):
            return _not_evaluated(out_dir, "moving_object", "No 2010 single-exposure frame epochs known", star, queries, metrics)
        bad = [vv for vv in sso.get("visits", []) if vv.get("status") != "ok"]
        if bad or not sso.get("visits"):
            return _not_evaluated(out_dir, "moving_object",
                                  f"{len(bad)} of {len(sso.get('visits', []))} visit ephemeris queries not usable: {[b.get('status') for b in bad][:3]}",
                                  star, queries, metrics)
        n_frames = len(sso["frame_mjds"])
        return make_vet_result(
            out_dir, check="moving_object", status="pass",
            note=f"No known Solar System body within {r_match:.0f}\" at any of {n_frames} 2010 frame epochs "
            f"({len(sso['visits'])} SkyBoT visit queries, observer {sso.get('observer')})",
            metrics=metrics,
        )
    mov = va.get("moving_objects_2010")
    if mov is None:
        return _not_evaluated(out_dir, "moving_object", "No ephemeris search input", star, [])
    ok = len(mov) == 0
    return make_vet_result(
        out_dir,
        check="moving_object",
        status="pass" if ok else "fail",
        note="No known Solar System bodies in ephemeris cone at 2010 WISE epochs" if ok else f"Asteroid match: {mov}",
        metrics={"moving_object_hits": mov},
    )


def check_second_instrument(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    sec = _flat_input(star, "second_instrument_ir", "second_instrument")
    if sec is None:
        return _not_evaluated(
            out_dir, "second_instrument",
            "No archival Spitzer, AKARI, IRAS, Herschel, or SPHEREx 22um-overlapping photometry at this position",
            star, ["second_instrument"] if _is_live(star) else [], {"second_instrument": None},
        )
    interp = sec.get("interpretation", "confirmed")
    if interp == "catalogue_match_unverified":
        return _not_evaluated(
            out_dir, "second_instrument",
            f"{sec.get('n_matches')} {sec.get('instrument')} catalogue row(s) within 10\"; photometry not yet compared with W4 (inspect)",
            star, ["second_instrument"], {"second_instrument": sec},
        )
    if interp == "no_excess_at_second_instrument":
        # Photometry at a comparable wavelength, deep enough to detect the W4 excess, sees only the photosphere.
        return make_vet_result(
            out_dir, check="second_instrument", status="fail",
            note=f"{sec.get('instrument')} photometry shows no excess where the W4 excess predicts a detection",
            metrics={"second_instrument": sec},
        )
    # Melis et al. 2012 case: instrument disagreement can be real secular variability
    return make_vet_result(
        out_dir,
        check="second_instrument",
        status="pass",
        note=f"Second-instrument archive {sec.get('instrument')} evaluated ({interp})",
        metrics={"second_instrument": sec},
    )


# Approximate inactive main-sequence locus of GALEX NUV (AB) - Gaia G (Vega) versus dereddened Gaia BP-RP.
# NOT taken from a published fit: anchored on a solar NUV-G ~= 5.4 (estimated from the solar UV spectral irradiance
# averaged over the NUV band, uncertain by ~0.5 mag), with the colour dependence taken as the Planck NUV-G colour
# difference relative to the Sun, amplified by 1.5 for stars redder than the Sun to mimic stronger NUV line
# blanketing (the factor is an assumption). The adopted uncertainty `youth_nuv_locus_sigma_mag` (0.5 mag) enters
# every comparison. Replace with an empirical GALEX x Gaia inactive-dwarf locus before relying on the UV sub-test.
NUV_G_LOCUS_BP_RP = (0.30, 0.55, 0.82, 1.00, 1.35, 1.80, 2.50)
NUV_G_LOCUS = (3.7, 4.6, 5.4, 6.5, 8.5, 10.9, 13.1)
# Extinction: A_NUV = 8.2 E(B-V) (Wyder et al. 2007, R_V = 3.1); Gaia bands from BAND_EXTINCTION_AV with A_V = 3.1 E(B-V)
R_NUV = 8.2
# ROSAT PSPC count rate -> 0.1-2.4 keV flux for coronal spectra: (5.30 HR1 + 8.31) x 1e-12 erg cm^-2 per count
# (Fleming et al. 1995); HR1 = 0 is assumed when HR1 is not reported.


def _youth_uv(star: dict[str, Any], galex: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    thr = float(cfg.vetting["youth_uv_excess_mag"])
    if galex.get("status") != "ok":
        return {"verdict": "not_queried", "status": galex.get("status")}
    match = galex.get("match")
    if match is None:
        return {"verdict": "non_constraining", "reason": "no GALEX AIS source (coverage gap and faint photosphere are indistinguishable)"}
    mags = star.get("mags") or {}
    if mags.get("G") is None or mags.get("BP") is None or mags.get("RP") is None:
        return {"verdict": "non_constraining", "reason": "no Gaia G/BP/RP photometry"}
    ebv = float(star.get("ebv_sfd") or 0.0)
    av = 3.1 * ebv
    bp_rp0 = float(mags["BP"]) - float(mags["RP"]) - av * (BAND_EXTINCTION_AV["BP"] - BAND_EXTINCTION_AV["RP"])
    g0 = float(mags["G"]) - av * BAND_EXTINCTION_AV["G"]
    nuv0 = float(match["nuv_mag"]) - R_NUV * ebv
    if not (NUV_G_LOCUS_BP_RP[0] <= bp_rp0 <= NUV_G_LOCUS_BP_RP[-1]):
        return {"verdict": "non_constraining", "reason": f"(BP-RP)0={bp_rp0:.2f} outside the locus range", "bp_rp0": bp_rp0}
    nuv_exp = g0 + float(np.interp(bp_rp0, NUV_G_LOCUS_BP_RP, NUV_G_LOCUS))
    excess = nuv_exp - nuv0
    sig = math.hypot(float(match.get("nuv_err") or cfg.vetting["youth_nuv_default_err_mag"]), float(cfg.vetting["youth_nuv_locus_sigma_mag"]))
    nsig = float(cfg.vetting["youth_uv_constraint_nsigma"])
    if excess > thr:
        verdict = "young"
    elif excess + nsig * sig < thr:
        verdict = "not_young"
    else:
        verdict = "non_constraining"
    return {"verdict": verdict, "nuv_mag": match["nuv_mag"], "nuv_err": match.get("nuv_err"), "sep_arcsec": match.get("sep_arcsec"),
            "nuv_expected_photosphere": nuv_exp, "nuv_excess_mag": excess, "excess_sigma_mag": sig, "bp_rp0": bp_rp0, "ebv_used": ebv}


def _youth_xray(star: dict[str, Any], youth: dict[str, Any], cfg: PipelineConfig) -> dict[str, Any]:
    lx_min = float(cfg.vetting["youth_xray_lx_min_erg_s"])
    d_pc = star.get("distance_pc")
    if d_pc is None:
        return {"verdict": "non_constraining", "reason": "no distance"}
    four_pi_d2 = 4.0 * math.pi * (float(d_pc) * PC_CM) ** 2
    surveys: dict[str, Any] = {}
    for key, lim_key in (("rosat_2rxs", "youth_rosat_flux_limit_erg_cm2_s"), ("erass1", "youth_erass1_flux_limit_erg_cm2_s")):
        sv = youth.get(key) or {}
        st = sv.get("status")
        if st is None:
            surveys[key] = {"verdict": "not_queried"}
            continue
        if st == "no_coverage":
            surveys[key] = {"verdict": "no_coverage", "note": sv.get("note")}
            continue
        if st != "ok":
            surveys[key] = {"verdict": "not_queried", "status": st}
            continue
        m = sv.get("match")
        if m is not None:
            if key == "rosat_2rxs":
                hr1 = m.get("hr1")
                ecf = (5.30 * (float(hr1) if hr1 is not None else 0.0) + 8.31) * 1e-12
                flux = float(m["rate_cts_s"]) * ecf
            else:
                flux = float(m["flux_erg_cm2_s"])
            lx = four_pi_d2 * flux
            surveys[key] = {"verdict": "young" if lx > lx_min else "not_young", "detected": True, "flux_erg_cm2_s": flux,
                            "lx_erg_s": lx, "sep_arcsec": m.get("sep_arcsec")}
        else:
            lim = float(cfg.vetting[lim_key])
            lx_lim = four_pi_d2 * lim
            surveys[key] = {"verdict": "not_young" if lx_lim < lx_min else "non_constraining", "detected": False,
                            "flux_limit_erg_cm2_s": lim, "lx_upper_limit_erg_s": lx_lim}
    verdicts = [s["verdict"] for s in surveys.values()]
    if "young" in verdicts:
        verdict = "young"
    elif "not_young" in verdicts:
        verdict = "not_young"
    elif "non_constraining" in verdicts:
        verdict = "non_constraining"
    else:
        verdict = "not_queried"
    return {"verdict": verdict, "lx_threshold_erg_s": lx_min, "surveys": surveys}


def _youth_association(assoc: dict[str, Any]) -> dict[str, Any]:
    if assoc.get("status") != "ok":
        return {"verdict": "not_queried", "status": assoc.get("status")}
    if assoc.get("young_groups") or assoc.get("young_object_type"):
        return {"verdict": "young", "young_groups": assoc.get("young_groups"), "simbad_object_type": assoc.get("simbad_object_type")}
    if assoc.get("clusters"):
        return {"verdict": "non_constraining", "reason": "cluster membership of unknown age", "clusters": assoc.get("clusters")}
    return {"verdict": "none_listed", "reason": "no SIMBAD-listed young group membership (catalogue incompleteness: not evidence of age)"}


def check_youth(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Youth indicators (UV excess, X-ray luminosity, young-group membership).

    `fail` on any positive indicator. `pass` only when BOTH the UV and the X-ray sub-tests are constraining and clean
    (a non-detection whose survey limit is above the youth threshold is not a pass). Association membership can only
    add a `fail`; the absence of a catalogued membership is never used as evidence of age.
    """
    lx_min = float(cfg.vetting["youth_xray_lx_min_erg_s"])
    thr_uv = float(cfg.vetting["youth_uv_excess_mag"])
    va = _va(star)
    if _is_live(star):
        youth = va.get("youth")
        queries = ["galex", "xray_2rxs", "xray_erass1", "young_association"]
        if youth is None:
            return _not_evaluated(out_dir, "youth", "UV/X-ray/association queries not run", star, queries)
        uv = _youth_uv(star, youth.get("galex") or {}, cfg)
        xr = _youth_xray(star, youth, cfg)
        assoc = _youth_association(youth.get("association") or {})
        metrics = {"uv": uv, "xray": xr, "association": assoc}
        summary = f"UV: {uv['verdict']}" + (f" (NUV excess {uv['nuv_excess_mag']:+.2f}+/-{uv['excess_sigma_mag']:.2f} mag)" if "nuv_excess_mag" in uv else "")
        xs = []
        for name, s in xr.get("surveys", {}).items():
            if "lx_erg_s" in s:
                xs.append(f"{name} Lx={s['lx_erg_s']:.1e}")
            elif "lx_upper_limit_erg_s" in s:
                xs.append(f"{name} Lx<{s['lx_upper_limit_erg_s']:.1e}")
            else:
                xs.append(f"{name} {s['verdict']}")
        summary += f"; X-ray: {xr['verdict']} ({', '.join(xs)} vs {lx_min:.0e} erg/s); association: {assoc['verdict']}"
        if "young" in (uv["verdict"], xr["verdict"], assoc["verdict"]):
            return make_vet_result(out_dir, check="youth", status="fail", note=f"Youth indicator present: {summary}", metrics=metrics)
        if uv["verdict"] == "not_young" and xr["verdict"] == "not_young":
            return make_vet_result(out_dir, check="youth", status="pass", note=f"No youth indicators at constraining depth: {summary}", metrics=metrics)
        return _not_evaluated(out_dir, "youth", f"Youth not constrained: {summary}", star, queries, metrics)

    uv_ex = va.get("galex_uv_excess_mag")
    xray = va.get("xray_luminosity_erg_s")
    assoc = va.get("young_association")
    metrics = {"young_association": assoc, "galex_uv_excess_mag": uv_ex, "xray_luminosity_erg_s": xray}
    young = (assoc is not None) or (uv_ex is not None and float(uv_ex) > thr_uv) or (xray is not None and float(xray) > lx_min)
    if young:
        return make_vet_result(
            out_dir, check="youth", status="fail",
            note=f"Youth indicators present (assoc={assoc}, UV={uv_ex}, Lx={xray})", metrics=metrics,
        )
    if uv_ex is None or xray is None:
        return _not_evaluated(out_dir, "youth", "UV excess or X-ray luminosity input missing", star, [], metrics)
    return make_vet_result(
        out_dir, check="youth", status="pass",
        note="No youth indicators (association, UV/X-ray excess)" + (" [synthetic fixture values]" if _is_fixture(star) else ""),
        metrics=metrics,
    )


def check_evolved(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    mg = float(star.get("abs_g_mag", 4.8))
    teff = float(star.get("teff_k", 5800.0))
    bright_lim = float(cfg.vetting["evolved_mg_bright_limit"])
    teff_default = str(star.get("teff_k_source", "")).startswith("default")
    metrics = {"abs_g_mag": mg, "teff_k": teff, "teff_k_source": star.get("teff_k_source")}
    if mg < bright_lim and teff_default:
        return _not_evaluated(out_dir, "evolved", f"M_G={mg:.2f} is giant-bright but Teff is a default (no GSP-Phot Teff)", star, [], metrics)
    is_evolved = (mg < bright_lim) and (teff < float(cfg.vetting["evolved_teff_max_k"]))
    metrics["is_evolved"] = is_evolved
    return make_vet_result(
        out_dir,
        check="evolved",
        status="fail" if is_evolved else "pass",
        note=f"HR position M_G={mg:.2f}, Teff={teff:.0f} K ({'giant branch' if is_evolved else 'main sequence'})",
        metrics=metrics,
    )


def check_binary(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    ruwe = float(star.get("ruwe", 1.0))
    va = _va(star)
    nss = _flat_input(star, "nss_flag", "gaia_meta")
    rv_sig = _flat_input(star, "rv_scatter_kms", "gaia_meta")
    rv_faint_unconstrained = bool(va.get("rv_faint_unconstrained", False))
    ruwe_max = float(cfg.vetting["binary_ruwe_max"])
    rv_max = float(cfg.vetting["binary_rv_scatter_max_kms"])
    metrics = {
        "ruwe": ruwe,
        "nss_flag": None if nss is None else int(nss),
        "rv_scatter_kms": None if rv_sig is None else float(rv_sig),
        "rv_faint_unconstrained": rv_faint_unconstrained,
    }
    nss_txt = "n/a" if nss is None else str(int(nss))
    rv_txt = "n/a" if rv_sig is None else f"{float(rv_sig):.2f} km/s"
    if (ruwe > ruwe_max) or (nss is not None and int(nss) > 0) or (rv_sig is not None and float(rv_sig) > rv_max and not rv_faint_unconstrained):
        return make_vet_result(out_dir, check="binary", status="fail", note=f"RUWE={ruwe:.2f}, NSS={nss_txt}, RV scatter={rv_txt}", metrics=metrics)
    if nss is None:
        return _not_evaluated(out_dir, "binary", f"RUWE={ruwe:.2f} but Gaia NSS/RV metadata not available", star, ["gaia_meta"], metrics)
    if rv_sig is None:
        return _not_evaluated(out_dir, "binary", f"RUWE={ruwe:.2f}, NSS={nss_txt}, no Gaia RV epoch statistic: close SB unconstrained", star, ["gaia_meta"], metrics)
    if float(rv_sig) > rv_max and rv_faint_unconstrained:
        return make_vet_result(
            out_dir, check="binary", status="not_evaluated",
            note=f"Faint-star Gaia RVS (RUWE={ruwe:.2f}, implied RV scatter={float(rv_sig):.2f} km/s without p-value): close SB1 binary unconstrained",
            metrics=metrics,
        )
    return make_vet_result(out_dir, check="binary", status="pass", note=f"RUWE={ruwe:.2f}, NSS={nss_txt}, RV scatter={rv_txt}", metrics=metrics)


def check_vision_rubric(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Real-pixel W3/W4 rubric: centroid offset (k sigma AND beam fraction), extendedness, brighter secondary peak.

    Synthetic fixtures (no real sky pixels) keep the synthetic renderer and their explicit offset/background inputs.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    cutout_path = out_dir / "wise_cutouts.png"
    va = _va(star)
    if _is_fixture(star):
        render_synthetic_fixture_cutouts(star, cutout_path)
        aw = star.get("allwise") or {}
        offset = float(aw.get("centroid_offset_arcsec", 0.0))
        bg_rms = float(aw.get("cirrus_bg_rms_mjy", 0.4))
        has_off_gal = va.get("deep_imaging_offset_galaxy") is not None
        limit = float(cfg.vetting["centroid_offset_max_beam_frac"]) * float(cfg.vetting["w4_beam_arcsec"])
        ok = (offset <= limit) and (bg_rms <= float(cfg.vetting["cirrus_bg_rms_max_mjy"])) and (not has_off_gal)
        return make_vet_result(
            out_dir,
            check="vision_rubric",
            status="pass" if ok else "fail",
            note="SYNTHETIC fixture cutout panel (no real pixels) graded against PSF/offset/cirrus rubric from explicit fixture inputs",
            metrics={"cutout_path": str(cutout_path), "offset_arcsec": offset, "cirrus_bg_rms_mjy": bg_rms, "synthetic_fixture": True},
            evidence_override=cutout_path,
        )

    bands = _cutout_bands(star)
    queries = ["cutouts_wise", "cutouts_ps1"]
    measured = {b: bands[b] for b in ("W3", "W4") if _band_ok(bands.get(b))}
    if not measured:
        return _not_evaluated(out_dir, "vision_rubric", "No real W3/W4 cutout pixels measured", star, queries,
                              {"bands": {b: (bands.get(b) or {}).get("status") for b in ("W3", "W4")}})
    v = cfg.vetting
    k = float(v["vision_centroid_max_sigma"])
    frac = float(v["centroid_offset_max_beam_frac"])
    ext_thr = float(v["vision_extended_fwhm_ratio"])
    min_snr = float(v["vision_min_fit_snr"])
    flags: list[str] = []
    low_snr: list[str] = []
    per_band: dict[str, Any] = {}
    for b, m in measured.items():
        beam = float(v["w3_beam_arcsec"] if b == "W3" else v["w4_beam_arcsec"])
        snr = m.get("flux_snr")
        off, sig = float(m["offset_arcsec"]), m.get("offset_sigma_arcsec")
        if snr is None or float(snr) < min_snr:
            low_snr.append(b)
        if off > frac * beam and sig is not None and off > k * float(sig):
            flags.append(f"{b} centroid offset {off:.2f}\" > {frac:.2f} beam and > {k:.0f} sigma ({float(sig):.2f}\")")
        if snr is not None and float(snr) >= min_snr and m.get("fwhm_ratio") is not None and float(m["fwhm_ratio"]) > ext_thr:
            flags.append(f"{b} extended: FWHM {m['fwhm_arcsec']:.1f}\" = {float(m['fwhm_ratio']):.2f} x reference {m['fwhm_ref_arcsec']:.1f}\"")
        if b == "W4" and m.get("brighter_secondary_peak"):
            flags.append(f"W4 secondary peak brighter than the target peak within the beam ({m.get('n_secondary_peaks')} secondary)")
        per_band[b] = {key: m.get(key) for key in ("offset_arcsec", "offset_sigma_arcsec", "fwhm_arcsec", "fwhm_ref_arcsec", "fwhm_ratio",
                                                    "flux_mjy", "flux_snr", "bg_rms_mjy_beam", "n_secondary_peaks", "brighter_secondary_peak",
                                                    "secondary_peaks", "sha256")}
    png = render_real_cutouts_panel(star, cutout_path)
    metrics = {"bands": per_band, "cutout_png": str(png) if png else None, "low_snr_bands": low_snr,
               "thresholds": {"k_sigma": k, "beam_frac": frac, "extended_fwhm_ratio": ext_thr, "min_fit_snr": min_snr},
               "epoch_jyear": (va.get("cutouts") or {}).get("epoch_jyear")}
    if flags:
        return make_vet_result(out_dir, check="vision_rubric", status="fail", note="Real-pixel rubric: " + "; ".join(flags), metrics=metrics)
    missing = [b for b in ("W3", "W4") if b not in measured]
    if missing or low_snr:
        return _not_evaluated(out_dir, "vision_rubric",
                              f"Real-pixel rubric incomplete: missing {missing}, S/N < {min_snr:.0f} in {low_snr}", star, queries, metrics)
    return make_vet_result(
        out_dir, check="vision_rubric", status="pass",
        note="Real-pixel W3/W4 rubric clean: centroids within limits, point-like FWHM, no brighter secondary peak in the W4 beam",
        metrics=metrics,
    )


def run_ir_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    """Run all 15 infrared checks + vision rubric on a target star."""
    checks = [
        check_xmatch_integrity,
        check_wise_quality,
        check_centroid_offset,
        check_neighbours,
        check_radio_counterpart,
        check_deep_imaging,
        check_chance_alignment,
        check_extragalactic_colours,
        check_cirrus,
        check_single_frames,
        check_moving_object,
        check_second_instrument,
        check_youth,
        check_evolved,
        check_binary,
        check_vision_rubric,
    ]
    return [fn(star, cfg, out_dir) for fn in checks]
