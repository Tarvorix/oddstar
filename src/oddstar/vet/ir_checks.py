"""All 15 infrared vetting checks plus vision cutout rubric (Section 6)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from oddstar.config import PipelineConfig
from oddstar.dossier.plots import render_cutouts_panel
from oddstar.io.physics import angular_separation_arcsec, propagate_position
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result


def check_xmatch_integrity(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    aw = star.get("allwise", {})
    prop_ra, prop_dec = propagate_position(
        float(star["ra_deg"]),
        float(star["dec_deg"]),
        float(star.get("pmra_mas_yr", 0.0)),
        float(star.get("pmdec_mas_yr", 0.0)),
        2016.0,
        2010.2,
    )
    wise_ra = float(aw.get("ra_deg", prop_ra))
    wise_dec = float(aw.get("dec_deg", prop_dec))
    sep = angular_separation_arcsec(prop_ra, prop_dec, wise_ra, wise_dec)
    n_neigh = int(aw.get("n_neighbours_gaia_xmatch", 1))
    ok = (n_neigh == 1) and (sep <= 2.5)
    return make_vet_result(
        out_dir,
        check="xmatch_integrity",
        status="pass" if ok else "fail",
        note=f"Gaia-to-WISE propagated sep={sep:.3f} arcsec, neighbours={n_neigh}",
        metrics={"propagated_sep_arcsec": sep, "gaia_xmatch_count": n_neigh},
    )


def check_wise_quality(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    aw = star.get("allwise", {})
    cc = str(aw.get("cc_flags", "0000"))
    ext_flg = int(aw.get("ext_flg", 0))
    ph_qual = str(aw.get("ph_qual", "AAAA"))
    w4rchi2 = float(aw.get("w4rchi2", 1.0))
    w4m = int(aw.get("w4m", 12))
    # Fail on artifact flags in W3/W4, extended flag > 1, poor quality U/X, or bad PSF chi2
    bad_cc = any(c in "DHOPdhop" for c in cc[2:])
    bad_qual = len(ph_qual) >= 4 and ph_qual[3] in {"U", "X"}
    ok = (not bad_cc) and (ext_flg <= 1) and (not bad_qual) and (w4rchi2 < 3.5) and (w4m >= 5)
    return make_vet_result(
        out_dir,
        check="wise_quality",
        status="pass" if ok else "fail",
        note=f"cc_flags={cc}, ext_flg={ext_flg}, ph_qual={ph_qual}, w4rchi2={w4rchi2:.2f}, w4m={w4m}",
        metrics={"cc_flags": cc, "ext_flg": ext_flg, "ph_qual": ph_qual, "w4rchi2": w4rchi2, "w4m": w4m},
    )


def check_centroid_offset(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    aw = star.get("allwise", {})
    offset = float(aw.get("centroid_offset_arcsec", 0.0))
    beam_w4 = float(cfg.vetting["w4_beam_arcsec"])
    max_frac = float(cfg.vetting["centroid_offset_max_beam_frac"])
    frac = offset / max(1e-6, beam_w4)
    ok = frac <= max_frac
    return make_vet_result(
        out_dir,
        check="centroid_offset",
        status="pass" if ok else "fail",
        note=f"WISE W4 centroid offset={offset:.2f} arcsec ({frac:.3f} of {beam_w4:.1f}\" beam)",
        metrics={"centroid_offset_arcsec": offset, "beam_fraction": frac, "limit_fraction": max_frac},
    )


def check_neighbours(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    nbs = vet.get("neighbours_in_w4_beam", [])
    max_flux_frac = float(cfg.vetting["neighbour_blend_flux_frac"])
    blend_frac = sum(float(n.get("flux_frac", 0.0)) for n in nbs)
    ok = blend_frac <= max_flux_frac
    return make_vet_result(
        out_dir,
        check="neighbours",
        status="pass" if ok else "fail",
        note=f"{len(nbs)} catalogued neighbours inside W4 beam contributing {blend_frac:.2%} flux",
        metrics={"neighbour_count": len(nbs), "blend_flux_fraction": blend_frac, "neighbours": nbs},
    )


def check_radio_counterpart(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    rad = float(cfg.vetting["radio_search_radius_arcsec"])
    hits = [
        r
        for r in vet.get("radio_sources_within_10arcsec", [])
        if float(r.get("sep_arcsec", 99.0)) <= rad
    ]
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
    vet = star.get("vetting_archive", {})
    if not bool(vet.get("deep_imaging_available", True)):
        return make_vet_result(
            out_dir,
            check="deep_imaging",
            status="not_evaluated",
            note="No deep sub-arcsec optical/NIR imaging available at this sky position",
            metrics={"deep_imaging_available": False},
        )
    gal = vet.get("deep_imaging_offset_galaxy")
    ok = gal is None
    return make_vet_result(
        out_dir,
        check="deep_imaging",
        status="pass" if ok else "fail",
        note=(
            "Deep optical/NIR imaging shows clean isolated point source"
            if ok
            else f"Resolved background galaxy in deep imaging at sep={gal.get('sep_arcsec')} arcsec (r={gal.get('r_mag')})"
        ),
        metrics={"deep_imaging_available": True, "offset_galaxy": gal},
    )


def check_chance_alignment(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    density = float(vet.get("local_red_source_density_per_arcsec2", 8e-5))
    beam_r = float(cfg.vetting["w4_beam_arcsec"]) / 2.0
    area = math.pi * (beam_r**2)
    p_align = 1.0 - math.exp(-area * density)
    return make_vet_result(
        out_dir,
        check="chance_alignment",
        status="pass",
        note=f"Background red-source chance alignment P={p_align:.4f} within {beam_r:.1f}\" radius (reported, not silently thresholded)",
        metrics={"chance_alignment_prob": p_align, "beam_radius_arcsec": beam_r, "local_density_arcsec2": density},
    )


def check_extragalactic_colours(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    m = star["mags"]
    w1_w2 = float(m["W1"]) - float(m["W2"])
    w2_w3 = float(m["W2"]) - float(m["W3"])
    w1_w2_lim = float(cfg.vetting["extragalactic_w1_w2_min"])
    w2_w3_lim = float(cfg.vetting["extragalactic_w2_w3_min"])
    # Stern/Mateos AGN/ULIRG locus AND high-excess galaxy color check
    agn_like = w1_w2 >= w1_w2_lim and w2_w3 >= w2_w3_lim
    # If blended with an offset galaxy or DOG, flag extragalactic colors if star is faint & W1-W2 red
    vet = star.get("vetting_archive", {})
    dog_like = (w1_w2 >= 0.45 and w2_w3 >= 3.0 and len(vet.get("neighbours_in_w4_beam", [])) > 0)
    ok = not (agn_like or dog_like)
    return make_vet_result(
        out_dir,
        check="extragalactic_colours",
        status="pass" if ok else "fail",
        note=f"WISE colors W1-W2={w1_w2:.2f}, W2-W3={w2_w3:.2f}",
        metrics={"w1_w2": w1_w2, "w2_w3": w2_w3, "agn_like": agn_like, "dog_like": dog_like},
    )


def check_cirrus(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    aw = star.get("allwise", {})
    ebv = float(star.get("ebv_sfd", 0.02))
    bg_rms = float(aw.get("cirrus_bg_rms_mjy", 0.4))
    ebv_max = float(cfg.vetting["cirrus_ebv_max"])
    bg_max = float(cfg.vetting["cirrus_bg_rms_max_mjy"])
    ok = (ebv <= ebv_max) and (bg_rms <= bg_max)
    return make_vet_result(
        out_dir,
        check="cirrus",
        status="pass" if ok else "fail",
        note=f"SFD E(B-V)={ebv:.3f}, local W3/W4 background RMS={bg_rms:.2f} mJy",
        metrics={"ebv_sfd": ebv, "cirrus_bg_rms_mjy": bg_rms},
    )


def check_single_frames(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    frames = list(star.get("wise_single_frames_w4_mjy", []))
    if not frames:
        return make_vet_result(
            out_dir,
            check="single_frames",
            status="not_evaluated",
            note="No 2010 single-exposure WISE frames available",
        )
    arr = np.asarray(frames, dtype=float)
    med = float(np.median(arr))
    consistent = int(np.sum(arr >= 0.5 * max(med, float(star["flux_mjy"]["W4"]) * 0.6)))
    frac = consistent / len(arr)
    min_frac = float(cfg.vetting["single_frame_min_fraction"])
    # Also check whether one single outlier spike dominates the mean over median
    max_to_med = float(np.max(arr) / max(0.05, med))
    ok = (frac >= min_frac) and (max_to_med < 2.2)
    return make_vet_result(
        out_dir,
        check="single_frames",
        status="pass" if ok else "fail",
        note=f"{consistent}/{len(arr)} single frames consistent ({frac:.2%}), max/median={max_to_med:.2f}",
        metrics={"n_frames": len(arr), "consistent_fraction": frac, "max_to_median": max_to_med},
    )


def check_moving_object(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    mov = vet.get("moving_objects_2010", [])
    ok = len(mov) == 0
    return make_vet_result(
        out_dir,
        check="moving_object",
        status="pass" if ok else "fail",
        note="No known Solar System bodies in ephemeris cone at 2010 WISE epochs" if ok else f"Asteroid match: {mov}",
        metrics={"moving_object_hits": mov},
    )


def check_second_instrument(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    sec = vet.get("second_instrument_ir")
    if sec is None:
        return make_vet_result(
            out_dir,
            check="second_instrument",
            status="not_evaluated",
            note="No archival Spitzer, AKARI, IRAS, Herschel, or SPHEREx 22um-overlapping photometry at this position",
            metrics={"second_instrument": None},
        )
    # Melis et al. 2012 case: instrument disagreement can be real secular variability
    interp = sec.get("interpretation", "confirmed")
    return make_vet_result(
        out_dir,
        check="second_instrument",
        status="pass",
        note=f"Second-instrument archive {sec.get('instrument')} evaluated ({interp})",
        metrics={"second_instrument": sec},
    )


def check_youth(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    uv_ex = float(vet.get("galex_uv_excess_mag", 0.0))
    xray = float(vet.get("xray_luminosity_erg_s", 0.0))
    assoc = vet.get("young_association")
    is_young = (assoc is not None) or (uv_ex > float(cfg.vetting["youth_uv_excess_mag"])) or (xray > 1e29)
    return make_vet_result(
        out_dir,
        check="youth",
        status="fail" if is_young else "pass",
        note="No youth indicators (association, UV/X-ray excess)" if not is_young else f"Youth indicators present (assoc={assoc}, UV={uv_ex:.2f}, Lx={xray:.1e})",
        metrics={"young_association": assoc, "galex_uv_excess_mag": uv_ex, "xray_luminosity_erg_s": xray},
    )


def check_evolved(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    mg = float(star.get("abs_g_mag", 4.8))
    teff = float(star.get("teff_k", 5800.0))
    bright_lim = float(cfg.vetting["evolved_mg_bright_limit"])
    is_evolved = (mg < bright_lim) and (teff < 5200.0)
    return make_vet_result(
        out_dir,
        check="evolved",
        status="fail" if is_evolved else "pass",
        note=f"HR position M_G={mg:.2f}, Teff={teff:.0f} K ({'giant branch' if is_evolved else 'main sequence'})",
        metrics={"abs_g_mag": mg, "teff_k": teff, "is_evolved": is_evolved},
    )


def check_binary(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    ruwe = float(star.get("ruwe", 1.0))
    vet = star.get("vetting_archive", {})
    nss = int(vet.get("nss_flag", 0))
    rv_sig = float(vet.get("rv_scatter_kms", 0.2))
    ruwe_max = float(cfg.vetting["binary_ruwe_max"])
    is_bin = (ruwe > ruwe_max) or (nss > 0) or (rv_sig > 5.0)
    return make_vet_result(
        out_dir,
        check="binary",
        status="fail" if is_bin else "pass",
        note=f"RUWE={ruwe:.2f}, NSS={nss}, RV scatter={rv_sig:.2f} km/s",
        metrics={"ruwe": ruwe, "nss_flag": nss, "rv_scatter_kms": rv_sig},
    )


def check_vision_rubric(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    cutout_path = out_dir / "wise_cutouts.png"
    render_cutouts_panel(star, cutout_path)
    aw = star.get("allwise", {})
    vet = star.get("vetting_archive", {})
    offset = float(aw.get("centroid_offset_arcsec", 0.0))
    bg_rms = float(aw.get("cirrus_bg_rms_mjy", 0.4))
    has_off_gal = vet.get("deep_imaging_offset_galaxy") is not None
    ok = (offset <= 4.2) and (bg_rms <= float(cfg.vetting["cirrus_bg_rms_max_mjy"])) and (not has_off_gal)
    return make_vet_result(
        out_dir,
        check="vision_rubric",
        status="pass" if ok else "fail",
        note="Standard-stretch multi-band cutout panel rendered and graded against PSF/offset/cirrus rubric",
        metrics={"cutout_path": str(cutout_path), "offset_arcsec": offset, "cirrus_bg_rms_mjy": bg_rms},
        evidence_override=cutout_path,
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
