"""Under-luminosity vetting checks (Section 6)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result


def check_astrometry(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    ruwe = float(star.get("ruwe", 1.0))
    plx_snr = float(star.get("parallax_snr", 50.0))
    nss = int(star.get("vetting_archive", {}).get("nss_flag", 0))
    ok = (ruwe <= float(cfg.sample.max_ruwe)) and (plx_snr >= float(cfg.sample.min_parallax_snr)) and (nss == 0)
    return make_vet_result(
        out_dir,
        check="astrometry",
        status="pass" if ok else "fail",
        note=f"RUWE={ruwe:.2f}, parallax_SNR={plx_snr:.1f}, NSS={nss}",
        metrics={"ruwe": ruwe, "parallax_snr": plx_snr, "nss_flag": nss},
    )


def check_extinction(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    sfd = float(star.get("ebv_sfd", 0.02))
    planck = float(star.get("ebv_planck", sfd))
    diff = abs(sfd - planck)
    max_diff = float(cfg.vetting["extinction_disagreement_mag"])
    # Heavily extinguished stars behind dust (E(B-V) > 0.3) are explained once extinction correction is applied
    high_dust = max(sfd, planck) > 0.30
    ok = (diff <= max_diff) and (not high_dust)
    return make_vet_result(
        out_dir,
        check="extinction",
        status="pass" if ok else "fail",
        note=f"E(B-V)_SFD={sfd:.3f}, E(B-V)_Planck={planck:.3f} (diff={diff:.3f}, heavy_dust={high_dust})",
        metrics={"ebv_sfd": sfd, "ebv_planck": planck, "disagreement": diff, "heavy_dust": high_dust},
    )


def check_metallicity(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    feh = float(star.get("feh", 0.0))
    stype = str(star.get("simbad_type", ""))
    is_subdwarf = (feh < -1.0) or ("Subdwarf" in stype)
    return make_vet_result(
        out_dir,
        check="metallicity",
        status="fail" if is_subdwarf else "pass",
        note=f"[Fe/H]={feh:.2f} dex ({'metal-poor subdwarf locus' if is_subdwarf else 'disk metallicity'})",
        metrics={"feh": feh, "is_subdwarf": is_subdwarf},
    )


def check_parameter_agreement(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    t_xp = float(star.get("teff_k", 5800.0))
    t_spec = float(star.get("teff_spec_k", t_xp))
    diff = abs(t_xp - t_spec)
    ok = diff <= 300.0
    return make_vet_result(
        out_dir,
        check="parameter_agreement",
        status="pass" if ok else "fail",
        note=f"Teff_XP={t_xp:.0f} K vs Teff_spec={t_spec:.0f} K (diff={diff:.0f} K)",
        metrics={"teff_xp_k": t_xp, "teff_spec_k": t_spec, "diff_k": diff},
    )


def check_companion(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    vet = star.get("vetting_archive", {})
    uv = float(vet.get("galex_uv_excess_mag", 0.0))
    ruwe = float(star.get("ruwe", 1.0))
    bad = (uv > float(cfg.vetting["youth_uv_excess_mag"])) or (ruwe > float(cfg.sample.max_ruwe))
    return make_vet_result(
        out_dir,
        check="companion",
        status="fail" if bad else "pass",
        note=f"GALEX UV excess={uv:.2f} mag, RUWE={ruwe:.2f}",
        metrics={"galex_uv_excess_mag": uv, "ruwe": ruwe},
    )


def run_underluminous_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [check_astrometry, check_extinction, check_metallicity, check_parameter_agreement, check_companion]
    return [fn(star, cfg, out_dir) for fn in checks]
