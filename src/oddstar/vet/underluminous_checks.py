"""Under-luminosity vetting checks (Section 6).

Input rule (tri-state honesty): a check may only `pass` on inputs that were actually measured.
* Live stars (`"live_enriched" in star`): an input that comes from a vetting query counts only when
  `star["vetting_archive"]["query_status"][<name>] == "ok"` (`galex`, `gaia_meta` for NSS/RV fields, `irsa_dust`);
  the 3D map counts only when `vetting_archive.dust3d.status == "ok"`.
* Fixtures and hand-built dicts: an explicit value counts; an absent key or None is missing.
A check whose available inputs all pass but which lacks one of its inputs returns `not_evaluated`, naming the
missing input; a failing available input still fails.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

from oddstar.config import PipelineConfig
from oddstar.io.physics import line_of_sight_ebv, map3d_ebv, star_applied_ebv
from oddstar.schemas import VettingResult
from oddstar.vet.common import make_vet_result


def _is_live(star: dict[str, Any]) -> bool:
    return "live_enriched" in star


def _query_ok(star: dict[str, Any], name: str) -> bool:
    qs = (star.get("vetting_archive") or {}).get("query_status") or {}
    return qs.get(name) == "ok"


def _vetting_input(star: dict[str, Any], key: str, query: str) -> Any:
    """Value of vetting_archive[key] if it counts as measured (see module docstring), else None."""
    va = star.get("vetting_archive") or {}
    if _is_live(star) and not _query_ok(star, query):
        return None
    return va.get(key)


def _star_value(star: dict[str, Any], key: str) -> float | None:
    val = star.get(key)
    if val is None:
        return None
    try:
        f = float(val)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def check_astrometry(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    ruwe = _star_value(star, "ruwe")
    plx_snr = _star_value(star, "parallax_snr")
    nss_raw = _vetting_input(star, "nss_flag", "gaia_meta")
    nss = None if nss_raw is None else int(nss_raw)
    metrics = {"ruwe": ruwe, "parallax_snr": plx_snr, "nss_flag": nss}
    bad = []
    if ruwe is not None and ruwe > float(cfg.sample.max_ruwe):
        bad.append(f"RUWE={ruwe:.2f}")
    if plx_snr is not None and plx_snr < float(cfg.sample.min_parallax_snr):
        bad.append(f"parallax_SNR={plx_snr:.1f}")
    if nss is not None and nss != 0:
        bad.append(f"NSS={nss}")
    missing = [n for n, v in (("ruwe", ruwe), ("parallax_snr", plx_snr), ("nss_flag (gaia_meta query)", nss)) if v is None]
    desc = f"RUWE={ruwe}, parallax_SNR={plx_snr}, NSS={nss}"
    if bad:
        return make_vet_result(out_dir, check="astrometry", status="fail", note=f"{desc}: {', '.join(bad)}", metrics=metrics)
    if missing:
        return make_vet_result(out_dir, check="astrometry", status="not_evaluated", note=f"Missing input: {', '.join(missing)} ({desc})", metrics=metrics)
    return make_vet_result(
        out_dir,
        check="astrometry",
        status="pass",
        note=f"RUWE={ruwe:.2f}, parallax_SNR={plx_snr:.1f}, NSS={nss}",
        metrics=metrics,
    )


def check_extinction(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    """Compare the reddening applied to the star with independent estimates of the dust towards it.

    Estimates: `ebv_star` (applied), the 3D-map value at the star's distance (`vetting_archive.dust3d`), and the
    total line-of-sight column (`vetting_archive.dust` or the top-level `ebv_sfd` map value; SF11 scale). Also
    compares two line-of-sight maps (`ebv_sfd` vs `ebv_planck`) when both exist.
    * fail: heavy dust (any estimate > `extinction_heavy_dust_ebv`: the correction dominates any luminosity
      deficit), applied value disagrees with the 3D value or exceeds the total column, or the two 2D maps disagree;
      tolerance `extinction_disagreement_mag` (in E(B-V)) or 2 sigma of the 3D value, whichever is larger.
    * not_evaluated: fewer than two estimates, or the total column greatly exceeds the applied value with no 3D map
      saying where the dust lies (intrinsic colour and luminosity uncertain).
    * pass: otherwise.
    """
    v = cfg.vetting
    tol = float(v["extinction_disagreement_mag"])
    heavy_thr = float(v["extinction_heavy_dust_ebv"])
    los_excess_thr = float(v["extinction_los_excess_summary_ebv"])
    factor = float(v["sfd_to_sf11_factor"])

    applied, applied_src = star_applied_ebv(star)
    # Line-of-sight from vetting enrichment counts for live stars only when the irsa_dust query succeeded
    star_los = star
    if _is_live(star) and not _query_ok(star, "irsa_dust"):
        va = dict(star.get("vetting_archive") or {})
        va.pop("dust", None)
        star_los = {**star, "vetting_archive": va}
    los = line_of_sight_ebv(star_los, sfd_to_sf11=factor)
    d3 = map3d_ebv(star)
    sfd_top = _star_value(star, "ebv_sfd")
    planck_top = _star_value(star, "ebv_planck")

    metrics: dict[str, Any] = {
        "ebv_star": applied,
        "ebv_star_source": applied_src,
        "ebv_los": None if los is None else los["ebv_los"],
        "ebv_los_source": None if los is None else los["source"],
        "ebv_3d": None if d3 is None else d3["ebv"],
        "ebv_3d_err": None if d3 is None else d3["ebv_err"],
        "ebv_sfd": sfd_top,
        "ebv_planck": planck_top,
        "tolerance": tol,
        "heavy_dust_threshold": heavy_thr,
    }
    estimates = [x for x in (applied, None if los is None else los["ebv_los"], None if d3 is None else d3["ebv"]) if x is not None]
    desc = (
        f"applied E(B-V)={'n/a' if applied is None else f'{applied:.3f}'} [{applied_src}], "
        f"line-of-sight={'n/a' if los is None else f'{los['ebv_los']:.3f}'}, "
        f"3D at distance={'n/a' if d3 is None else f'{d3['ebv']:.3f}+/-{d3['ebv_err']:.3f}'}"
    )
    if applied is None or len(estimates) < 2:
        return make_vet_result(
            out_dir,
            check="extinction",
            status="not_evaluated",
            note=f"Fewer than two reddening estimates ({desc}); need ebv_star plus a 3D-map or line-of-sight value",
            metrics=metrics,
        )

    problems: list[str] = []
    heavy = max(estimates) > heavy_thr
    metrics["heavy_dust"] = heavy
    if heavy:
        problems.append(f"heavy dust (max estimate {max(estimates):.3f} > {heavy_thr:.2f})")
    if d3 is not None:
        tol3 = max(tol, 2.0 * d3["ebv_err"])
        diff3 = abs(applied - d3["ebv"])
        metrics["disagreement_3d"] = diff3
        if diff3 > tol3:
            problems.append(f"applied vs 3D differ by {diff3:.3f} (> {tol3:.3f})")
    if los is not None and applied - los["ebv_los"] > tol:
        problems.append(f"applied exceeds the whole line-of-sight column by {applied - los['ebv_los']:.3f}")
    if sfd_top is not None and planck_top is not None and abs(sfd_top - planck_top) > tol:
        problems.append(f"SFD vs Planck line-of-sight maps differ by {abs(sfd_top - planck_top):.3f}")
    if problems:
        return make_vet_result(out_dir, check="extinction", status="fail", note=f"{desc}: {'; '.join(problems)}", metrics=metrics)

    if los is not None and d3 is None and (los["ebv_los"] - applied) > los_excess_thr:
        metrics["los_minus_applied"] = los["ebv_los"] - applied
        return make_vet_result(
            out_dir,
            check="extinction",
            status="not_evaluated",
            note=(
                f"{desc}: the line-of-sight column exceeds the applied reddening by {los['ebv_los'] - applied:.3f} and no 3D "
                "map says where that dust lies, so the star's intrinsic colour and luminosity are uncertain (see the "
                "two-case extinction summary)"
            ),
            metrics=metrics,
        )
    return make_vet_result(out_dir, check="extinction", status="pass", note=f"{desc}: consistent", metrics=metrics)


def check_metallicity(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    feh = _star_value(star, "feh")
    stype = str(star.get("simbad_type") or "")
    subdwarf_type = "Subdwarf" in stype
    is_subdwarf = subdwarf_type or (feh is not None and feh < -1.0)
    metrics = {"feh": feh, "is_subdwarf": is_subdwarf}
    if is_subdwarf:
        return make_vet_result(
            out_dir,
            check="metallicity",
            status="fail",
            note=f"[Fe/H]={feh} dex, SIMBAD type '{stype}' (metal-poor subdwarf locus)",
            metrics=metrics,
        )
    if feh is None:
        return make_vet_result(out_dir, check="metallicity", status="not_evaluated", note="Missing input: [Fe/H] (feh)", metrics=metrics)
    return make_vet_result(
        out_dir,
        check="metallicity",
        status="pass",
        note=f"[Fe/H]={feh:.2f} dex (disk metallicity)",
        metrics=metrics,
    )


def check_parameter_agreement(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    t_xp = _star_value(star, "teff_k")
    t_spec = _star_value(star, "teff_spec_k")
    metrics = {"teff_xp_k": t_xp, "teff_spec_k": t_spec}
    if _is_live(star) and not star.get("teff_spec_source"):
        # Live rows copy Gaia GSP-Phot Teff into teff_spec_k: not an independent spectroscopic determination
        return make_vet_result(
            out_dir,
            check="parameter_agreement",
            status="not_evaluated",
            note="Missing input: independent spectroscopic Teff (live teff_spec_k is a copy of Gaia GSP-Phot Teff)",
            metrics=metrics,
        )
    if t_xp is None or t_spec is None:
        missing = [n for n, val in (("teff_k", t_xp), ("teff_spec_k", t_spec)) if val is None]
        return make_vet_result(out_dir, check="parameter_agreement", status="not_evaluated", note=f"Missing input: {', '.join(missing)}", metrics=metrics)
    diff = abs(t_xp - t_spec)
    metrics["diff_k"] = diff
    ok = diff <= 300.0
    return make_vet_result(
        out_dir,
        check="parameter_agreement",
        status="pass" if ok else "fail",
        note=f"Teff_XP={t_xp:.0f} K vs Teff_spec={t_spec:.0f} K (diff={diff:.0f} K)",
        metrics=metrics,
    )


def check_companion(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> VettingResult:
    uv_raw = _vetting_input(star, "galex_uv_excess_mag", "galex")
    uv = None if uv_raw is None else float(uv_raw)
    ruwe = _star_value(star, "ruwe")
    metrics = {"galex_uv_excess_mag": uv, "ruwe": ruwe}
    bad = []
    if uv is not None and uv > float(cfg.vetting["youth_uv_excess_mag"]):
        bad.append(f"GALEX UV excess={uv:.2f} mag")
    if ruwe is not None and ruwe > float(cfg.sample.max_ruwe):
        bad.append(f"RUWE={ruwe:.2f}")
    if bad:
        return make_vet_result(out_dir, check="companion", status="fail", note="; ".join(bad), metrics=metrics)
    missing = [n for n, val in (("galex_uv_excess_mag (galex query)", uv), ("ruwe", ruwe)) if val is None]
    if missing:
        return make_vet_result(out_dir, check="companion", status="not_evaluated", note=f"Missing input: {', '.join(missing)}", metrics=metrics)
    return make_vet_result(
        out_dir,
        check="companion",
        status="pass",
        note=f"GALEX UV excess={uv:.2f} mag, RUWE={ruwe:.2f}",
        metrics=metrics,
    )


def run_underluminous_vetting(star: dict[str, Any], cfg: PipelineConfig, out_dir: Path) -> list[VettingResult]:
    checks = [check_astrometry, check_extinction, check_metallicity, check_parameter_agreement, check_companion]
    return [fn(star, cfg, out_dir) for fn in checks]
