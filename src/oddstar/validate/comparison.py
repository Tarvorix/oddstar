"""Regenerate `validation/certificates/dr3_vs_dr4_comparison.{json,md}` from the current certificates.

Run `python -m oddstar.validate.comparison` after `oddstar validate <channel>` has been run for every channel on both
releases. Only certificates that `search` would accept (current config hash, passed, current code fingerprint) are
used; a missing or stale certificate is an error, so the comparison can never describe certificates that are not in force.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from oddstar.channels import ALL_CHANNEL_NAMES
from oddstar.config import load_config, repo_root
from oddstar.errors import ArtifactNotFoundError
from oddstar.schemas import ValidationCertificate
from oddstar.validate.harness import certificate_status

OPEN_CASE_SOURCE_ID = "4582145718700356096"


def build_comparison(root: Path) -> dict[str, Any]:
    from oddstar.cli import resolve_star_record

    certs: dict[str, dict[str, ValidationCertificate]] = {}
    snr: dict[str, float] = {}
    hashes: dict[str, str] = {}
    for rel in ("dr3", "dr4"):
        cfg, short_hash, _ = load_config(root=root, gaia_release_override=rel)
        hashes[rel] = short_hash
        snr[rel] = float(resolve_star_record(OPEN_CASE_SOURCE_ID, cfg)["parallax_snr"])
        certs[rel] = {}
        for ch in ALL_CHANNEL_NAMES:
            row = certificate_status(ch, cfg, root)
            if not row["exists"] or row["stale"]:
                raise ArtifactNotFoundError(f"No current certificate for {ch} ({rel}): {'; '.join(row['reasons'])}")
            certs[rel][ch] = ValidationCertificate.model_validate_json((root / row["certificate"]).read_text(encoding="utf-8"))
    comparison = []
    for ch in ALL_CHANNEL_NAMES:
        c3, c4 = certs["dr3"][ch], certs["dr4"][ch]
        comparison.append(
            {
                "channel": ch,
                "dr3_passed": bool(c3.passed),
                "dr4_passed": bool(c4.passed),
                "dr3_blind_rate": c3.blind_recovery_rate,
                "dr4_blind_rate": c4.blind_recovery_rate,
                "dr3_far": c3.false_alarm_rate,
                "dr4_far": c4.false_alarm_rate,
                "dr3_inj_dp": c3.injection_completeness_design_point,
                "dr4_inj_dp": c4.injection_completeness_design_point,
            }
        )
    return {
        "channels_revalidated": len(comparison),
        "config_hash_dr3": hashes["dr3"],
        "config_hash_dr4": hashes["dr4"],
        "open_case_parallax_snr_dr3": snr["dr3"],
        "open_case_parallax_snr_dr4": snr["dr4"],
        "snr_improvement_factor": round(snr["dr4"] / snr["dr3"], 2),
        "comparison": comparison,
    }


def write_comparison(root: Path) -> tuple[Path, Path]:
    comp = build_comparison(root)
    out_dir = root / "validation" / "certificates"
    jp = out_dir / "dr3_vs_dr4_comparison.json"
    mp = out_dir / "dr3_vs_dr4_comparison.md"
    jp.write_text(json.dumps(comp, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Gaia DR3 vs DR4 Re-Validation Comparison Report",
        "",
        f"- Config hash: DR3 `{comp['config_hash_dr3']}`, DR4 `{comp['config_hash_dr4']}`",
        "- Validation runs on the offline fixture set (synthetic stars and literature analogues), not on live archive data.",
        f"- Open case {OPEN_CASE_SOURCE_ID} parallax SNR: DR3={comp['open_case_parallax_snr_dr3']:.1f} -> "
        f"DR4={comp['open_case_parallax_snr_dr4']:.1f} ({comp['snr_improvement_factor']:.2f}x)",
        "",
        "| Channel | DR3 Passed | DR4 Passed | DR3 Blind | DR4 Blind | DR3 Inj DP | DR4 Inj DP | DR3 FAR | DR4 FAR |",
        "| --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for r in comp["comparison"]:
        lines.append(
            f"| `{r['channel']}` | {r['dr3_passed']} | {r['dr4_passed']} | {r['dr3_blind_rate']:.2f} | {r['dr4_blind_rate']:.2f} "
            f"| {r['dr3_inj_dp']:.2f} | {r['dr4_inj_dp']:.2f} | {r['dr3_far']:.4f} | {r['dr4_far']:.4f} |"
        )
    mp.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return jp, mp


if __name__ == "__main__":
    for p in write_comparison(repo_root()):
        print(p)
