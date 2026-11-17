"""CI governance check: committed validation certificates must match the working tree (config hash + code fingerprint).

Wire into `oddstar.ci_checks.run_all_ci_checks` as
`"certificates": check_certificates_current(root),`.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path

from pydantic import ValidationError

from oddstar.config import compute_config_hash, load_config
from oddstar.schemas import ValidationCertificate
from oddstar.validate.fingerprint import compute_code_fingerprint, diff_file_hashes, source_root_for

CERT_DIR = Path("validation") / "certificates"
# JSON files in the certificate directory that are reports, not certificates
NON_CERTIFICATE_FILES = {"dr3_vs_dr4_comparison.json"}


def _certificate_files(root: Path) -> list[Path]:
    """Certificate JSON files tracked by git under validation/certificates (all on disk when git is unavailable)."""
    cert_dir = root / CERT_DIR
    try:
        res = subprocess.run(
            ["git", "ls-files", "-z", "--", CERT_DIR.as_posix()],
            cwd=str(root),
            capture_output=True,
            text=True,
            check=False,
        )
        tracked = [root / p for p in res.stdout.split("\0") if p] if res.returncode == 0 else None
    except OSError:
        tracked = None
    files = tracked if tracked is not None else sorted(cert_dir.glob("*.json")) if cert_dir.exists() else []
    return sorted(p for p in files if p.suffix == ".json" and p.name not in NON_CERTIFICATE_FILES and p.exists())


def check_certificates_current(root: Path) -> list[str]:
    """Report every committed certificate whose config hash or code fingerprint differs from the working tree.

    For each certificate the working-tree config (`configs/default.yaml` with the certificate's Gaia release) is
    hashed and the channel's code fingerprint is computed from `<root>/src`. A certificate without a fingerprint is
    reported too: it cannot show which code it validated. Superseded certificates (older config hashes) must be
    removed from the repository, since they no longer certify anything.
    """
    violations: list[str] = []
    files = _certificate_files(root)
    if not files:
        return violations
    src_root = source_root_for(root)
    cfg_path = root / "configs" / "default.yaml"
    hashes: dict[str, str] = {}
    fingerprints: dict[str, tuple[str, dict[str, str]]] = {}
    for path in files:
        rel = path.relative_to(root).as_posix()
        try:
            cert = ValidationCertificate.model_validate_json(path.read_text(encoding="utf-8"))
        except (ValidationError, json.JSONDecodeError, UnicodeDecodeError) as exc:
            violations.append(f"{rel}: not a readable ValidationCertificate ({type(exc).__name__})")
            continue
        if cert.gaia_release not in hashes:
            cfg, _, _ = load_config(cfg_path, root=root, gaia_release_override=cert.gaia_release)
            hashes[cert.gaia_release] = compute_config_hash(cfg.model_dump())[:12]
        cur_hash = hashes[cert.gaia_release]
        if cert.config_hash != cur_hash:
            violations.append(f"{rel}: config_hash {cert.config_hash} != working tree {cur_hash} ({cert.gaia_release})")
        if cert.channel not in fingerprints:
            try:
                fingerprints[cert.channel] = compute_code_fingerprint(cert.channel, src_root)
            except FileNotFoundError as exc:
                violations.append(f"{rel}: channel '{cert.channel}' has no source package ({exc})")
                continue
        cur_fp, cur_files = fingerprints[cert.channel]
        if cert.code_fingerprint is None:
            violations.append(f"{rel}: no code_fingerprint (issued before code fingerprinting); re-validate")
        elif cert.code_fingerprint != cur_fp:
            changed = diff_file_hashes(cert.details.get("code_fingerprint_files"), cur_files)
            detail = f" [{', '.join(changed)}]" if changed else ""
            violations.append(
                f"{rel}: code_fingerprint {cert.code_fingerprint[:12]} != working tree {cur_fp[:12]}{detail}"
            )
    return violations
