"""Code fingerprints for validation certificates (a certificate goes stale when the code it validated changes).

The fingerprint of a channel is a SHA-256 over the `oddstar` source files that determine its validated behaviour:

1. every module of the channel package (`oddstar/channels/<channel>/*.py`) and the transitive closure of the
   `oddstar` modules they import (module-level and function-level imports, found statically with `ast`);
2. the validation machinery: `oddstar.validate.harness` (which runs the four tests and imports the fixture generator
   `oddstar.io.fixtures_data`, `oddstar.io.physics`, the vetting engine and the tier assignment used by the impostor
   test) and the transitive closure of its imports, plus `oddstar.io.physics` and `oddstar.io.fixtures_data`
   explicitly.

Two rules keep one channel's fingerprint independent of the other channels: the registry `oddstar.channels`
(`build_channel`) is fingerprinted as a file but its imports (every channel) are not followed, and from the
validation roots no other `oddstar.channels.<name>` package is entered (the harness imports `ir_excess` and
`radio_narrowband` helpers that only their own branches use). Another channel's code is included only when the
channel's own package imports it (e.g. `generic_outlier` imports `underluminous.predict_ms_abs_g`).

Each file is read as UTF-8 with CRLF/CR normalised to LF and hashed; the fingerprint is the SHA-256 of the sorted
"<posix path relative to the source root>\\0<file sha256>\\n" lines, so it does not depend on the checkout location,
the platform's line endings or file-system order.
"""

from __future__ import annotations

import ast
import hashlib
from pathlib import Path

PACKAGE = "oddstar"
REGISTRY_MODULE = "oddstar.channels"
VALIDATION_ROOT_MODULES = ("oddstar.validate.harness", "oddstar.io.physics", "oddstar.io.fixtures_data")


def default_source_root() -> Path:
    """Directory that contains the imported `oddstar` package (the code that actually runs)."""
    import oddstar

    return Path(oddstar.__file__).resolve().parent.parent


def source_root_for(root: Path) -> Path:
    """Source root of a repository checkout (`<root>/src`), falling back to the imported package's location."""
    cand = root / "src"
    if (cand / PACKAGE / "__init__.py").exists():
        return cand
    return default_source_root()


def _module_file(module: str, src_root: Path) -> Path | None:
    base = src_root.joinpath(*module.split("."))
    if (base / "__init__.py").is_file():
        return base / "__init__.py"
    if base.with_suffix(".py").is_file():
        return base.with_suffix(".py")
    return None


def _module_name(path: Path, src_root: Path) -> str:
    rel = path.relative_to(src_root).with_suffix("")
    parts = list(rel.parts)
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def _imported_modules(path: Path, src_root: Path) -> set[str]:
    """`oddstar` modules imported anywhere in `path` (absolute and relative imports, including function bodies)."""
    mod = _module_name(path, src_root)
    is_pkg = path.name == "__init__.py"
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    found: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                if alias.name == PACKAGE or alias.name.startswith(PACKAGE + "."):
                    found.add(alias.name)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                pkg_parts = mod.split(".") if is_pkg else mod.split(".")[:-1]
                pkg_parts = pkg_parts[: len(pkg_parts) - (node.level - 1)]
                base = ".".join([*pkg_parts, node.module] if node.module else pkg_parts)
            else:
                base = node.module or ""
            if not (base == PACKAGE or base.startswith(PACKAGE + ".")):
                continue
            found.add(base)
            for alias in node.names:
                sub = f"{base}.{alias.name}"
                if _module_file(sub, src_root) is not None:
                    found.add(sub)
    return found


def _closure(roots: list[str], src_root: Path, *, allowed_channel: str | None) -> set[str]:
    """Modules reachable from `roots`; channel packages other than `allowed_channel` are not entered."""
    seen: set[str] = set()
    stack = list(roots)
    while stack:
        mod = stack.pop()
        if mod in seen:
            continue
        if mod.startswith(REGISTRY_MODULE + ".") and allowed_channel is not None:
            if not (mod == allowed_channel or mod.startswith(allowed_channel + ".")):
                continue
        path = _module_file(mod, src_root)
        if path is None:
            continue
        seen.add(mod)
        if mod == REGISTRY_MODULE:
            continue  # name -> class dispatch only; do not pull in every channel
        stack.extend(sorted(_imported_modules(path, src_root)))
    return seen


def channel_source_files(channel: str, src_root: Path | None = None) -> list[str]:
    """Sorted source files (posix paths relative to `src_root`) that determine `channel`'s validated behaviour."""
    src = src_root or default_source_root()
    pkg_dir = src / PACKAGE / "channels" / channel
    if not (pkg_dir / "__init__.py").is_file():
        raise FileNotFoundError(f"Channel package not found: {pkg_dir}")
    chan_mods = sorted(_module_name(p, src) for p in pkg_dir.rglob("*.py"))
    # 1. the channel package and everything it imports (other channels only when it imports them itself)
    mods = _closure(chan_mods, src, allowed_channel=None)
    # 2. the validation harness, fixture generator and physics, without entering other channel packages
    mods |= _closure(list(VALIDATION_ROOT_MODULES), src, allowed_channel=f"{REGISTRY_MODULE}.{channel}")
    files = {_module_file(m, src) for m in mods}
    return sorted(p.relative_to(src).as_posix() for p in files if p is not None)


def file_sha256_normalized(path: Path) -> str:
    """SHA-256 of a text file with CRLF and CR line endings normalised to LF."""
    text = path.read_bytes().decode("utf-8")
    norm = text.replace("\r\n", "\n").replace("\r", "\n")
    return hashlib.sha256(norm.encode("utf-8")).hexdigest()


def compute_code_fingerprint(channel: str, src_root: Path | None = None) -> tuple[str, dict[str, str]]:
    """(fingerprint, {relative path: file sha256}) for `channel`."""
    src = src_root or default_source_root()
    per_file = {rel: file_sha256_normalized(src / rel) for rel in channel_source_files(channel, src)}
    digest = hashlib.sha256()
    for rel in sorted(per_file):
        digest.update(f"{rel}\0{per_file[rel]}\n".encode("utf-8"))
    return digest.hexdigest(), per_file


def diff_file_hashes(recorded: dict[str, str] | None, current: dict[str, str]) -> list[str]:
    """Human-readable list of files added, removed or changed since `recorded`."""
    if not recorded:
        return []
    out: list[str] = []
    for rel in sorted(set(recorded) | set(current)):
        if rel not in current:
            out.append(f"removed {rel}")
        elif rel not in recorded:
            out.append(f"added {rel}")
        elif recorded[rel] != current[rel]:
            out.append(f"changed {rel}")
    return out
