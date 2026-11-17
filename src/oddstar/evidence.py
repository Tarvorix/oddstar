"""Committed evidence store for registered candidates (AGENTS.md: never mark a check `pass` without an evidence file).

Vetting evidence written by `make_vet_result` lives under git-ignored `runs/` or `dossiers/`. When `fuse` or
`inspect --register` registers a candidate, `export_candidate_evidence` copies a compact, size-guarded version into
the tracked directory `registry/evidence/<source_id>/<run_id>/`:

- `vet_<check>.json` per vetting check (check, status, note, reduced metrics, and the raw-cache manifest IDs, SHA-256
  values and query URLs found in the check's evidence and metrics),
- `measurement_<channel>.json` per channel measurement (values, uncertainties, caveats, reduced extras, score),
- small rendered evidence PNGs (above `evidence.max_image_bytes` only their SHA-256 and size are recorded),
- `index.json`: every file with its SHA-256, config hash, git commit, pipeline version, UTC time, tier, flagged
  channels and the exported vetting rows.

"Never commit raw data" still holds: `reduce_large_values` replaces every long array (time series, pixel data) and
every long string with a summary plus the SHA-256 and length of the full value. It applies to all exported JSON, not to
particular checks. An evidence directory referenced by a ledger entry is never overwritten.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

from oddstar.registry import REGISTRATION_EVENTS, clean_source_id, read_ledger
from oddstar.schemas import (
    ChannelScore,
    EvidenceFile,
    EvidenceIndex,
    OmittedEvidenceFile,
    VettingResult,
)

INDEX_NAME = "index.json"
IMAGE_SUFFIXES = {".png"}
_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_SAFE = re.compile(r"[^A-Za-z0-9._+-]+")


# ------------------------------------------------------------------------------------------------ hashing


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _jsonable(obj: Any) -> Any:
    """Plain JSON types: numpy scalars/arrays to Python, non-finite floats to None, tuples/sets to lists."""
    if hasattr(obj, "tolist") and not isinstance(obj, (str, bytes)):
        obj = obj.tolist()
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple, set)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, (str, int, bool)) or obj is None:
        return obj
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def canonical_json(obj: Any) -> bytes:
    return json.dumps(_jsonable(obj), sort_keys=True, separators=(",", ":")).encode("utf-8")


def _dump(obj: Any) -> bytes:
    return (json.dumps(obj, indent=2, sort_keys=True) + "\n").encode("utf-8")


# ------------------------------------------------------------------------------------------- size guard


def _count_scalars(value: Any, stop_after: int) -> int:
    """Number of scalar leaves in a nested list, counting at most `stop_after + 1`."""
    n = 0
    stack = [value]
    while stack:
        cur = stack.pop()
        if isinstance(cur, list):
            stack.extend(cur)
        elif isinstance(cur, dict):
            stack.extend(cur.values())
        else:
            n += 1
        if n > stop_after:
            break
    return n


def _flatten_numbers(value: Any) -> list[float] | None:
    out: list[float] = []
    stack = [value]
    while stack:
        cur = stack.pop()
        if isinstance(cur, list):
            stack.extend(cur)
        elif isinstance(cur, bool):
            return None
        elif isinstance(cur, (int, float)):
            out.append(float(cur))
        elif cur is None:
            continue
        else:
            return None
    return out


def _shape(value: list[Any]) -> list[int]:
    shape = [len(value)]
    cur: Any = value
    while isinstance(cur, list) and cur and isinstance(cur[0], list):
        cur = cur[0]
        shape.append(len(cur))
    return shape


def _array_summary(value: list[Any], max_elements: int, max_chars: int) -> dict[str, Any]:
    out: dict[str, Any] = {
        "_reduced": "array",
        "length": len(value),
        "shape": _shape(value),
        "n_scalars": _count_scalars(value, 1 << 62),
        "sha256": sha256_bytes(canonical_json(value)),
    }
    nums = _flatten_numbers(value)
    if nums is not None and nums:
        finite = sorted(x for x in nums if math.isfinite(x))
        if finite:
            mid = len(finite) // 2
            median = finite[mid] if len(finite) % 2 else 0.5 * (finite[mid - 1] + finite[mid])
            mean = sum(finite) / len(finite)
            var = sum((x - mean) ** 2 for x in finite) / len(finite)
            out["summary"] = {
                "n_finite": len(finite),
                "min": finite[0],
                "max": finite[-1],
                "mean": mean,
                "median": median,
                "std": math.sqrt(var),
            }
    else:
        out["head"] = [reduce_large_values(v, max_elements=max_elements, max_chars=max_chars) for v in value[:3]]
    return out


def reduce_large_values(value: Any, *, max_elements: int, max_chars: int) -> Any:
    """General size guard: long arrays and long strings become a summary with the SHA-256 and length of the full value.

    The SHA-256 is taken over the canonical JSON of the full value (sorted keys, no spaces), so it can be recomputed
    from the run artefacts or the raw cache. Dicts are walked recursively; nothing is dropped silently.
    """
    value = _jsonable(value)
    if isinstance(value, dict):
        return {k: reduce_large_values(v, max_elements=max_elements, max_chars=max_chars) for k, v in value.items()}
    if isinstance(value, list):
        if _count_scalars(value, max_elements) > max_elements:
            return _array_summary(value, max_elements, max_chars)
        return [reduce_large_values(v, max_elements=max_elements, max_chars=max_chars) for v in value]
    if isinstance(value, str) and len(value) > max_chars:
        return {"_reduced": "string", "length": len(value), "sha256": sha256_bytes(value.encode("utf-8")), "head": value[:200]}
    return value


def guard_payload(payload: dict[str, Any], *, max_elements: int, max_chars: int, max_bytes: int) -> dict[str, Any]:
    """Apply `reduce_large_values`, tightening the limits until the JSON fits in `max_bytes`."""
    me, mc = max_elements, max_chars
    for _ in range(6):
        reduced = reduce_large_values(payload, max_elements=me, max_chars=mc)
        if len(_dump(reduced)) <= max_bytes:
            return reduced
        me, mc = max(1, me // 4), max(16, mc // 4)
    full = canonical_json(payload)
    keep = {k: reduced[k] for k in ("check", "status", "channel") if k in reduced}
    keep["_reduced"] = {"kind": "object", "bytes": len(full), "sha256": sha256_bytes(full)}
    return keep


# -------------------------------------------------------------------------------------------- inputs


def collect_inputs(*objs: Any) -> dict[str, list[str]]:
    """Raw-cache manifest IDs, SHA-256 values and query URLs mentioned anywhere in the given objects."""
    manifests: set[str] = set()
    hashes: set[str] = set()
    urls: set[str] = set()
    stack: list[tuple[str, Any]] = [("", o) for o in objs]
    while stack:
        key, cur = stack.pop()
        if hasattr(cur, "tolist") and not isinstance(cur, (str, bytes)):
            continue  # numeric arrays carry no references
        if isinstance(cur, dict):
            stack.extend((str(k), v) for k, v in cur.items())
        elif isinstance(cur, (list, tuple)):
            stack.extend((key, v) for v in cur)
        elif isinstance(cur, str):
            k = key.lower()
            if "manifest" in k and cur:
                manifests.add(cur)
            if (k == "sha256" or k.endswith("_sha256")) and _HEX64.match(cur):
                hashes.add(cur)
            for m in re.finditer(r"https?://[^\s'\"<>]+", cur):
                urls.add(m.group(0).rstrip(".,;)"))
    return {"manifest_ids": sorted(manifests), "sha256": sorted(hashes), "query_urls": sorted(urls)}


# -------------------------------------------------------------------------------------------- export


@dataclass
class EvidenceExport:
    directory: Path
    rel_dir: str
    index_sha256: str
    index: EvidenceIndex
    vetting_rows: list[dict[str, Any]] = field(default_factory=list)


def _safe(name: str) -> str:
    return _SAFE.sub("_", name).strip("._") or "x"


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.resolve().relative_to(root.resolve()))
    except ValueError:
        return str(path)


def referenced_evidence_dirs(root: Path, registry_dir: str = "registry") -> set[str]:
    return {str(r["evidence_dir"]) for r in read_ledger(root, registry_dir) if r.get("evidence_dir")}


def _target_dir(root: Path, evidence_dir: str, registry_dir: str, source_id: str, run_id: str) -> tuple[Path, str]:
    base = root / evidence_dir / _safe(source_id)
    referenced = referenced_evidence_dirs(root, registry_dir)
    name = _safe(run_id)
    candidate = base / name
    n = 2
    while _rel(candidate, root) in referenced:
        candidate = base / f"{name}.{n}"
        n += 1
    return candidate, _rel(candidate, root)


def rejection_checks(vetting: list[VettingResult]) -> list[str]:
    """Failing checks that cause Rejected (artifact/integrity failures, as in `fuse.engine.assign_tier`)."""
    from oddstar.fuse.engine import ARTIFACT_INTEGRITY_CHECKS

    return sorted({v.check for v in vetting if v.status == "fail" and v.check in ARTIFACT_INTEGRITY_CHECKS})


def export_candidate_evidence(
    *,
    root: Path,
    cfg: Any,
    source_id: str,
    run_id: str,
    tier: str,
    vetting: list[VettingResult],
    measurements: Iterable[dict[str, Any]] = (),
    scores: Iterable[ChannelScore] = (),
    images: Iterable[Path] = (),
    config_hash: str,
    config_hash_full: str = "",
    git_commit: str,
    pipeline_version: str,
) -> EvidenceExport:
    """Write the committed evidence directory for one candidate registration and return its index hash.

    Tier A/B/C/Explained export every vetting check, every channel measurement and the given images. Rejected exports
    only the failing check(s) that caused the rejection.
    """
    sid = clean_source_id(source_id)
    ev = cfg.evidence
    lim = {"max_elements": int(ev.max_array_elements), "max_chars": int(ev.max_string_chars), "max_bytes": int(ev.max_json_bytes)}
    out_dir, rel_dir = _target_dir(root, cfg.storage.evidence_dir, cfg.storage.registry_dir, sid, run_id)
    if out_dir.exists():
        shutil.rmtree(out_dir)  # unreferenced earlier export of the same run (referenced ones get a new name)
    out_dir.mkdir(parents=True, exist_ok=True)

    scores = list(scores)
    flagged = sorted({s.channel for s in scores if s.flagged})
    scope = "rejection_only" if tier == "Rejected" else "full"
    if scope == "rejection_only":
        keep = set(rejection_checks(vetting))
        vet_rows = [v for v in vetting if v.check in keep]
        measurements, scores, images = [], [], []
    else:
        vet_rows = list(vetting)

    files: list[EvidenceFile] = []
    omitted: list[OmittedEvidenceFile] = []

    def _write(name: str, data: bytes, kind: str, **extra: Any) -> str:
        p = out_dir / name
        p.write_bytes(data)
        files.append(EvidenceFile(path=name, sha256=sha256_bytes(data), bytes=len(data), kind=kind, **extra))  # type: ignore[arg-type]
        return name

    def _image(src: Path, prefix: str = "") -> dict[str, Any]:
        size = src.stat().st_size
        digest = sha256_file(src)
        name = _safe(prefix + src.name)
        if size <= int(ev.max_image_bytes):
            if not any(f.path == name for f in files):
                _write(name, src.read_bytes(), "image")
            return {"file": name, "sha256": digest, "bytes": size}
        omitted.append(OmittedEvidenceFile(name=name, sha256=digest, bytes=size,
                                           reason=f"image larger than evidence.max_image_bytes={int(ev.max_image_bytes)}"))
        return {"file": None, "sha256": digest, "bytes": size, "omitted": "over size cap"}

    rows: list[dict[str, Any]] = []
    for v in sorted(vet_rows, key=lambda x: x.check):
        src_payload: Any = None
        source_ref: dict[str, Any] | None = None
        image_refs: list[dict[str, Any]] = []
        if v.evidence:
            src = Path(v.evidence)
            if not src.is_absolute():
                src = root / src
            if src.exists() and src.is_file():
                source_ref = {"path": _rel(src, root), "sha256": sha256_file(src), "bytes": src.stat().st_size}
                if src.suffix.lower() in IMAGE_SUFFIXES:
                    image_refs.append(_image(src, prefix=f"{v.check}_"))
                else:
                    try:
                        src_payload = json.loads(src.read_text(encoding="utf-8"))
                    except (UnicodeDecodeError, json.JSONDecodeError):
                        src_payload = None
            else:
                source_ref = {"path": str(v.evidence), "missing_at_export": True}
        cutout = (v.metrics or {}).get("cutout_path")
        if cutout and Path(str(cutout)).exists() and Path(str(cutout)).suffix.lower() in IMAGE_SUFFIXES:
            if not (source_ref and source_ref.get("path") == _rel(Path(str(cutout)), root)):
                image_refs.append(_image(Path(str(cutout)), prefix=f"{v.check}_"))
        payload = {
            "check": v.check,
            "status": v.status,
            "note": v.note,
            "metrics": v.metrics,
            "inputs": collect_inputs(v.metrics, src_payload, v.note),
            "source_evidence": source_ref,
            "images": image_refs,
        }
        payload = guard_payload(payload, max_elements=lim["max_elements"], max_chars=lim["max_chars"], max_bytes=lim["max_bytes"])
        fname = _write(f"vet_{_safe(v.check)}.json", _dump(payload), "vetting", check=v.check)
        rows.append({"check": v.check, "status": v.status, "evidence": fname})

    score_by_ch = {s.channel: s for s in scores}
    seen_channels: set[str] = set()
    for m in measurements:
        ch = str(m.get("channel", "unknown"))
        seen_channels.add(ch)
        sc = score_by_ch.get(ch)
        payload = {
            "channel": ch,
            "source_id": sid,
            "values": m.get("values", {}),
            "uncertainties": m.get("uncertainties", {}),
            "epochs": m.get("epochs", []),
            "caveats": m.get("caveats", []),
            "extra": m.get("extra", {}),
            "score": sc.model_dump(exclude={"evidence"}) if sc else None,
            "inputs": collect_inputs(m.get("extra", {})),
        }
        payload = guard_payload(payload, max_elements=lim["max_elements"], max_chars=lim["max_chars"], max_bytes=lim["max_bytes"])
        _write(f"measurement_{_safe(ch)}.json", _dump(payload), "measurement", channel=ch)
    for ch, sc in sorted(score_by_ch.items()):
        if ch in seen_channels:
            continue
        payload = {"channel": ch, "source_id": sid, "values": {}, "score": sc.model_dump(exclude={"evidence"}),
                   "note": "score only: the measurement was not available at export"}
        payload = guard_payload(payload, max_elements=lim["max_elements"], max_chars=lim["max_chars"], max_bytes=lim["max_bytes"])
        _write(f"measurement_{_safe(ch)}.json", _dump(payload), "measurement", channel=ch)

    for img in images:
        p = Path(img)
        if p.exists() and p.suffix.lower() in IMAGE_SUFFIXES:
            _image(p)

    index = EvidenceIndex(
        source_id=sid,
        run_id=run_id,
        tier=tier,  # type: ignore[arg-type]
        scope=scope,  # type: ignore[arg-type]
        flagged_channels=flagged,
        config_hash=config_hash,
        config_hash_full=config_hash_full,
        git_commit=git_commit,
        pipeline_version=pipeline_version,
        timestamp_utc=datetime.now(timezone.utc).isoformat(),
        vetting=rows,
        files=sorted(files, key=lambda f: f.path),
        omitted_files=omitted,
    )
    index_bytes = _dump(index.model_dump(mode="json"))
    (out_dir / INDEX_NAME).write_bytes(index_bytes)
    return EvidenceExport(directory=out_dir, rel_dir=rel_dir, index_sha256=sha256_bytes(index_bytes), index=index, vetting_rows=rows)


# -------------------------------------------------------------------------------------------- verify


def verify_evidence_dir(directory: Path) -> list[str]:
    """Re-hash every file listed in `index.json`; report missing, altered and unlisted files."""
    problems: list[str] = []
    idx_path = directory / INDEX_NAME
    if not idx_path.exists():
        return [f"{directory}: index.json missing"]
    try:
        index = EvidenceIndex.model_validate_json(idx_path.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - any parse failure is a verification failure
        return [f"{idx_path}: unreadable index ({exc})"]
    listed = set()
    for f in index.files:
        listed.add(f.path)
        p = directory / f.path
        if not p.exists():
            problems.append(f"{p}: listed in index.json but missing")
        elif sha256_file(p) != f.sha256:
            problems.append(f"{p}: SHA-256 mismatch (index {f.sha256[:12]}, file {sha256_file(p)[:12]})")
    for p in sorted(directory.iterdir()):
        if p.is_file() and p.name != INDEX_NAME and p.name not in listed:
            problems.append(f"{p}: file not listed in index.json")
    return problems


def verify_evidence(root: Path, cfg_storage: Any | None = None, source_id: str | None = None) -> dict[str, Any]:
    """Verify every evidence directory (optionally one source) and every ledger reference to evidence.

    Returns `problems` (empty when everything matches), `verified_dirs`, `ledger_entries_with_evidence` and
    `legacy_without_evidence` (registrations written before the evidence store existed; reported, never failed).
    """
    evidence_dir = getattr(cfg_storage, "evidence_dir", "registry/evidence")
    registry_dir = getattr(cfg_storage, "registry_dir", "registry")
    sid_filter = clean_source_id(source_id) if source_id else None
    problems: list[str] = []
    verified: list[str] = []
    base = root / evidence_dir
    if base.exists():
        for idx in sorted(base.rglob(INDEX_NAME)):
            d = idx.parent
            if sid_filter and d.parent.name != _safe(sid_filter):
                continue
            errs = verify_evidence_dir(d)
            problems.extend(errs)
            if not errs:
                verified.append(_rel(d, root))

    with_evidence = 0
    legacy: list[str] = []
    for n, row in enumerate(read_ledger(root, registry_dir), start=1):
        sid = str(row.get("source_id", ""))
        if sid_filter and sid != sid_filter:
            continue
        ref, digest = row.get("evidence_dir"), row.get("evidence_index_sha256")
        if not ref and not digest:
            if row.get("event_type") in REGISTRATION_EVENTS:
                legacy.append(str(row.get("entry_id")))
            continue
        with_evidence += 1
        tag = f"ledger line {n} ({row.get('entry_id')})"
        if not ref or not digest:
            problems.append(f"{tag}: evidence_dir and evidence_index_sha256 must both be set")
            continue
        d = root / str(ref)
        idx = d / INDEX_NAME
        if not idx.exists():
            problems.append(f"{tag}: evidence index {idx} missing")
            continue
        if sha256_file(idx) != digest:
            problems.append(f"{tag}: index.json SHA-256 differs from the ledger reference")
        if not base.exists() or not str(d.resolve()).startswith(str(base.resolve())):
            problems.append(f"{tag}: evidence directory outside {evidence_dir}")
        listed = {}
        try:
            listed = {f.path: f for f in EvidenceIndex.model_validate_json(idx.read_text(encoding="utf-8")).files}
        except Exception:  # noqa: BLE001
            problems.append(f"{tag}: index.json unreadable")
        rows = (row.get("details") or {}).get("vetting")
        if rows is None:
            problems.append(f"{tag}: registration with evidence must list its vetting rows in details.vetting")
            rows = []
        for vr in rows:
            if vr.get("status") not in ("pass", "fail"):
                continue
            ev = vr.get("evidence")
            if not ev or ev not in listed or not (d / ev).exists():
                problems.append(f"{tag}: {vr.get('status')} row '{vr.get('check')}' has no evidence file in {ref}")
    return {
        "problems": problems,
        "verified_dirs": verified,
        "ledger_entries_with_evidence": with_evidence,
        "legacy_without_evidence": len(legacy),
        "legacy_entry_ids": legacy,
    }
