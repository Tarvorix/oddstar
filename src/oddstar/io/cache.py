"""Content-addressed raw archive cache, manifest ledger, and 10 GB budget enforcer."""

from __future__ import annotations

import functools
import hashlib
import json
import os
import threading
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from oddstar.schemas import ManifestRecord


def _locked(method: Any) -> Any:
    """Run a RawCache method under the instance lock."""

    @functools.wraps(method)
    def wrapper(self: "RawCache", *args: Any, **kwargs: Any) -> Any:
        with self._lock:
            return method(self, *args, **kwargs)

    return wrapper


class RawCache:
    """Stores immutable raw archive blobs by SHA-256 and tracks manifest provenance."""

    def __init__(self, raw_dir: Path, max_cache_gb: float = 10.0) -> None:
        self.raw_dir = raw_dir
        self.max_cache_bytes = int(max_cache_gb * (1024**3))
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.raw_dir / "manifests.jsonl"
        self._index: dict[tuple[str, str, str], ManifestRecord] = {}
        self._by_sha: dict[str, ManifestRecord] = {}
        # Running byte total of the blobs on disk; None until first needed (one directory walk per process, not per put)
        self._cached_bytes: int | None = None
        self.recovered_truncated_manifest_line = False
        # Several fetch workers share one cache (survey NEOWISE ingest); index, byte total and manifest appends are
        # guarded so a manifest line is never interleaved and the budget count stays right.
        self._lock = threading.RLock()
        self._load_manifests()

    def _load_manifests(self) -> None:
        if not self.manifest_path.exists():
            return
        text = self.manifest_path.read_text(encoding="utf-8")
        lines = text.splitlines()
        for i, line in enumerate(lines):
            if not line.strip():
                continue
            try:
                rec = ManifestRecord.model_validate_json(line)
            except ValueError:
                # A run killed while appending leaves the last line cut short (no trailing newline). That request simply
                # was not recorded and is fetched again. A malformed line anywhere else is corruption and still raises.
                if i == len(lines) - 1 and not text.endswith("\n"):
                    self.recovered_truncated_manifest_line = True
                    continue
                raise
            key = (rec.archive, rec.release, rec.request_signature)
            self._index[key] = rec
            self._by_sha[rec.sha256] = rec
        if self.recovered_truncated_manifest_line:
            self._rewrite_manifests()

    def _rewrite_manifests(self) -> None:
        # One row per request: several requests can return byte-identical payloads (e.g. empty tables) and share a blob
        rows = [rec.model_dump_json() for rec in self._index.values()]
        tmp = self.manifest_path.with_suffix(".jsonl.tmp")
        tmp.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
        tmp.replace(self.manifest_path)

    @_locked
    def lookup(self, archive: str, release: str, request_signature: str) -> ManifestRecord | None:
        rec = self._index.get((archive, release, request_signature))
        if rec is None:
            return None
        blob = Path(rec.blob_path)
        if not rec.evicted and blob.exists():
            # Touch mtime for LRU tracking
            os.utime(blob, None)
            return rec
        return None

    def read_bytes(self, rec: ManifestRecord) -> bytes:
        return Path(rec.blob_path).read_bytes()

    @_locked
    def total_cached_bytes(self) -> int:
        total = 0
        for p in self.raw_dir.rglob("*"):
            if p.is_file() and p.name != "manifests.jsonl" and not p.name.endswith(".tmp"):
                total += p.stat().st_size
        self._cached_bytes = total
        return total

    def _mark_evicted(self, sha: str) -> None:
        for key, rec in list(self._index.items()):
            if rec.sha256 == sha and not rec.evicted:
                self._index[key] = rec.model_copy(update={"evicted": True})
        if sha in self._by_sha:
            self._by_sha[sha] = self._by_sha[sha].model_copy(update={"evicted": True})

    @_locked
    def enforce_budget(self, incoming_bytes: int = 0) -> list[str]:
        """Evict oldest raw blobs if cache exceeds max_cache_bytes; preserve manifest hashes."""
        evicted_shas: list[str] = []
        current = self._cached_bytes if self._cached_bytes is not None else self.total_cached_bytes()
        if current + incoming_bytes <= self.max_cache_bytes:
            return evicted_shas
        current = self.total_cached_bytes()  # re-measure from disk before evicting anything
        if current + incoming_bytes <= self.max_cache_bytes:
            return evicted_shas

        blobs = [
            p
            for p in self.raw_dir.rglob("*")
            if p.is_file() and p.name != "manifests.jsonl" and not p.name.endswith(".tmp")
        ]
        blobs.sort(key=lambda p: p.stat().st_mtime)
        for blob in blobs:
            if current + incoming_bytes <= self.max_cache_bytes:
                break
            size = blob.stat().st_size
            sha = blob.name
            blob.unlink(missing_ok=True)
            current -= size
            evicted_shas.append(sha)
            self._mark_evicted(sha)
        self._cached_bytes = current

        if evicted_shas:
            self._rewrite_manifests()
        return evicted_shas

    @_locked
    def evict_raw_blob(self, sha256: str) -> bool:
        """Process-and-evict helper for bulky radio/pixel payloads while retaining manifest provenance."""
        rec = self._by_sha.get(sha256)
        if rec is None:
            return False
        blob = Path(rec.blob_path)
        if blob.exists():
            size = blob.stat().st_size
            blob.unlink()
            if self._cached_bytes is not None:
                self._cached_bytes = max(0, self._cached_bytes - size)
        self._mark_evicted(sha256)
        self._rewrite_manifests()
        return True

    @_locked
    def put(
        self,
        *,
        archive: str,
        release: str,
        request_signature: str,
        query_text: str,
        payload: bytes,
        row_count: int,
        offline_fixture: bool = False,
    ) -> ManifestRecord:
        sha = hashlib.sha256(payload).hexdigest()
        existing = self.lookup(archive, release, request_signature)
        if existing is not None and existing.sha256 == sha and existing.offline_fixture == offline_fixture:
            return existing

        self.enforce_budget(incoming_bytes=len(payload))
        archive_dir = self.raw_dir / archive
        archive_dir.mkdir(parents=True, exist_ok=True)
        blob_path = archive_dir / sha
        if not blob_path.exists():
            blob_path.write_bytes(payload)
            if self._cached_bytes is not None:
                self._cached_bytes += len(payload)

        now = datetime.now(timezone.utc).isoformat()
        manifest_id = hashlib.sha256(f"{archive}:{release}:{request_signature}:{sha}".encode()).hexdigest()[:16]
        rec = ManifestRecord(
            manifest_id=manifest_id,
            archive=archive,
            release=release,
            request_signature=request_signature,
            query_text=query_text,
            utc_time=now,
            row_count=row_count,
            byte_size=len(payload),
            sha256=sha,
            blob_path=str(blob_path),
            evicted=False,
            offline_fixture=offline_fixture,
        )
        self._index[(archive, release, request_signature)] = rec
        self._by_sha[sha] = rec
        with self.manifest_path.open("a", encoding="utf-8") as fh:
            fh.write(rec.model_dump_json() + "\n")
        return rec

    def all_manifests(self) -> list[ManifestRecord]:
        return list(self._by_sha.values())

    def stats(self) -> dict[str, Any]:
        used_bytes = self.total_cached_bytes()
        return {
            "manifest_count": len(self._by_sha),
            "used_bytes": used_bytes,
            "used_gb": round(used_bytes / (1024**3), 6),
            "max_cache_gb": round(self.max_cache_bytes / (1024**3), 3),
        }
