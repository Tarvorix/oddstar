"""Content-addressed raw archive cache, manifest ledger, and 10 GB budget enforcer."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from oddstar.schemas import ManifestRecord


class RawCache:
    """Stores immutable raw archive blobs by SHA-256 and tracks manifest provenance."""

    def __init__(self, raw_dir: Path, max_cache_gb: float = 10.0) -> None:
        self.raw_dir = raw_dir
        self.max_cache_bytes = int(max_cache_gb * (1024**3))
        self.raw_dir.mkdir(parents=True, exist_ok=True)
        self.manifest_path = self.raw_dir / "manifests.jsonl"
        self._index: dict[tuple[str, str, str], ManifestRecord] = {}
        self._by_sha: dict[str, ManifestRecord] = {}
        self._load_manifests()

    def _load_manifests(self) -> None:
        if not self.manifest_path.exists():
            return
        for line in self.manifest_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rec = ManifestRecord.model_validate_json(line)
            key = (rec.archive, rec.release, rec.request_signature)
            self._index[key] = rec
            self._by_sha[rec.sha256] = rec

    def _rewrite_manifests(self) -> None:
        rows = [rec.model_dump_json() for rec in self._by_sha.values()]
        tmp = self.manifest_path.with_suffix(".jsonl.tmp")
        tmp.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
        tmp.replace(self.manifest_path)

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

    def total_cached_bytes(self) -> int:
        total = 0
        for p in self.raw_dir.rglob("*"):
            if p.is_file() and p.name != "manifests.jsonl" and not p.name.endswith(".tmp"):
                total += p.stat().st_size
        return total

    def enforce_budget(self, incoming_bytes: int = 0) -> list[str]:
        """Evict oldest raw blobs if cache exceeds max_cache_bytes; preserve manifest hashes."""
        evicted_shas: list[str] = []
        current = self.total_cached_bytes()
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
            if sha in self._by_sha:
                rec = self._by_sha[sha]
                updated = rec.model_copy(update={"evicted": True})
                self._by_sha[sha] = updated
                self._index[(rec.archive, rec.release, rec.request_signature)] = updated

        if evicted_shas:
            self._rewrite_manifests()
        return evicted_shas

    def evict_raw_blob(self, sha256: str) -> bool:
        """Process-and-evict helper for bulky radio/pixel payloads while retaining manifest provenance."""
        rec = self._by_sha.get(sha256)
        if rec is None:
            return False
        blob = Path(rec.blob_path)
        if blob.exists():
            blob.unlink()
        updated = rec.model_copy(update={"evicted": True})
        self._by_sha[sha256] = updated
        self._index[(rec.archive, rec.release, rec.request_signature)] = updated
        self._rewrite_manifests()
        return True

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
        if existing is not None and existing.sha256 == sha:
            return existing

        self.enforce_budget(incoming_bytes=len(payload))
        archive_dir = self.raw_dir / archive
        archive_dir.mkdir(parents=True, exist_ok=True)
        blob_path = archive_dir / sha
        if not blob_path.exists():
            blob_path.write_bytes(payload)

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
