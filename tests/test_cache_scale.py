"""Raw cache behaviour that matters for survey-scale runs: resume after a killed run, no per-put directory walk."""

from __future__ import annotations

from pathlib import Path

import pytest

from oddstar.io.cache import RawCache


def _put(cache: RawCache, sig: str, payload: bytes) -> None:
    cache.put(archive="neowise", release="r", request_signature=sig, query_text=sig, payload=payload, row_count=1)


def test_truncated_last_manifest_line_is_dropped_and_refetchable(tmp_path: Path) -> None:
    cache = RawCache(tmp_path)
    _put(cache, "a", b"aaaa")
    _put(cache, "b", b"bbbb")
    text = cache.manifest_path.read_text(encoding="utf-8")
    cache.manifest_path.write_text(text[: len(text) - 40], encoding="utf-8")  # killed mid-append: last line cut short

    resumed = RawCache(tmp_path)
    assert resumed.recovered_truncated_manifest_line is True
    assert resumed.lookup("neowise", "r", "a") is not None
    assert resumed.lookup("neowise", "r", "b") is None  # not recorded, so it is fetched again
    _put(resumed, "b", b"bbbb")
    assert RawCache(tmp_path).lookup("neowise", "r", "b") is not None


def test_corrupt_line_in_the_middle_still_raises(tmp_path: Path) -> None:
    cache = RawCache(tmp_path)
    _put(cache, "a", b"aaaa")
    _put(cache, "b", b"bbbb")
    lines = cache.manifest_path.read_text(encoding="utf-8").splitlines()
    cache.manifest_path.write_text(lines[0][:30] + "\n" + lines[1] + "\n", encoding="utf-8")
    with pytest.raises(ValueError):
        RawCache(tmp_path)


def test_put_does_not_walk_the_cache_each_time(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    cache = RawCache(tmp_path)
    calls = {"n": 0}
    original = cache.total_cached_bytes

    def counting() -> int:
        calls["n"] += 1
        return original()

    monkeypatch.setattr(cache, "total_cached_bytes", counting)
    for i in range(50):
        _put(cache, f"s{i}", f"payload-{i}".encode())
    assert calls["n"] == 1
    assert cache._cached_bytes == sum(len(f"payload-{i}".encode()) for i in range(50))


def test_identical_payloads_keep_one_manifest_row_per_request(tmp_path: Path) -> None:
    cache = RawCache(tmp_path, max_cache_gb=1e-9)  # one byte: every put evicts what came before
    _put(cache, "empty-1", b"no rows")
    _put(cache, "empty-2", b"no rows")
    _put(cache, "other", b"something else")
    sigs = {key[2] for key in RawCache(tmp_path)._index}  # every request keeps its provenance row on disk
    assert {"empty-1", "empty-2", "other"} <= sigs
