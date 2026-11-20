"""`vet --live --fetch-workers`: flagged stars are fetched a few at a time and vetted in the order of the search."""

from __future__ import annotations

import json
import shutil
import threading
import time
from pathlib import Path
from typing import Any

import pytest

from oddstar import cli
from oddstar.cli import _OrderedFetch, dispatch
from oddstar.config import repo_root
from oddstar.errors import ArchiveAccessError
from oddstar.io.cache import RawCache

RUN_ID = "test-vet-fetch-workers"
IDS = ["901", "902", "903", "904", "905", "906"]


def test_ordered_fetch_runs_ahead_and_keeps_the_callers_order() -> None:
    running, peak, lock = [0], [0], threading.Lock()

    def fetch(key: str) -> str:
        with lock:
            running[0] += 1
            peak[0] = max(peak[0], running[0])
        time.sleep(0.05)
        with lock:
            running[0] -= 1
        if key == "bad":
            raise ArchiveAccessError("archive down", retryable=True)
        return key.upper()

    keys = ["a", "b", "c", "d", "a"]
    f = _OrderedFetch(fetch, keys, 3)
    try:
        assert [f.get(k) for k in keys] == ["A", "B", "C", "D", "A"]
    finally:
        f.close()
    assert peak[0] == 3  # three at a time, never more

    # one worker, or one key: nothing runs ahead, each get fetches on the spot
    calls: list[str] = []
    one = _OrderedFetch(lambda k: calls.append(k) or k, ["x", "y"], 1)
    assert calls == [] and one.get("y") == "y" and calls == ["y"]
    one.close()

    # a failed fetch raises at its own place in the order, after the keys before it were handed out
    g = _OrderedFetch(fetch, ["a", "bad", "c"], 3)
    try:
        assert g.get("a") == "A"
        with pytest.raises(ArchiveAccessError, match="archive down"):
            g.get("bad")
    finally:
        g.close()


def _write_run(run_dir: Path) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    scores = [{"source_id": sid, "channel": "optical_dips", "raw_statistic": 9.0, "false_alarm_prob": 1e-6, "flagged": True,
               "rank_score": 6.0, "summary": "dip", "feature_contributions": {}} for sid in IDS]
    (run_dir / "search_optical_dips.json").write_text(json.dumps({
        "run_id": RUN_ID, "channel": "optical_dips", "scores": scores, "measurements": [], "flagged_ids": IDS}), encoding="utf-8")


def _vet(workers: int, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, fail: str | None = None) -> tuple[int, dict[str, Any], dict[str, Any], list[str]]:
    threads: list[str] = []

    def fake_fetch(sid: str, cfg: Any, cache: Any) -> tuple[dict[str, Any], list[Any]]:
        threads.append(threading.current_thread().name)
        time.sleep(0.02)
        if sid == fail:
            raise ArchiveAccessError(f"Gaia DR3 {sid} not reachable", retryable=True)
        return {"source_id": sid, "ra": 10.0 + int(sid) / 1000.0, "dec": 5.0, "live_enriched": True, "vetting_archive": {"query_status": {}}}, []

    seen: list[str] = []

    def fake_vet_target(star: dict[str, Any], scores: list[Any], cfg: Any, ev_dir: Path) -> list[Any]:
        seen.append(star["source_id"])
        star["touched"] = True  # a check that writes into its star must not reach another star's copy
        return []

    monkeypatch.setattr(cli, "fetch_single_star_live", fake_fetch)
    monkeypatch.setattr(cli, "vet_target", fake_vet_target)
    monkeypatch.setattr(cli, "RawCache", lambda *a, **kw: RawCache(tmp_path / f"raw{workers}"))
    run_dir = repo_root() / "runs" / RUN_ID
    _write_run(run_dir)
    code, payload = dispatch(["vet", RUN_ID, "--live", "--fetch-workers", str(workers), "--json"])
    out = json.loads((run_dir / "vetting.json").read_text(encoding="utf-8")) if (run_dir / "vetting.json").exists() else {}
    return code, payload, out, seen + ["|"] + threads


def test_vet_fetches_several_stars_at_once_and_vets_them_in_order(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    run_dir = repo_root() / "runs" / RUN_ID
    try:
        code, payload, out3, log3 = _vet(3, monkeypatch, tmp_path)
        assert code == 0, payload
        vetted, threads = log3[: log3.index("|")], log3[log3.index("|") + 1:]
        assert vetted == IDS and payload["data"]["vetted_count"] == len(IDS)
        assert len(threads) == len(IDS) and all(t.startswith("vetfetch") for t in threads) and len(set(threads)) > 1
        shutil.rmtree(run_dir)

        code, payload, out1, log1 = _vet(1, monkeypatch, tmp_path)
        assert code == 0, payload
        vetted, threads = log1[: log1.index("|")], log1[log1.index("|") + 1:]
        assert vetted == IDS and not any(t.startswith("vetfetch") for t in threads)
        assert out1 == out3  # the same vetting whatever the number of workers
        shutil.rmtree(run_dir)

        # a star that cannot be fetched stops the vetting with its own error, as it did one star at a time
        code, payload, out, log = _vet(3, monkeypatch, tmp_path, fail="903")
        assert code != 0 and "903" in payload["message"] and out == {}
        assert log[: log.index("|")] == ["901", "902"]
    finally:
        shutil.rmtree(run_dir, ignore_errors=True)
