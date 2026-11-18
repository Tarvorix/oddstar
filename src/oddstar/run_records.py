"""Committed run records: what each search or fusion run covered and found, small enough to live in the repository.

`runs/` is git-ignored (it holds per-star measurements and copies of the sample), so a statement such as "7,800
stars searched, one survivor" could not be checked against anything in the repository. A run record is the part of
a run that such statements rest on:

- `registry/runs/<run_id>.json`: for each search in the run its channel, release, config hash, git commit, seed,
  certificate, sample cuts, how many stars were in the sample / searched / not searched (and why) / flagged, the
  flagged stars with their statistics, the vetting table of each flagged star, the null-report numbers, and the
  SHA-256 of the full local search payload. For a fusion run: the input runs, tier counts and candidate rows.
- `registry/samples/<sha16>.json`: the sorted source_ids of a sample, stored once and shared by every run that
  searched the same stars.

Records contain no raw archive data and no per-star measurements. `verify_run_records` re-checks every record
against its sample file and is part of the CI governance checks.
"""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from oddstar.errors import ArtifactNotFoundError, PolicyViolationError
from oddstar.evidence import canonical_json, sha256_bytes, sha256_file

RECORD_VERSION = 1
RUNS_SUBDIR = "runs"
SAMPLES_SUBDIR = "samples"
MAX_NOTE_CHARS = 400
RUN_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,120}$")
COVERAGE_KEYS = ("n_input", "n_searched", "n_not_searched", "n_no_coverage", "n_failed", "n_not_fetched", "n_insufficient_visits",
                 "n_no_light_curve", "n_single_epoch", "n_no_spectrum", "n_too_faint", "n_near_bright_star",
                 "n_no_plate_record", "n_not_downloaded", "n_no_comparison")
NULL_REPORT_KEYS = ("n_searched", "n_flagged", "expected_false_alarms", "tier_breakdown", "upper_limits", "contour_50", "contour_90")
FLAGGED_KEYS = ("source_id", "raw_statistic", "false_alarm_prob", "rank_score", "summary")
FUSION_ROW_KEYS = ("source_id", "tier", "rank_score", "combined_fap", "flagged_channels", "patterns", "open_holes", "explained_reason")


def source_ids_sha256(source_ids: list[str]) -> str:
    """SHA-256 of the sorted, newline-joined source_ids (the identity of a sample)."""
    return sha256_bytes("\n".join(sorted(str(s) for s in source_ids)).encode("utf-8"))


def _dump(obj: Any) -> str:
    return json.dumps(obj, indent=2, sort_keys=True, ensure_ascii=True) + "\n"


def _clip(text: Any) -> str:
    s = " ".join(str(text or "").split())
    return s if len(s) <= MAX_NOTE_CHARS else s[: MAX_NOTE_CHARS - 1] + "…"


def _sample_ids(sdata: dict[str, Any]) -> list[str]:
    if "sample_stars" in sdata:
        return [str(s["source_id"]) for s in sdata["sample_stars"]]
    return [str(r["source_id"]) for r in sdata.get("scores", [])]


def _search_entry(search_file: Path, run_dir: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    sdata = json.loads(search_file.read_text(encoding="utf-8"))
    channel = str(sdata["channel"])
    ids = _sample_ids(sdata)
    ids_sha = source_ids_sha256(ids)
    sample = {"record_version": RECORD_VERSION, "n": len(ids), "source_ids_sha256": ids_sha, "source_ids": sorted(ids)}
    flagged = [{k: row.get(k) for k in FLAGGED_KEYS} for row in sdata.get("scores", []) if row.get("flagged")]
    flagged.sort(key=lambda r: str(r["source_id"]))
    entry: dict[str, Any] = {
        "channel": channel,
        "gaia_release": sdata.get("gaia_release"),
        "config_hash": sdata.get("config_hash"),
        "git_commit": sdata.get("git_commit"),
        "seed": sdata.get("seed"),
        "certificate": sdata.get("certificate"),
        "cuts": sdata.get("cuts", []),
        "sample": {"n": len(ids), "source_ids_sha256": ids_sha, "file": f"{SAMPLES_SUBDIR}/{ids_sha[:16]}.json"},
        "n_searched": int(sdata["n_searched"]),
        "n_flagged": int(sdata["n_flagged"]),
        "flagged": flagged,
        "search_payload_sha256": sha256_file(search_file),
    }
    if "n_not_searched" in sdata:
        entry["coverage"] = {k: sdata.get(k) for k in COVERAGE_KEYS}
        # The search payload groups these ids by reason ({category: [source_id, ...]}); the record keeps one flat list
        nsi = sdata.get("not_searched_ids", [])
        flat = [x for ids_ in nsi.values() for x in ids_] if isinstance(nsi, dict) else list(nsi)
        entry["not_searched_ids"] = sorted(str(s) for s in flat)
        for key in ("neowise_fetch", "tess_fetch", "spectra_fetch", "plate_fetch", "radio_fetch"):
            nf = sdata.get(key)
            if isinstance(nf, dict):
                entry[key] = nf
    vf = run_dir / "vetting.json"
    if vf.exists():
        vdata = json.loads(vf.read_text(encoding="utf-8")).get("vetting", {})
        flagged_ids = {str(r["source_id"]) for r in flagged}
        entry["vetting"] = {
            sid: [{"check": v.get("check"), "status": v.get("status"), "note": _clip(v.get("note"))} for v in rows]
            for sid, rows in sorted(vdata.items())
            if sid in flagged_ids
        }
        entry["vetting_payload_sha256"] = sha256_file(vf)
    nr = run_dir / f"null_report_{channel}.json"
    if nr.exists():
        ndata = json.loads(nr.read_text(encoding="utf-8"))
        entry["null_report"] = {k: ndata.get(k) for k in NULL_REPORT_KEYS}
    return entry, sample


def build_run_record(run_dir: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """The committed record of the run in `run_dir` and the sample records it refers to."""
    run_id = run_dir.name
    if not RUN_ID_RE.match(run_id):
        raise PolicyViolationError(f"Run id {run_id!r} is not a safe file name", hint="Use letters, digits, '.', '_' and '-' only")
    if not run_dir.is_dir():
        raise ArtifactNotFoundError(f"Run directory not found: {run_dir}")
    record: dict[str, Any] = {"record_version": RECORD_VERSION, "run_id": run_id}
    state_p = run_dir / "run_state.json"
    if state_p.exists():
        state = json.loads(state_p.read_text(encoding="utf-8"))
        record["timestamp_utc"] = state.get("timestamp_utc")
    samples: dict[str, dict[str, Any]] = {}
    searches: list[dict[str, Any]] = []
    for sf in sorted(run_dir.glob("search_*.json")):
        entry, sample = _search_entry(sf, run_dir)
        searches.append(entry)
        samples[sample["source_ids_sha256"]] = sample
    if searches:
        record["searches"] = searches
    fp = run_dir / "fusion.json"
    if fp.exists():
        fdata = json.loads(fp.read_text(encoding="utf-8"))
        inputs_p = run_dir / "fuse_inputs.json"
        cands = [{k: c.get(k) for k in FUSION_ROW_KEYS} for c in fdata.get("candidates", [])]
        tier_counts: dict[str, int] = {}
        for c in cands:
            tier_counts[str(c["tier"])] = tier_counts.get(str(c["tier"]), 0) + 1
        record["fusion"] = {
            "config_hash": fdata.get("config_hash"),
            "git_commit": fdata.get("git_commit"),
            "timestamp_utc": fdata.get("timestamp_utc"),
            "input_run_ids": json.loads(inputs_p.read_text(encoding="utf-8")).get("input_run_ids", []) if inputs_p.exists() else [],
            "n_stars_total": fdata.get("n_stars_total"),
            "n_trials_total": fdata.get("n_trials_total"),
            "expected_false_positives": fdata.get("expected_false_positives"),
            "consistent_with_noise": fdata.get("consistent_with_noise"),
            "tier_counts": dict(sorted(tier_counts.items())),
            "candidates": cands,
            "fusion_payload_sha256": sha256_file(fp),
        }
    if "searches" not in record and "fusion" not in record:
        raise ArtifactNotFoundError(f"Run {run_id} has no search_*.json and no fusion.json to record")
    return record, list(samples.values())


def export_run_record(root: Path, run_id: str, *, runs_dir: str = "runs", registry_dir: str = "registry", replace: bool = False) -> dict[str, Any]:
    """Write `registry/runs/<run_id>.json` and the sample files it needs. Returns the paths and the record.

    An existing record with different content is not overwritten unless `replace` (a committed record is cited by
    documents; a run that changed should normally get a new run id).
    """
    record, samples = build_run_record(root / runs_dir / run_id)
    reg = root / registry_dir
    rec_path = reg / RUNS_SUBDIR / f"{run_id}.json"
    text = _dump(record)
    if rec_path.exists() and rec_path.read_text(encoding="utf-8") != text and not replace:
        raise PolicyViolationError(
            f"{rec_path.relative_to(root)} exists with different content",
            hint="Use a new run id, or pass --replace if the record was never committed or the run gained vetting / a null report",
        )
    written: list[str] = []
    for sample in samples:
        sp = reg / SAMPLES_SUBDIR / f"{sample['source_ids_sha256'][:16]}.json"
        stext = _dump(sample)
        if sp.exists():
            if sp.read_text(encoding="utf-8") != stext:
                raise PolicyViolationError(f"{sp.relative_to(root)} exists with different content (hash prefix collision or edited file)")
        else:
            sp.parent.mkdir(parents=True, exist_ok=True)
            sp.write_text(stext, encoding="utf-8")
            written.append(str(sp.relative_to(root)))
    rec_path.parent.mkdir(parents=True, exist_ok=True)
    rec_path.write_text(text, encoding="utf-8")
    written.append(str(rec_path.relative_to(root)))
    return {"record": record, "record_path": str(rec_path.relative_to(root)), "written": written,
            "record_sha256": sha256_bytes(canonical_json(record))}


def load_run_records(root: Path, registry_dir: str = "registry") -> list[dict[str, Any]]:
    """Every committed run record, sorted by run id."""
    d = root / registry_dir / RUNS_SUBDIR
    return [json.loads(p.read_text(encoding="utf-8")) for p in sorted(d.glob("*.json"))] if d.is_dir() else []


def verify_run_records(root: Path, registry_dir: str = "registry") -> dict[str, Any]:
    """Check every committed run record against its sample files and its own counts."""
    reg = root / registry_dir
    problems: list[str] = []
    n_records = 0
    used_samples: set[str] = set()
    runs_d = reg / RUNS_SUBDIR
    for p in sorted(runs_d.glob("*.json")) if runs_d.is_dir() else []:
        rel = p.relative_to(root)
        n_records += 1
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
        except ValueError as exc:
            problems.append(f"{rel}: not valid JSON ({exc})")
            continue
        if rec.get("record_version") != RECORD_VERSION:
            problems.append(f"{rel}: record_version {rec.get('record_version')!r} (expected {RECORD_VERSION})")
        if rec.get("run_id") != p.stem:
            problems.append(f"{rel}: run_id {rec.get('run_id')!r} differs from the file name")
        if "searches" not in rec and "fusion" not in rec:
            problems.append(f"{rel}: neither a search nor a fusion record")
        for s in rec.get("searches", []):
            ch = s.get("channel")
            smp = s.get("sample") or {}
            if int(s.get("n_flagged", -1)) != len(s.get("flagged", [])):
                problems.append(f"{rel} [{ch}]: n_flagged {s.get('n_flagged')} but {len(s.get('flagged', []))} flagged rows")
            if int(s.get("n_searched", 0)) > int(smp.get("n", 0)):
                problems.append(f"{rel} [{ch}]: n_searched {s.get('n_searched')} exceeds the sample size {smp.get('n')}")
            cov = s.get("coverage")
            if cov is not None:
                if int(cov.get("n_searched", -1)) != int(s.get("n_searched", -2)):
                    problems.append(f"{rel} [{ch}]: coverage.n_searched differs from n_searched")
                if int(cov.get("n_not_searched", -1)) != len(s.get("not_searched_ids", [])):
                    problems.append(f"{rel} [{ch}]: coverage.n_not_searched differs from the number of not_searched_ids")
                if int(cov.get("n_searched", 0)) + int(cov.get("n_not_searched", 0)) != int(smp.get("n", -1)):
                    problems.append(f"{rel} [{ch}]: searched + not searched differs from the sample size")
            vet = s.get("vetting")
            if vet is not None:
                extra = sorted(set(vet) - {str(r.get("source_id")) for r in s.get("flagged", [])})
                if extra:
                    problems.append(f"{rel} [{ch}]: vetting rows for stars that were not flagged: {extra[:5]}")
            sp = reg / str(smp.get("file", ""))
            if not smp.get("file") or not sp.is_file():
                problems.append(f"{rel} [{ch}]: sample file {smp.get('file')!r} is missing")
                continue
            used_samples.add(sp.name)
            try:
                sample = json.loads(sp.read_text(encoding="utf-8"))
            except ValueError as exc:
                problems.append(f"{sp.relative_to(root)}: not valid JSON ({exc})")
                continue
            ids = [str(x) for x in sample.get("source_ids", [])]
            if source_ids_sha256(ids) != smp.get("source_ids_sha256") or sample.get("source_ids_sha256") != smp.get("source_ids_sha256"):
                problems.append(f"{rel} [{ch}]: sample file {smp.get('file')} does not hash to {str(smp.get('source_ids_sha256'))[:16]}")
            if len(ids) != int(smp.get("n", -1)) or len(set(ids)) != len(ids):
                problems.append(f"{rel} [{ch}]: sample file holds {len(ids)} ids ({len(set(ids))} distinct), record says {smp.get('n')}")
            idset = set(ids)
            stray = [r.get("source_id") for r in s.get("flagged", []) if str(r.get("source_id")) not in idset]
            stray += [x for x in s.get("not_searched_ids", []) if str(x) not in idset]
            if stray:
                problems.append(f"{rel} [{ch}]: {len(stray)} flagged or not-searched ids are not in the sample (e.g. {stray[0]})")
        fus = rec.get("fusion")
        if fus is not None:
            counts: dict[str, int] = {}
            for c in fus.get("candidates", []):
                counts[str(c.get("tier"))] = counts.get(str(c.get("tier")), 0) + 1
            if counts != {k: int(v) for k, v in (fus.get("tier_counts") or {}).items()}:
                problems.append(f"{rel}: fusion tier_counts differ from the candidate rows")
    samples_d = reg / SAMPLES_SUBDIR
    orphans = sorted(p.name for p in samples_d.glob("*.json") if p.name not in used_samples) if samples_d.is_dir() else []
    for name in orphans:
        problems.append(f"{registry_dir}/{SAMPLES_SUBDIR}/{name}: not referenced by any run record")
    return {"n_records": n_records, "n_samples": len(used_samples), "problems": problems}
