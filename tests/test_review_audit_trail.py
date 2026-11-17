"""Review ingest idempotence, R10 model-family independence, human checklists, review status and CI checks."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

from oddstar.ci_checks import check_review_independence_wording, check_review_ingestion
from oddstar.cli import dispatch
from oddstar.config import repo_root
from oddstar.registry import read_ledger
from oddstar.review.status import model_family, parse_checklists, review_status, review_status_line
from oddstar.schemas import ReviewRecord

SID = "5628404318169472"


def _workspace(tmp: Path, *, with_records: bool = True) -> Path:
    src = repo_root()
    root = tmp / "ws"
    (root / "configs").mkdir(parents=True)
    shutil.copy(src / "configs" / "default.yaml", root / "configs" / "default.yaml")
    (root / "pyproject.toml").write_text("[project]\nname='ws'\n", encoding="utf-8")
    (root / "registry").mkdir()
    (root / "registry" / "ledger.jsonl").write_text("", encoding="utf-8")
    (root / "review" / "records").mkdir(parents=True)
    (root / "review" / "packets").mkdir(parents=True)
    if with_records:
        for p in (src / "review" / "records").glob("REV-20261117-*.json"):
            shutil.copy(p, root / "review" / "records" / p.name)
        shutil.copy(src / "review" / "adversarial-2026-11-17.md", root / "review" / "adversarial-2026-11-17.md")
    return root


def _synthetic_r3(root: Path, reviewer_model: str, review_id: str = "REV-20261201-R3B-5628404318169472") -> Path:
    rec = ReviewRecord(
        review_id=review_id, round=3, source_id=SID, reviewer_model=reviewer_model, timestamp_utc="2026-12-01T00:00:00+00:00",
        verdict="reproduced", notes="synthetic test record",
    )
    p = root.parent / f"{review_id}.json"
    p.write_text(rec.model_dump_json(indent=2), encoding="utf-8")
    return p


def test_model_family_normalisation() -> None:
    assert model_family("OpenAI GPT-6 family") == "openai-gpt"
    assert model_family("gpt-5.1") == "openai-gpt"
    assert model_family("Claude Opus 5") == "anthropic-claude"
    assert model_family("Gemini 3 Pro (Google)") == "google-gemini"
    assert model_family("Llama-4") == "meta-llama"
    assert model_family("human astronomer") == "human"
    # Placeholder labels are never treated as a known family
    assert model_family("Independent-Model-Family-R3") == "unrecognized:independent-model-family"
    assert model_family("Argon-Blind-R1") == "unrecognized:argon-blind"


def test_ingest_all_is_idempotent_and_never_rewrites_stored_records(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    before = {p.name: p.read_bytes() for p in (root / "review" / "records").glob("*.json")}
    code, pay = dispatch(["review", "ingest", "--all", "--json"], root=root)
    assert code == 0, pay
    rows = read_ledger(root)
    types = [r["event_type"] for r in rows]
    assert types.count("review_ingested") == 8 and types.count("finding_opened") == 4
    assert {r["details"]["finding_id"] for r in rows if r["event_type"] == "finding_opened"} == {
        "F-5628404318169472-20261117-R2-PRIOR-IR",
        "F-5628404318169472-20261117-R2-OPTICAL",
        "F-5628404318169472-20261117-R3-UNRESOLVED-BLEND",
        "F-1111149103659931264-20261117-R2-KNOWN-SUBDWARF",
    }
    r3 = next(r for r in rows if r["entry_id"] == f"REG-REV-20261117-R3-{SID}")
    assert r3["provenance"]["model_family"] == "openai-gpt" and r3["provenance"]["independent_family"] is False
    # Second run (and single-file re-ingest of a stored record) appends nothing and leaves records byte-identical
    assert dispatch(["review", "ingest", "--all", "--json"], root=root)[1]["data"]["appended_entries"] == []
    code, one = dispatch(["review", "ingest", str(root / "review" / "records" / f"REV-20261117-R2-{SID}.json"), "--json"], root=root)
    assert code == 0 and one["data"]["already_ingested"] is True
    assert len(read_ledger(root)) == len(rows)
    assert {p.name: p.read_bytes() for p in (root / "review" / "records").glob("*.json")} == before
    assert check_review_ingestion(root) == []


def test_conflicting_record_under_same_id_is_refused(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    stored = root / "review" / "records" / f"REV-20261117-R3-{SID}.json"
    altered = json.loads(stored.read_text(encoding="utf-8"))
    altered["verdict"] = "explained"
    other = tmp_path / "altered.json"
    other.write_text(json.dumps(altered), encoding="utf-8")
    code, err = dispatch(["review", "ingest", str(other), "--json"], root=root)
    assert code != 0 and "already stored with different content" in err["message"]


def test_round3_from_a_different_family_makes_the_review_independent(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    dispatch(["review", "ingest", "--all", "--json"], root=root)
    st = review_status(root, SID)
    assert st.independent_family is False and st.rounds == [1, 2, 3]
    assert st.summary == "3 rounds, one model family, human checklist 1 of 4"
    assert review_status_line(root, SID) == st.summary

    # A further round-3 record from the same family changes nothing
    code, same = dispatch(["review", "ingest", str(_synthetic_r3(root, "OpenAI GPT-6.5", "REV-20261201-R3A-5628404318169472")), "--json"], root=root)
    assert code == 0 and same["data"]["independent_family"] is False
    assert review_status(root, SID).independent_family is False

    code, pay = dispatch(["review", "ingest", str(_synthetic_r3(root, "Claude Opus 5")), "--json"], root=root)
    assert code == 0, pay
    assert pay["data"]["model_family"] == "anthropic-claude" and pay["data"]["independent_family"] is True
    ev = read_ledger(root)[-1]
    assert ev["event_type"] == "review_ingested" and ev["provenance"]["independent_family"] is True
    assert ev["provenance"]["earlier_round_families"] == ["openai-gpt"]
    st2 = review_status(root, SID)
    assert st2.independent_family is True and st2.distinct_families == ["anthropic-claude", "openai-gpt"]
    assert st2.summary == "3 rounds (5 records), two model families (round 3 independent), human checklist 1 of 4"
    code, status = dispatch(["status", "--json"], root=root)
    assert status["data"]["review_status"][SID]["independent_family"] is True


def test_unrecognised_family_labels_never_count_as_independent(tmp_path: Path) -> None:
    root = _workspace(tmp_path, with_records=False)
    for rnd, model in ((1, "Argon-Blind-R1"), (2, "Argon-Adversarial-R2"), (3, "Independent-Model-Family-R3")):
        rec = ReviewRecord(review_id=f"REV-X-R{rnd}-1", round=rnd, source_id="111111111111111111", reviewer_model=model,
                           timestamp_utc=f"2026-11-1{rnd}T00:00:00+00:00", verdict="reproduced")
        (root / "review" / "records" / f"{rec.review_id}.json").write_text(rec.model_dump_json(), encoding="utf-8")
    dispatch(["review", "ingest", "--all", "--json"], root=root)
    st = review_status(root, "111111111111111111")
    assert st.independent_family is False and "unverifiable" in st.independence_basis
    assert st.summary.startswith("3 rounds, model family unverified")


def test_human_checklists_are_counted_per_source() -> None:
    per = parse_checklists(repo_root())
    assert (per[SID]["ticked"], per[SID]["total"]) == (1, 4)
    assert (per["1111149103659931264"]["ticked"], per["1111149103659931264"]["total"]) == (1, 2)
    assert (per["4582145718700356096"]["ticked"], per["4582145718700356096"]["total"]) == (0, 2)
    assert "review/adversarial-2026-11-17.md" in per[SID]["files"]


def test_status_lists_open_findings_per_source(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    code, before = dispatch(["status", "--json"], root=root)
    assert all(f.endswith("(not ingested)") for f in before["data"]["open_findings_by_source"][SID])
    dispatch(["review", "ingest", "--all", "--json"], root=root)
    code, st = dispatch(["status", "--json"], root=root)
    by = st["data"]["open_findings_by_source"]
    assert by[SID] == [
        "F-5628404318169472-20261117-R2-OPTICAL",
        "F-5628404318169472-20261117-R2-PRIOR-IR",
        "F-5628404318169472-20261117-R3-UNRESOLVED-BLEND",
    ]
    assert by["1111149103659931264"] == ["F-1111149103659931264-20261117-R2-KNOWN-SUBDWARF"]
    assert any(w.startswith(f"review {SID}: 3 rounds, one model family") for w in st["warnings"])


def test_ci_review_ingestion_flags_uningested_records_and_unlessoned_closures(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    errs = check_review_ingestion(root)
    assert sum("no review_ingested ledger event" in e for e in errs) == 8
    assert sum("has no finding_opened ledger event" in e for e in errs) == 4
    dispatch(["review", "ingest", "--all", "--json"], root=root)
    assert check_review_ingestion(root) == []
    bad_close = {"entry_id": "REG-CLOSE-x", "timestamp_utc": "2026-11-18T00:00:00+00:00", "source_id": SID,
                 "event_type": "finding_closed", "status": "closed", "provenance": {},
                 "details": {"finding_id": "F-5628404318169472-20261117-R2-PRIOR-IR"}}
    with (root / "registry" / "ledger.jsonl").open("a", encoding="utf-8") as fh:
        fh.write(json.dumps(bad_close) + "\n")
    assert any("PRIOR-IR closed in the ledger" in e and "without a lesson" in e for e in check_review_ingestion(root))


def test_independent_review_wording_check(tmp_path: Path) -> None:
    root = _workspace(tmp_path)
    dispatch(["review", "ingest", "--all", "--json"], root=root)
    (root / "docs").mkdir()
    (root / "docs" / "log.md").write_text(
        "# Log\n\n"
        f"Gaia DR3 {SID} stays at Tier C after three rounds of independent adversarial review.\n\n"
        f"- Gaia DR3 {SID}: the review was not independent review by a second family.\n"
        f"- VarWISE independently catalogs the variability of {SID}.\n\n"
        "Eight independent channels search the sky; independent review of other stars is planned.\n",
        encoding="utf-8",
    )
    errs = check_review_independence_wording(root)
    assert len(errs) == 1 and errs[0].startswith("docs/log.md:3:") and SID in errs[0]
    dispatch(["review", "ingest", str(_synthetic_r3(root, "Gemini 3")), "--json"], root=root)
    assert check_review_independence_wording(root) == []
