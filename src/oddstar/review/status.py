"""Honest review status per source (rules R5, R7, R10).

- `model_family` normalises a free-text `reviewer_model` to a family label. Only names that contain a known vendor or
  model-family token are recognised; anything else becomes `unrecognized:<slug>`, and independence involving an
  unrecognised family is never asserted.
- `parse_checklists` counts ticked / unticked items of the human checklists (`review/packets/*_human_checklist.md`
  and the `## Human checklist` section of `review/adversarial-*.md`) per source.
- `review_status` combines stored records, ledger events and checklists into a `ReviewStatus`, and
  `review_status_line` returns the one-line summary that documents must use (for example
  "3 rounds, one model family, human checklist 1 of 4").
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

from oddstar.registry import clean_source_id, read_ledger
from oddstar.schemas import ReviewRecord, ReviewStatus

# (family label, pattern) - first match wins. Patterns match whole tokens of the lower-cased reviewer_model.
_FAMILY_PATTERNS: list[tuple[str, str]] = [
    ("openai-gpt", r"\b(openai|gpt[-\s]?\d*(\.\d+)?|chatgpt|codex)\b"),
    ("anthropic-claude", r"\b(anthropic|claude)\b"),
    ("google-gemini", r"\b(google|gemini|gemma|palm|bard)\b"),
    ("meta-llama", r"\b(meta|llama)\b"),
    ("mistral", r"\b(mistral|mixtral|codestral)\b"),
    ("xai-grok", r"\b(xai|grok)\b"),
    ("deepseek", r"\bdeepseek\b"),
    ("alibaba-qwen", r"\bqwen\b"),
    ("human", r"\b(human|person|astronomer)\b"),
]
SOURCE_ID_RE = re.compile(r"\b\d{15,20}\b")
_ITEM_RE = re.compile(r"^\s*[-*]\s+\[([ xX])\]\s+(.*)$")
_NUMBER_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}


def model_family(reviewer_model: str) -> str:
    """Normalised model family of a `reviewer_model` string ("OpenAI GPT-6 family" -> "openai-gpt")."""
    text = (reviewer_model or "").strip().lower()
    for label, pat in _FAMILY_PATTERNS:
        if re.search(pat, text):
            return label
    slug = re.sub(r"[^a-z0-9]+", "-", text).strip("-")
    slug = re.sub(r"-r\d+$", "", slug)
    return f"unrecognized:{slug or 'empty'}"


def family_recognized(family: str) -> bool:
    return not family.startswith("unrecognized:")


def independence(records: list[ReviewRecord]) -> tuple[bool, str, dict[str, str]]:
    """Whether a round-3 record comes from a recognised family absent from rounds 1-2 (all recognised)."""
    fams = {f"R{r.round}:{r.review_id}": model_family(r.reviewer_model) for r in records}
    early = {model_family(r.reviewer_model) for r in records if r.round in (1, 2)}
    late = [model_family(r.reviewer_model) for r in records if r.round == 3]
    if not late:
        return False, "no round-3 record", fams
    if not early:
        return False, "no round-1/2 record to compare with", fams
    unknown = sorted(f for f in early | set(late) if not family_recognized(f))
    if unknown:
        return False, f"unverifiable: unrecognized model family {', '.join(unknown)}", fams
    new = sorted(set(late) - early)
    if new:
        return True, f"round 3 family {', '.join(new)} differs from rounds 1-2 ({', '.join(sorted(early))})", fams
    return False, f"round 3 uses the rounds 1-2 family ({', '.join(sorted(early))})", fams


# ----------------------------------------------------------------------------------------------- checklists


def _aliases_from_headings(text: str) -> dict[str, set[str]]:
    """`## <source_id> / alias / alias` headings -> {source_id: {aliases}}."""
    out: dict[str, set[str]] = {}
    for line in text.splitlines():
        m = re.match(r"^#{2,3}\s+(\d{15,20})\b(.*)$", line.strip())
        if not m:
            continue
        sid = m.group(1)
        names = {p.strip() for p in m.group(2).split("/") if p.strip()}
        names |= {n.split()[-1] for n in names if n.startswith("AllWISE ")}
        out.setdefault(sid, set()).update(names)
    return out


def _parse_items(lines: list[str]) -> list[tuple[bool, str]]:
    items = []
    for line in lines:
        m = _ITEM_RE.match(line)
        if m:
            items.append((m.group(1).lower() == "x", m.group(2)))
    return items


def _attribute(items: list[tuple[bool, str]], aliases: dict[str, set[str]], doc_sources: list[str]) -> list[tuple[bool, list[str]]]:
    """Explicit source ID / alias in the item; otherwise the previous item's sources; otherwise every doc source."""
    out: list[tuple[bool, list[str]]] = []
    prev: list[str] | None = None
    for ticked, text in items:
        found = sorted(set(SOURCE_ID_RE.findall(text)))
        for sid, names in aliases.items():
            if any(n and n in text for n in names):
                found.append(sid)
        found = sorted(set(found))
        if not found:
            found = list(prev) if prev else list(doc_sources)
        prev = found
        out.append((ticked, found))
    return out


def parse_checklists(root: Path) -> dict[str, dict[str, Any]]:
    """Per source: ticked/total human checklist items and the files they come from."""
    per: dict[str, dict[str, Any]] = {}

    def _add(sid: str, ticked: bool, rel: str) -> None:
        d = per.setdefault(sid, {"ticked": 0, "total": 0, "files": []})
        d["total"] += 1
        d["ticked"] += int(ticked)
        if rel not in d["files"]:
            d["files"].append(rel)

    pkt_dir = root / "review" / "packets"
    for p in sorted(pkt_dir.glob("*_human_checklist.md")) if pkt_dir.exists() else []:
        rel = str(p.relative_to(root))
        current: list[str] = []
        for line in p.read_text(encoding="utf-8").splitlines():
            if line.startswith("#"):
                ids = SOURCE_ID_RE.findall(line)
                if ids:
                    current = ids
                continue
            m = _ITEM_RE.match(line)
            if m and current:
                for sid in current:
                    _add(sid, m.group(1).lower() == "x", rel)

    rev_dir = root / "review"
    for p in sorted(rev_dir.glob("adversarial-*.md")) if rev_dir.exists() else []:
        rel = str(p.relative_to(root))
        text = p.read_text(encoding="utf-8")
        aliases = _aliases_from_headings(text)
        doc_sources = sorted(aliases)
        section: list[str] = []
        inside = False
        for line in text.splitlines():
            if re.match(r"^##\s+", line):
                inside = bool(re.match(r"^##\s+human checklist\b", line.strip(), flags=re.IGNORECASE))
                continue
            if inside:
                section.append(line)
        for ticked, sids in _attribute(_parse_items(section), aliases, doc_sources):
            for sid in sids:
                _add(sid, ticked, rel)
    return per


# ------------------------------------------------------------------------------------------------ status


def load_review_records(root: Path) -> list[tuple[Path, ReviewRecord]]:
    rec_dir = root / "review" / "records"
    out = []
    for p in sorted(rec_dir.glob("*.json")) if rec_dir.exists() else []:
        out.append((p, ReviewRecord.model_validate_json(p.read_text(encoding="utf-8"))))
    return out


def open_findings_by_source(rows: list[dict[str, Any]], records: list[tuple[Path, ReviewRecord]] | None = None) -> dict[str, list[str]]:
    """Findings opened in the ledger without a later finding_closed, plus open record findings never opened there."""
    state: dict[str, tuple[str, bool]] = {}
    for row in rows:
        et = row.get("event_type")
        fid = (row.get("details") or {}).get("finding_id")
        if not fid:
            continue
        if et == "finding_opened":
            state[fid] = (str(row.get("source_id")), True)
        elif et == "finding_closed":
            sid = state.get(fid, (str(row.get("source_id")), False))[0]
            state[fid] = (sid, False)
    out: dict[str, list[str]] = {}
    for fid, (sid, is_open) in state.items():
        if is_open:
            out.setdefault(sid, []).append(fid)
    for _, rec in records or []:
        sid = clean_source_id(rec.source_id)
        for f in rec.findings:
            if f.status == "open" and f.finding_id not in state:
                out.setdefault(sid, []).append(f"{f.finding_id} (not ingested)")
    return {k: sorted(v) for k, v in out.items()}


def _count_word(n: int) -> str:
    return _NUMBER_WORDS.get(n, str(n))


def summary_line(st: ReviewStatus) -> str:
    rounds = len(set(st.rounds))
    if not st.rounds:
        head = "no review rounds"
    else:
        head = f"{rounds} round{'s' if rounds != 1 else ''}"
        if len(st.rounds) != rounds:
            head += f" ({len(st.rounds)} records)"
    fams = st.distinct_families
    if not fams:
        fam = "no model family recorded"
    elif any(not family_recognized(f) for f in fams):
        fam = "model family unverified"
    elif len(fams) == 1:
        fam = "one model family"
    else:
        fam = f"{_count_word(len(fams))} model families" + (
            " (round 3 independent)" if st.independent_family else " (round 3 not independent)"
        )
    chk = f"human checklist {st.checklist_ticked} of {st.checklist_total}" if st.checklist_total else "no human checklist"
    return f"{head}, {fam}, {chk}"


def review_status(root: Path, source_id: str, *, registry_dir: str = "registry",
                  _cache: dict[str, Any] | None = None) -> ReviewStatus:
    sid = clean_source_id(source_id)
    cache = _cache if _cache is not None else {}
    if "records" not in cache:
        cache["records"] = load_review_records(root)
        cache["rows"] = read_ledger(root, registry_dir)
        cache["checklists"] = parse_checklists(root)
    recs = [r for _, r in cache["records"] if clean_source_id(r.source_id) == sid]
    rows = cache["rows"]
    ingested = {str(r.get("entry_id")) for r in rows if r.get("event_type") == "review_ingested"}
    indep, basis, fams = independence(recs)
    chk = cache["checklists"].get(sid, {"ticked": 0, "total": 0, "files": []})
    st = ReviewStatus(
        source_id=sid,
        rounds=sorted(r.round for r in recs),
        model_families=fams,
        distinct_families=sorted({model_family(r.reviewer_model) for r in recs}),
        independent_family=indep,
        independence_basis=basis,
        records_ingested=sum(1 for r in recs if f"REG-{r.review_id}" in ingested),
        records_total=len(recs),
        checklist_ticked=int(chk["ticked"]),
        checklist_total=int(chk["total"]),
        checklist_files=list(chk["files"]),
        open_findings=open_findings_by_source(rows, cache["records"]).get(sid, []),
    )
    return st.model_copy(update={"summary": summary_line(st)})


def review_status_line(root: Path, source_id: str, *, registry_dir: str = "registry") -> str:
    """The honest one-line review status of a source; document generation must use this wording."""
    return review_status(root, source_id, registry_dir=registry_dir).summary


def reviewed_source_ids(root: Path) -> list[str]:
    return sorted({clean_source_id(r.source_id) for _, r in load_review_records(root)})
