"""Pydantic v2 schemas for all oddstar pipeline artifacts and CLI contracts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Literal

from pydantic import BaseModel, Field


class CLIResponse(BaseModel):
    ok: bool = True
    summary: str
    artifacts: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    next: str = ""
    data: dict[str, Any] = Field(default_factory=dict)


class CLIError(BaseModel):
    ok: bool = False
    code: str
    message: str
    hint: str = ""
    retryable: bool = False


class ManifestRecord(BaseModel):
    manifest_id: str
    archive: str
    release: str
    request_signature: str
    query_text: str
    utc_time: str
    row_count: int
    byte_size: int
    sha256: str
    blob_path: str
    evicted: bool = False
    offline_fixture: bool = False


class SampleCutStep(BaseModel):
    cut_name: str
    description: str
    count_before: int
    count_after: int


class SampleSpec(BaseModel):
    channel: str
    gaia_release: str
    config_hash: str
    seed: int
    cuts: list[SampleCutStep] = Field(default_factory=list)
    total_count: int
    healpix_cells: list[int] = Field(default_factory=list)
    source_ids: list[str] = Field(default_factory=list)
    table_path: str = ""


class Measurement(BaseModel):
    source_id: str
    channel: str
    values: dict[str, float] = Field(default_factory=dict)
    uncertainties: dict[str, float] = Field(default_factory=dict)
    epochs: list[str] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)
    extra: dict[str, Any] = Field(default_factory=dict)


class EvidenceRef(BaseModel):
    kind: str
    path: str
    description: str
    sha256: str = ""


class ChannelScore(BaseModel):
    source_id: str
    channel: str
    raw_statistic: float
    false_alarm_prob: float
    flagged: bool
    rank_score: float
    summary: str
    feature_contributions: dict[str, float] = Field(default_factory=dict)
    evidence: list[EvidenceRef] = Field(default_factory=list)


class VettingResult(BaseModel):
    check: str
    status: Literal["pass", "fail", "not_evaluated"]
    evidence: str | None = None
    note: str = ""
    metrics: dict[str, Any] = Field(default_factory=dict)


class HypothesisRow(BaseModel):
    name: Literal[
        "instrument_artifact",
        "catalog_mismatch",
        "blend_neighbouring_star",
        "background_galaxy",
        "galactic_cirrus",
        "solar_system_object",
        "young_star_disk",
        "collision_debris",
        "stellar_companion",
        "evolved_star_mass_loss",
        "known_variable_class",
        "unexplained",
    ]
    status: Literal["disfavoured", "consistent", "untested"]
    evidence: list[str] = Field(default_factory=list)
    note: str = ""


class ValidationCertificate(BaseModel):
    channel: str
    gaia_release: str
    split: Literal["dev", "holdout", "all"]
    config_hash: str
    git_commit: str
    timestamp_utc: str
    passed: bool
    blind_recovery_rate: float
    blind_total: int
    blind_recovered: int
    impostor_pass: bool
    impostor_max_tier: str
    impostors_evaluated: int
    injection_completeness_design_point: float
    completeness_grid: list[dict[str, Any]] = Field(default_factory=list)
    false_alarm_rate: float
    false_alarm_ci_low: float
    false_alarm_ci_high: float
    expected_false_positives: float
    calibration_curve: list[dict[str, float]] = Field(default_factory=list)
    details: dict[str, Any] = Field(default_factory=dict)
    # SHA-256 over the source files that determine the channel's validated behaviour (oddstar.validate.fingerprint);
    # None on certificates issued before code fingerprinting, which therefore cannot authorise a search.
    code_fingerprint: str | None = None


class Coords(BaseModel):
    ra_deg: float
    dec_deg: float
    epoch: str = "J2016.0"


class Provenance(BaseModel):
    run_id: str
    git_commit: str
    config_hash: str
    certificates: list[str] = Field(default_factory=list)
    manifests: list[str] = Field(default_factory=list)


class StarMetadata(BaseModel):
    distance_pc: float | None = None
    teff_k: float | None = None
    luminosity_lsun: float | None = None
    ruwe: float | None = None
    simbad_type: str | None = None
    known_in_literature: bool = False
    gaia_release: str = "dr3"
    ebv_star: float | None = None
    ebv_star_source: str | None = None


class ChannelDossierEntry(BaseModel):
    name: str
    measurements: dict[str, float | None] = Field(default_factory=dict)
    raw_statistic: float | None = None
    false_alarm_prob: float | None = None
    flagged: bool = False
    epochs: list[str] = Field(default_factory=list)
    caveats: list[str] = Field(default_factory=list)
    evidence: list[str] = Field(default_factory=list)
    labels: dict[str, Any] = Field(default_factory=dict)


class DossierReviewSummary(BaseModel):
    round: int
    record: str
    verdict: str
    reviewer_model: str = ""


class Dossier(BaseModel):
    source_id: str
    coords: Coords
    tier: Literal["A", "B", "C", "Explained", "Rejected"]
    provenance: Provenance
    star: StarMetadata
    channels: list[ChannelDossierEntry] = Field(default_factory=list)
    vetting: list[VettingResult] = Field(default_factory=list)
    hypotheses: list[HypothesisRow] = Field(default_factory=list)
    cross_channel_patterns: list[str] = Field(default_factory=list)
    reviews: list[DossierReviewSummary] = Field(default_factory=list)
    what_would_settle_it: list[str] = Field(default_factory=list)
    lead_lines: list[str] = Field(default_factory=list)
    reproduce: str
    extinction_cases: dict[str, Any] | None = None


class ReviewFinding(BaseModel):
    finding_id: str
    severity: Literal["material", "minor", "info"]
    summary: str
    evidence: str
    reproduce_cmd: str
    status: Literal["open", "closed"] = "open"
    linked_test: str | None = None
    lesson_id: str | None = None
    unautomated_reason: str | None = None


class ReviewRecord(BaseModel):
    review_id: str
    round: int
    source_id: str
    reviewer_model: str
    timestamp_utc: str
    measurements: dict[str, float] = Field(default_factory=dict)
    uncertainties: dict[str, float] = Field(default_factory=dict)
    queries_run: list[str] = Field(default_factory=list)
    findings: list[ReviewFinding] = Field(default_factory=list)
    verdict: Literal["reproduced", "not_reproduced", "explained", "unresolved"]
    notes: str = ""


class ReviewPacketTarget(BaseModel):
    packet_target_id: str
    source_id: str
    ra_deg: float
    dec_deg: float
    claim: str | None = None


class ReviewPacket(BaseModel):
    packet_id: str
    round: int
    created_at_utc: str
    prompt_template_path: str
    prompt_text: str
    targets: list[ReviewPacketTarget] = Field(default_factory=list)


class LessonRecord(BaseModel):
    lesson_id: str
    finding_id: str
    source_id: str
    title: str
    what_happened: str
    generic_lesson: str
    automated_check: str
    linked_test: str
    added_to_impostors: bool = False
    skill_updated: str | None = None
    commit_sha: str = ""


class NullReport(BaseModel):
    channel: str
    run_id: str
    gaia_release: str
    config_hash: str
    git_commit: str
    timestamp_utc: str
    sample_cuts: list[SampleCutStep] = Field(default_factory=list)
    sky_coverage_healpix: list[int] = Field(default_factory=list)
    distance_distribution: dict[str, int] = Field(default_factory=dict)
    spectral_type_distribution: dict[str, int] = Field(default_factory=dict)
    completeness_grid: list[dict[str, Any]] = Field(default_factory=list)
    contour_50: list[dict[str, Any]] = Field(default_factory=list)
    contour_90: list[dict[str, Any]] = Field(default_factory=list)
    n_searched: int
    n_flagged: int
    expected_false_alarms: float
    tier_breakdown: dict[str, int] = Field(default_factory=dict)
    upper_limits: list[dict[str, Any]] = Field(default_factory=list)
    manifests: list[str] = Field(default_factory=list)
    reproduce: str


class CandidateFusionRow(BaseModel):
    source_id: str
    tier: Literal["A", "B", "C", "Explained", "Rejected"]
    rank_score: float
    combined_fap: float
    flagged_channels: list[str] = Field(default_factory=list)
    patterns: list[str] = Field(default_factory=list)
    open_holes: list[str] = Field(default_factory=list)
    explained_reason: str | None = None


class FusionResult(BaseModel):
    run_id: str
    timestamp_utc: str
    config_hash: str
    git_commit: str
    n_stars_total: int
    n_trials_total: int
    expected_false_positives: float
    consistent_with_noise: bool
    cross_channel_patterns: list[dict[str, Any]] = Field(default_factory=list)
    candidates: list[CandidateFusionRow] = Field(default_factory=list)


class RegistryEntry(BaseModel):
    entry_id: str
    timestamp_utc: str
    source_id: str
    event_type: Literal[
        "candidate_registered",
        "tier_updated",
        "review_ingested",
        "finding_opened",
        "finding_closed",
        "lesson_added",
        "explained_anomaly",
    ]
    tier: Literal["A", "B", "C", "Explained", "Rejected"] | None = None
    status: str = ""
    provenance: dict[str, Any] = Field(default_factory=dict)
    details: dict[str, Any] = Field(default_factory=dict)


class RunState(BaseModel):
    run_id: str
    stage: str
    channel: str
    gaia_release: str
    config_hash: str
    git_commit: str
    seed: int
    status: Literal["completed", "failed", "dry_run"]
    timestamp_utc: str
    artifacts: list[str] = Field(default_factory=list)
    summary: dict[str, Any] = Field(default_factory=dict)


SCHEMA_MODELS: dict[str, type[BaseModel]] = {
    "cli_response": CLIResponse,
    "cli_error": CLIError,
    "manifest_record": ManifestRecord,
    "sample_spec": SampleSpec,
    "measurement": Measurement,
    "channel_score": ChannelScore,
    "vetting_result": VettingResult,
    "hypothesis_row": HypothesisRow,
    "validation_certificate": ValidationCertificate,
    "dossier": Dossier,
    "review_record": ReviewRecord,
    "review_packet": ReviewPacket,
    "lesson_record": LessonRecord,
    "null_report": NullReport,
    "fusion_result": FusionResult,
    "registry_entry": RegistryEntry,
    "run_state": RunState,
}


def export_all_schemas(root: Path) -> list[Path]:
    """Export JSON Schemas for every artifact model to schemas/ and review/."""
    schema_dir = root / "schemas"
    schema_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, model_cls in SCHEMA_MODELS.items():
        path = schema_dir / f"{name}.schema.json"
        schema_dict = model_cls.model_json_schema()
        path.write_text(json.dumps(schema_dict, indent=2, sort_keys=True) + "\n")
        written.append(path)
    review_dir = root / "review"
    review_dir.mkdir(parents=True, exist_ok=True)
    rr_path = review_dir / "review_record.schema.json"
    rr_path.write_text(json.dumps(ReviewRecord.model_json_schema(), indent=2, sort_keys=True) + "\n")
    written.append(rr_path)
    return written
