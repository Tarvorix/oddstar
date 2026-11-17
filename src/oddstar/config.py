"""Versioned pipeline configuration loading, validation, and deterministic hashing."""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any

import yaml
from pydantic import BaseModel, Field, ValidationError

from oddstar.errors import ConfigError


class StorageConfig(BaseModel):
    max_cache_gb: float = Field(gt=0.0, le=50.0)
    healpix_nside: int = Field(ge=1)
    raw_dir: str = "data/raw"
    tables_dir: str = "data/tables"
    runs_dir: str = "runs"
    dossiers_dir: str = "dossiers"
    registry_dir: str = "registry"
    evidence_dir: str = "registry/evidence"


class GaiaConfig(BaseModel):
    release: str = "dr3"


class SampleConfig(BaseModel):
    seed: int
    field_sample_size: int = Field(ge=100)
    max_distance_pc: float = Field(gt=0.0)
    min_abs_gal_lat_deg: float = Field(ge=0.0)
    min_parallax_snr: float = Field(gt=0.0)
    max_ruwe: float = Field(gt=0.0)
    require_unique_allwise_xmatch: bool = True
    require_clean_2mass: bool = True


class ArchivePolicyConfig(BaseModel):
    enabled: bool = True
    release: str
    rate_limit_rps: float = Field(gt=0.0)
    requires_registration: bool = False


class ValidationPolicyConfig(BaseModel):
    blind_recovery_min: float = Field(ge=0.0, le=1.0)
    impostor_max_tier: str = "C"
    injection_completeness_min: float = Field(ge=0.0, le=1.0)
    false_alarm_confidence_level: float = Field(gt=0.0, lt=1.0)
    field_sample_size: int = Field(ge=100)
    injection_trials_per_cell: int = Field(ge=1)


class FusionPolicyConfig(BaseModel):
    tier_a_min_channels: int = 2
    tier_a_max_fap: float = 1e-4
    tier_b_max_fap: float = 1e-3
    # Since 1.1.0 (moved from oddstar.fuse.engine literals; values unchanged)
    fap_floor: float = Field(gt=0.0, lt=1.0)
    fallback_false_alarm_rate_per_star: float = Field(gt=0.0, lt=1.0)


class EvidenceStoreConfig(BaseModel):
    """Size guards for the committed evidence store (`storage.evidence_dir`); raw data never enters the repository."""

    max_array_elements: int = Field(ge=1)
    max_string_chars: int = Field(ge=16)
    max_json_bytes: int = Field(ge=1024)
    max_image_bytes: int = Field(ge=0)


class LiveIngestionConfig(BaseModel):
    """Live Gaia x AllWISE x 2MASS sample cuts (`oddstar.io.live_survey`) and live search filters (`oddstar.cli`).

    `regions` maps a region preset (halo, bulge, local) to its numeric Gaia cuts; the SQL clause structure of each
    preset stays in `live_survey._region_where_sql`.
    """

    allwise_number_of_neighbours: int
    allwise_number_of_mates: int
    tmass_number_of_neighbours: int
    tmass_number_of_mates: int
    regions: dict[str, dict[str, float]]
    ir_excess_w4_upper_limit_ph_qual: list[str]
    ir_excess_artifact_cc_flag_bands: list[int]


class PipelineConfig(BaseModel):
    version: str
    storage: StorageConfig
    gaia: GaiaConfig
    sample: SampleConfig
    channels: dict[str, dict[str, Any]]
    vetting: dict[str, float]
    validation: ValidationPolicyConfig
    fusion: FusionPolicyConfig
    archives: dict[str, ArchivePolicyConfig]
    evidence: EvidenceStoreConfig
    live_ingestion: LiveIngestionConfig


def repo_root(start: Path | None = None) -> Path:
    """Resolve repository root directory by locating pyproject.toml / configs."""
    cur = (start or Path.cwd()).resolve()
    for candidate in [cur, *cur.parents]:
        if (candidate / "pyproject.toml").exists() and (candidate / "configs").exists():
            return candidate
    return Path(__file__).resolve().parents[2]


def compute_config_hash(raw_dict: dict[str, Any]) -> str:
    """Compute canonical SHA-256 hex digest (16-char prefix and full hash) of config."""
    canonical = json.dumps(raw_dict, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(canonical.encode("utf-8")).hexdigest()


def current_git_commit(root: Path | None = None) -> str:
    """Return current HEAD commit SHA or 'uncommitted-workspace' when repo has no commits yet."""
    base = root or repo_root()
    try:
        res = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(base),
            capture_output=True,
            text=True,
            check=False,
        )
        if res.returncode == 0 and res.stdout.strip():
            return res.stdout.strip()
    except OSError:
        pass
    return "0000000000000000000000000000000000000000"


def load_config(
    config_path: Path | None = None,
    *,
    root: Path | None = None,
    gaia_release_override: str | None = None,
) -> tuple[PipelineConfig, str, str]:
    """Load and validate pipeline config, returning (config, short_hash, full_hash)."""
    base = root or repo_root()
    path = config_path or (base / "configs" / "default.yaml")
    if not path.exists():
        raise ConfigError(
            f"Config file not found: {path}",
            hint="Verify --config path or restore configs/default.yaml",
        )
    try:
        raw = yaml.safe_load(path.read_text(encoding="utf-8"))
        if not isinstance(raw, dict):
            raise ConfigError("Config YAML root must be a mapping")
        if gaia_release_override:
            raw.setdefault("gaia", {})["release"] = gaia_release_override
        cfg = PipelineConfig.model_validate(raw)
    except (yaml.YAMLError, ValidationError) as exc:
        raise ConfigError(f"Invalid pipeline config: {exc}", hint="Fix configs/default.yaml") from exc

    full_hash = compute_config_hash(cfg.model_dump())
    return cfg, full_hash[:12], full_hash
