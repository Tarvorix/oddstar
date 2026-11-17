"""Unified Channel protocol definition (Section 4)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol, runtime_checkable

from oddstar.config import PipelineConfig
from oddstar.schemas import ChannelScore, EvidenceRef, Measurement, SampleSpec


@runtime_checkable
class Channel(Protocol):
    name: str

    def sample(self, cfg: PipelineConfig, **kwargs: Any) -> SampleSpec: ...

    def measure(self, star: dict[str, Any], data: dict[str, Any]) -> Measurement: ...

    def score(self, m: Measurement) -> ChannelScore: ...

    def evidence(self, star: dict[str, Any], m: Measurement, out_dir: Path) -> list[EvidenceRef]: ...
