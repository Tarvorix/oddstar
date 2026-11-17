"""Registry of all eight search channels."""

from __future__ import annotations

from pathlib import Path

from oddstar.channels.base import Channel
from oddstar.channels.generic_outlier import GenericOutlierChannel
from oddstar.channels.ir_excess import IrExcessChannel
from oddstar.channels.ir_timedomain import IrTimedomainChannel
from oddstar.channels.optical_dips import OpticalDipsChannel
from oddstar.channels.radio_narrowband import RadioNarrowbandChannel
from oddstar.channels.spectral_lines import SpectralLinesChannel
from oddstar.channels.underluminous import UnderluminousChannel
from oddstar.channels.vanishing import VanishingChannel
from oddstar.config import PipelineConfig
from oddstar.errors import ArtifactNotFoundError

CHANNEL_CLASSES = {
    "ir_excess": IrExcessChannel,
    "ir_timedomain": IrTimedomainChannel,
    "optical_dips": OpticalDipsChannel,
    "radio_narrowband": RadioNarrowbandChannel,
    "underluminous": UnderluminousChannel,
    "generic_outlier": GenericOutlierChannel,
    "spectral_lines": SpectralLinesChannel,
    "vanishing": VanishingChannel,
}

ALL_CHANNEL_NAMES = tuple(CHANNEL_CLASSES.keys())


def build_channel(name: str, cfg: PipelineConfig, root: Path, **kwargs: object) -> Channel:
    if name not in CHANNEL_CLASSES:
        raise ArtifactNotFoundError(
            f"Unknown channel '{name}'",
            hint=f"Valid channels: {', '.join(ALL_CHANNEL_NAMES)}",
        )
    return CHANNEL_CLASSES[name](cfg, root, **kwargs)


__all__ = ["ALL_CHANNEL_NAMES", "CHANNEL_CLASSES", "Channel", "build_channel"]
