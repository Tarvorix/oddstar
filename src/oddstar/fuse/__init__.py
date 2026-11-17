"""Cross-channel fusion and ranking."""

from oddstar.fuse.engine import (
    CHANNEL_INSTRUMENT_GROUP,
    assign_tier,
    build_hypothesis_table,
    calibrate_fap,
    detect_cross_channel_patterns,
    fuse_results,
)

__all__ = [
    "CHANNEL_INSTRUMENT_GROUP",
    "assign_tier",
    "build_hypothesis_table",
    "calibrate_fap",
    "detect_cross_channel_patterns",
    "fuse_results",
]
