"""Infrared excess search channel."""

from oddstar.channels.ir_excess.channel import (
    IrExcessChannel,
    build_empirical_locus,
    fit_excess_blackbody_and_shell,
    fit_photosphere_model,
    predict_from_empirical_locus,
)

__all__ = [
    "IrExcessChannel",
    "build_empirical_locus",
    "fit_excess_blackbody_and_shell",
    "fit_photosphere_model",
    "predict_from_empirical_locus",
]
