"""Vetting checks and orchestration."""

from oddstar.vet.engine import vet_target
from oddstar.vet.ir_checks import run_ir_vetting
from oddstar.vet.lc_checks import run_lc_vetting
from oddstar.vet.radio_checks import run_radio_vetting
from oddstar.vet.spectral_checks import run_spectral_vetting
from oddstar.vet.underluminous_checks import run_underluminous_vetting
from oddstar.vet.vanishing_checks import run_vanishing_vetting

__all__ = [
    "run_ir_vetting",
    "run_lc_vetting",
    "run_radio_vetting",
    "run_spectral_vetting",
    "run_underluminous_vetting",
    "run_vanishing_vetting",
    "vet_target",
]
