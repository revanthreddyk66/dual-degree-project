"""
solar_halo_project.src
======================
Core physics modules for the Solar ICS Halo analysis.
"""

from .solar_modulation import (
    modulation_potential_model_I,
    modulation_potential_model_II,
    modulation_potential_model_III,
    apply_force_field,
    modulated_spectrum_at_distance,
)
from .local_interstellar_spectrum import CosmicRayLIS
from .halo_model import SolarHaloModel, build_models_for_phi_values
from .utils import (
    setup_paper_style,
    generate_mock_1year_data,
    MODEL_COLORS, MODEL_LINESTYLES,
    DATA_COLOR, DATA_MARKER, DATA_MS, DATA_CAPSIZE,
)

__all__ = [
    'modulation_potential_model_I',
    'modulation_potential_model_II',
    'modulation_potential_model_III',
    'apply_force_field',
    'modulated_spectrum_at_distance',
    'CosmicRayLIS',
    'SolarHaloModel',
    'build_models_for_phi_values',
    'setup_paper_style',
    'generate_mock_1year_data',
    'MODEL_COLORS',
    'MODEL_LINESTYLES',
    'DATA_COLOR',
    'DATA_MARKER',
    'DATA_MS',
    'DATA_CAPSIZE',
]
