"""
solar_halo_project.src
======================
Core physics and plotting utilities for the Solar ICS-halo analysis.

This package assumes you already possess *real, pre-filtered* Fermi-LAT data
(e.g. a counts cube, exposure cube, spectra, radial profiles, …).  The code
itself does **not** create synthetic or “mock” data anywhere.
"""

# ─────────────── Solar-modulation utilities ────────────────
from .solar_modulation import (
    modulation_potential_model_I,
    modulation_potential_model_II,
    modulation_potential_model_III,
    apply_force_field,
    modulated_spectrum_at_distance,
)

# ─────────────── Local interstellar spectra ────────────────
from .local_interstellar_spectrum import CosmicRayLIS

# ─────────────── Inverse-Compton halo model ────────────────
from .halo_model import SolarHaloModel, build_models_for_phi_values

# ─────────────── Plot styling helpers  ────────────────
from .utils import (
    setup_paper_style,
    MODEL_COLORS,
    MODEL_LINESTYLES,
    DATA_COLOR,
    DATA_MARKER,
    DATA_MS,
    DATA_CAPSIZE,
)

# (optional) IO helpers for *real* analysis products
from .io import load_results, save_results   # assumes you use results/*.pkl

__all__ = [
    # modulation
    "modulation_potential_model_I",
    "modulation_potential_model_II",
    "modulation_potential_model_III",
    "apply_force_field",
    "modulated_spectrum_at_distance",
    # LIS
    "CosmicRayLIS",
    # halo model
    "SolarHaloModel",
    "build_models_for_phi_values",
    # plotting / style
    "setup_paper_style",
    "MODEL_COLORS",
    "MODEL_LINESTYLES",
    "DATA_COLOR",
    "DATA_MARKER",
    "DATA_MS",
    "DATA_CAPSIZE",
    # I/O helpers for real data (no mock!)
    "load_results",
    "save_results",
]
