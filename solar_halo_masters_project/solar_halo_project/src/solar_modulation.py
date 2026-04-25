"""
solar_modulation.py
===================
Solar modulation models for cosmic-ray electrons and positrons.

Implements three force-field models from Linden et al. (2026):
  - Model I  : Radially-dependent potential (Fujii & McDonald 2005 parameterization)
  - Model II : No additional modulation inside 1 AU
  - Model III: Energy-dependent power-law correction

Physical basis:
    The force-field approximation (Gleeson & Axford 1968) is a steady-state
    diffusion-convection equation in the solar system frame. Under spherical
    symmetry, the modulated spectrum J(r, E) at heliocentric distance r is
    related to the Local Interstellar Spectrum J(∞, E) by:

        J(r, E) = J(∞, E) * [E² - E₀²] / [(E + eΦ(r))² - E₀²]

    where E₀ = 0.511 MeV is the electron rest mass energy and Φ(r) is the
    solar modulation potential in units of MV (= MeV/e).

References:
    - Gleeson & Axford (1968): ApJ 154, 1011
    - Fujii & McDonald (2005): Advances in Space Research 35, 611
    - Linden et al. (2026): arXiv:2505.04625
"""

import numpy as np
from typing import Union


# Physical constants
M_ELECTRON_MEV = 0.511          # Electron rest mass energy [MeV]
AU_TO_PC = 4.848e-6             # Astronomical unit to parsec
R_HELIOSPHERE_AU = 100.0        # Heliospheric boundary [AU]


def modulation_potential_model_I(r_AU: Union[float, np.ndarray],
                                  Phi0_MV: float) -> Union[float, np.ndarray]:
    """
    Model I: Radially-dependent solar modulation potential.

    Parameterization from Eq. (8) of Moskalenko et al. (2006), fitting
    data from IMP-8, Voyager 2, and Pioneer 10 during Solar Cycle 21 minimum.

    Formula:
        Φ₁(r) = Φ₀ * [r^(-0.1) - r_b^(-0.1)] / [(1 AU)^(-0.1) - r_b^(-0.1)]

    Parameters
    ----------
    r_AU : float or array
        Heliocentric distance in AU.
    Phi0_MV : float
        Modulation potential at 1 AU [MV = MeV/e].

    Returns
    -------
    float or array
        Modulation potential Φ₁(r) [MV].
    """
    r_b = R_HELIOSPHERE_AU
    r_1au = 1.0
    numerator = r_AU**(-0.1) - r_b**(-0.1)
    denominator = r_1au**(-0.1) - r_b**(-0.1)
    return Phi0_MV * numerator / denominator


def modulation_potential_model_II(r_AU: Union[float, np.ndarray],
                                   Phi0_MV: float) -> Union[float, np.ndarray]:
    """
    Model II: No modulation inside 1 AU, Model I outside.

    Motivated by Parker Solar Probe data (Li et al. 2022) suggesting
    weaker-than-expected modulation within 1 AU.

    Formula:
        Φ₂(r) = Φ₀              for r < 1 AU
        Φ₂(r) = Φ₁(r)           for r ≥ 1 AU

    Parameters
    ----------
    r_AU : float or array
        Heliocentric distance in AU.
    Phi0_MV : float
        Modulation potential at 1 AU [MV].

    Returns
    -------
    float or array
        Modulation potential Φ₂(r) [MV].
    """
    r_AU = np.atleast_1d(np.asarray(r_AU, dtype=float))
    Phi = np.where(r_AU < 1.0,
                   Phi0_MV,
                   modulation_potential_model_I(r_AU, Phi0_MV))
    return Phi.squeeze() if Phi.size == 1 else Phi


def modulation_potential_model_III(r_AU: Union[float, np.ndarray],
                                    Ee_GeV: float,
                                    Phi0_MV: float,
                                    alpha: float) -> Union[float, np.ndarray]:
    """
    Model III: Energy-dependent power-law correction to Model I.

    Adds a phenomenological energy dependence to probe energy-dependent
    modulation effects in the ICS halo.

    Formula:
        Φ₃(r, E) = Φ₁(r) * (E / 10 GeV)^(-α)    for E < 10 GeV
        Φ₃(r, E) = Φ₁(r)                          for E ≥ 10 GeV

    where α ≥ 0 increases modulation at lower energies.

    Parameters
    ----------
    r_AU : float or array
        Heliocentric distance in AU.
    Ee_GeV : float
        Cosmic-ray electron/positron energy [GeV].
    Phi0_MV : float
        Modulation potential at 1 AU [MV].
    alpha : float
        Power-law index (≥ 0). Best-fit: α = 0.0 ± 0.17 (Linden+2026).

    Returns
    -------
    float or array
        Modulation potential Φ₃(r, E) [MV].
    """
    Phi1 = modulation_potential_model_I(r_AU, Phi0_MV)
    E_ref_GeV = 10.0
    if Ee_GeV < E_ref_GeV:
        return Phi1 * (Ee_GeV / E_ref_GeV)**(-alpha)
    else:
        return Phi1


def apply_force_field(J_LIS: Union[float, np.ndarray],
                      Ee_GeV: Union[float, np.ndarray],
                      Phi_MV: float) -> Union[float, np.ndarray]:
    """
    Apply the force-field approximation to modulate the LIS spectrum.

    Gleeson & Axford (1968) force-field solution:

        J(r, E) = J(∞, E) * [E² - E₀²] / [(E + eΦ)² - E₀²]

    where energies are in GeV and Φ is in GV (= GeV/e).

    Parameters
    ----------
    J_LIS : float or array
        Local interstellar spectrum [cm⁻² s⁻¹ sr⁻¹ GeV⁻¹].
    Ee_GeV : float or array
        Cosmic-ray electron/positron kinetic energy [GeV].
    Phi_MV : float
        Modulation potential [MV = MeV/e].

    Returns
    -------
    float or array
        Modulated spectrum J(r, E) [cm⁻² s⁻¹ sr⁻¹ GeV⁻¹].
    """
    E0_GeV = M_ELECTRON_MEV * 1e-3          # 0.511 MeV → GeV
    Phi_GV = Phi_MV * 1e-3                  # MV → GV

    Ee_GeV = np.asarray(Ee_GeV, dtype=float)

    numerator   = Ee_GeV**2 - E0_GeV**2
    denominator = (Ee_GeV + Phi_GV)**2 - E0_GeV**2

    # Avoid division by zero / negative values at very low energies
    denominator = np.where(denominator <= 0, np.inf, denominator)
    numerator   = np.where(numerator   <= 0, 0.0,    numerator)

    return J_LIS * numerator / denominator


def modulated_spectrum_at_distance(J_LIS: np.ndarray,
                                    Ee_GeV: np.ndarray,
                                    r_AU: float,
                                    Phi0_MV: float,
                                    model: str = 'I',
                                    alpha: float = 0.0) -> np.ndarray:
    """
    Compute the modulated e+/e- spectrum at heliocentric distance r.

    Parameters
    ----------
    J_LIS : array
        Local interstellar spectrum at energies Ee_GeV.
    Ee_GeV : array
        Electron/positron energies [GeV].
    r_AU : float
        Heliocentric distance [AU].
    Phi0_MV : float
        Modulation potential at 1 AU [MV].
    model : str
        Which modulation model: 'I', 'II', or 'III'.
    alpha : float
        Energy-dependence index for Model III.

    Returns
    -------
    array
        Modulated spectrum [same units as J_LIS].
    """
    if model == 'I':
        Phi = modulation_potential_model_I(r_AU, Phi0_MV)
        return apply_force_field(J_LIS, Ee_GeV, Phi)
    elif model == 'II':
        Phi = modulation_potential_model_II(r_AU, Phi0_MV)
        return apply_force_field(J_LIS, Ee_GeV, Phi)
    elif model == 'III':
        result = np.zeros_like(J_LIS)
        for i, (J, E) in enumerate(zip(J_LIS, Ee_GeV)):
            Phi = modulation_potential_model_III(r_AU, E, Phi0_MV, alpha)
            result[i] = apply_force_field(J, E, Phi)
        return result
    else:
        raise ValueError(f"Unknown model '{model}'. Choose 'I', 'II', or 'III'.")


if __name__ == "__main__":
    # Quick sanity check
    import matplotlib.pyplot as plt

    r_arr = np.linspace(0.1, 10, 300)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Plot radial dependence of modulation potential
    for Phi0 in [250, 500, 750, 1000]:
        phi_I  = modulation_potential_model_I(r_arr, Phi0)
        phi_II = modulation_potential_model_II(r_arr, Phi0)
        axes[0].plot(r_arr, phi_I,  label=f'Model I,  Φ₀={Phi0} MV')
        axes[0].plot(r_arr, phi_II, '--', label=f'Model II, Φ₀={Phi0} MV')

    axes[0].set_xlabel('Heliocentric Distance [AU]')
    axes[0].set_ylabel('Φ(r) [MV]')
    axes[0].set_title('Solar Modulation Potential vs Distance')
    axes[0].legend(fontsize=7)
    axes[0].set_xlim(0, 5)
    axes[0].axvline(1.0, ls=':', color='gray', label='1 AU')

    # Plot force-field effect on a toy spectrum
    Ee = np.logspace(-2, 2, 200)           # 10 MeV – 100 GeV
    J_LIS_toy = Ee**(-3.3) * 1e4           # toy power-law LIS

    for Phi0 in [0, 250, 500, 1000]:
        J_mod = apply_force_field(J_LIS_toy, Ee, Phi0)
        axes[1].loglog(Ee, Ee**2 * J_mod, label=f'Φ={Phi0} MV')

    axes[1].set_xlabel('Energy [GeV]')
    axes[1].set_ylabel('E² × J(E)')
    axes[1].set_title('Force-Field Modulation Effect')
    axes[1].legend()

    plt.tight_layout()
    plt.savefig('/home/claude/solar_halo_project/tests/test_modulation.png', dpi=100)
    print("✓ solar_modulation.py sanity check passed.")
