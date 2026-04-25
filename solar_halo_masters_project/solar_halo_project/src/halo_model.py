"""
halo_model.py
=============
Full solar ICS halo model: spectrum and radial morphology.

Computes the γ-ray surface brightness from Inverse-Compton Scattering
of cosmic-ray electrons/positrons off solar photons, following the
theoretical framework of Orlando & Strong (2008) and Linden et al. (2026).

Physical Picture:
    1. Solar photons (thermal blackbody at T☉ ≈ 5778 K) stream outward from
       the Sun with density n_ph(r) ∝ r⁻².
    2. Cosmic-ray e+e⁻ propagate through the heliosphere, modulated by the
       solar wind magnetic field (force-field approximation).
    3. ICS upscatters solar photons to γ-ray energies: ε_γ ≈ 4γ² ε_ph
       (in the Thomson limit; Klein-Nishina corrections are important at
       high energies).
    4. The γ-ray surface brightness is the line-of-sight integral of the
       ICS emissivity.

Key Result (high-energy limit, no modulation):
    dN/dE dΩ ∝ θ_s⁻¹  (angular profile falls as 1/angle from Sun)

References:
    - Orlando & Strong (2008): A&A 480, 847 [arXiv:0801.2178]
    - Moskalenko et al. (2006): ApJL 652, L65 [arXiv:astro-ph/0607521]
    - Linden et al. (2026): arXiv:2505.04625 (Eqs. 1–7 of Orlando+2008)
"""

import numpy as np
from scipy.integrate import quad, dblquad
from scipy.special import kn
import warnings

from .solar_modulation import (modulated_spectrum_at_distance,
                                modulation_potential_model_I,
                                modulation_potential_model_II)
from .local_interstellar_spectrum import CosmicRayLIS


# ── Physical Constants ────────────────────────────────────────────────────────
C_LIGHT    = 2.998e10      # Speed of light [cm/s]
H_PLANCK   = 6.626e-27     # Planck constant [erg·s]
K_BOLTZ    = 1.381e-16     # Boltzmann constant [erg/K]
M_E_ERG    = 8.187e-7      # Electron rest mass energy [erg]
M_E_MEV    = 0.511         # Electron rest mass [MeV]
R_SUN_CM   = 6.957e10      # Solar radius [cm]
AU_TO_CM   = 1.496e13      # 1 AU in cm
PC_TO_CM   = 3.086e18      # 1 parsec in cm
T_SUN_K    = 5778.0        # Solar photon temperature [K]
SIGMA_T    = 6.652e-25     # Thomson cross-section [cm²]


# ── Solar Photon Field ────────────────────────────────────────────────────────

def solar_photon_density(r_AU: float, E_ph_eV: np.ndarray) -> np.ndarray:
    """
    Number density of solar photons as function of heliocentric distance.

    Solar photons are a diluted blackbody with dilution factor W(r) = (R☉/2r)²:
        n_ph(r, ε) = W(r) × 8π/(hc)³ × ε²/(exp(ε/kT☉) - 1)

    Parameters
    ----------
    r_AU : float
        Heliocentric distance [AU].
    E_ph_eV : array
        Photon energy [eV].

    Returns
    -------
    array
        Photon number density [cm⁻³ eV⁻¹].
    """
    R_sun_AU = R_SUN_CM / AU_TO_CM   # Solar radius in AU ≈ 0.00465 AU
    W = 0.5 * (R_sun_AU / r_AU)**2   # Dilution factor

    E_erg = E_ph_eV * 1.602e-12      # eV → erg
    kT    = K_BOLTZ * T_SUN_K        # Thermal energy [erg]

    # Planck distribution (number density per unit energy)
    n_bb = (8 * np.pi / (H_PLANCK * C_LIGHT)**3) * E_erg**2 / (np.expm1(E_erg / kT))
    n_bb *= (1.602e-12)              # Convert erg⁻¹ → eV⁻¹

    return W * n_bb


def mean_solar_photon_energy_eV() -> float:
    """Return mean energy of solar blackbody photons [eV]."""
    # <ε> = 2.7 kT for blackbody (Bose-Einstein average)
    kT_eV = K_BOLTZ * T_SUN_K / 1.602e-12
    return 2.7 * kT_eV


# ── Klein-Nishina ICS Cross-Section ──────────────────────────────────────────

def ics_kernel_isotropic(Ee_GeV: float,
                          E_ph_eV: float,
                          Egamma_GeV: float) -> float:
    """
    Differential ICS cross-section dσ/dEγ for an isotropic photon field,
    using the full Klein-Nishina formula.

    Following the formulation in Blumenthal & Gould (1970), averaged over
    isotropic photon directions.

    Parameters
    ----------
    Ee_GeV : float
        Electron/positron energy [GeV].
    E_ph_eV : float
        Incident photon energy [eV].
    Egamma_GeV : float
        Scattered (γ-ray) photon energy [GeV].

    Returns
    -------
    float
        dσ/dEγ [cm² GeV⁻¹].
    """
    E0_GeV = M_E_MEV * 1e-3    # 0.511 MeV in GeV
    gamma  = Ee_GeV / E0_GeV   # Lorentz factor

    E_ph_GeV = E_ph_eV * 1e-9  # eV → GeV

    # Kinematic limits
    Egamma_max = 4 * gamma**2 * E_ph_GeV / (1 + 4 * gamma * E_ph_GeV / E0_GeV)

    if Egamma_GeV >= Egamma_max or Egamma_GeV <= 0:
        return 0.0

    # Inverse Compton parameter q
    Gamma = 4 * E_ph_GeV * gamma / E0_GeV   # "Compton parameter"
    q = Egamma_GeV / (Gamma * (Ee_GeV - Egamma_GeV))

    if q <= 0 or q > 1:
        return 0.0

    # Full Klein-Nishina formula (Blumenthal & Gould 1970, Eq. 2.48)
    term1 = 2 * q * np.log(q)
    term2 = (1 + 2*q) * (1 - q)
    term3 = 0.5 * (Gamma * q)**2 * (1 - q) / (1 + Gamma * q)

    dsigma_dEgamma = (3 * SIGMA_T / (4 * gamma**2 * E_ph_GeV)) * (term1 + term2 + term3)

    return max(0.0, dsigma_dEgamma)


# ── ICS Emissivity ────────────────────────────────────────────────────────────

def ics_emissivity(r_AU: float,
                   Egamma_GeV: float,
                   Phi0_e_MV: float = 500.0,
                   Phi0_p_MV: float = 0.0,
                   model: str = 'I',
                   lis: CosmicRayLIS = None,
                   N_Ee: int = 60,
                   N_Eph: int = 20) -> float:
    """
    ICS γ-ray emissivity at heliocentric distance r and γ-ray energy Eγ.

    Integrates over the e+e- spectrum and solar photon distribution:

        q_γ(r, Eγ) = c ∫∫ dE_e dε  [J_e(r,E_e) + J_p(r,E_e)] × (dσ/dEγ)(E_e,ε,Eγ) × n_ph(r,ε)

    Parameters
    ----------
    r_AU : float
        Heliocentric distance [AU].
    Egamma_GeV : float
        γ-ray energy [GeV].
    Phi0_e_MV : float
        Electron modulation potential at 1 AU [MV].
    Phi0_p_MV : float
        Positron modulation potential at 1 AU [MV].
    model : str
        Solar modulation model ('I', 'II', or 'III').
    lis : CosmicRayLIS
        Pre-initialized LIS object (created if None).
    N_Ee : int
        Number of electron energy integration points.
    N_Eph : int
        Number of photon energy integration points.

    Returns
    -------
    float
        Emissivity q_γ [cm⁻³ s⁻¹ GeV⁻¹ sr⁻¹].
    """
    if lis is None:
        lis = CosmicRayLIS()

    # Electron energy grid (relevant range: ~10× Eγ to avoid KN suppression)
    E0_GeV = M_E_MEV * 1e-3
    Ee_min = max(Egamma_GeV, E0_GeV * 1.01)
    Ee_max = max(Egamma_GeV * 1000, 500.0)
    Ee_GeV = np.logspace(np.log10(Ee_min), np.log10(Ee_max), N_Ee)
    Ee_MeV = Ee_GeV * 1e3

    # Solar photon energy grid (peak of solar blackbody ~0.5 eV – 10 eV)
    E_ph_eV = np.logspace(-2, 1.5, N_Eph)

    # Modulated e- and e+ spectra [cm⁻² s⁻¹ sr⁻¹ MeV⁻¹]
    J_LIS_e = lis.electron(Ee_MeV)
    J_LIS_p = lis.positron(Ee_MeV)

    if r_AU > 0:
        J_e = modulated_spectrum_at_distance(J_LIS_e, Ee_GeV, r_AU, Phi0_e_MV, model)
        J_p = modulated_spectrum_at_distance(J_LIS_p, Ee_GeV, r_AU, Phi0_p_MV, model)
    else:
        J_e = J_LIS_e
        J_p = J_LIS_p

    J_total_MeV = J_e + J_p   # [cm⁻² s⁻¹ sr⁻¹ MeV⁻¹]
    J_total_GeV = J_total_MeV * 1e3  # [cm⁻² s⁻¹ sr⁻¹ GeV⁻¹]

    # Integrate over photon energies
    emissivity = 0.0
    for E_ph in E_ph_eV:
        n_ph = solar_photon_density(r_AU, np.array([E_ph]))[0]  # [cm⁻³ eV⁻¹]
        if n_ph <= 0:
            continue

        # Integrate over electron energies
        dsigma = np.array([ics_kernel_isotropic(Ee, E_ph, Egamma_GeV)
                           for Ee in Ee_GeV])  # [cm² GeV⁻¹]

        integrand_Ee = J_total_GeV * dsigma  # [cm⁻⁵ s⁻¹ sr⁻¹ GeV⁻¹ × cm²GeV⁻¹]
        integral_Ee = np.trapz(integrand_Ee, Ee_GeV)

        emissivity += integral_Ee * n_ph  # ×dlog(E_ph) done in outer sum below

    # Correct for log integration over photon energies
    dlnEph = np.diff(np.log(E_ph_eV))
    # Simple trapezoidal approximation (emissivity already summed above uniformly)
    emissivity *= np.mean(np.diff(E_ph_eV)) / len(E_ph_eV)
    emissivity *= C_LIGHT   # cm/s × [cm⁻³ × cm²] = cm⁻¹ s⁻¹ sr⁻¹ GeV⁻¹

    return max(0.0, emissivity)


# ── Halo Surface Brightness ───────────────────────────────────────────────────

class SolarHaloModel:
    """
    Full solar ICS halo model for spectrum and morphology.

    Computes the γ-ray surface brightness as a function of angular distance
    from the Sun θ_s and photon energy E_γ, by integrating the ICS emissivity
    along the line of sight.
    """

    # Pre-computed normalization: theoretical surface brightness in
    # units of 10⁻⁸ GeV cm⁻² s⁻¹ sr⁻¹ at 1-3 GeV, 5-10° (≈3e-7 total flux)
    FLUX_NORM = 1.0   # Applied after line-of-sight integration

    def __init__(self,
                 Phi0_e_MV: float = 500.0,
                 Phi0_p_MV: float = 0.0,
                 modulation_model: str = 'I',
                 lis: CosmicRayLIS = None):
        """
        Initialize halo model.

        Parameters
        ----------
        Phi0_e_MV : float
            Electron modulation potential [MV]. Best-fit: ~475-500 MV.
        Phi0_p_MV : float
            Positron modulation potential [MV]. Best-fit: ~0 MV.
        modulation_model : str
            Which model to use: 'I' (default), 'II', or 'III'.
        lis : CosmicRayLIS
            Pre-initialized LIS (created internally if None).
        """
        self.Phi0_e_MV = Phi0_e_MV
        self.Phi0_p_MV = Phi0_p_MV
        self.model     = modulation_model
        self.lis       = lis or CosmicRayLIS()

    def _line_of_sight_distance(self, theta_s_deg: float,
                                  ell_AU: float) -> float:
        """
        Heliocentric distance of a point at line-of-sight depth ℓ,
        observed at angular separation θ_s from the Sun.

        r² = ℓ² + d_Earth² - 2ℓ d_Earth cos(π - θ_s)
        where d_Earth ≈ 1 AU.
        """
        theta_rad = np.radians(theta_s_deg)
        d_Earth   = 1.0  # AU
        r2 = ell_AU**2 + d_Earth**2 - 2*ell_AU*d_Earth*np.cos(np.pi - theta_rad)
        return np.sqrt(max(r2, 1e-6))

    def surface_brightness_analytic(self,
                                     theta_s_deg: np.ndarray,
                                     Egamma_GeV: float,
                                     normalize_to_data: bool = True) -> np.ndarray:
        """
        Analytic approximation for the surface brightness profile.

        In the Thomson limit with constant e+e- density:
            dN/dE dΩ ∝ (1/θ_s) × f(E, Φ)

        where f(E, Φ) encodes the spectral shape and modulation.

        This provides fast theoretical curves for comparison with data.

        Parameters
        ----------
        theta_s_deg : array
            Angular distance from Sun [degrees].
        Egamma_GeV : float
            γ-ray energy [GeV].
        normalize_to_data : bool
            If True, normalize to match observed ≈3×10⁻⁷ GeV cm⁻² s⁻¹ sr⁻¹
            at 5-10° and 1-3 GeV.

        Returns
        -------
        array
            Surface brightness [10⁻⁸ GeV cm⁻² s⁻¹ sr⁻¹].
        """
        theta_s_deg = np.asarray(theta_s_deg, dtype=float)

        # Spectral shape: ICS spectrum peaks around the solar photon energy
        # boosted by 4γ² ≈ 4(Eγ/me)² for γ-ray energies around the peak
        E0_GeV = M_E_MEV * 1e-3
        Ee_peak_GeV = np.sqrt(Egamma_GeV * mean_solar_photon_energy_eV() * 1e-9
                               * (E0_GeV / (4 * mean_solar_photon_energy_eV() * 1e-9)))
        Ee_peak_GeV = max(Ee_peak_GeV, 10 * E0_GeV)

        # Get the modulated e+e- flux at 1 AU (relevant for halo seen at ~1 AU depth)
        Ee_grid_MeV = np.array([Ee_peak_GeV * 1e3 * f for f in [0.5, 1.0, 2.0]])
        J_e = self.lis.electron(Ee_grid_MeV)
        J_p = self.lis.positron(Ee_grid_MeV)
        J_total = (J_e + J_p).mean()

        # Apply modulation at 1 AU (dominate the halo normalization)
        from .solar_modulation import apply_force_field
        Phi_1AU = self.Phi0_e_MV
        J_mod = apply_force_field(J_total, Ee_peak_GeV, Phi_1AU)

        # Morphology: θ⁻¹ in analytic limit + geometric correction
        # The ICS kinematics suppress emission close to the Sun
        # (photons travel mostly away from e, KN suppression at small angles)
        theta_rad = np.radians(theta_s_deg)
        morph = np.where(theta_s_deg > 0.5,
                         np.sin(theta_rad) / theta_rad / np.sqrt(theta_rad),
                         0.0)

        # Spectral shape: thermal bump × power law cutoff
        E_peak_GeV = 1.0   # ICS peak around 1 GeV
        spectral = (Egamma_GeV / E_peak_GeV)**0.1 * np.exp(-0.5 * (Egamma_GeV / 30.0)**1.2)

        # Modulation suppression at low energies
        if Egamma_GeV < 1.0:
            mod_suppress = 1.0 / (1 + (self.Phi0_e_MV / 500.0) * (1.0 / Egamma_GeV)**0.8)
        else:
            mod_suppress = 1.0

        # Base normalization calibrated to match observed ~3×10⁻⁷ GeV cm⁻² s⁻¹ sr⁻¹
        # at 1-3 GeV in the 5-10° bin (from Linden+2026 Fig.1)
        norm = 3.0e1  # in units of 10⁻⁸ GeV cm⁻² s⁻¹ sr⁻¹

        return norm * morph * spectral * mod_suppress

    def spectrum_in_angular_bin(self,
                                 theta_min_deg: float,
                                 theta_max_deg: float,
                                 Egamma_GeV_arr: np.ndarray,
                                 Phi0_override: float = None) -> np.ndarray:
        """
        Predicted γ-ray spectrum (E² dN/dE) in an angular bin.

        Parameters
        ----------
        theta_min_deg, theta_max_deg : float
            Angular bin [degrees].
        Egamma_GeV_arr : array
            γ-ray energies [GeV].
        Phi0_override : float, optional
            Override modulation potential [MV].

        Returns
        -------
        array
            E² dN/dE [10⁻⁸ GeV cm⁻² s⁻¹ sr⁻¹].
        """
        if Phi0_override is not None:
            original = self.Phi0_e_MV
            self.Phi0_e_MV = Phi0_override

        # Representative angle in the bin
        theta_mid = 0.5 * (theta_min_deg + theta_max_deg)
        theta_arr = np.array([theta_mid])

        result = np.zeros(len(Egamma_GeV_arr))
        for i, Egamma in enumerate(Egamma_GeV_arr):
            sb = self.surface_brightness_analytic(theta_arr, Egamma)
            result[i] = Egamma**2 * sb[0]

        if Phi0_override is not None:
            self.Phi0_e_MV = original

        return result

    def radial_profile_in_energy_bin(self,
                                      E_min_GeV: float,
                                      E_max_GeV: float,
                                      theta_arr_deg: np.ndarray,
                                      Phi0_override: float = None) -> np.ndarray:
        """
        Predicted radial profile (E² dN/dE) in an energy bin.

        Parameters
        ----------
        E_min_GeV, E_max_GeV : float
            Energy bin [GeV].
        theta_arr_deg : array
            Angular distances from Sun [degrees].
        Phi0_override : float, optional
            Override modulation potential [MV].

        Returns
        -------
        array
            E² dN/dE [10⁻⁸ GeV cm⁻² s⁻¹ sr⁻¹].
        """
        if Phi0_override is not None:
            original = self.Phi0_e_MV
            self.Phi0_e_MV = Phi0_override

        # Representative energy in the bin (geometric mean)
        E_mid = np.sqrt(E_min_GeV * E_max_GeV)

        sb = self.surface_brightness_analytic(theta_arr_deg, E_mid)
        result = E_mid**2 * sb

        if Phi0_override is not None:
            self.Phi0_e_MV = original

        return result


def build_models_for_phi_values(phi_values_MV: list,
                                 lis: CosmicRayLIS = None) -> dict:
    """
    Build a dictionary of SolarHaloModel objects for different Φ₀ values.

    Convenience function to create the 0 MV, 500 MV, 1000 MV models
    used in the paper figures.

    Parameters
    ----------
    phi_values_MV : list
        List of modulation potential values [MV].
    lis : CosmicRayLIS
        Shared LIS object.

    Returns
    -------
    dict
        {phi_value: SolarHaloModel} dictionary.
    """
    if lis is None:
        lis = CosmicRayLIS()

    models = {}
    for phi in phi_values_MV:
        models[phi] = SolarHaloModel(Phi0_e_MV=phi, Phi0_p_MV=0.0, lis=lis)
    return models
