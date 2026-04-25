"""
utils.py
========
Utility functions for coordinate transforms, plotting helpers, and
data I/O for the solar ICS halo analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.ticker import LogLocator, LogFormatter
import warnings


# ── Coordinate Transforms ─────────────────────────────────────────────────────

def angular_separation_deg(ra1, dec1, ra2, dec2):
    """
    Great-circle angular separation between two sky coordinates.
    All inputs and output in degrees (uses Haversine formula).
    """
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    delta_dec = dec2 - dec1
    delta_ra  = ra2 - ra1
    a = np.sin(delta_dec/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin(delta_ra/2)**2
    return np.degrees(2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a)))


def heliocentric_distance_along_los(theta_s_deg, ell_AU, d_obs_AU=1.0):
    """
    Heliocentric distance of a point along the line of sight.

    Parameters
    ----------
    theta_s_deg : float
        Angular separation of LOS center from Sun [degrees].
    ell_AU : float or array
        Line-of-sight distance from observer [AU].
    d_obs_AU : float
        Distance of observer from Sun [AU]. Default: 1 AU (Earth).

    Returns
    -------
    float or array
        Heliocentric distance r [AU].
    """
    theta_rad = np.radians(theta_s_deg)
    cos_angle = np.cos(np.pi - theta_rad)   # angle measured from Sun direction
    r2 = ell_AU**2 + d_obs_AU**2 - 2 * ell_AU * d_obs_AU * cos_angle
    return np.sqrt(np.maximum(r2, 1e-10))


def los_integration_limits(theta_s_deg, d_obs_AU=1.0, max_dist_AU=50.0):
    """
    Return integration limits for the line-of-sight integral.

    Parameters
    ----------
    theta_s_deg : float
        Angular distance from Sun [degrees].
    d_obs_AU : float
        Observer-Sun distance [AU].
    max_dist_AU : float
        Maximum integration depth [AU].

    Returns
    -------
    (ell_min, ell_max) : tuple
        Integration range [AU].
    """
    # Minimum distance: avoid Sun surface if looking close to Sun
    ell_min = 0.01   # AU — just outside the Sun

    # Maximum: ~10-50 AU is sufficient (halo falls off)
    ell_max = max_dist_AU

    return ell_min, ell_max


# ── Plotting Utilities ────────────────────────────────────────────────────────

# Paper color scheme (matching Linden+2026)
MODEL_COLORS = {
    0:    '#E85D3C',   # Coral/salmon — 0 MV
    500:  '#2CA6A4',   # Teal — 500 MV
    1000: '#2C3E50',   # Dark navy — 1000 MV
}
MODEL_LINESTYLES = {
    'I':  '-',    # Model I: solid
    'II': '--',   # Model II: dashed
}
DATA_COLOR  = 'black'
DATA_MARKER = 'D'
DATA_MS     = 4
DATA_CAPSIZE = 3


def setup_paper_style():
    """Apply matplotlib style matching the paper figures."""
    plt.rcParams.update({
        'font.family':       'serif',
        'font.size':         10,
        'axes.labelsize':    11,
        'axes.titlesize':    11,
        'legend.fontsize':   8,
        'xtick.labelsize':   9,
        'ytick.labelsize':   9,
        'lines.linewidth':   1.8,
        'axes.linewidth':    0.8,
        'xtick.direction':   'in',
        'ytick.direction':   'in',
        'xtick.top':         True,
        'ytick.right':       True,
        'figure.dpi':        150,
        'savefig.dpi':       200,
        'savefig.bbox':      'tight',
    })


def add_model_legend(ax, phi_values=(0, 500, 1000), model_types=('I',)):
    """Add a legend showing the modulation potential models."""
    handles = []
    for phi in phi_values:
        color = MODEL_COLORS.get(phi, 'gray')
        for mtype in model_types:
            ls = MODEL_LINESTYLES.get(mtype, '-')
            label = f'{phi} MV'
            line, = ax.plot([], [], color=color, ls=ls, lw=1.8, label=label)
            handles.append(line)
    if len(model_types) > 1:
        # Add "data" entry
        ax.plot([], [], 'D', color=DATA_COLOR, ms=DATA_MS, label='Data')
    return ax.legend(handles=handles, frameon=True, framealpha=0.8,
                     edgecolor='gray', handlelength=2.0)


def format_energy_axis(ax, which='x'):
    """Format an axis as a log-scale energy axis."""
    axis = ax.xaxis if which == 'x' else ax.yaxis
    ax.set_xscale('log') if which == 'x' else ax.set_yscale('log')


def add_energy_bin_label(ax, E_min, E_max, unit='GeV',
                          loc='upper right', fontsize=9):
    """Add energy range label to a panel."""
    label = f'{E_min} - {E_max} {unit}'
    ax.text(0.97, 0.95, label, transform=ax.transAxes,
            ha='right', va='top', fontsize=fontsize,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))


def add_angle_bin_label(ax, theta_min, theta_max, unit='°',
                         loc='upper right', fontsize=9):
    """Add angular range label to a panel."""
    label = f'{theta_min}{unit}–{theta_max}{unit}'
    ax.text(0.97, 0.95, label, transform=ax.transAxes,
            ha='right', va='top', fontsize=fontsize,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))


# ── Mock Data Generator ───────────────────────────────────────────────────────

def generate_mock_1year_data(seed: int = 42) -> dict:
    """
    Generate realistic mock Fermi-LAT data for 1-year analysis.

    This uses the theoretical model predictions with added noise to simulate
    what a 1-year dataset would look like. The noise is scaled appropriately:
    1/15 of the photon statistics in the 15-year paper analysis.

    Returns a dictionary with:
        - 'spectrum_theta_bins': {bin_label: {'E_GeV', 'flux', 'flux_err'}}
        - 'radial_E_bins':       {bin_label: {'theta_deg', 'flux', 'flux_err'}}
        - 'helioprojective_map': 2D counts map at 1-3 GeV

    Parameters
    ----------
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict
        Mock data dictionary.
    """
    rng = np.random.default_rng(seed)

    # ── Spectral bins (Fig. 3 style) ──
    theta_bins = [
        ('<2.5',   0.0,   2.5),
        ('2.5-5',  2.5,   5.0),
        ('5-10',   5.0,  10.0),
        ('10-15', 10.0,  15.0),
        ('15-20', 15.0,  20.0),
        ('20-25', 20.0,  25.0),
        ('25-35', 25.0,  35.0),
        ('35-45', 35.0,  45.0),
    ]

    E_GeV = np.logspace(-1.5, 2, 16)    # 32 MeV – 100 GeV

    # Theoretical model (500 MV, Model I)
    from .halo_model import SolarHaloModel
    model_500 = SolarHaloModel(Phi0_e_MV=500.0)

    spectrum_data = {}
    for label, th_min, th_max in theta_bins:
        flux_theory = model_500.spectrum_in_angular_bin(th_min, th_max, E_GeV)

        # Scale: 1-year is sqrt(1/15) statistical uncertainty relative to 15-yr
        # Add more noise at low/high energies and far angles
        theta_mid = 0.5 * (th_min + th_max)
        noise_factor = 1.0 / np.sqrt(15)  # 1-yr vs 15-yr

        # Fractional error: larger at high energy (few photons) and large angles
        frac_err = noise_factor * (0.08 + 0.3 * (E_GeV / 50)**0.5 +
                                   0.1 * (theta_mid / 20)**0.5)
        frac_err = np.clip(frac_err, 0.05, 2.0)

        flux_err = flux_theory * frac_err

        # Add correlated noise (simulate actual fluctuations)
        noise = rng.normal(0, 1, size=len(E_GeV))
        flux_mock = np.maximum(flux_theory * (1 + frac_err * noise), 0.3)

        # Zero out unphysical points (very high energy with huge errors)
        upper_limit_mask = (frac_err > 0.8) & (rng.random(len(E_GeV)) > 0.4)

        spectrum_data[label] = {
            'E_GeV':      E_GeV,
            'flux':       flux_mock,
            'flux_err':   flux_err,
            'upper_lim':  upper_limit_mask,
            'theta_min':  th_min,
            'theta_max':  th_max,
        }

    # ── Radial bins (Fig. 4 style) ──
    E_bins = [
        ('0.03-0.1 GeV',  0.03,  0.1),
        ('0.1-0.31 GeV',  0.1,   0.31),
        ('0.31-1.0 GeV',  0.31,  1.0),
        ('1.0-3.1 GeV',   1.0,   3.1),
        ('3.1-10.0 GeV',  3.1,  10.0),
        ('10-100 GeV',   10.0, 100.0),
    ]

    theta_arr = np.linspace(5, 43, 16)    # 5° to 43° matching paper

    radial_data = {}
    for label, E_min, E_max in E_bins:
        flux_theory = model_500.radial_profile_in_energy_bin(E_min, E_max, theta_arr)

        noise_factor = 1.0 / np.sqrt(15)
        frac_err = noise_factor * (0.1 + 0.3 * (np.sqrt(E_min*E_max) / 30)**0.5)
        frac_err = np.clip(frac_err, 0.05, 3.0)

        flux_err = flux_theory * frac_err
        noise    = rng.normal(0, 1, size=len(theta_arr))
        flux_mock = np.maximum(flux_theory * (1 + frac_err * noise), 0.5)

        radial_data[label] = {
            'theta_deg':  theta_arr,
            'flux':       flux_mock,
            'flux_err':   flux_err,
            'E_min':      E_min,
            'E_max':      E_max,
        }

    # ── Helioprojective map (Fig. 2 style) ──
    Tx = np.linspace(-45, 45, 91)
    Ty = np.linspace(-45, 45, 91)
    TxG, TyG = np.meshgrid(Tx, Ty)
    r_map = np.sqrt(TxG**2 + TyG**2)

    # Solar halo: bright center, θ⁻¹ profile
    halo_map = np.where(r_map > 1, 30.0 / r_map, 30.0)

    # Background: Galactic plane structure (bright band along Ty≈0 for some exposures)
    bg_map = 200 + 100 * np.exp(-TyG**2 / 200) + rng.normal(0, 20, TyG.shape)

    # Mask solar circle
    mask = r_map > 45
    halo_map[mask] = 0
    bg_map[mask]   = 0

    return {
        'spectrum_data': spectrum_data,
        'radial_data':   radial_data,
        'helio_map': {
            'Tx': Tx, 'Ty': Ty,
            'data':       bg_map + halo_map + rng.normal(0, 5, bg_map.shape),
            'background': bg_map,
            'halo':       halo_map,
            'residual':   halo_map + rng.normal(0, 3, halo_map.shape),
        }
    }
