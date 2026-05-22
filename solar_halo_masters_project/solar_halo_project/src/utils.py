"""
utils.py
========
Utility functions for coordinate transforms, plotting helpers, and
data I/O for the solar ICS halo analysis.
"""

import numpy as np
import matplotlib.pyplot as plt

# ── Coordinate Transforms ────────────────────────────────────────────

def angular_separation_deg(ra1, dec1, ra2, dec2):
    """Great–circle separation between two (RA,Dec) points in degrees."""
    ra1, dec1, ra2, dec2 = map(np.radians, [ra1, dec1, ra2, dec2])
    ddec = dec2 - dec1
    dra  = ra2 - ra1
    a = (np.sin(ddec/2)**2 +
         np.cos(dec1)*np.cos(dec2)*np.sin(dra/2)**2)
    return np.degrees(2*np.arctan2(np.sqrt(a), np.sqrt(1-a)))


def heliocentric_distance_along_los(theta_s_deg, ell_AU, d_obs_AU=1.0):
    """
    Heliocentric distance of a point that is ell_AU along the line of
    sight seen at an angle theta_s_deg from the Sun.
    """
    theta = np.radians(theta_s_deg)
    r2 = (ell_AU**2 + d_obs_AU**2 -
          2*ell_AU*d_obs_AU*np.cos(np.pi - theta))
    return np.sqrt(np.maximum(r2, 1e-10))


def los_integration_limits(theta_s_deg, d_obs_AU=1.0, max_dist_AU=50.0):
    """Return (ell_min, ell_max) integration limits in AU."""
    return 0.01, max_dist_AU          # 0.01 AU ≃ just outside the Sun

# ── Plotting Utilities ───────────────────────────────────────────────

MODEL_COLORS = {0: '#E85D3C', 500: '#2CA6A4', 1000: '#2C3E50'}
MODEL_LINESTYLES = {'I': '-', 'II': '--'}
DATA_COLOR  = 'black'
DATA_MARKER = 'D'
DATA_MS     = 4
DATA_CAPSIZE = 3


def setup_paper_style():
    """Matplotlib rcParams matching the paper figures."""
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
    """Convenience function to add a model legend to an Axes."""
    handles = []
    for phi in phi_values:
        c = MODEL_COLORS.get(phi, 'gray')
        for m in model_types:
            ls = MODEL_LINESTYLES.get(m, '-')
            h, = ax.plot([], [], color=c, ls=ls, lw=1.8, label=f'{phi} MV')
            handles.append(h)
    if len(model_types) > 1:
        ax.plot([], [], DATA_MARKER, color=DATA_COLOR, ms=DATA_MS, label='Data')
    return ax.legend(handles=handles, frameon=True, framealpha=0.8,
                     edgecolor='gray', handlelength=2.0)


def add_energy_bin_label(ax, E_min, E_max, unit='GeV',
                          loc='upper right', fontsize=9):
    ax.text(0.97, 0.95, f'{E_min} – {E_max} {unit}', transform=ax.transAxes,
            ha='right', va='top', fontsize=fontsize,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))


def add_angle_bin_label(ax, theta_min, theta_max, unit='°',
                         loc='upper right', fontsize=9):
    ax.text(0.97, 0.95, f'{theta_min}{unit}–{theta_max}{unit}',
            transform=ax.transAxes,
            ha='right', va='top', fontsize=fontsize,
            bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1))
