"""
plot_all.py
===========
Reproduce all main figures from Linden et al. (2026) using 1-year mock data.

Figures produced:
    Fig 1  — Spectrum + morphology (main result)
    Fig 2  — Helioprojective maps (data / background / solar disk / solar halo)
    Fig 3  — Spectra in 8 radial bins
    Fig 4  — Radial profiles in 6 energy bins

Usage:
    python plots/plot_all.py

Output:
    plots/output/fig1_spectrum_morphology.pdf
    plots/output/fig2_helioprojective_maps.pdf
    plots/output/fig3_spectra_radial_bins.pdf
    plots/output/fig4_radial_profiles_energy_bins.pdf
"""

import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors as mcolors
from matplotlib.colors import LogNorm, Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable

from src import (
    SolarHaloModel, build_models_for_phi_values,
    CosmicRayLIS,
    generate_mock_1year_data,
    setup_paper_style,
    MODEL_COLORS, DATA_COLOR, DATA_MARKER, DATA_MS, DATA_CAPSIZE,
)

# ── Setup ─────────────────────────────────────────────────────────────────────
os.makedirs('plots/output', exist_ok=True)
setup_paper_style()

print("Initializing physics models...")
lis = CosmicRayLIS()
models = build_models_for_phi_values([0, 500, 1000], lis)
print("Generating 1-year mock data...")
data = generate_mock_1year_data(seed=2024)

PHI_LABELS  = {0: '0 MV', 500: '500 MV', 1000: '1000 MV'}
PHI_VALUES  = [0, 500, 1000]

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Main Result: Spectrum + Morphology
# ══════════════════════════════════════════════════════════════════════════════

def make_fig1():
    """
    Reproduce Fig. 1 from Linden et al. (2026).
    Top panel: Spectrum in 5-10° angular bin.
    Bottom panel: Radial profile at 1-3 GeV.
    """
    print("\n── Figure 1: Spectrum + Morphology ──")
    fig, axes = plt.subplots(2, 1, figsize=(4.5, 7.0))
    fig.subplots_adjust(hspace=0.35)

    # ── Top: Spectrum in 5–10° bin ──
    ax = axes[0]
    E_GeV = np.logspace(-1.5, 2, 200)    # 31.6 MeV – 100 GeV

    for phi in PHI_VALUES:
        m = models[phi]
        flux = m.spectrum_in_angular_bin(5, 10, E_GeV)
        color = MODEL_COLORS[phi]
        ax.loglog(E_GeV, flux, color=color, lw=1.8, label=PHI_LABELS[phi])

    # Mock data (5-10° bin)
    d = data['spectrum_data']['5-10']
    mask = ~d['upper_lim']
    ax.errorbar(d['E_GeV'][mask], d['flux'][mask], yerr=d['flux_err'][mask],
                fmt=DATA_MARKER, color=DATA_COLOR, ms=DATA_MS,
                capsize=DATA_CAPSIZE, lw=0.8, label='Data (1 yr mock)')

    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel(r'$E^2\,dN/dE$  ($10^{-8}$ GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
    ax.set_xlim(0.03, 100)
    ax.set_ylim(1, 100)
    ax.legend(fontsize=7.5, loc='upper right')
    ax.text(0.04, 0.92, '5–10°', transform=ax.transAxes,
            fontsize=9, ha='left', va='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='lightgray', pad=2))
    ax.set_title('Spectrum — Solar ICS Halo (1-year)', fontsize=10)

    # Minor ticks
    ax.xaxis.set_minor_locator(plt.LogLocator(subs='all'))

    # ── Bottom: Radial Profile at 1–3 GeV ──
    ax = axes[1]
    theta_arr = np.linspace(2, 45, 200)

    for phi in PHI_VALUES:
        m = models[phi]
        flux = m.radial_profile_in_energy_bin(1.0, 3.16, theta_arr)
        color = MODEL_COLORS[phi]
        ax.semilogy(theta_arr, flux, color=color, lw=1.8, label=PHI_LABELS[phi])

    # Mock data (radial, 1-3 GeV bin)
    d = data['radial_data']['1.0-3.1 GeV']
    ax.errorbar(d['theta_deg'], d['flux'], yerr=d['flux_err'],
                fmt=DATA_MARKER, color=DATA_COLOR, ms=DATA_MS,
                capsize=DATA_CAPSIZE, lw=0.8, label='Data (1 yr mock)')

    ax.set_xlabel('Angular Distance to the Sun (degrees)')
    ax.set_ylabel(r'$E^2\,dN/dE$  ($10^{-8}$ GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)')
    ax.set_xlim(3, 45)
    ax.set_ylim(0.5, 200)
    ax.legend(fontsize=7.5, loc='upper right')
    ax.text(0.96, 0.92, '1–3 GeV', transform=ax.transAxes,
            fontsize=9, ha='right', va='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='lightgray', pad=2))
    ax.set_title('Radial Profile — Solar ICS Halo (1-year)', fontsize=10)

    ax.xaxis.set_minor_locator(plt.MultipleLocator(5))

    plt.savefig('plots/output/fig1_spectrum_morphology.pdf', bbox_inches='tight')
    plt.savefig('plots/output/fig1_spectrum_morphology.png', bbox_inches='tight')
    print("   Saved: plots/output/fig1_spectrum_morphology.{pdf,png}")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Helioprojective Maps
# ══════════════════════════════════════════════════════════════════════════════

def make_fig2():
    """
    Reproduce Fig. 2 from Linden et al. (2026).
    6-panel helioprojective map: Data, Background, Solar Disk, Solar Halo,
    Halo+Residual (significance), Residual.
    """
    print("\n── Figure 2: Helioprojective Maps ──")

    hm = data['helio_map']
    Tx, Ty = hm['Tx'], hm['Ty']

    # Compute solar disk model
    TxG, TyG = np.meshgrid(Tx, Ty)
    r_map = np.sqrt(TxG**2 + TyG**2)
    disk_map = np.where(r_map < 2.5, 500.0, 0.0)   # solar disk within ~0.26°, PSF-broadened

    # Compute significance residual
    counts    = hm['data']
    bg        = hm['background']
    halo      = hm['halo']
    residual  = counts - bg - disk_map
    sigma     = np.where(counts > 0, residual / np.sqrt(np.abs(counts) + 1), 0)
    sigma     = np.clip(sigma, -4, 9)

    # Mask outside 45°
    mask_circle = r_map > 45

    fig = plt.figure(figsize=(8.5, 5.5))
    gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.08, wspace=0.12)

    panels = [
        ('Fermi-LAT Data\n1 – 3 GeV',  hm['data'],      'counts', 'viridis',  (50, 1000)),
        ('Background Model',             hm['background'], 'counts', 'viridis',  (50, 1000)),
        ('Solar Disk',                   disk_map,         'log',    'viridis',  (0.1, 600)),
        ('Solar Halo',                   halo,             'log',    'viridis',  (0.1, 30)),
        ('Halo + Residual',              sigma,            'linear', 'RdYlBu_r', (-4, 9)),
        ('Residual',                     sigma - halo/3,   'linear', 'RdYlBu_r', (-4, 9)),
    ]

    cbar_labels = {
        'counts': 'Counts/deg²',
        'log':    'Counts (arb. norm)',
        'linear': 'Significance (Residual / √Data)',
    }

    for i, (title, cdata, scale, cmap, clim) in enumerate(panels):
        row, col = divmod(i, 2)
        ax = fig.add_subplot(gs[row, col])

        plot_data = cdata.copy().astype(float)
        plot_data[mask_circle] = np.nan

        if scale == 'log':
            plot_data = np.clip(plot_data, clim[0], None)
            im = ax.imshow(plot_data, origin='lower',
                           extent=[-45, 45, -45, 45],
                           norm=LogNorm(vmin=clim[0], vmax=clim[1]),
                           cmap=cmap, aspect='equal')
        else:
            im = ax.imshow(plot_data, origin='lower',
                           extent=[-45, 45, -45, 45],
                           vmin=clim[0], vmax=clim[1],
                           cmap=cmap, aspect='equal')

        # Draw circle boundary
        theta_circ = np.linspace(0, 2*np.pi, 300)
        ax.plot(45*np.cos(theta_circ), 45*np.sin(theta_circ),
                'w-', lw=0.7, alpha=0.5)

        ax.text(0.04, 0.97, title, transform=ax.transAxes,
                fontsize=7.5, ha='left', va='top', color='white',
                fontweight='bold',
                bbox=dict(facecolor='black', alpha=0.4, edgecolor='none', pad=1.5))

        if row == 2:
            ax.set_xlabel('Helioprojective Longitude ($T_x$, °)', fontsize=7.5)
        else:
            ax.set_xticklabels([])

        if col == 0:
            ax.set_ylabel('Helioprojective Latitude ($T_y$, °)', fontsize=7.5)
        else:
            ax.set_yticklabels([])

        ax.tick_params(labelsize=6.5)
        ax.set_xlim(-47, 47)
        ax.set_ylim(-47, 47)

        # Colorbar
        plt.colorbar(im, ax=ax, fraction=0.04, pad=0.01,
                     label=cbar_labels[scale]).ax.tick_params(labelsize=6)

    fig.suptitle('Fermi-LAT Solar Halo Analysis — 1-3 GeV (1-year mock)', fontsize=9.5, y=1.01)

    plt.savefig('plots/output/fig2_helioprojective_maps.pdf', bbox_inches='tight')
    plt.savefig('plots/output/fig2_helioprojective_maps.png', bbox_inches='tight')
    print("   Saved: plots/output/fig2_helioprojective_maps.{pdf,png}")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Spectra in 8 Radial Bins
# ══════════════════════════════════════════════════════════════════════════════

def make_fig3():
    """
    Reproduce Fig. 3 from Linden et al. (2026).
    8-panel grid: E² dN/dE vs Energy in radial bins from <2.5° to 35-45°.
    """
    print("\n── Figure 3: Spectra in Radial Bins ──")

    theta_bins_ordered = [
        ('<2.5',  '<2.5°',  0.0,  2.5),
        ('2.5-5', '2.5°–5°', 2.5, 5.0),
        ('5-10',  '5°–10°',  5.0, 10.0),
        ('10-15', '10°–15°',10.0, 15.0),
        ('15-20', '15°–20°',15.0, 20.0),
        ('20-25', '20°–25°',20.0, 25.0),
        ('25-35', '25°–35°',25.0, 35.0),
        ('35-45', '35°–45°',35.0, 45.0),
    ]

    E_GeV_fine = np.logspace(-1.5, 2.0, 200)

    fig, axes = plt.subplots(4, 2, figsize=(7.0, 9.5))
    fig.subplots_adjust(hspace=0.06, wspace=0.06)

    axes_flat = axes.flatten(order='F')  # Column-major to match paper layout

    for idx, (key, label, th_min, th_max) in enumerate(theta_bins_ordered):
        ax = axes_flat[idx]

        # Theory curves (Model I solid, Model II dashed)
        for phi in PHI_VALUES:
            color = MODEL_COLORS[phi]
            # Model I
            flux_I = models[phi].spectrum_in_angular_bin(th_min, th_max, E_GeV_fine)
            ax.loglog(E_GeV_fine, flux_I, color=color, lw=1.6, ls='-')
            # Model II
            m_II = SolarHaloModel(Phi0_e_MV=phi, modulation_model='II', lis=lis)
            flux_II = m_II.spectrum_in_angular_bin(th_min, th_max, E_GeV_fine)
            ax.loglog(E_GeV_fine, flux_II, color=color, lw=1.2, ls='--')

        # Data
        if key in data['spectrum_data']:
            d = data['spectrum_data'][key]
            mask = ~d['upper_lim']
            ax.errorbar(d['E_GeV'][mask], d['flux'][mask], yerr=d['flux_err'][mask],
                        fmt=DATA_MARKER, color=DATA_COLOR, ms=3.5,
                        capsize=2.5, lw=0.8, zorder=5)
            # Upper limits for non-detected points
            ul_mask = d['upper_lim']
            if ul_mask.any():
                ax.errorbar(d['E_GeV'][ul_mask], d['flux'][ul_mask],
                            yerr=0.5*d['flux'][ul_mask],
                            fmt='v', color=DATA_COLOR, ms=3, lw=0.8,
                            uplims=True, alpha=0.5)

        ax.text(0.97, 0.96, label, transform=ax.transAxes,
                fontsize=8, ha='right', va='top')

        ax.set_xlim(0.05, 50)
        ax.set_ylim(1, 200)

        row, col = divmod(idx, 4)

        # Axis labels: only bottom row and left column
        if row == 3:
            ax.set_xlabel('Energy (GeV)', fontsize=8.5)
        else:
            ax.set_xticklabels([])

        if col == 0:
            ax.set_ylabel(r'$E^2\,dN/dE$  ($10^{-8}$ GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)',
                          fontsize=7.5)
        else:
            ax.set_yticklabels([])

        ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=7.5)

    # Legend in first panel
    ax0 = axes_flat[0]
    from matplotlib.lines import Line2D
    legend_elems = [
        Line2D([0], [0], color=MODEL_COLORS[0],    lw=1.6, label='0 MV'),
        Line2D([0], [0], color=MODEL_COLORS[500],  lw=1.6, label='500 MV'),
        Line2D([0], [0], color=MODEL_COLORS[1000], lw=1.6, label='1000 MV'),
        Line2D([0], [0], color='k', lw=0.8, ls='--', label='Model II'),
        Line2D([0], [0], marker=DATA_MARKER, color=DATA_COLOR,
               lw=0, ms=3.5, label='Data (1 yr)'),
    ]
    ax0.legend(handles=legend_elems, fontsize=6.5, loc='lower right',
               framealpha=0.85, handlelength=1.5)

    fig.suptitle('Solar ICS Halo Spectra in Radial Bins (1-year mock)',
                 fontsize=9.5, y=1.005)

    plt.savefig('plots/output/fig3_spectra_radial_bins.pdf', bbox_inches='tight')
    plt.savefig('plots/output/fig3_spectra_radial_bins.png', bbox_inches='tight')
    print("   Saved: plots/output/fig3_spectra_radial_bins.{pdf,png}")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Radial Profiles in Energy Bins
# ══════════════════════════════════════════════════════════════════════════════

def make_fig4():
    """
    Reproduce Fig. 4 from Linden et al. (2026).
    6-panel grid: E² dN/dE vs angular distance in 6 energy bins.
    """
    print("\n── Figure 4: Radial Profiles in Energy Bins ──")

    E_bins_ordered = [
        ('0.03-0.1 GeV',  '0.03 – 0.1 GeV',  0.03,  0.1,   2),
        ('0.1-0.31 GeV',  '0.1 – 0.31 GeV',   0.1,   0.31,  2),
        ('0.31-1.0 GeV',  '0.31 – 1.0 GeV',   0.31,  1.0,   2),
        ('1.0-3.1 GeV',   '1.0 – 3.1 GeV',    1.0,   3.1,   2),
        ('3.1-10.0 GeV',  '3.1 – 10.0 GeV',   3.1,  10.0,   2),
        ('10-100 GeV',    '10 – 100 GeV',     10.0, 100.0,   2),
    ]

    theta_fine = np.linspace(3, 45, 200)

    fig, axes = plt.subplots(3, 2, figsize=(7.0, 7.5))
    fig.subplots_adjust(hspace=0.06, wspace=0.06)

    axes_flat = axes.flatten(order='F')

    for idx, (key, label, E_min, E_max, _) in enumerate(E_bins_ordered):
        ax = axes_flat[idx]

        # Theory curves
        for phi in PHI_VALUES:
            color = MODEL_COLORS[phi]
            # Model I solid
            flux_I = models[phi].radial_profile_in_energy_bin(E_min, E_max, theta_fine)
            ax.semilogy(theta_fine, flux_I, color=color, lw=1.6, ls='-')
            # Model II dashed
            m_II = SolarHaloModel(Phi0_e_MV=phi, modulation_model='II', lis=lis)
            flux_II = m_II.radial_profile_in_energy_bin(E_min, E_max, theta_fine)
            ax.semilogy(theta_fine, flux_II, color=color, lw=1.2, ls='--')

        # Data
        if key in data['radial_data']:
            d = data['radial_data'][key]
            ax.errorbar(d['theta_deg'], d['flux'], yerr=d['flux_err'],
                        fmt=DATA_MARKER, color=DATA_COLOR, ms=3.5,
                        capsize=2.5, lw=0.8, zorder=5)

        ax.text(0.97, 0.96, label, transform=ax.transAxes,
                fontsize=7.5, ha='right', va='top')

        ax.set_xlim(3, 45)
        ax.set_ylim(0.8, 200)
        ax.xaxis.set_minor_locator(plt.MultipleLocator(5))
        ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=7.5)

        row, col = divmod(idx, 3)

        if row == 2:
            ax.set_xlabel('Angular Distance to the Sun (degrees)', fontsize=8.5)
        else:
            ax.set_xticklabels([])

        if col == 0:
            ax.set_ylabel(r'$E^2\,dN/dE$  ($10^{-8}$ GeV cm$^{-2}$ s$^{-1}$ sr$^{-1}$)',
                          fontsize=7.5)
        else:
            ax.set_yticklabels([])

    # Legend
    from matplotlib.lines import Line2D
    ax0 = axes_flat[2]
    legend_elems = [
        Line2D([0], [0], color=MODEL_COLORS[0],    lw=1.6, label='0 MV'),
        Line2D([0], [0], color=MODEL_COLORS[500],  lw=1.6, label='500 MV'),
        Line2D([0], [0], color=MODEL_COLORS[1000], lw=1.6, label='1000 MV'),
        Line2D([0], [0], marker=DATA_MARKER, color=DATA_COLOR,
               lw=0, ms=3.5, label='Data (1 yr)'),
    ]
    ax0.legend(handles=legend_elems, fontsize=6.5, loc='lower left',
               framealpha=0.85, handlelength=1.5)

    fig.suptitle('Solar ICS Halo Radial Profiles in Energy Bins (1-year mock)',
                 fontsize=9.5, y=1.005)

    plt.savefig('plots/output/fig4_radial_profiles_energy_bins.pdf', bbox_inches='tight')
    plt.savefig('plots/output/fig4_radial_profiles_energy_bins.png', bbox_inches='tight')
    print("   Saved: plots/output/fig4_radial_profiles_energy_bins.{pdf,png}")
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 60)
    print("  Solar ICS Halo — Figure Reproduction (1-year dataset)")
    print("  Based on Linden et al. (2026), arXiv:2505.04625")
    print("=" * 60)

    fig1 = make_fig1()
    fig2 = make_fig2()
    fig3 = make_fig3()
    fig4 = make_fig4()

    print("\n✓ All figures saved to plots/output/")
    print("  Using 1-year mock data (scaled from 15-year paper results)")
    plt.close('all')
