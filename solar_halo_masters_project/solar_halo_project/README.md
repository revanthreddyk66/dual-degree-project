# Solar ICS Halo — Fermi-LAT Analysis
### Master's Project Repository

> **Based on:** *"First Observations of Solar Halo Gamma Rays Over a Full Solar Cycle"*  
> Linden et al. (2026) — [arXiv:2505.04625](https://arxiv.org/abs/2505.04625)

---

## Overview

This repository contains a full, self-contained implementation of the solar Inverse-Compton Scattering (ICS) halo analysis using Fermi-LAT data. It is structured for a **1-year subset analysis** as a master's project, with all physics models, theoretical predictions, and plot-reproduction tools clearly documented.

The **solar ICS halo** arises when ambient cosmic-ray electrons and positrons (e⁺e⁻) scatter off sunlight photons, upscattering them to γ-ray energies via the Inverse-Compton process. The halo extends up to **45°** from the Sun and is detectable from **31.6 MeV to 100 GeV**.

---

##  Repository Structure

```
solar_halo_project/
│
├── README.md                        ← This file
├── requirements.txt                 ← Python dependencies
├── configs/
│   └── analysis_config.yaml         ← All tunable parameters
│
├── src/                             ← Core source code
│   ├── __init__.py
│   ├── solar_modulation.py          ← Force-field solar modulation models (I, II, III)
│   ├── ics_emissivity.py            ← ICS γ-ray emissivity calculations
│   ├── halo_model.py                ← Full solar halo model (spectrum + morphology)
│   ├── local_interstellar_spectrum.py  ← e+/e- LIS from Bisschoff+2019
│   ├── background_model.py          ← Data-driven astrophysical background
│   ├── statistical_fitting.py       ← Log-likelihood fitting with iminuit
│   └── utils.py                     ← Coordinate transforms, helper functions
│
├── models/                          ← Pre-computed model grids
│   └── README.md
│
│
├── notebooks/
│   └── 01_solar_halo_tutorial.ipynb 
│
├── data/
│   └── README.md                    ← Instructions for downloading Fermi-LAT data
│
├── docs/
│   ├── physics_derivations.md       ← Full mathematical derivations
│   └── pipeline_overview.md         ← Analysis pipeline explanation
│
```

---


### 1. Inverse-Compton Scattering (ICS)

The γ-ray emissivity from ICS is:

$$q_\gamma(r, E_\gamma) = \int dE_e \; J_{e^\pm}(r, E_e) \; \frac{d\sigma_\text{KN}}{dE_\gamma}(E_e, E_\gamma) \; n_\text{photon}(r)$$

where:
- $J_{e^\pm}(r, E_e)$ — cosmic-ray e⁺/e⁻ intensity at heliocentric distance $r$
- $\frac{d\sigma_\text{KN}}{dE_\gamma}$ — Klein-Nishina differential cross-section
- $n_\text{photon}(r) \propto r^{-2}$ — solar photon number density

### 2. Solar Modulation — Force-Field Approximation

The **Gleeson-Axford force-field model** (1968) relates the modulated spectrum at heliocentric distance $r$ to the Local Interstellar Spectrum (LIS):

$$J_{e^\pm}(r, E_e) = J_{e^\pm}(\infty, E_e) \cdot \frac{E_e^2 - E_0^2}{(E_e + e\Phi(r))^2 - E_0^2}$$

where:
- $J_{e^\pm}(\infty, E_e)$ — Local Interstellar Spectrum (from Bisschoff+2019/GALPROP)
- $E_0 = m_e c^2 = 0.511$ MeV — electron rest mass energy
- $e\Phi(r)$ — modulation potential at distance $r$

### 3. Modulation Potential Models

| Model | Description | Equation |
|-------|-------------|----------|
| **Model I** (default) | Radially-dependent, follows Solar Cycle 21 data | $\Phi_1(r) = \Phi_0 \left[ r^{-0.1} - r_b^{-0.1} \right] / \left[ (1\,\text{AU})^{-0.1} - r_b^{-0.1} \right]$ |
| **Model II** | No modulation inside 1 AU | $\Phi_2(r) = \Phi_0$ for $r < 1$ AU, else $\Phi_1(r)$ |
| **Model III** | Energy-dependent modulation | $\Phi_3(r, E_e) = \Phi_1(r) \times (E_e / 10\,\text{GeV})^{-\alpha}$ for $E_e < 10$ GeV |

where $r_b = 100$ AU is the heliopause boundary radius.

### 4. γ-ray Morphology

The line-of-sight integral of the ICS emissivity produces a surface brightness:

$$\frac{dN}{dE\,d\Omega}(\theta_s, E_\gamma) = \int_\text{l.o.s.} q_\gamma(r(\ell, \theta_s), E_\gamma) \, d\ell$$

At high energies with no modulation, this gives the simple analytic result:

$$\frac{dN}{dE\,d\Omega} \propto \theta_s^{-1}$$

---


### Installation

```bash
git clone https://github.com/YOUR_USERNAME/solar_halo_masters
cd solar_halo_masters
pip install -r requirements.txt
```

### Run theory plots (no Fermi data needed)

```bash
python plots/plot_fig3_spectral_bins.py --theory-only
python plots/plot_fig4_radial_bins.py --theory-only
```

### Run with 1-year Fermi data

```bash
# First configure your data paths:
nano configs/analysis_config.yaml

# Then run the full pipeline:
python plots/plot_all.py --config configs/analysis_config.yaml
```

---

## Dependencies

See `requirements.txt`. Key packages:
- `numpy`, `scipy`, `matplotlib` — numerical computing & plotting
- `astropy`, `healpy` — astronomical coordinates & HEALPix maps
- `iminuit` — likelihood minimization
- `fermitools` — Fermi-LAT data processing (requires separate install)
- `fermipy` — high-level Fermi analysis (requires separate install)
- `sunpy` — solar coordinate transforms

---




