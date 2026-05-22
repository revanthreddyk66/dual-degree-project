# Physics Derivations — Solar ICS Halo

## 1. Physical Setup

### The Solar ICS Halo

The Sun is surrounded by a diffuse glow of γ-rays produced by **Inverse-Compton Scattering** (ICS): ambient cosmic-ray electrons and positrons (e⁺e⁻) collide with solar photons and upscatter them to MeV–GeV energies.

**Two competing effects:**
- At **high energy** (>10 GeV): electrons free-stream; halo morphology ∝ θ⁻¹
- At **low energy** (<1 GeV): solar wind magnetic field modulates (suppresses) the e⁺e⁻ flux near the Sun → *solar modulation*

---

## 2. Inverse-Compton Scattering

### 2.1 Basic Kinematics

In the Thomson limit (when the photon energy in the electron rest frame ε' ≪ m_e c²), the upscattered photon energy is:

$$\varepsilon_\gamma \approx \frac{4}{3} \gamma^2 \varepsilon_\odot (1 - \beta\cos\alpha)$$

where:
- γ = E_e / m_e c² is the electron Lorentz factor
- ε⊙ ≈ 1 eV is the typical solar photon energy  
- α is the angle between the electron velocity and photon direction

**Peak γ-ray energy:** For a GeV electron (γ ≈ 2000), scattering a visible photon (ε⊙ ≈ 2 eV):

$$\varepsilon_\gamma^\text{peak} \approx 4\gamma^2 \varepsilon_\odot \sim 4 \times 10^6 \times 2 \text{ eV} \sim \text{few GeV}$$

This explains why the solar ICS halo peaks around 1–3 GeV.

### 2.2 Klein-Nishina Cross-Section

The full differential ICS cross-section (Klein & Nishina 1929) is:

$$\frac{d\sigma_\text{KN}}{dE_\gamma}(E_e, \varepsilon, E_\gamma) = \frac{3\sigma_T}{4\gamma^2 \varepsilon} \left[ 2q\ln q + (1+2q)(1-q) + \frac{(\Gamma q)^2(1-q)}{2(1+\Gamma q)} \right]$$

where:
- σ_T = 6.652 × 10⁻²⁵ cm² (Thomson cross-section)
- Γ = 4εγ / (m_e c²) is the "Compton parameter"
- q = E_γ / [Γ(E_e - E_γ)] ∈ [1/(1+Γ), 1]

**Thomson limit:** Γ ≪ 1 → reduces to σ_T × (kinematic factor)  
**Klein-Nishina regime:** Γ ≫ 1 → cross-section suppressed at high energies

### 2.3 Solar Photon Density

Solar photons are a diluted blackbody:

$$n_\odot(r, \varepsilon) = W(r) \cdot \frac{8\pi}{(hc)^3} \cdot \frac{\varepsilon^2}{e^{\varepsilon/kT_\odot} - 1}$$

where the **dilution factor** W(r) accounts for the finite size of the Sun:

$$W(r) = \frac{1}{2} \left(\frac{R_\odot}{r}\right)^2 \propto r^{-2}$$

with T⊙ = 5778 K and R⊙ = 6.957 × 10¹⁰ cm.

---

## 3. ICS Emissivity

The γ-ray emissivity (power per unit volume per unit solid angle per unit energy) is:

$$q_\gamma(r, E_\gamma) = c \int_0^\infty d\varepsilon \int_{E_{\min}}^\infty dE_e \; [J_{e^-}(r, E_e) + J_{e^+}(r, E_e)] \; \frac{d\sigma_\text{KN}}{dE_\gamma} \; n_\odot(r, \varepsilon)$$

**Key dependencies:**
- Scales as `r⁻²` from the solar photon density
- Modified near the Sun by solar modulation of J_e±

---

## 4. Solar Halo Surface Brightness

### 4.1 Line-of-Sight Integral

The observed γ-ray surface brightness at angular distance θ_s from the Sun:

$$\frac{dN}{dE_\gamma d\Omega}(\theta_s, E_\gamma) = \int_{-\infty}^{+\infty} q_\gamma\left(r(\ell, \theta_s), E_\gamma\right) d\ell$$

where r(ℓ, θ_s) is the heliocentric distance at depth ℓ along the line of sight:

$$r^2 = \ell^2 + d_\oplus^2 - 2\ell d_\oplus \cos(\pi - \theta_s)$$

with d⊕ = 1 AU (Earth-Sun distance).

### 4.2 Analytic Result (High-Energy Limit)

At E_γ ≫ 1 GeV, solar modulation is negligible and the e+e⁻ density is nearly constant in the heliosphere. Then:

$$n_\odot(r) \propto r^{-2} \quad \Rightarrow \quad q_\gamma(r) \propto r^{-2}$$

Integrating r⁻² over the line of sight at angle θ_s gives:

$$\frac{dN}{dE_\gamma d\Omega} \propto \theta_s^{-1}$$

**This is the famous 1/θ morphology of the solar ICS halo.**

---

## 5. Solar Modulation — Force-Field Approximation

### 5.1 Physical Picture

The **heliospheric magnetic field** (HMF), carried by the solar wind, creates a turbulent barrier that scatters and decelerates low-energy cosmic rays entering the heliosphere. The net effect is a reduction in flux and softening of the spectrum near the Sun.

### 5.2 Gleeson-Axford Force-Field Equation

The steady-state diffusion-convection transport equation in the heliosphere, under spherical symmetry, reduces to (Gleeson & Axford 1968):

$$\frac{\partial J}{\partial t} + \nabla \cdot \mathbf{S} = 0$$

where the streaming flux S depends on the diffusion coefficient κ and solar wind velocity V. The **force-field approximation** solves this exactly when κ ∝ βp (rigidity times velocity), giving:

$$\boxed{J(r, E) = J_\text{LIS}(\infty, E) \cdot \frac{E^2 - E_0^2}{(E + e\Phi(r))^2 - E_0^2}}$$

where:
- J_LIS is the Local Interstellar Spectrum (unmodulated)
- E is the total kinetic energy
- E₀ = m_e c² = 0.511 MeV is the rest mass energy
- eΦ(r) is the **modulation potential** in MeV (or MV for unit charge particles)

**Physical interpretation:** Each particle loses energy eΦ(r) adiabatically as it propagates from the heliospheric boundary (at r_b ≈ 100 AU) to heliocentric distance r. The modulation potential integrates all energy losses along the path.

### 5.3 Radial Dependence of Φ(r)

**Model I** (Fujii & McDonald 2005, from Solar Cycle 21 spacecraft data):

$$\Phi_1(r) = \Phi_0 \cdot \frac{r^{-0.1} - r_b^{-0.1}}{(1\,\text{AU})^{-0.1} - r_b^{-0.1}}$$

where r_b = 100 AU. This form comes from the radial dependence of the cosmic-ray mean free path λ ∝ r^{0.9} (implying κ ∝ r^{0.9} βp), consistent with a 1/f magnetic power spectrum.

**Model II** (conservative limit):

$$\Phi_2(r) = \begin{cases} \Phi_0 & r < 1\text{ AU} \\ \Phi_1(r) & r \geq 1\text{ AU} \end{cases}$$

Motivated by Parker Solar Probe data (Li et al. 2022) suggesting weak modulation inside 1 AU.

**Model III** (energy-dependent):

$$\Phi_3(r, E_e) = \begin{cases} \Phi_1(r) \cdot (E_e / 10\,\text{GeV})^{-\alpha} & E_e < 10\text{ GeV} \\ \Phi_1(r) & E_e \geq 10\text{ GeV} \end{cases}$$

Allows the modulation to be stronger at lower energies. **Best-fit from paper: α = 0.0 ± 0.17** (no evidence for energy dependence).

---

## 6. Local Interstellar Spectrum

The unmodulated e⁺/e⁻ spectrum outside the heliosphere, from Bisschoff, Potgieter & Aslam (2019), computed with GALPROP fitting Voyager 1 + PAMELA data. At GeV energies:

$$J_{e^-}(\infty, E) \approx C_- \left(\frac{E}{1\,\text{GeV}}\right)^{-3.28}$$

$$J_{e^+}(\infty, E) \approx C_+ \left(\frac{E}{1\,\text{GeV}}\right)^{-2.85}$$

The positron fraction is ~5-10% at GeV energies and rising at higher energies (AMS-02 positron excess).

---

## 7. Statistical Fitting

### 7.1 Log-Likelihood

The analysis minimizes the Poisson log-likelihood:

$$-2\ln\mathcal{L} = 2\sum_i \left[ \mu_i - n_i + n_i \ln(n_i/\mu_i) \right]$$

where n_i are observed counts and μ_i are predicted counts in each HEALPix pixel and energy bin.

### 7.2 Three-Component Model

For each HEALPix pixel and energy bin, the predicted counts are:

$$\mu_i = N_\text{bg} \cdot B_i + N_\text{disk} \cdot D_i + N_\text{halo}^{(1)} \cdot H_i^{(\Phi=0)} + N_\text{halo}^{(2)} \cdot H_i^{(\Phi=1000)}$$

where:
- B_i = data-driven astrophysical background
- D_i = solar disk template (PSF-convolved)  
- H_i^{(Φ)} = solar halo model templates at two extreme modulation potentials

### 7.3 Best-Fit Results (Linden+2026, 15 years)

| Parameter | Best-fit | 5σ range |
|-----------|----------|----------|
| Φ₀(e⁻) [MV] | 475 | 415 – 530 |
| Φ₀(e⁺) [MV] | 0 | 0 – 223 |
| Detection significance | >100σ | — |
| Energy range | 31.6 MeV – 100 GeV | — |
| Angular range | 0 – 45° | — |

---

## 8. Key Observational Results

### 8.1 Solar Modulation Time Variation

By splitting the 15-year dataset into yearly bins, the paper provides the **first γ-ray measurement of solar modulation as a function of time**, finding:
- Minimum modulation during 2008 and 2019 solar minima (~200-300 MV)
- Maximum modulation during 2014 solar maximum (~600-800 MV)
- Consistent with PAMELA and AMS-02 local electron measurements

### 8.2 Azimuthal Asymmetry

Evidence for **different modulation potentials** along (Tx) vs across (Ty) the ecliptic plane:
- Φ₀(Tx) ≈ 400 MV (ecliptic plane — heliospheric current sheet direction)
- Φ₀(Ty) ≈ 550 MV (perpendicular — polar drift direction)
- Preferred at ~4.7σ; possibly related to Solar Cycle 24 polarity

---


