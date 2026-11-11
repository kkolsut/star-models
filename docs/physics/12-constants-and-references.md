# 12. Constants, symbols and references

## Constants

Implemented, from `crates/stellar-core/src/constants.rs` (CGS).

| Constant | Symbol | Value | Unit |
|---|---|---|---|
| Gravitational constant | $`G`$ | $`6.67430 \times 10^{-8}`$ | cm³ g⁻¹ s⁻² |
| Speed of light | $`c`$ | $`2.99792458 \times 10^{10}`$ | cm s⁻¹ |
| Stefan–Boltzmann constant | $`\sigma`$ | $`5.670374419 \times 10^{-5}`$ | erg cm⁻² s⁻¹ K⁻⁴ |
| Radiation constant | $`a = 4\sigma/c`$ | $`7.56573 \times 10^{-15}`$ | erg cm⁻³ K⁻⁴ |
| Boltzmann constant | $`k`$ | $`1.380649 \times 10^{-16}`$ | erg K⁻¹ |
| Atomic mass unit | $`m_u`$ | $`1.66053907 \times 10^{-24}`$ | g |
| Solar mass | $`M_\odot`$ | $`1.98847 \times 10^{33}`$ | g |
| Solar radius (IAU 2015) | $`R_\odot`$ | $`6.957 \times 10^{10}`$ | cm |
| Solar luminosity (IAU 2015) | $`L_\odot`$ | $`3.828 \times 10^{33}`$ | erg s⁻¹ |
| Solar effective temperature | $`T_{\mathrm{eff},\odot}`$ | 5772 | K |
| Julian year | | $`3.15576 \times 10^{7}`$ | s |

Needed by planned upgrades (CODATA values, to be added to `constants.rs`):

| Constant | Symbol | Value | Unit | Used by |
|---|---|---|---|---|
| Planck constant | $`h`$ | $`6.62607015 \times 10^{-27}`$ | erg s | ionisation |
| Reduced Planck constant | $`\hbar`$ | $`1.054571817 \times 10^{-27}`$ | erg s | degeneracy |
| Electron mass | $`m_e`$ | $`9.1093837 \times 10^{-28}`$ | g | degeneracy, ionisation |
| Electronvolt | eV | $`1.602176634 \times 10^{-12}`$ | erg | ionisation, reaction energies |
| Elementary charge | $`e`$ | $`4.80320 \times 10^{-10}`$ | esu | Coulomb coupling diagnostic |

## Symbols

| Symbol | Meaning |
|---|---|
| $`m`$, $`M`$, $`q = m/M`$ | enclosed mass, total mass, mass fraction |
| $`r`$, $`R`$ | radius, radius of the photosphere |
| $`P`$, $`T`$, $`\rho`$ | pressure, temperature, density |
| $`L`$ | luminosity through the sphere of radius $`r`$ |
| $`X`$, $`Y`$, $`Z`$ | mass fractions of hydrogen, helium, heavier elements |
| $`X_{\mathrm{CNO}}`$, $`f_{\mathrm{CNO}}`$ | mass fraction of C+N+O; its share of $`Z`$ |
| $`\mu`$, $`\mu_e`$, $`\mu_i`$ | mean molecular weight: all particles, per electron, per ion |
| $`\beta`$ | gas pressure as a fraction of the total |
| $`\alpha`$, $`\delta`$, $`\varphi`$ | $`(\partial\ln\rho/\partial\ln P)_T`$, $`-(\partial\ln\rho/\partial\ln T)_P`$, $`(\partial\ln\rho/\partial\ln\mu)_{P,T}`$ |
| $`c_P`$, $`h`$, $`s`$ | specific heat at constant pressure, specific enthalpy, specific entropy |
| $`\Gamma_1`$ | first adiabatic exponent |
| $`\kappa`$ | Rosseland mean opacity |
| $`\varepsilon_{pp}`$, $`\varepsilon_{\mathrm{CNO}}`$, $`\varepsilon_{3\alpha}`$, $`\varepsilon_{\mathrm{grav}}`$, $`\varepsilon_\nu`$ | energy generation and loss rates per gram |
| $`Q`$ | energy released per gram of fuel |
| $`f_{Z_1Z_2}`$, $`\zeta`$ | screening factor and its composition parameter |
| $`\psi`$, $`\phi`$ | pp branching factors (energy, helium production) |
| $`\nabla`$, $`\nabla_{\mathrm{rad}}`$, $`\nabla_{\mathrm{ad}}`$, $`\nabla_{\mathrm{L}}`$, $`\nabla_\mu`$ | temperature gradients $`d\ln T/d\ln P`$: actual, radiative, adiabatic, Ledoux; and $`d\ln\mu/d\ln P`$ |
| $`g`$, $`H_P`$, $`\ell`$, $`\alpha_{\mathrm{MLT}}`$ | gravity, pressure scale height, mixing length and its parameter |
| $`U`$, $`W`$, $`\xi`$, $`\eta`$ | quantities of the mixing-length cubic (page 5 only) |
| $`v_{\mathrm{conv}}`$, $`D`$, $`\alpha_{\mathrm{ov}}`$ | convective velocity, diffusion coefficient, overshoot parameter |
| $`\tau`$, $`q(\tau)`$, $`T_{\mathrm{eff}}`$ | optical depth, Hopf function, effective temperature |
| $`\Omega`$, $`\omega = \Omega/\Omega_{\mathrm{crit}}`$, $`j`$, $`J`$, $`f_P`$ | rotation rate, its ratio to break-up, specific and total angular momentum, centrifugal factor (page 8) |
| $`\omega`$, $`\Omega^2`$, $`\Pi`$ | pulsation angular frequency, its dimensionless square, period (page 9) |
| $`\xi`$, $`p`$, $`\xi_r`$, $`p'`$ | pulsation displacement and pressure perturbations (page 9) |
| $`c_s`$, $`S_l`$, $`N`$ | sound speed, Lamb frequency, buoyancy frequency |
| $`\nu_{\max}`$, $`\Delta\nu`$, $`\Delta\Pi_l`$ | frequency of maximum power, large separation, g-mode period spacing |
| $`\dot M`$ | mass-loss rate (negative) |

Three symbols are reused and are distinguished by page: $`\xi`$ (mixing-length cubic, Lane–Emden
radius, pulsation displacement), $`\omega`$ and $`\Omega`$ (rotation, pulsation), and $`\alpha`$
(EOS derivative, pp branching parameter).

## References

Implemented physics and methods:

* **Kippenhahn & Weigert (1990)**, *Stellar Structure and Evolution*, Springer: nuclear rate fits in
  use.
* **Kippenhahn, Weigert & Weiss (2012)**, *Stellar Structure and Evolution*, 2nd ed., Springer:
  structure equations, EOS, mixing-length formulation, Henyey method, polytropes.
* **Carroll & Ostlie (2017)**, *An Introduction to Modern Astrophysics*, 2nd ed.: analytic opacity
  formulas, the StatStar model.
* **Böhm-Vitense (1958)**, Z. Astrophys. 46, 108: mixing-length theory.
* **Henyey, Forbes & Gould (1964)**, ApJ 139, 306: relaxation method.
* **Clayton (1968)**, *Principles of Stellar Evolution and Nucleosynthesis*: weak screening, pp
  branching.
* **Cox (1980)**, *Theory of Stellar Pulsation*: linear adiabatic wave equation.
* **Brown et al. (1991)**, ApJ 368, 599; **Kjeldsen & Bedding (1995)**, A&A 293, 87: scaling
  relations.
* **Asplund et al. (2009)**, ARA&A 47, 481; **IAU 2015 Resolution B3**: solar composition preset and
  nominal solar values.

Planned upgrades (each to be checked when used):

* **Paczyński (1983)**, ApJ 267, 315: degenerate electron pressure interpolation.
* **Krishna Swamy (1966)**, ApJ 145, 174: solar $`T(\tau)`$ relation.
* **de Jager, Nieuwenhuijzen & van der Hucht (1988)**, A&AS 72, 259; **Reimers (1975)**;
  **Vink, de Koter & Lamers (2001)**, A&A 369, 574: mass loss.
* **von Zeipel (1924)**, MNRAS 84, 665; **Kippenhahn & Thomas (1970)**; **Endal & Sofia (1976)**,
  ApJ 210, 184: rotation.
* **Unno, Osaki, Ando, Saio & Shibahashi (1989)**, *Nonradial Oscillations of Stars*, 2nd ed.:
  non-radial modes.
* **Herwig (2000)**, A&A 360, 952: exponential overshoot.
* **Thoul, Bahcall & Loeb (1994)**, ApJ 421, 828: element diffusion.
* **Itoh et al. (1996)**, ApJS 102, 411: neutrino losses.
* **Caughlan & Fowler (1988)**, At. Data Nucl. Data Tables 40, 283: reaction rates.
* **Iglesias & Rogers (1996)**, ApJ 464, 943 (OPAL); **Ferguson et al. (2005)**, ApJ 623, 585:
  opacity tables.
* **Eggleton (1971)**, MNRAS 151, 351: bounded scale height.
* **Acharya et al. (2025)**, Rev. Mod. Phys. ("Solar Fusion III"): current $`S`$ factors for the pp
  chain and CNO cycle.
* **Ledoux (1951)**, ApJ 114, 373: rotational splitting constant.
* **Pedersen et al. (2018)**, A&A 614, A128: comparison of overshoot prescriptions.
* **Vinyoles et al. (2017)**, ApJ 835, 202: standard solar model calibration.
* **Paxton et al. (2011–2019)**, ApJS (the MESA instrument papers): free-energy EOS blending,
  energy equation, rotation factors $`f_P`$, $`f_T`$.

* **Trampedach et al. (2014)**, MNRAS 442, 805: Hopf-corrected radiative gradient and flux
  relations (§2.5, eqs. 35–38); range of the Krishna Swamy relation (§5).
* **Henyey, Vardya & Bodenheimer (1965)**, ApJ 142, 841: convection with optically thin elements
  (eqs. 34–36, 39–40).
* **Chugunov & DeWitt (2009)**, Phys. Rev. C 80, 014611: screening in ionic mixtures, candidate for
  the strong regime.
* **Graboske et al. (1973)**, ApJ 181, 457: screening regimes (not accessed).
* **Reese (2015)**: rotational splitting, eq. 3.
* **Timmes & Arnett (1999)**, ApJS 125, 277; the pynucastro EOS documentation: Fermi–Dirac
  electron gas.

### Source comparisons recorded (external review, 2026-11-11)

| Item | Page | Result | Source location |
|---|---|---|---|
| Fermi–Dirac $`n_e`$, $`P_e`$, $`u_e`$ | 2 | integral definition and normalisation agree, including the $`\vartheta F_{3/2}`$ term; pressure and energy prefactors re-derived by change of variable | pynucastro EOS derivation (which notes a missing factor in the printed Timmes–Arnett expression) |
| EOS consistency residual | 2 | free-energy approach and fixed-composition identity supported | Timmes's thermodynamics reference, eqs. 1–8 |
| General buoyancy term $`B`$ | 5, 9 | equivalent; $`B = \delta\,B_{\mathrm{MESA}}`$ | Paxton et al. 2013, §3.3, eqs. 5–8 |
| Exponential overshoot | 5 | factor 2, scale, normalisation and cut-off supported; Herwig's cut-off is $`10^{-2}`$ cm² s⁻¹. The inward reference point $`f_0`$ is the MESA convention, not Herwig's. The value $`f_0 = 0.005`$ is unverified | Herwig 2000, §2, eq. 2; MESA controls documentation (`overshoot_f0`) |
| Hopf gradient and flux fractions | 6 | agree, within that construction's assumptions; does not cover the modified mixing-length cubic | Trampedach et al. 2014, eqs. 35–38 |
| Ionisation equilibrium | 2 | chemical-potential relation and classical limit re-derived; signs, mass ratio and ground-state factors consistent | derivation |
| Mass-coordinate diffusion equation | 5 | checked analytically | derivation |
| Pressure-ionisation potential and constants | 2 | agree with the reproduction; original Pols et al. 1995 not read. $`\mu`$, $`P`$, $`u`$ from it are derivations | Kovetz, Yaron & Prialnik 2009, §3.1, eq. 28 |
| Screening: plasma free energy, $`H_{ij}`$, quantum correction | 4 | agree. The electron-response extension is not from the source | Chugunov & DeWitt 2009, eq. 24, eq. A4, §§II, IV |
| Triple-α screening basis | 4 | agrees | Kushnir, Waxman & Chugunov 2019, eqs. 1, 7–8 |
| Thin-element efficiency | 6 | agrees; $`U_e`$ is a derivation from it | Henyey, Vardya & Bodenheimer 1965, eqs. 32, 34–36, 39–40 |
| Cowling equations | 9 | checked analytically from continuity, momentum and adiabatic relations; not compared with Unno et al. | derivation |
| Ledoux constant | 8 | integral agrees | Reese 2015, eq. 3 |
| Krishna Swamy coefficients | 6 | agree with a contemporaneous reproduction (eq. 4); original not read; range 0.02–10 in $`\tau`$ | Trampedach et al. 2014, §5 |
| pp branching: $`S`$ factors and $`C_\alpha \simeq 1.49\times10^{17}`$ | 4 | leading-order estimate only | Solar Fusion III, Table 1 |

Still unverified (gates stay): the Kippenhahn–Weigert(–Weiss) rate coefficients, the attribution of
$`C_\alpha`$ to Clayton, the Paczyński formulas, the Graboske screening prescription, the form of
pressure ionisation, the mass-loss fits, the gravity-darkening exponents, the Cowling equations as
printed in Unno et al., the neutrino fits, the bounded scale height, the opacity-table details, the
helium-enrichment parameters, and for Krishna Swamy's relation the original coefficients and its
optical-depth scale. The modified mixing-length cubic of page 6 is a derivation of these pages and
their reviews, with no published source compared.

None of the reviews looked at the code, re-ran the measured numbers or tested numerical
convergence.

Journal references in this list are quoted from memory; confirm volume and page before citing them
elsewhere.

---
Previous: [Validity, limits and order of work](11-validity-and-roadmap.md) · [Contents](README.md)
