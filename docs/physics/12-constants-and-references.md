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
* **Cox (1980)**, *Theory of Stellar Pulsation*: linear adiabatic wave equation; classical
  polytrope periods ($`\Omega^2 \approx 2.712`$ for $`n = 1.5`$, 9.26 for $`n = 3`$, quoted from memory
  for the test of plan 1.3, **[check source]**).
* **Pekeris (1938)**, ApJ 88, 189: radial modes of the uniform compressible sphere,
  $`\Omega_k^2 = \Gamma_1(k+1)(2k+3) - 4`$ (reproduced by the solver to $`10^{-3}`$, plan 1.3).
* **Brown et al. (1991)**, ApJ 368, 599; **Kjeldsen & Bedding (1995)**, A&A 293, 87: scaling
  relations.
* **Asplund et al. (2009)**, ARA&A 47, 481; **IAU 2015 Resolution B3**: solar composition preset and
  nominal solar values.

Planned upgrades (each to be checked when used):

* **Paczyński (1983)**, ApJ 267, 315–321: degenerate electron pressure interpolation (eq. 4,
  p. 317).
* **Pols, Tout, Eggleton & Han (1995)**, MNRAS 274, 964–974: ionisation and pressure ionisation
  (eqs. 14–18, 26–29).
* **Kovetz, Yaron & Prialnik (2009)**, MNRAS 395, 1857–1874: reproduction of the pressure-ionisation
  term (§3.1, eq. 28).
* **NIST Atomic Spectra Database** (ionisation energies, retrieved 2026-11-11).
* **Krishna Swamy (1966)**, ApJ 145, 174: solar $`T(\tau)`$ relation (eq. 33, p. 182; Rosseland
  mean optical depth).
* **de Jager, Nieuwenhuijzen & van der Hucht (1988)**, A&AS 72, 259–289; **Reimers (1975)**, Mém.
  Soc. R. Sci. Liège, 6e sér., 8, 369 (original not read; quoted in **Reimers (1977)**, A&A 61,
  217); **Vink, de Koter & Lamers (2001)**, A&A 369, 574–588: mass loss.
* **von Zeipel (1924)**, MNRAS 84, 665; **Lucy (1967)**, Z. Astrophys. 65, 89: gravity darkening.
  **Kippenhahn & Thomas (1970)**; **Endal & Sofia (1976)**, ApJ 210, 184: rotation.
* **Unno, Osaki, Ando, Saio & Shibahashi (1989)**, *Nonradial Oscillations of Stars*, 2nd ed.:
  non-radial modes (not obtained). **Christensen-Dalsgaard (2014)**, *Lecture Notes on Stellar
  Oscillations*, 5th ed., Aarhus University: Cowling equations as compared (eqs. 5.12–5.14).
* **Herwig (2000)**, A&A 360, 952–968: exponential overshoot (§3, eq. 2, p. 954).
* **Langer, El Eid & Fricke (1985)**, A&A 145, 179–191: convective mixing as diffusion, $`D = v\ell/3`$.
* **Thoul, Bahcall & Loeb (1994)**, ApJ 421, 828–842: element diffusion.
* **Itoh, Hayashi, Nishikawa & Kohyama (1996)**, ApJS 102, 411–424: neutrino losses.
* **Fowler, Caughlan & Zimmerman (1975)**, ARA&A 13, 69; **Harris, Fowler, Caughlan & Zimmerman
  (1983)**, ARA&A 21, 165: the rate tables behind the implemented fits.
* **Caughlan & Fowler (1988)**, At. Data Nucl. Data Tables 40, 283–334: reaction rates.
* **Angulo et al. (1999)**, Nucl. Phys. A 656, 3–183 (NACRE; not obtained).
* **Wang et al. (2021)**, Chin. Phys. C 45, 030003 (AME2020 tables): reaction Q values.
* **Kunz et al. (2002)**, ApJ 567, 643; **deBoer et al. (2017)**, Rev. Mod. Phys. 89, 035007:
  alternative and current ¹²C(α,γ) evaluations (not read; named as options).
* **MESA source code** (MESAHub/mesa, `main`, read 2026-11-12): `rates/private/ratelib.f90`
  (`rate_c12ag_fxt`, CF88 reproduction) and `net/private/net_screen.f90` (`eval_screen_triple`).
* **McDonald & Zijlstra (2015)**, MNRAS 448, 502–521: calibration of Reimers' $`\eta`$.
* **Iglesias & Rogers (1996)**, ApJ 464, 943–953 (OPAL); **Iglesias & Rogers (1993)**, ApJ 412, 752
  (C/O-enhanced tables); **Ferguson et al. (2005)**, ApJ 623, 585–596: opacity tables.
  **Grevesse & Noels (1993)**, in *Origin and Evolution of the Elements*, CUP: mixture of the
  ready-made tables.
* **Cassisi, Potekhin, Pietrinferni, Catelan & Salaris (2007)**, ApJ 661, 1094–1104;
  **Cassisi, Potekhin, Salaris & Pietrinferni (2021)**, A&A 654, A149; **Potekhin, Pons & Page
  (2015)**, Space Sci. Rev. 191, 239: electron conduction (Ioffe conductivity table).
* **Paxton et al. (2011)**, ApJS 192, 3, §5.1: bounded scale height. **Eggleton (1971)**, MNRAS 151,
  351–364, and **Eggleton (1972)**, MNRAS 156, 361–376: the Eggleton code; neither contains the
  bounded scale height (1972, p. 364, raises the problem without a formula).
* **Acharya et al. (2025)**, Rev. Mod. Phys. 97, 035002 ("Solar Fusion III"): current $`S`$ factors
  for the pp chain and CNO cycle.
* **Ledoux (1951)**, ApJ 114, 373: rotational splitting constant.
* **Pedersen et al. (2018)**, A&A 614, A128: comparison of overshoot prescriptions.
* **Vinyoles et al. (2017)**, ApJ 835, 202: standard solar model calibration.
* **Choi et al. (2016)**, ApJ 823, 102 (MIST); **Bressan et al. (2012)**, MNRAS 427, 127–145
  (PARSEC); **Planck Collaboration (2020)**, A&A 641, A6; **Aver et al. (2021)**, JCAP 2021(03),
  027: primordial helium and $`\Delta Y/\Delta Z`$.
* **Paxton et al. (2011–2019)**, ApJS (the MESA instrument papers; 2011: 192, 3; 2013: 208, 4):
  free-energy EOS blending, energy equation, rotation factors $`f_P`$, $`f_T`$.
* **Christensen-Dalsgaard et al. (1996)**, Science 272, 1286 (Model S): solar radial fundamental mode,
  page 13 of the About text.

* **Trampedach et al. (2014)**, MNRAS 442, 805–820: Hopf-corrected radiative gradient and flux
  relations (§2.5, eqs. 35–38); range of the Krishna Swamy relation (§5).
* **Henyey, Vardya & Bodenheimer (1965)**, ApJ 142, 841: convection with optically thin elements
  (eqs. 34–36, 39–40).
* **Chugunov & DeWitt (2009)**, Phys. Rev. C 80, 014611: screening in ionic mixtures, candidate for
  the strong regime.
* **Potekhin & Chabrier (2012)**, A&A 538, A115: electron screening, $`k_{\mathrm{TF}}`$ and the
  published electron-polarisation interpolation.
* **Kushnir, Waxman & Chugunov (2019)**, MNRAS 486, 449–452: triple-α screening basis.
* **Graboske et al. (1973)**, ApJ 181, 457: screening regimes (not accessed).
* **Reese (2015)**, A&A 578, A37: rotational splitting, eq. 3.
* **Timmes & Arnett (1999)**, ApJS 125, 277–294; the pynucastro EOS documentation: Fermi–Dirac
  electron gas.

### Source comparisons recorded

External review, 2026-11-11, and task 0.3 (same date; originals read unless stated).

| Item | Page | Result | Source location |
|---|---|---|---|
| Fermi–Dirac $`n_e`$, $`P_e`$, $`u_e`$ | 2 | integral definition and normalisation agree, including the $`\vartheta F_{3/2}`$ term; pressure and energy prefactors re-derived by change of variable | pynucastro EOS derivation (which notes a missing factor in the printed Timmes–Arnett expression) |
| EOS consistency residual | 2 | free-energy approach and fixed-composition identity supported | Timmes's thermodynamics reference, eqs. 1–8 |
| Ionisation equilibrium | 2 | chemical-potential relation and classical limit re-derived; signs, mass ratio and ground-state factors consistent | derivation |
| Ionisation energies and ground-state weights (H, He, He⁺) | 2 | agree (13.598, 24.587, 54.418 eV; $`g`$ from ²S₁/₂, ¹S₀); Pols et al. use the same weights with 13.60, 24.58, 54.40 eV | NIST ASD; Pols et al. 1995, eqs. 14–16, p. 966 |
| Pressure-ionisation potential and constants | 2 | agree: form of $`g_{\mathrm{PI}}`$, $`x = n_e m_{\mathrm{H}}`$, 13.60 eV, $`(c_1 \ldots c_4) = (3, 0.25, 2, 0.03)`$, subtraction at $`N_{e0}`$, $`\Delta\mu_e`$. Departure stated: these pages use the exact Fermi–Dirac $`\eta`$ inside $`g_{\mathrm{PI}}`$, Pols et al. the approximation of eq. 29. $`\mu`$, $`P`$, $`u`$ from it are derivations | Pols et al. 1995, §2.2.2, eqs. 18, 26–29, pp. 966–968; Kovetz, Yaron & Prialnik 2009, §3.1, eq. 28 |
| Paczyński interpolation | 2 | $`P_e`$, $`P_d`$, $`f`$ agree. **Corrected:** $`K_1 = 9.91\times10^{12}`$, $`K_2 = 1.231\times10^{15}`$ (they follow from $`m_{\mathrm{H}}`$, not $`m_u`$; the page had the $`m_u`$ values). The paper gives only $`(\partial U_e/\partial\ln T)_\rho`$; the page's $`u_e`$ has that derivative, its $`T`$-independent part is a choice of these pages | Paczyński 1983, eq. 4, p. 317 |
| Opacity-table variable $`R = \rho/T_6^3`$, OPAL grid | 3 | agrees; tables are not rectangular | Iglesias & Rogers 1996, §2.2, p. 944, Table 3; OPAL `GN93hz` header |
| Low-temperature table grid | 3 | agrees; the distributed files have 85 $`T`$ and 155 compositions against 75 and 120 in the paper | Ferguson et al. 2005, §3 p. 587, §6 p. 595 |
| Heavy-element mixture of the tables | 3 | **corrected:** ready-made OPAL (types 1 and 2) and the F05 base set use GN93, not A09 | Iglesias & Rogers 1996, Table 1 p. 944, §5 p. 953; Ferguson et al. 2005, §6 p. 595; OPAL tables page |
| OPAL type 2 (C/O) tables | 3 | agrees; $`X \le 0.7`$, 8 values of $`Z`$, GN93 base | OPAL tables page; Iglesias & Rogers 1993 |
| Radiative + conductive harmonic sum | 3 | agrees; the conduction fit is not chosen | Cassisi et al. 2007, §2, eq. 1 (arXiv version) |
| KW1990 pp, CNO, 3α fits (implemented) | 4 | exponents and polynomials agree; prefactors re-derived to within 0.3%; the book itself not read | Harris et al. 1983, Table 1, p. 168; Fowler, Caughlan & Zimmerman 1975, Table 1, pp. 72, 78 |
| Solar Fusion III $`S_{11}`$, $`S_{33}`$, $`S_{34}`$, $`S_{1,14}`$; $`\langle\sigma v\rangle`$ | 4 | agree. SF III gives $`S(0)`$, $`S'`$, $`S''`$, not thermal rates; $`C_\alpha = 1.45`$–$`1.47\times10^{17}`$ by integration (leading order 1.49). With constant $`S_{1,14}`$ the ¹⁴N(p,γ) rate is 0.52–0.56 of the rate behind the implemented CNO fit (rough estimate) | Acharya et al. 2025, Table I, eqs. 3–4 |
| Solar Fusion III $`S'(0)`$, $`S''(0)`$ (checked 2026-11-12, LaTeX source of arXiv:2405.06470) | 4 | Table I: $`S'_{11} = 4.5\times10^{-24}`$ b, $`S''_{11} = 9.9\times10^{-23}`$ b MeV⁻¹ (§III.F: 11.0 MeV⁻¹ and 242 MeV⁻² times $`S(0)`$); $`S'_{33} = -4.90`$ b, $`S''_{33} = 22.42`$ b MeV⁻¹ (§IV.B fit $`5.21 - 4.90E + 11.21E^2`$ confirms the ½ of eq. 4); $`S'_{34} = -3.03\times10^{-4}`$ b ($`-0.54`$ MeV⁻¹, valid 0–500 keV), $`S''_{34}`$ not given; $`S_{1,14}`$: $`S(0)`$ only. Coded in plan 1.6 | Acharya et al. 2025, Table I, eq. 4, §§III.F, IV.B, V.E, IX |
| ¹²C(α,γ) and 3α $`Q`$ values | 4 | agree | AME2020 (Wang et al. 2021) |
| CF88 ¹²C(α,γ) rate | 4 | transcribed from a reproduction only; gate stays | ORNL CF88 Fortran, reaction 62 |
| Screening: plasma free energy, $`H_{ij}`$, quantum correction $`t_{ij}`$ | 4 | agree | Chugunov & DeWitt 2009, eqs. 12, 24, 25, A4, §§II, IV; Potekhin & Chabrier 2012, eq. 22 |
| Electron-response extension | 4 | its pieces are published, the combination is not (gate stays); a published alternative is recorded | Potekhin & Chabrier 2012, eqs. 14–17, 20–21; Chugunov & DeWitt 2009, eq. A2 |
| Triple-α screening | 4 | basis agrees; the product form is checked in its limits only (derivation), unpublished at intermediate coupling (gate stays) | Kushnir, Waxman & Chugunov 2019, eqs. 1, 7–8 |
| Thermal-neutrino fits | 1 | reference, processes, units (per volume) and ranges recorded; coefficients not yet transcribed | Itoh et al. 1996, §§2–6, eqs. 2.1–5.28, Tables 1–4, 12 |
| Bounded scale height | 5 | formula agrees with MESA. **Attribution corrected:** not in Eggleton 1971 (read in full) | Paxton et al. 2011, §5.1; Eggleton 1971, App. eq. 57; Eggleton 1972, p. 364 |
| Exponential overshoot | 5 | factor 2, $`H_v = fH_P`$, $`D_0`$ near the edge, cut-off $`10^{-2}`$ cm² s⁻¹, $`f \simeq 0.016`$ agree. $`r_0`$ and $`H_{P,0}`$ are the MESA convention; Pedersen et al. use $`H_P`$ at the boundary. **Corrected:** $`f_0 = 0.005`$ withdrawn (no $`f_0`$ in Herwig; MESA default 0; Pedersen et al. 0.001) | Herwig 2000, §3, eq. 2, p. 954; Paxton et al. 2011, §5.2, eq. 2; Pedersen et al. 2018, §§2.1–2.2, §3.1, Table 2; MESA controls (`overshoot_f0`) |
| General buoyancy term $`B`$ | 5, 9 | equivalent; $`B = \delta\,B_{\mathrm{MESA}}`$ | Paxton et al. 2013, §3.3, eqs. 5–8 |
| Mass-coordinate diffusion equation, $`D = v\ell/3`$ | 5 | agree | Langer, El Eid & Fricke 1985, eq. 1 and text before eq. 2, p. 180 |
| Settling advection term | 5 | agrees after $`dm = 4\pi r^2\rho\,dr`$; the fits are for the solar interior only | Thoul, Bahcall & Loeb 1994, eqs. 1–2, 5–6, 19, 54, 56–58 |
| Hopf gradient and flux fractions | 6 | agree, within that construction's assumptions; does not cover the modified mixing-length cubic | Trampedach et al. 2014, eqs. 35–38 |
| Thin-element efficiency | 6 | agrees; $`U_e`$ is a derivation from it | Henyey, Vardya & Bodenheimer 1965, eqs. 32, 34–36, 39–40 |
| Krishna Swamy coefficients and optical-depth scale | 6 | coefficients agree; $`\tau`$ is the Rosseland mean depth (no $`\tau_{5000}`$ conversion); fit to Table 2, $`0 \le \tau \le 10`$ | Krishna Swamy 1966, eq. 33 p. 182, eq. 1 p. 174, p. 175, Table 2 p. 177; Trampedach et al. 2014, §5 |
| de Jager mass loss, linear form | 7 | coefficients agree; the paper's own linear reduction of its 20-term Chebyshev fit; valid for $`\log T_{\mathrm{eff}}`$ 3.3–4.8, $`\log L`$ 2.5–6.7 | de Jager et al. 1988, §4 remark 3 p. 264; eq. 1 pp. 262–263; Table VI |
| Metallicity factor of hot-star winds | 7 | **corrected:** $`(Z/0.02)^{0.7}`$ → $`(Z/0.019)^{0.69}`$ (O stars; 0.64 for B supergiants) | Vink et al. 2001, §7 eqs. 21–22, §4 (arXiv version) |
| Reimers coefficient | 7 | $`4\times10^{-13}`$ confirmed second-hand; original not read; $`\eta \approx 0.5`$ unverified (gate stays) | Reimers 1977, pp. 217, 223 |
| Gravity-darkening exponents | 8 | 0.25 agrees (flux ∝ $`g_{\mathrm{eff}}`$); **corrected:** 0.08 is Lucy's, not von Zeipel's | von Zeipel 1924, eqs. 37–38, pp. 675–676; Lucy 1967, eqs. 3, 6, pp. 90–91 |
| Ledoux constant | 8 | integral agrees | Reese 2015, eq. 3 |
| Cowling equations and boundary conditions | 9 | agree; not compared with Unno et al. (book not obtained) | Christensen-Dalsgaard 2014, eqs. 5.12–5.14 p. 75; 4.39, 4.65, 4.68 |
| $`Y_p`$, $`\Delta Y/\Delta Z`$ | 10 | **corrected** to the MIST pair 0.249 / 1.5 (no source gives 0.248 with 1.5); modern $`Y_p`$ 0.2453–0.2454 listed | Choi et al. 2016, §3.1, eqs. 1–3; Planck 2018 VI, §2.1; Aver et al. 2021, abstract |
| Solar radial fundamental mode (About §13, finding B4) | – | Model S: 258.0 μHz = 64.6 min (not 63); computed with the exact radial wave equation on the model mesh, the same code reproducing the BiSON $`l = 0`$ frequencies for $`n = 6`$–12 to 0.1–0.4 μHz | Christensen-Dalsgaard et al. 1996 (`fgong.l5bi.d.15`); Broomhall et al. 2009, MNRAS 396, L100, Table 1 |

### Items settled by decision (task 0.3, 2026-11-12)

The originals below were not available. Each gate is closed by the decision shown, not by
comparison with the original; the reason is recorded on the page named.

| Item | Page | Decision | Basis |
|---|---|---|---|
| Kippenhahn–Weigert–Weiss (2012) pp and CNO fits | 4 | not adopted; pp and ¹⁴N(p,γ) rates are built from Solar Fusion III by quadrature of its eq. 3 with eq. 4 and Table I, tabulated with C¹ interpolation | Acharya et al. 2025 (verified) |
| $`C_\alpha`$ attributed to Clayton (1968) | 4 | attribution dropped; the constant used is the one computed from Solar Fusion III | Acharya et al. 2025 |
| Caughlan & Fowler (1988) ¹²C(α,γ) fit | 4 | adopted as transcribed, with a multiplier parameter (default 1) | two independent reproductions agree term by term: ORNL CF88 Fortran; MESA `rate_c12ag_fxt` |
| Electron-response extension of the screening | 4 | not adopted; the published Potekhin & Chabrier (2012) eqs. 20–21 form is coded | Potekhin & Chabrier 2012 (verified) |
| Triple-α product form at intermediate coupling | 4 | adopted; exact in the weak and strong limits | derivation; MESA `eval_screen_triple` uses the same product of pair factors |
| Graboske et al. (1973) screening | 4 | not needed (Chugunov & DeWitt chosen) | – |
| Reimers (1975) and $`\eta`$ | 7 | form confirmed second-hand; $`\eta`$ a parameter, default 0.477 | Reimers 1977, p. 217; McDonald & Zijlstra 2015, abstract and eq. 1 |
| OP opacity tables | 3 | dropped from the plan; OPAL is the interior set | OPAL checked and has the type 2 tables |
| Conduction fit | 3 | Potekhin's conductivity table (Cassisi et al. 2007, 2021) | Ioffe page read 2026-11-12; Cassisi et al. 2007 §2, eq. 1 (verified) |
| Cowling and full non-radial equations as printed in Unno et al. | 9 | Christensen-Dalsgaard (2014) is the reference | lecture notes, eqs. 4.39, 4.61–4.68, 5.12–5.14 (verified) |
| Solar size of settling | 5 | now compared: $`\Delta Y_S = -0.029`$ / $`-0.030`$, $`Z_S`$ −9% / −10% | Vinyoles et al. 2017, Table 4 |

Still open, not gated: the modified mixing-length cubic of page 6 is a derivation of these pages and
their reviews, with no published source compared; the coefficients of the Itoh et al. neutrino
fits are to be transcribed from the paper (available) when task 9.6 is coded; the printed CF88
table and the books (KWW 2012, Clayton 1968, Unno et al. 1989) should be compared when they become
available.

None of the reviews looked at the code, re-ran the measured numbers or tested numerical
convergence, except where the table says so (the implemented rate fits).

Journal volumes and first pages in this list were checked on 2026-11-11 against the ADS scans
(papers before 2000) or Crossref and arXiv (later papers); Reimers (1975) and Grevesse & Noels
(1993) are quoted from citing papers.

---
Previous: [Validity, limits and order of work](11-validity-and-roadmap.md) · [Contents](README.md)
