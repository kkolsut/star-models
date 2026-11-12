# 11. Validity, limits and order of work

## What the implemented physics gets right and wrong

**[measured]** with $`X = 0.70`$, $`Z = 0.02`$, $`\alpha_{\mathrm{MLT}} = 1.8`$, 800 points:

| 1 M☉ ZAMS | This model | Detailed ZAMS models (approximate) | Suspected contributors (none isolated) |
|---|---|---|---|
| $`L`$ | 0.4195 L☉ | 0.7 L☉ | interior opacity; pp rate; composition |
| $`R`$ | 1.1146 R☉ | 0.89 R☉ | opacity, EOS without ionisation, atmosphere |
| $`T_{\mathrm{eff}}`$ | 4400 K | 5600 K | follows from $`L`$ and $`R`$ |
| Convective envelope base | 0.968 R | about 0.73 R | envelope opacity, ionisation, atmosphere entropy |
| Central $`\kappa`$ | 2.72 cm² g⁻¹ | roughly half | Kramers formulas |

How to read this table:

* The comparison is with **ZAMS** models. The present-day Sun (1 L☉, 1 R☉, 5772 K, envelope base
  0.713 R, at 4.57 Gyr, after calibration and with settling) is a different benchmark and is used
  only for the solar calibration below.
* The last column lists suspects, not a diagnosis. The only measured fact is the central opacity.
  $`L \propto 1/\kappa`$ is a homology argument for fixed structure; a real model readjusts
  $`T_c`$ and $`\rho_c`$, so "opacity twice too high" does not by itself explain "luminosity 40%
  too low". The contributions can be separated only by experiment: swap one ingredient at a time
  (opacity table, EOS, atmosphere, rates) and record the change. That experiment is part of
  phase 7.
* The depth of the convective envelope depends jointly on the opacity, the EOS, the atmosphere
  (which sets the entropy of the adiabat), the composition, $`\alpha_{\mathrm{MLT}}`$ and age. No
  single upgrade is promised to move it to the solar value.

Qualitatively right: radiative core and convective envelope in Sun-like stars, convective cores and
CNO burning above about 1.2 M☉, deep envelopes at low mass. Massive stars are expected to be closer
to detailed models because their opacity is mostly electron scattering; this has not been checked
against a published grid.

With $`Z = 0.008`$ the model reproduces the textbook StatStar Sun (0.8348 L☉, 5515 K), which checks
the implementation, not the physics.

## Where each part is valid

| Part | Valid | Breaks down |
|---|---|---|
| EOS (ideal gas + radiation) | 0.4–100 M☉ interiors on the main sequence | below 0.4 M☉; envelopes of cool stars; helium cores |
| Opacity (analytic) | hot interiors of massive stars | Sun-like stars (factor 2); below 2600 K |
| Nuclear rates | main-sequence hydrogen burning (pp II share included through $`\psi`$) | pp III share; late burning stages; screening beyond the weak regime (reported) |
| Mixing-length theory | deep, efficient convection, where $`\nabla \to \nabla_{\mathrm{ad}}`$ whatever the theory | superadiabatic surface layers (result depends on the free parameter $`\alpha_{\mathrm{MLT}}`$); boundaries (local theory; core overshoot only as an optional step in the mixing); anything time-dependent |
| Convective boundaries | Sun-like stars | stars with convective cores (overshoot only as an optional, uncalibrated step, off by default); layers with composition gradients (no Ledoux criterion; $`N^2`$ shows them as a diagnostic) |
| Atmosphere | hot stars | cool stars (one $`\kappa`$, grey) |
| Evolution | main sequence | after the main sequence (qualitative) |
| Radial pulsation period | fundamental mode and low overtones of a model in hydrostatic equilibrium, adiabatic, reflecting surface; as accurate as the model's $`R`$ and envelope | modes near the acoustic cut-off ($`\omega/\omega_{\mathrm{ac}}`$ reported, flagged above 0.5); says nothing about driving |

The mass ranges are rough guides. Whether a part of the physics holds is a local question, and the
following diagnostics (to be added as profile outputs with the EOS work) answer it point by point:

| Diagnostic | Definition | Meaning |
|---|---|---|
| Degeneracy | $`(\rho/\mu_e)\,/\,(2.4\times10^{-8}\,T^{3/2})`$ | $`\gtrsim 1`$: electron degeneracy pressure exceeds the ideal electron pressure |
| Coulomb coupling | $`\Gamma_{\mathrm{C}} = \dfrac{\bar Z^2 e^2}{a_i\,k\,T}`$, $`a_i = \left(\dfrac{3}{4\pi n_i}\right)^{1/3}`$ | $`\gtrsim 0.1`$: non-ideal corrections matter; $`\gtrsim 1`$: ideal gas fails |
| Departure from full ionisation | $`1 - x_{\mathrm{H}}`$, $`1 - x_{\mathrm{He}2}`$ from the ionisation equilibrium | not $`\approx 0`$ (partly ionised **or neutral**): $`\mu`$ and $`\rho`$ of the implemented EOS are wrong |
| Partial ionisation | $`x(1-x)`$ for each stage | not $`\approx 0`$: enhanced heat capacity; $`\nabla_{\mathrm{ad}}`$, $`c_P`$, $`\Gamma_1`$ of the implemented EOS are wrong |
| Screening regime | exponent of the screening factor (reported per model and run since plan 1.11) | $`\gtrsim 0.3`$: outside weak screening; at the cap: rates unreliable |
| Radiation pressure | $`1 - \beta`$ | $`\to 1`$: $`\Gamma_1 \to 4/3`$, model close to instability |

The degeneracy coefficient is where $`P_{\mathrm{nd}} = P_{\mathrm{nr}}`$ in the formulas of
[Equation of state](02-equation-of-state.md).

## Size of the effects that are missing

| Effect | Size | Decision |
|---|---|---|
| Tabulated opacity (D14) | central $`\kappa`$ twice the table value; expected to be the largest single error in solar $`L`$, size not isolated | planned, phase 7 |
| Ionisation in the EOS (D15a) | $`\nabla_{\mathrm{ad}}`$ 0.4 → about 0.1 in the ionisation zones; deepens the envelope, by how much depends on opacity and atmosphere | planned, phase 6 |
| Electron degeneracy (D15b) | enables < 0.4 M☉ and red giants | planned, phase 6 |
| Core overshoot (D16) | main sequence 10–25% longer above 1.2 M☉; measured +26% (3 M☉) and +23% (5 M☉) at $`\alpha_{\mathrm{ov}} = 0.2`$ | step overshoot done (phase 2, off by default); diffusive form phase 8 |
| Mass loss (E1) | 1% of $`M`$ at 20 M☉, 10–20% at 60 M☉ on the main sequence | planned, phase 4 |
| Rotation (E3) | a few % in structure, 20–30% in lifetime for fast rotators | planned, phases 5 and 8 |
| Integrated, non-grey atmosphere (D13, F3) | about 100 K in $`T_{\mathrm{eff}}`$ of cool stars | planned, phase 3 |
| pp branching and newer rates (D11) | pp rate 8% low, plus up to 50% from $`\psi`$; measured: ZAMS Sun $`L`$ −0.3%, 1 M☉ main sequence +2.6% | $`\psi`$ done (phase 1, hybrid with the 1990 pp rate); rates from Solar Fusion III with $`S'`$, $`S''`$ still planned |
| Separate $`Q`$ for CNO (D12) | CNO-burning lifetimes 4% too long; measured −4.3% (3 M☉), −4.4% (5 M☉) | done, phase 1 |
| Element settling (F1) | Sun: surface $`Y`$ −0.03, $`Z`$ −10% in 4.6 Gyr | planned, phase 8 |
| Molecular hydrogen (F2) | decides full convection at 0.1–0.5 M☉ | with phase 6 |
| Non-radial modes (E4) | a new output, not a correction | diagnostics ($`N^2`$, $`S_l`$, spacings) done in phase 2; modes planned, phase 9 |
| General relativity (E2) | $`10^{-6}`$ to $`10^{-5}`$ | not planned; redshift output done (phase 1) |
| Magnetic fields, spots (F6) | radii of active M dwarfs several % larger | note only |
| Turbulent pressure | shifts mode frequencies | not planned |
| Thermal neutrinos (F7) | red-giant helium cores from about $`10^8`$ K (plasma neutrinos); all later stages | with the red-giant work, phases 6–7 |
| Total energy conservation check | separates numerical error from physics in evolution runs; measured: Sun to 4.57 Gyr residual $`-4.5\times10^{-4}`$, of which $`10^{-6}`$ numerical | done, phase 1 |
| Helium flash, thermal pulses (F9) | beyond helium ignition | not planned |
| Binary interaction (F11) | most massive stars | out of scope |

Sizes in the "detailed models" column and in this table, other than the measured ones, are
literature values quoted from memory.

## Order of work

Three pieces of infrastructure unlock most of the list:

1. **Opacity tables and a better atmosphere** (D14, D13, F3): everything about cool stars.
2. **An EOS with ionisation, degeneracy and H₂** (D15a, D15b, F2): low-mass stars and red giants.
3. **A composition-diffusion step**: overshoot, semiconvection, rotational mixing and settling
   (D16, D17, E3 stage 3, F1).

| Phase | Content | Pages |
|---|---|---|
| 0 | Regression baseline, text fixes | – |
| 1 | Burning fixes, $`\psi`$, separate $`Q`$, scale-height cap, geometric means, redshift | 4, 5, 7, 9, 10 |
| 2 | Core overshoot; $`N^2`$, $`S_l`$ and asymptotic spacings | 5, 9 |
| 3 | Integrated atmosphere, non-grey $`T(\tau)`$ | 6 |
| 4 | Mass loss | 7 |
| 5 | Rotation: shape (5.1) now; structure (5.2) after phase 7 | 8 |
| 6 | EOS: degeneracy, then ionisation and H₂ | 2 |
| 7 | Opacity tables, conduction, thermal neutrinos (with 6), provisional solar calibration | 3, 1 |
| 8 | Composition-diffusion step and what it unlocks; final solar calibration | 5, 8 |
| 9 | Non-radial modes, pre-main sequence, carbon/oxygen, screening regimes | 9, 7, 4 |

Two points on this order, from the review of these pages on 2026-11-11:

* **Consistency before features.** The total-energy check
  ([Evolution in time](07-evolution.md)) and the EOS consistency test
  ([Equation of state](02-equation-of-state.md)) come before the features they judge: the energy
  check in phase 1 (plan 1.10), the EOS test at the start of phase 6 (plan 6.1), before any new
  EOS is written.
* **Phases 3, 6 and 7 are one piece of work for cool stars.** Atmosphere, EOS and opacity each
  change the envelope, and none can be validated alone. They can be coded in the listed order, but
  results for Sun-like stars are to be judged only when all three are in.
* **Structural rotation comes after that block** (agreed 2026-11-11): task 5.2 runs after the
  calibration at the end of phase 7, then mixing (phase 8). The display-only Roche shape and
  gravity darkening (5.1) stay independent.
* **Source verification is a separate prerequisite**: nothing marked **[check source]** is coded
  before it has been compared with its source (plan rule 7).

Order of execution: 0 → 1 → {2, 4, 5.1} → 3 → 6 → 7 → 5.2 → 8 → 9.

**Solar calibration** (provisional after phases 3, 6 and 7; repeated and final after settling
exists in phase 8). Three parameters are adjusted, the initial helium
$`Y_0`$, the initial metallicity $`Z_0`$ and $`\alpha_{\mathrm{MLT}}`$, so that a 1 M☉ model at
4.57 Gyr matches three observables: $`L_\odot`$, $`R_\odot`$ and the surface $`(Z/X)_\odot`$ of the
adopted abundance mixture. The mixture and the treatment of settling must be stated with the
result, because both change the calibrated values. Without settling (before phase 8) the surface
$`Z/X`$ equals the initial one, and the calibration reduces to $`Y_0`$ and
$`\alpha_{\mathrm{MLT}}`$ against $`L`$ and $`R`$.

Details, files and tests are in [`../implementation-plan.md`](../implementation-plan.md).

---
Previous: [Numerical method and checks](10-numerics-and-checks.md) · [Contents](README.md) · Next: [Constants, symbols and references](12-constants-and-references.md)
