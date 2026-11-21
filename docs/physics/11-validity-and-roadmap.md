# 11. Validity, limits and order of work

## What the implemented physics gets right and wrong

**[measured]** with $`X = 0.70`$, $`Z = 0.02`$, $`\alpha_{\mathrm{MLT}} = 1.8`$, 800 points (columns 2–3):

| 1 M☉ ZAMS | Analytic physics (default) | All stage-4 upgrades | Detailed ZAMS models (approximate) | Measured contributors (table below) |
|---|---|---|---|---|
| $`L`$ | 0.4165 L☉ | 0.6475 L☉ | 0.70 L☉ (Tout et al. 1996 fit: 0.698) | opacity tables (+63% alone); the EOS changes it by −7% to +62% depending on the version |
| $`R`$ | 1.1218 R☉ | 0.9046 R☉ | 0.89 R☉ (fit: 0.888) | EOS with ionisation (−27% alone); tables alone +24% |
| $`T_{\mathrm{eff}}`$ | 4378 K | 5444 K | about 5600 K | follows from $`L`$ and $`R`$ |
| Convective envelope base | 0.967 R | 0.727 R | about 0.73 R | ionisation (0.761 alone); tables and atmosphere move it by 0.01–0.04 |
| Central $`\kappa`$ | 2.74 cm² g⁻¹ | 1.57 cm² g⁻¹ | roughly half the analytic value | tables: 0.57 of the analytic value at the same state |

**Swap experiment (plan 7.7)** **[measured]**, 1 M☉ ZAMS, $`X = 0.70`$, $`Z = 0.02`$,
$`\alpha_{\mathrm{MLT}} = 1.8`$, 800 points; one ingredient at a time, then all
(`examples/physics_assessment.rs`; full output in
[`docs/measurements/physics-assessment-zams.md`](../measurements/physics-assessment-zams.md)). The
nuclear rates have one implementation (Solar Fusion III, plan 1.6) and are not swapped. The last
column is the computing time on the build machine (an emulated CPU, roughly 20× slower than a
desktop, shared with other jobs).

| Physics | $`L`$ [L☉] | $`R`$ [R☉] | $`T_{\mathrm{eff}}`$ [K] | envelope base [r/R] | $`T_c`$ [K] | $`\rho_c`$ | $`\kappa_c`$ | min $`\nabla_{\mathrm{ad}}`$, $`q > 0.5`$ | $`P_{\mathrm{phot}}`$ [dyn cm⁻²] | time [s] |
|---|---|---|---|---|---|---|---|---|---|---|
| analytic (default) | 0.4165 | 1.1218 | 4378 | 0.9668 | 1.2515e7 | 62.13 | 2.737 | 0.398 | 1.76e5 | 8 |
| atmosphere: integrated grey | 0.4165 | 1.1099 | 4402 | 0.9573 | 1.2515e7 | 62.13 | 2.737 | 0.398 | 3.61e5 | 5 |
| atmosphere: Krishna Swamy | 0.4165 | 1.1209 | 4380 | 0.9663 | 1.2515e7 | 62.13 | 2.737 | 0.398 | 2.31e5 | 7 |
| EOS: Fermi–Dirac, fully ionised | 0.3863 | 1.1082 | 4323 | 0.9653 | 1.2330e7 | 61.21 | 2.834 | 0.398 | 1.92e5 | 170 |
| EOS: + ionisation, no pressure ionisation | 0.6765 | 0.9462 | 5381 | 0.7934 | 1.3784e7 | 68.19 | 2.190 | 0.115 (log T = 4.11) | 5.42e4 | 171 |
| EOS: + pressure ionisation + H₂ (full) | 0.3942 | 0.8134 | 5071 | 0.7610 | 1.2369e7 | 61.33 | 2.811 | 0.128 (log T = 4.14) | 9.26e4 | 420 |
| opacity: tables | 0.6790 | 1.3914 | 4442 | 0.9567 | 1.3299e7 | 81.60 | 1.541 | 0.394 | 7.98e4 | 57 |
| all upgrades | 0.6475 | 0.9046 | 5444 | 0.7266 | 1.3144e7 | 80.72 | 1.569 | 0.120 (log T = 4.12) | 7.25e4 | 1040 |

How to read this table:

* The comparison is with **ZAMS** models. The present-day Sun (1 L☉, 1 R☉, 5772 K, envelope base
  0.713 R, at 4.57 Gyr, after calibration and with settling) is a different benchmark and is used
  only for the solar calibration below.
* **Luminosity: the interior opacity, as suspected.** The tables alone raise $`L`$ from 0.417 to
  0.679 L☉ (central κ 2.74 → 1.54, $`T_c`$ +6%, $`\rho_c`$ +31%). No other single ingredient moves
  $`L`$ the right way: Fermi–Dirac electrons lower it by 7% (degeneracy pressure at $`\eta_c = -1.5`$),
  and the atmosphere does not change it.
* **Radius and envelope: the EOS.** Ionisation lowers $`\nabla_{\mathrm{ad}}`$ to 0.12–0.13 near 10⁴ K
  (from 0.40) and deepens the envelope from 0.967 to 0.761 R; with the analytic opacity this
  shrinks $`R`$ to 0.81 R☉. The tables alone *raise* $`R`$ (2–3× more opacity at 10⁵–10⁶ K, page 3),
  so partial upgrades do not move the model monotonically towards the Sun.
* **Ionisation without pressure ionisation is not usable:** it recombines 26% of the central
  hydrogen, which raises $`\mu_c`$ and gives a spuriously bright, compact star (0.677 L☉, 0.946 R☉).
  With the Pols term the centre is fully ionised again ($`L`$ back to 0.394 with the analytic opacity).
* **The atmosphere** changes $`R`$ by 1% (grey integrated) or 0.1% (Krishna Swamy) and the envelope
  base by up to 0.01 R; its effect on the photospheric pressure is large (factor 2), on the
  interior small, as expected for the ZAMS Sun.
* **All together** the 1 M☉ ZAMS agrees with detailed models to 7% in $`L`$, 2% in $`R`$ and 0.003 R
  in the envelope base, without any calibration ($`X = 0.70`$, $`Z = 0.02`$, $`\alpha = 1.8`$).

**Defaults (decided 2026-11-17, plan 7.7).** The analytic physics stays the default of the app and
the CLI; the full set is offered as the preset "All upgrades" (`PhysicsOptions::all_upgrades()`),
calibrated below. Reasons: (1) cost: the full set takes about 130 times longer per model (a ZAMS
model takes tens of seconds to a minute on a desktop, an evolution run tens of minutes), too slow
for the app's auto-solve on every slider change; (2) the partial sets are not better than the
analytic one (see above), so there is no cheaper intermediate default; (3) the regression baseline
and every measurement quoted for the analytic physics stay valid. The decision is to be revisited
when the Fermi–Dirac EOS is tabulated (open item). The regression tests keep pinning the analytic
physics; the full set has its own reference test (`tests/cool_star.rs`, ignored by default for its
cost).

* The depth of the convective envelope depends jointly on the opacity, the EOS, the atmosphere
  (which sets the entropy of the adiabat), the composition, $`\alpha_{\mathrm{MLT}}`$ and age; the
  table separates the first three.

Qualitatively right: radiative core and convective envelope in Sun-like stars, convective cores and
CNO burning above about 1.2 M☉, deep envelopes at low mass. Massive stars are closer to detailed
models because their opacity is mostly electron scattering: within 10% in $`L`$ and $`R`$ for 5–20 M☉
against the published grid (plan 7.9, below).

With $`Z = 0.008`$ the model reproduces the textbook StatStar Sun (0.8348 L☉, 5515 K), which checks
the implementation, not the physics.

## Where each part is valid

| Part | Valid | Breaks down |
|---|---|---|
| EOS (ideal gas + radiation, default) | 0.4–100 M☉ interiors on the main sequence | below 0.4 M☉; envelopes of cool stars; helium cores |
| EOS, full (Fermi–Dirac, ionisation, pressure ionisation, H₂; selectable) | cool envelopes, low-mass stars, partly degenerate cores | Coulomb corrections (≈1% in the Sun) missing; metal ionisation crude (pool); pairs above 10⁹ K; stops runs at degenerate cores (stage 5) |
| Opacity (analytic, default) | hot interiors of massive stars | Sun-like stars (factor 2, plan 7.3 table on page 3); below 2600 K |
| Opacity tables (OPAL + Ferguson, GN93, conduction; selectable) | 2.7 ≤ log T ≤ 8.7, scaled-solar compositions | C/O-rich matter after helium burning (clamped, qualitative); outside the tables (analytic blend) |
| Nuclear rates | main-sequence hydrogen burning (Solar Fusion III rates; pp II share included through $`\psi`$) | pp III share; late burning stages; screening beyond the weak regime (reported) |
| Mixing-length theory | deep, efficient convection, where $`\nabla \to \nabla_{\mathrm{ad}}`$ whatever the theory | superadiabatic surface layers (result depends on the free parameter $`\alpha_{\mathrm{MLT}}`$); boundaries (local theory; core overshoot only as an optional step in the mixing); anything time-dependent |
| Convective boundaries | Sun-like stars | stars with convective cores (overshoot only as an optional, uncalibrated step, off by default); layers with composition gradients (no Ledoux criterion; $`N^2`$ shows them as a diagnostic) |
| Atmosphere (grey, default) | hot stars | cool stars (one $`\kappa`$, grey); the integrated Krishna Swamy atmosphere (selectable) is a solar fit, an extrapolation for other stars |
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
| Screening regime | exponent of the screening factor (reported per model and run since plan 1.11) | $`\gtrsim 0.3`$: outside weak screening (all-regime prescription since plan 9.7; helium ignition not validated) |
| Radiation pressure | $`1 - \beta`$ | $`\to 1`$: $`\Gamma_1 \to 4/3`$, model close to instability |

The degeneracy coefficient is where $`P_{\mathrm{nd}} = P_{\mathrm{nr}}`$ in the formulas of
[Equation of state](02-equation-of-state.md).

## Size of the effects that are missing

| Effect | Size | Decision |
|---|---|---|
| Tabulated opacity (D14) | measured: central $`\kappa`$ 0.57 of the analytic value; alone raises the ZAMS Sun's $`L`$ by 63% | done, stage 4 (selectable) |
| Ionisation in the EOS (D15a) | measured: $`\nabla_{\mathrm{ad}}`$ 0.40 → 0.12 near 10⁴ K; envelope base 0.967 → 0.761 R alone, 0.727 R with tables and atmosphere | done, stage 4 (selectable) |
| Electron degeneracy (D15b) | measured: ZAMS Sun $`\eta_c = -1.5`$, $`L`$ −7%; enables < 0.4 M☉ and red giants | done, stage 4 (selectable); red-giant cores need stage 5 |
| Core overshoot (D16) | main sequence 10–25% longer above 1.2 M☉; measured +26% (3 M☉) and +23% (5 M☉) at $`\alpha_{\mathrm{ov}} = 0.2`$ | step overshoot done (phase 2, off by default); diffusive form phase 8 |
| Mass loss (E1) | 1% of $`M`$ at 20 M☉, 10–20% at 60 M☉ on the main sequence; measured with the standard rates: 60 M☉ loses 14.3% on the main sequence, 20 M☉ 2.4%, 5 M☉ 0.10% | done, phase 4 (off by default) |
| Rotation (E3) | a few % in structure, 20–30% in lifetime for fast rotators | shape and gravity darkening for display done (5.1); structure planned, phases 5 and 8 |
| Integrated, non-grey atmosphere (D13, F3) | measured on the analytic ZAMS Sun: $`R`$ −1% (grey integrated), −0.1% (Krishna Swamy); $`P_{\mathrm{phot}}`$ up to 2× | done, stage 4 (selectable) |
| pp branching and newer rates (D11) | pp rate 5–7% low, CNO rate about 2× high, plus up to 50% from $`\psi`$; measured with the per-channel $`Q`$ and the other stage-1 changes: ZAMS Sun $`L`$ −1.0%, 1 M☉ main sequence 14.03 → 15.12 Gyr, 3 M☉ 285.9 → 268.6 Myr | done, phase 1: $`\psi`$, pp and ¹⁴N(p,γ) rates from Solar Fusion III |
| Separate $`Q`$ for CNO (D12) | CNO-burning lifetimes 4% too long; measured −4.3% (3 M☉), −4.4% (5 M☉) | done, phase 1 |
| Element settling (F1) | Sun: surface $`Y`$ −0.03, $`Z`$ −10% in 4.6 Gyr | planned, phase 8 |
| Molecular hydrogen (F2) | decides full convection at 0.1–0.5 M☉ | done, stage 4 (in the full EOS) |
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

**Provisional calibration (plan 7.8, done 2026-11-17)** **[measured]**, full physics set
(`PhysicsOptions::all_upgrades()`: Fermi–Dirac EOS with ionisation, pressure ionisation and H₂;
OPAL + Ferguson tables with conduction; Krishna Swamy atmosphere), mixture Grevesse & Noels (1993)
with $`(Z/X) = 0.0245`$, no diffusion, 400 points and coarse time steps (targets
$`\Delta\log = 0.05`$, $`\Delta X = 0.05`$; V.2 puts the mesh error at 400 points at 2×10⁻⁴ in
$`L`$) — **provisional** (`examples/solar_calibration.rs`, Newton with a difference Jacobian and
a Broyden update; log in [`docs/measurements/solar-calibration.md`](../measurements/solar-calibration.md)):

| $`Y_0`$ | $`X_0`$ | $`Z_0`$ | $`\alpha_{\mathrm{MLT}}`$ | $`L`$, $`R`$ at 4.57 Gyr | $`T_{\mathrm{eff}}`$ | envelope base | $`T_c`$ | $`X_c`$ | ZAMS $`L`$ |
|---|---|---|---|---|---|---|---|---|---|
| 0.2834 | 0.6994 | 0.01714 | 2.020 | 0.9999, 0.9999 | 5772 K | 0.733 R | 1.547×10⁷ K | 0.345 | 0.727 L☉ |

For comparison, the B16-GS98 standard solar model (Vinyoles et al. 2017, Table 4; with diffusion,
another code and MLT normalisation) has $`Y_{\mathrm{ini}} = 0.2718`$, $`Z_{\mathrm{ini}} = 0.0187`$,
$`\alpha_{\mathrm{MLT}} = 2.18`$, envelope base 0.7116 R (helioseismic 0.713), central
$`X_c = 1 - Y_C - Z_C = 0.347`$. Our $`Y_0`$ is 0.012 higher and the envelope 0.02 R shallower,
both in the direction expected without settling (settling lowers the surface $`Z/X`$ and deepens
the envelope; it is plan 8.5). The parameters are the preset `StarParams::sun_calibrated()` and
the app's *Calibrated Sun* button; `tests/cool_star.rs::calibrated_sun` (ignored, slow) re-checks
them. The final three-parameter calibration is plan 8.6.

**Validation against published models (plan 7.9)** **[measured]**. Different physics, unlike
track V. Reference set: the Eggleton-code models of Pols et al. (1998), represented by the ZAMS fits
of Tout et al. (1996; quoted fit errors at $`Z = 0.02`$ below 3% in $`L`$ and 1.2% in $`R`$) and the
main-sequence lifetime fits of Hurley, Pols & Tout (2000); the solar reference is the B16-GS98
standard solar model above. Our models: all upgrades, $`X = 0.70`$, $`Z = 0.02`$,
$`\alpha_{\mathrm{MLT}} = 2.02`$ (the calibrated value; Pols et al. used $`X = 0.70`$, $`Z = 0.02`$,
$`\alpha = 2.0`$), 800 points (`examples/validation.rs`; output in
[`docs/measurements/validation.md`](../measurements/validation.md)). Tolerances stated in advance:
5% in $`L`$, 3% in $`R`$ (the fit errors plus a margin for different input physics).

| $`M`$ [M☉] | $`L`$ [L☉] | fit | $`\Delta L`$ | $`R`$ [R☉] | fit | $`\Delta R`$ | $`T_{\mathrm{eff}}`$ [K] | $`T_c`$ [K] | $`\rho_c`$ | env. base | core $`m/M`$ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.5 | 0.0398 | 0.0379 | **+5.1%** | 0.4417 | 0.4588 | **−3.7%** | 3880 | 8.47e6 | 69.4 | 0.548 | – |
| 0.8 | 0.2241 | 0.2363 | **−5.2%** | 0.7180 | 0.7292 | −1.5% | 4686 | 1.10e7 | 76.1 | 0.679 | – |
| 1 | 0.6511 | 0.6977 | **−6.7%** | 0.8872 | 0.8882 | −0.1% | 5505 | 1.32e7 | 80.7 | 0.721 | – |
| 1.5 | 4.473 | 4.628 | −3.4% | 1.4494 | 1.4372 | +0.8% | 6972 | 1.81e7 | 80.0 | 0.969 | 0.071 |
| 2 | 14.98 | 15.89 | **−5.7%** | 1.6403 | 1.6102 | +1.9% | 8867 | 2.07e7 | 62.2 | 0.999 | 0.144 |
| 3 | 75.00 | 79.98 | **−6.2%** | 2.0322 | 1.9791 | +2.7% | 11915 | 2.35e7 | 38.1 | – | 0.199 |
| 5 | 512.0 | 530.1 | −3.4% | 2.7052 | 2.6356 | +2.6% | 16693 | 2.67e7 | 19.8 | – | 0.246 |
| 10 | 5505 | 5552 | −0.8% | 3.9951 | 3.9419 | +1.4% | 24875 | 3.10e7 | 8.71 | – | 0.317 |
| 20 | 42604 | 42845 | −0.6% | 5.8838 | 5.9990 | −1.9% | 34187 | 3.50e7 | 4.47 | – | 0.446 |

Differences beyond the tolerances (bold) and their named suspects (open items, none isolated):

* **$`L`$ 5–7% low for 0.8–3 M☉.** Suspects: the nuclear rates (Solar Fusion III against the
  Caughlan & Fowler 1988 rates of the Eggleton code; our CNO rate is about half the old one, which
  matters at 1.5–3 M☉), and the opacity mixture (GN93 here; Pols et al. used OPAL tables of their
  time). The 1 M☉ value is the same 7% found against the detailed ZAMS in the table at the top.
* **0.5 M☉: $`L`$ +5%, $`R`$ −4%.** Suspects: the EOS (no Coulomb corrections here; the Eggleton EOS
  has a different pressure-ionisation treatment) and the Krishna Swamy atmosphere used outside the
  solar regime.
* The Tout et al. fits interpolate a grid with their own errors (up to 3% / 1.2%), so differences
  of that size are not attributable.

*Main-sequence lifetimes* (to $`X_c < 10^{-3}`$; all upgrades with the time-step targets of the
calibration, $`\Delta\log = 0.05`$, $`\Delta X = 0.05`$; analytic values from the regression runs):

| $`M`$ [M☉] | all upgrades | analytic | Hurley et al. (2000) fit | $`\Delta`$ all | $`\Delta`$ analytic |
|---|---|---|---|---|---|
| 1 | 10.16 Gyr | 15.12 Gyr | 11.00 Gyr | **−7.7%** | +37% |
| 3 | 312.9 Myr | 268.6 Myr | 377.5 Myr | **−17.1%** | −29% |
| 5 | 84.6 Myr | 72.0 Myr | 104.0 Myr | **−18.7%** | −31% |

Tolerance 5% (the fit's $`t_{\mathrm{MS}} = \max(0.95, t_{\mathrm{hook}}/t_{\mathrm{BGB}})\,t_{\mathrm{BGB}}`$ is
itself a definition with a few per cent of slack). Suspects: **core overshoot** for 3 and 5 M☉ — the
Pols et al. (1998) models include convective overshooting, our default has none; the step overshoot
of plan 2.1 at $`\alpha_{\mathrm{ov}} = 0.2`$ lengthens our analytic 3 and 5 M☉ main sequence by 26%
and 23% (page 5), the size of the gap (to be measured with the full set; open item). For 1 M☉ the
fit's end point (95% of the time to the base of the giant branch) and our $`X_c < 10^{-3}`$ differ in
definition; the remaining 8% is an open item (the ZAMS $`L`$ is 7% *lower* than the fit, which
would lengthen, not shorten, the main sequence).

*The analytic default against the same fits* ($`\alpha = 1.8`$): $`\Delta L`$ = −74%, −47%, −40%, −26%,
−18%, −6%, +7%, +9%, +4% and $`\Delta R`$ = +53%, +40%, +26%, −14%, −14%, −11%, −8%, −6%, −7% for 0.5,
0.8, 1, 1.5, 2, 3, 5, 10, 20 M☉. This checks the statement at the top of the page that massive
stars are closer to detailed models: for the analytic physics they are within 10% (against 40–70%
at and below 1 M☉), not better; with all upgrades they agree to 1–2%.

*Validity ranges following from this table* (all upgrades): ZAMS $`L`$ to about 7% and $`R`$ to
about 4% for 0.5–20 M☉; main-sequence lifetimes 8% (1 M☉) to 19% (5 M☉) shorter than the grid,
which includes overshoot; the Sun after calibration to $`10^{-4}`$ in $`L`$, $`R`$ by construction, its
envelope base to 0.02 R (no settling). These replace the "qualitative" label for the full set on the
main sequence; the analytic default stays qualitative.

Details, files and tests are in [`../implementation-plan.md`](../implementation-plan.md).

---
Previous: [Numerical method and checks](10-numerics-and-checks.md) · [Contents](README.md) · Next: [Constants, symbols and references](12-constants-and-references.md)
