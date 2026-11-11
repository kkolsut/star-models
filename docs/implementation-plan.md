# Implementation plan for the physics review

This plan turns the proposals of [`physics-review.md`](physics-review.md) into ordered work. IDs
(B1, C7, D14, E3, F5 …) and equation numbers (§8.x) refer to that file. Nothing here has been
started.

> **Revised 2026-11-11** after a review of [`physics/`](physics/README.md). That folder is now the
> current statement of the equations; where it differs from `physics-review.md` §8, use it. Changes
> to this plan: structural rotation (5.2) moves behind phases 3, 6 and 7; two consistency checks
> are added (1.10, 6.1); the tasks for the atmosphere, mass loss, degeneracy, pp branching and
> neutrinos are corrected; source verification is a gate (rule 7). Task numbers are unchanged.

Effort: **S** = under a day, **M** = a few days, **L** = a week or more. These are rough estimates
for someone who knows the code.

## Overview

| Phase | Content | Changes default results? | Effort | Depends on |
|---|---|---|---|---|
| 0 | Baseline tests, text fixes | no | S | – |
| 1 | Small corrections and cheap outputs | slightly | M | 0 |
| 2 | Convective boundaries, seismic diagnostics | only with overshoot on | M | 1 |
| 3 | Atmosphere | yes, cool stars | M | 1 |
| 4 | Mass loss | only when enabled | M | 1 |
| 5.1 | Rotation: appearance only (Roche shape, gravity darkening) | display only | S–M | 1 |
| 6 | Equation of state | yes | L | 3 |
| 7 | Opacity tables, joint assessment, solar calibration | yes, strongly | L | 3, 6 |
| 5.2 | Rotation: structure | only when enabled | M | 7 (calibrated), 4 helps |
| 8 | Composition-diffusion step and what it unlocks | only when enabled | L | 2; 5.2 for rotational mixing |
| 9 | New modes and phases (non-radial solver, pre-main sequence, C/O) | no | L | 2, 6, 7 |

**Order of execution:** 0 → 1 → {2, 4, 5.1 in any order} → 3 → 6 → 7 → 5.2 → 8 → 9.

* Phases 3, 6 and 7 (atmosphere, EOS, opacity) are one block for cool stars. Each changes the
  envelope and none can be validated alone: results for Sun-like stars are judged, and the baseline
  for them re-pinned with a physical meaning, only at the end of phase 7.
* Structural rotation (5.2) waits for that block. Its corrections are a few percent, smaller than
  the present errors of the closures it would sit on. The display-only 5.1 is independent.
* Phases 2, 4 and 5.1 are independent of each other after phase 1.

## Rules for every task

1. **Physics lives in `stellar-core`.** New microphysics goes behind the existing traits (`Eos`,
   `Opacity`, `Nuclear` in `physics/mod.rs`). The core must keep building without file or thread
   access, because of the planned web build ([`web-version.md`](web-version.md)): tables are
   embedded with `include_bytes!`, not read at run time.
2. **New inputs are optional and default to the present behaviour.** A new field of `StarParams`
   needs, in this order:
   * the field with `#[serde(default)]` and a line in `validate` and `interpolate`
     (`structure/model.rs`);
   * a CLI flag (`crates/stellar-cli/src/main.rs`);
   * a slider or checkbox (`crates/stellar-app/src/app.rs`), `ZamsKey` in `worker.rs` if the ZAMS line
     depends on it, and a row in `views/compare.rs`;
   * the parameter tables in `views/about.rs` and `docs/about/03-parameters.md`.
3. **The Jacobian is built by finite differences with a step of $`10^{-6}`$.** Any nested iteration
   inside the physics (root search for $`\rho`$, Saha, atmosphere integration) must converge to
   about $`10^{-12}`$ relative, and any interpolation must have continuous first derivatives.
   Otherwise Newton stalls.
4. **Documentation moves with the code.** For each change:
   * edit the text in `crates/stellar-app/src/views/about.rs` and the matching page in `docs/about/`;
   * add or change equations in `crates/stellar-app/equations/equations.json`, then run
     `cd tools/equations && npm install && npm run render` (regenerates the SVGs and
     `src/equations.rs`);
   * `cargo test -p stellar-app` checks that every equation is present in `docs/about/`.
5. **Quoted numbers are re-measured** whenever default results change (list in task 1.9).
6. **Before merging:** `cargo fmt --all --check`, `cargo clippy --workspace --all-targets -- -D warnings`,
   `cargo test --workspace` (the CI in `.github/workflows/ci.yml` runs the same).
7. **Source verification is a gate.** No formula or coefficient marked **[check source]** in
   `physics/` (or [source] in `physics-review.md` §8) is coded until it has been compared with the
   cited source. The verification is its own step, done before the task starts: record the
   reference with equation or page number next to the formula in `physics/`, remove the marker, and
   repeat the reference in the code comment. If the source disagrees, fix `physics/` first.

## Phase 0: baseline and text

### 0.1 Baseline regression test (S)

* **Mostly done on the unmerged branch `quick-wins`:** `crates/stellar-core/tests/regression.rs` pins
  ZAMS and evolved values to $`10^{-4}`$. The task is therefore: merge `quick-wins` into `main`, then
  add 5 and 20 M☉ if they are not covered. The description below is what the test must provide.

* **Why:** later phases must change results on purpose only.
* **Do:** new `crates/stellar-core/tests/baseline.rs` that pins, to 0.1%, the present values for
  1, 5 and 20 M☉ ZAMS models ($`L`$, $`R`$, $`T_{\mathrm{eff}}`$, $`T_c`$, core mass, envelope
  base), the Sun at 4.57 Gyr, and the 5 M☉ main-sequence lifetime. Current values are in
  `physics-review.md` §2 and §7.
* **Done when:** the test passes on the unchanged code. Each later task that changes defaults
  updates these numbers in the same commit, with the reason in the commit message.

### 0.2 Text fixes B1–B6 (S)

* **Files:** `views/about.rs`, `docs/about/04`, `06`, `07`, `08`, `10`, `12`, `13`, `18`,
  `equations/equations.json` (ids `nuc_pp` and `nuc_psi`), `physics/nuclear.rs` (module
  comment and the name `psi_pp`).
* **Do:** rename $`\psi`$ to $`g_{11}`$ and state $`\psi = 1`$; cite Kippenhahn & Weigert for the
  rates; update the assumptions of section 6 and the Ledoux remark of section 8; reword the solar
  period; document the 40-iteration limit, the burning-limit floor and the use of $`Q_{\mathrm{H}}`$
  for CNO; qualify "exactly 1" for the virial check.
* **Done when:** equations re-rendered, `cargo test -p stellar-app` passes, no result changes.

## Phase 1: small corrections and cheap outputs

| Task | ID | Change | Files | Test | Effort |
|---|---|---|---|---|---|
| 1.1 | C7 | Burn mixed zones with the zone-averaged rate, then clip (§8.1) | `structure/evolve.rs` (`advance_composition`, reuse `hydrogen_rates`, `convective_runs`) | extend `tests/evolution.rs::hydrogen_burnt_matches_energy_radiated` to 3 M☉ to the end of the main sequence; ratio within 0.5% | S |
| 1.2 | C8 | Helium time-step limit (§8.1) | `evolve.rs` (`burn_limit_s`) | unit test: a model with $`\varepsilon_{3\alpha} > 0`$ gets a finite limit | S |
| 1.3 | C10 | Geometric means in the pulsation solver (§8.11) | `pulsation.rs` (`Lawe::interval`) | uniform-sphere test still passes; solar period changes less between 800 and 1600 points than before | S |
| 1.4 | C9 | Bounded scale height in the mixing length (§8.2) | `physics/convection.rs` (`mlt`) | new unit test at $`g \to 0`$; ZAMS $`L`$, $`R`$ unchanged to $`10^{-4}`$ | S |
| 1.5 | D12 | Separate $`Q_{pp}`$ and $`Q_{\mathrm{CNO}}`$ (§8.1) | `evolve.rs` (constants, `hydrogen_rates`, `advance_composition`), budget test | budget test with the per-channel energies | S |
| 1.6 | D11 | Branching factor $`\psi`$ and newer rate fits (§8.1) | `physics/nuclear.rs` | unit tests: $`\phi \to 1`$ at low $`T`$, $`\phi \le 2`$; solar-centre pp II share 10–40%; existing order-of-magnitude tests | S–M |
| 1.7 | F5 | Optional "helium follows metallicity" mode (§8.11) | `app.rs` (UI only; sets $`X`$ from $`Z`$), CLI flag | UI check; no core change | S |
| 1.8 | E2 | Gravitational redshift in the results (§8.8); add relativity to "Not modelled" | `structure/model.rs` (`Globals`, `build_model`), `output.rs`, app results panel, CLI print | Sun at 1 R☉ gives 636 m/s | S |
| 1.9 | – | Re-measure quoted numbers | `docs/about/12`, `13`, `15`, `views/about.rs`, `tests/baseline.rs` | – | S |
| 1.10 | – | Total-energy check for evolution runs (`physics/07`): accumulate $`\int L\,dt`$, $`\int\!\!\int\varepsilon_{\mathrm{nuc}}\,dm\,dt`$ and $`\Delta(E_{\mathrm{int}} + E_{\mathrm{grav}})`$; report the residual. Rename the existing check "fuel accounting" | `structure/evolve.rs`, `structure/model.rs` (`Globals`: $`E_{\mathrm{int}}`$, $`E_{\mathrm{grav}}`$), `output.rs`, results panel | residual below 1% for 1 and 5 M☉ on the main sequence; on halving the time step the residual converges to a limit, which need not be zero: what remains is the physical residual of the omitted composition term in $`\varepsilon_{\mathrm{grav}}`$, reported separately from the part that shrinks | S–M |

Notes:

* **1.10 comes early on purpose.** Every later phase is judged with it. It also measures the
  neglected composition term of $`\varepsilon_{\mathrm{grav}}`$ (D17).
* **1.6 changes hydrogen consumption too.** With $`\psi \ne 1`$ the burning step must use
  $`Q_{pp}^{\mathrm{eff}} = Q_{pp,\mathrm{I}}\,\psi/\phi`$ (`physics/04`); `Nuclear` therefore has
  to return $`\phi`$ or the hydrogen-consumption rate, not only $`\varepsilon_{pp}`$. Take the
  $`S`$ factors from Solar Fusion III; the KWW 2012 fits are an interim step.

* **1.6 needs the coefficients checked first.** §8.1 marks the newer fits and $`C_\alpha`$ as quoted
  from memory. Compute $`\alpha`$ from the definition with consistently defined
  $`\langle\sigma v\rangle`$ (`physics/04`); the leading-order
  $`C_\alpha \simeq 1.49\times10^{17}`$ from Solar Fusion III is the check value, not the
  implementation. Take the three $`S`$ factors from that one compilation and
  record the source in the code comment.
* **1.6 and 1.5 change default results.** Expect a slightly lower $`T_c`$ and a few % shorter
  lifetimes for CNO-burning stars. $`L`$ of the Sun is expected to stay near 0.42 L☉ (the opacity is the leading suspect for that value,
  not an established cause; see `physics/11`).
* **1.9 list:** 0.42 L☉, 1.11 R☉, 4400 K, 0.97 R, 17 Gyr, 1.43 Gyr (1.8 M☉), 76 Myr (5 M☉), 0.99 h,
  budget ratio 0.995, "26 steps", "about 1650 steps".

## Phase 2: convective boundaries and seismic diagnostics

### 2.1 Core overshoot, D16 (S–M)

* **Input:** `alpha_ov` in `StarParams`, default 0, range 0–0.5 (rule 2).
* **Do:** in `evolve.rs`, after `convective_runs`, extend a run that starts at the centre outward to
  $`r_{\mathrm{ov}}`$ (§8.2) before averaging rates and composition. The structure equations are
  untouched.
* **Outputs:** report the mixed-core mass next to the convective-core mass.
* **Tests:** `alpha_ov = 0` reproduces the baseline; with 0.2 the 3 M☉ main-sequence lifetime is
  longer and the helium core at the end is larger; a 1.1 M☉ star without a convective core is
  unaffected.
* **Docs:** sections 3, 8 and 12; move overshoot out of "Not included".

### 2.2 Buoyancy and Lamb frequencies, asymptotic spacings, E4 stage 1 (S–M)

* **Do:** add `n2` (formula of `physics/09`, with $`\nabla_\mu`$ from the stored composition; from
  phase 6 on with the general composition term $`B`$ of `physics/05`, which needs
  $`(\partial\ln\rho/\partial X_i)_{P,T}`$ from the EOS) and `cs` to `Profile`
  (`structure/model.rs`); a test that compares it with the definition
  $`g\,(\Gamma_1^{-1}\,d\ln P/dr - d\ln\rho/dr)`$ in the envelope; one $`\Delta\Pi_1`$ per connected
  radiative region, labelled as a diagnostic; add
  $`\Delta\nu_{\mathrm{asym}}`$ and $`\Delta\Pi_1`$ to `Globals`; add both columns to
  `output.rs`.
* **App:** a propagation diagram ($`N`$, $`S_1`$, $`S_2`$ against the chosen x axis) as a new
  profile panel; the two spacings in the results panel next to the scaling values.
* **Tests:** in homogeneous (ZAMS) models, $`N^2 < 0`$ where `convective` is set and $`> 0`$
  elsewhere, with a tolerance around zero for nearly adiabatic layers; in evolved models the
  formula is tested against its differential definition instead, and Schwarzschild and Ledoux
  instability are reported separately; for the ZAMS Sun
  $`\Delta\nu_{\mathrm{asym}}`$ is within 10% of the scaling value; $`\Delta\Pi_1`$ is finite only
  when a radiative region exists.
* **Docs:** sections 4, 5, 13.

## Phase 3: atmosphere (D13, F3)

* **Do (M):**
  1. `physics/atmosphere.rs`: a function integrating §8.4 in $`\ln\tau`$ (fixed-step RK4, about
     60 steps) and returning $`P(\tau_{\mathrm{fit}})`$; a `TTau` enum with `Eddington` and
     `KrishnaSwamy`.
  2. `structure/equations.rs` (`outer_bc`): replace `photosphere_pressure` by the integrated value;
     with a general $`q(\tau)`$ use the three surface relations of §8.4. $`T_{\mathrm{eff}}`$ is then
     derived from $`L`$ and $`R`$, not read from the last mesh point (`build_model` uses
     `pr.t[n - 1]` today, for `teff` and `nu_max_uhz`).
  3. `structure/initial.rs`: `solve_photosphere` uses the same function.
* **Corrections from `physics/06`:**
  * start from $`P(\tau_0) = P(0) + g\tau_0/\kappa_0`$ with $`P(0) = aT^4(0)/3`$, not from
    $`aT^4(\tau_0)/3`$; check that the result does not depend on $`\tau_0`$. Document this as
    the grey/LTE radiation-pressure approximation; the constant-$`\kappa`$ test verifies the
    Eddington limit only;
  * integrate the pair $`(P, T)`$ with the Hopf-corrected gradient
    $`\nabla_{\mathrm{rad,atm}} = A\,\nabla_{\mathrm{rad}}`$, $`A = 1 + q'(\tau)`$; where it
    exceeds $`\nabla_{\mathrm{ad}}`$ solve the **modified** cubic of `physics/06`
    ($`8U/9 \to 8U/(9A)`$, $`W \to W_A`$, element losses kept as $`U`$): `mlt` in
    `physics/convection.rs` gets $`A`$ as an argument, default 1, and a unit test that $`A = 1`$
    reproduces the present solution. Changing $`W`$ alone is wrong. Locate the photosphere from
    the integrated $`T(\tau)`$, not from $`\tau + q = 4/3`$, once convection is present; use the same corrected gradient in the
    interior down to where $`|q'|`$ is negligible (needs $`\tau`$ along the outer mesh points).
    Do not switch from an imposed $`T(\tau)`$ to the uncorrected gradient. Accepted first version,
    to be labelled approximate: uncorrected interior with $`\tau_{\mathrm{fit}}`$ deep enough that
    $`|q'| < 10^{-3}`$;
  * element losses in the atmosphere cubic use $`U_e = U\tau_e^2/(3 + \tau_e^2)`$,
    $`\tau_e = \kappa\rho\ell`$ (decided 2026-11-11, `physics/06`), coded in the form that does
    not multiply a diverging $`U`$ by a vanishing factor; tests: $`U_e \to U`$ for
    $`\tau_e \gg 1`$, $`U_e \to U\tau_e^2/3`$ for $`\tau_e \ll 1`$;
  * test: $`\nabla`$ is continuous across the onset of convection and across the fitting point;
    a radiative atmosphere reproduces $`T(\tau)`$;
  * state which kind of $`q(\tau)`$ is used (radiative, or a fit that already contains
    convection, as Krishna Swamy's is) and that it is extrapolated below $`\tau \approx 0.02`$;
  * keep the photospheric radius and the radius of the last mesh point distinct in the code, even
    if $`\Delta r`$ is set to zero in the first version, and state the approximation.
* **Assessment:** this phase is verified on its own only by the constant-$`\kappa`$ limit and on hot
  stars. Its effect on cool stars is judged with phases 6 and 7.
* **Risk:** the integration sits inside the finite-difference Jacobian (rule 3). Use a fixed number
  of steps, no adaptive control, so the result is a smooth function of $`(R, T_{\mathrm{eff}})`$.
* **Tests:** with constant $`\kappa`$ and Eddington's $`q`$ the integral equals the present formula
  to $`10^{-6}`$; ZAMS sequence 0.4–20 M☉ still converges; hot stars change by less than 0.1%.
* **Done when:** baseline updated for cool stars; section 9 rewritten.

## Phase 4: mass loss (E1)

1. **Carry the current mass (S).** `StarModel` gets a current mass distinct from `params.mass`;
   `solve_structure`, `pulsation::fundamental_mode`, `build_model` and the exports use it. No
   behaviour change yet.
2. **Rates (S).** New `physics/mass_loss.rs`: de Jager et al. (1988) for hot stars, Reimers for
   cool ones, a switch at about 8000 K, the metallicity factor of §8.7. Unit tests against the
   table in `physics-review.md` §5.
3. **Evolution step (M).** In `Evolution::step`: $`M \leftarrow M + \dot M\,\Delta t`$; rebuild
   $`m_k = q_k M`$; interpolate the composition **and the previous model's $`\ln T`$, $`\ln P`$**
   at fixed $`m`$ onto the new mesh (new helper in `structure/mesh.rs`, next to `remap`), so that
   $`\varepsilon_{\mathrm{grav}}`$ compares the same gas (`physics/07`). Comparing at fixed $`q`$
   is not an acceptable first version: the missing term is a rate and does not shrink with the
   time step. New limit $`|\Delta M|/M \le 10^{-3}`$ per step, for the interpolation near the
   surface.
4. **Inputs and outputs (S).** A wind on/off option and an efficiency factor in
   `EvolutionOptions`; current mass and $`\dot M`$ in the track CSV, the history panels and the
   timeline.
5. **Tests.** Wind off reproduces the baseline; the hydrogen budget of §8.7 closes to 1% for
   60 M☉; the 60 M☉ star ends the main sequence with 5–15% less mass; 5 M☉ loses less than 0.1%.
6. **Energy check.** The total-energy check of 1.10, extended by the energy carried off by the
   wind, closes to 1% for 60 M☉. Expect trouble in the outermost layers, where the advection term
   nearly cancels the time derivative.

## Phase 5: rotation (E3)

### 5.1 Appearance (S–M)

* **Do:** a display-only input $`\omega = \Omega/\Omega_{\mathrm{crit}}`$ (0–0.95). Compute the Roche
  shape and $`T_{\mathrm{eff}}(\theta)`$ from §8.9 in a small pure function (in `stellar-core`, so
  it is testable), and use it in `views/model3d.rs` for the surface mesh and its colour.
* **Tests:** $`\omega = 0`$ gives a sphere; $`\omega = 1`$ gives $`R_{\mathrm{eq}}/R_{\mathrm{pole}} = 1.5`$;
  the surface integral of $`\sigma T_{\mathrm{eff}}^4`$ equals $`L`$.
* **Label it clearly** as a display of a non-rotating model until 5.2 exists.

### 5.2 Structure (M): scheduled after phase 7

* **When:** after the solar calibration at the end of phase 7 (see the overview). Not before.
* **What it is:** the spherical approximation of `physics/08`, not the Kippenhahn–Thomas method;
  $`f_T = 1`$. The radiative gradient is **divided by $`f_P`$** (`local.rs`), as spherical
  radiative diffusion requires, and the convection formulas use $`g\,f_P`$.
* **With mass loss:** removed shells take their $`j`$ with them
  ($`dJ/dt = \tfrac23\Omega_s R^2 \dot M`$); no magnetic braking. The total-energy check includes
  $`E_{\mathrm{rot}}`$.
* **Input:** initial $`\omega`$ on the ZAMS in `StarParams` (default 0).
* **Do:**
  1. `StructureProblem` gets a per-point specific angular momentum $`j_k`$; `interval` multiplies
     `dlnp_dm` by $`f_P`$ (formula in `physics-review.md` §5, second form).
  2. ZAMS: solid rotation. $`j_k`$ depends on $`r_k`$, so iterate: solve, set
     $`j = \tfrac23\Omega r^2`$, solve again, until $`R`$ changes by less than $`10^{-6}`$ (two or
     three passes).
  3. Evolution: keep $`j_k`$ in radiative layers; enforce solid rotation in convective zones
     (§8.9), using the zones of the previous model as the mixing step does.
  4. `Profile` gets $`\Omega`$; `Globals` gets surface $`v_{\mathrm{eq}}`$,
     $`\Omega/\Omega_{\mathrm{crit}}`$ and total $`J`$; the virial check uses the rotating form
     (§8.9).
* **Tests:** $`\omega = 0`$ reproduces the baseline; total $`J`$ is conserved along a run to
  $`10^{-6}`$ with the wind off, and changes by the integrated wind loss with it on; the
  temperature equation reproduces the diffusion law for a given $`L`$ at any $`\omega`$; $`L`$
  decreases monotonically with $`\omega`$ at fixed mass; the rotating virial
  check stays within $`10^{-3}`$.
* **Validity is to be measured, not assumed:** compare $`L`$, $`R`$, $`T_c`$ with published
  rotating models at two or three rotation rates and record the differences in `physics/08`. Until
  then the feature is labelled illustrative. No rotational mixing (phase 8) and no braking.

## Phase 6: equation of state

Do this before the opacity tables: the tables need a realistic $`\rho(P, T)`$ in the envelope.

### 6.1 Generic derivative layer (S)

* **Do:** a helper in `physics/eos.rs` that turns any $`\rho(P,T)`$ and $`h(P,T)`$ into `EosState`
  using §8.5 (first block), by central differences in $`\ln P`$, $`\ln T`$.
* **Test:** applied to the present ideal-gas-plus-radiation EOS it reproduces the analytic
  $`c_P`$, $`\nabla_{\mathrm{ad}}`$, $`\Gamma_1`$ to $`10^{-6}`$.
* **Consistency test (new):** a reusable test that scans an `Eos` over its $`(P,T)`$ range and
  reports the largest dimensionless residual
  $`\mathcal{R} = \rho\,(\partial h/\partial P)_{T,X_i} - (1-\delta)`$ against a combined absolute
  and relative tolerance, and checks $`c_V > 0`$ and
  $`(\partial\ln\rho/\partial\ln P)_T > 0`$ (`physics/02`). The present EOS passes to round-off. Every EOS added later must run it and state
  its result in the docs.

### 6.2 Electron degeneracy, D15b (M)

* **Decided 2026-11-11:** Fermi–Dirac electrons with ideal ions and radiation (`physics/02`). The
  Paczyński interpolation is an optional comparison model only.
* **Do:**
  1. generalised Fermi–Dirac integrals $`F_{1/2}`$, $`F_{3/2}`$, $`F_{5/2}`$ by segmented
     Gauss quadrature with break points that move smoothly with $`\eta`$ (new
     `physics/fermi_dirac.rs`); $`10^{-12}`$ on a stated $`(\eta, \vartheta)`$ domain, tested
     against a high-precision reference; calls outside the domain, or with
     $`P \le P_{\mathrm{rad}}`$, are errors; tests against the non-degenerate limit and the two
     degenerate limits ($`x_F \ll 1`$, $`x_F \gg 1`$);
  2. new `Eos` implementation; root search in $`\eta`$ to $`10^{-13}`$; add $`h`$, $`\hbar`$,
     $`m_e`$ to `constants.rs`;
  3. last and optional: the Paczyński formulas as a second implementation, with a comparison of
     $`\rho`$, $`\nabla_{\mathrm{ad}}`$, $`\Gamma_1`$ recorded in `physics/02`.
* **Cost:** the EOS is called many times per Jacobian; measure the time per ZAMS model before and
  after, and tabulate $`F_k`$ on a $`C^1`$ grid only if it is needed.
* **Consistency:** passes the 6.1 test to quadrature accuracy.
* **Also:** add the local diagnostics of
  `physics/11` (degeneracy, Coulomb coupling) to `Profile`.
* **Also touch:** `EosState.beta` is used by `sanity_check` and shown in the app; define it as
  $`1 - P_{\mathrm{rad}}/P`$. `temperature_for_density` already works for any EOS.
* **Selection:** `Physics` gets a second constructor; a physics selector in the app and
  `--eos <ideal|degenerate>` in the CLI. Keep `ideal` as default until 6.3 and 7 are in.
* **Tests:** the ideal EOS is recovered to $`10^{-6}`$ at states chosen so that the physical
  difference is below that tolerance: $`y = n_e\lambda_e^3/2 \ll 1`$ and $`kT/m_ec^2 \ll 1`$ (the
  leading correction is $`P_e/n_ekT = 1 + y/2^{5/2}`$). Not at the solar centre, where the
  electrons are weakly degenerate: there, compare with an independent Fermi–Dirac reference and
  report the actual departure from the ideal gas. Limits
  $`P \propto \rho^{5/3}`$ and $`\rho^{4/3}`$; ZAMS sequence converges down to 0.1 M☉; evolution of
  1 M☉ past hydrogen exhaustion develops an isothermal, partly degenerate core.

### 6.3 Ionisation, D15a and F2 (L)

* **Do:** ionisation equilibrium for H and He from $`\mu_j = \mu_{j+1} + \mu_e`$ with the
  Fermi–Dirac $`\eta`$ (`physics/02`), solved together with $`\eta`$ and $`n_e`$; ionisation energy
  in $`h`$. Classical Saha is the test limit, not the implementation. Add eV to `constants.rs`. H₂
  (F2) as a further species with the same structure, last.
* **Pressure ionisation (decided 2026-11-11):** no threshold switch; the Pols-type free-energy
  term of `physics/02`, from which the shift of $`\chi_j`$ and the corrections to $`P`$ and
  $`u`$ all follow. Verify the potential and its four constants against the source first
  (rule 7; only a reproduction has been read). Part of the task: a scan of the declared
  $`(P,T)`$ domain for a unique stable branch; if it fails, restrict the domain and report, or
  fall back to an embedded free-energy table. Without this term the EOS is wrong in dense
  interiors, so 6.3 is not complete, and not the default, until it is in.
* **Checks:** the 6.1 residual and stability conditions over the whole range, including across the
  ionisation zones.
* **Consequences elsewhere:** `Composition::mu` is no longer the mean molecular weight of the gas;
  everything must read `EosState.mu`. The opacity formulas assume full ionisation and are only
  meaningful together with phase 7 at low temperature.
* **Risks:** convergence in the ionisation zones, where $`\nabla_{\mathrm{ad}}`$ drops to about 0.1
  and $`\Gamma_1`$ approaches 4/3; the mesh may need more points in $`\log P`$ there.
* **Tests:** full ionisation above $`10^6`$ K reproduces 6.2; $`\nabla_{\mathrm{ad}}`$ minimum
  between 0.1 and 0.2 near $`10^4`$ K at photospheric densities; the ZAMS Sun's envelope base moves
  inward.

## Phase 7: opacity tables (D14)

1. **Decide the data source and licence (blocking, not technical).** Candidates are the OPAL or OP
   Rosseland tables for the interior and a low-temperature set (for example Ferguson et al. 2005).
   Check that the files may be redistributed inside the repository and the binary; record the
   mixture they assume. If they may not be embedded, the app needs a "load tables" path, which
   conflicts with rule 1 for the web build.
2. **Reader and interpolation (M).** New `physics/opacity_table.rs`: parse at build time or at
   start-up from embedded bytes; interpolate $`\log\kappa`$ in $`(\log T, \log R)`$ with a
   $`C^1`$ scheme (rule 3); interpolate between tables in $`X`$ and $`Z`$; blend the two temperature
   ranges (§8.6).
3. **Edges (S).** Outside the table, fall back smoothly to the analytic formula (it stays as the
   `AnalyticOpacity` implementation and as the default until the tables are validated).
4. **Evolved compositions (M).** Hydrogen-free and carbon-rich mixtures need the corresponding
   tables; until then clamp to the nearest composition and say so.
5. **Conduction (S, after 6.2).** Harmonic sum of §8.6 with a published fit.
6. **Tests.** Table nodes are reproduced exactly; derivatives are continuous across cell edges;
   table and analytic formula are compared at identical $`(\rho, T, X, Z)`$ along the present
   ZAMS solar model, and the ratio is recorded (no pass value is prescribed: "half" was a
   comparison between centres of different models). The resulting $`L`$ of the ZAMS Sun is
   recorded, not required to fall in a range.
7. **Joint assessment.** With phases 3, 6 and 7 in, swap one ingredient at a time (opacity, EOS,
   atmosphere, rates) on the 1 M☉ ZAMS model and record what each contributes to $`L`$, $`R`$ and
   the envelope base. This replaces the "suspected contributors" column of `physics/11` by
   measurements.
8. **Then calibrate.** Adjust the initial helium, the initial metallicity and
   $`\alpha_{\mathrm{MLT}}`$ so that a 1 M☉ model at 4.57 Gyr has 1 L☉, 1 R☉ and the surface
   $`Z/X`$ of the adopted mixture. Before settling exists (phase 8) the surface $`Z/X`$ equals the
   initial one, so this reduces to helium and $`\alpha_{\mathrm{MLT}}`$ against $`L`$ and $`R`$;
   state the mixture and "no diffusion" with the result. Add it as a test; rewrite section 15 of
   the About pages. Repeat the calibration after phase 8.

## Phase 8: composition-diffusion step

* **Core (M–L).** New `structure/mixing.rs`: implicit solve of §8.3 on the mass mesh for each of
  $`X`$, $`Y`$, $`Z`$, with zero-flux boundaries and exact conservation; a small tridiagonal solver
  next to `solver/banded.rs`. Called from `advance_composition` in place of the averaging.
* **Test first:** with a very large $`D`$ in convective zones and zero elsewhere it must reproduce
  the instant-mixing results (baseline) to $`10^{-6}`$; total mass of each species is conserved to
  round-off.
* **Then, each as its own task:**
  1. exponential overshoot (replaces the step overshoot of 2.1 as an option), with $`D_0`$ and
     $`H_{P,0}`$ taken at $`r_0`$, $`f_0 H_P`$ inside the convective zone; $`f`$, $`f_0`$ and the
     cut-off are separate parameters (`physics/05`);
  2. semiconvection with the Ledoux criterion (D17, §8.2);
  3. rotational mixing (E3 stage 3), needs 5.2;
  4. gravitational settling (F1), needs the advection term and diffusion velocities.

## Phase 9: new modes and phases

| Task | ID | Needs | Notes | Effort |
|---|---|---|---|---|
| 9.1 Non-radial modes, Cowling | E4 stage 2 | 2.2 | new `BvpProblem` with §8.10; frequency scan using the banded solver directly; node counting; échelle diagram | M |
| 9.2 Non-radial modes, full | E4 stage 3 | 9.1 | five unknowns per point; validate against 9.1 at high order and against the radial solver at $`l = 0`$ | M–L |
| 9.3 Mode patterns in the 3D view | E4 | 9.1 | $`Y_l^m`$ surface displacement | S |
| 9.4 Pre-main-sequence contraction | F8 | 3, 7 for accuracy | $`n = 1.5`$ start in `structure/initial.rs`; evolution with $`\varepsilon_{\mathrm{grav}}`$ only; new phase labels | M |
| 9.5 Carbon and oxygen | D18 | 6.2 | two more abundances in `Composition` and `Profile`; first fix the definitions in the table of `physics/04` ($`Z_0`$, catalysts, produced C/O, $`\mu`$, screening) and the opacity for C/O-rich matter (7.4) | M |
| 9.7 Screening for all regimes | – | 6.2, 9.5 | decided 2026-11-11: Chugunov & DeWitt (2009) for pp and CNO, static two-stage form for triple-α, exponent cap removed (`physics/04`); do it with the red-giant work. Before coding, check the two parts that are not from a published source: the electron-response extension for the weak limit and the triple-α product form. Tests: weak limit equals the present Salpeter factor; strong ion-sphere limit; continuity. The helium-ignition gate is lifted only after these pass | M |
| 9.6 Thermal-neutrino losses | F7 | 6.2 | one more term in `interval` and `inner_bc`. Do it **with** the red-giant work (right after 6.2 and conduction 7.5), not last: plasma neutrinos already matter in helium cores near $`10^8`$ K. Keep it separate from the reaction neutrinos in the $`Q`$ values | S |

## Open decisions

| Decision | Needed before | Options |
|---|---|---|
| Opacity data source and licence | phase 7 | embed in the repository; or load at run time (no web build for that mode) |
| One physics set or a selector | 6.2 | keep "analytic" as a teaching mode and add "tabulated"; or replace |
| Where new sliders go | 2.1 | main panel; or the existing Advanced group |
| Formulas marked **[check source]** in `physics/` | each task that uses one (rule 7): 1.4, 1.6, 1.7, 3, 4, 5.1, 6.2, 6.3, 7, 8, 9.1, 9.6 | verify against the cited source; record equation or page number in `physics/` and in the code |
| ~~Electron EOS for degeneracy~~ | 6.2 | **decided 2026-11-11:** Fermi–Dirac electrons + ideal ions + radiation; Paczyński as optional comparison |
| ~~Pressure ionisation~~ | 6.3 | **decided 2026-11-11:** the Pols-type free-energy term of `physics/02`, conditional on a unique stable branch over the declared domain (scan is part of 6.3). Fallback if the scan fails: embedded free-energy table |
| ~~Screening beyond the weak regime~~ | before any claim about helium ignition | **decided 2026-11-11:** Chugunov & DeWitt (2009) for binary reactions, static two-stage form for triple-α (`physics/04`). Still to check before coding: the electron-response extension for the weak limit (the review's own, unpublished) and the triple-α product form at intermediate coupling |
| ~~Thin convective elements in the atmosphere~~ | 3 | **decided 2026-11-11:** $`U_e = U\tau_e^2/(3 + \tau_e^2)`$, $`\tau_e = \kappa\rho\ell`$ (`physics/06`); no new parameter |

## Suggested first three weeks

1. Phase 0 and tasks 1.1–1.5, 1.8, 1.10. In parallel, source verification for 1.6 and phase 3.
2. Tasks 1.6, 1.7, 1.9 and phase 2.
3. Phase 3, then phase 4 or 5.1 according to interest: mass loss for correct massive-star tracks,
   the Roche shape for the most visible result. After that the block 6 → 7; structural rotation
   (5.2) only after it.
