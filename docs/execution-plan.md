# Execution plan for the physics upgrades

The plan to work from. It turns every "Planned" block of [`physics/`](physics/README.md) into an
ordered task with files, tests and a ticket. Written 2026-11-11 against `main` at `53e86ec`.

* Equations: [`physics/`](physics/README.md). Page numbers below (p. 1 … p. 12) refer to it.
* Reasons and sizes of the effects: [`physics-review.md`](physics-review.md) (IDs B1 … F11).
* Earlier, longer plan with the history of decisions: [`implementation-plan.md`](implementation-plan.md).
  Task numbers are the same there and in `physics/` ("plan 1.6"); tasks that are new here are
  marked *new*.
* Tickets: Vikunja project "Star Models", one sub-project per stage. A ticket title starts with the
  task number in brackets, for example `[1.6]`.

Effort: **S** under a day, **M** a few days, **L** a week or more.

## Where the code stands

* `quick-wins` and the physics documents are merged; nothing of the planned physics is coded.
* `tests/regression.rs` pins 0.5, 1, 3 and 10 M☉ ZAMS models, the Sun at 4.57 Gyr and the 3 M☉
  main-sequence lifetime to 10⁻⁴.
* Baseline, 1 M☉ ZAMS: 0.4207 L☉, 1.1133 R☉, 4406 K, envelope base at 0.968 R.
* Two places where `physics/` is ahead of the code and must be corrected by the tasks named:
  * p. 9 says the pulsation solver reports ω/ω_ac. It does not. Task 1.3.
  * the README branch note still describes `main` at `3d48f6e` with `quick-wins` unmerged.
    Task 0.1.

## Stages

| Stage | Content | Tasks | Changes default results | Effort | Needs |
|---|---|---|---|---|---|
| 0 | Baseline, text fixes, source checks | 0.1–0.3 | no | S | – |
| 1 | Small corrections and consistency checks | 1.1–1.11 | slightly | M | 0 |
| 2 | Convective boundaries, seismic diagnostics | 2.1, 2.2 | no (overshoot off) | M | 1 |
| 3 | Mass loss, rotating-star display | 4.1–4.3, 5.1 | no (wind off) | M | 1 |
| 4 | Cool-star block: atmosphere, EOS, opacity, calibration | 3.1–3.3, 6.1–6.5, 7.1–7.4, 7.7–7.9 | yes, strongly | L | 1 |
| 5 | Red-giant physics | 5.R, 7.5, 9.6, 9.7, 9.5 | only past the main sequence | L | 4 |
| 6 | Rotation in the structure, mixing by diffusion | 5.2, 8.1–8.6 | no (off by default) | L | 4; 2.1 |
| 7 | New modes and phases | 9.1–9.4, 9.8 | no | L | 2.2, 4 |
| V | Verification with identical physics | V.1–V.3 | no | M | – (parallel) |

Order: 0 → 1 → {2, 3} → 4 → 5 → 6 → 7. Stages 2 and 3 are independent of each other. V.2 should
be finished before stage 4, so that later changes are not confused with mesh effects.

Inside stage 4 the order is 3.x → 6.x → 7.x. The three pieces each change the envelope and none
can be validated alone: no task before 7.7 has an acceptance range on solar L, R or the depth of
the convection zone.

## Rules for every task

1. **Physics lives in `stellar-core`**, behind the traits `Eos`, `Opacity`, `Nuclear`
   (`physics/mod.rs`). No file or thread access in the core (web build): tables are embedded with
   `include_bytes!`.
2. **New inputs are optional and default to the present behaviour.** A new `StarParams` field
   needs: `#[serde(default)]`, `validate`, `interpolate` (`structure/model.rs`); a CLI flag
   (`stellar-cli/src/main.rs`); a control in `stellar-app/src/app.rs`, `ZamsKey` in `worker.rs` if
   the ZAMS line depends on it, a row in `views/compare.rs`; the parameter tables in
   `views/about.rs` and `docs/about/03-parameters.md`.
3. **The Jacobian is built by finite differences with h = 10⁻⁶.** Any iteration inside the physics
   converges to about 10⁻¹²; any interpolation has continuous first derivatives; no adaptive step
   control inside a function the solver differentiates.
4. **Documentation moves with the code.** Move the block from "Planned" to "Implemented" in
   `physics/`; edit `views/about.rs` and the matching `docs/about/` page; change
   `stellar-app/equations/equations.json` and run `npm run render` in `tools/equations`.
5. **Every change of default results** updates `tests/regression.rs` in the same commit, with the
   reason in the commit message, and the quoted numbers (list in 1.9).
6. **Source gate.** Nothing marked **[check source]** in `physics/` is coded before it is compared
   with its source. Record the equation or page number in `physics/` and in the code comment,
   remove the marker. If the source disagrees, fix `physics/` first. Task 0.3 lists the slices.
7. **New physics is off by default** until the joint assessment 7.7 decides the defaults.
8. **Before merging:** `cargo fmt --all --check`,
   `cargo clippy --workspace --all-targets -- -D warnings`, `cargo test --workspace`.

## Stage 0: start

### 0.1 Baseline and stale statements (S)

* Add 5 and 20 M☉ ZAMS models and the 5 M☉ main-sequence lifetime to `tests/regression.rs`.
* Re-run every **[measured]** number of `physics/` on the current `main`; correct any that moved.
* Replace the branch note in `physics/README.md`.
* Done when `cargo test --workspace` passes on otherwise unchanged code.

### 0.2 About text: findings B1–B6 (S)

* Files: `views/about.rs`, `docs/about/04, 06, 07, 08, 10, 12, 13, 15, 18`, `equations.json`
  (`nuc_pp`, `nuc_psi`), `physics/nuclear.rs` (module comment, the name `psi_pp`).
* Rename ψ to g₁₁ and state ψ = 1; cite Kippenhahn & Weigert (1990) for the rates; limit "no
  burning history" to ZAMS models; say Schwarzschild is used throughout; reword the solar period;
  document the 40-iteration limit, the burning-limit floor and Q_H used for CNO; qualify "exactly
  1" for the virial check; opacity as the leading suspect only for 0.42 L☉.
* Done when the equations are re-rendered, `cargo test -p stellar-app` passes, no result changes.

### 0.3 Source verification gate (standing)

Worked in slices; the slice for a task is done before that task starts (rule 6).

| For task | To verify | Status (2026-11-12) |
|---|---|---|
| 1.4 | bounded scale height | done: formula is Paxton et al. 2011 §5.1, not in Eggleton (1971) |
| 1.6 | Kippenhahn–Weigert(–Weiss) pp and CNO fits; Solar Fusion III thermal rates and S factors; attribution of the old C_α | done: implemented KW1990 fits agree with Harris et al. 1983 / Fowler et al. 1975; SF III agrees. Decided: KWW 2012 fits not adopted, rates built from SF III; Clayton attribution dropped. Rough estimate: implemented CNO rate about 2× too high |
| 1.7 | Y_p, ΔY/ΔZ | done: corrected to MIST 0.249 / 1.5 |
| 3.3 | Krishna Swamy coefficients in the original; its optical-depth scale (τ₅₀₀₀ or Rosseland) | done: eq. 33 agrees; Rosseland depth |
| 4.1 | de Jager et al. and Reimers fits; metallicity exponent | done: de Jager agrees; exponent corrected to (Z/0.019)^0.69. Decided: η default 0.477 (McDonald & Zijlstra 2015) |
| 5.1 | gravity-darkening exponents | done: 0.25 von Zeipel, 0.08 Lucy (1967) |
| 6.3, 6.4 | ionisation energies and weights; Pols et al. (1995) potential and constants in the original | done: agree (NIST; Pols eqs. 14–18, 26–29) |
| 6.5 | Paczyński formulas, including the energy expression | done: K₁, K₂ corrected to the published m_H values |
| 7.1–7.5 | table formats, mixture, type 2 (C/O) tables, conduction fit | done: OPAL and Ferguson formats, type 2 agree; mixture corrected to GN93. Decided: OP dropped; conduction from Potekhin's table (Cassisi et al. 2007, 2021). Table licences unclear (7.1) |
| 9.6, 9.7 | Itoh et al. fits; Chugunov & DeWitt coefficients; the unpublished electron-response extension; the triple-α product form | done: Itoh reference, units, ranges; C&D t_ij agree. Decided: Potekhin & Chabrier 2012 form replaces the electron-response extension; triple-α product form kept (as in MESA) |
| 9.5 | ¹²C(α,γ)¹⁶O rate and energies | done: Q values agree (AME2020). Decided: CF88 as reproduced (ORNL and MESA agree), with a multiplier |
| 8.2, 8.5 | overshoot f₀; Thoul et al. diffusion velocities | done: f₀ = 0.005 withdrawn; Thoul agrees; solar size of settling checked (Vinyoles et al. 2017) |
| 9.1 | Cowling equations as printed in Unno et al. | done: agree with Christensen-Dalsgaard (2014), adopted as the reference (Unno not obtained) |
| all | journal volume and page numbers in p. 12 | done |

Details and source locations: the table "Source comparisons recorded" in
[`physics/12-constants-and-references.md`](physics/12-constants-and-references.md).

## Stage 1: small corrections and consistency checks

| Task | ID | Page | Change | Files | Test | Effort |
|---|---|---|---|---|---|---|
| 1.1 | C7 | 7 | Average the burning rate over a mixed zone, then burn and clip; same for helium | `structure/evolve.rs` (`advance_composition`, `hydrogen_rates`, `convective_runs`) | fuel accounting for 3 M☉ to the end of the main sequence within 0.5% | S |
| 1.2 | C8 | 7 | Helium time-step limit where ε₃α > 0 | `evolve.rs` (`burn_limit_s`) | a model with ε₃α > 0 gets a finite limit | S |
| 1.3 | C10 | 9 | Geometric means in the pulsation solver; report ω/ω_ac and flag above 0.5; verification of the mode | `pulsation.rs` (`Lawe::interval`), `RadialMode`, results panel | uniform-sphere test passes; solar period changes less between 800 and 1600 points; first overtones have the right node counts; polytrope periods (n = 1.5, 3) | S–M |
| 1.4 | C9 | 5 | Bounded scale height in the mixing length | `physics/convection.rs` (`mlt`) | unit test at g → 0; ZAMS L, R unchanged to 10⁻⁴ | S |
| 1.5 | D12 | 4 | Separate Q_pp and Q_CNO | `evolve.rs` | fuel accounting with per-channel energies | S |
| 1.10 | – | 7 | Total-energy check: ∫L dt, ∫∫ε_nuc dm dt, Δ(E_int + E_grav); rename the old check "fuel accounting" | `evolve.rs`, `structure/model.rs` (`Globals`), `output.rs`, results panel, CLI | residual below 1% for 1 and 5 M☉ on the main sequence; on halving the time step it converges to a limit, reported separately from the part that shrinks | S–M |
| 1.6 | D11 | 4 | pp branching (φ, ψ) with Q_pp^eff = Q_pp,I ψ/φ; pp and ¹⁴N(p,γ) rates from Solar Fusion III (tabulated quadrature) | `physics/nuclear.rs`, `physics/mod.rs` (`NuclearRates`), `evolve.rs` | φ → 1 at low T, φ ≤ 2; stable form for large α; fuel accounting closes; pp II share at the solar centre reported | S–M |
| 1.7 | F5 | 10 | Optional "helium follows metallicity" mode | `app.rs`, CLI flag | UI check; no core change | S |
| 1.8 | E2 | 10 | Gravitational redshift output; relativity under "Not modelled" | `structure/model.rs` (`Globals`), `output.rs`, results panel, CLI | 1 M☉ at 1 R☉ gives 636 m/s | S |
| 1.11 *new* | – | 4 | Screening diagnostic: largest exponent of a model and of a run; warn above 0.3 and at the cap | `physics/nuclear.rs`, `Globals`, results panel | exponent at the solar centre reported; a dense test state triggers each warning | S |
| 1.9 | – | all | Re-measure quoted numbers, re-pin the baseline | `docs/about/12, 13, 15`, `views/about.rs`, `physics/`, `tests/regression.rs` | – | S |

Order: 1.1, 1.2, 1.5 → 1.10 → 1.4, 1.3, 1.8, 1.11, 1.7 → 1.6 → 1.9.

*Status (2026-11-12, branch `agent/small-corrections`).* Done: 1.1–1.11. For 1.6 the S′(0), S″(0)
of Solar Fusion III Table I were checked against the source on 2026-11-12 (recorded in `physics/04`
and `12`); the pp and ¹⁴N(p,γ) rates and α are built from them by quadrature of eq. 3, tabulated
with C¹ interpolation, so no hybrid remains (¹⁴N(p,γ): Table I gives S(0) only). The measured
changes of default results are in `physics/04`, `07`, `11` and in the commit messages.

* **1.10 comes early on purpose.** Every later stage is judged with it. Conventions as on p. 7: u
  is thermal only, ε_nuc already excludes reaction neutrinos, ε_ν is thermal neutrinos only (zero
  for now). Its converged residual measures the neglected composition term of ε_grav (D17), for
  which no equation is planned.
* **1.6:** α is computed from its definition with consistently defined ⟨σv⟩; the leading-order
  C_α ≃ 1.49×10¹⁷ is a check value only. `Nuclear` must return φ or the hydrogen-consumption rate,
  not only ε_pp. Keeping the 1990 fit for the baseline with a modern α is allowed if it is called
  a hybrid and the 8% difference is stated.
* **1.5, 1.6 and 1.1 change default results** (a few % in CNO lifetimes, slightly lower T_c).
  The ¹⁴N(p,γ) rate from Solar Fusion III is roughly half the implemented one (0.3 estimate), so
  for CNO-burning stars 1.6 is expected to change T_c and lifetimes more than that; measure it.
* **1.9 list:** 0.42 L☉, 1.11 R☉, 4400 K, 0.97 R, 17 Gyr, 1.43 Gyr (1.8 M☉), 76 Myr (5 M☉),
  0.99 h, fuel ratio 0.995, the step counts quoted in About §12.

## Stage 2: convective boundaries and seismic diagnostics

### 2.1 Step overshoot of convective cores, D16 (S–M), p. 5

* Needs 1.1, 1.4.
* Input `alpha_ov` in `StarParams`, default 0, range 0–0.5 (rule 2).
* In `evolve.rs`, extend a convective run that starts at the centre to
  r_ov = r_cc + α_ov·min(H_P(r_cc), r_cc) before averaging rates and composition. Mixing only;
  the gradient stays radiative.
* Output: mixed-core mass next to the convective-core mass.
* Tests: `alpha_ov = 0` reproduces the baseline; main-sequence lifetime of 3 and 5 M☉ grows
  monotonically with α_ov; 1.1 M☉ without a convective core is unaffected; both energy checks
  close.

### 2.2 Buoyancy and Lamb frequencies, asymptotic spacings, E4 stage 1 (S–M), p. 9

* `Profile` gets `cs`, `n2`; `Globals` gets Δν_asym and one ΔΠ₁ per connected radiative region.
* N² from the production formula with φ∇_μ (exact for the present EOS); from 6.2 on with the
  general term B of p. 5. The definition g(Γ₁⁻¹ dlnP/dr − dlnρ/dr) is the cross-check only.
* Disconnected N² > 0 regions are not summed; values are labelled candidate cavities and
  asymptotic estimates.
* App: propagation diagram (N, S₁, S₂) as a profile panel; the spacings next to the scaling
  values. Columns in `output.rs`.
* Tests: in ZAMS models N² < 0 where `convective` is set, with a tolerance around zero in nearly
  adiabatic layers; in evolved models the formula agrees with the definition away from
  composition jumps, and Schwarzschild and Ledoux instability are reported separately; Δν_asym of
  the ZAMS Sun within 10% of the scaling value.

## Stage 3: mass loss and rotating-star display

### 4.1 Current mass and wind rates (S), p. 7

* `StarModel` carries a current mass distinct from `params.mass`; `solve_structure`,
  `pulsation::fundamental_mode`, `build_model` and the exports use it. No behaviour change.
* New `physics/mass_loss.rs`: de Jager et al. (hot), Reimers (cool, parameter η), switch near
  8000 K, factor (Z/0.019)^0.69 (Vink et al. 2001, eqs. 21–22). Unit tests against the table in
  `physics-review.md` §5.

### 4.2 Evolution with a wind (M), p. 7

* Needs 4.1, 1.10.
* In `Evolution::step`: M ← M + Ṁ Δt; m_k = q_k M; interpolate composition **and** the previous
  model's ln T, ln P at fixed m onto the new mesh (helper in `structure/mesh.rs`, next to
  `remap`). Comparing at fixed q is not an acceptable first version: the missing term is a rate.
* New limit |ΔM|/M ≤ 10⁻³ per step.
* Tests: wind off reproduces the baseline; M(t) equals the integrated rate to rounding; a
  constant-rate run converges on halving the step; fuel accounting with the X_s|Ṁ| term closes to
  1% for 60 M☉; 60 M☉ ends the main sequence 5–15% lighter; 5 M☉ loses less than 0.1%.
* Expect trouble in the outermost layers, where advection nearly cancels the time derivative.

### 4.3 Wind: inputs, outputs, energy check (S)

* Wind on/off and an efficiency factor in `EvolutionOptions`, app and CLI; M and Ṁ in the track
  CSV, history panels and timeline.
* Total-energy check extended by the energy carried off by the wind; closes to 1% for 60 M☉.

### 5.1 Rotation, display only: Roche shape and gravity darkening, E3 stage 1 (S–M), p. 8

* A pure function in `stellar-core` (testable): x(θ) for ω = Ω/Ω_crit in 0–0.95, g_eff(θ),
  T_eff(θ) ∝ g_eff^β with β = 0.25 (radiative envelope, von Zeipel) or about 0.08 (convective,
  Lucy 1967), normalised by L. Used by `views/model3d.rs` for the mesh and its colour.
* Label: "shape only, interior not recomputed".
* Tests: ω = 0 gives a sphere; ω = 1 gives R_eq/R_pole = 1.5; ∮σT_eff⁴ dA = L.

## Stage 4: cool-star block

### 3.1 Integrated grey atmosphere, D13 (M), p. 6

* `physics/atmosphere.rs`: integrate dP/dτ = g/κ(ρ(P,T), T) in ln τ with a fixed number of RK4
  steps (about 60) from τ₀ ≈ 10⁻⁴ to τ_fit; `TTau` enum with `Eddington`.
* Start: P(τ₀) = P(0) + gτ₀/κ₀ with P(0) = aT⁴(0)/3, **not** aT⁴(τ₀)/3; κ₀ by iteration to 10⁻¹².
* `structure/equations.rs` (`outer_bc`) and `structure/initial.rs` (`solve_photosphere`) use it.
* Tests: constant κ with q = 2/3 equals the present formula to 10⁻⁶; result independent of τ₀;
  boundary derivatives match finite differences; ZAMS sequence 0.4–20 M☉ converges; hot stars
  change by less than 0.1%.

### 3.2 Convection in the atmosphere: Hopf-corrected transport, thin elements (M), p. 6

* Needs 3.1.
* Integrate the pair (P, T) with ∇_rad,atm = A∇_rad, A = 1 + q′(τ). Where it exceeds ∇_ad solve
  η³ + (8U/9A)(η² + 2U_e η − W_A) = 0. `mlt` in `physics/convection.rs` gets A (default 1) and
  U_e. Changing W alone is wrong.
* U_e = 3acT³κ / [c_P(3 + τ_e²)] · √(8H_P/gδ), τ_e = κρℓ (the form that does not multiply a
  diverging U by a vanishing factor). U_e replaces U only in the two element-loss places.
* The same corrected gradient in the interior down to where |q′| is negligible (needs τ along the
  outer mesh points). Accepted first version, labelled approximate: uncorrected interior with
  τ_fit deep enough that |q′| < 10⁻³.
* Tests: A = 1 and τ_e → ∞ reproduce the present cubic; U_e → U for τ_e ≫ 1 and → Uτ_e²/3 for
  τ_e ≪ 1; ∇ continuous at the onset of convection and across the fitting point; a radiative
  atmosphere reproduces T(τ).

### 3.3 Non-grey T(τ) and surface conditions with a general q(τ), F3 (S–M), p. 6

* Needs 3.2 and the 0.3 slice (coefficients, optical-depth scale).
* `TTau::KrishnaSwamy`, selectable, default Eddington. State with it: a fit to the observed Sun
  (convection partly inside it), range 0.02 ≤ τ ≤ 10, extrapolated to τ = 0 for the start.
* T_eff is derived from L and R, not read from the last mesh point (`build_model` uses
  `pr.t[n - 1]` for `teff` and `nu_max_uhz`). The photosphere is located from the integrated
  T(τ). Photospheric radius and radius of the last mesh point stay distinct in the code even if
  Δr = 0 in the first version.
* Hot stars keep the grey Eddington atmosphere (range of the P_rad approximation, p. 6).
* No acceptance range on solar L or R; the change is reported.

### 6.1 EOS derivative layer and consistency test (S), p. 2

* `physics/eos.rs`: helper turning any ρ(P,T), h(P,T) into `EosState` by central differences in
  ln P, ln T. `EosState.beta` defined as 1 − P_rad/P.
* Reusable scan of an `Eos` over its (P,T) range: largest
  ℛ = ρ(∂h/∂P)_T − (1 − δ) against ε_abs + ε_rel|1 − δ|; c_V > 0; α > 0.
* Tests: the present EOS reproduces its analytic c_P, ∇_ad, Γ₁ to 10⁻⁶ and gives ℛ = 0 to
  round-off.

### 6.2 Fermi–Dirac electrons, D15b (M–L), p. 2

* Needs 6.1.
* New `physics/fermi_dirac.rs`: F₁/₂, F₃/₂, F₅/₂ by segmented Gauss quadrature with break points
  that move smoothly with η; 10⁻¹² on a stated (η, ϑ) domain against a high-precision reference;
  classical series below the lower η bound; errors outside the domain.
* New `Eos`: ideal ions + Fermi–Dirac electrons + radiation; bracketed root search in η to 10⁻¹³;
  P ≤ P_rad is an error, not a clamp. Add h, ħ, m_e, e to `constants.rs`.
* Selector: second `Physics` constructor, `--eos <ideal|degenerate>`, app selector; `ideal` stays
  the default.
* Local diagnostics of p. 11 in `Profile`: degeneracy, Coulomb coupling Γ_C, 1 − β.
* Composition derivatives (∂lnρ/∂X_i)_{P,T} for the buoyancy term B; 2.2 switches to B.
* Tests: ideal EOS recovered to 10⁻⁶ where y = n_eλ_e³/2 ≪ 1 and ϑ ≪ 1 (not at the solar centre:
  there compare with an independent reference and report the departure); limits ρ^{5/3} and
  ρ^{4/3}; ℛ to quadrature accuracy; ZAMS sequence converges down to 0.1 M☉; 1 M☉ past hydrogen
  exhaustion develops an isothermal, partly degenerate core.
* Measure the time per ZAMS model before and after; tabulate F_k on a C¹ grid only if needed.
* Until stage 5 is done: "as far as possible" runs stop when the core becomes degenerate, with a
  message.

### 6.3 Ionisation of H and He, D15a (L), p. 2

* Needs 6.2.
* Equilibrium from μ_j = μ_{j+1} + μ_e with the Fermi–Dirac η, solved together with η and n_e;
  ground-state weights; u_ion in h. Classical Saha is the test limit, not the implementation.
  Add eV to `constants.rs`.
* Everything reads `EosState.mu`; `Composition::mu` is no longer the molecular weight of the gas.
* Diagnostics: 1 − x_H, 1 − x_He2, x(1 − x) per stage in `Profile`.
* Tests: Saha limit; full ionisation at high T reproduces 6.2; ∇_ad minimum between 0.1 and 0.2
  near 10⁴ K at photospheric densities; neutral and ionised limits of μ; ℛ and stability across
  the ionisation zones.
* Not the default, and not complete, without 6.4: ideal equilibrium recombines in dense interiors.
* Risk: Newton convergence where Γ₁ → 4/3; the mesh may need more points in log P there.

### 6.4 Pressure ionisation: Pols-type free-energy term (M–L), p. 2

* Needs 6.3 and the 0.3 slice (original paper).
* F_PI with the four published constants; μ_e,PI, P_PI, u_PI all derived from it; ionisation
  energies shifted by μ_e,PI. No threshold switch.
* Part of the task: scan of the declared (P,T) domain for a unique stable branch. If it fails,
  restrict the domain and report; no branch switching or pressure floors. Fallback: an embedded
  free-energy table.
* Tests: dilute gas recovers 6.3; n = n₀ gives F_PI = P_PI = u_PI = 0; no recombination in a
  dense interior; equilibrium residuals, c_V > 0, α > 0 and ℛ over the domain.

### 6.5 Molecular hydrogen; Paczyński comparison model (M, optional pieces), p. 2

* H₂ (F2) as a further species with a dissociation equilibrium of the same structure. Decides
  full convection at 0.1–0.5 M☉.
* Paczyński interpolation as a second selectable `Eos`, comparison only: run the 6.1 scan on it
  and record ρ, ∇_ad, Γ₁ against 6.2 along a model in `physics/02`.
* Coulomb corrections are not planned.

### 7.1 Opacity data: source, licence, mixture (S, blocking), p. 3

* High-temperature set: OPAL (decided in 0.3; OP dropped). Low-temperature set: Ferguson et al.
  (2005). Conduction: Potekhin's table of Cassisi et al. (2007, 2021). Check that the files may be redistributed in the repository and the binary.
* Record the heavy-element mixture; the same mixture defines f_CNO in the rates and ⟨A⟩ in the
  EOS.
* If the data may not be embedded, a run-time "load tables" path is needed, which conflicts with
  rule 1 for the web build: decide before coding.

### 7.2 Table reader and C¹ interpolation, D14 (M), p. 3

* Needs 7.1.
* New `physics/opacity_table.rs` behind `Opacity`: embedded bytes; log κ in (log T, log R),
  R = ρ/T₆³, with continuous first derivatives; interpolation in X and Z.
* Tests: nodes reproduced exactly; derivatives continuous across cell edges and equal to finite
  differences.

### 7.3 Blending, edges, comparison at identical state (S–M), p. 3

* Needs 7.2.
* Blend the two temperature ranges over 3.75 ≤ log T ≤ 4.0; smooth fallback to the analytic
  formulas outside the tables; `AnalyticOpacity` stays selectable and the default.
* Table against analytic formula at identical (ρ, T, X, Z) along the present ZAMS solar model:
  record the ratio (`examples/opacity_map.rs`). No pass value. Record the resulting L of the ZAMS
  Sun.
* ZAMS sequence converges for the baseline masses with tables on.

### 7.4 Opacity for evolved compositions (M), p. 3

* Needs 7.2.
* Hydrogen-free tables; type 2 (excess C and O) tables once 9.5 exists. Until then clamp to the
  nearest composition and say so in the app; results after helium ignition stay qualitative.

### 7.7 Joint assessment of atmosphere, EOS and opacity (M), p. 11

* Needs 3.x, 6.1–6.4, 7.1–7.3, V.2.
* On the 1 M☉ ZAMS model, swap one ingredient at a time (opacity, EOS, atmosphere, rates) and all
  together. Table of L, R, T_eff, envelope base, T_c, ρ_c, central κ, ∇_ad in the ionisation
  zones, photospheric pressure.
* The table replaces the "suspected contributors" column of p. 11. Decide the new defaults from
  it (rule 7). Re-pin `tests/regression.rs`.

### 7.8 Provisional solar calibration (S–M), p. 11

* Needs 7.7.
* Adjust Y₀ and α_MLT so that 1 M☉ at 4.57 Gyr has L☉ and R☉ (without settling the surface Z/X
  equals the initial one). State the mixture and "no diffusion" with the result. Add it as a
  test; rewrite About §15 and the validity table of p. 11.

### 7.9 Validation against published models (M)

* Needs 7.8.
* Reference set, cited: a standard solar model; ZAMS values and main-sequence lifetimes for
  0.5–20 M☉ from a published grid.
* Table of relative differences in L, R, T_eff, T_c, ρ_c, envelope base, lifetimes. Each
  difference above a stated tolerance gets a named suspect or an open item. The validity ranges in
  the app follow from this table.

## Stage 5: red-giant physics

### 5.R Remeshing during evolution (M) *not in `physics/`; prerequisite*

* The mesh is fixed in q (p. 7, "Limits"); a hydrogen-burning shell is far thinner than the
  interior spacing. Reuses the fixed-m interpolation of 4.2.
* Refinement criteria on gradients of ln P, ln T, X, ε; conservative remapping of composition.
* Tests: remeshing a converged model changes L and R by less than the solver tolerance; species
  mass conserved to rounding; shell luminosity converges with the point count.

### 7.5 Electron conduction (S), p. 3

* Needs 6.2, 7.2. 1/κ = 1/κ_rad + 1/κ_cond with a published fit chosen with the tables.
* Tests: negligible at the solar centre; dominant in a degenerate helium core.

### 9.6 Thermal-neutrino losses, F7 (S), p. 1

* Needs 6.2. −ε_ν in `interval` and `inner_bc` (`structure/equations.rs`) and in the total-energy
  check; Itoh et al. (1996) fits, applied wherever they are valid, no temperature threshold. Kept
  separate from the reaction neutrinos in the Q values.
* Tests: negligible on the main sequence; the temperature maximum moves off-centre in a
  degenerate core; the energy residual closes.

### 9.7 Screening for all regimes (M), p. 4

* Needs 6.2 (∂n_e/∂μ_e). Settled in 0.3: code the electron screening factor of Potekhin &
  Chabrier (2012, eqs. 20–21), not the review's extension; keep the triple-α product of pair
  factors (as MESA's `eval_screen_triple`).
* Chugunov & DeWitt (2009) for pp and p + ¹⁴N with the actual mixture; static two-stage form for
  triple-α (t_ij = 1); H kept as a logarithm; the cap of 2 removed.
* Tests: weak limit equals the present Salpeter factor; triple-α weak limit 12e²/(kTD); strong
  limit 0.9Γ_e[Z_c^{5/3} − Z_i^{5/3} − Z_j^{5/3}]; continuous derivatives through intermediate
  coupling; solar-centre rates change by less than a stated tolerance.
* The gate on helium ignition is lifted only after these pass, and by a separate decision.

### 9.5 Composition vector, carbon and oxygen, D18 (M–L), p. 4

* Needs 6.2; 7.4 for the opacity.
* First fix the definitions of the table on p. 4: Z₀, X_CNO, produced X_C and X_O, μ/μ_e/μ_i from
  all species, screening ζ. One state used by EOS, opacity and rates.
* ¹²C(α,γ)¹⁶O with the abundance equations and Q_Cα of p. 4; time-step limit on C.
* Tests: mass fractions sum to 1 to rounding; fuel accounting per channel; C/O ratio of a
  helium-burning core reported.
* F4 (³He, ¹²C, ¹⁴N out of equilibrium) becomes possible; not scheduled.

Out of reach after this stage, and stated in the app: helium flash, thermal pulses.

## Stage 6: rotation in the structure, mixing by diffusion

### 5.2 Rotation in the structure equations, E3 stage 2 (M), p. 8

* After 7.8. Spherical approximation (f_T = 1), not the Kippenhahn–Thomas method.
* Input: initial ω on the ZAMS in `StarParams`, default 0.
* `StructureProblem` gets j_k; `interval` multiplies `dlnp_dm` by f_P = 1 − 3j²/(2Gmr);
  `local.rs` divides ∇_rad by f_P; convection uses g·f_P.
* ZAMS: solid rotation, iterate j = ⅔Ωr² with the solve until R changes by less than 10⁻⁶.
  Evolution: j conserved in radiative layers, solid rotation in convective zones (zones of the
  previous model). With a wind, removed shells take their j.
* Outputs: Ω in `Profile`; v_eq, Ω/Ω_crit, J in `Globals`; rotating virial check; E_rot in the
  total-energy check.
* Tests: ω = 0 reproduces the baseline; J conserved to 10⁻⁶ with the wind off and changing by the
  integrated loss with it on; the temperature equation is the diffusion law at any ω; L falls
  monotonically with ω; rotating virial check within 10⁻³.
* Validity is measured: compare L, R, T_c with published rotating models at two or three rates,
  record in `physics/08`; labelled illustrative until then.

### 8.1 Composition-diffusion step (M–L), p. 5

* Needs 2.1 (as the reference for overshoot), 1.4.
* New `structure/mixing.rs`: implicit solve of ∂X_i/∂t = ∂/∂m[(4πr²ρ)² D ∂X_i/∂m] per species,
  zero-flux boundaries, exact conservation; tridiagonal solver next to `solver/banded.rs`. Called
  from `advance_composition`; D_conv = ⅓ v_conv ℓ.
* Test first: very large D in convective zones and zero elsewhere reproduces instant mixing to
  10⁻⁶; species mass conserved to round-off.

### 8.2 Exponential overshoot (S–M), p. 5

* Needs 8.1. D₀ and H_P,0 at r₀, a distance f₀H_P inside the boundary; f, f₀ and the cut-off are
  separate stated parameters; bounded H_P for small cores. Option next to the step overshoot.
* Tests: f = 0 reproduces 8.1; lifetime grows monotonically with f; smooth profile edge.

### 8.3 Ledoux criterion and semiconvection, D17 (M), p. 5

* Needs 8.1; 6.2 for the general term B. ∇_L = ∇_ad + B/δ; semiconvective layers get a slow
  diffusion coefficient. Schwarzschild and Ledoux flags reported separately (shared with 2.2).

### 8.4 Rotational mixing, E3 stage 3 (M), p. 8

* Needs 8.1, 5.2. Circulation and shear as terms in D; the coefficients are chosen when the task
  starts (none are specified in `physics/`).

### 8.5 Gravitational settling, F1 (M), p. 5

* Needs 8.1 and the 0.3 slice. Advection term with the diffusion velocities of Thoul, Bahcall &
  Loeb (1994).
* Test: Sun over 4.6 Gyr, surface Y falls by about 0.03 and Z by about 10%.

### 8.6 Final solar calibration (S–M), p. 11

* Needs 8.5, 7.8. Three parameters (Y₀, Z₀, α_MLT) against L☉, R☉ and the surface (Z/X)☉ of the
  adopted mixture. Replaces the provisional values; re-run 7.9.

## Stage 7: new modes and phases

| Task | ID | Page | Needs | Content | Effort |
|---|---|---|---|---|---|
| 9.1 | E4 stage 2 | 9 | 2.2 | Non-radial modes in the Cowling approximation, l ≥ 1: new `BvpProblem`; frequency scan with the two-sided matching determinant, poles told from roots; labels by node count. Tests: polytrope eigenfrequencies from the literature; spacings approach the asymptotic values of 2.2 at high order | M |
| 9.2 | E4 stage 3 | 9 | 9.1 | Full fourth-order problem from Christensen-Dalsgaard (2014), eqs. 4.61–4.68 (Unno et al. not available); validate against 9.1 at high order and against the radial solver at l = 0 | M–L |
| 9.3 | E4 | 9 | 9.1 | Échelle diagram; Y_l^m surface patterns in the 3D view | S |
| 9.8 *new* | E3 | 8 | 9.1; 5.2 for Ω | Rotational splitting: Ledoux constant C_nl from the eigenfunctions, δν = m(1 − C_nl)Ω/2π. Tests: C_nl → 0 for high-order p modes, → 1/[l(l+1)] for high-order g modes | S |
| 9.4 | F8 | 7 | 3.x, 7.x for accuracy | Pre-main-sequence contraction: n = 1.5 start in `structure/initial.rs`; ε_grav only until hydrogen ignites; first step from t_KH; new phase labels. Tests: L = −d(E_int + E_grav)/dt; the track ends on the ZAMS model of the same mass | M |

Limits to document with 9.1–9.3: no excitation; composition steps from instant mixing exaggerate
mode trapping (better after 8.1); g modes need more points in the core.

## Track V: verification with identical physics

Independent of the physics upgrades; from the scientific review of 2026-11-05.

| Task | Content | Effort |
|---|---|---|
| V.1 | Independent solver (shooting to a fitting point) with the same analytic EOS, opacity, rates, convection and atmosphere; compare 1 and 3 M☉ ZAMS globals and profiles; store as reference data for CI | M |
| V.2 | Mesh and central-cut-off convergence: 400/800/1600/3200 points for ZAMS, mid-main-sequence and end-of-main-sequence models; local profiles, not only globals; observed order; m₀ from 10⁻⁸ to 10⁻⁴ M. Before stage 4 | S–M |
| V.3 | Solver robustness: finite-difference Jacobian against directional derivatives; effect of the kink at ∇_rad = ∇_ad on Newton convergence; positivity guards in trial steps; independence of the starting model | S–M |

## Open decisions

| Decision | Needed before | Options |
|---|---|---|
| Opacity data source and licence | 7.2 | embed in the repository; or load at run time (no web build for that mode) |
| One physics set or a selector | 6.2 | keep "analytic" as a teaching mode and add the new set; or replace |
| Where new controls go | 2.1 | main panel; or the Advanced group |
| Conduction fit | 7.5 | to be chosen with the tables |
| Rotational-mixing coefficients | 8.4 | to be chosen when 8.1 exists |
| Lifting the helium-ignition gate | after 9.7 | by decision, on the test results |

Decided 2026-11-11 and not reopened here: Fermi–Dirac electrons with ideal ions and radiation;
Pols-type pressure-ionisation term, conditional on the stability scan; Chugunov & DeWitt screening
with a static two-stage triple-α; thin-element term U_e in the atmosphere cubic; structural
rotation after the cool-star block.

## Not planned

| Item | Reason |
|---|---|
| General relativity in the structure (E2) | 10⁻⁶ to 10⁻⁵, below the discretisation error; redshift output only (1.8) |
| Other convection theories, turbulent pressure | change only the thin superadiabatic layer |
| Coulomb corrections in the EOS | about 1% in the Sun; a tabulated EOS would include them |
| Composition term in ε_grav (D17) | no equation yet; measured by the residual of 1.10 |
| Non-adiabatic pulsations (excitation) | needs the iron opacity peak and all of stages 4 and 7; ticket kept as "later" |
| Out-of-equilibrium ³He, ¹²C, ¹⁴N (F4) | possible after 9.5; not scheduled |
| Helium flash, thermal pulses (F9) | need a hydrodynamic term |
| Magnetic fields and spots (F6), binaries (F11) | out of scope |

## Coverage of `physics/`

| Page | Planned block | Task |
|---|---|---|
| 1 | thermal neutrinos / rotation factor / variable mass | 9.6 / 5.2 / 4.2 |
| 2 | derivative layer / Fermi–Dirac / Paczyński / ionisation / pressure ionisation / H₂ | 6.1 / 6.2 / 6.5 / 6.3 / 6.4 / 6.5 |
| 3 | tables / mixture / evolved compositions / conduction | 7.2, 7.3 / 7.1 / 7.4 / 7.5 |
| 4 | branching, newer fits / Q by channel / C and O / screening diagnostic / screening / F4 | 1.6 / 1.5 / 9.5 / 1.11 / 9.7 / not planned |
| 5 | bounded H_P / step overshoot / Ledoux / diffusion, exponential overshoot / settling / atmosphere cubic | 1.4 / 2.1 / 8.3 / 8.1, 8.2 / 8.5 / 3.2 |
| 6 | integrated atmosphere / Hopf transport, U_e / non-grey T(τ), surface conditions | 3.1 / 3.2 / 3.3 |
| 7 | zone burning / helium limit / total energy / mass loss / pre-main sequence | 1.1 / 1.2 / 1.10 / 4.1–4.3 / 9.4 |
| 8 | shape / structure / mixing / splitting | 5.1 / 5.2 / 8.4 / 9.8 |
| 9 | geometric means, cut-off / frequencies, spacings / Cowling / full / mode finding | 1.3 / 2.2 / 9.1 / 9.2 / 9.1 |
| 10 | redshift / rotating virial / helium tied to Z | 1.8 / 5.2 / 1.7 |
| 11 | local diagnostics / swap experiment / calibration | 6.2, 6.3, 1.11 / 7.7 / 7.8, 8.6 |
| 12 | new constants / unverified sources | 6.2, 6.3 / 0.3 |
