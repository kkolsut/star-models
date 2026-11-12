# Physics review: `docs/about/` against the implementation

Reviewed on 2026-11-10 at commit `3d48f6e`. This page compares the physics described in
[`docs/about/`](about/README.md) with the code in `crates/stellar-core`, and proposes improvements.
Nothing in the code or in the About pages was changed.

> **Branch note.** The review was made on `main`. The unmerged branch `quick-wins` (`origin/quick-wins`)
> already adds pinned regression tests (`crates/stellar-core/tests/regression.rs`), residual reporting,
> an opacity diagnostic and limit labels in the docs, and changes some quoted numbers (for example
> 17 Gyr becomes 14 Gyr). It does not touch findings B1–B4, C7–C10 or any item of sections D–F.
> Re-check B5 and B6 against that branch after it is merged.
 The equations needed for each proposal are
collected in [section 8](#8-equations-for-the-proposed-changes); the work breakdown is in
[`implementation-plan.md`](implementation-plan.md).

**Result in short.** The documented equations match the code and the textbooks they cite. There are
six places where the text is wrong or out of date, four small implementation issues (none changes a
result by more than about 0.1%), and a list of physics upgrades. The upgrade that matters most is the
opacity: it is the reason the model Sun has 0.42 L☉.

## 1. Scope and method

* **Compared:** every formula in sections 6–13 and 16 of the About pages, against
  `constants.rs`, `composition.rs`, `physics/*.rs`, `structure/*.rs`, `pulsation.rs` and
  `solver/henyey.rs`.
* **Re-derived by hand:** the thermodynamic identities ($`c_P`$, $`\nabla_{\mathrm{ad}}`$,
  $`\Gamma_1`$), the SI → CGS conversion of the opacity coefficients, the reaction $`Q`$ values, the
  photospheric pressure, and the pulsation equations in mass coordinates.
* **Re-run:** the numbers quoted in the About pages, with the existing release build of
  `stellar-cli` (built 2026-11-05, the day of the last commit).
* **Not done:** no code was changed and the test suite was not run. Reference values from
  detailed (table-based) models are quoted from the literature from memory and are approximate.

Each finding below is marked **[measured]** (CLI run), **[derived]** (hand calculation) or
**[read]** (from reading the code only).

## 2. Verdict by section

| About section | Docs vs code | Docs vs textbook | Notes |
|---|---|---|---|
| 6 Structure equations | match | match | Assumptions list is out of date (B2) |
| 7 Microphysics: EOS | match | match | $`c_P`$ equals KWW eq. 13.21 after expansion |
| 7 Microphysics: opacity | match | match | Carroll & Ostlie coefficients convert exactly |
| 7 Microphysics: nuclear | match | symbol and source wrong | B1 |
| 8 Energy transport | match | match | Ledoux remark out of date (B3) |
| 9 Boundary conditions | match | match | Constant-κ atmosphere (D13) |
| 10 Numerical method | match | n/a | Evolution's iteration limit not listed (B5) |
| 11 Initial model | match | match | |
| 12 Evolution | match, three details missing | match | B5, C7, C8 |
| 13 Pulsations | match | match | "observed" solar period (B4) |
| 16 Constants | match | CODATA 2018 / IAU 2015 | |

Numbers reproduced **[measured]**:

| Quantity | About pages | CLI |
|---|---|---|
| Sun ZAMS, $`L`$ | ≈ 0.42 L☉ | 0.4207 |
| Sun ZAMS, $`R`$ | ≈ 1.11 R☉ | 1.1133 |
| Sun ZAMS, $`T_{\mathrm{eff}}`$ | ≈ 4400 K | 4406 |
| Sun ZAMS, envelope base | 0.97 R | 0.9678 |
| $`Z = 0.008`$: $`L`$, $`T_{\mathrm{eff}}`$ | ≈ 0.83 L☉, ≈ 5500 K | 0.8348, 5515 |
| Sun to 4.57 Gyr, energy budget | 0.995 | 0.9947 |
| 5 M☉ main-sequence lifetime | 76 Myr | 76.1 Myr |

## 3. Findings

### A. Confirmed correct

* **Structure equations and discretisation.** `structure/equations.rs` implements the four
  equations in $`\ln r, \ln P, \ln T, L`$ at the interval midpoint. Band width is
  $`k_l = k_u = 5`$ as documented.
* **EOS [derived].** With $`c_V = \frac{k}{\mu m_u}\left[\tfrac32 + \tfrac{12(1-\beta)}{\beta}\right]`$
  and $`c_P - c_V = \frac{k}{\mu m_u}\frac{(4-3\beta)^2}{\beta^2}`$ the documented $`c_P`$ follows;
  $`\nabla_{\mathrm{ad}} = P\delta/(\rho T c_P)`$ then gives the documented expression.
* **Opacity [derived].** $`4.34\times10^{21}`$, $`3.68\times10^{18}`$ (SI) × $`10^{4}`$ give the CGS
  values; $`0.708 \times 1000^{0.2} = 2.82`$; $`7.9\times10^{-34} \times 10 \times \sqrt{1000} = 2.5\times10^{-31}`$.
* **Mixing-length theory.** The cubic, the substitution $`\eta = \xi - U`$, the starting bound and
  $`v_{\mathrm{conv}}`$ agree with Kippenhahn, Weigert & Weiss §7.
* **Photosphere [derived].** $`P_{\mathrm{rad}}(\tau) = \tfrac{a}{4} T_{\mathrm{eff}}^4 (\tau + \tfrac23)`$,
  so the gas pressure at $`\tau = 2/3`$ is $`\tfrac23 g/\kappa - a T_{\mathrm{eff}}^4/6`$, and the
  total is the documented value. It is consistent with the EOS radiation pressure at
  $`T = T_{\mathrm{eff}}`$.
* **Screening, $`Q_{\mathrm{He}}`$ [derived].** $`\zeta^2`$ terms are $`Z(Z+1)/A`$ for H, He and
  (approximately) CNO. $`7.275\ \mathrm{MeV} / 12\,m_u = 5.85\times10^{17}`$ erg/g.
* **Pulsations.** The two equations and three boundary conditions are the standard linear adiabatic
  wave equation (Cox 1980) rewritten with $`dr = dm / 4\pi r^2\rho`$.
* **Mesh, homology exponents, time-step rules, phases, stop conditions, constants.**

### B. Text that is wrong or out of date

Each of these needs the same edit in `docs/about/` and in `crates/stellar-app/src/views/about.rs`.

**B1. The pp-chain factor $`\psi`$ is mislabelled, and the source is misattributed [derived].**
Section 7 calls $`1 + 0.0123\,T_6^{1/3} + 0.0109\,T_6^{2/3} + 0.0009\,T_6`$ "$`\psi`$". In the
source formula this polynomial is $`g_{11}`$, a correction to the reaction rate. $`\psi`$ is a
different factor: it accounts for the pp II and pp III branches and lies between 1 and about 1.5.
The code (`physics/nuclear.rs`, variable `psi_pp`) contains only the polynomial, so it uses
$`\psi = 1`$. Also, the coefficients $`2.38\times10^{6}`$, $`8.67\times10^{27}`$ and
$`5.09\times10^{11}`$ are those of Kippenhahn & Weigert (1990), not Carroll & Ostlie as stated in
section 7, section 18 and the module comment (Carroll & Ostlie's pp coefficient converts to
$`2.41\times10^{6}`$).
*Fix:* rename the symbol to $`g_{11}`$ (equation id `nuc_psi` in `equations.json`), say that
$`\psi = 1`$ is assumed, and cite Kippenhahn & Weigert for the rates.

**B2. Section 6 "Assumptions" predates evolution.** It says the star is "chemically homogeneous …
(no burning history yet)" and that composition gradients "can be added later". Evolution already
stores and uses a composition per mesh point. *Fix:* limit the statement to ZAMS models and drop the
note.

**B3. Section 8, Ledoux criterion.** "The Ledoux criterion is not needed, because ZAMS models have
no composition gradients." Evolved models do have them (section 12 says so). *Fix:* say that the
Schwarzschild criterion is used throughout and that the stabilising effect of composition gradients
(Ledoux criterion, semiconvection) is neglected.

**B4. Section 13, solar period.** "The observed solar value is about 1 hour." The fundamental radial
mode is not observed in the Sun; about 63 minutes is the value computed for standard solar models.
The agreement of 0.99 h is also partly a coincidence, since the model has $`R = 1.11\,R_\odot`$.
*Fix:* "standard solar models give about 63 minutes".
*Checked 2026-11-11:* Model S (`fgong.l5bi.d.15`, users-phys.au.dk/jcd/solar_models) gives
258.0 μHz = 64.6 min for the radial fundamental (exact radial wave equation, linear finite elements on
the model mesh; the same code reproduces the BiSON l = 0 frequencies of Broomhall et al. 2009, Table 1,
to 0.1–0.4 μHz for n = 6–12). The About text quotes 64.6 min.

**B5. Details of section 12 that are in the code but not in the text [read].**
* Evolution steps allow 40 Newton iterations (`evolve.rs`, `Evolution::step`); section 10 lists only
  250 and 60.
* The burning limit is $`\min(0.01, \max(0.3X, 10^{-5}))`$ and skips shells with
  $`X \le 10^{-8}`$ (`burn_limit_s`).
* $`Q_{\mathrm{H}}`$ is the pp-chain value (26.2 MeV per helium nucleus) and is also used for the
  CNO cycle (see D12).

**B6. Section 4, virial check.** "Exactly 1" holds up to the surface term
$`4\pi R^3 P_{\mathrm{ph}}`$ and the cut at $`m_0`$. Both are far below the quoted
$`10^{-5}`$, so one clause is enough.

### C. Implementation issues

**C7. Burning is clipped per point before a convective zone is mixed [read; effect measured as
negligible].** In `evolve.rs`, `advance_composition` burns each point at its own rate, limits the
amount to that point's $`X`$, and only then averages the zone. The time-step limit, however, uses
the zone-averaged rate (`hydrogen_rates`). Near hydrogen exhaustion the centre of a convective core
can therefore be clipped, so slightly less hydrogen is removed than the energy released requires.
Measured budget $`Q_{\mathrm{H}}\,\Delta M_{\mathrm{H}} / \int L\,dt`$ at the end of the main
sequence: 0.9984 (1.5 M☉), 0.9991 (3 M☉), 0.9992 (5 M☉). The effect is at most 0.1%.
*Fix (optional, small):* average the rate over the zone first (reuse `hydrogen_rates`), then burn
and clip. Same for helium.

**C8. No time-step limit for helium burning [read].** `burn_limit_s` looks at $`X`$ only. Once
triple-α runs, only the structural limit (0.02 in $`\ln`$) controls the step. *Fix:* apply the same
rule to $`Y`$ where $`\varepsilon_{3\alpha} > 0`$, and document it.

**C9. Mixing length is not bounded near the centre [measured].** $`H_P \to \infty`$ as
$`g \to 0`$. In the 5 M☉ ZAMS core $`\alpha H_P / r`$ is 4900 at the innermost point, 50 at
$`q = 10^{-3}`$ and 5 at $`q = 0.03`$. The structure is unaffected because the core is adiabatic,
but the plotted $`v_{\mathrm{conv}}`$ there is formal. *Fix:* $`\ell = \min(\alpha H_P,\, r)`$ in
`physics/convection.rs` (`mlt` needs $`r`$), or a note in section 8.

**C10. Pulsation midpoints use arithmetic means [read].** `pulsation.rs` (`Lawe::interval`) averages
$`r`$, $`P`$, $`\rho`$ linearly, while the structure solver averages their logarithms. In the outer
layers, where $`P`$ falls by a large factor per interval, the two differ. *Fix:* geometric means.
Test: the period should change less under mesh refinement.

### D. Physics improvements

Effort: S = under a day, M = days, L = a week or more.

| # | Improvement | What it fixes | Effort |
|---|---|---|---|
| D11 | pp branching factor $`\psi(T, Y)`$; newer coefficients (KWW 2012: pp $`2.57\times10^{6}`$, CNO $`8.24\times10^{27}`$ with its high-$`T`$ cut-off) | pp rate about 8% low, plus up to 50% from $`\psi`$; small effect on $`L`$, visible in $`T_c`$ | S |
| D12 | Separate $`Q`$: pp $`6.27\times10^{18}`$, CNO $`\approx 5.99\times10^{18}`$ erg/g (25.0 MeV after neutrinos) | Lifetimes of CNO-burning stars about 4% too long | S |
| D16 | Step overshoot of convective cores by $`\alpha_{\mathrm{ov}} H_P`$ in the mixing step | Main-sequence lifetimes and core masses above 1.2 M☉ | S–M |
| D13 | Integrate $`dP/d\tau = g/\kappa(P, T(\tau))`$ through the atmosphere instead of using one κ | H⁻ opacity varies as $`T^9`$, so one value is a poor average in cool stars | M |
| D15b | Electron degeneracy with the Paczyński (1983) interpolation formula, behind the `Eos` trait | Stars below 0.4 M☉; helium cores, and with them a red giant branch | M |
| D14 | Tabulated opacities (OPAL/OP and low-temperature tables) behind the `Opacity` trait; electron conduction once D15b exists | See below | L |
| D15a | Saha ionisation of H and He in the EOS | Solar convection zone at 0.97 R instead of 0.71 R | L |
| D18 | $`^{12}\mathrm{C}(\alpha,\gamma)^{16}\mathrm{O}`$, with C and O tracked separately | Only needed if core helium burning is to be followed | M |
| D17 | Composition term in $`\varepsilon_{\mathrm{grav}}`$; Ledoux criterion and semiconvection | Small; list under "not modelled" | – |

**D14 in numbers [measured].** At the centre of the ZAMS Sun ($`T = 1.26\times10^{7}`$ K,
$`\rho = 64`$ g cm⁻³) the analytic formulas give $`\kappa = 2.73`$ cm² g⁻¹. Table-based solar models
have roughly half that. Since $`L \propto 1/\kappa`$ in a radiative interior, this alone accounts
for most of the factor between 0.42 and 0.7 L☉. A cheaper interim step is a better analytic fit for
the interior, but tables are the real fix and the trait is already in place.

## 4. Suggested order

1. **Text:** B1–B6 (About pages, `about.rs`, `equations.json`).
2. **Small code changes:** C8, C7, D11, D12. Check with the existing budget test in
   `crates/stellar-core/tests/evolution.rs`, extended to a 3 M☉ run to the end of the main sequence.
3. **Diagnostics:** C9, C10.
4. **Physics:** D16, D13, D15b, then D14 and D15a. The last two are what would make the solar model
   quantitative.

After D11, D12 or D14 the numbers quoted in sections 12, 13 and 15 of the About pages (0.42 L☉,
17 Gyr, 76 Myr, 0.99 h, 0.995) have to be recomputed.

## 5. Four effects not in the model: mass loss, curvature, rotation, non-radial modes

None of the four is implemented. The sizes below are estimates **[derived]** from the model's own
ZAMS values and from standard formulas quoted from memory; none was run.

| Effect | Size on the ZAMS / main sequence | Size after the main sequence | Worth adding? |
|---|---|---|---|
| Mass loss | negligible below about 15 M☉; 1% of $`M`$ at 20 M☉, 10–20% at 60 M☉ | large in every giant and supergiant phase | yes, for massive stars |
| Spacetime curvature | $`10^{-6}`$ to $`10^{-5}`$ | the same, until a compact remnant forms | no; show the redshift only |
| Rotation | negligible for the Sun; several % in structure and 20–30% in lifetime for fast-rotating A, B and O stars | cores spin up, envelopes slow down | yes, in stages |
| Non-radial modes | a property of every star, not a correction to the structure | mixed modes appear | yes, as a new output |

For Sun-like stars all four are far smaller than the opacity error (D14). For massive stars, where
the opacity is mostly electron scattering and the model is already close to detailed ones, mass loss
and rotation are the largest missing effects.

### E1. Mass loss

**Size.** With the de Jager, Nieuwenhuijzen & van der Hucht (1988) fit

```math
\log\left(-\dot M\right) = -8.158 + 1.769\,\log\frac{L}{L_\odot} - 1.676\,\log T_{\mathrm{eff}} \qquad [M_\odot\,\mathrm{yr}^{-1}]
```

the ZAMS models give:

| Mass | $`\dot M`$ [M☉/yr] | Lost on the main sequence |
|---|---|---|
| 5 M☉ | $`4\times10^{-11}`$ | 0.003 M☉ |
| 20 M☉ | $`3\times10^{-8}`$ | about 0.2 M☉ |
| 60 M☉ | $`1.3\times10^{-6}`$ (rising with $`L`$) | 5–10 M☉ |

For cool giants the usual formula is Reimers (1975),
$`\dot M = -4\times10^{-13}\,\eta\,L R / M`$ in solar units with $`\eta \approx 0.5`$. It removes
0.1–0.2 M☉ on the red giant branch of a 1 M☉ star, but this model does not reach that branch
(section 15 of the About pages), so it would do almost nothing until D15b is in.

**What has to change.**
* The mass is fixed in three places: `StructureProblem.mass` and `m` (built from $`q \times M`$),
  `solve_structure` (reads `params.mass`), and `pulsation.rs`. A current mass has to be carried in
  the model.
* The mesh is in $`q = m/M`$. When $`M`$ shrinks, a mesh point no longer follows the same gas, so:
  * composition must be interpolated at fixed $`m`$ onto the new mesh every step;
  * $`\varepsilon_{\mathrm{grav}}`$ needs the previous $`\ln P`$, $`\ln T`$ at the same $`m`$.
    A common simplification is to keep comparing at fixed $`q`$, which is accurate when
    $`|\Delta M|/M`$ per step is small.
* A new time-step limit, for example $`|\Delta M|/M \le 10^{-3}`$ per step.
* The budget check of section 12 needs the hydrogen carried away by the wind.

**Effort:** M. **Test:** with $`\dot M = 0`$ results are unchanged; total hydrogen (star plus wind)
is conserved; a 60 M☉ track ends the main sequence at lower $`L`$ than without mass loss.

### E2. Spacetime curvature (general relativity)

**Size.** The relativistic form of hydrostatic equilibrium (Tolman–Oppenheimer–Volkoff) multiplies
$`dP/dr`$ by

```math
\left(1 + \frac{P}{\rho c^2}\right)\left(1 + \frac{4\pi r^3 P}{m c^2}\right)\left(1 - \frac{2 G m}{r c^2}\right)^{-1}
```

Every factor differs from 1 by about $`GM/(Rc^2)`$:

| Star | $`GM/(Rc^2)`$ |
|---|---|
| Sun | $`2.1\times10^{-6}`$ |
| 60 M☉ ZAMS model ($`R = 10.5\,R_\odot`$) | $`1.2\times10^{-5}`$ |
| White dwarf | $`\sim10^{-4}`$ |
| Neutron star | $`\sim0.2`$ |

The model's own discretisation error (virial check) is $`10^{-5}`$ to $`10^{-4}`$, so the correction
is below what the solver resolves. It only decides anything for neutron stars, for white dwarfs near
the Chandrasekhar mass, and for supermassive stars above about $`10^5`$ M☉, where $`\Gamma_1 - 4/3`$
becomes as small as $`GM/(Rc^2)`$. None of these is in the 0.1–100 M☉ range, and runs stop at
$`T_c = 3\times10^8`$ K, long before a remnant forms.

**Proposal.** Do not put it in the solver. Two cheap additions instead (S):
* show the surface gravitational redshift $`z = GM/(Rc^2)`$ as a velocity (636 m/s for the Sun) in
  the results panel;
* add "general relativity" to the "Not modelled" list in section 15, with the number above.

### E3. Rotation

**Size.** The ratio of centrifugal to gravitational acceleration at the equator is
$`\Omega^2 R^3/(GM)`$. For the Sun (25.4 d) it is $`2\times10^{-5}`$. Stars above about 1.3 M☉ are
not slowed by a magnetised wind and typically rotate at 150–300 km/s, which is 30–60% of their
break-up speed (about 500 km/s for the 5 M☉ model). For them rotation:

1. supports part of the weight, so the star behaves like one of slightly lower mass (a few % in
   $`L`$);
2. flattens the star, up to $`R_{\mathrm{eq}}/R_{\mathrm{pole}} = 1.5`$ at break-up;
3. makes the poles hotter than the equator (von Zeipel: $`T_{\mathrm{eff}} \propto g_{\mathrm{eff}}^{1/4}`$),
   so $`L`$ and $`T_{\mathrm{eff}}`$ as seen depend on the viewing angle;
4. mixes fresh hydrogen into the core, which lengthens the main sequence by 20–30% in fast rotators
   and brings helium and nitrogen to the surface.

**Three stages.**

* **Stage 1, appearance (S–M, no solver change).** Compute the surface shape from the Roche model
  for a chosen $`\Omega/\Omega_{\mathrm{crit}}`$ and colour it by von Zeipel's law in the 3D view
  (`views/model3d.rs`). Show the range of apparent $`L`$, $`T_{\mathrm{eff}}`$ in the HR diagram.
* **Stage 2, structure (M).** The standard one-dimensional treatment (Kippenhahn & Thomas 1970;
  Endal & Sofia 1976) averages the centrifugal force over a sphere. To first order it is one factor
  in hydrostatic equilibrium:

  ```math
  \frac{dP}{dm} = -\frac{G m}{4\pi r^4}\, f_P, \qquad f_P = 1 - \frac{2\,\Omega^2 r^3}{3\,G m} = 1 - \frac{3\,j^2}{2\,G m\, r}
  ```

  where $`j = \tfrac23 \Omega r^2`$ is the specific angular momentum of a shell. The second form
  depends only on the local unknown $`r`$, so it fits `StructureProblem::interval` without changing
  the band structure. During evolution, keep $`j`$ per shell in radiative layers and enforce solid
  rotation in convective zones (reuse `convective_runs`). The core then spins up and the envelope
  slows down as they should. Valid up to roughly 70% of break-up.
* **Stage 3, mixing (L).** Rotational mixing is a diffusion process. The code only has instant
  mixing of convective zones, so this needs an implicit diffusion step for the composition. The same
  step would serve overshoot (D16) and semiconvection (D17), so it is worth building once.

**Test for stage 2:** $`\Omega = 0`$ reproduces the current models; the virial check, extended by
the rotational energy, stays at 1; $`L`$ falls as $`\Omega`$ rises.

### E4. Non-radial pulsations

**What exists.** `pulsation.rs` solves the radial ($`l = 0`$) fundamental mode only. Solar-like
oscillations come from scaling relations, not from the model's interior.

**What non-radial modes add.** Two families, set by two frequencies that the model already
determines:

```math
S_l^2 = \frac{l(l+1)\,c_s^2}{r^2}, \qquad N^2 = \frac{g^2 \rho}{P}\left[\delta\left(\nabla_{\mathrm{ad}} - \nabla\right) + \varphi\,\nabla_\mu\right], \qquad c_s^2 = \frac{\Gamma_1 P}{\rho}, \quad \nabla_\mu = \frac{d\ln\mu}{d\ln P}
```

* p modes (pressure-driven, $`\omega > S_l, N`$): the Sun's 5-minute oscillations, δ Scuti and
  β Cephei stars.
* g modes (buoyancy-driven, $`\omega < S_l, N`$), which exist only in radiative layers: γ Doradus
  stars (1.3–1.9 M☉) and slowly pulsating B stars (3–8 M☉).

On the ZAMS the two families are well separated in frequency. After the main sequence the core
contracts, $`N`$ rises there, and g modes move into the p-mode range. The result is mixed modes,
which are closely spaced and harder to find numerically.

**Three stages.**

* **Stage 1, asymptotic quantities (S).** From integrals over the existing profiles:

  ```math
  \Delta\nu = \left[2\int_0^R \frac{dr}{c_s}\right]^{-1}, \qquad \Delta\Pi_l = \frac{2\pi^2}{\sqrt{l(l+1)}}\left[\int_{N^2>0} N\,\frac{dr}{r}\right]^{-1}
  ```

  This replaces the $`\Delta\nu`$ scaling relation by a value that depends on the actual interior,
  adds the g-mode period spacing, and allows a propagation diagram ($`N`$ and $`S_l`$ against
  $`r`$) as a new profile view. $`N^2`$ should come from the formula above, not from differencing
  $`\rho`$.
* **Stage 2, modes in the Cowling approximation (M).** Neglecting the perturbation of the
  gravitational potential leaves two unknowns plus the eigenvalue: the same size as the present
  radial problem, so the existing `BvpProblem` and Henyey solver apply. It is accurate for high-order
  or high-$`l`$ modes and off by several % for the lowest ones.
* **Stage 3, full fourth-order problem (M–L).** Four unknowns plus the eigenvalue.

**The real difficulty is finding and labelling modes, not the equations.** The current method (start
Newton from a guess and accept the result if it has no nodes) does not extend to hundreds of modes.
The usual approach is to scan in frequency: for a fixed $`\omega`$ the problem without one boundary
condition is linear, one banded solve gives the residual of the remaining condition, and its sign
changes locate the modes. Modes are then labelled by counting nodes (Eckart–Scuflaire–Osaki scheme).

**Limits to state in the docs.**
* Which modes are excited still cannot be computed (no non-adiabatic analysis). The mechanism that
  drives β Cephei and slowly pulsating B stars is an opacity peak of iron near $`2\times10^5`$ K,
  which the Kramers formula does not contain. Empirical regions in the HR diagram would be used, as
  for the present strip.
* Evolved models have sharp composition steps because mixing is instant. They produce spikes in
  $`N^2`$ and exaggerated mode trapping until a diffusion step exists (E3 stage 3).
* g modes need more mesh points in the core than the structure does.

**Display.** Surface patterns $`Y_l^m`$ in the 3D view, and an échelle diagram of $`\nu_{n,l}`$.
With E3, rotational splitting $`\delta\nu \approx m\,\Omega/2\pi`$ follows at little cost.

### Order for these four

1. E4 stage 1 and the E2 redshift output (both S, no solver change).
2. E1 mass loss (M): the largest physical effect for the massive stars the model already handles well.
3. E3 stages 1 and 2.
4. E4 stages 2–3, then E3 stage 3 together with a composition-diffusion step.

## 6. Other effects to consider

This section lists what remains after sections 3 and 5. Sizes are literature values quoted from
memory, not results of this code.

### Inside the present mass range and phases

| # | Effect | Where it matters | Size | Depends on | Effort |
|---|---|---|---|---|---|
| F1 | **Element diffusion** (gravitational settling of He and metals) | Sun-like stars over Gyr | Sun: surface $`Y`$ falls by about 0.03 and $`Z`$ by about 10% in 4.6 Gyr; ages of old low-mass stars change by about 10% | composition-diffusion step (E3 stage 3) | M |
| F2 | **Molecular hydrogen and non-ideal gas** (H₂ formation, Coulomb forces, pressure ionisation) | 0.1–0.5 M☉, which the mass slider allows | H₂ lowers $`\nabla_{\mathrm{ad}}`$ in the envelope and sets where stars become fully convective; Coulomb terms are about 1% in the Sun | belongs with the EOS work (D15a, D15b) | L |
| F3 | **Non-grey atmosphere** (a $`T(\tau)`$ relation from model atmospheres or an empirical fit instead of Eddington's) | every star with a convective envelope | of order 100 K in $`T_{\mathrm{eff}}`$, a few % in $`R`$ | D13 | S after D13 |
| F4 | **Burning out of equilibrium at the ZAMS** (³He build-up; ¹²C → ¹⁴N conversion) | first 10–100 Myr; ³He takes Gyr below 0.5 M☉ | a short-lived convective core in 1 M☉ models; lower energy per reaction in M dwarfs | separate abundances for ³He, C, N | M |
| F5 | **Helium tied to metallicity** | the $`X`$ and $`Z`$ sliders | galactic enrichment gives roughly $`Y \approx 0.25 + 1.5\,Z`$; the sliders allow combinations no real star has | none | S (an optional "linked" mode) |
| F6 | **Magnetic fields and spots** | active M and K dwarfs | observed radii several % larger than standard models; also brakes rotation (E3) | none for a note | note only |

### Needed only if the model is extended to new phases

| # | Effect | Needed for | Depends on |
|---|---|---|---|
| F7 | **Thermal neutrino losses** (plasma, photo, pair) | degenerate helium cores and $`T > 5\times10^8`$ K; they cool the core and move helium ignition off-centre | D15b |
| F8 | **Pre-main-sequence contraction** (Hayashi and Henyey tracks, deuterium burning) | showing how a star arrives on the ZAMS | see below |
| F9 | **Helium flash, thermal pulses, later burning stages** | anything beyond core helium ignition | D15b, D18, F7, a hydrodynamic term for the flash |
| F10 | **Thermohaline mixing** | red giants after the luminosity bump | D15b, diffusion step |
| F11 | **Binary interaction** (mass transfer, tides, mergers) | most massive stars exchange mass with a companion | E1, E3; a different kind of program |

**F8 is the one new phase within reach.** The solver already has
$`\varepsilon_{\mathrm{grav}}`$, so a large, fully convective starting model (an $`n = 1.5`$
polytrope instead of $`n = 3`$ in `structure/initial.rs`) would contract to the ZAMS by itself. The
limit is accuracy: the position of the Hayashi track is set by the photosphere, where this model is
weakest (H⁻ formula, grey atmosphere, 2600 K floor). It would be qualitative until D13 and D14 are in.

### Already covered or not worth pursuing

* Electron conduction (with D14), turbulent pressure and time-dependent convection (small; listed in
  section 8 of the About pages), semiconvection and the composition term of
  $`\varepsilon_{\mathrm{grav}}`$ (D17).
* Departures from hydrostatic equilibrium: negligible in every phase the model reaches.
* Pair instability and the Eddington-limit envelope physics of stars above 100 M☉: outside the mass
  range.

### What this means for the order of work

Three pieces of infrastructure unlock most of sections 3, 5 and 6:

1. **Tabulated opacity and a better atmosphere** (D14, D13, F3): everything about cool stars.
2. **An EOS with ionisation, degeneracy and H₂** (D15a, D15b, F2): low-mass stars and all red-giant
   physics (F7, F9, F10).
3. **A composition-diffusion step** (E3 stage 3): overshoot, semiconvection, rotational mixing,
   element diffusion (D16, D17, F1), and smoother $`N^2`$ for E4.

F5 is independent of all three and takes an afternoon.

## 7. Is the convection model good enough?

The convection model has three parts, and the answer differs for each.

| Part | What the code does | Good enough? |
|---|---|---|
| Temperature gradient | Mixing-length theory, one fixed $`\alpha_{\mathrm{MLT}}`$ (`physics/convection.rs`) | Yes. Not the limiting factor |
| Where convection occurs | Local Schwarzschild criterion, no overshoot | No, for stars with convective cores |
| Mixing | Instant and complete inside each zone of the previous model (`evolve.rs`) | Yes for the mixing speed; the zone edges inherit the problem above |

### Temperature gradient: keep mixing-length theory

Mixing-length theory is what almost every one-dimensional stellar evolution code uses. Deep
convection is adiabatic whatever the theory, so the choice only matters in the thin superadiabatic
layer below the photosphere. **[measured]** on the ZAMS models:

* 5 M☉ core: $`\nabla - \nabla_{\mathrm{ad}} \le 1.4\times10^{-7}`$.
* Sun: $`\nabla - \nabla_{\mathrm{ad}}`$ peaks at 0.115 at $`r = 0.9997\,R`$ and is
  $`2\times10^{-5}`$ at the base of the zone.

The alternatives (full-spectrum turbulence models, $`\alpha`$ calibrated on 3D simulations as a
function of $`T_{\mathrm{eff}}`$ and $`\log g`$, time-dependent convection) change only that layer.
Time-dependent convection is needed to compute the red edge of the instability strip, which this
code takes from an empirical band anyway.

**The $`\alpha_{\mathrm{MLT}}`$ slider does almost nothing at present [measured].**

| Model | $`\alpha = 1.0`$ | $`\alpha = 1.8`$ | $`\alpha = 2.5`$ |
|---|---|---|---|
| 1 M☉: $`R`$ [R☉], $`T_{\mathrm{eff}}`$ [K] | 1.1207, 4391 | 1.1133, 4406 | 1.1106, 4411 |
| 1.3 M☉ | 1.1849, 6050 | 1.1840, 6052 | 1.1830, 6055 |
| 0.5 M☉ | 0.6996, 2180 | 0.6990, 2181 | 0.6988, 2181 |

The solar radius changes by 0.9% over the whole range. In detailed solar models the same range
changes it by roughly 10–20% (from memory). The reason is not the convection theory: the model's
solar convection zone holds only $`2\times10^{-6}`$ of the mass, because the opacity is too low at
$`10^5`$–$`10^6`$ K and there are no ionisation zones (D14, D15a). Section 3 of the About pages
("a larger $`\alpha`$ means … a smaller radius and a higher $`T_{\mathrm{eff}}`$") is correct in
sign but should give this size. A better convection theory would bring nothing until D14 and D15a
are in.

Small refinements of the existing formulation, none urgent:

* cap the mixing length near the centre (C9);
* the optically thin correction of Henyey, Vardya & Bodenheimer (1965) for blobs near the
  photosphere (S);
* turbulent pressure, about 10% of the pressure at the top of the solar convection zone; it shifts
  p-mode frequencies but hardly changes $`R`$ or $`L`$.

### Zone boundaries: the real weakness

* **No overshoot of convective cores (D16).** Rising gas does not stop where the Schwarzschild
  criterion says. Detailed models extend cores by about 0.2 $`H_P`$, calibrated on cluster diagrams
  and eclipsing binaries. This lengthens the main sequence by roughly 10–25% above about 1.2 M☉,
  widens the main-sequence band and raises the later luminosity. It is the largest
  convection-related error in the evolution results. Near 1.2 M☉ the extension must be limited by
  the core's own size, because $`H_P`$ diverges at the centre (the same issue as C9).
* **Composition jumps at the boundary (D17).** The criterion is evaluated point by point with the
  local composition. Where composition changes across the edge, opacity and $`\nabla_{\mathrm{rad}}`$
  jump too, and the Schwarzschild and Ledoux criteria disagree. This matters outside the shrinking
  cores of stars above about 15 M☉ (semiconvection) and for the growth of helium-burning cores.
  The effect is second-order next to overshoot, and the phases concerned are already qualitative.
* **No mixing below convective envelopes.** It only affects lithium and the fine structure of the
  boundary seen by helioseismology; neither is tracked.

### Instant mixing: justified

**[measured]** The convective velocity in the 5 M☉ core is about $`10^4`$ cm/s, so material crosses
the core in about a month. The main sequence lasts 76 Myr. Instant mixing is an excellent
approximation in every phase the model reaches; it fails only in late burning stages. Using the
zones of the previous model is acceptable with the present time-step limits. The side effect is that
the core edge moves in whole mesh intervals, which leaves small steps in the composition profile
(relevant for E4).

### Recommendation

1. Leave mixing-length theory as it is.
2. Add core overshoot (D16) and the mixing-length cap (C9).
3. State in the About pages that $`\alpha_{\mathrm{MLT}}`$ has little effect with the present
   opacity and EOS, with the numbers above.
4. Revisit the theory only after D14 and D15a, when $`\alpha`$ can be calibrated on the Sun.

## 8. Equations for the proposed changes

> **Superseded in part (2026-11-11).** [`physics/`](physics/README.md) is the current statement of
> these equations. After a review it corrects this section in: the atmosphere starting pressure
> (8.4), the status of the degenerate electron energy (8.5), the treatment of
> $`\varepsilon_{\mathrm{grav}}`$ with mass loss (8.7), the rotation approximation and rotational
> splitting (8.9, 8.10), the hydrogen consumption with pp branching (8.1) and the thermal-neutrino
> threshold (8.11). Later review rounds the same day also replaced: the Paczyński formulas as the
> degenerate EOS (now Fermi–Dirac electrons; 8.5), classical Saha with a pressure switch (now
> chemical-potential equilibrium; 8.5), the overshoot definition (8.3), the composition term in
> $`N^2`$ and the Ledoux criterion (8.2, 8.10), the atmosphere–convection matching (8.4) and the
> radiative gradient with rotation (8.9). This section is kept as a record of the first proposal
> only. Where the two differ, use `physics/`.

One block per proposal, in the notation of the About pages (CGS units). Status of each block:

* **[derived]**: worked out here from definitions; can be used as written.
* **[source]**: quoted from memory of the cited source. **Check every coefficient against the source
  before coding it.**

In the $`N^2`$ formula of E4, $`\varphi = (\partial\ln\rho/\partial\ln\mu)_{P,T}`$, which is 1 for
the present EOS.

### 8.1 Nuclear burning (B1, D11, D12, C7, C8, D18)

**pp branching factor [derived].** Let $`r_{11}`$, $`r_{33}`$, $`r_{34}`$ be the rates of
p + p, ³He + ³He and ³He + ⁴He. With ³He in equilibrium, $`r_{11} = 2 r_{33} + r_{34}`$, and the
helium production relative to the pp I branch alone is

```math
\phi = \frac{2\,(r_{33} + r_{34})}{r_{11}} = 1 - \alpha + \sqrt{\alpha\,(\alpha + 2)}, \qquad \alpha = \frac{\lambda_{34}^2}{\lambda_{11}\,\lambda_{33}}\left(\frac{Y}{4X}\right)^{2} = C_\alpha \left(\frac{Y}{4X}\right)^{2} \exp\!\left(-100\,T_6^{-1/3}\right)
```

$`\phi`$ runs from 1 (pp I only) to 2 (every ³He meets a ⁴He). The exponent 100 follows from the
three Gamow exponents ($`2 \times 128.1 - 33.7 - 122.5`$). **[source]** Clayton (1968) gives
$`C_\alpha \approx 1.2\times10^{17}`$ with the cross-sections of that time; with current ones it is
several times smaller, so it must be recomputed from the three $`S`$ factors.

The fraction of helium made through pp II, and the energy factor (pp II loses 25.66 MeV instead of
26.20 MeV to matter; pp III is neglected, which is safe below $`T_6 \approx 20`$):

```math
F_{\mathrm{II}} = \frac{2\,(\phi - 1)}{\phi}, \qquad \psi = \phi\,\left(1 - 0.021\,F_{\mathrm{II}}\right), \qquad \varepsilon_{pp} \;\to\; \psi\,\varepsilon_{pp}
```

**Newer rate fits [source]** (Kippenhahn, Weigert & Weiss 2012, §18.5; $`T_9 = T/10^9`$ K):

```math
\varepsilon_{pp} = 2.57\times10^{4}\;\psi\, f_{11}\, g_{11}\;\rho\, X^2\; T_9^{-2/3} \exp\!\left(-3.381\, T_9^{-1/3}\right), \qquad g_{11} = 1 + 3.82\,T_9 + 1.51\,T_9^2 + 0.144\,T_9^3 - 0.0114\,T_9^4
```

```math
\varepsilon_{\mathrm{CNO}} = 8.24\times10^{25}\; g_{14,1}\; X_{\mathrm{CNO}}\, X\, \rho\; T_9^{-2/3} \exp\!\left(-15.231\, T_9^{-1/3} - \left(\frac{T_9}{0.8}\right)^{2}\right), \qquad g_{14,1} = 1 - 2.00\,T_9 + 3.41\,T_9^2 - 2.43\,T_9^3
```

**Energy per gram, by channel [derived]** (26.731 MeV per ⁴He is $`6.398\times10^{18}`$ erg/g,
minus neutrinos):

```math
Q_{pp} = 6.27\times10^{18}, \qquad Q_{\mathrm{CNO}} = 5.99\times10^{18}\ \mathrm{erg\,g^{-1}}, \qquad \frac{\partial X}{\partial t} = -\frac{\varepsilon_{pp}}{Q_{pp}} - \frac{\varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{CNO}}}
```

**Burning in a mixed zone [derived]** (C7): average the rate first, then burn and clip.

```math
\langle \dot X \rangle_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} \dot X\, dm}{\int_{\mathrm{zone}} dm}, \qquad X_{\mathrm{zone}}^{\mathrm{new}} = \max\!\left(0,\; \bar X_{\mathrm{zone}} + \langle \dot X \rangle_{\mathrm{zone}}\,\Delta t\right)
```

**Helium time-step limit** (C8), wherever $`\varepsilon_{3\alpha} > 0`$:

```math
|\Delta Y| \le \min\!\left(0.01,\ \max(0.3\,Y,\ 10^{-5})\right)
```

**Carbon and oxygen [derived]** (D18). With $`Q_{3\alpha} = 5.85\times10^{17}`$ erg per gram of
helium and $`Q_{\mathrm{C}\alpha} = 7.162\ \mathrm{MeV}/12\,m_u = 5.76\times10^{17}`$ erg per gram
of carbon consumed:

```math
\frac{\partial X_{\mathrm{C}}}{\partial t} = \frac{\varepsilon_{3\alpha}}{Q_{3\alpha}} - \frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad \frac{\partial X_{\mathrm{O}}}{\partial t} = \frac{4}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad \frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{3\alpha}} - \frac{1}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}
```

The rate $`\varepsilon_{\mathrm{C}\alpha}`$ itself is to be taken from Kippenhahn, Weigert & Weiss
(2012, §18.5) or Caughlan & Fowler (1988); it is not reproduced here.

### 8.2 Convection (C9, D16, D17)

**Bounded scale height [source]** (the form used by Eggleton 1971 and in MESA). It needs no extra
input to `mlt`, and is finite at the centre:

```math
\tilde H_P = \min\!\left(\frac{P}{\rho\, g},\ \sqrt{\frac{P}{G\,\rho^{2}}}\right), \qquad \ell = \alpha_{\mathrm{MLT}}\,\tilde H_P
```

**Step overshoot of a convective core** (D16), applied in the mixing step only; the temperature
gradient in the extension stays radiative. $`r_{\mathrm{cc}}`$ is the Schwarzschild boundary:

```math
r_{\mathrm{ov}} = r_{\mathrm{cc}} + \alpha_{\mathrm{ov}}\,\min\!\left(H_P(r_{\mathrm{cc}}),\ r_{\mathrm{cc}}\right), \qquad \alpha_{\mathrm{ov}} \approx 0.2
```

**Ledoux criterion** (D17): a layer with $`\nabla_{\mathrm{ad}} < \nabla_{\mathrm{rad}} < \nabla_{\mathrm{L}}`$
is semiconvective.

```math
\nabla_{\mathrm{L}} = \nabla_{\mathrm{ad}} + \frac{\varphi}{\delta}\,\nabla_\mu, \qquad \nabla_\mu = \frac{d\ln\mu}{d\ln P}
```

### 8.3 Composition diffusion step (E3 stage 3, F1, diffusive overshoot)

**Mixing as diffusion [source]**, solved implicitly (one tridiagonal system per species):

```math
\frac{\partial X_i}{\partial t} = \frac{\partial}{\partial m}\left[\left(4\pi r^{2}\rho\right)^{2} D\,\frac{\partial X_i}{\partial m}\right]
```

Diffusion coefficients: convection, and exponential overshoot beyond a boundary at distance $`z`$
(Herwig 2000, $`f \approx 0.016`$):

```math
D_{\mathrm{conv}} = \tfrac{1}{3}\,v_{\mathrm{conv}}\,\ell, \qquad D_{\mathrm{ov}}(z) = D_0\,\exp\!\left(-\frac{2z}{f\,H_P}\right)
```

Instant mixing is the limit $`D \to \infty`$. Gravitational settling (F1) adds an advection term
with the diffusion velocity $`w_i`$ of each species, from Thoul, Bahcall & Loeb (1994):

```math
\frac{\partial X_i}{\partial t} = -\frac{\partial}{\partial m}\left(4\pi r^{2}\rho\, X_i\, w_i\right)
```

### 8.4 Atmosphere (D13, F3)

**Integration in optical depth [derived].** Start at small $`\tau_0`$ and integrate to the fitting
point $`\tau_{\mathrm{fit}}`$:

```math
\frac{dP}{d\tau} = \frac{g}{\kappa\!\left(\rho(P, T),\, T\right)}, \qquad T^{4}(\tau) = \tfrac{3}{4}\,T_{\mathrm{eff}}^{4}\left[\tau + q(\tau)\right], \qquad P(\tau_0) = \frac{a\,T^{4}(\tau_0)}{3} + \frac{g\,\tau_0}{\kappa_0}
```

$`\kappa_0`$ is found by iterating the last relation, as `solve_photosphere` does now. With
Eddington's $`q = 2/3`$ and constant $`\kappa`$ this reduces to the present boundary condition.

**Empirical solar relation [source]** (Krishna Swamy 1966):

```math
q(\tau) = 1.39 - 0.815\,e^{-2.54\,\tau} - 0.025\,e^{-30\,\tau}
```

With a general $`q(\tau)`$ the outermost mesh point is no longer at $`T = T_{\mathrm{eff}}`$. The
two surface conditions become

```math
L = 4\pi R^{2}\sigma T_{\mathrm{eff}}^{4}, \qquad T_N^{4} = \tfrac{3}{4}\,T_{\mathrm{eff}}^{4}\left[\tau_{\mathrm{fit}} + q(\tau_{\mathrm{fit}})\right], \qquad P_N = P_{\mathrm{atm}}(\tau_{\mathrm{fit}};\, g,\, T_{\mathrm{eff}})
```

### 8.5 Equation of state (D15a, D15b)

**Derived quantities for any EOS [derived].** Given $`\rho(P,T)`$ and the specific enthalpy
$`h(P,T)`$, by numerical differentiation:

```math
\alpha = \left(\frac{\partial\ln\rho}{\partial\ln P}\right)_{T}, \quad \delta = -\left(\frac{\partial\ln\rho}{\partial\ln T}\right)_{P}, \quad c_P = \left(\frac{\partial h}{\partial T}\right)_{P}, \quad \nabla_{\mathrm{ad}} = \frac{P\,\delta}{T\,\rho\,c_P}, \quad \Gamma_1 = \frac{1}{\alpha - \delta\,\nabla_{\mathrm{ad}}}
```

**Electron degeneracy [source]** (Paczyński 1983). $`P = P_{\mathrm{ion}} + P_e + P_{\mathrm{rad}}`$
with

```math
P_{\mathrm{ion}} = \frac{\rho\,k\,T}{\mu_i\,m_u}, \qquad P_e = \sqrt{P_{\mathrm{nd}}^{2} + P_{\mathrm{d}}^{2}}, \qquad P_{\mathrm{nd}} = \frac{\rho\,k\,T}{\mu_e\,m_u}, \qquad P_{\mathrm{d}} = \left(P_{\mathrm{nr}}^{-2} + P_{\mathrm{r}}^{-2}\right)^{-1/2}
```

```math
P_{\mathrm{nr}} = K_1\left(\frac{\rho}{\mu_e}\right)^{5/3}, \quad K_1 = \frac{(3\pi^2)^{2/3}\,\hbar^2}{5\,m_e\,m_u^{5/3}} = 1.0036\times10^{13}, \qquad P_{\mathrm{r}} = K_2\left(\frac{\rho}{\mu_e}\right)^{4/3}, \quad K_2 = \frac{(3\pi^2)^{1/3}\,\hbar c}{4\,m_u^{4/3}} = 1.2435\times10^{15}
```

```math
\frac{1}{\mu_i} = X + \frac{Y}{4} + \frac{Z}{\langle A\rangle}, \qquad \mu_e = \frac{2}{1 + X}
```

Internal energy of the electrons (needed for $`h`$):

```math
u_e = \frac{P_e}{(f - 1)\,\rho}, \qquad f = \frac{d\ln P_{\mathrm{d}}}{d\ln\rho} = \frac{5}{3}\left(\frac{P_{\mathrm{d}}}{P_{\mathrm{nr}}}\right)^{2} + \frac{4}{3}\left(\frac{P_{\mathrm{d}}}{P_{\mathrm{r}}}\right)^{2}
```

The `Eos` trait takes $`(P, T)`$, so $`\rho`$ is found by a one-dimensional root search; $`P`$ is
monotonic in $`\rho`$.

**Ionisation [source]** (Saha), for H → H⁺, He → He⁺, He⁺ → He²⁺:

```math
\frac{n_{j+1}\,n_e}{n_j} = \frac{2\,g_{j+1}}{g_j}\left(\frac{2\pi\,m_e\,k\,T}{h^{2}}\right)^{3/2} \exp\!\left(-\frac{\chi_j}{k\,T}\right)
```

| Step | $`\chi`$ [eV] | $`2 g_{j+1}/g_j`$ |
|---|---|---|
| H → H⁺ | 13.598 | 1 |
| He → He⁺ | 24.587 | 4 |
| He⁺ → He²⁺ | 54.418 | 1 |

With ionisation fractions $`x_{\mathrm{H}}`$, $`x_{\mathrm{He}1}`$ (singly) and $`x_{\mathrm{He}2}`$
(doubly), the mean molecular weight and the ionisation energy per gram (part of $`h`$) are

```math
\frac{1}{\mu} = X\,(1 + x_{\mathrm{H}}) + \frac{Y}{4}\left(1 + x_{\mathrm{He}1} + 2\,x_{\mathrm{He}2}\right) + \frac{Z}{2}, \qquad u_{\mathrm{ion}} = \frac{1}{m_u}\left[X\,x_{\mathrm{H}}\,\chi_{\mathrm{H}} + \frac{Y}{4}\left(x_{\mathrm{He}1}\,\chi_{\mathrm{He}} + x_{\mathrm{He}2}\,(\chi_{\mathrm{He}} + \chi_{\mathrm{He}^+})\right)\right]
```

Pure Saha fails at high density (pressure ionisation); the usual simple fix is to force full
ionisation above a density or pressure threshold.

### 8.6 Opacity tables (D14)

**[source]** Tables are given on a grid of $`\log T`$ and $`\log R`$, for a set of $`(X, Z)`$:

```math
R = \frac{\rho}{T_6^{3}}, \qquad \frac{1}{\kappa} = \frac{1}{\kappa_{\mathrm{rad}}} + \frac{1}{\kappa_{\mathrm{cond}}}
```

Interpolate $`\log\kappa`$ with a method that has continuous first derivatives (the Jacobian is
built by finite differences), and blend high- and low-temperature tables over an overlap, for
example $`3.75 \le \log T \le 4.0`$.

### 8.7 Mass loss (E1)

The two rate formulas are in section 5. Metallicity scaling of hot-star winds **[source]**:
$`\dot M \propto (Z/0.019)^{0.69}`$ (Vink, de Koter & Lamers 2001, eqs. 21–22; corrected in task 0.3 from $`(Z/0.02)^{0.7}`$).

**Moving mesh [derived].** With $`m_k = q_k\,M(t)`$, a quantity $`f`$ followed at fixed mass obeys

```math
\left(\frac{\partial f}{\partial t}\right)_{m} = \left(\frac{\partial f}{\partial t}\right)_{q} - q\,\frac{\dot M}{M}\,\left(\frac{\partial f}{\partial q}\right)_{t}
```

The second term is the correction to $`\varepsilon_{\mathrm{grav}}`$ and to the composition when
both are kept on the $`q`$ mesh. **Hydrogen budget with a wind**, $`X_s`$ being the surface value:

```math
\left[\int X\,dm\right]_{t}^{0} = \int_0^t \frac{L_{\mathrm{nuc,H}}}{Q_{\mathrm{H}}}\,dt' + \int_0^t X_s\,|\dot M|\,dt'
```

### 8.8 Gravitational redshift (E2)

```math
z = \frac{G M}{R\,c^{2}}, \qquad v_z = c\,z = \frac{G M}{R\,c}
```

### 8.9 Rotation (E3)

**Roche surface [derived].** With $`\omega = \Omega/\Omega_{\mathrm{crit}}`$ and
$`x(\theta) = R(\theta)/R_{\mathrm{pole}}`$ ($`\theta`$ measured from the rotation axis):

```math
\Omega_{\mathrm{crit}}^{2} = \frac{8\,G M}{27\,R_{\mathrm{pole}}^{3}}, \qquad \frac{1}{x} + \frac{4}{27}\,\omega^{2}\,x^{2}\sin^{2}\theta = 1, \qquad x(\theta) = \frac{3}{\omega\sin\theta}\,\cos\!\left[\frac{\pi + \arccos(\omega\sin\theta)}{3}\right]
```

**Effective gravity and gravity darkening [source]** (von Zeipel 1924; exponent 0.25 for radiative
envelopes, about 0.08 for convective ones). The constant is fixed by
$`\oint \sigma T_{\mathrm{eff}}^{4}\,dA = L`$:

```math
g_r = -\frac{G M}{R^{2}} + \Omega^{2} R \sin^{2}\theta, \quad g_\theta = \Omega^{2} R \sin\theta\cos\theta, \quad g_{\mathrm{eff}} = \sqrt{g_r^{2} + g_\theta^{2}}, \qquad T_{\mathrm{eff}}(\theta) \propto g_{\mathrm{eff}}^{1/4}
```

**Structure [derived].** $`f_P`$ is given in section 5. Angular momentum of a shell, solid rotation
of a convective zone, and the virial check with rotation:

```math
j = \tfrac{2}{3}\,\Omega\,r^{2}, \qquad \Omega_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} j\,dm}{\int_{\mathrm{zone}} \tfrac{2}{3}\,r^{2}\,dm}, \qquad \int_0^M \frac{3P}{\rho}\,dm + \int_0^M \tfrac{2}{3}\,\Omega^{2} r^{2}\,dm = \int_0^M \frac{G m}{r}\,dm
```

### 8.10 Non-radial pulsations (E4)

$`S_l`$, $`N`$ and the two asymptotic spacings are in section 5.

**Cowling approximation [source]** (Unno et al. 1989), for the radial displacement $`\xi_r`$ and
the Eulerian pressure perturbation $`p'`$:

```math
\frac{d\xi_r}{dr} = -\left(\frac{2}{r} - \frac{g}{c_s^{2}}\right)\xi_r + \frac{1}{\rho\,c_s^{2}}\left(\frac{S_l^{2}}{\omega^{2}} - 1\right)p', \qquad \frac{dp'}{dr} = \rho\left(\omega^{2} - N^{2}\right)\xi_r - \frac{g}{c_s^{2}}\,p'
```

```math
\frac{p'}{\rho} = \frac{\omega^{2}\,r\,\xi_r}{l} \;\; (\text{centre}), \qquad p' = \rho\,g\,\xi_r \;\; (\text{surface})
```

For the full fourth-order problem use the dimensionless variables $`y_1 \dots y_4`$ of Unno et al.
(1989, §18); they are not reproduced here.

**Frequency scan [derived].** For a trial $`\omega`$, impose the central condition and the
normalisation $`\xi_r(R) = 1`$, solve the linear system once, and evaluate

```math
\mathcal{D}(\omega) = \left[p' - \rho\,g\,\xi_r\right]_{r = R}
```

Each sign change of $`\mathcal{D}`$ brackets a mode; refine by bisection or by the existing Newton
solver with $`\omega^{2}`$ as an unknown.

### 8.11 Small items (C10, F5, F7, F8)

**Geometric means in the pulsation solver** (C10):
$`r_{k+1/2} = \sqrt{r_k\,r_{k+1}}`$, and the same for $`P`$ and $`\rho`$.

**Helium tied to metallicity [source]** (F5), with primordial $`Y_p`$:

```math
Y = Y_p + \frac{\Delta Y}{\Delta Z}\,Z, \qquad Y_p \approx 0.248, \quad \frac{\Delta Y}{\Delta Z} \approx 1.5, \qquad X = 1 - Y - Z
```

**Neutrino losses** (F7): $`dL/dm = \varepsilon_{\mathrm{nuc}} + \varepsilon_{\mathrm{grav}} - \varepsilon_\nu`$,
with $`\varepsilon_\nu(\rho, T)`$ from the fits of Itoh et al. (1996).

**Pre-main-sequence start [derived]** (F8): an $`n = 1.5`$ polytrope, for which the existing
`LaneEmden` type gives

```math
\xi_1 = 3.6538, \qquad \rho_c = 5.991\,\bar\rho, \qquad P_c = 0.770\,\frac{G M^{2}}{R^{4}}, \qquad t_{\mathrm{KH}} \approx \frac{G M^{2}}{R\,L}
```

$`t_{\mathrm{KH}}`$ sets the first time step.
