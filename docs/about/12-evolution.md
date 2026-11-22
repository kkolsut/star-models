# 12. Evolution in time

Starting from the ZAMS model, **Evolve** follows the star in time. Before every step the mesh is
refined where a burning shell or a composition step needs it, with the composition remapped so that
every species is conserved; the checkbox *remesh during evolution* (on by default) turns this off and
keeps the ZAMS mesh. The total mass
is constant unless the wind is switched on (see [Mass loss](#mass-loss); off by default). Each time step
has two parts:

1. The composition changes by nuclear burning and convective mixing.
2. The structure is solved again, including the energy released or absorbed as layers contract or
   expand.

## Composition

*Hydrogen burning* (pp chain and CNO cycle)

```math
\frac{\partial X}{\partial t} = -\frac{\varepsilon_{pp}}{Q_{pp}^{\mathrm{eff}}} - \frac{\varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{CNO}}}, \qquad Q_{pp}^{\mathrm{eff}} = 6.27\times10^{18}\,\left(1 - 0.021\,F_{\mathrm{II}}\right), \qquad Q_{\mathrm{CNO}} = 5.99\times10^{18}\ \mathrm{erg\,g^{-1}}
```

Each channel uses its own energy per gram of hydrogen burnt: 26.20 MeV per helium nucleus for the pp I
chain and 25.03 MeV for the CNO cycle, after neutrino losses. The pp II branch deposits 2.1% less than
pp I, so $`Q_{pp}^{\mathrm{eff}}`$ falls with the pp II share $`F_{\mathrm{II}}`$ (see
[Microphysics](07-microphysics.md)).

*Helium burning*: triple-α makes carbon, ¹²C(α,γ) turns it into oxygen; both are tracked apart from
the initial metals Z

```math
\frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}} - \frac{1}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad \frac{\partial X_{\mathrm{C}}}{\partial t} = \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}} - \frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad \frac{\partial X_{\mathrm{O}}}{\partial t} = \frac{4}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad Q_{\mathrm{He}} = 5.85\times10^{17}, \quad Q_{\mathrm{C}\alpha} = 5.76\times10^{17}\ \mathrm{erg\,g^{-1}}
```

*Instant mixing*

```math
X \;\leftarrow\; \bar X_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} X\, dm}{\int_{\mathrm{zone}} dm} \qquad \text{(same for } Y, Z, X_{\mathrm{C}}, X_{\mathrm{O}}\text{) in every convective zone}
```

After mixing, every mixed zone burns as one unit, at its mass-averaged rate:

```math
X_{\mathrm{zone}}^{\mathrm{new}} = \max\!\left(0,\; \bar X_{\mathrm{zone}} + \langle \dot X \rangle_{\mathrm{zone}}\,\Delta t\right), \qquad \langle \dot X \rangle_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} \dot X\, dm}{\int_{\mathrm{zone}} dm}
```

and the same for helium.

* Burning uses the energy-generation rates of the previous model. This is an explicit step, kept small
  by the time-step limits below. No zone burns more than it contains.
* The CNO catalyst abundance $`X_{\mathrm{CNO}} = f_{\mathrm{CNO}}\,Z`$ keeps its initial value.
* Mixing uses the Schwarzschild convective zones of the previous model. A convective core that shrinks
  during the main sequence therefore leaves a hydrogen gradient behind.

*Step overshoot* (input $`\alpha_{\mathrm{ov}}`$, default 0 = off). Before mixing, a convective core that
starts at the centre is extended from its Schwarzschild boundary $`r_{\mathrm{cc}}`$ to

```math
r_{\mathrm{ov}} = r_{\mathrm{cc}} + \alpha_{\mathrm{ov}}\,\min\!\left(H_P(r_{\mathrm{cc}}),\ r_{\mathrm{cc}}\right)
```

and every mesh point with $`r \le r_{\mathrm{ov}}`$ is mixed with the core.

* $`H_P(r_{\mathrm{cc}})`$ is interpolated between the mesh points around the boundary;
  $`\min(H_P, r_{\mathrm{cc}})`$ keeps the extension finite for small cores.
* Mixing only: the temperature gradient in the extension stays radiative.
* The mixed core feeds more hydrogen into the burning region. With $`\alpha_{\mathrm{ov}} = 0.2`$ the
  main sequence lasts 338 Myr instead of 269 Myr for 3 M☉ (+26%) and 88.6 Myr instead of 72.0 Myr
  for 5 M☉ (+23%). A 1.1 M☉ star, which has no convective core, is unaffected.
* $`\alpha_{\mathrm{ov}} \approx 0.2`$ is a typical calibration for intermediate-mass stars; it depends
  on mass and is not universal.

## Structure with gravitational energy

The energy equation gains the gravitational term:

```math
\frac{dL}{dm} = \varepsilon_{\mathrm{nuc}} + \varepsilon_{\mathrm{grav}}
```

```math
\varepsilon_{\mathrm{grav}} = -T\frac{ds}{dt} = -c_P T\left(\frac{d\ln T}{dt} - \nabla_{\mathrm{ad}}\frac{d\ln P}{dt}\right) \approx -\frac{c_P T}{\Delta t}\left(\Delta\ln T - \nabla_{\mathrm{ad}}\,\Delta\ln P\right)
```

$`\varepsilon_{\mathrm{grav}}`$ is evaluated implicitly, against the previous model at the same mass.
The Henyey solver therefore handles both:

* contraction ($`\varepsilon_{\mathrm{grav}} > 0`$: the layer heats up and releases gravitational
  energy);
* expansion ($`\varepsilon_{\mathrm{grav}} < 0`$).

Both happen on the thermal timescale after the main sequence. The previous model is the starting guess
for Newton, and the luminosity unknown is re-scaled to the current luminosity every step.

## Time steps

```math
|\Delta X| \le \min\!\left(0.01,\ \max(0.3\,X,\ 10^{-5})\right)\ \ (\text{same for } Y \text{ where } \varepsilon_{3\alpha} > 0), \qquad \Delta t_{\mathrm{next}} = \Delta t \cdot \min\!\left(2,\ \max\!\left(\tfrac{1}{2},\ \frac{0.02}{\max|\Delta \ln r, \Delta\ln P, \Delta\ln T, \Delta\ln L|}\right)\right)
```

* **Burning limit:** no shell may lose more than 0.01 of $`X`$, or more than 30% of its hydrogen, in
  one step; the same holds for helium wherever triple-α burns. A mixed zone counts as a single shell.
  The 30% rule has a floor of $`10^{-5}`$, and shells with $`X`$ (or $`Y`$) $`\le 10^{-8}`$ are
  skipped, so the step does not shrink without bound as the fuel runs out.
* **Structural limit:** the step grows or shrinks so that $`\ln r`$, $`\ln P`$, $`\ln T`$ and
  $`\ln L`$ change by about 0.02 per step. A step that changes them by more than 0.08 is rejected and
  repeated with half the time step.
* **Newton failure** (no convergence within 40 iterations): the step is divided by 3 and retried, up
  to 12 times. The smallest step is 1 year.
* **First step:** 10% of the burning limit, at most $`10^6`$ years.
* **Wind limit** (with mass loss): $`|\Delta M|/M \le 10^{-3}`$ per step.

## Mass loss

Off by default (**mass loss (wind)** in the Evolution section, `--wind` in the CLI). When on, the
star loses mass at the rate

```math
\log\left(-\dot M\right) = -8.158 + 1.769\,\log\frac{L}{L_\odot} - 1.676\,\log T_{\mathrm{eff}} + 0.69\,\log\frac{Z}{0.019} \qquad (T_{\mathrm{eff}} \ge 8000\ \mathrm{K};\ M_\odot\,\mathrm{yr}^{-1})
```

```math
-\dot M = 4\times10^{-13}\,\eta\,\frac{L\,R}{M} \qquad (T_{\mathrm{eff}} < 8000\ \mathrm{K};\ L, R, M \text{ in solar units};\ M_\odot\,\mathrm{yr}^{-1})
```

* **Hot stars:** the linear form of the de Jager, Nieuwenhuijzen & van der Hucht (1988) fit (§4,
  remark 3), times the metallicity factor of Vink, de Koter & Lamers (2001, eqs. 21–22). Multiplying
  the two is a choice made here, not a prescription of either paper. $`Z`$ is the surface value.
* **Cool stars:** Reimers (1975), with $`\eta = 0.477`$ by default (McDonald & Zijlstra 2015,
  calibrated on low-mass red giants in globular clusters; for other stars an extrapolation).
* The two meet as a step at 8000 K. **efficiency** (default 1) multiplies both rates.

Each step the mass changes explicitly, with the rate of the previous model. The mesh stays fixed in
$`q = m/M`$, so a mesh point no longer follows the same gas: before the structure is solved, the
composition and the previous model's $`\ln P`$, $`\ln T`$ are interpolated to the new mesh at fixed
mass $`m`$. $`\varepsilon_{\mathrm{grav}}`$ then compares the same gas. (Comparing at fixed $`q`$ would
miss a term proportional to $`\dot M`$ that does not shrink with the time step.)

```math
M \leftarrow M + \dot M\,\Delta t, \qquad m_k = q_k\,M, \qquad f^{\mathrm{prev}}_k = f^{\mathrm{old}}\!\left(q_k\,\frac{M + \dot M\,\Delta t}{M}\right), \qquad \frac{|\Delta M|}{M} \le 10^{-3}
```

Measured with the default composition ($`Z = 0.02`$) and 800 points, to the end of the main sequence:

| Mass | Ṁ on the ZAMS [M☉/yr] | Final mass | Lost |
|---|---:|---:|---:|
| 5 M☉ | $`3.9\times10^{-11}`$ | 4.995 M☉ | 0.10% |
| 20 M☉ | $`2.8\times10^{-8}`$ | 19.52 M☉ | 2.4% |
| 60 M☉ | $`1.3\times10^{-6}`$ | 51.39 M☉ | 14.3% |

The 60 M☉ main sequence ends at almost the same age (3.54 against 3.52 Myr) but at 16% lower
luminosity. Below about 15 M☉ the wind is negligible on the main sequence. The gas uncovered by the wind
decompresses; in the 60 M☉ run this absorbs about $`10^{-4}`$ of the luminosity in the outer 1% of the
mass. Halving every step limit changes the final 60 M☉ mass by $`5\times10^{-5}`$ of itself.

## Evolutionary phases

| Phase | Condition |
|---|---|
| Zero-age main sequence | age 0, homogeneous composition |
| Main sequence | core hydrogen burning, $`X_c \ge 10^{-3}`$ |
| H-shell burning (subgiant / Hertzsprung gap) | central hydrogen exhausted ($`X_c < 10^{-3}`$); hydrogen burns in a shell around the helium core |
| Red giant | shell burning, radius above 5 × the ZAMS radius and $`T_{\mathrm{eff}} < 5000`$ K |
| Core helium burning | triple-α supplies more than 10% of the nuclear energy |

## Stop conditions

| Choice | The run stops when |
|---|---|
| End of main sequence | central hydrogen is below $`10^{-3}`$ |
| Age | the chosen age is reached; the last step is shortened to land on it |
| As far as possible | $`T_{\mathrm{eff}}`$ drops below 2600 K (outside the range of the analytic opacities), $`T_c`$ exceeds $`3 \times 10^8`$ K, or the solver cannot continue even with tiny steps. Results after the main sequence are qualitative only (see [Accuracy & limitations](15-limitations.md#evolution)) |
| Safety limits (always) | 3000 steps, or an age of 30 Gyr |

## Checks

**Total-energy check.** Over a run, the energy radiated must equal the nuclear energy released minus
the change of the star's thermal and gravitational energy. $`u`$ is the thermal energy per gram; rest
mass and nuclear binding energy enter only through $`\varepsilon_{\mathrm{nuc}}`$.

```math
\mathcal{R} = \frac{\int L\,dt - \int\!\!\int\varepsilon_{\mathrm{nuc}}\,dm\,dt + \Delta(E_{\mathrm{int}} + E_{\mathrm{grav}})}{\int L\,dt}, \qquad E_{\mathrm{int}} = \int u\,dm, \qquad E_{\mathrm{grav}} = -\int \frac{G m}{r}\,dm
```

The residual has two parts, reported separately:

* The $`\varepsilon_{\mathrm{grav}}`$ of the solver leaves out the heat from the change of composition
  at fixed $`P`$ and $`T`$. The run adds that term up (the change of the enthalpy $`u + P/\rho`$); it
  does not shrink with the time step.
* The rest is numerical and shrinks with the time step.

For the Sun evolved to 4.57 Gyr the residual is $`-4.5\times10^{-4}`$: $`-4.45\times10^{-4}`$ is the
composition term and $`-9\times10^{-7}`$ numerical (halved by halving the time steps). For 5 M☉ to the
end of the main sequence they are $`-1.09\times10^{-3}`$ and $`2\times10^{-6}`$. Automated tests
require below 1% and check both parts on halving the steps.

With a wind, the gas removed carries off its thermal and gravitational energy, $`E_{\mathrm{wind}}`$
(negative: the star loses bound gas). It is reported as a third part:

```math
E_{\mathrm{wind}} = \sum_{\mathrm{steps}} \int_{M_{\mathrm{new}}}^{M_{\mathrm{old}}} \left(u - \frac{G m}{r}\right) dm, \qquad \mathcal{R} = \frac{\int L\,dt - \int\!\!\int\varepsilon_{\mathrm{nuc}}\,dm\,dt + \Delta(E_{\mathrm{int}} + E_{\mathrm{grav}}) + E_{\mathrm{wind}}}{\int L\,dt}
```

For 60 M☉ to the end of the main sequence $`E_{\mathrm{wind}}/\int L\,dt = -5.0\times10^{-4}`$, and the
residual is $`-3.16\times10^{-3}`$, nearly all of it the composition term, as without a wind
($`-3.24\times10^{-3}`$); the numerical part is $`-2.3\times10^{-5}`$.

**Fuel accounting.** The hydrogen burnt, each step with the effective energy per gram of the model it
burnt in, against the energy radiated:

```math
\sum_{\mathrm{steps}} Q_{\mathrm{H}}^{\mathrm{eff}}\,\Delta\!\left[\int X\, dm\right] \approx \int_0^t L\, dt', \qquad Q_{\mathrm{H}}^{\mathrm{eff}} = \frac{\int (\varepsilon_{pp} + \varepsilon_{\mathrm{CNO}})\,dm}{\int \left(\varepsilon_{pp}/Q_{pp}^{\mathrm{eff}} + \varepsilon_{\mathrm{CNO}}/Q_{\mathrm{CNO}}\right) dm}
```

With a wind, the hydrogen carried off, $`X_s|\Delta M|`$ with the surface value $`X_s`$, is not
counted as burnt:

```math
\Delta X_{\mathrm{burnt}} = -\Delta\!\left[\int X\,dm\right] - X_s\,|\Delta M|
```

It is close to 1 only on the main sequence, where hydrogen burning supplies nearly all the
luminosity. 60 M☉ with a wind gives 1.0012 (1.0003 without). For the Sun evolved to 4.57 Gyr it is 0.990 (0.995 with halved time steps: the difference
is the first-order error of the explicit burning step); for 3 M☉ to the end of the main sequence
0.996. Along every run:

* the virial check stays at 1 to about $`10^{-4}`$;
* the energy check, which includes $`\varepsilon_{\mathrm{grav}}`$, stays at 1.

The app shows the total-energy residual and the fuel accounting under the checks of an evolved model;
the CLI prints them, and the track and history CSV files have them as columns.

### Time-step convergence

The same tracks were run with the step targets halved repeatedly: structural change 2%, 1%, 0.5% and
0.25%, with the hydrogen limit $`\Delta X`$ halved alongside. The results converge at **first order**,
as expected from burning explicitly with the previous model's rates: each halving roughly halves the
difference.

| Quantity | 2% (default) | 1% | 0.5% | 0.25% |
|---|---:|---:|---:|---:|
| 1 M☉ at 4.57 Gyr: $`L`$ [L☉] | 0.51296 | 0.51326 | 0.51341 | 0.51349 |
| 1 M☉ at 4.57 Gyr: $`X_c`$ | 0.50958 | 0.50931 | 0.50917 | 0.50910 |
| 3 M☉ main-sequence lifetime [Myr] | 268.56 | 267.70 | 267.32 | 267.09 |
| 3 M☉ convective core at $`X_c = 0.3`$ [$`m/M`$] | 0.11221 | 0.11210 | 0.11215 | 0.11221 |
| Steps (3 M☉ to the end of the main sequence) | 137 | 276 | 556 | 1114 |

At the default settings the errors are small compared with the physics uncertainties:

* 1 M☉ at 4.57 Gyr: $`L`$, $`R`$ and $`T_c`$ within about 0.1%;
* 3 M☉: main-sequence lifetime within about 0.6% of the extrapolated value (about 266.9 Myr);
* values exactly at the end of the main sequence scatter by up to about 1%, because the last step lands
  anywhere below $`X_c = 10^{-3}`$.

A step is rejected when it changes the structure by more than $`4\times`$ the target. In these runs the
step-size control keeps every step well below that, so a stricter threshold ($`1.5\times`$) changes
nothing. Reproduce the table with:

```sh
cargo run --release -p stellar-core --example timestep_study
```

## In the app

* **Evolution section** (left panel): choose the stop condition and press **Evolve**. Tick **mass
  loss (wind)** for a run with a wind, with its **efficiency** and **Reimers η**.
  * **Stop** interrupts a run.
  * **Continue** resumes it from its latest model, for example with a later stop condition.
  * **Reset to ZAMS** discards the run.
  * **Sun's age** evolves to 4.57 Gyr.
* **Timeline** (above the views): drag to choose an age, step with **<** and **>**, **Play** through
  the run, or tick **follow latest** to watch it as it is computed. Every tab shows the selected
  model.
* **HR track and timeline models:** every step is added to the HR track. Full models are kept for the
  timeline only when one of these changes, which keeps memory bounded on runs of thousands of steps:
  * $`T_{\mathrm{eff}}`$ by more than 0.006 dex;
  * $`L`$ or $`R`$ by more than 0.015 dex;
  * $`X_c`$ by more than 0.02;
  * the phase or the pulsation state.
* **Exports:** **Pin current model** pins the selected age. **Track CSV** exports the global
  properties of every timeline model. The Time evolution tab exports every step (history CSV). Both
  have the current mass and $`\dot M`$ as columns (`M_msun`, `Mdot_msun_yr`).
* **Time evolution tab:** twelve linked time-series panels covering the whole run; a run with a wind
  adds two (mass, mass-loss rate), and the timeline shows $`M`$ and $`\dot M`$.
  * The x axis is age, or model number, which spreads out the fast phases after the main sequence.
  * Vertical lines mark phase changes, and a vertical marker shows the selected age.
  * Shaded intervals in the period panel mark the instability strip.
  * A click on any panel selects that age.
  * **fit all on screen** shrinks the panels to show all of them at once.
* Changing any input parameter discards the run, because the ZAMS star itself changes.
* **Scale:**
  * With **fixed for run**, the y axes, the $`\log P`$ axis and the colour scales span every model
    computed in the run, and are updated as models arrive.
    * Fractions use 0–1.
    * The energy-generation range starts 8 dex below each model's peak.
    * The HR diagram shows the whole diagram.
  * With **adapt to model**, each view fits the model on screen, and the HR diagram is zoomed on the
    track.
  * Switching the mode resets any manual zoom.

> Typical runs: the Sun to 4.57 Gyr takes 27 steps (0.3 s). A 5 M☉ star from the ZAMS across the
> Hertzsprung gap to 2600 K takes about 1650 steps (15 s), computed in the background while the
> timeline fills.

---
← Previous: [Initial model & continuation](11-initial-model.md) · [Contents](README.md) · Next: [Pulsations & oscillations](13-pulsations.md) →
