# 12. Evolution in time

Starting from the ZAMS model, **Evolve** follows the star in time on the same mass mesh. The total mass
is constant: there is no mass loss. Each time step has two parts:

1. The composition changes by nuclear burning and convective mixing.
2. The structure is solved again, including the energy released or absorbed as layers contract or
   expand.

## Composition

*Hydrogen burning* (pp chain and CNO cycle)

```math
\frac{\partial X}{\partial t} = -\frac{\varepsilon_{pp} + \varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{H}}}, \qquad Q_{\mathrm{H}} = 6.27\times10^{18}\ \mathrm{erg\,g^{-1}}
```

$`Q_{\mathrm{H}}`$ is the pp-chain value (26.2 MeV per helium nucleus formed). The same value is used
for the CNO cycle, whose energy per helium nucleus deposited in the star is slightly lower.

*Helium burning* (triple-α); the carbon it produces counts as metals

```math
\frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}}, \qquad \frac{\partial Z}{\partial t} = \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}}, \qquad Q_{\mathrm{He}} = 5.85\times10^{17}\ \mathrm{erg\,g^{-1}}
```

*Instant mixing*

```math
X \;\leftarrow\; \bar X_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} X\, dm}{\int_{\mathrm{zone}} dm} \qquad \text{(same for } Y, Z\text{) in every convective zone}
```

* Burning uses the energy-generation rates of the previous model. This is an explicit step, kept small
  by the time-step limits below.
* The CNO catalyst abundance $`X_{\mathrm{CNO}} = f_{\mathrm{CNO}}\,Z`$ keeps its initial value.
* Mixing uses the Schwarzschild convective zones of the previous model. A convective core that shrinks
  during the main sequence therefore leaves a hydrogen gradient behind.

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
|\Delta X| \le \min\!\left(0.01,\ \max(0.3\,X,\ 10^{-5})\right), \qquad \Delta t_{\mathrm{next}} = \Delta t \cdot \min\!\left(2,\ \max\!\left(\tfrac{1}{2},\ \frac{0.02}{\max|\Delta \ln r, \Delta\ln P, \Delta\ln T, \Delta\ln L|}\right)\right)
```

* **Burning limit:** no shell may lose more than 0.01 of $`X`$, or more than 30% of its hydrogen, in
  one step. A convective zone counts as a single shell. The 30% rule has a floor of
  $`\Delta X = 10^{-5}`$, and shells with $`X \le 10^{-8}`$ are skipped, so the step does not shrink
  without bound as hydrogen runs out.
* **Structural limit:** the step grows or shrinks so that $`\ln r`$, $`\ln P`$, $`\ln T`$ and
  $`\ln L`$ change by about 0.02 per step. A step that changes them by more than 0.08 is rejected and
  repeated with half the time step.
* **Newton failure** (no convergence within 40 iterations): the step is divided by 3 and retried, up
  to 12 times. The smallest step is 1 year.
* **First step:** 10% of the burning limit, at most $`10^6`$ years.

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

The hydrogen burnt, times $`Q_{\mathrm{H}}`$, must equal the energy radiated over the run:

```math
Q_{\mathrm{H}}\left[\int X\, dm\right]_{t}^{0} = \int_0^t L\, dt'
```

For the Sun evolved to 4.57 Gyr the ratio is 0.995; an automated test checks it. Along every run:

* the virial check stays at 1 to about $`10^{-4}`$;
* the energy check, which includes $`\varepsilon_{\mathrm{grav}}`$, stays at 1.

### Time-step convergence

The same tracks were run with the step targets halved repeatedly: structural change 2%, 1%, 0.5% and
0.25%, with the hydrogen limit $`\Delta X`$ halved alongside. The results converge at **first order**,
as expected from burning explicitly with the previous model's rates: each halving roughly halves the
difference.

| Quantity | 2% (default) | 1% | 0.5% | 0.25% |
|---|---:|---:|---:|---:|
| 1 M☉ at 4.57 Gyr: $`L`$ [L☉] | 0.52123 | 0.52155 | 0.52171 | 0.52180 |
| 1 M☉ at 4.57 Gyr: $`X_c`$ | 0.51311 | 0.51290 | 0.51279 | 0.51273 |
| 3 M☉ main-sequence lifetime [Myr] | 285.89 | 285.06 | 284.61 | 284.38 |
| 3 M☉ convective core at $`X_c = 0.3`$ [$`m/M`$] | 0.11280 | 0.11276 | 0.11289 | 0.11286 |
| Steps (3 M☉ to the end of the main sequence) | 138 | 279 | 561 | 1125 |

At the default settings the errors are small compared with the physics uncertainties:

* 1 M☉ at 4.57 Gyr: $`L`$, $`R`$ and $`T_c`$ within about 0.1%;
* 3 M☉: main-sequence lifetime within about 0.6% of the extrapolated value (about 284.2 Myr);
* values exactly at the end of the main sequence scatter by up to about 1%, because the last step lands
  anywhere below $`X_c = 10^{-3}`$.

A step is rejected when it changes the structure by more than $`4\times`$ the target. In these runs the
step-size control keeps every step well below that, so a stricter threshold ($`1.5\times`$) changes
nothing. Reproduce the table with:

```sh
cargo run --release -p stellar-core --example timestep_study
```

## In the app

* **Evolution section** (left panel): choose the stop condition and press **Evolve**.
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
  properties of every timeline model. The Time evolution tab exports every step (history CSV).
* **Time evolution tab:** twelve linked time-series panels covering the whole run.
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

> Typical runs: the Sun to 4.57 Gyr takes 26 steps (0.3 s). A 5 M☉ star from the ZAMS across the
> Hertzsprung gap to 2600 K takes about 1650 steps (15 s), computed in the background while the
> timeline fills.

---
← Previous: [Initial model & continuation](11-initial-model.md) · [Contents](README.md) · Next: [Pulsations & oscillations](13-pulsations.md) →
