# 4. Results: global quantities

These are shown in the left panel, the Compare table and the JSON export.

| Quantity | Meaning |
|---|---|
| $`R`$ [R☉] | Radius of the photosphere (optical depth $`\tau = 2/3`$), i.e. $`r`$ at the outermost mesh point |
| $`L`$ [L☉] | Luminosity emerging at the surface |
| $`T_{\mathrm{eff}}`$ [K], $`\log g`$ [cgs] | Effective temperature and surface gravity (see [below](#effective-temperature-and-surface-gravity)) |
| $`T_c`$, $`\rho_c`$, $`P_c`$ | Temperature, density and pressure at the innermost mesh point ($`m = 10^{-6} M`$), which are effectively the central values |
| Convective core [$`m/M`$] | Mass fraction of the convective region around the centre, or "none" if the centre is radiative. Its radius is shown in the views. |
| Conv. envelope base [$`r/R`$] | Inner edge of the outer convection zone, or "none" if the envelope is radiative. It is the zone that reaches within $`10^{-6}`$ of the total mass from the surface. |
| pp / CNO / 3α [%] | Share of the total nuclear energy production from the pp chain, the CNO cycle and helium burning |
| Virial check | Should be 1 (see [below](#virial-check)) |
| Energy check | Should be 1 (see [below](#energy-check)) |
| Newton iterations | Iterations in the last solver step (a continuation may take several steps) |
| Max residual | Largest residual of the discretised equations (see [below](#residuals)) |
| Age | Time since the ZAMS (0 unless the star has been evolved) |
| Phase | Evolutionary phase, see [Evolution in time](12-evolution.md#evolutionary-phases) |
| $`X_c`$ / He core [$`m/M`$] | Central hydrogen mass fraction, and the mass fraction of the hydrogen-exhausted core ($`X < 10^{-3}`$) |
| Fund. radial mode | Period of the fundamental radial pulsation mode, see [Pulsations & oscillations](13-pulsations.md) |
| $`\nu_{\max}`$ / $`\Delta\nu`$ [μHz] | Solar-like oscillation frequencies from the asteroseismic scaling relations (meaningful for stars with convective envelopes) |

## Effective temperature and surface gravity

```math
L = 4\pi R^2 \sigma T_{\mathrm{eff}}^4, \qquad \log g = \log_{10} \frac{G M}{R^2}
```

## Virial check

For any star in hydrostatic equilibrium the virial theorem makes this ratio exactly 1. Any deviation
measures the discretisation error, typically $`10^{-5}`$ to $`10^{-4}`$.

```math
\frac{\displaystyle\int_0^M \frac{3P}{\rho}\, dm}{\displaystyle\int_0^M \frac{G m}{r}\, dm} = 1 \qquad \text{(virial theorem)}
```

## Energy check

The energy generated inside (nuclear plus gravitational) must equal the energy leaving the surface.
$`L_0`$ is the luminosity of the innermost point.

```math
\frac{\displaystyle\int_{m_0}^M \varepsilon\, dm + L_0}{L} = 1 \qquad \text{(energy conservation)}
```

The integral uses the solver's own midpoint rule. It therefore integrates the same difference equation
the solver satisfies, and holds by construction for a converged model: it is a **discretisation
consistency check**, not an independent test of energy conservation. An evolution run has an
independent energy budget, described in [Evolution in time](12-evolution.md#checks).

## Residuals

The largest residual of each discretised equation (continuity, hydrostatic equilibrium, energy,
transport, and the boundary conditions), evaluated for the converged model. They are dimensionless and
typically about $`10^{-14}`$. Hover over the checks line in the app to see each one; the CLI prints
them, and the JSON export has them under `globals.residuals`.

> The virial, energy and residual checks measure **numerical** consistency. They say nothing about how
> close the model is to a real star; see [Accuracy & limitations](15-limitations.md).

## "Sun today" reference column

Present-day solar values, for orientation:

* **IAU 2015 nominal values:** $`R = 1\,R_\odot`$, $`L = 1\,L_\odot`$, $`T_{\mathrm{eff}} = 5772`$ K,
  $`\log g = 4.438`$.
* **Standard solar model centre:**
  * $`T_c \approx 1.57 \times 10^{7}`$ K
  * $`\rho_c \approx 150`$ g cm⁻³
  * $`P_c \approx 2.4 \times 10^{17}`$ dyn cm⁻²
* **Base of the convective envelope**, from helioseismology: $`0.713\,R`$.
* **Energy split:** about 99% pp, 1% CNO.

The present Sun has burnt about half of its central hydrogen. It is therefore brighter, larger and more
centrally condensed than a ZAMS model of the same mass.

> In the CLI the same quantities are printed by `stellar-cli` and written by `--json`.

---
← Previous: [Input parameters](03-parameters.md) · [Contents](README.md) · Next: [Profile quantities & export](05-profiles-and-export.md) →
