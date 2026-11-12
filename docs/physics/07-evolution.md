# 7. Evolution in time

## Implemented

Source: `crates/stellar-core/src/structure/evolve.rs`, `structure/equations.rs`.

Evolution is a sequence of structure models on a fixed mesh in $`q = m/M`$. Each time step
$`\Delta t`$ has two parts.

### Step 1: composition

Burning uses the rates of the previous model (an explicit step), with the equations of
[Nuclear burning](04-nuclear-burning.md):

```math
X^{\mathrm{new}} = X - \frac{\varepsilon_{pp} + \varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{H}}}\,\Delta t, \qquad Y^{\mathrm{new}} = Y + \frac{\varepsilon_{pp} + \varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{H}}}\,\Delta t - \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}}\,\Delta t
```

Each point burns at its own rate and the amount is limited to what the point contains. Then every
convective zone of the previous model is mixed to its mass average
([Convection and mixing](05-convection-and-mixing.md)). A convective core that shrinks leaves a
hydrogen gradient behind.

### Step 2: structure with gravitational energy

The structure equations are solved again for the new composition. The energy equation contains the
heat released or absorbed as each layer contracts or expands:

```math
\varepsilon_{\mathrm{grav}} = -T\frac{ds}{dt} = -c_P T\left(\frac{d\ln T}{dt} - \nabla_{\mathrm{ad}}\frac{d\ln P}{dt}\right) \approx -\frac{c_P T}{\Delta t}\left(\Delta\ln T - \nabla_{\mathrm{ad}}\,\Delta\ln P\right)
```

$`\Delta`$ is the change against the previous model at the same mass. The term is evaluated
implicitly (with the new $`c_P`$, $`T`$, $`\nabla_{\mathrm{ad}}`$), so contraction
($`\varepsilon_{\mathrm{grav}} > 0`$) and expansion ($`\varepsilon_{\mathrm{grav}} < 0`$) are both
handled. The contribution of the composition change to $`ds`$ is neglected.

### Time-step control

*Burning limit.* No shell may change its hydrogen by more than the following in one step; a
convective zone counts as one shell with its mass-averaged rate, and shells with
$`X \le 10^{-8}`$ are ignored.

```math
|\Delta X| \le \min\!\left(0.01,\ \max(0.3\,X,\ 10^{-5})\right)
```

*Structural limit.* The next step is scaled so that the largest change of $`\ln r`$, $`\ln P`$,
$`\ln T`$ (anywhere) and $`\ln L`$ (surface) is about 0.02:

```math
\Delta t_{\mathrm{next}} = \Delta t \cdot \min\!\left(2,\ \max\!\left(\tfrac{1}{2},\ \frac{0.02}{\max|\Delta \ln r, \Delta\ln P, \Delta\ln T, \Delta\ln L|}\right)\right)
```

* A step that changes the structure by more than 0.08 is repeated with half the time step.
* If Newton fails (limit 40 iterations), the step is divided by 3, up to 12 attempts; the smallest
  step is 1 year.
* The first step is 10% of the burning limit, at most $`10^6`$ years.

### Phases and stop conditions

| Phase | Condition |
|---|---|
| Zero-age main sequence | age 0 |
| Main sequence | $`X_c \ge 10^{-3}`$ |
| H-shell burning | $`X_c < 10^{-3}`$ |
| Red giant | shell burning, $`R > 5\,R_{\mathrm{ZAMS}}`$ and $`T_{\mathrm{eff}} < 5000`$ K |
| Core helium burning | triple-α gives more than 10% of the nuclear energy |

A run stops at the chosen condition (end of the main sequence, a given age, or as far as possible),
or when $`T_{\mathrm{eff}} < 2600`$ K, $`T_c > 3\times10^8`$ K, after 3000 steps, or at 30 Gyr.

### Fuel accounting (implemented check)

The code compares the hydrogen burnt with the energy radiated:

```math
Q_{\mathrm{H}}\left[\int X\, dm\right]_{t}^{0} \approx \int_0^t L\, dt'
```

This is a *fuel* check, not energy conservation. It is close to 1 only on the main sequence, where
hydrogen burning supplies nearly all the luminosity; contraction, expansion and helium burning break
it by design. What it does test exactly is the burning step itself:
$`Q_{\mathrm{H}}\,[\int X\,dm]^0_t = \int\!\!\int (\varepsilon_{pp} + \varepsilon_{\mathrm{CNO}})\,dm\,dt'`$.

**[measured]** ratio of the two sides: 0.9947 (1 M☉ to 4.57 Gyr); 0.9984, 0.9991, 0.9992 (1.5, 3 and
5 M☉ to the end of the main sequence). The difference from 1 mixes the change of thermal and
gravitational energy with numerical error; the check cannot separate them.

### Total energy conservation (not implemented as a check)

For constant mass and hydrostatic equilibrium:

```math
\int_0^t L\,dt' = \int_0^t\!\!\int_0^M \left(\varepsilon_{\mathrm{nuc}} - \varepsilon_\nu\right) dm\,dt' \;-\; \Delta\!\left(E_{\mathrm{int}} + E_{\mathrm{grav}}\right), \qquad E_{\mathrm{int}} = \int_0^M u\,dm, \quad E_{\mathrm{grav}} = -\int_0^M \frac{G m}{r}\,dm
```

Conventions that keep every term counted once:

* $`u`$ is the thermal energy only (gas, radiation and, later, ionisation and degeneracy energy). It
  contains no rest mass and no nuclear binding energy; those enter only through
  $`\varepsilon_{\mathrm{nuc}}`$.
* $`\varepsilon_{\mathrm{nuc}}`$ is the heat given to the gas, with reaction neutrinos already
  removed in the $`Q`$ values. $`\varepsilon_\nu`$ is the thermal-neutrino loss only (zero today).
* $`u`$ depends on composition through $`\mu`$. $`\Delta E_{\mathrm{int}}`$ computed from the two
  models includes that change; the implemented $`\varepsilon_{\mathrm{grav}}`$ does not (it drops
  the composition term of $`ds`$). The mismatch is a known residual of this equation, small on the
  main sequence.
* The work done at the surface, $`4\pi R^2 P_{\mathrm{ph}}\,dR/dt`$, is omitted; it is negligible
  for a photospheric boundary.

### Limits

* Results after the main sequence are qualitative. Without degeneracy and with the analytic opacity,
  low-mass stars do not climb the red giant branch.
* The mesh is fixed in $`q`$; a thin burning shell is resolved only by raising the number of points.

## Planned

### Burn a mixed zone as one unit (C7, plan 1.1)

Average the rate over the zone first, then burn and limit:

```math
\langle \dot X \rangle_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} \dot X\, dm}{\int_{\mathrm{zone}} dm}, \qquad X_{\mathrm{zone}}^{\mathrm{new}} = \max\!\left(0,\; \bar X_{\mathrm{zone}} + \langle \dot X \rangle_{\mathrm{zone}}\,\Delta t\right)
```

This removes a small inconsistency between the burning step and the time-step limit near hydrogen
exhaustion (effect at most 0.1% today).

### Helium time-step limit (C8, plan 1.2)

Wherever $`\varepsilon_{3\alpha} > 0`$:

```math
|\Delta Y| \le \min\!\left(0.01,\ \max(0.3\,Y,\ 10^{-5})\right)
```

### Mass loss (E1, plan phase 4)

*Hot stars* (de Jager, Nieuwenhuijzen & van der Hucht 1988), in M☉ yr⁻¹, with a metallicity factor
(Vink, de Koter & Lamers 2001):

```math
\log\left(-\dot M\right) = -8.158 + 1.769\,\log\frac{L}{L_\odot} - 1.676\,\log T_{\mathrm{eff}}, \qquad \dot M \propto \left(\frac{Z}{0.019}\right)^{0.69}
```

*Source comparison.* The three coefficients agree with de Jager et al. (1988), §4, remark 3,
p. 264. That linear form is the paper's own reduction of its fit to the linear terms; the fit
itself is a 20-term Chebyshev expansion in $`\log T_{\mathrm{eff}}`$ and $`\log(L/L_\odot)`$
(eq. 1, p. 262, normalised form p. 263, coefficients in Table VI), valid for
$`3.3 \le \log T_{\mathrm{eff}} \le 4.8`$ and $`2.5 \le \log(L/L_\odot) \le 6.7`$ (p. 263), and
covering O through M stars, not only hot ones. The linear form drops the higher terms; at the
paper's central point ($`\log T_{\mathrm{eff}} = 4.05`$, $`\log L/L_\odot = 4.6`$) it gives
$`\log(-\dot M) = -6.81`$, the value the paper quotes there (p. 264); its departure from the full
fit elsewhere is not quantified in the paper. The unit is not stated on the pages read; the values
of Table I (e.g. 4.90 for $`-\log(-\dot M)`$ of an O3 If star, p. 270) are in M☉ yr⁻¹.

The metallicity factor was compared with Vink et al. (2001), §7, eqs. 21–22 and the text after
them, and the abstract: $`\dot M \propto Z^{0.85}\,v_\infty^{p}`$ at fixed $`v_\infty/v_{\mathrm{esc}}`$;
with $`v_\infty \propto Z^{0.13}`$ this becomes $`\dot M \propto Z^{0.69}`$ for
$`T_{\mathrm{eff}} \gtrsim 25\,000`$ K and $`Z^{0.64}`$ for B supergiants below, for
$`1/30 \le Z/Z_\odot \le 3`$, with $`Z_\odot = 0.019`$ (§4). *Changed:* the exponent was 0.7 and the
reference metallicity 0.02; they are now the source's 0.69 and 0.019 (O-star side). Multiplying the
de Jager fit by Vink's factor is a combination made here, not a prescription of either paper; the
Vink factor belongs to Vink's own recipe (eq. 24, which has $`+0.85\log(Z/Z_\odot)`$ at fixed
$`v_\infty/v_{\mathrm{esc}}`$).

*Cool giants* (Reimers 1975), in solar units, default $`\eta = 0.477`$:

```math
\dot M = -4\times10^{-13}\,\eta\,\frac{L\,R}{M} \qquad [M_\odot\,\mathrm{yr}^{-1}]
```

The original (Reimers 1975, Mém. Soc. R. Sci. Liège, 6e sér., 8, 369) was not obtained. The
coefficient $`4\times10^{-13}`$ for $`L/(gR)`$ in solar units, i.e. $`L R / M`$, is confirmed by
Reimers (1977, A&A 61, 217), p. 217, who quotes it as his 1975 proposal and recalibrates it to
$`1.4\times10^{-13}`$ (p. 223, i.e. $`\eta \approx 0.35`$ in the form above). The factor $`\eta`$
does not appear there.

*Decided 2026-11-12:* $`\eta`$ is a user parameter with default 0.477, the median over 56 Galactic
globular clusters of McDonald & Zijlstra (2015, MNRAS 448, 502; arXiv:1501.00874 read 2026-11-12),
$`\eta_R = 0.477 \pm 0.070`$ (standard deviation) $`^{+0.050}_{-0.062}`$ (systematic), abstract and
§2.2 eq. 1, which is the formula above with the same $`4\times10^{-13}`$. They find $`\eta`$ varies
by less than about 30% over a factor 200 in iron abundance, so no metallicity factor is applied to
the cool-star branch. The calibration is for low-mass red giants; for other stars it is an
extrapolation, stated in the UI tooltip. The original Reimers (1975) paper was not read; its form
is confirmed by Reimers (1977) and by McDonald & Zijlstra.

Estimated loss on the main sequence from the ZAMS models: 0.003 M☉ (5 M☉), about 0.2 M☉ (20 M☉),
5–10 M☉ (60 M☉).

*Moving mesh.* With $`m_k = q_k\,M(t)`$, a quantity $`f`$ followed at fixed mass obeys

```math
\left(\frac{\partial f}{\partial t}\right)_{m} = \left(\frac{\partial f}{\partial t}\right)_{q} - q\,\frac{\dot M}{M}\,\left(\frac{\partial f}{\partial q}\right)_{t}
```

This applies to the composition **and** to the thermal variables. Comparing $`T`$ and $`P`$ at
fixed $`q`$ in $`\varepsilon_{\mathrm{grav}}`$ would drop a term
$`q\,(\dot M/M)\,T\,\partial s/\partial q`$. That term is a rate: a shorter time step does not make
it smaller, so it cannot be controlled by a step limit. With the transformation:

```math
\varepsilon_{\mathrm{grav}} = -\frac{c_P T}{\Delta t}\left[\left(\Delta_q \ln T - \nabla_{\mathrm{ad}}\,\Delta_q \ln P\right) - q\,\frac{\Delta M}{M}\left(\frac{\partial\ln T}{\partial q} - \nabla_{\mathrm{ad}}\,\frac{\partial\ln P}{\partial q}\right)\right]
```

$`\Delta_q`$ is the change at fixed $`q`$; the $`q`$-derivatives are taken on the new model. In
practice the previous model's $`\ln T`$, $`\ln P`$ and composition are interpolated to the new
$`m_k`$ before the step, which is the same thing. A limit on the mass lost per step is still used,
for the accuracy of that interpolation near the surface:

```math
\frac{|\Delta M|}{M} \le 10^{-3} \quad \text{per step}
```

In the outermost layers, which adjust thermally much faster than they are removed, the advection
term nearly cancels the time derivative; this cancellation is the main numerical difficulty of the
task.

*Fuel accounting with a wind*, $`X_s`$ being the surface hydrogen fraction (the total-energy
equation gains the thermal and gravitational energy carried off by the lost mass):

```math
\left[\int X\,dm\right]_{t}^{0} = \int_0^t\!\!\int_0^{M(t')} \left(\frac{\varepsilon_{pp}}{Q_{pp}^{\mathrm{eff}}} + \frac{\varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{CNO}}}\right) dm\,dt' + \int_0^t X_s\,|\dot M|\,dt'
```

### Pre-main-sequence contraction (F8, plan 9.4)

Start from a large, fully convective model (an $`n = 1.5`$ polytrope) and evolve with
$`\varepsilon_{\mathrm{grav}}`$ as the only energy source until hydrogen ignites:

```math
\rho_c = 5.991\,\bar\rho, \qquad P_c = 0.770\,\frac{G M^{2}}{R^{4}}, \qquad t_{\mathrm{KH}} \approx \frac{G M^{2}}{R\,L}
```

$`t_{\mathrm{KH}}`$ sets the first time step. The track is qualitative until the atmosphere and
opacity upgrades are in.

### Composition term in $`\varepsilon_{\mathrm{grav}}`$ (D17)

Listed as neglected; small on the main sequence. No equation is planned yet.

---
Previous: [Boundary conditions](06-boundary-conditions.md) · [Contents](README.md) · Next: [Rotation](08-rotation.md)
