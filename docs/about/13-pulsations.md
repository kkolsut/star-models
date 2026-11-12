# 13. Pulsations & oscillations

## Radial pulsation modes

The fundamental radial mode is computed for every model: the ZAMS model and every model of an evolution
run. It comes from the linear adiabatic wave equation, where $`\xi`$ is the relative radial
displacement and $`p`$ the relative pressure perturbation:

```math
\frac{d\xi}{dm} = -\frac{3\xi + p/\Gamma_1}{4\pi r^3\rho}, \qquad \frac{dp}{dm} = \frac{G m}{4\pi r^4 P}\left[\left(4 + \frac{\omega^2 r^3}{G m}\right)\xi + p\right], \qquad \xi = \frac{\delta r}{r},\; p = \frac{\delta P}{P}
```

```math
3\xi + \frac{p}{\Gamma_1} = 0 \;\; (\text{centre}), \qquad p + (4 + \Omega^2)\,\xi = 0,\;\; \xi = 1 \;\; (\text{surface}), \qquad \Omega^2 = \frac{\omega^2 R^3}{G M}
```

*First adiabatic exponent* (ideal gas + radiation)

```math
\Gamma_1 = \left(\frac{\partial\ln P}{\partial\ln\rho}\right)_{\!s} = \frac{32 - 24\beta - 3\beta^2}{24 - 21\beta}
```

* The dimensionless eigenvalue $`\Omega^2`$ is an extra unknown of the boundary-value problem.
* It is solved with the same Henyey method as the structure, with three unknowns per mesh point:
  $`(\xi, p, \Omega^2)`$.
* The fundamental mode is the solution without nodes in $`\xi`$.
* During evolution the previous model's mode is the starting guess, which keeps the solver on the
  fundamental branch.

Period and pulsation constant:

```math
\Pi = \frac{2\pi}{\omega}, \qquad Q = \Pi\,\sqrt{\frac{\bar\rho}{\bar\rho_\odot}}
```

**Validation.** For a uniform-density sphere the solver reproduces the exact result

```math
\Omega^2 = 3\Gamma_1 - 4 \qquad \text{(uniform-density sphere, } \xi = \text{const)}
```

and the ZAMS Sun's fundamental radial mode comes out at 0.99 hours. This mode is not observed in the
Sun; the standard solar model Model S (Christensen-Dalsgaard et al. 1996) gives 64.6 minutes
(258 μHz). The agreement is partly a coincidence, since the model Sun has $`R = 1.11\,R_\odot`$.

## What is calculated and what is not

| | |
|---|---|
| **Calculated** | the fundamental adiabatic radial mode: period, $`\Omega^2`$, $`Q`$ and the mode shape $`\xi(r)`$ |
| **Empirical inference** | membership of the instability strip, from the model's position in the HR diagram |
| **Scaling relations** | $`\nu_{\max}`$ and $`\Delta\nu`$, estimated from $`g`$, $`T_{\mathrm{eff}}`$ and the mean density |
| **Asymptotic estimates** | $`\Delta\nu`$ from the sound travel time, and $`\Delta\Pi_1`$ for every candidate g-mode cavity, from the profiles of $`c_s`$ and $`N`$ |
| **Not modelled** | driving and damping (growth rates), the excitation mechanism, amplitudes, overtones and non-radial mode frequencies |

δ Scuti stars in particular usually pulsate in several radial and non-radial modes at once, which a
fundamental radial-mode solver cannot represent.

## Which stars pulsate: the instability strip

Whether a mode is actually driven needs a non-adiabatic analysis, which is not done here. In real stars
the driving is the κ mechanism in the He II ionisation zone.

Instead, an approximate classical instability strip in the HR diagram decides it. Its edges are linear
fits to the observed strip:

* δ Scuti stars: 6300–8700 K at $`\log L \approx 0.5`$;
* classical Cepheids: about 5400–6800 K at $`\log L \approx 3.5`$.

```math
\log T_{\mathrm{blue}} = 3.94 - 0.037\,(\log L - 0.5), \qquad \log T_{\mathrm{red}} = 3.80 - 0.023\,(\log L - 0.5), \qquad \log L \ge 0.3
```

| Type | Where | Typical periods |
|---|---|---|
| **δ Scuti-type** | inside the strip with $`\log L < 1.7`$, on or near the main sequence | about 0.02–0.3 d |
| **Cepheid-type** | inside the strip with $`\log L \ge 1.7`$, for example while a 4–10 M☉ star crosses the Hertzsprung gap | days |

Inside the strip the 3D model pulsates with the computed fundamental mode:

* Each layer is displaced by $`\xi(r)`$: largest at the surface and small in the core.
* One cycle is shown per 1.6 s.
* The amplitude is exaggerated (6% of the radius).
* The surface is brightest near maximum compression.

The HR diagram shows the strip as a shaded band, and the timeline names the pulsator type and period.

## Solar-like oscillations

Stars with convective envelopes, like the Sun, oscillate in many p modes excited by turbulent
convection. Their frequency of maximum power and their large frequency separation follow the
asteroseismic scaling relations, shown in the results panel:

```math
\nu_{\max} = 3090\,\mu\mathrm{Hz}\;\frac{g}{g_\odot}\left(\frac{T_{\mathrm{eff}}}{T_{\mathrm{eff},\odot}}\right)^{-1/2}, \qquad \Delta\nu = 135.1\,\mu\mathrm{Hz}\,\sqrt{\frac{\bar\rho}{\bar\rho_\odot}}
```

## Characteristic frequencies and asymptotic spacings

Every model also carries the sound speed, the Lamb frequency and the buoyancy (Brunt–Väisälä)
frequency $`N`$:

```math
c_s^2 = \frac{\Gamma_1 P}{\rho}, \qquad S_l^2 = \frac{l(l+1)\,c_s^2}{r^2}
```

```math
N^2 = \frac{g^2 \rho}{P}\left[\delta\left(\nabla_{\mathrm{ad}} - \nabla\right) + \varphi\,\nabla_\mu\right], \qquad \nabla_\mu = \frac{d\ln\mu}{d\ln P}, \qquad \varphi = 1
```

* $`\nabla_\mu`$ is taken from the profile by differences; $`\varphi = 1`$ holds for the ideal-gas EOS.
* $`N^2 < 0`$ marks Ledoux instability. In a homogeneous model it coincides with the convective
  (Schwarzschild) flag; with a composition gradient a Schwarzschild-unstable layer can have
  $`N^2 > 0`$.
* The definition $`N^2 = g(\Gamma_1^{-1}\,d\ln P/dr - d\ln\rho/dr)`$ agrees with this formula to a
  median of $`2\times10^{-5}`$ in an evolved 3 M☉ model. It is a test only, because it differences
  nearly equal numbers.
* p modes propagate where the frequency is above both $`N`$ and $`S_l`$, g modes where it is below
  both. More profiles shows this propagation diagram.

```math
\Delta\nu = \left[2\int_0^R \frac{dr}{c_s}\right]^{-1}, \qquad \Delta\Pi_1 = \frac{2\pi^2}{\sqrt{2}}\left[\int_{r_1}^{r_2} N\,\frac{dr}{r}\right]^{-1}
```

* $`\Delta\nu`$ from the interior replaces the scaling relation by a value that depends on the
  structure. For the ZAMS Sun the two agree to 1% (116.0 against 115.0 μHz).
* $`\Delta\Pi_1`$ is computed separately for every connected $`N^2 > 0`$ region: a radiative core, the
  radiative zone between a convective core and envelope, or a radiative envelope. Regions are never
  summed. The value over the whole region is the low-frequency limit.
* Both are asymptotic estimates (high radial order), not computed spectra. $`\Delta\nu`$ carries the
  model's surface errors.

> Real Cepheids change their radius by a few percent, and δ Scuti stars by much less. Periods range from
> about half an hour to weeks. The animation is not to scale in time or amplitude.

---
← Previous: [Evolution in time](12-evolution.md) · [Contents](README.md) · Next: [How the views are made](14-views.md) →
