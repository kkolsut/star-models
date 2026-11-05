# 7. Microphysics

## Equation of state: fully ionised ideal gas plus radiation

```math
P = P_{\mathrm{gas}} + P_{\mathrm{rad}} = \frac{\rho\, k\, T}{\mu\, m_u} + \frac{a T^4}{3}
```

*Mean molecular weight (full ionisation)*

```math
\frac{1}{\mu} = 2X + \frac{3}{4}\,Y + \frac{1}{2}\,Z
```

```math
\beta = \frac{P_{\mathrm{gas}}}{P}, \qquad \delta = -\left(\frac{\partial \ln\rho}{\partial \ln T}\right)_{\!P} = \frac{4 - 3\beta}{\beta}
```

```math
c_P = \frac{k}{\mu\, m_u}\left[\frac{3}{2} + \frac{12\,(1-\beta)}{\beta} + \frac{(4 - 3\beta)^2}{\beta^2}\right]
```

```math
\nabla_{\mathrm{ad}} = \left(\frac{\partial \ln T}{\partial \ln P}\right)_{\!s} = \frac{2\,(4 - 3\beta)}{32 - 24\beta - 3\beta^2}
```

$`\nabla_{\mathrm{ad}}`$ is 0.4 for a pure gas and 0.25 for pure radiation.

Given $`P`$ and $`T`$, the density follows directly from the first equation. During Newton iterations
$`\beta`$ is limited to at least $`10^{-3}`$, which keeps unphysical trial states finite. Converged
models that reach this floor are rejected.

## Opacity: analytic Rosseland mean (CGS)

*Electron scattering*

```math
\kappa_{\mathrm{es}} = 0.2\,(1 + X)
```

*Kramers free–free*

```math
\kappa_{\mathrm{ff}} = 3.68\times10^{22}\,(1 - Z)(1 + X)\,\rho\,T^{-3.5}
```

*Kramers bound–free*, where $`t`$ is the guillotine factor

```math
\kappa_{\mathrm{bf}} = \frac{4.34\times10^{25}}{t}\,Z\,(1 + X)\,\rho\,T^{-3.5}, \qquad t = 2.82\,\bigl[\rho\,(1 + X)\bigr]^{0.2}
```

*H⁻ ion* (cool envelopes)

```math
\kappa_{\mathrm{H}^-} = 2.5\times10^{-31}\,\frac{Z}{0.02}\,\rho^{1/2}\,T^{9}
```

*Molecular floor*

```math
\kappa_{\mathrm{mol}} = 0.1\,Z
```

*Combined opacity*

```math
\kappa = \kappa_{\mathrm{mol}} + \left[\frac{1}{\kappa_{\mathrm{H}^-}} + \frac{1}{\kappa_{\mathrm{es}} + \kappa_{\mathrm{ff}} + \kappa_{\mathrm{bf}}}\right]^{-1}
```

The harmonic combination makes H⁻ dominate only in cool layers. Without it the steep $`T^9`$ term would
run away at high temperature. For $`Z = 0`$, the H⁻ term uses $`Z = 10^{-4}`$. These coefficients
follow Carroll & Ostlie (2017), converted to CGS.

## Nuclear energy generation [erg g⁻¹ s⁻¹]

```math
T_6 = \frac{T}{10^6\,\mathrm{K}}, \qquad T_8 = \frac{T}{10^8\,\mathrm{K}}
```

*pp chain*

```math
\varepsilon_{pp} = 2.38\times10^{6}\;\rho\, X^2\, f_{1}\, \psi\; T_6^{-2/3}\, \exp\!\left(-33.80\, T_6^{-1/3}\right)
```

```math
\psi = 1 + 0.0123\, T_6^{1/3} + 0.0109\, T_6^{2/3} + 0.0009\, T_6
```

*CNO cycle*

```math
\varepsilon_{\mathrm{CNO}} = 8.67\times10^{27}\;\rho\, X\, X_{\mathrm{CNO}}\, f_{7}\, C\; T_6^{-2/3}\, \exp\!\left(-152.28\, T_6^{-1/3}\right)
```

```math
C = 1 + 0.0027\, T_6^{1/3} - 0.00778\, T_6^{2/3} - 0.000149\, T_6, \qquad X_{\mathrm{CNO}} = f_{\mathrm{CNO}}\, Z
```

*Triple-α*

```math
\varepsilon_{3\alpha} = 5.09\times10^{11}\;\rho^2\, Y^3\, f_{4}\, f_{8}\; T_8^{-3}\, \exp\!\left(-\frac{44.027}{T_8}\right)
```

*Weak electron screening* (Salpeter / Clayton)

```math
f_{Z_1 Z_2} = \exp\!\left(0.188\, Z_1 Z_2\, \zeta\, \rho^{1/2}\, T_6^{-3/2}\right), \qquad \zeta^2 = 2X + 1.5\,Y + 4Z
```

The charge products are:

| Reaction | $`Z_1 Z_2`$ |
|---|---|
| p + p | 1 |
| p + ¹⁴N (the slowest CNO step) | 7 |
| α + α | 4 |
| α + ⁸Be | 8 |

The screening exponent is capped at 2, since weak screening is only valid for small values.

> EOS, opacity and nuclear rates sit behind interfaces (`Eos`, `Opacity`, `Nuclear`), so tabulated
> physics can replace them without touching the solver.

---
← Previous: [Structure equations](06-structure-equations.md) · [Contents](README.md) · Next: [Energy transport & convection](08-energy-transport.md) →
