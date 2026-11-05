# 9. Boundary conditions

## Centre

The equations are singular at $`m = 0`$. The first mesh point therefore sits at $`m_0 = 10^{-6} M`$ and
uses the leading terms of the series expansion about the centre:

```math
r_0 = \left(\frac{3\, m_0}{4\pi\, \rho_0}\right)^{1/3}, \qquad L_0 = \varepsilon_0\, m_0, \qquad m_0 = 10^{-6} M
```

## Surface: grey Eddington atmosphere

```math
T^4(\tau) = \frac{3}{4}\, T_{\mathrm{eff}}^4 \left(\tau + \frac{2}{3}\right)
```

The outermost mesh point is the photosphere, at optical depth $`\tau = 2/3`$, where
$`T = T_{\mathrm{eff}}`$. Two conditions hold there:

```math
L = 4\pi R^2 \sigma T^4 \qquad (T = T_{\mathrm{eff}} \text{ at } \tau = \tfrac{2}{3})
```

```math
P = \frac{2}{3}\, \frac{G M}{R^2\, \kappa} + \frac{a\, T_{\mathrm{eff}}^4}{6}
```

The pressure condition comes from integrating hydrostatic equilibrium in optical depth,
$`dP/d\tau = g/\kappa`$, from $`\tau = 0`$, where only the radiation pressure
$`a\,T_{\mathrm{eff}}^4/6`$ remains. Here $`\kappa`$ is evaluated self-consistently at the
photospheric $`(P, T)`$.

## Counting equations

$`2`$ central conditions $`+\ 2`$ surface conditions $`+\ 4(N - 1)`$ difference equations $`= 4N`$
equations, for the $`4N`$ unknowns.

---
← Previous: [Energy transport & convection](08-energy-transport.md) · [Contents](README.md) · Next: [Numerical method](10-numerical-method.md) →
