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

## Surface, optional: integrated atmosphere

*Physics › Atmosphere* offers an integrated atmosphere instead (plans 3.1–3.3). Hydrostatic
equilibrium $`dP/d\tau = g/\kappa`$ is integrated with the local opacity $`\kappa(\rho(P, T), T)`$
from $`\tau_0 = 10^{-6}`$ to a fitting point $`\tau_{\mathrm{fit}}`$ by 100 fixed RK4 steps in
$`\ln\tau`$, with $`T^4 = \tfrac34 T_{\mathrm{eff}}^4[\tau + q(\tau)]`$. The start is
$`P(\tau_0) = aT^4(0)/3 + g\tau_0/\kappa_0`$ (grey, LTE radiation pressure: adequate for cool stars).
With convection, $`P`$ and $`T`$ are integrated together and $`\nabla`$ comes from the mixing-length
cubic with the Hopf factor $`A = 1 + q'(\tau)`$ and radiative losses of optically thin elements
(Henyey, Vardya & Bodenheimer 1965).

* grey: $`q = 2/3`$, fit at $`\tau = 2/3`$ (the photosphere).
* Krishna Swamy (1966): $`q = 1.39 - 0.815\,e^{-2.54\tau} - 0.025\,e^{-30\tau}`$ on the Rosseland
  depth, a fit to the observed Sun; fit at $`\tau = 3.1`$ where $`|q'| < 10^{-3}`$ (the interior keeps
  the uncorrected gradient: approximate). Above $`T_{\mathrm{eff}} \approx 8000`$ K it blends into the
  grey relation.

The last mesh point is then the fitting point, a height $`\Delta r = \int d\tau/(\kappa\rho)`$ below
the photosphere; the model radius is $`R = r_N + \Delta r`$ and $`T_{\mathrm{eff}}`$ follows from
$`L = 4\pi R^2\sigma T_{\mathrm{eff}}^4`$. The atmosphere is plane parallel ($`g = GM/R^2`$; its mass is
neglected). Details: [`docs/physics/06`](../physics/06-boundary-conditions.md).

## Counting equations

$`2`$ central conditions $`+\ 2`$ surface conditions $`+\ 4(N - 1)`$ difference equations $`= 4N`$
equations, for the $`4N`$ unknowns.

---
← Previous: [Energy transport & convection](08-energy-transport.md) · [Contents](README.md) · Next: [Numerical method](10-numerical-method.md) →
