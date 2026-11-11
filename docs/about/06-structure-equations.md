# 6. Structure equations

The model solves the classical equations of stellar structure for a spherical star in hydrostatic and
thermal equilibrium. The independent variable is the mass $`m`$ enclosed within radius $`r`$ (the
Lagrangian form). It runs from the centre ($`m = 0`$) to the surface ($`m = M`$).

*Mass conservation*

```math
\frac{dr}{dm} = \frac{1}{4\pi r^2 \rho}
```

*Hydrostatic equilibrium*

```math
\frac{dP}{dm} = -\frac{G m}{4\pi r^4}
```

*Energy conservation*

```math
\frac{dL}{dm} = \varepsilon = \varepsilon_{pp} + \varepsilon_{\mathrm{CNO}} + \varepsilon_{3\alpha}
```

*Energy transport*

```math
\frac{dT}{dm} = -\frac{G m T}{4\pi r^4 P}\,\nabla, \qquad \nabla \equiv \frac{d\ln T}{d\ln P}
```

The temperature gradient $`\nabla`$ comes from radiative diffusion or from convection (see
[Energy transport & convection](08-energy-transport.md)).

The solver works with $`\ln r`$, $`\ln P`$, $`\ln T`$ and $`L`$, which keeps quantities spanning many
orders of magnitude well conditioned. $`\rho`$, $`\kappa`$, $`\varepsilon`$ and $`\nabla`$ are
functions of $`(P, T)`$ and the composition, given by the [microphysics](07-microphysics.md).

## Assumptions

* ZAMS models are static: the gravitational energy term $`\varepsilon_{\mathrm{grav}} = -T\,ds/dt`$ is
  zero. During evolution $`\varepsilon_{\mathrm{grav}}`$ is added to $`\varepsilon`$ (see
  [Evolution in time](12-evolution.md)).
* ZAMS models are chemically homogeneous: $`X`$, $`Y`$, $`Z`$ are the same everywhere. Evolved models
  store a composition per mesh point, set by burning and convective mixing (see
  [Evolution in time](12-evolution.md)).
* No neutrino losses, rotation, magnetic fields, mass loss or element diffusion.

---
← Previous: [Profile quantities & export](05-profiles-and-export.md) · [Contents](README.md) · Next: [Microphysics](07-microphysics.md) →
