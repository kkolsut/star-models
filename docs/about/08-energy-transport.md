# 8. Energy transport & convection

## Radiative transport

If all the energy is carried by photon diffusion, the temperature gradient must be the *radiative
gradient*:

```math
\nabla_{\mathrm{rad}} = \frac{3\,\kappa\, L\, P}{16\pi\, a\, c\, G\, m\, T^4}
```

## Convective instability: Schwarzschild criterion

```math
\nabla_{\mathrm{rad}} > \nabla_{\mathrm{ad}} \;\Rightarrow\; \text{convective}, \qquad \nabla_{\mathrm{rad}} \le \nabla_{\mathrm{ad}} \;\Rightarrow\; \text{radiative},\; \nabla = \nabla_{\mathrm{rad}}
```

In a convective layer a rising blob that expands adiabatically stays hotter, and so lighter, than its
surroundings, and keeps rising. The Schwarzschild criterion is used throughout, also in evolved
models that have composition gradients: their stabilising effect (Ledoux criterion, semiconvection) is
neglected.

## Convective transport: mixing-length theory (MLT)

The model follows Böhm-Vitense (1958) in the formulation of Kippenhahn, Weigert & Weiss (2012, §7).

*Gravity, pressure scale height, mixing length*

```math
g = \frac{G m}{r^2}, \qquad H_P = \frac{P}{\rho\, g}, \qquad \ell = \alpha_{\mathrm{MLT}}\, H_P
```

```math
U = \frac{3\, a\, c\, T^3}{c_P\, \rho^2\, \kappa\, \ell^2}\sqrt{\frac{8 H_P}{g\, \delta}}, \qquad W = \nabla_{\mathrm{rad}} - \nabla_{\mathrm{ad}}
```

*The MLT cubic*

```math
(\xi - U)^3 + \frac{8U}{9}\left(\xi^2 - U^2 - W\right) = 0
```

The cubic has exactly one root with $`\xi \ge U`$. It is solved for $`\eta = \xi - U`$ by Newton's
method, starting from the upper bound $`\sqrt{U^2 + W} - U`$, and converges monotonically to machine
precision. Then:

*Actual gradient*

```math
\nabla = \nabla_{\mathrm{ad}} + \eta\,(\eta + 2U), \qquad \eta = \xi - U \ge 0
```

*Convective velocity*

```math
v_{\mathrm{conv}} = \eta\, \ell \sqrt{\frac{g\, \delta}{8 H_P}}
```

*Convective flux fraction*

```math
\frac{F_{\mathrm{conv}}}{F} = \frac{\nabla_{\mathrm{rad}} - \nabla}{\nabla_{\mathrm{rad}}}
```

The two limits:

* **Deep interior:** convection is efficient ($`U \ll 1`$) and $`\nabla \approx \nabla_{\mathrm{ad}}`$.
  The layer is almost exactly adiabatic.
* **Just below the photosphere:** convection is inefficient ($`U \gg 1`$) and $`\nabla`$ rises towards
  $`\nabla_{\mathrm{rad}}`$. This superadiabatic layer is the only place where $`\alpha_{\mathrm{MLT}}`$
  matters.

> Not included: convective overshooting, semiconvection, turbulent pressure.

---
← Previous: [Microphysics](07-microphysics.md) · [Contents](README.md) · Next: [Boundary conditions](09-boundary-conditions.md) →
