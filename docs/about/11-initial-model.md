# 11. Initial model & continuation

## Cold start: an $`n = 3`$ polytrope

The Lane–Emden equation for $`n = 3`$ (the Eddington standard model) is integrated with 4th-order
Runge–Kutta, with step $`10^{-3}`$:

```math
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) = -\theta^n, \qquad \theta(0) = 1,\; \theta'(0) = 0, \qquad n = 3
```

The polytrope is scaled to homology estimates of the ZAMS radius and luminosity ($`M`$ in solar
masses):

```math
\frac{R}{R_\odot} \approx \begin{cases} 0.89\, M^{0.8} & M < 1.3 \\ 0.89 \cdot 1.3^{0.8}\, (M/1.3)^{0.6} & M \ge 1.3 \end{cases} \qquad \frac{L}{L_\odot} \approx \begin{cases} 0.7\, M^{4.5} & M < 2 \\ 0.7 \cdot 2^{4.5}\, (M/2)^{3.5} & M \ge 2 \end{cases}
```

Central values and the polytropic structure follow. $`T`$ comes from inverting the EOS at each point:

```math
\rho_c = 54.18\, \bar\rho, \qquad P_c = 11.05\, \frac{G M^2}{R^4}, \qquad \rho = \rho_c\, \theta^3, \quad P = P_c\, \theta^4
```

The outer layers are patched with a plane-parallel hydrostatic envelope; the larger of the two
pressures is used:

```math
P_{\mathrm{env}}(m) = \frac{G M\, (M - m)}{4\pi R^4} + P_{\mathrm{ph}}
```

For $`L(m)`$, $`\varepsilon`$ is integrated over the guessed structure and rescaled to the estimated
luminosity.

## Warm start and continuation

When parameters change, the solver starts from the model already on screen, scaled by homology. It is
interpolated onto the new mesh if $`N`$ changed:

```math
R \propto M^{0.75}, \qquad P \propto M^{-1}, \qquad T \propto M^{0.25}, \qquad L \propto M^{3.5}
```

* The solver first tries to reach the new parameters in a single step.
* If Newton fails, it halves the step along the path, down to 1/1024 of the path. The path is
  geometric in $`M`$ and linear in $`X`$, $`Z`$, $`\alpha_{\mathrm{MLT}}`$ and $`f_{\mathrm{CNO}}`$.
* It doubles the step again after each success.

If that fails too, the fallbacks are, in order:

1. a cold start from the polytrope;
2. a 1 M☉ reference model with the same composition, followed by continuation to the target;
3. a 1 M☉ reference model with the solar composition, followed by continuation to the target.

## Background solving

Solves run in a background thread, so the interface stays responsive.

* A new request cancels the one in progress, and the results of superseded requests are discarded.
* Edits are debounced by 0.15 s.
* The HR-diagram ZAMS line is computed in its own thread.

> The solver exposes single Newton steps (`HenyeyNewton::step`), so a future web build can iterate a
> few steps per frame (see [the web version notes](../web-version.md)).

---
← Previous: [Numerical method](10-numerical-method.md) · [Contents](README.md) · Next: [Evolution in time](12-evolution.md) →
