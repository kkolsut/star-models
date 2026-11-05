# 10. Numerical method

## Mesh

The mesh has $`N`$ points in $`q = m/M`$, from $`10^{-6}`$ to 1, uniformly spaced in the stretched
coordinate below. Its three terms place about 15% / 45% / 40% of the points in the centre / interior /
envelope:

```math
s(q) = a \ln q + b\, q - c \ln\!\left(1 - q + 10^{-12}\right), \qquad q = \frac{m}{M}
```

This gives:

* logarithmic spacing in $`q`$ near the centre;
* uniform spacing through the interior;
* logarithmic spacing in $`1 - q`$ towards the surface, down to about $`10^{-13}`$.

The outer layers contain very little mass, but the pressure there falls by many orders of magnitude,
so they need many points.

## Unknowns and discretisation

Each mesh point carries four unknowns. $`L_{\mathrm{scale}}`$ is the expected luminosity, which keeps
all unknowns of order 1:

```math
\mathbf{y}_k = \left(\ln r_k,\; \ln P_k,\; \ln T_k,\; L_k / L_{\mathrm{scale}}\right), \qquad k = 0, \dots, N-1
```

Between neighbouring points, the change of each unknown equals $`\Delta m`$ times its derivative
$`\mathbf{f}`$ (the right-hand sides of the structure equations) at the midpoint state. This scheme is
second-order accurate:

```math
\mathbf{y}_{k+1} - \mathbf{y}_k = \left(m_{k+1} - m_k\right)\, \mathbf{f}\!\left(\frac{m_k + m_{k+1}}{2},\; \frac{\mathbf{y}_k + \mathbf{y}_{k+1}}{2}\right)
```

## Henyey relaxation (Newton–Raphson on the whole star)

All $`4N`$ residuals $`\mathbf{F}(\mathbf{y})`$ (difference equations plus boundary conditions) are
driven to zero together by Newton's method (Henyey, Forbes & Gould 1964). This is the standard method
of stellar evolution codes, and it is also what the time evolution uses. Each iteration solves a linear
system and applies a damped correction:

```math
\mathsf{J}\, \delta\mathbf{y} = -\mathbf{F}(\mathbf{y}), \qquad \mathbf{y} \leftarrow \mathbf{y} + \lambda\, \delta\mathbf{y}, \qquad \lambda = \min\!\left(1,\; \min_j \frac{\Delta_j^{\max}}{\max_k |\delta y_{k,j}|}\right)
```

The damping limits $`\Delta_j^{\max}`$ are one per unknown:

| Unknown | $`\Delta_j^{\max}`$ |
|---|---|
| $`\ln r`$ | 0.2 |
| $`\ln P`$ | 0.5 |
| $`\ln T`$ | 0.2 |
| $`L/L_{\mathrm{scale}}`$ | 0.5 |

**Jacobian.** It is built by central finite differences, one mesh point at a time. Each equation
involves only neighbouring points, so no analytic derivatives of the physics are needed:

```math
\mathsf{J}_{ij} = \frac{\partial F_i}{\partial y_j} \approx \frac{F_i(\mathbf{y} + h\,\mathbf{e}_j) - F_i(\mathbf{y} - h\,\mathbf{e}_j)}{2h}, \qquad h = 10^{-6}
```

The step $`h`$ is absolute for $`\ln r`$, $`\ln P`$ and $`\ln T`$, and relative for $`L`$.

**Linear solve.** $`\mathsf{J}`$ is banded, with 5 sub- and 5 super-diagonals. It is solved by Gaussian
elimination with partial pivoting in banded storage, so the cost grows linearly with $`N`$.

**Convergence.** A model is converged when the largest correction is below $`10^{-8}`$ with an
undamped step ($`\lambda = 1`$).

* This typically takes 10–30 iterations, with quadratic convergence at the end.
* The limit is 250 iterations from a cold start and 60 from a warm start.
* A small correction alone does not prove the equations are satisfied, so the residuals of the
  converged solution are evaluated as well. They are dimensionless (differences of $`\ln r`$,
  $`\ln P`$, $`\ln T`$ and of $`L`$ in units of the luminosity scale). Converged models reach about
  $`10^{-14}`$; a model above $`10^{-6}`$ is rejected. The largest residual of each equation is reported
  with the model.

**Checks.** Every converged model is checked, and a model failing any check is rejected:

* all values finite;
* $`r`$ increasing outward;
* $`P`$ decreasing outward;
* $`L > 0`$;
* $`\beta`$ above the floor;
* every residual below $`10^{-6}`$.

## Post-processing

* **Profiles:** the full local physics (EOS, $`\kappa`$, $`\varepsilon`$, MLT) is evaluated at every
  mesh point.
* **Zone boundaries:** located where $`\nabla_{\mathrm{rad}} - \nabla_{\mathrm{ad}}`$ changes sign, by
  linear interpolation between mesh points.
* **Integrals** (energy split and virial check) use the trapezoidal rule in $`m`$. The energy check
  uses the solver's own midpoint rule (see [Results](04-results.md#energy-check)).

## Validation

The same generic solver reproduces the Lane–Emden polytrope zeros $`\xi_1 = \pi`$ ($`n = 1`$) and
$`\xi_1 = 6.89685`$ ($`n = 3`$). Automated tests also cover:

* EOS limits, MLT limits, and positivity of opacity and rates;
* the banded solver against a dense solver;
* mesh refinement, and warm start against cold start;
* qualitative trends along the ZAMS (convective cores, CNO dominance, monotonic $`L(M)`$).

> A solar model takes about 13 iterations and 40–60 ms from a cold start. After a slider change the
> warm start typically takes 20–100 ms.

---
← Previous: [Boundary conditions](09-boundary-conditions.md) · [Contents](README.md) · Next: [Initial model & continuation](11-initial-model.md) →
