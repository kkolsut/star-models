# 10. Numerical method and checks

How the equations of pages 1 to 9 are solved, and the quantities derived from a converged model.

## Implemented

Source: `crates/stellar-core/src/structure/mesh.rs`, `solver/henyey.rs`, `solver/banded.rs`,
`structure/initial.rs`, `structure/solve.rs`, `structure/model.rs`.

### Mesh

$`N`$ points in $`q = m/M`$ from $`10^{-6}`$ to 1, uniform in the stretched coordinate

```math
s(q) = a \ln q + b\, q - c \ln\!\left(1 - q + 10^{-12}\right), \qquad a = \frac{0.15}{\ln 10^{6}}, \quad b = 0.45, \quad c = \frac{0.40}{\ln 10^{12}}
```

This puts about 15% of the points near the centre (logarithmic in $`q`$), 45% in the interior
(uniform in $`q`$) and 40% in the envelope (logarithmic in $`1 - q`$).

### Unknowns and discretisation

```math
\mathbf{y}_k = \left(\ln r_k,\; \ln P_k,\; \ln T_k,\; L_k / L_{\mathrm{scale}}\right), \qquad k = 0, \dots, N-1
```

Between neighbouring points the equations are applied at the midpoint state (second-order
accurate), with $`\mathbf{f}`$ the right-hand sides of the
[structure equations](01-structure-equations.md):

```math
\mathbf{y}_{k+1} - \mathbf{y}_k = \left(m_{k+1} - m_k\right)\, \mathbf{f}\!\left(\frac{m_k + m_{k+1}}{2},\; \frac{\mathbf{y}_k + \mathbf{y}_{k+1}}{2}\right)
```

Counting: 2 central conditions, 2 surface conditions and $`4(N-1)`$ difference equations for
$`4N`$ unknowns.

### Henyey relaxation

All residuals $`\mathbf{F}(\mathbf{y})`$ are driven to zero together by damped Newton iteration:

```math
\mathsf{J}\, \delta\mathbf{y} = -\mathbf{F}(\mathbf{y}), \qquad \mathbf{y} \leftarrow \mathbf{y} + \lambda\, \delta\mathbf{y}, \qquad \lambda = \min\!\left(1,\; \min_j \frac{\Delta_j^{\max}}{\max_k |\delta y_{k,j}|}\right)
```

| Unknown | $`\Delta_j^{\max}`$ |
|---|---|
| $`\ln r`$ | 0.2 |
| $`\ln P`$ | 0.5 |
| $`\ln T`$ | 0.2 |
| $`L/L_{\mathrm{scale}}`$ | 0.5 |

The Jacobian is built by central finite differences and is banded (5 sub- and 5 super-diagonals):

```math
\mathsf{J}_{ij} \approx \frac{F_i(\mathbf{y} + h\,\mathbf{e}_j) - F_i(\mathbf{y} - h\,\mathbf{e}_j)}{2h}, \qquad h = 10^{-6}
```

A model is converged when the largest correction is below $`10^{-8}`$ with $`\lambda = 1`$. Iteration
limits: 250 from a cold start, 60 from a warm start, 40 in an evolution step.

> **Consequence for every upgrade.** Because derivatives are taken numerically with
> $`h = 10^{-6}`$, any iteration inside the physics (root search for $`\rho`$, ionisation,
> atmosphere integration) must converge to about $`10^{-12}`$, and any interpolation must have
> continuous first derivatives.

### Starting model

*Cold start*: an $`n = 3`$ polytrope from the Lane–Emden equation,

```math
\frac{1}{\xi^2}\frac{d}{d\xi}\left(\xi^2 \frac{d\theta}{d\xi}\right) = -\theta^n, \qquad \theta(0) = 1,\; \theta'(0) = 0
```

```math
\rho_c = 54.18\, \bar\rho, \qquad P_c = 11.05\, \frac{G M^2}{R^4}, \qquad \rho = \rho_c\, \theta^3, \quad P = P_c\, \theta^4
```

scaled to rough estimates of the ZAMS radius and luminosity ($`M`$ in solar masses):

```math
\frac{R}{R_\odot} \approx \begin{cases} 0.89\, M^{0.8} & M < 1.3 \\ 0.89 \cdot 1.3^{0.8}\, (M/1.3)^{0.6} & M \ge 1.3 \end{cases} \qquad \frac{L}{L_\odot} \approx \begin{cases} 0.7\, M^{4.5} & M < 2 \\ 0.7 \cdot 2^{4.5}\, (M/2)^{3.5} & M \ge 2 \end{cases}
```

The outer layers are patched with a plane-parallel envelope:

```math
P_{\mathrm{env}}(m) = \frac{G M\, (M - m)}{4\pi R^4} + P_{\mathrm{ph}}
```

*Warm start*: the previous model scaled by homology, then continuation in the parameters with step
halving if Newton fails.

```math
R \propto M^{0.75}, \qquad P \propto M^{-1}, \qquad T \propto M^{0.25}, \qquad L \propto M^{3.5}
```

### Derived quantities

```math
L = 4\pi R^2 \sigma T_{\mathrm{eff}}^4, \qquad \log g = \log_{10} \frac{G M}{R^2}, \qquad \bar\rho = \frac{3 M}{4\pi R^3}
```

Convective-zone boundaries are located where $`\nabla_{\mathrm{rad}} - \nabla_{\mathrm{ad}}`$ changes
sign, by linear interpolation between mesh points.

### Self-checks

*Virial check.* Equal to 1 for hydrostatic equilibrium, up to the surface term
$`4\pi R^3 P_{\mathrm{ph}}`$ and the cut at $`m_0`$ (both negligible). Deviations of $`10^{-5}`$ to
$`10^{-4}`$ measure the discretisation error.

```math
\frac{\displaystyle\int_0^M \frac{3P}{\rho}\, dm}{\displaystyle\int_0^M \frac{G m}{r}\, dm} = 1
```

*Energy check*, with the solver's own midpoint rule. It holds by construction for a converged model,
so it tests the discretisation, not energy conservation over time (for that see the budget in
[Evolution in time](07-evolution.md)).

```math
\frac{\displaystyle\int_{m_0}^M \left(\varepsilon_{\mathrm{nuc}} + \varepsilon_{\mathrm{grav}}\right) dm + L_0}{L} = 1
```

*Gravitational redshift* (E2, plan 1.8). An output only:

```math
z = \frac{G M}{R\,c^{2}}, \qquad v_z = c\,z = \frac{G M}{R\,c}
```

636 m/s for 1 M☉ at 1 R☉ (test); results panel, compare tab, CLI, summary JSON
(`grav_redshift_m_s`) and track CSV. Relativistic corrections to the structure are $`2\times10^{-6}`$
(Sun) to $`1.2\times10^{-5}`$ (60 M☉), below the discretisation error, so they stay out of the solver;
About §15 lists general relativity under "Not modelled".

*Sanity checks*: all values finite, $`r`$ increasing, $`P`$ decreasing, $`L > 0`$, $`\beta`$ above
its floor.

*Validation of the solver*: Lane–Emden zeros $`\xi_1 = \pi`$ ($`n = 1`$) and $`\xi_1 = 6.89685`$
($`n = 3`$); banded solver against a dense one; mesh refinement (800 → 1600 points changes $`L`$,
$`R`$, $`T_c`$ of the Sun by less than 0.5%).

### Helium tied to metallicity (F5, plan 1.7)

An optional input mode, not a change of the physics: the checkbox "Helium follows metallicity" in the
app and `--helium-from-z` in the CLI set $`X`$ from $`Z`$; `stellar-core` is unchanged. The parameters
are those of the MIST models (Choi et al. 2016, §3.1, eqs. 1–3):

```math
Y = Y_p + \frac{\Delta Y}{\Delta Z}\,Z, \qquad Y_p = 0.249, \quad \frac{\Delta Y}{\Delta Z} = \frac{Y_{\odot,\mathrm{proto}} - Y_p}{Z_{\odot,\mathrm{proto}}} = 1.5, \qquad X = 1 - Y - Z
```

with the Asplund et al. (2009) protosolar $`Y = 0.2703`$, $`Z = 0.0142`$ (Choi et al. 2016,
Table 1), so that at $`Z = 0.0142`$ the mode reproduces the protosolar preset of the code
($`X = 0.7154`$) to within rounding (the mode gives 0.7155). *Changed:* the page gave $`Y_p \approx 0.248`$ without a source;
no source read gives that value together with 1.5. The pairs found are:

| Source | $`Y_p`$ | $`\Delta Y/\Delta Z`$ | Location |
|---|---|---|---|
| Choi et al. 2016 (MIST), from Planck 2015 | 0.249 | 1.5 | §3.1, eqs. 1–2 |
| Bressan et al. 2012 (PARSEC), $`Y_p`$ from Komatsu et al. 2011 | 0.2485 | 1.78 | §3 (from its solar calibration), used in §4 |
| Planck 2018 VI, BBN (PArthENoPE) at the Planck baryon density | 0.2454 | – | §2.1; Table 2 caption |
| Aver et al. 2021, H II regions (measured) | 0.2453 ± 0.0034 | – | abstract |

Current determinations of $`Y_p`$ are lower than 0.249 by about 0.004. With $`Y_p = 0.2453`$ and the
same protosolar anchor the slope would be $`(0.2703 - 0.2453)/0.0142 = 1.76`$; the two pairs agree
at $`Z = 0.0142`$ and differ by at most 0.0037 in $`Y`$ (at $`Z = 0`$) for $`Z \le 0.0142`$. Either pair may be
offered; the default follows MIST so that both numbers come from one source.

## Planned

### Virial check with rotation (E3)

See [Rotation](08-rotation.md).

---
Previous: [Pulsations](09-pulsations.md) · [Contents](README.md) · Next: [Validity, limits and order of work](11-validity-and-roadmap.md)
