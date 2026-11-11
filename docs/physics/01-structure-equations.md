# 1. Structure equations

## Implemented

Source: `crates/stellar-core/src/structure/equations.rs`.

A star is treated as a sphere in hydrostatic equilibrium. The independent variable is the mass $`m`$
enclosed within radius $`r`$, from the centre ($`m = 0`$) to the surface ($`m = M`$).

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
\frac{dL}{dm} = \varepsilon_{\mathrm{nuc}} + \varepsilon_{\mathrm{grav}}, \qquad \varepsilon_{\mathrm{nuc}} = \varepsilon_{pp} + \varepsilon_{\mathrm{CNO}} + \varepsilon_{3\alpha}
```

*Energy transport*

```math
\frac{dT}{dm} = -\frac{G m T}{4\pi r^4 P}\,\nabla, \qquad \nabla \equiv \frac{d\ln T}{d\ln P}
```

The closure relations $`\rho(P,T)`$, $`\kappa(\rho,T)`$, $`\varepsilon_{\mathrm{nuc}}(\rho,T)`$ and
$`\nabla`$ are given in pages 2 to 5. $`\varepsilon_{\mathrm{grav}}`$ is zero for a static (ZAMS)
model and is defined in [Evolution in time](07-evolution.md).

### Form used by the solver

The unknowns are $`\ln r`$, $`\ln P`$, $`\ln T`$ and $`L`$:

```math
\frac{d\ln r}{dm} = \frac{1}{4\pi r^3 \rho}, \qquad \frac{d\ln P}{dm} = -\frac{G m}{4\pi r^4 P}, \qquad \frac{d\ln T}{dm} = \nabla\,\frac{d\ln P}{dm}, \qquad \frac{dL}{dm} = \varepsilon_{\mathrm{nuc}} + \varepsilon_{\mathrm{grav}}
```

### Assumptions

* Spherical symmetry; no rotation and no magnetic field.
* Hydrostatic equilibrium at all times (no acceleration term).
* Constant total mass.
* A ZAMS model is chemically homogeneous and in thermal equilibrium
  ($`\varepsilon_{\mathrm{grav}} = 0`$). An evolved model has a composition per mesh point.
* No energy loss by neutrinos other than what is already subtracted in the reaction energies.
* Newtonian gravity.

## Planned

### Thermal-neutrino losses (F7, plan 9.6, done with phases 6–7)

```math
\frac{dL}{dm} = \varepsilon_{\mathrm{nuc}} + \varepsilon_{\mathrm{grav}} - \varepsilon_\nu
```

$`\varepsilon_\nu(\rho, T, \text{composition})`$ is the loss by *thermal* neutrinos (plasma, photo,
pair, bremsstrahlung), from the fits of Itoh, Hayashi, Nishikawa & Kohyama (1996, ApJS 102, 411–424;
original read 2026-11-11). It is separate from the *reaction* neutrinos, which are already
subtracted in the $`Q`$ values of [Nuclear burning](04-nuclear-burning.md) and must not be counted
twice.

What the paper gives, for whoever codes it (no coefficients are copied here; take them from the
paper):

| Process | Location in Itoh et al. (1996) | Fit valid for |
|---|---|---|
| pair | §2, eqs. 2.1–2.10, Table 1 (fit of Itoh et al. 1989) | $`10^{0} \le \rho/\mu_e \le 10^{14}`$ g cm⁻³, $`10^{7} \le T \le 10^{11}`$ K |
| photo | §3, eqs. 3.1–3.8, Table 2 (fit of Itoh et al. 1989) | $`10^{0} \le \rho/\mu_e \le 10^{11}`$ g cm⁻³, $`10^{7} \le T \le 10^{11}`$ K |
| plasma | §4, eqs. 4.1–4.11 (the fit of Haft, Raffelt & Weiss 1994, their eqs. 23–27, reproduced) | range not stated in the paper |
| bremsstrahlung, gas | §5.1, eqs. 5.1–5.15 | $`10^{0} \le \rho/\mu_e \le 10^{12}`$ g cm⁻³, $`10^{8} \le T \le 10^{11}`$ K |
| bremsstrahlung, liquid | §5.2–5.3, eqs. 5.16–5.28, Tables 3–4 | $`10^{4} \le \rho \le 4.3\times10^{11}`$ g cm⁻³, $`1 \le \Gamma \le 160`$ |
| bremsstrahlung, crystal ($`\Gamma \ge 180`$) | later in §5 (lattice and phonon parts) | not recorded here |
| recombination | §6, Table 12 | not recorded here |

The rates are per unit volume (erg cm⁻³ s⁻¹, §1), so $`\varepsilon_\nu`$ is the sum divided by
$`\rho`$. The stated accuracy is about 10% (pair), 1% (photo) and 5% (plasma) where that process
dominates; it is poor outside those regions, where the paper recommends its numerical tables
(CD-ROM) instead (§§2–3). "Applied wherever they are valid" in plan 9.6 refers to these ranges; what
to do outside them is still to be decided.

There is no single temperature threshold: the rate depends on density and composition as well.
Plasma neutrinos already cool the helium cores of red giants near $`10^8`$ K, where they shift the
temperature maximum off-centre and set the core mass at helium ignition. The term therefore belongs
with the red-giant work (degenerate EOS and conduction, phases 6 and 7), not after it.

### Rotation (E3, plan 5.2)

A spherically averaged centrifugal force adds one factor to hydrostatic equilibrium; see
[Rotation](08-rotation.md):

```math
\frac{dP}{dm} = -\frac{G m}{4\pi r^4}\, f_P, \qquad f_P = 1 - \frac{2\,\Omega^2 r^3}{3\,G m}
```

### Variable mass (E1, plan phase 4)

With mass loss, $`M`$ depends on time and the mesh points $`m_k = q_k M(t)`$ no longer follow the
same gas; see [Evolution in time](07-evolution.md).

### Not planned: general relativity (E2)

The relativistic form of hydrostatic equilibrium multiplies $`dP/dr`$ by

```math
\left(1 + \frac{P}{\rho c^2}\right)\left(1 + \frac{4\pi r^3 P}{m c^2}\right)\left(1 - \frac{2 G m}{r c^2}\right)^{-1}
```

Each factor differs from 1 by about $`GM/(Rc^2)`$: $`2.1\times10^{-6}`$ for the Sun and
$`1.2\times10^{-5}`$ for the 60 M☉ ZAMS model. This is below the discretisation error of the solver
($`10^{-5}`$ to $`10^{-4}`$), so it stays out. Only the surface redshift is planned as an output; see
[Numerical method and checks](10-numerics-and-checks.md).

---
[Contents](README.md) · Next: [Equation of state](02-equation-of-state.md)
