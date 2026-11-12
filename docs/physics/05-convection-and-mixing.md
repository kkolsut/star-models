# 5. Convection and mixing

This page covers three things: where energy is carried by convection, what the temperature gradient
is there, and how convection mixes the composition.

## Implemented

Source: `crates/stellar-core/src/structure/local.rs` (radiative gradient),
`crates/stellar-core/src/physics/convection.rs` (mixing-length theory),
`crates/stellar-core/src/structure/evolve.rs` (mixing).

### Radiative gradient

The gradient needed to carry the whole luminosity by radiation:

```math
\nabla_{\mathrm{rad}} = \frac{3\,\kappa\, L\, P}{16\pi\, a\, c\, G\, m\, T^4}
```

### Schwarzschild criterion

```math
\nabla_{\mathrm{rad}} > \nabla_{\mathrm{ad}} \;\Rightarrow\; \text{convective}, \qquad \nabla_{\mathrm{rad}} \le \nabla_{\mathrm{ad}} \;\Rightarrow\; \text{radiative},\; \nabla = \nabla_{\mathrm{rad}}
```

It is applied locally, in ZAMS and evolved models alike. Composition gradients are not taken into
account (no Ledoux criterion).

### Mixing-length theory

Formulation of Kippenhahn, Weigert & Weiss (2012, §7), after Böhm-Vitense (1958).

*Gravity, pressure scale height, mixing length*

```math
g = \frac{G m}{r^2}, \qquad H_P = \frac{P}{\rho\, g}, \qquad \ell = \alpha_{\mathrm{MLT}}\, H_P
```

*Efficiency parameter and driving*

```math
U = \frac{3\, a\, c\, T^3}{c_P\, \rho^2\, \kappa\, \ell^2}\sqrt{\frac{8 H_P}{g\, \delta}}, \qquad W = \nabla_{\mathrm{rad}} - \nabla_{\mathrm{ad}}
```

*The cubic*, and the form solved in the code with $`\eta = \xi - U \ge 0`$

```math
(\xi - U)^3 + \frac{8U}{9}\left(\xi^2 - U^2 - W\right) = 0 \quad\Longleftrightarrow\quad \eta^3 + \frac{8U}{9}\left(\eta^2 + 2U\eta - W\right) = 0
```

Newton's method starts from the upper bound $`\eta_{\max} = \sqrt{U^2 + W} - U`$ and converges
monotonically.

*Gradient, velocity and flux fraction*

```math
\nabla = \nabla_{\mathrm{ad}} + \eta\,(\eta + 2U), \qquad v_{\mathrm{conv}} = \eta\, \ell \sqrt{\frac{g\, \delta}{8 H_P}}, \qquad \frac{F_{\mathrm{conv}}}{F} = \frac{\nabla_{\mathrm{rad}} - \nabla}{\nabla_{\mathrm{rad}}}
```

Limits: $`U \ll 1`$ gives $`\nabla \to \nabla_{\mathrm{ad}}`$ (efficient, deep interior);
$`U \gg 1`$ gives $`\nabla \to \nabla_{\mathrm{rad}}`$ (inefficient, just below the photosphere).

### Mixing

After each burning step, every convective zone of the previous model is mixed completely:

```math
X \;\leftarrow\; \bar X_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} X\, dm}{\int_{\mathrm{zone}} dm} \qquad \text{(same for } Y, Z\text{)}
```

### Step overshoot of convective cores (D16, plan 2.1)

Input $`\alpha_{\mathrm{ov}}`$ (`StarParams::alpha_ov`), default 0 (off), range 0–0.5. Before the
mixing step, a convective zone that starts at the innermost point is extended from its Schwarzschild
boundary $`r_{\mathrm{cc}}`$ (interpolated linearly in $`\nabla_{\mathrm{rad}} - \nabla_{\mathrm{ad}}`$) to

```math
r_{\mathrm{ov}} = r_{\mathrm{cc}} + \alpha_{\mathrm{ov}}\,\min\!\left(H_P(r_{\mathrm{cc}}),\ r_{\mathrm{cc}}\right), \qquad \alpha_{\mathrm{ov}} \approx 0.2
```

* Every mesh point with $`r \le r_{\mathrm{ov}}`$ is mixed with the core, and the burning rate is
  averaged over the same region for the time-step limit. A convective run reached by the extension is
  merged into it.
* $`H_P(r_{\mathrm{cc}})`$ is interpolated linearly in $`r`$ between the two points around the
  boundary. At the boundary $`H_P`$ is finite, so the bounded scale height of plan 1.4 is not needed
  here; $`\min(H_P, r_{\mathrm{cc}})`$ limits small cores.
* Mixing only: the temperature gradient in the extension stays radiative, and the structure
  equations do not see $`\alpha_{\mathrm{ov}}`$. A ZAMS model is identical for every value.
* Output: the mixed-core mass and radius (`mixed_core_q`, `mixed_core_r`) next to the convective core;
  equal to it for $`\alpha_{\mathrm{ov}} = 0`$.

Envelope overshoot and overshoot of convective shells are not included.

There is no other mixing outside the Schwarzschild boundaries.

### How good this is **[measured]**

* **Gradient.** Deep convection is adiabatic: $`\nabla - \nabla_{\mathrm{ad}} \le 1.4\times10^{-7}`$
  in the 5 M☉ core. In the Sun it peaks at 0.115 at $`r = 0.9997\,R`$.
* **$`\alpha_{\mathrm{MLT}}`$ has almost no effect today.** For 1 M☉, $`\alpha = 1.0, 1.8, 2.5`$
  gives $`R = 1.1207, 1.1133, 1.1106`$ R☉. The cause is the shallow convective envelope
  ($`2\times10^{-6}`$ of the mass), which comes from the opacity and the EOS, not from the
  convection theory.
* **Instant mixing is justified.** Gas crosses the 5 M☉ core in about a month
  ($`v_{\mathrm{conv}} \approx 10^4`$ cm/s); the main sequence lasts 76 Myr.
* **The mixing length is unbounded at the centre**, because $`H_P \to \infty`$ as $`g \to 0`$:
  $`\alpha H_P / r = 4900`$ at the innermost point of the 5 M☉ core. Only the plotted velocity is
  affected.
* **Overshoot lengthens the main sequence.** Main-sequence lifetime for
  $`\alpha_{\mathrm{ov}} = 0, 0.1, 0.2`$: 285.9, 323.6, 360.1 Myr for 3 M☉ and 76.1, 85.2, 93.7 Myr for
  5 M☉ (+26% and +23% at 0.2). The largest mixed core of the 3 M☉ run grows from 0.176 to 0.243 of the
  mass. 1.1 M☉, which keeps a radiative core, is bit-identical with and without overshoot to 5 Gyr.
  Fuel accounting closes to 0.15% in all six runs (`tests/overshoot.rs`). The default stays 0 until
  the joint assessment (plan 7.7).
* **The step is not a prediction.** Without overshoot (the default) the convection-related error of the
  evolution results is still dominated by its absence.

## Planned

### Bounded scale height (C9, plan 1.4)

Finite at the centre, and needs no extra input:

```math
\tilde H_P = \min\!\left(\frac{P}{\rho\, g},\ \sqrt{\frac{P}{G\,\rho^{2}}}\right), \qquad \ell = \alpha_{\mathrm{MLT}}\,\tilde H_P
```

Source: Paxton et al. (2011, ApJS 192, 3), §5.1, unnumbered (arXiv:1009.1622v1, p. 28): the
alternate definition $`\lambda_P' = (P/G\rho^2)^{1/2}`$ is used when $`\lambda_P' < \lambda_P = P/g\rho`$,
which is the minimum above. Paxton et al. attribute it to Eggleton (1971), but that paper
(MNRAS 151, 351–364, read in full) does not contain it: it uses the ordinary $`H_P`$ in its
mixing-length formulae (Appendix, eq. 57). Eggleton (1972, MNRAS 156, 361), p. 364, only notes
that $`H_P \to \infty`$ at the centre and "perhaps should be replaced by something which tends to
a finite value"; it gives no formula. Cite Paxton et al. 2011 for the formula.

### Ledoux criterion and semiconvection (D17, plan phase 8)

```math
\nabla_{\mathrm{L}} = \nabla_{\mathrm{ad}} + \frac{\varphi}{\delta}\,\nabla_\mu, \qquad \nabla_\mu = \frac{d\ln\mu}{d\ln P}, \qquad \varphi = \left(\frac{\partial\ln\rho}{\partial\ln\mu}\right)_{P,T}
```

A layer with $`\nabla_{\mathrm{ad}} < \nabla_{\mathrm{rad}} < \nabla_{\mathrm{L}}`$ is semiconvective:
unstable by Schwarzschild, stabilised by the composition gradient, and mixed slowly.
$`\varphi = 1`$ for the implemented EOS.

*General composition term.* $`\varphi\,\nabla_\mu`$ assumes the composition enters the EOS only
through one mean molecular weight. That holds for the implemented EOS. It does not for the planned
ones, where $`\rho`$ depends separately on the ion and electron abundances, on ionisation and later
on non-ideal terms. The general form, with the sum over independent species (one mass fraction is
eliminated by $`\sum_i X_i = 1`$):

```math
B = \sum_i \left(\frac{\partial\ln\rho}{\partial X_i}\right)_{P,T,X_{j\ne i}} \frac{dX_i}{d\ln P}, \qquad \nabla_{\mathrm{L}} = \nabla_{\mathrm{ad}} + \frac{B}{\delta}
```

Reading of the derivative: when an independent $`X_i`$ is varied, the eliminated species varies by
the opposite amount so that the sum stays 1, and all other independent $`X_j`$ are held fixed. It
is not a derivative with every other mass fraction fixed.

Convention: this $`B`$ equals $`\delta`$ times the $`B`$ of Paxton et al. (2013, §3.3, eqs. 5–8),
where $`N^2 = (g^2\rho/P)\,\delta\,(\nabla_{\mathrm{ad}} - \nabla + B_{\mathrm{MESA}})`$. The two
forms were compared and are equivalent.

$`B`$ reduces to $`\varphi\nabla_\mu`$ for an ideal mixture. The partial derivatives come from the
EOS (by differences at fixed $`P`$, $`T`$). From phase 6 on, $`B`$ is the production expression, in
the Ledoux criterion and in $`N^2`$ ([Pulsations](09-pulsations.md)).

### Mixing as diffusion (plan phase 8)

Replaces instant mixing; solved implicitly on the mass mesh for each species (Langer, El Eid &
Fricke 1985, A&A 145, 179, eq. 1, p. 180, without its nuclear term):

```math
\frac{\partial X_i}{\partial t} = \frac{\partial}{\partial m}\left[\left(4\pi r^{2}\rho\right)^{2} D\,\frac{\partial X_i}{\partial m}\right]
```

*Diffusion coefficients*: convection, and exponential overshoot beyond a reference point $`r_0`$
(Herwig 2000, $`f \approx 0.016`$, p. 954)

```math
D_{\mathrm{conv}} = \tfrac{1}{3}\,v_{\mathrm{conv}}\,\ell, \qquad D_{\mathrm{ov}}(r) = D_0\,\exp\!\left(-\frac{2\,|r - r_0|}{f\,H_{P,0}}\right)
```

Sources: $`D_{\mathrm{conv}}`$ is Langer et al. (1985), p. 180, text before eq. 2 ("$`D_c = 1/3\,v_c l`$",
after Weaver et al. 1978). $`D_{\mathrm{ov}}`$ is Herwig (2000), §3, eq. 2, p. 954
($`D_{\mathrm{OV}} = D_0 \exp(-2z/H_v)`$, $`H_v = f H_P`$, $`z`$ the geometric distance from the
edge of the convective zone), in the form with $`r_0`$ and $`H_{P,0}`$ of Paxton et al. (2011),
§5.2, eq. 2 ($`\lambda_{P,0}`$ "the pressure scale height at that location"). Pedersen et al.
(2018, A&A 614, A128, §2.2, eq. 2; compared in arXiv:1802.02051v1) write the same exponential with
$`H_{P}`$ at the Schwarzschild boundary instead of at $`r_0`$; the two differ at order $`f_0`$.

*Where $`D_0`$ is taken.* Not at the boundary itself: mixing-length theory gives
$`v_{\mathrm{conv}} \to 0`$ there, so $`D_0`$ would vanish and with it the overshoot. $`D_0`$ is
taken at a reference point $`r_0`$ inside the convective zone, a distance $`f_0 H_P`$ from the
Schwarzschild boundary $`r_{\mathrm{b}}`$ (on the convective side):

```math
|r_0 - r_{\mathrm{b}}| = f_0\,H_P(r_{\mathrm{b}}), \qquad D_0 = D_{\mathrm{conv}}(r_0), \qquad H_{P,0} = H_P(r_0)
```

The exponential of the formula above starts at $`r_0`$, with scale $`f\,H_{P,0}`$. Between $`r_0`$
and the far side it replaces $`D_{\mathrm{conv}}`$. $`f`$ and $`f_0`$ are two separate calibrated
parameters; $`f_0`$ does not modify the scale. *Value of $`f_0`$.* No source fixes it, and the
earlier "$`f_0 \approx 0.005`$ is a common choice" is withdrawn: Herwig (2000) has no $`f_0`$; the
MESA controls documentation (`overshoot_f0`, read 2026-11-11) gives default 0, requires the user
to set a positive value, and illustrates it with $`f_0 = 0.05`$, $`f = 0.25`$; Pedersen et al.
(2018, §3.1, arXiv v1 p. 4) use $`f_0 = 0.001`$ ("≈ 2 steps in resolution" into the core) with
$`f = 0.010`$–$`0.020`$ (their Table 2). Here $`f_0`$ is a user parameter with no sourced default;
it should be small compared with $`f`$ and a few mesh points wide. An exponential has no fixed
reach: how far the mixing is effective depends on the time available, and in the code on a
cut-off, below which $`D`$ is set to zero (a stated parameter; a relative one such as
$`D < 10^{-10} D_0`$ is a modelling choice of these pages, while Herwig 2000, p. 954, uses an
absolute $`10^{-2}`$ cm² s⁻¹). Two sources, cited separately: the exponential, its factor 2, a
normalisation near the boundary ($`D_0`$ "the diffusion coefficient in the convectively unstable
region near the convective boundary") and the absolute cut-off are Herwig's (2000, §3, eq. 2,
p. 954; earlier quoted here as §2) **[compared]**; the explicit reference point a distance
$`f_0 H_P`$ inside the zone, with the exponential starting there, is the convention of the MESA
code (Paxton et al. 2011, §5.2, eq. 2; `overshoot_f0` control; Pedersen et al. 2018, §2.1–2.2,
$`r_0 = r_{\mathrm{cc}} - f_0 H_{P,\mathrm{cc}}`$) **[compared]**, not Herwig's own
boundary-distance definition. $`H_P`$ is the bounded scale height $`\tilde H_P`$ for small cores.

*The parameters are not universal.* $`f \approx 0.016`$ and $`\alpha_{\mathrm{ov}} \approx 0.2`$ are
typical calibrations for intermediate-mass main-sequence stars; observations suggest a dependence on
mass. The two are not interchangeable by a fixed factor (a rough guide is
$`\alpha_{\mathrm{ov}} \sim 10 f`$), and they give different composition profiles: a step produces a
sharp edge, the exponential a smooth one. Both stay user parameters.

Instant mixing is the limit $`D \to \infty`$ inside convective zones and $`D = 0`$ elsewhere; the
new step must reproduce it. Rotational mixing ([Rotation](08-rotation.md)) adds further terms to
$`D`$.

### Gravitational settling (F1, plan phase 8)

Helium and metals sink slowly relative to hydrogen. An advection term with the diffusion velocity
$`w_i`$ of each species, from Thoul, Bahcall & Loeb (1994):

```math
\frac{\partial X_i}{\partial t} = -\frac{\partial}{\partial m}\left(4\pi r^{2}\rho\, X_i\, w_i\right)
```

Source: Thoul, Bahcall & Loeb (1994), eqs. 1–2, p. 829:
$`\partial X_s/\partial t = -(\rho r^2)^{-1}\,\partial_r[r^2 X_s T^{5/2}\xi_s]`$ with
$`\xi_s = w_s\rho/T^{5/2}`$ (Lagrangian time derivative); with $`dm = 4\pi r^2\rho\,dr`$ this is
the form above. $`w_s`$ is measured relative to the mean (mass-weighted) velocity (eqs. 5–6,
p. 830), so $`\sum_s \rho_s w_s = 0`$ (eq. 19, p. 831). The velocities come from Burgers'
equations with residual heat flows (eqs. 11–13, 17–19, pp. 830–831), solved as the linear system
eqs. 21–38 (pp. 831–832); radiative forces are neglected (p. 828). Their analytic fits
$`\xi_s = A_p\,\partial_r\ln p + A_T\,\partial_r\ln T + A_H\,\partial_r\ln C_H`$ (eq. 3, p. 829) with
quadratic coefficients in $`X`$ for H, O and Fe (eqs. 56–58, p. 839) use constant Coulomb
logarithms at their solar-centre values and are stated to hold for the solar interior only
(p. 840); $`w_{\mathrm{He}}`$ then follows from eq. 54 (p. 836). Outside the Sun the full system,
not the fits, is needed.

Size in the Sun over 4.6 Gyr: surface $`Y`$ falls by about 0.03 and $`Z`$ by about 10%. Checked
2026-11-12 against the standard solar models of Vinyoles et al. (2017, ApJ 835, 202, Table 4;
arXiv:1611.09867v4): B16-GS98 has $`Y_{\mathrm{ini}} = 0.2718`$, $`Y_S = 0.2426`$ ($`\Delta Y = -0.029`$)
and $`Z_{\mathrm{ini}} = 0.0187`$, $`Z_S = 0.0170`$ (−9%); B16-AGSS09met has 0.2613 → 0.2317
($`-0.030`$) and 0.0149 → 0.0134 (−10%). These are a test target for plan 8.5, not inputs.

### Mixing-length cubic in the atmosphere (D13, plan phase 3)

With a non-grey $`T(\tau)`$ relation the cubic above changes near the surface: its coefficient and
its driving term take the Hopf factor $`A = 1 + q'(\tau)`$, and optically thin elements need a
separate loss parameter. The equations are in
[Boundary conditions](06-boundary-conditions.md); they reduce to the cubic of this page for
$`A = 1`$.

### Not planned

* A different convection theory (full-spectrum turbulence, 3D-calibrated $`\alpha`$, time-dependent
  convection). They change only the thin superadiabatic layer and bring nothing before the opacity
  and EOS upgrades.
* Turbulent pressure. It is about 10% of the pressure at the top of the solar convection zone and
  shifts oscillation frequencies, but hardly changes $`R`$ or $`L`$.

---
Previous: [Nuclear burning](04-nuclear-burning.md) · [Contents](README.md) · Next: [Boundary conditions](06-boundary-conditions.md)
