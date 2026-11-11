# 8. Rotation

## Implemented

None. All models are non-rotating.

## Why it matters

The ratio of centrifugal to gravitational acceleration at the equator is

```math
\epsilon_{\mathrm{rot}} = \frac{\Omega^2 R^3}{G M}
```

* Sun (period 25.4 d): $`2\times10^{-5}`$, negligible.
* Stars above about 1.3 M☉ are not braked by a magnetised wind. They typically rotate at
  150–300 km/s, which is 30–60% of the break-up speed (about 500 km/s for the 5 M☉ model).

For those stars rotation supports part of the weight, flattens the star, makes the poles hotter than
the equator, and mixes fresh hydrogen into the core (main sequence 20–30% longer in fast rotators,
literature value).

## Planned

### Stage 1: shape and gravity darkening (E3, plan 5.1)

Display only; the structure stays non-rotating.

*Roche surface.* With $`\omega = \Omega/\Omega_{\mathrm{crit}}`$,
$`x(\theta) = R(\theta)/R_{\mathrm{pole}}`$ and $`\theta`$ measured from the rotation axis:

```math
\Omega_{\mathrm{crit}}^{2} = \frac{8\,G M}{27\,R_{\mathrm{pole}}^{3}}, \qquad \frac{1}{x} + \frac{4}{27}\,\omega^{2}\,x^{2}\sin^{2}\theta = 1
```

```math
x(\theta) = \frac{3}{\omega\sin\theta}\,\cos\!\left[\frac{\pi + \arccos(\omega\sin\theta)}{3}\right]
```

At break-up ($`\omega = 1`$) the equatorial radius is 1.5 times the polar radius.

*Effective gravity*

```math
g_r = -\frac{G M}{R^{2}} + \Omega^{2} R \sin^{2}\theta, \qquad g_\theta = \Omega^{2} R \sin\theta\cos\theta, \qquad g_{\mathrm{eff}} = \sqrt{g_r^{2} + g_\theta^{2}}
```

*Gravity darkening* (von Zeipel 1924) **[check source]**: exponent 0.25 for radiative envelopes,
about 0.08 for convective ones. The constant is fixed by the luminosity.

```math
T_{\mathrm{eff}}(\theta) \propto g_{\mathrm{eff}}^{1/4}, \qquad \oint \sigma\, T_{\mathrm{eff}}^{4}(\theta)\, dA = L
```

### Stage 2: structure (E3, plan 5.2)

**A spherical approximation**, not the Kippenhahn & Thomas (1970) / Endal & Sofia (1976) method.
Here the radial centrifugal force is averaged over a sphere
($`\langle\sin^2\theta\rangle = 2/3`$) and the layers stay spherical. The full method works on
distorted equipotential surfaces and corrects two equations with geometric factors: hydrostatic
equilibrium ($`f_P`$) and radiative transport ($`f_T`$, entering as
$`\nabla_{\mathrm{rad}} \to \nabla_{\mathrm{rad}}\,f_T/f_P`$). Both corrections start at order
$`\Omega^2`$.

*Hydrostatic equilibrium*

```math
\frac{dP}{dm} = -\frac{G m}{4\pi r^4}\, f_P, \qquad f_P = 1 - \frac{2\,\Omega^2 r^3}{3\,G m} = 1 - \frac{3\,j^2}{2\,G m\, r}
```

*Specific angular momentum of a shell.* The second form of $`f_P`$ depends only on the local unknown
$`r`$ once $`j`$ is given, so it fits the existing solver.

```math
j = \tfrac{2}{3}\,\Omega\,r^{2}
```

*Evolution of the angular momentum*

* Radiative layers: $`j`$ of each shell is conserved.
* Convective zones rotate as a solid body:

```math
\Omega_{\mathrm{zone}} = \frac{\int_{\mathrm{zone}} j\,dm}{\int_{\mathrm{zone}} \tfrac{2}{3}\,r^{2}\,dm}
```

* ZAMS model: solid rotation at the chosen $`\omega`$.
* Total $`J = \int j\,dm`$ is conserved along the run only without mass loss. With a wind
  (phase 4), the lost mass carries its angular momentum away. For an isotropic wind that co-rotates
  with a spherical surface ($`\Omega_s`$ the surface rotation rate, $`\dot M < 0`$):

```math
\frac{dJ}{dt} = \frac{2}{3}\,\Omega_s\,R^{2}\,\dot M
```

  This follows automatically when the removed shells take their own $`j`$ with them. It excludes
  magnetic lever arms (magnetic braking, far stronger in cool stars) and anisotropic winds.

*Energy.* The total-energy equation of [Evolution in time](07-evolution.md) gains the rotational
kinetic energy:

```math
E_{\mathrm{rot}} = \int_0^M \tfrac{1}{3}\,\Omega^{2} r^{2}\,dm = \int_0^M \tfrac{1}{2}\,j\,\Omega\,dm, \qquad \Delta\!\left(E_{\mathrm{int}} + E_{\mathrm{grav}}\right) \to \Delta\!\left(E_{\mathrm{int}} + E_{\mathrm{grav}} + E_{\mathrm{rot}}\right)
```

*Virial check with rotation*

```math
\int_0^M \frac{3P}{\rho}\,dm + \int_0^M \tfrac{2}{3}\,\Omega^{2} r^{2}\,dm = \int_0^M \frac{G m}{r}\,dm
```

*Transport closure.* Radiative diffusion on spherical layers is unchanged by rotation:

```math
\frac{dT}{dm} = -\frac{3\,\kappa\,L}{64\pi^{2}\,a\,c\,r^{4}\,T^{3}}
```

Since $`\nabla = d\ln T/d\ln P`$ and $`dP/dm`$ now carries $`f_P`$, the radiative gradient must
carry the inverse factor, or the temperature equation would no longer be this diffusion law:

```math
\nabla_{\mathrm{rad,rot}} = \frac{\nabla_{\mathrm{rad}}}{f_P} = \frac{3\,\kappa\,L\,P}{16\pi\,a\,c\,G\,m\,T^{4}\,f_P}
```

This is the spherical limit ($`f_T = 1`$) of the full method's $`f_T/f_P`$. The same effective
gravity enters the convection formulas: $`g \to g\,f_P`$ in $`H_P`$, $`U`$ and
$`v_{\mathrm{conv}}`$; $`\nabla_{\mathrm{ad}}`$ is unchanged. The Schwarzschild criterion compares
$`\nabla_{\mathrm{rad,rot}}`$ with $`\nabla_{\mathrm{ad}}`$.

*Range of validity: not established.* No error bound is claimed. The approximation is expected to
be adequate while $`1 - f_P`$ is a few percent at the surface, and to fail well before break-up,
where the surface is far from spherical. Before the feature is presented as more than
illustrative, compare $`L`$, $`R`$ and $`T_c`$ against published rotating models at two or three
rotation rates and record the differences here.

Not included at this stage: rotational mixing, magnetic braking, the geometric factors
$`f_P`$ and $`f_T`$ of the full method.

### Stage 3: rotational mixing (plan phase 8)

Rotation mixes composition by slow circulation and shear. It enters as additional terms in the
diffusion coefficient $`D`$ of the composition-diffusion step
([Convection and mixing](05-convection-and-mixing.md)). The coefficients are to be chosen when that
step exists.

### With non-radial modes

Rotation splits a mode of radial order $`n`$ and degree $`l`$ into $`2l+1`$ components of azimuthal
order $`m = -l \dots l`$. To first order in $`\Omega`$, for solid rotation, in the observer's frame:

```math
\delta\nu_{nlm} = \nu_{nlm} - \nu_{nl0} = m\,\left(1 - C_{nl}\right)\frac{\Omega}{2\pi}
```

Sign convention: perturbations vary as $`\exp[i(m\phi - \omega t)]`$, so $`m > 0`$ is a prograde
mode and has the higher frequency. $`C_{nl}`$ is the Ledoux constant, an integral over the
eigenfunctions (radial and horizontal displacements $`\xi_r`$, $`\xi_h`$):

```math
C_{nl} = \frac{\int_0^R \left(2\,\xi_r\,\xi_h + \xi_h^{2}\right)\rho\,r^{2}\,dr}{\int_0^R \left[\xi_r^{2} + l(l+1)\,\xi_h^{2}\right]\rho\,r^{2}\,dr}
```

(Ledoux 1951; the integral agrees with Reese 2015, eq. 3, **[compared]**). Limits: $`C_{nl} \to 0`$ for high-order p
modes; $`C_{nl} \to 1/[l(l+1)]`$ for high-order g modes, so $`l = 1`$ g modes are split by only
half the rotation frequency. Valid only for slow rotation ($`\Omega \ll \omega`$ and
$`\Omega \ll \Omega_{\mathrm{crit}}`$); for differential rotation $`\Omega`$ moves inside a
kernel-weighted integral.

---
Previous: [Evolution in time](07-evolution.md) · [Contents](README.md) · Next: [Pulsations](09-pulsations.md)
