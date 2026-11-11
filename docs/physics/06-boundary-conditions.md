# 6. Boundary conditions

Four first-order equations need four boundary conditions: two at the centre and two at the surface.

## Implemented

Source: `crates/stellar-core/src/structure/equations.rs` (`inner_bc`, `outer_bc`),
`crates/stellar-core/src/physics/atmosphere.rs`.

### Centre

The equations are singular at $`m = 0`$. The first mesh point sits at $`m_0 = 10^{-6} M`$ and uses
the leading terms of the expansion about the centre, with the local density and energy generation:

```math
r_0 = \left(\frac{3\, m_0}{4\pi\, \rho_0}\right)^{1/3}, \qquad L_0 = \left(\varepsilon_{\mathrm{nuc},0} + \varepsilon_{\mathrm{grav},0}\right) m_0
```

### Surface: grey Eddington atmosphere

```math
T^4(\tau) = \frac{3}{4}\, T_{\mathrm{eff}}^4 \left(\tau + \frac{2}{3}\right)
```

The outermost mesh point is the photosphere, at optical depth $`\tau = 2/3`$, where
$`T = T_{\mathrm{eff}}`$. The two conditions there:

```math
L = 4\pi R^2 \sigma T^4, \qquad P = \frac{2}{3}\, \frac{G M}{R^2\, \kappa} + \frac{a\, T^4}{6}
```

The pressure condition is hydrostatic equilibrium in optical depth, $`dP/d\tau = g/\kappa`$,
integrated from $`\tau = 0`$ with **constant** $`\kappa`$. At $`\tau = 0`$ only the radiation
pressure remains: $`P_{\mathrm{rad}}(\tau) = \tfrac{a}{4} T_{\mathrm{eff}}^4 (\tau + \tfrac23)`$,
which is $`a T_{\mathrm{eff}}^4/6`$ at the top. $`\kappa`$ is evaluated at the photospheric
$`(P, T)`$, so the condition is implicit and is solved together with the structure.

### Limits

* One value of $`\kappa`$ stands for the whole atmosphere. The H⁻ opacity varies as $`T^9`$, so this
  is a poor average in cool stars.
* The Eddington relation is grey. Real atmospheres of cool stars are a few hundred kelvin away from
  it at small $`\tau`$.

## Planned

### Integrated atmosphere (D13, plan phase 3)

Integrate from a small $`\tau_0`$ to the fitting point $`\tau_{\mathrm{fit}}`$ with the local
opacity:

```math
\frac{dP}{d\tau} = \frac{g}{\kappa\!\left(\rho(P, T),\, T\right)}, \qquad T^{4}(\tau) = \tfrac{3}{4}\,T_{\mathrm{eff}}^{4}\left[\tau + q(\tau)\right]
```

*Starting value (grey, LTE approximation for the radiation pressure).* $`P`$ is the total pressure.
At $`\tau = 0`$ the gas pressure vanishes; the radiation pressure there is taken as that of an
isotropic Planck field at the local temperature, and the first small step is taken with the opacity
at the top:

```math
P(0) = \frac{a\,T^{4}(0)}{3} = \frac{a}{4}\,T_{\mathrm{eff}}^{4}\,q(0), \qquad P(\tau_0) = P(0) + \frac{g\,\tau_0}{\kappa_0}
```

$`\kappa_0 = \kappa(P(\tau_0), T(\tau_0))`$ is found by iterating the last relation. The radiation
term is taken at $`\tau = 0`$, not at $`\tau_0`$: using $`a T^4(\tau_0)/3`$ would count the
radiation-pressure gradient twice (an excess of $`a T_{\mathrm{eff}}^4 \tau_0/4`$ in the Eddington
case). The start is exact for constant $`\kappa`$ and first order in $`\tau_0`$ otherwise, so
$`\tau_0`$ is chosen small enough ($`\sim10^{-4}`$) that the result does not depend on it.

With $`q = 2/3`$ and constant $`\kappa`$ the integration to $`\tau = 2/3`$ gives exactly the
implemented condition, which is the test for the new code. That test verifies the grey Eddington
limit only.

*What this approximation is.* $`P_{\mathrm{rad}} = aT^4/3`$ is not a general boundary condition for
a non-grey atmosphere. At small optical depth the radiation field is neither isotropic nor at the
local temperature, and a temperature Hopf function alone does not fix its pressure. With the
frequency-integrated moments $`J`$, $`H`$, $`K`$, the flux-mean and Rosseland opacities
$`\kappa_F`$, $`\kappa_R`$ and the Rosseland depth $`\tau_R`$:

```math
P_{\mathrm{rad}} = \frac{4\pi K}{c}, \qquad \frac{dP_{\mathrm{rad}}}{d\tau_R} = \frac{\kappa_F}{\kappa_R}\,\frac{F_{\mathrm{rad}}}{c}, \qquad \frac{dP_{\mathrm{mat}}}{d\tau_R} = \frac{g}{\kappa_R} - \frac{\kappa_F}{\kappa_R}\,\frac{F_{\mathrm{rad}}}{c}
```

A consistent treatment integrates the material pressure with a stated radiation boundary condition
and moment closure, and adds $`P_{\mathrm{rad}}`$ afterwards. These pages keep the simpler
prescription ($`\kappa_F = \kappa_R`$, $`K = J/3`$, $`J = B(T)`$) and state its range: it is
adequate where radiation pressure and radiative acceleration are small against gravity, which
covers the cool stars the non-grey relation is meant for. For hot stars, where they are not small,
the atmosphere stays the grey Eddington one.

### Convection in the atmosphere: Hopf-corrected transport

The $`T(\tau)`$ relation describes radiative equilibrium. In cool stars convection can reach above
the fitting point, and the transport equations used there must be consistent with the same
relation. Switching from the imposed $`T(\tau)`$ to the ordinary interior gradient of page 5 is
**not** consistent: it leaves a jump in $`\nabla`$ wherever $`q`$ still varies.

*Radiative gradient of the atmosphere.* Differentiating $`T^4 = \tfrac34 T_{\mathrm{eff}}^4[\tau + q(\tau)]`$
and using $`dP/d\tau = g/\kappa`$, with the total flux $`F_{\mathrm{tot}} = \sigma T_{\mathrm{eff}}^4`$:

```math
\nabla_{\mathrm{rad,atm}} = \frac{3\,\kappa\,F_{\mathrm{tot}}\,P}{16\,\sigma\,g\,T^{4}}\left[1 + \frac{dq}{d\tau}\right] = \nabla_{\mathrm{rad}}\left[1 + q'(\tau)\right]
```

The first factor is the interior $`\nabla_{\mathrm{rad}}`$ of page 5 (since
$`F_{\mathrm{tot}} = L/4\pi r^2`$ and $`\sigma = ac/4`$). The bracket is the Hopf correction; it
tends to 1 at depth, where $`q' \to 0`$ and the atmosphere joins the interior without a seam.

*Formulation.* The atmosphere is integrated as a pair of equations, not with an imposed
temperature:

```math
\frac{dP}{d\tau} = \frac{g}{\kappa}, \qquad \frac{d\ln T}{d\ln P} = \nabla, \qquad \nabla = \begin{cases} \nabla_{\mathrm{rad,atm}} & \nabla_{\mathrm{rad,atm}} \le \nabla_{\mathrm{ad}} \\ \text{modified mixing-length solution (below)} & \text{otherwise} \end{cases}
```

*The mixing-length cubic with the Hopf correction.* Replacing $`W`$ alone is not enough. In the
cubic of page 5, $`U`$ plays two roles: it measures the radiative losses of a convective element,
and (through the coefficient $`8U/9`$) it carries the conductivity of the background, which
converts the convective flux into a gradient difference. The Hopf correction changes the
background law only: the atmosphere transports the flux as if its conductivity were divided by

```math
A = 1 + q'(\tau) > 0, \qquad W_A = A\,\nabla_{\mathrm{rad}} - \nabla_{\mathrm{ad}} = \nabla_{\mathrm{rad,atm}} - \nabla_{\mathrm{ad}}
```

Flux balance, $`F_{\mathrm{rad}}/F_{\mathrm{tot}} = \nabla/(A\nabla_{\mathrm{rad}})`$, then gives
$`A\nabla_{\mathrm{rad}} - \nabla = A\,(9/8U)\,\eta^3`$ in place of page 5's
$`\nabla_{\mathrm{rad}} - \nabla = (9/8U)\,\eta^3`$, while the element relation
$`\nabla = \nabla_{\mathrm{ad}} + \eta^2 + 2U\eta`$ is unchanged:

```math
\eta^{3} + \frac{8U}{9A}\left(\eta^{2} + 2U\eta - W_A\right) = 0, \qquad \nabla = \nabla_{\mathrm{ad}} + \eta^{2} + 2U\eta
```

It reduces to page 5 for $`A = 1`$. This closure keeps page 5's **optically thick** prescription
for the element losses, which is an approximation in the atmosphere. With a separate element-loss
parameter $`U_e`$ (from a thin-element prescription, see below) the general form is

```math
\eta^{3} + \frac{8U}{9A}\left(\eta^{2} + 2U_e\,\eta - W_A\right) = 0, \qquad \nabla = \nabla_{\mathrm{ad}} + \eta^{2} + 2U_e\,\eta
```

$`U_e`$ is set by the element prescription; it is not obtained by rescaling every $`U`$. (Form
proposed by the Codex review of 2026-11-11 and re-derived here from the flux balance.)

In radiative layers the formulation reproduces the $`T(\tau)`$ relation exactly; in convective
layers $`\nabla`$ is continuous at the onset because the cubic gives
$`\nabla \to \nabla_{\mathrm{rad,atm}}`$ as $`W_A \to 0`$. The same corrected gradient is used in the interior equations down to the depth
where $`|q'|`$ is negligible, so there is one transport law on both sides of the fitting point.
The split of the flux then follows from the corrected gradient:

```math
\frac{F_{\mathrm{rad}}}{F_{\mathrm{tot}}} = \frac{\nabla}{\nabla_{\mathrm{rad,atm}}}, \qquad \frac{F_{\mathrm{conv}}}{F_{\mathrm{tot}}} = 1 - \frac{\nabla}{\nabla_{\mathrm{rad,atm}}}
```

Source: Trampedach et al. (2014), §2.5, eqs. (35)–(38) **[compared]** (gradient and flux fractions,
within that construction's assumptions; the comparison does not cover the cubic above).

*Fitting point.* Both values come from the integration:

```math
T_N = T_{\mathrm{atm}}\!\left(\tau_{\mathrm{fit}};\, g,\, T_{\mathrm{eff}}\right), \qquad P_N = P_{\mathrm{atm}}\!\left(\tau_{\mathrm{fit}};\, g,\, T_{\mathrm{eff}}\right)
```

*Simpler alternative, explicitly approximate.* Keep the uncorrected interior gradient and require
the fitting point to lie so deep that $`|q'(\tau_{\mathrm{fit}})|`$ is below a stated tolerance
(for the relation below, $`|q'| < 10^{-3}`$ needs $`\tau_{\mathrm{fit}} \gtrsim 3`$), with the
atmosphere above integrated as formulated here. This is acceptable as a first version if labelled.

*Which kind of $`q(\tau)`$.* What matters is not where a Hopf function comes from but what it
describes. A *radiative* Hopf function gives the temperature the layer would have if radiation
carried the whole flux; the formulation above is built for it. Such a function can be extracted
from a convective simulation by reducing the stratification to its radiative part (the Trampedach
construction does this), and using it causes no double counting. A fit to the *full* temperature
stratification of an atmosphere with convection (an observed solar $`T(\tau)`$, for example) is a
different object: part of the convective effect is already inside it, and applying mixing-length
theory on top counts that part twice to some degree. Which of the two a chosen relation is must be
stated with it.

### Thin convective elements: the loss parameter $`U_e`$ (decided)

**Decided 2026-11-11** (recommendation of the Codex review, accepted).

The formulas of page 5 treat the radiative losses of a convective element in the diffusion
approximation, valid when the element is optically thick. The controlling quantity is the optical
thickness of the element, not the depth in the atmosphere alone:

```math
\tau_e = \kappa\,\rho\,\ell
```

*Source.* Henyey, Vardya & Bodenheimer (1965, eqs. 32, 34–36, 39–40) **[compared]** define the
convective efficiency $`\gamma`$ and give its transparent and opaque limits and an interpolation,
with $`y`$ a number describing the temperature profile inside the element:

```math
\gamma = \frac{\nabla - \nabla_e}{\nabla_e - \nabla_{\mathrm{ad}}} = \frac{c_P\,\rho\,v}{8\,\sigma\,T^{3}}\;\frac{1 + y\,\tau_e^{2}}{\tau_e}
```

*Translation into the notation of page 5* (derived in the review, re-derived here). Page 5 has
$`\nabla - \nabla_e = \eta^2`$, $`\nabla_e - \nabla_{\mathrm{ad}} = 2U_e\eta`$ and
$`v = \eta\,\ell\sqrt{g\delta/8H_P}`$, so $`\gamma = \eta/(2U_e)`$ and

```math
U_e = U\,\frac{\tau_e^{2}}{3\,(1 + y\,\tau_e^{2})}
```

$`U_e \to U`$ for opaque elements requires $`y = 1/3`$. With that choice:

```math
U_e = U\,\frac{\tau_e^{2}}{3 + \tau_e^{2}} = \frac{3\,a\,c\,T^{3}\,\kappa}{c_P\left(3 + \tau_e^{2}\right)}\sqrt{\frac{8 H_P}{g\,\delta}}
```

The second form is the one to code: it does not multiply a diverging $`U`$ by a vanishing factor.
$`U_e`$ replaces $`U`$ only in the two element-loss places of the general cubic above; the
coefficient $`8U/(9A)`$ keeps $`U`$, and $`U_e`$ is not divided by $`A`$.

| Limit | Result |
|---|---|
| $`\tau_e \gg 1`$ | $`U_e \to U`$: page 5 (with $`A = 1`$) |
| $`\tau_e \ll 1`$ | $`U_e \to U\tau_e^2/3 = (a c T^3 \kappa/c_P)\sqrt{8H_P/(g\delta)}`$: losses vanish with the opacity, and the $`1/\kappa`$ divergence of the diffusion formula is gone |

No new free parameter: $`\alpha_{\mathrm{MLT}}`$ stays the calibration parameter. Henyey's
alternative $`y = 3/(4\pi^2)`$ would give $`U_e/U \to 4\pi^2/9`$ at depth, changing the
normalisation of deep convection; it is not used.

What stays approximate: local mixing-length theory with an interpolated radiative exchange; and one
Rosseland opacity used both for transport and for the emission of a transparent element, which
physically involves absorption (Planck-type) means.

### Non-grey $`T(\tau)`$ relation (F3) **[check source]**

Empirical solar relation of Krishna Swamy (1966):

```math
q(\tau) = 1.39 - 0.815\,e^{-2.54\,\tau} - 0.025\,e^{-30\,\tau}
```

The coefficients agree with a contemporaneous paper that reproduces the relation (its eq. 4); the
original has not been read. Its optical-depth scale is also unverified: if it is a monochromatic
depth ($`\tau_{5000}`$) and not the Rosseland depth used here, it must be converted with
$`d\tau_{5000}/d\tau_R = \kappa_{5000}/\kappa_R`$ before use; the two scales are not
interchangeable **[check source]**. It is a fit to the observed Sun, so it is of the second kind
above (convection is partly inside it), and its supported range is about $`0.02 \le \tau \le 10`$
(Trampedach et al. 2014, §5). Using it down to $`\tau = 0`$ for the starting pressure is an
extrapolation and is labelled as such; so is using it for stars unlike the Sun.

### Surface conditions with a general $`q(\tau)`$

The outermost mesh point is then no longer at $`T = T_{\mathrm{eff}}`$. The unknown
$`T_{\mathrm{eff}}`$ is defined by the luminosity, and the two conditions on the last point $`N`$
become (written for an atmosphere that is radiative down to the fitting point; in general $`T_N`$
comes from the integration of the previous section, which gives the same value in that case):

```math
L = 4\pi R^{2}\sigma T_{\mathrm{eff}}^{4}, \qquad T_N^{4} = \tfrac{3}{4}\,T_{\mathrm{eff}}^{4}\left[\tau_{\mathrm{fit}} + q(\tau_{\mathrm{fit}})\right], \qquad P_N = P_{\mathrm{atm}}\!\left(\tau_{\mathrm{fit}};\, g,\, T_{\mathrm{eff}}\right)
```

*Which radius.* $`R`$ in these relations is the photospheric radius, where $`T = T_{\mathrm{eff}}`$.
In a radiative atmosphere that is where $`\tau + q(\tau) = 4/3`$; once convection changes the
integrated temperature, the identity is not available and the depth $`\tau_{\mathrm{ph}}`$ is
located from the actual $`T(\tau)`$ of the integration. The last mesh point sits at $`r_N = R - \Delta r`$, with

```math
\Delta r = \int_{\tau_{\mathrm{ph}}}^{\tau_{\mathrm{fit}}} \frac{d\tau}{\kappa\,\rho}
```

from the atmosphere integration. For $`\tau_{\mathrm{fit}}`$ of order 1 in a dwarf,
$`\Delta r / R \sim 10^{-4}`$ and $`r_N = R`$ is acceptable (plane-parallel atmosphere, $`g`$
constant). For a deep fitting point or a giant it is not, and $`\Delta r`$ and the mass of the
atmosphere must be carried explicitly.

---
Previous: [Convection and mixing](05-convection-and-mixing.md) · [Contents](README.md) · Next: [Evolution in time](07-evolution.md)
