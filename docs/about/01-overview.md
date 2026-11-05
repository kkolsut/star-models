# 1. Overview

Star Models computes the internal structure of a star on the **zero-age main sequence (ZAMS)**. Such a
star:

* is chemically homogeneous;
* has just started burning hydrogen in its core;
* is in hydrostatic and thermal equilibrium.

You choose the mass and the composition. The app then solves the four differential equations of
stellar structure from the centre to the photosphere, with energy carried outward by radiation or by
convection. It shows how pressure, temperature, density and many other quantities vary through the
star.

Models are spherically symmetric, non-rotating and non-magnetic, without mass loss. The defaults are
those of the Sun:

* 1 M☉
* $`X = 0.70`$
* $`Z = 0.02`$
* $`\alpha_{\mathrm{MLT}} = 1.8`$

The default model is therefore the young Sun as it was about 4.6 billion years ago.

Press **Evolve** to follow the star in time:

* hydrogen burns to helium;
* the core contracts and the envelope expands;
* the star moves along an evolutionary track in the HR diagram.

Stars that cross the instability strip pulsate.

## Physics in use

| Component | Implementation |
|---|---|
| Equation of state | ideal gas + radiation (fully ionized) |
| Opacity | analytic (Kramers + H⁻ + e⁻ scattering) |
| Nuclear reactions | analytic pp + CNO + 3α (weak screening) |

The [Microphysics](07-microphysics.md) section gives every formula.

## What happens when you change a parameter

1. The parameters are checked: $`X + Z \le 1`$, and every value must be within its range.
2. A background thread solves the structure equations with a Newton–Raphson relaxation method (the
   Henyey method). It starts from the model currently on screen and moves step by step to the new
   parameters.
3. The converged solution is evaluated at every mesh point to give the profiles. The app locates the
   convective zones and computes global quantities and consistency checks.
4. Every view is redrawn from that model. The ZAMS line in the HR diagram is recomputed if the
   composition changed.

With **Evolve**, a background thread advances the star in time step by step (see
[Evolution in time](12-evolution.md)). It also computes the fundamental radial pulsation mode of every
model (see [Pulsations & oscillations](13-pulsations.md)).

> The following sections describe each of these steps in detail, with all formulas and numerical
> settings. In the app, right-click any equation to copy its LaTeX source.

---
[Contents](README.md) · Next: [Quick start & usage](02-usage.md) →
