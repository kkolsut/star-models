# 3. Input parameters

These parameters define a model. The ranges are those of the sliders in the left panel; the app also
shows each current value.

| Parameter | Symbol | Default | Range |
|---|---|---|---|
| Mass | $`M`$ [M☉] | 1.0 | 0.1 – 100 (logarithmic) |
| Hydrogen mass fraction | $`X`$ | 0.70 | 0 – $`(1 - Z)`$ |
| Metallicity | $`Z`$ | 0.02 | 0.0001 – 0.1 (logarithmic) |
| Helium mass fraction | $`Y`$ | derived | $`Y = 1 - X - Z`$ |
| Mixing-length parameter | $`\alpha_{\mathrm{MLT}}`$ | 1.8 | 0.3 – 3.0 |
| CNO fraction of $`Z`$ (Advanced) | $`f_{\mathrm{CNO}}`$ | 0.7 | 0 – 1 |
| Core overshoot (Advanced) | $`\alpha_{\mathrm{ov}}`$ | 0 (off) | 0 – 0.5 |
| Mesh points (Advanced) | $`N`$ | 800 | 200 – 3000 (logarithmic) |

## Mass $`M`$

The total mass of the star; it determines almost everything else.

* Low-mass stars are cool, dim and small, with deep convective envelopes. At the lowest masses they are
  fully convective.
* Above about 1.2 M☉ a convective core appears and the CNO cycle takes over from the pp chain.
* Massive stars are hot and luminous, and partly supported by radiation pressure ($`\beta < 1`$).

Below about 0.4 M☉ electron degeneracy and non-ideal gas effects matter. They are not modelled, so
results there are qualitative only.

## Hydrogen mass fraction $`X`$

The fraction of the mass in hydrogen. More hydrogen lowers the mean molecular weight $`\mu`$, which
makes the star dimmer and cooler at a given mass. $`X`$ also enters:

* the pp-chain rate ($`\propto X^2`$);
* the CNO rate ($`\propto X`$);
* the electron-scattering opacity $`0.2\,(1 + X)`$.

## Metallicity $`Z`$

The fraction of the mass in all elements heavier than helium.

* $`Z`$ controls the bound–free and H⁻ opacities. A higher $`Z`$ makes the star more opaque: dimmer,
  cooler and larger, with a deeper convective envelope.
* $`Z`$ also sets the abundance of the CNO catalysts.

## Helium mass fraction $`Y`$

Not set directly: it is whatever mass is left after hydrogen and metals, $`Y = 1 - X - Z`$.

## Mixing-length parameter $`\alpha_{\mathrm{MLT}}`$

In mixing-length theory a convective blob travels a distance $`\ell = \alpha_{\mathrm{MLT}}\,H_P`$
before it dissolves into its surroundings. $`H_P`$ is the pressure scale height.

* $`\alpha_{\mathrm{MLT}}`$ only matters where convection is inefficient: in the superadiabatic layers
  just below the photosphere of stars with convective envelopes.
* There, a larger $`\alpha_{\mathrm{MLT}}`$ means more efficient convection, a smaller radius and a
  higher $`T_{\mathrm{eff}}`$.
* Deep convection is almost exactly adiabatic, whatever the value of $`\alpha_{\mathrm{MLT}}`$.

Real solar models calibrate $`\alpha_{\mathrm{MLT}} \approx 1.7\text{–}2.2`$ to match the Sun.

## CNO fraction of $`Z`$: $`f_{\mathrm{CNO}}`$ (Advanced)

The fraction of the metals that is carbon, nitrogen and oxygen. The CNO-cycle rate uses
$`X_{\mathrm{CNO}} = f_{\mathrm{CNO}}\,Z`$. It only matters for stars hot enough to burn hydrogen via
the CNO cycle.

## Core overshoot $`\alpha_{\mathrm{ov}}`$ (Advanced)

Step overshoot of a convective core, used by evolution runs only. Before each mixing step the
convective core is extended from its Schwarzschild boundary $`r_{\mathrm{cc}}`$ to
$`r_{\mathrm{ov}} = r_{\mathrm{cc}} + \alpha_{\mathrm{ov}}\min(H_P(r_{\mathrm{cc}}), r_{\mathrm{cc}})`$, and
the whole region is mixed (see [Evolution in time](12-evolution.md#composition)).

* Mixing only: the temperature gradient in the extension stays radiative, so a ZAMS model does not
  change.
* A larger $`\alpha_{\mathrm{ov}}`$ feeds more hydrogen into the core and lengthens the main sequence of
  stars with convective cores. Stars without one are not affected.
* Calibrations for intermediate-mass stars give $`\alpha_{\mathrm{ov}} \approx 0.2`$, with a dependence
  on mass; the value is not universal.

## Mesh points $`N`$ (Advanced)

The number of points of the mass grid (see [Numerical method](10-numerical-method.md)). More points
give a more accurate but slower solution. Doubling from 800 to 1600 changes $`L`$, $`R`$ and $`T_c`$ by
less than 0.5% for the Sun.

## Presets

| Preset | Values |
|---|---|
| **Sun (X 0.70, Z 0.02)** | 1 M☉, the classic textbook solar composition, $`\alpha_{\mathrm{MLT}} = 1.8`$, $`f_{\mathrm{CNO}} = 0.7`$ |
| **Sun AGSS09** | the same, with the protosolar composition of Asplund et al. (2009): $`X = 0.7154`$, $`Z = 0.0142`$ |

## Evolution settings

**Stop at** offers three choices:

* the end of the main sequence ($`X_c < 10^{-3}`$, the default);
* a given age in Gyr (default 4.57);
* as far as possible.

The time-step controls are fixed (see [Evolution in time](12-evolution.md)). The mesh-points setting
also applies to evolution: more points resolve the thin hydrogen-burning shell after the main sequence
better.

## Display settings (they do not change the model)

| Setting | Effect |
|---|---|
| **x axis** | $`r/R`$ (fractional radius), $`m/M`$ (fractional mass, the solver's own coordinate) or $`\log P`$ (pressure, centre on the left) |
| **Overlay quantity** (Compare) | which profile is overlaid for all models |
| **3D:** cut faces, core as sphere, animate, speed, auto-rotate, zoom | see [How the views are made](14-views.md) |
| **Cross-section:** animate convection cells | turns the cell animation on or off |

---
← Previous: [Quick start & usage](02-usage.md) · [Contents](README.md) · Next: [Results: global quantities](04-results.md) →
