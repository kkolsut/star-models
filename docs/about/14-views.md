# 14. How the views are made

## Colour

* **Lines and points** use a colourblind-safe categorical palette, validated for colour-vision
  deficiencies.
  * The current model is always slot 1 (blue).
  * Each pinned model keeps the slot it received when pinned, so colours never move between models.
  * A legend is always shown; the Compare table lists exact values.
* **Gray bands** mark the convective zones of the current model. Gray is a neutral colour, never a
  series colour.
* **Scale:** during an evolution run, the timeline's Scale switch sets whether axes and colour scales
  stay fixed for the whole run (the default) or adapt to the shown model. See
  [Evolution in time](12-evolution.md#in-the-app).
* **Continuous fields** use single-hue, logarithmic colour scales: orange for temperature, blue for
  density.
  * In 2D views the low end fades towards the background.
  * In 3D, hotter or denser is always brighter.

## Profiles and More profiles

* Logarithmic quantities are plotted as $`\log_{10}`$, with ticks labelled $`10^n`$. When the range is
  narrow, $`2\cdot10^n`$ and $`5\cdot10^n`$ are labelled too.
* Energy generation is hidden where it falls below $`10^{-8}`$ of its peak. Convective velocity is
  drawn only in convective layers.
* The $`\log P`$ axis puts the centre on the left. It helps because the outer 10% of the radius spans
  almost half of the pressure range (46% for the Sun), and the outer 1% still spans about a fifth.

## Cross-section

A disc split in two halves, coloured by $`\log T`$ (left) and $`\log\rho`$ (right), between the
surface and central values.

* Convective zones are outlined and carry animated rolls whose direction alternates, as in real
  convection cells.
* The hover read-out is interpolated from the profile.

## 3D model

* **Geometry:** the sphere has the octant $`x, y, z > 0`$ removed. The three exposed quarter-discs are
  coloured by $`T`$ or $`\rho`$ from the radial profile.
* **Surface:** the approximate colour of a black body at $`T_{\mathrm{eff}}`$. When the envelope is
  convective, it is mottled by an animated granulation pattern.
* **Core ball:**
  * the convective core, if the star has one smaller than the whole star;
  * otherwise, the sphere that produces 90% of the luminosity.

  It "boils" when convective.
* **Radiative zones:** 180 photon sparks random-walk outward. The number per zone is proportional to
  its face area. Short trails show the zig-zag path of radiative diffusion.
* **Convective zones:** rows of circulating cells about $`0.16\,R`$ deep, at most 40 per row and face.
  Parcels rising (moving outward) are amber; sinking parcels are blue.
* **White arcs:** boundaries between radiative and convective zones.
* **Cut-face colouring by hydrogen $`X`$** (linear scale): hydrogen-rich gas is dark and hydrogen-poor
  gas is light. During evolution the burnt core shows up as a light region (or a light core ball).
* **Pulsation** (models in the instability strip): every layer is drawn at radius
  $`r\,\bigl(1 + A\,\xi(r)\sin\varphi\bigr)`$.
  * $`\xi(r)`$ is the computed fundamental-mode shape, equal to 1 at the surface.
  * $`A = 6\%`$, and one cycle takes 1.6 s.
  * The surface brightens near maximum compression.
  * The legend gives the real period.
* **Rendering:** orthographic projection on the CPU.
  * About 40,000 shaded triangles are depth-sorted every frame (painter's algorithm, about 2 ms).
  * They are drawn as one mesh with 4× multisampling.
  * Overlays use an exact ray test against the cut-away geometry, so they are hidden correctly behind
    the sphere and the core.

> The cross-section and the 3D view are visualisations of a **1D** model. Photon paths, convection
> cells, granulation and boiling are illustrative: their size, shape and speed are not computed, and
> this is not a 3D radiation-hydrodynamics simulation. Animation speeds are illustrative, not to scale. In the Sun a photon needs of order $`10^4`$–$`10^5`$
> years to diffuse out of the radiative zone. Convective turnover takes from minutes (granules) to about
> a month (base of the envelope).

## HR diagram

* **ZAMS line:** connects models at 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6.5, 8,
  10, 13, 16, 20, 25, 30 and 40 M☉.
  * They are computed in the background with the current $`X`$, $`Z`$, $`\alpha_{\mathrm{MLT}}`$ and
    $`f_{\mathrm{CNO}}`$.
  * Each has 400 mesh points and is reached by continuation from the current model.
* **Axes:** $`T_{\mathrm{eff}}`$ increases to the left, as is conventional.
* **Sun:** the ☉ symbol marks the present-day Sun (5772 K, 1 L☉).
* **Evolutionary track:** joins every step of the current run.
  * Labels mark the ZAMS, the end of the main sequence and the latest age.
  * Hovering over the track shows age, phase, $`X_c`$, radius and pulsation.
* **Shaded band:** the approximate classical instability strip (see
  [Pulsations & oscillations](13-pulsations.md)).

## Time evolution

The plots have one point per computed step, not only the timeline models.

* **Energy-source panel:** shows shares of the surface luminosity.
  * The nuclear shares are the pp/CNO/3α fractions times $`(1 - L_{\mathrm{grav}}/L)`$.
  * The gravitational share is $`\int \varepsilon_{\mathrm{grav}}\,dm \,/\, L`$, so the four add up
    to 1. It can be negative while the star expands.
* **Legends** are drawn next to each panel title.
* **Logarithmic axes** label round values when the range is less than a decade.

## Compare

A table of global properties for the current and pinned models, plus an overlay of one chosen profile
quantity on the selected x axis.

---
← Previous: [Pulsations & oscillations](13-pulsations.md) · [Contents](README.md) · Next: [Accuracy & limitations](15-limitations.md) →
