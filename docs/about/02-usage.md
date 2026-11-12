# 2. Quick start & usage

## Computing a model

* Pick a preset or move the sliders in the left panel.
  * With **Auto-solve** on, a new model is computed 0.15 s after you stop editing; otherwise press
    **Solve**.
  * **Cancel** stops a running solve.
* The status bar at the bottom shows the progress: the stage, the Newton iteration and the size of its
  correction. When the solve is done, it shows the number of mesh points, the iterations and the time
  taken.
* If the left panel says *"Parameters changed; showing the previous model"*, the views do not yet show
  the current slider values: a solve is pending or failed. Errors are shown in red in the status bar.
* A note appears below 0.4 M☉, where the ideal-gas equation of state is no longer adequate.

## Evolving a star

* In the **Evolution** section of the left panel, choose where to stop and press **Evolve**. The stop
  can be:
  * the end of the main sequence;
  * an age;
  * as far as possible.
* The run is computed in the background. The HR track grows and the views follow the latest model
  ("follow latest").
* Use the timeline above the views to pick any computed age, step with **<** and **>**, or **Play**
  through the run.
* **Stop** interrupts a run, **Continue** resumes it, **Reset to ZAMS** discards it. **Sun's age**
  evolves to 4.57 Gyr.
* **Scale** sets whether axes move while you scrub or play. It is in the Evolution section of the left
  panel and in the timeline, and can be set before a run starts.
  * **fixed for run** (default) keeps every axis and colour scale constant over the whole run, so only
    the star changes.
  * **adapt to model** fits each view to the model on screen for the most detail, and zooms the HR
    diagram in on the track.
* **True size** draws the star to scale relative to its largest radius in the run, with a faint
  outline of that size, so you see it expand.
  * It is in the left panel for both views, and each view also has its own checkbox.
  * After a giant phase the early star is small; zoom in on the 3D view to see it.
* To see a pulsating star, evolve one of these, then select a model inside the shaded instability
  strip and open the 3D model:
  * 1.5–2 M☉: a δ Scuti star on the main sequence;
  * 5 M☉ "as far as possible": a Cepheid crossing.

## Comparing and exporting

* **Pin current model** keeps a copy for comparison.
  * You can pin up to 7 models, and each keeps its colour.
  * Pinned models are overlaid in Profiles and Compare and plotted in the HR diagram.
  * In Compare, untick a model to hide it or press **×** to unpin it.
* **Export** offers three files:
  * **Profile CSV**: one row per mesh point (see [Profile quantities & export](05-profiles-and-export.md));
  * **Summary JSON**: parameters, global results and convective zones;
  * **Screenshot PNG**: the whole window as currently shown.
* The buttons in the top-right corner switch between light, dark and system theme.

## Tabs

| Tab | What it shows |
|---|---|
| **Profiles (P, T, ρ)** | Pressure, temperature and density against radius $`r/R`$ or mass $`m/M`$. The three plots share the x axis and the hover cursor. |
| **More profiles** | Luminosity and enclosed mass, opacity, nuclear energy generation by process, temperature gradients, chemical composition, gravitational energy release, convective velocity, convective flux fraction and gas-pressure fraction, the propagation diagram ($`N`$, $`S_1`$, $`S_2`$) with the asymptotic spacings, and the sound speed. The x axis can also be $`\log P`$, which spreads out the thin outer envelope. |
| **Cross-section** | A slice through the star: temperature on the left half and density on the right. Convective zones carry animated cells. |
| **3D model** | A cut-away sphere with animated radiative and convective energy transport. It can be coloured by hydrogen content to show the helium core, and pulsates when the star is in the instability strip. |
| **HR diagram** | The model on a computed ZAMS line, with the evolutionary track of the current run, the instability strip, the present-day Sun and the pinned models. |
| **Time evolution** | After an evolution run, linked plots of how the global parameters change with age (or model number): $`L`$, $`R`$, $`T_{\mathrm{eff}}`$, $`\log g`$, central temperature and density, central composition, core masses, energy sources, pulsation period and solar-like oscillation frequencies. Click a plot to show that age everywhere. |
| **Compare** | A table of global properties for the current and pinned models, plus an overlay plot of any profile quantity. |

## Plot controls (all 2D plots)

| Action | Control |
|---|---|
| Pan | drag |
| Zoom | ⌘ + scroll (Ctrl + scroll on Linux/Windows), or pinch |
| Zoom into a box | right-drag |
| Reset | double-click |

* **Hover** for values: the tooltip names the model and gives the x coordinate and the de-logged y
  value with units.
* **Gray bands** mark the convective zones of the current model.

## 3D model controls

| Action | Control |
|---|---|
| Rotate | drag |
| Zoom towards the cursor (0.3× to 400×) | scroll, or pinch |
| Pan | right-drag or shift-drag |
| Reset | double-click |

* **Zoom to surface layers** magnifies the outer envelope at the cut, where thin convection zones are
  otherwise too small to see.
* **Cut faces** chooses temperature or density colouring.
* **Core as sphere** shows the core as a ball (as in many textbook illustrations) or cuts it like the
  rest of the star.
* **Animate**, **speed** and **auto-rotate** control the motion.
* Hover over the cut faces, the core ball or the surface to read the local state.

> Display settings (axes, colouring, zoom, animation) never change the model itself.

---
← Previous: [Overview](01-overview.md) · [Contents](README.md) · Next: [Input parameters](03-parameters.md) →
