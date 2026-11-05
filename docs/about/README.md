# How Star Models works

This folder is the in-app **About** page in Markdown, for review on GitHub. It describes:

* how the models are computed: equations, physics and numerical method;
* every input parameter and output quantity;
* how to use the app and the command-line tool;
* the accuracy, the limitations and the references.

GitHub typesets the equations (LaTeX, rendered with MathJax). In other Markdown viewers they appear as
LaTeX source.

## Contents

| # | Section | What it covers |
|---|---|---|
| 1 | [Overview](01-overview.md) | What a model is, the physics in use, what happens when you change a parameter |
| 2 | [Quick start & usage](02-usage.md) | Computing, evolving, comparing and exporting models; the tabs and their controls |
| 3 | [Input parameters](03-parameters.md) | Mass, composition, mixing length, advanced settings, presets, evolution and display settings |
| 4 | [Results: global quantities](04-results.md) | R, L, $`T_{\mathrm{eff}}`$, central values, convective zones, self-checks, the "Sun today" column |
| 5 | [Profile quantities & export](05-profiles-and-export.md) | Every column of the profile CSV, and the summary JSON |
| 6 | [Structure equations](06-structure-equations.md) | The four equations of stellar structure and their assumptions |
| 7 | [Microphysics](07-microphysics.md) | Equation of state, opacity, nuclear energy generation |
| 8 | [Energy transport & convection](08-energy-transport.md) | Radiative gradient, Schwarzschild criterion, mixing-length theory |
| 9 | [Boundary conditions](09-boundary-conditions.md) | Centre expansions and the grey Eddington atmosphere |
| 10 | [Numerical method](10-numerical-method.md) | Mesh, discretisation, Henyey relaxation, post-processing, validation |
| 11 | [Initial model & continuation](11-initial-model.md) | The polytrope cold start, warm starts, continuation, background solving |
| 12 | [Evolution in time](12-evolution.md) | Burning, mixing, $`\varepsilon_{\mathrm{grav}}`$, time steps, phases, stop conditions, checks |
| 13 | [Pulsations & oscillations](13-pulsations.md) | Fundamental radial mode, instability strip, solar-like oscillations |
| 14 | [How the views are made](14-views.md) | Colour, profiles, cross-section, 3D model, HR diagram, time evolution, compare |
| 15 | [Accuracy & limitations](15-limitations.md) | What the analytic physics gets right and wrong; what is not modelled; plans |
| 16 | [Constants & units](16-constants.md) | Physical and solar constants |
| 17 | [Command-line tool](17-cli.md) | `stellar-cli` options and examples |
| 18 | [References & software](18-references.md) | Textbooks, papers and libraries |

## Keeping this in sync with the app

The app's About tab (`crates/stellar-app/src/views/about.rs`) is the original; these files follow
it section by section. If you change one, change the other.

* **Equations:** the display equations are copied verbatim from
  `crates/stellar-app/equations/equations.json`, the same LaTeX the app typesets. A test
  (`cargo test -p stellar-app`) fails if any equation in that file is missing from these pages.
* **Live values:** the app shows some values live, such as the current parameter values. Here they
  appear as defaults and ranges only.
