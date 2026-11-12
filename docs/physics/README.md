# Star Models: physics reference

This folder is the single reference for the physics of Star Models: what the code computes today, and
the upgrades agreed in the review of 2026-11-10. Each page covers one topic and has the same two parts:

* **Implemented**: the equations as they are in `crates/stellar-core` on `main`.
* **Planned**: the equations of each upgrade, with its ID from
  [`../physics-review.md`](../physics-review.md) and its phase in
  [`../implementation-plan.md`](../implementation-plan.md).

## Contents

| # | Page | Implemented | Planned |
|---|---|---|---|
| 1 | [Structure equations](01-structure-equations.md) | four equations in mass coordinate | neutrino term, rotation factor, variable mass |
| 2 | [Equation of state](02-equation-of-state.md) | ideal gas + radiation, fully ionised | degeneracy, ionisation, H₂ |
| 3 | [Opacity](03-opacity.md) | analytic Kramers + H⁻ + scattering | tables, conduction |
| 4 | [Nuclear burning](04-nuclear-burning.md) | pp, CNO, triple-α, weak screening | pp branching, newer fits, C/O, neutrinos |
| 5 | [Convection and mixing](05-convection-and-mixing.md) | Schwarzschild, mixing-length theory, instant mixing | overshoot, Ledoux, diffusion, settling |
| 6 | [Boundary conditions](06-boundary-conditions.md) | centre expansion, grey Eddington photosphere | integrated atmosphere, non-grey $`T(\tau)`$ |
| 7 | [Evolution in time](07-evolution.md) | burning, $`\varepsilon_{\mathrm{grav}}`$, time steps, phases | burning fixes, mass loss, pre-main sequence |
| 8 | [Rotation](08-rotation.md) | none | shape, gravity darkening, structure, mixing |
| 9 | [Pulsations](09-pulsations.md) | radial fundamental mode, strip, scaling relations | buoyancy frequency, asymptotics, non-radial modes |
| 10 | [Numerical method and checks](10-numerics-and-checks.md) | mesh, Henyey solver, starting model, self-checks | redshift output |
| 11 | [Validity, limits and order of work](11-validity-and-roadmap.md) | measured accuracy | what is deliberately left out |
| 12 | [Constants, symbols and references](12-constants-and-references.md) | | |

## Status markers

| Marker | Meaning |
|---|---|
| **Implemented** | In the code on `main`; the equation was compared with the source file named next to it |
| **Planned (ID)** | Not in the code. The ID (for example D14, E1) points to the review and the plan |
| **[check source]** | A planned formula or coefficient quoted from memory of the cited book or paper. Verify it against the source before coding it |
| **[compared]** | A planned formula that an external review (2026-11-11) compared with a named source; the comparison is listed in [Constants, symbols and references](12-constants-and-references.md). Whoever codes it still records the equation number in the code |
| **[measured]** | A number obtained by running the existing release build of `stellar-cli` |

These pages went through four rounds of external review on 2026-11-11 (theory only; the code was
not part of it). The source comparisons made there cover the listed items only and are not a
blanket verification of everything marked **[check source]**.

Task 0.3 (2026-11-11 and 2026-11-12) worked through every remaining marker: most were compared with
the original, and the rest were settled by a recorded decision where the original was not
available (a verified substitute source, independent reproductions, or dropping the item from the
plan). No marker is left. Which items rest on a decision rather than the original is listed at the
end of [page 12](12-constants-and-references.md). New planned formulas quoted from memory get the
marker again.

## How this relates to the other documents

* [`../about/`](../about/README.md) mirrors the in-app About page. It describes the app as it is
  (usage, views, parameters, physics) and must stay in sync with
  `crates/stellar-app/src/views/about.rs`. It contains only implemented physics.
* This folder contains the physics only, implemented and planned together, so that an upgrade can be
  read next to what it replaces.
* [`../physics-review.md`](../physics-review.md) explains why each upgrade is proposed and how large
  the effect is. [`../implementation-plan.md`](../implementation-plan.md) says in what order and in
  which files.
* [`../execution-plan.md`](../execution-plan.md) is the plan to work from: one task per planned
  block of these pages, with files, tests and the ticket number.

When an upgrade is implemented: move its block from "Planned" to "Implemented" here, and update the
About page as the plan describes.

## Conventions

* CGS units throughout; solar units only for display.
* The independent variable is the enclosed mass $`m`$, or $`q = m/M`$.
* $`T_6 = T/10^6`$ K, $`T_8 = T/10^8`$ K, $`T_9 = T/10^9`$ K.
* Equations use GitHub's LaTeX syntax: display equations in `math` code blocks, inline ones as
  `` $`…`$ ``.
* Symbols are listed in [Constants, symbols and references](12-constants-and-references.md).

> **Code note.** The implemented parts describe `main` at commit `954b6d3`, which includes the
> former `quick-wins` branch (diagnostics and tests, no change to the physics equations given here).
> Every **[measured]** number was re-run on that commit (task 0.1, 2026-11-11) and none moved.
