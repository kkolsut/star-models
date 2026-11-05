# Star Models

Interactive **zero-age main-sequence (ZAMS) stellar structure** models, written in Rust.

You choose the mass, composition (X, Y = 1 − X − Z, Z) and the mixing-length parameter. The app
solves the four stellar structure equations, with both **radiative and convective energy transport**,
and shows pressure, temperature and density as functions of radius. It also shows a cross-section of
the star, an HR diagram and side-by-side comparisons. Defaults are the Sun (M = 1 M☉, X = 0.70,
Z = 0.02, α_MLT = 1.8).

Runs natively on **Linux and macOS**, and builds on Windows too (egui/eframe).

## Quick start

```sh
# Rust toolchain (if needed)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

cargo run --release                 # the GUI (default member)
cargo run --release -p stellar-cli  # headless: prints a solar model summary
cargo test --workspace              # unit + integration tests
```

Linux system packages for the GUI (Debian/Ubuntu):

```sh
sudo apt-get install -y libxkbcommon-dev libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libwayland-dev
```

File dialogs on Linux use the XDG desktop portal (`xdg-desktop-portal` must be running, which it is
on standard GNOME/KDE desktops).

## The app

| Tab | What it shows |
|---|---|
| **Profiles (P, T, ρ)** | log P, log T, log ρ against r/R or m/M, with linked axes and a shared cursor. Convective zones are shaded. Pinned models are overlaid. |
| **More profiles** | L/L_s and m/M, opacity κ, ε split into pp, CNO and 3α, ∇_rad vs ∇_ad vs ∇ (shows *why* a layer is convective), convective velocity, F_conv/F and β = P_gas/P. The x axis can be r/R, m/M or **log P**; log P spreads out the thin outer envelope. |
| **Cross-section** | Split disc with temperature on the left and density on the right (log color scales). Convective zones are outlined and carry animated convection cells. Hover to read the local state. |
| **3D model** | A cut-away sphere: one octant removed, cut faces coloured by T or ρ, and the energy-generating (or convective) core as a ball in the corner. The surface has the black-body colour at T_eff. Animated transport: photon sparks random-walk outward through radiative zones, convection cells carry rising hot (amber) and sinking cool (blue) gas, the surface shows granulation for a convective envelope, and a convective core "boils". Drag to rotate, scroll or pinch to zoom at the cursor (up to 400×, e.g. **Zoom to surface layers** for the Sun's thin envelope), right- or shift-drag to pan, double-click to reset. Hover to read the local state. |
| **HR diagram** | The current model on a ZAMS line computed for the current composition, alongside the present-day Sun and the pinned models. |
| **Compare** | Table of global properties for the current and pinned models, plus an overlay of any profile quantity. |

The left panel holds the parameters (auto-solve while you edit), the global results next to
present-day solar values, **Pin current model**, and export buttons: **profile CSV** (25 columns),
**summary JSON** and a **PNG screenshot**.

When you change a parameter, the solver warm-starts from the model on screen, so a new solution
usually takes 20–100 ms.

### CLI

```sh
stellar-cli --mass 1 --x 0.70 --z 0.02 -o sun.csv --json sun.json
stellar-cli --sweep 0.5,1,2,5,10,20      # a ZAMS sequence by continuation
stellar-cli --mass 3 -v                  # print every Newton iteration
```

## Physics and numerics

Code layout: `crates/stellar-core` holds the physics and solver (no GUI dependencies),
`crates/stellar-app` the egui app, and `crates/stellar-cli` the command-line tool.

* **Structure equations** in mass coordinates, with unknowns (ln r, ln P, ln T, L) per mesh point:
  mass continuity, hydrostatic equilibrium, energy generation, and energy transport with
  d ln T / d ln P = ∇.
* **Energy transport:** ∇ = ∇_rad in radiative layers. Where the Schwarzschild criterion
  ∇_rad > ∇_ad holds, ∇ comes from **mixing-length theory** (Kippenhahn, Weigert & Weiss cubic,
  l = α_MLT H_P). It is ≈ ∇_ad deep inside and superadiabatic near the surface.
* **EOS:** fully ionized ideal gas plus radiation pressure (β, ∇_ad, c_P, δ).
* **Opacity:** analytic. Electron scattering + Kramers free–free and bound–free (with guillotine
  factor) + H⁻ + a molecular floor.
* **Nuclear:** analytic pp chain, CNO cycle and triple-α rates with weak screening.
* **Boundary conditions:** central expansions at m₀ = 10⁻⁶ M, and a grey Eddington photosphere at
  τ = ⅔ at the surface.
* **Solver:** Henyey relaxation, i.e. Newton–Raphson on the full discretized system.
  * The Jacobian is built from local finite differences and solved with a hand-written banded LU.
  * Corrections are damped.
  * Continuation in (log M, X, Z, α) is used for robustness.
  * The initial model is an n = 3 polytrope scaled by homology.
  * The mesh has 800 points, logarithmic in q near the centre and logarithmic in 1 − q down to
    ~10⁻¹³ near the surface.
* **Self-checks** reported for every model: virial theorem ∫3P/ρ dm / ∫Gm/r dm ≈ 1 and energy
  conservation ∫ε dm / L ≈ 1. Both hold to better than 10⁻³, typically ~10⁻⁵ to 10⁻⁴.

### Accuracy: what to expect

The analytic microphysics gets the **qualitative** structure right:

* radiative core and convective envelope in the Sun
* convective cores above ~1.2 M☉
* CNO dominance above ~1.5 M☉
* fully convective stars at low mass
* massive stars (where electron scattering dominates the opacity) close to detailed ZAMS models

It is **not quantitatively solar**. With X = 0.70 and Z = 0.02 the 1 M☉ model has L ≈ 0.42 L☉,
R ≈ 1.11 R☉ and a thin convective envelope (base at 0.97 R). A detailed ZAMS Sun has L ≈ 0.7 L☉,
R ≈ 0.89 R☉ and its envelope base near 0.71 R. Two things cause this:

* The textbook Kramers opacity is too high in the core and far too low at 10⁵–10⁶ K.
* The EOS has no H/He ionization zones.

With Z = 0.008 the same code reproduces the classic Carroll & Ostlie "StatStar" Sun
(L ≈ 0.83 L☉, T_eff ≈ 5500 K), which confirms the implementation matches the textbook physics.
Also, a ZAMS model is the *young* Sun; the present Sun has evolved for 4.6 Gyr.

## Roadmap

* **Tabulated microphysics:** OPAL/OP opacity tables plus Ferguson low-temperature tables,
  an EOS with Saha H/He ionization, and better reaction rates. These plug into the existing
  `Eos` / `Opacity` / `Nuclear` traits. The solver needs no analytic derivatives.
* **Evolution:** per-shell composition X(t), ε_grav, Ledoux criterion and an age parameter
  (e.g. the Sun at 4.57 Gyr). The Henyey solver and the per-shell composition storage are already
  in place for this.
* **Web build:** eframe compiles to WASM. The Newton solver already exposes a step-wise API
  (`HenyeyNewton::step`) so a browser build can iterate per frame instead of on threads.

## Development notes

* `STAR_MODELS_CAPTURE=<dir>` runs the app non-interactively. It visits every tab, saves a PNG of
  each, then exits. Optional settings for this mode:
  * `STAR_MODELS_THEME=light|dark`
  * `STAR_MODELS_MASS=<M☉>`
  * `STAR_MODELS_PIN=0.5,3` (pre-pinned models)
  * `STAR_MODELS_AXIS=radius|mass|pressure`
  * `STAR_MODELS_TABS=model3d,hr` (only capture these tabs)
  * `STAR_MODELS_3D_ZOOM=surface` (start the 3D view zoomed on the envelope)
* The 3D view is rendered on the CPU: rotation, shading and a depth sort (painter's algorithm) of about
  40k triangles per frame (~2 ms), handed to egui as one mesh with 4× MSAA. Overlays use an exact
  ray test against the cut-away geometry for occlusion. It needs no GPU-specific code, so it carries
  over to a WASM build.
* Chart colors use a validated colorblind-safe categorical palette (`crates/stellar-app/src/theme.rs`).
  The current model is always slot 1, and each pinned model keeps its slot. Sequential color
  scales are single-hue ramps.
