# Web version: design notes

Status: **ideas only, not implemented.** This note records what a browser version of Star Models
would take, the options for running the computation, and the decisions taken so far.

## Short answer

Everything runs **in the browser, on the client**. The physics and solver crate (`stellar-core`) and
the egui app (`stellar-app`) both compile to WebAssembly (`wasm32-unknown-unknown`). The web version
is a set of static files (HTML, `.wasm`, `.js`):

* no server-side code;
* no server costs;
* nothing the user computes leaves their machine.

Estimated effort:

| Approach | Developer work |
|---|---|
| Option A: time-sliced on the main thread | about 2.5–3.5 days |
| Option B: Web Worker | about 1–2 days more than A |

The native desktop app stays as it is; both builds come from the same code.

## What blocks a web build today

Most of the app is pure Rust and compiles to WASM unchanged:

* the physics and the solver;
* `egui_plot` plots;
* the SVG equations (`egui_extras` + `resvg`);
* the CPU 3D renderer;
* PNG encoding.

These are the places that need work:

| Where | What | Web replacement |
|---|---|---|
| `crates/stellar-app/src/worker.rs`: `Worker::submit`, `Worker::submit_zams`, `Worker::submit_evolution` | Each runs on a native thread (`std::thread::spawn`): structure solves, the HR-diagram ZAMS line, evolution runs. A standard WASM build has no threads. | Option A or B below |
| `stellar-core`: `solve_model` / `continuation` (`structure/solve.rs`), `Evolution::step` (`structure/evolve.rs`), the Newton loop in `solve_structure` | Run to completion in one call. On the browser's main thread a long run would freeze the page. | Yield points (A) or a worker (B) |
| `crates/stellar-app/src/worker.rs` | `std::time::Instant` panics in the browser | the `web-time` crate (same API, works on both) |
| `crates/stellar-app/src/app.rs`: `save_text`, `write_png`, Screenshot PNG | Native file dialogs (`rfd::FileDialog`) and `std::fs` writes | Browser downloads: build a `Blob`, click a temporary `<a download>` via `web-sys`; encode the PNG in memory |
| `crates/stellar-app/src/app.rs`: `StarApp::new` | The `STAR_MODELS_*` capture/test environment variables | Native only: compile out with `cfg(not(target_arch = "wasm32"))` |
| `crates/stellar-app/src/main.rs` | Native window options (`NativeOptions`, size, `multisampling: 4`) | A `#[cfg(target_arch = "wasm32")]` entry point using `eframe::WebRunner` on a `<canvas>`, with the web anti-aliasing option |
| Fonts | On the desktop, symbols the bundled egui fonts lack (☉, ∇, ⌘, ≈, ≪, √, superscript and subscript digits) may come from system fonts. A browser has no system fonts. | Test coverage in a browser; if needed, bundle a fallback font (e.g. a subset of Noto Sans Math / Noto Sans Symbols, ~0.2–0.6 MB) via `egui::FontDefinitions` |
| `StarModel.solution`, `StarModel.l_scale`, `RadialMode.dp` (`#[serde(skip)]`) | Not serialised today | Only matters for option B, where whole models cross the worker boundary |

## Option A: time-sliced main thread (simplest)

The solver and the UI share the browser's main thread. The solver gives control back between Newton
iterations, so the page keeps drawing.

* **Core change:** make the drivers `async` (`solve_model`, `continuation`, `cold_start`,
  `solve_structure`, `Evolution::step`). Add a yield point after each Newton iteration through a new
  `SolveControl::yield_now()` that returns a future. On the desktop that future is always ready.
  * One Newton iteration costs about 3–4 ms on the desktop (the solar model needs 13 iterations in
    40–60 ms) and is expected to take about 5–8 ms in WASM.
* **Desktop:** keep the worker threads and drive the futures with a tiny `block_on` (e.g.
  `pollster`). Keep synchronous wrappers (`solve_model_blocking`, …) so the CLI and the tests don't
  change.
* **Web:** run the futures with `wasm_bindgen_futures::spawn_local`. `yield_now()` waits for a
  `setTimeout(0)` (e.g. `gloo-timers`), so the browser can paint between slices.
  * The existing `Msg` channel (`std::sync::mpsc`) and the `AtomicBool` cancel flags work
    unchanged on a single thread.
  * Several jobs (a model solve and the ZAMS line) interleave naturally.

Pros:

* one bundle, and any static host works;
* Stop / Cancel take effect at the next slice;
* the smallest change to the code.

Cons:

* the frame rate dips while computing;
* runs are about 1.5–2× slower than on the desktop. For example, the full 5 M☉ "as far as possible"
  run (about 1650 steps, 15 s native) would take roughly 40–60 s;
* computation pauses or slows down when the tab is in the background (browser timer throttling).

## Option B: Web Worker

The solver runs in a dedicated browser worker, a separate WASM instance; the page only draws.

* **Bundling:** a second entry point, e.g. with `gloo-worker` or a Trunk `data-type="worker"`
  target.
* **Messages:** requests (solve, ZAMS line, evolve, cancel) and results (progress, models, evolved
  models, done) as serde messages. Models must be serialised whole, including the fields that are
  skipped today. A model is about 200 kB, and an evolution run streams one per step.
* **Cancelling:** a worker busy in a synchronous loop cannot read a "cancel" message. Either
  terminate and restart the worker (about 0.1 s), or use Option A's async drivers inside the worker,
  which gives it yield points.

Pros:

* the UI stays at full frame rate during long runs;
* runs keep going in background tabs;
* real parallelism between computing and drawing.

Cons:

* more plumbing: two entry points and message types;
* models are copied through `postMessage`;
* cancellation needs care.

A natural path is **A first, then B if needed**, reusing A's async drivers inside the worker.

**Not recommended:** WASM threads with `SharedArrayBuffer` (`wasm-bindgen-rayon` style).

* They need COOP/COEP HTTP headers, which simple static hosts such as GitHub Pages can't set.
* They need a nightly `build-std` toolchain.

## Work common to both options

* [ ] Add a `lib` target to `stellar-app` and a `#[cfg(target_arch = "wasm32")]` start function
      using `eframe::WebRunner`.
* [ ] Add `index.html` with a `<canvas>` and a Trunk configuration. `trunk serve` gives a local
      preview; `trunk build --release` writes the static files to `dist/`.
* [ ] Add web-only dependencies under `[target.'cfg(target_arch = "wasm32")'.dependencies]`:
      `wasm-bindgen`, `wasm-bindgen-futures`, `web-sys` (Blob, Url, HtmlAnchorElement), `web-time`,
      and `gloo-timers` for option A.
* [ ] Compile out `rfd`, `std::fs` and the `STAR_MODELS_*` variables on the web.
* [ ] Export CSV, JSON and PNG as browser downloads.
* [ ] Check font coverage in a browser and add a fallback font if needed.
* [ ] CI: add a `wasm32-unknown-unknown` build job (build only).
* [ ] README: add web build and preview instructions.

## Hosting (decided)

Start with a **local build and preview**:

```sh
rustup target add wasm32-unknown-unknown
cargo install trunk
trunk serve            # http://127.0.0.1:8080
trunk build --release  # static files in dist/
```

The `dist/` folder can later go to any static host (GitHub Pages, Netlify, Cloudflare Pages, a plain
web server). Options A and B need no server-side code and no special HTTP headers.

## Expectations and risks

* **Download size:** about 3–5 MB compressed, including the 61 embedded equation SVGs.
  `wasm-opt` and brotli/gzip help; it can be cached after the first load.
* **Browsers:** recent Chrome, Firefox and Safari through eframe's wgpu backend (WebGPU, or WebGL2
  as fallback). Safari needs the most testing.
* **Speed:** WASM numerics run about 1.3–2× slower than native. A single model takes 0.1–0.3 s; long
  evolution runs take tens of seconds (with option A, see above).
* **Memory:** timeline snapshots are already thinned (well under 100 MB), within WASM limits.
* **Clipboard:** "Copy LaTeX" needs a secure context: https, or localhost when previewing.
* **3D view:** the CPU renderer draws about 40 000 triangles per frame (about 2 ms natively). In WASM
  it should still fit a frame, but it competes with option A's solver slices.

## Verification (when implemented)

* `trunk serve`, then in Chrome, Firefox and Safari:
  * every tab renders, with no missing glyphs and the equations typeset;
  * the 3D animation and pulsation run;
  * the Sun to 4.57 Gyr and 5 M☉ "as far as possible" runs work, with Stop / Continue;
  * the timeline and the Time evolution tab work, including click-to-select;
  * every export arrives as a download.
* The native app and `cargo test --workspace` are unchanged.

## Open decision

**Option A or B** is to be chosen when this is implemented. The recommendation is A first: it's
smaller and keeps one bundle. Move to B if the UI must stay at full frame rate during long evolution
runs.
