//! "About" tab: how the models are computed, what every parameter and output
//! means, how to use each view, limitations and references.
//!
//! Numbers quoted here come from the code (stellar-core). Where a value is
//! available at run time (current parameters, physics in use, ZAMS masses,
//! CSV columns, version) it is read live, so the text cannot drift from the
//! implementation.
//!
//! The same text is in `docs/about/` as Markdown, for reading on GitHub. Keep
//! both in sync; a test checks that every equation appears there.

use egui::{Align, RichText, Ui};
use stellar_core::output::CSV_COLUMNS;
use stellar_core::{Physics, StarParams};

use crate::theme;
use crate::worker::ZAMS_MASSES;

#[derive(Default)]
pub struct AboutState {
    scroll_to: Option<usize>,
    /// Frames left to keep re-issuing the scroll request: the content height
    /// is only known after the first layout pass.
    scroll_frames: u8,
    current: usize,
}

impl AboutState {
    /// Open the page scrolled to section `i` (used by the capture mode).
    pub fn at_section(i: usize) -> Self {
        Self { scroll_to: Some(i.min(SECTIONS.len() - 1)), scroll_frames: 4, current: 0 }
    }
}

const SECTIONS: &[&str] = &[
    "Overview",
    "Quick start & usage",
    "Input parameters",
    "Results: global quantities",
    "Profile quantities & export",
    "Structure equations",
    "Microphysics",
    "Energy transport & convection",
    "Boundary conditions",
    "Numerical method",
    "Initial model & continuation",
    "Evolution in time",
    "Pulsations & oscillations",
    "How the views are made",
    "Accuracy & limitations",
    "Constants & units",
    "Command-line tool",
    "References & software",
];

// ---------------------------------------------------------------------------
// Small typographic helpers

fn heading(ui: &mut Ui, text: &str) -> egui::Response {
    ui.add_space(18.0);
    let r = ui.label(RichText::new(text).size(22.0).strong());
    ui.separator();
    r
}

fn sub(ui: &mut Ui, text: &str) {
    ui.add_space(8.0);
    ui.label(RichText::new(text).size(15.5).strong());
}

fn p(ui: &mut Ui, text: &str) {
    ui.add(egui::Label::new(text).wrap());
}

fn muted(ui: &mut Ui, dark: bool, text: &str) {
    ui.add(egui::Label::new(RichText::new(text).color(theme::secondary(dark))).wrap());
}

/// An equation (or code line), indented and in monospace. A trailing
/// annotation, separated by three or more spaces and containing no '=', is
/// shown in secondary text after it. Both wrap, so a line never widens the page.
fn eq(ui: &mut Ui, text: &str) {
    let (formula, note) = match text.rfind("   ") {
        Some(i) if !text[i..].contains('=') => (text[..i].trim_end(), text[i..].trim()),
        _ => (text, ""),
    };
    let dark = ui.visuals().dark_mode;
    ui.horizontal_wrapped(|ui| {
        ui.add_space(22.0);
        ui.add(egui::Label::new(RichText::new(formula).monospace()).wrap());
        if !note.is_empty() {
            ui.add_space(12.0);
            ui.add(egui::Label::new(RichText::new(note).color(theme::secondary(dark))).wrap());
        }
    });
}

/// Display size of the typeset equations relative to their natural size.
const MATH_SCALE: f32 = 0.92;

/// A typeset equation (pre-rendered from LaTeX, see `equations/equations.json`),
/// tinted with the text colour, with an optional annotation after it.
/// Right-click offers "Copy LaTeX".
fn math(ui: &mut Ui, id: &str, note: &str) {
    let Some(e) = crate::equations::get(id) else {
        eq(ui, &format!("[missing equation {id}]"));
        return;
    };
    let dark = ui.visuals().dark_mode;
    ui.add_space(3.0);
    ui.horizontal_wrapped(|ui| {
        ui.add_space(22.0);
        let max_w = (ui.available_width() - 8.0).max(60.0);
        let scale = MATH_SCALE.min(max_w / e.width);
        let image = egui::Image::from_bytes(format!("bytes://equations/{}.svg", e.id), e.svg)
            .tint(theme::primary(dark))
            .fit_to_exact_size(egui::vec2(e.width, e.height) * scale)
            .sense(egui::Sense::click());
        let response = ui.add(image);
        response.context_menu(|ui| {
            if ui.button("Copy LaTeX").clicked() {
                ui.ctx().copy_text(e.tex.to_string());
                ui.close();
            }
        });
        if !note.is_empty() {
            ui.add_space(14.0);
            ui.label(RichText::new(note).color(theme::secondary(dark)));
        }
    });
    ui.add_space(3.0);
}

fn bullet(ui: &mut Ui, text: &str) {
    ui.horizontal_top(|ui| {
        ui.add_space(8.0);
        ui.label("•");
        ui.add(egui::Label::new(text).wrap());
    });
}

/// A bullet whose first part is emphasised ("term — explanation").
fn term(ui: &mut Ui, name: &str, text: &str) {
    ui.horizontal_top(|ui| {
        ui.add_space(8.0);
        ui.label("•");
        ui.add(egui::Label::new(RichText::new(name).strong()));
        ui.add(egui::Label::new(text).wrap());
    });
}

/// A parameter block: name, metadata line, description.
fn param(ui: &mut Ui, dark: bool, name: &str, meta: &str, text: &str) {
    ui.add_space(6.0);
    ui.label(RichText::new(name).strong().size(14.5));
    ui.label(RichText::new(meta).monospace().color(theme::secondary(dark)));
    p(ui, text);
}

// ---------------------------------------------------------------------------

pub fn show(ui: &mut Ui, state: &mut AboutState, params: &StarParams, physics: &Physics) {
    let dark = ui.visuals().dark_mode;

    egui::Panel::left("about_toc").resizable(false).exact_size(230.0).show(ui, |ui| {
        ui.add_space(6.0);
        ui.label(RichText::new("Contents").strong());
        ui.add_space(4.0);
        for (i, name) in SECTIONS.iter().enumerate() {
            if ui.selectable_label(state.current == i, *name).clicked() {
                state.scroll_to = Some(i);
                state.scroll_frames = 3;
            }
        }
        ui.add_space(12.0);
        ui.label(
            RichText::new(format!("Star Models v{}", env!("CARGO_PKG_VERSION"))).small().color(theme::muted(dark)),
        );
    });

    // Top of the visible scroll viewport, in screen coordinates.
    let viewport_top = ui.cursor().top();
    egui::ScrollArea::vertical().auto_shrink(false).id_salt("about_scroll").show(ui, |ui| {
        // Comfortable line length on wide windows, full width on narrow ones.
        ui.set_max_width(ui.available_width().min(980.0));
        let mut current = 0;
        let mut section = |ui: &mut Ui, i: usize, body: &mut dyn FnMut(&mut Ui)| {
            let r = heading(ui, SECTIONS[i]);
            if state.scroll_to == Some(i) {
                // Jump straight there: long animated scrolls feel sluggish in a long document.
                r.scroll_to_me_animation(Some(Align::TOP), egui::style::ScrollAnimation::none());
                state.scroll_frames = state.scroll_frames.saturating_sub(1);
                if state.scroll_frames == 0 {
                    state.scroll_to = None;
                }
                ui.ctx().request_repaint();
            }
            if r.rect.top() <= viewport_top + 60.0 {
                current = i;
            }
            body(ui);
        };

        section(ui, 0, &mut |ui| overview(ui, dark, physics));
        section(ui, 1, &mut |ui| usage(ui, dark));
        section(ui, 2, &mut |ui| parameters(ui, dark, params));
        section(ui, 3, &mut |ui| results(ui, dark));
        section(ui, 4, &mut |ui| profiles(ui, dark));
        section(ui, 5, &mut |ui| equations(ui));
        section(ui, 6, &mut |ui| microphysics(ui, dark));
        section(ui, 7, &mut |ui| transport(ui, dark));
        section(ui, 8, &mut |ui| boundaries(ui, dark));
        section(ui, 9, &mut |ui| numerics(ui, dark));
        section(ui, 10, &mut |ui| initial(ui, dark));
        section(ui, 11, &mut |ui| evolution(ui, dark));
        section(ui, 12, &mut |ui| pulsation(ui, dark));
        section(ui, 13, &mut |ui| views(ui, dark));
        section(ui, 14, &mut |ui| limitations(ui, dark));
        section(ui, 15, &mut |ui| constants(ui, dark));
        section(ui, 16, &mut |ui| cli(ui, dark));
        section(ui, 17, &mut |ui| references(ui, dark));
        ui.add_space(40.0);
        state.current = current;
    });
}

// ---------------------------------------------------------------------------
// Sections

fn overview(ui: &mut Ui, dark: bool, physics: &Physics) {
    p(
        ui,
        "Star Models computes the internal structure of a star on the zero-age main sequence (ZAMS). \
         Such a star is chemically homogeneous, has just started burning hydrogen in its core, and is in \
         hydrostatic and thermal equilibrium. You choose the mass and the composition. The app then solves \
         the four differential equations of stellar structure from the centre to the photosphere, with \
         energy carried outward by radiation or by convection. It shows how pressure, temperature, density \
         and many other quantities vary through the star.",
    );
    p(
        ui,
        "Models are spherically symmetric, non-rotating and non-magnetic, without mass loss. The defaults are those \
         of the Sun (1 M☉, X = 0.70, Z = 0.02, α_MLT = 1.8). The default model is therefore the young Sun as it was \
         about 4.6 billion years ago. Press Evolve to follow the star in time: hydrogen burns to helium, the core \
         contracts, the envelope expands, and the star moves along an evolutionary track in the HR diagram. Stars \
         that cross the instability strip pulsate.",
    );
    muted(
        ui,
        dark,
        "What \"ZAMS\" means here: an idealised ZAMS, chemically homogeneous and static (ε_grav = 0), so nuclear burning \
         alone supplies the luminosity. Detailed evolution codes usually reach the ZAMS along a pre-main-sequence track and \
         define it by a criterion such as L_nuc/L ≈ 0.99. The two agree closely, but compare models at matched \
         evolutionary stages: a ZAMS model is not the present-day Sun.",
    );
    sub(ui, "Physics in use");
    p(ui, &physics.describe());
    sub(ui, "What happens when you change a parameter");
    bullet(ui, "The parameters are checked: X + Z ≤ 1, and every value must be within its range.");
    bullet(
        ui,
        "A background thread solves the structure equations with a Newton–Raphson relaxation method (Henyey method). \
         It starts from the model currently on screen and moves step by step to the new parameters.",
    );
    bullet(
        ui,
        "The converged solution is evaluated at every mesh point to give the profiles. Convective zones are located, \
         and global quantities and consistency checks are computed.",
    );
    bullet(
        ui,
        "Every view is redrawn from that model; the ZAMS line in the HR diagram is recomputed if the composition changed.",
    );
    bullet(
        ui,
        "With Evolve, a background thread advances the star in time step by step (see \"Evolution in time\"), and \
         computes the fundamental radial pulsation mode of every model (see \"Pulsations & oscillations\").",
    );
    muted(
        ui,
        dark,
        "The sections below describe each of these steps in detail, with all formulas and numerical settings. \
         Right-click any equation to copy its LaTeX source.",
    );
}

fn usage(ui: &mut Ui, dark: bool) {
    sub(ui, "Computing a model");
    bullet(
        ui,
        "Pick a preset or move the sliders in the left panel. With Auto-solve on, a new model is computed 0.15 s \
         after you stop editing; otherwise press Solve. Cancel stops a running solve.",
    );
    bullet(
        ui,
        "The status bar at the bottom shows progress (the stage, the Newton iteration and the size of its \
         correction). When the solve is done, it shows the number of mesh points, the iterations and the time taken.",
    );
    bullet(
        ui,
        "If the left panel says \"Parameters changed; showing the previous model\", the views do not yet show the \
         current slider values: a solve is pending or failed. Errors are shown in red in the status bar.",
    );
    bullet(ui, "A note appears below 0.4 M☉, where the ideal-gas equation of state is no longer adequate.");

    sub(ui, "Evolving a star");
    bullet(
        ui,
        "In the Evolution section of the left panel, choose where to stop (end of main sequence, an age, or as far as possible) and press Evolve.",
    );
    bullet(
        ui,
        "The run is computed in the background. The HR track grows and the views follow the latest model (\"follow latest\").",
    );
    bullet(
        ui,
        "Use the timeline above the views to pick any computed age, step with < and >, or Play through the run.",
    );
    bullet(
        ui,
        "Stop interrupts a run, Continue resumes it, Reset to ZAMS discards it. \"Sun's age\" evolves to 4.57 Gyr.",
    );
    bullet(
        ui,
        "Scale (Evolution section of the left panel, also in the timeline; can be set before a run starts): \
         \"fixed for run\" (default) keeps every axis and colour scale constant over the whole run, \
         so only the star changes while you scrub or play. \"adapt to model\" fits each view to the model on screen \
         for the most detail, and zooms the HR diagram in on the track.",
    );
    bullet(
        ui,
        "True size (left panel for both views, or each view's own checkbox): draws the star to scale relative to its \
         largest radius in the run, with a \
         faint outline of that size, so you see it expand. After a giant phase the early star is small; zoom in on the 3D view.",
    );
    bullet(
        ui,
        "For a pulsating star, evolve 1.5–2 M☉ (δ Scuti on the main sequence) or 5 M☉ \"as far as possible\" (Cepheid crossing), then select a model inside the shaded instability strip and open the 3D model.",
    );

    sub(ui, "Comparing and exporting");
    bullet(
        ui,
        "Pin current model keeps a copy for comparison. You can pin up to 7 models, and each keeps its colour. \
         Pinned models are overlaid in Profiles and Compare and plotted in the HR diagram. In Compare, untick a \
         model to hide it or press × to unpin it.",
    );
    bullet(
        ui,
        "Export: Profile CSV (one row per mesh point, see \"Profile quantities & export\"), Summary JSON (parameters, \
         global results, convective zones) and Screenshot PNG (the whole window as currently shown).",
    );
    bullet(ui, "The buttons in the top-right corner switch between light, dark and system theme.");

    sub(ui, "Tabs");
    term(
        ui,
        "Profiles (P, T, ρ):",
        "pressure, temperature and density against radius r/R or mass m/M. The three plots share the x axis and the hover cursor.",
    );
    term(
        ui,
        "More profiles:",
        "luminosity and enclosed mass, opacity, nuclear energy generation by process, temperature gradients, \
         chemical composition, gravitational energy release, convective velocity, convective flux fraction, \
         gas-pressure fraction, the propagation diagram (N, S₁, S₂) with the asymptotic spacings, and the sound speed. \
         The x axis can also be log P, which spreads out the thin outer envelope.",
    );
    term(
        ui,
        "Cross-section:",
        "a slice through the star, with temperature on the left half and density on the right. Convective zones carry animated cells.",
    );
    term(
        ui,
        "3D model:",
        "a cut-away sphere with animated radiative and convective energy transport. It can be coloured by hydrogen \
         content to show the helium core, and pulsates when the star is in the instability strip.",
    );
    term(
        ui,
        "HR diagram:",
        "the model on a computed ZAMS line, with the evolutionary track of the current run, the instability strip, \
         the present-day Sun and the pinned models.",
    );
    term(
        ui,
        "Time evolution:",
        "after an evolution run, linked plots of how the global parameters change with age (or model number): L, R, \
         T_eff, log g, central temperature and density, central composition, core masses, energy sources, pulsation \
         period and solar-like oscillation frequencies. Click a plot to show that age everywhere.",
    );
    term(
        ui,
        "Compare:",
        "a table of global properties for the current and pinned models, plus an overlay plot of any profile quantity.",
    );

    sub(ui, "Plot controls (all 2D plots)");
    bullet(
        ui,
        "Drag to pan. ⌘ + scroll (Ctrl + scroll on Linux/Windows) or pinch to zoom. Right-drag to zoom into a box. Double-click to reset.",
    );
    bullet(
        ui,
        "Hover for values: the tooltip names the model, gives the x coordinate and the de-logged y value with units.",
    );
    bullet(ui, "Gray bands mark the convective zones of the current model.");

    sub(ui, "3D model controls");
    bullet(
        ui,
        "Drag to rotate. Scroll or pinch to zoom towards the cursor (0.3× to 400×). Right-drag or shift-drag to pan. Double-click to reset.",
    );
    bullet(
        ui,
        "\"Zoom to surface layers\" magnifies the outer envelope at the cut, where thin convection zones are otherwise too small to see.",
    );
    bullet(
        ui,
        "\"Cut faces\" chooses temperature or density colouring. \"Core as sphere\" shows the core as a ball (as in many textbook \
         illustrations) or cuts it like the rest of the star. \"Animate\", \"speed\" and \"auto-rotate\" control the motion.",
    );
    bullet(ui, "Hover over the cut faces, the core ball or the surface to read the local state.");
    muted(ui, dark, "Display settings (axes, colouring, zoom, animation) never change the model itself.");
}

fn parameters(ui: &mut Ui, dark: bool, pr: &StarParams) {
    p(ui, "These parameters define a model. The current value is shown live; ranges are those of the sliders.");
    param(
        ui,
        dark,
        "Mass  M  [M☉]",
        &format!("default 1.0 · range 0.1 – 100 (logarithmic) · current {:.3}", pr.mass),
        "Total mass of the star; it determines almost everything else.\n\
         • Low-mass stars are cool, dim and small, with deep convective envelopes, and are fully convective \
         at the lowest masses.\n\
         • Above about 1.2 M☉ a convective core appears and the CNO cycle takes over from the pp chain.\n\
         • Massive stars are hot and luminous, and partly supported by radiation pressure (β < 1).\n\
         Below about 0.4 M☉ electron degeneracy and non-ideal gas effects matter; they are not modelled, \
         so results there are qualitative only.",
    );
    param(
        ui,
        dark,
        "Hydrogen mass fraction  X",
        &format!("default 0.70 · range 0 – (1 − Z) · current {:.4}", pr.x),
        "Fraction of the mass in hydrogen. More hydrogen lowers the mean molecular weight μ, which makes the \
         star dimmer and cooler at a given mass. X also enters the pp-chain rate (∝ X²), the CNO rate (∝ X) \
         and the electron-scattering opacity 0.2 (1 + X).",
    );
    param(
        ui,
        dark,
        "Metallicity  Z",
        &format!("default 0.02 · range 0.0001 – 0.1 (logarithmic) · current {:.5}", pr.z),
        "Fraction of the mass in all elements heavier than helium. Z controls the bound–free and H⁻ opacities. \
         A higher Z makes the star more opaque: dimmer, cooler, larger, with a deeper convective envelope. \
         Z also sets the abundance of the CNO catalysts.",
    );
    param(
        ui,
        dark,
        "Helium mass fraction  Y",
        &format!("derived: Y = 1 − X − Z · current {:.4}", pr.y()),
        "Not set directly; it is whatever mass is left after hydrogen and metals.\n\
         Helium follows metallicity (checkbox under the sliders; CLI --helium-from-z): galaxies enrich helium together \
         with metals. With the box ticked, Y follows Z as in the MIST models (Choi et al. 2016) and the X slider \
         follows Z (Y_p is the primordial helium abundance):",
    );
    math(ui, "par_helium", "");
    p(
        ui,
        "At Z = 0.0142 this gives the protosolar X = 0.7155 (the AGSS09 preset has 0.7154). It is an input mode only; \
         the physics does not change. Measurements of Y_p are now about 0.004 lower (0.245); with the same protosolar \
         anchor the slope would be 1.76.",
    );
    param(
        ui,
        dark,
        "Mixing-length parameter  α_MLT",
        &format!("default 1.8 · range 0.3 – 3.0 · current {:.2}", pr.alpha_mlt),
        "In mixing-length theory a convective blob travels a distance l = α_MLT · H_P (H_P is the pressure scale height, \
         bounded near the centre) before it dissolves into its surroundings. α_MLT only matters where convection is inefficient: in the \
         superadiabatic layers just below the photosphere of stars with convective envelopes. There, a larger α_MLT \
         means more efficient convection, a smaller radius and a higher T_eff. Deep convection is almost exactly \
         adiabatic whatever α_MLT is. Real solar models calibrate α_MLT ≈ 1.7–2.2 to match the Sun.",
    );
    param(
        ui,
        dark,
        "CNO fraction of Z  f_CNO  (Advanced)",
        &format!("default 0.7 · range 0 – 1 · current {:.2}", pr.f_cno),
        "Fraction of the metals that is carbon, nitrogen and oxygen. The CNO-cycle rate uses X_CNO = f_CNO · Z. \
         It only matters for stars hot enough to burn hydrogen via the CNO cycle.",
    );
    param(
        ui,
        dark,
        "Core overshoot  α_ov  (Advanced)",
        &format!("default 0 (off) · range 0 – 0.5 · current {:.2}", pr.alpha_ov),
        "Step overshoot of a convective core, used by evolution runs only. Before each mixing step the convective \
         core is extended from its Schwarzschild boundary r_cc to r_ov = r_cc + α_ov · min(H_P(r_cc), r_cc), and the \
         whole region is mixed. Mixing only: the temperature gradient in the extension stays radiative, so a ZAMS \
         model does not change. A larger α_ov feeds more hydrogen into the core and lengthens the main sequence of \
         stars with convective cores; it has no effect on stars without one. Calibrations for intermediate-mass \
         stars give α_ov ≈ 0.2, with a dependence on mass; the value is not universal.",
    );
    param(
        ui,
        dark,
        "Mesh points  N  (Advanced)",
        &format!("default 800 · range 200 – 3000 (logarithmic) · current {}", pr.n_points),
        "Number of points of the mass grid (see \"Numerical method\"). More points give a more accurate but slower \
         solution. Doubling from 800 to 1600 changes L, R and T_c by less than 0.5% for the Sun.",
    );
    sub(ui, "Presets");
    term(ui, "Sun (X 0.70, Z 0.02):", "1 M☉, the classic textbook solar composition, α_MLT = 1.8, f_CNO = 0.7.");
    term(
        ui,
        "Sun AGSS09:",
        "the same, with the protosolar composition of Asplund et al. (2009): X = 0.7154, Z = 0.0142.",
    );
    sub(ui, "Evolution settings");
    term(
        ui,
        "Stop at:",
        "end of main sequence (X_c < 10⁻³, default), a given age in Gyr (default 4.57), or as far as possible.",
    );
    p(
        ui,
        "The time-step controls are fixed (see \"Evolution in time\"). The mesh-points setting also applies to evolution: more points resolve the thin hydrogen-burning shell after the main sequence better.",
    );
    sub(ui, "Display settings (do not change the model)");
    term(
        ui,
        "x axis:",
        "r/R (fractional radius), m/M (fractional mass, the solver's own coordinate) or log P (pressure, centre on the left).",
    );
    term(ui, "Overlay quantity (Compare):", "which profile is overlaid for all models.");
    term(ui, "3D: cut faces, core as sphere, animate, speed, auto-rotate, zoom:", "see \"How the views are made\".");
    term(ui, "Cross-section: animate convection cells:", "turns the cell animation on or off.");
}

fn results(ui: &mut Ui, dark: bool) {
    p(ui, "Shown in the left panel, the Compare table and the JSON export.");
    term(ui, "R [R☉]:", "radius of the photosphere (optical depth τ = ⅔), i.e. r at the outermost mesh point.");
    term(ui, "L [L☉]:", "luminosity emerging at the surface.");
    term(ui, "T_eff [K] and log g [cgs]:", "effective temperature and surface gravity:");
    math(ui, "res_teff", "");
    term(
        ui,
        "T_c, ρ_c, P_c:",
        "temperature, density and pressure at the innermost mesh point (m = 10⁻⁶ M), which are effectively the central values.",
    );
    term(
        ui,
        "Convective core [m/M]:",
        "mass fraction of the convective region around the centre, or \"none\" if the centre is radiative. Its radius is shown in the views.",
    );
    term(
        ui,
        "Mixed core [m/M]:",
        "shown when α_ov > 0: the convective core extended by step overshoot, i.e. the region mixed in the next \
         evolution step (see \"Evolution in time\").",
    );
    term(
        ui,
        "Conv. envelope base [r/R]:",
        "inner edge of the outer convection zone, or \"none\" if the envelope is radiative. It is the zone that reaches \
         within 10⁻⁶ of the total mass from the surface.",
    );
    term(
        ui,
        "pp / CNO / 3α [%]:",
        "share of the total nuclear energy production from the pp chain, the CNO cycle and helium burning.",
    );
    term(
        ui,
        "Virial check:",
        "for any star in hydrostatic equilibrium the virial theorem makes this ratio 1, up to the surface pressure \
         term 4πR³P_ph and the central cut at m₀ (both far below 10⁻⁵), so any deviation measures the \
         discretisation error (typically 10⁻⁵ to 10⁻⁴):",
    );
    math(ui, "res_virial", "");
    term(
        ui,
        "Energy check:",
        "the energy generated inside (nuclear plus gravitational) must equal the energy leaving the surface. L₀ is the \
         luminosity of the innermost point. The integral uses the solver's own midpoint rule, so it integrates the same \
         difference equation the solver satisfies and holds by construction for a converged model: it is a \
         discretisation consistency check, not an independent test of energy conservation (the independent \
         total-energy check of an evolution run is described in \"Evolution in time\"):",
    );
    math(ui, "res_energy", "");
    term(ui, "Newton iterations:", "iterations in the last solver step (a continuation may take several steps).");
    term(
        ui,
        "Max residual:",
        "largest residual of the discretised equations (continuity, hydrostatic, energy, transport, boundary \
         conditions) of the converged model. They are dimensionless and typically about 10⁻¹⁴; hover over the checks \
         line to see each one.",
    );
    muted(
        ui,
        dark,
        "The virial, energy and residual checks measure numerical consistency. They say nothing about how close the \
         model is to a real star (see \"Accuracy & limitations\").",
    );
    term(ui, "Age:", "time since the ZAMS (0 unless the star has been evolved).");
    term(ui, "Phase:", "evolutionary phase, see \"Evolution in time\".");
    term(
        ui,
        "X_c / He core [m/M]:",
        "central hydrogen mass fraction, and the mass fraction of the hydrogen-exhausted core (X < 10⁻³).",
    );
    term(
        ui,
        "Fund. radial mode:",
        "period of the fundamental radial pulsation mode, see \"Pulsations & oscillations\".",
    );
    term(
        ui,
        "ω/ω_ac (fund. mode):",
        "frequency of that mode over the acoustic cut-off at the photosphere. The pulsation solver's surface condition \
         needs it to be small; a warning appears above 0.5 (see \"Pulsations & oscillations\").",
    );
    term(
        ui,
        "ν_max / Δν [μHz]:",
        "solar-like oscillation frequencies from the asteroseismic scaling relations (meaningful for stars with convective envelopes).",
    );
    term(
        ui,
        "Δν asymptotic [μHz]:",
        "large frequency separation from the interior, [2 ∫ dr / c_s]⁻¹ (see \"Pulsations & oscillations\").",
    );
    term(
        ui,
        "ΔΠ₁ [s]:",
        "asymptotic dipole g-mode period spacing of the largest candidate cavity. Every cavity is listed in More \
         profiles and in the JSON (g_cavities).",
    );
    term(
        ui,
        "Grav. redshift [m/s]:",
        "light leaving the surface loses energy climbing out of the star's gravity. Expressed as the Doppler velocity \
         that would give the same shift:",
    );
    math(ui, "res_redshift", "");
    p(
        ui,
        "It is an output only: general relativity is not modelled in the structure (see \"Accuracy & limitations\").",
    );
    term(
        ui,
        "Screening exponent:",
        "free electrons partly shield the charges of reacting nuclei, which speeds the reactions up by a factor \
         f = exp(…) (see \"Microphysics\"). The formula used holds while the exponent is small (weak screening). The \
         value shown is the largest exponent over the reactions that supply at least 1% of the nuclear energy, where \
         they burn. Above 0.3 a warning says the rates are approximate; at the cap of 2 that they are unreliable. The \
         ZAMS Sun has 0.23. Evolution runs also report the largest value of the run.",
    );
    term(
        ui,
        "Run: total-energy residual, fuel accounting:",
        "evolved models only: the energy checks of the run up to this model, see \"Evolution in time\".",
    );
    sub(ui, "\"Sun today\" reference column");
    p(
        ui,
        "Present-day solar values, for orientation:\n\
         • IAU 2015 nominal values: R = 1 R☉, L = 1 L☉, T_eff = 5772 K, log g = 4.438.\n\
         • Standard solar model centre: T_c ≈ 1.57 × 10⁷ K, ρ_c ≈ 150 g cm⁻³, P_c ≈ 2.4 × 10¹⁷ dyn cm⁻².\n\
         • Base of the convective envelope from helioseismology: 0.713 R.\n\
         • Energy split: about 99% pp, 1% CNO.\n\
         The present Sun has burnt about half of its central hydrogen, so it is brighter, larger and more centrally \
         condensed than a ZAMS model of the same mass.",
    );
    muted(ui, dark, "In the CLI the same quantities are printed by stellar-cli and written by --json.");
}

fn profiles(ui: &mut Ui, dark: bool) {
    p(
        ui,
        "Every model provides these quantities at every mesh point, from the centre to the photosphere. \
         The Profile CSV export writes one column per quantity, in this order:",
    );
    let desc: &[(&str, &str)] = &[
        ("q", "mass fraction m/M"),
        ("m_msun", "enclosed mass m [M☉]"),
        ("r_rsun", "radius r [R☉]"),
        ("r_cm", "radius r [cm]"),
        ("P_dyn_cm2", "total pressure P = P_gas + P_rad [dyn cm⁻²]"),
        ("T_K", "temperature T [K]"),
        ("rho_g_cm3", "density ρ [g cm⁻³]"),
        ("L_lsun", "luminosity L(r) flowing through the sphere of radius r [L☉]"),
        ("kappa_cm2_g", "Rosseland mean opacity κ [cm² g⁻¹]"),
        ("eps_pp", "energy generation by the pp chain [erg g⁻¹ s⁻¹]"),
        ("eps_cno", "energy generation by the CNO cycle [erg g⁻¹ s⁻¹]"),
        ("eps_3a", "energy generation by triple-α helium burning [erg g⁻¹ s⁻¹]"),
        ("eps_total", "total nuclear energy generation ε [erg g⁻¹ s⁻¹]"),
        ("eps_grav", "gravitational energy release ε_grav = −T ds/dt [erg g⁻¹ s⁻¹] (0 on the ZAMS)"),
        ("nabla_rad", "radiative gradient ∇_rad: the d ln T / d ln P needed to carry L by radiation alone"),
        ("nabla_ad", "adiabatic gradient ∇_ad = (∂ ln T / ∂ ln P) at constant entropy"),
        ("nabla", "actual gradient ∇ = d ln T / d ln P used in the model"),
        ("convective", "1 if convective by the Schwarzschild criterion (∇_rad > ∇_ad), else 0"),
        ("v_conv_cm_s", "convective velocity from mixing-length theory [cm s⁻¹] (0 in radiative layers)"),
        ("F_conv_frac", "fraction of the energy flux carried by convection"),
        ("H_P_cm", "pressure scale height H_P = P / (ρ g) [cm]"),
        ("mu", "mean molecular weight μ"),
        ("beta", "gas-pressure fraction β = P_gas / P"),
        ("Gamma1", "first adiabatic exponent Γ₁ (used for pulsations)"),
        ("X", "hydrogen mass fraction (uniform in ZAMS models, depleted by burning during evolution)"),
        ("Y", "helium mass fraction"),
        ("Z", "metal mass fraction"),
        ("c_s_cm_s", "adiabatic sound speed c_s = √(Γ₁ P / ρ) [cm s⁻¹]"),
        ("N2_s-2", "buoyancy frequency squared N² = (g²ρ/P)[δ(∇_ad − ∇) + ∇_μ] [s⁻²]; < 0 where Ledoux-unstable"),
    ];
    egui::Grid::new("csv_cols").num_columns(2).striped(true).spacing([16.0, 3.0]).show(ui, |ui| {
        for name in CSV_COLUMNS {
            let d = desc.iter().find(|(n, _)| n == name).map(|(_, d)| *d).unwrap_or("");
            ui.label(RichText::new(*name).monospace());
            ui.label(d);
            ui.end_row();
        }
    });
    sub(ui, "Summary JSON");
    p(
        ui,
        "params: the input parameters (mass, x, z, alpha_mlt, f_cno, n_points, alpha_ov).\n\
         physics: the description of the physics used.\n\
         globals: every quantity of \"Results: global quantities\", plus newton_iterations, max_correction and residuals.\n\
         convective_zones: one entry per zone, with q_inner, q_outer, r_inner and r_outer (r in units of R).",
    );
    muted(
        ui,
        dark,
        "All values are in CGS units unless the column name says otherwise. Numbers are written with 10 significant digits.",
    );
}

fn equations(ui: &mut Ui) {
    p(
        ui,
        "The model solves the classical equations of stellar structure for a spherical star in hydrostatic and thermal \
         equilibrium. The independent variable is the mass m enclosed within radius r (the Lagrangian form), running \
         from the centre (m = 0) to the surface (m = M):",
    );
    math(ui, "struct_dr", "mass conservation");
    math(ui, "struct_dp", "hydrostatic equilibrium");
    math(ui, "struct_dl", "energy conservation");
    math(ui, "struct_dt", "energy transport");
    p(
        ui,
        "The temperature gradient ∇ comes from radiative diffusion or from convection (next sections). The solver works \
         with ln r, ln P, ln T and L, which keeps quantities spanning many orders of magnitude well conditioned. ρ, κ, ε \
         and ∇ are functions of (P, T, composition) given by the microphysics.",
    );
    sub(ui, "Assumptions");
    bullet(
        ui,
        "ZAMS models are static: the gravitational energy term ε_grav = −T ds/dt is zero. During evolution ε_grav is \
         added to ε (see \"Evolution in time\").",
    );
    bullet(
        ui,
        "ZAMS models are chemically homogeneous: X, Y, Z are the same everywhere. Evolved models store a composition \
         per mesh point, set by burning and convective mixing (see \"Evolution in time\").",
    );
    bullet(ui, "No neutrino losses, rotation, magnetic fields, mass loss or element diffusion.");
}

fn microphysics(ui: &mut Ui, dark: bool) {
    sub(ui, "Equation of state: fully ionised ideal gas plus radiation");
    math(ui, "eos_p", "");
    math(ui, "eos_mu", "mean molecular weight (full ionisation)");
    math(ui, "eos_beta", "");
    math(ui, "eos_cp", "");
    math(ui, "eos_nabla_ad", "0.4 for pure gas, 0.25 for pure radiation");
    p(
        ui,
        "Given P and T, the density follows directly from the first equation. During Newton iterations β is limited to \
         at least 10⁻³, which keeps unphysical trial states finite; converged models are rejected if they reach this floor.",
    );

    sub(ui, "Opacity: analytic Rosseland mean (CGS)");
    math(ui, "kap_es", "electron scattering");
    math(ui, "kap_ff", "Kramers free–free");
    math(ui, "kap_bf", "Kramers bound–free, t = guillotine factor");
    math(ui, "kap_hm", "H⁻ ion (cool envelopes)");
    math(ui, "kap_mol", "molecular floor");
    math(ui, "kap_total", "combined opacity");
    p(
        ui,
        "The harmonic combination makes H⁻ dominate only in cool layers. Without it the steep T⁹ term would run away at \
         high temperature. For Z = 0, H⁻ uses Z = 10⁻⁴. These coefficients follow Carroll & Ostlie (2017), converted to CGS.",
    );
    p(
        ui,
        "Across regimes (plot in docs/about/07-microphysics.md):\n\
         • The H⁻/Kramers blend is smooth: where the two are equal (log T ≈ 4.18 at log R = −1.5) it gives half of each, \
         a dip of at most a factor of 2 with no discontinuity.\n\
         • Too high in the core: κ ≈ 2.7 at the centre of the 1 M☉ ZAMS model, against roughly 1.2–1.5 in standard solar \
         models with OPAL/OP opacities. Since L ∝ 1/κ in a radiative core (a homology argument), this is the leading \
         suspect for the underluminous analytic Sun; its share has not been isolated.\n\
         • Too low below the convection zone: about 9 at the present Sun's convection-zone base (2.2 × 10⁶ K, \
         0.19 g cm⁻³) against roughly 20, so the convective envelope ends too shallow.\n\
         • No line opacity or ionisation bumps (such as the Z bump near log T ≈ 5.3): Kramers laws are smooth power laws.",
    );

    sub(ui, "Nuclear energy generation [erg g⁻¹ s⁻¹]");
    p(
        ui,
        "The pp and CNO rates are built from the reaction rates of Solar Fusion III (Acharya et al. 2025). The rate of \
         each reaction is the thermal average of its cross section, written with the astrophysical factor S(E) and the \
         values S(0), S′(0), S″(0) recommended in their Table I (μ is the reduced mass, b sets the Coulomb barrier):",
    );
    math(ui, "nuc_sigmav", "thermal reaction rate");
    p(
        ui,
        "S(0) in MeV b: p + p 4.09 × 10⁻²⁵, ³He + ³He 5.21, ³He + ⁴He 5.61 × 10⁻⁴, p + ¹⁴N 1.68 × 10⁻³ (for the last \
         only S(0) is recommended). The integrals are computed once, tabulated against temperature and interpolated \
         smoothly. Compared with the formulas of Kippenhahn & Weigert (1990) used before, the p + p rate is 5–7% higher \
         and the p + ¹⁴N rate about half.",
    );
    math(ui, "nuc_t6", "");
    math(ui, "nuc_pp", "pp chain");
    p(
        ui,
        "ε_pp,I counts the pp I chain only: two p + p reactions per helium nucleus, 26.20 MeV deposited. Once helium is \
         abundant, ³He also burns with ⁴He (the pp II chain), which makes more helium per p + p reaction; ψ accounts for \
         it (pp III is neglected):",
    );
    math(ui, "nuc_branch", "pp II branching");
    p(
        ui,
        "• φ (1 to 2) is the number of helium nuclei made per pp I chain, F_II the share made through pp II, which \
         deposits 2.1% less energy per helium nucleus.\n\
         • α uses the same Solar Fusion III rates, with the screening factors.\n\
         • ψ is near 1 on the ZAMS and grows as helium builds up: φ = 1.02 at the ZAMS Sun's centre, 1.10 at 4.57 Gyr, where 17% of the helium is made through pp II.\n\
         • The hydrogen consumed per unit energy changes with it; see \"Evolution in time\".",
    );
    math(ui, "nuc_cno", "CNO cycle: the slowest step, p + ¹⁴N, sets the rate");
    math(ui, "nuc_c", "");
    math(ui, "nuc_3a", "triple-α (Kippenhahn & Weigert 1990)");
    p(ui, "Weak electron screening (Salpeter / Clayton):");
    math(ui, "nuc_screen", "");
    p(
        ui,
        "The charge products are Z₁Z₂ = 1 (p + p), 7 (p + ¹⁴N, the slowest CNO step), and 4 and 8 (α + α, α + ⁸Be). \
         The exponent is capped at 2, since weak screening is only valid for small values. Every model reports its \
         largest exponent, with a warning above 0.3 and at the cap (see \"Results\").",
    );
    muted(
        ui,
        dark,
        "EOS, opacity and nuclear rates sit behind interfaces (Eos, Opacity, Nuclear), so tabulated physics can replace them without touching the solver.",
    );
}

fn transport(ui: &mut Ui, dark: bool) {
    sub(ui, "Radiative transport");
    p(ui, "If all the energy is carried by photon diffusion, the temperature gradient must be");
    math(ui, "tr_nabla_rad", "radiative gradient");
    sub(ui, "Convective instability: Schwarzschild criterion");
    math(ui, "tr_schwarzschild", "");
    p(
        ui,
        "In a convective layer a rising blob that expands adiabatically stays hotter, and so lighter, than its surroundings, \
         and keeps rising. The Schwarzschild criterion is used throughout, also in evolved models that have composition \
         gradients: their stabilising effect (Ledoux criterion, semiconvection) is neglected.",
    );
    sub(ui, "Convective transport: mixing-length theory (MLT)");
    p(ui, "Following Böhm-Vitense (1958) in the formulation of Kippenhahn, Weigert & Weiss (2012, §7):");
    math(ui, "mlt_defs", "gravity, pressure scale height, mixing length");
    p(
        ui,
        "H_P grows without bound towards the centre, where g → 0. The mixing length therefore uses the smaller of H_P \
         and √(P/Gρ²), which stays finite (Paxton et al. 2011). This changes only the convective velocity near the \
         centre; the gradient there is adiabatic either way, and no global result moves.",
    );
    math(ui, "mlt_u", "");
    math(ui, "mlt_cubic", "the MLT cubic");
    p(
        ui,
        "The cubic has exactly one root with ξ ≥ U. It is solved for η = ξ − U by Newton's method from the upper bound \
         √(U² + W) − U, which converges monotonically to machine precision. Then:",
    );
    math(ui, "mlt_nabla", "actual gradient");
    math(ui, "mlt_v", "convective velocity");
    math(ui, "mlt_fconv", "convective flux fraction");
    bullet(ui, "Deep interior: convection is efficient (U ≪ 1) and ∇ ≈ ∇_ad. The layer is almost exactly adiabatic.");
    bullet(
        ui,
        "Just below the photosphere: convection is inefficient (U ≫ 1) and ∇ rises towards ∇_rad. This superadiabatic layer \
         is the only place where α_MLT matters.",
    );
    p(
        ui,
        "Overshoot: with α_ov > 0 a convective core is extended for the composition mixing only (see \"Evolution in \
         time\"); the temperature gradient in the extension stays radiative, so the structure equations are unchanged.",
    );
    muted(ui, dark, "Not included: overshoot of envelopes, semiconvection, turbulent pressure.");
}

fn boundaries(ui: &mut Ui, dark: bool) {
    sub(ui, "Centre");
    p(
        ui,
        "The equations are singular at m = 0, so the first mesh point sits at m₀ = 10⁻⁶ M and uses the leading terms of \
         the series expansion about the centre:",
    );
    math(ui, "bc_centre", "");
    sub(ui, "Surface: grey Eddington atmosphere");
    math(ui, "bc_eddington", "");
    p(
        ui,
        "The outermost mesh point is the photosphere, at optical depth τ = ⅔, where T = T_eff. Two conditions hold there:",
    );
    math(ui, "bc_lum", "");
    math(ui, "bc_pressure", "");
    p(
        ui,
        "The pressure condition comes from integrating hydrostatic equilibrium in optical depth, dP/dτ = g/κ, from τ = 0, \
         where only the radiation pressure a T_eff⁴ / 6 remains. Here κ is evaluated self-consistently at the photospheric \
         (P, T).",
    );
    muted(
        ui,
        dark,
        "2 central + 2 surface conditions + 4(N − 1) difference equations = 4N equations for the 4N unknowns.",
    );
}

fn numerics(ui: &mut Ui, dark: bool) {
    sub(ui, "Mesh");
    p(
        ui,
        "N points in q = m/M, from 10⁻⁶ to 1, uniformly spaced in the stretched coordinate below. Its three terms place \
         about 15% / 45% / 40% of the points in the centre / interior / envelope:",
    );
    math(ui, "num_mesh", "");
    p(
        ui,
        "This gives logarithmic spacing in q near the centre, uniform spacing through the interior, and logarithmic \
         spacing in 1 − q towards the surface, down to about 10⁻¹³. The outer layers contain very little mass, but the \
         pressure there falls by many orders of magnitude, so they need many points.",
    );
    sub(ui, "Unknowns and discretisation");
    p(
        ui,
        "Each mesh point carries four unknowns. L_scale is the expected luminosity, which keeps all unknowns of order 1:",
    );
    math(ui, "num_unknowns", "");
    p(
        ui,
        "Between neighbouring points, the change of each unknown equals Δm times its derivative f (the right-hand sides of \
         the structure equations) at the midpoint state. This scheme is second-order accurate:",
    );
    math(ui, "num_difference", "");
    sub(ui, "Henyey relaxation (Newton–Raphson on the whole star)");
    p(
        ui,
        "All 4N residuals F(y) (difference equations plus boundary conditions) are driven to zero together by Newton's \
         method (Henyey, Forbes & Gould 1964). This is the standard method of stellar evolution codes, and it also allows \
         time evolution to be added later. Each iteration solves a linear system and applies a damped correction:",
    );
    math(ui, "num_newton", "");
    bullet(
        ui,
        "The damping limits in λ (one per unknown) are 0.2 for ln r, 0.5 for ln P, 0.2 for ln T and 0.5 for L/L_scale.",
    );
    bullet(
        ui,
        "Jacobian: built by central finite differences, one mesh point at a time. Each equation involves only \
         neighbouring points, so no analytic derivatives of the physics are needed:",
    );
    math(ui, "num_jacobian", "for ln r, ln P, ln T; relative step for L");
    bullet(
        ui,
        "Linear solve: J is banded (5 sub- and 5 super-diagonals). It is solved by Gaussian elimination with partial \
         pivoting in banded storage, so the cost grows linearly with N.",
    );
    bullet(
        ui,
        "Convergence: the largest correction is below 10⁻⁸ with an undamped step (λ = 1). This typically takes 10–30 \
         iterations, with quadratic convergence at the end. The limit is 250 iterations from a cold start, 60 from a \
         warm start and 40 for each evolution step. A small correction alone does not prove the equations are satisfied, so the residuals of the \
         converged solution are also evaluated; converged models reach about 10⁻¹⁴.",
    );
    bullet(
        ui,
        "Every converged model is checked: all values finite, r increasing outward, P decreasing outward, L > 0, and β \
         above the floor, and every residual below 10⁻⁶. A model failing any check is rejected.",
    );
    sub(ui, "Post-processing");
    bullet(ui, "Profiles: the full local physics (EOS, κ, ε, MLT) is evaluated at every mesh point.");
    bullet(
        ui,
        "Zone boundaries: located where ∇_rad − ∇_ad changes sign, by linear interpolation between mesh points.",
    );
    bullet(ui, "Integrals (energy split, virial and energy checks) use the trapezoidal rule in m.");
    sub(ui, "Validation");
    p(
        ui,
        "The same generic solver reproduces the Lane–Emden polytrope zeros ξ₁ = π (n = 1) and 6.89685 (n = 3). \
         Automated tests also cover:\n\
         • EOS limits, MLT limits, and positivity of opacity and rates;\n\
         • the banded solver against a dense solver;\n\
         • mesh refinement, and warm start against cold start;\n\
         • qualitative trends along the ZAMS (convective cores, CNO dominance, monotonic L(M)).",
    );
    muted(
        ui,
        dark,
        "A solar model takes about 13 iterations and 40–60 ms from a cold start; after a slider change the warm start typically takes 20–100 ms.",
    );
}

fn initial(ui: &mut Ui, dark: bool) {
    sub(ui, "Cold start: an n = 3 polytrope");
    p(
        ui,
        "The Lane–Emden equation for n = 3 (the Eddington standard model) is integrated with 4th-order Runge–Kutta, step 10⁻³:",
    );
    math(ui, "init_lane_emden", "");
    p(ui, "It is scaled to homology estimates of the ZAMS radius and luminosity (M in solar masses):");
    math(ui, "init_homology", "");
    p(ui, "Central values and the polytropic structure; T comes from inverting the EOS at each point:");
    math(ui, "init_central", "");
    p(
        ui,
        "The outer layers are patched with a plane-parallel hydrostatic envelope (the larger of the two pressures is used):",
    );
    math(ui, "init_envelope", "");
    bullet(ui, "L(m): ε is integrated over the guessed structure and rescaled to the estimated luminosity.");
    sub(ui, "Warm start and continuation");
    p(
        ui,
        "When parameters change, the solver starts from the model already on screen, scaled by homology. It is \
         interpolated onto the new mesh if N changed:",
    );
    math(ui, "init_scaling", "");
    bullet(
        ui,
        "It first tries to reach the new parameters in a single step. If Newton fails, it halves the step along the path \
         (geometric in M, linear in X, Z, α_MLT, f_CNO), down to 1/1024 of the path. It doubles the step again after each success.",
    );
    bullet(
        ui,
        "Fallbacks: a cold start from the polytrope; then a 1 M☉ reference model (same composition, then the solar \
         composition), followed by continuation to the target.",
    );
    sub(ui, "Background solving");
    p(
        ui,
        "Solves run in a background thread, so the interface stays responsive. A new request cancels the one in progress, \
         and the results of superseded requests are discarded. Edits are debounced by 0.15 s. The HR-diagram ZAMS line is \
         computed in its own thread.",
    );
    muted(
        ui,
        dark,
        "The solver exposes single Newton steps (HenyeyNewton::step), so a future web build can iterate a few steps per frame.",
    );
}

fn evolution(ui: &mut Ui, dark: bool) {
    p(
        ui,
        "Starting from the ZAMS model, Evolve follows the star in time on the same mass mesh (the total mass is \
         constant: no mass loss). Each time step has two parts. First the composition changes by nuclear burning \
         and convective mixing. Then the structure is solved again, including the energy released or absorbed as \
         layers contract or expand.",
    );
    sub(ui, "Composition");
    math(ui, "evo_burn_h", "hydrogen burning (pp chain and CNO cycle)");
    p(
        ui,
        "Each channel uses its own energy per gram of hydrogen burnt: 26.20 MeV per helium nucleus for the pp I chain \
         and 25.03 MeV for the CNO cycle, after neutrino losses. The pp II branch deposits 2.1% less than pp I, so \
         Q_pp^eff falls with the pp II share F_II (see \"Microphysics\").",
    );
    math(ui, "evo_burn_he", "helium burning (triple-α); its carbon counts as metals");
    math(ui, "evo_mix", "instant mixing");
    p(
        ui,
        "Step overshoot (input α_ov, default 0 = off): before mixing, a convective core that starts at the centre is \
         extended from its Schwarzschild boundary r_cc to",
    );
    math(ui, "evo_overshoot", "every mesh point with r ≤ r_ov is mixed with the core");
    p(
        ui,
        "H_P(r_cc) is interpolated between the mesh points around the boundary; min(H_P, r_cc) keeps the extension finite \
         for small cores. The overshoot is applied to the mixing only: the temperature gradient in the extension stays \
         radiative. The mixed core feeds more hydrogen into the burning region. With α_ov = 0.2 the main sequence of \
         3 M☉ lasts 338 Myr instead of 269 Myr (+26%), and of 5 M☉ 88.6 Myr instead of 72.0 Myr (+23%). A 1.1 M☉ star, \
         which has no convective core, is unaffected. α_ov ≈ 0.2 is a typical calibration for intermediate-mass stars; \
         it depends on mass and is not universal.",
    );
    p(ui, "After mixing, every mixed zone burns as one unit, at its mass-averaged rate:");
    math(ui, "evo_zone_burn", "the same for helium");
    p(
        ui,
        "Burning uses the energy-generation rates of the previous model (an explicit step, kept small by the time-step \
         limits below), and no zone burns more than it contains. The CNO catalyst abundance X_CNO = f_CNO · Z keeps its \
         initial value. Mixing uses the Schwarzschild convective zones of the previous model, so a convective core \
         that shrinks during the main sequence leaves a hydrogen gradient behind.",
    );
    sub(ui, "Structure with gravitational energy");
    p(ui, "The energy equation gains the gravitational term:");
    math(ui, "evo_dl", "");
    math(ui, "evo_egrav", "");
    p(
        ui,
        "ε_grav is evaluated implicitly, against the previous model at the same mass. The Henyey solver therefore \
         handles contraction (ε_grav > 0: the layer heats up and releases gravitational energy) and expansion \
         (ε_grav < 0), which happen on the thermal timescale after the main sequence. The previous model is the \
         starting guess for Newton, and the luminosity unknown is re-scaled to the current luminosity every step.",
    );
    sub(ui, "Time steps");
    math(ui, "evo_dt", "");
    bullet(
        ui,
        "Burning limit: no shell may lose more than 0.01 of X, or more than 30% of its hydrogen, in one step; the \
         same holds for helium wherever triple-α burns. A mixed zone counts as a single shell. The 30% rule has a \
         floor of 10⁻⁵, and shells with X (or Y) ≤ 10⁻⁸ are skipped, so the step does not shrink without bound as \
         the fuel runs out.",
    );
    bullet(
        ui,
        "Structural limit: the step grows or shrinks so that ln r, ln P, ln T and ln L change by about 0.02 per step. \
         A step that changes them by more than 0.08 is rejected and repeated with half the time step.",
    );
    bullet(
        ui,
        "If Newton does not converge within 40 iterations, the step is divided by 3 and retried, up to 12 times; the \
         smallest step is 1 year.",
    );
    bullet(ui, "First step: 10% of the burning limit, at most 10⁶ years.");
    sub(ui, "Evolutionary phases");
    term(ui, "Zero-age main sequence:", "age 0, homogeneous composition.");
    term(ui, "Main sequence:", "core hydrogen burning, central X ≥ 10⁻³.");
    term(
        ui,
        "H-shell burning (subgiant / Hertzsprung gap):",
        "central hydrogen exhausted (X_c < 10⁻³); hydrogen burns in a shell around the helium core.",
    );
    term(ui, "Red giant:", "shell burning, radius above 5 × the ZAMS radius and T_eff below 5000 K.");
    term(ui, "Core helium burning:", "triple-α supplies more than 10% of the nuclear energy.");
    sub(ui, "Stop conditions");
    term(ui, "End of main sequence:", "central hydrogen below 10⁻³.");
    term(ui, "Age:", "the chosen age; the last step is shortened to land on it.");
    term(
        ui,
        "As far as possible:",
        "until T_eff drops below 2600 K (outside the range of the analytic opacities), T_c exceeds 3 × 10⁸ K, \
         or the solver cannot continue even with tiny steps. Results after the main sequence are qualitative only \
         (see \"Accuracy & limitations\").",
    );
    term(ui, "Safety limits (always):", "3000 steps and an age of 30 Gyr.");
    sub(ui, "Checks");
    p(
        ui,
        "Total-energy check. Over a run, the energy radiated must equal the nuclear energy released minus the change of \
         the star's thermal and gravitational energy (u is the thermal energy per gram; rest mass and nuclear binding \
         energy enter only through ε_nuc):",
    );
    math(ui, "evo_energy", "");
    p(
        ui,
        "The residual has two parts, reported separately. The ε_grav of the solver leaves out the heat from the change of \
         composition at fixed P and T; the run adds that term up (the change of the enthalpy u + P/ρ), and it does not \
         shrink with the time step. The rest is numerical and shrinks with the time step. For the Sun evolved to \
         4.57 Gyr the residual is −4.5 × 10⁻⁴, of which −4.45 × 10⁻⁴ is the composition term and −9 × 10⁻⁷ numerical \
         (halved by halving the time steps); for 5 M☉ to the end of the main sequence −1.09 × 10⁻³ and 2 × 10⁻⁶. \
         Automated tests require below 1% and check both parts on halving the steps.",
    );
    p(
        ui,
        "Fuel accounting. The hydrogen burnt, each step with the effective energy per gram of the model it burnt in, \
         against the energy radiated:",
    );
    math(ui, "evo_budget", "");
    p(
        ui,
        "It is close to 1 only on the main sequence, where hydrogen burning supplies nearly all the luminosity. For the \
         Sun evolved to 4.57 Gyr it is 0.990 (0.995 with halved time steps: the difference is the first-order error of \
         the explicit burning step); for 3 M☉ to the end of the main sequence 0.996. Along every run the virial check \
         stays at 1 to about 10⁻⁴, and the energy check, which includes ε_grav, at 1. The app shows the total-energy \
         residual and the fuel accounting under the checks of an evolved model; the CLI prints them.",
    );
    sub(ui, "Time-step convergence");
    p(
        ui,
        "The same tracks were run with the step targets halved repeatedly (structural change 2%, 1%, 0.5%, 0.25%, with \
         the hydrogen limit ΔX halved alongside). The results converge at first order, as expected from burning \
         explicitly with the previous model's rates: each halving roughly halves the difference. At the default \
         settings the errors are small compared with the physics uncertainties:\n\
         • 1 M☉ at 4.57 Gyr: L, R and T_c within about 0.1%.\n\
         • 3 M☉: main-sequence lifetime within about 0.6% (268.6 Myr against about 266.9 Myr extrapolated); \
         convective-core mass at X_c = 0.3 within 0.1%.\n\
         • Values exactly at the end of the main sequence scatter by up to about 1%, because the last step lands \
         anywhere below X_c = 10⁻³.\n\
         A step is rejected when it changes the structure by more than 4 × the target; in these runs the step-size \
         control keeps every step well below that, so a stricter rejection threshold (1.5 ×) changes nothing. \
         Run the study with: cargo run --release -p stellar-core --example timestep_study",
    );
    sub(ui, "In the app");
    bullet(
        ui,
        "Evolution section (left panel): choose the stop condition and press Evolve. Stop interrupts a run; Continue \
         resumes it from its latest model, for example with a later stop condition. Reset to ZAMS discards the run. \
         Sun's age evolves to 4.57 Gyr.",
    );
    bullet(
        ui,
        "Timeline (above the views): drag to choose an age, step with < and >, Play through the run, or tick \"follow \
         latest\" to watch it as it is computed. Every tab shows the selected model.",
    );
    bullet(
        ui,
        "Every step is added to the HR track. Full models are kept for the timeline when T_eff changes by more than \
         0.006 dex, L or R by more than 0.015 dex, X_c by more than 0.02, or when the phase or pulsation state changes. \
         This keeps memory bounded on runs of thousands of steps.",
    );
    bullet(
        ui,
        "Pin current model pins the selected age. Track CSV exports the global properties of every timeline model; \
         the Time evolution tab exports every step (history CSV).",
    );
    bullet(
        ui,
        "Time evolution tab: twelve linked time-series panels covering the whole run. The x axis is age, or model number, \
         which spreads out the fast phases after the main sequence. Vertical lines mark phase changes, a vertical \
         marker shows the selected age, shaded intervals in the period panel mark the instability strip, and a click \
         on any panel selects that age. \"fit all on screen\" shrinks the panels to show all of them at once.",
    );
    bullet(ui, "Changing any input parameter discards the run, because the ZAMS star itself changes.");
    bullet(
        ui,
        "Scale: with \"fixed for run\" the y axes, the log P axis and the colour scales span every model computed in the \
         run (updated as models arrive). Fractions use 0–1. The energy-generation range starts 8 dex below each model's \
         peak, and the HR diagram shows the whole diagram. With \"adapt to model\" each view fits the model on screen, and \
         the HR diagram is zoomed on the track. Switching the mode resets any manual zoom.",
    );
    muted(
        ui,
        dark,
        "Typical runs: the Sun to 4.57 Gyr in 27 steps (0.3 s). A 5 M☉ star from the ZAMS across the Hertzsprung gap \
         to 2600 K takes about 1650 steps (15 s), computed in the background while the timeline fills.",
    );
}

fn pulsation(ui: &mut Ui, dark: bool) {
    sub(ui, "Radial pulsation modes");
    p(
        ui,
        "For every model (the ZAMS model and every model of an evolution run), the fundamental radial mode is computed \
         from the linear adiabatic wave equation (ξ: relative radial displacement, p: relative pressure perturbation):",
    );
    math(ui, "puls_lawe", "");
    math(ui, "puls_bc", "");
    math(ui, "puls_gamma1", "ideal gas + radiation");
    p(
        ui,
        "The dimensionless eigenvalue Ω² is an extra unknown of the boundary-value problem. It is solved with the same \
         Henyey method as the structure, with three unknowns per mesh point (ξ, p, Ω²). The fundamental mode is the \
         solution without nodes in ξ. During evolution the previous model's mode is the starting guess, which keeps \
         the solver on the fundamental branch. At the midpoints between mesh points r, P and ρ are geometric means, \
         as in the structure solver. The surface condition treats the photosphere as a perfect reflector. That holds \
         when the mode frequency is well below the acoustic cut-off of the atmosphere; the ratio is reported with \
         every mode, and a warning appears above 0.5:",
    );
    math(ui, "puls_cutoff", "");
    p(ui, "For the ZAMS Sun ω/ω_ac = 0.095. Period and pulsation constant:");
    math(ui, "puls_period", "");
    p(ui, "Validation: for a uniform-density sphere the solver reproduces the exact result");
    math(ui, "puls_homogeneous", "");
    p(ui, "and its overtones (Pekeris 1938), which the library computes with a given number of nodes:");
    math(ui, "puls_overtones", "");
    p(
        ui,
        "to 10⁻³ for k = 0–3. Polytropes of index 1.5 and 3 give fundamental Ω² = 2.706 and 9.254, in agreement \
         to 10⁻⁴ with an independent calculation (the same equation solved by shooting in Lane–Emden variables). \
         Doubling the mesh from 800 to 1600 points changes the solar period by 4 × 10⁻⁵. The ZAMS Sun's fundamental \
         radial mode comes out at 1.00 hours. This mode is not observed in the Sun; the standard solar model Model S \
         gives 64.6 minutes (258 μHz). The agreement is partly a coincidence, since the model Sun has R = 1.12 R☉.",
    );
    sub(ui, "What is calculated and what is not");
    term(ui, "Calculated:", "the fundamental adiabatic radial mode: period, Ω², Q and the mode shape ξ(r).");
    term(
        ui,
        "Empirical inference:",
        "membership of the instability strip, from the model's position in the HR diagram.",
    );
    term(ui, "Scaling relations:", "ν_max and Δν, estimated from g, T_eff and the mean density.");
    term(
        ui,
        "Not modelled:",
        "driving and damping (growth rates), the excitation mechanism, amplitudes and non-radial modes. Radial overtones \
         are computed by the library (and tested) but not shown in the app. δ Scuti stars in particular usually \
         pulsate in several radial and non-radial modes at once.",
    );
    sub(ui, "Which stars pulsate: the instability strip");
    p(
        ui,
        "Whether a mode is actually driven (by the κ mechanism in the He II ionisation zone) needs a non-adiabatic \
         analysis, which is not done here. Instead, an approximate classical instability strip in the HR diagram decides \
         it. These are linear fits to the observed strip, from δ Scuti stars (6300–8700 K at log L ≈ 0.5) to classical \
         Cepheids (about 5400–6800 K at log L ≈ 3.5):",
    );
    math(ui, "puls_strip", "");
    term(
        ui,
        "δ Scuti-type:",
        "inside the strip with log L < 1.7, on or near the main sequence (periods of about 0.02–0.3 d).",
    );
    term(
        ui,
        "Cepheid-type:",
        "inside the strip with log L ≥ 1.7, for example while a 4–10 M☉ star crosses the Hertzsprung gap (periods of days).",
    );
    p(
        ui,
        "Inside the strip the 3D model pulsates with the computed fundamental mode. Each layer is displaced by ξ(r), \
         largest at the surface and small in the core. One cycle is shown per 1.6 s, the amplitude is exaggerated (6% \
         of the radius), and the surface is brightest near maximum compression. The HR diagram shows the strip as a \
         shaded band, and the timeline names the pulsator type and period.",
    );
    sub(ui, "Solar-like oscillations");
    p(
        ui,
        "Stars with convective envelopes, like the Sun, oscillate in many p modes excited by turbulent convection. \
         Their frequency of maximum power and their large frequency separation follow the asteroseismic scaling \
         relations, shown in the results panel:",
    );
    math(ui, "seis_scaling", "");
    sub(ui, "Characteristic frequencies and asymptotic spacings");
    p(ui, "Every model also carries the sound speed, the Lamb frequency and the buoyancy (Brunt–Väisälä) frequency N:");
    math(ui, "seis_cs_lamb", "sound speed and Lamb frequency");
    math(ui, "seis_n2", "buoyancy frequency; φ = 1 for the ideal-gas EOS");
    p(
        ui,
        "∇_μ is taken from the profile by differences. N² < 0 marks Ledoux instability. In a homogeneous model it \
         coincides with the convective (Schwarzschild) flag; with a composition gradient a Schwarzschild-unstable layer \
         can have N² > 0. The definition N² = g(Γ₁⁻¹ d ln P/dr − d ln ρ/dr) agrees with this formula to a median of \
         2 × 10⁻⁵ in an evolved 3 M☉ model; it serves as a test only, because it differences nearly equal numbers. \
         p modes propagate where the frequency is above both N and S_l, g modes where it is below both (the propagation \
         diagram in More profiles).",
    );
    math(ui, "seis_asymptotic", "asymptotic large separation and dipole period spacing");
    p(
        ui,
        "Δν from the interior replaces the scaling relation by a value that depends on the structure; for the ZAMS Sun \
         the two agree to 1% (114.8 against 113.7 μHz). ΔΠ₁ is computed separately for every connected N² > 0 region \
         (a candidate g-mode cavity: a radiative core, the radiative zone between a convective core and envelope, a \
         radiative envelope); regions are never summed. The value over the whole region is the low-frequency limit. \
         Both are asymptotic estimates (high radial order), not computed spectra, and Δν carries the model's surface \
         errors.",
    );
    muted(
        ui,
        dark,
        "Real Cepheids change their radius by a few percent and δ Scuti stars by much less; periods range from about \
         half an hour to weeks. The animation is not to scale in time or amplitude.",
    );
}

fn views(ui: &mut Ui, dark: bool) {
    sub(ui, "Colour");
    bullet(
        ui,
        "Lines and points use a colourblind-safe categorical palette, validated for colour-vision deficiencies. The current \
         model is always slot 1 (blue). Each pinned model keeps the slot it received when pinned, so colours never move \
         between models. A legend is always shown; the Compare table lists exact values.",
    );
    bullet(ui, "Gray bands: convective zones of the current model (a neutral colour, never a series colour).");
    bullet(
        ui,
        "During an evolution run, the timeline's Scale switch sets whether axes and colour scales stay fixed for the whole \
         run (default) or adapt to the shown model; see \"Evolution in time\".",
    );
    bullet(
        ui,
        "Continuous fields use single-hue, logarithmic colour scales: orange for temperature, blue for density. In 2D views \
         the low end fades towards the background. In 3D, hotter or denser is always brighter.",
    );
    sub(ui, "Profiles and More profiles");
    bullet(
        ui,
        "Logarithmic quantities are plotted as log₁₀, with ticks labelled 10ⁿ (plus 2·10ⁿ and 5·10ⁿ when the range is narrow).",
    );
    bullet(
        ui,
        "Energy generation is hidden where it falls below 10⁻⁸ of its peak. Convective velocity is drawn only in convective layers.",
    );
    bullet(
        ui,
        "log P axis: centre on the left. It is useful because the outer 10% of the radius spans almost half of the pressure range (46% for the Sun), and the outer 1% still spans about a fifth.",
    );
    bullet(
        ui,
        "Propagation diagram: N/2π where N² > 0, and the Lamb frequencies S₁/2π and S₂/2π, in μHz. Below the plots, \
         the asymptotic Δν and every candidate g-mode cavity with its ΔΠ₁.",
    );
    sub(ui, "Cross-section");
    p(
        ui,
        "A disc split in two halves, coloured by log T (left) and log ρ (right), between the surface and central values. \
         Convective zones are outlined and carry animated rolls whose direction alternates, as in real convection cells. \
         The hover read-out is interpolated from the profile.",
    );
    sub(ui, "3D model");
    bullet(
        ui,
        "The sphere has the octant x, y, z > 0 removed. The three exposed quarter-discs are coloured by T or ρ from the radial profile.",
    );
    bullet(
        ui,
        "The surface shows the approximate colour of a black body at T_eff. It is mottled by an animated granulation pattern when the envelope is convective.",
    );
    bullet(
        ui,
        "Core ball: the convective core, if the star has one that is smaller than the whole star; otherwise the sphere \
         that produces 90% of the luminosity. It \"boils\" when convective.",
    );
    bullet(
        ui,
        "Radiative zones: 180 photon sparks (number per zone ∝ its face area) random-walk outward, with short trails \
         showing the zig-zag path of radiative diffusion.",
    );
    bullet(
        ui,
        "Convective zones: rows of circulating cells about 0.16 R deep, at most 40 per row and face. Parcels rising \
         (moving outward) are amber; sinking parcels are blue.",
    );
    bullet(ui, "White arcs: boundaries between radiative and convective zones.");
    bullet(
        ui,
        "Cut-face colouring by hydrogen X (linear scale): hydrogen-rich gas dark, hydrogen-poor gas light. During \
         evolution the burnt core shows up as a light region (or a light core ball).",
    );
    bullet(
        ui,
        "Pulsation (models in the instability strip): every layer is drawn at r · (1 + A ξ(r) sin φ), with the computed \
         fundamental-mode shape ξ(r) (1 at the surface), A = 6% and one cycle per 1.6 s. The surface brightens near \
         maximum compression. The legend gives the real period.",
    );
    bullet(
        ui,
        "Rendering: orthographic projection on the CPU. About 40,000 shaded triangles are depth-sorted every frame \
         (painter's algorithm, about 2 ms) and drawn as one mesh with 4× multisampling. Overlays use an exact ray test \
         against the cut-away geometry, so they are hidden correctly behind the sphere and the core.",
    );
    muted(
        ui,
        dark,
        "The cross-section and the 3D view are visualisations of a 1D model. Photon paths, convection cells, granulation \
         and boiling are illustrative: their size, shape and speed are not computed, and this is not a 3D \
         radiation-hydrodynamics simulation. Animation speeds are illustrative, not to scale. In the Sun a photon needs of order 10⁴–10⁵ years to diffuse out \
         of the radiative zone, while convective turnover takes from minutes (granules) to about a month (base of the envelope).",
    );
    sub(ui, "HR diagram");
    p(
        ui,
        &format!(
            "The ZAMS line connects models at {} M☉. They are computed in the background with the current X, Z, α_MLT \
             and f_CNO, 400 mesh points each, by continuation from the current model. T_eff increases to the left, as is \
             conventional. The ☉ symbol marks the present-day Sun (5772 K, 1 L☉). The evolutionary track joins every \
             step of the current run; labels mark the ZAMS, the end of the main sequence and the latest age, and \
             hovering over the track shows age, phase, X_c, radius and pulsation. The shaded band is the approximate \
             classical instability strip (see \"Pulsations & oscillations\").",
            {
                let mut ms = ZAMS_MASSES.to_vec();
                ms.sort_by(f64::total_cmp);
                ms.iter().map(|m| format!("{m}")).collect::<Vec<_>>().join(", ")
            }
        ),
    );
    sub(ui, "Time evolution");
    p(
        ui,
        "One point per computed step (not only the timeline models). The energy-source panel shows shares of the surface \
         luminosity: the nuclear shares are the pp/CNO/3α fractions times (1 − L_grav/L), and the gravitational share is \
         ∫ε_grav dm / L, so the four add up to 1. It can be negative while the star expands. Legends are drawn next to \
         each panel title. Logarithmic axes label round values when the range is less than a decade.",
    );
    sub(ui, "Compare");
    p(
        ui,
        "A table of global properties for the current and pinned models, plus an overlay of one chosen profile quantity on the selected x axis.",
    );
}

fn limitations(ui: &mut Ui, dark: bool) {
    p(
        ui,
        "The analytic microphysics gets the structure qualitatively right:\n\
         • a radiative core and convective envelope in Sun-like stars;\n\
         • convective cores above about 1.2 M☉, with CNO dominance above about 1.5 M☉;\n\
         • deep convective envelopes in low-mass stars;\n\
         • massive stars close to detailed models, because their opacity is mostly electron scattering.",
    );
    p(
        ui,
        "It is not quantitatively solar. With X = 0.70 and Z = 0.02 the 1 M☉ model has L ≈ 0.42 L☉, R ≈ 1.12 R☉, \
         T_eff ≈ 4380 K and a thin convective envelope (base at 0.97 R). A detailed ZAMS Sun has L ≈ 0.7 L☉, \
         R ≈ 0.89 R☉ and its envelope base near 0.71 R. The leading suspects, whose shares have not been isolated: the \
         textbook Kramers opacity is too high in the core (the leading suspect for the low luminosity) and far too low \
         at 10⁵–10⁶ K, and the EOS has no ionisation zones.",
    );
    p(
        ui,
        "With Z = 0.008 the same physics reproduces the classic StatStar Sun of Carroll & Ostlie (L ≈ 0.83 L☉, \
         T_eff ≈ 5500 K), which shows that the implementation is consistent with the textbook formulas. This is a check \
         of the code (verification), not a validation against the real Sun.",
    );
    p(
        ui,
        "The numerical checks (virial, energy, residuals, mesh refinement, time-step convergence) show that the \
         equations are solved accurately. They say nothing about how realistic the physics is.",
    );
    sub(ui, "Evolution");
    bullet(
        ui,
        "Main-sequence lifetimes follow the model luminosities. The analytic Sun is underluminous (0.42 L☉ at the ZAMS), \
         so it burns slowly and leaves the main sequence after about 15.1 Gyr instead of about 10 Gyr. Intermediate-mass \
         and massive stars come out within about 25% of detailed models (for example 1.32 Gyr for 1.8 M☉ and 72 Myr \
         for 5 M☉, against roughly 1.2–1.5 Gyr and 90–100 Myr).",
    );
    bullet(
        ui,
        "After the main sequence the results are qualitative. Without electron degeneracy, hydrogen and helium \
         ionisation or molecular opacities, low-mass stars do not climb the red giant branch. They cool towards 2600 K at \
         roughly constant luminosity, the helium flash cannot occur, and runs stop there.",
    );
    bullet(
        ui,
        "The hydrogen-burning shell after the main sequence is thin in mass. Global results are converged at 800 points \
         (5 M☉ crossing: L within 0.2% of 3000 points), but more mesh points resolve the shell's profile better.",
    );
    bullet(
        ui,
        "Core overshoot is a step in the mixing only (α_ov, off by default); its value is a calibration, not a \
         prediction. Not included: diffusive (exponential) overshoot, envelope overshoot, semiconvection, mass loss, \
         rotational mixing, neutrino losses.",
    );
    bullet(
        ui,
        "Pulsation: the period is computed, but driving is not. The instability strip is an empirical band, so stars \
         are flagged by position, not by a stability analysis.",
    );
    sub(ui, "Not modelled");
    bullet(ui, "Partial ionisation of H and He (Saha equation): this affects the convective envelopes of cool stars.");
    bullet(ui, "Electron degeneracy and Coulomb interactions: these matter below about 0.4 M☉.");
    bullet(
        ui,
        "Tabulated opacities (OPAL/OP, low-temperature molecular tables) and detailed model atmospheres (a grey Eddington atmosphere is used).",
    );
    bullet(
        ui,
        "Neutrino losses; rotation; magnetic fields; mass loss; element diffusion; non-radial oscillation modes.",
    );
    bullet(
        ui,
        "Very massive stars (above about 50 M☉) approach the Eddington limit, where this simple treatment becomes crude.",
    );
    bullet(
        ui,
        "General relativity. The structure equations are Newtonian; the relativistic corrections are \
         GM/(Rc²) ≈ 2 × 10⁻⁶ for the Sun and 1.2 × 10⁻⁵ for 60 M☉, below the discretisation error. The gravitational \
         redshift of the surface light is reported as an output (636 m/s for 1 M☉ at 1 R☉).",
    );
    bullet(
        ui,
        "Screening beyond the weak regime. The weak-screening formula is used with its exponent capped at 2; each \
         model reports its largest exponent and warns above 0.3 and at the cap.",
    );
    sub(ui, "Planned");
    p(
        ui,
        "• Tabulated opacities and an EOS with ionisation and degeneracy, plugged into the existing physics interfaces. \
         This is the step that would make the solar model quantitative and allow the red giant branch.\n\
         • Adaptive remeshing around burning shells.\n\
         • A web version.",
    );
    muted(
        ui,
        dark,
        "Use the app to understand how stars are built and how their structure depends on mass and composition, not for precision work.",
    );
}

fn constants(ui: &mut Ui, dark: bool) {
    p(ui, "All calculations use CGS units internally. Solar units are used for display.");
    use stellar_core::constants as k;
    let rows: &[(&str, String)] = &[
        ("Gravitational constant G", format!("{:.5e} cm³ g⁻¹ s⁻²", k::G)),
        ("Speed of light c", format!("{:.8e} cm s⁻¹", k::C)),
        ("Stefan–Boltzmann σ", format!("{:.6e} erg cm⁻² s⁻¹ K⁻⁴", k::SIGMA_SB)),
        ("Radiation constant a = 4σ/c", format!("{:.5e} erg cm⁻³ K⁻⁴", k::A_RAD)),
        ("Boltzmann k", format!("{:.6e} erg K⁻¹", k::K_B)),
        ("Atomic mass unit m_u", format!("{:.6e} g", k::M_U)),
        ("Solar mass M☉", format!("{:.5e} g", k::M_SUN)),
        ("Solar radius R☉ (IAU 2015)", format!("{:.4e} cm", k::R_SUN)),
        ("Solar luminosity L☉ (IAU 2015)", format!("{:.4e} erg s⁻¹", k::L_SUN)),
        ("Solar T_eff (IAU 2015)", format!("{} K", k::TEFF_SUN)),
    ];
    egui::Grid::new("constants").num_columns(2).striped(true).spacing([20.0, 3.0]).show(ui, |ui| {
        for (name, value) in rows {
            ui.label(*name);
            ui.label(RichText::new(value).monospace());
            ui.end_row();
        }
    });
    muted(
        ui,
        dark,
        "Pressure in dyn cm⁻² (1 dyn cm⁻² = 0.1 Pa), density in g cm⁻³, energy generation in erg g⁻¹ s⁻¹ (1 erg g⁻¹ s⁻¹ = 10⁻⁴ W kg⁻¹).",
    );
}

fn cli(ui: &mut Ui, dark: bool) {
    p(ui, "The same solver is available headless, for scripting and batch runs:");
    eq(ui, "cargo run --release -p stellar-cli -- [options]");
    let opts: &[(&str, &str)] = &[
        ("--mass <M>", "stellar mass in M☉ (default 1)"),
        ("--x <X>", "hydrogen mass fraction (default 0.70)"),
        ("--z <Z>", "metallicity (default 0.02)"),
        ("--alpha <α>", "mixing-length parameter (default 1.8)"),
        ("--f-cno <f>", "CNO fraction of Z (default 0.7)"),
        ("--points <N>", "mesh points (default 800)"),
        ("-o, --output <file>", "write the profile CSV"),
        ("--json <file>", "write the summary JSON"),
        ("--sweep <M1,M2,…>", "solve a ZAMS sequence by continuation and print a table"),
        ("--evolve", "evolve the star from the ZAMS and print the track"),
        ("--until <tams|age|max>", "stop condition: end of main sequence (default), age, or as far as possible"),
        ("--max-age <Gyr>", "target age with --until age; otherwise a safety limit"),
        ("--track <file>", "write the evolutionary track (one row per model) as CSV"),
        ("-v, --verbose", "print every Newton iteration (every model with --evolve)"),
        ("-h, --help / -V, --version", "print help / version"),
    ];
    egui::Grid::new("cli_opts").num_columns(2).striped(true).spacing([20.0, 3.0]).show(ui, |ui| {
        for (o, d) in opts {
            ui.label(RichText::new(*o).monospace());
            ui.label(*d);
            ui.end_row();
        }
    });
    sub(ui, "Examples");
    eq(ui, "stellar-cli --mass 1 -o sun.csv --json sun.json");
    eq(ui, "stellar-cli --sweep 0.5,1,2,5,10,20");
    eq(ui, "stellar-cli --mass 3 --z 0.008 -v");
    eq(ui, "stellar-cli --mass 1 --evolve --until age --max-age 4.57");
    eq(ui, "stellar-cli --mass 5 --evolve --until max --track track.csv -o final_profile.csv");
    muted(
        ui,
        dark,
        "Both the CLI and the app use the stellar-core library, so identical parameters give the same model, to within the solver \
         tolerance. The app warm-starts from the previous model; the CLI starts cold unless --sweep is used.",
    );
}

fn references(ui: &mut Ui, dark: bool) {
    sub(ui, "Physics and methods");
    term(
        ui,
        "Kippenhahn, Weigert & Weiss (2012),",
        "Stellar Structure and Evolution, 2nd ed., Springer: the structure equations, mixing-length formulation, Henyey method and polytropes.",
    );
    term(
        ui,
        "Kippenhahn & Weigert (1990),",
        "Stellar Structure and Evolution, 1st ed., Springer: the triple-α rate formula (and the pp and CNO formulas \
         used before the Solar Fusion III rates).",
    );
    term(
        ui,
        "Carroll & Ostlie (2017),",
        "An Introduction to Modern Astrophysics, 2nd ed., Cambridge University Press: the analytic opacity formulas and the StatStar model.",
    );
    term(ui, "Hansen, Kawaler & Trimble (2004),", "Stellar Interiors, 2nd ed., Springer.");
    term(ui, "Henyey, Forbes & Gould (1964),", "ApJ 139, 306: the relaxation method for stellar structure.");
    term(ui, "Böhm-Vitense (1958),", "Zeitschrift für Astrophysik 46, 108: mixing-length theory.");
    term(
        ui,
        "Clayton (1968),",
        "Principles of Stellar Evolution and Nucleosynthesis, McGraw-Hill: thermonuclear rates and weak screening.",
    );
    term(ui, "Asplund, Grevesse, Sauval & Scott (2009),", "ARA&A 47, 481: the AGSS09 solar composition preset.");
    term(
        ui,
        "Acharya et al. (2025),",
        "Rev. Mod. Phys. 97, 035002 (\"Solar Fusion III\"): the astrophysical factors S(0), S′(0), S″(0) behind the \
         pp and CNO rates and the pp II branching factor ψ.",
    );
    term(
        ui,
        "Paxton et al. (2011),",
        "ApJS 192, 3 (MESA), §5.1: the bounded pressure scale height in the mixing length.",
    );
    term(
        ui,
        "Choi et al. (2016),",
        "ApJ 823, 102 (MIST), §3.1: Y_p and ΔY/ΔZ of the \"helium follows metallicity\" mode.",
    );
    term(ui, "IAU 2015 Resolution B3:", "nominal solar radius, luminosity and effective temperature.");
    term(
        ui,
        "Cox (1980),",
        "Theory of Stellar Pulsation, Princeton University Press: the linear adiabatic wave equation and its boundary conditions.",
    );
    term(
        ui,
        "Pekeris (1938),",
        "ApJ 88, 189: the radial modes of the uniform-density sphere, used to test the pulsation solver.",
    );
    term(
        ui,
        "Unno, Osaki, Ando, Saio & Shibahashi (1989),",
        "Nonradial Oscillations of Stars, 2nd ed., University of Tokyo Press.",
    );
    term(
        ui,
        "Brown, Gilliland, Noyes & Ramsey (1991),",
        "ApJ 368, 599, and Kjeldsen & Bedding (1995), A&A 293, 87: the ν_max and Δν scaling relations.",
    );
    sub(ui, "Software");
    term(
        ui,
        "Rust:",
        "programming language; the physics and solver are written from scratch in the stellar-core crate.",
    );
    term(
        ui,
        "egui / eframe 0.37 and egui_plot 0.38:",
        "user interface, plots and window (native on Linux, macOS and Windows).",
    );
    term(
        ui,
        "rfd:",
        "native file dialogs. png: screenshot encoding. serde / serde_json: JSON export. clap: CLI arguments.",
    );
    ui.add_space(6.0);
    ui.horizontal(|ui| {
        ui.hyperlink_to("Rust", "https://www.rust-lang.org");
        ui.label("·");
        ui.hyperlink_to("egui", "https://github.com/emilk/egui");
        ui.label("·");
        ui.hyperlink_to("egui_plot", "https://github.com/emilk/egui_plot");
    });
    muted(
        ui,
        dark,
        "Black-body surface colours in the 3D view come from an approximate temperature-to-sRGB table (illustrative only).",
    );
}

#[cfg(test)]
mod tests {
    #[test]
    fn every_referenced_equation_exists() {
        let src = include_str!("about.rs");
        let mut count = 0;
        for part in src.split("math(ui, \"").skip(1) {
            let id = part.split('"').next().unwrap();
            if id.is_empty() || id.contains(' ') {
                continue;
            }
            assert!(crate::equations::get(id).is_some(), "equation `{id}` missing: run tools/equations");
            count += 1;
        }
        assert!(count >= 40, "expected the About page to reference the equations, found {count}");
    }

    /// docs/about/ holds this page as Markdown for GitHub; every equation must
    /// appear there verbatim as a ```math block, so the two cannot drift apart.
    #[test]
    fn markdown_docs_contain_every_equation() {
        let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/about");
        let mut docs = String::new();
        for entry in std::fs::read_dir(&dir).expect("docs/about exists") {
            let path = entry.unwrap().path();
            if path.extension().is_some_and(|e| e == "md") {
                // Windows checkouts may convert to CRLF; compare line endings-agnostically.
                docs += &std::fs::read_to_string(&path).unwrap().replace("\r\n", "\n");
            }
        }
        for e in crate::equations::EQUATIONS {
            let block = format!("```math\n{}\n```", e.tex);
            assert!(
                docs.contains(&block),
                "equation `{}` is missing from docs/about (or differs from equations.json)",
                e.id
            );
        }
    }

    #[test]
    fn equation_svgs_are_self_contained() {
        for e in crate::equations::EQUATIONS {
            let svg = std::str::from_utf8(e.svg).unwrap();
            assert!(svg.starts_with("<svg"), "{}", e.id);
            assert!(!svg.contains("currentColor"), "{} must use white glyphs for tinting", e.id);
            assert!(!svg.contains("merror"), "{} contains a TeX error", e.id);
            assert!(e.width > 0.0 && e.height > 0.0);
        }
    }
}
