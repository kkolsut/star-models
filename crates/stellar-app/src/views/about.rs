//! "About" tab: how the models are computed, what every parameter and output
//! means, how to use each view, limitations and references.
//!
//! Numbers quoted here come from the code (stellar-core). Where a value is
//! available at run time (current parameters, physics in use, ZAMS masses,
//! CSV columns, version) it is read live, so the text cannot drift from the
//! implementation.

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
        section(ui, 5, &mut |ui| equations(ui, dark));
        section(ui, 6, &mut |ui| microphysics(ui, dark));
        section(ui, 7, &mut |ui| transport(ui, dark));
        section(ui, 8, &mut |ui| boundaries(ui, dark));
        section(ui, 9, &mut |ui| numerics(ui, dark));
        section(ui, 10, &mut |ui| initial(ui, dark));
        section(ui, 11, &mut |ui| views(ui, dark));
        section(ui, 12, &mut |ui| limitations(ui, dark));
        section(ui, 13, &mut |ui| constants(ui, dark));
        section(ui, 14, &mut |ui| cli(ui, dark));
        section(ui, 15, &mut |ui| references(ui, dark));
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
        "Every model is static (no time evolution), spherically symmetric, non-rotating and non-magnetic. \
         The defaults are those of the Sun (1 M☉, X = 0.70, Z = 0.02, α_MLT = 1.8). The default model is \
         therefore the young Sun as it was about 4.6 billion years ago, not today's Sun.",
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
         convective velocity, convective flux fraction and gas-pressure fraction. The x axis can also be log P, \
         which spreads out the thin outer envelope.",
    );
    term(
        ui,
        "Cross-section:",
        "a slice through the star, with temperature on the left half and density on the right. Convective zones carry animated cells.",
    );
    term(ui, "3D model:", "a cut-away sphere with animated radiative and convective energy transport.");
    term(
        ui,
        "HR diagram:",
        "the model on a computed ZAMS line, together with the present-day Sun and the pinned models.",
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
        "Not set directly; it is whatever mass is left after hydrogen and metals.",
    );
    param(
        ui,
        dark,
        "Mixing-length parameter  α_MLT",
        &format!("default 1.8 · range 0.3 – 3.0 · current {:.2}", pr.alpha_mlt),
        "In mixing-length theory a convective blob travels a distance l = α_MLT · H_P (H_P is the pressure scale height) \
         before it dissolves into its surroundings. α_MLT only matters where convection is inefficient: in the \
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
        "for any star in hydrostatic equilibrium the virial theorem makes this ratio exactly 1, so any deviation \
         measures the discretisation error (typically 10⁻⁵ to 10⁻⁴):",
    );
    math(ui, "res_virial", "");
    term(
        ui,
        "Energy check:",
        "the energy generated inside must equal the energy leaving the surface. L₀ = ε₀ m₀ is the luminosity of the innermost point:",
    );
    math(ui, "res_energy", "");
    term(ui, "Newton iterations:", "iterations in the last solver step (a continuation may take several steps).");
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
        ("nabla_rad", "radiative gradient ∇_rad: the d ln T / d ln P needed to carry L by radiation alone"),
        ("nabla_ad", "adiabatic gradient ∇_ad = (∂ ln T / ∂ ln P) at constant entropy"),
        ("nabla", "actual gradient ∇ = d ln T / d ln P used in the model"),
        ("convective", "1 if convective by the Schwarzschild criterion (∇_rad > ∇_ad), else 0"),
        ("v_conv_cm_s", "convective velocity from mixing-length theory [cm s⁻¹] (0 in radiative layers)"),
        ("F_conv_frac", "fraction of the energy flux carried by convection"),
        ("H_P_cm", "pressure scale height H_P = P / (ρ g) [cm]"),
        ("mu", "mean molecular weight μ"),
        ("beta", "gas-pressure fraction β = P_gas / P"),
        ("X", "hydrogen mass fraction (uniform in ZAMS models)"),
        ("Y", "helium mass fraction"),
        ("Z", "metal mass fraction"),
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
        "params: the input parameters (mass, x, z, alpha_mlt, f_cno, n_points).\n\
         physics: the description of the physics used.\n\
         globals: every quantity of \"Results: global quantities\", plus newton_iterations and max_correction.\n\
         convective_zones: one entry per zone, with q_inner, q_outer, r_inner and r_outer (r in units of R).",
    );
    muted(
        ui,
        dark,
        "All values are in CGS units unless the column name says otherwise. Numbers are written with 10 significant digits.",
    );
}

fn equations(ui: &mut Ui, dark: bool) {
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
        "Static: the gravitational energy term ε_grav = −T ds/dt is zero, as appropriate for a star in thermal equilibrium on the ZAMS.",
    );
    bullet(ui, "Chemically homogeneous: X, Y, Z are the same everywhere (no burning history yet).");
    bullet(ui, "No neutrino losses, rotation, magnetic fields, mass loss or element diffusion.");
    muted(
        ui,
        dark,
        "Composition is stored per mesh point, so composition gradients from evolution can be added later without changing the solver.",
    );
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

    sub(ui, "Nuclear energy generation [erg g⁻¹ s⁻¹]");
    math(ui, "nuc_t6", "");
    math(ui, "nuc_pp", "pp chain");
    math(ui, "nuc_psi", "");
    math(ui, "nuc_cno", "CNO cycle");
    math(ui, "nuc_c", "");
    math(ui, "nuc_3a", "triple-α");
    p(ui, "Weak electron screening (Salpeter / Clayton):");
    math(ui, "nuc_screen", "");
    p(
        ui,
        "The charge products are Z₁Z₂ = 1 (p + p), 7 (p + ¹⁴N, the slowest CNO step), and 4 and 8 (α + α, α + ⁸Be). \
         The exponent is capped at 2, since weak screening is only valid for small values.",
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
         and keeps rising. The Ledoux criterion is not needed, because ZAMS models have no composition gradients.",
    );
    sub(ui, "Convective transport: mixing-length theory (MLT)");
    p(ui, "Following Böhm-Vitense (1958) in the formulation of Kippenhahn, Weigert & Weiss (2012, §7):");
    math(ui, "mlt_defs", "gravity, pressure scale height, mixing length");
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
    muted(ui, dark, "Not included: convective overshooting, semiconvection, turbulent pressure.");
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
         iterations, with quadratic convergence at the end. The limit is 250 iterations from a cold start and 60 from a \
         warm start.",
    );
    bullet(
        ui,
        "Every converged model is checked: all values finite, r increasing outward, P decreasing outward, L > 0, and β \
         above the floor. A model failing any check is rejected.",
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
        "Rendering: orthographic projection on the CPU. About 40,000 shaded triangles are depth-sorted every frame \
         (painter's algorithm, about 2 ms) and drawn as one mesh with 4× multisampling. Overlays use an exact ray test \
         against the cut-away geometry, so they are hidden correctly behind the sphere and the core.",
    );
    muted(
        ui,
        dark,
        "Animation speeds are illustrative, not to scale. In the Sun a photon needs of order 10⁴–10⁵ years to diffuse out \
         of the radiative zone, while convective turnover takes from minutes (granules) to about a month (base of the envelope).",
    );
    sub(ui, "HR diagram");
    p(
        ui,
        &format!(
            "The ZAMS line connects models at {} M☉. They are computed in the background with the current X, Z, α_MLT \
             and f_CNO, 400 mesh points each, by continuation from the current model. T_eff increases to the left, as is \
             conventional. The ☉ symbol marks the present-day Sun (5772 K, 1 L☉).",
            {
                let mut ms = ZAMS_MASSES.to_vec();
                ms.sort_by(f64::total_cmp);
                ms.iter().map(|m| format!("{m}")).collect::<Vec<_>>().join(", ")
            }
        ),
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
        "It is not quantitatively solar. With X = 0.70 and Z = 0.02 the 1 M☉ model has L ≈ 0.42 L☉, R ≈ 1.11 R☉, \
         T_eff ≈ 4400 K and a thin convective envelope (base at 0.97 R). A detailed ZAMS Sun has L ≈ 0.7 L☉, \
         R ≈ 0.89 R☉ and its envelope base near 0.71 R. The textbook Kramers opacity is too high in the core and far \
         too low at 10⁵–10⁶ K, and the EOS has no ionisation zones.",
    );
    p(
        ui,
        "With Z = 0.008 the same physics reproduces the classic StatStar Sun of Carroll & Ostlie (L ≈ 0.83 L☉, \
         T_eff ≈ 5500 K), which confirms that the implementation matches the textbook formulas.",
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
        "Neutrino losses; time evolution (age, burnt composition, ε_grav); rotation; magnetic fields; mass loss; diffusion; overshooting.",
    );
    bullet(
        ui,
        "Very massive stars (above about 50 M☉) approach the Eddington limit, where this simple treatment becomes crude.",
    );
    sub(ui, "Planned");
    p(
        ui,
        "• Tabulated opacities and an EOS with ionisation, plugged into the existing physics interfaces.\n\
         • Time evolution, for example the Sun at 4.57 Gyr.\n\
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
        ("-v, --verbose", "print every Newton iteration"),
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
        "Carroll & Ostlie (2017),",
        "An Introduction to Modern Astrophysics, 2nd ed., Cambridge University Press: the analytic opacity and nuclear-rate formulas, and the StatStar model.",
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
    term(ui, "IAU 2015 Resolution B3:", "nominal solar radius, luminosity and effective temperature.");
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
