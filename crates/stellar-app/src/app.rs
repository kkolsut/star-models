//! Application state, side panel (parameters + results), tabs and status bar.

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use egui::{Color32, RichText};
use stellar_core::output::{profile_csv, summary_json};
use stellar_core::{NoControl, Physics, SolveError, SolveEvent, StarModel, StarParams, solve_model};

use crate::evolution::{EvoAction, EvolutionRun, EvolutionSettings, ScaleMode, StopChoice};
use crate::theme;
use crate::views::about::AboutState;
use crate::views::compare::{CompareAction, CompareState};
use crate::views::cross_section::CrossSectionState;
use crate::views::history::{HistoryAction, HistoryState};
use crate::views::hr::ZamsLine;
use crate::views::model3d::Model3dState;
use crate::views::{PinnedModel, ViewCtx, XAxis, model_label, swatch};
use crate::worker::{Msg, Worker, ZamsKey};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tab {
    Profiles,
    Extra,
    CrossSection,
    Model3d,
    Hr,
    History,
    Compare,
    About,
}

impl Tab {
    const ALL: [Tab; 8] =
        [Tab::Profiles, Tab::Extra, Tab::CrossSection, Tab::Model3d, Tab::Hr, Tab::History, Tab::Compare, Tab::About];

    fn title(self) -> &'static str {
        match self {
            Tab::Profiles => "Profiles (P, T, ρ)",
            Tab::Extra => "More profiles",
            Tab::CrossSection => "Cross-section",
            Tab::Model3d => "3D model",
            Tab::Hr => "HR diagram",
            Tab::History => "Time evolution",
            Tab::Compare => "Compare",
            Tab::About => "About",
        }
    }

    fn slug(self) -> &'static str {
        match self {
            Tab::Profiles => "profiles",
            Tab::Extra => "extra",
            Tab::CrossSection => "cross_section",
            Tab::Model3d => "model3d",
            Tab::Hr => "hr",
            Tab::History => "history",
            Tab::Compare => "compare",
            Tab::About => "about",
        }
    }
}

/// Automated screenshot mode (STAR_MODELS_CAPTURE=<dir>): visits every tab,
/// saves a PNG of each, then quits. Used to check the rendering.
struct Capture {
    dir: PathBuf,
    tabs: Vec<Tab>,
    step: usize,
    frames: u32,
    requested: bool,
    saved: Arc<AtomicBool>,
    started: f64,
}

/// Present-day Sun, for orientation next to the model values.
struct SunToday;
impl SunToday {
    const R: f64 = 1.0;
    const L: f64 = 1.0;
    const TEFF: f64 = 5772.0;
    const TC: f64 = 1.57e7;
    const RHOC: f64 = 150.0;
    const ENV_BASE: f64 = 0.713;
}

pub struct StarApp {
    params: StarParams,
    model: Option<Arc<StarModel>>,
    worker: Worker,
    auto_solve: bool,
    changed_at: f64,
    pending: Option<(u64, StarParams)>,
    failed: Option<StarParams>,
    stage: String,
    iteration: Option<(usize, f64)>,
    last_error: Option<String>,
    last_elapsed: Option<Duration>,
    pinned: Vec<PinnedModel>,
    tab: Tab,
    /// x axis of the P/T/ρ tab (radius or mass).
    profile_axis: XAxis,
    /// x axis of the detail and compare tabs (radius, mass or log P).
    detail_axis: XAxis,
    cross: CrossSectionState,
    model3d: Model3dState,
    about: AboutState,
    history: HistoryState,
    compare: CompareState,
    zams: ZamsLine,
    notice: Option<String>,
    capture: Option<Capture>,
    /// Current evolution run (track and timeline), if any.
    evolution: Option<EvolutionRun>,
    evo_settings: EvolutionSettings,
    /// How views scale while an evolution run is shown.
    scale_mode: ScaleMode,
    /// Scale mode of the previous frame (a change resets manual plot zoom).
    last_scale_mode: ScaleMode,
    /// Capture mode: evolution to run before taking screenshots, and which model to show.
    capture_evolve: Option<(StopChoice, f64, String)>,
    capture_evolve_started: bool,
}

impl StarApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
        // SVG loader for the typeset equations of the About page.
        egui_extras::install_image_loaders(&cc.egui_ctx);
        let physics = Arc::new(Physics::analytic());
        let mut params = StarParams::sun();
        if let Some(m) = std::env::var("STAR_MODELS_MASS").ok().and_then(|s| s.parse().ok()) {
            params.mass = m;
        }
        match std::env::var("STAR_MODELS_THEME").as_deref() {
            Ok("light") => cc.egui_ctx.set_theme(egui::Theme::Light),
            Ok("dark") => cc.egui_ctx.set_theme(egui::Theme::Dark),
            _ => {}
        }

        // Optional pre-pinned models (used by the capture mode and for demos).
        let mut pinned = Vec::new();
        if let Ok(list) = std::env::var("STAR_MODELS_PIN") {
            let sun = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).ok();
            for (i, m) in list.split(',').filter_map(|s| s.trim().parse::<f64>().ok()).enumerate() {
                let p = StarParams { mass: m, ..StarParams::sun() };
                if let Ok(model) = solve_model(&p, &physics, sun.as_ref(), &mut NoControl) {
                    pinned.push(PinnedModel {
                        label: model_label(&model),
                        model: Arc::new(model),
                        slot: i + 1,
                        visible: true,
                    });
                }
            }
        }

        let tabs: Vec<Tab> = match std::env::var("STAR_MODELS_TABS") {
            Ok(list) => Tab::ALL.into_iter().filter(|t| list.split(',').any(|s| s.trim() == t.slug())).collect(),
            Err(_) => Tab::ALL.to_vec(),
        };
        let capture = std::env::var("STAR_MODELS_CAPTURE").ok().map(|d| Capture {
            dir: PathBuf::from(d),
            tabs,
            step: 0,
            frames: 0,
            requested: false,
            saved: Arc::new(AtomicBool::new(false)),
            started: 0.0,
        });

        let detail_axis = match std::env::var("STAR_MODELS_AXIS").as_deref() {
            Ok("pressure") => XAxis::Pressure,
            Ok("mass") => XAxis::Mass,
            _ => XAxis::Radius,
        };

        let scale_mode = if std::env::var("STAR_MODELS_SCALE").as_deref() == Ok("adapt") {
            ScaleMode::AdaptToModel
        } else {
            ScaleMode::FixedForRun
        };

        let mut app = Self {
            params,
            model: None,
            worker: Worker::new(cc.egui_ctx.clone(), physics),
            auto_solve: true,
            changed_at: -1.0,
            pending: None,
            failed: None,
            stage: String::new(),
            iteration: None,
            last_error: None,
            last_elapsed: None,
            pinned,
            tab: Tab::Profiles,
            profile_axis: XAxis::Radius,
            detail_axis,
            cross: CrossSectionState::default(),
            model3d: Model3dState::default(),
            history: HistoryState::default(),
            about: std::env::var("STAR_MODELS_ABOUT_SECTION")
                .ok()
                .and_then(|v| v.parse().ok())
                .map(AboutState::at_section)
                .unwrap_or_default(),
            compare: CompareState::default(),
            zams: ZamsLine::default(),
            notice: None,
            capture,
            evolution: None,
            evo_settings: EvolutionSettings::default(),
            scale_mode,
            last_scale_mode: scale_mode,
            capture_evolve: std::env::var("STAR_MODELS_EVOLVE").ok().map(|v| {
                let select = std::env::var("STAR_MODELS_EVO_SELECT").unwrap_or_else(|_| "latest".into());
                match v.split_once(':') {
                    Some(("age", a)) => (StopChoice::Age, a.parse().unwrap_or(4.57), select),
                    _ if v == "max" => (StopChoice::AsFarAsPossible, 4.57, select),
                    _ => (StopChoice::EndOfMainSequence, 4.57, select),
                }
            }),
            capture_evolve_started: false,
        };
        match std::env::var("STAR_MODELS_3D_FACES").as_deref() {
            Ok("hydrogen") => app.model3d.quantity = crate::views::model3d::FaceQuantity::Hydrogen,
            Ok("density") => app.model3d.quantity = crate::views::model3d::FaceQuantity::Density,
            _ => {}
        }
        if std::env::var("STAR_MODELS_HISTORY_FIT").is_ok() {
            app.history.fit_all = true;
        }
        if std::env::var("STAR_MODELS_HISTORY_AXIS").as_deref() == Ok("step") {
            app.history.axis = crate::views::history::HistoryAxis::Step;
        }
        if std::env::var("STAR_MODELS_TRUE_SIZE").is_ok() {
            app.cross.true_size = true;
            app.model3d.true_size = true;
        }
        if let Some(t) = std::env::var("STAR_MODELS_3D_TIME").ok().and_then(|v| v.parse().ok()) {
            app.model3d.freeze_at(t);
        }
        if std::env::var("STAR_MODELS_3D_ZOOM").as_deref() == Ok("surface") {
            app.model3d.zoom_to_surface();
        }
        app.solve();
        app
    }

    fn solve(&mut self) {
        if let Err(e) = self.params.validate() {
            self.last_error = Some(e);
            return;
        }
        let job = self.worker.submit(self.params.clone(), self.model.clone());
        self.pending = Some((job, self.params.clone()));
        self.failed = None;
        self.stage = "starting".into();
        self.iteration = None;
    }

    fn handle_messages(&mut self) {
        for msg in self.worker.poll() {
            match msg {
                Msg::Progress { job, event } if self.worker.is_current(job) => match event {
                    SolveEvent::Stage(s) => {
                        self.stage = s;
                        self.iteration = None;
                    }
                    SolveEvent::Iteration { iter, max_correction, .. } => self.iteration = Some((iter, max_correction)),
                },
                Msg::Progress { .. } => {}
                Msg::Finished { job, result, elapsed } => {
                    if !self.worker.is_current(job) {
                        continue;
                    }
                    self.worker.finish(job);
                    let params = self.pending.take().map(|(_, p)| p);
                    match *result {
                        Ok(m) => {
                            if self.evolution.as_ref().is_some_and(|r| r.snapshots[0].params != m.params) {
                                self.evolution = None;
                            }
                            self.model = Some(Arc::new(m));
                            self.last_error = None;
                            self.last_elapsed = Some(elapsed);
                            self.ensure_zams();
                        }
                        Err(SolveError::Cancelled) => {}
                        Err(e) => {
                            self.last_error = Some(e.to_string());
                            self.failed = params;
                        }
                    }
                }
                Msg::Zams { key, point } => {
                    if self.zams.key == Some(key) {
                        self.zams.points.push(point);
                    }
                }
                Msg::ZamsDone { key } => {
                    if self.zams.key == Some(key) {
                        self.zams.done = true;
                    }
                }
                Msg::Evolved { job, model } => {
                    if let Some(run) = &mut self.evolution
                        && run.job == Some(job)
                    {
                        run.push(*model);
                    }
                }
                Msg::EvolutionDone { job, reason, dt_yr, r_zams } => {
                    self.worker.finish(job);
                    if let Some(run) = &mut self.evolution
                        && run.job == Some(job)
                    {
                        run.job = None;
                        run.finished = Some(reason);
                        run.resume = Some((dt_yr, r_zams));
                    }
                }
            }
        }
    }

    /// Start computing the ZAMS line if the composition changed.
    fn ensure_zams(&mut self) {
        let Some(model) = &self.model else { return };
        let key = ZamsKey::of(&model.params);
        if self.zams.key != Some(key) {
            self.zams = ZamsLine { key: Some(key), points: Vec::new(), done: false };
            self.worker.submit_zams(model.params.clone(), model.clone());
        }
    }

    fn maybe_auto_solve(&mut self, now: f64) {
        if !self.auto_solve || self.params.validate().is_err() {
            return;
        }
        let solved = self.model.as_ref().map(|m| &m.params) == Some(&self.params);
        let in_flight = self.pending.as_ref().map(|(_, p)| p) == Some(&self.params);
        let failed = self.failed.as_ref() == Some(&self.params);
        if !solved && !in_flight && !failed && now - self.changed_at > 0.15 {
            self.solve();
        }
    }

    fn pin_current(&mut self) {
        let Some(model) = &self.shown_model() else { return };
        let used: Vec<usize> = self.pinned.iter().map(|p| p.slot).collect();
        let Some(slot) = (1..theme::MAX_SLOTS).find(|s| !used.contains(s)) else {
            self.notice = Some("At most 7 pinned models (one per palette color). Unpin one first.".into());
            return;
        };
        self.pinned.push(PinnedModel { model: model.clone(), slot, label: model_label(model), visible: true });
    }

    /// The model shown in the views: the selected age of an evolution run,
    /// otherwise the ZAMS model.
    fn shown_model(&self) -> Option<Arc<StarModel>> {
        self.evolution.as_ref().map(|r| r.selected_model().clone()).or_else(|| self.model.clone())
    }

    fn evo_action(&mut self, action: EvoAction) {
        match action {
            EvoAction::SunAge => {
                self.evo_settings.stop = StopChoice::Age;
                self.evo_settings.target_age_gyr = 4.57;
                self.start_evolution();
            }
            EvoAction::Start => self.start_evolution(),
            EvoAction::Stop => self.worker.cancel(),
            EvoAction::Reset => self.evolution = None,
        }
    }

    fn history_action(&mut self, action: HistoryAction) {
        match action {
            HistoryAction::Evolve => self.start_evolution(),
            HistoryAction::SelectAge(age) => {
                if let Some(run) = &mut self.evolution {
                    let best = run
                        .snapshots
                        .iter()
                        .enumerate()
                        .min_by(|a, b| (a.1.globals.age_yr - age).abs().total_cmp(&(b.1.globals.age_yr - age).abs()))
                        .map(|(i, _)| i);
                    if let Some(i) = best {
                        run.selected = i;
                        run.follow = false;
                        run.playing = false;
                    }
                }
            }
            HistoryAction::Export => {
                if let Some(run) = &self.evolution {
                    let stem = format!("star_{:.2}Msun_history.csv", run.snapshots[0].params.mass);
                    self.notice = save_text(&stem, "CSV", "csv", &crate::evolution::history_csv(&run.track));
                }
            }
        }
    }

    /// Start a new run from the ZAMS model, or continue the current run from
    /// its latest model when possible.
    fn start_evolution(&mut self) {
        let Some(zams) = self.model.clone() else { return };
        if zams.params != self.params {
            self.notice = Some("Wait for the ZAMS model of the current parameters first.".into());
            return;
        }
        let opts = self.evo_settings.options();
        if let Some(run) = &mut self.evolution
            && !run.running()
            && let Some(resume) = run.resume
        {
            let latest_age = run.latest().globals.age_yr;
            let target_reached = matches!(opts.stop, stellar_core::StopCondition::Age(a) if a <= latest_age);
            if !target_reached {
                run.job = Some(self.worker.submit_evolution(run.latest().clone(), Some(resume), opts));
                run.finished = None;
                run.follow = true;
                return;
            }
        }
        let mut run = EvolutionRun::new(zams.clone());
        run.job = Some(self.worker.submit_evolution(zams, None, opts));
        self.evolution = Some(run);
    }

    // ------------------------------------------------------------------ panels

    fn params_panel(&mut self, ui: &mut egui::Ui) {
        let dark = ui.visuals().dark_mode;
        ui.add_space(4.0);
        ui.heading("Star Models");
        ui.label(RichText::new("Stellar structure and evolution").color(theme::secondary(dark)));
        ui.label(
            RichText::new(format!(
                "Physics: analytic Kramers + H⁻ opacity, fully ionized ideal gas + radiation, MLT (α = {:.2}), grey \
                 atmosphere. Qualitative, not quantitatively solar.",
                self.params.alpha_mlt
            ))
            .small()
            .color(theme::muted(dark)),
        )
        .on_hover_text(format!(
            "{}\n\nComposition: X, Y, Z with a fixed CNO fraction of Z. See About › Accuracy & limitations.",
            self.worker.physics().describe()
        ));
        ui.separator();

        let before = self.params.clone();
        ui.horizontal_wrapped(|ui| {
            ui.label("Preset:");
            if ui.button("Sun (X 0.70, Z 0.02)").clicked() {
                self.params = StarParams { n_points: self.params.n_points, ..StarParams::sun() };
            }
            if ui
                .button("Sun AGSS09")
                .on_hover_text("protosolar X = 0.7154, Z = 0.0142 (Asplund et al. 2009)")
                .clicked()
            {
                self.params = StarParams { n_points: self.params.n_points, ..StarParams::sun_agss09() };
            }
        });
        ui.add_space(4.0);
        egui::Grid::new("params").num_columns(2).spacing([8.0, 6.0]).show(ui, |ui| {
            ui.label("Mass").on_hover_text("total stellar mass");
            ui.add(
                egui::Slider::new(&mut self.params.mass, 0.1..=100.0).logarithmic(true).suffix(" M☉").max_decimals(3),
            );
            ui.end_row();

            ui.label("X (hydrogen)");
            let x_max = 1.0 - self.params.z;
            ui.add(egui::Slider::new(&mut self.params.x, 0.0..=x_max).max_decimals(4));
            ui.end_row();

            ui.label("Z (metals)");
            ui.add(egui::Slider::new(&mut self.params.z, 1e-4..=0.1).logarithmic(true).max_decimals(5));
            ui.end_row();

            ui.label("Y (helium)");
            ui.label(RichText::new(format!("{:.4}  = 1 − X − Z", self.params.y())).monospace());
            ui.end_row();

            ui.label("α_MLT").on_hover_text("mixing-length parameter: mixing length = α · pressure scale height");
            ui.add(egui::Slider::new(&mut self.params.alpha_mlt, 0.3..=3.0).max_decimals(2));
            ui.end_row();
        });
        egui::CollapsingHeader::new("Advanced").show(ui, |ui| {
            egui::Grid::new("adv").num_columns(2).spacing([8.0, 6.0]).show(ui, |ui| {
                ui.label("CNO fraction of Z")
                    .on_hover_text("mass fraction of metals in C, N and O (feeds the CNO cycle)");
                ui.add(egui::Slider::new(&mut self.params.f_cno, 0.0..=1.0).max_decimals(2));
                ui.end_row();
                ui.label("Mesh points");
                ui.add(egui::Slider::new(&mut self.params.n_points, 200..=3000).logarithmic(true));
                ui.end_row();
            });
        });
        // Keep X + Z ≤ 1 after Z changes.
        self.params.x = self.params.x.min(1.0 - self.params.z);
        if self.params != before {
            self.changed_at = ui.input(|i| i.time);
        }

        ui.add_space(4.0);
        ui.horizontal(|ui| {
            ui.checkbox(&mut self.auto_solve, "Auto-solve");
            if ui.button("Solve").clicked() {
                self.solve();
            }
            if self.worker.busy() && ui.button("Cancel").clicked() {
                self.worker.cancel();
                self.pending = None;
            }
        });
        if let Some(w) = self.params.ideal_gas_warning() {
            ui.label(RichText::new(format!("Note: {w}")).small().color(theme::secondary(dark)));
        }
        if let Err(e) = self.params.validate() {
            ui.colored_label(Color32::from_rgb(0xd0, 0x3b, 0x3b), format!("Invalid parameters: {e}"));
        }

        ui.separator();
        let prefs = crate::evolution::DisplayPrefs {
            scale: &mut self.scale_mode,
            true_size_cross: &mut self.cross.true_size,
            true_size_3d: &mut self.model3d.true_size,
        };
        if let Some(action) = crate::evolution::controls(ui, &mut self.evo_settings, self.evolution.as_ref(), prefs) {
            self.evo_action(action);
        }
        ui.separator();
        self.results(ui);
    }

    fn results(&mut self, ui: &mut egui::Ui) {
        let dark = ui.visuals().dark_mode;
        let Some(model) = self.shown_model() else {
            ui.label("No model yet.");
            return;
        };
        let g = &model.globals;
        let stale = model.params != self.params;
        ui.horizontal(|ui| {
            swatch(ui, theme::series(0, dark));
            ui.label(RichText::new(format!("Current model: {}", model_label(&model))).strong());
        });
        if stale {
            ui.label(
                RichText::new("Parameters changed; showing the previous model.").small().color(theme::secondary(dark)),
            );
        }
        ui.add_space(2.0);
        egui::Grid::new("results").num_columns(3).striped(true).spacing([10.0, 3.0]).show(ui, |ui| {
            ui.label("");
            ui.label(RichText::new("model").strong());
            ui.label(RichText::new("Sun today").color(theme::muted(dark)))
                .on_hover_text("present-day Sun (4.6 Gyr old) for orientation; a ZAMS model is the young star");
            ui.end_row();
            let mut row = |name: &str, tip: &str, val: String, sun: String| {
                let label = ui.label(name);
                if !tip.is_empty() {
                    label.on_hover_text(tip);
                }
                ui.monospace(val);
                ui.label(RichText::new(sun).monospace().color(theme::muted(dark)));
                ui.end_row();
            };
            row("R [R☉]", "", format!("{:.4}", g.radius), format!("{:.3}", SunToday::R));
            row("L [L☉]", "", crate::views::fmt_value(g.luminosity), format!("{:.3}", SunToday::L));
            row("T_eff [K]", "", format!("{:.0}", g.teff), format!("{:.0}", SunToday::TEFF));
            row("log g", "", format!("{:.3}", g.log_g), "4.438".into());
            row("T_c [K]", "", format!("{:.3e}", g.tc), format!("{:.3e}", SunToday::TC));
            row("ρ_c [g cm⁻³]", "", crate::views::fmt_value(g.rhoc), format!("{:.0}", SunToday::RHOC));
            row("P_c [dyn cm⁻²]", "", format!("{:.3e}", g.pc), "2.4e17".into());
            row(
                "Conv. core [m/M]",
                "",
                if g.conv_core_q > 0.0 { format!("{:.4}", g.conv_core_q) } else { "none".into() },
                "none".into(),
            );
            row(
                "Conv. env. base [r/R]",
                "",
                g.conv_env_base_r.map(|r| format!("{r:.4}")).unwrap_or_else(|| "none".into()),
                format!("{:.3}", SunToday::ENV_BASE),
            );
            row("pp / CNO [%]", "", format!("{:.1} / {:.1}", 100.0 * g.frac_pp, 100.0 * g.frac_cno), "99 / 1".into());
            row("Age", "", stellar_core::fmt_age(g.age_yr), "4.57 Gyr".into());
            row("Phase", "", short_phase(g.phase).into(), "main seq.".into());
            row("X_c / He core [m/M]", "", format!("{:.4} / {:.3}", g.x_c, g.he_core_q), "0.34 / 0".into());
            let period = model
                .pulsation
                .as_ref()
                .map(|p| crate::evolution::fmt_period(p.period_days()))
                .unwrap_or_else(|| "–".into());
            row(
                "Fund. radial mode",
                "calculated: period of the linear adiabatic fundamental radial mode (no driving or damping)",
                period,
                "≈ 1.1 h".into(),
            );
            row(
                "ν_max / Δν [μHz]",
                "scaling relations from g, T_eff and mean density: estimates, not a computed oscillation spectrum",
                format!("{:.0} / {:.1}", g.nu_max_uhz, g.delta_nu_uhz),
                "3090 / 135".into(),
            );
        });
        if let Some(kind) = stellar_core::pulsation::classify(g.teff, g.luminosity) {
            ui.label(
                RichText::new(format!("In the empirical instability strip: {} (by HR position)", kind.label()))
                    .small()
                    .color(theme::secondary(dark)),
            )
            .on_hover_text(
                "Calculated: the fundamental adiabatic radial period.\n\
                 Empirical inference: membership of the observed instability strip, from T_eff and L only.\n\
                 Not modelled: driving and damping (growth rates), the excitation mechanism, amplitudes, \
                 overtones and non-radial modes.",
            );
        }
        if matches!(
            g.phase,
            stellar_core::Phase::ShellHydrogen | stellar_core::Phase::RedGiant | stellar_core::Phase::CoreHelium
        ) {
            ui.label(
                RichText::new("After the main sequence: qualitative only (see About › Accuracy & limitations)")
                    .small()
                    .color(theme::secondary(dark)),
            )
            .on_hover_text(POST_MS_NOTE);
        }
        let r = &g.residuals;
        ui.label(
            RichText::new(format!(
                "Checks: virial {:.5}, energy {:.5}; {} Newton iterations, max residual {:.0e}",
                g.virial_ratio,
                g.energy_ratio,
                g.newton_iterations,
                r.max()
            ))
            .small()
            .color(theme::muted(dark)),
        )
        .on_hover_text(format!(
            "Numerical consistency of the solution, not agreement with real stars.\n\
             virial: global hydrostatic consistency (≈ 1)\n\
             energy: discretisation check; holds by construction for a converged model\n\
             largest residual of each discretised equation (dimensionless):\n\
             continuity {:.1e}, hydrostatic {:.1e}, energy {:.1e}, transport {:.1e}, boundary {:.1e}",
            r.continuity, r.hydrostatic, r.energy, r.transport, r.boundary
        ));

        ui.add_space(6.0);
        ui.horizontal_wrapped(|ui| {
            if ui
                .button("Pin current model")
                .on_hover_text("keep this model for comparison (Compare tab, overlays, HR diagram)")
                .clicked()
            {
                self.pin_current();
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("Export:");
            let stem = format!("star_{:.2}Msun", model.params.mass);
            if ui.button("Profile CSV").clicked() {
                self.notice = save_text(&format!("{stem}.csv"), "CSV", "csv", &profile_csv(&model));
            }
            if ui.button("Summary JSON").clicked() {
                self.notice = save_text(&format!("{stem}.json"), "JSON", "json", &summary_json(&model));
            }
            if let Some(run) = &self.evolution
                && ui.button("Track CSV").on_hover_text("global properties of every model on the timeline").clicked()
            {
                let refs: Vec<&StarModel> = run.snapshots.iter().map(|m| m.as_ref()).collect();
                self.notice =
                    save_text(&format!("{stem}_track.csv"), "CSV", "csv", &stellar_core::output::track_csv(&refs));
            }
            if ui.button("Screenshot PNG").clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .set_file_name(format!("{stem}_{}.png", self.tab.slug()))
                    .add_filter("PNG", &["png"])
                    .save_file()
            {
                let saved = Arc::new(AtomicBool::new(false));
                request_png(ui.ctx(), path.clone(), saved);
                self.notice = Some(format!("Screenshot saved to {}", path.display()));
            }
        });
        if let Some(n) = &self.notice {
            ui.label(RichText::new(n).small().color(theme::secondary(dark)));
        }
    }

    fn status_bar(&self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if let Some(run) = self.evolution.as_ref().filter(|r| r.running()) {
                ui.spinner();
                let g = &run.latest().globals;
                ui.label(format!(
                    "Evolving: {} steps, age {}, {}, L = {} L☉, T_eff = {:.0} K",
                    run.track.len() - 1,
                    stellar_core::fmt_age(g.age_yr),
                    g.phase.label(),
                    crate::views::fmt_value(g.luminosity),
                    g.teff
                ));
            } else if self.worker.busy() {
                ui.spinner();
                let it = self
                    .iteration
                    .map(|(i, c)| format!("  · Newton iteration {i}, max|Δy| = {c:.2e}"))
                    .unwrap_or_default();
                ui.label(format!("Solving: {}{it}", self.stage));
            } else if let Some(e) = &self.last_error {
                ui.colored_label(Color32::from_rgb(0xd0, 0x3b, 0x3b), "Error:");
                ui.label(e);
            } else if let (Some(m), Some(t)) = (&self.model, self.last_elapsed) {
                ui.label(format!(
                    "Converged: {} mesh points, {} Newton iterations in the last step, {:.0} ms",
                    m.profile.len(),
                    m.globals.newton_iterations,
                    t.as_secs_f64() * 1e3
                ));
            }
        });
    }

    fn top_bar(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            for tab in Tab::ALL {
                ui.selectable_value(&mut self.tab, tab, tab.title());
            }
            ui.separator();
            match self.tab {
                Tab::Profiles => {
                    ui.label("x axis:");
                    for axis in [XAxis::Radius, XAxis::Mass] {
                        ui.radio_value(&mut self.profile_axis, axis, axis.name());
                    }
                }
                Tab::Extra | Tab::Compare => {
                    ui.label("x axis:");
                    for axis in XAxis::ALL {
                        ui.radio_value(&mut self.detail_axis, axis, axis.name()).on_hover_text(match axis {
                            XAxis::Pressure => "log P spreads out the thin outer envelope, where opacity, gradients and convection change fastest",
                            _ => "",
                        });
                    }
                }
                _ => {}
            }
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                egui::widgets::global_theme_preference_buttons(ui);
            });
        });
    }

    fn central(&mut self, ui: &mut egui::Ui) {
        if self.tab == Tab::About {
            crate::views::about::show(ui, &mut self.about, &self.params, self.worker.physics());
            return;
        }
        if let Some(run) = &mut self.evolution {
            crate::evolution::timeline(ui, run, &mut self.scale_mode);
            ui.separator();
        }
        // The switch lives in the left panel and in the timeline; either may change it.
        let scale_changed = self.scale_mode != self.last_scale_mode;
        self.last_scale_mode = self.scale_mode;
        let Some(model) = self.shown_model() else {
            ui.centered_and_justified(|ui| {
                if let Some(e) = &self.last_error {
                    ui.label(format!("Could not compute a model: {e}"));
                } else {
                    ui.spinner();
                }
            });
            return;
        };
        let x_axis = if self.tab == Tab::Profiles { self.profile_axis } else { self.detail_axis };
        let track = self.evolution.as_ref().map(|r| r.track.as_slice()).unwrap_or(&[]);
        let fixed_mode = self.scale_mode == ScaleMode::FixedForRun;
        let fixed = self.evolution.as_ref().filter(|_| fixed_mode).map(|r| &r.ranges);
        let v = ViewCtx {
            model: &model,
            pinned: &self.pinned,
            x_axis,
            dark: ui.visuals().dark_mode,
            track,
            fixed,
            run: self.evolution.as_ref().map(|r| &r.ranges),
            zoom_track: self.evolution.is_some() && !fixed_mode,
            reset_plots: scale_changed,
        };
        match self.tab {
            Tab::Profiles => crate::views::profiles::show(ui, &v),
            Tab::Extra => crate::views::extra::show(ui, &v),
            Tab::CrossSection => crate::views::cross_section::show(ui, &v, &mut self.cross),
            Tab::Model3d => crate::views::model3d::show(ui, &v, &mut self.model3d),
            Tab::Hr => crate::views::hr::show(ui, &v, &self.zams),
            Tab::History => {
                if let Some(action) = crate::views::history::show(ui, &v, &mut self.history) {
                    self.history_action(action);
                }
            }
            Tab::About => {} // handled above
            Tab::Compare => {
                if let Some(action) = crate::views::compare::show(ui, &v, &mut self.compare) {
                    match action {
                        CompareAction::ToggleVisible(i) => self.pinned[i].visible ^= true,
                        CompareAction::Unpin(i) => {
                            self.pinned.remove(i);
                        }
                        CompareAction::UnpinAll => self.pinned.clear(),
                    }
                }
            }
        }
    }

    fn drive_capture(&mut self, ctx: &egui::Context, now: f64) {
        let Some(cap) = &mut self.capture else { return };
        ctx.request_repaint();
        if cap.started == 0.0 {
            cap.started = now;
        }
        let mut ready = self.model.is_some() && (self.zams.done || now - cap.started > 20.0);
        if ready && let Some((stop, age, select)) = self.capture_evolve.clone() {
            if !self.capture_evolve_started {
                self.capture_evolve_started = true;
                self.evo_settings.stop = stop;
                self.evo_settings.target_age_gyr = age;
                self.start_evolution();
                return;
            }
            let finished = self.evolution.as_ref().is_some_and(|r| !r.running());
            ready = finished || now - cap.started > 180.0;
            if ready && let Some(run) = &mut self.evolution {
                run.follow = false;
                run.playing = false;
                let n = run.snapshots.len();
                run.selected = match select.split_once(':') {
                    Some(("frac", f)) => {
                        ((f.parse::<f64>().unwrap_or(1.0) * (n - 1) as f64).round() as usize).min(n - 1)
                    }
                    _ if select == "pulsating" => run
                        .snapshots
                        .iter()
                        .position(|m| {
                            stellar_core::pulsation::classify(m.globals.teff, m.globals.luminosity)
                                == Some(stellar_core::pulsation::PulsatorType::Cepheid)
                        })
                        .or_else(|| {
                            run.snapshots.iter().position(|m| {
                                stellar_core::pulsation::classify(m.globals.teff, m.globals.luminosity).is_some()
                            })
                        })
                        .unwrap_or(n - 1),
                    _ => n - 1,
                };
            }
        }
        if !ready {
            return;
        }
        if cap.step >= cap.tabs.len() {
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            return;
        }
        let tab = cap.tabs[cap.step];
        self.tab = tab;
        cap.frames += 1;
        if cap.frames == 6 && !cap.requested {
            cap.requested = true;
            let path = cap.dir.join(format!("{}.png", tab.slug()));
            request_png(ctx, path, cap.saved.clone());
        }
        if cap.saved.swap(false, Ordering::SeqCst) {
            cap.step += 1;
            cap.frames = 0;
            cap.requested = false;
        }
    }
}

fn save_text(name: &str, kind: &str, ext: &str, content: &str) -> Option<String> {
    let path = rfd::FileDialog::new().set_file_name(name).add_filter(kind, &[ext]).save_file()?;
    Some(match std::fs::write(&path, content) {
        Ok(()) => format!("Saved {}", path.display()),
        Err(e) => format!("Could not save {}: {e}", path.display()),
    })
}

/// Ask egui for a screenshot of the next frame and write it as PNG.
fn request_png(ctx: &egui::Context, path: PathBuf, saved: Arc<AtomicBool>) {
    let ctx2 = ctx.clone();
    ctx.request_screenshot(move |image| {
        if let Err(e) = write_png(&path, &image) {
            eprintln!("screenshot failed: {e}");
        }
        saved.store(true, Ordering::SeqCst);
        ctx2.request_repaint();
    });
}

fn write_png(path: &std::path::Path, image: &egui::ColorImage) -> Result<(), Box<dyn std::error::Error>> {
    let [w, h] = image.size;
    let file = std::io::BufWriter::new(std::fs::File::create(path)?);
    let mut enc = png::Encoder::new(file, w as u32, h as u32);
    enc.set_color(png::ColorType::Rgba);
    enc.set_depth(png::BitDepth::Eight);
    let mut writer = enc.write_header()?;
    let bytes: Vec<u8> = image.pixels.iter().flat_map(|c| c.to_srgba_unmultiplied()).collect();
    writer.write_image_data(&bytes)?;
    Ok(())
}

impl eframe::App for StarApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        let now = ui.input(|i| i.time);
        self.handle_messages();
        self.maybe_auto_solve(now);
        let ctx = ui.ctx().clone();
        self.drive_capture(&ctx, now);
        if self.worker.busy() || (self.auto_solve && self.changed_at >= 0.0 && now - self.changed_at < 0.3) {
            ctx.request_repaint_after(Duration::from_millis(50));
        }
        if let Some(run) = &mut self.evolution {
            run.tick(ui.input(|i| i.stable_dt).min(0.1));
            if run.playing {
                ctx.request_repaint();
            }
        }

        egui::Panel::top("top").show(ui, |ui| {
            ui.add_space(3.0);
            self.top_bar(ui);
            ui.add_space(3.0);
        });
        egui::Panel::bottom("status").show(ui, |ui| self.status_bar(ui));
        egui::Panel::left("params").resizable(true).default_size(360.0).min_size(300.0).show(ui, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| self.params_panel(ui));
        });
        egui::CentralPanel::default().show(ui, |ui| self.central(ui));
    }
}

/// Compact phase name for the results table.
/// Why post-main-sequence models are only qualitative.
pub const POST_MS_NOTE: &str = "Not modelled: electron degeneracy, H/He ionization, conductive opacity, neutrino losses, \
     strong screening and molecular opacities. The mass mesh is fixed, so thin burning shells are not resolved.";

fn short_phase(p: stellar_core::Phase) -> &'static str {
    use stellar_core::Phase::*;
    match p {
        Zams => "ZAMS",
        MainSequence => "main seq.",
        ShellHydrogen => "H shell",
        RedGiant => "red giant",
        CoreHelium => "core He",
    }
}
