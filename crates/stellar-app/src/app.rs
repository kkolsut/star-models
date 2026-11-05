//! Application state, side panel (parameters + results), tabs and status bar.

use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Duration;

use egui::{Color32, RichText};
use stellar_core::output::{profile_csv, summary_json};
use stellar_core::{NoControl, Physics, SolveError, SolveEvent, StarModel, StarParams, solve_model};

use crate::theme;
use crate::views::compare::{CompareAction, CompareState};
use crate::views::cross_section::CrossSectionState;
use crate::views::hr::ZamsLine;
use crate::views::{PinnedModel, ViewCtx, XAxis, model_label, swatch};
use crate::worker::{Msg, Worker, ZamsKey};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Tab {
    Profiles,
    Extra,
    CrossSection,
    Hr,
    Compare,
}

impl Tab {
    const ALL: [Tab; 5] = [Tab::Profiles, Tab::Extra, Tab::CrossSection, Tab::Hr, Tab::Compare];

    fn title(self) -> &'static str {
        match self {
            Tab::Profiles => "Profiles (P, T, ρ)",
            Tab::Extra => "More profiles",
            Tab::CrossSection => "Cross-section",
            Tab::Hr => "HR diagram",
            Tab::Compare => "Compare",
        }
    }

    fn slug(self) -> &'static str {
        match self {
            Tab::Profiles => "profiles",
            Tab::Extra => "extra",
            Tab::CrossSection => "cross_section",
            Tab::Hr => "hr",
            Tab::Compare => "compare",
        }
    }
}

/// Automated screenshot mode (STAR_MODELS_CAPTURE=<dir>): visits every tab,
/// saves a PNG of each, then quits. Used to check the rendering.
struct Capture {
    dir: PathBuf,
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
    compare: CompareState,
    zams: ZamsLine,
    notice: Option<String>,
    capture: Option<Capture>,
}

impl StarApp {
    pub fn new(cc: &eframe::CreationContext<'_>) -> Self {
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

        let capture = std::env::var("STAR_MODELS_CAPTURE").ok().map(|d| Capture {
            dir: PathBuf::from(d),
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
            compare: CompareState::default(),
            zams: ZamsLine::default(),
            notice: None,
            capture,
        };
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
        let Some(model) = &self.model else { return };
        let used: Vec<usize> = self.pinned.iter().map(|p| p.slot).collect();
        let Some(slot) = (1..theme::MAX_SLOTS).find(|s| !used.contains(s)) else {
            self.notice = Some("At most 7 pinned models (one per palette color). Unpin one first.".into());
            return;
        };
        self.pinned.push(PinnedModel { model: model.clone(), slot, label: model_label(model), visible: true });
    }

    // ------------------------------------------------------------------ panels

    fn params_panel(&mut self, ui: &mut egui::Ui) {
        let dark = ui.visuals().dark_mode;
        ui.add_space(4.0);
        ui.heading("Star Models");
        ui.label(RichText::new("Zero-age main-sequence structure").color(theme::secondary(dark)));
        ui.label(
            RichText::new("Physics: analytic (Kramers + H⁻ opacity, ideal gas + radiation, MLT convection)")
                .small()
                .color(theme::muted(dark)),
        )
        .on_hover_text(self.worker.physics().describe());
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
        self.results(ui);
    }

    fn results(&mut self, ui: &mut egui::Ui) {
        let dark = ui.visuals().dark_mode;
        let Some(model) = self.model.clone() else {
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
            let mut row = |name: &str, val: String, sun: String| {
                ui.label(name);
                ui.monospace(val);
                ui.label(RichText::new(sun).monospace().color(theme::muted(dark)));
                ui.end_row();
            };
            row("R [R☉]", format!("{:.4}", g.radius), format!("{:.3}", SunToday::R));
            row("L [L☉]", crate::views::fmt_value(g.luminosity), format!("{:.3}", SunToday::L));
            row("T_eff [K]", format!("{:.0}", g.teff), format!("{:.0}", SunToday::TEFF));
            row("log g", format!("{:.3}", g.log_g), "4.438".into());
            row("T_c [K]", format!("{:.3e}", g.tc), format!("{:.3e}", SunToday::TC));
            row("ρ_c [g cm⁻³]", crate::views::fmt_value(g.rhoc), format!("{:.0}", SunToday::RHOC));
            row("P_c [dyn cm⁻²]", format!("{:.3e}", g.pc), "2.4e17".into());
            row(
                "Conv. core [m/M]",
                if g.conv_core_q > 0.0 { format!("{:.4}", g.conv_core_q) } else { "none".into() },
                "none".into(),
            );
            row(
                "Conv. env. base [r/R]",
                g.conv_env_base_r.map(|r| format!("{r:.4}")).unwrap_or_else(|| "none".into()),
                format!("{:.3}", SunToday::ENV_BASE),
            );
            row("pp / CNO [%]", format!("{:.1} / {:.1}", 100.0 * g.frac_pp, 100.0 * g.frac_cno), "99 / 1".into());
        });
        ui.label(
            RichText::new(format!(
                "Checks: virial {:.5}, energy {:.5}; {} Newton iterations",
                g.virial_ratio, g.energy_ratio, g.newton_iterations
            ))
            .small()
            .color(theme::muted(dark)),
        );

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
            if self.worker.busy() {
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
        let Some(model) = self.model.clone() else {
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
        let v = ViewCtx { model: &model, pinned: &self.pinned, x_axis, dark: ui.visuals().dark_mode };
        match self.tab {
            Tab::Profiles => crate::views::profiles::show(ui, &v),
            Tab::Extra => crate::views::extra::show(ui, &v),
            Tab::CrossSection => crate::views::cross_section::show(ui, &v, &mut self.cross),
            Tab::Hr => crate::views::hr::show(ui, &v, &self.zams),
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
        let ready = self.model.is_some() && (self.zams.done || now - cap.started > 20.0);
        if !ready {
            return;
        }
        if cap.step >= Tab::ALL.len() {
            ctx.send_viewport_cmd(egui::ViewportCommand::Close);
            return;
        }
        let tab = Tab::ALL[cap.step];
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
