//! Evolution run state in the app: the evolutionary track (every step), a
//! thinned list of full models for the timeline, playback, and the controls.
//!
//! Runs can take thousands of steps. Every step adds a small track point
//! (for the HR diagram). A full model (~200 kB) is kept only when it differs
//! noticeably from the previous kept one, so memory stays bounded while the
//! timeline still resolves every phase.

use std::sync::Arc;

use egui::RichText;
use stellar_core::pulsation::{self, PulsatorType};
use stellar_core::{EvolutionOptions, Globals, Phase, StarModel, StopCondition, fmt_age};

use crate::theme;

/// One point of the evolutionary track.
#[derive(Clone, Copy, Debug)]
pub struct TrackPoint {
    pub age_yr: f64,
    pub teff: f64,
    pub lum: f64,
    pub radius: f64,
    pub x_c: f64,
    pub phase: Phase,
    pub period_days: Option<f64>,
    pub pulsator: Option<PulsatorType>,
}

impl TrackPoint {
    pub fn of(m: &StarModel) -> Self {
        let g = &m.globals;
        Self {
            age_yr: g.age_yr,
            teff: g.teff,
            lum: g.luminosity,
            radius: g.radius,
            x_c: g.x_c,
            phase: g.phase,
            period_days: m.pulsation.as_ref().map(|p| p.period_days()),
            pulsator: pulsation::classify(g.teff, g.luminosity),
        }
    }
}

/// Is `b` different enough from `a` to keep both as timeline snapshots?
fn significant(a: &Globals, b: &Globals) -> bool {
    let dlog = |x: f64, y: f64| (x / y).log10().abs();
    dlog(a.teff, b.teff) > 0.006
        || dlog(a.luminosity, b.luminosity) > 0.015
        || dlog(a.radius, b.radius) > 0.015
        || (a.x_c - b.x_c).abs() > 0.02
        || a.phase != b.phase
        || pulsation::classify(a.teff, a.luminosity) != pulsation::classify(b.teff, b.luminosity)
}

pub struct EvolutionRun {
    /// Worker job id while running.
    pub job: Option<u64>,
    pub track: Vec<TrackPoint>,
    /// Thinned full models; the last element is always the latest model.
    pub snapshots: Vec<Arc<StarModel>>,
    pub selected: usize,
    pub follow: bool,
    pub playing: bool,
    /// Playback speed [snapshots per second].
    pub play_rate: f32,
    play_accum: f32,
    pub finished: Option<String>,
    /// (next Δt [yr], ZAMS radius [R☉]) to resume from the latest model.
    pub resume: Option<(f64, f64)>,
}

impl EvolutionRun {
    pub fn new(zams: Arc<StarModel>) -> Self {
        Self {
            job: None,
            track: vec![TrackPoint::of(&zams)],
            snapshots: vec![zams],
            selected: 0,
            follow: true,
            playing: false,
            play_rate: 12.0,
            play_accum: 0.0,
            finished: None,
            resume: None,
        }
    }

    pub fn running(&self) -> bool {
        self.job.is_some()
    }

    pub fn latest(&self) -> &Arc<StarModel> {
        self.snapshots.last().expect("a run always has its ZAMS model")
    }

    pub fn selected_model(&self) -> &Arc<StarModel> {
        &self.snapshots[self.selected.min(self.snapshots.len() - 1)]
    }

    /// Add a newly computed model.
    pub fn push(&mut self, model: StarModel) {
        self.track.push(TrackPoint::of(&model));
        let n = self.snapshots.len();
        let model = Arc::new(model);
        // The current last element is a candidate: keep it only if it differs
        // enough from the one before it, otherwise the new model replaces it.
        if n >= 2 && !significant(&self.snapshots[n - 2].globals, &self.snapshots[n - 1].globals) {
            self.snapshots[n - 1] = model;
        } else {
            self.snapshots.push(model);
        }
        if self.follow {
            self.selected = self.snapshots.len() - 1;
        }
    }

    /// Advance playback by `dt` seconds of wall time.
    pub fn tick(&mut self, dt: f32) {
        if !self.playing {
            return;
        }
        self.play_accum += dt * self.play_rate;
        while self.play_accum >= 1.0 {
            self.play_accum -= 1.0;
            if self.selected + 1 < self.snapshots.len() {
                self.selected += 1;
            } else if !self.running() {
                self.playing = false;
                break;
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StopChoice {
    EndOfMainSequence,
    Age,
    AsFarAsPossible,
}

pub struct EvolutionSettings {
    pub stop: StopChoice,
    pub target_age_gyr: f64,
}

impl Default for EvolutionSettings {
    fn default() -> Self {
        Self { stop: StopChoice::EndOfMainSequence, target_age_gyr: 4.57 }
    }
}

impl EvolutionSettings {
    pub fn options(&self) -> EvolutionOptions {
        let stop = match self.stop {
            StopChoice::EndOfMainSequence => StopCondition::EndOfMainSequence,
            StopChoice::Age => StopCondition::Age(self.target_age_gyr * 1e9),
            StopChoice::AsFarAsPossible => StopCondition::AsFarAsPossible,
        };
        EvolutionOptions { stop, ..Default::default() }
    }
}

/// Requests from the evolution controls.
pub enum EvoAction {
    Start,
    Stop,
    Reset,
    SunAge,
}

/// Evolution section of the left panel.
pub fn controls(ui: &mut egui::Ui, settings: &mut EvolutionSettings, run: Option<&EvolutionRun>) -> Option<EvoAction> {
    let dark = ui.visuals().dark_mode;
    let mut action = None;
    ui.label(RichText::new("Evolution").strong());
    ui.horizontal_wrapped(|ui| {
        ui.label("Stop at:");
        ui.radio_value(&mut settings.stop, StopChoice::EndOfMainSequence, "end of main sequence")
            .on_hover_text("stop when the central hydrogen is exhausted (X_c < 0.001)");
        ui.radio_value(&mut settings.stop, StopChoice::AsFarAsPossible, "as far as possible")
            .on_hover_text("continue after the main sequence until the analytic physics or the solver gives out");
    });
    ui.horizontal(|ui| {
        ui.radio_value(&mut settings.stop, StopChoice::Age, "age");
        ui.add(
            egui::DragValue::new(&mut settings.target_age_gyr)
                .range(0.001..=30.0)
                .speed(0.01)
                .max_decimals(3)
                .suffix(" Gyr"),
        );
    });
    ui.horizontal_wrapped(|ui| {
        let running = run.is_some_and(EvolutionRun::running);
        let can_continue = run.is_some_and(|r| r.resume.is_some() && !r.running());
        let label = if can_continue { "Continue" } else { "Evolve" };
        if ui
            .add_enabled(!running, egui::Button::new(label))
            .on_hover_text("evolve the star in time from the ZAMS (or continue the current run)")
            .clicked()
        {
            action = Some(EvoAction::Start);
        }
        if ui.add_enabled(running, egui::Button::new("Stop")).clicked() {
            action = Some(EvoAction::Stop);
        }
        if ui.add_enabled(run.is_some() && !running, egui::Button::new("Reset to ZAMS")).clicked() {
            action = Some(EvoAction::Reset);
        }
        if ui.button("Sun's age").on_hover_text("evolve to 4.57 Gyr").clicked() {
            action = Some(EvoAction::SunAge);
        }
    });
    if let Some(run) = run {
        let latest = &run.latest().globals;
        let text = if run.running() {
            format!("running: {} steps, age {}, {}", run.track.len() - 1, fmt_age(latest.age_yr), latest.phase.label())
        } else if let Some(reason) = &run.finished {
            format!("{} steps to {}: {reason}", run.track.len() - 1, fmt_age(latest.age_yr))
        } else {
            String::new()
        };
        ui.label(RichText::new(text).small().color(theme::secondary(dark)));
    }
    action
}

/// Timeline above the views: scrub, step, play through the computed ages.
pub fn timeline(ui: &mut egui::Ui, run: &mut EvolutionRun) {
    let dark = ui.visuals().dark_mode;
    let n = run.snapshots.len();
    ui.horizontal(|ui| {
        ui.label(RichText::new("Age").strong());
        if ui.small_button("<").on_hover_text("previous model").clicked() && run.selected > 0 {
            run.selected -= 1;
            run.follow = false;
        }
        let play = if run.playing { "Pause" } else { "Play" };
        if ui.small_button(play).on_hover_text("play through the computed models").clicked() {
            if !run.playing && run.selected + 1 >= n {
                run.selected = 0;
            }
            run.playing = !run.playing;
            run.follow = false;
        }
        if ui.small_button(">").on_hover_text("next model").clicked() && run.selected + 1 < n {
            run.selected += 1;
        }
        let ages: Vec<f64> = run.snapshots.iter().map(|m| m.globals.age_yr).collect();
        let mut idx = run.selected.min(n - 1);
        let width = (ui.available_width() - 150.0).clamp(160.0, 1200.0);
        ui.spacing_mut().slider_width = width;
        let resp = ui.add(
            egui::Slider::new(&mut idx, 0..=n.saturating_sub(1))
                .show_value(false)
                .custom_formatter(move |v, _| fmt_age(ages[(v as usize).min(ages.len() - 1)])),
        );
        if resp.changed() {
            run.selected = idx;
            run.follow = idx + 1 == n;
            run.playing = false;
        }
        ui.checkbox(&mut run.follow, "follow latest");
        if run.follow {
            run.selected = n - 1;
        }
    });
    ui.horizontal_wrapped(|ui| {
        let m = run.selected_model();
        let g = &m.globals;
        let mut text = format!("{} · {} · X_c = {:.3}", fmt_age(g.age_yr), g.phase.label(), g.x_c);
        if let Some(kind) = pulsation::classify(g.teff, g.luminosity) {
            let period =
                m.pulsation.as_ref().map(|p| format!(", P = {}", fmt_period(p.period_days()))).unwrap_or_default();
            text.push_str(&format!(" · {}{period}", kind.label()));
        }
        ui.label(RichText::new(text).color(theme::secondary(dark)));
    });
}

/// Period with a sensible unit.
pub fn fmt_period(days: f64) -> String {
    if days >= 1.0 {
        format!("{days:.2} d")
    } else if days * 24.0 >= 1.0 {
        format!("{:.2} h", days * 24.0)
    } else {
        format!("{:.1} min", days * 1440.0)
    }
}
