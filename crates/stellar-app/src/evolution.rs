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
    pub log_g: f64,
    pub tc: f64,
    pub rhoc: f64,
    pub x_c: f64,
    pub y_c: f64,
    pub he_core_q: f64,
    pub conv_core_q: f64,
    pub conv_env_base_r: Option<f64>,
    /// Shares of the surface luminosity from pp, CNO, 3α and gravitational energy.
    pub l_pp: f64,
    pub l_cno: f64,
    pub l_3a: f64,
    pub l_grav: f64,
    pub nu_max: f64,
    pub delta_nu: f64,
    pub phase: Phase,
    pub period_days: Option<f64>,
    pub pulsator: Option<PulsatorType>,
}

impl TrackPoint {
    pub fn of(m: &StarModel) -> Self {
        let g = &m.globals;
        // frac_pp/cno/3a are shares of the nuclear luminosity; L_nuc / L = 1 − L_grav / L.
        let nuc = 1.0 - g.frac_grav;
        Self {
            age_yr: g.age_yr,
            teff: g.teff,
            lum: g.luminosity,
            radius: g.radius,
            log_g: g.log_g,
            tc: g.tc,
            rhoc: g.rhoc,
            x_c: g.x_c,
            y_c: g.y_c,
            he_core_q: g.he_core_q,
            conv_core_q: g.conv_core_q,
            conv_env_base_r: g.conv_env_base_r,
            l_pp: g.frac_pp * nuc,
            l_cno: g.frac_cno * nuc,
            l_3a: g.frac_3a * nuc,
            l_grav: g.frac_grav,
            nu_max: g.nu_max_uhz,
            delta_nu: g.delta_nu_uhz,
            phase: g.phase,
            period_days: m.pulsation.as_ref().map(|p| p.period_days()),
            pulsator: pulsation::classify(g.teff, g.luminosity),
        }
    }
}

/// Every step of a run as CSV (the full history, not only the timeline models).
pub fn history_csv(track: &[TrackPoint]) -> String {
    use std::fmt::Write as _;
    let mut s = String::from(
        "step,age_yr,phase,L_lsun,R_rsun,Teff_K,log_g,Tc_K,rhoc_g_cm3,Xc,Yc,he_core_q,conv_core_q,conv_env_base_r,\
         L_pp_frac,L_cno_frac,L_3a_frac,L_grav_frac,P_fund_days,pulsator,nu_max_uHz,delta_nu_uHz\n",
    );
    for (i, t) in track.iter().enumerate() {
        let _ = writeln!(
            s,
            "{i},{:.6e},{:?},{:.6e},{:.6e},{:.2},{:.4},{:.6e},{:.6e},{:.6},{:.6},{:.6},{:.6},{},{:.5},{:.5},{:.5},{:.5},{},{},{:.4},{:.4}",
            t.age_yr,
            t.phase,
            t.lum,
            t.radius,
            t.teff,
            t.log_g,
            t.tc,
            t.rhoc,
            t.x_c,
            t.y_c,
            t.he_core_q,
            t.conv_core_q,
            t.conv_env_base_r.map(|r| format!("{r:.5}")).unwrap_or_default(),
            t.l_pp,
            t.l_cno,
            t.l_3a,
            t.l_grav,
            t.period_days.map(|p| format!("{p:.6e}")).unwrap_or_default(),
            t.pulsator.map(|p| p.label()).unwrap_or(""),
            t.nu_max,
            t.delta_nu,
        );
    }
    s
}

/// How plots, colour scales and the star's size respond as the shown model changes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScaleMode {
    /// Ranges cover every model of the run: axes and colours stay still, the star is drawn to scale.
    FixedForRun,
    /// Each view fits the model on screen (best detail; axes move during playback).
    AdaptToModel,
}

/// A value range [lo, hi]; empty while lo > hi.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Range {
    pub lo: f64,
    pub hi: f64,
}

impl Range {
    pub const EMPTY: Range = Range { lo: f64::INFINITY, hi: f64::NEG_INFINITY };
    pub const UNIT: Range = Range { lo: 0.0, hi: 1.0 };

    pub fn include(&mut self, v: f64) {
        if v.is_finite() {
            self.lo = self.lo.min(v);
            self.hi = self.hi.max(v);
        }
    }

    pub fn is_valid(&self) -> bool {
        self.lo <= self.hi
    }
}

/// Ranges covering every model of a run, in plotted units (log₁₀ for
/// logarithmic quantities). Updated incrementally as models arrive.
#[derive(Clone, Copy, Debug)]
pub struct RunRanges {
    pub log_p: Range,
    pub log_t: Range,
    pub log_rho: Range,
    pub log_kappa: Range,
    pub log_eps: Range,
    pub log_nabla: Range,
    pub log_vconv: Range,
    pub eps_grav: Range,
    /// Luminosity L(r) [erg s⁻¹] (linear).
    pub lum: Range,
    pub log_teff: Range,
    pub log_lum_sun: Range,
    /// Photospheric radius [R☉].
    pub radius: Range,
    /// Largest hydrogen mass fraction.
    pub x_max: f64,
}

impl RunRanges {
    pub fn of(m: &StarModel) -> Self {
        let mut r = Self {
            log_p: Range::EMPTY,
            log_t: Range::EMPTY,
            log_rho: Range::EMPTY,
            log_kappa: Range::EMPTY,
            log_eps: Range::EMPTY,
            log_nabla: Range::EMPTY,
            log_vconv: Range::EMPTY,
            eps_grav: Range::EMPTY,
            lum: Range::EMPTY,
            log_teff: Range::EMPTY,
            log_lum_sun: Range::EMPTY,
            radius: Range::EMPTY,
            x_max: 0.0,
        };
        r.include(m);
        r
    }

    pub fn include(&mut self, m: &StarModel) {
        let p = &m.profile;
        // Same floor as the energy-generation plot: 8 dex below the peak.
        let eps_floor = p.eps.iter().cloned().fold(0.0, f64::max) * 1e-8;
        let lg = |v: f64| if v > 0.0 { v.log10() } else { f64::NAN };
        for k in 0..p.len() {
            self.log_p.include(lg(p.p[k]));
            self.log_t.include(lg(p.t[k]));
            self.log_rho.include(lg(p.rho[k]));
            self.log_kappa.include(lg(p.kappa[k]));
            for e in [p.eps_pp[k], p.eps_cno[k], p.eps_3a[k], p.eps[k]] {
                if e > eps_floor {
                    self.log_eps.include(e.log10());
                }
            }
            for n in [p.nabla_rad[k], p.nabla_ad[k], p.nabla[k]] {
                self.log_nabla.include(lg(n));
            }
            if p.convective[k] {
                self.log_vconv.include(lg(p.v_conv[k]));
            }
            self.eps_grav.include(p.eps_grav[k]);
            self.lum.include(p.l[k]);
            self.x_max = self.x_max.max(p.x[k]);
        }
        self.eps_grav.include(0.0);
        let g = &m.globals;
        self.log_teff.include(g.teff.log10());
        self.log_lum_sun.include(g.luminosity.log10());
        self.radius.include(g.radius);
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
    /// Ranges over every model computed in this run.
    pub ranges: RunRanges,
}

impl EvolutionRun {
    pub fn new(zams: Arc<StarModel>) -> Self {
        Self {
            ranges: RunRanges::of(&zams),
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
        self.ranges.include(&model);
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
/// Display options for evolution runs, editable before and during a run.
pub struct DisplayPrefs<'a> {
    pub scale: &'a mut ScaleMode,
    /// "True size" of the Cross-section and 3D views (kept per view; set together here).
    pub true_size_cross: &'a mut bool,
    pub true_size_3d: &'a mut bool,
}

/// The two scale options as selectable buttons (shared by the left panel and the timeline).
fn scale_buttons(ui: &mut egui::Ui, scale: &mut ScaleMode) {
    ui.selectable_value(scale, ScaleMode::FixedForRun, "fixed for run").on_hover_text(
        "axes and colour scales cover every model of the run, so they stay still while the star changes; \
         the HR diagram shows the whole diagram",
    );
    ui.selectable_value(scale, ScaleMode::AdaptToModel, "adapt to model").on_hover_text(
        "every view fits the model on screen: most detail, but axes and colours change during playback; \
         the HR diagram zooms in on the track",
    );
}

pub fn controls(
    ui: &mut egui::Ui,
    settings: &mut EvolutionSettings,
    run: Option<&EvolutionRun>,
    prefs: DisplayPrefs<'_>,
) -> Option<EvoAction> {
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
    // Display options: can be set before a run starts and apply while it is shown.
    ui.horizontal_wrapped(|ui| {
        ui.label("Scale:");
        scale_buttons(ui, prefs.scale);
    });
    ui.horizontal(|ui| {
        let mixed = *prefs.true_size_cross != *prefs.true_size_3d;
        let mut both = *prefs.true_size_cross && *prefs.true_size_3d;
        let resp = ui.add(egui::Checkbox::new(&mut both, "true size").indeterminate(mixed)).on_hover_text(
            "Cross-section and 3D model: draw the star to scale relative to its largest radius in the run, \
                 so you see it expand (each view also has its own checkbox)",
        );
        if resp.changed() {
            *prefs.true_size_cross = both;
            *prefs.true_size_3d = both;
        }
        if run.is_none() {
            ui.label(RichText::new("(apply once a run exists)").small().color(theme::secondary(dark)));
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

/// Timeline above the views: scrub, step, play through the computed ages, and
/// choose how the views scale.
pub fn timeline(ui: &mut egui::Ui, run: &mut EvolutionRun, scale: &mut ScaleMode) {
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
    ui.horizontal(|ui| {
        ui.label("Scale:");
        scale_buttons(ui, scale);
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
