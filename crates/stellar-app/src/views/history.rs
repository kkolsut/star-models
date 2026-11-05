//! Time evolution of the star's global parameters over an evolution run:
//! linked time-series plots (luminosity, radius, temperatures, central
//! conditions, composition, core masses, energy sources, pulsation period,
//! asteroseismic frequencies). Every step of the run is plotted. Phase changes
//! are marked, the age shown in the other views is a vertical marker, and a
//! click on any plot selects that age.

use egui::{Align2, Id, RichText};
use egui_plot::{GridInput, GridMark, HoverPosition, Line, Plot, PlotPoint, PlotPoints, PlotUi, Span, Text, VLine};
use stellar_core::{Phase, fmt_age};

use super::{ViewCtx, fmt_value, log_spacer, log_tick};
use crate::evolution::{TrackPoint, fmt_period};
use crate::theme;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HistoryAxis {
    /// Age since the ZAMS (linear).
    Age,
    /// Model (step) number: spreads out the fast phases after the main sequence.
    Step,
}

pub struct HistoryState {
    pub axis: HistoryAxis,
    /// Shrink the panels so all of them fit on one screen.
    pub fit_all: bool,
}

impl Default for HistoryState {
    fn default() -> Self {
        Self { axis: HistoryAxis::Age, fit_all: false }
    }
}

/// Grid marks at 1-2-5 × 10ⁿ steps, about eight labelled marks across the range.
fn nice_spacer(input: GridInput) -> Vec<GridMark> {
    let (lo, hi) = input.bounds;
    let span = (hi - lo).max(1e-12);
    let raw = span / 8.0;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0].into_iter().map(|k| k * mag).find(|&s| s >= raw).unwrap_or(10.0 * mag);
    let minor = step / 5.0;
    let first = (lo / minor).floor() as i64;
    let last = (hi / minor).ceil() as i64;
    (first..=last)
        .map(|i| {
            let v = i as f64 * minor;
            let major = (v / step - (v / step).round()).abs() < 1e-6;
            GridMark { value: v, step_size: if major { step } else { minor } }
        })
        .filter(|m| m.value >= lo && m.value <= hi)
        .collect()
}

/// Labels for [`nice_spacer`] marks, with as many decimals as the step needs.
fn nice_label(mark: GridMark, _range: &std::ops::RangeInclusive<f64>) -> String {
    let step = mark.step_size;
    let decimals = if step >= 1.0 { 0 } else { (-step.log10().floor()) as usize };
    format!("{:.*}", decimals, mark.value)
}

pub enum HistoryAction {
    /// Show the model closest to this age [yr] in all views.
    SelectAge(f64),
    /// Start an evolution run (no run yet).
    Evolve,
    /// Export every step as CSV.
    Export,
}

/// A plotted series: name, categorical slot, value per track point.
struct Curve {
    name: &'static str,
    slot: usize,
    get: fn(&TrackPoint) -> Option<f64>,
}

/// Shared x data and annotations for all panels.
struct Frame<'a> {
    track: &'a [TrackPoint],
    xs: Vec<f64>,
    x_label: String,
    axis: HistoryAxis,
    /// Divisor turning years into the age-axis unit.
    unit: f64,
    selected_x: f64,
    /// (x, label) where the evolutionary phase changes.
    phase_marks: Vec<(f64, &'static str)>,
    link: Id,
    dark: bool,
    /// Full x range of the run (every panel shows all of it).
    x_range: (f64, f64),
}

impl Frame<'_> {
    fn nearest(&self, x: f64) -> usize {
        let i = self.xs.partition_point(|&v| v < x).min(self.xs.len() - 1);
        if i > 0 && (x - self.xs[i - 1]).abs() < (self.xs[i] - x).abs() { i - 1 } else { i }
    }

    fn age_of_x(&self, x: f64) -> f64 {
        match self.axis {
            HistoryAxis::Age => x * self.unit,
            HistoryAxis::Step => self.track[self.nearest(x)].age_yr,
        }
    }
}

/// Unit for an age axis covering up to `max_yr`.
fn age_unit(max_yr: f64) -> (f64, &'static str) {
    if max_yr >= 1e9 {
        (1e9, "Gyr")
    } else if max_yr >= 1e6 {
        (1e6, "Myr")
    } else if max_yr >= 1e3 {
        (1e3, "kyr")
    } else {
        (1.0, "yr")
    }
}

/// Split a curve into runs of valid values (positive when logarithmic).
fn segments(f: &Frame<'_>, get: fn(&TrackPoint) -> Option<f64>, log: bool) -> Vec<Vec<[f64; 2]>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    for (t, &x) in f.track.iter().zip(&f.xs) {
        match get(t).filter(|v| v.is_finite() && (!log || *v > 0.0)) {
            Some(v) => cur.push([x, if log { v.log10() } else { v }]),
            None if !cur.is_empty() => out.push(std::mem::take(&mut cur)),
            None => {}
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// One time-series panel. Returns the clicked x position, if any.
#[allow(clippy::too_many_arguments)]
fn panel(
    ui: &mut egui::Ui,
    f: &Frame<'_>,
    title: &str,
    y_label: &str,
    log: bool,
    curves: &[Curve],
    height: f32,
    first: bool,
    extra: &dyn Fn(&mut PlotUi<'_>),
    extra_legend: Option<&str>,
) -> Option<f64> {
    // Title with an inline legend: swatches never cover the data, even when the
    // panels are short.
    ui.horizontal_wrapped(|ui| {
        ui.label(RichText::new(title).strong());
        if curves.len() > 1 || extra_legend.is_some() {
            ui.add_space(6.0);
            for c in curves {
                super::swatch(ui, theme::series(c.slot, f.dark));
                ui.label(RichText::new(c.name).small());
            }
            if let Some(name) = extra_legend {
                super::swatch(ui, theme::convective_fill(f.dark).gamma_multiply(3.0));
                ui.label(RichText::new(name).small());
            }
        }
    });
    let mut plot = Plot::new(("history", title))
        .height(height)
        .x_axis_label(f.x_label.clone())
        .y_axis_label(y_label)
        .allow_scroll(false)
        .x_grid_spacer(nice_spacer)
        .x_axis_formatter(nice_label)
        .link_axis(f.link, [true, false])
        .link_cursor(f.link, [true, false])
        .label_formatter(move |pos: &HoverPosition<'_>| {
            let (name, p) = match pos {
                HoverPosition::NearDataPoint { plot_name, position, .. } => (*plot_name, *position),
                HoverPosition::Elsewhere { position } => ("", *position),
            };
            let t = &f.track[f.nearest(p.x)];
            let step = f.nearest(p.x);
            let head = format!("{} (model {step})\n{}", fmt_age(t.age_yr), t.phase.label());
            if name.is_empty() {
                return Some(head);
            }
            let y = if log { 10f64.powf(p.y) } else { p.y };
            Some(format!("{head}\n{name} = {}", fmt_value(y)))
        });
    plot = plot.include_x(f.x_range.0).include_x(f.x_range.1);
    plot = if log {
        plot.y_grid_spacer(log_spacer).y_axis_formatter(log_tick)
    } else {
        plot.y_axis_formatter(|mark, _| {
            let v = mark.value;
            let a = v.abs();
            if a == 0.0 || (1e-2..1e4).contains(&a) { format!("{v:.2}") } else { format!("{v:.1e}") }
        })
    };
    let muted = theme::muted(f.dark);
    let resp = plot.show(ui, |pu| {
        extra(pu);
        let top = pu.plot_bounds().max()[1];
        // Phase changes: a line each; labels only where they do not collide.
        let mut last_label_px = f32::NEG_INFINITY;
        for &(x, label) in &f.phase_marks {
            pu.add(VLine::new("", x).color(muted).width(1.0));
            let px = pu.screen_from_plot(PlotPoint::new(x, top)).x;
            if first && px - last_label_px > 130.0 {
                pu.add(
                    Text::new("", PlotPoint::new(x, top), format!(" {label}")).color(muted).anchor(Align2::LEFT_TOP),
                );
                last_label_px = px;
            }
        }
        for c in curves {
            let color = theme::series(c.slot, f.dark);
            for (i, seg) in segments(f, c.get, log).into_iter().enumerate() {
                pu.line(Line::new(c.name, PlotPoints::from(seg)).color(color).width(2.0).id((title, c.name, i)));
            }
        }
        // The age shown in the other views.
        pu.add(VLine::new("", f.selected_x).color(theme::primary(f.dark).gamma_multiply(0.7)).width(1.5));
        pu.pointer_coordinate()
    });
    if resp.response.clicked() { resp.inner.map(|p| p.x) } else { None }
}

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>, state: &mut HistoryState) -> Option<HistoryAction> {
    let dark = v.dark;
    if v.track.len() < 2 {
        let mut action = None;
        ui.add_space(40.0);
        ui.vertical_centered(|ui| {
            ui.label(RichText::new("No evolution run yet").size(18.0).strong());
            ui.add_space(6.0);
            ui.label(
                RichText::new(
                    "This tab shows how the star's global parameters change with age: luminosity, radius, temperatures, \
                     central conditions, composition, core masses, energy sources and pulsation period.",
                )
                .color(theme::secondary(dark)),
            );
            ui.add_space(10.0);
            if ui.button("Evolve (with the settings in the left panel)").clicked() {
                action = Some(HistoryAction::Evolve);
            }
        });
        return action;
    }

    let mut action = None;
    ui.horizontal(|ui| {
        ui.label("x axis:");
        ui.radio_value(&mut state.axis, HistoryAxis::Age, "age");
        ui.radio_value(&mut state.axis, HistoryAxis::Step, "model number").on_hover_text(
            "one unit per computed model: spreads out the fast phases after the main sequence, \
             which take only a few percent of the age",
        );
        ui.separator();
        ui.checkbox(&mut state.fit_all, "fit all on screen");
        ui.separator();
        if ui.button("Export history CSV").on_hover_text("every step of the run, all global quantities").clicked() {
            action = Some(HistoryAction::Export);
        }
        ui.label(
            RichText::new("click a plot to show that age in all views · vertical lines: phase changes")
                .small()
                .color(theme::muted(dark)),
        );
    });

    let track = v.track;
    let (unit, unit_name) = age_unit(track.last().map(|t| t.age_yr).unwrap_or(1.0));
    let xs: Vec<f64> = match state.axis {
        HistoryAxis::Age => track.iter().map(|t| t.age_yr / unit).collect(),
        HistoryAxis::Step => (0..track.len()).map(|i| i as f64).collect(),
    };
    let x_of_age = |age: f64| -> f64 {
        match state.axis {
            HistoryAxis::Age => age / unit,
            HistoryAxis::Step => {
                let i = track.partition_point(|t| t.age_yr < age).min(track.len() - 1);
                i as f64
            }
        }
    };
    let mut phase_marks = Vec::new();
    for w in 1..track.len() {
        let (a, b) = (track[w - 1].phase, track[w].phase);
        if a != b && a != Phase::Zams {
            let label = match b {
                Phase::ShellHydrogen => "end of main sequence",
                Phase::RedGiant => "red giant",
                Phase::CoreHelium => "core He burning",
                Phase::MainSequence => "main sequence",
                Phase::Zams => "",
            };
            phase_marks.push((xs[w], label));
        }
    }
    let f = Frame {
        track,
        selected_x: x_of_age(v.model.globals.age_yr),
        xs,
        x_label: match state.axis {
            HistoryAxis::Age => format!("age [{unit_name}]"),
            HistoryAxis::Step => "model number".into(),
        },
        axis: state.axis,
        unit,
        phase_marks,
        link: Id::unique("history_link"),
        dark,
        x_range: (0.0, 0.0),
    };
    let f = Frame { x_range: (f.xs[0], *f.xs.last().unwrap()), ..f };

    // Shaded intervals where the star lies in the instability strip.
    let strip_spans: Vec<(f64, f64)> = {
        let mut spans = Vec::new();
        let mut start: Option<f64> = None;
        for (t, &x) in track.iter().zip(&f.xs) {
            match (t.pulsator.is_some(), start) {
                (true, None) => start = Some(x),
                (false, Some(a)) => {
                    spans.push((a, x));
                    start = None;
                }
                _ => {}
            }
        }
        if let Some(a) = start {
            spans.push((a, *f.xs.last().unwrap()));
        }
        spans
    };
    let strip = |pu: &mut PlotUi<'_>| {
        for &(a, b) in &strip_spans {
            pu.add(Span::new("", a..=b).fill(theme::convective_fill(dark)).border_width(0.0));
        }
    };
    let none = |_: &mut PlotUi<'_>| {};

    // Three rows per screen (scroll for more), or all six rows at once.
    let rows = if state.fit_all { 6.0 } else { 3.0 };
    let h = ((ui.available_height() - rows * (ui.text_style_height(&egui::TextStyle::Body) + 10.0) - 30.0) / rows)
        .max(if state.fit_all { 90.0 } else { 150.0 });
    let mut clicked: Option<f64> = None;
    let mut take = |c: Option<f64>| {
        if c.is_some() {
            clicked = c;
        }
    };
    egui::ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        ui.columns(2, |cols| {
            take(panel(
                &mut cols[0],
                &f,
                "Luminosity",
                "L [L☉]",
                true,
                &[Curve { name: "L", slot: 0, get: |t| Some(t.lum) }],
                h,
                true,
                &none,
                None,
            ));
            take(panel(
                &mut cols[1],
                &f,
                "Radius",
                "R [R☉]",
                true,
                &[Curve { name: "R", slot: 0, get: |t| Some(t.radius) }],
                h,
                true,
                &none,
                None,
            ));
            take(panel(
                &mut cols[0],
                &f,
                "Effective temperature",
                "T_eff [K]",
                true,
                &[Curve { name: "T_eff", slot: 0, get: |t| Some(t.teff) }],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[1],
                &f,
                "Surface gravity",
                "log g [cgs]",
                false,
                &[Curve { name: "log g", slot: 0, get: |t| Some(t.log_g) }],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[0],
                &f,
                "Central temperature",
                "T_c [K]",
                true,
                &[Curve { name: "T_c", slot: 0, get: |t| Some(t.tc) }],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[1],
                &f,
                "Central density",
                "ρ_c [g cm⁻³]",
                true,
                &[Curve { name: "ρ_c", slot: 0, get: |t| Some(t.rhoc) }],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[0],
                &f,
                "Central composition",
                "mass fraction",
                false,
                &[
                    Curve { name: "X_c (hydrogen)", slot: 0, get: |t| Some(t.x_c) },
                    Curve { name: "Y_c (helium)", slot: 1, get: |t| Some(t.y_c) },
                ],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[1],
                &f,
                "Core masses",
                "m / M",
                false,
                &[
                    Curve { name: "helium core (X < 0.001)", slot: 0, get: |t| Some(t.he_core_q) },
                    Curve { name: "convective core", slot: 1, get: |t| Some(t.conv_core_q) },
                ],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[0],
                &f,
                "Energy sources (share of the surface luminosity)",
                "fraction of L",
                false,
                &[
                    Curve { name: "pp chain", slot: 0, get: |t| Some(t.l_pp) },
                    Curve { name: "CNO cycle", slot: 1, get: |t| Some(t.l_cno) },
                    Curve { name: "triple-α", slot: 2, get: |t| Some(t.l_3a) },
                    Curve { name: "gravitational (contraction > 0)", slot: 3, get: |t| Some(t.l_grav) },
                ],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[1],
                &f,
                "Fundamental radial pulsation period",
                "P [d]",
                true,
                &[Curve { name: "P (fundamental mode)", slot: 0, get: |t| t.period_days }],
                h,
                false,
                &strip,
                Some("in instability strip"),
            ));
            take(panel(
                &mut cols[0],
                &f,
                "Base of the convective envelope",
                "r / R",
                false,
                &[Curve { name: "envelope base", slot: 0, get: |t| t.conv_env_base_r }],
                h,
                false,
                &none,
                None,
            ));
            take(panel(
                &mut cols[1],
                &f,
                "Solar-like oscillations (scaling relations)",
                "frequency [μHz]",
                true,
                &[
                    Curve { name: "ν_max", slot: 0, get: |t| Some(t.nu_max) },
                    Curve { name: "Δν", slot: 1, get: |t| Some(t.delta_nu) },
                ],
                h,
                false,
                &none,
                None,
            ));
        });
        // Current values under the plots.
        let t = &track[f.nearest(f.selected_x)];
        ui.add_space(6.0);
        ui.label(
            RichText::new(format!(
                "Shown model: {} · {} · L = {} L☉ · R = {} R☉ · T_eff = {:.0} K · X_c = {:.4}{}",
                fmt_age(t.age_yr),
                t.phase.label(),
                fmt_value(t.lum),
                fmt_value(t.radius),
                t.teff,
                t.x_c,
                match (t.pulsator, t.period_days) {
                    (Some(k), Some(p)) => format!(" · {}, P = {}", k.label(), fmt_period(p)),
                    _ => String::new(),
                }
            ))
            .color(theme::secondary(dark)),
        );
    });
    if let Some(x) = clicked {
        action = Some(HistoryAction::SelectAge(f.age_of_x(x)));
    }
    action
}
