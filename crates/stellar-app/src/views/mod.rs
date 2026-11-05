//! Tab views and shared plotting helpers.

pub mod about;
pub mod compare;
pub mod cross_section;
pub mod extra;
pub mod hr;
pub mod model3d;
pub mod profiles;

use std::sync::Arc;

use egui::{Color32, Id};
use egui_plot::{Corner, GridInput, GridMark, HoverPosition, Legend, Line, Plot, PlotPoints, PlotUi, Span};
use stellar_core::{Profile, StarModel};

use crate::theme;

/// Which coordinate the profile plots use on the x axis.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum XAxis {
    /// r / R: the natural "profile" view.
    Radius,
    /// m / M: the Lagrangian coordinate the solver uses.
    Mass,
    /// log P, centre on the left. Spreads out the thin outer envelope, where
    /// opacity, gradients and convection change fastest.
    Pressure,
}

impl XAxis {
    pub const ALL: [XAxis; 3] = [XAxis::Radius, XAxis::Mass, XAxis::Pressure];

    /// Short name for selectors.
    pub fn name(self) -> &'static str {
        match self {
            XAxis::Radius => "radius r/R",
            XAxis::Mass => "mass m/M",
            XAxis::Pressure => "pressure log P",
        }
    }

    /// Axis title.
    pub fn label(self) -> &'static str {
        match self {
            XAxis::Radius => "r / R",
            XAxis::Mass => "m / M",
            XAxis::Pressure => "pressure P [dyn cm⁻²]  (centre ← → surface)",
        }
    }

    pub fn values(self, model: &StarModel) -> Vec<f64> {
        let p = &model.profile;
        match self {
            XAxis::Radius => {
                let r_surf = *p.r.last().unwrap_or(&1.0);
                p.r.iter().map(|r| r / r_surf).collect()
            }
            XAxis::Mass => p.q.clone(),
            XAxis::Pressure => p.p.iter().map(|v| v.log10()).collect(),
        }
    }

    /// Axis coordinate at normalised radius `s` (for convective-zone bounds).
    pub fn at_radius(self, model: &StarModel, s: f64) -> f64 {
        if self == XAxis::Radius {
            return s;
        }
        let p = &model.profile;
        let r_surf = *p.r.last().unwrap_or(&1.0);
        let rn: Vec<f64> = p.r.iter().map(|r| r / r_surf).collect();
        let xs = self.values(model);
        if s <= rn[0] {
            return xs[0];
        }
        let i = rn.partition_point(|&r| r < s).clamp(1, rn.len() - 1);
        let w = (s - rn[i - 1]) / (rn[i] - rn[i - 1]);
        xs[i - 1] + w * (xs[i] - xs[i - 1])
    }

    /// Hover text for an x value.
    pub fn format_x(self, x: f64) -> String {
        match self {
            XAxis::Radius => format!("r/R = {x:.4}"),
            XAxis::Mass => format!("m/M = {x:.6}"),
            XAxis::Pressure => format!("P = {} dyn cm⁻²", fmt_value(10f64.powf(x))),
        }
    }

    /// Apply axis-specific settings (orientation, ticks, full centre-to-surface
    /// extent of `model`) to a plot. A fixed extent keeps linked plots aligned,
    /// even when a series covers only part of the star.
    pub fn configure<'a>(self, plot: Plot<'a>, model: &StarModel) -> Plot<'a> {
        let plot = plot.x_axis_label(self.label());
        match self {
            XAxis::Radius | XAxis::Mass => plot.include_x(0.0).include_x(1.0),
            XAxis::Pressure => {
                let p = &model.profile.p;
                plot.invert_x(true)
                    .x_grid_spacer(log_spacer)
                    .x_axis_formatter(log_tick)
                    .include_x(p[0].log10())
                    .include_x(p[p.len() - 1].log10())
            }
        }
    }
}

/// A model pinned for comparison. Its categorical slot is assigned once at
/// pin time and never changes.
pub struct PinnedModel {
    pub model: Arc<StarModel>,
    pub slot: usize,
    pub label: String,
    pub visible: bool,
}

/// Everything a view needs to draw the current state.
pub struct ViewCtx<'a> {
    pub model: &'a Arc<StarModel>,
    pub pinned: &'a [PinnedModel],
    pub x_axis: XAxis,
    pub dark: bool,
}

impl ViewCtx<'_> {
    /// Short label for the current model.
    pub fn current_label(&self) -> String {
        model_label(self.model)
    }

    pub fn visible_pinned(&self) -> impl Iterator<Item = &PinnedModel> {
        self.pinned.iter().filter(|p| p.visible)
    }
}

pub fn model_label(m: &StarModel) -> String {
    format!("{:.2} M☉, X={:.3}, Z={:.4}", m.params.mass, m.params.x, m.params.z)
}

/// A plotted quantity: how to read it from a profile and how to label it.
#[derive(Clone, Copy)]
pub struct Quantity {
    pub title: &'static str,
    pub symbol: &'static str,
    pub unit: &'static str,
    pub log: bool,
    pub get: fn(&Profile) -> &Vec<f64>,
    /// Legend corner that the curves usually leave empty.
    pub legend: Corner,
}

pub const PRESSURE: Quantity =
    Quantity {
        title: "Pressure", symbol: "P", unit: "dyn cm⁻²", log: true, get: |p| &p.p, legend: Corner::LeftBottom
    };
pub const TEMPERATURE: Quantity =
    Quantity { title: "Temperature", symbol: "T", unit: "K", log: true, get: |p| &p.t, legend: Corner::LeftBottom };
pub const DENSITY: Quantity =
    Quantity {
        title: "Density", symbol: "ρ", unit: "g cm⁻³", log: true, get: |p| &p.rho, legend: Corner::LeftBottom
    };
pub const LUMINOSITY: Quantity = Quantity {
    title: "Luminosity",
    symbol: "L",
    unit: "erg s⁻¹",
    log: false,
    get: |p| &p.l,
    legend: Corner::RightBottom,
};
pub const OPACITY: Quantity =
    Quantity {
        title: "Opacity", symbol: "κ", unit: "cm² g⁻¹", log: true, get: |p| &p.kappa, legend: Corner::LeftTop
    };
pub const EPSILON: Quantity = Quantity {
    title: "Energy generation",
    symbol: "ε",
    unit: "erg g⁻¹ s⁻¹",
    log: true,
    get: |p| &p.eps,
    legend: Corner::RightTop,
};
pub const NABLA: Quantity = Quantity {
    title: "Temperature gradient ∇",
    symbol: "∇",
    unit: "",
    log: true,
    get: |p| &p.nabla,
    legend: Corner::LeftTop,
};
pub const V_CONV: Quantity = Quantity {
    title: "Convective velocity",
    symbol: "v_conv",
    unit: "cm s⁻¹",
    log: true,
    get: |p| &p.v_conv,
    legend: Corner::LeftTop,
};

pub const COMPARABLE: &[Quantity] = &[TEMPERATURE, PRESSURE, DENSITY, LUMINOSITY, OPACITY, EPSILON, NABLA, V_CONV];

/// (x, y) points; y is log10 for log quantities, and non-positive values are dropped.
pub fn points(xs: &[f64], ys: &[f64], log: bool) -> Vec<[f64; 2]> {
    xs.iter()
        .zip(ys)
        .filter(|(x, y)| x.is_finite() && y.is_finite() && (!log || **y > 0.0))
        .map(|(&x, &y)| [x, if log { y.log10() } else { y }])
        .collect()
}

/// Split a series into runs where `keep` holds, so gaps are not bridged.
pub fn segments(xs: &[f64], ys: &[f64], keep: &[bool], log: bool) -> Vec<Vec<[f64; 2]>> {
    let mut out = Vec::new();
    let mut cur = Vec::new();
    for i in 0..xs.len() {
        let ok = keep[i] && ys[i].is_finite() && (!log || ys[i] > 0.0);
        if ok {
            cur.push([xs[i], if log { ys[i].log10() } else { ys[i] }]);
        } else if !cur.is_empty() {
            out.push(std::mem::take(&mut cur));
        }
    }
    if !cur.is_empty() {
        out.push(cur);
    }
    out
}

/// Superscript form of an integer exponent, e.g. 17 → "¹⁷", −3 → "⁻³".
pub fn superscript(n: i64) -> String {
    const DIGITS: [char; 10] = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹'];
    let mut s = String::new();
    if n < 0 {
        s.push('⁻');
    }
    for c in n.unsigned_abs().to_string().chars() {
        s.push(DIGITS[c as usize - '0' as usize]);
    }
    s
}

/// Grid marks for an axis holding log10 values: decades (thinned to at most
/// ~8 labels), plus 2× and 5× marks when the range spans only a few decades.
pub fn log_spacer(input: GridInput) -> Vec<GridMark> {
    let (lo, hi) = input.bounds;
    let span = (hi - lo).max(1e-9);
    let step = [1.0, 2.0, 3.0, 5.0, 10.0, 20.0, 50.0].into_iter().find(|s| span / s <= 8.0).unwrap_or(100.0);
    let mut marks = Vec::new();
    let first = lo.floor() as i64;
    let last = hi.ceil() as i64;
    for e in first..=last {
        let v = e as f64;
        let major = (e as f64 / step).fract().abs() < 1e-9;
        if v >= lo && v <= hi {
            marks.push(GridMark { value: v, step_size: if major { step } else { 1.0_f64.min(step / 2.0) } });
        }
        if span <= 3.0 {
            for k in [2.0_f64, 5.0] {
                let m = v + k.log10();
                if m >= lo && m <= hi {
                    marks.push(GridMark { value: m, step_size: 0.3 });
                }
            }
        }
    }
    marks
}

/// Tick labels for an axis holding log10 values: 10ⁿ at decades and 2·10ⁿ /
/// 5·10ⁿ at the intermediate marks.
pub fn log_tick(mark: GridMark, _range: &std::ops::RangeInclusive<f64>) -> String {
    let v = mark.value;
    let e = v.floor();
    let frac = v - e;
    if frac.abs() < 1e-9 || (1.0 - frac).abs() < 1e-9 {
        format!("10{}", superscript(v.round() as i64))
    } else if (frac - 2f64.log10()).abs() < 1e-9 {
        format!("2·10{}", superscript(e as i64))
    } else if (frac - 5f64.log10()).abs() < 1e-9 {
        format!("5·10{}", superscript(e as i64))
    } else {
        String::new()
    }
}

/// Human-readable number.
pub fn fmt_value(v: f64) -> String {
    let a = v.abs();
    if a == 0.0 {
        "0".into()
    } else if (1e-2..1e4).contains(&a) {
        format!("{v:.4}")
    } else {
        format!("{v:.3e}")
    }
}

/// Hover label: series name, x and the de-logged y value with its unit.
pub fn hover_formatter(
    x_axis: XAxis,
    symbol: String,
    unit: &'static str,
    log: bool,
) -> impl Fn(&HoverPosition<'_>) -> Option<String> {
    move |pos| {
        let (name, p) = match pos {
            HoverPosition::NearDataPoint { plot_name, position, .. } => (*plot_name, *position),
            HoverPosition::Elsewhere { position } => ("", *position),
        };
        let y = if log { 10f64.powf(p.y) } else { p.y };
        let head = if name.is_empty() { String::new() } else { format!("{name}\n") };
        Some(format!("{head}{}\n{symbol} = {} {unit}", x_axis.format_x(p.x), fmt_value(y)))
    }
}

/// Shade the convective zones of `model` with neutral vertical spans. Only the
/// first plot of a group names them (one legend entry is enough).
pub fn shade_convective(plot_ui: &mut PlotUi<'_>, model: &StarModel, x_axis: XAxis, dark: bool, named: bool) {
    for z in &model.zones {
        let (a, b) = match x_axis {
            XAxis::Mass => (z.q_inner, z.q_outer),
            _ => (x_axis.at_radius(model, z.r_inner), x_axis.at_radius(model, z.r_outer)),
        };
        let (a, b) = if a <= b { (a, b) } else { (b, a) };
        plot_ui.add(
            Span::new(if named { "convective zone" } else { "" }, a..=b)
                .fill(theme::convective_fill(dark))
                .border_width(0.0),
        );
    }
}

/// A plot of one quantity: the current model in slot 0 plus visible pinned
/// models in their own slots, with the current model's convective zones shaded.
/// `legend` = false hides the legend (when a sibling plot already shows it).
pub fn quantity_plot(
    ui: &mut egui::Ui,
    v: &ViewCtx<'_>,
    q: Quantity,
    height: f32,
    link: Option<Id>,
    salt: &str,
    legend: bool,
) {
    let mut plot = Plot::new(("qplot", salt, q.symbol))
        .height(height)
        .y_axis_label(if q.unit.is_empty() { q.symbol.to_string() } else { format!("{} [{}]", q.symbol, q.unit) })
        .label_formatter(hover_formatter(v.x_axis, q.symbol.to_string(), q.unit, q.log))
        .allow_scroll(false);
    plot = v.x_axis.configure(plot, v.model);
    if legend {
        plot = plot.legend(Legend::default().position(q.legend));
    }
    if q.log {
        plot = plot.y_grid_spacer(log_spacer).y_axis_formatter(log_tick);
    }
    if let Some(id) = link {
        plot = plot.link_axis(id, [true, false]).link_cursor(id, [true, false]);
    }
    plot.show(ui, |pu| {
        shade_convective(pu, v.model, v.x_axis, v.dark, legend);
        for pin in v.visible_pinned() {
            let xs = v.x_axis.values(&pin.model);
            let pts = points(&xs, (q.get)(&pin.model.profile), q.log);
            pu.line(
                Line::new(pin.label.clone(), PlotPoints::from(pts)).color(theme::series(pin.slot, v.dark)).width(2.0),
            );
        }
        let xs = v.x_axis.values(v.model);
        let pts = points(&xs, (q.get)(&v.model.profile), q.log);
        pu.line(
            Line::new(format!("current: {}", v.current_label()), PlotPoints::from(pts))
                .color(theme::series(0, v.dark))
                .width(2.0),
        );
    });
}

/// Small colored swatch used next to text labels (identity is never color-alone).
pub fn swatch(ui: &mut egui::Ui, color: Color32) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
    ui.painter().rect_filled(rect.shrink(1.0), 2.0, color);
}

/// Caption explaining the gray convective bands, with a swatch.
pub fn convective_caption(ui: &mut egui::Ui, dark: bool) {
    ui.horizontal(|ui| {
        swatch(ui, theme::convective_fill(dark).gamma_multiply(3.0));
        ui.label(
            egui::RichText::new("gray bands: convective zones of the current model").color(theme::secondary(dark)),
        );
    });
}
