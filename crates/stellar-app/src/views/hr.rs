//! Hertzsprung–Russell diagram: log T_eff (increasing to the left) vs log L,
//! with the ZAMS line computed for the current composition, the present-day
//! Sun, the current model and the pinned models.

use egui::Align2;
use egui_plot::{
    GridInput, GridMark, HoverPosition, Legend, Line, MarkerShape, Plot, PlotPoint, PlotPoints, Points, Polygon, Text,
};
use stellar_core::constants::TEFF_SUN;
use stellar_core::pulsation::{STRIP_MIN_LOG_L, strip_edges};
use stellar_core::{Phase, fmt_age};

use super::{ViewCtx, log_spacer, log_tick};
use crate::theme;
use crate::worker::{ZamsKey, ZamsPoint};

#[derive(Default)]
pub struct ZamsLine {
    pub key: Option<ZamsKey>,
    pub points: Vec<ZamsPoint>,
    pub done: bool,
}

const LABELLED_MASSES: &[f64] = &[0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0];

/// T_eff grid lines at round temperatures (the axis holds log10 T_eff).
fn teff_spacer(input: GridInput) -> Vec<GridMark> {
    const MAJOR: &[f64] = &[1000.0, 2000.0, 5000.0, 10_000.0, 20_000.0, 50_000.0, 100_000.0];
    const MINOR: &[f64] = &[1500.0, 3000.0, 4000.0, 6000.0, 7000.0, 8000.0, 15_000.0, 30_000.0, 40_000.0, 70_000.0];
    let (lo, hi) = input.bounds;
    let mark = |t: f64, step: f64| GridMark { value: t.log10(), step_size: step };
    MAJOR
        .iter()
        .map(|&t| mark(t, 0.3))
        .chain(MINOR.iter().map(|&t| mark(t, 0.1)))
        .filter(|m| m.value >= lo && m.value <= hi)
        .collect()
}

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>, zams: &ZamsLine) {
    let dark = v.dark;
    ui.add(
        egui::Label::new(
            egui::RichText::new(
                "ZAMS line: models for 0.3–40 M☉ with the current X, Z and α_MLT. \
                 The present-day Sun has evolved for 4.6 Gyr, so it lies above its ZAMS position. \
                 Use Evolve in the left panel to compute an evolutionary track; the shaded band is the \
                 classical instability strip, where stars pulsate.",
            )
            .color(theme::secondary(dark)),
        )
        .wrap(),
    );
    if !zams.done {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label(format!("computing ZAMS… {}/{}", zams.points.len(), crate::worker::ZAMS_MASSES.len()));
        });
    }

    // Explicit bounds. Normally they contain every point (the whole diagram),
    // so the view fits the data even while the ZAMS line is still arriving.
    // With the adaptive scale during a run, they zoom in on the track.
    let (mut xs, mut ys): (Vec<f64>, Vec<f64>) = if v.zoom_track && !v.track.is_empty() {
        (Vec::new(), Vec::new())
    } else {
        let mut xs: Vec<f64> = zams.points.iter().map(|p| p.teff.log10()).collect();
        let mut ys: Vec<f64> = zams.points.iter().map(|p| p.luminosity.log10()).collect();
        xs.push(TEFF_SUN.log10());
        ys.push(0.0);
        (xs, ys)
    };
    for m in std::iter::once(&**v.model).chain(v.visible_pinned().map(|p| &*p.model)) {
        xs.push(m.globals.teff.log10());
        ys.push(m.globals.luminosity.log10());
    }
    for t in v.track {
        xs.push(t.teff.log10());
        ys.push(t.lum.log10());
    }
    let (x0, x1) =
        (xs.iter().cloned().fold(f64::INFINITY, f64::min), xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    let (y0, y1) =
        (ys.iter().cloned().fold(f64::INFINITY, f64::min), ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max));
    // Keep a sensible minimum extent when zoomed on a short track.
    let (x0, x1) = ((x0.min(0.5 * (x0 + x1) - 0.03)), x1.max(0.5 * (x0 + x1) + 0.03));
    let (y0, y1) = ((y0.min(0.5 * (y0 + y1) - 0.15)), y1.max(0.5 * (y0 + y1) + 0.15));
    let zoom = v.zoom_track && !v.track.is_empty();
    let pad_x = if zoom { 0.1 * (x1 - x0) } else { 0.04 };
    let pad_y = if zoom { 0.1 * (y1 - y0) } else { 0.3 };
    let target = [x0 - pad_x, x1 + pad_x, y0 - pad_y, y1 + pad_y];

    let mut plot = Plot::new("hr_diagram");
    if zoom {
        // Zoomed on the track: explicit default bounds (other items, such as
        // the full ZAMS line, do not widen the view). Reset whenever the track's
        // extent changes, so the view follows a run that is still computing.
        let key = egui::Id::unique("hr_zoom_target");
        let changed = ui.data_mut(|d| {
            let prev: Option<[f64; 4]> = d.get_temp(key);
            d.insert_temp(key, target);
            prev != Some(target)
        });
        plot = plot.default_x_bounds(target[0], target[1]).default_y_bounds(target[2], target[3]);
        if changed {
            plot = plot.reset();
        }
    } else {
        plot = plot.include_x(target[0]).include_x(target[1]).include_y(target[2]).include_y(target[3]);
    }
    if v.reset_plots {
        plot = plot.reset();
    }
    plot.invert_x(true)
        .legend(Legend::default().position(if zoom {
            egui_plot::Corner::LeftBottom
        } else {
            egui_plot::Corner::LeftTop
        }))
        .x_axis_label("effective temperature T_eff [K]")
        .y_axis_label("luminosity L [L☉]")
        .x_grid_spacer(teff_spacer)
        .x_axis_formatter(|mark, _| {
            let t = (10f64.powf(mark.value) / 100.0).round() * 100.0;
            if t >= 1e4 { format!("{:.0}k", t / 1e3) } else { format!("{t:.0}") }
        })
        .y_grid_spacer(log_spacer)
        .y_axis_formatter(log_tick)
        .label_formatter(|pos: &HoverPosition<'_>| {
            if let HoverPosition::NearDataPoint { plot_name: "evolutionary track", index, .. } = pos
                && let Some(t) = v.track.get(*index)
            {
                let puls = match (t.pulsator, t.period_days) {
                    (Some(k), Some(p)) => format!("\n{}, P = {}", k.label(), crate::evolution::fmt_period(p)),
                    _ => String::new(),
                };
                return Some(format!(
                    "{}\n{}\nT_eff = {:.0} K, L = {} L☉\nR = {} R☉, X_c = {:.4}{puls}",
                    fmt_age(t.age_yr),
                    t.phase.label(),
                    t.teff,
                    super::fmt_value(t.lum),
                    super::fmt_value(t.radius),
                    t.x_c
                ));
            }
            let (name, p) = match pos {
                HoverPosition::NearDataPoint { plot_name, position, .. } => (*plot_name, *position),
                HoverPosition::Elsewhere { position } => ("", *position),
            };
            let head = if name.is_empty() { String::new() } else { format!("{name}\n") };
            Some(format!("{head}T_eff = {:.0} K\nL = {} L☉", 10f64.powf(p.x), super::fmt_value(10f64.powf(p.y))))
        })
        .show(ui, |pu| {
            let muted = theme::muted(dark);
            let ink = theme::primary(dark);

            // Classical instability strip (approximate), as a neutral band.
            let (l_lo, l_hi) = (STRIP_MIN_LOG_L, 5.5);
            let edge = |blue: bool| {
                (0..=20).map(move |i| {
                    let ll = l_lo + (l_hi - l_lo) * i as f64 / 20.0;
                    let (b, r) = strip_edges(ll);
                    [if blue { b } else { r }, ll]
                })
            };
            let band: Vec<[f64; 2]> = edge(true).chain(edge(false).collect::<Vec<_>>().into_iter().rev()).collect();
            pu.add(
                Polygon::new("instability strip (approx.)", PlotPoints::from(band))
                    .fill_color(theme::convective_fill(dark))
                    .stroke(egui::Stroke::new(1.0, muted.gamma_multiply(0.6))),
            );

            let mut pts = zams.points.clone();
            pts.sort_by(|a, b| a.mass.total_cmp(&b.mass));
            if pts.len() >= 2 {
                let line: Vec<[f64; 2]> = pts.iter().map(|p| [p.teff.log10(), p.luminosity.log10()]).collect();
                pu.line(Line::new("ZAMS (computed)", PlotPoints::from(line)).color(muted).width(2.0));
            }
            for p in &pts {
                if LABELLED_MASSES.iter().any(|m| (m - p.mass).abs() < 1e-9) {
                    pu.add(
                        Text::new("", PlotPoint::new(p.teff.log10(), p.luminosity.log10()), format!("  {} M☉", p.mass))
                            .color(muted)
                            .anchor(Align2::LEFT_TOP),
                    );
                }
            }

            // Present-day Sun: open circle with a centre dot (the ☉ symbol).
            let sun = [[TEFF_SUN.log10(), 0.0]];
            pu.add(
                Points::new("Sun (today)", PlotPoints::from(sun.to_vec()))
                    .shape(MarkerShape::Circle)
                    .filled(false)
                    .radius(8.0)
                    .color(ink),
            );
            pu.add(
                Points::new("Sun (today)", PlotPoints::from(sun.to_vec()))
                    .shape(MarkerShape::Circle)
                    .radius(2.0)
                    .color(ink)
                    .id("sun_dot"),
            );

            for pin in v.visible_pinned() {
                let g = &pin.model.globals;
                pu.add(
                    Points::new(pin.label.clone(), PlotPoints::from(vec![[g.teff.log10(), g.luminosity.log10()]]))
                        .shape(MarkerShape::Circle)
                        .radius(6.0)
                        .color(theme::series(pin.slot, dark)),
                );
            }
            // Evolutionary track of the current run, with labelled ages.
            if v.track.len() >= 2 {
                let line: Vec<[f64; 2]> = v.track.iter().map(|t| [t.teff.log10(), t.lum.log10()]).collect();
                pu.line(
                    Line::new("evolutionary track", PlotPoints::from(line)).color(theme::series(0, dark)).width(2.0),
                );
                let mut marks = vec![0usize];
                if let Some(i) = v.track.iter().position(|t| t.phase != Phase::MainSequence && t.phase != Phase::Zams) {
                    marks.push(i);
                }
                marks.push(v.track.len() - 1);
                marks.dedup();
                for i in marks {
                    let t = &v.track[i];
                    let label = match (i, t.phase) {
                        (0, _) => "  ZAMS".to_string(),
                        (_, Phase::ShellHydrogen | Phase::RedGiant | Phase::CoreHelium) if i + 1 < v.track.len() => {
                            format!("  end of main sequence, {}", fmt_age(t.age_yr))
                        }
                        _ => format!("  {}", fmt_age(t.age_yr)),
                    };
                    pu.add(
                        Text::new("", PlotPoint::new(t.teff.log10(), t.lum.log10()), label)
                            .color(ink)
                            .anchor(Align2::LEFT_BOTTOM),
                    );
                }
            }

            let g = &v.model.globals;
            pu.add(
                Points::new(
                    format!("current: {}", v.current_label()),
                    PlotPoints::from(vec![[g.teff.log10(), g.luminosity.log10()]]),
                )
                .shape(MarkerShape::Diamond)
                .radius(8.0)
                .color(theme::series(0, dark)),
            );
        });
}
