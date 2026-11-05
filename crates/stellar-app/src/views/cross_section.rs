//! Cross-section of the star: a disc split in two halves, temperature on the
//! left (orange ramp) and density on the right (blue ramp), both on log scales.
//! Convective zones are outlined and overlaid with animated convection cells;
//! radiative zones are left plain. Hovering reads out the local state.

use std::f32::consts::{PI, TAU};

use egui::{Align2, Color32, FontId, Mesh, Pos2, Rect, Sense, Shape, Stroke, pos2, vec2};

use super::{ViewCtx, fmt_value, superscript};
use crate::theme::{self, Ramp};

pub struct CrossSectionState {
    pub animate: bool,
}

impl Default for CrossSectionState {
    fn default() -> Self {
        Self { animate: true }
    }
}

/// Interpolate `ys` at normalised radius `s` (profile radii normalised to R).
fn interp(rn: &[f64], ys: &[f64], s: f64) -> f64 {
    if s <= rn[0] {
        return ys[0];
    }
    let i = rn.partition_point(|&r| r < s).clamp(1, rn.len() - 1);
    let w = (s - rn[i - 1]) / (rn[i] - rn[i - 1]);
    ys[i - 1] + w * (ys[i] - ys[i - 1])
}

struct Field {
    log: Vec<f64>,
    lo: f64,
    hi: f64,
    ramp: Ramp,
}

impl Field {
    fn new(values: &[f64], ramp: Ramp) -> Self {
        let log: Vec<f64> = values.iter().map(|v| v.max(1e-300).log10()).collect();
        let lo = log.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = log.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Self { log, lo, hi, ramp }
    }
    fn color_at(&self, rn: &[f64], s: f64) -> Color32 {
        let v = interp(rn, &self.log, s);
        self.ramp.sample((v - self.lo) / (self.hi - self.lo).max(1e-12))
    }
}

/// Paint a half disc (`left` = x < 0) as a smooth mesh with colors per radius.
fn half_disc(mesh: &mut Mesh, c: Pos2, r_px: f32, rn: &[f64], field: &Field, left: bool) {
    const NR: usize = 180;
    const NA: usize = 90;
    let base = mesh.vertices.len() as u32;
    // Left half spans angles 90°..270° (screen y points down); right half −90°..90°.
    let a0 = if left { PI / 2.0 } else { -PI / 2.0 };
    for i in 0..=NR {
        let s = i as f64 / NR as f64;
        let color = field.color_at(rn, s);
        let rr = r_px * s as f32;
        for j in 0..=NA {
            let a = a0 + PI * j as f32 / NA as f32;
            mesh.colored_vertex(c + vec2(a.cos(), a.sin()) * rr, color);
        }
    }
    let row = (NA + 1) as u32;
    for i in 0..NR as u32 {
        for j in 0..NA as u32 {
            let v00 = base + i * row + j;
            let v01 = v00 + 1;
            let v10 = v00 + row;
            let v11 = v10 + 1;
            mesh.add_triangle(v00, v10, v11);
            mesh.add_triangle(v00, v11, v01);
        }
    }
}

/// An open circular arrow (a convection cell) centred at `c`.
fn cell(shapes: &mut Vec<Shape>, c: Pos2, radius: f32, phase: f32, clockwise: bool, stroke: Stroke) {
    let sweep = 1.55 * PI;
    let n = 20;
    let dir = if clockwise { 1.0 } else { -1.0 };
    let pts: Vec<Pos2> = (0..=n)
        .map(|k| {
            let a = phase + dir * sweep * k as f32 / n as f32;
            c + vec2(a.cos(), a.sin()) * radius
        })
        .collect();
    // Arrowhead at the end of the arc, pointing along the direction of motion.
    let end_a = phase + dir * sweep;
    let tip = c + vec2(end_a.cos(), end_a.sin()) * radius;
    let tangent = vec2(-end_a.sin(), end_a.cos()) * dir;
    let normal = vec2(end_a.cos(), end_a.sin());
    let hs = (radius * 0.45).clamp(2.5, 7.0);
    let head = vec![
        tip + tangent * hs,
        tip - tangent * hs * 0.2 + normal * hs * 0.6,
        tip - tangent * hs * 0.2 - normal * hs * 0.6,
    ];
    shapes.push(Shape::line(pts, stroke));
    shapes.push(Shape::convex_polygon(head, stroke.color, Stroke::NONE));
}

fn colorbar(painter: &egui::Painter, rect: Rect, field: &Field, title: &str, unit: &str, label_left: bool, dark: bool) {
    let (ink, muted) = (theme::primary(dark), theme::muted(dark));
    const N: usize = 64;
    let mut mesh = Mesh::default();
    for i in 0..=N {
        let t = i as f64 / N as f64;
        let y = rect.bottom() - rect.height() * t as f32;
        let color = field.ramp.sample(t);
        mesh.colored_vertex(pos2(rect.left(), y), color);
        mesh.colored_vertex(pos2(rect.right(), y), color);
    }
    for i in 0..N as u32 {
        let b = 2 * i;
        mesh.add_triangle(b, b + 1, b + 3);
        mesh.add_triangle(b, b + 3, b + 2);
    }
    painter.add(Shape::mesh(mesh));
    let font = FontId::proportional(12.0);
    let (x, align) =
        if label_left { (rect.left() - 6.0, Align2::RIGHT_CENTER) } else { (rect.right() + 6.0, Align2::LEFT_CENTER) };
    for e in (field.lo.ceil() as i64)..=(field.hi.floor() as i64) {
        let t = ((e as f64 - field.lo) / (field.hi - field.lo)) as f32;
        let y = rect.bottom() - rect.height() * t;
        let tick_x = if label_left { rect.left() - 3.0 } else { rect.right() + 3.0 };
        let edge = if label_left { rect.left() } else { rect.right() };
        painter.line_segment([pos2(edge, y), pos2(tick_x, y)], Stroke::new(1.0, muted));
        painter.text(pos2(x, y), align, format!("10{}", superscript(e)), font.clone(), muted);
    }
    painter.text(
        pos2(rect.center().x, rect.top() - 8.0),
        Align2::CENTER_BOTTOM,
        format!("{title} [{unit}]"),
        FontId::proportional(13.0),
        ink,
    );
}

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>, state: &mut CrossSectionState) {
    let dark = v.dark;
    let ink = theme::primary(dark);
    let p = &v.model.profile;
    let r_surf = *p.r.last().unwrap();
    let rn: Vec<f64> = p.r.iter().map(|r| r / r_surf).collect();

    ui.horizontal(|ui| {
        ui.checkbox(&mut state.animate, "Animate convection cells");
        ui.label(
            egui::RichText::new(
                "Left half: temperature · right half: density · both on log scales · hover to read values",
            )
            .color(theme::secondary(dark)),
        );
    });

    let (resp, painter) = ui.allocate_painter(ui.available_size(), Sense::hover());
    let area = resp.rect;
    let margin_x = 210.0_f32.min(area.width() * 0.25);
    let r_px = ((area.width() - 2.0 * margin_x).min(area.height() - 70.0) / 2.0).max(40.0);
    let c = pos2(area.center().x, area.center().y + 12.0);

    let t_field = Field::new(&p.t, Ramp::orange(dark));
    let rho_field = Field::new(&p.rho, Ramp::blue(dark));

    let mut mesh = Mesh::default();
    half_disc(&mut mesh, c, r_px, &rn, &t_field, true);
    half_disc(&mut mesh, c, r_px, &rn, &rho_field, false);
    painter.add(Shape::mesh(mesh));
    // 2px surface gap between the two halves.
    let surface = ui.visuals().panel_fill;
    painter.line_segment([c - vec2(0.0, r_px + 1.0), c + vec2(0.0, r_px + 1.0)], Stroke::new(2.0, surface));

    // Convective zones: outlines plus animated cells.
    let time = ui.input(|i| i.time) as f32;
    let cell_stroke = Stroke::new(1.6, ink.gamma_multiply(0.75));
    let mut shapes = Vec::new();
    for z in &v.model.zones {
        let (a, b) = (z.r_inner as f32 * r_px, z.r_outer as f32 * r_px);
        if z.r_inner > 1e-3 {
            painter.circle_stroke(c, a, Stroke::new(1.5, ink.gamma_multiply(0.8)));
        }
        if z.r_outer < 0.999 {
            painter.circle_stroke(c, b, Stroke::new(1.5, ink.gamma_multiply(0.8)));
        }
        let (rc, cr) = if a < 1.0 { (b * 0.5, b * 0.3) } else { (0.5 * (a + b), 0.38 * (b - a)) };
        if cr < 2.5 {
            continue;
        }
        let count = ((TAU * rc / (2.3 * cr)).floor() as usize).clamp(1, 80);
        let omega = 1.2;
        for k in 0..count {
            let ang = TAU * k as f32 / count as f32;
            let center = c + vec2(ang.cos(), ang.sin()) * if count == 1 { 0.0 } else { rc };
            let cw = k % 2 == 0;
            let phase = ang + if cw { time * omega } else { -time * omega };
            cell(&mut shapes, center, cr, phase, cw, cell_stroke);
        }
    }
    painter.extend(shapes);
    painter.circle_stroke(c, r_px, Stroke::new(1.5, ink.gamma_multiply(0.9)));

    // Half titles.
    painter.text(
        c + vec2(-r_px * 0.5, -r_px - 10.0),
        Align2::CENTER_BOTTOM,
        "Temperature",
        FontId::proportional(15.0),
        ink,
    );
    painter.text(c + vec2(r_px * 0.5, -r_px - 10.0), Align2::CENTER_BOTTOM, "Density", FontId::proportional(15.0), ink);

    // Region labels below the disc.
    let mut regions: Vec<String> = Vec::new();
    let g = &v.model.globals;
    if g.conv_core_q > 0.0 {
        regions.push(format!("Convective core: r/R < {:.3}  (m/M < {:.3})", g.conv_core_r, g.conv_core_q));
    }
    for z in &v.model.zones {
        let is_core = z.r_inner < 1e-3;
        let is_env = Some(z.r_inner) == g.conv_env_base_r;
        if !is_core && !is_env {
            regions.push(format!("Convective shell: r/R {:.3}–{:.3}", z.r_inner, z.r_outer));
        }
    }
    if let Some(rb) = g.conv_env_base_r {
        regions.push(format!("Convective envelope: r/R > {rb:.3}"));
    }
    let radiative = if v.model.zones.is_empty() {
        "Fully radiative".to_string()
    } else if v.model.zones.len() == 1 && v.model.zones[0].r_inner < 1e-3 && v.model.zones[0].r_outer > 0.999 {
        "Fully convective".to_string()
    } else {
        "Radiative zones: plain (no cells)".to_string()
    };
    regions.push(radiative);
    let text = regions.join("   ·   ");
    painter.text(
        pos2(c.x, (c.y + r_px + 14.0).min(area.bottom() - 4.0)),
        Align2::CENTER_TOP,
        text,
        FontId::proportional(13.0),
        theme::secondary(dark),
    );

    // Colorbars: T on the left, ρ on the right.
    let bar_h = r_px * 1.5;
    let bar_w = 14.0;
    let left_bar = Rect::from_min_size(pos2(c.x - r_px - 70.0, c.y - bar_h / 2.0), vec2(bar_w, bar_h));
    let right_bar = Rect::from_min_size(pos2(c.x + r_px + 56.0, c.y - bar_h / 2.0), vec2(bar_w, bar_h));
    colorbar(&painter, left_bar, &t_field, "T", "K", true, dark);
    colorbar(&painter, right_bar, &rho_field, "ρ", "g cm⁻³", false, dark);

    // Hover read-out.
    if let Some(ptr) = resp.hover_pos() {
        let d = ptr - c;
        let s = (d.length() / r_px) as f64;
        if s <= 1.0 {
            let k = rn.partition_point(|&r| r < s).min(rn.len() - 1);
            let at =
                |ys: &[f64]| 10f64.powf(interp(&rn, &ys.iter().map(|v| v.max(1e-300).log10()).collect::<Vec<_>>(), s));
            let zone = if p.convective[k] { "convective" } else { "radiative" };
            let text = format!(
                "r/R = {:.4}\nm/M = {:.4}\nT = {} K\nρ = {} g cm⁻³\nP = {} dyn cm⁻²\nL/L_s = {:.4}\n{zone}",
                s,
                interp(&rn, &p.q, s),
                fmt_value(at(&p.t)),
                fmt_value(at(&p.rho)),
                fmt_value(at(&p.p)),
                interp(&rn, &p.l, s) / p.l[p.len() - 1],
            );
            resp.on_hover_text_at_pointer(text);
        }
    }

    if state.animate && !v.model.zones.is_empty() {
        ui.ctx().request_repaint();
    }
}
