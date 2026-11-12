//! Additional profiles of the current model: luminosity and mass, opacity,
//! energy generation by process, temperature gradients (why each layer is
//! convective), composition, gravitational energy, convective velocity,
//! convective flux and β, the propagation diagram (N, S₁, S₂) and the sound speed.

use egui::{Id, RichText};
use egui_plot::{Corner, Legend, Line, Plot, PlotPoints};

use crate::evolution::{Range, RunRanges};

use super::{ViewCtx, convective_caption, hover_formatter, log_spacer, log_tick, points, segments, shade_convective};
use crate::theme;

struct Series {
    name: String,
    slot: usize,
    pts: Vec<Vec<[f64; 2]>>,
}

#[allow(clippy::too_many_arguments)]
fn panel(
    ui: &mut egui::Ui,
    v: &ViewCtx<'_>,
    title: &str,
    y_label: &str,
    log: bool,
    series: Vec<Series>,
    height: f32,
    link: Id,
    legend: Corner,
    range: fn(&RunRanges) -> Range,
) {
    ui.label(RichText::new(title).strong());
    let mut plot = Plot::new(("extra", title))
        .height(height)
        .legend(Legend::default().position(legend))
        .y_axis_label(y_label)
        .allow_scroll(false)
        .link_axis(link, [true, false])
        .link_cursor(link, [true, false])
        .label_formatter(hover_formatter(v.x_axis, "value".to_string(), "", log));
    plot = v.x_axis.configure(plot, v.model, v.fixed);
    plot = v.apply_scale(plot, range);
    if log {
        plot = plot.y_grid_spacer(log_spacer).y_axis_formatter(log_tick);
    } else {
        // Readable ticks for very small or large linear values (e.g. ε_grav).
        plot = plot.y_axis_formatter(|mark, _| {
            let v = mark.value;
            let a = v.abs();
            if a == 0.0 || (1e-2..1e4).contains(&a) { format!("{v:.2}") } else { format!("{v:.1e}") }
        });
    }
    plot.show(ui, |pu| {
        shade_convective(pu, v.model, v.x_axis, v.dark, false);
        for s in series {
            let color = theme::series(s.slot, v.dark);
            for (i, seg) in s.pts.into_iter().enumerate() {
                // Segments share the series name, so the legend shows one entry.
                pu.line(
                    Line::new(s.name.clone(), PlotPoints::from(seg)).color(color).width(2.0).id(("seg", s.slot, i)),
                );
            }
        }
    });
}

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>) {
    let p = &v.model.profile;
    let xs = v.x_axis.values(v.model);
    let n = p.len();
    let link = Id::unique("extra_link");
    let all = vec![true; n];

    // Characteristic frequencies as cyclic frequencies [μHz].
    let to_uhz = 1e6 / std::f64::consts::TAU;
    // (Models read from old files may lack c_s and N²: treat them as empty.)
    let n2 = |k: usize| p.n2.get(k).copied().unwrap_or(0.0);
    let nu_n: Vec<f64> = (0..n).map(|k| n2(k).max(0.0).sqrt() * to_uhz).collect();
    let stable: Vec<bool> = (0..n).map(|k| n2(k) > 0.0).collect();
    let lamb = |l: u32| -> Vec<f64> {
        p.cs.iter().zip(&p.r).map(|(&c, &r)| stellar_core::seismic::lamb_frequency(l, c, r) * to_uhz).collect()
    };
    let (nu_s1, nu_s2) = (lamb(1), lamb(2));

    let l_surf = p.l[n - 1];
    let l_norm: Vec<f64> = p.l.iter().map(|l| l / l_surf).collect();

    // Hide energy generation more than 8 dex below the peak.
    let eps_floor = p.eps.iter().cloned().fold(0.0, f64::max) * 1e-8;
    let above = |ys: &[f64]| -> Vec<bool> { ys.iter().map(|&e| e > eps_floor).collect() };

    convective_caption(ui, v.dark);
    let label_rows = 3.0 * (ui.text_style_height(&egui::TextStyle::Body) + ui.spacing().item_spacing.y);
    // Three rows fit the window; the rows below (convection, seismic) are reached by scrolling.
    let h = ((ui.available_height() - label_rows - 3.0 * ui.spacing().item_spacing.y) / 3.0).max(110.0);

    egui::ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        ui.columns(2, |cols| {
            panel(
                &mut cols[0],
                v,
                "Luminosity and enclosed mass",
                "fraction",
                false,
                vec![
                    Series { name: "L / L_surface".into(), slot: 0, pts: vec![points(&xs, &l_norm, false)] },
                    Series { name: "m / M".into(), slot: 1, pts: vec![points(&xs, &p.q, false)] },
                ],
                h,
                link,
                Corner::RightBottom,
                |_| Range::UNIT,
            );
            panel(
                &mut cols[1],
                v,
                "Rosseland opacity",
                "κ [cm² g⁻¹]",
                true,
                vec![Series { name: "κ".into(), slot: 0, pts: vec![points(&xs, &p.kappa, true)] }],
                h,
                link,
                Corner::LeftTop,
                |r| r.log_kappa,
            );
            panel(
                &mut cols[0],
                v,
                "Nuclear energy generation",
                "ε [erg g⁻¹ s⁻¹]",
                true,
                vec![
                    Series { name: "pp chain".into(), slot: 0, pts: segments(&xs, &p.eps_pp, &above(&p.eps_pp), true) },
                    Series {
                        name: "CNO cycle".into(),
                        slot: 1,
                        pts: segments(&xs, &p.eps_cno, &above(&p.eps_cno), true),
                    },
                    Series {
                        name: "triple-α".into(), slot: 2, pts: segments(&xs, &p.eps_3a, &above(&p.eps_3a), true)
                    },
                ],
                h,
                link,
                Corner::RightTop,
                |r| r.log_eps,
            );
            panel(
                &mut cols[1],
                v,
                "Temperature gradients (convective where ∇_rad > ∇_ad)",
                "∇ = d ln T / d ln P",
                true,
                vec![
                    Series { name: "∇_rad (radiative)".into(), slot: 0, pts: vec![points(&xs, &p.nabla_rad, true)] },
                    Series { name: "∇_ad (adiabatic)".into(), slot: 1, pts: vec![points(&xs, &p.nabla_ad, true)] },
                    Series { name: "∇ (actual, MLT)".into(), slot: 2, pts: vec![points(&xs, &p.nabla, true)] },
                ],
                h,
                link,
                Corner::LeftTop,
                |r| r.log_nabla,
            );
            panel(
                &mut cols[0],
                v,
                "Chemical composition (changes as the star evolves)",
                "mass fraction",
                false,
                vec![
                    Series { name: "X (hydrogen)".into(), slot: 0, pts: vec![points(&xs, &p.x, false)] },
                    Series { name: "Y (helium)".into(), slot: 1, pts: vec![points(&xs, &p.y, false)] },
                    Series { name: "Z (metals)".into(), slot: 2, pts: vec![points(&xs, &p.z, false)] },
                ],
                h,
                link,
                Corner::RightBottom,
                |_| Range::UNIT,
            );
            panel(
                &mut cols[1],
                v,
                "Gravitational energy release ε_grav = −T ds/dt (> 0: contraction)",
                "ε [erg g⁻¹ s⁻¹]",
                false,
                vec![Series { name: "ε_grav".into(), slot: 0, pts: vec![points(&xs, &p.eps_grav, false)] }],
                h,
                link,
                Corner::RightTop,
                |r| r.eps_grav,
            );
            panel(
                &mut cols[0],
                v,
                "Convective velocity (mixing-length theory)",
                "v_conv [cm s⁻¹]",
                true,
                vec![Series { name: "v_conv".into(), slot: 0, pts: segments(&xs, &p.v_conv, &p.convective, true) }],
                h,
                link,
                Corner::LeftTop,
                |r| r.log_vconv,
            );
            panel(
                &mut cols[1],
                v,
                "Convective flux fraction and gas-pressure fraction",
                "fraction",
                false,
                vec![
                    Series { name: "F_conv / F".into(), slot: 0, pts: segments(&xs, &p.f_conv, &all, false) },
                    Series { name: "β = P_gas / P".into(), slot: 1, pts: vec![points(&xs, &p.beta, false)] },
                ],
                h,
                link,
                Corner::LeftBottom,
                |_| Range::UNIT,
            );
            panel(
                &mut cols[0],
                v,
                "Propagation diagram: p modes above N and S_l, g modes below both",
                "cyclic frequency [μHz]",
                true,
                vec![
                    Series {
                        name: "N/2π (buoyancy, where N² > 0)".into(),
                        slot: 0,
                        pts: segments(&xs, &nu_n, &stable, true),
                    },
                    Series { name: "S₁/2π (Lamb, l = 1)".into(), slot: 1, pts: vec![points(&xs, &nu_s1, true)] },
                    Series { name: "S₂/2π (Lamb, l = 2)".into(), slot: 2, pts: vec![points(&xs, &nu_s2, true)] },
                ],
                h,
                link,
                Corner::RightTop,
                |_| Range::EMPTY,
            );
            panel(
                &mut cols[1],
                v,
                "Adiabatic sound speed c_s = √(Γ₁ P / ρ)",
                "c_s [cm s⁻¹]",
                true,
                vec![Series { name: "c_s".into(), slot: 0, pts: vec![points(&xs, &p.cs, true)] }],
                h,
                link,
                Corner::LeftBottom,
                |_| Range::EMPTY,
            );
        });
        cavities_table(ui, v);
    });
}

/// The candidate g-mode cavities of the current model, with their asymptotic ΔΠ₁.
fn cavities_table(ui: &mut egui::Ui, v: &ViewCtx<'_>) {
    let g = &v.model.globals;
    ui.add_space(4.0);
    ui.label(RichText::new("Asymptotic spacings").strong());
    ui.label(format!(
        "Δν = {:.2} μHz from [2 ∫ dr / c_s]⁻¹ (scaling relation: {:.2} μHz). Asymptotic, dominated by the outer layers.",
        g.delta_nu_asym_uhz, g.delta_nu_uhz
    ));
    if g.g_cavities.is_empty() {
        ui.label("No region with N² > 0: no candidate g-mode cavity.");
        return;
    }
    egui::Grid::new("g_cavities").num_columns(4).striped(true).spacing([16.0, 3.0]).show(ui, |ui| {
        ui.label(RichText::new("candidate cavity").strong());
        ui.label(RichText::new("r / R").strong());
        ui.label(RichText::new("m / M").strong());
        ui.label(RichText::new("ΔΠ₁ [s]").strong());
        ui.end_row();
        for c in &g.g_cavities {
            ui.label(c.label());
            ui.monospace(format!("{:.4} – {:.4}", c.r_inner, c.r_outer));
            ui.monospace(format!("{:.4} – {:.4}", c.q_inner, c.q_outer));
            ui.monospace(format!("{:.1}", c.delta_pi1_s));
            ui.end_row();
        }
    });
    ui.label(
        RichText::new(
            "Each connected N² > 0 region is a separate candidate cavity; they are not summed. ΔΠ₁ = 2π²/√2 · \
             [∫ N dr/r]⁻¹ over the whole region is the low-frequency asymptotic estimate, not a predicted spectrum.",
        )
        .color(theme::secondary(v.dark)),
    );
}
