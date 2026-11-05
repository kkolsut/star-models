//! Additional profiles of the current model: luminosity and mass, opacity,
//! energy generation by process, temperature gradients (why each layer is
//! convective), convective velocity, convective flux and β.

use egui::{Id, RichText};
use egui_plot::{Corner, Legend, Line, Plot, PlotPoints};

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
    plot = v.x_axis.configure(plot, v.model);
    if log {
        plot = plot.y_grid_spacer(log_spacer).y_axis_formatter(log_tick);
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

    let l_surf = p.l[n - 1];
    let l_norm: Vec<f64> = p.l.iter().map(|l| l / l_surf).collect();

    // Hide energy generation more than 8 dex below the peak.
    let eps_floor = p.eps.iter().cloned().fold(0.0, f64::max) * 1e-8;
    let above = |ys: &[f64]| -> Vec<bool> { ys.iter().map(|&e| e > eps_floor).collect() };

    convective_caption(ui, v.dark);
    let label_rows = 3.0 * (ui.text_style_height(&egui::TextStyle::Body) + ui.spacing().item_spacing.y);
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
            );
        });
    });
}
