//! Compare the current model with pinned models: a table of global properties
//! (the exact-value view) and an overlay plot of any profile quantity.

use egui::RichText;
use stellar_core::Globals;

use super::{COMPARABLE, ViewCtx, fmt_value, quantity_plot, swatch};
use crate::theme;

#[derive(Default)]
pub struct CompareState {
    pub quantity: usize,
}

pub enum CompareAction {
    ToggleVisible(usize),
    Unpin(usize),
    UnpinAll,
}

type Row = (&'static str, fn(&Globals, &stellar_core::StarParams) -> String);

const ROWS: &[Row] = &[
    ("Mass [M☉]", |g, _| format!("{:.3}", g.mass)),
    ("X / Y / Z", |_, p| format!("{:.3} / {:.3} / {:.4}", p.x, p.y(), p.z)),
    ("α_MLT", |_, p| format!("{:.2}", p.alpha_mlt)),
    ("Radius [R☉]", |g, _| format!("{:.4}", g.radius)),
    ("Luminosity [L☉]", |g, _| fmt_value(g.luminosity)),
    ("T_eff [K]", |g, _| format!("{:.0}", g.teff)),
    ("log g [cgs]", |g, _| format!("{:.3}", g.log_g)),
    ("T_c [K]", |g, _| format!("{:.4e}", g.tc)),
    ("ρ_c [g cm⁻³]", |g, _| fmt_value(g.rhoc)),
    ("P_c [dyn cm⁻²]", |g, _| format!("{:.4e}", g.pc)),
    ("Convective core [m/M]", |g, _| if g.conv_core_q > 0.0 { format!("{:.4}", g.conv_core_q) } else { "none".into() }),
    ("Conv. envelope base [r/R]", |g, _| g.conv_env_base_r.map(|r| format!("{r:.4}")).unwrap_or_else(|| "none".into())),
    ("pp / CNO / 3α [%]", |g, _| {
        format!("{:.1} / {:.1} / {:.2}", 100.0 * g.frac_pp, 100.0 * g.frac_cno, 100.0 * g.frac_3a)
    }),
    ("Virial check", |g, _| format!("{:.5}", g.virial_ratio)),
];

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>, state: &mut CompareState) -> Option<CompareAction> {
    let mut action = None;
    let dark = v.dark;

    if v.pinned.is_empty() {
        ui.label(
            RichText::new(
                "No pinned models yet. Use \"Pin current model\" in the left panel to keep a model for comparison.",
            )
            .color(theme::secondary(dark)),
        );
    }

    egui::ScrollArea::horizontal().id_salt("compare_table").show(ui, |ui| {
        egui::Grid::new("compare_grid").striped(true).spacing([18.0, 4.0]).show(ui, |ui| {
            ui.label(RichText::new("Quantity").strong());
            ui.horizontal(|ui| {
                swatch(ui, theme::series(0, dark));
                ui.label(RichText::new("current").strong());
            });
            for (i, pin) in v.pinned.iter().enumerate() {
                ui.horizontal(|ui| {
                    swatch(ui, theme::series(pin.slot, dark));
                    let mut vis = pin.visible;
                    if ui
                        .checkbox(&mut vis, RichText::new(&pin.label).strong())
                        .on_hover_text("show in plots")
                        .changed()
                    {
                        action = Some(CompareAction::ToggleVisible(i));
                    }
                    if ui.small_button("×").on_hover_text("unpin").clicked() {
                        action = Some(CompareAction::Unpin(i));
                    }
                });
            }
            ui.end_row();
            for (name, f) in ROWS {
                ui.label(*name);
                ui.monospace(f(&v.model.globals, &v.model.params));
                for pin in v.pinned {
                    ui.monospace(f(&pin.model.globals, &pin.model.params));
                }
                ui.end_row();
            }
        });
    });

    ui.add_space(6.0);
    ui.horizontal(|ui| {
        ui.label("Overlay quantity:");
        egui::ComboBox::from_id_salt("compare_q").selected_text(COMPARABLE[state.quantity].title).show_ui(ui, |ui| {
            for (i, q) in COMPARABLE.iter().enumerate() {
                ui.selectable_value(&mut state.quantity, i, q.title);
            }
        });
        if !v.pinned.is_empty() && ui.button("Unpin all").clicked() {
            action = Some(CompareAction::UnpinAll);
        }
    });
    let h = ui.available_height().max(200.0);
    quantity_plot(ui, v, COMPARABLE[state.quantity], h, None, "compare", true);
    action
}
