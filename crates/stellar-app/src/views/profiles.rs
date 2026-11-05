//! Main profile view: P, T and ρ as functions of radius (or mass),
//! stacked with linked x axes and a shared hover cursor.

use egui::Id;

use super::{DENSITY, PRESSURE, TEMPERATURE, ViewCtx, quantity_plot};

pub fn show(ui: &mut egui::Ui, v: &ViewCtx<'_>) {
    let link = Id::unique("profiles_link");
    let spacing = ui.spacing().item_spacing.y;
    let h = ((ui.available_height() - 2.0 * spacing) / 3.0).max(120.0);
    egui::ScrollArea::vertical().auto_shrink(false).show(ui, |ui| {
        quantity_plot(ui, v, PRESSURE, h, Some(link), "profiles", true);
        quantity_plot(ui, v, TEMPERATURE, h, Some(link), "profiles", false);
        quantity_plot(ui, v, DENSITY, h, Some(link), "profiles", false);
    });
}
