//! Star Models — interactive zero-age main-sequence stellar structure explorer.

mod app;
mod theme;
mod views;
mod worker;

fn main() -> eframe::Result {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("Star Models")
            .with_inner_size([1440.0, 920.0])
            .with_min_inner_size([960.0, 640.0]),
        // Anti-aliased edges for the software-rendered 3D mesh.
        multisampling: 4,
        ..Default::default()
    };
    eframe::run_native("Star Models", options, Box::new(|cc| Ok(Box::new(app::StarApp::new(cc)))))
}
