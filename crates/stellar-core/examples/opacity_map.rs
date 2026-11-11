//! Opacity diagnostic: the analytic Rosseland-mean opacity across regimes.
//!
//! Writes an SVG with two panels:
//! * κ against log T at fixed log R = log ρ − 3 log T₆ (the slicing used by
//!   OPAL/OP tables), with the solar ZAMS model's own κ(T) path overlaid;
//! * the individual contributions (electron scattering, free–free,
//!   bound–free, H⁻, molecular floor) at log R = −1.5, near the solar path.
//!
//! It also prints κ at a few reference points.
//!
//! ```text
//! cargo run --release -p stellar-core --example opacity_map -- docs/about/images/opacity.svg
//! ```

use std::fmt::Write as _;

use stellar_core::physics::Opacity;
use stellar_core::physics::opacity::AnalyticOpacity;
use stellar_core::{Composition, NoControl, Physics, StarParams, solve_model};

const W: f64 = 980.0;
const H: f64 = 470.0;
const PANEL_W: f64 = 400.0;
const PANEL_H: f64 = 330.0;
const TOP: f64 = 62.0;
const LEFT: [f64; 2] = [70.0, 560.0];
const LOG_T: (f64, f64) = (3.4, 8.0);
const LOG_K: (f64, f64) = (-3.0, 6.0);
const LOG_RS: [f64; 5] = [-5.0, -4.0, -3.0, -2.0, -1.0];
const LOG_R_COMPONENTS: f64 = -1.5;

fn rho_of(log_r: f64, log_t: f64) -> f64 {
    10f64.powf(log_r + 3.0 * (log_t - 6.0))
}

struct Panel {
    left: f64,
}

impl Panel {
    fn x(&self, log_t: f64) -> f64 {
        self.left + (log_t - LOG_T.0) / (LOG_T.1 - LOG_T.0) * PANEL_W
    }
    fn y(&self, log_k: f64) -> f64 {
        TOP + (LOG_K.1 - log_k) / (LOG_K.1 - LOG_K.0) * PANEL_H
    }
    fn path(&self, pts: &[(f64, f64)]) -> String {
        let mut d = String::new();
        for (i, (t, k)) in pts.iter().enumerate() {
            let _ = write!(d, "{}{:.1},{:.1} ", if i == 0 { "M" } else { "L" }, self.x(*t), self.y(*k));
        }
        d
    }
}

/// Curve sampled in log T.
fn curve(f: impl Fn(f64) -> f64) -> Vec<(f64, f64)> {
    let n = 460;
    (0..=n)
        .map(|i| {
            let lt = LOG_T.0 + (LOG_T.1 - LOG_T.0) * i as f64 / n as f64;
            (lt, f(lt).log10())
        })
        .collect()
}

fn axes(svg: &mut String, p: &Panel, clip: &str, title: &str, subtitle: &str) {
    let (x0, x1, y0, y1) = (p.left, p.left + PANEL_W, TOP, TOP + PANEL_H);
    let _ = writeln!(
        svg,
        r#"<clipPath id="{clip}"><rect x="{x0}" y="{y0}" width="{PANEL_W}" height="{PANEL_H}"/></clipPath>"#
    );
    let _ = writeln!(svg, r#"<text class="title" x="{x0}" y="22">{title}</text>"#);
    let _ = writeln!(svg, r#"<text class="sub" x="{x0}" y="40">{subtitle}</text>"#);
    for k in (LOG_K.0 as i32)..=(LOG_K.1 as i32) {
        let y = p.y(k as f64);
        let _ = writeln!(svg, r#"<line class="grid" x1="{x0}" x2="{x1}" y1="{y:.1}" y2="{y:.1}"/>"#);
        let _ =
            writeln!(svg, r#"<text class="tick" x="{:.1}" y="{:.1}" text-anchor="end">{k}</text>"#, x0 - 8.0, y + 4.0);
    }
    for t in (LOG_T.0.ceil() as i32)..=(LOG_T.1.floor() as i32) {
        let x = p.x(t as f64);
        let _ = writeln!(svg, r#"<text class="tick" x="{x:.1}" y="{:.1}" text-anchor="middle">{t}</text>"#, y1 + 18.0);
    }
    let _ = writeln!(svg, r#"<line class="axis" x1="{x0}" x2="{x1}" y1="{y1}" y2="{y1}"/>"#);
    let _ = writeln!(
        svg,
        r#"<text class="label" x="{:.1}" y="{:.1}" text-anchor="middle">log T [K]</text>"#,
        x0 + PANEL_W / 2.0,
        y1 + 40.0
    );
    let _ = writeln!(
        svg,
        r#"<text class="label" transform="translate({:.1},{:.1}) rotate(-90)" text-anchor="middle">log κ [cm² g⁻¹]</text>"#,
        x0 - 42.0,
        y0 + PANEL_H / 2.0
    );
}

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "opacity.svg".into());
    let comp = Composition::new(0.70, 0.02, 0.7).unwrap();
    let kappa = |rho: f64, t: f64| AnalyticOpacity.kappa(rho, t, &comp);

    let mut svg = String::new();
    let _ = writeln!(
        svg,
        r##"<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="system-ui, -apple-system, 'Segoe UI', sans-serif">
<style>
  .bg {{ fill: #fcfcfb; }}
  .title {{ font-size: 15px; font-weight: 600; fill: #0b0b0b; }}
  .sub, .label {{ font-size: 12px; fill: #52514e; }}
  .tick {{ font-size: 11px; fill: #898781; }}
  .note {{ font-size: 11px; fill: #52514e; }}
  .dl {{ font-size: 11px; fill: #52514e; }}
  .grid {{ stroke: #e1e0d9; stroke-width: 1; }}
  .axis {{ stroke: #898781; stroke-width: 1; }}
  .line {{ fill: none; stroke-width: 2; stroke-linejoin: round; stroke-linecap: round; }}
  .r0 {{ stroke: #86b6ef; }} .r1 {{ stroke: #5598e7; }} .r2 {{ stroke: #2a78d6; }} .r3 {{ stroke: #1c5cab; }} .r4 {{ stroke: #104281; }}
  .sun {{ stroke: #eb6834; stroke-width: 3; }}
  .c1 {{ stroke: #2a78d6; }} .c2 {{ stroke: #eb6834; }} .c3 {{ stroke: #1baf7a; }} .c4 {{ stroke: #eda100; }} .c5 {{ stroke: #e87ba4; }}
  .total {{ stroke: #0b0b0b; stroke-width: 2.5; }}
  .cross {{ stroke: #898781; stroke-width: 1; stroke-dasharray: 3 3; }}
  @media (prefers-color-scheme: dark) {{
    .bg {{ fill: #1a1a19; }}
    .title {{ fill: #f0efec; }}
    .sub, .label, .note, .dl {{ fill: #c3c2b7; }}
    .grid {{ stroke: #2c2c2a; }}
    .r0 {{ stroke: #184f95; }} .r1 {{ stroke: #256abf; }} .r2 {{ stroke: #3987e5; }} .r3 {{ stroke: #6da7ec; }} .r4 {{ stroke: #b7d3f6; }}
    .sun {{ stroke: #d95926; }}
    .c1 {{ stroke: #3987e5; }} .c2 {{ stroke: #d95926; }} .c3 {{ stroke: #199e70; }} .c4 {{ stroke: #c98500; }} .c5 {{ stroke: #d55181; }}
    .total {{ stroke: #f0efec; }}
  }}
</style>
<title>Analytic opacity across regimes</title>
<rect class="bg" width="{W}" height="{H}"/>"##
    );

    // Panel A: total κ at fixed log R, plus the solar ZAMS path.
    let a = Panel { left: LEFT[0] };
    axes(
        &mut svg,
        &a,
        "clipA",
        "Total κ at fixed log R = log ρ − 3 log T₆",
        "X = 0.70, Z = 0.02 · orange: the 1 M☉ ZAMS model along its own ρ–T path",
    );
    let _ = writeln!(svg, r#"<g clip-path="url(#clipA)">"#);
    let curves: Vec<Vec<(f64, f64)>> =
        LOG_RS.iter().map(|&lr| curve(|lt| kappa(rho_of(lr, lt), 10f64.powf(lt)))).collect();
    for (i, pts) in curves.iter().enumerate() {
        let _ = writeln!(svg, r#"<path class="line r{i}" d="{}"/>"#, a.path(pts));
    }
    let physics = Physics::analytic();
    let sun = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).expect("solar model");
    let pr = &sun.profile;
    let sun_pts: Vec<(f64, f64)> = (0..pr.len()).map(|k| (pr.t[k].log10(), pr.kappa[k].log10())).collect();
    let _ = writeln!(svg, r#"<path class="line sun" d="{}"/>"#, a.path(&sun_pts));
    let _ = writeln!(svg, "</g>");
    // Legend in the empty upper-right corner (the curves converge at high T,
    // so direct end labels would collide there).
    let mut entries: Vec<(String, String)> =
        LOG_RS.iter().enumerate().rev().map(|(i, lr)| (format!("r{i}"), format!("log R = {lr}"))).collect();
    entries.push(("sun".into(), format!("1 M☉ ZAMS model (centre κ = {:.1})", pr.kappa[0])));
    let (lx, mut ly) = (a.x(5.3), a.y(5.75));
    for (class, text) in &entries {
        let _ = writeln!(
            svg,
            r#"<line class="line {class}" x1="{lx:.1}" x2="{:.1}" y1="{ly:.1}" y2="{ly:.1}"/><text class="dl" x="{:.1}" y="{:.1}">{text}</text>"#,
            lx + 22.0,
            lx + 28.0,
            ly + 4.0
        );
        ly += 15.0;
    }

    // Panel B: components at log R = −1.5.
    let b = Panel { left: LEFT[1] };
    axes(
        &mut svg,
        &b,
        "clipB",
        &format!("Contributions at log R = {LOG_R_COMPONENTS}"),
        "dashed: H⁻ equals Kramers + scattering; the blend is half of each",
    );
    // Curve classes: es, ff, bf, H⁻ (indices into `components`), molecular floor, total.
    let classes = ["c1", "c2", "c3", "c4", "c5", "total"];
    let value = |i: usize, lt: f64| -> f64 {
        let (rho, t) = (rho_of(LOG_R_COMPONENTS, lt), 10f64.powf(lt));
        match i {
            0..=3 => AnalyticOpacity::components(rho, t, &comp)[i],
            4 => 0.1 * comp.z,
            _ => kappa(rho, t),
        }
    };
    // Where H⁻ equals the Kramers + scattering sum, the blend is half of each.
    let cross = (0..4000).map(|i| LOG_T.0 + (LOG_T.1 - LOG_T.0) * i as f64 / 4000.0).find(|&lt| {
        let c = AnalyticOpacity::components(rho_of(LOG_R_COMPONENTS, lt), 10f64.powf(lt), &comp);
        c[3] >= c[0] + c[1] + c[2]
    });
    if let Some(lt) = cross {
        let x = b.x(lt);
        let _ = writeln!(svg, r#"<line class="cross" x1="{x:.1}" x2="{x:.1}" y1="{TOP}" y2="{}"/>"#, TOP + PANEL_H);
    }
    let _ = writeln!(svg, r#"<g clip-path="url(#clipB)">"#);
    let comp_pts: Vec<Vec<(f64, f64)>> = (0..classes.len()).map(|i| curve(|lt| value(i, lt))).collect();
    for (class, pts) in classes.iter().zip(&comp_pts) {
        let _ = writeln!(svg, r#"<path class="line {class}" d="{}"/>"#, b.path(pts));
    }
    let _ = writeln!(svg, "</g>");
    // Direct labels at positions (log T, log κ) chosen to avoid the curves.
    for (name, lt, lk) in [
        ("electron scattering", 4.9, -0.25),
        ("free–free", 5.6, 0.45),
        ("bound–free", 5.6, 2.65),
        ("H⁻", 4.38, 5.3),
        ("molecular floor", 6.2, -2.5),
        ("total", 4.45, 3.65),
    ] {
        let _ = writeln!(svg, r#"<text class="dl" x="{:.1}" y="{:.1}">{name}</text>"#, b.x(lt), b.y(lk));
    }
    let _ = writeln!(
        svg,
        r#"<text class="dl" x="{:.1}" y="{:.1}">κ = κ_mol + [1/κ_H⁻ + 1/(κ_es + κ_ff + κ_bf)]⁻¹</text>"#,
        b.x(5.05),
        b.y(5.4)
    );

    let _ = writeln!(
        svg,
        r#"<text class="note" x="{}" y="{}">Analytic opacity of Star Models (crates/stellar-core/src/physics/opacity.rs). Generated by examples/opacity_map.rs.</text>"#,
        LEFT[0],
        H - 12.0
    );
    svg.push_str("</svg>\n");
    std::fs::write(&out, svg).expect("writing SVG");
    eprintln!("written {out}");

    println!("| point | T [K] | ρ [g cm⁻³] | log R | κ analytic [cm² g⁻¹] |");
    println!("|---|---:|---:|---:|---:|");
    let row = |name: &str, t: f64, rho: f64| {
        let lr = rho.log10() - 3.0 * (t / 1e6).log10();
        println!("| {name} | {t:.3e} | {rho:.3e} | {lr:.2} | {:.3} |", kappa(rho, t));
    };
    row("1 M☉ ZAMS model centre", pr.t[0], pr.rho[0]);
    row("present Sun, centre", 1.57e7, 150.0);
    row("present Sun, base of convection zone", 2.2e6, 0.19);
    row("present Sun, r ≈ 0.5 R", 3.9e6, 1.3);
    if let Some(lt) = cross {
        println!("\nH⁻ / Kramers crossover at log R = {LOG_R_COMPONENTS}: log T = {lt:.3}");
    }
}
