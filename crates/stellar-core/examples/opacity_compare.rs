//! Plan 7.3: the tabulated opacity against the analytic formulas at identical
//! (ρ, T, X, Z), along the analytic 1 M☉ ZAMS model (X = 0.70, Z = 0.02).
//! Prints a markdown table of κ_table/κ_analytic at selected depths; no pass
//! value is implied.
//!
//! Run: `cargo run --release -p stellar-core --example opacity_compare`

use stellar_core::physics::Opacity;
use stellar_core::physics::opacity::AnalyticOpacity;
use stellar_core::physics::opacity_table::TabulatedOpacity;
use stellar_core::{NoControl, Physics, StarParams, solve_model};

fn main() {
    let m = solve_model(&StarParams::sun(), &Physics::analytic(), None, &mut NoControl).unwrap();
    let p = &m.profile;
    let tab = TabulatedOpacity { conduction: true };
    println!("| r/R | log T | log ρ | log R | κ analytic | κ tables | ratio | table weight |");
    println!("|---|---|---|---|---|---|---|---|");
    let r_tot = *p.r.last().unwrap();
    let mut last_lt = f64::INFINITY;
    for k in 0..p.len() {
        let lt = p.t[k].log10();
        // One row per 0.25 dex in T, plus the centre and the surface.
        if !(k == 0 || k == p.len() - 1 || lt < last_lt - 0.25) {
            continue;
        }
        last_lt = lt;
        let c = p.composition_at(k, 0.7);
        let ka = AnalyticOpacity.kappa(p.rho[k], p.t[k], &c);
        let parts = tab.parts(p.rho[k], p.t[k], &c);
        let kt = tab.kappa(p.rho[k], p.t[k], &c);
        println!(
            "| {:.4} | {:.2} | {:.2} | {:.2} | {:.3e} | {:.3e} | {:.3} | {:.2} |",
            p.r[k] / r_tot,
            lt,
            p.rho[k].log10(),
            p.rho[k].log10() - 3.0 * (lt - 6.0),
            ka,
            kt,
            kt / ka,
            parts.table_weight
        );
    }
}
