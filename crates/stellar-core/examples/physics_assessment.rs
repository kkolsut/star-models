//! Plan 7.7: joint assessment of atmosphere, EOS and opacity on the 1 M☉
//! ZAMS model. One ingredient at a time is swapped from the analytic set,
//! then all together, and a table of the quantities named by the plan is
//! printed (markdown). The nuclear rates have a single implementation since
//! plan 1.6 (Solar Fusion III), so they are not swapped.
//!
//! Run: `cargo run --release -p stellar-core --example physics_assessment [mass] [X] [Z] [alpha]`

use stellar_core::physics::options::*;
use stellar_core::{NoControl, Physics, StarModel, StarParams, solve_model};

fn row(name: &str, m: &StarModel, secs: f64) -> String {
    let g = &m.globals;
    let p = &m.profile;
    let n = p.len();
    // Smallest ∇_ad outside the core (the ionisation zones) and where.
    let (k_min, nad_min) = (0..n)
        .filter(|&k| p.q[k] > 0.5)
        .map(|k| (k, p.nabla_ad[k]))
        .fold((0, f64::INFINITY), |a, b| if b.1 < a.1 { b } else { a });
    format!(
        "| {name} | {:.4} | {:.4} | {:.0} | {} | {:.4e} | {:.2} | {:.3} | {:.3} (log T = {:.2}) | {:.3e} | {:.1} |",
        g.luminosity,
        g.radius,
        g.teff,
        g.conv_env_base_r.map(|r| format!("{r:.4}")).unwrap_or_else(|| "–".into()),
        g.tc,
        g.rhoc,
        p.kappa[0],
        nad_min,
        p.t[k_min].log10(),
        g.surface.p_phot,
        secs
    )
}

fn main() {
    let args: Vec<f64> = std::env::args().skip(1).map(|s| s.parse().unwrap()).collect();
    let mass = args.first().copied().unwrap_or(1.0);
    let base_p = StarParams {
        mass,
        x: args.get(1).copied().unwrap_or(0.70),
        z: args.get(2).copied().unwrap_or(0.02),
        alpha_mlt: args.get(3).copied().unwrap_or(1.8),
        ..StarParams::sun()
    };
    let a = PhysicsOptions::analytic();
    let sets = [
        ("analytic (all defaults)", a),
        ("atmosphere: integrated grey", PhysicsOptions { atmosphere: AtmosphereChoice::IntegratedGrey, ..a }),
        ("atmosphere: Krishna Swamy", PhysicsOptions { atmosphere: AtmosphereChoice::KrishnaSwamy, ..a }),
        ("EOS: Fermi–Dirac, fully ionised", PhysicsOptions { eos: EosChoice::Degenerate, ..a }),
        ("EOS: + ionisation, no PI", PhysicsOptions { eos: EosChoice::Ionisation, ..a }),
        ("EOS: + pressure ionisation + H₂", PhysicsOptions { eos: EosChoice::Full, ..a }),
        ("opacity: tables", PhysicsOptions { opacity: OpacityChoice::Tables, ..a }),
        ("all upgrades", PhysicsOptions::all_upgrades()),
    ];
    println!(
        "{mass} M☉ ZAMS, X = {}, Z = {}, α_MLT = {}, {} points.\n",
        base_p.x, base_p.z, base_p.alpha_mlt, base_p.n_points
    );
    println!(
        "| Physics | L [L☉] | R [R☉] | T_eff [K] | envelope base [r/R] | T_c [K] | ρ_c [g cm⁻³] | κ_c [cm² g⁻¹] | min ∇_ad, q > 0.5 | P_phot [dyn cm⁻²] | time [s] |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|---|");
    for (name, o) in sets {
        let phys = Physics::from_options(&o);
        let t0 = std::time::Instant::now();
        match solve_model(&base_p, &phys, None, &mut NoControl) {
            Ok(m) => println!("{}", row(name, &m, t0.elapsed().as_secs_f64())),
            Err(e) => println!("| {name} | failed: {e} |"),
        }
    }
}
