//! Plan 7.8: provisional solar calibration. Adjusts the initial helium Y₀ and
//! α_MLT so that 1 M☉ at 4.57 Gyr has L☉ and R☉, with the surface Z/X equal to
//! the initial one (no diffusion) and fixed at the Grevesse & Noels (1993)
//! photospheric value (Z/X)☉ = 0.0245, the mixture of the opacity tables.
//!
//! Newton's method in (Y₀, α) with a finite-difference Jacobian on the first
//! iteration and Broyden updates afterwards; each evaluation is a ZAMS model
//! and an evolution to 4.57 Gyr with the default time-step controls.
//!
//! Run: `cargo run --release -p stellar-core --example solar_calibration [physics] [points] [dlog] [dx]`
//! with physics `all` (default; every stage-4 upgrade) or `analytic`; `dlog`
//! and `dx` are the time-step targets of the evolution (defaults 0.02, 0.01).

use stellar_core::constants::{L_SUN, R_SUN};
use stellar_core::mixture::GN93;
use stellar_core::physics::options::PhysicsOptions;
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

/// Grevesse & Noels (1993) photospheric (Z/X)☉.
const ZX_SUN: f64 = 0.0245;
const AGE: f64 = 4.57e9;

fn params(y: f64, alpha: f64, n: usize) -> StarParams {
    let x = (1.0 - y) / (1.0 + ZX_SUN);
    StarParams { mass: 1.0, x, z: ZX_SUN * x, alpha_mlt: alpha, f_cno: GN93.f_cno, n_points: n, alpha_ov: 0.0 }
}

/// (ln L/L☉, ln R/R☉) at 4.57 Gyr, and the model.
fn evaluate(
    p: &StarParams,
    physics: &Physics,
    warm: Option<&StarModel>,
    steps: (f64, f64),
) -> ([f64; 2], StarModel, StarModel) {
    let t0 = std::time::Instant::now();
    let zams = solve_model(p, physics, warm, &mut NoControl).expect("ZAMS");
    let opts = EvolutionOptions {
        stop: StopCondition::Age(AGE),
        pulsation: false,
        dlog_target: steps.0,
        dx_max: steps.1,
        ..Default::default()
    };
    let mut evo = Evolution::new(zams.clone(), opts);
    loop {
        match evo.step(physics, &mut NoControl).expect("evolution") {
            EvolutionStep::Advanced { .. } => {}
            EvolutionStep::Finished(r) => {
                assert!((evo.model.globals.age_yr / AGE - 1.0).abs() < 1e-6, "stopped early: {r}");
                break;
            }
        }
    }
    let g = &evo.model.globals;
    eprintln!("  ({} steps, {:.0} s)", evo.steps, t0.elapsed().as_secs_f64());
    ([g.luminosity.ln(), g.radius.ln()], zams, evo.model)
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let opts = match args.first().map(String::as_str) {
        Some("analytic") => PhysicsOptions::analytic(),
        _ => PhysicsOptions::all_upgrades(),
    };
    let n: usize = args.get(1).map(|s| s.parse().unwrap()).unwrap_or(800);
    let steps: (f64, f64) = (
        args.get(2).map(|s| s.parse().unwrap()).unwrap_or(0.02),
        args.get(3).map(|s| s.parse().unwrap()).unwrap_or(0.01),
    );
    let physics = Physics::from_options(&opts);
    eprintln!("physics: {}", physics.describe());
    let mut u = [0.27_f64, 1.9_f64]; // (Y₀, α)
    let (mut f, mut zams, mut last) = evaluate(&params(u[0], u[1], n), &physics, None, steps);
    eprintln!("Y = {:.5}, α = {:.4}: ln L = {:+.5}, ln R = {:+.5}", u[0], u[1], f[0], f[1]);
    // Finite-difference Jacobian.
    let h = [0.005, 0.05];
    let mut jac = [[0.0; 2]; 2];
    for j in 0..2 {
        let mut v = u;
        v[j] += h[j];
        let (fj, _, _) = evaluate(&params(v[0], v[1], n), &physics, Some(&zams), steps);
        for i in 0..2 {
            jac[i][j] = (fj[i] - f[i]) / h[j];
        }
        eprintln!("  ∂/∂{}: {:+.4} {:+.4}", ["Y", "α"][j], jac[0][j], jac[1][j]);
    }
    for it in 0..8 {
        if f[0].abs() < 2e-4 && f[1].abs() < 2e-4 {
            break;
        }
        let det = jac[0][0] * jac[1][1] - jac[0][1] * jac[1][0];
        let du = [(-f[0] * jac[1][1] + f[1] * jac[0][1]) / det, (-jac[0][0] * f[1] + jac[1][0] * f[0]) / det];
        let du = [du[0].clamp(-0.03, 0.03), du[1].clamp(-0.5, 0.5)];
        let u_new = [u[0] + du[0], u[1] + du[1]];
        let (f_new, z_new, m_new) = evaluate(&params(u_new[0], u_new[1], n), &physics, Some(&zams), steps);
        // Broyden update.
        let df = [f_new[0] - f[0], f_new[1] - f[1]];
        let jdu = [jac[0][0] * du[0] + jac[0][1] * du[1], jac[1][0] * du[0] + jac[1][1] * du[1]];
        let nn = du[0] * du[0] + du[1] * du[1];
        for i in 0..2 {
            for j in 0..2 {
                jac[i][j] += (df[i] - jdu[i]) * du[j] / nn;
            }
        }
        u = u_new;
        f = f_new;
        zams = z_new;
        last = m_new;
        eprintln!("it {it}: Y = {:.5}, α = {:.4}: ln L = {:+.2e}, ln R = {:+.2e}", u[0], u[1], f[0], f[1]);
    }
    let p = params(u[0], u[1], n);
    let g = &last.globals;
    println!(
        "Provisional solar calibration ({}, {} points, time-step targets dlog = {}, dX = {}, no diffusion):\n",
        physics.describe(),
        n,
        steps.0,
        steps.1
    );
    println!(
        "| Y₀ | X₀ | Z₀ | α_MLT | L [L☉] | R [R☉] | T_eff [K] | envelope base [r/R] | T_c [K] | X_c | ZAMS L [L☉] |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|---|");
    println!(
        "| {:.4} | {:.4} | {:.5} | {:.3} | {:.4} | {:.4} | {:.0} | {} | {:.4e} | {:.4} | {:.4} |",
        p.y(),
        p.x,
        p.z,
        p.alpha_mlt,
        g.luminosity,
        g.radius,
        g.teff,
        g.conv_env_base_r.map(|r| format!("{r:.4}")).unwrap_or_else(|| "–".into()),
        g.tc,
        g.x_c,
        zams.globals.luminosity
    );
    let _ = (L_SUN, R_SUN);
}
