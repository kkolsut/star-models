//! Track V.2: mesh and central-cut-off convergence.
//!
//! Solves the same star on nested meshes of 401, 801, 1601 and 3201 points
//! (uniform in the stretched coordinate, so every point of a coarse mesh is
//! also a point of the finer ones) and compares
//!
//! * global values (L, R, T_eff, T_c, ρ_c, convective boundaries),
//! * local profiles at the shared points: ln r, ln P, ln T, L/L_s, ε and the
//!   hydrogen gradient dX/dq,
//!
//! for ZAMS, mid-main-sequence and end-of-main-sequence models of 1 and 3 M☉.
//! The evolved models are reached with the same fixed sequence of time steps
//! on every mesh, so that only the spatial discretisation differs. The
//! observed order is p = log₂(|e_N| / |e_2N|) with e_N the difference to the
//! 3201-point model.
//!
//! The second part varies the central cut-off m₀ = q₀ M from 10⁻⁸ to 10⁻⁴ at
//! 1601 points.
//!
//! Run: `cargo run --release -p stellar-core --example mesh_convergence`

use stellar_core::solver::NewtonOptions;
use stellar_core::structure::initial::polytrope_guess;
use stellar_core::structure::mesh::Q_INNER;
use stellar_core::structure::mesh::make_mesh_from;
use stellar_core::structure::solve::newton;
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition,
};

const MESHES: [usize; 4] = [401, 801, 1601, 3201];

fn zams(mass: f64, n: usize, q0: f64, physics: &Physics) -> StarModel {
    let params = StarParams { mass, n_points: n, ..StarParams::sun() };
    let comp = params.composition().unwrap();
    let q = make_mesh_from(n, q0);
    let guess = polytrope_guess(physics, mass, &comp, &q, params.alpha_mlt);
    newton(&params, physics, &q, guess.y, guess.l_scale, NewtonOptions { max_iter: 250, tol: 1e-8 }, &mut NoControl)
        .unwrap_or_else(|e| panic!("M = {mass}, N = {n}, q0 = {q0:e}: {e}"))
}

fn opts(end: f64) -> EvolutionOptions {
    EvolutionOptions {
        stop: StopCondition::Age(end),
        dx_max: 1.0,
        reject_factor: 1e9,
        pulsation: false,
        ..Default::default()
    }
}

/// Time-step schedule: the accepted steps of a run on the coarsest mesh that
/// asks for `dt_yr` (the solver's 30% burning limit may shorten a step) and
/// stops exactly at each of `ages_yr`.
fn schedule(zams: StarModel, dt_yr: f64, ages_yr: &[f64], physics: &Physics) -> Vec<f64> {
    let mut evo = Evolution::new(zams, opts(*ages_yr.last().unwrap()));
    let mut out = Vec::new();
    for &age in ages_yr {
        while evo.model.globals.age_yr < age * (1.0 - 1e-12) {
            evo.dt_yr = dt_yr.min(age - evo.model.globals.age_yr);
            match evo.step(physics, &mut NoControl).unwrap() {
                EvolutionStep::Advanced { dt_yr: dt, .. } => out.push(dt),
                EvolutionStep::Finished(r) => panic!("finished early: {r}"),
            }
        }
    }
    out
}

/// Evolve with the given time steps, returning the models at `ages_yr` and
/// the number of steps the burning limit split in two (identical schedules
/// otherwise).
fn evolve(zams: StarModel, steps: &[f64], ages_yr: &[f64], physics: &Physics) -> (Vec<StarModel>, usize) {
    let mut evo = Evolution::new(zams, opts(steps.iter().sum::<f64>() * 2.0));
    let mut out = Vec::new();
    let mut next = 0;
    let mut split = 0;
    let mut target = 0.0;
    for &dt in steps {
        target += dt;
        let mut first = true;
        while evo.model.globals.age_yr < target * (1.0 - 1e-12) {
            evo.dt_yr = target - evo.model.globals.age_yr;
            match evo.step(physics, &mut NoControl).unwrap() {
                EvolutionStep::Advanced { .. } => {}
                EvolutionStep::Finished(r) => panic!("finished early: {r}"),
            }
            if !first {
                split += 1;
            }
            first = false;
        }
        if next < ages_yr.len() && (evo.model.globals.age_yr / ages_yr[next] - 1.0).abs() < 1e-9 {
            out.push(evo.model.clone());
            next += 1;
        }
    }
    assert_eq!(out.len(), ages_yr.len(), "schedule does not hit the output ages");
    (out, split)
}

fn globals(m: &StarModel) -> Vec<(&'static str, f64)> {
    let g = &m.globals;
    vec![
        ("L", g.luminosity),
        ("R", g.radius),
        ("T_eff", g.teff),
        ("T_c", g.tc),
        ("rho_c", g.rhoc),
        ("X_c", g.x_c),
        ("q_cc", g.conv_core_q),
        ("r_env", g.conv_env_base_r.unwrap_or(f64::NAN)),
    ]
}

/// Local quantities at every point of the coarse mesh (index stride `s` in the fine one).
fn locals(m: &StarModel, stride: usize) -> Vec<(&'static str, Vec<f64>)> {
    let p = &m.profile;
    let idx: Vec<usize> = (0..p.len()).step_by(stride).collect();
    let pick = |v: &[f64], f: &dyn Fn(f64) -> f64| idx.iter().map(|&k| f(v[k])).collect::<Vec<_>>();
    let l_s = p.l[p.len() - 1];
    // dX/dq between consecutive coarse points.
    let dxdq: Vec<f64> = idx.windows(2).map(|w| (p.x[w[1]] - p.x[w[0]]) / (p.q[w[1]] - p.q[w[0]])).collect();
    vec![
        ("ln r", pick(&p.r, &|v| v.ln())),
        ("ln P", pick(&p.p, &|v| v.ln())),
        ("ln T", pick(&p.t, &|v| v.ln())),
        ("L/L_s", pick(&p.l, &|v| v / l_s)),
        ("eps/eps_c", {
            let e0 = p.eps[0].max(1e-300);
            pick(&p.eps, &|v| v / e0)
        }),
        ("dX/dq", dxdq),
    ]
}

fn report(label: &str, models: &[StarModel]) {
    println!("\n### {label}\n");
    let fine = models.last().unwrap();
    // Globals.
    println!("| quantity | {} |", MESHES.iter().map(|n| format!("N = {n}")).collect::<Vec<_>>().join(" | "));
    println!("|---|{}", "---|".repeat(MESHES.len()));
    let gs: Vec<_> = models.iter().map(globals).collect();
    for (i, (name, _)) in gs[0].iter().enumerate() {
        let vals: Vec<String> = gs.iter().map(|g| format!("{:.6e}", g[i].1)).collect();
        println!("| {name} | {} |", vals.join(" | "));
    }
    println!("\nRelative difference to N = {} and observed order:\n", MESHES[3]);
    println!("| quantity | e(401) | e(801) | e(1601) | p(401/801) | p(801/1601) |");
    println!("|---|---|---|---|---|---|");
    for (i, (name, f)) in gs[3].iter().enumerate() {
        let e: Vec<f64> = (0..3).map(|j| (gs[j][i].1 - f) / f.abs().max(1e-300)).collect();
        let p = |a: f64, b: f64| if b.abs() > 0.0 && a.abs() > 0.0 { (a / b).abs().log2() } else { f64::NAN };
        println!(
            "| {name} | {:.2e} | {:.2e} | {:.2e} | {:.2} | {:.2} |",
            e[0],
            e[1],
            e[2],
            p(e[0], e[1]),
            p(e[1], e[2])
        );
    }
    // Local profiles on the 401-point mesh.
    println!(
        "\nLocal profiles: largest |difference| to N = {} over the 401 shared points (ε relative to its",
        MESHES[3]
    );
    println!("central value; dX/dq absolute), and observed order:\n");
    println!("| profile | e(401) | e(801) | e(1601) | p(401/801) | p(801/1601) |");
    println!("|---|---|---|---|---|---|");
    let strides = [1, 2, 4, 8];
    let lf = locals(fine, strides[3]);
    let ls: Vec<_> = (0..3).map(|j| locals(&models[j], strides[j])).collect();
    for (i, (name, fv)) in lf.iter().enumerate() {
        let e: Vec<f64> =
            (0..3).map(|j| ls[j][i].1.iter().zip(fv).map(|(a, b)| (a - b).abs()).fold(0.0, f64::max)).collect();
        let p = |a: f64, b: f64| if b > 0.0 && a > 0.0 { (a / b).log2() } else { f64::NAN };
        println!(
            "| {name} | {:.2e} | {:.2e} | {:.2e} | {:.2} | {:.2} |",
            e[0],
            e[1],
            e[2],
            p(e[0], e[1]),
            p(e[1], e[2])
        );
    }
    // Nesting check.
    let q_c = &models[0].profile.q;
    let q_f = &fine.profile.q;
    let nest = q_c.iter().enumerate().map(|(k, q)| (q - q_f[8 * k]).abs() / q).fold(0.0, f64::max);
    println!("\n(nesting: largest relative q mismatch {nest:.1e})");
}

fn main() {
    let physics = Physics::analytic();
    // (mass, fixed time step, mid-MS age, near-end-of-MS age)
    let cases = [(1.0, 2.0e8, 8.0e9, 1.5e10), (3.0, 5.0e6, 1.3e8, 2.55e8)];
    for &(mass, dt, mid, end) in &cases {
        let mut zs = Vec::new();
        let mut mids = Vec::new();
        let mut ends = Vec::new();
        let steps = schedule(zams(mass, MESHES[0], Q_INNER, &physics), dt, &[mid, end], &physics);
        eprintln!("M = {mass}: {} steps", steps.len());
        for &n in &MESHES {
            let t0 = std::time::Instant::now();
            let z = zams(mass, n, Q_INNER, &physics);
            let (ev, split) = evolve(z.clone(), &steps, &[mid, end], &physics);
            eprintln!(
                "M = {mass}, N = {n}: {:.1} s, {split} steps split by the burning limit",
                t0.elapsed().as_secs_f64()
            );
            zs.push(z);
            mids.push(ev[0].clone());
            ends.push(ev[1].clone());
        }
        report(&format!("{mass} M☉ ZAMS"), &zs);
        report(&format!("{mass} M☉ at {:.3} Gyr (X_c = {:.3})", mid / 1e9, mids[3].globals.x_c), &mids);
        report(&format!("{mass} M☉ at {:.3} Gyr (X_c = {:.4})", end / 1e9, ends[3].globals.x_c), &ends);
    }

    // Central cut-off.
    println!("\n### Central cut-off m₀ = q₀ M, 1601 points\n");
    for &mass in &[1.0, 3.0] {
        let reference = zams(mass, 1601, 1e-8, &physics);
        let fr = globals(&reference);
        println!("\n{mass} M☉ ZAMS, relative difference to q₀ = 10⁻⁸:\n");
        println!("| q₀ | {} |", fr.iter().map(|(n, _)| *n).collect::<Vec<_>>().join(" | "));
        println!("|---|{}", "---|".repeat(fr.len()));
        for q0 in [1e-7, 1e-6, 1e-5, 1e-4] {
            let m = zams(mass, 1601, q0, &physics);
            let g = globals(&m);
            let row: Vec<String> =
                g.iter().zip(&fr).map(|((_, a), (_, b))| format!("{:.1e}", (a - b) / b.abs().max(1e-300))).collect();
            println!("| {q0:.0e} | {} |", row.join(" | "));
        }
    }
}
