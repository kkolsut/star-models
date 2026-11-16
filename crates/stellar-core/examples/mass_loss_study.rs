//! Main-sequence evolution with a wind (plan 4.2): mass lost, fuel accounting
//! and total-energy check, with default and halved time steps.
//!
//! `cargo run --release -p stellar-core --example mass_loss_study [mass ...]`

use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, Wind,
    solve_model,
};

fn run(mass: f64, opts: EvolutionOptions, n_points: usize) -> (StarModel, usize, String) {
    let physics = Physics::analytic();
    let params = StarParams { mass, n_points, ..StarParams::sun() };
    let zams = solve_model(&params, &physics, None, &mut NoControl).unwrap();
    let mut evo = Evolution::new(zams, opts);
    loop {
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => {}
            EvolutionStep::Finished(reason) => return (evo.model, evo.steps, reason),
        }
    }
}

fn main() {
    let masses: Vec<f64> = std::env::args().skip(1).filter_map(|a| a.parse().ok()).collect();
    let masses = if masses.is_empty() { vec![5.0, 20.0, 60.0] } else { masses };
    for mass in masses {
        let base = EvolutionOptions {
            stop: StopCondition::EndOfMainSequence,
            pulsation: false,
            wind: Wind::standard(),
            ..Default::default()
        };
        let half = EvolutionOptions {
            dx_max: 0.5 * base.dx_max,
            dlog_target: 0.5 * base.dlog_target,
            dm_max: 0.5 * base.dm_max,
            ..base
        };
        let off = EvolutionOptions { wind: Wind::Off, ..base };
        for (name, opts) in [("wind off", off), ("wind", base), ("wind, half steps", half)] {
            let t0 = std::time::Instant::now();
            let (m, steps, reason) = run(mass, opts, 800);
            let g = &m.globals;
            let c = g.energy_check().unwrap();
            println!(
                "{mass:>5} M☉ {name:<17} {steps:>4} steps, {:.3} Myr, M = {:.4} M☉ ({:+.3}%), Ṁ = {:.3e} M☉/yr, \
                 L = {:.4e} L☉, T_eff = {:.0} K; residual {:.3e} (wind {:.3e}, comp {:.3e}, num {:.3e}); fuel {:.5}; {:.1} s; {reason}",
                g.age_yr / 1e6,
                g.mass,
                100.0 * (g.mass / mass - 1.0),
                g.mdot_msun_yr,
                g.luminosity,
                g.teff,
                c.residual,
                c.wind,
                c.composition,
                c.numerical,
                c.fuel_ratio,
                t0.elapsed().as_secs_f64(),
            );
        }
    }
}
