//! Time-step convergence study for stellar evolution.
//!
//! Runs the same tracks with progressively tighter time-step controls and
//! prints the quantities that should converge: main-sequence lifetime, the
//! state at fixed central hydrogen, and the present-day Sun.
//!
//! ```text
//! cargo run --release -p stellar-core --example timestep_study
//! ```

use std::time::Instant;

use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

struct Setting {
    name: &'static str,
    dlog_target: f64,
    dx_max: f64,
    reject_factor: f64,
}

const SETTINGS: &[Setting] = &[
    Setting { name: "default 2% / ΔX 1%", dlog_target: 0.02, dx_max: 0.01, reject_factor: 4.0 },
    Setting { name: "1% / ΔX 0.5%", dlog_target: 0.01, dx_max: 0.005, reject_factor: 4.0 },
    Setting { name: "0.5% / ΔX 0.25%", dlog_target: 0.005, dx_max: 0.0025, reject_factor: 4.0 },
    Setting { name: "0.25% / ΔX 0.125%", dlog_target: 0.0025, dx_max: 0.00125, reject_factor: 4.0 },
    Setting { name: "2% / ΔX 1%, reject 1.5×", dlog_target: 0.02, dx_max: 0.01, reject_factor: 1.5 },
];

fn run(zams: &StarModel, physics: &Physics, s: &Setting, stop: StopCondition) -> (Vec<StarModel>, f64) {
    let opts = EvolutionOptions {
        stop,
        dlog_target: s.dlog_target,
        dx_max: s.dx_max,
        reject_factor: s.reject_factor,
        pulsation: false,
        max_steps: 100_000,
        ..Default::default()
    };
    let t0 = Instant::now();
    let mut evo = Evolution::new(zams.clone(), opts);
    let mut models = vec![zams.clone()];
    loop {
        match evo.step(physics, &mut NoControl).expect("evolution step") {
            EvolutionStep::Advanced { .. } => models.push(evo.model.clone()),
            EvolutionStep::Finished(_) => return (models, t0.elapsed().as_secs_f64()),
        }
    }
}

/// Linear interpolation of `f` at the age where X_c first drops to `xc`.
fn at_xc(models: &[StarModel], xc: f64, f: impl Fn(&StarModel) -> f64) -> f64 {
    for w in models.windows(2) {
        let (a, b) = (w[0].globals.x_c, w[1].globals.x_c);
        if a >= xc && b < xc {
            let t = (a - xc) / (a - b);
            return f(&w[0]) + t * (f(&w[1]) - f(&w[0]));
        }
    }
    f64::NAN
}

fn main() {
    let physics = Physics::analytic();

    println!("## 1 M☉ to 4.57 Gyr");
    println!("| setting | steps | L [L☉] | R [R☉] | T_c [K] | X_c | time [s] |");
    println!("|---|---:|---:|---:|---:|---:|---:|");
    let sun = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).unwrap();
    for s in SETTINGS {
        let (m, secs) = run(&sun, &physics, s, StopCondition::Age(4.57e9));
        let g = &m.last().unwrap().globals;
        println!(
            "| {} | {} | {:.5} | {:.5} | {:.5e} | {:.5} | {:.1} |",
            s.name,
            m.len() - 1,
            g.luminosity,
            g.radius,
            g.tc,
            g.x_c,
            secs
        );
    }

    for mass in [1.0, 3.0] {
        println!("\n## {mass} M☉ to the end of the main sequence");
        println!(
            "| setting | steps | MS lifetime [Myr] | q_cc at X_c=0.3 | L at X_c=0.3 | T_eff at X_c=0.3 | L at TAMS | T_eff at TAMS | time [s] |"
        );
        println!("|---|---:|---:|---:|---:|---:|---:|---:|---:|");
        let zams = solve_model(&StarParams { mass, ..StarParams::sun() }, &physics, None, &mut NoControl).unwrap();
        for s in SETTINGS {
            let (m, secs) = run(&zams, &physics, s, StopCondition::EndOfMainSequence);
            let end = &m.last().unwrap().globals;
            println!(
                "| {} | {} | {:.2} | {:.5} | {:.5} | {:.1} | {:.5} | {:.1} | {:.1} |",
                s.name,
                m.len() - 1,
                end.age_yr / 1e6,
                at_xc(&m, 0.3, |m| m.globals.conv_core_q),
                at_xc(&m, 0.3, |m| m.globals.luminosity),
                at_xc(&m, 0.3, |m| m.globals.teff),
                end.luminosity,
                end.teff,
                secs
            );
        }
    }
}
