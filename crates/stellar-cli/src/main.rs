//! Headless ZAMS model builder.
//!
//! ```text
//! stellar-cli --mass 1 --x 0.70 --z 0.02 -o sun.csv --json sun.json
//! stellar-cli --sweep 0.5,1,2,5,10
//! ```

use std::time::Instant;

use anyhow::{Context, Result, bail};
use clap::Parser;
use stellar_core::output::{profile_csv, summary_json};
use stellar_core::{Physics, SolveControl, SolveEvent, StarModel, StarParams, solve_model};

#[derive(Parser, Debug)]
#[command(version, about = "Build a zero-age main-sequence stellar model")]
struct Args {
    /// Stellar mass [M☉]
    #[arg(long, default_value_t = 1.0)]
    mass: f64,
    /// Hydrogen mass fraction X
    #[arg(long, default_value_t = 0.70)]
    x: f64,
    /// Metallicity Z (Y = 1 − X − Z)
    #[arg(long, default_value_t = 0.02)]
    z: f64,
    /// Mixing-length parameter α_MLT
    #[arg(long, default_value_t = 1.8)]
    alpha: f64,
    /// Fraction of Z in C, N and O
    #[arg(long, default_value_t = 0.7)]
    f_cno: f64,
    /// Number of mesh points
    #[arg(long, default_value_t = 800)]
    points: usize,
    /// Write the radial profile to this CSV file
    #[arg(short, long)]
    output: Option<String>,
    /// Write parameters + global properties to this JSON file
    #[arg(long)]
    json: Option<String>,
    /// Comma-separated list of masses: solve a ZAMS sequence by continuation
    #[arg(long, value_delimiter = ',')]
    sweep: Option<Vec<f64>>,
    /// Print every Newton iteration
    #[arg(short, long)]
    verbose: bool,
}

struct Printer {
    verbose: bool,
}

impl SolveControl for Printer {
    fn report(&mut self, event: SolveEvent) {
        match event {
            SolveEvent::Stage(s) => eprintln!("» {s}"),
            SolveEvent::Iteration { iter, max_correction, residual, damping } if self.verbose => {
                eprintln!("   it {iter:3}  max|Δy| = {max_correction:9.3e}  max|F| = {residual:9.3e}  λ = {damping:.3}")
            }
            _ => {}
        }
    }
}

fn print_summary(m: &StarModel) {
    let g = &m.globals;
    println!(
        "M = {:.3} M☉   X = {:.4}  Y = {:.4}  Z = {:.4}   α_MLT = {}",
        g.mass,
        m.params.x,
        m.params.y(),
        m.params.z,
        m.params.alpha_mlt
    );
    println!("  R      = {:.4} R☉", g.radius);
    println!("  L      = {:.4} L☉", g.luminosity);
    println!("  T_eff  = {:.0} K   log g = {:.3}", g.teff, g.log_g);
    println!("  T_c    = {:.4e} K   ρ_c = {:.4e} g/cm³   P_c = {:.4e} dyn/cm²", g.tc, g.rhoc, g.pc);
    println!("  energy: pp {:.1}%  CNO {:.1}%  3α {:.2}%", 100.0 * g.frac_pp, 100.0 * g.frac_cno, 100.0 * g.frac_3a);
    if g.conv_core_q > 0.0 {
        println!("  convective core: m/M = {:.4}, r/R = {:.4}", g.conv_core_q, g.conv_core_r);
    } else {
        println!("  radiative core");
    }
    match (g.conv_env_base_r, g.conv_env_base_q) {
        (Some(r), Some(q)) => println!("  convective envelope base: r/R = {r:.4}, m/M = {q:.6}"),
        _ => println!("  no convective envelope"),
    }
    println!(
        "  zones: {}",
        m.zones.iter().map(|z| format!("[r/R {:.4}–{:.4}]", z.r_inner, z.r_outer)).collect::<Vec<_>>().join(" ")
    );
    println!(
        "  checks: virial {:.5}  energy {:.5}  Newton iters {}",
        g.virial_ratio, g.energy_ratio, g.newton_iterations
    );
}

fn main() -> Result<()> {
    let args = Args::parse();
    let physics = Physics::analytic();
    let base = StarParams {
        mass: args.mass,
        x: args.x,
        z: args.z,
        alpha_mlt: args.alpha,
        f_cno: args.f_cno,
        n_points: args.points,
    };
    let mut ctl = Printer { verbose: args.verbose };

    if let Some(masses) = &args.sweep {
        if masses.is_empty() {
            bail!("--sweep needs at least one mass");
        }
        let mut prev: Option<StarModel> = None;
        println!(
            "{:>8} {:>9} {:>11} {:>8} {:>11} {:>9} {:>7} {:>9} {:>9}",
            "M/M☉", "R/R☉", "L/L☉", "Teff", "Tc", "ρc", "CNO%", "q_cc", "r_env"
        );
        for &m in masses {
            let p = StarParams { mass: m, ..base.clone() };
            let t0 = Instant::now();
            let model = solve_model(&p, &physics, prev.as_ref(), &mut ctl).with_context(|| format!("M = {m}"))?;
            let g = &model.globals;
            println!(
                "{:8.3} {:9.4} {:11.4e} {:8.0} {:11.4e} {:9.3} {:7.2} {:9.4} {:>9}   ({:.2}s)",
                g.mass,
                g.radius,
                g.luminosity,
                g.teff,
                g.tc,
                g.rhoc,
                100.0 * g.frac_cno,
                g.conv_core_q,
                g.conv_env_base_r.map(|r| format!("{r:.4}")).unwrap_or_else(|| "-".into()),
                t0.elapsed().as_secs_f64()
            );
            prev = Some(model);
        }
        return Ok(());
    }

    let t0 = Instant::now();
    let model = solve_model(&base, &physics, None, &mut ctl)?;
    eprintln!("solved in {:.2}s", t0.elapsed().as_secs_f64());
    if let Some(w) = base.ideal_gas_warning() {
        eprintln!("warning: {w}");
    }
    print_summary(&model);
    if let Some(path) = &args.output {
        std::fs::write(path, profile_csv(&model)).with_context(|| format!("writing {path}"))?;
        eprintln!("profile written to {path}");
    }
    if let Some(path) = &args.json {
        std::fs::write(path, summary_json(&model)).with_context(|| format!("writing {path}"))?;
        eprintln!("summary written to {path}");
    }
    Ok(())
}
