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
    /// Helium follows metallicity: Y = 0.249 + 1.5 Z, X = 1 − Y − Z (MIST, Choi et al. 2016); --x is ignored
    #[arg(long)]
    helium_from_z: bool,
    /// Mixing-length parameter α_MLT
    #[arg(long, default_value_t = 1.8)]
    alpha: f64,
    /// Fraction of Z in C, N and O
    #[arg(long, default_value_t = 0.7)]
    f_cno: f64,
    /// Number of mesh points
    #[arg(long, default_value_t = 800)]
    points: usize,
    /// Step overshoot of a convective core, α_ov in units of min(H_P, r_cc) (0–0.5; mixing only, used by --evolve)
    #[arg(long, default_value_t = 0.0)]
    alpha_ov: f64,
    /// Write the radial profile to this CSV file
    #[arg(short, long)]
    output: Option<String>,
    /// Write parameters + global properties to this JSON file
    #[arg(long)]
    json: Option<String>,
    /// Comma-separated list of masses: solve a ZAMS sequence by continuation
    #[arg(long, value_delimiter = ',')]
    sweep: Option<Vec<f64>>,
    /// Evolve the star from the ZAMS (see --until and --max-age)
    #[arg(long)]
    evolve: bool,
    /// Stop condition for --evolve: tams (end of main sequence), age (needs --max-age) or max (as far as possible)
    #[arg(long, default_value = "tams")]
    until: String,
    /// Age to stop at [Gyr] (with --until age), or safety limit otherwise
    #[arg(long)]
    max_age: Option<f64>,
    /// Write the evolutionary track to this CSV file (with --evolve)
    #[arg(long)]
    track: Option<String>,
    /// Print every Newton iteration
    #[arg(short, long)]
    verbose: bool,
}

/// Primordial helium Y_p and enrichment ΔY/ΔZ of the "helium follows
/// metallicity" mode (MIST, Choi et al. 2016, §3.1; docs/physics/10).
const Y_PRIMORDIAL: f64 = 0.249;
const DY_DZ: f64 = 1.5;

/// X for a given Z with Y = Y_p + (ΔY/ΔZ) Z.
fn x_from_z(z: f64) -> f64 {
    1.0 - (Y_PRIMORDIAL + DY_DZ * z) - z
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
    println!("  gravitational redshift v_z = GM/(Rc) = {:.0} m/s", g.grav_redshift_m_s);
    println!("  energy: pp {:.1}%  CNO {:.1}%  3α {:.2}%", 100.0 * g.frac_pp, 100.0 * g.frac_cno, 100.0 * g.frac_3a);
    if g.conv_core_q > 0.0 {
        println!("  convective core: m/M = {:.4}, r/R = {:.4}", g.conv_core_q, g.conv_core_r);
        if m.params.alpha_ov > 0.0 {
            println!(
                "  mixed core (α_ov = {}): m/M = {:.4}, r/R = {:.4}",
                m.params.alpha_ov, g.mixed_core_q, g.mixed_core_r
            );
        }
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
        "  seismic: Δν = {:.2} μHz (asymptotic), {:.2} μHz (scaling); ν_max = {:.0} μHz (scaling)",
        g.delta_nu_asym_uhz, g.delta_nu_uhz, g.nu_max_uhz
    );
    if g.g_cavities.is_empty() {
        println!("  candidate g-mode cavities: none");
    }
    for c in &g.g_cavities {
        println!(
            "  candidate g-mode cavity ({}): r/R {:.4}–{:.4}, ΔΠ₁ = {:.1} s (asymptotic estimate)",
            c.label(),
            c.r_inner,
            c.r_outer,
            c.delta_pi1_s
        );
    }
    let s = &g.screening;
    if !s.reaction.is_empty() {
        println!("  screening: largest exponent {:.3} ({}, m/M = {:.3})", s.max_exponent, s.reaction, s.q);
    }
    if let Some(w) = s.warning() {
        println!("  warning: {w}");
    }
    if let Some(mode) = &m.pulsation {
        println!(
            "  fundamental radial mode: P = {:.4} h, ω/ω_ac = {:.3} at the photosphere",
            mode.period_s / 3600.0,
            mode.omega_over_ac
        );
        if let Some(w) = mode.cutoff_warning() {
            println!("  warning: {w}");
        }
    }
    if let Some(c) = g.energy_check() {
        println!(
            "  run: total-energy residual {:.2e} (composition term {:.2e}, numerical {:.2e}); fuel accounting {:.5}",
            c.residual, c.composition, c.numerical, c.fuel_ratio
        );
        if let Some(run) = &g.run {
            println!("  run: largest screening exponent {:.3}", run.max_screening);
        }
    }
    let r = &g.residuals;
    println!(
        "  checks: virial {:.5}  energy {:.5} (discretisation)  Newton iters {}  max|F| {:.1e}",
        g.virial_ratio,
        g.energy_ratio,
        g.newton_iterations,
        r.max()
    );
    println!(
        "  residuals: continuity {:.1e}  hydrostatic {:.1e}  energy {:.1e}  transport {:.1e}  boundary {:.1e}",
        r.continuity, r.hydrostatic, r.energy, r.transport, r.boundary
    );
    println!("  physics: {}; MLT α = {}", m.physics, m.params.alpha_mlt);
}

fn main() -> Result<()> {
    let args = Args::parse();
    let physics = Physics::analytic();
    let x = if args.helium_from_z { x_from_z(args.z) } else { args.x };
    let base = StarParams {
        mass: args.mass,
        x,
        z: args.z,
        alpha_mlt: args.alpha,
        f_cno: args.f_cno,
        n_points: args.points,
        alpha_ov: args.alpha_ov,
    };
    base.validate().map_err(anyhow::Error::msg)?;
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

    if args.evolve {
        return evolve(&args, &base, &physics, &mut ctl);
    }

    let t0 = Instant::now();
    let mut model = solve_model(&base, &physics, None, &mut ctl)?;
    eprintln!("solved in {:.2}s", t0.elapsed().as_secs_f64());
    model.pulsation = stellar_core::pulsation::fundamental_mode(&model, None).ok();
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

fn evolve(args: &Args, base: &StarParams, physics: &Physics, ctl: &mut Printer) -> Result<()> {
    use stellar_core::{Evolution, EvolutionOptions, EvolutionStep, StopCondition, fmt_age};
    let stop = match args.until.as_str() {
        "tams" => StopCondition::EndOfMainSequence,
        "age" => StopCondition::Age(args.max_age.context("--until age needs --max-age <Gyr>")? * 1e9),
        "max" => StopCondition::AsFarAsPossible,
        other => bail!("unknown --until {other} (use tams, age or max)"),
    };
    let mut opts = EvolutionOptions { stop, ..Default::default() };
    if let (Some(a), false) = (args.max_age, matches!(stop, StopCondition::Age(_))) {
        opts.max_age_yr = a * 1e9;
    }
    let t0 = Instant::now();
    let zams = solve_model(base, physics, None, ctl)?;
    let mut evo = Evolution::new(zams.clone(), opts);
    let mut track = vec![zams];
    println!(
        "{:>5} {:>11} {:>10} {:>10} {:>7} {:>8} {:>10} {:>8} {:>7} {:>9} {:>11} {:>8} {:>9}  phase",
        "step", "age", "dt [yr]", "L/L☉", "Teff", "R/R☉", "Tc", "Xc", "q_He", "P_fund", "pulsator", "energy", "E-resid"
    );
    let print = |step: usize, m: &StarModel, dt: f64| {
        let g = &m.globals;
        let period = m.pulsation.as_ref().map(|p| format!("{:.4}d", p.period_days())).unwrap_or_else(|| "-".into());
        let puls = stellar_core::pulsation::classify(g.teff, g.luminosity)
            .map(|p| match p {
                stellar_core::pulsation::PulsatorType::DeltaScuti => "δ Sct",
                stellar_core::pulsation::PulsatorType::Cepheid => "Cepheid",
            })
            .unwrap_or("");
        println!(
            "{:5} {:>11} {:10.3e} {:10.4e} {:7.0} {:8.4} {:10.4e} {:8.5} {:7.4} {:>9} {:>11} {:8.5} {:>9}  {}",
            step,
            fmt_age(g.age_yr),
            dt,
            g.luminosity,
            g.teff,
            g.radius,
            g.tc,
            g.x_c,
            g.he_core_q,
            period,
            puls,
            g.energy_ratio,
            g.energy_check().map(|c| format!("{:.2e}", c.residual)).unwrap_or_else(|| "-".into()),
            g.phase.label()
        );
    };
    print(0, &track[0], 0.0);
    loop {
        match evo.step(physics, ctl)? {
            EvolutionStep::Advanced { dt_yr, .. } => {
                if args.verbose || evo.steps.is_multiple_of(10) {
                    print(evo.steps, &evo.model, dt_yr);
                }
                track.push(evo.model.clone());
            }
            EvolutionStep::Finished(reason) => {
                print(evo.steps, &evo.model, 0.0);
                eprintln!("finished after {} steps in {:.1}s: {reason}", evo.steps, t0.elapsed().as_secs_f64());
                println!();
                print_summary(&evo.model);
                break;
            }
        }
    }
    if let Some(path) = &args.track {
        let refs: Vec<&StarModel> = track.iter().collect();
        std::fs::write(path, stellar_core::output::track_csv(&refs)).with_context(|| format!("writing {path}"))?;
        eprintln!("track written to {path}");
    }
    if let Some(path) = &args.output {
        std::fs::write(path, profile_csv(&evo.model)).with_context(|| format!("writing {path}"))?;
        eprintln!("final profile written to {path}");
    }
    Ok(())
}
