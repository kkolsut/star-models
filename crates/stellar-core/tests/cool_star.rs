//! Stage 4 (cool-star block): solver-level tests of the integrated atmosphere
//! (plans 3.1–3.3), the opacity tables (7.2–7.4) and the Fermi–Dirac EOS
//! (6.2–6.5). Unit tests of each ingredient live next to its code.
//!
//! The tests that solve models with a Fermi–Dirac EOS are marked `#[ignore]`
//! because they take minutes; run them with `cargo test --release -- --ignored`.

use stellar_core::physics::options::{AtmosphereChoice, EosChoice, OpacityChoice, PhysicsOptions};
use stellar_core::solver::BvpProblem;
use stellar_core::structure::equations::{I_L, I_LNP, I_LNR, I_LNT, NVAR, StructureProblem};
use stellar_core::structure::mesh::make_mesh;
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

fn options(atmosphere: AtmosphereChoice, eos: EosChoice, opacity: OpacityChoice) -> PhysicsOptions {
    PhysicsOptions { eos, opacity, atmosphere }
}

fn sequence(physics: &Physics, masses: &[f64]) -> Vec<StarModel> {
    let mut out: Vec<StarModel> = Vec::new();
    for &m in masses {
        let p = StarParams { mass: m, ..StarParams::sun() };
        let model = solve_model(&p, physics, out.last(), &mut NoControl).unwrap_or_else(|e| panic!("M = {m}: {e}"));
        out.push(model);
    }
    out
}

/// Plan 3.1: the ZAMS sequence 0.4–20 M☉ converges with the integrated
/// atmosphere, and hot stars change by less than 0.1% against the closed form.
#[test]
fn integrated_atmosphere_sequence() {
    let masses = [1.0, 0.7, 0.4, 2.0, 5.0, 10.0, 20.0];
    let grey = sequence(&Physics::analytic(), &masses);
    for atm in [AtmosphereChoice::IntegratedGrey, AtmosphereChoice::KrishnaSwamy] {
        let phys = Physics::from_options(&options(atm, EosChoice::Ideal, OpacityChoice::Analytic));
        let models = sequence(&phys, &masses);
        for (a, b) in models.iter().zip(&grey) {
            let (ga, gb) = (&a.globals, &b.globals);
            eprintln!(
                "{atm:?} M = {}: L {:.5} vs {:.5}, R {:.5} vs {:.5}, T_eff {:.0} vs {:.0}",
                ga.mass, ga.luminosity, gb.luminosity, ga.radius, gb.radius, ga.teff, gb.teff
            );
            if gb.teff > 10_000.0 {
                for (name, x, y) in
                    [("L", ga.luminosity, gb.luminosity), ("R", ga.radius, gb.radius), ("T_eff", ga.teff, gb.teff)]
                {
                    assert!((x / y - 1.0).abs() < 1e-3, "{atm:?} M = {}: {name} {x} vs {y}", ga.mass);
                }
            }
            // T_eff is derived from L and R.
            let teff = stellar_core::physics::atmosphere::effective_temperature(
                ga.luminosity * stellar_core::constants::L_SUN,
                ga.radius * stellar_core::constants::R_SUN,
            );
            assert!((teff / ga.teff - 1.0).abs() < 1e-10);
        }
    }
}

/// Plan 3.1: the derivatives of the outer boundary conditions are smooth:
/// central differences with h and h/10 agree.
#[test]
fn boundary_derivatives_are_smooth() {
    let phys = Physics::from_options(&options(AtmosphereChoice::KrishnaSwamy, EosChoice::Ideal, OpacityChoice::Tables));
    let m = solve_model(&StarParams::sun(), &phys, None, &mut NoControl).unwrap();
    let q = make_mesh(m.params.n_points);
    let comp = m.params.composition().unwrap();
    let problem = StructureProblem::new(
        &phys,
        m.globals.mass * stellar_core::constants::M_SUN,
        &q,
        comp,
        m.params.alpha_mlt,
        m.l_scale,
    );
    let n = q.len();
    let yn: Vec<f64> = m.solution[(n - 1) * NVAR..].to_vec();
    // Plan 3.2: ∇ is continuous across the fitting point (interior A = 1
    // against the atmosphere's A = 1 + q'(τ_fit), |q'| < 1e-3).
    let atm = problem.atmosphere_at(&yn).unwrap();
    let (inside, outside) = (m.profile.nabla[n - 1], *atm.nabla.last().unwrap());
    eprintln!("∇ at the fitting point: interior {inside:.6}, atmosphere {outside:.6}");
    assert!((inside / outside - 1.0).abs() < 2e-3, "∇ jumps at the fitting point: {inside} vs {outside}");
    for var in [I_LNR, I_LNP, I_LNT, I_L] {
        let d = |h: f64| {
            let (mut a, mut b) = (yn.clone(), yn.clone());
            a[var] += h;
            b[var] -= h;
            let (mut fa, mut fb) = ([0.0; 2], [0.0; 2]);
            problem.outer_bc(&a, &mut fa);
            problem.outer_bc(&b, &mut fb);
            [(fa[0] - fb[0]) / (2.0 * h), (fa[1] - fb[1]) / (2.0 * h)]
        };
        let (d1, d2) = (d(1e-4), d(1e-5));
        for i in 0..2 {
            assert!((d1[i] - d2[i]).abs() < 1e-4 * d1[i].abs().max(1.0), "var {var}, eq {i}: {} vs {}", d1[i], d2[i]);
        }
    }
}

/// Plan 7.3: the solver converges for the baseline masses with the tables
/// (other physics analytic). Below about 0.75 M☉ the tables with the fully
/// ionised ideal gas have no converged ZAMS model (molecular opacities at
/// T_eff ≈ 3200 K with μ = 0.6 gas): the low-mass end needs the full EOS
/// (`fully_convective_low_mass_stars`).
#[test]
fn tables_sequence_converges() {
    let phys = Physics::from_options(&options(AtmosphereChoice::Grey, EosChoice::Ideal, OpacityChoice::Tables));
    let models = sequence(&phys, &[1.0, 0.85, 3.0, 5.0, 10.0, 20.0]);
    let sun = &models[0].globals;
    eprintln!(
        "1 M☉ with tables: L = {:.4}, R = {:.4}, κ_c = {:.3}",
        sun.luminosity, sun.radius, models[0].profile.kappa[0]
    );
    // At the centre of the ZAMS Sun the tables give about half the analytic κ.
    assert!(models[0].profile.kappa[0] < 2.0);
}

/// Plan 6.2: ZAMS sequence with Fermi–Dirac electrons down to 0.1 M☉.
#[test]
#[ignore = "slow: Fermi–Dirac EOS"]
fn degenerate_eos_low_mass_sequence() {
    let phys = Physics::from_options(&options(AtmosphereChoice::Grey, EosChoice::Degenerate, OpacityChoice::Analytic));
    let models = sequence(&phys, &[1.0, 0.6, 0.4, 0.25, 0.15, 0.1]);
    let last = &models.last().unwrap().globals;
    eprintln!("0.1 M☉: η_c = {:.2}, R = {:.4}, L = {:.3e}", last.eta_c, last.radius, last.luminosity);
    assert!(last.eta_c > models[0].globals.eta_c);
}

/// Plan 6.5: with H₂ and ionisation, low-mass stars become fully convective
/// (reported: the mass below which no radiative core is left).
#[test]
#[ignore = "slow: full EOS"]
fn fully_convective_low_mass_stars() {
    let phys = Physics::from_options(&PhysicsOptions::all_upgrades());
    let models = sequence(&phys, &[1.0, 0.7, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2]);
    for m in &models {
        let g = &m.globals;
        let full = g.conv_core_q > 0.0 && g.conv_env_base_r.is_none() && m.zones.len() == 1;
        eprintln!(
            "M = {:.2}: zones {:?}, fully convective: {full}",
            g.mass,
            m.zones.iter().map(|z| (z.q_inner, z.q_outer)).collect::<Vec<_>>()
        );
    }
    let m02 = models.last().unwrap();
    assert!(m02.zones.len() == 1 && m02.zones[0].q_inner <= m02.profile.q[0] && m02.zones[0].q_outer >= 1.0 - 1e-9);
}

/// Plan 6.2: 1 M☉ past hydrogen exhaustion develops an isothermal, partly
/// degenerate helium core; "as far as possible" stops at a degenerate core.
#[test]
#[ignore = "slow: Fermi–Dirac EOS evolution"]
fn degenerate_helium_core() {
    let phys = Physics::from_options(&options(AtmosphereChoice::Grey, EosChoice::Degenerate, OpacityChoice::Analytic));
    let zams = solve_model(&StarParams { n_points: 400, ..StarParams::sun() }, &phys, None, &mut NoControl).unwrap();
    let opts = EvolutionOptions { stop: StopCondition::AsFarAsPossible, pulsation: false, ..Default::default() };
    let mut evo = Evolution::new(zams, opts);
    let reason = loop {
        match evo.step(&phys, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => {}
            EvolutionStep::Finished(r) => break r,
        }
    };
    let m = &evo.model;
    let p = &m.profile;
    let core = p.x.iter().take_while(|x| **x < 1e-3).count();
    eprintln!("stopped: {reason}; η_c = {:.2}, He core q = {:.3}", m.globals.eta_c, m.globals.he_core_q);
    assert!(core > 5);
    // Nearly isothermal core: T varies by less than 15% across it.
    let t = &p.t[..core];
    let (lo, hi) = t.iter().fold((f64::INFINITY, 0.0f64), |(a, b), v| (a.min(*v), b.max(*v)));
    assert!(hi / lo < 1.15, "core T from {lo:e} to {hi:e}");
    assert!(m.globals.eta_c > 0.0);
    let _ = (I_LNP, I_LNT);
}

/// Plan 7.8: the provisional solar calibration of the full set reproduces L☉
/// and R☉ at 4.57 Gyr (400 points, the time-step targets of the calibration).
#[test]
#[ignore = "slow: full physics evolution (about an hour on a desktop)"]
fn calibrated_sun() {
    let phys = Physics::from_options(&PhysicsOptions::all_upgrades());
    let p = StarParams { n_points: 400, ..StarParams::sun_calibrated() };
    let zams = solve_model(&p, &phys, None, &mut NoControl).unwrap();
    let opts = EvolutionOptions {
        stop: StopCondition::Age(4.57e9),
        pulsation: false,
        dlog_target: 0.05,
        dx_max: 0.05,
        ..Default::default()
    };
    let mut evo = Evolution::new(zams, opts);
    while let EvolutionStep::Advanced { .. } = evo.step(&phys, &mut NoControl).unwrap() {}
    let g = &evo.model.globals;
    eprintln!(
        "4.57 Gyr: L = {:.5}, R = {:.5}, T_eff = {:.0}, envelope base {:?}",
        g.luminosity, g.radius, g.teff, g.conv_env_base_r
    );
    // Rounding of the published parameters (4 digits) moves L by about 2e-3.
    assert!((g.luminosity - 1.0).abs() < 5e-3 && (g.radius - 1.0).abs() < 2e-3);
}

/// Plan 6.5: the Paczyński interpolation is run through the consistency scan
/// (it is not derived from a free energy, so ℛ is reported, not required to
/// vanish) and compared with the fully ionised Fermi–Dirac EOS along the
/// analytic ZAMS Sun.
#[test]
fn paczynski_comparison() {
    use stellar_core::physics::Eos;
    use stellar_core::physics::eos::{FermiDiracEos, Paczynski};
    use stellar_core::physics::eos_fd::stellar_scan;
    let c = StarParams::sun().composition().unwrap();
    let scan = stellar_scan(&Paczynski, &c, 10);
    eprintln!("Paczyński scan: {scan:?}");
    assert!(scan.unstable == 0, "c_V or α not positive somewhere");
    let m = solve_model(&StarParams { n_points: 200, ..StarParams::sun() }, &Physics::analytic(), None, &mut NoControl)
        .unwrap();
    let fd = FermiDiracEos::fully_ionised();
    let p = &m.profile;
    let (mut d_rho, mut d_nad, mut d_g1) = (0.0f64, 0.0f64, 0.0f64);
    for k in 0..p.len() {
        let (a, b) = (Paczynski.eval(p.p[k], p.t[k], &c), fd.eval(p.p[k], p.t[k], &c));
        d_rho = d_rho.max((a.rho / b.rho - 1.0).abs());
        d_nad = d_nad.max((a.nabla_ad / b.nabla_ad - 1.0).abs());
        d_g1 = d_g1.max((a.gamma1 / b.gamma1 - 1.0).abs());
    }
    eprintln!("along the ZAMS Sun, largest relative differences: ρ {d_rho:.2e}, ∇_ad {d_nad:.2e}, Γ₁ {d_g1:.2e}");
    assert!(d_rho < 0.05 && d_nad < 0.05 && d_g1 < 0.05);
}
