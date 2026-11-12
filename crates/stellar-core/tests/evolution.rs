//! Stellar evolution: main-sequence trends, energy budget, end of the main
//! sequence, and radial pulsation periods.

use stellar_core::constants::YEAR;
use stellar_core::structure::evolve::{burning_rates, hydrogen_mass};
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Phase, Physics, StarModel, StarParams, StopCondition,
    solve_model,
};

fn run(mass: f64, stop: StopCondition) -> (Vec<StarModel>, String) {
    run_with(mass, EvolutionOptions { stop, ..Default::default() })
}

fn run_with(mass: f64, opts: EvolutionOptions) -> (Vec<StarModel>, String) {
    let physics = Physics::analytic();
    let params = StarParams { mass, ..StarParams::sun() };
    let zams = solve_model(&params, &physics, None, &mut NoControl).unwrap();
    let mut evo = Evolution::new(zams.clone(), opts);
    let mut models = vec![zams];
    loop {
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => models.push(evo.model.clone()),
            EvolutionStep::Finished(reason) => return (models, reason),
        }
    }
}

#[test]
fn sun_evolves_along_the_main_sequence() {
    let (models, reason) = run(1.0, StopCondition::Age(4.57e9));
    assert!(reason.contains("target age"), "{reason}");
    let (zams, now) = (&models[0], models.last().unwrap());
    assert!((now.globals.age_yr - 4.57e9).abs() < 1e3);
    assert_eq!(now.globals.phase, Phase::MainSequence);
    // Brighter, larger, hotter centre, less central hydrogen.
    assert!(now.globals.luminosity > 1.1 * zams.globals.luminosity);
    assert!(now.globals.radius > zams.globals.radius);
    assert!(now.globals.tc > zams.globals.tc);
    assert!(now.globals.x_c < 0.65 && now.globals.x_c > 0.2, "X_c = {}", now.globals.x_c);
    // The surface is untouched by burning.
    assert!((now.profile.x.last().unwrap() - 0.70).abs() < 1e-6);
    for m in &models {
        assert!((m.globals.virial_ratio - 1.0).abs() < 1e-3, "virial {}", m.globals.virial_ratio);
        assert!((m.globals.energy_ratio - 1.0).abs() < 1e-3, "energy {}", m.globals.energy_ratio);
    }
}

#[test]
fn hydrogen_burnt_matches_energy_radiated() {
    // Fuel accounting: hydrogen burnt × Q_H^eff against the energy radiated.
    let (models, _) = run(1.0, StopCondition::Age(4.57e9));
    let check = models.last().unwrap().globals.energy_check().unwrap();
    eprintln!("1 M☉ to 4.57 Gyr: fuel ratio {:.5}", check.fuel_ratio);
    assert!((check.fuel_ratio - 1.0).abs() < 0.03, "fuel ratio {}", check.fuel_ratio);
}

/// Plans 1.1 and 1.5: the burning step removes exactly the hydrogen implied by
/// the per-channel rates (ε_pp/Q_pp^eff + ε_CNO/Q_CNO, zone-averaged), and the
/// fuel accounting of a 3 M☉ star closes to 0.5% at the end of the main sequence.
#[test]
fn fuel_accounting_with_per_channel_energies() {
    let (models, _) = run(3.0, StopCondition::EndOfMainSequence);
    let mut worst: f64 = 0.0;
    for w in models.windows(2) {
        if w[1].globals.x_c < 0.05 {
            break; // clipping at exhaustion is allowed to remove less
        }
        let dt = (w[1].globals.age_yr - w[0].globals.age_yr) * YEAR;
        let (h, _) = burning_rates(&w[0]);
        let p = &w[0].profile;
        let expected: f64 = (0..p.len())
            .map(|k| {
                let lo = if k == 0 { 0.0 } else { 0.5 * (p.m[k] + p.m[k - 1]) };
                let hi = if k + 1 == p.len() { p.m[k] } else { 0.5 * (p.m[k] + p.m[k + 1]) };
                (hi - lo) * h[k] * dt
            })
            .sum();
        let burnt = hydrogen_mass(&w[0]) - hydrogen_mass(&w[1]);
        worst = worst.max((burnt / expected - 1.0).abs());
    }
    assert!(worst < 1e-9, "burning step against the rates: worst relative difference {worst:.2e}");
    let check = models.last().unwrap().globals.energy_check().unwrap();
    eprintln!("3 M☉ to the end of the main sequence: fuel ratio {:.5}", check.fuel_ratio);
    assert!((check.fuel_ratio - 1.0).abs() < 0.005, "fuel ratio {}", check.fuel_ratio);
}

/// Plan 1.10: the total-energy check closes to 1% on the main sequence, and on
/// halving the time step its numerical part shrinks while the composition term
/// (D17) converges to a limit.
#[test]
fn total_energy_check_on_the_main_sequence() {
    for (mass, stop) in [(1.0, StopCondition::Age(4.57e9)), (5.0, StopCondition::EndOfMainSequence)] {
        let base = EvolutionOptions { stop, pulsation: false, ..Default::default() };
        let half = EvolutionOptions { dx_max: 0.5 * base.dx_max, dlog_target: 0.5 * base.dlog_target, ..base };
        let a = run_with(mass, base).0.last().unwrap().globals.energy_check().unwrap();
        let b = run_with(mass, half).0.last().unwrap().globals.energy_check().unwrap();
        eprintln!(
            "M = {mass}: residual {:.3e} → {:.3e} (composition {:.3e} → {:.3e}, numerical {:.3e} → {:.3e}); fuel {:.5} → {:.5}",
            a.residual, b.residual, a.composition, b.composition, a.numerical, b.numerical, a.fuel_ratio, b.fuel_ratio
        );
        for c in [a, b] {
            assert!(c.residual.abs() < 0.01, "M = {mass}: residual {}", c.residual);
        }
        assert!(b.numerical.abs() < 0.75 * a.numerical.abs(), "M = {mass}: numerical part does not shrink");
        assert!(
            (b.composition - a.composition).abs() < 0.1 * a.composition.abs(),
            "M = {mass}: composition term does not converge"
        );
    }
}

#[test]
fn intermediate_mass_star_reaches_end_of_main_sequence() {
    let (models, reason) = run(3.0, StopCondition::EndOfMainSequence);
    assert!(reason.contains("hydrogen exhausted"), "{reason}");
    let (zams, end) = (&models[0], models.last().unwrap());
    let age_myr = end.globals.age_yr / 1e6;
    assert!(age_myr > 150.0 && age_myr < 800.0, "main-sequence lifetime {age_myr} Myr");
    // The convective core shrinks during the main sequence of a 3 M☉ star.
    let q_cc_mid = models[models.len() / 2].globals.conv_core_q;
    assert!(q_cc_mid < zams.globals.conv_core_q, "{q_cc_mid} vs {}", zams.globals.conv_core_q);
    assert!(end.globals.x_c < 1e-3);
}

#[test]
fn solar_fundamental_radial_mode_is_about_an_hour() {
    let physics = Physics::analytic();
    let sun = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).unwrap();
    let mode = stellar_core::pulsation::fundamental_mode(&sun, None).unwrap();
    assert_eq!(mode.nodes, 0);
    let hours = mode.period_s / 3600.0;
    eprintln!("solar ZAMS fundamental radial mode: {hours:.3} h, Ω² = {:.3}", mode.omega2);
    assert!(hours > 0.6 && hours < 1.6, "P = {hours} h");
    // ξ grows outward (the surface moves most).
    assert!(mode.xi[0].abs() < mode.xi.last().unwrap().abs());
}
