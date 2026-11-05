//! Stellar evolution: main-sequence trends, energy budget, end of the main
//! sequence, and radial pulsation periods.

use stellar_core::constants::{L_SUN, YEAR};
use stellar_core::structure::evolve::Q_H;
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Phase, Physics, StarModel, StarParams, StopCondition,
    solve_model,
};

fn run(mass: f64, stop: StopCondition) -> (Vec<StarModel>, String) {
    let physics = Physics::analytic();
    let params = StarParams { mass, ..StarParams::sun() };
    let zams = solve_model(&params, &physics, None, &mut NoControl).unwrap();
    let mut evo = Evolution::new(zams.clone(), EvolutionOptions { stop, ..Default::default() });
    let mut models = vec![zams];
    loop {
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => models.push(evo.model.clone()),
            EvolutionStep::Finished(reason) => return (models, reason),
        }
    }
}

/// ∫ X dm over the star [g].
fn hydrogen_mass(m: &StarModel) -> f64 {
    let p = &m.profile;
    (1..p.len()).map(|k| 0.5 * (p.x[k] + p.x[k - 1]) * (p.m[k] - p.m[k - 1])).sum()
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
    let (models, _) = run(1.0, StopCondition::Age(4.57e9));
    let burnt = hydrogen_mass(&models[0]) - hydrogen_mass(models.last().unwrap());
    let radiated: f64 = models
        .windows(2)
        .map(|w| {
            let dt = (w[1].globals.age_yr - w[0].globals.age_yr) * YEAR;
            0.5 * (w[0].globals.luminosity + w[1].globals.luminosity) * L_SUN * dt
        })
        .sum();
    let ratio = burnt * Q_H / radiated;
    eprintln!("nuclear energy released / energy radiated = {ratio:.5}");
    assert!((ratio - 1.0).abs() < 0.03, "nuclear energy / radiated energy = {ratio}");
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
