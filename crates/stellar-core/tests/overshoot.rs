//! Step overshoot of convective cores (plan 2.1): off by default, mixing only,
//! longer main sequences for stars with convective cores, no effect on stars
//! without one, and the fuel accounting still closes.

use stellar_core::constants::{L_SUN, YEAR};
use stellar_core::structure::evolve::{Q_H, core_overshoot};
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

fn zams(mass: f64, alpha_ov: f64) -> StarModel {
    let p = StarParams { mass, alpha_ov, ..StarParams::sun() };
    solve_model(&p, &Physics::analytic(), None, &mut NoControl).unwrap()
}

/// ∫ X dm over the star [g].
fn hydrogen_mass(m: &StarModel) -> f64 {
    let p = &m.profile;
    (1..p.len()).map(|k| 0.5 * (p.x[k] + p.x[k - 1]) * (p.m[k] - p.m[k - 1])).sum()
}

struct Run {
    end: StarModel,
    /// Hydrogen burnt × Q_H / energy radiated.
    fuel_ratio: f64,
    /// Largest |energy_ratio − 1| of any model.
    worst_energy: f64,
    /// Largest mixed-core mass of any model [m/M], and whether it always contained the convective core.
    max_mixed_q: f64,
    mixed_contains_convective: bool,
}

fn run(mass: f64, alpha_ov: f64, stop: StopCondition) -> Run {
    let physics = Physics::analytic();
    let z = zams(mass, alpha_ov);
    let h0 = hydrogen_mass(&z);
    let mut evo = Evolution::new(z, EvolutionOptions { stop, pulsation: false, ..Default::default() });
    let (mut radiated, mut worst_energy, mut max_mixed_q, mut contains) = (0.0, 0.0_f64, 0.0_f64, true);
    loop {
        let (age0, l0) = (evo.model.globals.age_yr, evo.model.globals.luminosity);
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => {
                let g = &evo.model.globals;
                radiated += 0.5 * (l0 + g.luminosity) * L_SUN * (g.age_yr - age0) * YEAR;
                worst_energy = worst_energy.max((g.energy_ratio - 1.0).abs());
                max_mixed_q = max_mixed_q.max(g.mixed_core_q);
                contains &= g.mixed_core_q >= g.conv_core_q;
            }
            EvolutionStep::Finished(_) => break,
        }
    }
    let fuel_ratio = (h0 - hydrogen_mass(&evo.model)) * Q_H / radiated;
    Run { end: evo.model, fuel_ratio, worst_energy, max_mixed_q, mixed_contains_convective: contains }
}

#[test]
fn overshoot_is_off_by_default_and_leaves_the_zams_unchanged() {
    assert_eq!(StarParams::sun().alpha_ov, 0.0);
    let (plain, ov) = (zams(3.0, 0.0), zams(3.0, 0.2));
    // Mixing only: the static structure does not depend on α_ov.
    assert_eq!(plain.globals.luminosity, ov.globals.luminosity);
    assert_eq!(plain.globals.radius, ov.globals.radius);
    assert_eq!(plain.globals.conv_core_q, ov.globals.conv_core_q);
    // Without overshoot the mixed core is the convective core.
    assert!(core_overshoot(&plain.profile, &plain.zones, 0.0).is_none());
    assert_eq!(plain.globals.mixed_core_q, plain.globals.conv_core_q);
    // With it, the mixed core reaches r_cc + α_ov min(H_P, r_cc).
    let o = core_overshoot(&ov.profile, &ov.zones, 0.2).unwrap();
    assert!(o.r_ov > o.r_cc);
    assert!((ov.globals.mixed_core_r - o.r_ov / ov.profile.r.last().unwrap()).abs() < 1e-12);
    assert!(
        ov.globals.mixed_core_q > ov.globals.conv_core_q,
        "{} vs {}",
        ov.globals.mixed_core_q,
        ov.globals.conv_core_q
    );
    // Out-of-range values are rejected.
    assert!(StarParams { alpha_ov: 0.6, ..StarParams::sun() }.validate().is_err());
    assert!(StarParams { alpha_ov: -0.1, ..StarParams::sun() }.validate().is_err());
}

#[test]
fn main_sequence_lifetime_grows_with_overshoot() {
    let cases: Vec<(f64, f64)> = [3.0, 5.0].iter().flat_map(|&m| [0.0, 0.1, 0.2].map(|a| (m, a))).collect();
    let runs: Vec<Run> = std::thread::scope(|s| {
        let handles: Vec<_> =
            cases.iter().map(|&(m, a)| s.spawn(move || run(m, a, StopCondition::EndOfMainSequence))).collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for (i, (&(m, a), r)) in cases.iter().zip(&runs).enumerate() {
        let age = r.end.globals.age_yr / 1e6;
        eprintln!(
            "M = {m} M☉, α_ov = {a}: main sequence {age:.2} Myr, max mixed core {:.4}, fuel ratio {:.5}, worst energy {:.1e}",
            r.max_mixed_q, r.fuel_ratio, r.worst_energy
        );
        assert!(r.end.globals.x_c < 1e-3, "M = {m}, α_ov = {a} did not reach the end of the main sequence");
        // Fuel accounting and the discretisation energy check close with overshoot too.
        assert!((r.fuel_ratio - 1.0).abs() < 0.03, "M = {m}, α_ov = {a}: fuel ratio {}", r.fuel_ratio);
        assert!(r.worst_energy < 1e-3, "M = {m}, α_ov = {a}: energy ratio off by {}", r.worst_energy);
        assert!(r.mixed_contains_convective);
        if i % 3 > 0 {
            let prev = &runs[i - 1];
            assert!(
                r.end.globals.age_yr > 1.05 * prev.end.globals.age_yr,
                "M = {m}: lifetime does not grow from α_ov = {} to {a}",
                cases[i - 1].1
            );
            assert!(r.max_mixed_q > prev.max_mixed_q);
        }
    }
}

#[test]
fn star_without_convective_core_is_unaffected() {
    let stop = StopCondition::Age(5e9);
    let (plain, ov) = std::thread::scope(|s| {
        let a = s.spawn(|| run(1.1, 0.0, stop));
        let b = s.spawn(|| run(1.1, 0.2, stop));
        (a.join().unwrap(), b.join().unwrap())
    });
    assert_eq!(plain.max_mixed_q, 0.0, "1.1 M☉ is expected to keep a radiative core");
    assert_eq!(ov.max_mixed_q, 0.0);
    let (a, b) = (&plain.end.globals, &ov.end.globals);
    assert_eq!(a.age_yr, b.age_yr);
    assert_eq!(a.x_c, b.x_c);
    assert_eq!(a.luminosity, b.luminosity);
    assert_eq!(a.radius, b.radius);
}
