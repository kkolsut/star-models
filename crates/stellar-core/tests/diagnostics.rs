//! Output-only diagnostics of stage 1: gravitational redshift (plan 1.8) and the
//! screening diagnostic (plan 1.11).

use stellar_core::Composition;
use stellar_core::constants::{M_SUN, R_SUN};
use stellar_core::physics::Nuclear;
use stellar_core::physics::nuclear::AnalyticNuclear;
use stellar_core::structure::model::{SCREENING_WEAK_LIMIT, grav_redshift_m_s, screening_warning};
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarParams, StopCondition, solve_model,
};

#[test]
fn gravitational_redshift() {
    let v = grav_redshift_m_s(M_SUN, R_SUN);
    assert!((v - 636.0).abs() < 0.5, "v_z = {v} m/s");
    let sun = solve_model(&StarParams::sun(), &Physics::analytic(), None, &mut NoControl).unwrap();
    let g = &sun.globals;
    assert!((g.grav_redshift_m_s - 636.0 / g.radius).abs() < 0.5, "{} m/s at {} R☉", g.grav_redshift_m_s, g.radius);
}

#[test]
fn screening_exponent_at_the_solar_centre() {
    let physics = Physics::analytic();
    let zams = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).unwrap();
    let s = &zams.globals.screening;
    let p = &zams.profile;
    let c = p.composition_at(0, 0.7);
    let centre = AnalyticNuclear::default().screening(p.rho[0], p.t[0], &c).unwrap();
    eprintln!("ZAMS Sun, centre: exponents pp chain (largest) {:.4}, p + ¹⁴N {:.4}", centre.pp, centre.cno);
    eprintln!("ZAMS Sun: largest screening exponent {:.4} ({}) at m/M = {:.4}", s.max_exponent, s.reaction, s.q);
    assert!(s.max_exponent > 0.0 && s.max_exponent < SCREENING_WEAK_LIMIT);
    assert!(s.reaction.contains("pp"), "{}", s.reaction);
    // Plan 9.7: the all-regime screening has no cap.
    assert_eq!(s.cap, None);
    assert!(s.warning().is_none());
    // Along a run the largest value is kept.
    let mut evo = Evolution::new(
        zams,
        EvolutionOptions { stop: StopCondition::Age(4.57e9), pulsation: false, ..Default::default() },
    );
    while let EvolutionStep::Advanced { .. } = evo.step(&physics, &mut NoControl).unwrap() {}
    let g = &evo.model.globals;
    let run = g.run.unwrap();
    eprintln!(
        "Sun at 4.57 Gyr: largest screening exponent {:.4} ({}); largest of the run {:.4}",
        g.screening.max_exponent, g.screening.reaction, run.max_screening
    );
    assert!(run.max_screening >= g.screening.max_exponent);
}

#[test]
fn dense_states_trigger_each_warning() {
    let c = Composition::new(0.0, 0.02, 0.7).unwrap();
    let cap = AnalyticNuclear::default().screening_cap();
    assert_eq!(cap, None);
    // Helium at 10⁸ K: weak at 10³ g/cm³, beyond 0.3 at 10⁵ and 10⁷ (no cap).
    let e = |rho: f64| AnalyticNuclear::default().screening(rho, 1e8, &c).unwrap().triple_alpha;
    assert!(screening_warning(e(1e3), cap).is_none(), "{}", e(1e3));
    for rho in [1e5, 1e7] {
        let w = screening_warning(e(rho), cap).unwrap();
        assert!(w.contains("weak regime"), "{} → {w}", e(rho));
    }
    // A capped implementation still gets the cap warning.
    assert!(screening_warning(2.5, Some(2.0)).unwrap().contains("cap"));
    assert!(e(1e7) > e(1e5) && e(1e5) > e(1e3));
}
