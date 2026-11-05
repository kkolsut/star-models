//! End-to-end checks of the solar ZAMS model built with analytic physics.
//!
//! The analytic (Kramers + H⁻) opacities are only qualitatively right, so the
//! bands on L and R are wide. Structural features (radiative core, convective
//! envelope, pp dominance) and the global consistency checks are strict.

use stellar_core::{NoControl, Physics, StarParams, solve_model};

#[test]
fn solar_model_structure() {
    let physics = Physics::analytic();
    let m = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).expect("Sun must converge");
    let g = &m.globals;

    // Global consistency of the discretised solution.
    assert!((g.virial_ratio - 1.0).abs() < 1e-3, "virial ratio {}", g.virial_ratio);
    assert!((g.energy_ratio - 1.0).abs() < 1e-3, "energy ratio {}", g.energy_ratio);

    // Loose physical bands (analytic opacity).
    assert!(g.luminosity > 0.3 && g.luminosity < 1.0, "L = {}", g.luminosity);
    assert!(g.radius > 0.8 && g.radius < 1.2, "R = {}", g.radius);
    assert!(g.tc > 1.1e7 && g.tc < 1.5e7, "Tc = {}", g.tc);
    assert!(g.rhoc > 40.0 && g.rhoc < 150.0, "ρc = {}", g.rhoc);

    // Structure: radiative core, convective envelope, pp chain dominates.
    assert_eq!(g.conv_core_q, 0.0, "the Sun has no convective core");
    let base = g.conv_env_base_r.expect("the Sun has a convective envelope");
    assert!(base > 0.5 && base < 1.0, "envelope base r/R = {base}");
    assert!(g.frac_pp > 0.9, "pp fraction {}", g.frac_pp);

    // Profiles: P, T, ρ decrease outward; L increases to the surface value.
    let p = &m.profile;
    assert!(p.p.windows(2).all(|w| w[1] < w[0]));
    assert!(p.t.windows(2).all(|w| w[1] <= w[0] * (1.0 + 1e-9)));
    assert!(p.l.windows(2).all(|w| w[1] >= w[0] * (1.0 - 1e-9)));
}

#[test]
fn mesh_refinement_changes_little() {
    let physics = Physics::analytic();
    let coarse = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).unwrap();
    let fine_params = StarParams { n_points: 1600, ..StarParams::sun() };
    let fine = solve_model(&fine_params, &physics, None, &mut NoControl).unwrap();
    let rel = |a: f64, b: f64| (a / b - 1.0).abs();
    assert!(rel(coarse.globals.luminosity, fine.globals.luminosity) < 5e-3);
    assert!(rel(coarse.globals.radius, fine.globals.radius) < 5e-3);
    assert!(rel(coarse.globals.tc, fine.globals.tc) < 5e-3);
}

#[test]
fn warm_start_matches_cold_start() {
    let physics = Physics::analytic();
    let sun = solve_model(&StarParams::sun(), &physics, None, &mut NoControl).unwrap();
    let target = StarParams { mass: 1.3, x: 0.72, z: 0.018, ..StarParams::sun() };
    let warm = solve_model(&target, &physics, Some(&sun), &mut NoControl).unwrap();
    let cold = solve_model(&target, &physics, None, &mut NoControl).unwrap();
    assert_eq!(warm.params, target);
    let rel = |a: f64, b: f64| (a / b - 1.0).abs();
    assert!(rel(warm.globals.luminosity, cold.globals.luminosity) < 1e-6);
    assert!(rel(warm.globals.radius, cold.globals.radius) < 1e-6);
}

#[test]
fn invalid_parameters_are_rejected() {
    let physics = Physics::analytic();
    let bad = StarParams { x: 0.9, z: 0.2, ..StarParams::sun() };
    assert!(solve_model(&bad, &physics, None, &mut NoControl).is_err());
}
