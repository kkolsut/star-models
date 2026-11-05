//! Regression tests: pin the current numerical results of the analytic physics.
//!
//! These are **not** validation against real stars (the analytic microphysics
//! is known to be off, see docs/about/15-limitations.md). They catch any
//! unintended change of the numbers: a failing test means a code change moved
//! the results, and the reference values must be updated deliberately, with a
//! reason, if the change is intended.
//!
//! The qualitative physics checks (radiative core, pp dominance, loose bands
//! around the Sun, ...) live in the other test files.

use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

/// Relative tolerance for static (ZAMS) models, converged to 1e-8 in the unknowns.
const ZAMS_RTOL: f64 = 1e-4;
/// Relative tolerance for evolved models. Tight enough to catch a change of the
/// time-step controls (halving the step moves X_c at 4.57 Gyr by ~4e-4).
const EVO_RTOL: f64 = 1e-4;

/// Compare (name, actual, expected) triples and report every mismatch at once.
fn check(label: &str, rtol: f64, values: &[(&str, f64, f64)]) {
    let bad: Vec<String> = values
        .iter()
        .filter(|(_, a, e)| {
            let rel = ((a - e) / e).abs();
            rel.is_nan() || rel > rtol
        })
        .map(|(n, a, e)| format!("  {n}: got {a:.8e}, expected {e:.8e} (rel. diff {:.2e})", (a - e) / e))
        .collect();
    assert!(bad.is_empty(), "{label}: results changed (rtol {rtol:.0e}):\n{}", bad.join("\n"));
}

fn zams(mass: f64) -> StarModel {
    let p = StarParams { mass, ..StarParams::sun() };
    solve_model(&p, &Physics::analytic(), None, &mut NoControl).unwrap_or_else(|e| panic!("M = {mass}: {e}"))
}

fn evolve(zams: StarModel, stop: StopCondition) -> StarModel {
    let physics = Physics::analytic();
    let opts = EvolutionOptions { stop, pulsation: false, ..Default::default() };
    let mut evo = Evolution::new(zams, opts);
    loop {
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => {}
            EvolutionStep::Finished(_) => return evo.model,
        }
    }
}

#[test]
fn solar_zams() {
    let m = zams(1.0);
    let g = &m.globals;
    check(
        "1 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 0.42066039),
            ("R", g.radius, 1.11332042),
            ("T_eff", g.teff, 4405.5425),
            ("T_c", g.tc, 1.26184629e7),
            ("rho_c", g.rhoc, 63.923655),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.96782054),
            ("frac_pp", g.frac_pp, 0.99883063),
        ],
    );
    assert_eq!(g.conv_core_q, 0.0);
}

#[test]
fn low_mass_zams() {
    let g = zams(0.5).globals;
    check(
        "0.5 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 9.9597925e-3),
            ("R", g.radius, 0.69900385),
            ("T_eff", g.teff, 2180.9684),
            ("T_c", g.tc, 7.2513022e6),
            ("rho_c", g.rhoc, 46.739435),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.75669948),
        ],
    );
}

#[test]
fn intermediate_and_massive_zams() {
    let g = zams(3.0).globals;
    check(
        "3 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 73.915913),
            ("R", g.radius, 1.8190833),
            ("T_eff", g.teff, 12548.296),
            ("T_c", g.tc, 2.2862448e7),
            ("rho_c", g.rhoc, 34.287181),
            ("conv. core q", g.conv_core_q, 0.17551975),
            ("frac_cno", g.frac_cno, 0.89994527),
        ],
    );
    let g = zams(10.0).globals;
    check(
        "10 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 6032.5176),
            ("R", g.radius, 3.8504842),
            ("T_eff", g.teff, 25923.495),
            ("T_c", g.tc, 3.0172342e7),
            ("rho_c", g.rhoc, 7.9662620),
            ("conv. core q", g.conv_core_q, 0.32945805),
        ],
    );
}

#[test]
fn sun_at_present_age() {
    let g = evolve(zams(1.0), StopCondition::Age(4.57e9)).globals;
    check(
        "1 M☉ at 4.57 Gyr",
        EVO_RTOL,
        &[
            ("L", g.luminosity, 0.52122620),
            ("R", g.radius, 1.15126839),
            ("T_c", g.tc, 1.36616900e7),
            ("X_c", g.x_c, 0.51311136),
        ],
    );
}

#[test]
fn three_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(3.0), StopCondition::EndOfMainSequence).globals;
    check("3 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 285.88946)]);
}
