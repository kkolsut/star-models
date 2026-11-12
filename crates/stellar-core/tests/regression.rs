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
            ("L", g.luminosity, 0.41952023),
            ("R", g.radius, 1.11459091),
            ("T_eff", g.teff, 4400.0444),
            ("T_c", g.tc, 1.26075225e7),
            ("rho_c", g.rhoc, 63.502542),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.96764385),
            ("frac_pp", g.frac_pp, 0.99886455),
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
            ("L", g.luminosity, 73.854749),
            ("R", g.radius, 1.8225059),
            ("T_eff", g.teff, 12533.914),
            ("T_c", g.tc, 2.2822366e7),
            ("rho_c", g.rhoc, 34.119270),
            ("conv. core q", g.conv_core_q, 0.17468611),
            ("frac_cno", g.frac_cno, 0.86886396),
        ],
    );
    let g = zams(5.0).globals;
    check(
        "5 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 558.42314),
            ("R", g.radius, 2.5192201),
            ("T_eff", g.teff, 17678.051),
            ("T_c", g.tc, 2.5997958e7),
            ("rho_c", g.rhoc, 18.085324),
            ("conv. core q", g.conv_core_q, 0.23080913),
            ("frac_cno", g.frac_cno, 0.97134011),
        ],
    );
    let g = zams(10.0).globals;
    check(
        "10 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 6032.4536),
            ("R", g.radius, 3.8509253),
            ("T_eff", g.teff, 25921.942),
            ("T_c", g.tc, 3.0169210e7),
            ("rho_c", g.rhoc, 7.9638479),
            ("conv. core q", g.conv_core_q, 0.32941480),
        ],
    );
    let g = zams(20.0).globals;
    check(
        "20 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 44392.543),
            ("R", g.radius, 5.7511193),
            ("T_eff", g.teff, 34936.380),
            ("T_c", g.tc, 3.3979180e7),
            ("rho_c", g.rhoc, 4.0590849),
            ("conv. core q", g.conv_core_q, 0.47361224),
            ("frac_cno", g.frac_cno, 0.99892792),
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
            ("L", g.luminosity, 0.51702692),
            ("R", g.radius, 1.15453137),
            ("T_c", g.tc, 1.36448091e7),
            ("X_c", g.x_c, 0.50705607),
        ],
    );
}

#[test]
fn three_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(3.0), StopCondition::EndOfMainSequence).globals;
    check("3 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 273.30198)]);
}

#[test]
fn five_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(5.0), StopCondition::EndOfMainSequence).globals;
    check("5 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 72.714814)]);
}
