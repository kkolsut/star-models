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
            ("L", g.luminosity, 0.41653940),
            ("R", g.radius, 1.12184982),
            ("T_eff", g.teff, 4377.9746),
            ("T_c", g.tc, 1.25150713e7),
            ("rho_c", g.rhoc, 62.125456),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.96683445),
            ("frac_pp", g.frac_pp, 0.99950072),
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
            ("L", g.luminosity, 9.9058747e-3),
            ("R", g.radius, 0.70106649),
            ("T_eff", g.teff, 2174.8043),
            ("T_c", g.tc, 7.2082182e6),
            ("rho_c", g.rhoc, 45.859153),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.75515909),
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
            ("L", g.luminosity, 75.388663),
            ("R", g.radius, 1.7682680),
            ("T_eff", g.teff, 12790.250),
            ("T_c", g.tc, 2.3498645e7),
            ("rho_c", g.rhoc, 37.243291),
            ("conv. core q", g.conv_core_q, 0.17332463),
            ("frac_cno", g.frac_cno, 0.83122113),
        ],
    );
    let g = zams(5.0).globals;
    check(
        "5 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 566.23297),
            ("R", g.radius, 2.4363188),
            ("T_eff", g.teff, 18038.828),
            ("T_c", g.tc, 2.6829218e7),
            ("rho_c", g.rhoc, 19.851672),
            ("conv. core q", g.conv_core_q, 0.23173091),
            ("frac_cno", g.frac_cno, 0.96249364),
        ],
    );
    let g = zams(10.0).globals;
    check(
        "10 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 6073.6732),
            ("R", g.radius, 3.7194083),
            ("T_eff", g.teff, 26421.198),
            ("T_c", g.tc, 3.1164189e7),
            ("rho_c", g.rhoc, 8.7710502),
            ("conv. core q", g.conv_core_q, 0.33063148),
        ],
    );
    let g = zams(20.0).globals;
    check(
        "20 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 44531.538),
            ("R", g.radius, 5.5525146),
            ("T_eff", g.teff, 35583.499),
            ("T_c", g.tc, 3.5120723e7),
            ("rho_c", g.rhoc, 4.4807706),
            ("conv. core q", g.conv_core_q, 0.47453142),
            ("frac_cno", g.frac_cno, 0.99859599),
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
            ("L", g.luminosity, 0.51295845),
            ("R", g.radius, 1.16160019),
            ("T_c", g.tc, 1.35327828e7),
            ("X_c", g.x_c, 0.50958056),
        ],
    );
}

#[test]
fn three_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(3.0), StopCondition::EndOfMainSequence).globals;
    check("3 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 268.56162)]);
}

#[test]
fn five_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(5.0), StopCondition::EndOfMainSequence).globals;
    check("5 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 72.003796)]);
}
