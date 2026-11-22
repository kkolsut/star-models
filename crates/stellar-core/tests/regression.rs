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
            ("L", g.luminosity, 0.41679912),
            ("R", g.radius, 1.1211539),
            ("T_eff", g.teff, 4380.0156),
            ("T_c", g.tc, 1.2523748e7),
            ("rho_c", g.rhoc, 62.251719),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.96690359),
            ("frac_pp", g.frac_pp, 0.99955309),
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
            ("L", g.luminosity, 0.0099180191),
            ("R", g.radius, 0.70049584),
            ("T_eff", g.teff, 2176.3565),
            ("T_c", g.tc, 7.2203005e6),
            ("rho_c", g.rhoc, 46.079184),
            ("env. base r/R", g.conv_env_base_r.unwrap_or(f64::NAN), 0.75555019),
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
            ("L", g.luminosity, 75.462632),
            ("R", g.radius, 1.7656348),
            ("T_eff", g.teff, 12802.923),
            ("T_c", g.tc, 2.3532408e7),
            ("rho_c", g.rhoc, 37.403350),
            ("conv. core q", g.conv_core_q, 0.17330211),
            ("frac_cno", g.frac_cno, 0.82973341),
        ],
    );
    let g = zams(5.0).globals;
    check(
        "5 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 566.42449),
            ("R", g.radius, 2.4343004),
            ("T_eff", g.teff, 18047.830),
            ("T_c", g.tc, 2.6850125e7),
            ("rho_c", g.rhoc, 19.897476),
            ("conv. core q", g.conv_core_q, 0.23176244),
            ("frac_cno", g.frac_cno, 0.96231021),
        ],
    );
    let g = zams(10.0).globals;
    check(
        "10 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 6074.0872),
            ("R", g.radius, 3.7180951),
            ("T_eff", g.teff, 26426.314),
            ("T_c", g.tc, 3.1174468e7),
            ("rho_c", g.rhoc, 8.7796603),
            ("conv. core q", g.conv_core_q, 0.33064479),
        ],
    );
    let g = zams(20.0).globals;
    check(
        "20 M☉ ZAMS",
        ZAMS_RTOL,
        &[
            ("L", g.luminosity, 44532.184),
            ("R", g.radius, 5.5515949),
            ("T_eff", g.teff, 35586.575),
            ("T_c", g.tc, 3.5126196e7),
            ("rho_c", g.rhoc, 4.4828596),
            ("conv. core q", g.conv_core_q, 0.47453577),
            ("frac_cno", g.frac_cno, 0.99859465),
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
            ("L", g.luminosity, 0.51340759),
            ("R", g.radius, 1.1608888),
            ("T_c", g.tc, 1.3542437e7),
            ("X_c", g.x_c, 0.50957201),
        ],
    );
}

#[test]
fn three_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(3.0), StopCondition::EndOfMainSequence).globals;
    check("3 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 268.37341)]);
}

#[test]
fn five_solar_masses_main_sequence_lifetime() {
    let g = evolve(zams(5.0), StopCondition::EndOfMainSequence).globals;
    check("5 M☉ end of main sequence", EVO_RTOL, &[("age [Myr]", g.age_yr / 1e6, 72.002996)]);
}
