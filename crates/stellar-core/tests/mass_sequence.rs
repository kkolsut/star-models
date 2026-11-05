//! Qualitative behaviour along the ZAMS: convective cores and CNO burning in
//! massive stars, deep convective envelopes in low-mass stars, and a monotonic
//! mass–luminosity relation.

use stellar_core::{NoControl, Physics, StarModel, StarParams, solve_model};

fn model(mass: f64, warm: Option<&StarModel>) -> StarModel {
    let p = StarParams { mass, ..StarParams::sun() };
    solve_model(&p, &Physics::analytic(), warm, &mut NoControl).unwrap_or_else(|e| panic!("M = {mass}: {e}"))
}

#[test]
fn three_solar_masses_has_convective_core_and_cno() {
    let m = model(3.0, None);
    assert!(m.globals.conv_core_q > 0.05, "core q = {}", m.globals.conv_core_q);
    assert!(m.globals.frac_cno > 0.5, "CNO fraction {}", m.globals.frac_cno);
}

#[test]
fn low_mass_star_has_deeper_envelope() {
    let sun = model(1.0, None);
    let half = model(0.5, Some(&sun));
    let r_sun = sun.globals.conv_env_base_r.unwrap();
    let r_half = half.globals.conv_env_base_r.unwrap();
    assert!(r_half < r_sun, "0.5 M☉ base {r_half} vs 1 M☉ base {r_sun}");
    assert_eq!(half.globals.conv_core_q, 0.0);
}

#[test]
fn zams_sequence_converges_and_is_monotonic() {
    let masses = [0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0];
    let mut prev: Option<StarModel> = None;
    let mut lums = Vec::new();
    for &m in &masses {
        let cur = model(m, prev.as_ref());
        assert!((cur.globals.virial_ratio - 1.0).abs() < 1e-3);
        lums.push(cur.globals.luminosity);
        prev = Some(cur);
    }
    assert!(lums.windows(2).all(|w| w[1] > w[0]), "L(M) must increase: {lums:?}");
}

#[test]
fn cold_start_works_across_masses() {
    for &m in &[0.3, 0.7, 2.5, 10.0, 40.0] {
        let cur = model(m, None);
        assert!(cur.globals.luminosity > 0.0);
    }
}
