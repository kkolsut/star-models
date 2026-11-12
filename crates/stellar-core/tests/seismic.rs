//! Seismic diagnostics (plan 2.2): buoyancy frequency against the convective
//! flag and against its definition, asymptotic Δν against the scaling relation,
//! and the candidate g-mode cavities.

use stellar_core::seismic::{self, CavityLocation};
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

fn zams(mass: f64) -> StarModel {
    let p = StarParams { mass, ..StarParams::sun() };
    solve_model(&p, &Physics::analytic(), None, &mut NoControl).unwrap()
}

/// Tolerance on the dimensionless buoyancy A = N² P / (g² ρ) around zero in
/// nearly adiabatic convective layers.
const A_TOL: f64 = 1e-6;

#[test]
fn zams_buoyancy_matches_the_convective_flag() {
    for mass in [0.5, 1.0, 1.5, 3.0, 5.0] {
        let m = zams(mass);
        let p = &m.profile;
        let g = p.gravity();
        assert_eq!(p.n2.len(), p.len());
        assert_eq!(p.cs.len(), p.len());
        for k in 0..p.len() {
            let a = seismic::dimensionless_buoyancy(p, &g, k);
            if p.convective[k] {
                assert!(a < A_TOL, "M = {mass}: convective point {k} has A = {a:.3e}");
            } else if k > 0 && k + 1 < p.len() {
                assert!(a > 0.0, "M = {mass}: radiative point {k} has A = {a:.3e}");
            }
            // Homogeneous: Schwarzschild and Ledoux agree.
            if k + 1 < p.len() && a.abs() > A_TOL {
                assert_eq!(p.schwarzschild_unstable(k), p.ledoux_unstable(k), "M = {mass}, k = {k}");
            }
        }
    }
}

#[test]
fn evolved_buoyancy_agrees_with_its_definition() {
    // 3 M☉ in the middle of the main sequence: a μ gradient above the shrinking core.
    let physics = Physics::analytic();
    let opts = EvolutionOptions { stop: StopCondition::Age(2.0e8), pulsation: false, ..Default::default() };
    let mut evo = Evolution::new(zams(3.0), opts);
    while let EvolutionStep::Advanced { .. } = evo.step(&physics, &mut NoControl).unwrap() {}
    let m = &evo.model;
    let p = &m.profile;
    let n = p.len();
    let g = p.gravity();
    let def = seismic::n2_from_definition(p, &g);
    let nmu = seismic::nabla_mu(&p.mu, &p.p);
    assert!(nmu.iter().any(|&v| v > 1e-2), "expected a composition gradient");

    // Compare in stably stratified layers, away from composition jumps (no ∇_μ
    // within two points), where finite differences are meaningful. In the
    // nearly adiabatic radiative layers just below the photosphere (A ~ 10⁻³)
    // the definition cancels badly, so the relative test needs A > 10⁻²; the
    // absolute difference in A is tested everywhere.
    let a_of = |n2: f64, k: usize| n2 * p.p[k] / (g[k] * g[k] * p.rho[k]);
    let pts: Vec<usize> = (2..n - 2)
        .filter(|&k| a_of(p.n2[k], k) > 1e-3)
        .filter(|&k| (k - 2..=k + 2).all(|i| nmu[i].abs() < 1e-3))
        .collect();
    assert!(pts.len() > 200, "{} points", pts.len());
    let worst_abs = pts.iter().map(|&k| (a_of(def[k], k) - a_of(p.n2[k], k)).abs()).fold(0.0, f64::max);
    let mut rel: Vec<f64> =
        pts.iter().filter(|&&k| a_of(p.n2[k], k) > 1e-2).map(|&k| ((def[k] - p.n2[k]) / p.n2[k]).abs()).collect();
    rel.sort_by(f64::total_cmp);
    let (median, worst) = (rel[rel.len() / 2], *rel.last().unwrap());
    eprintln!(
        "N²: production vs definition over {} points: median {median:.2e}, worst {worst:.2e} (A > 10⁻²); \
         worst |ΔA| {worst_abs:.2e}",
        pts.len()
    );
    assert!(median < 1e-3, "median relative difference {median}");
    assert!(worst < 0.05, "worst relative difference {worst}");
    assert!(worst_abs < 2e-3, "worst difference in A {worst_abs}");

    // Schwarzschild and Ledoux instability are evaluated separately. With
    // ∇_μ ≥ 0 every Ledoux-unstable layer is Schwarzschild-unstable, not the reverse.
    let schwarzschild = (0..n).filter(|&k| p.schwarzschild_unstable(k)).count();
    let ledoux = (0..n).filter(|&k| p.ledoux_unstable(k)).count();
    eprintln!("Schwarzschild-unstable points: {schwarzschild}, Ledoux-unstable points: {ledoux}");
    for k in 0..n - 1 {
        if p.ledoux_unstable(k) && seismic::dimensionless_buoyancy(p, &g, k) < -A_TOL {
            assert!(p.schwarzschild_unstable(k), "k = {k}: Ledoux- but not Schwarzschild-unstable");
        }
    }
    assert!(ledoux <= schwarzschild);
}

#[test]
fn solar_asymptotic_delta_nu_is_close_to_the_scaling_value() {
    let g = zams(1.0).globals;
    let ratio = g.delta_nu_asym_uhz / g.delta_nu_uhz;
    eprintln!("ZAMS Sun: Δν asymptotic {:.2} μHz, scaling {:.2} μHz", g.delta_nu_asym_uhz, g.delta_nu_uhz);
    assert!((ratio - 1.0).abs() < 0.10, "ratio {ratio}");
}

#[test]
fn g_mode_cavities_follow_the_radiative_regions() {
    // The Sun: a radiative core up to the base of the convective envelope.
    let sun = zams(1.0);
    let c = seismic::main_cavity(&sun.globals.g_cavities).expect("radiative core");
    assert_eq!(c.location, CavityLocation::Core);
    let env = sun.globals.conv_env_base_r.unwrap();
    assert!((c.r_outer - env).abs() < 0.01, "cavity ends at {} vs envelope base {env}", c.r_outer);
    assert!(c.delta_pi1_s.is_finite() && c.delta_pi1_s > 0.0);

    // 1.5 M☉: convective core and a thin convective envelope; the main cavity lies between them.
    let m = zams(1.5);
    let c = seismic::main_cavity(&m.globals.g_cavities).unwrap();
    assert_eq!(c.location, CavityLocation::Intermediate);
    assert!(c.r_inner >= m.globals.conv_core_r - 0.01 && c.r_outer <= m.globals.conv_env_base_r.unwrap() + 0.01);

    // Every listed cavity lies in a radiative region and has a finite spacing.
    for model in [&sun, &m] {
        let r_surf = *model.profile.r.last().unwrap();
        for cav in &model.globals.g_cavities {
            assert!(cav.delta_pi1_s.is_finite() && cav.delta_pi1_s > 0.0);
            let k = model.profile.r.partition_point(|&r| r / r_surf < 0.5 * (cav.r_inner + cav.r_outer));
            assert!(!model.profile.convective[k], "cavity {cav:?} centred in a convective layer");
        }
    }
}
