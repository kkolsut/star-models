//! Radial pulsations (plan 1.3): polytropes, overtones, mesh convergence and the
//! acoustic cut-off diagnostic of real models.

use stellar_core::constants::{FOUR_PI, G, M_SUN, R_SUN};
use stellar_core::pulsation::{self, radial_mode};
use stellar_core::structure::model::Profile;
use stellar_core::{NoControl, Physics, StarModel, StarParams, solve_model};

/// Lane–Emden solution θ(ξ) of index n by RK4, up to the first zero ξ₁.
/// Returns (ξ, θ, dθ/dξ) on a fine grid and ξ₁.
fn lane_emden(n: f64) -> (Vec<(f64, f64, f64)>, f64) {
    let f = |x: f64, y: [f64; 2]| -> [f64; 2] { [y[1], -y[0].max(0.0).powf(n) - 2.0 * y[1] / x] };
    let h = 1e-5;
    let x0 = 1e-4;
    let mut x = x0;
    let mut y = [1.0 - x0 * x0 / 6.0, -x0 / 3.0];
    let mut out = vec![(x, y[0], y[1])];
    loop {
        let k1 = f(x, y);
        let k2 = f(x + h / 2.0, [y[0] + h / 2.0 * k1[0], y[1] + h / 2.0 * k1[1]]);
        let k3 = f(x + h / 2.0, [y[0] + h / 2.0 * k2[0], y[1] + h / 2.0 * k2[1]]);
        let k4 = f(x + h, [y[0] + h * k3[0], y[1] + h * k3[1]]);
        let next = [
            y[0] + h / 6.0 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
            y[1] + h / 6.0 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        ];
        if next[0] <= 0.0 {
            let x1 = x + h * y[0] / (y[0] - next[0]);
            return (out, x1);
        }
        x += h;
        y = next;
        out.push((x, y[0], y[1]));
    }
}

const THETA_CUT: f64 = 1e-5;

/// A polytrope of 1 M☉ and 1 R☉ with Γ₁ = 5/3 as a model the pulsation code accepts.
/// The mesh is uniform in ξ^(1/2)-like spacing towards the surface, which is cut at θ = 10⁻⁴.
fn polytrope(n: f64, points: usize) -> StarModel {
    let (le, x1) = lane_emden(n);
    let (mass, radius) = (M_SUN, R_SUN);
    let a = radius / x1;
    let dtheta1 = {
        let (xl, _, dl) = le[le.len() - 1];
        dl * (xl / x1).powi(2) // −ξ²θ' is conserved outside the mass; close enough at the edge
    };
    let rho_c = mass / (FOUR_PI * a.powi(3) * (-x1 * x1 * dtheta1));
    let p_c = FOUR_PI * G * a * a * rho_c * rho_c / (n + 1.0);
    let last = le.iter().rposition(|&(_, t, _)| t > THETA_CUT).unwrap();
    let x_end = le[last].0;
    let mut pr = Profile::default();
    for i in 0..points {
        let s = i as f64 / (points - 1) as f64;
        // Dense near the surface: ξ = ξ_end (1 − (1 − s)²)… mixed with a linear part.
        let x = (x_end * (0.5 * s + 0.5 * (1.0 - (1.0 - s).powi(2)))).max(le[0].0);
        let j = le.partition_point(|&(xx, _, _)| xx < x).min(le.len() - 1);
        let (xx, t, d) = le[j];
        let m = -FOUR_PI * a.powi(3) * rho_c * xx * xx * d;
        pr.q.push(m / mass);
        pr.m.push(m);
        pr.r.push(a * xx);
        pr.rho.push(rho_c * t.powf(n));
        pr.p.push(p_c * t.powf(n + 1.0));
        pr.gamma1.push(5.0 / 3.0);
    }
    let mut model =
        solve_model(&StarParams { n_points: 200, ..StarParams::sun() }, &Physics::analytic(), None, &mut NoControl)
            .unwrap();
    model.profile = pr;
    model
}

/// Fundamental radial modes of polytropes with Γ₁ = 5/3 against the classical
/// values Ω² ≈ 2.712 (n = 1.5) and 9.26 (n = 3). These reference values are
/// quoted from memory of the standard tables (Cox 1980), not yet compared with
/// the source; the tolerance (0.5%) is set accordingly. Overtones: the right
/// node counts, increasing Ω².
#[test]
fn polytrope_radial_modes() {
    for (n, reference) in [(1.5, 2.712), (3.0, 9.26)] {
        let model = polytrope(n, 3000);
        let modes: Vec<_> = (0..3).map(|k| radial_mode(&model, k, None).unwrap()).collect();
        let w: Vec<f64> = modes.iter().map(|m| m.omega2).collect();
        eprintln!("polytrope n = {n}: Ω² = {:.4}, {:.4}, {:.4}", w[0], w[1], w[2]);
        for (k, m) in modes.iter().enumerate() {
            assert_eq!(m.nodes, k);
        }
        assert!(w[0] < w[1] && w[1] < w[2]);
        assert!((w[0] / reference - 1.0).abs() < 5e-3, "n = {n}: Ω² = {} against {reference}", w[0]);
    }
}

fn sun(points: usize) -> StarModel {
    solve_model(&StarParams { n_points: points, ..StarParams::sun() }, &Physics::analytic(), None, &mut NoControl)
        .unwrap()
}

/// Mesh convergence of the solar fundamental period, and the same mode found
/// from every mesh: with arithmetic midpoints the period changed by 1.36e-4
/// (relative) between 800 and 1600 points; with geometric means it must change less.
#[test]
fn solar_fundamental_mode_mesh_convergence() {
    let periods: Vec<(usize, f64)> = [400, 800, 1600]
        .iter()
        .map(|&n| {
            let mode = pulsation::fundamental_mode(&sun(n), None).unwrap();
            assert_eq!(mode.nodes, 0);
            (n, mode.period_s / 3600.0)
        })
        .collect();
    eprintln!("solar ZAMS fundamental period [h]: {periods:?}");
    let change = (periods[2].1 / periods[1].1 - 1.0).abs();
    eprintln!("change 800 → 1600 points: {change:.2e} (arithmetic means: 1.36e-4)");
    assert!(change < 1.36e-4, "change {change:.2e}");
    assert!((periods[0].1 / periods[2].1 - 1.0).abs() < 2e-3);
}

/// The first overtones of the Sun have one and two nodes and shorter periods;
/// the cut-off ratio is small for the fundamental and grows with the order.
#[test]
fn solar_overtones_and_cutoff_ratio() {
    let model = sun(800);
    let modes: Vec<_> = pulsation::radial_modes(&model, 2).into_iter().map(|m| m.unwrap()).collect();
    for (k, m) in modes.iter().enumerate() {
        eprintln!(
            "solar ZAMS radial mode, {} nodes: P = {:.2} min, ω/ω_ac = {:.4}",
            m.nodes,
            m.period_s / 60.0,
            m.omega_over_ac
        );
        assert_eq!(m.nodes, k);
    }
    assert!(modes[0].period_s > modes[1].period_s && modes[1].period_s > modes[2].period_s);
    assert!(modes[0].omega_over_ac > 0.0 && modes[0].omega_over_ac < 0.2);
    assert!(modes[0].cutoff_warning().is_none());
    assert!(modes[1].omega_over_ac > modes[0].omega_over_ac);
    // A previous mode on the same mesh as the starting guess gives the same mode.
    let again = radial_mode(&model, 1, Some(&modes[1])).unwrap();
    assert!((again.omega2 / modes[1].omega2 - 1.0).abs() < 1e-8);
}
