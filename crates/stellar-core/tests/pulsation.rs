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

/// Independent reference for the polytrope modes: the same adiabatic wave
/// equation written in Lane–Emden variables and solved by shooting, with no
/// mesh and no surface cut. With r = a z, θ(z) the Lane–Emden function of index
/// n, z₁ its first zero and w₁ = −z₁²θ′(z₁), the equations of `pulsation.rs`
/// become
///
///   dξ/dz = −(3ξ + p/Γ₁)/z,
///   dp/dz = −(n + 1)(θ′/θ) [(4 + Ω² (z/z₁)³ w₁/(−z²θ′)) ξ + p],
///
/// regular at the centre (3ξ + p/Γ₁ = 0) and at the surface (p + (4 + Ω²)ξ = 0).
/// Both regular solutions are integrated (RK4) to z₁/2, from the centre outward
/// and from the surface inward: the irregular solutions decay in both
/// directions. Eigenvalues are the zeros of ξ_out p_in − ξ_in p_out there.
struct Shooter {
    n: f64,
    gamma1: f64,
    z1: f64,
    dtheta1: f64,
}

impl Shooter {
    fn new(n: f64, gamma1: f64) -> Self {
        let (le, z1) = lane_emden(n);
        let (zl, _, dl) = le[le.len() - 1];
        // θ'' = −2θ'/z − θⁿ at the surface, where θⁿ is 0, or 1 for n = 0.
        let surf = if n == 0.0 { 1.0 } else { 0.0 };
        let dtheta1 = dl + (-2.0 * dl / zl - surf) * (z1 - zl);
        Self { n, gamma1, z1, dtheta1 }
    }

    fn derivs(&self, om2: f64, z: f64, y: [f64; 4]) -> [f64; 4] {
        let [th, dth, xi, p] = y;
        let w1 = -self.z1 * self.z1 * self.dtheta1;
        let ddth = -th.max(0.0).powf(self.n) - 2.0 * dth / z;
        let dxi = -(3.0 * xi + p / self.gamma1) / z;
        let coeff = 4.0 + om2 * (z / self.z1).powi(3) * w1 / (-z * z * dth);
        let dp = -(self.n + 1.0) * (dth / th) * (coeff * xi + p);
        [dth, ddth, dxi, dp]
    }

    /// RK4 from z to z_end; the step adapts near the surface, where θ′/θ grows.
    fn integrate(&self, om2: f64, mut z: f64, mut y: [f64; 4], z_end: f64) -> [f64; 4] {
        let dir = (z_end - z).signum();
        while (z_end - z) * dir > 0.0 {
            let h = (1e-4_f64).min(0.01 * (self.z1 - z).abs().max(1e-12)).min((z_end - z).abs()) * dir;
            let add = |y: [f64; 4], k: [f64; 4], f: f64| [0, 1, 2, 3].map(|i| y[i] + f * k[i]);
            let k1 = self.derivs(om2, z, y);
            let k2 = self.derivs(om2, z + h / 2.0, add(y, k1, h / 2.0));
            let k3 = self.derivs(om2, z + h / 2.0, add(y, k2, h / 2.0));
            let k4 = self.derivs(om2, z + h, add(y, k3, h));
            y = [0, 1, 2, 3].map(|i| y[i] + h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
            z += h;
        }
        y
    }

    fn mismatch(&self, om2: f64) -> f64 {
        let (n, zm) = (self.n, 0.5 * self.z1);
        let z0: f64 = 1e-4;
        let centre =
            [1.0 - z0 * z0 / 6.0 + n * z0.powi(4) / 120.0, -z0 / 3.0 + n * z0.powi(3) / 30.0, 1.0, -3.0 * self.gamma1];
        let out = self.integrate(om2, z0, centre, zm);
        let d = 1e-8 * self.z1;
        let ddth1 = -2.0 * self.dtheta1 / self.z1 - if n == 0.0 { 1.0 } else { 0.0 };
        let surface = [-self.dtheta1 * d, self.dtheta1 - ddth1 * d, 1.0, -(4.0 + om2)];
        let inn = self.integrate(om2, self.z1 - d, surface, zm);
        out[2] * inn[3] - inn[2] * out[3]
    }

    /// The first `count` eigenvalues Ω², by scanning and bisection.
    fn eigenvalues(&self, count: usize) -> Vec<f64> {
        let mut roots = Vec::new();
        let (mut lo, mut f_lo) = (0.05, self.mismatch(0.05));
        while roots.len() < count && lo < 200.0 {
            let hi = lo * 1.04;
            let f_hi = self.mismatch(hi);
            if f_lo.signum() != f_hi.signum() {
                let (mut a, mut b, mut fa) = (lo, hi, f_lo);
                for _ in 0..60 {
                    let m = 0.5 * (a + b);
                    let fm = self.mismatch(m);
                    if fm.signum() == fa.signum() {
                        a = m;
                        fa = fm;
                    } else {
                        b = m;
                    }
                }
                roots.push(0.5 * (a + b));
            }
            (lo, f_lo) = (hi, f_hi);
        }
        roots
    }
}

/// The shooting reference reproduces the uniform sphere (n = 0) exactly:
/// Ω²_k = Γ₁ (k + 1)(2k + 3) − 4 (Pekeris 1938).
#[test]
fn shooting_reference_reproduces_the_uniform_sphere() {
    let g1 = 5.0 / 3.0;
    let w = Shooter::new(0.0, g1).eigenvalues(3);
    for (k, w) in w.iter().enumerate() {
        let exact = g1 * ((k + 1) * (2 * k + 3)) as f64 - 4.0;
        assert!((w / exact - 1.0).abs() < 1e-6, "k = {k}: Ω² = {w}, exact {exact}");
    }
}

/// Radial modes of the n = 1.5 and n = 3 polytropes (Γ₁ = 5/3): the Henyey
/// solver on a 3000-point model (surface cut at θ = 10⁻⁵) against the shooting
/// reference above (agreement 4e-5 for n = 1.5 and 1.1e-4 for n = 3), for the fundamental and the first two overtones, with the
/// right node counts.
///
/// Literature values for comparison: Ω² ≈ 2.712 (n = 1.5) and 9.26 (n = 3) for
/// the fundamental, quoted from memory of the classical tables (Hurley, Roberts
/// & Wright 1966; Cox 1980; for n = 3 Schwarzschild 1941; precise values in
/// Christensen-Dalsgaard & Mullan 1994) and not yet compared with a source
/// [check source]. They are printed, not tested; the shooting reference gives
/// 2.70587 and 9.25471, so the remembered 2.712 is probably wrong.
#[test]
fn polytrope_radial_modes() {
    for (n, literature) in [(1.5, 2.712), (3.0, 9.26)] {
        let reference = Shooter::new(n, 5.0 / 3.0).eigenvalues(3);
        let model = polytrope(n, 3000);
        let modes: Vec<_> = (0..3).map(|k| radial_mode(&model, k, None).unwrap()).collect();
        eprintln!(
            "polytrope n = {n}: shooting Ω² = {:.5}, {:.5}, {:.5}; Henyey {:.5}, {:.5}, {:.5}; literature (unchecked) {literature}",
            reference[0], reference[1], reference[2], modes[0].omega2, modes[1].omega2, modes[2].omega2
        );
        for (k, m) in modes.iter().enumerate() {
            assert_eq!(m.nodes, k);
            assert!((m.omega2 / reference[k] - 1.0).abs() < 3e-4, "n = {n}, k = {k}: {} vs {}", m.omega2, reference[k]);
        }
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
