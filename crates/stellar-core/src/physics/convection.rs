//! Convection: Schwarzschild criterion plus mixing-length theory (MLT), in the
//! formulation of Kippenhahn, Weigert & Weiss (2012, §7).
//!
//! With mixing length l = α H_P:
//!
//! U = 3 a c T³ / (c_P ρ² κ l²) · √(8 H_P / (g δ)),   W = ∇_rad − ∇_ad
//!
//! the cubic (ξ − U)³ + (8U/9)(ξ² − U² − W) = 0 is solved for ξ, which gives
//! ∇ = ξ² − U² + ∇_ad. Here we solve it for η = ξ − U ≥ 0:
//!
//! η³ + (8U/9)(η² + 2Uη − W) = 0,   ∇ = ∇_ad + η(η + 2U),   ∇ − ∇_e = η².

use crate::constants::{A_RAD, C};

/// Local quantities needed by MLT.
#[derive(Clone, Copy, Debug)]
pub struct MltInput {
    pub nabla_rad: f64,
    pub nabla_ad: f64,
    pub t: f64,
    pub rho: f64,
    pub p: f64,
    pub kappa: f64,
    pub cp: f64,
    pub delta: f64,
    /// Local gravity g = G m / r².
    pub g: f64,
    /// Mixing-length parameter α_MLT.
    pub alpha: f64,
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct MltResult {
    /// Actual temperature gradient d ln T / d ln P.
    pub nabla: f64,
    /// Whether the layer is convectively unstable (Schwarzschild).
    pub convective: bool,
    /// Convective velocity [cm/s].
    pub v_conv: f64,
    /// Fraction of the energy flux carried by convection.
    pub f_conv: f64,
}

pub fn mlt(i: &MltInput) -> MltResult {
    if !(i.nabla_rad > i.nabla_ad) {
        return MltResult { nabla: i.nabla_rad, convective: false, v_conv: 0.0, f_conv: 0.0 };
    }
    let hp = i.p / (i.rho * i.g);
    let lm = i.alpha * hp;
    let u = (3.0 * A_RAD * C * i.t.powi(3) / (i.cp * i.rho * i.rho * i.kappa * lm * lm)
        * (8.0 * hp / (i.g * i.delta)).sqrt())
    .clamp(0.0, 1e100);
    let w = i.nabla_rad - i.nabla_ad;
    let eta = solve_mlt_cubic(u, w);
    let nabla = i.nabla_ad + eta * (eta + 2.0 * u);
    let v_conv = eta * lm * (i.g * i.delta / (8.0 * hp)).sqrt();
    MltResult {
        nabla,
        convective: true,
        v_conv: if v_conv.is_finite() { v_conv } else { 0.0 },
        f_conv: ((i.nabla_rad - nabla) / i.nabla_rad).clamp(0.0, 1.0),
    }
}

/// Solve η³ + (8U/9)(η² + 2Uη − W) = 0 for the unique root η ≥ 0 (W > 0).
///
/// On η ≥ 0 the function is increasing and convex, and it is non-negative at
/// the upper bound η_max = √(U² + W) − U. Newton's method started from η_max
/// therefore converges monotonically from above, without overshooting.
pub fn solve_mlt_cubic(u: f64, w: f64) -> f64 {
    if !(w > 0.0) {
        return 0.0;
    }
    let c = 8.0 * u / 9.0;
    let mut eta = w / ((u * u + w).sqrt() + u);
    for _ in 0..400 {
        let f = eta * eta * eta + c * (eta * eta + 2.0 * u * eta - w);
        let df = 3.0 * eta * eta + c * (2.0 * eta + 2.0 * u);
        if !(df > 0.0) {
            break;
        }
        let new = (eta - f / df).max(0.0);
        let done = (eta - new).abs() <= 1e-15 * new.abs();
        eta = new;
        if done {
            break;
        }
    }
    eta
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cubic_root_satisfies_equation() {
        for &u in &[1e-8, 1e-3, 0.1, 1.0, 10.0, 1e4] {
            for &w in &[1e-6, 1e-2, 0.5, 10.0, 1e4] {
                let eta = solve_mlt_cubic(u, w);
                let c = 8.0 * u / 9.0;
                let f = eta.powi(3) + c * (eta * eta + 2.0 * u * eta - w);
                let scale = eta.powi(3) + c * (eta * eta + 2.0 * u * eta + w);
                assert!(f.abs() <= 1e-10 * scale, "u={u} w={w} eta={eta} f={f}");
                // ∇ must lie between ∇_ad and ∇_rad.
                let excess = eta * (eta + 2.0 * u);
                assert!(excess >= 0.0 && excess <= w * (1.0 + 1e-12));
            }
        }
    }

    fn input(alpha: f64, rho: f64) -> MltInput {
        MltInput {
            nabla_rad: 1.0,
            nabla_ad: 0.4,
            t: 1e6,
            rho,
            p: 1e14,
            kappa: 10.0,
            cp: 3e8,
            delta: 1.0,
            g: 1e5,
            alpha,
        }
    }

    #[test]
    fn efficient_convection_is_adiabatic() {
        let r = mlt(&input(1.8, 0.1));
        assert!(r.convective);
        assert!((r.nabla - 0.4).abs() < 1e-4, "∇ = {}", r.nabla);
        // F_conv/F = (∇_rad − ∇)/∇_rad = (1 − 0.4)/1.
        assert!((r.f_conv - 0.6).abs() < 1e-3, "f_conv = {}", r.f_conv);
    }

    #[test]
    fn inefficient_convection_is_radiative() {
        let r = mlt(&input(1e-4, 1e-9));
        assert!(r.convective);
        assert!((r.nabla - 1.0).abs() < 1e-2, "∇ = {}", r.nabla);
    }

    #[test]
    fn radiative_layer() {
        let mut i = input(1.8, 0.1);
        i.nabla_rad = 0.2;
        let r = mlt(&i);
        assert!(!r.convective);
        assert_eq!(r.nabla, 0.2);
        assert_eq!(r.v_conv, 0.0);
    }
}
