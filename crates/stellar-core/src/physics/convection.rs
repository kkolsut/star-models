//! Convection: Schwarzschild criterion plus mixing-length theory (MLT), in the
//! formulation of Kippenhahn, Weigert & Weiss (2012, §7).
//!
//! With mixing length l = α H̃_P, where the scale height is bounded so that it
//! stays finite at the centre (Paxton et al. 2011, ApJS 192, 3, §5.1):
//!
//! H̃_P = min(P / (ρ g), √(P / (G ρ²)))
//!
//! The ordinary H_P = P/(ρg) is kept wherever it comes from the temperature
//! gradient itself (in U and in the convective velocity below).
//!
//!
//! U = 3 a c T³ / (c_P ρ² κ l²) · √(8 H_P / (g δ)),   W = ∇_rad − ∇_ad
//!
//! the cubic (ξ − U)³ + (8U/9)(ξ² − U² − W) = 0 is solved for ξ, which gives
//! ∇ = ξ² − U² + ∇_ad. Here we solve it for η = ξ − U ≥ 0:
//!
//! η³ + (8U/9)(η² + 2Uη − W) = 0,   ∇ = ∇_ad + η(η + 2U),   ∇ − ∇_e = η².
//!
//! Generalised form used in the atmosphere (`docs/physics/06`, plan 3.2): with
//! the Hopf factor A = 1 + q'(τ) of the background radiative transport and the
//! element-loss parameter U_e of optically thin elements,
//!
//! η³ + (8U/9A)(η² + 2U_e η − W_A) = 0,   ∇ = ∇_ad + η² + 2U_e η,   W_A = A∇_rad − ∇_ad,
//!
//! U_e = 3acT³κ / [c_P (3 + τ_e²)] · √(8H_P/(gδ)),   τ_e = κρℓ
//!
//! (Henyey, Vardya & Bodenheimer 1965 with y = 1/3). A = 1 and opaque elements
//! (U_e = U) give the cubic above.

use crate::constants::{A_RAD, C, G};

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
    /// Hopf factor A = 1 + q'(τ) of the background radiative transport
    /// (1 outside the atmosphere).
    pub hopf_a: f64,
    /// Radiative losses of optically thin elements (U_e in place of U in the
    /// element relation); false = the opaque-element formula of page 5.
    pub thin_elements: bool,
}

impl MltInput {
    /// Input with A = 1 and opaque elements (the interior formula).
    #[allow(clippy::too_many_arguments)]
    pub fn interior(
        nabla_rad: f64,
        nabla_ad: f64,
        t: f64,
        rho: f64,
        p: f64,
        kappa: f64,
        cp: f64,
        delta: f64,
        g: f64,
        alpha: f64,
    ) -> Self {
        Self { nabla_rad, nabla_ad, t, rho, p, kappa, cp, delta, g, alpha, hopf_a: 1.0, thin_elements: false }
    }
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

/// Scale height for the mixing length, min(P/(ρg), √(P/(Gρ²))): finite as g → 0.
pub fn bounded_scale_height(p: f64, rho: f64, g: f64) -> f64 {
    (p / (rho * g)).min((p / (G * rho * rho)).sqrt())
}

pub fn mlt(i: &MltInput) -> MltResult {
    let nabla_rad = i.hopf_a * i.nabla_rad;
    if !(nabla_rad > i.nabla_ad) {
        return MltResult { nabla: nabla_rad, convective: false, v_conv: 0.0, f_conv: 0.0 };
    }
    let hp = i.p / (i.rho * i.g);
    let lm = i.alpha * bounded_scale_height(i.p, i.rho, i.g);
    let root = (8.0 * hp / (i.g * i.delta)).sqrt();
    let u = (3.0 * A_RAD * C * i.t.powi(3) / (i.cp * i.rho * i.rho * i.kappa * lm * lm) * root).clamp(0.0, 1e100);
    let ue = if i.thin_elements { element_loss(i.t, i.rho, i.kappa, i.cp, lm, root) } else { u };
    let w = nabla_rad - i.nabla_ad;
    let eta = solve_mlt_cubic_general(u / i.hopf_a, ue, w);
    let nabla = i.nabla_ad + eta * (eta + 2.0 * ue);
    let v_conv = eta * lm * (i.g * i.delta / (8.0 * hp)).sqrt();
    MltResult {
        nabla,
        convective: true,
        v_conv: if v_conv.is_finite() { v_conv } else { 0.0 },
        f_conv: ((nabla_rad - nabla) / nabla_rad).clamp(0.0, 1.0),
    }
}

/// U_e = 3acT³κ / [c_P (3 + τ_e²)] · √(8H_P/(gδ)) with τ_e = κρℓ; `root` is
/// √(8H_P/(gδ)). Written so that it never multiplies a diverging U by a
/// vanishing factor.
pub fn element_loss(t: f64, rho: f64, kappa: f64, cp: f64, lm: f64, root: f64) -> f64 {
    let tau_e = kappa * rho * lm;
    (3.0 * A_RAD * C * t.powi(3) * kappa / (cp * (3.0 + tau_e * tau_e)) * root).clamp(0.0, 1e100)
}

/// Solve η³ + (8U/9)(η² + 2Uη − W) = 0 for the unique root η ≥ 0 (W > 0).
///
/// On η ≥ 0 the function is increasing and convex, and it is non-negative at
/// the upper bound η_max = √(U² + W) − U. Newton's method started from η_max
/// therefore converges monotonically from above, without overshooting.
pub fn solve_mlt_cubic(u: f64, w: f64) -> f64 {
    solve_mlt_cubic_general(u, u, w)
}

/// Solve η³ + (8U_c/9)(η² + 2U_e η − W) = 0 for the root η ≥ 0 (W > 0), with
/// U_c = U/A the coefficient and U_e the element-loss parameter. The same
/// argument as for [`solve_mlt_cubic`] applies, with η_max = √(U_e² + W) − U_e.
pub fn solve_mlt_cubic_general(uc: f64, ue: f64, w: f64) -> f64 {
    if !(w > 0.0) {
        return 0.0;
    }
    let c = 8.0 * uc / 9.0;
    let u = ue;
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
            hopf_a: 1.0,
            thin_elements: false,
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

    /// Plan 1.4: at g → 0 the mixing length stays finite (√(P/Gρ²)) and MLT
    /// returns finite results with ∇_ad ≤ ∇ ≤ ∇_rad. (As the buoyancy vanishes,
    /// convection becomes inefficient and ∇ → ∇_rad; the innermost mesh point
    /// sits at m₀ > 0, where it is still adiabatic.)
    #[test]
    fn bounded_scale_height_at_the_centre() {
        let (p, rho) = (1e17, 100.0);
        let cap = (p / (G * rho * rho)).sqrt();
        assert_eq!(bounded_scale_height(p, rho, 0.0), cap);
        assert_eq!(bounded_scale_height(p, rho, 1e-3), cap);
        // Far from the centre the ordinary scale height is used.
        assert_eq!(bounded_scale_height(p, rho, 1e6), p / (rho * 1e6));
        for g in [1e-30, 1e-10, 1e-3, 1.0] {
            let mut i = input(1.8, 0.1);
            i.g = g;
            let r = mlt(&i);
            assert!(r.convective && r.nabla.is_finite() && r.v_conv.is_finite(), "g = {g}: {r:?}");
            assert!(r.nabla >= 0.4 && r.nabla <= 1.0, "g = {g}: ∇ = {}", r.nabla);
        }
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

    /// Plan 3.2: A = 1 and τ_e → ∞ reproduce the cubic of page 5.
    #[test]
    fn general_form_reduces_to_page_5() {
        for (alpha, rho) in [(1.8, 0.1), (1e-2, 1e-7), (1.0, 1e-5)] {
            let a = mlt(&input(alpha, rho));
            let mut i = input(alpha, rho);
            i.thin_elements = true;
            i.kappa = 1e8; // τ_e = κρℓ ≫ 1
            let mut j = input(alpha, rho);
            j.kappa = 1e8;
            let (b, c) = (mlt(&i), mlt(&j));
            assert!((b.nabla - c.nabla).abs() < 1e-10 * c.nabla, "{b:?} vs {c:?}");
            assert!(a.convective);
        }
    }

    /// U_e → U for τ_e ≫ 1 and → U τ_e²/3 for τ_e ≪ 1.
    #[test]
    fn element_loss_limits() {
        let (t, cp, root): (f64, f64, f64) = (6000.0, 3e8, 1e-3);
        let u_of =
            |kappa: f64, rho: f64, lm: f64| 3.0 * A_RAD * C * t.powi(3) / (cp * rho * rho * kappa * lm * lm) * root;
        let lm = 1e7;
        // Opaque: κρℓ = 1e4.
        let (kappa, rho) = (10.0, 1e-4);
        let r = element_loss(t, rho, kappa, cp, lm, root) / u_of(kappa, rho, lm);
        assert!((r - 1.0).abs() < 1e-7, "{r}");
        // Transparent: κρℓ = 1e-3.
        let (kappa, rho) = (1e-3, 1e-7);
        let tau_e = kappa * rho * lm;
        let r = element_loss(t, rho, kappa, cp, lm, root) / (u_of(kappa, rho, lm) * tau_e * tau_e / 3.0);
        assert!((r - 1.0).abs() < 1e-6, "{r}");
    }

    /// ∇ is continuous at the onset of convection: as A∇_rad → ∇_ad from above,
    /// ∇ → A∇_rad, for any A and U_e.
    #[test]
    fn continuous_at_onset_with_hopf_factor() {
        for a in [0.8, 1.0, 1.7] {
            for thin in [false, true] {
                let mut i = input(1.8, 1e-6);
                i.hopf_a = a;
                i.thin_elements = thin;
                i.kappa = 1e-2;
                i.nabla_rad = 0.4 / a * (1.0 + 1e-9);
                let r = mlt(&i);
                assert!(r.convective);
                assert!((r.nabla - a * i.nabla_rad).abs() < 1e-8, "A = {a}: {} vs {}", r.nabla, a * i.nabla_rad);
                i.nabla_rad = 0.4 / a * (1.0 - 1e-9);
                let r = mlt(&i);
                assert!(!r.convective && (r.nabla - a * i.nabla_rad).abs() < 1e-15);
            }
        }
    }

    #[test]
    fn general_cubic_root() {
        for &uc in &[1e-6, 0.1, 10.0] {
            for &ue in &[1e-8, 1e-3, 1.0, 100.0] {
                for &w in &[1e-4, 0.3, 50.0] {
                    let eta = solve_mlt_cubic_general(uc, ue, w);
                    let c = 8.0 * uc / 9.0;
                    let f = eta.powi(3) + c * (eta * eta + 2.0 * ue * eta - w);
                    let scale = eta.powi(3) + c * (eta * eta + 2.0 * ue * eta + w);
                    assert!(f.abs() <= 1e-10 * scale, "uc={uc} ue={ue} w={w}");
                }
            }
        }
    }
}
