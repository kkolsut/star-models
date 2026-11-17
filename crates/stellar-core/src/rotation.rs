//! Rotation, display only (plan 5.1, `docs/physics/08-rotation.md`, "Stage 1"):
//! the Roche surface of a uniformly rotating star and gravity darkening.
//!
//! This is post-processing of a non-rotating model. The interior is not
//! recomputed; the polar radius is taken equal to the radius of the
//! non-rotating model, and the whole mass is treated as a point mass (Roche
//! model).
//!
//! With ω = Ω/Ω_crit, x(θ) = R(θ)/R_pole and θ measured from the rotation axis:
//!
//! * Ω_crit² = 8 G M / (27 R_pole³)
//! * 1/x + (4/27) ω² x² sin²θ = 1, solved by
//!   x(θ) = 3/(ω sin θ) · cos[(π + arccos(ω sin θ))/3]
//! * g_r = −G M/R² + Ω² R sin²θ, g_θ = Ω² R sin θ cos θ, g_eff = √(g_r² + g_θ²)
//! * T_eff(θ) ∝ g_eff^β with β = 0.25 for radiative envelopes (von Zeipel 1924,
//!   eq. 37) or 0.08 for convective ones (Lucy 1967, eqs. 3 and 6), the constant
//!   fixed by ∮ σ T_eff⁴ dA = L.

use crate::constants::{G, L_SUN, M_SUN, R_SUN, SIGMA_SB};
use crate::structure::model::StarModel;

/// Gravity-darkening exponent of a radiative envelope (von Zeipel 1924).
pub const BETA_RADIATIVE: f64 = 0.25;
/// Gravity-darkening exponent of a convective envelope (Lucy 1967).
pub const BETA_CONVECTIVE: f64 = 0.08;
/// Largest ω = Ω/Ω_crit offered for display.
pub const OMEGA_MAX: f64 = 0.95;
/// A convective envelope whose base lies below this fraction of the radius
/// counts as convective for the gravity-darkening exponent.
pub const DEEP_ENVELOPE_R: f64 = 0.99;

/// x(θ) = R(θ)/R_pole of the Roche surface for ω = Ω/Ω_crit ∈ [0, 1] and
/// colatitude θ [rad].
pub fn roche_x(omega: f64, theta: f64) -> f64 {
    let a = omega.clamp(0.0, 1.0) * theta.sin().abs();
    if a < 1e-6 {
        // Series of the closed form: x = 1 + 4a²/27 + O(a⁴).
        return 1.0 + 4.0 * a * a / 27.0;
    }
    let mut x = 3.0 / a * ((std::f64::consts::PI + a.min(1.0).acos()) / 3.0).cos();
    // The closed form loses digits for small a (the cosine is near zero):
    // polish with Newton on 1/x + (4/27) a² x² − 1 = 0. Not near break-up,
    // where the root is double.
    if a < 0.5 {
        let k = 4.0 / 27.0 * a * a;
        for _ in 0..2 {
            let f = 1.0 / x + k * x * x - 1.0;
            let df = -1.0 / (x * x) + 2.0 * k * x;
            x -= f / df;
        }
    }
    x
}

/// Effective gravity (g_r, g_θ) in units of G M / R_pole² at the Roche surface.
pub fn roche_gravity(omega: f64, theta: f64) -> (f64, f64) {
    let x = roche_x(omega, theta);
    let w = 8.0 / 27.0 * omega * omega * x;
    let (s, c) = theta.sin_cos();
    (-1.0 / (x * x) + w * s * s, w * s * c)
}

/// Gravity-darkening exponent for a model: 0.08 (Lucy) if it has a convective
/// envelope reaching below 0.99 R, else 0.25 (von Zeipel).
pub fn beta_for(model: &StarModel) -> f64 {
    match model.globals.conv_env_base_r {
        Some(r) if r < DEEP_ENVELOPE_R => BETA_CONVECTIVE,
        _ => BETA_RADIATIVE,
    }
}

/// A rotating star for display: Roche shape and gravity-darkened T_eff.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RocheStar {
    /// ω = Ω/Ω_crit.
    pub omega: f64,
    /// Mass [g].
    pub mass: f64,
    /// Polar radius [cm].
    pub r_pole: f64,
    /// Luminosity [erg s⁻¹].
    pub luminosity: f64,
    /// Gravity-darkening exponent β.
    pub beta: f64,
    /// Constant of T_eff = c · g_eff^β (g_eff in units of G M / R_pole²) [K].
    c: f64,
}

/// Number of intervals in colatitude for the surface integrals (pole to equator).
const N_THETA: usize = 2000;

impl RocheStar {
    /// Rotating version of a star of mass M [g], polar radius R_pole [cm] and
    /// luminosity L [erg s⁻¹]. ω is clamped to [0, 1].
    pub fn new(mass: f64, r_pole: f64, luminosity: f64, omega: f64, beta: f64) -> Self {
        let mut star = Self { omega: omega.clamp(0.0, 1.0), mass, r_pole, luminosity, beta, c: 1.0 };
        // ∮ σ T⁴ dA = σ c⁴ ∮ g^{4β} dA = L
        let flux_integral = star.surface_integral(|g| g.powf(4.0 * beta));
        star.c = (luminosity / (SIGMA_SB * flux_integral)).powf(0.25);
        star
    }

    /// Rotating version of a model: R_pole = R, the model's M and L, β from
    /// its envelope ([`beta_for`]) unless given.
    pub fn of_model(model: &StarModel, omega: f64, beta: Option<f64>) -> Self {
        let g = &model.globals;
        Self::new(
            g.mass * M_SUN,
            g.radius * R_SUN,
            g.luminosity * L_SUN,
            omega,
            beta.unwrap_or_else(|| beta_for(model)),
        )
    }

    /// G M / R_pole² [cm s⁻²].
    fn g0(&self) -> f64 {
        G * self.mass / (self.r_pole * self.r_pole)
    }

    /// Critical angular velocity Ω_crit = √(8 G M / (27 R_pole³)) [rad s⁻¹].
    pub fn omega_crit(&self) -> f64 {
        (8.0 * G * self.mass / (27.0 * self.r_pole.powi(3))).sqrt()
    }

    /// Angular velocity Ω [rad s⁻¹].
    pub fn angular_velocity(&self) -> f64 {
        self.omega * self.omega_crit()
    }

    /// x(θ) = R(θ)/R_pole.
    pub fn x(&self, theta: f64) -> f64 {
        roche_x(self.omega, theta)
    }

    /// Radius at colatitude θ [cm].
    pub fn radius(&self, theta: f64) -> f64 {
        self.r_pole * self.x(theta)
    }

    /// Equatorial radius [cm].
    pub fn r_equator(&self) -> f64 {
        self.radius(std::f64::consts::FRAC_PI_2)
    }

    /// Equatorial velocity Ω R_eq [cm s⁻¹].
    pub fn v_equator(&self) -> f64 {
        self.angular_velocity() * self.r_equator()
    }

    /// Effective gravity |g_eff| at colatitude θ [cm s⁻²].
    pub fn g_eff(&self, theta: f64) -> f64 {
        let (gr, gt) = roche_gravity(self.omega, theta);
        gr.hypot(gt) * self.g0()
    }

    /// Local effective temperature T_eff(θ) [K].
    pub fn teff(&self, theta: f64) -> f64 {
        let (gr, gt) = roche_gravity(self.omega, theta);
        self.c * gr.hypot(gt).powf(self.beta)
    }

    /// ∮ f(g) dA over the surface [cm²] for g = g_eff / (G M / R_pole²): the
    /// midpoint rule in θ from pole to equator, doubled. The area element is
    /// dA = 2π R sin θ √(R² + (dR/dθ)²) dθ with dx/dθ = −x g_θ / g_r (the
    /// surface is an equipotential, normal to g_eff).
    fn surface_integral(&self, f: impl Fn(f64) -> f64) -> f64 {
        let h = std::f64::consts::FRAC_PI_2 / N_THETA as f64;
        let mut sum = 0.0;
        for i in 0..N_THETA {
            let theta = (i as f64 + 0.5) * h;
            let x = self.x(theta);
            let (gr, gt) = roche_gravity(self.omega, theta);
            let dx = if gr != 0.0 { -x * gt / gr } else { 0.0 };
            let da = std::f64::consts::TAU * x * theta.sin() * (x * x + dx * dx).sqrt();
            sum += f(gr.hypot(gt)) * da * h;
        }
        2.0 * sum * self.r_pole * self.r_pole
    }

    /// Surface area [cm²].
    pub fn area(&self) -> f64 {
        self.surface_integral(|_| 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn no_rotation_is_a_sphere() {
        for i in 0..=20 {
            let theta = PI * i as f64 / 20.0;
            assert!((roche_x(0.0, theta) - 1.0).abs() < 1e-15);
        }
        let star = RocheStar::new(M_SUN, R_SUN, L_SUN, 0.0, BETA_RADIATIVE);
        assert!((star.area() / (4.0 * PI * R_SUN * R_SUN) - 1.0).abs() < 1e-6);
        // Uniform T_eff = (L / 4πR²σ)^¼.
        let t = (L_SUN / (4.0 * PI * R_SUN * R_SUN * SIGMA_SB)).powf(0.25);
        for theta in [0.0, 0.7, FRAC_PI_2] {
            assert!((star.teff(theta) / t - 1.0).abs() < 1e-6);
        }
    }

    #[test]
    fn break_up_gives_three_halves() {
        assert!((roche_x(1.0, FRAC_PI_2) - 1.5).abs() < 1e-12);
        assert!((roche_x(1.0, 0.0) - 1.0).abs() < 1e-15);
        // At break-up the effective gravity vanishes at the equator.
        let (gr, gt) = roche_gravity(1.0, FRAC_PI_2);
        assert!(gr.abs() < 1e-12 && gt.abs() < 1e-12);
    }

    #[test]
    fn closed_form_solves_the_roche_equation() {
        for &omega in &[1e-7, 1e-3, 0.3, 0.7, 0.95, 1.0] {
            for i in 0..=30 {
                let theta = PI * i as f64 / 30.0;
                let x = roche_x(omega, theta);
                let s = theta.sin();
                let res = 1.0 / x + 4.0 / 27.0 * omega * omega * x * x * s * s - 1.0;
                assert!(res.abs() < 1e-12, "ω = {omega}, θ = {theta}: residual {res}");
            }
        }
    }

    /// ∮ σ T_eff⁴ dA = L, checked with an independent quadrature: Simpson's
    /// rule in μ = cos θ, with dR/dθ from finite differences.
    #[test]
    fn luminosity_is_conserved() {
        for &(omega, beta) in &[(0.0, 0.25), (0.5, 0.25), (0.8, 0.08), (0.95, 0.25), (0.95, 0.08)] {
            let star = RocheStar::new(20.0 * M_SUN, 5.0 * R_SUN, 4e4 * L_SUN, omega, beta);
            let n = 4000; // even
            let mut sum = 0.0;
            for i in 0..=n {
                let mu = i as f64 / n as f64; // pole (1) to equator (0), one hemisphere
                let theta = mu.acos().clamp(1e-9, FRAC_PI_2);
                let h = 1e-6;
                let r = star.radius(theta);
                let (lo, hi) = (theta - h, (theta + h).min(FRAC_PI_2));
                let dr = (star.radius(hi) - star.radius(lo)) / (hi - lo);
                // dA = 2π R √(R² + R'²) sin θ dθ = 2π R √(R² + R'²) dμ
                let f = SIGMA_SB * star.teff(theta).powi(4) * std::f64::consts::TAU * r * (r * r + dr * dr).sqrt();
                let w = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                sum += w * f;
            }
            let l = 2.0 * sum / (3.0 * n as f64);
            assert!((l / star.luminosity - 1.0).abs() < 2e-4, "ω = {omega}, β = {beta}: {}", l / star.luminosity);
        }
    }

    #[test]
    fn poles_are_hotter_than_the_equator() {
        let star = RocheStar::new(3.0 * M_SUN, 2.0 * R_SUN, 80.0 * L_SUN, 0.9, BETA_RADIATIVE);
        assert!(star.teff(0.0) > star.teff(FRAC_PI_2) * 1.2);
        // ω = 0.9: x(π/2) = (3/0.9) cos[(π + arccos 0.9)/3] = 1.21550.
        assert!((star.r_equator() / star.r_pole - 1.21550).abs() < 1e-5, "{}", star.r_equator() / star.r_pole);
        // Ω_crit² R_pole³ = 8GM/27.
        let w = star.omega_crit();
        assert!((w * w * star.r_pole.powi(3) / (8.0 * G * star.mass / 27.0) - 1.0).abs() < 1e-12);
    }
}
