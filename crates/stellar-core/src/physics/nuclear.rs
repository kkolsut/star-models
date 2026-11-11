//! Analytic nuclear energy generation rates (Kippenhahn & Weigert 1990 forms, CGS).
//!
//! * pp chain:  ε_pp  = 2.38e6 ρ X² f₁₁ g₁₁ T₆^(−2/3) exp(−33.80 T₆^(−1/3)),
//!   g₁₁ = 1 + 0.0123 T₆^(1/3) + 0.0109 T₆^(2/3) + 0.0009 T₆ (rate correction).
//!   The pp II/III branching factor ψ of the source formula is taken as ψ = 1.
//! * CNO cycle: ε_CNO = 8.67e27 ρ X X_CNO f₁₄ C T₆^(−2/3) exp(−152.28 T₆^(−1/3))
//! * triple-α:  ε_3α  = 5.09e11 ρ² Y³ f T₈^(−3) exp(−44.027 / T₈)
//!
//! Each includes a weak-screening factor f = exp(0.188 Z₁Z₂ ζ ρ^½ T₆^(−3/2)).

use super::{Nuclear, NuclearRates};
use crate::composition::Composition;

#[derive(Debug, Clone, Copy, Default)]
pub struct AnalyticNuclear;

/// Cap on the weak-screening exponent. The weak-screening approximation is
/// only valid for small values; the cap also guards against overflow in
/// unphysical intermediate states.
const MAX_SCREENING_EXPONENT: f64 = 2.0;

impl Nuclear for AnalyticNuclear {
    fn name(&self) -> &'static str {
        "analytic pp + CNO + 3α (weak screening)"
    }

    fn eps(&self, rho: f64, t: f64, c: &Composition) -> NuclearRates {
        if !(t > 0.0 && rho > 0.0) {
            return NuclearRates::default();
        }
        let t6 = t * 1e-6;
        let t613 = t6.cbrt();
        let t623 = t613 * t613;

        // ζ² = Σ X_i Z_i (Z_i + 1) / A_i (metals approximated by CNO).
        let zeta = (2.0 * c.x + 1.5 * c.y + 4.0 * c.z).sqrt();
        let screen_base = 0.188 * zeta * rho.sqrt() * t6.powf(-1.5);
        let screen = |z1z2: f64| (z1z2 * screen_base).min(MAX_SCREENING_EXPONENT).exp();

        // g₁₁: correction to the p + p rate (not the branching factor ψ, which is 1 here).
        let g11 = 1.0 + 0.0123 * t613 + 0.0109 * t623 + 0.0009 * t6;
        let pp = 2.38e6 * rho * c.x * c.x * screen(1.0) * g11 * (-33.80 / t613).exp() / t623;

        let c_cno = (1.0 + 0.0027 * t613 - 0.00778 * t623 - 0.000149 * t6).max(0.0);
        let cno = 8.67e27 * rho * c.x * c.x_cno * screen(7.0) * c_cno * (-152.28 / t613).exp() / t623;

        let t8 = t * 1e-8;
        let triple_alpha =
            5.09e11 * rho * rho * c.y.powi(3) * screen(4.0) * screen(8.0) * (-44.027 / t8).exp() / t8.powi(3);

        NuclearRates { pp, cno, triple_alpha }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn present_solar_centre_order_of_magnitude() {
        // Present-day solar centre: T ≈ 15.7 MK, ρ ≈ 150 g/cm³, X ≈ 0.34.
        let c = Composition::new(0.34, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear.eps(150.0, 1.57e7, &c);
        assert!(e.pp > 3.0 && e.pp < 50.0, "ε_pp = {}", e.pp);
        assert!(e.cno < e.pp, "pp should dominate in the Sun");
        assert!(e.triple_alpha < 1e-20);
    }

    #[test]
    fn cno_dominates_in_hot_cores() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear.eps(20.0, 2.8e7, &c);
        assert!(e.cno > e.pp, "pp {} vs cno {}", e.pp, e.cno);
    }

    #[test]
    fn finite_and_nonnegative() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        for lt in 30..=95 {
            for lr in -120..=60 {
                let t = 10f64.powf(lt as f64 / 10.0);
                let rho = 10f64.powf(lr as f64 / 10.0);
                let e = AnalyticNuclear.eps(rho, t, &c);
                assert!(e.total().is_finite() && e.total() >= 0.0);
            }
        }
    }
}
