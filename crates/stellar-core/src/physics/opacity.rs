//! Analytic Rosseland-mean opacity.
//!
//! Combines:
//! * electron scattering: κ_es = 0.2 (1 + X)
//! * Kramers free–free:   κ_ff = 3.68e22 (1 − Z)(1 + X) ρ T^-3.5
//! * Kramers bound–free:  κ_bf = 4.34e25 / t · Z (1 + X) ρ T^-3.5, with the
//!   guillotine factor t = 2.82 [ρ (1 + X)]^0.2
//! * H⁻ (cool envelopes): κ_H⁻ = 2.5e-31 (Z / 0.02) ρ^½ T⁹
//! * a small molecular floor κ_mol = 0.1 Z
//!
//! Combination: κ = κ_mol + [1/κ_H⁻ + 1/(κ_es + κ_ff + κ_bf)]⁻¹, which stops
//! the steep H⁻ term from running away at high temperature.
//! Coefficients follow Carroll & Ostlie (2017), converted to CGS. The result is
//! qualitatively right only; tabulated opacities replace it later.

use super::Opacity;
use crate::composition::Composition;

#[derive(Debug, Clone, Copy, Default)]
pub struct AnalyticOpacity;

impl AnalyticOpacity {
    /// Individual contributions (es, ff, bf, H⁻) for diagnostics.
    pub fn components(rho: f64, t: f64, c: &Composition) -> [f64; 4] {
        let k_es = 0.2 * (1.0 + c.x);
        let t35 = t.powf(-3.5);
        let k_ff = 3.68e22 * (1.0 - c.z) * (1.0 + c.x) * rho * t35;
        let guillotine = 2.82 * (rho * (1.0 + c.x)).powf(0.2);
        let k_bf = 4.34e25 / guillotine * c.z * (1.0 + c.x) * rho * t35;
        let k_hm = 2.5e-31 * (c.z.max(1e-4) / 0.02) * rho.sqrt() * t.powi(9);
        [k_es, k_ff, k_bf, k_hm]
    }
}

impl Opacity for AnalyticOpacity {
    fn name(&self) -> &'static str {
        "analytic (Kramers + H⁻ + e⁻ scattering)"
    }

    fn kappa(&self, rho: f64, t: f64, c: &Composition) -> f64 {
        let [k_es, k_ff, k_bf, k_hm] = Self::components(rho, t, c);
        let k_mol = 0.1 * c.z;
        k_mol + 1.0 / (1.0 / k_hm + 1.0 / (k_es + k_ff + k_bf))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn positive_and_finite_on_grid() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        for lt in 30..=90 {
            for lr in -120..=30 {
                let t = 10f64.powf(lt as f64 / 10.0);
                let rho = 10f64.powf(lr as f64 / 10.0);
                let k = AnalyticOpacity.kappa(rho, t, &c);
                assert!(k.is_finite() && k > 0.0, "κ({rho}, {t}) = {k}");
            }
        }
    }

    #[test]
    fn electron_scattering_floor_in_hot_interior() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        let k = AnalyticOpacity.kappa(1.0, 5e7, &c);
        assert!((k / (0.2 * 1.7) - 1.0).abs() < 0.1, "κ = {k}");
    }
}
