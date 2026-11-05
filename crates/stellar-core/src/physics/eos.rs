//! Equation of state: fully ionized ideal gas plus radiation.

use super::{Eos, EosState};
use crate::composition::Composition;
use crate::constants::{A_RAD, R_GAS};

/// Lower bound on β = P_gas/P. It keeps the EOS finite for the unphysical
/// (radiation-dominated) states a Newton iteration can visit.
pub const MIN_BETA: f64 = 1e-3;

/// Fully ionized ideal gas + black-body radiation:
/// P = ρ k T / (μ m_u) + a T⁴ / 3.
#[derive(Debug, Clone, Copy, Default)]
pub struct IdealGasRadiation;

impl Eos for IdealGasRadiation {
    fn name(&self) -> &'static str {
        "ideal gas + radiation (fully ionized)"
    }

    fn eval(&self, p: f64, t: f64, comp: &Composition) -> EosState {
        let mu = comp.mu();
        let p_rad = A_RAD * t.powi(4) / 3.0;
        let beta = (1.0 - p_rad / p).max(MIN_BETA);
        let p_gas = beta * p;
        let rho = p_gas * mu / (R_GAS * t);
        let delta = (4.0 - 3.0 * beta) / beta;
        let b2 = beta * beta;
        let cp = R_GAS / mu * (1.5 + 12.0 * (1.0 - beta) / beta + (4.0 - 3.0 * beta).powi(2) / b2);
        let nabla_ad = 2.0 * (4.0 - 3.0 * beta) / (32.0 - 24.0 * beta - 3.0 * b2);
        let gamma1 = (32.0 - 24.0 * beta - 3.0 * b2) / (24.0 - 21.0 * beta);
        EosState { rho, p_gas, p_rad, beta, mu, nabla_ad, cp, delta, gamma1 }
    }
}

/// Invert an EOS for temperature: find T such that ρ(P, T) = ρ.
/// This works for any EOS where ρ decreases monotonically with T at fixed P.
pub fn temperature_for_density(eos: &dyn Eos, p: f64, rho: f64, comp: &Composition) -> f64 {
    let (mut lo, mut hi) = (1e2_f64.ln(), 1e10_f64.ln());
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if eos.eval(p, mid.exp(), comp).rho > rho {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-13 {
            break;
        }
    }
    (0.5 * (lo + hi)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sun() -> Composition {
        Composition::new(0.70, 0.02, 0.7).unwrap()
    }

    #[test]
    fn gas_dominated_limit() {
        let s = IdealGasRadiation.eval(1e17, 1.5e7, &sun());
        assert!(s.beta > 0.99);
        assert!((s.nabla_ad - 0.4).abs() < 2e-3);
        assert!((s.gamma1 - 5.0 / 3.0).abs() < 5e-3);
        let cp_ideal = 2.5 * R_GAS / sun().mu();
        assert!((s.cp / cp_ideal - 1.0).abs() < 0.05);
    }

    #[test]
    fn radiation_dominated_limit() {
        // β → MIN_BETA gives ∇_ad → 1/4.
        let t: f64 = 1e8;
        let p = A_RAD * t.powi(4) / 3.0 / (1.0 - 1e-3);
        let s = IdealGasRadiation.eval(p, t, &sun());
        assert!((s.nabla_ad - 0.25).abs() < 1e-3, "∇_ad = {}", s.nabla_ad);
        assert!((s.gamma1 - 4.0 / 3.0).abs() < 2e-3, "Γ₁ = {}", s.gamma1);
    }

    #[test]
    fn density_round_trip() {
        let comp = sun();
        let (p, t) = (2.3e17, 1.4e7);
        let s = IdealGasRadiation.eval(p, t, &comp);
        let p_back = s.rho * R_GAS * t / comp.mu() + A_RAD * t.powi(4) / 3.0;
        assert!((p_back / p - 1.0).abs() < 1e-12);
        let t_back = temperature_for_density(&IdealGasRadiation, p, s.rho, &comp);
        assert!((t_back / t - 1.0).abs() < 1e-9);
    }
}
