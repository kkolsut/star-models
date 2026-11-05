//! Outer boundary condition: a grey Eddington photosphere.
//!
//! With T⁴(τ) = ¾ T_eff⁴ (τ + ⅔), the photosphere lies at τ = ⅔, where T = T_eff.
//! Integrating hydrostatic equilibrium dP/dτ = g/κ from τ = 0 (where only the
//! radiation pressure a T_eff⁴ / 6 remains) gives
//!
//! P_ph = ⅔ g / κ + a T_eff⁴ / 6.

use super::Physics;
use crate::composition::Composition;
use crate::constants::{A_RAD, FOUR_PI, G, SIGMA_SB};

/// Total photospheric pressure for mass m, radius r, temperature t and opacity κ.
pub fn photosphere_pressure(m: f64, r: f64, t: f64, kappa: f64) -> f64 {
    2.0 / 3.0 * G * m / (r * r * kappa) + A_RAD * t.powi(4) / 6.0
}

/// Effective temperature from L = 4π R² σ T_eff⁴.
pub fn effective_temperature(l: f64, r: f64) -> f64 {
    (l.max(0.0) / (FOUR_PI * SIGMA_SB * r * r)).powf(0.25)
}

/// Solve P_ph = ⅔ g / κ(ρ(P_ph, T), T) + a T⁴/6 by damped fixed-point
/// iteration in ln P. The damping makes the iteration converge for opacity
/// laws κ ∝ ρ^a with |a| ≤ 1.
pub fn solve_photosphere(phys: &Physics, m: f64, r: f64, teff: f64, comp: &Composition) -> f64 {
    let mut lnp = 1e5_f64.ln();
    for _ in 0..500 {
        let p = lnp.exp();
        let rho = phys.eos.eval(p, teff, comp).rho;
        let kappa = phys.opacity.kappa(rho, teff, comp);
        let target = photosphere_pressure(m, r, teff, kappa).ln();
        let new = 0.5 * (lnp + target);
        let done = (new - lnp).abs() < 1e-12;
        lnp = new;
        if done {
            break;
        }
    }
    lnp.exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{L_SUN, M_SUN, R_SUN, TEFF_SUN};

    #[test]
    fn solar_teff() {
        let t = effective_temperature(L_SUN, R_SUN);
        assert!((t - TEFF_SUN).abs() < 5.0, "T_eff = {t}");
    }

    #[test]
    fn solar_photosphere_pressure_is_reasonable() {
        let phys = Physics::analytic();
        let comp = Composition::new(0.7, 0.02, 0.7).unwrap();
        let p = solve_photosphere(&phys, M_SUN, R_SUN, TEFF_SUN, &comp);
        // The real solar photosphere has P ~ 1e5 dyn/cm²; the analytic opacity is crude.
        assert!(p > 1e3 && p < 1e7, "P_ph = {p}");
        let rho = phys.eos.eval(p, TEFF_SUN, &comp).rho;
        let kappa = phys.opacity.kappa(rho, TEFF_SUN, &comp);
        let resid = p / photosphere_pressure(M_SUN, R_SUN, TEFF_SUN, kappa) - 1.0;
        assert!(resid.abs() < 1e-8);
    }
}
