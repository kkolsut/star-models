//! Equation of state: fully ionized ideal gas plus radiation; the generic
//! derivative layer and thermodynamic-consistency scan (plan 6.1). The
//! Fermi–Dirac EOS with ionisation is in [`super::eos_fd`].

pub use super::eos_fd::{FermiDiracEos, Paczynski};
use super::{Eos, EosState};
use crate::composition::Composition;
use crate::constants::{A_RAD, H_PLANCK, K_B, M_E, M_U, R_GAS};

/// Classical electron degeneracy parameter η = ln(n_e λ_e³ / 2), λ_e the
/// thermal de Broglie wavelength (a diagnostic for EOSs without Fermi–Dirac
/// electrons).
pub fn classical_eta(rho: f64, t: f64, mu_e: f64) -> f64 {
    let n_e = rho / (mu_e * M_U);
    let lambda3 = (H_PLANCK * H_PLANCK / (2.0 * std::f64::consts::PI * M_E * K_B * t)).powf(1.5);
    (0.5 * n_e * lambda3).ln()
}

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
        // u = (3/2) P_gas/ρ + 3 P_rad/ρ = (3/2 + 3(1 − β)/β) P_gas/ρ (gas plus radiation).
        let u = (1.5 + 3.0 * (1.0 - beta) / beta) * p_gas / rho;
        EosState {
            rho,
            p_gas,
            p_rad,
            beta,
            mu,
            nabla_ad,
            cp,
            delta,
            gamma1,
            u,
            // Filled in by the profile output (see `structure::model`), not on
            // every call: the ideal EOS is evaluated millions of times.
            eta: f64::NAN,
            x_h: 1.0,
            x_he1: 0.0,
            x_he2: 1.0,
            x_h2: 0.0,
        }
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

/// Step in ln P and ln T of the derivative layer.
pub const DERIV_STEP: f64 = 1e-5;

/// Generic derivative layer (plan 6.1): the [`EosState`] of any EOS given as
/// `f(ln P, ln T) = (ln ρ, h, μ)`, by central differences in ln P and ln T.
/// `None` if `f` fails at any of the five points.
pub fn derive_state(f: &dyn Fn(f64, f64) -> Option<(f64, f64, f64)>, p: f64, t: f64) -> Option<EosState> {
    let (lp, lt) = (p.ln(), t.ln());
    let hs = DERIV_STEP;
    let (lr, h, mu) = f(lp, lt)?;
    let (lr_pp, h_pp, _) = f(lp + hs, lt)?;
    let (lr_pm, h_pm, _) = f(lp - hs, lt)?;
    let (lr_tp, h_tp, _) = f(lp, lt + hs)?;
    let (lr_tm, h_tm, _) = f(lp, lt - hs)?;
    let alpha = (lr_pp - lr_pm) / (2.0 * hs);
    let delta = -(lr_tp - lr_tm) / (2.0 * hs);
    let cp = (h_tp - h_tm) / (2.0 * hs) / t;
    let _dh_dlnp = (h_pp - h_pm) / (2.0 * hs);
    Some(EosState::from_derivatives(p, t, lr.exp(), h, alpha, delta, cp, A_RAD * t.powi(4) / 3.0, mu))
}

/// Result of a thermodynamic-consistency and stability scan (plan 6.1).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ConsistencyScan {
    /// Points evaluated and points where the EOS is undefined.
    pub points: usize,
    pub undefined: usize,
    /// Largest |ℛ| / (ε_abs + ε_rel |1 − δ|) and where (ln P, ln T).
    pub worst_ratio: f64,
    pub worst_at: (f64, f64),
    /// Largest |ℛ|.
    pub max_r: f64,
    /// Smallest c_V and α = (∂ln ρ/∂ln P)_T, and the number of points where
    /// either is not positive.
    pub min_cv: f64,
    pub min_alpha: f64,
    pub unstable: usize,
}

impl ConsistencyScan {
    /// Whether every defined point satisfies |ℛ| ≤ ε_abs + ε_rel|1 − δ|, c_V > 0, α > 0.
    pub fn passed(&self) -> bool {
        self.worst_ratio <= 1.0 && self.unstable == 0
    }
}

/// The consistency residual ℛ = ρ(∂h/∂P)_T − (1 − δ), c_V and α at one point,
/// from [`Eos::eval`] alone (differences of h = u + P/ρ at fixed T).
pub fn consistency_at(eos: &dyn Eos, p: f64, t: f64, comp: &Composition) -> Option<(f64, f64, f64, f64)> {
    let s = eos.try_eval(p, t, comp).ok()?;
    let hs: f64 = 1e-5;
    let a = eos.try_eval(p * hs.exp(), t, comp).ok()?;
    let b = eos.try_eval(p * (-hs).exp(), t, comp).ok()?;
    let dh_dp = (a.h() - b.h()) / (p * hs.exp() - p * (-hs).exp());
    let r = s.rho * dh_dp - (1.0 - s.delta);
    // α from Γ₁ = 1/(α − δ∇_ad).
    let alpha = 1.0 / s.gamma1 + s.delta * s.nabla_ad;
    let cv = s.cp - p * s.delta * s.delta / (s.rho * t * alpha);
    Some((r, cv, alpha, s.delta))
}

/// Scan `eos` over a grid of `n × n` points in ln P ∈ `lnp`, ln T ∈ `lnt`,
/// with tolerances ε_abs, ε_rel on ℛ. Points with P ≤ P_rad are skipped.
pub fn consistency_scan(
    eos: &dyn Eos,
    comp: &Composition,
    lnp: (f64, f64),
    lnt: (f64, f64),
    n: usize,
    eps_abs: f64,
    eps_rel: f64,
) -> ConsistencyScan {
    let mut out = ConsistencyScan { min_cv: f64::INFINITY, min_alpha: f64::INFINITY, ..Default::default() };
    for i in 0..n {
        for j in 0..n {
            let lp = lnp.0 + (lnp.1 - lnp.0) * i as f64 / (n - 1) as f64;
            let lt = lnt.0 + (lnt.1 - lnt.0) * j as f64 / (n - 1) as f64;
            let (p, t) = (lp.exp(), lt.exp());
            if !(p > 1.01 * A_RAD * t.powi(4) / 3.0) {
                continue;
            }
            out.points += 1;
            let Some((r, cv, alpha, delta)) = consistency_at(eos, p, t, comp) else {
                out.undefined += 1;
                continue;
            };
            let ratio = r.abs() / (eps_abs + eps_rel * (1.0 - delta).abs());
            if ratio > out.worst_ratio {
                out.worst_ratio = ratio;
                out.worst_at = (lp, lt);
            }
            out.max_r = out.max_r.max(r.abs());
            out.min_cv = out.min_cv.min(cv);
            out.min_alpha = out.min_alpha.min(alpha);
            if !(cv > 0.0 && alpha > 0.0) {
                out.unstable += 1;
            }
        }
    }
    out
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
        let u_ideal = 1.5 * R_GAS * 1.5e7 / sun().mu();
        assert!((s.u / u_ideal - 1.0).abs() < 0.05, "u = {}", s.u);
    }

    /// u = (3/2) P_gas/ρ + a T⁴/ρ, independently of how β is used.
    #[test]
    fn internal_energy_gas_plus_radiation() {
        let (p, t) = (1e16, 3e7);
        let s = IdealGasRadiation.eval(p, t, &sun());
        let expected = 1.5 * s.p_gas / s.rho + A_RAD * t.powi(4) / s.rho;
        assert!((s.u / expected - 1.0).abs() < 1e-12, "{} vs {expected}", s.u);
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

    /// Plan 6.1: the derivative layer applied to ρ(P,T), h(P,T) of the ideal
    /// EOS reproduces its closed-form c_P, ∇_ad, Γ₁ to 10⁻⁶.
    #[test]
    fn derivative_layer_reproduces_the_closed_forms() {
        let comp = sun();
        let f = |lp: f64, lt: f64| {
            let s = IdealGasRadiation.eval(lp.exp(), lt.exp(), &comp);
            Some((s.rho.ln(), s.h(), s.mu))
        };
        for (p, t) in [(1e17, 1.5e7), (1e16, 3e7), (1e5, 6000.0), (1e18, 1e8)] {
            let a = IdealGasRadiation.eval(p, t, &comp);
            let b = derive_state(&f, p, t).unwrap();
            for (name, x, y) in [
                ("ρ", a.rho, b.rho),
                ("c_P", a.cp, b.cp),
                ("∇_ad", a.nabla_ad, b.nabla_ad),
                ("Γ₁", a.gamma1, b.gamma1),
                ("δ", a.delta, b.delta),
                ("u", a.u, b.u),
            ] {
                assert!((x / y - 1.0).abs() < 1e-6, "{name} at ({p:e}, {t:e}): {x} vs {y}");
            }
        }
    }

    /// The ideal EOS gives ℛ = 0 (to the round-off of the differences), c_V > 0, α > 0.
    #[test]
    fn ideal_eos_is_consistent() {
        let s = consistency_scan(&IdealGasRadiation, &sun(), (5.0, 45.0), (8.0, 18.5), 25, 1e-9, 1e-6);
        assert!(s.points > 200, "{s:?}");
        assert!(s.passed(), "{s:?}");
    }
}
