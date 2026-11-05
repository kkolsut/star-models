//! Initial model for the Henyey solver: an n = 3 polytrope (the Eddington
//! standard model), scaled to homology estimates of R and L. The near-surface
//! layers are patched with a plane-parallel hydrostatic envelope.

use super::equations::{I_L, I_LNP, I_LNR, I_LNT, NVAR};
use crate::composition::Composition;
use crate::constants::{FOUR_PI, G, L_SUN, M_SUN, R_SUN};
use crate::physics::Physics;
use crate::physics::atmosphere::{effective_temperature, solve_photosphere};
use crate::physics::eos::temperature_for_density;

/// Tabulated solution of the Lane–Emden equation
/// θ'' + (2/ξ) θ' + θⁿ = 0, with θ(0) = 1 and θ'(0) = 0.
#[derive(Clone, Debug)]
pub struct LaneEmden {
    pub n: f64,
    pub xi: Vec<f64>,
    pub theta: Vec<f64>,
    pub dtheta: Vec<f64>,
    /// First zero of θ.
    pub xi1: f64,
    /// θ'(ξ₁).
    pub dtheta1: f64,
}

impl LaneEmden {
    /// Integrate with fixed-step RK4 (step `h`), starting from the series
    /// expansion about the centre.
    pub fn solve(n: f64, h: f64) -> Self {
        let rhs = |xi: f64, th: f64, ph: f64| -> (f64, f64) { (ph, -th.max(0.0).powf(n) - 2.0 * ph / xi) };
        let mut xi: f64 = 1e-4;
        let mut th = 1.0 - xi * xi / 6.0 + n * xi.powi(4) / 120.0;
        let mut ph = -xi / 3.0 + n * xi.powi(3) / 30.0;
        let (mut xs, mut ts, mut ps) = (vec![xi], vec![th], vec![ph]);
        loop {
            let (k1t, k1p) = rhs(xi, th, ph);
            let (k2t, k2p) = rhs(xi + 0.5 * h, th + 0.5 * h * k1t, ph + 0.5 * h * k1p);
            let (k3t, k3p) = rhs(xi + 0.5 * h, th + 0.5 * h * k2t, ph + 0.5 * h * k2p);
            let (k4t, k4p) = rhs(xi + h, th + h * k3t, ph + h * k3p);
            let th_new = th + h / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
            let ph_new = ph + h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
            if th_new <= 0.0 {
                let w = th / (th - th_new);
                let xi1 = xi + w * h;
                let dtheta1 = ph + w * (ph_new - ph);
                return Self { n, xi: xs, theta: ts, dtheta: ps, xi1, dtheta1 };
            }
            xi += h;
            th = th_new;
            ph = ph_new;
            xs.push(xi);
            ts.push(th);
            ps.push(ph);
            assert!(xi < 1e3, "Lane–Emden integration did not terminate (n ≥ 5?)");
        }
    }

    /// ρ_c / ρ̄ = ξ₁ / (3 |θ'(ξ₁)|).
    pub fn central_condensation(&self) -> f64 {
        self.xi1 / (3.0 * self.dtheta1.abs())
    }

    /// P_c = W_n G M² / R⁴ with W_n = 1 / (4π (n+1) θ'(ξ₁)²).
    pub fn pressure_constant(&self) -> f64 {
        1.0 / (FOUR_PI * (self.n + 1.0) * self.dtheta1 * self.dtheta1)
    }

    /// Mass fraction q(ξ) for every tabulated point.
    pub fn mass_fractions(&self) -> Vec<f64> {
        let m1 = self.xi1 * self.xi1 * self.dtheta1.abs();
        self.xi.iter().zip(&self.dtheta).map(|(&x, &p)| x * x * (-p) / m1).collect()
    }
}

/// Homology estimate of the ZAMS radius [R☉] for mass `m` [M☉].
pub fn radius_guess(m: f64) -> f64 {
    if m < 1.3 { 0.89 * m.powf(0.8) } else { 0.89 * 1.3_f64.powf(0.8) * (m / 1.3).powf(0.6) }
}

/// Homology estimate of the ZAMS luminosity [L☉] for mass `m` [M☉].
pub fn luminosity_guess(m: f64) -> f64 {
    if m < 2.0 { 0.7 * m.powf(4.5) } else { 0.7 * 2.0_f64.powf(4.5) * (m / 2.0).powf(3.5) }
}

pub struct InitialGuess {
    pub y: Vec<f64>,
    pub l_scale: f64,
}

/// Build the starting model on mesh `q` for total mass `mass_msun`.
pub fn polytrope_guess(physics: &Physics, mass_msun: f64, comp: &Composition, q: &[f64]) -> InitialGuess {
    let le = LaneEmden::solve(3.0, 1e-3);
    let qs = le.mass_fractions();

    let mass = mass_msun * M_SUN;
    let radius = radius_guess(mass_msun) * R_SUN;
    let lum = luminosity_guess(mass_msun) * L_SUN;

    let rho_mean = mass / (FOUR_PI / 3.0 * radius.powi(3));
    let rho_c = le.central_condensation() * rho_mean;
    let p_c = le.pressure_constant() * G * mass * mass / radius.powi(4);
    let t_c = temperature_for_density(physics.eos.as_ref(), p_c, rho_c, comp);
    let teff = effective_temperature(lum, radius);
    let p_ph = solve_photosphere(physics, mass, radius, teff, comp);

    let interp = |qv: f64, ys: &[f64], y_end: f64| -> f64 {
        if qv <= qs[0] {
            return ys[0];
        }
        let last = qs.len() - 1;
        if qv >= qs[last] {
            let w = ((qv - qs[last]) / (1.0 - qs[last])).clamp(0.0, 1.0);
            return ys[last] + w * (y_end - ys[last]);
        }
        let i = qs.partition_point(|&v| v < qv).max(1);
        let w = (qv - qs[i - 1]) / (qs[i] - qs[i - 1]);
        ys[i - 1] + w * (ys[i] - ys[i - 1])
    };

    let mut r = Vec::with_capacity(q.len());
    let mut p = Vec::with_capacity(q.len());
    let mut t = Vec::with_capacity(q.len());
    for &qv in q {
        let (ri, pi, ti) = if qv >= 1.0 {
            (radius, p_ph, teff)
        } else {
            let xi = interp(qv, &le.xi, le.xi1);
            let theta = interp(qv, &le.theta, 0.0).max(0.0);
            let r_poly = if qv <= qs[0] {
                (3.0 * qv * mass / (FOUR_PI * rho_c)).cbrt()
            } else {
                (radius * xi / le.xi1).min(radius)
            };
            // Plane-parallel envelope: P ≈ G M (M − m) / (4π R⁴) + P_ph.
            let p_env = G * mass * mass * (1.0 - qv) / (FOUR_PI * radius.powi(4)) + p_ph;
            let pi = (p_c * theta.powi(4)).max(p_env);
            let ti = (t_c * (pi / p_c).powf(0.25)).max(teff);
            (r_poly, pi, ti)
        };
        r.push(ri);
        p.push(pi);
        t.push(ti);
    }
    // The radius must increase strictly outward.
    for k in 1..r.len() {
        if r[k] <= r[k - 1] {
            r[k] = r[k - 1] * (1.0 + 1e-9);
        }
    }

    // Luminosity: integrate ε over the guess, then rescale to the homology L.
    let mut cum = vec![0.0; q.len()];
    let eps: Vec<f64> = (0..q.len())
        .map(|k| {
            let rho = physics.eos.eval(p[k], t[k], comp).rho;
            physics.nuclear.eps(rho, t[k], comp).total()
        })
        .collect();
    cum[0] = eps[0] * q[0] * mass;
    for k in 1..q.len() {
        cum[k] = cum[k - 1] + 0.5 * (eps[k] + eps[k - 1]) * (q[k] - q[k - 1]) * mass;
    }
    let total = *cum.last().unwrap();
    let l_of = |k: usize| -> f64 { if total > 0.0 && total.is_finite() { lum * cum[k] / total } else { lum * q[k] } };

    let mut y = vec![0.0; q.len() * NVAR];
    for k in 0..q.len() {
        y[k * NVAR + I_LNR] = r[k].ln();
        y[k * NVAR + I_LNP] = p[k].ln();
        y[k * NVAR + I_LNT] = t[k].ln();
        y[k * NVAR + I_L] = l_of(k) / lum;
    }
    InitialGuess { y, l_scale: lum }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lane_emden_known_zeros() {
        let n1 = LaneEmden::solve(1.0, 1e-3);
        assert!((n1.xi1 - std::f64::consts::PI).abs() < 1e-6, "ξ₁(n=1) = {}", n1.xi1);
        let n3 = LaneEmden::solve(3.0, 1e-3);
        assert!((n3.xi1 - 6.896_85).abs() < 1e-4, "ξ₁(n=3) = {}", n3.xi1);
        assert!((n3.central_condensation() - 54.18).abs() < 0.05);
        assert!((n3.pressure_constant() - 11.05).abs() < 0.02);
    }
}
