//! Linear adiabatic radial pulsations, and where stars pulsate.
//!
//! The radial displacement ξ = δr/r and the Lagrangian pressure perturbation
//! p = δP/P of a mode with angular frequency ω obey (in mass coordinates)
//!
//!   dξ/dm = −(3ξ + p/Γ₁) / (4π r³ ρ)
//!   dp/dm = G m / (4π r⁴ P) · [ (4 + ω² r³/(G m)) ξ + p ]
//!
//! with regularity at the centre, 3ξ + p/Γ₁ = 0, and at the surface
//! p + (4 + Ω²) ξ = 0 and ξ = 1 (normalisation), where Ω² = ω² R³/(G M).
//! The eigenvalue Ω² is solved for as an extra unknown with the generic
//! Henyey solver. The fundamental mode is the solution without nodes.
//!
//! Whether a star actually pulsates (κ mechanism) needs a non-adiabatic
//! analysis; here an approximate classical instability strip in the HR
//! diagram decides it.

use serde::{Deserialize, Serialize};

use crate::constants::{FOUR_PI, G, M_SUN, R_SUN};
use crate::solver::{BvpProblem, HenyeyNewton, NewtonOptions};
use crate::structure::model::StarModel;

/// A radial pulsation mode.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RadialMode {
    /// Period [s].
    pub period_s: f64,
    /// Dimensionless eigenvalue Ω² = ω² R³ / (G M).
    pub omega2: f64,
    /// Number of nodes of ξ (0 for the fundamental mode).
    pub nodes: usize,
    /// ξ = δr/r on the model mesh, normalised to 1 at the surface.
    pub xi: Vec<f64>,
    /// δP/P on the model mesh.
    #[serde(skip)]
    pub dp: Vec<f64>,
    /// Pulsation constant Q = Π √(ρ̄/ρ̄☉) [days].
    pub q_days: f64,
}

impl RadialMode {
    pub fn period_days(&self) -> f64 {
        self.period_s / 86400.0
    }
}

struct Lawe<'a> {
    m: &'a [f64],
    r: &'a [f64],
    p: &'a [f64],
    rho: &'a [f64],
    gamma1: &'a [f64],
    mass: f64,
    radius: f64,
}

impl BvpProblem for Lawe<'_> {
    fn nvar(&self) -> usize {
        3
    }
    fn n_inner(&self) -> usize {
        1
    }
    fn npts(&self) -> usize {
        self.m.len()
    }
    fn inner_bc(&self, y: &[f64], out: &mut [f64]) {
        out[0] = 3.0 * y[0] + y[1] / self.gamma1[0];
    }
    fn outer_bc(&self, y: &[f64], out: &mut [f64]) {
        out[0] = y[1] + (4.0 + y[2]) * y[0];
        out[1] = y[0] - 1.0;
    }
    fn interval(&self, k: usize, a: &[f64], b: &[f64], out: &mut [f64]) {
        let dm = self.m[k + 1] - self.m[k];
        let mm = 0.5 * (self.m[k] + self.m[k + 1]);
        let rm = 0.5 * (self.r[k] + self.r[k + 1]);
        let pm = 0.5 * (self.p[k] + self.p[k + 1]);
        let rhom = 0.5 * (self.rho[k] + self.rho[k + 1]);
        let g1 = 0.5 * (self.gamma1[k] + self.gamma1[k + 1]);
        let (xi, dp, w) = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]));
        out[0] = (b[0] - a[0]) + dm * (3.0 * xi + dp / g1) / (FOUR_PI * rm.powi(3) * rhom);
        let coeff = 4.0 + w * (self.mass / mm) * (rm / self.radius).powi(3);
        out[1] = (b[1] - a[1]) - dm * G * mm / (FOUR_PI * rm.powi(4) * pm) * (coeff * xi + dp);
        out[2] = b[2] - a[2];
    }
    fn max_correction(&self, var: usize) -> f64 {
        if var == 2 { 3.0 } else { 2.0 }
    }
}

/// Solve for a radial mode on given profiles (CGS), starting from a guess.
#[allow(clippy::too_many_arguments)]
pub fn solve_mode(
    m: &[f64],
    r: &[f64],
    p: &[f64],
    rho: &[f64],
    gamma1: &[f64],
    guess_xi: &[f64],
    guess_dp: &[f64],
    guess_omega2: f64,
) -> Result<(Vec<f64>, Vec<f64>, f64), String> {
    let n = m.len();
    let mass = m[n - 1];
    let radius = r[n - 1];
    let problem = Lawe { m, r, p, rho, gamma1, mass, radius };
    let y0: Vec<f64> = (0..n).flat_map(|k| [guess_xi[k], guess_dp[k], guess_omega2]).collect();
    let mut nr = HenyeyNewton::new(&problem, y0, NewtonOptions { max_iter: 60, tol: 1e-9 });
    nr.run()?;
    let y = nr.into_solution();
    let xi: Vec<f64> = (0..n).map(|k| y[3 * k]).collect();
    let dp: Vec<f64> = (0..n).map(|k| y[3 * k + 1]).collect();
    let omega2 = y[2];
    if !(omega2 > 0.0) || xi.iter().chain(&dp).any(|v| !v.is_finite()) {
        return Err("no physical solution".into());
    }
    Ok((xi, dp, omega2))
}

fn count_nodes(xi: &[f64]) -> usize {
    xi.windows(2).filter(|w| w[0].signum() != w[1].signum() && w[0] != 0.0).count()
}

/// Fundamental radial mode of a model. A previous mode (from the preceding
/// model of an evolution run, on the same mesh) is used as the starting guess
/// when available; otherwise a few homologous guesses are tried.
pub fn fundamental_mode(model: &StarModel, previous: Option<&RadialMode>) -> Result<RadialMode, String> {
    let pr = &model.profile;
    let n = pr.len();
    let mut guesses: Vec<(Vec<f64>, Vec<f64>, f64)> = Vec::new();
    if let Some(prev) = previous.filter(|m| m.xi.len() == n && m.dp.len() == n) {
        guesses.push((prev.xi.clone(), prev.dp.clone(), prev.omega2));
    }
    for &w in &[3.0, 8.0, 15.0, 1.5, 30.0] {
        let xi = vec![1.0; n];
        let dp: Vec<f64> = (0..n).map(|k| -3.0 * pr.gamma1[k] * (1.0 - pr.q[k]) - (4.0 + w) * pr.q[k]).collect();
        guesses.push((xi, dp, w));
    }
    let mut last = String::from("no guess converged");
    for (gx, gp, gw) in guesses {
        match solve_mode(&pr.m, &pr.r, &pr.p, &pr.rho, &pr.gamma1, &gx, &gp, gw) {
            Ok((xi, dp, omega2)) => {
                let nodes = count_nodes(&xi);
                if nodes != 0 {
                    last = format!("found a mode with {nodes} nodes instead of the fundamental");
                    continue;
                }
                let mass = model.params.mass * M_SUN;
                let radius = pr.r[n - 1];
                let omega = (omega2 * G * mass / radius.powi(3)).sqrt();
                let period_s = std::f64::consts::TAU / omega;
                let rho_mean = mass / (FOUR_PI / 3.0 * radius.powi(3));
                let rho_sun = M_SUN / (FOUR_PI / 3.0 * R_SUN.powi(3));
                let q_days = period_s / 86400.0 * (rho_mean / rho_sun).sqrt();
                return Ok(RadialMode { period_s, omega2, nodes, xi, dp, q_days });
            }
            Err(e) => last = e,
        }
    }
    Err(last)
}

/// Classical pulsators that can be identified from the HR position.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum PulsatorType {
    /// δ Scuti star: on or near the main sequence in the instability strip.
    DeltaScuti,
    /// Cepheid-like: luminous star crossing the instability strip.
    Cepheid,
}

impl PulsatorType {
    pub fn label(self) -> &'static str {
        match self {
            PulsatorType::DeltaScuti => "δ Scuti-type pulsator",
            PulsatorType::Cepheid => "Cepheid-type pulsator",
        }
    }
}

/// Lowest luminosity [log L/L☉] at which the strip is applied.
pub const STRIP_MIN_LOG_L: f64 = 0.3;
/// Luminosity separating δ Scuti from Cepheid-like pulsators [log L/L☉].
pub const STRIP_CEPHEID_LOG_L: f64 = 1.7;

/// Approximate blue and red edges of the classical instability strip, as
/// log₁₀ T_eff at a given log₁₀ L/L☉. Linear fits to the observed strip:
/// δ Scuti stars (6300–8700 K at log L ≈ 0.5) up to classical Cepheids
/// (≈ 5400–6800 K at log L ≈ 3.5).
pub fn strip_edges(log_l: f64) -> (f64, f64) {
    let blue = 3.94 - 0.037 * (log_l - 0.5);
    let red = 3.80 - 0.023 * (log_l - 0.5);
    (blue, red)
}

/// Pulsator type if (T_eff, L) lies inside the instability strip.
pub fn classify(teff: f64, lum_lsun: f64) -> Option<PulsatorType> {
    let log_l = lum_lsun.log10();
    if log_l < STRIP_MIN_LOG_L {
        return None;
    }
    let (blue, red) = strip_edges(log_l);
    let lt = teff.log10();
    if lt <= blue && lt >= red {
        Some(if log_l < STRIP_CEPHEID_LOG_L { PulsatorType::DeltaScuti } else { PulsatorType::Cepheid })
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A uniform-density sphere with constant Γ₁ has the homologous
    /// fundamental mode ξ = const, with Ω² = 3Γ₁ − 4.
    #[test]
    fn homogeneous_sphere_fundamental_mode() {
        let n = 400;
        let (mass, radius) = (M_SUN, R_SUN);
        let rho = mass / (FOUR_PI / 3.0 * radius.powi(3));
        let q: Vec<f64> = (0..n).map(|i| 1e-6 + (1.0 - 1e-6) * (i as f64 / (n - 1) as f64).powf(1.5)).collect();
        let m: Vec<f64> = q.iter().map(|q| q * mass).collect();
        let r: Vec<f64> = q.iter().map(|q| radius * q.cbrt()).collect();
        let p: Vec<f64> = r
            .iter()
            .map(|r| (2.0 * std::f64::consts::PI / 3.0) * G * rho * rho * (radius * radius - r * r) + 1.0)
            .collect();
        let rho_v = vec![rho; n];
        for &g1 in &[5.0 / 3.0, 1.5] {
            let gam = vec![g1; n];
            let dp0: Vec<f64> = vec![-3.0 * g1; n];
            let (xi, _, w) = solve_mode(&m, &r, &p, &rho_v, &gam, &vec![1.0; n], &dp0, 2.0).unwrap();
            assert!((w - (3.0 * g1 - 4.0)).abs() < 2e-3, "Γ₁ = {g1}: Ω² = {w}");
            assert_eq!(count_nodes(&xi), 0);
            assert!(xi.iter().all(|x| (x - 1.0).abs() < 1e-2));
        }
    }

    #[test]
    fn instability_strip() {
        // A 1.8 M☉ main-sequence star at ~7500 K, 10 L☉: δ Scuti.
        assert_eq!(classify(7500.0, 10.0), Some(PulsatorType::DeltaScuti));
        // A 6000 K star at 3000 L☉: Cepheid.
        assert_eq!(classify(6000.0, 3000.0), Some(PulsatorType::Cepheid));
        // The Sun and hot B stars are outside.
        assert_eq!(classify(5772.0, 1.0), None);
        assert_eq!(classify(20000.0, 1000.0), None);
    }
}
