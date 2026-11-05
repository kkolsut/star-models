//! The stellar structure equations as a Henyey boundary-value problem.
//!
//! Independent variable: mass m. Unknowns per point:
//! `y = (ln r, ln P, ln T, L / L_scale)`.
//!
//! * d ln r / dm = 1 / (4π r³ ρ)                    (mass continuity)
//! * d ln P / dm = −G m / (4π r⁴ P)                  (hydrostatic equilibrium)
//! * dL / dm     = ε_nuc                             (energy generation)
//! * d ln T / dm = ∇ · d ln P / dm                   (energy transport)
//!
//! Here ∇ = ∇_rad in radiative layers. In convective layers
//! (Schwarzschild: ∇_rad > ∇_ad), ∇ comes from mixing-length theory.
//!
//! Boundary conditions:
//! * inner (first point, m₀ ≪ M): r₀³ = 3 m₀ / (4π ρ₀) and L₀ = ε₀ m₀
//! * outer (photosphere): L = 4π R² σ T⁴ and P = ⅔ G M / (R² κ) + a T⁴ / 6

use super::local::{LocalState, evaluate};
use crate::composition::Composition;
use crate::constants::{FOUR_PI, G, SIGMA_SB};
use crate::physics::Physics;
use crate::physics::atmosphere::photosphere_pressure;
use crate::solver::BvpProblem;

pub const NVAR: usize = 4;
pub const I_LNR: usize = 0;
pub const I_LNP: usize = 1;
pub const I_LNT: usize = 2;
pub const I_L: usize = 3;

pub struct StructureProblem<'a> {
    pub physics: &'a Physics,
    /// Total mass [g].
    pub mass: f64,
    /// Mass coordinate of every mesh point [g].
    pub m: Vec<f64>,
    /// Composition at every mesh point (uniform for ZAMS models).
    pub comp: Vec<Composition>,
    pub alpha_mlt: f64,
    /// Luminosity scale [erg/s], used to keep the L unknown of order unity.
    pub l_scale: f64,
}

impl<'a> StructureProblem<'a> {
    pub fn new(physics: &'a Physics, mass: f64, q: &[f64], comp: Composition, alpha_mlt: f64, l_scale: f64) -> Self {
        Self { physics, mass, m: q.iter().map(|&q| q * mass).collect(), comp: vec![comp; q.len()], alpha_mlt, l_scale }
    }

    /// Evaluate the full local physics at mesh point `k` for unknowns `y`.
    pub fn local_at(&self, k: usize, y: &[f64]) -> LocalState {
        evaluate(
            self.physics,
            self.alpha_mlt,
            &self.comp[k],
            self.m[k],
            y[I_LNR].exp(),
            y[I_LNP].exp(),
            y[I_LNT].exp(),
            y[I_L] * self.l_scale,
        )
    }
}

impl BvpProblem for StructureProblem<'_> {
    fn nvar(&self) -> usize {
        NVAR
    }

    fn n_inner(&self) -> usize {
        2
    }

    fn npts(&self) -> usize {
        self.m.len()
    }

    fn inner_bc(&self, y0: &[f64], out: &mut [f64]) {
        let comp = &self.comp[0];
        let (p, t) = (y0[I_LNP].exp(), y0[I_LNT].exp());
        let rho = self.physics.eos.eval(p, t, comp).rho;
        let m0 = self.m[0];
        out[0] = y0[I_LNR] - (3.0 * m0 / (FOUR_PI * rho)).ln() / 3.0;
        let eps = self.physics.nuclear.eps(rho, t, comp).total();
        out[1] = y0[I_L] - eps * m0 / self.l_scale;
    }

    fn interval(&self, k: usize, a: &[f64], b: &[f64], out: &mut [f64]) {
        let dm = self.m[k + 1] - self.m[k];
        let m = 0.5 * (self.m[k] + self.m[k + 1]);
        let comp = self.comp[k].lerp(&self.comp[k + 1], 0.5);
        let mid = |i: usize| 0.5 * (a[i] + b[i]);
        let s = evaluate(
            self.physics,
            self.alpha_mlt,
            &comp,
            m,
            mid(I_LNR).exp(),
            mid(I_LNP).exp(),
            mid(I_LNT).exp(),
            mid(I_L) * self.l_scale,
        );
        let r = s.r;
        let dlnp_dm = -G * m / (FOUR_PI * r.powi(4) * s.p);
        out[0] = (b[I_LNR] - a[I_LNR]) - dm / (FOUR_PI * r.powi(3) * s.rho());
        out[1] = (b[I_LNP] - a[I_LNP]) - dm * dlnp_dm;
        out[2] = (b[I_L] - a[I_L]) - dm * s.eps() / self.l_scale;
        out[3] = (b[I_LNT] - a[I_LNT]) - dm * s.nabla() * dlnp_dm;
    }

    fn outer_bc(&self, yn: &[f64], out: &mut [f64]) {
        let k = self.m.len() - 1;
        let comp = &self.comp[k];
        let (r, p, t) = (yn[I_LNR].exp(), yn[I_LNP].exp(), yn[I_LNT].exp());
        let rho = self.physics.eos.eval(p, t, comp).rho;
        let kappa = self.physics.opacity.kappa(rho, t, comp);
        out[0] = yn[I_LNP] - photosphere_pressure(self.mass, r, t, kappa).ln();
        let l = yn[I_L] * self.l_scale;
        out[1] = l / (FOUR_PI * SIGMA_SB * r * r * t.powi(4)) - 1.0;
    }

    fn max_correction(&self, var: usize) -> f64 {
        match var {
            I_LNR => 0.2,
            I_LNP => 0.5,
            I_LNT => 0.2,
            _ => 0.5,
        }
    }

    fn fd_step(&self, var: usize, value: f64) -> f64 {
        match var {
            I_L => 1e-6 * value.abs().max(1e-6),
            _ => 1e-6,
        }
    }
}
