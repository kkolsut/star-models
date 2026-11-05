//! Evaluate all the local physics at one (m, r, P, T, L) state.

use crate::composition::Composition;
use crate::constants::{A_RAD, C, FOUR_PI, G};
use crate::physics::convection::{MltInput, MltResult, mlt};
use crate::physics::{EosState, NuclearRates, Physics};

#[derive(Clone, Copy, Debug)]
pub struct LocalState {
    pub m: f64,
    pub r: f64,
    pub p: f64,
    pub t: f64,
    pub l: f64,
    pub eos: EosState,
    pub kappa: f64,
    pub nuc: NuclearRates,
    /// Local gravity G m / r².
    pub g: f64,
    /// Pressure scale height P / (ρ g).
    pub hp: f64,
    /// Radiative gradient ∇_rad = 3 κ L P / (16 π a c G m T⁴).
    pub nabla_rad: f64,
    pub mlt: MltResult,
}

impl LocalState {
    pub fn rho(&self) -> f64 {
        self.eos.rho
    }
    pub fn eps(&self) -> f64 {
        self.nuc.total()
    }
    /// Actual temperature gradient d ln T / d ln P.
    pub fn nabla(&self) -> f64 {
        self.mlt.nabla
    }
}

/// Evaluate EOS, opacity, nuclear rates and convection at one point.
#[allow(clippy::too_many_arguments)]
pub fn evaluate(
    phys: &Physics,
    alpha_mlt: f64,
    comp: &Composition,
    m: f64,
    r: f64,
    p: f64,
    t: f64,
    l: f64,
) -> LocalState {
    let eos = phys.eos.eval(p, t, comp);
    let kappa = phys.opacity.kappa(eos.rho, t, comp);
    let nuc = phys.nuclear.eps(eos.rho, t, comp);
    let g = G * m / (r * r);
    let hp = p / (eos.rho * g);
    let nabla_rad = 3.0 * kappa * l * p / (4.0 * FOUR_PI * A_RAD * C * G * m * t.powi(4));
    let mlt = mlt(&MltInput {
        nabla_rad,
        nabla_ad: eos.nabla_ad,
        t,
        rho: eos.rho,
        p,
        kappa,
        cp: eos.cp,
        delta: eos.delta,
        g,
        alpha: alpha_mlt,
    });
    LocalState { m, r, p, t, l, eos, kappa, nuc, g, hp, nabla_rad, mlt }
}
