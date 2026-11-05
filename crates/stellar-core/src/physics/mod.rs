//! Microphysics: equation of state, opacity, nuclear energy generation,
//! convection (mixing-length theory) and the photospheric boundary condition.
//!
//! EOS, opacity and nuclear rates sit behind traits, so tabulated physics
//! (OPAL/OP opacities, Saha EOS, …) can be swapped in without touching the solver.

pub mod atmosphere;
pub mod convection;
pub mod eos;
pub mod nuclear;
pub mod opacity;

use crate::composition::Composition;

/// Thermodynamic state returned by an [`Eos`] for given (P, T).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct EosState {
    /// Density [g cm^-3].
    pub rho: f64,
    /// Gas pressure [dyn cm^-2].
    pub p_gas: f64,
    /// Radiation pressure [dyn cm^-2].
    pub p_rad: f64,
    /// β = P_gas / P.
    pub beta: f64,
    /// Mean molecular weight.
    pub mu: f64,
    /// Adiabatic temperature gradient ∇_ad = (∂ln T/∂ln P)_s.
    pub nabla_ad: f64,
    /// Specific heat at constant pressure [erg g^-1 K^-1].
    pub cp: f64,
    /// δ = −(∂ln ρ/∂ln T)_P.
    pub delta: f64,
    /// First adiabatic exponent Γ₁ = (∂ln P/∂ln ρ)_s (used for pulsations).
    pub gamma1: f64,
}

/// Equation of state: density and thermodynamic derivatives as functions of (P, T).
pub trait Eos: Send + Sync {
    fn name(&self) -> &'static str;
    fn eval(&self, p: f64, t: f64, comp: &Composition) -> EosState;
}

/// Rosseland mean opacity κ(ρ, T) [cm^2 g^-1].
pub trait Opacity: Send + Sync {
    fn name(&self) -> &'static str;
    fn kappa(&self, rho: f64, t: f64, comp: &Composition) -> f64;
}

/// Nuclear energy generation rates [erg g^-1 s^-1], split by process.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct NuclearRates {
    pub pp: f64,
    pub cno: f64,
    pub triple_alpha: f64,
}

impl NuclearRates {
    pub fn total(&self) -> f64 {
        self.pp + self.cno + self.triple_alpha
    }
}

/// Nuclear energy generation ε(ρ, T).
pub trait Nuclear: Send + Sync {
    fn name(&self) -> &'static str;
    fn eps(&self, rho: f64, t: f64, comp: &Composition) -> NuclearRates;
}

/// Bundle of microphysics used by the structure solver.
pub struct Physics {
    pub eos: Box<dyn Eos>,
    pub opacity: Box<dyn Opacity>,
    pub nuclear: Box<dyn Nuclear>,
}

impl Physics {
    /// Self-contained analytic physics: ideal gas + radiation EOS, Kramers +
    /// H⁻ + electron-scattering opacity, analytic pp / CNO / 3α rates.
    pub fn analytic() -> Self {
        Self {
            eos: Box::new(eos::IdealGasRadiation),
            opacity: Box::new(opacity::AnalyticOpacity),
            nuclear: Box::new(nuclear::AnalyticNuclear),
        }
    }

    /// Short human-readable description of the physics in use.
    pub fn describe(&self) -> String {
        format!(
            "EOS: {}; opacity: {}; nuclear: {}; atmosphere: grey Eddington (T = T_eff at τ = 2/3)",
            self.eos.name(),
            self.opacity.name(),
            self.nuclear.name()
        )
    }
}

impl Default for Physics {
    fn default() -> Self {
        Self::analytic()
    }
}
