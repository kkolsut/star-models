//! Microphysics: equation of state, opacity, nuclear energy generation,
//! convection (mixing-length theory) and the photospheric boundary condition.
//!
//! EOS, opacity and nuclear rates sit behind traits, so tabulated physics
//! (OPAL/OP opacities, Saha EOS, …) can be swapped in without touching the solver.

pub mod atmosphere;
pub mod convection;
pub mod eos;
pub mod eos_fd;
pub mod fermi_dirac;
pub mod interp;
pub mod mass_loss;
pub mod nuclear;
pub mod opacity;
pub mod opacity_table;
pub mod options;

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
    /// Specific internal (thermal) energy u [erg g^-1]: no rest mass, no
    /// nuclear binding energy (used by the total-energy check of evolution runs).
    /// With ionisation it includes the ionisation energy (zero: neutral atoms).
    pub u: f64,
    /// Electron degeneracy parameter η = μ_e^chem/kT (kinetic chemical
    /// potential). NaN from EOSs without Fermi–Dirac electrons; the profile
    /// output then uses the classical value ln(n_e λ_e³/2).
    pub eta: f64,
    /// Ionisation fractions: H⁺/H, He⁺/He, He²⁺/He (1, 0, 1 when fully ionised).
    pub x_h: f64,
    pub x_he1: f64,
    pub x_he2: f64,
    /// Fraction of the hydrogen nuclei bound in H₂.
    pub x_h2: f64,
}

impl EosState {
    /// Thermodynamic state from ρ, h and the first derivatives at fixed
    /// composition (`docs/physics/02`, plan 6.1):
    /// α = (∂ln ρ/∂ln P)_T, δ = −(∂ln ρ/∂ln T)_P, c_P = (∂h/∂T)_P.
    /// The remaining fields (η, ionisation) are set to the fully ionised values.
    #[allow(clippy::too_many_arguments)]
    pub fn from_derivatives(
        p: f64,
        t: f64,
        rho: f64,
        h: f64,
        alpha: f64,
        delta: f64,
        cp: f64,
        p_rad: f64,
        mu: f64,
    ) -> Self {
        let nabla_ad = p * delta / (t * rho * cp);
        let gamma1 = 1.0 / (alpha - delta * nabla_ad);
        let beta = 1.0 - p_rad / p;
        EosState {
            rho,
            p_gas: p - p_rad,
            p_rad,
            beta,
            mu,
            nabla_ad,
            cp,
            delta,
            gamma1,
            u: h - p / rho,
            eta: f64::NAN,
            x_h: 1.0,
            x_he1: 0.0,
            x_he2: 1.0,
            x_h2: 0.0,
        }
    }

    /// A state of NaNs: what an EOS returns outside its domain (the solver then
    /// rejects the trial state; nothing is clamped).
    pub fn invalid() -> Self {
        let n = f64::NAN;
        EosState {
            rho: n,
            p_gas: n,
            p_rad: n,
            beta: n,
            mu: n,
            nabla_ad: n,
            cp: n,
            delta: n,
            gamma1: n,
            u: n,
            eta: n,
            x_h: n,
            x_he1: n,
            x_he2: n,
            x_h2: n,
        }
    }

    /// Specific enthalpy h = u + P/ρ.
    pub fn h(&self) -> f64 {
        self.u + (self.p_gas + self.p_rad) / self.rho
    }
}

/// Equation of state: density and thermodynamic derivatives as functions of (P, T).
pub trait Eos: Send + Sync {
    fn name(&self) -> &'static str;
    /// State at (P, T); outside the EOS's domain every field is NaN
    /// ([`EosState::invalid`]).
    fn eval(&self, p: f64, t: f64, comp: &Composition) -> EosState;
    /// Whether the EOS treats electron degeneracy (Fermi–Dirac electrons).
    fn degenerate_electrons(&self) -> bool {
        false
    }
    /// Like [`Eos::eval`], with the reason when (P, T) lies outside the domain.
    fn try_eval(&self, p: f64, t: f64, comp: &Composition) -> Result<EosState, String> {
        let s = self.eval(p, t, comp);
        if s.rho.is_finite() { Ok(s) } else { Err(format!("EOS undefined at P = {p:e}, T = {t:e}")) }
    }
}

/// Rosseland mean opacity κ(ρ, T) [cm^2 g^-1].
pub trait Opacity: Send + Sync {
    fn name(&self) -> &'static str;
    fn kappa(&self, rho: f64, t: f64, comp: &Composition) -> f64;
}

/// Nuclear energy generation rates [erg g^-1 s^-1], split by process.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NuclearRates {
    /// pp chain, including the branching factor ψ.
    pub pp: f64,
    pub cno: f64,
    pub triple_alpha: f64,
    /// pp branching: helium nuclei made per pp I chain, φ ∈ [1, 2] (1 = pp I only).
    /// Hydrogen consumption follows from it (see [`nuclear::q_pp_eff`]).
    pub pp_phi: f64,
}

impl Default for NuclearRates {
    fn default() -> Self {
        Self { pp: 0.0, cno: 0.0, triple_alpha: 0.0, pp_phi: 1.0 }
    }
}

impl NuclearRates {
    pub fn total(&self) -> f64 {
        self.pp + self.cno + self.triple_alpha
    }
}

/// Exponents of the screening factors f = exp(exponent), per channel, before
/// any cap is applied. For a channel with several reactions the largest is given.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ScreeningExponents {
    pub pp: f64,
    pub cno: f64,
    pub triple_alpha: f64,
}

/// Nuclear energy generation ε(ρ, T).
pub trait Nuclear: Send + Sync {
    fn name(&self) -> &'static str;
    fn eps(&self, rho: f64, t: f64, comp: &Composition) -> NuclearRates;
    /// Screening exponents at (ρ, T), for the screening diagnostic; `None` if
    /// the rates are unscreened or the implementation does not report them.
    fn screening(&self, _rho: f64, _t: f64, _comp: &Composition) -> Option<ScreeningExponents> {
        None
    }
    /// Value at which the implementation caps a screening exponent, if any.
    fn screening_cap(&self) -> Option<f64> {
        None
    }
}

/// Bundle of microphysics used by the structure solver.
pub struct Physics {
    pub eos: Box<dyn Eos>,
    pub opacity: Box<dyn Opacity>,
    pub nuclear: Box<dyn Nuclear>,
    /// Outer boundary condition.
    pub atmosphere: atmosphere::Atmosphere,
}

impl Physics {
    /// Self-contained analytic physics: ideal gas + radiation EOS, Kramers +
    /// H⁻ + electron-scattering opacity, analytic pp / CNO / 3α rates, grey
    /// Eddington photosphere.
    pub fn analytic() -> Self {
        Self {
            eos: Box::new(eos::IdealGasRadiation),
            opacity: Box::new(opacity::AnalyticOpacity),
            nuclear: Box::new(nuclear::AnalyticNuclear),
            atmosphere: atmosphere::Atmosphere::GreyClosedForm,
        }
    }

    /// Physics assembled from the selectable options.
    pub fn from_options(o: &options::PhysicsOptions) -> Self {
        Self {
            eos: o.eos.build(),
            opacity: o.opacity.build(),
            nuclear: Box::new(nuclear::AnalyticNuclear),
            atmosphere: o.atmosphere.build(),
        }
    }

    /// Short human-readable description of the physics in use.
    pub fn describe(&self) -> String {
        format!(
            "EOS: {}; opacity: {}; nuclear: {}; atmosphere: {}",
            self.eos.name(),
            self.opacity.name(),
            self.nuclear.name(),
            self.atmosphere.describe()
        )
    }
}

impl Default for Physics {
    fn default() -> Self {
        Self::analytic()
    }
}
