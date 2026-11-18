//! Selectable physics (plans 3.x, 6.x, 7.x). The defaults reproduce the
//! analytic physics of [`super::Physics::analytic`]; the joint assessment
//! (plan 7.7) records what each choice changes.

use serde::{Deserialize, Serialize};

use super::atmosphere::{Atmosphere, IntegratedAtmosphere, TTau};
use super::{Eos, Opacity, eos, opacity, opacity_table};

/// Equation of state.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum EosChoice {
    /// Fully ionised ideal gas + radiation.
    #[default]
    Ideal,
    /// Fermi–Dirac electrons, ideal ions, radiation; full ionisation (plan 6.2).
    Degenerate,
    /// As `Degenerate` with H and He ionisation equilibrium (plan 6.3). Not
    /// complete without pressure ionisation: recombines in dense interiors.
    Ionisation,
    /// Ionisation with the Pols-type pressure-ionisation term (plan 6.4) and
    /// molecular hydrogen (plan 6.5).
    Full,
    /// Paczyński (1983) interpolation, for comparison only (plan 6.5).
    Paczynski,
}

impl EosChoice {
    pub const ALL: [EosChoice; 5] =
        [EosChoice::Ideal, EosChoice::Degenerate, EosChoice::Ionisation, EosChoice::Full, EosChoice::Paczynski];

    pub fn build(self) -> Box<dyn Eos> {
        match self {
            EosChoice::Ideal => Box::new(eos::IdealGasRadiation),
            EosChoice::Degenerate => Box::new(eos::FermiDiracEos::fully_ionised()),
            EosChoice::Ionisation => Box::new(eos::FermiDiracEos::ionising(false, false)),
            EosChoice::Full => Box::new(eos::FermiDiracEos::ionising(true, true)),
            EosChoice::Paczynski => Box::new(eos::Paczynski),
        }
    }

    /// Command-line name.
    pub fn key(self) -> &'static str {
        match self {
            EosChoice::Ideal => "ideal",
            EosChoice::Degenerate => "degenerate",
            EosChoice::Ionisation => "ionisation",
            EosChoice::Full => "full",
            EosChoice::Paczynski => "paczynski",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            EosChoice::Ideal => "ideal gas + radiation (fully ionised)",
            EosChoice::Degenerate => "Fermi–Dirac electrons (fully ionised)",
            EosChoice::Ionisation => "Fermi–Dirac + H/He ionisation (no pressure ionisation)",
            EosChoice::Full => "Fermi–Dirac + ionisation + pressure ionisation + H₂",
            EosChoice::Paczynski => "Paczyński interpolation (comparison only)",
        }
    }

    pub fn from_key(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|c| c.key() == s)
    }
}

/// Opacity.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum OpacityChoice {
    /// Kramers + H⁻ + electron scattering formulas.
    #[default]
    Analytic,
    /// OPAL GN93 + Ferguson et al. (2005) tables with Potekhin conduction.
    Tables,
}

impl OpacityChoice {
    pub const ALL: [OpacityChoice; 2] = [OpacityChoice::Analytic, OpacityChoice::Tables];

    pub fn build(self) -> Box<dyn Opacity> {
        match self {
            OpacityChoice::Analytic => Box::new(opacity::AnalyticOpacity),
            OpacityChoice::Tables => Box::new(opacity_table::TabulatedOpacity::default()),
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            OpacityChoice::Analytic => "analytic",
            OpacityChoice::Tables => "tables",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            OpacityChoice::Analytic => "analytic (Kramers + H⁻ + e⁻ scattering)",
            OpacityChoice::Tables => "OPAL + Ferguson tables, conduction (GN93 mixture)",
        }
    }

    pub fn from_key(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|c| c.key() == s)
    }
}

/// Atmosphere (outer boundary condition).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum AtmosphereChoice {
    /// Grey Eddington photosphere with one κ.
    #[default]
    Grey,
    /// Integrated grey atmosphere with convection (plans 3.1, 3.2).
    IntegratedGrey,
    /// Integrated atmosphere with the Krishna Swamy T(τ) and convection (plan 3.3).
    KrishnaSwamy,
}

impl AtmosphereChoice {
    pub const ALL: [AtmosphereChoice; 3] =
        [AtmosphereChoice::Grey, AtmosphereChoice::IntegratedGrey, AtmosphereChoice::KrishnaSwamy];

    pub fn build(self) -> Atmosphere {
        match self {
            AtmosphereChoice::Grey => Atmosphere::GreyClosedForm,
            AtmosphereChoice::IntegratedGrey => {
                Atmosphere::Integrated(IntegratedAtmosphere::new(TTau::Eddington, true))
            }
            AtmosphereChoice::KrishnaSwamy => {
                Atmosphere::Integrated(IntegratedAtmosphere::new(TTau::KrishnaSwamy, true))
            }
        }
    }

    pub fn key(self) -> &'static str {
        match self {
            AtmosphereChoice::Grey => "grey",
            AtmosphereChoice::IntegratedGrey => "integrated",
            AtmosphereChoice::KrishnaSwamy => "krishna-swamy",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            AtmosphereChoice::Grey => "grey Eddington, one κ (τ = 2/3)",
            AtmosphereChoice::IntegratedGrey => "integrated grey, with convection",
            AtmosphereChoice::KrishnaSwamy => "integrated Krishna Swamy T(τ), with convection",
        }
    }

    pub fn from_key(s: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|c| c.key() == s)
    }
}

/// The selectable parts of the physics.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct PhysicsOptions {
    #[serde(default)]
    pub eos: EosChoice,
    #[serde(default)]
    pub opacity: OpacityChoice,
    #[serde(default)]
    pub atmosphere: AtmosphereChoice,
}

impl PhysicsOptions {
    /// The analytic physics (the default).
    pub fn analytic() -> Self {
        Self::default()
    }

    /// Every upgrade of stage 4 on: degenerate, ionising EOS with pressure
    /// ionisation and H₂, opacity tables, Krishna Swamy atmosphere.
    pub fn all_upgrades() -> Self {
        Self { eos: EosChoice::Full, opacity: OpacityChoice::Tables, atmosphere: AtmosphereChoice::KrishnaSwamy }
    }
}
