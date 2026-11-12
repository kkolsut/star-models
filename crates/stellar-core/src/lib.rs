//! # stellar-core
//!
//! Computes static, chemically homogeneous zero-age main-sequence (ZAMS) stellar
//! models. It solves the four stellar structure equations in Lagrangian (mass)
//! form with a Henyey-type Newton–Raphson relaxation solver.
//!
//! The crate has no GUI dependencies, so it can be used from the CLI, the egui
//! desktop app, or (later) a WASM build.
//!
//! Units are CGS internally. Results carry solar-normalised helpers.

// `!(a > b)` is used on purpose throughout the numerics: it also takes the safe
// branch when a value is NaN, which `a <= b` would not.
#![allow(clippy::neg_cmp_op_on_partial_ord)]

pub mod composition;
pub mod constants;
pub mod output;
pub mod physics;
pub mod pulsation;
pub mod seismic;
pub mod solver;
pub mod structure;

pub use composition::Composition;
pub use physics::Physics;
pub use structure::evolve::{Evolution, EvolutionOptions, EvolutionStep, StopCondition, fmt_age};
pub use structure::model::{ConvectiveZone, Globals, Phase, Profile, Residuals, StarModel, StarParams};
pub use structure::solve::{NoControl, SolveControl, SolveError, SolveEvent, solve_model};
