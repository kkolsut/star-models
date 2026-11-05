//! Numerical machinery: a banded linear solver and a generic Henyey relaxation
//! (Newton–Raphson) solver for two-point boundary-value problems.

pub mod banded;
pub mod henyey;

pub use henyey::{BvpProblem, HenyeyNewton, IterInfo, NewtonOptions, StepOutcome};
