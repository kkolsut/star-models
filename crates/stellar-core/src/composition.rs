//! Chemical composition (mass fractions) and derived mean molecular weights.

use serde::{Deserialize, Serialize};

/// Mass fractions of hydrogen (X), helium (Y) and metals (Z), plus the
/// combined C+N+O mass fraction used by the CNO-cycle rate.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Composition {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub x_cno: f64,
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum CompositionError {
    #[error("hydrogen mass fraction X must be in [0, 1], got {0}")]
    InvalidX(f64),
    #[error("metallicity Z must be in [0, 1), got {0}")]
    InvalidZ(f64),
    #[error("X + Z must not exceed 1 (got X + Z = {0})")]
    Overfull(f64),
    #[error("CNO fraction of metals must be in [0, 1], got {0}")]
    InvalidCno(f64),
}

impl Composition {
    /// Build a composition from X and Z; Y = 1 − X − Z. `f_cno` is the fraction
    /// of Z that is carbon, nitrogen and oxygen.
    pub fn new(x: f64, z: f64, f_cno: f64) -> Result<Self, CompositionError> {
        if !(0.0..=1.0).contains(&x) {
            return Err(CompositionError::InvalidX(x));
        }
        if !(0.0..1.0).contains(&z) {
            return Err(CompositionError::InvalidZ(z));
        }
        if x + z > 1.0 + 1e-12 {
            return Err(CompositionError::Overfull(x + z));
        }
        if !(0.0..=1.0).contains(&f_cno) {
            return Err(CompositionError::InvalidCno(f_cno));
        }
        Ok(Self { x, y: (1.0 - x - z).max(0.0), z, x_cno: f_cno * z })
    }

    /// Mean molecular weight of a fully ionized gas:
    /// 1/μ ≈ 2X + 3Y/4 + Z/2.
    pub fn mu(&self) -> f64 {
        1.0 / (2.0 * self.x + 0.75 * self.y + 0.5 * self.z)
    }

    /// Mean molecular weight per free electron, μ_e = 2 / (1 + X).
    pub fn mu_e(&self) -> f64 {
        2.0 / (1.0 + self.x)
    }

    /// Linear interpolation between two compositions.
    pub fn lerp(&self, other: &Self, t: f64) -> Self {
        let l = |a: f64, b: f64| a + (b - a) * t;
        Self { x: l(self.x, other.x), y: l(self.y, other.y), z: l(self.z, other.z), x_cno: l(self.x_cno, other.x_cno) }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn solar_mu() {
        let c = Composition::new(0.70, 0.02, 0.7).unwrap();
        assert!((c.y - 0.28).abs() < 1e-12);
        // 1/μ = 1.4 + 0.21 + 0.01 = 1.62
        assert!((c.mu() - 1.0 / 1.62).abs() < 1e-12);
        assert!((c.mu_e() - 2.0 / 1.7).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_input() {
        assert!(Composition::new(0.9, 0.2, 0.7).is_err());
        assert!(Composition::new(-0.1, 0.02, 0.7).is_err());
        assert!(Composition::new(0.7, 0.02, 1.5).is_err());
    }
}
