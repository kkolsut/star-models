//! Chemical composition (mass fractions) and derived mean molecular weights.
//!
//! The composition vector (plan 9.5, `docs/physics/04`, "Carbon and oxygen"):
//!
//! | field | meaning |
//! |---|---|
//! | `x`, `y` | hydrogen, helium |
//! | `z` | Z₀: the initial heavy elements, a fixed scaled mixture (GN93) |
//! | `x_c`, `x_o` | carbon and oxygen produced by helium burning, in addition to Z₀ |
//! | `x_cno` | catalysts of the CNO cycle, f_CNO Z₀ plus produced carbon |
//!
//! X + Y + Z₀ + X_C + X_O = 1. Produced carbon enters `x_cno`
//! ([`Composition::cno_catalysts`]); it matters for the CNO rate only where it
//! is mixed into hydrogen, since that rate is proportional to X. The mean
//! molecular weights and the screening count every species (C: A = 12, six
//! electrons; O: A = 16, eight electrons).

use serde::{Deserialize, Serialize};

use crate::mixture::GN93;

/// Mass fractions of hydrogen (X), helium (Y), the initial metals (Z₀), and
/// carbon and oxygen made by helium burning, plus the combined C+N+O mass
/// fraction used by the CNO-cycle rate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Composition {
    pub x: f64,
    pub y: f64,
    /// Initial heavy elements Z₀ (scaled GN93 mixture).
    pub z: f64,
    pub x_cno: f64,
    /// Carbon produced by helium burning (in addition to the carbon in Z₀).
    #[serde(default)]
    pub x_c: f64,
    /// Oxygen produced by ¹²C(α,γ)¹⁶O (in addition to the oxygen in Z₀).
    #[serde(default)]
    pub x_o: f64,
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
        Ok(Self { x, y: (1.0 - x - z).max(0.0), z, x_cno: f_cno * z, x_c: 0.0, x_o: 0.0 })
    }

    /// CNO catalysts: f_CNO Z₀ plus the produced carbon (`docs/physics/04`).
    pub fn cno_catalysts(f_cno: f64, z0: f64, x_c: f64) -> f64 {
        f_cno * z0 + x_c
    }

    /// All elements heavier than helium, Z₀ + X_C + X_O (what the opacity
    /// tables and formulas, which know one Z, are given).
    pub fn z_all(&self) -> f64 {
        self.z + self.x_c + self.x_o
    }

    /// Sum of all mass fractions (1 up to rounding).
    pub fn sum(&self) -> f64 {
        self.x + self.y + self.z + self.x_c + self.x_o
    }

    /// Nuclei per atomic mass unit, Σ X_i/A_i (A = 1, 4, ⟨A⟩ of GN93, 12, 16).
    pub fn ions_per_mu(&self) -> f64 {
        self.x + self.y / 4.0 + self.z / GN93.mean_a + self.x_c / 12.0 + self.x_o / 16.0
    }

    /// Electrons per atomic mass unit at full ionisation, Σ X_i Z_i/A_i.
    pub fn electrons_per_mu(&self) -> f64 {
        self.x + self.y / 2.0 + self.z * GN93.z_over_a + 0.5 * (self.x_c + self.x_o)
    }

    /// Σ X_i Z_i²/A_i per atomic mass unit (for ⟨Z²⟩ and the Debye length).
    pub fn z2_per_mu(&self) -> f64 {
        self.x + self.y + self.z * GN93.mean_z2 / GN93.mean_a + 3.0 * self.x_c + 4.0 * self.x_o
    }

    /// Mean molecular weight of a fully ionized gas:
    /// 1/μ ≈ 2X + 3Y/4 + Z/2 + 7X_C/12 + 9X_O/16 (the metals Z₀ in the
    /// usual approximation (1 + Z_i)/A_i ≈ 1/2).
    pub fn mu(&self) -> f64 {
        1.0 / (2.0 * self.x + 0.75 * self.y + 0.5 * self.z + 7.0 / 12.0 * self.x_c + 9.0 / 16.0 * self.x_o)
    }

    /// Mean molecular weight per free electron, μ_e = 2 / (1 + X) (He, C, O
    /// and the metals all give ½ electron per atomic mass unit).
    pub fn mu_e(&self) -> f64 {
        2.0 / (1.0 + self.x)
    }

    /// Linear interpolation between two compositions.
    pub fn lerp(&self, other: &Self, t: f64) -> Self {
        let l = |a: f64, b: f64| a + (b - a) * t;
        Self {
            x: l(self.x, other.x),
            y: l(self.y, other.y),
            z: l(self.z, other.z),
            x_cno: l(self.x_cno, other.x_cno),
            x_c: l(self.x_c, other.x_c),
            x_o: l(self.x_o, other.x_o),
        }
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
        assert_eq!(c.z_all(), c.z);
    }

    #[test]
    fn carbon_and_oxygen_count_as_species() {
        let c = Composition { x: 0.0, y: 0.2, z: 0.02, x_cno: 0.014, x_c: 0.5, x_o: 0.28 };
        assert!((c.sum() - 1.0).abs() < 1e-15);
        assert!((c.z_all() - 0.8).abs() < 1e-15);
        // Pure carbon: μ = 12/7, μ_e = 2; pure oxygen: μ = 16/9.
        let pure_c = Composition { x_c: 1.0, ..Default::default() };
        assert!((pure_c.mu() - 12.0 / 7.0).abs() < 1e-12);
        assert!((pure_c.mu_e() - 2.0).abs() < 1e-12);
        assert!((pure_c.electrons_per_mu() - 0.5).abs() < 1e-15);
        assert!((pure_c.ions_per_mu() - 1.0 / 12.0).abs() < 1e-15);
        assert!((pure_c.z2_per_mu() - 3.0).abs() < 1e-15);
        let pure_o = Composition { x_o: 1.0, ..Default::default() };
        assert!((pure_o.mu() - 16.0 / 9.0).abs() < 1e-12);
        // Mean molecular weight from the ion and electron counts agrees for C and O.
        for v in [pure_c, pure_o] {
            assert!((1.0 / v.mu() - v.ions_per_mu() - v.electrons_per_mu()).abs() < 1e-12);
        }
    }

    #[test]
    fn rejects_bad_input() {
        assert!(Composition::new(0.9, 0.2, 0.7).is_err());
        assert!(Composition::new(-0.1, 0.02, 0.7).is_err());
        assert!(Composition::new(0.7, 0.02, 1.5).is_err());
    }
}
