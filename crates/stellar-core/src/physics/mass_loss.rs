//! Stellar winds: mass-loss rates (plan 4.1, `docs/physics/07-evolution.md`, "Mass loss").
//!
//! * Hot stars (T_eff ≥ 8000 K): the linear form of the de Jager, Nieuwenhuijzen
//!   & van der Hucht (1988) fit, §4 remark 3, p. 264,
//!   log(−Ṁ) = −8.158 + 1.769 log(L/L☉) − 1.676 log T_eff  [M☉ yr⁻¹],
//!   times the metallicity factor (Z/0.019)^0.69 of Vink, de Koter & Lamers
//!   (2001), §7 eqs. 21–22 (O-star side). Combining the two is a choice made
//!   here, not a prescription of either paper.
//! * Cool stars (T_eff < 8000 K): Reimers (1975) as confirmed by Reimers (1977,
//!   p. 217) and McDonald & Zijlstra (2015, eq. 1),
//!   −Ṁ = 4 × 10⁻¹³ η L R / M  [M☉ yr⁻¹, solar units], default η = 0.477
//!   (McDonald & Zijlstra 2015, median over 56 globular clusters). No
//!   metallicity factor (they find η varies by less than about 30% over a
//!   factor 200 in iron abundance).
//!
//! The switch between the two is a step at 8000 K.
//!
//! All functions return the *loss rate* −Ṁ ≥ 0 in M☉ yr⁻¹; [`Wind::mdot`] returns
//! the signed Ṁ ≤ 0.

use serde::{Deserialize, Serialize};

/// T_eff [K] at and above which the hot-star (de Jager) rate is used.
pub const T_SWITCH: f64 = 8000.0;
/// Default Reimers η (McDonald & Zijlstra 2015).
pub const REIMERS_ETA_DEFAULT: f64 = 0.477;
/// Reference metallicity of the Vink et al. (2001) factor (their Z☉).
pub const Z_REF: f64 = 0.019;
/// Metallicity exponent of the Vink et al. (2001) factor (T_eff ≳ 25 000 K).
pub const Z_EXPONENT: f64 = 0.69;

/// de Jager et al. (1988), linear form: −Ṁ [M☉ yr⁻¹] for L [L☉] and T_eff [K].
pub fn de_jager(l_lsun: f64, teff: f64) -> f64 {
    10f64.powf(-8.158 + 1.769 * l_lsun.log10() - 1.676 * teff.log10())
}

/// Metallicity factor (Z / 0.019)^0.69 of Vink et al. (2001).
pub fn metallicity_factor(z: f64) -> f64 {
    (z.max(0.0) / Z_REF).powf(Z_EXPONENT)
}

/// Reimers: −Ṁ = 4 × 10⁻¹³ η L R / M [M☉ yr⁻¹] for L [L☉], R [R☉], M [M☉].
pub fn reimers(l_lsun: f64, r_rsun: f64, m_msun: f64, eta: f64) -> f64 {
    4e-13 * eta * l_lsun * r_rsun / m_msun
}

/// Which formula applies at a given T_eff.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WindRegime {
    /// de Jager et al. (1988) with the Vink et al. (2001) metallicity factor.
    Hot,
    /// Reimers (1975).
    Cool,
}

impl WindRegime {
    pub fn of(teff: f64) -> Self {
        if teff >= T_SWITCH { WindRegime::Hot } else { WindRegime::Cool }
    }

    pub fn label(self) -> &'static str {
        match self {
            WindRegime::Hot => "de Jager et al. (1988) × (Z/0.019)^0.69",
            WindRegime::Cool => "Reimers (1975)",
        }
    }
}

/// Combined rate: −Ṁ [M☉ yr⁻¹] for a star with L [L☉], R [R☉], M [M☉],
/// T_eff [K] and surface metallicity Z, with Reimers parameter η.
pub fn loss_rate(l_lsun: f64, r_rsun: f64, m_msun: f64, teff: f64, z: f64, eta: f64) -> f64 {
    match WindRegime::of(teff) {
        WindRegime::Hot => de_jager(l_lsun, teff) * metallicity_factor(z),
        WindRegime::Cool => reimers(l_lsun, r_rsun, m_msun, eta),
    }
}

/// Mass loss of an evolution run.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Wind {
    /// Constant mass (the default).
    #[default]
    Off,
    /// de Jager (hot) / Reimers (cool) rates, multiplied by `efficiency`.
    Standard {
        /// Factor applied to both rates (1 = the published rates).
        efficiency: f64,
        /// Reimers η.
        reimers_eta: f64,
    },
    /// A fixed loss rate −Ṁ [M☉ yr⁻¹], for tests and experiments.
    Constant { loss_msun_yr: f64 },
}

impl Wind {
    /// The published rates (efficiency 1, η = 0.477).
    pub fn standard() -> Self {
        Wind::Standard { efficiency: 1.0, reimers_eta: REIMERS_ETA_DEFAULT }
    }

    pub fn is_off(&self) -> bool {
        match *self {
            Wind::Off => true,
            Wind::Standard { efficiency, .. } => !(efficiency > 0.0),
            Wind::Constant { loss_msun_yr } => !(loss_msun_yr > 0.0),
        }
    }

    /// Signed rate Ṁ ≤ 0 [M☉ yr⁻¹] for a star with L [L☉], R [R☉], M [M☉],
    /// T_eff [K] and surface metallicity Z.
    pub fn mdot(&self, l_lsun: f64, r_rsun: f64, m_msun: f64, teff: f64, z: f64) -> f64 {
        let loss = match *self {
            Wind::Off => 0.0,
            Wind::Standard { efficiency, reimers_eta } => {
                efficiency.max(0.0) * loss_rate(l_lsun, r_rsun, m_msun, teff, z, reimers_eta)
            }
            Wind::Constant { loss_msun_yr } => loss_msun_yr.max(0.0),
        };
        if loss.is_finite() { -loss } else { 0.0 }
    }

    /// Check the parameters.
    pub fn validate(&self) -> Result<(), String> {
        match *self {
            Wind::Off => Ok(()),
            Wind::Standard { efficiency, reimers_eta } => {
                if !(0.0..=10.0).contains(&efficiency) {
                    return Err("wind efficiency must be in [0, 10]".into());
                }
                if !(0.0..=5.0).contains(&reimers_eta) {
                    return Err("Reimers η must be in [0, 5]".into());
                }
                Ok(())
            }
            Wind::Constant { loss_msun_yr } => {
                if !(loss_msun_yr >= 0.0 && loss_msun_yr.is_finite()) {
                    return Err("constant mass-loss rate must be ≥ 0".into());
                }
                Ok(())
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn de_jager_central_point_of_the_paper() {
        // de Jager et al. (1988), p. 264: log(−Ṁ) = −6.81 at log T_eff = 4.05, log L = 4.6.
        let r = de_jager(10f64.powf(4.6), 10f64.powf(4.05));
        assert!((r.log10() + 6.81).abs() < 0.005, "{}", r.log10());
    }

    #[test]
    fn metallicity_factor_is_one_at_reference() {
        assert!((metallicity_factor(Z_REF) - 1.0).abs() < 1e-15);
        assert!((metallicity_factor(0.0019) - 0.1f64.powf(0.69)).abs() < 1e-12);
    }

    #[test]
    fn reimers_solar_units() {
        // A 1 M☉ giant with L = 1000 L☉, R = 100 R☉, η = 0.477.
        let r = reimers(1000.0, 100.0, 1.0, 0.477);
        assert!((r - 4e-13 * 0.477 * 1e5).abs() < 1e-20);
    }

    #[test]
    fn switch_at_8000_k() {
        assert_eq!(WindRegime::of(7999.0), WindRegime::Cool);
        assert_eq!(WindRegime::of(8000.0), WindRegime::Hot);
        let hot = loss_rate(1e4, 10.0, 10.0, 8000.0, 0.019, 0.5);
        assert!((hot - de_jager(1e4, 8000.0)).abs() < 1e-20);
        let cool = loss_rate(1e4, 10.0, 10.0, 7999.0, 0.019, 0.5);
        assert!((cool - reimers(1e4, 10.0, 10.0, 0.5)).abs() < 1e-24);
    }

    #[test]
    fn wind_signs_and_off() {
        assert_eq!(Wind::Off.mdot(1e5, 10.0, 30.0, 3e4, 0.02), 0.0);
        assert!(Wind::standard().mdot(1e5, 10.0, 30.0, 3e4, 0.02) < 0.0);
        let half = Wind::Standard { efficiency: 0.5, reimers_eta: 0.477 }.mdot(1e5, 10.0, 30.0, 3e4, 0.02);
        assert!((half / Wind::standard().mdot(1e5, 10.0, 30.0, 3e4, 0.02) - 0.5).abs() < 1e-14);
        assert_eq!(Wind::Constant { loss_msun_yr: 1e-7 }.mdot(1.0, 1.0, 1.0, 5800.0, 0.02), -1e-7);
        assert!(Wind::Standard { efficiency: 0.0, reimers_eta: 0.477 }.is_off());
    }
}
