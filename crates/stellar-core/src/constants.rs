//! Physical and astronomical constants in CGS units.

/// Gravitational constant [cm^3 g^-1 s^-2].
pub const G: f64 = 6.674_30e-8;
/// Speed of light [cm s^-1].
pub const C: f64 = 2.997_924_58e10;
/// Stefan–Boltzmann constant [erg cm^-2 s^-1 K^-4].
pub const SIGMA_SB: f64 = 5.670_374_419e-5;
/// Radiation density constant a = 4σ/c [erg cm^-3 K^-4].
pub const A_RAD: f64 = 4.0 * SIGMA_SB / C;
/// Boltzmann constant [erg K^-1].
pub const K_B: f64 = 1.380_649e-16;
/// Atomic mass unit [g].
pub const M_U: f64 = 1.660_539_066_60e-24;
/// Gas constant per atomic mass unit, k_B / m_u [erg g^-1 K^-1].
pub const R_GAS: f64 = K_B / M_U;

/// Solar mass [g].
pub const M_SUN: f64 = 1.988_47e33;
/// Nominal solar radius [cm] (IAU 2015).
pub const R_SUN: f64 = 6.957e10;
/// Nominal solar luminosity [erg s^-1] (IAU 2015).
pub const L_SUN: f64 = 3.828e33;
/// Nominal solar effective temperature [K] (IAU 2015).
pub const TEFF_SUN: f64 = 5772.0;

/// 4π, used everywhere in the structure equations.
pub const FOUR_PI: f64 = 4.0 * std::f64::consts::PI;

/// Julian year [s].
pub const YEAR: f64 = 3.155_76e7;
