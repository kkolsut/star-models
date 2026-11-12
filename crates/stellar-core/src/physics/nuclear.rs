//! Analytic nuclear energy generation rates (Kippenhahn & Weigert 1990 forms, CGS).
//!
//! * pp chain:  ε_pp  = ψ ε_pp,I, ε_pp,I = 2.38e6 ρ X² f₁₁ g₁₁ T₆^(−2/3) exp(−33.80 T₆^(−1/3)),
//!   g₁₁ = 1 + 0.0123 T₆^(1/3) + 0.0109 T₆^(2/3) + 0.0009 T₆ (rate correction).
//!   ψ = φ (1 − 0.021 F_II) is the pp II branching factor (pp III neglected),
//!   with φ = 1 + 2/(√(1 + 2/α) + 1), F_II = 2(φ − 1)/φ and α from its definition
//!   α = (Y/4X)² ⟨σv⟩₃₄² / (⟨σv⟩₁₁ ⟨σv⟩₃₃), the thermal rates integrated from the
//!   Solar Fusion III S(0) values (Acharya et al. 2025, Table I, eq. 3; constant S).
//!   This is a hybrid: the pp I rate is the 1990 fit, α is built from Solar
//!   Fusion III (the two p + p rates differ by about 8%).
//! * CNO cycle: ε_CNO = 8.67e27 ρ X X_CNO f₁₄ C T₆^(−2/3) exp(−152.28 T₆^(−1/3))
//! * triple-α:  ε_3α  = 5.09e11 ρ² Y³ f T₈^(−3) exp(−44.027 / T₈)
//!
//! Each includes a weak-screening factor f = exp(0.188 Z₁Z₂ ζ ρ^½ T₆^(−3/2)).

use std::sync::OnceLock;

use super::{Nuclear, NuclearRates, ScreeningExponents};
use crate::composition::Composition;

/// Energy deposited per gram of hydrogen burnt through the pp I chain
/// [erg g⁻¹] (26.20 MeV per ⁴He nucleus after neutrino losses).
pub const Q_PP_I: f64 = 6.27e18;
/// Energy deposited per gram of hydrogen burnt through the CNO cycle
/// [erg g⁻¹] (25.03 MeV per ⁴He nucleus after neutrino losses).
pub const Q_CNO: f64 = 5.99e18;
/// Energy released per gram of helium burnt to carbon by the triple-α
/// reaction [erg g⁻¹] (7.275 MeV per ¹²C nucleus).
pub const Q_HE: f64 = 5.85e17;

/// Fraction of helium made through pp II, F_II = 2(φ − 1)/φ.
pub fn pp2_fraction(phi: f64) -> f64 {
    let phi = phi.max(1.0);
    2.0 * (phi - 1.0) / phi
}

/// Energy deposited per gram of hydrogen burnt by the pp chain with branching
/// factor φ: Q_pp^eff = Q_pp,I ψ/φ = Q_pp,I (1 − 0.021 F_II) [erg g⁻¹].
pub fn q_pp_eff(phi: f64) -> f64 {
    Q_PP_I * (1.0 - 0.021 * pp2_fraction(phi))
}

/// φ = 1 − α + √(α(α + 2)) in the form that is stable for large α.
pub fn pp_phi(alpha: f64) -> f64 {
    if !(alpha > 0.0) {
        return 1.0;
    }
    1.0 + 2.0 / ((1.0 + 2.0 / alpha).sqrt() + 1.0)
}

/// Thermal rates for the pp branching α from Solar Fusion III.
pub mod sf3 {
    /// Boltzmann constant [keV K⁻¹].
    const K_B_KEV: f64 = 8.617_333_262e-8;
    /// Atomic mass unit [keV].
    const M_U_KEV: f64 = 931_494.102_42;
    /// Fine-structure constant.
    const ALPHA_FS: f64 = 1.0 / 137.035_999_084;

    /// A non-resonant reaction with constant astrophysical factor S(0).
    #[derive(Clone, Copy, Debug)]
    pub struct Reaction {
        pub z1z2: f64,
        /// Reduced mass [u] (integer mass numbers).
        pub mu: f64,
        /// S(0) [MeV b] (Acharya et al. 2025, Table I).
        pub s0: f64,
    }

    /// ¹H(p,e⁺ν)²H.
    pub const PP: Reaction = Reaction { z1z2: 1.0, mu: 0.5, s0: 4.09e-25 };
    /// ³He(³He,2p)⁴He.
    pub const HE3_HE3: Reaction = Reaction { z1z2: 4.0, mu: 1.5, s0: 5.21 };
    /// ³He(⁴He,γ)⁷Be.
    pub const HE3_HE4: Reaction = Reaction { z1z2: 4.0, mu: 12.0 / 7.0, s0: 5.61e-4 };

    /// ln of ∫₀^∞ S(E) exp(−E/kT − √(E_G/E)) dE / S(0) [keV], eq. 3 of
    /// Acharya et al. (2025) without its constant prefactor, by composite
    /// Simpson quadrature around the Gamow peak. E_G = 2 μ c² (π α Z₁Z₂)².
    pub fn ln_gamow_integral(r: &Reaction, t: f64) -> f64 {
        let kt = K_B_KEV * t;
        let b = std::f64::consts::PI * ALPHA_FS * r.z1z2 * (2.0 * r.mu * M_U_KEV).sqrt();
        let e0 = (0.5 * b * kt).powf(2.0 / 3.0);
        let width = 4.0 / 3f64.sqrt() * (e0 * kt).sqrt();
        let f = |e: f64| -e / kt - b / e.sqrt();
        let f0 = f(e0);
        let (lo, hi) = ((e0 - 12.0 * width).max(1e-9 * e0), e0 + 12.0 * width);
        let n = 2000;
        let h = (hi - lo) / n as f64;
        let sum: f64 = (0..=n)
            .map(|i| {
                let w = if i == 0 || i == n {
                    1.0
                } else if i % 2 == 1 {
                    4.0
                } else {
                    2.0
                };
                w * (f(lo + i as f64 * h) - f0).exp()
            })
            .sum();
        f0 + (sum * h / 3.0).ln()
    }

    /// ln[α (4X/Y)²] without screening, by direct quadrature:
    /// α (4X/Y)² = ⟨σv⟩₃₄² / (⟨σv⟩₁₁ ⟨σv⟩₃₃), with ⟨σv⟩ ∝ μ^(−1/2) (kT)^(−3/2) S(0) J.
    pub fn ln_alpha_coefficient(t: f64) -> f64 {
        let term = |r: &Reaction| r.s0.ln() - 0.5 * r.mu.ln() + ln_gamow_integral(r, t);
        2.0 * term(&HE3_HE4) - term(&PP) - term(&HE3_HE3)
    }
}

/// Table of ln[α (4X/Y)²] against x = T₆^(−1/3), interpolated with cubic Hermite
/// polynomials (C¹), extrapolated linearly outside [X_MIN, X_MAX].
struct AlphaTable {
    x: Vec<f64>,
    v: Vec<f64>,
    d: Vec<f64>,
}

/// x = T₆^(−1/3) range of the table: T₆ from 1 to 200.
const X_MIN: f64 = 0.170_997_594_667_669_7;
const X_MAX: f64 = 1.0;
const ALPHA_TABLE_POINTS: usize = 161;

fn alpha_table() -> &'static AlphaTable {
    static TABLE: OnceLock<AlphaTable> = OnceLock::new();
    TABLE.get_or_init(|| {
        let n = ALPHA_TABLE_POINTS;
        let x: Vec<f64> = (0..n).map(|i| X_MIN + (X_MAX - X_MIN) * i as f64 / (n - 1) as f64).collect();
        let at = |x: f64| sf3::ln_alpha_coefficient(1e6 / x.powi(3));
        let v: Vec<f64> = x.iter().map(|&x| at(x)).collect();
        let h = 1e-4;
        let d: Vec<f64> = x.iter().map(|&x| (at(x + h) - at(x - h)) / (2.0 * h)).collect();
        AlphaTable { x, v, d }
    })
}

/// ln[α (4X/Y)²] at temperature T [K], unscreened, from the table.
pub fn ln_alpha_coefficient(t: f64) -> f64 {
    let tab = alpha_table();
    let x = 1.0 / (t * 1e-6).cbrt();
    let n = tab.x.len();
    if x <= tab.x[0] {
        return tab.v[0] + tab.d[0] * (x - tab.x[0]);
    }
    if x >= tab.x[n - 1] {
        return tab.v[n - 1] + tab.d[n - 1] * (x - tab.x[n - 1]);
    }
    let h = tab.x[1] - tab.x[0];
    let i = (((x - tab.x[0]) / h) as usize).min(n - 2);
    let s = (x - tab.x[i]) / h;
    let (s2, s3) = (s * s, s * s * s);
    tab.v[i] * (2.0 * s3 - 3.0 * s2 + 1.0)
        + h * tab.d[i] * (s3 - 2.0 * s2 + s)
        + tab.v[i + 1] * (-2.0 * s3 + 3.0 * s2)
        + h * tab.d[i + 1] * (s3 - s2)
}

#[derive(Debug, Clone, Copy, Default)]
pub struct AnalyticNuclear;

/// Cap on the weak-screening exponent. The weak-screening approximation is
/// only valid for small values; the cap also guards against overflow in
/// unphysical intermediate states.
const MAX_SCREENING_EXPONENT: f64 = 2.0;

/// ζ ρ^½ T₆^(−3/2) × 0.188: the weak-screening exponent per unit Z₁Z₂.
fn screening_base(rho: f64, t: f64, c: &Composition) -> f64 {
    // ζ² = Σ X_i Z_i (Z_i + 1) / A_i (metals approximated by CNO).
    let zeta = (2.0 * c.x + 1.5 * c.y + 4.0 * c.z).sqrt();
    0.188 * zeta * rho.sqrt() * (t * 1e-6).powf(-1.5)
}

impl Nuclear for AnalyticNuclear {
    fn name(&self) -> &'static str {
        "analytic pp (with pp II branching) + CNO + 3α (weak screening)"
    }

    fn screening(&self, rho: f64, t: f64, c: &Composition) -> Option<ScreeningExponents> {
        if !(t > 0.0 && rho > 0.0) {
            return None;
        }
        let b = screening_base(rho, t, c);
        // pp: p + p (³He + ³He and ³He + ⁴He enter only through α, with 4b);
        // CNO: p + ¹⁴N; triple-α: α + ⁸Be.
        Some(ScreeningExponents { pp: 4.0 * b, cno: 7.0 * b, triple_alpha: 8.0 * b })
    }

    fn screening_cap(&self) -> Option<f64> {
        Some(MAX_SCREENING_EXPONENT)
    }

    fn eps(&self, rho: f64, t: f64, c: &Composition) -> NuclearRates {
        if !(t > 0.0 && rho > 0.0) {
            return NuclearRates::default();
        }
        let t6 = t * 1e-6;
        let t613 = t6.cbrt();
        let t623 = t613 * t613;

        let screen_base = screening_base(rho, t, c);
        let screen = |z1z2: f64| (z1z2 * screen_base).min(MAX_SCREENING_EXPONENT).exp();

        // g₁₁: correction to the p + p rate (not the branching factor ψ).
        let g11 = 1.0 + 0.0123 * t613 + 0.0109 * t623 + 0.0009 * t6;
        let pp_1 = 2.38e6 * rho * c.x * c.x * screen(1.0) * g11 * (-33.80 / t613).exp() / t623;

        // pp II branching: α with the screening factors f₃₄² / (f₁₁ f₃₃).
        let pp_phi = if c.x > 0.0 && c.y > 0.0 {
            let n4_np = c.y / (4.0 * c.x);
            pp_phi(n4_np * n4_np * ln_alpha_coefficient(t).exp() * screen(4.0) / screen(1.0))
        } else {
            1.0
        };
        let psi = pp_phi * (1.0 - 0.021 * pp2_fraction(pp_phi));
        let pp = psi * pp_1;

        let c_cno = (1.0 + 0.0027 * t613 - 0.00778 * t623 - 0.000149 * t6).max(0.0);
        let cno = 8.67e27 * rho * c.x * c.x_cno * screen(7.0) * c_cno * (-152.28 / t613).exp() / t623;

        let t8 = t * 1e-8;
        let triple_alpha =
            5.09e11 * rho * rho * c.y.powi(3) * screen(4.0) * screen(8.0) * (-44.027 / t8).exp() / t8.powi(3);

        NuclearRates { pp, cno, triple_alpha, pp_phi }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn present_solar_centre_order_of_magnitude() {
        // Present-day solar centre: T ≈ 15.7 MK, ρ ≈ 150 g/cm³, X ≈ 0.34.
        let c = Composition::new(0.34, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear.eps(150.0, 1.57e7, &c);
        assert!(e.pp > 3.0 && e.pp < 50.0, "ε_pp = {}", e.pp);
        assert!(e.cno < e.pp, "pp should dominate in the Sun");
        assert!(e.triple_alpha < 1e-20);
    }

    #[test]
    fn cno_dominates_in_hot_cores() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear.eps(20.0, 2.8e7, &c);
        assert!(e.cno > e.pp, "pp {} vs cno {}", e.pp, e.cno);
    }

    /// α (4X/Y)² against the leading-order closed form C_α exp(−100 T₆^(−1/3)):
    /// the quadrature with constant S(0) gives C_α = 1.47, 1.46, 1.45 × 10¹⁷ at
    /// T₆ = 10, 15, 30 (physics/04, recomputed 2026-11-11).
    #[test]
    fn alpha_matches_recorded_check_values() {
        for (t6, c) in [(10.0, 1.47e17), (15.0, 1.46e17), (30.0, 1.45e17)] {
            let direct = sf3::ln_alpha_coefficient(t6 * 1e6).exp() / (-100.0 / f64::cbrt(t6)).exp();
            assert!((direct / c - 1.0).abs() < 0.01, "T6 = {t6}: C_α = {direct:.4e}");
        }
    }

    /// The table reproduces the direct quadrature to 10⁻⁴ (rule 3).
    #[test]
    fn alpha_table_reproduces_quadrature() {
        let mut worst: f64 = 0.0;
        for i in 0..=400 {
            let t6 = 10f64.powf(0.0 + 2.3 * i as f64 / 400.0);
            let d = ln_alpha_coefficient(t6 * 1e6) - sf3::ln_alpha_coefficient(t6 * 1e6);
            worst = worst.max(d.abs());
        }
        assert!(worst < 1e-4, "largest |Δ ln α| = {worst:.2e}");
    }

    #[test]
    fn branching_limits() {
        assert_eq!(pp_phi(0.0), 1.0);
        assert!((pp_phi(1e-12) - 1.0).abs() < 1e-5);
        assert!(pp_phi(1e12) <= 2.0 && pp_phi(1e12) > 2.0 - 1e-5);
        for a in [1e-6_f64, 1e-3, 0.1, 1.0, 10.0, 1e3] {
            let naive = 1.0 - a + (a * (a + 2.0)).sqrt();
            assert!((pp_phi(a) - naive).abs() < 1e-9, "α = {a}");
        }
        // φ → 1 at low temperature for any composition.
        let c = Composition::new(0.3, 0.02, 0.7).unwrap();
        assert!((AnalyticNuclear.eps(100.0, 3e6, &c).pp_phi - 1.0).abs() < 1e-6);
        // Q_pp^eff: pp I value without branching, 2.1% lower for pure pp II.
        assert_eq!(q_pp_eff(1.0), Q_PP_I);
        assert!((q_pp_eff(2.0) / Q_PP_I - 0.979).abs() < 1e-12);
    }

    /// pp II share at the centre of the present-day Sun (reported, no pass range).
    #[test]
    fn pp2_share_at_solar_centre() {
        let c = Composition::new(0.34, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear.eps(150.0, 1.57e7, &c);
        let f2 = pp2_fraction(e.pp_phi);
        eprintln!("present solar centre: φ = {:.4}, F_II = {:.3}", e.pp_phi, f2);
        assert!(e.pp_phi > 1.0 && e.pp_phi < 2.0);
    }

    #[test]
    fn finite_and_nonnegative() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        for lt in 30..=95 {
            for lr in -120..=60 {
                let t = 10f64.powf(lt as f64 / 10.0);
                let rho = 10f64.powf(lr as f64 / 10.0);
                let e = AnalyticNuclear.eps(rho, t, &c);
                assert!(e.total().is_finite() && e.total() >= 0.0);
            }
        }
    }
}
