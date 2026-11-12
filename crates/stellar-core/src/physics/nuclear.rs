//! Nuclear energy generation rates (CGS).
//!
//! * pp chain:  ε_pp = ψ ε_pp,I with ε_pp,I = (r₁₁/2) E_I / ρ, r₁₁ = ½ n_p² ⟨σv⟩₁₁ f₁₁,
//!   E_I = 26.20 MeV deposited per completed pp I chain.
//!   ψ = φ (1 − 0.021 F_II) is the pp II branching factor (pp III neglected),
//!   with φ = 1 + 2/(√(1 + 2/α) + 1), F_II = 2(φ − 1)/φ and α from its definition
//!   α = (n₄/n_p)² ⟨σv⟩₃₄² / (⟨σv⟩₁₁ ⟨σv⟩₃₃) × f₃₄²/(f₁₁ f₃₃).
//! * CNO cycle: ε_CNO = n_p n₁₄ ⟨σv⟩₁,₁₄ f₁₄ E_CNO / ρ, the p + ¹⁴N bottleneck
//!   in CN equilibrium (all catalysts counted as ¹⁴N), E_CNO = 25.03 MeV per
//!   cycle deposited.
//! * triple-α:  ε_3α  = 5.09e11 ρ² Y³ f T₈^(−3) exp(−44.027 / T₈)
//!   (Kippenhahn & Weigert 1990).
//!
//! The thermal rates ⟨σv⟩ of p + p, ³He + ³He, ³He + ⁴He and p + ¹⁴N are eq. 3
//! of Solar Fusion III (Acharya et al. 2025, Rev. Mod. Phys. 97, 035002)
//! integrated numerically with S(E) = S(0) + S′(0) E + ½ S″(0) E² (their eq. 4)
//! and the S(0), S′(0), S″(0) of their Table I; tabulated once and interpolated
//! with C¹ cubic Hermite polynomials.
//!
//! Each includes a weak-screening factor f = exp(0.188 Z₁Z₂ ζ ρ^½ T₆^(−3/2)).

use std::sync::OnceLock;

use super::{Nuclear, NuclearRates, ScreeningExponents};
use crate::composition::Composition;
use crate::constants::M_U;

/// Energy deposited per gram of hydrogen burnt through the pp I chain
/// [erg g⁻¹] (26.20 MeV per ⁴He nucleus after neutrino losses).
pub const Q_PP_I: f64 = 6.27e18;
/// Energy deposited per gram of hydrogen burnt through the CNO cycle
/// [erg g⁻¹] (25.03 MeV per ⁴He nucleus after neutrino losses).
pub const Q_CNO: f64 = 5.99e18;
/// Energy released per gram of helium burnt to carbon by the triple-α
/// reaction [erg g⁻¹] (7.275 MeV per ¹²C nucleus).
pub const Q_HE: f64 = 5.85e17;

/// Atomic masses [u] (AME): ¹H, ⁴He, ¹⁴N.
const A_H1: f64 = 1.007_825;
const A_HE4: f64 = 4.002_603;
const A_N14: f64 = 14.003_074;

/// Energy per completed pp I chain and per CNO cycle [erg], from Q_PP_I and
/// Q_CNO, so that the burning step and the rates count hydrogen identically.
const E_PP_I: f64 = Q_PP_I * 4.0 * A_H1 * M_U;
const E_CNO: f64 = Q_CNO * 4.0 * A_H1 * M_U;

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

/// Thermal reaction rates from Solar Fusion III (Acharya et al. 2025).
pub mod sf3 {
    use crate::constants::{C, K_B, M_U};

    /// Fine-structure constant.
    const ALPHA_FS: f64 = 1.0 / 137.035_999_084;
    /// 1 MeV [erg]; 1 b [cm²].
    const MEV: f64 = 1.602_176_634e-6;
    const BARN: f64 = 1e-24;

    /// A non-resonant reaction with S(E) = S(0) + S′(0) E + ½ S″(0) E².
    #[derive(Clone, Copy, Debug)]
    pub struct Reaction {
        pub z1z2: f64,
        /// Reduced mass [u].
        pub mu: f64,
        /// S(0) [MeV b], S′(0) [b], S″(0) [b MeV⁻¹] (Acharya et al. 2025, Table I).
        pub s0: f64,
        pub s1: f64,
        pub s2: f64,
    }

    const fn reduced(a: f64, b: f64) -> f64 {
        a * b / (a + b)
    }

    /// ¹H(p,e⁺ν)²H. Table I; S′/S = 11.0 MeV⁻¹, S″/S = 242 MeV⁻² (§III.F).
    pub const PP: Reaction =
        Reaction { z1z2: 1.0, mu: reduced(1.007_825, 1.007_825), s0: 4.09e-25, s1: 4.5e-24, s2: 9.9e-23 };
    /// ³He(³He,2p)⁴He. Table I; S(E) = 5.21 − 4.90 E + 11.21 E² MeV b (§IV.B).
    pub const HE3_HE3: Reaction =
        Reaction { z1z2: 4.0, mu: reduced(3.016_029, 3.016_029), s0: 5.21, s1: -4.90, s2: 22.42 };
    /// ³He(α,γ)⁷Be. Table I; S′/S = −0.54 MeV⁻¹ for E = 0–500 keV, S″ consistent
    /// with zero (§V.E).
    pub const HE3_HE4: Reaction =
        Reaction { z1z2: 4.0, mu: reduced(3.016_029, 4.002_603), s0: 5.61e-4, s1: -3.03e-4, s2: 0.0 };
    /// ¹⁴N(p,γ)¹⁵O. Table I gives S(0) only (sum of the partial S(0), Table XI);
    /// resonances are not included (they matter above about T₆ = 100).
    pub const N14_P: Reaction =
        Reaction { z1z2: 7.0, mu: reduced(1.007_825, 14.003_074), s0: 1.68e-3, s1: 0.0, s2: 0.0 };

    /// ln ⟨σv⟩ [cm³ s⁻¹] at temperature T [K], eq. 3 of Solar Fusion III:
    /// ⟨σv⟩ = √(8/(πμ)) (kT)^(−3/2) ∫ S(E) exp(−E/kT − b/√E) dE,
    /// b = π α Z₁Z₂ √(2μc²), by composite Simpson quadrature around the Gamow
    /// peak (±12 widths, 2000 intervals). S(E) is kept non-negative.
    pub fn ln_sigma_v(r: &Reaction, t: f64) -> f64 {
        let mu = r.mu * M_U;
        let kt = K_B * t;
        let b = std::f64::consts::PI * ALPHA_FS * r.z1z2 * (2.0 * mu * C * C).sqrt();
        let e0 = (0.5 * b * kt).powf(2.0 / 3.0);
        let width = 4.0 / 3f64.sqrt() * (e0 * kt).sqrt();
        let f = |e: f64| -e / kt - b / e.sqrt();
        let f0 = f(e0);
        let s = |e: f64| {
            let em = e / MEV;
            (r.s0 + r.s1 * em + 0.5 * r.s2 * em * em).max(0.0) / r.s0
        };
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
                let e = lo + i as f64 * h;
                w * s(e) * (f(e) - f0).exp()
            })
            .sum();
        0.5 * (8.0 / (std::f64::consts::PI * mu)).ln() - 1.5 * kt.ln()
            + (r.s0 * MEV * BARN).ln()
            + f0
            + (sum * h / 3.0).ln()
    }
}

/// Table of a smooth function of x = T₆^(−1/3), interpolated with cubic
/// Hermite polynomials (C¹), extrapolated linearly outside [X_MIN, X_MAX].
struct Table {
    x: Vec<f64>,
    v: Vec<f64>,
    d: Vec<f64>,
}

/// x = T₆^(−1/3) range of the tables: T₆ from 1 to 200.
const X_MIN: f64 = 0.170_997_594_667_669_7;
const X_MAX: f64 = 1.0;
const TABLE_POINTS: usize = 161;

impl Table {
    fn build(f: impl Fn(f64) -> f64) -> Self {
        let n = TABLE_POINTS;
        let x: Vec<f64> = (0..n).map(|i| X_MIN + (X_MAX - X_MIN) * i as f64 / (n - 1) as f64).collect();
        let at = |x: f64| f(1e6 / x.powi(3));
        let v: Vec<f64> = x.iter().map(|&x| at(x)).collect();
        let h = 1e-4;
        let d: Vec<f64> = x.iter().map(|&x| (at(x + h) - at(x - h)) / (2.0 * h)).collect();
        Self { x, v, d }
    }

    fn at_t(&self, t: f64) -> f64 {
        let x = 1.0 / (t * 1e-6).cbrt();
        let n = self.x.len();
        if x <= self.x[0] {
            return self.v[0] + self.d[0] * (x - self.x[0]);
        }
        if x >= self.x[n - 1] {
            return self.v[n - 1] + self.d[n - 1] * (x - self.x[n - 1]);
        }
        let h = self.x[1] - self.x[0];
        let i = (((x - self.x[0]) / h) as usize).min(n - 2);
        let s = (x - self.x[i]) / h;
        let (s2, s3) = (s * s, s * s * s);
        self.v[i] * (2.0 * s3 - 3.0 * s2 + 1.0)
            + h * self.d[i] * (s3 - 2.0 * s2 + s)
            + self.v[i + 1] * (-2.0 * s3 + 3.0 * s2)
            + h * self.d[i + 1] * (s3 - s2)
    }
}

/// Tabulated ln ⟨σv⟩ of the four reactions: [pp, ³He³He, ³He⁴He, p¹⁴N].
fn tables() -> &'static [Table; 4] {
    static TABLES: OnceLock<[Table; 4]> = OnceLock::new();
    TABLES.get_or_init(|| {
        [sf3::PP, sf3::HE3_HE3, sf3::HE3_HE4, sf3::N14_P].map(|r| Table::build(move |t| sf3::ln_sigma_v(&r, t)))
    })
}

/// ln ⟨σv⟩ [cm³ s⁻¹] of p + p, ³He + ³He, ³He + ⁴He and p + ¹⁴N at T [K], from the tables.
pub fn ln_sigma_v(t: f64) -> [f64; 4] {
    let tab = tables();
    [tab[0].at_t(t), tab[1].at_t(t), tab[2].at_t(t), tab[3].at_t(t)]
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
        "Solar Fusion III pp (with pp II branching) + CNO, analytic 3α (weak screening)"
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
        let screen_base = screening_base(rho, t, c);
        let screen = |z1z2: f64| (z1z2 * screen_base).min(MAX_SCREENING_EXPONENT).exp();
        let [l11, l33, l34, l114] = ln_sigma_v(t);

        // Number densities [cm⁻³].
        let n_p = rho * c.x / (A_H1 * M_U);
        let n_14 = rho * c.x_cno / (A_N14 * M_U);

        // pp I: r₁₁ = ½ n_p² ⟨σv⟩₁₁ f₁₁ reactions per cm³ s, one chain per two reactions.
        let r11 = 0.5 * n_p * n_p * l11.exp() * screen(1.0);
        let pp_1 = 0.5 * r11 * E_PP_I / rho;

        // pp II branching: α with the screening factors f₃₄² / (f₁₁ f₃₃) = f₄/f₁.
        let pp_phi = if c.x > 0.0 && c.y > 0.0 {
            let n4_np = c.y * A_H1 / (c.x * A_HE4);
            pp_phi(n4_np * n4_np * (2.0 * l34 - l11 - l33).exp() * screen(4.0) / screen(1.0))
        } else {
            1.0
        };
        let psi = pp_phi * (1.0 - 0.021 * pp2_fraction(pp_phi));
        let pp = psi * pp_1;

        // CNO: one cycle per p + ¹⁴N reaction in CN equilibrium.
        let cno = n_p * n_14 * l114.exp() * screen(7.0) * E_CNO / rho;

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

    /// The recorded check of physics/04: with constant S(0) and integer mass
    /// numbers, α (4X/Y)² / exp(−100 T₆^(−1/3)) = C_α = 1.47, 1.46, 1.45 × 10¹⁷ at
    /// T₆ = 10, 15, 30 (computed for those pages 2026-11-11).
    #[test]
    fn alpha_matches_recorded_check_values() {
        let flat = |r: sf3::Reaction, mu: f64| sf3::Reaction { mu, s1: 0.0, s2: 0.0, ..r };
        let (pp, h33, h34) = (flat(sf3::PP, 0.5), flat(sf3::HE3_HE3, 1.5), flat(sf3::HE3_HE4, 12.0 / 7.0));
        for (t6, c) in [(10.0, 1.47e17), (15.0, 1.46e17), (30.0, 1.45e17)] {
            let t = t6 * 1e6;
            let ln_a = 2.0 * sf3::ln_sigma_v(&h34, t) - sf3::ln_sigma_v(&pp, t) - sf3::ln_sigma_v(&h33, t);
            let c_alpha = ln_a.exp() / (-100.0 / f64::cbrt(t6)).exp();
            assert!((c_alpha / c - 1.0).abs() < 0.01, "T6 = {t6}: C_α = {c_alpha:.4e}");
        }
    }

    /// The tables reproduce the direct quadrature to 10⁻⁴ in ln⟨σv⟩ (rule 3).
    #[test]
    fn tables_reproduce_quadrature() {
        let reactions = [sf3::PP, sf3::HE3_HE3, sf3::HE3_HE4, sf3::N14_P];
        let mut worst = [0.0_f64; 4];
        for i in 0..=400 {
            let t = 1e6 * 10f64.powf(2.3 * i as f64 / 400.0);
            let tab = ln_sigma_v(t);
            for j in 0..4 {
                worst[j] = worst[j].max((tab[j] - sf3::ln_sigma_v(&reactions[j], t)).abs());
            }
        }
        assert!(worst.iter().all(|w| *w < 1e-4), "largest |Δ ln⟨σv⟩| = {worst:.2?}");
    }

    /// Against the 1990 fits that were replaced (unscreened, pp I only):
    /// physics/04 estimated the Solar Fusion III p + p rate about 8% above the
    /// fit and the ¹⁴N(p,γ) rate 0.52–0.56 of it at T₆ = 15–30. Reported; loose bounds.
    #[test]
    fn comparison_with_the_1990_fits() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        for t6 in [10.0_f64, 15.0, 20.0, 30.0] {
            let (t613, t623) = (t6.cbrt(), t6.cbrt().powi(2));
            let rho = 1.0;
            let g11 = 1.0 + 0.0123 * t613 + 0.0109 * t623 + 0.0009 * t6;
            let pp_fit = 2.38e6 * rho * c.x * c.x * g11 * (-33.80 / t613).exp() / t623;
            let c_cno = 1.0 + 0.0027 * t613 - 0.00778 * t623 - 0.000149 * t6;
            let cno_fit = 8.67e27 * rho * c.x * c.x_cno * c_cno * (-152.28 / t613).exp() / t623;
            let [l11, _, _, l114] = ln_sigma_v(t6 * 1e6);
            let n_p = rho * c.x / (A_H1 * M_U);
            let pp_new = 0.25 * n_p * n_p * l11.exp() * E_PP_I / rho;
            let cno_new = n_p * rho * c.x_cno / (A_N14 * M_U) * l114.exp() * E_CNO / rho;
            let (rp, rc) = (pp_new / pp_fit, cno_new / cno_fit);
            eprintln!("T6 = {t6}: SF III / 1990 fit: pp {rp:.4}, CNO {rc:.4}");
            assert!(rp > 1.0 && rp < 1.2, "pp ratio {rp}");
            assert!(rc > 0.45 && rc < 0.65, "CNO ratio {rc}");
        }
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
