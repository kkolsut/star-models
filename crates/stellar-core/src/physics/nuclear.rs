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
//! * ¹²C(α,γ)¹⁶O: ε_Cα = f_Cα ρ X_C Y N_A⟨σv⟩ f Q_Cα / 4 with the Caughlan &
//!   Fowler (1988) fit, on the carbon made by helium burning (plan 9.5).
//!
//! The thermal rates ⟨σv⟩ of p + p, ³He + ³He, ³He + ⁴He and p + ¹⁴N are eq. 3
//! of Solar Fusion III (Acharya et al. 2025, Rev. Mod. Phys. 97, 035002)
//! integrated numerically with S(E) = S(0) + S′(0) E + ½ S″(0) E² (their eq. 4)
//! and the S(0), S′(0), S″(0) of their Table I; tabulated once and interpolated
//! with C¹ cubic Hermite polynomials.
//!
//! Each includes a screening factor f = e^H from [`super::screening`]: Chugunov &
//! DeWitt (2009) with the electron factor of Potekhin & Chabrier (2012) for the
//! binary reactions, the static two-stage form for triple-α (plan 9.7). H is
//! not capped. Its weak limit is the former Salpeter factor
//! exp(0.188 Z₁Z₂ ζ ρ^½ T₆^(−3/2)).

use std::sync::OnceLock;

use super::screening::{Pair, Plasma};
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
/// Energy released per gram of carbon burnt to oxygen by ¹²C(α,γ)¹⁶O
/// [erg g⁻¹] (7.162 MeV per reaction, AME2020; `docs/physics/04`).
pub const Q_C_ALPHA: f64 = 5.759e17;

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

/// ¹²C(α,γ)¹⁶O: N_A⟨σv⟩ [cm³ mol⁻¹ s⁻¹] of Caughlan & Fowler (1988), as
/// reproduced term by term in the ORNL Fortran version of the compilation
/// (reaction 62) and in MESA's `rate_c12ag_fxt` (`docs/physics/04`, decided
/// 2026-11-12), without MESA's factor 1.7.
pub fn c12_alpha_na_sigma_v(t: f64) -> f64 {
    let t9 = t * 1e-9;
    if !(t9 > 0.0) {
        return 0.0;
    }
    let t913 = t9.cbrt();
    let t923 = t913 * t913;
    let t9m23 = 1.0 / t923;
    let t92 = t9 * t9;
    let gamow = (-32.120 / t913).exp();
    1.04e8 / (t92 * (1.0 + 0.0489 * t9m23).powi(2)) * gamow * (-(t9 / 3.496).powi(2)).exp()
        + 1.76e8 / (t92 * (1.0 + 0.2654 * t9m23).powi(2)) * gamow
        + 1.25e3 * t9.powf(-1.5) * (-27.499 / t9).exp()
        + 1.43e-2 * t9.powi(5) * (-15.541 / t9).exp()
}

/// The analytic nuclear network: Solar Fusion III pp chain (with pp II
/// branching) and CNO bottleneck, triple-α, and ¹²C(α,γ)¹⁶O, screened by the
/// all-regime prescription of [`super::screening`] (plan 9.7).
#[derive(Debug, Clone, Copy)]
pub struct AnalyticNuclear {
    /// Multiplier f_Cα of the ¹²C(α,γ)¹⁶O rate (1: CF88; MESA's
    /// `rate_c12ag_fxt` uses 1.7). The spread is the present uncertainty.
    pub f_c_alpha: f64,
}

impl Default for AnalyticNuclear {
    fn default() -> Self {
        Self { f_c_alpha: 1.0 }
    }
}

/// Reduced masses [u] of the screened pairs.
const MU_PP: f64 = A_H1 * A_H1 / (A_H1 + A_H1);
const MU_33: f64 = 3.016_029 * 3.016_029 / (2.0 * 3.016_029);
const MU_34: f64 = 3.016_029 * A_HE4 / (3.016_029 + A_HE4);
const MU_P14: f64 = A_H1 * A_N14 / (A_H1 + A_N14);
const MU_C_ALPHA: f64 = 12.0 * A_HE4 / (12.0 + A_HE4);

/// Ceiling on any one channel's ε [erg g⁻¹ s⁻¹], applied to ln ε. Not a cap on
/// the screening exponent (plan 9.7 removed that): a guard against overflow
/// reached only in unphysical trial states of the Newton iteration (for
/// example 10³ K at 10⁻² g cm⁻³ with classical electrons, where H ≈ 10³),
/// many orders of magnitude above any rate of a real model.
pub const EPS_CEILING: f64 = 1e30;

/// e^{ln ε}, limited to [`EPS_CEILING`]; 0 for ln ε = −∞.
fn bounded(ln_eps: f64) -> f64 {
    ln_eps.min(EPS_CEILING.ln()).exp()
}

/// Screening exponents ln f of every reaction at one state.
#[derive(Clone, Copy, Debug, Default)]
struct Screens {
    pp: f64,
    he33: f64,
    he34: f64,
    cno: f64,
    triple_alpha: f64,
    c_alpha: f64,
}

/// Below this temperature the triple-α rate is under 10⁻⁵⁰ erg g⁻¹ s⁻¹ at any
/// stellar density, and its screening factor is not evaluated (f = 1).
const T_MIN_TRIPLE_ALPHA: f64 = 3e7;

/// The pairs: p + p, ³He + ³He, ³He + ⁴He, p + ¹⁴N, α + α and ⁸Be + α
/// (static), ¹²C + α.
fn pairs() -> &'static [Pair; 7] {
    static PAIRS: OnceLock<[Pair; 7]> = OnceLock::new();
    PAIRS.get_or_init(|| {
        [
            Pair::new(1.0, 1.0, Some(MU_PP)),
            Pair::new(2.0, 2.0, Some(MU_33)),
            Pair::new(2.0, 2.0, Some(MU_34)),
            Pair::new(1.0, 7.0, Some(MU_P14)),
            Pair::new(2.0, 2.0, None),
            Pair::new(4.0, 2.0, None),
            Pair::new(2.0, 6.0, Some(MU_C_ALPHA)),
        ]
    })
}

/// Last states evaluated: the finite-difference Jacobian repeats (ρ, T)
/// whenever it perturbs r or L.
type ScreenKey = [u64; 8];
thread_local! {
    static SCREEN_MEMO: std::cell::RefCell<Vec<(ScreenKey, Screens)>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Screening exponents of the reactions that can run in composition `c`
/// (hydrogen burning needs X > 0, ¹²C(α,γ) produced carbon, triple-α
/// T ≥ 3 × 10⁷ K and Y > 0); the others are 0.
fn screens(rho: f64, t: f64, eta: f64, c: &Composition) -> Screens {
    let key = [rho, t, eta, c.x, c.y, c.z, c.x_c, c.x_o].map(f64::to_bits);
    if let Some(s) = SCREEN_MEMO.with(|m| m.borrow().iter().find(|(k, _)| *k == key).map(|(_, s)| *s)) {
        return s;
    }
    let s = match Plasma::new(rho, t, eta, c) {
        None => Screens::default(),
        Some(pl) => {
            let p = pairs();
            let (hydrogen, helium) = (c.x > 0.0, c.y > 0.0);
            let h = |on: bool, i: usize| if on { pl.h_pair(&p[i]) } else { 0.0 };
            let hot = helium && t >= T_MIN_TRIPLE_ALPHA;
            Screens {
                pp: h(hydrogen, 0),
                he33: h(hydrogen && helium, 1),
                he34: h(hydrogen && helium, 2),
                cno: h(hydrogen, 3),
                triple_alpha: if hot { pl.h_pair(&p[4]) + pl.h_pair(&p[5]) } else { 0.0 },
                c_alpha: h(helium && c.x_c > 0.0, 6),
            }
        }
    };
    SCREEN_MEMO.with(|m| {
        let mut m = m.borrow_mut();
        if m.len() >= 8 {
            m.remove(0);
        }
        m.push((key, s));
    });
    s
}

impl Nuclear for AnalyticNuclear {
    fn name(&self) -> &'static str {
        "Solar Fusion III pp (with pp II branching) + CNO, analytic 3α, CF88 ¹²C(α,γ) (Chugunov–DeWitt screening)"
    }

    fn screening_eta(&self, rho: f64, t: f64, eta: f64, c: &Composition) -> Option<ScreeningExponents> {
        if !(t > 0.0 && rho > 0.0) {
            return None;
        }
        // Every reaction, also those that cannot run in this composition.
        let p = pairs();
        let pl = Plasma::new(rho, t, eta, c)?;
        let h = |i: usize| pl.h_pair(&p[i]);
        let s = Screens { pp: h(0), he33: h(1), he34: h(2), cno: h(3), triple_alpha: h(4) + h(5), c_alpha: h(6) };
        // pp chain: the largest of p + p and the ³He reactions entering φ.
        Some(ScreeningExponents {
            pp: s.pp.max(s.he33).max(s.he34),
            cno: s.cno,
            triple_alpha: s.triple_alpha,
            c_alpha: s.c_alpha,
        })
    }

    fn eps_eta(&self, rho: f64, t: f64, eta: f64, c: &Composition) -> NuclearRates {
        if !(t > 0.0 && rho > 0.0) {
            return NuclearRates::default();
        }
        let s = screens(rho, t, eta, c);
        let [l11, l33, l34, l114] = ln_sigma_v(t);

        // Number densities [cm⁻³].
        let n_p = rho * c.x / (A_H1 * M_U);
        let n_14 = rho * c.x_cno / (A_N14 * M_U);

        // pp I: r₁₁ = ½ n_p² ⟨σv⟩₁₁ f₁₁ reactions per cm³ s, one chain per two
        // reactions; ε_pp,I = r₁₁ E_I/(2ρ). Rates are formed as logarithms.
        let ln_pp_1 = (0.25 * n_p * n_p * E_PP_I / rho).ln() + l11 + s.pp;

        // pp II branching: α with the screening factors f₃₄² / (f₁₁ f₃₃).
        let pp_phi = if c.x > 0.0 && c.y > 0.0 {
            let n4_np = c.y * A_H1 / (c.x * A_HE4);
            pp_phi(n4_np * n4_np * (2.0 * (l34 + s.he34) - (l11 + s.pp) - (l33 + s.he33)).exp())
        } else {
            1.0
        };
        let psi = pp_phi * (1.0 - 0.021 * pp2_fraction(pp_phi));
        let pp = bounded(psi.ln() + ln_pp_1);

        // CNO: one cycle per p + ¹⁴N reaction in CN equilibrium.
        let cno = bounded((n_p * n_14 * E_CNO / rho).ln() + l114 + s.cno);

        let t8 = t * 1e-8;
        let triple_alpha =
            bounded((5.09e11 * rho * rho * c.y.powi(3) / t8.powi(3)).ln() + s.triple_alpha - 44.027 / t8);

        // ¹²C(α,γ)¹⁶O: ε = ρ X_C Y N_A⟨σv⟩ Q_Cα / 4 (Q_Cα per gram of carbon),
        // on the carbon produced by helium burning.
        let c_alpha = if c.x_c > 0.0 && c.y > 0.0 {
            bounded((self.f_c_alpha * rho * c.x_c * c.y * c12_alpha_na_sigma_v(t) * Q_C_ALPHA / 4.0).ln() + s.c_alpha)
        } else {
            0.0
        };

        NuclearRates { pp, cno, triple_alpha, c_alpha, pp_phi }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn present_solar_centre_order_of_magnitude() {
        // Present-day solar centre: T ≈ 15.7 MK, ρ ≈ 150 g/cm³, X ≈ 0.34.
        let c = Composition::new(0.34, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear::default().eps(150.0, 1.57e7, &c);
        assert!(e.pp > 3.0 && e.pp < 50.0, "ε_pp = {}", e.pp);
        assert!(e.cno < e.pp, "pp should dominate in the Sun");
        assert!(e.triple_alpha < 1e-20);
    }

    #[test]
    fn cno_dominates_in_hot_cores() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear::default().eps(20.0, 2.8e7, &c);
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
        assert!((AnalyticNuclear::default().eps(100.0, 3e6, &c).pp_phi - 1.0).abs() < 1e-6);
        // Q_pp^eff: pp I value without branching, 2.1% lower for pure pp II.
        assert_eq!(q_pp_eff(1.0), Q_PP_I);
        assert!((q_pp_eff(2.0) / Q_PP_I - 0.979).abs() < 1e-12);
    }

    /// pp II share at the centre of the present-day Sun (reported, no pass range).
    #[test]
    fn pp2_share_at_solar_centre() {
        let c = Composition::new(0.34, 0.02, 0.7).unwrap();
        let e = AnalyticNuclear::default().eps(150.0, 1.57e7, &c);
        let f2 = pp2_fraction(e.pp_phi);
        eprintln!("present solar centre: φ = {:.4}, F_II = {:.3}", e.pp_phi, f2);
        assert!(e.pp_phi > 1.0 && e.pp_phi < 2.0);
    }

    /// Plan 9.7: at the centre of the present Sun (classical electrons, as with
    /// the ideal EOS) the all-regime screening changes the rates from the
    /// former weak (Salpeter) factors, ζ² = 2X + 1.5Y + 4Z, by less than the
    /// stated tolerance: 1% for pp, 15% for CNO. Measured −0.8% and −12%. The
    /// CNO change is the intermediate coupling of ¹⁴N itself: Γ_N = 7^{5/3} Γ_e
    /// = 1.7 at Γ_e = 0.067, where the free energy f₀ is well below its
    /// Debye–Hückel form; the rigid-background exponent of p + ¹⁴N is 0.156
    /// against 0.269 in the weak limit, the electron factor 1.33 against 1.27.
    #[test]
    fn solar_centre_rates_against_weak_screening() {
        let c = Composition::new(0.34, 0.02, 0.7).unwrap();
        let (rho, t) = (150.0, 1.57e7);
        let new = AnalyticNuclear::default().eps(rho, t, &c);
        let b = 0.188 * (2.0 * c.x + 1.5 * c.y + 4.0 * c.z).sqrt() * rho.sqrt() * (t * 1e-6).powf(-1.5);
        let [l11, l33, l34, l114] = ln_sigma_v(t);
        let n_p = rho * c.x / (A_H1 * M_U);
        let n4_np = c.y * A_H1 / (c.x * A_HE4);
        let phi = pp_phi(n4_np * n4_np * (2.0 * l34 - l11 - l33).exp() * (4.0 * b).exp() / b.exp());
        let pp_old = phi * (1.0 - 0.021 * pp2_fraction(phi)) * 0.25 * n_p * n_p * (l11 + b).exp() * E_PP_I / rho;
        let cno_old = n_p * rho * c.x_cno / (A_N14 * M_U) * (l114 + 7.0 * b).exp() * E_CNO / rho;
        let (dpp, dcno) = (new.pp / pp_old - 1.0, new.cno / cno_old - 1.0);
        eprintln!("solar centre, new/old − 1: pp {dpp:+.4e}, CNO {dcno:+.4e}");
        assert!(dpp.abs() < 0.01 && dcno.abs() < 0.15);
    }

    /// ¹²C(α,γ)¹⁶O: positive, increasing with T over helium burning, zero
    /// without carbon; ε/Q_Cα is the carbon burnt per gram and second,
    /// 12 m_u n_C n_α ⟨σv⟩ / ρ.
    #[test]
    fn carbon_alpha_rate() {
        let ts = [5e7, 1e8, 2e8, 3e8, 5e8, 1e9];
        let r: Vec<f64> = ts.iter().map(|&t| c12_alpha_na_sigma_v(t)).collect();
        assert!(r.iter().all(|v| *v > 0.0 && v.is_finite()));
        assert!(r.windows(2).all(|w| w[1] > w[0]), "{r:?}");
        let pure_he = Composition { x: 0.0, y: 0.98, z: 0.02, x_cno: 0.014, x_c: 0.0, x_o: 0.0 };
        assert_eq!(AnalyticNuclear::default().eps(1e4, 2e8, &pure_he).c_alpha, 0.0);
        let c = Composition { x_c: 0.4, y: 0.58, ..pure_he };
        let (rho, t) = (1e4, 2e8);
        let nuc = AnalyticNuclear { f_c_alpha: 1.0 };
        let eps = nuc.eps(rho, t, &c).c_alpha;
        let h = Plasma::new(rho, t, f64::NAN, &c).unwrap().h(2.0, 6.0, MU_C_ALPHA);
        let n_c = rho * c.x_c / (12.0 * M_U);
        let n_a = rho * c.y / (4.0 * M_U);
        let sv = c12_alpha_na_sigma_v(t) * M_U; // N_A⟨σv⟩ / N_A
        let want = 12.0 * M_U * n_c * n_a * sv * h.exp() / rho;
        assert!((eps / Q_C_ALPHA / want - 1.0).abs() < 1e-12, "{} vs {want}", eps / Q_C_ALPHA);
        // The multiplier scales the rate.
        let doubled = AnalyticNuclear { f_c_alpha: 2.0 }.eps(rho, t, &c).c_alpha;
        assert!((doubled / eps - 2.0).abs() < 1e-12);
        // Q_Cα from the AME2020 Q value, 7.1619 MeV per ¹²C.
        assert!((7.1619e6 * crate::constants::EV / (12.0 * M_U) / Q_C_ALPHA - 1.0).abs() < 2e-4);
    }

    #[test]
    fn finite_and_nonnegative() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        for lt in 30..=95 {
            for lr in -120..=60 {
                let t = 10f64.powf(lt as f64 / 10.0);
                let rho = 10f64.powf(lr as f64 / 10.0);
                let e = AnalyticNuclear::default().eps(rho, t, &c);
                assert!(e.total().is_finite() && e.total() >= 0.0);
            }
        }
    }
}
