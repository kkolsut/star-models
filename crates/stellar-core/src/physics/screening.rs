//! Electron screening for all regimes (plan 9.7, `docs/physics/04`,
//! "Screening prescription for all regimes", decided 2026-11-11/12).
//!
//! *Binary reactions.* Chugunov & DeWitt (2009, eq. A4) with the one-component
//! plasma free energy of Potekhin & Chabrier (their eq. 24 as printed in
//! `docs/physics/04`):
//!
//! H_ij = (C_ij + Γ_ij²)/(1 + Γ_ij²) · [f₀(Γ_i/t_ij) + f₀(Γ_j/t_ij) − f₀(Γ_c/t_ij)],
//! C_ij = 3 Z_iZ_j √(⟨Z²⟩/⟨Z⟩) / (Z_c^{5/2} − Z_i^{5/2} − Z_j^{5/2}),
//!
//! with Z_c = Z_i + Z_j, Γ_s = Z_s^{5/3} Γ_e, Γ_ij = 2Z_iZ_j/(Z_i^{1/3} + Z_j^{1/3}) Γ_e
//! and the quantum correction t_ij of Chugunov & DeWitt (2009, eq. 25; the
//! same as eq. 22 of Potekhin & Chabrier 2012). ⟨Z⟩, ⟨Z²⟩ are averages over
//! the actual ions (H, He, the GN93 metals, produced C and O).
//!
//! *Electrons.* The rigid-background exponent is multiplied by the factor of
//! Potekhin & Chabrier (2012, eqs. 20–21),
//! √(1 + k_TF² D̃_ion²), D̃_ion² = D_ion² + (0.6 a_ij)², k_TF² = 4πe² ∂n_e/∂μ_e.
//! Potekhin & Chabrier used the factor with their own rigid-background
//! exponent; applying it to eq. A4 of Chugunov & DeWitt (2009) is a small
//! extrapolation of the source, decided in `docs/physics/04` and checked by
//! the tests below in both limits. ∂n_e/∂μ_e comes from the electron model of
//! the EOS: classical (n_e/kT) when the EOS has no Fermi–Dirac electrons (η
//! not given), the Fermi–Dirac integrals at the EOS's η otherwise.
//!
//! *Triple-α.* Static two-stage form, t_ij = 1: ln f_3α = H_{2,2} + H_{4,2}
//! (the product of pair factors, as MESA's `eval_screen_triple`). Exact in the
//! weak limit (12e²/(kTD)) and in the strong limit; an approximation at
//! intermediate coupling.
//!
//! H is returned as a logarithm and not capped.

use super::eos_fd::electron_gas;
use crate::composition::Composition;
use crate::constants::{C, E_CHARGE, FOUR_PI, HBAR, K_B, M_E, M_U};

/// Coefficients of the free-energy fit f₀(Γ); A₃ = −√3/2 − A₁/√A₂.
const A1: f64 = -0.907;
const A2: f64 = 0.629_54;
const A3: f64 = 0.277_104_595_662_215_8;
const B1: f64 = 0.004_56;
const B2: f64 = 211.6;
const B3: f64 = -1e-4;
const B4: f64 = 0.004_62;

/// Free energy of a one-component plasma per ion in units of kT, f₀(Γ)
/// (negative; −Γ^{3/2}/√3 for Γ → 0, ≈ −0.907 Γ for Γ ≫ 1).
pub fn f0(g: f64) -> f64 {
    if !(g > 0.0) {
        return 0.0;
    }
    // For small Γ the bracketed differences cancel; use their series.
    let s = g.sqrt();
    let first = if g < 1e-4 * A2 {
        // √(Γ(A₂+Γ)) − A₂ asinh √(Γ/A₂) = (2/3) Γ^{3/2}/√A₂ (1 − 3Γ/(10A₂) + …)
        2.0 / 3.0 * g * s / A2.sqrt() * (1.0 - 0.3 * g / A2)
    } else {
        (g * (A2 + g)).sqrt() - A2 * (g / A2).sqrt().asinh()
    };
    let second = if s < 1e-4 { g * s / 3.0 * (1.0 - 0.6 * g) } else { s - s.atan() };
    let third =
        if g < 1e-6 * B2 { g * g / (2.0 * B2) * (1.0 - 2.0 * g / (3.0 * B2)) } else { g - B2 * (g / B2).ln_1p() };
    A1 * first + 2.0 * A3 * second + B1 * third + 0.5 * B3 * (g * g / B4).ln_1p()
}

/// The charge-dependent constants of one reacting pair, computed once.
#[derive(Clone, Copy, Debug)]
pub struct Pair {
    /// Z_i^{5/3}, Z_j^{5/3}, Z_c^{5/3}.
    z53: [f64; 3],
    /// C_ij / √(⟨Z²⟩/⟨Z⟩) = 3Z_iZ_j / (Z_c^{5/2} − Z_i^{5/2} − Z_j^{5/2}).
    c_ij: f64,
    /// Γ_ij / Γ_e = 2Z_iZ_j / (Z_i^{1/3} + Z_j^{1/3}).
    g_ij: f64,
    /// y_ij = 4Z_iZ_j/Z_c² and its powers 0.14, 0.19.
    y: f64,
    y014: f64,
    y019: f64,
    /// τ_ij T^{1/3} [K^{1/3}] (None: static form, t_ij = 1).
    tau_t13: Option<f64>,
    /// a_ij / a_e = (Z_i^{1/3} + Z_j^{1/3})/2.
    a_ij: f64,
}

impl Pair {
    /// Pair of charges `zi`, `zj`; `mu_amu` is the reduced mass [u] for the
    /// quantum correction t_ij, `None` for the static form.
    pub fn new(zi: f64, zj: f64, mu_amu: Option<f64>) -> Self {
        let zc = zi + zj;
        let y = 4.0 * zi * zj / (zc * zc);
        let tau_t13 = mu_amu.map(|mu| {
            (27.0 * std::f64::consts::PI.powi(2) * mu * M_U * zi * zi * zj * zj * E_CHARGE.powi(4)
                / (2.0 * K_B * HBAR * HBAR))
                .cbrt()
        });
        Self {
            z53: [zi, zj, zc].map(|z: f64| z.powf(5.0 / 3.0)),
            c_ij: 3.0 * zi * zj / (zc.powf(2.5) - zi.powf(2.5) - zj.powf(2.5)),
            g_ij: 2.0 * zi * zj / (zi.cbrt() + zj.cbrt()),
            y,
            y014: y.powf(0.14),
            y019: y.powf(0.19),
            tau_t13,
            a_ij: 0.5 * (zi.cbrt() + zj.cbrt()),
        }
    }
}

/// The plasma at (ρ, T) for one composition: the quantities shared by all
/// reactions.
#[derive(Clone, Copy, Debug)]
pub struct Plasma {
    /// kT [erg].
    pub kt: f64,
    /// T^{1/3} [K^{1/3}].
    t13: f64,
    /// Electron-sphere radius a_e = (3/(4π n_e))^{1/3} [cm].
    pub a_e: f64,
    /// Γ_e = e²/(a_e kT).
    pub gamma_e: f64,
    /// ⟨Z²⟩/⟨Z⟩ over the ions.
    pub z2_over_z: f64,
    /// D_ion⁻² = (4πe²/kT) Σ n_s Z_s² [cm⁻²].
    pub d_ion_inv2: f64,
    /// k_TF² = 4πe² ∂n_e/∂μ_e [cm⁻²].
    pub k_tf2: f64,
    /// (kT/n_e) ∂n_e/∂μ_e: 1 for classical electrons, below 1 with degeneracy.
    pub degeneracy_factor: f64,
}

impl Plasma {
    /// The plasma at density `rho`, temperature `t`, full ionisation, with
    /// the electron degeneracy parameter `eta` of the EOS (NaN: classical
    /// electrons). `None` outside the physical range.
    pub fn new(rho: f64, t: f64, eta: f64, c: &Composition) -> Option<Self> {
        if !(rho > 0.0 && t > 0.0) {
            return None;
        }
        let ye = c.electrons_per_mu();
        let z2 = c.z2_per_mu();
        if !(ye > 0.0 && z2 > 0.0) {
            return None;
        }
        let kt = K_B * t;
        let n_e = rho * ye / M_U;
        let a_e = (3.0 / (FOUR_PI * n_e)).cbrt();
        let e2 = E_CHARGE * E_CHARGE;
        let degeneracy_factor = degeneracy_factor(eta, kt);
        Some(Self {
            kt,
            t13: t.cbrt(),
            a_e,
            gamma_e: e2 / (a_e * kt),
            z2_over_z: z2 / ye,
            d_ion_inv2: FOUR_PI * e2 / kt * rho * z2 / M_U,
            k_tf2: FOUR_PI * e2 * n_e * degeneracy_factor / kt,
            degeneracy_factor,
        })
    }

    /// Debye–Hückel screening length with electrons, D⁻² = D_ion⁻² + k_TF² [cm⁻²].
    pub fn d_inv2(&self) -> f64 {
        self.d_ion_inv2 + self.k_tf2
    }

    /// Weak-screening (Debye–Hückel) exponent Z_iZ_j e²/(kT D), with electrons.
    pub fn h_weak(&self, zi: f64, zj: f64) -> f64 {
        zi * zj * E_CHARGE * E_CHARGE / self.kt * self.d_inv2().sqrt()
    }

    /// Rigid-background exponent of Chugunov & DeWitt (2009, eq. A4) of a pair.
    pub fn h_ion_pair(&self, p: &Pair) -> f64 {
        let ge = self.gamma_e;
        let gij = p.g_ij * ge;
        let tij = match p.tau_t13 {
            None => 1.0,
            Some(tt) => {
                let z = 3.0 * gij * self.t13 / tt;
                (1.0 + 0.013 * p.y * p.y * z + 0.406 * p.y014 * z * z + (0.062 * p.y019 + 1.8 / gij) * z * z * z).cbrt()
            }
        };
        let g = ge / tij;
        let g2 = gij * gij;
        let cij = p.c_ij * self.z2_over_z.sqrt();
        (cij + g2) / (1.0 + g2) * (f0(p.z53[0] * g) + f0(p.z53[1] * g) - f0(p.z53[2] * g))
    }

    /// Electron-screening factor of Potekhin & Chabrier (2012, eqs. 20–21) of
    /// a pair, √(1 + k_TF² D̃_ion²) with D̃_ion² = D_ion² + (0.6 a_ij)².
    pub fn electron_factor_pair(&self, p: &Pair) -> f64 {
        let a_ij = p.a_ij * self.a_e;
        (1.0 + self.k_tf2 * (1.0 / self.d_ion_inv2 + 0.36 * a_ij * a_ij)).sqrt()
    }

    /// Screening exponent H_ij = ln f_ij of a pair, with electrons.
    pub fn h_pair(&self, p: &Pair) -> f64 {
        self.h_ion_pair(p) * self.electron_factor_pair(p)
    }

    /// Rigid-background exponent between charges `zi`, `zj`; with `mu_amu`
    /// (reduced mass [u]) the quantum correction t_ij is applied, `None` gives
    /// the static form (t_ij = 1).
    pub fn h_ion(&self, zi: f64, zj: f64, mu_amu: Option<f64>) -> f64 {
        self.h_ion_pair(&Pair::new(zi, zj, mu_amu))
    }

    /// Electron-screening factor between charges `zi`, `zj`.
    pub fn electron_factor(&self, zi: f64, zj: f64) -> f64 {
        self.electron_factor_pair(&Pair::new(zi, zj, None))
    }

    /// Screening exponent H_ij = ln f_ij of a non-resonant binary reaction
    /// between charges `zi`, `zj` with reduced mass `mu_amu` [u].
    pub fn h(&self, zi: f64, zj: f64, mu_amu: f64) -> f64 {
        self.h_pair(&Pair::new(zi, zj, Some(mu_amu)))
    }

    /// Static exponent (t_ij = 1) of a pair, with electrons.
    pub fn h_static(&self, zi: f64, zj: f64) -> f64 {
        self.h_pair(&Pair::new(zi, zj, None))
    }

    /// ln f_3α = H_{2,2} + H_{4,2}, both static: ⁸Be formation and the capture
    /// of the third α.
    pub fn h_triple_alpha(&self) -> f64 {
        self.h_static(2.0, 2.0) + self.h_static(4.0, 2.0)
    }
}

/// (kT/n_e) (∂n_e/∂μ_e)_T for electrons with degeneracy parameter `eta`
/// (classical, 1, if `eta` is not finite or very negative).
fn degeneracy_factor(eta: f64, kt: f64) -> f64 {
    if !eta.is_finite() || eta < -30.0 {
        return 1.0;
    }
    let theta = kt / (M_E * C * C);
    match electron_gas(eta, theta) {
        Ok(g) if g.n > 0.0 && g.dn_deta > 0.0 => g.dn_deta / g.n,
        _ => 1.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn helium() -> Composition {
        Composition { x: 0.0, y: 0.98, z: 0.02, x_cno: 0.014, x_c: 0.0, x_o: 0.0 }
    }

    #[test]
    fn free_energy_limits() {
        assert!((A3 - (-(3f64.sqrt()) / 2.0 - A1 / A2.sqrt())).abs() < 1e-15);
        for g in [1e-8_f64, 1e-6, 1e-4] {
            let dh = -g.powf(1.5) / 3f64.sqrt();
            assert!((f0(g) / dh - 1.0).abs() < 1e-3, "Γ = {g}: {} vs {dh}", f0(g));
        }
        // Strong coupling: −0.907Γ, less B₁Γ and the √Γ terms (1.1% at Γ = 10⁴).
        let g = 1e4;
        assert!((f0(g) / (A1 * g) - 1.0).abs() < 0.02, "f₀({g}) = {}", f0(g));
        // Continuous across the series switches.
        for g0 in [1e-4 * A2, 1e-8, 1e-6 * B2] {
            let (a, b) = (f0(g0 * (1.0 - 1e-9)), f0(g0 * (1.0 + 1e-9)));
            assert!((a - b).abs() < 1e-6 * a.abs().max(1e-300), "Γ = {g0}: {a} vs {b}");
        }
    }

    /// Weak coupling: H_ij → Z_iZ_j e²/(kTD); with classical electrons that is
    /// the Salpeter formula 0.188 Z_iZ_j ζ ρ^½ T₆^{−3/2}, ζ² = Σ X_i Z_i(Z_i + 1)/A_i.
    #[test]
    fn weak_limit_is_salpeter() {
        let c = Composition::new(0.7, 0.02, 0.7).unwrap();
        let (rho, t) = (1e-6, 1e8);
        let pl = Plasma::new(rho, t, f64::NAN, &c).unwrap();
        let zeta = (c.z2_per_mu() + c.electrons_per_mu()).sqrt();
        for (zi, zj, mu) in [(1.0, 1.0, 0.5), (1.0, 7.0, 0.93), (2.0, 2.0, 2.0), (2.0, 6.0, 3.0)] {
            let salpeter = 0.188 * zi * zj * zeta * rho.sqrt() * (t * 1e-6).powf(-1.5);
            let h = pl.h(zi, zj, mu);
            assert!((h / salpeter - 1.0).abs() < 2e-3, "{zi}+{zj}: H = {h:e}, Salpeter {salpeter:e}");
            assert!((h / pl.h_weak(zi, zj) - 1.0).abs() < 1e-3);
        }
        // The present (pre-9.7) formula used ζ² = 2X + 1.5Y + 4Z; same for Z = 0.
        let c0 = Composition::new(0.7, 0.0, 0.0).unwrap();
        let pl0 = Plasma::new(rho, t, f64::NAN, &c0).unwrap();
        let old = 0.188 * (2.0 * c0.x + 1.5 * c0.y).sqrt() * rho.sqrt() * (t * 1e-6).powf(-1.5);
        assert!((pl0.h(1.0, 1.0, 0.5) / old - 1.0).abs() < 2e-3);
    }

    /// Triple-α, weak: H_3α → 12 e²/(kTD).
    #[test]
    fn triple_alpha_weak_limit() {
        let pl = Plasma::new(1e-4, 1e8, f64::NAN, &helium()).unwrap();
        let h = pl.h_triple_alpha();
        let want = pl.h_weak(3.0, 4.0); // 12 e²/(kTD)
        assert!((h / want - 1.0).abs() < 2e-3, "{h:e} vs {want:e}");
    }

    /// Strong coupling, classical liquid: H_ij ≈ 0.9 Γ_e [Z_c^{5/3} − Z_i^{5/3} − Z_j^{5/3}]
    /// for the rigid background with t_ij = 1; the static triple-α telescopes
    /// to 3f₀(Γ₂) − f₀(Γ₆) there.
    #[test]
    fn strong_limit() {
        // Γ_e ≈ 300 for a cold, dense helium plasma.
        let pl = Plasma::new(1e9, 3e5, f64::NAN, &helium()).unwrap();
        assert!(pl.gamma_e > 100.0, "Γ_e = {}", pl.gamma_e);
        let p53 = |z: f64| z.powf(5.0 / 3.0);
        for (zi, zj) in [(1.0, 1.0), (2.0, 2.0), (2.0, 6.0), (4.0, 2.0)] {
            let want = 0.9 * pl.gamma_e * (p53(zi + zj) - p53(zi) - p53(zj));
            let h = pl.h_ion(zi, zj, None);
            assert!((h / want - 1.0).abs() < 0.03, "{zi}+{zj}: {h} vs {want}");
        }
        let h3 = pl.h_ion(2.0, 2.0, None) + pl.h_ion(4.0, 2.0, None);
        let g = |z: f64| p53(z) * pl.gamma_e;
        let tele = 3.0 * f0(g(2.0)) - f0(g(6.0));
        assert!((h3 / tele - 1.0).abs() < 0.01, "{h3} vs {tele}");
    }

    /// Dilute limit: all exponents → 0; ln f is smooth through intermediate
    /// coupling (no jumps in d ln H/d ln ρ), from 10⁻⁶ to 10⁷ g cm⁻³ at 10⁸ K.
    #[test]
    fn dilute_limit_and_smoothness() {
        let c = helium();
        let t = 1e8;
        let h = |lr: f64| Plasma::new(10f64.powf(lr), t, f64::NAN, &c).unwrap().h(2.0, 6.0, 3.0);
        assert!(h(-12.0) < 1e-6 && h(-12.0) > 0.0);
        let n = 600;
        let step = 13.0 / n as f64;
        let slope: Vec<f64> = (0..=n)
            .map(|i| {
                let lr = -6.0 + i as f64 * step;
                (h(lr + 1e-4).ln() - h(lr - 1e-4).ln()) / (2e-4 * std::f64::consts::LN_10)
            })
            .collect();
        for w in slope.windows(2) {
            assert!((w[1] - w[0]).abs() < 0.05, "slope jumps {} → {}", w[0], w[1]);
        }
        // Weak: slope ½ (H ∝ ρ^½).
        assert!((slope[0] - 0.5).abs() < 0.02, "{}", slope[0]);
    }

    /// The electron factor gives the Debye length with electrons in the weak
    /// limit, and degenerate electrons screen less than classical ones.
    #[test]
    fn electron_factor() {
        let c = helium();
        let pl = Plasma::new(1e-6, 1e8, f64::NAN, &c).unwrap();
        let want = (pl.d_inv2() / pl.d_ion_inv2).sqrt();
        assert!((pl.electron_factor(2.0, 2.0) / want - 1.0).abs() < 1e-4);
        let degenerate = Plasma::new(1e6, 1e8, 10.0, &c).unwrap();
        assert!(degenerate.degeneracy_factor < 0.2, "{}", degenerate.degeneracy_factor);
        let classical = Plasma::new(1e6, 1e8, f64::NAN, &c).unwrap();
        assert!(degenerate.h(2.0, 2.0, 2.0) < classical.h(2.0, 2.0, 2.0));
        // Non-degenerate η gives the classical factor.
        assert!((Plasma::new(1e-3, 1e8, -15.0, &c).unwrap().degeneracy_factor - 1.0).abs() < 1e-6);
    }
}
