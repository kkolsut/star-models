//! Fermi–Dirac electron EOS with ideal ions and radiation (plan 6.2), with
//! the ionisation equilibrium of H and He (6.3), the Pols-type
//! pressure-ionisation term (6.4) and molecular hydrogen (6.5); and the
//! Paczyński (1983) interpolation as a comparison model (6.5).
//! Equations: `docs/physics/02-equation-of-state.md`.
//!
//! *Structure.* Every quantity is an explicit function of the electron
//! degeneracy parameter η and T: the Fermi–Dirac integrals give n_e, P_e, u_e;
//! the equilibrium conditions give the ionisation (and dissociation) ratios at
//! that η; charge neutrality and conservation of nuclei then give ρ in closed
//! form. For given (P, T) the EOS solves P(η, T) = P for η (safeguarded
//! Newton, |Δη| < 10⁻¹³). The derivatives follow from central differences of
//! the explicit map in (η, ln T) and the chain rule (no further root solve).
//!
//! *Composition.* Ions per atomic mass unit: X (hydrogen), Y/4, Z/⟨A⟩ (the
//! metal mixture of [`crate::mixture`], GN93). Electrons of the metals: Z⟨Z/A⟩
//! per atomic mass unit; with ionisation they are treated as a pool of bound
//! electrons that ionises with the hydrogen energy and weights (an
//! approximation, stated: no metal ionisation stages are modelled, the metals
//! neither supply spurious free electrons in the atmosphere nor go missing in
//! the interior).
//!
//! *Approximations stated with 6.3:* ground-state statistical weights only,
//! ion mass ratios 1. H₂ (6.5): rigid rotor in the high-temperature limit and
//! harmonic oscillator, dissociation energy D₀ = 4.478 eV from v = 0. Pols
//! term (6.4): η in the term is the ideal Fermi–Dirac inversion (a departure
//! from Pols et al. 1995, eq. 29). Coulomb corrections and electron–positron
//! pairs are not included.

use super::eos::{ConsistencyScan, derive_state};
use super::fermi_dirac::{FdError, FdIntegrals, fermi_dirac};
use super::{Eos, EosState};
use crate::composition::Composition;
use crate::constants::{A_RAD, C, EV, H_PLANCK, K_B, M_E, M_H, M_U};
use crate::mixture::GN93;

/// Ionisation energies [eV] (NIST ASD, retrieved 2026-11-11).
pub const CHI_H: f64 = 13.598_434_599_702;
pub const CHI_HE1: f64 = 24.587_389_011;
pub const CHI_HE2: f64 = 54.417_765_528_2;
/// H₂: dissociation energy from v = 0 [eV], rotational and vibrational
/// temperatures [K] (Huber & Herzberg 1979, as in the NIST Chemistry WebBook:
/// B₀ = B_e − α_e/2 = 59.32 cm⁻¹ → Θ_r = 85.3 K, ω_e = 4401.2 cm⁻¹ → Θ_v = 6332 K,
/// D₀ = 36118 cm⁻¹).
pub const D0_H2: f64 = 4.4781;
pub const THETA_ROT_H2: f64 = 85.3;
pub const THETA_VIB_H2: f64 = 6332.0;
/// Pols et al. (1995) pressure-ionisation constants (c₁ [g cm⁻³], c₂, c₃, c₄ [g cm⁻³]) and energy scale.
pub const POLS_C: (f64, f64, f64, f64) = (3.0, 0.25, 2.0, 0.03);
pub const POLS_CHI: f64 = 13.60;

/// Prefactor of n_e: 8π√2 (m_e c / h)³.
fn k_n() -> f64 {
    8.0 * std::f64::consts::PI * std::f64::consts::SQRT_2 * (M_E * C / H_PLANCK).powi(3)
}

/// Electron gas at (η, ϑ): n_e, P_e, ρu_e (energy per volume, no rest mass),
/// ∂n_e/∂η, ∂n_e/∂ϑ.
#[derive(Clone, Copy, Debug)]
pub struct ElectronGas {
    pub n: f64,
    pub p: f64,
    pub e: f64,
    pub dn_deta: f64,
    pub dn_dtheta: f64,
}

pub fn electron_gas(eta: f64, theta: f64) -> Result<ElectronGas, FdError> {
    let FdIntegrals { f, d_eta, d_theta } = fermi_dirac(eta, theta)?;
    let kn = k_n();
    let mc2 = M_E * C * C;
    let t15 = theta.powf(1.5);
    let t25 = t15 * theta;
    let n = kn * t15 * (f[0] + theta * f[1]);
    let p = 2.0 / 3.0 * kn * mc2 * t25 * (f[1] + 0.5 * theta * f[2]);
    let e = kn * mc2 * t25 * (f[1] + theta * f[2]);
    let dn_deta = kn * t15 * (d_eta[0] + theta * d_eta[1]);
    let dn_dtheta = kn * (1.5 * theta.sqrt() * (f[0] + theta * f[1]) + t15 * (d_theta[0] + f[1] + theta * d_theta[1]));
    Ok(ElectronGas { n, p, e, dn_deta, dn_dtheta })
}

/// Ideal Fermi–Dirac inversion: η with n_e(η, ϑ) = n (Newton in η on ln n,
/// which is increasing and nearly linear for η ≪ 0, to |ln n_e − ln n| < 10⁻¹⁴),
/// starting from η = `guess` whose electron gas is `start` if known.
pub fn invert_density_from(
    n: f64,
    theta: f64,
    guess: f64,
    start: Option<ElectronGas>,
) -> Result<(f64, ElectronGas), FdError> {
    let mut eta = guess;
    let mut gas = match start {
        Some(g) => g,
        None => electron_gas(eta, theta)?,
    };
    for _ in 0..200 {
        let f = gas.n.ln() - n.ln();
        if f.abs() < 1e-14 {
            break;
        }
        let df = gas.dn_deta / gas.n;
        let step = (-f / df).clamp(-20.0, 20.0);
        eta += step;
        gas = electron_gas(eta, theta)?;
        if step.abs() < 1e-15 * eta.abs().max(1.0) {
            break;
        }
    }
    Ok((eta, gas))
}

/// [`invert_density_from`] without a known starting gas.
pub fn invert_density(n: f64, theta: f64, guess: f64) -> Result<(f64, ElectronGas), FdError> {
    invert_density_from(n, theta, guess, None)
}

/// Starting value of η for electron density `n` at ϑ (non-degenerate or
/// degenerate estimate).
fn eta_guess(n: f64, t: f64, theta: f64) -> f64 {
    let lambda3 = (H_PLANCK * H_PLANCK / (2.0 * std::f64::consts::PI * M_E * K_B * t)).powf(1.5);
    let nd = (0.5 * n * lambda3).ln();
    if nd < 0.0 {
        nd
    } else {
        let xf = H_PLANCK / (M_E * C) * (3.0 * n / (8.0 * std::f64::consts::PI)).cbrt();
        ((1.0 + xf * xf).sqrt() - 1.0) / theta
    }
}

const MEMO_SIZE: usize = 8;

/// Last converged η: (flags, ln P, ln T, X, Z, η).
type LastEta = Option<(u8, f64, f64, f64, f64, f64)>;

thread_local! {
    /// Last converged η: the starting value of the next root search when the
    /// state is close (purely a speed-up; the result is converged to the same
    /// tolerance either way).
    static LAST_ETA: std::cell::Cell<LastEta> = const { std::cell::Cell::new(None) };
    /// The last few states (exact repeats of (P, T, X, Y, Z)).
    static MEMO: std::cell::RefCell<Vec<([u64; 6], EosState)>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// Fermi–Dirac EOS (see the module documentation).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FermiDiracEos {
    /// H/He ionisation equilibrium; false = full ionisation (plan 6.2).
    pub ionisation: bool,
    /// Pols-type pressure ionisation (needs `ionisation`).
    pub pressure_ionisation: bool,
    /// Molecular hydrogen (needs `ionisation`).
    pub molecules: bool,
}

/// Everything at one (η, T).
#[derive(Clone, Copy, Debug, Default)]
pub struct Point {
    pub ln_p: f64,
    pub ln_rho: f64,
    pub h: f64,
    pub mu: f64,
    pub x_h: f64,
    pub x_he1: f64,
    pub x_he2: f64,
    pub x_h2: f64,
    /// Pressure-ionisation shift of the electron chemical potential, μ_PI/kT.
    pub psi: f64,
    pub p_pi: f64,
}

/// The Pols term: Q(n, T) = n g(n, T) and its partial derivatives.
struct PolsQ {
    q: f64,
    qn: f64,
    qt: f64,
}

/// Q, ∂Q/∂n, ∂Q/∂T at free-electron density `n` with its ideal η and gas.
fn pols_q(n: f64, t: f64, theta: f64, eta: f64, gas: &ElectronGas) -> PolsQ {
    let (c1, c2, c3, c4) = POLS_C;
    let x = M_H * n;
    let a = (c1 / x).powf(c2);
    let e = (-a).exp();
    let kt = K_B * t;
    let b = POLS_CHI * EV / kt + eta + c3 * (x / c4).ln_1p();
    let g = e * b;
    // ∂η/∂n and ∂η/∂T at fixed n.
    let eta_n = 1.0 / gas.dn_deta;
    let eta_t = -gas.dn_dtheta * (theta / t) / gas.dn_deta;
    let de_dn = e * c2 * a / n;
    let db_dn = eta_n + c3 * (M_H / c4) / (1.0 + x / c4);
    let gn = de_dn * b + e * db_dn;
    let gt = e * (-POLS_CHI * EV / (kt * t) + eta_t);
    PolsQ { q: n * g, qn: g + n * gn, qt: n * gt }
}

impl FermiDiracEos {
    /// Fermi–Dirac electrons with full ionisation (plan 6.2).
    pub fn fully_ionised() -> Self {
        Self { ionisation: false, pressure_ionisation: false, molecules: false }
    }

    /// With the ionisation equilibrium (6.3), optionally pressure ionisation
    /// (6.4) and H₂ (6.5).
    pub fn ionising(pressure_ionisation: bool, molecules: bool) -> Self {
        Self { ionisation: true, pressure_ionisation, molecules }
    }

    /// The explicit map (η, ln T) → state.
    pub fn point(&self, eta: f64, ln_t: f64, c: &Composition) -> Result<Point, FdError> {
        let t = ln_t.exp();
        let kt = K_B * t;
        let theta = kt / (M_E * C * C);
        let gas = electron_gas(eta, theta)?;
        let n_e = gas.n;
        let za = GN93.z_over_a;
        let p_rad = A_RAD * t.powi(4) / 3.0;

        // Pressure-ionisation shift ψ = μ_PI/kT = −Q_n(n_e).
        let psi = if self.pressure_ionisation { -pols_q(n_e, t, theta, eta, &gas).qn } else { 0.0 };

        // Ionisation ratios (log space), equilibrium with the full η.
        let (x_h, x_he1, x_he2, ln_r_h) = if self.ionisation {
            let ln_r_h = 0.5f64.ln() - eta - CHI_H * EV / kt - psi;
            let z1 = 2.0f64.ln() - eta - CHI_HE1 * EV / kt - psi;
            let z2 = z1 + 0.5f64.ln() - eta - CHI_HE2 * EV / kt - psi;
            let m = 0.0f64.max(z1).max(z2);
            let (e0, e1, e2) = ((-m).exp(), (z1 - m).exp(), (z2 - m).exp());
            let s = e0 + e1 + e2;
            (1.0 / (1.0 + (-ln_r_h).exp()), e1 / s, e2 / s, ln_r_h)
        } else {
            (1.0, 0.0, 1.0, f64::INFINITY)
        };
        let x_pool = x_h;
        // Free electrons per m_u from He and the metal pool, b·m_u.
        let nu_b = c.y / 4.0 * (x_he1 + 2.0 * x_he2) + c.z * za * x_pool;

        // Hydrogen: atoms n0, protons r_H n0, molecules n0²/K_D.
        let (rho, n_h0, n_hp, n_h2) = if self.ionisation && self.molecules && c.x > 0.0 {
            let r_h = ln_r_h.min(700.0).exp();
            let ln_kd = 1.5 * (std::f64::consts::PI * M_H * kt / (H_PLANCK * H_PLANCK)).ln() + 4f64.ln()
                - ((t / (2.0 * THETA_ROT_H2)).ln())
                + (-(-THETA_VIB_H2 / t).exp()).ln_1p()
                - D0_H2 * EV / kt;
            let inv_kd = (-ln_kd).min(700.0).exp();
            let b = nu_b / M_U;
            let a2 = 2.0 * b * inv_kd;
            let b2 = (1.0 + r_h) * b + c.x * r_h / M_U;
            let c2 = c.x * n_e / M_U;
            let n0 = 2.0 * c2 / (b2 + (b2 * b2 + 4.0 * a2 * c2).sqrt());
            let rho = M_U / c.x * (2.0 * n0 * n0 * inv_kd + (1.0 + r_h) * n0);
            (rho, n0, r_h * n0, n0 * n0 * inv_kd)
        } else {
            let nu = c.x * x_h + nu_b;
            let rho = n_e * M_U / nu;
            let n_h = c.x * rho / M_U;
            (rho, n_h * (1.0 - x_h), n_h * x_h, 0.0)
        };
        if !(rho > 0.0 && rho.is_finite()) {
            return Err(FdError::OutOfDomain { eta, theta });
        }
        let n_he = c.y * rho / (4.0 * M_U);
        let n_z = c.z * rho / (GN93.mean_a * M_U);
        let n_heavy = n_h0 + n_hp + n_h2 + n_he + n_z;

        // Pressure-ionisation pressure and energy.
        let (p_pi, u_pi) = if self.pressure_ionisation {
            let n_tot = rho / M_U * (c.x + c.y / 2.0 + c.z * za);
            let (eta0, gas0) = invert_density_from(n_tot, theta, eta, Some(gas))?;
            let q0 = pols_q(n_tot, t, theta, eta0, &gas0);
            let q = pols_q(n_e, t, theta, eta, &gas);
            (kt * (n_tot * q0.qn - q0.q - n_e * q.qn + q.q), -kt * t / rho * (q0.qt - q.qt))
        } else {
            (0.0, 0.0)
        };

        let p = n_heavy * kt + gas.p + p_rad + p_pi;
        if !(p > 0.0) {
            return Err(FdError::OutOfDomain { eta, theta });
        }
        // Energy per volume; zero: neutral atoms (with ionisation) or fully
        // ionised gas (without).
        let mut e = 1.5 * kt * n_heavy + gas.e + A_RAD * t.powi(4);
        if self.ionisation {
            e += CHI_H * EV * (n_hp + c.z * za * rho / M_U * x_pool)
                + n_he * (CHI_HE1 * EV * x_he1 + (CHI_HE1 + CHI_HE2) * EV * x_he2);
            if n_h2 > 0.0 {
                let vib = K_B * THETA_VIB_H2 / ((THETA_VIB_H2 / t).exp() - 1.0);
                e += n_h2 * (kt + vib - D0_H2 * EV);
            }
        }
        let u = e / rho + u_pi;
        let h = u + p / rho;
        let mu = rho / (M_U * (n_heavy + n_e));
        let n_hnuc = n_h0 + n_hp + 2.0 * n_h2;
        let x_h2 = if n_hnuc > 0.0 { 2.0 * n_h2 / n_hnuc } else { 0.0 };
        Ok(Point { ln_p: p.ln(), ln_rho: rho.ln(), h, mu, x_h, x_he1, x_he2, x_h2, psi, p_pi })
    }

    /// Solve P(η, T) = P for η.
    pub fn solve_eta(&self, p: f64, t: f64, c: &Composition) -> Result<(f64, Point), String> {
        let p_rad = A_RAD * t.powi(4) / 3.0;
        if !(t > 0.0 && p > p_rad) {
            return Err(format!("P = {p:e} not above the radiation pressure {p_rad:e} at T = {t:e}"));
        }
        let lt = t.ln();
        let lp = p.ln();
        let theta = K_B * t / (M_E * C * C);
        // Start from the fully ionised ideal gas.
        let mu0 = 1.0 / (2.0 * c.x + 0.75 * c.y + c.z * (GN93.z_over_a + 1.0 / GN93.mean_a));
        let rho0 = (p - p_rad) * mu0 * M_U / (K_B * t);
        let n0 = rho0 / M_U * (c.x + c.y / 2.0 + c.z * GN93.z_over_a);
        let flags = self.flags();
        let warm = LAST_ETA.with(|cell| {
            cell.get().and_then(|(f, lp0, lt0, x0, z0, eta0)| {
                ((f == flags)
                    && (lp - lp0).abs() < 0.1
                    && (lt - lt0).abs() < 0.1
                    && (c.x - x0).abs() < 0.02
                    && (c.z - z0).abs() < 0.01)
                    .then_some(eta0)
            })
        });
        let mut eta = warm.unwrap_or_else(|| eta_guess(n0, t, theta).clamp(-90.0, 5e3));
        let f = |eta: f64| -> Result<(f64, Point), String> {
            let pt = self.point(eta, lt, c).map_err(|e| e.to_string())?;
            Ok((pt.ln_p - lp, pt))
        };
        let (mut lo, mut hi): (Option<f64>, Option<f64>) = (None, None);
        let (mut r, mut pt) = f(eta)?;
        // Previous iterate for the secant slope.
        let mut prev: Option<(f64, f64)> = None;
        for _ in 0..200 {
            if r == 0.0 {
                break;
            }
            if r < 0.0 {
                lo = Some(eta);
            } else {
                hi = Some(eta);
            }
            let d = match prev {
                Some((e0, r0)) if e0 != eta && r0 != r => (r - r0) / (eta - e0),
                _ => {
                    let h = 1e-6 * eta.abs().max(1.0);
                    let (rp, _) = f(eta + h)?;
                    (rp - r) / h
                }
            };
            let mut next = if d > 0.0 { eta - r / d } else { eta + if r < 0.0 { 5.0 } else { -5.0 } };
            // Keep inside the bracket; bisect if the step leaves it.
            if let (Some(a), Some(b)) = (lo, hi) {
                if !(next > a.min(b) && next < a.max(b)) {
                    next = 0.5 * (a + b);
                }
            } else {
                next = next.clamp(eta - 20.0, eta + 20.0 * eta.abs().max(1.0));
            }
            let step = next - eta;
            prev = Some((eta, r));
            eta = next;
            (r, pt) = f(eta)?;
            if step.abs() < 1e-13 * eta.abs().max(1.0) {
                break;
            }
        }
        if !(r.abs() < 1e-11) {
            return Err(format!("η iteration did not converge at P = {p:e}, T = {t:e}"));
        }
        LAST_ETA.with(|cell| cell.set(Some((flags, lp, lt, c.x, c.z, eta))));
        Ok((eta, pt))
    }

    fn flags(&self) -> u8 {
        self.ionisation as u8 | (self.pressure_ionisation as u8) << 1 | (self.molecules as u8) << 2
    }

    /// Full state with derivatives at (P, T).
    pub fn state(&self, p: f64, t: f64, c: &Composition) -> Result<EosState, String> {
        let (eta, pt) = self.solve_eta(p, t, c)?;
        let lt = t.ln();
        let he = 1e-5 * eta.abs().max(1.0);
        let ht = 1e-5;
        let e = |r: Result<Point, FdError>| r.map_err(|e| e.to_string());
        let ep = e(self.point(eta + he, lt, c))?;
        let em = e(self.point(eta - he, lt, c))?;
        let tp = e(self.point(eta, lt + ht, c))?;
        let tm = e(self.point(eta, lt - ht, c))?;
        let d = |a: f64, b: f64, h: f64| (a - b) / (2.0 * h);
        let (p_e, r_e, h_e) = (d(ep.ln_p, em.ln_p, he), d(ep.ln_rho, em.ln_rho, he), d(ep.h, em.h, he));
        let (p_t, r_t, h_t) = (d(tp.ln_p, tm.ln_p, ht), d(tp.ln_rho, tm.ln_rho, ht), d(tp.h, tm.h, ht));
        if !(p_e > 0.0) {
            return Err(format!("P not increasing with η at P = {p:e}, T = {t:e}"));
        }
        let alpha = r_e / p_e;
        let delta = -(r_t - r_e * p_t / p_e);
        let cp = (h_t - h_e * p_t / p_e) / t;
        let rho = pt.ln_rho.exp();
        let mut s = EosState::from_derivatives(p, t, rho, pt.h, alpha, delta, cp, A_RAD * t.powi(4) / 3.0, pt.mu);
        s.eta = eta;
        s.x_h = pt.x_h;
        s.x_he1 = pt.x_he1;
        s.x_he2 = pt.x_he2;
        s.x_h2 = pt.x_h2;
        Ok(s)
    }
}

impl Eos for FermiDiracEos {
    fn name(&self) -> &'static str {
        match (self.ionisation, self.pressure_ionisation, self.molecules) {
            (false, _, _) => "Fermi–Dirac electrons + ideal ions + radiation (fully ionised)",
            (true, false, false) => "Fermi–Dirac electrons, H/He ionisation (no pressure ionisation)",
            (true, true, false) => "Fermi–Dirac electrons, H/He ionisation + Pols pressure ionisation",
            (true, false, true) => "Fermi–Dirac electrons, H/He ionisation + H₂",
            (true, true, true) => "Fermi–Dirac electrons, H/He ionisation + Pols pressure ionisation + H₂",
        }
    }

    fn eval(&self, p: f64, t: f64, comp: &Composition) -> EosState {
        self.try_eval(p, t, comp).unwrap_or_else(|_| EosState::invalid())
    }

    fn degenerate_electrons(&self) -> bool {
        true
    }

    fn try_eval(&self, p: f64, t: f64, comp: &Composition) -> Result<EosState, String> {
        // The finite-difference Jacobian repeats (P, T) whenever it perturbs
        // r or L: remember the last few states.
        let key = [self.flags() as u64, p.to_bits(), t.to_bits(), comp.x.to_bits(), comp.y.to_bits(), comp.z.to_bits()];
        if let Some(s) = MEMO.with(|m| m.borrow().iter().find(|(k, _)| *k == key).map(|(_, s)| *s)) {
            return Ok(s);
        }
        let s = self.state(p, t, comp)?;
        MEMO.with(|m| {
            let mut m = m.borrow_mut();
            if m.len() >= MEMO_SIZE {
                m.remove(0);
            }
            m.push((key, s));
        });
        Ok(s)
    }
}

/// Paczyński (1983, eq. 4) interpolation between non-degenerate, degenerate
/// non-relativistic and degenerate relativistic electrons, with ideal ions
/// and radiation; u_e = P_e/((f − 1)ρ). Comparison model only: not derived
/// from a free energy (see `docs/physics/02`). Published constants
/// K₁ = 9.91×10¹², K₂ = 1.231×10¹⁵ (cgs, ρ/μ_e in g cm⁻³).
#[derive(Clone, Copy, Debug, Default)]
pub struct Paczynski;

impl Paczynski {
    pub const K1: f64 = 9.91e12;
    pub const K2: f64 = 1.231e15;

    fn mu(c: &Composition) -> (f64, f64) {
        let mu_i = 1.0 / (c.x + c.y / 4.0 + c.z / GN93.mean_a);
        let mu_e = 1.0 / (c.x + c.y / 2.0 + c.z * GN93.z_over_a);
        (mu_i, mu_e)
    }

    /// (P, u) at (ρ, T).
    pub fn p_u(rho: f64, t: f64, c: &Composition) -> (f64, f64) {
        let (mu_i, mu_e) = Self::mu(c);
        let p_ion = rho * K_B * t / (mu_i * M_U);
        let p_nd = rho * K_B * t / (mu_e * M_U);
        let r = rho / mu_e;
        let p_nr = Self::K1 * r.powf(5.0 / 3.0);
        let p_r = Self::K2 * r.powf(4.0 / 3.0);
        let p_d = (p_nr.powi(-2) + p_r.powi(-2)).powf(-0.5);
        let p_e = (p_nd * p_nd + p_d * p_d).sqrt();
        let f = 5.0 / 3.0 * (p_d / p_nr).powi(2) + 4.0 / 3.0 * (p_d / p_r).powi(2);
        let p_rad = A_RAD * t.powi(4) / 3.0;
        let u = 1.5 * p_ion / rho + p_e / ((f - 1.0) * rho) + 3.0 * p_rad / rho;
        (p_ion + p_e + p_rad, u)
    }

    /// ρ at (P, T): bisection in ln ρ (P increases with ρ), to 10⁻¹⁴.
    pub fn rho(p: f64, t: f64, c: &Composition) -> Option<f64> {
        if !(p > A_RAD * t.powi(4) / 3.0) {
            return None;
        }
        let (mut lo, mut hi) = (-40.0_f64, 25.0_f64);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if Self::p_u(mid.exp(), t, c).0 < p {
                lo = mid;
            } else {
                hi = mid;
            }
            if hi - lo < 1e-14 {
                break;
            }
        }
        Some((0.5 * (lo + hi)).exp())
    }
}

impl Eos for Paczynski {
    fn name(&self) -> &'static str {
        "Paczyński (1983) interpolation (comparison only)"
    }

    fn eval(&self, p: f64, t: f64, c: &Composition) -> EosState {
        let (mu_i, mu_e) = Self::mu(c);
        let mu = 1.0 / (1.0 / mu_i + 1.0 / mu_e);
        let f = |lp: f64, lt: f64| {
            let (p, t) = (lp.exp(), lt.exp());
            let rho = Self::rho(p, t, c)?;
            let (_, u) = Self::p_u(rho, t, c);
            Some((rho.ln(), u + p / rho, mu))
        };
        derive_state(&f, p, t).unwrap_or_else(EosState::invalid)
    }
}

/// Consistency scan of an EOS over the stellar range used in the tests and
/// recorded in `docs/physics/02`: ln P from 10² to 10¹⁹ dyn cm⁻², T from 3000 K
/// to 10⁸ K.
pub fn stellar_scan(eos: &dyn Eos, comp: &Composition, n: usize) -> ConsistencyScan {
    super::eos::consistency_scan(
        eos,
        comp,
        (2.0 * std::f64::consts::LN_10, 19.0 * std::f64::consts::LN_10),
        (3000f64.ln(), 1e8f64.ln()),
        n,
        1e-6,
        1e-4,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::R_GAS;
    use crate::physics::eos::IdealGasRadiation;

    fn sun() -> Composition {
        Composition::new(0.7, 0.02, GN93.f_cno).unwrap()
    }

    /// Ideal gas with the same ions and electrons as the Fermi–Dirac EOS.
    fn ideal_rho(p: f64, t: f64, c: &Composition) -> f64 {
        let inv_mu = c.x + c.y / 4.0 + c.z / GN93.mean_a + c.x + c.y / 2.0 + c.z * GN93.z_over_a;
        (p - A_RAD * t.powi(4) / 3.0) / (inv_mu * R_GAS * t)
    }

    /// Plan 6.2: the ideal EOS is recovered where the electrons are
    /// non-degenerate (n_e λ³/2 = e^η ≪ 1) and non-relativistic (ϑ = kT/m_ec² ≪ 1):
    /// to 10⁻⁶ where both are below 10⁻⁶, otherwise to their size.
    #[test]
    fn ideal_limit() {
        let c = sun();
        let eos = FermiDiracEos::fully_ionised();
        for (p, t) in [(1e3, 3e3), (1e2, 5e3), (1e4, 1e4), (1e7, 1e5), (3e9, 3e5)] {
            let s = eos.state(p, t, &c).unwrap();
            let theta = K_B * t / (M_E * C * C);
            let tol = 1e-6 + s.eta.exp() + 3.0 * theta;
            let want = ideal_rho(p, t, &c);
            assert!((s.rho / want - 1.0).abs() < tol, "ρ at ({p:e}, {t:e}), η = {}: {} vs {want}", s.eta, s.rho);
            // Derived quantities against the closed forms of the ideal EOS
            // (radiation included), which depend on β only.
            let id = IdealGasRadiation.eval(p, t, &c);
            for (name, a, b) in
                [("∇_ad", s.nabla_ad, id.nabla_ad), ("Γ₁", s.gamma1, id.gamma1), ("δ", s.delta, id.delta)]
            {
                assert!((a / b - 1.0).abs() < tol + 1e-6, "{name} at ({p:e}, {t:e}): {a} vs {b}");
            }
        }
    }

    /// The solar centre: compare with the ideal gas and report the departure
    /// (η ≈ −1.5 there: degeneracy pressure lowers ρ at given P by about 2%).
    #[test]
    fn solar_centre_departure() {
        let c = Composition::new(0.35, 0.02, GN93.f_cno).unwrap();
        let s = FermiDiracEos::fully_ionised().state(2.3e17, 1.55e7, &c).unwrap();
        let d = s.rho / ideal_rho(2.3e17, 1.55e7, &c) - 1.0;
        eprintln!("solar centre: η = {:.3}, ρ_FD/ρ_ideal − 1 = {d:.4}", s.eta);
        assert!(s.eta > -3.0 && s.eta < 0.0);
        assert!(d < 0.0 && d > -0.05);
    }

    /// Completely degenerate electron pressure (Chandrasekhar):
    /// P = (π m_e⁴c⁵/3h³)[x(2x² − 3)√(1 + x²) + 3 asinh x], x = p_F/m_e c.
    fn chandrasekhar(n_e: f64) -> f64 {
        let x = H_PLANCK / (M_E * C) * (3.0 * n_e / (8.0 * std::f64::consts::PI)).cbrt();
        std::f64::consts::PI * M_E.powi(4) * C.powi(5) / (3.0 * H_PLANCK.powi(3))
            * (x * (2.0 * x * x - 3.0) * (1.0 + x * x).sqrt() + 3.0 * x.asinh())
    }

    /// Degenerate limits: cold gas against the Chandrasekhar formula, and its
    /// limits P_e = K₁(ρ/μ_e)^{5/3}, K₂(ρ/μ_e)^{4/3} with the m_u constants
    /// 1.0036×10¹³ and 1.2435×10¹⁵.
    #[test]
    fn degenerate_limits() {
        // Pure helium: μ_e = 2.
        for (rho_mu, t) in [(10.0, 1e4), (1e4, 1e5), (1e6, 1e6), (5e8, 1e7)] {
            let n_e = rho_mu / M_U;
            let theta = K_B * t / (M_E * C * C);
            let (eta, gas) = invert_density(n_e, theta, eta_guess(n_e, t, theta)).unwrap();
            assert!(eta > 30.0, "η = {eta}");
            let want = chandrasekhar(n_e);
            // Finite-temperature correction ∝ (kT/E_F)².
            assert!((gas.p / want - 1.0).abs() < 2e-3, "ρ/μ_e = {rho_mu}: {} vs {want}", gas.p);
        }
        let k1 = chandrasekhar(10.0 / M_U) / 10f64.powf(5.0 / 3.0);
        assert!((k1 / 1.0036e13 - 1.0).abs() < 1e-3, "K₁ = {k1:e}");
        let r: f64 = 1e14;
        let k2 = chandrasekhar(r / M_U) / r.powf(4.0 / 3.0);
        assert!((k2 / 1.2435e15 - 1.0).abs() < 1e-3, "K₂ = {k2:e}");
    }

    #[test]
    fn eta_root_and_round_trip() {
        let c = sun();
        let eos = FermiDiracEos::fully_ionised();
        for (p, t) in [(1e17, 1.5e7), (1e22, 1e8), (1e8, 1e5)] {
            let (eta, pt) = eos.solve_eta(p, t, &c).unwrap();
            assert!((pt.ln_p - p.ln()).abs() < 1e-12, "η = {eta}");
        }
        // P ≤ P_rad is an error, not a clamp.
        let t: f64 = 1e7;
        assert!(eos.state(0.5 * A_RAD * t.powi(4) / 3.0, t, &c).is_err());
        assert!(eos.eval(0.5 * A_RAD * t.powi(4) / 3.0, t, &c).rho.is_nan());
    }

    /// Thermodynamic consistency of every variant (plan 6.1 test).
    #[test]
    fn consistency_scans() {
        let c = sun();
        for eos in [
            FermiDiracEos::fully_ionised(),
            FermiDiracEos::ionising(false, false),
            FermiDiracEos::ionising(true, false),
            FermiDiracEos::ionising(true, true),
        ] {
            let s = stellar_scan(&eos, &c, 12);
            eprintln!("{}: {s:?}", eos.name());
            assert!(s.passed(), "{}: {s:?}", eos.name());
        }
    }

    /// Plan 6.3: Saha limit, full ionisation at high T, ∇_ad minimum in the
    /// hydrogen ionisation zone, limits of μ.
    #[test]
    fn ionisation() {
        let c = sun();
        let eos = FermiDiracEos::ionising(false, false);
        // Saha for hydrogen at a dilute state: n_H⁺ n_e / n_H⁰ = (2πm_e kT/h²)^{3/2} e^{−χ/kT}.
        let (p, t) = (1e2, 9000.0);
        let s = eos.state(p, t, &c).unwrap();
        assert!(s.x_h > 0.01 && s.x_h < 0.99, "x_H = {}", s.x_h);
        let theta = K_B * t / (M_E * C * C);
        let gas = electron_gas(s.eta, theta).unwrap();
        let saha = (2.0 * std::f64::consts::PI * M_E * K_B * t / (H_PLANCK * H_PLANCK)).powf(1.5)
            * (-CHI_H * EV / (K_B * t)).exp();
        let lhs = s.x_h / (1.0 - s.x_h) * gas.n;
        assert!((lhs / saha - 1.0).abs() < 1e-3, "Saha: {lhs:e} vs {saha:e}");
        // Nearly fully ionised at high T and low density reproduces 6.2, to
        // the size of the remaining He⁺ fraction.
        let (p, t) = (1e14, 1e7);
        let a = eos.state(p, t, &c).unwrap();
        let b = FermiDiracEos::fully_ionised().state(p, t, &c).unwrap();
        let rest = 1.0 - a.x_he2;
        assert!(rest < 1e-3, "1 − x_He2 = {rest}");
        assert!((a.rho / b.rho - 1.0).abs() < 3.0 * rest + 1e-6, "{} vs {}", a.rho, b.rho);
        assert!((a.nabla_ad / b.nabla_ad - 1.0).abs() < 30.0 * rest + 1e-6);
        // ∇_ad has a minimum between 0.1 and 0.2 near 10⁴ K at photospheric densities.
        let min = (0..60)
            .map(|k| 10f64.powf(3.7 + 0.01 * k as f64))
            .map(|t| eos.state(1e5, t, &c).unwrap().nabla_ad)
            .fold(1.0, f64::min);
        assert!(min > 0.05 && min < 0.2, "∇_ad,min = {min}");
        // Neutral and ionised limits of μ.
        let neutral = eos.state(1e4, 3000.0, &c).unwrap();
        let want = 1.0 / (c.x + c.y / 4.0 + c.z / GN93.mean_a);
        assert!((neutral.mu / want - 1.0).abs() < 1e-3, "μ = {} vs {want}", neutral.mu);
        let ionised = eos.state(1e14, 1e7, &c).unwrap();
        let want = 1.0 / (2.0 * c.x + 0.75 * c.y + c.z * (GN93.z_over_a + 1.0 / GN93.mean_a));
        assert!((ionised.mu / want - 1.0).abs() < 1e-3, "μ = {} vs {want}", ionised.mu);
    }

    /// Plan 6.4: dilute gas recovers 6.3; no recombination in a dense
    /// interior; n = n₀ gives a vanishing term.
    #[test]
    fn pressure_ionisation() {
        let c = sun();
        let pols = FermiDiracEos::ionising(true, false);
        let ideal = FermiDiracEos::ionising(false, false);
        let (a, b) = (pols.state(1e3, 8000.0, &c).unwrap(), ideal.state(1e3, 8000.0, &c).unwrap());
        assert!((a.rho / b.rho - 1.0).abs() < 1e-6 && (a.x_h - b.x_h).abs() < 1e-6);
        // Dense, comparatively cool interior (ρ ≈ 10³ g/cm³, T = 10⁶ K): ideal
        // equilibrium recombines hydrogen, the Pols term keeps it ionised.
        let (p, t) = (5e17, 1e6);
        let (a, b) = (pols.state(p, t, &c).unwrap(), ideal.state(p, t, &c).unwrap());
        assert!(a.x_h > 0.99 && a.x_he2 > 0.99, "Pols: x_H = {}, x_He2 = {}", a.x_h, a.x_he2);
        assert!(b.x_h < a.x_h, "ideal: x_H = {}", b.x_h);
        // Full ionisation (n → n₀): P_PI → 0.
        let (eta, pt) = pols.solve_eta(1e18, 1e8, &c).unwrap();
        assert!(pt.p_pi.abs() < 1e-5 * 1e18, "P_PI = {}, η = {eta}", pt.p_pi);
    }

    /// Plan 6.5: molecular limit of μ at low T.
    #[test]
    fn molecular_hydrogen() {
        let c = sun();
        let eos = FermiDiracEos::ionising(true, true);
        let s = eos.state(1e6, 1500.0, &c).unwrap();
        assert!(s.x_h2 > 0.99, "x_H₂ = {}", s.x_h2);
        let want = 1.0 / (c.x / 2.0 + c.y / 4.0 + c.z / GN93.mean_a);
        assert!((s.mu / want - 1.0).abs() < 1e-2, "μ = {} vs {want}", s.mu);
        // Dissociated at 10⁴ K.
        assert!(eos.state(1e4, 1e4, &c).unwrap().x_h2 < 1e-3);
    }

    #[test]
    fn paczynski_limits() {
        let c = sun();
        let s = Paczynski.eval(1e15, 1e7, &c);
        let id = FermiDiracEos::fully_ionised().state(1e15, 1e7, &c).unwrap();
        assert!((s.rho / id.rho - 1.0).abs() < 1e-3, "{} vs {}", s.rho, id.rho);
    }
}
