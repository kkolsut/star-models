//! Outer boundary condition: the atmosphere (`docs/physics/06`).
//!
//! Two kinds, selected by [`Atmosphere`]:
//!
//! * [`Atmosphere::GreyClosedForm`] (the default): the grey Eddington
//!   photosphere with one value of κ, T⁴(τ) = ¾T_eff⁴(τ + ⅔), the last mesh
//!   point at τ = ⅔ where T = T_eff and
//!
//!   P_ph = ⅔ g / κ + a T_eff⁴ / 6.
//!
//! * [`Atmosphere::Integrated`] (plans 3.1–3.3): hydrostatic equilibrium
//!   dP/dτ = g/κ(ρ(P,T), T) integrated with the local opacity from τ₀ to the
//!   fitting point τ_fit by a fixed number of RK4 steps in ln τ, with a
//!   T(τ) relation T⁴ = ¾T_eff⁴[τ + q(τ)] ([`TTau`]: grey Eddington or the
//!   solar relation of Krishna Swamy 1966). With convection on, the pair
//!   (P, T) is integrated and ∇ comes from the mixing-length cubic with the
//!   Hopf factor A = 1 + q'(τ) and optically thin elements
//!   ([`super::convection`]).
//!
//! *Start (grey, LTE approximation for the radiation pressure).* P is the total
//! pressure; P(0) = aT⁴(0)/3 = (a/4)T_eff⁴ q(0) and P(τ₀) = P(0) + gτ₀/κ₀, with
//! κ₀ = κ at (P(τ₀), T(τ₀)) found by iteration. The radiation term is taken at
//! τ = 0, not τ₀. This is adequate where radiation pressure and radiative
//! acceleration are small against gravity (cool stars); it is not a general
//! boundary condition for a non-grey atmosphere.
//!
//! *Surface conditions.* T_eff is defined by L = 4πR²σT_eff⁴, where R is the
//! photospheric radius (T = T_eff). The last mesh point is the fitting point,
//! at r_N = R − Δr with Δr = ∫_{τ_ph}^{τ_fit} dτ/(κρ) from the integration
//! (plane parallel: g = GM/R² and m = M throughout the atmosphere). The two
//! conditions on the last point are T_N = T(τ_fit) and P_N = P(τ_fit).
//!
//! *Interior.* The interior equations keep the uncorrected gradient (A = 1);
//! with the Krishna Swamy relation the fitting point is placed where
//! |q'(τ_fit)| < 10⁻³ (τ_fit = 3.1), the accepted first version of plan 3.2,
//! labelled approximate. Thin-element losses are used on both sides of the
//! fitting point when the atmosphere has convection, so the transport law is
//! the same there.

use serde::{Deserialize, Serialize};

use super::Physics;
use super::convection::{MltInput, mlt};
use super::interp::smoothstep;
use crate::composition::Composition;
use crate::constants::{A_RAD, FOUR_PI, G, SIGMA_SB};

/// Total photospheric pressure for mass m, radius r, temperature t and opacity κ
/// (grey closed form).
pub fn photosphere_pressure(m: f64, r: f64, t: f64, kappa: f64) -> f64 {
    2.0 / 3.0 * G * m / (r * r * kappa) + A_RAD * t.powi(4) / 6.0
}

/// Effective temperature from L = 4π R² σ T_eff⁴.
pub fn effective_temperature(l: f64, r: f64) -> f64 {
    (l.max(0.0) / (FOUR_PI * SIGMA_SB * r * r)).powf(0.25)
}

/// T(τ) relation T⁴ = ¾ T_eff⁴ [τ + q(τ)].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum TTau {
    /// Grey Eddington: q = 2/3.
    #[default]
    Eddington,
    /// Krishna Swamy (1966, eq. 33): q = 1.39 − 0.815 e^{−2.54τ} − 0.025 e^{−30τ},
    /// τ the Rosseland depth. A fit to the observed Sun (convection partly
    /// inside it), supported for 0.02 ≤ τ ≤ 10; its use at τ → 0 and for stars
    /// unlike the Sun is an extrapolation. Above T_eff ≈ 8000 K it is blended
    /// into the grey relation (see [`TTau::hot_weight`]).
    KrishnaSwamy,
}

impl TTau {
    pub fn label(self) -> &'static str {
        match self {
            TTau::Eddington => "grey Eddington",
            TTau::KrishnaSwamy => "Krishna Swamy (1966), solar",
        }
    }

    /// Weight of the grey relation in hot stars: 0 below log T_eff = 3.9
    /// (7943 K), 1 above 4.0, C¹ in between. Hot stars keep the grey Eddington
    /// atmosphere (range of the radiation-pressure approximation, page 6).
    pub fn hot_weight(self, teff: f64) -> f64 {
        match self {
            TTau::Eddington => 1.0,
            TTau::KrishnaSwamy => smoothstep((teff.log10() - 3.9) / 0.1),
        }
    }

    /// q(τ) and q'(τ) for effective temperature `teff`.
    pub fn q(self, tau: f64, teff: f64) -> (f64, f64) {
        let w = self.hot_weight(teff);
        let (ks, dks) = (
            1.39 - 0.815 * (-2.54 * tau).exp() - 0.025 * (-30.0 * tau).exp(),
            0.815 * 2.54 * (-2.54 * tau).exp() + 0.025 * 30.0 * (-30.0 * tau).exp(),
        );
        (w * (2.0 / 3.0) + (1.0 - w) * ks, (1.0 - w) * dks)
    }

    /// Fitting depth: τ = ⅔ for the grey relation; for Krishna Swamy the depth
    /// where |q'| < 10⁻³.
    pub fn default_tau_fit(self) -> f64 {
        match self {
            TTau::Eddington => 2.0 / 3.0,
            TTau::KrishnaSwamy => 3.1,
        }
    }
}

/// Options of the integrated atmosphere.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct IntegratedAtmosphere {
    pub t_tau: TTau,
    /// Optical depth of the fitting point (the last mesh point).
    pub tau_fit: f64,
    /// Starting optical depth τ₀.
    pub tau0: f64,
    /// Number of RK4 steps in ln τ.
    pub steps: usize,
    /// Integrate (P, T) with mixing-length convection (plan 3.2); otherwise T
    /// follows the T(τ) relation (radiative atmosphere, plan 3.1).
    pub convection: bool,
}

impl IntegratedAtmosphere {
    pub fn new(t_tau: TTau, convection: bool) -> Self {
        Self { t_tau, tau_fit: t_tau.default_tau_fit(), tau0: 1e-6, steps: 100, convection }
    }
}

/// Outer boundary condition.
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum Atmosphere {
    /// Grey Eddington photosphere with one κ: P = ⅔ g/κ + aT⁴/6 at τ = ⅔.
    #[default]
    GreyClosedForm,
    /// Integrated atmosphere (plans 3.1–3.3).
    Integrated(IntegratedAtmosphere),
}

impl Atmosphere {
    pub fn describe(&self) -> String {
        match self {
            Atmosphere::GreyClosedForm => "grey Eddington (T = T_eff at τ = 2/3, one κ)".into(),
            Atmosphere::Integrated(a) => format!(
                "integrated, {} T(τ), fit at τ = {}{}",
                a.t_tau.label(),
                a.tau_fit,
                if a.convection { ", with convection" } else { ", radiative" }
            ),
        }
    }

    /// Whether convective elements use the thin-element losses (only with an
    /// integrated atmosphere with convection, on both sides of the fitting point).
    pub fn thin_elements(&self) -> bool {
        matches!(self, Atmosphere::Integrated(a) if a.convection)
    }
}

/// Result of an atmosphere integration.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct AtmosphereSolution {
    /// Pressure and temperature at the fitting point.
    pub p_fit: f64,
    pub t_fit: f64,
    pub teff: f64,
    /// Surface gravity used [cm s⁻²].
    pub g: f64,
    /// Photospheric radius R [cm] and its height above the fitting point.
    pub r_phot: f64,
    pub dr: f64,
    /// Optical depth of the photosphere (T = T_eff).
    pub tau_ph: f64,
    /// Whether any layer of the atmosphere is convective.
    pub convective: bool,
    /// Samples (τ, P, T, ∇, convective) at the RK4 nodes, for tests and output.
    pub tau: Vec<f64>,
    pub p: Vec<f64>,
    pub t: Vec<f64>,
    pub nabla: Vec<f64>,
    pub conv: Vec<bool>,
}

struct Ctx<'a> {
    phys: &'a Physics,
    a: &'a IntegratedAtmosphere,
    comp: &'a Composition,
    alpha: f64,
    g: f64,
    teff: f64,
}

impl Ctx<'_> {
    /// T from the T(τ) relation.
    fn t_relation(&self, tau: f64) -> f64 {
        let (q, _) = self.a.t_tau.q(tau, self.teff);
        (0.75 * (tau + q)).powf(0.25) * self.teff
    }

    /// Local state: (d ln P/ds, ∇, convective, κρ) at s = ln τ, given (ln P, ln T).
    fn rhs(&self, s: f64, lnp: f64, lnt: f64) -> (f64, f64, bool, f64) {
        let tau = s.exp();
        let (p, t) = (lnp.exp(), lnt.exp());
        let eos = self.phys.eos.eval(p, t, self.comp);
        let kappa = self.phys.opacity.kappa(eos.rho, t, self.comp);
        let dlnp = tau * self.g / (kappa * p);
        let (nabla, conv) = if self.a.convection {
            let (_, dq) = self.a.t_tau.q(tau, self.teff);
            let nabla_rad = 3.0 * kappa * p * self.teff.powi(4) / (16.0 * self.g * t.powi(4));
            let r = mlt(&MltInput {
                nabla_rad,
                nabla_ad: eos.nabla_ad,
                t,
                rho: eos.rho,
                p,
                kappa,
                cp: eos.cp,
                delta: eos.delta,
                g: self.g,
                alpha: self.alpha,
                hopf_a: 1.0 + dq,
                thin_elements: true,
            });
            (r.nabla, r.convective)
        } else {
            (0.0, false)
        };
        (dlnp, nabla, conv, kappa * eos.rho)
    }
}

/// Integrate the atmosphere for gravity `g` and effective temperature `teff`.
pub fn integrate(
    phys: &Physics,
    a: &IntegratedAtmosphere,
    alpha: f64,
    comp: &Composition,
    g: f64,
    teff: f64,
) -> AtmosphereSolution {
    let ctx = Ctx { phys, a, comp, alpha, g, teff };
    let tau0 = a.tau0.min(0.5 * a.tau_fit);
    let t0 = ctx.t_relation(tau0);
    let (q0, _) = a.t_tau.q(0.0, teff);
    let p_top = 0.25 * A_RAD * teff.powi(4) * q0;
    // P(τ₀) = P(0) + gτ₀/κ₀, κ₀ at (P(τ₀), T(τ₀)): Newton in ln P with a
    // numerical derivative, to 10⁻¹³.
    let resid = |lnp: f64| -> f64 {
        let p = lnp.exp();
        let rho = phys.eos.eval(p, t0, comp).rho;
        let kappa = phys.opacity.kappa(rho, t0, comp);
        lnp - (p_top + g * tau0 / kappa).ln()
    };
    let mut lnp = (p_top + g * tau0 / phys.opacity.kappa(1e-10, t0, comp)).ln().max((2.0 * p_top).ln());
    for _ in 0..100 {
        let f = resid(lnp);
        let h = 1e-6;
        let df = (resid(lnp + h) - resid(lnp - h)) / (2.0 * h);
        let step = if df > 0.1 { -f / df } else { -0.5 * f };
        let step = step.clamp(-1.0, 1.0);
        lnp += step;
        if step.abs() < 1e-13 {
            break;
        }
    }

    let n = a.steps.max(2);
    let s0 = tau0.ln();
    let s1 = a.tau_fit.ln();
    let h = (s1 - s0) / n as f64;
    // State y = (ln P, ln T, z), z = ∫ dτ/(κρ) from τ₀.
    let deriv = |s: f64, y: [f64; 3]| -> ([f64; 3], f64, bool) {
        // Radiative atmosphere: T from the relation, not from the state.
        let lnt = if a.convection { y[1] } else { ctx.t_relation(s.exp()).ln() };
        let (dlnp, nabla, conv, kr) = ctx.rhs(s, y[0], lnt);
        let dlnt = if a.convection { nabla * dlnp } else { 0.0 };
        ([dlnp, dlnt, s.exp() / kr], nabla, conv)
    };
    let mut y = [lnp, t0.ln(), 0.0];
    let mut sol = AtmosphereSolution { teff, g, ..Default::default() };
    // Samples (s, ln T, d ln T/ds, z, dz/ds) for locating the photosphere.
    let mut samples: Vec<[f64; 5]> = Vec::with_capacity(n + 1);
    let record = |s: f64, y: &[f64; 3], sol: &mut AtmosphereSolution, samples: &mut Vec<[f64; 5]>| {
        let tau = s.exp();
        let (d, nabla, conv) = deriv(s, *y);
        let (lnt, dlnt) = if a.convection {
            (y[1], d[1])
        } else {
            // Radiative: T from the relation; d ln T/ds by differentiating it.
            let (q, dq) = a.t_tau.q(tau, teff);
            (ctx.t_relation(tau).ln(), 0.25 * tau * (1.0 + dq) / (tau + q))
        };
        sol.tau.push(tau);
        sol.p.push(y[0].exp());
        sol.t.push(lnt.exp());
        sol.nabla.push(if a.convection { nabla } else { dlnt / d[0] });
        sol.conv.push(conv);
        samples.push([s, lnt, dlnt, y[2], d[2]]);
    };
    record(s0, &y, &mut sol, &mut samples);
    let mut s = s0;
    for _ in 0..n {
        let add = |y: &[f64; 3], k: &[f64; 3], f: f64| [y[0] + f * k[0], y[1] + f * k[1], y[2] + f * k[2]];
        let (k1, _, _) = deriv(s, y);
        let (k2, _, _) = deriv(s + 0.5 * h, add(&y, &k1, 0.5 * h));
        let (k3, _, _) = deriv(s + 0.5 * h, add(&y, &k2, 0.5 * h));
        let (k4, _, _) = deriv(s + h, add(&y, &k3, h));
        for i in 0..3 {
            y[i] += h / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
        }
        s += h;
        record(s, &y, &mut sol, &mut samples);
    }
    sol.p_fit = y[0].exp();
    sol.t_fit = *sol.t.last().unwrap();
    sol.convective = sol.conv.iter().any(|&c| c);

    // Photosphere: where ln T = ln T_eff, by cubic Hermite interpolation
    // between the RK4 nodes; beyond the fitting point, by extrapolation with
    // the last slopes.
    let target = teff.ln();
    let last = samples[samples.len() - 1];
    let (s_ph, z_ph) = match samples.windows(2).find(|w| (w[0][1] - target) * (w[1][1] - target) <= 0.0) {
        Some(w) => {
            let (a0, b0) = (w[0], w[1]);
            let hh = b0[0] - a0[0];
            let cubic = |t: f64, f0: f64, d0: f64, f1: f64, d1: f64| {
                let (t2, t3) = (t * t, t * t * t);
                (2.0 * t3 - 3.0 * t2 + 1.0) * f0
                    + (t3 - 2.0 * t2 + t) * hh * d0
                    + (-2.0 * t3 + 3.0 * t2) * f1
                    + (t3 - t2) * hh * d1
            };
            // Bisection on the cubic (monotonic in practice), fixed 60 steps.
            let (mut lo, mut hi) = (0.0, 1.0);
            let sign0 = a0[1] - target;
            for _ in 0..60 {
                let mid = 0.5 * (lo + hi);
                if (cubic(mid, a0[1], a0[2], b0[1], b0[2]) - target) * sign0 > 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let t = 0.5 * (lo + hi);
            (a0[0] + t * hh, cubic(t, a0[3], a0[4], b0[3], b0[4]))
        }
        None => {
            let s_ph = last[0] + (target - last[1]) / last[2];
            (s_ph, last[3] + (s_ph - last[0]) * last[4])
        }
    };
    sol.tau_ph = s_ph.exp();
    sol.dr = last[3] - z_ph;
    sol
}

/// Surface conditions for a star of mass `m` whose last mesh point is at
/// radius `r_n` with luminosity `l`: the photospheric radius R = r_n + Δr is
/// found by a fixed number of substitutions (Δr/R ≲ 10⁻³ for dwarfs).
pub fn surface(
    phys: &Physics,
    a: &IntegratedAtmosphere,
    alpha: f64,
    comp: &Composition,
    m: f64,
    r_n: f64,
    l: f64,
) -> AtmosphereSolution {
    let mut r = r_n;
    let mut sol = AtmosphereSolution::default();
    for _ in 0..6 {
        let teff = effective_temperature(l, r);
        sol = integrate(phys, a, alpha, comp, G * m / (r * r), teff);
        let r_new = r_n + sol.dr;
        let done = (r_new - r).abs() <= 1e-13 * r;
        r = r_new;
        if done {
            break;
        }
    }
    sol.r_phot = r;
    sol
}

/// Solve the grey closed form P_ph = ⅔ g / κ(ρ(P_ph, T), T) + a T⁴/6 by damped
/// fixed-point iteration in ln P. The damping makes the iteration converge for
/// opacity laws κ ∝ ρ^a with |a| ≤ 1.
pub fn solve_photosphere(phys: &Physics, m: f64, r: f64, teff: f64, comp: &Composition) -> f64 {
    let mut lnp = 1e5_f64.ln();
    for _ in 0..500 {
        let p = lnp.exp();
        let rho = phys.eos.eval(p, teff, comp).rho;
        let kappa = phys.opacity.kappa(rho, teff, comp);
        let target = photosphere_pressure(m, r, teff, kappa).ln();
        let new = 0.5 * (lnp + target);
        let done = (new - lnp).abs() < 1e-12;
        lnp = new;
        if done {
            break;
        }
    }
    lnp.exp()
}

/// Pressure and temperature of the outermost mesh point of a starting model
/// for mass `m`, radius `r` and luminosity `l`, for any atmosphere.
pub fn starting_surface(phys: &Physics, alpha: f64, m: f64, r: f64, l: f64, comp: &Composition) -> (f64, f64) {
    match &phys.atmosphere {
        Atmosphere::GreyClosedForm => {
            let teff = effective_temperature(l, r);
            (solve_photosphere(phys, m, r, teff, comp), teff)
        }
        Atmosphere::Integrated(a) => {
            let s = surface(phys, a, alpha, comp, m, r, l);
            (s.p_fit, s.t_fit)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{L_SUN, M_SUN, R_SUN, TEFF_SUN};
    use crate::physics::{Opacity, eos::IdealGasRadiation, nuclear::AnalyticNuclear};

    fn sun() -> Composition {
        Composition::new(0.7, 0.02, 0.7).unwrap()
    }

    #[test]
    fn solar_teff() {
        let t = effective_temperature(L_SUN, R_SUN);
        assert!((t - TEFF_SUN).abs() < 5.0, "T_eff = {t}");
    }

    #[test]
    fn solar_photosphere_pressure_is_reasonable() {
        let phys = Physics::analytic();
        let comp = sun();
        let p = solve_photosphere(&phys, M_SUN, R_SUN, TEFF_SUN, &comp);
        // The real solar photosphere has P ~ 1e5 dyn/cm²; the analytic opacity is crude.
        assert!(p > 1e3 && p < 1e7, "P_ph = {p}");
        let rho = phys.eos.eval(p, TEFF_SUN, &comp).rho;
        let kappa = phys.opacity.kappa(rho, TEFF_SUN, &comp);
        let resid = p / photosphere_pressure(M_SUN, R_SUN, TEFF_SUN, kappa) - 1.0;
        assert!(resid.abs() < 1e-8);
    }

    struct ConstKappa(f64);
    impl Opacity for ConstKappa {
        fn name(&self) -> &'static str {
            "constant"
        }
        fn kappa(&self, _rho: f64, _t: f64, _c: &Composition) -> f64 {
            self.0
        }
    }

    fn const_physics(k: f64) -> Physics {
        Physics {
            eos: Box::new(IdealGasRadiation),
            opacity: Box::new(ConstKappa(k)),
            nuclear: Box::new(AnalyticNuclear),
            atmosphere: Atmosphere::GreyClosedForm,
        }
    }

    /// Plan 3.1: constant κ with q = 2/3 integrated to τ = 2/3 equals the
    /// closed form to 10⁻⁶, independently of τ₀.
    #[test]
    fn constant_kappa_reproduces_the_closed_form() {
        let g = G * M_SUN / (R_SUN * R_SUN);
        for k in [1e-3, 0.1, 10.0] {
            let phys = const_physics(k);
            let want = photosphere_pressure(M_SUN, R_SUN, TEFF_SUN, k);
            for tau0 in [1e-3, 1e-4, 1e-6] {
                let a = IntegratedAtmosphere { tau0, ..IntegratedAtmosphere::new(TTau::Eddington, false) };
                let s = integrate(&phys, &a, 1.8, &sun(), g, TEFF_SUN);
                assert!((s.p_fit / want - 1.0).abs() < 1e-6, "κ = {k}, τ₀ = {tau0}: {} vs {want}", s.p_fit);
                assert!((s.t_fit / TEFF_SUN - 1.0).abs() < 1e-14);
                assert!(s.dr.abs() < 1e-6 * R_SUN);
            }
        }
    }

    /// Result independent of τ₀ with the real (analytic) opacity, and of the
    /// number of steps.
    #[test]
    fn independent_of_tau0_and_steps() {
        let phys = Physics::analytic();
        let g = G * M_SUN / (R_SUN * R_SUN);
        for t_tau in [TTau::Eddington, TTau::KrishnaSwamy] {
            let base = IntegratedAtmosphere::new(t_tau, false);
            let r = integrate(&phys, &base, 1.8, &sun(), g, TEFF_SUN).p_fit;
            for (tau0, steps) in [(1e-7, 100), (1e-6, 200), (1e-8, 400)] {
                let a = IntegratedAtmosphere { tau0, steps, ..base };
                let p = integrate(&phys, &a, 1.8, &sun(), g, TEFF_SUN).p_fit;
                assert!((p / r - 1.0).abs() < 2e-6, "{t_tau:?} τ₀ = {tau0}, {steps} steps: {p} vs {r}");
            }
        }
    }

    /// Plan 3.2: a radiative atmosphere integrated as a (P, T) pair reproduces
    /// the T(τ) relation.
    #[test]
    fn radiative_pair_reproduces_t_tau() {
        // Hot enough that the analytic opacity keeps the atmosphere radiative.
        let phys = Physics::analytic();
        let g = 1e4;
        let teff = 9000.0;
        for t_tau in [TTau::Eddington, TTau::KrishnaSwamy] {
            let a = IntegratedAtmosphere::new(t_tau, true);
            let s = integrate(&phys, &a, 1.8, &sun(), g, teff);
            // (The Krishna Swamy relation has A = 1 + q' up to 3.8 at the top,
            // where A∇_rad can exceed ∇_ad; convection there is inefficient and
            // ∇ ≈ A∇_rad, so the comparison still holds.)
            assert!(t_tau == TTau::KrishnaSwamy || !s.convective);
            let b = IntegratedAtmosphere::new(t_tau, false);
            let r = integrate(&phys, &b, 1.8, &sun(), g, teff);
            for k in 0..s.tau.len() {
                assert!((s.t[k] / r.t[k] - 1.0).abs() < 2e-6, "{t_tau:?} τ = {}: {} vs {}", s.tau[k], s.t[k], r.t[k]);
            }
            assert!((s.p_fit / r.p_fit - 1.0).abs() < 2e-6);
        }
    }

    #[test]
    fn krishna_swamy_relation() {
        // At τ = 0: T/T_eff = (¾ q(0))^¼ = 0.801 (paper's fit; tabulated 0.800).
        let (q0, _) = TTau::KrishnaSwamy.q(0.0, 5000.0);
        assert!(((0.75 * q0).powf(0.25) - 0.801).abs() < 5e-4);
        // |q'| < 1e-3 at the default fitting point.
        let (_, dq) = TTau::KrishnaSwamy.q(3.1, 5000.0);
        assert!(dq < 1e-3);
        // Hot stars: grey.
        assert_eq!(TTau::KrishnaSwamy.q(0.3, 12000.0), (2.0 / 3.0, 0.0));
    }

    /// The photosphere of a radiative Krishna Swamy atmosphere lies where
    /// τ + q(τ) = 4/3, above the fitting point.
    #[test]
    fn photosphere_location() {
        let phys = Physics::analytic();
        let g = G * M_SUN / (R_SUN * R_SUN);
        let a = IntegratedAtmosphere::new(TTau::KrishnaSwamy, false);
        let s = integrate(&phys, &a, 1.8, &sun(), g, TEFF_SUN);
        let (q, _) = TTau::KrishnaSwamy.q(s.tau_ph, TEFF_SUN);
        assert!((s.tau_ph + q - 4.0 / 3.0).abs() < 1e-6, "τ_ph = {}", s.tau_ph);
        assert!(s.dr > 0.0 && s.dr < 1e-3 * R_SUN, "Δr = {}", s.dr);
    }
}
