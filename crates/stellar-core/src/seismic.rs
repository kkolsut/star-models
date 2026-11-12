//! Seismic diagnostics of a model: sound speed, Lamb and buoyancy (Brunt–Väisälä)
//! frequencies, and the asymptotic spacings Δν and ΔΠ₁.
//!
//! Everything here is evaluated from the profiles of a converged model; nothing
//! feeds back into the structure. See docs/physics/09-pulsations.md,
//! "Characteristic frequencies" and "Asymptotic spacings".

use serde::{Deserialize, Serialize};

use crate::structure::model::Profile;

/// φ = (∂ln ρ/∂ln μ)_{P,T}. Exactly 1 for the implemented EOS (ideal gas plus
/// radiation, composition entering only through μ). Replaced by the general
/// composition term B of docs/physics/05 when a non-ideal EOS is added (plan 6.2).
pub const PHI: f64 = 1.0;

/// Smallest dimensionless buoyancy A = N² P / (g² ρ) counted as radiative
/// (stably stratified) when locating g-mode cavities. Efficient convection has
/// |A| ≲ 10⁻⁶ (∇ − ∇_ad ≤ 1.4 × 10⁻⁷ in the 5 M☉ core); radiative zones have
/// A ~ 0.1. The tolerance keeps numerical noise of either sign in nearly
/// adiabatic layers from opening spurious cavities.
pub const A_RADIATIVE: f64 = 1e-4;

/// A connected region with N² > 0: a *candidate* g-mode cavity.
///
/// The value of ΔΠ₁ uses the whole region and is the low-frequency asymptotic
/// estimate; it is not a predicted spectrum (see docs/physics/09-pulsations.md).
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct GCavity {
    /// Inner and outer edge [r/R].
    pub r_inner: f64,
    pub r_outer: f64,
    /// Inner and outer edge [m/M].
    pub q_inner: f64,
    pub q_outer: f64,
    /// ∫ N dr / r over the region [s⁻¹].
    pub n_integral: f64,
    /// Asymptotic dipole period spacing ΔΠ₁ = 2π² / √2 / ∫ N dr/r [s].
    pub delta_pi1_s: f64,
    /// Where the region lies.
    pub location: CavityLocation,
    /// A thin layer (less than 0.1 in ln r, not at the centre): listed, but not a
    /// g-mode cavity in the asymptotic sense.
    pub thin: bool,
}

/// Position of a candidate g-mode cavity in the star.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum CavityLocation {
    /// Starts at the innermost mesh point (radiative core).
    Core,
    /// Between two convective regions.
    Intermediate,
    /// Reaches the surface: a radiative envelope, or the thin radiative layer
    /// above a convective envelope.
    Surface,
}

impl GCavity {
    /// Short description of where the region lies.
    pub fn label(&self) -> &'static str {
        match (self.location, self.thin) {
            (CavityLocation::Core, _) => "radiative core",
            (CavityLocation::Intermediate, false) => "radiative zone",
            (CavityLocation::Intermediate, true) => "thin radiative layer",
            (CavityLocation::Surface, false) => "radiative envelope",
            (CavityLocation::Surface, true) => "thin surface layer",
        }
    }
}

/// Sound speed c_s = √(Γ₁ P / ρ) [cm s⁻¹].
pub fn sound_speed(gamma1: f64, p: f64, rho: f64) -> f64 {
    (gamma1 * p / rho).sqrt()
}

/// Lamb frequency S_l = √(l(l+1)) c_s / r [rad s⁻¹].
pub fn lamb_frequency(l: u32, cs: f64, r: f64) -> f64 {
    let l = f64::from(l);
    (l * (l + 1.0)).sqrt() * cs / r
}

/// ∇_μ = d ln μ / d ln P along the profile, by centred differences (one-sided
/// at the ends).
pub fn nabla_mu(mu: &[f64], p: &[f64]) -> Vec<f64> {
    let n = mu.len();
    if n < 2 {
        return vec![0.0; n];
    }
    (0..n)
        .map(|k| {
            let (a, b) = (k.saturating_sub(1), (k + 1).min(n - 1));
            let dlnp = (p[b] / p[a]).ln();
            if dlnp.abs() > 0.0 { (mu[b] / mu[a]).ln() / dlnp } else { 0.0 }
        })
        .collect()
}

/// Buoyancy frequency squared from the production formula
/// N² = (g² ρ / P) [δ (∇_ad − ∇) + φ ∇_μ] [s⁻²], at every point.
///
/// `g` is the magnitude of the local gravity G m / r², `nabla` the actual
/// temperature gradient of the model.
pub fn buoyancy_n2(
    g: &[f64],
    rho: &[f64],
    p: &[f64],
    delta: &[f64],
    nabla_ad: &[f64],
    nabla: &[f64],
    mu: &[f64],
) -> Vec<f64> {
    let nmu = nabla_mu(mu, p);
    (0..g.len()).map(|k| g[k] * g[k] * rho[k] / p[k] * (delta[k] * (nabla_ad[k] - nabla[k]) + PHI * nmu[k])).collect()
}

/// N² from its definition, N² = g (Γ₁⁻¹ d ln P/dr − d ln ρ/dr), by centred
/// differences on the mesh. A cross-check of [`buoyancy_n2`] only: it
/// differences nearly equal numbers and is noisy in the deep interior.
pub fn n2_from_definition(pr: &Profile, gravity: &[f64]) -> Vec<f64> {
    let n = pr.len();
    (0..n)
        .map(|k| {
            let (a, b) = (k.saturating_sub(1), (k + 1).min(n - 1));
            let dr = pr.r[b] - pr.r[a];
            if !(dr > 0.0) {
                return 0.0;
            }
            let dlnp = (pr.p[b] / pr.p[a]).ln() / dr;
            let dlnrho = (pr.rho[b] / pr.rho[a]).ln() / dr;
            gravity[k] * (dlnp / pr.gamma1[k] - dlnrho)
        })
        .collect()
}

/// Dimensionless buoyancy A = N² P / (g² ρ) at point `k` of a profile.
pub fn dimensionless_buoyancy(pr: &Profile, gravity: &[f64], k: usize) -> f64 {
    let g2 = gravity[k] * gravity[k];
    if g2 > 0.0 { pr.n2[k] * pr.p[k] / (g2 * pr.rho[k]) } else { 0.0 }
}

/// Asymptotic large frequency separation Δν = [2 ∫₀^R dr / c_s]⁻¹ [μHz].
/// The integral from the centre to the innermost point uses c_s there.
pub fn delta_nu_asymptotic_uhz(r: &[f64], cs: &[f64]) -> f64 {
    let n = r.len();
    if n < 2 {
        return f64::NAN;
    }
    let mut t = r[0] / cs[0];
    for k in 1..n {
        t += 0.5 * (1.0 / cs[k] + 1.0 / cs[k - 1]) * (r[k] - r[k - 1]);
    }
    1e6 / (2.0 * t)
}

/// Candidate g-mode cavities: connected regions with A > [`A_RADIATIVE`],
/// centre to surface. Disconnected regions are kept apart, never summed.
pub fn g_cavities(pr: &Profile, gravity: &[f64]) -> Vec<GCavity> {
    let n = pr.len();
    if n < 2 {
        return Vec::new();
    }
    let radius = pr.r[n - 1];
    let stable: Vec<bool> = (0..n).map(|k| dimensionless_buoyancy(pr, gravity, k) > A_RADIATIVE).collect();
    let mut out = Vec::new();
    let mut k = 0;
    while k < n {
        if !stable[k] {
            k += 1;
            continue;
        }
        let a = k;
        while k + 1 < n && stable[k + 1] {
            k += 1;
        }
        let b = k;
        k += 1;
        if b == a {
            continue;
        }
        let nf = |i: usize| pr.n2[i].max(0.0).sqrt() / pr.r[i];
        let integral: f64 = (a + 1..=b).map(|i| 0.5 * (nf(i) + nf(i - 1)) * (pr.r[i] - pr.r[i - 1])).sum();
        if !(integral > 0.0) {
            continue;
        }
        let location = if a == 0 {
            CavityLocation::Core
        } else if b == n - 1 {
            CavityLocation::Surface
        } else {
            CavityLocation::Intermediate
        };
        out.push(GCavity {
            r_inner: if a == 0 { 0.0 } else { pr.r[a] / radius },
            r_outer: pr.r[b] / radius,
            q_inner: if a == 0 { 0.0 } else { pr.q[a] },
            q_outer: pr.q[b],
            n_integral: integral,
            delta_pi1_s: 2.0 * std::f64::consts::PI.powi(2) / 2f64.sqrt() / integral,
            location,
            thin: a != 0 && (pr.r[b] / pr.r[a]).ln() < 0.1,
        });
    }
    out
}

/// The cavity with the largest ∫ N dr/r among those that are not thin layers
/// (the one with the smallest ΔΠ₁), if any. A convenience for one-number
/// displays; the full list is in `Globals::g_cavities`.
pub fn main_cavity(cavities: &[GCavity]) -> Option<&GCavity> {
    cavities.iter().filter(|c| !c.thin).max_by(|a, b| a.n_integral.total_cmp(&b.n_integral))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn nabla_mu_of_power_law() {
        // μ ∝ P^0.1 gives ∇_μ = 0.1 everywhere.
        let p: Vec<f64> = (0..20).map(|i| 10f64.powf(17.0 - 0.5 * i as f64)).collect();
        let mu: Vec<f64> = p.iter().map(|p| 0.6 * (p / 1e17).powf(0.1)).collect();
        for v in nabla_mu(&mu, &p) {
            assert!((v - 0.1).abs() < 1e-12, "{v}");
        }
    }

    #[test]
    fn delta_nu_of_uniform_sound_speed() {
        // c_s constant: Δν = c_s / (2R).
        let r: Vec<f64> = (1..=100).map(|i| i as f64 * 1e9).collect();
        let cs = vec![1e7; r.len()];
        let expect = 1e7 / (2.0 * 1e11) * 1e6;
        assert!((delta_nu_asymptotic_uhz(&r, &cs) / expect - 1.0).abs() < 1e-12);
    }

    fn profile(n2: &[f64]) -> (Profile, Vec<f64>) {
        let n = n2.len();
        let r: Vec<f64> = (1..=n).map(|i| i as f64 * 1e9).collect();
        let pr = Profile {
            q: (1..=n).map(|i| i as f64 / n as f64).collect(),
            p: (0..n).map(|i| 1e17 * (1.0 - 0.5 * i as f64 / n as f64)).collect(),
            rho: vec![1.0; n],
            n2: n2.to_vec(),
            r,
            ..Default::default()
        };
        let g = (0..n).map(|i| 1e5 + 1e3 * i as f64).collect();
        (pr, g)
    }

    #[test]
    fn no_cavity_without_a_stable_region() {
        let (pr, g) = profile(&[-1e-8; 50]);
        assert!(g_cavities(&pr, &g).is_empty());
        // Noise of either sign around zero (nearly adiabatic) opens no cavity either.
        let noise: Vec<f64> = (0..50).map(|i| if i % 2 == 0 { 1e-16 } else { -1e-16 }).collect();
        let (pr, g) = profile(&noise);
        assert!(g_cavities(&pr, &g).is_empty());
    }

    #[test]
    fn disconnected_regions_are_kept_apart() {
        // Stable – convective – stable – convective.
        let n2: Vec<f64> = (0..60).map(|i| if i < 15 || (30..45).contains(&i) { 1e-6 } else { -1e-6 }).collect();
        let (pr, g) = profile(&n2);
        let c = g_cavities(&pr, &g);
        assert_eq!(c.len(), 2);
        assert_eq!(c[0].location, CavityLocation::Core);
        assert_eq!(c[1].location, CavityLocation::Intermediate);
        for cav in &c {
            assert!(cav.delta_pi1_s.is_finite() && cav.delta_pi1_s > 0.0);
        }
        // ΔΠ₁ is per cavity: the two values are not combined.
        assert!(c[0].delta_pi1_s != c[1].delta_pi1_s);
    }

    #[test]
    fn lamb_frequency_l1() {
        assert!((lamb_frequency(1, 2.0, 1.0) - 2.0 * 2f64.sqrt()).abs() < 1e-15);
    }
}
