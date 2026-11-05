//! Input parameters and the converged stellar model (profiles + global properties).

use serde::{Deserialize, Serialize};

use super::equations::{I_L, I_LNP, I_LNR, I_LNT, NVAR, StructureProblem};
use crate::composition::{Composition, CompositionError};
use crate::constants::{G, L_SUN, M_SUN, R_SUN};

/// User-facing model parameters.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct StarParams {
    /// Total mass [M☉].
    pub mass: f64,
    /// Hydrogen mass fraction X.
    pub x: f64,
    /// Metallicity Z (Y = 1 − X − Z).
    pub z: f64,
    /// Mixing-length parameter α_MLT.
    pub alpha_mlt: f64,
    /// Fraction of Z in C, N and O (feeds the CNO cycle).
    pub f_cno: f64,
    /// Number of mesh points.
    pub n_points: usize,
}

impl Default for StarParams {
    fn default() -> Self {
        Self::sun()
    }
}

impl StarParams {
    pub const MASS_RANGE: (f64, f64) = (0.1, 100.0);
    pub const POINTS_RANGE: (usize, usize) = (100, 5000);

    /// The Sun with the classic composition X = 0.70, Y = 0.28, Z = 0.02.
    pub fn sun() -> Self {
        Self { mass: 1.0, x: 0.70, z: 0.02, alpha_mlt: 1.8, f_cno: 0.7, n_points: 800 }
    }

    /// The Sun with AGSS09 protosolar composition (Asplund et al. 2009).
    pub fn sun_agss09() -> Self {
        Self { x: 0.7154, z: 0.0142, ..Self::sun() }
    }

    pub fn y(&self) -> f64 {
        1.0 - self.x - self.z
    }

    pub fn composition(&self) -> Result<Composition, CompositionError> {
        Composition::new(self.x, self.z, self.f_cno)
    }

    pub fn validate(&self) -> Result<(), String> {
        let (mlo, mhi) = Self::MASS_RANGE;
        if !(mlo..=mhi).contains(&self.mass) {
            return Err(format!("mass must be in [{mlo}, {mhi}] M☉"));
        }
        self.composition().map_err(|e| e.to_string())?;
        if !(self.alpha_mlt > 0.0 && self.alpha_mlt <= 5.0) {
            return Err("α_MLT must be in (0, 5]".into());
        }
        let (plo, phi) = Self::POINTS_RANGE;
        if !(plo..=phi).contains(&self.n_points) {
            return Err(format!("mesh points must be in [{plo}, {phi}]"));
        }
        Ok(())
    }

    /// Path through parameter space for continuation: geometric in mass,
    /// linear in everything else.
    pub fn interpolate(&self, other: &Self, t: f64) -> Self {
        let l = |a: f64, b: f64| a + (b - a) * t;
        Self {
            mass: (self.mass.ln() + (other.mass.ln() - self.mass.ln()) * t).exp(),
            x: l(self.x, other.x),
            z: l(self.z, other.z),
            alpha_mlt: l(self.alpha_mlt, other.alpha_mlt),
            f_cno: l(self.f_cno, other.f_cno),
            n_points: other.n_points,
        }
    }

    /// Mass below which the ideal-gas EOS is a poor approximation.
    pub fn ideal_gas_warning(&self) -> Option<&'static str> {
        (self.mass < 0.4).then_some(
            "below ~0.4 M☉ electron degeneracy and non-ideal effects matter; \
             the ideal-gas EOS is inaccurate",
        )
    }
}

/// Radial profiles at every mesh point, in CGS units.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Profile {
    /// m / M.
    pub q: Vec<f64>,
    pub m: Vec<f64>,
    pub r: Vec<f64>,
    pub p: Vec<f64>,
    pub t: Vec<f64>,
    pub rho: Vec<f64>,
    pub l: Vec<f64>,
    pub kappa: Vec<f64>,
    pub eps_pp: Vec<f64>,
    pub eps_cno: Vec<f64>,
    pub eps_3a: Vec<f64>,
    pub eps: Vec<f64>,
    pub nabla_rad: Vec<f64>,
    pub nabla_ad: Vec<f64>,
    pub nabla: Vec<f64>,
    pub convective: Vec<bool>,
    pub v_conv: Vec<f64>,
    pub f_conv: Vec<f64>,
    pub hp: Vec<f64>,
    pub mu: Vec<f64>,
    pub beta: Vec<f64>,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl Profile {
    pub fn len(&self) -> usize {
        self.q.len()
    }

    pub fn is_empty(&self) -> bool {
        self.q.is_empty()
    }
}

/// A convective region. Radii are in units of the stellar radius.
#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct ConvectiveZone {
    pub q_inner: f64,
    pub q_outer: f64,
    pub r_inner: f64,
    pub r_outer: f64,
}

/// Global properties of a converged model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Globals {
    /// [M☉]
    pub mass: f64,
    /// [R☉]
    pub radius: f64,
    /// [L☉]
    pub luminosity: f64,
    /// [K]
    pub teff: f64,
    /// log₁₀ g [cgs]
    pub log_g: f64,
    /// Central temperature [K].
    pub tc: f64,
    /// Central density [g cm^-3].
    pub rhoc: f64,
    /// Central pressure [dyn cm^-2].
    pub pc: f64,
    /// Mass fraction of the convective core (0 if radiative core).
    pub conv_core_q: f64,
    /// Radius of the convective core [r/R].
    pub conv_core_r: f64,
    /// Base of the convective envelope [r/R], if any.
    pub conv_env_base_r: Option<f64>,
    /// Mass coordinate of the base of the convective envelope [m/M], if any.
    pub conv_env_base_q: Option<f64>,
    /// Fractions of the luminosity from pp, CNO, 3α.
    pub frac_pp: f64,
    pub frac_cno: f64,
    pub frac_3a: f64,
    /// Virial check: ∫3P/ρ dm / ∫G m/r dm (should be ≈ 1).
    pub virial_ratio: f64,
    /// Energy check: ∫ε dm / L_surface (should be ≈ 1).
    pub energy_ratio: f64,
    pub newton_iterations: usize,
    pub max_correction: f64,
}

/// A converged stellar model.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct StarModel {
    pub params: StarParams,
    pub physics: String,
    pub globals: Globals,
    pub zones: Vec<ConvectiveZone>,
    pub profile: Profile,
    /// Raw solution vector (ln r, ln P, ln T, L/L_scale per point), kept for warm starts.
    #[serde(skip)]
    pub solution: Vec<f64>,
    #[serde(skip)]
    pub l_scale: f64,
}

/// Evaluate the profiles and global quantities from a converged solution.
pub(crate) fn build_model(
    params: &StarParams,
    problem: &StructureProblem<'_>,
    q: &[f64],
    y: Vec<f64>,
    newton_iterations: usize,
    max_correction: f64,
) -> StarModel {
    let n = q.len();
    let mut pr = Profile::default();
    for k in 0..n {
        let s = problem.local_at(k, &y[k * NVAR..(k + 1) * NVAR]);
        let c = &problem.comp[k];
        pr.q.push(q[k]);
        pr.m.push(problem.m[k]);
        pr.r.push(s.r);
        pr.p.push(s.p);
        pr.t.push(s.t);
        pr.rho.push(s.rho());
        pr.l.push(s.l);
        pr.kappa.push(s.kappa);
        pr.eps_pp.push(s.nuc.pp);
        pr.eps_cno.push(s.nuc.cno);
        pr.eps_3a.push(s.nuc.triple_alpha);
        pr.eps.push(s.eps());
        pr.nabla_rad.push(s.nabla_rad);
        pr.nabla_ad.push(s.eos.nabla_ad);
        pr.nabla.push(s.nabla());
        pr.convective.push(s.mlt.convective);
        pr.v_conv.push(s.mlt.v_conv);
        pr.f_conv.push(s.mlt.f_conv);
        pr.hp.push(s.hp);
        pr.mu.push(s.eos.mu);
        pr.beta.push(s.eos.beta);
        pr.x.push(c.x);
        pr.y.push(c.y);
        pr.z.push(c.z);
    }
    // The outermost point sits on the photosphere; there the local ∇ is
    // dominated by the boundary condition, not by transport, so classify it
    // like its neighbour.
    if n >= 2 {
        pr.convective[n - 1] = pr.convective[n - 2];
    }

    let radius = pr.r[n - 1];
    let lum = pr.l[n - 1];
    let mass = params.mass * M_SUN;

    let zones = find_zones(&pr, radius);
    let (conv_core_q, conv_core_r) =
        zones.first().filter(|z| z.q_inner <= q[0]).map(|z| (z.q_outer, z.r_outer)).unwrap_or((0.0, 0.0));
    let env = zones.last().filter(|z| z.q_outer >= 1.0 - 1e-6 && z.q_inner > q[0]);

    // Integrals with the trapezoid rule over mass.
    let trap = |f: &dyn Fn(usize) -> f64| -> f64 {
        (1..n).map(|k| 0.5 * (f(k) + f(k - 1)) * (pr.m[k] - pr.m[k - 1])).sum::<f64>()
    };
    let e_pp = trap(&|k| pr.eps_pp[k]);
    let e_cno = trap(&|k| pr.eps_cno[k]);
    let e_3a = trap(&|k| pr.eps_3a[k]);
    let e_tot = (e_pp + e_cno + e_3a).max(f64::MIN_POSITIVE);
    let thermal = trap(&|k| 3.0 * pr.p[k] / pr.rho[k]);
    let grav = trap(&|k| G * pr.m[k] / pr.r[k]);

    let globals = Globals {
        mass: params.mass,
        radius: radius / R_SUN,
        luminosity: lum / L_SUN,
        teff: pr.t[n - 1],
        log_g: (G * mass / (radius * radius)).log10(),
        tc: pr.t[0],
        rhoc: pr.rho[0],
        pc: pr.p[0],
        conv_core_q,
        conv_core_r,
        conv_env_base_r: env.map(|z| z.r_inner),
        conv_env_base_q: env.map(|z| z.q_inner),
        frac_pp: e_pp / e_tot,
        frac_cno: e_cno / e_tot,
        frac_3a: e_3a / e_tot,
        virial_ratio: thermal / grav,
        energy_ratio: (e_tot + pr.l[0]) / lum,
        newton_iterations,
        max_correction,
    };

    StarModel {
        params: params.clone(),
        physics: problem.physics.describe(),
        globals,
        zones,
        profile: pr,
        solution: y,
        l_scale: problem.l_scale,
    }
}

/// Locate convective zones and interpolate their boundaries linearly in
/// ∇_rad − ∇_ad.
fn find_zones(pr: &Profile, radius: f64) -> Vec<ConvectiveZone> {
    let n = pr.len();
    let w = |k: usize| pr.nabla_rad[k] - pr.nabla_ad[k];
    let boundary = |k: usize| -> (f64, f64) {
        // Transition lies between k−1 and k.
        let (w0, w1) = (w(k - 1), w(k));
        let f = if (w1 - w0).abs() > 0.0 { (-w0 / (w1 - w0)).clamp(0.0, 1.0) } else { 0.5 };
        let qb = pr.q[k - 1] + f * (pr.q[k] - pr.q[k - 1]);
        let rb = pr.r[k - 1] + f * (pr.r[k] - pr.r[k - 1]);
        (qb, rb / radius)
    };
    let mut zones = Vec::new();
    let mut start: Option<(f64, f64)> = None;
    for k in 0..n {
        let c = pr.convective[k];
        match (c, start) {
            (true, None) => {
                start = Some(if k == 0 { (pr.q[0], 0.0) } else { boundary(k) });
            }
            (false, Some((qi, ri))) => {
                let (qo, ro) = boundary(k);
                zones.push(ConvectiveZone { q_inner: qi, q_outer: qo, r_inner: ri, r_outer: ro });
                start = None;
            }
            _ => {}
        }
    }
    if let Some((qi, ri)) = start {
        zones.push(ConvectiveZone { q_inner: qi, q_outer: 1.0, r_inner: ri, r_outer: 1.0 });
    }
    zones
}

/// Basic physical sanity checks on a converged model.
pub(crate) fn sanity_check(model: &StarModel) -> Result<(), String> {
    let p = &model.profile;
    let finite = p.r.iter().chain(&p.p).chain(&p.t).chain(&p.rho).chain(&p.l).all(|v| v.is_finite());
    if !finite {
        return Err("non-finite values in converged model".into());
    }
    if !p.r.windows(2).all(|w| w[1] > w[0]) {
        return Err("radius is not monotonic".into());
    }
    if !p.p.windows(2).all(|w| w[1] < w[0]) {
        return Err("pressure does not decrease outward".into());
    }
    if !(model.globals.luminosity > 0.0) {
        return Err("non-positive surface luminosity".into());
    }
    if p.beta.iter().any(|&b| b <= crate::physics::eos::MIN_BETA * 1.0001) {
        return Err("radiation pressure exceeds total pressure somewhere".into());
    }
    Ok(())
}

/// Extract (ln r, ln P, ln T, l) from the raw solution of mesh point `k`.
pub fn unknowns_at(model: &StarModel, k: usize) -> [f64; 4] {
    let y = &model.solution[k * NVAR..(k + 1) * NVAR];
    [y[I_LNR], y[I_LNP], y[I_LNT], y[I_L]]
}
