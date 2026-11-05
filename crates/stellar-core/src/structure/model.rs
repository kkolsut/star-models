//! Input parameters and the converged stellar model (profiles + global properties).

use serde::{Deserialize, Serialize};

use super::equations::{I_L, I_LNP, I_LNR, I_LNT, NVAR, StructureProblem};
use crate::composition::{Composition, CompositionError};
use crate::constants::{FOUR_PI, G, L_SUN, M_SUN, R_SUN, TEFF_SUN};
use crate::solver::ResidualSummary;

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
    /// Gravitational energy release ε_grav = −T ds/dt [erg g⁻¹ s⁻¹] (0 on the ZAMS).
    pub eps_grav: Vec<f64>,
    pub nabla_rad: Vec<f64>,
    pub nabla_ad: Vec<f64>,
    pub nabla: Vec<f64>,
    pub convective: Vec<bool>,
    pub v_conv: Vec<f64>,
    pub f_conv: Vec<f64>,
    pub hp: Vec<f64>,
    pub mu: Vec<f64>,
    pub beta: Vec<f64>,
    /// First adiabatic exponent Γ₁.
    pub gamma1: Vec<f64>,
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

/// Evolutionary phase of a model.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum Phase {
    /// Zero-age main sequence: homogeneous, just started burning hydrogen.
    #[default]
    Zams,
    /// Core hydrogen burning.
    MainSequence,
    /// Hydrogen exhausted in the centre, burning in a shell (subgiant / Hertzsprung gap).
    ShellHydrogen,
    /// Shell burning, expanded to more than 5× the ZAMS radius and cooler than 5000 K.
    RedGiant,
    /// Helium burning in the core.
    CoreHelium,
}

impl Phase {
    pub fn label(self) -> &'static str {
        match self {
            Phase::Zams => "zero-age main sequence",
            Phase::MainSequence => "main sequence (core H burning)",
            Phase::ShellHydrogen => "H-shell burning (subgiant / Hertzsprung gap)",
            Phase::RedGiant => "red giant (H-shell burning)",
            Phase::CoreHelium => "core He burning",
        }
    }
}

/// Global properties of a converged model.
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Globals {
    /// Age since the ZAMS [yr].
    #[serde(default)]
    pub age_yr: f64,
    #[serde(default)]
    pub phase: Phase,
    /// Central hydrogen and helium mass fractions.
    #[serde(default)]
    pub x_c: f64,
    #[serde(default)]
    pub y_c: f64,
    /// Mass fraction of the hydrogen-exhausted core (X < 10⁻³), 0 on the main sequence.
    #[serde(default)]
    pub he_core_q: f64,
    /// Mean density 3M / (4πR³) [g cm⁻³].
    #[serde(default)]
    pub mean_density: f64,
    /// Asteroseismic scaling: frequency of maximum oscillation power [μHz].
    #[serde(default)]
    pub nu_max_uhz: f64,
    /// Asteroseismic scaling: large frequency separation [μHz].
    #[serde(default)]
    pub delta_nu_uhz: f64,
    /// Fraction of the luminosity supplied by gravitational contraction (∫ε_grav dm / L).
    #[serde(default)]
    pub frac_grav: f64,
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
    /// Energy check: (Σ Δm (ε_nuc + ε_grav) + L₀) / L_surface with the solver's midpoint rule (≈ 1).
    /// This is a discretisation consistency check: it integrates the same
    /// difference equation the solver satisfies, so it holds by construction
    /// for a converged model. It is not an independent energy-conservation test.
    pub energy_ratio: f64,
    pub newton_iterations: usize,
    /// Largest Newton correction in the last iteration (the convergence test).
    pub max_correction: f64,
    /// Largest residuals of the discretised equations of the converged model.
    #[serde(default)]
    pub residuals: Residuals,
}

/// Largest absolute residual of each discretised structure equation over the
/// mesh. All residuals are dimensionless (differences of ln r, ln P, ln T, and
/// of L in units of the luminosity scale).
#[derive(Clone, Copy, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct Residuals {
    pub continuity: f64,
    pub hydrostatic: f64,
    pub energy: f64,
    pub transport: f64,
    /// Inner and outer boundary conditions.
    pub boundary: f64,
}

impl Residuals {
    pub fn from_summary(s: &ResidualSummary) -> Self {
        let at = |i: usize| s.interval.get(i).copied().unwrap_or(0.0);
        Self {
            continuity: at(I_LNR),
            hydrostatic: at(I_LNP),
            energy: at(I_L),
            transport: at(I_LNT),
            boundary: s.boundary,
        }
    }

    /// Largest of all residuals.
    pub fn max(&self) -> f64 {
        [self.continuity, self.hydrostatic, self.energy, self.transport, self.boundary].into_iter().fold(0.0, f64::max)
    }
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
    /// Fundamental radial pulsation mode, if computed.
    #[serde(default)]
    pub pulsation: Option<crate::pulsation::RadialMode>,
}

/// Evaluate the profiles and global quantities from a converged solution.
pub(crate) fn build_model(
    params: &StarParams,
    problem: &StructureProblem<'_>,
    q: &[f64],
    y: Vec<f64>,
    newton_iterations: usize,
    max_correction: f64,
    residuals: Residuals,
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
        pr.eps_grav.push(problem.eps_grav_at(k, &y[k * NVAR..(k + 1) * NVAR]));
        pr.nabla_rad.push(s.nabla_rad);
        pr.nabla_ad.push(s.eos.nabla_ad);
        pr.nabla.push(s.nabla());
        pr.convective.push(s.mlt.convective);
        pr.v_conv.push(s.mlt.v_conv);
        pr.f_conv.push(s.mlt.f_conv);
        pr.hp.push(s.hp);
        pr.mu.push(s.eos.mu);
        pr.beta.push(s.eos.beta);
        pr.gamma1.push(s.eos.gamma1);
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
    // Energy generated, with the solver's own midpoint rule: this verifies
    // that the converged solution conserves energy globally. (A point-based
    // quadrature would mostly measure how well thin burning shells are sampled.)
    let (mut e_nuc_mid, mut e_grav) = (0.0, 0.0);
    for k in 0..n - 1 {
        let (s, eg) = problem.interval_state(k, &y[k * NVAR..(k + 1) * NVAR], &y[(k + 1) * NVAR..(k + 2) * NVAR]);
        let dm = problem.m[k + 1] - problem.m[k];
        e_nuc_mid += dm * s.eps();
        e_grav += dm * eg;
    }
    let thermal = trap(&|k| 3.0 * pr.p[k] / pr.rho[k]);
    let grav = trap(&|k| G * pr.m[k] / pr.r[k]);

    // Hydrogen-exhausted core: contiguous points from the centre with X < 10⁻³.
    let he_core_q = pr.x.iter().zip(&pr.q).take_while(|(x, _)| **x < 1e-3).last().map(|(_, q)| *q).unwrap_or(0.0);
    let mean_density = mass / (FOUR_PI / 3.0 * radius.powi(3));
    let g_surf = G * mass / (radius * radius);
    let g_sun = G * M_SUN / (R_SUN * R_SUN);
    let rho_sun = M_SUN / (FOUR_PI / 3.0 * R_SUN.powi(3));

    let globals = Globals {
        age_yr: 0.0,
        phase: Phase::Zams,
        x_c: pr.x[0],
        y_c: pr.y[0],
        he_core_q,
        mean_density,
        nu_max_uhz: 3090.0 * (g_surf / g_sun) * (pr.t[n - 1] / TEFF_SUN).powf(-0.5),
        delta_nu_uhz: 135.1 * (mean_density / rho_sun).sqrt(),
        frac_grav: e_grav / lum,
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
        energy_ratio: (e_nuc_mid + e_grav + pr.l[0]) / lum,
        newton_iterations,
        max_correction,
        residuals,
    };

    StarModel {
        params: params.clone(),
        physics: problem.physics.describe(),
        globals,
        zones,
        profile: pr,
        solution: y,
        l_scale: problem.l_scale,
        pulsation: None,
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
