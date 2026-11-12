//! Stellar evolution: a sequence of structure models in time.
//!
//! Each time step Δt:
//! 1. Composition: each convective zone is mixed instantly to its
//!    mass-averaged composition (a convective core is first extended by step
//!    overshoot, r_ov = r_cc + α_ov min(H_P(r_cc), r_cc); mixing only, the
//!    temperature gradient there stays radiative). Then hydrogen burns to
//!    helium, X ← X − (ε_pp/Q_pp^eff + ε_CNO/Q_CNO) Δt, and helium to carbon,
//!    Y ← Y − ε_3α Δt / Q_He. A mixed zone burns as one unit at its
//!    mass-averaged rate, limited to what it contains.
//! 2. Structure: the four structure equations are solved for the new
//!    composition, with the gravitational energy term
//!    ε_grav = −c_P T (Δ ln T − ∇_ad Δ ln P) / Δt evaluated implicitly
//!    against the previous model on the same mass mesh (Henyey method,
//!    warm-started from the previous model).
//! 3. Time-step control: Δt is limited by the hydrogen and helium consumption
//!    (ΔX, ΔY ≤ dx_max per shell, and ≤ 30% of the local X, Y) and by the
//!    structural change per step (target `dlog_target` in ln r, ln P, ln T and
//!    ln L). A rejected or failed step is retried with a smaller Δt.
//!
//! Along the run the totals of the total-energy check are accumulated
//! ([`RunTotals`]): ∫L dt, ∫∫ε_nuc dm dt and the composition term omitted from
//! ε_grav; with Δ(E_int + E_grav) they close to the discretisation error.

use super::equations::{I_L, I_LNP, I_LNR, I_LNT, NVAR, PrevStep};
use super::model::{ConvectiveZone, Phase, Profile, RunTotals, StarModel};
use super::solve::{SolveControl, SolveError, solve_structure};
use crate::composition::Composition;
use crate::constants::{L_SUN, YEAR};
use crate::physics::Physics;
pub use crate::physics::nuclear::{Q_CNO, Q_HE, Q_PP_I, q_pp_eff};
use crate::pulsation;
use crate::solver::NewtonOptions;

/// When an evolution run should stop.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StopCondition {
    /// Stop when the central hydrogen is exhausted (X_c < 10⁻³).
    EndOfMainSequence,
    /// Stop at the given age [yr].
    Age(f64),
    /// Continue until the solver can no longer follow the star, or a safety limit is reached.
    AsFarAsPossible,
}

#[derive(Clone, Copy, Debug)]
pub struct EvolutionOptions {
    pub stop: StopCondition,
    /// Largest change of the hydrogen mass fraction in any shell per step.
    pub dx_max: f64,
    /// Target largest change of ln r, ln P, ln T and ln L per step.
    pub dlog_target: f64,
    /// A step whose largest change exceeds `reject_factor × dlog_target` is
    /// rejected and retried with half the time step.
    pub reject_factor: f64,
    /// Safety limits.
    pub max_steps: usize,
    pub max_age_yr: f64,
    /// Compute the fundamental radial pulsation mode of every model.
    pub pulsation: bool,
}

impl Default for EvolutionOptions {
    fn default() -> Self {
        Self {
            stop: StopCondition::EndOfMainSequence,
            dx_max: 0.01,
            dlog_target: 0.02,
            reject_factor: 4.0,
            max_steps: 3000,
            max_age_yr: 3.0e10,
            pulsation: true,
        }
    }
}

/// Result of one call to [`Evolution::step`].
#[derive(Clone, Debug, PartialEq)]
pub enum EvolutionStep {
    /// A new model was computed (see [`Evolution::model`]).
    Advanced { dt_yr: f64, attempts: usize },
    /// The run is over, with the reason.
    Finished(String),
}

/// Central hydrogen fraction below which hydrogen counts as exhausted.
pub const X_EXHAUSTED: f64 = 1e-3;
/// Smallest time step tried before giving up [yr].
const MIN_DT_YR: f64 = 1.0;

pub struct Evolution {
    /// Latest converged model.
    pub model: StarModel,
    /// Time step to try next [yr].
    pub dt_yr: f64,
    pub steps: usize,
    pub opts: EvolutionOptions,
    r_zams: f64,
    finished: Option<String>,
}

impl Evolution {
    /// Start evolving from a converged ZAMS model.
    pub fn new(mut zams: StarModel, opts: EvolutionOptions) -> Self {
        zams.globals.age_yr = 0.0;
        zams.globals.phase = Phase::Zams;
        let r_zams = zams.globals.radius;
        zams.globals.run = Some(RunTotals::start(&zams.globals));
        let dt = (0.1 * burn_limit_s(&zams, opts.dx_max) / YEAR).clamp(MIN_DT_YR, 1e6);
        Self { model: zams, dt_yr: dt, steps: 0, opts, r_zams, finished: None }
    }

    /// Continue an earlier run from its latest model.
    pub fn resume(mut model: StarModel, dt_yr: f64, r_zams: f64, opts: EvolutionOptions) -> Self {
        if model.globals.run.is_none() {
            model.globals.run = Some(RunTotals::start(&model.globals));
        }
        Self { model, dt_yr: dt_yr.max(MIN_DT_YR), steps: 0, opts, r_zams, finished: None }
    }

    pub fn r_zams(&self) -> f64 {
        self.r_zams
    }

    fn stop_reason(&self) -> Option<String> {
        let g = &self.model.globals;
        if self.steps >= self.opts.max_steps {
            return Some(format!("reached the limit of {} steps", self.opts.max_steps));
        }
        if g.age_yr >= self.opts.max_age_yr * (1.0 - 1e-9) {
            return Some(format!("reached the age limit of {:.1} Gyr", self.opts.max_age_yr / 1e9));
        }
        match self.opts.stop {
            StopCondition::EndOfMainSequence if g.x_c < X_EXHAUSTED => {
                Some("central hydrogen exhausted: end of the main sequence".into())
            }
            StopCondition::Age(a) if g.age_yr >= a * (1.0 - 1e-9) => {
                Some(format!("reached the target age of {}", fmt_age(a)))
            }
            StopCondition::AsFarAsPossible if g.teff < 2600.0 => {
                Some("T_eff dropped below 2600 K, outside the range of the analytic opacities".into())
            }
            StopCondition::AsFarAsPossible if g.tc > 3e8 => Some("central temperature above 3 × 10⁸ K".into()),
            _ => None,
        }
    }

    /// Advance by one time step (with retries at smaller Δt if needed).
    pub fn step(&mut self, physics: &Physics, ctl: &mut dyn SolveControl) -> Result<EvolutionStep, SolveError> {
        if let Some(r) = &self.finished {
            return Ok(EvolutionStep::Finished(r.clone()));
        }
        if let Some(r) = self.stop_reason() {
            self.finished = Some(r.clone());
            return Ok(EvolutionStep::Finished(r));
        }
        let old = &self.model;
        let age = old.globals.age_yr;
        let mut dt = self.dt_yr.min(burn_limit_s(old, self.opts.dx_max) / YEAR).max(MIN_DT_YR);
        let target_age = match self.opts.stop {
            StopCondition::Age(a) => a.min(self.opts.max_age_yr),
            _ => self.opts.max_age_yr,
        };
        dt = dt.min(target_age - age).max(1e-6);

        // Re-normalise the luminosity unknown to the current luminosity.
        let l_scale = old.globals.luminosity * L_SUN;
        let mut y0 = old.solution.clone();
        for k in 0..old.profile.len() {
            y0[k * NVAR + I_L] *= old.l_scale / l_scale;
        }
        let ln_p: Vec<f64> = (0..old.profile.len()).map(|k| old.solution[k * NVAR + I_LNP]).collect();
        let ln_t: Vec<f64> = (0..old.profile.len()).map(|k| old.solution[k * NVAR + I_LNT]).collect();

        let mut attempts = 0;
        loop {
            attempts += 1;
            let comps = advance_composition(old, dt * YEAR);
            let comps_used = comps.clone();
            let prev = PrevStep { ln_p: ln_p.clone(), ln_t: ln_t.clone(), dt: dt * YEAR };
            let opts = NewtonOptions { max_iter: 40, tol: 1e-8 };
            match solve_structure(
                &old.params,
                physics,
                &old.profile.q,
                comps,
                Some(prev),
                y0.clone(),
                l_scale,
                opts,
                ctl,
            ) {
                Ok(mut m) => {
                    let change = max_change(old, &m);
                    if change > self.opts.reject_factor * self.opts.dlog_target && dt > MIN_DT_YR && attempts < 12 {
                        dt *= 0.5;
                        continue;
                    }
                    m.globals.age_yr = age + dt;
                    m.globals.phase = classify(&m, self.r_zams);
                    m.globals.run = Some(accumulate_run(old, &m, &comps_used, physics, dt * YEAR));
                    if self.opts.pulsation {
                        m.pulsation = pulsation::fundamental_mode(&m, old.pulsation.as_ref()).ok();
                    }
                    let grow = (self.opts.dlog_target / change.max(1e-9)).clamp(0.5, 2.0);
                    self.dt_yr = dt * grow;
                    self.steps += 1;
                    self.model = m;
                    return Ok(EvolutionStep::Advanced { dt_yr: dt, attempts });
                }
                Err(SolveError::Cancelled) => return Err(SolveError::Cancelled),
                Err(e) => {
                    dt /= 3.0;
                    if dt < MIN_DT_YR || attempts >= 12 {
                        let reason = format!("the solver could not continue ({e})");
                        self.finished = Some(reason.clone());
                        return Ok(EvolutionStep::Finished(reason));
                    }
                }
            }
        }
    }
}

/// Human-readable age.
pub fn fmt_age(yr: f64) -> String {
    if yr >= 1e9 {
        format!("{:.3} Gyr", yr / 1e9)
    } else if yr >= 1e6 {
        format!("{:.2} Myr", yr / 1e6)
    } else if yr >= 1e3 {
        format!("{:.1} kyr", yr / 1e3)
    } else {
        format!("{yr:.0} yr")
    }
}

/// Mass associated with each mesh point (half of each neighbouring interval).
fn point_masses(m: &[f64]) -> Vec<f64> {
    let n = m.len();
    (0..n)
        .map(|k| {
            let lo = if k == 0 { 0.0 } else { 0.5 * (m[k] + m[k - 1]) };
            let hi = if k + 1 == n { m[k] } else { 0.5 * (m[k] + m[k + 1]) };
            (hi - lo).max(0.0)
        })
        .collect()
}

/// Consumption rates per point [s⁻¹]: hydrogen, (ε_pp/Q_pp^eff + ε_CNO/Q_CNO),
/// and helium by the triple-α reaction, ε_3α/Q_He. Each is averaged over every
/// mixed zone, which burns as one unit.
pub fn burning_rates(model: &StarModel) -> (Vec<f64>, Vec<f64>) {
    let p = &model.profile;
    let mut h: Vec<f64> = (0..p.len()).map(|k| p.eps_pp[k] / q_pp_eff(p.phi_at(k)) + p.eps_cno[k] / Q_CNO).collect();
    let mut he: Vec<f64> = p.eps_3a.iter().map(|e| e / Q_HE).collect();
    let dm = point_masses(&p.m);
    for (a, b) in mixed_runs(model) {
        let mass: f64 = dm[a..=b].iter().sum::<f64>().max(f64::MIN_POSITIVE);
        for rate in [&mut h, &mut he] {
            let avg = (a..=b).map(|k| rate[k] * dm[k]).sum::<f64>() / mass;
            rate[a..=b].iter_mut().for_each(|r| *r = avg);
        }
    }
    (h, he)
}

/// Largest change of a mass fraction allowed in one step: min(dx_max, max(0.3 v, 10⁻⁵)).
fn fraction_limit(v: f64, dx_max: f64) -> f64 {
    dx_max.min((0.3 * v).max(1e-5))
}

/// Largest time step [s] keeping the hydrogen and helium changes per shell within limits.
fn burn_limit_s(model: &StarModel, dx_max: f64) -> f64 {
    let (h, he) = burning_rates(model);
    let p = &model.profile;
    let limit = |rates: &[f64], v: &[f64]| {
        rates
            .iter()
            .zip(v)
            .filter(|(rate, v)| **rate > 0.0 && **v > 1e-8)
            .map(|(rate, v)| fraction_limit(*v, dx_max) / rate)
            .fold(f64::INFINITY, f64::min)
    };
    limit(&h, &p.x).min(limit(&he, &p.y))
}

/// ∫X dm with the point masses of the composition step [g].
pub fn hydrogen_mass(model: &StarModel) -> f64 {
    let p = &model.profile;
    point_masses(&p.m).iter().zip(&p.x).map(|(dm, x)| dm * x).sum()
}

/// Effective energy per gram of hydrogen burnt in a model [erg g⁻¹]: the
/// hydrogen-burning power over the hydrogen consumption rate.
pub fn q_h_eff(model: &StarModel) -> f64 {
    let p = &model.profile;
    let dm = point_masses(&p.m);
    let (h, _) = burning_rates(model);
    let power: f64 = (0..p.len()).map(|k| dm[k] * (p.eps_pp[k] + p.eps_cno[k])).sum();
    let rate: f64 = (0..p.len()).map(|k| dm[k] * h[k]).sum();
    if rate > 0.0 { power / rate } else { Q_PP_I }
}

/// Totals of the run after the step from `old` to `new` (time step `dt_s`,
/// composition `comps` given to the structure solver).
fn accumulate_run(old: &StarModel, new: &StarModel, comps: &[Composition], physics: &Physics, dt_s: f64) -> RunTotals {
    let mut run = old.globals.run.unwrap_or_else(|| RunTotals::start(&old.globals));
    run.l_dt += new.globals.luminosity * L_SUN * dt_s;
    run.nuc_dt += new.globals.e_nuc * dt_s;
    run.fuel += (hydrogen_mass(old) - hydrogen_mass(new)) * q_h_eff(old);
    // Composition term of ε_grav, −(∂h/∂X_i)_{P,T} dX_i/dt, evaluated at the
    // previous (P, T): Σ Δm [h(P, T, X_new) − h(P, T, X_old)], h = u + P/ρ.
    let p = &old.profile;
    let dm = point_masses(&p.m);
    let h = |k: usize, c: &Composition| {
        let s = physics.eos.eval(p.p[k], p.t[k], c);
        s.u + p.p[k] / s.rho
    };
    run.comp_term += (0..p.len())
        .map(|k| {
            let c_old = Composition { x: p.x[k], y: p.y[k], z: p.z[k], x_cno: comps[k].x_cno };
            dm[k] * (h(k, &comps[k]) - h(k, &c_old))
        })
        .sum::<f64>();
    run.max_screening = run.max_screening.max(new.globals.screening.max_exponent);
    run
}

/// Index ranges [a, b] of consecutive convective points.
fn convective_runs(conv: &[bool]) -> Vec<(usize, usize)> {
    let mut runs = Vec::new();
    let mut start = None;
    for (k, &c) in conv.iter().enumerate() {
        match (c, start) {
            (true, None) => start = Some(k),
            (false, Some(a)) => {
                runs.push((a, k - 1));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(a) = start {
        runs.push((a, conv.len() - 1));
    }
    runs
}

/// Extent of the step-overshoot region above a convective core.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoreOvershoot {
    /// Schwarzschild boundary of the core [cm].
    pub r_cc: f64,
    /// Outer edge of the mixed region, r_cc + α_ov min(H_P(r_cc), r_cc) [cm].
    pub r_ov: f64,
    /// m/M at `r_ov` (linear interpolation in r).
    pub q_ov: f64,
    /// Last mesh point with r ≤ r_ov: the mixed run is 0..=k_ov.
    pub k_ov: usize,
}

/// Step overshoot of a convective core that starts at the centre, for
/// α_ov > 0. `None` without a convective core or with α_ov = 0.
///
/// H_P(r_cc) is interpolated linearly in r between the two mesh points around
/// the boundary; min(H_P, r_cc) keeps the extension finite for small cores.
pub fn core_overshoot(p: &Profile, zones: &[ConvectiveZone], alpha_ov: f64) -> Option<CoreOvershoot> {
    if !(alpha_ov > 0.0) || p.is_empty() || !p.convective[0] {
        return None;
    }
    let n = p.len();
    let radius = p.r[n - 1];
    let core = zones.first().filter(|z| z.q_inner <= p.q[0])?;
    let r_cc = core.r_outer * radius;
    let at = |v: &[f64], r: f64| -> f64 {
        let i = p.r.partition_point(|&x| x < r).clamp(1, n - 1);
        let w = ((r - p.r[i - 1]) / (p.r[i] - p.r[i - 1])).clamp(0.0, 1.0);
        v[i - 1] + w * (v[i] - v[i - 1])
    };
    let hp = at(&p.hp, r_cc);
    let r_ov = (r_cc + alpha_ov * hp.min(r_cc)).min(radius);
    let last_conv = p.convective.iter().take_while(|&&c| c).count() - 1;
    let k_ov = p.r.partition_point(|&r| r <= r_ov).saturating_sub(1).max(last_conv);
    Some(CoreOvershoot { r_cc, r_ov, q_ov: at(&p.q, r_ov), k_ov })
}

/// Index ranges [a, b] mixed as one unit: the convective runs, with a
/// convective core extended by step overshoot (runs it reaches are merged).
fn mixed_runs(model: &StarModel) -> Vec<(usize, usize)> {
    let p = &model.profile;
    let mut runs = convective_runs(&p.convective);
    if let Some(ov) = core_overshoot(p, &model.zones, model.params.alpha_ov) {
        let mut b = runs[0].1.max(ov.k_ov);
        let mut i = 1;
        while i < runs.len() && runs[i].0 <= b + 1 {
            b = b.max(runs[i].1);
            i += 1;
        }
        runs.splice(0..i, [(0, b)]);
    }
    runs
}

/// Mix every mixed zone of `model` (a convective core including its overshoot
/// region) to its mass average, then burn for Δt [s]: a mixed zone burns as one
/// unit at its mass-averaged rate, and no point burns more than it contains.
pub fn advance_composition(model: &StarModel, dt_s: f64) -> Vec<Composition> {
    let p = &model.profile;
    let n = p.len();
    let (mut x, mut y, mut z) = (p.x.clone(), p.y.clone(), p.z.clone());
    let dm = point_masses(&p.m);
    for (a, b) in mixed_runs(model) {
        let mass: f64 = dm[a..=b].iter().sum::<f64>().max(f64::MIN_POSITIVE);
        let avg = |v: &[f64]| (a..=b).map(|k| v[k] * dm[k]).sum::<f64>() / mass;
        let (xa, ya, za) = (avg(&x), avg(&y), avg(&z));
        for k in a..=b {
            x[k] = xa;
            y[k] = ya;
            z[k] = za;
        }
    }
    let (h_rate, he_rate) = burning_rates(model);
    for k in 0..n {
        let burnt_h = (h_rate[k] * dt_s).min(x[k]).max(0.0);
        x[k] -= burnt_h;
        y[k] += burnt_h;
        let burnt_he = (he_rate[k] * dt_s).min(y[k]).max(0.0);
        y[k] -= burnt_he;
        z[k] += burnt_he;
    }
    // CNO catalysts keep their initial abundance (He-burning products are not CNO catalysts for H burning here).
    let x_cno = model.params.f_cno * model.params.z;
    (0..n)
        .map(|k| {
            let s = x[k] + y[k] + z[k];
            Composition { x: x[k] / s, y: y[k] / s, z: z[k] / s, x_cno }
        })
        .collect()
}

/// Largest change of ln r, ln P, ln T over the star, and of ln L at the surface.
fn max_change(a: &StarModel, b: &StarModel) -> f64 {
    let mut d: f64 = 0.0;
    for k in 0..a.profile.len() {
        for j in [I_LNR, I_LNP, I_LNT] {
            d = d.max((b.solution[k * NVAR + j] - a.solution[k * NVAR + j]).abs());
        }
    }
    d.max((b.globals.luminosity / a.globals.luminosity).ln().abs())
}

/// Evolutionary phase from the model's state.
pub fn classify(m: &StarModel, r_zams: f64) -> Phase {
    let g = &m.globals;
    if g.frac_3a > 0.1 {
        Phase::CoreHelium
    } else if g.x_c >= X_EXHAUSTED {
        if g.age_yr > 0.0 { Phase::MainSequence } else { Phase::Zams }
    } else if g.radius > 5.0 * r_zams && g.teff < 5000.0 {
        Phase::RedGiant
    } else {
        Phase::ShellHydrogen
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{NoControl, StarParams, solve_model};

    fn sun() -> StarModel {
        solve_model(&StarParams { n_points: 300, ..StarParams::sun() }, &Physics::analytic(), None, &mut NoControl)
            .unwrap()
    }

    /// Plan 1.2: wherever ε_3α > 0 the helium consumption limits the step.
    #[test]
    fn helium_burning_gets_a_finite_limit() {
        let mut m = sun();
        let n = m.profile.len();
        m.profile.eps_pp = vec![0.0; n];
        m.profile.eps_cno = vec![0.0; n];
        m.profile.eps_3a = vec![0.0; n];
        assert_eq!(burn_limit_s(&m, 0.01), f64::INFINITY);
        m.profile.eps_3a[0] = 1e3;
        let limit = burn_limit_s(&m, 0.01);
        assert!(limit.is_finite() && limit > 0.0);
        // ΔY = min(0.01, 0.3 Y) at the rate ε_3α / Q_He.
        let expected = 0.01 / (1e3 / Q_HE);
        assert!((limit / expected - 1.0).abs() < 1e-12, "{limit} vs {expected}");
    }

    /// Plan 1.1: a mixed zone burns at its averaged rate and stays uniform; the
    /// hydrogen removed equals the averaged rate times Δt.
    #[test]
    fn mixed_zone_burns_as_one_unit() {
        let p = StarParams { mass: 3.0, n_points: 300, ..StarParams::sun() };
        let m = solve_model(&p, &Physics::analytic(), None, &mut NoControl).unwrap();
        let runs = mixed_runs(&m);
        let (a, b) = runs[0];
        assert_eq!(a, 0);
        let (h, _) = burning_rates(&m);
        let dt = 0.5 * m.profile.x[0] / h[0];
        let comps = advance_composition(&m, dt);
        for k in a..=b {
            assert!((comps[k].x - comps[0].x).abs() < 1e-15);
        }
        assert!((comps[0].x - (m.profile.x[0] - h[0] * dt)).abs() < 1e-12);
        // A step that would burn more than the zone holds empties it, uniformly.
        let comps = advance_composition(&m, 3.0 * dt);
        assert!((a..=b).all(|k| comps[k].x == 0.0));
    }
}
