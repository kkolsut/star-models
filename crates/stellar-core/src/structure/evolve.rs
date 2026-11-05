//! Stellar evolution: a sequence of structure models in time.
//!
//! Each time step Δt:
//! 1. Composition: hydrogen burns to helium, X ← X − ε_H Δt / Q_H, and helium
//!    to carbon, Y ← Y − ε_3α Δt / Q_He. Each convective zone is then mixed
//!    instantly to its mass-averaged composition.
//! 2. Structure: the four structure equations are solved for the new
//!    composition, with the gravitational energy term
//!    ε_grav = −c_P T (Δ ln T − ∇_ad Δ ln P) / Δt evaluated implicitly
//!    against the previous model on the same mass mesh (Henyey method,
//!    warm-started from the previous model).
//! 3. Time-step control: Δt is limited by the hydrogen consumption
//!    (ΔX ≤ dx_max per shell, and ≤ 30% of the local X) and by the structural
//!    change per step (target `dlog_target` in ln r, ln P, ln T and ln L).
//!    A rejected or failed step is retried with a smaller Δt.

use super::equations::{I_L, I_LNP, I_LNR, I_LNT, NVAR, PrevStep};
use super::model::{Phase, StarModel};
use super::solve::{SolveControl, SolveError, solve_structure};
use crate::composition::Composition;
use crate::constants::{L_SUN, YEAR};
use crate::physics::Physics;
use crate::pulsation;
use crate::solver::NewtonOptions;

/// Energy deposited per gram of hydrogen converted to helium [erg g⁻¹]
/// (26.2 MeV per ⁴He nucleus, net of neutrino losses).
pub const Q_H: f64 = 6.27e18;
/// Energy released per gram of helium burnt to carbon by the triple-α
/// reaction [erg g⁻¹] (7.275 MeV per ¹²C nucleus).
pub const Q_HE: f64 = 5.85e17;

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
        let dt = (0.1 * burn_limit_s(&zams, opts.dx_max) / YEAR).clamp(MIN_DT_YR, 1e6);
        Self { model: zams, dt_yr: dt, steps: 0, opts, r_zams, finished: None }
    }

    /// Continue an earlier run from its latest model.
    pub fn resume(model: StarModel, dt_yr: f64, r_zams: f64, opts: EvolutionOptions) -> Self {
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

/// Rate of change of X per point [s⁻¹], averaged over each convective zone
/// (which is mixed, so it burns as one unit).
fn hydrogen_rates(model: &StarModel) -> Vec<f64> {
    let p = &model.profile;
    let mut rate: Vec<f64> = (0..p.len()).map(|k| (p.eps_pp[k] + p.eps_cno[k]) / Q_H).collect();
    let dm = point_masses(&p.m);
    for (a, b) in convective_runs(&p.convective) {
        let mass: f64 = dm[a..=b].iter().sum();
        let avg = (a..=b).map(|k| rate[k] * dm[k]).sum::<f64>() / mass.max(f64::MIN_POSITIVE);
        rate[a..=b].iter_mut().for_each(|r| *r = avg);
    }
    rate
}

/// Largest time step [s] keeping the hydrogen change per shell within limits.
fn burn_limit_s(model: &StarModel, dx_max: f64) -> f64 {
    let rates = hydrogen_rates(model);
    rates
        .iter()
        .zip(&model.profile.x)
        .filter(|(rate, x)| **rate > 0.0 && **x > 1e-8)
        .map(|(rate, x)| dx_max.min((0.3 * x).max(1e-5)) / rate)
        .fold(f64::INFINITY, f64::min)
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

/// Burn for Δt [s], then mix every convective zone of `model`.
pub fn advance_composition(model: &StarModel, dt_s: f64) -> Vec<Composition> {
    let p = &model.profile;
    let n = p.len();
    let (mut x, mut y, mut z) = (p.x.clone(), p.y.clone(), p.z.clone());
    for k in 0..n {
        let burnt_h = ((p.eps_pp[k] + p.eps_cno[k]) * dt_s / Q_H).min(x[k]).max(0.0);
        x[k] -= burnt_h;
        y[k] += burnt_h;
        let burnt_he = (p.eps_3a[k] * dt_s / Q_HE).min(y[k]).max(0.0);
        y[k] -= burnt_he;
        z[k] += burnt_he;
    }
    let dm = point_masses(&p.m);
    for (a, b) in convective_runs(&p.convective) {
        let mass: f64 = dm[a..=b].iter().sum::<f64>().max(f64::MIN_POSITIVE);
        let avg = |v: &[f64]| (a..=b).map(|k| v[k] * dm[k]).sum::<f64>() / mass;
        let (xa, ya, za) = (avg(&x), avg(&y), avg(&z));
        for k in a..=b {
            x[k] = xa;
            y[k] = ya;
            z[k] = za;
        }
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
