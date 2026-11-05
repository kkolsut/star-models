//! Solve driver: cold start from a polytrope, warm start / continuation from a
//! previously converged model, with step halving when Newton fails.

use super::equations::{I_LNP, I_LNR, I_LNT, NVAR, StructureProblem};
use super::initial::polytrope_guess;
use super::mesh::{make_mesh, remap};
use super::model::{StarModel, StarParams, build_model, sanity_check};
use crate::physics::Physics;
use crate::solver::{HenyeyNewton, NewtonOptions, StepOutcome};

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum SolveError {
    #[error("invalid parameters: {0}")]
    InvalidParams(String),
    #[error("solver did not converge: {0}")]
    NoConvergence(String),
    #[error("cancelled")]
    Cancelled,
}

/// Progress events emitted while solving.
#[derive(Debug, Clone, PartialEq)]
pub enum SolveEvent {
    Stage(String),
    Iteration { iter: usize, max_correction: f64, residual: f64, damping: f64 },
}

/// Hook for progress reporting and cancellation.
pub trait SolveControl {
    fn cancelled(&self) -> bool {
        false
    }
    fn report(&mut self, _event: SolveEvent) {}
}

/// A [`SolveControl`] that never cancels and ignores progress.
pub struct NoControl;
impl SolveControl for NoControl {}

const COLD_START: NewtonOptions = NewtonOptions { max_iter: 250, tol: 1e-8 };
const WARM_START: NewtonOptions = NewtonOptions { max_iter: 60, tol: 1e-8 };
/// Smallest continuation step (as a fraction of the full path).
const MIN_STEP: f64 = 1.0 / 1024.0;

/// Compute a ZAMS model for `params`.
///
/// If `warm` is given, the solver continues from that converged model, which
/// is much faster and more robust for small parameter changes. Otherwise (or
/// if that fails) it starts from an n = 3 polytrope. If that also fails, it
/// solves a 1 M☉ reference model and continues from there.
pub fn solve_model(
    params: &StarParams,
    physics: &Physics,
    warm: Option<&StarModel>,
    ctl: &mut dyn SolveControl,
) -> Result<StarModel, SolveError> {
    params.validate().map_err(SolveError::InvalidParams)?;

    if let Some(w) = warm {
        match continuation(w, params, physics, ctl) {
            Ok(m) => return Ok(m),
            Err(SolveError::Cancelled) => return Err(SolveError::Cancelled),
            Err(_) => {} // fall back to a cold start
        }
    }

    ctl.report(SolveEvent::Stage(format!("cold start from n=3 polytrope (M = {:.3} M☉)", params.mass)));
    let mut last_err = match cold_start(params, physics, ctl) {
        Ok(m) => return Ok(m),
        Err(SolveError::Cancelled) => return Err(SolveError::Cancelled),
        Err(e) => e,
    };

    // Fall back to a 1 M☉ reference model (first with the requested
    // composition, then the default solar one) and continue from it.
    let mut references = vec![StarParams { mass: 1.0, ..params.clone() }];
    if params.x != StarParams::sun().x || params.z != StarParams::sun().z {
        references.push(StarParams { n_points: params.n_points, ..StarParams::sun() });
    }
    for reference in references {
        if reference == *params {
            continue;
        }
        ctl.report(SolveEvent::Stage(format!(
            "fallback: reference model M = 1 M☉, X = {}, Z = {}",
            reference.x, reference.z
        )));
        match cold_start(&reference, physics, ctl) {
            Ok(ref_model) => match continuation(&ref_model, params, physics, ctl) {
                Ok(m) => return Ok(m),
                Err(SolveError::Cancelled) => return Err(SolveError::Cancelled),
                Err(e) => last_err = e,
            },
            Err(SolveError::Cancelled) => return Err(SolveError::Cancelled),
            Err(e) => last_err = e,
        }
    }
    Err(last_err)
}

fn cold_start(params: &StarParams, physics: &Physics, ctl: &mut dyn SolveControl) -> Result<StarModel, SolveError> {
    let comp = params.composition().map_err(|e| SolveError::InvalidParams(e.to_string()))?;
    let q = make_mesh(params.n_points);
    let guess = polytrope_guess(physics, params.mass, &comp, &q);
    newton(params, physics, &q, guess.y, guess.l_scale, COLD_START, ctl)
}

/// Run Newton iterations for fixed parameters from initial guess `y0`.
pub fn newton(
    params: &StarParams,
    physics: &Physics,
    q: &[f64],
    y0: Vec<f64>,
    l_scale: f64,
    opts: NewtonOptions,
    ctl: &mut dyn SolveControl,
) -> Result<StarModel, SolveError> {
    let comp = params.composition().map_err(|e| SolveError::InvalidParams(e.to_string()))?;
    let mass = params.mass * crate::constants::M_SUN;
    let problem = StructureProblem::new(physics, mass, q, comp, params.alpha_mlt, l_scale);
    let mut nr = HenyeyNewton::new(&problem, y0, opts);
    loop {
        if ctl.cancelled() {
            return Err(SolveError::Cancelled);
        }
        match nr.step() {
            StepOutcome::Continue(info) => ctl.report(SolveEvent::Iteration {
                iter: info.iter,
                max_correction: info.max_correction,
                residual: info.residual,
                damping: info.damping,
            }),
            StepOutcome::Converged(info) => {
                ctl.report(SolveEvent::Iteration {
                    iter: info.iter,
                    max_correction: info.max_correction,
                    residual: info.residual,
                    damping: info.damping,
                });
                let iters = nr.iterations();
                let model = build_model(params, &problem, q, nr.into_solution(), iters, info.max_correction);
                sanity_check(&model).map_err(SolveError::NoConvergence)?;
                return Ok(model);
            }
            StepOutcome::Failed(msg) => return Err(SolveError::NoConvergence(msg)),
        }
    }
}

/// Homology-scaled guess for `target` built from the converged `model`.
fn homology_guess(model: &StarModel, target: &StarParams) -> (Vec<f64>, Vec<f64>, f64) {
    let q_new = make_mesh(target.n_points);
    let mut y = if model.profile.q.len() == q_new.len() {
        model.solution.clone()
    } else {
        remap(&model.profile.q, &model.solution, NVAR, &q_new)
    };
    // R ∝ M^a, P ∝ M²/R⁴, T ∝ M/R, L ∝ M^b. The unknown L/L_scale is kept;
    // only the luminosity scale changes.
    let (a, b) = (0.75, 3.5);
    let d = (target.mass / model.params.mass).ln();
    for k in 0..q_new.len() {
        y[k * NVAR + I_LNR] += a * d;
        y[k * NVAR + I_LNP] += (2.0 - 4.0 * a) * d;
        y[k * NVAR + I_LNT] += (1.0 - a) * d;
    }
    (q_new, y, model.l_scale * (b * d).exp())
}

/// Continue from a converged model to the target parameters. The full step is
/// tried first; on failure the step is halved repeatedly.
fn continuation(
    from: &StarModel,
    target: &StarParams,
    physics: &Physics,
    ctl: &mut dyn SolveControl,
) -> Result<StarModel, SolveError> {
    let start = from.params.clone();
    if start == *target && !from.solution.is_empty() {
        return Ok(from.clone());
    }
    let mut current = from.clone();
    let (mut t, mut dt) = (0.0_f64, 1.0_f64);
    while t < 1.0 {
        let t_try = (t + dt).min(1.0);
        let p = start.interpolate(target, t_try);
        ctl.report(SolveEvent::Stage(format!(
            "continuation {:>3.0}% (M = {:.3} M☉, X = {:.3}, Z = {:.4})",
            t_try * 100.0,
            p.mass,
            p.x,
            p.z
        )));
        let (q, y0, l_scale) = homology_guess(&current, &p);
        match newton(&p, physics, &q, y0, l_scale, WARM_START, ctl) {
            Ok(m) => {
                current = m;
                t = t_try;
                dt = (dt * 2.0).min(1.0);
            }
            Err(SolveError::Cancelled) => return Err(SolveError::Cancelled),
            Err(e) => {
                dt *= 0.5;
                if dt < MIN_STEP {
                    return Err(e);
                }
            }
        }
    }
    // The last step was made with the interpolated parameters; make sure the
    // stored parameters are exactly the target.
    current.params = target.clone();
    Ok(current)
}
