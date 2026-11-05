//! Generic Henyey-type relaxation solver for two-point boundary-value problems.
//!
//! The problem is discretised on `npts` mesh points with `nvar` unknowns per
//! point. The equations are ordered as
//!
//! ```text
//! [ inner BC (n_inner rows) | interval 0 (nvar rows) | … | interval npts−2 | outer BC (nvar − n_inner rows) ]
//! ```
//!
//! Interval k couples only points k and k+1, so the Jacobian is banded with
//! kl = n_inner + nvar − 1 and ku = 2 nvar − 1 − n_inner. Each Newton step
//! builds the Jacobian by local central finite differences (no analytic
//! derivatives are needed), solves it with a banded LU, and applies a damped
//! correction.
//!
//! [`HenyeyNewton::step`] performs one Newton iteration at a time. Native
//! callers loop over it on a worker thread; a future WASM build can run a few
//! steps per frame.

use super::banded::BandedMatrix;

/// A discretised two-point boundary-value problem.
pub trait BvpProblem {
    /// Unknowns per mesh point.
    fn nvar(&self) -> usize;
    /// Number of boundary conditions at the inner (first) point.
    fn n_inner(&self) -> usize;
    /// Number of mesh points.
    fn npts(&self) -> usize;
    /// Inner boundary residuals (`n_inner` values) from the first point.
    fn inner_bc(&self, y0: &[f64], out: &mut [f64]);
    /// Outer boundary residuals (`nvar − n_inner` values) from the last point.
    fn outer_bc(&self, yn: &[f64], out: &mut [f64]);
    /// Difference-equation residuals (`nvar` values) between points k and k+1.
    fn interval(&self, k: usize, yk: &[f64], yk1: &[f64], out: &mut [f64]);
    /// Largest allowed absolute correction of variable `var` per iteration
    /// (used for damping).
    fn max_correction(&self, _var: usize) -> f64 {
        f64::INFINITY
    }
    /// Finite-difference step for variable `var` at its current value.
    fn fd_step(&self, _var: usize, value: f64) -> f64 {
        1e-7 * value.abs().max(1.0)
    }
}

#[derive(Clone, Copy, Debug)]
pub struct NewtonOptions {
    pub max_iter: usize,
    /// Convergence threshold on the largest (undamped) correction.
    pub tol: f64,
}

impl Default for NewtonOptions {
    fn default() -> Self {
        Self { max_iter: 100, tol: 1e-9 }
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct IterInfo {
    pub iter: usize,
    /// Largest absolute correction proposed by Newton (before damping).
    pub max_correction: f64,
    /// Largest absolute residual before the correction was applied.
    pub residual: f64,
    /// Damping factor applied (1 = full Newton step).
    pub damping: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum StepOutcome {
    Continue(IterInfo),
    Converged(IterInfo),
    Failed(String),
}

pub struct HenyeyNewton<'a, P: BvpProblem + ?Sized> {
    problem: &'a P,
    y: Vec<f64>,
    opts: NewtonOptions,
    iter: usize,
    jac: BandedMatrix,
    res: Vec<f64>,
    nv: usize,
    ni: usize,
    np: usize,
}

impl<'a, P: BvpProblem + ?Sized> HenyeyNewton<'a, P> {
    pub fn new(problem: &'a P, y0: Vec<f64>, opts: NewtonOptions) -> Self {
        let nv = problem.nvar();
        let ni = problem.n_inner();
        let np = problem.npts();
        assert!(np >= 2, "need at least two mesh points");
        assert!(ni <= nv, "too many inner boundary conditions");
        assert_eq!(y0.len(), nv * np, "initial guess has wrong length");
        let n = nv * np;
        let kl = ni + nv - 1;
        let ku = 2 * nv - 1 - ni;
        Self { problem, y: y0, opts, iter: 0, jac: BandedMatrix::new(n, kl, ku), res: vec![0.0; n], nv, ni, np }
    }

    pub fn solution(&self) -> &[f64] {
        &self.y
    }

    pub fn into_solution(self) -> Vec<f64> {
        self.y
    }

    pub fn iterations(&self) -> usize {
        self.iter
    }

    #[inline]
    fn row_interval(&self, k: usize) -> usize {
        self.ni + self.nv * k
    }

    #[inline]
    fn row_outer(&self) -> usize {
        self.ni + self.nv * (self.np - 1)
    }

    /// Evaluate all residuals for state `y`.
    pub fn residual(&self, y: &[f64], f: &mut [f64]) {
        let (nv, ni, np, p) = (self.nv, self.ni, self.np, self.problem);
        p.inner_bc(&y[0..nv], &mut f[0..ni]);
        for k in 0..np - 1 {
            let r = self.row_interval(k);
            p.interval(k, &y[nv * k..nv * (k + 1)], &y[nv * (k + 1)..nv * (k + 2)], &mut f[r..r + nv]);
        }
        let r = self.row_outer();
        p.outer_bc(&y[nv * (np - 1)..], &mut f[r..r + (nv - ni)]);
    }

    fn assemble_jacobian(&mut self) {
        let (nv, ni, np, p) = (self.nv, self.ni, self.np, self.problem);
        self.jac.zero();
        let mut plus = vec![0.0; nv];
        let mut minus = vec![0.0; nv];
        let mut fp = vec![0.0; nv];
        let mut fm = vec![0.0; nv];
        for k in 0..np {
            for j in 0..nv {
                let col = nv * k + j;
                let v = self.y[col];
                let h = p.fd_step(j, v);
                plus.copy_from_slice(&self.y[nv * k..nv * (k + 1)]);
                minus.copy_from_slice(&plus);
                plus[j] = v + h;
                minus[j] = v - h;
                let inv = 1.0 / (plus[j] - minus[j]);

                if k == 0 {
                    p.inner_bc(&plus, &mut fp[..ni]);
                    p.inner_bc(&minus, &mut fm[..ni]);
                    for i in 0..ni {
                        self.jac.set(i, col, (fp[i] - fm[i]) * inv);
                    }
                }
                if k >= 1 {
                    let prev = &self.y[nv * (k - 1)..nv * k];
                    p.interval(k - 1, prev, &plus, &mut fp);
                    p.interval(k - 1, prev, &minus, &mut fm);
                    let r0 = self.row_interval(k - 1);
                    for i in 0..nv {
                        self.jac.set(r0 + i, col, (fp[i] - fm[i]) * inv);
                    }
                }
                if k + 1 < np {
                    let next = &self.y[nv * (k + 1)..nv * (k + 2)];
                    p.interval(k, &plus, next, &mut fp);
                    p.interval(k, &minus, next, &mut fm);
                    let r0 = self.row_interval(k);
                    for i in 0..nv {
                        self.jac.set(r0 + i, col, (fp[i] - fm[i]) * inv);
                    }
                }
                if k == np - 1 {
                    let no = nv - ni;
                    p.outer_bc(&plus, &mut fp[..no]);
                    p.outer_bc(&minus, &mut fm[..no]);
                    let r0 = self.row_outer();
                    for i in 0..no {
                        self.jac.set(r0 + i, col, (fp[i] - fm[i]) * inv);
                    }
                }
            }
        }
    }

    /// Perform one damped Newton iteration.
    pub fn step(&mut self) -> StepOutcome {
        if self.iter >= self.opts.max_iter {
            return StepOutcome::Failed(format!("no convergence after {} Newton iterations", self.opts.max_iter));
        }
        self.iter += 1;

        let mut res = std::mem::take(&mut self.res);
        self.residual(&self.y, &mut res);
        if let Some(i) = res.iter().position(|v| !v.is_finite()) {
            self.res = res;
            return StepOutcome::Failed(format!(
                "non-finite residual at equation {i} (mesh point ~{})",
                i.saturating_sub(self.ni) / self.nv
            ));
        }
        let residual = res.iter().fold(0.0_f64, |a, v| a.max(v.abs()));

        self.assemble_jacobian();
        let mut dy: Vec<f64> = res.iter().map(|v| -v).collect();
        self.res = res;
        if let Err(e) = self.jac.solve_in_place(&mut dy) {
            return StepOutcome::Failed(format!("linear solve failed: {e}"));
        }
        if dy.iter().any(|v| !v.is_finite()) {
            return StepOutcome::Failed("non-finite Newton correction".into());
        }

        // Damping: limit the largest correction of each variable.
        let mut damping: f64 = 1.0;
        for j in 0..self.nv {
            let limit = self.problem.max_correction(j);
            if limit.is_finite() {
                let biggest = (0..self.np).fold(0.0_f64, |a, k| a.max(dy[self.nv * k + j].abs()));
                if biggest > limit {
                    damping = damping.min(limit / biggest);
                }
            }
        }
        let max_correction = dy.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        for (y, d) in self.y.iter_mut().zip(&dy) {
            *y += damping * d;
        }

        let info = IterInfo { iter: self.iter, max_correction, residual, damping };
        if max_correction < self.opts.tol && damping == 1.0 {
            StepOutcome::Converged(info)
        } else {
            StepOutcome::Continue(info)
        }
    }

    /// Iterate until convergence or failure.
    pub fn run(&mut self) -> Result<IterInfo, String> {
        loop {
            match self.step() {
                StepOutcome::Continue(_) => {}
                StepOutcome::Converged(info) => return Ok(info),
                StepOutcome::Failed(msg) => return Err(msg),
            }
        }
    }
}
