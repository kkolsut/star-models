//! Validate the generic Henyey solver on the Lane–Emden equation. It is solved
//! as a boundary-value problem with the eigenvalue λ = ξ₁² as an extra
//! unknown, and the zeros ξ₁ are known accurately.
//!
//! With x = ξ/ξ₁ ∈ [x₀, 1], the unknowns are (θ, φ = dθ/dx, λ):
//!   θ' = φ,   φ' = −2φ/x − λ θⁿ,   λ' = 0
//!   inner: θ(x₀) = 1 − λx₀²/6,  φ(x₀) = −λx₀/3      outer: θ(1) = 0

use stellar_core::solver::{BvpProblem, HenyeyNewton, NewtonOptions};

struct LaneEmdenBvp {
    n: f64,
    x: Vec<f64>,
}

fn pow_signed(v: f64, n: f64) -> f64 {
    v.signum() * v.abs().powf(n)
}

impl BvpProblem for LaneEmdenBvp {
    fn nvar(&self) -> usize {
        3
    }
    fn n_inner(&self) -> usize {
        2
    }
    fn npts(&self) -> usize {
        self.x.len()
    }
    fn inner_bc(&self, y: &[f64], out: &mut [f64]) {
        let x0 = self.x[0];
        out[0] = y[0] - (1.0 - y[2] * x0 * x0 / 6.0);
        out[1] = y[1] + y[2] * x0 / 3.0;
    }
    fn outer_bc(&self, y: &[f64], out: &mut [f64]) {
        out[0] = y[0];
    }
    fn interval(&self, k: usize, a: &[f64], b: &[f64], out: &mut [f64]) {
        let dx = self.x[k + 1] - self.x[k];
        let xm = 0.5 * (self.x[k] + self.x[k + 1]);
        let (th, ph, lam) = (0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]));
        out[0] = b[0] - a[0] - dx * ph;
        out[1] = b[1] - a[1] - dx * (-2.0 * ph / xm - lam * pow_signed(th, self.n));
        out[2] = b[2] - a[2];
    }
    fn max_correction(&self, var: usize) -> f64 {
        if var == 2 { 5.0 } else { 0.3 }
    }
}

/// Initial guess shapes for θ(x).
enum Guess {
    Parabola,
    /// The n = 1 solution sin(πx)/(πx), a closer shape for centrally condensed polytropes.
    Sinc,
}

fn solve_xi1(n: f64, lambda_guess: f64, guess: Guess) -> f64 {
    let npts = 4000;
    let x0 = 1e-4;
    let x: Vec<f64> = (0..npts).map(|i| x0 + (1.0 - x0) * i as f64 / (npts - 1) as f64).collect();
    let pi = std::f64::consts::PI;
    let y0: Vec<f64> = x
        .iter()
        .flat_map(|&xi| match guess {
            Guess::Parabola => [1.0 - xi * xi, -2.0 * xi, lambda_guess],
            Guess::Sinc => {
                let a = pi * xi;
                [a.sin() / a, (a * a.cos() - a.sin()) / (a * xi), lambda_guess]
            }
        })
        .collect();
    let problem = LaneEmdenBvp { n, x };
    let mut nr = HenyeyNewton::new(&problem, y0, NewtonOptions { max_iter: 100, tol: 1e-11 });
    nr.run().expect("Lane–Emden BVP should converge");
    nr.solution()[2].sqrt()
}

#[test]
fn polytrope_n1_first_zero_is_pi() {
    let xi1 = solve_xi1(1.0, 8.0, Guess::Parabola);
    assert!((xi1 - std::f64::consts::PI).abs() < 1e-5, "ξ₁ = {xi1}");
}

#[test]
fn polytrope_n3_first_zero() {
    let xi1 = solve_xi1(3.0, 40.0, Guess::Sinc);
    assert!((xi1 - 6.896_849).abs() < 1e-4, "ξ₁ = {xi1}");
}
