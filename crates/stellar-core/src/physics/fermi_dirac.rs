//! Generalised Fermi–Dirac integrals (plan 6.2, `docs/physics/02`)
//!
//! F_k(η, ϑ) = ∫₀^∞ x^k (1 + ϑx/2)^{1/2} / (e^{x−η} + 1) dx,   k = 1/2, 3/2, 5/2,
//!
//! with their derivatives with respect to η and ϑ, all from one quadrature
//! pass.
//!
//! *Quadrature.* The range is split at break points that are smooth functions
//! of η: with x_b = ln(1 + e^η) (≈ η for η ≫ 1, → 0 for η ≪ −1),
//!
//! [0, b₁], [b₁, b₂], [b₂, x_b], [x_b, x_b + 12], [x_b + 12, ∞),
//! b₁ = min(ln(1 + e^{x_b − 25}), x_b), b₂ = ln(1 + e^{x_b − 5}) limited to [b₁, x_b].
//!
//! The first four are integrated in v = √x (which removes the x^{1/2}
//! behaviour at the origin, also when a segment starts just above it) with
//! fixed Gauss–Legendre rules; the tail is integrated with a Gauss–Laguerre
//! rule. Fixed nodes, no adaptive control: the result is a smooth function of
//! (η, ϑ) for the finite-difference Jacobian.
//!
//! *Domain and accuracy* (test `matches_high_precision_reference`, against
//! values computed with mpmath at 40 digits): relative error below 10⁻¹² for
//! −100 ≤ η ≤ 10⁴, 10⁻⁸ ≤ ϑ ≤ 1 and ϑ·x_b ≤ 10 (the electron Fermi energy
//! up to 5 m_e c²). Outside that domain [`fermi_dirac`] returns an error; below
//! η = −100 the classical (non-degenerate) limit F_k ∝ e^η is used.

use std::sync::OnceLock;

/// Lower bound of the quadrature domain in η; below it the classical series.
pub const ETA_MIN: f64 = -100.0;
pub const ETA_MAX: f64 = 1e4;
pub const THETA_MIN: f64 = 1e-8;
pub const THETA_MAX: f64 = 1.0;
/// Largest ϑ·x_b (relativistic degeneracy) inside the domain.
pub const THETA_XB_MAX: f64 = 10.0;

/// F_k, ∂F_k/∂η and ∂F_k/∂ϑ for k = 1/2, 3/2, 5/2 (index 0, 1, 2).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FdIntegrals {
    pub f: [f64; 3],
    pub d_eta: [f64; 3],
    pub d_theta: [f64; 3],
}

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum FdError {
    #[error("Fermi–Dirac integral outside its domain: η = {eta}, ϑ = {theta}")]
    OutOfDomain { eta: f64, theta: f64 },
}

struct Rules {
    leg_a: (Vec<f64>, Vec<f64>),
    leg_b: (Vec<f64>, Vec<f64>),
    lag: (Vec<f64>, Vec<f64>),
}

fn rules() -> &'static Rules {
    static R: OnceLock<Rules> = OnceLock::new();
    R.get_or_init(|| Rules { leg_a: gauss_legendre(32), leg_b: gauss_legendre(40), lag: gauss_laguerre(32) })
}

/// Gauss–Legendre nodes and weights on [0, 1].
pub fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    // Golub–Welsch: Jacobi matrix with b_k = k / √(4k² − 1).
    let diag = vec![0.0; n];
    let off: Vec<f64> = (1..n).map(|k| k as f64 / ((4 * k * k - 1) as f64).sqrt()).collect();
    let (x, w0) = golub_welsch(&diag, &off, 2.0);
    // Map [−1, 1] → [0, 1].
    let mut pairs: Vec<(f64, f64)> = x.iter().zip(&w0).map(|(x, w)| (0.5 * (x + 1.0), 0.5 * w)).collect();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Polish the nodes with Newton on P_n and recompute weights.
    for p in &mut pairs {
        let mut t = 2.0 * p.0 - 1.0;
        for _ in 0..3 {
            let (pn, dpn) = legendre(n, t);
            t -= pn / dpn;
        }
        let (_, dpn) = legendre(n, t);
        p.0 = 0.5 * (t + 1.0);
        p.1 = 1.0 / ((1.0 - t * t) * dpn * dpn); // 2/((1−t²)P'²) on [−1,1], halved
    }
    pairs.into_iter().unzip()
}

fn legendre(n: usize, x: f64) -> (f64, f64) {
    let (mut p0, mut p1) = (1.0, x);
    for k in 2..=n {
        let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
        p0 = p1;
        p1 = p2;
    }
    let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
    (p1, dp)
}

/// Gauss–Laguerre nodes and weights for ∫₀^∞ e^{−x} f(x) dx.
pub fn gauss_laguerre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let diag: Vec<f64> = (0..n).map(|k| (2 * k + 1) as f64).collect();
    let off: Vec<f64> = (1..n).map(|k| k as f64).collect();
    let (x, w) = golub_welsch(&diag, &off, 1.0);
    let mut pairs: Vec<(f64, f64)> = x.into_iter().zip(w).collect();
    pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
    // Polish nodes with Newton on L_n; weights x / ((n+1) L_{n+1}(x))².
    for p in &mut pairs {
        let mut t = p.0;
        for _ in 0..3 {
            let (l, dl, _) = laguerre(n, t);
            t -= l / dl;
        }
        let (_, _, ln1) = laguerre(n, t);
        p.0 = t;
        p.1 = t / (((n + 1) as f64) * ln1).powi(2);
    }
    pairs.into_iter().unzip()
}

/// L_n(x), L_n'(x), L_{n+1}(x).
fn laguerre(n: usize, x: f64) -> (f64, f64, f64) {
    let (mut l0, mut l1) = (1.0, 1.0 - x);
    for k in 1..n {
        let l2 = ((2 * k + 1) as f64 - x) * l1 / (k + 1) as f64 - k as f64 * l0 / (k + 1) as f64;
        l0 = l1;
        l1 = l2;
    }
    // l1 = L_n, l0 = L_{n−1}
    let ln1 = ((2 * n + 1) as f64 - x) * l1 / (n + 1) as f64 - n as f64 * l0 / (n + 1) as f64;
    let dl = n as f64 * (l1 - l0) / x;
    (l1, dl, ln1)
}

/// Eigenvalues and first components of the eigenvectors of a symmetric
/// tridiagonal matrix (implicit QL), giving Gauss nodes and weights μ₀ v₀².
fn golub_welsch(diag: &[f64], off: &[f64], mu0: f64) -> (Vec<f64>, Vec<f64>) {
    let n = diag.len();
    let mut d = diag.to_vec();
    let mut e = off.to_vec();
    e.push(0.0);
    let mut z = vec![0.0; n];
    z[0] = 1.0;
    for l in 0..n {
        let mut iter = 0;
        loop {
            let mut m = l;
            while m < n - 1 {
                let dd = d[m].abs() + d[m + 1].abs();
                if e[m].abs() <= f64::EPSILON * dd {
                    break;
                }
                m += 1;
            }
            if m == l {
                break;
            }
            iter += 1;
            assert!(iter < 60, "QL did not converge");
            let mut g = (d[l + 1] - d[l]) / (2.0 * e[l]);
            let mut r = g.hypot(1.0);
            g = d[m] - d[l] + e[l] / (g + r.copysign(g));
            let (mut s, mut c, mut p) = (1.0, 1.0, 0.0);
            let mut i = m;
            while i > l {
                i -= 1;
                let mut f = s * e[i];
                let b = c * e[i];
                r = f.hypot(g);
                e[i + 1] = r;
                if r == 0.0 {
                    d[i + 1] -= p;
                    e[m] = 0.0;
                    break;
                }
                s = f / r;
                c = g / r;
                g = d[i + 1] - p;
                r = (d[i] - g) * s + 2.0 * c * b;
                p = s * r;
                d[i + 1] = g + p;
                g = c * r - b;
                f = z[i + 1];
                z[i + 1] = s * z[i] + c * f;
                z[i] = c * z[i] - s * f;
            }
            if r == 0.0 && i >= l {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
    }
    let w = z.iter().map(|v| mu0 * v * v).collect();
    (d, w)
}

/// Smooth positive part ln(1 + e^x), without overflow.
#[inline]
pub fn softplus(x: f64) -> f64 {
    if x > 30.0 { x + (-x).exp().ln_1p() } else { x.exp().ln_1p() }
}

/// Generalised Fermi–Dirac integrals at (η, ϑ).
pub fn fermi_dirac(eta: f64, theta: f64) -> Result<FdIntegrals, FdError> {
    let xb = softplus(eta);
    if !(eta <= ETA_MAX && (THETA_MIN..=THETA_MAX).contains(&theta) && theta * xb <= THETA_XB_MAX) {
        return Err(FdError::OutOfDomain { eta, theta });
    }
    if eta < ETA_MIN {
        return Ok(classical(eta, theta));
    }
    let r = rules();
    let mut acc = FdIntegrals::default();
    let b1 = softplus(xb - 25.0).min(xb);
    let b2 = softplus(xb - 5.0).clamp(b1, xb);
    let segs = [(0.0, b1, &r.leg_a), (b1, b2, &r.leg_a), (b2, xb, &r.leg_a), (xb, xb + 12.0, &r.leg_b)];
    for (a, b, rule) in segs {
        // Every finite segment in v = √x: the integrand 2v·x^k·√(1 + ϑx/2)·f
        // is smooth in v, wherever the segment starts.
        let (va, vb) = (a.sqrt(), b.sqrt());
        let dv = vb - va;
        if !(dv > 0.0) {
            continue;
        }
        for (u, w) in rule.0.iter().zip(&rule.1) {
            let v = va + dv * u;
            accumulate(&mut acc, v * v, 2.0 * v * dv * w, eta, theta, None);
        }
    }
    let x0 = xb + 12.0;
    for (t, w) in r.lag.0.iter().zip(&r.lag.1) {
        let x = x0 + t;
        // e^{−t} is in the weight: the Fermi factor is e^{η−x} / (1 + e^{η−x})
        // = e^{−t} · e^{η−x₀} / (1 + e^{η−x}).
        accumulate(&mut acc, x, *w, eta, theta, Some(eta - x0));
    }
    Ok(acc)
}

/// Add one node. With `tail = Some(η − x₀)`, the weight already contains e^{−t}.
#[inline]
fn accumulate(acc: &mut FdIntegrals, x: f64, w: f64, eta: f64, theta: f64, tail: Option<f64>) {
    let e = (x - eta).min(700.0);
    // f = 1/(e^{x−η} + 1); f(1 − f) = ∂f/∂η.
    let (f, df) = match tail {
        None => {
            let f = 1.0 / (e.exp() + 1.0);
            (f, f * (1.0 - f))
        }
        Some(shift) => {
            let q = (-e).exp(); // e^{η−x}
            let f = shift.exp() / (1.0 + q);
            (f, f / (1.0 + q))
        }
    };
    let g = (1.0 + 0.5 * theta * x).sqrt();
    let sx = x.sqrt();
    let mut xk = sx; // x^{1/2}
    for k in 0..3 {
        let base = w * xk * g;
        acc.f[k] += base * f;
        acc.d_eta[k] += base * df;
        acc.d_theta[k] += w * xk * x * f / (4.0 * g);
        xk *= x;
    }
}

/// Classical (non-degenerate) limit below η = ETA_MIN: F_k(η) = e^{η − η_min} F_k(η_min),
/// exact up to terms of relative size e^{η_min} ≈ 10⁻⁴⁴ (the series
/// Σ (−1)^{j+1} e^{jη} I_k(j) to its first term).
fn classical(eta: f64, theta: f64) -> FdIntegrals {
    let mut a = fermi_dirac(ETA_MIN, theta).expect("inside the domain");
    let s = (eta - ETA_MIN).exp();
    for k in 0..3 {
        a.f[k] *= s;
        a.d_eta[k] = a.f[k];
        a.d_theta[k] *= s;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss_rules_integrate_polynomials() {
        let (x, w) = gauss_legendre(32);
        let s: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(20)).sum();
        assert!((s - 1.0 / 21.0).abs() < 1e-15);
        let (x, w) = gauss_laguerre(32);
        let s: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(10)).sum();
        assert!((s / 3_628_800.0 - 1.0).abs() < 1e-13, "{s}");
    }

    /// High-precision reference values (mpmath, 40 digits; script in
    /// `tools/reference/fermi_dirac_reference.py`): (η, ϑ, F_{1/2}, F_{3/2}, F_{5/2}).
    #[allow(clippy::excessive_precision)]
    const REFERENCE: &[(f64, f64, [f64; 3])] = &include!("fermi_dirac_reference.in");

    #[test]
    fn matches_high_precision_reference() {
        let mut worst = 0.0_f64;
        for &(eta, theta, want) in REFERENCE {
            let got = fermi_dirac(eta, theta).unwrap();
            for k in 0..3 {
                let rel = (got.f[k] / want[k] - 1.0).abs();
                worst = worst.max(rel);
                assert!(rel < 1e-12, "η = {eta}, ϑ = {theta}, k = {k}: {} vs {} ({rel:.1e})", got.f[k], want[k]);
            }
        }
        eprintln!("largest relative error {worst:.1e} over {} points", REFERENCE.len());
    }

    #[test]
    fn derivatives_match_differences() {
        for &(eta, theta) in &[(-30.0, 1e-3), (-2.0, 0.01), (0.5, 0.1), (12.0, 0.02), (300.0, 0.01), (5000.0, 1e-3)] {
            let a = fermi_dirac(eta, theta).unwrap();
            let h = 1e-5 * eta.abs().max(1.0);
            let p = fermi_dirac(eta + h, theta).unwrap();
            let m = fermi_dirac(eta - h, theta).unwrap();
            let ht = 1e-5 * theta;
            let pt = fermi_dirac(eta, theta + ht).unwrap();
            let mt = fermi_dirac(eta, theta - ht).unwrap();
            for k in 0..3 {
                let fd = (p.f[k] - m.f[k]) / (2.0 * h);
                assert!((fd / a.d_eta[k] - 1.0).abs() < 1e-6, "η = {eta}, k = {k}: {fd} vs {}", a.d_eta[k]);
                let fdt = (pt.f[k] - mt.f[k]) / (2.0 * ht);
                assert!((fdt / a.d_theta[k] - 1.0).abs() < 1e-6, "ϑ, η = {eta}, k = {k}: {fdt} vs {}", a.d_theta[k]);
            }
        }
    }

    #[test]
    fn non_degenerate_and_degenerate_limits() {
        // η → −∞, ϑ → 0: F_k → Γ(k + 1) e^η.
        let g = [0.886_226_925_452_758, 1.329_340_388_179_137, 3.323_350_970_447_843];
        let a = fermi_dirac(-40.0, 1e-6).unwrap();
        for k in 0..3 {
            assert!((a.f[k] / ((-40.0f64).exp() * g[k]) - 1.0).abs() < 1e-6);
        }
        // η ≫ 1, ϑ → 0: F_{1/2} → (2/3) η^{3/2} (1 + π²/(8η²)) (Sommerfeld), and
        // the relativistic factor √(1 + ϑx/2) adds (3/20)ϑη to first order.
        let (eta, th) = (2000.0_f64, 1e-8);
        let a = fermi_dirac(eta, th).unwrap();
        let want = 2.0 / 3.0
            * eta.powf(1.5)
            * (1.0 + std::f64::consts::PI.powi(2) / (8.0 * eta * eta))
            * (1.0 + 0.15 * th * eta);
        assert!((a.f[0] / want - 1.0).abs() < 1e-9, "{} vs {want}", a.f[0]);
        // Continuity at the switch to the classical series.
        let e = 1e-9;
        let (l, r) = (fermi_dirac(ETA_MIN - e, 0.01).unwrap(), fermi_dirac(ETA_MIN + e, 0.01).unwrap());
        // F ∝ e^η here: compare F e^{−η}.
        assert!((l.f[0] / r.f[0] * (2.0 * e).exp() - 1.0).abs() < 1e-12);
        assert!(fermi_dirac(1.0, 2.0).is_err());
        assert!(fermi_dirac(1e5, 1e-3).is_err());
    }
}
