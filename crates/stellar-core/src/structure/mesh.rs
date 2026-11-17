//! Mass mesh in q = m/M.
//!
//! The points are uniform in a stretched coordinate
//!
//! s(q) = a ln q + b q − c ln(1 − q + ε)
//!
//! which gives roughly:
//! * logarithmic spacing in q near the centre (q ≥ 1e-6),
//! * uniform spacing in q through the interior,
//! * logarithmic spacing in 1 − q towards the surface, down to 1 − q ≈ 1e-13.
//!
//! The outer envelope above the convection zone holds very little mass, but
//! P, T and ρ change by orders of magnitude there, so it needs this resolution.

/// Innermost mesh point, q₀ = m₀/M.
pub const Q_INNER: f64 = 1e-6;
const SURFACE_SOFTENING: f64 = 1e-12;
// Roughly 15% / 45% / 40% of the points in the core / interior / envelope.
const W_CORE: f64 = 0.15 / 13.8155;
const W_MID: f64 = 0.45;
const W_ENV: f64 = 0.40 / 27.631;

/// Stretched mesh coordinate s(q); monotonically increasing on (0, 1].
pub fn coordinate(q: f64) -> f64 {
    W_CORE * q.ln() + W_MID * q - W_ENV * (1.0 - q + SURFACE_SOFTENING).ln()
}

/// Build a mesh of `n` points from q = Q_INNER to q = 1 (inclusive).
pub fn make_mesh(n: usize) -> Vec<f64> {
    assert!(n >= 10, "mesh needs at least 10 points");
    let s0 = coordinate(Q_INNER);
    let s1 = coordinate(1.0);
    (0..n)
        .map(|i| {
            if i == 0 {
                return Q_INNER;
            }
            if i == n - 1 {
                return 1.0;
            }
            let s = s0 + (s1 - s0) * i as f64 / (n - 1) as f64;
            let (mut lo, mut hi) = (Q_INNER, 1.0_f64);
            for _ in 0..200 {
                let mid = 0.5 * (lo + hi);
                if mid <= lo || mid >= hi {
                    break;
                }
                if coordinate(mid) < s {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        })
        .collect()
}

/// Interpolate a solution vector (`nv` variables per point) from one mesh to
/// another, linearly in the stretched coordinate.
pub fn remap(q_old: &[f64], y_old: &[f64], nv: usize, q_new: &[f64]) -> Vec<f64> {
    assert_eq!(y_old.len(), q_old.len() * nv);
    let s_old: Vec<f64> = q_old.iter().map(|&q| coordinate(q)).collect();
    let mut out = Vec::with_capacity(q_new.len() * nv);
    for &q in q_new {
        let s = coordinate(q);
        let i = s_old.partition_point(|&v| v < s).clamp(1, s_old.len() - 1);
        let (s0, s1) = (s_old[i - 1], s_old[i]);
        let w = if s1 > s0 { ((s - s0) / (s1 - s0)).clamp(0.0, 1.0) } else { 0.0 };
        for j in 0..nv {
            let a = y_old[(i - 1) * nv + j];
            let b = y_old[i * nv + j];
            out.push(a + w * (b - a));
        }
    }
    out
}

/// Interpolate a solution vector (`nv` variables per point) given on the mesh
/// `q` of a star of mass `mass_old` to the same mesh `q` of a star of mass
/// `mass_new`, at fixed mass coordinate m: point k of the new star, at
/// m = q_k · mass_new, takes the old values at q = q_k · mass_new / mass_old.
///
/// This follows the same gas when the total mass changes (a wind removes the
/// outermost layers). The interpolation is linear in the stretched coordinate,
/// as in [`remap`]; points outside the old mesh take its end values.
pub fn remap_fixed_mass(q: &[f64], y: &[f64], nv: usize, mass_old: f64, mass_new: f64) -> Vec<f64> {
    let ratio = mass_new / mass_old;
    let targets: Vec<f64> = q.iter().map(|&qk| (qk * ratio).clamp(q[0], 1.0)).collect();
    remap(q, y, nv, &targets)
}

/// ∫_{m_cut}^{m_n} f dm with the trapezoid rule on the points `m` (increasing),
/// `f` interpolated linearly at `m_cut`. Zero if `m_cut ≥ m_n`.
pub fn outer_integral(m: &[f64], f: &[f64], m_cut: f64) -> f64 {
    let n = m.len();
    if !(m_cut < m[n - 1]) {
        return 0.0;
    }
    let i = m.partition_point(|&v| v <= m_cut).clamp(1, n - 1);
    let w = ((m_cut - m[i - 1]) / (m[i] - m[i - 1])).clamp(0.0, 1.0);
    let f_cut = f[i - 1] + w * (f[i] - f[i - 1]);
    let m_lo = m_cut.max(m[i - 1]);
    let mut sum = 0.5 * (f_cut + f[i]) * (m[i] - m_lo);
    for k in i + 1..n {
        sum += 0.5 * (f[k] + f[k - 1]) * (m[k] - m[k - 1]);
    }
    sum
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn remap_fixed_mass_follows_the_gas() {
        let q = make_mesh(400);
        // A profile that is a function of m alone, with M_old = 1.
        let f = |m: f64| (1.0 - m + 1e-14).ln() + 3.0 * m;
        let y: Vec<f64> = q.iter().map(|&qk| f(qk)).collect();
        let m_new = 0.999;
        let out = remap_fixed_mass(&q, &y, 1, 1.0, m_new);
        for (k, &qk) in q.iter().enumerate().take(q.len() - 1) {
            let want = f(qk * m_new);
            assert!((out[k] - want).abs() < 2e-3 * want.abs().max(1.0), "k = {k}: {} vs {want}", out[k]);
        }
        // Unchanged mass: identity.
        let same = remap_fixed_mass(&q, &y, 1, 1.0, 1.0);
        assert!(same.iter().zip(&y).all(|(a, b)| (a - b).abs() < 1e-12));
    }

    #[test]
    fn outer_integral_of_a_linear_function() {
        let m: Vec<f64> = (0..=10).map(|i| i as f64).collect();
        let f: Vec<f64> = m.iter().map(|&v| 2.0 * v).collect();
        // ∫_{7.5}^{10} 2m dm = 100 − 56.25
        assert!((outer_integral(&m, &f, 7.5) - 43.75).abs() < 1e-12);
        assert!((outer_integral(&m, &f, 7.0) - 51.0).abs() < 1e-12);
        assert_eq!(outer_integral(&m, &f, 10.0), 0.0);
    }

    #[test]
    fn mesh_is_monotonic_and_resolves_surface() {
        let q = make_mesh(800);
        assert_eq!(q[0], Q_INNER);
        assert_eq!(*q.last().unwrap(), 1.0);
        assert!(q.windows(2).all(|w| w[1] > w[0]), "mesh must increase strictly");
        let last_interior = 1.0 - q[q.len() - 2];
        assert!(last_interior < 1e-11, "outermost interval 1−q = {last_interior}");
        // Interior spacing should be modest.
        let max_dq = q.windows(2).map(|w| w[1] - w[0]).fold(0.0, f64::max);
        assert!(max_dq < 0.01, "max Δq = {max_dq}");
    }

    #[test]
    fn remap_identity() {
        let q = make_mesh(100);
        let y: Vec<f64> = q.iter().flat_map(|&v| [v, 2.0 * v]).collect();
        let back = remap(&q, &y, 2, &q);
        for (a, b) in y.iter().zip(&back) {
            assert!((a - b).abs() < 1e-12);
        }
    }
}
