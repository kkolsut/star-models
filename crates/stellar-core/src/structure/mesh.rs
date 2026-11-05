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

#[cfg(test)]
mod tests {
    use super::*;

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
