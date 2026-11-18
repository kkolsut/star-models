//! C¹ interpolation on rectilinear grids: piecewise cubic Hermite in one and
//! two dimensions, with node derivatives from the parabola through each node
//! and its neighbours (non-uniform grids allowed).
//!
//! The interpolants reproduce the node values exactly and have continuous
//! first derivatives everywhere (the derivative at a node is shared by the two
//! cells that meet there), as the finite-difference Jacobian of the solver
//! requires.

/// Derivatives at the nodes of a 1-D table (`x` strictly increasing): the
/// derivative of the parabola through each node and its two neighbours
/// (one-sided parabola at the ends, linear for two nodes).
pub fn node_slopes(x: &[f64], f: &[f64]) -> Vec<f64> {
    let n = x.len();
    assert_eq!(n, f.len());
    if n == 1 {
        return vec![0.0];
    }
    if n == 2 {
        let s = (f[1] - f[0]) / (x[1] - x[0]);
        return vec![s, s];
    }
    let mut d = vec![0.0; n];
    for i in 1..n - 1 {
        let (h0, h1) = (x[i] - x[i - 1], x[i + 1] - x[i]);
        let (s0, s1) = ((f[i] - f[i - 1]) / h0, (f[i + 1] - f[i]) / h1);
        d[i] = (h1 * s0 + h0 * s1) / (h0 + h1);
    }
    let (h0, h1) = (x[1] - x[0], x[2] - x[1]);
    let (s0, s1) = ((f[1] - f[0]) / h0, (f[2] - f[1]) / h1);
    d[0] = ((2.0 * h0 + h1) * s0 - h0 * s1) / (h0 + h1);
    let m = n - 1;
    let (h0, h1) = (x[m - 1] - x[m - 2], x[m] - x[m - 1]);
    let (s0, s1) = ((f[m - 1] - f[m - 2]) / h0, (f[m] - f[m - 1]) / h1);
    d[m] = ((2.0 * h1 + h0) * s1 - h1 * s0) / (h0 + h1);
    d
}

/// Index i of the cell [x_i, x_{i+1}] containing `v` (clamped to the table).
pub fn cell(x: &[f64], v: f64) -> usize {
    let n = x.len();
    debug_assert!(n >= 2);
    x.partition_point(|&a| a <= v).clamp(1, n - 1) - 1
}

/// Cubic Hermite basis on [0, 1]: (h00, h10, h01, h11) and their derivatives.
#[inline]
fn hermite(t: f64) -> ([f64; 4], [f64; 4]) {
    let t2 = t * t;
    let t3 = t2 * t;
    (
        [2.0 * t3 - 3.0 * t2 + 1.0, t3 - 2.0 * t2 + t, -2.0 * t3 + 3.0 * t2, t3 - t2],
        [6.0 * t2 - 6.0 * t, 3.0 * t2 - 4.0 * t + 1.0, -6.0 * t2 + 6.0 * t, 3.0 * t2 - 2.0 * t],
    )
}

/// 1-D C¹ cubic Hermite interpolant.
#[derive(Clone, Debug)]
pub struct Hermite1 {
    pub x: Vec<f64>,
    pub f: Vec<f64>,
    pub d: Vec<f64>,
}

impl Hermite1 {
    pub fn new(x: Vec<f64>, f: Vec<f64>) -> Self {
        let d = node_slopes(&x, &f);
        Self { x, f, d }
    }

    /// Value and derivative at `v`, which must lie inside the table (it is
    /// clamped to the end cells and the cubic is evaluated there).
    pub fn eval(&self, v: f64) -> (f64, f64) {
        if self.x.len() == 1 {
            return (self.f[0], 0.0);
        }
        let i = cell(&self.x, v);
        let h = self.x[i + 1] - self.x[i];
        let t = (v - self.x[i]) / h;
        let (b, db) = hermite(t);
        let val = b[0] * self.f[i] + b[1] * h * self.d[i] + b[2] * self.f[i + 1] + b[3] * h * self.d[i + 1];
        let der = (db[0] * self.f[i] + db[2] * self.f[i + 1]) / h + db[1] * self.d[i] + db[3] * self.d[i + 1];
        (val, der)
    }
}

/// Weights of the C¹ cubic Hermite interpolant at `v` with respect to the
/// node values: f(v) = Σ w_j f_j over the returned (index, weight) pairs.
/// The node slopes are those of [`node_slopes`], which are linear in f, so
/// the interpolant is a fixed linear combination of at most four nodes.
pub fn hermite_weights(x: &[f64], v: f64) -> Vec<(usize, f64)> {
    let n = x.len();
    if n == 1 {
        return vec![(0, 1.0)];
    }
    let i = cell(x, v);
    let h = x[i + 1] - x[i];
    let t = ((v - x[i]) / h).clamp(0.0, 1.0);
    let (b, _) = hermite(t);
    let mut w: Vec<(usize, f64)> = vec![(i, b[0]), (i + 1, b[2])];
    // Slopes as linear combinations of node values.
    let slope = |k: usize| -> Vec<(usize, f64)> {
        if n == 2 {
            let s = 1.0 / (x[1] - x[0]);
            return vec![(0, -s), (1, s)];
        }
        if k == 0 {
            let (h0, h1) = (x[1] - x[0], x[2] - x[1]);
            let a = (2.0 * h0 + h1) / (h0 + h1);
            let c = h0 / (h0 + h1);
            // a s0 − c s1
            return vec![(0, -a / h0), (1, a / h0 + c / h1), (2, -c / h1)];
        }
        if k == n - 1 {
            let m = n - 1;
            let (h0, h1) = (x[m - 1] - x[m - 2], x[m] - x[m - 1]);
            let a = (2.0 * h1 + h0) / (h0 + h1);
            let c = h1 / (h0 + h1);
            // a s1 − c s0
            return vec![(m - 2, c / h0), (m - 1, -a / h1 - c / h0), (m, a / h1)];
        }
        let (h0, h1) = (x[k] - x[k - 1], x[k + 1] - x[k]);
        let (a, c) = (h1 / (h0 + h1), h0 / (h0 + h1));
        // a s0 + c s1
        vec![(k - 1, -a / h0), (k, a / h0 - c / h1), (k + 1, c / h1)]
    };
    for (k, coef) in [(i, b[1] * h), (i + 1, b[3] * h)] {
        for (j, s) in slope(k) {
            w.push((j, coef * s));
        }
    }
    // Merge duplicates.
    w.sort_by_key(|p| p.0);
    let mut out: Vec<(usize, f64)> = Vec::with_capacity(4);
    for (j, v) in w {
        match out.last_mut() {
            Some(l) if l.0 == j => l.1 += v,
            _ => out.push((j, v)),
        }
    }
    out
}

/// 2-D C¹ bicubic Hermite interpolant on a rectilinear grid, values stored row
/// by row (`f[i * ny + j]` at (x_i, y_j)).
#[derive(Clone, Debug)]
pub struct Hermite2 {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    f: Vec<f64>,
    fx: Vec<f64>,
    fy: Vec<f64>,
    fxy: Vec<f64>,
}

impl Hermite2 {
    pub fn new(x: Vec<f64>, y: Vec<f64>, f: Vec<f64>) -> Self {
        let (nx, ny) = (x.len(), y.len());
        assert_eq!(f.len(), nx * ny);
        assert!(nx >= 2 && ny >= 2);
        let mut fx = vec![0.0; nx * ny];
        let mut fy = vec![0.0; nx * ny];
        let mut fxy = vec![0.0; nx * ny];
        for i in 0..nx {
            let row: Vec<f64> = (0..ny).map(|j| f[i * ny + j]).collect();
            let d = node_slopes(&y, &row);
            for j in 0..ny {
                fy[i * ny + j] = d[j];
            }
        }
        for j in 0..ny {
            let col: Vec<f64> = (0..nx).map(|i| f[i * ny + j]).collect();
            let d = node_slopes(&x, &col);
            let colfy: Vec<f64> = (0..nx).map(|i| fy[i * ny + j]).collect();
            let dd = node_slopes(&x, &colfy);
            for i in 0..nx {
                fx[i * ny + j] = d[i];
                fxy[i * ny + j] = dd[i];
            }
        }
        Self { x, y, f, fx, fy, fxy }
    }

    /// Value at node (i, j).
    pub fn node(&self, i: usize, j: usize) -> f64 {
        self.f[i * self.y.len() + j]
    }

    /// Value and the two partial derivatives at (u, v). Outside the grid the
    /// interpolant is continued C¹ by its tangent plane at the nearest
    /// boundary point (linear extrapolation); callers decide whether to use it.
    pub fn eval(&self, u: f64, v: f64) -> (f64, f64, f64) {
        let (x0, x1) = (self.x[0], *self.x.last().unwrap());
        let (y0, y1) = (self.y[0], *self.y.last().unwrap());
        let (uc, vc) = (u.clamp(x0, x1), v.clamp(y0, y1));
        let (f, fu, fv) = self.eval_inside(uc, vc);
        (f + fu * (u - uc) + fv * (v - vc), fu, fv)
    }

    fn eval_inside(&self, u: f64, v: f64) -> (f64, f64, f64) {
        let ny = self.y.len();
        let i = cell(&self.x, u);
        let j = cell(&self.y, v);
        let hx = self.x[i + 1] - self.x[i];
        let hy = self.y[j + 1] - self.y[j];
        let (bx, dbx) = hermite((u - self.x[i]) / hx);
        let (by, dby) = hermite((v - self.y[j]) / hy);
        let mut val = 0.0;
        let mut du = 0.0;
        let mut dv = 0.0;
        for (a, ii) in [(0usize, i), (1, i + 1)] {
            for (b, jj) in [(0usize, j), (1, j + 1)] {
                let k = ii * ny + jj;
                // Basis indices: value at corner a → bx[2a], slope → bx[2a+1].
                let (ax, sx) = (2 * a, 2 * a + 1);
                let (ay, sy) = (2 * b, 2 * b + 1);
                let terms = [
                    (self.f[k], ax, ay, 1.0),
                    (self.fx[k], sx, ay, hx),
                    (self.fy[k], ax, sy, hy),
                    (self.fxy[k], sx, sy, hx * hy),
                ];
                for (c, px, py, scale) in terms {
                    let c = c * scale;
                    val += c * bx[px] * by[py];
                    du += c * dbx[px] * by[py] / hx;
                    dv += c * bx[px] * dby[py] / hy;
                }
            }
        }
        (val, du, dv)
    }
}

/// C¹ step from 0 (at t ≤ 0) to 1 (at t ≥ 1): 3t² − 2t³.
pub fn smoothstep(t: f64) -> f64 {
    let t = t.clamp(0.0, 1.0);
    t * t * (3.0 - 2.0 * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hermite1_reproduces_nodes_and_cubics_on_uniform_grids() {
        let x: Vec<f64> = (0..8).map(|i| i as f64 * 0.5).collect();
        let f: Vec<f64> = x.iter().map(|v| v * v).collect();
        let h = Hermite1::new(x.clone(), f.clone());
        for (xi, fi) in x.iter().zip(&f) {
            assert!((h.eval(*xi).0 - fi).abs() < 1e-14);
        }
        // Parabola slopes are exact for quadratics.
        for k in 0..70 {
            let v = k as f64 * 0.05;
            let (val, der) = h.eval(v);
            assert!((val - v * v).abs() < 1e-12, "{v}");
            assert!((der - 2.0 * v).abs() < 1e-12);
        }
    }

    #[test]
    fn hermite_weights_match_the_interpolant() {
        let x = vec![0.0, 0.3, 1.0, 1.2, 2.5, 4.0];
        let f = vec![1.0, -0.5, 2.0, 0.3, 4.0, -1.0];
        let h = Hermite1::new(x.clone(), f.clone());
        for k in 0..=80 {
            let v = k as f64 * 0.05;
            let w = hermite_weights(&x, v);
            let s: f64 = w.iter().map(|(j, c)| c * f[*j]).sum();
            assert!((s - h.eval(v).0).abs() < 1e-12, "v = {v}: {s} vs {}", h.eval(v).0);
        }
    }

    #[test]
    fn hermite2_is_c1_and_exact_at_nodes() {
        let x = vec![0.0, 0.1, 0.25, 0.5, 0.8, 1.0];
        let y = vec![-1.0, -0.5, 0.0, 0.5, 1.5];
        let g = |a: f64, b: f64| (2.0 * a).sin() + a * b * b - 0.3 * b;
        let f: Vec<f64> = x.iter().flat_map(|&a| y.iter().map(move |&b| g(a, b))).collect();
        let h = Hermite2::new(x.clone(), y.clone(), f);
        for (i, &a) in x.iter().enumerate() {
            for (j, &b) in y.iter().enumerate() {
                assert!((h.eval(a, b).0 - g(a, b)).abs() < 1e-14);
                assert_eq!(h.node(i, j), g(a, b));
            }
        }
        // Derivatives continuous across a cell edge, and equal to finite differences.
        for &b in &[-0.7, 0.2, 1.1] {
            let e = 1e-9;
            let l = h.eval(0.25 - e, b);
            let r = h.eval(0.25 + e, b);
            assert!((l.1 - r.1).abs() < 1e-6 && (l.2 - r.2).abs() < 1e-6);
            let (_, fu, fv) = h.eval(0.4, b);
            let d = 1e-6;
            let fdu = (h.eval(0.4 + d, b).0 - h.eval(0.4 - d, b).0) / (2.0 * d);
            let fdv = (h.eval(0.4, b + d).0 - h.eval(0.4, b - d).0) / (2.0 * d);
            assert!((fu - fdu).abs() < 1e-7 && (fv - fdv).abs() < 1e-7);
        }
        for &a in &[0.05, 0.6] {
            let e = 1e-9;
            let l = h.eval(a, 0.0 - e);
            let r = h.eval(a, 0.0 + e);
            assert!((l.1 - r.1).abs() < 1e-6 && (l.2 - r.2).abs() < 1e-6);
        }
    }
}
