//! Banded linear solver: Gaussian elimination with partial (row) pivoting.
//!
//! Storage is row-major. Each row `i` keeps a window of columns
//! `[i − kl, i + kl + ku]`, which leaves room for the fill-in that row
//! pivoting creates (the same idea as LAPACK's `gbtrf`).

// Explicit index loops mirror the band arithmetic and are clearer here.
#![allow(clippy::needless_range_loop)]

#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[error("singular matrix (zero pivot in column {column})")]
pub struct SingularMatrix {
    pub column: usize,
}

#[derive(Debug, Clone)]
pub struct BandedMatrix {
    n: usize,
    kl: usize,
    ku: usize,
    width: usize,
    data: Vec<f64>,
}

impl BandedMatrix {
    /// A zero n×n matrix with `kl` sub-diagonals and `ku` super-diagonals.
    pub fn new(n: usize, kl: usize, ku: usize) -> Self {
        let width = 2 * kl + ku + 1;
        Self { n, kl, ku, width, data: vec![0.0; n * width] }
    }

    pub fn n(&self) -> usize {
        self.n
    }

    pub fn zero(&mut self) {
        self.data.fill(0.0);
    }

    #[inline]
    fn in_window(&self, i: usize, j: usize) -> bool {
        j + self.kl >= i && j <= i + self.kl + self.ku
    }

    #[inline]
    fn idx(&self, i: usize, j: usize) -> usize {
        debug_assert!(self.in_window(i, j), "({i}, {j}) outside band");
        i * self.width + (j + self.kl - i)
    }

    pub fn get(&self, i: usize, j: usize) -> f64 {
        if self.in_window(i, j) { self.data[self.idx(i, j)] } else { 0.0 }
    }

    /// Set element (i, j). It must lie inside the original band
    /// (−kl ≤ j − i ≤ ku).
    pub fn set(&mut self, i: usize, j: usize, v: f64) {
        assert!(j + self.kl >= i && j <= i + self.ku, "element ({i}, {j}) outside band kl={} ku={}", self.kl, self.ku);
        let k = self.idx(i, j);
        self.data[k] = v;
    }

    /// Solve A x = b in place, where `b` becomes x. This destroys the matrix.
    pub fn solve_in_place(&mut self, b: &mut [f64]) -> Result<(), SingularMatrix> {
        let n = self.n;
        assert_eq!(b.len(), n);
        let span = self.kl + self.ku;
        for k in 0..n {
            let rmax = (k + self.kl).min(n - 1);
            let cmax = (k + span).min(n - 1);
            // Partial pivoting.
            let mut p = k;
            let mut best = self.data[self.idx(k, k)].abs();
            for r in k + 1..=rmax {
                let v = self.data[self.idx(r, k)].abs();
                if v > best {
                    best = v;
                    p = r;
                }
            }
            if !(best > 0.0) || !best.is_finite() {
                return Err(SingularMatrix { column: k });
            }
            if p != k {
                for c in k..=cmax {
                    let (ik, ip) = (self.idx(k, c), self.idx(p, c));
                    self.data.swap(ik, ip);
                }
                b.swap(k, p);
            }
            let pivot = self.data[self.idx(k, k)];
            for r in k + 1..=rmax {
                let irk = self.idx(r, k);
                let f = self.data[irk] / pivot;
                if f == 0.0 {
                    continue;
                }
                self.data[irk] = 0.0;
                for c in k + 1..=cmax {
                    let ikc = self.idx(k, c);
                    let irc = self.idx(r, c);
                    self.data[irc] -= f * self.data[ikc];
                }
                b[r] -= f * b[k];
            }
        }
        // Back substitution.
        for k in (0..n).rev() {
            let cmax = (k + span).min(n - 1);
            let mut s = b[k];
            for c in k + 1..=cmax {
                s -= self.data[self.idx(k, c)] * b[c];
            }
            b[k] = s / self.data[self.idx(k, k)];
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Small deterministic PRNG so the tests need no extra dependency.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64 / (1u64 << 53) as f64) * 2.0 - 1.0
        }
    }

    fn dense_solve(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Vec<f64> {
        let n = b.len();
        for k in 0..n {
            let p = (k..n).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs())).unwrap();
            a.swap(k, p);
            b.swap(k, p);
            for r in k + 1..n {
                let f = a[r][k] / a[k][k];
                for c in k..n {
                    a[r][c] -= f * a[k][c];
                }
                b[r] -= f * b[k];
            }
        }
        let mut x = vec![0.0; n];
        for k in (0..n).rev() {
            let s: f64 = (k + 1..n).map(|c| a[k][c] * x[c]).sum();
            x[k] = (b[k] - s) / a[k][k];
        }
        x
    }

    #[test]
    fn matches_dense_solver() {
        let mut rng = Lcg(42);
        for &(n, kl, ku) in &[(10, 1, 1), (40, 5, 5), (37, 3, 7), (64, 5, 5)] {
            let mut band = BandedMatrix::new(n, kl, ku);
            let mut dense = vec![vec![0.0; n]; n];
            for i in 0..n {
                for j in i.saturating_sub(kl)..=(i + ku).min(n - 1) {
                    // A weak diagonal forces the solver to pivot.
                    let v = rng.next() + if i == j { 0.05 } else { 0.0 };
                    band.set(i, j, v);
                    dense[i][j] = v;
                }
            }
            let b: Vec<f64> = (0..n).map(|_| rng.next()).collect();
            let expected = dense_solve(dense, b.clone());
            let mut x = b;
            band.solve_in_place(&mut x).unwrap();
            for (a, e) in x.iter().zip(&expected) {
                assert!((a - e).abs() < 1e-8 * (1.0 + e.abs()), "{a} vs {e}");
            }
        }
    }

    #[test]
    fn detects_singular() {
        let mut m = BandedMatrix::new(3, 1, 1);
        m.set(0, 0, 1.0);
        m.set(1, 1, 0.0);
        m.set(2, 2, 1.0);
        let mut b = vec![1.0, 1.0, 1.0];
        assert!(m.solve_in_place(&mut b).is_err());
    }
}
