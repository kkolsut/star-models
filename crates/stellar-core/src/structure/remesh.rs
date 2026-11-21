//! Remeshing during evolution (plan 5.R).
//!
//! The base mesh is fixed in q = m/M ([`super::mesh`]). A hydrogen-burning
//! shell, or the composition step left behind by a receding convective core,
//! can be far thinner than the interior spacing. Before every evolution step
//! the mesh is rebuilt from the base mesh: an interval is bisected (in q) while
//! the previous model changes across it by more than
//!
//! * `dlnp` in ln P or `dlnt` in ln T,
//! * `dx` in X,
//! * `dl` × L_nuc in the nuclear luminosity ∫ε_nuc dm (every channel),
//!
//! the previous model being interpolated linearly in m. Intervals are split
//! worst first, at most `max_depth` times, until the point budget is used. As
//! every step starts again from the base mesh, refinement follows a moving
//! shell and disappears where it is no longer needed; the split points are
//! dyadic subdivisions of the base mesh, so a point that is still needed keeps
//! its place from step to step.
//!
//! *Composition* is remapped conservatively ([`remap_composition`]): every point
//! stands for the mass between the midpoints to its neighbours (the point
//! masses of the composition step), with uniform composition, and a new cell
//! takes the mass-weighted mean of the old cells it overlaps. The mass of every
//! species is conserved to rounding. Where a new cell straddles two old ones
//! this mixes them (numerical diffusion over one cell), which is the price of
//! conservation.
//!
//! *Structure.* The previous model's ln P, ln T (for ε_grav) and the starting
//! solution are interpolated to the new points at fixed m, linearly in the
//! stretched mesh coordinate ([`super::mesh::remap`]).

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use super::model::StarModel;
use crate::composition::Composition;

/// Refinement criteria and limits.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RemeshOptions {
    /// Largest change of ln P across an interval.
    pub dlnp: f64,
    /// Largest change of ln T across an interval.
    pub dlnt: f64,
    /// Largest change of X across an interval.
    pub dx: f64,
    /// Largest share of the nuclear luminosity generated in one interval.
    pub dl: f64,
    /// Largest number of bisections of a base interval.
    pub max_depth: u32,
    /// Largest number of points inserted, as a multiple of the base points.
    pub max_extra: f64,
}

impl Default for RemeshOptions {
    fn default() -> Self {
        Self { dlnp: 0.1, dlnt: 0.05, dx: 0.05, dl: 0.02, max_depth: 8, max_extra: 1.0 }
    }
}

/// The quantities the criteria look at, as functions of q on the old mesh.
struct Fields {
    q: Vec<f64>,
    ln_p: Vec<f64>,
    ln_t: Vec<f64>,
    x: Vec<f64>,
    /// Cumulative nuclear luminosity ∫₀^m ε_nuc dm / L_nuc (0 to 1).
    l_nuc: Vec<f64>,
}

impl Fields {
    fn from_model(model: &StarModel) -> Self {
        let p = &model.profile;
        let n = p.len();
        let mut l = vec![0.0; n];
        l[0] = p.eps[0].max(0.0) * p.m[0];
        for k in 1..n {
            l[k] = l[k - 1] + 0.5 * (p.eps[k].max(0.0) + p.eps[k - 1].max(0.0)) * (p.m[k] - p.m[k - 1]);
        }
        let total = l[n - 1];
        let l_nuc = if total > 0.0 { l.iter().map(|v| v / total).collect() } else { vec![0.0; n] };
        Self {
            q: p.q.clone(),
            ln_p: p.p.iter().map(|v| v.ln()).collect(),
            ln_t: p.t.iter().map(|v| v.ln()).collect(),
            x: p.x.clone(),
            l_nuc,
        }
    }

    /// (ln P, ln T, X, L_nuc) at q, linear in q between the old points.
    fn at(&self, q: f64) -> [f64; 4] {
        let n = self.q.len();
        let i = self.q.partition_point(|&v| v < q).clamp(1, n - 1);
        let w = ((q - self.q[i - 1]) / (self.q[i] - self.q[i - 1])).clamp(0.0, 1.0);
        let l = |v: &[f64]| v[i - 1] + w * (v[i] - v[i - 1]);
        [l(&self.ln_p), l(&self.ln_t), l(&self.x), l(&self.l_nuc)]
    }
}

/// Interval of the candidate mesh, ordered by how badly it violates the criteria.
struct Interval {
    lo: f64,
    hi: f64,
    v_lo: [f64; 4],
    v_hi: [f64; 4],
    depth: u32,
    excess: f64,
}

impl PartialEq for Interval {
    fn eq(&self, o: &Self) -> bool {
        self.excess == o.excess
    }
}
impl Eq for Interval {}
impl PartialOrd for Interval {
    fn partial_cmp(&self, o: &Self) -> Option<Ordering> {
        Some(self.cmp(o))
    }
}
impl Ord for Interval {
    fn cmp(&self, o: &Self) -> Ordering {
        self.excess.total_cmp(&o.excess)
    }
}

/// Largest ratio of a change across an interval to its limit.
fn excess(a: &[f64; 4], b: &[f64; 4], o: &RemeshOptions) -> f64 {
    let lim = [o.dlnp, o.dlnt, o.dx, o.dl];
    (0..4).map(|i| (b[i] - a[i]).abs() / lim[i]).fold(0.0, f64::max)
}

/// The mesh for the next step: `base` with the intervals that the previous
/// model `model` does not resolve bisected. Returns `base` itself when
/// nothing needs refining.
pub fn refine_mesh(model: &StarModel, base: &[f64], opts: &RemeshOptions) -> Vec<f64> {
    let f = Fields::from_model(model);
    let budget = (opts.max_extra * base.len() as f64) as usize;
    let mut heap = BinaryHeap::new();
    let mut done: Vec<(f64, f64)> = Vec::with_capacity(base.len() + budget);
    let vals: Vec<[f64; 4]> = base.iter().map(|&q| f.at(q)).collect();
    for k in 0..base.len() - 1 {
        let e = excess(&vals[k], &vals[k + 1], opts);
        let iv = Interval { lo: base[k], hi: base[k + 1], v_lo: vals[k], v_hi: vals[k + 1], depth: 0, excess: e };
        if e > 1.0 { heap.push(iv) } else { done.push((iv.lo, iv.hi)) }
    }
    let mut inserted = 0;
    while let Some(iv) = heap.pop() {
        let mid = 0.5 * (iv.lo + iv.hi);
        if inserted >= budget || iv.depth >= opts.max_depth || !(mid > iv.lo && mid < iv.hi) {
            done.push((iv.lo, iv.hi));
            continue;
        }
        inserted += 1;
        let vm = f.at(mid);
        for (lo, hi, a, b) in [(iv.lo, mid, iv.v_lo, vm), (mid, iv.hi, vm, iv.v_hi)] {
            let e = excess(&a, &b, opts);
            let child = Interval { lo, hi, v_lo: a, v_hi: b, depth: iv.depth + 1, excess: e };
            if e > 1.0 { heap.push(child) } else { done.push((lo, hi)) }
        }
    }
    if inserted == 0 {
        return base.to_vec();
    }
    let mut q: Vec<f64> = done.iter().map(|d| d.0).collect();
    q.push(*base.last().unwrap());
    q.sort_by(f64::total_cmp);
    q.dedup();
    q
}

/// Boundaries of the cells of the points `m`: 0, the midpoints, m_{n−1}.
fn cell_edges(m: &[f64]) -> Vec<f64> {
    let n = m.len();
    let mut e = Vec::with_capacity(n + 1);
    e.push(0.0);
    for k in 1..n {
        e.push(0.5 * (m[k - 1] + m[k]));
    }
    e.push(m[n - 1]);
    e
}

/// Conservative remap of the composition from the points `m_old` to `m_new`
/// (same total mass): each new cell gets the mass-weighted mean of the old
/// cells it overlaps. X, Y, Z₀, X_C, X_O are remapped; the CNO catalysts are
/// recomputed as f_CNO Z₀ + X_C.
pub fn remap_composition(m_old: &[f64], comps: &[Composition], m_new: &[f64], f_cno: f64) -> Vec<Composition> {
    let (eo, en) = (cell_edges(m_old), cell_edges(m_new));
    let mut out = Vec::with_capacity(m_new.len());
    let mut i = 0;
    for j in 0..m_new.len() {
        let (lo, hi) = (en[j], en[j + 1]);
        let mut acc = [0.0; 5];
        let mut mass = 0.0;
        while i + 1 < eo.len() && eo[i + 1] <= lo {
            i += 1;
        }
        let mut k = i;
        while k < comps.len() && eo[k] < hi {
            let w = (eo[k + 1].min(hi) - eo[k].max(lo)).max(0.0);
            let c = &comps[k];
            for (a, v) in acc.iter_mut().zip([c.x, c.y, c.z, c.x_c, c.x_o]) {
                *a += w * v;
            }
            mass += w;
            k += 1;
        }
        let c = if mass > 0.0 {
            let [x, y, z, x_c, x_o] = acc.map(|a| a / mass);
            Composition { x, y, z, x_cno: Composition::cno_catalysts(f_cno, z, x_c), x_c, x_o }
        } else {
            // A cell of zero width (cannot happen for a strictly increasing mesh).
            comps[i.min(comps.len() - 1)]
        };
        out.push(c);
    }
    out
}

/// Mass of each species in the cells of `m`: ∫X dm, ∫Y dm, ∫Z₀ dm, ∫X_C dm, ∫X_O dm [g].
pub fn species_masses(m: &[f64], comps: &[Composition]) -> [f64; 5] {
    let e = cell_edges(m);
    let mut s = [0.0; 5];
    for (k, c) in comps.iter().enumerate() {
        let dm = e[k + 1] - e[k];
        for (a, v) in s.iter_mut().zip([c.x, c.y, c.z, c.x_c, c.x_o]) {
            *a += dm * v;
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn comp(x: f64, x_c: f64) -> Composition {
        let z = 0.02;
        Composition { x, y: 1.0 - x - z - x_c, z, x_cno: 0.7 * z + x_c, x_c, x_o: 0.0 }
    }

    #[test]
    fn remap_conserves_every_species() {
        let m_old: Vec<f64> = (0..200).map(|i| (i as f64 + 1.0).powf(1.3)).collect();
        let total = *m_old.last().unwrap();
        let comps: Vec<Composition> = (0..200)
            .map(|i| {
                let x = if i < 70 {
                    0.0
                } else if i < 75 {
                    0.14 * (i - 70) as f64
                } else {
                    0.7
                };
                comp(x, if i < 30 { 0.3 } else { 0.0 })
            })
            .collect();
        // A new mesh: some points removed, some inserted, same total mass.
        let mut m_new: Vec<f64> = m_old.iter().step_by(3).copied().collect();
        m_new.extend(m_old.windows(2).skip(60).take(30).map(|w| 0.5 * (w[0] + w[1])));
        m_new.push(total);
        m_new.sort_by(f64::total_cmp);
        m_new.dedup();
        let new = remap_composition(&m_old, &comps, &m_new, 0.7);
        let (a, b) = (species_masses(&m_old, &comps), species_masses(&m_new, &new));
        for i in 0..5 {
            assert!((a[i] - b[i]).abs() <= 1e-12 * total, "species {i}: {} vs {}", a[i], b[i]);
        }
        for c in &new {
            assert!((c.sum() - 1.0).abs() < 1e-12);
        }
        // The same mesh maps onto itself.
        let same = remap_composition(&m_old, &comps, &m_old, 0.7);
        for (u, v) in same.iter().zip(&comps) {
            assert!((u.x - v.x).abs() < 1e-14 && (u.x_c - v.x_c).abs() < 1e-14);
        }
    }
}
