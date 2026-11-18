//! Tabulated Rosseland mean opacities (plan 7.1–7.4) and electron conduction
//! (plan 7.5, brought forward with the tables).
//!
//! Data, embedded in the binary (`crates/stellar-core/data/`, see the README
//! there for sources, checksums and the licence status):
//!
//! * interior, 3.75 ≤ log T ≤ 8.70: OPAL type 1 tables `GN93hz`
//!   (Iglesias & Rogers 1996), 126 (X, Z) compositions;
//! * low temperature, 2.70 ≤ log T ≤ 4.50: Ferguson et al. (2005), set
//!   `f05.g93`, 155 (X, Z) compositions, molecules and grains included;
//! * electron conduction: the table `condtab21wd.dat` of A. Y. Potekhin
//!   (Cassisi et al. 2007, with the 2021 update of Cassisi et al. 2021;
//!   reference Potekhin, Pons & Page 2015).
//!
//! All radiative tables use the Grevesse & Noels (1993) heavy-element mixture
//! ([`crate::mixture`]).
//!
//! log κ is interpolated in (log T, log R), R = ρ / T₆³, by bicubic Hermite
//! interpolation with node derivatives from local parabolas: C¹, and exact at
//! the nodes. Between tables it is interpolated in X and in ζ = ln(Z + 10⁻⁵)
//! with the same C¹ cubic (composition is not a solver unknown, but evolution
//! changes it smoothly). The two sets are blended in log κ with a C¹ weight
//! over 3.75 ≤ log T ≤ 4.0. Outside the region where a table has data (the
//! high-density corner of the OPAL tables, the hydrogen-free low-R corner,
//! beyond the grid edges) the table is continued C¹ by its tangent plane and
//! blended over 0.5 dex into the analytic formulas of [`super::opacity`]:
//! nothing is extrapolated silently beyond that band.
//!
//! Electron conduction adds as a conductance, 1/κ = 1/κ_rad + 1/κ_cond, with
//! κ_cond = 16σT³/(3ρK) (Cassisi et al. 2007, eq. 1) and K from Potekhin's
//! table for a pure plasma of charge Z_eff = ⟨Z²⟩^{1/2} averaged over the ions
//! by number (Cassisi et al. 2007, §2.2.3), interpolated in log Z_eff.
//!
//! Compositions outside the tables are clamped to the nearest one: X to the
//! tabulated range at each Z, Z to at most 0.1. After helium burning the
//! metals are mostly carbon and oxygen, which a scaled-solar table does not
//! describe (OPAL type 2 tables would; plan 7.4 needs the composition vector
//! of plan 9.5). [`TabulatedOpacity::composition_note`] reports the clamp.

use std::sync::OnceLock;

use super::Opacity;
use super::interp::{Hermite1, Hermite2, hermite_weights, smoothstep};
use super::opacity::AnalyticOpacity;
use crate::composition::Composition;
use crate::constants::SIGMA_SB;
use crate::mixture;

const OPAL_GN93: &str = include_str!("../../data/opal/GN93hz");
const POTEKHIN: &str = include_str!("../../data/potekhin/condtab21wd.dat");

/// The Ferguson et al. (2005) files, embedded one by one.
macro_rules! ferguson {
    ($($name:literal),* $(,)?) => {
        &[$(include_str!(concat!("../../data/ferguson2005_gn93/", $name))),*]
    };
}
const FERGUSON_GN93: &[&str] = ferguson!(
    "g0.0.tron",
    "g0.00001.tron",
    "g0.00003.tron",
    "g0.0001.tron",
    "g0.0003.tron",
    "g0.001.tron",
    "g0.002.tron",
    "g0.004.tron",
    "g0.01.tron",
    "g0.02.tron",
    "g0.03.tron",
    "g0.04.tron",
    "g0.05.tron",
    "g0.06.tron",
    "g0.08.tron",
    "g0.1.tron",
    "g1.0.tron",
    "g1.00001.tron",
    "g1.00003.tron",
    "g1.0001.tron",
    "g1.0003.tron",
    "g1.001.tron",
    "g1.002.tron",
    "g1.004.tron",
    "g1.01.tron",
    "g1.02.tron",
    "g1.03.tron",
    "g1.04.tron",
    "g1.05.tron",
    "g1.06.tron",
    "g1.08.tron",
    "g1.1.tron",
    "g10.0.tron",
    "g2.0.tron",
    "g2.00001.tron",
    "g2.00003.tron",
    "g2.0001.tron",
    "g2.0003.tron",
    "g2.001.tron",
    "g2.002.tron",
    "g2.004.tron",
    "g2.01.tron",
    "g2.02.tron",
    "g2.03.tron",
    "g2.04.tron",
    "g2.05.tron",
    "g2.06.tron",
    "g2.08.tron",
    "g2.1.tron",
    "g35.0.tron",
    "g35.00001.tron",
    "g35.00003.tron",
    "g35.0001.tron",
    "g35.0003.tron",
    "g35.001.tron",
    "g35.002.tron",
    "g35.004.tron",
    "g35.01.tron",
    "g35.02.tron",
    "g35.03.tron",
    "g35.04.tron",
    "g35.05.tron",
    "g35.06.tron",
    "g35.08.tron",
    "g35.1.tron",
    "g5.0.tron",
    "g5.00001.tron",
    "g5.00003.tron",
    "g5.0001.tron",
    "g5.0003.tron",
    "g5.001.tron",
    "g5.002.tron",
    "g5.004.tron",
    "g5.01.tron",
    "g5.02.tron",
    "g5.03.tron",
    "g5.04.tron",
    "g5.05.tron",
    "g5.06.tron",
    "g5.08.tron",
    "g5.1.tron",
    "g7.0.tron",
    "g7.00001.tron",
    "g7.00003.tron",
    "g7.0001.tron",
    "g7.0003.tron",
    "g7.001.tron",
    "g7.002.tron",
    "g7.004.tron",
    "g7.01.tron",
    "g7.02.tron",
    "g7.03.tron",
    "g7.04.tron",
    "g7.05.tron",
    "g7.06.tron",
    "g7.08.tron",
    "g7.1.tron",
    "g8.0.tron",
    "g8.00001.tron",
    "g8.00003.tron",
    "g8.0001.tron",
    "g8.0003.tron",
    "g8.001.tron",
    "g8.002.tron",
    "g8.004.tron",
    "g8.01.tron",
    "g8.02.tron",
    "g8.03.tron",
    "g8.04.tron",
    "g8.05.tron",
    "g8.06.tron",
    "g8.08.tron",
    "g8.1.tron",
    "g9.0.tron",
    "g9.00001.tron",
    "g9.00003.tron",
    "g9.0001.tron",
    "g9.0003.tron",
    "g9.001.tron",
    "g9.002.tron",
    "g9.004.tron",
    "g9.01.tron",
    "g9.02.tron",
    "g9.03.tron",
    "g9.04.tron",
    "g9.05.tron",
    "g9.06.tron",
    "g9.08.tron",
    "g9.1.tron",
    "g92.08.tron",
    "g94.06.tron",
    "g95.0.tron",
    "g95.00001.tron",
    "g95.00003.tron",
    "g95.0001.tron",
    "g95.0003.tron",
    "g95.001.tron",
    "g95.002.tron",
    "g95.004.tron",
    "g95.01.tron",
    "g95.02.tron",
    "g95.03.tron",
    "g95.04.tron",
    "g95.05.tron",
    "g96.04.tron",
    "g97.03.tron",
    "g98.02.tron",
    "g99.01.tron",
    "g996.004.tron",
    "g998.002.tron",
    "g999.001.tron",
    "g9997.0003.tron",
    "g9999.0001.tron",
    "g99997.00003.tron",
    "g99999.00001.tron",
);

/// Width [dex] of the band over which a table is blended into the analytic
/// formulas outside the region where it has data.
pub const EDGE_BAND: f64 = 0.5;
/// log T range of the blend between the low-temperature and the OPAL tables.
pub const BLEND_LOG_T: (f64, f64) = (3.75, 4.0);
/// Offset in ζ = ln(Z + Z_OFFSET), the variable of the interpolation in Z.
const Z_OFFSET: f64 = 1e-5;
/// Largest tabulated metallicity.
pub const Z_MAX: f64 = 0.1;

/// One composition: log κ on a (log T, log R) grid.
#[derive(Clone, Debug)]
pub struct KappaTable {
    pub x: f64,
    pub z: f64,
    grid: Hermite2,
    /// First and last log R with data in each log T row, interpolated in log T.
    r_lo: Hermite1,
    r_hi: Hermite1,
    /// Number of entries without data (filled for the interpolation).
    pub missing: usize,
}

impl KappaTable {
    /// Build from rows (log T ascending) of optional values on the log R grid.
    fn new(x: f64, z: f64, log_t: Vec<f64>, log_r: Vec<f64>, rows: Vec<Vec<Option<f64>>>) -> Self {
        let ny = log_r.len();
        let mut f = Vec::with_capacity(log_t.len() * ny);
        let (mut lo, mut hi) = (Vec::new(), Vec::new());
        let mut missing = 0;
        for row in &rows {
            let valid: Vec<usize> = (0..ny).filter(|&j| row.get(j).copied().flatten().is_some()).collect();
            assert!(valid.len() >= 2, "opacity table row with fewer than two values");
            let (a, b) = (valid[0], *valid.last().unwrap());
            lo.push(log_r[a]);
            hi.push(log_r[b]);
            let v = |j: usize| row[j].unwrap();
            for j in 0..ny {
                let val = match row.get(j).copied().flatten() {
                    Some(v) => v,
                    None => {
                        missing += 1;
                        if j < a {
                            v(a) + (v(a + 1) - v(a)) * (log_r[j] - log_r[a]) / (log_r[a + 1] - log_r[a])
                        } else if j > b {
                            v(b) + (v(b) - v(b - 1)) * (log_r[j] - log_r[b]) / (log_r[b] - log_r[b - 1])
                        } else {
                            // Interior gap: linear between the neighbouring values.
                            let l = (0..j).rev().find(|&k| row[k].is_some()).unwrap();
                            let r = (j + 1..ny).find(|&k| row[k].is_some()).unwrap();
                            v(l) + (v(r) - v(l)) * (log_r[j] - log_r[l]) / (log_r[r] - log_r[l])
                        }
                    }
                };
                f.push(val);
            }
        }
        Self {
            x,
            z,
            r_lo: Hermite1::new(log_t.clone(), lo),
            r_hi: Hermite1::new(log_t.clone(), hi),
            grid: Hermite2::new(log_t, log_r, f),
            missing,
        }
    }

    /// log κ at (log T, log R) and the weight of the table (1 where it has
    /// data, falling C¹ to 0 at [`EDGE_BAND`] outside).
    pub fn eval(&self, lt: f64, lr: f64) -> (f64, f64) {
        let (f, _, _) = self.grid.eval(lt, lr);
        let t0 = self.grid.x[0];
        let t1 = *self.grid.x.last().unwrap();
        let s = |d: f64| smoothstep((d + EDGE_BAND) / EDGE_BAND);
        let w = s(lt - t0) * s(t1 - lt) * s(lr - self.r_lo.eval(lt).0) * s(self.r_hi.eval(lt).0 - lr);
        (f, w)
    }

    /// log T and log R grid.
    pub fn axes(&self) -> (&[f64], &[f64]) {
        (&self.grid.x, &self.grid.y)
    }

    /// Tabulated value at node (i, j) (after filling).
    pub fn node(&self, i: usize, j: usize) -> f64 {
        self.grid.node(i, j)
    }
}

/// A set of tables for many (X, Z), interpolated between compositions.
#[derive(Clone, Debug)]
pub struct TableSet {
    /// ζ = ln(Z + Z_OFFSET) of the distinct metallicities, ascending.
    zeta: Vec<f64>,
    /// For each metallicity: the tables sorted by X.
    by_z: Vec<Vec<KappaTable>>,
}

impl TableSet {
    fn new(mut tables: Vec<KappaTable>) -> Self {
        tables.sort_by(|a, b| a.z.total_cmp(&b.z).then(a.x.total_cmp(&b.x)));
        let mut by_z: Vec<Vec<KappaTable>> = Vec::new();
        for t in tables {
            match by_z.last_mut() {
                Some(g) if (g[0].z - t.z).abs() < 1e-9 => g.push(t),
                _ => by_z.push(vec![t]),
            }
        }
        let zeta = by_z.iter().map(|g| (g[0].z + Z_OFFSET).ln()).collect();
        Self { zeta, by_z }
    }

    /// Number of tables.
    pub fn len(&self) -> usize {
        self.by_z.iter().map(Vec::len).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn tables(&self) -> impl Iterator<Item = &KappaTable> {
        self.by_z.iter().flatten()
    }

    /// Whether (X, Z) lies inside the tabulated compositions.
    pub fn covers(&self, x: f64, z: f64) -> bool {
        let zmax = self.by_z.last().unwrap()[0].z;
        if !(z <= zmax + 1e-12) {
            return false;
        }
        // X range of the two bracketing metallicities.
        let zeta = (z + Z_OFFSET).ln();
        let i = super::interp::cell(&self.zeta, zeta);
        [i, (i + 1).min(self.by_z.len() - 1)].iter().all(|&k| {
            let g = &self.by_z[k];
            x >= g[0].x - 1e-12 && x <= g.last().unwrap().x + 1e-12
        })
    }

    /// log κ and table weight at (X, Z, log T, log R).
    pub fn eval(&self, x: f64, z: f64, lt: f64, lr: f64) -> (f64, f64) {
        let zeta = (z.clamp(0.0, Z_MAX) + Z_OFFSET).ln();
        let (mut f, mut w) = (0.0, 0.0);
        for (iz, wz) in hermite_weights(&self.zeta, zeta.clamp(self.zeta[0], *self.zeta.last().unwrap())) {
            let g = &self.by_z[iz];
            let xs: Vec<f64> = g.iter().map(|t| t.x).collect();
            let xc = x.clamp(xs[0], *xs.last().unwrap());
            for (ix, wx) in hermite_weights(&xs, xc) {
                let (fi, wi) = g[ix].eval(lt, lr);
                f += wz * wx * fi;
                w += wz * wx * wi;
            }
        }
        (f, w.clamp(0.0, 1.0))
    }
}

/// Potekhin's conductivity table: log K [erg cm⁻¹ s⁻¹ K⁻¹] on (log ρ, log T)
/// for pure plasmas of ion charge Z.
#[derive(Clone, Debug)]
pub struct ConductionTable {
    log_z: Vec<f64>,
    grids: Vec<Hermite2>,
}

impl ConductionTable {
    fn parse(text: &str) -> Self {
        let mut lines = text.lines().skip(1).filter(|l| !l.trim().is_empty());
        let (mut log_z, mut grids) = (Vec::new(), Vec::new());
        while let Some(head) = lines.next() {
            let h: Vec<f64> = head.split_whitespace().map(|v| v.parse().unwrap()).collect();
            let (z, log_t) = (h[0], h[1..].to_vec());
            let mut log_rho = Vec::new();
            let mut vals = Vec::new();
            for _ in 0..64 {
                let r: Vec<f64> = lines.next().unwrap().split_whitespace().map(|v| v.parse().unwrap()).collect();
                log_rho.push(r[0]);
                vals.extend_from_slice(&r[1..]);
            }
            log_z.push(z.log10());
            grids.push(Hermite2::new(log_rho, log_t, vals));
        }
        Self { log_z, grids }
    }

    /// log₁₀ of the electron thermal conductivity for charge `z_ion`.
    pub fn log_k(&self, z_ion: f64, log_rho: f64, log_t: f64) -> f64 {
        let lz = z_ion.log10().clamp(self.log_z[0], *self.log_z.last().unwrap());
        hermite_weights(&self.log_z, lz).into_iter().map(|(i, w)| w * self.grids[i].eval(log_rho, log_t).0).sum()
    }

    /// Conductive opacity κ_cond = 16σT³/(3ρK) [cm² g⁻¹].
    pub fn kappa(&self, z_ion: f64, rho: f64, t: f64) -> f64 {
        let lk = self.log_k(z_ion, rho.log10(), t.log10());
        16.0 * SIGMA_SB * t.powi(3) / (3.0 * rho) / 10f64.powf(lk)
    }
}

/// Effective ion charge ⟨Z²⟩^{1/2}, averaged over the ions by number
/// (Cassisi et al. 2007, §2.2.3), for H, He and the metal mixture.
pub fn effective_charge(c: &Composition) -> f64 {
    let n_h = c.x;
    let n_he = c.y / 4.0;
    let n_z = c.z / mixture::GN93.mean_a;
    let n = n_h + n_he + n_z;
    if !(n > 0.0) {
        return 1.0;
    }
    ((n_h + 4.0 * n_he + mixture::GN93.mean_z2 * n_z) / n).sqrt()
}

fn parse_opal(text: &str) -> Vec<KappaTable> {
    let lines: Vec<&str> = text.lines().collect();
    let mut tables = Vec::new();
    let mut k = 0;
    // Skip the summary; tables start after the "Tables" banner.
    while k < lines.len() && !lines[k].contains("************************************ Tables") {
        k += 1;
    }
    while k < lines.len() {
        let l = lines[k];
        if l.starts_with("TABLE #") {
            let x = field(l, "X=");
            let z = field(l, "Z=");
            // Find the header row "logT ...".
            let mut h = k + 1;
            while !lines[h].trim_start().starts_with("logT") {
                h += 1;
            }
            let log_r: Vec<f64> = lines[h].split_whitespace().skip(1).map(|v| v.parse().unwrap()).collect();
            let mut log_t = Vec::new();
            let mut rows = Vec::new();
            let mut r = h + 1;
            while r < lines.len() && !lines[r].starts_with("TABLE #") {
                let toks: Vec<&str> = lines[r].split_whitespace().collect();
                if !toks.is_empty() {
                    log_t.push(toks[0].parse::<f64>().unwrap());
                    let mut row: Vec<Option<f64>> =
                        toks[1..].iter().map(|v| v.parse::<f64>().ok().filter(|v| (v - 9.999).abs() > 1e-9)).collect();
                    row.resize(log_r.len(), None);
                    rows.push(row);
                }
                r += 1;
            }
            tables.push(KappaTable::new(x, z, log_t, log_r.clone(), rows));
            k = r;
        } else {
            k += 1;
        }
    }
    tables
}

fn field(line: &str, key: &str) -> f64 {
    let i = line.find(key).unwrap() + key.len();
    line[i..].split_whitespace().next().unwrap().parse().unwrap()
}

/// Ferguson et al. (2005) file: header with X and Z, the log R row, then rows
/// of log T (descending) and 19 values in fixed columns of 7 characters (the
/// values can touch, so they are split by column, not by white space).
fn parse_ferguson(text: &str) -> KappaTable {
    let mut lines = text.lines();
    let head = lines.next().unwrap();
    let x = field(head, "X=");
    let z = field(head, "Z=");
    let mut log_r = Vec::new();
    let mut rows: Vec<(f64, Vec<Option<f64>>)> = Vec::new();
    for l in lines {
        let t = l.trim_start();
        if let Some(rest) = t.strip_prefix("log T") {
            log_r = rest.split_whitespace().map(|v| v.parse().unwrap()).collect();
            continue;
        }
        if t.is_empty() || log_r.is_empty() || !t.chars().next().unwrap().is_ascii_digit() {
            continue;
        }
        let lt: f64 = l[..6].trim().parse().unwrap();
        let vals = (0..log_r.len())
            .map(|j| l.get(6 + 7 * j..6 + 7 * (j + 1)).and_then(|s| s.trim().parse::<f64>().ok()))
            .collect();
        rows.push((lt, vals));
    }
    rows.sort_by(|a, b| a.0.total_cmp(&b.0));
    let log_t = rows.iter().map(|r| r.0).collect();
    KappaTable::new(x, z, log_t, log_r, rows.into_iter().map(|r| r.1).collect())
}

/// The parsed data, shared by every [`TabulatedOpacity`].
#[derive(Debug)]
pub struct OpacityData {
    pub opal: TableSet,
    pub low_t: TableSet,
    pub conduction: ConductionTable,
}

static DATA: OnceLock<OpacityData> = OnceLock::new();

/// Parse the embedded tables (once per process).
pub fn data() -> &'static OpacityData {
    DATA.get_or_init(|| OpacityData {
        opal: TableSet::new(parse_opal(OPAL_GN93)),
        low_t: TableSet::new(FERGUSON_GN93.iter().map(|t| parse_ferguson(t)).collect()),
        conduction: ConductionTable::parse(POTEKHIN),
    })
}

/// OPAL + Ferguson et al. (2005) radiative opacities (GN93 mixture), analytic
/// formulas outside the tables, optional electron conduction.
#[derive(Debug, Clone, Copy)]
pub struct TabulatedOpacity {
    /// Add electron conduction (Potekhin's table).
    pub conduction: bool,
}

impl Default for TabulatedOpacity {
    fn default() -> Self {
        Self { conduction: true }
    }
}

/// Parts of a tabulated opacity, for diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KappaParts {
    /// Radiative opacity [cm² g⁻¹] (tables blended with the analytic formulas).
    pub radiative: f64,
    /// Conductive opacity [cm² g⁻¹] (∞ without conduction).
    pub conductive: f64,
    /// Weight of the tables in the radiative opacity (1 inside the tables).
    pub table_weight: f64,
    /// Weight of the OPAL set in the blend with the low-temperature set.
    pub opal_weight: f64,
}

impl TabulatedOpacity {
    pub fn parts(&self, rho: f64, t: f64, c: &Composition) -> KappaParts {
        let d = data();
        let lt = t.log10();
        let lr = rho.log10() - 3.0 * (lt - 6.0);
        let wo = smoothstep((lt - BLEND_LOG_T.0) / (BLEND_LOG_T.1 - BLEND_LOG_T.0));
        let (mut lk, mut w) = (0.0, 0.0);
        if wo > 0.0 {
            let (f, wt) = d.opal.eval(c.x, c.z, lt, lr);
            lk += wo * f;
            w += wo * wt;
        }
        if wo < 1.0 {
            let (f, wt) = d.low_t.eval(c.x, c.z, lt, lr);
            lk += (1.0 - wo) * f;
            w += (1.0 - wo) * wt;
        }
        let radiative = if w >= 1.0 {
            10f64.powf(lk)
        } else {
            let la = AnalyticOpacity.kappa(rho, t, c).log10();
            10f64.powf(w * lk + (1.0 - w) * la)
        };
        let conductive = if self.conduction { d.conduction.kappa(effective_charge(c), rho, t) } else { f64::INFINITY };
        KappaParts { radiative, conductive, table_weight: w, opal_weight: wo }
    }

    /// Note for the app when the composition lies outside the tables (X or Z
    /// clamped, or Z > 0.1 as after helium burning); `None` inside.
    pub fn composition_note(c: &Composition) -> Option<&'static str> {
        let d = data();
        if c.z > Z_MAX + 1e-12 {
            Some(
                "Z above 0.1 (carbon and oxygen from helium burning): opacity taken from the Z = 0.1 \
                 scaled-solar tables; results qualitative (type 2 C/O tables not used)",
            )
        } else if !d.opal.covers(c.x, c.z) || !d.low_t.covers(c.x, c.z) {
            Some("composition outside the opacity tables: clamped to the nearest tabulated X, Z")
        } else {
            None
        }
    }
}

impl Opacity for TabulatedOpacity {
    fn name(&self) -> &'static str {
        if self.conduction {
            "OPAL GN93 + Ferguson et al. 2005 tables, Potekhin conduction"
        } else {
            "OPAL GN93 + Ferguson et al. 2005 tables (no conduction)"
        }
    }

    fn kappa(&self, rho: f64, t: f64, c: &Composition) -> f64 {
        let p = self.parts(rho, t, c);
        if p.conductive.is_finite() { 1.0 / (1.0 / p.radiative + 1.0 / p.conductive) } else { p.radiative }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn comp(x: f64, z: f64) -> Composition {
        Composition::new(x, z, mixture::GN93.f_cno).unwrap()
    }

    #[test]
    fn tables_parse() {
        let d = data();
        assert_eq!(d.opal.len(), 126);
        assert_eq!(d.low_t.len(), 155);
        let t = d.opal.tables().find(|t| t.x == 0.7 && t.z == 0.02).unwrap();
        let (lt, lr) = t.axes();
        assert_eq!(lt.len(), 70);
        assert_eq!(lr.len(), 19);
        // Table 3 of Iglesias & Rogers (1996) / the GN93hz file: X = 0.7, Z = 0.02.
        let f = d.low_t.tables().find(|t| t.x == 0.7 && t.z == 0.02).unwrap();
        assert_eq!(f.axes().0.len(), 85);
    }

    /// Every tabulated node is reproduced exactly by the interpolation, for
    /// compositions on the grid (spot checks against the file contents).
    #[test]
    fn nodes_reproduced() {
        let d = data();
        // GN93hz, table X = 0.7, Z = 0.02, log T = 3.75 and 7.00 rows (file lines).
        let cases = [(0.7, 0.02, 3.75, -8.0), (0.7, 0.02, 7.0, -3.0), (0.35, 0.001, 5.5, -2.5), (0.0, 0.02, 6.0, 0.5)];
        for (x, z, lt, lr) in cases {
            let t = d.opal.tables().find(|t| t.x == x && t.z == z).unwrap();
            let (ts, rs) = t.axes();
            let i = ts.iter().position(|v| (v - lt).abs() < 1e-9).unwrap();
            let j = rs.iter().position(|v| (v - lr).abs() < 1e-9).unwrap();
            let want = t.node(i, j);
            let (got, w) = d.opal.eval(x, z, lt, lr);
            assert!((got - want).abs() < 1e-12, "({x}, {z}, {lt}, {lr}): {got} vs {want}");
            assert_eq!(w, 1.0);
        }
        // Ferguson node.
        let t = d.low_t.tables().find(|t| t.x == 0.7 && t.z == 0.02).unwrap();
        let (got, _) = d.low_t.eval(0.7, 0.02, 3.5, -3.0);
        let (ts, rs) = t.axes();
        let i = ts.iter().position(|v| (v - 3.5).abs() < 1e-9).unwrap();
        let j = rs.iter().position(|v| (v + 3.0).abs() < 1e-9).unwrap();
        assert!((got - t.node(i, j)).abs() < 1e-12);
    }

    /// Spot values read from the files by eye.
    #[test]
    fn file_values() {
        let d = data();
        // GN93hz: X = 0, Z = 0, log T = 3.75, log R = −4.0 → −6.025.
        assert!((d.opal.eval(0.0, 0.0, 3.75, -4.0).0 + 6.025).abs() < 1e-12);
        // f05.g93 g0.02.tron (X = 0, Z = 0.02): log T = 4.5, log R = −8 → −0.644.
        assert!((d.low_t.eval(0.0, 0.02, 4.5, -8.0).0 + 0.644).abs() < 1e-12);
    }

    /// C¹: log κ and its derivatives are continuous across cell edges and
    /// across the blend, and the derivatives equal finite differences.
    #[test]
    fn continuity_and_derivatives() {
        let op = TabulatedOpacity { conduction: true };
        let c = comp(0.7, 0.02);
        let lk = |lrho: f64, lt: f64| op.kappa(10f64.powf(lrho), 10f64.powf(lt), &c).log10();
        // Points on cell edges (log T = 5.0 is a node; log R = −3 at ρ = 10^(−3 + 3(lt − 6))).
        for &(lt, dlr) in &[(5.0, 0.0), (3.75, 0.1), (4.0, -0.2), (3.875, 0.3), (6.0, 0.0), (7.3, 0.25)] {
            let lrho = -3.0 + dlr + 3.0 * (lt - 6.0);
            let h = 1e-5;
            let e = 1e-9;
            let d_lt = |l: f64| (lk(lrho, l + h) - lk(lrho, l - h)) / (2.0 * h);
            let d_lr = |l: f64| (lk(l + h, lt) - lk(l - h, lt)) / (2.0 * h);
            let jump_t = (d_lt(lt + e) - d_lt(lt - e)).abs();
            let jump_r = (d_lr(lrho + e) - d_lr(lrho - e)).abs();
            assert!(jump_t < 1e-4 && jump_r < 1e-4, "lt = {lt}: jumps {jump_t:.1e} {jump_r:.1e}");
            assert!((lk(lrho, lt + e) - lk(lrho, lt - e)).abs() < 1e-6);
        }
    }

    #[test]
    fn solar_interior_value() {
        // Centre of a solar model: T = 1.57e7 K, ρ = 150 g/cm³, X = 0.35, Z = 0.02 → κ ≈ 1.2 (OPAL).
        let op = TabulatedOpacity { conduction: true };
        let k = op.kappa(150.0, 1.57e7, &comp(0.35, 0.02));
        assert!(k > 0.9 && k < 1.6, "κ = {k}");
        let p = op.parts(150.0, 1.57e7, &comp(0.35, 0.02));
        assert!(p.radiative / p.conductive < 0.01, "conduction negligible at the solar centre");
    }

    #[test]
    fn conduction_dominates_in_a_degenerate_helium_core() {
        let op = TabulatedOpacity { conduction: true };
        let c = comp(0.0, 0.02);
        let p = op.parts(1e6, 5e7, &c);
        assert!(p.conductive < 0.1 * p.radiative, "κ_cond = {}, κ_rad = {}", p.conductive, p.radiative);
        // Effective charge of helium.
        assert!((effective_charge(&Composition::new(0.0, 0.0, 0.0).unwrap()) - 2.0).abs() < 1e-12);
    }

    #[test]
    fn analytic_fallback_outside() {
        let op = TabulatedOpacity { conduction: false };
        let c = comp(0.7, 0.02);
        // log T = 9.5: far above the tables.
        let k = op.parts(1.0, 10f64.powf(9.5), &c);
        assert_eq!(k.table_weight, 0.0);
        assert!((k.radiative / AnalyticOpacity.kappa(1.0, 10f64.powf(9.5), &c) - 1.0).abs() < 1e-12);
        // Positive and finite on a wide grid.
        for lt in 25..=95 {
            for lrho in -14..=8 {
                let k = op.kappa(10f64.powi(lrho), 10f64.powf(lt as f64 / 10.0), &c);
                assert!(k.is_finite() && k > 0.0);
            }
        }
    }

    /// Plan 7.4: hydrogen-free tables are used for helium cores, continuously
    /// in X → 0; compositions beyond the tables are clamped and reported.
    #[test]
    fn hydrogen_free_and_clamping() {
        let d = data();
        for (lt, lr) in [(7.5, -2.0), (6.0, -4.0), (5.0, -3.5)] {
            let (k0, _) = d.opal.eval(0.0, 0.02, lt, lr);
            let (k1, _) = d.opal.eval(1e-6, 0.02, lt, lr);
            assert!((k1 - k0).abs() < 1e-5, "X → 0 at ({lt}, {lr}): {k1} vs {k0}");
        }
        assert!(TabulatedOpacity::composition_note(&comp(0.7, 0.02)).is_none());
        assert!(TabulatedOpacity::composition_note(&comp(0.0, 0.02)).is_none());
        let co_rich = Composition { x: 0.0, y: 0.5, z: 0.5, x_cno: 0.5 };
        assert!(TabulatedOpacity::composition_note(&co_rich).is_some());
        // Clamped: Z = 0.5 gives the Z = 0.1 value.
        let (a, _) = d.opal.eval(0.0, 0.5, 7.0, -3.0);
        let (b, _) = d.opal.eval(0.0, 0.1, 7.0, -3.0);
        assert_eq!(a, b);
    }
}
