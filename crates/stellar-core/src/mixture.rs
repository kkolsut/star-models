//! Heavy-element mixture (plan 7.1).
//!
//! The opacity tables (OPAL `GN93hz`, Ferguson et al. 2005 `f05.g93`) are
//! computed for the solar mixture of Grevesse & Noels (1993). The same mixture
//! defines the CNO fraction of the burning rates and the mean ion of the
//! metals in the EOS, so all three use one set of numbers.
//!
//! Number and mass fractions of the metals from the header of `GN93hz`
//! (lines 31–60); atomic masses from the IUPAC 2013 standard atomic weights.
//! The averages below are computed from them (see the unit test).

/// Averages over the metals of a mixture.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mixture {
    pub name: &'static str,
    /// Mass fraction of C + N + O in the metals.
    pub f_cno: f64,
    /// Mean atomic mass ⟨A⟩ (by number).
    pub mean_a: f64,
    /// Mean nuclear charge ⟨Z⟩ (by number).
    pub mean_z: f64,
    /// Mean squared charge ⟨Z²⟩ (by number).
    pub mean_z2: f64,
    /// Electrons per atomic mass unit per unit metal mass fraction, ⟨Z⟩/⟨A⟩
    /// (1/2 in the approximation of the analytic EOS).
    pub z_over_a: f64,
}

/// (symbol, Z, A, number fraction among the metals) for GN93.
pub const GN93_ELEMENTS: [(&str, u32, f64, f64); 19] = [
    ("C", 6, 12.011, 0.245518),
    ("N", 7, 14.007, 0.064578),
    ("O", 8, 15.999, 0.512966),
    ("Ne", 10, 20.180, 0.083210),
    ("Na", 11, 22.990, 0.001479),
    ("Mg", 12, 24.305, 0.026308),
    ("Al", 13, 26.982, 0.002042),
    ("Si", 14, 28.086, 0.024552),
    ("P", 15, 30.974, 0.000195),
    ("S", 16, 32.06, 0.011222),
    ("Cl", 17, 35.45, 0.000219),
    ("Ar", 18, 39.948, 0.002291),
    ("K", 19, 39.098, 0.000091),
    ("Ca", 20, 40.078, 0.001586),
    ("Ti", 22, 47.867, 0.000075),
    ("Cr", 24, 51.996, 0.000329),
    ("Mn", 25, 54.938, 0.000170),
    ("Fe", 26, 55.845, 0.021877),
    ("Ni", 28, 58.693, 0.001293),
];

/// Grevesse & Noels (1993), as used by the OPAL and Ferguson tables.
pub const GN93: Mixture = Mixture {
    name: "Grevesse & Noels (1993)",
    f_cno: 0.708_70,
    mean_a: 17.0175,
    mean_z: 8.4429,
    mean_z2: 82.802,
    z_over_a: 0.496_13,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn averages_follow_from_the_table() {
        let n: f64 = GN93_ELEMENTS.iter().map(|e| e.3).sum();
        let a = GN93_ELEMENTS.iter().map(|e| e.2 * e.3).sum::<f64>() / n;
        let z = GN93_ELEMENTS.iter().map(|e| e.1 as f64 * e.3).sum::<f64>() / n;
        let z2 = GN93_ELEMENTS.iter().map(|e| (e.1 as f64).powi(2) * e.3).sum::<f64>() / n;
        let cno = GN93_ELEMENTS[..3].iter().map(|e| e.2 * e.3).sum::<f64>() / (a * n);
        assert!((n - 1.0).abs() < 1e-5);
        assert!((a / GN93.mean_a - 1.0).abs() < 1e-5, "⟨A⟩ = {a}");
        assert!((z / GN93.mean_z - 1.0).abs() < 1e-5, "⟨Z⟩ = {z}");
        assert!((z2 / GN93.mean_z2 - 1.0).abs() < 1e-5, "⟨Z²⟩ = {z2}");
        assert!((cno / GN93.f_cno - 1.0).abs() < 1e-5, "f_CNO = {cno}");
        assert!((z / a / GN93.z_over_a - 1.0).abs() < 1e-5);
        // The CNO mass fraction of the GN93hz header: 0.173285 + 0.053152 + 0.482273.
        assert!((GN93.f_cno - 0.708_71).abs() < 2e-5);
    }
}
