//! Export of models to CSV (profiles) and JSON (parameters + global properties).

use std::fmt::Write as _;

use serde::Serialize;

use crate::constants::{L_SUN, M_SUN, R_SUN};
use crate::structure::model::{ConvectiveZone, Globals, StarModel, StarParams};

/// CSV column names, in output order.
pub const CSV_COLUMNS: &[&str] = &[
    "q",
    "m_msun",
    "r_rsun",
    "r_cm",
    "P_dyn_cm2",
    "T_K",
    "rho_g_cm3",
    "L_lsun",
    "kappa_cm2_g",
    "eps_pp",
    "eps_cno",
    "eps_3a",
    "eps_total",
    "nabla_rad",
    "nabla_ad",
    "nabla",
    "convective",
    "v_conv_cm_s",
    "F_conv_frac",
    "H_P_cm",
    "mu",
    "beta",
    "X",
    "Y",
    "Z",
];

/// Profile as CSV (one row per mesh point, centre to surface).
pub fn profile_csv(model: &StarModel) -> String {
    let p = &model.profile;
    let mut s = String::with_capacity(p.len() * 400);
    s.push_str(&CSV_COLUMNS.join(","));
    s.push('\n');
    for k in 0..p.len() {
        let row = [
            p.q[k],
            p.m[k] / M_SUN,
            p.r[k] / R_SUN,
            p.r[k],
            p.p[k],
            p.t[k],
            p.rho[k],
            p.l[k] / L_SUN,
            p.kappa[k],
            p.eps_pp[k],
            p.eps_cno[k],
            p.eps_3a[k],
            p.eps[k],
            p.nabla_rad[k],
            p.nabla_ad[k],
            p.nabla[k],
            if p.convective[k] { 1.0 } else { 0.0 },
            p.v_conv[k],
            p.f_conv[k],
            p.hp[k],
            p.mu[k],
            p.beta[k],
            p.x[k],
            p.y[k],
            p.z[k],
        ];
        for (i, v) in row.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            if i == 16 {
                let _ = write!(s, "{}", *v as u8);
            } else {
                let _ = write!(s, "{v:.10e}");
            }
        }
        s.push('\n');
    }
    s
}

#[derive(Serialize)]
struct Summary<'a> {
    params: &'a StarParams,
    physics: &'a str,
    globals: &'a Globals,
    convective_zones: &'a [ConvectiveZone],
}

/// Parameters, global properties and convective zones as pretty JSON.
pub fn summary_json(model: &StarModel) -> String {
    serde_json::to_string_pretty(&Summary {
        params: &model.params,
        physics: &model.physics,
        globals: &model.globals,
        convective_zones: &model.zones,
    })
    .expect("model summary is always serialisable")
}
