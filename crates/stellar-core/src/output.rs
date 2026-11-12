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
    "eps_grav",
    "nabla_rad",
    "nabla_ad",
    "nabla",
    "convective",
    "v_conv_cm_s",
    "F_conv_frac",
    "H_P_cm",
    "mu",
    "beta",
    "Gamma1",
    "X",
    "Y",
    "Z",
    "c_s_cm_s",
    "N2_s-2",
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
            p.eps_grav[k],
            p.nabla_rad[k],
            p.nabla_ad[k],
            p.nabla[k],
            if p.convective[k] { 1.0 } else { 0.0 },
            p.v_conv[k],
            p.f_conv[k],
            p.hp[k],
            p.mu[k],
            p.beta[k],
            p.gamma1[k],
            p.x[k],
            p.y[k],
            p.z[k],
            p.cs.get(k).copied().unwrap_or(f64::NAN),
            p.n2.get(k).copied().unwrap_or(f64::NAN),
        ];
        for (i, v) in row.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            if CSV_COLUMNS[i] == "convective" {
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

/// Evolutionary track as CSV: one row per model.
pub fn track_csv(models: &[&StarModel]) -> String {
    let mut s = String::from(
        "age_yr,phase,L_lsun,Teff_K,R_rsun,log_g,Tc_K,rhoc_g_cm3,Xc,Yc,he_core_q,conv_core_q,mixed_core_q,conv_env_base_r,\
         frac_pp,frac_cno,frac_3a,frac_grav,P_fund_days,Q_days,omega_over_ac,pulsator,nu_max_uHz,delta_nu_uHz,\
         delta_nu_asym_uHz,delta_Pi1_s,v_z_m_s,screening_max,energy_residual,energy_residual_comp,fuel_ratio\n",
    );
    for m in models {
        let g = &m.globals;
        let (period, q, ratio_ac) = m
            .pulsation
            .as_ref()
            .map(|p| (p.period_days(), p.q_days, p.omega_over_ac))
            .unwrap_or((f64::NAN, f64::NAN, f64::NAN));
        let check = g.energy_check();
        let opt = |v: Option<f64>| v.map(|v| format!("{v:.4e}")).unwrap_or_default();
        let pulsator = crate::pulsation::classify(g.teff, g.luminosity).map(|p| p.label()).unwrap_or("");
        let _ = writeln!(
            s,
            "{:.6e},{:?},{:.6e},{:.2},{:.6e},{:.4},{:.6e},{:.6e},{:.6},{:.6},{:.6},{:.6},{:.6},{},{:.5},{:.5},{:.5},{:.5},{:.6e},{:.5},{:.4},{},{:.4},{:.4},{:.4},{},{:.2},{:.4},{},{},{}",
            g.age_yr,
            g.phase,
            g.luminosity,
            g.teff,
            g.radius,
            g.log_g,
            g.tc,
            g.rhoc,
            g.x_c,
            g.y_c,
            g.he_core_q,
            g.conv_core_q,
            g.mixed_core_q,
            g.conv_env_base_r.map(|r| format!("{r:.5}")).unwrap_or_default(),
            g.frac_pp,
            g.frac_cno,
            g.frac_3a,
            g.frac_grav,
            period,
            q,
            ratio_ac,
            pulsator,
            g.nu_max_uhz,
            g.delta_nu_uhz,
            g.delta_nu_asym_uhz,
            crate::seismic::main_cavity(&g.g_cavities).map(|c| format!("{:.2}", c.delta_pi1_s)).unwrap_or_default(),
            g.grav_redshift_m_s,
            g.screening.max_exponent,
            opt(check.map(|c| c.residual)),
            opt(check.map(|c| c.composition)),
            check.map(|c| format!("{:.6}", c.fuel_ratio)).unwrap_or_default(),
        );
    }
    s
}
