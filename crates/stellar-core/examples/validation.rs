//! Plan 7.9: comparison with published models (different physics, unlike
//! track V). The published grid is represented by the fitting formulae of
//! the Eggleton-code models (Pols et al. 1998):
//!
//! * ZAMS L and R: Tout, Pols, Eggleton & Han (1996), MNRAS 281, 257 (errors
//!   below 3% in L and 1.2% in R at Z = 0.02);
//! * main-sequence lifetime t_MS = max(x, t_hook/t_BGB) t_BGB: Hurley, Pols &
//!   Tout (2000), MNRAS 315, 543, eqs. 4–7.
//!
//! Coefficients transcribed from the SSE/BSE source as distributed with COSMIC
//! (`src/cosmic/src/zdata.h`, `SSE/SSE_zcnsts.f`, `SSE/SSE_zfuncs.f`), with
//! Z☉ = 0.02 in ζ = log₁₀(Z/0.02).
//!
//! Run: `cargo run --release -p stellar-core --example validation -- [physics] [X] [Z] [alpha] [lifetimes]`
//! with physics `all` or `analytic`; `lifetimes` also evolves 1, 3 and 5 M☉ to
//! the end of the main sequence.

use stellar_core::physics::options::PhysicsOptions;
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Physics, StarModel, StarParams, StopCondition, solve_model,
};

/// Coefficients of zdata.h, xz(1..76) (L_ZAMS 1–35, R_ZAMS 36–76).
const XZ: [f64; 76] = [
    3.970417e-01,
    -3.2913574e-01,
    3.4776688e-01,
    3.7470851e-01,
    9.011915e-02,
    8.527626e+00,
    -2.441225973e+01,
    5.643597107e+01,
    3.706152575e+01,
    5.4562406e+00,
    2.5546e-04,
    -1.23461e-03,
    -2.3246e-04,
    4.5519e-04,
    1.6176e-04,
    5.432889e+00,
    -8.62157806e+00,
    1.344202049e+01,
    1.451584135e+01,
    3.39793084e+00,
    5.563579e+00,
    -1.032345224e+01,
    1.944322980e+01,
    1.897361347e+01,
    4.16903097e+00,
    7.8866060e-01,
    -2.90870942e+00,
    6.54713531e+00,
    4.05606657e+00,
    5.3287322e-01,
    5.86685e-03,
    -1.704237e-02,
    3.872348e-02,
    2.570041e-02,
    3.83376e-03,
    1.715359e+00,
    6.2246212e-01,
    -9.2557761e-01,
    -1.16996966e+00,
    -3.0631491e-01,
    6.597788e+00,
    -4.2450044e-01,
    -1.213339427e+01,
    -1.073509484e+01,
    -2.51487077e+00,
    1.008855000e+01,
    -7.11727086e+00,
    -3.167119479e+01,
    -2.424848322e+01,
    -5.33608972e+00,
    1.012495e+00,
    3.2699690e-01,
    -9.23418e-03,
    -3.876858e-02,
    -4.12750e-03,
    7.490166e-02,
    2.410413e-02,
    7.233664e-02,
    3.040467e-02,
    1.97741e-03,
    1.077422e-02,
    3.082234e+00,
    9.447205e-01,
    -2.15200882e+00,
    -2.49219496e+00,
    -6.3848738e-01,
    1.784778e+01,
    -7.4534569e+00,
    -4.896066856e+01,
    -4.005386135e+01,
    -9.09331816e+00,
    2.2582e-04,
    -1.86899e-03,
    3.88783e-03,
    1.42402e-03,
    -7.671e-05,
];
/// xt(1..31): T_BGB 1–17, t_hook 18–31.
const XT: [f64; 31] = [
    1.593890e+03,
    2.053038e+03,
    1.231226e+03,
    2.327785e+02,
    2.706708e+03,
    1.483131e+03,
    5.772723e+02,
    7.411230e+01,
    1.466143e+02,
    -1.048442e+02,
    -6.795374e+01,
    -1.391127e+01,
    4.141960e-02,
    4.564888e-02,
    2.958542e-02,
    5.571483e-03,
    3.426349e-01,
    1.949814e+01,
    1.758178e+00,
    -6.008212e+00,
    -4.470533e+00,
    4.903830e+00,
    5.212154e-02,
    3.166411e-02,
    -2.750074e-03,
    -2.271549e-03,
    1.312179e+00,
    -3.294936e-01,
    9.231860e-02,
    2.610989e-02,
    8.073972e-01,
];

struct Fits {
    a: [f64; 27],
    x_ms: f64,
}

impl Fits {
    fn new(z: f64) -> Self {
        let lz = (z / 0.02).log10();
        let p5 = |i: usize| XZ[i] + lz * (XZ[i + 1] + lz * (XZ[i + 2] + lz * (XZ[i + 3] + lz * XZ[i + 4])));
        let p4 = |i: usize| XT[i] + lz * (XT[i + 1] + lz * (XT[i + 2] + lz * XT[i + 3]));
        let mut a = [0.0; 27];
        for k in 0..7 {
            a[1 + k] = p5(5 * k);
        }
        // R: msp(8..12) from xz(36..60), msp(13) = xz(61), msp(14..16) from xz(62..76).
        for k in 0..5 {
            a[8 + k] = p5(35 + 5 * k);
        }
        a[13] = XZ[60];
        for k in 0..3 {
            a[14 + k] = p5(61 + 5 * k);
        }
        for k in 0..4 {
            a[17 + k] = p4(4 * k);
        }
        a[21] = XT[16];
        a[22] = p4(17);
        a[23] = XT[21];
        a[24] = p4(22);
        a[25] = p4(26);
        a[26] = XT[30];
        let x_ms = 0.95f64.max((0.95 - 10.0 / 3.0 * (z - 0.01)).max(0.99f64.min(0.98 - 100.0 / 7.0 * (z - 0.001))));
        Self { a, x_ms }
    }

    fn l_zams(&self, m: f64) -> f64 {
        let a = &self.a;
        (a[1] * m.powf(5.5) + a[2] * m.powi(11))
            / (a[3] + m.powi(3) + a[4] * m.powi(5) + a[5] * m.powi(7) + a[6] * m.powi(8) + a[7] * m.powf(9.5))
    }

    fn r_zams(&self, m: f64) -> f64 {
        let a = &self.a;
        let mx = m.sqrt();
        ((a[8] * m.powi(2) + a[9] * m.powi(6)) * mx + a[10] * m.powi(11) + (a[11] + a[12] * mx) * m.powi(19))
            / (a[13] + a[14] * m.powi(2) + (a[15] * m.powi(8) + m.powi(18) + a[16] * m.powi(19)) * mx)
    }

    /// Main-sequence lifetime [Myr].
    fn t_ms(&self, m: f64) -> f64 {
        let a = &self.a;
        let tbgb =
            (a[17] + a[18] * m.powi(4) + a[19] * m.powf(5.5) + m.powi(7)) / (a[20] * m.powi(2) + a[21] * m.powi(7));
        let thook = (1.0 - 0.01 * (a[22] / m.powf(a[23])).max(a[24] + a[25] / m.powf(a[26]))).max(0.5);
        self.x_ms.max(thook) * tbgb
    }
}

fn lifetime(zams: StarModel, physics: &Physics) -> f64 {
    let opts = EvolutionOptions { stop: StopCondition::EndOfMainSequence, pulsation: false, ..Default::default() };
    let mut evo = Evolution::new(zams, opts);
    loop {
        match evo.step(physics, &mut NoControl).expect("evolution") {
            EvolutionStep::Advanced { .. } => {}
            EvolutionStep::Finished(_) => return evo.model.globals.age_yr / 1e6,
        }
    }
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let opts = match args.first().map(String::as_str) {
        Some("analytic") => PhysicsOptions::analytic(),
        _ => PhysicsOptions::all_upgrades(),
    };
    let num = |i: usize, d: f64| args.get(i).map(|s| s.parse().unwrap()).unwrap_or(d);
    let (x, z, alpha) = (num(1, 0.70), num(2, 0.02), num(3, 1.8));
    let lifetimes = args.iter().any(|a| a == "lifetimes");
    let physics = Physics::from_options(&opts);
    let fits = Fits::new(z);
    println!("{}; X = {x}, Z = {z}, α_MLT = {alpha}\n", physics.describe());
    println!(
        "| M [M☉] | L | L (Tout+96) | ΔL | R | R (Tout+96) | ΔR | T_eff [K] | T_c [K] | ρ_c | env. base [r/R] | core [m/M] |"
    );
    println!("|---|---|---|---|---|---|---|---|---|---|---|---|");
    let base = StarParams { x, z, alpha_mlt: alpha, f_cno: stellar_core::mixture::GN93.f_cno, ..StarParams::sun() };
    let mut prev: Option<StarModel> = None;
    let mut models = Vec::new();
    for m in [1.0, 0.8, 0.5, 1.5, 2.0, 3.0, 5.0, 10.0, 20.0] {
        let p = StarParams { mass: m, ..base.clone() };
        // Warm start from the previous mass when it is within a factor of two.
        let warm = prev.as_ref().filter(|w| (w.globals.mass - m).abs() / m < 1.0);
        let model = match solve_model(&p, &physics, warm, &mut NoControl) {
            Ok(model) => model,
            Err(e) => {
                println!("| {m} | failed: {e} |");
                continue;
            }
        };
        let g = &model.globals;
        let (lt, rt) = (fits.l_zams(m), fits.r_zams(m));
        println!(
            "| {m} | {:.4} | {:.4} | {:+.1}% | {:.4} | {:.4} | {:+.1}% | {:.0} | {:.3e} | {:.2} | {} | {:.3} |",
            g.luminosity,
            lt,
            100.0 * (g.luminosity / lt - 1.0),
            g.radius,
            rt,
            100.0 * (g.radius / rt - 1.0),
            g.teff,
            g.tc,
            g.rhoc,
            g.conv_env_base_r.map(|r| format!("{r:.3}")).unwrap_or_else(|| "–".into()),
            g.conv_core_q
        );
        prev = Some(model.clone());
        models.push((m, model));
    }
    if lifetimes {
        println!("\n| M [M☉] | t_MS [Myr] | t_MS (Hurley+00) [Myr] | Δ |");
        println!("|---|---|---|---|");
        for (m, model) in models.into_iter().filter(|(m, _)| [1.0, 3.0, 5.0].contains(m)) {
            let t = lifetime(model, &physics);
            let h = fits.t_ms(m);
            println!("| {m} | {t:.1} | {h:.1} | {:+.1}% |", 100.0 * (t / h - 1.0));
        }
    }
}
