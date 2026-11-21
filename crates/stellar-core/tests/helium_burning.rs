//! The composition vector (H, He, Z₀, C, O) and ¹²C(α,γ)¹⁶O (plan 9.5), on a
//! pure-helium star: a 1 M☉ helium main-sequence star burns helium in its
//! convective core from the start, so helium burning can be followed without
//! the red-giant phases.

use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Phase, Physics, StarModel, StarParams, StopCondition,
    solve_model,
};

fn helium_star() -> (Vec<StarModel>, String) {
    let physics = Physics::analytic();
    let params = StarParams { mass: 1.0, x: 0.0, z: 0.02, ..StarParams::sun() };
    let zams = solve_model(&params, &physics, None, &mut NoControl).unwrap();
    let opts = EvolutionOptions {
        stop: StopCondition::AsFarAsPossible,
        pulsation: false,
        max_steps: 400,
        ..Default::default()
    };
    let mut evo = Evolution::new(zams.clone(), opts);
    let mut models = vec![zams];
    loop {
        if models.last().unwrap().globals.y_c < 0.1 {
            return (models, "Y_c below 0.1".into());
        }
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => models.push(evo.model.clone()),
            EvolutionStep::Finished(r) => return (models, r),
        }
    }
}

#[test]
fn helium_burning_makes_carbon_and_oxygen() {
    let (models, reason) = helium_star();
    let (first, last) = (&models[0], models.last().unwrap());
    let g = &last.globals;
    eprintln!(
        "1 M☉ helium star: {} models ({reason}), age {:.1} Myr; centre Y {:.3}, C {:.3}, O {:.3}, C/O {:.3}; \
         L {:.1} L☉; energy 3α {:.3}, ¹²C(α,γ) {:.3}",
        models.len(),
        g.age_yr / 1e6,
        g.y_c,
        g.c_c,
        g.o_c,
        g.c_c / g.o_c,
        g.luminosity,
        g.frac_3a,
        g.frac_ca
    );
    assert_eq!(first.globals.phase, Phase::Zams);
    assert!(models[1..].iter().all(|m| m.globals.phase == Phase::CoreHelium));
    assert!(g.y_c < 0.5, "helium burning did not get far: Y_c = {}", g.y_c);
    // Carbon first, oxygen later: the C/O ratio falls as helium is used up.
    assert!(g.c_c > 0.1 && g.o_c > 0.01, "C {} O {}", g.c_c, g.o_c);
    let mid = models.iter().find(|m| m.globals.y_c < 0.75).unwrap();
    assert!(mid.globals.c_c / mid.globals.o_c > g.c_c / g.o_c);
    assert!(g.frac_ca > 0.0 && g.frac_3a > g.frac_ca);
    // Z₀ is not changed by burning; mass fractions sum to 1 at every point.
    for m in &models {
        let p = &m.profile;
        for k in 0..p.len() {
            let (c, o) = p.c_o_at(k);
            let s = p.x[k] + p.y[k] + p.z[k] + c + o;
            assert!((s - 1.0).abs() < 1e-12, "sum {s} at k = {k}");
            assert!((p.z[k] - 0.02).abs() < 1e-12);
        }
    }
}

/// Fuel accounting per channel: carbon made × Q_He against ∫∫ε_3α dm dt, and
/// carbon burnt to oxygen × Q_Cα against ∫∫ε_Cα dm dt (to the time
/// discretisation: the burning step uses the rates of the previous model).
#[test]
fn fuel_accounting_per_channel() {
    let (models, _) = helium_star();
    let c = models.last().unwrap().globals.energy_check().unwrap();
    let [h, he, ca] = c.fuel_by_channel;
    eprintln!(
        "1 M☉ helium star: fuel × Q / ∫∫ε: H {h:?}, 3α {he:?}, ¹²C(α,γ) {ca:?}; all channels / ∫L dt {:.4}; \
         energy residual {:.2e}",
        c.fuel_ratio, c.residual
    );
    assert!(h.is_none(), "no hydrogen to burn");
    let (he, ca) = (he.unwrap(), ca.unwrap());
    assert!((he - 1.0).abs() < 0.02, "3α {he}");
    assert!((ca - 1.0).abs() < 0.05, "¹²C(α,γ) {ca}");
    assert!((c.fuel_ratio - 1.0).abs() < 0.05, "{}", c.fuel_ratio);
}
