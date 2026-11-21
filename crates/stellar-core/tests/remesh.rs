//! Remeshing during evolution (plan 5.R): a converged model is unchanged by
//! remeshing to within the discretisation error, every species is conserved by
//! the remap, and the luminosity of a hydrogen-burning shell converges as the
//! refinement criteria are tightened.

use stellar_core::solver::NewtonOptions;
use stellar_core::structure::equations::NVAR;
use stellar_core::structure::mesh::remap;
use stellar_core::structure::remesh::{RemeshOptions, refine_mesh, remap_composition, species_masses};
use stellar_core::structure::solve::solve_structure;
use stellar_core::{
    Evolution, EvolutionOptions, EvolutionStep, NoControl, Phase, Physics, StarModel, StarParams, StopCondition,
    solve_model,
};

/// Criteria tight enough to refine a main-sequence model nearly everywhere.
const TIGHT: RemeshOptions =
    RemeshOptions { dlnp: 0.02, dlnt: 0.01, dx: 0.01, dl: 0.005, max_depth: 4, max_extra: 1.0 };

fn compositions(m: &StarModel) -> Vec<stellar_core::Composition> {
    (0..m.profile.len()).map(|k| m.profile.composition_at(k, m.params.f_cno)).collect()
}

/// A converged ZAMS model re-solved on a remeshed grid: L and R change by less
/// than 10⁻⁴, the size of the mesh error at 800 points (V.2: 5 × 10⁻⁵).
#[test]
fn remeshing_a_converged_model() {
    let physics = Physics::analytic();
    for mass in [1.0, 3.0] {
        let params = StarParams { mass, ..StarParams::sun() };
        let m = solve_model(&params, &physics, None, &mut NoControl).unwrap();
        let q_old = &m.profile.q;
        // Loose criteria leave the mesh alone.
        let loose = RemeshOptions { dlnp: 10.0, dlnt: 10.0, dx: 10.0, dl: 10.0, ..TIGHT };
        assert_eq!(refine_mesh(&m, q_old, &loose), *q_old);
        let q_new = refine_mesh(&m, q_old, &TIGHT);
        assert!(q_new.len() > q_old.len() + 200, "{} → {} points", q_old.len(), q_new.len());
        assert!(q_new.windows(2).all(|w| w[1] > w[0]));
        assert!(q_old.iter().all(|q| q_new.contains(q)), "base points kept");
        let mass_g = m.mass_g();
        let m_old: Vec<f64> = q_old.iter().map(|q| q * mass_g).collect();
        let m_new: Vec<f64> = q_new.iter().map(|q| q * mass_g).collect();
        let comps = remap_composition(&m_old, &compositions(&m), &m_new, params.f_cno);
        let y0 = remap(q_old, &m.solution, NVAR, &q_new);
        let opts = NewtonOptions { max_iter: 40, tol: 1e-8 };
        let r =
            solve_structure(&params, mass, &physics, &q_new, comps, None, y0, m.l_scale, opts, &mut NoControl).unwrap();
        let dl = r.globals.luminosity / m.globals.luminosity - 1.0;
        let dr = r.globals.radius / m.globals.radius - 1.0;
        eprintln!("{mass} M☉ ZAMS, {} → {} points: ΔL/L = {dl:.2e}, ΔR/R = {dr:.2e}", q_old.len(), q_new.len());
        assert!(dl.abs() < 1e-4 && dr.abs() < 1e-4);
    }
}

/// Every species is conserved to rounding by the remap, on an evolved 3 M☉
/// model with the composition gradient left by its receding core.
#[test]
fn remap_conserves_species_of_an_evolved_model() {
    let physics = Physics::analytic();
    let zams = solve_model(&StarParams { mass: 3.0, ..StarParams::sun() }, &physics, None, &mut NoControl).unwrap();
    let opts = EvolutionOptions { stop: StopCondition::Age(1.3e8), pulsation: false, ..Default::default() };
    let mut evo = Evolution::new(zams, opts);
    while let EvolutionStep::Advanced { .. } = evo.step(&physics, &mut NoControl).unwrap() {}
    let m = evo.model;
    let mass_g = m.mass_g();
    let q_new = refine_mesh(&m, &m.profile.q, &TIGHT);
    let m_old: Vec<f64> = m.profile.q.iter().map(|q| q * mass_g).collect();
    let m_new: Vec<f64> = q_new.iter().map(|q| q * mass_g).collect();
    let old = compositions(&m);
    let new = remap_composition(&m_old, &old, &m_new, m.params.f_cno);
    let (a, b) = (species_masses(&m_old, &old), species_masses(&m_new, &new));
    for i in 0..5 {
        assert!((a[i] - b[i]).abs() <= 1e-13 * mass_g, "species {i}: {:e} vs {:e}", a[i], b[i]);
    }
    assert!(new.iter().all(|c| (c.sum() - 1.0).abs() < 1e-13));
}

/// Total hydrogen-burning luminosity [erg s⁻¹].
fn l_hydrogen(m: &StarModel) -> f64 {
    m.globals.e_nuc * (m.globals.frac_pp + m.globals.frac_cno)
}

/// Shell burning: from one model of a 3 M☉ star past the main sequence (400
/// base points), the same time step on the base mesh and with remeshing of
/// increasing fineness. The shell luminosity converges: the default criteria
/// are closer to the finest than the base mesh is, and within 1% of it.
#[test]
fn shell_luminosity_converges_with_refinement() {
    let physics = Physics::analytic();
    let params = StarParams { mass: 3.0, n_points: 400, ..StarParams::sun() };
    let zams = solve_model(&params, &physics, None, &mut NoControl).unwrap();
    let opts = EvolutionOptions { stop: StopCondition::AsFarAsPossible, pulsation: false, ..Default::default() };
    let mut evo = Evolution::new(zams, opts);
    let r_zams = evo.r_zams();
    // Into the Hertzsprung gap: a helium core of 5% of the mass, burning in a shell.
    while !(evo.model.globals.phase == Phase::ShellHydrogen && evo.model.globals.he_core_q > 0.05) {
        match evo.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { .. } => {}
            EvolutionStep::Finished(r) => panic!("finished before the shell phase: {r}"),
        }
    }
    let start = evo.model.clone();
    let dt = 0.25 * evo.dt_yr;
    let fine = RemeshOptions { dlnp: 0.05, dlnt: 0.025, dx: 0.0125, dl: 0.005, max_depth: 10, max_extra: 4.0 };
    let finer = RemeshOptions { dlnp: 0.025, dlnt: 0.0125, dx: 0.00625, dl: 0.0025, max_depth: 12, max_extra: 8.0 };
    let step = |remesh: Option<RemeshOptions>| {
        let mut e = Evolution::resume(start.clone(), dt, r_zams, EvolutionOptions { remesh, ..opts });
        e.dt_yr = dt;
        match e.step(&physics, &mut NoControl).unwrap() {
            EvolutionStep::Advanced { dt_yr, .. } => assert!((dt_yr / dt - 1.0).abs() < 1e-9, "step split"),
            EvolutionStep::Finished(r) => panic!("{r}"),
        }
        (l_hydrogen(&e.model), e.model.profile.len())
    };
    let (l_base, n_base) = step(None);
    let (l_def, n_def) = step(Some(RemeshOptions::default()));
    let (l_fine, n_fine) = step(Some(fine));
    let (l_finer, n_finer) = step(Some(finer));
    let rel = |l: f64| l / l_finer - 1.0;
    eprintln!(
        "3 M☉, q_He = {:.3}: L_H relative to the finest ({n_finer} points): base mesh ({n_base}) {:+.2e}, default \
         ({n_def}) {:+.2e}, fine ({n_fine}) {:+.2e}",
        start.globals.he_core_q,
        rel(l_base),
        rel(l_def),
        rel(l_fine)
    );
    assert!(n_base < n_def && n_def <= n_fine && n_fine <= n_finer);
    assert!(rel(l_def).abs() < rel(l_base).abs());
    assert!(rel(l_fine).abs() <= rel(l_def).abs() + 1e-4);
    assert!(rel(l_def).abs() < 0.01);
}
