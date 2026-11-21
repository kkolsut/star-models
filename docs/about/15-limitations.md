# 15. Accuracy & limitations

The analytic microphysics gets the structure **qualitatively** right:

* a radiative core and convective envelope in Sun-like stars;
* convective cores above about 1.2 M☉, with CNO dominance above about 1.5 M☉;
* deep convective envelopes in low-mass stars;
* massive stars close to detailed models, because their opacity is mostly electron scattering.

It is **not quantitatively solar**. With $`X = 0.70`$ and $`Z = 0.02`$:

| 1 M☉ ZAMS | $`L`$ | $`R`$ | Convective envelope base |
|---|---|---|---|
| This model (analytic physics) | ≈ 0.42 L☉ | ≈ 1.12 R☉ | $`0.97\,R`$ (thin) |
| This model, *Physics › All upgrades* | ≈ 0.65 L☉ | ≈ 0.90 R☉ | $`0.73\,R`$ |
| Detailed ZAMS Sun | ≈ 0.7 L☉ | ≈ 0.89 R☉ | ≈ $`0.73\,R`$ |

This model also has $`T_{\mathrm{eff}} \approx 4380`$ K. The causes were measured by swapping one ingredient at a time (plan 7.7): the opacity tables alone raise L to 0.68 L☉ (the analytic central κ is 1.8 times the OPAL value), the EOS with ionisation alone moves the envelope base to 0.76 R, and the integrated atmosphere changes R by at most 1%. Only all of them together come close: with Physics › All upgrades the same 1 M☉ model has L = 0.648 L☉, R = 0.905 R☉, T_eff = 5444 K and its envelope base at 0.727 R, uncalibrated. That set is about 130 times slower and is not the default (see [`docs/physics/11`](../physics/11-validity-and-roadmap.md)).

*Provisional solar calibration* (Physics › Calibrated Sun): with all upgrades, $`Y_0 = 0.2834`$,
$`Z_0 = 0.01714`$ ($`Z/X = 0.0245`$, Grevesse & Noels 1993) and $`\alpha_{\mathrm{MLT}} = 2.020`$ give $`L_\odot`$
and $`R_\odot`$ at 4.57 Gyr, $`T_{\mathrm{eff}} = 5772`$ K and the envelope base at 0.733 R (helioseismic 0.713 R;
no element settling yet, which deepens it). Against the published Eggleton-code grid (Tout et al. 1996, Hurley et al. 2000; plan 7.9) the full set gives ZAMS luminosities within 7% and radii within 4% for 0.5–20 M☉, and main-sequence lifetimes 8–19% shorter (that grid includes core overshoot, which is off here by default).

With $`Z = 0.008`$ the same physics reproduces the classic StatStar Sun of Carroll & Ostlie
($`L \approx 0.83\,L_\odot`$, $`T_{\mathrm{eff}} \approx 5500`$ K). This confirms that the
implementation is consistent with the textbook formulas. This is a check of the code (verification),
not a validation against the real Sun.

The numerical checks (virial, energy, residuals, mesh refinement, time-step convergence) show that the
equations are solved accurately. They say nothing about how realistic the physics is.

## Evolution

* **Lifetimes:** main-sequence lifetimes follow the model luminosities.
  * The analytic Sun is underluminous (0.42 L☉ at the ZAMS), so it burns slowly. It leaves the main
    sequence after about 15.1 Gyr instead of about 10 Gyr.
  * Intermediate-mass and massive stars come out within about 25% of detailed models:

    | Mass | This model | Detailed models |
    |---|---|---|
    | 1.8 M☉ | 1.32 Gyr | roughly 1.2–1.5 Gyr |
    | 5 M☉ | 72 Myr | roughly 90–100 Myr |
* **After the main sequence** the results are qualitative.
  * Without electron degeneracy, hydrogen and helium ionisation or molecular opacities, low-mass stars
    do not climb the red giant branch.
  * Instead they cool towards 2600 K at roughly constant luminosity. The helium flash cannot occur, and
    runs stop there.
* **Burning shell:** the hydrogen-burning shell after the main sequence is thin in mass.
  * Global results are converged at 800 points: for the 5 M☉ crossing, $`L`$ is within 0.2% of the
    3000-point value.
  * More mesh points resolve the shell's profile better.
* **Core overshoot** is a step in the mixing only ($`\alpha_{\mathrm{ov}}`$, off by default); its value is
  a calibration, not a prediction.
* **Mass loss** (off by default) uses two fitted rates: de Jager et al. (1988) in its linear form,
  times a metallicity factor taken from Vink et al. (2001), above 8000 K, and Reimers below. The
  combination, the step at 8000 K and the efficiency factor are choices, not predictions; real
  hot-star winds depend on more than $`L`$, $`T_{\mathrm{eff}}`$ and $`Z`$. Reimers' $`\eta`$ is
  calibrated on low-mass red giants, which this model does not reach. No Wolf–Rayet,
  luminous-blue-variable or dust-driven winds.
* **Not included:** diffusive (exponential) overshoot, envelope overshoot, semiconvection, rotational
  mixing, neutrino losses.
* **Pulsation:** the period is computed, but driving is not. The instability strip is an empirical band,
  so stars are flagged by position, not by a stability analysis.

## Not modelled

* Partial ionisation of H and He (Saha equation), which affects the convective envelopes of cool stars.
* Electron degeneracy and Coulomb interactions, which matter below about 0.4 M☉.
* Tabulated opacities (OPAL/OP, low-temperature molecular tables) and detailed model atmospheres; a grey
  Eddington atmosphere is used instead.
* Neutrino losses, rotation in the structure, magnetic fields, element diffusion and non-radial
  oscillation modes. The 3D view can draw a rotating star (Roche shape, gravity darkening), but the
  interior is not recomputed.
* Very massive stars (above about 50 M☉) approach the Eddington limit, where this simple treatment
  becomes crude.
* General relativity. The structure equations are Newtonian; the relativistic corrections are
  $`GM/(Rc^2) \approx 2\times10^{-6}`$ for the Sun and $`1.2\times10^{-5}`$ for 60 M☉, below the
  discretisation error. The gravitational redshift of the surface light is reported as an output
  (636 m/s for 1 M☉ at 1 R☉).
* Helium ignition. Screening covers all regimes (Chugunov & DeWitt), but thermal-neutrino losses are
  missing, so runs with Fermi–Dirac electrons still stop when the core becomes degenerate; each model
  reports its largest screening exponent and warns above 0.3.

## Planned

* Tabulated opacities and an EOS with ionisation and degeneracy, plugged into the existing physics
  interfaces. This is the step that would make the solar model quantitative and allow the red giant
  branch.
* A web version (see [the web version notes](../web-version.md)).

> Use the app to understand how stars are built and how their structure depends on mass and
> composition, not for precision work.

---
← Previous: [How the views are made](14-views.md) · [Contents](README.md) · Next: [Constants & units](16-constants.md) →
