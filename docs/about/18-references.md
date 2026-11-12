# 18. References & software

## Physics and methods

* **Kippenhahn, Weigert & Weiss (2012)**, *Stellar Structure and Evolution*, 2nd ed., Springer: the
  structure equations, mixing-length formulation, Henyey method and polytropes.
* **Kippenhahn & Weigert (1990)**, *Stellar Structure and Evolution*, 1st ed., Springer: the analytic
  triple-α rate formula (and the pp and CNO formulas used before the Solar Fusion III rates).
* **Carroll & Ostlie (2017)**, *An Introduction to Modern Astrophysics*, 2nd ed., Cambridge University
  Press: the analytic opacity formulas and the StatStar model.
* **Hansen, Kawaler & Trimble (2004)**, *Stellar Interiors*, 2nd ed., Springer.
* **Henyey, Forbes & Gould (1964)**, ApJ 139, 306: the relaxation method for stellar structure.
* **Böhm-Vitense (1958)**, Zeitschrift für Astrophysik 46, 108: mixing-length theory.
* **Clayton (1968)**, *Principles of Stellar Evolution and Nucleosynthesis*, McGraw-Hill: thermonuclear
  rates and weak screening.
* **Asplund, Grevesse, Sauval & Scott (2009)**, ARA&A 47, 481: the AGSS09 solar composition preset.
* **Acharya et al. (2025)**, Rev. Mod. Phys. 97, 035002 ("Solar Fusion III"): the astrophysical
  factors $`S(0)`$, $`S'(0)`$, $`S''(0)`$ behind the pp and CNO rates and the pp II branching factor
  $`\psi`$.
* **Paxton et al. (2011)**, ApJS 192, 3 (MESA), §5.1: the bounded pressure scale height in the mixing
  length.
* **Choi et al. (2016)**, ApJ 823, 102 (MIST), §3.1: $`Y_p`$ and $`\Delta Y/\Delta Z`$ of the "helium follows
  metallicity" mode.
* **IAU 2015 Resolution B3:** nominal solar radius, luminosity and effective temperature.
* **Cox (1980)**, *Theory of Stellar Pulsation*, Princeton University Press: the linear adiabatic wave
  equation and its boundary conditions.
* **Pekeris (1938)**, ApJ 88, 189: the radial modes of the uniform-density sphere, used to test the
  pulsation solver.
* **Unno, Osaki, Ando, Saio & Shibahashi (1989)**, *Nonradial Oscillations of Stars*, 2nd ed.,
  University of Tokyo Press.
* **Brown, Gilliland, Noyes & Ramsey (1991)**, ApJ 368, 599, and **Kjeldsen & Bedding (1995)**, A&A 293,
  87: the $`\nu_{\max}`$ and $`\Delta\nu`$ scaling relations.

## Software

* **[Rust](https://www.rust-lang.org):** the programming language. The physics and solver are written
  from scratch in the `stellar-core` crate.
* **[egui / eframe](https://github.com/emilk/egui) 0.37** and
  **[egui_plot](https://github.com/emilk/egui_plot) 0.38:** user interface, plots and window (native on
  Linux, macOS and Windows).
* **rfd:** native file dialogs.
* **png:** screenshot encoding.
* **serde / serde_json:** JSON export.
* **clap:** CLI arguments.
* **MathJax:** typesets the equations. The app embeds them pre-rendered as SVG; GitHub renders these
  pages with MathJax too.

> Black-body surface colours in the 3D view come from an approximate temperature-to-sRGB table
> (illustrative only).

---
← Previous: [Command-line tool](17-cli.md) · [Contents](README.md)
