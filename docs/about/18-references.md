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
* **de Jager, Nieuwenhuijzen & van der Hucht (1988)**, A&AS 72, 259: the hot-star mass-loss rate
  (linear form, §4 remark 3).
* **Vink, de Koter & Lamers (2001)**, A&A 369, 574: the metallicity factor $`(Z/0.019)^{0.69}`$ of the
  hot-star wind.
* **Reimers (1975)**, Mém. Soc. R. Sci. Liège, 6e sér., 8, 369, as quoted by **Reimers (1977)**, A&A
  61, 217, and **McDonald & Zijlstra (2015)**, MNRAS 448, 502: the cool-star mass-loss rate and its
  default $`\eta = 0.477`$.
* **von Zeipel (1924)**, MNRAS 84, 665, and **Lucy (1967)**, Z. Astrophys. 65, 89: gravity darkening
  of a rotating star (display only).
* **IAU 2015 Resolution B3:** nominal solar radius, luminosity and effective temperature.
* **Cox (1980)**, *Theory of Stellar Pulsation*, Princeton University Press: the linear adiabatic wave
  equation and its boundary conditions.
* **Pekeris (1938)**, ApJ 88, 189: the radial modes of the uniform-density sphere, used to test the
  pulsation solver.
* **Unno, Osaki, Ando, Saio & Shibahashi (1989)**, *Nonradial Oscillations of Stars*, 2nd ed.,
  University of Tokyo Press.
* **Brown, Gilliland, Noyes & Ramsey (1991)**, ApJ 368, 599, and **Kjeldsen & Bedding (1995)**, A&A 293,
  87: the $`\nu_{\max}`$ and $`\Delta\nu`$ scaling relations.

## Cool-star physics: data and methods (optional physics)

* **Iglesias & Rogers (1996)**, ApJ 464, 943: the OPAL opacity tables (file GN93hz, embedded unchanged).
* **Ferguson, Alexander, Allard, Barman, Bodnarik, Hauschildt, Heffner-Wong & Tamanai (2005)**, ApJ 623, 585: the low-temperature opacity tables with molecules and grains (set f05.g93, embedded unchanged).
* **Grevesse & Noels (1993)**, in Origin and Evolution of the Elements, eds. Prantzos, Vangioni-Flam & Cassé, Cambridge University Press, p. 15: the heavy-element mixture of both opacity tables.
* **Potekhin, Pons & Page (2015)**, Space Sci. Rev. 191, 239; Cassisi, Potekhin, Pietrinferni, Catelan & Salaris (2007), ApJ 661, 1094; Cassisi, Potekhin, Salaris & Pietrinferni (2021), A&A 654, A149: the electron thermal conductivity table of A. Y. Potekhin (condtab21wd.dat, Ioffe Institute, embedded unchanged).
* **Krishna Swamy (1966)**, ApJ 145, 174: the solar T(τ) relation (eq. 33).
* **Henyey, Vardya & Bodenheimer (1965)**, ApJ 142, 841: convection with optically thin elements.
* **Trampedach, Stein, Christensen-Dalsgaard, Nordlund & Asplund (2014)**, MNRAS 442, 805: the Hopf-corrected radiative gradient in the atmosphere and the range of the Krishna Swamy relation.
* **Pols, Tout, Eggleton & Han (1995)**, MNRAS 274, 964: ionisation and the pressure-ionisation free-energy term.
* **Kovetz, Yaron & Prialnik (2009)**, MNRAS 395, 1857: a reproduction of the pressure-ionisation term.
* **Paczyński (1983)**, ApJ 267, 315: the interpolation formula for the degenerate electron pressure (comparison EOS).
* **Timmes & Arnett (1999)**, ApJS 125, 277: the Fermi–Dirac electron gas of stellar EOSs.
* **Kramida, Ralchenko, Reader & NIST ASD Team**, NIST Atomic Spectra Database: ionisation energies of H and He (retrieved 2026-11-11).
* **Huber & Herzberg (1979)**, Molecular Spectra and Molecular Structure IV: Constants of Diatomic Molecules, Van Nostrand Reinhold: the H₂ constants.
* **Golub & Welsch (1969)**, Math. Comp. 23, 221: Gauss quadrature nodes and weights.
* **Vinyoles, Serenelli, Villante, Basu, Bergström, Gonzalez-Garcia, Maltoni, Peña-Garay & Song (2017)**, ApJ 835, 202: the B16 standard solar models, used for comparison with the provisional solar calibration.
* **Tout, Pols, Eggleton & Han (1996)**, MNRAS 281, 257: ZAMS luminosity and radius fits of the Eggleton-code models, used for validation (plan 7.9).
* **Hurley, Pols & Tout (2000)**, MNRAS 315, 543: main-sequence lifetime fits of the same models (coefficients as in the SSE/BSE source distributed with COSMIC).

The three tables are embedded unchanged; they remain the work of their authors and are not
covered by the licence of the code. Sources, checksums and licence status:
[`crates/stellar-core/data/README.md`](../../crates/stellar-core/data/README.md).

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
