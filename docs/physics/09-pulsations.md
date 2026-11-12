# 9. Pulsations

## Implemented

Source: `crates/stellar-core/src/pulsation.rs`, `structure/model.rs` (scaling relations).

### Radial fundamental mode

Linear adiabatic wave equation in mass coordinate, for the relative displacement $`\xi`$ and the
relative Lagrangian pressure perturbation $`p`$ of a mode with angular frequency $`\omega`$:

```math
\xi = \frac{\delta r}{r}, \qquad p = \frac{\delta P}{P}
```

```math
\frac{d\xi}{dm} = -\frac{3\xi + p/\Gamma_1}{4\pi r^3\rho}, \qquad \frac{dp}{dm} = \frac{G m}{4\pi r^4 P}\left[\left(4 + \frac{\omega^2 r^3}{G m}\right)\xi + p\right]
```

*Boundary conditions and eigenvalue*

```math
3\xi + \frac{p}{\Gamma_1} = 0 \;\; (\text{centre}), \qquad p + (4 + \Omega^2)\,\xi = 0,\;\; \xi = 1 \;\; (\text{surface}), \qquad \Omega^2 = \frac{\omega^2 R^3}{G M}
```

$`\Omega^2`$ is solved as an extra unknown with the Henyey solver. The fundamental mode is the
solution without nodes in $`\xi`$.

*What the surface condition assumes.* It is the regularity condition of the wave equation at a
surface where the pressure goes to zero ($`dp/dm`$ stays finite only if the bracket vanishes),
applied at the photosphere with $`m = M`$, $`r = R`$. It treats the surface as a perfect reflector
and ignores the atmosphere above it. That is justified when the mode frequency is well below the
acoustic cut-off frequency of the atmosphere, $`\omega_{\mathrm{ac}} \approx c_s/(2H_P)`$, which
is usually, but not always, true of the fundamental mode. It is not assumed: the solver reports
$`\omega/\omega_{\mathrm{ac}}`$ at the photosphere and flags the result when the ratio is not small
(above about 0.5). The condition is not adequate for high overtones.

*What the calculation does not say.* The equations are adiabatic: they give a period for every
model, whether or not the star pulsates. Adiabaticity itself fails in the outer layers, where the
thermal time is shorter than the period; the period is little affected, the growth rate entirely.

*Period and pulsation constant*

```math
\Pi = \frac{2\pi}{\omega}, \qquad Q = \Pi\,\sqrt{\frac{\bar\rho}{\bar\rho_\odot}}
```

*Validation*: for a uniform-density sphere with constant $`\Gamma_1`$,

```math
\Omega^2 = 3\Gamma_1 - 4, \qquad \xi = \text{const}
```

The interval midpoints of this solver use arithmetic means of $`r`$, $`P`$, $`\rho`$.

### Instability strip

Whether a mode is driven is not computed. An empirical band in the HR diagram decides it, for
$`\log L/L_\odot \ge 0.3`$:

```math
\log T_{\mathrm{blue}} = 3.94 - 0.037\,(\log L - 0.5), \qquad \log T_{\mathrm{red}} = 3.80 - 0.023\,(\log L - 0.5)
```

Inside the band a star is labelled δ Scuti-type for $`\log L < 1.7`$ and Cepheid-type above.

The label is a lookup in observed boundaries. It is independent of the computed mode and is not a
stability result: nothing in the adiabatic calculation can place a model inside or outside the
strip. Since the model's $`T_{\mathrm{eff}}`$ is itself off for cool stars
([Validity](11-validity-and-roadmap.md)), the label inherits that error.

### Solar-like oscillations: scaling relations

```math
\nu_{\max} = 3090\,\mu\mathrm{Hz}\;\frac{g}{g_\odot}\left(\frac{T_{\mathrm{eff}}}{T_{\mathrm{eff},\odot}}\right)^{-1/2}, \qquad \Delta\nu = 135.1\,\mu\mathrm{Hz}\,\sqrt{\frac{\bar\rho}{\bar\rho_\odot}}
```

These use only $`M`$, $`R`$ and $`T_{\mathrm{eff}}`$, not the interior.

### Characteristic frequencies (E4 stage 1, plan 2.2)

Source: `crates/stellar-core/src/seismic.rs`; `Profile::cs`, `Profile::n2`; `Globals::delta_nu_asym_uhz`,
`Globals::g_cavities`. Diagnostics only: nothing feeds back into the structure.

*Sound speed, Lamb frequency and buoyancy (Brunt–Väisälä) frequency*

```math
c_s^2 = \frac{\Gamma_1 P}{\rho}, \qquad S_l^2 = \frac{l(l+1)\,c_s^2}{r^2}, \qquad N^2 = \frac{g^2 \rho}{P}\left[\delta\left(\nabla_{\mathrm{ad}} - \nabla\right) + \varphi\,\nabla_\mu\right]
```

$`g > 0`$ is the magnitude of the gravitational acceleration.

*Relation to the convective flag.* In a chemically homogeneous model $`N^2 < 0`$ where the
Schwarzschild criterion gives convection. The two do not coincide in general: the flag ignores
composition, while $`N^2`$ contains it, so a stabilising composition gradient gives $`N^2 > 0`$
inside a Schwarzschild-unstable (semiconvective) layer. $`N^2 < 0`$ is Ledoux instability. In
efficient convection $`\nabla \to \nabla_{\mathrm{ad}}`$ and $`N^2`$ is zero to within numerical
noise, of either sign.

The term $`\varphi\nabla_\mu`$ is the ideal-mixture form and is exact for the implemented EOS. For
the planned EOS it is replaced by the general composition term $`B`$ of
[Convection and mixing](05-convection-and-mixing.md):

```math
N^2 = \frac{g^2 \rho}{P}\left[\delta\left(\nabla_{\mathrm{ad}} - \nabla\right) + B\right]
```

This is the production expression: it does not difference nearly equal numbers. $`\nabla_\mu`$ is
taken by centred differences of $`\ln\mu`$ and $`\ln P`$ along the mesh (one-sided at the ends);
$`\delta`$, $`\nabla_{\mathrm{ad}}`$ and the actual $`\nabla`$ come from the EOS and the MLT of the
model. The definition
itself serves as a cross-check on the same profile (it is noisy in the deep interior, so it is a
test, not the output):

```math
N^2 = g\left(\frac{1}{\Gamma_1}\frac{d\ln P}{dr} - \frac{d\ln\rho}{dr}\right)
```

* p modes (pressure is the restoring force) propagate where $`\omega > S_l`$ and $`\omega > N`$.
* g modes (buoyancy) propagate where $`\omega < S_l`$ and $`\omega < N`$, so only in radiative
  layers.

### Asymptotic spacings (E4 stage 1, plan 2.2)

*Large frequency separation of p modes, and period spacing of g modes*

```math
\Delta\nu = \left[2\int_0^R \frac{dr}{c_s}\right]^{-1}, \qquad \Delta\Pi_l = \frac{2\pi^2}{\sqrt{l(l+1)}}\left[\int_{r_1}^{r_2} N\,\frac{dr}{r}\right]^{-1}
```

Both come from integrals over the existing profiles. The first replaces the $`\Delta\nu`$ scaling
relation by a value that depends on the interior.

*Implementation.* The $`\Delta\nu`$ integral uses the trapezoid rule in $`r`$, with
$`r_0/c_{s,0}`$ for the part inside the innermost point. A cavity is a run of consecutive points with
$`A = N^2 P/(g^2\rho) > 10^{-4}`$; the tolerance keeps the noise of nearly adiabatic layers (of either
sign, $`|A| \lesssim 10^{-6}`$) from opening spurious cavities, while radiative zones have
$`A \sim 0.1`$. Each cavity is reported with its edges in $`r/R`$ and $`m/M`$, $`\int N\,dr/r`$, its
own $`\Delta\Pi_1`$, a location (core, between convective regions, reaching the surface) and a flag
for thin layers ($`\ln(r_2/r_1) < 0.1`$ away from the centre, for example the radiative layer above a
convective envelope). The single number shown in the results panel is the cavity with the largest
$`\int N\,dr/r`$ that is not thin.

*Limits of the g-mode integral.* $`r_1`$ to $`r_2`$ is one connected g-mode cavity, a radiative
region with $`N^2 > 0`$. Disconnected regions are **not** summed: modes trapped in different
cavities form separate sequences. The output is therefore one value per radiative region, labelled
by where it lies. For a star with a convective core and a thin convective envelope the relevant
cavity is the radiative zone between them; for a red giant it is the radiative core. A narrow
$`N^2 > 0`$ layer (for example between two ionisation zones) is listed but is not a g-mode cavity
in the asymptotic sense unless many wavelengths fit in it.

A connected $`N^2 > 0`$ region is a *candidate* cavity. Where a g mode of a given frequency and
degree actually propagates is set by $`\omega < N`$ **and** $`\omega < S_l`$, so the true limits
$`r_1`$, $`r_2`$ depend on $`\omega`$ and $`l`$; the tabulated value uses the whole region and is
the low-frequency limit. Cavities separated by a thin evanescent layer are coupled, and the
observed spectrum is then not two independent sequences (mixed modes in subgiants and red giants
are the best-known case). The per-cavity values are diagnostics, not predicted spectra.

Both formulas are asymptotic (high radial order). The integrand $`1/c_s`$ is dominated by the outer
layers, so $`\Delta\nu`$ carries the model's surface errors.

### How good this is **[measured]**

* **N² against the convective flag** (ZAMS, 0.5–5 M☉): $`N^2 \le 0`$ at every convective point
  ($`A \le -10^{-11}`$ at most) and $`A > 10^{-3}`$ at every interior radiative point.
* **N² against its definition** (3 M☉ at 200 Myr, $`X_c \approx 0.3`$, away from composition jumps):
  median relative difference $`2.4\times10^{-5}`$; worst 2.2% where $`A > 10^{-2}`$; worst
  $`|\Delta A| = 7\times10^{-4}`$ in the nearly adiabatic radiative layers just below the photosphere.
* **Schwarzschild and Ledoux** differ in evolved models: in the same model one layer at the edge of the
  shrinking core is Schwarzschild-unstable but has $`N^2 > 0`$ (149 against 148 points).
* **Δν.** ZAMS Sun: 115.96 μHz asymptotic against 115.01 μHz scaling (+0.8%). The ratio falls with
  mass, to 0.91 for 3 and 5 M☉ ZAMS models.
* **ΔΠ₁.** ZAMS Sun: one cavity, the radiative core to 0.967 R (the envelope base), ΔΠ₁ = 3206 s,
  plus a thin surface layer. 1.5 M☉: the zone between core and envelope, 3512 s. 3 M☉: the radiative
  envelope, 5646 s. All tests in `tests/seismic.rs`.

## Planned

### Geometric means (C10, plan 1.3)

For consistency with the structure solver, at interval midpoints:

```math
r_{k+1/2} = \sqrt{r_k\,r_{k+1}}, \qquad P_{k+1/2} = \sqrt{P_k\,P_{k+1}}, \qquad \rho_{k+1/2} = \sqrt{\rho_k\,\rho_{k+1}}
```

### Non-radial modes, Cowling approximation (E4 stage 2, plan 9.1)

Neglecting the perturbation of the gravitational potential, for the radial displacement $`\xi_r`$
and the Eulerian pressure perturbation $`p'`$ of a mode of degree $`l`$:

```math
\frac{d\xi_r}{dr} = -\left(\frac{2}{r} - \frac{g}{c_s^{2}}\right)\xi_r + \frac{1}{\rho\,c_s^{2}}\left(\frac{S_l^{2}}{\omega^{2}} - 1\right)p'
```

```math
\frac{dp'}{dr} = \rho\left(\omega^{2} - N^{2}\right)\xi_r - \frac{g}{c_s^{2}}\,p'
```

*Source.* Compared (task 0.3, 2026-11-11) with Christensen-Dalsgaard, *Lecture Notes on Stellar
Oscillations*, 5th ed. (2014), §5.2.2, eqs. 5.12–5.14, p. 75, which give the same pair with the
coefficient $`\Gamma_1^{-1} H_p^{-1}`$, $`H_p^{-1} = -d\ln p/dr`$, in place of $`g/c_s^2`$. The two are
equal by hydrostatic equilibrium and $`c_s^2 = \Gamma_1 P/\rho`$. The notes obtain them from the full
equations (eqs. 4.61–4.62, pp. 65–66) by dropping $`\Phi'`$; their $`N^2`$ (eq. 4.63) and $`S_l^2`$
(eq. 4.60) are the definitions used above. *Decided 2026-11-12:* these lecture notes are the reference for
the Cowling equations and their boundary conditions. Unno et al. (1989), which the notes cite for
the boundary conditions (§18.1), could not be obtained and is not needed for this task.

*Boundary conditions*

```math
\frac{p'}{\rho} = \frac{\omega^{2}\,r\,\xi_r}{l} \;\; (\text{centre}), \qquad p' = \rho\,g\,\xi_r \;\; (\text{surface})
```

Both follow from the same notes: the centre from $`\xi_r \simeq l\,\xi_h`$ (eq. 4.65, p. 67) with
$`\xi_h = p'/(\rho r \omega^2)`$ (eq. 4.39, p. 62, with $`\Phi' = 0`$); the surface from
$`\delta p = p' + \xi_r\,dp/dr = 0`$ (eq. 4.68, p. 67) with $`dp/dr = -\rho g`$.

The central condition holds for $`l \ge 1`$ only (it divides by $`l`$). Radial modes ($`l = 0`$)
have their own regularity condition and are the subject of the implemented solver; a comparison of
the two at $`l = 0`$ must use that condition.

Two unknowns plus the eigenvalue: the same size as the radial problem. Accurate for high-order or
high-$`l`$ modes, where the notes cite numerical confirmation (§5.2.1, p. 75); for low-order,
low-$`l`$ modes it can fail badly (the notes single out the $`l = 1`$ f mode, "drastically
different" from the full equations). The earlier statement "several percent off for the lowest
ones" was not from a source and is withdrawn; the size of the error is to be measured against
plan 9.2.

### Non-radial modes, full problem (E4 stage 3, plan 9.2)

Four unknowns plus the eigenvalue, including the perturbed gravitational potential. *Decided
2026-11-12:* take the equations and boundary conditions from the same, available lecture notes
(Christensen-Dalsgaard 2014, eqs. 4.61–4.62, with the central and surface conditions around
eqs. 4.65–4.68, p. 67; the dimensionless form $`y_1 \dots y_4`$ is set up from them when the task
starts) instead of Unno et al. (1989), which was not available. Check the result against the Cowling limit above and against
the radial solver at $`l = 0`$.

### Finding the modes

For a trial $`\omega`$, impose the central condition and the normalisation $`\xi_r(R) = 1`$, solve
the linear system once, and evaluate the remaining surface condition:

```math
\mathcal{D}(\omega) = \left[p' - \rho\,g\,\xi_r\right]_{r = R}
```

A sign change of $`\mathcal{D}`$ brackets a mode only if $`\mathcal{D}`$ is continuous there. With
the normalisation $`\xi_r(R) = 1`$, $`\mathcal{D}`$ has a pole wherever the solution has a node at
the surface, and it changes sign across the pole as well. Two safeguards, of which the second is
preferred:

* accept a bracket only if $`|\mathcal{D}|`$ decreases towards the root on refinement (it grows
  towards a pole);
* use a two-sided matching discriminant: integrate from the centre and from the surface to an
  interior fitting point $`r_f`$, and take the determinant of the two solutions there (subscripts
  $`c`$ and $`s`$ for the integrations from centre and surface):

```math
\mathcal{D}(\omega) = \det\begin{pmatrix} \xi_{r,c}(r_f) & \xi_{r,s}(r_f) \\ p'_c(r_f) & p'_s(r_f) \end{pmatrix}
```

  The starting vectors at the two boundaries, and any rescaling applied during the integration,
  must be continuous functions of $`\omega`$, finite, non-zero, and rescaled only by positive
  factors (a sign flip introduces a discontinuity). Then $`\mathcal{D}`$ is continuous for
  $`\omega > 0`$ and has no artificial poles. No stronger claim is made: the equations contain
  $`1/\omega^2`$, so the search interval excludes $`\omega = 0`$; continuity is checked numerically
  when a bracket is refined; and a scan for sign changes can still miss modes (two between
  neighbouring scan points), so the node count of consecutive modes is used to detect gaps. Modes are labelled by counting nodes
(Eckart–Scuflaire–Osaki scheme).

### Not planned yet: excitation

A non-adiabatic analysis would say which modes grow. It also needs an opacity with the iron peak
near $`2\times10^5`$ K, which drives β Cephei and slowly pulsating B stars. Until then, driving stays
empirical.

---
Previous: [Rotation](08-rotation.md) · [Contents](README.md) · Next: [Numerical method and checks](10-numerics-and-checks.md)
