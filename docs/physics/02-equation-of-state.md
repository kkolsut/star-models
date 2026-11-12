# 2. Equation of state

The equation of state (EOS) gives the density and the thermodynamic derivatives as functions of
$`(P, T)`$ and the composition. In the code it is the `Eos` trait, which returns `EosState`
(`crates/stellar-core/src/physics/mod.rs`).

## Implemented: fully ionised ideal gas plus radiation

Source: `crates/stellar-core/src/physics/eos.rs`, `crates/stellar-core/src/composition.rs`.

```math
P = P_{\mathrm{gas}} + P_{\mathrm{rad}} = \frac{\rho\, k\, T}{\mu\, m_u} + \frac{a T^4}{3}
```

*Mean molecular weight (full ionisation) and mean molecular weight per electron*

```math
\frac{1}{\mu} = 2X + \frac{3}{4}\,Y + \frac{1}{2}\,Z, \qquad \mu_e = \frac{2}{1 + X}
```

*Gas-pressure fraction and thermal expansion coefficient*

```math
\beta = \frac{P_{\mathrm{gas}}}{P}, \qquad \delta = -\left(\frac{\partial \ln\rho}{\partial \ln T}\right)_{\!P} = \frac{4 - 3\beta}{\beta}
```

*Specific heat at constant pressure*

```math
c_P = \frac{k}{\mu\, m_u}\left[\frac{3}{2} + \frac{12\,(1-\beta)}{\beta} + \frac{(4 - 3\beta)^2}{\beta^2}\right]
```

*Adiabatic gradient and first adiabatic exponent*

```math
\nabla_{\mathrm{ad}} = \left(\frac{\partial \ln T}{\partial \ln P}\right)_{\!s} = \frac{2\,(4 - 3\beta)}{32 - 24\beta - 3\beta^2}, \qquad \Gamma_1 = \left(\frac{\partial\ln P}{\partial\ln\rho}\right)_{\!s} = \frac{32 - 24\beta - 3\beta^2}{24 - 21\beta}
```

Limits: $`\nabla_{\mathrm{ad}} = 0.4`$, $`\Gamma_1 = 5/3`$ for a pure gas; $`\nabla_{\mathrm{ad}} = 0.25`$,
$`\Gamma_1 = 4/3`$ for pure radiation.

Given $`P`$ and $`T`$, the density follows directly from the first equation. During Newton
iterations $`\beta`$ is limited to at least $`10^{-3}`$; a converged model that reaches this floor is
rejected.

### What this EOS leaves out

* Partial ionisation of H and He: $`\nabla_{\mathrm{ad}}`$ stays at 0.4 in the envelope, where it
  really drops to about 0.1.
* Electron degeneracy, which matters below about 0.4 M☉ and in helium cores.
* Molecular hydrogen, Coulomb forces and pressure ionisation.

## Planned

### Derived quantities for any EOS (plan 6.1)

Every new EOS provides $`\rho(P,T)`$ and the specific enthalpy $`h(P,T)`$. The rest follows by
numerical differentiation:

```math
\alpha = \left(\frac{\partial\ln\rho}{\partial\ln P}\right)_{T}, \qquad \delta = -\left(\frac{\partial\ln\rho}{\partial\ln T}\right)_{P}, \qquad c_P = \left(\frac{\partial h}{\partial T}\right)_{P}
```

```math
\nabla_{\mathrm{ad}} = \frac{P\,\delta}{T\,\rho\,c_P}, \qquad \Gamma_1 = \frac{1}{\alpha - \delta\,\nabla_{\mathrm{ad}}}
```

Applied to the implemented EOS (for which $`\alpha = 1/\beta`$), these reproduce the closed forms
above.

**Thermodynamic consistency.** $`\rho(P,T)`$ and $`h(P,T)`$ cannot be chosen independently. Both
must follow from one potential (a free energy), which is equivalent to the identity

```math
\left(\frac{\partial h}{\partial P}\right)_{T} = \frac{1 - \delta}{\rho}
```

If it fails, $`\nabla_{\mathrm{ad}}`$, $`c_P`$ and $`\Gamma_1`$ computed from the formulas above are
not mutually consistent, and $`\varepsilon_{\mathrm{grav}}`$ does not conserve energy. The derivatives are at fixed
composition $`X_i`$.

*The test.* Both sides vanish for an ideal gas without radiation ($`\delta = 1`$), so a relative
difference is undefined there. Every new EOS is scanned over its $`(P,T)`$ range with the
dimensionless residual

```math
\mathcal{R} = \rho\left(\frac{\partial h}{\partial P}\right)_{T, X_i} - (1 - \delta), \qquad |\mathcal{R}| \le \epsilon_{\mathrm{abs}} + \epsilon_{\mathrm{rel}}\,|1 - \delta|
```

and the largest $`|\mathcal{R}|`$ is reported. The implemented EOS gives $`\mathcal{R} = 0`$.
Together with it, two stability conditions are checked at every point:

```math
c_V > 0, \qquad \alpha = \left(\frac{\partial\ln\rho}{\partial\ln P}\right)_{T} > 0
```

What the test covers: thermodynamic integrability at fixed composition. It does not verify
chemical (ionisation) equilibrium, the composition derivatives used for buoyancy, or anything
about accuracy against nature.

### Electron degeneracy (D15b, plan 6.2): Fermi–Dirac electrons **[compared]**

**Decided 2026-11-11:** the degenerate EOS is a Fermi–Dirac electron gas with ideal ions and
radiation. It treats partial degeneracy and relativity consistently and needs no table. The
Paczyński interpolation further down is kept only as an optional comparison model.

```math
P = P_{\mathrm{ion}} + P_e + P_{\mathrm{rad}}, \qquad P_{\mathrm{ion}} = \frac{\rho\,k\,T}{\mu_i\,m_u}, \qquad \frac{1}{\mu_i} = X + \frac{Y}{4} + \frac{Z}{\langle A\rangle}, \qquad P_{\mathrm{rad}} = \frac{a T^4}{3}
```

*Electrons.* With the degeneracy parameter $`\eta = \mu_e^{\mathrm{chem}}/(kT)`$ (electron chemical
potential without rest mass; not the molecular weight $`\mu_e`$), the relativity parameter $`\vartheta = kT/(m_e c^2)`$ and the generalised
Fermi–Dirac integrals

```math
F_k(\eta, \vartheta) = \int_0^\infty \frac{x^{k}\,\left(1 + \tfrac{1}{2}\vartheta x\right)^{1/2}}{e^{\,x - \eta} + 1}\,dx
```

the number density, pressure and energy density of the electrons are

```math
n_e = \frac{8\pi\sqrt{2}}{h^{3}}\,m_e^{3}c^{3}\,\vartheta^{3/2}\left[F_{1/2} + \vartheta\,F_{3/2}\right]
```

```math
P_e = \frac{16\pi\sqrt{2}}{3\,h^{3}}\,m_e^{4}c^{5}\,\vartheta^{5/2}\left[F_{3/2} + \tfrac{1}{2}\vartheta\,F_{5/2}\right], \qquad \rho\,u_e = \frac{8\pi\sqrt{2}}{h^{3}}\,m_e^{4}c^{5}\,\vartheta^{5/2}\left[F_{3/2} + \vartheta\,F_{5/2}\right]
```

*Closure and enthalpy*

```math
n_e = \frac{\rho}{\mu_e\,m_u}, \qquad h = \frac{5}{2}\,\frac{k T}{\mu_i\,m_u} + u_e + \frac{P_e}{\rho} + \frac{4\,a T^{4}}{3\,\rho}
```

All three electron quantities derive from one grand potential, so the consistency identity above
holds exactly (up to quadrature error). Limits, to be used as tests:

| Limit | Result |
|---|---|
| $`\eta \to -\infty`$ (non-degenerate), $`\vartheta \ll 1`$ | $`P_e = n_e kT`$, $`u_e = \tfrac32 P_e/\rho`$: the implemented EOS |
| $`\eta \gg 1`$ (degenerate), $`x_F \ll 1`$ | $`P_e = K_1(\rho/\mu_e)^{5/3}`$, $`u_e = \tfrac32 P_e/\rho`$ |
| $`\eta \gg 1`$ (degenerate), $`x_F \gg 1`$ | $`P_e = K_2(\rho/\mu_e)^{4/3}`$, $`u_e = 3P_e/\rho`$ |

Large $`\eta`$ says only that the gas is degenerate. Whether it is relativistic is set by the Fermi
momentum:

```math
x_F = \frac{p_F}{m_e c} = \frac{h}{m_e c}\left(\frac{3\,n_e}{8\pi}\right)^{1/3}
```

$`K_1`$, $`K_2`$ are the closed forms given below for the Paczyński model, but with $`m_u`$ in place of
$`m_{\mathrm{H}}`$ (consistent with $`n_e = \rho/\mu_e m_u`$): $`1.0036\times10^{13}`$ and $`1.2435\times10^{15}`$ (cgs).

*Admissible input.* $`T > 0`$ and $`P > P_{\mathrm{rad}}(T) = aT^4/3`$. Outside it there is no
solution; the call returns an error, not a clamped state.

*Solving for the density.* The `Eos` trait takes $`(P,T)`$. The unknown is $`\eta`$: for a trial
$`\eta`$ the formulas give $`n_e`$, hence $`\rho`$, hence the total $`P`$, which is monotonic in
$`\eta`$. A bracketed one-dimensional root search converges to $`10^{-13}`$.

*Numerical requirements.* The result must be a smooth function of $`(P,T)`$ for the
finite-difference Jacobian, so there is no adaptive error control. Fixed nodes on a fixed interval
are not enough either: the integrand has a step of width $`\sim 1`$ at $`x = \eta`$, which moves
without bound as $`\eta`$ grows. The integration range is therefore split into segments whose
break points are smooth functions of $`\eta`$ (around $`x = \eta`$ for $`\eta > 0`$), each with a
fixed Gauss rule and a variable change that removes the $`x^{1/2}`$ behaviour at the origin.

The accuracy is claimed only on a stated domain, to be fixed and tested when implemented, for
example $`-100 \le \eta \le 10^4`$, $`10^{-6} \le \vartheta \le 1`$ at relative accuracy
$`10^{-12}`$ against a high-precision reference. Below the lower $`\eta`$ bound the classical
series is used; a call beyond the domain is an error. Derivatives come from the generic layer above, or
analytically from $`\partial F_k/\partial\eta`$ and $`\partial F_k/\partial\vartheta`$.

*Not included:* electron–positron pairs (irrelevant below about $`10^9`$ K), Coulomb corrections,
partial ionisation, pressure ionisation. The last three need a separate extension that is itself
derived from a free energy (next sections); adding them as independent corrections to $`P`$ or
$`h`$ would break the consistency this choice was made for.

### Comparison model: Paczyński interpolation

Optional, selectable for comparison; not used as the basis of evolution runs. Interpolation formula
of Paczyński (1983, ApJ 267, 315, eq. 4 on p. 317; compared with the original in task 0.3). The
paper writes it for a hydrogen-free He/C/Mg mixture, $`1/\mu_i = X_4/4 + X_{12}/12 + X_{24}/24`$
and $`\mu_e = 2`$, with a mass unit $`H`$ in $`P_i`$ and $`P_{\mathrm{nd}}`$ that it does not
define. The general $`1/\mu_i`$ below and the use of $`m_u`$ are generalisations of these pages:

```math
P = P_{\mathrm{ion}} + P_e + P_{\mathrm{rad}}, \qquad P_{\mathrm{ion}} = \frac{\rho\,k\,T}{\mu_i\,m_u}, \qquad \frac{1}{\mu_i} = X + \frac{Y}{4} + \frac{Z}{\langle A\rangle}
```

```math
P_e = \sqrt{P_{\mathrm{nd}}^{2} + P_{\mathrm{d}}^{2}}, \qquad P_{\mathrm{nd}} = \frac{\rho\,k\,T}{\mu_e\,m_u}, \qquad P_{\mathrm{d}} = \left(P_{\mathrm{nr}}^{-2} + P_{\mathrm{r}}^{-2}\right)^{-1/2}
```

*Non-relativistic and relativistic limits of complete degeneracy*

```math
P_{\mathrm{nr}} = K_1\left(\frac{\rho}{\mu_e}\right)^{5/3}, \qquad K_1 = \frac{(3\pi^2)^{2/3}\,\hbar^2}{5\,m_e\,m_{\mathrm{H}}^{5/3}} = 9.91\times10^{12}
```

```math
P_{\mathrm{r}} = K_2\left(\frac{\rho}{\mu_e}\right)^{4/3}, \qquad K_2 = \frac{(3\pi^2)^{1/3}\,\hbar c}{4\,m_{\mathrm{H}}^{4/3}} = 1.231\times10^{15}
```

The numbers are the published ones (cgs; the paper writes $`K_1 = 9.91\times10^{12}\,\mu_e^{-5/3}`$,
$`K_2 = 1.231\times10^{15}\,\mu_e^{-4/3}`$ with $`\rho`$ in place of $`\rho/\mu_e`$). The
paper gives no closed form; the closed forms reproduce the published numbers (9.906×10¹²,
1.2306×10¹⁵) with the mass of the hydrogen atom $`m_{\mathrm{H}}`$, and give 1.0036×10¹³ and
1.2435×10¹⁵ with $`m_u`$. *Corrected in task 0.3:* this page previously quoted the $`m_u`$ values
as Paczyński's. A reproduction of the paper uses the published values; with $`m_u`$ the
degenerate pressure is 1.3% (non-relativistic) and 1.0% (relativistic) higher at given
$`\rho/\mu_e`$.

*Internal energy of the electrons* (enters $`h`$)

The paper gives only the temperature derivative of the specific electron energy at fixed density,
together with $`f`$:

```math
\left(\frac{\partial u_e}{\partial\ln T}\right)_{\!\rho} = \frac{1}{f - 1}\,\frac{P_e}{\rho}\left(\frac{\partial\ln P_e}{\partial\ln T}\right)_{\!\rho}, \qquad f = \frac{d\ln P_{\mathrm{d}}}{d\ln\rho} = \frac{5}{3}\left(\frac{P_{\mathrm{d}}}{P_{\mathrm{nr}}}\right)^{2} + \frac{4}{3}\left(\frac{P_{\mathrm{d}}}{P_{\mathrm{r}}}\right)^{2}
```

and says that eq. 4 "combined with standard thermodynamic relations" gives all quantities needed.
This page uses

```math
u_e = \frac{P_e}{(f - 1)\,\rho}
```

Since $`f`$ depends on $`\rho`$ only, this has exactly the published temperature derivative. Its
temperature-independent part, $`P_{\mathrm{d}}/((f-1)\rho)`$ at $`T = 0`$, is a choice of these
pages and is not in the paper. The other energies of eq. 4 are $`u_{\mathrm{ion}} = \tfrac32
P_{\mathrm{ion}}/\rho`$ and $`u_{\mathrm{rad}} = 3P_{\mathrm{rad}}/\rho`$. The paper states
worst errors of 5% in total pressure or adiabatic gradient in small areas of the
$`(\rho, T)`$ plane, and below 1% in most cases, for $`T < 10^{9.5}`$ K.

The `Eos` trait takes $`(P,T)`$, so $`\rho`$ is found by a one-dimensional root search ($`P`$ is
monotonic in $`\rho`$).

**Status of this prescription.** It is an interpolation between limits, not an EOS derived from a
free energy. The pressure formula and the energy formula are each reasonable in the limits
(non-degenerate: $`f \to 5/3`$, $`u_e = \tfrac32 P_e/\rho`$; relativistic degenerate:
$`f \to 4/3`$, $`u_e = 3P_e/\rho`$), but between them the pair does not satisfy the consistency
identity above, and $`f`$ is taken from the degenerate branch even where the thermal pressure
dominates. (That regime is not itself an error: thermal pressure can dominate a relativistically
degenerate branch only if $`kT \gtrsim m_e c^2`$, where $`u_e = 3P_e/\rho`$ is the correct limit.
The concern is the lack of a common potential in between.) The energy expression is checked against
the paper only through its temperature derivative (above).

Its role is therefore limited to comparison: run the consistency test on it, and compare
$`\rho`$, $`\nabla_{\mathrm{ad}}`$ and $`\Gamma_1`$ with the Fermi–Dirac EOS along a model. The
differences are recorded here when measured.

### Ionisation of hydrogen and helium (D15a, plan 6.3)

*General condition.* The classical Saha equation assumes non-degenerate electrons and cannot be
combined as it stands with the Fermi–Dirac gas above. Ionisation equilibrium for each step
$`j \to j+1`$ is the equality of chemical potentials, obtained by minimising the combined free
energy of ions, electrons and radiation:

```math
\mu_j^{\mathrm{chem}} = \mu_{j+1}^{\mathrm{chem}} + \mu_e^{\mathrm{chem}}
```

The superscript marks chemical potentials; $`\mu`$, $`\mu_e`$, $`\mu_i`$ without it remain the
dimensionless molecular weights. With classical ions, the electron chemical potential
$`\mu_e^{\mathrm{chem}} = \eta\,kT`$ of the Fermi–Dirac gas (kinetic part, the same convention as
above) and the ionisation energy $`\chi_j`$ counted as the energy of the higher stage, the general
result is

```math
\frac{n_{j+1}}{n_j} = \frac{\mathcal{Z}_{j+1}}{\mathcal{Z}_j}\left(\frac{m_{j+1}}{m_j}\right)^{3/2} \exp\!\left(-\eta - \frac{\chi_j}{kT}\right)
```

with the internal partition functions $`\mathcal{Z}`$ (energies measured from each stage's ground
state) and the masses of the two stages from their translational partition functions. Two
approximations are made here and must be stated with the implementation:

* $`\mathcal{Z}_j \to g_j`$, the statistical weight of the ground state. Excited states are
  dropped; a full partition function diverges for an isolated atom and needs the same cut-off
  physics as pressure ionisation.
* $`m_{j+1}/m_j \to 1`$ (error of order $`m_e/m_u`$).

```math
\frac{n_{j+1}}{n_j} \simeq \frac{g_{j+1}}{g_j}\,\exp\!\left(-\eta - \frac{\chi_j}{kT}\right)
```

$`\eta`$, the ionisation fractions and $`n_e`$ are solved together: the free electrons released by
ionisation are the $`n_e(\eta, \vartheta)`$ of the electron gas.

*Classical limit.* For non-degenerate **and** non-relativistic electrons, $`\eta \ll -1`$ and
$`\vartheta \ll 1`$, $`e^{\eta} = \tfrac12 n_e\,(h^2/2\pi m_e kT)^{3/2}`$ and the condition becomes
the Saha equation:

```math
\frac{n_{j+1}\,n_e}{n_j} = \frac{2\,g_{j+1}}{g_j}\left(\frac{2\pi\,m_e\,k\,T}{h^{2}}\right)^{3/2} \exp\!\left(-\frac{\chi_j}{k\,T}\right)
```

| Step | $`\chi`$ [eV] | Ground levels | $`2 g_{j+1}/g_j`$ |
|---|---|---|---|
| H → H⁺ | 13.598 | H I ²S₁/₂ ($`g = 2`$) → bare ($`g = 1`$) | 1 |
| He → He⁺ | 24.587 | He I ¹S₀ ($`g = 1`$) → He II ²S₁/₂ ($`g = 2`$) | 4 |
| He⁺ → He²⁺ | 54.418 | He II ²S₁/₂ ($`g = 2`$) → bare ($`g = 1`$) | 1 |

Energies and ground levels compared with the NIST Atomic Spectra Database, ionisation energies
(retrieved 2026-11-11): 13.598434599702, 24.587389011 and 54.4177655282 eV; $`g = 2J+1`$ of the
ground level, electronic only (nuclear spin cancels in the ratios). The same weight factors
($`\tfrac12`$, 2, $`\tfrac12`$ for $`n_{j+1}/n_j`$) appear in Pols et al. (1995, eqs. 14–16,
p. 966), which use the rounded energies 13.60, 24.58 and 54.40 eV.

With ionisation fractions $`x_{\mathrm{H}}`$, $`x_{\mathrm{He}1}`$ (singly ionised) and
$`x_{\mathrm{He}2}`$ (doubly ionised):

```math
\frac{1}{\mu} = X\,(1 + x_{\mathrm{H}}) + \frac{Y}{4}\left(1 + x_{\mathrm{He}1} + 2\,x_{\mathrm{He}2}\right) + \frac{Z}{2}
```

```math
u_{\mathrm{ion}} = \frac{1}{m_u}\left[X\,x_{\mathrm{H}}\,\chi_{\mathrm{H}} + \frac{Y}{4}\left(x_{\mathrm{He}1}\,\chi_{\mathrm{He}} + x_{\mathrm{He}2}\,(\chi_{\mathrm{He}} + \chi_{\mathrm{He}^+})\right)\right]
```

$`u_{\mathrm{ion}}`$ is part of $`h`$; it is what lowers $`\nabla_{\mathrm{ad}}`$ in the ionisation
zones. The extension must stay consistent: $`\mu`$ and $`u_{\mathrm{ion}}`$ come from the same
chemical-equilibrium solution, the free electrons it yields are the $`n_e`$ of the electron gas above, and the result is
checked with the consistency test.

*Pressure ionisation.* Ideal ionisation equilibrium fails at high density: it predicts
recombination in stellar interiors, where atoms are in fact crushed. An abrupt switch to full
ionisation above a pressure threshold (the earlier proposal) is withdrawn: it makes the EOS
discontinuous and breaks the consistency identity. The remedy must be a term in the free energy (an
interaction or excluded-volume term, which appears in the equilibrium condition as a lowering of
$`\chi_j`$ and, from the same term, as corrections to $`P`$ and $`u`$). The form chosen
follows; a tabulated free-energy EOS is the fallback if its stability condition fails.

### Pressure ionisation: Pols-type free-energy term (decided)

**Decided 2026-11-11** (recommendation of the Codex review, accepted). The decision stands on the
condition stated at the end of this section: the stability scan must pass.

*Why this form.* It is analytic, needs no table, and reuses the Fermi–Dirac inversion already
required. It is empirical.

*The term.* From Pols et al. (1995, MNRAS 274, 964, §2.2.2, eqs. 26–28, p. 967), compared with
the original in task 0.3; it agrees with the reproduction of Kovetz, Yaron & Prialnik (2009, §3.1,
eq. 28). Pols et al. write it per unit mass, $`\Delta F_{\mathrm{PI}} = -N_e kT\,g_{\mathrm{PI}}`$,
and add $`\Delta F_{\mathrm{PI}}(N_e) - \Delta F_{\mathrm{PI}}(N_{e0})`$ to the free energy; the
form below is the same per volume. Let $`n`$ be the number density of free electrons and $`n_0`$
that of all electrons, bound and free; $`m_{\mathrm{H}}`$ is the hydrogen mass, as in the source:

```math
x = m_{\mathrm{H}}\,n, \qquad g(n, T) = \exp\!\left[-\left(\frac{c_1}{x}\right)^{c_2}\right]\left[\frac{13.60\ \mathrm{eV}}{kT} + \eta(n, T) + c_3\,\ln\!\left(1 + \frac{x}{c_4}\right)\right]
```

```math
F_{\mathrm{PI}} = kT\,V\left[n_0\,g(n_0, T) - n\,g(n, T)\right], \qquad (c_1, c_2, c_3, c_4) = (3\ \mathrm{g\,cm^{-3}},\ 0.25,\ 2,\ 0.03\ \mathrm{g\,cm^{-3}})
```

$`\eta(n,T)`$ is the ideal Fermi–Dirac degeneracy parameter that belongs to the density $`n`$ (the
inversion of the $`n_e`$ formula above). This is a departure from the source: Pols et al. (eq. 29,
p. 968) evaluate $`\eta`$ in this term, for both $`n`$ and $`n_0`$, from a cruder approximation to
the Fermi–Dirac integral, valid at large degeneracy, to avoid inverting it for $`n_0`$; their
equilibrium $`\eta`$ elsewhere is the full one. The four constants are kept at their published
values; Pols et al. call them adjustable and chose them for agreement with the MHD EOS on the
location of the ionisation zones.

*What follows from it* (derived in the review, re-derived here from
$`\mu = \partial F/\partial N_e`$, $`P = -\partial F/\partial V`$,
$`U = -T^2\,\partial(F/T)/\partial T`$). With $`Q(n,T) = n\,g(n,T)`$ and subscripts for partial
derivatives:

```math
\mu_{e,\mathrm{PI}}^{\mathrm{chem}} = -kT\,Q_n(n, T)
```

```math
P_{\mathrm{PI}} = kT\left[n_0\,Q_n(n_0, T) - Q(n_0, T) - n\,Q_n(n, T) + Q(n, T)\right], \qquad u_{\mathrm{PI}} = -\frac{kT^{2}}{\rho}\left[Q_T(n_0, T) - Q_T(n, T)\right]
```

$`Q_T`$ includes the temperature dependence of the inversion $`\eta(n,T)`$. Ionising one bound
electron raises $`n`$ and leaves $`n_0`$ unchanged, so the equilibrium condition becomes

```math
\frac{n_{j+1}}{n_j} = \frac{g_{j+1}}{g_j}\,\exp\!\left[-\eta - \frac{\chi_j + \mu_{e,\mathrm{PI}}^{\mathrm{chem}}}{kT}\right]
```

The shift of the ionisation energy is $`\mu_{e,\mathrm{PI}}^{\mathrm{chem}}`$ (negative at high
density). The $`\eta`$ in this condition and in the Fermi–Dirac integrals stays the ideal one; the
interaction enters only through the separate term.

*The EOS with it.* Charge neutrality, conservation of nuclei, the equilibrium conditions and the
pressure equation are solved together:

```math
P = P_{\mathrm{ion}} + P_e + P_{\mathrm{rad}} + P_{\mathrm{PI}}, \qquad u = u_{\mathrm{ion,th}} + u_e + u_{\mathrm{ion}} + \frac{a T^{4}}{\rho} + u_{\mathrm{PI}}, \qquad h = u + \frac{P}{\rho}
```

($`u_{\mathrm{ion,th}} = \tfrac32 kT/\mu_i m_u`$ is the thermal energy of the nuclei and atoms,
$`u_{\mathrm{ion}}`$ the ionisation energy defined above.) Everything comes from the one
equilibrated potential, including its $`n_0`$ term.

| Test | Expected |
|---|---|
| Dilute gas | the exponential suppresses $`F_{\mathrm{PI}}`$: ideal equilibrium and the Saha limit are recovered |
| Full ionisation, $`n = n_0`$ | $`F_{\mathrm{PI}} = P_{\mathrm{PI}} = u_{\mathrm{PI}} = 0`$: the Fermi–Dirac EOS above |
| Dense interior | no spurious recombination |
| Whole declared domain | equilibrium residuals; $`c_V > 0`$; positive compressibility; the consistency residual $`\mathcal{R}`$ |

*Condition attached to the decision.* A smooth potential does not guarantee a smooth EOS: the equilibrium can
have competing minima or lose stability. A unique stable branch is required over the supported
domain. If the scan finds a failure, it is reported and the domain is restricted; it is not hidden
by branch switching or pressure floors.

*What stays approximate.* An empirical model: no excited states, no realistic atomic interactions,
one hydrogen-like energy scale (13.60 eV) for a hydrogen–helium mixture. It does not provide
Coulomb corrections or H₂.

*Alternatives not chosen.* An excluded-volume term needs radii for every species, mixture rules
and a treatment of the electrons, which are further modelling decisions. An embedded free-energy
table is the more accurate route but adds data, redistribution questions, composition coverage and
interpolation; note that for a table, $`C^1`$ interpolation of the potential is not enough, since
pressure is already a first derivative.

### Molecular hydrogen and non-ideal terms (F2)

H₂ is a further species treated the same way (dissociation equilibrium in place of ionisation). It
matters for 0.1–0.5 M☉. Coulomb corrections are about 1% in the Sun and are not planned separately;
a tabulated EOS would include them.

---
Previous: [Structure equations](01-structure-equations.md) · [Contents](README.md) · Next: [Opacity](03-opacity.md)
