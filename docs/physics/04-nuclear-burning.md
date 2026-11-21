# 4. Nuclear burning

Energy generation rates $`\varepsilon(\rho, T)`$ in erg g⁻¹ s⁻¹, and the composition changes they
cause. In the code the rates are the `Nuclear` trait.

## Implemented

Source: `crates/stellar-core/src/physics/nuclear.rs` (rates),
`crates/stellar-core/src/physics/screening.rs` (screening, plan 9.7),
`crates/stellar-core/src/structure/evolve.rs` (composition).

The pp and CNO rates are built from Solar Fusion III (Acharya et al. 2025, Rev. Mod. Phys. 97,
035002; plan 1.6, decided 2026-11-12). The triple-α rate is the fit of Kippenhahn & Weigert (1990),
¹²C(α,γ)¹⁶O the Caughlan & Fowler (1988) fit (plan 9.5, [Carbon and oxygen](#carbon-and-oxygen-d18-plan-95-implemented)).

*Source check of Solar Fusion III* (LaTeX source of arXiv:2405.06470 read 2026-11-12): eq. 3 is the
thermal rate, eq. 4 the expansion $`S(E) \simeq S(0) + S'(0)E + \tfrac12 S''(0)E^2`$ (E in MeV), and
Table I lists:

| Reaction | $`S(0)`$ [MeV b] | $`S'(0)`$ [b] | $`S''(0)`$ [b MeV⁻¹] | Text |
|---|---|---|---|---|
| ¹H(p,e⁺ν)²H | $`4.09\times10^{-25}`$ | $`4.5\times10^{-24}`$ | $`9.9\times10^{-23}`$ | §III.F: $`S'/S = 11.0 \pm 0.2`$ MeV⁻¹, $`S''/S = 242 \pm 72`$ MeV⁻² (corrections of about 7% and 0.5% at the solar Gamow peak) |
| ³He(³He,2p)⁴He | 5.21 | −4.90 | 22.42 | §IV.B: $`S(E) = 5.21 - 4.90E + 11.21E^2`$ MeV b, so the ½ of eq. 4 is confirmed |
| ³He(α,γ)⁷Be | $`5.61\times10^{-4}`$ | $`-3.03\times10^{-4}`$ | – | §V.E: $`S'/S = -0.54 \pm 0.07`$ MeV⁻¹ for E = 0–500 keV; "second derivative consistent with zero" |
| ¹⁴N(p,γ)¹⁵O | $`1.68\times10^{-3}`$ | – | – | §IX: sum of the partial $`S(0)`$ (Table XI), 1.68 ± 0.14 keV b; no derivatives recommended |

The same review confirms weak (Salpeter) screening for $`Z_1 Z_2 < 10`$ in the solar core (§X.B).
The all-regime screening coded in plan 9.7 departs from it by −12% for p + ¹⁴N at the solar centre
(see [the measurement](#screening-prescription-for-all-regimes-plan-97-implemented)).

### Rates

*Thermal rates* (eq. 3 of Solar Fusion III), for p + p, ³He + ³He, ³He + ⁴He and p + ¹⁴N:

```math
\langle\sigma v\rangle = \sqrt{\frac{8}{\pi\mu}}\,(kT)^{-3/2} \int_0^\infty S(E)\, \exp\!\left(-\frac{E}{kT} - \frac{b}{\sqrt{E}}\right) dE, \qquad b = \pi\alpha_{\mathrm{fs}} Z_1 Z_2 \sqrt{2\mu c^2}, \qquad S(E) = S(0) + S'(0)\,E + \tfrac12 S''(0)\,E^2
```

* Reduced masses from atomic masses (¹H 1.007825, ³He 3.016029, ⁴He 4.002603, ¹⁴N 14.003074 u), as
  in the 1990 fits (their exponent 33.80 for p + p needs $`\mu = 0.5039`$; with $`\mu = 1/2`$ it
  would be 33.72, a 4% change of the rate).
* The integral is computed by composite Simpson quadrature over ±12 widths of the Gamow peak,
  tabulated once at first use on 161 points in $`T_6^{-1/3}`$ for $`1 \le T_6 \le 200`$, and
  interpolated with cubic Hermite polynomials (C¹, rule 3); outside the table it is extrapolated
  linearly in $`T_6^{-1/3}`$. The table reproduces the quadrature to better than $`10^{-4}`$ in
  $`\ln\langle\sigma v\rangle`$ (unit test).
* $`S(E)`$ is kept non-negative. The ³He + ⁴He slope is valid to 500 keV and the ¹⁴N(p,γ) rate has
  no resonances: both hold for hydrogen burning, below about $`T_6 = 100`$.

*pp chain.* The pp I rate times the branching factor $`\psi`$ (next section):

```math
\varepsilon_{pp} = \psi\,\varepsilon_{pp,\mathrm{I}}, \qquad \varepsilon_{pp,\mathrm{I}} = \frac{r_{11}}{2}\,\frac{E_{\mathrm{I}}}{\rho}, \qquad r_{11} = \tfrac12\, n_p^2\, \langle\sigma v\rangle_{11}\, f_1, \qquad n_p = \frac{\rho X}{m_{\mathrm{H}}}
```

$`E_{\mathrm{I}} = 4 m_{\mathrm{H}}\, Q_{pp,\mathrm{I}} = 26.20`$ MeV is the energy deposited per completed
pp I chain (two p + p reactions).

*CNO cycle.* The slowest step, p + ¹⁴N, sets the rate; in CN equilibrium all catalysts are counted
as ¹⁴N:

```math
\varepsilon_{\mathrm{CNO}} = n_p\, n_{14}\, \langle\sigma v\rangle_{1,14}\, f_7\, \frac{E_{\mathrm{CNO}}}{\rho}, \qquad n_{14} = \frac{\rho X_{\mathrm{CNO}}}{m_{14}}, \qquad X_{\mathrm{CNO}} = f_{\mathrm{CNO}}\, Z
```

$`E_{\mathrm{CNO}} = 4 m_{\mathrm{H}}\, Q_{\mathrm{CNO}} = 25.03`$ MeV per cycle. (Fowler et al. 1975 quote
24.970 MeV; the 0.25% difference comes from the neutrino energies. Taking both energies from the
$`Q`$ values of the burning step keeps the fuel accounting exact.)

*Triple-α* (Kippenhahn & Weigert 1990)

```math
\varepsilon_{3\alpha} = 5.09\times10^{11}\;\rho^2\, Y^3\, f_{4}\, f_{8}\; T_8^{-3}\, \exp\!\left(-\frac{44.027}{T_8}\right)
```

**[measured]** against the replaced 1990 fits (unscreened, pp I only; unit test): the Solar Fusion III
p + p rate is 1.050, 1.054, 1.057, 1.065 times the fit at $`T_6 = 10`$, 15, 20, 30 (estimated "about
8%" before $`S'`$ and $`S''`$ were included); the ¹⁴N(p,γ) rate is 0.520, 0.527, 0.534, 0.546 times
the fit (estimated 0.52–0.56 with constant $`S`$, below).

*The replaced 1990 fits, for reference.* Before plan 1.6 the code used, as given by Kippenhahn &
Weigert (1990),
$`\varepsilon_{pp} = 2.38\times10^{6}\rho X^2 f_1 g_{11} T_6^{-2/3}\exp(-33.80\,T_6^{-1/3})`$ with
$`g_{11} = 1 + 0.0123 T_6^{1/3} + 0.0109 T_6^{2/3} + 0.0009 T_6`$, and
$`\varepsilon_{\mathrm{CNO}} = 8.67\times10^{27}\rho X X_{\mathrm{CNO}} f_7 g_{14,1} T_6^{-2/3}\exp(-152.28\,T_6^{-1/3})`$
with $`g_{14,1} = 1 + 0.0027 T_6^{1/3} - 0.00778 T_6^{2/3} - 0.000149 T_6`$. The book was not
available for the source check of 2026-11-11; the fits were compared with the compilations they are
built on (originals read on ADS; $`T_9^{1/3} = 0.1\,T_6^{1/3}`$, $`T_9 = 0.1\,T_8`$):

| Fit | Compilation entry | Exponent and correction polynomial | Prefactor re-derived |
|---|---|---|---|
| pp | Harris et al. (1983), ARA&A 21, 165, Table 1 (p. 168), H1(P,E+NU)H2: $`3.82\times10^{-15}`$ (same polynomial and exponent in Fowler et al. 1975, Table 1, p. 72, with 4.21e-15) | 33.80 and 0.0123, 0.0109, 0.000938 → 0.0009: agree | $`2.377\times10^{6}`$ with $`n_p = \rho X/m_{\mathrm{H}}`$ and 13.1 MeV per p + p reaction (half of 26.2 MeV); code 2.38e6 |
| CNO | Fowler, Caughlan & Zimmerman (1975), ARA&A 13, 69, Table 1 (p. 78), N14(P,G)O15: $`5.08\times10^{7}`$, "Q = 24.970 per N14(P,G)O15 reaction for full CN-cycle" | 152.28 and 0.0027, −0.00778, −0.000149: agree; the source's further terms ($`T_9^{4/3}`$, $`T_9^{5/3}`$, $`\exp[-(T_9/3.090)^2]`$, resonances) are dropped | $`8.674\times10^{27}`$ with $`n_p = \rho X/m_{\mathrm{H}}`$, $`n_{14} = \rho X_{\mathrm{CNO}}/(14\,m_u)`$ and 24.97 MeV; code 8.67e27 |
| triple-α | Harris et al. (1983), Table 1 (p. 168), HE4(2A,G)C12: $`2.79\times10^{-8}\,T_9^{-3}\exp(-4.4027/T_9)`$, Q = 7.275 (the same high-temperature term in the ORNL reproduction of Caughlan & Fowler 1988); the 1975 entry (p. 74) had 2.49e-8 and 4.4109 | 44.027: agrees | $`5.100\times10^{11}`$ with $`\rho^2Y^3N_A/384`$ and 7.275 MeV; code 5.09e11 (still used) |

### Electron screening

All regimes (plan 9.7): Chugunov & DeWitt (2009) for the binary reactions, with the electron factor
of Potekhin & Chabrier (2012) and the static two-stage form for triple-α; see
[Screening prescription for all regimes](#screening-prescription-for-all-regimes-plan-97-implemented).
Its weak limit is the formula used before plan 9.7 (Salpeter), which had its exponent capped at 2:

```math
f_{Z_1 Z_2} = \exp\!\left(0.188\, Z_1 Z_2\, \zeta\, \rho^{1/2}\, T_6^{-3/2}\right), \qquad \zeta^2 = 2X + 1.5\,Y + 4Z
```

| Reaction | $`Z_1 Z_2`$ |
|---|---|
| p + p | 1 |
| p + ¹⁴N | 7 |
| α + α | 4 |
| α + ⁸Be | 8 |
| ¹²C + α | 12 |

### Screening prescription for all regimes (plan 9.7, implemented)

**Decided 2026-11-11** (recommendation of the Codex review, accepted). Two parts of it were not
from a published source. *Settled 2026-11-12:* the electron-response extension is replaced by the
published form of Potekhin & Chabrier (2012, eqs. 20–21); the product form for triple-α is kept,
with MESA's implementation as precedent (details below). Chosen over Graboske et al. (1973) because
one analytic interpolation is easier to keep smooth than joined branches; Graboske's prescription
is therefore not needed and was not compared.

*Binary reactions:* Chugunov & DeWitt (2009, eq. A4, with eq. 24 for the plasma free energy)
**[compared]**. With number fractions $`x_s = n_s/\sum n_s`$ of the actual ion mixture and
$`\langle Z^p\rangle = \sum_s x_s Z_s^p`$:

```math
a_e = \left(\frac{3}{4\pi n_e}\right)^{1/3}, \qquad \Gamma_e = \frac{e^{2}}{a_e\,kT}, \qquad \Gamma_s = Z_s^{5/3}\,\Gamma_e, \qquad \Gamma_{ij} = \frac{2\,Z_i Z_j}{Z_i^{1/3} + Z_j^{1/3}}\,\Gamma_e
```

*Free energy of a one-component plasma* (fit of Potekhin & Chabrier):

```math
f_0(\Gamma) = A_1\left[\sqrt{\Gamma\,(A_2 + \Gamma)} - A_2\,\operatorname{asinh}\sqrt{\Gamma/A_2}\right] + 2A_3\left[\sqrt{\Gamma} - \arctan\sqrt{\Gamma}\right] + B_1\left[\Gamma - B_2\ln\!\left(1 + \frac{\Gamma}{B_2}\right)\right] + \frac{B_3}{2}\ln\!\left(1 + \frac{\Gamma^{2}}{B_4}\right)
```

| $`A_1`$ | $`A_2`$ | $`A_3`$ | $`B_1`$ | $`B_2`$ | $`B_3`$ | $`B_4`$ |
|---|---|---|---|---|---|---|
| −0.907 | 0.62954 | $`-\sqrt{3}/2 - A_1/\sqrt{A_2}`$ | 0.00456 | 211.6 | $`-10^{-4}`$ | 0.00462 |

*Screening exponent.* With the compound charge $`Z_c = Z_i + Z_j`$:

```math
H_{ij} = \frac{C_{ij} + \Gamma_{ij}^{2}}{1 + \Gamma_{ij}^{2}}\left[f_0\!\left(\frac{\Gamma_i}{t_{ij}}\right) + f_0\!\left(\frac{\Gamma_j}{t_{ij}}\right) - f_0\!\left(\frac{\Gamma_c}{t_{ij}}\right)\right], \qquad f_{ij} = e^{H_{ij}}
```

```math
C_{ij} = \frac{3\,Z_i Z_j\,\sqrt{\langle Z^{2}\rangle/\langle Z\rangle}}{Z_c^{5/2} - Z_i^{5/2} - Z_j^{5/2}}
```

*Quantum correction for non-resonant reactions.* With the reduced mass $`\mu_{ij}`$ of the two
nuclei:

```math
\tau_{ij} = \left(\frac{27\pi^{2}\,\mu_{ij}\,Z_i^{2}Z_j^{2}\,e^{4}}{2\,kT\,\hbar^{2}}\right)^{1/3}, \qquad z_{ij} = \frac{3\,\Gamma_{ij}}{\tau_{ij}}, \qquad y_{ij} = \frac{4\,Z_i Z_j}{Z_c^{2}}
```

```math
t_{ij} = \left[1 + 0.013\,y_{ij}^{2}\,z_{ij} + 0.406\,y_{ij}^{0.14}\,z_{ij}^{2} + \left(0.062\,y_{ij}^{0.19} + \frac{1.8}{\Gamma_{ij}}\right)z_{ij}^{3}\right]^{1/3}
```

$`z_{ij}`$ is $`\zeta_{ij}`$ of Chugunov & DeWitt (2009, eq. 12) and $`t_{ij}`$ their eq. 25; the same
expression is printed as eq. 22 of Potekhin & Chabrier (2012, A&A 538, A115), whose footnote 1 on
the substitution $`\Gamma_j \to \Gamma_j/t_{ij}`$ corrects a typo of the 2009 paper. Both agree with
the form above (checked 2026-11-11).

*Electron response in the weak limit* (an extension proposed by the review, **not** in the
published formula; no publication of it was found). *Decided 2026-11-12: not adopted;* the
published alternative of Potekhin & Chabrier (2012) below is the one to code. The extension is kept
here only as the reasoning that led to it. The published $`C_{ij}`$
assumes a rigid electron background (Chugunov & DeWitt 2009, App. A; Potekhin & Chabrier 2012, §1),
so it does not reproduce the electron part of the implemented Salpeter formula. To recover the
correct weak limit for classical or degenerate electrons, replace $`C_{ij}`$ by
$`H_{\mathrm{DH}}/H_{\mathrm{LM,DH}}`$:

```math
D^{-2} = 4\pi e^{2}\left[\frac{\sum_s n_s Z_s^{2}}{kT} + \left(\frac{\partial n_e}{\partial \mu_e^{\mathrm{chem}}}\right)_{T}\right], \qquad H_{\mathrm{DH}} = \frac{Z_i Z_j\,e^{2}}{kT\,D}, \qquad H_{\mathrm{LM,DH}} = \frac{\Gamma_e^{3/2}}{\sqrt{3}}\left[Z_c^{5/2} - Z_i^{5/2} - Z_j^{5/2}\right]
```

$`\partial n_e/\partial\mu_e^{\mathrm{chem}}`$ comes from the Fermi–Dirac EOS. This fades out at
strong coupling and does not model electron polarisation there.

What is published and agrees: $`D^{-2}`$ with the electron term $`k_{\mathrm{TF}}^{2} = 4\pi e^{2}\,\partial n_e/\partial\mu_e`$
and $`H_{\mathrm{DH}}`$ are eqs. 14–17 of Potekhin & Chabrier (2012); $`H_{\mathrm{LM,DH}}`$ is eq. A2 of
Chugunov & DeWitt (2009). What is not published: using their ratio in place of $`C_{ij}`$.

*Published alternative* (Potekhin & Chabrier 2012, eqs. 20–21, read 2026-11-11). Multiply a
rigid-background exponent $`h_{\mathrm{ii}}`$ by an electron-screening factor:

```math
h_{\mathrm{appr}} = h_{\mathrm{ii}}\,\sqrt{1 + k_{\mathrm{TF}}^{2}\,\tilde D_{\mathrm{ion}}^{2}}, \qquad \tilde D_{\mathrm{ion}} = \sqrt{D_{\mathrm{ion}}^{2} + (0.6\,a_{ij})^{2}}, \qquad D_{\mathrm{ion}}^{-2} = \frac{4\pi e^{2}}{kT}\sum_s n_s Z_s^{2}, \qquad a_{ij} = \frac{a_i + a_j}{2}
```

with $`a_s = Z_s^{1/3} a_e`$ and 0.6 their only fitted parameter. In the weak limit
$`h_{\mathrm{ii}} \to Z_iZ_je^{2}/(kT D_{\mathrm{ion}})`$ and the factor is $`D_{\mathrm{ion}}/D`$, so it
gives the same limit $`H_{\mathrm{DH}}`$; unlike the review's extension it keeps an electron
polarisation term at strong coupling. Potekhin & Chabrier used it with their own rigid-background
$`h_{\mathrm{ii}}`$ (their eq. 19 or the fit of Chugunov & DeWitt 2009b), not with eq. A4 of
Chugunov & DeWitt (2009a); applying it to $`H_{ij}`$ above is a small extrapolation of the source.
*Decided 2026-11-12:* this published form is coded in plan 9.7, applied to $`H_{ij}`$ above; the
extrapolation (CD2009 eq. A4 instead of PC2012's own $`h_{\mathrm{ii}}`$) is stated in the code
comment and tested against both limits in the table below.

*Use.* $`f_{11}`$ for pp and $`f_{1,7}`$ for the CNO bottleneck, with the actual mixture
(including produced carbon, plan 9.5). The exponent cap of 2 is removed; $`H`$ is kept as
a logarithm internally.

*Triple-α* (resonant). The static form, $`t_{ij} = 1`$, applied once to the unscreened rate:

```math
f_{3\alpha} \approx \exp\!\left[H_{2,2}^{\mathrm{static}} + H_{4,2}^{\mathrm{static}}\right]
```

The first term shifts the equilibrium abundance of ⁸Be, the second the capture of the third α
through the Hoyle state. The basis is the difference of Coulomb chemical potentials,
$`\ln f_{3\alpha} = (3\mu_\alpha^{C} - \mu_{\mathrm{C}}^{C})/kT`$ (Kushnir, Waxman & Chugunov 2019,
MNRAS 486, 449, eqs. 1, 7–8 **[compared]**; they credit this triple-α factor to Yakovlev &
Shalybkov 1989, Astrophys. Space Phys. Rev. 7, 311, not read), valid for non-degenerate ions with
$`\Gamma/\tau \ll 1`$ and a small Coulomb contribution to the electron free energy (their
assumptions i–iii, §2). In the strong limit the two stages telescope exactly to
$`3f_0(\Gamma_2) - f_0(\Gamma_6)`$, and in the weak limit they add to $`(4 + 8)\,e^2/(kTD)`$, which
is the Debye–Hückel value of $`(3\mu_\alpha^{C} - \mu_{\mathrm{C}}^{C})/kT`$ (derivation). At
intermediate coupling the product of two interpolation factors is an approximation. *Decided
2026-11-12: adopted.* No paper printing the product form was found, but it is what MESA does:
`eval_screen_triple` in `net/private/net_screen.f90` (MESAHub/mesa `main`, read 2026-11-12) screens
⁴He + ⁴He + ⁴He as the product of the pair factors for (⁴He, ⁴He) and (⁴He, $`Z = 4`$, $`A = 8`$),
with any pair prescription, including Chugunov's. It is exact in both limits (above). The factors
$`(C_{ij} + \Gamma_{ij}^2)/(1 + \Gamma_{ij}^2)`$ of the two stages use different $`\Gamma_{ij}`$, so the
sum is not exactly a difference of chemical potentials there. Evaluating the chemical-potential
difference directly from one free-energy expression would avoid this; that, too, is not from a
publication. The
non-resonant correction $`t_{ij}`$ must not be applied to this chain.

| Test | Expected |
|---|---|
| Weak coupling | $`H_{ij} \to Z_iZ_je^2/(kTD)`$; with classical electrons, the implemented Salpeter formula |
| Triple-α, weak | $`H_{3\alpha} \to 12\,e^2/(kTD)`$ |
| Strong, classical liquid | $`H_{ij} \approx 0.9\,\Gamma_e\left[Z_c^{5/3} - Z_i^{5/3} - Z_j^{5/3}\right]`$ |
| Dilute | all factors → 1; derivatives continuous through intermediate coupling |

*What stays approximate.* No tunnelling correction for the resonant chain, no plasma change of
resonance widths, no pycnonuclear burning. Screening of this kind does not make the ideal-ion EOS
adequate at strong coupling: Coulomb terms in the EOS are a separate requirement. Whether the gate
on helium ignition can be lifted with this prescription is decided after it is implemented and
tested, not by adopting it.

*Implementation* (`crates/stellar-core/src/physics/screening.rs`, plan 9.7, 2026-11-21). As
above, with these choices:

* **Electrons.** $`k_{\mathrm{TF}}^2 = 4\pi e^2\,\partial n_e/\partial\mu_e`$ uses the electron model
  of the EOS: $`n_e/kT`$ (classical) for an EOS without Fermi–Dirac electrons, so that with the ideal
  EOS the weak limit is exactly the Salpeter formula; otherwise the Fermi–Dirac integrals at the
  EOS's $`\eta`$, $`(\partial n_e/\partial\eta)/kT`$ (`Nuclear::eps_eta`). $`n_e`$, $`\Gamma_e`$ and
  $`\langle Z\rangle`$, $`\langle Z^2\rangle`$ are for full ionisation of H, He, the GN93 metals
  and the produced C and O.
* **Species.** p + p and p + ¹⁴N with $`t_{ij}`$ (reduced masses from the atomic masses); ³He + ³He
  and ³He + ⁴He enter $`\alpha`$ of the pp branching as $`2H_{34} - H_{11} - H_{33}`$; triple-α
  static, $`H_{2,2} + H_{4,2}`$; ¹²C + α with $`t_{ij}`$.
* **No cap on $`H`$.** Rates are formed as logarithms, $`\ln\varepsilon = \ln(\ldots) + \ln\langle\sigma v\rangle + H`$.
  Each channel's $`\varepsilon`$ is limited to $`10^{30}`$ erg g⁻¹ s⁻¹, far above any real model
  (a guard against overflow in unphysical Newton trial states, such as $`10^3`$ K at
  $`10^{-2}`$ g cm⁻³ with classical electrons, where $`H \approx 10^3`$).
* **Cost.** The charge-dependent constants of each pair are computed once, the last eight states
  are remembered (the finite-difference Jacobian repeats $`(\rho, T)`$ when it perturbs $`r`$ or
  $`L`$), and only reactions that can run are screened (hydrogen burning needs $`X > 0`$,
  ¹²C(α,γ) produced carbon, triple-α $`T \ge 3\times10^7`$ K, where its rate is below
  $`10^{-50}`$ erg g⁻¹ s⁻¹). **[measured]** A ZAMS model costs about 1.6× the time of the weak
  formula (1 M☉: 1.8 s against 1.1 s; 3 M☉: 10.9 s against 6.1 s, release build).

*Tests* (unit tests of `screening.rs` and `nuclear.rs`): the free energy has its Debye–Hückel and
strong limits ($`f_0(1) = -0.4373`$, the Monte Carlo value of the one-component plasma is −0.436);
weak limit equal to the Salpeter formula with $`\zeta`$ from the same species list, and to the
pre-9.7 formula for $`Z = 0`$ (to 2×10⁻³ at $`\Gamma_e \approx 2\times10^{-5}`$); triple-α weak
limit $`12e^2/(kTD)`$; strong limit $`0.9\,\Gamma_e[Z_c^{5/3} - Z_i^{5/3} - Z_j^{5/3}]`$ within 3% and
the static triple-α telescoping to $`3f_0(\Gamma_2) - f_0(\Gamma_6)`$ within 1%; $`d\ln H/d\log\rho`$
continuous from $`10^{-6}`$ to $`10^{7}`$ g cm⁻³ at $`10^8`$ K; degenerate electrons screen less than
classical ones.

**[measured]** *Solar centre* (present Sun, $`\rho = 150`$, $`T = 1.57\times10^7`$ K, $`X = 0.34`$,
classical electrons): the rates change from the weak formula by −0.8% (pp) and −12% (CNO); stated
tolerances 1% and 15% (unit test). The CNO change is the intermediate coupling of ¹⁴N itself,
$`\Gamma_{\mathrm N} = 7^{5/3}\Gamma_e = 1.7`$ at $`\Gamma_e = 0.067`$, where $`f_0`$ is well below its
Debye–Hückel form (−0.87 against −1.29); the rigid-background exponent of p + ¹⁴N is 0.156 against
0.269 in the weak limit, the electron factor 1.33 against 1.27. This is larger than the agreement
with Salpeter screening for $`Z_1Z_2 < 10`$ in the solar core that Solar Fusion III reports (§X.B,
above): the linear-mixing form of Chugunov & DeWitt (2009) treats a trace ion by its own ion-sphere
coupling. **Open item:** compare p + ¹⁴N at the solar centre with a published intermediate-screening
calculation for trace ions before the CNO rate of the Sun is quoted to better than 15%.
*ZAMS models* (1 M☉): $`L`$ +0.07%, $`R`$ −0.05%, $`T_c`$ +0.07% against the weak formula.

### Helium ignition gate

The weak (Debye–Hückel) formula used before plan 9.7 was valid only while its exponent was small,
and its cap at 2 was a numerical safeguard. In the dense, degenerate helium core of a red giant the
triple-α reactions are outside the weak regime. Plan 9.7 replaced it by the prescription of the
next section, which covers weak, intermediate and strong screening and takes the electron
degeneracy from the EOS; its tests pass (`physics/screening.rs`).

The gate itself is **not lifted** (lifting it is a separate decision, plan 9.7): with a Fermi–Dirac
EOS an "as far as possible" run still stops when the core becomes degenerate ($`\eta_c > 1`$), and
the model must not be read as predicting helium ignition (core mass, luminosity at the tip, the
flash). Thermal-neutrino losses (plan 9.6) are still missing, and they set where a degenerate
core ignites. The largest screening exponent of every model and run is still reported, with a
warning above 0.3 ([Screening diagnostic](#screening-diagnostic-plan-111)).

### pp branching factor (D11, plan 1.6)

Let $`r_{11}`$, $`r_{33}`$, $`r_{34}`$ be the rates of p + p, ³He + ³He and ³He + ⁴He. With ³He in
equilibrium, $`r_{11} = 2 r_{33} + r_{34}`$. Helium production relative to the pp I branch alone:

```math
\phi = \frac{2\,(r_{33} + r_{34})}{r_{11}} = 1 + \frac{2}{\sqrt{1 + 2/\alpha} + 1}, \qquad \alpha = \left(\frac{n_4}{n_p}\right)^{2} \frac{\langle\sigma v\rangle_{34}^{\,2}}{\langle\sigma v\rangle_{11}\,\langle\sigma v\rangle_{33}}, \qquad \frac{n_4}{n_p} = \frac{Y}{4X}
```

The form of $`\phi`$ is the stable equivalent of $`1 - \alpha + \sqrt{\alpha(\alpha + 2)}`$; $`\phi`$
runs from 1 to 2. The fraction of helium made through pp II, and the energy factor (pp II gives
25.66 MeV instead of 26.20 MeV to the gas; pp III is neglected, safe below $`T_6 \approx 20`$):

```math
F_{\mathrm{II}} = \frac{2\,(\phi - 1)}{\phi}, \qquad \psi = \phi\,\left(1 - 0.021\,F_{\mathrm{II}}\right), \qquad \varepsilon_{pp} = \psi\,\varepsilon_{pp,\mathrm{I}}
```

*How $`\alpha`$ is computed.* From its definition, with the thermal rates of the previous section
(one evaluation, Solar Fusion III with $`S'`$ and $`S''`$), the number ratio from the atomic masses,
$`n_4/n_p = (Y/m_4)/(X/m_{\mathrm{H}}) \approx Y/4X`$, and the screening factors
$`f_{34}^2/(f_{11} f_{33})`$.

* **[measured]** pp II share in the models: centre of the ZAMS Sun $`\phi = 1.025`$
  ($`F_{\mathrm{II}} = 0.048`$); centre of the Sun at 4.57 Gyr $`\phi = 1.096`$ ($`F_{\mathrm{II}} = 0.175`$).
  Reported, not tested against a range (the analytic Sun is not a solar model).
*Check value.* The leading-order approximation of $`\alpha`$ (constant $`S`$, no screening, leading
Gamow term), with $`\tau_{ij} = 42.49\,(Z_i^2 Z_j^2\,\mu_{ij})^{1/3}\,T_6^{-1/3}`$ and
$`2\tau_{34} - \tau_{11} - \tau_{33} = 99.98\,T_6^{-1/3}`$, is

```math
\alpha \simeq C_\alpha \left(\frac{Y}{4X}\right)^{2} \exp\!\left(-100\,T_6^{-1/3}\right), \qquad C_\alpha \simeq \frac{S_{34}^{2}}{S_{11}\,S_{33}}\left(\frac{4\,\mu_{11}\,\mu_{33}}{\mu_{34}^{2}}\right)^{1/3} \simeq 1.49\times10^{17}
```

It is used only as a check of the quadrature: with constant $`S(0)`$ and integer mass numbers the
quadrature gives $`C_\alpha = 1.47`$, 1.46, $`1.45\times10^{17}`$ at $`T_6 = 10`$, 15, 30 (unit test,
as computed for these pages on 2026-11-11). The value $`1.2\times10^{17}`$ once quoted from
memory of Clayton (1968) is not used (decided 2026-11-12, attribution dropped).

### Composition changes

The hydrogen consumption follows each channel's energy per gram (D12, plan 1.5; with pp branching,
plan 1.6). $`\phi`$ counts helium nuclei and $`\psi`$ energy, so the energy delivered per gram of
hydrogen burnt by the pp chain depends on the branch:

```math
\frac{\partial X}{\partial t} = -\frac{\varepsilon_{pp}}{Q_{pp}^{\mathrm{eff}}} - \frac{\varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{CNO}}}, \qquad Q_{pp}^{\mathrm{eff}} = Q_{pp,\mathrm{I}}\,\frac{\psi}{\phi} = Q_{pp,\mathrm{I}}\left(1 - 0.021\,F_{\mathrm{II}}\right)
```

```math
Q_{pp,\mathrm{I}} = 6.27\times10^{18}, \qquad Q_{\mathrm{CNO}} = 5.99\times10^{18}\ \mathrm{erg\,g^{-1}}
```

```math
\frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}} - \frac{1}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad Q_{\mathrm{He}} = 5.85\times10^{17}\ \mathrm{erg\,g^{-1}}
```

and the carbon and oxygen of [Carbon and oxygen](#carbon-and-oxygen-d18-plan-95-implemented);
$`Z = Z_0`$ is not changed by burning.

* 26.731 MeV per ⁴He is $`6.398\times10^{18}`$ erg/g. After neutrino losses the pp I chain deposits
  26.20 MeV ($`Q_{pp,\mathrm{I}}`$) and the CNO cycle 25.03 MeV ($`Q_{\mathrm{CNO}}`$).
* $`Q_{\mathrm{He}}`$ is 7.275 MeV per carbon nucleus.
* Keeping the fixed $`Q_{pp,\mathrm{I}}`$ with the branching-corrected rate would underestimate the
  hydrogen consumption by up to 2.1% ($`F_{\mathrm{II}} = 1`$).
* Carbon and oxygen from helium burning are tracked separately (plan 9.5; before, carbon was counted
  in $`Z`$ and the catalysts kept their initial value). The catalysts are
  $`X_{\mathrm{CNO}} = f_{\mathrm{CNO}}Z_0 + X_{\mathrm C}`$.
* The fits assume that ³He and the CNO isotopes are in equilibrium.
* **[measured]** effect of the per-channel $`Q`$ alone: main sequence 4.3% shorter for 3 M☉ (285.89 →
  273.52 Myr) and 4.4% for 5 M☉ (76.09 → 72.71 Myr).

### Carbon and oxygen (D18, plan 9.5, implemented)

With $`Q_{3\alpha} = Q_{\mathrm{He}}`$ per gram of helium and
$`Q_{\mathrm{C}\alpha} = 7.162\ \mathrm{MeV}/12\,m_u = 5.76\times10^{17}`$ erg per gram of carbon
consumed by ¹²C(α,γ)¹⁶O:

```math
\frac{\partial X_{\mathrm{C}}}{\partial t} = \frac{\varepsilon_{3\alpha}}{Q_{3\alpha}} - \frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad \frac{\partial X_{\mathrm{O}}}{\partial t} = \frac{4}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}
```

```math
\frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{3\alpha}} - \frac{1}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}
```

*Energies (checked 2026-11-11).* With the AME2020 mass excesses (Wang et al. 2021, file
`mass_1.mas20` of the IAEA AMDC: ⁴He 2424.916 keV, ¹²C 0, ¹⁶O −4737.002 keV), $`Q`$ of ¹²C(α,γ)¹⁶O is
7161.9 keV and $`Q_{3\alpha} = 3 \times 2424.916 = 7274.7`$ keV, as used above. Per gram:
$`7.162\ \mathrm{MeV}/12\,m_u = 5.759\times10^{17}`$ erg g⁻¹ (carbon) and
$`7.275\ \mathrm{MeV}/(3 \times 4.0026\,m_u) = 5.846\times10^{17}`$ erg g⁻¹ (helium), consistent with
the rounded values. The mass factors 4/3 and 1/3 use mass numbers (12 + 4 → 16); the mass defect is
below 10⁻³ of the mass.

The rate $`\varepsilon_{\mathrm{C}\alpha}`$ is the Caughlan & Fowler (1988, At. Data Nucl. Data
Tables 40, 283–334) fit (*decided 2026-11-12*). The printed table was not available; the fit is
taken from two independent reproductions that agree term by term: the ORNL Fortran version of the
compilation and MESA's `rate_c12ag_fxt` (`rates/private/ratelib.f90`, MESAHub/mesa `main`, read
2026-11-12; its $`q_1 = 1/12.222016 = 1/3.496^2`$). The CF88 fit as reproduced in the ORNL Fortran version of the
compilation (reaction 62, C12(A,G)O16, Q = 7.162; read via the Internet Archive copy of
`www.phy.ornl.gov/astrophysics/data/cf88/fortran_rates.html`) is

```math
N_A\langle\sigma v\rangle = \frac{1.04\times10^{8}}{T_9^{2}\,(1 + 0.0489\,T_9^{-2/3})^{2}}\,e^{-32.120\,T_9^{-1/3} - (T_9/3.496)^{2}} + \frac{1.76\times10^{8}}{T_9^{2}\,(1 + 0.2654\,T_9^{-2/3})^{2}}\,e^{-32.120\,T_9^{-1/3}} + 1.25\times10^{3}\,T_9^{-3/2}\,e^{-27.499/T_9} + 1.43\times10^{-2}\,T_9^{5}\,e^{-15.541/T_9}
```

in cm³ mol⁻¹ s⁻¹. With $`n_{\mathrm{C}} = \rho X_{\mathrm{C}}/(12\,m_u)`$ and $`n_\alpha = \rho Y/(4\,m_u)`$ this gives
$`\varepsilon_{\mathrm{C}\alpha} = \rho\,X_{\mathrm{C}}\,Y\,N_A\,\big(N_A\langle\sigma v\rangle\big)\,Q/48`$ (derivation). This is a transcription from a reproduction, not from the
printed table; two reproductions agreeing makes a transcription error unlikely, and the printed
table should still be compared when it becomes available. Harris et al. (1983, p. 167) already
noted that the first term of their own fit should be changed to 1.58e8, a reminder that this rate
has been revised repeatedly.

*Uncertainty, as a parameter.* The rate is coded with a multiplier $`f_{\mathrm{C}\alpha}`$ (default
1, i.e. CF88). MESA's `rate_c12ag_fxt` applies 1.7 to the same fit (comment "1.7 times cf88
value"), and MESA also carries the Kunz et al. (2002, ApJ 567, 643) fit as `rate_c12ag_kunz`; the
spread between these is the size of the present uncertainty. A modern evaluation (for example
deBoer et al. 2017, Rev. Mod. Phys. 89, 035007) can replace the default after 9.5 without changing
the interface.

*What the composition vector must then mean.* Tracking C and O is not only two more equations. The
following have to be fixed together, or the pieces contradict each other:

| Quantity | Definition |
|---|---|
| $`Z_0`$ | initial heavy elements, fixed scaled mixture; sets the opacity table and the initial CNO catalysts |
| $`X_{\mathrm{CNO}}`$ | catalysts of the CNO cycle: $`f_{\mathrm{CNO}} Z_0`$. Carbon made by helium burning is **not** added while it sits in hydrogen-free layers; it counts only where it is mixed into hydrogen |
| $`X_{\mathrm{C}}`$, $`X_{\mathrm{O}}`$ | carbon and oxygen produced by helium burning, in addition to $`Z_0`$ |
| $`\mu`$, $`\mu_e`$, $`\mu_i`$ | computed from all species (C: $`A=12`$, 6 electrons; O: $`A=16`$, 8 electrons), used by the EOS and by $`\nabla_\mu`$ in the Ledoux criterion and $`N^2`$ |
| Screening $`\zeta`$ | from the same species list |
| Opacity | see the mixture note in [Opacity](03-opacity.md) |

*Implementation* (plan 9.5, 2026-11-21). `Composition` carries `x_c`, `x_o` besides X, Y and
$`Z = Z_0`$; `Profile` carries them per point (`x_c`, `x_o`, and `eps_ca`); `Globals` the central
values (`c_c`, `o_c`) and the share of ¹²C(α,γ) in the nuclear luminosity (`frac_ca`). The CLI
prints the central C/O ratio.

* $`X_{\mathrm{CNO}} = f_{\mathrm{CNO}}Z_0 + X_{\mathrm C}`$ at every point
  (`Composition::cno_catalysts`). As the CNO rate is proportional to $`X`$, the produced carbon
  matters only where it is mixed into hydrogen, as the table above requires.
* EOS: ions $`X + Y/4 + Z_0/\langle A\rangle + X_{\mathrm C}/12 + X_{\mathrm O}/16`$ and electrons
  $`X + Y/2 + Z_0\langle Z/A\rangle + (X_{\mathrm C} + X_{\mathrm O})/2`$ per $`m_u`$; the ideal EOS
  $`1/\mu = 2X + 3Y/4 + Z_0/2 + 7X_{\mathrm C}/12 + 9X_{\mathrm O}/16`$. Conduction uses
  $`\langle Z^2\rangle^{1/2}`$ of the same list.
* Opacity: tables and formulas are given $`Z_0 + X_{\mathrm C} + X_{\mathrm O}`$ (what they were
  given before plan 9.5, when carbon counted in $`Z`$); C/O-enhanced tables are still not used.
* Composition step: after mixing (all five fractions), H → He, then He → C by $`\varepsilon_{3\alpha}/Q_{3\alpha}`$,
  then C → O by $`\varepsilon_{\mathrm C\alpha}/Q_{\mathrm C\alpha}`$, each limited to what the point
  (or mixed zone) contains; a gram of carbon takes ⅓ g of helium. Time step: $`\Delta X_{\mathrm C}`$
  from ¹²C(α,γ) is limited like $`\Delta X`$ and $`\Delta Y`$.
* Rate multiplier $`f_{\mathrm{C}\alpha}`$: `AnalyticNuclear::f_c_alpha`, default 1 (CF88).
* Fuel accounting per channel (`EnergyCheck::fuel_by_channel`): carbon made by triple-α,
  $`\Delta\!\int\! X_{\mathrm C}dm + \tfrac34\Delta\!\int\! X_{\mathrm O}dm`$, times $`Q_{3\alpha}`$, and carbon
  burnt, $`\tfrac34\Delta\!\int\! X_{\mathrm O}dm`$, times $`Q_{\mathrm C\alpha}`$, each against
  $`\int\!\!\int\varepsilon\,dm\,dt`$ of the channel; `fuel_ratio` now adds all channels. The
  helium channels integrate the rates of the previous model, those the composition step uses:
  $`\varepsilon_{\mathrm C\alpha}`$ grows with $`X_{\mathrm C}`$ by up to 10% per step, and the rates of
  the new model left a lag of 5%.

*Test: a 1 M☉ helium star* (`tests/helium_burning.rs`, analytic physics, 2026-11-21)
**[measured]**. A ZAMS model with $`X = 0`$ does not converge from the default starting guess; it
is reached by continuation in $`X`$ from the solar composition (half-decades to $`10^{-9}`$, then
decades). Below $`X \approx 10^{-10}`$ the CNO cycle can no longer supply the luminosity and the core
heats to helium ignition: $`T_c = 1.34\times10^8`$ K, $`L = 222\,L_\odot`$, $`R = 0.20\,R_\odot`$,
all of it from triple-α. Evolved (with remeshing) until $`Y_c = 0.1`$: 89 models, 9.4 Myr,
central $`X_{\mathrm C} = 0.69`$, $`X_{\mathrm O} = 0.19`$ (C/O = 3.6), $`L = 324\,L_\odot`$;
¹²C(α,γ) supplies 52% of the nuclear luminosity at the end (triple-α dominates at
$`Y_c = 0.75`$). Mass fractions sum to 1 to $`10^{-12}`$ at every point; $`Z_0`$ is unchanged.
Fuel accounting: triple-α 1.0003, ¹²C(α,γ) 1.0003, all channels against $`\int L\,dt`$ 1.0012;
total-energy residual $`-6.6\times10^{-3}`$. The C/O ratio of the analytic physics is not
validated (no neutrino losses, no semiconvection or overshoot at the core boundary, CF88 rate).

### Screening diagnostic (plan 1.11)

Every model reports its largest screening exponent $`H`$ (since plan 9.7 of the all-regime
prescription, with the EOS's electrons), over the channels that supply at least 1%
of the nuclear luminosity and the points where such a channel burns at more than 1% of its peak
rate; an evolution run also keeps the largest value of any of its models. For the pp chain the
exponent of ³He + ³He ($`Z_1 Z_2 = 4`$, through $`\alpha`$) is the largest. The app and the CLI
warn

* above 0.3: outside the weak regime (intermediate or strong screening; helium ignition is not
  validated, see the [gate](#helium-ignition-gate));
* at a cap, for rates that cap the exponent: the rates are unreliable. The all-regime screening
  of plan 9.7 has no cap.

**[measured]** ZAMS Sun (plan 9.7, 2026-11-21): 0.178 (³He + ³He, largest at $`m/M = 0.54`$;
before plan 9.7, with the weak formula, 0.23). 0.5 M☉ ZAMS: 0.33 (³He + ³He; before 0.45), so cool
dense cores of low-mass stars still get the "outside the weak regime" warning; 3, 5, 10, 20 M☉:
0.080, 0.051, 0.029, 0.018 (p + ¹⁴N; before 0.12, 0.07, 0.04, 0.02).

## Planned

### Kippenhahn–Weigert–Weiss (2012) fits (not adopted)

The pp and ¹⁴N(p,γ) rates are now built from Solar Fusion III (implemented, plan 1.6, see
[Rates](#rates)). The plan asked for an embedded table (`include_bytes!` or a `const` array); the
code instead tabulates the quadrature once at first use (`OnceLock`), which needs no file access
either (rule 1) and keeps the coefficients of Table I in the source.

Why not the interim fits: the Kippenhahn, Weigert & Weiss (2012) coefficients below could not be
checked (book and NACRE not available on 2026-11-11 or 2026-11-12), and the plan already required
them to be compared with rates built from Solar Fusion III before adoption. Building those rates
directly removes the step. The coefficients are kept below for reference only. **They are not
adopted and must not be coded.**

*Status of the check (2026-11-11).* The book was not available (no legitimate copy found), so the
section number and every coefficient below are still unverified. The forms (polynomial in $`T_9`$,
the $`\exp[-(T_9/0.8)^2]`$ cut-off) look like the NACRE analytic approximations (Angulo et al. 1999,
Nucl. Phys. A 656, 3–183), but that paper was not obtained either (publisher access only), so the
attribution is a guess to be checked.

For the comparison asked for above: Solar Fusion III gives $`S_{1,14}(0) = 1.68\times10^{-3}`$ MeV b
(Table I). With constant $`S`$ this gives a ¹⁴N(p,γ)¹⁵O rate 0.52–0.56 times the Fowler–Caughlan
rates (1975, behind the implemented CNO fit, and the 1988 reproduction) at $`T_6 = 15`$–30 (rough: $`S'`$, $`S''`$ and
resonances neglected; computed 2026-11-11). The implemented CNO rate is therefore probably about a
factor 2 too high at main-sequence temperatures. The rates built from Solar Fusion III (above)
contain the reduction by construction; the quadrature replaced this rough estimate in plan 1.6
(measured ratio 0.520–0.546 at $`T_6 = 10`$–30, see [Rates](#rates)).

```math
\varepsilon_{pp} = 2.57\times10^{4}\;\psi\, f_{11}\, g_{11}\;\rho\, X^2\; T_9^{-2/3} \exp\!\left(-3.381\, T_9^{-1/3}\right)
```

```math
g_{11} = 1 + 3.82\,T_9 + 1.51\,T_9^2 + 0.144\,T_9^3 - 0.0114\,T_9^4
```

```math
\varepsilon_{\mathrm{CNO}} = 8.24\times10^{25}\; g_{14,1}\; X_{\mathrm{CNO}}\, X\, \rho\; T_9^{-2/3} \exp\!\left(-15.231\, T_9^{-1/3} - \left(\frac{T_9}{0.8}\right)^{2}\right)
```

```math
g_{14,1} = 1 - 2.00\,T_9 + 3.41\,T_9^2 - 2.43\,T_9^3
```

### Burning out of equilibrium (F4)

Separate abundances for ³He, ¹²C and ¹⁴N would follow the approach to equilibrium on the ZAMS. Not
scheduled; needs the composition vector of D18.

---
Previous: [Opacity](03-opacity.md) · [Contents](README.md) · Next: [Convection and mixing](05-convection-and-mixing.md)
