# 4. Nuclear burning

Energy generation rates $`\varepsilon(\rho, T)`$ in erg g⁻¹ s⁻¹, and the composition changes they
cause. In the code the rates are the `Nuclear` trait.

## Implemented

Source: `crates/stellar-core/src/physics/nuclear.rs` (rates),
`crates/stellar-core/src/structure/evolve.rs` (composition).

The rate fits are those of Kippenhahn & Weigert (1990).

### Rates

*pp chain.* The pp I branch only; the branching factor $`\psi`$ is taken as 1.

```math
\varepsilon_{pp} = 2.38\times10^{6}\;\rho\, X^2\, f_{1}\, g_{11}\; T_6^{-2/3}\, \exp\!\left(-33.80\, T_6^{-1/3}\right)
```

```math
g_{11} = 1 + 0.0123\, T_6^{1/3} + 0.0109\, T_6^{2/3} + 0.0009\, T_6
```

> The About page and `equations.json` call $`g_{11}`$ "$`\psi`$" (finding B1). The name used here is
> the correct one.

*CNO cycle.* The slowest step, p + ¹⁴N, sets the rate.

```math
\varepsilon_{\mathrm{CNO}} = 8.67\times10^{27}\;\rho\, X\, X_{\mathrm{CNO}}\, f_{7}\, g_{14,1}\; T_6^{-2/3}\, \exp\!\left(-152.28\, T_6^{-1/3}\right)
```

```math
g_{14,1} = 1 + 0.0027\, T_6^{1/3} - 0.00778\, T_6^{2/3} - 0.000149\, T_6, \qquad X_{\mathrm{CNO}} = f_{\mathrm{CNO}}\, Z
```

*Triple-α*

```math
\varepsilon_{3\alpha} = 5.09\times10^{11}\;\rho^2\, Y^3\, f_{4}\, f_{8}\; T_8^{-3}\, \exp\!\left(-\frac{44.027}{T_8}\right)
```

### Electron screening

Weak screening (Salpeter), with the exponent capped at 2:

```math
f_{Z_1 Z_2} = \exp\!\left(0.188\, Z_1 Z_2\, \zeta\, \rho^{1/2}\, T_6^{-3/2}\right), \qquad \zeta^2 = 2X + 1.5\,Y + 4Z
```

| Reaction | $`Z_1 Z_2`$ |
|---|---|
| p + p | 1 |
| p + ¹⁴N | 7 |
| α + α | 4 |
| α + ⁸Be | 8 |

### Composition changes

```math
\frac{\partial X}{\partial t} = -\frac{\varepsilon_{pp} + \varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{H}}}, \qquad Q_{\mathrm{H}} = 6.27\times10^{18}\ \mathrm{erg\,g^{-1}}
```

```math
\frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}}, \qquad \frac{\partial Z}{\partial t} = \frac{\varepsilon_{3\alpha}}{Q_{\mathrm{He}}}, \qquad Q_{\mathrm{He}} = 5.85\times10^{17}\ \mathrm{erg\,g^{-1}}
```

* $`Q_{\mathrm{H}}`$ is 26.2 MeV per helium nucleus, the pp-chain value after neutrino losses. It is
  used for the CNO cycle as well.
* $`Q_{\mathrm{He}}`$ is 7.275 MeV per carbon nucleus.
* Carbon from helium burning is counted in $`Z`$. The catalyst abundance $`X_{\mathrm{CNO}}`$ keeps
  its initial value.
* The fits assume that ³He and the CNO isotopes are in equilibrium.

## Planned

### pp branching factor (D11, plan 1.6)

Let $`r_{11}`$, $`r_{33}`$, $`r_{34}`$ be the rates of p + p, ³He + ³He and ³He + ⁴He. With ³He in
equilibrium, $`r_{11} = 2 r_{33} + r_{34}`$. Helium production relative to the pp I branch alone:

```math
\phi = \frac{2\,(r_{33} + r_{34})}{r_{11}} = 1 - \alpha + \sqrt{\alpha\,(\alpha + 2)}, \qquad \alpha = C_\alpha \left(\frac{Y}{4X}\right)^{2} \exp\!\left(-100\,T_6^{-1/3}\right)
```

$`\phi`$ runs from 1 to 2. The fraction of helium made through pp II and the energy factor (pp II
gives 25.66 MeV instead of 26.20 MeV to the gas; pp III is neglected, safe below
$`T_6 \approx 20`$):

```math
F_{\mathrm{II}} = \frac{2\,(\phi - 1)}{\phi}, \qquad \psi = \phi\,\left(1 - 0.021\,F_{\mathrm{II}}\right), \qquad \varepsilon_{pp} = \psi\,\varepsilon_{pp,\mathrm{I}}
```

*Fuel consumption must change with it.* $`\phi`$ counts helium nuclei, $`\psi`$ counts energy, so the
energy delivered per gram of hydrogen burnt depends on the branch:

```math
\frac{\partial X}{\partial t}\bigg|_{pp} = -\frac{\varepsilon_{pp}}{Q_{pp}^{\mathrm{eff}}}, \qquad Q_{pp}^{\mathrm{eff}} = Q_{pp,\mathrm{I}}\,\frac{\psi}{\phi} = Q_{pp,\mathrm{I}}\left(1 - 0.021\,F_{\mathrm{II}}\right), \qquad Q_{pp,\mathrm{I}} = 6.27\times10^{18}\ \mathrm{erg\,g^{-1}}
```

Here $`\varepsilon_{pp,\mathrm{I}}`$ is the baseline rate (pp I only) and $`\varepsilon_{pp}`$ the
corrected one. The baseline must use the same rate convention as $`\alpha`$ below. In a fully
consistent calculation it is built from the same $`\langle\sigma v\rangle_{11}`$, with $`E_{\mathrm{I}}`$
the energy deposited per completed pp I chain:

```math
r_{11} = \tfrac{1}{2}\,n_p^{2}\,\langle\sigma v\rangle_{11}, \qquad \varepsilon_{pp,\mathrm{I}} = \frac{r_{11}}{2\rho}\,E_{\mathrm{I}}
```

Using the implemented (1990) fit for $`\varepsilon_{pp,\mathrm{I}}`$ with a modern $`\alpha`$ is a
hybrid; it is allowed as an interim step if it is called that and the difference between the two
p + p rates is stated (about 8%, see [Validity](11-validity-and-roadmap.md)).

For large $`\alpha`$ the expression for $`\phi`$ subtracts nearly equal numbers; evaluate the
equivalent form

```math
\phi = 1 + \frac{2}{\sqrt{1 + 2/\alpha} + 1}
``` Most of the extra helium production is already carried by
$`\psi`$ in $`\varepsilon_{pp}`$. Keeping the fixed $`Q_{pp,\mathrm{I}}`$ with the corrected rate
would therefore underestimate hydrogen consumption by at most 2.1% ($`F_{\mathrm{II}} = 1`$),
within the pp I/II approximation.

*How $`\alpha`$ is obtained.* The specification is the definition, evaluated with consistently
defined thermal rates $`\langle\sigma v\rangle`$ from one evaluation (Solar Fusion III), including
their energy dependence and screening:

```math
\alpha = \left(\frac{n_4}{n_p}\right)^{2} \frac{\langle\sigma v\rangle_{34}^{\,2}}{\langle\sigma v\rangle_{11}\,\langle\sigma v\rangle_{33}}, \qquad \frac{n_4}{n_p} = \frac{Y}{4X}
```

The closed form with $`C_\alpha`$ above is its leading-order approximation (constant $`S`$ factors,
no screening, leading Gamow term), with reduced masses $`\mu_{ij}`$ in atomic mass units:

```math
C_\alpha \simeq \frac{S_{34}^{2}}{S_{11}\,S_{33}}\left(\frac{4\,\mu_{11}\,\mu_{33}}{\mu_{34}^{2}}\right)^{1/3}
```

With the Solar Fusion III values (Table 1, in MeV b) $`S_{11}(0) = 4.09\times10^{-25}`$,
$`S_{33}(0) = 5.21`$, $`S_{34}(0) = 5.61\times10^{-4}`$ and $`\mu_{11} = 1/2`$,
$`\mu_{33} = 3/2`$, $`\mu_{34} = 12/7`$, this gives $`C_\alpha \simeq 1.49\times10^{17}`$ (computed
by the reviewer of these pages and re-computed here; the exponent 100 follows from the same Gamow
terms). It is an estimate for checking the implementation, not a precision fit.

The value $`1.2\times10^{17}`$ quoted earlier from memory of Clayton (1968) is of the same size;
its attribution is unverified **[check source]**. The earlier expectation that a modern value would
be "several times smaller" was wrong and is withdrawn.

### Newer rate fits (D11) **[check source]**

Kippenhahn, Weigert & Weiss (2012, §18.5). These are an interim step, not the final reference: the
coefficients are to be compared with rates built from the Solar Fusion III $`S`$ factors before
they are adopted.

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

### Energy per gram by channel (D12, plan 1.5)

26.731 MeV per ⁴He is $`6.398\times10^{18}`$ erg/g. After neutrino losses:

```math
Q_{pp,\mathrm{I}} = 6.27\times10^{18}, \qquad Q_{\mathrm{CNO}} = 5.99\times10^{18}\ \mathrm{erg\,g^{-1}}, \qquad \frac{\partial X}{\partial t} = -\frac{\varepsilon_{pp}}{Q_{pp}^{\mathrm{eff}}} - \frac{\varepsilon_{\mathrm{CNO}}}{Q_{\mathrm{CNO}}}
```

$`Q_{pp}^{\mathrm{eff}} = Q_{pp,\mathrm{I}}`$ as long as $`\psi = 1`$; with branching it is the
expression given above.

### Carbon and oxygen (D18, plan 9.5)

With $`Q_{3\alpha} = Q_{\mathrm{He}}`$ per gram of helium and
$`Q_{\mathrm{C}\alpha} = 7.162\ \mathrm{MeV}/12\,m_u = 5.76\times10^{17}`$ erg per gram of carbon
consumed by ¹²C(α,γ)¹⁶O:

```math
\frac{\partial X_{\mathrm{C}}}{\partial t} = \frac{\varepsilon_{3\alpha}}{Q_{3\alpha}} - \frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}, \qquad \frac{\partial X_{\mathrm{O}}}{\partial t} = \frac{4}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}
```

```math
\frac{\partial Y}{\partial t} = -\frac{\partial X}{\partial t} - \frac{\varepsilon_{3\alpha}}{Q_{3\alpha}} - \frac{1}{3}\,\frac{\varepsilon_{\mathrm{C}\alpha}}{Q_{\mathrm{C}\alpha}}
```

The rate $`\varepsilon_{\mathrm{C}\alpha}`$ is to be taken from Kippenhahn, Weigert & Weiss (2012)
or Caughlan & Fowler (1988).

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

### Screening beyond the weak regime (gate for helium ignition)

The implemented formula is weak (Debye–Hückel) screening, valid while its exponent is small. The
cap at 2 is a numerical safeguard, not a model of intermediate or strong screening. In the dense,
degenerate helium core of a red giant the triple-α reactions are outside the weak regime, and the
formula also ignores electron degeneracy in the screening cloud.

Consequence: the degenerate EOS makes a helium core appear, but the model must not be read as
predicting helium ignition (core mass, luminosity at the tip, the flash) until screening has a
prescription covering weak, intermediate and strong regimes. The prescription chosen is given in
the next section; the gate stays until it is implemented and tested.
Until then: report the largest screening exponent reached in a run, and warn when it exceeds about
0.3 (still approximately weak) and when it hits the cap.

### Screening prescription for all regimes (decided)

**Decided 2026-11-11** (recommendation of the Codex review, accepted). Two parts of it are not from
a published source and must be checked before they are coded: the electron-response extension for
the weak limit, and the product form for triple-α at intermediate coupling. Chosen over
Graboske et al. (1973) because one analytic interpolation is easier to keep smooth than joined
branches; Graboske's paper was not opened.

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

*Electron response in the weak limit* (an extension proposed by the review, **not** in the
published formula, and not verified). The published $`C_{ij}`$ assumes a rigid electron
background, so it does not reproduce the electron part of the implemented Salpeter formula. To
recover the correct weak limit for classical or degenerate electrons, replace $`C_{ij}`$ by
$`H_{\mathrm{DH}}/H_{\mathrm{LM,DH}}`$:

```math
D^{-2} = 4\pi e^{2}\left[\frac{\sum_s n_s Z_s^{2}}{kT} + \left(\frac{\partial n_e}{\partial \mu_e^{\mathrm{chem}}}\right)_{T}\right], \qquad H_{\mathrm{DH}} = \frac{Z_i Z_j\,e^{2}}{kT\,D}, \qquad H_{\mathrm{LM,DH}} = \frac{\Gamma_e^{3/2}}{\sqrt{3}}\left[Z_c^{5/2} - Z_i^{5/2} - Z_j^{5/2}\right]
```

$`\partial n_e/\partial\mu_e^{\mathrm{chem}}`$ comes from the Fermi–Dirac EOS. This fades out at
strong coupling and does not model electron polarisation there.

*Use.* $`f_{11}`$ for pp and $`f_{1,7}`$ for the CNO bottleneck, with the actual mixture
(including produced carbon once it is tracked). The exponent cap of 2 is removed; $`H`$ is kept as
a logarithm internally.

*Triple-α* (resonant). The static form, $`t_{ij} = 1`$, applied once to the unscreened rate:

```math
f_{3\alpha} \approx \exp\!\left[H_{2,2}^{\mathrm{static}} + H_{4,2}^{\mathrm{static}}\right]
```

The first term shifts the equilibrium abundance of ⁸Be, the second the capture of the third α
through the Hoyle state. The basis is the difference of Coulomb chemical potentials,
$`\ln f_{3\alpha} = (3\mu_\alpha^{C} - \mu_{\mathrm{C}}^{C})/kT`$ (Kushnir, Waxman & Chugunov 2019,
eqs. 1, 7–8 **[compared]**), valid for classical ions; electron degeneracy is allowed. In the
strong limit the two stages telescope exactly to $`3f_0(\Gamma_2) - f_0(\Gamma_6)`$. At
intermediate coupling the product of two interpolation factors is an approximation. The
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

### Burning out of equilibrium (F4)

Separate abundances for ³He, ¹²C and ¹⁴N would follow the approach to equilibrium on the ZAMS. Not
scheduled; needs the composition vector of D18.

---
Previous: [Opacity](03-opacity.md) · [Contents](README.md) · Next: [Convection and mixing](05-convection-and-mixing.md)
