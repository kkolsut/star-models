# 3. Opacity

The opacity is the Rosseland mean $`\kappa(\rho, T)`$ in cm² g⁻¹. In the code it is the `Opacity`
trait.

## Implemented: analytic formulas

Source: `crates/stellar-core/src/physics/opacity.rs`. Coefficients follow Carroll & Ostlie (2017),
converted from SI to CGS.

*Electron scattering*

```math
\kappa_{\mathrm{es}} = 0.2\,(1 + X)
```

*Kramers free–free*

```math
\kappa_{\mathrm{ff}} = 3.68\times10^{22}\,(1 - Z)(1 + X)\,\rho\,T^{-3.5}
```

*Kramers bound–free*, with the guillotine factor $`t`$

```math
\kappa_{\mathrm{bf}} = \frac{4.34\times10^{25}}{t}\,Z\,(1 + X)\,\rho\,T^{-3.5}, \qquad t = 2.82\,\bigl[\rho\,(1 + X)\bigr]^{0.2}
```

*H⁻ ion* (cool envelopes). For $`Z < 10^{-4}`$ the formula uses $`Z = 10^{-4}`$.

```math
\kappa_{\mathrm{H}^-} = 2.5\times10^{-31}\,\frac{Z}{0.02}\,\rho^{1/2}\,T^{9}
```

*Molecular floor*

```math
\kappa_{\mathrm{mol}} = 0.1\,Z
```

*Combined opacity*

```math
\kappa = \kappa_{\mathrm{mol}} + \left[\frac{1}{\kappa_{\mathrm{H}^-}} + \frac{1}{\kappa_{\mathrm{es}} + \kappa_{\mathrm{ff}} + \kappa_{\mathrm{bf}}}\right]^{-1}
```

The harmonic combination lets H⁻ dominate only in cool layers; without it the $`T^9`$ term would run
away at high temperature.

### Accuracy

* At the centre of the ZAMS Sun ($`T = 1.26\times10^{7}`$ K, $`\rho = 64`$ g cm⁻³) the formulas give
  $`\kappa = 2.73`$ **[measured]**. Detailed solar models have roughly half that at
  their own centre (literature value). This is a comparison between centres of different models,
  not at the same state: the test still to be made is the analytic formula against a table at
  identical $`(\rho, T, X, Z)`$. The high opacity is the leading suspect for the model Sun's
  0.42 L☉, but its share has not been isolated (see [Validity](11-validity-and-roadmap.md)).
* Between $`10^5`$ and $`10^6`$ K the formulas are too low, which makes the solar convective
  envelope too shallow.
* There are no opacity peaks from metal lines (for example the iron peak near $`2\times10^5`$ K).
* Below 2600 K the formulas are outside their range; evolution runs stop there.

## Planned

### Tabulated opacity (D14, plan phase 7)

Tables (OPAL for the interior, a low-temperature set for the envelope) are given on a grid of
$`\log T`$ and $`\log R`$ for a set of $`(X, Z)`$:

```math
R = \frac{\rho}{T_6^{3}}, \qquad T_6 = 10^{-6}\,T
```

This definition agrees with Iglesias & Rogers (1996, §2.2, p. 944), the header of the OPAL file
`GN93hz`, and Ferguson et al. (2005, §3, p. 587). Source check (task 0.3, 2026-11-11) of the two
sets named in task 7.1. *Decided 2026-11-12:* OPAL is the interior set and OP is dropped from the
plan. OPAL is the set that was checked, it alone provides the C/O-enhanced type 2 tables needed
after helium burning (below), and its grid and the Ferguson grid overlap where the blend is
planned. OP can be added later as an alternative behind the same `Opacity` trait; it would need its
own source check then:

| | OPAL (Iglesias & Rogers 1996) | Low temperature (Ferguson et al. 2005) |
|---|---|---|
| $`\log T`$ | 3.75 to 8.70, 70 values (step 0.05, coarser at the top) | 2.70 to 4.50; step 0.05, 0.01 between 2.9 and 3.5 |
| $`\log R`$ | −8.0 to +1.0, 19 values (step 0.5) | −8.0 to +1.0, 19 values (step 0.5) |
| compositions per set | 126 $`(X, Z)`$ tables | 120 in the paper; the distributed `f05.g93` archive has 155 files and 85 temperatures |
| missing entries | yes: tables are not rectangular. Blank at the highest densities and temperatures (degenerate, conduction-dominated, cut near $`\eta = 8`$); 9.999 at low $`T`$ and low $`R`$ in some hydrogen-free tables. In the file, the $`X = 0.7`$, $`Z = 0.02`$ rows checked ($`\log T`$ = 3.75, 7.00, 8.70) agree with Table 3 of the paper (pp. 945–946) | none in the 155 files of `f05.g93` |
| file layout | ASCII, one 77-line block per composition after a 240-line header | ASCII, one file per composition; fixed-width columns of 7 characters, and at $`Z = 0`$ values below −10 touch their neighbours, so the reader must split by column, not by white space |
| where stated | paper §2.2, p. 944; `GN93hz` header, lines 1–60 | paper §6, p. 595; the files themselves |

The two sets overlap in $`3.75 \le \log T \le 4.5`$, so the blending range below lies inside both.

Requirements on the implementation:

* interpolate $`\log\kappa`$ with continuous first derivatives, because the Jacobian is built by
  finite differences;
* blend the high- and low-temperature tables over an overlap, for example
  $`3.75 \le \log T \le 4.0`$;
* fall back smoothly to the analytic formulas outside the tables, including the missing
  high-density corner of the OPAL tables (where conduction, below, carries the heat);
* keep the analytic formulas selectable.

### Heavy-element mixture

Tables are computed for one fixed mixture of heavy elements, scaled by $`Z`$. The same mixture must
define $`f_{\mathrm{CNO}}`$ in the burning rates and $`\langle A\rangle`$ in the EOS.

Which mixtures exist ready-made (checked 2026-11-11; this corrects the earlier example "Asplund et
al. 2009", which no ready-made OPAL table uses):

* OPAL: the precomputed solar Type 1 set and all Type 2 sets use Grevesse & Noels (1993)
  (Iglesias & Rogers 1996, Table 1, p. 944; OPAL "existing tables" page; `GN93hz` header). Other
  mixtures, Asplund et al. (2009) among them, need tables generated with the OPAL web form
  (Iglesias & Rogers 1996, §5, p. 953).
* Ferguson et al. (2005): base mixture Grevesse & Noels (1993) (§6, p. 595). The Wichita State page
  now also offers sets for Grevesse & Sauval 1998, Asplund et al. 2004, 2009 and 2021, Caffau et al.
  2011 and others.

So the only matching pair available without generating tables is Grevesse & Noels (1993). The app's
solar composition preset is Asplund et al. (2009); choosing between the two is part of task 7.1.

A scaled-solar $`(X, Z)`$ table is wrong for helium-burning material, where the "metals" are almost
pure carbon and oxygen. If C and O are tracked (D18), the opacity needs tables with excess C and O
as extra parameters: the OPAL "Type 2" tables. They are tabulated in $`\Delta X_{\mathrm{C}}`$ and
$`\Delta X_{\mathrm{O}}`$, the C and O mass fractions beyond those contained in $`Z`$, at
$`X = 0, 0.03, 0.1, 0.35, 0.7`$ (none above 0.7) and $`Z = 0, 0.001, 0.004, 0.01, 0.02, 0.03,
0.05, 0.1`$, with the Grevesse & Noels (1993) base mixture (OPAL "existing tables" page; method in
Iglesias & Rogers 1993, ApJ 412, 752; the 1996 paper, p. 943 and §5, p. 953, lists C and O as
interpolation variables). Until then, results after helium ignition stay qualitative, whatever else
is upgraded.

### Electron conduction (with D15b)

In degenerate matter heat is also carried by electrons. The two channels add as conductances:

```math
\frac{1}{\kappa} = \frac{1}{\kappa_{\mathrm{rad}}} + \frac{1}{\kappa_{\mathrm{cond}}}
```

This agrees with Cassisi et al. (2007, §2, the sentence before eq. 1; read in the arXiv version
astro-ph/0703011), where $`\kappa_{\mathrm{cond}} = 16\sigma T^3/(3\rho K)`$ (eq. 1) with
$`K`$ the electron thermal conductivity.

*Decided 2026-11-12:* $`K`$, and with it $`\kappa_{\mathrm{cond}}(\rho, T)`$, comes from the table of
thermal conductivities of A. Potekhin at `www.ioffe.ru/astro/conduct/` (page read 2026-11-12):
pure, non-magnetic electron–ion plasmas, with a Fortran interpolation tool, last updated
2021-05-18. It implements Cassisi et al. (2007) (ρ from 10⁻⁶ to 10⁹ g cm⁻³, T from 10³ to 10⁹ K,
ion charge 1 to 60, per §2.4 of that paper) with the later corrections listed on the page,
including the bridge between moderate and strong degeneracy of Cassisi et al. (2021, A&A 654,
A149). It is the same source MESA uses. The site asks users to cite Potekhin, Pons & Page (2015,
Space Sci. Rev. 191, 239) and states no licence, so embedding it falls under the same permission
question as the radiative tables (plan 7.1). How the pure-plasma values are combined for a mixture
is taken from Cassisi et al. (2007, §2) when 7.5 is coded; that paper is available (arXiv
astro-ph/0703011).

---
Previous: [Equation of state](02-equation-of-state.md) · [Contents](README.md) · Next: [Nuclear burning](04-nuclear-burning.md)
