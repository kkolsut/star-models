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

### Tabulated opacity (D14, plan phase 7) **[check source]**

Tables (OPAL or OP for the interior, a low-temperature set for the envelope) are given on a grid of
$`\log T`$ and $`\log R`$ for a set of $`(X, Z)`$:

```math
R = \frac{\rho}{T_6^{3}}
```

Requirements on the implementation:

* interpolate $`\log\kappa`$ with continuous first derivatives, because the Jacobian is built by
  finite differences;
* blend the high- and low-temperature tables over an overlap, for example
  $`3.75 \le \log T \le 4.0`$;
* fall back smoothly to the analytic formulas outside the tables;
* keep the analytic formulas selectable.

### Heavy-element mixture

Tables are computed for one fixed mixture of heavy elements (to be named with the table set, for
example Asplund et al. 2009), scaled by $`Z`$. The same mixture must define $`f_{\mathrm{CNO}}`$ in
the burning rates and $`\langle A\rangle`$ in the EOS.

A scaled-solar $`(X, Z)`$ table is wrong for helium-burning material, where the "metals" are almost
pure carbon and oxygen. If C and O are tracked (D18), the opacity needs tables with excess C and O
as extra parameters (the OPAL "type 2" form) **[check source]**. Until then, results after helium
ignition stay qualitative, whatever else is upgraded.

### Electron conduction (with D15b)

In degenerate matter heat is also carried by electrons. The two channels add as conductances:

```math
\frac{1}{\kappa} = \frac{1}{\kappa_{\mathrm{rad}}} + \frac{1}{\kappa_{\mathrm{cond}}}
```

$`\kappa_{\mathrm{cond}}(\rho, T)`$ from a published fit, to be chosen with the tables.

---
Previous: [Equation of state](02-equation-of-state.md) · [Contents](README.md) · Next: [Nuclear burning](04-nuclear-burning.md)
