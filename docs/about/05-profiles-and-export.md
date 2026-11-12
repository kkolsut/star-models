# 5. Profile quantities & export

Every model provides these quantities at every mesh point, from the centre to the photosphere. The
**Profile CSV** export writes one column per quantity, in this order:

| # | Column | Quantity |
|---|---|---|
| 1 | `q` | mass fraction $`m/M`$ |
| 2 | `m_msun` | enclosed mass $`m`$ [M☉] |
| 3 | `r_rsun` | radius $`r`$ [R☉] |
| 4 | `r_cm` | radius $`r`$ [cm] |
| 5 | `P_dyn_cm2` | total pressure $`P = P_{\mathrm{gas}} + P_{\mathrm{rad}}`$ [dyn cm⁻²] |
| 6 | `T_K` | temperature $`T`$ [K] |
| 7 | `rho_g_cm3` | density $`\rho`$ [g cm⁻³] |
| 8 | `L_lsun` | luminosity $`L(r)`$ flowing through the sphere of radius $`r`$ [L☉] |
| 9 | `kappa_cm2_g` | Rosseland mean opacity $`\kappa`$ [cm² g⁻¹] |
| 10 | `eps_pp` | energy generation by the pp chain [erg g⁻¹ s⁻¹] |
| 11 | `eps_cno` | energy generation by the CNO cycle [erg g⁻¹ s⁻¹] |
| 12 | `eps_3a` | energy generation by triple-α helium burning [erg g⁻¹ s⁻¹] |
| 13 | `eps_total` | total nuclear energy generation $`\varepsilon`$ [erg g⁻¹ s⁻¹] |
| 14 | `eps_grav` | gravitational energy release $`\varepsilon_{\mathrm{grav}} = -T\,ds/dt`$ [erg g⁻¹ s⁻¹] (0 on the ZAMS) |
| 15 | `nabla_rad` | radiative gradient $`\nabla_{\mathrm{rad}}`$: the $`d\ln T/d\ln P`$ needed to carry $`L`$ by radiation alone |
| 16 | `nabla_ad` | adiabatic gradient $`\nabla_{\mathrm{ad}} = (\partial\ln T/\partial\ln P)_s`$, at constant entropy |
| 17 | `nabla` | actual gradient $`\nabla = d\ln T/d\ln P`$ used in the model |
| 18 | `convective` | 1 if convective by the Schwarzschild criterion ($`\nabla_{\mathrm{rad}} > \nabla_{\mathrm{ad}}`$), else 0 |
| 19 | `v_conv_cm_s` | convective velocity from mixing-length theory [cm s⁻¹] (0 in radiative layers) |
| 20 | `F_conv_frac` | fraction of the energy flux carried by convection |
| 21 | `H_P_cm` | pressure scale height $`H_P = P/(\rho g)`$ [cm] |
| 22 | `mu` | mean molecular weight $`\mu`$ |
| 23 | `beta` | gas-pressure fraction $`\beta = P_{\mathrm{gas}}/P`$ |
| 24 | `Gamma1` | first adiabatic exponent $`\Gamma_1`$ (used for pulsations) |
| 25 | `X` | hydrogen mass fraction (uniform in ZAMS models, depleted by burning during evolution) |
| 26 | `Y` | helium mass fraction |
| 27 | `Z` | metal mass fraction |
| 28 | `c_s_cm_s` | adiabatic sound speed $`c_s = \sqrt{\Gamma_1 P/\rho}`$ [cm s⁻¹] |
| 29 | `N2_s-2` | buoyancy frequency squared $`N^2`$ [s⁻²]; negative where the layer is Ledoux-unstable |

## Summary JSON

| Key | Contents |
|---|---|
| `params` | the input parameters: `mass`, `x`, `z`, `alpha_mlt`, `f_cno`, `n_points`, `alpha_ov` |
| `physics` | the description of the physics used |
| `globals` | every quantity of [Results: global quantities](04-results.md), plus `newton_iterations`, `max_correction` and `residuals` |
| `convective_zones` | one entry per zone, with `q_inner`, `q_outer`, `r_inner` and `r_outer` ($`r`$ in units of $`R`$) |

> All values are in CGS units unless the column name says otherwise. Numbers are written with 10
> significant digits.

---
← Previous: [Results: global quantities](04-results.md) · [Contents](README.md) · Next: [Structure equations](06-structure-equations.md) →
