# 16. Constants & units

All calculations use CGS units internally. Solar units are used for display. The values are those in
`crates/stellar-core/src/constants.rs`.

| Constant | Symbol | Value | Unit |
|---|---|---|---|
| Gravitational constant | $`G`$ | $`6.67430 \times 10^{-8}`$ | cm³ g⁻¹ s⁻² |
| Speed of light | $`c`$ | $`2.99792458 \times 10^{10}`$ | cm s⁻¹ |
| Stefan–Boltzmann constant | $`\sigma`$ | $`5.670374 \times 10^{-5}`$ | erg cm⁻² s⁻¹ K⁻⁴ |
| Radiation constant | $`a = 4\sigma/c`$ | $`7.56573 \times 10^{-15}`$ | erg cm⁻³ K⁻⁴ |
| Boltzmann constant | $`k`$ | $`1.380649 \times 10^{-16}`$ | erg K⁻¹ |
| Atomic mass unit | $`m_u`$ | $`1.660539 \times 10^{-24}`$ | g |
| Solar mass | $`M_\odot`$ | $`1.98847 \times 10^{33}`$ | g |
| Solar radius (IAU 2015) | $`R_\odot`$ | $`6.957 \times 10^{10}`$ | cm |
| Solar luminosity (IAU 2015) | $`L_\odot`$ | $`3.828 \times 10^{33}`$ | erg s⁻¹ |
| Solar effective temperature (IAU 2015) | $`T_{\mathrm{eff},\odot}`$ | 5772 | K |

Units used in the outputs:

| Quantity | Unit | SI equivalent |
|---|---|---|
| Pressure | dyn cm⁻² | 1 dyn cm⁻² = 0.1 Pa |
| Density | g cm⁻³ | 1 g cm⁻³ = 1000 kg m⁻³ |
| Energy generation | erg g⁻¹ s⁻¹ | 1 erg g⁻¹ s⁻¹ = 10⁻⁴ W kg⁻¹ |

---
← Previous: [Accuracy & limitations](15-limitations.md) · [Contents](README.md) · Next: [Command-line tool](17-cli.md) →
