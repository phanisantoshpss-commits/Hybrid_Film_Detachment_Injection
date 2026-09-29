# HybridFilmDetachmentInjection

A custom runtime-selectable OpenFOAM surface-film `injectionModel` for combined **gravity driven dripping** and **shear driven film stripping**.

The model was developed as part of the **OpenFOAM Programming** project to study detachment of a liquid film and its conversion into Lagrangian droplets.

**Author:** Satya Sai Phani Santosh Pilaka  
**OpenFOAM version:** v2306

---

## Model concept

Film detachment is allowed through two mechanisms:

1. **Gravity driven dripping**, described using a Bond number.
2. **Shear stripping**, described using a film Weber number.

The model calculates one total detached mass. The relative gravity and shear contributions are used to determine the characteristic droplet diameter.

---

## Gravity-driven detachment

The local Bond number is

```math
Bo =
\frac{(\rho_l-\rho_g)\,g_n\,\delta^2}
{\sigma}
```

where

```math
g_n = \mathbf{g}\cdot\mathbf{n}
```

is the gravity component normal to the film surface.

Gravity detachment is active when

```math
Bo > Bo_{\mathrm{crit}}
```

with the exceedance factor

```math
F_{Bo}
=
\max
\left(
\frac{Bo}{Bo_{\mathrm{crit}}}-1,\,
0
\right).
```

---

## Gas-shear stripping

The relative gas-film velocity is

```math
\mathbf{U}_{rel}
=
\mathbf{U}_g-\mathbf{U}_f .
```

Only its tangential component is used:

```math
\mathbf{U}_{rel,t}
=
\mathbf{U}_{rel}
-
(\mathbf{U}_{rel}\cdot\mathbf{n})\mathbf{n}.
```

The film Weber number is

```math
We_f
=
\frac{\rho_g
|\mathbf{U}_{rel,t}|^2
\delta}
{\sigma}.
```

Shear stripping is active when

```math
We_f > We_{\mathrm{crit}}
```

with

```math
F_{We}
=
\max
\left(
\frac{We_f}{We_{\mathrm{crit}}}-1,\,
0
\right).
```

---

## Detached mass

The total threshold exceedance is

```math
E = F_{Bo}+F_{We}.
```

A local capillary timescale is defined as

```math
\tau_\sigma
=
\sqrt{
\frac{\rho_l\delta^3}{\sigma}
}.
```

The detachment rate is

```math
\lambda
=
\frac{C_m E}{\tau_\sigma}.
```

The fraction of locally available film mass detached during one timestep is

```math
f_{\mathrm{det}}
=
1-\exp(-\lambda\Delta t).
```

Therefore,

```math
m_{\mathrm{detach}}
=
f_{\mathrm{det}}\,
m_{\mathrm{available}}.
```

The same mass is removed from the Eulerian film and transferred to the Lagrangian cloud.

---

## Droplet diameter

### Gravity only

For gravity-driven dripping, the capillary length is

```math
l_c
=
\sqrt{
\frac{\sigma}{\rho_l |\mathbf{g}|}
}.
```

The droplet diameter is

```math
d_{Bo}
=
C_d l_c
```

with the default value

```math
C_d = 3.3.
```

This follows the capillary-length approach used by OpenFOAM's `BrunDrippingInjection`.

### Shear only

For aerodynamic stripping, the Sauter mean diameter is estimated using the Sattelmayer prefilming-airblast correlation:

```math
d_{We}
=
0.67
\frac{\sigma^{0.75}}
{U_g^{1.57}}.
```

This is an empirical correlation and its use here represents a modelling approximation for wall-film stripping.

### Gravity and shear simultaneously

If both mechanisms are active,

```math
W_{Bo}
=
\frac{F_{Bo}}{F_{Bo}+F_{We}}
```

and

```math
W_{We}
=
\frac{F_{We}}{F_{Bo}+F_{We}}.
```

The hybrid diameter is calculated using an SMD-inspired weighted harmonic relation:

```math
d_{\mathrm{hybrid}}
=
\left(
\frac{W_{Bo}}{d_{Bo}}
+
\frac{W_{We}}{d_{We}}
\right)^{-1}.
```

This hybrid relation is a closure introduced in the present model.

---

## Runtime coefficients

Example `surfaceFilmProperties` configuration:

```cpp
injectionModels
(
    hybridFilmDetachmentInjection
);

hybridFilmDetachmentInjectionCoeffs
{
    bondNumCrit     1.0;
    weberNumCrit    5.0;
    cm              0.25;
    dCoeff          3.3;
}
```

---

## Source structure

```text
HybridFilmDetachmentInjection/
├── hybridFilmDetachmentInjection.H
├── hybridFilmDetachmentInjection.C
└── Make/
    ├── files
    └── options
```

The model is registered through the OpenFOAM runtime-selection mechanism and can therefore be selected directly from `surfaceFilmProperties`.

---

## Validation

The model is tested using four characteristic conditions:

| Case | Film condition | Gas-film condition | Expected behaviour |
|---|---|---|---|
| A | Thin film | Low relative velocity | No detachment |
| B | Thick film | Low relative velocity | Gravity dripping |
| C | Thin film | High relative velocity | Shear stripping |
| D | Thick film | High relative velocity | Hybrid detachment |

The planned parametric study evaluates gas velocity and monitors film thickness, detached mass, detachment location and generated droplet diameter.

---

## Limitations of the Validation Cases

- The initial validation cases are designed mainly to verify the model logic and branch activation rather than reproduce a specific experiment.
- Constant or simplified film and gas properties are assumed in the first tests.
- The critical values \(Bo_{\mathrm{crit}}\) and \(We_{\mathrm{crit}}\) are treated as model parameters and are not yet calibrated against experimental data.
- The shear-droplet diameter correlation is based on prefilming airblast atomization and may not be directly transferable to all wall-film geometries.
- Secondary breakup, droplet coalescence, evaporation, and detailed interface-wave effects are not considered during the first validation stage.
- The model is tested under high velocity regimes even though when considered as Laminar to clearly distinguish between the mechanisms involved and avoid high computational costs.

---

## Possible Improvements

- Validate and calibrate the model against experimental film-detachment data.
- Introduce viscosity effects through additional dimensionless groups such as the Ohnesorge number.
- Replace the simplified shear-diameter correlation with a correlation specifically developed for wall-film stripping.
- Add diagnostic output fields for \(Bo\), \(We_f\), \(F_{Bo}\), \(F_{We}\), inorder for better understanding of the mechanism involved.
- Extend the model to distinguish different detachment regimes and include more advanced droplet size distributions instead of a single characteristic diameter.
- Perform timestep and mesh sensitivity studies to verify numerical independence of the predicted detachment behaviour.

---

## References

1. Brun, P. T., Damiano, A., Rieu, P., Balestra, G., & Gallaire, F. (2015). *Rayleigh–Taylor instability under an inclined plane*. Physics of Fluids, 27, 084107.

2. Lefebvre, A. H. (1988). *Atomization and Sprays*. Hemisphere Publishing / CRC Press.

3. Sattelmayer, T., & Wittig, S. (1986). *Internal Flow Effects in Prefilming Airblast Atomizers: Mechanisms of Atomization and Droplet Spectra*. Journal of Engineering for Gas Turbines and Power, 108(3), 465–472.

4. OpenFOAM — `BrunDrippingInjection`, `kinematicSingleLayer`, and surface-film modelling framework.

---

## Status

This project is under development and is intended primarily as an OpenFOAM programming and surface-film modelling study.