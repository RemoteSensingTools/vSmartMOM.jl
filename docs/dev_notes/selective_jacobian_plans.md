# Retrieval-selected Jacobian plans

## Status and scope

This note documents the retrieval-selective analytic-linearization machinery
introduced for the synthetic OCO RRS/XCO2 inversions. The first concrete
retrieval flavor is `OCO_RRS_synth`; the underlying `JacobianPlan` and
`ActiveParameterLayout` types are retrieval-independent.

The feature changes which physical tangent columns are constructed and
propagated. It does **not** change the forward radiative-transfer equations,
the elemental/doubling/interaction derivative formulas, or the physical
meaning of vSmartMOM's native aerosol parameters.

Linearized Raman scattering remains unsupported. Both corrected and
uncorrected inversions therefore use the elastic `noRS` linearized forward
model. RRS enters only through the uncorrected truth measurement.

## Motivation

The historical linearized path allocated a native basis containing

```text
surface pressure
+ 7 parameters per aerosol
+ every layer of every variable gas (including the q-driven H2O block)
+ surface parameters
+ source parameters
```

For the 16-layer, three-aerosol OCO experiment this is much larger than the OE
state. More importantly, the derivative dimension multiplies the large
spectral phase matrices and all elemental, doubling, and interaction operator
workspaces. Removing unused columns only from the returned Jacobian would save
almost no computation or peak memory.

The selective path instead compiles a small band-local physical basis and
applies it while assembling the moment-invariant aerosol/gas cache. Therefore
inactive columns never enter the combined phase-Jacobian tensor or the MOM
operator workspaces.

## Public types and responsibilities

### `ParameterKey`

A `ParameterKey` is the stable identity of one physical parameter:

- `kind`: parameter family, such as `:atmosphere`, `:gas`, `:aerosol`,
  `:surface`, or `:source`;
- `component`: species/component identity, or `nothing`;
- `field`: physical field such as `:surface_pressure`, `:vmr`, `:tau_ref`,
  `:z0`, `:P0`, `:SIF760`, or `:mSIF`;
- `layer`: TOA-to-BOA layer index, with zero meaning not layer-resolved;
- `band`: spectral-band index, with zero meaning shared across bands.

Keys, rather than integer offsets, connect heterogeneous band-local
Jacobians to one global state.

### `ActiveParameterLayout`

One `ActiveParameterLayout` describes the actual trailing derivative axis of
one RT call. It stores:

- ordered local keys and human-readable names;
- the ordered global keys and names;
- `local_to_global`, used to scatter a local Jacobian into the global basis;
- `native_layer_columns`, which identify the historical physical-optics
  columns required to construct the active atmospheric basis;
- the number of atmospheric columns;
- contiguous surface, SIF, and canopy ranges.

Atmospheric columns must retain native component order: pressure, aerosol 1,
aerosol 2, ..., then gas. Retrieval ordering may differ because
`local_to_global` performs the final reordering. Native columns must be
positive, unique, and sorted. Surface/source blocks follow the atmospheric
block and must be contiguous because their analytic seed routines consume
unit ranges.

### `JacobianPlan`

A `JacobianPlan` owns the global named basis and one
`ActiveParameterLayout` per band. Retrieval-specific policy ends at this
compiled object. Optical-property assembly and the RT kernels consume the
generic plan without dispatching on retrieval type.

### `PlannedRTModelLin`

`PlannedRTModelLin` pairs the ordinary `RTModelLin` with a `JacobianPlan`.
Unknown properties are forwarded to the underlying model so existing
diagnostic code can continue reading fields such as `tau_dot_abs` and
`lin_aerosol_optics`.

## Construction and execution sequence

```text
model_from_parameters(OCO_RRS_synth(), params)
  |
  +-- build ordinary forward + linearized optical model
  |     +-- retain all forward absorption and aerosol optics
  |     +-- skip q-H2O tangent generation
  |     +-- skip ForwardDiff Mie/microphysics derivatives
  |          (zero structural tangents preserve internal type/shape contracts)
  |
  +-- jacobian_plan(OCO_RRS_synth(), ...)
        +-- identify CO2 layer centers below 10 km
        +-- create the global key/name basis
        +-- create each band-local layout and native-column map

rt_run_lin(model, planned_lin; i_band, sources)
  |
  +-- validate the band surface and SIF source column counts
  +-- allocate RT tangent operators with the local parameter count
  +-- split native atmospheric selection by component
  +-- cache only selected aerosol tau/SSA and gas tangent columns
  +-- if no selected aerosol Mie parameter exists:
  |     compute forward aerosol phase blocks only
  +-- attach compact phase tangents and mix optical properties
  +-- propagate the compact basis through elemental/doubling/interaction
  +-- return ObserverRTResultLin with result.layout = band layout

globalize_jacobian(K_local, result.layout)
  +-- allocate the global trailing dimension
  +-- initialize inactive band parameters to exact zero
  +-- scatter local columns through local_to_global
```

Forward aerosol optics and H2O absorption are never removed. Only their
derivatives are omitted when the retrieval declares them fixed.

## `OCO_RRS_synth` state mappings

The current 16-layer prior fixes layers 1--4 because their centers are above
10 km. CO2 layers 5--16 are active. The global 30-column physical precursor
is:

| Global columns | Parameters |
|---:|:---|
| 1 | surface pressure |
| 2:13 | CO2 VMR in layers 5:16, TOA-to-BOA |
| 14:16 | sulfate, organic-carbon, and UTLS-sulfate native `tau_ref` |
| 17:19 | sulfate, organic-carbon, and UTLS-sulfate physical `z0` |
| 20:22 | O2 A-band surface `P0`, `P1`, `P2` |
| 23:25 | weak-CO2 surface `P0`, `P1`, `P2` |
| 26:28 | strong-CO2 surface `P0`, `P1`, `P2` |
| 29:30 | `SIF760`, `mSIF` |

Band-local order follows the optical assembly contract, not the global OE
order:

| Local columns | O2 A band (12) | Weak/strong CO2 band (22) |
|---:|:---|:---|
| 1 | surface pressure | surface pressure |
| 2:7 | `(tau_ref,z0)` for each of 3 aerosols | same |
| 8:10 | band surface `P0:P2` | CO2 VMR layers 5:7 |
| 11:12 | `SIF760,mSIF` | CO2 VMR layers 8:9 |
| 13:19 | -- | CO2 VMR layers 10:16 |
| 20:22 | -- | band surface `P0:P2` |

`globalize_jacobian` maps these local positions to the global table. Thus the
O2 A-band CO2 entries and both CO2-band SIF entries are exact zeros without
ever being propagated.

The flavor currently validates all of the following:

- exactly three spectral bands;
- exactly three aerosol species;
- exactly three Legendre surface coefficients in every band;
- CO2 configured as a variable molecule;
- at least one CO2 layer center below 10 km.

## Native physical basis versus OE coordinates

The returned compact/globalized Jacobian is still expressed in vSmartMOM's
native physical coordinates. The inversion layer must apply the following
transformations.

### AOD reference wavelength

Let `tau_ref` be optical depth at the scattering configuration's native
reference wavelength and let

```text
AOD760 = s760 * tau_ref,
```

where `s760` is the fixed aerosol extinction-ratio interpolation evaluated at
760 nm. At fixed microphysics,

```text
dF/dAOD760 = (dF/dtau_ref) / s760.
```

For the positive OE coordinate `u = log(AOD760)`, this simplifies to

```text
dF/du = AOD760 * dF/dAOD760
      = tau_ref * dF/dtau_ref.
```

This simplification depends on fixed microphysics. If refractive index or size
distribution later enters the state, the derivative of `s760` must also be
included.

### Aerosol height

The altitude-form model profile is `LogNormal(log(z0), sigma0)`, but its native
linearization is with respect to physical `z0` in km. For the positive OE
coordinate `v = log(z0)`:

```text
dF/dv = z0 * dF/dz0.
```

The same diagonal chain-rule matrix transforms the associated prior/posterior
covariances.

### Other units

- CO2 columns are derivatives per unit VMR, not per ppm; multiply by `1e-6`
  for a ppm-coordinate Jacobian.
- Surface pressure follows the model input unit, hPa.
- `SIF760` is spectral radiance per cm-1 at 760 nm.
- `mSIF` is the derivative of that radiance with respect to wavenumber in
  cm-1. Wavelength-domain SIF priors must be converted before retrieval.

## Adding another retrieval flavor

1. Define a zero-field subtype of `AbstractJacobianFlavor`.
2. Implement
   `jacobian_plan(::MyFlavor, params, model, lin_model)::JacobianPlan`.
3. Build stable global `ParameterKey`s and names.
4. For every band, place atmospheric keys first in native component order,
   followed by contiguous surface/source blocks.
5. Supply the corresponding sorted native atmospheric column indices.
6. Override either upstream-work trait when appropriate:

   ```julia
   requires_aerosol_microphysics_jacobians(::MyFlavor) = false
   requires_h2o_jacobians(::MyFlavor) = false
   ```

7. Test the local/global maps, compact allocation dimensions, exact equality
   against selected full-model columns, and central finite differences for
   every physical parameter class.

No retrieval-specific method should be added to elemental, doubling, or
interaction kernels.

## Validation record

The regression is `test/test_selective_jacobians.jl`. Run it from `test/`:

```bash
julia --project=. -e '
    using Test, vSmartMOM
    include("test_selective_jacobians.jl")
'
```

The test covers:

1. exact OCO plan dimensions, names, local/global scattering, and zero fill;
2. exact forward and Jacobian equality between compact propagation and the
   corresponding columns of a full native run;
3. proof that the cached aerosol and gas derivative dimensions are compact;
4. independent two-sided finite differences.

Central finite-difference results recorded on 2026-08-27:

| Parameter class | Maximum absolute error | Maximum relative error |
|:---|---:|---:|
| aerosol `tau_ref` | `1.84e-10` | `2.77e-9` |
| aerosol physical `z0` | `3.39e-12` | `1.06e-6` |
| selected layer absorption | `1.08e-11` | `2.50e-8` |
| surface pressure | `8.84e-12` | `8.79e-6` |
| surface `P0` | `1.17e-12` | `1.78e-7` |
| surface `P1` | `1.17e-12` | `1.78e-7` |
| surface `P2` | `1.17e-12` | `1.78e-7` |
| `SIF760` | `1.69e-11` | `8.74e-7` |
| `mSIF` | `1.32e-9` | `4.58e-7` |

Finite differences run on CPU/Float64. The implementation remains
architecture-generic through the existing array-type abstraction; GPU
execution requires the ordinary CUDA regression on a host with a functional
driver.

## Files introduced or changed

- `src/CoreRT/parameter_layout.jl`: abstract/full/active layouts, stable keys,
  plan/flavor types, validation and accessors.
- `src/CoreRT/tools/jacobian_plans.jl`: extension function, upstream-work
  traits, `OCO_RRS_synth` compiler, model-construction dispatch, globalizer.
- `src/CoreRT/types_lin.jl`: planned-model wrapper and abstract layout field on
  linearized observer results.
- `src/CoreRT/tools/lin_model_from_parameters.jl`: optional suppression of
  H2O and aerosol-microphysics tangent generation while preserving forward
  physics.
- `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties_lin.jl`:
  component-local native selection, compact invariant cache, and forward-only
  aerosol phase path when no Mie tangent is active.
- `src/CoreRT/rt_run_lin.jl`: plan-aware public entry point, source/layout
  validation, compact allocations, and compact-cache plumbing.
- `src/CoreRT/CoreRT.jl` and `src/vSmartMOM.jl`: include/export surface.
- `test/test_selective_jacobians.jl` and `test/runtests.jl`: mapping,
  allocation, equivalence, and finite-difference regressions.
- `docs/src/pages/jacobians.md`,
  `docs/src/pages/concepts/06_linearization.md`, `AGENTS.md`, and the RRS-XCO2
  inversion guides: user/developer documentation and retrieval-boundary rules.
