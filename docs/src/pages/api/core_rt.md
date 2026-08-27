# CoreRT API

`CoreRT` owns the adding-doubling solver, the RT model containers, and the
linearized/Jacobian-facing model layout.

## Forward and Linearized Modes

```@docs
vSmartMOM.FwdMode
vSmartMOM.LinMode
vSmartMOM.model_from_parameters_lin
vSmartMOM.rt_run_lin
```

## Model Types

```@docs
vSmartMOM.AbstractRTModel
vSmartMOM.RTModel
vSmartMOM.SolverConfig
vSmartMOM.Atmosphere
vSmartMOM.RayleighScattering
vSmartMOM.AerosolState
vSmartMOM.Optics
vSmartMOM.OpticsLin
vSmartMOM.CoreRT.vSmartMOM_Parameters
vSmartMOM.CoreRT.AtmosphericProfile
vSmartMOM.CoreRT.ObsGeometry
vSmartMOM.CoreRT.QuadPoints
vSmartMOM.CoreRT.CompositeLayer
vSmartMOM.CoreRT.AddedLayer
vSmartMOM.CoreRT.RTNumericalParameters
```

## Height-Aware Forward and Linearized Output

`rt_run` returns a named `ObserverRTResult`. Endpoint fields are `nothing`
when not selected by `geometry.obs_alt`; strict-interior outputs are stored as
`LevelRadiance` records in `result.levels`. See the
[`geometry` schema](../IO/Schema/geometry.md) for the scalar/vector selection
convention.

The analytic linearized solver returns `ObserverRTResultLin`. It stores
endpoint Jacobians alongside the endpoint radiances and uses
`LevelRadianceLin` records for co-located forward fields and Jacobians at each
strict-interior height. Both result types retain their historical tuple
iteration conventions.

```@docs
vSmartMOM.ObserverRTResult
vSmartMOM.LevelRadiance
vSmartMOM.ObserverRTResultLin
vSmartMOM.LevelRadianceLin
vSmartMOM.total_downwelling
vSmartMOM.total_downwelling_jacobian
```

> Source-term types (`SolarBeam`, `BlackbodySource`, `SurfaceSIF`,
> `SourceSet`, `NoSource`, AD-mode traits, `prepare_source` /
> `surface_source_contribute!` dispatchers) are documented on the
> [Source terms](../extending/sources.md#api-reference) page.

## Batch Scene Processing

`BatchContext` enables efficient multi-scene loop workflows (e.g. ensemble retrievals,
parameter sweeps) by caching expensive one-time setup work (Mie, Fourier decomposition,
HITRAN parsing) and exposing cheap per-scene update functions.

```@docs
vSmartMOM.BatchContext
vSmartMOM.update_model!
vSmartMOM.update_aerosol_loading!
vSmartMOM.update_aerosol_microphysics!
```

## Stream-Level RT Output

```@docs
vSmartMOM.CoreRT.StreamRTResult
vSmartMOM.CoreRT.rt_run_streams(::Any)
```

## Jacobian Parameter Layout

```@docs
vSmartMOM.CoreRT.AbstractParameterLayout
vSmartMOM.CoreRT.ParameterLayout
vSmartMOM.CoreRT.ParameterKey
vSmartMOM.CoreRT.ActiveParameterLayout
vSmartMOM.CoreRT.n_total
vSmartMOM.CoreRT.aerosol_range
vSmartMOM.CoreRT.gas_range
vSmartMOM.CoreRT.gas_profile_range
vSmartMOM.CoreRT.gas_layer_index
vSmartMOM.CoreRT.atmosphere_range
vSmartMOM.CoreRT.surface_range
vSmartMOM.CoreRT.surface_index
vSmartMOM.CoreRT.sif_range
vSmartMOM.CoreRT.n_layer_params
vSmartMOM.CoreRT.canopy_range
vSmartMOM.CoreRT.surface_parameter_count
```

## Retrieval-Selected Jacobian Plans

The policy types below compile a retrieval definition into a small band-local
tangent space before the layer optical properties and MOM operators are
allocated. See [Compute Jacobians](../jacobians.md#Retrieval-Selected-Jacobians)
for the `OCO_RRS_synth` state ordering, units, and chain rules.

```@docs
vSmartMOM.CoreRT.AbstractJacobianFlavor
vSmartMOM.CoreRT.OCO_RRS_synth
vSmartMOM.CoreRT.JacobianPlan
vSmartMOM.CoreRT.PlannedRTModelLin
vSmartMOM.CoreRT.jacobian_plan
vSmartMOM.CoreRT.requires_aerosol_microphysics_jacobians
vSmartMOM.CoreRT.requires_h2o_jacobians
vSmartMOM.CoreRT.globalize_jacobian
vSmartMOM.sif_reference_state
```
