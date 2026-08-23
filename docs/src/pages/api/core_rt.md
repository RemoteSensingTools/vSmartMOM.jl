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
vSmartMOM.CoreRT.AbstractFourierConvergence
vSmartMOM.CoreRT.AllFourierMoments
vSmartMOM.CoreRT.IntensityConvergence
```

## Single-Scattering Correction

```@docs
vSmartMOM.CoreRT.AbstractSingleScatteringCorrection
vSmartMOM.CoreRT.NoSSCorrection
vSmartMOM.CoreRT.TMSCorrection
vSmartMOM.CoreRT.rt_run_ss_exact
```

## Fourier Loop Bound

The azimuthal Fourier loop bound is derived per band at model build: each
active component (surface, sources, Rayleigh, aerosol modes) declares the
highest Fourier order it can contribute via the `component_m_max` trait, and
the band bound is their maximum, clamped by the user's `max_m`. At run time,
the [`IntensityConvergence`](@ref vSmartMOM.CoreRT.IntensityConvergence)
strategy above may terminate the loop below that ceiling once the accumulated
Stokes-I contribution at the user view angles stabilises — the two mechanisms
compose.

```@docs
vSmartMOM.CoreRT.component_m_max
vSmartMOM.CoreRT._aggregate_m_max
```

## Column Core and TOA Fast Path

Every public forward entry point funnels into the internal single-column
solver `_rt_run_column`; `rt_run_toa` is the lean TOA-only wrapper used by
retrieval loops that never consume BOA or interior fields.

```@docs
vSmartMOM.CoreRT._rt_run_column
vSmartMOM.CoreRT.rt_run_toa
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

## Atmosphere/Surface Split

Run the atmosphere phase once, replay the surface phase per BRDF — see the
[Fast Re-runs & Batch Processing](../batch_processing.md) guide.

```@docs
vSmartMOM.rt_run_atmosphere
vSmartMOM.rt_run_surface
vSmartMOM.rt_run_multi_surface
vSmartMOM.AtmosphereRTCache
```

## Analytic Lambertian Albedo Closure

```@docs
vSmartMOM.LambertianClosure
vSmartMOM.lambertian_closure
vSmartMOM.albedo_jacobian
vSmartMOM.invert_albedo
```

## Scenario Sweeps

```@docs
vSmartMOM.ScenarioSweep
vSmartMOM.SweepResult
vSmartMOM.run_sweep
vSmartMOM.remake_geometry
```

## Stream-Level RT Output

```@docs
vSmartMOM.CoreRT.StreamRTResult
vSmartMOM.CoreRT.rt_run_streams(::Any)
```

## Jacobian Parameter Layout

```@docs
vSmartMOM.CoreRT.ParameterLayout
vSmartMOM.CoreRT.n_total
vSmartMOM.CoreRT.aerosol_range
vSmartMOM.CoreRT.gas_range
vSmartMOM.CoreRT.gas_profile_range
vSmartMOM.CoreRT.gas_layer_index
vSmartMOM.CoreRT.surface_range
vSmartMOM.CoreRT.surface_index
vSmartMOM.CoreRT.surface_parameter_count
vSmartMOM.CoreRT.sif_range
vSmartMOM.CoreRT.atmosphere_range
vSmartMOM.CoreRT.n_layer_params
vSmartMOM.CoreRT.canopy_range
vSmartMOM.sif_reference_state
```

## Quadrature Types

```@docs
vSmartMOM.CoreRT.AbstractQuadratureType
vSmartMOM.CoreRT.RadauQuad
vSmartMOM.CoreRT.GaussLegQuad
```

## Quadrature Construction

```@docs
vSmartMOM.CoreRT.rt_set_streams
```
