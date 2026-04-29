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
```

## Jacobian Parameter Layout

```@docs
vSmartMOM.CoreRT.ParameterLayout
vSmartMOM.CoreRT.n_total
vSmartMOM.CoreRT.aerosol_range
vSmartMOM.CoreRT.gas_range
vSmartMOM.CoreRT.surface_range
vSmartMOM.CoreRT.surface_index
vSmartMOM.CoreRT.n_layer_params
vSmartMOM.CoreRT.canopy_range
```
