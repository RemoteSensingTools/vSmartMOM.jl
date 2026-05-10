# IO Overview

The IO submodule centralizes loading and validation of simulation inputs, decoupled from the CoreRT solvers.

Use `read_parameters` as the user-facing entry point. It dispatches on the input type:

- `read_parameters(path::AbstractString)` loads a registered file format.
- `read_parameters(dict::Dict)` converts an in-memory config.
- `read_parameters(src::IOSource)` loads a typed source such as GEOS-Chem NetCDF.
- `read_parameters(params::vSmartMOM_Parameters)` passes an already parsed parameter object through.

The explicit aliases `parameters_from_file`, `parameters_from_dict`, and `parameters_from_source` remain available when code wants the call site to state the input kind. `parameters_from_yaml` is also supported for YAML files.

- `read_atmos_profile(path|dict)` loads an atmospheric profile only.
- Format registry: `IO.Formats` selects a loader based on source/extension (YAML and TOML supported).
- Safe parsing: enums and types are parsed with explicit maps. Surfaces (BRDF) and spectral band ranges are parsed without `eval`.

## Quick start

```julia
using vSmartMOM

# 1) YAML or TOML file
params = read_parameters(joinpath(pkgdir(vSmartMOM), "src", "CoreRT", "DefaultParameters.yaml"))
model  = model_from_parameters(params)
R, T   = rt_run(model)

# 2) Dict in memory
cfg = Dict(
  "radiative_transfer" => Dict(
    "spec_bands" => ["(1e7/777):0.015:(1e7/757)"],
    "surface" => ["LambertianSurfaceScalar(0.15)"],
    "quadrature_type" => "GaussLegQuad()",
    "polarization_type" => "Stokes_I()",
    "max_m" => 3,
    "Δ_angle" => 2.0,
    "l_trunc" => 20,
    "depol" => 0.0,
    "float_type" => "Float64",
    "architecture" => "default_architecture",
  ),
  "geometry" => Dict("sza"=>60, "vza"=>[60,45,30,15,0,15,30,45,60], "vaz"=>[180,180,180,180,0,0,0,0,0], "obs_alt"=>1000.0),
  "atmospheric_profile" => Dict("T"=>fill(260.0,34), "p"=>collect(range(0.14, stop=1005.0, length=35)), "profile_reduction"=>-1),
)
params2 = read_parameters(cfg)
```

## Formats

- YAML and TOML: registered by default. Add more by calling `IO.Formats.register_format`.
- GEOS-Chem and generic NetCDF sources use typed `IOSource` objects and the same `read_parameters` entry point.

## Safety maps

- Quadrature: `RadauQuad`, `GaussLegQuad`
- Polarization: `Stokes_I`, `Stokes_IQ`, `Stokes_IQU`, `Stokes_IQUV`
- Architecture: `CPU`, `GPU`, `default_architecture`
- Decomposition: `NAI2`, `PCW`
- Broadening: `Voigt`, `Lorentz`, `Doppler`
- Complex error function: `HumlicekWeidemann32SDErrorFunction`
- Surfaces (BRDF): LambertianSurfaceScalar, LambertianSurfaceSpectrum, LambertianSurfaceLegendre, LambertianSurfaceSpline, rpvSurfaceScalar, RossLiSurfaceScalar, CoxMunkSurface, CanopySurface

All are passed in as strings in the YAML and mapped safely without `eval`.
