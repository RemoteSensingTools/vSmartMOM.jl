# IO Overview

The IO submodule centralizes loading and validation of simulation inputs, decoupled from the CoreRT solvers.

- Multiple dispatch entry points:
  - `read_parameters(path::AbstractString)` → load from YAML (via formats registry)
  - `read_parameters(dict::Dict)` → convert an in-memory config
  - `read_atmos_profile(path|dict)` → load an atmospheric profile only
- Format registry: `IO.Formats` selects a loader based on source/extension (YAML supported).
- Safe parsing: enums and types are parsed with explicit maps. Surfaces (BRDF) are parsed safely.
- **Exception**: `spec_bands` uses `eval()` to support flexible Julia range expressions and Unitful conversions (e.g., `"(1e7/777):0.015:(1e7/757)"`). All other fields avoid eval.

## Quick start

```julia
using vSmartMOM

# 1) YAML file
params = read_parameters(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))
model  = model_from_parameters(params)
R, T   = rt_run(model)

# 2) Dict in memory
cfg = Dict(
  "radiative_transfer" => Dict(
    "spec_bands" => ["(1e7/777):0.015:(1e7/757)"],
    "surface" => ["LambertianSurfaceScalar(0.15)"],
    "quadrature_type" => "GaussQuadFullSphere()",
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

- YAML: registered by default. Add more by calling `IO.Formats.register_format`.

## Safety maps

- Quadrature: `RadauQuad`, `GaussQuadHemisphere`, `GaussQuadFullSphere`
- Polarization: `Stokes_I`, `Stokes_IQ`, `Stokes_IQU`, `Stokes_IQUV`
- Architecture: `CPU`, `GPU`, `default_architecture`
- Decomposition: `NAI2`, `PCW`
- Broadening: `Voigt`, `Lorentz`, `Doppler`
- Complex error function: `HumlicekWeidemann32SDErrorFunction`
- Surfaces (BRDF): LambertianSurfaceScalar, LambertianSurfaceSpectrum, LambertianSurfaceLegendre, LambertianSurfaceSpline, rpvSurfaceScalar, RossLiSurfaceScalar, CoxMunkSurface, CanopySurface

All are passed in as strings in the YAML and mapped safely (without eval, except for spec_bands range expressions).
