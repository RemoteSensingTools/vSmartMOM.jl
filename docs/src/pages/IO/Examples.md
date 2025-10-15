# IO examples

## Load from YAML file

```julia
using vSmartMOM
params = read_parameters(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))
model  = model_from_parameters(params)
R, T   = rt_run(model)
```

## Build from a Dict

```julia
cfg = Dict(
  "radiative_transfer" => Dict(
    "spec_bands" => ["(1e7/777):0.02:(1e7/757)"],
    "surface" => ["LambertianSurfaceScalar(0.2)"],
    "quadrature_type" => "RadauQuad()",
    "polarization_type" => "Stokes_IQU()",
    "max_m" => 3,
    "Δ_angle" => 2.0,
    "l_trunc" => 20,
    "depol" => 0.0,
    "float_type" => "Float32",
    "architecture" => "CPU()",
  ),
  "geometry" => Dict("sza"=>30, "vza"=>[0, 30, 60], "vaz"=>[0, 0, 0], "obs_alt"=>1000.0),
  "atmospheric_profile" => Dict("T"=>fill(255.0,3), "p"=>[0.1, 500.0, 900.0, 1005.0], "profile_reduction"=>-1),
)
params = read_parameters(cfg)
```

## Surfaces (BRDF)

You can specify different surface types per band, e.g.,

```yaml
radiative_transfer:
  surface:
    - LambertianSurfaceScalar(0.1)
    - rpvSurfaceScalar(0.08, 0.3, -0.2, 30.0)
```

## Atmospheric profile only

```julia
using vSmartMOM
ap = read_atmos_profile(joinpath(dirname(pathof(vSmartMOM)), "CoreRT", "DefaultParameters.yaml"))
```
