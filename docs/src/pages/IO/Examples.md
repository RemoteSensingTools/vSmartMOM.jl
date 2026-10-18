# IO examples

## Load from YAML file

```julia
using vSmartMOM
params = read_parameters(joinpath(pkgdir(vSmartMOM), "src", "CoreRT", "DefaultParameters.yaml"))
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
    "l_trunc" => 20,
    "depol" => 0.0,
    "float_type" => "Float32",
    "architecture" => "CPU()",
  ),
  "geometry" => Dict("sza"=>30, "vza"=>[0, 30, 60], "vaz"=>[0, 0, 0], "obs_alt"=>[0.0]),
  "atmospheric_profile" => Dict("T"=>fill(255.0,3), "p"=>[0.1, 500.0, 900.0, 1005.0], "profile_reduction"=>-1),
)
params = read_parameters(cfg)
```

Here `[0.0]` requests TOA upwelling and BOA downwelling. To add radiances at
5 km above BOA, use `"obs_alt" => [0.0, 5.0]` and read the interior pair from
`rt_run(model).levels`.

## Surfaces (BRDF)

You can specify different surface types per band, e.g.,

```yaml
radiative_transfer:
  surface:
    - LambertianSurfaceScalar(0.1)
    - rpvSurfaceScalar(0.08, 0.3, -0.2, 30.0)
```

## Atmospheric profile only

The standalone profile readers take a mapping with `T`, `p_half` (interface
pressures in hPa, TOA first) and optional `q` and `vmr` at its top level; a
scene file such as `DefaultParameters.yaml` (profile nested under
`atmospheric_profile`, interfaces named `p`) is not accepted.

```julia
using vSmartMOM
ap = read_atmos_profile_dict(Dict(
  "T"      => [220.0, 250.0, 280.0],
  "p_half" => [0.1, 300.0, 700.0, 1000.0],
  "q"      => [1e-5, 1e-3, 5e-3],
  "vmr"    => Dict("CO2" => 4.0e-4),
))
# or, from a YAML file holding the same keys:
# ap = read_atmos_profile("profile.yaml")
```
