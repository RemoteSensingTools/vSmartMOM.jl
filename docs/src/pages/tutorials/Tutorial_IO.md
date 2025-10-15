```@meta
EditURL = "Tutorial_IO.jl"
```

# IO quick tutorial

Goal: Load parameters from YAML and Dict, construct a model, and run RT.

This tutorial demonstrates the IO submodule entry points and the minimal
set of required fields.

## 1) From YAML
```julia
using vSmartMOM
const VS = vSmartMOM
params = VS.read_parameters(joinpath(dirname(pathof(VS)), "CoreRT", "DefaultParameters.yaml"))
model  = VS.model_from_parameters(params)
R, T   = VS.rt_run(model)
size(R)
```

## 2) From a Dict
```julia
using Unitful
cfg = Dict(
  "radiative_transfer" => Dict(
    "spec_bands" => ["(1e7/777):0.02:(1e7/757)"],
    "surface" => ["LambertianSurfaceScalar(0.2)"],
    "quadrature_type" => "RadauQuad()",
    "polarization_type" => "Stokes_I()",
    "max_m" => 3,
    "Δ_angle" => 2.0,
    "l_trunc" => 20,
    "depol" => 0.0,
    "float_type" => "Float64",
    "architecture" => "default_architecture",
  ),
  "geometry" => Dict("sza"=>30, "vza"=>[0, 30, 60], "vaz"=>[0, 0, 0], "obs_alt"=>1000.0),
  "atmospheric_profile" => Dict(
      "T"=>fill(255.0,3),
      "p"=>[0.1, 500.0, 900.0, 1005.0],
      "profile_reduction"=>-1,
  ),
)
params2 = VS.read_parameters(cfg)
model2  = VS.model_from_parameters(params2)
R2, T2  = VS.rt_run(model2)
size(R2)
```

## Mandatory vs optional
Mandatory top-level sections: `radiative_transfer`, `geometry`, `atmospheric_profile`.
Optional: `absorption`, `scattering`.
See the IO → Schema doc page for all keys and constraints.

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

