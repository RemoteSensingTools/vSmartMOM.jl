# # IO Quick Tutorial
#
# Goal: Load parameters from YAML and Dict, construct a model, and run RT.
#
# This tutorial demonstrates the IO submodule entry points and the minimal
# set of required fields.

# ## 1) From YAML

using vSmartMOM

params = vSmartMOM.read_parameters(joinpath(pkgdir(vSmartMOM), "src", "CoreRT", "DefaultParameters.yaml"))
model  = vSmartMOM.model_from_parameters(params)
R, T   = vSmartMOM.rt_run(model)
println("Reflectance shape: ", size(R))

# ## 2) From a Dict

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
      "T"=>fill(255.0, 3),
      "p"=>[0.1, 500.0, 900.0, 1005.0],
      "profile_reduction"=>-1,
  ),
)
params2 = vSmartMOM.read_parameters(cfg)
model2  = vSmartMOM.model_from_parameters(params2)
R2, T2  = vSmartMOM.rt_run(model2)
println("Dict-based reflectance shape: ", size(R2))

# ## Mandatory vs optional
# Mandatory top-level sections: `radiative_transfer`, `geometry`, `atmospheric_profile`.
# Optional: `absorption`, `scattering`.
# See the IO Schema doc page for all keys and constraints.
