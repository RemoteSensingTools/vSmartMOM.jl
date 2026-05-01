# # IO Quick Tutorial
#
# Goal: Load parameters from YAML and Dict, construct a model, and run RT.
#
# This tutorial demonstrates the IO submodule entry points and the minimal
# set of required fields.

# ## 1) From YAML

using vSmartMOM

params = vSmartMOM.read_parameters(joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml"))
model  = vSmartMOM.model_from_parameters(params)
R, T   = vSmartMOM.rt_run(model)
println("Reflectance shape: ", size(R))

# ## 2) From a Dict

cfg = Dict(
  "radiative_transfer" => Dict(
    "spec_bands" => ["[12987.0]"],
    "surface" => ["LambertianSurfaceScalar(0.2)"],
    "quadrature_type" => "GaussQuadHemisphere()",
    "polarization_type" => "Stokes_I()",
    "max_m" => 1,
    "Δ_angle" => 2.0,
    "l_trunc" => 1,
    "depol" => 0.0,
    "float_type" => "Float64",
    "architecture" => "CPU()",
  ),
  "geometry" => Dict("sza"=>30, "vza"=>[0], "vaz"=>[0], "obs_alt"=>1000.0),
  "atmospheric_profile" => Dict(
      "T"=>[260.0, 280.0],
      "p"=>[100.0, 600.0, 1000.0],
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
