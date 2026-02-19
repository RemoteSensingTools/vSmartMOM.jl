# # Quick Start: vSmartMOM Radiative Transfer
#
# vSmartMOM is a polarized radiative transfer solver using the adding-doubling method.
# This 5-minute tutorial runs a forward RT simulation and computes Jacobians.
#
# ## 1) Load the package

using vSmartMOM

# ## 2) Load parameters from YAML
#
# Parameters define the atmosphere, surface, geometry, and spectral bands.
# See the IO Schema documentation for the full YAML format.

yaml_path = joinpath(dirname(dirname(pathof(vSmartMOM))),
                     "test", "test_parameters", "ParamsEMIT_fast.yaml")
params = parameters_from_yaml(yaml_path)
params.architecture = vSmartMOM.Architectures.CPU()

# ## 3) Build the model and run forward RT

model = model_from_parameters(params)
R, T = rt_run(model)

# ## 4) Interpret results
#
# - **R** = reflectance at top-of-atmosphere (TOA) `[nVZA × nStokes × nSpec]`
# - **T** = transmittance at bottom-of-atmosphere (BOA)
#
# Dimensions: view zenith angles × Stokes components (I,Q,U,V) × spectral points.

println("R shape: ", size(R))
println("T shape: ", size(T))
println("R(nadir, I, λ₁) = ", R[1, 1, 1])

# ## 5) Linearized RT for Jacobians
#
# Use `model_from_parameters(LinMode(), params)` and `rt_run` to get analytic
# derivatives of R and T with respect to atmospheric and surface parameters.
# No finite differences required — Jacobians come from a single RT pass.

using vSmartMOM.CoreRT

model_lin, lin_model = model_from_parameters(LinMode(), params)
NAer = length(params.scattering_params.rt_aerosols)
NGas = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R_lin, T_lin, dR, dT = rt_run(model_lin, lin_model, NAer, NGas, NSurf)

# `dR` has shape `[nParams × nVZA × nStokes × nSpec]` — derivatives w.r.t. aerosols,
# gas VMRs, and surface albedo.

println("dR shape (Jacobian): ", size(dR))

# ## 6) GPU support
#
# For NVIDIA GPUs with CUDA.jl, set `params.architecture = vSmartMOM.Architectures.GPU()`
# before building the model. The same `model_from_parameters` and `rt_run` calls apply.
# See the GPU tutorial for more.
