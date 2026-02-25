# # Jacobians: Linearized Radiative Transfer
#
# This tutorial shows how to compute analytic Jacobians (derivatives of
# top-of-atmosphere radiance with respect to atmospheric and surface parameters)
# using the linearized RT mode in vSmartMOM.
#
# The linearized mode propagates derivatives through the adding-doubling
# solver alongside the forward radiance, producing exact Jacobians in a
# single RT pass — no finite differences needed.
#
# ## 1) Setup: load a test case

using vSmartMOM
using vSmartMOM.CoreRT

params = parameters_from_yaml(
    joinpath(dirname(dirname(pathof(vSmartMOM))),
             "test", "test_parameters", "ParamsEMIT_fast.yaml"))
params.architecture = vSmartMOM.Architectures.CPU()

# ## 2) Build the linearized model
#
# `LinMode()` tells the constructor to compute both the forward model and
# the derivative containers (`vSmartMOM_Lin`).

model, lin_model = model_from_parameters(LinMode(), params)

# `model` is the same forward model as `model_from_parameters(params)`.
# `lin_model` holds pre-computed ∂τ/∂x arrays for gases, aerosols, and
# aerosol optics.

println("τ̇_abs bands: ", length(lin_model.τ̇_abs),
        "  shape[1]: ", size(lin_model.τ̇_abs[1]))
println("τ̇_aer bands: ", length(lin_model.τ̇_aer),
        "  shape[1]: ", size(lin_model.τ̇_aer[1]))

# ## 3) Run linearized RT
#
# The linearized `rt_run` returns `(R, T, dR, dT)`:
# - `R`  — reflected Stokes field  `[nVZA × nStokes × nSpec]`
# - `T`  — transmitted field
# - `dR` — Jacobian of R w.r.t. all parameters  `[nVZA × nStokes × nSpec × nParams]`
# - `dT` — Jacobian of T

NAer  = length(params.scattering_params.rt_aerosols)
NGas  = size(lin_model.τ̇_abs[1], 1)
NSurf = 1

R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)

println("R  shape: ", size(R))
println("dR shape: ", size(dR))

# ## 4) Interpret the Jacobian layout
#
# The derivative dimension of `dR` is ordered as:
#
# | Index range                | Parameter                        |
# |:---------------------------|:---------------------------------|
# | `1 : NAer*7`               | Aerosol properties (7 per type)  |
# | `NAer*7+1 : NAer*7+NGas`   | Gas VMRs (per variable molecule) |
# | `NAer*7+NGas+1 : end`      | Surface albedo                   |
#
# For this test case:

Nparams = NAer * 7 + NGas + NSurf
println("Total Jacobian parameters: ", Nparams,
        " (NAer×7=", NAer*7, ", NGas=", NGas, ", NSurf=", NSurf, ")")

# ## 5) Inspect individual derivatives
#
# Surface albedo Jacobian (last parameter, first VZA, Stokes-I):

dR_albedo = dR[:, 1, :, end]
println("dR/d(albedo) at nadir, first 5 spectral points: ",
        round.(dR_albedo[1, 1:min(5,end)], digits=6))

# Gas VMR Jacobians:

if NGas > 0
    igas_start = NAer * 7 + 1
    dR_gas1 = dR[:, 1, :, igas_start]
    println("dR/d(gas₁ VMR) at nadir, first 5 points: ",
            round.(dR_gas1[1, 1:min(5,end)], digits=6))
end

# ## 6) Verify against finite differences (optional)
#
# For validation, one can perturb a parameter by ε and compare:
#
# ```julia
# ε = 1e-4
# params_pert = deepcopy(params)
# # e.g. perturb surface albedo
# # ... modify params_pert ...
# model_pert = model_from_parameters(params_pert)
# R_pert, = rt_run(model_pert)
# K_fd = (R_pert .- R) ./ ε
# # Compare K_fd with dR[:, :, :, end]
# ```
#
# The test suite in `test/test_forward_lin.jl` performs this comparison
# systematically for all parameter types.
