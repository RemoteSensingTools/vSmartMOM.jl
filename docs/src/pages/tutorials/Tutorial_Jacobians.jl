# # Jacobians: Linearized RT Tutorial
#
# ### Introduction
# This tutorial demonstrates how to compute analytic Jacobians
# (derivatives of radiance with respect to physical state parameters)
# using vSmartMOM's linearized RT solver.
#
# You will learn:
# 1. How to set up a linearized model
# 2. How to extract Jacobians from the output
# 3. How to interpret the parameter ordering
# 4. How to compare analytic Jacobians with finite differences
#
# ---
#
# ### Load packages

using vSmartMOM
using vSmartMOM.CoreRT

# ---
#
# ### 1. Build a linearized model
#
# The linearized RT requires a `LinMode()` marker when building the model.
# This tells `model_from_parameters` to also precompute the derivative
# quantities (∂τ/∂x, ∂ω/∂x, ∂Z/∂x) needed by the chain rule.

params = parameters_from_yaml("test/test_parameters/JacobianTestFast.yaml")
model, lin_model = model_from_parameters(LinMode(), params)

# `model` is the usual forward model (optical properties, geometry, etc.).
# `lin_model` is a `vSmartMOM_Lin` struct holding the linearized optical
# properties w.r.t. state-vector elements.

# ---
#
# ### 2. Run the linearized solver
#
# Count the number of aerosol, gas, and surface parameters:

NAer = length(params.scattering_params.rt_aerosols)   # number of aerosol types
NGas = size(lin_model.τ̇_abs[1], 1)                     # number of gas species
NSurf = 1                                              # surface albedo

R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)

# `R` is the reflected radiance: `(nVZA, nStokes, nSpec)`
# `dR` contains the Jacobians: `(nVZA, nStokes, nSpec, Nparams)`
# where `Nparams = NAer * 7 + NGas + NSurf`

println("Reflected radiance shape: ", size(R))
println("Jacobian shape:           ", size(dR))

# ---
#
# ### 3. Parameter ordering
#
# The last dimension of `dR` corresponds to the following parameters
# (for each aerosol type, 7 derivatives):
#
# | Index range           | Parameter                                |
# |-----------------------|------------------------------------------|
# | 1..7                  | Aerosol 1: τ_ref, nᵣ, nᵢ, rₘ, σ_g, p₀, σ_p |
# | 8..7*NAer             | Aerosol 2..NAer (same 7 per aerosol)     |
# | 7*NAer+1..7*NAer+NGas | Gas VMR scaling factors                  |
# | 7*NAer+NGas+1         | Surface albedo                           |
#

Nparams = NAer * 7 + NGas + NSurf
@assert size(dR, 4) == Nparams

# Extract Jacobian w.r.t. surface albedo (last parameter):
dR_dalb = dR[:, :, :, end]
println("dR/d(albedo) max value: ", maximum(abs.(dR_dalb)))

# Extract Jacobian w.r.t. aerosol reference optical depth (first parameter):
dR_dtau = dR[:, :, :, 1]
println("dR/d(τ_ref)  max value: ", maximum(abs.(dR_dtau)))

# ---
#
# ### 4. Quick finite-difference comparison
#
# To verify analytic Jacobians, perturb a parameter and compare.
# Here we perturb the surface albedo.

ε = 1e-5
alb_base = params.brdf[1].albedo

# Perturbed run (forward only — no need for linearized solver)
params_pert = deepcopy(params)
params_pert.brdf[1] = vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}(alb_base + ε)

model_pert = model_from_parameters(params_pert)
R_pert, _, _, _, _, _, _ = rt_run(model_pert)

# Finite-difference Jacobian
dR_fd = (R_pert .- R) ./ ε

# Compare at first VZA, Stokes I, all wavelengths
println("\n--- Surface albedo Jacobian comparison (VZA 1, Stokes I) ---")
for iλ in axes(R, 3)
    analytic = dR_dalb[1, 1, iλ]
    fd       = dR_fd[1, 1, iλ]
    relerr   = abs(analytic - fd) / max(abs(fd), 1e-15)
    println("  λ[$iλ]: analytic=$(round(analytic, sigdigits=6)), FD=$(round(fd, sigdigits=6)), rel_err=$(round(relerr, sigdigits=3))")
end
