# Scattering Module Example

This example shows a complete NAI2 workflow, optional PCW setup, and AD usage.

```julia
using vSmartMOM.Scattering
using Distributions
using FastGaussQuadrature
using Parameters

#
# STEP 1: Define aerosol size distribution and refractive index
#
r_m = 0.30               # median radius [um]
sigma_r = 2.0            # geometric stddev [-]
n_r = 1.30               # real refractive index
n_i = 0.001              # imaginary refractive index (use positive convention)
r_max = 30.0             # [um]
nquad_radius = 1500      # quadrature points for size integration

size_distribution = LogNormal(log(r_m), log(sigma_r))
aero = Aerosol(size_distribution, n_r, n_i)

#
# STEP 2: Build Mie model settings
#
lambda_um = 0.55
polarization = Stokes_IQUV()
truncation = δBGE(20, 2.0)   # l_max = 20, exclusion angle = 2 deg

model_NAI2 = make_mie_model(
    NAI2(),
    aero,
    lambda_um,
    polarization,
    truncation,
    r_max,
    nquad_radius,
)

#
# STEP 3: Compute aerosol optical properties
#
aerosol_optics = compute_aerosol_optical_properties(model_NAI2)

println("omega_tilde = ", aerosol_optics.ω̃)
println("k_ext       = ", aerosol_optics.k)
println("f_t         = ", aerosol_optics.fᵗ)

@unpack α, β, γ, δ, ϵ, ζ = aerosol_optics.greek_coefs
println("Greek coefficient length = ", length(β))

#
# STEP 4: Reconstruct phase-matrix elements on a custom angular grid
#
mu_quad, _ = gausslegendre(500)
scattering_matrix = reconstruct_phase(aerosol_optics.greek_coefs, mu_quad)
println("f11 at forward angle = ", scattering_matrix.f₁₁[end])

#
# STEP 5: Optional AD run (returns ONE AerosolOptics object; Jacobian in `.derivs`)
#
aerosol_optics_ad = compute_aerosol_optical_properties(model_NAI2; autodiff=true)
println("AD derivative array size = ", size(aerosol_optics_ad.derivs))

#
# STEP 6: Convenience APIs for cross-sections / phase function
#
mu_pf, w_mu_pf, p11, C_ext, C_sca, g = phase_function(aero, lambda_um, r_max, nquad_radius)
println("phase_function: C_ext = ", C_ext, ", C_sca = ", C_sca, ", g = ", g)

XS_ext, XS_sca, Cext_eff, Csca_eff = compute_aerosol_XS(aero, lambda_um, r_max, nquad_radius)
println("compute_aerosol_XS: XS_ext = ", XS_ext, ", XS_sca = ", XS_sca)
println("compute_aerosol_XS: Cext_eff = ", Cext_eff, ", Csca_eff = ", Csca_eff)

#
# OPTIONAL: PCW model setup
# - Provide precomputed Wigner matrices from file:
#     model_PCW = make_mie_model(PCW(), aero, lambda_um, polarization, truncation,
#                                r_max, nquad_radius, "path/to/wigner_values.jld")
#
# - Or compute Wigner matrices directly (expensive for large N_max):
#     N_max = 120
#     wigner_A, wigner_B = compute_wigner_values(N_max)
#     model_PCW = make_mie_model(PCW(), aero, lambda_um, polarization, truncation,
#                                r_max, nquad_radius, wigner_A, wigner_B)
```
