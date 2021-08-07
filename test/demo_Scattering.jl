## 
## Using RadiativeTransfer.Scattering to perform Mie computations
## 

using Revise
using RadiativeTransfer
using RadiativeTransfer.Scattering

# ? Should we make a wrapper around LogNormal?
using Distributions
using FastGaussQuadrature
# using BenchmarkTools

## 
## STEP 1: Create the Aerosol
## 

# Aerosol particle distribution and properties 
μ  = (π/3) * 0.55 /(2π)                # Log mean radius
σ  = 1.1               # Log stddev of radius
nᵣ = 1.5                # Real part of refractive index
nᵢ = 0.0                # Imag part of refractive index
r_max = 30.0            # Maximum radius
nquad_radius = 2500     # Number of quadrature points for integrating of size dist.

size_distribution = make_log_normal_size_dist(μ, σ)

# Create the aerosol
aero = Aerosol(size_distribution, nᵣ, nᵢ)

## 
## STEP 2: Create the Mie Calculations model
## 

λ = 0.55   # Incident wavelength
polarization_type = Stokes_I()
truncation_type = δBGE(10, 10)
wigner_file_path = "/home/rjeyaram/RadiativeTransfer/src/Scattering/Mie/wigner_values.jld"

model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius)
model_PCW = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_file_path)

## 
## STEP 3: Perform the Mie Calculations
## 

# @btime aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
# @btime aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

# Compute the derivatives!
aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2, autodiff=true);
aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

μ_quad, w_μ = gausslegendre(1000)
f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = Scattering.reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ_quad);

# aerosol_optics_NAI2.greek_coefs.α ≈ aerosol_optics_PCW.greek_coefs.α

# aerosol_optics_NAI2.greek_coefs ≈ aerosol_optics_PCW.greek_coefs
aerosol_optics_NAI2 ≈ aerosol_optics_PCW
