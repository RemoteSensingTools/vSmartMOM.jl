## 
## Using RadiativeTransfer.Scattering to perform Mie computations
## 

using Revise
using RadiativeTransfer
using RadiativeTransfer.Scattering

# ? Should we make a wrapper around LogNormal?
using Distributions
# using BenchmarkTools

## 
## STEP 1: Create the Aerosol
## 

# Aerosol particle distribution and properties 
μ  = 0.3                # Log mean radius
σ  = 6.82               # Log stddev of radius
r_max = 30.0            # Maximum radius
nquad_radius = 2500     # Number of quadrature points for integrating of size dist.
nᵣ = 1.3                # Real part of refractive index
nᵢ = 0.001                # Imag part of refractive index

size_distribution = LogNormal(log(μ), log(σ))

# Create the aerosol
aero = make_univariate_aerosol(size_distribution, r_max, nquad_radius, nᵣ, nᵢ)

## 
## STEP 2: Create the Mie Calculations model
## 

λ = 0.55   # Incident wavelength
polarization_type = Stokes_IQUV()
truncation_type = δBGE(10, 10)
wigner_file_path = "/home/rjeyaram/RadiativeTransfer/src/Scattering/Mie/wigner_values.jld"

model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type)
model_PCW = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, wigner_file_path)

## 
## STEP 3: Perform the Mie Calculations
## 

# @btime aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
# @btime aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

# Compute the derivatives!
aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

# aerosol_optics_NAI2.greek_coefs.α ≈ aerosol_optics_PCW.greek_coefs.α

# aerosol_optics_NAI2.greek_coefs ≈ aerosol_optics_PCW.greek_coefs
aerosol_optics_NAI2 ≈ aerosol_optics_PCW
