# Scattering Module Example

```julia
using vSmartMOM
using vSmartMOM.Scattering
using Distributions
using BenchmarkTools

## 
## STEP 1: Create the Aerosol
## 

# Aerosol particle distribution and properties 
μ  = 0.3                # Log mean radius (μm)
σ  = 6.82               # Log stddev of radius (μm)
r_max = 30.0            # Maximum radius (μm)
nquad_radius = 2500     # Number of quadrature points for integrating of size dist.
nᵣ = 1.3                # Real part of refractive index
nᵢ = 0.0                # Imag part of refractive index

size_distribution = LogNormal(log(μ), log(σ))

# Create the aerosol
aero = Aerosol(size_distribution, nᵣ, nᵢ)

## 
## STEP 2: Create the Mie Calculations model
## 

λ = 0.55        # Incident wavelength (μm)
polarization_type = Stokes_IQUV()  
l_max = 10      # Trunction length for legendre terms
Δ_angle = 2     # Exclusion angle for forward peak (in fitting procedure) `[degrees]`
truncation_type = δBGE(l_max, Δ_angle)

# NAI2 Method
model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius)

# PCW Method with saved/loaded Wigner tables
wigner_file_path = "PATH_TO_SAVED_WIGNER_MATRIX"
model_PCW = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_file_path)

# PCW Method with newly computed Wigner tables
wigner_A, wigner_B = compute_wigner_values(600) # Specify N_max
model_PCW_computed_wigner = make_mie_model(PCW(), aero, λ, polarization_type, truncation_type, wigner_A, wigner_B)

## 
## STEP 3: Perform the Mie Calculations
## 

aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
aerosol_optics_PCW = compute_aerosol_optical_properties(model_PCW);

# To perform auto-differentiation w.r.t. μ, σ, nᵣ, nᵢ
aerosol_optics_NAI2_AD, derivs = compute_aerosol_optical_properties(model_NAI2, autodiff=true);

## 
## STEP 4: Obtain Scattering Matrix
## 

μ_quad, w_μ = gausslegendre(1000); # Quadrature points/weights
scattering_matrix = reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ_quad);

```