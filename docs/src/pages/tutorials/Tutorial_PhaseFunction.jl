# # Simple Mie Phase Function Tutorial
# This is just the first test, more to follow. 
# 
# ---
# In the following, we will just walk through how to compute greek coefficients (which will be require for radiative transfer calculations) as well as actual phase matrices from the Mie theory.

# ---
# First let's use the required packages

using RadiativeTransfer.PhaseFunction
using Distributions
using Plots
#----------------------------------------------------------------------------

## Aerosol particle distribution and properties 
μ  = 0.3        ## Log-normal median radius [μm]
σ  = 2.0        ## Log-normal stddev of radius
r_max = 30.0    ## Maximum radius [μm]
n  = 2500       ## Number of quadrature points for integrating of size dist.
nᵣ = 1.3        ## Real part of refractive index
nᵢ = 0.0        ## Imag part of refractive index (sign changed, use only + here)

## Create a Size Distribution (from Julia's Distributions package)
size_distribution = LogNormal(log(μ), log(σ))

## Create the aerosol
aero = make_univariate_aerosol(size_distribution, r_max, n, nᵣ, nᵢ)
#----------------------------------------------------------------------------

λ = 0.55                             ## Incident wavelength [μm]
polarization_type = Stokes_IQUV()    ## Polarization type 
truncation_type   = δBGE(20, 2)      ## Truncation type
#----------------------------------------------------------------------------

## Create a Mie model, using the Siewert method NAI2
model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type)
#----------------------------------------------------------------------------

## Compute aerosol optical properties:
aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
#----------------------------------------------------------------------------

# --- 
# simple example of how to use the Documentation, add `?` in front and get the DocStrings of aerosol_optics_NAI2 :
# 

?aerosol_optics_NAI2
#----------------------------------------------------------------------------

# ---
# ### Let's plot the greek coefficients
# which are basically giving us the legendre decomposition of the phase matrix components:

using Parameters
@unpack α,β,γ,δ,ϵ,ζ = aerosol_optics_NAI2.greek_coefs
p1 = plot(α,  title="α")
p2 = plot(β,  title="β")
p3 = plot(γ,  title="γ")
p4 = plot(δ,  title="δ")
p5 = plot(ϵ,  title="ϵ")
p6 = plot(ζ,  title="ζ")
plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), legend=false)
xlims!(0,100)
#----------------------------------------------------------------------------

using FastGaussQuadrature
μ, w_μ = gausslegendre(1000)
## Reconstruct Phase Functions from greek coefficients (overkill for Siewert, mostly for Wigner method)
f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = PhaseFunction.reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ);
#----------------------------------------------------------------------------

# ---
# #### Plot only phase function for I (f₁₁) and the I -> Q transition in the phase matrix (f₁₂) for the Stokes Vector [I,Q,U,V]


p1 = plot(μ, f₁₁, yscale=:log10, title="f₁₁")
p2 = plot(μ, f₁₂ ./ f₁₁,  title="f₁₂/f₁₁")

plot(p1, p2, layout=(2, 1), legend=false)
xlabel!("cos(Θ)")
#----------------------------------------------------------------------------

anim = @animate for r = 0.3:0.2:5
    local size_distribution = LogNormal(log(r), log(σ))
    ## Create the aerosol
    local aero       = make_univariate_aerosol(size_distribution, r_max, n, nᵣ, nᵢ)
    local model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type)
    local aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
    local f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = PhaseFunction.reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ);
    ## @show f₁₁[1]
    p1 = plot(μ, f₁₁, yscale=:log10, title="f₁₁", label="r(μm)=$r")
    ylims!(1e-3, 1e3)
    p2 = plot(μ, f₁₂ ./ f₁₁,  title="f₁₂/f₁₁", label="Q/I")
    ylims!(-1.1, 1.1)
    plot(p1, p2, layout=(2, 1))
end
#----------------------------------------------------------------------------

gif(anim, fps = 5)
#----------------------------------------------------------------------------


#----------------------------------------------------------------------------
