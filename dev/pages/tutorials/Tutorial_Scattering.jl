# # Scattering: Mie Phase Function Tutorial
# 
# ### Introduction
# In the following tutorial, we will walk through how to compute phase function matrices from 
# Mie theory.

# Along the way, we will also calculate greek coefficients (which are essentially the 
# legendre decomposition of the phase matrix components), which will be required for radiative 
# transfer calculations. 

# ---

# ### Performing the phase function calculations

# First let's use the required packages

using vSmartMOM.Scattering
using Distributions
using Plots
#----------------------------------------------------------------------------

# Now, we define the aerosol size distribution and properties. We only support univariate
# aerosols for now, but will add multivariate aerosols soon. 

μ     = 0.3;             ## Log-normal median radius [μm]
σ     = 2.0;             ## Log-normal stddev of radius
nᵣ    = 1.3;             ## Real part of refractive index
nᵢ    = 0.0;             ## Imag part of refractive index (sign changed, use only + here)
r_max = 30.0;            ## Maximum radius [μm]
nquad_radius  = 2500;    ## Number of quadrature points for integrating of size dist.

## Create a Size Distribution (using Julia's Distributions package)
size_distribution = LogNormal(log(μ), log(σ));

#----------------------------------------------------------------------------

# Create the aerosol
aero = Aerosol(size_distribution, nᵣ, nᵢ)
#----------------------------------------------------------------------------

# Now, we define some scattering and truncation properties: 

λ = 0.55;                             ## Incident wavelength [μm]
polarization_type = Stokes_IQUV();    ## Polarization type 

## Trunction length for legendre terms
l_max = 20;

## Exclusion angle for forward peak (in fitting procedure) `[degrees]`
Δ_angle = 2;

## Truncation type
truncation_type   = δBGE(l_max, Δ_angle);      
#----------------------------------------------------------------------------

## Create a Mie model, using the Siewert method NAI2
model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius);
#----------------------------------------------------------------------------

## Compute aerosol optical properties:
aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
#----------------------------------------------------------------------------

# ### Plotting the greek coefficients
# (These are basically giving us the legendre decomposition of the phase matrix components)

using Parameters
@unpack α,β,γ,δ,ϵ,ζ = aerosol_optics_NAI2.greek_coefs;
p1 = plot(α,  title="α");
p2 = plot(β,  title="β");
p3 = plot(γ,  title="γ");
p4 = plot(δ,  title="δ");
p5 = plot(ϵ,  title="ϵ");
p6 = plot(ζ,  title="ζ");
plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), legend=false)
xlims!(0,100)
#----------------------------------------------------------------------------

# ### Reconstructing the Phase Functions from the greek coefficients

using FastGaussQuadrature
μ_quad, w_μ = gausslegendre(1000)
scattering_matrix = Scattering.reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ_quad);
@unpack f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = scattering_matrix
#----------------------------------------------------------------------------

# Plot only phase function for I (f₁₁) and the I -> Q transition in the phase matrix (f₁₂) for the Stokes Vector [I,Q,U,V]

p1 = plot(μ_quad, f₁₁, yscale=:log10, title="f₁₁")
p2 = plot(μ_quad, f₁₂ ./ f₁₁,  title="f₁₂/f₁₁")

plot(p1, p2, layout=(2, 1), legend=false)
xlabel!("cos(Θ)")
#----------------------------------------------------------------------------

p1 = plot([acos.(μ_quad); -reverse(acos.(μ_quad))], log10.([f₁₁ ; reverse(f₁₁)]), proj=:polar, yscale=:log10, title="f₁₁", lims=(-3,4.2), yaxis=false)
p2 = plot([acos.(μ_quad); -reverse(acos.(μ_quad))], [f₁₂ ./ f₁₁ ; reverse(f₁₂ ./ f₁₁)], proj=:polar, title="f₁₂/f₁₁", lims=(0,0.6), yaxis=false)

plot(p1, p2, layout=(1, 2), legend=false)
xlabel!("Θ")

#----------------------------------------------------------------------------

anim = @animate for r = 0.03:0.05:4.3
    @show r
    local size_distribution = LogNormal(log(r), log(σ))
    ## Create the aerosol
    local aero       = Aerosol(size_distribution, nᵣ, nᵢ)
    local model_NAI2 = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type, r_max, nquad_radius)
    local aerosol_optics_NAI2 = compute_aerosol_optical_properties(model_NAI2);
    local scattering_matrix = Scattering.reconstruct_phase(aerosol_optics_NAI2.greek_coefs, μ_quad);
    @unpack f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ = scattering_matrix
    ## @show f₁₁[1]

    p1 = plot(μ_quad, f₁₁, yscale=:log10, title="f₁₁", label="r(μm)=$r")
    ylims!(1e-3, 1e3)
    p2 = plot(μ_quad, f₁₂ ./ f₁₁,  title="f₁₂/f₁₁", label="Q/I")
    ylims!(-1.1, 1.1)
    plot(p1, p2, layout=(2, 1))

    p3 = plot([acos.(μ_quad); -reverse(acos.(μ_quad))], log10.([f₁₁ ; reverse(f₁₁)]), proj=:polar, yscale=:log10, title="f₁₁, r = $(r)μm", lims=(-3,4.2), yaxis=false)
    p4 = plot([acos.(μ_quad); -reverse(acos.(μ_quad))], [f₁₂ ./ f₁₁ ; reverse(f₁₂ ./ f₁₁)], proj=:polar, title="f₁₂/f₁₁, r = $(r)μm", lims=(0,1.0), yaxis=false)

    plot(plot(p1, p2, layout=(2,1), legend=false), plot(p3, p4, layout=(1,2), legend=false), layout=(2, 1), legend=false)

    plot!(size=(600,800))

end
gif(anim, fps = 5)
