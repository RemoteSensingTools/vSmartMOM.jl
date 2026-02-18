using Revise
using RadiativeTransfer, RadiativeTransfer.vSmartMOM, RadiativeTransfer.Scattering
using Statistics

# Load YAML files into parameter struct
parameters = parameters_from_yaml("test/test_parameters/O2Parameters.yaml");
# Create model struct (precomputes optical properties) from parameters
model      = model_from_parameters(parameters);
m = 0
aerosol_optics = model.aerosol_optics[1];
RaylZ⁺⁺, RaylZ⁻⁺ = Scattering.compute_Z_moments(model.params.polarization_type, Array(model.quad_points.qp_μ), model.greek_rayleigh, m);
AerZ⁺⁺, AerZ⁻⁺   = Scattering.compute_Z_moments(model.params.polarization_type, Array(model.quad_points.qp_μ), aerosol_optics[1].greek_coefs, m)

iz = 5
τRayl = model.τ_rayl[1][iz];
τAer  = model.τ_aer[1][iz]
ϖRayl = 1.0

FT = Float64
τ = FT(0)
ϖ = FT(0)
A = FT(0)
Z⁺⁺ = similar(RaylZ⁺⁺); 
Z⁻⁺ = similar(RaylZ⁺⁺);
if (τRayl + sum(τAer)) < eps(FT)
    fill!(Z⁺⁺, 0); fill!(Z⁻⁺, 0);
    return FT(0), FT(1), Z⁺⁺, Z⁻⁺
end

τ += τRayl
ϖ += τRayl * ϖRayl
A += τRayl * ϖRayl
Z⁺⁺ = τRayl * ϖRayl * RaylZ⁺⁺
Z⁻⁺ = τRayl * ϖRayl * RaylZ⁻⁺
for i = 1:1
    #@show τ, ϖ , A, τAer[i]
    τ   += τAer
    ϖ   += τAer * aerosol_optics[i].ω̃
    A   += τAer * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ)
    Z⁺⁺ += τAer * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ) * AerZ⁺⁺[:,:]
    Z⁻⁺ += τAer * aerosol_optics[i].ω̃ * (1 - aerosol_optics[i].fᵗ) * AerZ⁻⁺[:,:]
    #@show τ, ϖ , A
end

Z⁺⁺ /= A
Z⁻⁺ /= A
A /= ϖ
ϖ /= τ

# Rescaling composite SSPs according to Eqs. A.3 of Sanghavi et al. (2013) or Eqs.(8) of Sanghavi & Stephens (2015)
#@show τRayl, τ,A,  ϖ
τ *= (FT(1) - (FT(1) - A) * ϖ)
ϖ *= A / (FT(1) - (FT(1) - A) * ϖ)#Suniti
#@show τRayl, τ
fscattRayl = τRayl/τ

i = 1;
Ray  = RadiativeTransfer.vSmartMOM.CoreScatteringOpticalProperties(τRayl,ϖRayl,RaylZ⁺⁺, RaylZ⁻⁺ )

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    @unpack fᵗ, ω̃ = aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    RadiativeTransfer.vSmartMOM.CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
end

Aero = createAero(τAer, aerosol_optics[1], AerZ⁺⁺, AerZ⁻⁺)
combo = Ray + Aero