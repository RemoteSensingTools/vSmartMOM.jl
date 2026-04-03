const c₂                 = 1.4387769
const cMassMol           = 1.66053873e-27 #grams per molecule for unit molec. mass
const cSqrtLn2divSqrtPi  = 0.469718639319144059835 #√(ln2/π)
const cLn2               = 0.6931471805599 #ln2
const cSqrtLn2           = 0.8325546111577 #√(ln2)
const cSqrt2Ln2          = 1.1774100225 #√(2ln2)
const cc_                = 2.99792458e8 #speed of light [m/s]
const cBolts_            = 1.3806503e-23 #Boltzmann const. [J/K]
const p_ref              = 1013.25  # reference pressure [hPa]
const t_ref              = 296.0    # reference temperature [K]
const nm_per_cm          = 1.0e7

function get_n₀_n₁(ieJ₁⁺, Δ)
    nSpec = size(ieJ₁⁺, 3)
    n₁_start = max(1, 1 - Δ)
    n₁_end   = min(nSpec, nSpec - Δ)
    n₁ = n₁_start:n₁_end
    n₀ = (n₁_start + Δ):(n₁_end + Δ)
    return n₀, n₁
end
# Currently assuming same T for all vertical atmospheric layers (so that a uniform Raman wavelength grid can be assumed for rt_interactions)
function getRamanAtmoConstants(ν̃::AbstractFloat, T::AbstractFloat)
    ν̃, T = promote(ν̃, T)
    n2 = InelasticScattering.getMolecularConstants(InelasticScattering.N₂(), (0.8));
    compute_effective_coefficents!(ν̃, T, n2)
    compute_energy_levels!(n2)
    compute_σ_Rayl_coeff!(n2)
    compute_σ_Rayl_VibRaman_coeff_hires!(T, n2)
    compute_σ_VibRaman_coeff!(T, n2)
    compute_σ_RoVibRaman_coeff!(T, n2)

    o2 = InelasticScattering.getMolecularConstants(InelasticScattering.O₂(), (0.2));
    compute_effective_coefficents!(ν̃, T, o2)
    compute_energy_levels!(o2)
    compute_σ_Rayl_coeff!(o2)
    compute_σ_Rayl_VibRaman_coeff_hires!(T, o2)
    compute_σ_VibRaman_coeff!(T, o2)
    compute_σ_RoVibRaman_coeff!(T, o2)

    return n2,o2
end

function getRamanAtmoConstants(ν̃::AbstractFloat, T::AbstractFloat, vmr_n2::AbstractFloat, vmr_o2::AbstractFloat)
    ν̃, T, vmr_n2, vmr_o2 = promote(ν̃, T, vmr_n2, vmr_o2)
    @assert 0<=vmr_n2+vmr_o2<=1
    n2 = InelasticScattering.getMolecularConstants(InelasticScattering.N₂(), vmr_n2);
    compute_effective_coefficents!(ν̃, T, n2)
    compute_energy_levels!(n2)
    compute_σ_Rayl_coeff!(n2)
    compute_σ_Rayl_VibRaman_coeff_hires!(T, n2)
    compute_σ_VibRaman_coeff!(T, n2)
    compute_σ_RoVibRaman_coeff!(T, n2)

    o2 = InelasticScattering.getMolecularConstants(InelasticScattering.O₂(), vmr_o2);
    compute_effective_coefficents!(ν̃, T, o2)
    compute_energy_levels!(o2)
    compute_σ_Rayl_coeff!(o2)
    compute_σ_Rayl_VibRaman_coeff_hires!(T, o2)
    compute_σ_VibRaman_coeff!(T, o2)
    compute_σ_RoVibRaman_coeff!(T, o2)

    return n2,o2
end

function compute_ϖ_Cabannes(RS_type::noRS, depol, λ₀)
    RS_type.ϖ_Cabannes = 1.0;
    return RS_type.ϖ_Cabannes;
end

#=function compute_ϖ_Cabannes(
            RS_type::Union{RRS, VS_0to1, VS_1to0, RRS_plus, VS_0to1_plus, VS_1to0_plus}, 
            depol, λ₀)
    ν₀ = 1e7/λ₀;

    
    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    ϖ_Cabannes = (4-3*depol)/(4+2*depol)
    #ϖ_Cabannes = σ_elastic/(σ_RRS+σ_elastic);
    return ϖ_Cabannes;
end
=#

function compute_ϖ_Cabannes(
            RS_type::Union{RRS, RRS_plus}, 
            λ₀)
    ν₀ = 1e7/λ₀;

    σ_Rayl =  RS_type.n2.vmr * RS_type.n2.effCoeff.σ_Rayl_coeff + 
                RS_type.o2.vmr * RS_type.o2.effCoeff.σ_Rayl_coeff 
    σ_Rayl *= ν₀^4 # I had initially assumed this to be the purely elastic (Cabannes) scattering cross-section. However, it is actually the Rayleigh cross-section, which becomes purely elastic for solids and liquids due to the absence of a rotational component. The Rayleigh cross-section is the sum of purely elastic (Cabannes) and inelastic (rotational Raman) scattering cross-sections.
    
    σ_RRS =  RS_type.n2.vmr * 
            ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' *
            RS_type.n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += RS_type.n2.vmr * 
            ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * 
            RS_type.n2.effCoeff.σ_RoRaman_coeff_JtoJm2

    σ_RRS += RS_type.o2.vmr * 
            ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * 
            RS_type.o2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += RS_type.o2.vmr * 
            ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * 
            RS_type.o2.effCoeff.σ_RoRaman_coeff_JtoJm2

    #σ_RVRS =  RS_type.n2.vmr * 
    #        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * 
    #        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    #σ_RVRS += RS_type.n2.vmr * 
    #        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * 
    #        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    #σ_RVRS += RS_type.n2.vmr * 
    #        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * 
    #        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    #σ_RVRS += RS_type.n2.vmr * 
    #        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * 
    #        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2
    
    #σ_RVRS += RS_type.o2.vmr * 
    #        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * 
    #        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    #σ_RVRS += RS_type.o2.vmr * 
    #        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * 
    #        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    #σ_RVRS += RS_type.o2.vmr * 
    #        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * 
    #        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    #σ_RVRS += RS_type.o2.vmr * 
    #        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * 
    #        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2

    #σ_VRS =  RS_type.n2.vmr * 
    #        ((ν₀.+RS_type.n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * 
    #        RS_type.n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    #σ_VRS += RS_type.n2.vmr * 
    #        ((ν₀.+RS_type.n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * 
    #        RS_type.n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    #σ_VRS += RS_type.o2.vmr * 
    #        ((ν₀.+RS_type.o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * 
    #        RS_type.o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    #σ_VRS += RS_type.o2.vmr * 
    #        ((ν₀.+RS_type.o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * 
    #        RS_type.o2.effCoeff.σ_VibRaman_coeff_1to0_hires    

    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    #ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    ϖ_Cabannes = 1.0 - σ_RRS/σ_Rayl;
    return ϖ_Cabannes;
end

function compute_ϖ_Cabannes(
    RS_type::Union{VS_0to1_plus, VS_1to0_plus}, 
    λ₀)
    ν₀ = 1e7/λ₀;

    σ_Rayl =  RS_type.n2.vmr * RS_type.n2.effCoeff.σ_Rayl_coeff + 
            RS_type.o2.vmr * RS_type.o2.effCoeff.σ_Rayl_coeff 
    σ_Rayl *= ν₀^4

    #=
    σ_RRS =  RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' *
        RS_type.n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * 
        RS_type.n2.effCoeff.σ_RoRaman_coeff_JtoJm2

    σ_RRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * 
        RS_type.o2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * 
        RS_type.o2.effCoeff.σ_RoRaman_coeff_JtoJm2
    =#
    σ_RVRS =  RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * 
        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    σ_RVRS += RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * 
        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    σ_RVRS += RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * 
        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    σ_RVRS += RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * 
        RS_type.n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2

    σ_RVRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * 
        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    σ_RVRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * 
        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    σ_RVRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * 
        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    σ_RVRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * 
        RS_type.o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2

    σ_VRS =  RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * 
        RS_type.n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    σ_VRS += RS_type.n2.vmr * 
        ((ν₀.+RS_type.n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * 
        RS_type.n2.effCoeff.σ_VibRaman_coeff_1to0_hires

    σ_VRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * 
        RS_type.o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    σ_VRS += RS_type.o2.vmr * 
        ((ν₀.+RS_type.o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * 
        RS_type.o2.effCoeff.σ_VibRaman_coeff_1to0_hires    

    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    #ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    #ϖ_Cabannes = σ_elastic/(σ_RRS+σ_elastic);
    ϖ_Cabannes = ϖ_Cabannes = 1.0 - (σ_VRS+σ_RVRS)/σ_Rayl; #σ_elastic/(σ_VRS+σ_RVRS+σ_elastic);
    return ϖ_Cabannes;
end

function compute_ϖ_Cabannes(λ₀, n2, o2)
    ν₀ = 1e7/λ₀;

    σ_Rayl =  n2.vmr * n2.effCoeff.σ_Rayl_coeff + o2.vmr * o2.effCoeff.σ_Rayl_coeff 
    σ_Rayl *= ν₀^4

    σ_RRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJm2

    σ_RRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * o2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * o2.effCoeff.σ_RoRaman_coeff_JtoJm2

    #σ_RVRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    #σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    #σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    #σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2
    #σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    #σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    #σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    #σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2

    #σ_VRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    #σ_VRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    #σ_VRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    #σ_VRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * o2.effCoeff.σ_VibRaman_coeff_1to0_hires    

    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    #ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    ϖ_Cabannes = 1.0 - σ_RRS/σ_Rayl;
    #@show σ_elastic
    #@show σ_RRS+σ_elastic
    return ϖ_Cabannes;
end

function compute_ϖ_Cabannes_VS(λ₀, n2, o2)
    ν₀ = 1e7/λ₀;

    σ_Rayl =  n2.vmr * n2.effCoeff.σ_Rayl_coeff + o2.vmr * o2.effCoeff.σ_Rayl_coeff 
    σ_Rayl *= ν₀^4

    #=
    σ_RRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJm2

    σ_RRS += o2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += o2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * n2.effCoeff.σ_RoRaman_coeff_JtoJm2
    =#
    
    σ_RVRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    σ_RVRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    σ_RVRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2

    σ_VRS =  n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    σ_VRS += n2.vmr * ((ν₀.+n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    σ_VRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    σ_VRS += o2.vmr * ((ν₀.+o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * o2.effCoeff.σ_VibRaman_coeff_1to0_hires    

    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    #ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    ϖ_Cabannes_VS = 1.0 - (σ_RVRS+σ_VRS)/σ_Rayl;
    #@show σ_elastic
    #@show σ_RVRS+σ_VRS+σ_elastic
    return ϖ_Cabannes_VS;
end

function compute_ϖ_Cabannes(λ₀, mol)
    ν₀ = 1e7/λ₀;

    σ_Rayl =  mol.effCoeff.σ_Rayl_coeff
    σ_Rayl *= ν₀^4

    σ_RRS =  ((ν₀.+mol.effCoeff.Δν̃_RoRaman_coeff_JtoJp2).^4)' * mol.effCoeff.σ_RoRaman_coeff_JtoJp2
    σ_RRS += ((ν₀.+mol.effCoeff.Δν̃_RoRaman_coeff_JtoJm2).^4)' * mol.effCoeff.σ_RoRaman_coeff_JtoJm2

    #σ_RVRS =  ((ν₀.+mol.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2).^4)' * mol.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2
    #σ_RVRS += ((ν₀.+mol.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2).^4)' * mol.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2
    #σ_RVRS += ((ν₀.+mol.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2).^4)' * mol.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2
    #σ_RVRS += ((ν₀.+mol.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2).^4)' * mol.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2
    
    #σ_VRS =  ((ν₀.+mol.effCoeff.Δν̃_VibRaman_coeff_0to1_hires).^4)' * mol.effCoeff.σ_VibRaman_coeff_0to1_hires
    #σ_VRS += ((ν₀.+mol.effCoeff.Δν̃_VibRaman_coeff_1to0_hires).^4)' * mol.effCoeff.σ_VibRaman_coeff_1to0_hires


    #RS_type.ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    #ϖ_Cabannes = σ_elastic/(σ_VRS+σ_RVRS+σ_RRS+σ_elastic);
    ϖ_Cabannes = 1.0 - σ_RRS/σ_Rayl;
    #@show σ_elastic
    #@show σ_RRS+σ_elastic
    return ϖ_Cabannes;
end

# New: γ_air_Rayleigh is computed from the effective Rayleigh depolarization ratios of N₂ and O₂, which are then combined to give the effective depolarization ratio of air. 
function compute_γ_air_Rayleigh!(λ₀::FT, RS_type::Union{RRS, VS_0to1, VS_1to0, RRS_plus, VS_0to1_plus, VS_1to0_plus}) where FT

    γN₂ = RS_type.n2.effCoeff.γ_C_Rayl
    σ0N₂ = RS_type.n2.effCoeff.σ_Rayl_coeff
    σN₂ = σ0N₂ * (3-4γN₂)/(1+2γN₂) # 128π^5 * α̅^2 
    VMR_N₂ = RS_type.n2.vmr

    γO₂ = RS_type.o2.effCoeff.γ_C_Rayl
    σ0O₂ = RS_type.o2.effCoeff.σ_Rayl_coeff
    σO₂ = σ0O₂ * (3-4γO₂)/(1+2γO₂) # 128π^5 * α̅^2 
    VMR_O₂ = RS_type.o2.vmr

    tmp1 = σN₂ * VMR_N₂ + σO₂ * VMR_O₂
    tmp2 = (σN₂ * VMR_N₂ * γN₂)/(3-4γN₂) + (σO₂ * VMR_O₂ * γO₂)/(3-4γO₂)
    
    γ_air_Rayleigh = 3/(4 + tmp1 / tmp2)
    σ_air_Rayleigh = (σ0N₂ * VMR_N₂ + σ0O₂ * VMR_O₂)*(nm_per_cm/λ₀)^4/(VMR_N₂ + VMR_O₂) # cm^{2}/molecule

    return γ_air_Rayleigh, σ_air_Rayleigh
end

# New: γ_air_Cabannes is computed by assuming the same form of the scattering equation (Is = σ*((1+2γ)/(3-4γ))*P*Ii*dz, see Eq.14 of Sanghavi(2022) for exact form) as for Rayleigh scattering      
function compute_γ_air_Cabannes!(λ₀::FT, RS_type::Union{RRS, VS_0to1, VS_1to0, RRS_plus, VS_0to1_plus, VS_1to0_plus}) where FT

    γN₂ = compute_γ_mol_Cabannes!(λ₀, RS_type.n2)[2]
    ϖN₂ = compute_ϖ_Cabannes(λ₀, RS_type.n2)
    σ0N₂ = RS_type.n2.effCoeff.σ_Rayl_coeff 
    σN₂ = ϖN₂ * σ0N₂ * (3-4γN₂)/(1+2γN₂) # 128π^5 * α̅^2 
    VMR_N₂ = RS_type.n2.vmr

    γO₂ = compute_γ_mol_Cabannes!(λ₀, RS_type.o2)[2]
    ϖO₂ = compute_ϖ_Cabannes(λ₀, RS_type.o2)
    σ0O₂ = RS_type.o2.effCoeff.σ_Rayl_coeff
    σO₂ = ϖO₂ * σ0O₂ * (3-4γO₂)/(1+2γO₂) # 128π^5 * α̅^2 
    VMR_O₂ = RS_type.o2.vmr

    tmp1 = σN₂ * VMR_N₂ + σO₂ * VMR_O₂
    tmp2 = (σN₂ * VMR_N₂ * γN₂)/(3-4γN₂) + (σO₂ * VMR_O₂ * γO₂)/(3-4γO₂)

    γ_air_Cabannes = 3/(4 + tmp1 / tmp2)

    ϖ_air_Cabannes = (ϖN₂ * σ0N₂ * VMR_N₂ + ϖO₂ * σ0O₂ * VMR_O₂) / (σ0N₂ * VMR_N₂ + σ0O₂ * VMR_O₂)
    
    return γ_air_Cabannes, ϖ_air_Cabannes;
end

# New, revised
function compute_γ_air_Cabannes!(λ₀::FT,
                        n2, o2) where FT
    γN₂ = compute_γ_mol_Cabannes!(λ₀, n2)[2]
    ϖN₂ = compute_ϖ_Cabannes(λ₀, n2)
    σ0N₂ = n2.effCoeff.σ_Rayl_coeff 
    σN₂ = ϖN₂ * σ0N₂ * (3-4γN₂)/(1+2γN₂) # 128π^5 * α̅^2 
    VMR_N₂ = n2.vmr

    γO₂ = compute_γ_mol_Cabannes!(λ₀, o2)[2]
    ϖO₂ = compute_ϖ_Cabannes(λ₀, o2)
    σ0O₂ = o2.effCoeff.σ_Rayl_coeff
    σO₂ = ϖO₂ * σ0O₂ * (3-4γO₂)/(1+2γO₂) # 128π^5 * α̅^2 
    VMR_O₂ = o2.vmr

    tmp1 = σN₂ * VMR_N₂ + σO₂ * VMR_O₂
    tmp2 = (σN₂ * VMR_N₂ * γN₂)/(3-4γN₂) + (σO₂ * VMR_O₂ * γO₂)/(3-4γO₂)

    γ_air_Cabannes = 3/(4 + tmp1 / tmp2)

    ϖ_air_Cabannes = (ϖN₂ * σ0N₂ * VMR_N₂ + ϖO₂ * σ0O₂ * VMR_O₂) / (σ0N₂ * VMR_N₂ + σ0O₂ * VMR_O₂)
    
    return γ_air_Cabannes, ϖ_air_Cabannes;
end

# New, revised
function compute_γ_air_Rayleigh!(λ₀::FT, n2, o2) where FT
    γN₂ = n2.effCoeff.γ_C_Rayl
    σ0N₂ = n2.effCoeff.σ_Rayl_coeff
    σN₂ = σ0N₂ * (3-4γN₂)/(1+2γN₂) # 128π^5 * α̅^2 
    VMR_N₂ = n2.vmr

    γO₂ = o2.effCoeff.γ_C_Rayl
    σ0O₂ = o2.effCoeff.σ_Rayl_coeff
    σO₂ = σ0O₂ * (3-4γO₂)/(1+2γO₂) # 128π^5 * α̅^2 
    VMR_O₂ = o2.vmr

    tmp1 = σN₂ * VMR_N₂ + σO₂ * VMR_O₂
    tmp2 = (σN₂ * VMR_N₂ * γN₂)/(3-4γN₂) + (σO₂ * VMR_O₂ * γO₂)/(3-4γO₂)
    
    γ_air_Rayleigh = 3/(4 + tmp1 / tmp2)
    σ_air_Rayleigh = (σ0N₂ * VMR_N₂ + σ0O₂ * VMR_O₂)*(nm_per_cm/λ₀)^4/(VMR_N₂ + VMR_O₂) # cm^{2}/molecule

    return γ_air_Rayleigh, σ_air_Rayleigh
end

# New: In the old version, mol.effCoeff.γ_C_Rayl was assumed to be γ_mol_Cabannes, but it actually was indeed γ_mol_Rayleigh. This assumption has now been corrected, and the correct γ_mol_Cabannes is computed in the function below. 
function compute_γ_mol_Cabannes!(λ₀::FT, mol) where FT
    ν̃ =  1.e7/λ₀
    effT  =  300.; #K assumed constant for Earth atmospheres

    #n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
    #greek_raman = get_greek_raman(RS_type, n2, o2);
    ϖ_Cabannes = compute_ϖ_Cabannes(λ₀, mol)
    γ_mol_Rayleigh = mol.effCoeff.γ_C_Rayl

    tmp1 = 1+2*γ_mol_Rayleigh
    tmp2 = 2+3*ϖ_Cabannes
    tmp3 = 1-ϖ_Cabannes
    tmpN = tmp1*tmp2-5
    tmpD = tmp1*tmp3+5
    γ_mol_Cabannes = 0.5*tmpN/tmpD


    #@show γ_mol_Rayleigh, tmp_chk
    return ϖ_Cabannes, γ_mol_Cabannes, γ_mol_Rayleigh
end

# for VS
function compute_γ_air_Rayleigh_VS!(λ₀::FT) where FT
    ν̃ =  1.e7/λ₀
    effT  =  300.; #K assumed constant for Earth atmospheres

    n2,o2 = InelasticScattering.getRamanAtmoConstants(ν̃,effT);
    #greek_raman = get_greek_raman(RS_type, n2, o2);
    ϖ_Cabannes_VS = compute_ϖ_Cabannes_VS(λ₀, n2, o2)
    γ_air_Cabannes = compute_γ_air_Cabannes!(n2, o2)

    tmp = (1/ϖ_Cabannes_VS)*((1+2γ_air_Cabannes)/(3-4γ_air_Cabannes))
    γ_air_Rayleigh = 0.5*((3*tmp-1)/(2*tmp+1))

    return ϖ_Cabannes_VS, γ_air_Cabannes, γ_air_Rayleigh
end

# Note: ν stands for wavenumber in the following (NOT frequency)
function apply_lineshape!(Δνᵢ, σᵢ,  # discrete transitions
                λ₀,                 # incident  wavelength [nm]
                # Note:  λ₀ can either be used to denote an individual 
                # wavelength or as a representative (e.g. central) 
                # wavelength for a specific band
                Δν_out,             # Output grid (equidistant)
                σ_out,              # σ at output grid
                pressure::Real,     # actual pressure [hPa]
                temperature::Real,  # Temperature (K)
                molMass; 
                wavelength_flag::Bool=false)
    σ_out .= 0;
    # Notify user of wavelength grid
    if (wavelength_flag)
        @info """
        Note: Rayleigh/Raman Cross-section reported to wavelength grid (nm)
        """  maxlog = 5
    end

    # Max min range (ignoring wing cutoff here)
    grid_max = maximum(Δν_out) 
    grid_min = minimum(Δν_out) 

    # Interpolators from grid bounds to index values
    grid_idx_interp_low  = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=1)
    grid_idx_interp_high = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=length(Δν_out))

    S_sum=0.0
    # Loop through all transition lines:
    for j in eachindex(Δνᵢ)
        # Test that this ν lies within the grid
        if grid_min < Δνᵢ[j] < grid_max
            
            ν = Δνᵢ[j] + nm_per_cm/λ₀#13500.0 #Dummy for now #Suniti

            # Compute Doppler HWHM, ν still needs to be supplied, @Suniti?:
            γ_d = ((cSqrt2Ln2 / cc_) * sqrt(cBolts_ / cMassMol) * sqrt(temperature) * ν / sqrt(molMass))

            # line intensity 
            S = σᵢ[j] *  ν^4 #Suniti
            S_sum += S
            #@show γ_d, Δνᵢ[j], S, ν, λ₀, cMassMol

            wing_cutoff = 2γ_d 

            # Calculate index range that this transition impacts
            ind_start = Int64(floor(grid_idx_interp_low(Δνᵢ[j] - wing_cutoff)))
            ind_stop  = Int64(ceil(grid_idx_interp_high(Δνᵢ[j] + wing_cutoff)))
            
            # Create views from the result and grid arrays
            result_view   = view(σ_out,  ind_start:ind_stop);
            grid_view     = view(Δν_out, ind_start:ind_stop);

            # Just Doppler broadening for now, can be modified with any line-shape later (using a kernel ideally)
            for I in eachindex(grid_view)
                # If we undersample the line-width, we have to make sure the integral is conserved (almost), TBD
                @inbounds result_view[I] += S * cSqrtLn2divSqrtPi * exp(-cLn2 * (((grid_view[I]) - Δνᵢ[j]) / γ_d)^2) / γ_d
                #@show grid_view[I], result_view[I]
            end
        end
    end
    #dν = Δν_out[2]-Δν_out[1] 
    nothing  
    #@show S_sum, sum(σ_out)*dν
end

function apply_gridlines!(Δνᵢ, σᵢ,  # discrete transitions
    λ₀,                 # incident  wavelength [nm]
    # Note:  λ₀ can either be used to denote an individual 
    # wavelength or as a representative (e.g. central) 
    # wavelength for a specific band
    ν_in,             # Output grid (equidistant)
    σ_out;              # σ at output grid
    
    #pressure::Real,     # actual pressure [hPa]
    #temperature::Real,  # Temperature (K)
    #molMass; 
    wavelength_flag::Bool=false)    
    σ_out .= 0;
    # Notify user of wavelength grid
    if (wavelength_flag)
        @info """
        Note: Rayleigh/Raman Cross-section reported to wavelength grid (nm)
        """  maxlog = 5
    end
    Δν_in = ν_in .- nm_per_cm/λ₀

    # Max min range (ignoring wing cutoff here)
    grid_max = maximum(Δν_in) 
    grid_min = minimum(Δν_in) 

    # Interpolators from grid bounds to index values
    #grid_idx_interp_low  = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=1) 
    #grid_idx_interp_high = LinearInterpolation(Δν_out, 1:1:length(Δν_out), extrapolation_bc=length(Δν_out))
    S_sum = 0.0
    # Loop through all transition lines:
    for j in eachindex(Δνᵢ)
        #@show(grid_min, Δνᵢ[j], grid_max)
        # Test that this ν lies within the grid
        if grid_min < Δνᵢ[j] < grid_max
            
            ν = Δνᵢ[j] + nm_per_cm/λ₀#13500.0 #Dummy for now #Suniti

            # Compute Doppler HWHM, ν still needs to be supplied, @Suniti?:
            #γ_d = ((cSqrt2Ln2 / cc_) * sqrt(cBolts_ / cMassMol) * sqrt(temperature) * ν / sqrt(molMass))
            # line intensity 
            S = σᵢ[j] *  ν^4 #Suniti
            #@show γ_d, Δνᵢ[j], S, ν, λ₀, cMassMol
            S_sum += S
            #wing_cutoff = 2γ_d 
            i=argmin(abs.(Δνᵢ[j].-Δν_in))
            #@show i, Δνᵢ[j]-Δν_in[i]    
            if Δν_in[i]<Δνᵢ[j]
                ind_start = i
                ind_stop  = i+1
            else
                ind_start = i-1
                ind_stop  = i
            end         
            # Calculate index range that this transition impacts
            #ind_start = Int64(floor(grid_idx_interp_low(Δνᵢ[j] - wing_cutoff)))
            #ind_stop  = Int64(ceil(grid_idx_interp_high(Δνᵢ[j] + wing_cutoff)))

            # Create views from the result and grid arrays
            result_view   = view(σ_out,  ind_start:ind_stop);
            grid_view     = view(Δν_in, ind_start:ind_stop);

            # Just Doppler broadening for now, can be modified with any line-shape later (using a kernel ideally)
            for I in eachindex(grid_view)
                # If we undersample the line-width, we have to make sure the integral is conserved (almost), TBD
                @inbounds result_view[I] += S/2.# * cSqrtLn2divSqrtPi * exp(-cLn2 * (((grid_view[I]) - Δνᵢ[j]) / γ_d)^2) / γ_d
                #@show grid_view[I], result_view[I]
            end
        end
    end
    #@show S_sum, sum(σ_out)
    nothing
end


#function compute_optical_Rayl!(grid_out,atmo_σ_Rayl, λ₀, n2, o2)
function compute_optical_Rayl(λ₀, n2, o2)
    atmo_σ_Rayl = n2.vmr * n2.effCoeff.σ_Rayl_coeff #σ_out #cross section in cm^2
    atmo_σ_Rayl += o2.vmr * o2.effCoeff.σ_Rayl_coeff #cross section in cm^2
    atmo_σ_Rayl *= (nm_per_cm/λ₀)^4
    #plot(1.e7/λ₀ .+ grid_out,atmo_σ_Rayl*1.e40)
    return atmo_σ_Rayl;
end

function compute_optical_RS!(RS_type::Union{RRS, RRS_plus}, grid_in, λ₀, n2, o2)
    #plotly()
    # grid_in is a uniform wavenumber grid covering the entire band spectrum 
    # TMP: grid_in = nm_per_cm/λ₀.+collect(-250:0.002:250) #this is a wavenumber grid
    # get_greek_raman!(RS_type, n2, o2)
    
    grid_in_collected = collect(grid_in)
    σ_out = similar(grid_in_collected);
    atmo_σ_RRS_JtoJp2 = similar(grid_in_collected);
    atmo_σ_RRS_JtoJm2 = similar(grid_in_collected);
    σ_tmp = similar(grid_in_collected);

    # N2
    apply_gridlines!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, n2.effCoeff.σ_RoRaman_coeff_JtoJp2,  λ₀, grid_in_collected, σ_out);
    atmo_σ_RRS_JtoJp2 = n2.vmr * σ_out #cross section in cm^2

    apply_gridlines!(n2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, n2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ₀, grid_in_collected, σ_out);
    atmo_σ_RRS_JtoJm2 = n2.vmr * σ_out #cross section in cm^2
    # O2
    apply_gridlines!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJp2, o2.effCoeff.σ_RoRaman_coeff_JtoJp2, λ₀, grid_in_collected, σ_out);
    atmo_σ_RRS_JtoJp2 += o2.vmr * σ_out #cross section in cm^2

    apply_gridlines!(o2.effCoeff.Δν̃_RoRaman_coeff_JtoJm2, o2.effCoeff.σ_RoRaman_coeff_JtoJm2, λ₀, grid_in_collected, σ_out);
    atmo_σ_RRS_JtoJm2 += o2.vmr * σ_out #cross section in cm^2

    σ_tmp .= atmo_σ_RRS_JtoJm2 .+ atmo_σ_RRS_JtoJp2
    index_ramangrid_out = findall(>(0), σ_tmp)
    atmo_σ_RRS = σ_tmp[index_ramangrid_out]
    if (nm_per_cm/λ₀>grid_in[1] && nm_per_cm/λ₀<grid_in[end])
        index_ramangrid_out .-= argmin(abs.(grid_in_collected .- nm_per_cm/λ₀))
    end 
    #for I in eachindex(atmo_σ_RRS)
    #    @show grid_in[argmin(abs.(grid_in .- nm_per_cm/λ₀))+index_ramangrid_out[I]], index_ramangrid_out[I], atmo_σ_RRS[I]
    #end    
    return index_ramangrid_out, atmo_σ_RRS;
    #plot(grid_out, atmo_σ_RRS_JtoJp2, yscale=:log10)
    #plot(1.e7/λ₀ .+ grid_out, atmo_σ_RRS_plot*1.e40)
end

function compute_optical_RS!(RS_type::Union{VS_0to1, VS_0to1_plus}, grid_in, λ₀, n2, o2)
    #plotly()
    #get_greek_raman(RS_type, n2, o2)
    #compute_ϖ_Cabannes!(RS_type, λ₀, n2, o2)

    #@show n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0], o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0]
    #νᵣ = 0.5*(n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0] + o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires[0])
    
    # TMP: grid_in = nm_per_cm/λ₀ .+ collect((νᵣ-750):0.002:(νᵣ+750))
    σ_out = similar(grid_in);
    #atmo_σ_VRS_0to1 = similar(grid_in);
    #atmo_σ_RVRS_0to1 = similar(grid_in);
    σ_tmpVRS = similar(grid_in);
    σ_tmpRVRS = similar(grid_in);
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2.parent; n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2.parent]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2.parent; n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_tmpRVRS = n2.vmr * σ_out #cross section in cm^2
    xin = n2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = n2.effCoeff.σ_VibRaman_coeff_0to1_hires
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_tmpVRS = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2.parent; o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2.parent]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2.parent; o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_tmpRVRS += o2.vmr * σ_out #cross section in cm^2
    xin = o2.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = o2.effCoeff.σ_VibRaman_coeff_0to1_hires
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_tmpVRS += o2.vmr * σ_out #cross section in cm^2

    index_VRSgrid_out = findall(>(0), σ_tmpVRS)
    atmo_σ_VRS_0to1 = σ_tmpVRS[index_VRSgrid_out]

    index_RVRSgrid_out = findall(>(0), σ_tmpRVRS)
    atmo_σ_RVRS_0to1 = σ_tmpRVRS[index_RVRSgrid_out]
    #plot(grid_out, atmo_σ_RRS_JtoJp2, yscale=:log10)
    #plot(grid_in, σ_tmpRVRS*1.e40)
    #plot!(grid_in, σ_tmpVRS*1.e40)
    return index_VRSgrid_out, atmo_σ_VRS_0to1, index_RVRSgrid_out, atmo_σ_RVRS_0to1;
end


function compute_optical_RS!(RS_type::Union{VS_1to0, VS_1to0_plus}, grid_in, λ₀, n2, o2)
    #plotly()
    get_greek_raman(RS_type, n2, o2)
    compute_ϖ_Cabannes!(RS_type, λ₀, n2, o2)
    #@show n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0], o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0]
    νᵣ = 0.5*(n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0] + o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires[0])
        
    # TMP: grid_in = nm_per_cm/λ₀ + collect((νᵣ-750):0.002:(νᵣ+750))
    σ_out = similar(grid_in);
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2.parent; n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2.parent]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2.parent; n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_RVRStmp = n2.vmr * σ_out #cross section in cm^2

    xin = n2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = n2.effCoeff.σ_VibRaman_coeff_1to0_hires
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_VRStmp = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2.parent; o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2.parent]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2.parent; o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_RVRStmp += o2.vmr * σ_out #cross section in cm^2

    xin = o2.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = o2.effCoeff.σ_VibRaman_coeff_1to0_hires
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_VRStmp += o2.vmr * σ_out #cross section in cm^2

    index_VRSgrid_out = findall(>(0), σ_VRStmp)
    atmo_σ_VRS_1to0 = σ_VRStmp[index_VRSgrid_out]

    index_RVRSgrid_out = findall(>(0), σ_RVRStmp)
    atmo_σ_RVRS_1to0 = σ_RVRStmp[index_RVRSgrid_out]

    return index_VRSgrid_out, atmo_σ_VRS_1to0, index_RVRSgrid_out, atmo_σ_RVRS_1to0;
    #plot(grid_out, atmo_σ_RRS_JtoJp2, yscale=:log10)
    #plot(grid_in, σ_RVRStmp*1.e40)
    #plot!(grid_in, σ_VRStmp*1.e40)
end

#===========For target grids==========#
function compute_optical_RVRS!(RS_type::Union{VS_0to1, VS_0to1_plus}, grid_in, λ₀, n2, o2)
    
    σ_out = similar(grid_in);
    #atmo_σ_VRS_0to1 = similar(grid_in);
    #atmo_σ_RVRS_0to1 = similar(grid_in);
    #σ_tmpVRS = similar(grid_in);
    σ_tmpRVRS = similar(grid_in);
    
    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2.parent; n2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2.parent]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2.parent; n2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_tmpRVRS = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJp2.parent; o2.effCoeff.Δν̃_RoVibRaman_coeff_0to1_JtoJm2.parent]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJp2.parent; o2.effCoeff.σ_RoVibRaman_coeff_0to1_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_tmpRVRS += o2.vmr * σ_out #cross section in cm^2

    index_RVRSgrid_out = findall(>(0), σ_tmpRVRS)
    atmo_σ_RVRS_0to1 = σ_tmpRVRS[index_RVRSgrid_out]

    return index_RVRSgrid_out, atmo_σ_RVRS_0to1;
end


function compute_optical_RVRS!(RS_type::Union{VS_1to0, VS_1to0_plus}, grid_in, λ₀, n2, o2)
    σ_out = similar(grid_in);

    # N2
    xin = [n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2.parent; n2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2.parent]
    yin = [n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2.parent; n2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_RVRStmp = n2.vmr * σ_out #cross section in cm^2

    # O2
    xin = [o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJp2.parent; o2.effCoeff.Δν̃_RoVibRaman_coeff_1to0_JtoJm2.parent]
    yin = [o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJp2.parent; o2.effCoeff.σ_RoVibRaman_coeff_1to0_JtoJm2.parent];
    nz_mask = abs.(xin) .> 0
    apply_gridlines!(xin[nz_mask], yin[nz_mask], λ₀, grid_in, σ_out);
    σ_RVRStmp += o2.vmr * σ_out #cross section in cm^2

    index_RVRSgrid_out = findall(>(0), σ_RVRStmp)
    atmo_σ_RVRS_1to0 = σ_RVRStmp[index_RVRSgrid_out]

    return index_RVRSgrid_out, atmo_σ_RVRS_1to0;
end

function compute_optical_VRS!(RS_type::Union{VS_0to1, VS_0to1_plus}, grid_in, λ₀, mol)
    σ_out = similar(grid_in);
    σ_tmpVRS = similar(grid_in);

    # mol: N2 OR O2
    xin = mol.effCoeff.Δν̃_VibRaman_coeff_0to1_hires
    yin = mol.effCoeff.σ_VibRaman_coeff_0to1_hires
    #apply_lineshape!(xin, yin, λ₀, collect(grid_out), σ_out, 1, 300.0, 40);
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_tmpVRS = mol.vmr * σ_out #cross section in cm^2

    index_VRSgrid_out = findall(>(0), σ_tmpVRS)
    atmo_σ_VRS_0to1 = σ_tmpVRS[index_VRSgrid_out]

    return index_VRSgrid_out, atmo_σ_VRS_0to1;
end

function compute_optical_VRS!(RS_type::Union{VS_1to0, VS_1to0_plus}, grid_in, λ₀, mol)
    σ_out = similar(grid_in);
    σ_tmpVRS = similar(grid_in);

    # mol: N2 OR O2
    xin = mol.effCoeff.Δν̃_VibRaman_coeff_1to0_hires
    yin = mol.effCoeff.σ_VibRaman_coeff_1to0_hires
    apply_gridlines!(xin, yin, λ₀, grid_in, σ_out);
    σ_tmpVRS = mol.vmr * σ_out #cross section in cm^2

    index_VRSgrid_out = findall(>(0), σ_tmpVRS)
    atmo_σ_VRS_1to0 = σ_tmpVRS[index_VRSgrid_out]

    return index_VRSgrid_out, atmo_σ_VRS_1to0;
end


#=======================================#

"""
    $(FUNCTIONNAME)(depol)
Returns the greek coefficients (as [`GreekCoefs`](@ref)) of the Rayleigh phase function given 
depolarization value. 
- `depol` Depolarization (best use 0 as default )
"""
function get_greek_raman!(RS_type::noRS, n2, o2)
    return nothing
end

# the following applies to both rovibrational and rotational Raman scattering (by both N2 and O2)
function get_greek_raman(RS_type::Union{RRS, RRS_plus, VS_0to1, VS_0to1_plus, VS_1to0, VS_1to0_plus}, 
                            n2, o2)
    depol = n2.effCoeff.rho_depol_RotRaman
    FT = eltype(depol)

    # Rayleigh Greek Parameters
    dpl_p = (1 - depol)  / (1 + depol/2)
    #dpl_q = (1 + depol)  / (1 - depol)
    dpl_r = (1 - 2depol) / (1 - depol)
  
    α  =  FT[0.0, 0.0,             3dpl_p]
    β  =  FT[1.0, 0.0,             0.5 * dpl_p]
    γ  =  FT[0.0, 0.0,             dpl_p * sqrt(1.5)] 
    δ  =  FT[0.0, dpl_p * dpl_r * 1.5, 0.0] 
    ϵ  =  FT[0.0, 0.0,             0.0] 
    ζ  =  FT[0.0, 0.0,             0.0]
    return GreekCoefs(α, β, γ, δ, ϵ, ζ);
    #return nothing
end

function get_greek_raman_VS(RS_type::Union{VS_0to1, VS_0to1_plus, VS_1to0, VS_1to0_plus}, 
                            in_molec)
    
    depol = in_molec.effCoeff.rho_depol_VibRaman
    
    FT = eltype(depol)
    # Rayleigh Greek Parameters
    dpl_p = (1 - depol)  / (1 + depol / 2)
    #dpl_q = (1 + depol)  / (1 - depol)
    dpl_r = (1 - 2depol) / (1 - depol)
  
    α  =  FT[0.0, 0.0,             3dpl_p]
    β  =  FT[1.0, 0.0,             0.5 * dpl_p]
    γ  =  FT[0.0, 0.0,             dpl_p * sqrt(1.5)] 
    δ  =  FT[0.0, dpl_p * dpl_r * 1.5, 0.0] 
    ϵ  =  FT[0.0, 0.0,             0.0] 
    ζ  =  FT[0.0, 0.0,             0.0]
    return GreekCoefs(α, β, γ, δ, ϵ, ζ);
    #return nothing
end

function compute_Rayl_depol(n2, o2)
    depol = (n2.vmr*n2.effCoeff.rho_depol_Rayl + o2.vmr*o2.effCoeff.rho_depol_Rayl)/(n2.vmr+o2.vmr)
    return depol
end


function computeRamanZλ!(RS_type::Union{RRS_plus,RRS}, pol_type, qp_μ, m, arr_type)
    RS_type.Z⁺⁺_λ₁λ₀, RS_type.Z⁻⁺_λ₁λ₀ =  Scattering.compute_Z_moments(pol_type, 
                                        qp_μ, 
                                        RS_type.greek_raman, 
                                        m, 
                                        arr_type = arr_type);
    nothing
end

function computeRamanZλ!(RS_type::Union{noRS_plus, noRS}, pol_type, qp_μ, m, arr_type)
    nothing
end

function computeRamanZλ!(RS_type::AbstractRamanType, pol_type, qp_μ, m, arr_type)
    RS_type.Z⁺⁺_λ₁λ₀, RS_type.Z⁻⁺_λ₁λ₀ = Scattering.compute_Z_moments(pol_type, 
                                        qp_μ, 
                                        RS_type.greek_raman, 
                                        m, 
                                        arr_type = arr_type);
    RS_type.Z⁺⁺_λ₁λ₀_VS_n2, RS_type.Z⁻⁺_λ₁λ₀_VS_n2 = 
                    Scattering.compute_Z_moments(pol_type, 
                                            Array(qp_μ), 
                                            RS_type.greek_raman_VS_n2, 
                                            m, 
                                            arr_type = arr_type);
    RS_type.Z⁺⁺_λ₁λ₀_VS_o2, RS_type.Z⁻⁺_λ₁λ₀_VS_o2 = 
                    Scattering.compute_Z_moments(pol_type, 
                                        Array(qp_μ), 
                                        RS_type.greek_raman_VS_o2, 
                                        m, 
                                        arr_type = arr_type);      
    nothing
end




