abstract type AbstractMolecule end

struct N₂ <: AbstractMolecule end
struct O₂ <: AbstractMolecule end
struct H₂ <: AbstractMolecule end

"""
    struct PolTensor{FT}

A struct which provides all polarizability tensor elements 
    relevant to Rayleigh, and rotational+vibrational Raman cross-sections  
    Ref: Buldakov et al.,(1999), Molecular Spectroscopy

    The polarizability at arbitrary wavenumbers [m⁻¹] and temperatures [K] can be computed as 
    
    α̅(2πcν, T) = α̅₀₀(1 + α_b⋅T + α_c⋅T²)/(1-(2πcν/ω₀)²

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct PolarizationTensor{FT}
    "Average Polarizability `[m³]` at an angular frequency of ω₀ `[Hz]` and T=0K"
    α̅₀₀::FT        
    "First derivative of `α̅₀₀` wrt the internuclear distance rₑ `[m²]`"
    α₀₀_prime::FT  
    "Reference frequency `[Hz]`"
    ω₀::FT     
    "Linear T-dependence coefficient `[K⁻¹]`"
    α_b::FT      
    "Quadratic T-dependence coefficient `[K⁻²]`"
    α_c::FT
    "Polarizability anisotropy `[m³]` " 
    γ̅₀₀ ::FT
    "First derivative of `γ̅` with respect to the internuclear distance `rₑ` `[m²]`"
    γ₀₀_prime::FT       
end

"""
    struct DunhamCoefficients{FT}

A struct which provides 
the nuclear spin multiplicity for a (currently diatomic homonuclear) molecule for odd 
and even values of initial rotational quantum states J_initial 

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct NuclearSpinMultiplicity
    "Nuclear Spin Multiplicity (Array of length 2 {odd J_initial, even J_initial})"
    gₛ::AbstractArray{Int}    
end

"""
    struct DunhamCoefficients{FT}

A struct which provides 
the Dunham expansion coefficients to compute the molecular energy levels 
corresponding to different rotational (J) and vibrational (v) states. The overall state 
of the molecule is given by its electronic ground state. #Constants from "Molecular 
Spectra and Molecular Structure IV. Constants of diatomic molecules" by K.P. Huber 
and G. Herzberg, 1978 (https://webbook.nist.gov/cgi/cbook.cgi?ID=C7727379&Mask=1000)

Needs description of Y elements

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct DunhamCoefficients{FT}
    "Dunham Expansion Coefficients (5x5 Matrix)"
    Y::AbstractArray{FT}    
end

Base.@kwdef mutable struct EffectiveCoefficients{FT}
    "Ambient Temperature `[K]` at which following effective coefficients are computed"
    T::FT 
    "Average Polarizability `[m³]` at wavenumber ν and temperature T"
    α̅::FT      #  Assumes a central wavelength effective over a broadband range of wavelengths
    γ̅::FT
    α_prime::FT
    γ_prime::FT
    #"King correction factor, See Table 1 and Eq 4 of Chance & Spurr (1997)"
    #F_King::FT #  See Table 1 and Eq 4 of Chance & Spurr (1997)
    "α/γ"
    ϵ::FT
    #"King correction factor for vibrational Raman scattering, see Tables H, I (CRS 13, D. A. Long (1977)) "
    #F_King_prime::FT
    "α_prime/γ_prime"
    ϵ_prime::FT
    γ_C_Rayl::FT # 3/(45(ϵ)^2+4)
    γ_C_RotRaman::FT # 3/4
    γ_C_VibRaman::FT # 3/(45(ϵ_prime)^2+4)
    γ_C_RoVibRaman::FT # 3/4
    rho_depol_Rayl::FT # 2γ_C_Rayl/(1+γ_C_Rayl)
    rho_depol_RotRaman::FT # 2γ_C_RotRaman/(1+γ_C_RotRaman)
    rho_depol_VibRaman::FT # 2γ_C_VibRaman/(1+γ_C_VibRaman)
    rho_depol_RoVibRaman::FT # 2γ_C_RoVibRaman/(1+γ_C_RoVibRaman)
    "Rayl cross-section σ(ν) = σ_Rayl_coeff * (ν)⁴"
    σ_Rayl_coeff::FT  #Cross-section = σ_Rayl_coeff * ν⁴
    σ_Rayl_coeff_hires::AbstractArray{FT}
    Δν̃_Rayl_coeff_hires::AbstractArray{FT}
    
    #"Wavenumber shift in emitted light after vᵢ, Jᵢ --> v_f, J_f transition, νₛ = ν + Δνₛ"
    #Δνₛ::AbstractArray{FT}  
    "Raman scattering coefficient for pure vibrational scattering"
    σ_VibRaman_coeff_0to1::FT
    σ_VibRaman_coeff_1to0::FT
    σ_VibRaman_coeff_0to1_hires::AbstractArray{FT} 
    σ_VibRaman_coeff_1to0_hires::AbstractArray{FT}
    
    "Vibrational Raman scattering wavenumber displacement Δν̃"
    Δν̃_VibRaman_coeff_0to1::FT 
    Δν̃_VibRaman_coeff_1to0::FT
    Δν̃_VibRaman_coeff_0to1_hires::AbstractArray{FT}
    Δν̃_VibRaman_coeff_1to0_hires::AbstractArray{FT}
    
    "Raman scattering coefficient for rovibrational scattering"
    σ_RoVibRaman_coeff_0to1_JtoJm2::AbstractArray{FT} 
    σ_RoVibRaman_coeff_0to1_JtoJp2::AbstractArray{FT} 
    σ_RoVibRaman_coeff_1to0_JtoJm2::AbstractArray{FT}      
    σ_RoVibRaman_coeff_1to0_JtoJp2::AbstractArray{FT}  #Cross-section = σ_VibRaman_coeff * ν⁴
    "Rovibrational wavenumber displacement Δν̃"
    Δν̃_RoVibRaman_coeff_0to1_JtoJm2::AbstractArray{FT} 
    Δν̃_RoVibRaman_coeff_0to1_JtoJp2::AbstractArray{FT} 
    Δν̃_RoVibRaman_coeff_1to0_JtoJm2::AbstractArray{FT} 
    Δν̃_RoVibRaman_coeff_1to0_JtoJp2::AbstractArray{FT} 
    "Rotational Raman cross-section σ(ν, vᵢ, Jᵢ --> v_f, J_f) = σ_Raman_coeff * (ν+Δνₛ)⁴"
    σ_RoRaman_coeff_JtoJm2::AbstractArray{FT}
    σ_RoRaman_coeff_JtoJp2::AbstractArray{FT} #Cross-section = σ_Raman_coeff * ν⁴
    "Rotational Raman wavenumber displacement Δν̃"
    Δν̃_RoRaman_coeff_JtoJm2::AbstractArray{FT} 
    Δν̃_RoRaman_coeff_JtoJp2::AbstractArray{FT} #Cross-section = σ_Raman_coeff * ν⁴
    "Energy levels (in wavenumbers [cm^{-1}]) for a given v, J state"
    E_vJ::AbstractArray{FT}  
end


Base.@kwdef struct MolecularConstants{FT}
    "Volume mixing ratio"
    vmr::FT
    "Dunham Expansion Coefficients (5x5 Matrix)"
    Y::AbstractArray{FT}  
    "Nuclear Spin Multiplicity (Array of length 2 {odd J_initial, even J_initial})"
    gₛ::AbstractArray{Int} 
    PolTensor::PolarizationTensor{FT}
    effCoeff::EffectiveCoefficients{FT}
end







