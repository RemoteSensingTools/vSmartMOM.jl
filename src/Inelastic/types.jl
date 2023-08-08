"""
    type AbstractRamanType
Abstract Raman type 
"""
abstract type AbstractRamanType{FT}  end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RRS{FT<:AbstractFloat} <: AbstractRamanType{FT} 

    "Molecular Constants for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constants for O2"
    o2::InelasticScattering.MolecularConstants{Float64}
    "Greek coeffs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    # ramanAtmoProp::RamanAtmosphereProperties
    fscattRayl#::Array{FT,1}
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    i_ref::Int
    n_Raman::Int
    bandSpecLim = []
    iBand = 1
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VS_0to1{FT<:AbstractFloat} <: AbstractRamanType{FT}  
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    #ramanAtmoProp::RamanAtmosphereProperties
    fscattRayl::FT
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    dτ₀::FT
    dτ₀_λ::FT
    k_Rayl_scatt::FT #σ_Rayl(λ_scatt)/σ_Rayl(λ_incident)
    n_Raman::Int
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VS_1to0{FT<:AbstractFloat} <: AbstractRamanType{FT}  
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    fscattRayl::FT
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    dτ₀::FT
    dτ₀_λ::FT
    k_Rayl_scatt::FT #σ_Rayl(λ_scatt)/σ_Rayl(λ_incident)
    n_Raman::Int
    #ramanAtmoProp::RamanAtmosphereProperties
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RVRS{FT<:AbstractFloat} <: AbstractRamanType{FT}  
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}

    bandSpecLim = []
    iBand::Array{Int,1} = [1]

    "Pre-computed optical properties"
    #ramanAtmoProp::RamanAtmosphereProperties
    fscattRayl::Array{FT,1}
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,2} #last dimension -> band index
    i_λ₁λ₀::Array{Int,2} #last dimension -> band index
    Z⁻⁺_λ₁λ₀::Array{FT,3} #last dimension -> band index
    Z⁺⁺_λ₁λ₀::Array{FT,3} #last dimension -> band index
    τ₀::FT
    n_Raman::Int
end

Base.@kwdef mutable struct noRS{FT} <: AbstractRamanType{FT} 
    fscattRayl::Array{FT,1} = [0.0]
    ϖ_Cabannes::Array{FT,1} = [1.0,1.0,1.0] #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    bandSpecLim = []
    iBand::Array{Int,1} = [1]
end

############################################################
############### Types for Concatenated Mode ################
############################################################


"""
    struct RRS_plus{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters for the concatenated mode of vSmartMOM
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RRS_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}  

    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = [1]
    grid_in::Array{StepRangeLen{FT},1} 

    "Molecular Constants for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constants for O2"
    o2::InelasticScattering.MolecularConstants{Float64}
    "Greek coeffs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    # ramanAtmoProp::RamanAtmosphereProperties
    #values for each band
    fscattRayl::Array{FT,1}
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    
    ϖ_λ₁λ₀::Array{FT,2} #last index represents the band iB
    i_λ₁λ₀::Array{Int,2} #last index represents the band iB

    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    i_ref::Int
    n_Raman::Int
    
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}  
    "Concatenated indices of band limits"
    bandSpecLim::Vector{Any} = []#Array{UnitRange{Int64},1}
    iBand::Vector{Any}       = []   #Array{Int,1} 
    grid_in::Vector{Any}     = []#Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}#AbstractRange{<:Real}#Vector{StepRangeLen{FT}} 

    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}
    
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs       = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    greek_raman_VS_n2::GreekCoefs = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    greek_raman_VS_o2::GreekCoefs = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])

    "Pre-computed optical properties"
    #ramanAtmoProp::RamanAtmosphereProperties
    #values for each band
    fscattRayl::Array{FT,1} = zeros(FT,1)
    ϖ_Cabannes::Array{FT,1} = zeros(FT,1)#elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    
    ϖ_λ₁λ₀::Array{FT,1}   = zeros(FT,1)#last index represents the band iB
    i_λ₁λ₀::Array{Int,1}  = zeros(Int,1)#last index represents the band iB
    Z⁻⁺_λ₁λ₀::Array{FT,2} = zeros(FT,1,1)#independent of band
    Z⁺⁺_λ₁λ₀::Array{FT,2} = zeros(FT,1,1) #independent of band
    
    ϖ_λ₁λ₀_VS_n2::Array{FT,1}   = zeros(FT,1)#last index represents the band iB
    i_λ₁λ₀_VS_n2::Array{Int,1}  = zeros(Int,1)#last index represents the band iB
    Z⁻⁺_λ₁λ₀_VS_n2::Array{FT,2} = zeros(FT,1,1)#independent of band
    Z⁺⁺_λ₁λ₀_VS_n2::Array{FT,2} = zeros(FT,1,1)#independent of band

    ϖ_λ₁λ₀_VS_o2::Array{FT,1}   = zeros(FT,1) #last index represents the band iB
    i_λ₁λ₀_VS_o2::Array{Int,1}  = zeros(Int,1) #last index represents the band iB
    Z⁻⁺_λ₁λ₀_VS_o2::Array{FT,2} = zeros(FT,1,1) #independent of band
    Z⁺⁺_λ₁λ₀_VS_o2::Array{FT,2} = zeros(FT,1,1) #independent of band

    i_λ₁λ₀_all::Array{Int,1}    = zeros(Int,1)
    i_ref::Int                  = 1
    n_Raman::Int                = 1
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct VS_1to0_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}  
    
    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = []
    grid_in::Array{StepRangeLen{FT},1} 

    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}
    
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs
    greek_raman_VS_n2::GreekCoefs
    greek_raman_VS_o2::GreekCoefs

    "Pre-computed optical properties"
    #values for each band
    fscattRayl::Array{FT,1} 
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering

    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}

    ϖ_λ₁λ₀_VS_n2::Array{FT,1} 
    i_λ₁λ₀_VS_n2::Array{Int,1} 
    Z⁻⁺_λ₁λ₀_VS_n2::Array{FT,2} #independent of band
    Z⁺⁺_λ₁λ₀_VS_n2::Array{FT,2} #independent of band

    ϖ_λ₁λ₀_VS_o2::Array{FT,1} 
    i_λ₁λ₀_VS_o2::Array{Int,1} 
    Z⁻⁺_λ₁λ₀_VS_o2::Array{FT,2} #independent of band
    Z⁺⁺_λ₁λ₀_VS_o2::Array{FT,2} #independent of band

    i_λ₁λ₀_all::Array{Int,1}
    i_ref::Int
    #dτ₀::FT
    #dτ₀_λ::FT
    #k_Rayl_scatt::FT #σ_Rayl(λ_scatt)/σ_Rayl(λ_incident)
    n_Raman::Int
    #ramanAtmoProp::RamanAtmosphereProperties
end

Base.@kwdef mutable struct noRS_plus{FT} <: AbstractRamanType{FT} 
    fscattRayl::Array{FT,1} = [0.0]
    ϖ_Cabannes::FT = 1.0 #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    bandSpecLim = []
    iBand::Array{Int,1} = []
end





