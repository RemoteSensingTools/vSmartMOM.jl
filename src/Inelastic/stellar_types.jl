"""
    type AbstractRamanType
Abstract Raman type 
"""
abstract type AbstractRamanType  end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_RRS{FT<:AbstractFloat} <: AbstractRamanType 

    "Molecular Constants for H2"
    h2::InelasticScattering.MolecularConstants{Float64}
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
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct sol_VS_0to1{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constant for H2"
    h2::InelasticScattering.MolecularConstants{FT}
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    #ramanAtmoProp::RamanAtmosphereProperties
    fscattRayl#::FT
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    dτ₀::FT
    dτ₀_λ::FT
    k_Rayl_scatt::FT #σ_Rayl(λ_scatt)/σ_Rayl(λ_incident)
    n_Raman::Int
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct sol_VS_1to0{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constant for H2"
    h2::InelasticScattering.MolecularConstants{FT}
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    fscattRayl#::FT
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    dτ₀::FT
    dτ₀_λ::FT
    k_Rayl_scatt::FT #σ_Rayl(λ_scatt)/σ_Rayl(λ_incident)
    n_Raman::Int
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #ramanAtmoProp::RamanAtmosphereProperties
end
#=
"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RVRS{FT<:AbstractFloat} <: AbstractRamanType 
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
=#



############################################################
############### Types for Concatenated Mode ################
############################################################


"""
    struct RRS_plus{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters for the concatenated mode of vSmartMOM
# Fields
$(DocStringExtensions.FIELDS)
"""
#=Base.@kwdef mutable struct RRS_plus{FT<:AbstractFloat} <: AbstractRamanType 

    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = [1]
    grid_in::Array{StepRangeLen{FT},1} 

    "Molecular Constants for H2"
    h2::InelasticScattering.MolecularConstants{Float64}
    
    "Greek coeffs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    # ramanAtmoProp::RamanAtmosphereProperties
    #values for each band
    fscattRayl#::Array{FT,1}
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    
    ϖ_λ₁λ₀::Array{FT,2} #last index represents the band iB
    i_λ₁λ₀::Array{Int,2} #last index represents the band iB

    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    i_ref::Int
    n_Raman::Int
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end
=#
"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType 
    "Concatenated indices of band limits"
    bandSpecLim::Vector{Any} = []#Array{UnitRange{Int64},1}
    iBand::Vector{Any}       = []   #Array{Int,1} 
    grid_in::Vector{Any}     = []#Vector{StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}}#AbstractRange{<:Real}#Vector{StepRangeLen{FT}} 

    "Molecular Constant for H2"
    h2::InelasticScattering.MolecularConstants{FT}
        
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs       = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    greek_raman_VS::GreekCoefs = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    
    "Pre-computed optical properties"
    #ramanAtmoProp::RamanAtmosphereProperties
    #values for each band
    fscattRayl::Array{FT,1} = zeros(FT,1)
    ϖ_Cabannes::Array{FT,1} = zeros(FT,1)#elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    
    ϖ_λ₁λ₀::Array{FT,1}   = zeros(FT,1)#last index represents the band iB
    i_λ₁λ₀::Array{Int,1}  = zeros(Int,1)#last index represents the band iB
    Z⁻⁺_λ₁λ₀::Array{FT,2} = zeros(FT,1,1)#independent of band
    Z⁺⁺_λ₁λ₀::Array{FT,2} = zeros(FT,1,1) #independent of band
    
    ϖ_λ₁λ₀_VS::Array{FT,1}   = zeros(FT,1)#last index represents the band iB
    i_λ₁λ₀_VS::Array{Int,1}  = zeros(Int,1)#last index represents the band iB
    Z⁻⁺_λ₁λ₀_VS::Array{FT,2} = zeros(FT,1,1)#independent of band
    Z⁺⁺_λ₁λ₀_VS::Array{FT,2} = zeros(FT,1,1)#independent of band

    i_λ₁λ₀_all::Array{Int,1}    = zeros(Int,1)
    i_ref::Int                  = 1
    n_Raman::Int                = 1
    F₀::Array{FT,2}             = zeros(FT,1,1) # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_VS_1to0_plus{FT<:AbstractFloat} <: AbstractRamanType 
    
    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = []
    grid_in::Array{StepRangeLen{FT},1} 

    "Molecular Constant for H2"
    h2::InelasticScattering.MolecularConstants{FT}
        
    "Greek coefs in Raman calculations" 
    greek_raman::GreekCoefs
    greek_raman_VS::GreekCoefs
    
    "Pre-computed optical properties"
    #values for each band
    fscattRayl::Array{FT,1} 
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering

    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}

    ϖ_λ₁λ₀_VS::Array{FT,1} 
    i_λ₁λ₀_VS::Array{Int,1} 
    Z⁻⁺_λ₁λ₀_VS::Array{FT,2} #independent of band
    Z⁺⁺_λ₁λ₀_VS::Array{FT,2} #independent of band

    i_λ₁λ₀_all::Array{Int,1}
    i_ref::Int
    #dτ₀::FT
    #dτ₀_λ::FT
    #k_Rayl_scatt::FT #σ_Rayl(λ_scatt)/σ_Rayl(λ_incident)
    n_Raman::Int
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #ramanAtmoProp::RamanAtmosphereProperties
end


# Replace all struct names <name> with sol_<name>
#=
Base.@kwdef mutable struct sol_RRS{FT<:AbstractFloat} <: AbstractRamanType 
    h2::InelasticScattering.MolecularConstants{Float64}
    
    greek_raman::GreekCoefs
    fscattRayl#::Array{FT,1}
    ϖ_Cabannes::Array{FT,1}
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    i_ref::Int
    n_Raman::Int
    bandSpecLim = []
    iBand = 1
    F₀::Array{FT,2}
    #SIF₀::Array{FT,2}
end

Base.@kwdef struct sol_VS_0to1{FT<:AbstractFloat} <: AbstractRamanType 
    h2::InelasticScattering.MolecularConstants{FT}
    greek_raman::GreekCoefs
    fscattRayl#::FT
    ϖ_Cabannes::FT
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    dτ₀::FT
    dτ₀_λ::FT
    k_Rayl_scatt::FT
    n_Raman::Int
    F₀::Array{FT,2}
    #SIF₀::Array{FT,2}
end

Base.@kwdef struct sol_VS_1to0{FT<:AbstractFloat} <: AbstractRamanType 
    h2::InelasticScattering.MolecularConstants{FT}
    greek_raman::GreekCoefs
    fscattRayl#::FT
    ϖ_Cabannes::FT
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    dτ₀::FT
    dτ₀_λ::FT
    k_Rayl_scatt::FT
    n_Raman::Int
    F₀::Array{FT,2}
    #SIF₀::Array{FT,2}
end

#=
Base.@kwdef mutable struct sol_RRS_plus{FT<:AbstractFloat} <: AbstractRamanType 
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = [1]
    grid_in::Array{StepRangeLen{FT},1}
    h2::InelasticScattering.MolecularConstants{Float64}
    greek_raman::GreekCoefs
    fscattRayl#::Array{FT,1}
    ϖ_Cabannes::Array{FT,1}
    ϖ_λ₁λ₀::Array{FT,2}
    i_λ₁λ₀::Array{Int,2}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    i_ref::Int
    n_Raman::Int
    F₀::Array{FT,2}
    #SIF₀::Array{FT,2}
end
=#
Base.@kwdef mutable struct sol_VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType 
    bandSpecLim::Vector{Any} = []
    iBand::Vector{Any} = []
    grid_in::Vector{Any} = []
    h2::InelasticScattering.MolecularConstants{FT}
    greek_raman::GreekCoefs = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    greek_raman_VS::GreekCoefs = GreekCoefs([FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)], [FT(1)])
    fscattRayl::Array{FT,1} = zeros(FT,1)
    ϖ_Cabannes::Array{FT,1} = zeros(FT,1)
    ϖ_λ₁λ₀::Array{FT,1} = zeros(FT,1)
    i_λ₁λ₀::Array{Int,1} = zeros(Int,1)
    Z⁻⁺_λ₁λ₀::Array{FT,2} = zeros(FT,1,1)
    Z⁺⁺_λ₁λ₀::Array{FT,2} = zeros(FT,1,1)
    ϖ_λ₁λ₀_VS::Array{FT,1} = zeros(FT,1)
    i_λ₁λ₀_VS::Array{Int,1} = zeros(Int,1)
    Z⁻⁺_λ₁λ₀_VS::Array{FT,2} = zeros(FT,1,1)
    Z⁺⁺_λ₁λ₀_VS::Array{FT,2} = zeros(FT,1,1)
    i_λ₁λ₀_all::Array{Int,1} = zeros(Int,1)
    i_ref::Int = 1
    n_Raman::Int = 1
    F₀::Array{FT,2} = zeros(FT,1,1)
    #SIF₀::Array{FT,2}
end

Base.@kwdef mutable struct sol_VS_1to0_plus{FT<:AbstractFloat} <: AbstractRamanType 
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = []
    grid_in::Array{StepRangeLen{FT},1}
    h2::InelasticScattering.MolecularConstants{FT}
    greek_raman::GreekCoefs
    greek_raman_VS::GreekCoefs
    fscattRayl::Array{FT,1}
    ϖ_Cabannes::Array{FT,1}
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    ϖ_λ₁λ₀_VS::Array{FT,1}
    i_λ₁λ₀_VS::Array{Int,1}
    Z⁻⁺_λ₁λ₀_VS::Array{FT,2}
    Z⁺⁺_λ₁λ₀_VS::Array{FT,2}
    i_λ₁λ₀_all::Array{Int,1}
    i_ref::Int
    n_Raman::Int
    F₀::Array{FT,2}
    #SIF₀::Array{FT,2}
end
=#