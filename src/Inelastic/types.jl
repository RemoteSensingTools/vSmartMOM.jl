"""
    type AbstractRamanType
Abstract Raman type 
"""
abstract type AbstractRamanType{FT}  end

"""
    RRS{FT} <: AbstractRamanType{FT}

Rotational Raman Scattering (RRS) parameters for atmospheric N₂ and O₂.
RRS redistributes photons from Fraunhofer absorption lines into the
surrounding continuum ("filling-in" effect), which is important for
high-resolution spectral retrievals (e.g., O₂ A-band).

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
    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand = 1
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    VS_0to1{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering: Stokes transition (v=0 → v=1).
Models the energy transfer from incident photons to vibrational
excitation of N₂ and O₂ molecules, shifting scattered photons to
longer wavelengths (lower wavenumbers).

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
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    VS_1to0{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering: anti-Stokes transition (v=1 → v=0).
Models the energy transfer from vibrationally excited N₂ and O₂
molecules back to incident photons, shifting scattered photons to
shorter wavelengths (higher wavenumbers).

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
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #ramanAtmoProp::RamanAtmosphereProperties
end
#=
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

    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand::Vector{Int} = Int[1]

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
"""
    noRS{FT} <: AbstractRamanType{FT}

No Raman Scattering (elastic-only RT).  This is the default Raman type
used by `rt_run(model)`.  All scattering is treated as elastic
(Cabannes fraction = 1).
"""
Base.@kwdef mutable struct noRS{FT} <: AbstractRamanType{FT}
    fscattRayl::Array{FT,1} = [0.0]
    ϖ_Cabannes::Array{FT,1} = [1.0, 1.0, 1.0] #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand::Vector{Int} = Int[1]
    F₀ = zeros(Float64, 1, 1)  # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀ = zeros(Float64, 1, 1) # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

############################################################
############### Types for Concatenated Mode ################
############################################################


"""
    RRS_plus{FT} <: AbstractRamanType{FT}

Rotational Raman Scattering parameters for the concatenated (multi-band)
mode of vSmartMOM.  Extends [`RRS`](@ref) with band-concatenation bookkeeping.

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
    fscattRayl#::Array{FT,1}
    ϖ_Cabannes::Array{FT,1} #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering

    ϖ_λ₁λ₀::Array{FT,2} #last index represents the band iB
    i_λ₁λ₀::Array{Int,2} #last index represents the band iB

    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    i_ref::Int
    n_Raman::Int
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    VS_0to1_plus{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering (v=0→1, Stokes) for the concatenated
(multi-band) mode.  Extends [`VS_0to1`](@ref) with band-concatenation
bookkeeping and separate N₂/O₂ vibrational channels.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}
    "Concatenated indices of band limits"
    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand::Vector{Int}       = Int[]
    grid_in::Vector{AbstractRange{Float64}} = AbstractRange{Float64}[] 

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
    F₀::Array{FT,2}             = zeros(FT,1,1) # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    VS_1to0_plus{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering (v=1→0, anti-Stokes) for the concatenated
(multi-band) mode.  Extends [`VS_1to0`](@ref) with band-concatenation
bookkeeping and separate N₂/O₂ vibrational channels.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct VS_1to0_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}
    
    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Vector{Int} = Int[]
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
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #ramanAtmoProp::RamanAtmosphereProperties
end

"""
    noRS_plus{FT} <: AbstractRamanType{FT}

No Raman Scattering (elastic-only) for the concatenated (multi-band) mode.
Extends [`noRS`](@ref) with band-concatenation bookkeeping.
"""
Base.@kwdef mutable struct noRS_plus{FT} <: AbstractRamanType{FT}
    fscattRayl::Array{FT,1} = [0.0]
    ϖ_Cabannes::FT = 1.0 #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand::Vector{Int} = Int[]
    F₀::Array{FT,2} = zeros(FT, 1, 1) # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end


#==============================================================#
# ############### Types for Stellar Mode #######################
#==============================================================#

"""
    sol_RRS{FT} <: AbstractRamanType{FT}

Rotational Raman Scattering for stellar (non-solar) atmospheres,
using H₂ molecular constants instead of N₂/O₂.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_RRS{FT<:AbstractFloat} <: AbstractRamanType{FT}

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
    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand = 1
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    #SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    sol_VS_0to1{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering (v=0→1, Stokes) for stellar atmospheres,
using H₂ molecular constants.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct sol_VS_0to1{FT<:AbstractFloat} <: AbstractRamanType{FT}
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
    sol_VS_1to0{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering (v=1→0, anti-Stokes) for stellar atmospheres,
using H₂ molecular constants.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct sol_VS_1to0{FT<:AbstractFloat} <: AbstractRamanType{FT}
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
Base.@kwdef mutable struct RVRS{FT<:AbstractFloat} <: AbstractRamanType{FT} 
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}

    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand::Vector{Int} = Int[1]

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
#=Base.@kwdef mutable struct RRS_plus{FT<:AbstractFloat} <: AbstractRamanType{FT} 

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
    sol_VS_0to1_plus{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering (v=0→1, Stokes) for stellar atmospheres
in the concatenated (multi-band) mode.  Uses H₂ molecular constants.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}
    "Concatenated indices of band limits"
    bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]
    iBand::Vector{Int}       = Int[]
    grid_in::Vector{AbstractRange{Float64}} = AbstractRange{Float64}[] 

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
    sol_VS_1to0_plus{FT} <: AbstractRamanType{FT}

Vibrational Raman Scattering (v=1→0, anti-Stokes) for stellar atmospheres
in the concatenated (multi-band) mode.  Uses H₂ molecular constants.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_VS_1to0_plus{FT<:AbstractFloat} <: AbstractRamanType{FT}
    
    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Vector{Int} = Int[]
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