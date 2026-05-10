"""
    type AbstractRamanType
Abstract Raman type 
"""
abstract type AbstractRamanType  end

# Molecular constants stay in Float64 even when FT is Float32. The Raman
# cross-section prefactors can sit near Float32's subnormal range; widening
# these constants keeps the molecular state stable while RT work arrays remain
# parameterized by FT.

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct RRS{FT<:AbstractFloat} <: AbstractRamanType 

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
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VS_0to1{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{Float64}
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
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VS_1to0{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{Float64}
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
Base.@kwdef mutable struct RVRS{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{Float64}

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
"""
    noRS(; kwargs...)
    noRS{FT}(; kwargs...)

Pure elastic-scattering mode. This mode carries no inelastic redistribution
state and is the default Raman type used by `rt_run(model)`. The untyped
keyword constructor infers `FT` from supplied floating-point arrays and falls
back to `Float64` when no typed arrays are provided.
"""
mutable struct noRS{FT} <: AbstractRamanType
    fscattRayl::Array{FT,1}
    # ϖ_Cabannes is indexed by iBand; default sized for up to 3 bands so
    # zero-arg `noRS()` works with multi-band smoke configs (e.g. EMIT 2-band).
    ϖ_Cabannes::Array{FT,1} # elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    bandSpecLim
    iBand::Array{Int,1}
    # F₀/SIF₀ placeholders: rt_run resizes to (pol_type.n, nSpec) before first use.
    # Defaults keep zero-arg `noRS()` ergonomic for pure-elastic callers.
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

_maybe_float_type(::Nothing) = nothing
_maybe_float_type(x::AbstractFloat) = typeof(x)
function _maybe_float_type(x)
    T = try
        eltype(x)
    catch
        return nothing
    end
    return T <: AbstractFloat ? T : nothing
end

function _configured_float_type(default::Type{<:AbstractFloat}, values...)
    FT = default
    has_configured_type = false
    for value in values
        T = _maybe_float_type(value)
        T === nothing && continue
        FT = has_configured_type ? promote_type(FT, T) : T
        has_configured_type = true
    end
    return FT
end

_config_vector(::Type{FT}, value, default) where {FT} =
    value === nothing ? default : collect(FT, value)
_config_matrix(::Type{FT}, value, default) where {FT} =
    value === nothing ? default : Array{FT,2}(value)
_config_scalar(::Type{FT}, value, default) where {FT} =
    value === nothing ? default : FT(value)

function noRS{FT}(; fscattRayl = FT[0],
                    ϖ_Cabannes = FT[1, 1, 1],
                    bandSpecLim = [],
                    iBand = [1],
                    F₀ = zeros(FT, 1, 1),
                    SIF₀ = zeros(FT, 1, 1)) where {FT}
    return noRS{FT}(
        _config_vector(FT, fscattRayl, FT[0]),
        _config_vector(FT, ϖ_Cabannes, FT[1, 1, 1]),
        bandSpecLim,
        collect(Int, iBand),
        _config_matrix(FT, F₀, zeros(FT, 1, 1)),
        _config_matrix(FT, SIF₀, zeros(FT, 1, 1)),
    )
end

function noRS(; fscattRayl = nothing,
                ϖ_Cabannes = nothing,
                bandSpecLim = [],
                iBand = [1],
                F₀ = nothing,
                SIF₀ = nothing)
    FT = _configured_float_type(Float64, fscattRayl, ϖ_Cabannes, F₀, SIF₀)
    return noRS{FT}(
        _config_vector(FT, fscattRayl, FT[0]),
        _config_vector(FT, ϖ_Cabannes, FT[1, 1, 1]),
        bandSpecLim,
        collect(Int, iBand),
        _config_matrix(FT, F₀, zeros(FT, 1, 1)),
        _config_matrix(FT, SIF₀, zeros(FT, 1, 1)),
    )
end

############################################################
############### Types for Concatenated Mode ################
############################################################


"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType 
    "Concatenated indices of band limits"
    bandSpecLim::Vector{UnitRange{Int64}} = UnitRange{Int64}[]
    iBand::Vector{Int}       = Int[]
    grid_in::Vector{StepRangeLen{FT}} = StepRangeLen{FT}[]

    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{Float64}

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
    SIF₀::Array{FT,2}           = zeros(FT,1,1)# Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

"""
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct VS_1to0_plus{FT<:AbstractFloat} <: AbstractRamanType 
    
    "Concatenated indices of band limits"
    bandSpecLim = Array{UnitRange{Int64},1}
    iBand::Array{Int,1} = []
    grid_in::Array{StepRangeLen{FT},1} 

    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{Float64}
    
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
    noRS_plus(; kwargs...)
    noRS_plus{FT}(; kwargs...)

Concatenated-grid pure elastic-scattering mode used by plus-mode code paths
that share the Raman interface but do not include inelastic redistribution.
The untyped keyword constructor infers `FT` from supplied floating-point
inputs and falls back to `Float64` when no typed inputs are provided.
"""
mutable struct noRS_plus{FT} <: AbstractRamanType
    fscattRayl::Array{FT,1}
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    bandSpecLim
    iBand::Array{Int,1}
    # F₀/SIF₀ placeholders: rt_run resizes to (pol_type.n, nSpec) before first use.
    F₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
    SIF₀::Array{FT,2} # Solar/Stellar irradiation Stokes vector of size (pol_type.n, nSpec)
end

function noRS_plus{FT}(; fscattRayl = FT[0],
                         ϖ_Cabannes = one(FT),
                         bandSpecLim = [],
                         iBand = Int[],
                         F₀ = zeros(FT, 1, 1),
                         SIF₀ = zeros(FT, 1, 1)) where {FT}
    return noRS_plus{FT}(
        _config_vector(FT, fscattRayl, FT[0]),
        _config_scalar(FT, ϖ_Cabannes, one(FT)),
        bandSpecLim,
        collect(Int, iBand),
        _config_matrix(FT, F₀, zeros(FT, 1, 1)),
        _config_matrix(FT, SIF₀, zeros(FT, 1, 1)),
    )
end

function noRS_plus(; fscattRayl = nothing,
                     ϖ_Cabannes = nothing,
                     bandSpecLim = [],
                     iBand = Int[],
                     F₀ = nothing,
                     SIF₀ = nothing)
    FT = _configured_float_type(Float64, fscattRayl, ϖ_Cabannes, F₀, SIF₀)
    return noRS_plus{FT}(
        _config_vector(FT, fscattRayl, FT[0]),
        _config_scalar(FT, ϖ_Cabannes, one(FT)),
        bandSpecLim,
        collect(Int, iBand),
        _config_matrix(FT, F₀, zeros(FT, 1, 1)),
        _config_matrix(FT, SIF₀, zeros(FT, 1, 1)),
    )
end


#==============================================================#
# ############### Types for Stellar Mode #######################
#==============================================================#

"""
    sol_RRS{FT}

Stellar/solar rotational Raman-scattering mode for H₂-dominated atmospheres.
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
    h2::InelasticScattering.MolecularConstants{Float64}
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
    h2::InelasticScattering.MolecularConstants{Float64}
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
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{Float64}

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
    struct RRS{FT<:AbstractFloat}
A struct which defines Rotational Raman Scattering parameters
# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef mutable struct sol_VS_0to1_plus{FT<:AbstractFloat} <: AbstractRamanType 
    "Concatenated indices of band limits"
    bandSpecLim::Vector{UnitRange{Int64}} = UnitRange{Int64}[]
    iBand::Vector{Int}       = Int[]
    grid_in::Vector{StepRangeLen{FT}} = StepRangeLen{FT}[]

    "Molecular Constant for H2"
    h2::InelasticScattering.MolecularConstants{Float64}
        
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
    h2::InelasticScattering.MolecularConstants{Float64}
        
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

"""
    has_inelastic(rs::AbstractRamanType) -> Bool

Return whether `rs` represents a Raman-aware mode that carries inelastic
source/scattering state. This is the public trait used by CoreRT boundary code
instead of spelling out `isa noRS` checks at each call site.
"""
has_inelastic(::AbstractRamanType) = true
has_inelastic(::Union{noRS, noRS_plus}) = false

"""
    uses_cabannes_phase(rs::AbstractRamanType) -> Bool

Return whether the elastic Rayleigh path should use Cabannes phase coefficients
because Raman redistribution is handled explicitly by `rs`.
"""
uses_cabannes_phase(rs::AbstractRamanType) = has_inelastic(rs)

"""
    needs_interaction_workspace(rs::AbstractRamanType) -> Bool

Return whether CoreRT should allocate the staged inelastic interaction
workspace for this Raman mode.
"""
needs_interaction_workspace(rs::AbstractRamanType) = has_inelastic(rs)

"""
    needs_rayleigh_expansion(rs::AbstractRamanType) -> Bool

Return whether per-layer Rayleigh scattering fractions must be expanded onto
the active Raman spectral grid before entering the core kernel.
"""
needs_rayleigh_expansion(rs::AbstractRamanType) = has_inelastic(rs)

"""
    normalize_raman_weights!(rs, model, iBand)

Normalize Raman redistribution weights for modes whose `ϖ_λ₁λ₀` represents a
rotational Raman fraction. Modes that do not require this normalization are
left unchanged.
"""
normalize_raman_weights!(::AbstractRamanType, model, iBand) = nothing

function normalize_raman_weights!(rs::RRS, model, iBand)
    iB = iBand[1]
    rs.ϖ_λ₁λ₀ .*= (1 - model.ϖ_Cabannes[iB]) / sum(rs.ϖ_λ₁λ₀)
    return nothing
end
