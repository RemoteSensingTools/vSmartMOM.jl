#=
 
This file contains the linearization of all relevant types that are used in the vSmartMOM module:

- `AtmosphericProfile` stores all relevant atmospheric profile information 
- `AbstractObsGeometry` specifies the RT geometry
- `RT_Aerosol` holds an Aerosol with additional RT parameters
- `AbstractQuadratureType` specifies the quadrature type to use
- `AbstractSourceType` specifies the source type
- `CompositeLayer` and `AddedLayer` specify the layer properties
- `AbstractScatteringInterface` specifies the scattering interface type
- `AbstractSurfaceType` specify the type of surface in the RT simulation
- `AbsorptionParameters`, `ScatteringParameters`, and `vSmartMOM_Model` hold model parameters
- `QuadPoints` holds quadrature points, weights, etc. 
- `ComputedAtmosphereProperties` and `ComputedLayerProperties` hold intermediate computed properties

=#

"Struct for an atmospheric profile"
struct linAtmosphericProfile{FT}#{FT, VMR <: Union{Real, Vector}}
    "Temperature Profile"
    T::Array{FT,1}
    "Pressure Profile (Full)"
    p_full::Array{FT,1}
    "Specific humidity profile"
    q::Array{FT,1}
    "Pressure Levels"
    p_half::Array{FT,1}
    "H2O Volume Mixing Ratio Profile"
    vmr_h2o::Array{FT,1}
    "Vertical Column Density (Dry)"
    vcd_dry::Array{FT,1}
    "Vertical Column Density (H2O)"
    vcd_h2o::Array{FT,1}
    "Volume Mixing Ratio of Constituent Gases"
    vmr_co2::Array{FT,1} #Dict{String, VMR}
    "Layer height (meters)"
    Δz::Array{FT,1}
    "Layer altitude (meters)"
    z::Array{FT,1}
    dzdpsurf::Array{FT,1}
    "Derivative of vmr_h2o with respect to multiplicative factors in the boundary layer and above"
    dVMR_H2O::Array{FT,2}
    "Derivative of vmr_co2 with respect to psurf, 3 multiplicative factors for uniform, exponential, and lognormal profiles, and their respective parameters H, z₀, and σ"
    dVMR_CO2::Array{FT,2}
end
#=
"Types for describing atmospheric parameters"
abstract type AbstractObsGeometry end

"Observation Geometry (basics)" 
Base.@kwdef struct ObsGeometry{FT} <: AbstractObsGeometry
    "Solar Zenith Angle `[Degree]`"
    sza::FT
    "Viewing Zenith Angle(s) `[Degree]`" 
    vza::Array{FT,1}
    "Viewing Azimuth Angle(s) `[Degree]`" 
    vaz::Array{FT,1}
    "Altitude of Observer `[Pa]`"
    obs_alt::FT
end
=#
#=
mutable struct RT_Aerosol{}#FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "Aerosol"
    aerosol::Aerosol#{FT}
    "Reference τ"
    τ_ref#::FT
    "Vertical distribution as function of p (using Distributions.jl)"
    profile::Distribution#::FT
end
=#
#=
"Quadrature Types for RT streams"
abstract type AbstractQuadratureType end

"""
    struct RadauQuad

Use Gauss Radau quadrature scheme, which includes the SZA as point (see Sanghavi vSmartMOM)

"""
struct RadauQuad <:AbstractQuadratureType end

"""
    struct GaussQuadHemisphere

Use Gauss quadrature scheme, define interval [-1,1] within an hemisphere (90⁰), repeat for both

"""
struct GaussQuadHemisphere <:AbstractQuadratureType end

"""
    struct GaussQuadFullSphere

Use Gauss quadrature scheme, define interval [-1,1] for full sphere (180⁰), take half of it (less points near horizon compared to GaussQuadHemisphere)

"""
struct GaussQuadFullSphere <:AbstractQuadratureType end

"Abstract Type for Source Function Integration"
abstract type AbstractSourceType end

"Use Dummy Node Integration (DNI), SZA has to be a full node, Radau required"
struct DNI <:AbstractSourceType end

"Use Source Function Integration (SFI), Solar beam embedded in source term, can work with all quadrature schemes (faster)"
struct SFI <:AbstractSourceType end

"Abstract Type for Layer R,T and J matrices"
abstract type AbstractLayer end
=#
"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct linCompositeLayer{FT} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    dR⁻⁺::AbstractArray{FT,4}
    "Composite layer Reflectance matrix R (from - -> +)"
    dR⁺⁻::AbstractArray{FT,4}
    "Composite layer transmission matrix T (from + -> +)"
    dT⁺⁺::AbstractArray{FT,4}
    "Composite layer transmission matrix T (from - -> -)"
    dT⁻⁻::AbstractArray{FT,4}
    "Composite layer source matrix J (in + direction)"
    dJ₀⁺::AbstractArray{FT,4}
    "Composite layer source matrix J (in - direction)"
    dJ₀⁻::AbstractArray{FT,4}
end

"Added (Single) Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct linAddedLayer{FT} <: AbstractLayer 
    "Added layer Reflectance matrix R (from + -> -)"
    dr⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix T (from + -> +)"
    dt⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from - -> +)"
    dr⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix T (from - -> -)"
    dt⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix J (in + direction)"
    dj₀⁺::AbstractArray{FT,4}
    "Added layer source matrix J (in - direction)"
    dj₀⁻::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from + -> -)"
    dxr⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix T (from + -> +)"
    dxt⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix R (from - -> +)"
    dxr⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix T (from - -> -)"
    dxt⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix J (in + direction)"
    dxj₀⁺::AbstractArray{FT,4}
    "Added layer source matrix J (in - direction)"
    dxj₀⁻::AbstractArray{FT,4}
end

"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
struct linCompositeLayerRS{FT} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    dR⁻⁺::AbstractArray{FT,3}
    "Composite layer Reflectance matrix R (from - -> +)"
    dR⁺⁻::AbstractArray{FT,3}
    "Composite layer transmission matrix T (from + -> +)"
    dT⁺⁺::AbstractArray{FT,3}
    "Composite layer transmission matrix T (from - -> -)"
    dT⁻⁻::AbstractArray{FT,3}
    "Composite layer source matrix J (in + direction)"
    dJ₀⁺::AbstractArray{FT,3}
    "Composite layer source matrix J (in - direction)"
    dJ₀⁻::AbstractArray{FT,3}

    # Additional Arrays for Raman scattering
    "Composite layer Reflectance matrix ieR (from + -> -)"
    dieR⁻⁺::AbstractArray{FT,4}
    "Composite layer Reflectance matrix ieR (from - -> +)"
    dieR⁺⁻::AbstractArray{FT,4}
    "Composite layer transmission matrix ieT (from + -> +)"
    dieT⁺⁺::AbstractArray{FT,4}
    "Composite layer transmission matrix ieT (from - -> -)"
    dieT⁻⁻::AbstractArray{FT,4}
    "Composite layer source matrix ieJ (in + direction)"
    dieJ₀⁺::AbstractArray{FT,4}
    "Composite layer source matrix ieJ (in - direction)"
    dieJ₀⁻::AbstractArray{FT,4}
end

"Added (Single) Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
struct linAddedLayerRS{FT} <: AbstractLayer 
    "Added layer Reflectance matrix R (from + -> -)"
    dr⁻⁺::AbstractArray{FT,3}
    "Added layer transmission matrix T (from + -> +)"
    dt⁺⁺::AbstractArray{FT,3}
    "Added layer Reflectance matrix R (from - -> +)"
    dr⁺⁻::AbstractArray{FT,3}
    "Added layer transmission matrix T (from - -> -)"
    dt⁻⁻::AbstractArray{FT,3}
    "Added layer source matrix J (in + direction)"
    dj₀⁺::AbstractArray{FT,3}
    "Added layer source matrix J (in - direction)"
    dj₀⁻::AbstractArray{FT,3}

    # Additional Arrays for Raman scattering
    "Added layer Reflectance matrix ieR (from + -> -)"
    dier⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix ieT (from + -> +)"
    diet⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix ieR (from - -> +)"
    dier⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix ieT (from - -> -)"
    diet⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix ieJ (in + direction)"
    diej₀⁺::AbstractArray{FT,4}
    "Added layer source matrix ieJ (in - direction)"
    diej₀⁻::AbstractArray{FT,4}

    "Added layer Reflectance matrix R (from + -> -)"
    dxr⁻⁺::AbstractArray{FT,3}
    "Added layer transmission matrix T (from + -> +)"
    dxt⁺⁺::AbstractArray{FT,3}
    "Added layer Reflectance matrix R (from - -> +)"
    dxr⁺⁻::AbstractArray{FT,3}
    "Added layer transmission matrix T (from - -> -)"
    dxt⁻⁻::AbstractArray{FT,3}
    "Added layer source matrix J (in + direction)"
    dxj₀⁺::AbstractArray{FT,3}
    "Added layer source matrix J (in - direction)"
    dxj₀⁻::AbstractArray{FT,3}

    # Additional Arrays for Raman scattering
    "Added layer Reflectance matrix ieR (from + -> -)"
    diexr⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix ieT (from + -> +)"
    diext⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix ieR (from - -> +)"
    diexr⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix ieT (from - -> -)"
    diext⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix ieJ (in + direction)"
    diexj₀⁺::AbstractArray{FT,4}
    "Added layer source matrix ieJ (in - direction)"
    diexj₀⁻::AbstractArray{FT,4}
end
# Multisensor Composite layers 
# Elastic
"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct linCompositeLayerMS{M} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    dtopR⁻⁺::M#AbstractArray{FT,4}
    "Composite layer Reflectance matrix R (from - -> +)"
    dtopR⁺⁻::M
    "Composite layer transmission matrix T (from + -> +)"
    dtopT⁺⁺::M
    "Composite layer transmission matrix T (from - -> -)"
    dtopT⁻⁻::M
    "Composite layer source matrix J (in + direction)"
    dtopJ₀⁺::M
    "Composite layer source matrix J (in - direction)"
    dtopJ₀⁻::M
    "Composite layer Reflectance matrix R (from + -> -)"
    dbotR⁻⁺::M
    "Composite layer Reflectance matrix R (from - -> +)"
    dbotR⁺⁻::M
    "Composite layer transmission matrix T (from + -> +)"
    dbotT⁺⁺::M
    "Composite layer transmission matrix T (from - -> -)"
    dbotT⁻⁻::M
    "Composite layer source matrix J (in + direction)"
    dbotJ₀⁺::M
    "Composite layer source matrix J (in - direction)"
    dbotJ₀⁻::M
end

# Inelastic
"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct linCompositeLayerMSRS{M1, M2} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    dtopR⁻⁺::M1
    "Composite layer Reflectance matrix R (from - -> +)"
    dtopR⁺⁻::M1
    "Composite layer transmission matrix T (from + -> +)"
    dtopT⁺⁺::M1
    "Composite layer transmission matrix T (from - -> -)"
    dtopT⁻⁻::M1
    "Composite layer source matrix J (in + direction)"
    dtopJ₀⁺::M1 
    "Composite layer source matrix J (in - direction)"
    dtopJ₀⁻::M1

    # Additional Arrays for Raman scattering
    "Composite layer Reflectance matrix ieR (from + -> -)"
    dtopieR⁻⁺::M2
    "Composite layer Reflectance matrix ieR (from - -> +)"
    dtopieR⁺⁻::M2
    "Composite layer transmission matrix ieT (from + -> +)"
    dtopieT⁺⁺::M2
    "Composite layer transmission matrix ieT (from - -> -)"
    dtopieT⁻⁻::M2
    "Composite layer source matrix ieJ (in + direction)"
    dtopieJ₀⁺::M2
    "Composite layer source matrix ieJ (in - direction)"
    dtopieJ₀⁻::M2

    "Composite layer Reflectance matrix R (from + -> -)"
    dbotR⁻⁺::M1
    "Composite layer Reflectance matrix R (from - -> +)"
    dbotR⁺⁻::M1
    "Composite layer transmission matrix T (from + -> +)"
    dbotT⁺⁺::M1
    "Composite layer transmission matrix T (from - -> -)"
    dbotT⁻⁻::M1
    "Composite layer source matrix J (in + direction)"
    dbotJ₀⁺::M1
    "Composite layer source matrix J (in - direction)"
    dbotJ₀⁻::M1

    # Additional Arrays for Raman scattering
    "Composite layer Reflectance matrix ieR (from + -> -)"
    dbotieR⁻⁺::M2
    "Composite layer Reflectance matrix ieR (from - -> +)"
    dbotieR⁺⁻::M2
    "Composite layer transmission matrix ieT (from + -> +)"
    dbotieT⁺⁺::M2
    "Composite layer transmission matrix ieT (from - -> -)"
    dbotieT⁻⁻::M2
    "Composite layer source matrix ieJ (in + direction)"
    dbotieJ₀⁺::M2
    "Composite layer source matrix ieJ (in - direction)"
    dbotieJ₀⁻::M2
end
#=
"Abstract Type for Scattering Interfaces" 
abstract type AbstractScatteringInterface end
"No scattering in either the added layer or the composite layer" 
struct ScatteringInterface_00 <: AbstractScatteringInterface end
"No scattering in inhomogeneous composite layer; Scattering in homogeneous layer, added to bottom of the composite layer." 
struct ScatteringInterface_01 <: AbstractScatteringInterface end
"Scattering in inhomogeneous composite layer; No scattering in homogeneous layer which is added to bottom of the composite layer." 
struct ScatteringInterface_10 <: AbstractScatteringInterface end
"Scattering in inhomogeneous composite layer; Scattering in homogeneous layer, added to bottom of the composite layer." 
struct ScatteringInterface_11 <: AbstractScatteringInterface end
"Scattering Interface between Surface and Composite Layer" 
struct ScatteringInterface_AtmoSurf <: AbstractScatteringInterface end

"Abstract Type for Surface Types" 
abstract type AbstractSurfaceType end

"Lambertian Surface (scalar per band)"
struct LambertianSurfaceScalar{FT} <: AbstractSurfaceType
    "Albedo (scalar)"
    albedo::FT
end

"Defined as Array (has to have the same length as the band!)"
struct LambertianSurfaceSpectrum{FT} <: AbstractSurfaceType
    "Albedo (vector)"
    albedo::AbstractArray{FT,1}
end

"Defined as Array (has to have the same length as the band!)"
struct rpvSurfaceScalar{FT} <: AbstractSurfaceType
    "Overall reflectance level parameter (scalar)"
    ρ₀::FT
    "Hotspot function parameter (1.0 = no hotspot)"
    ρ_c::FT
    "Anisotropy shape parameter. k < 1.0 (> 1.0) corresponds to a bowl (bell) shape."
    k::FT
    "Asymmetry parameter, Θ < 0.0 (> 0.0) corresponds to a predominantly backward (forward) scattering."
    Θ::FT
end

struct RossLiSurfaceScalar{FT} <: AbstractSurfaceType
    "Volumetric RossThick  fraction"
    fvol::FT
    "Geometric LiSparse fraction"
    fgeo::FT
    "Isotropic reflectance fraction"
    fiso::FT
end

"Defined by Legendre polynomial terms as function of spectral grid, which is scaled to [-1,1] (degree derived from length of `a_coeff`)"
struct LambertianSurfaceLegendre{FT} <: AbstractSurfaceType
    "albedo = legendre_coeff[1] * P₀ + legendre_coeff[2]*P₁ + legendre_coeff[3]*P₂ + ... "
    legendre_coeff::AbstractArray{FT,1}
end
=#
#=
"""
    struct AbsorptionParameters

A struct which holds all absorption-related parameters (before any computations)
"""
mutable struct AbsorptionParameters
    "Molecules to use for absorption calculations (`nBand, nMolecules`)"
    molecules::AbstractArray
    "Volume-Mixing Ratios"
    vmr::Dict
    "Type of broadening function (Doppler/Lorentz/Voigt)"
    broadening_function::AbstractBroadeningFunction
    "Complex Error Function to use in Voigt calculations"
    CEF::AbstractComplexErrorFunction
    "Wing cutoff to use in cross-section calculation (cm⁻¹)"
    wing_cutoff::Integer
    "Lookup table type"
    luts::AbstractArray 
end
=#
#=
"""
    struct ScatteringParameters

A struct which holds all scattering-related parameters (before any computations)
"""
mutable struct ScatteringParameters{FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "List of scattering aerosols and their properties"
    rt_aerosols::Vector{RT_Aerosol}
    "Maximum aerosol particle radius for quadrature points/weights (µm)"
    r_max::FT
    "Number of quadrature points for integration of size distribution"
    nquad_radius::Integer
    "Reference wavelength (µm)"
    λ_ref::FT
    "Reference refractive index"
    n_ref::Complex{FT}
    "Algorithm to use for fourier decomposition (NAI2/PCW)"
    decomp_type::AbstractFourierDecompositionType
end
=#
#=
"""
    struct vSmartMOM_Parameters

A struct which holds all initial model parameters (before any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Parameters{FT<:Union{AbstractFloat, ForwardDiff.Dual}} 

    # radiative_transfer group
    "Spectral bands (`nBand`)"
    spec_bands::AbstractArray
    "Surface (Bidirectional Reflectance Distribution Function)"
    brdf::AbstractArray
    "Quadrature type for RT streams (RadauQuad/GaussQuadHemisphere/GaussQuadFullSphere)"
    quadrature_type::AbstractQuadratureType
    "Type of polarization (I/IQ/IQU/IQUV)"
    polarization_type::AbstractPolarizationType
    "Hard cutoff for maximum number of Fourier moments to loop over"
    max_m::Integer
    "Exclusion angle for forward peak [deg]"
    Δ_angle::FT
    "Truncation length for legendre terms (scalar for now, can do `nBand` later)"
    l_trunc::Integer
    "Depolarization factor"
    depol::FT    
    "Float type to use in the RT (Float64/Float32)"
    float_type::DataType
    "Architecture to use for calculations (CPU/GPU)"
    architecture::AbstractArchitecture

    # geometry group
    "Solar zenith angle [deg]"
    sza::FT
    "Viewing zenith angles [deg]"
    vza::AbstractArray{FT}
    "Viewing azimuthal angles [deg]"
    vaz::AbstractArray{FT}
    "Altitude of observer [Pa]"
    obs_alt::FT

    # atmospheric_profile group
    "Temperature Profile [K]"
    T::AbstractArray{FT}
    "Pressure Profile [hPa]"
    p::AbstractArray{FT}
    "Specific humidity profile"
    q::AbstractArray{FT}
    "Length of profile reduction"
    profile_reduction_n::Integer

    # absorption group
    "Optional struct that holds all absorption-related parameters"
    absorption_params::Union{AbsorptionParameters, Nothing}

    # scattering group
    "Optional struct that holds all aerosol scattering-related parameters"
    scattering_params::Union{ScatteringParameters, Nothing}
    
end
=#
#=
"""
    struct QuadPoints{FT} 

A struct which holds Quadrature points, weights, etc

# Fields
$(DocStringExtensions.FIELDS)
"""
struct QuadPoints{FT} 
    "μ₀, cos(SZA)"
    μ₀::FT
    "Index in quadrature points with sun"
    iμ₀::Int
    "Index in quadrature points with sun (in qp_μN)"
    iμ₀Nstart::Int
    "Quadrature points"
    qp_μ::AbstractArray{FT,1}
    "Weights of quadrature points"
    wt_μ::AbstractArray{FT,1}
    "Quadrature points (repeated for polarizations)"
    qp_μN::AbstractArray{FT,1}
    "Weights of quadrature points (repeated for polarizations)"
    wt_μN::AbstractArray{FT,1}
    "Number of quadrature points"
    Nquad::Int
end
=#
"""
    struct vSmartMOM_lin

A struct which holds all derived model parameters (including any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_lin #<:Union{Vector{Vector{dAerosolOptics}}, Vector{Array{FT},1}, Vector{Array{FT},3}} where FT#{FT} <: AbstractFloat#{linAE, linTAE, linTAE}#, linTAB, linTR, linTAE, linPRO}

    #"Struct with all individual parameters"
    #params::PA # vSmartMOM_Parameters
    
    "Truncated aerosol optics"
    #lin_aerosol_optics::dAerosolOptics{FT}  #::linAE # AbstractArray{AbstractArray{AerosolOptics}}
    lin_aerosol_optics::Any #Vector{Vector{dAerosolOptics}}
    #"Greek coefs in Rayleigh calculations" 
    #greek_rayleigh::GR # GreekCoefs
    #"Quadrature points/weights, etc"
    #quad_points::QP # QuadPoints

    #"Array to hold cross-sections over entire atmospheric profile"
    
    #"Rayleigh optical thickness"
    lin_τ_rayl::Any #Vector{Array{FT},1} #Array{Matrix{FT},1}#linTR # AbstractArray{AbstractArray}
    #"Molecular absorption"
    lin_τ_abs::Any #Vector{Array{FT},3} #linTAB # AbstractArray{AbstractArray} 
    #"Aerosol optical thickness"
    lin_τ_aer_psurf::Any #Vector{Array{FT},1} #Array{Matrix{FT},1}
    lin_τ_aer_z₀::Any #Vector{Array{FT},1}#Array{Matrix{FT},1}#::linTAE # AbstractArray{AbstractArray}
    lin_τ_aer_σz::Any #Vector{Array{FT},1}#Array{Matrix{FT},1}#::linTAE
    #"Observational Geometry (includes sza, vza, vaz)"
    #obs_geom::Ogeom # ObsGeometry
    #"Atmospheric profile to use"
    #lin_profile::linPRO #AtmosphericProfile
end

"""
    struct ComputedAtmosphereProperties

A struct which holds (for the entire atmosphere) all key layer optical properties required for the RT core solver

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct linComputedAtmosphereProperties

    "Absorption optical depth vectors (wavelength dependent)"
    lin_τ_λ_all
    "Albedo vectors (wavelength dependent)"
    lin_ϖ_λ_all
    "Absorption optical depth scalars (not wavelength dependent)"
    lin_τ_all
    "Albedo scalars (not wavelength dependent)"
    lin_ϖ_all
    "Combined Z moments (forward)"
    lin_Z⁺⁺_all
    "Combined Z moments (backward)"
    lin_Z⁻⁺_all
    "Maximum dτs"
    lin_dτ_max_all
    "dτs"
    lin_dτ_all
    #"Number of doublings (for all layers)"
    #ndoubl_all
    "dτs (wavelength dependent)"
    lin_dτ_λ_all
    "All expk"
    lin_expk_all
    #"Scattering flags"
    #scatter_all
    "Sum of optical thicknesses of all layers above the current layer"
    lin_τ_sum_all
    #"elastic (Cabannes) scattering fraction of Rayleigh (Cabannes+Raman) scattering per layer"
    #ϖ_Cabannes_all
    "Rayleigh fraction of scattering cross section per layer"
    lin_fscattRayl_all
    #"Scattering interface type for each layer"
    #scattering_interfaces_all
end



"""
    struct ComputedLayerProperties

A struct which holds all key layer optical properties required for the RT core solver

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct linComputedLayerProperties

    "Absorption optical depth vector (wavelength dependent)"
    lin_τ_λ 
    "Albedo vector (wavelength dependent)"
    lin_ϖ_λ 
    "Absorption optical depth scalar (not wavelength dependent)"
    lin_τ 
    "Albedo scalar (not wavelength dependent)"
    lin_ϖ  
    "Combined Z moment (forward)"
    lin_Z⁺⁺ 
    "Combined Z moment (backward)"
    lin_Z⁻⁺ 
    "Maximum dτ"
    lin_dτ_max 
    "dτ"
    lin_dτ     
    #"Number of doublings"
    #ndoubl
    "dτ (wavelength dependent)"
    lin_dτ_λ 
    "expk"
    lin_expk 
    #"Scattering flag"
    #scatter 
    "Sum of optical thicknesses of all layers above the current layer"
    lin_τ_sum
    "Fraction of scattering caused by Rayleigh"
    lin_fscattRayl
    #"Elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering"
    #ϖ_Cabannes 
    #"Scattering interface type for current layer"
    #scattering_interface
end

abstract type linAbstractOpticalProperties end

# Core optical Properties COP
Base.@kwdef struct linCoreScatteringOpticalProperties{FT,FT2,FT3} <:  linAbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    lin_τ::FT 
    "Single scattering albedo"
    lin_ϖ::FT2   
    "Z scattering matrix (forward)"
    lin_Z⁺⁺::FT3 
    "Z scattering matrix (backward)"
    lin_Z⁻⁺::FT3
end

# Core optical Properties COP with directional cross section 
Base.@kwdef struct linCoreDirectionalScatteringOpticalProperties{FT,FT2,FT3,FT4} <:  linAbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    lin_τ::FT 
    "Single scattering albedo"
    lin_ϖ::FT2   
    "Z scattering matrix (forward)"
    lin_Z⁺⁺::FT3 
    "Z scattering matrix (backward)"
    lin_Z⁻⁺::FT3
    "Ross kernel; cross section projection factor along µ (G ∈ [0,1], 1 for isotropic σ)"
    lin_G::FT4
end

Base.@kwdef struct linCoreAbsorptionOpticalProperties{FT} <:  linAbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    lin_τ::FT 
end

struct paired_xdx
    x::AbstractOpticalProperties #{xFT, xFT2, xFT3}
    dx::linAbstractOpticalProperties #{dxFT, dxFT2, dxFT3}
    #x::Vector{AbstractOpticalProperties} #{xFT, xFT2, xFT3}
    #dx::Vector{linAbstractOpticalProperties} #{dxFT, dxFT2, dxFT3}
end

struct scatt_paired_xdx
    x::CoreScatteringOpticalProperties #{xFT, xFT2, xFT3}
    dx::linCoreScatteringOpticalProperties #{dxFT, dxFT2, dxFT3}
    #x::Vector{CoreScatteringOpticalProperties} #{xFT, xFT2, xFT3}
    #dx::Vector{linCoreScatteringOpticalProperties} #{dxFT, dxFT2, dxFT3}
end

struct dirscatt_paired_xdx
    x::CoreDirectionalScatteringOpticalProperties #{xFT, xFT2, xFT3}
    dx::linCoreDirectionalScatteringOpticalProperties #{dxFT, dxFT2, dxFT3}
    #x::Vector{CoreDirectionalScatteringOpticalProperties} #{xFT, xFT2, xFT3}
    #dx::Vector{linCoreDirectionalScatteringOpticalProperties} #{dxFT, dxFT2, dxFT3}
end

struct abs_paired_xdx
    x::CoreAbsorptionOpticalProperties #{xFT, xFT2, xFT3}
    dx::linCoreAbsorptionOpticalProperties #{dxFT, dxFT2, dxFT3}
    #x::Vector{CoreAbsorptionOpticalProperties} #{xFT, xFT2, xFT3}
    #dx::Vector{linCoreAbsorptionOpticalProperties} #{dxFT, dxFT2, dxFT3}
end

# Adding Core Optical Properties, can have mixed dimensions!
function Base.:+(z1::scatt_paired_xdx, 
    z2::scatt_paired_xdx)# where {xFT, xFT2, xFT3, yFT, yFT2, yFT3}     
    #( (x::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
    #              dx::linCoreScatteringOpticalProperties{dxFT, dxFT2, dxFT3}), 
    #              (y::CoreScatteringOpticalProperties{yFT, yFT2, yFT3}, 
    #              dy::linCoreScatteringOpticalProperties{dyFT, dyFT2, dyFT3}) 
    #            ) where {xFT, xFT2, xFT3, yFT, yFT2, yFT3} 
    # Predefine some arrays:
    x = z1.x
    dx = z1.dx 
    y = z2.x
    dy = z2.dx 

    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺
    yZ⁺⁺ = y.Z⁺⁺
    yZ⁻⁺ = y.Z⁻⁺

    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ
    wy = y.τ .* y.ϖ  
    @show size(wx), size(wy), size(τ)
    ϖ  =  (wx .+ wy) ./ τ #(wx + wy) ./ τ

    @show size(dx.lin_τ), size(dy.lin_τ)
    nparams = size(dx.lin_τ,2)
    lin_τ = similar(dx.lin_τ)
    lin_ϖ = similar(dx.lin_τ)
    #for i = 1:nparams
    lin_τ  = dx.lin_τ .+ dy.lin_τ #dx.lin_τ + dy.lin_τ
    @show size(lin_ϖ), size(dx.lin_τ), size(x.ϖ), size(dx.lin_ϖ * x.τ')
    @show size(dy.lin_τ * y.ϖ), size(dy.lin_ϖ * y.τ')
    @show size(lin_τ .* ϖ')
    lin_ϖ  = ((dx.lin_τ * x.ϖ + 
                dx.lin_ϖ * x.τ') .+
            (dy.lin_τ * y.ϖ + 
            dy.lin_ϖ * y.τ') .- 
            lin_τ .* ϖ')./τ'
    #end
    
    #@show xFT, xFT2, xFT3
    #all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, dy.lin_Z⁺⁺, dy.lin_Z⁻⁺)) : nothing
    #all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, dx.lin_Z⁺⁺, dx.lin_Z⁻⁺)) : nothing

    n = length(τ);
    
    wx = reshape(wx,1,1,n)
    #wy = reshape(wy,1,1,n)
        
    Z⁺⁺ = (xZ⁺⁺ .* wx .+ wy .* yZ⁺⁺)./(wx .+ wy)
    Z⁻⁺ = (xZ⁻⁺ .* wx .+ wy .* yZ⁻⁺)./(wx .+ wy)

    tmpx1 = (dx.lin_τ .* x.ϖ .+ dx.lin_ϖ .* x.τ')
    tmpy1 = (dy.lin_τ .* y.ϖ .+ dy.lin_ϖ .* y.τ')
    tmpx2 = reshape(wx, 1,1,1,n)#x.τ * x.ϖ #wx
    tmpy2 = wy #y.τ * y.ϖ #wy
    @show size(Z⁺⁺)
    lin_Z⁺⁺ = zeros(nparams,size(Z⁺⁺,1), size(Z⁺⁺,2), size(Z⁺⁺,3))
    @show '1'
    #ylin_Z⁺⁺ = zeros(nparams,size(Z⁺⁺,1), size(Z⁺⁺,2), size(Z⁺⁺,3))
    #@show '2'
    lin_Z⁻⁺ = zeros(nparams,size(Z⁺⁺,1), size(Z⁺⁺,2), size(Z⁺⁺,3))
    @show '3'
    #ylin_Z⁻⁺ = zeros(nparams,size(Z⁺⁺,1), size(Z⁺⁺,2), size(Z⁺⁺,3))
    #@show '4'
    for i = 1:nparams
        lin_Z⁺⁺[i,:,:,:] = (x.Z⁺⁺ .- Z⁺⁺).*reshape(tmpx1[i,:],1,1,n) .+ 
                        dx.dZ⁺⁺[i,:,:].*tmpx2 .+ 
                        y.Z⁺⁺.*reshape(tmpy1[i,:],1,1,n) .+ 
                        dy.dZ⁺⁺[i,:,:].*tmpy2 
        lin_Z⁻⁺[i,:,:,:] = (x.Z⁻⁺ .- Z⁻⁺).*reshape(tmpx1[i,:],1,1,n) .+ 
                        dx.dZ⁻⁺[i,:,:].*tmpx2 .+
                        y.Z⁻⁺.*reshape(tmpx1[i,:],1,1,n) .+ 
                        dy.dZ⁻⁺[i,:,:].*tmpy2 
    end
    lin_Z⁺⁺ = lin_Z⁺⁺./reshape(wx .+ wy, 1,1,1,n)
    lin_Z⁻⁺ = lin_Z⁻⁺./reshape(wx .+ wy, 1,1,1,n)

    return scatt_paired_xdx(CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, lin_Z⁺⁺, lin_Z⁻⁺)) 
end

# Concatenate Core Optical Properties, can have mixed dimensions!
# Check with Christian what exactly this is supposed to do
#=
function Base.:*(x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties, dx::linCoreScatteringOpticalProperties, dy::linCoreScatteringOpticalProperties) 
    arr_type  = array_type(architecture(x.τ))
    x = expandOpticalProperties(x,arr_type);
    y = expandOpticalProperties(y,arr_type);
    
    dx = expandOpticalProperties(dx,arr_type);
    dy = expandOpticalProperties(dy,arr_type);
    
    CoreScatteringOpticalProperties([x.τ; y.τ],[x.ϖ; y.ϖ],cat(x.Z⁺⁺,y.Z⁺⁺, dims=3), cat(x.Z⁻⁺,y.Z⁻⁺, dims=3) ), 
    linCoreScatteringOpticalProperties(cat(dx.dτ, dy.dτ, dims=2),cat(dx.dϖ, dy.dϖ, dims=2), cat(dx.dZ⁺⁺,dy.dZ⁺⁺, dims=4), cat(dx.dZ⁻⁺,dy.dZ⁻⁺, dims=4) )
end
=#

function Base.:+(z1::scatt_paired_xdx, z2::abs_paired_xdx)
    #(x::CoreScatteringOpticalProperties, dx::linCoreScatteringOpticalProperties, y::CoreAbsorptionOpticalProperties, dy::linCoreAbsorptionOpticalProperties) 
    # Predefine some arrays: 
    x = z1.x
    dx = z1.dx 
    y = z2.x
    dy = z2.dx 

    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ
    ϖ  =  wx ./ τ

    lin_τ  = dx.lin_τ + dy.lin_τ
    lin_ϖ  = (dx.lin_τ .* x.ϖ + x.τ .* dx.lin_ϖ - 
            ϖ .* lin_τ)./τ
    
    
    #@show xFT, xFT2, xFT3
    #all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, dy.lin_Z⁺⁺, dy.lin_Z⁻⁺)) : nothing
    #all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, dx.lin_Z⁺⁺, dx.lin_Z⁻⁺)) : nothing
    
    return scatt_paired_xdx(CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, dx.lin_Z⁺⁺, dx.lin_Z⁻⁺))
end

function Base.:+(z2::abs_paired_xdx, z1::scatt_paired_xdx) #(y::CoreAbsorptionOpticalProperties, dy::linCoreAbsorptionOpticalProperties, x::CoreScatteringOpticalProperties, dx::linCoreScatteringOpticalProperties) 
    x = z1.x
    dx = z1.dx 
    y = z2.x
    dy = z2.dx 
    
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ
    ϖ  =  wx ./ τ

    lin_τ  = dx.lin_τ + dy.lin_τ
    lin_ϖ  = (dx.lin_τ .* x.ϖ + x.τ .* dx.lin_ϖ - 
            ϖ .* lin_τ)./τ

    return scatt_paired_xdx(CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺), linCoreScatteringOpticalProperties(lin_τ, lin_ϖ, dx.lin_Z⁺⁺, dx.lin_Z⁻⁺))
end

#=
function Base.:*( x::FT, y::CoreScatteringOpticalProperties{FT}, dy::linCoreScatteringOpticalProperties{FT} ) where FT
    CoreScatteringOpticalProperties(y.τ * x, y.ϖ, y.Z⁺⁺, y.Z⁻⁺, y.G), linCoreScatteringOpticalProperties(dy.lin_τ * x, dy.lin_ϖ, dy.linZ⁺⁺, dy.linZ⁻⁺, dy.linG)
end
=#
