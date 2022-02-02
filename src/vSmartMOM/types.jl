#=
 
This file contains all types that are used in the vSmartMOM module:

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
struct AtmosphericProfile{FT, VMR <: Union{Real, Vector}}
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
    vmr::Dict{String, VMR}
end

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

mutable struct RT_Aerosol{}#FT<:Union{AbstractFloat, ForwardDiff.Dual}}
    "Aerosol"
    aerosol::Aerosol#{FT}
    "Reference τ"
    τ_ref#::FT
    "Pressure peak (Pa)"
    p₀#::FT
    "Pressure peak width (Pa)"
    σp#::FT
end

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

"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct CompositeLayer{FT} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    R⁻⁺::AbstractArray{FT,3}
    "Composite layer Reflectance matrix R (from - -> +)"
    R⁺⁻::AbstractArray{FT,3}
    "Composite layer transmission matrix T (from + -> +)"
    T⁺⁺::AbstractArray{FT,3}
    "Composite layer transmission matrix T (from - -> -)"
    T⁻⁻::AbstractArray{FT,3}
    "Composite layer source matrix J (in + direction)"
    J₀⁺::AbstractArray{FT,3}
    "Composite layer source matrix J (in - direction)"
    J₀⁻::AbstractArray{FT,3}
end

"Added (Single) Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct AddedLayer{FT} <: AbstractLayer 
    "Added layer Reflectance matrix R (from + -> -)"
    r⁻⁺::AbstractArray{FT,3}
    "Added layer transmission matrix T (from + -> +)"
    t⁺⁺::AbstractArray{FT,3}
    "Added layer Reflectance matrix R (from - -> +)"
    r⁺⁻::AbstractArray{FT,3}
    "Added layer transmission matrix T (from - -> -)"
    t⁻⁻::AbstractArray{FT,3}
    "Added layer source matrix J (in + direction)"
    J₀⁺::AbstractArray{FT,3}
    "Added layer source matrix J (in - direction)"
    J₀⁻::AbstractArray{FT,3}
end

"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
struct CompositeLayerRS{FT} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    R⁻⁺::AbstractArray{FT,3}
    "Composite layer Reflectance matrix R (from - -> +)"
    R⁺⁻::AbstractArray{FT,3}
    "Composite layer transmission matrix T (from + -> +)"
    T⁺⁺::AbstractArray{FT,3}
    "Composite layer transmission matrix T (from - -> -)"
    T⁻⁻::AbstractArray{FT,3}
    "Composite layer source matrix J (in + direction)"
    J₀⁺::AbstractArray{FT,3}
    "Composite layer source matrix J (in - direction)"
    J₀⁻::AbstractArray{FT,3}

    # Additional Arrays for Raman scattering
    "Composite layer Reflectance matrix ieR (from + -> -)"
    ieR⁻⁺::AbstractArray{FT,4}
    "Composite layer Reflectance matrix ieR (from - -> +)"
    ieR⁺⁻::AbstractArray{FT,4}
    "Composite layer transmission matrix ieT (from + -> +)"
    ieT⁺⁺::AbstractArray{FT,4}
    "Composite layer transmission matrix ieT (from - -> -)"
    ieT⁻⁻::AbstractArray{FT,4}
    "Composite layer source matrix ieJ (in + direction)"
    ieJ₀⁺::AbstractArray{FT,4}
    "Composite layer source matrix ieJ (in - direction)"
    ieJ₀⁻::AbstractArray{FT,4}
end

"Added (Single) Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
struct AddedLayerRS{FT} <: AbstractLayer 
    "Added layer Reflectance matrix R (from + -> -)"
    r⁻⁺::AbstractArray{FT,3}
    "Added layer transmission matrix T (from + -> +)"
    t⁺⁺::AbstractArray{FT,3}
    "Added layer Reflectance matrix R (from - -> +)"
    r⁺⁻::AbstractArray{FT,3}
    "Added layer transmission matrix T (from - -> -)"
    t⁻⁻::AbstractArray{FT,3}
    "Added layer source matrix J (in + direction)"
    J₀⁺::AbstractArray{FT,3}
    "Added layer source matrix J (in - direction)"
    J₀⁻::AbstractArray{FT,3}

    # Additional Arrays for Raman scattering
    "Added layer Reflectance matrix ieR (from + -> -)"
    ier⁻⁺::AbstractArray{FT,4}
    "Added layer transmission matrix ieT (from + -> +)"
    iet⁺⁺::AbstractArray{FT,4}
    "Added layer Reflectance matrix ieR (from - -> +)"
    ier⁺⁻::AbstractArray{FT,4}
    "Added layer transmission matrix ieT (from - -> -)"
    iet⁻⁻::AbstractArray{FT,4}
    "Added layer source matrix ieJ (in + direction)"
    ieJ₀⁺::AbstractArray{FT,4}
    "Added layer source matrix ieJ (in - direction)"
    ieJ₀⁻::AbstractArray{FT,4}
end

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

"Defined by Legendre polynomial terms as function of spectral grid, which is scaled to [-1,1] (degree derived from length of `a_coeff`)"
struct LambertianSurfaceLegendre{FT} <: AbstractSurfaceType
    "albedo = legendre_coeff[1] * P₀ + legendre_coeff[2]*P₁ + legendre_coeff[3]*P₂ + ... "
    legendre_coeff::AbstractArray{FT,1}
end

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
end

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
    "Algorithm to use for fourier decomposition (NAI2/PCW)"
    decomp_type::AbstractFourierDecompositionType
end

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

"""
    struct vSmartMOM_Model

A struct which holds all derived model parameters (including any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Model

    "Struct with all individual parameters"
    params::vSmartMOM_Parameters
    
    "Truncated aerosol optics"
    aerosol_optics::AbstractArray{AbstractArray{AerosolOptics}}
    "Greek coefs in Rayleigh calculations" 
    greek_rayleigh::GreekCoefs
    "Quadrature points/weights, etc"
    quad_points::QuadPoints

    "Array to hold cross-sections over entire atmospheric profile"
    τ_abs::AbstractArray{AbstractArray}
    "Rayleigh optical thickness"
    τ_rayl::AbstractArray{AbstractArray}
    "Aerosol optical thickness"
    τ_aer::AbstractArray{AbstractArray}

    "Observational Geometry (includes sza, vza, vaz)"
    obs_geom::ObsGeometry
    "Atmospheric profile to use"
    profile::AtmosphericProfile
end

"""
    struct ComputedAtmosphereProperties

A struct which holds (for the entire atmosphere) all key layer optical properties required for the RT core solver

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ComputedAtmosphereProperties

    "Absorption optical depth vectors (wavelength dependent)"
    τ_λ_all
    "Albedo vectors (wavelength dependent)"
    ϖ_λ_all
    "Absorption optical depth scalars (not wavelength dependent)"
    τ_all
    "Albedo scalars (not wavelength dependent)"
    ϖ_all
    "Combined Z moments (forward)"
    Z⁺⁺_all
    "Combined Z moments (backward)"
    Z⁻⁺_all
    "Maximum dτs"
    dτ_max_all
    "dτs"
    dτ_all
    "Number of doublings (for all layers)"
    ndoubl_all
    "dτs (wavelength dependent)"
    dτ_λ_all
    "All expk"
    expk_all
    "Scattering flags"
    scatter_all
    "Sum of optical thicknesses of all layers above the current layer"
    τ_sum_all
    #"elastic (Cabannes) scattering fraction of Rayleigh (Cabannes+Raman) scattering per layer"
    #ϖ_Cabannes_all
    "Rayleigh fraction of scattering cross section per layer"
    fscattRayl_all
    "Scattering interface type for each layer"
    scattering_interfaces_all
end

# TODO SUNITI: write a function to compute these properties and create this structure
#Base.@kwdef struct RamanAtmosphereProperties
    #"band spectral grid"
    #grid_in
    #"inelastic scattering SSA"
    #ϖ_λ₀λ₁
    #"inelastic scattering index"
    #i_λ₀λ₁
    #"inelastic (vibrational) scattering SSA: split later for each molecule"
    #ϖ_vib_λ₀λ₁
    #"inelastic (vibrational) scattering index: split later for each molecule"
    #i_vib_λ₀λ₁
    #"Greek coefs in Rayleigh calculations" 
    #greek_raman::GreekCoefs
    #"Combined o2 and n2 Z moments for rotational/rovibrational RS  (forward)"
    #Z⁺⁺_RRS #same for rotational and rovibrational scattering
    #"Combined o2 and n2 Z moments for rotational/rovibrational RS  (backward)"
    #Z⁻⁺_RRS #same for rotational and rovibrational scattering
    #"Combined o2 and n2 Z moments for vibrational RS (forward): split later for each molecule"
    #Z⁺⁺_VRS #same for rotational and rovibrational scattering
    #"Combined o2 and n2 Z moments for vibrational RS (backward): split later for each molecule"
    #Z⁻⁺_VRS #same for rotational and rovibrational scattering
#end

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
Base.@kwdef mutable struct RRS{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constants for N2"
    n2::InelasticScattering.MolecularConstants{Float64}
    "Molecular Constants for O2"
    o2::InelasticScattering.MolecularConstants{Float64}
    "Greek coeffs in Raman calculations" 
    greek_raman::GreekCoefs
    "Pre-computed optical properties"
    # ramanAtmoProp::RamanAtmosphereProperties
    fscattRayl::FT
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
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
Base.@kwdef struct VS_0to1{FT<:AbstractFloat} <: AbstractRamanType 
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
Base.@kwdef struct VS_1to0{FT<:AbstractFloat} <: AbstractRamanType 
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
Base.@kwdef struct RVRS{FT<:AbstractFloat} <: AbstractRamanType 
    "Molecular Constant for N2"
    n2::InelasticScattering.MolecularConstants{FT}
    "Molecular Constant for O2"
    o2::InelasticScattering.MolecularConstants{FT}
    "Pre-computed optical properties"
    #ramanAtmoProp::RamanAtmosphereProperties
    fscattRayl::FT
    ϖ_Cabannes::FT #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
    ϖ_λ₁λ₀::Array{FT,1}
    i_λ₁λ₀::Array{Int,1}
    Z⁻⁺_λ₁λ₀::Array{FT,2}
    Z⁺⁺_λ₁λ₀::Array{FT,2}
    τ₀::FT
    n_Raman::Int
end

Base.@kwdef mutable struct noRS{FT} <: AbstractRamanType
    fscattRayl::FT = 0.0
    ϖ_Cabannes::FT = 1.0 #elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering
end



"""
    struct ComputedLayerProperties

A struct which holds all key layer optical properties required for the RT core solver

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ComputedLayerProperties

    "Absorption optical depth vector (wavelength dependent)"
    τ_λ 
    "Albedo vector (wavelength dependent)"
    ϖ_λ 
    "Absorption optical depth scalar (not wavelength dependent)"
    τ 
    "Albedo scalar (not wavelength dependent)"
    ϖ  
    "Combined Z moment (forward)"
    Z⁺⁺ 
    "Combined Z moment (backward)"
    Z⁻⁺ 
    "Maximum dτ"
    dτ_max 
    "dτ"
    dτ     
    "Number of doublings"
    ndoubl
    "dτ (wavelength dependent)"
    dτ_λ 
    "expk"
    expk 
    "Scattering flag"
    scatter 
    "Sum of optical thicknesses of all layers above the current layer"
    τ_sum
    "Fraction of scattering caused by Rayleigh"
    fscattRayl
    #"Elastic fraction (Cabannes) of Rayleigh (Cabannes+Raman) scattering"
    #ϖ_Cabannes 
    "Scattering interface type for current layer"
    scattering_interface
end

abstract type AbstractOpticalProperties end

# Core optical Properties COP
Base.@kwdef struct CoreScatteringOpticalProperties{FT} <:  AbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    τ::Union{FT, AbstractArray{FT,1}} 
    "Single scattering albedo"
    ϖ::Union{FT, AbstractArray{FT,1}}   
    "Z scattering matrix (forward)"
    Z⁺⁺::Union{AbstractArray{FT,2}, AbstractArray{FT,3}, Nothing}  
    "Z scattering matrix (backward)"
    Z⁻⁺::Union{AbstractArray{FT,2}, AbstractArray{FT,3}, Nothing}
end

Base.@kwdef struct CoreAbsorptionOpticalProperties{FT} <:  AbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    τ::Union{FT, AbstractArray{FT,1}} 
end

# Adding Core Optical Properties, can have mixed dimensions!
function Base.:+( x::CoreScatteringOpticalProperties{FT}, y::CoreScatteringOpticalProperties{FT} ) where FT
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    wy = y.τ .* y.ϖ  
    w  = wx .+ wy
    ϖ  =  (w) ./ τ
    
    all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)) : nothing
    all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)) : nothing

    # A bit more tediuous for Z matrices:
    length(unique(w))  == 1 ? w = unique(w) : nothing
    length(unique(wx)) == 1 ? wx = unique(wx) : nothing
    length(unique(wy)) == 1 ? wy = unique(wy) : nothing
    n = length(w);
    if n == 1
        Z⁺⁺ = ((wx .* x.Z⁺⁺) .+ (wy .* y.Z⁺⁺)) ./ w 
        Z⁻⁺ = ((wx .* x.Z⁻⁺) .+ (wy .* y.Z⁻⁺)) ./ w
    else
        nx = length(wx)
        ny = length(wy)
        nx > 1 ? wx = reshape(wx,1,1,n) : nothing
        ny > 1 ? wy = reshape(wy,1,1,n) : nothing
        w = reshape(w,1,1,n)
        size(x.Z⁺⁺,3) == 1 ? xZ⁺⁺ = repeat(x.Z⁺⁺,1,1,n) : xZ⁺⁺ = x.Z⁺⁺
        size(x.Z⁻⁺,3) == 1 ? xZ⁻⁺ = repeat(x.Z⁻⁺,1,1,n) : xZ⁻⁺ = x.Z⁻⁺
        size(y.Z⁺⁺,3) == 1 ? yZ⁺⁺ = repeat(y.Z⁺⁺,1,1,n) : yZ⁺⁺ = y.Z⁺⁺
        size(y.Z⁻⁺,3) == 1 ? yZ⁻⁺ = repeat(y.Z⁻⁺,1,1,n) : yZ⁻⁺ = y.Z⁻⁺
        Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) ./ w 
        Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺) ./ w
    end
    CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺)  
end

# Concatenate Core Optical Properties, can have mixed dimensions!
function Base.:*( x::CoreScatteringOpticalProperties{FT}, y::CoreScatteringOpticalProperties{FT} ) where FT
    arr_type  = array_type(architecture(x.τ))

    x = expandOpticalProperties(x,arr_type);
    y = expandOpticalProperties(y,arr_type);
    CoreScatteringOpticalProperties([x.τ; y.τ],[x.ϖ; y.ϖ],cat(x.Z⁺⁺,y.Z⁺⁺, dims=3), cat(x.Z⁻⁺,y.Z⁻⁺, dims=3))
end

function Base.:+( x::CoreScatteringOpticalProperties{FT}, y::CoreAbsorptionOpticalProperties{FT} ) where FT
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    ϖ  = (wx) ./ τ
    CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)
end

function Base.:+(  y::CoreAbsorptionOpticalProperties{FT}, x::CoreScatteringOpticalProperties{FT} ) where FT
    x + y
end


function Base.:*( x::FT, y::CoreScatteringOpticalProperties{FT} ) where FT
    CoreScatteringOpticalProperties(y.τ * x, y.ϖ, y.Z⁺⁺, y.Z⁻⁺)
end

