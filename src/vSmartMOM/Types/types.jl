
"Struct for an atmospheric profile"
struct AtmosphericProfile{FT}
    "Latitude of the AtmosphericProfile"
    lat
    "Longitude of the AtmosphericProfile"
    lon
    "Surface Pressure"
    psurf::FT
    "Temperature Profile"
    T::Array{FT,1}
    "Specific humidity profile"
    q::Array{FT,1}
    "Pressure Profile"
    p::Array{FT,1}
    "Pressure Levels"
    p_levels::Array{FT,1}
    "H2O Volume Mixing Ratio Profile"
    vmr_h2o::Array{FT,1}
    "Vertical Column Density (Dry)"
    vcd_dry::Array{FT,1}
    "Vertical Column Density (H2O)"
    vcd_h2o::Array{FT,1}
end

"Types for describing atmospheric parameters"
abstract type AbstractObsGeometry end

"Observation Geometry (basics)" 
@with_kw struct ObsGeometry{FT} <: AbstractObsGeometry
    "Altitude of Observer `[Pa]`"
    obs_alt::FT
    "Solar Zenith Angle `[Degree]`"
    sza::FT
    "Viewing Zenith Angle(s) `[Degree]`" 
    vza::Array{FT,1}
    "Viewing Azimuth Angle(s) `[Degree]`" 
    vaz::Array{FT,1}
end

"Quadrature Types for RT streams"
abstract type AbstractQuadratureType end

# <<Christian>>
struct RadauQuad <:AbstractQuadratureType end
struct GaussQuadHemisphere <:AbstractQuadratureType end
struct GaussQuadFullSphere <:AbstractQuadratureType end

"Abstract Type for Source Function Integration"
abstract type AbstractSourceType end
"Dummy Node Integration"
struct DNI <:AbstractSourceType end
"Source Function Integration"
struct SFI <:AbstractSourceType end

"Abstract Type for Layer R and T matrices"
abstract type AbstractLayer end

"Composite Layer Matrices"
@with_kw mutable struct CompositeLayer{FT} <: AbstractLayer 

    "Composite R (Forward) <<Suniti>>"
    R⁻⁺::AbstractArray{FT,3}
    "Composite R (Backward) <<Suniti>>"
    R⁺⁻::AbstractArray{FT,3}
    "Composite T (Forward) <<Suniti>>"
    T⁺⁺::AbstractArray{FT,3}
    "Composite T (Backward) <<Suniti>>"
    T⁻⁻::AbstractArray{FT,3}
    "Composite J (Forward) <<Suniti>>"
    J₀⁺::AbstractArray{FT,3}
    "Composite J (Backward) <<Suniti>>"
    J₀⁻::AbstractArray{FT,3}
end

"Added (Single) Layer Matrices"
@with_kw mutable struct AddedLayer{FT} <: AbstractLayer 
    
    "Added R (Forward) <<Suniti>>"
    r⁻⁺::AbstractArray{FT,3}
    "Added T (Forward) <<Suniti>>"
    t⁺⁺::AbstractArray{FT,3}
    "Added R (Backward) <<Suniti>>"
    r⁺⁻::AbstractArray{FT,3}
    "Added T (Backward) <<Suniti>>"
    t⁻⁻::AbstractArray{FT,3}
    "Added J (Forward) <<Suniti>>"
    J₀⁺::AbstractArray{FT,3}
    "Added J (Backward) <<Suniti>>"
    J₀⁻::AbstractArray{FT,3}
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

"Defined by polynomial terms as function of λ (degree derived from length of `a_coeff`)"
struct LambertianSurfacePolyFit{FT} <: AbstractSurfaceType
    "albedo(λ) = a_coeff[1] + a_coeff[2]*λ + a_coeff[3]*λ² + ... "
    a_coeff::AbstractArray{FT,1}
end

"""
    struct vSmartMOM_Parameters

A struct which holds all initial model parameters (before any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Parameters{FT<:Union{AbstractFloat, ForwardDiff.Dual}} 

    "Float type to use in the RT (Float64/Float32)"
    float_type::DataType
    "band center wavelength (`nBand`)"
    λ_band::AbstractArray{FT}
    "reference wavelength"
    λ_ref::FT
    "Depolarization value"
    depol::FT
    "Truncation length for legendre terms (scalar for now, can do `nBand` later)"
    l_trunc::Integer
    "Exclusion angle for forward peak"
    Δ_angle::FT
    "Hard cutoff for maximum number of m iterations"
    max_m::Integer
    "Type of polarization (I/IQ/IQU/IQUV)"
    polarization_type::AbstractPolarizationType
    "Altitude of observer `[Pa]`"
    obs_alt::FT
    "Solar zenith angle [deg]"
    sza::FT
    "Viewing zenith angles [deg]"
    vza::AbstractArray{FT}
    "Viewing azimuthal angles [deg]"
    vaz::AbstractArray{FT}
    "Number of aerosol species"
    nAer::Integer
    "AOD at Reference wavelength (`nAer`)"
    τAer_ref::AbstractArray#{FT} #Suniti
    "Log mean radius (`nAer`)"
    μ::AbstractArray{FT} #Suniti
    "Log stddev of radius (`nAer`)"
    σ::AbstractArray{FT} #Suniti
    "Maximum radius"
    r_max::FT
    "Number of quadrature points for integration of size distribution"
    nquad_radius::Integer
    "Real part of refractive index (`nBand,nAer`)"
    nᵣ::AbstractArray#{FT}
    "Imag part of refractive index (`nBand,nAer`)"
    nᵢ::AbstractArray{FT}
    "Pressure peak [Pa] (`nAer`)"
    p₀::AbstractArray#{FT}
    "Pressure peak width [Pa] (`nAer`)"
    σp::AbstractArray{FT}
    "Path to atmospheric profile file"
    file::AbstractString
    "Time index to retrieve from file "
    timeIndex::Union{Integer, Nothing}
    "Latitude of atmospheric profile"
    lat::Union{Number, Nothing}
    "Longitude of atmospheric profile"
    lon::Union{Number, Nothing}
    "Length of profile reduction"
    profile_reduction_n::Integer
    "Starting wavelength for absorption grid (`nBand`)"
    spec_grid_start::AbstractArray{FT}
    "Ending wavelength for absorption grid (`nBand`)"
    spec_grid_end::AbstractArray{FT}
    "Number of spectral points in absorption grid (`nBand`)"
    spec_grid_n::AbstractArray{Integer}
    "Type of broadening function (Doppler/Lorentz/Voigt)"
    broadening_function::AbstractBroadeningFunction
    "Wing cutoff to use in cross-section calculation"
    wing_cutoff::Integer
    "Complex Error Function to use in Voigt calculations"
    CEF::AbstractComplexErrorFunction
    "Volume mixing ratio (needs to be per profile, gas species and band later)"
    vmr::FT
    "Surface reflectance type"
    brdf#::AbstractArray{AbstractSurfaceType}
    "Algorithm to use for fourier decomposition (NAI2/PCW)"
    decomp_type::AbstractFourierDecompositionType
    "Quadrature type for RT streams (RadauQuad/GaussQuadHemisphere/GaussQuadFullSphere)"
    quadrature_type::AbstractQuadratureType
    "Architecture to use for calculations (CPU/GPU)"
    architecture::AbstractArchitecture
end

"Quadrature points, weights, etc"
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
    aerosol_optics::AbstractArray #Array{Array{AerosolOptics}}
    "Greek coefs in Rayleigh calculations" 
    greek_rayleigh::GreekCoefs
    "Quadrature points/weights, etc"
    quad_points::QuadPoints
    "Array to hold cross-sections over entire atmospheric profile"
    τ_abs::AbstractArray
    "Rayleigh optical thickness"
    τRayl::AbstractArray
    "Aerosol optical thickness"
    τAer::AbstractArray
    "Observational Geometry (includes sza, vza, vaz)"
    obs_geom::ObsGeometry
    "Atmospheric profile to use"
    profile::AtmosphericProfile
    "Absorption model for computing absorption cross sections (will be an array later, for nGasSpecies,nBand)"
    #absorption_model::AbstractCrossSectionModel
    "Surface BRDF"
    brdf::AbstractSurfaceType
end

"A struct to internally hold the computed atmosphere properties <<Suniti>>"
abstract type ComputedProperties end
@with_kw struct ComputedAtmosphereProperties

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
    "Scattering interface type for each layer"
    scattering_interfaces_all
end

"A struct to internally hold the computed single-layer properties <<Suniti>>"
@with_kw struct ComputedLayerProperties

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
    "Scattering interface type for current layer"
    scattering_interface
end