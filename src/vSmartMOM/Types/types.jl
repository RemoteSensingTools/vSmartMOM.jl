
"Struct for an atmospheric profile"
struct AtmosphericProfile{FT}
    lat::FT
    lon::FT
    psurf::FT
    T::Array{FT,1}
    q::Array{FT,1}
    p::Array{FT,1}
    p_levels::Array{FT,1}
    vmr_h2o::Array{FT,1}
    vcd_dry::Array{FT,1}
    vcd_h2o::Array{FT,1}
end

# Types for describing atmospheric parameters #
abstract type AbstractObsGeometry end

"Observation Geometry (basics)" 
@with_kw struct ObsGeometry{FT} <: AbstractObsGeometry
    "altitude of observer `[Pa]`"
    obs_alt::FT
    "solar zenith angle `[Degree]`"
    sza::FT
    "viewing zenith angle(s) `[Degree]`" 
    vza::Array{FT,1}
    "viewing azimuth angle(s) `[Degree]`" 
    vaz::Array{FT,1}
end

"Quadrature Types for RT streams"
abstract type AbstractQuadratureType end

struct RadauQuad <:AbstractQuadratureType end
struct GaussQuadHemisphere <:AbstractQuadratureType end
struct GaussQuadFullSphere <:AbstractQuadratureType end


abstract type AbstractLayer end

@with_kw struct CompositeLayer{FT} <: AbstractLayer 

    # Composite layer R and T matrices
    R⁻⁺::AbstractArray{FT,3}
    R⁺⁻::AbstractArray{FT,3}
    T⁺⁺::AbstractArray{FT,3}
    T⁻⁻::AbstractArray{FT,3}

end

@with_kw struct AddedLayer{FT} <: AbstractLayer 
    
    # Added layer R and T matrices
    r⁻⁺::AbstractArray{FT,3}
    t⁺⁺::AbstractArray{FT,3}
    r⁺⁻::AbstractArray{FT,3}
    t⁻⁻::AbstractArray{FT,3}

end

abstract type AbstractScatteringInterface end
struct ScatteringInterface_00 <: AbstractScatteringInterface end
struct ScatteringInterface_01 <: AbstractScatteringInterface end
struct ScatteringInterface_10 <: AbstractScatteringInterface end
struct ScatteringInterface_11 <: AbstractScatteringInterface end

abstract type AbstractvSmartMOMModel end

"""
    struct vSmartMOM_Parameters

A struct which holds all initial model parameters (before any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Parameters{FT} #<: AbstractvSmartMOMModel

    "Float type to use in the RT (Float64/Float32)"
    float_type::DataType
    "Incident wavelength"
    λ::FT
    "Depolarization value"
    depol::FT
    "Truncation length for legendre terms"
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
    "Log mean radius"
    μ::FT
    "Log stddev of radius"
    σ::FT
    "Maximum radius"
    r_max::FT
    "Number of quadrature points for integration of size distribution"
    nquad_radius::Integer
    "Real part of refractive index"
    nᵣ::FT
    "Imag part of refractive index"
    nᵢ::FT
    "Pressure peak [Pa]"
    p₀::FT
    "Pressure peak width [Pa]"
    σp::FT
    "Path to atmospheric profile file"
    file::AbstractString
    "Time index to retrieve from file "
    timeIndex::Integer
    "Latitude of atmospheric profile"
    lat::Number
    "Longitude of atmospheric profile"
    lon::Number
    "Length of profile reduction"
    profile_reduction_n::Integer
    "Starting wavenumber for absorption grid"
    ν_min::FT
    "Ending wavenumber for absorption grid"
    ν_max::FT
    "Number of spectral points in absorption grid"
    ν_step::FT
    "Type of broadening function (Doppler/Lorentz/Voigt)"
    broadening_function::AbstractBroadeningFunction
    "Wing cutoff to use in cross-section calculation"
    wing_cutoff::Integer
    "Complex Error Function to use in Voigt calculations"
    CEF::AbstractComplexErrorFunction
    "Volume mixing ratio"
    vmr::FT
    "Algorithm to use for fourier decomposition (NAI2/PCW)"
    decomp_type::AbstractFourierDecompositionType
    "Quadrature type for RT streams (RadauQuad/GaussQuadHemisphere/GaussQuadFullSphere)"
    quadrature_type::AbstractQuadratureType
    "Architecture to use for calculations (CPU/GPU)"
    architecture::AbstractArchitecture
end

"""
    struct vSmartMOM_Model

A struct which holds all derived model parameters (including any computations)

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Model <: AbstractvSmartMOMModel

    "Struct with all individual parameters"
    params::vSmartMOM_Parameters
    "Abstract greek coefficient truncation type "
    truncation_type::AbstractTruncationType
    "Mie model to use in phase function calculations"
    mie_model::MieModel
    "Computed aerosol optics"
    aerosol_optics::AerosolOptics
    "Computed aerosol optics (truncated)"
    aerosol_optics_truncated::AerosolOptics
    "Greek coefs in Rayleigh calculations" 
    greek_rayleigh::GreekCoefs
    "Number of quadrature points"
    Nquad::Integer
    "Quadrature points"
    qp_μ::AbstractArray
    "Quadrature weights"
    wt_μ::AbstractArray
    "Rayleigh "
    τRayl::AbstractArray
    "Rayleigh optical thickness"
    ϖRayl::AbstractArray
    "Aerosols optical thickness"
    τAer::AbstractArray
    "Aerosols albedo"
    ϖAer::AbstractArray
    "Observational Geometry (includes sza, vza, vaz)"
    obs_geom::ObsGeometry
    "Aerosol type"
    aerosol::AbstractAerosolType
    "Atmospheric profile to use"
    profile::AtmosphericProfile
    "Absorption model for computing absorption cross sections"
    absorption_model::AbstractCrossSectionModel
    "Array to hold cross-sections over entire atmospheric profile"
    τ_abs::AbstractArray

end