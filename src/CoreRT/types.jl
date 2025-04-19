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
    "Layer height (meters)"
    Δz::Array{FT,1}
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
    "Vertical distribution as function of p (using Distributions.jl)"
    profile::Distribution#::FT
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
    j₀⁺::AbstractArray{FT,3}
    "Added layer source matrix J (in - direction)"
    j₀⁻::AbstractArray{FT,3}
    "Added layer temporary space to avoid allocations"
    temp1::Union{AbstractArray{FT,3}, Nothing}
    "Added layer temporary space to avoid allocations"
    temp2::Union{AbstractArray{FT,3}, Nothing}
    "Pointer to temporary space to avoid allocations"
    temp1_ptr::Union{CuArray{CuPtr{FT}, 1, CUDA.DeviceMemory}, Nothing}
    "Pointer to temporary space to avoid allocations"
    temp2_ptr::Union{CuArray{CuPtr{FT}, 1, CUDA.DeviceMemory}, Nothing}
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
# Multisensor Composite layers 
# Elastic
"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct CompositeLayerMS{M} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    topR⁻⁺::M#AbstractArray{FT,4}
    "Composite layer Reflectance matrix R (from - -> +)"
    topR⁺⁻::M
    "Composite layer transmission matrix T (from + -> +)"
    topT⁺⁺::M
    "Composite layer transmission matrix T (from - -> -)"
    topT⁻⁻::M
    "Composite layer source matrix J (in + direction)"
    topJ₀⁺::M
    "Composite layer source matrix J (in - direction)"
    topJ₀⁻::M
    "Composite layer Reflectance matrix R (from + -> -)"
    botR⁻⁺::M
    "Composite layer Reflectance matrix R (from - -> +)"
    botR⁺⁻::M
    "Composite layer transmission matrix T (from + -> +)"
    botT⁺⁺::M
    "Composite layer transmission matrix T (from - -> -)"
    botT⁻⁻::M
    "Composite layer source matrix J (in + direction)"
    botJ₀⁺::M
    "Composite layer source matrix J (in - direction)"
    botJ₀⁻::M
end

# Inelastic
"Composite Layer Matrices (`-/+` defined in τ coordinates, i.e. `-`=outgoing, `+`=incoming"
Base.@kwdef struct CompositeLayerMSRS{M1, M2} <: AbstractLayer 
    "Composite layer Reflectance matrix R (from + -> -)"
    topR⁻⁺::M1
    "Composite layer Reflectance matrix R (from - -> +)"
    topR⁺⁻::M1
    "Composite layer transmission matrix T (from + -> +)"
    topT⁺⁺::M1
    "Composite layer transmission matrix T (from - -> -)"
    topT⁻⁻::M1
    "Composite layer source matrix J (in + direction)"
    topJ₀⁺::M1 
    "Composite layer source matrix J (in - direction)"
    topJ₀⁻::M1

    # Additional Arrays for Raman scattering
    "Composite layer Reflectance matrix ieR (from + -> -)"
    topieR⁻⁺::M2
    "Composite layer Reflectance matrix ieR (from - -> +)"
    topieR⁺⁻::M2
    "Composite layer transmission matrix ieT (from + -> +)"
    topieT⁺⁺::M2
    "Composite layer transmission matrix ieT (from - -> -)"
    topieT⁻⁻::M2
    "Composite layer source matrix ieJ (in + direction)"
    topieJ₀⁺::M2
    "Composite layer source matrix ieJ (in - direction)"
    topieJ₀⁻::M2

    "Composite layer Reflectance matrix R (from + -> -)"
    botR⁻⁺::M1
    "Composite layer Reflectance matrix R (from - -> +)"
    botR⁺⁻::M1
    "Composite layer transmission matrix T (from + -> +)"
    botT⁺⁺::M1
    "Composite layer transmission matrix T (from - -> -)"
    botT⁻⁻::M1
    "Composite layer source matrix J (in + direction)"
    botJ₀⁺::M1
    "Composite layer source matrix J (in - direction)"
    botJ₀⁻::M1

    # Additional Arrays for Raman scattering
    "Composite layer Reflectance matrix ieR (from + -> -)"
    botieR⁻⁺::M2
    "Composite layer Reflectance matrix ieR (from - -> +)"
    botieR⁺⁻::M2
    "Composite layer transmission matrix ieT (from + -> +)"
    botieT⁺⁺::M2
    "Composite layer transmission matrix ieT (from - -> -)"
    botieT⁻⁻::M2
    "Composite layer source matrix ieJ (in + direction)"
    botieJ₀⁺::M2
    "Composite layer source matrix ieJ (in - direction)"
    botieJ₀⁻::M2
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

"Defined by a simple spline from Interpolations.jl"
struct LambertianSurfaceSpline{FT} <: AbstractSurfaceType
    interpolator::AbstractInterpolation{FT}
    wlGrid::AbstractArray{FT,1} # Has to be added here as it won't otherwise be available in the Lambertian Surface Routine
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
    "Lookup table type"
    luts::AbstractArray 
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
    "Reference refractive index"
    n_ref::Complex{FT}
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
mutable struct vSmartMOM_Model{PA, AE, GR, QP, TAB, TR, TAE, Ogeom, PRO}

    "Struct with all individual parameters"
    params::PA # vSmartMOM_Parameters
    
    "Truncated aerosol optics"
    aerosol_optics::AE # AbstractArray{AbstractArray{AerosolOptics}}
    "Greek coefs in Rayleigh calculations" 
    greek_rayleigh::GR # GreekCoefs
    "Quadrature points/weights, etc"
    quad_points::QP # QuadPoints

    "Array to hold cross-sections over entire atmospheric profile"
    τ_abs::TAB # AbstractArray{AbstractArray}
    "Rayleigh optical thickness"
    τ_rayl::TR # AbstractArray{AbstractArray}
    "Aerosol optical thickness"
    τ_aer::TAE # AbstractArray{AbstractArray}

    "Observational Geometry (includes sza, vza, vaz)"
    obs_geom::Ogeom # ObsGeometry
    "Atmospheric profile to use"
    profile::PRO #AtmosphericProfile
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
Base.@kwdef struct CoreScatteringOpticalProperties{FT,FT2,FT3} <:  AbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    τ::FT 
    "Single scattering albedo"
    ϖ::FT2   
    "Z scattering matrix (forward)"
    Z⁺⁺::FT3 
    "Z scattering matrix (backward)"
    Z⁻⁺::FT3
end

# Core optical Properties COP with directional cross section 
Base.@kwdef struct CoreDirectionalScatteringOpticalProperties{FT,FT2,FT3,FT4} <:  AbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    τ::FT 
    "Single scattering albedo"
    ϖ::FT2   
    "Z scattering matrix (forward)"
    Z⁺⁺::FT3 
    "Z scattering matrix (backward)"
    Z⁻⁺::FT3
    "Ross kernel; cross section projection factor along µ (G ∈ [0,1], 1 for isotropic σ)"
    G::FT4
end

Base.@kwdef struct CoreAbsorptionOpticalProperties{FT} <:  AbstractOpticalProperties
    "Absorption optical depth (scalar or wavelength dependent)"
    τ::FT 
end

# Adding Core Optical Properties, can have mixed dimensions!
function Base.:+( x::CoreScatteringOpticalProperties{xFT, xFT2, xFT3}, 
                  y::CoreScatteringOpticalProperties{yFT, yFT2, yFT3} 
                ) where {xFT, xFT2, xFT3, yFT, yFT2, yFT3} 
    # Predefine some arrays:            
    xZ⁺⁺ = x.Z⁺⁺
    xZ⁻⁺ = x.Z⁻⁺
    yZ⁺⁺ = y.Z⁺⁺
    yZ⁻⁺ = y.Z⁻⁺

    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    wy = y.τ .* y.ϖ  
    w  = wx .+ wy
    ϖ  =  w ./ τ
    
    #@show xFT, xFT2, xFT3
    all(wx .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, y.Z⁺⁺, y.Z⁻⁺)) : nothing
    all(wy .== 0.0) ? (return CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)) : nothing

    n = length(w);
    
    wy = wy ./ w
    wx = wx ./ w
    wx = reshape(wx,1,1,n)
    wy = reshape(wy,1,1,n)
        
    Z⁺⁺ = (wx .* xZ⁺⁺ .+ wy .* yZ⁺⁺) 
    Z⁻⁺ = (wx .* xZ⁻⁺ .+ wy .* yZ⁻⁺)

    CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺)  
end

# Concatenate Core Optical Properties, can have mixed dimensions!
function Base.:*( x::CoreScatteringOpticalProperties, y::CoreScatteringOpticalProperties ) 
    arr_type  = array_type(architecture(x.τ))
    x = expandOpticalProperties(x,arr_type);
    y = expandOpticalProperties(y,arr_type);
    CoreScatteringOpticalProperties([x.τ; y.τ],[x.ϖ; y.ϖ],cat(x.Z⁺⁺,y.Z⁺⁺, dims=3), cat(x.Z⁻⁺,y.Z⁻⁺, dims=3) )
end

function Base.:+( x::CoreScatteringOpticalProperties, y::CoreAbsorptionOpticalProperties ) 
    τ  = x.τ .+ y.τ
    wx = x.τ .* x.ϖ 
    #@show size(wx), size(τ)
    ϖ  = (wx) ./ τ
    CoreScatteringOpticalProperties(τ, ϖ, x.Z⁺⁺, x.Z⁻⁺)
end

function Base.:+(  y::CoreAbsorptionOpticalProperties, x::CoreScatteringOpticalProperties ) 
    x + y
end


function Base.:*( x::FT, y::CoreScatteringOpticalProperties{FT} ) where FT
    CoreScatteringOpticalProperties(y.τ * x, y.ϖ, y.Z⁺⁺, y.Z⁻⁺, y.G)
end

# From https://gist.github.com/mcabbott/80ac43cca3bee8f57809155a5240519f
function _repeat(x::AbstractArray, counts::Integer...)
    N = max(ndims(x), length(counts))
    size_y = ntuple(d -> size(x,d) * get(counts, d, 1), N)
    size_x2 = ntuple(d -> isodd(d) ? size(x, 1+d÷2) : 1, 2*N)

    ## version without mutation
    # ignores = ntuple(d -> reshape(Base.OneTo(counts[d]), ntuple(_->1, 2d-1)..., :), length(counts))
    # y = reshape(broadcast(first∘tuple, reshape(x, size_x2), ignores...), size_y)

    # ## version with mutation
    size_y2 = ntuple(d -> isodd(d) ? size(x, 1+d÷2) : get(counts, d÷2, 1), 2*N)
    y = similar(x, size_y)
    reshape(y, size_y2) .= reshape(x, size_x2)
    y
end