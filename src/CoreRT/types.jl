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
- `AbsorptionParameters`, `ScatteringParameters`, and `RTModel` hold model parameters
- `QuadPoints` holds quadrature points, weights, etc. 
- `ComputedAtmosphereProperties` and `ComputedLayerProperties` hold intermediate computed properties

=#

# Conditional type for CUDA pointer arrays (only needed when CUDA is loaded)
# When CUDA is not loaded, this will just be Nothing
# When CUDA is loaded, CUDAExt will properly handle CuArray types
const MaybeCuPtrArray = Union{AbstractArray, Nothing}

"""
    AtmosphericProfile{FT, VMR}

Vertical atmospheric state on a pressure grid, ordered from TOA to BOA.

Stores temperature, pressure (full levels and half levels), humidity,
dry and wet vertical column densities, and volume mixing ratios for
trace gases.  Constructed internally by [`model_from_parameters`](@ref)
from the raw arrays in [`vSmartMOM_Parameters`](@ref).

`p_full` has `N` levels; `p_half` has `N-1` layer-boundary pressures.
`T`, `q`, `vmr_h2o`, `vcd_dry`, `vcd_h2o`, and `Δz` are layer quantities
of length `N`.

# Fields
$(DocStringExtensions.FIELDS)
"""
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

"""
    RT_Aerosol{FT}

An aerosol species combined with its RT-relevant parameters: reference optical
depth (`τ_ref`) and the vertical pressure distribution (`profile`).
"""
mutable struct RT_Aerosol{FT<:Real}
    "Aerosol optical/microphysical properties"
    aerosol::Aerosol
    "Reference optical depth at λ_ref"
    τ_ref::FT
    "Vertical distribution as function of pressure (from Distributions.jl)"
    profile::Distribution
    "Optional analytic phase function; `nothing` selects the Mie path"
    phase_function::Union{Nothing, Scattering.AbstractAnalyticPhaseFunction}
    "Single-scattering albedo for analytic phase-function aerosols"
    ϖ::FT
end

RT_Aerosol(aerosol::Aerosol, τ_ref::FT,
           profile::Distribution) where {FT<:Real} =
    RT_Aerosol{FT}(aerosol, τ_ref, profile, nothing, one(FT))

RT_Aerosol(aerosol::Aerosol, τ_ref::FT, profile::Distribution,
           phase_function::Scattering.AbstractAnalyticPhaseFunction;
           ϖ::Real = one(FT)) where {FT<:Real} =
    RT_Aerosol{FT}(aerosol, τ_ref, profile, phase_function, convert(FT, ϖ))

"Quadrature Types for RT streams"
abstract type AbstractQuadratureType end

"""
    struct RadauQuad

Use Gauss Radau quadrature scheme, which includes the SZA as point (see Sanghavi vSmartMOM)

"""
struct RadauQuad <:AbstractQuadratureType end

"""
    struct GaussLegQuad

Half-space Gauss-Legendre quadrature on `[0, 1]` for the upper hemisphere.
The default for plane-parallel RT.
"""
struct GaussLegQuad <: AbstractQuadratureType end

"Abstract Type for Source Function Integration"
abstract type AbstractSourceType end

"Use Dummy Node Integration (DNI), SZA has to be a full node, Radau required"
struct DNI <:AbstractSourceType end

"Use Source Function Integration (SFI), Solar beam embedded in source term, can work with all quadrature schemes (faster)"
struct SFI <:AbstractSourceType end

"Abstract Type for Layer R,T and J matrices"
abstract type AbstractLayer end

"""
    CompositeLayer{FT} <: AbstractLayer

Accumulated (composite) layer matrices produced by the interaction step of the
adding/doubling RT method.  Stores the reflectance (`R`), transmission (`T`),
and source (`J`) matrices for the combined atmospheric slab from the top of
atmosphere down to the current layer.

Sign convention: `−` = outgoing (upward, decreasing τ), `+` = incoming
(downward, increasing τ).

# Fields
- `R⁻⁺::AbstractArray{FT,3}`: reflectance from incoming (+) to outgoing (−)
- `R⁺⁻::AbstractArray{FT,3}`: reflectance from outgoing (−) to incoming (+)
- `T⁺⁺::AbstractArray{FT,3}`: transmission in the incoming (+) direction
- `T⁻⁻::AbstractArray{FT,3}`: transmission in the outgoing (−) direction
- `J₀⁺::AbstractArray{FT,3}`: source in the incoming (+) direction
- `J₀⁻::AbstractArray{FT,3}`: source in the outgoing (−) direction
"""
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
    # v0.7 Phase A.2a — per-source J₀ slots. Each non-solar source carries its
    # own J₀⁺/J₀⁻ accumulator so doubling and interaction can apply per-source
    # math (e.g. thermal's expk = 1 vs solar's exp(-dτ/μ₀)). Empty NT for
    # solar-only runs → no extra memory, no behavior change.
    "Per-source composite-layer J₀ slots (v0.7 Phase A.2a, NamedTuple{(:thermal, …)})"
    J₀_by_src::NamedTuple = (;)
end

"""
    AddedLayer{FT} <: AbstractLayer

Single homogeneous layer matrices produced by the elemental/doubling steps of
the RT solver.  After doubling, these represent a full atmospheric layer that
is subsequently combined with the composite layer via the interaction step.

Lower-case field names (`r`, `t`, `j`) distinguish added-layer quantities from
composite-layer quantities (`R`, `T`, `J`).

# Fields
- `r⁻⁺::AbstractArray{FT,3}`: reflectance from incoming (+) to outgoing (−)
- `t⁺⁺::AbstractArray{FT,3}`: transmission in the incoming (+) direction
- `r⁺⁻::AbstractArray{FT,3}`: reflectance from outgoing (−) to incoming (+)
- `t⁻⁻::AbstractArray{FT,3}`: transmission in the outgoing (−) direction
- `j₀⁺::AbstractArray{FT,3}`: source in the incoming (+) direction
- `j₀⁻::AbstractArray{FT,3}`: source in the outgoing (−) direction
- `temp1`, `temp2`: pre-allocated workspace arrays to avoid allocations
- `temp1_ptr`, `temp2_ptr`: CUDA pointer arrays (ignored on CPU)
"""
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
    "Pointer to temporary space to avoid allocations (CUDA-specific, ignored on CPU)"
    temp1_ptr::MaybeCuPtrArray
    "Pointer to temporary space to avoid allocations (CUDA-specific, ignored on CPU)"
    temp2_ptr::MaybeCuPtrArray
    "Doubling workspace: geometric-progression reflectance (same shape as t⁺⁺)"
    dbl_gp_refl::Union{AbstractArray{FT,3}, Nothing} = nothing
    "Doubling workspace: source temp +"
    dbl_j₁⁺::Union{AbstractArray{FT,3}, Nothing} = nothing
    "Doubling workspace: source temp -"
    dbl_j₁⁻::Union{AbstractArray{FT,3}, Nothing} = nothing
    # v0.7 Phase A.2a — per-source j₀ slots. Each non-solar source contributes
    # its own [`SourceSlot`](@ref) (j₀⁺/j₀⁻ + per-source `expk` for doubling +
    # j₁ workspace). The doubling/interaction routines iterate this NT to
    # apply per-source math (e.g. thermal's expk = 1 vs solar's exp(-dτ/μ₀)).
    # Empty NT in solar-only runs → bit-equal to pre-A.2a behaviour.
    "Per-source j₀ slots (v0.7 Phase A.2a, NamedTuple{(:thermal, …)} of SourceSlot)"
    j₀_by_src::NamedTuple = (;)
end

"""
    SourceSlot{FT} <: AbstractLayer

Per-source carrier for the j₀ doubling state. One slot per non-solar source
type in the active `prepared_sources` (Phase A.2a). Holds the source-
specific elemental j₀⁺/j₀⁻ (filled by `contribute!`), the doubling-
workspace j₁ buffers, and the per-source `expk` attenuation factor used by
the doubling recurrence (`exp(-dτ/μ₀)` for solar, `ones` for thermal, etc.).

The composite-layer counterpart is the `J₀_by_src` NamedTuple on
[`CompositeLayer`](@ref); each slot key (`:thermal`, …) maps to a plain
`(J₀⁺, J₀⁻)` pair (the composite layer has no doubling workspace, since
doubling already produced the full per-source J₀ before interaction).

# Fields
- `j₀⁺::AbstractArray{FT,3}` — incoming-direction source for this source
- `j₀⁻::AbstractArray{FT,3}` — outgoing-direction source
- `dbl_j₁⁺::AbstractArray{FT,3}` — doubling-loop workspace
- `dbl_j₁⁻::AbstractArray{FT,3}` — doubling-loop workspace
- `expk::AbstractArray{FT,1}` — per-spectral attenuation factor squared each
  doubling step; for thermal this stays `1.0` because thermal is self-
  generated (no incident beam to attenuate). For solar this is the
  legacy `exp(-dτ/μ₀)` — but the solar slot lives in the legacy `j₀⁺/j₀⁻`
  fields, not in `j₀_by_src`.
"""
Base.@kwdef struct SourceSlot{FT} <: AbstractLayer
    j₀⁺::AbstractArray{FT,3}
    j₀⁻::AbstractArray{FT,3}
    dbl_j₁⁺::AbstractArray{FT,3}
    dbl_j₁⁻::AbstractArray{FT,3}
    expk::AbstractArray{FT,1}
end

"""
    CompositeSourceSlot{FT}

Per-source composite-layer carrier for the J₀ vectors after doubling +
interaction. Holds J₀⁺ and J₀⁻ on the active architecture; no workspace
buffers are needed (the interaction step does the combination in place).
"""
Base.@kwdef struct CompositeSourceSlot{FT} <: AbstractLayer
    J₀⁺::AbstractArray{FT,3}
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
    "Added layer source matrix j (in + direction)"
    j₀⁺::AbstractArray{FT,3}
    "Added layer source matrix j (in - direction)"
    j₀⁻::AbstractArray{FT,3}

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

"""
    LambertianSurfaceScalar{FT} <: AbstractSurfaceType

Isotropic Lambertian surface with a single scalar albedo `α` shared across
all spectral points in a band.  The BRDF is `ρ = α/π`, independent of
viewing and illumination angles.  Only the `m = 0` Fourier moment is nonzero.

Supports analytical Jacobians with respect to albedo in linearized RT.

# Fields
- `albedo::FT`: hemispherical albedo ∈ [0, 1].
"""
struct LambertianSurfaceScalar{FT} <: AbstractSurfaceType
    "Albedo (scalar)"
    albedo::FT
end

"""
    LambertianSurfaceSpectrum{FT} <: AbstractSurfaceType

Isotropic Lambertian surface with a per-spectral-point albedo vector.
The vector length must match the number of spectral points in the band.

# Fields
- `albedo::AbstractArray{FT,1}`: spectral albedo vector, one value per grid point.
"""
struct LambertianSurfaceSpectrum{FT} <: AbstractSurfaceType
    "Albedo (vector)"
    albedo::AbstractArray{FT,1}
end

"""
    rpvSurfaceScalar{FT} <: AbstractSurfaceType

Rahman-Pinty-Verstraete (1993) four-parameter empirical BRDF model.
Combines a Minnaert angular shape, a Henyey-Greenstein-like hot-spot
function, and a geometric bowl/bell correction.  Scalar only (Stokes-I);
no built-in linearization.

# Fields
$(DocStringExtensions.FIELDS)
"""
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

"""
    RossLiSurfaceScalar{FT} <: AbstractSurfaceType

Ross-Li kernel-based BRDF (Lucht, Schaaf & Strahler, 2000).  Decomposes
the surface reflectance into three semi-physical kernels: isotropic,
RossThick (volumetric canopy scattering), and LiSparse (geometric
mutual shadowing).  Standard MODIS/RAMI parameterization.  Scalar only
(Stokes-I); no built-in linearization.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RossLiSurfaceScalar{FT} <: AbstractSurfaceType
    "Volumetric RossThick  fraction"
    fvol::FT
    "Geometric LiSparse fraction"
    fgeo::FT
    "Isotropic reflectance fraction"
    fiso::FT
end

"""
    CoxMunkSurface{FT} <: AbstractSurfaceType

Cox-Munk (1954) ocean surface: Fresnel reflection from wind-roughened wave facets.
Supports full polarization (Mueller matrix) and optional Lambertian whitecap contribution.

Isotropic slope distribution: σ² = 0.003 + 0.00512·U.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct CoxMunkSurface{FT} <: AbstractSurfaceType
    "10-m wind speed [m/s]"
    wind_speed::FT
    "Complex refractive index of water; nothing → built-in Segelstein 1981 lookup"
    n_water::Union{Nothing, Complex{FT}, Vector{Complex{FT}}} = nothing
    "Whitecap Lambertian albedo (Koepke 1984)"
    whitecap_albedo::FT = 0.22
    "Include whitecap contribution"
    include_whitecaps::Bool = true
    "Include Smith (1967) shadow masking"
    shadowing::Bool = true
end

"Defined by Legendre polynomial terms as function of spectral grid, which is scaled to `[-1,1]` (degree derived from length of `a_coeff`)"
struct LambertianSurfaceLegendre{FT} <: AbstractSurfaceType
    "albedo = legendre_coeff[1] * P₀ + legendre_coeff[2]*P₁ + legendre_coeff[3]*P₂ + ... "
    legendre_coeff::AbstractArray{FT,1}
end

"Defined by a simple spline from Interpolations.jl"
struct LambertianSurfaceSpline{FT} <: AbstractSurfaceType
    interpolator::AbstractInterpolation{FT}
    wlGrid::AbstractArray{FT,1}
end

"""
    CanopySurface{FT} <: AbstractSurfaceType

Composite lower boundary: a vegetation canopy (one or more scattering layers)
backed by a soil BRDF. When used as the surface in `rt_run`, the canopy
sub-layers are processed internally via the adding-doubling method before
interacting with the soil, producing an effective canopy+soil reflectance.

Requires `CanopyOptics.jl` for leaf-angle distribution and scattering models.
"""
mutable struct CanopySurface{FT} <: AbstractSurfaceType
    "Soil BRDF (any AbstractSurfaceType, e.g. LambertianSurfaceScalar)"
    soil::AbstractSurfaceType
    "Total leaf area index"
    LAI::FT
    "Number of canopy sub-layers (1 = big-leaf)"
    n_layers::Int
    "Leaf angle distribution (from CanopyOptics, e.g. spherical_leaves())"
    LAD
    "Canopy scattering model (e.g. BiLambertianCanopyScattering)"
    canopy_scattering
    "CanopyOptics Z-matrix integration controls"
    canopy_quadrature
    "Canopy clumping model; affects effective G for propagation, not Z normalization"
    canopy_clumping
    "Leaf reflectance (scalar or spectral vector on leaf_optics_grid)"
    leaf_reflectance::Union{FT, Vector{FT}}
    "Leaf transmittance (scalar or spectral vector on leaf_optics_grid)"
    leaf_transmittance::Union{FT, Vector{FT}}
    "Wavelength/wavenumber grid for spectral leaf R/T (nothing if scalar)"
    leaf_optics_grid::Union{Nothing, Vector{FT}}
    "Unit of leaf_optics_grid: :nm (wavelength) or :cm_inv (wavenumber)"
    grid_unit::Symbol
    "Include within-canopy atmospheric absorption between sub-layers"
    include_atm::Bool
    "Pressure thickness of within-canopy air column in hPa (e.g. 3 hPa ≈ 30 m canopy); nothing = no canopy atmosphere"
    canopy_dp::Union{Nothing, FT}
    "Per-sub-layer LAI fractions (nothing = uniform); must sum to 1"
    lai_fractions::Union{Nothing, Vector{FT}}
    "Within-canopy atmospheric τ (spectral, set by rt_run for include_atm=true)"
    _within_canopy_τ::Union{Nothing, Vector{FT}}
    "Lazily-initialized cache (Zup, Zdown, working arrays); set on first use"
    _cache  # Union{Nothing, CanopyCache} — typed at first init
end

function CanopySurface(; soil::AbstractSurfaceType,
                        LAI, n_layers::Int=1, LAD=nothing,
                        canopy_scattering=nothing,
                        canopy_quadrature=nothing,
                        canopy_clumping=nothing,
                        leaf_reflectance=0.4, leaf_transmittance=0.05,
                        leaf_optics_grid=nothing, grid_unit::Symbol=:nm,
                        include_atm::Bool=false,
                        canopy_dp=nothing,
                        lai_fractions=nothing)
    FT = typeof(float(LAI))
    if LAD === nothing
        LAD = CanopyOptics.spherical_leaves()
    end
    if canopy_scattering === nothing
        R_scalar = leaf_reflectance isa Number ? FT(leaf_reflectance) :
                   FT(sum(leaf_reflectance) / length(leaf_reflectance))
        T_scalar = leaf_transmittance isa Number ? FT(leaf_transmittance) :
                   FT(sum(leaf_transmittance) / length(leaf_transmittance))
        canopy_scattering = CanopyOptics.BiLambertianCanopyScattering(
            R=R_scalar, T=T_scalar)
    end
    if canopy_quadrature === nothing
        canopy_quadrature = CanopyOptics.CanopyQuadrature()
    end
    if canopy_clumping === nothing
        canopy_clumping = CanopyOptics.NoClumping{FT}()
    end
    lr = leaf_reflectance isa Number ? FT(leaf_reflectance) : convert(Vector{FT}, leaf_reflectance)
    lt = leaf_transmittance isa Number ? FT(leaf_transmittance) : convert(Vector{FT}, leaf_transmittance)

    # Accept Unitful grid (e.g. collect(400:2500) .* u"nm") or plain floats + grid_unit
    if leaf_optics_grid !== nothing && eltype(leaf_optics_grid) <: Unitful.Quantity
        u_grid = unit(first(leaf_optics_grid))
        if u_grid == u"nm"
            grid_unit = :nm
        elseif u_grid == u"cm^-1"
            grid_unit = :cm_inv
        else
            error("Unsupported leaf_optics_grid unit: $u_grid; expected u\"nm\" or u\"cm^-1\"")
        end
        lg = convert(Vector{FT}, ustrip.(leaf_optics_grid))
    else
        lg = leaf_optics_grid === nothing ? nothing : convert(Vector{FT}, leaf_optics_grid)
    end
    @assert grid_unit in (:nm, :cm_inv) "grid_unit must be :nm or :cm_inv"

    dp = canopy_dp === nothing ? nothing : FT(canopy_dp)
    lf = lai_fractions === nothing ? nothing : convert(Vector{FT}, lai_fractions)
    CanopySurface{FT}(soil, FT(LAI), n_layers, LAD, canopy_scattering,
                      canopy_quadrature, canopy_clumping,
                      lr, lt, lg, grid_unit, include_atm, dp, lf, nothing, nothing)
end

"""
    struct AbsorptionParameters

A struct which holds all absorption-related parameters (before any computations).
The species list per band is the union of `fixed_molecules` (no Jacobian) and
`variable_molecules` (Jacobian computed). H₂O is handled separately: when
`atmospheric_profile.q` is provided, line absorption is computed from
`profile.vmr_h2o` using `h2o_lut` (if present) or HITRAN-on-the-fly otherwise,
and is automatically treated as variable in the linearised path.
"""
mutable struct AbsorptionParameters{FM,VM,V,BF,CE,LT,HL}
    "Fixed-abundance molecules per band (no Jacobian computed). Must NOT include H2O."
    fixed_molecules::FM
    "Variable-abundance molecules per band (Jacobian computed w.r.t. VMR). Must NOT include H2O."
    variable_molecules::VM
    "Volume-Mixing Ratios"
    vmr::V
    "Type of broadening function (Doppler/Lorentz/Voigt)"
    broadening_function::BF
    "Complex Error Function to use in Voigt calculations"
    CEF::CE
    "Wing cutoff to use in cross-section calculation (cm⁻¹)"
    wing_cutoff::Real
    "Lookup table per (band, position-in-fixed_molecules ∪ variable_molecules)"
    luts::LT
    "Optional H₂O LUT (one per band). `nothing` means HITRAN-on-the-fly fallback."
    h2o_lut::HL
    "Optional list of HITRAN CIA file paths (one per collision pair)"
    cia_files::Vector{String}
    "Optional path to AER MT_CKD water-vapor continuum NetCDF (e.g. absco-ref_wv-mt-ckd.nc)"
    mtckd_file::String
end

"""
    struct ScatteringParameters

A struct which holds all scattering-related parameters (before any computations)
"""
mutable struct ScatteringParameters{FT<:Real}
    "List of scattering aerosols and their properties"
    rt_aerosols::Vector{RT_Aerosol{FT}}
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
    RTNumericalParameters{FT}

Centralised home for tunable numerical knobs that were previously hardcoded
in the RT kernels. Lives once on `vSmartMOM_Parameters` / `RTModel`, threaded
through to the kernels where needed. Designed to grow as more constants are
lifted (Glint Legendre expansion cap, etc.).

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct RTNumericalParameters{FT<:AbstractFloat}
    "Doubling resolution: cap on the elemental-layer optical depth
    `dτ_max = threshold · μ_min`, where `μ_min` is the smallest TRUE
    (positive-weight) quadrature stream cosine. Smaller values produce
    finer initial layers (more doublings → more accurate single-scatter
    at the cost of more roundoff through doubling). Default `0.001` is
    conservative legacy; raising to `~0.01` reduces doubling iterations."
    dτ_max_threshold::FT = FT(0.001)

    "Absolute floor on `dτ_max` (and hence on `dτ_initial`), expressed
    as a multiple of `eps(FT)`. Prevents the doubling discretisation
    from collapsing below FT precision regardless of geometry. Default
    `1024·eps(FT)` ≈ 1.2e-4 for Float32 / 2.3e-13 for Float64 — the
    F32 value caps `ndoubl` around 13–14 for τ_total≤1 (well within
    Float32 representability); the F64 value is far below any sensible
    `threshold·μ_min` so it never activates in practice."
    dτ_min_floor::FT = FT(1024) * eps(FT)

    "BLAS thread cap applied at every `rt_run` invocation. `nothing`
    means leave the BLAS thread setting alone (use whatever
    `LinearAlgebra.BLAS.get_num_threads()` returned at session start).
    Why cap: the batched-GEMM call sites in `cpu_batched.jl` and the
    elemental/doubling/interaction kernels operate on small matrices
    (NSTREAMS·n_stokes ≈ 12-48 per side) tiled across a wide spectral
    batch. Multi-threaded BLAS coordination dominates the work at that
    matrix size and serializes against the spectral-batch parallelism.
    The VLIDORT baseline harness empirically picked 8 as the sweet spot
    on a 128-thread machine; `1` is the most defensive choice when a
    single process owns all spectral points and threading happens at
    the batched outer level. Set via the `numerics.blas_threads` YAML
    key, or programmatically with
    `RTNumericalParameters{FT}(blas_threads = 8)`."
    blas_threads::Union{Nothing, Int} = nothing

    "Verbose output flag. When `false` (default) `rt_run`/`rt_run_lin`
    suppress the per-call `print_timer()` dump at the end of the run.
    Set `true` (YAML key `numerics.verbose: true` or programmatically
    `RTNumericalParameters{FT}(verbose = true)`) to bring back the
    timing tree — useful for profiling, noisy for production loops.
    Per-band `Computing profile for X...` and `AOD at band X...`
    messages were demoted to `@debug` and are silent unless you set
    `ENV[\"JULIA_DEBUG\"] = \"vSmartMOM\"`."
    verbose::Bool = false
end

"""
    vSmartMOM_Parameters{FT<:Real}

Top-level container for all user-specified model parameters **before** any
derived quantities are computed.  Groups radiative-transfer settings
(spectral bands, BRDF, quadrature, polarization), observation geometry
(SZA, VZA, VAZ), the atmospheric profile (T, p, q), and optional
absorption/scattering sub-parameter structs.

# Fields
$(DocStringExtensions.FIELDS)
"""
mutable struct vSmartMOM_Parameters{FT<:Real} 

    # radiative_transfer group
    "Spectral bands — Vector of wavenumber grids, one per band"
    spec_bands::Vector{Vector{FT}}
    "Surface (Bidirectional Reflectance Distribution Function)"
    brdf::Vector{<:AbstractSurfaceType}
    "Quadrature type for RT streams (RadauQuad/GaussLegQuad)"
    quadrature_type::AbstractQuadratureType
    "Type of polarization (I/IQ/IQU/IQUV)"
    polarization_type::AbstractPolarizationType
    "Hard cutoff for maximum number of Fourier moments to loop over"
    max_m::Integer
    "Exclusion angle for forward peak in degrees (legacy — see `truncation`)"
    Δ_angle::FT
    "Truncation length for legendre terms (scalar for now, can do `nBand` later)"
    l_trunc::Integer
    "Phase-function truncation method.
    Defaults to `δBGE(l_trunc, Δ_angle)` for backward compatibility;
    set explicitly to `NoTruncation()` for canopy-only or
    smooth-phase-function runs (Sanghavi & Stephens 2015 §2 — the
    `f_tr → 0` limit). The `Δ_angle` field above is then ignored."
    truncation::Scattering.AbstractTruncationType
    "Depolarization factor"
    depol::FT
    "Numerical-knob parameters (doubling threshold, etc.) — defaults to
    `RTNumericalParameters{FT}()` if not specified in YAML."
    numerics::RTNumericalParameters
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

    # === Phase D — `nstreams`-first schema (v0.7) ============================
    # New-schema fields appended after the existing positional fields so the
    # legacy positional constructor call in `parameters_from_dict` and any
    # downstream consumers using positional construction keep working unchanged.
    # All default to `nothing` — populated by the parser only when the user
    # opts into the new schema (presence of `nstreams` in YAML).

    "Phase D — primary user-facing resolution knob: weighted streams per
    hemisphere. Public contract: `stream_l_cap = 2·nstreams - 1`. When
    `nothing`, the parser populated `max_m`/`l_trunc` from a legacy YAML;
    new-schema configs set this and may omit `max_m`/`l_trunc`."
    nstreams::Union{Int, Nothing}

    "Phase D — explicit per-band Fourier loop bound (order, optional). When
    set, clamps the trait-aggregator output. Lower hard cap on top of
    `stream_l_cap`."
    m_max_override::Union{Int, Nothing}

    "Phase D — derived projection cap: `2·nstreams - 1` for new schema, or
    `l_trunc` for legacy. Internally consumed by the trait aggregator and
    the truncation-resolver. Always populated; defaults to `l_trunc` until
    `nstreams` is explicitly set."
    stream_l_cap::Int

    "Phase D — legacy `l_trunc` value retained verbatim from YAML when the
    user explicitly set it. `nothing` for new-schema configs that derive
    `stream_l_cap` from `nstreams` instead."
    legacy_l_cap_override::Union{Int, Nothing}

end

"""
    QuadPoints{FT}

Gauss or Gauss-Radau quadrature points and weights used for the angular
discretisation of the radiative transfer equation.  Also stores the cosine
of the solar zenith angle (`μ₀`) and its index within the quadrature grid.

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
    "Total number of quadrature points (weighted streams + zero-weight SZA/VZA output nodes)"
    Nquad::Int
    "Number of weighted streams per hemisphere (count of nonzero weights). Public contract: stream_l_cap = 2·Nstreams - 1."
    Nstreams::Int
end

# ============================================================================
# New hierarchical model types (Oceananigans-inspired)
# ============================================================================

"Abstract base type for all RT models"
abstract type AbstractRTModel{ARCH<:AbstractArchitecture, FT<:AbstractFloat} end

"""
    SolverConfig{FT, PT, QT}

Immutable RT solver configuration: polarization mode, quadrature scheme,
Fourier truncation, and Legendre truncation. Fixed after model construction;
never differentiated.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SolverConfig{FT<:AbstractFloat, PT<:AbstractPolarizationType, QT<:AbstractQuadratureType}
    "Type of polarization (Stokes_I / IQU / IQUV)"
    polarization_type::PT
    "Quadrature type (RadauQuad / GaussLegQuad)"
    quadrature_type::QT
    "Per-band Fourier loop bound, **order semantics**: loop runs `m = 0:m_max_bands[iBand]`. Equals `n_fourier_moments_bands .- 1`."
    m_max_bands::Vector{Int}
    "Per-band Fourier moment count: `n_fourier_moments_bands[iBand] = m_max_bands[iBand] + 1`. Provided so consumers don't have to remember the count↔order convention."
    n_fourier_moments_bands::Vector{Int}
    "Per-band max truncated Legendre index"
    l_max::Vector{Int}
    "Legendre truncation order (user-specified)"
    l_trunc::Int
    "Exclusion angle for forward peak [deg]"
    Δ_angle::FT
    "Depolarization factor"
    depol::FT
    "Phase C flag: when `true` (current default), per-band `m_max_bands` are derived from `component_m_max(c, ctx)` traits across active components — Cox-Munk / RPV / RossLi / canopy run to their full `user_l_cap` instead of the historical half-truncated aggregator. Flip to `false` to fall back to the legacy `min(ceil((l_max+1)/2), params.max_m)` aggregator (bit-equal to Phase B)."
    use_component_traits::Bool
end

"""
    Atmosphere{FT, VMR}

The atmospheric column: profile data plus spectral grid definitions.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Atmosphere{FT<:AbstractFloat, VMR}
    "Atmospheric profile (T, p, q, vmr, vcd, etc.)"
    profile::AtmosphericProfile{FT, VMR}
    "Spectral bands — Vector of wavenumber grids, one per band"
    spec_bands::Vector{Vector{FT}}
end

"""
    RayleighScattering{FT, GC}

Precomputed Rayleigh and Cabannes (elastic) scattering properties.
Derived from the depolarization ratio; fixed for a given spectral band.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RayleighScattering{FT<:AbstractFloat, GC<:Scattering.GreekCoefs}
    "Greek coefficients for total Rayleigh scattering (single set or per-band)"
    greek_rayleigh::Union{GC, Vector{GC}}
    "Greek coefficients for Cabannes (pure elastic) Rayleigh per band"
    greek_cabannes::Vector{GC}
    "Pure elastic (Cabannes) fraction of Rayleigh scattering per band"
    ϖ_Cabannes::Vector{FT}
end

"""
    AerosolState{FT, AO}

Per-band aerosol scattering optics and optical depth profiles.
Primary differentiable state for aerosol retrievals.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AerosolState{FT<:AbstractFloat, AO, AT<:AbstractArray{FT}}
    "Truncated aerosol optics: aerosol_optics[iBand][iAer]"
    aerosol_optics::Vector{Vector{AO}}
    "Aerosol optical depth profiles: `τ_aer[iBand][iAer, iLayer]` (or `[iAer, nSpec, iLayer]` for lin)"
    τ_aer::Vector{AT}
end

"""
    Optics{FT, RS, AS}

Container for all precomputed optical properties that feed into the RT solver.
Groups Rayleigh scattering, aerosol scattering, absorption optical depths,
and Rayleigh optical depths.

**AD boundary**: `aerosols` and `τ_abs` hold differentiable state;
`rayleigh` and `τ_rayl` are typically fixed.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Optics{FT<:AbstractFloat, RS<:RayleighScattering, AS<:AerosolState}
    "Rayleigh/Cabannes scattering properties"
    rayleigh::RS
    "Aerosol scattering optics and optical depths"
    aerosols::AS
    "Absorption optical depth: τ_abs[iBand][nSpec × nLayers]"
    τ_abs::Vector{Matrix{FT}}
    "Rayleigh optical depth: τ_rayl[iBand][nSpec × nLayers]"
    τ_rayl::Vector{Matrix{FT}}
end

"""
    RTModel{ARCH, FT} <: AbstractRTModel{ARCH, FT}

Central model state for vSmartMOM radiative transfer.

Only two type parameters:
- `ARCH`: compute architecture (`CPU` or `GPU`)
- `FT`: floating-point precision (`Float32` or `Float64`)

All physics sub-components are organized hierarchically:
- `solver`: RT solver configuration (polarization, quadrature, truncation)
- `geometry`: observation geometry (SZA, VZA, VAZ, observer altitude)
- `quad_points`: precomputed quadrature points and weights
- `atmosphere`: atmospheric state (profile + spectral bands)
- `optics`: all precomputed optical properties (Rayleigh, aerosol, absorption)
- `surfaces`: per-band surface BRDF models

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RTModel{ARCH<:AbstractArchitecture, FT<:AbstractFloat} <: AbstractRTModel{ARCH, FT}
    "Compute architecture (CPU/GPU)"
    architecture::ARCH
    "RT solver configuration"
    solver::SolverConfig{FT}
    "Numerical-knob parameters (doubling threshold, etc.) — see
    `RTNumericalParameters`"
    numerics::RTNumericalParameters{FT}
    "Observation geometry"
    geometry::ObsGeometry{FT}
    "Quadrature points and weights"
    quad_points::QuadPoints{FT}
    "Atmospheric state (profile + spectral bands)"
    atmosphere::Atmosphere{FT}
    "Precomputed optical properties"
    optics::Optics{FT}
    "Surface models, one per spectral band"
    surfaces::Vector{<:AbstractSurfaceType}
    "Source-term scene (v0.6 source-term refactor). Defaults to a single
    [`SolarBeam`](@ref) — i.e. the historical unit-Stokes-I direct solar
    illumination — and can be overridden at `rt_run` call time via the
    `sources=` kwarg. Mirrors the abstract-typed `surfaces` slot."
    sources::AbstractSource
end

# ── Accessor functions (Oceananigans-style) ───────────────────────────────

"Return the compute architecture of the model"
Architectures.architecture(m::AbstractRTModel) = m.architecture
"Return the float type of the model"
float_type(::AbstractRTModel{ARCH, FT}) where {ARCH, FT} = FT
"Return the array constructor for the model's architecture"
Architectures.array_type(m::AbstractRTModel) = Architectures.array_type(m.architecture)
"Return the polarization type"
polarization_type(m::AbstractRTModel) = m.solver.polarization_type
"Number of aerosol species"
n_aerosols(m::RTModel) = isempty(m.optics.aerosols.aerosol_optics) ? 0 :
    length(first(m.optics.aerosols.aerosol_optics))
"""
    m_max_bands(m::RTModel) -> Vector{Int}

Per-band Fourier loop bound (order). Loop runs `m = 0:m_max_bands(m)[iBand]`.
"""
m_max_bands(m::RTModel) = m.solver.m_max_bands

"""
    n_fourier_moments_bands(m::RTModel) -> Vector{Int}

Per-band Fourier moment count (= `m_max_bands(m) .+ 1`). Provided for
consumers that explicitly want count rather than order.
"""
n_fourier_moments_bands(m::RTModel) = m.solver.n_fourier_moments_bands

# ── RTModel accessors ─────────────────────────────────────────────────────
get_surface(m::RTModel, iBand) = m.surfaces[iBand]
get_surfaces(m::RTModel) = m.surfaces
get_spec_bands(m::RTModel) = m.atmosphere.spec_bands

# ── Convenience property forwarding ──────────────────────────────────────
# Allows shorthand access (model.τ_abs, model.profile, model.obs_geom, etc.)
# on RTModel, forwarding to the appropriate sub-struct field.

function Base.getproperty(m::RTModel, s::Symbol)
    s === :τ_abs && return m.optics.τ_abs
    s === :τ_rayl && return m.optics.τ_rayl
    s === :τ_aer && return m.optics.aerosols.τ_aer
    s === :aerosol_optics && return m.optics.aerosols.aerosol_optics
    s === :greek_rayleigh && return m.optics.rayleigh.greek_rayleigh
    s === :greek_cabannes && return m.optics.rayleigh.greek_cabannes
    s === :ϖ_Cabannes && return m.optics.rayleigh.ϖ_Cabannes
    s === :obs_geom && return m.geometry
    s === :profile && return m.atmosphere.profile
    s === :l_max && return m.solver.l_max
    return getfield(m, s)
end

# ── Linearized optics ────────────────────────────────────────────────────

"""
    OpticsLin{FT}

Linearized (Jacobian) counterpart of [`Optics`](@ref). Each field stores
derivatives of the corresponding forward-model field with respect to the
physical state vector.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct OpticsLin{A, B, C}
    "∂τ_abs/∂x per band: Vector of arrays [NGas × nSpec × nLayers]"
    τ̇_abs::A
    "∂τ_aer/∂x per band: Vector of arrays [NAer × 7 × nSpec × nLayers]"
    τ̇_aer::B
    "Linearized aerosol optics per band per aerosol"
    lin_aerosol_optics::C
end

# =========================================================================

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

"""
    CoreScatteringOpticalProperties{FT,FT2,FT3} <: AbstractOpticalProperties

Core optical properties for a single atmospheric layer used by the RT solver.
Bundles the optical depth, single-scattering albedo, and the forward/backward
scattering phase-matrix expansion (`Z` matrices).

Multiple `CoreScatteringOpticalProperties` can be combined with `+` (mixing
Rayleigh and aerosol) or `*` (vertical concatenation).

# Fields
- `τ::FT`: layer optical depth (scalar or per-wavelength vector)
- `ϖ::FT2`: single-scattering albedo
- `Z⁺⁺::FT3`: forward-scattering Z matrix
- `Z⁻⁺::FT3`: backward-scattering Z matrix
"""
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
    ϖ  = w ./ ifelse.(τ .> zero.(τ), τ, one.(τ))
    
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
    ϖ  = wx ./ ifelse.(τ .> zero.(τ), τ, one.(τ))
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

    #     ## version with mutation
    size_y2 = ntuple(d -> isodd(d) ? size(x, 1+d÷2) : get(counts, d÷2, 1), 2*N)
    y = similar(x, size_y)
    reshape(y, size_y2) .= reshape(x, size_x2)
    y
end

# ===========================================================
# GPU Memory Pre-allocation Workspace (added for unified branch)
# ===========================================================

"""
    RTWorkspace{FT, AT3}

Pre-allocated workspace for RT computations to avoid repeated GPU memory allocations.
Used in the doubling and interaction kernels.
"""
mutable struct RTWorkspace{FT, AT3<:AbstractArray{FT,3}}
    "Geometric progression: (I - R⁺⁻R⁻⁺)⁻¹"
    gp_refl::AT3
    "T⁺⁺ × gp_refl"
    tt_gp::AT3
    "Temporary source J₁⁺"
    J₁⁺::AT3
    "Temporary source J₁⁻"
    J₁⁻::AT3
    # Temporaries for batch_inv! (pivot and info arrays)
    "Pivot array for CUBLAS LU factorization"
    pivot::AbstractMatrix{Cint}
    "Info array for CUBLAS LU factorization"
    info::AbstractVector{Cint}
    # Temporaries for interaction (3D: NquadN × NquadN × nSpec)
    "Interaction temporary: (I - R⁺⁻r⁻⁺)⁻¹"
    tmp_inv::AT3
    "Interaction temporary: T × tmp_inv"
    T_inv::AT3
    "General purpose 3D temporary"
    tmp3d_a::AT3
    "General purpose 3D temporary"
    tmp3d_b::AT3
end

"""
    make_rt_workspace(FT, arr_type, NquadN, nSpec)

Create an RTWorkspace with CPU Array-backed temporaries.
"""
function make_rt_workspace(FT::Type, arr_type, NquadN::Int, nSpec::Int)
    dims3 = (NquadN, NquadN, nSpec)
    dims_J = (NquadN, 1, nSpec)
    
    RTWorkspace(
        arr_type(zeros(FT, dims3)),     # gp_refl
        arr_type(zeros(FT, dims3)),     # tt_gp
        arr_type(zeros(FT, dims_J)),    # J₁⁺
        arr_type(zeros(FT, dims_J)),    # J₁⁻
        zeros(Cint, NquadN, nSpec),     # pivot (replaced by GPU version if needed)
        zeros(Cint, nSpec),             # info  (replaced by GPU version if needed)
        arr_type(zeros(FT, dims3)),     # tmp_inv
        arr_type(zeros(FT, dims3)),     # T_inv
        arr_type(zeros(FT, dims3)),     # tmp3d_a
        arr_type(zeros(FT, dims3)),     # tmp3d_b
    )
end

"""
    make_gpu_rt_workspace(FT, NquadN, nSpec)

Create an RTWorkspace with GPU-backed temporaries. 
Stub function -- overwritten by the CUDA extension when CUDA is loaded.
"""
function make_gpu_rt_workspace end
