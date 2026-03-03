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

"Defined by Legendre polynomial terms as function of spectral grid, which is scaled to [-1,1] (degree derived from length of `a_coeff`)"
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
                      lr, lt, lg, grid_unit, include_atm, dp, lf, nothing, nothing)
end

"""
    struct AbsorptionParameters

A struct which holds all absorption-related parameters (before any computations).
`molecules` is the full list per band; `fixed_molecules` and `variable_molecules`
partition it for analytic Jacobians (fixed = no derivative, variable = derivative computed).
"""
mutable struct AbsorptionParameters{M,FM,VM,V,BF,CE,LT}
    "Molecules to use for absorption calculations (`nBand, nMolecules`)"
    molecules::M
    "Fixed-abundance molecules per band (no Jacobian computed)"
    fixed_molecules::FM
    "Variable-abundance molecules per band (Jacobian computed w.r.t. VMR)"
    variable_molecules::VM
    "Volume-Mixing Ratios"
    vmr::V
    "Type of broadening function (Doppler/Lorentz/Voigt)"
    broadening_function::BF
    "Complex Error Function to use in Voigt calculations"
    CEF::CE
    "Wing cutoff to use in cross-section calculation (cm⁻¹)"
    wing_cutoff::Real
    "Lookup table type"
    luts::LT
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
    "Number of quadrature points"
    Nquad::Int
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
    "Quadrature type (RadauQuad / GaussQuadHemisphere / GaussQuadFullSphere)"
    quadrature_type::QT
    "Hard cutoff for maximum number of Fourier moments (scalar, user-specified)"
    max_m::Int
    "Per-band adjusted max Fourier iterations (derived from aerosol greek coef lengths)"
    max_m_bands::Vector{Int}
    "Per-band max truncated Legendre index"
    l_max::Vector{Int}
    "Legendre truncation order (user-specified)"
    l_trunc::Int
    "Exclusion angle for forward peak [deg]"
    Δ_angle::FT
    "Depolarization factor"
    depol::FT
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
    "Aerosol optical depth profiles: τ_aer[iBand][iAer, iLayer] (or [iAer, nSpec, iLayer] for lin)"
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
"Per-band max Fourier moments"
max_m_bands(m::RTModel, iBand::Int) = m.solver.max_m_bands[iBand]

# ── RTModel accessors ─────────────────────────────────────────────────────
get_surface(m::RTModel, iBand) = m.surfaces[iBand]
get_surfaces(m::RTModel) = m.surfaces
get_spec_bands(m::RTModel) = m.atmosphere.spec_bands
get_max_m(m::RTModel) = m.solver.max_m

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
    s === :max_m && return m.solver.max_m_bands
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