"""
    Abstract type hierarchy for the next-generation aerosol pipeline.

These types describe per-bin physics (size distribution, refractive index,
mixing rule, bin integration strategy) in a way that the Mie kernel
replacement (deferred to a follow-up PR) can dispatch on without churning
IO or scene-loading code. The existing config-driven types in `types.jl`
are unchanged and still drive the legacy paths.
"""

# ============================================================================
# Size distributions
# ============================================================================

abstract type AbstractSizeDistribution{FT} end

"""
    LogNormalSizeDistribution{FT} <: AbstractSizeDistribution{FT}

Lognormal number-size distribution parameterized by geometric mean radius
`rg` (μm) and geometric standard deviation `σg` (dimensionless).
"""
struct LogNormalSizeDistribution{FT} <: AbstractSizeDistribution{FT}
    rg::FT
    σg::FT
end

"""
    GammaSizeDistribution{FT} <: AbstractSizeDistribution{FT}

Gamma size distribution with shape `α` and effective radius `β` (μm).
"""
struct GammaSizeDistribution{FT} <: AbstractSizeDistribution{FT}
    α::FT
    β::FT
end

"""
    DiscreteBinSizeDistribution{FT} <: AbstractSizeDistribution{FT}

Tabulated size spectrum used by sectional schemes (e.g., TOMAS).
`r_centers` and `r_widths` are bin-center radii and dlogD widths (μm);
`dN_dlogD` is the spectrum.
"""
struct DiscreteBinSizeDistribution{FT} <: AbstractSizeDistribution{FT}
    r_centers::Vector{FT}
    r_widths::Vector{FT}
    dN_dlogD::Vector{FT}
end

# ============================================================================
# Bin integration strategies (how a sectional scheme is consumed by Mie)
# ============================================================================

abstract type AbstractBinIntegration end

"""Fit a single lognormal to all bins, then run Mie once."""
struct LogNormalFit <: AbstractBinIntegration end

"""
    ConstantIntegrationPerBin(nquad=160)

Integrate Mie quantities over each sectional bin with `nquad`
Gauss-Legendre points in log-size space, assuming constant `dN/dlogD`
inside each bin. The integrated number in each bin is conserved exactly.
"""
struct ConstantIntegrationPerBin <: AbstractBinIntegration
    nquad::Int
    function ConstantIntegrationPerBin(nquad::Integer = 160)
        nquad > 0 || throw(ArgumentError("nquad must be positive"))
        return new(Int(nquad))
    end
end

"""
    LinearIntegrationPerBin(nquad=160)

Integrate Mie quantities over each sectional bin with `nquad`
Gauss-Legendre points in log-size space, using a conservative piecewise
linear reconstruction of `dN/dlogD` from neighboring bins. The bin-integrated
number is conserved, while resolving within-bin Mie oscillations better than
a bin-center approximation.
"""
struct LinearIntegrationPerBin <: AbstractBinIntegration
    nquad::Int
    function LinearIntegrationPerBin(nquad::Integer = 160)
        nquad > 0 || throw(ArgumentError("nquad must be positive"))
        return new(Int(nquad))
    end
end

"""
    DirectBinSum()

Compatibility/diagnostic mode that evaluates Mie only at each bin center.
This is fast but not fine enough for production Mie over wide sectional bins;
prefer [`LinearIntegrationPerBin`](@ref) or [`ConstantIntegrationPerBin`](@ref).
"""
struct DirectBinSum <: AbstractBinIntegration end

# ============================================================================
# Refractive index hierarchy (parallel to RefractiveIndexLUT/Database from
# types.jl — those remain the IO format; these are the per-bin handles)
# ============================================================================

abstract type AbstractRefractiveIndex{FT} end

"""Wavelength-independent complex refractive index `n = nᵣ + i·nᵢ`."""
struct ConstantRI{FT} <: AbstractRefractiveIndex{FT}
    nᵣ::FT
    nᵢ::FT
end

"""
    WavelengthDependentRI{FT, I} <: AbstractRefractiveIndex{FT}

Tabulated RI vs λ (μm), with `nᵣ`/`nᵢ` callables (typically interpolators)
returning the real and imaginary parts at a given wavelength.
"""
struct WavelengthDependentRI{FT,I} <: AbstractRefractiveIndex{FT}
    λ::Vector{FT}
    nᵣ::I
    nᵢ::I
end

# ============================================================================
# Species composition + mixing rules
# ============================================================================

"""
    SpeciesComposition{FT}

Composition of one aerosol bin: a vector of species names, their fractions
(in either volume or mass space), and a `basis` symbol (`:volume` or
`:mass`). Fractions are expected to sum to 1 within numerical tolerance.
"""
struct SpeciesComposition{FT}
    species::Vector{String}
    fractions::Vector{FT}
    basis::Symbol
end

abstract type AbstractMixingRule end

"""External mixing: each species occupies its own particle population."""
struct ExternalMixing <: AbstractMixingRule end

"""Volume-weighted effective complex RI (simple internal mixing)."""
struct VolumeWeightedMixing <: AbstractMixingRule end

"""Maxwell-Garnett effective-medium mixing (stub — not implemented in this PR)."""
struct MaxwellGarnettMixing <: AbstractMixingRule end

"""Bruggeman effective-medium mixing (stub — not implemented in this PR)."""
struct BruggemanMixing <: AbstractMixingRule end

# ============================================================================
# Sectional-aerosol abstractions
#
# The contract these types enforce (any future TOMAS variant, MAAM, MOSAIC,
# custom sectional scheme) is: every sectional payload owns a size grid,
# a per-bin × per-layer number array, a species vector + per-bin × per-layer
# mass array, plus opt-in RI / mixing / integration policies. Downstream
# RT / Mie code reads ONLY through the generic accessors below — TOMAS-specific
# constants (IBINS, MOLWT, species order) live inside the TOMAS reader.
# ============================================================================

# `AbstractAerosolScheme{FT} <: AerosolScheme` is declared in types.jl so
# that the concrete `TOMASScheme{FT}` / `TwoMomentScheme{FT}` structs in
# the same file can subtype it. This file only adds the bin-data side of
# the hierarchy.

"""
    AerosolSizeGrid{FT}

Sectional size grid for any binned aerosol scheme. Holds bin edges
(optional — not every scheme produces explicit edges), bin centers, and
bin widths in the chosen coordinate.

# Fields
- `edges::Union{Nothing, Vector{FT}}`: Bin-edge values, length `nbins+1`,
  in `units`. `nothing` when the scheme reports centers only.
- `centers::Vector{FT}`: Bin-center values, length `nbins`, in `units`.
- `widths::Vector{FT}`: Bin widths, length `nbins`. For log-spaced bins
  this is `Δlog₁₀`; for linear bins it is the linear width; for mass-bins
  it is the mass width. The reader records what is natural for the scheme.
- `coordinate::Symbol`: One of `:diameter`, `:radius`, `:mass`. Tells
  downstream code how to interpret `centers`/`edges`.
- `spacing::Symbol`: `:linear`, `:log10`, `:mass_quadruple`, `:mass_double`,
  etc. Free-form — schemes document what they emit.
- `units::String`: Unit string for the coordinate (`"m"`, `"μm"`, `"nm"`,
  `"kg"`, …).
"""
struct AerosolSizeGrid{FT}
    edges::Union{Nothing, Vector{FT}}
    centers::Vector{FT}
    widths::Vector{FT}
    coordinate::Symbol
    spacing::Symbol
    units::String
end

Base.length(g::AerosolSizeGrid) = length(g.centers)

"""
    AbstractAerosolBinData{FT}

Base type for sectional bin payloads (one per scene cell). Concrete
subtypes carry the grid, prognostic number + per-species mass, and
opt-in RI / mixing / integration policies. The Mie kernel (deferred)
dispatches on these.

Downstream code uses the generic accessors below — never field access:
[`nbins`](@ref), [`nlayers`](@ref), [`species_list`](@ref),
[`bin_composition`](@ref).
"""
abstract type AbstractAerosolBinData{FT} end

"""
    SectionalAerosolData{FT, S} <: AbstractAerosolBinData{FT}

Generic sectional-aerosol payload for one scene cell. `S <: AbstractAerosolScheme{FT}`
identifies the originating scheme (TOMAS variant, MAAM, …). The fields
are agnostic to that scheme — the scheme parameter only carries metadata
(densities, hydrophilicity, GCHP variable prefixes) that downstream code
can consult when needed.

# Fields
- `scheme::S`: The originating scheme (e.g. `TOMASScheme{FT}`).
- `size_grid::AerosolSizeGrid{FT}`: Bin geometry (centers / edges / widths).
- `number::Matrix{FT}`: Particle number concentration (`#/cm³`), shape
  `(nbins, nlayers)`. Orientation TOA→BOA.
- `species::Vector{String}`: Species labels in storage order.
- `species_mass::Array{FT,3}`: Mass concentration (`μg/m³`), shape
  `(nbins, nlayers, nspecies)`.
- `metadata::NamedTuple`: Scheme-specific provenance (e.g. `(; variant=:tomas15)`).
- `ri_database::Union{Nothing, RefractiveIndexDatabase{FT}}`: RI vs λ,
  optional. Required for Mie computation.
- `mixing_rule::AbstractMixingRule`: How species RIs combine within a bin.
- `integration::AbstractBinIntegration`: Bin → Mie consumption strategy.
"""
struct SectionalAerosolData{FT, S<:AbstractAerosolScheme{FT}} <: AbstractAerosolBinData{FT}
    scheme::S
    size_grid::AerosolSizeGrid{FT}
    number::Matrix{FT}
    species::Vector{String}
    species_mass::Array{FT,3}
    metadata::NamedTuple
    ri_database::Union{Nothing, RefractiveIndexDatabase{FT}}
    mixing_rule::AbstractMixingRule
    integration::AbstractBinIntegration
end

# Generic accessors ---------------------------------------------------------

"""
    nbins(d::AbstractAerosolBinData) -> Int

Number of size bins. Generic accessor — works on any sectional payload.
"""
nbins(d::SectionalAerosolData) = length(d.size_grid)

"""
    nlayers(d::AbstractAerosolBinData) -> Int

Number of atmospheric layers represented in `d.number` / `d.species_mass`.
"""
nlayers(d::SectionalAerosolData) = size(d.number, 2)

"""
    species_list(d::AbstractAerosolBinData) -> Vector{String}

Species label vector (in storage order).
"""
species_list(d::SectionalAerosolData) = d.species

"""
    densities(scheme::AbstractAerosolScheme, species_codes) -> Vector{FT}

Per-species bulk density (kg/m³) under the scheme's density policy. Concrete
schemes provide this; the fallback errors out so consumers don't silently
inherit a wrong default.
"""
function densities end
function refractive_index_key end
function read_aerosol_cell end

_densities_error(scheme) = error("densities not implemented for aerosol scheme $(typeof(scheme))")
densities(scheme::AbstractAerosolScheme, _species_codes) = _densities_error(scheme)

refractive_index_key(::AbstractAerosolScheme, species_code::AbstractString) =
    String(species_code)

read_aerosol_cell(scheme::AbstractAerosolScheme, _args...; _kwargs...) =
    error("read_aerosol_cell not implemented for aerosol scheme $(typeof(scheme))")

"""
    bin_composition(d::AbstractAerosolBinData, ibin, ilev; basis=:volume)

Per-bin per-layer [`SpeciesComposition`](@ref). The `basis` kwarg selects
mass- or volume-fraction reporting; `:volume` (default) converts per-species
mass to volume using the scheme's [`densities`](@ref) before normalizing.

Zero-mass bins return a composition with all-zero fractions; the Mie
boundary treats this as "no aerosol contribution".
"""
function bin_composition(d::SectionalAerosolData{FT}, ibin::Integer, ilev::Integer;
                         basis::Symbol = :volume) where {FT}
    1 ≤ ibin ≤ nbins(d)   || throw(BoundsError(d, (ibin, ilev)))
    1 ≤ ilev ≤ nlayers(d) || throw(BoundsError(d, (ibin, ilev)))
    m = view(d.species_mass, ibin, ilev, :)

    if basis === :mass
        total = sum(m)
        fracs = total > zero(FT) ? FT.(collect(m ./ total)) :
                                   fill(zero(FT), length(d.species))
        return SpeciesComposition{FT}(copy(d.species), fracs, :mass)
    elseif basis === :volume
        ρ = densities(d.scheme, d.species)              # Vector{FT}, kg/m³
        v = m ./ ρ
        total = sum(v)
        fracs = total > zero(FT) ? FT.(collect(v ./ total)) :
                                   fill(zero(FT), length(d.species))
        return SpeciesComposition{FT}(copy(d.species), fracs, :volume)
    else
        throw(ArgumentError("bin_composition basis must be :mass or :volume, got :$basis"))
    end
end

# ============================================================================
# Effective-RI computation (the Mie boundary)
# ============================================================================

"""
    effective_ri(bin, db, λ, rule)

Effective complex refractive index for one bin under an internal-mixing
rule. The composition `bin` provides species fractions; `db` provides
per-species `Complex{FT}` RI at wavelength `λ` (μm).

- [`VolumeWeightedMixing`](@ref): volume-weighted complex average (the
  fractions must already be in `:volume` basis — checked at runtime).
- [`MaxwellGarnettMixing`](@ref) / [`BruggemanMixing`](@ref): not implemented
  in this PR.

For [`ExternalMixing`](@ref), use [`external_ri_components`](@ref) instead —
it returns per-species RI + weight pairs that the caller runs Mie over.
"""
function effective_ri(bin::SpeciesComposition{FT},
                      db::RefractiveIndexDatabase{FT},
                      λ::Real,
                      ::VolumeWeightedMixing) where {FT}
    bin.basis == :volume || throw(ArgumentError(
        "VolumeWeightedMixing requires SpeciesComposition with basis=:volume, " *
        "got :$(bin.basis)"))
    n = Complex{FT}(zero(FT), zero(FT))
    @inbounds for (sp, frac) in zip(bin.species, bin.fractions)
        n += frac * get_refractive_index(db, sp, λ)
    end
    return n
end

function effective_ri(bin::SpeciesComposition{FT},
                      db::RefractiveIndexDatabase{FT},
                      λ::Real,
                      rule::VolumeWeightedMixing,
                      scheme::AbstractAerosolScheme{FT}) where {FT}
    bin.basis == :volume || throw(ArgumentError(
        "VolumeWeightedMixing requires SpeciesComposition with basis=:volume, " *
        "got :$(bin.basis)"))
    n = Complex{FT}(zero(FT), zero(FT))
    @inbounds for (sp, frac) in zip(bin.species, bin.fractions)
        key = refractive_index_key(scheme, sp)
        n += frac * get_refractive_index(db, key, λ)
    end
    return n
end

function effective_ri(d::SectionalAerosolData{FT},
                      ibin::Integer, ilev::Integer,
                      db::RefractiveIndexDatabase{FT},
                      λ::Real,
                      rule::VolumeWeightedMixing = VolumeWeightedMixing()) where {FT}
    bin = bin_composition(d, ibin, ilev; basis=:volume)
    return effective_ri(bin, db, λ, rule, d.scheme)
end

effective_ri(::SpeciesComposition, ::RefractiveIndexDatabase, ::Real,
             ::MaxwellGarnettMixing) =
    error("MaxwellGarnettMixing is not implemented in this PR (deferred)")

effective_ri(::SpeciesComposition, ::RefractiveIndexDatabase, ::Real,
             ::BruggemanMixing) =
    error("BruggemanMixing is not implemented in this PR (deferred)")

"""
    external_ri_components(bin, db, λ) -> Vector{Tuple{Complex{FT}, FT}}

For external mixing, return per-species (RI, fraction) pairs at wavelength
`λ` (μm). The caller is expected to run Mie once per species and sum the
per-species optics weighted by `fraction`.
"""
function external_ri_components(bin::SpeciesComposition{FT},
                                db::RefractiveIndexDatabase{FT},
                                λ::Real) where {FT}
    return [(get_refractive_index(db, sp, λ), frac)
            for (sp, frac) in zip(bin.species, bin.fractions)]
end

function external_ri_components(bin::SpeciesComposition{FT},
                                db::RefractiveIndexDatabase{FT},
                                λ::Real,
                                scheme::AbstractAerosolScheme{FT}) where {FT}
    return [(get_refractive_index(db, refractive_index_key(scheme, sp), λ), frac)
            for (sp, frac) in zip(bin.species, bin.fractions)]
end

# ============================================================================
# Optics-computation entry point (stub — Mie kernel replacement deferred)
# ============================================================================

"""
    compute_aerosol_optics(data, λ_grid, pol_type, decomp_type) -> AerosolOptics

Entry point for converting an [`AbstractAerosolBinData`](@ref) payload into
single-scattering optics over a wavelength grid. In this PR the Mie kernel
is **not** implemented; this method emits a `@warn` and delegates to the
legacy placeholder so callers see a clear "deferred" signal.
"""
function compute_aerosol_optics(::SectionalAerosolData,
                                ::AbstractVector,
                                _pol_type, _decomp_type)
    @warn "compute_aerosol_optics(::SectionalAerosolData, ...) is a stub — " *
          "Mie kernel replacement is deferred to a follow-up PR" maxlog=1
    error("compute_aerosol_optics(::SectionalAerosolData, ...) not implemented")
end
