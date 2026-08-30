"Common interface for full and retrieval-selected Jacobian layouts."
abstract type AbstractParameterLayout end

"""
    ParameterLayout

Describes the ordering of physical parameters in the Jacobian derivative dimension.

Instead of hardcoding `1 + 7*NAer + NGas + NSurf` throughout the codebase, all index
arithmetic goes through this struct. Each aerosol carries `aerosol_params`
sub-parameters (currently 7: `τ_ref, nᵣ, nᵢ, rₘ, σ_g`, profile location,
profile width), followed by the layer-resolved gas VMR columns and surface parameters.
Gas columns are flattened species-major: all TOA-to-BOA layers for gas 1,
then all layers for gas 2, and so on. Surface parameters are followed by the
optional two-column SIF block `[SIF755, slope]`, then canopy parameters.
The final two columns are `(p₀, σ_p)` for a pressure-form `Normal` profile and
`(z₀, σ₀)` for an altitude-form `LogNormal` profile.

# Example
```julia
layout = ParameterLayout(n_aerosols=1, n_gases=6, n_surface=1) # 2 gases × 3 layers
psurf_index(layout)        # 1
aerosol_range(layout, 1)   # 2:8
gas_profile_range(layout, 1, 3) # 9:11
gas_layer_index(layout, 2, 3, 3) # 14
surface_range(layout)      # 15:15
n_total(layout)            # 15
```
"""
struct ParameterLayout <: AbstractParameterLayout
    n_atmosphere::Int
    aerosol_params::Int
    n_aerosols::Int
    n_gases::Int
    n_surface::Int
    n_sif::Int
    n_canopy::Int
end

"""
    ParameterKey

Stable, retrieval-independent identity for one physical Jacobian parameter.
`band == 0` denotes a parameter shared by all bands and `layer == 0` denotes
a parameter without a vertical-layer index.  Retrieval flavours use these
keys to compile small band-local tangent spaces into one global state vector.

`kind` identifies the parameter family; `component` identifies an aerosol,
gas, surface, or source within that family; `field` identifies the physical
quantity; and `layer`/`band` locate it vertically and spectrally. The key is
metadata only: unit and coordinate transformations remain the responsibility
of the retrieval boundary.
"""
struct ParameterKey
    kind::Symbol
    component::Union{Nothing,Int,String}
    field::Symbol
    layer::Int
    band::Int
end

ParameterKey(kind::Symbol, field::Symbol;
             component::Union{Nothing,Int,String}=nothing,
             layer::Integer=0, band::Integer=0) =
    ParameterKey(kind, component, field, Int(layer), Int(band))

"""
    ActiveParameterLayout

Compact parameter layout used by one spectral band of a selective
linearized calculation. `keys` are ordered exactly as the trailing tangent
dimension propagated by the RT kernels. `local_to_global` scatters those
columns into `global_keys`; `native_layer_columns` selects the corresponding
columns from the historical full physical-optics layout before the expensive
adding-doubling propagation begins.

Layer parameters must precede surface/source parameters. Native layer columns
must be positive, unique, sorted, and remain in the historical component
order (pressure, aerosols, gas); `local_to_global` supplies any desired
retrieval reordering. Surface and SIF columns must each be contiguous because
their analytic seed routines accept unit ranges.
"""
struct ActiveParameterLayout <: AbstractParameterLayout
    keys::Vector{ParameterKey}
    parameter_names::Vector{String}
    global_keys::Vector{ParameterKey}
    global_parameter_names::Vector{String}
    local_to_global::Vector{Int}
    native_layer_columns::Vector{Int}
    n_layer::Int
    surface_columns::UnitRange{Int}
    sif_columns::UnitRange{Int}
    canopy_columns::UnitRange{Int}
end

function ActiveParameterLayout(keys::Vector{ParameterKey},
                               parameter_names::Vector{String},
                               global_keys::Vector{ParameterKey},
                               global_parameter_names::Vector{String},
                               native_layer_columns::Vector{Int},
                               n_layer::Integer;
                               surface_columns::UnitRange{Int}=1:0,
                               sif_columns::UnitRange{Int}=1:0,
                               canopy_columns::UnitRange{Int}=1:0)
    length(keys) == length(parameter_names) || throw(DimensionMismatch(
        "active parameter keys and names must have equal lengths"))
    length(global_keys) == length(global_parameter_names) || throw(DimensionMismatch(
        "global parameter keys and names must have equal lengths"))
    length(native_layer_columns) == n_layer || throw(DimensionMismatch(
        "native_layer_columns must contain one entry per active layer parameter"))
    all(>(0), native_layer_columns) || throw(ArgumentError(
        "native layer columns must be positive"))
    length(unique(native_layer_columns)) == length(native_layer_columns) ||
        throw(ArgumentError("native layer columns must be unique"))
    issorted(native_layer_columns) || throw(ArgumentError(
        "native layer columns must retain native component order; use " *
        "local_to_global to impose a different retrieval-state order"))
    0 <= n_layer <= length(keys) || throw(ArgumentError(
        "n_layer must lie between zero and the local parameter count"))
    length(unique(keys)) == length(keys) || throw(ArgumentError(
        "active parameter keys must be unique within a band"))
    length(unique(global_keys)) == length(global_keys) || throw(ArgumentError(
        "global parameter keys must be unique"))
    lookup = Dict(key => i for (i, key) in enumerate(global_keys))
    local_to_global = Int[]
    for key in keys
        haskey(lookup, key) || throw(ArgumentError(
            "band-local parameter key $key is absent from the global layout"))
        push!(local_to_global, lookup[key])
    end
    covered = vcat(collect(1:Int(n_layer)), collect(surface_columns),
                   collect(sif_columns), collect(canopy_columns))
    sort(covered) == collect(1:length(keys)) || throw(ArgumentError(
        "layer, surface, SIF, and canopy blocks must cover each local column exactly once"))
    return ActiveParameterLayout(keys, parameter_names, global_keys,
        global_parameter_names, local_to_global, native_layer_columns,
        Int(n_layer), surface_columns, sif_columns, canopy_columns)
end

"Base type for retrieval-specific Jacobian selection policies."
abstract type AbstractJacobianFlavor end

"""
    OCO_RRS_synth()

Selective Jacobian flavour for the synthetic three-band OCO RRS/XCO2
retrieval. It selects the native physical columns corresponding to surface
pressure, the non-fixed CO2 layers, aerosol reference optical depth and
profile location, three Legendre surface coefficients per band, and the two
O2-A-band SIF parameters. Conversion from native aerosol `tau_ref`/`z0` to
the retrieval's `log(AOD760)`/`log(z0)` coordinates is intentionally applied
at the retrieval boundary. Aerosol microphysics and profile widths, H2O, and
the fixed upper CO2 layers remain forward-model inputs only.
"""
struct OCO_RRS_synth <: AbstractJacobianFlavor end

"""
    JacobianPlan

A compiled global named physical basis plus one compact
[`ActiveParameterLayout`](@ref) per spectral band. `global_keys` determine
cross-band identity and order; `bands[i]` determines the derivative dimension
actually allocated and propagated by band `i`. Retrieval-specific policy is
fully compiled into this object before the RT kernels are entered.
"""
struct JacobianPlan{F<:AbstractJacobianFlavor}
    flavor::F
    global_keys::Vector{ParameterKey}
    global_parameter_names::Vector{String}
    bands::Vector{ActiveParameterLayout}
end

@inline band_layout(plan::JacobianPlan, i_band::Integer) = plan.bands[i_band]
@inline n_global(plan::JacobianPlan) = length(plan.global_keys)
@inline parameter_names(layout::ActiveParameterLayout) = layout.parameter_names
@inline parameter_names(plan::JacobianPlan) = plan.global_parameter_names
@inline local_to_global(layout::ActiveParameterLayout) = layout.local_to_global
@inline native_layer_columns(layout::ActiveParameterLayout) = layout.native_layer_columns

ParameterLayout(; n_atmosphere::Int=1, aerosol_params::Int=7, n_aerosols::Int=0, n_gases::Int=0,
                  n_surface::Int=1, n_sif::Int=0, n_canopy::Int=0) =
    ParameterLayout(n_atmosphere, aerosol_params, n_aerosols, n_gases, n_surface, n_sif, n_canopy)

"Total number of retrieval parameters."
@inline n_total(pl::ParameterLayout) =
    pl.n_atmosphere + pl.aerosol_params * pl.n_aerosols + pl.n_gases + pl.n_surface + pl.n_sif + pl.n_canopy

"Atmospheric-state range; column 1 is surface pressure."
@inline atmosphere_range(pl::ParameterLayout) = 1:pl.n_atmosphere
@inline psurf_index(::ParameterLayout) = 1

"Index range for aerosol `iaer` (1-based)."
@inline function aerosol_range(pl::ParameterLayout, iaer::Int)
    offset = pl.n_atmosphere + (iaer - 1) * pl.aerosol_params
    return (offset + 1):(offset + pl.aerosol_params)
end

"Index range for all gas VMR parameters."
@inline function gas_range(pl::ParameterLayout)
    start = pl.n_atmosphere + pl.aerosol_params * pl.n_aerosols + 1
    return start:(start + pl.n_gases - 1)
end

"Index range for all vertical layers of gas species `igas` (species-major ordering)."
@inline function gas_profile_range(pl::ParameterLayout, igas::Int, n_layers::Int)
    n_layers > 0 || throw(ArgumentError("n_layers must be positive"))
    start = first(gas_range(pl)) + (igas - 1) * n_layers
    stop = start + n_layers - 1
    stop <= last(gas_range(pl)) || throw(BoundsError(gas_range(pl), start:stop))
    return start:stop
end

"Jacobian index for gas species `igas` in atmospheric layer `iz`."
@inline gas_layer_index(pl::ParameterLayout, igas::Int, iz::Int, n_layers::Int) =
    gas_profile_range(pl, igas, n_layers)[iz]

"Index range for surface parameters."
@inline function surface_range(pl::ParameterLayout)
    start = pl.n_atmosphere + pl.aerosol_params * pl.n_aerosols + pl.n_gases + 1
    return start:(start + pl.n_surface - 1)
end

"Index of a specific surface parameter (1-based within the surface block)."
@inline surface_index(pl::ParameterLayout, i::Int=1) =
    pl.n_atmosphere + pl.aerosol_params * pl.n_aerosols + pl.n_gases + i

"Number of retrieval parameters contributed by a surface model."
surface_parameter_count(::AbstractSurfaceType) = 1
surface_parameter_count(s::LambertianSurfaceLegendre) = length(s.legendre_coeff)

"Index range for SIF parameters: per-wavenumber radiance at 760 nm, then dSIF/dν."
@inline function sif_range(pl::ParameterLayout)
    start = last(surface_range(pl)) + 1
    return start:(start + pl.n_sif - 1)
end

@inline sif755_index(pl::ParameterLayout) = first(sif_range(pl))
@inline sif_slope_index(pl::ParameterLayout) = first(sif_range(pl)) + 1
@inline sif760_index(pl::ParameterLayout) = first(sif_range(pl))
@inline msif_index(pl::ParameterLayout) = first(sif_range(pl)) + 1

"Number of layer-level parameters (aerosol + gas, excluding surface and canopy)."
@inline n_layer_params(pl::ParameterLayout) =
    pl.n_atmosphere + pl.aerosol_params * pl.n_aerosols + pl.n_gases

"""Index range for canopy parameters (LAI, leaf_R, leaf_T, ...)."""
@inline function canopy_range(pl::ParameterLayout)
    start = pl.n_atmosphere + pl.aerosol_params * pl.n_aerosols + pl.n_gases + pl.n_surface + pl.n_sif + 1
    return start:(start + pl.n_canopy - 1)
end

# Active layouts deliberately expose only operations that remain meaningful
# after heterogeneous parameter selection. Retrieval code should use the
# stable ParameterKey metadata rather than infer semantics from offsets.
@inline n_total(pl::ActiveParameterLayout) = length(pl.keys)
@inline n_layer_params(pl::ActiveParameterLayout) = pl.n_layer
@inline surface_range(pl::ActiveParameterLayout) = pl.surface_columns
@inline surface_index(pl::ActiveParameterLayout, i::Int=1) = pl.surface_columns[i]
@inline sif_range(pl::ActiveParameterLayout) = pl.sif_columns
@inline sif755_index(pl::ActiveParameterLayout) = pl.sif_columns[1]
@inline sif_slope_index(pl::ActiveParameterLayout) = pl.sif_columns[2]
@inline sif760_index(pl::ActiveParameterLayout) = pl.sif_columns[1]
@inline msif_index(pl::ActiveParameterLayout) = pl.sif_columns[2]
@inline canopy_range(pl::ActiveParameterLayout) = pl.canopy_columns

function atmosphere_range(pl::ActiveParameterLayout)
    indices = findall(k -> k.kind === :atmosphere, pl.keys)
    isempty(indices) && return 1:0
    return first(indices):last(indices)
end

function psurf_index(pl::ActiveParameterLayout)
    index = findfirst(k -> k.kind === :atmosphere &&
                          k.field === :surface_pressure, pl.keys)
    index === nothing && throw(ArgumentError(
        "surface pressure is inactive in this Jacobian layout"))
    return index
end

@inline function _active_kind_indices(pl::ActiveParameterLayout, kind::Symbol)
    return findall(k -> k.kind === kind, pl.keys)
end

gas_range(pl::ActiveParameterLayout) = _active_kind_indices(pl, :gas)

function aerosol_range(pl::ActiveParameterLayout, iaer::Int)
    return findall(k -> k.kind === :aerosol && k.component == iaer, pl.keys)
end
