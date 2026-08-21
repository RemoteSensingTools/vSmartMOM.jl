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
struct ParameterLayout
    n_atmosphere::Int
    aerosol_params::Int
    n_aerosols::Int
    n_gases::Int
    n_surface::Int
    n_sif::Int
    n_canopy::Int
end

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
