"""
    ParameterLayout

Describes the ordering of physical parameters in the Jacobian derivative dimension.

Instead of hardcoding `7*NAer + NGas + NSurf` throughout the codebase, all index
arithmetic goes through this struct.  Each aerosol carries `aerosol_params`
sub-parameters (currently 7: `τ_ref, nᵣ, nᵢ, rₘ, σ_g, p₀, σ_p`), followed by
one slot per gas VMR and one per surface parameter.

# Example
```julia
layout = ParameterLayout(n_aerosols=1, n_gases=2, n_surface=1)
aerosol_range(layout, 1)   # 1:7
gas_range(layout)          # 8:9
surface_range(layout)      # 10:10
n_total(layout)            # 10
```
"""
struct ParameterLayout
    aerosol_params::Int
    n_aerosols::Int
    n_gases::Int
    n_surface::Int
end

ParameterLayout(; aerosol_params::Int=7, n_aerosols::Int=0, n_gases::Int=0,
                  n_surface::Int=1) =
    ParameterLayout(aerosol_params, n_aerosols, n_gases, n_surface)

"Total number of retrieval parameters."
@inline n_total(pl::ParameterLayout) =
    pl.aerosol_params * pl.n_aerosols + pl.n_gases + pl.n_surface

"Index range for aerosol `iaer` (1-based)."
@inline function aerosol_range(pl::ParameterLayout, iaer::Int)
    offset = (iaer - 1) * pl.aerosol_params
    return (offset + 1):(offset + pl.aerosol_params)
end

"Index range for all gas VMR parameters."
@inline function gas_range(pl::ParameterLayout)
    start = pl.aerosol_params * pl.n_aerosols + 1
    return start:(start + pl.n_gases - 1)
end

"Index range for surface parameters."
@inline function surface_range(pl::ParameterLayout)
    start = pl.aerosol_params * pl.n_aerosols + pl.n_gases + 1
    return start:(start + pl.n_surface - 1)
end

"Index of a specific surface parameter (1-based within the surface block)."
@inline surface_index(pl::ParameterLayout, i::Int=1) =
    pl.aerosol_params * pl.n_aerosols + pl.n_gases + i

"Number of layer-level parameters (aerosol + gas, excluding surface)."
@inline n_layer_params(pl::ParameterLayout) =
    pl.aerosol_params * pl.n_aerosols + pl.n_gases
