#=
Vertical-grid conventions for surface-pressure perturbations.

A retrieval that adjusts surface pressure p_surf must say what happens to the
rest of the vertical grid when p_surf changes. The two conventions below
describe the same atmosphere at the reference p_surf; they differ only in how
the grid is rebuilt at a trial p_surf and, consistently, in the analytic
surface-pressure tangents used by the linearized model.
=#

"""
    AbstractVerticalGridConvention

How the model's vertical grid responds to a change in surface pressure
`p_surf` while the layer temperatures, specific humidities, and trace-gas
mixing ratios are held fixed.

Concrete conventions:

- [`FixedPressureGrid`](@ref) (default): every pressure interface above the
  surface is fixed; only the bottom boundary `p_half[end] = p_surf` moves.
- [`FixedAltitudeGrid`](@ref): every interface is fixed in geometric altitude
  above the surface; a change in `p_surf` rescales every interface pressure.

A convention enters in two places, which must agree for the analytic Jacobian
to be the derivative of the forward model actually evaluated:

1. [`surface_pressure_interfaces`](@ref) rebuilds the interface pressures for
   a trial `p_surf`;
2. `model_from_parameters(LinMode(), params; vertical_grid=...)` selects the
   matching surface-pressure tangents (`psurf_profile_tangents`).
"""
abstract type AbstractVerticalGridConvention end

"""
    FixedPressureGrid()

Default convention. The interface pressures `p_half[1:end-1]` are fixed and
only the bottom boundary `p_half[end] = p_surf` moves. Consequently only the
bottom layer changes: its dry and H₂O columns, its geometric thickness, and
its mid-layer pressure `p_full[end]`, with `dp_full[end]/dp_surf = 1/2`.
Geometric altitudes above the bottom layer shift with the bottom-layer
thickness.

The tangent describes the final model grid. With `profile_reduction > 0` the
reduced interfaces are spaced from the top pressure to `p_surf`, so a rebuilt
reduced grid moves more than its bottom interface; the bottom-only tangent is
then the derivative of a different (fixed-final-grid) forward model.
"""
struct FixedPressureGrid <: AbstractVerticalGridConvention end

"""
    FixedAltitudeGrid()

Every interface is fixed in geometric altitude above the surface, with the
layer temperatures `T`, specific humidities `q`, and mixing ratios fixed on
those altitude layers. This is the natural convention for a fixed-height
observer (flux tower, aircraft at a fixed altitude), around which the
layering must not move when `p_surf` changes.

vSmartMOM's layer hydrostatics (`compute_atmos_profile_fields`) give each
layer the thickness

    Δzᵢ = (R Tᵢ / (g₀ Mᵢ)) · ln(p_half[i+1] / p_half[i]),

with constant `g₀` and a moist-air molar mass `Mᵢ` fixed by `qᵢ`. Holding
every `Δzᵢ`, `Tᵢ`, and `Mᵢ` fixed therefore fixes every pressure ratio
`p_half[i+1]/p_half[i]`, and the whole grid scales with the surface pressure:

    p_half[i](p_surf) = p_half_ref[i] · p_surf / p_surf_ref.

This is exact for the model's isothermal-layer, constant-gravity
hydrostatics; it is not a first-order approximation. It follows that, for
every layer `i`,

    dp_half[i]/dp_surf = p_half[i]/p_surf,   dp_full[i]/dp_surf = p_full[i]/p_surf,
    dN_dry,i/dp_surf   = N_dry,i/p_surf,     dN_H₂O,i/dp_surf   = N_H₂O,i/p_surf,
    dΔzᵢ/dp_surf       = 0,

because each hydrostatic column is proportional to its pressure thickness.
Every layer's cross sections therefore respond to `p_surf`, not only the
bottom layer's.

Aerosol profiles defined in geometric altitude above the surface
(`LogNormal(log(z₀), σ₀)`) keep their layer allocation, because no interface
moves in altitude. Aerosol profiles defined in pressure (`Normal(p₀, σp)`)
keep `p₀` and `σp` fixed in pressure and are redistributed as the layer
pressures scale.
"""
struct FixedAltitudeGrid <: AbstractVerticalGridConvention end

"""
    surface_pressure_interfaces(convention, p_half, p_surf)

Return the interface pressures [hPa, TOA to BOA] for a trial surface pressure
`p_surf`, starting from the reference interfaces `p_half` (whose last entry
is the reference surface pressure). The input is not modified.

At `p_surf == p_half[end]` both conventions return an exact copy of
`p_half`, so the reference atmosphere is identical bit for bit.
"""
function surface_pressure_interfaces end

function surface_pressure_interfaces(::FixedPressureGrid,
                                     p_half::AbstractVector{FT},
                                     p_surf::Real) where {FT}
    p_new = collect(p_half)
    p_new[end] = FT(p_surf)
    p_new[end] > p_new[end - 1] || throw(ArgumentError(
        "surface pressure $p_surf hPa is not below the final layer top " *
        "$(p_new[end - 1]) hPa"))
    return p_new
end

function surface_pressure_interfaces(::FixedAltitudeGrid,
                                     p_half::AbstractVector{FT},
                                     p_surf::Real) where {FT}
    ps = FT(p_surf)
    isfinite(ps) && ps > zero(FT) || throw(ArgumentError(
        "surface pressure must be finite and positive; got $p_surf hPa"))
    # One common factor for every interface: all layer pressure ratios, and
    # hence all hydrostatic layer thicknesses, are unchanged. The factor is
    # exactly one at the reference surface pressure.
    scale = ps / p_half[end]
    p_new = collect(p_half) .* scale
    p_new[end] = ps
    return p_new
end
