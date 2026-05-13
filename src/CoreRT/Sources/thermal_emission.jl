#=

`ThermalEmission` — first-class per-layer Planck volume source for the
v0.6/v0.7 source-term framework. Represents isothermal in-atmosphere
thermal-IR emission `B(T_iz)` at every layer, contributing to both the
upwelling (`j₀⁻`) and downwelling (`j₀⁺`) source vectors of each
elemental added layer before doubling propagates the layer.

Mirrors the [`SurfaceSIF`](@ref) dispatch pattern:

  - Isotropic source → only the m=0 Fourier moment is non-zero.
  - The injection multiplies by `2π` so the uniform `0.5/π` weighting
    applied downstream by `postprocessing_vza!` (rt_run.jl:239,503) brings
    the contribution back to `B(T)` in radiance units. Compare to
    `SurfaceSIF` (surface_sif.jl:177) which uses a factor of `2` because
    `SIF₀` is *irradiance* and the `1/π` in the combined `(2π)·(1/π) = 2`
    converts irradiance → Lambertian radiance. `B(T)` is *already*
    radiance, so the conversion factor is bare `2π`.

**TIR weight rule (overrides REFACTOR_SPEC_v6 §2.11; corrected 2026-05-11):**
TIR is isotropic and does not undergo Fourier decomposition. The
`weight = 1.0; if(m==0) weight = 0.5` rule in legacy Fortran
`rt_atmo_path.f90:132-134` (and the uniform `0.5/π` in Julia
`rt_run.jl:239,503`) applies to azimuthally-varying terms, not to thermal.
Thermal contributes with weight 1 at m=0 and 0 at m>0; the gating
`m == 0 || return nothing` plus the explicit `2π` factor at injection
together produce that semantics through the existing reconstruction loop.

The legacy `BD_scripts/molec_opac.py:12` "multiplied by 2" comment is *not*
about this Fourier weight — it refers to the disk-integration normalization
(half-ellipsoid vs full-ellipsoid surface area) which is an ExoOptics
geometry fix, not an RT-side fix.

# Per-layer thermal source formula

For an isothermal absorbing-emitting elemental slab at temperature `T_iz`,
total optical depth `dτ_λ` per spectral point and SSA `ϖ_λ`, the exact
thermal source contribution to outgoing stream `μᵢ` (both upward and
downward) is:

```
j_thermal[i, n] = (1 - ϖ_λ[n]) · B(T_iz, ν[n]) · (1 - exp(-dτ_λ[n] / μᵢ))
```

This matches the Fell (1997, FU Berlin PhD thesis) finite-δ thermal source
with the L'Hôpital limit handled by `-expm1(-x)`, consistent with the
elastic kernel's exact-finite-δ design (elemental.jl Eqs. 1.52–1.56). In
the linear-thin-layer limit (`dτ/μᵢ → 0`) this reduces to the legacy
Fortran formula `j_elt = dτ · (1-ϖ) · B / μ` from `rt_elemental.f90:226`.

The source is unpolarized — only the I (Stokes-1) component is filled;
Q/U/V remain zero. Doubling propagates the per-elemental contributions to
the full layer source, and interaction stitches the layers together.

=#

"""
    ThermalEmission(; B_layer = nothing) <: AbstractSource

User-facing per-layer thermal volume source. Carries a pre-computed Planck
radiance matrix `B_layer` of shape `(nLayers, nSpec)` in
**mW · m⁻² · sr⁻¹ · cm⁻¹** (matching [`vSmartMOM.SolarModel.planck_spectrum_wn`](@ref)).

Two construction patterns:

```julia
# Pre-computed B_layer
ThermalEmission(B_layer = my_planck_matrix)

# T_layers + spectral grid → computes B_layer internally
ThermalEmission(T_layers, ν_grid)
```

Compose with other sources via `+`:

```julia
sources = SolarBeam() + ThermalEmission(T_layers, ν_grid)
```

# Fields
- `B_layer :: Union{Nothing, AbstractMatrix}`: Planck radiance per layer
  per spectral point with shape `(nLayers, nSpec)`. Stored without an
  `eltype` constraint — `prepare_source` does the precision conversion
  against the model's `FT` at solve time, matching [`SolarBeam`](@ref)'s
  FT-deferred design.

# AD mode
[`source_ad_mode`](@ref) returns [`AnalyticSourceJacobian`](@ref); a future
`source_tangent!` body provides the analytic
`(∂j/∂T_layer, ∂j/∂τ, ∂j/∂ϖ)` derivatives (Phase A.2).
"""
struct ThermalEmission <: AbstractSource
    B_layer :: Union{Nothing, AbstractMatrix}
end

ThermalEmission(; B_layer=nothing) = ThermalEmission(B_layer)

"""
    ThermalEmission(T_layers::AbstractVector, ν_grid::AbstractVector) -> ThermalEmission

Convenience constructor: compute `B_layer[iz, n] = B(T_layers[iz], ν_grid[n])`
via [`vSmartMOM.SolarModel.planck_spectrum_wn`](@ref), in
mW · m⁻² · sr⁻¹ · cm⁻¹.
"""
function ThermalEmission(T_layers::AbstractVector{<:Real},
                         ν_grid::AbstractVector{<:Real})
    nLayers = length(T_layers)
    ν = collect(float.(ν_grid))
    nSpec = length(ν)
    B_layer = zeros(Float64, nLayers, nSpec)
    @inbounds for iz in 1:nLayers
        B_layer[iz, :] .= vSmartMOM.SolarModel.planck_spectrum_wn(float(T_layers[iz]), ν)
    end
    return ThermalEmission(B_layer)
end

source_ad_mode(::ThermalEmission) = AnalyticSourceJacobian()

Base.show(io::IO, s::ThermalEmission) =
    print(io, "ThermalEmission(B_layer=",
          s.B_layer === nothing ? "zeros" : summary(s.B_layer), ")")

"""
    PreparedThermalEmission{FT, AT} <: AbstractPreparedSource

Kernel-ready per-layer Planck source. `B_layer` is materialised on the
model's array type at the right `(nLayers, nSpec)` shape and `FT`
precision.
"""
struct PreparedThermalEmission{FT<:AbstractFloat, AT<:AbstractMatrix} <: AbstractPreparedSource
    B_layer :: AT
end

source_ad_mode(::PreparedThermalEmission) = AnalyticSourceJacobian()

# v0.7 Phase A.2a — thermal lives in its own j₀_by_src slot keyed `:thermal`.
# `source_expk_init` defaults to `ones` (inherited from `AbstractPreparedSource`
# in Sources/types.jl), which is correct for thermal — the bottom sub-layer's
# emission is not pre-attenuated by the top, so the doubling recurrence
# `j₁⁺ = j₀⁺ · expk` collapses to `j₁⁺ = j₀⁺` (Fortran rt_doubling.f90:196).
source_key(::PreparedThermalEmission) = :thermal

Base.show(io::IO, p::PreparedThermalEmission) =
    print(io, "PreparedThermalEmission(B_layer=", summary(p.B_layer), ")")

"""
    prepare_source(s::ThermalEmission, FT, pol_n, nSpec, arr_type) -> PreparedThermalEmission

Resolve a [`ThermalEmission`](@ref) into a kernel-ready
[`PreparedThermalEmission`](@ref). The default (`B_layer === nothing`)
materialises a zero `(1, nSpec)` matrix on the active architecture (no-op
source, useful as a placeholder). A user-supplied `B_layer` is
precision-converted; the leading `nLayers` dimension is taken from the
matrix itself rather than from `pol_n` since thermal is unpolarized.

Note that `pol_n` is unused — included only to keep the
[`prepare_source`](@ref) contract uniform across source types.
"""
function prepare_source(s::ThermalEmission, FT::Type{<:AbstractFloat},
                        pol_n::Integer, nSpec::Integer, arr_type)
    if s.B_layer === nothing
        B = zeros(FT, 1, nSpec)
        return PreparedThermalEmission{FT, typeof(arr_type(B))}(arr_type(B))
    else
        size(s.B_layer, 2) == nSpec || error(
            "ThermalEmission: B_layer column count $(size(s.B_layer, 2)) does not " *
            "match nSpec=$nSpec. The matrix must be (nLayers, nSpec).")
        B_dev = arr_type(convert(Array{FT, 2}, s.B_layer))
        return PreparedThermalEmission{FT, typeof(B_dev)}(B_dev)
    end
end

# ============================================================================
# Volume-source dispatch — `contribute!` per-layer
# ============================================================================
#
# The volume source `contribute!` is called once per layer during the
# rt_kernel layer loop, between `elemental!` and `doubling!`. For
# `NoSource` / `SourceSet` / non-volume sources (`PreparedSolarBeam`,
# `PreparedSurfaceSIF`) it is a no-op. For `PreparedThermalEmission` it
# fills the per-elemental thermal contribution.
#
# Dispatch signature (intentionally explicit rather than var-args to make
# the `iz`/`m`/`pol_type`/`quad_points`/`architecture` arguments part of
# the named API):
#
#   contribute!(prep, added_layer, ϖ_λ, dτ_λ,
#               iz, m, pol_type, quad_points, architecture)
#
# ============================================================================

# Default no-op for non-volume prepared sources. SolarBeam already
# contributes inside the elemental kernel via the SFI path; SurfaceSIF
# contributes at the surface via `surface_source_contribute!`.
contribute!(::AbstractPreparedSource, _added_layer,
            _ϖ_λ, _dτ_λ, _iz, _m, _pol_type, _quad_points, _architecture) = nothing

# NoSource and SourceSet dispatchers — these replace the older var-args
# stubs at the bottom of solar_beam.jl for the volume-source seam.
contribute!(::NoSource, _added_layer,
            _ϖ_λ, _dτ_λ, _iz, _m, _pol_type, _quad_points, _architecture) = nothing

function contribute!(s::SourceSet, added_layer,
                     ϖ_λ, dτ_λ, iz::Integer, m::Integer,
                     pol_type, quad_points, architecture)
    @inbounds for src in s.sources
        contribute!(src, added_layer, ϖ_λ, dτ_λ, iz, m, pol_type, quad_points, architecture)
    end
    return nothing
end

"""
    contribute!(prep::PreparedThermalEmission, added_layer,
                ϖ_λ, dτ_λ, iz, m, pol_type, quad_points, architecture)

Inject the per-layer Planck volume source into the elemental added layer's
`j₀⁺` and `j₀⁻` accumulators, in the I (Stokes-1) component only, at the
m=0 Fourier moment.

For each stream `μᵢ` and spectral point `n`:

```
δ = 2π · (1 - ϖ_λ[n]) · B_layer[iz, n] · (1 - exp(-dτ_λ[n] / μᵢ))
j₀⁺[i_I, 1, n] += δ
j₀⁻[i_I, 1, n] += δ
```

where `i_I = (iμ - 1) · pol_type.n + 1` indexes the Stokes-I row for
quadrature stream `iμ`. Q/U/V rows are untouched (thermal Planck is
unpolarized).

The `2π` factor undoes the uniform `0.5/π` weight applied by
`postprocessing_vza!` so that `B(T)` survives at full strength after
reconstruction. This is exactly analogous to the factor of `2` used in
`surface_sif.jl:177`, with the extra `π` accounting for `B(T)` being a
radiance rather than a Lambertian-disk irradiance.

For `m > 0` this is a no-op — thermal is isotropic and contributes only at
m=0.
"""
function contribute!(prep::PreparedThermalEmission, added_layer,
                     ϖ_λ::AbstractVector, dτ_λ::AbstractVector,
                     iz::Integer, m::Integer,
                     pol_type, quad_points, architecture)
    m == 0 || return nothing
    # Skip on placeholder zero-matrix (B_layer === zeros) to avoid the
    # Nquad·nSpec write loop entirely when the user passed no thermal data.
    iszero(prep.B_layer) && return nothing
    # When the prepared source was built with a placeholder default
    # (1 × nSpec), the layer index iz can exceed `size(B_layer, 1)` —
    # in that case treat the layer as having zero Planck contribution.
    iz <= size(prep.B_layer, 1) || return nothing
    # v0.7 Phase A.2a — write to the per-source `:thermal` slot in
    # `added_layer.j₀_by_src` (not the legacy solar j₀⁺/j₀⁻ fields). The slot
    # must exist; if `make_added_layer` was called without a thermal source
    # in `prepared_sources`, the slot won't be there and we return silently
    # (defensive — should be reachable only via test/stub callers).
    slot = get(added_layer.j₀_by_src, :thermal, nothing)
    slot === nothing && return nothing

    FT = eltype(slot.j₀⁺)
    qp_μ = quad_points.qp_μ
    Nquad = length(qp_μ)
    nStokes = pol_type.n

    # Pull B_layer for this layer onto the host (small slice, copy-out is
    # cheap on GPU; the broadcast below keeps everything else on device).
    B_iz = Array(view(prep.B_layer, iz, :))           # (nSpec,)
    ϖ_h  = Array(ϖ_λ)                                  # (nSpec,)
    dτ_h = Array(dτ_λ)                                 # (nSpec,)
    μ_h  = Array(qp_μ)                                 # (Nquad,)

    nSpec = length(B_iz)
    # Build the (Nquad, nSpec) thermal contribution slab on host then
    # transfer to the active arch in one shot.
    j_th = zeros(FT, Nquad, nSpec)
    twoπ = FT(2π)
    @inbounds for n in 1:nSpec
        omega_abs = one(FT) - FT(ϖ_h[n])
        Bn = FT(B_iz[n])
        coeff = twoπ * omega_abs * Bn
        for iμ in 1:Nquad
            μi = FT(μ_h[iμ])
            μi > eps(FT) || continue
            j_th[iμ, n] = coeff * (-expm1(-FT(dτ_h[n]) / μi))
        end
    end

    arr_type = array_type(architecture)
    j_th_dev = arr_type(j_th)

    # NB: the slot is zeroed by `rt_kernel!` immediately before this
    # contribute! call — see the per-source slot-reset block at the top
    # of the noRS{FT} `rt_kernel!` body. Writing only the Stokes-I rows
    # here is therefore safe (Q/U/V stay at zero from that reset).
    @inbounds for iμ in 1:Nquad
        i_I = (iμ - 1) * nStokes + 1
        slot.j₀⁺[i_I, 1, :] .= j_th_dev[iμ, :]
        slot.j₀⁻[i_I, 1, :] .= j_th_dev[iμ, :]
    end
    return nothing
end

# ============================================================================
# `prepare_sources` already delegates per-source via `map(...)` in
# solar_beam.jl. No changes needed there — `prepare_source(::ThermalEmission, ...)`
# above plugs into that walk automatically.
# ============================================================================

# Predicate: does this prepared-source set carry a thermal volume source?
# Mirrors `has_solar_beam`/`has_surface_sif` in the sibling source files.
has_thermal_emission(::NoSource) = false
has_thermal_emission(::PreparedThermalEmission) = true
has_thermal_emission(::AbstractPreparedSource) = false

function has_thermal_emission(s::SourceSet)
    @inbounds for src in s.sources
        src isa PreparedThermalEmission && return true
    end
    return false
end
