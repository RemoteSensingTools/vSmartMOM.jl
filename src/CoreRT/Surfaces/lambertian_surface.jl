#=
============================================================================
Lambertian (isotropic) surface BRDF
============================================================================

A Lambertian surface scatters incident flux equally in every outgoing
direction.  Its bidirectional reflectance is constant:

    ρ(μᵢ, μᵣ, Δϕ)  =  a / π

where `a` ∈ [0, 1] is the spectrally-varying albedo.  Because the BRDF has
no angular structure, only the m = 0 Fourier moment is non-zero — the
`m > 0` branch sets the surface added-layer to identity / zero.  The
factor of 2 in `ρ = 2·a` here is the `1/π · 2π` from converting the
hemispheric flux `a` to a Lambertian radiance × azimuthal-weight
compensation (the same trick is used in `inject_surface_SIF!` below).

Four flavours share the same equation, differing ONLY in how the spectral
albedo vector `a` is produced:
  • `LambertianSurfaceScalar`    — single-band scalar albedo (broadcast)
  • `LambertianSurfaceSpectrum`  — user-supplied per-grid-point vector
  • `LambertianSurfaceLegendre`  — albedo = Σ cₙ · Pₙ(λ̃) (spectral)
  • `LambertianSurfaceSpline`    — albedo from a spline interpolator

The albedo lives behind `surface_albedo(brdf, τ_sum)` (one method per
flavour); the surface-layer scaffold (`create_surface_layer!`) is written
ONCE for all four below.

────────────────────────────────────────────────────────────────────────
Surface-layer conventions (uniform across all four flavours)
────────────────────────────────────────────────────────────────────────
The surface is the *bottom* layer of the column, so only `r⁻⁺`, `r⁺⁻`
and `j₀⁻` can change the TOA reflectance `R` (see
`CoreKernel/interaction.jl`): the surface's `t⁻⁻` updates the composite
`T⁻⁻` only *after* `R⁻⁺` is formed in the final interaction and is never
read again, and `t⁺⁺` / `j₀⁺` flow into the composite `T⁺⁺` / `J₀⁺`
(the below-surface downward field). Consequently the choice of `t` and
`j₀⁺` is invisible to `R`, which is why the Scalar, Legendre/Spline and
linearized paths historically disagreed without breaking the golden TOA
suites. We unify them as follows:

  • `t⁺⁺ = t⁻⁻ = I` (identity) — matches the long-standing Scalar/generic
    forward path; keeps the bit-pinned Scalar TOA results unchanged.
  • `j₀⁺ = attenuated direct beam` (μ-selector · exp(-τ/μ₀)). This is the
    Scalar/generic convention and is REQUIRED: at m=0 the surface's `j₀⁺`
    is consumed by `interaction_hdrf!` (the BHR-downward flux diagnostic,
    `bhr_J₀⁺ += j₀⁺[iμ₀]·μ₀`). The old Legendre/Spline path zeroed it
    (a "Suniti double-check" stub) and silently broke that diagnostic.
  • `j₀⁻ = μ₀ · R_surf · I₀ · exp(-τ/μ₀)` (the reflected direct beam).
  • m > 0: everything zero (`r⁻⁺ = r⁺⁻ = 0`, `j₀± = 0`) with `t = I`.
    The pre-consolidation code zeroed `r⁻⁺` twice and never `r⁺⁻`
    (benign only because m=0 runs first and leaves r⁺⁻ already zero) —
    fixed here.
============================================================================
=#

# ---------------------------------------------------------------------------
# Per-flavour spectral albedo
# ---------------------------------------------------------------------------
"""
    surface_albedo(brdf, τ_sum) -> FT or Vector{FT}

Spectral (bare, un-doubled) Lambertian albedo `a` for each spectral point.
The surface scaffold multiplies this by 2 to form the m=0 Fourier moment
(`ρ = 2a`, the `1/π · 2π` flux→radiance/azimuth normalization). `τ_sum`
supplies the spectral axis length (one entry per spectral point).

Methods:
  • `LambertianSurfaceScalar`   — constant `albedo` (caller broadcasts).
  • `LambertianSurfaceSpectrum` — the stored vector (length must == nSpec).
  • `LambertianSurfaceLegendre` — Σ cₙ Pₙ(x), x ∈ range(-1, 1, nSpec).
    NOTE: the Legendre basis is evaluated over an *index*-space grid
    `range(-1, 1, length=nSpec)`, NOT over a physical wavelength grid.
    This is a historical quirk preserved bit-for-bit; the polynomial maps
    the spectral *index* to [-1, 1], so the same coefficients give a
    different curve for a different number of spectral points.
  • `LambertianSurfaceSpline`   — `interpolator(wlGrid)`.
"""
surface_albedo(b::LambertianSurfaceScalar{FT}, τ_sum) where {FT} = b.albedo

function surface_albedo(b::LambertianSurfaceSpectrum{FT}, τ_sum) where {FT}
    nSpec = length(τ_sum)
    length(b.albedo) == nSpec || throw(DimensionMismatch(
        "LambertianSurfaceSpectrum albedo length $(length(b.albedo)) ≠ " *
        "number of spectral points $nSpec"))
    return b.albedo
end

function surface_albedo(b::LambertianSurfaceLegendre{FT}, τ_sum) where {FT}
    n = length(b.legendre_coeff)
    # Constant term: P₀(x) ≡ 1, so the albedo is just the coefficient.
    # (compute_legendre_poly asserts nmax > 1, so the n == 1 case must be
    # handled analytically — it errored for as long as this surface existed.)
    n == 1 && return fill(b.legendre_coeff[1], length(τ_sum))
    # Legendre basis evaluated over INDEX space range(-1,1,nSpec) — see note above.
    x = collect(range(FT(-1), FT(1), length=length(τ_sum)))
    P = Scattering.compute_legendre_poly(x, n)[1]
    return P * b.legendre_coeff
end

surface_albedo(b::LambertianSurfaceSpline{FT}, τ_sum) where {FT} =
    b.interpolator(b.wlGrid)

# Float type carried by each flavour's albedo data.
_albedo_eltype(::LambertianSurfaceScalar{FT}) where {FT}   = FT
_albedo_eltype(::LambertianSurfaceSpectrum{FT}) where {FT} = FT
_albedo_eltype(::LambertianSurfaceLegendre{FT}) where {FT} = FT
_albedo_eltype(::LambertianSurfaceSpline{FT}) where {FT}   = FT

# ---------------------------------------------------------------------------
# Shared scaffold helpers (also used by the generic AbstractSurfaceType
# method in rpv_surface.jl — keeps the source fill / zero logic in one place)
# ---------------------------------------------------------------------------
"""
    _surface_source!(added_layer, R_surf, τ_sum, quad_points, pol_type, FT)

Fill the surface added-layer source vectors (`j₀⁺`, `j₀⁻`) for the m=0
moment under Source-Function-Integration. `R_surf` is the *dense* m=0
reflectance matrix (already carrying the `2a` factor — pre the
quadrature-weight multiply). `j₀⁺` is the attenuated direct beam (consumed
by `interaction_hdrf!`); `j₀⁻` is its reflection.
"""
function _surface_source!(added_layer, R_surf, τ_sum, quad_points, pol_type, ::Type{FT}) where {FT}
    (; qp_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points
    I₀_NquadN = similar(qp_μN)
    I₀_NquadN[:] .= zero(FT)
    I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀
    atten = exp.(-τ_sum / μ₀)'
    added_layer.j₀⁺[:, 1, :] .= I₀_NquadN .* atten
    added_layer.j₀⁻[:, 1, :] .= μ₀ * (R_surf * I₀_NquadN) .* atten
    return nothing
end

"""
    _fill_surface_layer!(added_layer, R_surf_weighted, T_surf)

Write the m=0 surface reflection/transmission into the added layer:
`r⁻⁺ = R_surf_weighted`, `r⁺⁻ = 0`, `t⁺⁺ = t⁻⁻ = T_surf` (identity).
`R_surf_weighted` already includes the `·Diag(qp_μN·wt_μN)` quadrature
weighting (and any spectral broadcast).
"""
function _fill_surface_layer!(added_layer, R_surf_weighted, T_surf)
    FT = eltype(added_layer.r⁻⁺)
    added_layer.r⁻⁺ .= R_surf_weighted
    added_layer.r⁺⁻ .= zero(FT)
    added_layer.t⁺⁺ .= T_surf
    added_layer.t⁻⁻ .= T_surf
    return nothing
end

"""
    _zero_surface_layer!(added_layer, T_surf)

m > 0 surface layer: a Lambertian/isotropic surface contributes nothing to
non-zero Fourier moments. Zero both reflectances and both source vectors,
leave transmission at identity. (Zeroes `r⁺⁻` explicitly — the old code
zeroed `r⁻⁺` twice and never `r⁺⁻`.)
"""
function _zero_surface_layer!(added_layer, T_surf)
    FT = eltype(added_layer.r⁻⁺)
    added_layer.r⁻⁺ .= zero(FT)
    added_layer.r⁺⁻ .= zero(FT)
    added_layer.t⁺⁺ .= T_surf
    added_layer.t⁻⁻ .= T_surf
    added_layer.j₀⁺ .= zero(FT)
    added_layer.j₀⁻ .= zero(FT)
    return nothing
end

# ---------------------------------------------------------------------------
# Unified Lambertian surface-layer scaffold (all four flavours)
# ---------------------------------------------------------------------------
"""
    create_surface_layer!(lambertian, added_layer, SFI, m, pol_type,
                          quad_points, τ_sum, architecture)

Computes (in place) surface optical properties for any Lambertian surface
flavour ([`LambertianSurfaceScalar`](@ref), `LambertianSurfaceSpectrum`,
`LambertianSurfaceLegendre`, `LambertianSurfaceSpline`) as an
[`AddedLayer`](@ref). The four differ only in how the spectral albedo is
produced (`surface_albedo`); the scaffold — including the m>0 / j₀⁺
conventions documented at the top of this file — is shared.

    - `lambertian` a Lambertian surface struct (one of the four flavours)
    - `SFI` bool if SFI (Source Function Integration) is used
    - `m` Fourier moment (starting at 0)
    - `pol_type` Polarization type struct
    - `quad_points` Quadrature points struct
    - `τ_sum` total optical thickness from TOA to the surface
    - `architecture` Compute architecture (GPU,CPU)
"""
function create_surface_layer!(lambertian::Union{LambertianSurfaceScalar,
                                                 LambertianSurfaceSpectrum,
                                                 LambertianSurfaceLegendre,
                                                 LambertianSurfaceSpline},
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI,
                               m::Int,
                               pol_type,
                               quad_points,
                               τ_sum,
                               architecture)

    FT = _albedo_eltype(lambertian)

    (; qp_μN, wt_μN) = quad_points
    arr_type = array_type(architecture)
    Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n
    T_surf = arr_type(Diagonal(ones(FT, pol_type.n * Nquad)))

    if m == 0
        # Bare albedo per spectral point, then the factor-2 m=0 normalization.
        # Spectral flavours (Spectrum/Legendre/Spline) return a HOST vector;
        # upload it BEFORE broadcasting against the device-resident τ_sum —
        # a mixed host/device broadcast errors on GPU. Scalar albedos broadcast
        # against device arrays directly.
        ρ₀ = surface_albedo(lambertian, τ_sum)
        ρ_dev = ρ₀ isa AbstractVector ? arr_type(FT.(ρ₀)) : FT(ρ₀)
        # `.+ zero.(τ_sum)` lifts a scalar albedo to a per-spectral vector.
        ρ = FT(2) .* (ρ_dev .+ zero.(τ_sum))

        # Unit Stokes-I selector replicated across quadrature directions
        # (albedo == 1); the spectral ρ is broadcast in afterwards — the
        # cheap pattern (no nSpec separate matrices materialized).
        R_unit = Matrix(Diagonal(vcat(FT(1), zeros(FT, pol_type.n - 1))))
        R_unit = repeat(R_unit', Nquad)
        R_unit = arr_type(repeat(R_unit', Nquad))

        # Source function of surface (m=0):
        #   j₀⁺ = attenuated direct beam  (consumed by interaction_hdrf!)
        #   j₀⁻ = μ₀ · (ρ · (R_unit·I₀)) · exp(-τ/μ₀)   — ρ broadcast per spec pt.
        # The grouping `μ₀ .* ((R_unit·I₀) .* ρ') .* exp'` is chosen so that
        # for a scalar albedo it is bit-identical to the historical Scalar
        # path `μ₀*(R_surf·I₀) .* exp'` (R_surf carried ρ on its diagonal):
        # IEEE multiplication is commutative, so baking ρ into the matvec
        # vs. broadcasting it in afterwards gives identical bits.
        if SFI
            (; iμ₀Nstart, iμ₀, μ₀) = quad_points
            I₀_NquadN = similar(qp_μN)
            I₀_NquadN[:] .= zero(FT)
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀
            added_layer.j₀⁺[:, 1, :] .= I₀_NquadN .* exp.(-τ_sum / μ₀)'
            added_layer.j₀⁻[:, 1, :] .= μ₀ .* ((R_unit * I₀_NquadN) .* ρ') .* exp.(-τ_sum / μ₀)'
        end

        # Spline pattern: weight the unit selector once, broadcast ρ in.
        R_weighted = (R_unit * Diagonal(qp_μN .* wt_μN)) .* reshape(ρ, 1, 1, :)
        _fill_surface_layer!(added_layer, R_weighted, T_surf)
    else
        _zero_surface_layer!(added_layer, T_surf)
    end
    return nothing
end

"""
    inject_surface_SIF!(brdf, added_layer, m, pol_type, SIF₀, architecture)

Add isotropic solar-induced fluorescence (SIF) surface emission to
`added_layer.j₀⁻`. SIF is Lambertian — only the m=0 Fourier moment carries
it — so higher moments are untouched. The factor 2 comes from
(1/π) × 2π: (1/π) normalizes the hemispheric SIF flux `SIF₀` into a
Lambertian radiance, and 2π compensates the `weight = 0.5/π` azimuthal
weighting applied downstream in `postprocessing_vza!` (SIF is isotropic
and must not be azimuthally weighted).

Ported from sanghavi `lambertian_surface.jl` (injection sites at
sanghavi lines 67-68 and 157-158). Non-Lambertian surfaces fall through
to a no-op.
"""
inject_surface_SIF!(::AbstractSurfaceType, _added_layer, _m, _pol_type, _SIF₀, _architecture) = nothing

function inject_surface_SIF!(
    ::Union{LambertianSurfaceScalar, LambertianSurfaceSpectrum,
            LambertianSurfaceLegendre, LambertianSurfaceSpline},
    added_layer::Union{AddedLayer, AddedLayerRS},
    m::Int,
    pol_type,
    SIF₀::AbstractArray,
    architecture,
)
    m == 0 || return nothing
    iszero(SIF₀) && return nothing
    FT = eltype(added_layer.j₀⁻)
    arr_type = array_type(architecture)
    Nquad = size(added_layer.j₀⁻, 1) ÷ pol_type.n
    added_layer.j₀⁻[:, 1, :] .+= FT(2) .* arr_type(repeat(FT.(SIF₀), Nquad))
    return nothing
end

"""
    _sif_source(RS_type)

Return `RS_type.SIF₀` if the field is declared, else `nothing`. Used by
`rt_run` / `rt_run_ss` to thread SIF into `inject_surface_SIF!` without
requiring every `AbstractRamanType` concrete to carry the field (e.g.
`_plus` variants have it commented out).
"""
_sif_source(RS_type) = hasproperty(RS_type, :SIF₀) ? RS_type.SIF₀ : nothing
