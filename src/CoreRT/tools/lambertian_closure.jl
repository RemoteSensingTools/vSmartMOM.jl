#=
============================================================================
LambertianClosure — analytic O(1)-per-albedo Lambertian surface closure
============================================================================

`rt_run_surface` (`rt_run_split.jl`) replays a full `create_surface_layer!`
+ `interaction!` per surface — exact for any BRDF, but it still pays an
`O(NquadN^3 · nSpec)` matrix solve per Fourier moment. For a Lambertian
surface the m=0 reflectance is rank-one
(`r⁻⁺ = 2a·u·vᵀ`, `Surfaces/lambertian_surface.jl:264-307`), so the final
adding step (`CoreKernel/interaction.jl`'s
`J₀⁻ ← J₀⁻ + T⁻⁻·(E − r⁻⁺R⁺⁻)⁻¹·(r⁻⁺J₀⁺ + j₀⁻)`) collapses via
Sherman-Morrison to a scalar rational function of the albedo `a` — once
three scalars/vectors are read from the cached m=0 atmosphere block, no
matrix algebra is needed at all. See `proposals/surface_split_albedo_sweep.md`
§3 for the full derivation; this file implements it.

Derivation summary (all quantities per spectral point `s`; `u` = Stokes-I
slot indicator over the NquadN streams, indices `1:pol_n:NquadN`;
`v = Diagonal(qp_μN .* wt_μN)·u`):

    S̄    = 2 · vᵀ R⁺⁻ u                    atmosphere's spherical albedo,
                                           viewed from below
    E_dw  = vᵀ J₀⁺ + μ₀ F₀ᴵ e^(−τ_surf/μ₀)  total downward flux at BOA
    t_up  = T⁻⁻ u                          upward transmission of an
                                           isotropic BOA source

    J₀⁻(a) = J₀⁻_atm + [2a·E_dw / (1 − a·S̄)] · t_up             (†)

All m>0 moments are exactly zero for a Lambertian surface layer
(`_zero_surface_layer!` in `lambertian_surface.jl`), so (†) is the WHOLE
albedo dependence; everything else in the TOA sum is captured by the a=0
("black-sky") replay, reused as-is via `rt_run_surface`.
============================================================================
=#

"""
    LambertianClosure{FT}

Analytic closure for the TOA Lambertian-surface reflectance as a function of
albedo `a`, built once from an [`AtmosphereRTCache`](@ref)'s `m=0` blocks.
Evaluating `c(a)` for a new albedo (scalar or per-spectral-point vector) is
an `O(n_vza · pol_n · nSpec)` broadcast — no surface build, no matrix
inversion, no Fourier loop.

# Fields
- `R_black :: Array{FT,3}` — TOA reflectance at `a = 0` (`(n_vza, pol_n,
  nSpec)`), i.e. `rt_run_surface(cache, LambertianSurfaceScalar(0))[1]`.
  Includes every Fourier moment (only `m=0` depends on `a`; `m>0` is
  albedo-independent and folded in here).
- `t_up :: Array{FT,3}` — `T⁻⁻u` restricted to the view-direction rows and
  pre-weighted by the exact `m=0` azimuthal/Stokes factor that
  `postprocessing_vza!` applies (`w₀ · Diagonal([1,1,0,0][1:pol_n])`, i.e.
  `cos(0)=1` on I/Q, `sin(0)=0` on U/V) — so evaluating the closure is a
  plain broadcast against this array, not a matrix multiply.
- `S̄ :: Vector{FT}` — atmosphere's spherical albedo viewed from below,
  `2·vᵀR⁺⁻u`, one value per spectral point.
- `E_dw :: Vector{FT}` — total (diffuse + direct) downward flux at BOA,
  `vᵀJ₀⁺ + μ₀F₀ᴵe^(−τ_surf/μ₀)`.
- `w₀ :: FT` — the `m=0` Fourier weight (`0.5/π`); stored for reference —
  already folded into `t_up`.

See [`lambertian_closure`](@ref) (constructor), the functor
`(c::LambertianClosure)(a)`, [`albedo_jacobian`](@ref), and
[`invert_albedo`](@ref).
"""
struct LambertianClosure{FT}
    R_black :: Array{FT, 3}
    t_up    :: Array{FT, 3}
    S̄       :: Vector{FT}
    E_dw    :: Vector{FT}
    w₀      :: FT
end

"""
    lambertian_closure(cache::AtmosphereRTCache) -> LambertianClosure

Build the analytic Lambertian-albedo closure (proposal §3.1) from an
[`AtmosphereRTCache`](@ref). Works on both `:full` and `:slim` caches (§4)
— only the `m=0` blocks are read, and both cache modes always keep those in
full.

Solar-only: throws `ArgumentError` if the cache carries per-source `J₀`
slots (e.g. thermal emission, which declares a `source_key` and therefore
gets its own composite slot) — the closure term for those sources is a
follow-up (§3 notes this is analogous to "no closure exists in general" for
non-Lambertian BRDFs). Also throws for an active `SurfaceSIF` source: SIF
is injected directly into the surface's legacy `j₀⁻` slot (not a per-source
composite slot, so the check above wouldn't catch it), but it still bounces
between surface and atmosphere in an albedo-dependent way that (†) does not
account for.

GPU note: the cache's device-resident arrays are read exactly once here
(`Array(...)`-converted up front), so the closure itself is built on the
host regardless of the cache's architecture; the returned `LambertianClosure`
is always host-resident, and evaluating it is plain CPU broadcasting.

Like the cache it is built from, the closure describes the scene *as of*
`rt_run_atmosphere` — any subsequent [`update_model!`](@ref) /
`update_aerosol_*!` scene change requires a fresh cache and a fresh closure.
"""
function lambertian_closure(cache::AtmosphereRTCache{FT}) where {FT}
    isempty(cache.J₀_by_src_per_m[1]) || throw(ArgumentError(
        "lambertian_closure: cache carries per-source J₀ slots (e.g. thermal " *
        "emission) — the analytic closure only covers the solar-only path " *
        "today (proposals/surface_split_albedo_sweep.md §3); extending (†) to " *
        "the extra source terms is a follow-up."))
    has_surface_sif(cache.prepared_sources) && throw(ArgumentError(
        "lambertian_closure: cache carries a SurfaceSIF source — SIF is " *
        "injected into the surface's legacy j₀⁻ slot and its albedo-dependent " *
        "reflection between surface and atmosphere is not captured by (†); " *
        "build the cache without a SurfaceSIF source, or use rt_run_surface " *
        "directly for this scene."))

    pol_type    = cache.pol_type
    quad_points = cache.quad_points
    pol_n  = pol_type.n
    NquadN = cache.NquadN
    nSpec  = cache.nSpec
    μ₀     = cache.μ₀

    # Host copies of the m=0 atmosphere-only blocks (read once — cheap even
    # when the cache is device-resident; see the GPU note above).
    R⁺⁻0   = Array(cache.R⁺⁻_per_m[1])
    J₀⁺0   = Array(cache.J₀⁺_per_m[1])
    T⁻⁻0   = Array(cache.T⁻⁻_per_m[1])
    τ_surf = Array(cache.τ_sum_surf)
    F₀ᴵ    = Array(cache.surface_F₀)[1, :]
    qp_μN  = Array(quad_points.qp_μN)
    wt_μN  = Array(quad_points.wt_μN)
    qp_μ   = Array(quad_points.qp_μ)

    # u = Stokes-I slot indicator (indices 1:pol_n:NquadN); v = D·u with
    # D = Diagonal(qp_μN .* wt_μN). Both vanish outside the Stokes-I slots
    # (u is a 0/1 indicator, v inherits its support), so every matvec below
    # only ever touches the `iS`-indexed sub-block — no need to materialize
    # the full length-NquadN u/v vectors.
    iS  = 1:pol_n:NquadN
    qwI = qp_μN[iS] .* wt_μN[iS]
    ones_iS = ones(FT, length(iS))

    S̄    = Vector{FT}(undef, nSpec)
    E_dw = Vector{FT}(undef, nSpec)
    t_up = Matrix{FT}(undef, NquadN, nSpec)
    @inbounds for s in 1:nSpec
        S̄[s]    = 2 * dot(qwI, view(R⁺⁻0, iS, iS, s) * ones_iS)
        E_dw[s] = dot(qwI, view(J₀⁺0, iS, 1, s)) + μ₀ * F₀ᴵ[s] * exp(-τ_surf[s] / μ₀)
        t_up[:, s] .= vec(sum(view(T⁻⁻0, :, iS, s), dims = 2))
    end

    # a=0 ("black-sky") replay — reuses the already-tested `rt_run_surface`
    # path (every Fourier moment, hemispherical diagnostics, ...) instead of
    # re-deriving postprocessing here.
    R_black = rt_run_surface(cache, LambertianSurfaceScalar(zero(FT)))[1]

    # View-row extraction + m=0 Stokes/azimuthal weighting, mirroring
    # `postprocessing_vza!` EXACTLY — same nearest-stream lookup, same index
    # helper, same weight structure: `_precompute_vza_weights` IS the
    # production helper, not a reimplementation.
    n_vza = length(cache.vza)
    w₀    = FT(0.5 / π)
    vza_info0 = _precompute_vza_weights(cache.vza, cache.vaz, qp_μ, pol_type, 0, w₀)
    t_up_weighted = zeros(FT, n_vza, pol_n, nSpec)
    for i in 1:n_vza
        istart, iend, w = vza_info0[i]
        for s in 1:nSpec
            t_up_weighted[i, :, s] .= w * t_up[istart:iend, s]
        end
    end

    return LambertianClosure{FT}(R_black, t_up_weighted, S̄, E_dw, w₀)
end

"""
    _check_albedo_domain(a, S̄)

Validate `a ≥ 0` and `a·S̄ < 1` (the rank-one adding-series factor
`1/(1−a·S̄)` diverges or flips sign otherwise) for either a scalar or a
per-spectral-point albedo vector. Internal helper shared by the functor and
`albedo_jacobian`.
"""
function _check_albedo_domain(a::Real, S̄::AbstractVector)
    a < 0 && throw(ArgumentError("LambertianClosure: albedo must be ≥ 0, got $a"))
    m = maximum(a .* S̄)
    m >= 1 && throw(ArgumentError("LambertianClosure: a·S̄ must be < 1 (adding series must converge); max a·S̄ = $m"))
    return nothing
end

function _check_albedo_domain(a::AbstractVector{<:Real}, S̄::AbstractVector)
    length(a) == length(S̄) || throw(DimensionMismatch(
        "LambertianClosure: albedo vector length $(length(a)) ≠ nSpec $(length(S̄))"))
    any(<(0), a) && throw(ArgumentError("LambertianClosure: albedo must be ≥ 0"))
    m = maximum(a .* S̄)
    m >= 1 && throw(ArgumentError("LambertianClosure: a·S̄ must be < 1 (adding series must converge); max a·S̄ = $m"))
    return nothing
end

"""
    (c::LambertianClosure)(a) -> Array{FT,3}

Evaluate the exact TOA reflectance at Lambertian albedo `a` (proposal §3.1,
Eq. †). `a` may be a scalar (spectrally-flat albedo, broadcasting against
every spectral point like [`LambertianSurfaceScalar`](@ref)) or a
length-`nSpec` vector (like [`LambertianSurfaceSpectrum`](@ref)). Returns an
array shaped `(n_vza, pol_n, nSpec)`, matching `rt_run_surface`'s `R_SFI`.
"""
function (c::LambertianClosure{FT})(a::Union{Real, AbstractVector{<:Real}}) where {FT}
    _check_albedo_domain(a, c.S̄)
    coeff = @. FT(2) * a * c.E_dw / (1 - a * c.S̄)
    return c.R_black .+ reshape(coeff, 1, 1, :) .* c.t_up
end

"""
    albedo_jacobian(c::LambertianClosure, a) -> Array{FT,3}

Exact analytic `∂R_SFI/∂a` at albedo `a` (same shape/broadcast rules as the
functor): `w₀ · 2·E_dw/(1−a·S̄)² · t_up` (proposal §3.2). Differentiating (†)
directly — `R_black` doesn't depend on `a`, so it drops out.
"""
function albedo_jacobian(c::LambertianClosure{FT}, a::Union{Real, AbstractVector{<:Real}}) where {FT}
    _check_albedo_domain(a, c.S̄)
    coeff = @. FT(2) * c.E_dw / (1 - a * c.S̄)^2
    return reshape(coeff, 1, 1, :) .* c.t_up
end

"""
    invert_albedo(c::LambertianClosure, R_meas; i_vza::Int = 1) -> Vector

Closed-form albedo retrieval from a measured TOA Stokes-I reflectance
`R_meas` (same `(n_vza, pol_n, nSpec)` shape as the functor's output — e.g.
from a real instrument, or from `rt_run`/`rt_run_surface`) at view index
`i_vza`. No iteration, no linearized RT run: invert (†) for `a` per spectral
point,

    y = (R_meas − R_black) / t_upᴵ;  a = y / (2·E_dw + S̄·y)

where `t_upᴵ = c.t_up[i_vza, 1, :]` already carries the `w₀`/Stokes-I
weighting.

When `|t_upᴵ[s]|` or `|E_dw[s]|` is at the machine-epsilon floor (the
Stokes-I upward transmission or the downward flux effectively vanishes — an
ill-conditioned spectral point, e.g. a fully opaque path), the inversion is
undefined and `a[s] = NaN` (rather than an arbitrarily large or garbage
value from dividing by ~0).
"""
function invert_albedo(c::LambertianClosure{FT}, R_meas::AbstractArray{<:Real, 3};
                       i_vza::Int = 1) where {FT}
    n_vza = size(c.R_black, 1)
    nSpec = length(c.S̄)
    1 <= i_vza <= n_vza || throw(ArgumentError("invert_albedo: i_vza=$i_vza out of range 1:$n_vza"))
    # Validate every axis the (@inbounds) loop indexes: with bounds checking
    # off, a too-short spectral axis would read garbage rather than error.
    size(R_meas, 1) >= i_vza || throw(DimensionMismatch(
        "invert_albedo: R_meas has $(size(R_meas, 1)) VZA rows, need at least $i_vza"))
    size(R_meas, 2) >= 1 || throw(DimensionMismatch(
        "invert_albedo: R_meas has no Stokes/polarization axis (size(R_meas, 2) == 0)"))
    size(R_meas, 3) == nSpec || throw(DimensionMismatch(
        "invert_albedo: R_meas has $(size(R_meas, 3)) spectral points, closure has $nSpec"))

    a = Vector{FT}(undef, nSpec)
    @inbounds for s in 1:nSpec
        t = c.t_up[i_vza, 1, s]
        if abs(t) <= eps(FT) || abs(c.E_dw[s]) <= eps(FT)
            a[s] = FT(NaN)
            continue
        end
        y = (R_meas[i_vza, 1, s] - c.R_black[i_vza, 1, s]) / t
        a[s] = y / (2 * c.E_dw[s] + c.S̄[s] * y)
    end
    return a
end
