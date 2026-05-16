# =============================================================================
# rt_run_atmosphere / rt_run_surface / rt_run_multi_surface
#
# Atmosphere/surface split of `rt_run`. Builds an `AtmosphereRTCache` once
# per scene and replays the surface phase per BRDF — the m-loop layer
# accumulation runs once across many surfaces.
#
# Scope (this PR):
#   * `noRS` only — Raman 4-D `ie*` snapshots are deferred.
#   * Single iBand (the first one). Multi-band multi-surface is deferred.
#   * Lambertian / CoxMunk / RossLi / RPV surfaces. Canopy is rejected
#     with a clear ArgumentError (the canopy preload would need to move
#     into `_prepare_canopy_for_run!`; deferred to a follow-up).
# =============================================================================

"""
    rt_run_atmosphere(model::RTModel;
                      target_brdfs = nothing,
                      i_band = 1,
                      sources = nothing) -> AtmosphereRTCache

Run the atmosphere phase of `rt_run` once and cache the per-Fourier-moment
composite-layer snapshot. The returned [`AtmosphereRTCache`](@ref) can be
consumed by [`rt_run_surface`](@ref) repeatedly to swap surfaces over the
same atmosphere without paying the m-loop layer accumulation cost.

# Keyword arguments
- `target_brdfs :: Vector{<:AbstractSurfaceType}` — surfaces the caller
  intends to swap in via `rt_run_surface`. The cache m-max is sized for
  `maximum(component_m_max(b, ctx) for b in target_brdfs)`. Defaults to
  `model.surfaces` (least surprising), so the cache always covers the
  model's own surface at minimum.
- `i_band :: Integer = 1` — spectral band index.
- `sources :: Union{Nothing, AbstractSource}` — source set, as in `rt_run`.

Only `noRS` is supported in this PR; the Raman variants raise
`ArgumentError`. Canopy surfaces are also rejected (lift the canopy
preload first — deferred).
"""
function rt_run_atmosphere(model;
                           target_brdfs = nothing,
                           i_band::Integer = 1,
                           sources::Union{Nothing, AbstractSource} = nothing)
    iBand = [Int(i_band)]
    brdfs = target_brdfs === nothing ? model.surfaces : Vector(target_brdfs)
    for b in brdfs
        if b isa CanopySurface
            throw(ArgumentError("rt_run_atmosphere does not yet support CanopySurface; canopy preload needs to be lifted into a helper first (deferred)."))
        end
    end

    # Compute the cache-wide m_max: max over (current per-band, target BRDFs).
    # `l_trunc` (Legendre truncation order from the SolverConfig) is the
    # upper-bound proxy for `user_l_cap` — Cox-Munk / RPV / RossLi return
    # this value, Lambertian returns 0.
    user_l_cap = model.solver.l_trunc
    ctx = (; user_l_cap, stream_l_cap = user_l_cap,
             truncation = nothing, m_max_override = nothing)
    m_max_cache = max(model.solver.m_max_bands[iBand[1]],
                      maximum(component_m_max(b, ctx) for b in brdfs))

    # Collection slots — filled by the callback. We hold concrete arrays
    # (one per m); their element types are pinned from the first push
    # (the callback's `cl.R⁻⁺` etc. carry the architecture-specific
    # concrete 3-D array type). The Fourier-loop bound is widened
    # *non-mutatively* via the `m_max_override` kwarg passed to `rt_run`
    # below — `model.solver.m_max_bands` is never touched.
    FT_pre = float_type(model)
    R⁻⁺ = Vector{Any}()           # narrowed via `[FT_arr...]` after fill
    R⁺⁻ = Vector{Any}()
    T⁺⁺ = Vector{Any}()
    T⁻⁻ = Vector{Any}()
    J₀⁺ = Vector{Any}()
    J₀⁻ = Vector{Any}()
    J0_by_src = NamedTuple[]
    sci_surf = AbstractScatteringInterface[]
    weights = FT_pre[]

    captured_ctx = Ref{Union{Nothing, NamedTuple}}(nothing)

    function cb(ctx)
        cl = ctx.composite_layer
        push!(R⁻⁺, deepcopy(cl.R⁻⁺))
        push!(R⁺⁻, deepcopy(cl.R⁺⁻))
        push!(T⁺⁺, deepcopy(cl.T⁺⁺))
        push!(T⁻⁻, deepcopy(cl.T⁻⁻))
        push!(J₀⁺, deepcopy(cl.J₀⁺))
        push!(J₀⁻, deepcopy(cl.J₀⁻))
        push!(J0_by_src, _deepcopy_J0_by_src(cl.J₀_by_src))
        push!(sci_surf, ctx.scattering_interface_surf)
        push!(weights, FT_pre(ctx.weight))
        # Last-fired context wins — only `m`/`weight`/composite_layer vary
        # across calls; the immutable surroundings are stable.
        captured_ctx[] = ctx
        return nothing
    end

    rt_run(model; i_band = i_band, sources = sources,
                  stop_after_atmosphere = true,
                  atm_snapshot_callback = cb,
                  m_max_override = m_max_cache)

    captured_ctx[] === nothing &&
        error("rt_run_atmosphere: callback never fired — atmosphere phase produced no Fourier moments.")
    c = captured_ctx[]

    FT  = float_type(model)
    AT3 = typeof(R⁻⁺[1])
    τ_surf_arr = deepcopy(c.τ_sum_surf)
    AT1 = typeof(τ_surf_arr)
    # Narrow the scattering-interface vector. The surface-bottom interface
    # is uniform across m for any single scene; the eltype reduces to one
    # concrete `<: AbstractScatteringInterface`.
    SCI = typeof(sci_surf[1])

    return AtmosphereRTCache{FT, AT3, AT1, c.arr_type,
                             typeof(c.pol_type), typeof(c.quad_points),
                             typeof(c.arch), typeof(c.RS_type),
                             typeof(c.prepared_sources), typeof(c.I_static),
                             SCI}(
        c.pol_type, c.quad_points, c.iμ₀, FT(c.μ₀),
        FT.(c.vza), FT.(c.vaz),
        c.nSpec, m_max_cache, c.NquadN, iBand,
        weights,
        c.arr_type, c.arch, c.SFI, c.RS_type, c.prepared_sources, c.I_static,
        AT3[R⁻⁺...], AT3[R⁺⁻...],
        AT3[T⁺⁺...], AT3[T⁻⁻...],
        AT3[J₀⁺...], AT3[J₀⁻...],
        J0_by_src, SCI[sci_surf...], τ_surf_arr)
end

# Per-source J₀ NamedTuple deepcopy helper. Each slot is a
# `CompositeSourceSlot{FT}` with `J₀⁺` / `J₀⁻` 3D arrays.
function _deepcopy_J0_by_src(nt::NamedTuple)
    isempty(nt) && return nt
    keys_ = keys(nt)
    vals = map(keys_) do k
        slot = nt[k]
        CompositeSourceSlot{eltype(slot.J₀⁺)}(
            J₀⁺ = deepcopy(slot.J₀⁺),
            J₀⁻ = deepcopy(slot.J₀⁻),
        )
    end
    return NamedTuple{keys_}(vals)
end

# =============================================================================

"""
    rt_run_surface(cache::AtmosphereRTCache, brdf;
                   verbose = false) -> (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)

Replay the surface phase of `rt_run` against a cached atmosphere. Returns
the same tuple shape as `rt_run(model)` under SFI.

Throws `ArgumentError` when `component_m_max(brdf, ctx) > cache.m_max` —
rebuild the cache with the BRDF in `target_brdfs`.
"""
function rt_run_surface(cache::AtmosphereRTCache{FT}, brdf;
                        verbose::Bool = false) where {FT}
    if brdf isa CanopySurface
        throw(ArgumentError("rt_run_surface does not yet support CanopySurface (deferred)."))
    end

    # Guard: BRDF must fit within the cache's m_max budget
    ctx_brdf = (; user_l_cap = cache.m_max,
                  stream_l_cap = cache.m_max,
                  truncation = nothing,
                  m_max_override = nothing)
    m_max_needed = component_m_max(brdf, ctx_brdf)
    if m_max_needed > cache.m_max
        throw(ArgumentError("rt_run_surface: BRDF requires m_max=$(m_max_needed) but cache was built with m_max=$(cache.m_max). Pass this BRDF (or one of the same family) to `target_brdfs` on `rt_run_atmosphere`."))
    end

    pol_type    = cache.pol_type
    quad_points = cache.quad_points
    iμ₀         = cache.iμ₀
    μ₀          = cache.μ₀
    vza         = cache.vza
    vaz         = cache.vaz
    nSpec       = cache.nSpec
    NquadN      = cache.NquadN
    arr_type    = cache.arr_type
    arch        = cache.arch
    SFI         = cache.SFI
    RS_type     = cache.RS_type
    I_static    = cache.I_static
    prepared_sources = cache.prepared_sources
    qp_μ        = quad_points.qp_μ

    dims = (NquadN, NquadN)

    # Fresh composite_layer + added_layer_surface that we'll repopulate per m.
    composite_layer     = make_composite_layer(RS_type, FT, arr_type, dims, nSpec;
                                               prepared_sources=prepared_sources)
    added_layer_surface = make_added_layer(RS_type, FT, arr_type, dims, nSpec;
                                           prepared_sources=prepared_sources)

    # Output arrays (same shape as rt_run)
    R       = zeros(FT, length(vza), pol_type.n, nSpec)
    T       = zeros(FT, length(vza), pol_type.n, nSpec)
    R_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    T_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    ieR_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    ieT_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    hdr     = zeros(FT, length(vza), pol_type.n, nSpec)
    bhr_dw  = zeros(FT, pol_type.n, nSpec)
    bhr_uw  = zeros(FT, pol_type.n, nSpec)

    _interaction_ws = _interaction_workspace(RS_type, composite_layer, added_layer_surface, arch)

    # Pre-allocated scratch for the HDRF source term — reused across m
    # so the `m_max+1` allocations the previous code did per `rt_run_surface`
    # call become one allocation (julia-style review finding).
    hdr_J₀⁻ = similar(composite_layer.J₀⁻)

    for m in 0:cache.m_max
        weight = cache.weight_per_m[m + 1]

        # Restore atmosphere snapshot into composite_layer
        copyto!(composite_layer.R⁻⁺, cache.R⁻⁺_per_m[m + 1])
        copyto!(composite_layer.R⁺⁻, cache.R⁺⁻_per_m[m + 1])
        copyto!(composite_layer.T⁺⁺, cache.T⁺⁺_per_m[m + 1])
        copyto!(composite_layer.T⁻⁻, cache.T⁻⁻_per_m[m + 1])
        copyto!(composite_layer.J₀⁺, cache.J₀⁺_per_m[m + 1])
        copyto!(composite_layer.J₀⁻, cache.J₀⁻_per_m[m + 1])
        # Per-source slots
        for k in keys(composite_layer.J₀_by_src)
            slot_src = cache.J₀_by_src_per_m[m + 1][k]
            slot_dst = composite_layer.J₀_by_src[k]
            copyto!(slot_dst.J₀⁺, slot_src.J₀⁺)
            copyto!(slot_dst.J₀⁻, slot_src.J₀⁻)
        end

        # Surface step (mirrors rt_run.jl:401-473 exactly for non-canopy)
        create_surface_layer!(brdf, added_layer_surface, SFI, m, pol_type,
                              quad_points, arr_type(cache.τ_sum_surf), arch)

        inject_surface_SIF!(brdf, added_layer_surface, m, pol_type,
                            _sif_source(RS_type), arch)
        surface_source_contribute!(prepared_sources, brdf,
                                   added_layer_surface, m, pol_type, arch)

        interaction!(RS_type, cache.sci_surf_per_m[m + 1], SFI,
                     composite_layer, added_layer_surface, I_static;
                     workspace = _interaction_ws)

        # `hdr_J₀⁻` is overwritten by `interaction_hdrf!` each iteration;
        # allocated once above.
        interaction_hdrf!(SFI, composite_layer, added_layer_surface, m,
                          pol_type, quad_points, hdr_J₀⁻, bhr_uw, bhr_dw)

        postprocessing_vza!(RS_type, iμ₀, pol_type, composite_layer,
                            vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI,
                            R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
        postprocessing_vza_hdrf!(RS_type, iμ₀, pol_type, hdr_J₀⁻,
                                 vza, qp_μ, m, vaz, μ₀, weight, nSpec, hdr)
    end

    # Cox-Munk SS correction (per-BRDF, not per-cache)
    if brdf isa CoxMunkSurface && SFI
        apply_ss_correction!(R_SFI, brdf, pol_type, vza, vaz, μ₀,
                             Array(cache.τ_sum_surf),
                             cache.m_max, nSpec)
    end

    verbose && print_timer()

    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr,
                  bhr_uw[1, :], bhr_dw[1, :]) : (R, T)
end

# =============================================================================

"""
    rt_run_multi_surface(model::RTModel, brdfs::Vector{<:AbstractSurfaceType};
                         i_band = 1, sources = nothing)

Sweep a set of BRDFs over a single atmosphere. The atmosphere is built
once via [`rt_run_atmosphere`](@ref); each BRDF is then run through
[`rt_run_surface`](@ref).

Returns a `Vector` of `rt_run`-style tuples, one per BRDF in `brdfs`.

# Acceptance
`rt_run_multi_surface(model, [model.surfaces[1]])[1] ≡ rt_run(model)`
(bit-exact for Lambertian surfaces).
"""
function rt_run_multi_surface(model, brdfs::AbstractVector;
                              i_band::Integer = 1,
                              sources::Union{Nothing, AbstractSource} = nothing)
    cache = rt_run_atmosphere(model;
                              target_brdfs = brdfs,
                              i_band = i_band,
                              sources = sources)
    return [rt_run_surface(cache, brdf) for brdf in brdfs]
end
