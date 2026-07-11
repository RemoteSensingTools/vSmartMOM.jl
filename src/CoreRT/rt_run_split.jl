# =============================================================================
# rt_run_atmosphere / rt_run_surface / rt_run_multi_surface
#
# Atmosphere/surface split of `rt_run`. In the Fourier loop (`rt_run.jl`),
# the surface enters the solve at exactly one point per moment `m`: after
# the z-loop (elemental → doubling → interaction over `Nz` layers) settles,
# a surface `AddedLayer` is built and fused into the composite with one
# final `interaction!`. Everything before that — the expensive part,
# `Nz` layers × (elemental + doubling + per-layer interaction), per moment,
# per spectral point — is completely surface-independent. This file lets
# callers pay that cost once (`rt_run_atmosphere`) and replay only the
# surface phase per BRDF (`rt_run_surface`), which is
# `O(NquadN^3 · nSpec)` per moment instead of `O(Nz · NquadN^3 · nSpec)`.
#
# Ported from the `gchp-io` branch prototype (see
# `proposals/surface_split_albedo_sweep.md` §1-2) and re-mirrored against
# the current `rt_run.jl` surface block. The two blocks must stay in
# lockstep: when `rt_run.jl`'s surface block changes, update the mirrored
# calls in `rt_run_surface` below (same call order, same kwargs).
#
# Scope (this PR):
#   * `noRS` only. `rt_run_atmosphere` delegates to the public
#     `rt_run(model; ...)` entry point, which always constructs `noRS`
#     internally — there is no way to reach this code path with a Raman
#     `RS_type`, so no separate runtime check is needed. Raman 4-D `ie*`
#     snapshots are deferred to a follow-up.
#   * Single `i_band` (an `Integer`, matching `rt_run`'s own kwarg — no
#     vector form here). Multi-band multi-surface is deferred.
#   * Lambertian / Cox-Munk / RossLi / RPV surfaces. `CanopySurface` is
#     rejected with `ArgumentError` in both `rt_run_atmosphere` (the
#     model's own surface, or any `target_brdfs` entry) and
#     `rt_run_surface` (the replayed `brdf`) — the canopy preload
#     (`_init_canopy_cache!` / `_compute_canopy_atm_tau!`, run once before
#     `rt_run`'s m loop) would need to move into a reusable helper first.
# =============================================================================

"""
    rt_run_atmosphere(model::RTModel;
                      target_brdfs = nothing,
                      i_band::Integer = 1,
                      sources = nothing,
                      cache_mode::Symbol = :auto) -> AtmosphereRTCache

Run the atmosphere phase of [`rt_run`](@ref) once and cache the
per-Fourier-moment composite-layer snapshot. The returned
[`AtmosphereRTCache`](@ref) can be consumed by [`rt_run_surface`](@ref)
repeatedly to swap surfaces over the same atmosphere without paying the
m-loop layer-accumulation cost again.

# Keyword arguments
- `target_brdfs` — surfaces the caller intends to swap in via
  `rt_run_surface`. The cache's Fourier loop bound is sized for
  `max(model's own m_max, maximum(component_m_max(b, ctx) for b in
  target_brdfs))`, so e.g. a Lambertian model can still serve a Cox-Munk
  replay. Defaults to `model.surfaces` (least surprising — the cache
  always covers the model's own surface at minimum).
- `i_band::Integer = 1` — spectral band index.
- `sources` — source set, as in `rt_run`.
- `cache_mode::Symbol = :auto` — `:full` (all six composite blocks at every
  Fourier moment), `:slim` (proposal §4: six blocks only at `m=0`; `J₀⁺`/
  `J₀⁻`/per-source slots only at `m>0` — valid only when every `target_brdfs`
  entry is Lambertian-family), or `:auto` (picks `:slim` iff every entry in
  `brdfs` — `target_brdfs`, or `model.surfaces` when left at its default —
  has `component_m_max(b, ctx) == 0`, else `:full`). See the memory-footprint
  note on [`AtmosphereRTCache`](@ref).

Only `noRS` is supported (see the module-level scope note above);
`CanopySurface` in `target_brdfs` (or the model's own surface, when
`target_brdfs` is left at its default) raises `ArgumentError`.

!!! warning "Caches do not survive scene updates"
    The cache snapshots the model's optics *at build time*. Any subsequent
    [`update_model!`](@ref), [`update_aerosol_loading!`](@ref) /
    [`update_aerosol_microphysics!`](@ref), or direct mutation of
    `model.optics` silently invalidates it (and any
    [`LambertianClosure`](@ref) derived from it) — the replay would keep
    answering for the *old* scene. Per-scenario order must always be:
    update the scene → `rt_run_atmosphere` → sweep surfaces. See the
    [Fast Re-runs & Batch Processing](@ref) guide.
"""
function rt_run_atmosphere(model;
                           target_brdfs = nothing,
                           i_band::Integer = 1,
                           sources::Union{Nothing, AbstractSource} = nothing,
                           cache_mode::Symbol = :auto)
    cache_mode in (:auto, :full, :slim) || throw(ArgumentError(
        "rt_run_atmosphere: cache_mode must be :auto, :full, or :slim (got $(repr(cache_mode)))."))

    iBand = [Int(i_band)]
    brdfs = target_brdfs === nothing ? model.surfaces : collect(target_brdfs)
    for b in brdfs
        if b isa CanopySurface
            throw(ArgumentError("rt_run_atmosphere does not support CanopySurface; the canopy preload (_init_canopy_cache!/_compute_canopy_atm_tau!) needs to move into a reusable helper first (deferred)."))
        end
    end

    # Cache-wide Fourier loop bound: max over (the model's own per-band
    # m_max, the target BRDFs' declared Fourier support). `model.solver.l_trunc`
    # is used as the `user_l_cap` proxy for `component_m_max` — exact for
    # aerosol-free scenes (this PR's test scope); with aerosols present the
    # true per-band cap can be tighter (see `_derive_m_max_bands_via_traits`
    # in `tools/model_from_parameters.jl`), so a widened cache could
    # over-allocate relative to what's strictly needed. `stream_l_cap` is
    # computed exactly from the quadrature (`2·Nstreams - 1`).
    user_l_cap = model.solver.l_trunc
    stream_l_cap = 2 * model.quad_points.Nstreams - 1
    ctx = (; user_l_cap, stream_l_cap, truncation = nothing, m_max_override = nothing)
    m_max_cache = max(model.solver.m_max_bands[iBand[1]],
                      maximum(component_m_max(b, ctx) for b in brdfs))

    # :auto -> :slim iff every target BRDF is Lambertian-family
    # (component_m_max == 0 — the only case where the m>0 surface layer is
    # provably zero, see proposal §4). Otherwise :full.
    resolved_cache_mode = cache_mode == :auto ?
        (all(==(0), (component_m_max(b, ctx) for b in brdfs)) ? :slim : :full) :
        cache_mode
    slim = resolved_cache_mode == :slim

    # Collection slots — filled by the callback below, one push per Fourier
    # moment (m_max_cache + 1 calls total; a handful for typical configs).
    # Held as `Vector{Any}` and narrowed to concrete element types once the
    # run completes — simpler than threading the concrete array type through
    # the callback closure, and the narrowing cost is negligible next to one
    # `rt_run` call.
    R⁻⁺ = Vector{Any}()
    R⁺⁻ = Vector{Any}()
    T⁺⁺ = Vector{Any}()
    T⁻⁻ = Vector{Any}()
    J₀⁺ = Vector{Any}()
    J₀⁻ = Vector{Any}()
    J0_by_src = NamedTuple[]
    sci_surf  = AbstractScatteringInterface[]
    FT_pre  = float_type(model)
    weights = FT_pre[]

    captured_ctx = Ref{Union{Nothing, NamedTuple}}(nothing)

    function cb(state)
        cl = state.composite_layer
        if !slim || state.m == 0
            push!(R⁻⁺, deepcopy(cl.R⁻⁺))
            push!(R⁺⁻, deepcopy(cl.R⁺⁻))
            push!(T⁺⁺, deepcopy(cl.T⁺⁺))
            push!(T⁻⁻, deepcopy(cl.T⁻⁻))
        else
            # Slim m>0: a Lambertian-family surface layer is exactly zero
            # here, so the final interaction is the identity on J₀⁺/J₀⁻ and
            # rt_run_surface never reads R/T at m>0 (see its slim branch) —
            # store zero-size placeholders of the SAME concrete array type
            # (so `AT3[R⁻⁺...]` below still narrows to one concrete type)
            # instead of a full (NquadN,NquadN,nSpec) deepcopy.
            push!(R⁻⁺, similar(cl.R⁻⁺, 0, 0, 0))
            push!(R⁺⁻, similar(cl.R⁺⁻, 0, 0, 0))
            push!(T⁺⁺, similar(cl.T⁺⁺, 0, 0, 0))
            push!(T⁻⁻, similar(cl.T⁻⁻, 0, 0, 0))
        end
        push!(J₀⁺, deepcopy(cl.J₀⁺))
        push!(J₀⁻, deepcopy(cl.J₀⁻))
        push!(J0_by_src, _deepcopy_J0_by_src(cl.J₀_by_src))
        push!(sci_surf, state.scattering_interface_surf)
        push!(weights, FT_pre(state.weight))
        # Only `m` / `weight` / `composite_layer` / `scattering_interface_surf`
        # vary across calls within one `rt_run` invocation; the surroundings
        # (pol_type, quad_points, arch, ...) are stable, so the last-fired
        # context is representative of all of them.
        captured_ctx[] = state
        return nothing
    end

    rt_run(model; i_band = i_band, sources = sources,
                  stop_after_atmosphere = true,
                  atm_snapshot_callback = cb,
                  m_max_override = m_max_cache)

    captured_ctx[] === nothing &&
        error("rt_run_atmosphere: callback never fired — the atmosphere phase produced no Fourier moments.")
    c = captured_ctx[]

    FT  = float_type(model)
    AT3 = typeof(R⁻⁺[1])
    τ_surf_arr = deepcopy(c.τ_sum_surf)
    F₀_arr     = deepcopy(c.surface_F₀)
    AT1 = typeof(τ_surf_arr)
    AT2 = typeof(F₀_arr)
    # The surface-bottom scattering interface is uniform across m for any
    # single scene, so the eltype narrows to one concrete
    # `<: AbstractScatteringInterface`.
    SCI = typeof(sci_surf[1])

    return AtmosphereRTCache{FT, AT3, AT1, AT2, c.arr_type,
                             typeof(c.pol_type), typeof(c.quad_points),
                             typeof(c.arch), typeof(c.RS_type),
                             typeof(c.prepared_sources), typeof(c.I_static),
                             SCI}(
        c.pol_type, c.quad_points, c.iμ₀, FT(c.μ₀),
        FT.(collect(c.vza)), FT.(collect(c.vaz)),
        c.nSpec, m_max_cache, user_l_cap, c.NquadN, resolved_cache_mode,
        iBand, weights,
        c.arr_type, c.arch, c.SFI, c.RS_type, c.prepared_sources, c.I_static,
        AT3[R⁻⁺...], AT3[R⁺⁻...],
        AT3[T⁺⁺...], AT3[T⁻⁻...],
        AT3[J₀⁺...], AT3[J₀⁻...],
        J0_by_src, SCI[sci_surf...],
        τ_surf_arr, F₀_arr)
end

# Per-source J₀ NamedTuple deepcopy helper. Each slot is a
# `CompositeSourceSlot{FT}` with `J₀⁺` / `J₀⁻` 3-D arrays.
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
                   verbose::Bool = false)
        -> (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)

Replay the surface phase of [`rt_run`](@ref) against a cached atmosphere.
Returns the same tuple shape as `rt_run(model)` under SFI (`cache.SFI` is
always `true` — `rt_run` hardcodes SFI).

Mirrors `rt_run.jl`'s surface block (`create_surface_layer!` →
`surface_source_contribute!` → `interaction!` → `interaction_hdrf!` →
`postprocessing_vza!`/`_hdrf!`) call-for-call, so results are bit-exact
with a monolithic `rt_run` call using the same `brdf` and the same Fourier
loop bound.

Throws `ArgumentError` for `CanopySurface` (not yet supported — see the
module-level scope note in this file), and when `component_m_max(brdf,
ctx)` — evaluated against the *model's* Fourier-support cap
(`cache.user_l_cap`), not the cache's own width (`cache.m_max`); see
[`AtmosphereRTCache`](@ref) — exceeds `cache.m_max`: rebuild the cache with
this BRDF (or one of the same family) in `target_brdfs`. Also throws
`ArgumentError` when `cache.cache_mode == :slim` and `brdf` is not
Lambertian-family (`component_m_max(brdf, ctx) != 0`) — the cache's m>0
blocks were discarded on the assumption that only a Lambertian-family
surface would ever be replayed; rebuild with `cache_mode = :full`.
"""
function rt_run_surface(cache::AtmosphereRTCache{FT}, brdf;
                        verbose::Bool = false) where {FT}
    if brdf isa CanopySurface
        throw(ArgumentError("rt_run_surface does not support CanopySurface (deferred — see the module-level scope note in rt_run_split.jl)."))
    end

    # Guard: does this BRDF fit within the cache's Fourier loop bound?
    ctx_brdf = (; user_l_cap = cache.user_l_cap, stream_l_cap = cache.user_l_cap,
                  truncation = nothing, m_max_override = nothing)
    m_max_needed = component_m_max(brdf, ctx_brdf)
    if m_max_needed > cache.m_max
        throw(ArgumentError("rt_run_surface: BRDF requires m_max=$(m_max_needed) but cache was built with m_max=$(cache.m_max). Pass this BRDF (or one of the same family) to `target_brdfs` on `rt_run_atmosphere`."))
    end

    # Guard: a :slim cache only kept m>0 atmosphere blocks on the assumption
    # that the replayed surface is Lambertian-family (m>0 surface layer
    # exactly zero, see proposal §4). A non-Lambertian BRDF needs the
    # discarded R/T blocks at m>0 — reject rather than silently replaying
    # against stale/empty data.
    slim = cache.cache_mode === :slim
    if slim && m_max_needed != 0
        throw(ArgumentError("rt_run_surface: cache was built with cache_mode=:slim (valid only for Lambertian-family surfaces — component_m_max==0); $(typeof(brdf)) needs m_max=$(m_max_needed) support at m>0, which this cache discarded. Rebuild with cache_mode=:full (or cache_mode=:auto with this BRDF in `target_brdfs`)."))
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

    # Fresh composite_layer + added_layer_surface, repopulated from the
    # snapshot every m below (mirrors rt_run.jl's single allocation before
    # the m loop — here it happens once per `rt_run_surface` call instead
    # of once per full `rt_run`).
    composite_layer     = make_composite_layer(RS_type, FT, arr_type, dims, nSpec;
                                               prepared_sources = prepared_sources)
    added_layer_surface = make_added_layer(RS_type, FT, arr_type, dims, nSpec;
                                           prepared_sources = prepared_sources)

    # Output arrays — same shapes as rt_run.
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

    # hdr_J₀⁻ is overwritten by `interaction_hdrf!` each moment; pre-allocate
    # once (mirrors rt_run.jl's pre-loop hoist — avoids one `similar()` per m).
    hdr_J₀⁻ = similar(composite_layer.J₀⁻)

    for m in 0:cache.m_max
        weight = cache.weight_per_m[m + 1]

        if slim && m > 0
            # Slim replay at m>0 (proposal §4): the cache discarded R/T here
            # because a Lambertian-family surface layer is EXACTLY zero at
            # m>0 (`_zero_surface_layer!`), so `create_surface_layer!` →
            # `surface_source_contribute!` → `interaction!` →
            # `interaction_hdrf!` reduce to an identity on J₀⁺/J₀⁻ and an
            # exact zero on `hdr_J₀⁻` (`hdr_J₀⁻ = r⁻⁺⊠J₀⁺ + j₀⁻ = 0⊠J₀⁺ + 0`)
            # — verified bit-exact against the full replay in
            # test/test_lambertian_closure.jl. Restore J₀⁺/J₀⁻ (+ per-source
            # slots) directly and skip straight to postprocessing.
            @timeit "Restore snapshot (slim)" begin
                copyto!(composite_layer.J₀⁺, cache.J₀⁺_per_m[m + 1])
                copyto!(composite_layer.J₀⁻, cache.J₀⁻_per_m[m + 1])
                for k in keys(composite_layer.J₀_by_src)
                    slot_src = cache.J₀_by_src_per_m[m + 1][k]
                    slot_dst = composite_layer.J₀_by_src[k]
                    copyto!(slot_dst.J₀⁺, slot_src.J₀⁺)
                    copyto!(slot_dst.J₀⁻, slot_src.J₀⁻)
                end
            end
            hdr_J₀⁻ .= zero(FT)
        else
            # Restore the atmosphere-only snapshot into composite_layer. The
            # cached blocks are already device-resident copies (built via
            # `deepcopy` in `rt_run_atmosphere`), so no `arr_type(...)` re-wrap
            # is needed here — `copyto!` handles same-device same-shape arrays
            # directly.
            @timeit "Restore snapshot" begin
                copyto!(composite_layer.R⁻⁺, cache.R⁻⁺_per_m[m + 1])
                copyto!(composite_layer.R⁺⁻, cache.R⁺⁻_per_m[m + 1])
                copyto!(composite_layer.T⁺⁺, cache.T⁺⁺_per_m[m + 1])
                copyto!(composite_layer.T⁻⁻, cache.T⁻⁻_per_m[m + 1])
                copyto!(composite_layer.J₀⁺, cache.J₀⁺_per_m[m + 1])
                copyto!(composite_layer.J₀⁻, cache.J₀⁻_per_m[m + 1])
                for k in keys(composite_layer.J₀_by_src)
                    slot_src = cache.J₀_by_src_per_m[m + 1][k]
                    slot_dst = composite_layer.J₀_by_src[k]
                    copyto!(slot_dst.J₀⁺, slot_src.J₀⁺)
                    copyto!(slot_dst.J₀⁻, slot_src.J₀⁻)
                end
            end

            # Surface step — mirrors rt_run.jl's surface block (non-canopy
            # branch) call-for-call. NOTE: current `rt_run` does NOT call
            # `inject_surface_SIF!` here (legacy `RS_type.SIF₀` is intentionally
            # ignored; SIF flows through `surface_source_contribute!` via the
            # `SurfaceSIF` source) — do not add it back; see rt_run.jl's own
            # comment at its surface block for the rationale.
            @timeit "Create Surface" create_surface_layer!(brdf, added_layer_surface,
                                SFI, m, pol_type, quad_points,
                                cache.τ_sum_surf, arch;
                                F₀ = cache.surface_F₀)

            surface_source_contribute!(prepared_sources, brdf, added_layer_surface,
                                       m, pol_type, arch)

            @timeit "interaction" interaction!(RS_type, cache.sci_surf_per_m[m + 1], SFI,
                                composite_layer, added_layer_surface, I_static;
                                workspace = _interaction_ws)

            @timeit "interaction_HDRF" interaction_hdrf!(SFI, composite_layer,
                                added_layer_surface, m, pol_type, quad_points,
                                hdr_J₀⁻, bhr_uw, bhr_dw)
        end

        @timeit "Postprocessing VZA" postprocessing_vza!(RS_type, iμ₀, pol_type,
                            composite_layer, vza, qp_μ, m, vaz, μ₀, weight,
                            nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)
        @timeit "Postprocessing HDRF" postprocessing_vza_hdrf!(RS_type, iμ₀, pol_type,
                            hdr_J₀⁻, vza, qp_μ, m, vaz, μ₀, weight, nSpec, hdr)
    end

    # Single-scattering correction for Cox-Munk specular hotspot (TMS) —
    # replayed per BRDF, using the cache's τ_sum_surf / m_max (mirrors
    # rt_run.jl's post-loop block; see its LIMITATION note re: non-unit F₀).
    if brdf isa CoxMunkSurface && SFI
        @timeit "SS Correction" apply_ss_correction!(R_SFI, brdf, pol_type, vza, vaz, μ₀,
                             Array(cache.τ_sum_surf), cache.m_max, nSpec)
    end

    verbose && print_timer()

    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr,
                  bhr_uw[1, :], bhr_dw[1, :]) : (R, T)
end

# =============================================================================

"""
    rt_run_multi_surface(model::RTModel, brdfs::AbstractVector;
                         i_band::Integer = 1, sources = nothing)

Sweep a set of BRDFs over a single atmosphere. The atmosphere is built once
via [`rt_run_atmosphere`](@ref) (with `target_brdfs = brdfs`, so the cache
is sized for every BRDF in the sweep); each BRDF is then run through
[`rt_run_surface`](@ref).

Returns a `Vector` of `rt_run`-style tuples, one per BRDF in `brdfs`, in the
same order.

# Acceptance
`rt_run_multi_surface(model, [model.surfaces[1]])[1]` is bit-exact with
`rt_run(model)` — same kernels, same values, no separate code path for the
"just replay the model's own surface" case.
"""
function rt_run_multi_surface(model, brdfs::AbstractVector;
                              i_band::Integer = 1,
                              sources::Union{Nothing, AbstractSource} = nothing)
    cache = rt_run_atmosphere(model; target_brdfs = brdfs,
                              i_band = i_band, sources = sources)
    return [rt_run_surface(cache, brdf) for brdf in brdfs]
end
