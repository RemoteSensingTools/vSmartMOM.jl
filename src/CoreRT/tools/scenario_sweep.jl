# =============================================================================
# Scenario sweep API (gchp-io PR — Phase D)
#
# Declarative description of the cross product (SZA × view-pair × BRDF) to
# evaluate for one scene, plus a driver that executes it efficiently:
#
#   * VZA + VAZ are PAIRED (one viewing direction per index), matching
#     `postprocessing_vza!` (postprocessing_vza.jl:23). The constructor
#     accepts paired or crossed forms and explicitly `vec`s the latter.
#   * Surfaces are free per-SZA: `rt_run_multi_surface` builds the
#     atmosphere once per SZA and replays the surface step per BRDF.
#   * SZA reuse depth is the documented "minimal" form: per-SZA we
#     re-run `model_from_parameters`. The expensive Mie + τ_abs portion
#     of model construction can be amortised across SZAs in a follow-up
#     (see the plan's `SceneOptics` design notes).
# =============================================================================

"""
    ScenarioSweep{FT}

Declarative cross product (SZA × view-pair × BRDF) over a single scene.
- `sza :: Vector{FT}` — crossed
- `view_pairs :: Vector{Tuple{FT,FT}}` — paired (vza, vaz). `postprocessing_vza!`
  produces one result per pair, not a Cartesian product.
- `brdfs :: Vector{AbstractSurfaceType}` — crossed

Build with the keyword constructor [`ScenarioSweep(; sza, vza, vaz, brdfs)`](@ref),
which accepts either paired (`length(vza) == length(vaz)`) or crossed
(`length(vza) != length(vaz)`) inputs.
"""
struct ScenarioSweep{FT}
    sza::Vector{FT}
    view_pairs::Vector{Tuple{FT, FT}}
    brdfs::Vector{<:AbstractSurfaceType}
end

"""
    ScenarioSweep(; sza, view_pairs, brdfs, FT=Float64) -> ScenarioSweep
    ScenarioSweep(; sza, vza, vaz, brdfs, view_mode::Symbol, FT=Float64) -> ScenarioSweep

Convenience constructors:

- Pass `view_pairs` as a `Vector{Tuple{FT,FT}}` to specify viewing directions
  directly — unambiguous, recommended for benchmark drivers.

- Or pass `vza` and `vaz` separately with an **explicit** `view_mode`:
  - `view_mode = :paired` — `vza[i]` paired with `vaz[i]` (requires
    `length(vza) == length(vaz)`).
  - `view_mode = :cross` — Cartesian product, row-major over (vza, vaz).

The constructor never infers paired-vs-crossed from vector lengths
(historical footgun: a 5-vza × 5-vaz sweep would silently collapse to
5 paired directions). Callers must commit.
"""
function ScenarioSweep(; sza::AbstractVector,
                         view_pairs::Union{Nothing, AbstractVector} = nothing,
                         vza::Union{Nothing, AbstractVector} = nothing,
                         vaz::Union{Nothing, AbstractVector} = nothing,
                         view_mode::Union{Nothing, Symbol} = nothing,
                         brdfs::AbstractVector{<:AbstractSurfaceType},
                         FT::DataType = Float64)
    sza_v = FT.(sza)

    pairs::Vector{Tuple{FT, FT}} = if view_pairs !== nothing
        (vza !== nothing || vaz !== nothing || view_mode !== nothing) &&
            throw(ArgumentError("ScenarioSweep: pass either `view_pairs=` OR `(vza, vaz, view_mode=)`, not both."))
        [(FT(p[1]), FT(p[2])) for p in view_pairs]

    elseif vza !== nothing || vaz !== nothing
        (vza === nothing || vaz === nothing) &&
            throw(ArgumentError("ScenarioSweep: when using `vza`/`vaz`, both are required."))
        view_mode === nothing &&
            throw(ArgumentError("ScenarioSweep: `view_mode` is required when passing `(vza, vaz)`. Use :paired or :cross."))
        vza_v = FT.(vza); vaz_v = FT.(vaz)
        if view_mode === :paired
            length(vza_v) == length(vaz_v) ||
                throw(ArgumentError("ScenarioSweep(view_mode=:paired): length(vza)=$(length(vza_v)) ≠ length(vaz)=$(length(vaz_v))."))
            [(vza_v[i], vaz_v[i]) for i in eachindex(vza_v)]
        elseif view_mode === :cross
            # Row-major Cartesian product: outer over vza, inner over vaz
            [(v, a) for v in vza_v for a in vaz_v]
        else
            throw(ArgumentError("ScenarioSweep: view_mode must be :paired or :cross, got :$view_mode"))
        end

    else
        throw(ArgumentError("ScenarioSweep: pass either `view_pairs=` or `(vza=, vaz=, view_mode=)`."))
    end

    # Preserve the concrete BRDF element type when all surfaces are uniform —
    # `Vector{AbstractSurfaceType}` would force dynamic dispatch in
    # `create_surface_layer!` inside the sweep's hot path.
    return ScenarioSweep{FT}(sza_v, pairs, collect(brdfs))
end

Base.show(io::IO, s::ScenarioSweep{FT}) where {FT} = print(io,
    "ScenarioSweep{", FT, "}(", length(s.sza), " SZA × ",
    length(s.view_pairs), " view-pair × ", length(s.brdfs), " BRDF)")

"""
    SweepResult{FT}

Output of [`run_sweep`](@ref). `R` and `T` are 5-D arrays of shape
`(N_sza, N_view_pair, N_brdf, pol_n, nSpec)`. The `sweep` field carries
the input axes for downstream introspection.
"""
struct SweepResult{FT}
    R::Array{FT, 5}
    T::Array{FT, 5}
    sweep::ScenarioSweep{FT}
end

Base.show(io::IO, r::SweepResult{FT}) where {FT} = print(io,
    "SweepResult{", FT, "}(R/T size=", size(r.R), ", sweep=", r.sweep, ")")

# =============================================================================
# SceneOptics — minimal wrapper (SZA reuse depth: rebuild per SZA)
# =============================================================================

"""
    SceneOptics{ARCH, FT}

Lightweight wrapper carrying the base parameters + a pre-built model for
one scene. Used by [`run_sweep`](@ref) to iterate over SZAs without the
caller having to re-supply the full parameter set. The SceneOptics
abstraction is kept thin in this PR: each SZA re-runs `model_from_parameters`.
A future optimisation will cache the μ₀-independent optics (`τ_abs`,
Mie greek coefs); see the plan notes.
"""
struct SceneOptics{ARCH<:AbstractArchitecture, FT<:AbstractFloat}
    model::RTModel{ARCH, FT}
    base_params::vSmartMOM_Parameters{FT}
end

"""
    scene_optics(model, base_params) -> SceneOptics

Bundle a built model with its source parameters so [`model_for_sza`](@ref)
can produce per-SZA models from the same scene.
"""
scene_optics(model::RTModel{ARCH, FT}, params::vSmartMOM_Parameters{FT}) where {ARCH, FT} =
    SceneOptics{ARCH, FT}(model, params)

"""
    model_for_sza(so::SceneOptics, sza::Real;
                  vza=nothing, vaz=nothing) -> RTModel

Build a new model for `so`'s scene at a different SZA (and optionally
different viewing geometry). Currently rebuilds the entire model — Mie
+ τ_abs reuse is a documented follow-up.
"""
function model_for_sza(so::SceneOptics{ARCH, FT}, sza::Real;
                       vza::Union{Nothing, AbstractVector}=nothing,
                       vaz::Union{Nothing, AbstractVector}=nothing) where {ARCH, FT}
    p = deepcopy(so.base_params)
    p.sza = FT(sza)
    vza !== nothing && (p.vza = FT.(collect(vza)))
    vaz !== nothing && (p.vaz = FT.(collect(vaz)))
    return model_from_parameters(p)
end

# =============================================================================
# Sweep driver
# =============================================================================

"""
    run_sweep(so::SceneOptics, sweep::ScenarioSweep;
              i_band=1, sources=nothing) -> SweepResult

Execute the (SZA × view-pair × BRDF) sweep. For each SZA, builds a model
via [`model_for_sza`](@ref), then calls [`rt_run_multi_surface`](@ref) to
amortise the atmosphere phase across all BRDFs.

Returns a [`SweepResult`](@ref) with `R`/`T` shaped
`(N_sza, N_view_pair, N_brdf, pol_n, nSpec)`.
"""
function run_sweep(so::SceneOptics{ARCH, FT}, sweep::ScenarioSweep{FT};
                   i_band::Integer = 1,
                   sources::Union{Nothing, AbstractSource} = nothing) where {ARCH, FT}
    N_sza  = length(sweep.sza)
    N_vp   = length(sweep.view_pairs)
    N_brdf = length(sweep.brdfs)
    vza_vec = FT[first(p) for p in sweep.view_pairs]
    vaz_vec = FT[last(p)  for p in sweep.view_pairs]

    # Output shapes are fully known from the scene before the loop runs —
    # preallocate so the loop body is type-stable end-to-end.
    pol_n = CoreRT.polarization_type(so.model).n
    nSpec = size(so.model.τ_abs[i_band], 1)
    R_all = zeros(FT, N_sza, N_vp, N_brdf, pol_n, nSpec)
    T_all = zeros(FT, N_sza, N_vp, N_brdf, pol_n, nSpec)

    for (i_sza, sza_val) in enumerate(sweep.sza)
        model = model_for_sza(so, sza_val; vza=vza_vec, vaz=vaz_vec)
        out = rt_run_multi_surface(model, sweep.brdfs;
                                   i_band=i_band, sources=sources)
        # `out[i_brdf]` is (R, T, ieR, ieT, hdr, bhr_uw, bhr_dw);
        # `out[i_brdf][1]` is R shaped (N_vp, pol_n, nSpec).
        for i_brdf in 1:N_brdf
            R_all[i_sza, :, i_brdf, :, :] .= out[i_brdf][1]
            T_all[i_sza, :, i_brdf, :, :] .= out[i_brdf][2]
        end
    end

    return SweepResult{FT}(R_all, T_all, sweep)
end
