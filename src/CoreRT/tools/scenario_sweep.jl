# =============================================================================
# Scenario sweeps: SZA × view-pair × BRDF over one scene
# (proposals/surface_split_albedo_sweep.md §6/§7, PR 3)
#
# Declarative description of the cross product (SZA × view-pair × BRDF) to
# evaluate for one scene, plus a driver that executes it efficiently:
#
#   * VZA + VAZ are PAIRED (one viewing direction per index), matching
#     `postprocessing_vza!` (tools/postprocessing_vza.jl) — it produces one
#     result per (vza, vaz) pair, not a Cartesian product. The constructor
#     accepts paired or crossed forms and requires an EXPLICIT `view_mode`
#     when passing separate `vza`/`vaz` vectors — ported verbatim from the
#     `gchp-io` prototype's footgun guard: a 5-vza × 5-vaz sweep must not
#     silently collapse to 5 paired directions just because the vectors
#     happen to be the same length as a Cartesian expansion.
#   * Surfaces are free per-SZA: `rt_run_multi_surface` builds the
#     atmosphere once per SZA and replays the surface step per BRDF
#     (surface_split_albedo_sweep.md §0/§2).
#   * SZA reuse depth: the `gchp-io` prototype rebuilt the ENTIRE model per
#     SZA via `model_from_parameters` (re-running HITRAN absorption + Mie
#     every SZA — see the prototype's `SceneOptics`/`model_for_sza`). On
#     `main` today, [`remake_geometry`](@ref) (tools/model_from_parameters.jl)
#     replaces that: only `geometry::ObsGeometry` and `quad_points::QuadPoints`
#     depend on `(sza, vza, vaz)` (verified bit-exact against a full
#     `model_from_parameters` rebuild in `test/test_scenario_sweep.jl`; see
#     `remake_geometry`'s docstring for the argument). `run_sweep` below
#     therefore shares one `atmosphere`/`optics`/`solver` across every SZA in
#     the sweep and only rebuilds the cheap geometry/quadrature per SZA —
#     the per-SZA HITRAN/Mie re-run the prototype paid is gone.
# =============================================================================

"""
    ScenarioSweep{FT}

Declarative cross product (SZA × view-pair × BRDF) over a single scene.
- `sza :: Vector{FT}` — crossed
- `view_pairs :: Vector{Tuple{FT,FT}}` — paired (vza, vaz). `postprocessing_vza!`
  produces one result per pair, not a Cartesian product.
- `brdfs :: Vector{<:AbstractSurfaceType}` — crossed

Build with the keyword constructor [`ScenarioSweep(; sza, vza, vaz, brdfs)`](@ref),
which accepts either paired (`view_mode = :paired`) or crossed
(`view_mode = :cross`) `(vza, vaz)` inputs, or an explicit `view_pairs` list.
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

- Pass `view_pairs` as a `Vector{<:Tuple}` to specify viewing directions
  directly — unambiguous, recommended for benchmark drivers.

- Or pass `vza` and `vaz` separately with an **explicit** `view_mode`:
  - `view_mode = :paired` — `vza[i]` paired with `vaz[i]` (requires
    `length(vza) == length(vaz)`).
  - `view_mode = :cross` — Cartesian product, row-major over (vza, vaz).

The constructor never infers paired-vs-crossed from vector lengths
(historical footgun, carried over from the `gchp-io` prototype: a
5-vza × 5-vaz sweep would silently collapse to 5 paired directions if
`view_mode` were inferred instead of required). Callers must commit.
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
# Sweep driver
# =============================================================================

"""
    run_sweep(model::RTModel, params::vSmartMOM_Parameters, sweep::ScenarioSweep;
             i_band::Integer = 1, sources = nothing) -> SweepResult

Execute the (SZA × view-pair × BRDF) sweep declared by `sweep` over the
scene carried by `model`/`params`. For each SZA, [`remake_geometry`](@ref)
builds a cheap per-geometry `RTModel` — sharing `atmosphere`/`optics`/
`solver`/`surfaces` with `model`, no HITRAN/Mie re-run — then
[`rt_run_multi_surface`](@ref) amortises the atmosphere phase (elemental →
doubling → interaction over `Nz` layers) across every BRDF in
`sweep.brdfs` within that SZA.

`model` and `params` must satisfy [`remake_geometry`](@ref)'s contract:
`model` was built via `model_from_parameters(params)` (or is otherwise
geometry-consistent with it — same `quadrature_type`/`l_trunc`/
`polarization_type`/`architecture`). `sweep.view_pairs` supplies ONE
shared `(vza, vaz)` list reused at every SZA — `rt_run_multi_surface`
runs every view pair in a single model call (they are just entries in the
model's `vza`/`vaz` vectors), so per-SZA cost is `O(N_brdf)` surface
replays, not `O(N_view_pair · N_brdf)`.

Returns a [`SweepResult`](@ref) with `R`/`T` shaped
`(N_sza, N_view_pair, N_brdf, pol_n, nSpec)`.
"""
function run_sweep(model::RTModel{ARCH, FT}, params::vSmartMOM_Parameters,
                   sweep::ScenarioSweep{FT};
                   i_band::Integer = 1,
                   sources::Union{Nothing, AbstractSource} = nothing) where {ARCH, FT}
    N_sza  = length(sweep.sza)
    N_vp   = length(sweep.view_pairs)
    N_brdf = length(sweep.brdfs)
    vza_vec = FT[first(p) for p in sweep.view_pairs]
    vaz_vec = FT[last(p)  for p in sweep.view_pairs]

    # Output shapes are fully known from the scene before the loop runs —
    # preallocate so the loop body is type-stable end-to-end.
    pol_n = polarization_type(model).n
    nSpec = size(model.τ_abs[i_band], 1)
    R_all = zeros(FT, N_sza, N_vp, N_brdf, pol_n, nSpec)
    T_all = zeros(FT, N_sza, N_vp, N_brdf, pol_n, nSpec)

    for (i_sza, sza_val) in enumerate(sweep.sza)
        model_sza = remake_geometry(model, params; sza = sza_val, vza = vza_vec, vaz = vaz_vec)
        out = rt_run_multi_surface(model_sza, sweep.brdfs; i_band = i_band, sources = sources)
        # `out[i_brdf]` is (R, T, ieR, ieT, hdr, bhr_uw, bhr_dw);
        # `out[i_brdf][1]` is R shaped (N_vp, pol_n, nSpec).
        for i_brdf in 1:N_brdf
            R_all[i_sza, :, i_brdf, :, :] .= out[i_brdf][1]
            T_all[i_sza, :, i_brdf, :, :] .= out[i_brdf][2]
        end
    end

    return SweepResult{FT}(R_all, T_all, sweep)
end
