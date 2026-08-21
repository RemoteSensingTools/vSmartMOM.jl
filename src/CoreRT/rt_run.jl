#=

This file contains the entry point for running the RT simulation, rt_run. 

There are two implementations: one that accepts the raw parameters, and one that accepts
the model. The latter should generally be used by users. 

=#




"""
    rt_run(model::RTModel; i_band=1) -> (R, T, ...)

Run the forward radiative transfer solver for one or more spectral bands.

Performs polarized adding-doubling RT through the atmosphere defined in `model`,
computing top-of-atmosphere (TOA) reflectance and bottom-of-atmosphere (BOA)
transmittance.  The solver iterates over azimuthal Fourier moments
`m = 0, …, max_m-1`, building layer R/T/J matrices via elemental → doubling →
interaction steps, then applies surface coupling and postprocessing.

Equivalent to `rt_run(noRS(), model, i_band)` (no Raman scattering).

# Arguments
- `model::RTModel`: Pre-built model from [`model_from_parameters`](@ref).
- `i_band::Integer=1`: Spectral band index (or vector of indices) to compute.

# Returns
A tuple `(R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)` when using
source-function integration (default), where:
- `R_SFI::Array{FT,3}`: TOA reflectance `[nVZA × nStokes × nSpec]`
- `T_SFI::Array{FT,3}`: BOA transmittance `[nVZA × nStokes × nSpec]`

For most use cases, only the first two elements are needed:
```julia
R, T = rt_run(model)
```

# Example
```julia
params = parameters_from_yaml("config/my_scene.yaml")
model = model_from_parameters(params)
R, T = rt_run(model)
R[1, 1, :]  # Stokes-I reflectance at first VZA across the spectrum
```

# See also
- [`rt_run(model, lin_model, NAer, NGas, NSurf)`](@ref) for linearized RT with Jacobians.
- [`model_from_parameters`](@ref) to build the model from parameters.
- [`rt_run_atmosphere`](@ref) / [`rt_run_surface`](@ref) / [`rt_run_multi_surface`](@ref)
  to amortize the layer-accumulation cost across many surfaces over the same
  atmosphere. They are implemented on top of the `atm_snapshot_callback` /
  `stop_after_atmosphere` / `m_max_override` keywords below — most callers
  should use the higher-level functions rather than these directly.
"""
function rt_run(model; i_band::Integer = 1,
                sources::Union{Nothing, AbstractSource} = nothing,
                streams_callback::Union{Nothing, Function} = nothing,
                atm_snapshot_callback::Union{Nothing, Function} = nothing,
                stop_after_atmosphere::Bool = false,
                m_max_override::Union{Nothing, Int} = nothing)
    rt_run(InelasticScattering.noRS{float_type(model)}(), model, i_band;
           sources, streams_callback,
           atm_snapshot_callback, stop_after_atmosphere, m_max_override)
end

"""
    StreamRTResult{FT}

Per-Fourier-moment radiative-transfer result at all internal quadrature
streams — the output of [`rt_run_streams`](@ref). Downstream consumers
(e.g. ExoOptics disk integration) reconstruct the Stokes vector at
arbitrary `(μ_v, μ_0, Δφ)` by interpolating over `(μ_v, μ_0)` and
Fourier-summing over `m` with `cos(m·Δφ) / sin(m·Δφ)` weights, instead
of one `rt_run` per pixel.

# Fields
- `qp_μ :: Vector{FT}` — full quadrature mu nodes (length `Nquad`).
- `iμ₀ :: Int` — index into `qp_μ` of the SZA stream (the closest to `μ₀`).
- `μ₀ :: FT` — cosine of the chosen incident-beam SZA.
- `pol_n :: Int` — number of Stokes components (`pol_type.n`).
- `weight :: Vector{FT}` — Fourier weight `(0.5/π for m=0; 1/π otherwise)`
  used for each moment by the internal post-processor.
- `R⁻⁺_per_m :: Vector{Array{FT, 3}}` — `(NquadN, NquadN, nSpec)` per
  Fourier moment, where `NquadN = Nquad · pol_n`. This is the
  full-stream reflection matrix (stokes-out blocks × stokes-in blocks).
- `J⁻_per_m :: Vector{Array{FT, 3}}` — `(NquadN, 1, nSpec)` per Fourier
  moment, the combined per-source SFI Stokes at all output streams
  (legacy slot + all per-source slots summed). Use this when you want
  the SFI-style output without re-running RT.
- `J⁺_per_m :: Vector{Array{FT, 3}}` — analog for BOA transmittance.

# Recovering `rt_run`'s output
```julia
streams = rt_run_streams(model; sources = sources)
# At a single (vza, vaz):
iμ = nearest_point(streams.qp_μ, cosd(vza))
_, istart, iend = get_indices(iμ, pol_type)
R_recovered = zero(rt_run(model; sources)[1][1, :, :])  # (n_stokes, nSpec)
for (mi, m) in enumerate(0:length(streams.J⁻_per_m)-1)
    cosmφ = cosd(m * vaz); sinmφ = sind(m * vaz)
    w_stokes = pol_n == 1 ? streams.weight[mi] * cosmφ :
               streams.weight[mi] *
               Diagonal([cosmφ, cosmφ, sinmφ, sinmφ][1:pol_n])
    for s in 1:nSpec
        R_recovered[:, s] .+= w_stokes * streams.J⁻_per_m[mi][istart:iend, 1, s]
    end
end
```

A unit test pins the bit-exact agreement of this reconstruction against
`rt_run` for one published-figure geometry.
"""
struct StreamRTResult{FT}
    qp_μ      :: Vector{FT}
    iμ₀       :: Int
    μ₀        :: FT
    pol_n     :: Int
    weight    :: Vector{FT}
    R⁻⁺_per_m :: Vector{Array{FT, 3}}     # TOA reflection block (NquadN, NquadN, nSpec)
    T⁺⁺_per_m :: Vector{Array{FT, 3}}     # BOA transmission block (NquadN, NquadN, nSpec)
    J⁻_per_m  :: Vector{Array{FT, 3}}     # SFI source-convolved upwelling (NquadN, 1, nSpec)
    J⁺_per_m  :: Vector{Array{FT, 3}}     # SFI source-convolved downwelling
    # Optical-depth profile (m-independent, captured once on m=0). Shape
    # `(nSpec, nLayer)`. Useful for ExoOptics post-processing — photon
    # budget, weighting-function decompositions, retrieval Jacobians.
    τ_total   :: Array{FT, 2}             # τ_rayl + τ_abs + (τ_aer if active)
    τ_rayl    :: Array{FT, 2}
    τ_abs     :: Array{FT, 2}
end

"""
    rt_run_streams(model; i_band=1, sources=nothing) -> StreamRTResult

Run the RT solver and return per-Fourier-moment Stokes matrices at all
quadrature streams instead of post-processed `(vza, vaz)` outputs. See
[`StreamRTResult`](@ref) for the data layout and a worked recovery example.

Internally just calls [`rt_run`](@ref) with a `streams_callback` that
copies `composite_layer.R⁻⁺`, `composite_layer.J₀⁺/⁻`, and the combined
per-source-slot SFI contributions out of the live layer accumulators
once per Fourier moment.
"""
function rt_run_streams(model; i_band::Integer = 1,
                         sources::Union{Nothing, AbstractSource} = nothing)
    FT = float_type(model)
    R⁻⁺_list = Vector{Array{FT, 3}}()
    T⁺⁺_list = Vector{Array{FT, 3}}()
    J⁻_list  = Vector{Array{FT, 3}}()
    J⁺_list  = Vector{Array{FT, 3}}()
    weights  = FT[]
    qp_μ_save = Vector{FT}()
    iμ₀_save = Ref(0)
    μ₀_save  = Ref(zero(FT))
    pol_n_save = Ref(0)

    cb = function (state)
        if isempty(qp_μ_save)
            append!(qp_μ_save, FT.(state.qp_μ))
            iμ₀_save[]   = state.iμ₀
            μ₀_save[]    = FT(state.μ₀)
            pol_n_save[] = state.pol_type.n
        end
        push!(weights, FT(state.weight))
        push!(R⁻⁺_list, Array{FT, 3}(_to_cpu(state.composite_layer.R⁻⁺)))
        push!(T⁺⁺_list, Array{FT, 3}(_to_cpu(state.composite_layer.T⁺⁺)))
        # Sum legacy slot + per-source slots so consumers see the same
        # combined SFI output that postprocessing_vza! accumulates into
        # R_SFI / T_SFI. The deepcopy on the legacy slot decouples our
        # stored copy from the live layer matrices (which the next m
        # iteration overwrites).
        J⁻_combined = Array{FT, 3}(_to_cpu(state.composite_layer.J₀⁻))
        J⁺_combined = Array{FT, 3}(_to_cpu(state.composite_layer.J₀⁺))
        for cslot in values(state.composite_layer.J₀_by_src)
            J⁻_combined .+= _to_cpu(cslot.J₀⁻)
            J⁺_combined .+= _to_cpu(cslot.J₀⁺)
        end
        push!(J⁻_list, J⁻_combined)
        push!(J⁺_list, J⁺_combined)
    end

    rt_run(model; i_band, sources, streams_callback = cb)

    # Capture the τ profile from the model. These are m-independent
    # quantities the user wants for post-processing (photon budget,
    # weighting functions, retrieval Jacobians). Shape: (nSpec, nLayer).
    # Aerosol contribution is NOT included here — `model.optics.aerosols.τ_aer`
    # is stored as (n_aer, nLayer) which doesn't broadcast directly. A
    # spectral aerosol contribution would need the per-band σ_ext applied
    # at the same time; defer to a follow-up when aerosols are wired.
    τ_rayl  = Array{FT, 2}(_to_cpu(model.τ_rayl[i_band]))
    τ_abs   = Array{FT, 2}(_to_cpu(model.τ_abs[i_band]))
    τ_total = τ_rayl .+ τ_abs

    return StreamRTResult{FT}(
        qp_μ_save, iμ₀_save[], μ₀_save[], pol_n_save[],
        weights, R⁻⁺_list, T⁺⁺_list, J⁻_list, J⁺_list,
        τ_total, τ_rayl, τ_abs,
    )
end

"""
    rt_run_test(RS_type, model, iBand; kwargs...)

Test entry point for RT calculations with explicit Raman type.
"""
function rt_run_test(RS_type::AbstractRamanType, model, iBand; kwargs...)
    rt_run(RS_type, model, iBand; kwargs...)
end

"""
    _interaction_workspace(rs, composite_layer, added_layer, arch)

Allocate the inelastic interaction workspace only for Raman-aware modes. Pure
elastic modes dispatch to `nothing`.
"""
_interaction_workspace(::Union{noRS, noRS_plus}, composite_layer, added_layer, arch) = nothing
_interaction_workspace(::AbstractRamanType, composite_layer, added_layer, arch) =
    InteractionWorkspace(composite_layer, added_layer; staged = _use_staged_interaction(arch, composite_layer))

# Whether interaction!(ScatteringInterface_11) pages its inelastic outputs through
# CPU (`staged=true`: memory-frugal but PCIe-bound) or runs GPU-only
# (`staged=false`: ~6× faster, but needs ~2 extra big 4-D temps in VRAM).
# Default: non-staged. The fused per-Δn kernels cut transient GPU memory ~7×, so
# non-staged now fits in the common case and is much faster (measured 87→39 s end-
# to-end @ nSpec=10000). The CUDA extension overrides the GPU method to fall back to
# staged ONLY when the non-staged footprint would not fit the device's free memory.
_use_staged_interaction(arch, composite_layer) = false

"""
    _expand_layer_rayleigh!(rs, fScattRayleigh, iz)

Populate per-layer Rayleigh scattering fractions for Raman-aware modes. Pure
elastic modes are a no-op.
"""
_expand_layer_rayleigh!(::Union{noRS, noRS_plus}, fScattRayleigh, iz) = nothing
function _expand_layer_rayleigh!(RS_type::AbstractRamanType, fScattRayleigh, iz)
    RS_type.fscattRayl = expandBandScalars(RS_type, fScattRayleigh[iz])
    return nothing
end

function _warn_explicit_depol_raman(RS_type::AbstractRamanType, model)
    if InelasticScattering.has_inelastic(RS_type) && model.solver.depol >= 0
        @warn "Raman-active RT is using an explicit nonnegative depol; current model construction applies this same depol to Rayleigh, Cabannes, and tau_rayl. Use depol: -1 for molecular Rayleigh/Cabannes auto mode." depol=model.solver.depol maxlog=1
    end
    return nothing
end

"""
    rt_run(RS_type, model::RTModel, iBand; sources=nothing)

Perform Radiative Transfer calculations with explicit Raman type.

Accepts an optional `sources::AbstractSource` argument (v0.6 source-term
refactor, Phase 2). When `nothing`, the legacy F₀ default is preserved
bit-for-bit. When a [`SourceSet`](@ref) (or a single
[`AbstractSource`](@ref)) is supplied, the F₀ carried by the first
[`SolarBeam`](@ref) is routed into `RS_type.F₀` ahead of the existing
SFI kernel call. Phase 5 will remove the `RS_type.F₀` indirection.

Also accepts the atmosphere/surface-split hooks consumed by
[`rt_run_atmosphere`](@ref): `atm_snapshot_callback` (fired once per Fourier
moment right after the layer loop, before the surface step),
`stop_after_atmosphere` (skip the surface step + postprocessing for this
call), and `m_max_override` (widen, never narrow, the Fourier loop bound).
All three default to a no-op / bit-exact-compatible value.
"""
function rt_run(RS_type::AbstractRamanType, model, iBand;
                sources::Union{Nothing, AbstractSource} = nothing,
                streams_callback::Union{Nothing, Function} = nothing,
                atm_snapshot_callback::Union{Nothing, Function} = nothing,
                stop_after_atmosphere::Bool = false,
                m_max_override::Union{Nothing, Int} = nothing)
    _warn_explicit_depol_raman(RS_type, model)

    # Apply the per-model BLAS thread cap once per `rt_run` invocation
    # (no-op when `model.numerics.blas_threads === nothing`). Lives here
    # so swapping models with different caps "just works" — the caller
    # doesn't have to remember to re-pin BLAS between runs.
    if model.numerics.blas_threads !== nothing
        LinearAlgebra.BLAS.set_num_threads(model.numerics.blas_threads)
    end

    (; obs_alt, sza, vza, vaz) = model.obs_geom   # Observational geometry properties
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀, iμ₀, Nquad) = model.quad_points # All quadrature points
    pol_type = CoreRT.polarization_type(model)
    # Per-band Fourier loop bound (order). Multi-band runs loop to the
    # max across bands; per-component skipping inside the loop is a
    # Phase C concern.
    m_max = maximum(m_max_bands(model)[ib] for ib in iBand)
    # Atmosphere/surface split (`rt_run_atmosphere`) — allow callers to widen
    # the Fourier loop bound without mutating `model.solver`. Only raises the
    # cap: silently lowering it would truncate scattering that the model was
    # built to resolve.
    if m_max_override !== nothing
        m_max = max(m_max, m_max_override)
    end
    (; quad_points) = model
    FT       = CoreRT.float_type(model)
    dτ_max_threshold = model.numerics.dτ_max_threshold   # numerical knob → rt_kernel!
    dτ_min_floor     = model.numerics.dτ_min_floor

    n_aer = CoreRT.n_aerosols(model)

    # Also to be changed if more than 1 band is used!!
    brdf = get_surface(model, iBand[1])
    if length(iBand) > 1
        @info "More than one band has been chosen, be aware that multiple BRDFs are not yet implemented and only the first one will be used!"
    end

    (; ϖ_Cabannes) = RS_type

    # Normalize ϖ_λ₁λ₀ so its sum equals the Raman fraction of scattering
    # (1 - ϖ_Cabannes) for this band. Missing on unified → inelastic ieR/ieT
    # was off by ~4× relative to sanghavi reference (see plans/PHASE_1B_STAGING.md §8).
    # Ported from sanghavi/src/CoreRT/rt_run.jl:293.
    InelasticScattering.normalize_raman_weights!(RS_type, model, iBand)

    Nz = length(model.profile.p_full)   # Number of vertical slices

    RS_type.bandSpecLim = UnitRange{Int}[]
    nSpec = 0;
    for iB in iBand
        nSpec0 = nSpec+1;
        nSpec += size(model.τ_abs[iB], 1); # Number of spectral points
        push!(RS_type.bandSpecLim,nSpec0:nSpec);
    end

    arr_type = CoreRT.array_type(model) # Type of array to use
    arch     = CoreRT.architecture(model)
    SFI = true                          # SFI flag
    NquadN =  Nquad * pol_type.n         # Nquad (multiplied by Stokes n)
    dims   = (NquadN,NquadN)              # nxn dims

    # Output arrays for reflected and transmitted solar irradiation at TOA and BOA
    @timeit "Arrays"  R       = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  T       = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  R_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  T_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  ieR_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  ieT_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  hdr     = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  bhr_dw     = zeros(FT, pol_type.n, nSpec)
    @timeit "Arrays"  bhr_uw     = zeros(FT, pol_type.n, nSpec)
    @timeit "Arrays"  hdr_J₀⁻    = zeros(FT, length(vza), pol_type.n, nSpec)
    # Notify user of processing parameters
    msg =
    """
    Processing on: $(arch)
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((NquadN, NquadN, nSpec))
    """
    @info msg

    # Resolve sources EARLY so the layer allocators can size per-source slots
    # (v0.7 Phase A.2a). Resolution order is unchanged: explicit `sources=`
    # kwarg > `model.sources` > legacy `RS_type.F₀` (only when the latter is
    # already user-shaped, to preserve the historical
    # `rs.F₀ = ...; rt_run(rs, model, ...)` test pattern).
    # Legacy `RS_type.SIF₀` is no longer consumed by rt_run; use SurfaceSIF.
    # `RS_type.F₀` remains a compatibility fallback when no SolarBeam is supplied.
    effective_sources = sources === nothing ? model.sources : sources
    prepared_sources = prepare_sources(effective_sources, FT, pol_type.n, nSpec, arr_type)

    # Create arrays — pass `prepared_sources` so per-source j₀ / J₀ slots
    # (e.g. `:thermal`) get allocated alongside the legacy solar buffers.
    @timeit "Creating layers" added_layer         =
        make_added_layer(RS_type, FT, arr_type, dims, nSpec;
                         prepared_sources=prepared_sources)
    # Just for now, only use noRS here. The surface added-layer needs the
    # same per-source slots so per-source j₀⁻ injection works at the surface
    # (Phase A.2c will wire surface-emission contributions here).
    @timeit "Creating layers" added_layer_surface =
        make_added_layer(RS_type, FT, arr_type, dims, nSpec;
                         prepared_sources=prepared_sources)
    @timeit "Creating layers" composite_layer     =
        make_composite_layer(RS_type, FT, arr_type, dims, nSpec;
                             prepared_sources=prepared_sources)
    @timeit "Creating arrays" I_static =
        Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))));

        # Note: Raman SS properties (ϖ_λ₁λ₀, Z matrices, etc.) are set up
    # in model_from_parameters when RS_type != noRS.

    # Build concatenated wavenumber grid for canopy spectral features
    _canopy_spec_wn = nothing
    if brdf isa CanopySurface
        _canopy_spec_wn = vcat([get_spec_bands(model)[iB] for iB in iBand]...)
    end

    # Pre-initialize canopy cache before the Fourier loop (Zazi precomputation is expensive)
    if brdf isa CanopySurface && brdf._cache === nothing
        @timeit "Canopy cache init" _init_canopy_cache!(
            brdf, added_layer_surface, pol_type, quad_points, arch;
            spec_bands_wn=_canopy_spec_wn, m_max=m_max)
    end

    # Pre-compute within-canopy atmospheric optical depth if requested
    if brdf isa CanopySurface && brdf.include_atm && brdf.canopy_dp !== nothing
        @timeit "Canopy atm tau" _compute_canopy_atm_tau!(brdf, model, _canopy_spec_wn)
    end

    if sources === nothing && size(RS_type.F₀) == (pol_type.n, nSpec)
        # User has pre-set RS_type.F₀; leave it alone for back-compat.
    else
        F₀_dev = extract_solar_F₀(prepared_sources, FT, pol_type.n, nSpec, arr_type)
        # `RS_type.F₀` historically lives on host memory; keep that contract by
        # converting back through Array. (Cheap: F₀ is (pol_n, nSpec).)
        RS_type.F₀ = Array{FT, 2}(F₀_dev)
    end
    surface_F₀ = arr_type(FT.(RS_type.F₀))

    # Phase 4: allocate the InteractionWorkspace once before the layer loop for
    # RRS/VS runs to avoid the per-call GPU allocation (sanghavi reported ~7 GB
    # per call on high-spec Float64 configs). noRS skips the workspace — the
    # elastic interaction! path doesn't consume it. `staged=true` pages large
    # 4-D ie buffers through CPU between passes; only beneficial on GPU where
    # device memory pressure matters. CPU runs use staged=false (CPU↔CPU copies
    # would be pure overhead).
    _interaction_ws = _interaction_workspace(RS_type, composite_layer, added_layer, arch)

    # --- Pre-loop hoists (m-independent work) ---
    # Build the m-invariant cache: device uploads of τ_rayl/τ_abs per layer,
    # pre-applied δ-M aerosol corrections, fScattRayleigh (requires one m=0 Z
    # pass whose Z matrices are then discarded), and the host copy of qp_μ.
    # Cost: equivalent to one constructCoreOpticalProperties call; paid once
    # instead of once per Fourier moment.
    @timeit "CacheInit" _m_inv_cache = build_m_invariant_cache(RS_type, iBand, model)

    # Cumulative optical depth (m-independent, saved for TMS correction).
    # Computed once from m=0 optical props (Z doesn't enter τ_sum_all).
    τ_sum_all = nothing

    # Device-resident τ_sum_all[:,end] slice used by create_surface_layer!.
    # Hoisted out of the m loop — τ_sum_all is m-independent, so the upload
    # only needs to happen once (after τ_sum_all is first computed below).
    τ_sum_end_dev = nothing

    # Scattering interface sequence (m-independent: detection uses τ·ϖ only).
    # Computed once on m=0 and reused for all subsequent moments.
    scattering_interfaces_all = nothing

    # Pre-allocate hdr_J₀⁻ once; interaction_hdrf! writes into it each moment.
    # The original `similar(composite_layer.J₀⁻)` inside the loop allocated a
    # fresh array each moment (and the pre-allocated hdr_J₀⁻ at line 303 was
    # never used because it has shape (nVZA, pol_n, nSpec), not (NquadN,1,nSpec)).
    hdr_J₀⁻ = similar(composite_layer.J₀⁻)

    # Per-layer max(τ·ϖ), computed ONCE on m = 0 and reused for every Fourier
    # moment: τ and ϖ are m-independent (only the Z matrices change with m),
    # yet the scatter test + dτ/ndoubl derivation historically re-reduced
    # τ·ϖ on the device 3× per (layer, m) — each a broadcast + reduction +
    # BLOCKING scalar read, ~3·Nz·(m_max+1) GPU round trips per solve
    # (nsys-measured as the dominant sync-bound term). Same scalar, same
    # downstream arithmetic ⇒ bitwise-identical results.
    max_τϖ_all = _M_INVARIANT_DTAU_CACHE_ENABLED[] ? fill(FT(NaN), Nz) : nothing

    # Loop over fourier moments
    for m = 0:m_max

        # Azimuthal weighting
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)
        # Set the Zλᵢλₒ interaction parameters for Raman (or nothing for noRS)
        @timeit "IE"  InelasticScattering.computeRamanZλ!(RS_type, pol_type,collect(qp_μ), m, arr_type)
        # Compute the core layer optical properties (cache-aware: only Z moments
        # are recomputed; τ/ϖ uploads and fScattRayleigh come from _m_inv_cache).
        @timeit "OpticalProps" layer_opt_props, fScattRayleigh   =
            constructCoreOpticalProperties(RS_type, iBand, m, model, _m_inv_cache);
        # τ_sum_all and scattering_interfaces_all are m-independent (extractEffectiveProps
        # only uses τ·ϖ from layer_opt_props, not Z). Compute on m=0, reuse thereafter.
        if τ_sum_all === nothing
            @timeit "Extract Optical Properties" scattering_interfaces_all, τ_sum_all =
                extractEffectiveProps(layer_opt_props, quad_points)
            # Upload the BOA τ_sum slice to device once.
            τ_sum_end_dev = arr_type(τ_sum_all[:,end])
        end
        # On m>0 we reuse the cached scattering_interfaces_all and τ_sum_all.

        # Loop over vertical layers:
        @showprogress 1 "Looping over layers ..." for iz = 1:Nz  # Count from TOA to BOA

            # Construct the atmospheric layer
            @timeit "Expand Bands" _expand_layer_rayleigh!(RS_type, fScattRayleigh, iz)

            # Expand all layer optical properties to their full dimension:
            @timeit "OpticalProps" layer_opt =
                expandOpticalProperties(layer_opt_props[iz], arr_type)

            # m-invariant max(τ·ϖ) for this layer (see the cache's comment
            # above the m-loop): one device reduction per layer on m = 0,
            # host-cached scalar thereafter.
            if max_τϖ_all !== nothing
                m == 0 && (max_τϖ_all[iz] = FT(maximum(layer_opt.τ .* layer_opt.ϖ)))
                mτϖ_iz = max_τϖ_all[iz]
            else
                mτϖ_iz = nothing
            end

            # Perform Core RT (doubling/elemental/interaction)
            @timeit "RT Kernel" rt_kernel!(RS_type, pol_type, SFI,
                        added_layer, composite_layer,
                        layer_opt,
                        scattering_interfaces_all[iz],
                        τ_sum_all[:,iz],
                        m, quad_points,
                        I_static,
                        arch,
                        qp_μN, iz;
                        workspace=_interaction_ws,
                        prepared_sources=prepared_sources,
                        dτ_max_threshold=dτ_max_threshold,
                        dτ_min_floor=dτ_min_floor,
                        max_τϖ=mτϖ_iz)
        end

        # Atmosphere/surface split (`rt_run_atmosphere` / `rt_run_surface`).
        # Fire the optional callback once the iz-loop has settled and the
        # composite layer holds the atmosphere-only state for this m — BEFORE
        # the surface step below. `rt_run_atmosphere` uses this to snapshot
        # the pre-surface composite per Fourier moment; `stop_after_atmosphere`
        # then skips the surface block + postprocessing + SS correction
        # entirely (saving real work, not just discarding outputs). When
        # `atm_snapshot_callback` is `nothing` and `stop_after_atmosphere =
        # false` (the default), this is a single branch per m, no
        # allocations — bit-exact backwards compat with pre-split behaviour.
        if atm_snapshot_callback !== nothing
            atm_snapshot_callback((;
                m, weight, pol_type, qp_μ, iμ₀, μ₀,
                vza, vaz, nSpec, NquadN, iBand, SFI,
                composite_layer,
                scattering_interface_surf = scattering_interfaces_all[end],
                # Device-resident — the same array the monolithic surface
                # step consumes below. Snapshotting the host slice
                # τ_sum_all[:, end] instead would hand rt_run_surface a CPU
                # Vector to broadcast against device arrays on GPU models.
                τ_sum_surf = τ_sum_end_dev,
                surface_F₀,
                arr_type, arch, RS_type, prepared_sources, I_static,
                quad_points,
            ))
        end
        if stop_after_atmosphere
            continue
        end

        # Create surface matrices:
        if brdf isa CanopySurface
            @timeit "Create Surface" create_surface_layer!(brdf,
                                added_layer_surface,
                                SFI, m,
                                pol_type,
                                quad_points,
                                τ_sum_end_dev,
                                arch;
                                spec_bands_wn=_canopy_spec_wn,
                                m_max=m_max,
                                F₀=surface_F₀)
        else
            @timeit "Create Surface" create_surface_layer!(brdf,
                                added_layer_surface,
                                SFI, m,
                                pol_type,
                                quad_points,
                                τ_sum_end_dev,
                                arch;
                                F₀=surface_F₀)
        end

        # Surface-emission sources are injected only through the source system.
        # The legacy `RS_type.SIF₀` path is intentionally ignored here to avoid
        # double-counting when a `SurfaceSIF` source is also present.
        surface_source_contribute!(prepared_sources, brdf, added_layer_surface, m, pol_type, arch)

        # One last interaction with surface:
        @timeit "interaction" interaction!(RS_type,
                                    scattering_interfaces_all[end],
                                    SFI,
                                    composite_layer,
                                    added_layer_surface,
                                    I_static;
                                    workspace=_interaction_ws)
        # hdr_J₀⁻ is pre-allocated before the m loop (avoids one similar() per moment).
        @timeit "interaction_HDRF" interaction_hdrf!(#RS_type,
                                    #bandSpecLim,
                                    #scattering_interfaces_all[end], 
                                    SFI, 
                                    composite_layer, 
                                    added_layer_surface, 
                                    m, pol_type, quad_points,
                                    hdr_J₀⁻, bhr_uw, bhr_dw)
        
        # Postprocess and weight according to vza
        @timeit "Postprocessing VZA" postprocessing_vza!(RS_type, 
                            iμ₀, pol_type, 
                            composite_layer, 
                            vza, qp_μ, m, vaz, μ₀, 
                            weight, nSpec, 
                            SFI, 
                            R, R_SFI, 
                            T, T_SFI,
                            ieR_SFI, ieT_SFI)

        @timeit "Postprocessing HDRF" postprocessing_vza_hdrf!(RS_type,
            iμ₀, pol_type,
            hdr_J₀⁻,
            vza, qp_μ, m, vaz, μ₀,
            weight, nSpec,
            hdr)

        # Phase H — per-moment streams export hook (v0.7+).
        # Optional callback called once per Fourier moment AFTER the layer
        # accumulators have settled but BEFORE the moment-aware quantities
        # are overwritten by the next m iteration. The callback receives
        # everything it needs to reproduce the post-processing offline:
        # the moment index `m`, the Fourier weight `w` that the internal
        # post-processor used, the polarization type (so it knows the
        # stokes layout), the quadrature mu nodes + iμ₀ (the SZA stream
        # index), μ₀ itself, and the composite_layer. Downstream
        # consumers (ExoOptics disk integration, custom post-processors)
        # use the per-moment R⁻⁺ / J₀⁻ stored in `composite_layer` to do
        # their own Fourier sum + interpolation onto arbitrary
        # (μ_v, μ_0, Δφ) without needing one rt_run per (sza, vza, vaz).
        #
        # Cost when the callback is `nothing`: a single branch in the m
        # loop, no allocations — bit-exact backwards compat.
        if streams_callback !== nothing
            streams_callback((;
                m, weight, pol_type, qp_μ, iμ₀, μ₀,
                composite_layer, nSpec))
        end
    end

    # Single-scattering correction for Cox-Munk specular hotspot (TMS).
    # LIMITATION: this TMS correction assumes a unit incident beam per spectral
    # point. With a non-unit surface F₀ (SolarBeam / RS_type.F₀) the Fourier
    # surface layer is F₀-scaled but this correction is not, so Cox-Munk SFI with
    # a real solar spectrum is inconsistent. The default (unit F₀) path is
    # unaffected. TODO (with Sanghavi): pass surface_F₀ here and scale the glint
    # correction per spectral point to complete the F₀ surface-source feature.
    if brdf isa CoxMunkSurface && SFI && !stop_after_atmosphere
        @timeit "SS Correction" apply_ss_correction!(
            R_SFI, brdf, pol_type, vza, vaz, μ₀,
            Array(τ_sum_all[:,end]), m_max, nSpec)
    end

    # Show timing statistics (only when the user asked — verbose flag in
    # `RTNumericalParameters`; default false to keep production loops quiet).
    model.numerics.verbose && print_timer()
    reset_timer!()

    # Return R_SFI or R, depending on the flag
    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw[1,:], bhr_dw[1,:]) : (R, T)
end


# =========================================================================
# Single-scatter approximation driver
# =========================================================================

"""
    rt_run_ss(model::RTModel; i_band=1) -> (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hem_R, hem_T)

Run the single-scatter approximation forward RT solver for one or more bands.

Computes only the single-scattering component of the TOA reflectance /
BOA transmittance by replacing the layer-doubling kernel with `rt_kernel_ss!`
and using `interaction_ss!` for the final
surface coupling. Useful as a fast/debug reference path; not a physically
complete RT solution (multiple scattering is dropped).

Equivalent to `rt_run_ss(noRS(), model, i_band)`.

Ported from sanghavi-branch `rt_run_ss` (see `plans/IMPLEMENTATION_PLAN_v2.md`
Phase 1c).

# Arguments
- `model::RTModel`: Pre-built model from [`model_from_parameters`](@ref).
- `i_band::Integer=1`: Spectral band index (or vector of indices).

# Returns
6-tuple `(R_SFI, T_SFI, ieR_SFI, ieT_SFI, hem_R, hem_T)`:
- `R_SFI`, `T_SFI`, `ieR_SFI`, `ieT_SFI`: shape `(nVza, nStokes, nSpec)` — same
  per-direction outputs as [`rt_run`](@ref).
- `hem_R`, `hem_T`: hemispherical-integrated TOA reflectance / BOA transmittance,
  shape `(nSpec,)`. Computed from the m=0 Fourier coefficient of the Stokes-I
  source function over the full upper hemisphere (Σⱼ J₀[j, 1, λ] · μⱼ · wⱼ).
"""
function rt_run_ss(model; i_band::Integer = 1, sources::Union{Nothing, AbstractSource} = nothing)
    rt_run_ss(InelasticScattering.noRS{float_type(model)}(), model, i_band; sources)
end

"""
    rt_run_test_ss(RS_type, model, iBand; kwargs...)

Test entry point for single-scatter RT with explicit Raman type.
"""
function rt_run_test_ss(RS_type::AbstractRamanType, model, iBand; kwargs...)
    rt_run_ss(RS_type, model, iBand; kwargs...)
end

"""
    rt_run_ss(RS_type, model::RTModel, iBand)

Single-scatter approximation driver with explicit Raman type. See
[`rt_run_ss`](@ref) for the user-facing entry point.
"""
function rt_run_ss(RS_type::AbstractRamanType, model, iBand;
                   sources::Union{Nothing, AbstractSource} = nothing)
    _warn_explicit_depol_raman(RS_type, model)

    # Per-model BLAS thread cap (see `rt_run` body for rationale).
    if model.numerics.blas_threads !== nothing
        LinearAlgebra.BLAS.set_num_threads(model.numerics.blas_threads)
    end
    # Numerics knobs threaded into rt_kernel_ss! (matches rt_run forward path).
    dτ_max_threshold = model.numerics.dτ_max_threshold
    dτ_min_floor     = model.numerics.dτ_min_floor

    (; vza, vaz) = model.obs_geom
    (; qp_μ, wt_μ, qp_μN, μ₀, iμ₀, Nquad) = model.quad_points
    pol_type = CoreRT.polarization_type(model)
    # Per-band Fourier loop bound (order). Multi-band runs loop to the
    # max across bands; per-component skipping is a Phase C concern.
    m_max = maximum(m_max_bands(model)[ib] for ib in iBand)
    (; quad_points) = model
    FT       = CoreRT.float_type(model)

    brdf = get_surface(model, iBand[1])
    if length(iBand) > 1
        @info "More than one band has been chosen; single-scatter path uses only the first BRDF."
    end

    # Same ϖ_λ₁λ₀ normalization as rt_run (ported from sanghavi rt_run_ss:466).
    InelasticScattering.normalize_raman_weights!(RS_type, model, iBand)

    Nz = length(model.profile.p_full)

    RS_type.bandSpecLim = UnitRange{Int}[]
    nSpec = 0
    for iB in iBand
        nSpec0 = nSpec + 1
        nSpec += size(model.τ_abs[iB], 1)
        push!(RS_type.bandSpecLim, nSpec0:nSpec)
    end

    arr_type = CoreRT.array_type(model)
    arch     = CoreRT.architecture(model)
    SFI = true
    NquadN = Nquad * pol_type.n
    dims   = (NquadN, NquadN)

    @timeit "Arrays"  R       = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  T       = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  R_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  T_SFI   = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  ieR_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  ieT_SFI = zeros(FT, length(vza), pol_type.n, nSpec)
    @timeit "Arrays"  hem_R   = zeros(FT, nSpec)
    @timeit "Arrays"  hem_T   = zeros(FT, nSpec)

    msg = """
    Single-scatter mode — processing on: $(arch)
    With FT: $(FT)
    Source Function Integration: $(SFI)
    Dimensions: $((NquadN, NquadN, nSpec))
    """
    @info msg

    @timeit "Creating layers" added_layer         = make_added_layer(RS_type, FT, arr_type, dims, nSpec)
    @timeit "Creating layers" added_layer_surface = make_added_layer(RS_type, FT, arr_type, dims, nSpec)
    @timeit "Creating layers" composite_layer     = make_composite_layer(RS_type, FT, arr_type, dims, nSpec)
    @timeit "Creating arrays" I_static            = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))))

    # Resolve sources (v0.6 source-term refactor). Resolution: kwarg >
    # model.sources > pre-set `RS_type.F₀` for back-compat. See `rt_run`.
    effective_sources = sources === nothing ? model.sources : sources
    prepared_sources = prepare_sources(effective_sources, FT, pol_type.n, nSpec, arr_type)
    if sources === nothing && size(RS_type.F₀) == (pol_type.n, nSpec)
        # User-pre-set F₀ honored.
    else
        F₀_dev = extract_solar_F₀(prepared_sources, FT, pol_type.n, nSpec, arr_type)
        RS_type.F₀ = Array{FT, 2}(F₀_dev)
    end
    surface_F₀ = arr_type(FT.(RS_type.F₀))

    τ_sum_all = nothing

    for m = 0:m_max
        weight = m == 0 ? FT(0.5/π) : FT(1.0/π)

        @timeit "IE" InelasticScattering.computeRamanZλ!(RS_type, pol_type, collect(qp_μ), m, arr_type)

        @timeit "OpticalProps" layer_opt_props, fScattRayleigh =
            constructCoreOpticalProperties(RS_type, iBand, m, model)

        @timeit "Extract Optical Properties" scattering_interfaces_all, τ_sum_all =
            extractEffectiveProps(layer_opt_props, quad_points)

        @showprogress 1 "SS looping over layers ..." for iz = 1:Nz
            @timeit "Expand Bands" _expand_layer_rayleigh!(RS_type, fScattRayleigh, iz)

            @timeit "OpticalProps" layer_opt = expandOpticalProperties(layer_opt_props[iz], arr_type)

            @timeit "RT SS Kernel" rt_kernel_ss!(RS_type, pol_type, SFI,
                        added_layer, composite_layer,
                        layer_opt,
                        scattering_interfaces_all[iz],
                        τ_sum_all[:, iz],
                        m, quad_points,
                        I_static,
                        arch,
                        qp_μN, iz;
                        dτ_max_threshold=dτ_max_threshold,
                        dτ_min_floor=dτ_min_floor)
        end

        @timeit "Create Surface" create_surface_layer!(brdf,
                            added_layer_surface,
                            SFI, m,
                            pol_type,
                            quad_points,
                            arr_type(τ_sum_all[:, end]),
                            arch;
                            F₀=surface_F₀)

        # Surface-emission sources are injected only through the source system.
        # The legacy `RS_type.SIF₀` path is intentionally ignored here to avoid
        # double-counting when a `SurfaceSIF` source is also present.
        surface_source_contribute!(prepared_sources, brdf, added_layer_surface, m, pol_type, arch)

        # SS mode uses interaction_ss! (no multiple-scattering doubling with surface).
        τsurf = zeros(FT, length(τ_sum_all[:, Nz + 1]))
        @timeit "interaction_ss" interaction_ss!(SFI,
                            composite_layer,
                            added_layer_surface,
                            τ_sum_all[:, Nz + 1],
                            τsurf,
                            quad_points,
                            arch)

        @timeit "Postprocessing VZA" postprocessing_vza!(RS_type,
                            iμ₀, pol_type,
                            composite_layer,
                            vza, qp_μ, m, vaz, μ₀,
                            weight, nSpec,
                            SFI,
                            R, R_SFI,
                            T, T_SFI,
                            ieR_SFI, ieT_SFI)

        # Hemispherical integration — only m=0 contributes under 2π azimuthal
        # integration (higher Fourier moments vanish). Weight factor
        # (0.5/π) × 2π = 1 makes the raw Σⱼ J₀[j,1,λ] · μⱼ · wⱼ direct.
        if m == 0 && SFI
            J₀⁻_cpu = Array(composite_layer.J₀⁻)
            J₀⁺_cpu = Array(composite_layer.J₀⁺)
            μ_arr = Array(qp_μ)
            w_arr = Array(wt_μ)
            nStokes = pol_type.n
            @inbounds for s = 1:nSpec, j = 1:Nquad
                j_I = (j - 1) * nStokes + 1
                hem_R[s] += J₀⁻_cpu[j_I, 1, s] * μ_arr[j] * w_arr[j]
                hem_T[s] += J₀⁺_cpu[j_I, 1, s] * μ_arr[j] * w_arr[j]
            end
        end
    end

    model.numerics.verbose && print_timer()
    reset_timer!()

    return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hem_R, hem_T) : (R, T, hem_R, hem_T)
end
