"""
    _rayleigh_greek_source(rs, greek_rayleigh, greek_cabannes)

Select the Rayleigh phase-coefficient source for the elastic path. Pure
elastic modes use the full Rayleigh coefficients; Raman-aware modes use
Cabannes coefficients because Raman redistribution is handled explicitly.
"""
_rayleigh_greek_source(::Union{noRS, noRS_plus}, greek_rayleigh, greek_cabannes) = greek_rayleigh
_rayleigh_greek_source(::AbstractRamanType, greek_rayleigh, greek_cabannes) = greek_cabannes

_same_greek_coefficients(a, b) =
    all(name -> getfield(a, name) == getfield(b, name), fieldnames(typeof(a)))

"Compute diffuse phase operators and, for external solar, direct μ₀ blocks."
function _compute_phase_blocks(model, greek, m, arr_type; tables=nothing, Π_pair=nothing)
    q = model.quad_points
    pol = CoreRT.polarization_type(model)
    if q.external_solar
        μdiff = collect(q.qp_μ)
        # Tabulated P/R/T are valid only on the grid they were built for.
        t = (tables !== nothing && tables.μ == μdiff) ? tables : nothing
        Z⁺⁺, Z⁻⁺ = Scattering.compute_Z_moments(pol, μdiff, greek, m;
                                                  arr_type=arr_type, tables=t,
                                                  Π_pair=(t === nothing ? nothing : Π_pair))
        Z₀⁺, Z₀⁻ = Scattering.compute_Z_source_moments(
            pol, μdiff, q.μ₀, greek, m; arr_type=arr_type)
        return Z⁺⁺, Z⁻⁺, Z₀⁺, Z₀⁻
    end
    μphase = collect(q.phase_qp_μ)
    t = (tables !== nothing && tables.μ == μphase) ? tables : nothing
    Z⁺⁺, Z⁻⁺ = Scattering.compute_Z_moments(
        pol, μphase, greek, m; arr_type=arr_type, tables=t,
        Π_pair=(t === nothing ? nothing : Π_pair))
    return Z⁺⁺, Z⁻⁺, nothing, nothing
end

"Evaluate aerosol phase operators only at retained optical nodes, then interpolate in ν."
function _compute_aerosol_phase_blocks(model, optics, ν_spec, m, arr_type;
                                       tables=nothing, Π_pair=nothing)
    optics.phase_ν === nothing && return _compute_phase_blocks(
        model, optics.greek_coefs, m, arr_type; tables, Π_pair)
    blocks = [_compute_phase_blocks(model, g, m, Array; tables, Π_pair)
              for g in optics.phase_greek]
    interp(k) = arr_type(_interpolate_phase_nodes(ν_spec, optics.phase_ν,
                                                   [b[k] for b in blocks]))
    Z⁺⁺, Z⁻⁺ = interp(1), interp(2)
    if blocks[1][3] === nothing
        return Z⁺⁺, Z⁻⁺, nothing, nothing
    end
    return Z⁺⁺, Z⁻⁺, interp(3), interp(4)
end

"""
    _rayleigh_fraction_of_total_extinction(rayleigh_layer, total_layer)

Fraction of the total spectral extinction that participates in Rayleigh-driven
inelastic scattering. The denominator must include gas absorption because the
RRS elemental kernels attenuate with the total `dτ_λ`.
"""
_rayleigh_fraction_of_total_extinction(rayleigh_layer, total_layer) =
    Array(rayleigh_layer.τ ./ total_layer.τ)

"""
    MInvariantCache

Cache of m-independent optical quantities computed once before the Fourier
loop in [`rt_run`](@ref).  Passed as `cache=` to
[`constructCoreOpticalProperties`](@ref) to skip redundant uploads,
host↔device copies, and fScattRayleigh recomputation on every moment.

# Fields
- `μ_host`: host copy of the phase-evaluation cosines. In external-solar
  SFI this is the diffuse/operator grid augmented by exact `μ₀`.
- `rayl_τ_dev[iBi][iz]`: device τ_rayl slice per band-index/layer — avoids re-upload.
- `rayl_ϖ_Cabannes[iBi]`: ϖ_Cabannes per band (scalar or device array) — avoids
  re-upload.
- `aer_τ_mod[iBi][iAer][iz]`: δ-M-modified aerosol τ per band-index/aerosol/layer —
  precomputed `(1-fᵗ·ω̃)·τ_aer`; a device nSpec vector for multi-λ bands, or a scalar
  for single-λ bands.
- `aer_ϖ_mod[iBi][iAer]`: δ-M-modified aerosol ϖ — precomputed `(1-fᵗ)·ω̃/(1-fᵗ·ω̃)`;
  an nSpec vector for multi-λ bands, or a scalar for single-λ bands.
- `τ_abs_dev[iBi][iz]`: device τ_abs slice per band-index/layer — avoids re-upload.
- `fScattRayleigh`: per-layer Rayleigh scattering fraction (host arrays) — avoids
  the per-m `Array(rayl.τ ./ combo.τ)` device→host copy.

Note: indices `iBi` run 1:length(iBand), not over the raw band numbers in `iBand`.
"""
struct MInvariantCache{FT}
    μ_host          :: Vector{FT}
    rayl_τ_dev      :: Vector{Vector}         # [iBi][iz]      → device nSpec vector
    rayl_ϖ_Cabannes :: Vector                 # [iBi]          → ϖ_Cabannes (scalar or arr)
    aer_τ_mod       :: Vector{Vector{Vector}} # [iBi][iAer][iz] → plain FT scalar
    aer_ϖ_mod       :: Vector{Vector}         # [iBi][iAer]    → plain FT scalar
    τ_abs_dev       :: Vector{Vector}         # [iBi][iz]      → device nSpec vector
    fScattRayleigh  :: Vector{Vector{Vector}} # [iBi][iz]      → host FT vector
    "Per band: `mode_layers[i]` = the ONLY layer aerosol mode `i` occupies (layer-
    resolved diagonal τ structure), 0 for an all-zero mode; `nothing` when modes
    genuinely overlap layers (dense combine required). Detected from the host τ_aer
    once — enables the O(nAero) sparse combine instead of O(nAero × nZ)."
    mode_layers     :: Vector{Union{Nothing, Vector{Int}}}
    "Scatterer-independent Z-moment precompute (generalised spherical functions
    P/R/T on ±μ up to the global l_max) — see `Scattering.ZMomentTables`. Built
    once per solve; `nothing` when `_Z_TABLES_ENABLED[] == false`."
    z_tables        :: Union{Nothing, Scattering.ZMomentTables}
end

"""
Runtime kill-switch for the sparse (layer-resolved) aerosol combine. Default
`true`. NOTE the sparse path is tolerance-equivalent, not bitwise: the dense
path's zero-τ combines each round through `(ϖ·τ + 0)/τ`, which the sparse path
skips (making it marginally MORE exact).
"""
const _SPARSE_AERO_COMBINE_ENABLED = Ref(true)

"""
Runtime kill-switch for the tabulated Z-moment path (per-solve `ZMomentTables`
+ per-m shared Π lists instead of rebuilding the generalised spherical
functions for every (scatterer, m) pair). Default `true`. The tabulated path
is bitwise-identical to the plain one — same sums in the same order.
"""
const _Z_TABLES_ENABLED = Ref(true)

"Detect the layer-resolved diagonal structure of `τ_aer_band` [iAer, nSpec, iz]."
function _detect_mode_layers(τ_aer_band)
    τ_host = τ_aer_band isa Array ? τ_aer_band : Array(τ_aer_band)
    nAero, _, nZ = size(τ_host)
    out = zeros(Int, nAero)
    for i in 1:nAero
        layers = [iz for iz in 1:nZ if any(!iszero, @view τ_host[i, :, iz])]
        if length(layers) > 1
            return nothing            # overlapping modes → dense combine
        end
        out[i] = isempty(layers) ? 0 : layers[1]
    end
    return out
end

"""
    build_m_invariant_cache(RS_type, iBand, model) -> MInvariantCache

Pre-compute all Fourier-moment-independent quantities used by
[`constructCoreOpticalProperties`](@ref): device uploads of τ_rayl and τ_abs,
pre-applied δ-M aerosol corrections, the Rayleigh scattering fraction
`fScattRayleigh` (computed via a single m=0 Z pass whose Z matrices are then
discarded), and the host copy of the quadrature μ nodes.

Called once before the `m = 0:m_max` Fourier loop in `rt_run`; the resulting
`MInvariantCache` is passed into each `constructCoreOpticalProperties` call via
the `cache=` keyword, eliminating O(Nz × Nbands × m_max) redundant device
uploads/downloads and one `collect(qp_μ)` per moment.
"""
function build_m_invariant_cache(RS_type::AbstractRamanType, iBand, model)
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh, greek_cabannes) = model
    arr_type = CoreRT.array_type(model)
    pol_type = CoreRT.polarization_type(model)
    FT       = CoreRT.float_type(model)

    μ_host = collect(model.quad_points.phase_qp_μ)
    nAero  = size(τ_aer[iBand[1]], 1)
    nZ     = size(τ_rayl[1], 2)
    nBands = length(iBand)

    rayl_τ_dev      = Vector{Vector}(undef, nBands)
    rayl_ϖ_Cabannes = Vector(undef, nBands)
    aer_τ_mod       = Vector{Vector{Vector}}(undef, nBands)
    aer_ϖ_mod       = Vector{Vector}(undef, nBands)
    τ_abs_dev       = Vector{Vector}(undef, nBands)
    fScattRayleigh  = Vector{Vector{Vector}}(undef, nBands)
    mode_layers     = Vector{Union{Nothing, Vector{Int}}}(undef, nBands)

    # Z-moment precompute: P/R/T depend only on (μ, l), Π only on (l, m) —
    # neither on the scatterer — so one table at the global l_max (longest
    # Greek expansion over Rayleigh + every aerosol mode, all bands) serves
    # every compute_Z_moments call of this solve.
    z_tables = nothing
    if _Z_TABLES_ENABLED[]
        l_max_global = 0
        for iB in iBand
            gr_source = _rayleigh_greek_source(RS_type, greek_rayleigh, greek_cabannes)
            gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source
            l_max_global = max(l_max_global, length(gr.β))
            for i in 1:size(τ_aer[iB], 1)
                l_max_global = max(l_max_global,
                                   length(aerosol_optics[iB][i].greek_coefs.β))
            end
        end
        # Build on the grid the diffuse Z blocks are evaluated on: qp_μ under
        # external solar (μ₀ handled separately via compute_Z_source_moments),
        # phase_qp_μ otherwise. _compute_phase_blocks re-checks the grid match.
        μ_tab = model.quad_points.external_solar ?
                collect(model.quad_points.qp_μ) : μ_host
        z_tables = Scattering.ZMomentTables(μ_tab, l_max_global)
    end
    Π_m0 = z_tables === nothing ? nothing :
           Scattering.make_Π_lists(pol_type, z_tables, 0)

    for (iBi, iB) in enumerate(iBand)
        gr_source = _rayleigh_greek_source(RS_type, greek_rayleigh, greek_cabannes)
        gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source

        # -- Upload m-independent τ/ϖ per layer --
        rayl_τ_dev_iB      = [_to_device(arr_type, τ_rayl[iB][:,iz]) for iz=1:nZ]
        rayl_ϖ_Cabannes_iB = RS_type.ϖ_Cabannes[iB]  # scalar or existing device array

        # Aerosol: pre-apply δ-M correction.
        # τ_aer[iB] is now 3-D [iAer, nSpec, iLayer], so τ_aer[iB][i,:,iz] is an
        # nSpec vector. ω̃ and fᵗ may be scalars (single-λ band) or vectors (multi-λ band);
        # broadcasting handles both cases.
        # aer_τ_mod[iBi][iAer][iz] → nSpec vector (or scalar for single-λ bands)
        # aer_ϖ_mod[iBi][iAer]     → nSpec vector (or scalar for single-λ bands)
        aer_τ_mod_iB = Vector{Vector}(undef, nAero)
        aer_ϖ_mod_iB = Vector(undef, nAero)
        for i = 1:nAero
            (; fᵗ, ω̃) = aerosol_optics[iB][i]
            # Device-type ω̃/fᵗ when they are host Vectors (multi-λ band on GPU):
            # aer_ϖ_mod is later passed into CoreScatteringOpticalProperties on the
            # device path (cache overload), so it must live on the device.
            # Scalars (single-λ) are unchanged — _to_device only acts on AbstractArray.
            ω̃_d  = ω̃  isa AbstractArray ? _to_device(arr_type, ω̃)  : ω̃
            fᵗ_d = fᵗ isa AbstractArray ? _to_device(arr_type, fᵗ) : fᵗ
            aer_ϖ_mod_iB[i] = (1 .- fᵗ_d) .* ω̃_d ./ (1 .- fᵗ_d .* ω̃_d)
            aer_τ_mod_iB[i] = [_to_device(arr_type, (1 .- fᵗ .* ω̃) .* τ_aer[iB][i, :, iz]) for iz=1:nZ]
        end

        τ_abs_dev_iB = [_to_device(arr_type, τ_abs[iB][:,iz]) for iz=1:nZ]

        # Compute fScattRayleigh once using m=0 Z moments.
        # τ/ϖ mixing in the `+` operator is Z-independent, so this is valid for all m.
        Rayl𝐙⁺⁺_m0, Rayl𝐙⁻⁺_m0, RaylZ₀⁺_m0, RaylZ₀⁻_m0 =
            _compute_phase_blocks(model, gr, 0, arr_type;
                                  tables=z_tables, Π_pair=Π_m0)
        rayl_m0 = [CoreScatteringOpticalProperties(rayl_τ_dev_iB[iz],
                    rayl_ϖ_Cabannes_iB, Rayl𝐙⁺⁺_m0, Rayl𝐙⁻⁺_m0,
                    RaylZ₀⁺_m0, RaylZ₀⁻_m0) for iz=1:nZ]
        ml = _SPARSE_AERO_COMBINE_ENABLED[] ? _detect_mode_layers(τ_aer[iB]) : nothing
        mode_layers[iBi] = ml
        combo_m0 = Vector{CoreScatteringOpticalProperties}(rayl_m0)
        for i = 1:nAero
            ml !== nothing && ml[i] == 0 && continue   # all-zero mode
            AerZ⁺⁺_m0, AerZ⁻⁺_m0, AerZ₀⁺_m0, AerZ₀⁻_m0 =
                _compute_aerosol_phase_blocks(model, aerosol_optics[iB][i],
                                              get_spec_bands(model)[iB], 0, arr_type;
                                              tables=z_tables, Π_pair=Π_m0)
            if ml === nothing
                # Dense combine: modes may overlap layers.
                aer_m0 = [CoreScatteringOpticalProperties(aer_τ_mod_iB[i][iz],
                             aer_ϖ_mod_iB[i], AerZ⁺⁺_m0, AerZ⁻⁺_m0,
                             AerZ₀⁺_m0, AerZ₀⁻_m0) for iz=1:nZ]
                combo_m0 = combo_m0 .+ aer_m0
            else
                # Sparse combine: layer-resolved diagonal τ — mode i lives ONLY
                # in layer ml[i]; adding its zero τ to every other layer is a
                # no-op the dense path pays O(nZ) device broadcasts for.
                iz = ml[i]
                combo_m0[iz] = combo_m0[iz] + CoreScatteringOpticalProperties(
                    aer_τ_mod_iB[i][iz], aer_ϖ_mod_iB[i], AerZ⁺⁺_m0, AerZ⁻⁺_m0,
                    AerZ₀⁺_m0, AerZ₀⁻_m0)
            end
        end
        # Device→host: needed for _expand_layer_rayleigh! (Raman path). Done once.
        # Denominator is the TOTAL extinction (incl. gas absorption), matching the
        # non-cache constructCoreOpticalProperties path: the RRS elemental kernels
        # attenuate with the total dτ_λ. Omitting τ_abs here (the gas-free combo_m0)
        # over-weighted Raman source terms inside absorption lines — the LuT bug.
        combo_m0_with_absorption =
            combo_m0 .+ [CoreAbsorptionOpticalProperties(τ_abs_dev_iB[iz]) for iz=1:nZ]
        fScattRayleigh_iB =
            [_rayleigh_fraction_of_total_extinction(rayl_m0[iz], combo_m0_with_absorption[iz])
             for iz=1:nZ]

        rayl_τ_dev[iBi]      = rayl_τ_dev_iB
        rayl_ϖ_Cabannes[iBi] = rayl_ϖ_Cabannes_iB
        aer_τ_mod[iBi]       = aer_τ_mod_iB
        aer_ϖ_mod[iBi]       = aer_ϖ_mod_iB
        τ_abs_dev[iBi]       = τ_abs_dev_iB
        fScattRayleigh[iBi]  = fScattRayleigh_iB
    end

    return MInvariantCache{FT}(μ_host, rayl_τ_dev, rayl_ϖ_Cabannes,
                                aer_τ_mod, aer_ϖ_mod, τ_abs_dev, fScattRayleigh,
                               mode_layers, z_tables)
end

function constructCoreOpticalProperties(RS_type::AbstractRamanType, iBand, m, model)
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh, greek_cabannes) = model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"

    arr_type = CoreRT.array_type(model)

    pol_type = CoreRT.polarization_type(model)

    # Quadrature points:
    # Phase matrices may carry one extra solar column/row beyond the square
    # diffuse operator. The elemental R/T kernel consumes the leading diffuse
    # block; the SFI source kernel consumes the exact μ₀ phase column.
    μ = collect(model.quad_points.phase_qp_μ)
    N = length(model.quad_points.qp_μN)
    # Number of Aerosols:
    nAero = size(τ_aer[iBand[1]],1)
    nZ    = size(τ_rayl[1],2)
    # Rayleigh Z matrix per band — pick greek_rayleigh for noRS (pure-elastic,
    # includes rotational Raman in the effective depol) or greek_cabannes for
    # any Raman-aware RS_type (rotational Raman is handled separately via the
    # RRS interaction, so the elastic path should use the lower-depol Cabannes
    # phase matrix). Mismatch here produces a ~1% bias on Stokes I and ~3% on Q
    # because the Cabannes depol (~0.007) vs Rayleigh depol (~0.028) shifts the
    # polarization-sensitive greek coefficients (β, δ) by ~3%.
    # See sanghavi src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:31-42.
    band_layer_props    = Vector{Vector{CoreScatteringOpticalProperties}}()
    band_fScattRayleigh = Vector{Vector}()
    for iB in iBand
        gr_source = _rayleigh_greek_source(RS_type, greek_rayleigh, greek_cabannes)
        gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻ =
            _compute_phase_blocks(model, gr, m, arr_type)
        rayl = [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), RS_type.ϖ_Cabannes[iB],
                Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻) for i=1:nZ]

        combo = rayl

        for i=1:nAero
            AerZ⁺⁺, AerZ⁻⁺, AerZ₀⁺, AerZ₀⁻ = _compute_aerosol_phase_blocks(
                model, aerosol_optics[iB][i], get_spec_bands(model)[iB], m, arr_type)
            # τ_aer[iB] is now 3-D [iAer, nSpec, iLayer]; extract per-spectral
            # τ vector for layer iz matching the Rayleigh pattern: arr_type(τ_rayl[iB][:,iz]).
            aer = [createAero(arr_type(τ_aer[iB][i,:,iz]),
                               aerosol_optics[iB][i],
                               AerZ⁺⁺, AerZ⁻⁺, AerZ₀⁺, AerZ₀⁻) for iz=1:nZ]
            combo = combo .+ aer
        end

        combo_with_absorption =
            combo .+ [CoreAbsorptionOpticalProperties(arr_type((τ_abs[iB][:,i]))) for i=1:nZ]

        fScattRayleigh =
            [_rayleigh_fraction_of_total_extinction(rayl[i], combo_with_absorption[i])
             for i=1:nZ]
        push!(band_layer_props, combo_with_absorption)
        push!(band_fScattRayleigh, fScattRayleigh)
    end

    layer_opt = [prod([band_layer_props[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    fscat_opt = [[band_fScattRayleigh[i][iz] for i=1:length(iBand)] for iz=1:nZ]
    return layer_opt, fscat_opt
end

"""
    constructCoreOpticalProperties(RS_type, iBand, m, model; cache=nothing)

Cache-aware overload: when `cache::MInvariantCache` is supplied, skip all
m-independent work (τ/ϖ device uploads, `collect(qp_μ)`, fScattRayleigh
computation) and only compute the Fourier-moment-dependent Z matrices.

When `cache === nothing` the function falls through to the original method
that recomputes everything — this preserves bit-exact backward compatibility
for callers that do not supply a cache (the lin path, rt_run_ss, etc.).
"""
function constructCoreOpticalProperties(RS_type::AbstractRamanType, iBand, m, model,
                                         cache::MInvariantCache)
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh, greek_cabannes) = model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"

    arr_type = CoreRT.array_type(model)
    pol_type = CoreRT.polarization_type(model)

    # Use cached host μ — avoids collect(qp_μ) per moment
    μ = cache.μ_host
    nAero = size(τ_aer[iBand[1]], 1)
    nZ    = size(τ_rayl[1], 2)

    band_layer_props = Vector{Vector{CoreScatteringOpticalProperties}}()

    # Π matrices for this m are scatterer-independent — assemble once from the
    # cached P/R/T tables and share across Rayleigh + every aerosol mode below.
    Π_m = cache.z_tables === nothing ? nothing :
          Scattering.make_Π_lists(pol_type, cache.z_tables, m)

    gr_source = _rayleigh_greek_source(RS_type, greek_rayleigh, greek_cabannes)
    selected_rayleigh = [gr_source isa AbstractVector ? gr_source[iB] : gr_source
                         for iB in iBand]
    # Rayleigh has no intrinsic wavelength dependence once its Greek
    # coefficients are fixed. Auto depolarization and update_model! may still
    # make those coefficients band-specific, so share only after an explicit
    # equality check.
    shared_rayleigh = length(selected_rayleigh) > 1 && all(
        g -> _same_greek_coefficients(g, selected_rayleigh[1]), selected_rayleigh)
    shared_blocks = shared_rayleigh ?
        _compute_phase_blocks(model, selected_rayleigh[1], m, arr_type;
                              tables=cache.z_tables, Π_pair=Π_m) : nothing

    for (iBi, iB) in enumerate(iBand)
        gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source

        # Compute Z moments for this m (only m-dependent work)
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺, RaylZ₀⁺, RaylZ₀⁻ = shared_blocks === nothing ?
            _compute_phase_blocks(model, gr, m, arr_type;
                                  tables=cache.z_tables, Π_pair=Π_m) : shared_blocks

        # Use cached device τ/ϖ uploads — avoids arr_type(τ_rayl[iB][:,i]) per layer
        rayl = [CoreScatteringOpticalProperties(cache.rayl_τ_dev[iBi][iz],
                    cache.rayl_ϖ_Cabannes[iBi], Rayl𝐙⁺⁺, Rayl𝐙⁻⁺,
                    RaylZ₀⁺, RaylZ₀⁻) for iz=1:nZ]

        ml = cache.mode_layers[iBi]
        combo = Vector{CoreScatteringOpticalProperties}(rayl)
        for i = 1:nAero
            ml !== nothing && ml[i] == 0 && continue   # all-zero mode
            AerZ⁺⁺, AerZ⁻⁺, AerZ₀⁺, AerZ₀⁻ = _compute_aerosol_phase_blocks(
                model, aerosol_optics[iB][i], get_spec_bands(model)[iB], m, arr_type;
                tables=cache.z_tables, Π_pair=Π_m)
            # Use cached pre-corrected δ-M values — avoids createAero recomputing fᵗ/ω̃.
            # aer_τ_mod[iBi][i][iz] is an nSpec device vector (or scalar for single-λ bands);
            # aer_ϖ_mod[iBi][i] is an nSpec vector (or scalar); both broadcast in CoreScatteringOpticalProperties.
            if ml === nothing
                # Dense combine (overlapping modes): every layer gets every mode.
                aer = [CoreScatteringOpticalProperties(cache.aer_τ_mod[iBi][i][iz],
                           cache.aer_ϖ_mod[iBi][i], AerZ⁺⁺, AerZ⁻⁺,
                           AerZ₀⁺, AerZ₀⁻) for iz=1:nZ]
                combo = combo .+ aer
            else
                # Sparse combine (layer-resolved diagonal τ): mode i occupies ONLY
                # layer ml[i]; the dense path would add its zero τ to the other
                # nZ−1 layers — O(nAero × nZ) wasted device broadcasts per moment
                # (measured 83% of the whole solve at nAero = nZ = 72).
                iz = ml[i]
                combo[iz] = combo[iz] + CoreScatteringOpticalProperties(
                    cache.aer_τ_mod[iBi][i][iz], cache.aer_ϖ_mod[iBi][i],
                    AerZ⁺⁺, AerZ⁻⁺, AerZ₀⁺, AerZ₀⁻)
            end
        end

        # Add cached device τ_abs uploads — avoids arr_type(τ_abs[iB][:,i]) per layer
        push!(band_layer_props,
              combo .+ [CoreAbsorptionOpticalProperties(cache.τ_abs_dev[iBi][iz]) for iz=1:nZ])
    end

    layer_opt = [prod([band_layer_props[i][iz] for i=1:length(iBand)]) for iz=1:nZ]
    # Return the pre-computed (cached) fScattRayleigh — avoids the per-m
    # Array(rayl.τ ./ combo.τ) device→host download.
    fscat_opt = [[cache.fScattRayleigh[iBi][iz] for iBi=1:length(iBand)] for iz=1:nZ]
    return layer_opt, fscat_opt
end

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺,
                    AerZ₀⁺=nothing, AerZ₀⁻=nothing)
    (; fᵗ, ω̃) = aerosol_optics
    # Device-type ω̃ and fᵗ to match τAer when they are host Vectors (multi-λ band
    # on GPU: τAer is already a device array, but ω̃/fᵗ remain host Vector{FT}).
    # Scalars (single-λ band) are left unchanged — `oftype` only acts on AbstractArray.
    ω̃d  = ω̃  isa AbstractArray ? oftype(τAer, ω̃)  : ω̃
    fᵗd = fᵗ isa AbstractArray ? oftype(τAer, fᵗ) : fᵗ
    τ_mod = (1 .- fᵗd .* ω̃d) .* τAer
    ϖ_mod = (1 .- fᵗd) .* ω̃d ./ (1 .- fᵗd .* ω̃d)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod, AerZ⁺⁺, AerZ⁻⁺,
                                    AerZ₀⁺, AerZ₀⁻)
end

# Extract scattering definitions and integrated absorptions for the source function!
function extractEffectiveProps(
            lods::Array,#{CoreScatteringOpticalProperties{FT},1}
            quad_points::QuadPoints{FT}
            ) where FT

    nSpec = length(lods[1].τ)
    nZ    = length(lods)
    scattering_interface = ScatteringInterface_00()
    scattering_interfaces_all = AbstractScatteringInterface[]
    τ_sum_all = similar(lods[1].τ, (nSpec, nZ+1))
    τ_sum_all[:,1] .= 0
    for iz = 1:nZ
        scatter = maximum(lods[iz].τ .* lods[iz].ϖ) > 2eps(FT)
        scattering_interface = get_scattering_interface(scattering_interface, scatter, iz)
        push!(scattering_interfaces_all, scattering_interface)
        @views τ_sum_all[:,iz+1] = τ_sum_all[:,iz] + getG_atSun(lods[iz], quad_points) * lods[iz].τ 
    end
    return scattering_interfaces_all, τ_sum_all
end

function getG_atSun(lod::CoreScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    FT(1)
end

function getG_atSun(lod::CoreDirectionalScatteringOpticalProperties,quad_points::QuadPoints{FT}) where FT
    (; iμ₀) = quad_points
    gfct = collect(lod.G)[iμ₀]
    return gfct
end


"""
    _to_device(arr_type, x)

Convert array `x` to the target device array type without copying when it is
already the correct concrete type.  This avoids the `Array(::Array)` copy that
the plain `arr_type(x)` call would trigger on CPU.
"""
@inline function _to_device(arr_type, x::AbstractArray)
    x isa arr_type ? x : arr_type(x)
end

"""
    _ensure_3d(Z)

Reshape a 2-D phase matrix `(NquadN, NquadN)` to 3-D `(NquadN, NquadN, 1)`
so that the elemental kernels — which index `Z[i,j,n2]` with `n2=1` when
`size(Z,3)==1` — always receive a 3-D array regardless of how the Z was built.
When `Z` is already 3-D the array is returned unchanged (no allocation).
"""
@inline _ensure_3d(Z::AbstractArray{T,2}) where {T} = reshape(Z, size(Z,1), size(Z,2), 1)
@inline _ensure_3d(Z::AbstractArray{T,3}) where {T} = Z

function expandOpticalProperties(in::CoreScatteringOpticalProperties, arr_type;
                                  expand_Z::Bool = false)
    (; τ, ϖ, Z⁺⁺, Z⁻⁺, Z₀⁺, Z₀⁻) = in
    solar_plus = Z₀⁺ === nothing ? nothing : _to_device(arr_type, _ensure_3d(Z₀⁺))
    solar_minus = Z₀⁻ === nothing ? nothing : _to_device(arr_type, _ensure_3d(Z₀⁻))
    if expand_Z && solar_plus !== nothing && size(solar_plus, 3) == 1
        solar_plus = _repeat(solar_plus, 1, 1, length(τ))
        solar_minus = _repeat(solar_minus, 1, 1, length(τ))
    end
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺, 3) == 1
        if expand_Z
            # Caller needs fully-expanded (NquadN, NquadN, nSpec) arrays,
            # e.g. Base.:* which concatenates multiple layers along dim=3.
            Z⁺⁺ = _repeat(Z⁺⁺, 1, 1, length(τ))
            Z⁻⁺ = _repeat(Z⁻⁺, 1, 1, length(τ))
            return CoreScatteringOpticalProperties(
                _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺),
                solar_plus, solar_minus)
        else
            # Fast path: all elemental/doubling kernels branch on size(Z,3)
            # and use index n2=1 when size==1, so no replication is needed.
            # Ensure Z is at least 3-D (reshape from 2-D if necessary) so
            # that Z[i,j,1] is always valid inside kernel code.
            Z3⁺⁺ = _to_device(arr_type, _ensure_3d(Z⁺⁺))
            Z3⁻⁺ = _to_device(arr_type, _ensure_3d(Z⁻⁺))
            return CoreScatteringOpticalProperties(
                _to_device(arr_type, τ), _to_device(arr_type, ϖ), Z3⁺⁺, Z3⁻⁺,
                solar_plus, solar_minus)
        end
    else
        @assert size(Z⁺⁺, 3) == length(τ) "Z and τ dimensions need to match"
        return CoreScatteringOpticalProperties(
            _to_device(arr_type, τ), _to_device(arr_type, ϖ),
            _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺),
            solar_plus, solar_minus)
    end
end

function expandOpticalProperties(in::CoreDirectionalScatteringOpticalProperties, arr_type;
                                  expand_Z::Bool = false)
    (; τ, ϖ, Z⁺⁺, Z⁻⁺, G) = in
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺, 3) == 1
        if expand_Z
            Z⁺⁺ = _repeat(Z⁺⁺, 1, 1, length(τ))
            Z⁻⁺ = _repeat(Z⁻⁺, 1, 1, length(τ))
            return CoreDirectionalScatteringOpticalProperties(
                _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺),
                _to_device(arr_type, G))
        else
            Z3⁺⁺ = _to_device(arr_type, _ensure_3d(Z⁺⁺))
            Z3⁻⁺ = _to_device(arr_type, _ensure_3d(Z⁻⁺))
            return CoreDirectionalScatteringOpticalProperties(
                _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                Z3⁺⁺, Z3⁻⁺, _to_device(arr_type, G))
        end
    else
        @assert size(Z⁺⁺, 3) == length(τ) "Z and τ dimensions need to match"
        return CoreDirectionalScatteringOpticalProperties(
            _to_device(arr_type, τ), _to_device(arr_type, ϖ),
            _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺),
            _to_device(arr_type, G))
    end
end

function expandBandScalars(RS_type, x)
    # `constructCoreOpticalProperties` returns one `x` entry per requested band,
    # in requested-band order. `RS_type.iBand` stores model band numbers, which
    # may be non-contiguous for single-band calls such as `iBand = [2]`.
    out = zeros(eltype(x[1]), sum(length.(RS_type.bandSpecLim)))
    for (ib_local, _) in enumerate(RS_type.iBand)
        band_range = RS_type.bandSpecLim[ib_local]
        out[band_range] .= expandScalar(x[ib_local], length(band_range))
    end
    return out
end

expandScalar(x,n) = x.*ones(n);
