"""
    _rayleigh_greek_source(rs, greek_rayleigh, greek_cabannes)

Select the Rayleigh phase-coefficient source for the elastic path. Pure
elastic modes use the full Rayleigh coefficients; Raman-aware modes use
Cabannes coefficients because Raman redistribution is handled explicitly.
"""
_rayleigh_greek_source(::Union{noRS, noRS_plus}, greek_rayleigh, greek_cabannes) = greek_rayleigh
_rayleigh_greek_source(::AbstractRamanType, greek_rayleigh, greek_cabannes) = greek_cabannes

"""
    MInvariantCache

Cache of m-independent optical quantities computed once before the Fourier
loop in [`rt_run`](@ref).  Passed as `cache=` to
[`constructCoreOpticalProperties`](@ref) to skip redundant uploads,
host↔device copies, and fScattRayleigh recomputation on every moment.

# Fields
- `μ_host`: host copy of quadrature cosines — avoids `collect(qp_μ)` per moment.
- `rayl_τ_dev[iBi][iz]`: device τ_rayl slice per band-index/layer — avoids re-upload.
- `rayl_ϖ_Cabannes[iBi]`: ϖ_Cabannes per band (scalar or device array) — avoids
  re-upload.
- `aer_τ_mod[iBi][iAer][iz]`: δ-M-modified aerosol τ (plain FT scalar) per
  band-index/aerosol/layer — precomputed `(1-fᵗ·ω̃)·τ_aer`.
- `aer_ϖ_mod[iBi][iAer]`: δ-M-modified aerosol ϖ (plain FT scalar) — precomputed
  `(1-fᵗ)·ω̃/(1-fᵗ·ω̃)`.
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

    μ_host = collect(model.quad_points.qp_μ)
    nAero  = size(τ_aer[iBand[1]], 1)
    nZ     = size(τ_rayl[1], 2)
    nBands = length(iBand)

    rayl_τ_dev      = Vector{Vector}(undef, nBands)
    rayl_ϖ_Cabannes = Vector(undef, nBands)
    aer_τ_mod       = Vector{Vector{Vector}}(undef, nBands)
    aer_ϖ_mod       = Vector{Vector}(undef, nBands)
    τ_abs_dev       = Vector{Vector}(undef, nBands)
    fScattRayleigh  = Vector{Vector{Vector}}(undef, nBands)

    for (iBi, iB) in enumerate(iBand)
        gr_source = _rayleigh_greek_source(RS_type, greek_rayleigh, greek_cabannes)
        gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source

        # -- Upload m-independent τ/ϖ per layer --
        rayl_τ_dev_iB      = [_to_device(arr_type, τ_rayl[iB][:,iz]) for iz=1:nZ]
        rayl_ϖ_Cabannes_iB = RS_type.ϖ_Cabannes[iB]  # scalar or existing device array

        # Aerosol: pre-apply δ-M correction (scalars — τ_aer[iB][i,iz] is a FT scalar)
        aer_τ_mod_iB = Vector{Vector}(undef, nAero)
        aer_ϖ_mod_iB = Vector(undef, nAero)
        for i = 1:nAero
            (; fᵗ, ω̃) = aerosol_optics[iB][i]
            aer_ϖ_mod_iB[i] = (1 - fᵗ) * ω̃ / (1 - fᵗ * ω̃)
            aer_τ_mod_iB[i] = [(1 - fᵗ * ω̃) * τ_aer[iB][i, iz] for iz=1:nZ]
        end

        τ_abs_dev_iB = [_to_device(arr_type, τ_abs[iB][:,iz]) for iz=1:nZ]

        # Compute fScattRayleigh once using m=0 Z moments.
        # τ/ϖ mixing in the `+` operator is Z-independent, so this is valid for all m.
        Rayl𝐙⁺⁺_m0, Rayl𝐙⁻⁺_m0 = Scattering.compute_Z_moments(pol_type, μ_host, gr, 0,
                                                                  arr_type=arr_type)
        rayl_m0 = [CoreScatteringOpticalProperties(rayl_τ_dev_iB[iz],
                    rayl_ϖ_Cabannes_iB, Rayl𝐙⁺⁺_m0, Rayl𝐙⁻⁺_m0) for iz=1:nZ]
        combo_m0 = rayl_m0
        for i = 1:nAero
            AerZ⁺⁺_m0, AerZ⁻⁺_m0 = Scattering.compute_Z_moments(
                pol_type, μ_host, aerosol_optics[iB][i].greek_coefs, 0, arr_type=arr_type)
            aer_m0 = [CoreScatteringOpticalProperties(aer_τ_mod_iB[i][iz],
                         aer_ϖ_mod_iB[i], AerZ⁺⁺_m0, AerZ⁻⁺_m0) for iz=1:nZ]
            combo_m0 = combo_m0 .+ aer_m0
        end
        # Device→host: needed for _expand_layer_rayleigh! (Raman path). Done once.
        fScattRayleigh_iB = [Array(rayl_m0[iz].τ ./ combo_m0[iz].τ) for iz=1:nZ]

        rayl_τ_dev[iBi]      = rayl_τ_dev_iB
        rayl_ϖ_Cabannes[iBi] = rayl_ϖ_Cabannes_iB
        aer_τ_mod[iBi]       = aer_τ_mod_iB
        aer_ϖ_mod[iBi]       = aer_ϖ_mod_iB
        τ_abs_dev[iBi]       = τ_abs_dev_iB
        fScattRayleigh[iBi]  = fScattRayleigh_iB
    end

    return MInvariantCache{FT}(μ_host, rayl_τ_dev, rayl_ϖ_Cabannes,
                                aer_τ_mod, aer_ϖ_mod, τ_abs_dev, fScattRayleigh)
end

function constructCoreOpticalProperties(RS_type::AbstractRamanType, iBand, m, model)
    (; τ_rayl, τ_aer, τ_abs, aerosol_optics, greek_rayleigh, greek_cabannes) = model
    @assert all(iBand .≤ length(τ_rayl)) "iBand exceeded number of bands"

    arr_type = CoreRT.array_type(model)

    pol_type = CoreRT.polarization_type(model)

    # Quadrature points:
    μ = collect(model.quad_points.qp_μ)
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
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, gr, m,
                                                        arr_type = arr_type)
        rayl = [CoreScatteringOpticalProperties(arr_type(τ_rayl[iB][:,i]), RS_type.ϖ_Cabannes[iB],
                Rayl𝐙⁺⁺, Rayl𝐙⁻⁺) for i=1:nZ]

        combo = rayl

        for i=1:nAero
            AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(
                                pol_type, μ,
                                aerosol_optics[iB][i].greek_coefs,
                                m, arr_type=arr_type)
            aer = createAero.(τ_aer[iB][i,:],
                              [aerosol_optics[iB][i]],
                              [AerZ⁺⁺], [AerZ⁻⁺])
            combo = combo .+ aer
        end

        fScattRayleigh = [Array(rayl[i].τ ./ combo[i].τ) for i=1:nZ]
        push!(band_layer_props,
              combo .+ [CoreAbsorptionOpticalProperties(arr_type((τ_abs[iB][:,i]))) for i=1:nZ])
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

    for (iBi, iB) in enumerate(iBand)
        gr_source = _rayleigh_greek_source(RS_type, greek_rayleigh, greek_cabannes)
        gr = gr_source isa AbstractVector ? gr_source[iB] : gr_source

        # Compute Z moments for this m (only m-dependent work)
        Rayl𝐙⁺⁺, Rayl𝐙⁻⁺ = Scattering.compute_Z_moments(pol_type, μ, gr, m,
                                                           arr_type=arr_type)

        # Use cached device τ/ϖ uploads — avoids arr_type(τ_rayl[iB][:,i]) per layer
        rayl = [CoreScatteringOpticalProperties(cache.rayl_τ_dev[iBi][iz],
                    cache.rayl_ϖ_Cabannes[iBi], Rayl𝐙⁺⁺, Rayl𝐙⁻⁺) for iz=1:nZ]

        combo = rayl
        for i = 1:nAero
            AerZ⁺⁺, AerZ⁻⁺ = Scattering.compute_Z_moments(
                pol_type, μ, aerosol_optics[iB][i].greek_coefs, m, arr_type=arr_type)
            # Use cached pre-corrected δ-M scalars — avoids createAero recomputing fᵗ/ω̃
            aer = [CoreScatteringOpticalProperties(cache.aer_τ_mod[iBi][i][iz],
                       cache.aer_ϖ_mod[iBi][i], AerZ⁺⁺, AerZ⁻⁺) for iz=1:nZ]
            combo = combo .+ aer
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

function createAero(τAer, aerosol_optics, AerZ⁺⁺, AerZ⁻⁺)
    (; fᵗ, ω̃) = aerosol_optics
    τ_mod = (1-fᵗ * ω̃ ) * τAer;
    ϖ_mod = (1-fᵗ) * ω̃/(1-fᵗ * ω̃)
    CoreScatteringOpticalProperties(τ_mod, ϖ_mod,AerZ⁺⁺, AerZ⁻⁺)
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
    (; τ, ϖ, Z⁺⁺, Z⁻⁺) = in
    @assert length(τ) == length(ϖ) "τ and ϖ sizes need to match"
    if size(Z⁺⁺, 3) == 1
        if expand_Z
            # Caller needs fully-expanded (NquadN, NquadN, nSpec) arrays,
            # e.g. Base.:* which concatenates multiple layers along dim=3.
            Z⁺⁺ = _repeat(Z⁺⁺, 1, 1, length(τ))
            Z⁻⁺ = _repeat(Z⁻⁺, 1, 1, length(τ))
            return CoreScatteringOpticalProperties(
                _to_device(arr_type, τ), _to_device(arr_type, ϖ),
                _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺))
        else
            # Fast path: all elemental/doubling kernels branch on size(Z,3)
            # and use index n2=1 when size==1, so no replication is needed.
            # Ensure Z is at least 3-D (reshape from 2-D if necessary) so
            # that Z[i,j,1] is always valid inside kernel code.
            Z3⁺⁺ = _to_device(arr_type, _ensure_3d(Z⁺⁺))
            Z3⁻⁺ = _to_device(arr_type, _ensure_3d(Z⁻⁺))
            return CoreScatteringOpticalProperties(
                _to_device(arr_type, τ), _to_device(arr_type, ϖ), Z3⁺⁺, Z3⁻⁺)
        end
    else
        @assert size(Z⁺⁺, 3) == length(τ) "Z and τ dimensions need to match"
        return CoreScatteringOpticalProperties(
            _to_device(arr_type, τ), _to_device(arr_type, ϖ),
            _to_device(arr_type, Z⁺⁺), _to_device(arr_type, Z⁻⁺))
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
    out = zeros(eltype(x[1]), sum([length(RS_type.bandSpecLim[iB]) for iB in RS_type.iBand]))
    for iB in RS_type.iBand
        out[RS_type.bandSpecLim[iB]] .= expandScalar(x[iB],length(RS_type.bandSpecLim[iB]))
    end
    return out
end

expandScalar(x,n) = x.*ones(n);
