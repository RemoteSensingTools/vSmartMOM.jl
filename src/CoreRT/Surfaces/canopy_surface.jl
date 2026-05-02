#=

Canopy-coupled surface: create_surface_layer! dispatch for CanopySurface.

Internally runs canopy sub-layers through the adding-doubling method, 
combines with the soil BRDF, and presents the effective canopy+soil 
reflectance as a surface boundary condition to the atmospheric RT.

=#

"""
    CanopyCache{FT}

Pre-allocated workspace for canopy RT inside `create_surface_layer!`.
Created lazily on first call and stored in `CanopySurface._cache`.

The cache also holds canopy Z matrices for every Fourier moment. For scalar
leaf optics this is a single spectral node; for spectral leaf optics these are
coarse-grid Z matrices used for spectral interpolation.
"""
struct CanopyCache{FT, AG, AL<:AddedLayer{FT}, CL<:CompositeLayer{FT}, IS}
    G::AG
    canopy_added::AL
    canopy_composite::CL
    soil_added::AL
    I_static::IS
    "Interpolated leaf reflectance on computation grid (nothing if scalar)"
    R_interp::Union{Nothing, Vector{FT}}
    "Interpolated leaf transmittance on computation grid (nothing if scalar)"
    T_interp::Union{Nothing, Vector{FT}}
    "Coarse wavenumber grid for Z-matrix interpolation (nothing if scalar)"
    Z_coarse_wn::Union{Nothing, Vector{FT}}
    "Z⁺⁺ on coarse grid [Nq, Nq, n_coarse, n_m] for each Fourier moment"
    Z_coarse_pp::Union{Nothing, Array{FT,4}}
    "Z⁻⁺ on coarse grid [Nq, Nq, n_coarse, n_m] for each Fourier moment"
    Z_coarse_mp::Union{Nothing, Array{FT,4}}
end

@inline _array_with_eltype(::Type{FT}, A::Array{FT,N}) where {FT,N} = A
@inline _array_with_eltype(::Type{FT}, A::AbstractArray{<:Any,N}) where {FT,N} =
    Array{FT,N}(A)

function _compute_canopy_Z_stack(canopy::CanopySurface{FT},
                                 scat_model, qp_μN,
                                 max_m::Int,
                                 ::Type{FT}) where {FT}
    n_m = max(max_m, 1)
    Zpp, Zmp = CanopyOptics.compute_Z_matrices_aniso_analytic(
        scat_model, qp_μN, canopy.LAD, n_m - 1;
        quadrature = canopy.canopy_quadrature)
    return _array_with_eltype(FT, Zpp), _array_with_eltype(FT, Zmp)
end

"""
    CanopySurface_from_prospect(leaf_prospect, prospect_grid; soil, LAI, kwargs...)

Convenience constructor: compute spectral leaf R/T from PROSPECT and build
a `CanopySurface` with the spectral leaf optics on the PROSPECT wavelength grid.

# Arguments
- `leaf_prospect`: a `CanopyOptics.LeafProspectProProperties` (N, Ccab, Ccar, ...)
- `prospect_grid`: wavelength range in nm (default 400:1:2500)
- Remaining keyword arguments are forwarded to `CanopySurface(...)`.
"""
function CanopySurface_from_prospect(leaf_prospect, prospect_grid=400.0:1.0:2500.0;
                                     soil::AbstractSurfaceType,
                                     LAI, n_layers::Int=1, LAD=nothing,
                                     canopy_quadrature=nothing,
                                     canopy_clumping=nothing,
                                     include_atm::Bool=false,
                                     canopy_dp=nothing,
                                     lai_fractions=nothing)
    opti = CanopyOptics.createLeafOpticalStruct(prospect_grid)
    R, T = CanopyOptics.prospect(leaf_prospect, opti)
    λ_nm = [Float64(v.val) for v in opti.λ]
    CanopySurface(; soil=soil, LAI=LAI, n_layers=n_layers, LAD=LAD,
                   leaf_reflectance=R, leaf_transmittance=T,
                   leaf_optics_grid=λ_nm, grid_unit=:nm,
                   canopy_quadrature=canopy_quadrature,
                   canopy_clumping=canopy_clumping,
                   include_atm=include_atm, canopy_dp=canopy_dp,
                   lai_fractions=lai_fractions)
end

"""
    _init_canopy_cache!(canopy, added_layer, pol_type, quad_points, architecture;
                        spec_bands_wn=nothing, max_m=3)

Compute the effective Ross `G` factor, canopy Z-matrix cache, and working
arrays for the canopy RT. Stored in `canopy._cache`.

When spectral leaf optics are provided (`canopy.leaf_optics_grid !== nothing`),
also interpolates leaf R/T to the computation grid and pre-computes Z-matrices
on a coarse spectral grid for per-wavenumber interpolation.
"""
function _init_canopy_cache!(canopy::CanopySurface{FT},
                             added_layer, pol_type,
                             quad_points, architecture;
                             spec_bands_wn::Union{Nothing, Vector{FT}}=nothing,
                             max_m::Int=3) where {FT}
    arr_type = array_type(architecture)
    (; qp_μN) = quad_points

    qμ = collect(qp_μN)
    G_raw = CanopyOptics.G(qμ, canopy.LAD)
    G = arr_type(FT.(CanopyOptics.effective_G(canopy.canopy_clumping, G_raw, qμ)))

    dims = (size(added_layer.r⁻⁺, 1), size(added_layer.r⁻⁺, 2))
    nSpec = size(added_layer.r⁻⁺, 3)

    rs_internal = noRS{FT}()
    canopy_added     = make_added_layer(rs_internal, FT, arr_type, dims, nSpec)
    canopy_composite = make_composite_layer(rs_internal, FT, arr_type, dims, nSpec)
    soil_added       = make_added_layer(rs_internal, FT, arr_type, dims, nSpec)
    I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))))

    R_interp = nothing
    T_interp = nothing
    Z_coarse_wn = nothing
    Z_coarse_pp = nothing
    Z_coarse_mp = nothing

    if canopy.leaf_optics_grid !== nothing && spec_bands_wn !== nothing
        R_interp, T_interp, Z_coarse_wn, Z_coarse_pp, Z_coarse_mp =
            _build_spectral_canopy_cache(canopy, qμ, collect(quad_points.wt_μN),
                                         spec_bands_wn, max_m, FT)
    else
        Z_stack_pp, Z_stack_mp =
            _compute_canopy_Z_stack(canopy, canopy.canopy_scattering, qμ, max_m, FT)
        Z_coarse_pp = reshape(Z_stack_pp, size(Z_stack_pp, 1), size(Z_stack_pp, 2), 1, size(Z_stack_pp, 3))
        Z_coarse_mp = reshape(Z_stack_mp, size(Z_stack_mp, 1), size(Z_stack_mp, 2), 1, size(Z_stack_mp, 3))
        _check_Z_flux_conservation(Z_coarse_pp, Z_coarse_mp, collect(quad_points.wt_μN), 1)
    end

    canopy._cache = CanopyCache(G,
                                canopy_added, canopy_composite,
                                soil_added, I_static,
                                R_interp, T_interp,
                                Z_coarse_wn, Z_coarse_pp, Z_coarse_mp)
    return nothing
end

"""
    _leaf_grid_to_wn(grid, grid_unit) -> Vector{FT}

Convert leaf optics grid to wavenumber (cm⁻¹) using Unitful.
Does NOT sort; call `sortperm` on the result if you need a sorted order.
"""
function _leaf_grid_to_wn(grid::Vector{FT}, grid_unit::Symbol) where {FT}
    if grid_unit == :cm_inv
        return copy(grid)
    else
        u_src = grid_unit == :nm ? u"nm" : error("Unknown grid_unit: $grid_unit")
        qty   = grid .* u_src
        wn    = ustrip.(uconvert.(u"cm^-1", qty, Spectral()))
        return convert(Vector{FT}, wn)
    end
end

"""
    _check_Z_flux_conservation(Z_pp, Z_mp, qp_μN, wt_μN, n_coarse; tol=0.05)

Verify the m=0 single-scattering normalisation of the canopy Z matrices
against the convention used by `Scattering.compute_Z_moments` for the
atmospheric path (Rayleigh / Mie). Both that routine and the canonical
`CanopyOptics.compute_Z_matrices` (line ~178 in leafAngleRoutines.jl)
divide Γ by `ϖ·G`, so Z is independent of leaf ω. The expected
hemispheric column integral is therefore

    Σᵢ wᵢ · (Z⁺⁺[i,j] + Z⁻⁺[i,j])  ≈  2     (for every incidence j)

`elemental!` multiplies Z by ϖ explicitly, which is what restores the
ω scaling at the layer level. If this check sees a mean of ≈ 1 instead
of 2 then the canopy Z source is missing the `/G(μ_in)` normalisation
and canopy reflectance will be biased low.
"""
function _check_Z_flux_conservation(Z_pp::Array{FT,4}, Z_mp::Array{FT,4},
                                     wt_μN, n_coarse;
                                     tol::FT=FT(5e-2)) where {FT}
    Nq = size(Z_pp, 1)
    w  = collect(wt_μN)
    target = FT(2)
    n_warn = 0
    flux_min = FT(Inf); flux_max = FT(-Inf); flux_sum = FT(0)
    for ic in 1:n_coarse, j in 1:Nq
        flux = FT(0)
        @inbounds for i in 1:Nq
            flux += w[i] * (Z_pp[i, j, ic, 1] + Z_mp[i, j, ic, 1])
        end
        flux_min = min(flux_min, flux)
        flux_max = max(flux_max, flux)
        flux_sum += flux
        if abs(flux - target) > tol * target
            n_warn += 1
        end
    end
    flux_mean = flux_sum / (n_coarse * Nq)
    if n_warn > 0
        @warn """Canopy Z-matrix m=0 flux normalisation off vs atmospheric convention.
                 Σᵢ wᵢ (Z⁺⁺ + Z⁻⁺)[:,j] expected ≈ 2, got mean=$(round(flux_mean, digits=4))
                 (range $(round(flux_min, digits=4))–$(round(flux_max, digits=4)))
                 at $n_warn / $(n_coarse * Nq) (μ', spectral) points (tol=$tol).
                 A mean ≈ 1 indicates the canopy Z source is missing the
                 /G(μ) normalisation that compute_Z_matrices applies."""
    end
    return nothing
end

"""
    _build_spectral_canopy_cache(canopy, qp_μN, spec_bands_wn, max_m, FT)

Interpolate leaf R/T from user grid to computation grid and pre-compute
Z-matrices on a coarse spectral grid for later per-wavenumber interpolation.
"""
function _build_spectral_canopy_cache(canopy::CanopySurface{FT},
                                      qp_μN, wt_μN, spec_bands_wn::Vector{FT},
                                      max_m::Int, ::Type{FT}) where {FT}
    leaf_wn = _leaf_grid_to_wn(canopy.leaf_optics_grid, canopy.grid_unit)
    leaf_R = canopy.leaf_reflectance isa Number ?
             fill(FT(canopy.leaf_reflectance), length(leaf_wn)) :
             convert(Vector{FT}, canopy.leaf_reflectance)
    leaf_T = canopy.leaf_transmittance isa Number ?
             fill(FT(canopy.leaf_transmittance), length(leaf_wn)) :
             convert(Vector{FT}, canopy.leaf_transmittance)

    sp = sortperm(leaf_wn)
    leaf_wn = leaf_wn[sp]
    leaf_R  = leaf_R[sp]
    leaf_T  = leaf_T[sp]

    itp_R = LinearInterpolation(leaf_wn, leaf_R)
    itp_T = LinearInterpolation(leaf_wn, leaf_T)

    R_interp = FT[clamp(itp_R(wn), FT(0), FT(1)) for wn in spec_bands_wn]
    T_interp = FT[clamp(itp_T(wn), FT(0), FT(1)) for wn in spec_bands_wn]

    wn_min, wn_max = extrema(spec_bands_wn)
    λ_max_nm = FT(1e7) / wn_min
    λ_min_nm = FT(1e7) / wn_max
    coarse_step_nm = FT(5)
    λ_coarse = collect(λ_min_nm:coarse_step_nm:λ_max_nm)
    if isempty(λ_coarse) || λ_coarse[end] < λ_max_nm
        push!(λ_coarse, λ_max_nm)
    end
    # LinearInterpolation needs at least 2 distinct knots; pad the coarse
    # grid when the spectral band is a single point or very narrow.
    if length(λ_coarse) < 2
        push!(λ_coarse, λ_coarse[end] + coarse_step_nm)
    end
    wn_coarse = sort(FT(1e7) ./ λ_coarse)

    R_coarse = FT[clamp(itp_R(wn), FT(0), FT(1)) for wn in wn_coarse]
    T_coarse = FT[clamp(itp_T(wn), FT(0), FT(1)) for wn in wn_coarse]

    Nq = length(qp_μN)
    n_coarse = length(wn_coarse)
    n_m = max(max_m, 1)

    Z_coarse_pp = zeros(FT, Nq, Nq, n_coarse, n_m)
    Z_coarse_mp = zeros(FT, Nq, Nq, n_coarse, n_m)

    for ic in 1:n_coarse
        scat_model = CanopyOptics.BiLambertianCanopyScattering(
            R=R_coarse[ic], T=T_coarse[ic])
        Zpp_stack, Zmp_stack = _compute_canopy_Z_stack(canopy, scat_model, qp_μN, n_m, FT)
        Z_coarse_pp[:, :, ic, :] .= Zpp_stack
        Z_coarse_mp[:, :, ic, :] .= Zmp_stack
    end

    _check_Z_flux_conservation(Z_coarse_pp, Z_coarse_mp, wt_μN, n_coarse)

    return R_interp, T_interp, wn_coarse, Z_coarse_pp, Z_coarse_mp
end

"""
    _interpolate_Z_spectral(Z_coarse_pp, Z_coarse_mp, Z_coarse_wn, spec_wn, m, FT)

Interpolate coarse-grid Z-matrices to the full computation grid for Fourier moment `m`.
Returns `(Z_pp_spec, Z_mp_spec)` as 3D arrays `[Nq, Nq, nSpec]`.
"""
function _interpolate_Z_spectral(Z_coarse_pp::Array{FT,4}, Z_coarse_mp::Array{FT,4},
                                  Z_coarse_wn::Vector{FT}, spec_wn::Vector{FT},
                                  m::Int, ::Type{FT}) where {FT}
    Nq = size(Z_coarse_pp, 1)
    nSpec = length(spec_wn)
    m_idx = m + 1
    if m_idx > size(Z_coarse_pp, 4)
        m_idx = size(Z_coarse_pp, 4)
    end

    Z_pp_out = zeros(FT, Nq, Nq, nSpec)
    Z_mp_out = zeros(FT, Nq, Nq, nSpec)

    for i in 1:Nq, j in 1:Nq
        itp_pp = LinearInterpolation(Z_coarse_wn, Z_coarse_pp[i, j, :, m_idx])
        itp_mp = LinearInterpolation(Z_coarse_wn, Z_coarse_mp[i, j, :, m_idx])
        for s in 1:nSpec
            wn = clamp(spec_wn[s], Z_coarse_wn[1], Z_coarse_wn[end])
            Z_pp_out[i, j, s] = itp_pp(wn)
            Z_mp_out[i, j, s] = itp_mp(wn)
        end
    end
    return Z_pp_out, Z_mp_out
end

"""
    create_surface_layer!(canopy::CanopySurface, added_layer, SFI, m, pol_type,
                          quad_points, τ_sum, architecture;
                          spec_bands_wn=nothing, max_m=3)

Compute the effective canopy+soil surface layer for Fourier moment `m`.

When spectral leaf optics are provided, `spec_bands_wn` (the concatenated
computation wavenumber grid) is required for initializing the spectral cache
on first call.

Internally:
1. Compute Z-matrices from the cached analytic stack (or coarse-grid interpolation).
2. For each canopy sub-layer: build `CoreDirectionalScatteringOpticalProperties`,
   run `elemental!` (canopy kernel) + `doubling!`.
3. Accumulate sub-layers via `interaction!`.
4. Create soil surface layer and interact with the canopy composite.
5. Copy the resulting R⁻⁺ and J₀ into the output `added_layer`.

# TODO: Canopy + Raman (RRS/VS) coupling is not currently implemented. The
# canopy BRDF code is noRS-only; calling rt_run with CanopySurface + RRS/VS
# has undefined behavior (the inelastic ieR/ieT branches in the kernels
# never see a canopy-adapted contribution and will return garbage or crash).
# Future work: couple the canopy path to the inelastic path — this is a
# separate project, not part of the sanghavi-unified merge.
"""
function create_surface_layer!(canopy::CanopySurface{FT},
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI, m::Int, pol_type,
                               quad_points, τ_sum,
                               architecture;
                               spec_bands_wn::Union{Nothing, Vector{FT}}=nothing,
                               max_m::Int=3) where {FT}
    arr_type = array_type(architecture)
    (; qp_μN, iμ₀) = quad_points

    if canopy._cache === nothing
        _init_canopy_cache!(canopy, added_layer, pol_type, quad_points, architecture;
                            spec_bands_wn=spec_bands_wn, max_m=max_m)
    end
    cache = canopy._cache::CanopyCache

    nSpec = size(added_layer.r⁻⁺, 3)
    Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n

    has_spectral = cache.R_interp !== nothing

    if has_spectral && cache.Z_coarse_pp !== nothing
        𝐙⁺⁺, 𝐙⁻⁺ = _interpolate_Z_spectral(cache.Z_coarse_pp, cache.Z_coarse_mp,
                                               cache.Z_coarse_wn, spec_bands_wn, m, FT)
    elseif cache.Z_coarse_pp !== nothing
        m_idx = min(m + 1, size(cache.Z_coarse_pp, 4))
        𝐙⁺⁺ = cache.Z_coarse_pp[:, :, 1, m_idx]
        𝐙⁻⁺ = cache.Z_coarse_mp[:, :, 1, m_idx]
    else
        error("Canopy Z cache was not initialized")
    end

    G_iμ₀ = Array(cache.G)[iμ₀]

    ϖ_canopy = _canopy_ssa(canopy, cache, nSpec, FT)

    τ_sum_curr = Array(τ_sum)

    for ilay in 1:canopy.n_layers
        LAI_sub = _sublayer_lai(canopy, ilay, FT)

        canopy_core = CoreDirectionalScatteringOpticalProperties(
            arr_type(LAI_sub * ones(FT, nSpec)),
            arr_type(ϖ_canopy),
            arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺), cache.G)

        τ_sum_dev = arr_type(τ_sum_curr)
        dτ, ndoubl, expk = init_layer(canopy_core, quad_points, pol_type, architecture)

        @timeit "canopy elemental" elemental!(
            pol_type, SFI, τ_sum_dev, dτ, canopy_core,
            m, ndoubl, true, quad_points, cache.canopy_added, architecture)
        @timeit "canopy doubling" doubling!(
            pol_type, SFI, expk, ndoubl, cache.canopy_added,
            cache.I_static, architecture)

        if ilay == 1
            copy_added_to_composite!(cache.canopy_composite, cache.canopy_added)
        else
            interaction!(ScatteringInterface_11(), SFI,
                         cache.canopy_composite, cache.canopy_added, cache.I_static)
        end

        @inbounds @. τ_sum_curr += G_iμ₀ * LAI_sub

        if canopy.include_atm && canopy._within_canopy_τ !== nothing && ilay < canopy.n_layers
            _interleave_atm_layer!(canopy, ilay, nSpec, τ_sum_curr, SFI,
                                   pol_type, quad_points, cache, arr_type, architecture)
        end
    end

    τ_sum_soil = arr_type(τ_sum_curr)
    create_surface_layer!(canopy.soil, cache.soil_added, SFI, m, pol_type,
                          quad_points, τ_sum_soil, architecture)

    @timeit "canopy-soil interaction" interaction!(
        ScatteringInterface_11(), SFI,
        cache.canopy_composite, cache.soil_added, cache.I_static)

    T_surf = arr_type(Diagonal(ones(FT, pol_type.n * Nquad)))
    added_layer.r⁻⁺ .= cache.canopy_composite.R⁻⁺
    added_layer.r⁺⁻ .= 0
    added_layer.t⁺⁺ .= T_surf
    added_layer.t⁻⁻ .= T_surf
    if SFI
        added_layer.j₀⁺ .= 0
        added_layer.j₀⁻ .= cache.canopy_composite.J₀⁻
    end
    return nothing
end

"""
    _canopy_ssa(canopy, cache, nSpec, FT) -> Vector{FT}

Return the canopy single-scattering albedo as a vector of length `nSpec`.
Uses interpolated spectral R/T when available, otherwise broadcasts a scalar.
"""
@inline function _canopy_ssa(canopy::CanopySurface{FT}, cache, nSpec::Int, ::Type{FT}) where {FT}
    if cache.R_interp !== nothing && cache.T_interp !== nothing
        return cache.R_interp .+ cache.T_interp
    elseif canopy.leaf_reflectance isa Number && canopy.leaf_transmittance isa Number
        return fill(FT(canopy.leaf_reflectance + canopy.leaf_transmittance), nSpec)
    else
        lr = canopy.leaf_reflectance isa Number ? canopy.leaf_reflectance :
             sum(canopy.leaf_reflectance) / length(canopy.leaf_reflectance)
        lt = canopy.leaf_transmittance isa Number ? canopy.leaf_transmittance :
             sum(canopy.leaf_transmittance) / length(canopy.leaf_transmittance)
        return fill(FT(lr + lt), nSpec)
    end
end

@inline function _sublayer_lai(canopy::CanopySurface{FT}, ilay::Int, ::Type{FT}) where {FT}
    if canopy.lai_fractions !== nothing
        return FT(canopy.LAI * canopy.lai_fractions[ilay])
    else
        return FT(canopy.LAI / canopy.n_layers)
    end
end

"""
    _compute_canopy_atm_tau!(canopy, model, spec_bands_wn)

Compute the total within-canopy atmospheric optical depth from `canopy.canopy_dp` (hPa).

Uses bottom-of-atmosphere conditions (T, p, VMR) from `model.profile` and scales
the bottom-layer absorption from the pre-computed `model.τ_abs`.

The Avogadro-based VCD calculation mirrors `compute_atmos_profile_fields`:
  `vcd = N_A × Δp / (M_dry × g × 100²) × 100`  [molec/cm²]
"""
function _compute_canopy_atm_tau!(canopy::CanopySurface{FT},
                                  model, spec_bands_wn::Vector{FT}) where {FT}
    dp = canopy.canopy_dp
    dp === nothing && return nothing

    profile = model.profile
    Nz = length(profile.p_full)

    Nₐ  = FT(6.0221408e23)   # Avogadro [mol⁻¹]
    g₀  = FT(9.8067)         # gravity [m/s²]
    dry_mass = FT(0.028964)  # dry air molar mass [kg/mol]

    vmr_h2o_bot = profile.vmr_h2o[Nz]
    vmr_dry = FT(1) - vmr_h2o_bot
    wet_mass = FT(0.018015)
    M = vmr_dry * dry_mass + vmr_h2o_bot * wet_mass
    vcd_canopy = Nₐ * FT(dp) / (M * g₀ * FT(100^2)) * FT(100)

    nSpec = length(spec_bands_wn)
    τ_canopy = zeros(FT, nSpec)

    offset = 0
    for (iB, τ_abs_band) in enumerate(model.τ_abs)
        nB = size(τ_abs_band, 1)
        τ_bot_layer = τ_abs_band[:, Nz]
        vcd_bot = profile.vcd_dry[Nz]
        if vcd_bot > FT(0)
            scale = vcd_canopy / vcd_bot
            for j in 1:nB
                @inbounds τ_canopy[offset + j] += τ_bot_layer[j] * scale
            end
        end
        offset += nB
    end

    canopy._within_canopy_τ = τ_canopy
    return nothing
end

"""
    _interleave_atm_layer!(canopy, ilay, nSpec, τ_sum_curr, SFI,
                           pol_type, quad_points, cache, arr_type, architecture)

Insert a thin atmospheric absorption layer between canopy sub-layers `ilay` and `ilay+1`.
Uses the pre-set `canopy._within_canopy_τ` (absorption-only, no scattering).
"""
function _interleave_atm_layer!(canopy::CanopySurface{FT}, ilay, nSpec, τ_sum_curr, SFI,
                                pol_type, quad_points, cache, arr_type, architecture) where {FT}
    τ_atm = canopy._within_canopy_τ
    if τ_atm === nothing || length(τ_atm) != nSpec || all(iszero, τ_atm)
        return nothing
    end
    τ_gap = arr_type(τ_atm ./ max(1, canopy.n_layers - 1))

    # Wrap the pure-absorption sub-layer as a CoreScatteringOpticalProperties
    # with ϖ = 0 and zero Z so the standard elemental! dispatch applies.
    # The ϖ = 0 multiplier zeroes out the scattering kernel exactly, leaving
    # t⁺⁺ = exp(-τ/μ) (the correct Beer–Lambert transmission) and r⁻⁺ = 0.
    NquadN = size(cache.canopy_added.r⁻⁺, 1)
    ϖ_zero = arr_type(zeros(FT, nSpec))
    Z_zero = arr_type(zeros(FT, NquadN, NquadN, nSpec))
    atm_scat = CoreScatteringOpticalProperties(τ_gap, ϖ_zero, Z_zero, Z_zero)

    # F₀ shape (nStokes, nSpec) — only its shape matters here since ϖ = 0
    # zeroes the SFI source contributions inside the elemental kernel.
    F₀_local = arr_type(zeros(FT, pol_type.n, nSpec))
    @views F₀_local[1, :] .= 1

    τ_sum_dev = arr_type(τ_sum_curr)
    ndoubl_abs = 0

    elemental!(pol_type, SFI, τ_sum_dev, τ_gap, F₀_local, atm_scat,
               0, ndoubl_abs, true, quad_points, cache.canopy_added, architecture)
    interaction!(ScatteringInterface_10(), SFI,
                 cache.canopy_composite, cache.canopy_added, cache.I_static)

    τ_gap_host = Array(τ_gap)
    @inbounds @. τ_sum_curr += τ_gap_host
    return nothing
end

"""
    create_surface_layer!(RS_type, canopy::CanopySurface, added_layer, added_layer_lin,
                          iparam, SFI, m, pol_type, quad_points, τ_sum, τ̇_sum, F₀, architecture)

Linearized surface dispatch for `CanopySurface`.  Computes the forward
canopy+soil reflectance (same as the non-linearized dispatch) and the
derivative of the soil albedo parameter at `iparam` via finite differences.

Canopy-specific parameter derivatives (LAI, leaf R/T) are not yet
implemented in the linearized path; their Jacobian slots should be populated
via the hybrid AD boundary (ForwardDiff through canopy optical properties).
"""
function create_surface_layer!(RS_type::noRS,
                               canopy::CanopySurface{FT},
                               added_layer::AddedLayer,
                               added_layer_lin,
                               iparam::Int,
                               SFI, m::Int, pol_type, quad_points,
                               τ_sum, τ̇_sum, F₀,
                               architecture) where {FT}

    create_surface_layer!(canopy, added_layer, SFI, m, pol_type,
                          quad_points, τ_sum, architecture)

    added_layer_lin.ap_ṙ⁻⁺ .= 0
    added_layer_lin.ap_ṙ⁺⁻ .= 0
    added_layer_lin.ap_ṫ⁺⁺ .= 0
    added_layer_lin.ap_ṫ⁻⁻ .= 0
    added_layer_lin.ap_J̇₀⁺ .= 0
    added_layer_lin.ap_J̇₀⁻ .= 0

    r_ref = Array(added_layer.r⁻⁺)
    j_ref = SFI ? Array(added_layer.j₀⁻) : nothing

    ε = FT(1e-7)
    perturbed = _perturb_soil_albedo(canopy, ε)
    invalidate_canopy_cache!(perturbed)
    arr_type = array_type(architecture)
    cache_p = canopy._cache
    if cache_p !== nothing
        _init_canopy_cache!(perturbed, added_layer, pol_type, quad_points, architecture)
    end

    temp_added = cache_p !== nothing ? cache_p.soil_added : 
        make_added_layer(noRS{FT}(), FT, arr_type,
            (size(added_layer.r⁻⁺,1), size(added_layer.r⁻⁺,2)),
            size(added_layer.r⁻⁺,3))
    create_surface_layer!(perturbed, temp_added, SFI, m, pol_type,
                          quad_points, τ_sum, architecture)

    ṙ = (Array(temp_added.r⁻⁺) .- r_ref) ./ ε
    added_layer_lin.ap_ṙ⁻⁺[iparam,:,:,:] .= arr_type(ṙ)

    if SFI && j_ref !== nothing
        j̇ = (Array(temp_added.j₀⁻) .- j_ref) ./ ε
        added_layer_lin.ap_J̇₀⁻[iparam,:,:,:] .= arr_type(j̇)
    end

    synchronize_if_gpu()
    return nothing
end

function _perturb_soil_albedo(canopy::CanopySurface{FT}, ε::FT) where {FT}
    soil = canopy.soil
    if soil isa LambertianSurfaceScalar
        perturbed_soil = LambertianSurfaceScalar(soil.albedo + ε)
    else
        perturbed_soil = soil
    end
    CanopySurface{FT}(perturbed_soil, canopy.LAI, canopy.n_layers, canopy.LAD,
                      canopy.canopy_scattering, canopy.canopy_quadrature,
                      canopy.canopy_clumping,
                      canopy.leaf_reflectance,
                      canopy.leaf_transmittance, canopy.leaf_optics_grid,
                      canopy.grid_unit, canopy.include_atm, canopy.canopy_dp,
                      canopy.lai_fractions, canopy._within_canopy_τ, nothing)
end

"""
    invalidate_canopy_cache!(canopy::CanopySurface)

Force recomputation of the canopy cache on the next `create_surface_layer!` call.
Use after changing canopy parameters (LAI, LAD, etc.) between runs.
"""
function invalidate_canopy_cache!(canopy::CanopySurface)
    canopy._cache = nothing
    return nothing
end
