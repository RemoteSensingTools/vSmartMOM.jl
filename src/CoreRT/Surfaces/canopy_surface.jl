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
"""
struct CanopyCache{FT, AZ, AG, AL<:AddedLayer{FT}, CL<:CompositeLayer{FT}, IS}
    Zup::AZ
    Zdown::AZ
    G::AG
    canopy_added::AL
    canopy_composite::CL
    soil_added::AL
    I_static::IS
end

"""
    _init_canopy_cache!(canopy, added_layer, pol_type, quad_points, architecture)

Compute azimuthal pre-integration (`Zup`, `Zdown`), the Ross `G` factor,
and allocate working arrays for the canopy RT.  Stored in `canopy._cache`.
"""
function _init_canopy_cache!(canopy::CanopySurface{FT},
                             added_layer, pol_type,
                             quad_points, architecture) where {FT}
    arr_type = array_type(architecture)
    (; qp_μN) = quad_points

    Zup, Zdown = CanopyOptics.precompute_Zazi(
        canopy.canopy_scattering, collect(qp_μN), canopy.LAD)
    G = arr_type(CanopyOptics.G(collect(qp_μN), canopy.LAD))

    dims = (size(added_layer.r⁻⁺, 1), size(added_layer.r⁻⁺, 2))
    nSpec = size(added_layer.r⁻⁺, 3)

    rs_internal = noRS{FT}()
    canopy_added     = make_added_layer(rs_internal, FT, arr_type, dims, nSpec)
    canopy_composite = make_composite_layer(rs_internal, FT, arr_type, dims, nSpec)
    soil_added       = make_added_layer(rs_internal, FT, arr_type, dims, nSpec)
    I_static = Diagonal(arr_type(Diagonal{FT}(ones(dims[1]))))

    canopy._cache = CanopyCache(Zup, Zdown, G,
                                canopy_added, canopy_composite,
                                soil_added, I_static)
    return nothing
end

"""
    create_surface_layer!(canopy::CanopySurface, added_layer, SFI, m, pol_type,
                          quad_points, τ_sum, architecture)

Compute the effective canopy+soil surface layer for Fourier moment `m`.

Internally:
1. Compute Z-matrices from cached azimuthal pre-integration.
2. For each canopy sub-layer: build `CoreDirectionalScatteringOpticalProperties`,
   run `elemental!` (canopy kernel) + `doubling!`.
3. Accumulate sub-layers via `interaction!`.
4. Create soil surface layer and interact with the canopy composite.
5. Copy the resulting R⁻⁺ and J₀ into the output `added_layer`.
"""
function create_surface_layer!(canopy::CanopySurface{FT},
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI, m::Int, pol_type,
                               quad_points, τ_sum,
                               architecture) where {FT}
    arr_type = array_type(architecture)
    (; qp_μN, iμ₀) = quad_points

    if canopy._cache === nothing
        _init_canopy_cache!(canopy, added_layer, pol_type, quad_points, architecture)
    end
    cache = canopy._cache::CanopyCache

    nSpec = size(added_layer.r⁻⁺, 3)
    Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n

    𝐙⁺⁺, 𝐙⁻⁺ = CanopyOptics.compute_Z_matrices_aniso(
        canopy.canopy_scattering, collect(qp_μN), canopy.LAD,
        cache.Zup, cache.Zdown, m)

    G_iμ₀ = Array(cache.G)[iμ₀]

    ϖ_canopy = _canopy_ssa(canopy, FT)

    τ_sum_curr = Array(τ_sum)

    for ilay in 1:canopy.n_layers
        LAI_sub = _sublayer_lai(canopy, ilay, FT)

        canopy_core = CoreDirectionalScatteringOpticalProperties(
            arr_type(LAI_sub * ones(FT, nSpec)),
            arr_type(ϖ_canopy * ones(FT, nSpec)),
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

@inline function _canopy_ssa(canopy::CanopySurface{FT}, ::Type{FT}) where {FT}
    if canopy.leaf_reflectance isa Number && canopy.leaf_transmittance isa Number
        return FT(canopy.leaf_reflectance + canopy.leaf_transmittance)
    else
        lr = canopy.leaf_reflectance isa Number ? canopy.leaf_reflectance : sum(canopy.leaf_reflectance)/length(canopy.leaf_reflectance)
        lt = canopy.leaf_transmittance isa Number ? canopy.leaf_transmittance : sum(canopy.leaf_transmittance)/length(canopy.leaf_transmittance)
        return FT(lr + lt)
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
    _interleave_atm_layer!(canopy, ilay, nSpec, τ_sum_curr, SFI,
                           pol_type, quad_points, cache, arr_type, architecture)

Insert a thin atmospheric absorption layer between canopy sub-layers `ilay` and `ilay+1`.
Uses the pre-set `canopy._within_canopy_τ` (absorption-only, no scattering).
"""
function _interleave_atm_layer!(canopy::CanopySurface{FT}, ilay, nSpec, τ_sum_curr, SFI,
                                pol_type, quad_points, cache, arr_type, architecture) where {FT}
    τ_atm = canopy._within_canopy_τ
    if τ_atm === nothing || length(τ_atm) != nSpec
        return nothing
    end
    τ_gap = arr_type(τ_atm ./ max(1, canopy.n_layers - 1))
    atm_abs = CoreAbsorptionOpticalProperties(τ_gap)

    τ_sum_dev = arr_type(τ_sum_curr)
    (; qp_μN, μ₀) = quad_points
    expk_abs = arr_type(exp.(-Array(τ_gap) ./ μ₀))
    ndoubl_abs = 0

    elemental!(pol_type, SFI, τ_sum_dev, τ_gap, atm_abs,
               0, ndoubl_abs, false, quad_points, cache.canopy_added, architecture)
    interaction!(ScatteringInterface_10(), SFI,
                 cache.canopy_composite, cache.canopy_added, cache.I_static)

    @inbounds @. τ_sum_curr += Array(τ_gap)
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
                      canopy.canopy_scattering, canopy.leaf_reflectance,
                      canopy.leaf_transmittance, canopy.include_atm,
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
