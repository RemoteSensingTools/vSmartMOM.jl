#=
Azimuthal weighting for multi-sensor RT postprocessing.
Handles interlayer flux coupling before the VZA loop.
=#

"""
    _set_unscattered_downwelling!(direct_dwJ, sensor_levels, τ_sum_all,
                                   F₀, vza, qp_μ, iμ₀, μ₀, pol_type)

Populate the separately reported collimated solar beam at strict-interior
interfaces. A requested VZA is on the solar ordinate when it resolves to the
same explicit quadrature node as `μ₀`. The direct carrier follows the
existing m=0 reconstruction convention: normalization `1/(2π)`, I/Q retained,
and U/V zero. It is deliberately independent of VAZ, matching the historical
BOA direct-beam carrier and published-reference corrections.

For a boundary with `k` atmospheric layers above it, cumulative optical depth
is `τ_sum_all[:, k+1]`.
"""
function _set_unscattered_downwelling!(direct_dwJ, sensor_levels, τ_sum_all,
                                        F₀, vza, qp_μ, iμ₀, μ₀, pol_type)
    FT = eltype(F₀)
    μ₀ > zero(FT) || return direct_dwJ

    qp_μ_cpu = _to_cpu(qp_μ)
    solar_views = findall(i -> nearest_point(qp_μ_cpu, cosd(vza[i])) == iμ₀,
                          eachindex(vza))
    isempty(solar_views) && return direct_dwJ

    F₀_cpu = _to_cpu(F₀)
    nSpec = size(F₀_cpu, 2)
    size(τ_sum_all, 1) == nSpec || throw(DimensionMismatch(
        "cumulative optical depth has $(size(τ_sum_all, 1)) spectral points, " *
        "but F₀ has $nSpec"))

    # This is the same m=0 Stokes selector used by _precompute_vza_weights.
    selector = FT[one(FT), one(FT), zero(FT), zero(FT)][1:pol_type.n]
    F₀_selected = F₀_cpu .* reshape(selector, :, 1)
    weight_m0 = FT(0.5 / π)

    for ims in eachindex(sensor_levels)
        boundary = sensor_levels[ims]
        boundary == 0 && continue
        1 <= boundary < size(τ_sum_all, 2) ||
            throw(BoundsError(τ_sum_all, (:, boundary + 1)))
        # Transfer only the requested cumulative-τ column. Copying the full
        # nSpec×(Nz+1) profile here is unnecessarily expensive on a GPU.
        τ_above = vec(_to_cpu(@view τ_sum_all[:, boundary + 1]))
        attenuation = exp.(-τ_above ./ μ₀)
        direct = weight_m0 .* F₀_selected .* reshape(attenuation, 1, :)
        for i in solar_views
            @views direct_dwJ[ims][i, :, :] .= direct
        end
    end
    return direct_dwJ
end

"Multi-sensor post-processing: elastic (noRS)"
function postprocessing_vza_ms!(RS_type::noRS,
        sensor_levels,
        iμ₀, pol_type,
        composite_layer,
        vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI,
        uwJ, dwJ, uwieJ, dwieJ,
        I_static,
        arr_type)

    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    topR⁺⁻ = composite_layer.topR⁺⁻
    botR⁻⁺ = composite_layer.botR⁻⁺
    topJ₀⁺ = composite_layer.topJ₀⁺
    topJ₀⁻ = composite_layer.topJ₀⁻
    botJ₀⁺ = composite_layer.botJ₀⁺
    botJ₀⁻ = composite_layer.botJ₀⁻

    nμ = size(topJ₀⁺[1], 1)
    el_type = eltype(topJ₀⁺[1])
    # Interlayer solves must stay on the model architecture. Convert each
    # completed sensor field to the host once before VZA post-processing.
    tuwJ = [arr_type(zeros(el_type, (nμ, 1, nSpec))) for _ in sensor_levels]
    tdwJ = [arr_type(zeros(el_type, (nμ, 1, nSpec))) for _ in sensor_levels]

    for ims in eachindex(sensor_levels)
        if sensor_levels[ims] == 0
            tuwJ[ims] .= botJ₀⁻[ims]
            tdwJ[ims] .= botJ₀⁺[ims]
        else
            compute_interlayer_flux!(RS_type, I_static,
                topR⁺⁻[ims], botR⁻⁺[ims],
                topJ₀⁺[ims], botJ₀⁻[ims],
                tdwJ[ims], tuwJ[ims], arr_type)
        end
    end

    tuwJ_cpu = _to_cpu.(tuwJ)
    tdwJ_cpu = _to_cpu.(tdwJ)

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for ims in eachindex(sensor_levels)
            for s = 1:nSpec
                uwJ[ims][i,:,s] .+= w * tuwJ_cpu[ims][istart:iend, 1, s]
                dwJ[ims][i,:,s] .+= w * tdwJ_cpu[ims][istart:iend, 1, s]
            end
        end
    end
end

"Multi-sensor post-processing: Raman / inelastic"
function postprocessing_vza_ms!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus},
        sensor_levels,
        iμ₀, pol_type,
        composite_layer,
        vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI,
        uwJ, dwJ, uwieJ, dwieJ,
        I_static,
        arr_type)

    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    topR⁺⁻  = composite_layer.topR⁺⁻
    botR⁻⁺  = composite_layer.botR⁻⁺
    topJ₀⁺  = composite_layer.topJ₀⁺
    topJ₀⁻  = composite_layer.topJ₀⁻
    botJ₀⁺  = composite_layer.botJ₀⁺
    botJ₀⁻  = composite_layer.botJ₀⁻
    botieR⁻⁺  = composite_layer.botieR⁻⁺
    topieR⁺⁻  = composite_layer.topieR⁺⁻
    topieJ₀⁺  = composite_layer.topieJ₀⁺
    botieJ₀⁺  = composite_layer.botieJ₀⁺
    botieJ₀⁻  = composite_layer.botieJ₀⁻

    nμ = size(topJ₀⁺[1], 1)
    el_type = eltype(topJ₀⁺[1])
    n_raman = size(topieJ₀⁺[1], 4)

    tuwJ   = [arr_type(zeros(el_type, (nμ, 1, nSpec)))          for _ in sensor_levels]
    tdwJ   = [arr_type(zeros(el_type, (nμ, 1, nSpec)))          for _ in sensor_levels]
    tuwieJ = [arr_type(zeros(el_type, (nμ, 1, nSpec, n_raman))) for _ in sensor_levels]
    tdwieJ = [arr_type(zeros(el_type, (nμ, 1, nSpec, n_raman))) for _ in sensor_levels]

    for ims in eachindex(sensor_levels)
        if sensor_levels[ims] == 0
            tuwJ[ims]   .= botJ₀⁻[ims]
            tdwJ[ims]   .= botJ₀⁺[ims]
            tuwieJ[ims] .= botieJ₀⁻[ims]
            tdwieJ[ims] .= botieJ₀⁺[ims]
        else
            compute_interlayer_flux!(RS_type, I_static,
                topR⁺⁻[ims], botR⁻⁺[ims],
                topJ₀⁺[ims], botJ₀⁻[ims],
                tdwJ[ims], tuwJ[ims],
                topieR⁺⁻[ims], botieR⁻⁺[ims],
                topieJ₀⁺[ims], botieJ₀⁻[ims],
                tdwieJ[ims], tuwieJ[ims], arr_type)
        end
    end

    # Pre-reduce the inelastic sources over the Raman dimension once per sensor
    # (VZA-independent), mirroring the single-sensor postprocessing_vza! — instead
    # of re-summing n_raman terms inside the vza×spec loop (the old hot loop).
    tuwJ_cpu   = _to_cpu.(tuwJ)
    tdwJ_cpu   = _to_cpu.(tdwJ)
    tuwieJr = [_to_cpu(dropdims(sum(t, dims=4), dims=4)) for t in tuwieJ]
    tdwieJr = [_to_cpu(dropdims(sum(t, dims=4), dims=4)) for t in tdwieJ]

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for ims in eachindex(sensor_levels)
            for s = 1:nSpec
                uwJ[ims][i,:,s]   .+= w * tuwJ_cpu[ims][istart:iend, 1, s]
                dwJ[ims][i,:,s]   .+= w * tdwJ_cpu[ims][istart:iend, 1, s]
                uwieJ[ims][i,:,s] .+= w * tuwieJr[ims][istart:iend, 1, s]
                dwieJ[ims][i,:,s] .+= w * tdwieJr[ims][istart:iend, 1, s]
            end
        end
    end
end

"Multi-sensor canopy post-processing: elastic (noRS)"
function postprocessing_vza_ms_canopy!(RS_type::noRS,
        sensor_levels,
        quad_points,
        iμ₀, pol_type,
        composite_layer,
        solJ₀,
        vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI,
        uwJ, dwJ, uwieJ, dwieJ,
        hdr_J₀⁻, bhr_uw, bhr_dw,
        I_static,
        arr_type)

    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    topR⁺⁻ = composite_layer.topR⁺⁻
    botR⁻⁺ = composite_layer.botR⁻⁺
    topJ₀⁺ = composite_layer.topJ₀⁺
    topJ₀⁻ = composite_layer.topJ₀⁻
    botJ₀⁺ = composite_layer.botJ₀⁺
    botJ₀⁻ = composite_layer.botJ₀⁻

    el_type = eltype(topJ₀⁺[1])
    nμ = size(topJ₀⁺[1], 1)
    tuwJ = [arr_type(zeros(el_type, (nμ, 1, nSpec))) for _ in sensor_levels]
    tdwJ = [arr_type(zeros(el_type, (nμ, 1, nSpec))) for _ in sensor_levels]

    for ims in eachindex(sensor_levels)
        if sensor_levels[ims] == 0
            tuwJ[ims] .= botJ₀⁻[ims]
            tdwJ[ims] .= botJ₀⁺[ims]
        else
            compute_interlayer_flux!(RS_type, I_static,
                topR⁺⁻[ims], botR⁻⁺[ims],
                topJ₀⁺[ims], botJ₀⁻[ims],
                tdwJ[ims], tuwJ[ims], arr_type)
        end
    end

    interaction_hdrf_canopy!(SFI, tdwJ[2], tuwJ[2], solJ₀,
                            m, pol_type, quad_points,
                            hdr_J₀⁻, bhr_uw, bhr_dw)

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for ims in eachindex(sensor_levels)
            _tuwJ = _to_cpu(tuwJ[ims])
            _tdwJ = _to_cpu(tdwJ[ims])
            for s = 1:nSpec
                uwJ[ims][i,:,s] .+= w * _tuwJ[istart:iend, 1, s]
                dwJ[ims][i,:,s] .+= w * _tdwJ[istart:iend, 1, s]
            end
        end
    end
end
