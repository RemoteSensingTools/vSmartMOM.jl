#=
Azimuthal weighting for multi-sensor RT postprocessing.
Handles interlayer flux coupling before the VZA loop.
=#

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
    el_type = typeof(topJ₀⁺[1][1,1,1])
    tuwJ = [zeros(el_type, (nμ, 1, nSpec)) for _ in sensor_levels]
    tdwJ = [zeros(el_type, (nμ, 1, nSpec)) for _ in sensor_levels]

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

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for ims in eachindex(sensor_levels)
            for s = 1:nSpec
                uwJ[ims][i,:,s] .+= w * tuwJ[ims][istart:iend, 1, s]
                dwJ[ims][i,:,s] .+= w * tdwJ[ims][istart:iend, 1, s]
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
    el_type = typeof(topJ₀⁺[1][1,1,1])
    n_raman = size(topieJ₀⁺[1], 4)

    tuwJ   = [zeros(el_type, (nμ, 1, nSpec))          for _ in sensor_levels]
    tdwJ   = [zeros(el_type, (nμ, 1, nSpec))          for _ in sensor_levels]
    tuwieJ = [zeros(el_type, (nμ, 1, nSpec, n_raman)) for _ in sensor_levels]
    tdwieJ = [zeros(el_type, (nμ, 1, nSpec, n_raman)) for _ in sensor_levels]

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

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for ims in eachindex(sensor_levels)
            for s = 1:nSpec
                uwJ[ims][i,:,s] .+= w * tuwJ[ims][istart:iend, 1, s]
                dwJ[ims][i,:,s] .+= w * tdwJ[ims][istart:iend, 1, s]
                for t = 1:n_raman
                    uwieJ[ims][i,:,s] .+= w * tuwieJ[ims][istart:iend, 1, s, t]
                    dwieJ[ims][i,:,s] .+= w * tdwieJ[ims][istart:iend, 1, s, t]
                end
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
