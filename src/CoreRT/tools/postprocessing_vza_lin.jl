#=
Azimuthal weighting for linearized (Jacobian) RT matrices.
Uses the same precomputed VZA weights as the forward postprocessing.
=#

"""
    postprocessing_vza!(RS_type::noRS, iμ₀, pol_type, composite_layer, composite_layer_lin,
                        vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI, R_SFI, T_SFI, Ṙ_SFI, Ṫ_SFI)

Azimuthally-weight linearized RT matrices (Jacobians).

Accumulates forward source terms `J₀⁺`, `J₀⁻` into `R_SFI`, `T_SFI`, and their
per-parameter derivatives `J̇₀⁺`, `J̇₀⁻` into `Ṙ_SFI`, `Ṫ_SFI`. Uses the same
`_precompute_vza_weights` as the forward postprocessing.

Architecture-aware: no-copy on CPU, minimal GPU→CPU transfer on GPU.
"""
function postprocessing_vza!(RS_type::noRS,
                    iμ₀, pol_type,
                    composite_layer,
                    composite_layer_lin,
                    vza, qp_μ, m, vaz, μ₀,
                    weight, nSpec,
                    SFI,
                    R_SFI, T_SFI,
                    Ṙ_SFI, Ṫ_SFI)

    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    J₀⁺ = _to_cpu(composite_layer.J₀⁺)
    J₀⁻ = _to_cpu(composite_layer.J₀⁻)
    J̇₀⁺ = _to_cpu(composite_layer_lin.J̇₀⁺)
    J̇₀⁻ = _to_cpu(composite_layer_lin.J̇₀⁻)

    Nparams = size(J̇₀⁻, 4)

    @inbounds for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        # Vectorise over the spectral dimension with one matrix product per VZA
        # (and, for the Jacobian, per parameter) instead of a scalar `for s` loop
        # that slice-allocated per spectral point × per parameter. `w` is the
        # azimuthal Stokes weight (Diagonal for n>1, scalar for n=1), so this is the
        # SAME `w * column` Stokes mix on every spectral column — identical result.
        # NOTE: matrix multiply `*`, not elementwise `.*`.
        @views R_SFI[i, :, :] .+= w * J₀⁻[istart:iend, 1, :]
        @views T_SFI[i, :, :] .+= w * J₀⁺[istart:iend, 1, :]
        for iparam = 1:Nparams
            @views Ṙ_SFI[i, :, :, iparam] .+= w * J̇₀⁻[istart:iend, 1, :, iparam]
            @views Ṫ_SFI[i, :, :, iparam] .+= w * J̇₀⁺[istart:iend, 1, :, iparam]
        end
    end
end

"""
    postprocessing_vza_ms_lin!(RS_type, top_layers, top_layers_lin,
                               bottom_layers, bottom_layers_lin, ...)

Solve and azimuthally reconstruct the diffuse radiance and analytic Jacobian
at each strict-interior observer interface.  Each top/bottom pair represents
the independently accumulated atmospheric subcolumns on either side of that
interface.  The coupled interlayer field is solved on the model architecture;
only the completed source vectors are transferred to the host outputs.
"""
function postprocessing_vza_ms_lin!(RS_type::noRS,
                    top_layers, top_layers_lin,
                    bottom_layers, bottom_layers_lin,
                    pol_type, vza, qp_μ, m, vaz, weight, nSpec,
                    uwJ, dwJ, uwJ_lin, dwJ_lin,
                    I_static, arr_type)

    isempty(top_layers) && return nothing
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    nμ = size(first(top_layers).J₀⁺, 1)
    Nparams = size(first(top_layers_lin).J̇₀⁺, 4)
    FT = eltype(first(top_layers).J₀⁺)

    tuwJ = [arr_type(zeros(FT, nμ, 1, nSpec)) for _ in top_layers]
    tdwJ = [arr_type(zeros(FT, nμ, 1, nSpec)) for _ in top_layers]
    tuwJ_lin = [arr_type(zeros(FT, nμ, 1, nSpec, Nparams)) for _ in top_layers]
    tdwJ_lin = [arr_type(zeros(FT, nμ, 1, nSpec, Nparams)) for _ in top_layers]

    for ims in eachindex(top_layers)
        top = top_layers[ims]
        top_lin = top_layers_lin[ims]
        bottom = bottom_layers[ims]
        bottom_lin = bottom_layers_lin[ims]
        compute_interlayer_flux!(RS_type, I_static,
            top.R⁺⁻, bottom.R⁻⁺,
            top.J₀⁺, bottom.J₀⁻,
            tdwJ[ims], tuwJ[ims],
            top_lin.Ṙ⁺⁻, bottom_lin.Ṙ⁻⁺,
            top_lin.J̇₀⁺, bottom_lin.J̇₀⁻,
            tdwJ_lin[ims], tuwJ_lin[ims], arr_type)
    end

    tuwJ_cpu = _to_cpu.(tuwJ)
    tdwJ_cpu = _to_cpu.(tdwJ)
    tuwJ_lin_cpu = _to_cpu.(tuwJ_lin)
    tdwJ_lin_cpu = _to_cpu.(tdwJ_lin)

    @inbounds for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for ims in eachindex(top_layers)
            @views uwJ[ims][i, :, :] .+= w * tuwJ_cpu[ims][istart:iend, 1, :]
            @views dwJ[ims][i, :, :] .+= w * tdwJ_cpu[ims][istart:iend, 1, :]
            for iparam in 1:Nparams
                @views uwJ_lin[ims][i, :, :, iparam] .+=
                    w * tuwJ_lin_cpu[ims][istart:iend, 1, :, iparam]
                @views dwJ_lin[ims][i, :, :, iparam] .+=
                    w * tdwJ_lin_cpu[ims][istart:iend, 1, :, iparam]
            end
        end
    end
    return nothing
end
