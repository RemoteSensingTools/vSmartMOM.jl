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
