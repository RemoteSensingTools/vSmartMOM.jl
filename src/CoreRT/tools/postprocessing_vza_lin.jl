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
        for s = 1:nSpec
            R_SFI[i,:,s] .+= w * J₀⁻[istart:iend, 1, s]
            T_SFI[i,:,s] .+= w * J₀⁺[istart:iend, 1, s]
            for iparam = 1:Nparams
                Ṙ_SFI[i,:,s,iparam] .+= w * J̇₀⁻[istart:iend, 1, s, iparam]
                Ṫ_SFI[i,:,s,iparam] .+= w * J̇₀⁺[istart:iend, 1, s, iparam]
            end
        end
    end
end
