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

    J₀⁺ = Array(composite_layer.J₀⁺)
    J₀⁻ = Array(composite_layer.J₀⁻)
    J̇₀⁺ = Array(composite_layer_lin.J̇₀⁺)
    J̇₀⁻ = Array(composite_layer_lin.J̇₀⁻)

    Nparams = size(J̇₀⁻, 1)

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for s = 1:nSpec
            R_SFI[i,:,s] .+= w * J₀⁻[istart:iend, 1, s]
            T_SFI[i,:,s] .+= w * J₀⁺[istart:iend, 1, s]
            for iparam = 1:Nparams
                Ṙ_SFI[iparam,i,:,s] .+= w * J̇₀⁻[iparam, istart:iend, 1, s]
                Ṫ_SFI[iparam,i,:,s] .+= w * J̇₀⁺[iparam, istart:iend, 1, s]
            end
        end
    end
end
