#=
Azimuthal weighting of RT matrices after all Fourier-moment kernel calculations.
Accumulates cos(m·φ) / sin(m·φ) weighted source terms into the output arrays
R_SFI, T_SFI (and their inelastic counterparts for Raman).
=#

"""
    _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

Precompute per-VZA index ranges and azimuthal weight matrices/scalars.
Returns a vector of `(istart, iend, w)` tuples where `w` is either a
scalar (n_stokes==1) or a Diagonal matrix.
"""
function _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)
    n = pol_type.n
    map(eachindex(vza)) do i
        iμ = nearest_point(qp_μ, cosd(vza[i]))
        _, istart, iend = get_indices(iμ, pol_type)
        cos_m_phi = cosd(m * vaz[i])
        sin_m_phi = sind(m * vaz[i])
        w = if n == 1
            weight * cos_m_phi
        else
            weight * Diagonal([cos_m_phi, cos_m_phi, sin_m_phi, sin_m_phi][1:n])
        end
        (istart, iend, w)
    end
end

"Perform post-processing to azimuthally-weight RT matrices (elastic, no Raman)"
function postprocessing_vza!(RS_type::noRS, iμ₀, pol_type,
        composite_layer, vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)

    _, istart0, iend0 = get_indices(iμ₀, pol_type)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    R⁻⁺ = collect(composite_layer.R⁻⁺)
    T⁺⁺ = collect(composite_layer.T⁺⁺)
    J₀⁺ = collect(composite_layer.J₀⁺)
    J₀⁻ = collect(composite_layer.J₀⁻)

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for s = 1:nSpec
            if SFI
                R_SFI[i,:,s] .+= w * J₀⁻[istart:iend, 1, s]
                T_SFI[i,:,s] .+= w * J₀⁺[istart:iend, 1, s]
            else
                R[i,:,s] .+= w * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀
                T[i,:,s] .+= w * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀
            end
        end
    end
end

"RAMI: Perform post-processing to azimuthally-weight hdr matrices"
function postprocessing_vza_hdrf!(RS_type::noRS, iμ₀, pol_type,
        hdr_J₀⁻, vza, qp_μ, m, vaz, μ₀, weight, nSpec, hdr)

    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)
    hdr_J₀⁻ = collect(hdr_J₀⁻)

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for s = 1:nSpec
            hdr[i,:,s] .+= w * hdr_J₀⁻[istart:iend, 1, s]
        end
    end
end

"Perform post-processing to azimuthally-weight RT matrices (Raman / inelastic)"
function postprocessing_vza!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus},
        iμ₀, pol_type, composite_layer,
        vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)

    _, istart0, iend0 = get_indices(iμ₀, pol_type)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    R⁻⁺ = collect(composite_layer.R⁻⁺)
    T⁺⁺ = collect(composite_layer.T⁺⁺)
    J₀⁺ = collect(composite_layer.J₀⁺)
    J₀⁻ = collect(composite_layer.J₀⁻)
    ieJ₀⁺ = collect(composite_layer.ieJ₀⁺)
    ieJ₀⁻ = collect(composite_layer.ieJ₀⁻)

    n_raman = size(ieJ₀⁺, 4)

    for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        for s = 1:nSpec
            if SFI
                R_SFI[i,:,s]  .+= w * J₀⁻[istart:iend, 1, s]
                T_SFI[i,:,s]  .+= w * J₀⁺[istart:iend, 1, s]
                for t = 1:n_raman
                    ieR_SFI[i,:,s] .+= w * ieJ₀⁻[istart:iend, 1, s, t]
                    ieT_SFI[i,:,s] .+= w * ieJ₀⁺[istart:iend, 1, s, t]
                end
            else
                R[i,:,s] .+= w * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀
                T[i,:,s] .+= w * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀
            end
        end
    end
end
