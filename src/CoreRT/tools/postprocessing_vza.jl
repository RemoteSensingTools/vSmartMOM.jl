#=
Azimuthal weighting of RT matrices after all Fourier-moment kernel calculations.
Accumulates cos(m·φ) / sin(m·φ) weighted source terms into the output arrays
R_SFI, T_SFI (and their inelastic counterparts for Raman).
=#

"""
    _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

Precompute per-VZA index ranges and azimuthal weight matrices for postprocessing.

For each view zenith angle (VZA), finds the nearest quadrature point, Stokes indices
`(istart, iend)`, and azimuthal weight ``w = \\text{weight} \\cdot \\cos(m \\phi)`` (or
``\\sin(m \\phi)`` for Stokes 3,4). For scalar RT, `w` is a scalar; for polarized RT,
`w` is a Diagonal matrix over Stokes components.

# Returns
- Vector of `(istart, iend, w)` tuples, one per VZA.
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

"""
    postprocessing_vza!(RS_type::noRS, iμ₀, pol_type, composite_layer, vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)

Azimuthally-weight RT matrices for elastic (no Raman) scattering.

Accumulates cos(m·φ)-weighted reflectance/transmittance from the composite layer
into `R`, `T` (collimated) or `R_SFI`, `T_SFI` (source function integration).
Uses quadrature indices for solar direction `iμ₀` and view directions `vza`, `vaz`.
"""
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

"""
    postprocessing_vza_hdrf!(RS_type, iμ₀, pol_type, hdr_J₀⁻, vza, qp_μ, m, vaz, μ₀, weight, nSpec, hdr)

RAMI benchmark: azimuthally-weight hemispherical-directional reflectance (HDR) matrices.

Accumulates weighted `hdr_J₀⁻` (downwelling source at surface) into `hdr` for each VZA.
"""
function postprocessing_vza_hdrf!(RS_type, iμ₀, pol_type,
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

"""
    postprocessing_vza!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, ...)

Azimuthally-weight RT matrices for Raman/inelastic scattering.

Same as elastic `postprocessing_vza!` but also accumulates inelastic source terms
`ieJ₀⁺`, `ieJ₀⁻` into `ieR_SFI`, `ieT_SFI` for each Raman shift.
"""
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
