#=
Azimuthal weighting of RT matrices after all Fourier-moment kernel calculations.
Accumulates cos(m·φ) / sin(m·φ) weighted source terms into the output arrays
R_SFI, T_SFI (and their inelastic counterparts for Raman).
=#

_to_cpu(x::Array) = x
_to_cpu(x) = Array(x)

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
    postprocessing_vza_toa!(composite_layer, vza, qp_μ, m, vaz,
                             pol_type, weight, R_TOA)

Accumulate only the upwelling SFI source at TOA. This is the lean output path
for disk-integrated exoplanet calculations: no BOA directional field and no
HDR/BHR diagnostics are constructed or postprocessed.
"""
function postprocessing_vza_toa!(composite_layer, vza, qp_μ, m, vaz,
                                  pol_type, weight, R_TOA)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)
    J₀⁻ = _to_cpu(composite_layer.J₀⁻)
    @inbounds for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        @views R_TOA[i, :, :] .+= w * J₀⁻[istart:iend, 1, :]
    end
    for cslot in values(composite_layer.J₀_by_src)
        J⁻_src = _to_cpu(cslot.J₀⁻)
        @inbounds for i in eachindex(vza)
            istart, iend, w = vza_info[i]
            @views R_TOA[i, :, :] .+= w * J⁻_src[istart:iend, 1, :]
        end
    end
    return nothing
end

"""TOA-only elastic and Raman source accumulation for an `RRS` composite."""
function postprocessing_vza_toa!(composite_layer::CompositeLayerRS,
                                  vza, qp_μ, m, vaz, pol_type, weight,
                                  R_TOA, ieR_TOA)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)
    J₀⁻ = _to_cpu(composite_layer.J₀⁻)
    ieJ₀⁻ = _to_cpu(dropdims(sum(composite_layer.ieJ₀⁻, dims=4), dims=4))
    @inbounds for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        @views R_TOA[i,:,:] .+= w * J₀⁻[istart:iend,1,:]
        @views ieR_TOA[i,:,:] .+= w * ieJ₀⁻[istart:iend,1,:]
    end
    nothing
end

"""
    postprocessing_vza!(RS_type::noRS, iμ₀, pol_type, composite_layer, vza, qp_μ, m, vaz, μ₀, weight, nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)

Azimuthally-weight RT matrices for elastic (no Raman) scattering.

Accumulates cos(m·φ)-weighted reflectance/transmittance from the composite layer
into `R`, `T` (collimated) or `R_SFI`, `T_SFI` (source function integration).
Uses quadrature indices for solar direction `iμ₀` and view directions `vza`, `vaz`.

Architecture-aware: no-copy on CPU, minimal GPU→CPU transfer on GPU.
"""
function postprocessing_vza!(RS_type::noRS, iμ₀, pol_type,
        composite_layer, vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)

    _, istart0, iend0 = get_indices(iμ₀, pol_type)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    if SFI
        # Legacy solar slot accumulation
        J₀⁺ = _to_cpu(composite_layer.J₀⁺)
        J₀⁻ = _to_cpu(composite_layer.J₀⁻)
        @inbounds for i in eachindex(vza)
            istart, iend, w = vza_info[i]
            # Vectorise over the spectral dimension: one matrix product per VZA
            # instead of a scalar `for s` loop that slice-allocated two tiny arrays
            # per spectral point (millions of allocations for nSpec≈40k). `w` is the
            # azimuthal Stokes weight (a Diagonal for n>1, scalar for n=1), so this
            # is the SAME `w * column` Stokes mix applied to every spectral column at
            # once — identical result. NOTE: matrix multiply `*`, not elementwise `.*`.
            @views R_SFI[i, :, :] .+= w * J₀⁻[istart:iend, 1, :]
            @views T_SFI[i, :, :] .+= w * J₀⁺[istart:iend, 1, :]
        end
        # v0.7 Phase A.2a — per-source J₀ slots. RT reconstruction is linear
        # in sources so each slot's contribution adds into the same R_SFI/T_SFI
        # output. The same Fourier `weight` and per-VZA cosine apply (thermal
        # only contributed at m=0 anyway — see `thermal_emission.jl`).
        for cslot in values(composite_layer.J₀_by_src)
            J⁺_src = _to_cpu(cslot.J₀⁺)
            J⁻_src = _to_cpu(cslot.J₀⁻)
            @inbounds for i in eachindex(vza)
                istart, iend, w = vza_info[i]
                @views R_SFI[i, :, :] .+= w * J⁻_src[istart:iend, 1, :]
                @views T_SFI[i, :, :] .+= w * J⁺_src[istart:iend, 1, :]
            end
        end
    else
        R⁻⁺ = _to_cpu(composite_layer.R⁻⁺)
        T⁺⁺ = _to_cpu(composite_layer.T⁺⁺)
        @inbounds for i in eachindex(vza)
            istart, iend, w = vza_info[i]
            for s = 1:nSpec
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
    hdr_J₀⁻ = _to_cpu(hdr_J₀⁻)

    @inbounds for i in eachindex(vza)
        istart, iend, w = vza_info[i]
        @views hdr[i, :, :] .+= w * hdr_J₀⁻[istart:iend, 1, :]
    end
end

"""
    postprocessing_vza!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus}, ...)

Azimuthally-weight RT matrices for Raman/inelastic scattering.

Same as elastic `postprocessing_vza!` but also accumulates the inelastic source
terms `ieJ₀⁺`, `ieJ₀⁻` (summed over all Raman shifts) into `ieR_SFI`, `ieT_SFI`.
"""
function postprocessing_vza!(RS_type::Union{RRS, VS_0to1_plus, VS_1to0_plus},
        iμ₀, pol_type, composite_layer,
        vza, qp_μ, m, vaz, μ₀, weight,
        nSpec, SFI, R, R_SFI, T, T_SFI, ieR_SFI, ieT_SFI)

    _, istart0, iend0 = get_indices(iμ₀, pol_type)
    vza_info = _precompute_vza_weights(vza, vaz, qp_μ, pol_type, m, weight)

    if SFI
        J₀⁺   = _to_cpu(composite_layer.J₀⁺)
        J₀⁻   = _to_cpu(composite_layer.J₀⁻)
        # The inelastic output sums over ALL Raman shifts:
        #   ieR_SFI[i,:,s] = Σₜ w·ieJ₀⁻[…,s,t] = w·(Σₜ ieJ₀⁻[…,s,t]).
        # That Raman-dim sum is VZA-independent, so reduce it ONCE on the device
        # (a single memory-bound pass) — collapsing ieJ₀± to a 3-D source with the
        # same shape as the elastic J₀± — instead of re-summing n_raman terms
        # inside the vza×spec loop (the old hot loop; ~400× slower at nSpec~10⁴).
        # Only the reduced (NquadN,1,nSpec) array crosses PCIe, not the full 4-D.
        ieJ₀⁺ = _to_cpu(dropdims(sum(composite_layer.ieJ₀⁺, dims=4), dims=4))
        ieJ₀⁻ = _to_cpu(dropdims(sum(composite_layer.ieJ₀⁻, dims=4), dims=4))

        @inbounds for i in eachindex(vza)
            istart, iend, w = vza_info[i]
            for s = 1:nSpec
                R_SFI[i,:,s]   .+= w * J₀⁻[istart:iend, 1, s]
                T_SFI[i,:,s]   .+= w * J₀⁺[istart:iend, 1, s]
                ieR_SFI[i,:,s] .+= w * ieJ₀⁻[istart:iend, 1, s]
                ieT_SFI[i,:,s] .+= w * ieJ₀⁺[istart:iend, 1, s]
            end
        end
    else
        R⁻⁺ = _to_cpu(composite_layer.R⁻⁺)
        T⁺⁺ = _to_cpu(composite_layer.T⁺⁺)
        @inbounds for i in eachindex(vza)
            istart, iend, w = vza_info[i]
            for s = 1:nSpec
                R[i,:,s] .+= w * (R⁻⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀
                T[i,:,s] .+= w * (T⁺⁺[istart:iend, istart0:iend0, s] / μ₀) * pol_type.I₀
            end
        end
    end
end
