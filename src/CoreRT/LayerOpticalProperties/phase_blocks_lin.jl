# Forward/tangent phase blocks, including spectral-node interpolation.
function _compute_phase_blocks_lin(model, greek, lin_greek, m, arr_type)
    q = model.quad_points
    pol = CoreRT.polarization_type(model)
    if q.external_solar
        μ = collect(q.qp_μ)
        Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺ = Scattering.compute_Z_moments(
            pol, μ, greek, lin_greek, m; arr_type=arr_type)
        Z₀⁺, Z₀⁻, Ż₀⁺, Ż₀⁻ = Scattering.compute_Z_source_moments(
            pol, μ, q.μ₀, greek, lin_greek, m; arr_type=arr_type)
        return Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺, Z₀⁺, Z₀⁻, Ż₀⁺, Ż₀⁻
    end
    Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺ = Scattering.compute_Z_moments(
        pol, collect(q.qp_μ), greek, lin_greek, m; arr_type=arr_type)
    return Z⁺⁺, Z⁻⁺, Ż⁺⁺, Ż⁻⁺, nothing, nothing, nothing, nothing
end

"Evaluate forward/tangent aerosol phase blocks at 2–3 retained nodes and interpolate."
function _compute_aerosol_phase_blocks_lin(model, optics, lin_optics, ν_spec,
                                            m, arr_type)
    # Both routes receive derivatives of the TRUNCATED Greek coefficients:
    # lin_model_from_parameters truncates each value/tangent pair before
    # storing it (or its spectral nodes). truncate_phase_lin differentiates
    # the fit and all six families' 1/(1-fᵗ) normalization. Consequently these
    # Z tangents contain the phase part of dfᵗ, not only the raw Mie dZ.
    # At fixed interpolation knots, angular evaluation and node interpolation
    # are linear in those coefficients, so the same maps act on the tangents.
    optics.phase_ν === nothing && return _compute_phase_blocks_lin(
        model, optics.greek_coefs, lin_optics.lin_greek_coefs, m, arr_type)
    length(optics.phase_ν) == length(optics.phase_greek) ==
        length(lin_optics.phase_lin_greek) || throw(DimensionMismatch(
        "aerosol phase-node values and tangents must have matching lengths"))
    blocks = [_compute_phase_blocks_lin(model, g, lg, m, Array) for
              (g, lg) in zip(optics.phase_greek, lin_optics.phase_lin_greek)]
    interp(k) = interpolate_phase_blocks(ν_spec, optics.phase_ν,
                                          [b[k] for b in blocks], arr_type)
    first_four = ntuple(interp, 4)
    if blocks[1][5] === nothing
        return first_four..., nothing, nothing, nothing, nothing
    end
    return first_four..., interp(5), interp(6), interp(7), interp(8)
end

