#=
Linearized Cox-Munk ocean surface.

Provides Jacobians of the surface reflectance matrix and SFI source terms
with respect to wind speed, following the pattern in lambertian_surface_lin.jl.

The derivative chain is:
  dR_surf/dU = dR_surf/dσ² · dσ²/dU
where dσ²/dU = 0.00256 for the per-axis slope variance σ² = ½(0.003 + 0.00512U)
(Cox & Munk 1954; see `wind_to_sigma2`), plus whitecap derivatives.

Uses the analytical `reflectance_and_deriv` function from coxmunk_surface.jl,
which shares all wind-speed-independent work (geometry, Fresnel reflection, Stokes basis changes)
between the forward and derivative evaluations.
=#

"""
    create_surface_layer!(RS_type::noRS, surf::CoxMunkSurface, added_layer, added_layer_lin,
                          iparam, SFI, m, pol_type, quad_points, τ_sum, τ̇_sum, F₀, architecture)

Compute forward + linearized surface reflection for a Cox-Munk surface.

The surface parameter is wind speed.  The derivative `∂R_surf/∂U` is computed
analytically via the derivative chain through the slope PDF, shadow factor,
and whitecap fraction — all sharing the wind-speed-independent Fresnel and
geometry computations with the forward evaluation.

A spectrally varying `n_water` is handled as in the forward builder: the
reflectance and its wind derivative are evaluated at the knots of
`_n_water_knots` and interpolated with the same weights, so forward
and Jacobian use the same per-point operator.
"""
function create_surface_layer!(RS_type::noRS,
                                surf::CoxMunkSurface{FT},
                                added_layer::AddedLayer,
                                added_layer_lin::AddedLayerLin,
                                iparam::Int,
                                SFI,
                                m::Int,
                                pol_type,
                                quad_points,
                                τ_sum, τ̇_sum, F₀,
                                architecture) where {FT}

    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points
    (; n) = pol_type
    arr_type = array_type(architecture)

    nparams = size(τ̇_sum, 2)
    nspec = length(τ_sum)
    Nquad = size(added_layer.r⁻⁺, 1) ÷ n
    tmp = arr_type(ones(FT, n * Nquad))
    T_surf = Diagonal(tmp)
    i₀ = iμ₀Nstart:(iμ₀Nstart + n - 1)
    n_water = _spectral_n_water(surf, nspec)

    n_water isa Complex || return _create_spectral_surface_layer_lin!(
        surf, n_water, added_layer, added_layer_lin, iparam, SFI, m, pol_type,
        quad_points, τ_sum, τ̇_sum, F₀, arr_type, T_surf)

    # ── Forward reflectance + analytical derivative in one pass ──
    μ_arr = collect(qp_μ)
    if m == 0
        ρ_raw, dρ_raw = reflectance_and_deriv(surf, pol_type, μ_arr, m; n_water)
        ρ          = FT(2) * ρ_raw
        dρ_dU_raw  = FT(2) * dρ_raw
    else
        ρ, dρ_dU_raw = reflectance_and_deriv(surf, pol_type, μ_arr, m; n_water)
    end
    R_surf = arr_type(ρ)
    Ṙ_surf = arr_type(dρ_dU_raw)

    # ── SFI source terms + derivatives ──
    if SFI
        F₀_NquadN = arr_type(zeros(FT, length(qp_μN), nspec))
        Ḟ₀_NquadN = arr_type(zeros(FT, length(qp_μN), nspec, nparams + 1))

        tmpF = F₀ .* arr_type(exp.(-τ_sum / μ₀))'
        F₀_NquadN[i₀, :] .= tmpF
        Ḟ₀_NquadN[i₀, :, 1:nparams] .= -reshape(tmpF, n, nspec, 1) .*
            reshape(τ̇_sum, 1, nspec, nparams) / μ₀

        added_layer.j₀⁺[:, :, :] .= zero(FT)
        added_layer.j₀⁻[:, 1, :] .= μ₀ * (R_surf * F₀_NquadN)

        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        for ii in 1:nspec
            for ctr in 1:nparams
                added_layer_lin.ap_J̇₀⁻[:, 1, ii, ctr] .= μ₀ * R_surf * Ḟ₀_NquadN[:, ii, ctr]
            end
            # Derivative of source w.r.t. wind speed (analytical)
            added_layer_lin.ap_J̇₀⁻[:, 1, ii, iparam] .= μ₀ * Ṙ_surf * F₀_NquadN[:, ii]
        end
    end

    # ── Apply quadrature weights ──
    R_surf_weighted = R_surf * Diagonal(qp_μN .* wt_μN)
    Ṙ_surf_weighted = Ṙ_surf * Diagonal(qp_μN .* wt_μN)

    # ── Store in AddedLayer ──
    added_layer.r⁻⁺ .= R_surf_weighted
    added_layer.r⁺⁻ .= zero(FT)
    added_layer.t⁺⁺ .= T_surf
    added_layer.t⁻⁻ .= zero(FT)

    # ── Store derivatives in AddedLayerLin ──
    added_layer_lin.ap_ṙ⁻⁺[:, :, :, iparam] .= Ṙ_surf_weighted
    added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
    added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
    added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)

    synchronize_if_gpu()
end

# Spectrally varying refractive index: per-point reflectance and wind
# derivative, interpolated between the operator knots with common weights.
function _create_spectral_surface_layer_lin!(surf::CoxMunkSurface{FT}, n_water,
                                             added_layer, added_layer_lin, iparam,
                                             SFI, m, pol_type, quad_points,
                                             τ_sum, τ̇_sum, F₀, arr_type,
                                             T_surf) where {FT}
    (; qp_μ, qp_μN, wt_μN, iμ₀Nstart, μ₀) = quad_points
    (; n) = pol_type
    nparams = size(τ̇_sum, 2)
    nspec = length(τ_sum)
    NquadN = size(added_layer.r⁻⁺, 1)
    i₀ = iμ₀Nstart:(iμ₀Nstart + n - 1)
    μ_arr = collect(qp_μ)
    fac = m == 0 ? FT(2) : one(FT)
    Rk, dRk, lo, hi, t = _spectral_knot_operators(n_water, arr_type) do nw
        ρ, dρ = reflectance_and_deriv(surf, pol_type, μ_arr, m; n_water = nw)
        (fac * ρ, fac * dρ)
    end
    # Unweighted per-point operators first (for the direct-beam sources).
    R_surf = added_layer.r⁻⁺
    _interpolate_knots!(R_surf, Rk, lo, hi, t)
    Ṙ_surf = similar(R_surf)
    _interpolate_knots!(Ṙ_surf, dRk, lo, hi, t)

    if SFI
        F₀_NquadN = arr_type(zeros(FT, NquadN, nspec))
        Ḟ₀_NquadN = arr_type(zeros(FT, NquadN, nspec, nparams + 1))
        tmpF = F₀ .* arr_type(exp.(-τ_sum / μ₀))'
        F₀_NquadN[i₀, :] .= tmpF
        Ḟ₀_NquadN[i₀, :, 1:nparams] .= -reshape(tmpF, n, nspec, 1) .*
            reshape(τ̇_sum, 1, nspec, nparams) / μ₀
        beam = reshape(F₀_NquadN, NquadN, 1, nspec)

        added_layer.j₀⁺ .= zero(FT)
        added_layer.j₀⁻ .= μ₀ .* batched_mul(R_surf, beam)
        added_layer_lin.ap_J̇₀⁺ .= zero(FT)
        for ctr in 1:nparams
            added_layer_lin.ap_J̇₀⁻[:, :, :, ctr] .= μ₀ .* batched_mul(
                R_surf, reshape(Ḟ₀_NquadN[:, :, ctr], NquadN, 1, nspec))
        end
        # Derivative of source w.r.t. wind speed (analytical)
        added_layer_lin.ap_J̇₀⁻[:, :, :, iparam] .= μ₀ .* batched_mul(Ṙ_surf, beam)
    end

    # ── Apply quadrature weights (column scaling, in place) ──
    weights = reshape(qp_μN .* wt_μN, 1, :)
    R_surf .*= weights
    Ṙ_surf .*= weights
    added_layer.r⁺⁻ .= zero(FT)
    added_layer.t⁺⁺ .= T_surf
    added_layer.t⁻⁻ .= zero(FT)

    added_layer_lin.ap_ṙ⁻⁺[:, :, :, iparam] .= Ṙ_surf
    added_layer_lin.ap_ṙ⁺⁻ .= zero(FT)
    added_layer_lin.ap_ṫ⁺⁺ .= zero(FT)
    added_layer_lin.ap_ṫ⁻⁻ .= zero(FT)

    synchronize_if_gpu()
end

# ──────────────────────────────────────────────────────────────────────────────
# Finite-difference version (kept for validation / comparison)
# ──────────────────────────────────────────────────────────────────────────────

"""
    reflectance_fd_deriv(surf, pol_type, μ, m; n_water, ΔU)

Compute Fourier reflectance and its finite-difference derivative w.r.t. wind speed.
Useful for validating the analytical derivative from `reflectance_and_deriv`.

Returns `(R, dR_dU_fd)` — both `[Nμ·n_stokes, Nμ·n_stokes]`.
"""
function reflectance_fd_deriv(surf::CoxMunkSurface{FT}, pol_type,
                               μ::AbstractArray{FT}, m::Int;
                               n_water::Complex{FT} = _single_n_water(surf),
                               ΔU::FT = max(FT(1e-4), FT(1e-4) * surf.wind_speed)) where FT
    surf_plus  = CoxMunkSurface{FT}(wind_speed = surf.wind_speed + ΔU,
                                      n_water = surf.n_water,
                                      whitecap_albedo = surf.whitecap_albedo,
                                      include_whitecaps = surf.include_whitecaps,
                                      shadowing = surf.shadowing)
    surf_minus = CoxMunkSurface{FT}(wind_speed = max(FT(0), surf.wind_speed - ΔU),
                                      n_water = surf.n_water,
                                      whitecap_albedo = surf.whitecap_albedo,
                                      include_whitecaps = surf.include_whitecaps,
                                      shadowing = surf.shadowing)

    R       = reflectance(surf,       pol_type, μ, m; n_water=n_water)
    R_plus  = reflectance(surf_plus,  pol_type, μ, m; n_water=n_water)
    R_minus = reflectance(surf_minus, pol_type, μ, m; n_water=n_water)

    actual_ΔU = (surf.wind_speed + ΔU) - max(FT(0), surf.wind_speed - ΔU)
    dR_fd = (R_plus - R_minus) / actual_ΔU

    return R, dR_fd
end
