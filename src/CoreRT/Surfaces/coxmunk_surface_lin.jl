#=
Linearized Cox-Munk ocean surface.

Provides Jacobians of the surface reflectance matrix and SFI source terms
with respect to wind speed, following the pattern in lambertian_surface_lin.jl.

The derivative chain is:
  dR_surf/dU = dR_surf/dσ² · dσ²/dU
where dσ²/dU = 0.00512 (Cox & Munk 1954), plus whitecap derivatives.

Uses the analytical `reflectance_and_deriv` function from coxmunk_surface.jl,
which shares all wind-speed-independent work (geometry, Fresnel, rotations)
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

    # ── Forward reflectance + analytical derivative in one pass ──
    μ_arr = collect(qp_μ)
    if m == 0
        ρ_raw, dρ_raw = reflectance_and_deriv(surf, pol_type, μ_arr, m)
        ρ          = FT(2) * ρ_raw
        dρ_dU_raw  = FT(2) * dρ_raw
    else
        ρ, dρ_dU_raw = reflectance_and_deriv(surf, pol_type, μ_arr, m)
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
                               n_water::Complex{FT} = _get_n_water(surf, FT(550)),
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
