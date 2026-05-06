#=
Cox-Munk ocean surface BRDF with full polarization support.

Implements the Cox-Munk (1954) wind-roughened ocean surface with:
- Fresnel reflection from tilted wave facets (full Mueller matrix)
- Isotropic Gaussian slope distribution parameterized by wind speed
- Optional Lambertian whitecap contribution (Monahan & O'Muircheartaigh 1980)
- Smith (1967) shadow masking function
- TMS (Truncated Multiple Scattering) single-scattering correction for the specular hotspot

References:
- Cox & Munk (1954), JOSA 44(11), doi:10.1364/JOSA.44.000838
- Mishchenko & Travis (1997), JGR 102(D14), doi:10.1029/97JD01084
- Zhai et al. (2010), Opt. Express 18(9), doi:10.1364/OE.18.009613
=#

# ──────────────────────────────────────────────────────────────────────
# Helper functions
# ──────────────────────────────────────────────────────────────────────

"""Isotropic slope variance σ² as a function of 10-m wind speed U [m/s] (Cox & Munk 1954)."""
wind_to_sigma2(U::FT) where FT = FT(0.003) + FT(0.00512) * U

"""Isotropic Gaussian slope PDF for slope components (zx, zy)."""
function cox_munk_pdf(zx::FT, zy::FT, σ²::FT) where FT
    return exp(-(zx^2 + zy^2) / (2σ²)) / (2 * FT(π) * σ²)
end

"""Whitecap fractional coverage (Monahan & O'Muircheartaigh 1980)."""
function whitecap_fraction(U::FT) where FT
    U <= zero(FT) && return zero(FT)
    return FT(2.95e-6) * U^FT(3.52)
end

"""
    smith_shadowing(μ, σ²)

Smith (1967) monostatic shadowing function Λ(μ) for a Gaussian slope distribution.
Returns the shadowing factor S = 1/(1 + Λ), where Λ is the shadowing integral.
"""
function smith_Λ(μ::FT, σ²::FT) where FT
    μ <= zero(FT) && return FT(1e10)  # fully shadowed at grazing
    σ = sqrt(σ²)
    cot_θ = μ / sqrt(max(FT(1e-30), one(FT) - μ^2))
    ν = cot_θ / (sqrt(FT(2)) * σ)
    # Λ = (exp(-ν²) / (√(2π) ν) - erfc(ν)) / 2
    # Use the approximation valid for all ν (Tsang et al. 1985):
    Λ = (exp(-ν^2) / (sqrt(FT(2π)) * ν) - erfc(ν)) / 2
    return max(zero(FT), Λ)
end

"""Bistatic shadowing factor for incident (μᵢ) and reflected (μᵣ) directions."""
function shadow_factor(μᵢ::FT, μᵣ::FT, σ²::FT) where FT
    return one(FT) / (one(FT) + smith_Λ(μᵢ, σ²) + smith_Λ(μᵣ, σ²))
end

# ──────────────────────────────────────────────────────────────────────
# Analytical derivatives w.r.t. σ² and wind speed
# ──────────────────────────────────────────────────────────────────────

"""Derivative of the slope PDF w.r.t. σ²: ``∂P/∂σ² = P (Z² − 2σ²)/(2σ⁴)``."""
function cox_munk_pdf_dσ²(zx::FT, zy::FT, σ²::FT) where FT
    Z² = zx^2 + zy^2
    P = cox_munk_pdf(zx, zy, σ²)
    return P * (Z² - FT(2) * σ²) / (FT(2) * σ²^2)
end

"""Derivative of Smith Λ(μ) w.r.t. σ²."""
function smith_Λ_dσ²(μ::FT, σ²::FT) where FT
    μ <= zero(FT) && return zero(FT)
    σ = sqrt(σ²)
    cot_θ = μ / sqrt(max(FT(1e-30), one(FT) - μ^2))
    ν = cot_θ / (sqrt(FT(2)) * σ)

    # Check if Λ itself is zero (derivative of max(0,·) is 0 there)
    Λ_raw = (exp(-ν^2) / (sqrt(FT(2π)) * ν) - erfc(ν)) / 2
    Λ_raw <= zero(FT) && return zero(FT)

    # dΛ/dν via product rule on exp(-ν²)/(√(2π) ν) and chain rule on erfc(ν):
    #   d/dν[exp(-ν²)/(√(2π)ν)] = exp(-ν²)(-2ν² − 1)/(√(2π)ν²)
    #   d/dν[erfc(ν)]            = −(2/√π) exp(-ν²)
    exp_neg_ν² = exp(-ν^2)
    dΛ_dν = (exp_neg_ν² * (-FT(2) * ν^2 - one(FT)) / (sqrt(FT(2π)) * ν^2) +
             FT(2) / sqrt(FT(π)) * exp_neg_ν²) / 2

    # dν/dσ² = −ν / (2σ²)
    dν_dσ² = -ν / (FT(2) * σ²)

    return dΛ_dν * dν_dσ²
end

"""Derivative of the bistatic shadow factor w.r.t. σ²: ``∂S/∂σ² = −S² (∂Λᵢ/∂σ² + ∂Λᵣ/∂σ²)``."""
function shadow_factor_dσ²(μᵢ::FT, μᵣ::FT, σ²::FT) where FT
    Λ_i = smith_Λ(μᵢ, σ²)
    Λ_r = smith_Λ(μᵣ, σ²)
    S = one(FT) / (one(FT) + Λ_i + Λ_r)
    dΛ_i = smith_Λ_dσ²(μᵢ, σ²)
    dΛ_r = smith_Λ_dσ²(μᵣ, σ²)
    return -S^2 * (dΛ_i + dΛ_r)
end

"""Derivative of whitecap fraction w.r.t. wind speed U: ``∂f_{wc}/∂U = 2.95×10⁻⁶ · 3.52 · U^{2.52}``."""
function whitecap_fraction_deriv(U::FT) where FT
    U <= zero(FT) && return zero(FT)
    return FT(2.95e-6) * FT(3.52) * U^FT(2.52)
end

# ──────────────────────────────────────────────────────────────────────
# Core BRDF Mueller matrix
# ──────────────────────────────────────────────────────────────────────

"""
    coxmunk_geometry(μᵢ, μᵣ, dϕ)

Compute the facet-reflection geometry for given incident (μᵢ, azimuth=0)
and reflected (μᵣ, azimuth=dϕ) directions.

Returns a NamedTuple with:
- `cos_β`: cosine of facet tilt angle
- `cos_θ_local`: cosine of the local incidence angle on the facet
- `zx, zy`: slope components
- `α₁, α₂`: rotation angles (scattering plane ↔ facet incidence plane)
"""
function coxmunk_geometry(μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    sin_θᵢ = sqrt(max(FT(0), one(FT) - μᵢ^2))
    sin_θᵣ = sqrt(max(FT(0), one(FT) - μᵣ^2))
    cos_dϕ = cos(dϕ)
    sin_dϕ = sin(dϕ)

    # Facet normal = bisector of (-incident) and reflected directions
    # Incident direction: (sin_θᵢ, 0, -μᵢ)  [downward]
    # Reflected direction: (sin_θᵣ cos(dϕ), sin_θᵣ sin(dϕ), μᵣ)  [upward]
    # Half-vector (unnormalized facet normal):
    nx = -sin_θᵢ + sin_θᵣ * cos_dϕ
    ny = sin_θᵣ * sin_dϕ
    nz = μᵢ + μᵣ
    norm_n = sqrt(nx^2 + ny^2 + nz^2)

    # Guard against degenerate geometry
    if norm_n < FT(1e-15)
        return (; cos_β=one(FT), cos_θ_local=one(FT),
                  zx=zero(FT), zy=zero(FT), α₁=zero(FT), α₂=zero(FT))
    end

    nx /= norm_n
    ny /= norm_n
    nz /= norm_n

    # Tilt angle
    cos_β = max(FT(1e-10), nz)

    # Local incidence angle on facet: cos_θ_local = |k_i · n|
    cos_θ_local = clamp((μᵢ + μᵣ) / (FT(2) * cos_β), FT(0), one(FT))

    # Slope components
    zx = -nx / cos_β
    zy = -ny / cos_β

    # ── Rotation angles (Mishchenko & Travis 1997, Eqs. 4-5) ──
    # α₁: rotation from scattering plane to the plane of incidence on the facet
    # α₂: rotation from the plane of reflection on the facet back to scattering plane

    # Scattering angle cosine
    cos_Θ = -μᵢ * μᵣ + sin_θᵢ * sin_θᵣ * cos_dϕ
    sin_Θ² = max(FT(0), one(FT) - cos_Θ^2)
    sin_Θ  = sqrt(sin_Θ²)

    if sin_Θ < FT(1e-12)
        # Forward or backward scattering: rotation angles are zero
        α₁ = zero(FT)
        α₂ = zero(FT)
    else
        # sin(α₁) and cos(α₁) from spherical geometry:
        # cos(α₁) = (μᵢ cos_Θ + μᵣ) / (sin_θᵢ · sin_Θ)  (not quite right for reflection)
        # We use the facet-based formulas from Zhai et al. (2010), Appendix A:

        # The local incidence plane is defined by (k_i, n).
        # We need the angle between the scattering plane (k_i, k_r) and this incidence plane.
        # Both planes contain k_i, so the angle between them is the dihedral angle.

        # Cross products to get plane normals:
        # Scattering plane normal: k_i × k_r  (defines the reference for Stokes)
        # Incidence plane normal:  k_i × n

        # k_i = (sin_θᵢ, 0, -μᵢ)  [in our convention: propagation direction]
        # k_r = (-sin_θᵣ cos_dϕ, -sin_θᵣ sin_dϕ, μᵣ)  [outgoing]
        # n   = (nx, ny, nz)

        # Scattering plane normal (k_i × k_r):
        sp_x = FT(0) * μᵣ - (-μᵢ) * (-sin_θᵣ * sin_dϕ)
        sp_y = (-μᵢ) * (-sin_θᵣ * cos_dϕ) - sin_θᵢ * μᵣ
        sp_z = sin_θᵢ * (-sin_θᵣ * sin_dϕ) - FT(0) * (-sin_θᵣ * cos_dϕ)
        # Simplify:
        sp_x = -μᵢ * sin_θᵣ * sin_dϕ
        sp_y = μᵢ * sin_θᵣ * cos_dϕ - sin_θᵢ * μᵣ
        sp_z = -sin_θᵢ * sin_θᵣ * sin_dϕ

        # Incidence plane normal (k_i × n):
        ip_x = FT(0) * nz - (-μᵢ) * ny
        ip_y = (-μᵢ) * nx - sin_θᵢ * nz
        ip_z = sin_θᵢ * ny - FT(0) * nx
        # Simplify:
        ip_x = μᵢ * ny
        ip_y = -μᵢ * nx - sin_θᵢ * nz
        ip_z = sin_θᵢ * ny

        # cos(α₁) = (sp · ip) / (|sp| |ip|)
        dot_si = sp_x * ip_x + sp_y * ip_y + sp_z * ip_z
        mag_sp = sqrt(sp_x^2 + sp_y^2 + sp_z^2)
        mag_ip = sqrt(ip_x^2 + ip_y^2 + ip_z^2)

        if mag_sp < FT(1e-15) || mag_ip < FT(1e-15)
            α₁ = zero(FT)
        else
            cos_α₁ = clamp(dot_si / (mag_sp * mag_ip), -one(FT), one(FT))
            # Sign from triple product (k_i · (sp × ip))
            cross_x = sp_y * ip_z - sp_z * ip_y
            cross_y = sp_z * ip_x - sp_x * ip_z
            cross_z = sp_x * ip_y - sp_y * ip_x
            sign_α₁ = sin_θᵢ * cross_x + FT(0) * cross_y + (-μᵢ) * cross_z
            α₁ = sign_α₁ >= 0 ? acos(cos_α₁) : -acos(cos_α₁)
        end

        # For α₂: angle between the reflection plane on the facet and the scattering plane
        # Reflection plane normal (k_r × n):
        rp_x = (-sin_θᵣ * sin_dϕ) * nz - μᵣ * ny
        rp_y = μᵣ * nx - (-sin_θᵣ * cos_dϕ) * nz
        rp_z = (-sin_θᵣ * cos_dϕ) * ny - (-sin_θᵣ * sin_dϕ) * nx

        dot_sr = sp_x * rp_x + sp_y * rp_y + sp_z * rp_z
        mag_rp = sqrt(rp_x^2 + rp_y^2 + rp_z^2)

        if mag_sp < FT(1e-15) || mag_rp < FT(1e-15)
            α₂ = zero(FT)
        else
            cos_α₂ = clamp(dot_sr / (mag_sp * mag_rp), -one(FT), one(FT))
            # Sign from triple product (k_r · (sp × rp))
            cross2_x = sp_y * rp_z - sp_z * rp_y
            cross2_y = sp_z * rp_x - sp_x * rp_z
            cross2_z = sp_x * rp_y - sp_y * rp_x
            sign_α₂ = (-sin_θᵣ * cos_dϕ) * cross2_x +
                       (-sin_θᵣ * sin_dϕ) * cross2_y +
                       μᵣ * cross2_z
            α₂ = sign_α₂ >= 0 ? acos(cos_α₂) : -acos(cos_α₂)
        end
    end

    return (; cos_β, cos_θ_local, zx, zy, α₁, α₂)
end

"""
    coxmunk_brdf_mueller(surf, n_stokes, μᵢ, μᵣ, dϕ; n_water)

Evaluate the Cox-Munk BRDF Mueller matrix at a single geometry.

Returns an `SMatrix{n_stokes, n_stokes}` with the BRDF value
(units: sr⁻¹), including whitecap contribution if enabled.
"""
function coxmunk_brdf_mueller(surf::CoxMunkSurface{FT}, n_stokes::Int,
                               μᵢ::FT, μᵣ::FT, dϕ::FT;
                               n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    σ² = wind_to_sigma2(surf.wind_speed)

    # Geometry
    geom = coxmunk_geometry(μᵢ, μᵣ, dϕ)
    (; cos_β, cos_θ_local, zx, zy, α₁, α₂) = geom

    # Slope PDF
    P = cox_munk_pdf(zx, zy, σ²)

    # Fresnel reflection
    r_s, r_p = fresnel_coefficients(n_water, cos_θ_local)
    M_F = fresnel_mueller(r_s, r_p, n_stokes)

    # Rotation matrices
    L₁ = stokes_rotation_matrix(-α₁, n_stokes)  # scattering plane → incidence plane
    L₂ = stokes_rotation_matrix(α₂, n_stokes)   # reflection plane → scattering plane

    # Full Mueller matrix for the facet contribution
    M_facet = L₂ * M_F * L₁

    # Cox-Munk BRDF weighting
    prefactor = P / (FT(4) * μᵢ * μᵣ * cos_β^4)

    # Shadow masking
    if surf.shadowing
        S = shadow_factor(μᵢ, μᵣ, σ²)
        prefactor *= S
    end

    ρ_glint = prefactor * M_facet

    # Whitecap contribution (unpolarized Lambertian)
    if surf.include_whitecaps
        f_wc = whitecap_fraction(surf.wind_speed)
        ρ_wc = _whitecap_mueller(surf.whitecap_albedo, n_stokes)
        return (one(FT) - f_wc) * ρ_glint + f_wc * ρ_wc
    else
        return ρ_glint
    end
end

"""
    coxmunk_brdf_mueller_and_deriv(surf, n_stokes, μᵢ, μᵣ, dϕ; n_water)

Evaluate the Cox-Munk BRDF Mueller matrix **and** its analytical derivative
w.r.t. wind speed U at a single geometry, sharing all U-independent work
(geometry, Fresnel, rotations).

Returns `(M, dM_dU)` — both `SMatrix{n_stokes, n_stokes}`.
"""
function coxmunk_brdf_mueller_and_deriv(surf::CoxMunkSurface{FT}, n_stokes::Int,
                                         μᵢ::FT, μᵣ::FT, dϕ::FT;
                                         n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    U  = surf.wind_speed
    σ² = wind_to_sigma2(U)
    dσ²_dU = FT(0.00512)

    # ── Geometry, Fresnel, rotations (all independent of U) ──
    geom = coxmunk_geometry(μᵢ, μᵣ, dϕ)
    (; cos_β, cos_θ_local, zx, zy, α₁, α₂) = geom

    r_s, r_p = fresnel_coefficients(n_water, cos_θ_local)
    M_F  = fresnel_mueller(r_s, r_p, n_stokes)
    L₁   = stokes_rotation_matrix(-α₁, n_stokes)
    L₂   = stokes_rotation_matrix(α₂, n_stokes)
    M_facet = L₂ * M_F * L₁

    geom_weight = one(FT) / (FT(4) * μᵢ * μᵣ * cos_β^4)

    # ── Forward + derivative of slope PDF ──
    P      = cox_munk_pdf(zx, zy, σ²)
    dP_dσ² = cox_munk_pdf_dσ²(zx, zy, σ²)

    # ── Forward + derivative of shadow factor ──
    if surf.shadowing
        S      = shadow_factor(μᵢ, μᵣ, σ²)
        dS_dσ² = shadow_factor_dσ²(μᵢ, μᵣ, σ²)
        prefactor      = P * S * geom_weight
        dprefactor_dσ² = (dP_dσ² * S + P * dS_dσ²) * geom_weight
    else
        prefactor      = P * geom_weight
        dprefactor_dσ² = dP_dσ² * geom_weight
    end

    ρ_glint      = prefactor * M_facet
    dρ_glint_dU  = (dprefactor_dσ² * dσ²_dU) * M_facet

    # ── Whitecap ──
    if surf.include_whitecaps
        f_wc     = whitecap_fraction(U)
        df_wc_dU = whitecap_fraction_deriv(U)
        ρ_wc     = _whitecap_mueller(surf.whitecap_albedo, n_stokes)

        M    = (one(FT) - f_wc) * ρ_glint  + f_wc * ρ_wc
        dM   = (one(FT) - f_wc) * dρ_glint_dU + df_wc_dU * (ρ_wc - ρ_glint)
        return M, dM
    else
        return ρ_glint, dρ_glint_dU
    end
end

"""
    reflectance_and_deriv(surf, pol_type, μ, m; n_water)

Fourier moment `m` of the Cox-Munk BRDF reflectance matrix **and** its
analytical derivative w.r.t. wind speed U, computed in a single pass over
the azimuthal quadrature.

Returns `(R, dR_dU)` — both `[Nμ·n_stokes, Nμ·n_stokes]` dense matrices.
"""
function reflectance_and_deriv(surf::CoxMunkSurface{FT}, pol_type,
                                μ::AbstractArray{FT}, m::Int;
                                n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    n  = pol_type.n
    Nμ = length(μ)
    nn = Nμ * n
    Rsurf  = zeros(FT, nn, nn)
    dRsurf = zeros(FT, nn, nn)

    nQuad_ϕ = 100
    ϕ, w = CanopyOptics.gauleg(nQuad_ϕ, FT(0), FT(π))

    for iϕ in eachindex(ϕ)
        dϕ = ϕ[iϕ]
        wϕ = w[iϕ]
        for j in 1:Nμ, i in 1:Nμ
            M, dM = coxmunk_brdf_mueller_and_deriv(surf, n, μ[i], μ[j], dϕ;
                                                     n_water=n_water)
            @inbounds for sj in 1:n, si in 1:n
                az  = _azimuthal_kernel(si, sj, m, dϕ)
                idx = (i-1)*n + si
                jdx = (j-1)*n + sj
                Rsurf[idx, jdx]  += wϕ * M[si, sj]  * az
                dRsurf[idx, jdx] += wϕ * dM[si, sj] * az
            end
        end
    end

    ff = m == 0 ? FT(1.0) : FT(2.0)
    return ff * Rsurf / FT(π), ff * dRsurf / FT(π)
end

"""Whitecap Mueller matrix: unpolarized Lambertian (only M[1,1] = albedo/π)."""
function _whitecap_mueller(albedo::FT, n_stokes::Int) where FT
    val = albedo / FT(π)
    if n_stokes == 1
        return SMatrix{1,1,FT}(val)
    elseif n_stokes == 2
        z = zero(FT)
        return SMatrix{2,2,FT}(val, z, z, z)
    elseif n_stokes == 3
        z = zero(FT)
        return SMatrix{3,3,FT}(val, z, z, z, z, z, z, z, z)
    else
        z = zero(FT)
        return SMatrix{4,4,FT}(val, z, z, z, z, z, z, z, z, z, z, z, z, z, z, z)
    end
end

"""Get refractive index from surface spec, falling back to built-in table."""
function _get_n_water(surf::CoxMunkSurface{FT}, λ_nm::FT) where FT
    if surf.n_water === nothing
        n = water_refractive_index(λ_nm)
        return Complex{FT}(real(n), imag(n))
    elseif surf.n_water isa Complex
        return surf.n_water::Complex{FT}
    else
        # Vector case: caller must select the right index
        return surf.n_water[1]::Complex{FT}
    end
end

# ──────────────────────────────────────────────────────────────────────
# Fourier decomposition: reflectance(surf, pol_type, μ, m)
# ──────────────────────────────────────────────────────────────────────

"""
    _azimuthal_kernel(si, sj, m, dϕ)

Azimuthal Fourier kernel for Stokes element (si, sj) at Fourier moment m.
Matches the postprocessing reconstruction weight `Diagonal([cos,cos,sin,sin])`.
"""
@inline function _azimuthal_kernel(si::Int, sj::Int, m::Int, dϕ::FT) where FT
    row_is_iq = (si <= 2)
    col_is_iq = (sj <= 2)
    if row_is_iq == col_is_iq
        return cos(m * dϕ)
    else
        return sin(m * dϕ)
    end
end

"""
    reflectance(surf::CoxMunkSurface, pol_type, μ, m)

Fourier moment `m` of the Cox-Munk BRDF reflectance matrix for quadrature
directions `μ`.

Returns `[Nμ·n_stokes, Nμ·n_stokes]` matrix including all Stokes coupling
(I-Q, U-V off-diagonal blocks).  Uses Gauss-Legendre quadrature over
azimuth [0, π] with 100 points.
"""
function reflectance(surf::CoxMunkSurface{FT}, pol_type, μ::AbstractArray{FT}, m::Int;
                     n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    n = pol_type.n
    Nμ = length(μ)
    nn = Nμ * n
    Rsurf = zeros(FT, nn, nn)

    nQuad_ϕ = 100
    ϕ, w = CanopyOptics.gauleg(nQuad_ϕ, FT(0), FT(π))

    for iϕ in eachindex(ϕ)
        dϕ = ϕ[iϕ]
        wϕ = w[iϕ]
        for j in 1:Nμ       # incident direction
            for i in 1:Nμ   # reflected direction
                M = coxmunk_brdf_mueller(surf, n, μ[i], μ[j], dϕ; n_water=n_water)
                # Place each Stokes sub-element with the correct azimuthal kernel
                @inbounds for sj in 1:n, si in 1:n
                    az = _azimuthal_kernel(si, sj, m, dϕ)
                    Rsurf[(i-1)*n + si, (j-1)*n + sj] += wϕ * M[si, sj] * az
                end
            end
        end
    end

    ff = m == 0 ? FT(1.0) : FT(2.0)
    return ff * Rsurf / FT(π)
end

# ──────────────────────────────────────────────────────────────────────
# Single-scattering (TMS) correction for the specular hotspot
# ──────────────────────────────────────────────────────────────────────

"""
    apply_ss_correction!(R_SFI, surf, pol_type, vza, vaz, μ₀, τ_total, max_m, nSpec;
                         n_water)

Truncated Multiple Scattering (TMS) correction for the specular sun glint peak.

After the Fourier loop, adds the difference between the exact single-scattering
surface contribution and the truncated Fourier reconstruction at each viewing geometry.

Modifies `R_SFI` in place.
"""
function apply_ss_correction!(R_SFI::AbstractArray{FT,3},
                               surf::CoxMunkSurface{FT},
                               pol_type, vza, vaz, μ₀::FT,
                               τ_total::AbstractVector{FT},
                               max_m::Int, nSpec::Int;
                               n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    n = pol_type.n

    for iv in eachindex(vza)
        μ_v = FT(cosd(vza[iv]))
        dϕ  = FT(deg2rad(vaz[iv]))

        # ── Exact single-scattering surface BRDF ──
        M_exact = coxmunk_brdf_mueller(surf, n, μ_v, μ₀, dϕ; n_water=n_water)

        # ── Fourier-reconstructed BRDF at this geometry ──
        M_fourier = zeros(FT, n, n)
        for m in 0:(max_m - 1)
            weight_m = m == 0 ? FT(0.5) : FT(1.0)
            for si in 1:n, sj in 1:n
                az = _azimuthal_kernel(si, sj, m, dϕ)
                # Re-evaluate the Fourier coefficient at the specific (μ_v, μ₀) pair
                # via the same quadrature used in reflectance()
                M_fourier[si, sj] += weight_m * az * _fourier_coeff_element(
                    surf, n, si, sj, μ_v, μ₀, m; n_water=n_water)
            end
        end

        # Apply the azimuthal weight (matching postprocessing_vza.jl)
        cos_m0_ϕ = one(FT)  # reconstruction at the actual azimuth is already done above
        # The Fourier coefficients were already weighted by cos/sin(m*dϕ) in the loop

        # Correction for each spectral point
        for s in 1:nSpec
            atten = μ₀ * exp(-τ_total[s] / μ₀)
            for si in 1:n
                correction = atten * (M_exact[si, 1] - M_fourier[si, 1])
                R_SFI[iv, si, s] += correction
            end
        end
    end
end

"""
Compute one element of the Fourier coefficient of the BRDF at a specific (μ_v, μ₀) pair.
Uses the same Gauss-Legendre quadrature as `reflectance`.
"""
function _fourier_coeff_element(surf::CoxMunkSurface{FT}, n_stokes::Int,
                                 si::Int, sj::Int, μᵢ::FT, μⱼ::FT, m::Int;
                                 n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    nQuad = 100
    ϕ, w = CanopyOptics.gauleg(nQuad, FT(0), FT(π))
    result = zero(FT)
    ff = m == 0 ? FT(1.0) : FT(2.0)
    for iϕ in eachindex(ϕ)
        M = coxmunk_brdf_mueller(surf, n_stokes, μᵢ, μⱼ, ϕ[iϕ]; n_water=n_water)
        az = _azimuthal_kernel(si, sj, m, ϕ[iϕ])
        result += w[iϕ] * M[si, sj] * az
    end
    return ff * result / FT(π)
end
