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

"""
    wind_to_sigma2(U)

Per-axis slope variance σ² of the isotropic Gaussian facet-slope distribution
for 10-m wind speed `U` in m/s. Cox & Munk (1954) give the *total* mean-square
slope ⟨zx² + zy²⟩ = 0.003 + 0.00512U; each of the two slope components carries
half of it, σ² = ½(0.003 + 0.00512U), which is the variance used by
`cox_munk_pdf` (and by VLIDORT's GISS Cox–Munk kernel).
"""
wind_to_sigma2(U::FT) where FT = (FT(0.003) + FT(0.00512) * U) / 2

"""dσ²/dU of `wind_to_sigma2`: ½ · 0.00512 = 0.00256 (m/s)⁻¹."""
wind_to_sigma2_dU(::Type{FT}) where FT = FT(0.00512) / 2

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
    smith_Λ(μ, σ²)

Smith (1967) shadowing integral Λ(μ), μ = cos θ, for a Gaussian distribution
of the slope component in the plane of the ray with variance σ² (the
per-axis variance, see `wind_to_sigma2`):

```math
Λ = \\int_c^∞ \\left(\\frac{s}{c} - 1\\right) f_σ(s)\\,ds
  = \\frac{1}{2}\\left[\\frac{e^{-ν^2}}{\\sqrt{π}\\,ν} - \\operatorname{erfc}(ν)\\right],
\\qquad c = \\cot θ,\\; ν = \\frac{c}{\\sqrt{2}\\,σ}.
```

`shadow_factor` combines Λ for the incident and reflected directions.
"""
function smith_Λ(μ::FT, σ²::FT) where FT
    μ <= zero(FT) && return FT(1e10)  # fully shadowed at grazing
    σ = sqrt(σ²)
    cot_θ = μ / sqrt(max(FT(1e-30), one(FT) - μ^2))
    ν = cot_θ / (sqrt(FT(2)) * σ)
    Λ = (exp(-ν^2) / (sqrt(FT(π)) * ν) - erfc(ν)) / 2
    # Λ ≥ 0 analytically; the guard only removes roundoff for large ν.
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
    Λ_raw = (exp(-ν^2) / (sqrt(FT(π)) * ν) - erfc(ν)) / 2
    Λ_raw <= zero(FT) && return zero(FT)

    # dΛ/dν via product rule on exp(-ν²)/(√π ν) and chain rule on erfc(ν):
    #   d/dν[exp(-ν²)/(√π ν)] = exp(-ν²)(-2ν² − 1)/(√π ν²)
    #   d/dν[erfc(ν)]         = −(2/√π) exp(-ν²)
    # The exp(-ν²)-proportional terms cancel, leaving
    dΛ_dν = -exp(-ν^2) / (2 * sqrt(FT(π)) * ν^2)

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
    coxmunk_geometry(μᵢ, μᵣ, Δϕ)

Facet geometry for specular reflection of the incident propagation direction
`kᵢ = (sinθᵢ, 0, −μᵢ)` (downward, azimuth 0) into `kᵣ = (sinθᵣ cosΔϕ,
sinθᵣ sinΔϕ, μᵣ)` (upward, azimuth Δϕ). This is the azimuth convention of the
atmospheric phase matrix (`Scattering.phase_matrix_first_column`): forward
reflection, and hence the sun glint, lies at Δϕ = 0.

The reflecting facet has normal `n̂ ∝ kᵣ − kᵢ`, tilt `cos β = n̂_z`, local
incidence `cos θ = kᵣ·n̂`, and slopes `zx = −n̂ₓ/n̂_z`, `zy = −n̂_y/n̂_z`.

Fresnel reflection is diagonal in the facet basis `(ŝ, p̂ = k × ŝ)`, where
`ŝ ∝ kᵢ × kᵣ` is normal to the plane of reflection. vSmartMOM Stokes vectors
refer to the meridian basis `(φ̂, θ̂)` of each beam, with
`Q = |E_φ|² − |E_θ|²`, `U = 2Re(E_φ E_θ*)`, `V = 2Im(E_φ E_θ*)`. The field
components transform as

    (E_s, E_p) = [cos χ  sin χ; sin χ  −cos χ] (E_φ, E_θ),   ŝ = cos χ φ̂ + sin χ θ̂,

an improper (reflecting) change of basis that is its own inverse. Its Stokes
matrix, `meridian_facet_stokes`, depends on `cos 2χ` and `sin 2χ`; `χ₁`
belongs to the incident and `χ₂` to the reflected beam.
"""
function coxmunk_geometry(μᵢ::FT, μᵣ::FT, Δϕ::FT) where FT
    sinθᵢ = sqrt(max(zero(FT), one(FT) - μᵢ^2))
    sinθᵣ = sqrt(max(zero(FT), one(FT) - μᵣ^2))
    cosΔϕ, sinΔϕ = cos(Δϕ), sin(Δϕ)

    # Propagation directions and their meridian bases (φ̂, θ̂)
    kᵢ  = SVector(sinθᵢ, zero(FT), -μᵢ)
    kᵣ  = SVector(sinθᵣ * cosΔϕ, sinθᵣ * sinΔϕ, μᵣ)
    φ̂ᵢ = SVector(zero(FT), one(FT), zero(FT))
    θ̂ᵢ = SVector(-μᵢ, zero(FT), -sinθᵢ)
    φ̂ᵣ = SVector(-sinΔϕ, cosΔϕ, zero(FT))
    θ̂ᵣ = SVector(μᵣ * cosΔϕ, μᵣ * sinΔϕ, -sinθᵣ)

    # Facet normal (half vector), tilt, local incidence, slopes. For unit
    # vectors |kᵣ − kᵢ| = 2 cos θ_local; using the norm avoids the cancellation
    # of dot(kᵣ, n̂) for nearly parallel grazing rays in Float32.
    h = kᵣ - kᵢ
    norm_h = norm(h)
    n̂ = h / norm_h
    cos_β = max(FT(1e-10), n̂[3])
    cos_θ_local = clamp(norm_h / 2, zero(FT), one(FT))
    zx = -n̂[1] / cos_β
    zy = -n̂[2] / cos_β

    # Normal to the plane of reflection; at exact backscatter (kᵣ = −kᵢ, i.e.
    # μᵣ = μᵢ and Δϕ = π: normal incidence on a facet tilted by θᵢ) every
    # ŝ ⟂ kᵢ is valid.
    s_vec = cross(kᵢ, kᵣ)
    ŝ = norm(s_vec) > eps(FT) ? normalize(s_vec) : φ̂ᵢ

    # cos χ = ŝ·φ̂ and sin χ = ŝ·θ̂ for each beam
    cos_χ₁, sin_χ₁ = dot(ŝ, φ̂ᵢ), dot(ŝ, θ̂ᵢ)
    cos_χ₂, sin_χ₂ = dot(ŝ, φ̂ᵣ), dot(ŝ, θ̂ᵣ)
    return (; cos_β, cos_θ_local, zx, zy,
              cos_2χ₁ = cos_χ₁^2 - sin_χ₁^2, sin_2χ₁ = 2 * cos_χ₁ * sin_χ₁,
              cos_2χ₂ = cos_χ₂^2 - sin_χ₂^2, sin_2χ₂ = 2 * cos_χ₂ * sin_χ₂)
end

"""
    meridian_facet_stokes(cos_2χ, sin_2χ)

Stokes matrix of the field-basis change
`(E_s, E_p) = [cos χ  sin χ; sin χ  −cos χ](E_φ, E_θ)` (see `coxmunk_geometry`):

    [1     0        0      0
     0   cos 2χ   sin 2χ   0
     0   sin 2χ  −cos 2χ   0
     0     0        0     −1]

It is symmetric and involutory, so the same matrix maps meridian → facet
Stokes vectors for the incident beam and facet → meridian for the
reflected beam.
"""
@inline function meridian_facet_stokes(cos_2χ::FT, sin_2χ::FT) where FT
    o, z = one(FT), zero(FT)
    # Column-major; the matrix is symmetric.
    return SMatrix{4,4,FT}(o, z,       z,       z,
                           z, cos_2χ,  sin_2χ,  z,
                           z, sin_2χ, -cos_2χ,  z,
                           z, z,       z,      -o)
end

"Leading `n_stokes × n_stokes` block of a 4×4 Stokes matrix."
@inline function _leading_stokes_block(M::SMatrix{4,4,FT}, n_stokes::Int) where FT
    n_stokes == 4 && return M
    n_stokes == 3 && return M[SOneTo(3), SOneTo(3)]
    n_stokes == 2 && return M[SOneTo(2), SOneTo(2)]
    return M[SOneTo(1), SOneTo(1)]
end

"""
    coxmunk_facet_mueller(geom, n_water, n_stokes)

Polarized Fresnel reflection of one facet in the meridian bases of the
incident and reflected beams: `L(χ₂) · M_F · L(χ₁)`, truncated to
`n_stokes` after the full 4×4 product (truncating the factors first would
drop the U-mediated part of the I–Q coupling).
"""
@inline function coxmunk_facet_mueller(geom, n_water::Complex{FT}, n_stokes::Int) where FT
    r_s, r_p = fresnel_coefficients(n_water, geom.cos_θ_local)
    M = meridian_facet_stokes(geom.cos_2χ₂, geom.sin_2χ₂) *
        fresnel_mueller(r_s, r_p, 4) *
        meridian_facet_stokes(geom.cos_2χ₁, geom.sin_2χ₁)
    return _leading_stokes_block(M, n_stokes)
end

"""
    coxmunk_brdf_mueller(surf, n_stokes, μᵢ, μᵣ, Δϕ; n_water)

Cox–Munk BRDF Mueller matrix (sr⁻¹) for light incident with cosine `μᵢ`
(azimuth 0) and reflected into cosine `μᵣ` at relative azimuth `Δϕ` (radians;
the glint is at Δϕ = 0, see `coxmunk_geometry`). Element `[i, j]` maps
incident Stokes component `j` to reflected component `i`, in vSmartMOM's
meridian-plane convention. Includes Smith shadowing and the whitecap
contribution when enabled.
"""
function coxmunk_brdf_mueller(surf::CoxMunkSurface{FT}, n_stokes::Int,
                               μᵢ::FT, μᵣ::FT, dϕ::FT;
                               n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    σ² = wind_to_sigma2(surf.wind_speed)

    # Facet geometry and its polarized Fresnel reflection
    geom = coxmunk_geometry(μᵢ, μᵣ, dϕ)
    M_facet = coxmunk_facet_mueller(geom, n_water, n_stokes)

    # Cox-Munk BRDF weighting
    P = cox_munk_pdf(geom.zx, geom.zy, σ²)
    # At extreme grazing P underflows to 0 while 4μᵢμᵣcos⁴β can underflow too
    # (Float32); the glint term is then exactly 0, not 0·Inf = NaN.
    prefactor = iszero(P) ? zero(FT) : P / (FT(4) * μᵢ * μᵣ * geom.cos_β^4)

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
(geometry, Fresnel reflection, Stokes basis changes).

Returns `(M, dM_dU)` — both `SMatrix{n_stokes, n_stokes}`.
"""
function coxmunk_brdf_mueller_and_deriv(surf::CoxMunkSurface{FT}, n_stokes::Int,
                                         μᵢ::FT, μᵣ::FT, dϕ::FT;
                                         n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    U  = surf.wind_speed
    σ² = wind_to_sigma2(U)
    dσ²_dU = wind_to_sigma2_dU(FT)

    # ── Geometry and polarized Fresnel reflection (independent of U) ──
    geom = coxmunk_geometry(μᵢ, μᵣ, dϕ)
    (; cos_β, zx, zy) = geom
    M_facet = coxmunk_facet_mueller(geom, n_water, n_stokes)

    # ── Forward + derivative of slope PDF ──
    P      = cox_munk_pdf(zx, zy, σ²)
    dP_dσ² = cox_munk_pdf_dσ²(zx, zy, σ²)
    # P and dP/dσ² vanish together at extreme grazing, where the weight
    # 1/(4μᵢμᵣcos⁴β) can overflow (Float32): the glint terms are then 0.
    geom_weight = iszero(P) ? zero(FT) : one(FT) / (FT(4) * μᵢ * μᵣ * cos_β^4)

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

The RT operator convention is `π * BRDF`, matching the generic surface
builder and its final `1/π` radiance reconstruction.
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
        for j in 1:Nμ, i in 1:Nμ       # incident j → reflected i
            M, dM = coxmunk_brdf_mueller_and_deriv(surf, n, μ[j], μ[i], dϕ;
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
    return ff * Rsurf, ff * dRsurf
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

Azimuthal Fourier kernel for Stokes element (si, sj) at Fourier moment m, in
the storage convention of the atmospheric phase-matrix moments
(`Scattering.compute_Z_moments`), so that surface and atmospheric operators
combine consistently in the adding method:

| block (row ← column) | kernel        |
|----------------------|---------------|
| IQ ← IQ, UV ← UV     | `cos(mϕ)`     |
| UV ← IQ              | `sin(mϕ)`     |
| IQ ← UV              | `−sin(mϕ)`    |

The same kernel synthesizes the real-space matrix from its moments, and the
UV ← IQ column matches the postprocessing weights for I, Q (cos) and U, V (sin).
"""
@inline function _azimuthal_kernel(si::Int, sj::Int, m::Int, dϕ::FT) where FT
    row_is_iq = (si <= 2)
    col_is_iq = (sj <= 2)
    if row_is_iq == col_is_iq
        return cos(m * dϕ)
    elseif row_is_iq
        return -sin(m * dϕ)
    else
        return sin(m * dϕ)
    end
end

"""
    reflectance(surf::CoxMunkSurface, pol_type, μ, m)

Fourier moment `m` of the Cox-Munk BRDF reflectance matrix for quadrature
directions `μ`.

Returns the Fourier moments of `π * BRDF` as an
`[Nμ·n_stokes, Nμ·n_stokes]` matrix, including all Stokes coupling
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
                M = coxmunk_brdf_mueller(surf, n, μ[j], μ[i], dϕ; n_water=n_water)
                # Place each Stokes sub-element with the correct azimuthal kernel
                @inbounds for sj in 1:n, si in 1:n
                    az = _azimuthal_kernel(si, sj, m, dϕ)
                    Rsurf[(i-1)*n + si, (j-1)*n + sj] += wϕ * M[si, sj] * az
                end
            end
        end
    end

    ff = m == 0 ? FT(1.0) : FT(2.0)
    return ff * Rsurf
end

# ──────────────────────────────────────────────────────────────────────
# Single-scattering (TMS) correction for the specular hotspot
# ──────────────────────────────────────────────────────────────────────

"""
    _coxmunk_glint_residual(surf, n, μ_v, μ₀, dϕ, m_max, n_water)

Exact minus truncated direct-solar BRDF column and its analytic wind tangent.
The surface builder doubles the zeroth moment and postprocessing halves it;
there is no additional half-weight in this physical BRDF reconstruction.
"""
function _coxmunk_glint_residual(surf::CoxMunkSurface{FT}, n, μ_v, μ₀,
                                  dϕ, m_max, n_water) where FT
    exact, dexact = coxmunk_brdf_mueller_and_deriv(
        surf, n, μ₀, μ_v, dϕ; n_water)
    truncated = zeros(FT, n)
    dtruncated = zeros(FT, n)
    ϕ, w = CanopyOptics.gauleg(100, zero(FT), FT(π))
    for i in eachindex(ϕ)
        M, dM = coxmunk_brdf_mueller_and_deriv(
            surf, n, μ₀, μ_v, ϕ[i]; n_water)
        for m in 0:m_max, si in 1:n
            scale = (m == 0 ? one(FT) : FT(2)) * w[i] / FT(π)
            az = _azimuthal_kernel(si, 1, m, dϕ) *
                 _azimuthal_kernel(si, 1, m, ϕ[i])
            truncated[si] += scale * az * M[si, 1]
            dtruncated[si] += scale * az * dM[si, 1]
        end
    end
    return exact[:, 1] - truncated, dexact[:, 1] - dtruncated
end

"""
    apply_ss_correction!(R_SFI, surf, pol_type, vza, vaz, μ₀, τ_total, m_max, nSpec;
                         F₀, n_water)

Add exact minus Fourier-truncated Cox–Munk direct-solar reflection at TOA.
`m_max` is the last included Fourier order. `F₀` is the incident Stokes spectrum
with shape `(n_stokes, nSpec)`; the direct-solar correction currently requires
unpolarized illumination. Both downward and upward paths use `τ_total`, the
solver's cumulative optical depth. A zero incident spectrum adds zero glint.

The linearized overload additionally accepts `(Ṙ, τ̇_total, wind_column)`
before `surf`. Its atmospheric columns follow `τ̇_total` and `wind_column`
selects the surface wind-speed column in `Ṙ`.
"""
function apply_ss_correction!(R_SFI::AbstractArray{FT,3},
                               surf::CoxMunkSurface{FT},
                               pol_type, vza, vaz, μ₀::FT,
                               τ_total::AbstractVector{FT},
                               m_max::Int, nSpec::Int;
                               F₀,
                               n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    return apply_ss_correction!(R_SFI, nothing, nothing, 0, surf, pol_type,
        vza, vaz, μ₀, τ_total, m_max, nSpec; F₀, n_water)
end

function apply_ss_correction!(R_SFI::AbstractArray{FT,3}, Ṙ,
                               τ̇_total, wind_column::Int,
                               surf::CoxMunkSurface{FT},
                               pol_type, vza, vaz, μ₀::FT,
                               τ_total::AbstractVector{FT},
                               m_max::Int, nSpec::Int;
                               F₀,
                               n_water::Complex{FT} = _get_n_water(surf, FT(550))) where FT
    _require_unpolarized_solar(F₀, "Cox–Munk direct-solar correction")
    for iv in eachindex(vza)
        μ_v = FT(cosd(vza[iv]))
        residual, dresidual = _coxmunk_glint_residual(
            surf, pol_type.n, μ_v, μ₀, FT(deg2rad(vaz[iv])), m_max, n_water)
        airmass = inv(μ₀) + inv(μ_v)
        for s in 1:nSpec, si in 1:pol_type.n
            atten = μ₀ * F₀[1,s] * exp(-τ_total[s] * airmass)
            correction = atten * residual[si]
            R_SFI[iv,si,s] += correction
            if Ṙ !== nothing
                for p in axes(τ̇_total, 2)
                    Ṙ[iv,si,s,p] -= airmass * τ̇_total[s,p] * correction
                end
                Ṙ[iv,si,s,wind_column] += atten * dresidual[si]
            end
        end
    end
    return nothing
end
