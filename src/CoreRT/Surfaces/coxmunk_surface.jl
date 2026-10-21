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
                               n_water::Complex{FT} = _single_n_water(surf)) where FT
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
                                         n_water::Complex{FT} = _single_n_water(surf)) where FT
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
                                n_water::Complex{FT} = _single_n_water(surf)) where FT
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

# ──────────────────────────────────────────────────────────────────────
# Water refractive index per spectral point
# ──────────────────────────────────────────────────────────────────────

"""
    _resolve_spectral_surface(surface, ν)

The surface as the RT solves it on the spectral grid `ν` (wavenumbers in
cm⁻¹, one per spectral point, in solver order). A `CoxMunkSurface` with
`n_water = nothing` is returned with a vector `n_water` holding the
Segelstein (1981) index at λ = 10⁷/ν nm of every point
([`water_refractive_index`](@ref)); a vector `n_water` must have one entry
per spectral point (`ArgumentError` otherwise); a scalar `n_water` and all
other surfaces are returned unchanged. The RT drivers call this before the
Fourier loop (a canopy resolves its soil per moment), so the surface layers,
the exact-glint correction and the Jacobians all see the same index.
"""
_resolve_spectral_surface(surface::AbstractSurfaceType, ν) = surface

function _resolve_spectral_surface(surf::CoxMunkSurface{FT}, ν) where {FT}
    n_water = surf.n_water
    if n_water === nothing
        λ_nm = 1e7 ./ Float64.(Array(ν))
        return CoxMunkSurface{FT}(wind_speed = surf.wind_speed,
            n_water = Complex{FT}.(water_refractive_index.(λ_nm)),
            whitecap_albedo = surf.whitecap_albedo,
            include_whitecaps = surf.include_whitecaps,
            shadowing = surf.shadowing)
    end
    n_water isa AbstractVector && _check_n_water_length(n_water, length(ν))
    return surf
end

function _check_n_water_length(n_water::AbstractVector, n_spec::Integer)
    length(n_water) == n_spec || throw(ArgumentError(
        "CoxMunkSurface: vector n_water has $(length(n_water)) entries, but the " *
        "band has $n_spec spectral points (one refractive index per spectral point)"))
    return nothing
end

"""
    _spectral_n_water(surf, n_spec)

Water refractive index the surface layer uses on `n_spec` spectral points:
one `Complex{FT}` when it is the same at every point (scalar `n_water`, or a
vector with equal entries), otherwise the per-point vector. A surface with
`n_water = nothing` must first be resolved on the spectral grid
(`_resolve_spectral_surface`).
"""
function _spectral_n_water(surf::CoxMunkSurface{FT}, n_spec::Integer) where {FT}
    n_water = surf.n_water
    n_water === nothing && throw(ArgumentError(
        "CoxMunkSurface(n_water = nothing) reached the surface builder without the " *
        "wavelengths of the spectral points (the RT drivers resolve it with " *
        "`_resolve_spectral_surface`). Pass n_water explicitly here: a scalar, or a " *
        "vector with one entry per spectral point (e.g. from water_refractive_index)."))
    n_water isa Complex && return n_water
    _check_n_water_length(n_water, n_spec)
    return all(==(first(n_water)), n_water) ? first(n_water) : n_water
end

"""
    _single_n_water(surf)

Default `n_water` of the single-wavelength BRDF functions
(`coxmunk_brdf_mueller`, `reflectance`, …): the scalar `n_water`, or the
common value of a vector with equal entries. Other vectors and `nothing`
carry no single index and throw an `ArgumentError`; pass `n_water` for the
wavelength of interest.
"""
function _single_n_water(surf::CoxMunkSurface)
    n_water = surf.n_water
    n_water isa Complex && return n_water
    n_water !== nothing && !isempty(n_water) && all(==(first(n_water)), n_water) &&
        return first(n_water)
    throw(ArgumentError(
        "CoxMunkSurface.n_water = $(n_water === nothing ? "nothing" : "a spectral vector") " *
        "has no single value; pass the n_water keyword for the wavelength of interest " *
        "(e.g. water_refractive_index(λ_nm))"))
end

"""
    _n_water_knots(n_water) -> (knots, lo, hi, t)

Knots for the spectral Cox–Munk operator. The BRDF depends on the spectral
point only through the Fresnel coefficients of `n_water[s]`, which vary
slowly across a band. The operator is evaluated exactly at the knot points
`knots` (indices into `n_water`, first and last point included) and, at every
other point `s`, interpolated linearly between the knots `lo[s]` and `hi[s]`
(indices into `knots`) with weight `t[s]`:

    R(n_water[s]) ≈ (1 − t[s]) R(n_water[knots[lo[s]]]) + t[s] R(n_water[knots[hi[s]]]).

`t[s]` projects `n_water[s]` on the complex segment between the two knot
indices; at knots the operator is exact. A segment between two knots is
accepted when, at every point on it, the same interpolation reproduces the
point's own Fresnel Mueller elements — (|r_s|² ± |r_p|²)/2 and the real and
imaginary parts of r_s r_p* — at the local incidence cosines
`_N_WATER_PROXY_COS` to `_N_WATER_TOL` times the local reflectance
(|r_s|² + |r_p|²)/2. The BRDF is a positive slope-weighted sum of these
elements over local incidence angles, so its intensity element is then
interpolated to about the same relative accuracy (including n close to 1,
where the reflectance is most sensitive). An index with real part ≤ 1 admits
total internal reflection, whose critical angle makes the operator
non-smooth in n; such points are always knots (exact). The check runs in
Float64; in Float32 the tolerance is at least 4 eps(Float32), the rounding of
the stored index. Knots are placed greedily (doubling, then bisection, from
the previous knot); every accepted segment is checked, though it need not be
the longest possible one. Adjacent points always form a segment. A constant
`n_water` gives one knot; an index that is not smooth from point to point
(e.g. noisy measured values) makes most points knots, i.e. one full
reflectance evaluation per spectral point.
"""
function _n_water_knots(n_water::AbstractVector{Complex{FT}}) where {FT}
    N = length(n_water)
    all(==(first(n_water)), n_water) && return [1], ones(Int, N), ones(Int, N), zeros(FT, N)
    tol = max(_N_WATER_TOL, 4 * eps(FT))
    n = ComplexF64.(n_water)
    F = [_fresnel_elements(n[s], c) for c in _N_WATER_PROXY_COS, s in 1:N]
    # Projection weight of point s on the segment a → e.
    function project(a, e, s)
        local Δ = n[e] - n[a]
        local δ = n[s] - n[a]
        return iszero(Δ) ? 0.0 : clamp(real(δ * conj(Δ)) / abs2(Δ), 0.0, 1.0)
    end
    function valid(a, e)
        any(s -> real(n[s]) <= 1, a:e) && return false
        for s in (a + 1):(e - 1)
            w = project(a, e, s)
            for j in axes(F, 1)
                fa, fe, fs = F[j, a], F[j, e], F[j, s]
                bound = tol * fs[1]
                for k in 1:4
                    abs((1 - w) * fa[k] + w * fe[k] - fs[k]) <= bound || return false
                end
            end
        end
        return true
    end
    knots = [1]
    a = 1
    while a < N
        good, probe = a + 1, a + 2
        while probe <= N && valid(a, probe)
            good, probe = probe, a + 2 * (probe - a)
        end
        bad = min(probe, N + 1)
        while bad - good > 1
            mid = (good + bad) ÷ 2
            valid(a, mid) ? (good = mid) : (bad = mid)
        end
        push!(knots, good)
        a = good
    end
    lo = ones(Int, N); hi = ones(Int, N); t = zeros(FT, N)
    for q in 1:(length(knots) - 1)
        a, e = knots[q], knots[q + 1]
        for s in a:e
            lo[s], hi[s] = q, q + 1
            t[s] = s == a ? zero(FT) : s == e ? one(FT) : FT(project(a, e, s))
        end
    end
    return knots, lo, hi, t
end

# The knots depend only on the values of `n_water`; within a run the same
# values reach the surface builder at every Fourier moment (and the glint
# correction), so the plan of the last values is kept, with a copy of those
# values, and reused while an equal vector is passed.
const _N_WATER_KNOT_CACHE = Ref{Any}(nothing)
const _N_WATER_KNOT_LOCK = ReentrantLock()

"`_n_water_knots(n_water)`, reused while vectors equal to the last one are passed. The returned arrays are shared with the cache: do not mutate them."
function _n_water_knot_plan(n_water::AbstractVector{Complex{FT}}) where {FT}
    lock(_N_WATER_KNOT_LOCK) do
        cached = _N_WATER_KNOT_CACHE[]
        if !(cached !== nothing && cached[1] isa Vector{Complex{FT}} && isequal(cached[1], n_water))
            cached = (Vector(n_water), _n_water_knots(n_water))
            _N_WATER_KNOT_CACHE[] = cached
        end
        return cached[2]::Tuple{Vector{Int}, Vector{Int}, Vector{Int}, Vector{FT}}
    end
end

"Fresnel Mueller elements ((|r_s|²+|r_p|²)/2, (|r_s|²−|r_p|²)/2, Re r_s r_p*, Im r_s r_p*) at incidence cosine `c`."
function _fresnel_elements(n::ComplexF64, c::Float64)
    r_s, r_p = fresnel_coefficients(n, c)
    r_sp = r_s * conj(r_p)
    return ((abs2(r_s) + abs2(r_p)) / 2, (abs2(r_s) - abs2(r_p)) / 2, real(r_sp), imag(r_sp))
end

# Interpolation tolerance of `_n_water_knots` (relative to the local Fresnel
# reflectance) and the local incidence cosines where it is checked;
# test_coxmunk_spectral.jl checks the resulting operator error across
# 0.2–4 µm and for other indices.
const _N_WATER_TOL = 1e-7
const _N_WATER_PROXY_COS = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.03)

"""
    _interpolate_knots!(R, Rk, lo, hi, t)

`R[:, :, s] = (1 − t[s]) Rk[:, :, lo[s]] + t[s] Rk[:, :, hi[s]]` for every
spectral point `s` (see `_n_water_knots`). All arguments live on the
device of `R`.
"""
function _interpolate_knots!(R::AbstractArray{FT,3}, Rk::AbstractArray{FT,3},
                             lo::AbstractVector{<:Integer}, hi::AbstractVector{<:Integer},
                             t::AbstractVector{FT}) where {FT}
    device = devi(architecture(R))
    kernel! = _interpolate_knots_kernel!(device)
    kernel!(R, Rk, lo, hi, t; ndrange = size(R))
    synchronize_if_gpu()
    return R
end

@kernel function _interpolate_knots_kernel!(R, @Const(Rk), @Const(lo), @Const(hi), @Const(t))
    i, j, s = @index(Global, NTuple)
    @inbounds begin
        w = t[s]
        R[i, j, s] = (one(w) - w) * Rk[i, j, lo[s]] + w * Rk[i, j, hi[s]]
    end
end

"""
    _spectral_knot_operators(f, n_water, arr_type) -> (Rk..., lo, hi, t)

Evaluate `f(n)`, which returns one host matrix or a tuple of host matrices of
equal size, at the knots of `n_water` (see `_n_water_knots`) and
upload the stacked `(rows, cols, nKnots)` arrays and the interpolation
indices and weights with `arr_type`.
"""
function _spectral_knot_operators(f, n_water::AbstractVector{Complex{FT}}, arr_type) where {FT}
    knots, lo, hi, t = _n_water_knot_plan(n_water)
    values = [f(n_water[k]) for k in knots]
    stacked = values[1] isa Tuple ?
        ntuple(i -> arr_type(stack(v[i] for v in values)), length(values[1])) :
        (arr_type(stack(values)),)
    return (stacked..., arr_type(lo), arr_type(hi), arr_type(t))
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
                     n_water::Complex{FT} = _single_n_water(surf)) where FT
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
# Surface layer
# ──────────────────────────────────────────────────────────────────────

"""
    create_surface_layer!(surf::CoxMunkSurface, added_layer, SFI, m, pol_type,
                          quad_points, τ_sum, architecture; F₀ = nothing)

Cox–Munk surface layer for Fourier moment `m`. With one refractive index for
the whole band (scalar `n_water`, or a vector with equal entries) this is the
shared analytic-surface assembly: one reflectance matrix for every spectral
point. A spectrally varying `n_water` (one entry per spectral point; resolve
`n_water = nothing` first with `_resolve_spectral_surface`) gives a
reflectance matrix per spectral point: the Fourier reflectance is evaluated
at the knots of `_n_water_knots` and interpolated on the device to
the other points, and the reflected direct beam `j₀⁻` uses each point's
matrix.
"""
function create_surface_layer!(surf::CoxMunkSurface{FT},
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI,
                               m::Int,
                               pol_type,
                               quad_points,
                               τ_sum,
                               architecture;
                               F₀=nothing) where {FT}
    (; qp_μ, qp_μN, wt_μN) = quad_points
    n_spec = length(τ_sum)
    Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n
    arr_type = array_type(architecture)
    T_surf = arr_type(Diagonal(ones(FT, pol_type.n * Nquad)))
    μ = collect(qp_μ)
    # The surface-layer convention carries an additional factor two at m=0.
    fac = m == 0 ? 2 : 1
    n_water = _spectral_n_water(surf, n_spec)

    if n_water isa Complex
        # Same operations as the analytic-surface assembly (analytic_surface.jl).
        ρ = fac * reflectance(surf, pol_type, μ, m; n_water)
        R_surf = arr_type(ρ)
        if SFI
            _surface_source!(added_layer, R_surf, τ_sum, quad_points, pol_type,
                             FT, architecture; F₀)
        end
        _fill_surface_layer!(added_layer, R_surf * Diagonal(qp_μN .* wt_μN), T_surf)
        return nothing
    end

    Rk, lo, hi, t = _spectral_knot_operators(
        n -> fac * reflectance(surf, pol_type, μ, m; n_water = n), n_water, arr_type)
    # Unweighted per-point reflectance first (for the direct-beam source),
    # then the quadrature-weight column scaling, in place.
    R_surf = added_layer.r⁻⁺
    _interpolate_knots!(R_surf, Rk, lo, hi, t)
    if SFI
        quad_points.external_solar && throw(ArgumentError(
            "external-solar SFI currently supports Lambertian surfaces only"))
        (; μ₀) = quad_points
        beam = _surface_beam_at_surface(F₀, FT, pol_type, quad_points, τ_sum,
                                        architecture)
        added_layer.j₀⁺[:, 1, :] .= beam
        added_layer.j₀⁻ .= μ₀ .* batched_mul(R_surf, reshape(beam, size(beam, 1), 1, n_spec))
    end
    R_surf .*= reshape(qp_μN .* wt_μN, 1, :)
    added_layer.r⁺⁻ .= zero(FT)
    added_layer.t⁺⁺ .= T_surf
    added_layer.t⁻⁻ .= zero(FT)
    return nothing
end

# ──────────────────────────────────────────────────────────────────────
# Single-scattering (TMS) correction for the specular hotspot
# ──────────────────────────────────────────────────────────────────────

"""
    _coxmunk_glint_truncated(surf, n, μ_v, μ₀, dϕ, m_max, n_water)

Fourier series through order `m_max` of the direct-solar BRDF column
(incident `μ₀`, reflected `μ_v` at relative azimuth `dϕ`) and of its wind
tangent, built from the same 100-node azimuth quadrature as the surface
operator. The surface builder doubles the zeroth moment and postprocessing
halves it; there is no additional half-weight in this physical BRDF
reconstruction.
"""
function _coxmunk_glint_truncated(surf::CoxMunkSurface{FT}, n, μ_v, μ₀,
                                  dϕ, m_max, n_water) where FT
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
    return truncated, dtruncated
end

"""
    _coxmunk_glint_residual(surf, n, μ_v, μ₀, dϕ, m_max, n_water)

Exact minus truncated direct-solar BRDF column and its analytic wind tangent
(see `_coxmunk_glint_truncated`).
"""
function _coxmunk_glint_residual(surf::CoxMunkSurface{FT}, n, μ_v, μ₀,
                                  dϕ, m_max, n_water) where FT
    exact, dexact = coxmunk_brdf_mueller_and_deriv(
        surf, n, μ₀, μ_v, dϕ; n_water)
    truncated, dtruncated = _coxmunk_glint_truncated(
        surf, n, μ_v, μ₀, dϕ, m_max, n_water)
    return exact[:, 1] - truncated, dexact[:, 1] - dtruncated
end

"""
    apply_ss_correction!(R_SFI, surf, pol_type, vza, vaz, μ₀, τ_total, m_max, nSpec;
                         F₀, n_water = nothing)

Add exact minus Fourier-truncated Cox–Munk direct-solar reflection at TOA.
`m_max` is the last included Fourier order. `F₀` is the incident Stokes spectrum
with shape `(n_stokes, nSpec)`; the direct-solar correction currently requires
unpolarized illumination. Both downward and upward paths use `τ_total`, the
solver's cumulative optical depth. A zero incident spectrum adds zero glint.

The refractive index is the surface's own, per spectral point (`n_water =
nothing`; a surface with `surf.n_water = nothing` must first be resolved with
`_resolve_spectral_surface`, as the RT drivers do); a `Complex`
keyword applies one index to every point. With a spectrally varying index the
exact term uses the index of each point and the truncated term is
interpolated between the knots of the surface operator
(`_n_water_knots`), so that it cancels the Fourier-series glint the
solver carried: the corrected single-bounce glint is exact at every point.

The correction is applied at TOA only. Interior upwelling levels (`obs_alt`
strictly inside the atmosphere) carry the Fourier-series glint, in the
forward and in the linearized solver.

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
                               n_water::Union{Nothing,Complex{FT}} = nothing) where FT
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
                               n_water::Union{Nothing,Complex{FT}} = nothing) where FT
    _require_unpolarized_solar(F₀, "Cox–Munk direct-solar correction")
    n_spec_water = n_water === nothing ? _spectral_n_water(surf, nSpec) : n_water
    n_spec_water isa Complex || return _apply_ss_correction_spectral!(
        R_SFI, Ṙ, τ̇_total, wind_column, surf, pol_type, vza, vaz, μ₀,
        τ_total, m_max, nSpec, F₀, n_spec_water)
    for iv in eachindex(vza)
        μ_v = FT(cosd(vza[iv]))
        residual, dresidual = _coxmunk_glint_residual(
            surf, pol_type.n, μ_v, μ₀, FT(deg2rad(vaz[iv])), m_max, n_spec_water)
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

# Spectrally varying refractive index: exact glint at every point, truncated
# series interpolated between the operator knots (see `apply_ss_correction!`).
function _apply_ss_correction_spectral!(R_SFI, Ṙ, τ̇_total, wind_column,
                                        surf::CoxMunkSurface{FT}, pol_type, vza, vaz,
                                        μ₀::FT, τ_total, m_max, nSpec, F₀,
                                        n_water::AbstractVector{Complex{FT}}) where FT
    n = pol_type.n
    knots, lo, hi, t = _n_water_knot_plan(n_water)
    for iv in eachindex(vza)
        μ_v = FT(cosd(vza[iv]))
        dϕ = FT(deg2rad(vaz[iv]))
        truncated = [_coxmunk_glint_truncated(surf, n, μ_v, μ₀, dϕ, m_max, n_water[k])
                     for k in knots]
        airmass = inv(μ₀) + inv(μ_v)
        for s in 1:nSpec
            exact, dexact = coxmunk_brdf_mueller_and_deriv(
                surf, n, μ₀, μ_v, dϕ; n_water = n_water[s])
            (T_lo, dT_lo), (T_hi, dT_hi) = truncated[lo[s]], truncated[hi[s]]
            w = t[s]
            atten = μ₀ * F₀[1,s] * exp(-τ_total[s] * airmass)
            for si in 1:n
                residual = exact[si, 1] - ((one(FT) - w) * T_lo[si] + w * T_hi[si])
                correction = atten * residual
                R_SFI[iv,si,s] += correction
                if Ṙ !== nothing
                    dresidual = dexact[si, 1] -
                                ((one(FT) - w) * dT_lo[si] + w * dT_hi[si])
                    for p in axes(τ̇_total, 2)
                        Ṙ[iv,si,s,p] -= airmass * τ̇_total[s,p] * correction
                    end
                    Ṙ[iv,si,s,wind_column] += atten * dresidual
                end
            end
        end
    end
    return nothing
end
