#=
Fresnel reflection utilities for polarized surface models.

Provides:
- `fresnel_coefficients`: complex amplitude reflection coefficients (r_s, r_p)
- `fresnel_mueller`: 4×4 (or n×n) Mueller matrix for Fresnel reflection
- `stokes_rotation_matrix`: rotation matrix L(φ) for Stokes reference plane

All functions use StaticArrays.SMatrix for performance in inner loops.
=#

"""
    fresnel_coefficients(n_rel, cos_θᵢ)

Complex Fresnel amplitude reflection coefficients for s- and p-polarization.

# Arguments
- `n_rel::Complex`: relative refractive index (n_transmitted / n_incident)
- `cos_θᵢ::Real`: cosine of the incidence angle

# Returns
- `(r_s, r_p)`: complex amplitude coefficients
"""
function fresnel_coefficients(n_rel::Complex{FT}, cos_θᵢ::FT) where FT
    sin_θᵢ² = max(FT(0), one(FT) - cos_θᵢ^2)
    # Snell's law: n_i sin(θ_i) = n_t sin(θ_t)  ⟹  cos(θ_t) = √(1 - (sin θ_i / n_rel)²)
    cos_θₜ = sqrt(one(Complex{FT}) - sin_θᵢ² / n_rel^2)

    r_s = (cos_θᵢ - n_rel * cos_θₜ) / (cos_θᵢ + n_rel * cos_θₜ)
    r_p = (n_rel * cos_θᵢ - cos_θₜ) / (n_rel * cos_θᵢ + cos_θₜ)
    return r_s, r_p
end

"""
    fresnel_mueller(r_s, r_p, n_stokes)

Mueller matrix for Fresnel reflection given complex amplitude coefficients.

For `n_stokes = 4` (IQUV):
```
M = [(|rs|²+|rp|²)/2   (|rs|²-|rp|²)/2    0              0           ]
    [(|rs|²-|rp|²)/2   (|rs|²+|rp|²)/2    0              0           ]
    [0                  0                   Re(rs·rp*)     Im(rs·rp*)  ]
    [0                  0                  -Im(rs·rp*)     Re(rs·rp*)  ]
```
"""
function fresnel_mueller(r_s::Complex{FT}, r_p::Complex{FT}, n_stokes::Int) where FT
    rs2 = abs2(r_s)
    rp2 = abs2(r_p)

    if n_stokes == 1
        return SMatrix{1,1,FT}((rs2 + rp2) / 2)
    end

    rsp = r_s * conj(r_p)
    re_rsp = real(rsp)
    im_rsp = imag(rsp)

    half_sum  = (rs2 + rp2) / 2
    half_diff = (rs2 - rp2) / 2

    if n_stokes == 2
        return SMatrix{2,2,FT}(
            half_sum,  half_diff,
            half_diff, half_sum
        )
    elseif n_stokes == 3
        return SMatrix{3,3,FT}(
            half_sum,  half_diff, FT(0),
            half_diff, half_sum,  FT(0),
            FT(0),     FT(0),    re_rsp
        )
    else  # n_stokes == 4
        return SMatrix{4,4,FT}(
            half_sum,  half_diff, FT(0),    FT(0),
            half_diff, half_sum,  FT(0),    FT(0),
            FT(0),     FT(0),    re_rsp,  -im_rsp,
            FT(0),     FT(0),    im_rsp,   re_rsp
        )
    end
end

"""
    stokes_rotation_matrix(φ, n_stokes)

Rotation matrix L(φ) for the Stokes reference plane.

Rotates the Stokes vector reference plane by angle φ (radians).
Convention matches vSmartMOM's D = [1, 1, -1, -1]:
```
L = [1    0       0      0]
    [0   cos2φ   sin2φ   0]
    [0  -sin2φ   cos2φ   0]
    [0    0       0      1]
```
"""
function stokes_rotation_matrix(φ::FT, n_stokes::Int) where FT
    if n_stokes == 1
        return SMatrix{1,1,FT}(one(FT))
    end

    c2 = cos(2φ)
    s2 = sin(2φ)

    if n_stokes == 2
        return SMatrix{2,2,FT}(
            one(FT), FT(0),
            FT(0),   c2
        )
    elseif n_stokes == 3
        return SMatrix{3,3,FT}(
            one(FT), FT(0), FT(0),
            FT(0),   c2,   -s2,
            FT(0),   s2,    c2
        )
    else  # n_stokes == 4
        return SMatrix{4,4,FT}(
            one(FT), FT(0), FT(0), FT(0),
            FT(0),   c2,   -s2,    FT(0),
            FT(0),   s2,    c2,    FT(0),
            FT(0),   FT(0), FT(0), one(FT)
        )
    end
end
