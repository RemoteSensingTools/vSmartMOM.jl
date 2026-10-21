#=
============================================================================
RPV (Rahman–Pinty–Verstraete) BRDF surface
============================================================================

A semi-empirical bidirectional reflectance for vegetated and bare-soil
canopies (Rahman, Pinty & Verstraete, JGR 1993).  The full BRDF is the
product of four physically-named factors:

    ρ(μᵢ, μᵣ, Δϕ) =  ρ₀ · M(μᵢ, μᵣ, k) · F(Θ, cos g) · H(ρ_c, G)
                     ↑     ↑               ↑            ↑
                     |     |               |            |
                     |     Minnaert        Henyey-      "bowl-shape"
                     overall   limb-       Greenstein   geometric
                     amplitude darkening   hot-spot     correction

Geometric quantities (computed in `reflectance(::rpvSurfaceScalar, …)`):
    θᵢ = acos(μᵢ),   θᵣ = acos(μᵣ)        viewing/illumination polar angles
    cos g = -μᵢμᵣ + sin θᵢ sin θᵣ cos Δϕ   scattering (phase) angle
    G = √(tan²θᵢ + tan²θᵣ + 2 tanθᵢ tanθᵣ cosΔϕ)   geometric scale

Free parameters (in `rpvSurfaceScalar`):
    ρ₀    overall albedo                 ρ_c    bowl-shape amplitude
    k     Minnaert exponent              Θ     hot-spot width

`reflectance(rpv, μ, m)` uses the common analytic-surface quadrature in
`analytic_surface.jl` to obtain the Fourier moments used by adding-doubling.

Polarized RT: only the I → I block is non-zero (RPV is a scalar model),
so for n_stokes > 1 the function returns 0 in the off-diagonal Stokes
slots.  See `rossli_surface.jl` and `coxmunk_surface.jl` for the
polarization-aware analogues.
============================================================================
=#

"""
    reflectance(rpv::rpvSurfaceScalar, n, μᵢ, μᵣ, dϕ)

RPV (Rahman-Pinty-Verstraete) BRDF model for scalar (n=1) reflectance.

The RPV model (Rahman et al., 1993) parameterizes the bidirectional reflectance as:
``\\rho = \\rho_0 \\, M(\\mu_i, \\mu_r, k) \\, F(\\Theta, \\cos g) \\, H(\\rho_c, G)``

- **ρ₀**: Overall amplitude (isotropic scaling).
- **k**: Minnaert limb-darkening exponent (controls angular distribution).
- **Θ**: Hot-spot parameter (controls backscatter peak width).
- **ρ_c**: Geometric term amplitude (controls bowl shape).

For polarized RT (n>1), returns zero. See Rahman, Pinty & Verstraete (1993), JGR.
"""
function reflectance(rpv::rpvSurfaceScalar{FT},  n, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    (; ρ₀, ρ_c, k, Θ) = rpv
    # Convert cosines to angles for RPV formula
    if n==1
        θᵢ   = acos(clamp(μᵢ, FT(-1), FT(1))) #assert 0<=θᵢ<=π/2 (ulp guard)
        θᵣ   = acos(clamp(μᵣ, FT(-1), FT(1))) #assert 0<=θᵣ<=π/2
        cosg = -μᵢ*μᵣ + sin(θᵢ)*sin(θᵣ)*cos(dϕ) #RAMI form: μᵢ*μᵣ + sin(θᵢ)*sin(θᵣ)*cos(dϕ) (vSmartMOM sign convention is compatible with that of Rahman, Pinty, Verstraete, 1993) 
        # G² ≥ 0 in ℝ but rounds negative in F32 on the tanθᵢ≈tanθᵣ,
        # cosΔφ≈−1 diagonal → NaN (the round-1 review's non-finite RPV)
        G    = sqrt(max(tan(θᵢ)^2 + tan(θᵣ)^2 + 2*tan(θᵢ)*tan(θᵣ)*cos(dϕ), FT(0))) #RAMI form: (tan(θᵢ)^2 + tan(θᵣ)^2 - 2*tan(θᵢ)*tan(θᵣ)*cos(dϕ))^FT(0.5)
        return ρ₀ * rpvM(μᵢ, μᵣ, k) * rpvF(Θ, cosg) * rpvH(ρ_c, G)
    else
        return FT(0)
    end
end

"""Minnaert term: ``M = (\\mu_i \\mu_r)^{k-1} / (\\mu_i + \\mu_r)^{1-k}``"""
function rpvM(μᵢ::FT, μᵣ::FT, k::FT) where FT
    return (μᵢ * μᵣ)^(k -1) /  (μᵢ + μᵣ)^(1 - k)
end

"""Geometric term: ``H = 1 + (1 - \\rho_c) / (1 + G)``, with ``G`` the phase angle."""
function rpvH(ρ_c::FT, G::FT) where FT
    return 1 + (1 - ρ_c) / (1 + G)
end

"""Hot-spot term: ``F = (1 - \\Theta^2) / (1 + \\Theta^2 + 2\\Theta \\cos g)^{1.5}``"""
function rpvF(θ::FT, cosg::FT) where FT
    θ = -θ #for RAMI only
    return (1 - θ^2) /  (1 + θ^2 + 2θ * cosg)^FT(1.5) #RAMI form: (1 - Θ^2) /  (1 + Θ^2 + 2Θ * cosg)^FT(1.5)
end
