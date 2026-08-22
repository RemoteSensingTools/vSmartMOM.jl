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

`reflectance(rpv, μ, m)` does the azimuthal Fourier integral over Δϕ to
get the m-th moment used in the adding-doubling solver; the generic
`reflectance(::AbstractSurfaceType, pol_type, μ, m)` below performs the
same Gauss-Legendre quadrature for any analytic BRDF.

Polarized RT: only the I → I block is non-zero (RPV is a scalar model),
so for n_stokes > 1 the function returns 0 in the off-diagonal Stokes
slots.  See `rossli_surface.jl` and `coxmunk_surface.jl` for the
polarization-aware analogues.
============================================================================
=#

"""
    create_surface_layer!(brdf::AbstractSurfaceType, added_layer, SFI, m,
                          pol_type, quad_points, τ_sum, architecture)

Generic surface-layer builder for any analytic BRDF (RPV, Ross-Li, …) via
its Fourier-moment `reflectance(brdf, pol_type, μ, m)`. Computes (in place)
the surface optical properties as an [`AddedLayer`](@ref), sharing the
source-fill (`_surface_source!`) and r/t-fill (`_fill_surface_layer!`)
helpers with the Lambertian scaffold (`lambertian_surface.jl`).

    - `brdf` any [`AbstractSurfaceType`](@ref) with a `reflectance` method
    - `SFI` bool if SFI (Source Function Integration) is used
    - `m` Fourier moment (starting at 0; factor-2 normalization applied at m=0)
    - `pol_type` Polarization type struct
    - `quad_points` Quadrature points struct
    - `τ_sum` total optical thickness from TOA to the surface
    - `architecture` Compute architecture (GPU,CPU)
"""
function create_surface_layer!(brdf::AbstractSurfaceType,
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI,
                               m::Int,
                               pol_type,
                               quad_points,
                               τ_sum,
                               architecture;
                               F₀=nothing)

    (; qp_μ, wt_μ, qp_μN, wt_μN) = quad_points
    FT = eltype(qp_μN)
    # Get size of added layer
    Nquad = size(added_layer.r⁻⁺,1) ÷ pol_type.n
    arr_type = array_type(architecture)
    T_surf = arr_type(Diagonal(ones(FT, pol_type.n*Nquad)))

    # Fourier-m reflectance from the analytic BRDF, with the factor-2 m=0
    # normalization (Albedo normalized by π → 1/π·2π for the 0th moment).
    ρ = (m == 0 ? 2 : 1) * vSmartMOM.CoreRT.reflectance(brdf, pol_type, collect(qp_μ), m)
    # Dense reflectance matrix on device.
    R_surf = arr_type(ρ)

    # Source function of surface (m=0 source fill shared with the Lambertian
    # scaffold; uses the attenuated direct beam — see lambertian_surface.jl).
    if SFI
        _surface_source!(added_layer, R_surf, τ_sum, quad_points, pol_type, FT, architecture; F₀=F₀)
    end

    # Quadrature-weight, then fill r/t via the shared scaffold helper.
    _fill_surface_layer!(added_layer, R_surf * Diagonal(qp_μN.*wt_μN), T_surf)

end

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

"""
    reflectance(brdf::AbstractSurfaceType, pol_type, μ, m)

Fourier moment `m` of the BRDF reflectance matrix for quadrature directions `μ`.

Computes ``R_{ij} = (f/\\pi) \\int_0^\\pi \\rho(n, \\mu_i, \\mu_j, \\phi) \\cos(m\\phi) \\, d\\phi``
with `f = 2` for m=0, `f = 1` otherwise. Returns `[nμ·n_stokes, nμ·n_stokes]` matrix.
"""
function reflectance(brdf::AbstractSurfaceType, pol_type, μ::AbstractArray{FT}, m::Int) where FT
    # Hardcoded nQuad for now, needs to go into brdf in the future!
    nQuad = 100

    ty = array_type(architecture(μ))
    # Size of Matrix
    nn = length(μ) * pol_type.n
    Rsurf = ty(zeros(FT,nn,nn)) 
    ff = m==0 ? FT(1.0) : FT(2.0)
    #@show ff
    for n = 1:pol_type.n
        f(x) = reflectance.((brdf,),n, μ, μ', x) * cos(m*x)
        ϕ,w  = CanopyOptics.gauleg(nQuad,   FT(0),  FT(π));
        # more clumsy now with quadrature:
        # TODO clean this up a bit
        b = f.(ϕ)
        c = similar(Rsurf[n:pol_type.n:end,n:pol_type.n:end]) * FT(0)
        for i in eachindex(b)
            c += w[i] * b[i]
        end

        #@show size(c), size(b[1]),  size(reflectance.((brdf,),n, μ, μ', 1.0))
        Rsurf[n:pol_type.n:end,n:pol_type.n:end] .= c/π
        # use quadgk before, was a bit slow 
        # quadgk(f, 0, π, rtol=1e-4)[1] / π
    end
    return ff * Rsurf
end

