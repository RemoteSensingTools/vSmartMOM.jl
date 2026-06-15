#=
============================================================================
Lambertian (isotropic) surface BRDF
============================================================================

A Lambertian surface scatters incident flux equally in every outgoing
direction.  Its bidirectional reflectance is constant:

    ρ(μᵢ, μᵣ, Δϕ)  =  a / π

where `a` ∈ [0, 1] is the spectrally-varying albedo.  Because the BRDF has
no angular structure, only the m = 0 Fourier moment is non-zero — the
`m > 0` branch sets the surface added-layer to identity / zero.  The
factor of 2 in `ρ = 2·a` here is the `1/π · 2π` from converting the
hemispheric flux `a` to a Lambertian radiance × azimuthal-weight
compensation (the same trick is used in `inject_surface_SIF!` below).

Three flavours share the same equation:
  • `LambertianSurfaceScalar`    — single-band scalar albedo
  • `LambertianSurfaceLegendre`  — albedo = Σ cₙ · Pₙ(λ̃) (spectral)
  • `LambertianSurfaceSpline`    — albedo from a spline interpolator

Each fills `added_layer.r⁻⁺ = a/π · μ_quad·w_quad` and leaves `r⁺⁻ = 0`.
The lower boundary is opaque, so `t⁻⁻ = 0`: an upward field entering from
below the surface would represent a sub-surface source, which is outside the
atmospheric RT problem. `t⁺⁺` is kept as the identity because the current
surface-as-layer coupling uses it as a pass-through for the downwelling field
at the surface (e.g. diagnostics/HDRF), not as physical subsurface
transmission.
============================================================================
=#

function _surface_solar_F₀(F₀, FT, pol_type, nSpec::Integer, arr_type)
    if F₀ === nothing
        F₀_host = zeros(FT, pol_type.n, nSpec)
        @views F₀_host[1, :] .= one(FT)
        return arr_type(F₀_host)
    end
    size(F₀) == (pol_type.n, nSpec) || error(
        "Surface solar source F₀ shape $(size(F₀)) does not match " *
        "(pol_type.n, nSpec) = ($(pol_type.n), $nSpec).")
    return arr_type(FT.(F₀))
end

function _surface_beam_at_surface(F₀, FT, pol_type, quad_points, τ_sum, architecture)
    (; qp_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points
    arr_type = array_type(architecture)
    nSpec = length(τ_sum)
    beam = arr_type(zeros(FT, length(qp_μN), nSpec))
    F₀_dev = _surface_solar_F₀(F₀, FT, pol_type, nSpec, arr_type)
    beam[iμ₀Nstart:pol_type.n*iμ₀, :] .= F₀_dev .* exp.(-τ_sum / μ₀)'
    return beam
end

"""
    $(FUNCTIONNAME)(lambertian::LambertianSurfaceScalar{FT})

Computes (in place) surface optical properties for a (scalar) lambertian albedo as [`AddedLayer`](@ref) 

    - `lambertian` a [`LambertianSurfaceScalar`](@ref) struct that defines albedo as scalar
    - `SFI` bool if SFI is used
    - `m` Fourier moment (starting at 0)
    - `pol_type` Polarization type struct
    - `quad_points` Quadrature points struct
    - `τ_sum` total optical thickness from TOA to the surface
    - `architecture` Compute architecture (GPU,CPU)
""" 
function create_surface_layer!(lambertian::LambertianSurfaceScalar{FT}, 
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI,
                               m::Int,
                               pol_type,
                               quad_points,
                               τ_sum,
                               architecture;
                               F₀=nothing) where {FT}
    
    (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points
    j₀⁺ = added_layer.j₀⁺
    j₀⁻ = added_layer.j₀⁻
    # Get size of added layer
    Nquad = size(added_layer.r⁻⁺,1) ÷ pol_type.n
    tmp    = ones(pol_type.n*Nquad)
    arr_type = array_type(architecture)
    T_surf = arr_type(Diagonal(tmp))
    if m == 0
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        ρ = FT(2) * lambertian.albedo#/FT(π)
        
        # Construct dense surface reflectance matrix and move to device
        R_surf = Matrix(Diagonal(vcat(ρ, zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)
        R_surf = arr_type(R_surf)

        
        # Source function of surface:
        if SFI
            beam_at_surface = _surface_beam_at_surface(F₀, FT, pol_type,
                                                       quad_points, τ_sum,
                                                       architecture)
            j₀⁺[:,1,:] .= beam_at_surface;
            j₀⁻[:,1,:] .= μ₀ * (R_surf * beam_at_surface);
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        

        #@show size(added_layer.r⁻⁺), size(R_surf), size(j₀⁻)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= zero(FT);

    else
        added_layer.r⁻⁺ .= zero(FT);
        added_layer.r⁻⁺ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= zero(FT);
        j₀⁺ .= zero(FT);
        j₀⁻ .= zero(FT);
    end
end

function create_surface_layer!(lambertian::LambertianSurfaceLegendre{FT}, 
    added_layer::Union{AddedLayer,AddedLayerRS},
    SFI,
    m::Int,
    pol_type,
    quad_points,
    τ_sum,
    architecture;
    F₀=nothing) where {FT}
    j₀⁺ = added_layer.j₀⁺
    j₀⁻ = added_layer.j₀⁻
    if m == 0
        (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points
        legendre_coeff = lambertian.legendre_coeff
        arr_type = array_type(architecture)
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        # a) Define range for legendre polynomial:
        x = collect(range(FT(-1), FT(1), length=length(τ_sum)));
        # Legendre Polynomial basis functions:
        P = Scattering.compute_legendre_poly(x,length(legendre_coeff))[1]
        # Evaluate Polynomial (as matrix multiplication)
        albedo = P * legendre_coeff
        ρ = arr_type(FT(2) .* albedo)
        # Get size of added layer
        dim = size(added_layer.r⁻⁺)
        Nquad = dim[1] ÷ pol_type.n

        R_surf = Matrix(Diagonal(vcat(FT(1), zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)

        # Move to architecture:
        R_surf = arr_type(R_surf)

        # Source function of surface:
        if SFI
            beam_at_surface = _surface_beam_at_surface(F₀, FT, pol_type,
                                                       quad_points, τ_sum,
                                                       architecture)
            j₀⁺[:] .= zero(FT)
            # Suniti double-check
            j₀⁻[:,1,:] = μ₀*(R_surf*beam_at_surface) .* ρ';
        end
        R_surf   = R_surf * Diagonal(qp_μN.*wt_μN)
        siz = size(added_layer.r⁻⁺)
        R_surf3D = reshape(reduce(hcat,[i*R_surf for i in collect(ρ)]), siz...);
        tmp    = ones(pol_type.n*Nquad)
        T_surf = arr_type(Diagonal(tmp))

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.j₀⁻)
        added_layer.r⁻⁺ .= R_surf3D;
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= zero(FT);

    else
        Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n
        arr_type = array_type(architecture)
        T_surf = arr_type(Diagonal(ones(FT, pol_type.n * Nquad)))
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.t⁺⁺[:] .= T_surf;
        added_layer.t⁻⁻[:] .= zero(FT);
        j₀⁺[:] .= zero(FT);
        j₀⁻[:] .= zero(FT);
    end
end

function create_surface_layer!(lambertian::LambertianSurfaceSpline{FT}, 
    added_layer::Union{AddedLayer,AddedLayerRS},
    SFI,
    m::Int,
    pol_type,
    quad_points,
    τ_sum,
    architecture;
    F₀=nothing) where {FT}
    j₀⁺ = added_layer.j₀⁺
    j₀⁻ = added_layer.j₀⁻
    if m == 0
        (; qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀) = quad_points

        arr_type = array_type(architecture)
        
        # Evaluate spline
        albedo = lambertian.interpolator(lambertian.wlGrid)
        ρ = arr_type(FT(2) .* albedo)
        # Get size of added layer
        dim = size(added_layer.r⁻⁺)
        Nquad = dim[1] ÷ pol_type.n

        R_surf = Matrix(Diagonal(vcat(FT(1), zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)

        # Move to architecture:
        R_surf = arr_type(R_surf)

        # Source function of surface:
        if SFI
            beam_at_surface = _surface_beam_at_surface(F₀, FT, pol_type,
                                                       quad_points, τ_sum,
                                                       architecture)
            j₀⁺[:] .= zero(FT)
            # Suniti double-check
            j₀⁻[:,1,:] = μ₀*(R_surf*beam_at_surface) .* ρ';
        end
        R_surf   = R_surf * Diagonal(qp_μN.*wt_μN)
        
        tmp    = ones(pol_type.n*Nquad)
        T_surf = arr_type(Diagonal(tmp))
        added_layer.r⁻⁺ .= R_surf .* reshape(ρ, 1, 1, :)
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= zero(FT);

    else
        Nquad = size(added_layer.r⁻⁺, 1) ÷ pol_type.n
        arr_type = array_type(architecture)
        T_surf = arr_type(Diagonal(ones(FT, pol_type.n * Nquad)))
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.t⁺⁺[:] .= T_surf;
        added_layer.t⁻⁻[:] .= zero(FT);
        j₀⁺[:] .= zero(FT);
        j₀⁻[:] .= zero(FT);
    end
end

function reflectance(sur::LambertianSurfaceScalar{FT}, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    return sur.albedo
end

"""
    inject_surface_SIF!(brdf, added_layer, m, pol_type, SIF₀, architecture)

Add isotropic solar-induced fluorescence (SIF) surface emission to
`added_layer.j₀⁻`. SIF is Lambertian — only the m=0 Fourier moment carries
it — so higher moments are untouched. The factor 2 comes from
(1/π) × 2π: (1/π) normalizes the hemispheric SIF flux `SIF₀` into a
Lambertian radiance, and 2π compensates the `weight = 0.5/π` azimuthal
weighting applied downstream in `postprocessing_vza!` (SIF is isotropic
and must not be azimuthally weighted).

Ported from sanghavi `lambertian_surface.jl` (injection sites at
sanghavi lines 67-68 and 157-158). Non-Lambertian surfaces fall through
to a no-op.
"""
inject_surface_SIF!(::AbstractSurfaceType, _added_layer, _m, _pol_type, _SIF₀, _architecture) = nothing

function inject_surface_SIF!(
    ::Union{LambertianSurfaceScalar, LambertianSurfaceLegendre, LambertianSurfaceSpline},
    added_layer::Union{AddedLayer, AddedLayerRS},
    m::Int,
    pol_type,
    SIF₀::AbstractArray,
    architecture,
)
    m == 0 || return nothing
    iszero(SIF₀) && return nothing
    FT = eltype(added_layer.j₀⁻)
    arr_type = array_type(architecture)
    Nquad = size(added_layer.j₀⁻, 1) ÷ pol_type.n
    added_layer.j₀⁻[:, 1, :] .+= FT(2) .* arr_type(repeat(FT.(SIF₀), Nquad))
    return nothing
end

