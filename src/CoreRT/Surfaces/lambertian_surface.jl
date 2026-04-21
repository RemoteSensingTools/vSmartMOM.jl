#=

This file specifies how to create surface layers, given the surface type, and related info

=#

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
                               architecture) where {FT}
    
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
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .= zero(FT);
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            
            j₀⁺[:,1,:] .= I₀_NquadN .* exp.(-τ_sum/μ₀)';
            j₀⁻[:,1,:] .= μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)';
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        

        #@show size(added_layer.r⁻⁺), size(R_surf), size(j₀⁻)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

    else
        added_layer.r⁻⁺ .= zero(FT);
        added_layer.r⁻⁺ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;
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
    architecture) where {FT}
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
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .= zero(FT);
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            j₀⁺[:] .= zero(FT)
            # Suniti double-check
            j₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* (ρ .* exp.(-τ_sum/μ₀))';
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
        added_layer.t⁻⁻ .= T_surf;

    else
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.t⁺⁺[:] .= zero(FT);
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
    architecture) where {FT}
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
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .= zero(FT);
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            j₀⁺[:] .= zero(FT)
            # Suniti double-check
            j₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* (ρ .* exp.(-τ_sum/μ₀))';
        end
        R_surf   = R_surf * Diagonal(qp_μN.*wt_μN)
        
        tmp    = ones(pol_type.n*Nquad)
        T_surf = arr_type(Diagonal(tmp))
        added_layer.r⁻⁺ .= R_surf .* reshape(ρ, 1, 1, :)
        added_layer.r⁺⁻ .= zero(FT);
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

    else
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.r⁻⁺[:] .= zero(FT);
        added_layer.t⁺⁺[:] .= zero(FT);
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

"""
    _sif_source(RS_type)

Return `RS_type.SIF₀` if the field is declared, else `nothing`. Used by
`rt_run` / `rt_run_ss` to thread SIF into `inject_surface_SIF!` without
requiring every `AbstractRamanType` concrete to carry the field (e.g.
`_plus` variants have it commented out).
"""
_sif_source(RS_type) = hasproperty(RS_type, :SIF₀) ? RS_type.SIF₀ : nothing



