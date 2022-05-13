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
                               τ_sum, F₀,
                               architecture) where {FT}
    
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀ = quad_points
    # Get size of added layer
    Nquad = size(added_layer.r⁻⁺,1) ÷ pol_type.n
    tmp    = ones(pol_type.n*Nquad)
    arr_type = array_type(architecture)
    T_surf = arr_type(Diagonal(tmp))
    if m == 0
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        ρ = 2lambertian.albedo#/FT(π)
        
        
        
        R_surf = Array(Diagonal(vcat(ρ, zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)
        
        # Move to architecture:
        R_surf = arr_type(R_surf)

        
        # Source function of surface:
        if SFI
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .=0;
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            
            added_layer.J₀⁺[iμ₀Nstart:pol_type.n*iμ₀,1,:] .= 
                F₀ .* exp.(-τ_sum/μ₀)';
            added_layer.J₀⁻[iμ₀Nstart:pol_type.n*iμ₀,1,:] = 
                μ₀*(R_surf[iμ₀Nstart:pol_type.n*iμ₀, iμ₀Nstart:pol_type.n*iμ₀]*
                F₀) .* exp.(-τ_sum/μ₀)';
            # added_layer.J₀⁺[:,1,:] .= I₀_NquadN .* exp.(-τ_sum/μ₀)';
            # added_layer.J₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)';
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.J₀⁻)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

    else
        added_layer.r⁻⁺ .= 0;
        added_layer.r⁻⁺ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;
        added_layer.J₀⁺ .= 0;
        added_layer.J₀⁻ .= 0;
    end
end

function create_surface_layer!(lambertian::LambertianSurfaceLegendre{FT}, 
    added_layer::Union{AddedLayer,AddedLayerRS},
    SFI,
    m::Int,
    pol_type,
    quad_points,
    τ_sum, F₀,
    architecture) where {FT}
    FT2 = Float64
    if m == 0
        @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀ = quad_points
        legendre_coeff = lambertian.legendre_coeff
        arr_type = array_type(architecture)
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        # a) Define range for legendre polymonial:
        x = collect(range(FT2(-1), FT2(1), length=length(τ_sum)));
        # Legendre Polynomial basis functions:
        P = Scattering.compute_legendre_poly(x,length(legendre_coeff))[1]
        # Evaluate Polynomial (as matrix multiplication)
        albedo = P * legendre_coeff
        ρ = arr_type(2albedo)
        # Get size of added layer
        dim = size(added_layer.r⁻⁺)
        Nquad = dim[1] ÷ pol_type.n

        R_surf = Array(Diagonal(vcat(FT(1), zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)

        # Move to architecture:
        R_surf = arr_type(R_surf)

        # Source function of surface:
        if SFI
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .=0;
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            added_layer.J₀⁺[:] .= 0
            # Suniti double-check
            # added_layer.J₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* (ρ .* exp.(-τ_sum/μ₀))';
            added_layer.J₀⁻[iμ₀Nstart:pol_type.n*iμ₀,1,:] = 
                μ₀*(R_surf[iμ₀Nstart:pol_type.n*iμ₀, iμ₀Nstart:pol_type.n*iμ₀]*F₀) .* (ρ .* exp.(-τ_sum/μ₀))';
        end
        R_surf   = R_surf * Diagonal(qp_μN.*wt_μN)
        siz = size(added_layer.r⁻⁺)
        R_surf3D = reshape(reduce(hcat,[i*R_surf for i in Array(ρ)]), siz...);
        tmp    = ones(pol_type.n*Nquad)
        T_surf = arr_type(Diagonal(tmp))

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.J₀⁻)
        added_layer.r⁻⁺ .= R_surf3D;
        added_layer.r⁺⁻ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

    else
        added_layer.r⁻⁺[:] .= 0;
        added_layer.r⁻⁺[:] .= 0;
        added_layer.t⁺⁺[:] .= 0;
        added_layer.t⁻⁻[:] .= 0;
        added_layer.J₀⁺[:] .= 0;
        added_layer.J₀⁻[:] .= 0;
    end
end