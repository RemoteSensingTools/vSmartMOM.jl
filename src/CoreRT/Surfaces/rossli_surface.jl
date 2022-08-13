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
function create_surface_layer!(rpv::RossLiSurfaceScalar{FT}, 
                               added_layer::Union{AddedLayer,AddedLayerRS},
                               SFI,
                               m::Int,
                               pol_type,
                               quad_points,
                               τ_sum,
                               architecture) where {FT}
    
    @unpack qp_μ, wt_μ, qp_μN, wt_μN, iμ₀Nstart, iμ₀, μ₀ = quad_points
    # Get size of added layer
    Nquad = size(added_layer.r⁻⁺,1) ÷ pol_type.n
    tmp    = ones(pol_type.n*Nquad)
    arr_type = array_type(architecture)
    T_surf = arr_type(Diagonal(tmp))
    #if m == 0 
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        #ρ = 2lambertian.albedo#/FT(π)
        
        R_surf = Array(Diagonal(vcat(1, zeros(FT,pol_type.n-1))))
        R_surf = repeat(R_surf',Nquad)
        R_surf = repeat(R_surf',Nquad)
        
        # Move to architecture:
        R_surf = arr_type(R_surf)

        
        # Source function of surface:
        if SFI
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .=0;
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            
            added_layer.j₀⁺[:,1,:] .= I₀_NquadN .* exp.(-τ_sum/μ₀)';
            added_layer.j₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)';
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.j₀⁻)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

    """
    else
        added_layer.r⁻⁺ .= 0;
        added_layer.r⁻⁺ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;
        added_layer.j₀⁺ .= 0;
        added_layer.j₀⁻ .= 0;
    end
    """
end

# Ross Li
function reflectance(RossLi::RossLiSurfaceScalar{FT}, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    @unpack fiso, fvol, fgeo = RossLi
    # TODO: Suniti, stupid calculations here:
    θᵢ   = acos(μᵢ) #assert 0<=θᵢ<=π/2
    θᵣ   = acos(μᵣ) #assert 0<=θᵣ<=π/2
    
    return fiso * RossLi_K_iso() + 
            fvol * RossLi_K_vol(Θᵢ, θᵣ, dϕ) + 
            fgeo * RossLi_K_geo(Θᵢ, θᵣ, dϕ)
end
# Isotropic kernel
function RossLi_K_iso() where FT
    return 1.0;
end

# Volumetric (Ross Thick) kernel
function RossLi_K_vol(θᵢ::FT, θᵣ::FT, dϕ::FT) where FT
    ξ = acos(cos(θᵢ)*cos(θᵣ) + sin(θᵢ)*sin(θᵣ)*cos(dϕ))
    K_vol = ((π/2 - ξ)*cos(ξ) + sin(ξ))/(cos(θᵢ)+cos(θᵣ)) - (π/4)
    return K_vol
end

#Geometric (Li Sparse) kernel
function RossLi_K_geo(θᵢ::FT, θᵣ::FT, dϕ::FT) where FT
    h_by_b = 2 #for RAMI
    b_by_r = 1 #for RAMI

    θᵢᵖ = atan(tan(θᵢ)*b_by_r)
    θᵣᵖ = atan(tan(θᵣ)*b_by_r)

    ξᵖ = acos(cos(θᵢᵖ)*cos(θᵣᵖ) + sin(θᵢᵖ)*sin(θᵣᵖ)*cos(dϕ))
    D = sqrt(tan^2(θᵢᵖ) + tan^2(θᵣᵖ) - 2tan(θᵢᵖ)*tan(θᵣᵖ)*cos(dϕ))
    t = acos(h_by_b * sqrt(D^2 + (tan(θᵢᵖ)*tan(θᵣᵖ)*sin(dϕ))^2)/(sec(θᵢᵖ)+sec(θᵣᵖ)))
    O = (1/π)*(t-sin(t)*cos(t))*(sec(θᵢᵖ)+sec(θᵣᵖ))
    return O - (sec(θᵢᵖ)+sec(θᵣᵖ)) + 0.5*(1+cos(ξᵖ))*sec(θᵢᵖ)*sec(θᵣᵖ)
end

function reflectance(RossLi::AbstractSurfaceType, pol_type, μ::AbstractArray{FT}, m::Int) where FT
    ty = array_type(architecture(μ))
    # Size of Matrix
    nn = length(μ) * pol_type.n
    Rsurf = ty(zeros(FT,nn,nn)) 
    for n = 1:pol_type.n
        f(x) = reflectance.((rpv,),n, μ, μ', x) * cos(m*x)
        Rsurf[n:pol_type.n:end,n:pol_type.n:end] .= quadgk(f, 0, π, rtol=1e-6)[1] / π
    end
    return Rsurf
end


@kernel function applyExpansion!(Rsurf,n_stokes::Int,  v)
    i, j = @index(Global, NTuple)
    # get indices:
    ii = (i-1)*n_stokes 
    jj = (j-1)*n_stokes 
    
    # Fill values:
    for i_n = 1:n_stokes
        Rsurf[ii+i_n, jj+i_n] = v[i,j][i_n]
    end
end

# 2D expansion:
function expandSurface!(Rsurf::AbstractArray{FT,2}, n_stokes::Int, v) where {FT}
    device = devi(architecture(Rsurf))
    applyExpansion_! = applyExpansion!(device)
    event = applyExpansion_!(Rsurf, n_stokes, v, ndrange=size(v));
    wait(device, event);
    synchronize_if_gpu();
    return nothing
end


