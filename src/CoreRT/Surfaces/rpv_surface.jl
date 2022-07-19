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
function create_surface_layer!(rpv::rpvSurfaceScalar{FT}, 
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
            
            added_layer.j₀⁺[:,1,:] .= I₀_NquadN .* exp.(-τ_sum/μ₀)';
            added_layer.j₀⁻[:,1,:] = μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)';
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.j₀⁻)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;

    else
        added_layer.r⁻⁺ .= 0;
        added_layer.r⁻⁺ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;
        added_layer.j₀⁺ .= 0;
        added_layer.j₀⁻ .= 0;
    end
end

function reflectance(rpv::rpvSurfaceScalar{FT}, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    @unpack ρ₀, ρ_c, k, Θ = rpv
    # TODO: Suniti, stupid calculations here:
    θᵢ   = acos(μᵢ)
    θᵣ   = acos(μᵣ)
    cosg = μᵢ*μᵣ + sin(θᵢ)*sin(θᵣ)*cos(dϕ)
    G    = (tan(θᵢ)^2 + tan(θᵣ)^2 - 2*tan(θᵢ)*tan(θᵣ)*cos(dϕ))^FT(0.5)
    ρ₀ * rpvM(μᵢ, μᵣ, k) * rpvF(Θ, cosg) * rpvH(ρ_c, G)
end

function rpvM(μᵢ::FT, μᵣ::FT, k::FT) where FT
    (μᵢ * μᵣ)^(k -1) /  (μᵢ + μᵣ)^(1 -k)
end

function rpvH(ρ_c::FT, G::FT) where FT
    1 + (1 - ρ_c) / (1 + G)
end

function rpvF(Θ::FT, cosg::FT) where FT
    (1 - Θ^2) /  (1 + Θ^2 + 2Θ * cosg)^FT(1.5)
end

function reflectance(rpv::rpvSurfaceScalar{FT}, μ::Array{FT}, m::Int) where FT
    @unpack ρ₀, ρ_c, k, Θ = rpv
    f(x) = reflectance.([rpv], μ, μ', [x]) * cos(m*x)
    quadgk(f, 0, π, rtol=1e-2)[1]
end


