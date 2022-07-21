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
    #if m == 0
        # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
        @show rpv
        if m==0
            ρ = 2*vSmartMOM.CoreRT.reflectance(rpv, pol_type.n, Array(qp_μ), m)
        else
            ρ = vSmartMOM.CoreRT.reflectance(rpv, pol_type.n, Array(qp_μ), m)
        end
        
        # Move to architecture:
        R_surf = arr_type(ρ) #arr_type(expandStokes(ρ, pol_type.n))

        
        # Source function of surface:
        if SFI
            I₀_NquadN = similar(qp_μN);
            I₀_NquadN[:] .= 0;
            I₀_NquadN[iμ₀Nstart:pol_type.n*iμ₀] = pol_type.I₀;
            
            added_layer.j₀⁺[:,1,:] .= I₀_NquadN .* exp.(-τ_sum/μ₀)';
            added_layer.j₀⁻[:,1,:] .= μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)';
        end
        R_surf = R_surf * Diagonal(qp_μN.*wt_μN)
        

        #@show size(added_layer.r⁻⁺), size(R_surf), size(added_layer.j₀⁻)
        added_layer.r⁻⁺ .= R_surf;
        added_layer.r⁺⁻ .= 0;
        added_layer.t⁺⁺ .= T_surf;
        added_layer.t⁻⁻ .= T_surf;
end

function reflectance(rpv::rpvSurfaceScalar{FT}, n, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    @unpack ρ₀, ρ_c, k, Θ = rpv
    #@unpack n = pol_type
    # TODO: Suniti, check stupid calculations here:
    θᵢ   = acos(μᵢ)
    θᵣ   = acos(μᵣ)
    cosg = μᵢ*μᵣ + sin(θᵢ)*sin(θᵣ)*cos(dϕ)
    G    = (tan(θᵢ)^2 + tan(θᵣ)^2 - 2*tan(θᵢ)*tan(θᵣ)*cos(dϕ))^FT(0.5)
    BRF = ρ₀ * rpvM(μᵢ, μᵣ, k) * rpvF(Θ, cosg) * rpvH(ρ_c, G)
    #vcat(BRF, zeros(eltype(BRF),n-1))
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

function reflectance(rpv::AbstractSurfaceType, n::Int, μ::AbstractArray{FT}, m::Int) where FT
    f(x) = reflectance.((rpv,),n, μ, μ', x) * cos(m*x)
    quadgk(f, 0, 2π, rtol=1e-6)[1]
end

function expandStokes(v::AbstractArray)
    d1,d2 = size(v)
    n = size(v[1,1])
    dn = n - 1
    M = Matrix{eltype(eltype(v))}(undef, d1*n, d1*n)
    for i in 1:d1, j in 1:d2
        ii = (i-1)*n + 1
        jj = (j-1)*n + 1 
        M[ii:ii+dn, jj:jj+dn] .= Diagonal(v[i,j])
    end
    return M
end

function expandStokes2(v)
    return Diagonal(v)
end
