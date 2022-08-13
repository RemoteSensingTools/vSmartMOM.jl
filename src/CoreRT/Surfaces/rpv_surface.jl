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
    
    # Albedo normalized by π (and factor 2 for 0th Fourier Moment)
    @show rpv
    if m==0
        ρ = 2*vSmartMOM.CoreRT.reflectance(rpv, pol_type, Array(qp_μ), m)
    else
        ρ = vSmartMOM.CoreRT.reflectance(rpv, pol_type, Array(qp_μ), m)
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

#Rahman Pinty Verstraete model
function reflectance(rpv::rpvSurfaceScalar{FT}, μᵢ::FT, μᵣ::FT, dϕ::FT) where FT
    @unpack ρ₀, ρ_c, k, Θ = rpv
    # TODO: Suniti, stupid calculations here:
    θᵢ   = acos(μᵢ) #assert 0<=θᵢ<=π/2
    θᵣ   = acos(μᵣ) #assert 0<=θᵣ<=π/2
    cosg = -μᵢ*μᵣ + sin(θᵢ)*sin(θᵣ)*cos(dϕ) #RAMI form: μᵢ*μᵣ + sin(θᵢ)*sin(θᵣ)*cos(dϕ) (vSmartMOM sign convention is compatible with that of Rahman, Pinty, Verstraete, 1993) 
    G    = (tan(θᵢ)^2 + tan(θᵣ)^2 + 2*tan(θᵢ)*tan(θᵣ)*cos(dϕ))^FT(0.5) #RAMI form: (tan(θᵢ)^2 + tan(θᵣ)^2 - 2*tan(θᵢ)*tan(θᵣ)*cos(dϕ))^FT(0.5)
    return ρ₀ * rpvM(μᵢ, μᵣ, k) * rpvF(Θ, cosg) * rpvH(ρ_c, G)
end

function rpvM(μᵢ::FT, μᵣ::FT, k::FT) where FT
    return (μᵢ * μᵣ)^(k -1) /  (μᵢ + μᵣ)^(1 - k)
end

function rpvH(ρ_c::FT, G::FT) where FT
    return 1 + (1 - ρ_c) / (1 + G)
end

function rpvF(Θ::FT, cosg::FT) where FT
    θ = -θ #for RAMI only
    return (1 - Θ^2) /  (1 + Θ^2 - 2Θ * cos(π-g))^FT(1.5) #RAMI form: (1 - Θ^2) /  (1 + Θ^2 + 2Θ * cosg)^FT(1.5)
end

function reflectance(rpv::rpvSurfaceScalar{FT}, μ::Array{FT}, m::Int) where FT
    @unpack ρ₀, ρ_c, k, Θ = rpv
    f(x) = reflectance.([rpv], μ, μ', [x]) * cos(m*x)
    return 2*quadgk(f, 0, π, rtol=1e-4)[1]
end

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
    _O = (1/π)*(t-sin(t)*cos(t))*(sec(θᵢᵖ)+sec(θᵣᵖ))
    return _O - (sec(θᵢᵖ)+sec(θᵣᵖ)) + 0.5*(1+cos(ξᵖ))*sec(θᵢᵖ)*sec(θᵣᵖ)
end

function reflectance(rpv::AbstractSurfaceType, pol_type, μ::AbstractArray{FT}, m::Int) where FT
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

