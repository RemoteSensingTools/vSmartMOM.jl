
"""
$(FUNCTIONNAME)(size_parameter,refractive_idx::Number,an,bn,Dn)
Computes Mie coefficients `an` and `bn` as a function of size parameter and complex refractive index. See eq 4.88 in Bohren and Huffman
- `size_parameter` size parameter of the aerosol (2πr/λ)
- `refractive_idx` refractive index of the aerosol (complex number)
- `an` and `bn` pre-allocated arrays, need to match at least n_max for the given size parameter
- `Dn` pre-allocated array for the logarithmic derivative (see BH, eq 4.88) (need to check whether it can be created internally without causing too many allocations)

The function returns a rounded integer, following conventions by BH, Rooj/Stap, Siewert 
"""
function compute_mie_ab!(size_param, refractive_idx::Number,an,bn,Dn)
    FT = typeof(refractive_idx)

    y = size_param * refractive_idx
    # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
    n_max = get_n_max(size_param)

    # Make sure downward recurrence starts higher up (at least 15, check eq. A9 in de Rooij and Stap, 1984, may need to check what is needed)
    nmx = round(Int, max(n_max, abs(y))+50 )
    @assert size(an)[1]>=n_max
    @assert size(an) == size(bn)
    fill!(Dn,0);
    # Dn as in eq 4.88, Bohren and Huffman, to calculate an and bn
    #Dn = zeros(FT, nmx)
    # Downward Recursion, eq. 4.89, Bohren and Huffman
    for n = nmx-1:-1:1
        rn = n+1
        #@show n, (rn/y) - (1 / (Dn[n+1] + rn/y))
        Dn[n] = (rn/y) - (1 / (Dn[n+1] + rn/y))
        #@show n, Dn[n]
    end

    # Get recursion for bessel functions ψ and ξ
    ψ₀ =  cos(size_param)
    ψ₁ =  sin(size_param)
    χ₀ = -sin(size_param)
    χ₁ =  cos(size_param)

    ξ₁ = FT(ψ₁, -χ₁)

    # This solves Bohren and Huffman eq. 4.88 for an and bn, computing updated ψ and ξ on the fly
    for n = 1:n_max  
        fn = (2n+1) / (n*(n+1))
        ψ  = (2n-1) * ψ₁/size_param - ψ₀
        χ  = (2n-1) * χ₁/size_param - χ₀

        ξ   = FT(ψ, -χ)
        t_a = Dn[n] / refractive_idx + n/size_param
        t_b = Dn[n] * refractive_idx + n/size_param
         
        an[n] = (t_a * ψ - ψ₁) / (t_a * ξ - ξ₁)
        bn[n] = (t_b * ψ - ψ₁) / (t_b * ξ - ξ₁)
        #@show n, ψ, ψ₁, ξ,  ξ₁, real(an[n])
        ψ₀ = ψ₁
        ψ₁ = ψ
        χ₀ = χ₁
        χ₁ = χ
        ξ₁ = FT(ψ₁, -χ₁)
    end
    return nothing
end


"""
$(FUNCTIONNAME)(size_parameter)
Computes the number of required Legendre functions  for a given size parameter. See eq 6 in Sanghavi 2014
- `size_parameter` size parameter of the aerosol (2πr/λ)
The function returns a rounded integer, following conventions by BH, Rooj/Stap, Siewert 
"""
function get_n_max(size_parameter)
    FT = eltype(size_parameter)
    round(Int,size_parameter + FT(4.05)*size_parameter^(1/3) + FT(10))
end


"""
$(FUNCTIONNAME)(an, bn, π_, τ_, S₁, S₂)
Determines the amplitude functions `S₁`,`S₂` in Mie theory
- `an` and `bn` pre-calculated Mie coefficients `an` and `bn`, see [`compute_mie_ab!`](@ref) function
- `π` and `τ` pre-calculated associated Legendre functions `π` and `τ`, see [`compute_mie_π_τ!`](@ref) function 
The function returns `S₁`,`S₂` as a function of the cosine of the scattering angle `ξ`. Users need to make sure `an` and `bn`, `π` and `τ` are pre-computed.
"""
function compute_mie_S₁S₂!(an, bn, π_, τ_, S₁, S₂)
    FT = eltype(an)
    nmax = size(an,1);
    nμ   = size(π_,1);
    @assert size(S₁) == size(S₂)
    @assert length(S₁) == nμ

    for l=1:nmax
        for iμ=1:nμ 
            S₁[iμ] += (2l + 1) / (l*(l+1)) * (an[l] * τ_[iμ,l] + bn[l] * π_[iμ,l])
            S₂[iμ] += (2l + 1) / (l*(l+1)) * (an[l] * π_[iμ,l] + bn[l] * τ_[iμ,l])
        end
    end
    return nothing
end


"""
$(FUNCTIONNAME)(n,xmin,xmax; norm=false)
Returns the `n` Gauss-Legendre quadrature points and weights with a change of interval between xmin and xmax
- `n` number of quadrature points
- `xmin`,`xmax` lower and upper bound of integral
- `norm`: if `true`, normalizes the weights so that a mean can be computed instead of full integration
The function returns `n` quadrature points ξ within [xmin,xmax] with associated weightes `w` 
"""
function gauleg(n,xmin,xmax; norm=false)
    ξ,w = gausslegendre( n )
    ξ = (xmax-xmin)/2 * ξ .+ (xmin+xmax)/2
    if norm
        w /= sum(w)
    else
        w *= (xmax-xmin)/2
    end
    return ξ,w
end

"""
$(FUNCTIONNAME)(greek_coefs, μ; returnLeg = false)
Returns the reconstructed elements of the 4x4 scattering matrix at positions f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ from the greek coefficients
- `greek_coefs` greek coefficients (Domke Type)
- `returnLeg` if `false` (default), just return `f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄`, if `true`, return `f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄, P, P²` (i.e. also the two legendre polynomials as matrices)
"""
function reconstruct_phase(greek_coefs, μ; returnLeg = false)
    FT = eltype(greek_coefs.α)
    #@assert length(μ) == length(α)
    lMax = length(greek_coefs.α);
    nμ = length(μ)
    P, P², R², T² = compute_legendre_poly(μ, lMax)
    # To stay general, we also don't assume f₂₂=f₁₁ or f₄₄=f₃₃
    # which only holds for spherical
    f₁₁   = zeros(FT, nμ)
    f₃₃   = zeros(FT, nμ)
    f₁₂   = zeros(FT, nμ)
    f₃₄   = zeros(FT, nμ)
    f₂₂   = zeros(FT, nμ)
    f₄₄   = zeros(FT, nμ)

    fac = zeros(lMax);
    for l=2:lMax-1
        fac[l+1] = sqrt(1 / ( ( l-1) * l * (l+1) * (l+2) ));
    end
    # In matrix form:
    f₁₁[:] = P * greek_coefs.β                               # a₁ in Rooij notation
    f₄₄[:] = P * greek_coefs.δ                               # a₄ in Rooij notation
    f₁₂[:] = P² * (fac .* greek_coefs.γ)                     # b₁ in Rooij notation
    f₃₄[:] = P² * (fac .* greek_coefs.ϵ)                     # b₂ in Rooij notation
    f₂₂[:] = R² * (fac .* greek_coefs.α) .+ T² * (fac .* greek_coefs.ζ)  # a₂ in Rooij notation
    f₃₃[:] = R² * (fac .* greek_coefs.ζ) .+ T² * (fac .* greek_coefs.α)  # a₃ in Rooij notation

    # For truncation in δ-BGE, we need P and P² as well, convenient to return here:
    if returnLeg
        return f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄, P, P²
    else
        return f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄
    end
end

function get_greek_rayleigh(depol)
    # Rayleigh Greek Parameters
    dpl_p = (1 - depol)  / (1 + depol/2)
    dpl_q = (1 + depol)  / (1 - depol)
    dpl_r = (1 - 2depol) / (1 - depol)
  
    α  =  [0.0, 0.0,             3dpl_p]
    β  =  [1.0, 0.0,             0.5*dpl_p]
    γ  =  [0.0, 0.0,             dpl_p*sqrt(1.5)] 
    δ  =  [0.0, dpl_p*dpl_r*1.5, 0.0] 
    ϵ  =  [0.0, 0.0,             0.0] 
    ζ  =  [0.0, 0.0,             0.0]
    return α, β, γ, δ, ϵ, ζ 
end

function construct_Π_matrix(mo::Stokes_IQUV, P,R,T,l::Int,m::Int; sign_change=false)
    if sign_change # (basically gets it for -μ due to symmetries on P,R,T)
        if mod(l-m,2) == 1
            Π = [SMatrix{4,4}([-P[i,l,m] 0 0 0 ; 0 -R[i,l,m] -T[i,l,m] 0; 0 -T[i,l,m] -R[i,l,m] 0; 0 0 0 -P[i,l,m]]) for i in 1:size(P,1)] 
        else
            Π = [SMatrix{4,4}([P[i,l,m] 0 0 0 ; 0 R[i,l,m] T[i,l,m] 0; 0 T[i,l,m] R[i,l,m] 0; 0 0 0 P[i,l,m]]) for i in 1:size(P,1)]
        end
    else
        Π = [SMatrix{4,4}([P[i,l,m] 0 0 0 ; 0 R[i,l,m] -T[i,l,m] 0; 0 -T[i,l,m] R[i,l,m] 0; 0 0 0 P[i,l,m]]) for i in 1:size(P,1)]
    end
    return Π
end

function construct_Π_matrix(mo::Stokes_IQU, P,R,T,l::Int,m::Int; sign_change=false)
    if sign_change # (basically gets it for -μ due to symmetries on P,R,T)
        if mod(l-m,2) == 1
            Π = [SMatrix{3,3}([-P[i,l,m] 0 0  ; 0 -R[i,l,m] -T[i,l,m] ; 0 -T[i,l,m] -R[i,l,m] ]) for i in 1:size(P,1)] 
        else
            Π = [SMatrix{3,3}([P[i,l,m] 0 0  ; 0 R[i,l,m] T[i,l,m] ; 0 T[i,l,m] R[i,l,m] ]) for i in 1:size(P,1)]
        end
    else
        Π = [SMatrix{3,3}([P[i,l,m] 0 0  ; 0 R[i,l,m] -T[i,l,m] ; 0 -T[i,l,m] R[i,l,m] ]) for i in 1:size(P,1)]
    end
    return Π
end

function construct_Π_matrix(mod::Stokes_I, P,R,T,l::Int,m::Int; sign_change=false)
    if sign_change # (basically gets it for -μ due to symmetries on P,R,T)
        Π = -P[:,l,m]
    else
        Π = P[:,l,m]
    end        
end

function construct_B_matrix(mod::Stokes_IQUV, α, β, γ, δ, ϵ, ζ,l::Int)
    𝐁 = SMatrix{4,4}([β[l] γ[l] 0 0 ; γ[l] α[l] 0 0; 0 0 ζ[l] ϵ[l]; 0 0 -ϵ[l] δ[l]])
end

function construct_B_matrix(mod::Stokes_IQU, α, β, γ, δ, ϵ, ζ,l::Int)
    𝐁 = SMatrix{3,3}([β[l] γ[l] 0 ; γ[l] α[l] 0 ; 0 0 ζ[l]])
end

function construct_B_matrix(mod::Stokes_I, α, β, γ, δ, ϵ, ζ,l::Int)
    𝐁 = β[l]
end


function compute_Z_moments(mod::AbstractPolarizationType, μ, α, β, γ, δ, ϵ, ζ, m::Int)
    FT = eltype(β)
    n = length(μ)
    
    # Set prefactor for moments (note 1-notation for `m` here):
    if m==1
        fact=0.5
    else
        fact = 1.0
    end

    # get Lmax just from length of array:
    Lmax = length(β)
    # Check that all μ are positive here ([0,1])
    @assert all(0 .< μ .≤ 1) "all μ's within compute_Z_moments have to be ∈ ]0,1]"
    # Compute legendre Polynomials at μ and up to lmax
    P,R,T    = PhaseFunction.compute_associated_legendre_PRT(μ,Lmax)
    P⁻,R⁻,T⁻ = PhaseFunction.compute_associated_legendre_PRT(-μ,Lmax)
  
    # Pre-compute all required B matrices
    𝐁_all = [construct_B_matrix(mod,α, β, γ, δ, ϵ, ζ,i) for i in 1:Lmax]
    # Get dimension of square matrix (easier for Scalar/Stokes dimensions)
    B_dim = Int(sqrt(length(𝐁_all[1])))
    
    # Create matrices:
    nb = B_dim*n
    𝐙⁺⁺ = zeros(FT,nb,nb)
    𝐙⁺⁻ = zeros(FT,nb,nb)
    A⁺⁺ = zeros(FT,B_dim,B_dim,n,n)
    A⁺⁻ = zeros(FT,B_dim,B_dim,n,n)

    # Iterate over l
    for l = m:Lmax
        #@show l
        # B matrix for l
        𝐁 = 𝐁_all[l];
        # Construct Π matrix for l,m pair (change to in place later!)
        # See eq. 15 in Sanghavi 2014, note that P,R,T are already normalized
        Π  = construct_Π_matrix(mod,P,R,T,l,m)
        #Π⁻ = construct_Π_matrix(mod,P,R,T,l,m; sign_change=true)
        Π⁻ = construct_Π_matrix(mod,P⁻,R⁻,T⁻,l,m)
        # Iterate over angles
        for i in eachindex(μ), j in eachindex(μ)
            if B_dim==1
                A⁺⁺[B_dim,B_dim,i,j] += Π[i] * 𝐁 * Π[j]
                A⁺⁻[B_dim,B_dim,i,j] += Π[i] * 𝐁 * Π⁻[j]
            else
                A⁺⁺[:,:,i,j] += Π[i] * 𝐁 * Π[j]
                A⁺⁻[:,:,i,j] += Π[i] * 𝐁 * Π⁻[j]
            end
        end
    end
    # Now get to the Z part:
    for imu in eachindex(μ), jmu in eachindex(μ)
        # Indices adjusted for size of A
        ii=(imu-1)*B_dim
        jj=(jmu-1)*B_dim
        
        # This is equivalent to Z̄ = 1/(1+δ) * C̄m+S̄m = 1/(1+δ) * (A+DAD+AD-DA) (see eq 11 in Sanghavi et al, 2013)
        for i=1:B_dim, j=1:B_dim
            𝐙⁺⁺[ii+i,jj+j] = 2fact*A⁺⁺[i,j,imu,jmu]
            if i<=2 && j>=3
                𝐙⁺⁻[ii+i,jj+j] = -2fact*A⁺⁻[i,j,imu,jmu]
            elseif i>=3 && j<=2
                𝐙⁺⁻[ii+i,jj+j] = -2fact*A⁺⁻[i,j,imu,jmu]
            else
                𝐙⁺⁻[ii+i,jj+j] = 2fact*A⁺⁻[i,j,imu,jmu]
            end
        end
    end
    return 𝐙⁺⁺,𝐙⁺⁻
end
    
