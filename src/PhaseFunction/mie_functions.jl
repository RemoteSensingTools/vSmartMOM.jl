

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
    nmx = round(Int, max(n_max, abs(y))+20 )
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
$(FUNCTIONNAME)(μ, nmax, π, τ)
Computes the associated Legendre functions  amplitude functions `π` and `τ` in Mie theory (stored internally). See eq 6 in Sanghavi 2014
- `μ` cosine of the scattering angle
- `nmax` max number of legendre terms (depends on size parameter, see [`get_n_max`](@ref))
- `π` and `τ` pre-allocated arrays, values will be stored there (need to be of size (nmax,length(μ)) )
"""
function compute_mie_π_τ!(μ, nmax, π_, τ_)
    @assert size(π_) == size(τ_) == (nmax,length(μ))
    # BH book, pages 94-96:
    π_[1,:] .= 1.0;
    π_[2,:] .= 3μ;
    τ_[1,:] .= μ;
    # This is equivalent to 3*cos(2*acos(μ))
    τ_[2,:] .= 6μ.^2 .-3;
    for n=2:nmax-1
        for i in eachindex(μ)
            π_[n+1,i] = ((2n + 1) * μ[i] * π_[n,i] - (n+1) * π_[n-1,i]) / n 
            τ_[n+1,i] = (n+1) * μ[i] * π_[n+1,i] - (n+2)*π_[n,i]
            # @show n+1,μ[i], π_[n+1,i], τ_[n+1,i], π_[n,i]
        end
    end
end

# Or use this? Not yet sure what is the best in the long run (pre-allocating? Looks a bit more cumbersome at times...)
function compute_mie_π_τ(μ, nmax)
    π_ = zeros(nmax,length(μ))
    τ_ = zeros(nmax,length(μ))
    #@assert size(π_) == size(τ_) == (nmax,length(μ))
    # BH book, pages 94-96:
    π_[1,:] .= 1.0;
    π_[2,:] .= 3μ;
    τ_[1,:] .= μ;
    # This is equivalent to 3*cos(2*acos(μ))
    τ_[2,:] .= 6μ.^2 .-3;
    for n=2:nmax-1
        for i in eachindex(μ)
            π_[n+1,i] = ((2n + 1) * μ[i] * π_[n,i] - (n+1) * π_[n-1,i]) / n 
            τ_[n+1,i] = (n+1) * μ[i] * π_[n+1,i] - (n+2)*π_[n,i]
            # @show n+1,μ[i], π_[n+1,i], τ_[n+1,i], π_[n,i]
        end
    end
    return π_,τ_
end


"""
$(FUNCTIONNAME)(an, bn, π_, τ_)
Returns the amplitude functions `S₁`,`S₂` in Mie theory
- `an` and `bn` pre-calculated Mie coefficients `an` and `bn`, see [`compute_mie_ab!`](@ref) function
- `π` and `τ` pre-calculated associated Legendre functions `π` and `τ`, see [`compute_mie_π_τ!`](@ref) function 
The function returns `S₁`,`S₂` as a function of the cosine of the scattering angle `ξ`. Users need to make sure `an` and `bn`, `π` and `τ` are pre-computed.
"""
function compute_mie_S1S2(an, bn, π_, τ_)
    nmax = size(an)[1];
    nμ   = size(π_)[2];
    S₁   = zeros(Complex{Float64}, nμ);
    S₂   = zeros(Complex{Float64}, nμ);
    for l=1:nmax
        for iμ=1:nμ 
            S₁[iμ] += (2l + 1) / (l*(l+1)) * (an[l] * τ_[l,iμ] + bn[l] * π_[l,iμ])
            S₂[iμ] += (2l + 1) / (l*(l+1)) * (an[l] * π_[l,iμ] + bn[l] *  τ_[l,iμ])
        end
    end

    return S₁,S₂
end


# DEBUG stage:
function eval_legendre(x,nmax)
    @assert nmax > 1
    #@assert size(P) == (nmax,length(x))
    P = zeros(nmax,length(x));
    # 0th Legendre polynomial, a constant
    P[1,:] .= 1;
    # 1st Legendre polynomial, x
    
    P[2,:] = x;
    for n=2:nmax-1
        for i in eachindex(x)
            P[n+1,i] = ((2n + 1) * x[i] * P[n,i] - n * P[n-1,i])/(n+1) 
        end
    end
    return P
end  

# DEBUG stage:
function average_anbn(an,bn,w,k)
    # Compute <an+bn>
    anbn = sum(w .* (an + bn),dims=2)
    # Compute <|an|²+|bn|²>
    abs2an_abs2bn = sum(w .* (abs2(an) + abs2(bn)),dims=2) 

    return anbn, abs2an_abs2bn
    #Cₑₓ = sum()
end

# DEBUG stage: Kernels, I don't think this will run well on a GPU, might have to spell everything out
@kernel function comp_ab!(@Const(grid),an,bn,Dn,n)
    I = @index(Global, Linear)
    compute_mie_ab!(grid[I],n,view(an,:,I),view(bn,:,I),Dn)
end

# DEBUG stage: Sνν
@kernel function compute_Sl_νν!(@Const(wignerA),@Const(wignerB),an,bn,w,lMax)
    # Indices over n and m
    n, m = @index(Global, NTuple)
    Sνν = zeros(eltype(an),lMax)
    # Outer loop over l
    for l = 1:lMax
        Sνν[l] += real(sum(w .* (an[n,:]' + bn[n,:]') .* (an[m,:] + bn[m,:]))) * wignerA[l,n,m]^2
        Sνν[l] += 1/2 *  sum(w .* abs2(an[n,:] + bn[n,:])) * wignerA[l,n,m]^2
    end
    
end