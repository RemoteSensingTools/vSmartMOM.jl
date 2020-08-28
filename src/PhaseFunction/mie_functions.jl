using KernelAbstractions


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
            
        ψ₀ = ψ₁
        ψ₁ = ψ
        χ₀ = χ₁
        χ₁ = χ
        ξ₁ = FT(ψ₁, -χ₁)
    end
    return nothing
end

function get_n_max(size_parameter)
    FT = eltype(size_parameter)
    round(Int,size_parameter + FT(4.05)*size_parameter^(1/3) + FT(10))
end

function average_anbn(an,bn,w,k)
    # Compute <an+bn>
    anbn = sum(w .* (an + bn),dims=2)
    # Compute <|an|²+|bn|²>
    abs2an_abs2bn = sum(w .* (abs2(an) + abs2(bn)),dims=2) 

    return anbn, abs2an_abs2bn
    #Cₑₓ = sum()
end

# Kernels, I don't think this will run well on a GPU, might have to spell everything out
@kernel function comp_ab!(@Const(grid),an,bn,Dn,n)
    I = @index(Global, Linear)
    compute_mie_ab!(grid[I],n,view(an,:,I),view(bn,:,I),Dn)
end

# Sνν
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

function compute_pi_tau!(μ, nmax, π_, τ_)
    #@assert length(π_) == length(τ_) == nmax
    π_[1,:] .= 1.0;
    π_[2,:] .= 3μ;
    τ_[1,:] .= μ;
    τ_[2,:] .= 6μ.^2 .-3;
    for n=2:nmax-1
        for i in eachindex(μ)
            π_[n+1,i] = ((2n + 1) * μ[i] * π_[n,i] - (n+1) * π_[n-1,1]) / n 
            τ_[n+1,i] = (n+1) * μ[i] * π_[n+1,i] - (n+2)*π_[n,i]
            # @show n+1,μ[i], π_[n+1,i], τ_[n+1,i], π_[n,i]
        end
    end
end

function eval_legendre!(x,nmax,P)
    @assert nmax > 1
    @assert size(P) == (nmax,length(x))
    # 0th Legendre polynomial, a constant
    P[1,:] .= 1;
    # 1st Legendre polynomial, x
    P[2,:] = x;
    for n=2:nmax-1
        for i in eachindex(x)
            P[n+1,:] = ((2n + 1) * x[i] * P[n,i] - n * P[n-1,1])/(n+1) 
        end
    end
end  
