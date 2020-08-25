using KernelAbstractions

function compute_mie_ab!(size_param, refractive_idx::Number,an,bn; N=350)
    FT = typeof(refractive_idx)

    y = size_param * refractive_idx
    # Maximum expansion (see eq. ...)
    n_max = round(Int,size_param + 4size_param^(1/3) + 2.0)
    # Make sure downward recurrence starts higher up (at least 15)
    nmx = round(Int, max(n_max, abs(y))+17 )
    @assert N>n_max
    
    # Dn as in eq 4.88, Bohren and Huffman, to calculate an and bn
    Dn = zeros(FT, nmx)
    # Downward Recursion, eq. 4.89, Bohren and Huffman
    for n = nmx-1:-1:1
        rn = n+1
        Dn[n] = (rn/y) - (1 / (Dn[n+1] + rn/y))
        #@show Dn[n]
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

# Kernels, I don't think this will run well on a GPU, might have to spell everything out
@kernel function comp_ab!(@Const(grid),an,bn,n)
    I = @index(Global, Linear)
    compute_mie_ab!(grid[I],n,view(an,:,I),view(bn,:,I))
end
