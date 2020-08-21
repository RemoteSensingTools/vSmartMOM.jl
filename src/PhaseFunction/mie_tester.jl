function compute_mie_ab(size_param::Real, refractive_idx::Number)
    y = size_param * refractive_idx
    # Maximum expansion (see eq. ...)
    n_max = round(Int,size_param + 4size_param^(1/3) + 2.0)
    # Not sure where this comes from ?
    nmx = round(Int, max(n_max, abs(y)) + 15)

    # Dn as in eq 4.88, Bohren and Huffman, to calculate an and bn
    Dn = zeros(ComplexF64, nmx)
    an = zeros(ComplexF64, n_max)
    bn = zeros(ComplexF64, n_max)
    #

    # Downward Recursion, eq. 4.89, Bohren and Huffman
    for n = nmx-1:-1:1
        rn = n+1
        Dn[n] = (rn/y) - (1 / (Dn[n+1] + rn/y))
    end

    # Get recursion for bessel functions ψ and ξ
    ψ₀ =  cos(size_param)
    ψ₁ =  sin(size_param)
    χ₀ = -sin(size_param)
    χ₁ =  cos(size_param)

    ξ₁ = ComplexF64(ψ₁, -χ₁)
    
    for n = 1:n_max
        fn = (2n+1) / (n*(n+1))
        ψ = (2n-1) * ψ₁/size_param - ψ₀
        χ = (2n-1) * χ₁/size_param - χ₀

        ξ = ComplexF64(ψ, -χ)
        t_a = Dn[n] / refractive_idx + n/size_param
        t_b = Dn[n] * refractive_idx + n/size_param
        if n==1
            @show ψ - ψ₁
            @show ξ - ξ₁
            @show t_a
            @show t_b
        end
        an[n] = (t_a * ψ - ψ₁) / (t_a * ξ - ξ₁)
        bn[n] = (t_b * ψ - ψ₁) / (t_b * ξ - ξ₁)
        
        ψ₀ = ψ₁
        ψ₁ = ψ
        χ₀ = χ₁
        χ₁ = χ
        ξ₁ = ComplexF64(ψ₁, -χ₁)
    end

    return an,bn
end