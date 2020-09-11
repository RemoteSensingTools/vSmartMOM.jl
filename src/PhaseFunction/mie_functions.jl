
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
$(FUNCTIONNAME)(α, β, γ, δ, ϵ, ζ, μ; returnLeg = false)
Returns the reconstructed elements of the 4x4 scattering matrix at positions f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄ from the greek coefficients α, β, γ, δ, ϵ, ζ
- `α, β, γ, δ, ϵ, ζ` greek coefficients (arrays)
- `returnLeg` if `false` (default), just return `f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄`, if `true`, return `f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄, P, P²` (i.e. also the two legendre polynomials as matrices)
"""
function reconstruct_phase(α, β, γ, δ, ϵ, ζ, μ; returnLeg = false)
    FT = eltype(α)
    #@assert length(μ) == length(α)
    lMax = length(α);
    nμ = length(μ)
    P, P², R², T² = compute_legendre_poly(μ,lMax)
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
    f₁₁[:] = P * β                               # a₁ in Rooij notation
    f₄₄[:] = P * δ                               # a₄ in Rooij notation
    f₁₂[:] = P² * (fac .* γ)                     # b₁ in Rooij notation
    f₃₄[:] = P² * (fac .* ϵ)                     # b₂ in Rooij notation
    f₂₂[:] = R² * (fac .* α) .+ T² * (fac .* ζ)  # a₂ in Rooij notation
    f₃₃[:] = R² * (fac .* ζ) .+ T² * (fac .* α)  # a₃ in Rooij notation

    # For truncation in δ-BGE, we need P and P² as well, convenient to return here:
    if returnLeg
        return f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄, P, P²
    else
        return f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄
    end
end