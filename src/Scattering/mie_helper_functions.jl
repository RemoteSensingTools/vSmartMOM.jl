#=

This file contains helper functions that are used throughout the module

=#

"""
    exp_m1(x)

Convenience function returning ``(-1)^x`` using the parity of `x` (1 if `x` is even, -1 if odd).
Used in Mie/PCW summations for sign alternation. Name reflects common notation in scattering literature.
"""
exp_m1(x) = iseven(x) ? 1 : -1

"""
    $(FUNCTIONNAME)(size_parameter)
Computes the number of required Legendre functions  for a given size parameter. 
See eq 6 in Sanghavi 2014
- `size_parameter` size parameter of the aerosol (2πr/λ)
The function returns a rounded integer, following conventions by BH, Rooj/Stap, Siewert 
"""
get_n_max(size_parameter) = (size_parameter>8.0) ? round(Int, size_parameter + 4.05 * size_parameter^(1/3) + 10) : round(Int, size_parameter + 4.0 * size_parameter^(1/3) + 1)

@doc raw"""
    $(FUNCTIONNAME)(size_param, refractive_idx::Number, an, bn, Dn)

Compute Mie coefficients ``a_n`` and ``b_n`` from the size parameter and complex
refractive index (Bohren & Huffman, eq. 4.88).

The logarithmic derivative ``D_n(y)`` is obtained via downward recursion
(BH eq. 4.89; de Rooij & Stap 1984, eq. A9).  This recursion is **always
performed in Float64** regardless of the element type of `Dn`, because the
cancellation-prone continued-fraction form loses significant digits in Float32
for ``|y| \gtrsim 50``.

# Arguments
- `size_param`: size parameter ``x = 2\pi r/\lambda``
- `refractive_idx`: complex refractive index of the particle
- `an`, `bn`: pre-allocated arrays (length ≥ `n_max(x)`) for the Mie coefficients
- `Dn`: pre-allocated array for the logarithmic derivative (length chosen by caller)
"""
function compute_mie_ab!(size_param, refractive_idx::Number, an, bn, Dn)
    # Compute y
    y = size_param * refractive_idx

    # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
    n_max = get_n_max(size_param)

    # Make sure downward recurrence starts higher up
    # (at least 15, check eq. A9 in de Rooij and Stap, 1984, may need to check what is needed)
    nmx = length(Dn)
    @assert size(an)[1] >= n_max
    @assert size(an) == size(bn)
    fill!(Dn, 0);

    # Dn as in eq 4.88, Bohren and Huffman, to calculate an and bn
    # Downward Recursion, eq. 4.89, Bohren and Huffman
    # Always performed in Float64 for numerical stability
    y64 = Complex{Float64}(y)
    Dn_prev = Complex{Float64}(0)
    @inbounds for n = (nmx - 1):-1:1
        ratio = (n + 1) / y64
        Dn_prev = ratio - 1 / (Dn_prev + ratio)
        Dn[n] = Dn_prev
    end

    # Get recursion for bessel functions ψ and ξ
    ψ₀, ψ₁, χ₀, χ₁ =  (cos(size_param), sin(size_param), -sin(size_param), cos(size_param))
    ξ₁ = ψ₁ + χ₁*im

    # This solves Bohren and Huffman eq. 4.88 for an and bn, computing updated ψ and ξ on the fly
    @inbounds for n = 1:n_max  
        ψ  = (2n - 1) * ψ₁ / size_param - ψ₀
        χ  = (2n - 1) * χ₁ / size_param - χ₀
        ξ   = ψ + χ*im
        t_a = Dn[n] / refractive_idx + n / size_param
        t_b = Dn[n] * refractive_idx + n / size_param
         
        an[n] = (t_a * ψ - ψ₁) / (t_a * ξ - ξ₁)
        bn[n] = (t_b * ψ - ψ₁) / (t_b * ξ - ξ₁)

        ψ₀ = ψ₁
        ψ₁ = ψ
        χ₀ = χ₁
        χ₁ = χ
        ξ₁ = ψ₁ + χ₁*im
    end
end

"""
    $(FUNCTIONNAME)(size_param, refractive_idx, an, bn, Dn)

Alternative implementation of Mie coefficients ``a_n``, ``b_n`` using upward recursion for ψ and χ
and Lentz's continued-fraction method for the logarithmic derivative. See Bohren & Huffman eq. 4.88.
Modifies `an`, `bn`, and `Dn` in-place. Same interface as [`compute_mie_ab!`](@ref).
"""
function compute_mie_ab_new!(size_param, refractive_idx::Number, an, bn, Dn)
    # Compute y
    y = size_param * refractive_idx

    # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
    n_max = get_n_max(size_param)

    # Make sure downward recurrence starts higher up 
    # (at least 15, check eq. A9 in de Rooij and Stap, 1984, may need to check what is needed)
    nmx = length(Dn)
    @assert size(an)[1] >= n_max
    @assert size(an) == size(bn)
    fill!(Dn, 0);


    
    #Computing ψ using downward recursion
    N_ = n_max+60
    
    ψ = zeros(N_)
    ψ[end]   = 0.0
    ψ[end-1] = 1.0
    for idx=N_-2:-1:1
        ψ[idx] = (2idx+1)*ψ[idx+1]/size_param - ψ[idx+2];    
    end
    
    
    #Computing ψ using upward recursion
    N_ = n_max
    ψ = zeros(N_)
    ψ[1] = sin(size_param);
    if N_>1
        ψ[2]  = (sin(size_param)/size_param)-cos(size_param);
        for idx = 3:N_
            ψ[idx] = (2idx-3)*ψ[idx-1]/size_param - ψ[idx-2];
        end
    end
    #computing χ using upward recursion
    N_ = n_max
    χ = zeros(N_)     
    χ[1] = cos(size_param)
    if N_>1
        χ[2] = cos(size_param)/size_param + sin(size_param)            
        for idx=3:N_
            χ[idx] = (2idx-3)*χ[idx-1]/size_param - χ[idx-2];
        end
    end

    #computing An (Lentz)
    result = zeros(Complex{FT}, nmx);
    z=y
    for n=1:nmx
        zinv   = 2/z
        alpha_ = (n + 0.5)*zinv
        aj     =-(n + 1.5)*zinv
        alpha_j1 = aj+1/alpha_
        alpha_j2 = aj
      
        ratio = alpha_j1/alpha_j2
        runratio = alpha_*ratio
      
        while abs(abs(ratio)-1) > 1e-20
            aj=zinv-aj
            alpha_j1=1/alpha_j1+aj
            alpha_j2=1/alpha_j2+aj
            ratio=alpha_j1/alpha_j2
            
            epsilon1 = 1.0e-2
            compare_1 = abs(alpha_j1/aj)
            compare_2 = abs(alpha_j2/aj)

            if abs(compare_1)<=epsilon1 || abs(compare_2)<=epsilon1   
                zinv *= -1
                aj = zinv - aj
                ratio = (1+aj*alpha_j1)/(1+aj*alpha_j2)
                alpha_j1 = aj + 1/alpha_j1
                alpha_j2 = aj + 1/alpha_j2
            end
            zinv *= -1;
            runratio=ratio*runratio;
        end 
        result[n] = -n/z;
        result[n] += runratio;
    end

    # Dn as in eq 4.88, Bohren and Huffman, to calculate an and bn
    # Downward Recursion, eq. 4.89, Bohren and Huffman
    @inbounds for n = (nmx - 1):-1:1
        Dn[n] = ((n+1) / y) - (1 / (Dn[n+1] + (n+1) / y))
    end

    # Get recursion for bessel functions ψ and ξ
    ψ₀, ψ₁, χ₀, χ₁ =  (cos(size_param), sin(size_param), -sin(size_param), cos(size_param))
    ξ₁ = ψ₁ -χ₁*im

    ψ0 = zeros(N_)
    χ0 = zeros(N_)

    # This solves Bohren and Huffman eq. 4.88 for an and bn, computing updated ψ and ξ on the fly
    @inbounds for n = 1:n_max  
        ψ  = (2n - 1) * ψ₁ / size_param - ψ₀
        χ  = (2n - 1) * χ₁ / size_param - χ₀
        ξ   = ψ -χ*im
        t_a = Dn[n] / refractive_idx + n / size_param
        t_b = Dn[n] * refractive_idx + n / size_param
        ψ0[n] = ψ₁
        χ0[n] = χ₁
        an[n] = (t_a * ψ - ψ₁) / (t_a * ξ - ξ₁)
        bn[n] = (t_b * ψ - ψ₁) / (t_b * ξ - ξ₁)

        ψ₀ = ψ₁
        ψ₁ = ψ
        χ₀ = χ₁
        χ₁ = χ
        ξ₁ = ψ₁ -χ₁*im
    end
end


""" 
    $(FUNCTIONNAME)(model::MieModel, λ, radius)
Compute all an, bn using compute_mie_ab!
Input: MieModel, wavelength (λ), radius
Output: an, bn. Both of shape (aerosol.nquad_radius, N_max) (N_max from aerosol.r_max)
"""
function compute_anbn(model::MieModel, λ, radius)
    
    (; computation_type, aerosol, r_max, nquad_radius, λ, polarization_type, truncation_type, wigner_A, wigner_B) = model
    (; size_distribution, nᵣ, nᵢ) = aerosol

    FT = eltype(λ)
    FT2 = eltype(nᵣ)

    # Find overall N_max from the maximum radius
    N_max = Scattering.get_n_max(2 * π * r_max / λ)

    # Where to store an, bn, computed over size distribution
    an = zeros(Complex{FT2}, nquad_radius, N_max)
    bn = zeros(Complex{FT2}, nquad_radius, N_max)

    # Pre-allocate Dn buffer for max size parameter (from r_max)
    y_max = (2 * π * r_max / λ) * abs(nᵣ - nᵢ)
    nmx_max = round(Int, max(N_max, y_max) + 51)
    Dn = zeros(Complex{FT2}, nmx_max)

    # Loop over the size distribution, and compute an, bn, for each size
    for i in 1:nquad_radius

        # Get current radius and size parameter
        r = radius[i] 
        size_param = 2 * π * r / λ

        # Compute an, bn (Dn is zeroed inside compute_mie_ab!)
        Scattering.compute_mie_ab!(size_param, nᵣ - nᵢ * im,
                                      view(an, i, :),
                                      view(bn, i, :), Dn)
    end

    return an, bn;
end

"""
    $(FUNCTIONNAME)(an, bn, ab_pairs, w, Nmax, N_max_)

Precompute size-averaged Mie coefficient products for PCW Greek-coefficient evaluation.

Fills `ab_pairs` with ``\\langle a_n^* a_m \\rangle``, ``\\langle a_n^* b_m \\rangle``, etc.,
enabling fast evaluation of ``(a_n^* \\pm b_n^*)(a_m \\pm b_m)`` in Sanghavi (2014) Eq. 22.

# Arguments
- `an`, `bn`: Mie coefficients, shape `(nquad_radius, N_max)`
- `ab_pairs`: tuple of 4 matrices `(mat_anam, mat_anbm, mat_bnam, mat_bnbm)` to fill
- `w`: probability weights for each radius
- `Nmax`: maximum expansion order
- `N_max_`: per-radius truncation indices (length `nquad_radius`)
"""
function compute_avg_anbns!(an, bn, ab_pairs, w, Nmax, N_max_)
    FT2 = eltype(an)

    # Unpack ab_pairs
    mat_anam, mat_anbm, mat_bnam, mat_bnbm = ab_pairs

    # Fill all matrices with 0
    [fill!(mat, 0) for mat in [mat_anam, mat_bnbm, mat_anbm, mat_bnam]]

    @inbounds for n in 1:Nmax, m in n:Nmax
        anam = bnbm = anbm = bnam = FT2(0);
        @inbounds for i = 1:size(an, 1)
            if m < N_max_[i] && n < N_max_[i]
                anam += w[i] * an[i,n]' * an[i,m]
                bnbm += w[i] * bn[i,n]' * bn[i,m]
                anbm += w[i] * an[i,n]' * bn[i,m]
                bnam += w[i] * bn[i,n]' * an[i,m]
            end
        end 
        @inbounds mat_anam[m,n] = anam;
        @inbounds mat_bnbm[m,n] = bnbm;
        @inbounds mat_anbm[m,n] = anbm;
        @inbounds mat_bnam[m,n] = bnam;
    end
end

"""
    $(FUNCTIONNAME)(an, bn, π_, τ_, S₁, S₂)
Determines the amplitude functions `S₁`,`S₂` in Mie theory
- `an` and `bn` pre-calculated Mie coefficients `an` and `bn`, see [`compute_mie_ab!`](@ref) function
- `π` and `τ` pre-calculated associated Legendre functions `π` and `τ`, see [`compute_mie_π_τ`](@ref) function
The function returns `S₁`,`S₂` as a function of the cosine of the scattering angle `ξ`. 
Users need to make sure `an` and `bn`, `π` and `τ` are pre-computed.
"""
function compute_mie_S₁S₂!(an, bn, π_, τ_, S₁, S₂)
    FT = eltype(an)
    nmax = size(an, 1);
    nμ   = size(π_, 1);

    # Verify sizes
    @assert size(S₁) == size(S₂)
    @assert length(S₁) == nμ

    @inbounds for l in 1:nmax, iμ in 1:nμ 
            S₁[iμ] += (2l + 1) / (l * (l + 1)) * (an[l] * τ_[iμ,l] + bn[l] * π_[iμ,l])
            S₂[iμ] += (2l + 1) / (l * (l + 1)) * (an[l] * π_[iμ,l] + bn[l] * τ_[iμ,l])
    end
end

"""
    $(FUNCTIONNAME)(n,xmin,xmax; norm=false)
Returns the `n` Gauss-Legendre quadrature points and weights with a change of interval between xmin and xmax
- `n` number of quadrature points
- `xmin`,`xmax` lower and upper bound of integral
- `norm`: if `true`, normalizes the weights so that a mean can be computed instead of full integration
The function returns `n` quadrature points ξ within [xmin,xmax] with associated weightes `w` 
"""
function gauleg(n, xmin, xmax; norm=false)
    ξ, w = gausslegendre(n)
    ξ = (xmax - xmin) / 2 * ξ .+ (xmin + xmax) / 2
    norm ? w /= sum(w) : w *= (xmax - xmin) / 2
    return ξ, w
end

@doc raw"""
    gauleg_log(n, r_min, r_max; norm=false)

Gauss-Legendre quadrature with nodes equidistant in ``\ln r``.

For a log-normal size distribution

```math
n(r) = \frac{1}{r\,\ln\sigma\,\sqrt{2\pi}}
       \exp\!\Bigl[-\frac{(\ln r - \ln r_m)^2}{2\ln^2\sigma}\Bigr],
```

the substitution ``u = \ln r`` gives ``dr = r\,du`` and transforms the
integrand into a Gaussian in ``u``.  Placing Gauss-Legendre nodes in
``[\ln r_{\min},\,\ln r_{\max}]`` concentrates points where the
distribution has significant mass.

The returned weights include the Jacobian ``r = e^u``:

```math
w_i^{\text{log}} = w_i^{\text{GL}}\;\cdot r_i\;\cdot
                    \frac{\ln r_{\max} - \ln r_{\min}}{2}.
```

# Arguments
- `n`: number of quadrature points
- `r_min`, `r_max`: integration bounds in radius space (both > 0)
- `norm`: if `true`, normalize weights to sum to 1 (for computing means)

# Returns
`(r, w)` — radius nodes and corresponding weights
"""
function gauleg_log(n, r_min, r_max; norm=false)
    ξ, w = gausslegendre(n)
    ln_min, ln_max = log(r_min), log(r_max)
    ln_r = (ln_max - ln_min) / 2 * ξ .+ (ln_max + ln_min) / 2
    r = exp.(ln_r)
    w = w .* r .* (ln_max - ln_min) / 2
    if norm
        w ./= sum(w)
    end
    return r, w
end

@doc raw"""
    reconstruct_phase(greek_coefs, μ; returnLeg = false)

Reconstruct angle-space scattering-matrix elements from Greek coefficients.

Reference: Sanghavi (2014), Fourier/Greek framework around Eq. (17).

`f₁₁` represents the scalar phase function and is normalized as:
```math
\frac{1}{4\pi}\int_0^{2\pi}d\phi \int_{-1}^1 p(\mu) d\mu  = 1
```

Using Legendre basis matrices computed at `μ`, this function evaluates:

```math
f_{11} = P\beta,\quad f_{44}=P\delta,\quad
f_{12}=P^2(\mathrm{fac}\odot\gamma),\quad
f_{34}=P^2(\mathrm{fac}\odot\epsilon),
```

```math
f_{22}=R^2(\mathrm{fac}\odot\alpha)+T^2(\mathrm{fac}\odot\zeta),\quad
f_{33}=R^2(\mathrm{fac}\odot\zeta)+T^2(\mathrm{fac}\odot\alpha).
```

# Arguments
- `greek_coefs`: [`GreekCoefs`](@ref) coefficients.
- `μ`: cosine of scattering angles where the matrix should be reconstructed.
- `returnLeg`: if `true`, also return `(P, P²)`.

# Returns
- `ScatteringMatrix` when `returnLeg=false`.
- `(ScatteringMatrix, P, P²)` when `returnLeg=true`.
"""
function reconstruct_phase(greek_coefs, μ; returnLeg=false)

    FT = eltype(greek_coefs.α)
    l_max = length(greek_coefs.α);
    nμ = length(μ)

    # Compute legendre polynomials
    P, P², R², T² = compute_legendre_poly(μ, l_max)

    # To stay general, we also don't assume f₂₂=f₁₁ or f₄₄=f₃₃
    # which only holds for spherical
    f₁₁, f₃₃, f₁₂, f₃₄, f₂₂, f₄₄ = (zeros(FT, nμ), zeros(FT, nμ), zeros(FT, nμ), 
                                    zeros(FT, nμ), zeros(FT, nμ), zeros(FT, nμ))

    # Compute prefactor
    fac = zeros(l_max);
    [fac[l + 1] = sqrt(1 / ((l-1) * l * (l+1) * (l+2))) for l = 2:(l_max-1)]

    # In matrix form:
    f₁₁[:] = P * greek_coefs.β                                           # a₁ in Rooij notation
    f₄₄[:] = P * greek_coefs.δ                                           # a₄ in Rooij notation
    f₁₂[:] = P² * (fac .* greek_coefs.γ)                                 # b₁ in Rooij notation
    f₃₄[:] = P² * (fac .* greek_coefs.ϵ)                                 # b₂ in Rooij notation
    f₂₂[:] = R² * (fac .* greek_coefs.α) .+ T² * (fac .* greek_coefs.ζ)  # a₂ in Rooij notation
    f₃₃[:] = R² * (fac .* greek_coefs.ζ) .+ T² * (fac .* greek_coefs.α)  # a₃ in Rooij notation

    # Put elements into a struct
    scattering_matrix = ScatteringMatrix(f₁₁, f₁₂, f₂₂, f₃₃, f₃₄, f₄₄)

    # For truncation in δ-BGE, we need P and P² as well, convenient to return here:
    return returnLeg ? (scattering_matrix, P, P²) : scattering_matrix
end

"""
    $(FUNCTIONNAME)(depol)
Returns the greek coefficients (as [`GreekCoefs`](@ref)) of the Rayleigh phase function given 
depolarization value. 
- `depol` Depolarization (best use 0 as default )
"""
function get_greek_rayleigh(depol::Number)
    FT = eltype(depol)
    # Rayleigh Greek Parameters
    dpl_p = (1 - depol)  / (1 + depol / 2)
    #dpl_q = (1 + depol)  / (1 - depol)
    dpl_r = (1 - 2depol) / (1 - depol)
  
    α  =  FT[0.0, 0.0,             3dpl_p]
    β  =  FT[1.0, 0.0,             0.5 * dpl_p]
    γ  =  FT[0.0, 0.0,             dpl_p * sqrt(1.5)] 
    δ  =  FT[0.0, dpl_p * dpl_r * 1.5, 0.0] 
    ϵ  =  FT[0.0, 0.0,             0.0] 
    ζ  =  FT[0.0, 0.0,             0.0]
    return GreekCoefs(α, β, γ, δ, ϵ, ζ)
end

"""
    $(FUNCTIONNAME)(k, an, bn, w)

Compute size-averaged scattering and extinction cross sections from Mie coefficients.

Uses Bohren & Huffman eq. 4.61: ``C_{\\mathrm{sca}} = \\frac{2\\pi}{k^2}\\sum_n (2n+1)(|a_n|^2+|b_n|^2)``
and eq. 4.62 for extinction. Averaging is over the size distribution with weights `w`.

# Arguments
- `k`: wavenumber (2π/λ)
- `an`, `bn`: Mie coefficients, shape `(nquad_radius, N_max)`
- `w`: probability weights for each radius

# Returns
- `(C_sca, C_ext)`: tuple of scattering and extinction cross sections
"""
function compute_avg_C_scatt_ext(k, an, bn, w)
    n_ = collect(1:size(an)[2]);
    n_ = 2n_ .+ 1
    coef = 2π / k^2 * n_'
    return (coef * (w' * (abs2.(an') + abs2.(bn'))')', coef * (w' * real(an + bn))')
end

"""
    $(FUNCTIONNAME)(size_distribution, wᵣ, r, r_max)

Compute probability weights for radii in the size distribution quadrature.

Combines PDF values at quadrature points with Gauss-Legendre weights `wᵣ`, then normalizes.
Used for averaging Mie properties over the size distribution. Logs the fraction of the
distribution cut by `r_max`.

# Arguments
- `size_distribution`: distribution object (e.g. from Distributions.jl)
- `wᵣ`: quadrature weights for radius points
- `r`: radius quadrature points
- `r_max`: maximum radius (truncation)

# Returns
- Normalized weights summing to 1
"""
function compute_wₓ(size_distribution, wᵣ, r, r_max) 
    
    wₓ = pdf.(size_distribution,r)      # Weights from distribution
    wₓ .*= wᵣ                           # pre multiply with wᵣ to get proper means eventually:

    # normalize (could apply a check whether cdf.(size_distribution,r_max) is larger than 0.99:
    #println("Test")
    @debug "Fraction of size distribution cut by max radius: $((1-cdf.(size_distribution,r_max))*100) %"  
    wₓ /= sum(wₓ)
    return wₓ
end

#####
##### Π-matrix construction methods (Sanghavi 2014, eq. 15)
#####

"""
    $(FUNCTIONNAME)(mo::Stokes_IQUV, P, R, T, l::Int, m::Int; sign_change=false)
Compute Π matrix for all stokes vector elements used in computations of the phase matrix 
See Sanghavi 2014, eq. 15
"""
function construct_Π_matrix(mo::Stokes_IQUV, P, R, T, l::Int, m::Int; sign_change=false)
    if sign_change # (basically gets it for -μ due to symmetries on P,R,T)
        if isodd(l-m)
            Π = [SMatrix{4,4}([-P[i,l,m] 0 0 0 ; 0 -R[i,l,m] -T[i,l,m] 0; 0 -T[i,l,m] -R[i,l,m] 0; 0 0 0 -P[i,l,m]]) for i in 1:size(P, 1)] 
        else
            Π = [SMatrix{4,4}([P[i,l,m] 0 0 0 ; 0 R[i,l,m] T[i,l,m] 0; 0 T[i,l,m] R[i,l,m] 0; 0 0 0 P[i,l,m]]) for i in 1:size(P, 1)]
        end
    else
        Π = [SMatrix{4,4}([P[i,l,m] 0 0 0 ; 0 R[i,l,m] -T[i,l,m] 0; 0 -T[i,l,m] R[i,l,m] 0; 0 0 0 P[i,l,m]]) for i in 1:size(P, 1)]
    end
    return Π
end

"""
    $(FUNCTIONNAME)(mo::Stokes_IQU, P, R, T, l::Int, m::Int; sign_change=false)
Compute Π matrix for  stokes vector elements I,Q,U used in computations of the phase matrix
See Sanghavi 2014, eq. 15
"""
function construct_Π_matrix(mo::Stokes_IQU, P, R, T, l::Int, m::Int; sign_change=false)
    if sign_change # (basically gets it for -μ due to symmetries on P,R,T)
        if isodd(l-m)
            Π = [SMatrix{3,3}([-P[i,l,m] 0 0  ; 0 -R[i,l,m] -T[i,l,m] ; 0 -T[i,l,m] -R[i,l,m] ]) for i in 1:size(P, 1)] 
        else
            Π = [SMatrix{3,3}([P[i,l,m] 0 0  ; 0 R[i,l,m] T[i,l,m] ; 0 T[i,l,m] R[i,l,m] ]) for i in 1:size(P, 1)]
        end
    else
        Π = [SMatrix{3,3}([P[i,l,m] 0 0  ; 0 R[i,l,m] -T[i,l,m] ; 0 -T[i,l,m] R[i,l,m] ]) for i in 1:size(P, 1)]
    end
    return Π
end

"""
    $(FUNCTIONNAME)(mo::Stokes_IQ, P, R, T, l::Int, m::Int; sign_change=false)
Compute Π matrix for Stokes vector elements I,Q used in computations of the
phase matrix. This is the I/Q block of Sanghavi 2014, eq. 15.
"""
function construct_Π_matrix(mo::Stokes_IQ, P, R, T, l::Int, m::Int; sign_change=false)
    if sign_change && isodd(l - m)
        return [SMatrix{2,2}([-P[i, l, m] 0;
                              0 -R[i, l, m]]) for i in 1:size(P, 1)]
    end
    return [SMatrix{2,2}([P[i, l, m] 0;
                          0 R[i, l, m]]) for i in 1:size(P, 1)]
end

"""
    $(FUNCTIONNAME)(mo::Stokes_I, P, R, T, l::Int, m::Int; sign_change=false)
Compute Π matrix for  stokes vector elements I used in computations of the phase matrix 

"""
construct_Π_matrix(mod::Stokes_I, P, R, T, l::Int, m::Int; sign_change=false) = sign_change ? -P[:,l,m] : P[:,l,m]

#####
##### B-matrix construction methods (Sanghavi 2014, eq. 16)
#####

"""
    $(FUNCTIONNAME)(mod::Stokes_IQUV, α, β, γ, δ, ϵ, ζ, l::Int)
Compute B matrix for all stokes vector elements used in computations of the phase matrix 
See Sanghavi 2014, eq. 16 
"""
construct_B_matrix(mod::Stokes_IQUV, α, β, γ, δ, ϵ, ζ, l::Int) = SMatrix{4,4}([β[l] γ[l] 0 0 ; γ[l] α[l] 0 0; 0 0 ζ[l] ϵ[l]; 0 0 -ϵ[l] δ[l]])

"""
    $(FUNCTIONNAME)(mod::Stokes_IQU, α, β, γ, δ, ϵ, ζ, l::Int)
Compute B matrix for stokes vector elements I,Q,U used in computations of the phase matrix
    See Sanghavi 2014, eq. 16 
"""
construct_B_matrix(mod::Stokes_IQU, α, β, γ, δ, ϵ, ζ, l::Int) = SMatrix{3,3}([β[l] γ[l] 0 ; γ[l] α[l] 0 ; 0 0 ζ[l]])

"""
    $(FUNCTIONNAME)(mod::Stokes_IQ, α, β, γ, δ, ϵ, ζ, l::Int)
Compute B matrix for Stokes vector elements I,Q used in computations of the
phase matrix. This is the I/Q block of Sanghavi 2014, eq. 16.
"""
construct_B_matrix(mod::Stokes_IQ, α, β, γ, δ, ϵ, ζ, l::Int) = SMatrix{2,2}([β[l] γ[l]; γ[l] α[l]])

"""
$(FUNCTIONNAME)(mod::Stokes_I, α, β, γ, δ, ϵ, ζ, l::Int)
Compute Π matrix for stokes vector elements I used in computations of the phase matrix
See Sanghavi 2014, eq. 16 
"""
construct_B_matrix(mod::Stokes_I, α, β, γ, δ, ϵ, ζ, l::Int) = β[l]


#=
"""
    $(FUNCTIONNAME)(mod::AbstractPolarizationType, μ, α, β, γ, δ, ϵ, ζ, m::Int)
Compute moments of the phase matrix 
"""
function compute_Z_moments(mod::AbstractPolarizationType, μ, greek_coefs::GreekCoefs, m::Int ; arr_type = Array)
    @unpack α, β, γ, δ, ϵ, ζ = greek_coefs
    FT = eltype(β)
    n = length(μ)

    # Change from 0-index to 1-index (i.e. the lowest m is 0 ), 
    # make more logical later to avoid confusion later (m=0 has a meaning!!)
    m = m+1
    
    # Set prefactor for moments (note 1-notation for `m` here):
    fact = (m == 1) ? 0.5 : 1.0

    # get l_max just from length of array:
    l_max = length(β)

    # Check that all μ are positive here ([0,1])
    # @show μ
    @assert all(0 .< μ .≤ 1) "all μ's within compute_Z_moments have to be ∈ ]0,1]"

    # Compute legendre Polynomials at μ and up to lmax
    P, R, T    = Scattering.compute_associated_legendre_PRT(μ, l_max)
    P⁻, R⁻, T⁻ = Scattering.compute_associated_legendre_PRT(-μ, l_max)
  
    # Pre-compute all required B matrices
    𝐁_all = [construct_B_matrix(mod, α, β, γ, δ, ϵ, ζ, i) for i in 1:l_max]

    # Get dimension of square matrix (easier for Scalar/Stokes dimensions)
    B_dim = Int(sqrt(length(𝐁_all[1])))
    
    # Create matrices:
    nb = B_dim * n
    𝐙⁺⁺, 𝐙⁻⁺ = (zeros(FT, nb, nb), zeros(FT, nb, nb))
    A⁺⁺, A⁻⁺ = (zeros(FT, B_dim, B_dim, n, n), zeros(FT, B_dim, B_dim, n, n))

    # Iterate over l
    @inbounds for l = m:l_max

        # B matrix for l
        𝐁 = 𝐁_all[l];

        # Construct Π matrix for l,m pair (change to in place later!)
        # See eq. 15 in Sanghavi 2014, note that P,R,T are already normalized
        Π  = construct_Π_matrix(mod, P, R, T, l, m)
        Π⁻ = construct_Π_matrix(mod, P⁻, R⁻, T⁻, l, m)

        # Iterate over angles
        for i in eachindex(μ), j in eachindex(μ)
            if B_dim == 1
                A⁺⁺[B_dim,B_dim,i,j] += Π[i] * 𝐁 * Π[j]
                A⁻⁺[B_dim,B_dim,i,j] += Π[i] * 𝐁 * Π⁻[j]
            else
                A⁺⁺[:,:,i,j] += Π[i] * 𝐁 * Π[j]
                A⁻⁺[:,:,i,j] += Π[i] * 𝐁 * Π⁻[j]
            end
        end
    end

    # Now get to the Z part:
    @inbounds for imu in eachindex(μ), jmu in eachindex(μ)
        
        # Indices adjusted for size of A
        ii, jj = ((imu - 1) * B_dim, (jmu - 1) * B_dim)
            
        # This is equivalent to Z̄ = 1/(1+δ) * C̄m+S̄m = 1/(1+δ) * (A+DAD+AD-DA) 
        # (see eq 11 in Sanghavi et al, 2013)
        for i in 1:B_dim, j in 1:B_dim
            𝐙⁺⁺[ii + i,jj + j] = 2fact * A⁺⁺[i,j,imu,jmu]
            if i <= 2 && j >= 3
                𝐙⁻⁺[ii + i,jj + j] = -2fact * A⁻⁺[i,j,imu,jmu]
            elseif i >= 3 && j <= 2
                𝐙⁻⁺[ii + i,jj + j] = -2fact * A⁻⁺[i,j,imu,jmu]
            else
                𝐙⁻⁺[ii + i,jj + j] = 2fact * A⁻⁺[i,j,imu,jmu]
            end
        end
    end

    # Return Z-moments
    return arr_type(𝐙⁺⁺), arr_type(𝐙⁻⁺)
end
=#
