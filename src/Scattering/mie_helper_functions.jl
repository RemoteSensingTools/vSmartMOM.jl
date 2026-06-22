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

"""
    _mie_dn_recursion!(Dn, y, nmx)

Perform the downward Dₙ recursion (BH eq. 4.89) for `y = size_param * refractive_idx`.

Two dispatch methods:
- `Complex{<:AbstractFloat}` (plain floats): stabilised in Float64 to avoid
  catastrophic cancellation for `|y| ≳ 50` in Float32.
- `Complex{<:ForwardDiff.Dual}`: uses native arithmetic so derivative tangents
  propagate correctly through the recursion; no Float64 cast is possible.
"""
function _mie_dn_recursion!(Dn, y::Complex{<:AbstractFloat}, nmx)
    y64 = Complex{Float64}(y)
    Dn_prev = Complex{Float64}(0)
    @inbounds for n = (nmx - 1):-1:1
        ratio = (n + 1) / y64
        Dn_prev = ratio - 1 / (Dn_prev + ratio)
        Dn[n] = Dn_prev
    end
end

function _mie_dn_recursion!(Dn, y, nmx)
    # Generic path: used for ForwardDiff Dual numbers and any other numeric type
    # where conversion to Complex{Float64} is not defined.
    Dn_prev = zero(y)
    @inbounds for n = (nmx - 1):-1:1
        ratio = (n + 1) / y
        Dn_prev = ratio - 1 / (Dn_prev + ratio)
        Dn[n] = Dn_prev
    end
end

@doc raw"""
    $(FUNCTIONNAME)(size_param, refractive_idx::Number, an, bn, Dn)

Compute Mie coefficients ``a_n`` and ``b_n`` from the size parameter and complex
refractive index (Bohren & Huffman, eq. 4.88).

The logarithmic derivative ``D_n(y)`` is obtained via downward recursion
(BH eq. 4.89; de Rooij & Stap 1984, eq. A9).  For plain
`Complex{<:AbstractFloat}` arguments the recursion is promoted to Float64
for numerical stability (cancellation-prone continued-fraction form loses
significant digits in Float32 for ``|y| \gtrsim 50``).  For
`Complex{<:ForwardDiff.Dual}` arguments native arithmetic is used so that
derivative tangents propagate correctly.

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
    # For plain floating-point types the recursion is performed in Float64
    # for numerical stability (cancellation-prone continued-fraction form
    # loses significant digits in Float32 for |y| ≳ 50).
    # For ForwardDiff Dual types we must keep native arithmetic so that
    # derivatives ∂aₙ/∂n propagate correctly through Dₙ.
    _mie_dn_recursion!(Dn, y, nmx)

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


