#=
=====================================================================
Siewert NAI-2 aerosol optical properties (cross-sections + Greek
coefficients of the polarized phase matrix).
=====================================================================
For a given Mie aerosol size distribution, this file computes:

 1. Bulk extinction Cₑₓₜ and single-scattering Cₛ cross-sections, and
    the single-scattering albedo ϖ̃ = Cₛ / Cₑₓₜ.
 2. The Greek expansion coefficients α(l), β(l), γ(l), δ(l), ϵ(l),
    ζ(l) of the polarized scattering matrix in the Hovenier basis. β is
    the same as the scalar Legendre expansion of the phase function;
    α…ζ together encode the full 4×4 Mueller matrix.
 3. The δ-M truncation factor fᵗ when a δBGE-style truncation is
    requested (otherwise fᵗ = 0). With NoTruncation(), the raw
    coefficients are returned as-is.

NAI-2 = "Numerical Algorithm for Inelastic-2" in Siewert's notation.
The companion `compute_PCW.jl` uses precomputed Wigner symbol tables
for the Greek expansion — same physics, different numerical strategy.

Reference: Siewert (1982) ApJ 245, 357; (2000a) JQSRT 64, 109.
=====================================================================
=#

@doc raw"""
    compute_aerosol_optical_properties(model::MieModel{<:NAI2,FT}, FT_out::Type) -> AerosolOptics

Compute bulk aerosol optical properties with the Siewert NAI-2 formulation
on the CPU.

This method is the internal two-argument form called by the architecture/method
router (`_dispatch_aerosol_optics`) in `phase_function_autodiff.jl`. The router
supplies `FT_out = _mie_output_type(model)` (the model's own `FT` type parameter)
so that a Float32 model produces Float32 greek coefficients. For
`architecture = GPU()` the GPU KernelAbstractions pipeline is dispatched
instead (see [`compute_aerosol_optical_properties_gpu`](@ref)).

The second argument `FT_out` sets the **output** float type of the returned
`AerosolOptics` (greek coefficients, ω̃, k). It is supplied by the router and
must be passed explicitly — there is no default.

**Internal precision**: the Mie cross-section sums (S₁, S₂, Cₑₓₜ, Cₛ arrays)
and the Greek-coefficient angular quadrature always run in `Float64` on the CPU
path, regardless of `FT`. This is intentional — the Dₙ downward continued-fraction
recursion and the size-integral sums are catastrophically unstable in Float32 for
large or absorbing particles (|y| ≳ 50). The final `convert.(FT_out, …)` call
narrows back to the requested output precision.

Reference:
- S. Sanghavi, *Revisiting the Fourier expansion of Mie scattering matrices in generalized spherical functions*, JQSRT 136 (2014), 16-27. https://doi.org/10.1016/j.jqsrt.2013.12.015
- Cross-section expressions are the standard Mie-series forms (paper Eq. (1)); Greek-coefficient Fourier expansion follows the Siewert/Domke framework discussed around paper Eq. (17).

For each radius quadrature point, Mie coefficients are used to compute
extinction and scattering cross-sections:

```math
C_{\mathrm{sca}}(r) = \frac{2\pi}{k^2}\sum_{n=1}^{N}(2n+1)\left(|a_n|^2 + |b_n|^2\right),
\qquad
C_{\mathrm{ext}}(r) = \frac{2\pi}{k^2}\sum_{n=1}^{N}(2n+1)\Re(a_n+b_n),
```

with ``k=2\pi/\lambda``. Bulk values are quadrature-weighted means over the
size distribution. The returned optical scalars are:

```math
\tilde{\omega} = \frac{\bar{C}_{\mathrm{sca}}}{\bar{C}_{\mathrm{ext}}},
\qquad
k = \bar{C}_{\mathrm{ext}}.
```

Greek coefficients `(α, β, γ, δ, ϵ, ζ)` are then derived from the bulk phase
matrix moments (Sanghavi, 2014).

# Returns
- [`AerosolOptics`](@ref) with `greek_coefs`, `ω̃`, `k`, and `fᵗ = 1`.

# Notes
- Uses the convention `nᵢ >= 0`.
- Radius and wavelength units must be consistent.
"""
function compute_aerosol_optical_properties(model::MieModel{FDT,FT}, FT_out::Type) where {FDT <: NAI2, FT}

    # Unpack the model
    (; computation_type, aerosol, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B) = model

    # Extract variables from aerosol struct:
    (; size_distribution, nᵣ, nᵢ) = aerosol

    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0

    # Internal computation precision for the CPU Mie path.
    #
    # The Dₙ downward continued-fraction recursion and the size-distribution
    # quadrature sums (S₁, S₂, Cₑₓₜ, Cₛ, Greek-coef integrals) are performed
    # in Float64 regardless of the model's FT.  This is the historically
    # deliberate requirement: Float32 Dₙ recursion loses ≳ 4 significant digits
    # for |y| ≳ 50 (absorbing/large particles), causing catastrophic error in aₙ/bₙ.
    # For ForwardDiff Dual types we keep native arithmetic (Float64 promotion is
    # not meaningful for Dual; the _mie_dn_recursion! generic path handles it).
    IC = FT <: AbstractFloat ? Float64 : FT  # IC = Internal Computation type

    # Get radius quadrature points and weights using log-space quadrature
    # (equidistant in ln(r) for efficient integration of log-normal distributions)
    r_min = max(quantile(size_distribution, 1e-8), 1e-6 * IC(r_max))
    r, wᵣ = gauleg_log(nquad_radius, IC(r_min), IC(r_max); norm=false)

    # Wavenumber (in IC precision)
    k = IC(2π) / IC(λ)

    # Size parameter
    x_size_param = k .* r # (2πr/λ)

    # Compute n_max for largest size:
    n_max = get_n_max(maximum(x_size_param))

    # Determine max amount of Gaussian quadrature points for angle dependence of
    # phase functions:
    n_mu = 2n_max - 1;

    # Obtain Gauss-Legendre quadrature points and weights for phase function
    μ, w_μ = gausslegendre(n_mu)
    # Promote μ/w_μ to IC (they come out as Float64 from gausslegendre, which is
    # already IC when FT<:AbstractFloat; no-op for Dual).
    μ  = IC.(μ)
    w_μ = IC.(w_μ)

    # Compute π and τ functions
    leg_π, leg_τ = compute_mie_π_τ(μ, n_max)

    # Pre-allocate internal arrays in IC (Float64 for plain-float models).
    # These carry the full-precision Mie amplitudes and cross-section integrands;
    # they are converted to FT_out only at the very end.
    S₁    = zeros(Complex{IC}, n_mu, nquad_radius)
    S₂    = zeros(Complex{IC}, n_mu, nquad_radius)
    f₁₁   = zeros(IC, n_mu, nquad_radius)
    f₃₃   = zeros(IC, n_mu, nquad_radius)
    f₁₂   = zeros(IC, n_mu, nquad_radius)
    f₃₄   = zeros(IC, n_mu, nquad_radius)
    C_ext = zeros(IC, nquad_radius)
    C_sca = zeros(IC, nquad_radius)

    # Standardized weights for the size distribution (in IC):
    wₓ = compute_wₓ(size_distribution, wᵣ, r, IC(r_max))

    # Pre-allocate buffers for the inner loop (sized for the largest particle).
    # m_ref is built in IC so the Complex arithmetic in compute_mie_ab! stays
    # in full precision regardless of the aerosol's nᵣ/nᵢ element type.
    m_ref = IC(nᵣ) - IC(nᵢ) * im
    y_max = maximum(x_size_param) * abs(m_ref)
    nmx_max = round(Int, max(n_max, y_max) + 51)
    an  = zeros(Complex{IC}, n_max)
    bn  = zeros(Complex{IC}, n_max)
    Dₙ  = zeros(Complex{IC}, nmx_max)
    n_  = IC.(2 .* collect(1:n_max) .+ 1)

    # Loop over size parameters
    for i = 1:length(x_size_param)

        n_max_i = get_n_max(x_size_param[i])

        # Zero the views that will be used
        an_v = view(an, 1:n_max_i)
        bn_v = view(bn, 1:n_max_i)
        fill!(an_v, 0)
        fill!(bn_v, 0)
        fill!(Dₙ, 0)

        # Compute aₙ,bₙ and S₁,S₂
        compute_mie_ab!(x_size_param[i], m_ref, an_v, bn_v, Dₙ)
        compute_mie_S₁S₂!(an_v, bn_v, leg_π, leg_τ, view(S₁, :, i), view(S₂, :, i))
        
        # Compute Extinction and scattering cross sections using pre-allocated n_
        n_v = view(n_, 1:n_max_i)
        @inbounds C_sca[i] = 2π / k^2 * dot(n_v, abs2.(an_v) .+ abs2.(bn_v))
        @inbounds C_ext[i] = 2π / k^2 * dot(n_v, real.(an_v .+ bn_v))
       
        # Compute scattering matrix components per size parameter (in-place).
        # inv_x2 is computed in IC (Float64) so all products stay in full precision.
        inv_x2 = IC(0.5) / x_size_param[i]^2
        @inbounds for iμ in 1:n_mu
            s1 = S₁[iμ, i]; s2 = S₂[iμ, i]
            f₁₁[iμ, i] =  inv_x2 * real(abs2(s1) + abs2(s2))
            f₃₃[iμ, i] =  inv_x2 * real(s1 * conj(s2) + s2 * conj(s1))
            f₁₂[iμ, i] = -inv_x2 * real(abs2(s1) - abs2(s2))
            f₃₄[iμ, i] = -inv_x2 * imag(s1 * conj(s2) - s2 * conj(s1))
        end
    end

    # Calculate bulk scattering and extinction cross-sections
    bulk_C_sca =  sum(wₓ .* C_sca)
    bulk_C_ext =  sum(wₓ .* C_ext)
    
    # Compute bulk scattering 
    wr = (4π * r.^2 .*  wₓ) 
    bulk_f₁₁   =  f₁₁ * wr
    bulk_f₃₃   =  f₃₃ * wr
    bulk_f₁₂   =  f₁₂ * wr
    bulk_f₃₄   =  f₃₄ * wr

    # Normalize Phase function with bulk scattering cross section.
    bulk_f₁₁ /= bulk_C_sca 
    bulk_f₃₃ /= bulk_C_sca
    bulk_f₁₂ /= bulk_C_sca
    bulk_f₃₄ /= bulk_C_sca
    
    # Range of l-values
    l_max = length(μ);

    # Get legendre polynomials for l-max 
    P, P², R², T² = compute_legendre_poly(μ, l_max)

    # Compute Greek coefficients in IC (Float64 for plain-float models).
    # All angular quadrature products (w_μ, bulk_fXX, P, P², R², T²) are already IC.
    α = zeros(IC, l_max)
    β = zeros(IC, l_max)
    δ = zeros(IC, l_max)
    γ = zeros(IC, l_max)
    ϵ = zeros(IC, l_max)
    ζ = zeros(IC, l_max)

    # Compute Greek coefficients from bulk scattering matrix elements (spherical only here!)
    for l = 0:length(β) - 1

        # pre-factor:
        fac = l ≥ 2 ? (2l + 1) / 2 * sqrt(1 / ((l - 1) * (l) * (l + 1) * (l + 2))) : 0

        # Compute Greek coefficients
        # Eq 17 in Sanghavi 2014

        δ[l + 1] = (2l + 1) / 2 * w_μ' * (bulk_f₃₃ .* P[:,l + 1])
        β[l + 1] = (2l + 1) / 2 * w_μ' * (bulk_f₁₁ .* P[:,l + 1])
        γ[l + 1] = fac      * w_μ' * (bulk_f₁₂ .* P²[:,l + 1])
        ϵ[l + 1] = fac      * w_μ' * (bulk_f₃₄ .* P²[:,l + 1])
        ζ[l + 1] = fac      * w_μ' * (bulk_f₃₃ .* R²[:,l + 1] + bulk_f₁₁ .* T²[:,l + 1])
        α[l + 1] = fac      * w_μ' * (bulk_f₁₁ .* R²[:,l + 1] + bulk_f₃₃ .* T²[:,l + 1])
    end

    # Convert internal Float64 results to the requested output type FT_out.
    # For ForwardDiff Dual types (IC == FT, not AbstractFloat), skip conversion
    # to preserve derivative tangents through the computation.
    if IC <: AbstractFloat && FT_out <: AbstractFloat
        greek_coefs = GreekCoefs(convert.(FT_out, α),
                                 convert.(FT_out, β),
                                 convert.(FT_out, γ),
                                 convert.(FT_out, δ),
                                 convert.(FT_out, ϵ),
                                 convert.(FT_out, ζ))
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=FT_out(bulk_C_sca / bulk_C_ext), k=FT_out(bulk_C_ext), fᵗ=FT_out(1))
    else
        greek_coefs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=(bulk_C_sca / bulk_C_ext), k=(bulk_C_ext), fᵗ=one(eltype(α)))
    end
end

@doc raw"""
    compute_ref_aerosol_extinction(model::MieModel{<:NAI2,FT}, FT_out::Type) -> Real

Compute only the bulk extinction cross-section for the aerosol model:

```math
\bar{C}_{\mathrm{ext}} = \int C_{\mathrm{ext}}(r)\,p(r)\,dr,
```

where ``C_{\mathrm{ext}}(r)`` is evaluated from the Mie series and ``p(r)`` is
represented with quadrature over `[0, r_max]`.

The second argument `FT_out` sets the return type. It is supplied by the caller
and must be passed explicitly — there is no default.  The internal cross-section
computation always runs in Float64 for plain floats (same numerical stability
rationale as `compute_aerosol_optical_properties`).

This helper avoids phase-function and Greek-coefficient reconstruction.

Reference: Sanghavi (2014), Eq. (1) for `C_ext(r)`.
"""
function compute_ref_aerosol_extinction(model::MieModel{FDT,FT}, FT_out::Type) where {FDT <: NAI2, FT}

    # Unpack the model
    (; computation_type, aerosol, λ, polarization_type, r_max, nquad_radius) = model

    # Extract variables from aerosol struct:
    (; size_distribution, nᵣ, nᵢ) = aerosol

    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (definition)"

    # Internal computation type: Float64 for plain floats (stability), native for Dual.
    # See compute_aerosol_optical_properties for the detailed rationale.
    IC = FT <: AbstractFloat ? Float64 : FT

    # Get radius quadrature points and weights using log-space quadrature
    r_min = max(quantile(size_distribution, 1e-8), 1e-6 * IC(r_max))
    r, wᵣ = gauleg_log(nquad_radius, IC(r_min), IC(r_max); norm=false)

    # Wavenumber (IC precision)
    k = IC(2π) / IC(λ)

    # Size parameter
    x_size_param = k .* r # (2πr/λ)

    # Compute n_max for largest size:
    n_max = get_n_max(maximum(x_size_param))

    # Pre-allocate C_ext in IC (Float64 for plain-float models)
    C_ext = zeros(IC, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, IC(r_max))

    # m_ref in IC so the Mie recursion stays numerically stable
    m_ref = IC(nᵣ) - IC(nᵢ) * im

    # Pre-allocate buffers outside the loop (sized for largest particle)
    y_max = maximum(x_size_param) * abs(m_ref)
    nmx_max = round(Int, max(n_max, y_max) + 51)
    an = zeros(Complex{IC}, n_max)
    bn = zeros(Complex{IC}, n_max)
    Dₙ = zeros(Complex{IC}, nmx_max)
    n_ = IC.(2 .* collect(1:n_max) .+ 1)

    # Loop over size parameters
    for i = 1:length(x_size_param)

        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max_i = get_n_max(x_size_param[i])

        an_v = view(an, 1:n_max_i)
        bn_v = view(bn, 1:n_max_i)
        fill!(an_v, 0); fill!(bn_v, 0)
        fill!(Dₙ, 0)

        n_v = view(n_, 1:n_max_i)

        # Compute aₙ,bₙ (Dₙ recursion runs in Float64 via _mie_dn_recursion!)
        compute_mie_ab!(x_size_param[i], m_ref, an_v, bn_v, Dₙ)

        # Compute Extinction cross section:
        C_ext[i] = IC(2π) / k^2 * (n_v' * real.(an_v .+ bn_v))
    end

    # Calculate bulk extinction cross-section (IC precision)
    bulk_C_ext = sum(wₓ .* C_ext)

    # Convert to requested output type (narrow Float64 → FT_out for plain floats)
    return FT_out <: AbstractFloat ? FT_out(bulk_C_ext) : bulk_C_ext
end

@doc raw"""
    phase_function(aerosol::Aerosol, λ, r_max, nquad_radius)
        -> (μ, w_μ, p11, C_ext, C_sca, g)

Compute scalar phase-function output for a size-distributed aerosol.

The returned `p11(μ)` is normalized by the bulk scattering cross-section and
satisfies:

```math
\frac{1}{4\pi}\int_0^{2\pi}\int_{-1}^{1} p_{11}(\mu)\,d\mu\,d\phi = 1.
```

Asymmetry factor:

```math
g = \frac{1}{2}\int_{-1}^{1}\mu\,p_{11}(\mu)\,d\mu.
```

Reference: Sanghavi (2014), Mie-amplitude/cross-section setup from Eq. (1), with NAI-style angular integration discussed in Secs. 3-4.
"""
function phase_function(aerosol::Aerosol, λ, r_max, nquad_radius) 

    # Extract variables from aerosol struct:
    (; size_distribution, nᵣ, nᵢ) = aerosol
    
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (definition)"

    # Get the refractive index's real part type
    FT = eltype(nᵣ);
    # @assert FT == Float64 "Aerosol computations require 64bit"
    # Get radius quadrature points and weights using log-space quadrature
    r_min = max(quantile(size_distribution, 1e-8), 1e-6 * r_max)
    r, wᵣ = gauleg_log(nquad_radius, r_min, r_max; norm=true)

    # Wavenumber
    k = 2π / λ

    # Size parameter
    x_size_param = k * r # (2πr/λ)

    # Compute n_max for largest size:
    n_max = get_n_max(maximum(x_size_param))

    # Determine max amount of Gaussian quadrature points for angle dependence of 
    # phase functions:
    n_mu = 2n_max - 1;

    # Obtain Gauss-Legendre quadrature points and weights for phase function
    μ, w_μ = gausslegendre(n_mu)

    # Compute π and τ functions
    leg_π, leg_τ = compute_mie_π_τ(μ, n_max)

    # Pre-allocate arrays:
    S₁    = zeros(Complex{FT}, n_mu, nquad_radius)
    S₂    = zeros(Complex{FT}, n_mu, nquad_radius)
    f₁₁   = zeros(FT, n_mu, nquad_radius)
    C_ext = zeros(FT, nquad_radius)
    C_sca = zeros(FT, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    
    # Loop over size parameters
    for i = 1:length(x_size_param)

        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = get_n_max(x_size_param[i])

        # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
        an = (zeros(Complex{FT}, n_max))
        bn = (zeros(Complex{FT}, n_max))

        # Weighting for sums of 2n+1
        n_ = collect(FT, 1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dn:
        y = x_size_param[i] * (aerosol.nᵣ - aerosol.nᵢ);
        nmx = round(Int, max(n_max, abs(y)) + 51)
        Dn = zeros(Complex{FT}, nmx)

        # Compute an,bn and S₁,S₂
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ - aerosol.nᵢ * im, an, bn, Dn)
        compute_mie_S₁S₂!(an, bn, leg_π, leg_τ, view(S₁, :, i), view(S₂, :, i))
        
        # Compute Extinction and scattering cross sections: 
        C_sca[i] = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))
        # @show r[i], x_size_param[i], C_ext[i], C_ext[i]/(4π*r[i]^2), C_ext[i]*1e-8
        # Compute scattering matrix components per size parameter (might change column/row ordering):
        f₁₁[:,i] =  0.5 / x_size_param[i]^2  * real(abs2.(S₁[:,i]) + abs2.(S₂[:,i]));
    end

    # Calculate bulk scattering and extinction coeffitientcs
    bulk_C_sca =  sum(wₓ .* C_sca)
    bulk_C_ext =  sum(wₓ .* C_ext)
    
    # Compute bulk scattering 
    wr = (4π * r.^2 .*  wₓ) 
    bulk_f₁₁   =  f₁₁ * wr

    # Normalize Phase function with bulk scattering cross section.
    bulk_f₁₁ /= bulk_C_sca 

    # Assymetry factor g = 0.5\int_{-1}^1 μ P(μ) dμ
    g = 1/2 * w_μ'*(μ .*  bulk_f₁₁ )
    return μ, w_μ, bulk_f₁₁, bulk_C_ext, bulk_C_sca, g
end

"""
    phase_function(r::FT, λ::FT, nᵣ::FT, nᵢ::FT) where {FT<:AbstractFloat}
        -> (μ, w_μ, p11, C_ext, C_sca, g)

Monodisperse version of [`phase_function`](@ref), for a single radius `r`.
Returns the same tuple as the size-distribution overload.
"""
function phase_function(r::FT, λ::FT, nᵣ::FT, nᵢ::FT) where {FT<:AbstractFloat} 
    # Imaginary part of the refractive index must be ≥ 0 (definition)
    @assert nᵢ ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (definition)"
    # Wavenumber
    k = 2π / λ  

    # Size parameter
    size_param = k * r # (2πr/λ)
    
    # Compute n_max (doubled here to make plots smoother):
    n_max = 2get_n_max(size_param)

    # Determine max amount of Gaussian quadrature points for angle dependence of 
    # phase functions:
    n_mu = 2n_max - 1;

    # Obtain Gauss-Legendre quadrature points and weights for phase function
    μ, w_μ = gausslegendre(n_mu)

    # Compute π and τ functions
    leg_π, leg_τ = compute_mie_π_τ(μ, n_max)

    # Pre-allocate arrays:
    S₁    = zeros(Complex{FT}, n_mu)
    S₂    = zeros(Complex{FT}, n_mu)

    # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
    an = (zeros(Complex{FT}, n_max))
    bn = (zeros(Complex{FT}, n_max))

    # Weighting for sums of 2n+1
    n_ = collect(FT, 1:n_max);
    n_ = 2n_ .+ 1

    # Pre-allocate Dn:
    y = size_param * (nᵣ - nᵢ);
    nmx = round(Int, max(n_max, abs(y)) + 51)
    Dn = zeros(Complex{FT}, nmx)

    # Compute an,bn and S₁,S₂
    compute_mie_ab!(size_param, nᵣ - nᵢ * im, an, bn, Dn)
    compute_mie_S₁S₂!(an, bn, leg_π, leg_τ, S₁, S₂)
        
    # Compute Extinction and scattering cross sections: 
    C_sca = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
    C_ext = 2π / k^2 * (n_' * real(an + bn))

    # Compute scattering matrix components per size parameter (might change column/row ordering):
    f₁₁ =  0.5 / size_param^2  * real(abs2.(S₁) + abs2.(S₂));
    f₁₁ *= 4π * r.^2
    f₁₁ /= C_sca

    # Asymmetry factor g
    g = 1/2 * w_μ'*(μ .*  f₁₁ )
    return μ, w_μ, f₁₁, C_ext, C_sca, g
end

@doc raw"""
    compute_aerosol_XS(aerosol::Aerosol, λ, r_max, nquad_radius)
        -> (XS_ext, XS_sca, Cext_eff, Csca_eff)

Compute bulk extinction/scattering cross-sections for a size-distributed aerosol
without reconstructing phase matrices.

# Returns
- `XS_ext`: bulk extinction cross-section ``\bar{C}_{\mathrm{ext}}``
- `XS_sca`: bulk scattering cross-section ``\bar{C}_{\mathrm{sca}}``
- `Cext_eff`: area-normalized extinction
- `Csca_eff`: area-normalized scattering

with

```math
C_{\mathrm{eff}} = \frac{\bar{C}}{\pi\int r^2 p(r)\,dr}.
```

Reference: Sanghavi (2014), Eq. (1) for particle-level `C_ext`/`C_sca`.
"""
function compute_aerosol_XS(aerosol::Aerosol, λ::FT, r_max::FT, nquad_radius::Int64) where {FT<:AbstractFloat} 
    # Extract variables from aerosol struct:
    (; size_distribution, nᵣ, nᵢ) = aerosol
    
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0

    # Get the refractive index's real part type
    #@show size_distribution.σ  
    #FT = eltype(size_distribution.σ);
    # @assert FT == Float64 "Aerosol computations require 64bit"
    # Get radius quadrature points and weights using log-space quadrature
    r_min = max(quantile(size_distribution, 1e-8), 1e-6 * r_max)
    r, wᵣ = gauleg_log(nquad_radius, r_min, r_max; norm=true)

    # Wavenumber
    k = 2π / λ

    # Size parameter
    x_size_param = k * r # (2πr/λ)

    # Compute n_max for largest size:
    n_max = get_n_max(maximum(x_size_param))

    # Determine max amount of Gaussian quadrature points for angle dependence of 
    # phase functions:
    #n_mu = 2n_max - 1;

    # Obtain Gauss-Legendre quadrature points and weights for phase function
    #μ, w_μ = gausslegendre(n_mu)

    # Compute π and τ functions
    #leg_π, leg_τ = compute_mie_π_τ(μ, n_max)

    # Pre-allocate arrays:
    #S₁    = zeros(Complex{FT}, n_mu, nquad_radius)
    #S₂    = zeros(Complex{FT}, n_mu, nquad_radius)
    #f₁₁   = zeros(FT, n_mu, nquad_radius)
    #f₃₃   = zeros(FT, n_mu, nquad_radius)
    #f₁₂   = zeros(FT, n_mu, nquad_radius)
    #f₃₄   = zeros(FT, n_mu, nquad_radius)
    C_ext = zeros(FT, nquad_radius)
    C_sca = zeros(FT, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    
    # Loop over size parameters
    for i = 1:length(x_size_param)

        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = get_n_max(x_size_param[i])

        # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
        an = (zeros(Complex{FT}, n_max))
        bn = (zeros(Complex{FT}, n_max))

        # Weighting for sums of 2n+1
        n_ = collect(FT, 1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dₙ  :
        y = x_size_param[i] * (aerosol.nᵣ - aerosol.nᵢ);
        nmx = round(Int, max(n_max, abs(y)) + 51)
        Dₙ  = zeros(Complex{FT}, nmx)

        # Compute aₙ,bₙ and S₁,S₂
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ - aerosol.nᵢ * im, an, bn, Dₙ )
        #compute_mie_S₁S₂!(an, bn, leg_π, leg_τ, view(S₁, :, i), view(S₂, :, i))
        
        # Compute Extinction and scattering cross sections: 
        C_sca[i] = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))
       
        # Compute scattering matrix components per size parameter:
        #f₁₁[:,i] =  0.5 / x_size_param[i]^2  * real(abs2.(S₁[:,i]) + abs2.(S₂[:,i]));
        #f₃₃[:,i] =  0.5 / x_size_param[i]^2  * real(S₁[:,i] .* conj(S₂[:,i]) + S₂[:,i] .* conj(S₁[:,i]));
        #f₁₂[:,i] = -0.5 / x_size_param[i]^2  * real(abs2.(S₁[:,i]) - abs2.(S₂[:,i]));
        #f₃₄[:,i] = -0.5 / x_size_param[i]^2  * imag(S₁[:,i] .* conj(S₂[:,i]) - S₂[:,i] .* conj(S₁[:,i]));

    end

    # Calculate bulk scattering and extinction cross-sections
    bulk_XS_sca = sum(wₓ .* C_sca)
    bulk_XS_ext = sum(wₓ .* C_ext)

    bulk_C_sca =  bulk_XS_sca/(π * sum(wₓ .* r.^2))
    bulk_C_ext =  bulk_XS_ext/(π * sum(wₓ .* r.^2))
    
    #wr = (4π * r.^2 .*  wₓ) 
    return bulk_XS_ext, bulk_XS_sca, bulk_C_ext, bulk_C_sca 

end
