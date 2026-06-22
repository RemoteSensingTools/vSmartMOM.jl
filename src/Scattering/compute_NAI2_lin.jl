#=
 
This file specifies how to compute aerosol optical properties using the Siewert-NAI2 method
 
=#

"""
    $(FUNCTIONNAME)(lin::LinMode, model::MieModel{FDT,FT}, FT2::Type) where {FDT<:NAI2, FT}

Reference: Suniti Sanghavi 2014, https://doi.org/10.1016/j.jqsrt.2013.12.015

Compute the aerosol optical properties and their derivatives (w.r.t. nᵣ, nᵢ, rₚ, σₚ)
using the Siewert-NAI2 linearized method.

The second argument `FT2` sets the **output** float type for both the forward greek
coefficients and the derivative (lin) greek coefficients. It must be passed
explicitly — there is no default.  The internal computation always runs in Float64 for
numerical stability of the Dₙ recursion (same rationale as the forward path).

Input: LinMode, MieModel
Output: (AerosolOptics, linAerosolOptics)
"""
function compute_aerosol_optical_properties(lin::LinMode, model::MieModel{FDT,FT}, FT2::Type) where {FDT <: NAI2, FT}

    # Unpack the model
    (; computation_type, aerosol, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B) = model

    # Extract variables from aerosol struct:
    (; size_distribution, nᵣ, nᵢ) = aerosol
    
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0

    # Internal computation precision: Float64 for plain floats (numerical stability of
    # the Dₙ downward recursion and S₁/S₂/Cₑₓₜ sums), native arithmetic for Dual types.
    # See compute_aerosol_optical_properties (forward path) for the detailed rationale.
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

    # Determine max amount of Gaussian quadrature points for angle dependence of
    # phase functions:
    n_mu = 2n_max - 1;

    # Obtain Gauss-Legendre quadrature points and weights for phase function
    μ, w_μ = gausslegendre(n_mu)
    μ   = IC.(μ)
    w_μ = IC.(w_μ)

    # Compute π and τ functions
    leg_π, leg_τ = compute_mie_π_τ(μ, n_max)

    # Pre-allocate internal arrays in IC (Float64 for plain-float models).
    # All Mie amplitudes, cross-section integrands, and their linearized derivatives
    # live in IC to avoid precision loss in the continued-fraction recursion.
    S₁    = zeros(Complex{IC}, n_mu, nquad_radius)
    S₂    = zeros(Complex{IC}, n_mu, nquad_radius)
    f₁₁   = zeros(IC, n_mu, nquad_radius)
    f₃₃   = zeros(IC, n_mu, nquad_radius)
    f₁₂   = zeros(IC, n_mu, nquad_radius)
    f₃₄   = zeros(IC, n_mu, nquad_radius)
    C_ext = zeros(IC, nquad_radius)
    C_sca = zeros(IC, nquad_radius)

    # derivatives with respect to nᵣ, nᵢ, respectively
    Ṡ₁    = zeros(Complex{IC}, 2, n_mu, nquad_radius)
    Ṡ₂    = zeros(Complex{IC}, 2, n_mu, nquad_radius)

    ḟ₁₁   = zeros(IC, 2, n_mu, nquad_radius)
    ḟ₃₃   = zeros(IC, 2, n_mu, nquad_radius)
    ḟ₁₂   = zeros(IC, 2, n_mu, nquad_radius)
    ḟ₃₄   = zeros(IC, 2, n_mu, nquad_radius)
    Ċ_ext = zeros(IC, 2, nquad_radius)
    Ċ_sca = zeros(IC, 2, nquad_radius)

    # derivatives with respect to nᵣ, nᵢ, rₚ, σₚ, respectively
    bulk_Ċ_ext = zeros(IC, 4)
    bulk_Ċ_sca = zeros(IC, 4)
    bulk_ϖ̇     = zeros(IC, 4)
    bulk_ḟ₁₁   = zeros(IC, 4, n_mu)
    bulk_ḟ₃₃   = zeros(IC, 4, n_mu)
    bulk_ḟ₁₂   = zeros(IC, 4, n_mu)
    bulk_ḟ₃₄   = zeros(IC, 4, n_mu)

    # Standardized weights for the size distribution (in IC):
    wₓ, ẇₓ = compute_wₓ(lin, size_distribution, wᵣ, r, IC(r_max))

    # Pre-allocate buffers for the inner loop (sized for the largest particle).
    # m_ref in IC keeps the Mie recursion in full precision.
    m_ref = IC(nᵣ) - IC(nᵢ) * im
    y_max = maximum(x_size_param) * abs(m_ref)
    nmx_max = round(Int, max(n_max, y_max) + 51)
    an  = zeros(Complex{IC}, n_max)
    bn  = zeros(Complex{IC}, n_max)
    ȧn  = zeros(Complex{IC}, 2, n_max)
    ḃn  = zeros(Complex{IC}, 2, n_max)
    Dₙ  = zeros(Complex{IC}, nmx_max)
    Ḋₙ  = zeros(Complex{IC}, 2, nmx_max)
    n_  = IC.(2 .* collect(1:n_max) .+ 1)

    # Loop over size parameters
    for i = 1:length(x_size_param)

        n_max_i = get_n_max(x_size_param[i])

        # Zero the views that will be used
        an_v = view(an, 1:n_max_i)
        bn_v = view(bn, 1:n_max_i)
        ȧn_v = view(ȧn, :, 1:n_max_i)
        ḃn_v = view(ḃn, :, 1:n_max_i)
        fill!(an_v, 0); fill!(bn_v, 0)
        fill!(ȧn_v, 0); fill!(ḃn_v, 0)
        fill!(Dₙ, 0);   fill!(Ḋₙ, 0)
        
        # Compute aₙ,bₙ and S₁,S₂
        compute_mie_ab!(x_size_param[i], m_ref, an_v, bn_v, Dₙ, ȧn_v, ḃn_v, Ḋₙ)
        compute_mie_S₁S₂!(an_v, bn_v, 
            ȧn_v, ḃn_v, 
            leg_π, leg_τ, 
            view(S₁, :, i), view(S₂, :, i), 
            view(Ṡ₁, :, :, i), view(Ṡ₂, :, :, i))

        # Compute Extinction and scattering cross sections using pre-allocated n_
        n_v = view(n_, 1:n_max_i)
        coef = 2π / k^2
        @inbounds C_sca[i] = coef * dot(n_v, abs2.(an_v) .+ abs2.(bn_v))
        @inbounds C_ext[i] = coef * dot(n_v, real.(an_v .+ bn_v))

        @inbounds for ctr=1:2
            cs = zero(IC)
            ce = zero(IC)
            for n in 1:n_max_i
                cs += n_[n] * real(an_v[n] * conj(ȧn[ctr,n]) + ȧn[ctr,n] * conj(an_v[n]) +
                                  bn_v[n] * conj(ḃn[ctr,n]) + ḃn[ctr,n] * conj(bn_v[n]))
                ce += n_[n] * real(ȧn[ctr,n] + ḃn[ctr,n])
            end
            Ċ_sca[ctr,i] = coef * cs
            Ċ_ext[ctr,i] = coef * ce
        end

        # Compute scattering matrix components per size parameter (in-place)
        # inv_x2 in IC so all phase-matrix products stay in full precision
        inv_x2 = IC(0.5) / x_size_param[i]^2
        @inbounds for iμ in 1:n_mu
            s1 = S₁[iμ, i]; s2 = S₂[iμ, i]
            f₁₁[iμ, i] =  inv_x2 * real(abs2(s1) + abs2(s2))
            f₃₃[iμ, i] =  inv_x2 * real(s1 * conj(s2) + s2 * conj(s1))
            f₁₂[iμ, i] = -inv_x2 * real(abs2(s1) - abs2(s2))
            f₃₄[iμ, i] = -inv_x2 * imag(s1 * conj(s2) - s2 * conj(s1))

            for ctr in 1:2
                ds1 = Ṡ₁[ctr, iμ, i]; ds2 = Ṡ₂[ctr, iμ, i]
                ḟ₁₁[ctr,iμ,i] =  inv_x2 * real(ds1*conj(s1) + s1*conj(ds1) +
                                                  ds2*conj(s2) + s2*conj(ds2))
                ḟ₃₃[ctr,iμ,i] =  inv_x2 * real(ds1*conj(s2) + s1*conj(ds2) +
                                                  ds2*conj(s1) + s2*conj(ds1))
                ḟ₁₂[ctr,iμ,i] = -inv_x2 * real(ds1*conj(s1) + s1*conj(ds1) -
                                                  ds2*conj(s2) - s2*conj(ds2))
                ḟ₃₄[ctr,iμ,i] = -inv_x2 * imag(ds1*conj(s2) + s1*conj(ds2) -
                                                  ds2*conj(s1) - s2*conj(ds1))
            end
        end
    end

    # Calculate bulk scattering and extinction cross-sections
    bulk_C_sca =  sum(wₓ .* C_sca)
    bulk_C_ext =  sum(wₓ .* C_ext)

    @views for ctr=1:2
        bulk_Ċ_sca[ctr] =  sum(wₓ .* Ċ_sca[ctr,:])
        bulk_Ċ_ext[ctr] =  sum(wₓ .* Ċ_ext[ctr,:])
    end
    @views for ctr=3:4
        bulk_Ċ_sca[ctr] =  sum(ẇₓ[ctr-2, :] .* C_sca)
        bulk_Ċ_ext[ctr] =  sum(ẇₓ[ctr-2, :] .* C_ext)
    end
    
    bulk_ϖ = bulk_C_sca / bulk_C_ext
    for ctr=1:4
        bulk_ϖ̇[ctr] = (bulk_Ċ_sca[ctr] - bulk_ϖ * bulk_Ċ_ext[ctr]) / bulk_C_ext
    end
    # Compute bulk scattering 
    wr = (4π * r.^2 .*  wₓ) 
    ẇr = zeros(FT, 2, length(r))
    
    @views for ctr=1:2
        ẇr[ctr,:] = (4π * r.^2 .*  ẇₓ[ctr,:]) 
    end
    bulk_f₁₁   =  f₁₁ * wr
    bulk_f₃₃   =  f₃₃ * wr
    bulk_f₁₂   =  f₁₂ * wr
    bulk_f₃₄   =  f₃₄ * wr

    # Normalize Phase function with bulk scattering cross section.
    bulk_f₁₁ /= bulk_C_sca 
    bulk_f₃₃ /= bulk_C_sca
    bulk_f₁₂ /= bulk_C_sca
    bulk_f₃₄ /= bulk_C_sca
    @views for ctr=1:2
        bulk_ḟ₁₁[ctr,:] = (ḟ₁₁[ctr,:,:] * wr - bulk_f₁₁ * bulk_Ċ_sca[ctr]) / bulk_C_sca
        bulk_ḟ₃₃[ctr,:] = (ḟ₃₃[ctr,:,:] * wr - bulk_f₃₃ * bulk_Ċ_sca[ctr]) / bulk_C_sca
        bulk_ḟ₁₂[ctr,:] = (ḟ₁₂[ctr,:,:] * wr - bulk_f₁₂ * bulk_Ċ_sca[ctr]) / bulk_C_sca
        bulk_ḟ₃₄[ctr,:] = (ḟ₃₄[ctr,:,:] * wr - bulk_f₃₄ * bulk_Ċ_sca[ctr]) / bulk_C_sca
    end
    @views for ctr=3:4
        bulk_ḟ₁₁[ctr,:] = (f₁₁ * ẇr[ctr-2,:] - bulk_f₁₁ * bulk_Ċ_sca[ctr]) / bulk_C_sca
        bulk_ḟ₃₃[ctr,:] = (f₃₃ * ẇr[ctr-2,:] - bulk_f₃₃ * bulk_Ċ_sca[ctr]) / bulk_C_sca
        bulk_ḟ₁₂[ctr,:] = (f₁₂ * ẇr[ctr-2,:] - bulk_f₁₂ * bulk_Ċ_sca[ctr]) / bulk_C_sca
        bulk_ḟ₃₄[ctr,:] = (f₃₄ * ẇr[ctr-2,:] - bulk_f₃₄ * bulk_Ċ_sca[ctr]) / bulk_C_sca
    end
    # Range of l-values
    l_max = length(μ);

    # Get legendre polynomials for l-max 
    P, P², R², T² = compute_legendre_poly(μ, l_max)

    # Compute Greek coefficients in IC (Float64 for plain-float models).
    # All angular quadrature products are already IC; narrowing to FT2 happens below.
    α = zeros(IC, l_max)
    β = zeros(IC, l_max)
    δ = zeros(IC, l_max)
    γ = zeros(IC, l_max)
    ϵ = zeros(IC, l_max)
    ζ = zeros(IC, l_max)

    α̇ = zeros(IC, 4, l_max)
    β̇ = zeros(IC, 4, l_max)
    δ̇ = zeros(IC, 4, l_max)
    γ̇ = zeros(IC, 4, l_max)
    ϵ̇ = zeros(IC, 4, l_max)
    ζ̇ = zeros(IC, 4, l_max)
    
    # Compute Greek coefficients from bulk scattering matrix elements (spherical only here!)
    for l = 0:length(β) - 1

        # pre-factor:
        fac = l ≥ 2 ? (2l + 1) / 2 * sqrt(1 / ((l - 1) * (l) * (l + 1) * (l + 2))) : 0

        # Compute Greek coefficients 
        # Eq 17 in Sanghavi 2014

        δ[l + 1] = (2l + 1) / 2 * w_μ' * (bulk_f₃₃ .* P[:,l + 1])
        β[l + 1] = (2l + 1) / 2 * w_μ' * (bulk_f₁₁ .* P[:,l + 1])
        γ[l + 1] = fac          * w_μ' * (bulk_f₁₂ .* P²[:,l + 1])
        ϵ[l + 1] = fac          * w_μ' * (bulk_f₃₄ .* P²[:,l + 1])
        ζ[l + 1] = fac          * w_μ' * (bulk_f₃₃ .* R²[:,l + 1] + bulk_f₁₁ .* T²[:,l + 1]) 
        α[l + 1] = fac          * w_μ' * (bulk_f₁₁ .* R²[:,l + 1] + bulk_f₃₃ .* T²[:,l + 1]) 

        @views for ctr=1:4
            δ̇[ctr,l + 1] = (2l + 1) / 2 * dot(w_μ, bulk_ḟ₃₃[ctr,:] .* P[:,l + 1])
            β̇[ctr,l + 1] = (2l + 1) / 2 * dot(w_μ, bulk_ḟ₁₁[ctr,:] .* P[:,l + 1])
            γ̇[ctr,l + 1] = fac          * dot(w_μ, bulk_ḟ₁₂[ctr,:] .* P²[:,l + 1])
            ϵ̇[ctr,l + 1] = fac          * dot(w_μ, bulk_ḟ₃₄[ctr,:] .* P²[:,l + 1])
            ζ̇[ctr,l + 1] = fac          * dot(w_μ, bulk_ḟ₃₃[ctr,:] .* R²[:,l + 1] .+ bulk_ḟ₁₁[ctr,:] .* T²[:,l + 1])
            α̇[ctr,l + 1] = fac          * dot(w_μ, bulk_ḟ₁₁[ctr,:] .* R²[:,l + 1] .+ bulk_ḟ₃₃[ctr,:] .* T²[:,l + 1])
        end    
    end

    # Convert internal IC (Float64) results to the requested output type FT2.
    # Symmetric conversion: both value greek coefficients AND derivative (lin)
    # greek coefficients are narrowed to FT2. For Dual types (IC != AbstractFloat),
    # skip conversion to preserve derivative tangents.
    if FT2 <: AbstractFloat
        # Create GreekCoefs object with α, β, γ, δ, ϵ, ζ
        greek_coefs = GreekCoefs(convert.(FT2, α), 
                                 convert.(FT2, β), 
                                 convert.(FT2, γ), 
                                 convert.(FT2, δ), 
                                 convert.(FT2, ϵ), 
                                 convert.(FT2, ζ))
        lin_greek_coefs = linGreekCoefs(convert.(FT2, α̇), 
                                         convert.(FT2, β̇), 
                                         convert.(FT2, γ̇), 
                                         convert.(FT2, δ̇), 
                                         convert.(FT2, ϵ̇), 
                                         convert.(FT2, ζ̇))
        # Return the packaged AerosolOptics object
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=FT2(bulk_ϖ), k=FT2(bulk_C_ext), fᵗ=FT2(1)),
            linAerosolOptics(lin_greek_coefs=lin_greek_coefs, ω̃̇=convert.(FT2, bulk_ϖ̇), k̇=convert.(FT2, bulk_Ċ_ext), ḟᵗ=FT2(0.0) .*convert.(FT2, bulk_Ċ_ext))# zeros(FT2,4))
    else
        # Dual / non-AbstractFloat path: keep IC arrays as-is (no narrowing)
        greek_coefs = GreekCoefs(α, β, γ, δ, ϵ, ζ)
        lin_greek_coefs = linGreekCoefs(α̇, β̇, γ̇, δ̇, ϵ̇, ζ̇)
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=(bulk_C_sca / bulk_C_ext), k=(bulk_C_ext), fᵗ=one(IC)),
            linAerosolOptics(lin_greek_coefs=lin_greek_coefs, ω̃̇=bulk_ϖ̇, k̇=bulk_Ċ_ext, ḟᵗ=zero(IC) .* bulk_Ċ_ext)
    end
end

function compute_ref_aerosol_extinction(lin::LinMode, model::MieModel{FDT,FT}, FT2::Type) where {FDT <: NAI2, FT}

    # Unpack the model
    (; computation_type, aerosol, λ, polarization_type, r_max, nquad_radius) = model

    # Extract variables from aerosol struct:
    (; size_distribution, nᵣ, nᵢ) = aerosol

    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (definition)"

    # Internal computation precision: Float64 for plain floats, native for Dual.
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

    # Pre-allocate arrays in IC
    C_ext = zeros(IC, nquad_radius)
    Ċ_ext = zeros(IC, 2, nquad_radius)
    bulk_Ċ_ext = zeros(IC, 4)

    # Standardized weights for the size distribution:
    wₓ, ẇₓ = compute_wₓ(lin, size_distribution, wᵣ, r, IC(r_max))

    # Pre-allocate buffers outside the loop (sized for largest particle)
    m_ref = IC(nᵣ) - IC(nᵢ) * im
    y_max = maximum(x_size_param) * abs(m_ref)
    nmx_max = round(Int, max(n_max, y_max) + 51)
    an  = zeros(Complex{IC}, n_max)
    bn  = zeros(Complex{IC}, n_max)
    ȧn  = zeros(Complex{IC}, 2, n_max)
    ḃn  = zeros(Complex{IC}, 2, n_max)
    Dₙ  = zeros(Complex{IC}, nmx_max)
    Ḋₙ  = zeros(Complex{IC}, 2, nmx_max)
    n_  = IC.(2 .* collect(1:n_max) .+ 1)

    # Loop over size parameters
    for i = 1:length(x_size_param)

        n_max_i = get_n_max(x_size_param[i])

        # Zero the views
        an_v = view(an, 1:n_max_i)
        bn_v = view(bn, 1:n_max_i)
        ȧn_v = view(ȧn, :, 1:n_max_i)
        ḃn_v = view(ḃn, :, 1:n_max_i)
        fill!(an_v, 0); fill!(bn_v, 0)
        fill!(ȧn_v, 0); fill!(ḃn_v, 0)
        fill!(Dₙ, 0); fill!(Ḋₙ, 0)

        compute_mie_ab!(x_size_param[i], m_ref, an_v, bn_v, Dₙ, ȧn_v, ḃn_v, Ḋₙ)

        # Compute Extinction cross sections
        n_v = view(n_, 1:n_max_i)
        coef = IC(2π) / k^2
        @inbounds C_ext[i] = coef * dot(n_v, real.(an_v .+ bn_v))
        @inbounds for ctr=1:2
            ce = zero(IC)
            for n in 1:n_max_i
                ce += n_[n] * real(ȧn[ctr,n] + ḃn[ctr,n])
            end
            Ċ_ext[ctr,i] = coef * ce
        end
    end

    # Calculate bulk extinction cross-section (IC precision)
    bulk_C_ext = sum(wₓ .* C_ext)
    for ctr=1:2
        bulk_Ċ_ext[ctr] = sum(wₓ .* Ċ_ext[ctr,:])
    end
    for ctr=3:4
        bulk_Ċ_ext[ctr] = sum(ẇₓ[ctr-2,:] .* C_ext)
    end

    # Return bulk extinction (converted to FT2) and derivative vector (converted to FT2).
    # Symmetric: both value and derivative are narrowed to the same output type.
    if FT2 <: AbstractFloat
        return FT2(bulk_C_ext), convert.(FT2, bulk_Ċ_ext)
    else
        return bulk_C_ext, bulk_Ċ_ext
    end
end
