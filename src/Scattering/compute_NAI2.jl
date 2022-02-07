#=
 
This file specifies how to compute aerosol optical properties using the Siewert-NAI2 method
 
=#

"""
    $(FUNCTIONNAME)(model::MieModel{FDT}) where FDT<:NAI2

Reference: Suniti Sanghavi 2014, https://doi.org/10.1016/j.jqsrt.2013.12.015

Compute the aerosol optical properties using the Siewert-NAI2 method
Input: MieModel, holding all computation and aerosol properties 
Output: AerosolOptics, holding all Greek coefficients and Cross-Sectional information
"""
function compute_aerosol_optical_properties(model::MieModel{FDT}, FT2::Type=Float64) where FDT <: NAI2

    # Unpack the model
    @unpack computation_type, aerosol, λ, polarization_type, truncation_type, r_max, nquad_radius, wigner_A, wigner_B = model

    # Extract variables from aerosol struct:
    @unpack size_distribution, nᵣ, nᵢ = aerosol
    
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0

    # Get the refractive index's real part type
    FT = eltype(nᵣ);
    # @assert FT == Float64 "Aerosol computations require 64bit"
    # Get radius quadrature points and weights (for mean, thus normalized):
    # 
    
    # Just sample from 0.25%ile to 99.75%ile:
    start,stop = quantile(size_distribution,[0.0025,0.9975])
    #r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true) 
    r, wᵣ = gauleg(nquad_radius, start, min(stop,r_max) ; norm=true) 
    
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
    f₃₃   = zeros(FT, n_mu, nquad_radius)
    f₁₂   = zeros(FT, n_mu, nquad_radius)
    f₃₄   = zeros(FT, n_mu, nquad_radius)
    C_ext = zeros(FT, nquad_radius)
    C_sca = zeros(FT, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    
    # Loop over size parameters
    @showprogress 1 "Computing PhaseFunctions Siewert NAI-2 style ..." for i = 1:length(x_size_param)

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
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ + aerosol.nᵢ * im, an, bn, Dₙ )
        compute_mie_S₁S₂!(an, bn, leg_π, leg_τ, view(S₁, :, i), view(S₂, :, i))
        
        # Compute Extinction and scattering cross sections: 
        C_sca[i] = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))
       
        # Compute scattering matrix components per size parameter:
        f₁₁[:,i] =  0.5 / x_size_param[i]^2  * real(abs2.(S₁[:,i]) + abs2.(S₂[:,i]));
        f₃₃[:,i] =  0.5 / x_size_param[i]^2  * real(S₁[:,i] .* conj(S₂[:,i]) + S₂[:,i] .* conj(S₁[:,i]));
        f₁₂[:,i] = -0.5 / x_size_param[i]^2  * real(abs2.(S₁[:,i]) - abs2.(S₂[:,i]));
        f₃₄[:,i] = -0.5 / x_size_param[i]^2  * imag(S₁[:,i] .* conj(S₂[:,i]) - S₂[:,i] .* conj(S₁[:,i]));

    end

    # Calculate bulk scattering and extinction coeffitientcs
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

    # Compute Greek coefficients:
    α = zeros(FT, l_max)
    β = zeros(FT, l_max)
    δ = zeros(FT, l_max)
    γ = zeros(FT, l_max)
    ϵ = zeros(FT, l_max)
    ζ = zeros(FT, l_max)
    
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

    # Check whether this is a Dual number (if so, don't do any conversions)
    if FT2 <: AbstractFloat
        # Create GreekCoefs object with α, β, γ, δ, ϵ, ζ
        greek_coefs = GreekCoefs(convert.(FT2, α), 
                                 convert.(FT2, β), 
                                 convert.(FT2, γ), 
                                 convert.(FT2, δ), 
                                 convert.(FT2, ϵ), 
                                 convert.(FT2, ζ))
        # Return the packaged AerosolOptics object
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=FT2(bulk_C_sca / bulk_C_ext), k=FT2(bulk_C_ext), fᵗ=FT2(1))

    else
        greek_coefs = GreekCoefs(α, β, γ,δ,ϵ,ζ)
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=(bulk_C_sca / bulk_C_ext), k=(bulk_C_ext), fᵗ=FT(1))
    end
end

function compute_ref_aerosol_extinction(model::MieModel{FDT}, FT2::Type=Float64) where FDT <: NAI2

    # Unpack the model
    @unpack computation_type, aerosol, λ, polarization_type, r_max, nquad_radius = model

    # Extract variables from aerosol struct:
    @unpack size_distribution, nᵣ, nᵢ = aerosol
    
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (definition)"

    # Get the refractive index's real part type
    FT = eltype(nᵣ);
    #@show FT
    #@assert FT == Float64 "Aerosol computations require 64bit"
    # Get radius quadrature points and weights (for mean, thus normalized):
    r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true) 
    
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
    C_ext = zeros(FT, nquad_radius)
    C_sca = zeros(FT, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    
    # Loop over size parameters
    @showprogress 1 "Computing extinction XS at reference wavelength ..." for i = 1:length(x_size_param)

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
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ + aerosol.nᵢ * im, an, bn, Dn)
        
        # Compute Extinction and scattering cross sections: 
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))
    end

    # Calculate bulk extinction coeffitient
    bulk_C_ext =  sum(wₓ .* C_ext)
    
    # Return the bulk extinction coeffitient
    return bulk_C_ext
end

"""
    $(FUNCTIONNAME)(aerosol::Aerosol, λ)

Compute phase function from aerosol distribution with log-normal mean μ [µm] and σ
Output: μ, w_μ, P, C_ext, C_sca, g
"""
function phase_function(aerosol::Aerosol, λ, r_max, nquad_radius) 

    # Extract variables from aerosol struct:
    @unpack size_distribution, nᵣ, nᵢ = aerosol
    
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0 "Imaginary part of the refractive index must be ≥ 0 (definition)"

    # Get the refractive index's real part type
    FT = eltype(nᵣ);
    # @assert FT == Float64 "Aerosol computations require 64bit"
    # Get radius quadrature points and weights (for mean, thus normalized):
    r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true) 
    
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
    @showprogress 1 "Computing PhaseFunctions Siewert NAI-2 style ..." for i = 1:length(x_size_param)

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
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ + aerosol.nᵢ * im, an, bn, Dn)
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
    $(FUNCTIONNAME)(r::FT, λ::FT, nᵣ::FT, nᵢ::FT)

Compute phase function from mono-modal aerosol with radius `r` at wavelength `λ`, both in `[μm]`
Output: μ, w_μ, P, C_ext, C_sca, g
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
    compute_mie_ab!(size_param, nᵣ + nᵢ * im, an, bn, Dn)
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