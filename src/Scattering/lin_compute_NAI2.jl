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
#TODO: define linAerosolOptics 
function compute_aerosol_optical_properties(model::MieModel{FDT}, FT2::Type=Float64) where FDT <: NAI2

    # Unpack the model
    @unpack aerosol, λ, r_max, nquad_radius = model
    @show nquad_radius, λ
    # Extract variables from aerosol struct:
    @unpack size_distribution, nᵣ, nᵢ = aerosol
    #@show typeof(size_distribution.σ), typeof(nᵣ)
    # Imaginary part of the refractive index must be ≥ 0
    @assert nᵢ ≥ 0

    # Get the refractive index's real part type
    #@show size_distribution.σ  
    # TODO: This is still very clumsy, the FT conversions are not good here.
    FT = eltype(nᵣ);

    #@show FT, ForwardDiff.valtype(size_distribution.σ)
    vFT = ForwardDiff.valtype(nᵣ)
    #@assert FT == Float64 "Aerosol computations require 64bit"
    # Get radius quadrature points and weights (for mean, thus normalized):
    # 
    
    # Just sample from 0.25%ile to 99.75%ile:
    #start,stop = quantile(size_distribution,[vFT(0.0001),vFT(0.99999999)])
    r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true) 
    #r, wᵣ = gauleg(nquad_radius, 0.0, min(stop,r_max) ; norm=true) 
    #@show start, stop
    #@show typeof(r), typeof(wᵣ)
    #@show typeof(convert.(FT, r))
    r  = convert.(FT, r)
    wᵣ = convert.(FT, wᵣ)
    # Wavenumber
    k = vFT(2π) / λ  

    # Size parameter
    x_size_param = k * r # (2πr/λ)

    # Compute n_max for largest size:
    n_max = get_n_max(maximum(x_size_param))

    # Determine max number of Gaussian quadrature points for angle dependence of 
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

    dS₁    = zeros(Complex{FT}, n_mu, nquad_radius)
    dS₂    = zeros(Complex{FT}, n_mu, nquad_radius)
    df₁₁   = zeros(FT, 2, n_mu, nquad_radius)
    df₃₃   = zeros(FT, 2, n_mu, nquad_radius)
    df₁₂   = zeros(FT, 2, n_mu, nquad_radius)
    df₃₄   = zeros(FT, 2, n_mu, nquad_radius)
    dC_ext = zeros(FT, 2, nquad_radius) #derivatives with respect to nᵣ, nᵢ
    dC_sca = zeros(FT, 2, nquad_radius) #derivatives with respect to nᵣ, nᵢ

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    dwₓ = zeros(FT, 2, length(wₓ))
    dwₓ[1,:] = wₓ.*(log.(r)-size_distribution.μ)/(size_distribution.σ)^2
    dwₓ[2,:] = wₓ.*(((log.(r)-size_distribution.μ)/size_distribution.σ)^2 .- 1)/size_distribution.σ 
    # Loop over size parameters
    @showprogress 1 "Computing PhaseFunctions Siewert NAI-2 style ..." for i = 1:length(x_size_param)

        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = get_n_max(x_size_param[i])

        # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
        an = (zeros(Complex{FT}, n_max))
        bn = (zeros(Complex{FT}, n_max))
        dan = (zeros(Complex{FT}, n_max))
        dbn = (zeros(Complex{FT}, n_max))
        
        # Weighting for sums of 2n+1
        n_ = collect(FT, 1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dₙ  :
        y = x_size_param[i] * (aerosol.nᵣ - aerosol.nᵢ);
        nmx = round(Int, max(n_max, abs(y)) + 51)
        Dₙ  = zeros(Complex{FT}, nmx)
        dDₙ  = zeros(Complex{FT}, nmx)

        # Compute aₙ,bₙ and S₁,S₂
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ - aerosol.nᵢ * im, an, bn, Dₙ, dan, dbn, dDₙ)
        compute_mie_S₁S₂!(an, bn, dan, dbn, leg_π, leg_τ, view(S₁, :, i), view(S₂, :, i), view(dS₁, :, i), view(dS₂, :, i))
        
        # Compute Extinction and scattering cross sections: 
        C_sca[i] = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))

        dC_sca[1,i] = 2π / k^2 * 2*(n_' * (real(dan).*real(an) +
                imag(dan).*imag(an) + 
                real(dbn).*real(bn) +
                imag(dbn).*imag(bn)))
        dC_sca[2,i] = 2π / k^2 * 2*(n_' * (real(an).*imag(dan) -
                real(dan).*imag(an) + 
                real(bn).*imag(dbn) -
                reak(dbn).*imag(bn)))

        dC_ext[1,i] = 2π / k^2 * (n_' * real(dan + dbn))
        dC_ext[2,i] = 2π / k^2 * (n_' * imag(dan + dbn))

        # Compute scattering matrix components per size parameter:
        f₁₁[:,i] =  vFT(0.5) / x_size_param[i]^2   * (abs2.(S₁[:,i]) + abs2.(S₂[:,i]));
        f₃₃[:,i] =  vFT(0.5) / x_size_param[i]^2   * (S₁[:,i] .* conj(S₂[:,i]) + S₂[:,i] .* conj(S₁[:,i]));
        f₁₂[:,i] = -vFT(0.5) / x_size_param[i]^2   * (abs2.(S₁[:,i]) - abs2.(S₂[:,i]));
        f₃₄[:,i] = -vFT(0.5) / x_size_param[i]^2   * imag(S₁[:,i] .* conj(S₂[:,i]) - S₂[:,i] .* conj(S₁[:,i]));

        df₁₁[1,:,i] =  vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*real(dS₁[:,i]) + 
                real(S₂[:,i]).*real(dS₂[:,i]) + 
                imag(S₁[:,i]).*imag(dS₁[:,i]) + 
                imag(S₂[:,i]).*imag(dS₂[:,i]));
        df₁₁[2,:,i] = vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*imag(dS₁[:,i]) - 
               real(dS₁[:,i]).*imag(S₁[:,i]) + 
               real(S₂[:,i]).*imag(dS₂[:,i]) - 
               real(dS₂[:,i]).*imag(S₂[:,i]));
        
        df₃₃[1,:,i] =  vFT(0.5) / x_size_param[i]^2   * 
            2*( real(dS₁[:,i]) .* real(S₂[:,i]) + 
                real(S₁[:,i]) .* real(dS₂[:,i]) +
                imag(dS₁[:,i]) .* imag(S₂[:,i]) + 
                imag(S₁[:,i]) .* imag(dS₂[:,i]));       
        df₃₃[2,:,i] =  vFT(0.5) / x_size_param[i]^2   * 
                2*( real(S₁[:,i]) .* imag(dS₂[:,i]) - 
                    real(dS₁[:,i]) .* imag(S₂[:,i]) -
                    imag(S₁[:,i]) .* real(dS₂[:,i]) + 
                    imag(dS₁[:,i]) .* real(S₂[:,i]));
            
        
        df₁₂[1,:,i] = -vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*real(dS₁[:,i]) - 
                real(S₂[:,i]).*real(dS₂[:,i]) + 
                imag(S₁[:,i]).*imag(dS₁[:,i]) - 
                imag(S₂[:,i]).*imag(dS₂[:,i]));
        df₁₂[2,:,i] = -vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*imag(dS₁[:,i]) - 
            real(dS₁[:,i]).*imag(S₁[:,i]) - 
            real(S₂[:,i]).*imag(dS₂[:,i]) + 
            real(dS₂[:,i]).*imag(S₂[:,i]));
    

        df₃₄[1,:,i] = -vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₂[:,i]).*imag(dS₁[:,i]) + 
                real(dS₂[:,i]).*imag(S₁[:,i]) -
                imag(S₂[:,i]).*real(dS₁[:,i]) -
                imag(dS₂[:,i]).*real(S₁[:,i]))
        df₃₄[2,:,i] = -vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*real(dS₂[:,i]) + 
               imag(S₁[:,i]).*imag(dS₂[:,i]) -
               real(dS₁[:,i]).*real(S₂[:,i]) -
               imag(dS₁[:,i]).*imag(S₂[:,i]))
    end

    # Calculate bulk scattering and extinction coeffitientcs
    bulk_C_sca =  sum(wₓ .* C_sca)
    bulk_C_ext =  sum(wₓ .* C_ext)
    d_bulk_C_sca   = zeros(FT, 4)
    d_bulk_C_ext   = zeros(FT, 4)
    dω̃ = zeros(FT, 4)
    
    for i=1:2
        d_bulk_C_sca[i] =  sum(wₓ .* dC_sca[i,:])
        d_bulk_C_ext[i] =  sum(wₓ .* dC_ext][i,:])
    end
    for i=3:4
        ii=i-2
        d_bulk_C_sca[i] =  sum(dwₓ[ii,:] .* C_sca)
        d_bulk_C_ext[i] =  sum(dwₓ[ii,:] .* C_ext])
    end
    ω̃ = bulk_C_sca/bulk_C_ext
    for i=1:4
        dω̃[i] = d_bulk_C_sca[i]/bulk_C_ext - bulk_C_sca*d_bulk_C_ext[i]/bulk_C_ext^2
    end
    
    # Compute bulk scattering 
    wr = (4π * r.^2 .*  wₓ) 
    dwr = zeros(FT, 4, length(wr))
    for i=1:4
        dwr[i,:] = (4π * r.^2 .*  dwₓ[i,:]) 
    end
    bulk_f₁₁   =  f₁₁ * wr
    bulk_f₃₃   =  f₃₃ * wr
    bulk_f₁₂   =  f₁₂ * wr
    bulk_f₃₄   =  f₃₄ * wr
    d_bulk_f₁₁ = zeros(FT, 4, n_mu)
    d_bulk_f₁₂ = zeros(FT, 4, n_mu)
    d_bulk_f₃₃ = zeros(FT, 4, n_mu)
    d_bulk_f₃₄ = zeros(FT, 4, n_mu)

    for i=1:2
        for j=1:n_mu
            d_bulk_f₁₁[i,j]   =  df₁₁[i,j,:] * wr
            d_bulk_f₃₃[i,j]   =  df₃₃[i,j,:] * wr
            d_bulk_f₁₂[i,j]   =  df₁₂[i,j,:] * wr
            d_bulk_f₃₄[i,j]   =  df₃₄[i,j,:] * wr
        end
    end
    for i=3:4
        ii=i-2
        for j=1:n_mu
            d_bulk_f₁₁[i,j]   =  f₁₁[j,:] * dwr[ii,:]
            d_bulk_f₃₃[i,j]   =  f₃₃[j,:] * dwr[ii,:]
            d_bulk_f₁₂[i,j]   =  f₁₂[j,:] * dwr[ii,:]
            d_bulk_f₃₄[i,j]   =  f₃₄[j,:] * dwr[ii,:]
        end
    end

    # Normalize Phase function with bulk scattering cross section.
    bulk_f₁₁ /= bulk_C_sca 
    bulk_f₃₃ /= bulk_C_sca
    bulk_f₁₂ /= bulk_C_sca
    bulk_f₃₄ /= bulk_C_sca

    d_bulk_f₁₁ /= bulk_C_sca 
    d_bulk_f₃₃ /= bulk_C_sca
    d_bulk_f₁₂ /= bulk_C_sca
    d_bulk_f₃₄ /= bulk_C_sca
    
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

    dα = zeros(FT, 4, l_max)
    dβ = zeros(FT, 4, l_max)
    dδ = zeros(FT, 4, l_max)
    dγ = zeros(FT, 4, l_max)
    dϵ = zeros(FT, 4, l_max)
    dζ = zeros(FT, 4, l_max)
    
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
        
        for i=1:4
            dδ[i,l + 1] = (2l + 1) / 2 * w_μ' * (d_bulk_f₃₃[i,:] .* P[:,l + 1])
            dβ[i,l + 1] = (2l + 1) / 2 * w_μ' * (d_bulk_f₁₁[i,:] .* P[:,l + 1])
            dγ[i,l + 1] = fac      * w_μ' * (d_bulk_f₁₂[i,:] .* P²[:,l + 1])
            dϵ[i,l + 1] = fac      * w_μ' * (d_bulk_f₃₄[i,:] .* P²[:,l + 1])
            dζ[i,l + 1] = fac      * w_μ' * (d_bulk_f₃₃[i,:] .* R²[:,l + 1] + d_bulk_f₁₁[i,:] .* T²[:,l + 1]) 
            dα[i,l + 1] = fac      * w_μ' * (d_bulk_f₁₁[i,:] .* R²[:,l + 1] + d_bulk_f₃₃[i,:] .* T²[:,l + 1]) 
        end
    end

    # Check whether this is a Dual number (if so, don't do any conversions)
    # TODO: Equally clumsy, needs to be fixed.
    if FT <: AbstractFloat
        #@show "Convert greek", FT2
        # Create GreekCoefs object with α, β, γ, δ, ϵ, ζ
        greek_coefs = GreekCoefs(convert.(FT2, α), 
                                 convert.(FT2, β), 
                                 convert.(FT2, γ), 
                                 convert.(FT2, δ), 
                                 convert.(FT2, ϵ), 
                                 convert.(FT2, ζ))
        d_greek_coefs=[];
        for i=1:4
            tmp_greek_coefs = GreekCoefs(convert.(FT2, dα[i,:]), 
                                        convert.(FT2, dβ[i,:]), 
                                        convert.(FT2, dγ[i,:]), 
                                        convert.(FT2, dδ[i,:]), 
                                        convert.(FT2, dϵ[i,:]), 
                                        convert.(FT2, dζ[i,:]))
            push!(d_greek_coefs, tmp_greek_coefs);
        end
        #@show typeof(convert.(FT2, β)), typeof(greek_coefs)
        # Return the packaged AerosolOptics object
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=FT2(ω̃), k=FT2(bulk_C_ext), fᵗ=FT2(1)), linAerosolOptics(d_greek_coefs=d_greek_coefs, dω̃=FT2(dω̃), dk=FT2(d_bulk_C_ext), dfᵗ=FT2(1))

    else
        greek_coefs = GreekCoefs(α,β,γ,δ,ϵ,ζ)
        d_greek_coefs=[];
        for i=1:4
            tmp_greek_coefs = GreekCoefs(dα[i,:], 
                                        dβ[i,:], 
                                        dγ[i,:], 
                                        dδ[i,:], 
                                        dϵ[i,:], 
                                        dζ[i,:])
            push!(d_greek_coefs, tmp_greek_coefs);
        end
        return AerosolOptics(greek_coefs=greek_coefs, ω̃=(bulk_C_sca / bulk_C_ext), k=(bulk_C_ext), fᵗ=FT(1)), , linAerosolOptics(d_greek_coefs=d_greek_coefs, dω̃=dω̃, dk=d_bulk_C_ext, dfᵗ=FT(1))
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
    # start,stop = quantile(size_distribution,[0.0001, 0.99999999])
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
    #leg_π, leg_τ = compute_mie_π_τ(μ, n_max)

    # Pre-allocate arrays:
    C_ext = zeros(FT, nquad_radius)
    dC_ext = zeros(FT, 2, nquad_radius)
    d_bulk_C_ext   = zeros(FT, 4)
    #C_sca = zeros(FT, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    dwₓ = zeros(FT, 2, length(wₓ))
    dwₓ[1,:] = wₓ.*(log.(r)-size_distribution.μ)/(size_distribution.σ)^2
    dwₓ[2,:] = wₓ.*(((log.(r)-size_distribution.μ)/size_distribution.σ)^2 .- 1)/size_distribution.σ 
    
    # Loop over size parameters
    @showprogress 1 "Computing extinction XS at reference wavelength ..." for i = 1:length(x_size_param)

        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = get_n_max(x_size_param[i])

        # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
        an = zeros(Complex{FT}, n_max)
        bn = zeros(Complex{FT}, n_max)
        dan = zeros(Complex{FT}, n_max)
        dbn = zeros(Complex{FT}, n_max)

        # Weighting for sums of 2n+1
        n_ = collect(FT, 1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dn:
        y = x_size_param[i] * (aerosol.nᵣ - aerosol.nᵢ);
        nmx = round(Int, max(n_max, abs(y)) + 51)
        Dn = zeros(Complex{FT}, nmx)
        dDn = zeros(Complex{FT}, nmx)

        # Compute an,bn and S₁,S₂
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ - aerosol.nᵢ * im, an, bn, Dn, dan, dbn, dDn)
        
        # Compute Extinction and scattering cross sections: 
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))
        dC_ext[1,i] = 2π / k^2 * (n_' * real(dan + dbn))
        dC_ext[2,i] = 2π / k^2 * (n_' * imag(dan + dbn))    
    end

    # Calculate bulk extinction coeffitient
    bulk_C_ext =  sum(wₓ .* C_ext)
    for i=1:2
        d_bulk_C_ext[i] =  sum(wₓ .* dC_ext[i,:])
    end
    for i=3:4
        ii=i-2
        d_bulk_C_ext[i] =  sum(dwₓ[ii,:] .* C_ext)
    end
    # Return the bulk extinction coeffitient
    return bulk_C_ext, d_bulk_C_ext
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

    dS₁    = zeros(Complex{FT}, n_mu, nquad_radius)
    dS₂    = zeros(Complex{FT}, n_mu, nquad_radius)
    df₁₁   = zeros(FT, 2, n_mu, nquad_radius)
    dC_ext = zeros(FT, 2, nquad_radius)
    dC_sca = zeros(FT, 2, nquad_radius)

    # Standardized weights for the size distribution:
    wₓ = compute_wₓ(size_distribution, wᵣ, r, r_max) 
    dwₓ = zeros(FT, 2, length(wₓ))
    dwₓ[1,:] = wₓ.*(log.(r)-size_distribution.μ)/(size_distribution.σ)^2
    dwₓ[2,:] = wₓ.*(((log.(r)-size_distribution.μ)/size_distribution.σ)^2 .- 1)/size_distribution.σ 
    
    # Loop over size parameters
    @showprogress 1 "Computing PhaseFunctions Siewert NAI-2 style ..." for i = 1:length(x_size_param)

        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = get_n_max(x_size_param[i])

        # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
        an = (zeros(Complex{FT}, n_max))
        bn = (zeros(Complex{FT}, n_max))
        dan = (zeros(Complex{FT}, n_max))
        dbn = (zeros(Complex{FT}, n_max))
       
        # Weighting for sums of 2n+1
        n_ = collect(FT, 1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dn:
        y = x_size_param[i] * (aerosol.nᵣ - aerosol.nᵢ);
        nmx = round(Int, max(n_max, abs(y)) + 51)
        Dn = zeros(Complex{FT}, nmx)
        dDn = zeros(Complex{FT}, nmx)
        
        # Compute an,bn and S₁,S₂
        compute_mie_ab!(x_size_param[i], aerosol.nᵣ + aerosol.nᵢ * im, an, bn, Dn, dan, dbn, dDn)
        compute_mie_S₁S₂!(an, bn, dan, dbn, leg_π, leg_τ, view(S₁, :, i), view(S₂, :, i), view(dS₁, :, i), view(dS₂, :, i))
        
        # Compute Extinction and scattering cross sections: 
        C_sca[i] = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
        C_ext[i] = 2π / k^2 * (n_' * real(an + bn))
        dC_sca[1,i] = 2π / k^2 * 2*(n_' * (real(dan).*real(an) +
                imag(dan).*imag(an) + 
                real(dbn).*real(bn) +
                imag(dbn).*imag(bn)))
        dC_sca[2,i] = 2π / k^2 * 2*(n_' * (real(an).*imag(dan) -
                real(dan).*imag(an) + 
                real(bn).*imag(dbn) -
                reak(dbn).*imag(bn)))

        dC_ext[1,i] = 2π / k^2 * (n_' * real(dan + dbn))
        dC_ext[2,i] = 2π / k^2 * (n_' * imag(dan + dbn))
        # @show r[i], x_size_param[i], C_ext[i], C_ext[i]/(4π*r[i]^2), C_ext[i]*1e-8
        # Compute scattering matrix components per size parameter (might change column/row ordering):
        f₁₁[:,i] =  0.5 / x_size_param[i]^2  * real(abs2.(S₁[:,i]) + abs2.(S₂[:,i]));
        df₁₁[1,:,i] =  vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*real(dS₁[:,i]) + 
                real(S₂[:,i]).*real(dS₂[:,i]) + 
                imag(S₁[:,i]).*imag(dS₁[:,i]) + 
                imag(S₂[:,i]).*imag(dS₂[:,i]));
        df₁₁[2,:,i] = vFT(0.5) / x_size_param[i]^2   * 
            2*(real(S₁[:,i]).*imag(dS₁[:,i]) - 
                real(dS₁[:,i]).*imag(S₁[:,i]) + 
                real(S₂[:,i]).*imag(dS₂[:,i]) - 
                real(dS₂[:,i]).*imag(S₂[:,i]));
    end

    # Calculate bulk scattering and extinction coeffitientcs
    bulk_C_sca =  sum(wₓ .* C_sca)
    bulk_C_ext =  sum(wₓ .* C_ext)
    d_bulk_C_sca   = zeros(FT, 4)
    d_bulk_C_ext   = zeros(FT, 4)
    #dω̃ = zeros(FT, 4)
    
    for i=1:2
        d_bulk_C_sca[i] =  sum(wₓ .* dC_sca[i,:])
        d_bulk_C_ext[i] =  sum(wₓ .* dC_ext][i,:])
    end
    for i=3:4
        ii=i-2
        d_bulk_C_sca[i] =  sum(dwₓ[ii,:] .* C_sca)
        d_bulk_C_ext[i] =  sum(dwₓ[ii,:] .* C_ext])
    end
    #ω̃ = bulk_C_sca/bulk_C_ext
    #for i=1:4
    #    dω̃[i] = d_bulk_C_sca[i]/bulk_C_ext - bulk_C_sca*d_bulk_C_ext[i]/bulk_C_ext^2
    #end
    
    
    # Compute bulk scattering 
    wr = (4π * r.^2 .*  wₓ) 
    dwr = zeros(FT, 4, length(wr))
    for i=1:4
        dwr[i,:] = (4π * r.^2 .*  dwₓ[i,:]) 
    end

    bulk_f₁₁   =  f₁₁ * wr
    d_bulk_f₁₁ = zeros(FT, 4, n_mu)
    for i=1:2
        for j=1:n_mu
            d_bulk_f₁₁[i,j]   =  df₁₁[i,j,:] * wr
        end
    end
    for i=3:4
        ii=i-2
        for j=1:n_mu
            d_bulk_f₁₁[i,j]   =  f₁₁[j,:] * dwr[ii,:]
        end
    end
 
    # Normalize Phase function with bulk scattering cross section.
    bulk_f₁₁ /= bulk_C_sca 
    d_bulk_f₁₁ /= bulk_C_sca 

    # Assymetry factor g = 0.5\int_{-1}^1 μ P(μ) dμ
    g = 1/2 * w_μ'*(μ .*  bulk_f₁₁ )
    dg = zeros(FT,4)
    for i=1:4
        dg[i] = 1/2 * w_μ'*(μ .*  bulk_f₁₁[i,:] )
    end
    return μ, w_μ, bulk_f₁₁, bulk_C_ext, bulk_C_sca, g,
                d_bulk_f₁₁, d_bulk_C_ext, d_bulk_C_sca, dg
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
    dS₁    = zeros(Complex{FT}, n_mu)
    dS₂    = zeros(Complex{FT}, n_mu)

    # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
    an = (zeros(Complex{FT}, n_max))
    bn = (zeros(Complex{FT}, n_max))
    dan = (zeros(Complex{FT}, n_max))
    dbn = (zeros(Complex{FT}, n_max))

    # Weighting for sums of 2n+1
    n_ = collect(FT, 1:n_max);
    n_ = 2n_ .+ 1

    # Pre-allocate Dn:
    y = size_param * (nᵣ - nᵢ);
    nmx = round(Int, max(n_max, abs(y)) + 51)
    Dn = zeros(Complex{FT}, nmx)
    dDn = zeros(Complex{FT}, nmx)
    
    # Compute an,bn and S₁,S₂
    compute_mie_ab!(size_param, nᵣ + nᵢ * im, an, bn, Dn)
    compute_mie_ab!(size_param, nᵣ + nᵢ * im, an, bn, Dn, dan, dbn, dDn)
    compute_mie_S₁S₂!(an, bn, leg_π, leg_τ, S₁, S₂)
    compute_mie_S₁S₂!(an, bn, dan, dbn, leg_π, leg_τ, S₁, S₂, dS₁, dS₂)
        
    # Compute Extinction and scattering cross sections: 
    C_sca = 2π / k^2 * (n_' * (abs2.(an) + abs2.(bn)))
    C_ext = 2π / k^2 * (n_' * real(an + bn))
    dC_sca = zeros(FT, 2)
    dC_ext = zeros(FT, 2)
    df₁₁ = zeros(FT, 2, n_mu)
    dg = zeros(FT, 2)

    dC_sca[1,i] = 2π / k^2 * 2*(n_' * (real(dan).*real(an) +
                imag(dan).*imag(an) + 
                real(dbn).*real(bn) +
                imag(dbn).*imag(bn)))
    dC_sca[2,i] = 2π / k^2 * 2*(n_' * (real(an).*imag(dan) -
                real(dan).*imag(an) + 
                real(bn).*imag(dbn) -
                reak(dbn).*imag(bn)))
    dC_ext[1,i] = 2π / k^2 * (n_' * real(dan + dbn))
    dC_ext[2,i] = 2π / k^2 * (n_' * imag(dan + dbn))

    # Compute scattering matrix components per size parameter (might change column/row ordering):
    f₁₁ =  0.5 / size_param^2  * real(abs2.(S₁) + abs2.(S₂));
    f₁₁ *= 4π * r.^2
    df₁₁[1,:] =  vFT(0.5) / size_param^2   * (4π * r.^2) *
        2*(real.(S₁).*real.(dS₁) + 
        real.(S₂).*real.(dS₂) + 
        imag.(S₁).*imag.(dS₁) + 
        imag.(S₂).*imag.(dS₂));
    df₁₁[2,:] = vFT(0.5) / x_size_param[i]^2   * (4π * r.^2) *
        2*(real.(S₁).*imag.(dS₁) - 
        real.(dS₁).*imag.(S₁) + 
        real.(S₂).*imag.(dS₂) - 
        real.(dS₂).*imag.(S₂));

    f₁₁ /= C_sca
    df₁₁ /= C_sca

    # Asymmetry factor g
    g = 1/2 * w_μ'*(μ .*  f₁₁ )
    for i=1:2
        dg[i] = 1/2 * w_μ'*(μ .*  df₁₁[i,:] )
    end
    return μ, w_μ, f₁₁, C_ext, C_sca, g, df₁₁, dC_ext, dC_sca, dg
end