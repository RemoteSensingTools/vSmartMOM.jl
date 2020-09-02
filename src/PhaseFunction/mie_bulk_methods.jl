


# ToDo: Enable arrays of aerosols (for μ̄, σ, nᵣ, nᵢ)
#function calc_aer_opt_prop(mod::NAI2, r, λ::Number, μ̄::Number, σ::Number, nᵣ, nᵢ; nquad_radius=2500, r_max = 30.0)
function calc_aer_opt_prop(mod::NAI2, aero::AbstractAerosolType, λ::Number)
    # Extract variables from struct:
    @unpack nquad_radius, nᵣ, nᵢ,r_max =  aero

    FT = eltype(nᵣ);
    # Get radius quadrature points and weights (for mean, thus normalized):
    r, wᵣ = gauleg(nquad_radius, 0.0, r_max ; norm=true) 
    
    # Size parameter
    x_sizeParam = 2π * r/λ
    # Compute Nmax for largest size:
    n_max = get_n_max(maximum(x_sizeParam))
    # Determine max amount of Gaussian quadrature points for angle dependence of phas functions:
    n_mu = 2n_max-1;
    # Obtain Gauss-Legendre quadrature points and weights for phase function
    μ, w_μ = gausslegendre( n_mu )
    # Compute π and τ functions
    leg_π, leg_τ = compute_mie_π_τ(μ, n_max)
    # Wavenumber:
    k = 2π/λ

    # Pre-allocate arrays:
    S₁    = zeros(Complex{FT},nquad_radius,n_mu)
    S₂    = zeros(Complex{FT},nquad_radius,n_mu)
    f₁₁   = zeros(FT, nquad_radius,n_mu)
    f₃₃   = zeros(FT, nquad_radius,n_mu)
    f₁₂   = zeros(FT, nquad_radius,n_mu)
    f₃₄   = zeros(FT, nquad_radius,n_mu)
    C_ext = zeros(FT, nquad_radius)
    C_sca = zeros(FT, nquad_radius)
    # Weights for the size distribution:
    wₓ = pdf.(aero.size_distribution,r)
    # pre multiply with wᵣ to get proper means eventually:
    wₓ .*= wᵣ
    # normalize (could apply a check whether cdf.(aero.size_distribution,r_max) is larger than 0.99:
    wₓ /= sum(wₓ)
    
    @showprogress 1 "Computing PhaseFunctions Siewert NAI-2 style ..." for i = 1:length(x_sizeParam)
        #println(i, " ", x_sizeParam[i])
        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = get_n_max(x_sizeParam[i])

        # In Domke methods, we want to pre-allocate these as 2D outside of this loop.
        an = (zeros(Complex{FT},n_max))
        bn = (zeros(Complex{FT},n_max))

        # Weighting for sums of 2n+1
        n_ = collect(FT,1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dn:
        y = x_sizeParam[i] * (aero.nᵣ-aero.nᵢ);
        nmx = round(Int, max(n_max, abs(y))+51 )
        Dn = zeros(Complex{FT},nmx)
        #@show typeof(Dn)
        #PhaseFunction.compute_mie_ab!(x_sizeParam[i],aero1.nᵣ-aero1.nᵢ*im,view(an,1:n_max),view(bn,1:n_max),Dn)
        #S1[i,:],S2[i,:] = PhaseFunction.compute_mie_S1S2(view(an,1:n_max), view(bn,1:n_max), leg_π, leg_τ)
        compute_mie_ab!(x_sizeParam[i],aero.nᵣ+aero.nᵢ*im,an,bn,Dn)
        compute_mie_S₁S₂!(an, bn, leg_π, leg_τ, view(S₁,i,:), view(S₂,i,:))
        
        # Compute Extinction and scattering cross sections: 
        C_sca[i] = 2pi/k^2 * (n_' * (abs2.(an) + abs2.(bn)))
        C_ext[i] = 2pi/k^2 * (n_' * real(an + bn))

        # Compute scattering matrix components per size parameter:
        f₁₁[i,:] = 0.5/x_sizeParam[i]^2  * real(abs2.(S₁[i,:]) + abs2.(S₂[i,:]));
        f₃₃[i,:] = 0.5/x_sizeParam[i]^2  * real(S₁[i,:] .* conj(S₂[i,:]) + S₂[i,:] .* conj(S₁[i,:]));
        f₁₂[i,:] = -0.5/x_sizeParam[i]^2  * real(abs2.(S₁[i,:]) - abs2.(S₂[i,:]));
        f₃₄[i,:] = -0.5/x_sizeParam[i]^2 * imag(S₁[i,:] .* conj(S₂[i,:]) - S₂[i,:] .* conj(S₁[i,:]));

    end
    bulk_C_sca =  sum(wₓ .* C_sca)
    bulk_C_ext =  sum(wₓ .* C_ext)
    bulk_f₁₁   =  sum(4π*r.^2 .*  wₓ .* f₁₁,dims=1)
    bulk_f₃₃   =  sum(4π*r.^2 .*  wₓ .* f₃₃,dims=1)
    bulk_f₁₂   =  sum(4π*r.^2 .*  wₓ .* f₁₂,dims=1)
    bulk_f₃₄   =  sum(4π*r.^2 .*  wₓ .* f₃₄,dims=1)

    # Normalize Phase function with bulk scattering cross section.
    bulk_f₁₁ /= bulk_C_sca 
    bulk_f₃₃ /= bulk_C_sca
    bulk_f₁₂ /= bulk_C_sca
    bulk_f₃₄ /= bulk_C_sca

    lMax = length(μ);
    P, P², R², T² = compute_legendre_poly(μ,lMax)

    # Compute Greek coefficients:
    α = zeros(FT,lMax)
    β = zeros(FT,lMax)
    δ = zeros(FT,lMax)
    γ = zeros(FT,lMax)
    ϵ = zeros(FT,lMax)
    ζ = zeros(FT,lMax)
    #@show size(avg_f11), size(P), size(w_μ), size(f11), size(wₓ)
    for l=0:length(β)-1
        # Factorial can be smarter through recursion I guess! 
        fac = (2l+1)/2 * sqrt(factorial(big(l-2))/factorial(big(l+2)))
        δ[l+1] = (2l+1)/2 * sum(w_μ' * (bulk_f₃₃[1,:] .* P[l+1,:]))
        β[l+1] = (2l+1)/2 * sum(w_μ' * (bulk_f₁₁[1,:] .* P[l+1,:]))
        γ[l+1] = fac * sum(w_μ' * (bulk_f₁₂[1,:] .* P²[l+1,:]))
        ϵ[l+1] = fac * sum(w_μ' * (bulk_f₃₄[1,:] .* P²[l+1,:]))
        ζ[l+1] = fac * sum(w_μ' * (bulk_f₃₃[1,:] .* R²[l+1,:] + bulk_f₁₁[1,:] .* T²[l+1,:]) )
        α[l+1] = fac * sum(w_μ' * (bulk_f₁₁[1,:] .* R²[l+1,:] + bulk_f₃₃[1,:] .* T²[l+1,:]) )
    end
    # return bulk_f₁₁, bulk_f₁₂, f₁₁, f₃₃, f₁₂,f₃₄,  C_ext, C_sca,bulk_C_sca, bulk_C_ext, α, β, γ, δ, ϵ, ζ
    return α, β, γ, δ, ϵ, ζ, bulk_C_sca, bulk_C_ext
end

function f_test(x)
    aero = UnivariateAerosol(LogNormal(log(x[1]), log(x[2])), 30.0, 2500, x[3], x[4])
    α, β, γ, δ, ϵ, ζ, bulk_C_sca, bulk_C_ext = calc_aer_opt_prop(NAI2(), aero, 0.55)
    return [α; β; γ; δ; ϵ; ζ; bulk_C_sca; bulk_C_ext]
          
end