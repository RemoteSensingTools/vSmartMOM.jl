# Quick testing function:
using Revise
using RadiativeTransfer.Scattering
using Distributions
using FastGaussQuadrature
using Plots
using KernelAbstractions
using CUDA
using ProgressMeter

# Define median radius, sigma and wavelength (all in micron!)
μ  = 0.3
σ  = 6.82
wl = 0.55

FT = Float64

# Generate aerosol:
aero1 = Scattering.UnivariateAerosol(LogNormal(log(μ), log(σ)), 30.0, 10000,1.3,-1e-6)

# Obtain Gauss Legendre Quadrature Points:
x,w = gausslegendre( aero1.nquad_radius )
maxSizeParam = 2π * aero1.r_max/wl
n_max = Scattering.get_n_max(maxSizeParam)
x_sizeParam  = x*maxSizeParam/2 .+ maxSizeParam/2

#x_sizeParam=
# 1) Gauss Legendre quadrature over μ (requires 2Nmax-1 to be accurate!)
n_mu = 2n_max-1;

μ, w_μ = gausslegendre( n_mu )

leg_π = zeros(n_max,n_mu)
leg_τ = zeros(n_max,n_mu)

# Pre-compute π and τ:
Scattering.compute_mie_π_τ!(μ, n_max, leg_π, leg_τ)

function testSiewert(x_sizeParam,leg_π,leg_τ, aero1, λ,μ,w_μ )
    size1 = size(leg_π)[2]
    # Wavenumber:
    k = 2π/λ
    # Radius
    r = x_sizeParam/2π*wl

    S₁  = zeros(Complex{FT},length(x_sizeParam),size1)
    S₂  = zeros(Complex{FT},length(x_sizeParam),size1)
    f₁₁ = zeros(length(x_sizeParam),size1)
    f₃₃ = zeros(length(x_sizeParam),size1)
    f₁₂ = zeros(length(x_sizeParam),size1)
    f₃₄ = zeros(length(x_sizeParam),size1)
    C_ext = zeros(length(x_sizeParam))
    C_sca = zeros(length(x_sizeParam))
    # Weights for the size distribution:
    wₓ = pdf.(aero1.size_distribution,r)
    # Try Gauss-Legendre later (add to wx)

    # normalize:
    wₓ /= sum(wₓ)
    #wₓ /= maximum(x_sizeParam/2π*wl)

    @showprogress 1 "Computing PhaseFunctions Siewert NAI-2 style ..." for i = 1:length(x_sizeParam)
        #println(i, " ", x_sizeParam[i])
        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = Scattering.get_n_max(x_sizeParam[i])
        an = (zeros(Complex{FT},n_max))
        bn = (zeros(Complex{FT},n_max))
        # Weighting for sums of 2n+1
        n_ = collect(1:n_max);
        n_ = 2n_ .+ 1

        # Pre-allocate Dn:
        y = x_sizeParam[i] * (aero1.nᵣ-aero1.nᵢ);
        nmx = round(Int, max(n_max, abs(y))+25 )
        Dn = zeros(Complex{FT},nmx)

        #Scattering.compute_mie_ab!(x_sizeParam[i],aero1.nᵣ-aero1.nᵢ*im,view(an,1:n_max),view(bn,1:n_max),Dn)
        #S1[i,:],S2[i,:] = Scattering.compute_mie_S1S2(view(an,1:n_max), view(bn,1:n_max), leg_π, leg_τ)
        Scattering.compute_mie_ab!(x_sizeParam[i],aero1.nᵣ-aero1.nᵢ*im,an,bn,Dn)
        S₁[i,:],S₂[i,:] = Scattering.compute_mie_S1S2(an, bn, leg_π, leg_τ)
        
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
    P, P², R², T² = Scattering.compute_legendre_poly(μ,lMax)
    # Compute Greek coefficients:
    α = zeros(lMax)
    β = zeros(lMax)
    δ = zeros(lMax)
    γ = zeros(lMax)
    ϵ = zeros(lMax)
    ζ = zeros(lMax)
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
    return bulk_f₁₁, bulk_f₁₂, f₁₁, f₃₃, f₁₂,f₃₄,  C_ext, C_sca,bulk_C_sca, bulk_C_ext, α, β, γ, δ, ϵ, ζ
end


@time bulk_f₁₁, bulk_f₁₂, f₁₁, f₃₃, f₁₂,f₃₄,  C_ext, C_sca,avg_C_sca, avg_C_ext, α, β, γ, δ, ϵ, ζ= testSiewert(x_sizeParam, leg_π,leg_τ, aero1,wl,μ,w_μ)
    
l = @layout [a ; b; c; d; e; f]
p1 = plot(α,label="α")
p2 = plot(β,label="β")
p3 = plot(γ,label="γ")
p4 = plot(δ,label="δ")
p5 = plot(ϵ,label="ϵ")
p6 = plot(ζ,label="ζ")

plot(p1, p2, p3,p4,p5,p6, layout = l)
plot!(size=(400,600))

# Test single particle:
function testSingle(index)
    lMax = 750;
    P⁰, P², R², T² = Scattering.compute_legendre_poly(μ,lMax)
    ϵ_s = zeros(lMax);
    #index = 100
    for l=0:lMax-1
        fac = (2l+1)/2 * sqrt(factorial(big(l-2))/factorial(big(l+2)))
        ϵ_s[l+1] = fac * sum(w_μ' * (f₃₄[index,:] .* P²[l+1,:]))
    end
    return ϵ_s
end

plot(testSingle(1000))


