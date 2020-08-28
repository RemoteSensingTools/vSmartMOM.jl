# Quick testing function:
using Revise
using RadiativeTransfer.PhaseFunction
using Distributions
using FastGaussQuadrature
using Plots
using KernelAbstractions
using CUDA

# Define median radius, sigma and wavelength (all in micron!)
μ  = 0.3
σ  = 6.82
wl = 0.55


# Generate aerosol:
aero1 = PhaseFunction.UnivariateAerosol(LogNormal(log(μ), log(σ)), 30.0, 10000,1.3,0.01)

# Obtain Gauss Legendre Quadrature Points:
x,w = gausslegendre( aero1.nquad_radius )
maxSizeParam = 2π * aero1.r_max/wl
n_max = PhaseFunction.get_n_max(maxSizeParam)
x_sizeParam  = x*maxSizeParam/2 .+ maxSizeParam/2

# 1) Gauss Legendre quadrature over μ
n_mu = 1000;
μ, w_μ = gausslegendre( n_mu )

leg_π = zeros(n_max,n_mu)
leg_τ = zeros(n_max,n_mu)
PhaseFunction.compute_mie_π_τ!(μ, n_max, leg_π, leg_τ)
# Choose just one for testing:


function testSiewert(x_sizeParam,leg_π,leg_τ, aero1)
    size1 = size(leg_π)[2]
    S1 = zeros(Complex{Float32},length(x_sizeParam),size1)
    S2 = zeros(Complex{Float32},length(x_sizeParam),size1)
    @showprogress 1 "Computing PhaseFunctions ..." for i = 1:length(x_sizeParam)
        #println(i, " ", x_sizeParam[i])
        # Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
        n_max = PhaseFunction.get_n_max(x_sizeParam[i])
        an = (zeros(Complex{Float32},n_max))
        bn = (zeros(Complex{Float32},n_max))

        # Pre-allocate Dn:
        y = x_sizeParam[i] * (aero1.nᵣ-aero1.nᵢ);
        nmx = round(Int, max(n_max, abs(y))+25 )
        Dn = zeros(Complex{Float32},nmx)

        #PhaseFunction.compute_mie_ab!(x_sizeParam[i],aero1.nᵣ-aero1.nᵢ*im,view(an,1:n_max),view(bn,1:n_max),Dn)
        #S1[i,:],S2[i,:] = PhaseFunction.compute_mie_S1S2(view(an,1:n_max), view(bn,1:n_max), leg_π, leg_τ)
        PhaseFunction.compute_mie_ab!(x_sizeParam[i],aero1.nᵣ-aero1.nᵢ*im,an,bn,Dn)
        S1[i,:],S2[i,:] = PhaseFunction.compute_mie_S1S2(an, bn, leg_π, leg_τ)
    end
    return S1,S2
end

@time testSiewert(x_sizeParam, leg_π,leg_τ, aero1)
    
# Choose index and plot phase functions:
i = 500
plot(μ, log.(real(S1[i,:] .* conj(S1[i,:]) + S2[i,:] .* conj(S2[i,:]))))
i = 2000
plot!(μ, log.(real(S1[i,:] .* conj(S1[i,:]) + S2[i,:] .* conj(S2[i,:]))))
i = 10000
plot!(μ, log.(real(S1[i,:] .* conj(S1[i,:]) + S2[i,:] .* conj(S2[i,:]))))






