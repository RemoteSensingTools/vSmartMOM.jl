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
# Maximum expansion (see eq. A17 from de Rooij and Stap, 1984)
n_max = PhaseFunction.get_n_max(maxSizeParam)
x_sizeParam  = x*maxSizeParam/2 .+ maxSizeParam/2

an = (zeros(Complex{Float32},n_max,length(x)))
bn = (zeros(Complex{Float32},n_max,length(x)))

# Pre-allocate Dn:
y = maxSizeParam * (aero1.nᵣ-aero1.nᵢ);

nmx = round(Int, max(n_max, abs(y))+25 )
Dn = zeros(Complex{Float32},nmx)
# This is still not really geared towards GPU, will have to think about how we could do that
function run()
    kernel! = PhaseFunction.comp_ab!(CPU(),4)
    event = kernel!(x_sizeParam,an,bn,Dn,aero1.nᵣ-aero1.nᵢ*im,ndrange=length(x)); wait(event)
end

x_sizeParam_CU = CuArray(x_sizeParam)
an_CU = CuArray(an)
bn_CU = CuArray(bn)
Dn_CU = CuArray(Dn)

function runGPU()
    kernel! = PhaseFunction.comp_ab!(CUDADevice())
    event = kernel!(x_sizeParam_CU,an_CU,bn_CU,Dn_CU,aero1.nᵣ-aero1.nᵢ*im,ndrange=length(x)); wait(event)
end
@time run()
# Testing Siewert method now (for one size parameter first):
# 1) Gauss Legendre quadrature over μ
n_mu = 1000;
μ, w_μ = gausslegendre( n_mu )
leg_π = zeros(n_max,n_mu)
leg_τ = zeros(n_max,n_mu)
PhaseFunction.compute_mie_π_τ!(μ, n_max, leg_π, leg_τ)






l = @layout [a; b; c; d]
p1 = plot(x_sizeParam, real(an[1,:]), label="Real(a1)")
p2 = plot(x_sizeParam, imag(an[1,:]), label="Imag(a1)")
p3 = plot(x_sizeParam, real(bn[1,:]), label="Real(b1)")
p4 = plot(x_sizeParam, imag(bn[1,:]), label="Imag(b1)")
plot(p1, p2, p3,p4, layout = l)
xlabel!("Size Parameter")






