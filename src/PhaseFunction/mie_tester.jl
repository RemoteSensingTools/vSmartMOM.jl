# Quick testing function:
using RadiativeTransfer.PhaseFunction
using Distributions
using FastGaussQuadrature
using Plots

# Define median radius, sigma and wavelength (all in micron!)
μ  = 0.3
σ  = 6.82
wl = 0.55


# Generate aerosol:
aero1 = PhaseFunction.UnivariateAerosol(LogNormal(log(μ), log(σ)), 30.0, 10000,1.3,0.01)
# Obtain Gauss Legendre Quadrature Points:
x,w = gausslegendre( aero1.nquad_radius )
maxSizeParam = 2π * aero1.r_max
x_sizeParam  = x*maxSizeParam/2 .+ maxSizeParam/2

an = (zeros(Complex{Float32},350,length(x)))
bn = (zeros(Complex{Float32},350,length(x)))

# This is still not really geared towards GPU, will have to think about how we could do that
function run()
    kernel! = PhaseFunction.comp_ab!(CPU(),4)
    event = kernel!(x_sizeParam,an,bn,aero1.nᵣ-aero1.nᵢ*im,ndrange=length(x)); wait(event)
end

@time run()

l = @layout [a; b; c; d]
p1 = plot(x_sizeParam, real(an[1,:]), label="Real(a1)")
p2 = plot(x_sizeParam, imag(an[1,:]), label="Imag(a1)")
p3 = plot(x_sizeParam, real(bn[1,:]), label="Real(b1)")
p4 = plot(x_sizeParam, imag(bn[1,:]), label="Imag(b1)")
plot(p1, p2, p3,p4, layout = l)
xlabel!("Size Parameter")






