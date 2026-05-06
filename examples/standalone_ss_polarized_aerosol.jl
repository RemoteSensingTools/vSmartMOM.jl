#!/usr/bin/env julia

using Printf
using vSmartMOM
using vSmartMOM.Scattering
using vSmartMOM.StandaloneSS

FT = Float64

geometry = SSGeometry(
    μ₀ = FT(0.82),
    μv = FT[0.35, 0.55, 0.75],
    Δϕ = deg2rad.(FT[0, 60, 120]))

phase = SyntheticPolarizedHenyeyGreensteinPhaseFunction(
    g = FT(0.45),
    polarization_fraction = FT(0.6))
greek = greek_coefficients(phase; l_max = 32, nquad = 96)

aerosol = GreekCoefsSSContributor(
    greek_coefs = greek,
    ϖ = FT(0.92),
    τ = reshape(FT[0.03, 0.04, 0.02], 3, 1))
absorption = AbsorptionSSContributor(
    τ = reshape(FT[0.01, 0.01, 0.01], 3, 1))

surface = CoxMunkSSSurface(
    wind_speed = FT(5.0),
    n_water = Complex{FT}(FT(1.34), FT(1e-8)),
    include_whitecaps = false,
    shadowing = true)

config = ExactSSConfig(
    geometry = geometry,
    surface = surface,
    contributors = (aerosol, absorption),
    I0 = FT[1],
    polarization_type = Stokes_IQU{FT}())

result = run_exact_ss(config; paths = :paths_1_2)

@assert any(abs.(result.path1[:, 2:3, :]) .> 0)
@assert any(abs.(result.path2[:, 2:3, :]) .> 0)

println("StandaloneSS polarizing aerosol example")
println("view  mu_v   dphi_deg    path1_I      path1_Q      path1_U      path2_I      path2_Q      path2_U")
for iv in eachindex(geometry.μv)
    @printf("%4d  %.2f  %8.1f  %11.4e %11.4e %11.4e %11.4e %11.4e %11.4e\n",
            iv,
            geometry.μv[iv],
            rad2deg(geometry.Δϕ[iv]),
            result.path1[iv, 1, 1],
            result.path1[iv, 2, 1],
            result.path1[iv, 3, 1],
            result.path2[iv, 1, 1],
            result.path2[iv, 2, 1],
            result.path2[iv, 3, 1])
end
