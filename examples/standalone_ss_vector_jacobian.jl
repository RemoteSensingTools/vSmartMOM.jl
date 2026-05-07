#!/usr/bin/env julia

using ForwardDiff
using Printf
using vSmartMOM
using vSmartMOM.StandaloneSS

FT = Float64

geometry = SSGeometry(
    μ₀ = FT(0.79),
    μv = FT[0.41, 0.73],
    Δϕ = FT[0.2, 1.1])

absorber = AbsorptionSSContributor(
    τ = FT[0.06 0.02;
           0.03 0.05])

wind0 = FT(4.0)
n_water = Complex{FT}(FT(1.34), FT(1e-8))

surface_from_wind(wind) = CoxMunkSSSurface(
    wind_speed = wind,
    n_water = n_water,
    include_whitecaps = false,
    shadowing = true)

config_from_wind(wind) = ExactSSConfig(
    geometry = geometry,
    surface = surface_from_wind(wind),
    contributors = (absorber,),
    I0 = FT[1.0, 0.8],
    polarization_type = vSmartMOM.Scattering.Stokes_IQ{FT}())

selector = SSMeasurementSelector(paths = :path2, stokes_indices = 1:2)

with_jacobians = run_exact_ss_with_jacobians(
    config_from_wind(wind0);
    paths = :path2,
    selector = selector)

dρ_dwind = surface_brdf_wind_jacobian(config_from_wind(wind0))
J_wind = chain_rule_combine_surface_brdf(
    with_jacobians.jacobians.path2.surface_brdf,
    dρ_dwind,
    selector)

fd_wind = ForwardDiff.jacobian(FT[wind0]) do wind
    result = run_exact_ss(config_from_wind(wind[1]); paths = :path2)
    selected_measurements(result, selector)
end

@assert J_wind ≈ fd_wind rtol=2e-10 atol=1e-12

println("StandaloneSS vector Cox-Munk wind Jacobian")
println("measurements: ", size(with_jacobians.measurements))
println("J_wind:      ", size(J_wind))
println("max |analytic - FD| = ", maximum(abs.(J_wind .- fd_wind)))
println()
println("row  radiance       d radiance / d wind")
for i in eachindex(with_jacobians.measurements)
    @printf("%3d  %12.6e  %12.6e\n",
            i, with_jacobians.measurements[i], J_wind[i, 1])
end
