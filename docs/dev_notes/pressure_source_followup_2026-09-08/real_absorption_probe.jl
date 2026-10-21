# Run from test/ with --project=. and an optional native ABSCO CO2 file argument.
# Synthetic regression fixtures are public; spectroscopy is read-only and is
# never copied into this repository. Includes the complete pressure regression
# first so the same construction/coordinate helpers are exercised.
using Test, LinearAlgebra, vSmartMOM
import AtmosphericAbsorption
include(normpath(joinpath(@__DIR__, "../../../test/test_absorption_pressure.jl")))

@testset "Real spectroscopy: fixed-grid surface pressure" begin
    for kind in (:line, :absco), species in (:fixed, :variable)
        kind === :absco && isempty(ARGS) && continue
        absorber = kind === :line ? _pressure_test_line_model(Float64) :
            AtmosphericAbsorption.read_absco(only(ARGS); FT=Float64,
                architecture=AtmosphericAbsorption.CPU(),
                wavenumber_range=(6199.7, 6200.3), broadener_vmr=:all)
        params = _pressure_test_parameters(Float64, absorber; species)
        if species === :variable
            params.absorption_params.variable_molecules[1] .= "CO2"
        end
        # Empty LUT storage selects the genuine configured HITRAN LBL branch,
        # rather than injecting a preconstructed line model through the LUT API.
        kind === :line && empty!(params.absorption_params.luts)
        model, lin = model_from_parameters(LinMode(), params)
        @test maximum(model.τ_abs[1]) > 0
        h = 0.02
        plus, minus = deepcopy(params), deepcopy(params)
        plus.p[end] += h
        minus.p[end] -= h
        mplus, mminus = model_from_parameters(plus), model_from_parameters(minus)
        fdτ = (mplus.τ_abs[1] - mminus.τ_abs[1]) / (2h)
        optical_error = norm(lin.τ̇_abs_psurf[1] - fdτ) / norm(fdτ)
        @test optical_error < 1e-5
        analytic = rt_run(model, lin, 0, size(lin.τ̇_abs[1], 1), 1)
        fdR = (rt_run(mplus).toa - rt_run(mminus).toa) / (2h)
        dR = analytic.toa_jacobian[:, :, :, 1]
        radiance_error = norm(dR - fdR) / norm(fdR)
        @test radiance_error < 1e-4
        weights = [1.0 2 1 0 0; 0 0 1 2 1] ./ 4
        convolved_error = norm(weights * vec(dR - fdR)) / norm(weights * vec(fdR))
        @test convolved_error < 1e-4
        println((; kind, species, optical_error, radiance_error, convolved_error))
    end
end
