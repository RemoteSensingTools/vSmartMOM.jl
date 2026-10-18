#!/usr/bin/env julia
# Retrieval-state plumbing of the vertical-grid convention (CPU, no RT solve).
# Run from test/: julia --project=.. \
#   ../sandbox/workflows/RRS_XCO2/inversion/test_vertical_grid_forward_mapping.jl
# Needs the spectroscopy/solar environment that VSmartMOMForward loads.

using LinearAlgebra
using Test

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "VSmartMOMForward.jl"))
using .VSmartMOMForward
using vSmartMOM

function _state(evaluator, psurf)
    x = zeros(30)
    x[1] = psurf
    x[VSmartMOMForward.CO2_RANGE] .= 400e-6
    x[last(VSmartMOMForward.CO2_RANGE)] = 360e-6
    x[VSmartMOMForward.LOG_AOD_RANGE] .= log.([0.2, 0.05, 0.01])
    x[VSmartMOMForward.LOG_HEIGHT_RANGE] .=
        [a.profile.μ for a in evaluator.base_parameters.scattering_params.rt_aerosols]
    x[VSmartMOMForward.SURFACE_RANGE] .= repeat([0.3, 0.0, 0.0], 3)
    return x
end

@testset "retrieval surface pressure follows the vertical-grid convention" begin
    default = OCOForwardEvaluator(; architecture=:CPU, float_type=Float32, nstreams=9)
    @test default.vertical_grid isa FixedPressureGrid
    altitude = OCOForwardEvaluator(; architecture=:CPU, float_type=Float32, nstreams=9,
                                   vertical_grid=FixedAltitudeGrid())
    reference = copy(default.base_parameters.p)
    @test reference[end] == 1000f0
    for (evaluator, grid) in ((default, FixedPressureGrid()), (altitude, FixedAltitudeGrid()))
        # At the reference pressure both conventions keep the template exactly.
        params = deepcopy(evaluator.base_parameters)
        apply_retrieval_state!(params, _state(evaluator, 1000.0), evaluator.tau_ref_scale;
                               vertical_grid=evaluator.vertical_grid)
        @test params.p == reference
        params = deepcopy(evaluator.base_parameters)
        apply_retrieval_state!(params, _state(evaluator, 985.0), evaluator.tau_ref_scale;
                               vertical_grid=evaluator.vertical_grid)
        @test params.p == surface_pressure_interfaces(grid, reference, 985f0)
        if grid isa FixedPressureGrid
            @test params.p[1:end-1] == reference[1:end-1]
        else
            @test params.p ./ reference ≈ fill(0.985f0, length(reference)) rtol=4eps(Float32)
        end
        # The template itself is never modified by a trial state.
        @test evaluator.base_parameters.p == reference
    end
    # XCO2 weights: bottom layer only (pressure grid) vs. none (altitude grid).
    x0, x1 = _state(default, 1000.0), _state(default, 985.0)
    @test column_averaged_co2_ppm(default, x0) == column_averaged_co2_ppm(altitude, x0)
    # Lower psurf thins the 360-ppm bottom layer, so pressure-grid XCO2 rises.
    @test column_averaged_co2_ppm(default, x1) > column_averaged_co2_ppm(default, x0)
    @test column_averaged_co2_ppm(altitude, x1) ≈ column_averaged_co2_ppm(altitude, x0) rtol=1e-6
end
