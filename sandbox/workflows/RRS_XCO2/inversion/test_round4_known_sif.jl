#!/usr/bin/env julia

using Random
using Test

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "Round4KnownSIF.jl"))
using .OptimalEstimation: ForwardEvaluation
using .Round4KnownSIF

@testset "round-4 known-SIF state map" begin
    @test NU_759_CM1 > NU_760_CM1
    @test DELTA_NU_760_MINUS_759_CM1 < 0

    off = Round4SIFMap(false)
    xoff = collect(1.0:28.0)
    core_off = expand_round4_state(off, xoff)
    @test active_state_count(off) == 28
    @test active_core_indices(off) == collect(1:28)
    @test core_off[1:28] == xoff
    @test core_off[29:30] == [0.0, 0.0]

    known = 4.7e-3
    slope = 1.2e-5
    on = Round4SIFMap(true; Lν759=known)
    xon = vcat(xoff, slope)
    core_on = expand_round4_state(on, xon)
    @test active_state_count(on) == 29
    @test active_core_indices(on) == vcat(collect(1:28), 30)
    @test core_on[30] == slope
    @test core_on[29] ==
        known + slope * DELTA_NU_760_MINUS_759_CM1
    @test core_on[29] + core_on[30] * (NU_759_CM1 - NU_760_CM1) ≈
        known atol=2e-16 rtol=0

    @test_throws ArgumentError Round4SIFMap(false; Lν759=1e-3)
    @test_throws ArgumentError Round4SIFMap(true; Lν759=-1e-3)
    @test_throws DimensionMismatch expand_round4_state(on, xoff)
end

@testset "round-4 Jacobian chain rule and forward wrapper" begin
    Random.seed!(0x759760)
    nmeasurement = 17
    operator = randn(nmeasurement, CORE_STATE_COUNT)
    offset = randn(nmeasurement)
    ranges = [1:5, 6:11, 12:17]
    core_evaluator(x) = ForwardEvaluation(
        offset + operator * x, operator, ranges;
        timing=(forward_seconds=1.0, jacobian_seconds=0.5,
                instrument_seconds=0.25))

    for map in (Round4SIFMap(false),
                Round4SIFMap(true; Lν759=4.7e-3))
        nstate = active_state_count(map)
        state = randn(nstate)
        wrapped = Round4ForwardEvaluator(core_evaluator, map)
        evaluation = wrapped(state)
        expected_core = expand_round4_state(map, state)

        @test evaluation.measurement == offset + operator * expected_core
        @test evaluation.jacobian == reduce_core_jacobian(map, operator)
        @test evaluation.band_ranges == ranges
        @test evaluation.timing.forward_seconds == 1.0

        # Check every reduced column independently, including the composite
        # mSIF column, against a central perturbation through the full map.
        step = 1e-6
        for index in 1:nstate
            plus = copy(state); plus[index] += step
            minus = copy(state); minus[index] -= step
            finite_difference =
                (wrapped(plus).measurement - wrapped(minus).measurement) /
                (2step)
            @test isapprox(
                finite_difference, evaluation.jacobian[:, index];
                atol=2e-8, rtol=2e-8)
        end
    end

    @test_throws DimensionMismatch reduce_core_jacobian(
        Round4SIFMap(false), zeros(4, 29))
end
