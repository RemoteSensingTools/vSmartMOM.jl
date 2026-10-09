#!/usr/bin/env julia

using Random
using Test

include(joinpath(@__DIR__, "OptimalEstimation.jl"))
include(joinpath(@__DIR__, "Round6FixedSIF.jl"))
using .OptimalEstimation: ForwardEvaluation
using .Round6FixedSIF

@testset "round-6 fixed-SIF state map" begin
    off = Round6SIFMap(false)
    state = collect(1.0:28.0)
    core_off = expand_round6_state(off, state)
    @test active_state_count(off) == 28
    @test active_core_indices(off) == collect(1:28)
    @test core_off[1:28] == state
    @test core_off[29:30] == [0.0, 0.0]

    anchor = 0.004818031987713776
    slope = 1.2291230681458325e-5
    on = Round6SIFMap(true; Lnu759=anchor, mSIF=slope)
    core_on = expand_round6_state(on, state)
    @test core_on[29] ==
        anchor + slope * DELTA_NU_760_MINUS_759_CM1
    @test core_on[30] == slope
    @test core_on[29] + core_on[30] *
        (NU_759_CM1 - NU_760_CM1) ≈ anchor atol=2e-16 rtol=0

    @test_throws ArgumentError Round6SIFMap(false; Lnu759=1e-3)
    @test_throws ArgumentError Round6SIFMap(false; mSIF=1e-5)
    @test_throws ArgumentError Round6SIFMap(true; Lnu759=-1e-3)
    @test_throws DimensionMismatch expand_round6_state(on, state[1:27])
end

@testset "round-6 fixed-column Jacobian" begin
    Random.seed!(0x500759)
    nmeasurement = 17
    operator = randn(nmeasurement, CORE_STATE_COUNT)
    offset = randn(nmeasurement)
    ranges = [1:5, 6:11, 12:17]
    core_evaluator(x) = ForwardEvaluation(
        offset + operator * x, operator, ranges;
        timing=(forward_seconds=1.0, jacobian_seconds=0.5,
                instrument_seconds=0.25))

    for map in (Round6SIFMap(false),
                Round6SIFMap(true;
                    Lnu759=0.004818031987713776,
                    mSIF=1.2291230681458325e-5))
        state = randn(ROUND6_STATE_COUNT)
        wrapped = Round6ForwardEvaluator(core_evaluator, map)
        evaluation = wrapped(state)
        @test evaluation.measurement ==
            offset + operator * expand_round6_state(map, state)
        @test evaluation.jacobian == operator[:, 1:ROUND6_STATE_COUNT]
        @test size(evaluation.jacobian) == (nmeasurement, 28)

        step = 1e-6
        for index in eachindex(state)
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
        Round6SIFMap(false), zeros(4, 29))
end
