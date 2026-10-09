#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "RetrievalCases.jl"))
using .RetrievalCases

const CAMPAIGN_ROOT = normpath(joinpath(
    @__DIR__, "..", "bottom_layer_XCO2_retrievals"))
const TRUTH_TABLE = joinpath(CAMPAIGN_ROOT, "truth", "true_states.dat")
const CURRY_LAUNCHER = joinpath(
    @__DIR__, "run_bottom_layer_retrieval_partition.sh")
const WURST_LAUNCHER = joinpath(
    @__DIR__, "run_bottom_layer_aerosol_retrievals.sh")

function parse_state_spec(value)
    states = Int[]
    for token in split(value, ',')
        bounds = parse.(Int, split(token, '-'))
        length(bounds) == 1 ? push!(states, only(bounds)) :
            append!(states, bounds[1]:bounds[2])
    end
    return states
end

function global_state_spec(source, variable)
    matched = match(Regex(variable * "=\"([^\"]+)\""), source)
    isnothing(matched) && error("missing $variable")
    return parse_state_spec(only(matched.captures))
end

function partition_blocks(source, variable)
    matches = collect(eachmatch(
        Regex("(?m)^\\s*" * variable * "=\\(([^)]*)\\)"), source))
    length(matches) == 2 || error(
        "expected two $variable assignments, found $(length(matches))")
    return map(matches) do matched
        states = Int[]
        for block in eachmatch(r"\"([0-9]+):([0-9]+)\"",
                               only(matched.captures))
            append!(states, parse(Int, block.captures[1]):
                            parse(Int, block.captures[2]))
        end
        states
    end
end

function assert_ascending_co2_blocks(states, truth_by_index)
    @test length(states) % 5 == 0
    for offset in 1:5:length(states)
        block = states[offset:offset + 4]
        @test [truth_by_index[state].bottom_co2_ppm for state in block] ==
            [360.0, 380.0, 400.0, 420.0, 440.0]
    end
end

@testset "bottom-layer cross-host phase schedule" begin
    curry_source = read(CURRY_LAUNCHER, String)
    wurst_source = read(WURST_LAUNCHER, String)
    truth_by_index = Dict(
        truth.state_index => truth for truth in read_truth_cases(TRUTH_TABLE))

    curry_nosif = partition_blocks(curry_source, "NOSIF_BLOCKS")
    curry_sif = partition_blocks(curry_source, "SIF_BLOCKS")
    wurst_nosif = partition_blocks(wurst_source, "NOSIF_BLOCKS")
    wurst_sif = partition_blocks(wurst_source, "SIF_BLOCKS")

    all_nosif = global_state_spec(curry_source, "ALL_NOSIF_SPEC")
    all_sif = global_state_spec(curry_source, "ALL_SIF_SPEC")
    @test all_nosif == global_state_spec(wurst_source, "ALL_NOSIF_SPEC")
    @test all_sif == global_state_spec(wurst_source, "ALL_SIF_SPEC")
    @test length(all_nosif) == length(unique(all_nosif)) == 40
    @test length(all_sif) == length(unique(all_sif)) == 40
    @test isempty(intersect(all_nosif, all_sif))
    @test sort(vcat(all_nosif, all_sif)) == collect(1:80)

    curry_nosif_all = vcat(curry_nosif...)
    curry_sif_all = vcat(curry_sif...)
    wurst_nosif_all = vcat(wurst_nosif...)
    wurst_sif_all = vcat(wurst_sif...)
    @test isempty(intersect(curry_nosif_all, wurst_nosif_all))
    @test isempty(intersect(curry_sif_all, wurst_sif_all))
    @test sort(vcat(curry_nosif_all, wurst_nosif_all)) == sort(all_nosif)
    @test sort(vcat(curry_sif_all, wurst_sif_all)) == sort(all_sif)

    @test all(state -> truth_by_index[state].aerosol_case == :none &&
                       truth_by_index[state].sif_case == :off,
              curry_nosif_all)
    @test all(state -> truth_by_index[state].aerosol_case == :none &&
                       truth_by_index[state].sif_case != :off,
              curry_sif_all)
    @test all(state -> truth_by_index[state].aerosol_case != :none &&
                       truth_by_index[state].sif_case == :off,
              wurst_nosif_all)
    @test all(state -> truth_by_index[state].aerosol_case != :none &&
                       truth_by_index[state].sif_case != :off,
              wurst_sif_all)

    # Curry owns complete five-CO2-value clear blocks. Wurst's retained
    # assignments deliberately split forest no-SIF states 051 and 052:055
    # across its two devices.
    for partition in vcat(curry_nosif, curry_sif)
        assert_ascending_co2_blocks(partition, truth_by_index)
    end

    expected_noise_order = repeat(vcat(11, collect(1:10)), inner=2)
    expected_class_order = repeat([:corrected, :uncorrected], 11)
    for state in 1:80
        experiments = build_experiments(
            [truth_by_index[state]];
            measurement_directory="/unused",
            noise_directory="/unused",
            validate_inputs=false)
        @test getfield.(experiments, :noise_index) == expected_noise_order
        @test getfield.(experiments, :measurement_class) == expected_class_order
    end

    for source in (curry_source, wurst_source)
        no_sif_work = findfirst("run_owned_blocks nosif off", source)
        no_sif_barrier = findfirst(
            "wait_for_global_phase nosif \"\${ALL_NOSIF_SPEC}\"", source)
        sif_work = findfirst("run_owned_blocks sif on", source)
        sif_barrier = findfirst(
            "wait_for_global_phase sif \"\${ALL_SIF_SPEC}\"", source)
        @test !isnothing(no_sif_work)
        @test !isnothing(no_sif_barrier)
        @test !isnothing(sif_work)
        @test !isnothing(sif_barrier)
        @test first(no_sif_work) < first(no_sif_barrier) <
              first(sif_work) < first(sif_barrier)
    end
end
