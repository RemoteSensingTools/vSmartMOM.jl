#!/usr/bin/env julia
# Compare the frozen Round7 correction against the unchanged Round5 instrument.
using NCDatasets
using Test
const BASE = joinpath(dirname(Base.active_project()), "sandbox", "workflows", "RRS_XCO2", "inversion")
include(joinpath(BASE, "instrument", "SyntheticOCO2.jl"))
using .SyntheticOCO2

function validate_round7_instrument(directory)
    files = sort(filter(name -> occursin(r"^OCO2round7_\d{3}\.nc$", name), readdir(directory)))
    @test length(files) == 80
    largest = 0.0
    @testset "Round7 instrument parity and exact correction identities" begin
        for file in files
            NCDataset(joinpath(directory, file), "r") do ds
                nu = Float64.(ds["correction_wavenumber"][:])
                stokes = permutedims(Float64.(ds["correction_stokes"][:, :]))
                analyzer = Float64.(ds["o2a_analyzer_coefficients"][:])
                result = process_stokes_spectrum(1e7 ./ nu, stokes, analyzer, band_spec(:o2a))
                correction = Float64.(ds["correction_observation"][:])
                difference = maximum(abs.(result .- correction[1:934]))
                largest = max(largest, difference)
                @test isapprox(result, correction[1:934]; rtol=2e-13, atol=2e-13)
                @test all(iszero, correction[935:end])
                unc = Float64.(ds["measurement_uncorrected"][:])
                imperfect = Float64.(ds["measurement_imperfectly_corrected"][:])
                @test imperfect == unc - correction
                noise = Float64.(ds["injected_measurement_noise"][:, :])
                @test Float64.(ds["imperfectly_corrected_perturbed"][:, :]) == imperfect .+ noise
                @test Float64.(ds["uncorrected_perturbed"][:, :]) == unc .+ noise
                @test all(iszero, noise[:, 11])
                @test all(iszero, Float64.(ds["normalized_noise_draw"][:, 11]))
            end
        end
    end
    println("maximum Python/Julia instrument absolute difference = ", largest)
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && validate_round7_instrument(only(ARGS))
