#!/usr/bin/env julia

using Test

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
include(joinpath(@__DIR__, "OCO2Noise.jl"))
using .SyntheticOCO2
using .OCO2Noise

@testset "OCO-2 L1B ATBD noise model" begin
    @test noise_band_spec(:o2a).max_ms == 7.00e20
    @test noise_band_spec(:o2a).min_ms == 7.50e16
    @test noise_band_spec(:weak_co2).max_ms == 2.45e20
    @test noise_band_spec(:strong_co2).max_ms == 1.25e20
    @test_throws ArgumentError noise_band_spec(:unknown)

    wavelength = [765.0, 1606.5, 2062.0]
    radiance = [1.0, 2.0, 3.0]
    photon = energy_to_photon_radiance(wavelength, radiance)
    @test photon_to_energy_radiance(wavelength, photon) ≈ radiance rtol=2e-15

    max_ms = 7.00e20
    c_photon = [0.01, 0.02]
    c_background = [0.005, 0.006]
    signal = [0.0, max_ms]
    expected = max_ms / 100 .* sqrt.([
        c_background[1]^2,
        100 * c_photon[2]^2 + c_background[2]^2,
    ])
    @test noise_equivalent_radiance_photon(
        signal, max_ms, c_photon, c_background) ≈ expected rtol=2e-15
    @test noise_equivalent_radiance_photon(
        -signal, max_ms, c_photon, c_background) ≈ expected rtol=2e-15

    spec = noise_band_spec(:o2a)
    test_photon = [spec.min_ms / 2, spec.max_ms * 2]
    test_wavelength = [760.0, 770.0]
    test_energy = photon_to_energy_radiance(test_wavelength, test_photon)
    statistics = noise_statistics(
        test_wavelength, test_energy, fill(0.01, 2), fill(0.005, 2), spec)
    @test statistics.below_min_ms == [true, false]
    @test statistics.above_max_ms == [false, true]
    @test statistics.variance_energy == statistics.nen_energy .^ 2
    @test all(statistics.snr .> 0)

    coefficient_path = joinpath(@__DIR__, "representative_snr_coefficients.nc")
    if isfile(coefficient_path)
        coefficients = read_representative_snr_coefficients(coefficient_path)
        @test Set(keys(coefficients)) == Set(getfield.(BAND_SPECS, :name))
        for spec in BAND_SPECS
            coefficient = coefficients[spec.name]
            @test coefficient.wavelength == synthetic_grid(spec)
            @test all(coefficient.c_photon .> 0)
            @test all(coefficient.c_background .> 0)
            @test all(coefficient.extrapolated_source_count .>= 0)
        end
        @test all(coefficients[:o2a].extrapolated_source_count .== 0)
        @test all(coefficients[:weak_co2].extrapolated_source_count .== 0)
        @test count(>(0), coefficients[:strong_co2].extrapolated_source_count) == 29
    end
end
