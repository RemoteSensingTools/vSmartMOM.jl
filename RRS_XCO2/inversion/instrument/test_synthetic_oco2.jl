#!/usr/bin/env julia

using NCDatasets
using Test

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
using .SyntheticOCO2

@testset "synthetic OCO-2 grids" begin
    @test length(synthetic_grid(band_spec(:o2a))) == 934
    @test length(synthetic_grid(band_spec(:weak_co2))) == 807
    @test length(synthetic_grid(band_spec(:strong_co2))) == 1001
    @test first(synthetic_grid(band_spec(:o2a))) == 758.0
    @test last(synthetic_grid(band_spec(:o2a))) == 771.995
    @test last(synthetic_grid(band_spec(:strong_co2))) == 2082.0
end

@testset "OCO analyzer projection" begin
    coefficients = [0.5, -0.3, -0.4, 0.0]
    stokes = [2.0 4.0; 0.1 -0.2; -0.3 0.4]
    expected_raw = 0.5 .* stokes[1, :] .-
                   (-0.3) .* stokes[2, :] .+
                   (-0.4) .* stokes[3, :]
    @test project_oco_analyzer(stokes, coefficients) == expected_raw
    @test project_oco_analyzer(stokes, coefficients; normalize=:m11) ==
          expected_raw ./ 0.5
    @test_throws ArgumentError project_oco_analyzer(
        stokes, coefficients; normalize=:determinant)
end

@testset "spectral-density Jacobian" begin
    wavelength = [500.0, 1000.0]
    per_cm = [1.0, 2.0]
    @test per_wavenumber_to_per_wavelength(wavelength, per_cm) ==
          per_cm .* 1e7 ./ wavelength .^ 2
end

@testset "Gaussian convolution and resampling" begin
    # Descending, mildly nonuniform source wavelengths exercise the truth-map
    # uniform-wavenumber ordering and the trapezoidal quadrature weights.
    source_wavelength = reverse([757.0 + 0.0048i + 2e-7i^2 for i in 0:3400])
    source_constant = fill(7.25, length(source_wavelength))
    target = synthetic_grid(band_spec(:o2a))
    convolved = gaussian_convolve_resample(
        source_wavelength, source_constant, target, 0.04)
    @test convolved ≈ source_constant[1] .* ones(length(target)) atol=2e-12

    # A normalized Gaussian leaves an affine spectrum unchanged at its center
    # on a well-resolved symmetric grid, apart from tiny quadrature error.
    affine = 2 .+ 0.3 .* source_wavelength
    convolved_affine = gaussian_convolve_resample(
        source_wavelength, affine, target, 0.04)
    @test convolved_affine ≈ 2 .+ 0.3 .* target atol=2e-8

    # The fixed instrument operator must act on every linearized parameter
    # column exactly as it acts on a forward perturbation.
    coefficients = [0.5, -0.47, -0.17, 0.0]
    nsource = length(source_wavelength)
    base = vcat(reshape(2 .+ 0.01 .* source_wavelength, 1, nsource),
                reshape(0.1 .* sin.(source_wavelength), 1, nsource),
                reshape(0.05 .* cos.(source_wavelength), 1, nsource))
    perturbations = Array{Float64}(undef, 3, nsource, 2)
    perturbations[:, :, 1] = 0.02 .* base
    perturbations[:, :, 2] = reverse(base; dims=2) .* 0.01
    processed_jacobian = process_stokes_jacobian(
        source_wavelength, perturbations, coefficients, band_spec(:o2a))
    epsilon = 1e-5
    processed_base = process_stokes_spectrum(
        source_wavelength, base, coefficients, band_spec(:o2a))
    for parameter_index in 1:2
        processed_perturbed = process_stokes_spectrum(
            source_wavelength,
            base .+ epsilon .* perturbations[:, :, parameter_index],
            coefficients,
            band_spec(:o2a))
        finite_difference = (processed_perturbed .- processed_base) ./ epsilon
        @test finite_difference ≈ processed_jacobian[:, parameter_index] rtol=2e-9 atol=2e-9
    end
end

@testset "strong-band source-support guard" begin
    # Retrievals recreate the truth product's 987-point base plus its
    # eight-point short-wavelength shoulder as one regular 0.1 cm^-1 grid.
    retrieval_wavenumber = collect((1e7 / 2084):0.1:(1e7 / 2041.7))
    @test length(retrieval_wavenumber) == 995
    retrieval_extension = required_wavenumber_extensions(
        retrieval_wavenumber, band_spec(:strong_co2))
    @test isempty(retrieval_extension.short)
    @test isempty(retrieval_extension.long)

    wavelength_file = normpath(joinpath(@__DIR__, "..", "..", "truth_map",
                                        "aerosol_chunked", "sim_wavelength.nc"))
    if isfile(wavelength_file)
        source_wavelength = NCDataset(wavelength_file) do dataset
            Float64.(dataset["strong_co2_wavelength"][:])
        end
        target = synthetic_grid(band_spec(:strong_co2))
        integrated = gaussian_convolve_resample(
            source_wavelength, ones(length(source_wavelength)), target,
            band_spec(:strong_co2).fwhm_nm)
        @test length(integrated) == length(target)
        @test all(isfinite, integrated)

        # The first 987 points reproduce the pre-merge basis and must still
        # fail the strict support guard without its eight-point shoulder.
        base_wavelength = source_wavelength[1:987]
        @test_throws ArgumentError gaussian_convolve_resample(
            base_wavelength, ones(length(base_wavelength)), target,
            band_spec(:strong_co2).fwhm_nm)

        source_wavenumber = NCDataset(wavelength_file) do dataset
            Float64.(dataset["strong_co2_wavenumber"][:])
        end
        integrated_extension = required_wavenumber_extensions(
            source_wavenumber, band_spec(:strong_co2))
        @test isempty(integrated_extension.short)
        @test isempty(integrated_extension.long)

        extension = required_wavenumber_extensions(
            source_wavenumber[1:987], band_spec(:strong_co2))
        @test length(extension.short) == 8
        @test isempty(extension.long)
        @test minimum(1e7 ./ extension.short) <= extension.required_short_nm
    end
end
