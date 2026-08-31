#!/usr/bin/env julia

"""
Derive wavelength-dependent representative OCO-2 noise coefficients.

The two `InstrumentHeader/snr_coef` spectra are read for every footprint in
the four nadir L1B files used by `oco_gain.jl`. Each coefficient spectrum is
mapped to the synthetic OCO wavelength grid using its own dispersion
polynomial. The component-wise median across the resulting 32 spectra is used
because both coefficients are non-negative and do not have the signed
cancellation issue of the Stokes analyzer terms.

Samples outside a source footprint's dispersion range use its nearest edge
coefficient. The output records how many of the 32 source spectra required
this extrapolation at every synthetic sample. This currently affects only the
short-wavelength edge of the synthetic strong-CO2 grid.

Environment variables:

- `OCO_L1B_ROOT`: L1B directory (default `/home/sanghavi/data/OCO2_jacobians`)
- `OCO_SNR_COEFFICIENTS_OUT`: output NetCDF path
"""

using Dates
using NCDatasets
using Printf
using Statistics

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
include(joinpath(@__DIR__, "OCO2Noise.jl"))
using .SyntheticOCO2
using .OCO2Noise

const DEFAULT_L1B_ROOT = "/home/sanghavi/data/OCO2_jacobians"
const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RRS_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const DEFAULT_OUTPUT = joinpath(
    RRS_ROOT, "inversion", "instrument", "representative_snr_coefficients.nc")
const ATBD_PATH = normpath(joinpath(
    @__DIR__, "..", "..", "resources", "OCO_L1B_ATBD.pdf"))
const L1B_FILENAMES = (
    "oco2_L1bScND_37214a_210630_B11006r_220326081005.h5",
    "oco2_L1bScND_42610a_220706_B11008r_220831171738.h5",
    "oco2_L1bScND_41856a_220515_B11008r_220804213809.h5",
    "oco2_L1bScND_45442a_230116_B11008r_230209231005.h5",
)

function dispersion_wavelength_nm(coefficients::AbstractVector)
    length(coefficients) == 6 || error(
        "expected six OCO dispersion coefficients; found $(length(coefficients))")
    sample = collect(1.0:1016.0)
    wavelength_um = zeros(Float64, length(sample))
    for power in 0:5
        wavelength_um .+= Float64(coefficients[power + 1]) .* sample .^ power
    end
    all(diff(wavelength_um) .> 0) || error(
        "OCO dispersion polynomial is not strictly increasing")
    return 1000 .* wavelength_um
end

function interpolate_edge(source_x::AbstractVector,
                          source_y::AbstractVector,
                          target_x::AbstractVector)
    length(source_x) == length(source_y) ||
        throw(DimensionMismatch("source coordinate and value lengths differ"))
    all(diff(source_x) .> 0) ||
        throw(ArgumentError("source coordinate must be strictly increasing"))
    values = Vector{Float64}(undef, length(target_x))
    extrapolated = falses(length(target_x))
    for (target_index, target) in pairs(target_x)
        if target <= first(source_x)
            values[target_index] = source_y[1]
            extrapolated[target_index] = target < first(source_x)
        elseif target >= last(source_x)
            values[target_index] = source_y[end]
            extrapolated[target_index] = target > last(source_x)
        else
            left = searchsortedlast(source_x, target)
            fraction = (target - source_x[left]) /
                       (source_x[left + 1] - source_x[left])
            values[target_index] = (1 - fraction) * source_y[left] +
                                   fraction * source_y[left + 1]
        end
    end
    return values, extrapolated
end

function read_l1b_coefficients(path::AbstractString)
    NCDataset(path) do dataset
        group = dataset.group["InstrumentHeader"]
        snr = Float64.(Array(group["snr_coef"][:, :, :, :]))
        dispersion = Float64.(Array(group["dispersion_coef_samp"][:, :, :]))
        size(snr) == (2, 1016, 8, 3) || error(
            "$path has unexpected snr_coef size $(size(snr))")
        size(dispersion) == (6, 8, 3) || error(
            "$path has unexpected dispersion size $(size(dispersion))")
        return snr, dispersion
    end
end

function derive(files)
    photon_columns = Dict(spec.name => Vector{Vector{Float64}}()
                          for spec in BAND_SPECS)
    background_columns = Dict(spec.name => Vector{Vector{Float64}}()
                              for spec in BAND_SPECS)
    extrapolated_columns = Dict(spec.name => Vector{BitVector}()
                                for spec in BAND_SPECS)

    for path in files
        @info "reading OCO SNR coefficients" path
        snr, dispersion = read_l1b_coefficients(path)
        for (band_index, spec) in pairs(BAND_SPECS)
            grid = synthetic_grid(spec)
            for footprint in 1:8
                source_wavelength = dispersion_wavelength_nm(
                    @view dispersion[:, footprint, band_index])
                photon, extrapolated = interpolate_edge(
                    source_wavelength,
                    @view(snr[1, :, footprint, band_index]), grid)
                background, extrapolated_background = interpolate_edge(
                    source_wavelength,
                    @view(snr[2, :, footprint, band_index]), grid)
                extrapolated == extrapolated_background || error(
                    "coefficient extrapolation masks disagree")
                push!(photon_columns[spec.name], photon)
                push!(background_columns[spec.name], background)
                push!(extrapolated_columns[spec.name], extrapolated)
            end
        end
    end

    return Dict(spec.name => begin
        photon = hcat(photon_columns[spec.name]...)
        background = hcat(background_columns[spec.name]...)
        extrapolated = hcat(extrapolated_columns[spec.name]...)
        sample_count = size(photon, 2)
        (
            wavelength=synthetic_grid(spec),
            c_photon=[median(@view photon[index, :]) for index in axes(photon, 1)],
            c_photon_q05=[quantile(@view(photon[index, :]), 0.05)
                          for index in axes(photon, 1)],
            c_photon_q95=[quantile(@view(photon[index, :]), 0.95)
                          for index in axes(photon, 1)],
            c_background=[median(@view background[index, :])
                          for index in axes(background, 1)],
            c_background_q05=[quantile(@view(background[index, :]), 0.05)
                              for index in axes(background, 1)],
            c_background_q95=[quantile(@view(background[index, :]), 0.95)
                              for index in axes(background, 1)],
            extrapolated_source_count=Int16.(vec(sum(extrapolated; dims=2))),
            sample_count=sample_count,
        )
    end for spec in BAND_SPECS)
end

function define_coefficient(output, name, dimension, values, long_name)
    variable = defVar(output, name, Float64, (dimension,))
    variable.attrib["long_name"] = long_name
    variable.attrib["units"] = "1"
    variable[:] = values
end

function write_netcdf(path, results, files)
    mkpath(dirname(path))
    isfile(path) && rm(path)
    NCDataset(path, "c") do output
        defDim(output, "band", length(BAND_SPECS))
        max_ms = defVar(output, "max_ms", Float64, ("band",))
        max_ms.attrib["units"] = "photons s-1 m-2 sr-1 um-1"
        max_ms.attrib["source"] = "OCO L1B ATBD Table 3-5"
        max_ms[:] = [noise_band_spec(spec.name).max_ms for spec in BAND_SPECS]
        min_ms = defVar(output, "min_ms", Float64, ("band",))
        min_ms.attrib["units"] = "photons s-1 m-2 sr-1 um-1"
        min_ms.attrib["source"] = "OCO L1B ATBD Table 3-6, OCO-2 column"
        min_ms.attrib["usage"] = "dynamic-range validation only; not an Eq. 3-8 floor"
        min_ms[:] = [noise_band_spec(spec.name).min_ms for spec in BAND_SPECS]

        for spec in BAND_SPECS
            result = results[spec.name]
            name = String(spec.name)
            dimension = name
            defDim(output, dimension, length(result.wavelength))
            wavelength = defVar(output, "$(name)_wavelength", Float64, (dimension,))
            wavelength.attrib["units"] = "nm"
            wavelength.attrib["long_name"] =
                "synthetic OCO sample-center wavelength for $(spec.long_name)"
            wavelength[:] = result.wavelength

            define_coefficient(output, "c_photon_$(name)", dimension,
                               result.c_photon,
                               "median OCO L1B photon-noise coefficient")
            define_coefficient(output, "c_photon_q05_$(name)", dimension,
                               result.c_photon_q05,
                               "5th percentile OCO L1B photon-noise coefficient")
            define_coefficient(output, "c_photon_q95_$(name)", dimension,
                               result.c_photon_q95,
                               "95th percentile OCO L1B photon-noise coefficient")
            define_coefficient(output, "c_background_$(name)", dimension,
                               result.c_background,
                               "median OCO L1B background-noise coefficient")
            define_coefficient(output, "c_background_q05_$(name)", dimension,
                               result.c_background_q05,
                               "5th percentile OCO L1B background-noise coefficient")
            define_coefficient(output, "c_background_q95_$(name)", dimension,
                               result.c_background_q95,
                               "95th percentile OCO L1B background-noise coefficient")

            extrapolated = defVar(
                output, "extrapolated_source_count_$(name)", Int16, (dimension,))
            extrapolated.attrib["long_name"] =
                "number of 32 source coefficient spectra using nearest-edge extrapolation"
            extrapolated.attrib["valid_range"] = Int16[0, result.sample_count]
            extrapolated[:] = result.extrapolated_source_count
        end

        output.attrib["band_order"] = join(string.(getfield.(BAND_SPECS, :name)), " ")
        output.attrib["source_dataset"] = "InstrumentHeader/snr_coef"
        output.attrib["source_dispersion_dataset"] =
            "InstrumentHeader/dispersion_coef_samp"
        output.attrib["source_files"] = join(basename.(files), ";")
        output.attrib["source_scope"] = "four nadir OCO-2 L1B files used by oco_gain.jl"
        output.attrib["estimator"] =
            "component-wise median after wavelength interpolation; 32 spectra per band"
        output.attrib["edge_policy"] =
            "nearest source coefficient; per-sample source count is stored"
        output.attrib["noise_equation"] =
            "OCO-2 Level 1B ATBD v3.0 rev0 Eq. 3-8"
        output.attrib["atbd_path"] = abspath(ATBD_PATH)
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["representative_snr_coefficients_complete"] = 1
    end
end

function write_table(path, results, files)
    open(path, "w") do io
        println(io, "# Representative OCO-2 L1B SNR coefficients")
        println(io, "# estimator = pointwise median of 4 files x 8 footprints")
        println(io, "# source_files = $(join(basename.(files), ';'))")
        println(io, "# MaxMS/MinMS = ATBD Tables 3-5/3-6; MinMS is validation only")
        println(io, "# band MaxMS MinMS Cphoton_min Cphoton_max Cbackground_min Cbackground_max any_extrapolated_samples all_extrapolated_samples source_spectra")
        for spec in BAND_SPECS
            result = results[spec.name]
            noise_spec = noise_band_spec(spec.name)
            any_extrapolated = count(>(0), result.extrapolated_source_count)
            all_extrapolated = count(==(result.sample_count),
                                     result.extrapolated_source_count)
            @printf(io, "%s %.9g %.9g %.12g %.12g %.12g %.12g %d %d %d\n",
                    spec.name, noise_spec.max_ms, noise_spec.min_ms,
                    extrema(result.c_photon)...,
                    extrema(result.c_background)...,
                    any_extrapolated, all_extrapolated, result.sample_count)
        end
    end
end

function main()
    root = get(ENV, "OCO_L1B_ROOT", DEFAULT_L1B_ROOT)
    files = joinpath.(root, L1B_FILENAMES)
    all(isfile, files) || error(
        "one or more required OCO L1B files are missing from $root")
    isfile(ATBD_PATH) || error("missing OCO L1B ATBD: $ATBD_PATH")
    output = get(ENV, "OCO_SNR_COEFFICIENTS_OUT", DEFAULT_OUTPUT)
    results = derive(files)
    write_netcdf(output, results, files)
    table_path = splitext(output)[1] * ".dat"
    write_table(table_path, results, files)
    for spec in BAND_SPECS
        result = results[spec.name]
        @info "representative OCO SNR coefficients" band=spec.name c_photon=extrema(result.c_photon) c_background=extrema(result.c_background) extrapolated_samples=count(>(0), result.extrapolated_source_count) fully_extrapolated_samples=count(==(result.sample_count), result.extrapolated_source_count)
    end
    @info "wrote representative OCO SNR products" output table_path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
