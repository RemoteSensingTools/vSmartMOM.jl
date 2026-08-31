#!/usr/bin/env julia

"""
Create synthetic OCO-2 measurement files from completed high-resolution truth
scenes.

For each selected band and radiative-transfer component, the processing order
is fixed and linear:

1. project `(I,Q,U)` through the representative OCO Stokes analyzer row;
2. retain the raw OCO analyzer response (no additional normalization);
3. convert radiance density from per cm^-1 to per nm;
4. convolve with the Table-1 Gaussian FWHM in wavelength space;
5. evaluate at the requested synthetic OCO sample centers.

The O2 A-band output retains processed Rayleigh, Cabannes, and RRS components
for quality control. The corrected measurement is Rayleigh; the uncorrected
measurement is Cabannes + RRS. Both CO2 measurement classes are identical and
use the truth-map Rayleigh simulations.

Usage:

    julia --project=. RRS_XCO2/inversion/instrument/process_truth_map.jl

Environment variables:

- `TRUTH_INPUT_DIRS`: colon-separated scene directories. Defaults to the
  top-level truth-map directory and its `aerosol_chunked` subdirectory.
- `SYNTHETIC_OCO_OUT`: output directory (default
  `RRS_XCO2/truth_map/OCO_radiances`).
- `OCO_STOKES_COEFFICIENTS`: representative coefficient NetCDF.
- `SUPPLEMENTAL_SHOULDER_DIR`: optional legacy/staging convolution-shoulder
  directory (default `RRS_XCO2/truth_map/convolution_shoulders`). Merged truth
  files already contain these samples and need no auxiliary directory.
- `BANDS`: comma-separated subset of `o2a,weak_co2,strong_co2`.
- `FIRST_STATE`, `LAST_STATE`: inclusive scene-index limits.
- `ALLOW_INCOMPLETE=1`: development-only override of the completed-scene
  attribute check. It does not bypass fill-value or Gaussian-support checks.
- `ALLOW_CLIPPED_GAUSSIAN=1`: development-only edge-kernel override. Never
  use this for production; generate supplemental shoulders instead.
- `FORCE=1`: replace existing synthetic output files.
"""

using Dates
using NCDatasets
using Printf

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
using .SyntheticOCO2

const INVERSION_ROOT = normpath(joinpath(@__DIR__, ".."))
const RRS_ROOT = normpath(joinpath(INVERSION_ROOT, ".."))
const DEFAULT_TRUTH_ROOT = joinpath(RRS_ROOT, "truth_map")
const DEFAULT_OUTPUT = joinpath(DEFAULT_TRUTH_ROOT, "OCO_radiances")
const DEFAULT_COEFFICIENTS = joinpath(@__DIR__, "representative_stokes_coefficients.nc")
const DEFAULT_SUPPLEMENTAL_SHOULDERS = joinpath(DEFAULT_TRUTH_ROOT, "convolution_shoulders")
const SUPPORT_SIGMA = 6.0

env_flag(name) = lowercase(get(ENV, name, "0")) in ("1", "true", "yes", "on")

function selected_bands()
    requested = Symbol.(strip.(split(get(ENV, "BANDS", "o2a,weak_co2,strong_co2"), ',')))
    isempty(requested) && error("BANDS selected no bands")
    for name in requested
        band_spec(name)
    end
    return requested
end

function state_index(path)
    match_result = match(r"hiressim_(\d+)\.nc$", basename(path))
    isnothing(match_result) && error("cannot parse state index from $path")
    return parse(Int, only(match_result.captures))
end

function discover_scenes(input_directories, first_state, last_state)
    scenes = Dict{Int,String}()
    for directory in input_directories
        isdir(directory) || continue
        for name in sort(readdir(directory))
            occursin(r"^hiressim_\d+\.nc$", name) || continue
            path = joinpath(directory, name)
            index = state_index(path)
            first_state <= index <= last_state || continue
            haskey(scenes, index) && error(
                "duplicate truth scene $index in $(scenes[index]) and $path")
            scenes[index] = path
        end
    end
    isempty(scenes) && error("no truth scenes found in $(join(input_directories, ", "))")
    return sort!(collect(scenes); by=first)
end

function wavelength_path_for(scene_path)
    local_path = joinpath(dirname(scene_path), "sim_wavelength.nc")
    isfile(local_path) || error("missing wavelength file beside scene: $local_path")
    return local_path
end

function read_wavelengths(path, bands)
    NCDataset(path) do dataset
        Dict(name => Float64.(nomissing(dataset["$(name)_wavelength"][:], NaN))
             for name in bands)
    end
end

function read_stokes(dataset, variable_name)
    haskey(dataset, variable_name) || error(
        "truth scene does not contain $variable_name")
    values = Array(dataset[variable_name][:, :])
    stokes = Float64.(nomissing(values, NaN))
    size(stokes, 1) == 3 || error(
        "$variable_name has $(size(stokes, 1)) Stokes rows; expected 3")
    all(isfinite, stokes) || error(
        "$variable_name contains missing, fill, or non-finite values")
    # NetCDF's default Float32 fill value is finite, so reject it explicitly.
    maximum(abs, stokes) < 1e30 || error(
        "$variable_name still contains unwritten NetCDF fill values")
    return stokes
end

function append_supplemental_shoulders(source_wavelength,
                                       source_stokes,
                                       supplemental_directory,
                                       index,
                                       band::Symbol,
                                       component::Symbol)
    path = joinpath(supplemental_directory,
                    @sprintf("hires_shoulders_%03d.nc", index))
    isfile(path) || return source_wavelength, source_stokes

    wavelength_parts = Vector{Float64}[Float64.(source_wavelength)]
    stokes_parts = Matrix{Float64}[Float64.(source_stokes)]
    NCDataset(path) do dataset
        if haskey(dataset.attrib, "supplemental_convolution_shoulders_complete")
            Int(dataset.attrib["supplemental_convolution_shoulders_complete"]) == 1 ||
                error("supplemental shoulder file is not complete: $path")
        end
        for side in (:short, :long)
            stem = "$(band)_$(side)_shoulder"
            variable_name = "radiance_$(component)_$(stem)"
            haskey(dataset, variable_name) || continue
            wavelength_name = "$(stem)_wavelength"
            haskey(dataset, wavelength_name) || error(
                "$path contains $variable_name but no $wavelength_name")
            supplemental_wavelength =
                Float64.(nomissing(dataset[wavelength_name][:], NaN))
            supplemental_stokes = read_stokes(dataset, variable_name)
            matches = map(supplemental_wavelength) do wavelength
                findfirst(value -> isapprox(value, wavelength; atol=1e-10, rtol=0),
                          source_wavelength)
            end
            if all(index -> !isnothing(index), matches)
                indices = Int[something(index) for index in matches]
                source_stokes[:, indices] == supplemental_stokes || error(
                    "$path duplicates embedded shoulder wavelengths but its " *
                    "radiances differ from the truth scene")
                continue
            elseif any(index -> !isnothing(index), matches)
                error("$path partially overlaps the embedded truth grid")
            end
            push!(wavelength_parts, supplemental_wavelength)
            push!(stokes_parts, supplemental_stokes)
        end
    end
    return vcat(wavelength_parts...), hcat(stokes_parts...)
end

function require_completed_scene(dataset, allow_incomplete)
    allow_incomplete && return
    generic_complete = haskey(dataset.attrib, "simulation_complete") &&
        Int(dataset.attrib["simulation_complete"]) == 1
    chunked_complete = haskey(dataset.attrib, "chunked_simulation_complete") &&
        Int(dataset.attrib["chunked_simulation_complete"]) == 1
    generic_complete || chunked_complete || error(
        "scene is not marked simulation_complete=1 or " *
        "chunked_simulation_complete=1")
end

function process_scene(scene_path, output_path, coefficient_path,
                       supplemental_directory, bands;
                       allow_incomplete=false,
                       require_full_support=true)
    coefficients = read_representative_coefficients(coefficient_path)
    wavelengths = read_wavelengths(wavelength_path_for(scene_path), bands)
    processed = Dict{Symbol,Dict{Symbol,Vector{Float64}}}()
    source_attributes = Dict{String,Any}()
    index = state_index(scene_path)

    NCDataset(scene_path) do dataset
        require_completed_scene(dataset, allow_incomplete)
        for key in keys(dataset.attrib)
            source_attributes[String(key)] = dataset.attrib[key]
        end

        if :o2a in bands
            spec = band_spec(:o2a)
            lambda = wavelengths[:o2a]
            rayleigh_lambda, rayleigh_stokes = append_supplemental_shoulders(
                lambda, read_stokes(dataset, "radiance_rayleigh_o2a"),
                supplemental_directory, index, :o2a, :rayleigh)
            cabannes_lambda, cabannes_stokes = append_supplemental_shoulders(
                lambda, read_stokes(dataset, "radiance_cabannes_o2a"),
                supplemental_directory, index, :o2a, :cabannes)
            rrs_lambda, rrs_stokes = append_supplemental_shoulders(
                lambda, read_stokes(dataset, "radiance_rrs_o2a"),
                supplemental_directory, index, :o2a, :rrs)
            rayleigh = process_stokes_spectrum(
                rayleigh_lambda, rayleigh_stokes,
                coefficients[:o2a], spec;
                support_sigma=SUPPORT_SIGMA,
                require_full_support=require_full_support)
            cabannes = process_stokes_spectrum(
                cabannes_lambda, cabannes_stokes,
                coefficients[:o2a], spec;
                support_sigma=SUPPORT_SIGMA,
                require_full_support=require_full_support)
            rrs = process_stokes_spectrum(
                rrs_lambda, rrs_stokes,
                coefficients[:o2a], spec;
                support_sigma=SUPPORT_SIGMA,
                require_full_support=require_full_support)
            processed[:o2a] = Dict(
                :rayleigh => rayleigh,
                :cabannes => cabannes,
                :rrs => rrs,
                :corrected => rayleigh,
                :uncorrected => cabannes .+ rrs,
            )
        end

        for name in (:weak_co2, :strong_co2)
            name in bands || continue
            spec = band_spec(name)
            source_wavelength, source_stokes = append_supplemental_shoulders(
                wavelengths[name],
                read_stokes(dataset, "radiance_rayleigh_$(name)"),
                supplemental_directory, index, name, :rayleigh)
            rayleigh = process_stokes_spectrum(
                source_wavelength,
                source_stokes,
                coefficients[name], spec;
                support_sigma=SUPPORT_SIGMA,
                require_full_support=require_full_support)
            processed[name] = Dict(
                :rayleigh => rayleigh,
                :corrected => rayleigh,
                :uncorrected => rayleigh,
            )
        end
    end

    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)
    NCDataset(output_path, "c") do output
        for name in bands
            spec = band_spec(name)
            grid = synthetic_grid(spec)
            dimension_name = String(name)
            defDim(output, dimension_name, length(grid))
            wavelength = defVar(output, "$(name)_wavelength", Float64, (dimension_name,))
            wavelength.attrib["units"] = "nm"
            wavelength.attrib["long_name"] = "synthetic OCO-2 sample-center wavelength"
            wavelength.attrib["sampling_interval_nm"] = spec.sampling_interval_nm
            wavelength.attrib["gaussian_fwhm_nm"] = spec.fwhm_nm
            wavelength.attrib["gaussian_sigma_nm"] = fwhm_to_sigma(spec.fwhm_nm)
            wavelength[:] = grid

            for (component, values) in sort!(collect(processed[name]); by=first)
                variable_name = "I_OCO_$(component)_$(name)"
                variable = defVar(output, variable_name, Float64, (dimension_name,))
                variable.attrib["units"] = "mW m-2 sr-1 nm-1"
                variable.attrib["long_name"] =
                    "Mueller-analyzer processed OCO intensity: $(component), $(spec.long_name)"
                variable.attrib["source_component"] = String(component)
                variable.attrib["instrument_processing"] =
                    "raw OCO analyzer; per-cm-1 to per-nm; Gaussian convolution; resampling"
                variable[:] = values
            end
        end

        for key in ("state_index", "surface", "aerosol_case", "sif_case",
                    "xco2_ppm", "psurf_hpa", "sza_deg", "vza_deg",
                    "atmospheric_layers", "spectroscopy_database",
                    "spectroscopy_version", "o2_absco_lut",
                    "o2_h2o_lut", "o2_truth_reused",
                    "o2_truth_reuse_source", "o2_truth_spectroscopy",
                    "o2_truth_grid_note", "co2_absco_regeneration_complete",
                    "weak_h2o_absco_lut", "weak_co2_absco_lut",
                    "strong_h2o_absco_lut", "strong_co2_absco_lut",
                    "h2o_line_absorption_by_band", "o2_vmr",
                    "atmospheric_profile_configuration",
                    "profile_preparation")
            haskey(source_attributes, key) &&
                (output.attrib[key] = source_attributes[key])
        end
        output.attrib["source_truth_scene"] = abspath(scene_path)
        output.attrib["representative_stokes_coefficients"] = abspath(coefficient_path)
        shoulder_provenance = if get(
            source_attributes, "strong_co2_short_shoulder_merged", 0) == 1
            "embedded in source truth strong_co2 grid"
        elseif isdir(supplemental_directory)
            abspath(supplemental_directory)
        else
            "none"
        end
        output.attrib["supplemental_convolution_shoulders"] = shoulder_provenance
        output.attrib["analyzer_projection"] = "M11*I - M12*Q + M13*U"
        output.attrib["analyzer_normalization"] =
            "none beyond the coefficients themselves; matches oco_gain.jl"
        output.attrib["spectral_density_conversion"] =
            "L_per_nm = L_per_cm-1 * 1e7/lambda_nm^2"
        output.attrib["ils_shape"] = "Gaussian in wavelength"
        output.attrib["gaussian_support_sigma"] = SUPPORT_SIGMA
        output.attrib["band_definitions"] =
            "o2a 758:0.015:772 nm FWHM=0.04 nm; " *
            "weak_co2 1594:0.031:1619 nm FWHM=0.08 nm; " *
            "strong_co2 2042:0.04:2082 nm FWHM=0.10 nm"
        output.attrib["corrected_measurement"] = "Rayleigh"
        output.attrib["uncorrected_measurement"] =
            "O2A: Cabannes+RRS; weak/strong CO2: Rayleigh"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["instrument_processing_complete"] = 1
    end
    return output_path
end

function main()
    input_default = join((DEFAULT_TRUTH_ROOT,
                          joinpath(DEFAULT_TRUTH_ROOT, "aerosol_chunked")), ':')
    input_directories = filter(path -> !isempty(path),
                               split(get(ENV, "TRUTH_INPUT_DIRS", input_default), ':'))
    output_directory = get(ENV, "SYNTHETIC_OCO_OUT", DEFAULT_OUTPUT)
    coefficient_path = get(ENV, "OCO_STOKES_COEFFICIENTS", DEFAULT_COEFFICIENTS)
    supplemental_directory = get(
        ENV, "SUPPLEMENTAL_SHOULDER_DIR", DEFAULT_SUPPLEMENTAL_SHOULDERS)
    bands = selected_bands()
    first_state = parse(Int, get(ENV, "FIRST_STATE", "1"))
    last_state = parse(Int, get(ENV, "LAST_STATE", "64"))
    allow_incomplete = env_flag("ALLOW_INCOMPLETE")
    require_full_support = !env_flag("ALLOW_CLIPPED_GAUSSIAN")
    force = env_flag("FORCE")

    scenes = discover_scenes(input_directories, first_state, last_state)
    mkpath(output_directory)
    for (index, scene_path) in scenes
        output_path = joinpath(output_directory, @sprintf("OCO2sims_%03d.nc", index))
        if isfile(output_path) && !force
            @info "skip existing synthetic OCO scene" index output_path
            continue
        end
        @info "process synthetic OCO scene" index scene_path bands
        process_scene(scene_path, output_path, coefficient_path,
                      supplemental_directory, bands;
                      allow_incomplete=allow_incomplete,
                      require_full_support=require_full_support)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
