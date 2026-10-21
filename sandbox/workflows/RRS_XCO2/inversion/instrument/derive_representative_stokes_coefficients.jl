#!/usr/bin/env julia

"""
Derive one sign-preserving OCO Stokes analyzer row per spectral band.

The OCO L1B field
`FootprintGeometry/footprint_stokes_coefficients` has Julia-order dimensions
`(coefficient, band, footprint, frame)`. For the present nadir truth study,
all available `oco2_L1bScND_*.h5` files in `OCO_L1B_ROOT` are pooled and every
valid sounding/footprint is given equal weight.

Direct component-wise means of M12 and M13 are not used: analyzer directions
rotate between soundings, so signed terms can cancel and produce an
unphysical weak polarizer. Instead, this script computes a circular median of
the analyzer direction atan(M13/M11, M12/M11), then selects the actual L1B
coefficient row nearest that median. The output therefore preserves signs and
the observed analyzer magnitude exactly.

Environment variables:

- `OCO_L1B_ROOT`: directory containing OCO L1B sounding files
  (default `/home/sanghavi/data/OCO2_jacobians`)
- `OCO_L1B_MODE`: `nadir` (default) or `all`
- `OCO_STOKES_OUT`: output NetCDF path
"""

using Dates
using NCDatasets
using Printf
using Statistics

include(joinpath(@__DIR__, "SyntheticOCO2.jl"))
using .SyntheticOCO2

const DEFAULT_L1B_ROOT = "/home/sanghavi/data/OCO2_jacobians"
const WORKFLOW_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RRS_ROOT = normpath(get(ENV, "RRS_XCO2_DATA_ROOT", WORKFLOW_ROOT))
const DEFAULT_OUTPUT = joinpath(
    RRS_ROOT, "inversion", "instrument", "representative_stokes_coefficients.nc")
const DATASET_PATH = "FootprintGeometry/footprint_stokes_coefficients"

wrap_angle(angle) = atan(sin(angle), cos(angle))

function selected_files(root::AbstractString, mode::AbstractString)
    isdir(root) || error("OCO L1B root does not exist: $root")
    normalized_mode = lowercase(mode)
    prefix = if normalized_mode == "nadir"
        "oco2_L1bScND_"
    elseif normalized_mode == "all"
        "oco2_L1bSc"
    else
        error("OCO_L1B_MODE must be 'nadir' or 'all'; received '$mode'")
    end
    files = sort(filter(name -> startswith(name, prefix) && endswith(name, ".h5"),
                        readdir(root)))
    isempty(files) && error("no matching OCO L1B files found in $root")
    return joinpath.(root, files)
end

function read_coefficients(path::AbstractString)
    NCDataset(path) do dataset
        group = dataset.group["FootprintGeometry"]
        haskey(group, "footprint_stokes_coefficients") ||
            error("$path does not contain $DATASET_PATH")
        variable = group["footprint_stokes_coefficients"]
        size(variable, 1) == 4 || error(
            "$path has $(size(variable, 1)) Stokes coefficients; expected 4")
        size(variable, 2) == length(BAND_SPECS) || error(
            "$path has $(size(variable, 2)) bands; expected $(length(BAND_SPECS))")
        return Array(variable[:, :, :, :])
    end
end

function collect_by_band(files)
    collected = [[Float64[] for _ in 1:4] for _ in BAND_SPECS]
    for path in files
        @info "reading OCO Stokes coefficients" path
        coefficients = read_coefficients(path)
        for band_index in eachindex(BAND_SPECS), coefficient_index in 1:4
            append!(collected[band_index][coefficient_index],
                    vec(Float64.(coefficients[coefficient_index, band_index, :, :])))
        end
    end
    return collected
end

function representative_row(columns)
    length(columns) == 4 || error("expected four coefficient columns")
    valid = isfinite.(columns[1]) .& isfinite.(columns[2]) .&
            isfinite.(columns[3]) .& isfinite.(columns[4]) .&
            (columns[1] .> 0)
    count(valid) > 0 || error("no valid OCO analyzer rows")
    coefficients = hcat((column[valid] for column in columns)...)

    angle = atan.(coefficients[:, 3] ./ coefficients[:, 1],
                  coefficients[:, 2] ./ coefficients[:, 1])
    mean_angle = atan(mean(sin.(angle)), mean(cos.(angle)))
    unwrapped = mean_angle .+ wrap_angle.(angle .- mean_angle)
    median_angle = wrap_angle(median(unwrapped))
    angular_distance = abs.(wrap_angle.(angle .- median_angle))
    representative_index = argmin(angular_distance)

    representative = vec(coefficients[representative_index, :])
    raw_mean = vec(mean(coefficients; dims=1))
    resultant_length = hypot(mean(cos.(angle)), mean(sin.(angle)))
    analyzer_degree = hypot(representative[2], representative[3]) /
                      representative[1]
    return (representative=representative,
            raw_mean=raw_mean,
            median_angle_deg=rad2deg(median_angle),
            resultant_length=resultant_length,
            analyzer_degree=analyzer_degree,
            sample_count=size(coefficients, 1),
            angular_residual_deg=rad2deg(angular_distance[representative_index]))
end

function write_netcdf(path, results, files, mode)
    mkpath(dirname(path))
    isfile(path) && rm(path)
    NCDataset(path, "c") do dataset
        defDim(dataset, "stokes_coefficient", 4)
        defDim(dataset, "band", length(BAND_SPECS))

        representative = defVar(dataset, "representative_stokes_coefficients",
                                Float64, ("stokes_coefficient", "band"))
        representative.attrib["long_name"] =
            "Observed OCO L1B analyzer row nearest the circular-median direction"
        representative.attrib["coefficient_order"] = "I Q U V"
        representative[:, :] = hcat((result.representative for result in results)...)

        raw_mean = defVar(dataset, "raw_signed_mean_stokes_coefficients",
                          Float64, ("stokes_coefficient", "band"))
        raw_mean.attrib["warning"] =
            "Diagnostic only; rotating signed Q/U terms cancel in this mean"
        raw_mean[:, :] = hcat((result.raw_mean for result in results)...)

        for (name, values, description) in (
            ("circular_median_angle", [r.median_angle_deg for r in results],
             "Circular median of atan(M13/M11,M12/M11)"),
            ("circular_resultant_length", [r.resultant_length for r in results],
             "Magnitude of the raw mean normalized analyzer direction"),
            ("representative_analyzer_degree", [r.analyzer_degree for r in results],
             "hypot(M12,M13)/M11 for the selected observed row"),
            ("representative_angular_residual", [r.angular_residual_deg for r in results],
             "Angular distance between selected row and circular median"),
        )
            variable = defVar(dataset, name, Float64, ("band",))
            variable.attrib["long_name"] = description
            if occursin("angle", name) || occursin("residual", name)
                variable.attrib["units"] = "degree"
            end
            variable[:] = values
        end

        sample_count = defVar(dataset, "sample_count", Int64, ("band",))
        sample_count[:] = [r.sample_count for r in results]

        dataset.attrib["band_order"] = join(string.(getfield.(BAND_SPECS, :name)), " ")
        dataset.attrib["band_long_names"] = join(getfield.(BAND_SPECS, :long_name), "; ")
        dataset.attrib["source_dataset"] = DATASET_PATH
        dataset.attrib["source_files"] = join(basename.(files), ";")
        dataset.attrib["source_scope"] = mode
        dataset.attrib["estimator"] =
            "circular median of analyzer direction followed by nearest observed row"
        dataset.attrib["reason_for_estimator"] =
            "preserves signed Q/U direction and analyzer magnitude; avoids component cancellation"
        dataset.attrib["projection_convention"] = "M11*I - M12*Q + M13*U"
        dataset.attrib["normalization"] =
            "none beyond the OCO coefficients themselves; matches oco_gain.jl"
        dataset.attrib["created_utc"] = string(now(UTC))
    end
end

function write_table(path, results, files, mode)
    open(path, "w") do io
        println(io, "# Representative OCO L1B Stokes analyzer rows")
        println(io, "# source_scope = $mode")
        println(io, "# estimator = circular median direction; nearest observed row")
        println(io, "# projection = M11*I - M12*Q + M13*U; no subsequent normalization")
        println(io, "# source_files = $(join(basename.(files), ';'))")
        println(io, "# band M11 M12 M13 M14 circular_median_angle_deg resultant_length sample_count")
        for (spec, result) in zip(BAND_SPECS, results)
            @printf(io, "%s %.12g %.12g %.12g %.12g %.9f %.9f %d\n",
                    spec.name, result.representative...,
                    result.median_angle_deg, result.resultant_length,
                    result.sample_count)
        end
    end
end

function main()
    root = get(ENV, "OCO_L1B_ROOT", DEFAULT_L1B_ROOT)
    mode = lowercase(get(ENV, "OCO_L1B_MODE", "nadir"))
    output = get(ENV, "OCO_STOKES_OUT", DEFAULT_OUTPUT)
    files = selected_files(root, mode)
    collected = collect_by_band(files)
    results = representative_row.(collected)

    for (spec, result) in zip(BAND_SPECS, results)
        @info "representative OCO analyzer" band=spec.name coefficients=result.representative circular_median_angle_deg=result.median_angle_deg raw_mean=result.raw_mean resultant_length=result.resultant_length sample_count=result.sample_count
    end

    write_netcdf(output, results, files, mode)
    table_path = splitext(output)[1] * ".dat"
    write_table(table_path, results, files, mode)
    @info "wrote representative OCO analyzer products" output table_path
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
