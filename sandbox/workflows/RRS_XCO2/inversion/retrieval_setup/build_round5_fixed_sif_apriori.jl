#!/usr/bin/env julia

module Round5FixedSIFPriorBuilder

"""
Construct the round-5 prior from the approved round-3 tapered prior.

Round 5 makes exactly three scientific changes:

1. both SIF coordinates are fixed and removed from the active state;
2. the UTLS/stratospheric aerosol `ln(AOD760)` and `ln(z0/km)` standard
   deviations are each divided by 10;
The CO2 mean and complete covariance are copied without modification so round
5 can be compared directly with round 4. A precision-matrix decorrelation is
reserved for a prospective round 6.
"""

using Dates
using LinearAlgebra
using NCDatasets
using Printf
using SHA

include(joinpath(@__DIR__, "..", "Round4SIFTruthConvention.jl"))
using .Round4SIFTruthConvention

export ROUND5_MODEL,
       ROUND5_ACTIVE_TO_FULL,
       ROUND5_ACTIVE_TO_CORE,
       STRATOSPHERIC_SIGMA_SCALE,
       default_output_directory,
       output_filename,
       file_sha256,
       build_round5_prior,
       write_round5_prior,
       write_round5_summary,
       main

const FULL_STATE_COUNT = 34
const ROUND5_ACTIVE_STATE_COUNT = 28
const CO2_FULL_INDICES = 2:17
const STRATOSPHERIC_AOD_FULL_INDEX = 20
const STRATOSPHERIC_Z0_FULL_INDEX = 23
const SIF760_FULL_INDEX = 33
const MSIF_FULL_INDEX = 34
const ROUND5_ACTIVE_TO_FULL = vcat(1, collect(6:32))
const ROUND5_ACTIVE_TO_CORE = collect(1:28)
const STRATOSPHERIC_SIGMA_SCALE = 0.1
const ROUND5_MODEL =
    "round5_fixed_sif759_msif_tight_utls_acos_mapped_tapered_vertical_correlation"
const DEFAULT_SOURCE_PRIOR = normpath(joinpath(
    @__DIR__, "..", "..", "bottom_layer_XCO2_retrievals",
    "retrieval_setup",
    "apriori_states_acos_mapped_tapered_vertical_correlation.nc"))
const DEFAULT_OUTPUT_DIRECTORY = normpath(joinpath(
    @__DIR__, "..", "..", "bottom_layer_XCO2_retrievals",
    "round5_fixed_sif", "retrieval_setup"))
const EXPECTED_SOURCE_ACTIVE = vcat(1, collect(6:34))
const EXPECTED_SOURCE_CO2_MODEL =
    "acos_mapped_tapered_vertical_correlation"
const SPECIES_ORDER = ("sulfate", "organic_carbon", "utls_sulfate")

default_output_directory() = DEFAULT_OUTPUT_DIRECTORY

file_sha256(path::AbstractString) = open(path, "r") do stream
    bytes2hex(sha256(stream))
end

function output_filename(mode::Symbol; extension::AbstractString="nc")
    mode in (:off, :on) || throw(ArgumentError(
        "round-5 SIF mode must be :off or :on"))
    extension in ("nc", "dat") || throw(ArgumentError(
        "round-5 prior extension must be nc or dat"))
    return "apriori_states_round5_fixed_sif_$(mode)_" *
           "tight_utls_acos_mapped_tapered_vertical_correlation.$extension"
end

function _attributes(attributes)
    Dict{String,Any}(String(key) => value for (key, value) in pairs(attributes))
end

function _read_source(path::AbstractString)
    isfile(path) || throw(ArgumentError(
        "missing approved round-3 tapered prior: $path"))
    sha = file_sha256(path)
    return NCDataset(path, "r") do dataset
        get(dataset.attrib, "apriori_complete", 0) == 1 || error(
            "source prior is not marked complete")
        String(get(dataset.attrib, "co2_covariance_model", "")) ==
            EXPECTED_SOURCE_CO2_MODEL || error(
                "round 5 requires the approved tapered round-3 CO2 prior")
        active = Int.(dataset["active_parameter_index"][:])
        active == EXPECTED_SOURCE_ACTIVE || error(
            "source prior does not have the expected 30-coordinate map")
        xa = Float64.(dataset["xa"][:, :])
        Sa = Float64.(dataset["Sa"][:, :, :])
        size(xa) == (FULL_STATE_COUNT, 4) || error(
            "source xa must be 34 by 4")
        size(Sa) == (FULL_STATE_COUNT, FULL_STATE_COUNT, 4) || error(
            "source Sa must be 34 by 34 by 4")
        return (;
            path=abspath(path), sha256=sha,
            attributes=_attributes(dataset.attrib), xa, Sa,
            surface_order=split(String(dataset.attrib["surface_order"])),
            parameter_names=split(String(dataset.attrib["parameter_names"])),
            parameter_units=split(
                String(dataset.attrib["parameter_units"]), " | "),
            layer_height=Float64.(dataset["co2_layer_center_height"][:]),
            layer_pressure=Float64.(dataset["co2_layer_center_pressure"][:]),
            taper_retention=Float64.(dataset[
                "co2_correlation_taper_adjacent_retention"][:]),
            base_adjacent=Float64.(dataset[
                "co2_base_adjacent_correlation"][:]),
            selected_adjacent=Float64.(dataset[
                "co2_selected_adjacent_correlation"][:]),
        )
    end
end

function _adjacent_correlations(covariance::AbstractMatrix)
    sigma = sqrt.(diag(covariance))
    return [covariance[index, index + 1] /
            (sigma[index] * sigma[index + 1])
            for index in 1:length(sigma) - 1]
end

function _fixed_sif(mode::Symbol)
    mode in (:off, :on) || throw(ArgumentError(
        "round-5 SIF mode must be :off or :on"))
    mode == :off && return (;
        Lnu759=0.0, Llambda759=0.0, SIF760=0.0, mSIF=0.0)
    convention = validate_round4_sif_truth_convention()
    Lnu759 = Float64(convention.diagnostic.Lnu)
    mSIF = Float64(convention.mSIF)
    nu759 = 1.0e7 / 759.0
    nu760 = 1.0e7 / 760.0
    SIF760 = Lnu759 + mSIF * (nu760 - nu759)
    return (;
        Lnu759,
        Llambda759=Float64(convention.diagnostic.Llambda),
        SIF760,
        mSIF,
    )
end

"""Build the four-surface round-5 prior and its complete audit record."""
function build_round5_prior(source_path::AbstractString, mode::Symbol)
    source = _read_source(source_path)
    fixed_sif = _fixed_sif(mode)
    xa = copy(source.xa)
    Sa = copy(source.Sa)

    for surface in axes(Sa, 3)
        # A congruence transform divides the two stratospheric coordinate
        # sigmas by 10 and would preserve correlations if any were added later.
        for index in (STRATOSPHERIC_AOD_FULL_INDEX,
                      STRATOSPHERIC_Z0_FULL_INDEX)
            Sa[index, :, surface] .*= STRATOSPHERIC_SIGMA_SCALE
            Sa[:, index, surface] .*= STRATOSPHERIC_SIGMA_SCALE
        end

        Sa[SIF760_FULL_INDEX:MSIF_FULL_INDEX, :, surface] .= 0.0
        Sa[:, SIF760_FULL_INDEX:MSIF_FULL_INDEX, surface] .= 0.0
        xa[SIF760_FULL_INDEX, surface] = fixed_sif.SIF760
        xa[MSIF_FULL_INDEX, surface] = fixed_sif.mSIF
    end

    active_mask = falses(FULL_STATE_COUNT)
    active_mask[ROUND5_ACTIVE_TO_FULL] .= true
    Sa_active = Array{Float64}(undef,
        ROUND5_ACTIVE_STATE_COUNT, ROUND5_ACTIVE_STATE_COUNT, 4)
    selected_adjacent = Vector{Float64}(undef, 11)
    for surface in axes(Sa, 3)
        active = Matrix(Symmetric(
            Sa[ROUND5_ACTIVE_TO_FULL, ROUND5_ACTIVE_TO_FULL, surface]))
        isposdef(Symmetric(active)) || error(
            "round-5 active covariance is not positive definite for surface $surface")
        Sa_active[:, :, surface] .= active
        Sa[CO2_FULL_INDICES, CO2_FULL_INDICES, surface] ==
            source.Sa[CO2_FULL_INDICES, CO2_FULL_INDICES, surface] || error(
                "round 5 changed the round-4 CO2 covariance")
        xa[CO2_FULL_INDICES, surface] ==
            source.xa[CO2_FULL_INDICES, surface] || error(
                "round 5 changed the round-4 CO2 mean")
        adjacent = copy(source.selected_adjacent)
        surface == 1 ? (selected_adjacent .= adjacent) :
            selected_adjacent == adjacent || error(
                "round-5 CO2 covariance differs between surfaces")
    end

    source_aod_sigma = sqrt.(diag(source.Sa[:, :, 1]))[18:20]
    source_z0_sigma = sqrt.(diag(source.Sa[:, :, 1]))[21:23]
    aerosol_aod_sigma = sqrt.(diag(Sa[:, :, 1]))[18:20]
    aerosol_z0_sigma = sqrt.(diag(Sa[:, :, 1]))[21:23]
    aerosol_aod_sigma[1:2] == source_aod_sigma[1:2] || error(
        "round 5 changed a tropospheric AOD prior")
    aerosol_z0_sigma[1:2] == source_z0_sigma[1:2] || error(
        "round 5 changed a tropospheric height prior")
    aerosol_aod_sigma[3] ==
        STRATOSPHERIC_SIGMA_SCALE * source_aod_sigma[3] || error(
            "round-5 stratospheric AOD sigma is not exactly 10x tighter")
    aerosol_z0_sigma[3] ==
        STRATOSPHERIC_SIGMA_SCALE * source_z0_sigma[3] || error(
            "round-5 stratospheric height sigma is not exactly 10x tighter")

    return (; mode, source, xa, Sa, Sa_active, active_mask,
            active_to_full=copy(ROUND5_ACTIVE_TO_FULL),
            active_to_core=copy(ROUND5_ACTIVE_TO_CORE), fixed_sif,
            source_aod_sigma, source_z0_sigma,
            aerosol_aod_sigma, aerosol_z0_sigma,
            selected_adjacent)
end

function write_round5_prior(prior;
                            output_path::AbstractString,
                            overwrite::Bool=false)
    isfile(output_path) && !overwrite && throw(ArgumentError(
        "round-5 prior already exists: $output_path"))
    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)

    NCDataset(output_path, "c") do output
        defDim(output, "parameter", FULL_STATE_COUNT)
        defDim(output, "parameter_2", FULL_STATE_COUNT)
        defDim(output, "active_parameter", ROUND5_ACTIVE_STATE_COUNT)
        defDim(output, "active_parameter_2", ROUND5_ACTIVE_STATE_COUNT)
        defDim(output, "surface", 4)
        defDim(output, "aerosol_species", 3)
        defDim(output, "co2_adjacent_layer_pair", 11)

        defVar(output, "xa", Float64,
               ("parameter", "surface"))[:, :] = prior.xa
        defVar(output, "Sa", Float64,
               ("parameter", "parameter_2", "surface"))[:, :, :] = prior.Sa
        defVar(output, "prior_sigma", Float64,
               ("parameter", "surface"))[:, :] =
            sqrt.(max.(0.0, [prior.Sa[index, index, surface]
                             for index in 1:FULL_STATE_COUNT,
                                 surface in 1:4]))
        defVar(output, "Sa_active", Float64,
               ("active_parameter", "active_parameter_2", "surface"))[:, :, :] =
            prior.Sa_active
        defVar(output, "active_parameter_index", Int16,
               ("active_parameter",))[:] = Int16.(prior.active_to_full)
        defVar(output, "active_core_parameter_index", Int16,
               ("active_parameter",))[:] = Int16.(prior.active_to_core)
        defVar(output, "active_mask", Int8,
               ("parameter",))[:] = Int8.(prior.active_mask)
        defVar(output, "co2_layer_center_height", Float64,
               ("parameter",))[:] = prior.source.layer_height
        defVar(output, "co2_layer_center_pressure", Float64,
               ("parameter",))[:] = prior.source.layer_pressure
        defVar(output, "co2_correlation_taper_adjacent_retention", Float64,
               ("co2_adjacent_layer_pair",))[:] =
            prior.source.taper_retention
        defVar(output, "co2_base_adjacent_correlation", Float64,
               ("co2_adjacent_layer_pair",))[:] =
            prior.source.base_adjacent
        defVar(output, "co2_selected_adjacent_correlation", Float64,
               ("co2_adjacent_layer_pair",))[:] = prior.selected_adjacent
        defVar(output, "aerosol_ln_aod_sigma_by_species", Float64,
               ("aerosol_species",))[:] = prior.aerosol_aod_sigma
        defVar(output, "aerosol_ln_z0_sigma_by_species", Float64,
               ("aerosol_species",))[:] = prior.aerosol_z0_sigma

        output.attrib["surface_order"] = join(prior.source.surface_order, " ")
        output.attrib["parameter_names"] =
            join(prior.source.parameter_names, " ")
        output.attrib["parameter_units"] =
            join(prior.source.parameter_units, " | ")
        output.attrib["state_order"] =
            "psurf; 16 CO2 TOA-to-BOA; 3 ln(AOD760); 3 ln(aerosol z0/km); 9 surface; fixed SIF760; fixed mSIF"
        output.attrib["retrieval_state_model"] = "round5_fixed_sif"
        output.attrib["round5_prior_model"] = ROUND5_MODEL
        output.attrib["round5_prior_definition_version"] = Int32(1)
        output.attrib["round5_sif_case"] = String(prior.mode)
        output.attrib["active_state_count"] = Int32(ROUND5_ACTIVE_STATE_COUNT)
        output.attrib["round5_active_to_full"] =
            join(prior.active_to_full, " ")
        output.attrib["round5_active_to_core"] =
            join(prior.active_to_core, " ")
        output.attrib["known_sif_wavelength_nm"] = 759.0
        output.attrib["known_sif_Lnu_mW_m-2_sr-1_per_cm-1"] =
            prior.fixed_sif.Lnu759
        output.attrib["known_sif_Llambda_mW_m-2_sr-1_nm-1"] =
            prior.fixed_sif.Llambda759
        output.attrib["fixed_SIF760_mW_m-2_sr-1_per_cm-1"] =
            prior.fixed_sif.SIF760
        output.attrib["fixed_mSIF_mW_m-2_sr-1_per_cm-2"] =
            prior.fixed_sif.mSIF
        output.attrib["sif_parameter_status"] =
            "SIF759 anchor and mSIF fixed exactly; both omitted from active state"
        output.attrib["co2_covariance_model"] =
            EXPECTED_SOURCE_CO2_MODEL
        output.attrib["co2_covariance_base_model"] =
            String(get(prior.source.attributes,
                       "co2_covariance_base_model", "acos_mapped"))
        output.attrib["co2_covariance_construction"] =
            String(prior.source.attributes["co2_covariance_construction"])
        output.attrib["round5_co2_prior_status"] =
            "xa[2:17] and Sa[2:17,2:17] copied exactly from source prior"
        output.attrib["stratospheric_aerosol_species"] = "utls_sulfate"
        output.attrib["stratospheric_aerosol_sigma_scale"] =
            STRATOSPHERIC_SIGMA_SCALE
        output.attrib["stratospheric_aod_sigma_source"] =
            prior.source_aod_sigma[3]
        output.attrib["stratospheric_aod_sigma_round5"] =
            prior.aerosol_aod_sigma[3]
        output.attrib["stratospheric_z0_sigma_source"] =
            prior.source_z0_sigma[3]
        output.attrib["stratospheric_z0_sigma_round5"] =
            prior.aerosol_z0_sigma[3]
        output.attrib["source_prior_path"] = prior.source.path
        output.attrib["source_prior_sha256"] = prior.source.sha256
        output.attrib["full_covariance_status"] =
            "positive semidefinite; four upper CO2 and both SIF coordinates fixed"
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["apriori_complete"] = Int32(1)
    end
    return abspath(output_path)
end

function write_round5_summary(prior;
                              output_path::AbstractString,
                              overwrite::Bool=false)
    isfile(output_path) && !overwrite && throw(ArgumentError(
        "round-5 prior summary already exists: $output_path"))
    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)
    open(output_path, "w") do io
        println(io, "# Round-5 fixed-SIF a priori: ", ROUND5_MODEL)
        println(io, "# SIF case: ", prior.mode)
        println(io, "# Source prior: ", prior.source.path)
        println(io, "# Source prior SHA-256: ", prior.source.sha256)
        println(io, "# Active coordinates: 28; full indices ",
                join(prior.active_to_full, " "))
        @printf(io, "# Fixed SIF: Lnu759=%.15e SIF760=%.15e mSIF=% .15e\n",
                prior.fixed_sif.Lnu759, prior.fixed_sif.SIF760,
                prior.fixed_sif.mSIF)
        println(io, "# CO2 prior: mean and full 16x16 covariance copied exactly from round 4.")
        @printf(io, "# ln(AOD760) sigma, sulfate/organic/UTLS: %.15e %.15e %.15e\n",
                prior.aerosol_aod_sigma...)
        @printf(io, "# ln(z0/km) sigma, sulfate/organic/UTLS: %.15e %.15e %.15e\n",
                prior.aerosol_z0_sigma...)
        println(io, "# CO2 adjacent correlations (unchanged in round 5)")
        for pair in eachindex(prior.selected_adjacent)
            @printf(io, "# layers %d-%d  %.15e\n",
                    pair + 4, pair + 5,
                    prior.selected_adjacent[pair])
        end
        println(io, "# surface full_index active_index parameter mean sigma units")
        active_position = Dict(index => position for
            (position, index) in enumerate(prior.active_to_full))
        for (surface_index, surface) in enumerate(prior.source.surface_order)
            for full_index in 1:FULL_STATE_COUNT
                @printf(io, "%-7s %2d %2d %-29s % .15e %.15e %s\n",
                        surface, full_index,
                        get(active_position, full_index, 0),
                        prior.source.parameter_names[full_index],
                        prior.xa[full_index, surface_index],
                        sqrt(max(prior.Sa[full_index, full_index,
                                          surface_index], 0.0)),
                        prior.source.parameter_units[full_index])
            end
        end
    end
    return abspath(output_path)
end

function main(args=ARGS)
    isempty(args) || error(
        "usage: build_round5_fixed_sif_apriori.jl (configure paths with environment variables)")
    source = get(ENV, "ROUND5_SOURCE_PRIOR", DEFAULT_SOURCE_PRIOR)
    output_directory = get(
        ENV, "ROUND5_PRIOR_OUTPUT_DIRECTORY", DEFAULT_OUTPUT_DIRECTORY)
    overwrite = get(ENV, "ROUND5_PRIOR_OVERWRITE", "0") == "1"
    for mode in (:off, :on)
        prior = build_round5_prior(source, mode)
        nc = joinpath(output_directory, output_filename(mode; extension="nc"))
        dat = joinpath(output_directory, output_filename(mode; extension="dat"))
        write_round5_prior(prior; output_path=nc, overwrite)
        write_round5_summary(prior; output_path=dat, overwrite)
        println(nc)
        println(dat)
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()

end # module Round5FixedSIFPriorBuilder
