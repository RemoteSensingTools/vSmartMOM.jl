#!/usr/bin/env julia

"""
Validate and catalogue the exact Wurst aerosol retrieval partition for the
bottom-layer CO2 campaign.

This preflight is intentionally separate from the radiance and noise-product
validators. The launcher runs those first for the complete 80-state dataset,
then uses this script to verify the retrieval-specific truth selection,
campaign-local paths, a-priori covariance, and execution order. No retrieval
is started by this script.
"""

using LinearAlgebra
using NCDatasets

include(joinpath(@__DIR__, "RetrievalCases.jl"))
using .RetrievalCases

const DEFAULT_CAMPAIGN_ROOT = normpath(joinpath(
    @__DIR__, "..", "bottom_layer_XCO2_retrievals"))
const NOSIF_AEROSOL_STATES = vcat(collect(11:15), collect(31:35),
                                  collect(51:55), collect(71:75))
const SIF_AEROSOL_STATES = vcat(collect(16:20), collect(36:40),
                                collect(56:60), collect(76:80))
const ORDERED_AEROSOL_STATES = vcat(NOSIF_AEROSOL_STATES,
                                    SIF_AEROSOL_STATES)
const AOD_PARAMETER_INDICES = collect(18:20)
const EXPECTED_AEROSOL_LN_AOD_SIGMA = 0.75
const EXPECTED_AOD_VARIANCE = EXPECTED_AEROSOL_LN_AOD_SIGMA^2

function validate_prior(path)
    isfile(path) || error("missing campaign-local a-priori file: $path")
    NCDataset(path) do dataset
        get(dataset.attrib, "apriori_complete", 0) == 1 || error(
            "a-priori file is not marked complete: $path")
        size(dataset["xa"]) == (34, 4) || error(
            "a-priori state array has the wrong shape")
        size(dataset["Sa_active"]) == (30, 30, 4) || error(
            "active a-priori covariance array has the wrong shape")
        active_parameter_indices =
            Int.(dataset["active_parameter_index"][:])
        active_parameter_indices == vcat(1, collect(6:34)) || error(
            "active-parameter layout is not the expected 30-state layout")
        aerosol_ln_aod_sigma = Float64(get(
            dataset.attrib, "aerosol_ln_aod_sigma", NaN))
        aerosol_ln_aod_sigma == EXPECTED_AEROSOL_LN_AOD_SIGMA || error(
            "aerosol_ln_aod_sigma is $aerosol_ln_aod_sigma, expected " *
            "$(EXPECTED_AEROSOL_LN_AOD_SIGMA) in the campaign-local prior")
        aod_active_indices = map(AOD_PARAMETER_INDICES) do parameter_index
            active_index = findfirst(==(parameter_index),
                                     active_parameter_indices)
            isnothing(active_index) && error(
                "AOD parameter $parameter_index is not active in the " *
                "campaign-local prior")
            active_index
        end
        for surface_index in 1:4
            covariance = Float64.(dataset["Sa_active"][:, :, surface_index])
            all(isfinite, covariance) || error(
                "surface $surface_index active covariance is non-finite")
            isposdef(Symmetric(covariance)) || error(
                "surface $surface_index active covariance is not positive definite")
            aod_variances = diag(covariance)[aod_active_indices]
            all(==(EXPECTED_AOD_VARIANCE), aod_variances) || error(
                "surface $surface_index active AOD covariance diagonal is " *
                "$(join(aod_variances, ", ")), expected " *
                "$(EXPECTED_AOD_VARIANCE) for all three ln(AOD760) states")
        end
        for band in ("o2a", "weak_co2", "strong_co2"), order in ("p1", "p2")
            key = "surface_$(order)_sigma_$(band)"
            Float64(dataset.attrib[key]) == 2e-3 || error(
                "$key is not 2e-3 in the campaign-local prior")
        end
        Float64(dataset.attrib["sif_wavelength_slope_prior_mw_m2_sr_nm2"]) == 0 ||
            error("SIF wavelength-slope prior is not centered on zero")
        Float64(dataset.attrib["sif_wavelength_slope_sigma_mw_m2_sr_nm2"]) ==
            0.002625 || error("SIF wavelength-slope prior sigma is not 0.002625")
    end
    return nothing
end

function validate_truth_selection(cases)
    length(cases) == 80 || error(
        "bottom-layer truth table must contain 80 states; found $(length(cases))")
    by_index = Dict(truth.state_index => truth for truth in cases)
    sort!(collect(keys(by_index))) == collect(1:80) || error(
        "bottom-layer truth indices are not exactly 1:80")
    ordered = TruthCase[]
    for state in ORDERED_AEROSOL_STATES
        haskey(by_index, state) || error("missing aerosol truth state $state")
        truth = by_index[state]
        truth.campaign == :bottom_layer_XCO2 || error(
            "state $state belongs to the wrong campaign")
        truth.co2_profile_mode == :bottom_layer || error(
            "state $state does not use the bottom-layer CO2 profile")
        truth.aerosol_case != :none || error(
            "state $state is not an aerosol scene")
        truth.fixed_upper_co2_ppm == truth.background_co2_ppm == 400.0 || error(
            "state $state does not use the fixed 400 ppm background")
        truth.bottom_layer_index == 16 || error(
            "state $state does not perturb layer 16")
        push!(ordered, truth)
    end
    all(truth -> truth.sif_case == :off, ordered[1:20]) || error(
        "the first 20 selected aerosol states are not all SIF-off")
    all(truth -> truth.sif_case != :off, ordered[21:40]) || error(
        "the final 20 selected aerosol states are not all SIF-on")
    return ordered
end

function validate_experiment_order(experiments)
    length(experiments) == 40 * 11 * 2 || error(
        "aerosol partition must contain 880 paired retrievals")
    state_order = unique(
        experiment.truth.state_index for experiment in experiments)
    state_order == ORDERED_AEROSOL_STATES || error(
        "aerosol truth-state order changed unexpectedly")
    for truth_state in ORDERED_AEROSOL_STATES
        state_experiments = filter(
            experiment -> experiment.truth.state_index == truth_state,
            experiments)
        length(state_experiments) == 22 || error(
            "state $truth_state does not contain 22 retrievals")
        paired_indices = [state_experiments[index].noise_index
                          for index in 1:2:length(state_experiments)]
        paired_indices == vcat(11, collect(1:10)) || error(
            "state $truth_state is not ordered perturbation 11 then 1:10")
        for index in 1:2:length(state_experiments)
            pair = state_experiments[index:index + 1]
            getfield.(pair, :measurement_class) == [:corrected, :uncorrected] ||
                error("state $truth_state has a non-adjacent measurement pair")
            pair[1].noise_index == pair[2].noise_index || error(
                "state $truth_state has mismatched pair perturbations")
            pair[1].random_seed == pair[2].random_seed || error(
                "state $truth_state has mismatched pair random seeds")
        end
    end
    return nothing
end

function main()
    campaign_root = abspath(get(
        ENV, "BOTTOM_RETRIEVAL_CAMPAIGN_ROOT", DEFAULT_CAMPAIGN_ROOT))
    truth_table = joinpath(campaign_root, "truth", "true_states.dat")
    measurement_directory = joinpath(campaign_root, "truth", "OCO_radiances")
    noise_directory = joinpath(measurement_directory, "noise_covariances")
    prior_path = joinpath(campaign_root, "retrieval_setup", "apriori_states.nc")
    output_root = joinpath(campaign_root, "retrievals")
    manifest_path = get(
        ENV, "BOTTOM_AEROSOL_MANIFEST_PATH",
        joinpath(output_root, "retrieval_manifest_wurst0_aerosol_all.dat"))

    isfile(truth_table) || error("missing campaign truth table: $truth_table")
    ordered_cases = validate_truth_selection(read_truth_cases(truth_table))
    validate_prior(prior_path)
    experiments = build_experiments(
        ordered_cases; measurement_directory, noise_directory)
    validate_experiment_order(experiments)
    write_manifest = lowercase(get(
        ENV, "BOTTOM_AEROSOL_WRITE_MANIFEST", "1")) in
        ("1", "true", "yes", "on")
    if write_manifest
        write_experiment_manifest(
            experiments; output_path=manifest_path, inversion_root=output_root)
    elseif !isfile(manifest_path)
        error("secondary worker cannot find authoritative aerosol manifest: " *
              manifest_path)
    end

    println("bottom-layer aerosol retrieval preflight passed")
    println("states_no_sif=$(join(lpad.(string.(NOSIF_AEROSOL_STATES), 3, '0'), ','))")
    println("states_sif=$(join(lpad.(string.(SIF_AEROSOL_STATES), 3, '0'), ','))")
    println("retrievals=$(length(experiments)) perturbation_order=11,1:10")
    println("manifest=$(abspath(manifest_path)) write_manifest=$write_manifest")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
