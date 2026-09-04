#!/usr/bin/env julia

"""
Validate paired bottom-layer CO2 retrieval products selected by the environment.

`VALIDATION_MODE=smoke` applies the smoke gate selected by
`BOTTOM_RETRIEVAL_SCENE_CLASS`:

- `clear_nosif` (default) requires the paired noiseless desert control
  (state 043).
- `aerosol_all_sif` accepts either or both paired noiseless 400-ppm aerosol
  controls: state 013 (SIF off) and state 018 (SIF on). This permits a
  no-SIF-first production queue to gate its two stages independently.
- `all` is reserved for structural product barriers spanning the clear and
  aerosol workers (and, after the no-SIF barrier, both SIF phases).

Scientific closure is mandatory for each corrected smoke retrieval; the
intentionally RRS-mismatched uncorrected retrieval must still be a complete,
finite, internally consistent OE result. `VALIDATION_MODE=products` performs
the latter structural validation for arbitrary selections in the selected
scene class. The default `auto` recognizes the canonical smoke-state set and
treats all other selections as product audits.
"""

using NCDatasets
using Printf
using LinearAlgebra

include(joinpath(@__DIR__, "RetrievalCases.jl"))
using .RetrievalCases
include(joinpath(@__DIR__, "RetrievalState.jl"))
using .RetrievalState

const DEFAULT_CAMPAIGN_ROOT = normpath(joinpath(
    @__DIR__, "..", "bottom_layer_XCO2_retrievals"))
const ACTIVE_LOG_AOD_INDICES = 14:16
const REQUIRED_LOG_AOD_VARIANCE = 0.75^2
const ACTIVE_LOG_AOD_NAMES = [
    "ln_sulfate_aod760",
    "ln_organic_carbon_aod760",
    "ln_utls_sulfate_aod760",
]

function parse_indices(value::AbstractString)
    indices = Int[]
    for token in split(value, ',')
        fields = split(strip(token), '-')
        if length(fields) == 1
            push!(indices, parse(Int, only(fields)))
        elseif length(fields) == 2
            append!(indices, parse(Int, fields[1]):parse(Int, fields[2]))
        else
            error("invalid index selection: $token")
        end
    end
    isempty(indices) && error("index selection is empty")
    length(unique(indices)) == length(indices) || error(
        "index selection contains duplicates")
    return indices
end

function finite_array(dataset, name)
    haskey(dataset, name) || error("missing variable $name")
    variable = dataset[name]
    # A single `[:]` is linear indexing in Julia and therefore flattens
    # matrices.  Preserve the stored NetCDF rank so covariance/Jacobian shape
    # and sub-block checks operate on their intended axes.
    indices = ntuple(_ -> Colon(), ndims(variable))
    values = Float64.(nomissing(variable[indices...], NaN))
    all(isfinite, values) || error("$name contains non-finite values")
    return values
end

"""Require a retrieval product to embed the selected campaign prior exactly."""
function validate_embedded_prior(parameter_names, prior_state,
                                 prior_covariance,
                                 expected_prior::RetrievalPrior, path)
    parameter_names == expected_prior.parameter_names || error(
        "embedded parameter names do not exactly match the " *
        "$(expected_prior.surface) campaign prior in $path")
    expected_prior.parameter_names[ACTIVE_LOG_AOD_INDICES] ==
        ACTIVE_LOG_AOD_NAMES || error(
            "campaign prior active indices 14:16 are not the three " *
            "ln(AOD760) parameters")
    prior_state == expected_prior.xa || error(
        "embedded a-priori state does not exactly match the " *
        "$(expected_prior.surface) campaign prior in $path")
    prior_covariance == expected_prior.Sa || error(
        "embedded active a-priori covariance does not exactly match the " *
        "$(expected_prior.surface) campaign prior in $path")

    for index in ACTIVE_LOG_AOD_INDICES
        expected_prior.Sa[index, index] == REQUIRED_LOG_AOD_VARIANCE || error(
            "campaign prior ln(AOD) variance at active index $index is " *
            "$(expected_prior.Sa[index, index]); expected " *
            "$REQUIRED_LOG_AOD_VARIANCE")
        prior_covariance[index, index] == REQUIRED_LOG_AOD_VARIANCE || error(
            "embedded ln(AOD) variance at active index $index is not " *
            "$REQUIRED_LOG_AOD_VARIANCE in $path")
    end
    return nothing
end

function validate_outcome_attributes(attributes, path;
                                     require_scientific_closure=false)
    Int(get(attributes, "retrieval_complete", 0)) == 1 || error(
        "retrieval is not marked complete: $path")
    outcome = Int(get(attributes, "outcome", -1))
    outcome in 1:4 || error("retrieval has illegal OE outcome $outcome: $path")
    converged = Int(get(attributes, "converged", -1))
    converged in (0, 1) || error("retrieval has invalid converged flag: $path")
    converged == Int(outcome in (1, 2)) || error(
        "retrieval outcome/converged metadata are inconsistent: $path")
    fit_quality_ok = Int(get(attributes, "fit_quality_ok", -1))
    fit_quality_ok in (0, 1) || error(
        "retrieval has invalid fit_quality_ok flag: $path")
    if require_scientific_closure
        outcome == 1 || error(
            "corrected desert-control closure requires OE outcome 1; " *
            "received $outcome: $path")
        converged == 1 || error(
            "corrected desert-control retrieval did not converge: $path")
        fit_quality_ok == 1 || error(
            "corrected desert-control retrieval failed fit quality: $path")
    end
    return (; outcome, converged=Bool(converged),
            fit_quality_ok=Bool(fit_quality_ok))
end

function validate_desert_control_selection(states, perturbations,
                                            truth_by_index)
    states == [43] || error(
        "smoke validation requires only desert control state 043; " *
        "received $(join(states, ','))")
    perturbations == [UNPERTURBED_INDEX] || error(
        "smoke validation requires only noiseless perturbation 11")
    haskey(truth_by_index, 43) || error("bottom truth table lacks state 043")
    truth = truth_by_index[43]
    truth.surface == :desert || error("state 043 is not the desert surface")
    truth.aerosol_case == :none || error("state 043 is not clear-sky")
    truth.sif_case == :off || error("state 043 is not SIF-off")
    truth.bottom_layer_index == 16 || error(
        "state 043 does not use bottom layer 16")
    truth.bottom_co2_ppm == truth.fixed_upper_co2_ppm == 400.0 || error(
        "state 043 is not the 400 ppm bottom-layer control")
    return nothing
end

function validate_aerosol_controls_selection(states, perturbations,
                                              truth_by_index)
    !isempty(states) && all(in((13, 18)), states) || error(
        "aerosol smoke validation accepts only states 013 and 018; " *
        "received $(join(states, ','))")
    length(unique(states)) == length(states) || error(
        "aerosol smoke validation received duplicate states")
    perturbations == [UNPERTURBED_INDEX] || error(
        "aerosol smoke validation requires only noiseless perturbation 11")
    for state in states
        haskey(truth_by_index, state) || error(
            "bottom truth table lacks aerosol smoke state $state")
        truth = truth_by_index[state]
        truth.surface == :urban || error(
            "aerosol smoke state $state is not the urban surface")
        truth.aerosol_case != :none || error(
            "aerosol smoke state $state does not contain aerosols")
        truth.bottom_layer_index == 16 || error(
            "aerosol smoke state $state does not use bottom layer 16")
        truth.bottom_co2_ppm == truth.fixed_upper_co2_ppm == 400.0 || error(
            "aerosol smoke state $state is not the 400 ppm control")
    end
    13 in states && truth_by_index[13].sif_case != :off && error(
        "aerosol smoke state 013 is not SIF-off")
    18 in states && truth_by_index[18].sif_case == :off && error(
        "aerosol smoke state 018 is not SIF-on")
    return nothing
end

function validate_scene_membership(truth, state, scene_class)
    if scene_class == "clear_nosif"
        truth.sif_case == :off || error("state $state is not SIF-off")
        truth.aerosol_case == :none || error("state $state is not clear-sky")
    elseif scene_class == "aerosol_all_sif"
        truth.aerosol_case != :none || error(
            "state $state is not an aerosol scene")
    elseif scene_class != "all"
        error("unsupported bottom-layer scene class: $scene_class")
    end
    return nothing
end

function validate_one(path, truth, perturbation, measurement_class, prior_path;
                      require_scientific_closure=false)
    isfile(path) || error("missing retrieval output: $path")
    expected_prior = load_retrieval_prior(truth.surface; path=prior_path)
    NCDataset(path) do dataset
        outcome = validate_outcome_attributes(
            dataset.attrib, path; require_scientific_closure)
        Int(dataset.attrib["truth_state_index"]) == truth.state_index || error(
            "truth-state mismatch in $path")
        Int(dataset.attrib["perturbation_index"]) == perturbation || error(
            "perturbation mismatch in $path")
        String(dataset.attrib["measurement_class"]) == String(measurement_class) || error(
            "measurement-class mismatch in $path")
        String(dataset.attrib["campaign"]) == "bottom_layer_XCO2" || error(
            "campaign mismatch in $path")
        String(dataset.attrib["co2_profile_mode"]) == "bottom_layer" || error(
            "CO2-profile mode mismatch in $path")
        Float64(dataset.attrib["fixed_upper_co2_ppm"]) == 400.0 || error(
            "fixed upper CO2 is not 400 ppm in $path")
        Int(dataset.attrib["truth_bottom_co2_layer_index"]) == 16 || error(
            "bottom-layer index mismatch in $path")
        Float64(dataset.attrib["truth_bottom_co2_ppm"]) == truth.bottom_co2_ppm || error(
            "bottom-layer VMR mismatch in $path")
        isapprox(Float64(dataset.attrib["truth_xco2_ppm"]), truth.xco2_ppm;
                 atol=1e-12, rtol=0) || error("truth XCO2 mismatch in $path")
        abspath(String(dataset.attrib["source_apriori"])) == prior_path || error(
            "a-priori provenance mismatch in $path")
        Int(dataset.attrib["nstreams"]) == 9 || error(
            "retrieval did not use nine streams in $path")
        String(dataset.attrib["jacobian_flavor"]) == "OCO_RRS_synth" || error(
            "Jacobian flavor mismatch in $path")
        String(dataset.attrib["retrieval_forward_scattering"]) == "noRS" || error(
            "forward-scattering mode mismatch in $path")

        final_chi_squared = finite_array(
            dataset, "final_band_reduced_chi_squared")
        maximum_chi_squared = Float64(dataset.attrib["maximum_band_chi_squared"])
        computed_fit_quality = all(<(maximum_chi_squared), final_chi_squared)
        computed_fit_quality == outcome.fit_quality_ok || error(
            "fit_quality_ok disagrees with terminal per-band reduced " *
            "chi-squared in $path")

        measurement = finite_array(dataset, "measurement_perturbed")
        length(measurement) == 2742 || error(
            "measurement length is $(length(measurement)); expected 2742")
        all(>(0), finite_array(dataset, "Se_diagonal")) || error(
            "nonpositive measurement variance in $path")
        size(dataset["final_jacobian"]) == (2742, 30) || error(
            "terminal Jacobian shape mismatch in $path")
        size(dataset["gain_matrix"]) == (30, 2742) || error(
            "gain-matrix shape mismatch in $path")
        size(dataset["a_priori_covariance"]) == (30, 30) || error(
            "a-priori covariance shape mismatch in $path")
        haskey(dataset.attrib, "parameter_names") || error(
            "missing embedded parameter names in $path")
        parameter_names = split(String(dataset.attrib["parameter_names"]))
        prior_state = finite_array(dataset, "a_priori_state")
        prior_covariance = finite_array(dataset, "a_priori_covariance")
        validate_embedded_prior(parameter_names, prior_state,
                                prior_covariance, expected_prior, path)
        for index in (21, 22, 24, 25, 27, 28)
            isapprox(prior_covariance[index, index], (2e-3)^2;
                     atol=0, rtol=1e-14) || error(
                "surface P1/P2 prior sigma is not 2e-3 at active index $index in $path")
        end
        lambda = 760.0
        conversion = 1e7
        transform = [lambda^2 / conversion 0.0;
                     -2lambda^3 / conversion^2 -lambda^4 / conversion^2]
        expected_sif_state = transform * [0.1, 0.0]
        expected_sif_covariance = transform *
            Diagonal([0.25^2, 0.002625^2]) * transform'
        prior_state[29:30] ≈ expected_sif_state || error(
            "zero-centered SIF slope prior mean mismatch in $path")
        prior_covariance[29:30, 29:30] ≈ expected_sif_covariance || error(
            "three-times-wider SIF slope prior covariance mismatch in $path")

        noise_expected = perturbation != UNPERTURBED_INDEX
        Bool(Int(dataset.attrib["noise_injected"])) == noise_expected || error(
            "noise-injection metadata mismatch in $path")
        if !noise_expected
            all(iszero, finite_array(dataset, "normalized_noise_draw")) || error(
                "noiseless retrieval has a nonzero normalized draw in $path")
            all(iszero, finite_array(dataset, "injected_measurement_noise")) || error(
                "noiseless retrieval has injected noise in $path")
        end
    end
    return nothing
end

function main()
    campaign_root = abspath(get(
        ENV, "BOTTOM_RETRIEVAL_CAMPAIGN_ROOT", DEFAULT_CAMPAIGN_ROOT))
    output_root = abspath(get(
        ENV, "RETRIEVAL_OUTPUT_ROOT", joinpath(campaign_root, "retrievals")))
    prior_path = abspath(get(
        ENV, "RETRIEVAL_PRIOR_PATH",
        joinpath(campaign_root, "retrieval_setup", "apriori_states.nc")))
    truth_table = joinpath(campaign_root, "truth", "true_states.dat")
    truth_by_index = Dict(
        truth.state_index => truth for truth in read_truth_cases(truth_table))
    states = parse_indices(get(ENV, "EXPECTED_STATES", "43"))
    perturbations = parse_indices(get(ENV, "EXPECTED_PERTURBATIONS", "11"))
    validation_mode = lowercase(get(ENV, "VALIDATION_MODE", "auto"))
    validation_mode in ("auto", "smoke", "products") || error(
        "VALIDATION_MODE must be auto, smoke, or products")
    scene_class = lowercase(get(
        ENV, "BOTTOM_RETRIEVAL_SCENE_CLASS", "clear_nosif"))
    scene_class in ("clear_nosif", "aerosol_all_sif", "all") || error(
        "BOTTOM_RETRIEVAL_SCENE_CLASS must be clear_nosif, " *
        "aerosol_all_sif, or all")
    canonical_smoke_states = scene_class == "clear_nosif" ? [43] : [13, 18]
    smoke_mode = validation_mode == "smoke" ||
        (validation_mode == "auto" && states == canonical_smoke_states &&
         perturbations == [UNPERTURBED_INDEX])
    if smoke_mode
        scene_class == "clear_nosif" ?
            validate_desert_control_selection(
                states, perturbations, truth_by_index) :
            validate_aerosol_controls_selection(
                states, perturbations, truth_by_index)
    end

    count = 0
    for state in states, perturbation in perturbations,
            measurement_class in (:corrected, :uncorrected)
        haskey(truth_by_index, state) || error("unknown bottom truth state $state")
        truth = truth_by_index[state]
        validate_scene_membership(truth, state, scene_class)
        path = joinpath(
            output_root, String(measurement_class),
            @sprintf("retrieval_state%03d_perturbation%02d.nc",
                     state, perturbation))
        require_scientific_closure = smoke_mode &&
            measurement_class == :corrected
        validate_one(
            path, truth, perturbation, measurement_class, prior_path;
            require_scientific_closure)
        count += 1
    end
    println("validated bottom-layer retrievals: files=$count " *
            "states=$(join(states, ',')) perturbations=$(join(perturbations, ','))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
