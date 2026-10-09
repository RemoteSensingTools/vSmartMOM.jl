#!/usr/bin/env julia

"""
Compare the isolated sigma(ln(AOD760)) = 0.75 noiseless retrievals with the
unchanged sigma = 2 production baselines for states 001 and 013.

The script is intentionally all-or-nothing. It validates all four sensitivity
products and all four baselines before creating either report. Retrieval files
are opened read-only and are never changed.

Usage:

    julia --project=. \
        RRS_XCO2/inversion/compare_ln_aod_sigma_0p75_sensitivity.jl

Outputs:

    RRS_XCO2/bottom_layer_XCO2_retrievals/sensitivity_tests/
        ln_aod_sigma_0p75/sigma2_vs_sigma0p75.dat
        ln_aod_sigma_0p75/sigma2_vs_sigma0p75.md
"""

using Dates
using NCDatasets
using Printf
using SHA

const HERE = abspath(@__DIR__)
const CAMPAIGN_ROOT = normpath(joinpath(
    HERE, "..", "bottom_layer_XCO2_retrievals"))
const TEST_ROOT = joinpath(
    CAMPAIGN_ROOT, "sensitivity_tests", "ln_aod_sigma_0p75")
const BASELINE_ROOT = joinpath(CAMPAIGN_ROOT, "retrievals")
const TEST_OUTPUT_ROOT = joinpath(TEST_ROOT, "retrievals")
const BASELINE_PRIOR = joinpath(
    CAMPAIGN_ROOT, "retrieval_setup", "apriori_states.nc")
const TEST_PRIOR = joinpath(TEST_ROOT, "retrieval_setup", "apriori_states.nc")
const SCENE_COMPONENTS = joinpath(CAMPAIGN_ROOT, "truth", "scene_components.dat")
const DAT_REPORT = joinpath(TEST_ROOT, "sigma2_vs_sigma0p75.dat")
const MARKDOWN_REPORT = joinpath(TEST_ROOT, "sigma2_vs_sigma0p75.md")

const STATES = (1, 13)
const CLASSES = (:corrected, :uncorrected)
const PERTURBATION = 11
const BAND_NAMES = ("o2a", "weak_co2", "strong_co2")
const AOD_NAMES = (
    "ln_sulfate_aod760",
    "ln_organic_carbon_aod760",
    "ln_utls_sulfate_aod760",
)
const DAT_ROW_FORMAT = Printf.Format(
    "%03d %s %.12g %s %d %s %d %d %d %d %d %d %.12g %.12g " *
    "%.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g %.12g " *
    "%.12g %.12g %.12g %.12g %.12g %.12g %s\n")

retrieval_path(root, state, class) = joinpath(
    root, String(class),
    @sprintf("retrieval_state%03d_perturbation%02d.nc", state, PERTURBATION))

function file_sha256(path)
    open(path, "r") do stream
        return bytes2hex(sha256(stream))
    end
end

outcome_name(value) = value == 1 ? "converged_fit_pass" :
    value == 2 ? "converged_fit_fail" :
    value == 3 ? "maximum_iterations" :
    value == 4 ? "maximum_divergences" : "unknown_$(value)"

yesno(value) = value ? "yes" : "no"
fmt(value) = @sprintf("%.7g", value)
fmt_signed(value) = @sprintf("%+.7g", value)
fmt_chi(values) = join(fmt.(values), " / ")

function require_inputs()
    missing = String[]
    for path in (BASELINE_PRIOR, TEST_PRIOR, SCENE_COMPONENTS)
        isfile(path) || push!(missing, path)
    end
    for state in STATES, class in CLASSES
        for root in (BASELINE_ROOT, TEST_OUTPUT_ROOT)
            path = retrieval_path(root, state, class)
            isfile(path) || push!(missing, path)
        end
    end
    isempty(missing) || error(
        "comparison requires every completed baseline and sensitivity product; " *
        "missing:\n  " * join(missing, "\n  "))
end

function read_truth_aod760(path, aerosol_case)
    in_aerosol_section = false
    for raw_line in eachline(path)
        line = strip(raw_line)
        isempty(line) && continue
        if line == "[AEROSOL_CASES]"
            in_aerosol_section = true
            continue
        elseif startswith(line, "[")
            in_aerosol_section = false
            continue
        end
        (in_aerosol_section && !startswith(line, '#')) || continue
        fields = split(line)
        first(fields) == aerosol_case || continue
        length(fields) >= 9 || error("malformed aerosol row in $path: $line")
        species = parse.(Float64, fields[6:8])
        total = parse(Float64, fields[9])
        isapprox(sum(species), total; atol=2e-9, rtol=0) || error(
            "species AOD760 values do not sum to the catalog total for $aerosol_case")
        return (; species, total)
    end
    error("aerosol case $aerosol_case was not found in $path")
end

function validate_prior_sigma(covariance, names, expected_sigma, path)
    indices = map(AOD_NAMES) do name
        index = findfirst(==(name), names)
        isnothing(index) && error("missing $name in $path")
        index
    end
    variances = [Float64(covariance[index, index]) for index in indices]
    expected_variance = expected_sigma^2
    all(==(expected_variance), variances) || error(
        "unexpected log-AOD prior variances in $path: $(variances); " *
        "expected $expected_variance")
    return collect(indices)
end

function read_prior_slice(path, surface)
    NCDataset(path, "r") do dataset
        Int(get(dataset.attrib, "apriori_complete", 0)) == 1 || error(
            "prior file is not marked complete: $path")
        surfaces = split(String(dataset.attrib["surface_order"]))
        surface_index = findfirst(==(surface), surfaces)
        isnothing(surface_index) && error("surface $surface is absent from $path")
        active_to_full = Int.(dataset["active_parameter_index"][:])
        full_names = split(String(dataset.attrib["parameter_names"]))
        maximum(active_to_full) <= length(full_names) || error(
            "active-state index exceeds the parameter-name list in $path")
        xa_full = Float64.(dataset["xa"][:, surface_index])
        covariance = Float64.(dataset["Sa_active"][:, :, surface_index])
        return (;
            names=full_names[active_to_full],
            state=xa_full[active_to_full],
            covariance,
        )
    end
end

function read_result(path, expected_state, expected_class, expected_sigma,
                     expected_prior)
    NCDataset(path, "r") do dataset
        attributes = dataset.attrib
        Int(get(attributes, "retrieval_complete", 0)) == 1 || error(
            "sensitivity comparison refuses incomplete retrieval: $path")
        Int(attributes["truth_state_index"]) == expected_state || error(
            "state mismatch in $path")
        Int(attributes["perturbation_index"]) == PERTURBATION || error(
            "perturbation mismatch in $path")
        String(attributes["measurement_class"]) == String(expected_class) || error(
            "measurement-class mismatch in $path")
        Int(attributes["noise_injected"]) == 0 || error(
            "comparison accepts only noiseless perturbation 11: $path")
        all(iszero, dataset["normalized_noise_draw"][:]) || error(
            "stored normalized noise is nonzero in $path")
        normpath(String(attributes["source_apriori"])) == normpath(expected_prior) ||
            error("unexpected source_apriori in $path")

        surface = String(attributes["surface"])
        names = split(String(attributes["parameter_names"]))
        prior_state = Float64.(dataset["a_priori_state"][:])
        prior_covariance = Float64.(dataset["a_priori_covariance"][:, :])
        prior_reference = read_prior_slice(expected_prior, surface)
        names == prior_reference.names || error(
            "embedded parameter names differ from $expected_prior in $path")
        isequal(prior_state, prior_reference.state) || error(
            "embedded a-priori state differs from $expected_prior in $path")
        isequal(prior_covariance, prior_reference.covariance) || error(
            "embedded a-priori covariance differs from $expected_prior in $path")
        aod_indices = validate_prior_sigma(
            prior_covariance, names, expected_sigma, path)

        final_state = Float64.(dataset["final_state"][:])
        state_history = Float64.(dataset["state_at_trial"][:, :])
        size(state_history, 1) == length(names) || error(
            "state-history dimension mismatch in $path")
        final_aod = exp.(final_state[aod_indices])
        evaluated_aod = exp.(state_history[aod_indices, :])
        all(isfinite, final_aod) && all(isfinite, evaluated_aod) || error(
            "non-finite physical AOD in $path")
        # A converged solve's separately evaluated terminal state is not a
        # trial in state_at_trial; include it in the evaluated-state maxima.
        all_evaluated_aod = hcat(evaluated_aod, final_aod)

        trial_index = Int.(dataset["trial_index"][:])
        accepted = Bool.(dataset["trial_accepted"][:])
        divergent = Bool.(dataset["trial_divergent"][:])
        d_sigma = Float64.(dataset["d_sigma_sq_scaled"][:])
        evaluation_seconds = Float64.(dataset["evaluation_seconds"][:])
        ntrial = length(trial_index)
        all(length(values) == ntrial for values in
            (accepted, divergent, d_sigma, evaluation_seconds)) || error(
            "trial-vector length mismatch in $path")
        trial_index == collect(1:ntrial) || error(
            "trial indices are not contiguous in $path")
        terminal_record = findlast(accepted)
        isnothing(terminal_record) && error("retrieval has no accepted trial: $path")
        divergence_count = Int(attributes["divergence_count"])
        divergence_count == count(divergent) || error(
            "divergence count disagrees with trial flags in $path")

        outcome = Int(attributes["outcome"])
        converged = Bool(attributes["converged"])
        fit_quality_ok = Bool(attributes["fit_quality_ok"])
        converged == (outcome in (1, 2)) || error(
            "outcome/converged mismatch in $path")
        convergence_threshold = Float64(attributes["convergence_threshold"])
        terminal_d_sigma_sq_scaled = d_sigma[terminal_record]
        converged == (terminal_d_sigma_sq_scaled < convergence_threshold) || error(
            "terminal accepted convergence metric disagrees with outcome in $path")

        chi_squared = Float64.(dataset["final_band_reduced_chi_squared"][:])
        length(chi_squared) == length(BAND_NAMES) || error(
            "expected three final band chi-squares in $path")
        maximum_band_chi_squared = Float64(attributes["maximum_band_chi_squared"])
        fit_quality_ok == all(<(maximum_band_chi_squared), chi_squared) || error(
            "fit-quality flag disagrees with final chi-squares in $path")

        final_evaluation_seconds = Float64(attributes["final_evaluation_seconds"])
        trial_evaluation_seconds = sum(evaluation_seconds)
        # A converged solve performs an additional evaluation at its proposed
        # terminal state. Outcomes 3/4 reuse an evaluated accepted trial.
        total_evaluation_seconds = trial_evaluation_seconds +
            (converged ? final_evaluation_seconds : 0.0)

        return (;
            path=abspath(path), state=expected_state, class=expected_class,
            sigma_ln_aod=expected_sigma,
            surface,
            aerosol_case=String(attributes["aerosol_case"]),
            outcome, outcome_name=outcome_name(outcome), converged,
            fit_quality_ok, trial_count=ntrial,
            accepted_count=count(accepted), divergence_count,
            terminal_accepted_trial=trial_index[terminal_record],
            terminal_d_sigma_sq_scaled, convergence_threshold, chi_squared,
            truth_xco2=Float64(attributes["truth_xco2_ppm"]),
            xco2=Float64(dataset["XCO2"][]), final_aod,
            final_total_aod=sum(final_aod),
            max_evaluated_sulfate_aod=maximum(all_evaluated_aod[1, :]),
            max_evaluated_total_aod=maximum(vec(sum(all_evaluated_aod; dims=1))),
            trial_evaluation_seconds, final_evaluation_seconds,
            total_evaluation_seconds,
        )
    end
end

function validate_pair(baseline, test)
    baseline.state == test.state || error("paired state mismatch")
    baseline.class == test.class || error("paired class mismatch")
    baseline.surface == test.surface || error(
        "paired surface mismatch for state $(baseline.state)")
    baseline.aerosol_case == test.aerosol_case || error(
        "paired aerosol-case mismatch for state $(baseline.state)")
    baseline.truth_xco2 == test.truth_xco2 || error(
        "paired truth-XCO2 mismatch for state $(baseline.state)")
end

function atomic_write(writer, path)
    mkpath(dirname(path))
    temporary, stream = mktemp(dirname(path))
    try
        writer(stream)
        close(stream)
        mv(temporary, path; force=true)
    catch
        isopen(stream) && close(stream)
        isfile(temporary) && rm(temporary; force=true)
        rethrow()
    end
end

function write_dat(path, results, baseline_hash, test_hash)
    atomic_write(path) do stream
        println(stream, "# sigma(ln(AOD760)) sensitivity comparison; perturbation 11 is exactly noiseless.")
        println(stream, "# baseline_prior_sha256 $baseline_hash")
        println(stream, "# test_prior_sha256 $test_hash")
        println(stream, "# total_evaluation_seconds includes the separate terminal evaluation only for converged outcomes.")
        println(stream,
            "state class sigma_ln_aod configuration outcome outcome_name converged " *
            "fit_quality_ok trial_count accepted_count divergence_count " *
            "terminal_accepted_trial terminal_accepted_d_sigma_sq_scaled " *
            "convergence_threshold chi2_o2a chi2_weak_co2 chi2_strong_co2 " *
            "truth_xco2_ppm xco2_ppm xco2_error_ppm " *
            "final_sulfate_aod760 final_organic_carbon_aod760 " *
            "final_utls_sulfate_aod760 final_total_aod760 " *
            "max_evaluated_sulfate_aod760 max_evaluated_total_aod760 " *
            "trial_evaluation_seconds final_evaluation_seconds " *
            "total_evaluation_seconds source_file")
        for result in results
            configuration = result.sigma_ln_aod == 2.0 ? "baseline" : "test"
            Printf.format(stream, DAT_ROW_FORMAT,
                result.state, String(result.class), result.sigma_ln_aod,
                configuration, result.outcome, result.outcome_name,
                Int(result.converged), Int(result.fit_quality_ok),
                result.trial_count, result.accepted_count,
                result.divergence_count, result.terminal_accepted_trial,
                result.terminal_d_sigma_sq_scaled, result.convergence_threshold,
                result.chi_squared..., result.truth_xco2, result.xco2,
                result.xco2 - result.truth_xco2, result.final_aod...,
                result.final_total_aod, result.max_evaluated_sulfate_aod,
                result.max_evaluated_total_aod,
                result.trial_evaluation_seconds,
                result.final_evaluation_seconds,
                result.total_evaluation_seconds, result.path)
        end
    end
end

function arrow(baseline, test; signed_delta=true)
    delta = test - baseline
    suffix = signed_delta ? " (Delta $(fmt_signed(delta)))" : ""
    return "$(fmt(baseline)) -> $(fmt(test))$suffix"
end

function write_markdown(path, pairs, truth_aod, baseline_hash, test_hash)
    atomic_write(path) do stream
        println(stream, "# Log-AOD prior sensitivity: sigma 2 vs 0.75")
        println(stream)
        println(stream, "Generated: `$(Dates.format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"))`")
        println(stream)
        println(stream,
            "All rows use perturbation 11 (exactly noiseless). Arrows show " *
            "`sigma=2 -> sigma=0.75`; Delta is test minus baseline. The state-step " *
            "convergence threshold is 2 and the per-band fit threshold is 1.4.")
        println(stream)
        println(stream, "- Baseline prior SHA-256: `$baseline_hash`")
        println(stream, "- Test prior SHA-256: `$test_hash`")
        println(stream)
        println(stream, "## Convergence and timing")
        println(stream)
        println(stream,
            "| State | Class | Outcome | Converged | Fit pass | Terminal accepted d_sigma_sq_scaled | Trials | Accepted | Divergences | Total evaluation seconds |")
        println(stream,
            "|---:|:---|:---|:---:|:---:|---:|---:|---:|---:|---:|")
        for (baseline, test) in pairs
            @printf(stream,
                "| %03d | %s | %d (%s) -> %d (%s) | %s -> %s | %s -> %s | %s | %d -> %d | %d -> %d | %d -> %d | %s |\n",
                baseline.state, String(baseline.class),
                baseline.outcome, baseline.outcome_name,
                test.outcome, test.outcome_name,
                yesno(baseline.converged), yesno(test.converged),
                yesno(baseline.fit_quality_ok), yesno(test.fit_quality_ok),
                arrow(baseline.terminal_d_sigma_sq_scaled,
                      test.terminal_d_sigma_sq_scaled),
                baseline.trial_count, test.trial_count,
                baseline.accepted_count, test.accepted_count,
                baseline.divergence_count, test.divergence_count,
                arrow(baseline.total_evaluation_seconds,
                      test.total_evaluation_seconds))
        end
        println(stream)
        println(stream,
            "Timing sums all evaluated trials and, for converged outcomes only, " *
            "the additional terminal-state evaluation. Outcomes 3/4 reuse an " *
            "already evaluated accepted trial and are not double-counted.")
        println(stream)
        println(stream, "## Final spectral fit and XCO2")
        println(stream)
        println(stream,
            "| State | Class | Reduced chi-squared (O2A / weak / strong) | XCO2 ppm | XCO2 error ppm |")
        println(stream, "|---:|:---|:---|---:|---:|")
        for (baseline, test) in pairs
            @printf(stream,
                "| %03d | %s | %s -> %s | %s | %s |\n",
                baseline.state, String(baseline.class),
                fmt_chi(baseline.chi_squared), fmt_chi(test.chi_squared),
                arrow(baseline.xco2, test.xco2),
                arrow(baseline.xco2 - baseline.truth_xco2,
                      test.xco2 - test.truth_xco2))
        end
        println(stream)
        println(stream, "## Aerosol AOD760")
        println(stream)
        println(stream,
            "| State | Class | Final sulfate | Final organic carbon | Final UTLS sulfate | Final total | Maximum evaluated sulfate | Maximum evaluated total | Truth total |")
        println(stream, "|---:|:---|---:|---:|---:|---:|---:|---:|---:|")
        for (baseline, test) in pairs
            truth = truth_aod[baseline.state]
            @printf(stream,
                "| %03d | %s | %s | %s | %s | %s | %s | %s | %s |\n",
                baseline.state, String(baseline.class),
                arrow(baseline.final_aod[1], test.final_aod[1]),
                arrow(baseline.final_aod[2], test.final_aod[2]),
                arrow(baseline.final_aod[3], test.final_aod[3]),
                arrow(baseline.final_total_aod, test.final_total_aod),
                arrow(baseline.max_evaluated_sulfate_aod,
                      test.max_evaluated_sulfate_aod),
                arrow(baseline.max_evaluated_total_aod,
                      test.max_evaluated_total_aod), fmt(truth.total))
        end
        println(stream)
        println(stream,
            "Maximum evaluated AOD includes rejected LM trial states; final AOD " *
            "is the terminal accepted state (or the separately evaluated " *
            "terminal proposal for a converged solve).")
    end
end

function main()
    isempty(ARGS) || error("this fixed campaign comparison accepts no arguments")
    # Do not create report paths until every expected retrieval exists.
    require_inputs()

    baseline_hash = file_sha256(BASELINE_PRIOR)
    test_hash = file_sha256(TEST_PRIOR)
    results = NamedTuple[]
    pairs = Tuple{NamedTuple,NamedTuple}[]
    truth_aod = Dict{Int,NamedTuple}()

    for state in STATES, class in CLASSES
        baseline = read_result(
            retrieval_path(BASELINE_ROOT, state, class),
            state, class, 2.0, BASELINE_PRIOR)
        test = read_result(
            retrieval_path(TEST_OUTPUT_ROOT, state, class),
            state, class, 0.75, TEST_PRIOR)
        validate_pair(baseline, test)
        push!(results, baseline, test)
        push!(pairs, (baseline, test))
        truth_aod[state] = read_truth_aod760(
            SCENE_COMPONENTS, baseline.aerosol_case)
    end

    write_dat(DAT_REPORT, results, baseline_hash, test_hash)
    write_markdown(MARKDOWN_REPORT, pairs, truth_aod,
                   baseline_hash, test_hash)
    println("wrote $(abspath(DAT_REPORT))")
    println("wrote $(abspath(MARKDOWN_REPORT))")
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
