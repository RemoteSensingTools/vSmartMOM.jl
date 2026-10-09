#!/usr/bin/env julia

"""
Compare the completed state-001 (no aerosol) and state-009 (aerosol) retrieval
ensembles.

The comparison is deliberately strict: both corrected and uncorrected files
must be complete for every perturbation 01:11 in both truth states.  Noisy
perturbations 01:10 define the empirical mean and sample variance;
perturbation 11 is the exact noiseless retrieval and is reported separately.

Usage:

    julia --project=. RRS_XCO2/inversion/compare_state001_state009_ensembles.jl

    # Inventory only; useful while state 009 is still running.
    julia --project=. RRS_XCO2/inversion/compare_state001_state009_ensembles.jl \
        --status

Optional paths:

    --inversion-root PATH
    --output-prefix PATH

The output prefix defaults to
`RRS_XCO2/inversion/state001_vs_state009_retrieval_comparison` and produces:

* `<prefix>_runs.dat`: one row per retrieval, including convergence and timing;
* `<prefix>_state_statistics.dat`: physical-state ensemble and posterior
  statistics; and
* `<prefix>_summary.md`: a human-readable comparison.

Layer CO2 is collapsed to the saved XCO2 diagnostic.  Its posterior variance
uses the tangent of the same pressure-dependent dry-column mapping independently
implemented and regression-checked by `plot_corrected_vs_uncorrected_errors.py`.
Log-AOD and log-height posterior covariances are mapped with the delta method.
SIF760 is displayed per wavelength at 760 nm, as in the existing Python tools;
mSIF remains in its native retrieval coordinate.
"""

using Dates
using LinearAlgebra: dot
using NCDatasets
using Printf
using Statistics

const HERE = @__DIR__
const TRUTH_TABLE = normpath(joinpath(HERE, "..", "truth_map", "true_states.dat"))
const SCENE_COMPONENTS = normpath(joinpath(
    HERE, "..", "truth_map", "scene_components.dat"))
const STATES = (1, 9)
const CLASSES = (:corrected, :uncorrected)
const PERTURBED = 1:10
const UNPERTURBED = 11
const EXPECTED_PERTURBATIONS = 1:11
const BAND_NAMES = ("o2a", "weak_co2", "strong_co2")

# These constants duplicate the independent XCO2 regression mapping in
# plot_corrected_vs_uncorrected_errors.py.  They originate from the same
# reduced 16-layer, 1000-hPa profile as the retrieval forward model.  Only the
# bottom-layer dry column changes with retrieved surface pressure.
const DRY_COLUMN_FRACTIONS_1000 = [
    0.06271458366696699, 0.06267999123579230,
    0.06264827279449964, 0.06261831034665068,
    0.06258999065374835, 0.06255915670625688,
    0.06253127935109044, 0.06250621030757231,
    0.06248208275638748, 0.06245969764265345,
    0.06243817296080489, 0.06241161154141432,
    0.06238468284683572, 0.06235227383324737,
    0.06232590314728601, 0.06229778020879323,
]
const REFERENCE_SURFACE_PRESSURE_HPA = 1000.0
const BOTTOM_LAYER_TOP_PRESSURE_HPA = 937.50875
const SIF_WAVENUMBER_TO_WAVELENGTH_760 = 1.0e7 / 760.0^2

const PARAMETER_SPECS = [
    (key="XCO2", unit="ppm"),
    (key="psurf", unit="hPa"),
    (key="sulfate_aod760", unit="1"),
    (key="organic_carbon_aod760", unit="1"),
    (key="utls_sulfate_aod760", unit="1"),
    (key="sulfate_z0", unit="km"),
    (key="organic_carbon_z0", unit="km"),
    (key="utls_sulfate_z0", unit="km"),
    (key="o2a_surface_P0", unit="1"),
    (key="o2a_surface_P1", unit="1"),
    (key="o2a_surface_P2", unit="1"),
    (key="weak_co2_surface_P0", unit="1"),
    (key="weak_co2_surface_P1", unit="1"),
    (key="weak_co2_surface_P2", unit="1"),
    (key="strong_co2_surface_P0", unit="1"),
    (key="strong_co2_surface_P1", unit="1"),
    (key="strong_co2_surface_P2", unit="1"),
    (key="SIF760", unit="mW_m-2_sr-1_nm-1"),
    (key="mSIF", unit="native"),
]

function usage(io::IO=stdout)
    println(io, "Usage: julia --project=. $(relpath(@__FILE__, pwd())) [options]")
    println(io, "  --status                 Inventory inputs without writing reports")
    println(io, "  --inversion-root PATH    corrected/ and uncorrected/ parent")
    println(io, "  --output-prefix PATH     Output prefix without an extension")
    println(io, "  -h, --help               Show this help")
end

function parse_arguments(args)
    inversion_root = HERE
    output_prefix = nothing
    status = false
    index = 1
    while index <= length(args)
        argument = args[index]
        if argument == "--status"
            status = true
        elseif argument in ("-h", "--help")
            usage()
            return nothing
        elseif argument in ("--inversion-root", "--output-prefix")
            index == length(args) && error("$argument requires a path")
            value = args[index + 1]
            argument == "--inversion-root" ?
                (inversion_root = abspath(value)) : (output_prefix = abspath(value))
            index += 1
        elseif startswith(argument, "--inversion-root=")
            inversion_root = abspath(split(argument, '='; limit=2)[2])
        elseif startswith(argument, "--output-prefix=")
            output_prefix = abspath(split(argument, '='; limit=2)[2])
        else
            error("unknown argument: $argument")
        end
        index += 1
    end
    isnothing(output_prefix) && (output_prefix = joinpath(
        inversion_root, "state001_vs_state009_retrieval_comparison"))
    return (; inversion_root, output_prefix, status)
end

retrieval_path(root, state, class, perturbation) = joinpath(
    root, String(class), @sprintf(
        "retrieval_state%03d_perturbation%02d.nc", state, perturbation))

function completion_status(path, state, class, perturbation)
    isfile(path) || return (; ready=false, reason="missing")
    try
        return NCDataset(path) do dataset
            get(dataset.attrib, "retrieval_complete", 0) == 1 ||
                return (; ready=false, reason="not_marked_complete")
            Int(dataset.attrib["truth_state_index"]) == state ||
                return (; ready=false, reason="wrong_truth_state")
            Symbol(String(dataset.attrib["measurement_class"])) == class ||
                return (; ready=false, reason="wrong_measurement_class")
            Int(dataset.attrib["perturbation_index"]) == perturbation ||
                return (; ready=false, reason="wrong_perturbation")
            return (; ready=true, reason="complete")
        end
    catch exception
        return (; ready=false,
                reason="unreadable_$(nameof(typeof(exception)))")
    end
end

function inventory(root; io=stdout)
    all_ready = true
    for state in STATES, class in CLASSES
        complete = Int[]
        pending = String[]
        for perturbation in EXPECTED_PERTURBATIONS
            path = retrieval_path(root, state, class, perturbation)
            status = completion_status(path, state, class, perturbation)
            if status.ready
                push!(complete, perturbation)
            else
                all_ready = false
                push!(pending, @sprintf("%02d:%s", perturbation, status.reason))
            end
        end
        complete_text = isempty(complete) ? "none" :
            join((@sprintf("%02d", value) for value in complete), ",")
        pending_text = isempty(pending) ? "none" : join(pending, ",")
        @printf(io, "state=%03d class=%-11s complete=%s pending=%s\n",
                state, String(class), complete_text, pending_text)
    end
    println(io, "ready=", all_ready)
    return all_ready
end

function require_complete_inputs(root)
    inventory(root) && return nothing
    error("state 001/009 comparison requires matched complete corrected and " *
          "uncorrected perturbations 01:11; use --status while runs continue")
end

function table_row(path, state_index)
    isfile(path) || error("missing truth table: $path")
    names = String[]
    for line in eachline(path)
        stripped = strip(line)
        if startswith(stripped, "# index ")
            names = split(strip(stripped[2:end]))
        elseif !isempty(stripped) && !startswith(stripped, '#')
            isempty(names) && error("truth-table header was not found in $path")
            values = split(stripped)
            length(values) == length(names) || error(
                "truth-table row has $(length(values)) fields; expected $(length(names))")
            row = Dict(zip(names, values))
            parse(Int, row["index"]) == state_index && return row
        end
    end
    error("truth state $state_index was not found in $path")
end

function aerosol_aod760(path, aerosol_case)
    isfile(path) || error("missing scene-component table: $path")
    active = false
    names = String[]
    for line in eachline(path)
        stripped = strip(line)
        if stripped == "[AEROSOL_CASES]"
            active = true
            continue
        elseif active && startswith(stripped, '[')
            break
        elseif !active || isempty(stripped)
            continue
        elseif startswith(stripped, "# case")
            names = split(strip(stripped[2:end]))
        elseif !startswith(stripped, '#') && !isempty(names)
            values = split(stripped)
            length(values) == length(names) || error(
                "malformed aerosol-case row in $path")
            row = Dict(zip(names, values))
            row["case"] == aerosol_case || continue
            return Dict(
                "sulfate_aod760" => parse(Float64, row["sulfate_AOD760"]),
                "organic_carbon_aod760" =>
                    parse(Float64, row["organic_AOD760"]),
                "utls_sulfate_aod760" =>
                    parse(Float64, row["utls_sulfate_AOD760"]),
            )
        end
    end
    error("aerosol case '$aerosol_case' was not found in $path")
end

function parameter_indices(names)
    length(unique(names)) == length(names) || error("duplicate parameter names")
    return Dict(name => index for (index, name) in enumerate(names))
end

function physical_state(names, state, saved_xco2)
    indices = parameter_indices(names)
    native(name) = state[indices[name]]
    values = Dict{String,Float64}(
        "XCO2" => Float64(saved_xco2),
        "psurf" => native("psurf"),
    )
    for species in ("sulfate", "organic_carbon", "utls_sulfate")
        values["$(species)_aod760"] = exp(native("ln_$(species)_aod760"))
        values["$(species)_z0"] = exp(native("ln_$(species)_z0"))
    end
    for band in BAND_NAMES, order in 0:2
        key = "$(band)_surface_P$order"
        values[key] = native(key)
    end
    values["SIF760"] = native("SIF760") *
        SIF_WAVENUMBER_TO_WAVELENGTH_760
    values["mSIF"] = native("mSIF")
    return values
end

"""Return independently evaluated XCO2 and its native-state gradient."""
function xco2_value_and_gradient(names, state, fixed_upper_ppm)
    indices = parameter_indices(names)
    psurf = state[indices["psurf"]]
    bottom_thickness = psurf - BOTTOM_LAYER_TOP_PRESSURE_HPA
    reference_thickness = REFERENCE_SURFACE_PRESSURE_HPA -
        BOTTOM_LAYER_TOP_PRESSURE_HPA
    bottom_thickness > 0 || error(
        "surface pressure $psurf hPa lies above the bottom-layer top")

    weights = copy(DRY_COLUMN_FRACTIONS_1000)
    weights[end] *= bottom_thickness / reference_thickness
    total = sum(weights)
    co2 = fill(Float64(fixed_upper_ppm) * 1e-6, 16)
    for layer in 5:16
        co2[layer] = state[indices[@sprintf("co2_vmr_layer%02d", layer)]]
    end
    mean_vmr = sum(co2 .* weights) / total

    gradient = zeros(Float64, length(state))
    for layer in 5:16
        gradient[indices[@sprintf("co2_vmr_layer%02d", layer)]] =
            1.0e6 * weights[layer] / total
    end
    bottom_weight_derivative = DRY_COLUMN_FRACTIONS_1000[end] /
        reference_thickness
    gradient[indices["psurf"]] = 1.0e6 * bottom_weight_derivative *
        (co2[end] - mean_vmr) / total
    return 1.0e6 * mean_vmr, gradient
end

function physical_gradients(names, state, saved_xco2, fixed_upper_ppm)
    indices = parameter_indices(names)
    nstate = length(state)
    gradients = Dict{String,Vector{Float64}}()
    xco2, gradients["XCO2"] = xco2_value_and_gradient(
        names, state, fixed_upper_ppm)
    isapprox(xco2, saved_xco2; rtol=3e-7, atol=2e-5) || error(
        "independent XCO2 map does not reproduce saved XCO2 " *
        "($xco2 versus $saved_xco2 ppm)")

    function single_gradient(name, scale=1.0)
        gradient = zeros(Float64, nstate)
        gradient[indices[name]] = scale
        return gradient
    end
    gradients["psurf"] = single_gradient("psurf")
    for species in ("sulfate", "organic_carbon", "utls_sulfate")
        aod_name = "ln_$(species)_aod760"
        height_name = "ln_$(species)_z0"
        gradients["$(species)_aod760"] = single_gradient(
            aod_name, exp(state[indices[aod_name]]))
        gradients["$(species)_z0"] = single_gradient(
            height_name, exp(state[indices[height_name]]))
    end
    for band in BAND_NAMES, order in 0:2
        key = "$(band)_surface_P$order"
        gradients[key] = single_gradient(key)
    end
    gradients["SIF760"] = single_gradient(
        "SIF760", SIF_WAVENUMBER_TO_WAVELENGTH_760)
    gradients["mSIF"] = single_gradient("mSIF")
    return gradients
end

function physical_posterior_variances(names, state, saved_xco2,
                                      fixed_upper_ppm, posterior)
    size(posterior) == (length(state), length(state)) || error(
        "posterior covariance shape does not match the state")
    all(isfinite, posterior) || error("posterior covariance contains non-finite values")
    isapprox(posterior, posterior'; rtol=1e-8, atol=1e-15) || error(
        "posterior covariance is not symmetric")
    gradients = physical_gradients(
        names, state, saved_xco2, fixed_upper_ppm)
    variances = Dict{String,Float64}()
    for spec in PARAMETER_SPECS
        gradient = gradients[spec.key]
        variance = dot(gradient, posterior * gradient)
        tolerance = 100eps(Float64) * max(1.0, maximum(abs, posterior)) *
            max(1.0, sum(abs, gradient)^2)
        variance >= -tolerance || error(
            "negative transformed posterior variance for $(spec.key): $variance")
        variances[spec.key] = max(0.0, variance)
    end
    return variances
end

function truth_values(row, aerosol_truth, prior_physical)
    truth = copy(prior_physical)
    merge!(truth, aerosol_truth)
    truth["XCO2"] = parse(Float64, row["xco2_ppm"])
    truth["psurf"] = parse(Float64, row["psurf_hpa"])
    truth["SIF760"] = parse(Float64, row["SIF760"]) *
        SIF_WAVENUMBER_TO_WAVELENGTH_760
    truth["mSIF"] = parse(Float64, row["mSIF"])
    for (state_band, truth_band) in zip(
            BAND_NAMES, ("o2a", "weak", "strong")), order in 0:2
        truth["$(state_band)_surface_P$order"] =
            parse(Float64, row["$(truth_band)_P$order"])
    end
    return truth
end

outcome_label(value) = value == 1 ? "converged_fit_pass" :
    value == 2 ? "converged_fit_fail" :
    value == 3 ? "maximum_iterations" :
    value == 4 ? "maximum_divergences" : "unknown_$value"

function read_run(path, expected_state, expected_class, expected_perturbation)
    return NCDataset(path) do dataset
        get(dataset.attrib, "retrieval_complete", 0) == 1 ||
            error("retrieval is not marked complete: $path")
        state_index = Int(dataset.attrib["truth_state_index"])
        class = Symbol(String(dataset.attrib["measurement_class"]))
        perturbation = Int(dataset.attrib["perturbation_index"])
        state_index == expected_state || error("truth-state mismatch in $path")
        class == expected_class || error("measurement-class mismatch in $path")
        perturbation == expected_perturbation || error(
            "perturbation-index mismatch in $path")

        names = split(String(dataset.attrib["parameter_names"]))
        final_state = Float64.(dataset["final_state"][:])
        length(names) == length(final_state) || error(
            "parameter-name/state mismatch in $path")
        final_xco2 = Float64(dataset["XCO2"][])
        prior_state = Float64.(dataset["a_priori_state"][:])
        prior_xco2 = Float64(dataset["a_priori_XCO2"][])
        posterior = Float64.(dataset["posterior_covariance"][:, :])
        fixed_upper_ppm = Float64(dataset.attrib["fixed_upper_co2_ppm"])

        trials = Int.(dataset["trial_index"][:])
        iterations = Int.(dataset["iteration_index"][:])
        accepted = Bool.(dataset["trial_accepted"][:])
        divergent = Bool.(dataset["trial_divergent"][:])
        d_sigma = Float64.(dataset["d_sigma_sq_scaled"][:])
        evaluation_seconds = Float64.(dataset["evaluation_seconds"][:])
        linear_seconds = Float64.(dataset["linear_algebra_seconds"][:])
        ntrial = length(trials)
        all(length(values) == ntrial for values in
            (iterations, accepted, divergent, d_sigma,
             evaluation_seconds, linear_seconds)) || error(
            "trial-vector lengths disagree in $path")
        trials == collect(1:ntrial) || error("trial indices are not contiguous in $path")
        any(accepted) || error("retrieval has no accepted trial in $path")

        divergence_count = Int(dataset.attrib["divergence_count"])
        divergence_count == count(divergent) || error(
            "stored divergence count disagrees with trial flags in $path")
        rejected_count = count(!, accepted)
        rejected_count == count(divergent) || error(
            "rejected/divergent trial flags disagree in $path")
        terminal_record = findlast(accepted)
        isnothing(terminal_record) && error("no accepted trial in $path")
        terminal_d_sigma = d_sigma[terminal_record]
        convergence_threshold = Float64(dataset.attrib["convergence_threshold"])
        maximum_iterations = Int(dataset.attrib["maximum_iterations"])
        outcome = Int(dataset.attrib["outcome"])
        converged = Bool(dataset.attrib["converged"])
        fit_quality_ok = Bool(dataset.attrib["fit_quality_ok"])
        converged == (outcome in (1, 2)) || error(
            "outcome/convergence flag mismatch in $path")
        final_chi = Float64.(dataset["final_band_reduced_chi_squared"][:])
        length(final_chi) == 3 || error("expected three final band chi-squares")
        chi_threshold = Float64(dataset.attrib["maximum_band_chi_squared"])
        fit_quality_ok == all(<(chi_threshold), final_chi) || error(
            "fit flag disagrees with final band chi-square in $path")
        step_condition_met = terminal_d_sigma < convergence_threshold
        converged == step_condition_met || error(
            "convergence flag disagrees with terminal accepted step in $path")

        final_evaluation_seconds = Float64(
            dataset.attrib["final_evaluation_seconds"])
        # A converged solve evaluates its accepted proposal once more at the
        # terminal state. Outcomes 3/4 reuse an already tabulated trial
        # evaluation, so adding final_evaluation_seconds there would double
        # count it.
        total_evaluation_seconds = sum(evaluation_seconds) +
            (converged ? final_evaluation_seconds : 0.0)
        total_linear_seconds = sum(linear_seconds)

        normalized_draw = Float64.(dataset["normalized_noise_draw"][:])
        if perturbation == UNPERTURBED
            all(iszero, normalized_draw) || error(
                "perturbation 11 is not noiseless in $path")
        end

        final_physical = physical_state(names, final_state, final_xco2)
        prior_physical = physical_state(names, prior_state, prior_xco2)
        posterior_variance = physical_posterior_variances(
            names, final_state, final_xco2, fixed_upper_ppm, posterior)

        return (;
            path=abspath(path), state_index, class, perturbation, names,
            final_state, final_physical, prior_physical, posterior_variance,
            fixed_upper_ppm,
            surface=String(dataset.attrib["surface"]),
            aerosol_case=String(dataset.attrib["aerosol_case"]),
            truth_xco2=Float64(dataset.attrib["truth_xco2_ppm"]),
            normalized_draw,
            trial_count=ntrial,
            accepted_iterations=count(accepted),
            maximum_accepted_iteration=maximum(iterations[accepted]),
            maximum_iterations,
            rejected_trials=rejected_count,
            divergent_trials=count(divergent),
            divergence_count,
            outcome,
            outcome_name=outcome_label(outcome),
            converged,
            fit_quality_ok,
            terminal_d_sigma,
            convergence_threshold,
            step_condition_met,
            final_chi,
            chi_threshold,
            band_fit_condition_met=all(<(chi_threshold), final_chi),
            total_evaluation_seconds,
            total_linear_seconds,
            final_evaluation_seconds,
        )
    end
end

function load_runs(root)
    runs = Dict{Tuple{Int,Symbol,Int},Any}()
    for state in STATES, class in CLASSES, perturbation in EXPECTED_PERTURBATIONS
        path = retrieval_path(root, state, class, perturbation)
        run = read_run(path, state, class, perturbation)
        runs[(state, class, perturbation)] = run
    end

    reference_names = runs[(first(STATES), first(CLASSES), 1)].names
    for run in values(runs)
        run.names == reference_names || error(
            "active parameter layout differs across retrieval products")
        run.truth_xco2 == Float64(run.fixed_upper_ppm) || error(
            "fixed upper CO2 does not match truth XCO2 in $(run.path)")
    end
    for state in STATES, perturbation in EXPECTED_PERTURBATIONS
        corrected = runs[(state, :corrected, perturbation)]
        uncorrected = runs[(state, :uncorrected, perturbation)]
        corrected.surface == uncorrected.surface || error(
            "surface mismatch in state $state perturbation $perturbation")
        corrected.aerosol_case == uncorrected.aerosol_case || error(
            "aerosol-case mismatch in state $state perturbation $perturbation")
        corrected.truth_xco2 == uncorrected.truth_xco2 || error(
            "truth-XCO2 mismatch in state $state perturbation $perturbation")
        corrected.normalized_draw == uncorrected.normalized_draw || error(
            "corrected/uncorrected normalized noise draws differ for " *
            "state $state perturbation $perturbation")
    end
    return runs
end

function state_truths(runs)
    truths = Dict{Int,Dict{String,Float64}}()
    for state in STATES
        reference = runs[(state, :corrected, 1)]
        row = table_row(TRUTH_TABLE, state)
        row["surface"] == reference.surface || error(
            "truth-table/retrieval surface mismatch for state $state")
        row["aerosol_case"] == reference.aerosol_case || error(
            "truth-table/retrieval aerosol mismatch for state $state")
        aerosol = aerosol_aod760(SCENE_COMPONENTS, reference.aerosol_case)
        truth = truth_values(row, aerosol, reference.prior_physical)
        truth["XCO2"] == reference.truth_xco2 || error(
            "truth-table/retrieval XCO2 mismatch for state $state")
        truths[state] = truth
    end
    return truths
end

function write_run_table(path, runs)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# State-001 versus state-009 retrieval run diagnostics.")
        println(io, "# total_evaluation_seconds counts each forward+Jacobian " *
                    "evaluation once, including the separate terminal evaluation " *
                    "only for converged outcomes.")
        println(io, "# state class perturbation noise_case trials " *
                    "accepted_iterations maximum_accepted_iteration " *
                    "maximum_iterations_allowed rejected_trials divergent_trials " *
                    "convergence outcome outcome_name fit_quality_ok " *
                    "terminal_d_sigma_scaled convergence_threshold " *
                    "step_condition_met o2a_reduced_chi2 weak_co2_reduced_chi2 " *
                    "strong_co2_reduced_chi2 band_chi2_threshold " *
                    "band_fit_condition_met total_evaluation_seconds " *
                    "total_linear_algebra_seconds final_evaluation_seconds path")
        for state in STATES, class in CLASSES,
                perturbation in EXPECTED_PERTURBATIONS
            run = runs[(state, class, perturbation)]
            noise_case = perturbation == UNPERTURBED ? "noiseless" : "perturbed"
            @printf(io,
                "%03d %-11s %02d %-9s %d %d %d %d %d %d %d %d %-22s %d %.12e %.12e %d %.12e %.12e %.12e %.12e %d %.12e %.12e %.12e %s\n",
                state, String(class), perturbation, noise_case,
                run.trial_count, run.accepted_iterations,
                run.maximum_accepted_iteration, run.maximum_iterations,
                run.rejected_trials, run.divergent_trials,
                Int(run.converged), run.outcome, run.outcome_name,
                Int(run.fit_quality_ok), run.terminal_d_sigma,
                run.convergence_threshold, Int(run.step_condition_met),
                run.final_chi..., run.chi_threshold,
                Int(run.band_fit_condition_met), run.total_evaluation_seconds,
                run.total_linear_seconds, run.final_evaluation_seconds,
                run.path)
        end
    end
    return path
end

function state_statistics(runs, truths)
    rows = NamedTuple[]
    for state in STATES, class in CLASSES, spec in PARAMETER_SPECS
        noisy_runs = [runs[(state, class, perturbation)] for perturbation in PERTURBED]
        values = [run.final_physical[spec.key] for run in noisy_runs]
        posterior_variances = [
            run.posterior_variance[spec.key] for run in noisy_runs]
        posterior_stds = sqrt.(posterior_variances)
        noiseless = runs[(state, class, UNPERTURBED)]
        truth = truths[state][spec.key]
        push!(rows, (;
            state, class, parameter=spec.key, unit=spec.unit, truth,
            sample_count=length(values),
            ensemble_mean=mean(values),
            ensemble_bias=mean(values) - truth,
            sample_variance=var(values; corrected=true),
            sample_std=std(values; corrected=true),
            mean_posterior_variance=mean(posterior_variances),
            mean_posterior_std=mean(posterior_stds),
            sqrt_mean_posterior_variance=sqrt(mean(posterior_variances)),
            noiseless_value=noiseless.final_physical[spec.key],
            noiseless_bias=noiseless.final_physical[spec.key] - truth,
            noiseless_posterior_variance=noiseless.posterior_variance[spec.key],
            noiseless_posterior_std=sqrt(noiseless.posterior_variance[spec.key]),
        ))
    end
    return rows
end

function write_statistics_table(path, rows)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# Physical compact-state statistics. Perturbations 01:10 " *
                    "define empirical statistics; perturbation 11 is separate.")
        println(io, "# Sample variance/std describe terminal retrieval spread, " *
                    "not OE posterior uncertainty. Posterior quantities use " *
                    "delta-method physical-coordinate transforms.")
        println(io, "# state class parameter units truth n_perturbed " *
                    "ensemble_mean ensemble_bias sample_variance sample_std " *
                    "mean_terminal_posterior_variance " *
                    "mean_terminal_posterior_std " *
                    "sqrt_mean_terminal_posterior_variance noiseless_value " *
                    "noiseless_bias noiseless_posterior_variance " *
                    "noiseless_posterior_std")
        for row in rows
            @printf(io,
                "%03d %-11s %-31s %-22s % .12e %d % .12e % .12e % .12e % .12e % .12e % .12e % .12e % .12e % .12e % .12e % .12e\n",
                row.state, String(row.class), row.parameter, row.unit, row.truth,
                row.sample_count, row.ensemble_mean, row.ensemble_bias,
                row.sample_variance, row.sample_std,
                row.mean_posterior_variance, row.mean_posterior_std,
                row.sqrt_mean_posterior_variance, row.noiseless_value,
                row.noiseless_bias, row.noiseless_posterior_variance,
                row.noiseless_posterior_std)
        end
    end
    return path
end

function range_text(values; digits=2)
    formatter(value) = @sprintf("%.*f", digits, value)
    return "$(formatter(minimum(values)))-$(formatter(maximum(values)))"
end

function outcome_counts(runs_for_group)
    counts = Dict{Int,Int}()
    for run in runs_for_group
        counts[run.outcome] = get(counts, run.outcome, 0) + 1
    end
    return join(("$outcome:$(counts[outcome])" for outcome in sort(collect(keys(counts)))), ",")
end

function stat_lookup(rows, state, class, parameter)
    index = findfirst(row -> row.state == state && row.class == class &&
                             row.parameter == parameter, rows)
    isnothing(index) && error("missing state statistic")
    return rows[index]
end

fmt(value) = @sprintf("%.6g", value)

function write_summary(path, runs, truths, rows, run_table, statistics_table)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# State 001 versus state 009 retrieval comparison")
        println(io)
        println(io, "Generated: `$(now(UTC))` UTC")
        println(io)
        println(io, "This report compares the matched completed retrieval ensembles for " *
                    "urban, 380-ppm, no-SIF truth states 001 (no aerosol) and " *
                    "009 (AOD760 = 0.28). Corrected and uncorrected products " *
                    "are both required for all eleven perturbation indices.")
        println(io)
        println(io, "Perturbations 01:10 define the empirical mean, bias, sample " *
                    "variance, and sample standard deviation. Perturbation 11 " *
                    "is noiseless and is never included in those ensemble statistics. " *
                    "Empirical spread and terminal OE posterior uncertainty are " *
                    "reported separately.")
        println(io)
        println(io, "Noise draws are paired between corrected and uncorrected " *
                    "retrievals of the same truth state and perturbation. They are " *
                    "not paired between states 001 and 009, so comparisons between " *
                    "the two state ensembles are independent; this report performs " *
                    "no paired cross-state test.")
        println(io)
        println(io, "## Convergence and timing")
        println(io)
        println(io, "The state-step test is strict: `d_sigma_sq_scaled < " *
                    "convergence_threshold`. Fit quality requires every band's " *
                    "terminal reduced chi-square to be below the stored per-band " *
                    "threshold; it does not define state-step convergence.")
        println(io)
        println(io, "| State | Class | Noisy converged | Noisy fit pass | Outcomes " *
                    "| Trials mean [range] | Accepted mean [range] / max | " *
                    "Rejected / divergent | Terminal step mean [max] / threshold | " *
                    "Mean final chi2 (O2A, weak, strong) | Eval time mean [range] s | " *
                    "Noiseless: conv/fit, trials, accepted, step, chi2, time s |")
        println(io, "|---:|:---|---:|---:|:---|:---|:---|:---|:---|:---|:---|:---|")
        for state in STATES, class in CLASSES
            noisy = [runs[(state, class, perturbation)] for perturbation in PERTURBED]
            noiseless = runs[(state, class, UNPERTURBED)]
            trial_counts = [run.trial_count for run in noisy]
            accepted_counts = [run.accepted_iterations for run in noisy]
            steps = [run.terminal_d_sigma for run in noisy]
            times = [run.total_evaluation_seconds for run in noisy]
            mean_chi = [mean(run.final_chi[band] for run in noisy) for band in 1:3]
            @printf(io,
                "| %03d | %s | %d/10 | %d/10 | %s | %.2f [%s] | %.2f [%s] / %d | %d / %d | %.4g [%.4g] / %.4g | %.4g, %.4g, %.4g | %.1f [%s] | %d/%d, %d, %d, %.4g, %.4g/%.4g/%.4g, %.1f |\n",
                state, String(class), count(run -> run.converged, noisy),
                count(run -> run.fit_quality_ok, noisy), outcome_counts(noisy),
                mean(trial_counts), range_text(trial_counts; digits=0),
                mean(accepted_counts), range_text(accepted_counts; digits=0),
                first(noisy).maximum_iterations,
                sum(run.rejected_trials for run in noisy),
                sum(run.divergent_trials for run in noisy),
                mean(steps), maximum(steps), first(noisy).convergence_threshold,
                mean_chi..., mean(times), range_text(times; digits=1),
                Int(noiseless.converged), Int(noiseless.fit_quality_ok),
                noiseless.trial_count, noiseless.accepted_iterations,
                noiseless.terminal_d_sigma, noiseless.final_chi...,
                noiseless.total_evaluation_seconds)
        end
        println(io)
        println(io, "## Compact physical-state statistics")
        println(io)
        println(io, "Values below show `mean +/- sample SD` for noisy perturbations, " *
                    "the mean terminal posterior SD, and the noiseless terminal value. " *
                    "The machine-readable table additionally includes sample " *
                    "variances, biases, posterior variances, and noiseless posterior " *
                    "uncertainties for every parameter.")
        for class in CLASSES
            println(io)
            println(io, "### $(uppercasefirst(String(class)))")
            println(io)
            println(io, "| Parameter | Units | Truth 001 | State 001 noisy mean +/- SD | " *
                        "Mean posterior SD | Noiseless 001 | Truth 009 | " *
                        "State 009 noisy mean +/- SD | Mean posterior SD | Noiseless 009 |")
            println(io, "|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|")
            for spec in PARAMETER_SPECS
                one = stat_lookup(rows, 1, class, spec.key)
                nine = stat_lookup(rows, 9, class, spec.key)
                println(io, "| `$(spec.key)` | $(spec.unit) | $(fmt(one.truth)) | " *
                            "$(fmt(one.ensemble_mean)) +/- $(fmt(one.sample_std)) | " *
                            "$(fmt(one.mean_posterior_std)) | " *
                            "$(fmt(one.noiseless_value)) | $(fmt(nine.truth)) | " *
                            "$(fmt(nine.ensemble_mean)) +/- $(fmt(nine.sample_std)) | " *
                            "$(fmt(nine.mean_posterior_std)) | " *
                            "$(fmt(nine.noiseless_value)) |")
            end
        end
        println(io)
        println(io, "## Machine-readable outputs")
        println(io)
        println(io, "- [`$(basename(run_table))`]($(basename(run_table)))")
        println(io, "- [`$(basename(statistics_table))`]($(basename(statistics_table)))")
        println(io)
        println(io, "Aerosol-off profile-center heights are still listed because they " *
                    "remain retrieval coordinates, but they are spectrally " *
                    "non-identifiable when their associated AOD is zero.")
    end
    return path
end

function main(args=ARGS)
    options = parse_arguments(args)
    isnothing(options) && return nothing
    if options.status
        inventory(options.inversion_root)
        return nothing
    end

    require_complete_inputs(options.inversion_root)
    runs = load_runs(options.inversion_root)
    truths = state_truths(runs)
    rows = state_statistics(runs, truths)
    run_table = options.output_prefix * "_runs.dat"
    statistics_table = options.output_prefix * "_state_statistics.dat"
    summary = options.output_prefix * "_summary.md"
    write_run_table(run_table, runs)
    write_statistics_table(statistics_table, rows)
    write_summary(summary, runs, truths, rows, run_table, statistics_table)
    println(run_table)
    println(statistics_table)
    println(summary)
    return nothing
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
