module RetrievalOutput

using Dates
using NCDatasets
using Printf

using ..OptimalEstimation
using ..RetrievalCases

export retrieval_output_path, write_retrieval_result, write_xco2_diagnostics!

const INVERSION_ROOT = RetrievalCases.INVERSION_ROOT

function retrieval_output_path(experiment::RetrievalExperiment;
                               inversion_root::AbstractString=INVERSION_ROOT)
    1 <= experiment.noise_index <= 10 || throw(ArgumentError(
        "production perturbation index must lie in 1:10"))
    directory = joinpath(inversion_root, String(experiment.measurement_class))
    filename = @sprintf(
        "retrieval_state%03d_perturbation%02d.nc",
        experiment.truth.state_index,
        experiment.noise_index)
    return joinpath(directory, filename)
end

function _define_vector(output, name, values, dimension; units="1", long_name="")
    variable = defVar(output, name, Float64, (dimension,))
    variable.attrib["units"] = units
    isempty(long_name) || (variable.attrib["long_name"] = long_name)
    variable[:] = values
    return variable
end

function _timing(record, name)
    haskey(record.model_timing, name) || return NaN
    value = record.model_timing[name]
    return value isa Real ? Float64(value) : NaN
end

function _validate_xco2_diagnostics(diagnostics, ntrial::Integer)
    diagnostics isa NamedTuple || throw(ArgumentError(
        "XCO2 diagnostics must be a named tuple"))
    required = (:a_priori_ppm, :trial_ppm, :final_ppm)
    all(name -> hasproperty(diagnostics, name), required) || throw(ArgumentError(
        "XCO2 diagnostics require fields " * join(string.(required), ", ")))
    trials = Float64.(diagnostics.trial_ppm)
    length(trials) == ntrial || throw(DimensionMismatch(
        "XCO2 trial history has $(length(trials)) entries; expected $ntrial"))
    values = vcat(Float64(diagnostics.a_priori_ppm), trials,
                  Float64(diagnostics.final_ppm))
    all(x -> isfinite(x) && 0 < x < 10_000, values) || throw(ArgumentError(
        "XCO2 diagnostics must contain finite values in (0,10000) ppm"))
    return (; a_priori_ppm=values[1], trial_ppm=trials,
            final_ppm=values[end])
end

"""Add the standard scalar and trial-resolved XCO2 variables to an open dataset."""
function write_xco2_diagnostics!(output, diagnostics, ntrial::Integer)
    values = _validate_xco2_diagnostics(diagnostics, ntrial)

    prior = defVar(output, "a_priori_XCO2", Float64, ())
    prior.attrib["units"] = "ppm"
    prior.attrib["long_name"] =
        "dry-air-column-weighted a priori CO2 volume-mixing ratio"
    prior[] = values.a_priori_ppm

    history = defVar(output, "XCO2_at_trial", Float64, ("trial",))
    history.attrib["units"] = "ppm"
    history.attrib["long_name"] =
        "dry-air-column-weighted CO2 volume-mixing ratio at every evaluated LM trial"
    history[:] = values.trial_ppm

    final = defVar(output, "XCO2", Float64, ())
    final.attrib["units"] = "ppm"
    final.attrib["long_name"] =
        "terminal dry-air-column-weighted CO2 volume-mixing ratio"
    final.attrib["weighting"] =
        "sum(CO2_VMR_layer * dry_air_VCD_layer) / sum(dry_air_VCD_layer)"
    final.attrib["surface_pressure_dependent_weights"] = 1
    final[] = values.final_ppm
    return values
end

"""Write one complete retrieval record to its corrected/uncorrected folder."""
function write_retrieval_result(experiment::RetrievalExperiment,
                                realization::MeasurementRealization,
                                result::OEResult,
                                xa::AbstractVector,
                                Sa::AbstractMatrix,
                                parameter_names::AbstractVector{<:AbstractString};
                                output_path::AbstractString=
                                    retrieval_output_path(experiment),
                                settings::OESettings=OESettings(),
                                solar_spectrum_path::AbstractString="",
                                xco2_diagnostics=nothing,
                                provenance::AbstractDict=Dict{String,Any}(),
                                overwrite::Bool=false)
    isfile(output_path) && !overwrite && throw(ArgumentError(
        "retrieval output already exists: $output_path"))
    nstate = length(result.final_state)
    nmeasurement = length(realization.perturbed)
    ntrial = length(result.records)
    nband = length(result.final_band_chi_squared)
    length(xa) == nstate || throw(DimensionMismatch("xa length differs from result state"))
    size(Sa) == (nstate, nstate) || throw(DimensionMismatch("Sa shape differs from state"))
    length(parameter_names) == nstate || throw(DimensionMismatch(
        "parameter-name count differs from state"))
    size(result.final_jacobian) == (nmeasurement, nstate) || throw(DimensionMismatch(
        "terminal Jacobian has an unexpected shape"))
    size(result.gain_matrix) == (nstate, nmeasurement) || throw(DimensionMismatch(
        "gain matrix has an unexpected shape"))

    mkpath(dirname(output_path))
    isfile(output_path) && rm(output_path)
    NCDataset(output_path, "c") do output
        defDim(output, "state", nstate)
        defDim(output, "state_2", nstate)
        defDim(output, "measurement", nmeasurement)
        defDim(output, "trial", ntrial)
        defDim(output, "band", nband)

        state_history = defVar(output, "state_at_trial", Float64, ("state", "trial"))
        state_history.attrib["long_name"] =
            "active retrieval-coordinate state at every evaluated LM trial"
        state_history[:, :] = hcat((record.state for record in result.records)...)
        isnothing(xco2_diagnostics) ||
            write_xco2_diagnostics!(output, xco2_diagnostics, ntrial)
        step_history = defVar(output, "proposed_state_step", Float64, ("state", "trial"))
        step_history.attrib["long_name"] = "LM step proposed from each evaluated trial"
        step_history[:, :] = hcat((record.proposed_step for record in result.records)...)

        _define_vector(output, "a_priori_state", xa, "state";
                       long_name="active a priori state")
        _define_vector(output, "final_state", result.final_state, "state";
                       long_name="terminal active retrieval state")
        prior_covariance = defVar(output, "a_priori_covariance", Float64,
                                  ("state", "state_2"))
        prior_covariance[:, :] = Sa
        posterior = defVar(output, "posterior_covariance", Float64,
                           ("state", "state_2"))
        posterior[:, :] = result.posterior_covariance
        averaging = defVar(output, "averaging_kernel", Float64,
                           ("state", "state_2"))
        averaging[:, :] = result.averaging_kernel

        jacobian = defVar(output, "final_jacobian", Float64,
                          ("measurement", "state"))
        jacobian.attrib["long_name"] = "instrument-processed terminal Jacobian K"
        jacobian[:, :] = result.final_jacobian
        gain = defVar(output, "gain_matrix", Float64, ("state", "measurement"))
        gain.attrib["long_name"] = "terminal OE gain matrix"
        gain[:, :] = result.gain_matrix

        _define_vector(output, "wavelength", realization.wavelength_nm, "measurement";
                       units="nm")
        _define_vector(output, "measurement_noiseless", realization.noiseless,
                       "measurement"; units="mW m-2 sr-1 nm-1")
        _define_vector(output, "measurement_perturbed", realization.perturbed,
                       "measurement"; units="mW m-2 sr-1 nm-1")
        _define_vector(output, "normalized_noise_draw", realization.normalized_draw,
                       "measurement"; units="1")
        _define_vector(output, "noise_standard_deviation", realization.noise_std,
                       "measurement"; units="mW m-2 sr-1 nm-1")
        variance = _define_vector(output, "Se_diagonal", realization.variance,
                                  "measurement";
                                  units="(mW m-2 sr-1 nm-1)^2")
        variance.attrib["frozen_during_retrieval"] = 1
        _define_vector(output, "final_forward_model", result.final_measurement,
                       "measurement"; units="mW m-2 sr-1 nm-1")
        _define_vector(output, "final_residual",
                       result.final_measurement - realization.perturbed,
                       "measurement"; units="mW m-2 sr-1 nm-1")

        band_index = zeros(Int8, nmeasurement)
        band_start = Int32[]
        band_end = Int32[]
        for (iband, range) in enumerate(realization.band_ranges)
            band_index[range] .= iband
            push!(band_start, first(range))
            push!(band_end, last(range))
        end
        band_variable = defVar(output, "band_index", Int8, ("measurement",))
        band_variable.attrib["flag_values"] = Int8[1, 2, 3]
        band_variable.attrib["flag_meanings"] = "o2a weak_co2 strong_co2"
        band_variable[:] = band_index
        start_variable = defVar(output, "band_start_index", Int32, ("band",))
        start_variable.attrib["index_convention"] = "one-based inclusive"
        start_variable[:] = band_start
        end_variable = defVar(output, "band_end_index", Int32, ("band",))
        end_variable.attrib["index_convention"] = "one-based inclusive"
        end_variable[:] = band_end

        function trial_vector(name, values; long_name="")
            variable = defVar(output, name, Float64, ("trial",))
            isempty(long_name) || (variable.attrib["long_name"] = long_name)
            variable[:] = values
            return variable
        end
        trial_number = defVar(output, "trial_index", Int16, ("trial",))
        trial_number[:] = Int16[record.trial for record in result.records]
        iteration_number = defVar(output, "iteration_index", Int16, ("trial",))
        iteration_number[:] = Int16[record.iteration for record in result.records]
        accepted = defVar(output, "trial_accepted", Int8, ("trial",))
        accepted.attrib["flag_values"] = Int8[0, 1]
        accepted.attrib["flag_meanings"] = "rejected accepted"
        accepted[:] = Int8[record.accepted for record in result.records]
        divergent = defVar(output, "trial_divergent", Int8, ("trial",))
        divergent[:] = Int8[record.divergent for record in result.records]
        trial_vector("gamma", [record.gamma for record in result.records])
        trial_vector("reduction_ratio",
                     [record.reduction_ratio for record in result.records])
        trial_vector("measurement_cost",
                     [record.measurement_cost for record in result.records])
        trial_vector("prior_cost", [record.prior_cost for record in result.records])
        trial_vector("total_cost", [record.total_cost for record in result.records])
        trial_vector("forecast_cost", [record.forecast_cost for record in result.records])
        trial_vector("d_sigma_sq", [record.d_sigma_sq for record in result.records])
        trial_vector("d_sigma_sq_scaled",
                     [record.d_sigma_sq_scaled for record in result.records])
        trial_vector("evaluation_seconds",
                     [record.evaluation_seconds for record in result.records])
        trial_vector("linear_algebra_seconds",
                     [record.linear_algebra_seconds for record in result.records])
        trial_vector("forward_seconds",
                     [_timing(record, :forward_seconds) for record in result.records])
        trial_vector("jacobian_seconds",
                     [_timing(record, :jacobian_seconds) for record in result.records])
        trial_vector("instrument_seconds",
                     [_timing(record, :instrument_seconds) for record in result.records])
        trial_vector("model_build_seconds",
                     [_timing(record, :model_build_seconds) for record in result.records])
        trial_vector("rt_linearized_seconds",
                     [_timing(record, :rt_linearized_seconds) for record in result.records])
        effective_count = defVar(output, "effective_parameter_count", Int16,
                                 ("trial",))
        effective_count[:] = Int16[
            record.effective_parameter_count for record in result.records]

        chi_history = defVar(output, "band_reduced_chi_squared", Float64,
                             ("band", "trial"))
        chi_history[:, :] = hcat((record.band_chi_squared for record in result.records)...)
        final_chi = defVar(output, "final_band_reduced_chi_squared", Float64,
                           ("band",))
        final_chi[:] = result.final_band_chi_squared

        truth = experiment.truth
        output.attrib["truth_state_index"] = truth.state_index
        output.attrib["perturbation_index"] = experiment.noise_index
        output.attrib["pair_index"] = experiment.pair_index
        output.attrib["retrieval_index"] = experiment.retrieval_index
        output.attrib["measurement_class"] = String(experiment.measurement_class)
        output.attrib["surface"] = String(truth.surface)
        output.attrib["aerosol_case"] = String(truth.aerosol_case)
        output.attrib["truth_xco2_ppm"] = truth.xco2_ppm
        output.attrib["fixed_upper_co2_layers"] = "1:4"
        output.attrib["fixed_upper_co2_ppm"] = truth.xco2_ppm
        output.attrib["fixed_upper_co2_source"] = "uniform truth-scene XCO2"
        output.attrib["random_seed_uint64"] = string(experiment.random_seed)
        output.attrib["random_distribution"] = "Uniform(-sqrt(3),sqrt(3))"
        output.attrib["parameter_names"] = join(parameter_names, " ")
        output.attrib["state_dimension"] = nstate
        output.attrib["nstreams"] = 9
        output.attrib["jacobian_flavor"] = "OCO_RRS_synth"
        output.attrib["retrieval_forward_scattering"] = "noRS"
        output.attrib["solar_source"] =
            "solar.out high-resolution transmission * Planck(5777 K)"
        isempty(solar_spectrum_path) ||
            (output.attrib["source_solar_spectrum"] = abspath(solar_spectrum_path))
        output.attrib["surface_legendre_coordinate"] =
            "truth complete band before supplemental convolution shoulders"
        output.attrib["convergence_threshold"] = settings.convergence_threshold
        output.attrib["maximum_iterations"] = settings.maximum_iterations
        output.attrib["maximum_divergences"] = settings.maximum_divergences
        output.attrib["maximum_band_chi_squared"] = settings.maximum_band_chi_squared
        output.attrib["initial_gamma"] = settings.initial_gamma
        output.attrib["outcome"] = result.outcome
        output.attrib["converged"] = Int(result.converged)
        output.attrib["fit_quality_ok"] = Int(result.fit_quality_ok)
        output.attrib["divergence_count"] = result.divergence_count
        output.attrib["final_evaluation_seconds"] = result.final_evaluation_seconds
        output.attrib["source_measurement"] = abspath(experiment.measurement_path)
        output.attrib["source_noise_covariance"] = abspath(experiment.noise_path)
        for (key, value) in provenance
            key_string = String(key)
            haskey(output.attrib, key_string) && throw(ArgumentError(
                "provenance key would replace retrieval metadata: $key_string"))
            output.attrib[key_string] = value
        end
        output.attrib["created_utc"] = string(now(UTC))
        output.attrib["retrieval_complete"] = 1
    end
    return output_path
end

end # module RetrievalOutput
