module Round7RetrievalCampaign

using NCDatasets
using Printf
using SHA
using ..OptimalEstimation: OESettings
using ..RetrievalCases
using ..Round5FixedSIF
import ..Round5RetrievalCampaign as R5

export ROUND7_STATE_MODEL, ROUND7_BASE_CHECKPOINT, ROUND7_CODE_FILES,
       round7_identity_values, round7_campaign_identity,
       validate_round7_manifest, validate_round7_prior,
       load_round7_observation, round7_realization, round7_experiments,
       validate_round7_output_root, round7_output_provenance,
       validate_round7_existing_output, round7_smoke_gate, file_sha256

const ROUND7_STATE_MODEL = "round7_fixed_sif_imperfect_correction"
const ROUND7_BASE_CHECKPOINT = "7acab57000dae259207a6760faae156cfa1734f6"
const ROUND7_CODE_FILES = ("Round7RetrievalCampaign.jl",
    "run_round7_imperfect_correction_retrievals.jl", "test_round7_retrieval_campaign.jl")
const NMEASUREMENT = 2742
const NSCENE = 80
const NPERTURBATION = 11
const DEFINITION = Dict{String,Any}(
    "round7_definition_version" => 1,
    "round7_correction_definition" => "Cabannes+RRS-Rayleigh",
    "round7_correction_operation" => "subtract",
    "round7_pressure_geometry_source" => "independent_scene_metadata",
    "round7_noise_policy" => "reuse_round5_uncorrected_noise_and_covariance")

file_sha256(path::AbstractString) = R5.file_sha256(path)
_digest_lines(values) = bytes2hex(sha256(join(values, "\n") * "\n"))
_observation_name(index) = @sprintf("OCO2round7_%03d.nc", index)

function _equal(attributes, name, expected, source)
    haskey(attributes, name) || error("$source is missing $name")
    actual = attributes[name]
    matches = actual isa Real && expected isa Real ?
        isapprox(actual, expected; rtol=2e-12, atol=2e-14) :
        string(actual) == string(expected)
    matches || error("$source: $name=$(repr(actual)); expected $(repr(expected))")
    return actual
end

function _hash(value, name; length=64)
    value isa AbstractString && occursin(Regex("^[0-9a-f]{$length}\$"), value) ||
        error("$name must be $length lowercase hexadecimal characters")
    return String(value)
end

"""
Portable deployment identity, independent of absolute paths. `code_paths` is
the three ROUND7_CODE_FILES in that order. `input_paths` is the ordered list
documented in the runner. Each set digest hashes newline-terminated file
digests; the codeset additionally starts with the pinned git tree digest.
The campaign digest hashes model, checkpoint, codeset, input set, and
observation-manifest digest, each terminated by a newline.
"""
function round7_identity_values(; checkpoint, base_tree_sha256, code_paths,
                                 input_paths, observation_manifest)
    _hash(checkpoint, "checkpoint"; length=40)
    _hash(base_tree_sha256, "base_tree_sha256")
    code = _digest_lines(vcat(base_tree_sha256, file_sha256.(code_paths)))
    inputs = _digest_lines(file_sha256.(input_paths))
    observations = file_sha256(observation_manifest)
    campaign = _digest_lines([ROUND7_STATE_MODEL, checkpoint, code, inputs, observations])
    return Dict(
        "ROUND7_CODE_CHECKPOINT_SHA" => checkpoint,
        "ROUND7_CODESET_SHA256" => code,
        "ROUND7_INPUT_SET_SHA256" => inputs,
        "ROUND7_CAMPAIGN_IDENTITY_SHA256" => campaign,
        "ROUND7_OBSERVATION_MANIFEST_SHA256" => observations)
end

function round7_campaign_identity(environment=ENV; expected=nothing)
    result = Dict{String,Any}("round7_campaign_identity_status" => "complete")
    for name in ("ROUND7_CODE_CHECKPOINT_SHA", "ROUND7_CODESET_SHA256",
                 "ROUND7_INPUT_SET_SHA256", "ROUND7_CAMPAIGN_IDENTITY_SHA256",
                 "ROUND7_OBSERVATION_MANIFEST_SHA256")
        # Accept the unprefixed spelling in the observation-generator contract.
        source = name == "ROUND7_OBSERVATION_MANIFEST_SHA256" &&
            !haskey(environment, name) ? "OBSERVATION_MANIFEST_SHA256" : name
        haskey(environment, source) || error("$source is required")
        value = _hash(environment[source], source;
                      length=endswith(name, "CHECKPOINT_SHA") ? 40 : 64)
        if name == "ROUND7_OBSERVATION_MANIFEST_SHA256" &&
                haskey(environment, "OBSERVATION_MANIFEST_SHA256")
            value == environment["OBSERVATION_MANIFEST_SHA256"] ||
                error("observation manifest identity aliases disagree")
        end
        isnothing(expected) || value == expected[name] ||
            error("$source differs from the actual code/input identity")
        result[lowercase(name)] = value
    end
    return result
end

"""Verify SHA256SUMS, all 80 immutable observations, and provenance.json."""
function validate_round7_manifest(directory; expected_sha256)
    root = realpath(directory)
    manifest = joinpath(root, "SHA256SUMS")
    file_sha256(manifest) == _hash(expected_sha256, "observation manifest SHA256") ||
        error("Round-7 observation manifest checksum mismatch")
    entries = Dict{String,String}()
    for line in eachline(manifest)
        isempty(strip(line)) && continue
        matched = match(r"^([0-9a-f]{64}) [ *](?:\./)?(OCO2round7_[0-9]{3}\.nc|provenance\.json)$", line)
        isnothing(matched) && error("invalid Round-7 SHA256SUMS entry: $line")
        digest, name = matched.captures
        haskey(entries, name) && error("duplicate Round-7 manifest entry: $name")
        entries[name] = digest
    end
    Set(keys(entries)) == Set(vcat(_observation_name.(1:NSCENE), "provenance.json")) ||
        error("Round-7 SHA256SUMS must list exactly states 001:080 plus provenance.json")
    for (name, digest) in entries
        path = joinpath(root, name)
        isfile(path) && !islink(path) || error("missing or symlink observation: $path")
        file_sha256(path) == digest || error("Round-7 observation checksum mismatch: $path")
    end
    return (; directory=root, path=manifest, sha256=expected_sha256, entries)
end

"""Keep the Round-5 validator unchanged; additionally check full fixed fields."""
function validate_round7_prior(prior, map; path, source_prior_path)
    R5.validate_round5_prior(prior, map; path)
    NCDataset(path) do ds
        _equal(ds.attrib, "source_prior_sha256", file_sha256(source_prior_path), path)
        fixed = expand_round5_state(map, zeros(ROUND5_STATE_COUNT))[29:30]
        full_mean = Float64.(ds["xa"][:, :])
        full_cov = Float64.(ds["Sa"][:, :, :])
        size(full_mean) == (34, 4) && size(full_cov) == (34, 34, 4) ||
            error("Round-5 full prior has the wrong shape")
        all(isfinite, full_mean) && all(isfinite, full_cov) || error("nonfinite Round-5 prior")
        all(iszero, full_cov[33:34, :, :]) && all(iszero, full_cov[:, 33:34, :]) ||
            error("fixed SIF prior covariance must be exactly zero")
        for surface in 1:4
            full_mean[33:34, surface] ≈ fixed || error("fixed SIF mean/map mismatch")
            active = R5.expected_active_to_full(map)
            ds["Sa_active"][:, :, surface] == full_cov[active, active, surface] ||
                error("active/full Round-5 prior covariance mismatch")
        end
        NCDataset(source_prior_path) do source
            full_mean[1:32, :] == source["xa"][1:32, :] ||
                error("non-SIF prior means differ from Round-5 controls")
            expected = Float64.(source["Sa"][:, :, :])
            for surface in 1:4, index in (20, 23)
                expected[index, :, surface] .*= 0.1
                expected[:, index, surface] .*= 0.1
            end
            full_cov[1:32, 1:32, :] == expected[1:32, 1:32, :] ||
                error("non-SIF covariance differs from Round-5 tight-UTLS controls")
        end
    end
    return prior
end

function _array(ds, name, shape; dimensions=nothing)
    haskey(ds, name) || error("Round-7 observation is missing $name")
    variable = ds[name]
    size(variable) == shape || error("$name has size $(size(variable)); expected $shape")
    isnothing(dimensions) || dimnames(variable) == dimensions ||
        error("$name has incorrect dimension order $(dimnames(variable))")
    values = Float64.(nomissing(Array(variable), NaN))
    all(isfinite, values) || error("$name has missing/nonfinite values")
    return values
end

# Tight elementwise tolerances, scaled by the operands (including cancellation).
# Inputs and the unchanged writer use Float64 even though RT runs in Float32.
function _arithmetic(actual, expected, label; scale=abs.(expected))
    all(abs.(actual .- expected) .<=
        32eps(Float64) .* max.(abs.(actual), abs.(expected), scale, floatmin(Float64))) ||
        error("Round-7 $label arithmetic mismatch")
end

"""
Read all eleven stored realizations without redrawing or reconstructing noise.
SIF provenance is the version-2 record copied from Round-5 measurements and
controls, not a new SIF release receipt. Matrices are (measurement,perturbation)
in Julia, i.e. (perturbation,measurement) in Python/on disk.
"""
function load_round7_observation(manifest, truth, map; truth_table_sha256,
                                 prior_sha256, analyzer_sha256)
    path = joinpath(manifest.directory, _observation_name(truth.state_index))
    digest = manifest.entries[basename(path)]
    file_sha256(path) == digest || error("Round-7 observation changed after manifest validation")
    observation = NCDataset(path) do ds
        for (key, value) in DEFINITION
            _equal(ds.attrib, key, value, path)
        end
        _equal(ds.attrib, "round7_observation_complete", 1, path)
        for key in (:state_index, :surface, :aerosol_case, :sif_case)
            _equal(ds.attrib, string(key), getproperty(truth, key), path)
        end
        _equal(ds.attrib, "source_truth_table_sha256", truth_table_sha256, path)
        _equal(ds.attrib, "round5_prior_sha256", prior_sha256, path)
        _equal(ds.attrib, "source_analyzer_sha256", analyzer_sha256, path)
        for key in ("source_generator_sha256", "source_grid_sha256")
            _hash(get(ds.attrib, key, ""), key)
        end
        # The pinned Round-5 forward model uses this independent scene geometry;
        # TruthCase does not carry geometry, and no retrieved psurf is substituted.
        for (key, value) in ("psurf_hpa" => 1000.0, "sza_deg" => 30.0,
                             "vza_deg" => 0.0, "relative_azimuth_deg" => 0.0)
            _equal(ds.attrib, key, value, path)
        end
        get(ds.dim, "measurement", 0) == NMEASUREMENT || error("measurement dimension must be 2742")
        get(ds.dim, "perturbation", 0) == NPERTURBATION || error("perturbation dimension must be 11")
        provenance = validated_sif_provenance(ds.attrib, truth.sif_case; source=path)
        # Optional extra generator hashes are retained, never interpreted as instructions.
        attributes = Dict{String,Any}(string(k) => v for (k, v) in pairs(ds.attrib))
        vector(name) = _array(ds, name, (NMEASUREMENT,); dimensions=("measurement",))
        matrix(name) = _array(ds, name, (NMEASUREMENT, NPERTURBATION);
                              dimensions=("measurement", "perturbation"))
        wavelength = vector("wavelength")
        all(>(0), wavelength) || error("nonpositive wavelength")
        starts = _array(ds, "band_start_index", (3,))
        stops = _array(ds, "band_end_index", (3,))
        starts == [1, 935, 1742] && stops == [934, 1741, 2742] ||
            error("Round-7 requires the unchanged Round-5 OCO band grid")
        ranges = UnitRange{Int}[Int(a):Int(b) for (a, b) in zip(starts, stops)]
        all(!isempty, ranges) && reduce(vcat, collect.(ranges)) == collect(1:NMEASUREMENT) ||
            error("band indices must partition all 2742 measurements (one-based)")
        all(r -> all(>(0), diff(wavelength[r])), ranges) || error("wavelengths not increasing within bands")
        _array(ds, "perturbation_index", (NPERTURBATION,);
               dimensions=("perturbation",)) == collect(1:NPERTURBATION) ||
            error("perturbation indices must be 1:11")
        haskey(ds, "random_seed_uint64") || error("missing stored random_seed_uint64")
        seed_variable = ds["random_seed_uint64"]
        size(seed_variable) == (NPERTURBATION,) && dimnames(seed_variable) == ("perturbation",) ||
            error("random_seed_uint64 must be a perturbation vector")
        raw_seeds = Array(seed_variable)
        all(v -> v isa UInt64, raw_seeds) || error("stored seeds must be native UInt64, never Float64")
        seeds = UInt64.(raw_seeds)
        iszero(seeds[11]) || error("perturbation 11 seed must be zero")
        coefficients = _array(ds, "mean_o2a_surface_coefficients", (3,))
        uncorrected = vector("measurement_uncorrected")
        ideal = vector("measurement_ideal_corrected")
        imperfect = vector("measurement_imperfectly_corrected")
        correction = vector("correction_observation")
        noise_std = vector("noise_standard_deviation")
        variance = vector("Se_diagonal")
        all(>(0), noise_std) && all(>(0), variance) || error("noise covariance must be positive")
        _arithmetic(variance, noise_std .^ 2, "covariance")
        _arithmetic(imperfect, uncorrected - correction, "subtraction";
                    scale=abs.(uncorrected) + abs.(correction))
        co2 = 935:NMEASUREMENT
        all(iszero, correction[co2]) && imperfect[co2] == uncorrected[co2] &&
            ideal[co2] == uncorrected[co2] || error("Round-7 correction must leave both CO2 bands exactly unchanged")
        uncorrected_perturbed = matrix("uncorrected_perturbed")
        perturbed = matrix("imperfectly_corrected_perturbed")
        draw = matrix("normalized_noise_draw")
        injected = matrix("injected_measurement_noise")
        _arithmetic(injected, noise_std .* draw, "injected noise")
        _arithmetic(uncorrected_perturbed, uncorrected .+ injected, "uncorrected noise addition";
                    scale=abs.(uncorrected) .+ abs.(injected))
        _arithmetic(perturbed, uncorrected_perturbed .- correction, "perturbed subtraction";
                    scale=abs.(uncorrected_perturbed) .+ abs.(correction))
        # The base writer requires exact equality. Never silently repair a stored draw.
        perturbed == imperfect .+ noise_std .* draw ||
            error("stored corrected realizations violate exact base-writer noise addition")
        perturbed[co2, :] == uncorrected_perturbed[co2, :] ||
            error("Round-7 noisy CO2 bands must remain exactly unchanged")
        all(iszero, draw[:, 11]) && all(iszero, injected[:, 11]) &&
            perturbed[:, 11] == imperfect && uncorrected_perturbed[:, 11] == uncorrected ||
            error("perturbation 11 must contain exact zero noise")
        (; path, sha256=digest, attributes, provenance, wavelength, ranges,
           coefficients, uncorrected, ideal, imperfect, correction,
           noise_std, variance, uncorrected_perturbed, perturbed, draw, injected, seeds)
    end
    R5.validate_round5_realization(map, truth, observation; source=path)
    return observation
end

function round7_realization(observation, index::Integer)
    1 <= index <= NPERTURBATION || throw(ArgumentError("perturbation must be 1:11"))
    return MeasurementRealization(copy(observation.imperfect),
        observation.perturbed[:, index], copy(observation.noise_std),
        copy(observation.variance), observation.draw[:, index],
        copy(observation.wavelength), copy(observation.ranges), copy(observation.provenance))
end

"""Use the base IDs/seeds only as metadata; the runner never generates noise."""
function round7_experiments(truths, directory)
    experiments = build_experiments(truths; validate_inputs=false)
    return [RetrievalExperiment(e.retrieval_index, e.pair_index, e.truth,
        e.noise_index, :corrected, e.random_seed,
        joinpath(directory, _observation_name(e.truth.state_index)),
        joinpath(directory, _observation_name(e.truth.state_index)))
        for e in experiments if e.measurement_class == :corrected]
end

function validate_round7_output_root(path; sif_on)
    root = abspath(path)
    expected = sif_on ? "retrievals_sif" : "retrievals_nosif"
    occursin(r"round7(?:[^0-9]|$)", lowercase(root)) && basename(root) == expected ||
        throw(ArgumentError("output must be an isolated round7 namespace ending in $expected"))
    # A new path can still have an existing symlink ancestor into older products.
    ancestor = root
    while !ispath(ancestor)
        ancestor = dirname(ancestor)
    end
    realpath(ancestor) == ancestor || error("Round-7 output root has a symlink ancestor")
    return root
end

function round7_output_provenance(map; prior_path, observation, campaign_identity,
                                  truth_table_sha256)
    result = R5.round5_output_provenance(map; prior_path)
    result["retrieval_state_model"] = ROUND7_STATE_MODEL
    merge!(result, DEFINITION, campaign_identity)
    merge!(result, Dict{String,Any}(
        "round7_observation_kind" => "imperfectly_corrected",
        "round7_underlying_state_model" => R5.ROUND5_STATE_MODEL,
        "round7_underlying_jacobian_flavor" => R5.ROUND5_JACOBIAN_FLAVOR,
        "round7_prior_policy" => "unchanged_round5_fixed_sif_tight_utls",
        "round7_source_observation" => basename(observation.path),
        "round7_source_observation_sha256" => observation.sha256,
        "source_truth_table_sha256" => truth_table_sha256,
        "round7_noise_realization_source" => "stored_round5_uncorrected_realization",
        "round7_float_type" => "Float32"))
    for name in ("source_analyzer_sha256", "source_grid_sha256", "source_generator_sha256")
        result["round7_observation_" * name] = observation.attributes[name]
    end
    NCDataset(prior_path) do ds
        result["round7_source_prior_sha256"] = ds.attrib["source_prior_sha256"]
    end
    return result
end

function round7_smoke_gate(converged, fit_quality_ok; source="retrieval")
    converged && fit_quality_ok || error(
        "Round-7 smoke gate failed for $source: converged=$converged fit_ok=$fit_quality_ok. " *
        "Imperfect correction can legitimately fail the fit gate; result retained, thresholds unchanged.")
    return true
end

"""Missing is resumable; any existing incomplete/mismatched product is fatal."""
function validate_round7_existing_output(path, experiment, realization, prior, map;
                                         provenance, smoke=false)
    ispath(path) || islink(path) || return false
    isfile(path) && !islink(path) || error("refusing nonregular output: $path")
    NCDataset(path) do ds
        for (key, value) in provenance
            _equal(ds.attrib, key, value, path)
        end
        for (key, value) in ("retrieval_complete" => 1,
                "truth_state_index" => experiment.truth.state_index,
                "perturbation_index" => experiment.noise_index,
                "random_seed_uint64" => string(experiment.random_seed),
                "surface" => experiment.truth.surface,
                "aerosol_case" => experiment.truth.aerosol_case,
                "sif_case" => experiment.truth.sif_case,
                "measurement_class" => "corrected", "state_dimension" => 28,
                "nstreams" => 9, "jacobian_flavor" => R5.ROUND5_JACOBIAN_FLAVOR)
            _equal(ds.attrib, key, value, path)
        end
        matching_sif_provenance(realization.provenance,
            validated_sif_provenance(ds.attrib, experiment.truth.sif_case; source=path))
        for (name, values) in (
                "active_core_parameter_index" => active_core_indices(map),
                "measurement_noiseless" => realization.noiseless,
                "measurement_perturbed" => realization.perturbed,
                "noise_standard_deviation" => realization.noise_std,
                "normalized_noise_draw" => realization.normalized_draw,
                "injected_measurement_noise" => realization.noise_std .* realization.normalized_draw,
                "Se_diagonal" => realization.variance,
                "wavelength" => realization.wavelength_nm,
                "band_start_index" => first.(realization.band_ranges),
                "band_end_index" => last.(realization.band_ranges),
                "a_priori_state" => prior.xa,
                "a_priori_covariance" => prior.Sa)
            haskey(ds, name) && Array(ds[name]) == values ||
                error("$path mismatches current $name; refusing overwrite")
        end
        settings = OESettings()
        for key in (:convergence_threshold, :maximum_iterations,
                    :maximum_divergences, :maximum_band_chi_squared, :initial_gamma)
            _equal(ds.attrib, string(key), getproperty(settings, key), path)
        end
        _equal(ds.attrib, "parameter_names", join(prior.parameter_names, " "), path)
        for (name, shape) in (("final_state", (28,)),
                ("final_forward_model", (NMEASUREMENT,)),
                ("final_jacobian", (NMEASUREMENT, 28)),
                ("final_band_reduced_chi_squared", (3,)))
            _array(ds, name, shape)
        end
        smoke && round7_smoke_gate(get(ds.attrib, "converged", 0) == 1,
            get(ds.attrib, "fit_quality_ok", 0) == 1; source=path)
    end
    return true
end

end # module
