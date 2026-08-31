#!/usr/bin/env julia

"""
Monitor the clean ABSCO truth-map campaign and email when retrieval inputs are
ready.

"Ready" is deliberately stricter than producer-process completion.  This
driver waits for all 64 high-resolution scenes, validates their ABSCO/profile
provenance and finite spectral arrays, then runs the shared closure and
retrieval-input pipeline:

1. truth/retrieval forward and Jacobian closure;
2. instrument-operator and retrieval-state unit tests;
3. Mueller processing, Gaussian convolution, and OCO-grid resampling;
4. validation of all 64 synthetic OCO radiance files;
5. noise-model unit tests and generation of all 64 covariance files;
6. validation of every covariance and the optimal-estimation output schema.

An email is sent only after every stage succeeds.  A failed producer or
validation stage sends a separate NOT READY notification.

Environment controls:

- `NOTIFY_EMAIL` (default: `git config user.email`)
- `POLL_SECONDS` (default: 300)
- `STALE_HOURS` (default: 12)
- `CUDA_DEVICE` (default: 1; used by the closure rerun)
- `FORCE_MONITOR=1` ignores an existing ready/failed sentinel
"""

using Dates
using JLD2
using NCDatasets
using Printf
using Sockets: gethostname

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const TRUTH_ROOT = joinpath(ROOT, "RRS_XCO2", "truth_map")
const AEROSOL_ROOT = joinpath(TRUTH_ROOT, "aerosol_chunked")
const INSTRUMENT_ROOT = joinpath(ROOT, "RRS_XCO2", "inversion", "instrument")
const INVERSION_ROOT = joinpath(ROOT, "RRS_XCO2", "inversion")
const STATUS_PATH = joinpath(TRUTH_ROOT, "ABSCO_PIPELINE_STATUS.md")
const READY_SENTINEL = joinpath(TRUTH_ROOT, ".absco_retrieval_ready")
const FAILED_SENTINEL = joinpath(TRUTH_ROOT, ".absco_pipeline_failed")
const PRODUCER_LOGS = (
    joinpath(AEROSOL_ROOT, "o2_exact_grid_wurst.log"),
    joinpath(TRUTH_ROOT, "o2_exact_grid_curry.log"),
)
const CHECKPOINTS = (
    joinpath(AEROSOL_ROOT, "o2_exact_grid_aerosol_1_64_checkpoint.jld2"),
    joinpath(TRUTH_ROOT, "o2_exact_grid_none_1_64_checkpoint.jld2"),
)
const EXPECTED_CLEAR_UNITS = 8
const EXPECTED_AEROSOL_UNITS = 8 * cld(2735, 64)
const POLL_SECONDS = parse(Int, get(ENV, "POLL_SECONDS", "300"))
const STALE_HOURS = parse(Float64, get(ENV, "STALE_HOURS", "12"))
const FORCE_MONITOR = lowercase(get(ENV, "FORCE_MONITOR", "0")) in
                      ("1", "true", "yes", "on")

const NO_AEROSOL_STATES = vcat(collect(1:8), collect(17:24),
                               collect(33:40), collect(49:56))
const AEROSOL_STATES = vcat(collect(9:16), collect(25:32),
                            collect(41:48), collect(57:64))
const BAND_VARIABLES = (
    "radiance_rayleigh_o2a",
    "radiance_cabannes_o2a",
    "radiance_rrs_o2a",
    "radiance_rayleigh_weak_co2",
    "radiance_rayleigh_strong_co2",
)

function configured_email()
    configured = strip(get(ENV, "NOTIFY_EMAIL", ""))
    !isempty(configured) && return configured
    email = try
        strip(read(`git -C $ROOT config --get user.email`, String))
    catch
        ""
    end
    isempty(email) && error("NOTIFY_EMAIL is unset and git user.email is unavailable")
    return email
end

timestamp() = string(now(UTC), " UTC")

function write_status(status::AbstractString, details::AbstractString="")
    open(STATUS_PATH, "w") do io
        println(io, "# ABSCO truth-to-retrieval pipeline status")
        println(io)
        println(io, "- Updated: `$(timestamp())`")
        println(io, "- Status: **$status**")
        isempty(details) || begin
            println(io)
            println(io, details)
        end
    end
end

function send_email(subject::AbstractString, body::AbstractString, recipient)
    @info "sending email notification" recipient subject
    command = `/usr/bin/mail -s $subject $recipient`
    run(pipeline(command; stdin=IOBuffer(body)))
end

scene_path(index::Integer, aerosol::Bool) = joinpath(
    aerosol ? AEROSOL_ROOT : TRUTH_ROOT,
    @sprintf("hiressim_%03d.nc", index))

function scene_is_complete(path)
    isfile(path) || return false
    try
        return NCDataset(path) do dataset
            get(dataset.attrib, "simulation_complete", 0) == 1 &&
            get(dataset.attrib, "o2_absco_regeneration_complete", 0) == 1 &&
            get(dataset.attrib, "o2_truth_regenerated", 0) == 1
        end
    catch
        # Scene files are created before their arrays are filled. An open/read
        # collision during an active write is therefore incomplete, not yet a
        # validation failure.
        return false
    end
end

function finite_truth_array(dataset, name)
    haskey(dataset, name) || error("completed truth scene is missing $name")
    values = Float64.(nomissing(dataset[name][:, :], NaN))
    size(values, 1) == 3 || error("$name does not contain I,Q,U rows")
    all(isfinite, values) || error("$name contains non-finite values")
    maximum(abs, values) < 1e30 || error("$name contains unwritten fill values")
    return values
end

function validate_scene(path, expected_index; aerosol::Bool)
    NCDataset(path) do dataset
        get(dataset.attrib, "simulation_complete", 0) == 1 ||
            error("scene is not marked complete: $path")
        Int(dataset.attrib["state_index"]) == expected_index ||
            error("state index mismatch in $path")
        String(dataset.attrib["spectroscopy_database"]) == "ABSCO" ||
            error("non-ABSCO spectroscopy in $path")
        String(dataset.attrib["spectroscopy_version"]) == "5.2" ||
            error("wrong ABSCO version in $path")
        Float64(dataset.attrib["o2_vmr"]) == 0.21 ||
            error("wrong O2 VMR in $path")
        Float64(dataset.attrib["psurf_hpa"]) == 1000.0 ||
            error("wrong surface pressure in $path")
        Int(dataset.attrib["atmospheric_layers"]) == 16 ||
            error("wrong atmospheric layer count in $path")
        for key in ("o2_absco_lut", "o2_h2o_lut",
                    "weak_h2o_absco_lut", "weak_co2_absco_lut",
                    "strong_h2o_absco_lut", "strong_co2_absco_lut")
            lut = String(dataset.attrib[key])
            isfile(lut) || error("recorded spectroscopy table is missing: $lut")
        end
        get(dataset.attrib, "co2_absco_regeneration_complete", 0) == 1 ||
            error("scene lacks completed CO2-only ABSCO regeneration: $path")
        get(dataset.attrib, "o2_absco_regeneration_complete", 0) == 1 ||
            error("scene lacks completed exact-grid O2 regeneration: $path")
        get(dataset.attrib, "o2_truth_regenerated", 0) == 1 ||
            error("scene lacks regenerated-O2 provenance: $path")
        get(dataset.attrib, "o2_truth_reused", 1) == 0 ||
            error("scene is still marked as reused O2 truth: $path")
        Int(dataset.attrib["o2_core_grid_version"]) == 2 ||
            error("scene predates exact O2 core-grid construction: $path")
        if aerosol
            get(dataset.attrib, "chunked_simulation_complete", 0) == 1 ||
                error("aerosol scene lacks chunked completion marker: $path")
        end
        for name in BAND_VARIABLES
            finite_truth_array(dataset, name)
        end
    end
    return true
end

function validate_wavelength_file(path)
    isfile(path) || error("missing truth wavelength file: $path")
    NCDataset(path) do dataset
        String(dataset.attrib["spectroscopy_database"]) == "ABSCO" ||
            error("wavelength file lacks ABSCO provenance: $path")
        for (band, expected) in (("o2a", 2735), ("weak_co2", 1281),
                                 ("strong_co2", 995))
            nu = Float64.(dataset["$(band)_wavenumber"][:])
            wavelength = Float64.(dataset["$(band)_wavelength"][:])
            length(nu) == expected || error("wrong $band grid length in $path")
            length(wavelength) == expected || error(
                "wrong $band wavelength length in $path")
            all(isfinite, nu) && all(isfinite, wavelength) ||
                error("non-finite $band grid in $path")
            all(diff(nu) .> 0) || error("non-monotonic $band wavenumber grid in $path")
        end
    end
end

function truth_progress()
    completed_clear = 0
    completed_aerosol = 0
    for index in NO_AEROSOL_STATES
        path = scene_path(index, false)
        if scene_is_complete(path)
            validate_scene(path, index; aerosol=false)
            completed_clear += 1
        end
    end
    for index in AEROSOL_STATES
        path = scene_path(index, true)
        if scene_is_complete(path)
            validate_scene(path, index; aerosol=true)
            completed_aerosol += 1
        end
    end
    return completed_clear, completed_aerosol
end

function checkpoint_progress()
    function completed_units(path)
        isfile(path) || return 0
        try
            saved = load(path)
            return length(get(saved, "completed", String[]))
        catch
            # Atomic producer writes use a temporary file followed by rename,
            # so a transient read failure should be rare. Treat it as no new
            # progress and retry at the next polling interval.
            return 0
        end
    end
    return (completed_units(CHECKPOINTS[2]), completed_units(CHECKPOINTS[1]))
end

function tail_text(path; maximum_bytes=200_000)
    isfile(path) || return ""
    open(path, "r") do io
        seek(io, max(filesize(path) - maximum_bytes, 0))
        return read(io, String)
    end
end

function check_producer_errors()
    pattern = r"(?m)^(ERROR:|LoadError:|OutOfMemoryError|CUDA error)"
    for path in PRODUCER_LOGS
        text = tail_text(path)
        occursin(pattern, text) && error(
            "truth producer reported an error in $path\n" * last(text, min(4000, length(text))))
    end
end

function latest_producer_activity()
    paths = String[PRODUCER_LOGS..., CHECKPOINTS...]
    append!(paths, (scene_path(index, false) for index in NO_AEROSOL_STATES))
    append!(paths, (scene_path(index, true) for index in AEROSOL_STATES))
    existing = filter(isfile, paths)
    isempty(existing) && return 0.0
    return maximum(stat(path).mtime for path in existing)
end

function run_stage(name, relative_script; environment=Dict{String,String}())
    script = joinpath(ROOT, relative_script)
    isfile(script) || error("missing readiness-stage script: $script")
    @info "starting retrieval-readiness stage" name script
    write_status("VALIDATING", "Current stage: `$name`.")
    command = `$(Base.julia_cmd()) --project=$ROOT $script`
    run(addenv(command, environment))
    @info "completed retrieval-readiness stage" name
end

function run_readiness_pipeline()
    validate_wavelength_file(joinpath(TRUTH_ROOT, "sim_wavelength.nc"))
    validate_wavelength_file(joinpath(AEROSOL_ROOT, "sim_wavelength.nc"))

    device = get(ENV, "CUDA_DEVICE", "1")
    run_stage("exact-grid A-band regeneration and CO2 completion",
              "RRS_XCO2/scripts/validate_regenerated_truth.jl")
    run_stage("truth/retrieval ABSCO closure",
              "RRS_XCO2/inversion/validate_truth_forward_closure.jl";
              environment=Dict("CUDA_DEVICE" => device))
    run_stage("instrument forward/Jacobian operator tests",
              "RRS_XCO2/inversion/instrument/test_synthetic_oco2.jl")
    run_stage("retrieval state/profile/source mapping tests",
              "RRS_XCO2/inversion/test_forward_state_mapping.jl")
    run_stage("synthetic OCO radiance generation",
              "RRS_XCO2/inversion/instrument/process_truth_map.jl";
              environment=Dict("FORCE" => "1"))
    run_stage("synthetic OCO radiance validation",
              "RRS_XCO2/inversion/instrument/validate_oco_radiances.jl")
    run_stage("OCO noise-model unit tests",
              "RRS_XCO2/inversion/instrument/test_oco2_noise.jl")
    run_stage("measurement covariance generation",
              "RRS_XCO2/inversion/instrument/generate_noise_covariances.jl";
              environment=Dict("FORCE" => "1"))
    run_stage("measurement covariance validation",
              "RRS_XCO2/inversion/instrument/validate_noise_covariances.jl")
    run_stage("optimal-estimation and output-schema tests",
              "RRS_XCO2/inversion/test_optimal_estimation.jl")
end

function main()
    recipient = configured_email()
    if !FORCE_MONITOR
        isfile(READY_SENTINEL) && error(
            "ready notification was already sent; use FORCE_MONITOR=1 to rerun")
        isfile(FAILED_SENTINEL) && error(
            "a previous monitor failed; inspect it and use FORCE_MONITOR=1 to rerun")
    end
    FORCE_MONITOR && foreach(path -> isfile(path) && rm(path),
                             (READY_SENTINEL, FAILED_SENTINEL))
    mkpath(AEROSOL_ROOT)
    last_counts = (-1, -1)
    last_units = (-1, -1)
    write_status("WAITING FOR TRUTH", "Notification recipient: `$recipient`.")

    try
        while true
            check_producer_errors()
            counts = truth_progress()
            units = checkpoint_progress()
            if counts != last_counts || units != last_units
                @info "truth progress" no_aerosol=counts[1] aerosol=counts[2] clear_units=units[1] aerosol_units=units[2]
                write_status(
                    "WAITING FOR TRUTH",
                    "Completed scenes: no aerosol `$(counts[1])/32`; " *
                    "aerosol `$(counts[2])/32`.\n\n" *
                    "Atomic production units: clear `$(units[1])/$EXPECTED_CLEAR_UNITS`; " *
                    "aerosol `$(units[2])/$EXPECTED_AEROSOL_UNITS`.\n\n" *
                    "Notification recipient: `$recipient`.")
                last_counts = counts
                last_units = units
            end
            counts == (32, 32) && break

            latest = latest_producer_activity()
            if latest > 0 && time() - latest > STALE_HOURS * 3600
                error("truth production has shown no file activity for " *
                      "more than $STALE_HOURS hours")
            end
            sleep(POLL_SECONDS)
        end

        @info "all truth scenes complete; starting retrieval-readiness pipeline"
        run_readiness_pipeline()
        completion = timestamp()
        body = """
The clean RRS-XCO2 ABSCO truth and retrieval-input pipeline completed successfully at $completion.

All 64 high-resolution truth scenes, 64 Mueller/convolved/resampled OCO radiance files, and 64 diagonal noise-covariance files passed validation. The truth/retrieval closure and retrieval-facing unit tests also passed.

New corrected and uncorrected retrieval runs are ready to start.

Closure report:
$INVERSION_ROOT/retrieval_setup/ABSCO_CLOSURE.md

Pipeline status:
$STATUS_PATH
"""
        send_email("RRS-XCO2 ready for new retrieval runs", body, recipient)
        open(READY_SENTINEL, "w") do io
            println(io, completion)
            println(io, recipient)
        end
        write_status("READY FOR RETRIEVALS", body)
        @info "retrieval-ready notification sent" recipient
    catch exception
        failure = sprint(showerror, exception, catch_backtrace())
        body = """
The RRS-XCO2 ABSCO truth-to-retrieval pipeline is NOT ready.

Failure time: $(timestamp())
Host: $(gethostname())

$failure

Inspect the monitor and producer logs under:
$TRUTH_ROOT
"""
        write_status("FAILED — NOT READY", "```text\n$failure\n```")
        try
            send_email("RRS-XCO2 pipeline failed — not ready", body, recipient)
        catch mail_exception
            @error "failure email could not be sent" exception=mail_exception
        end
        open(FAILED_SENTINEL, "w") do io
            println(io, timestamp())
            println(io, failure)
        end
        rethrow()
    end
end

abspath(PROGRAM_FILE) == abspath(@__FILE__) && main()
