# Local-only test suite — NOT part of the CI unit suite (test/runtests.jl).
#
# Run locally with:
#     julia --project=test test/local/runtests.jl
#
# Covers everything that needs resources the GitHub CI runners lack:
#   - GPU / Metal tests (test/local/gpu/) — Mie kernels run on the KA CPU
#     backend always; the CUDA-only tests are guarded and skip without a device.
#   - data/feature-dependent configs (test/local/test_parameters/):
#       * H2O-override scenes — declare H2O in the molecule list + an explicit
#         VMR, a config-level override NOT yet implemented (H2O is driven by
#         atmospheric_profile.q). Parked until the override lands.
#       * *_newLUT* / *SIF_grid* / O2Parameters* — need external ABSCO/HITRAN
#         LUTs (e.g. VSMARTMOM_HITRAN_LUT_DIR or /…/data paths) not in the repo.
#       * legacy reference scenes.
# Every config either runs, is skipped for a declared absent resource, or is a
# documented known failure checked with `@test_broken`. Skips:
#   * an absolute file path or `${ENV:VAR}` reference that does not resolve
#     (decided before parsing);
#   * a GPU architecture when no functional CUDA device is present (decided
#     after parsing, before model construction, so parse errors still fail).
# Any other parse/build error (MethodError, invalid configuration, solver
# failure) is a test error. Set VSMARTMOM_LOCAL_DATA_REQUIRED=true to require
# that at least one config actually ran.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

const _LOCAL_CUDA_OK = try
    @eval using CUDA
    CUDA.functional()
catch
    false
end

# Documented in local/README.md: the config-level H2O override is parked, and
# the legacy reference scenes predate the current schema. Each entry pins the
# exact current failure; an unexpected pass or a different error fails.
const _LOCAL_KNOWN_FAILURES = Dict(
    "2BandParameters.yaml"             => (:parse, ArgumentError, "H2O must not appear in"),
    "3BandParameters.yaml"             => (:parse, ArgumentError, "H2O must not appear in"),
    "3BandParameters_canopy.yaml"      => (:parse, ArgumentError, "H2O must not appear in"),
    "S2BandParameters.yaml"            => (:parse, ArgumentError, "H2O must not appear in"),
    "ParamsEMIT_MODTRANcomp.yaml"      => (:parse, ArgumentError, "H2O must not appear in"),
    "ParamsEMIT_MODTRANcomp_alb02.yaml"=> (:parse, ArgumentError, "H2O must not appear in"),
    "ThreeBandsParameters.yaml"        => (:parse, ArgumentError, "radiative_transfer/spec_bands"),
    "O2ParametersVS.yaml"              => (:build, BoundsError, "AbstractSurfaceType} at index [2]"),
)

"""Unresolved absolute paths and `\${ENV:VAR}` references in non-comment YAML text."""
function _missing_local_resources(path)
    absent = String[]
    for line in eachline(path)
        code = first(split(line, '#'; limit = 2))
        for m in eachmatch(r"\$\{ENV:([A-Za-z_][A-Za-z0-9_]*)\}([^\"'\s,\]]*)", code)
            var, rest = m.captures
            if !haskey(ENV, var)
                push!(absent, "ENV[$var]")
            elseif !isempty(rest) && !ispath(ENV[var] * rest)
                push!(absent, ENV[var] * rest)
            end
        end
        for m in eachmatch(r"[\"'\s\[](/[^\"'\s,\]]+)", code)
            ispath(m.captures[1]) || push!(absent, m.captures[1])
        end
    end
    return unique(absent)
end

_error_text(err) = sprint(showerror, err)

function _check_known_failure(path, (stage, errtype, fragment))
    # A build-stage entry must parse; a parse error there is a test error.
    params = stage === :build ? parameters_from_yaml(path) : nothing
    err = try
        stage === :parse ? parameters_from_yaml(path) : model_from_parameters(params)
        nothing
    catch e
        e
    end
    @test_broken err === nothing
    @test err === nothing || (err isa errtype && occursin(fragment, _error_text(err)))
end

@testset "local (GPU + data-dependent)" begin
    # ── GPU / Metal (the included runner manages its own cwd) ──────────────
    include(joinpath(@__DIR__, "gpu", "runtests.jl"))

    # ── data/feature-dependent config parity ───────────────────────────────
    params_dir = joinpath(@__DIR__, "test_parameters")
    orig_cwd = pwd()
    cd(normpath(joinpath(@__DIR__, "..")))   # test/ root, so config data paths resolve
    try
        @testset "data-dependent configs (forward/lin parity where buildable)" begin
            n_ok = 0
            skipped = String[]
            for path in filter(e -> endswith(e, ".yaml"), readdir(params_dir; join = true))
                cfg = basename(path)
                @testset "$cfg" begin
                    absent = _missing_local_resources(path)
                    known = get(_LOCAL_KNOWN_FAILURES, cfg, nothing)
                    # Parse-stage known failures do not read data files.
                    if known !== nothing && (known[1] === :parse || isempty(absent))
                        _check_known_failure(path, known)
                    elseif !isempty(absent)
                        @info "skipped: missing local data" config = cfg missing = absent
                        push!(skipped, cfg)
                    else
                        params = parameters_from_yaml(path)
                        if params.architecture isa vSmartMOM.Architectures.GPU && !_LOCAL_CUDA_OK
                            @info "skipped: GPU config without a functional CUDA device" config = cfg
                            push!(skipped, cfg)
                        else
                            fwd = model_from_parameters(params)
                            lin_fwd, _ = model_from_parameters(LinMode(), params)
                            @test CoreRT.m_max_bands(fwd) == CoreRT.m_max_bands(lin_fwd)
                            n_ok += 1
                        end
                    end
                end
            end
            @info "local configs run in this environment" n_ok skipped
            if get(ENV, "VSMARTMOM_LOCAL_DATA_REQUIRED", "false") == "true"
                @test n_ok > 0
            end
        end
    finally
        cd(orig_cwd)
    end
end
