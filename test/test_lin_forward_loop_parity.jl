# Phase B — forward / lin Fourier loop-bound parity
# =========================================================================
#
# Asserts that `model_from_parameters(...)` and
# `model_from_parameters(LinMode(), ...)` produce identical
# `m_max_bands` (per-band Fourier loop bound, order semantics) for
# every YAML in `test/test_parameters/`. Without this guard, the
# forward and lin paths could silently disagree on which Fourier
# moments to include — exactly the latent bug the
# `_derive_m_max_bands` helper unifies.
#
# Every config under test/test_parameters/ MUST parse, build forward, AND
# build linearized — a failure in any of those is now a HARD test failure,
# not a silent skip (silent skips let stale/broken configs hide; see the
# 14 configs that had silently rotted before the test/setup split). Configs
# that need external data, external LUTs (VSMARTMOM_HITRAN_LUT_DIR), or the
# not-yet-implemented H2O-override feature live under
# test/setup/test_parameters/ and run via test/setup/runtests.jl instead.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

const _PARAMS_DIR = joinpath(@__DIR__, "test_parameters")

_yaml_paths() = filter(p -> endswith(p, ".yaml"), readdir(_PARAMS_DIR; join = true))

@testset "Phase B — forward/lin m_max_bands parity" begin
    n_compared = 0
    for path in _yaml_paths()
        cfg = basename(path)
        params = try
            parameters_from_yaml(path)
        catch err
            @error "Config in test_parameters/ failed to parse — fix it or move it to test/setup/test_parameters/" config=cfg exception=(err, catch_backtrace())
            @test false
            continue
        end

        fwd = try
            model_from_parameters(params)
        catch err
            @error "Config in test_parameters/ failed to build forward — fix it or move it to test/setup/test_parameters/" config=cfg exception=(err, catch_backtrace())
            @test false
            continue
        end

        # lin-build failures are real bugs — not wrapped.
        lin_model, _ = model_from_parameters(LinMode(), params)

        fwd_m = CoreRT.m_max_bands(fwd)
        lin_m = CoreRT.m_max_bands(lin_model)
        @test fwd_m == lin_m
        @test CoreRT.n_fourier_moments_bands(fwd) == fwd_m .+ 1
        @test CoreRT.n_fourier_moments_bands(lin_model) == lin_m .+ 1
        n_compared += 1
    end
    @info "Phase B parity coverage" n_compared
    @test n_compared > 0
end
