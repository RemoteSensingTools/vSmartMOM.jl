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
# Parse failures and forward-build failures are still skipped (RRS/VS
# configs need the explicit RS_type argument, etc.) — those are tested
# elsewhere. But once a config builds forward, the lin path MUST also
# build cleanly: we no longer silently swallow lin-build errors,
# because masking them lets latent regressions slip through CI (see
# the Rayleigh-only LinMode FieldError that surfaced in the v2.1
# Codex review).

using vSmartMOM
using vSmartMOM.CoreRT
using Test

const _PARAMS_DIR = joinpath(@__DIR__, "test_parameters")

_yaml_paths() = filter(p -> endswith(p, ".yaml"), readdir(_PARAMS_DIR; join = true))

@testset "Phase B — forward/lin m_max_bands parity" begin
    n_compared = 0
    for path in _yaml_paths()
        params = try
            parameters_from_yaml(path)
        catch err
            @info "Skipping (parse failed)" path err
            continue
        end

        fwd = try
            model_from_parameters(params)
        catch err
            @info "Skipping (forward build failed)" path = basename(path) err
            continue
        end

        # No longer wrap in try/catch — lin-build failures are real bugs.
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
