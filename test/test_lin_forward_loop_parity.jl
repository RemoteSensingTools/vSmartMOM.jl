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
# `test/test_parameters/` holds ONLY CI-buildable configs, so any parse or
# forward/lin build failure here is a HARD failure — not a silent skip. (A
# silent skip on a shared-FS dev box is misleading: it can "see" /home/.../data
# LUTs that CI's clean runners cannot, so configs that need local data must NOT
# live here.) Scenes that need external LUTs / ABSCO data / a GPU / the
# unimplemented H2O-override live under `test/local/test_parameters/` and run
# via `test/local/runtests.jl` on a machine that has the data. If a config
# trips this guard on CI, either fix it or move it to test/local/.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

const _PARAMS_DIR = joinpath(@__DIR__, "test_parameters")

_yaml_paths() = filter(p -> endswith(p, ".yaml"), readdir(_PARAMS_DIR; join = true))

@testset "Phase B — forward/lin m_max_bands parity" begin
    # aarch64 macOS CI runners overflow the default task stack deep inside
    # the HITRAN line-by-line absorption build for the Raman configs
    # (x86-64 survives — the recursion depth sits near the limit). Run each
    # model build on a 64 MiB task stack; behaviorally identical elsewhere.
    # Root-causing the deep recursion in the ARM absorption path is tracked
    # as a follow-up.
    function _build_big_stack(f)
        t = Task(f, 64 << 20)
        schedule(t)
        return fetch(t)
    end

    n_compared = 0
    for path in _yaml_paths()
        cfg = basename(path)
        params = try
            parameters_from_yaml(path)
        catch err
            @error "CI config failed to parse — fix it, or move it to test/local/test_parameters/ if it needs local data/LUTs" config=cfg exception=(err, catch_backtrace())
            @test false
            continue
        end

        fwd = try
            _build_big_stack(() -> model_from_parameters(params))
        catch err
            @error "CI config failed to build forward — fix it, or move it to test/local/test_parameters/ if it needs local data/LUTs" config=cfg exception=(err, catch_backtrace())
            @test false
            continue
        end

        # lin-build failures are real bugs — not wrapped in @test-false.
        lin_model, _ = _build_big_stack(() -> model_from_parameters(LinMode(), params))

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
