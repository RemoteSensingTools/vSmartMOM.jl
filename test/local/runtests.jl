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
# Any config that *does* build in your environment is checked for forward/lin
# Fourier-loop parity; the rest report what they need and are skipped.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

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
            for path in filter(e -> endswith(e, ".yaml"), readdir(params_dir; join = true))
                cfg = basename(path)
                params = try
                    parameters_from_yaml(path)
                catch err
                    @info "needs local data/setup (parse)" config = cfg reason = sprint(showerror, err)
                    continue
                end
                fwd = try
                    model_from_parameters(params)
                catch err
                    @info "needs local data/setup (build)" config = cfg reason = sprint(showerror, err)
                    continue
                end
                lin_model, _ = model_from_parameters(LinMode(), params)
                @test CoreRT.m_max_bands(fwd) == CoreRT.m_max_bands(lin_model)
                n_ok += 1
            end
            @info "local configs buildable in this environment" n_ok
        end
    finally
        cd(orig_cwd)
    end
end
