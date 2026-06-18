# Setup-dependent test suite — NOT part of the CI unit suite (test/runtests.jl).
#
# These configs need something the CI environment does not have:
#   - H2O-override configs (2Band/3Band/3Band_canopy/S2Band/ParamsEMIT_MODTRANcomp*):
#     they declare H2O in the molecule list + an explicit H2O VMR — a config-level
#     override that is NOT yet implemented (H2O is currently driven by
#     atmospheric_profile.q). Parked here until the override lands.
#   - *_newLUT* configs: need the ABSCO LUTs via the VSMARTMOM_HITRAN_LUT_DIR
#     environment variable (data not shipped in the repo).
#   - ThreeBandsParameters / O2ParametersVS / O2_parameters2_SIF_grid /
#     WCO2_parameters_SIF_grid: legacy example scenes kept for reference.
#
# Run locally with:
#     julia --project=test test/setup/runtests.jl
# Any config that *does* build in your environment is checked for forward/lin
# Fourier-loop parity; the rest report what they need and are skipped.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

const _SETUP_PARAMS = joinpath(@__DIR__, "test_parameters")
const _SETUP_ORIG_CWD = pwd()
# Run from the test/ root so any relative data paths inside the configs resolve.
cd(normpath(joinpath(@__DIR__, "..")))
try
    @testset "setup-dependent configs (forward/lin parity where buildable)" begin
        n_ok = 0
        for path in filter(e -> endswith(e, ".yaml"), readdir(_SETUP_PARAMS; join = true))
            cfg = basename(path)
            params = try
                parameters_from_yaml(path)
            catch err
                @info "needs setup (parse)" config = cfg reason = sprint(showerror, err)
                continue
            end
            fwd = try
                model_from_parameters(params)
            catch err
                @info "needs setup (build)" config = cfg reason = sprint(showerror, err)
                continue
            end
            lin_model, _ = model_from_parameters(LinMode(), params)
            @test CoreRT.m_max_bands(fwd) == CoreRT.m_max_bands(lin_model)
            n_ok += 1
        end
        @info "setup configs buildable in this environment" n_ok
    end
finally
    cd(_SETUP_ORIG_CWD)
end
