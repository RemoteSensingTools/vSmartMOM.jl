# =================================================================
# Minimal RAMI Smoke Test
#
# Verifies the core RT engine works for a RAMI-style scenario by
# calling rt_run directly with a RAMI YAML config (no gas absorption).
# Bypasses sandbox/rami/rami_tools.jl to avoid world-age issues.
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using Test

yaml = "../sandbox/rami/RamiNoGas.yaml"

@testset "RAMI Smoke Test" begin
    params = parameters_from_yaml(yaml)
    params.architecture = vSmartMOM.Architectures.CPU()
    params.max_m = 3
    params.l_trunc = 20

    @testset "Model construction" begin
        model = model_from_parameters(params)
        @test !isnothing(model)
        @test all(isfinite.(model.profile.T))
        @test all(model.profile.p_full .> 0)
    end

    @testset "RT run" begin
        model = model_from_parameters(params)
        result = rt_run(model)
        R = result isa Tuple ? result[1] : result

        @test ndims(R) >= 2
        @test all(isfinite.(R))

        I_vals = R[:, 1, :]
        @test all(I_vals .>= 0)

        nVza = length(model.obs_geom.vza)
        @test size(R, 1) == nVza
    end

    @testset "Multi-band consistency" begin
        model = model_from_parameters(params)
        Nbands = length(params.spec_bands)
        Nsurf  = length(params.brdf)
        for ib in 1:min(Nbands, Nsurf)
            result = rt_run(model; i_band=ib)
            R = result isa Tuple ? result[1] : result
            @test all(isfinite.(R))
            @test all(R[:, 1, :] .>= 0)
        end
    end
end
