using Test
using vSmartMOM
using vSmartMOM.CoreRT

function _run_height_lin(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = isnothing(params.scattering_params) ? 0 :
           length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    return rt_run(model, lin_model, NAer, NGas, 1)
end

@testset "interior-height multisensor GPU" begin
    config = joinpath("..", "config", "quickstart.yaml")

    for FT in (Float64, Float32)
        cpu_params = parameters_from_yaml(config)
        cpu_params.float_type = FT
        cpu_params.obs_alt = [0.0, 5.0]
        cpu_params.vza = [60.0, 45.0]
        cpu_params.vaz = [180.0, 180.0]
        cpu_params.architecture = vSmartMOM.Architectures.CPU()
        cpu_result = rt_run(model_from_parameters(cpu_params))

        gpu_params = parameters_from_yaml(config)
        gpu_params.float_type = FT
        gpu_params.obs_alt = [0.0, 5.0]
        gpu_params.vza = [60.0, 45.0]
        gpu_params.vaz = [180.0, 180.0]
        gpu_params.architecture = vSmartMOM.Architectures.GPU()
        gpu_result = rt_run(model_from_parameters(gpu_params))

        rtol = FT === Float64 ? 1e-10 : 2e-5
        @test gpu_result.toa ≈ cpu_result.toa rtol=rtol
        @test gpu_result.boa ≈ cpu_result.boa rtol=rtol
        @test length(gpu_result.levels) == 1
        @test gpu_result.levels[1].height_km == FT(5)
        @test gpu_result.levels[1].upwelling ≈
              cpu_result.levels[1].upwelling rtol=rtol
        @test gpu_result.levels[1].downwelling ≈
              cpu_result.levels[1].downwelling rtol=rtol
        @test gpu_result.levels[1].unscattered_downwelling ≈
              cpu_result.levels[1].unscattered_downwelling rtol=rtol
        @test gpu_result.levels[1].unscattered_downwelling[1, 1, 1] > 0
        @test iszero(gpu_result.levels[1].unscattered_downwelling[2, 1, 1])

        cpu_lin_params = parameters_from_yaml(config)
        cpu_lin_params.float_type = FT
        cpu_lin_params.obs_alt = [0.0, 5.0]
        cpu_lin_params.vza = [60.0, 45.0]
        cpu_lin_params.vaz = [180.0, 180.0]
        cpu_lin_params.architecture = vSmartMOM.Architectures.CPU()
        cpu_lin = _run_height_lin(cpu_lin_params)

        gpu_lin_params = parameters_from_yaml(config)
        gpu_lin_params.float_type = FT
        gpu_lin_params.obs_alt = [0.0, 5.0]
        gpu_lin_params.vza = [60.0, 45.0]
        gpu_lin_params.vaz = [180.0, 180.0]
        gpu_lin_params.architecture = vSmartMOM.Architectures.GPU()
        gpu_lin = _run_height_lin(gpu_lin_params)

        lin_rtol = FT === Float64 ? 2e-9 : 5e-4
        @test gpu_lin.toa ≈ cpu_lin.toa rtol=lin_rtol
        @test gpu_lin.boa ≈ cpu_lin.boa rtol=lin_rtol
        @test gpu_lin.toa_jacobian ≈ cpu_lin.toa_jacobian rtol=lin_rtol
        @test gpu_lin.boa_jacobian ≈ cpu_lin.boa_jacobian rtol=lin_rtol
        gpu_level = only(gpu_lin.levels)
        cpu_level = only(cpu_lin.levels)
        @test gpu_level.upwelling ≈ cpu_level.upwelling rtol=lin_rtol
        @test gpu_level.downwelling ≈ cpu_level.downwelling rtol=lin_rtol
        @test gpu_level.unscattered_downwelling ≈
              cpu_level.unscattered_downwelling rtol=lin_rtol
        @test gpu_level.upwelling_jacobian ≈
              cpu_level.upwelling_jacobian rtol=lin_rtol
        @test gpu_level.downwelling_jacobian ≈
              cpu_level.downwelling_jacobian rtol=lin_rtol
        @test gpu_level.unscattered_downwelling_jacobian ≈
              cpu_level.unscattered_downwelling_jacobian rtol=lin_rtol
    end
end
