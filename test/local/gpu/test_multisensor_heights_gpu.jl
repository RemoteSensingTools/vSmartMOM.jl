using Test
using vSmartMOM
using vSmartMOM.CoreRT

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
    end
end
