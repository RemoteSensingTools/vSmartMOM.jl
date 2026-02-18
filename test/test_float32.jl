using vSmartMOM
using Test

@testset "Float32 Consistency" begin

    yaml64 = "test_parameters/PureRayleighParameters.yaml"
    yaml32 = "test_parameters/PureRayleighParameters_f32.yaml"

    @testset "Float64 baseline" begin
        params64 = vSmartMOM.parameters_from_yaml(yaml64)
        params64.max_m = 2
        params64.l_trunc = 10
        params64.architecture = vSmartMOM.CPU()
        model64 = model_from_parameters(params64)
        R64, _, _, _, _, _, _ = rt_run(model64)
        @test eltype(R64) == Float64
        @test all(isfinite, R64)
        @test size(R64, 1) == length(params64.vza)
    end

    @testset "Float32 forward" begin
        params32 = vSmartMOM.parameters_from_yaml(yaml32)
        params32.max_m = 2
        params32.l_trunc = 10
        params32.architecture = vSmartMOM.CPU()
        model32 = model_from_parameters(params32)
        R32, _, _, _, _, _, _ = rt_run(model32)
        @test eltype(R32) == Float32
        @test all(isfinite, R32)
        @test size(R32, 1) == length(params32.vza)
    end

    @testset "Float32 vs Float64 agreement" begin
        params64 = vSmartMOM.parameters_from_yaml(yaml64)
        params64.max_m = 2
        params64.l_trunc = 10
        params64.architecture = vSmartMOM.CPU()
        model64 = model_from_parameters(params64)
        R64, _, _, _, _, _, _ = rt_run(model64)
        
        params32 = vSmartMOM.parameters_from_yaml(yaml32)
        params32.max_m = 2
        params32.l_trunc = 10
        params32.architecture = vSmartMOM.CPU()
        model32 = model_from_parameters(params32)
        R32, _, _, _, _, _, _ = rt_run(model32)
        
        R64_f32 = Float32.(R64)
        rel_diff = abs.(R64_f32 .- R32) ./ max.(abs.(R64_f32), Float32(1e-10))
        max_rel = maximum(rel_diff)
        mean_rel = sum(rel_diff) / length(rel_diff)
        
        println("  Float32 vs Float64: max_rel=$max_rel, mean_rel=$mean_rel")
        @test max_rel < 0.01  # 1% relative tolerance
    end
end
