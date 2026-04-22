@testset "IO exports" begin
    @test isdefined(Main, :read_parameters)
    @test read_parameters === vSmartMOM.read_parameters

    params_path = joinpath(@__DIR__, "test_parameters", "PureRayleighParameters.yaml")
    @test read_parameters(params_path) isa vSmartMOM.CoreRT.vSmartMOM_Parameters
end
