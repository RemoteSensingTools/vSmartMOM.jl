using Test
using vSmartMOM

@testset "docs quickstart" begin
    scene = joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml")
    @test isfile(scene)

    params = read_parameters(scene)
    @test params.architecture isa vSmartMOM.Architectures.CPU

    model = model_from_parameters(params)
    R, T = rt_run(model)

    @test size(R) == (length(params.vza), params.polarization_type.n, length(params.spec_bands[1]))
    @test size(T) == size(R)
    @test all(isfinite, R)
    @test all(isfinite, T)
end
