using Test, vSmartMOM, YAML, Logging
using vSmartMOM.CoreRT

@testset "Shared parsed HITRAN data" begin
    clear_spectroscopy_cache!()
    tasks = [Threads.@spawn CoreRT._hitran_lines("O2",Float64) for _ in 1:4]
    cached = fetch.(tasks)
    @test all(x -> x === first(cached), cached)
    @test length(CoreRT._HITRAN_LINES) == 1
    fresh = CoreRT.AtmosphericAbsorption.load_lines(
        CoreRT.AtmosphericAbsorption.HitranPort(artifact("O2")); FT=Float64)
    @test first(cached).ν0 == fresh.ν0
    @test first(cached).S == fresh.S
    @test CoreRT._hitran_lines("O2",Float32) !== first(cached)
    @test length(CoreRT._HITRAN_LINES) == 2
    with_logger(NullLogger()) do
        cfg = YAML.load_file("test_parameters/JacobianTestFast.yaml")
        delete!(cfg,"scattering")
        params = read_parameters(cfg)
        params.spec_bands[1] = collect(range(first(params.spec_bands[1]),last(params.spec_bands[1]);length=5))
        forward = model_from_parameters(params)
        linearized,_ = model_from_parameters(LinMode(),params)
        @test forward.τ_abs[1] ≈ linearized.τ_abs[1] rtol=1e-12
        @test CoreRT._hitran_lines("O2",Float64) === first(cached)
        @test first(cached).S == fresh.S
        context = BatchContext(params)
        @test only(only(context.absorption_models)).lines === first(cached)
        @test context.model.τ_abs[1] ≈ forward.τ_abs[1] rtol=1e-12
    end
    clear_spectroscopy_cache!()
    @test isempty(CoreRT._HITRAN_LINES)
    @test first(cached).S == fresh.S # clearing ownership does not mutate users
end
