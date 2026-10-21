# The linearized solver handles one band per call: a one-element band vector
# is that band, and several bands are rejected with a clear message.
using Test, vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.InelasticScattering

@testset "Linearized RT band selection" begin
    p = read_parameters(joinpath(@__DIR__, "test_parameters", "JacobianTestRayleigh.yaml"))
    model, lin = model_from_parameters(LinMode(), p)
    ngas = size(lin.τ̇_abs[1], 1)
    naer = CoreRT.n_aerosols(model)
    nsurf = CoreRT.surface_parameter_count(CoreRT.get_surface(model, 1))
    rs = noRS{Float64}()
    by_index = rt_run(rs, model, lin, naer, ngas, nsurf, 1)
    by_vector = rt_run(rs, model, lin, naer, ngas, nsurf, [1])
    @test by_vector[1] == by_index[1]
    @test by_vector[3] == by_index[3]
    @test_throws r"one band per call" rt_run(rs, model, lin, naer, ngas, nsurf, [1, 1])
end
