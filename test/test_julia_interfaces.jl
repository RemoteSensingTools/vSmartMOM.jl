using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering

@testset "RTModel property introspection" begin
    params = read_parameters(Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[12999.9, 13000.0, 13000.1]"],
            "surface" => ["LambertianSurfaceScalar(0.1)"],
            "nstreams" => 3,
            "polarization_type" => "Stokes_I()",
            "truncation" => "NoTruncation()",
            "depol" => -1,
            "float_type" => "Float64",
            "architecture" => "CPU()"),
        "geometry" => Dict(
            "sza" => 30.0, "vza" => [20.0], "vaz" => [0.0], "obs_alt" => [0]),
        "atmospheric_profile" => Dict(
            "T" => [270.0], "p" => [100.0, 1000.0], "profile_reduction" => -1)))
    model = model_from_parameters(params)
    aliases = (:τ_abs, :τ_rayl, :τ_aer, :aerosol_optics, :greek_rayleigh,
               :greek_cabannes, :ϖ_Cabannes, :obs_geom, :profile, :l_max)
    @test all(hasproperty(model, name) for name in aliases)
    @test all(getproperty(model, name) !== nothing for name in aliases)
end

@testset "GreekCoefs tolerance keywords" begin
    a = Scattering.get_greek_rayleigh(0.03)
    b = deepcopy(a)
    b.β[1] += 1e-8
    @test isapprox(a, b; rtol=0, atol=1e-7)
    @test !isapprox(a, b; rtol=0, atol=1e-9)
end

@testset "elastic Cabannes compatibility overload" begin
    elastic = vSmartMOM.InelasticScattering.noRS{
        Float32}(ϖ_Cabannes=Float32[0.2, 0.3])
    result = vSmartMOM.InelasticScattering.compute_ϖ_Cabannes(
        elastic, 0.03f0, 760f0)
    @test result === elastic.ϖ_Cabannes
    @test result == Float32[1, 1]
end
