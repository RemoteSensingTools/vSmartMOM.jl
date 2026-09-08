using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.StandaloneSS

function docs_minimal_parameter_dict()
    return Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[12987.0]"],
            "surface" => ["LambertianSurfaceScalar(0.1)"],
            "quadrature_type" => "GaussLegQuad()",
            "polarization_type" => "Stokes_I()",
            "max_m" => 1,
            "Δ_angle" => 2.0,
            "l_trunc" => 1,
            "depol" => 0.0,
            "float_type" => "Float64",
            "architecture" => "CPU()",
        ),
        "geometry" => Dict(
            "sza" => 30.0,
            "vza" => [0.0],
            "vaz" => [0.0],
            "obs_alt" => [0.0],
        ),
        "atmospheric_profile" => Dict(
            "T" => [260.0, 280.0],
            "p" => [100.0, 600.0, 1000.0],
            "profile_reduction" => -1,
        ),
    )
end

@testset "published source configuration example" begin
    path = joinpath(pkgdir(vSmartMOM), "docs", "src", "pages", "IO", "Schema", "sources.md")
    code = match(r"```julia\n(.*?)```"s, read(path, String)).captures[1]
    example = Module(gensym(:SourceExample))
    Base.include_string(example, code, path)
    R = Base.invokelatest(getfield, example, :R)
    R_sif = Base.invokelatest(getfield, example, :R_sif)
    @test all(isfinite, R)
    @test all(isfinite, R_sif)
    @test all(R_sif .> R)
end

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

@testset "docs IO loaders" begin
    cfg = docs_minimal_parameter_dict()
    params = read_parameters(cfg)

    @test params.architecture isa vSmartMOM.Architectures.CPU
    @test params.brdf[1] isa vSmartMOM.CoreRT.LambertianSurfaceScalar{Float64}
    @test params.spec_bands == [[12987.0]]

    model = model_from_parameters(params)
    R, T = rt_run(model)
    @test size(R) == (1, 1, 1)
    @test size(T) == size(R)
end

@testset "docs Jacobian layout" begin
    layout = CoreRT.ParameterLayout(aerosol_params = 7,
                                    n_aerosols = 2,
                                    n_gases = 3,
                                    n_surface = 1)

    @test CoreRT.n_total(layout) == 19
    @test CoreRT.aerosol_range(layout, 1) == 2:8
    @test CoreRT.aerosol_range(layout, 2) == 9:15
    @test CoreRT.gas_range(layout) == 16:18
    @test CoreRT.surface_index(layout) == 19
end

@testset "docs StandaloneSS vector Jacobian example" begin
    include(joinpath(pkgdir(vSmartMOM), "examples",
                     "standalone_ss_vector_jacobian.jl"))
    @test true
end
