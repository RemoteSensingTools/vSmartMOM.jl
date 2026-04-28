using Test
using vSmartMOM
using vSmartMOM.IO
using vSmartMOM.CoreRT

function _minimal_parameter_dict(; surface = ["LambertianSurfaceScalar(0.1)"],
                                   float_type = "Float64",
                                   architecture = "CPU()",
                                   quadrature_type = "GaussQuadHemisphere()")
    return Dict(
        "radiative_transfer" => Dict(
            "spec_bands" => ["[13000.0, 13000.1]"],
            "surface" => surface,
            "quadrature_type" => quadrature_type,
            "polarization_type" => "Stokes_I()",
            "max_m" => 1,
            "Δ_angle" => 0.0,
            "l_trunc" => 20,
            "depol" => 0.0,
            "float_type" => float_type,
            "architecture" => architecture,
        ),
        "geometry" => Dict(
            "sza" => 30.0,
            "vza" => [0.0],
            "vaz" => [0.0],
            "obs_alt" => 1000.0,
        ),
        "atmospheric_profile" => Dict(
            "T" => [290.0],
            "p" => [1000.0, 0.0],
            "profile_reduction" => 1,
        ),
    )
end

@testset "valid minimal dict" begin
    params = parameters_from_dict(_minimal_parameter_dict())

    @test params.float_type === Float64
    @test params.architecture isa vSmartMOM.Architectures.CPU
    @test params.brdf[1] isa LambertianSurfaceScalar{Float64}
    @test params.spec_bands == [[13000.0, 13000.1]]
end

@testset "TOML file input" begin
    mktempdir() do dir
        path = joinpath(dir, "params.toml")
        open(path, "w") do io
            write(io, """
                [radiative_transfer]
                spec_bands = ["[13000.0, 13000.1]"]
                surface = ["LambertianSurfaceScalar(0.1)"]
                quadrature_type = "GaussQuadHemisphere()"
                polarization_type = "Stokes_I()"
                max_m = 1
                "Δ_angle" = 0.0
                l_trunc = 20
                depol = 0.0
                float_type = "Float64"
                architecture = "CPU()"

                [geometry]
                sza = 30.0
                vza = [0.0]
                vaz = [0.0]
                obs_alt = 1000.0

                [atmospheric_profile]
                T = [290.0]
                p = [1000.0, 0.0]
                profile_reduction = 1
            """)
        end

        params = parameters_from_file(path)
        @test params.float_type === Float64
        @test params.brdf[1] isa LambertianSurfaceScalar{Float64}
        @test read_parameters(path).spec_bands == [[13000.0, 13000.1]]
        @test_throws ArgumentError parameters_from_yaml(path)
    end
end

@testset "surface dispatch parser" begin
    surf = vSmartMOM.IO.parse_surface_str("CoxMunkSurface(wind_speed=5.0)", Float32)

    @test surf isa CoxMunkSurface{Float32}
    @test surf.wind_speed === Float32(5)
    @test_throws ArgumentError vSmartMOM.IO.parse_surface_str("UnknownSurface(1.0)", Float64)
    @test_throws ArgumentError vSmartMOM.IO.parse_surface_str("LambertianSurfaceScalar()", Float64)
    @test_throws ArgumentError vSmartMOM.IO.parse_surface_str("CoxMunkSurface()", Float64)
end

@testset "environment path expansion" begin
    env_key = "VSMARTMOM_TEST_LUT_DIR"
    old_value = get(ENV, env_key, nothing)

    try
        ENV[env_key] = "/tmp/vsmartmom-luts"
        @test vSmartMOM.IO._expand_env_path("\${ENV:$(env_key)}/O2.jld2") == normpath("/tmp/vsmartmom-luts/O2.jld2")
        @test vSmartMOM.IO._expand_env_path("relative/O2.jld2") == "relative/O2.jld2"

        delete!(ENV, env_key)
        @test_throws ArgumentError vSmartMOM.IO._expand_env_path("\${ENV:$(env_key)}/O2.jld2")
    finally
        if old_value === nothing
            delete!(ENV, env_key)
        else
            ENV[env_key] = old_value
        end
    end
end

@testset "invalid configuration raises ArgumentError" begin
    cfg = _minimal_parameter_dict(float_type = "Float16")
    @test_throws ArgumentError parameters_from_dict(cfg)

    cfg = _minimal_parameter_dict(architecture = "TPU()")
    @test_throws ArgumentError parameters_from_dict(cfg)

    cfg = _minimal_parameter_dict(quadrature_type = "BadQuad()")
    @test_throws ArgumentError parameters_from_dict(cfg)

    cfg = _minimal_parameter_dict()
    delete!(cfg["geometry"], "sza")
    @test_throws ArgumentError parameters_from_dict(cfg)
end
