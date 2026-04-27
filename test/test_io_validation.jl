using Test
using vSmartMOM
using vSmartMOM.IO
using vSmartMOM.IO.Formats
using vSmartMOM.IO.Sources

struct UnsupportedSourceForTest <: Formats.IOSource end

@testset "format registry dispatch" begin
    data = Dict("radiative_transfer" => Dict())

    @test Formats.load_config(Formats.DictSource(data)) === data
    @test_throws ArgumentError Formats.load_config(Formats.FileSource("scene.txt"))
    @test_throws ArgumentError Formats.load_config(UnsupportedSourceForTest())
end

@testset "NetCDF source constructors" begin
    mktemp() do path, io
        close(io)

        @test Sources.GeosChemSource(path, 1, 2, 3) isa Sources.GeosChemSource
        @test Sources.NetCDFGridSource(path, 1, 2) isa Sources.NetCDFGridSource

        @test_throws ArgumentError Sources.GeosChemSource(path, 0, 2, 3)
        @test_throws ArgumentError Sources.GeosChemSource(path, 1, 0, 3)
        @test_throws ArgumentError Sources.GeosChemSource(path, 1, 2, 7)
        @test_throws ArgumentError Sources.NetCDFGridSource(path, 0, 2)
        @test_throws ArgumentError Sources.NetCDFGridSource(path, 1, 0)
    end

    @test_throws ArgumentError Sources.GeosChemSource("missing-geoschem.nc4", 1, 1, 1)
    @test_throws ArgumentError Sources.NetCDFGridSource("missing-grid.nc", 1, 1)
end

@testset "atmospheric profile validation" begin
    @test_throws ArgumentError read_atmos_profile("profile.txt")
end
