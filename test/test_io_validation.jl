using Test
using vSmartMOM
using vSmartMOM.IO
using vSmartMOM.IO.Formats
using vSmartMOM.IO.Sources
using NCDatasets

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

    @test vSmartMOM.IO._geoschem_specific_humidity([0.0, 1.0, 10.0]) ==
          [0.0, 0.001, 0.01]
end

@testset "GEOS-Chem conversion contract" begin
    mktempdir() do dir
        path = joinpath(dir, "geoschem.nc4")
        ds = NCDataset(path, "c")
        try
            for (name, n) in (("Xdim", 1), ("Ydim", 1), ("nf", 1),
                              ("lev", 2), ("time", 1))
                defDim(ds, name, n)
            end

            defVar(ds, "lats", Float64, ("Xdim", "Ydim", "nf"))[:, :, :] .= 34.0
            defVar(ds, "lons", Float64, ("Xdim", "Ydim", "nf"))[:, :, :] .= -118.0
            defVar(ds, "time", Float64, ("time",))[:] .= 0.0
            defVar(ds, "Met_DELP", Float64,
                   ("Xdim", "Ydim", "nf", "lev", "time"))[:, :, :, :, :] .=
                reshape([400.0, 500.0], 1, 1, 1, 2, 1)
            defVar(ds, "Met_PS2WET", Float64,
                   ("Xdim", "Ydim", "nf", "time"))[:, :, :, :] .= 1000.0

            T_var = defVar(ds, "Met_T", Float64,
                           ("Xdim", "Ydim", "nf", "lev", "time"))
            T_var.attrib["units"] = "K"
            T_var[:, :, :, :, :] .= reshape([280.0, 240.0], 1, 1, 1, 2, 1)

            q_var = defVar(ds, "Met_SPHU", Float64,
                           ("Xdim", "Ydim", "nf", "lev", "time"))
            q_var.attrib["units"] = "g kg-1"
            q_var[:, :, :, :, :] .= reshape([10.0, 1.0], 1, 1, 1, 2, 1)

            co2_var = defVar(ds, "SpeciesConcVV_CO2", Float64,
                             ("Xdim", "Ydim", "nf", "lev", "time"))
            co2_var.attrib["units"] = "mol mol-1 dry"
            co2_var[:, :, :, :, :] .= reshape([4.2e-4, 4.0e-4], 1, 1, 1, 2, 1)
        finally
            close(ds)
        end

        src = Sources.GeosChemSource(path, 1, 1, 1)
        config = vSmartMOM.IO.geoschem_to_dict(src)
        @test haskey(config, "radiative_transfer")
        @test haskey(config, "geometry")
        @test config["atmospheric_profile"]["T"] == [240.0, 280.0]
        @test config["atmospheric_profile"]["p"] == [100.0, 600.0, 1000.0]
        @test config["atmospheric_profile"]["q"] == [0.001, 0.01]
        @test config["absorption"]["fixed_molecules"] == [["O2", "CO2"]]
        @test config["absorption"]["vmr"]["O2"] == 0.2095
        @test !haskey(config["absorption"]["vmr"], "H2O")

        params = parameters_from_source(src)
        @test params.T == [240.0, 280.0]
        @test params.q == [0.001, 0.01]
        @test params.absorption_params.vmr["CO2"] == [4.0e-4, 4.2e-4]
    end
end

@testset "atmospheric profile validation" begin
    @test_throws ArgumentError read_atmos_profile("profile.txt")

    profile = read_atmos_profile_dict(Dict(
        "T" => [240, 280],
        "p_half" => [100, 500, 1000],
        "q" => [0.0, 0.01],
        "vmr" => Dict("CO2" => [4.0e-4, 4.1e-4], "O2" => 0.2095),
    ))
    @test profile.T == [240.0, 280.0]
    @test profile.q == [0.0, 0.01]
    @test profile.p_full == [300.0, 750.0]
    @test profile.p_half == [100.0, 500.0, 1000.0]
    @test profile.vmr["CO2"] == [4.0e-4, 4.1e-4]
    @test profile.vmr["O2"] == 0.2095
    @test length(profile.Δz) == 2
    @test all(profile.Δz .> 0)

    interface_profile = read_atmos_profile_dict(Dict(
        "T" => [240, 280],
        "p_half" => [100, 500, 1000],
        "q" => [0.0, 0.01, 0.02],
        "vmr" => Dict("CO2" => [4.0e-4, 4.1e-4, 4.3e-4]),
    ))
    @test interface_profile.q == vSmartMOM.CoreRT._layer_centered_input(
        "specific humidity", [0.0, 0.01, 0.02], [100.0, 500.0, 1000.0])
    @test interface_profile.vmr["CO2"] ==
          vSmartMOM.CoreRT._layer_centered_input(
              "VMR CO2", [4.0e-4, 4.1e-4, 4.3e-4],
              [100.0, 500.0, 1000.0])

    dry_profile = read_atmos_profile_dict(Dict(
        "T" => [250.0], "p_half" => [100.0, 1000.0]))
    @test dry_profile.q == [0.0]
    @test isempty(dry_profile.vmr)

    @test_throws ArgumentError read_atmos_profile_dict(Dict(
        "T" => [250.0], "p_half" => [1000.0, 100.0]))
    @test_throws ArgumentError read_atmos_profile_dict(Dict(
        "T" => [250.0], "p_half" => [100.0, 1000.0], "q" => [1.0]))
end
