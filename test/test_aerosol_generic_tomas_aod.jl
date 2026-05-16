using Test
using NCDatasets
using vSmartMOM
using vSmartMOM.Aerosols
using YAML

const _AOD_RI_DATABASE = joinpath(dirname(@__DIR__), "data",
                                  "refractive_indices_database.yaml")
const _AOD_GCHP_FILE = get(ENV, "VSMARTMOM_GCHP_FILE",
                           "/home/cfranken/data/GEOSChem.Custom.20190702_0000z.nc4")
if !isdefined(Main, :AerosolMieExtinctionLUT)
    include(joinpath(dirname(@__DIR__), "src", "Aerosols", "LUT", "mie_extinction_lut.jl"))
end

@testset "Generic TOMAS schemes" begin
    expected_bins = Dict(:tomas12 => 12, :tomas15 => 15, :tomas30 => 30, :tomas40 => 40)
    expected_spacing = Dict(:tomas12 => :mass_quadruple,
                            :tomas15 => :mass_quadruple,
                            :tomas30 => :mass_double,
                            :tomas40 => :mass_double)

    for variant in keys(expected_bins)
        scheme = TOMASScheme(variant)
        @test scheme isa TOMASScheme
        @test TOMAS15Scheme(variant) isa TOMASScheme
        @test scheme.variant === variant
        @test scheme.n_bins == expected_bins[variant]
        @test nbins(SectionalAerosolData(
            scheme, scheme.size_grid,
            zeros(Float64, scheme.n_bins, 1),
            copy(scheme.species),
            zeros(Float64, scheme.n_bins, 1, length(scheme.species)),
            (; variant),
            nothing, VolumeWeightedMixing(), LinearIntegrationPerBin(8))) == expected_bins[variant]
        @test scheme.size_grid.spacing === expected_spacing[variant]
        @test length(scheme.bin_edges) == scheme.n_bins + 1
        @test length(scheme.bin_centers) == scheme.n_bins
    end
end

@testset "TOMAS legacy shim" begin
    config = YAML.load_file(joinpath(dirname(@__DIR__), "examples",
                                     "aerosol_config_tomas15.yaml"))
    scheme = TOMAS15Scheme(config)
    @test scheme isa TOMASScheme
    @test scheme.variant === :legacy_config
    @test scheme.n_bins == 15
    @test scheme.diam_min == 10.0
    @test scheme.diam_max == 10000.0
end

@testset "TOMAS RI mapping" begin
    db = load_refractive_index_database(_AOD_RI_DATABASE)
    scheme = TOMASScheme(:tomas15)
    comp = SpeciesComposition{Float64}(["SF", "DUST"], [0.25, 0.75], :volume)
    n_eff = effective_ri(comp, db, 0.760, VolumeWeightedMixing(), scheme)
    expected = 0.25 * get_refractive_index(db, "sulfate_suso", 0.760) +
               0.75 * get_refractive_index(db, "dust_opac", 0.760)
    @test n_eff ≈ expected
end

@testset "Bin integration modes" begin
    @test LogNormalFit() isa AbstractBinIntegration
    @test ConstantIntegrationPerBin().nquad == 160
    @test LinearIntegrationPerBin().nquad == 160
    @test ConstantIntegrationPerBin(32).nquad == 32
    @test LinearIntegrationPerBin(32).nquad == 32
    @test DirectBinSum() isa AbstractBinIntegration
    @test_throws ArgumentError ConstantIntegrationPerBin(0)
    @test_throws ArgumentError LinearIntegrationPerBin(0)
end

@testset "AOD diagnostic chunk indexing" begin
    @test vSmartMOM.IO._range_chunks(1:2:9, 2) == [[1, 3], [5, 7], [9]]
    @test vSmartMOM.IO._range_chunks(4:6, 2) == [[4, 5], [6]]
    span, locs = vSmartMOM.IO._chunk_span([1, 3, 5])
    @test span == 1:5
    @test locs == [1, 3, 5]
end

@testset "Column AOD diagnostic" begin
    scheme = TOMASScheme(:tomas15; include_species=["SF"])
    number = zeros(Float64, scheme.n_bins, 2)
    number[4, :] .= 100.0
    mass = zeros(Float64, scheme.n_bins, 2, 1)
    mass[4, :, 1] .= 1.0
    data = SectionalAerosolData(
        scheme, scheme.size_grid, number, ["SF"], mass,
        (; variant=:tomas15),
        load_refractive_index_database(_AOD_RI_DATABASE),
        VolumeWeightedMixing(), LinearIntegrationPerBin(8))

    aod = compute_column_aod(data, [300.0, 700.0, 1000.0],
                             [260.0, 290.0], [0.0, 0.01],
                             [0.760, 1.600, 2.200])
    aod_constant = compute_column_aod(data, [300.0, 700.0, 1000.0],
                                      [260.0, 290.0], [0.0, 0.01],
                                      [0.760, 1.600, 2.200];
                                      integration=ConstantIntegrationPerBin(8))
    aod_direct = compute_column_aod(data, [300.0, 700.0, 1000.0],
                                    [260.0, 290.0], [0.0, 0.01],
                                    [0.760, 1.600, 2.200];
                                    integration=DirectBinSum())
    @test size(aod) == (3,)
    @test all(isfinite.(aod))
    @test all(aod .>= 0)
    @test any(aod .> 0)
    @test all(isfinite.(aod_constant))
    @test all(isfinite.(aod_direct))
end

@testset "Mie extinction LUT" begin
    lut = Main.AerosolMieExtinctionLUT.build_mie_extinction_lut(
        [0.001, 1.0, 100.0], [1.4, 1.5], [0.0, 0.01]; FT=Float32)
    q_mid = Main.AerosolMieExtinctionLUT.evaluate_extinction_lut(lut, 0.1, 1.45, 0.005)
    @test isfinite(q_mid)
    @test q_mid >= 0

    scheme = TOMASScheme(:tomas15; include_species=["SF"])
    number = zeros(Float64, scheme.n_bins, 1)
    number[4, 1] = 100.0
    mass = zeros(Float64, scheme.n_bins, 1, 1)
    mass[4, 1, 1] = 1.0
    data = SectionalAerosolData(
        scheme, scheme.size_grid, number, ["SF"], mass,
        (; variant=:tomas15),
        load_refractive_index_database(_AOD_RI_DATABASE),
        VolumeWeightedMixing(), DirectBinSum())

    aod = compute_column_aod(data, [500.0, 1000.0], [280.0], [0.0], [0.760];
                             mie_lut=lut)
    @test size(aod) == (1,)
    @test isfinite(aod[1])
    @test aod[1] >= 0
end

@testset "GCHP full-file AOD diagnostic" begin
    if get(ENV, "VSMARTMOM_RUN_FULL_GCHP_AOD", "0") == "1" && isfile(_AOD_GCHP_FILE)
        mktempdir() do dir
            out = write_gchp_aod_diagnostic(joinpath(dir, "gchp_aod.nc4"), _AOD_GCHP_FILE)
            NCDataset(out) do ds
                @test haskey(ds, "aod")
                aod = ds["aod"][:, :, :, :]
                @test NCDatasets.dimnames(ds["aod"]) == ("wavelength", "face", "y", "x")
                @test size(aod, 1) == 3
                @test all(isfinite.(aod))
                @test all(aod .>= 0)
                @test maximum(aod) > minimum(aod)
            end
        end
    else
        @info "Skipping full GCHP AOD diagnostic; set VSMARTMOM_RUN_FULL_GCHP_AOD=1 and VSMARTMOM_GCHP_FILE to enable"
        @test true
    end
end
