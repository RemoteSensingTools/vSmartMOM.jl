module IOExportSmoke

using Test
using vSmartMOM

const PARAMS_PATH = joinpath(@__DIR__, "test_parameters", "PureRayleighParameters.yaml")

@testset "top-level IO exports" begin
    exported_names = (
        :read_parameters,
        :parameters_from_file,
        :parameters_from_source,
        :parameters_from_yaml,
        :parameters_from_dict,
        :read_atmos_profile,
        :read_atmos_profile_dict,
        :GeosChemSource,
        :NetCDFGridSource,
        :NetCDFSource,
        :geoschem_to_dict,
        :read_geoschem_profile,
    )

    for name in exported_names
        @test isdefined(@__MODULE__, name)
        @test getfield(@__MODULE__, name) === getfield(vSmartMOM, name)
    end

    params = read_parameters(PARAMS_PATH)
    @test params isa vSmartMOM.CoreRT.vSmartMOM_Parameters
    @test read_parameters(params) === params
    @test parameters_from_file(PARAMS_PATH).spec_bands == params.spec_bands
    @test parameters_from_yaml(PARAMS_PATH).spec_bands == params.spec_bands
end

end
