using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using vSmartMOM.Architectures

const VLIDORT_SOLAR_TESTER_NSTREAMS = 8
const VLIDORT_SOLAR_TESTER_NMOMENTS = 2 * VLIDORT_SOLAR_TESTER_NSTREAMS - 1

function solar_tester_stream_counts(config_path::AbstractString)
    params = parameters_from_yaml(config_path)
    FT = params.float_type
    obs_geom = CoreRT.ObsGeometry{FT}(params.sza, params.vza, params.vaz, params.obs_alt)
    quad_points = CoreRT.rt_set_streams(params.quadrature_type,
                                        params.l_trunc,
                                        obs_geom,
                                        params.polarization_type,
                                        Architectures.array_type(params.architecture))
    weights = Array(quad_points.wt_μ)
    return (
        params = params,
        total_nodes = length(weights),
        weighted_nodes = count(abs.(weights) .> 10eps(FT)),
        zero_weight_nodes = count(abs.(weights) .<= 10eps(FT)),
    )
end

@testset "VLIDORT baseline: solar_tester definition checks" begin
    cases = (
        "scalar" => joinpath(CONFIG_DIR, "solar_tester.yaml"),
        "vector" => joinpath(CONFIG_DIR, "solar_tester_vector.yaml"),
    )

    for (label, path) in cases
        diag = solar_tester_stream_counts(path)
        params = diag.params

        @test params.quadrature_type isa CoreRT.GaussQuadHemisphere
        @test diag.weighted_nodes == VLIDORT_SOLAR_TESTER_NSTREAMS
        @test params.l_trunc == VLIDORT_SOLAR_TESTER_NMOMENTS
        @test params.max_m >= VLIDORT_SOLAR_TESTER_NMOMENTS + 1
        @test params.truncation isa Scattering.NoTruncation

        @info "solar_tester discretization [$label]" weighted_nodes=diag.weighted_nodes total_nodes=diag.total_nodes zero_weight_nodes=diag.zero_weight_nodes max_m=params.max_m l_trunc=params.l_trunc truncation=typeof(params.truncation)
    end
end
