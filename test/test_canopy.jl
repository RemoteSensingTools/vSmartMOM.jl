# =================================================================
# Canopy Surface Tests
#
# Tests the CanopySurface integration into the standard rt_run() flow:
#   1. CanopySurface struct construction
#   2. YAML-based construction with canopy: section
#   3. Forward RT with CanopySurface vs Lambertian baseline
#   4. Multi-layer canopy consistency
#   5. Comparison with legacy rt_run_canopy() (when available)
# =================================================================

using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.Scattering
using vSmartMOM.InelasticScattering
using CanopyOptics
using Test
using LinearAlgebra

println("="^60)
println("Canopy Surface Tests")
println("="^60)

@testset "CanopySurface type construction" begin
    soil = LambertianSurfaceScalar(0.1)
    canopy = CanopySurface(; soil=soil, LAI=3.0)
    @test canopy.LAI == 3.0
    @test canopy.n_layers == 1
    @test canopy.soil === soil
    @test canopy.include_atm == false
    @test canopy._cache === nothing

    canopy_ml = CanopySurface(; soil=soil, LAI=4.0, n_layers=4,
                               leaf_reflectance=0.45, leaf_transmittance=0.1)
    @test canopy_ml.n_layers == 4
    @test canopy_ml.leaf_reflectance == 0.45
    @test canopy_ml.leaf_transmittance == 0.1

    canopy_frac = CanopySurface(; soil=soil, LAI=3.0, n_layers=3,
                                 lai_fractions=[0.5, 0.3, 0.2])
    @test canopy_frac.lai_fractions == [0.5, 0.3, 0.2]
end

@testset "CanopySurface YAML parsing" begin
    params = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    @test params.brdf[1] isa CanopySurface
    canopy = params.brdf[1]
    @test canopy.LAI == 3.0
    @test canopy.n_layers == 1
    @test canopy.soil isa LambertianSurfaceScalar
    @test canopy.soil.albedo ≈ 0.1
    @test canopy.leaf_reflectance ≈ 0.4
    @test canopy.leaf_transmittance ≈ 0.05
end

@testset "CanopySurface forward RT" begin
    params = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()

    model = model_from_parameters(params)
    @test !isnothing(model)
    @test model.params.brdf[1] isa CanopySurface

    result = rt_run(model; i_band=1)
    R = result isa Tuple ? result[1] : result

    nVza = length(model.obs_geom.vza)
    nStokes = model.params.polarization_type.n
    nSpec = length(model.params.spec_bands[1])

    @test ndims(R) == 3
    @test size(R, 1) == nVza
    @test size(R, 2) == nStokes
    @test size(R, 3) == nSpec
    @test all(isfinite.(R))

    I_vals = R[:, 1, :]
    @test all(I_vals .> 0)
    @test maximum(I_vals) < 1.0

    println("  Canopy RT I: min=$(minimum(I_vals)), max=$(maximum(I_vals))")
end

@testset "Canopy vs Lambertian comparison" begin
    params_canopy = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    params_canopy.architecture = vSmartMOM.Architectures.CPU()

    params_lamb = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    params_lamb.architecture = vSmartMOM.Architectures.CPU()
    params_lamb.brdf[1] = LambertianSurfaceScalar(0.1)

    model_canopy = model_from_parameters(params_canopy)
    model_lamb   = model_from_parameters(params_lamb)

    R_canopy = rt_run(model_canopy; i_band=1)[1]
    invalidate_canopy_cache!(model_canopy.params.brdf[1])
    R_lamb   = rt_run(model_lamb; i_band=1)[1]

    @test all(isfinite.(R_canopy))
    @test all(isfinite.(R_lamb))

    @test !isapprox(R_canopy, R_lamb, rtol=0.01)
    println("  Canopy I[1]: $(R_canopy[1,1,1])")
    println("  Lambertian I[1]: $(R_lamb[1,1,1])")
end

@testset "Multi-layer canopy consistency" begin
    params = parameters_from_yaml("test_parameters/CanopyTest.yaml")
    params.architecture = vSmartMOM.Architectures.CPU()

    soil = params.brdf[1].soil
    LAI = params.brdf[1].LAI

    params.brdf[1] = CanopySurface(; soil=soil, LAI=LAI, n_layers=1)
    model1 = model_from_parameters(params)
    R1 = rt_run(model1; i_band=1)[1]
    invalidate_canopy_cache!(model1.params.brdf[1])

    params.brdf[1] = CanopySurface(; soil=soil, LAI=LAI, n_layers=4)
    model4 = model_from_parameters(params)
    R4 = rt_run(model4; i_band=1)[1]
    invalidate_canopy_cache!(model4.params.brdf[1])

    @test all(isfinite.(R1))
    @test all(isfinite.(R4))

    println("  1-layer I[1]: $(R1[1,1,1])")
    println("  4-layer I[1]: $(R4[1,1,1])")
end

@testset "ParameterLayout canopy extension" begin
    layout = ParameterLayout(n_aerosols=1, n_gases=2, n_surface=1, n_canopy=3)
    @test n_total(layout) == 7 + 2 + 1 + 3
    @test canopy_range(layout) == 11:13
    @test surface_range(layout) == 10:10

    layout0 = ParameterLayout(n_aerosols=0, n_gases=0, n_surface=1, n_canopy=0)
    @test n_total(layout0) == 1
    @test canopy_range(layout0) == 2:1  # empty range
end

println("\n" * "="^60)
println("Canopy surface tests complete.")
println("="^60)
