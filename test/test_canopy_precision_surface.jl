module CanopySurfacePrecisionTests

using Test, TOML, LinearAlgebra, vSmartMOM, CanopyOptics
using vSmartMOM.CoreRT, vSmartMOM.InelasticScattering

const SETTINGS = TOML.parsefile(joinpath(@__DIR__, "canopy_fixtures/precision_surface.toml"))
const ORIGINAL_BLAS_THREADS = BLAS.get_num_threads()
BLAS.set_num_threads(SETTINGS["numerics"]["blas_threads"])

function conservative_reflection(::Type{FT}, streams, nlayer, architecture) where {FT}
    c, n = SETTINGS["canopy"], SETTINGS["numerics"]
    arr = vSmartMOM.Architectures.array_type(architecture)
    pol = CoreRT.Stokes_I{FT}()
    obs = CoreRT.ObsGeometry(FT(c["solar_zenith_degrees"]), FT[0], FT[0], FT[0])
    qp = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(), 2streams-1, obs, pol, arr)
    dim = length(qp.qp_μ)
    added = CoreRT.make_added_layer(noRS{FT}(), FT, arr, (dim,dim), 1)
    canopy = CoreRT.CanopySurface(;
        soil=CoreRT.LambertianSurfaceScalar(FT(c["soil_albedo"])),
        LAI=FT(c["lai"]), n_layers=nlayer,
        LAD=beta_leaves(FT(c["leaf_beta_alpha"]), FT(c["leaf_beta_beta"])),
        leaf_reflectance=FT(c["reflectance"]), leaf_transmittance=FT(c["transmittance"]),
        canopy_quadrature=CanopyQuadrature(n_leaf=n["leaf_nodes"], n_azimuth=n["azimuth_nodes"]))
    CoreRT.create_surface_layer!(canopy, added, true, 0, pol, qp, arr(zeros(FT,1)), architecture)
    mu, w = collect(qp.qp_μ), Array(qp.wt_μ)
    direct = sum(mu .* w .* Array(added.j₀⁻)[:,1,1]) / qp.μ₀
    diffuse = FT(2)*sum(mu .* w .* (Array(added.r⁻⁺)[:,:,1]*ones(FT,dim)))
    (; direct, diffuse)
end

@testset "Canopy numerical controls" begin
    n = SETTINGS["numerics"]
    for FT in (Float32, Float64)
        controls = CoreRT.CanopyNumericalParameters(FT;
            dτ_max_threshold=n["custom_threshold"], dτ_min_floor=n["custom_floor"])
        @test controls.dτ_max_threshold isa FT
        @test controls.dτ_min_floor isa FT
        soil = CoreRT.LambertianSurfaceScalar(FT(SETTINGS["canopy"]["soil_albedo"]))
        surfaces = CoreRT.AbstractSurfaceType[soil]
        dict = Dict("LAI" => SETTINGS["canopy"]["lai"], "numerics" => Dict(
            "dτ_max_threshold" => n["custom_threshold"], "dτ_min_floor" => n["custom_floor"]))
        vSmartMOM.IO._parse_canopy_section(dict, FT, surfaces)
        @test surfaces[1].canopy_numerics === controls
        perturbed = CoreRT._perturb_soil_albedo(surfaces[1], zero(FT))
        @test perturbed.canopy_numerics === controls
        @test_throws ArgumentError CoreRT.CanopyNumericalParameters(FT; dτ_max_threshold=zero(FT))
        @test_throws ArgumentError CoreRT.CanopyNumericalParameters(FT; dτ_min_floor=-one(FT))
        @test_throws ArgumentError CoreRT.CanopyNumericalParameters(FT; dτ_max_threshold=FT(NaN))
    end
end

@testset "Production canopy conservative photon budget" begin
    n, v = SETTINGS["numerics"], SETTINGS["validation"]
    for layers in n["physical_layers"]
        previous_error = Inf
        for streams in n["streams"]
            reference = conservative_reflection(Float64, streams, layers, vSmartMOM.Architectures.CPU())
            single = conservative_reflection(Float32, streams, layers, vSmartMOM.Architectures.CPU())
            for key in (:direct, :diffuse)
                if streams in n["converged_streams"]
                    @test abs(getproperty(reference,key)-1) <= v["conservation_atol"]
                    @test abs(getproperty(single,key)-1) <= v["conservation_atol"]
                end
                @test abs(getproperty(single,key)-getproperty(reference,key)) <= v["precision_atol"]
            end
            angular_error = max(abs(reference.direct-1), abs(reference.diffuse-1))
            @test angular_error < previous_error
            previous_error = angular_error
            println("Conservative canopy: ", (;streams,layers,Float32=single,Float64=reference))
        end
    end
end

if get(ENV,"VSMART_CANOPY_INCREMENT_GPU","false")=="true"
    using CUDA
    CUDA.functional() || error("CUDA is required for the requested canopy test")
    CUDA.allowscalar(false)
    @testset "Production canopy photon budget: CUDA" begin
        n, v = SETTINGS["numerics"], SETTINGS["validation"]
        for FT in (Float32,Float64), layers in n["physical_layers"]
            result = conservative_reflection(FT, last(n["streams"]), layers, vSmartMOM.Architectures.GPU())
            @test abs(result.direct-1) <= v["conservation_atol"]
            @test abs(result.diffuse-1) <= v["conservation_atol"]
        end
    end
end

BLAS.set_num_threads(ORIGINAL_BLAS_THREADS)
end
