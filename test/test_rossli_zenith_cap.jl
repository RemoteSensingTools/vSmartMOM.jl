using Test
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Scattering

# Ross–Li kernel zenith cap (decision 10 of the v2.2 decision log): a
# parameter, default 80°; `nothing` evaluates the kernels at the true angles
# (MODIS/VLIDORT). The VLIDORT comparison is Case D of test/vlidort_baseline.
@testset "Ross–Li zenith cap" begin
    @testset "construction" begin
        s = CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12)       # positional, as before
        @test s isa CoreRT.RossLiSurfaceScalar{Float64}
        @test (s.fvol, s.fgeo, s.fiso, s.zenith_cap_deg) === (0.05, 0.02, 0.12, 80.0)
        @test CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12; zenith_cap_deg = nothing).zenith_cap_deg === nothing
        @test CoreRT.RossLiSurfaceScalar(0.05f0, 0.02f0, 0.12f0; zenith_cap_deg = 85).zenith_cap_deg === 85f0
        @test CoreRT.RossLiSurfaceScalar{Float32}(0.05, 0.02, 0.12).zenith_cap_deg === 80f0
        @test_throws ArgumentError CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12; zenith_cap_deg = 0)
        @test_throws ArgumentError CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12; zenith_cap_deg = 95)
        @test_throws ArgumentError CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12; zenith_cap_deg = true)
        @test sprint(show, s) == "RossLi(fvol=0.05, fgeo=0.02, fiso=0.12)"
        @test occursin("zenith_cap_deg=nothing",
            sprint(show, CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12; zenith_cap_deg = nothing)))
    end

    @testset "kernel angles" begin
        capped = CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12)
        free = CoreRT.RossLiSurfaceScalar(0.05, 0.02, 0.12; zenith_cap_deg = nothing)
        μ80 = cosd(80.0)
        for (μᵢ, μᵣ, Δϕ) in ((0.9, 0.8, 0.4), (0.5, 0.3, 2.0), (μ80 + 0.01, 0.7, 1.0))
            # Below the cap both definitions agree exactly.
            @test CoreRT.reflectance(capped, 1, μᵢ, μᵣ, Δϕ) == CoreRT.reflectance(free, 1, μᵢ, μᵣ, Δϕ)
        end
        μ_grazing = cosd(89.7)     # the most grazing 16-stream Gauss node
        for (μᵣ, Δϕ) in ((0.8, 0.0), (0.5, π / 2))
            # With the cap a grazing direction behaves as θ = 80°…
            @test CoreRT.reflectance(capped, 1, μ_grazing, μᵣ, Δϕ) ==
                  CoreRT.reflectance(capped, 1, cosd(80.0), μᵣ, Δϕ)
            # …without it the LiSparse kernel is evaluated at 89.7°.
            @test CoreRT.reflectance(free, 1, μ_grazing, μᵣ, Δϕ) !=
                  CoreRT.reflectance(capped, 1, μ_grazing, μᵣ, Δϕ)
        end
        # The MODIS-definition BRDF turns negative at grazing (VLIDORT: −2.55
        # at SZA 30°, θ = 89.7°, RAZ 0 for these weights); the capped one does not.
        @test CoreRT.reflectance(free, 1, cosd(30.0), μ_grazing, 0.0) < 0
        @test CoreRT.reflectance(capped, 1, cosd(30.0), μ_grazing, 0.0) > 0
        @test CoreRT.reflectance(free, 2, 0.5, 0.5, 0.0) == 0     # scalar BRDF
    end

    @testset "YAML/TOML surface string" begin
        surf(s) = vSmartMOM.IO.parse_surface_str(s, Float64)
        @test surf("RossLiSurfaceScalar(0.05, 0.02, 0.12)").zenith_cap_deg === 80.0
        for spelling in ("nothing", "null", "None")
            s = surf("RossLiSurfaceScalar(0.05, 0.02, 0.12, zenith_cap_deg=$spelling)")
            @test s.zenith_cap_deg === nothing
            @test (s.fvol, s.fgeo, s.fiso) == (0.05, 0.02, 0.12)
        end
        @test surf("RossLiSurfaceScalar(0.05, 0.02, 0.12, zenith_cap_deg=85)").zenith_cap_deg === 85.0
        kw = surf("RossLiSurfaceScalar(fvol=0.05, fgeo=0.02, fiso=0.12, zenith_cap_deg=null)")
        @test (kw.fvol, kw.fgeo, kw.fiso, kw.zenith_cap_deg) == (0.05, 0.02, 0.12, nothing)
        @test vSmartMOM.IO.parse_surface_str("RossLiSurfaceScalar(0.05, 0.02, 0.12)", Float32) isa
              CoreRT.RossLiSurfaceScalar{Float32}
        @test_throws ArgumentError surf("RossLiSurfaceScalar(0.05, 0.02)")
        @test_throws ArgumentError surf("RossLiSurfaceScalar(0.05, 0.02, 0.12, cap=80)")
        @test_throws ArgumentError surf("RossLiSurfaceScalar(0.05, 0.02, 0.12, zenith_cap_deg=0)")
        @test_throws ArgumentError surf("RossLiSurfaceScalar(0.05, 0.02, 0.12, zenith_cap_deg=true)")
        @test_throws ArgumentError surf("RossLiSurfaceScalar(fvol=0.05, fgeo=0.02)")
        @test_throws ArgumentError surf("RossLiSurfaceScalar(0.05, 0.02, 0.12, fvol=0.05)")
        # Other surfaces parse as before.
        @test surf("CoxMunkSurface(wind_speed=5.0)").n_water === nothing
        @test surf("CoxMunkSurface(wind_speed=5.0, n_water=null)").n_water === nothing
        @test_throws ArgumentError surf("CoxMunkSurface(5.0, shadowing=false)")
    end

    @testset "RT with and without the cap" begin
        params = read_parameters(joinpath(pkgdir(vSmartMOM), "config", "vegetation_rossli.yaml"))
        rossli = only(unique(params.brdf))
        @test rossli.zenith_cap_deg == 80
        R_default = rt_run(model_from_parameters(params)).toa
        params.brdf .= Ref(CoreRT.RossLiSurfaceScalar(rossli.fvol, rossli.fgeo, rossli.fiso;
                                                       zenith_cap_deg = 80))
        @test rt_run(model_from_parameters(params)).toa == R_default
        params.brdf .= Ref(CoreRT.RossLiSurfaceScalar(rossli.fvol, rossli.fgeo, rossli.fiso;
                                                       zenith_cap_deg = nothing))
        R_free = rt_run(model_from_parameters(params)).toa
        @test all(isfinite, R_free)
        @test R_free != R_default
    end
end
