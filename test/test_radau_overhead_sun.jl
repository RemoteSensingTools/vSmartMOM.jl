# =========================================================================
# Radau quadrature with the sun at zenith (SZA = 0, μ₀ = 1).
#
# μ₀ = 1 coincides with the fixed Radau node, so `rt_set_streams` takes the
# single-interval branch. That branch must map the weights from [−1, 1] to
# [0, 1] (factor ½) exactly like the two-interval branch, and the solar test
# must compare μ₀ with the MAPPED nodes. Before the fix the weights summed to
# 2 and the radiances were negative. The SZA = 0 solution must be continuous
# with a slightly oblique sun and agree with Gauss quadrature. A slightly
# oblique sun can make the Radau sub-interval [μ₀, 1] narrower than the node
# spacing representable in FT; coincident nodes must then be merged (they
# used to make `rt_set_streams` throw), and nodes a few ulps apart must not
# cost the elemental kernels accuracy.
# =========================================================================

using Test
using Logging
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering: Stokes_IQU

const _RADAU_YAML = joinpath(@__DIR__, "test_parameters", "PureRayleighParameters.yaml")

@testset "Radau quadrature at SZA = 0" begin
    pol_type = Stokes_IQU()
    arr = Array

    @testset "weights and nodes ($FT)" for FT in (Float64, Float32)
        rtol = FT == Float64 ? 1e-12 : 1e-6
        for Ltrunc in (5, 10, 21, 127)
            # Small SZA makes the sub-interval [μ₀, 1] narrower than the FT
            # node spacing (in Float32 up to SZA ≈ 1°): coincident nodes are
            # merged, μ₀ itself is kept.
            for sza in (0.0, 1e-6, 1e-4, 0.1, 1.0, 30.0)
                obs = CoreRT.ObsGeometry(FT(sza), FT[0, 30], FT[0, 0], FT(0))
                qp = CoreRT.rt_set_streams(CoreRT.RadauQuad(), Ltrunc, obs, pol_type, arr)
                wt = Array(qp.wt_μ)
                μ = Array(qp.qp_μ)
                @test qp.μ₀ == cosd(FT(sza))
                @test allunique(μ)
                @test all(0 .< μ .<= 1)
                # Half-space weights integrate dμ over [0, 1] ...
                @test sum(wt) ≈ 1 rtol=rtol
                # ... and ∫₀¹ μ dμ = 1/2 exactly.
                @test sum(wt .* μ) ≈ 0.5 rtol=rtol
                # μ₀ is a weighted ordinate.
                @test μ[qp.iμ₀] == qp.μ₀
                @test wt[qp.iμ₀] > 0
            end
        end
    end

    @testset "radiance continuity at SZA = 0" begin
        function run(quad, sza)
            with_logger(ConsoleLogger(stderr, Logging.Warn)) do
                p = parameters_from_yaml(_RADAU_YAML)
                p.architecture = vSmartMOM.Architectures.CPU()
                p.polarization_type = pol_type
                p.quadrature_type = quad
                p.sza = sza
                p.vza = [0.0, 20.0, 45.0, 70.0]
                p.vaz = [0.0, 60.0, 120.0, 180.0]
                p.brdf = [LambertianSurfaceScalar(0.2)]
                model = model_from_parameters(p; external_solar=false)
                model.τ_rayl[1] .= 0.5 / size(model.τ_rayl[1], 2)
                model.τ_abs[1] .= 0
                r = rt_run(model)
                (r.toa, r.boa)
            end
        end
        R0, T0 = run(CoreRT.RadauQuad(), 0.0)
        Rδ, Tδ = run(CoreRT.RadauQuad(), 1e-6)
        Rν, Tν = run(CoreRT.RadauQuad(), 1e-5)
        Rε, Tε = run(CoreRT.RadauQuad(), 1e-4)
        Rg, Tg = run(CoreRT.GaussLegQuad(), 0.0)
        @test all(>(0), R0[:, 1, :])
        @test all(>(0), T0[:, 1, :])
        # vza = 0 resolves to the solar ordinate only at SZA = 0, where BOA
        # also carries the unscattered beam; compare the other views. At
        # SZA = 1e-6, μ₀ is within 4 eps of the fixed node μ = 1, which is
        # moved onto μ₀ (single-interval grid). At SZA = 1e-5 and 1e-4 the
        # sub-interval [μ₀, 1] is 1.5e-14 and 1.5e-12 wide, so its weighted
        # nodes are one to a few ulps apart (merged only when exactly equal).
        # The elemental factors do not divide by the node difference (limit
        # form, decision 14e), so the nadir view (μ = 1 node) stays continuous:
        # its change from SZA = 0 is ≈ 1e-11 here, against 2e-6 at SZA = 1e-5
        # with the former 1/(μᵢ - μⱼ) expression.
        for (R, T) in ((Rδ, Tδ), (Rν, Tν), (Rε, Tε))
            @test R[2:end, :, :] ≈ R0[2:end, :, :] rtol=1e-6
            @test T[2:end, :, :] ≈ T0[2:end, :, :] rtol=1e-6
            @test R ≈ R0 rtol=1e-6
            @test R[1, :, :] ≈ R0[1, :, :] rtol=1e-9
        end
        # Different quadrature, same physics.
        @test R0[:, 1, :] ≈ Rg[:, 1, :] rtol=2e-3
        @test T0[:, 1, :] ≈ Tg[:, 1, :] rtol=2e-3
    end
end
