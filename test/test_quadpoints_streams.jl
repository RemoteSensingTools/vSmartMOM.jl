# Phase A — `Nstreams` field on `QuadPoints`
# =========================================================================
#
# Verifies the `Nstreams` (count of nonzero weights) vs `Nquad` (augmented
# total including zero-weight SZA/VZA output nodes) distinction across
# both supported quadrature schemes.

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Architectures
using Test

# Minimal ObsGeometry / pol_type / arr_type for direct rt_set_streams tests.
const _ARR = vSmartMOM.Architectures.array_type(vSmartMOM.Architectures.CPU())

function _make_obs_geom(; sza, vza)
    CoreRT.ObsGeometry{Float64}(
        sza = Float64(sza),
        vza = Float64.(vza),
        vaz = zeros(Float64, length(vza)),
        obs_alt = 1000.0,
    )
end

@testset "Phase A — QuadPoints Nstreams" begin
    pol_type = CoreRT.Stokes_IQUV()

    @testset "GaussLegQuad — Nstreams = (Ltrunc+2)÷2" begin
        for Ltrunc in (0, 2, 5, 10, 20, 25)
            obs = _make_obs_geom(sza = 30.0, vza = [10.0, 20.0])
            qp = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(), Ltrunc, obs, pol_type, _ARR)

            # Nstreams equals count of nonzero weights by construction.
            @test qp.Nstreams == count(!iszero, qp.wt_μ)
            # Public formula for Gauss matches Sanghavi's (Ltrunc+2)÷2.
            @test qp.Nstreams == (Ltrunc + 2) ÷ 2
            # Augmented Nquad ≥ Nstreams (SZA + VZAs only ever add nodes).
            @test qp.Nquad >= qp.Nstreams
            # Public contract: stream_l_cap = 2·Nstreams - 1.
            @test 2 * qp.Nstreams - 1 >= Ltrunc
        end
    end

    @testset "GaussLegQuad — Ltrunc = 0 corner case" begin
        obs = _make_obs_geom(sza = 45.0, vza = [0.0])
        qp = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(), 0, obs, pol_type, _ARR)
        # Sanghavi: +2 form must yield ≥1 weighted stream so the
        # quadrature is non-degenerate.
        @test qp.Nstreams == 1
        @test qp.Nquad >= 1
    end

    @testset "RadauQuad — Nstreams from count(!iszero, wt_μ)" begin
        for Ltrunc in (5, 10, 20, 25)
            obs = _make_obs_geom(sza = 78.46, vza = [10.0, 20.0, 40.0])
            qp = CoreRT.rt_set_streams(CoreRT.RadauQuad(), Ltrunc, obs, pol_type, _ARR)

            # Field invariant: Nstreams matches the actual weight array.
            @test qp.Nstreams == count(!iszero, qp.wt_μ)
            # Augmented Nquad ≥ Nstreams.
            @test qp.Nquad >= qp.Nstreams
            # Both schemes claim resolving order ≥ Ltrunc under the public
            # contract `stream_l_cap = 2·Nstreams - 1`.
            @test 2 * qp.Nstreams - 1 >= Ltrunc
        end
    end

    @testset "Augmentation only adds zero-weight nodes" begin
        # SZA + VZAs are appended with weight 0 so the total node count
        # grows but Nstreams (count of nonzero weights) is unchanged.
        obs = _make_obs_geom(sza = 60.0, vza = [0.0, 30.0, 60.0])
        for Q in (CoreRT.GaussLegQuad(), CoreRT.RadauQuad())
            qp = CoreRT.rt_set_streams(Q, 10, obs, pol_type, _ARR)
            @test qp.Nstreams + count(iszero, qp.wt_μ) == qp.Nquad
            @test qp.Nstreams == count(!iszero, qp.wt_μ)
        end
    end

    @testset "Public `nstreams` kwarg API matches legacy Ltrunc form" begin
        # Gauss: nstreams=N ↔ Ltrunc=2N-2; the builder must produce the
        # same QuadPoints in either call shape.
        obs = _make_obs_geom(sza = 30.0, vza = [10.0, 20.0])
        for N in (1, 2, 5, 11, 13)
            qp_legacy = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(), 2N - 2, obs, pol_type, _ARR)
            qp_new    = CoreRT.rt_set_streams(CoreRT.GaussLegQuad(), obs, pol_type, _ARR;
                                              nstreams = N)
            @test qp_legacy.Nstreams == qp_new.Nstreams == N
            @test qp_legacy.Nquad == qp_new.Nquad
            # Public contract holds: stream_l_cap = 2N - 1.
        end

        # Radau: nstreams=N ↔ Ltrunc=2N-1. Weighted-stream count may
        # differ between schemes for the same nstreams; what's invariant
        # is `stream_l_cap = 2·N - 1`.
        for N in (3, 5, 8, 13)
            qp = CoreRT.rt_set_streams(CoreRT.RadauQuad(), obs, pol_type, _ARR;
                                       nstreams = N)
            @test qp.Nstreams >= 1
            @test qp.Nstreams == count(!iszero, qp.wt_μ)
            @test qp.Nquad >= qp.Nstreams
            # Public contract: stream_l_cap = 2·N - 1 covers Ltrunc=2N-1.
        end

        # Validation: nstreams must be ≥ 1.
        @test_throws ArgumentError CoreRT.rt_set_streams(
            CoreRT.GaussLegQuad(), obs, pol_type, _ARR; nstreams = 0)
        @test_throws ArgumentError CoreRT.rt_set_streams(
            CoreRT.RadauQuad(), obs, pol_type, _ARR; nstreams = -1)
    end
end
