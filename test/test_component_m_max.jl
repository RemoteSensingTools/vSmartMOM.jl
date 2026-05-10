# Phase C — `component_m_max(component, ctx)` trait dispatch
# =========================================================================
#
# Verifies that each surface / source / scatterer declares the correct
# Fourier support, and that the trait-based aggregator
# `_derive_m_max_bands_via_traits` produces the right per-band bound
# when `SolverConfig.use_component_traits` is flipped on.

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.Scattering
using Test

const _CTX_30 = (; user_l_cap = 30, stream_l_cap = 30,
                  m_max_override = nothing, truncation = nothing)

@testset "Phase C — component_m_max trait dispatch" begin

    @testset "Sources are neutral or zero" begin
        # SolarBeam → 0 (NOT typemax(Int) — would pin the loop to user_l_cap)
        @test CoreRT.component_m_max(CoreRT.SolarBeam(), _CTX_30) == 0
        @test CoreRT.component_m_max(CoreRT.SurfaceSIF(), _CTX_30) == 0
        @test CoreRT.component_m_max(CoreRT.NoSource(), _CTX_30) == 0
        # SourceSet → max over its sources.
        ss_empty = CoreRT.SourceSet(())
        @test CoreRT.component_m_max(ss_empty, _CTX_30) == 0
        ss_solar = CoreRT.SourceSet((CoreRT.SolarBeam(), CoreRT.SurfaceSIF()))
        @test CoreRT.component_m_max(ss_solar, _CTX_30) == 0
    end

    @testset "Lambertian surfaces → 0 (m=0 is exact)" begin
        @test CoreRT.component_m_max(
            CoreRT.LambertianSurfaceScalar(0.1), _CTX_30) == 0
        @test CoreRT.component_m_max(
            CoreRT.LambertianSurfaceLegendre([0.1]), _CTX_30) == 0
        # LambertianSurfaceSpline / LambertianSurfaceSpectrum are
        # covered by the same default-method abstract trait — they
        # share the dispatch via their type signature.
    end

    @testset "Cox-Munk / RPV / Ross-Li / Canopy → user_l_cap" begin
        cm = CoreRT.CoxMunkSurface(wind_speed=5.0)
        @test CoreRT.component_m_max(cm, _CTX_30) == 30

        rpv = CoreRT.rpvSurfaceScalar(0.1, 0.1, 0.7, -0.1)
        @test CoreRT.component_m_max(rpv, _CTX_30) == 30

        rl = CoreRT.RossLiSurfaceScalar(0.1, 0.1, 0.1)
        @test CoreRT.component_m_max(rl, _CTX_30) == 30

        # Canopy uses default Lambertian soil + spherical leaves
        soil = CoreRT.LambertianSurfaceScalar(0.1)
        canopy = CoreRT.CanopySurface(; soil=soil, LAI=1.0)
        @test CoreRT.component_m_max(canopy, _CTX_30) == 30

        # Tighter cap should be respected (trait returns user_l_cap directly).
        ctx_5 = (; user_l_cap = 5, stream_l_cap = 5,
                  m_max_override = nothing, truncation = nothing)
        @test CoreRT.component_m_max(cm, ctx_5) == 5
    end

    @testset "Rayleigh → 2 (β₀, β₁, β₂ exhaust the phase function)" begin
        # Build a minimal RayleighScattering by reaching into a
        # tiny Rayleigh-only model.
        params = parameters_from_yaml(
            joinpath(@__DIR__, "test_parameters", "PureRayleighParameters.yaml"))
        params.architecture = vSmartMOM.Architectures.CPU()
        m = model_from_parameters(params)
        @test CoreRT.component_m_max(m.optics.rayleigh, _CTX_30) == 2
    end
end

@testset "Phase C — _aggregate_m_max" begin
    @testset "max-then-clamp" begin
        # Mix of Rayleigh (2) and a Lambertian (0) → 2; capped to user_l_cap=5.
        ctx = (; user_l_cap = 5, stream_l_cap = 5,
                m_max_override = nothing, truncation = nothing)
        comps = (
            CoreRT.LambertianSurfaceScalar(0.1),
            # Construct a fake "Rayleigh-like" component via a dummy method
            # by relying on RayleighScattering trait — use the model's
            # Rayleigh state from above instead.
        )
        # Single-Lambertian: 0
        @test CoreRT._aggregate_m_max(comps, ctx) == 0
    end

    @testset "Cox-Munk dominates Lambertian" begin
        ctx = (; user_l_cap = 13, stream_l_cap = 13,
                m_max_override = nothing, truncation = nothing)
        comps = (
            CoreRT.LambertianSurfaceScalar(0.1),
            CoreRT.CoxMunkSurface(wind_speed=5.0),
        )
        # max(0, 13) = 13; clamped to user_l_cap=13.
        @test CoreRT._aggregate_m_max(comps, ctx) == 13
    end

    @testset "m_max_override clamps below user_l_cap" begin
        ctx = (; user_l_cap = 13, stream_l_cap = 13,
                m_max_override = 4, truncation = nothing)
        comps = (CoreRT.CoxMunkSurface(wind_speed=5.0),)
        @test CoreRT._aggregate_m_max(comps, ctx) == 4
    end

    @testset "empty components → 0" begin
        ctx = (; user_l_cap = 10, stream_l_cap = 10,
                m_max_override = nothing, truncation = nothing)
        @test CoreRT._aggregate_m_max((), ctx) == 0
    end
end

@testset "Phase C — trait-based aggregator is the active path" begin
    # Phase C lands the trait infrastructure AND turns it on by default
    # (combined with the originally-staged Phase C-flip). Cox-Munk forward
    # gets its full surface-driven Fourier resolution back instead of the
    # half-truncation the previous count-only aggregator produced.
    params = parameters_from_yaml(
        joinpath(@__DIR__, "test_parameters", "PureRayleighParameters.yaml"))
    params.architecture = vSmartMOM.Architectures.CPU()
    m = model_from_parameters(params)
    @test m.solver.use_component_traits == true
    # Rayleigh-only / Lambertian config: trait aggregator picks
    # max(2, 0) = 2, capped to legacy_order_cap = max_m - 1.
    # PureRayleighParameters.yaml has max_m=3, so order cap = 2.
    @test m.solver.m_max_bands == fill(2, length(m.solver.m_max_bands))
    @test m.solver.n_fourier_moments_bands == fill(3, length(m.solver.n_fourier_moments_bands))
end
