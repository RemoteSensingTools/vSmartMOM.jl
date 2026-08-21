using Test
using GCHPAerosolOptics

# Regression baseline: port these from the gchp-io branch as the first green tests.
#   gchp-io:test/test_aerosol_generic_tomas_aod.jl   — TOMAS-15 column AOD
#   gchp-io:test/test_layer_resolved_aerosol_optics.jl — layer-resolved optics

@testset "GCHPAerosolOptics" begin
    @testset "Stage 1 — column AOD (port from gchp-io)" begin
        # TODO: read a small GCHP fixture, compute column_aod, assert vs GCHP's AOD diagnostic.
        @test_skip false
    end

    @testset "Stage 2 — representations" begin
        # TODO: lognormal fit recovers a known synthetic distribution; spline matches bins
        #       in the fine-bin limit.
        @test_skip false
    end

    @testset "Stage 3 — layer optics → AerosolOptics" begin
        # TODO: per-bin-summed optics conserve τ; ϖ ∈ [0,1]; Greek β[1] normalization holds.
        @test_skip false
    end

    @testset "ML reduction" begin
        # TODO: eof_expand(eof_reduce(p)) ≈ p within retained-variance tolerance.
        @test_skip false
    end
end
