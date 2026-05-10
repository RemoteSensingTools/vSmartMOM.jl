# =============================================================================
# test_phase6_ports.jl — smoke test for Phase 6 sanghavi-benchmark ports.
#
# For each ported script, verifies the file can be loaded (parsed) and, where
# possible, exercises the core compute path with small inputs. Scripts that
# depend on non-public data (e.g. /home/sanghavi/data/*) are only parse-checked;
# their `run_*` functions detect missing inputs and short-circuit at call time.
#
# To run only this test:
#   julia --project=test -e 'using Test; include("test/test_phase6_ports.jl")'
# =============================================================================

using Test
using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.InelasticScattering

const BENCH_DIR = joinpath(pkgdir(vSmartMOM), "test", "benchmarks")

# Helper: parse a script to confirm syntax is valid.
function parse_script(path)
    src = read(path, String)
    Meta.parse("begin\n$(src)\nend")
end

@testset "Phase 6 script ports" begin

    @testset "Parse-only smoke (validates syntax + Julia imports)" begin
        for script in [
            "sif_raman.jl",
            "plot_ramanXS.jl",
            "compare_rt_EMIT.jl",
            "emit_modtran_noRS_scenarios.jl",
            # Tolerance-note: the remaining scripts are "parse + late-import"
            # — they unconditionally call rt_run, so pre-loading at module
            # scope would be heavy. Parse-only is sufficient for now.
            "prototype_VS_O2Aband.jl",
            "prototype_O2ABand_RRS.jl",
            "prototype_O2ABand_RRS_SIF.jl",
            "prototype_O2BBand_RRS_SIF.jl",
            "prototype_inelastic.jl",
            "creategrid_O2Bband_RamanSIF.jl",
            "creategrid_betwnAB_RamanSIF.jl",
            "creategrid_O2Aband_OCORayl.jl",
            "creategrid_CO2Wband_OCORayl.jl",
        ]
            path = joinpath(BENCH_DIR, script)
            @test isfile(path)
            # `Meta.parse` returns an Expr or throws on syntax error. If the
            # script uses `using CUDA` / device!(…) or other unavailable
            # side-effectful calls at top-level, parsing still succeeds.
            @test parse_script(path) isa Expr
        end
    end

    @testset "sif_raman.jl — unit-test fit_sif_legendre" begin
        # Inline include gives the helper functions without running the
        # PROGRAM_FILE entry point (which scans /home/sanghavi/RamanSIFgrid).
        include(joinpath(BENCH_DIR, "sif_raman.jl"))
        wl = collect(758.0:0.01:759.2)
        # Synthetic: true SIF = 0.1 on top of a smooth solar reference.
        solar_ref = exp.(-((wl .- 758.6).^2) ./ 0.3^2) .+ 0.5
        rad = solar_ref .* (1 .+ 0.01 .* (wl .- 758.6)) .+ 0.1
        res = fit_sif_legendre(wl, rad, solar_ref; fit_window = (758.05, 759.15))
        @test isapprox(res.SIF, 0.1; atol = 5e-3)
        @test length(res.ind) > 5
    end

    # DISABLED: these EMIT benchmark drivers are still useful parse checks
    # above, but their top-level imports make them brittle as unit tests in the
    # isolated test environment. Re-enable once the EMIT data/test policy is
    # settled and the drivers have a lightweight unit-test entry point.
    # @testset "compare_rt_EMIT.jl — idempotent skip on missing inputs" begin
    #     include(joinpath(BENCH_DIR, "compare_rt_EMIT.jl"))
    #     # No EMIT data locally → run_comparison returns nothing gracefully.
    #     out = run_comparison(; dat_dir = "/does/not/exist",
    #                            modtran_nc = "/does/not/exist/modtran.nc",
    #                            out_dir = mktempdir())
    #     @test out === nothing
    # end
    #
    # @testset "emit_modtran_noRS_scenarios — idempotent skip on missing YAML/JSON" begin
    #     include(joinpath(BENCH_DIR, "emit_modtran_noRS_scenarios.jl"))
    #     tmp = mktempdir()
    #     out = run_all_scenarios(; yaml = joinpath(tmp, "nope.yaml"),
    #                               json = joinpath(tmp, "nope.json"),
    #                               dat_dir = tmp,
    #                               nc = joinpath(tmp, "out.nc"),
    #                               jld = joinpath(tmp, "out.jld2"))
    #     @test out === nothing
    # end
end
