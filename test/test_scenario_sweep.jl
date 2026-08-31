# =========================================================================
# Scenario sweeps: SZA × view-pair × BRDF over one scene
# (see `docs/dev_notes/proposals/surface_split_albedo_sweep.md` §6/§7, PR 3)
# =========================================================================
#
# Two things are under test here:
#
#   1. `remake_geometry` (tools/model_from_parameters.jl) — the cheap
#      geometry-only `RTModel` rebuild that `run_sweep` uses instead of a
#      full `model_from_parameters` re-run per SZA. Only `geometry` and
#      `quad_points` depend on (sza, vza, vaz); `atmosphere`/`optics`/
#      `solver`/`surfaces` are shared by reference. Acceptance criterion:
#      `rt_run(remake_geometry(model, params; sza=θ))` is BIT-EXACT
#      (`==`, not `≈`) with `rt_run(model_from_parameters(p′))` for
#      `p′ = deepcopy(params)` with `p′.sza = θ` — same kernels, same
#      inputs, so any mismatch means a hidden geometry dependency was
#      missed (see `remake_geometry`'s docstring for why `quad_points.
#      Nstreams`/`solver.m_max_bands` are expected to be SZA-invariant).
#   2. `ScenarioSweep`/`SweepResult`/`run_sweep` — the declarative
#      (SZA × view-pair × BRDF) driver built on top of `remake_geometry`
#      (outer loop) + `rt_run_multi_surface` (inner loop, amortizing the
#      atmosphere phase across BRDFs within one SZA). Every `(i,j,k)`
#      slice of `SweepResult.R` must equal a direct monolithic `rt_run`
#      call using the SAME vza/vaz vectors (all view pairs live in one
#      model call — see the module docstring on why this matters for
#      bit-exactness) and the matching BRDF swapped in.
#
# Reuses the cheap Rayleigh-only config from test_surface_split.jl
# (2 spectral points, nstreams=11 ⇒ m_max=2 for a Lambertian surface).

using vSmartMOM
using vSmartMOM.CoreRT
using Test

@testset "Scenario sweep" begin

    arch = vSmartMOM.Architectures.CPU()

    function _base_params()
        p = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        p.architecture = arch
        p.brdf[1] = CoreRT.LambertianSurfaceScalar(0.3)
        return p
    end

    # ---------------------------------------------------------------------
    # 1. remake_geometry equivalence (deliverable B).
    # ---------------------------------------------------------------------
    @testset "remake_geometry equivalence" begin
        params = _base_params()
        model  = model_from_parameters(params)

        for θ in (10.0, 60.0)
            model_geo = CoreRT.remake_geometry(model, params; sza = θ)

            p_ref = deepcopy(params)
            p_ref.sza = θ
            model_ref = model_from_parameters(p_ref)

            # The invariant that makes bit-exactness possible (see
            # `remake_geometry`'s docstring): Nstreams — hence
            # solver.m_max_bands, which `remake_geometry` reuses verbatim
            # instead of recomputing — does not depend on SZA. Assert it
            # explicitly so a future `rt_set_streams` change that breaks
            # this fails loudly here rather than silently degrading.
            @test model_geo.quad_points.Nstreams == model_ref.quad_points.Nstreams
            @test model_geo.solver.m_max_bands == model_ref.solver.m_max_bands

            r_geo = rt_run(model_geo)
            r_ref = rt_run(model_ref)
            @test length(r_geo) == length(r_ref)
            for i in eachindex(r_geo)
                @test r_geo[i] == r_ref[i]
            end
        end

        # `remake_geometry` also accepts vza/vaz overrides (used internally
        # by `run_sweep`) — same bit-exactness contract.
        @testset "remake_geometry with vza/vaz override" begin
            vza_new = [40.0, 15.0, 0.0]
            vaz_new = [0.0, 180.0, 0.0]
            θ = 45.0
            model_geo = CoreRT.remake_geometry(model, params; sza = θ, vza = vza_new, vaz = vaz_new)

            p_ref = deepcopy(params)
            p_ref.sza = θ; p_ref.vza = vza_new; p_ref.vaz = vaz_new
            model_ref = model_from_parameters(p_ref)

            r_geo = rt_run(model_geo)
            r_ref = rt_run(model_ref)
            for i in eachindex(r_geo)
                @test r_geo[i] == r_ref[i]
            end
        end
    end

    # ---------------------------------------------------------------------
    # 2. Small sweep: 2 SZA × 2 paired view-pairs × 2 BRDFs, checked slice-
    #    by-slice against a direct monolithic `rt_run` reference built with
    #    the SAME vza/vaz vectors (run_sweep runs all view pairs in one
    #    model call, so the reference must too, for bit-exactness).
    # ---------------------------------------------------------------------
    @testset "Small sweep vs monolithic reference" begin
        params = _base_params()
        model  = model_from_parameters(params)

        vza_pairs = [40.0, 15.0]
        vaz_pairs = [0.0, 180.0]
        sza_list  = [10.0, 45.0]
        lamb_a = CoreRT.LambertianSurfaceScalar(0.1)
        lamb_b = CoreRT.LambertianSurfaceScalar(0.4)

        sweep = ScenarioSweep(; sza = sza_list, vza = vza_pairs, vaz = vaz_pairs,
                              view_mode = :paired, brdfs = [lamb_a, lamb_b])

        result = run_sweep(model, params, sweep)

        pol_n = CoreRT.polarization_type(model).n
        nSpec = size(model.τ_abs[1], 1)
        @test size(result.R) == (2, 2, 2, pol_n, nSpec)
        @test size(result.T) == (2, 2, 2, pol_n, nSpec)

        for (i_sza, sza_val) in enumerate(sweep.sza), (i_brdf, brdf) in enumerate(sweep.brdfs)
            p_ref = deepcopy(params)
            p_ref.sza = sza_val
            p_ref.vza = vza_pairs
            p_ref.vaz = vaz_pairs
            p_ref.brdf[1] = brdf
            model_ref = model_from_parameters(p_ref)
            ref = rt_run(model_ref)   # (R, T, ieR, ieT, hdr, bhr_uw, bhr_dw), R/T shaped (N_vp, pol_n, nSpec)

            @test result.R[i_sza, :, i_brdf, :, :] == ref[1]
            @test result.T[i_sza, :, i_brdf, :, :] == ref[2]
        end
    end

    # ---------------------------------------------------------------------
    # 3. `view_mode` guard.
    # ---------------------------------------------------------------------
    @testset "view_mode guard" begin
        lamb = CoreRT.LambertianSurfaceScalar(0.1)

        # vza/vaz passed without view_mode -> throws (the footgun guard;
        # NEVER inferred from vector lengths).
        @test_throws ArgumentError ScenarioSweep(; sza = [30.0], vza = [10.0, 20.0],
                                                 vaz = [0.0, 180.0], brdfs = [lamb])

        # Passing both view_pairs AND vza/vaz/view_mode is also rejected.
        @test_throws ArgumentError ScenarioSweep(; sza = [30.0], view_pairs = [(10.0, 0.0)],
                                                 vza = [10.0], vaz = [0.0], view_mode = :paired,
                                                 brdfs = [lamb])

        # Neither view_pairs nor vza/vaz given.
        @test_throws ArgumentError ScenarioSweep(; sza = [30.0], brdfs = [lamb])

        # :paired requires equal lengths.
        @test_throws ArgumentError ScenarioSweep(; sza = [30.0], vza = [10.0, 20.0], vaz = [0.0],
                                                 view_mode = :paired, brdfs = [lamb])

        # Invalid view_mode symbol.
        @test_throws ArgumentError ScenarioSweep(; sza = [30.0], vza = [10.0], vaz = [0.0],
                                                 view_mode = :bogus, brdfs = [lamb])

        # :cross expands to the full Cartesian product, row-major over (vza, vaz).
        sweep_cross = ScenarioSweep(; sza = [30.0], vza = [10.0, 20.0, 30.0], vaz = [0.0, 180.0],
                                    view_mode = :cross, brdfs = [lamb])
        @test length(sweep_cross.view_pairs) == 3 * 2
        @test sweep_cross.view_pairs == [(10.0, 0.0), (10.0, 180.0), (20.0, 0.0),
                                         (20.0, 180.0), (30.0, 0.0), (30.0, 180.0)]

        # :paired keeps the pairing as given.
        sweep_paired = ScenarioSweep(; sza = [30.0], vza = [10.0, 20.0], vaz = [0.0, 180.0],
                                     view_mode = :paired, brdfs = [lamb])
        @test sweep_paired.view_pairs == [(10.0, 0.0), (20.0, 180.0)]
    end

    # ---------------------------------------------------------------------
    # 4. SweepResult shape + axis ordering.
    # ---------------------------------------------------------------------
    @testset "SweepResult shape + axis ordering" begin
        params = _base_params()
        model  = model_from_parameters(params)
        pol_n  = CoreRT.polarization_type(model).n
        nSpec  = size(model.τ_abs[1], 1)

        lamb_a = CoreRT.LambertianSurfaceScalar(0.0)
        lamb_b = CoreRT.LambertianSurfaceScalar(0.2)
        lamb_c = CoreRT.LambertianSurfaceScalar(0.5)

        sweep = ScenarioSweep(; sza = [20.0, 40.0, 55.0], vza = [30.0], vaz = [0.0],
                              view_mode = :paired, brdfs = [lamb_a, lamb_b, lamb_c])
        result = run_sweep(model, params, sweep)

        @test size(result.R) == (3, 1, 3, pol_n, nSpec)   # (N_sza, N_view_pair, N_brdf, pol_n, nSpec)
        @test size(result.T) == size(result.R)
        @test result.sweep === sweep

        # Axis ordering: bump ONE axis at a time and confirm the slice that
        # should change is the only one that does, relative to the (1,1,1)
        # corner — this pins down which array dimension is SZA vs BRDF
        # (the view-pair axis has only one entry here, so it's covered by
        # the small-sweep test above instead).
        @test result.R[1, 1, 1, :, :] != result.R[2, 1, 1, :, :]   # SZA axis moves
        @test result.R[1, 1, 1, :, :] != result.R[1, 1, 2, :, :]   # BRDF axis moves
        @test result.R[1, 1, 1, :, :] != result.R[1, 1, 3, :, :]

        # a=0 Lambertian at every SZA reproduces the Rayleigh-only "black
        # sky" reflectance exactly (cross-check against a direct monolithic
        # run, independent of the small-sweep testset above).
        for (i_sza, sza_val) in enumerate(sweep.sza)
            p_ref = deepcopy(params)
            p_ref.sza = sza_val
            p_ref.vza = [30.0]; p_ref.vaz = [0.0]
            p_ref.brdf[1] = lamb_a
            model_ref = model_from_parameters(p_ref)
            ref = rt_run(model_ref)[1]
            @test result.R[i_sza, :, 1, :, :] == ref
        end
    end
end
