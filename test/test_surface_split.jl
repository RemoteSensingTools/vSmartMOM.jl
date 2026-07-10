# =========================================================================
# Atmosphere/surface split — rt_run_atmosphere / rt_run_surface /
# rt_run_multi_surface (see `proposals/surface_split_albedo_sweep.md`)
# =========================================================================
#
# `rt_run`'s Fourier loop touches the surface at exactly one point per
# moment `m`: after the z-loop (elemental → doubling → interaction over
# `Nz` layers) settles, a surface `AddedLayer` is built and fused into the
# composite with one final `interaction!`. `rt_run_atmosphere` snapshots
# the pre-surface composite per moment (via `atm_snapshot_callback` +
# `stop_after_atmosphere`); `rt_run_surface` replays only the surface step
# against that snapshot. Since both paths call the *same* kernels with the
# *same* inputs, results must be bit-exact (`==`), not merely `isapprox` —
# any non-equality here means the split has drifted from `rt_run.jl`.
#
# Test map:
#   1. Identity — rt_run_multi_surface(model, [model.surfaces[1]])[1] ==
#      rt_run(model), on a small Lambertian Rayleigh-only scene.
#   2. Surface swap — cache built once, replayed against a different
#      Lambertian albedo, vs a monolithic run with `model.surfaces[1]`
#      mutated to the same albedo (surfaces is a plain
#      `Vector{AbstractSurfaceType}` — swapping an element is legal).
#   3. RPV — (a) a model whose own surface IS RPV, replayed against
#      itself; (b) the m_max-widening path: a cache built from a
#      Lambertian model but sized for an RPV `target_brdfs`, replayed and
#      checked against a monolithic RPV run at the SAME Fourier loop
#      bound (asserted first, since a mismatched m_max would silently
#      change the Fourier sum and invalidate the bit-exactness check).
#   4. Guards — CanopySurface rejected (both as the model's own surface
#      and as a `target_brdfs` entry, and at replay time); an
#      m_max-exceeding BRDF replayed against an under-sized cache throws.
#   5. Multi-surface sweep — two Lambertian albedos in one
#      `rt_run_multi_surface` call, each checked against its own
#      monolithic reference.
#   6. Cox-Munk identity — exercises the replayed single-scattering (TMS)
#      correction; more expensive (`m_max=21` for the ocean_coxmunk.yaml
#      config) but still within budget (~30s combined).
#
# NOTE on RPV `l_trunc`: `config/vegetation_rpv.yaml` uses `nstreams: 11`
# (⇒ `l_trunc = 21`), and running RPV to its native `m_max = 21` hits a
# PRE-EXISTING numerical fragility in `interaction!`'s batched `inv`
# (`matrix contains Infs or NaNs`) that reproduces identically on
# unmodified `main` via a plain monolithic `rt_run` — i.e. it is a
# property of the RPV Fourier-reflectance / adding-doubling combination
# at that resolution, not something this PR's split introduces. Verified
# non-monotonic in `l_trunc` (5 and 13 are clean, 7 and 9 are not) —
# consistent with a conditioning issue rather than a "high m is always
# bad" one. Out of scope for this PR (which ports the split mechanism,
# not RPV numerics); the tests below mutate `l_trunc`/`max_m` down to a
# known-stable, cheap value (5) so they exercise real multi-moment RPV
# support without tripping the unrelated fragility.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

@testset "Atmosphere/surface split" begin

    arch = vSmartMOM.Architectures.CPU()

    # ---------------------------------------------------------------------
    # 1. Identity: rt_run_multi_surface(model, [model.surfaces[1]])[1] ==
    #    rt_run(model), bit-exact.
    # ---------------------------------------------------------------------
    @testset "Identity (Lambertian, Rayleigh-only)" begin
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = arch
        params.brdf[1] = CoreRT.LambertianSurfaceScalar(0.3)
        model = model_from_parameters(params)

        ref = rt_run(model)
        res = rt_run_multi_surface(model, [model.surfaces[1]])[1]

        @test length(ref) == length(res)
        for i in eachindex(ref)
            @test ref[i] == res[i]
        end
    end

    # ---------------------------------------------------------------------
    # 2. Surface swap: build the cache once (over the model's own 0.1
    #    albedo), replay against a DIFFERENT albedo (0.3) — Lambertian's
    #    component_m_max is always 0, so any Rayleigh-driven cache can
    #    replay any Lambertian albedo, whether or not it was in
    #    `target_brdfs`. Reference: mutate `model.surfaces[1]` in place
    #    (a plain `Vector{AbstractSurfaceType}`) and run monolithically.
    # ---------------------------------------------------------------------
    @testset "Surface swap (Lambertian albedo)" begin
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = arch
        params.brdf[1] = CoreRT.LambertianSurfaceScalar(0.1)
        model = model_from_parameters(params)

        swapped = CoreRT.LambertianSurfaceScalar(0.3)
        cache = rt_run_atmosphere(model)
        res   = rt_run_surface(cache, swapped)

        model_swapped = deepcopy(model)
        model_swapped.surfaces[1] = swapped
        ref = rt_run(model_swapped)

        for i in eachindex(ref)
            @test ref[i] == res[i]
        end
    end

    # ---------------------------------------------------------------------
    # 3. RPV — same-model replay, then the m_max-widening path.
    # ---------------------------------------------------------------------
    @testset "RPV" begin
        # Small, known-stable Fourier width (see the module-level NOTE
        # above for why 21 — the config's native nstreams-derived value —
        # is avoided).
        function _rpv_params()
            p = parameters_from_yaml("../config/vegetation_rpv.yaml")
            p.architecture = arch
            p.l_trunc = 5
            p.max_m   = 6
            return p
        end

        params_rpv = _rpv_params()
        model_rpv  = model_from_parameters(params_rpv)
        rpv_surf   = model_rpv.surfaces[1]
        @test rpv_surf isa CoreRT.rpvSurfaceScalar
        @test model_rpv.solver.m_max_bands[1] == 5   # sanity: RPV drives the loop bound here

        ref_rpv = rt_run(model_rpv)

        @testset "Same-model replay" begin
            cache_a = rt_run_atmosphere(model_rpv)   # target_brdfs defaults to model.surfaces
            @test cache_a.m_max == model_rpv.solver.m_max_bands[1]
            res_a = rt_run_surface(cache_a, rpv_surf)
            for i in eachindex(ref_rpv)
                @test ref_rpv[i] == res_a[i]
            end
        end

        @testset "m_max widening from a Lambertian-only cache" begin
            params_lamb = _rpv_params()
            params_lamb.brdf[1] = CoreRT.LambertianSurfaceScalar(0.1)
            model_lamb = model_from_parameters(params_lamb)
            @test model_lamb.solver.m_max_bands[1] < model_rpv.solver.m_max_bands[1]

            cache_b = rt_run_atmosphere(model_lamb; target_brdfs = [rpv_surf])
            # Apples-to-apples check FIRST: the cache must have widened to
            # exactly what `component_m_max` demands for RPV — otherwise
            # the replay below would silently sum a truncated Fourier
            # series and the bit-exactness check would be meaningless.
            @test cache_b.m_max == model_rpv.solver.m_max_bands[1]

            res_b = rt_run_surface(cache_b, rpv_surf)
            for i in eachindex(ref_rpv)
                @test ref_rpv[i] == res_b[i]
            end
        end
    end

    # ---------------------------------------------------------------------
    # 4. Guards.
    # ---------------------------------------------------------------------
    @testset "Guards" begin
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = arch
        params.brdf[1] = CoreRT.LambertianSurfaceScalar(0.1)
        model = model_from_parameters(params)

        soil   = CoreRT.LambertianSurfaceScalar(0.1)
        canopy = CoreRT.CanopySurface(; soil = soil, LAI = 1.0)

        # CanopySurface in `target_brdfs`.
        @test_throws ArgumentError rt_run_atmosphere(model; target_brdfs = [canopy])

        # CanopySurface as the model's own surface (default target_brdfs).
        model_canopy = deepcopy(model)
        model_canopy.surfaces[1] = canopy
        @test_throws ArgumentError rt_run_atmosphere(model_canopy)

        # CanopySurface at replay time.
        cache_lamb = rt_run_atmosphere(model)   # sized for Lambertian only (m_max=2)
        @test_throws ArgumentError rt_run_surface(cache_lamb, canopy)

        # An m_max-exceeding BRDF (RPV needs `user_l_cap`, here 21) against
        # a cache sized only for Lambertian (m_max=2) must throw — not
        # silently truncate the Fourier sum.
        rpv = CoreRT.rpvSurfaceScalar(0.1, 0.1, 0.7, -0.1)
        @test cache_lamb.m_max == 2
        @test_throws ArgumentError rt_run_surface(cache_lamb, rpv)
    end

    # ---------------------------------------------------------------------
    # 5. Multi-surface sweep.
    # ---------------------------------------------------------------------
    @testset "Multi-surface sweep" begin
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = arch
        params.brdf[1] = CoreRT.LambertianSurfaceScalar(0.1)
        model = model_from_parameters(params)

        lamb_a = CoreRT.LambertianSurfaceScalar(0.1)
        lamb_b = CoreRT.LambertianSurfaceScalar(0.4)
        results = rt_run_multi_surface(model, [lamb_a, lamb_b])
        @test length(results) == 2

        model_a = deepcopy(model); model_a.surfaces[1] = lamb_a
        model_b = deepcopy(model); model_b.surfaces[1] = lamb_b
        ref_a = rt_run(model_a)
        ref_b = rt_run(model_b)

        for i in eachindex(ref_a)
            @test results[1][i] == ref_a[i]
            @test results[2][i] == ref_b[i]
        end
    end

    # ---------------------------------------------------------------------
    # 6. Cox-Munk identity — exercises the replayed TMS single-scattering
    #    correction (`apply_ss_correction!`, called per-BRDF inside
    #    `rt_run_surface`, mirroring rt_run.jl's post-loop block). Costs
    #    ~30s combined (nstreams=11 ⇒ m_max=21 for Cox-Munk's
    #    `component_m_max = user_l_cap`) but stays within the suite's
    #    per-file budget.
    # ---------------------------------------------------------------------
    @testset "Cox-Munk identity" begin
        params = parameters_from_yaml("../config/ocean_coxmunk.yaml")
        params.architecture = arch
        model = model_from_parameters(params)

        ref = rt_run(model)
        res = rt_run_multi_surface(model, [model.surfaces[1]])[1]

        for i in eachindex(ref)
            @test ref[i] == res[i]
        end
    end
end
