# =========================================================================
# LambertianClosure — analytic O(1)-per-albedo Lambertian surface closure
# (docs/dev_notes/proposals/surface_split_albedo_sweep.md §3) + the :slim AtmosphereRTCache
# (§4).
# =========================================================================
#
# `lambertian_closure(cache)` reads only the cached m=0 atmosphere blocks
# (§3.1's S̄, E_dw, t_up) and reuses the already-tested `rt_run_surface`
# replay for the a=0 ("black-sky") term. Evaluating the closure at a new
# albedo is then a plain broadcast — no surface build, no adding-step
# matrix inverse — so it should agree with the *exact* `rt_run_surface`
# replay to within the gap between Sherman-Morrison (closure) and batched-LU
# (replay) floating-point paths, not bit-for-bit.
#
# Test map:
#   1. Closure vs replay — scalar albedo, several values.
#   2. Closure vs replay — spectral (per-point) albedo vector.
#   3. Jacobian — vs central FD of the closure itself (tight) and vs central
#      FD of the full replay (looser: it's now also comparing against the
#      Sherman-Morrison/LU gap).
#   4. Inversion round-trip (invert_albedo ∘ evaluate ≈ identity), scalar
#      and spectral.
#   5. a=0 sanity — the closure term must vanish IDENTICALLY (`==`, not
#      `≈`) since coeff(0) = 0 exactly.
#   6. Guards — albedo domain, per-source J₀ slots (thermal), an active
#      SurfaceSIF source (flows through the legacy j₀⁻ slot, not a
#      per-source slot, so it needs its own check — see the closure's
#      docstring), invalid `cache_mode`, `i_vza` out of range, and the
#      `invert_albedo` NaN branch on a synthetic degenerate closure.
#   7. Slim cache (§4) — replay bit-exact vs `:full`; `:auto` picks `:slim`
#      for an all-Lambertian target; a non-Lambertian BRDF replayed against
#      a `:slim` cache throws; the closure builds fine on a `:slim` cache.
#
# Reuses the cheap Rayleigh-only config from test_surface_split.jl
# (2 spectral points, 9 VZAs, nstreams=11 ⇒ m_max=2 for a Lambertian
# surface) to keep wall time low.

using vSmartMOM
using vSmartMOM.CoreRT
using Test

@testset "Lambertian closure" begin

    arch = vSmartMOM.Architectures.CPU()

    function _rayleigh_model(; albedo = 0.3)
        params = parameters_from_yaml("test_parameters/PureRayleighParameters.yaml")
        params.architecture = arch
        params.brdf[1] = CoreRT.LambertianSurfaceScalar(albedo)
        return model_from_parameters(params)
    end

    model = _rayleigh_model()
    nSpec = size(model.τ_abs[1], 1)
    @test nSpec == 2   # sanity: the spectral-vector test below assumes this

    # `cache_mode` defaults to `:auto`, which resolves to `:slim` here (every
    # target BRDF — the model's own Lambertian surface — has
    # component_m_max == 0). This means every test below except the explicit
    # "Slim cache" testset is ALSO exercising the :slim replay path — a
    # deliberate extra check that :slim is truly transparent to
    # `lambertian_closure` and `rt_run_surface`.
    cache = rt_run_atmosphere(model)
    @test cache.cache_mode == :slim
    c = lambertian_closure(cache)

    # ---------------------------------------------------------------------
    # 1. Closure vs replay — scalar albedo.
    # ---------------------------------------------------------------------
    @testset "Closure vs replay (scalar albedo)" begin
        max_rel_dev = 0.0
        for a in (0.0, 0.05, 0.3, 0.8)
            R_c   = c(a)
            R_ref = rt_run_surface(cache, CoreRT.LambertianSurfaceScalar(a))[1]
            @test isapprox(R_c, R_ref; rtol = 1e-12, atol = 1e-12)
            dev = maximum(abs.(R_c .- R_ref) ./ max.(abs.(R_ref), eps()))
            max_rel_dev = max(max_rel_dev, dev)
        end
        # Loud tripwire for a convention/constant-factor error (proposal §3.1's
        # "IMPORTANT constants caveat") — the isapprox above alone could pass
        # at a systematically-relaxed tolerance without this second check.
        @test max_rel_dev < 1e-10
    end

    # ---------------------------------------------------------------------
    # 2. Closure vs replay — spectral (per-point) albedo vector.
    # ---------------------------------------------------------------------
    @testset "Closure vs replay (spectral albedo)" begin
        a_vec = [0.0, 0.9]   # spans [0, 0.9]; nSpec == 2 for this config
        @test length(a_vec) == cache.nSpec
        R_c   = c(a_vec)
        R_ref = rt_run_surface(cache, CoreRT.LambertianSurfaceSpectrum(a_vec))[1]
        @test isapprox(R_c, R_ref; rtol = 1e-12, atol = 1e-12)
        dev = maximum(abs.(R_c .- R_ref) ./ max.(abs.(R_ref), eps()))
        @test dev < 1e-10
    end

    # ---------------------------------------------------------------------
    # 3. Jacobian.
    # ---------------------------------------------------------------------
    @testset "Jacobian" begin
        a0 = 0.3
        h  = 1e-6
        J  = albedo_jacobian(c, a0)

        # Tight: central FD of the closure itself — same analytic function,
        # so the only disagreement is O(h²) truncation + O(eps/h) roundoff.
        Jfd_closure = (c(a0 + h) .- c(a0 - h)) ./ (2h)
        @test isapprox(J, Jfd_closure; rtol = 1e-5)

        # Looser: central FD of the full replay — now also crossing the
        # Sherman-Morrison (closure) vs batched-LU (replay) floating-point gap.
        Jfd_replay = (rt_run_surface(cache, CoreRT.LambertianSurfaceScalar(a0 + h))[1] .-
                      rt_run_surface(cache, CoreRT.LambertianSurfaceScalar(a0 - h))[1]) ./ (2h)
        @test isapprox(J, Jfd_replay; rtol = 1e-4)
    end

    # ---------------------------------------------------------------------
    # 4. Inversion round-trip.
    # ---------------------------------------------------------------------
    @testset "Inversion round-trip" begin
        for a in (0.0, 0.05, 0.3, 0.8)
            R_meas = c(a)
            a_rec  = invert_albedo(c, R_meas)
            @test all(isapprox.(a_rec, a; atol = 1e-9))
        end

        a_vec     = [0.2, 0.7]
        R_meas    = c(a_vec)
        a_rec_vec = invert_albedo(c, R_meas)
        @test isapprox(a_rec_vec, a_vec; atol = 1e-9)
    end

    # ---------------------------------------------------------------------
    # 5. a=0 sanity — the closure term must vanish IDENTICALLY.
    # ---------------------------------------------------------------------
    @testset "a=0 sanity" begin
        R_black_direct = rt_run_surface(cache, CoreRT.LambertianSurfaceScalar(0.0))[1]
        @test c(0.0) == R_black_direct
        @test c(0.0) == c.R_black
    end

    # ---------------------------------------------------------------------
    # 6. Guards.
    # ---------------------------------------------------------------------
    @testset "Guards" begin
        @test_throws ArgumentError c(-0.1)
        @test_throws ArgumentError albedo_jacobian(c, -0.1)

        # a·S̄ ≥ 1: the rank-one adding series must converge.
        a_bad = 1 / minimum(c.S̄) + 1
        @test_throws ArgumentError c(a_bad)

        # Spectral-vector length mismatch.
        @test_throws DimensionMismatch c([0.1, 0.2, 0.3])

        # Invalid cache_mode.
        @test_throws ArgumentError rt_run_atmosphere(model; cache_mode = :bogus)

        # i_vza out of range.
        @test_throws ArgumentError invert_albedo(c, c(0.3); i_vza = 1000)

        # Per-source J₀ slot present (thermal emission declares a source_key
        # and gets its own composite slot) — the closure's solar-only guard.
        pol_type = CoreRT.polarization_type(model)
        Nz = length(model.profile.p_full)
        FT = CoreRT.float_type(model)
        B = fill(FT(1.0), Nz, nSpec)
        cache_thermal = rt_run_atmosphere(model; sources = SolarBeam() + ThermalEmission(B_layer = B))
        @test !isempty(cache_thermal.J₀_by_src_per_m[1])
        @test_throws ArgumentError lambertian_closure(cache_thermal)

        # An active SurfaceSIF source flows through the legacy j₀⁻ slot (not
        # a per-source slot), so it needs its own guard — see the closure's
        # docstring for why (†) doesn't capture its albedo-dependent bounce.
        sif_spec = zeros(FT, pol_type.n, nSpec); sif_spec[1, :] .= FT(0.01)
        cache_sif = rt_run_atmosphere(model; sources = SolarBeam() + SurfaceSIF(SIF₀ = sif_spec))
        @test isempty(cache_sif.J₀_by_src_per_m[1])   # confirms the OTHER guard wouldn't catch this
        @test_throws ArgumentError lambertian_closure(cache_sif)

        # invert_albedo's NaN branch (near-zero t_upᴵ / E_dw) — exercised on
        # a synthetic degenerate closure since the real scene never produces
        # exactly-zero transmission.
        R_black0 = zeros(FT, 1, 1, 2)
        t_up0    = zeros(FT, 1, 1, 2)   # t_upᴵ == 0 everywhere -> NaN guard
        c_degenerate = CoreRT.LambertianClosure{FT}(R_black0, t_up0, FT[0.1, 0.1], FT[0.5, 0.5], FT(0.5 / π))
        @test all(isnan, invert_albedo(c_degenerate, R_black0))

        # Shape guards: a wrong spectral length must raise a clear
        # DimensionMismatch, NOT silently index out of bounds (the loop is
        # @inbounds), and an out-of-range i_vza an ArgumentError.
        @test_throws DimensionMismatch invert_albedo(c_degenerate, zeros(FT, 1, 1, 3))
        @test_throws ArgumentError invert_albedo(c_degenerate, R_black0; i_vza = 2)
    end

    # ---------------------------------------------------------------------
    # 7. Slim cache (§4).
    # ---------------------------------------------------------------------
    @testset "Slim cache" begin
        cache_full = rt_run_atmosphere(model; cache_mode = :full)
        cache_slim = rt_run_atmosphere(model; cache_mode = :slim)
        @test cache_full.cache_mode == :full
        @test cache_slim.cache_mode == :slim

        # Bit-exact replay agreement, :full vs :slim, at two albedos.
        for a in (0.1, 0.4)
            r_full = rt_run_surface(cache_full, CoreRT.LambertianSurfaceScalar(a))
            r_slim = rt_run_surface(cache_slim, CoreRT.LambertianSurfaceScalar(a))
            @test length(r_full) == length(r_slim)
            for i in eachindex(r_full)
                @test r_full[i] == r_slim[i]
            end
        end

        # :auto already checked at the top of this file (`cache.cache_mode
        # == :slim`) — this is the same introspection, restated here next to
        # the rest of the slim-cache checks for locality.
        @test rt_run_atmosphere(model).cache_mode == :slim

        # Non-Lambertian BRDF against a :slim cache throws (the cache's own
        # m_max is too narrow for RPV here, so this already throws via the
        # m_max guard — see the widened-m_max case just below for the
        # slim-mode-specific guard in isolation).
        rpv = CoreRT.rpvSurfaceScalar(0.1, 0.1, 0.7, -0.1)
        @test_throws ArgumentError rt_run_surface(cache_slim, rpv)

        # Force cache_mode=:slim while widening m_max via target_brdfs=[rpv],
        # so the m_max guard alone would NOT catch it — isolates the new
        # slim-mode-specific guard in `rt_run_surface`.
        cache_forced_slim = rt_run_atmosphere(model; target_brdfs = [rpv], cache_mode = :slim)
        @test cache_forced_slim.cache_mode == :slim
        @test cache_forced_slim.m_max >= CoreRT.component_m_max(rpv,
            (; user_l_cap = cache_forced_slim.user_l_cap,
               stream_l_cap = 2 * cache_forced_slim.quad_points.Nstreams - 1,
               truncation = nothing, m_max_override = nothing))
        @test_throws ArgumentError rt_run_surface(cache_forced_slim, rpv)

        # The closure works on a :slim cache (already exercised by every
        # other testset above via `cache`/`c`, since :auto resolved to
        # :slim) — also check it against an explicitly-built :slim cache and
        # against the :full cache for numerical agreement.
        c_slim = lambertian_closure(cache_slim)
        c_full = lambertian_closure(cache_full)
        @test c_slim(0.3) == c_full(0.3)
        @test c_slim.R_black == c_full.R_black
        @test c_slim.t_up == c_full.t_up
        @test c_slim.S̄ == c_full.S̄
        @test c_slim.E_dw == c_full.E_dw
    end
end
