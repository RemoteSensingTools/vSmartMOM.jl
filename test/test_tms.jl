# =========================================================================
# TMS single-scattering correction (AbstractSingleScatteringCorrection)
# =========================================================================
#
# Covers: strategy defaults, the Taylor-free path factor, the f = 0 no-op
# (memo §3.5 — mask and exact addition cancel at machine precision for
# Rayleigh, where truncated ≡ exact), the parity-policy rejections, the
# exact-SS reference, and the truncation-order convergence diagnostic
# (memo §5.7 — TMS-corrected results at two δBGE orders agree better than
# uncorrected ones; the remaining gap bounds the out-of-scope IMS term).

using vSmartMOM
using vSmartMOM.CoreRT
using vSmartMOM.CoreRT: NoSSCorrection, TMSCorrection, RTNumericalParameters,
                        _ss_path_factor, rt_run_ss_exact, LambertianSurfaceScalar
using vSmartMOM: Scattering
using vSmartMOM.Scattering: δBGE, Aerosol
using Distributions
using ForwardDiff
using Test

const _TMS_TD = @__DIR__

_tms_set_numerics!(p, sc) = (p.numerics = RTNumericalParameters{p.float_type}(
    dτ_max_threshold    = p.numerics.dτ_max_threshold,
    dτ_min_floor        = p.numerics.dτ_min_floor,
    blas_threads        = p.numerics.blas_threads,
    verbose             = p.numerics.verbose,
    fourier_convergence = p.numerics.fourier_convergence,
    ss_correction       = sc); p)

@testset "strategy defaults" begin
    @test RTNumericalParameters{Float64}().ss_correction isa NoSSCorrection
    @test TMSCorrection() isa CoreRT.AbstractSingleScatteringCorrection
end

@testset "path factor: expm1 limit and Dual safety" begin
    μv, μ₀ = 0.8, 0.6
    s = 1 / μv + 1 / μ₀
    # Small-Δτ limit: pf → exp(−τa·s)·μ₀/(μv+μ₀)·Δτ·s, exact via expm1.
    τa, Δτ = 0.3, 1e-12
    pf = _ss_path_factor(τa, Δτ, μv, μ₀)
    @test pf ≈ exp(-τa * s) * μ₀ / (μv + μ₀) * Δτ * s rtol = 1e-10
    # Float32 does not collapse to 0/garbage at tiny Δτ.
    pf32 = _ss_path_factor(0.3f0, 1f-6, 0.8f0, 0.6f0)
    @test isfinite(pf32) && pf32 > 0
    # ForwardDiff derivative w.r.t. Δτ matches the analytic limit at Δτ→0.
    d = ForwardDiff.derivative(x -> _ss_path_factor(τa, x, μv, μ₀), 0.0)
    @test d ≈ exp(-τa * s) * μ₀ / (μv + μ₀) * s rtol = 1e-12
end

@testset "TMS f = 0 no-op (Rayleigh-only)" begin
    params = parameters_from_yaml(joinpath(_TMS_TD,
        "test_parameters/PureRayleighParameters.yaml"))
    params.brdf[1] = LambertianSurfaceScalar(0.2)
    m_off = model_from_parameters(_tms_set_numerics!(deepcopy(params), NoSSCorrection()))
    m_on  = model_from_parameters(_tms_set_numerics!(deepcopy(params), TMSCorrection()))
    R_off, T_off = rt_run(m_off)
    R_on,  T_on  = rt_run(m_on)
    # The downwelling field is untouched by design (mask acts on j₀⁻ only).
    @test T_on == T_off
    # Rayleigh: truncated ≡ exact, so the masked beam term and the exact
    # addition cancel — the doubling recursion composes the first-order
    # term analytically, so this holds at machine precision (memo §3.5).
    @test R_on[:, 1, :] ≈ R_off[:, 1, :] rtol = 1e-10
    @test maximum(abs.(R_on .- R_off)) < 1e-9 * maximum(abs.(R_off))
end

@testset "parity-policy rejections" begin
    params = parameters_from_yaml(joinpath(_TMS_TD,
        "test_parameters/PureRayleighParameters.yaml"))
    params.brdf[1] = LambertianSurfaceScalar(0.2)
    _tms_set_numerics!(params, TMSCorrection())

    # Linearized driver: matched-or-rejected, never silent.
    model_fwd, model_lin = model_from_parameters(LinMode(), deepcopy(params))
    @test_throws ArgumentError rt_run(model_fwd, model_lin, 0,
                                      size(model_lin.τ̇_abs[1], 1), 1; i_band = 1)

    # Interior sensor heights (multisensor path bypasses the column core).
    params_ms = deepcopy(params)
    params_ms.obs_alt = [0.0, 5.0]
    model_ms = model_from_parameters(params_ms; external_solar = false)
    @test_throws ArgumentError rt_run(model_ms)
end

@testset "exact-SS reference is finite, positive, bounded by the full field" begin
    params = parameters_from_yaml(joinpath(_TMS_TD,
        "test_parameters/PureRayleighParameters.yaml"))
    params.brdf[1] = LambertianSurfaceScalar(0.0)   # black surface
    model = model_from_parameters(params)
    ss = rt_run_ss_exact(model)
    R, _ = rt_run(model)
    @test all(isfinite, ss)
    @test all(ss[:, 1, :] .> 0)
    # First order can never exceed the full multiple-scattering field over
    # a black surface.
    @test all(ss[:, 1, :] .<= R[:, 1, :])
end

@testset "TMS through the split and the Lambertian closure" begin
    # The exact-SS term is surface-independent, computed once at cache-build
    # time, added by rt_run_surface — and therefore folded into the
    # closure's R_black automatically. All three routes must agree with the
    # monolithic TMS run bit-exactly (identical mask, identical addition).
    params = parameters_from_yaml(joinpath(_TMS_TD,
        "test_parameters/PureRayleighParameters.yaml"))
    params.brdf[1] = LambertianSurfaceScalar(0.2)
    _tms_set_numerics!(params, TMSCorrection())
    model = model_from_parameters(deepcopy(params))

    ref = rt_run(model)                                   # monolithic TMS
    cache = rt_run_atmosphere(model)
    @test cache.tms_ss !== nothing
    res = rt_run_surface(cache, model.surfaces[1])        # replay TMS
    for i in eachindex(ref)
        @test ref[i] == res[i]
    end

    clos = lambertian_closure(cache)                      # closure TMS
    R_clos = clos(0.2)
    @test R_clos ≈ ref[1] rtol = 1e-12

    # A NoSSCorrection cache stores no term and stays historical.
    params_off = deepcopy(params)
    _tms_set_numerics!(params_off, NoSSCorrection())
    cache_off = rt_run_atmosphere(model_from_parameters(params_off))
    @test cache_off.tms_ss === nothing
end

@testset "exact engine vs full-kernel Fourier SS (aerosol)" begin
    # Independent validation of the exact-SS engine (rotations, scattering
    # angle, normalization, path algebra): for a moderately peaked kernel
    # the Fourier-space SS solver with NoTruncation and ample streams
    # converges to the real-space exact SS. Black surface isolates the
    # first-order field.
    p = parameters_from_yaml(joinpath(_TMS_TD,
        "test_parameters/JacobianTestFast.yaml"))
    p.brdf[1] = LambertianSurfaceScalar(0.0)
    p.scattering_params.rt_aerosols[1].aerosol =
        Aerosol(LogNormal(log(1.0), log(1.6)), 1.45, 0.001)
    p.scattering_params.rt_aerosols[1].τ_ref = 0.3
    p.sza = 70.0; p.vza = [65.0, 30.0]; p.vaz = [0.0, 0.0]
    p.truncation = Scattering.NoTruncation()
    p.l_trunc = 40
    p.max_m   = 41
    m = model_from_parameters(p)
    R_fourier = rt_run_ss(m)[1]
    R_exact   = rt_run_ss_exact(m)
    @test maximum(abs.(R_fourier[:, 1, :] .- R_exact[:, 1, :]) ./
                  abs.(R_exact[:, 1, :])) < 1e-3   # measured 1.5e-4
end

@testset "truncation-order convergence (memo §5.7)" begin
    # Coarse aerosol with a real forward peak; FIXED streams (l_trunc) so
    # only the δBGE truncation order changes between the two builds. At
    # backscatter-side geometries (δ-M's home turf) the TMS-corrected pair
    # must agree better than the uncorrected pair. NOTE: toward the forward
    # peak (small Θ) the order-dependent N–T attenuation residual — the
    # out-of-scope IMS term (Nakajima & Tanaka 1988) — grows and can invert
    # this ordering; that regime is deliberately NOT asserted here (see
    # docs/dev_notes/TMS_DISCOVERY.md, scope limits).
    base = parameters_from_yaml(joinpath(_TMS_TD,
        "test_parameters/JacobianTestFast.yaml"))
    base.brdf[1] = LambertianSurfaceScalar(0.1)
    base.scattering_params.rt_aerosols[1].aerosol =
        Aerosol(LogNormal(log(1.0), log(1.6)), 1.45, 0.001)
    base.scattering_params.rt_aerosols[1].τ_ref = 0.3

    function run_at(l_tr, sc)
        p = deepcopy(base)
        p.truncation = δBGE{p.float_type}(l_tr)
        _tms_set_numerics!(p, sc)
        m = model_from_parameters(p)
        R, _ = rt_run(m)
        return Array(R[:, 1, :])
    end

    R_off_lo = run_at(6,  NoSSCorrection())
    R_off_hi = run_at(14, NoSSCorrection())
    R_on_lo  = run_at(6,  TMSCorrection())
    R_on_hi  = run_at(14, TMSCorrection())

    gap_off = maximum(abs.(R_off_lo .- R_off_hi))
    gap_on  = maximum(abs.(R_on_lo .- R_on_hi))
    @info "TMS truncation-order convergence" gap_uncorrected = gap_off gap_tms = gap_on
    @test gap_on < gap_off
    # The correction must be genuinely active for this truncating aerosol
    # (the Rayleigh no-op test above pins the exact-cancellation case; here
    # exact ≠ truncated, so on/off must differ measurably). The magnitude
    # is geometry- and fit-dependent; assert activity, not a fixed ratio.
    activity = maximum(abs.(R_on_lo .- R_off_lo))
    @info "TMS activity |on - off| (l_max = 6)" activity
    @test activity > 1e-5
end
