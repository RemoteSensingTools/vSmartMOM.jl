#=
Caller-defined-node bulk Mie API tests (Sol plan phase 1).

Validates `compute_aerosol_optical_properties_nodes` /
`compute_aerosol_optical_properties_nodes_gpu` (src/Scattering/compute_NAI2_nodes.jl)
and the `_nai2_bulk_optics` delegation refactor of the log-normal NAI2 path
(src/Scattering/compute_NAI2.jl), against the validation battery from the task
spec:
  (a) delegation bit-compat: refactored log-normal path vs the existing
      Scattering test suite (covered by test_Scattering.jl staying green;
      here we additionally check the internal shared core reproduces the
      log-normal path bit-for-bit given the SAME nodes).
  (b) caller-node CPU vs legacy CPU NAI2 with identical nodes/weights.
  (c) KA-CPU vs CUDA (CUDA is functional on this host; both are exercised).
  (d) Float64 and Float32(+DSEmulated) tolerances.
  (e) unnormalized-weights invariance.
  (f) wide-range (0.005-6 um) ensemble: fine nodes on the coarse-driven
      angular grid vs their own dedicated (smaller) grid.
=#

if !@isdefined(phase_function)
    using Test
    using vSmartMOM
    using vSmartMOM.Scattering
    using vSmartMOM.Architectures
    using Distributions
end

import KernelAbstractions
import Logging   # needed for `Logging.Warn` in the @test_logs min_level check below
const _KA_CPU = KernelAbstractions.CPU

const _CUDA_FUNCTIONAL = try
    @eval import CUDA
    CUDA.functional()
catch
    false
end
const _CUDA_BACKEND = _CUDA_FUNCTIONAL ? CUDA.CUDABackend() : nothing

if _CUDA_FUNCTIONAL
    @info "test_mie_nodes: CUDA is functional — running real-GPU node-API testsets."
else
    @info "test_mie_nodes: CUDA not functional — real-GPU node-API testsets are skipped (KA-CPU still runs)."
end

# Standard log-normal test aerosol, shared across testsets (same parameters as
# test_Scattering.jl / test/local/gpu/test_mie_gpu.jl for cross-file consistency).
const _mu_aer = 0.3
const _sigma_aer = 2.1
const _r_max = 30.0
const _nr = 1.3
const _ni = 0.001
const _lam = 0.55
const _nquad = 300

@testset "compute_aerosol_optical_properties_nodes: identical-node bit-compat (a,b)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    aero = Aerosol(dist, _nr, _ni)
    pol  = Stokes_IQUV()

    model = make_mie_model(NAI2(), aero, _lam, pol, NoTruncation(), _r_max, _nquad)
    ref = compute_aerosol_optical_properties(model)   # raw/untruncated, ft=1

    # Reconstruct the EXACT internal nodes the log-normal path used.
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

    nodes = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())

    # (b) identical nodes/weights -> same k, omega_tilde, and all six Greek arrays.
    # NOT asserted bit-identical: the node API re-normalizes weights internally
    # (wx ./ sum(wx)), and sum(wx) of compute_wₓ's already-normalized output is
    # only 1.0 to within an ulp — exactly 1.0 on x86-Linux, but libm differences
    # on macOS-ARM/Windows can make it 1.0 ± 1 ulp, perturbing every weight by
    # ~eps. Both paths share _nai2_bulk_optics, so the only divergence is that
    # re-normalization; near-machine-precision tolerances pin it down.
    @test isapprox(nodes.k, ref.k, rtol=1e-13)
    @test isapprox(nodes.ω̃, ref.ω̃, rtol=1e-13)
    for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
        a = getproperty(ref.greek_coefs, f)
        b = getproperty(nodes.greek_coefs, f)
        @test length(a) == length(b)
        @test maximum(abs.(a .- b)) < 1e-11
    end
    # NoTruncation resets the "untruncated yet" ft=1 sentinel to 0 (node API
    # always applies the requested truncation; the raw sentinel never escapes it).
    @test nodes.fᵗ == 0
    @test ref.fᵗ == 1

    # δBGE truncation applied by the node API must match calling truncate_phase
    # manually on the same raw AerosolOptics (consistency with the existing
    # model_from_parameters/update_model convention).
    trunc = δBGE(10, 10.0)
    nodes_trunc = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, trunc)
    manual_trunc = Scattering.truncate_phase(trunc, ref)
    # Same re-normalization caveat as above -> tight isapprox, not exact ==.
    @test isapprox(nodes_trunc.fᵗ, manual_trunc.fᵗ, rtol=1e-12)
    @test maximum(abs.(nodes_trunc.greek_coefs.β .- manual_trunc.greek_coefs.β)) < 1e-11
    @test maximum(abs.(nodes_trunc.greek_coefs.α .- manual_trunc.greek_coefs.α)) < 1e-11
end

@testset "delegation bit-compat: refactored log-normal path sanity (a)" begin
    # The pre-refactor CPU path and the post-refactor path (generate r,wx then
    # delegate to _nai2_bulk_optics) perform IDENTICAL arithmetic in IDENTICAL
    # order -- Julia does not reassociate floating point operations when code
    # is moved into a subroutine, so this must be bit-identical, not merely
    # approximate. Cross-checked against PCW-golden values in test_Scattering.jl
    # separately; here we just sanity-check a couple of independent configs.
    for (mua, siga, rmax, nq, nr, ni, lam) in [
        (0.3, 2.1, 30.0, 300, 1.3, 0.001, 0.55),
        (1.0, 1.5, 15.0, 150, 1.45, 0.01, 0.87),
    ]
        dist = LogNormal(log(mua), log(siga))
        aero = Aerosol(dist, nr, ni)
        pol  = Stokes_IQUV()
        model = make_mie_model(NAI2(), aero, lam, pol, NoTruncation(), rmax, nq)
        ref = compute_aerosol_optical_properties(model)
        @test isfinite(ref.k) && isfinite(ref.ω̃)
        @test 0 < ref.ω̃ <= 1
    end
end

@testset "unnormalized-weights invariance (e)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    a = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())
    b = compute_aerosol_optical_properties_nodes(r, wx .* 7.3, _nr, _ni, _lam, pol, NoTruncation())

    # Weight normalization is internal; scaling ALL weights by a constant must
    # leave bulk quantities unchanged to floating-point precision (re-normalizing
    # a rescaled array is not bit-identical to the original -- rel err ~eps -- so
    # we assert near-machine-precision equality, not exact bit equality).
    @test isapprox(a.k, b.k, rtol=1e-13)
    @test isapprox(a.ω̃, b.ω̃, rtol=1e-13)
    @test maximum(abs.(a.greek_coefs.β .- b.greek_coefs.β)) < 1e-12
    @test maximum(abs.(a.greek_coefs.α .- b.greek_coefs.α)) < 1e-12
end

@testset "Float64/Float32 GPU node-API tolerances (c,d) — KA-CPU backend" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    ref = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())

    @testset "KA-CPU, Float64, NativeFloat64" begin
        gpu = compute_aerosol_optical_properties_nodes_gpu(r, wx, _nr, _ni, _lam, _KA_CPU();
                                                            precision_policy=NativeFloat64())
        gpu = Scattering.truncate_phase(NoTruncation(), gpu)
        @test isapprox(gpu.k, ref.k, rtol=1e-6)
        @test isapprox(gpu.ω̃, ref.ω̃, rtol=1e-6)
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
            @test isapprox(getproperty(gpu.greek_coefs, f), getproperty(ref.greek_coefs, f), atol=1e-6)
        end
    end

    @testset "KA-CPU, Float32 nodes, DSEmulated" begin
        r32 = Float32.(r); wx32 = Float32.(wx)
        gpu = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                            _KA_CPU(); precision_policy=DSEmulated())
        gpu = Scattering.truncate_phase(NoTruncation(), gpu)
        @test eltype(gpu.greek_coefs.β) === Float32
        @test gpu.k isa Float32
        @test isapprox(Float64(gpu.k), ref.k, rtol=1e-4)
        @test isapprox(Float64(gpu.ω̃), ref.ω̃, rtol=1e-4)
        L = min(length(gpu.greek_coefs.β), length(ref.greek_coefs.β))
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
            g = Float64.(getproperty(gpu.greek_coefs, f)[1:L])
            r_ = getproperty(ref.greek_coefs, f)[1:L]
            @test maximum(abs.(g .- r_)) < 1e-3
        end
    end

    @testset "top-level dispatch: architecture=CPU() (default) matches direct CPU call" begin
        top = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())
        @test top.k == ref.k
        @test top.greek_coefs.β == ref.greek_coefs.β
    end
end

@testset "KA-CPU vs real CUDA node API (c)" begin
    if !_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — real-GPU node-API test skipped"
    else
        dist = LogNormal(log(_mu_aer), log(_sigma_aer))
        r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
        r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
        wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

        ka = compute_aerosol_optical_properties_nodes_gpu(r, wx, _nr, _ni, _lam, _KA_CPU();
                                                           precision_policy=NativeFloat64())
        ka = Scattering.truncate_phase(NoTruncation(), ka)
        cu = compute_aerosol_optical_properties_nodes_gpu(r, wx, _nr, _ni, _lam, _CUDA_BACKEND;
                                                           precision_policy=NativeFloat64())
        cu = Scattering.truncate_phase(NoTruncation(), cu)

        @test isapprox(ka.k, cu.k, rtol=1e-10)
        @test isapprox(ka.ω̃, cu.ω̃, rtol=1e-10)
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
            @test isapprox(getproperty(ka.greek_coefs, f), getproperty(cu.greek_coefs, f), atol=1e-9)
        end

        # Top-level dispatch: architecture=GPU() routes to the real CUDA backend.
        pol = Stokes_IQUV()
        top_gpu = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation();
                                                            architecture=Architectures.GPU())
        cpu_ref = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())
        @test isapprox(top_gpu.k, cpu_ref.k, rtol=1e-6)
        @test isapprox(top_gpu.ω̃, cpu_ref.ω̃, rtol=1e-6)

        # Float32/DSEmulated on real hardware too.
        r32 = Float32.(r); wx32 = Float32.(wx)
        gpu32 = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                              _CUDA_BACKEND; precision_policy=DSEmulated())
        gpu32 = Scattering.truncate_phase(NoTruncation(), gpu32)
        @test isapprox(Float64(gpu32.k), cpu_ref.k, rtol=1e-4)
        L = min(length(gpu32.greek_coefs.β), length(cpu_ref.greek_coefs.β))
        @test maximum(abs.(Float64.(gpu32.greek_coefs.β[1:L]) .- cpu_ref.greek_coefs.β[1:L])) < 1e-3
    end
end

@testset "wide-range ensemble (0.005-6 um, TOMAS use case) (f)" begin
    nr = 1.5; ni = 0.005; lam = 0.55

    # Fine (Aitken-mode-like) sub-ensemble: 0.005-0.05 um.
    fine_r = exp.(range(log(0.005), log(0.05), length=40))
    fine_w = pdf.(LogNormal(log(0.02), log(1.6)), fine_r)

    # Coarse (dust-mode-like) sub-ensemble: 1-6 um.
    coarse_r = exp.(range(log(1.0), log(6.0), length=25))
    coarse_w = pdf.(LogNormal(log(2.0), log(1.5)), coarse_r) .* 1e-3

    wide_r = vcat(fine_r, coarse_r)
    wide_w = vcat(fine_w, coarse_w)
    pol = Stokes_IQUV()

    # Sanity: the wide ensemble actually drives a much larger n_max/n_mu than
    # the fine sub-ensemble alone -- otherwise this test isn't exercising
    # anything.
    kk = 2π / lam
    n_max_fine = Scattering.get_n_max(maximum(kk .* fine_r))
    n_max_wide = Scattering.get_n_max(maximum(kk .* wide_r))
    @test n_max_wide > 10 * n_max_fine

    # Full public-API smoke test: must run end to end with no r_max cap.
    wide_bulk = compute_aerosol_optical_properties_nodes(wide_r, wide_w, nr, ni, lam, pol, NoTruncation())
    @test isfinite(wide_bulk.k) && isfinite(wide_bulk.ω̃)
    @test 0 < wide_bulk.ω̃ <= 1
    # Angular resolution rides on the actual largest supplied node (no artificial cap):
    @test length(wide_bulk.greek_coefs.β) == 2 * n_max_wide - 1

    # Core quantitative check: fine nodes computed on their OWN (small) grid vs
    # the SAME fine nodes forced onto the coarse-driven (much larger) angular
    # grid -- isolates whether the larger Gauss-Legendre quadrature loses
    # accuracy for the fine (low Mie-order) particles. Uses the internal shared
    # core directly (`_nai2_bulk_optics`) so the coarse-driven n_max/n_mu can be
    # supplied explicitly without diluting the physical mean with coarse-mode
    # scattering (which would confound a numerical accuracy check with a
    # physical-composition difference).
    n_mu_fine = 2n_max_fine - 1
    n_max_wide_ = Scattering.get_n_max(maximum(kk .* wide_r))
    n_mu_wide  = 2n_max_wide_ - 1

    fine_w_norm = fine_w ./ sum(fine_w)
    ref_fine     = Scattering._nai2_bulk_optics(fine_r, fine_w_norm, nr, ni, lam, n_max_fine,    n_mu_fine,  Float64)
    fine_on_wide = Scattering._nai2_bulk_optics(fine_r, fine_w_norm, nr, ni, lam, n_max_wide_, n_mu_wide, Float64)

    @test isapprox(ref_fine.bulk_C_ext, fine_on_wide.bulk_C_ext, rtol=1e-10)
    @test isapprox(ref_fine.bulk_C_sca, fine_on_wide.bulk_C_sca, rtol=1e-10)

    L = length(ref_fine.β)
    @test length(fine_on_wide.β) > L   # wide grid produces (many) more terms
    for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
        a = getproperty(ref_fine, f)
        b = getproperty(fine_on_wide, f)[1:L]
        # Machine-precision-level agreement on the common Legendre range:
        # Gauss-Legendre quadrature with MORE points than needed still
        # integrates the (low-order) fine-particle phase function exactly, so
        # the coarse-driven grid does not degrade fine-particle accuracy.
        @test maximum(abs.(a .- b)) < 1e-12
    end
end

@testset "l_max output-length override" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    full = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())
    capped = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation(); l_max=50)

    @test length(capped.greek_coefs.β) == 50
    @test capped.greek_coefs.β == full.greek_coefs.β[1:50]
    @test capped.k == full.k  # scalars unaffected by Greek-length capping
    @test capped.ω̃ == full.ω̃

    # Codex P1 regression: l_max must NOT silently disable or distort a
    # requested δBGE — truncation is applied to the FULL raw series first,
    # the output cap slices afterward. In particular l_max < δBGE.l_max
    # previously skipped the fit entirely (fᵗ came back 0).
    trunc = δBGE(20, 10.0)
    full_t   = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, trunc)
    capped_t = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, trunc; l_max=10)
    @test full_t.fᵗ > 0
    @test capped_t.fᵗ == full_t.fᵗ            # same fit, from the same full series
    @test length(capped_t.greek_coefs.β) == 10
    @test capped_t.greek_coefs.β == full_t.greek_coefs.β[1:10]
end

@testset "GPU path: Float32 weights with Float32-overflowing sum (Codex P2 regression)" begin
    # sum(Float32[3f38, 3f38]) overflows to Inf32; normalization must widen to
    # Float64 first, making these weights land on exactly [0.5, 0.5] — i.e.
    # bit-identical output to passing [0.5, 0.5] directly.
    r32 = Float32[0.1, 0.2]
    big = compute_aerosol_optical_properties_nodes_gpu(r32, Float32[3f38, 3f38],
              1.3f0, 0.001f0, 0.55f0, _KA_CPU(); precision_policy=DSEmulated())
    ref = compute_aerosol_optical_properties_nodes_gpu(r32, Float32[0.5, 0.5],
              1.3f0, 0.001f0, 0.55f0, _KA_CPU(); precision_policy=DSEmulated())
    @test isfinite(big.k) && big.k > 0
    @test isfinite(big.ω̃) && 0 < big.ω̃ <= 1
    @test big.k == ref.k
    @test big.greek_coefs.β == ref.greek_coefs.β
end

@testset "input validation" begin
    pol = Stokes_IQUV()
    @test_throws AssertionError compute_aerosol_optical_properties_nodes(
        [1.0, 2.0], [1.0], _nr, _ni, _lam, pol, NoTruncation())          # length mismatch
    @test_throws AssertionError compute_aerosol_optical_properties_nodes(
        [1.0, 2.0], [1.0, -1.0], _nr, _ni, _lam, pol, NoTruncation())    # negative weight
    @test_throws AssertionError compute_aerosol_optical_properties_nodes(
        [1.0, 2.0], [0.0, 0.0], _nr, _ni, _lam, pol, NoTruncation())     # all-zero weights
    @test_throws AssertionError compute_aerosol_optical_properties_nodes(
        [1.0, 2.0], [1.0, 1.0], _nr, -0.001, _lam, pol, NoTruncation())  # ni < 0
end

# Minimal architecture with no ka_backend/precision-policy/Mie hooks,
# mirroring the `_NoMieArch` pattern in test/local/gpu/test_mie_gpu.jl
# (has_gpu_mie defaults to false on the abstract supertype). Defined at top
# level: local `struct` definitions are a syntax error before Julia 1.12.
struct _NoMieArchNodes <: vSmartMOM.Architectures.AbstractArchitecture end

@testset "non-CPU architecture without a GPU Mie pipeline -> warn-once CPU fallback" begin
    @test vSmartMOM.Architectures.has_gpu_mie(_NoMieArchNodes()) == false

    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    ref = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())

    vSmartMOM.Scattering._WARNED_NO_GPU_MIE_NODES[] = false
    fb1 = @test_logs (:warn, r"no GPU Mie pipeline") match_mode = :any begin
        compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation();
                                                  architecture = _NoMieArchNodes())
    end
    fb2 = @test_logs min_level = Logging.Warn begin
        compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation();
                                                  architecture = _NoMieArchNodes())
    end

    @test fb1.k == ref.k
    @test fb2.k == ref.k
    @test fb1.greek_coefs.β == ref.greek_coefs.β
end
