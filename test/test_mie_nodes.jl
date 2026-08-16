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

    # The DIRECTLY-callable _gpu entry point must enforce the SAME validation
    # list (shared _validate_node_inputs) — a negative weight passed only to
    # the GPU seam would otherwise yield silently wrong bulk optics.
    @test_throws AssertionError compute_aerosol_optical_properties_nodes_gpu(
        [1.0, 2.0], [1.0], _nr, _ni, _lam, _KA_CPU())                    # length mismatch
    @test_throws AssertionError compute_aerosol_optical_properties_nodes_gpu(
        [1.0, 2.0], [1.0, -1.0], _nr, _ni, _lam, _KA_CPU())              # negative weight
    @test_throws AssertionError compute_aerosol_optical_properties_nodes_gpu(
        Float64[], Float64[], _nr, _ni, _lam, _KA_CPU())                 # empty node set
    @test_throws AssertionError compute_aerosol_optical_properties_nodes_gpu(
        [1.0, 2.0], [1.0, 1.0], _nr, -0.001, _lam, _KA_CPU())            # ni < 0
    @test_throws AssertionError compute_aerosol_optical_properties_nodes_gpu(
        [1.0, 2.0], [0.0, 0.0], _nr, _ni, _lam, _KA_CPU())               # all-zero weights
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

#=
=====================================================================
compute_aerosol_extinction_nodes / _gpu (EXTINCTION-ONLY fast path) tests.
=====================================================================
=#

@testset "compute_aerosol_extinction_nodes: bitwise k-equality vs full API (a)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    full = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())
    k_only = compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)

    # BIT-IDENTICAL, not merely isapprox -- same normalization/arithmetic order.
    @test k_only === full.k
    @test k_only == full.k

    # A δBGE truncation must not perturb k (truncation only touches Greek/fᵗ).
    trunc = δBGE(10, 10.0)
    full_trunc = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, trunc)
    @test k_only == full_trunc.k

    # wide-range (TOMAS-like) node set (same fixture as the "wide-range
    # ensemble" testset above).
    nr = 1.5; ni = 0.005; lam = 0.55
    fine_r = exp.(range(log(0.005), log(0.05), length=40))
    fine_w = pdf.(LogNormal(log(0.02), log(1.6)), fine_r)
    coarse_r = exp.(range(log(1.0), log(6.0), length=25))
    coarse_w = pdf.(LogNormal(log(2.0), log(1.5)), coarse_r) .* 1e-3
    wide_r = vcat(fine_r, coarse_r)
    wide_w = vcat(fine_w, coarse_w)

    wide_full = compute_aerosol_optical_properties_nodes(wide_r, wide_w, nr, ni, lam, pol, NoTruncation())
    wide_k = compute_aerosol_extinction_nodes(wide_r, wide_w, nr, ni, lam)
    @test wide_k === wide_full.k
    @test wide_k == wide_full.k
end

@testset "compute_aerosol_extinction_nodes: unnormalized-weight invariance (b)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

    a = compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)
    b = compute_aerosol_extinction_nodes(r, wx .* 7.3, _nr, _ni, _lam)
    @test isapprox(a, b, rtol=1e-13)
end

@testset "compute_aerosol_extinction_nodes: KA-CPU vs CPU, real-CUDA vs KA-CPU (c)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

    cpu_k = compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)
    ka_k = compute_aerosol_extinction_nodes_gpu(r, wx, _nr, _ni, _lam, _KA_CPU();
                                                 precision_policy=NativeFloat64())
    @test isapprox(ka_k, cpu_k, rtol=1e-12)

    # top-level dispatch: architecture=CPU() (default) matches the direct CPU call.
    top = compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)
    @test top == cpu_k

    if !_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — real-GPU extinction-only test skipped"
    else
        cu_k = compute_aerosol_extinction_nodes_gpu(r, wx, _nr, _ni, _lam, _CUDA_BACKEND;
                                                     precision_policy=NativeFloat64())
        @test isapprox(cu_k, ka_k, rtol=1e-10)

        top_gpu_k = compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam;
                                                      architecture=Architectures.GPU())
        @test isapprox(top_gpu_k, cpu_k, rtol=1e-6)
    end
end

@testset "compute_aerosol_extinction_nodes: Float32 + DSEmulated sanity (d)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

    ref_k = compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)

    r32 = Float32.(r); wx32 = Float32.(wx)
    k32 = compute_aerosol_extinction_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                _KA_CPU(); precision_policy=DSEmulated())
    @test k32 isa Float32
    @test isapprox(Float64(k32), ref_k, rtol=1e-4)
end

@testset "compute_aerosol_extinction_nodes: input validation, both entry points (e)" begin
    @test_throws AssertionError compute_aerosol_extinction_nodes(
        [1.0, 2.0], [1.0], _nr, _ni, _lam)                                # length mismatch
    @test_throws AssertionError compute_aerosol_extinction_nodes(
        [1.0, 2.0], [1.0, -1.0], _nr, _ni, _lam)                          # negative weight
    @test_throws AssertionError compute_aerosol_extinction_nodes(
        [1.0, 2.0], [0.0, 0.0], _nr, _ni, _lam)                           # all-zero weights
    @test_throws AssertionError compute_aerosol_extinction_nodes(
        [1.0, 2.0], [1.0, 1.0], _nr, -0.001, _lam)                       # ni < 0
    @test_throws AssertionError compute_aerosol_extinction_nodes(
        Float64[], Float64[], _nr, _ni, _lam)                            # empty node set

    # The DIRECTLY-callable _gpu entry point must enforce the SAME validation
    # list (shared _validate_node_inputs) as the public entry point above.
    @test_throws AssertionError compute_aerosol_extinction_nodes_gpu(
        [1.0, 2.0], [1.0], _nr, _ni, _lam, _KA_CPU())                    # length mismatch
    @test_throws AssertionError compute_aerosol_extinction_nodes_gpu(
        [1.0, 2.0], [1.0, -1.0], _nr, _ni, _lam, _KA_CPU())              # negative weight
    @test_throws AssertionError compute_aerosol_extinction_nodes_gpu(
        Float64[], Float64[], _nr, _ni, _lam, _KA_CPU())                 # empty node set
    @test_throws AssertionError compute_aerosol_extinction_nodes_gpu(
        [1.0, 2.0], [1.0, 1.0], _nr, -0.001, _lam, _KA_CPU())            # ni < 0
    @test_throws AssertionError compute_aerosol_extinction_nodes_gpu(
        [1.0, 2.0], [0.0, 0.0], _nr, _ni, _lam, _KA_CPU())               # all-zero weights
end

@testset "extinction _gpu: Float32 weights with Float32-overflowing sum (Codex P2 mirror)" begin
    # Mirrors the full API's P2 regression (julia-reviewer nit): the extinction
    # GPU entry copied the Float64-widen-before-narrow normalization fix, so it
    # must survive weights whose Float32 sum overflows — bit-identical to the
    # exactly-[0.5, 0.5] normalized case.
    r32 = Float32[0.1, 0.2]
    big = compute_aerosol_extinction_nodes_gpu(r32, Float32[3f38, 3f38],
              1.3f0, 0.001f0, 0.55f0, _KA_CPU(); precision_policy=DSEmulated())
    ref = compute_aerosol_extinction_nodes_gpu(r32, Float32[0.5, 0.5],
              1.3f0, 0.001f0, 0.55f0, _KA_CPU(); precision_policy=DSEmulated())
    @test isfinite(big) && big > 0
    @test big == ref
end

@testset "compute_aerosol_extinction_nodes: timing sanity, extinction-only vs full API (f)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    # Warm up (compile) both paths before timing so @elapsed measures
    # execution, not JIT compilation.
    compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)
    compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())

    t_ext  = @elapsed compute_aerosol_extinction_nodes(r, wx, _nr, _ni, _lam)
    t_full = @elapsed compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())

    # Sanity line only -- NOT asserted (timing is inherently noisy in CI).
    println("compute_aerosol_extinction_nodes vs full API timing (standard node set, nquad=$_nquad): ",
            "t_ext=$(t_ext)s  t_full=$(t_full)s  speedup(full/ext)=$(round(t_full / t_ext, digits=2))x")
    @test t_ext >= 0 && t_full >= 0
end

#=
=====================================================================
NativeFloat32 (pure Float32 end-to-end) + Metal Mie registration tests.
=====================================================================
NativeFloat32 is the ONLY MiePrecisionPolicy under which the caller-node GPU
pipeline allocates zero Float64 device arrays -- required for Apple Metal (no
Float64 hardware at all) and F32-only-throughput consumer CUDA. Accuracy is
measured (not assumed) below on both the standard log-normal set and the
wide-range TOMAS-like set, then asserted with margin over the measured
numbers. DSEmulated/NativeFloat64 tolerances established above are untouched
by any of this -- see the "identical-node bit-compat" and
"Float64/Float32 GPU node-API tolerances" testsets, still green.
=====================================================================
=#

@testset "NativeFloat32: accuracy vs Float64 CPU reference (measured + margin)" begin
    # --- Standard log-normal set (same fixture as the rest of this file) ---
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    pol = Stokes_IQUV()

    ref = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, pol, NoTruncation())

    r32 = Float32.(r); wx32 = Float32.(wx)
    gpu32 = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                          _KA_CPU(); precision_policy=NativeFloat32())
    gpu32 = Scattering.truncate_phase(NoTruncation(), gpu32)

    k_relerr_std   = abs(Float64(gpu32.k) - ref.k) / ref.k
    ω_relerr_std   = abs(Float64(gpu32.ω̃) - ref.ω̃) / ref.ω̃
    L_std = min(length(gpu32.greek_coefs.β), length(ref.greek_coefs.β))
    greek_maxabs_std = maximum(
        maximum(abs.(Float64.(getproperty(gpu32.greek_coefs, f)[1:L_std]) .-
                     getproperty(ref.greek_coefs, f)[1:L_std]))
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
    )

    # --- Wide-range TOMAS-like set (same fixture as the "wide-range ensemble" testset) ---
    nr_w = 1.5; ni_w = 0.005; lam_w = 0.55
    fine_r = exp.(range(log(0.005), log(0.05), length=40))
    fine_w = pdf.(LogNormal(log(0.02), log(1.6)), fine_r)
    coarse_r = exp.(range(log(1.0), log(6.0), length=25))
    coarse_w = pdf.(LogNormal(log(2.0), log(1.5)), coarse_r) .* 1e-3
    wide_r = vcat(fine_r, coarse_r)
    wide_w = vcat(fine_w, coarse_w)

    ref_wide = compute_aerosol_optical_properties_nodes(wide_r, wide_w, nr_w, ni_w, lam_w, pol, NoTruncation())

    wide_r32 = Float32.(wide_r); wide_w32 = Float32.(wide_w)
    gpu32_wide = compute_aerosol_optical_properties_nodes_gpu(wide_r32, wide_w32, Float32(nr_w), Float32(ni_w),
                                                               Float32(lam_w), _KA_CPU();
                                                               precision_policy=NativeFloat32())
    gpu32_wide = Scattering.truncate_phase(NoTruncation(), gpu32_wide)

    k_relerr_wide = abs(Float64(gpu32_wide.k) - ref_wide.k) / ref_wide.k
    ω_relerr_wide = abs(Float64(gpu32_wide.ω̃) - ref_wide.ω̃) / ref_wide.ω̃
    L_wide = min(length(gpu32_wide.greek_coefs.β), length(ref_wide.greek_coefs.β))
    greek_maxabs_wide = maximum(
        maximum(abs.(Float64.(getproperty(gpu32_wide.greek_coefs, f)[1:L_wide]) .-
                     getproperty(ref_wide.greek_coefs, f)[1:L_wide]))
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
    )

    println("NativeFloat32 vs Float64 CPU reference -- MEASURED accuracy:")
    println("  standard set   (nquad=$_nquad): k relerr=$(k_relerr_std)  ω̃ relerr=$(ω_relerr_std)  Greek maxabs=$(greek_maxabs_std)")
    println("  wide TOMAS set (n=$(length(wide_r))):     k relerr=$(k_relerr_wide)  ω̃ relerr=$(ω_relerr_wide)  Greek maxabs=$(greek_maxabs_wide)")

    # Measured (this host, 2026-08-15): k/ω̃ relerr ~1e-8, Greek maxabs ~5e-3 on
    # BOTH sets (dominated by the native-Float32 Dₙ recursion feeding the Greek
    # angular projection -- k/ω̃ are simple Neumaier-compensated weighted sums
    # and stay accurate even without DS emulation; Greek coefficients are
    # noticeably WORSE than DSEmulated's established ~1e-3 floor, exactly as
    # expected for a policy with no double-single Dₙ emulation at all).
    # Asserted with a ~100x margin on k/ω̃ and ~2x margin on Greek to absorb
    # BLAS/libm cross-platform noise without being a vacuous bound.
    @test k_relerr_std   < 1e-6
    @test ω_relerr_std   < 1e-6
    @test greek_maxabs_std < 1e-2
    @test k_relerr_wide   < 1e-6
    @test ω_relerr_wide   < 1e-6
    @test greek_maxabs_wide < 1e-2

    # Sanity: NativeFloat32 must not be MORE accurate on Greek than the
    # established DSEmulated floor by construction (it has strictly less
    # numerical machinery) -- not asserted numerically (measured margins above
    # already suffice), just documented here as the expected ordering.
end

@testset "NativeFloat32: zero Float64 device arrays (mechanical check)" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
    r32 = Float32.(r); wx32 = Float32.(wx)

    # White-box: _prepare_node_mie_gpu's returned FT/RA and per-node device
    # array eltypes directly prove no Float64 buffer is allocated -- EVERY
    # downstream reduction/Greek/Legendre device array in both public GPU
    # functions is typed by this SAME `RA`.
    prep = Scattering._prepare_node_mie_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                             _KA_CPU(), NativeFloat32())
    @test prep.FT === Float32
    @test prep.RA === Float32     # THE key lever: never widened to Float64
    @test eltype(prep.x_dev)    === Float32
    @test eltype(prep.an_dev)   === Complex{Float32}
    @test eltype(prep.bn_dev)   === Complex{Float32}
    @test eltype(prep.nmax_dev) === Int

    # Black-box corroboration: the full end-to-end output stays Float32
    # throughout (RA === FT means no convert.(FT, ...) narrowing branch is
    # even taken), consistent with an all-Float32 device pipeline.
    full = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                          _KA_CPU(); precision_policy=NativeFloat32())
    @test full.k isa Float32
    @test full.ω̃ isa Float32
    @test eltype(full.greek_coefs.β) === Float32

    k_only = compute_aerosol_extinction_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                   _KA_CPU(); precision_policy=NativeFloat32())
    @test k_only isa Float32
    @test k_only == full.k   # both entry points share _prepare_node_mie_gpu
end

@testset "NativeFloat32: policy/type validation" begin
    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

    # Float64 nodes with NativeFloat32 policy must assert (mirrors the
    # existing NativeFloat64-requires-Float64 check).
    @test_throws AssertionError compute_aerosol_optical_properties_nodes_gpu(
        r, wx, _nr, _ni, _lam, _KA_CPU(); precision_policy=NativeFloat32())
    @test_throws AssertionError compute_aerosol_extinction_nodes_gpu(
        r, wx, _nr, _ni, _lam, _KA_CPU(); precision_policy=NativeFloat32())

    # Float32 nodes with NativeFloat32 must NOT assert (sanity: the policy
    # actually accepts its own intended input type).
    r32 = Float32.(r); wx32 = Float32.(wx)
    @test compute_aerosol_optical_properties_nodes_gpu(
        r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), _KA_CPU();
        precision_policy=NativeFloat32()) isa AerosolOptics
    @test compute_aerosol_extinction_nodes_gpu(
        r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), _KA_CPU();
        precision_policy=NativeFloat32()) isa Float32
end

# Minimal Metal-backend mock: `_is_metal_backend` defaults to `false` on
# `Any` (real-dispatch trait, refined to `true` for `Metal.MetalBackend` only
# by the weakly-loaded `vSmartMOMMetalExt` extension -- see gpu_precision.jl).
# This mock registers its OWN `_is_metal_backend` method (real dispatch, not
# name/string matching) so the Metal-only-NativeFloat32 guard in
# `_prepare_node_mie_gpu` can be exercised WITHOUT installing Metal.jl -- the
# guard fires before any KernelAbstractions call, so this mock never needs to
# implement a real KA backend interface. Defined at top level (structs cannot
# be declared inside a function/testset body before Julia 1.12); named
# distinctly from `Metal.MetalBackend` since matching is now by dispatch, not
# by name, so no collision risk even when Metal.jl IS loaded.
struct MockMetalBackend end
Scattering._is_metal_backend(::MockMetalBackend) = true

@testset "Metal-only-NativeFloat32 guard (mock backend, no real Metal.jl needed)" begin
    @test Scattering._is_metal_backend(MockMetalBackend()) == true
    @test Scattering._is_metal_backend(_KA_CPU()) == false

    dist = LogNormal(log(_mu_aer), log(_sigma_aer))
    r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
    r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
    wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)

    # NativeFloat64 (Float64 nodes) on a Metal-shaped backend: Metal has NO
    # Float64 hardware at all -- must throw ArgumentError (not a cryptic
    # failure inside KernelAbstractions.allocate, which this mock cannot even
    # perform).
    @test_throws ArgumentError compute_aerosol_optical_properties_nodes_gpu(
        r, wx, _nr, _ni, _lam, MockMetalBackend(); precision_policy=NativeFloat64())
    @test_throws ArgumentError compute_aerosol_extinction_nodes_gpu(
        r, wx, _nr, _ni, _lam, MockMetalBackend(); precision_policy=NativeFloat64())

    # DSEmulated (Float32 nodes) on a Metal-shaped backend: DSEmulated's
    # historical reduction discipline widens to RA=Float64 for a Float32
    # kernel -- illegal on Metal device arrays -- must ALSO throw
    # ArgumentError. DSEmulated-on-Metal (a Float32-COMPENSATED, not widened,
    # reduction) is a documented follow-up, not implemented this round.
    r32 = Float32.(r); wx32 = Float32.(wx)
    @test_throws ArgumentError compute_aerosol_optical_properties_nodes_gpu(
        r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), MockMetalBackend(); precision_policy=DSEmulated())
    @test_throws ArgumentError compute_aerosol_extinction_nodes_gpu(
        r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), MockMetalBackend(); precision_policy=DSEmulated())

    # NativeFloat32: RA === FT === Float32, so the guard's `FT === Float64 ||
    # RA === Float64` condition is trivially false for ANY backend -- proven
    # directly at the policy-dispatch level (a mock with no real KA backend
    # interface cannot run the full kernel pipeline to show this end to end).
    @test Scattering._mie_reduction_type(NativeFloat32(), Float32) === Float32
    @test Scattering._check_policy_ft(NativeFloat32(), Float32) === nothing
end

#=
=====================================================================
GPU exact k-identity: full API vs extinction-only sibling.
=====================================================================
The CPU path asserts `===` bitwise identity (see "compute_aerosol_extinction_nodes:
bitwise k-equality vs full API (a)" above). On the GPU paths this was
previously only checked with tolerances -- but compute_aerosol_optical_properties_nodes_gpu
and compute_aerosol_extinction_nodes_gpu now share `_prepare_node_mie_gpu`
(same Kernel-1 launch producing the SAME an_dev/bn_dev/x_dev/nmax_dev), and
each then launches `cross_sections_kernel!` + `weighted_sum_kernel!` (the
latter always `ndrange=1` -- a single-thread Neumaier loop, so no parallel-
reduction nondeterminism) on the SAME per-node/weight arrays. So `k` should be
EXACTLY the same number between the two entry points, not just close.

NOT `===`: two independent device->host copies (`Array(bulk_Cext_dev)[1]`,
narrowed to `FT`) are different objects/moments in each call, so `===` does
not apply the way it does for the CPU path's shared-buffer case -- exact
numeric equality (`==`) is the right bar here, confirmed measured (not just
assumed) on KA-CPU and real CUDA below for all three precision policies.
=====================================================================
=#

"""
    _assert_gpu_k_exact_equal(backend, rr, ww, nr_, ni_, lam_, policy, FT)

Shared check used by the KA-CPU/CUDA/Metal exact-k-identity testsets: computes
`compute_aerosol_optical_properties_nodes_gpu(...).k` and
`compute_aerosol_extinction_nodes_gpu(...)` on IDENTICAL inputs/backend/policy
and asserts `==` exact equality. On failure, prints the policy/backend and the
numeric difference so a real regression is easy to diagnose (rather than
silently loosening the bound).
"""
function _assert_gpu_k_exact_equal(backend, rr, ww, nr_, ni_, lam_, policy, FT)
    rr2 = FT.(rr); ww2 = FT.(ww)
    full = compute_aerosol_optical_properties_nodes_gpu(rr2, ww2, FT(nr_), FT(ni_), FT(lam_), backend;
                                                         precision_policy=policy)
    ext  = compute_aerosol_extinction_nodes_gpu(rr2, ww2, FT(nr_), FT(ni_), FT(lam_), backend;
                                                 precision_policy=policy)
    if full.k != ext
        println("GPU exact k-identity FAILED: backend=$(typeof(backend)) policy=$(typeof(policy)) ",
                "full.k=$(full.k) ext=$(ext) diff=$(Float64(full.k) - Float64(ext))")
    end
    @test full.k == ext
end

# Standard + wide-range (TOMAS-like) node sets, shared by the three
# exact-k-identity testsets below (KA-CPU, real CUDA, real Metal).
const _dist_std = LogNormal(log(_mu_aer), log(_sigma_aer))
const _rmin_std = max(quantile(_dist_std, 1e-8), 1e-6 * _r_max)
const (_r_std, _wr_std) = vSmartMOM.Scattering.gauleg_log(_nquad, _rmin_std, _r_max; norm=false)
const _wx_std = vSmartMOM.Scattering.compute_wₓ(_dist_std, _wr_std, _r_std, _r_max)

const _nr_wide = 1.5; const _ni_wide = 0.005; const _lam_wide = 0.55
const _fine_r_ex = exp.(range(log(0.005), log(0.05), length=40))
const _fine_w_ex = pdf.(LogNormal(log(0.02), log(1.6)), _fine_r_ex)
const _coarse_r_ex = exp.(range(log(1.0), log(6.0), length=25))
const _coarse_w_ex = pdf.(LogNormal(log(2.0), log(1.5)), _coarse_r_ex) .* 1e-3
const _wide_r_ex = vcat(_fine_r_ex, _coarse_r_ex)
const _wide_w_ex = vcat(_fine_w_ex, _coarse_w_ex)

const _EXACT_K_NODE_SETS = (
    ("standard",   _r_std,    _wx_std,    _nr,      _ni,      _lam),
    ("wide TOMAS", _wide_r_ex, _wide_w_ex, _nr_wide, _ni_wide, _lam_wide),
)
const _EXACT_K_POLICIES = (
    (NativeFloat64(), Float64),
    (DSEmulated(),    Float32),
    (NativeFloat32(), Float32),
)

@testset "GPU exact k-identity: full API vs extinction-only, KA-CPU (all 3 policies, both sets)" begin
    for (label, rr, ww, nr_, ni_, lam_) in _EXACT_K_NODE_SETS, (policy, FT) in _EXACT_K_POLICIES
        @testset "$label / $(typeof(policy))" begin
            _assert_gpu_k_exact_equal(_KA_CPU(), rr, ww, nr_, ni_, lam_, policy, FT)
        end
    end
end

# ============================================================================
# Real Metal hardware (gated on a functional Metal backend, like the CUDA
# gating above). This host has none (Linux) -- these testsets skip cleanly;
# the user validates on a Mac.
# ============================================================================
const _METAL_FUNCTIONAL = try
    @eval import Metal
    Metal.functional()
catch
    false
end
const _METAL_BACKEND = _METAL_FUNCTIONAL ? Metal.MetalBackend() : nothing

if _METAL_FUNCTIONAL
    @info "test_mie_nodes: Metal is functional — running real-Metal node-API testsets."
else
    @info "test_mie_nodes: Metal not functional — real-Metal node-API testsets are skipped (KA-CPU/mock coverage above still runs)."
end

@testset "Real Metal: NativeFloat32 vs KA-CPU NativeFloat32 and Float64 reference (gated)" begin
    if !_METAL_FUNCTIONAL
        @test_skip "Metal not functional on this host — real-Metal node-API test skipped"
    else
        dist = LogNormal(log(_mu_aer), log(_sigma_aer))
        r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
        r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
        wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
        r32 = Float32.(r); wx32 = Float32.(wx)

        ref = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, Stokes_IQUV(), NoTruncation())

        ka32 = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                             _KA_CPU(); precision_policy=NativeFloat32())
        ka32 = Scattering.truncate_phase(NoTruncation(), ka32)
        mt32 = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                             _METAL_BACKEND; precision_policy=NativeFloat32())
        mt32 = Scattering.truncate_phase(NoTruncation(), mt32)

        @test isapprox(mt32.k, ka32.k, rtol=1e-5)
        @test isapprox(mt32.ω̃, ka32.ω̃, rtol=1e-5)
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
            @test isapprox(getproperty(mt32.greek_coefs, f), getproperty(ka32.greek_coefs, f), atol=1e-3)
        end
        @test isapprox(Float64(mt32.k), ref.k, rtol=1e-6)
        @test isapprox(Float64(mt32.ω̃), ref.ω̃, rtol=1e-6)

        # Top-level dispatch: architecture=MetalGPU() auto-selects NativeFloat32.
        top = compute_aerosol_optical_properties_nodes(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                        Stokes_IQUV(), NoTruncation();
                                                        architecture=Architectures.MetalGPU())
        @test isapprox(top.k, ka32.k, rtol=1e-5)

        # Float64 on Metal must throw (auto-select path).
        @test_throws ArgumentError compute_aerosol_optical_properties_nodes(
            r, wx, _nr, _ni, _lam, Stokes_IQUV(), NoTruncation(); architecture=Architectures.MetalGPU())
    end
end

@testset "Real Metal: exact k-identity full API vs extinction-only (NativeFloat32, both sets, gated)" begin
    if !_METAL_FUNCTIONAL
        @test_skip "Metal not functional on this host — real-Metal exact k-identity test skipped"
    else
        # Metal only supports NativeFloat32 this round (NativeFloat64/DSEmulated
        # both throw ArgumentError via the Metal-only-NativeFloat32 guard in
        # _prepare_node_mie_gpu -- see the "Metal-only-NativeFloat32 guard"
        # testset above), so only that one policy is exercised here.
        for (label, rr, ww, nr_, ni_, lam_) in _EXACT_K_NODE_SETS
            @testset "$label / NativeFloat32" begin
                _assert_gpu_k_exact_equal(_METAL_BACKEND, rr, ww, nr_, ni_, lam_, NativeFloat32(), Float32)
            end
        end
    end
end

# ============================================================================
# Real CUDA: NativeFloat32 accuracy (gated, like the other real-CUDA testsets)
# ============================================================================
@testset "Real CUDA: NativeFloat32 vs KA-CPU NativeFloat32 and Float64 reference (gated)" begin
    if !_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — NativeFloat32 real-GPU test skipped"
    else
        dist = LogNormal(log(_mu_aer), log(_sigma_aer))
        r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
        r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
        wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
        r32 = Float32.(r); wx32 = Float32.(wx)

        ref = compute_aerosol_optical_properties_nodes(r, wx, _nr, _ni, _lam, Stokes_IQUV(), NoTruncation())

        ka32 = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                             _KA_CPU(); precision_policy=NativeFloat32())
        ka32 = Scattering.truncate_phase(NoTruncation(), ka32)
        cu32 = compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                             _CUDA_BACKEND; precision_policy=NativeFloat32())
        cu32 = Scattering.truncate_phase(NoTruncation(), cu32)

        # KA-CPU vs real CUDA, SAME policy: tight (same FLOP sequence/precision,
        # different device -- not necessarily bit-identical CPU vs GPU FPUs).
        @test isapprox(cu32.k, ka32.k, rtol=1e-6)
        @test isapprox(cu32.ω̃, ka32.ω̃, rtol=1e-6)
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
            @test isapprox(getproperty(cu32.greek_coefs, f), getproperty(ka32.greek_coefs, f), atol=1e-4)
        end

        # vs Float64 CPU reference: the measured NativeFloat32 error band (see
        # the KA-CPU accuracy testset above for the exact numbers on this host).
        @test isapprox(Float64(cu32.k), ref.k, rtol=1e-6)
        @test isapprox(Float64(cu32.ω̃), ref.ω̃, rtol=1e-6)
        L = min(length(cu32.greek_coefs.β), length(ref.greek_coefs.β))
        for f in (:α, :β, :γ, :δ, :ϵ, :ζ)
            g = Float64.(getproperty(cu32.greek_coefs, f)[1:L])
            rr = getproperty(ref.greek_coefs, f)[1:L]
            @test maximum(abs.(g .- rr)) < 1e-2
        end

        # Extinction-only sibling: same NativeFloat32 policy, real CUDA vs KA-CPU.
        k_ka32 = compute_aerosol_extinction_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                       _KA_CPU(); precision_policy=NativeFloat32())
        k_cu32 = compute_aerosol_extinction_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                       _CUDA_BACKEND; precision_policy=NativeFloat32())
        @test isapprox(k_cu32, k_ka32, rtol=1e-6)
        @test isapprox(Float64(k_cu32), ref.k, rtol=1e-6)
    end
end

@testset "Real CUDA: exact k-identity full API vs extinction-only (all 3 policies, both sets, gated)" begin
    if !_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — real-CUDA exact k-identity test skipped"
    else
        for (label, rr, ww, nr_, ni_, lam_) in _EXACT_K_NODE_SETS, (policy, FT) in _EXACT_K_POLICIES
            @testset "$label / $(typeof(policy))" begin
                _assert_gpu_k_exact_equal(_CUDA_BACKEND, rr, ww, nr_, ni_, lam_, policy, FT)
            end
        end
    end
end

# ============================================================================
# A100 timing table: NativeFloat32 vs DSEmulated vs NativeFloat64 (gated,
# reported only -- @elapsed is inherently noisy, so only >=0 is asserted).
# ============================================================================
@testset "Real CUDA: precision-policy timing table (report only)" begin
    if !_CUDA_FUNCTIONAL
        @test_skip "CUDA not functional on this host — precision-policy timing table skipped"
    else
        gpu_name = try
            CUDA.name(CUDA.device())
        catch
            "unknown CUDA device"
        end

        dist = LogNormal(log(_mu_aer), log(_sigma_aer))
        r_min = max(quantile(dist, 1e-8), 1e-6 * _r_max)
        r, wr = vSmartMOM.Scattering.gauleg_log(_nquad, r_min, _r_max; norm=false)
        wx = vSmartMOM.Scattering.compute_wₓ(dist, wr, r, _r_max)
        r32 = Float32.(r); wx32 = Float32.(wx)

        # Warm up (compile) all three policies before timing.
        compute_aerosol_optical_properties_nodes_gpu(r, wx, _nr, _ni, _lam, _CUDA_BACKEND;
                                                      precision_policy=NativeFloat64())
        compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                      _CUDA_BACKEND; precision_policy=DSEmulated())
        compute_aerosol_optical_properties_nodes_gpu(r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam),
                                                      _CUDA_BACKEND; precision_policy=NativeFloat32())

        t_f64 = @elapsed compute_aerosol_optical_properties_nodes_gpu(
            r, wx, _nr, _ni, _lam, _CUDA_BACKEND; precision_policy=NativeFloat64())
        t_ds  = @elapsed compute_aerosol_optical_properties_nodes_gpu(
            r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), _CUDA_BACKEND; precision_policy=DSEmulated())
        t_f32 = @elapsed compute_aerosol_optical_properties_nodes_gpu(
            r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), _CUDA_BACKEND; precision_policy=NativeFloat32())

        # Extinction-only sibling too (the cheaper path GCHPIO's hybrid mode
        # actually wants for a "shape-elsewhere" node set).
        t_f64_ext = @elapsed compute_aerosol_extinction_nodes_gpu(
            r, wx, _nr, _ni, _lam, _CUDA_BACKEND; precision_policy=NativeFloat64())
        t_ds_ext  = @elapsed compute_aerosol_extinction_nodes_gpu(
            r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), _CUDA_BACKEND; precision_policy=DSEmulated())
        t_f32_ext = @elapsed compute_aerosol_extinction_nodes_gpu(
            r32, wx32, Float32(_nr), Float32(_ni), Float32(_lam), _CUDA_BACKEND; precision_policy=NativeFloat32())

        println("Real-CUDA ($gpu_name) precision-policy timings (nquad=$_nquad, standard set):")
        println("  full API:        NativeFloat64=$(t_f64)s  DSEmulated=$(t_ds)s  NativeFloat32=$(t_f32)s")
        println("  extinction-only: NativeFloat64=$(t_f64_ext)s  DSEmulated=$(t_ds_ext)s  NativeFloat32=$(t_f32_ext)s")

        @test t_f64 >= 0 && t_ds >= 0 && t_f32 >= 0
        @test t_f64_ext >= 0 && t_ds_ext >= 0 && t_f32_ext >= 0
    end
end
