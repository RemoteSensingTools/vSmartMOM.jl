#=
GPU Mie scattering tests -- accuracy and performance evaluation.

Tests three precision strategies against Float64 CPU reference:
  1. CPU baseline (existing code)
  2. GPU/CPU-backend with NativeFloat64 Dn recursion
  3. GPU/CPU-backend with DSEmulated Dn recursion (Float32 DoubleSingle pairs)

Accuracy metrics:
  - DS arithmetic unit tests (vs Float64)
  - Neumaier summation unit tests
  - Dn recursion accuracy for challenging size parameters
  - Full NAI2 pipeline: Greek coefficients, SSA, extinction

Performance metrics:
  - Wall-clock timing for each scheme at multiple nquad_radius values
=#

using Test
using vSmartMOM
using Distributions
using LinearAlgebra
using Printf
using Logging

# KernelAbstractions is re-exported through vSmartMOM's Scattering module
using vSmartMOM.Scattering
using vSmartMOM.Architectures: CPU, GPU
using vSmartMOM.Scattering: DoubleSingle, ComplexDS, NeumaierAccum,
    TwoSum, TwoProd, ds_add, ds_sub, ds_mul, ds_div, ds_inv,
    cds_add, cds_sub, cds_mul, cds_div, cds_inv, cds_from_real, cds_mul_real, cds_complex, to_complex,
    neumaier_add, neumaier_sum, ComplexNeumaier, cneumaier_add, cneumaier_sum,
    NativeFloat64, DSEmulated,
    compute_aerosol_optical_properties_gpu, compute_aerosol_optical_properties,
    get_n_max, compute_mie_ab!, MieModel, NAI2, Aerosol, AerosolOptics,
    GreekCoefs, Stokes_IQUV, δBGE, make_mie_model

import KernelAbstractions
const KA_CPU = KernelAbstractions.CPU

# Real-GPU tests are gated on a functional CUDA device. On CPU-only machines
# `CUDA_FUNCTIONAL` stays false and every GPU testset below is skipped, so the
# file passes everywhere.
const CUDA_FUNCTIONAL = try
    @eval import CUDA
    CUDA.functional()
catch
    false
end
const CUDA_BACKEND = CUDA_FUNCTIONAL ? CUDA.CUDABackend() : nothing

if CUDA_FUNCTIONAL
    @info "test_mie_gpu: CUDA is functional — running real-GPU Mie testsets."
else
    @info "test_mie_gpu: CUDA not functional — skipping real-GPU Mie testsets."
end

# ============================================================================
# Unit Tests: DoubleSingle Arithmetic
# ============================================================================
@testset "DoubleSingle arithmetic" begin

    @testset "TwoSum error-free" begin
        # TwoSum should be exact: a + b = hi + lo
        for _ in 1:10000
            a = randn(Float32)
            b = randn(Float32)
            ds = TwoSum(a, b)
            # Check exactness in Float64
            exact = Float64(a) + Float64(b)
            reconstructed = Float64(ds.hi) + Float64(ds.lo)
            @test reconstructed == exact
        end
    end

    @testset "TwoProd error-free" begin
        for _ in 1:10000
            a = randn(Float32)
            b = randn(Float32)
            ds = TwoProd(a, b)
            exact = Float64(a) * Float64(b)
            reconstructed = Float64(ds.hi) + Float64(ds.lo)
            @test reconstructed == exact
        end
    end

    @testset "DS add/sub/mul/div accuracy" begin
        max_rel_err_add = 0.0
        max_rel_err_mul = 0.0
        max_rel_err_div = 0.0

        for _ in 1:100000
            a64 = randn() * 1000
            b64 = randn() * 1000

            a_ds = DoubleSingle{Float32}(a64)
            b_ds = DoubleSingle{Float32}(b64)

            # Add. Guard against catastrophic cancellation (a ≈ -b): the
            # RELATIVE error of any finite-precision add is unbounded there
            # (tiny exact sum in the denominator), so unlucky draws made this
            # stochastic bound flaky. DS accuracy under cancellation is
            # covered deterministically by the Dn-recursion testset below.
            sum_ds = ds_add(a_ds, b_ds)
            sum_exact = a64 + b64
            if abs(sum_exact) > 1e-3 * max(abs(a64), abs(b64))
                rel = abs(convert(Float64, sum_ds) - sum_exact) / abs(sum_exact)
                max_rel_err_add = max(max_rel_err_add, rel)
            end

            # Multiply
            prod_ds = ds_mul(a_ds, b_ds)
            prod_exact = a64 * b64
            if abs(prod_exact) > 1e-30
                rel = abs(convert(Float64, prod_ds) - prod_exact) / abs(prod_exact)
                max_rel_err_mul = max(max_rel_err_mul, rel)
            end

            # Divide (avoid division by near-zero)
            if abs(b64) > 0.01
                div_ds = ds_div(a_ds, b_ds)
                div_exact = a64 / b64
                rel = abs(convert(Float64, div_ds) - div_exact) / abs(div_exact)
                max_rel_err_div = max(max_rel_err_div, rel)
            end
        end

        @printf("  DS add max relative error: %.2e\n", max_rel_err_add)
        @printf("  DS mul max relative error: %.2e\n", max_rel_err_mul)
        @printf("  DS div max relative error: %.2e\n", max_rel_err_div)

        # DS provides ~44-48 mantissa bits; addition's renormalization step
        # accumulates Float32-level error on the lo terms, so worst-case add
        # error is ~eps(Float32)^2 ≈ 1e-9 to 1e-10.  Still vastly better than
        # Float32 alone (~6e-8).  The Dn recursion test below validates that
        # this precision is more than sufficient for Mie computations.
        @test max_rel_err_add < 1e-9
        @test max_rel_err_mul < 1e-13
        @test max_rel_err_div < 1e-12
    end

    @testset "ComplexDS arithmetic" begin
        max_rel_err = 0.0
        for _ in 1:10000
            a = randn(ComplexF64) * 100
            b = randn(ComplexF64) * 100

            a_cds = ComplexDS{Float32}(a)
            b_cds = ComplexDS{Float32}(b)

            # Complex multiply
            prod_cds = cds_mul(a_cds, b_cds)
            prod_exact = a * b
            prod_approx = convert(Complex{Float64}, prod_cds)

            if abs(prod_exact) > 1e-20
                rel = abs(prod_approx - prod_exact) / abs(prod_exact)
                max_rel_err = max(max_rel_err, rel)
            end

            # Complex divide
            if abs(b) > 0.01
                div_cds = cds_div(a_cds, b_cds)
                div_exact = a / b
                div_approx = convert(Complex{Float64}, div_cds)
                rel = abs(div_approx - div_exact) / abs(div_exact)
                max_rel_err = max(max_rel_err, rel)
            end
        end

        @printf("  ComplexDS max relative error: %.2e\n", max_rel_err)
        @test max_rel_err < 1e-12
    end
end

# ============================================================================
# Unit Tests: Neumaier Compensated Summation
# ============================================================================
@testset "Neumaier compensated summation" begin

    @testset "Harmonic series" begin
        # Sum 1/n for n=1..10000 in Float32
        # Naive Float32 loses significant digits
        n = 10000
        exact = sum(1.0 / k for k in 1:n)  # Float64 reference

        # Naive Float32
        naive = Float32(0)
        for k = 1:n
            naive += Float32(1) / Float32(k)
        end

        # Neumaier Float32
        acc = NeumaierAccum{Float32}()
        for k = 1:n
            acc = neumaier_add(acc, Float32(1) / Float32(k))
        end
        compensated = Float64(neumaier_sum(acc))

        naive_err = abs(Float64(naive) - exact)
        comp_err = abs(compensated - exact)

        @printf("  Harmonic(10000): exact=%.10f\n", exact)
        @printf("  Naive Float32 error:    %.2e\n", naive_err)
        @printf("  Neumaier Float32 error: %.2e\n", comp_err)

        # Neumaier should be significantly better
        @test comp_err < naive_err
        @test comp_err < 1e-3  # Should be within ~Float32 precision of exact
    end

    @testset "Alternating series" begin
        # Sum (-1)^k / k for k=1..10000
        n = 10000
        exact = sum((-1.0)^k / k for k in 1:n)

        acc = NeumaierAccum{Float32}()
        naive = Float32(0)
        for k = 1:n
            term = Float32((-1)^k) / Float32(k)
            naive += term
            acc = neumaier_add(acc, term)
        end

        naive_err = abs(Float64(naive) - exact)
        comp_err = abs(Float64(neumaier_sum(acc)) - exact)

        @printf("  Alternating(10000): exact=%.10f\n", exact)
        @printf("  Naive error:    %.2e\n", naive_err)
        @printf("  Neumaier error: %.2e\n", comp_err)

        @test comp_err < naive_err || comp_err < 1e-3
    end

    @testset "Wide dynamic range (size distribution weights)" begin
        # Simulate the 10-order-of-magnitude weight range in Mie size distributions
        n = 2500
        values = Float32.(exp.(range(log(1e-10), log(1e0), length=n)))
        weights = Float32.(randn(n).^2 .+ 0.1)

        exact = sum(Float64.(values) .* Float64.(weights))

        naive = Float32(0)
        acc = NeumaierAccum{Float32}()
        for i = 1:n
            naive += values[i] * weights[i]
            acc = neumaier_add(acc, values[i] * weights[i])
        end

        naive_err = abs(Float64(naive) - exact) / abs(exact)
        comp_err = abs(Float64(neumaier_sum(acc)) - exact) / abs(exact)

        @printf("  Wide range: naive rel err = %.2e, Neumaier rel err = %.2e\n",
                naive_err, comp_err)
        @test comp_err < 1e-5  # Much better than naive
    end
end

# ============================================================================
# Dn Recursion Validation
# ============================================================================
@testset "Dn recursion accuracy" begin

    # Test cases spanning the range of challenging size parameters
    test_cases = [
        (10.0, complex(1.5, -0.1)),
        (50.0, complex(1.5, -0.1)),
        (100.0, complex(1.5, -0.1)),
        (200.0, complex(1.33, 0.0)),
        (100.0, complex(2.0, -1.0)),
        (500.0, complex(1.5, -0.1)),
    ]

    for (x, m) in test_cases
        y = x * m
        n_max = get_n_max(x)
        nmx = round(Int, max(n_max, abs(y)) + 51)

        # CPU Float64 reference
        Dn_ref = zeros(ComplexF64, nmx)
        compute_mie_ab!(x, m, zeros(ComplexF64, n_max), zeros(ComplexF64, n_max), Dn_ref)

        # DS emulated (Float32 pairs)
        y_ds = ComplexDS{Float32}(Complex{Float64}(y))
        Dn_prev = cds_complex(Float32(0), Float32(0))
        Dn_ds = zeros(ComplexF64, n_max)

        for n = (nmx - 1):-1:1
            n_plus_1 = cds_from_real(DoubleSingle{Float32}(Float32(n + 1)))
            ratio = cds_div(n_plus_1, y_ds)
            sum_val = cds_add(Dn_prev, ratio)
            inv_sum = cds_inv(sum_val)
            Dn_prev = cds_sub(ratio, inv_sum)
            if n <= n_max
                Dn_ds[n] = convert(Complex{Float64}, Dn_prev)
            end
        end

        # Compare
        max_rel_err = 0.0
        for n = 1:n_max
            if abs(Dn_ref[n]) > 1e-30
                rel = abs(Dn_ds[n] - Dn_ref[n]) / abs(Dn_ref[n])
                max_rel_err = max(max_rel_err, rel)
            end
        end

        @printf("  Dn(x=%.0f, m=%s): max rel error = %.2e (n_max=%d, nmx=%d)\n",
                x, m, max_rel_err, n_max, nmx)
        @test max_rel_err < 1e-7  # DS provides ~48 bits → should be well within this
    end
end

# ============================================================================
# Full NAI2 Pipeline: GPU vs CPU Reference
# ============================================================================
@testset "NAI2 GPU vs CPU pipeline" begin

    # Standard test aerosol (same as test_Scattering.jl)
    μ_aer = 0.3
    σ_aer = 2.1
    r_max = 30.0
    nᵣ = 1.3
    nᵢ = 0.001
    λ = 0.55

    size_distribution = LogNormal(log(μ_aer), log(σ_aer))
    aero = Aerosol(size_distribution, nᵣ, nᵢ)

    polarization_type = Stokes_IQUV()
    truncation_type = δBGE(10, 10.0)

    for nquad in [100, 500, 1000]
        model = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type,
                               r_max, nquad)

        # CPU reference (existing code)
        ref = compute_aerosol_optical_properties(model)

        # GPU pipeline on CPU backend (NativeFloat64 Dn)
        gpu_f64 = compute_aerosol_optical_properties_gpu(model, KA_CPU();
                    precision_policy=NativeFloat64())

        # GPU pipeline on CPU backend (DSEmulated Dn)
        gpu_ds = compute_aerosol_optical_properties_gpu(model, KA_CPU();
                    precision_policy=DSEmulated())

        @testset "nquad=$nquad, NativeFloat64" begin
            @test isapprox(gpu_f64.ω̃, ref.ω̃, rtol=1e-6)
            @test isapprox(gpu_f64.k, ref.k, rtol=1e-6)

            for (fname, field) in [(:α, :α), (:β, :β), (:γ, :γ),
                                   (:δ, :δ), (:ϵ, :ϵ), (:ζ, :ζ)]
                ref_vals = getproperty(ref.greek_coefs, field)
                gpu_vals = getproperty(gpu_f64.greek_coefs, field)
                @test isapprox(gpu_vals, ref_vals, atol=1e-6)
            end
        end

        @testset "nquad=$nquad, DSEmulated" begin
            # DS provides ~48 bits → slightly less accurate than Float64
            @test isapprox(gpu_ds.ω̃, ref.ω̃, rtol=1e-4)
            @test isapprox(gpu_ds.k, ref.k, rtol=1e-4)

            for (fname, field) in [(:α, :α), (:β, :β), (:γ, :γ),
                                   (:δ, :δ), (:ϵ, :ϵ), (:ζ, :ζ)]
                ref_vals = getproperty(ref.greek_coefs, field)
                gpu_vals = getproperty(gpu_ds.greek_coefs, field)
                @test isapprox(gpu_vals, ref_vals, atol=1e-3)
            end
        end
    end
end

# ============================================================================
# Real CUDA GPU: full NAI2 pipeline on actual hardware (CUDABackend)
# ============================================================================
# Mirrors the KA-CPU accuracy testset above but runs on the real device. Guarded
# by CUDA_FUNCTIONAL so CPU-only machines skip this entirely.
@testset "NAI2 real-GPU (CUDABackend) vs CPU" begin
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — real-GPU pipeline test skipped"
    else
        μ_aer = 0.3; σ_aer = 2.1; r_max = 30.0
        nᵣ = 1.3; nᵢ = 0.001; λ = 0.55
        aero = Aerosol(LogNormal(log(μ_aer), log(σ_aer)), nᵣ, nᵢ)
        polarization_type = Stokes_IQUV()
        truncation_type   = δBGE(10, 10.0)

        for nquad in [100, 500, 1000]
            model = make_mie_model(NAI2(), aero, λ, polarization_type,
                                   truncation_type, r_max, nquad)
            ref = compute_aerosol_optical_properties(model)  # CPU reference

            # Backend-explicit real-GPU call (NativeFloat64 Dn recursion)
            gpu = compute_aerosol_optical_properties_gpu(model, CUDA_BACKEND;
                        precision_policy = NativeFloat64())

            @testset "real-GPU nquad=$nquad, NativeFloat64" begin
                @test isapprox(gpu.ω̃, ref.ω̃, rtol=1e-6)
                @test isapprox(gpu.k,  ref.k,  rtol=1e-6)
                for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    @test isapprox(getproperty(gpu.greek_coefs, field),
                                   getproperty(ref.greek_coefs, field), atol=1e-6)
                end
            end
        end
    end
end

# ============================================================================
# Single-verb dispatch: make_mie_model(...; architecture=GPU()) on real hardware
# ============================================================================
# Exercises the multiple-dispatch API surface end to end:
#   make_mie_model(...; architecture=GPU())  +  compute_aerosol_optical_properties(model)
# must route to the GPU pipeline and match the CPU reference within NativeFloat64
# tolerances (1e-6 SSA/k).
@testset "Single-verb GPU dispatch" begin
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — single-verb GPU dispatch test skipped"
    else
        μ_aer = 0.3; σ_aer = 2.1; r_max = 30.0
        nᵣ = 1.3; nᵢ = 0.001; λ = 0.55
        aero = Aerosol(LogNormal(log(μ_aer), log(σ_aer)), nᵣ, nᵢ)
        pol = Stokes_IQUV(); trunc = δBGE(10, 10.0)

        for nquad in [100, 500]
            m_cpu = make_mie_model(NAI2(), aero, λ, pol, trunc, r_max, nquad)
            m_gpu = make_mie_model(NAI2(), aero, λ, pol, trunc, r_max, nquad;
                                   architecture = GPU())

            # Type-parameter sanity: architecture rides on the model type
            @test m_gpu isa MieModel{NAI2, Float64, GPU}
            @test m_cpu isa MieModel{NAI2, Float64, CPU}

            ref = compute_aerosol_optical_properties(m_cpu)  # CPU dispatch
            gpu = compute_aerosol_optical_properties(m_gpu)  # GPU dispatch (auto NativeFloat64)

            @testset "single-verb nquad=$nquad" begin
                @test isapprox(gpu.ω̃, ref.ω̃, rtol=1e-6)
                @test isapprox(gpu.k,  ref.k,  rtol=1e-6)
                for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    @test isapprox(getproperty(gpu.greek_coefs, field),
                                   getproperty(ref.greek_coefs, field), atol=1e-6)
                end
            end
        end
    end
end

# ============================================================================
# Float32 GPU model: compensated host reduction keeps F32 accurate
# ============================================================================
# A genuinely Float32 model (FT === Float32) runs the device kernels in Float32
# (DSEmulated Dn recursion) — this is the only option on F32-only GPUs (e.g.
# Metal). The host-side size-distribution reduction and Greek-coefficient
# quadrature are widened to Float64 with Neumaier compensation, so the Greek
# coefficients land at the Float32 representational floor (~1e-4) instead of the
# ~5e-3 they would hit if the reduction itself ran in Float32. SSA/k reach ~1e-7.
@testset "Float32 GPU compensated-reduction accuracy" begin
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — Float32 GPU accuracy test skipped"
    else
        # Float64 reference at identical physical parameters
        aero64 = Aerosol(LogNormal(log(0.3),  log(2.1)),  1.3,   0.001)
        m64 = make_mie_model(NAI2(), aero64, 0.55, Stokes_IQUV(), δBGE(10, 10.0), 30.0, 500)
        ref = compute_aerosol_optical_properties(m64)

        # True Float32 model on the GPU with DSEmulated (Float32) Dn recursion
        aero32 = Aerosol(LogNormal(log(0.3f0), log(2.1f0)), 1.3f0, 0.001f0)
        m32 = make_mie_model(NAI2(), aero32, 0.55f0, Stokes_IQUV(), δBGE(10, 10.0f0),
                             30.0f0, 500; architecture = GPU(),
                             precision_policy = DSEmulated())
        gpu32 = compute_aerosol_optical_properties(m32)

        # Output type rides on the model's FT (Float32 here)
        @test eltype(gpu32.greek_coefs.β) === Float32
        @test gpu32.ω̃ isa Float32

        # SSA / k accurate to ~1e-7 thanks to compensated bulk reductions
        @test abs(Float64(gpu32.ω̃) - ref.ω̃) / abs(ref.ω̃) < 1e-6
        @test abs(Float64(gpu32.k)  - ref.k)  / abs(ref.k)  < 1e-6

        # Greek coefficients at the Float32 floor (compensated reduction → ~1e-4,
        # not the ~5e-3 of an all-Float32 reduction). Compare on the common length.
        L = min(length(gpu32.greek_coefs.β), length(ref.greek_coefs.β))
        for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
            g = Float64.(getproperty(gpu32.greek_coefs, field)[1:L])
            r = getproperty(ref.greek_coefs, field)[1:L]
            @test maximum(abs.(g .- r)) < 1e-3
        end
    end
end

# ============================================================================
# Float32 GPU DEFAULT precision policy (precision_policy = nothing → DSEmulated)
# ============================================================================
# Regression for the P1 finding: a Float32 GPU model built WITHOUT an explicit
# precision_policy (nothing → auto) must NOT resolve to NativeFloat64 (whose GPU
# kernel asserts FT === Float64) — it must auto-select DSEmulated, the
# Float32-native path. Several existing GPU YAMLs use float_type: Float32 and
# never set a policy, so this is the path they actually hit. The DSEmulated host
# reduction is Float64-widened + Neumaier-compensated, so accuracy is at the
# Float32 representational floor (1e-4 SSA/k, established DS tolerances).
@testset "Float32 GPU default-policy auto-selects DSEmulated" begin
    # FT-aware policy selection is testable WITHOUT a GPU device: the policy
    # function is provided by the loaded CUDA extension regardless of
    # CUDA.functional(). Guard the actual compute on CUDA_FUNCTIONAL.
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — Float32 GPU default-policy test skipped"
    else
        # 1) The FT-aware default policy itself (provided by the CUDA extension)
        @test vSmartMOM.Architectures.default_mie_precision_policy(GPU(), Float32) isa DSEmulated
        @test vSmartMOM.Architectures.default_mie_precision_policy(GPU(), Float64) isa NativeFloat64

        # 2) Float64 reference at identical physical parameters
        aero64 = Aerosol(LogNormal(log(0.3), log(2.1)), 1.3, 0.001)
        m64 = make_mie_model(NAI2(), aero64, 0.55, Stokes_IQUV(), δBGE(10, 10.0), 30.0, 500)
        ref = compute_aerosol_optical_properties(m64)

        # 3) Float32 GPU model with NO explicit policy (precision_policy = nothing).
        #    FT === Float32 (derived from the Float32 r_max). This MUST run — the
        #    router resolves nothing → DSEmulated() via the FT-aware default — not
        #    trip the NativeFloat64 FT === Float64 assert.
        aero32 = Aerosol(LogNormal(log(0.3f0), log(2.1f0)), 1.3f0, 0.001f0)
        m32 = make_mie_model(NAI2(), aero32, 0.55f0, Stokes_IQUV(), δBGE(10, 10.0f0),
                             30.0f0, 500; architecture = GPU())
        @test m32 isa MieModel{NAI2, Float32, GPU}
        @test m32.precision_policy === nothing   # auto

        gpu32 = compute_aerosol_optical_properties(m32)   # would assert if NativeFloat64

        # Output type rides on the model's FT (Float32)
        @test eltype(gpu32.greek_coefs.β) === Float32
        @test gpu32.ω̃ isa Float32

        # Matches the CPU Float64 reference within the established DS tolerances.
        @test abs(Float64(gpu32.ω̃) - ref.ω̃) / abs(ref.ω̃) < 1e-4
        @test abs(Float64(gpu32.k)  - ref.k)  / abs(ref.k)  < 1e-4

        L = min(length(gpu32.greek_coefs.β), length(ref.greek_coefs.β))
        for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
            g = Float64.(getproperty(gpu32.greek_coefs, field)[1:L])
            r = getproperty(ref.greek_coefs, field)[1:L]
            @test maximum(abs.(g .- r)) < 1e-3
        end
    end
end

# ============================================================================
# P2 regression: Float32 aerosol params + Float64 r_max/λ → auto-selects DSEmulated
# ============================================================================
# The reviewer's exact case: Aerosol(dist, 1.3f0, 0.001f0) combined with a
# Float64 r_max (30.0) and Float64 λ (0.55). In this model:
#   - _mie_output_type  → Float64  (r_max's type drives the MieModel FT parameter)
#   - _mie_kernel_type  → Float32  (eltype(nᵣ) — how compute_NAI2_gpu derives FT)
# Before the fix the router passed FT2 (Float64) to default_mie_precision_policy,
# returning NativeFloat64(), whose GPU kernel asserts FT === Float64 — but
# FT = eltype(nᵣ) = Float32 → assertion failure. The fix uses _mie_kernel_type
# (Float32) for the policy decision, returning DSEmulated(), which is the
# Float32-native path and must run without error and match the CPU reference.
@testset "P2: Float32 aerosol + Float64 r_max — auto-selects DSEmulated" begin
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — P2 mixed-precision policy test skipped"
    else
        # The reviewer's exact construction: Float32 nᵣ/nᵢ, Float64 r_max/λ.
        aero_mixed = Aerosol(LogNormal(log(0.3), log(2.1)), 1.3f0, 0.001f0)
        m_mixed = make_mie_model(NAI2(), aero_mixed, 0.55, Stokes_IQUV(),
                                 δBGE(10, 10.0), 30.0, 500;
                                 architecture = GPU())

        # Confirm the type-level situation that triggered the bug:
        #   FT (r_max type) = Float64  →  MieModel{NAI2, Float64, GPU}
        #   eltype(nᵣ) = Float32       →  kernel must run Float32 (DSEmulated)
        @test m_mixed isa MieModel{NAI2, Float64, GPU}
        @test eltype(m_mixed.aerosol.nᵣ) === Float32   # kernel precision axis
        @test m_mixed.precision_policy === nothing       # nothing → auto

        # Before the fix this would @assert-fail with "NativeFloat64 … requires Float64".
        gpu_mixed = compute_aerosol_optical_properties(m_mixed)

        # The auto policy must have resolved to DSEmulated (Float32-native path).
        # Verify indirectly: output type rides on FT2 = Float64 (r_max type).
        @test gpu_mixed.ω̃ isa Float64
        @test eltype(gpu_mixed.greek_coefs.β) === Float64

        # Float64 CPU reference at identical physical parameters.
        aero64 = Aerosol(LogNormal(log(0.3), log(2.1)), 1.3, 0.001)
        m64    = make_mie_model(NAI2(), aero64, 0.55, Stokes_IQUV(), δBGE(10, 10.0), 30.0, 500)
        ref    = compute_aerosol_optical_properties(m64)

        # DSEmulated tolerance (Float32 kernel, Float64-widened host reduction):
        # SSA/k within 1e-4, Greek coefficients within 1e-3.
        @test abs(gpu_mixed.ω̃ - ref.ω̃) / abs(ref.ω̃) < 1e-4
        @test abs(gpu_mixed.k  - ref.k)  / abs(ref.k)  < 1e-4

        L = min(length(gpu_mixed.greek_coefs.β), length(ref.greek_coefs.β))
        for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
            g = getproperty(gpu_mixed.greek_coefs, field)[1:L]
            r = getproperty(ref.greek_coefs,       field)[1:L]
            @test maximum(abs.(g .- r)) < 1e-3
        end
    end
end

# ============================================================================
# Non-CPU architecture WITHOUT a GPU Mie pipeline → warn-once CPU fallback
# ============================================================================
# Regression for the P1 finding: a non-CPU architecture that has no GPU Mie
# pipeline (the trait `has_gpu_mie` is false) must NOT route through the GPU
# branch (which would call ka_backend on an architecture with no backend hook
# and error). Instead the router must warn ONCE and compute Mie on the CPU,
# returning results bit-identical to the CPU path (RT arrays are unaffected;
# only Mie falls back).
#
# NOTE: as of the NativeFloat32 + Metal Mie registration work, `MetalGPU()`
# itself is NO LONGER an example of this situation once `vSmartMOMMetalExt` is
# loaded -- the Metal extension now registers `has_gpu_mie(::MetalGPU) = true`
# (NativeFloat32 route; see ext/vSmartMOMMetalExt.jl). This testset instead
# uses a minimal LOCAL architecture type that genuinely has no registered GPU
# Mie pipeline, to keep exercising the fallback trait mechanism itself.
#
# We test the trait mechanism directly with a minimal local AbstractArchitecture
# subtype: `has_gpu_mie` defaults to false on the abstract supertype, so the
# router takes the CPU branch. This needs no GPU device, so it runs everywhere.
#
# Minimal architecture with NO ka_backend / precision-policy / Mie hooks, defined
# at top level (structs cannot be declared inside the function `@testset` wraps).
struct _NoMieArch <: vSmartMOM.Architectures.AbstractArchitecture end

@testset "Non-CPU no-GPU-Mie architecture → warn-once CPU fallback" begin
    # Trait default is false (refinement is opt-in, like ka_backend).
    @test vSmartMOM.Architectures.has_gpu_mie(_NoMieArch()) == false

    aero = Aerosol(LogNormal(log(0.3), log(2.1)), 1.3, 0.001)
    pol  = Stokes_IQUV(); trunc = δBGE(10, 10.0)

    m_cpu  = make_mie_model(NAI2(), aero, 0.55, pol, trunc, 30.0, 200)
    m_fall = make_mie_model(NAI2(), aero, 0.55, pol, trunc, 30.0, 200;
                            architecture = _NoMieArch())

    ref = compute_aerosol_optical_properties(m_cpu)

    # Reset the module-global warn-once flag so this test is self-contained
    # regardless of whether an earlier testset already tripped the fallback.
    vSmartMOM.Scattering._WARNED_NO_GPU_MIE[] = false

    # First call warns exactly once; @test_logs asserts a matching @warn fires.
    fb1 = @test_logs (:warn, r"no GPU Mie pipeline") match_mode = :any begin
        compute_aerosol_optical_properties(m_fall)
    end
    # Subsequent calls must NOT warn again (warn-once flag): assert no
    # Warn-level (or above) records are produced.
    fb2 = @test_logs min_level = Logging.Warn begin
        compute_aerosol_optical_properties(m_fall)
    end

    # Fallback results are IDENTICAL to the CPU path (not just approximate).
    @test fb1.ω̃ == ref.ω̃
    @test fb1.k  == ref.k
    for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
        @test getproperty(fb1.greek_coefs, field) == getproperty(ref.greek_coefs, field)
        @test getproperty(fb2.greek_coefs, field) == getproperty(ref.greek_coefs, field)
    end
end

# ============================================================================
# YAML-level smoke: model_from_parameters with architecture=GPU()
# ============================================================================
# Build a full RTModel from a real aerosol YAML, once on CPU and once with the
# params' architecture flipped to GPU(). The per-band truncated aerosol_optics
# must match within NativeFloat64 tolerances (1e-6 SSA/k). This validates the
# production wiring in model_from_parameters.jl (single-verb call +
# architecture=params.architecture).
@testset "YAML model_from_parameters GPU vs CPU aerosol_optics" begin
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — YAML GPU model_from_parameters test skipped"
    else
        params_cpu = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
        params_gpu = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
        params_gpu.architecture = GPU()

        model_cpu = model_from_parameters(params_cpu)
        model_gpu = model_from_parameters(params_gpu)

        ao_cpu = model_cpu.aerosol_optics   # [iBand][iAer]
        ao_gpu = model_gpu.aerosol_optics

        @test length(ao_gpu) == length(ao_cpu)
        for iBand in eachindex(ao_cpu)
            for iAer in eachindex(ao_cpu[iBand])
                a = ao_cpu[iBand][iAer]
                b = ao_gpu[iBand][iAer]
                @test isapprox(b.ω̃, a.ω̃, rtol=1e-6)
                @test isapprox(b.k,  a.k,  rtol=1e-6)
                @test isapprox(b.fᵗ, a.fᵗ, atol=1e-6)
                for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
                    @test isapprox(getproperty(b.greek_coefs, field),
                                   getproperty(a.greek_coefs, field), atol=1e-6)
                end
            end
        end
    end
end

# Minimal mock backend registering `Scattering._is_metal_backend(...) = true`
# (real dispatch, not name matching -- see gpu_precision.jl), so the
# Metal-only-NativeFloat32 guard in compute_aerosol_optical_properties_gpu can
# be exercised WITHOUT installing Metal.jl. Defined at top level (structs
# cannot be declared inside a function/testset body before Julia 1.12).
struct _FakeMetalBackendMie end
vSmartMOM.Scattering._is_metal_backend(::_FakeMetalBackendMie) = true

# ============================================================================
# NativeFloat32 (log-normal MieModel GPU path): genuine kernel-1 dispatch
# ============================================================================
# Metal-capability-leak fix follow-up: `has_gpu_mie` is a single
# architecture-level trait shared by this log-normal path
# (compute_aerosol_optical_properties_gpu) and the caller-node GPU path, so
# registering Metal for one registers it for both. This path's Kernel 1 now
# dispatches via the SAME `Scattering._mie_kernel1(policy, backend)` helper
# the caller-node path uses (NativeFloat64/NativeFloat32 -> native kernel,
# DSEmulated -> double-single kernel), and carries the same Metal-only
# Float64-device-array guard (`Scattering._is_metal_backend`). UNLIKE the
# caller-node path, this path's host-side size-distribution reduction is
# ALWAYS Float64-widened regardless of policy (see the `RA` comment in
# compute_NAI2_gpu.jl), so `NativeFloat32` here is expected to land close to
# `DSEmulated`'s accuracy (not the node path's ~5e-3 Greek floor) -- only the
# Dₙ-recursion kernel differs (native Float32 vs Float32 double-single
# pairs), and that recursion is not the accuracy bottleneck once the
# reduction itself is Float64-compensated.
@testset "NativeFloat32 (log-normal MieModel GPU path): KA-CPU accuracy + kernel dispatch" begin
    aero64 = Aerosol(LogNormal(log(0.3), log(2.1)), 1.3, 0.001)
    m64 = make_mie_model(NAI2(), aero64, 0.55, Stokes_IQUV(), δBGE(10, 10.0), 30.0, 500)
    ref = compute_aerosol_optical_properties(m64)

    aero32 = Aerosol(LogNormal(log(0.3f0), log(2.1f0)), 1.3f0, 0.001f0)
    m32 = make_mie_model(NAI2(), aero32, 0.55f0, Stokes_IQUV(), δBGE(10, 10.0f0), 30.0f0, 500)
    @test m32 isa MieModel{NAI2, Float32, CPU}

    # FT2 defaults to Float64 on this low-level entry point (unlike the
    # single-verb compute_aerosol_optical_properties(model), which supplies
    # FT2 = _mie_output_type(model) automatically) -- pass it explicitly so
    # the Float32 model's output stays Float32, matching the sibling
    # "Float32 GPU compensated-reduction accuracy" testset's expectations.
    gpu_native = compute_aerosol_optical_properties_gpu(m32, KA_CPU(); precision_policy=NativeFloat32(), FT2=Float32)
    gpu_ds     = compute_aerosol_optical_properties_gpu(m32, KA_CPU(); precision_policy=DSEmulated(), FT2=Float32)

    @test eltype(gpu_native.greek_coefs.β) === Float32
    @test gpu_native.ω̃ isa Float32

    k_relerr_native = abs(Float64(gpu_native.k) - ref.k) / ref.k
    ω_relerr_native = abs(Float64(gpu_native.ω̃) - ref.ω̃) / ref.ω̃
    @printf("  NativeFloat32 (log-normal, KA-CPU): k relerr=%.2e  ω̃ relerr=%.2e\n", k_relerr_native, ω_relerr_native)

    # Same Float64-widened host reduction regardless of policy => NativeFloat32
    # lands close to (not necessarily better/worse than) DSEmulated, both at
    # the ~1e-6/1e-4 established DS tolerances -- NOT the node path's looser
    # ~1e-2 Greek floor (that path's reduction is never widened for
    # NativeFloat32, by design; this path's is, always).
    @test k_relerr_native < 1e-6
    @test ω_relerr_native < 1e-6
    L = min(length(gpu_native.greek_coefs.β), length(ref.greek_coefs.β))
    for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
        g = Float64.(getproperty(gpu_native.greek_coefs, field)[1:L])
        r_ = getproperty(ref.greek_coefs, field)[1:L]
        @test maximum(abs.(g .- r_)) < 1e-3
    end

    # Metal-only-NativeFloat32 guard (mock backend, no real Metal.jl needed):
    # NativeFloat64 (or any FT===Float64 model) on a Metal-shaped backend must
    # throw ArgumentError -- proves the guard added to compute_NAI2_gpu.jl
    # actually fires on this path too, not just the caller-node path.
    @test_throws ArgumentError compute_aerosol_optical_properties_gpu(
        m64, _FakeMetalBackendMie(); precision_policy=NativeFloat64())
end

@testset "NativeFloat32 (log-normal MieModel GPU path): real CUDA (gated)" begin
    if !CUDA_FUNCTIONAL
        @test_skip "CUDA not functional — NativeFloat32 log-normal real-GPU test skipped"
    else
        aero64 = Aerosol(LogNormal(log(0.3), log(2.1)), 1.3, 0.001)
        m64 = make_mie_model(NAI2(), aero64, 0.55, Stokes_IQUV(), δBGE(10, 10.0), 30.0, 500)
        ref = compute_aerosol_optical_properties(m64)

        aero32 = Aerosol(LogNormal(log(0.3f0), log(2.1f0)), 1.3f0, 0.001f0)
        m32 = make_mie_model(NAI2(), aero32, 0.55f0, Stokes_IQUV(), δBGE(10, 10.0f0), 30.0f0, 500)

        ka_native = compute_aerosol_optical_properties_gpu(m32, KA_CPU(); precision_policy=NativeFloat32(), FT2=Float32)
        cu_native = compute_aerosol_optical_properties_gpu(m32, CUDA_BACKEND; precision_policy=NativeFloat32(), FT2=Float32)

        @test isapprox(cu_native.k, ka_native.k, rtol=1e-6)
        @test isapprox(cu_native.ω̃, ka_native.ω̃, rtol=1e-6)
        for field in (:α, :β, :γ, :δ, :ϵ, :ζ)
            @test isapprox(getproperty(cu_native.greek_coefs, field),
                           getproperty(ka_native.greek_coefs, field), atol=1e-4)
        end

        k_relerr = abs(Float64(cu_native.k) - ref.k) / ref.k
        ω_relerr = abs(Float64(cu_native.ω̃) - ref.ω̃) / ref.ω̃
        @printf("  NativeFloat32 (log-normal, real CUDA): k relerr=%.2e  ω̃ relerr=%.2e\n", k_relerr, ω_relerr)
        @test k_relerr < 1e-6
        @test ω_relerr < 1e-6
    end
end

# ============================================================================
# Performance Benchmarks
# ============================================================================
# The benchmark sweep and the nquad=2500 accuracy table below take several
# minutes on the KA-CPU backend, so they are opt-in: run with
# VSMARTMOM_MIE_BENCH=1 to include them.
if get(ENV, "VSMARTMOM_MIE_BENCH", "0") == "1"

@testset "Performance comparison" begin

    μ_aer = 0.3
    σ_aer = 2.1
    r_max = 30.0
    nᵣ = 1.3
    nᵢ = 0.001
    λ = 0.55

    size_distribution = LogNormal(log(μ_aer), log(σ_aer))
    aero = Aerosol(size_distribution, nᵣ, nᵢ)
    polarization_type = Stokes_IQUV()
    truncation_type = δBGE(10, 10.0)

    println("\n" * "="^70)
    println("Performance Benchmarks: CPU baseline vs GPU kernels (CPU backend)")
    println("="^70)
    @printf("%-12s  %12s  %12s  %12s  %8s\n",
            "nquad", "CPU (s)", "GPU-F64 (s)", "GPU-DS (s)", "Speedup")
    println("-"^70)

    for nquad in [100, 500, 1000, 2500]
        model = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type,
                               r_max, nquad)

        # Warmup
        compute_aerosol_optical_properties(model)
        compute_aerosol_optical_properties_gpu(model, KA_CPU(); precision_policy=NativeFloat64())
        compute_aerosol_optical_properties_gpu(model, KA_CPU(); precision_policy=DSEmulated())

        # Time CPU baseline
        t_cpu = @elapsed for _ in 1:3
            compute_aerosol_optical_properties(model)
        end
        t_cpu /= 3

        # Time GPU-F64 (CPU backend)
        t_gpu_f64 = @elapsed for _ in 1:3
            compute_aerosol_optical_properties_gpu(model, KA_CPU(); precision_policy=NativeFloat64())
        end
        t_gpu_f64 /= 3

        # Time GPU-DS (CPU backend)
        t_gpu_ds = @elapsed for _ in 1:3
            compute_aerosol_optical_properties_gpu(model, KA_CPU(); precision_policy=DSEmulated())
        end
        t_gpu_ds /= 3

        @printf("%-12d  %12.4f  %12.4f  %12.4f  %7.2fx\n",
                nquad, t_cpu, t_gpu_f64, t_gpu_ds, t_cpu / t_gpu_f64)
    end
    println("="^70)
    println("Note: GPU-backend benchmarks on actual CUDA hardware require")
    println("      CUDA.jl loaded. Run with: julia --project -e 'using CUDA; ...'")
end

# ============================================================================
# Accuracy Summary Table
# ============================================================================
@testset "Accuracy summary" begin
    μ_aer = 0.3
    σ_aer = 2.1
    r_max = 30.0
    nᵣ = 1.3
    nᵢ = 0.001
    λ = 0.55
    nquad = 2500

    size_distribution = LogNormal(log(μ_aer), log(σ_aer))
    aero = Aerosol(size_distribution, nᵣ, nᵢ)
    polarization_type = Stokes_IQUV()
    truncation_type = δBGE(10, 10.0)

    model = make_mie_model(NAI2(), aero, λ, polarization_type, truncation_type,
                           r_max, nquad)

    ref = compute_aerosol_optical_properties(model)
    gpu_f64 = compute_aerosol_optical_properties_gpu(model, KA_CPU(); precision_policy=NativeFloat64())
    gpu_ds  = compute_aerosol_optical_properties_gpu(model, KA_CPU(); precision_policy=DSEmulated())

    println("\n" * "="^70)
    println("Accuracy Summary (nquad=$nquad)")
    println("="^70)
    @printf("%-12s  %15s  %15s\n", "Metric", "GPU-F64 err", "GPU-DS err")
    println("-"^70)

    # SSA
    ssa_err_f64 = abs(gpu_f64.ω̃ - ref.ω̃) / abs(ref.ω̃)
    ssa_err_ds  = abs(gpu_ds.ω̃ - ref.ω̃) / abs(ref.ω̃)
    @printf("%-12s  %15.2e  %15.2e\n", "SSA (rel)", ssa_err_f64, ssa_err_ds)

    # Extinction
    k_err_f64 = abs(gpu_f64.k - ref.k) / abs(ref.k)
    k_err_ds  = abs(gpu_ds.k - ref.k) / abs(ref.k)
    @printf("%-12s  %15.2e  %15.2e\n", "k_ext (rel)", k_err_f64, k_err_ds)

    # Greek coefficients
    for (name, field) in [("alpha", :α), ("beta", :β), ("gamma", :γ),
                          ("delta", :δ), ("epsilon", :ϵ), ("zeta", :ζ)]
        ref_v = getproperty(ref.greek_coefs, field)
        f64_v = getproperty(gpu_f64.greek_coefs, field)
        ds_v  = getproperty(gpu_ds.greek_coefs, field)

        max_abs_f64 = maximum(abs.(f64_v .- ref_v))
        max_abs_ds  = maximum(abs.(ds_v .- ref_v))
        @printf("%-12s  %15.2e  %15.2e\n", "$name (abs)", max_abs_f64, max_abs_ds)
    end
    println("="^70)

    # Assertions
    @test ssa_err_f64 < 1e-6
    @test ssa_err_ds  < 1e-4
    @test k_err_f64   < 1e-6
    @test k_err_ds    < 1e-4
end

end # if VSMARTMOM_MIE_BENCH
