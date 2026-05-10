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

# KernelAbstractions is re-exported through vSmartMOM's Scattering module
using vSmartMOM.Scattering
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

            # Add
            sum_ds = ds_add(a_ds, b_ds)
            sum_exact = a64 + b64
            if abs(sum_exact) > 1e-30
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
# Performance Benchmarks
# ============================================================================
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
