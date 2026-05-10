#!/usr/bin/env julia
# ==========================================================================
# Mie Code Float32 vs Float64 Stability Benchmark
# ==========================================================================
# Tests whether Float32 produces accurate Mie coefficients and optical
# properties across a range of size parameters (x = 2πr/λ).
# The downward recursion for Dₙ is the most numerically sensitive part.
# ==========================================================================

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using vSmartMOM
using vSmartMOM.Scattering
using Distributions
using Printf
using Statistics
using LinearAlgebra

println("=" ^ 72)
println("  Float32 vs Float64 Mie Stability Test")
println("=" ^ 72)

# ── Part 1: Single-particle Mie coefficient comparison ───────────────────

println("\n─── Part 1: Single-particle aₙ, bₙ accuracy ───\n")
println("  Testing compute_mie_ab! for individual size parameters")
println("  Refractive index: m = 1.3 - 1e-8i\n")

nᵣ = 1.3
nᵢ = 1e-8
m = nᵣ - nᵢ * im

size_params = [0.1, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0,
               100.0, 150.0, 200.0, 300.0, 500.0]

@printf("  %8s  %6s  %12s  %12s  %12s  %12s\n",
        "x", "n_max", "max|Δaₙ|/|aₙ|", "max|Δbₙ|/|bₙ|", "max|ΔDₙ|/|Dₙ|", "status")
println("  " * "─" ^ 72)

for x in size_params
    n_max = Scattering.get_n_max(x)
    y = x * m
    nmx = round(Int, max(n_max, abs(y)) + 51)

    # Float64 reference
    an64 = zeros(ComplexF64, n_max)
    bn64 = zeros(ComplexF64, n_max)
    Dn64 = zeros(ComplexF64, nmx)
    Scattering.compute_mie_ab!(Float64(x), ComplexF64(m), an64, bn64, Dn64)

    # Float32
    an32 = zeros(ComplexF32, n_max)
    bn32 = zeros(ComplexF32, n_max)
    Dn32 = zeros(ComplexF32, nmx)
    Scattering.compute_mie_ab!(Float32(x), ComplexF32(m), an32, bn32, Dn32)

    # Relative errors
    mask_a = abs.(an64) .> 1e-30
    mask_b = abs.(bn64) .> 1e-30
    mask_D = abs.(Dn64) .> 1e-30

    err_a = any(mask_a) ? maximum(abs.(an64[mask_a] .- ComplexF64.(an32[mask_a])) ./ abs.(an64[mask_a])) : 0.0
    err_b = any(mask_b) ? maximum(abs.(bn64[mask_b] .- ComplexF64.(bn32[mask_b])) ./ abs.(bn64[mask_b])) : 0.0
    err_D = any(mask_D) ? maximum(abs.(Dn64[mask_D] .- ComplexF64.(Dn32[mask_D])) ./ abs.(Dn64[mask_D])) : 0.0

    status = max(err_a, err_b) < 1e-3 ? "OK" :
             max(err_a, err_b) < 1e-1 ? "WARN" : "FAIL"

    @printf("  %8.1f  %6d  %12.2e  %12.2e  %12.2e  %12s\n",
            x, n_max, err_a, err_b, err_D, status)
end

# ── Part 2: Detailed breakdown at critical size parameters ───────────────

println("\n─── Part 2: Coefficient-by-coefficient analysis at x=50, 100, 200 ───\n")

for x in [50.0, 100.0, 200.0]
    n_max = Scattering.get_n_max(x)
    y = x * m
    nmx = round(Int, max(n_max, abs(y)) + 51)

    an64 = zeros(ComplexF64, n_max)
    bn64 = zeros(ComplexF64, n_max)
    Dn64 = zeros(ComplexF64, nmx)
    Scattering.compute_mie_ab!(Float64(x), ComplexF64(m), an64, bn64, Dn64)

    an32 = zeros(ComplexF32, n_max)
    bn32 = zeros(ComplexF32, n_max)
    Dn32 = zeros(ComplexF32, nmx)
    Scattering.compute_mie_ab!(Float32(x), ComplexF32(m), an32, bn32, Dn32)

    rel_err_a = abs.(an64 .- ComplexF64.(an32)) ./ max.(abs.(an64), 1e-30)
    rel_err_b = abs.(bn64 .- ComplexF64.(bn32)) ./ max.(abs.(bn64), 1e-30)

    # Find where errors exceed threshold
    bad_a = findall(rel_err_a .> 1e-3)
    bad_b = findall(rel_err_b .> 1e-3)

    @printf("  x = %.0f (n_max = %d):\n", x, n_max)
    if isempty(bad_a) && isempty(bad_b)
        @printf("    All coefficients within 0.1%% accuracy\n")
    else
        if !isempty(bad_a)
            @printf("    aₙ > 0.1%% error at n = %s (first 5)\n",
                    join(string.(bad_a[1:min(5, length(bad_a))]), ", "))
            @printf("    Max aₙ error: %.2e at n=%d\n",
                    maximum(rel_err_a), argmax(rel_err_a))
        end
        if !isempty(bad_b)
            @printf("    bₙ > 0.1%% error at n = %s (first 5)\n",
                    join(string.(bad_b[1:min(5, length(bad_b))]), ", "))
            @printf("    Max bₙ error: %.2e at n=%d\n",
                    maximum(rel_err_b), argmax(rel_err_b))
        end
    end

    # Dₙ error distribution (first n_max elements)
    Dn_rel = abs.(Dn64[1:n_max] .- ComplexF64.(Dn32[1:n_max])) ./
             max.(abs.(Dn64[1:n_max]), 1e-30)
    @printf("    Dₙ errors (first %d): median=%.2e, 95th=%.2e, max=%.2e\n",
            n_max, median(Dn_rel), quantile(Dn_rel, 0.95), maximum(Dn_rel))
    println()
end

# ── Part 3: Full size-distribution-averaged properties ───────────────────

println("─── Part 3: Full aerosol optics comparison (Float32 vs Float64) ───\n")

test_configs = [
    (r_max = 5.0,  nq = 200, label = "Small (r_max=5)"),
    (r_max = 10.0, nq = 200, label = "Medium (r_max=10)"),
    (r_max = 30.0, nq = 200, label = "Large (r_max=30)"),
    (r_max = 50.0, nq = 200, label = "Very large (r_max=50)"),
]

λ = 0.55
sd64 = LogNormal(log(0.3), log(1.5))

for cfg in test_configs
    println("  $(cfg.label), λ=$(λ) μm, nquad=$(cfg.nq)")

    # Float64 reference
    aerosol64 = Aerosol(sd64, nᵣ, nᵢ)
    pol64 = Stokes_IQU{Float64}()
    trunc64 = δBGE(20, 2.0)
    model64 = make_mie_model(NAI2(), aerosol64, λ, pol64, trunc64, cfg.r_max, cfg.nq)
    t64 = @timed begin
        optics64 = compute_aerosol_optical_properties(model64)
    end

    # Float32 -- need Float32 types throughout
    sd32 = LogNormal(log(Float32(0.3)), log(Float32(1.5)))
    aerosol32 = Aerosol(sd32, Float32(nᵣ), Float32(nᵢ))
    pol32 = Stokes_IQU{Float32}()
    trunc32 = δBGE(20, Float32(2.0))
    model32 = make_mie_model(NAI2(), aerosol32, Float32(λ), pol32, trunc32,
                             Float32(cfg.r_max), cfg.nq)
    local optics32, t32_time
    try
        t32 = @timed begin
            optics32 = compute_aerosol_optical_properties(model32, Float32)
        end
        t32_time = t32.time
    catch e
        @printf("    Float32 FAILED: %s\n\n", sprint(showerror, e))
        continue
    end

    optics64_result = t64.value

    # Compare scalar quantities
    err_ω = abs(Float64(optics32.ω̃) - optics64_result.ω̃) /
            max(abs(optics64_result.ω̃), 1e-30)
    err_k = abs(Float64(optics32.k) - optics64_result.k) /
            max(abs(optics64_result.k), 1e-30)

    # Compare Greek coefficients
    β64 = optics64_result.greek_coefs.β
    β32 = Float64.(optics32.greek_coefs.β)
    mask_β = abs.(β64) .> 1e-10
    err_β = any(mask_β) ? maximum(abs.(β64[mask_β] .- β32[mask_β]) ./ abs.(β64[mask_β])) : 0.0

    α64 = optics64_result.greek_coefs.α
    α32 = Float64.(optics32.greek_coefs.α)
    mask_α = abs.(α64) .> 1e-10
    err_α = any(mask_α) ? maximum(abs.(α64[mask_α] .- α32[mask_α]) ./ abs.(α64[mask_α])) : 0.0

    status = max(err_ω, err_k, err_β) < 1e-3 ? "OK" :
             max(err_ω, err_k, err_β) < 1e-1 ? "WARN" : "FAIL"

    @printf("    ω̃:  F64=%.8f  F32=%.8f  rel_err=%.2e\n",
            optics64_result.ω̃, optics32.ω̃, err_ω)
    @printf("    k:   F64=%.8e  F32=%.8e  rel_err=%.2e\n",
            optics64_result.k, optics32.k, err_k)
    @printf("    β:   max_rel_err=%.2e (over %d nonzero coeffs)\n",
            err_β, sum(mask_β))
    @printf("    α:   max_rel_err=%.2e (over %d nonzero coeffs)\n",
            err_α, sum(mask_α))
    @printf("    Speed: F64=%.3fs  F32=%.3fs  (%.1fx)\n",
            t64.time, t32_time, t64.time / t32_time)
    @printf("    Status: %s\n\n", status)
end

# ── Part 4: Absorbing aerosol (higher nᵢ) ───────────────────────────────

println("─── Part 4: Absorbing aerosol (nᵢ=0.01, r_max=20) ───\n")

nᵢ_abs = 0.01
m_abs = nᵣ - nᵢ_abs * im

for x in [10.0, 50.0, 100.0, 200.0]
    n_max = Scattering.get_n_max(x)
    y = x * m_abs
    nmx = round(Int, max(n_max, abs(y)) + 51)

    an64 = zeros(ComplexF64, n_max)
    bn64 = zeros(ComplexF64, n_max)
    Dn64 = zeros(ComplexF64, nmx)
    Scattering.compute_mie_ab!(Float64(x), ComplexF64(m_abs), an64, bn64, Dn64)

    an32 = zeros(ComplexF32, n_max)
    bn32 = zeros(ComplexF32, n_max)
    Dn32 = zeros(ComplexF32, nmx)
    Scattering.compute_mie_ab!(Float32(x), ComplexF32(m_abs), an32, bn32, Dn32)

    mask_a = abs.(an64) .> 1e-30
    mask_b = abs.(bn64) .> 1e-30
    err_a = any(mask_a) ? maximum(abs.(an64[mask_a] .- ComplexF64.(an32[mask_a])) ./ abs.(an64[mask_a])) : 0.0
    err_b = any(mask_b) ? maximum(abs.(bn64[mask_b] .- ComplexF64.(bn32[mask_b])) ./ abs.(bn64[mask_b])) : 0.0

    status = max(err_a, err_b) < 1e-3 ? "OK" :
             max(err_a, err_b) < 1e-1 ? "WARN" : "FAIL"

    @printf("  x=%6.0f (m=1.3-0.01i): aₙ err=%.2e  bₙ err=%.2e  [%s]\n",
            x, err_a, err_b, status)
end

# ── Part 5: Cross-section accuracy ──────────────────────────────────────

println("\n─── Part 5: Cross-section accuracy (C_ext, C_sca) ───\n")

@printf("  %8s  %12s  %12s  %12s  %12s  %8s\n",
        "x", "ΔC_ext/C_ext", "ΔC_sca/C_sca", "Δω̃/ω̃", "Δg/g", "status")
println("  " * "─" ^ 68)

for x in [1.0, 5.0, 10.0, 50.0, 100.0, 200.0, 300.0]
    n_max = Scattering.get_n_max(x)
    y = x * m
    nmx = round(Int, max(n_max, abs(y)) + 51)
    n_mu = 2 * n_max - 1
    k_val = 2π / 0.55

    # Float64
    an64 = zeros(ComplexF64, n_max)
    bn64 = zeros(ComplexF64, n_max)
    Dn64 = zeros(ComplexF64, nmx)
    Scattering.compute_mie_ab!(Float64(x), ComplexF64(m), an64, bn64, Dn64)
    n_ = Float64.(2 .* collect(1:n_max) .+ 1)
    C_sca64 = 2π / k_val^2 * dot(n_, abs2.(an64) .+ abs2.(bn64))
    C_ext64 = 2π / k_val^2 * dot(n_, real.(an64 .+ bn64))
    ω64 = C_sca64 / C_ext64

    # Float32
    an32 = zeros(ComplexF32, n_max)
    bn32 = zeros(ComplexF32, n_max)
    Dn32 = zeros(ComplexF32, nmx)
    Scattering.compute_mie_ab!(Float32(x), ComplexF32(m), an32, bn32, Dn32)
    n32 = Float32.(2 .* collect(1:n_max) .+ 1)
    C_sca32 = 2π / Float32(k_val)^2 * dot(n32, abs2.(an32) .+ abs2.(bn32))
    C_ext32 = 2π / Float32(k_val)^2 * dot(n32, real.(an32 .+ bn32))
    ω32 = C_sca32 / C_ext32

    err_ext = abs(Float64(C_ext32) - C_ext64) / max(abs(C_ext64), 1e-30)
    err_sca = abs(Float64(C_sca32) - C_sca64) / max(abs(C_sca64), 1e-30)
    err_ω = abs(Float64(ω32) - ω64) / max(abs(ω64), 1e-30)

    # Asymmetry factor from aₙ, bₙ (simplified, single particle)
    g64 = 0.0
    g32 = Float32(0.0)
    for n in 1:n_max-1
        fac1 = n * (n + 2) / (n + 1)
        fac2 = (2n + 1) / (n * (n + 1))
        g64 += fac1 * real(an64[n] * conj(an64[n+1]) + bn64[n] * conj(bn64[n+1])) +
               fac2 * real(an64[n] * conj(bn64[n]))
        g32 += Float32(fac1) * real(an32[n] * conj(an32[n+1]) + bn32[n] * conj(bn32[n+1])) +
               Float32(fac2) * real(an32[n] * conj(bn32[n]))
    end
    g64 *= 4π / (k_val^2 * C_sca64)
    g32 *= Float32(4π) / (Float32(k_val)^2 * C_sca32)
    err_g = abs(g64) > 1e-10 ? abs(Float64(g32) - g64) / abs(g64) : 0.0

    status = max(err_ext, err_sca, err_ω) < 1e-3 ? "OK" :
             max(err_ext, err_sca, err_ω) < 1e-1 ? "WARN" : "FAIL"

    @printf("  %8.0f  %12.2e  %12.2e  %12.2e  %12.2e  %8s\n",
            x, err_ext, err_sca, err_ω, err_g, status)
end

# ── Summary ─────────────────────────────────────────────────────────────

println("\n" * "=" ^ 72)
println("  Summary")
println("=" ^ 72)
println("""
  The downward recursion for Dₙ(y) with y = x·m is the critical path.
  For Float32 (7 significant digits), precision loss occurs when:
    - The recursion length nmx = max(n_max, |y|) + 51 grows large
    - |y| > ~100 starts showing significant Dₙ errors
    - These propagate into aₙ, bₙ and then into C_ext, C_sca, ω̃

  Recommendation: use Float32 only for x < ~50 (r < ~4 μm at λ=0.55 μm).
  For typical atmospheric aerosols (r_max=50, λ=0.3-2.5 μm), Float64 is
  essential for the Mie recursions, though downstream operations (Greek
  coefficient integration, RT kernel) could potentially use Float32.
""")
println("=" ^ 72)
