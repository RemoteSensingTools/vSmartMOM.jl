#!/usr/bin/env julia
# ==========================================================================
# Mie Code GPU Feasibility Benchmark
# ==========================================================================
# Tests GPU-accelerated Mie coefficient computation using CUDA kernels.
# The per-radius Mie computation is embarrassingly parallel -- each radius
# point computes its own Dₙ recursion + aₙ, bₙ independently.
# ==========================================================================

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using vSmartMOM
using vSmartMOM.Scattering
using Distributions
using Printf
using Statistics
using LinearAlgebra
using CUDA

CUDA.allowscalar(true)

println("=" ^ 72)
println("  Mie GPU Feasibility Benchmark")
println("=" ^ 72)
println("  GPU: $(CUDA.name(CUDA.device()))")
println()

# ==========================================================================
# GPU kernel: batched compute_mie_ab!
# ==========================================================================
# Each thread handles one radius point. The recursion is serial per thread
# but independent across threads.

function gpu_compute_mie_ab_kernel!(
    an_out,     # (n_max_global, nquad)  Complex
    bn_out,     # (n_max_global, nquad)  Complex
    x_sp,       # (nquad,) size parameters
    m_re,       # real part of refractive index (scalar)
    m_im,       # imag part of refractive index (scalar)
    n_maxs,     # (nquad,) per-radius n_max values
    n_max_global, # max n_max across all radii
    nmx_global  # max nmx across all radii (for Dₙ workspace)
)
    nquad = length(x_sp)
    # Shared workspace for Dₙ is too large for shared mem; use global
    # Each thread gets a column in the pre-allocated Dₙ buffer
    Dn = CUDA.zeros(ComplexF64, nmx_global, nquad)

    # Launch kernel
    threads_per_block = 128
    nblocks = cld(nquad, threads_per_block)

    @cuda threads=threads_per_block blocks=nblocks _mie_ab_kernel!(
        an_out, bn_out, Dn, x_sp,
        m_re, m_im, n_maxs, n_max_global, nmx_global, nquad)

    CUDA.synchronize()
    return nothing
end

function _mie_ab_kernel!(an_out, bn_out, Dn, x_sp,
                         m_re, m_im, n_maxs, n_max_global, nmx_global, nquad)
    idx = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    if idx > nquad
        return nothing
    end

    x = x_sp[idx]
    n_max = n_maxs[idx]
    m = Complex(m_re, -m_im)
    y = x * m

    nmx = nmx_global  # use global max to keep it simple

    # Downward recursion for Dₙ
    for n = (nmx - 1):-1:1
        Dn[n, idx] = ((n + 1) / y) - (1.0 / (Dn[n + 1, idx] + (n + 1) / y))
    end

    # Bessel function recursion
    ψ₀ = cos(x)
    ψ₁ = sin(x)
    χ₀ = -sin(x)
    χ₁ = cos(x)
    ξ₁ = Complex(ψ₁, χ₁)

    for n = 1:n_max
        ψ = (2n - 1) * ψ₁ / x - ψ₀
        χ = (2n - 1) * χ₁ / x - χ₀
        ξ = Complex(ψ, χ)
        t_a = Dn[n, idx] / m + n / x
        t_b = Dn[n, idx] * m + n / x

        an_out[n, idx] = (t_a * ψ - ψ₁) / (t_a * ξ - ξ₁)
        bn_out[n, idx] = (t_b * ψ - ψ₁) / (t_b * ξ - ξ₁)

        ψ₀ = ψ₁
        ψ₁ = ψ
        χ₀ = χ₁
        χ₁ = χ
        ξ₁ = Complex(ψ₁, χ₁)
    end

    return nothing
end

# ==========================================================================
# CPU reference (using existing Scattering functions)
# ==========================================================================
function cpu_compute_all_ab(x_sp, nᵣ, nᵢ, n_max_global, nmx_global)
    nquad = length(x_sp)
    m = nᵣ - nᵢ * im
    an = zeros(ComplexF64, n_max_global, nquad)
    bn = zeros(ComplexF64, n_max_global, nquad)
    Dn_buf = zeros(ComplexF64, nmx_global)

    for i in 1:nquad
        fill!(Dn_buf, 0)
        nm_i = Scattering.get_n_max(x_sp[i])
        Scattering.compute_mie_ab!(x_sp[i], m,
                                    view(an, 1:nm_i, i),
                                    view(bn, 1:nm_i, i),
                                    Dn_buf)
    end
    return an, bn
end

# ==========================================================================
# Benchmark loop
# ==========================================================================

nᵣ = 1.3
nᵢ = 1e-8
λ = 0.55

println("─── Benchmark: CPU vs GPU compute_mie_ab! ───\n")
@printf("  %6s  %6s  %6s  %8s  %8s  %8s  %12s  %12s\n",
        "r_max", "nquad", "n_max", "CPU (s)", "GPU (s)", "speedup",
        "max|Δaₙ|rel", "max|Δbₙ|rel")
println("  " * "─" ^ 80)

for (r_max, nquad) in [(10.0, 100), (10.0, 500), (10.0, 1000),
                        (30.0, 100), (30.0, 500), (30.0, 1000), (30.0, 2500),
                        (50.0, 500), (50.0, 1000), (50.0, 2500)]
    k = 2π / λ
    r, _ = Scattering.gauleg(nquad, 0.0, r_max; norm = false)
    x_sp = k .* r

    n_max_per = [Scattering.get_n_max(x) for x in x_sp]
    n_max_global = maximum(n_max_per)
    nmx_global = round(Int, max(n_max_global,
                     maximum(x_sp) * abs(nᵣ - nᵢ * im)) + 51)

    # CPU benchmark
    cpu_compute_all_ab(x_sp, nᵣ, nᵢ, n_max_global, nmx_global)  # warmup
    GC.gc()
    t_cpu = @elapsed begin
        an_cpu, bn_cpu = cpu_compute_all_ab(x_sp, nᵣ, nᵢ, n_max_global, nmx_global)
    end

    # GPU benchmark
    x_sp_d = CuArray(x_sp)
    n_maxs_d = CuArray(Int32.(n_max_per))
    an_d = CUDA.zeros(ComplexF64, n_max_global, nquad)
    bn_d = CUDA.zeros(ComplexF64, n_max_global, nquad)

    # warmup
    gpu_compute_mie_ab_kernel!(an_d, bn_d, x_sp_d, nᵣ, nᵢ,
                                n_maxs_d, n_max_global, nmx_global)

    GC.gc(); CUDA.reclaim()
    CUDA.synchronize()
    t_gpu = @elapsed begin
        fill!(an_d, 0); fill!(bn_d, 0)
        gpu_compute_mie_ab_kernel!(an_d, bn_d, x_sp_d, nᵣ, nᵢ,
                                    n_maxs_d, n_max_global, nmx_global)
    end

    # Accuracy check
    an_gpu = Array(an_d)
    bn_gpu = Array(bn_d)

    mask_a = abs.(an_cpu) .> 1e-30
    mask_b = abs.(bn_cpu) .> 1e-30
    err_a = any(mask_a) ? maximum(abs.(an_cpu[mask_a] .- an_gpu[mask_a]) ./ abs.(an_cpu[mask_a])) : 0.0
    err_b = any(mask_b) ? maximum(abs.(bn_cpu[mask_b] .- bn_gpu[mask_b]) ./ abs.(bn_cpu[mask_b])) : 0.0

    speedup = t_cpu / t_gpu
    @printf("  %6.0f  %6d  %6d  %8.4f  %8.4f  %8.1fx  %12.2e  %12.2e\n",
            r_max, nquad, n_max_global, t_cpu, t_gpu, speedup, err_a, err_b)

    CUDA.unsafe_free!(x_sp_d)
    CUDA.unsafe_free!(n_maxs_d)
    CUDA.unsafe_free!(an_d)
    CUDA.unsafe_free!(bn_d)
end

# ==========================================================================
# GPU S₁S₂ computation
# ==========================================================================

println("\n─── Benchmark: GPU S₁S₂ computation ───\n")

function gpu_compute_S1S2!(S1, S2, an, bn, leg_pi, leg_tau, n_maxs, nquad)
    n_mu = size(leg_pi, 1)
    threads = 256
    blocks = cld(nquad * n_mu, threads)

    @cuda threads=threads blocks=blocks _S1S2_kernel!(
        S1, S2, an, bn, leg_pi, leg_tau, n_maxs, n_mu, nquad)
    CUDA.synchronize()
end

function _S1S2_kernel!(S1, S2, an, bn, leg_pi, leg_tau, n_maxs, n_mu, nquad)
    tid = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    total = n_mu * nquad
    if tid > total
        return nothing
    end

    imu = ((tid - 1) % n_mu) + 1
    iq  = ((tid - 1) ÷ n_mu) + 1

    nmax = n_maxs[iq]
    s1 = Complex(0.0, 0.0)
    s2 = Complex(0.0, 0.0)

    for l in 1:nmax
        fac = (2l + 1) / (l * (l + 1))
        s1 += fac * (an[l, iq] * leg_tau[imu, l] + bn[l, iq] * leg_pi[imu, l])
        s2 += fac * (an[l, iq] * leg_pi[imu, l] + bn[l, iq] * leg_tau[imu, l])
    end

    S1[imu, iq] = s1
    S2[imu, iq] = s2
    return nothing
end

# Test S₁S₂ at a single configuration
r_max_test = 30.0
nquad_test = 500

k = 2π / λ
r, _ = Scattering.gauleg(nquad_test, 0.0, r_max_test; norm = false)
x_sp = k .* r
n_max_per = [Scattering.get_n_max(x) for x in x_sp]
n_max_global = maximum(n_max_per)
n_mu = 2 * n_max_global - 1
nmx_global = round(Int, max(n_max_global,
                 maximum(x_sp) * abs(nᵣ - nᵢ * im)) + 51)

μ, w_μ = Scattering.gausslegendre(n_mu)
leg_π, leg_τ = Scattering.compute_mie_π_τ(μ, n_max_global)

# CPU reference: ab + S1S2
an_cpu, bn_cpu = cpu_compute_all_ab(x_sp, nᵣ, nᵢ, n_max_global, nmx_global)
S1_cpu = zeros(ComplexF64, n_mu, nquad_test)
S2_cpu = zeros(ComplexF64, n_mu, nquad_test)
for i in 1:nquad_test
    nm = n_max_per[i]
    Scattering.compute_mie_S₁S₂!(view(an_cpu, 1:nm, i),
                                   view(bn_cpu, 1:nm, i),
                                   leg_π, leg_τ,
                                   view(S1_cpu, :, i),
                                   view(S2_cpu, :, i))
end

# GPU: ab + S1S2
x_sp_d = CuArray(x_sp)
n_maxs_d = CuArray(Int32.(n_max_per))
an_d = CUDA.zeros(ComplexF64, n_max_global, nquad_test)
bn_d = CUDA.zeros(ComplexF64, n_max_global, nquad_test)
gpu_compute_mie_ab_kernel!(an_d, bn_d, x_sp_d, nᵣ, nᵢ,
                            n_maxs_d, n_max_global, nmx_global)

leg_π_d = CuArray(leg_π)
leg_τ_d = CuArray(leg_τ)
S1_d = CUDA.zeros(ComplexF64, n_mu, nquad_test)
S2_d = CUDA.zeros(ComplexF64, n_mu, nquad_test)

# warmup
gpu_compute_S1S2!(S1_d, S2_d, an_d, bn_d, leg_π_d, leg_τ_d, n_maxs_d, nquad_test)

GC.gc(); CUDA.reclaim(); CUDA.synchronize()

# Timed GPU S1S2
fill!(S1_d, 0); fill!(S2_d, 0)
t_gpu_s = @elapsed begin
    gpu_compute_S1S2!(S1_d, S2_d, an_d, bn_d, leg_π_d, leg_τ_d, n_maxs_d, nquad_test)
end

# Timed CPU S1S2 (just the S1S2 part, ab already computed)
t_cpu_s = @elapsed begin
    S1_tmp = zeros(ComplexF64, n_mu, nquad_test)
    S2_tmp = zeros(ComplexF64, n_mu, nquad_test)
    for i in 1:nquad_test
        nm = n_max_per[i]
        Scattering.compute_mie_S₁S₂!(view(an_cpu, 1:nm, i),
                                       view(bn_cpu, 1:nm, i),
                                       leg_π, leg_τ,
                                       view(S1_tmp, :, i),
                                       view(S2_tmp, :, i))
    end
end

S1_gpu = Array(S1_d)
S2_gpu = Array(S2_d)
mask_s1 = abs.(S1_cpu) .> 1e-20
err_s1 = any(mask_s1) ? maximum(abs.(S1_cpu[mask_s1] .- S1_gpu[mask_s1]) ./ abs.(S1_cpu[mask_s1])) : 0.0

@printf("  S₁S₂ (r_max=%.0f, nquad=%d, n_mu=%d):\n", r_max_test, nquad_test, n_mu)
@printf("    CPU: %.4f s\n", t_cpu_s)
@printf("    GPU: %.4f s  (%.1fx speedup)\n", t_gpu_s, t_cpu_s / t_gpu_s)
@printf("    Max S₁ relative error: %.2e\n", err_s1)

# ==========================================================================
# Full pipeline timing: CPU forward vs GPU (ab + S1S2 only)
# ==========================================================================
println("\n─── Full pipeline: CPU forward Mie vs GPU (ab+S₁S₂) ───\n")

for (r_max_t, nq_t) in [(10.0, 500), (30.0, 500), (30.0, 1000)]
    sd = LogNormal(log(0.3), log(1.5))
    aerosol = Aerosol(sd, nᵣ, nᵢ)
    pol = Stokes_IQU{Float64}()
    trunc = δBGE(20, 2.0)
    model = make_mie_model(NAI2(), aerosol, λ, pol, trunc, r_max_t, nq_t)

    # CPU full forward
    compute_aerosol_optical_properties(model)  # warmup
    GC.gc()
    t_cpu_full = @elapsed compute_aerosol_optical_properties(model)

    # GPU: just the ab + S1S2 part
    k_t = 2π / λ
    r_t, _ = Scattering.gauleg(nq_t, 0.0, r_max_t; norm = false)
    x_sp_t = k_t .* r_t
    nmp_t = [Scattering.get_n_max(x) for x in x_sp_t]
    nmg_t = maximum(nmp_t)
    nmu_t = 2 * nmg_t - 1
    nmx_t = round(Int, max(nmg_t, maximum(x_sp_t) * abs(nᵣ - nᵢ * im)) + 51)

    μ_t, _ = Scattering.gausslegendre(nmu_t)
    lp_t, lt_t = Scattering.compute_mie_π_τ(μ_t, nmg_t)

    x_d = CuArray(x_sp_t)
    nm_d = CuArray(Int32.(nmp_t))
    a_d = CUDA.zeros(ComplexF64, nmg_t, nq_t)
    b_d = CUDA.zeros(ComplexF64, nmg_t, nq_t)
    lp_d = CuArray(lp_t)
    lt_d = CuArray(lt_t)
    s1_d = CUDA.zeros(ComplexF64, nmu_t, nq_t)
    s2_d = CUDA.zeros(ComplexF64, nmu_t, nq_t)

    # warmup
    gpu_compute_mie_ab_kernel!(a_d, b_d, x_d, nᵣ, nᵢ, nm_d, nmg_t, nmx_t)
    gpu_compute_S1S2!(s1_d, s2_d, a_d, b_d, lp_d, lt_d, nm_d, nq_t)

    GC.gc(); CUDA.reclaim(); CUDA.synchronize()

    fill!(a_d, 0); fill!(b_d, 0); fill!(s1_d, 0); fill!(s2_d, 0)
    t_gpu_full = @elapsed begin
        gpu_compute_mie_ab_kernel!(a_d, b_d, x_d, nᵣ, nᵢ, nm_d, nmg_t, nmx_t)
        gpu_compute_S1S2!(s1_d, s2_d, a_d, b_d, lp_d, lt_d, nm_d, nq_t)
    end

    @printf("  r_max=%4.0f, nquad=%4d: CPU_full=%.3fs  GPU_ab+S₁S₂=%.4fs  (%.1fx)\n",
            r_max_t, nq_t, t_cpu_full, t_gpu_full, t_cpu_full / t_gpu_full)

    # cleanup
    for arr in [x_d, nm_d, a_d, b_d, lp_d, lt_d, s1_d, s2_d]
        CUDA.unsafe_free!(arr)
    end
end

# ==========================================================================
# GPU memory analysis
# ==========================================================================
println("\n─── GPU Memory Requirements ───\n")

for (r_max_t, nq_t) in [(30.0, 1000), (50.0, 2500)]
    k_t = 2π / λ
    x_max = k_t * r_max_t
    nmg = Scattering.get_n_max(x_max)
    nmu = 2 * nmg - 1
    nmx = round(Int, max(nmg, x_max * abs(nᵣ - nᵢ * im)) + 51)

    # Memory per array
    mem_Dn = nmx * nq_t * 16  # ComplexF64
    mem_an = nmg * nq_t * 16
    mem_bn = nmg * nq_t * 16
    mem_S1 = nmu * nq_t * 16
    mem_S2 = nmu * nq_t * 16
    mem_pi = nmu * nmg * 8
    mem_tau = nmu * nmg * 8
    total = mem_Dn + mem_an + mem_bn + mem_S1 + mem_S2 + mem_pi + mem_tau

    @printf("  r_max=%4.0f, nquad=%4d: n_max=%d, n_mu=%d\n", r_max_t, nq_t, nmg, nmu)
    @printf("    Dₙ buffer:  %6.1f MB\n", mem_Dn / 1e6)
    @printf("    aₙ, bₙ:     %6.1f MB each\n", mem_an / 1e6)
    @printf("    S₁, S₂:     %6.1f MB each\n", mem_S1 / 1e6)
    @printf("    π, τ:       %6.1f MB each\n", mem_pi / 1e6)
    @printf("    Total:      %6.1f MB (%.1f%% of GPU)\n",
            total / 1e6, total / CUDA.totalmem(CUDA.device()) * 100)
    println()
end

println("=" ^ 72)
println("  GPU Feasibility Benchmark Complete")
println("=" ^ 72)
