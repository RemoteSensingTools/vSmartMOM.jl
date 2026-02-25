#!/usr/bin/env julia
# ==========================================================================
# Mie Code Performance Benchmark
# ==========================================================================
# Compares forward, analytic-linearized, and ForwardDiff-AD Mie performance.
# Isolates bottlenecks: per-radius allocation, S₁S₂ computation, Greek coefs.
# ==========================================================================

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", ".."))

using vSmartMOM
using vSmartMOM.Scattering
using Distributions
using Printf
using Statistics

# ── helpers ──────────────────────────────────────────────────────────────────
function build_model(; λ = 0.55, nᵣ = 1.3, nᵢ = 1e-8,
                       r_median = 0.3, σ_g = 1.5,
                       r_max = 30.0, nquad_radius = 500,
                       l_trunc = 20, Δ_angle = 2.0)
    sd = LogNormal(log(r_median), log(σ_g))
    aerosol = Aerosol(sd, nᵣ, nᵢ)
    pol = Stokes_IQU{Float64}()
    trunc = δBGE(l_trunc, Δ_angle)
    make_mie_model(NAI2(), aerosol, λ, pol, trunc, r_max, nquad_radius)
end

function timed_run(f, label; nruns = 3)
    # warmup
    f()
    GC.gc()

    times = Float64[]
    allocs = Int[]
    for _ in 1:nruns
        GC.gc()
        stats = @timed f()
        push!(times, stats.time)
        push!(allocs, Int(stats.bytes))
    end
    t_med = median(times)
    a_med = median(allocs)
    @printf("  %-28s  %8.3f s  %8.1f MB alloc  (n=%d)\n",
            label, t_med, a_med / 1e6, nruns)
    return t_med, a_med
end

# ── header ───────────────────────────────────────────────────────────────────
println("=" ^ 72)
println("  Mie Code Performance Benchmark")
println("=" ^ 72)

# ==========================================================================
# Part 1: Sweep nquad_radius and r_max
# ==========================================================================
println("\n─── Part 1: Scaling with nquad_radius and r_max ───\n")
println("  λ = 0.55 μm, nᵣ = 1.3, nᵢ = 1e-8, LogNormal(0.3, 1.5)")

results = []

for r_max in [10.0, 30.0, 50.0]
    for nq in [100, 500, 1000]
        model = build_model(r_max = r_max, nquad_radius = nq)
        x_max = 2π * r_max / 0.55
        n_max_est = Scattering.get_n_max(x_max)
        n_mu_est = 2 * n_max_est - 1
        @printf("\n  r_max=%4.0f  nquad=%4d  →  x_max≈%.0f  n_max≈%d  n_mu≈%d\n",
                r_max, nq, x_max, n_max_est, n_mu_est)

        t_fwd, a_fwd = timed_run("Forward") do
            compute_aerosol_optical_properties(model)
        end

        t_lin, a_lin = timed_run("Linearized (analytic)") do
            compute_aerosol_optical_properties(LinMode(), model)
        end

        t_ad, a_ad = timed_run("ForwardDiff AD") do
            compute_aerosol_optical_properties(model; autodiff = true)
        end

        push!(results, (; r_max, nq, t_fwd, t_lin, t_ad, a_fwd, a_lin, a_ad))
    end
end

println("\n─── Summary table ───\n")
@printf("  %5s  %5s  %10s  %10s  %10s  %8s  %8s\n",
        "r_max", "nquad", "Fwd (s)", "Lin (s)", "AD (s)", "Lin/Fwd", "AD/Fwd")
println("  " * "─" ^ 65)
for r in results
    @printf("  %5.0f  %5d  %10.4f  %10.4f  %10.4f  %8.2fx  %8.2fx\n",
            r.r_max, r.nq, r.t_fwd, r.t_lin, r.t_ad,
            r.t_lin / r.t_fwd, r.t_ad / r.t_fwd)
end

# ==========================================================================
# Part 2: Isolate individual bottlenecks at a single config
# ==========================================================================
println("\n─── Part 2: Bottleneck isolation (r_max=30, nquad=500) ───\n")

model = build_model(r_max = 30.0, nquad_radius = 500)
(; aerosol, λ, r_max, nquad_radius) = model
(; size_distribution, nᵣ, nᵢ) = aerosol

k = 2π / λ
r, wᵣ = Scattering.gauleg(nquad_radius, 0.0, r_max; norm = false)
x_sp = k .* r
n_max_global = Scattering.get_n_max(maximum(x_sp))
n_mu = 2 * n_max_global - 1
μ, w_μ = Scattering.gausslegendre(n_mu)
leg_π, leg_τ = Scattering.compute_mie_π_τ(μ, n_max_global)

# --- 2a: per-radius compute_mie_ab! timing ---
begin
    # Pre-allocate for largest size
    y_max = maximum(x_sp) * abs(nᵣ - nᵢ * im)
    nmx_max = round(Int, max(n_max_global, abs(y_max)) + 51)
    an_buf = zeros(ComplexF64, n_max_global)
    bn_buf = zeros(ComplexF64, n_max_global)
    Dn_buf = zeros(ComplexF64, nmx_max)

    function bench_mie_ab_loop()
        for i in 1:length(x_sp)
            fill!(an_buf, 0)
            fill!(bn_buf, 0)
            fill!(Dn_buf, 0)
            nm_i = Scattering.get_n_max(x_sp[i])
            Scattering.compute_mie_ab!(x_sp[i], nᵣ - nᵢ * im,
                                        view(an_buf, 1:nm_i),
                                        view(bn_buf, 1:nm_i),
                                        Dn_buf)
        end
    end
    timed_run(bench_mie_ab_loop, "compute_mie_ab! (all radii)"; nruns = 5)
end

# --- 2b: per-radius compute_mie_S₁S₂! timing ---
begin
    S₁ = zeros(ComplexF64, n_mu, nquad_radius)
    S₂ = zeros(ComplexF64, n_mu, nquad_radius)

    function bench_S1S2_loop()
        fill!(S₁, 0)
        fill!(S₂, 0)
        for i in 1:length(x_sp)
            fill!(an_buf, 0); fill!(bn_buf, 0); fill!(Dn_buf, 0)
            nm_i = Scattering.get_n_max(x_sp[i])
            Scattering.compute_mie_ab!(x_sp[i], nᵣ - nᵢ * im,
                                        view(an_buf, 1:nm_i),
                                        view(bn_buf, 1:nm_i),
                                        Dn_buf)
            Scattering.compute_mie_S₁S₂!(view(an_buf, 1:nm_i),
                                           view(bn_buf, 1:nm_i),
                                           leg_π, leg_τ,
                                           view(S₁, :, i),
                                           view(S₂, :, i))
        end
    end
    timed_run(bench_S1S2_loop, "ab! + S₁S₂! (all radii)"; nruns = 5)
end

# --- 2c: Greek coefficient loop timing ---
begin
    f₁₁ = rand(n_mu)
    f₃₃ = rand(n_mu)
    f₁₂ = rand(n_mu)
    f₃₄ = rand(n_mu)
    P, P², R², T² = Scattering.compute_legendre_poly(μ, n_mu)
    α_ = zeros(n_mu)
    β_ = zeros(n_mu)
    γ_ = zeros(n_mu)
    δ_ = zeros(n_mu)
    ϵ_ = zeros(n_mu)
    ζ_ = zeros(n_mu)

    function bench_greek()
        for l = 0:n_mu - 1
            fac = l ≥ 2 ? (2l + 1) / 2 * sqrt(1.0 / ((l - 1) * l * (l + 1) * (l + 2))) : 0.0
            δ_[l + 1] = (2l + 1) / 2 * (w_μ' * (f₃₃ .* P[:, l + 1]))
            β_[l + 1] = (2l + 1) / 2 * (w_μ' * (f₁₁ .* P[:, l + 1]))
            γ_[l + 1] = fac * (w_μ' * (f₁₂ .* P²[:, l + 1]))
            ϵ_[l + 1] = fac * (w_μ' * (f₃₄ .* P²[:, l + 1]))
            ζ_[l + 1] = fac * (w_μ' * (f₃₃ .* R²[:, l + 1] .+ f₁₁ .* T²[:, l + 1]))
            α_[l + 1] = fac * (w_μ' * (f₁₁ .* R²[:, l + 1] .+ f₃₃ .* T²[:, l + 1]))
        end
    end
    timed_run(bench_greek, "Greek coef loop (n_mu=$n_mu)"; nruns = 5)
end

# --- 2d: Allocation measurement in original code ---
println("\n  Allocation analysis (single full call):")

stats_fwd = @timed compute_aerosol_optical_properties(model)
@printf("    Forward:     %8.3f s  %8.1f MB  %d allocs\n",
        stats_fwd.time, stats_fwd.bytes / 1e6, 0)  # alloc count not in @timed

stats_lin = @timed compute_aerosol_optical_properties(LinMode(), model)
@printf("    Linearized:  %8.3f s  %8.1f MB  %d allocs\n",
        stats_lin.time, stats_lin.bytes / 1e6, 0)

stats_ad = @timed compute_aerosol_optical_properties(model; autodiff = true)
@printf("    ForwardDiff: %8.3f s  %8.1f MB  %d allocs\n",
        stats_ad.time, stats_ad.bytes / 1e6, 0)

# ==========================================================================
# Part 3: Accuracy cross-check (AD vs analytic linearized)
# ==========================================================================
println("\n─── Part 3: AD vs Analytic Linearized accuracy ───\n")

model_check = build_model(r_max = 30.0, nquad_radius = 500)
aero_fwd = compute_aerosol_optical_properties(model_check)
aero_ad = compute_aerosol_optical_properties(model_check; autodiff = true)
aero_an, lin_an = compute_aerosol_optical_properties(LinMode(), model_check)

# AD derivs shape: (6*greek_length+2, 4), columns = [r_m, σ, nᵣ, nᵢ]
ad_derivs = aero_ad.derivs
L = length(aero_fwd.greek_coefs.β)
param_names = ["r_m", "σ", "nᵣ", "nᵢ"]

println("  Forward values match:")
@printf("    ω̃: AD=%.8f  Fwd=%.8f  |diff|=%.2e\n",
        aero_ad.ω̃, aero_fwd.ω̃, abs(aero_ad.ω̃ - aero_fwd.ω̃))
@printf("    k:  AD=%.8e  Fwd=%.8e  |diff|=%.2e\n",
        aero_ad.k, aero_fwd.k, abs(aero_ad.k - aero_fwd.k))

# Compare analytic derivatives (lin) with AD derivatives
# lin_an has: lin_greek_coefs (α̇,β̇,...) with shape (4, L) for [nᵣ,nᵢ,μ,σ]
# AD has columns [r_m, σ, nᵣ, nᵢ]  
# Mapping: AD col 3 (nᵣ) ↔ lin row 1, AD col 4 (nᵢ) ↔ lin row 2

println("\n  Derivative comparison (analytic lin vs AD, for nᵣ and nᵢ):")
for (ad_col, lin_row, name) in [(3, 1, "nᵣ"), (4, 2, "nᵢ")]
    ad_dk = ad_derivs[end, ad_col]
    an_dk = lin_an.k̇[lin_row]
    rdiff_k = abs(ad_dk) > 1e-15 ? abs(ad_dk - an_dk) / abs(ad_dk) : abs(ad_dk - an_dk)
    @printf("    d(k)/d(%s):  AD=%.6e  Analytic=%.6e  rel_diff=%.2e\n",
            name, ad_dk, an_dk, rdiff_k)

    ad_dw = ad_derivs[end - 1, ad_col]
    an_dw = lin_an.ω̃̇[lin_row]
    rdiff_w = abs(ad_dw) > 1e-15 ? abs(ad_dw - an_dw) / abs(ad_dw) : abs(ad_dw - an_dw)
    @printf("    d(ω̃)/d(%s): AD=%.6e  Analytic=%.6e  rel_diff=%.2e\n",
            name, ad_dw, an_dw, rdiff_w)

    ad_beta = ad_derivs[L+1:2L, ad_col]
    an_beta = lin_an.lin_greek_coefs.β̇[lin_row, :]
    mask = abs.(ad_beta) .> 1e-15
    if any(mask)
        re = maximum(abs.(ad_beta[mask] .- an_beta[mask]) ./ abs.(ad_beta[mask]))
        @printf("    d(β)/d(%s):  max_rel_err=%.2e over %d nonzero coeffs\n",
                name, re, sum(mask))
    end
end

println("\n" * "=" ^ 72)
println("  Benchmark complete")
println("=" ^ 72)
