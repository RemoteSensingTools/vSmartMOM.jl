#!/usr/bin/env julia
"""
Quick test of Distributions.jl for log-normal size distributions
"""

using Distributions
using Printf

println("="^70)
println("Testing Distributions.jl for Log-Normal Size Distributions")
println("="^70)
println()

"""
    lognormal(r, params)

Log-normal size distribution using Distributions.jl

Parameters:
- r: Particle radius [μm]
- params: [N_total, r_med, sigma_g]

Returns:
- dN/dlogr: Number per log-radius bin [#/cm³]
"""
function lognormal(r::AbstractVector, p::AbstractVector)
    N_total, r_med, sigma_g = p
    # Create LogNormal distribution with μ=log(r_med), σ=log(sigma_g)
    dist = LogNormal(log(r_med), log(sigma_g))
    # LogNormal pdf is per unit r, convert to per log10(r):
    # dN/dlog10(r) = dN/dr × dr/dlog10(r) = dN/dr × r × ln(10)
    return N_total .* pdf.(dist, r) .* r .* log(10)
end

# Test parameters
N_total = 1000.0  # #/cm³
r_med = 0.1       # μm
sigma_g = 2.0     # geometric standard deviation

println("Test parameters:")
@printf("  N_total = %.1f #/cm³\n", N_total)
@printf("  r_med   = %.3f μm\n", r_med)
@printf("  σ_g     = %.2f\n", sigma_g)
println()

# Create distribution
dist = LogNormal(log(r_med), log(sigma_g))

# Test points
r_test = [0.01, 0.05, 0.1, 0.2, 0.5]

println("Log-normal distribution properties:")
@printf("  Mean radius: %.4f μm\n", mean(dist))
@printf("  Median radius: %.4f μm\n", median(dist))
@printf("  Mode radius: %.4f μm\n", mode(dist))
@printf("  Std deviation: %.4f\n", std(dist))
println()

println("PDF at test points:")
for r in r_test
    dN = lognormal([r], [N_total, r_med, sigma_g])[1]
    @printf("  r = %.3f μm: dN/dlogr = %.2f #/cm³\n", r, dN)
end

println()

# Verify manual calculation matches Distributions.jl
println("Verification (manual vs Distributions.jl):")
r_verify = 0.1

# Original manual formula (per natural log):
manual_ln = (N_total / (sqrt(2π) * log(sigma_g))) * 
            exp(-0.5 * (log(r_verify / r_med) / log(sigma_g))^2)

# Convert to per log10:
manual = manual_ln * log(10)

with_dist = lognormal([r_verify], [N_total, r_med, sigma_g])[1]

@printf("  At r = %.3f μm:\n", r_verify)
@printf("    Manual formula (ln): %.2f #/cm³\n", manual_ln)
@printf("    Manual formula (log10): %.2f #/cm³\n", manual)
@printf("    Distributions.jl: %.2f #/cm³\n", with_dist)
@printf("    Difference: %.2e (%.2f%%)\n", abs(manual - with_dist), 
        abs(manual - with_dist) / manual * 100)

println()
println("="^70)
println("✓ Test complete!")
println("="^70)
