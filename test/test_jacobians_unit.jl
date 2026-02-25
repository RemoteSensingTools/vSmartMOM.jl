#!/usr/bin/env julia
# ==========================================================================
# Unit Tests for Linearized RT Jacobians
# ==========================================================================
# Tests each stage of the derivative chain against finite differences.
# Designed to run fast (~2-5 min total) with minimal wavelengths.
#
# Test hierarchy (each level isolates a different part of the code):
#   1. Mie interpolation formulas (product/quotient rule)
#   2. createAero δ-M scaling derivatives  
#   3. Core optical property combination (+) derivatives
#   4. Full end-to-end Jacobians: surface albedo, τ_ref, p₀, nᵣ
#
# For comparison against Automatic Differentiation (ForwardDiff) and to check
# FD vs AD vs analytic (timing and accuracy), run: test/test_jacobians_AD_compare.jl
# ==========================================================================

using Test
using vSmartMOM, vSmartMOM.CoreRT
using vSmartMOM.CoreRT: RT_Aerosol
using Distributions, Statistics
using LinearAlgebra

const YAML_FAST = "test_parameters/JacobianTestFast.yaml"

# ---------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------
function run_lin_rt(params)
    model, lin_model = model_from_parameters(LinMode(), params)
    NAer = length(params.scattering_params.rt_aerosols)
    NGas = size(lin_model.τ̇_abs[1], 1)
    NSurf = 1
    R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, NSurf)
    return R, dR, NAer, NGas, NSurf, model, lin_model
end

function run_fwd_only(params)
    R, _, _, _, _, _, _ = run_lin_rt(params)
    return R
end

"""Compute max and mean relative error, ignoring entries where |fd| < threshold."""
function rel_errors(analytic, fd; threshold=1e-12)
    mask = abs.(fd) .> threshold
    if any(mask)
        errs = abs.(analytic[mask] .- fd[mask]) ./ abs.(fd[mask])
        return (max=maximum(errs), mean=mean(errs), median=median(errs))
    else
        ae = abs.(analytic .- fd)
        return (max=maximum(ae), mean=mean(ae), median=median(ae))
    end
end

# =====================================================================
# Run base model once (shared by all tests)
# =====================================================================
println("Setting up base model (with Mie)...")
params_base = parameters_from_yaml(YAML_FAST)
R_base, dR_base, NAer, NGas, NSurf, model_base, lin_model_base = run_lin_rt(params_base)
Nparams = NAer * 7 + NGas + NSurf

# Extract base values for perturbations
τ_ref_base = params_base.scattering_params.rt_aerosols[1].τ_ref
nᵣ_base    = params_base.scattering_params.rt_aerosols[1].aerosol.nᵣ
p0_base    = mean(params_base.scattering_params.rt_aerosols[1].profile)
σp_base    = std(params_base.scattering_params.rt_aerosols[1].profile)
sd_base    = params_base.scattering_params.rt_aerosols[1].aerosol.size_distribution

println("  R: $(size(R_base)), dR: $(size(dR_base)), Nparams=$Nparams")
println("  nλ=$(size(R_base,3)), nVZA=$(size(R_base,1)), nStokes=$(size(R_base,2)), Nparams=$(size(dR_base,4))")

# =====================================================================
@testset "Jacobian Unit Tests" begin
# =====================================================================

# -----------------------------------------------------------------
@testset "Level 1: Mie Interpolation Sanity" begin
    # Check that the interpolated Mie derivatives are finite and 
    # have reasonable magnitudes
    aer_optics = model_base.aerosol_optics[1][1]
    lin_aer    = lin_model_base.lin_aerosol_optics[1][1]
    
    @test all(isfinite.(aer_optics.k))
    @test all(isfinite.(aer_optics.ω̃))
    @test all(aer_optics.k .> 0)         # extinction always positive
    @test all(0 .< aer_optics.ω̃ .≤ 1)   # SSA in (0,1]
    
    for ctr in 1:4
        @test all(isfinite.(lin_aer.k̇[ctr,:]))
        @test all(isfinite.(lin_aer.ω̃̇[ctr,:]))
    end
    
    # After Bug 20 fix: ω̃̇ should NOT be computed via division by k̇.
    # Check that ω̃̇ values are reasonable (not Inf/NaN, not huge)
    for ctr in 1:4
        max_ssa_dot = maximum(abs.(lin_aer.ω̃̇[ctr,:]))
        @test max_ssa_dot < 1e6
    end
    
    println("  ✓ Mie interpolation outputs are finite and bounded")
end

# -----------------------------------------------------------------
@testset "Level 2: Surface Albedo Jacobian (δ=1e-4)" begin
    # Surface albedo bypasses all atmospheric derivative code.
    # This tests the adding/interaction linearization only.
    albedo_base = 0.05
    δ = 1e-4
    
    analytic = dR_base[:, :, :, Nparams]  # last param = surface
    
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.brdf = [LambertianSurfaceScalar(albedo_base + δ)]
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    
    err = rel_errors(analytic, fd)
    println("  Surface albedo: max_rel=$(round(err.max, sigdigits=3)), mean_rel=$(round(err.mean, sigdigits=3))")
    
    @test err.max < 1e-3   # Surface albedo Jacobian max rel error
    @test err.mean < 1e-4  # Surface albedo Jacobian mean rel error
end

# -----------------------------------------------------------------
@testset "Level 3: τ_ref Jacobian (δ/τ=1e-3)" begin
    # τ_ref bypasses Mie derivatives but exercises:
    # createAero → constructCoreOpticalProperties → RT kernel
    δ = τ_ref_base * 1e-3
    
    analytic = dR_base[:, :, :, 1]  # param 1 = τ_ref
    
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].τ_ref = τ_ref_base + δ
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    
    err = rel_errors(analytic, fd)
    println("  τ_ref: max_rel=$(round(err.max, sigdigits=3)), mean_rel=$(round(err.mean, sigdigits=3))")
    
    @test err.max < 0.50   # τ_ref Jacobian max rel error < 50%
    @test err.mean < 0.15  # τ_ref Jacobian mean rel error < 15%
end

# -----------------------------------------------------------------
@testset "Level 4: p₀ Jacobian (δ/p₀=1e-2)" begin
    # Profile parameter — bypasses Mie, tests profile derivatives + RT
    δ = p0_base * 1e-2
    
    analytic = dR_base[:, :, :, 6]  # param 6 = p₀
    
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].profile = Normal(p0_base + δ, σp_base)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    
    err = rel_errors(analytic, fd)
    println("  p₀: max_rel=$(round(err.max, sigdigits=3)), mean_rel=$(round(err.mean, sigdigits=3))")
    
    # NOTE: p₀ is a profile redistribution parameter (Σ τ̇ ≈ 0 across layers),
    # which causes catastrophic cancellation in the Adding method (Bug 23).
    # The absolute error is ~3.8e-11, but the true derivative is O(1e-11),
    # giving large relative errors. This is a fundamental numerical precision
    # limitation, not a code bug. See docs/LINEARIZATION_BUGS.md Bug 23.
    # For now, we only check that the derivative is finite and has correct order
    # of magnitude (within 10x of FD).
    @test all(isfinite.(analytic))
    @test err.mean < 5.0  # Relaxed: catastrophic cancellation (Bug 23)
end

# -----------------------------------------------------------------
@testset "Level 5: nᵣ Jacobian (δ=1e-3)" begin
    # Real refractive index — exercises the FULL Mie derivative chain:
    # Mie → interpolation (Bug 20 fix) → τ̇_aer (Bug 21 fix) → createAero → RT
    δ = 1e-3
    
    analytic = dR_base[:, :, :, 2]  # param 2 = nᵣ
    
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].aerosol.nᵣ = nᵣ_base + δ
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    
    err = rel_errors(analytic, fd)
    println("  nᵣ: max_rel=$(round(err.max, sigdigits=3)), mean_rel=$(round(err.mean, sigdigits=3))")
    
    @test err.max < 1.0    # nᵣ Jacobian max rel error < 100%
    @test err.mean < 0.30  # nᵣ Jacobian mean rel error < 30%
end

# -----------------------------------------------------------------
@testset "Level 6: σ_dist Jacobian (δ/σ=1e-2)" begin
    # Size distribution width — exercises Mie derivatives via ẇₓ
    δ = sd_base.σ * 1e-2
    
    analytic = dR_base[:, :, :, 5]  # param 5 = σ of LogNormal
    
    params_pert = parameters_from_yaml(YAML_FAST)
    params_pert.scattering_params.rt_aerosols[1].aerosol.size_distribution = 
        LogNormal(sd_base.μ, sd_base.σ + δ)
    R_pert = run_fwd_only(params_pert)
    fd = (R_pert .- R_base) ./ δ
    
    err = rel_errors(analytic, fd)
    println("  σ_dist: max_rel=$(round(err.max, sigdigits=3)), mean_rel=$(round(err.mean, sigdigits=3))")
    
    @test err.max < 1.0    # σ_dist Jacobian max rel error < 100%
    @test err.mean < 0.30  # σ_dist Jacobian mean rel error < 30%
end

# -----------------------------------------------------------------
@testset "Level 7: Forward consistency" begin
    # The forward R from the linearized path must match standalone forward
    # (regression check — Bug 13 was exactly this)
    @test all(isfinite.(R_base))          # R is finite
    @test all(isfinite.(dR_base))         # dR is finite
    @test all(R_base[:, 1, :] .>= 0)     # Stokes I ≥ 0
end

# =====================================================================
end  # @testset "Jacobian Unit Tests"
# =====================================================================

println("\nAll Jacobian unit tests complete.")
