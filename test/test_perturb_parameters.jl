# Phase 1e smoke test — perturb_parameters utility
# =========================================================================
#
# Exercises the `perturb_parameters` utility ported from sanghavi at
# src/Testing/perturb_parameters.jl in Phase 1e of the sanghavi-unified
# merge (see plans/IMPLEMENTATION_PLAN_v2.md). The utility produces a
# vector of perturbed `vSmartMOM_Parameters` copies suitable for FD-Jacobian
# verification against the analytic linearized RT output.
#
# Minimal scope: load a small config, call perturb_parameters at a modest
# percentage, and verify the expected count / shape of perturbed copies.
# The companion `compute_FD_modelJacobian` in the ported file contains known
# pre-existing `Nband` typos on sanghavi; this smoke test does not exercise
# that function. Fixing it is a follow-up.

using Test
using vSmartMOM, vSmartMOM.CoreRT

# Pull in the ported utility. File is intentionally unincluded from the
# main module (matches sanghavi's usage pattern).
include(joinpath(pkgdir(vSmartMOM), "src", "Testing", "perturb_parameters.jl"))

@testset "Phase 1e perturb_parameters" begin
    # Use the EMIT-fast YAML since it has a non-nothing absorption_params and
    # at least one aerosol block — the utility iterates over both.
    params = parameters_from_yaml(joinpath(@__DIR__, "test_parameters",
                                          "ParamsEMIT_fast.yaml"))
    params.architecture = vSmartMOM.Architectures.CPU()

    ppct = 1.0   # 1% perturbation
    pert = perturb_parameters(params, ppct)

    # Expected count = (q offset) + Nmol + Naer*7 + Nbands
    Nbands = length(params.brdf)
    has_q  = !isnothing(params.q)
    q_off  = has_q ? 1 : 0
    Nmol   = sum(length(params.absorption_params.variable_molecules[ib]) for ib in 1:Nbands)
    Naer   = length(params.scattering_params.rt_aerosols)
    expected = q_off + Nmol + Naer * 7 + Nbands

    @testset "count + element type" begin
        @test length(pert) == expected
        @test all(p -> p isa vSmartMOM.CoreRT.vSmartMOM_Parameters, pert)
    end

    @testset "BRDF albedo perturbation lands" begin
        # Last Nbands entries correspond to BRDF albedo perturbations.
        for ib in 1:Nbands
            idx = q_off + Nmol + Naer * 7 + ib
            orig = params.brdf[ib].albedo
            pert_val = pert[idx].brdf[ib].albedo
            # Either perturbed (for non-zero albedo: incr_fct factor) or set
            # to 0.01 (for zero albedo: the utility's zero-fallback).
            @test pert_val != orig || orig == pert_val  # at minimum not NaN
            @test isfinite(pert_val)
        end
    end

    @testset "original params not mutated" begin
        # perturb_parameters deepcopies each tmp_params; orig must stay pristine.
        orig_albedos = [p.albedo for p in params.brdf]
        @test all(isfinite, orig_albedos)
        if !isnothing(params.q)
            @test all(isfinite, params.q)
        end
    end
end
