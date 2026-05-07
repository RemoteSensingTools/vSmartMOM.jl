using vSmartMOM
using vSmartMOM.Architectures
using vSmartMOM.Absorption
using vSmartMOM.Scattering
using vSmartMOM.CoreRT
using vSmartMOM.SolarModel
using Test
using DelimitedFiles
using Statistics
using LinearAlgebra
using ProgressMeter
using WignerSymbols
using Distributions
using JLD2

const _VSMARTMOM_TEST_ORIGINAL_CWD = pwd()
cd(@__DIR__)
try

# Core module tests
@testset "Absorption" begin include("test_Absorption.jl") end
@testset "Scattering" begin include("test_Scattering.jl") end
@testset "CoreRT" begin include("test_CoreRT.jl") end
@testset "Batched Kernels" begin include("test_batched_kernels.jl") end
@testset "SolarModel" begin include("test_SolarModel.jl") end

# Forward model tests (these require YAML parameter files + data)
@testset "Forward noRS" begin include("test_forward_noRS.jl") end

# Jacobian unit tests
@testset "Jacobian Unit" begin include("test_jacobians_unit.jl") end

# Type stability tests (no external data)
@testset "Type Stability" begin include("test_type_stability.jl") end

# Float32 consistency tests
@testset "Float32 Consistency" begin include("test_float32.jl") end

# Quality gates
@testset "Quality Gates" begin include("test_quality.jl") end

# Top-level public IO API
@testset "IO Exports" begin include("test_io_exports.jl") end

# Parameter parser regression tests
@testset "Parameter Parser" begin include("test_parameters_parser.jl") end

# IO validation regression tests
@testset "IO Validation" begin include("test_io_validation.jl") end

@testset "Canopy Surface" begin include("test_canopy.jl") end

# Phase A — `Nstreams` field on `QuadPoints` (count of nonzero weights;
# distinct from augmented `Nquad`). See docs/src/pages/conventions.md §6.
@testset "QuadPoints streams" begin include("test_quadpoints_streams.jl") end

# Phase-function truncation invariants (Sanghavi & Stephens 2015)
@testset "Truncation" begin include("test_truncation.jl") end

# Cox-Munk ocean surface tests
@testset "Cox-Munk Surface" begin include("test_coxmunk.jl") end

# GPU-specific tests (conditional on CUDA availability)
CUDA_AVAILABLE = try
    using CUDA
    CUDA.functional()
catch
    false
end

if CUDA_AVAILABLE
    @testset "Raman GPU" begin include("test_forward_raman_gpu.jl") end
end

# Phase 1b regression gate — RRS forward model vs frozen reference.
# Skipped by default on CPU-only machines (run takes ~3 min); set
# PHASE1B_CPU=1 or have CUDA available to actually run.
@testset "Phase 1b RRS regression" begin include("test_forward_raman_phase1b.jl") end

# Phase 1c single-scatter driver smoke test.
@testset "Phase 1c SS driver" begin include("test_forward_ss.jl") end

# Standalone exact single-scatter Phase 1.
@testset "StandaloneSS" begin include("StandaloneSS/runtests.jl") end

# Phase 1d — Aerosols module wire-in (close `using vSmartMOM.Aerosols` export gap).
@testset "Aerosols" begin include("test_Aerosols.jl") end

# Phase 1e — perturb_parameters utility ported from sanghavi.
@testset "Phase 1e perturb_parameters" begin include("test_perturb_parameters.jl") end

# Phase 3a — SIF injection + data loaders (Lambertian surface + sif_loader.jl).
# Loader subtests auto-skip when their fixtures
# (`src/SIF_emission/{sif-spectra.csv, ficus_refl_600to800nm.dat}`) are absent;
# the inject_surface_SIF! kernel coverage runs unconditionally.
@testset "Phase 3a SIF" begin include("test_sif.jl") end

# Phase 6 — sanghavi test/benchmarks/*.jl script ports (parse + light-unit).
@testset "Phase 6 script ports" begin include("test_phase6_ports.jl") end

# v0.6 source-term refactor — AbstractSource vocabulary, SolarBeam, BlackbodySource,
# SurfaceSIF, surface_source_contribute! double-dispatch, prepared_sources flow.
# Includes the full Phase 1 → 5b regression assertions and end-to-end bit-equality
# checks via small CPU rt_run scenarios on PureRayleighParameters.
@testset "Sources (v0.6)" begin include("test_sources.jl") end

# VLIDORT 2.8.3 baseline validation suite — Siewert 2000 PROBLEM_IIA Stokes-I,
# solar_tester scalar (Task 1), solar_tester vector (Task 1, IQU). Reference
# tables are committed under test/vlidort_baseline/reference_data/; no PyVLIDORT
# / Fortran runtime needed. The included file already wraps its own
# `@testset "VLIDORT baseline"`, so we don't double-nest.
include("vlidort_baseline/runtests.jl")

finally
    cd(_VSMARTMOM_TEST_ORIGINAL_CWD)
end
