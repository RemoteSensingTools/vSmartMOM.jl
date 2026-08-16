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

# Caller-defined-node bulk Mie API (Sol plan phase 1) -- caller-node
# compute_aerosol_optical_properties_nodes/_gpu; delegation refactor of the
# log-normal NAI2 path via the shared _nai2_bulk_optics core.
@testset "Mie caller-node API" begin include("test_mie_nodes.jl") end

# Batched caller-node Mie seam (proposals/batched_mie_nodes_seam.md v2) --
# exact-nmax grouped GPU batching + reusable single-owner workspace.
@testset "Mie caller-node batched API" begin include("test_mie_nodes_batched.jl") end
@testset "CoreRT" begin include("test_CoreRT.jl") end
@testset "Batched Kernels" begin include("test_batched_kernels.jl") end
@testset "SolarModel" begin include("test_SolarModel.jl") end

# Forward model tests (these require YAML parameter files + data)
@testset "Forward noRS" begin include("test_forward_noRS.jl") end

# Molecular Rayleigh/Cabannes/RRS optical-property identities.
@testset "Rayleigh/Cabannes Raman" begin include("test_ray_cab_raman.jl") end

# Jacobian unit tests
@testset "Jacobian Unit" begin include("test_jacobians_unit.jl") end

# Type stability tests (no external data)
@testset "Type Stability" begin include("test_type_stability.jl") end

# Float32 consistency tests
@testset "Float32 Consistency" begin include("test_float32.jl") end

# Quality gates
@testset "Quality Gates" begin include("test_quality.jl") end

# Package health (Aqua) and static analysis snapshot (JET)
@testset "Aqua" begin include("test_aqua.jl") end
@testset "JET"  begin include("test_jet.jl")  end

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

# Phase B — forward and lin paths must derive the same per-band Fourier
# loop bound. Regression for the previous lin-only precedence bug at
# even `l_max`; both paths now flow through `_derive_m_max_bands`.
@testset "Forward/lin m_max parity" begin include("test_lin_forward_loop_parity.jl") end

# Phase C — `component_m_max(c, ctx)` trait dispatch + flag-gated
# aggregator. **Default flag is `true`** in v2.1: Cox-Munk / RossLi /
# RPV / canopy run to their full `user_l_cap`. The flag exists as an
# escape hatch (`use_component_traits = false`) so the legacy
# half-truncated aggregator can still be selected for byte-equality
# regression.
@testset "Component m_max traits" begin include("test_component_m_max.jl") end

# Phase D — schema-doc invariants. Ensures every top-level YAML/TOML
# block has a per-block doc page and the JSON Schema covers the v0.7
# resolution knobs.
@testset "Schema docs" begin include("test_schema_docs.jl") end

# Phase-function truncation invariants (Sanghavi & Stephens 2015)
@testset "Truncation" begin include("test_truncation.jl") end

# VS Raman path smoke tests — guards the previously-broken VS_0to1/VS_1to0
# call chain (undefined compute_optical_Rayl!/compute_ϖ_Cabannes! calls).
@testset "Inelastic VS smoke" begin include("test_inelastic_vs_smoke.jl") end

# ForwardDiff hybrid-AD path — guards the Dual-safe Dₙ recursion in
# compute_mie_ab! (Float64-stabilized for plain floats, native for Duals)
# and the AD albedo Jacobian. Runtime ~2 min.
@testset "Hybrid AD" begin include("test_hybrid_ad.jl") end

# Cox-Munk ocean surface tests
@testset "Cox-Munk Surface" begin include("test_coxmunk.jl") end

# Lambertian surface scaffold — cross-flavor consistency (Scalar/Legendre/
# Spline/Spectrum share one create_surface_layer!), m>0 conventions, and the
# LambertianSurfaceSpectrum end-to-end regression.
@testset "Lambertian surfaces" begin include("test_lambertian_surfaces.jl") end

# Atmosphere/surface split — rt_run_atmosphere / rt_run_surface /
# rt_run_multi_surface. Bit-exactness regression against monolithic rt_run
# (same kernels, same inputs); see proposals/surface_split_albedo_sweep.md.
@testset "Surface split" begin include("test_surface_split.jl") end

# Analytic Lambertian-albedo closure (O(1) per albedo) over a cached
# atmosphere, + the :slim AtmosphereRTCache mode; see
# proposals/surface_split_albedo_sweep.md §3-4.
@testset "Lambertian closure" begin include("test_lambertian_closure.jl") end

# SZA × view-pair × BRDF scenario sweeps: `remake_geometry`'s bit-exactness
# vs a full `model_from_parameters` rebuild, plus `run_sweep`'s output
# against a monolithic `rt_run` reference; see
# proposals/surface_split_albedo_sweep.md §6-7 (PR 3).
@testset "Scenario sweep" begin include("test_scenario_sweep.jl") end

# NOTE: GPU/Metal tests live in test/local/gpu/ and run via
# `julia --project=test test/local/runtests.jl` (the local-only suite); they are
# NOT part of this CI suite. Likewise, configs needing external LUTs/ABSCO data
# or the unimplemented H2O override live in test/local/test_parameters/.

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

# Batch-processing API: BatchContext + update_model!
# Validates bit-equality of optical depths and radiances across scene updates,
# round-trip consistency, and guard rails.
@testset "BatchContext update_model!" begin include("test_update_model.jl") end

# v0.6 source-term refactor — AbstractSource vocabulary, SolarBeam, BlackbodySource,
# SurfaceSIF, surface_source_contribute! double-dispatch, prepared_sources flow.
# Includes the full Phase 1 → 5b regression assertions and end-to-end bit-equality
# checks via small CPU rt_run scenarios on PureRayleighParameters.
@testset "Sources (v0.6)" begin include("test_sources.jl") end

# v0.7 Phase A — ThermalEmission per-layer Planck volume source.
# Pins the corrected TIR weight rule (overrides REFACTOR_SPEC_v6 §2.11):
# isotropic source fires only at m=0 with 2π factor to undo downstream 0.5/π.
@testset "ThermalEmission (v0.7 Phase A)" begin include("test_thermal_emission.jl") end

# VLIDORT 2.8.3 baseline validation suite — Siewert 2000 PROBLEM_IIA Stokes-I,
# solar_tester scalar (Task 1), solar_tester vector (Task 1, IQU). Reference
# tables are committed under test/vlidort_baseline/reference_data/; no PyVLIDORT
# / Fortran runtime needed. The included file already wraps its own
# `@testset "VLIDORT baseline"`, so we don't double-nest.
include("vlidort_baseline/runtests.jl")

finally
    cd(_VSMARTMOM_TEST_ORIGINAL_CWD)
end
