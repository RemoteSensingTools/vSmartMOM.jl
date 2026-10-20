# vSmartMOM tests

Run the default CPU suite from this directory:

```bash
julia --project=. -e 'using Pkg; Pkg.develop(path=".."); Pkg.instantiate()'
CUDA_VISIBLE_DEVICES='' julia --project=. runtests.jl
```

Other runners: [`local/runtests.jl`](local/README.md) (GPU and external data),
`vlidort_baseline/runtests.jl` and `StandaloneSS/runtests.jl` (both also
included by `runtests.jl`).

## Test-file inventory

A file being present does not mean it runs. Each root `test_*.jl` file is
listed below with one of these statuses:

- **active**: included by `runtests.jl`, so it runs in the default CPU suite.
- **active + flagged**: active, but some cases run only when an environment
  flag or a functional CUDA device is present (noted per file).
- **flagged**: included by `runtests.jl`, but every case needs a flag or CUDA.
- **optional**: not in the default runner; it needs external data or another
  environment.
- **benchmark**: reports timings; it is not a correctness test.
- **support**: helper code for other tests, not a suite.
- **historical**: kept for reference; its thresholds are not current accuracy
  gates.

Flags: `GPU` = `VSMARTMOM_JACOBIAN_GPU_TEST=true` (and a CUDA device);
`CUDA` = runs automatically when a CUDA device is functional; `Metal` = runs
automatically on macOS with Metal.jl installed and `Metal.functional()`.
Backend selectors (`VSMARTMOM_SOURCE_GPU_TEST=true`, `EXTERNAL_SOLAR_TEST_GPU=1`,
`VSMARTMOM_COXMUNK_ARCH`) run the same cases on another backend instead of the
CPU; they do not add cases.

| File | Status | What it covers |
| --- | --- | --- |
| `test_Absorption.jl` | active + flagged | HITRAN parsing, line shapes and cross sections; CUDA arm when `Architectures.has_cuda()` |
| `test_Aerosols.jl` | active | Aerosols module: refractive-index database, TOMAS-15 and two-moment schemes |
| `test_CoreRT.jl` | active | 6SV1 and Natraj reference benchmarks (fixtures in `benchmarks/`) |
| `test_Scattering.jl` | active + flagged | Mie and Greek coefficients; slow Wigner 3-j cases only with `VSMARTMOM_FULL_TESTS=true` |
| `test_SolarModel.jl` | active | Planck spectra |
| `test_absorption_batched.jl` | optional | Batched vs per-layer absorption on a real ABSCO table; activates the package environment, needs an AtmosphericAbsorption API newer than the pinned 0.1.2, and exits 0 without testing when that API or the table is absent |
| `test_absorption_column_integration.jl` | active | Column-exact CIA and MT_CKD integration; external-file cases are `@test_skip` when `VSMARTMOM_TEST_O2_O2_CIA` / `VSMARTMOM_TEST_MTCKD` files are absent |
| `test_absorption_pressure.jl` | active + flagged | Pressure derivatives of line/LUT cross sections; CUDA arm under `CUDA` |
| `test_abundance_jacobian_contract.jl` | active | Gas-abundance Jacobian contract: columns whose forward absorption also responds to the abundance (line-by-line self broadening, MT_CKD, CIA) are rejected unless `abundance_jacobians=:frozen_cross_section`; exact ABSCO H₂O-broadener term (FD-checked) |
| `test_aerosol_reference.jl` | active | Normalized size-distribution tangents |
| `test_aqua.jl` | active | Aqua package health |
| `test_batched_kernels.jl` | active + flagged | Batched solve/inverse mutation contract; smoke arm under `Metal` |
| `test_canopy.jl` | active | Canopy surface construction and YAML |
| `test_component_m_max.jl` | active | Per-component Fourier `m_max` traits |
| `test_copy_parameters.jl` | active | Independent trial parameters and read-only LUT sharing |
| `test_coxmunk.jl` | active | Cox–Munk Fresnel and BRDF units |
| `test_coxmunk_polarization.jl` | active | Cox–Munk Stokes convention anchors (uses `coxmunk_jones_reference.jl`) |
| `test_coxmunk_rt_contract.jl` | active | Cox–Munk radiance and Jacobian contract (`VSMARTMOM_COXMUNK_ARCH` selects the backend) |
| `test_elemental_boundary_jacobians.jl` | active + flagged | Elemental derivatives at boundary inputs; CUDA arm under `GPU` |
| `test_elemental_limit_form.jl` | active + flagged | Elemental t⁺⁺/J₀⁺ factors near coincident ordinates: BigFloat accuracy (elastic, linearized, canopy, Raman kernels), merge-window continuity, Jacobian vs FD near μ₀ ≈ μᵢ; CUDA arm under `VSMARTMOM_JACOBIAN_GPU_TEST` |
| `test_external_solar_raman.jl` | active | External-solar RRS columns (`EXTERNAL_SOLAR_TEST_GPU=1` for GPU) |
| `test_external_solar_sfi.jl` | active | External-μ₀ SFI regression (`EXTERNAL_SOLAR_TEST_GPU=1` for GPU) |
| `test_float32.jl` | active | Float32/Float64 consistency |
| `test_forward_lin.jl` | historical | Legacy EMIT analytic-vs-FD Jacobians; 1–50 % thresholds with `isnan` allowances, so a pass is not accuracy evidence (see file header) |
| `test_forward_noRS.jl` | active + flagged | Elastic forward-model smoke test; CUDA arm under `CUDA` |
| `test_forward_raman.jl` | optional | Legacy RRS forward workflow; needs `test_parameters/O2Parameters.yaml` (CPU) or `O2Parameters_GPU.yaml` (GPU), which are not in `test_parameters/` (local copies live in `local/test_parameters/`) |
| `test_forward_raman_phase1b.jl` | flagged | RRS forward vs frozen reference (inelastic part rescaled in the test for the linear-interpolation line placement); runs only with `PHASE1B_CPU=1` or CUDA, otherwise one `@test_skip` |
| `test_forward_ss.jl` | active | Single-scatter driver smoke test |
| `test_fourier_convergence.jl` | active | Fourier early-exit strategies |
| `test_fused_doubling.jl` | active + flagged | Fused doubling kernels on KA-CPU; CUDA arm under `CUDA` when device `BENCH_DEV` (default 0) can be selected |
| `test_fused_interaction.jl` | active + flagged | Fused interaction kernels on KA-CPU; CUDA arm under `CUDA` when device `BENCH_DEV` (default 0) can be selected |
| `test_h2o_self_broadening.jl` | active | H2O self-broadening with moist mole fraction |
| `test_helpers.jl` | support | `VSmartMOMTestSupport`: shared linearized/forward run helpers and `rel_errors` |
| `test_hitran_cache.jl` | active | HITRAN download-cache identity and coverage (offline fixtures) |
| `test_hybrid_ad.jl` | active + flagged | CPU linearized FD check, ForwardDiff Mie; GPU-vs-CPU linearized and forward arms under `CUDA` |
| `test_inelastic_vs_smoke.jl` | active | Vibrational Raman call-chain smoke test |
| `test_io_exports.jl` | active | Top-level IO exports |
| `test_io_validation.jl` | active | IO format registry and validation |
| `test_jacobian_batched.jl` | active + flagged | Batched Jacobian propagation on CPU; CUDA arm under `GPU` |
| `test_jacobian_blocked.jl` | active + flagged | Blocked tangent products on CPU; CUDA arm under `GPU` |
| `test_jacobian_tiled.jl` | flagged | Shared-memory/tiled CUDA Jacobian products; included once by `runtests.jl` under `GPU` |
| `test_jacobian_upstream_contract.jl` | active | Selected derivatives require upstream availability |
| `test_jacobian_vendor.jl` | flagged | CUDA vendor Jacobian pointer plans; all cases under `GPU` |
| `test_jacobians_unit.jl` | active | Linearized-RT Jacobians vs finite differences, plus an ordinary-forward consistency check |
| `test_jet.jl` | active + flagged | JET snapshot; advisory by default, the baseline-count assertion runs only with `VSMARTMOM_JET_STRICT=1` |
| `test_julia_interfaces.jl` | active | RTModel property introspection, GreekCoefs tolerance keywords, elastic Cabannes overload |
| `test_lambertian_closure.jl` | active | Analytic Lambertian-albedo closure |
| `test_lambertian_jacobian_batched.jl` | active + flagged | Batched Lambertian Jacobians on CPU; CUDA arm under `GPU` |
| `test_lambertian_surfaces.jl` | active | Lambertian surface flavours |
| `test_lin_band_selection.jl` | active | Linearized RT solves one band per call: a one-element band vector matches the integer band, and several bands are rejected |
| `test_lin_forward_loop_parity.jl` | active | Forward/linearized Fourier loop-bound parity |
| `test_local_jacobian.jl` | active + flagged | Local optical-basis Jacobians on CPU (uses `local_jacobian_fixture.jl`); CUDA arm under `GPU` |
| `test_m_invariant_dtau.jl` | active | m-invariant dτ cache, bitwise A/B |
| `test_mie_nodes.jl` | active + flagged | Caller-defined-node Mie API; CUDA arm under `CUDA`, Metal arms under `Metal` |
| `test_mie_nodes_batched.jl` | active + flagged | Batched caller-node Mie API; CUDA arm under `CUDA` |
| `test_molecular_rayleigh.jl` | active | Molecular Rayleigh optical properties; N₂/O₂ NIR dispersion fits, dispersive tensor, H₂ refractivity check |
| `test_rayleigh_anisotropy.jl` | active | `rayleigh_anisotropy` parsing; `depol: -1` Rayleigh depolarization in the forward, linearized and batch builders |
| `test_multisensor_heights.jl` | active | Observer heights and interior levels |
| `test_nai2_fixed_radius_grid.jl` | active | Fixed radius grid of the log-normal Mie paths: nodes independent of μ/σ, size tangents vs central differences (optical properties and TOA radiance at the JacobianTestFast radius settings), `scattering.r_min`, PCW on the same nodes |
| `test_nonscattering_layers.jl` | active | Layers with zero scattering optical depth: equality with vanishing scattering and with the merged column at TOA and BOA, and with vanishing scattering at interior levels (IQU, albedo 0 and 0.3), the pure-absorber thermal closed form, Float32/Float64 thermal parity, and linearized parity |
| `test_parameters_parser.jl` | active | Parameter parser regressions |
| `test_performance.jl` | benchmark | Forward vs linearized timing report (no `@testset`) |
| `test_perturb_parameters.jl` | active | `perturb_parameters` utility |
| `test_phase0_inplace.jl` | active | In-place doubling/interaction updates, bitwise |
| `test_phase6_ports.jl` | active | Ported benchmark scripts parse (light units) |
| `test_phase_interpolation.jl` | active + flagged | Phase interpolation on CPU; CUDA arm under `GPU` |
| `test_quadpoints_streams.jl` | active | `QuadPoints.Nstreams` |
| `test_quality.jl` | active | Taplo schema configuration |
| `test_radau_overhead_sun.jl` | active | RadauQuad near an overhead sun: weights and geometry at SZA 0 in Float32/Float64, and continuity from SZA 0 to 1e-6, 1e-5, 1e-4 (sub-interval nodes ulps apart) and Gauss |
| `test_ray_cab_raman.jl` | active | Rayleigh/Cabannes/RRS identities |
| `test_raman_placzek_teller.jl` | active | Placzek–Teller sum rule of the per-J Q/S/O lines; per-J Cabannes depolarization |
| `test_raman_rrs_lines.jl` | active | RRS line weights on the band grid: a band without rotational lines is rejected; linear-interpolation line placement; offsets from the reference grid point; removed `getRamanSSProp!` overloads |
| `test_raman_vs_plus.jl` | active | `VS_0to1_plus` forward run over concatenated bands (pure Rayleigh, one layer); source-band depletion of both VS `_plus` modes; VS `_plus` builder against the standard builder for both `rayleigh_anisotropy` settings, `params` untouched, Float32 and line-by-line absorption |
| `test_rt_entry_guards.jl` | active | Solve-boundary guards: reused-RS F₀, `sources=` applies to its call only, HDR/BHR direct beam, TMS streams/bands, multisensor Cox–Munk TOA, U/V beams, thermal paths, linearized surfaces |
| `test_scenario_sweep.jl` | active | Scenario sweeps vs monolithic `rt_run` |
| `test_schema_docs.jl` | active | Schema documentation invariants |
| `test_selective_jacobians.jl` | active | Retrieval-selective Jacobians |
| `test_sif.jl` | active | SIF injection and loaders (loader cases skip without their fixtures) |
| `test_source_adding.jl` | active | Equivalent-source adding (`VSMARTMOM_SOURCE_GPU_TEST=true` runs it on GPU) |
| `test_source_adding_sif.jl` | active | Equivalent-source adding with SIF (`VSMARTMOM_SOURCE_GPU_TEST=true` runs it on GPU) |
| `test_source_validation.jl` | active | Source request validation |
| `test_sources.jl` | active | Source-term vocabulary and routing |
| `test_sparse_combine.jl` | active | Sparse vs dense aerosol combine, bitwise |
| `test_surface_split.jl` | active | Atmosphere/surface split vs monolithic `rt_run` |
| `test_thermal_emission.jl` | active | Thermal emission volume source |
| `test_tms.jl` | active | TMS single-scattering correction |
| `test_truncation.jl` | active | Phase-function truncation invariants |
| `test_type_stability.jl` | active | Type stability |
| `test_update_model.jl` | active | BatchContext `update_model!` |
| `test_update_model_transaction.jl` | active | BatchContext failure-state integrity |
| `test_z_jacobian_tables.jl` | active + flagged | Tabulated phase Jacobians on CPU; CUDA arm under `GPU` |
| `test_z_tables.jl` | active | Tabulated Z moments, bitwise |

Independent references in the default suite: 6SV1 and Natraj (`test_CoreRT.jl`,
fixtures in `benchmarks/`), and Siewert 2000 plus VLIDORT 2.8.3 solar_tester
(`vlidort_baseline/`, tables in `vlidort_baseline/reference_data/`). Other root
scripts (`benchmark_rt.jl`, `benchmark_fwd_vs_lin.jl`) and most of
`benchmarks/` are benchmarks, not tests.
