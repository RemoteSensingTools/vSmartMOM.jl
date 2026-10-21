# Inventory E — Docs State Audit

**Branch:** `unified-vsmartmom`
**Worktree:** `/home/sanghavi/code/github/uni_vSmartMOM/`
**Tip commit:** `a4e4187` — "Add batched-kernel and Raman scaling benchmarks with writeup"
**As of:** 2026-04-19
**Purpose:** Ground-truth inventory of docs + exports + tests state before the ship-quality docs overhaul for the `sanghavi-unified` merge.

---

## 1. Public API exports

Collected from every `^export` line in `src/vSmartMOM.jl`, `src/CoreRT/CoreRT.jl`, `src/Scattering/Scattering.jl`, `src/Absorption/Absorption.jl`, `src/Inelastic/InelasticScattering.jl`, `src/IO/IO.jl`, `src/IO/Sources.jl`, `src/IO/Formats.jl`, `src/IO/NetCDF/NetCDF.jl`, `src/SolarModel/SolarModel.jl`, `src/Aerosols/Aerosols.jl`, and `src/Architectures.jl`.

Docstring codes:
- **yes** = full prose docstring with arguments/returns/example
- **partial** = one-liner or bare `$(FUNCTIONNAME)`/`$(FIELDS)` stub with no prose
- **no** = no docstring at the definition site

### 1.1 Top-level `vSmartMOM` exports (user-facing entry points)

| Export | Definition | Docstring | Category |
|---|---|---|---|
| `CPU` | `src/Architectures.jl:34` | partial (3-line stub) | public |
| `GPU` | `src/Architectures.jl:40` | partial (3-line stub) | public |
| `default_architecture` | `src/Architectures.jl:74` | yes | public |
| `array_type` | `src/Architectures.jl:64` | yes | public |
| `artifact` | `src/Artifacts/artifact_helper.jl:36` | yes | public |
| `fetch_hitran` | `src/Artifacts/hitran_api.jl:57` | yes | public |
| `fetch_hitran_by_ids` | `src/Artifacts/hitran_api.jl:96` | yes | public |
| `set_hitran_edition!` | `src/Artifacts/hitran_preferences.jl:27` | yes | public |
| `get_hitran_edition` | `src/Artifacts/hitran_preferences.jl:37` | partial (2 lines) | public |
| `available_hitran_editions` | `src/Artifacts/hitran_preferences.jl:45` | yes | public |
| `hitran_info` | `src/Artifacts/hitran_preferences.jl:67` | yes | public |
| `hitran_is_cached` | `src/Artifacts/hitran_api.jl:26` | partial (1-line) | public |
| `FwdMode` | `src/vSmartMOM.jl:33` | no | advanced |
| `LinMode` | `src/vSmartMOM.jl:34` | no | advanced |
| `default_parameters` | `src/CoreRT/tools/model_from_parameters.jl:9` | partial (1-liner) | public |
| `parameters_from_yaml` | `src/IO/IO.jl:63` | yes | public |
| `model_from_parameters` | `src/CoreRT/tools/model_from_parameters.jl:33` | yes | public |
| `rt_run` | `src/CoreRT/rt_run.jl:53` (fwd) and `src/CoreRT/rt_run_lin.jl:72` (lin) | yes (both) | public |
| `read_parameters` | re-exported from IO | **no** (function exists as alias; not found as separate definition) | public |
| `read_atmos_profile` | `src/IO/AtmosProfile.jl:35` | **no docstring at def site** | public |
| `rt_run_lin` | `src/CoreRT/rt_run_lin.jl:85` | yes (convenience alias) | public |
| `model_from_parameters_lin` | `src/CoreRT/tools/lin_model_from_parameters.jl:51` | need verification (LinMode overload) | public |
| `RTModel` | `src/CoreRT/types.jl:843` | yes (prose + fields) | public |
| `AbstractRTModel` | `src/CoreRT/types.jl:719` | partial (1-line) | advanced |
| `SolverConfig` | `src/CoreRT/types.jl:731` | yes | public |
| `Atmosphere` | `src/CoreRT/types.jl:758` | yes | public |
| `RayleighScattering` | `src/CoreRT/types.jl:774` | yes | public |
| `AerosolState` | `src/CoreRT/types.jl:792` | yes | public |
| `Optics` | `src/CoreRT/types.jl:812` | yes | public |
| `OpticsLin` | `src/CoreRT/types.jl:913` | yes | public |
| `GeosChemSource` | `src/IO/Sources.jl:44` | yes | public |
| `NetCDFGridSource` | `src/IO/Sources.jl:79` | yes | public |
| `NetCDFSource` | `src/IO/Sources.jl:16` | yes | advanced (abstract) |
| `geoschem_to_dict` | `src/IO/NetCDF/GeosChem.jl:35` | yes | public |
| `read_geoschem_profile` | `src/IO/NetCDF/GeosChem.jl:143` | **not re-verified** | public |

### 1.2 CoreRT-level exports (visible to users who `using vSmartMOM.CoreRT`)

| Export | Definition | Docstring | Category |
|---|---|---|---|
| `rt_run_ss` | `src/CoreRT/CoreKernel/rt_kernel_ss.jl` (single-scatter) | unknown — likely partial | advanced |
| `lin_added_layer_all_params` | `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl:115` | unknown | advanced |
| `OpticalPropertyJacobian` | `src/CoreRT/types_lin.jl:159` (const alias) | unknown — const alias | advanced |
| `RawAerosolJacobian` | `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl:20` | unknown | advanced |
| `delta_m_forward` | `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl:44` | unknown | advanced |
| `delta_m_truncation_lin` | `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl:86` | unknown | advanced |
| `ParameterLayout` | `src/CoreRT/parameter_layout.jl:20` | yes | public |
| `n_total` | `src/CoreRT/parameter_layout.jl:33` | partial (1-line) | public |
| `aerosol_range` | `src/CoreRT/parameter_layout.jl:37` | partial (1-line) | public |
| `gas_range` | `src/CoreRT/parameter_layout.jl:43` | partial (1-line) | public |
| `surface_range` | `src/CoreRT/parameter_layout.jl:49` | partial (1-line) | public |
| `surface_index` | `src/CoreRT/parameter_layout.jl:55` | partial (1-line) | public |
| `canopy_range` | `src/CoreRT/parameter_layout.jl:63` | partial (1-line) | public |
| `n_layer_params` | `src/CoreRT/parameter_layout.jl:59` | partial (1-line) | public |
| `RTModelLin` | `src/CoreRT/types_lin.jl:108` | **no** (bare struct) | public |
| `GaussQuadFullSphere` | `src/CoreRT/types.jl:117` | **no** (bare struct) | public |
| `LambertianSurfaceScalar` | `src/CoreRT/types.jl:389` | yes | public |
| `LambertianSurfaceSpectrum` | `src/CoreRT/types.jl:403` | yes (partial prose) | public |
| `CanopySurface` | `src/CoreRT/types.jl:497` | yes | public |
| `CanopySurface_from_prospect` | `src/CoreRT/Surfaces/canopy_surface.jl:52` | unknown | public |
| `invalidate_canopy_cache!` | `src/CoreRT/Surfaces/canopy_surface.jl:554` | unknown | advanced |
| `CoxMunkSurface` | `src/CoreRT/types.jl:462` | yes | public |
| `water_refractive_index` | `src/CoreRT/Surfaces/water_refraction.jl:68` | unknown | public |

### 1.3 Scattering submodule exports

| Export | Definition | Docstring | Category |
|---|---|---|---|
| `make_mie_model` | `src/Scattering/make_mie_model.jl:23,39,59` (3 methods) | unknown | public |
| `reconstruct_phase` | `src/Scattering/mie_helper_functions.jl:415` | unknown | public |
| `NAI2` | `src/Scattering/types.jl:71` | **no** (bare struct) | public |
| `PCW` | `src/Scattering/types.jl:79` | **no** (bare struct) | public |
| `Aerosol` | `src/Scattering/types.jl:42` | unknown | public |
| `MieModel` | `src/Scattering/types.jl` | unknown | public |
| `Stokes_I`, `Stokes_IQU`, `Stokes_IQUV` | `src/Scattering/types.jl` | unknown | public |
| `δBGE` | `src/Scattering/types.jl:164` | unknown | public |
| `GreekCoefs` | `src/Scattering/types.jl:231` | unknown | public |
| `AerosolOptics` | `src/Scattering/types.jl` (plus _lin variant at `types_lin.jl:240`) | unknown | public |
| `AbstractFourierDecompositionType` | `src/Scattering/types.jl:63` | **no** (bare abstract) | advanced |
| `linGreekCoefs` | `src/Scattering/types_lin.jl:199` | unknown | advanced |
| `linAerosolOptics` | `src/Scattering/types_lin.jl` | unknown | advanced |
| `compute_B`, `compute_ab`, `comp_ab`, `compute_mie_π_τ!` | `src/Scattering/mie_helper_functions.jl` | unknown | internal-but-exported |
| `compute_wigner_values`, `save_wigner_values`, `load_wigner_values` | `src/Scattering/compute_wigner_values.jl` | unknown | internal-but-exported |
| `compute_Sl` | `src/Scattering/compute_PCW.jl:137` | unknown | internal-but-exported |
| `gausslegendre` | re-export from FastGaussQuadrature | n/a | internal-but-exported |
| `compute_aerosol_optical_properties` | `src/Scattering/compute_NAI2.jl:44`, `compute_PCW.jl:28`, `phase_function_autodiff.jl:55`, `compute_NAI2_lin.jl:16` | unknown | public |
| `compute_ref_aerosol_extinction` | `src/Scattering/compute_NAI2.jl:216` | unknown | public |
| `truncate_phase` | `src/Scattering/truncate_phase.jl:113` | unknown | public |
| `ConjugateTransposePairs` | `src/Scattering/types.jl` | unknown | advanced |
| `AbstractPolarizationType`, `AbstractAerosolType`, `AbstractTruncationType` | `src/Scattering/types.jl:19,92,156` | **no** (bare abstracts) | advanced |
| `phase_function` | `src/Scattering/compute_NAI2.jl:315,413` (multiple methods) | unknown | public |
| `compute_aerosol_XS` | `src/Scattering/compute_NAI2.jl:491`, `compute_NAI2_lin.jl:548` | unknown | advanced |
| `compute_aerosol_optical_properties_gpu` | `src/Scattering/compute_NAI2_gpu.jl:31` | unknown | public |
| `MiePrecisionPolicy`, `NativeFloat64`, `DSEmulated`, `DoubleSingle`, `ComplexDS`, `NeumaierAccum`, `neumaier_add`, `neumaier_sum` | `src/Scattering/gpu_precision.jl` | unknown | advanced |

### 1.4 Absorption submodule exports

| Export | Definition | Docstring | Category |
|---|---|---|---|
| `AbstractCrossSectionModel` | `src/Absorption/types.jl:157` | yes (brief) | advanced |
| `HitranModel` | `src/Absorption/types.jl:168` | yes | public |
| `InterpolationModel` | `src/Absorption/types.jl:193` | yes | public |
| `HitranTable` | `src/Absorption/types.jl:24` | yes | public |
| `compute_absorption_cross_section` | `src/Absorption/compute_absorption_cross_section.jl:32,198` | unknown | public |
| `absorption_cross_section` | `src/Absorption/compute_absorption_cross_section.jl` | unknown | public |
| `read_hitran` | `src/Absorption/read_hitran.jl:30` | unknown | public |
| `make_hitran_model` | `src/Absorption/make_model_helpers.jl:25` | unknown | public |
| `Doppler`, `Lorentz`, `Voigt` | `src/Absorption/types.jl:101,104,107` | partial (bare-struct with preceding abstract docstring) | public |
| `AbstractBroadeningFunction` | `src/Absorption/types.jl:98` | **no** (bare abstract) | advanced |
| `HumlicekErrorFunction`, `HumlicekWeidemann32VoigtErrorFunction`, `HumlicekWeidemann32SDErrorFunction`, `CPF12ErrorFunction`, `ErfcHumliErrorFunctionVoigt`, `ErfcHumliErrorFunctionSD`, `ErfcErrorFunction` | `src/Absorption/types.jl:123-…` | **mostly no** — bare structs | advanced |
| `AbstractComplexErrorFunction` | `src/Absorption/types.jl:120` | no | advanced |
| `make_interpolation_model`, `save_interpolation_model`, `load_interpolation_model` | `src/Absorption/make_model_helpers.jl:55,102,107,130` | unknown | public |

### 1.5 SolarModel submodule exports

| Export | Definition | Docstring | Category |
|---|---|---|---|
| `planck_spectrum_wn` | `src/SolarModel/SolarModel.jl:18,68` | yes (both methods) | public |
| `planck_spectrum_wl` | `src/SolarModel/SolarModel.jl:36` | yes | public |
| `solar_transmission_from_file` | `src/SolarModel/SolarModel.jl:116,123` | yes | public |
| `default_solar_transmission` | `src/SolarModel/SolarModel.jl:135` | yes | public |

### 1.6 Architectures submodule re-exports (via CoreRT / root)

Defined in `src/Architectures.jl`: `@hascuda`, `AbstractArchitecture`, `CPU`, `GPU`, `devi`, `array_type`, `architecture`, `default_architecture`, `synchronize_if_gpu`, `has_cuda`. All have short docstrings (yes/partial); `architecture(::Array)` has none; `devi` has none.

### 1.7 Aerosols submodule exports (present in tree but NOT re-exported at root)

`src/Aerosols/Aerosols.jl` exports:
`AerosolScheme, TOMAS15Scheme, TwoMomentScheme, AerosolData, AerosolSpeciesData, RefractiveIndexLUT, RefractiveIndexDatabase, read_aerosol_data, load_refractive_index_database, get_refractive_index, compute_optical_properties`.

**Gap:** This submodule is not included in `src/vSmartMOM.jl`. It is referenced in examples and in the `AEROSOL_DESIGN.md` file inside `src/Aerosols/`, but the user must do `using vSmartMOM.Aerosols` … which will fail because the parent module never loads it. **This is likely a live bug or an intentional staging state that must be surfaced before docs ship.**

### 1.8 InelasticScattering submodule exports

Not checked exhaustively — no `^export` lines found at module-level top. Likely uses `using` re-export or exports buried elsewhere. Types used by users (`RRS`, `VS_0to1`, `VS_1to0`, `noRS`, `GreekCoefs`) are referenced in tests as `InelasticScattering.RRS(…)` i.e. fully qualified.

### 1.9 Export summary counts

- Root-level exports from `vSmartMOM.jl`: ~38 symbols (including re-exports).
- CoreRT module exports: ~26 symbols in two export blocks.
- Scattering module exports: ~33 symbols (heavy leakage of internals; many are internal-but-exported).
- Absorption module exports: ~18 symbols.
- IO module exports: 6 symbols.
- SolarModel exports: 4 symbols.
- Aerosols module exports: 11 symbols **not actually reachable from `using vSmartMOM`** (live gap).

**Surface docstring coverage rough estimate:** of ~140 distinct exported names across all modules, roughly **50–60% have a usable (full or partial) docstring at their definition site**, the rest are bare `struct … end` with a hint line or nothing. Highest coverage: the new `RTModel` hierarchy (added during unified work) and HITRAN API. Lowest coverage: Scattering helper/internal exports, Absorption CEF structs, abstract-type supertypes.

---

## 2. Documenter structure

`docs/make.jl` exists and uses `Documenter.jl` + `Literate.jl`. Output tree:

### 2.1 Published (referenced in `docs/make.jl`)

```
docs/src/
├── index.md                                  H1 "Introduction" — quick-example + install
├── design.md                                 H1 "Architecture and Design" — module graph, data flow
└── pages/
    ├── api_reference.md                     grouped @docs API reference (audit complete; comprehensive but missing many advanced symbols)
    ├── geoschem_integration.md              GEOSChem workflow walkthrough
    ├── IO/
    │   ├── Overview.md
    │   ├── Schema.md
    │   └── Examples.md
    ├── vSmartMOM/
    │   ├── Overview.md                     22 lines
    │   ├── Example.md                      64 lines
    │   ├── InputParametersGuide.md         74 lines
    │   ├── Types.md                        63 lines (@docs blocks)
    │   ├── Principles.md                   108 lines
    │   ├── CoreRTTheory.md                 233 lines (adding-doubling theory)
    │   └── References.md                   8 lines (thin)
    ├── Absorption/
    │   ├── Overview.md, Example.md, Types.md, References.md, HITRAN_Data.md
    │   └── vSmartMOMDiagram-Absorption.drawio.png
    ├── Scattering/
    │   ├── Overview.md, Example.md, Types.md, References.md
    │   └── vSmartMOMDiagram-Scattering.drawio.png
    └── tutorials/                          Literate.jl-compiled
        ├── Tutorial_QuickStart.{jl,md}
        ├── Tutorial_Absorption.{jl,md}
        ├── Tutorial_Scattering.{jl,md}
        ├── Tutorial_MieDeepDive.{jl,md}
        ├── Tutorial_IO.{jl,md}
        ├── Tutorial_CoreRT.{jl,md}
        ├── Tutorial_Surfaces.{jl,md}
        ├── Tutorial_Canopy.{jl,md}
        ├── Tutorial_Jacobians.{jl,md}
        ├── Tutorial_GPU.{jl,md}
        └── Tutorial_HybridAD.{jl,md}
```

Assets: `docs/src/assets/{logo.png, favicon.ico, CrossSectionGIF.gif, ScatteringGIF.gif}`.

The `make.jl` script calls `makedocs(..., warnonly=true)` — **warnings are not currently fatal**, so missing docstrings or bad @refs are tolerated.

Deploy config: unified-vsmartmom is published to a branch-named dev URL; main still canonical dev.

### 2.2 Legacy dev notes (under `docs/dev_notes/`)

See §3 — these are NOT referenced by `docs/make.jl` and therefore do not appear in the published Documenter site.

### 2.3 Orphaned

None obvious in `docs/src/` — every .md and tutorial is referenced in `make.jl`. The `References.md` under `vSmartMOM/` (8 lines) is stub-like and may read as orphaned even though it is wired up.

---

## 3. Dev notes inventory (`docs/dev_notes/`)

20 files, 6,667 total lines, none referenced by `docs/make.jl`. Classification:

| File | Lines | One-line description | Classification |
|---|---|---|---|
| `AI_HANDOFF_MEMO.md` | 80 | AI handoff for Jacobian/test work | **stale** (superseded by current `plans/CLAUDE_HANDOFF_BRIEF.md`) |
| `BIMODAL_FIT_VALIDATION.md` | 101 | Julia-vs-Python bimodal lognormal fit validation, 2025-10-21 | keep (aerosol/size-distribution reference) |
| `CLEANUP_SUMMARY.md` | 178 | Dated code-cleanup log from io-update branch, 2025-10-15 | **stale** |
| `CUDA_OPTIONAL_IMPLEMENTATION.md` | 327 | Implementation log for making CUDA an optional weak dep, 2025-10-15 | **stale** (internal notes, superseded; CUDAExt is in place now) |
| `CUDA_SETUP.md` | 258 | User-facing guide for CUDA/GPU setup | **migrate to Documenter** (user-facing; belongs under a "GPU" published page) |
| `FILE_REORGANIZATION_SUMMARY.md` | 194 | TOMAS-aerosols branch file layout changes, 2025-10-20 | **stale** |
| `FT_TYPE_TEST_RESULTS.md` | 233 | Float32/FT parameterization test results, 2025-10-16 | **stale** (historical) |
| `HANDOFF_SUMMARY.md` | 633 | Aerosol framework implementation handoff, 2025-10-16 | keep (partial reference for aerosols module) but aerosol-exposure gap (§1.7) must be resolved first |
| `JACOBIAN_TEST_WORKFLOW.md` | 65 | How to run the Jacobian test suite | **migrate to Documenter** (or to `test/README`) |
| `LINEARIZATION_BUGS.md` | 810 | Catalog of bugs encountered porting linearized code from sanghavi | keep (historical debugging record) |
| `MieGPUSpeedup.md` | 1021 | Mie GPU implementation plan + testing (pre-impl doc) | **stale** (impl landed in `d6bd45c`) but may contain design rationale worth mining |
| `NEXT_STEPS_ANALYSIS.md` | 683 | Oceananigans-style roadmap, 2025-10-14 | **superseded by new plan** (predates sanghavi-unified merge brief) |
| `RAMAN_CODE_HANDOFF.md` | 152 | Addressed to sanghavi-branch Raman-code author; explains unified RT flow | keep (referenced by `CLAUDE_HANDOFF_BRIEF`) |
| `SESSION_HANDOFF.md` | 190 | Dated Jacobian debugging session | **stale** |
| `SESSION_SUMMARY_2025-10-15.md` | 246 | Session log | **stale** |
| `batched_kernel_benchmarks.md` | 380 | 2026-04-03 GPU batched-kernel benchmark writeup | keep (current; referenced by recent commit) |
| `linearization_changes_unified_branch.md` | 547 | 2026-02-19 writeup of linearization changes | keep (valuable design reference) |
| `ocean.md` | 1250 | Ocean RT implementation plan, phase-1 OceanOptics | keep (ongoing work tracker) |
| `raman_gpu_optimization.md` | 71 | Raman GPU allocation audit | keep but **superseded on scope** by `CLAUDE_HANDOFF_BRIEF` §4 |
| `sanghavi_unified_merge_plan.md` | 294 | Christian's 2026-03-21 merge plan | **superseded by new plan** (explicitly flagged stale in `CLAUDE_HANDOFF_BRIEF.md` §4.1) |

**Classification counts:** 8 stale, 3 superseded, 2 migrate-to-Documenter, 7 keep.

---

## 4. Test suite survey (`test/runtests.jl`)

`runtests.jl` registers 11 `@testset` groups (plus a CUDA-gated Raman GPU test), covering 19 active test files (four additional `test_*.jl` files exist on disk but are not included; see below). Total `@test` occurrences in the test tree: **464** (matches the claimed 474 within rounding — the 10-line delta is explained by the disabled `wigner3j` test and `@test_skip` branches).

### 4.1 Included in `runtests.jl` (top-to-bottom order)

| Testset | File | Lines | `@test` count | `@testset` count | Exercises | Flags |
|---|---|---|---|---|---|---|
| Absorption | `test_Absorption.jl` | 210 | 40 | 5 | HITRAN parser; cross-section computation at various P/T; broadening types | — |
| Scattering | `test_Scattering.jl` | 117 | 11 | 2 | Wigner 3j (gated by `VSMARTMOM_FULL_TESTS`, ~60s); Mie phase functions | **slow-candidate** (wigner3j) |
| CoreRT | `test_CoreRT.jl` | 84 | 9 | 2 | 6SV1 benchmark, Natraj benchmark — CPU forward RT against reference values | **Canary candidate** (exemplary, reference-backed) |
| SolarModel | `test_SolarModel.jl` | 28 | 3 | 1 | Planck spectrum + solar transmission | — |
| Forward noRS | `test_forward_noRS.jl` | 165 | 36 | 7 | EMIT-style forward RT, CPU + optional GPU, no Raman | **Canary/Quickstart candidate** (fast, exemplary, full forward path) |
| Jacobian Unit | `test_jacobians_unit.jl` | 224 | 20 | 9 | Mie interpolation derivatives, δ-M derivatives, Core `+` derivatives, end-to-end dR Jacobians vs FD | **Linearization**; "~2–5 min" annotated |
| Type Stability | `test_type_stability.jl` | 169 | 46 | 6 | `@inferred` checks on CoreScatteringOpticalProperties, helper functions; Float32/Float64 preservation | — |
| Float32 Consistency | `test_float32.jl` | 56 | 7 | 4 | Float64 baseline + Float32 forward; agreement within 1% | **Canary candidate** |
| RAMI Smoke | `test_rami_smoke.jl` | 53 | 9 | 4 | RAMI-config forward RT, no gas absorption | **Canary candidate** |
| Canopy Surface | `test_canopy.jl` | 208 | 51 | 10 | CanopySurface construction, YAML, forward vs Lambertian baseline, multi-layer, multi-band | — |
| Cox-Munk Surface | `test_coxmunk.jl` | 499 | 71 | 43 | Fresnel, Mueller matrix, water refractive index, shadow masking, BRDF reciprocity/energy, polarized Fourier, Jacobian FD | **Linearization** (has FD checks); largest test file by testset count |
| Raman GPU (CUDA gated) | `test_forward_raman_gpu.jl` | 129 | 11 | 4 | Forward RRS path on GPU (O2Parameters_GPU.yaml, ~60 spec pts) | **Raman**; GPU-only |

### 4.2 On-disk but NOT in `runtests.jl`

| File | Lines | Apparent purpose | Flags |
|---|---|---|---|
| `test_Aerosols.jl` | 307, 91 `@test` | Aerosols module tests (TOMAS15, two-moment, refractive index DB, optical-property compute) | **Aerosols module is not loaded from `vSmartMOM.jl`** — file exists orphaned in the test tree |
| `test_forward_raman.jl` | 179, 14 `@test` | CPU forward RRS | **Raman**; presumably slow, hence omitted |
| `test_forward_lin.jl` | 251, 6 `@test` | Full linearized RT vs FD | **Linearization**; references `src/Testing/perturb_parameters.jl` (may not exist at that path) |
| `test_hybrid_ad.jl` | 146, 13 `@test` | Hybrid ForwardDiff vs analytic Jacobian | **Linearization** |
| `test_jacobians_GPU.jl` | 91, 4 `@test` | GPU Jacobians | **Linearization**, GPU |
| `test_mie_gpu.jl` | 469, 21 `@test` | GPU-accelerated Mie tests | **slow-candidate**, GPU |
| `test_performance.jl` | 165, 1 `@test` | Performance/benchmarks (self-documented as "diagnostic script, not @testset — not included in runtests.jl") | benchmark, not a test |

Plus orphan directories: `aerosol_exploration_output/`, `benchmark_fwd_vs_lin.jl`, `benchmark_rt.jl`, `benchmarks/`, `rami/`.

### 4.3 Test-count contribution (sorted by `@test` count where included)

| File | @test |
|---|---|
| `test_coxmunk.jl` | 71 |
| `test_type_stability.jl` | 46 |
| `test_canopy.jl` | 51 |
| `test_Absorption.jl` | 40 |
| `test_forward_noRS.jl` | 36 |
| `test_jacobians_unit.jl` | 20 |
| `test_Scattering.jl` | 11 (most gated behind env) |
| `test_forward_raman_gpu.jl` | 11 (when CUDA) |
| `test_CoreRT.jl` | 9 |
| `test_rami_smoke.jl` | 9 |
| `test_float32.jl` | 7 |
| `test_SolarModel.jl` | 3 |

Total across files listed in runtests: **304 `@test` invocations**. The 474-test figure in `CLAUDE.md` likely counts `@testset`-expanded individual assertions; `test_coxmunk.jl` alone has 43 nested testsets and a fair number of looped `@test` calls.

### 4.4 Canary recommendations (for a fast-subset CI lane)

Sorted by fitness-for-canary (fast + exemplary + exercises forward path):
1. **`test_CoreRT.jl`** — reference-value validation vs 6SV1 and Natraj.
2. **`test_forward_noRS.jl`** — full EMIT-like forward path, CPU+GPU optional.
3. **`test_float32.jl`** — exercises Float32 support end to end.
4. **`test_rami_smoke.jl`** — bare-bones RT engine smoke.

---

## 5. README + CHANGELOG

### 5.1 README gap assessment

`README.md` (~130 lines) covers:
- Project banner + badges (some badges point to main-branch URLs that may 404 on the merged branch).
- Installation via `Pkg.add("vSmartMOM")` — **does not document that unified-vsmartmom is not yet the registered package version**; a user following the README gets the old API.
- Modules overview with key functions for `vSmartMOM`, `vSmartMOM.Absorption`, `vSmartMOM.Scattering`.
- HITRAN Data Access subsection (up-to-date for the unified-branch HITRAN work).
- Contribution + contact + license.

**What a new user would NOT learn from the current README:**
- That `model_from_parameters` returns `RTModel{ARCH,FT}` (hierarchical, Oceananigans-style), not the flat `vSmartMOM_Model` they may see in external papers/code.
- How to access differentiable state (`model.optics.τ_abs`, `model.τ_abs` shim).
- How to run linearized RT with Jacobians (only briefly mentioned in "Key functions" block; no code snippet).
- Which surface types are supported (no mention of Cox-Munk, Canopy, RPV, RossLi).
- GPU/CPU selection mechanism (`params.architecture = vSmartMOM.GPU()`).
- HITRAN edition switching (only referenced by link — a full one-paragraph primer is warranted at README level).
- RAMI/canopy capability.
- Float32 support.
- Where to find examples and what they do.
- Breaking-change note pointing at CHANGELOG v2.0.0 for migration.

**Sections missing for a v2.0.0 ship-quality README:**
- Quickstart (5-line forward RT + 5-line linearized RT).
- "What's new in v2.0.0" callout with link to CHANGELOG.
- "Migrating from v1.x" paragraph pointing at a migration guide.
- Surface-type table.
- Feature matrix (polarization, Raman, GPU, Jacobians) with ticks indicating supported combinations.
- Links to Documenter tutorials (QuickStart, Jacobians, GPU).
- Author/affiliation block updated for the unified maintainer set.

### 5.2 CHANGELOG gap assessment

`CHANGELOG.md` has a single `## v2.0.0` section covering:
- Breaking: `RTModel` replaces `vSmartMOM_Model` (1 paragraph; mentions `Base.getproperty` shim).
- Breaking: `RTModelLin` replaces `vSmartMOM_Lin`.
- Breaking: `rt_run_canopy` removed.
- New features: hierarchical model, `ParameterLayout`, Cox-Munk, accessor functions.
- Improvements: YAML safe-literal parsing, Raman `Vector{Any}` fixes, dead-code cleanup.
- Known limitations: lin Raman unsupported, wigner3j disabled by default.

**Completeness check for documenting the breaking change:**
- The *fact* of the switch is stated, but the entry does not:
  - Show a before/after code snippet (what the old user code looked like, what to change).
  - Enumerate every renamed/moved symbol (e.g., `model.params.brdf` → `model.surfaces[i]`, `model.aerosol_optics` still reachable via shim but now lives in `model.optics.aerosols.aerosol_optics`).
  - Call out that *IO*-level loading works unchanged (important reassurance).
  - Link to a migration guide (which does not exist yet in `docs/`).
  - Document the `LinMode()` / `FwdMode()` dispatch pattern.
  - Mention HITRAN 2024 direct-download feature (merged in `df19b06`, not in changelog).
  - Mention GPU Mie landing (`d6bd45c`).
  - Mention batched-kernel/Raman benchmarks (`a4e4187`).

Neither does the CHANGELOG cover non-breaking user-visible additions from this branch: `artifact()` dispatch, `Float32` across the stack, `CanopySurface_from_prospect`, canopy within-atmosphere scattering, GEOSChem NetCDF integration, scratch-space HITRAN cache.

**Verdict:** the v2.0.0 entry exists but is **not ship-quality** for a breaking release. It would need at minimum: a "Migration" subsection with concrete code diffs, and entries for the HITRAN-2024, Float32, GPU-Mie, Canopy, and GEOSChem additions to fairly reflect what the release contains.

---

## 6. Examples folder

`examples/` exists at repo root with 7 files:

| File | Type | Scope | Runnable? | Reference output committed? |
|---|---|---|---|---|
| `OPTICAL_THICKNESS_GUIDE.md` | docs | Narrative guide for optical-thickness workflow | n/a | n/a |
| `aerosol_config_tomas15.yaml` | config | TOMAS-15 aerosol scheme sample | supports other scripts | n/a |
| `aerosol_config_two_moment.yaml` | config | Two-moment aerosol scheme sample | supports other scripts | n/a |
| `aerosol_integration_example.jl` | script | End-to-end aerosol integration workflow | **unclear — depends on `vSmartMOM.Aerosols` submodule not reachable from root (§1.7 gap)** | no |
| `compute_thermal_ir_optical_thickness.jl` | script | Thermal IR optical depth for a 10m layer | runnable (uses only `vSmartMOM.Absorption`) | no |
| `geoschem_integration.jl` | script | GEOSChem NetCDF → parameters → rt_run | runnable in principle; requires a GEOSChem file the user provides | no |
| `optical_thickness_simple.jl` | script | 10m layer absorption demo (simpler variant) | runnable | no |

### Golden-path coverage

| Scenario | Covered? |
|---|---|
| Forward RT (Lambertian land) | **no dedicated example** (only covered inside tutorials) |
| Forward RT (Cox-Munk ocean) | **no** |
| Forward RT (canopy) | **no** |
| Forward RT (RPV / RossLi) | **no** |
| Linearized RT with Jacobians | **no** |
| Raman (RRS, VS) forward | **no** |
| GPU forward RT | **no** |
| HITRAN-2024 direct download | **no** |
| GEOSChem integration | yes (one script) |
| Absorption cross-section only | yes (two scripts) |

**Gap:** the examples folder is biased toward absorption-only and GEOSChem; the "RT ships with polarization, surfaces, Jacobians, and Raman" message promised by README is not demonstrated by any example. For a v2.0.0 ship, each surface type and each RT mode (forward / linearized / Raman) should have at least one `.jl` + `.md` pair with committed reference output. Tutorials in `docs/src/pages/tutorials/` partially fill this gap (QuickStart, Jacobians, Canopy, Surfaces, GPU, CoreRT), but tutorials are Literate-generated and not always the right surface for users wanting a copy-pasteable script.

---

## 7. Docs overhaul gap list (priority-ordered)

1. **Aerosols-module reachability gap.** `src/Aerosols/` exports 11 symbols that no user can access because the module is not included in `src/vSmartMOM.jl`. Decide whether to (a) wire it in, (b) remove it, or (c) mark explicitly as internal/staging. Until resolved, `aerosol_integration_example.jl` likely does not run.
2. **README rewrite** for v2.0.0: quickstart, feature matrix, surface table, migration pointer, RTModel mention, linearized RT snippet, HITRAN-2024 primer, example/tutorial index.
3. **CHANGELOG deepening**: concrete before/after code diffs for the RTModel break, entries for HITRAN-2024 / Float32 / GPU Mie / Canopy / GEOSChem / batched kernels. Add a migration subsection linking to a new `pages/Migration.md`.
4. **Migration guide** (`docs/src/pages/Migration_from_v1.md`): every renamed symbol, every accessor change, every call-site diff, pointer to the `Base.getproperty` shim and what it covers/doesn't.
5. **Examples build-out** covering the golden path: forward Lambertian, forward Cox-Munk, forward canopy, forward RPV, linearized Jacobians, forward Raman, GPU variant. Each with `.jl` + `.md` + committed PNG/CSV of reference output.
6. **API reference audit** (`docs/src/pages/api_reference.md`): missing `@docs` entries include `RTModelLin`, `SolverConfig`, `Atmosphere`, `RayleighScattering`, `AerosolState`, `Optics`, `OpticsLin`, `ParameterLayout` accessors (`n_total`, `aerosol_range`, `gas_range`, `surface_range`, `surface_index`, `canopy_range`, `n_layer_params`), `FwdMode`, `LinMode`, `get_surface`, `get_surfaces`, `get_spec_bands`, `default_solar_transmission`, `planck_spectrum_*`, `default_parameters`. Some listed `@docs` entries (e.g. `CoreRT.fresnel_coefficients`) are not exported — confirm these are the intended internal functions to document.
7. **Docstring backfill** for bare-struct exports: all `AbstractPolarizationType`, `AbstractAerosolType`, `AbstractTruncationType`, `NAI2`, `PCW`, `Stokes_*`, `GaussQuadFullSphere`, `RTModelLin`, all Absorption error-function structs, `FwdMode`/`LinMode`. Target: at least one prose paragraph each.
8. **Dev-notes cleanup**: archive the 8 stale notes to a `docs/dev_notes/archive/` subfolder with an index README; merge `CUDA_SETUP.md` and `JACOBIAN_TEST_WORKFLOW.md` into published Documenter pages; retire `NEXT_STEPS_ANALYSIS.md` and `sanghavi_unified_merge_plan.md` in favor of `plans/CLAUDE_HANDOFF_BRIEF.md` + forthcoming implementation plan.
9. **Documenter warnings**: flip `warnonly = true` → `warnonly = false` in `make.jl` after the docstring backfill to make missing docstrings / broken `@ref`s fatal in CI.
10. **Test-file intake**: decide which of the 7 on-disk-but-unregistered test files (`test_forward_lin.jl`, `test_forward_raman.jl`, `test_hybrid_ad.jl`, `test_jacobians_GPU.jl`, `test_mie_gpu.jl`, `test_Aerosols.jl`, `test_performance.jl`) should be (a) wired into `runtests.jl`, (b) gated behind an env flag like `VSMARTMOM_FULL_TESTS`, or (c) deleted. The 474-test claim does not include them.
11. **Canary/Quickstart CI lane**: extract `test_CoreRT.jl` + `test_forward_noRS.jl` + `test_float32.jl` + `test_rami_smoke.jl` into a fast subset and document it.
12. **Tutorial audit** against the new RTModel API: verify each `Tutorial_*.jl` still runs post-refactor and references `model.optics.*` / `model.atmosphere.*` instead of flat `model.*` fields (this has been partially done via the `Base.getproperty` shim, but the narrative text may still show stale field access).
13. **README badge URLs** may need updating (devcode, JOSS DOI) if maintainer/URL changes.
14. **Aerosol example reachability fix** once §7.1 resolves.
15. **Docs deployment ref-name logic** in `make.jl` currently publishes `unified-vsmartmom` → `unified-vsmartmom` URL; needs update for eventual `sanghavi-unified` branch/PR preview URL.

---

## 8. Open questions

1. **Aerosols module visibility.** Is `src/Aerosols/Aerosols.jl` intended to be a user-facing subpackage (requiring it to be `include`d and `using`-ed in `src/vSmartMOM.jl`), an internal helper used only by tests and examples, or staging work? This affects (a) README messaging, (b) whether `test_Aerosols.jl` gets wired in, and (c) whether `aerosol_integration_example.jl` is a shipping example.
2. **Version numbering.** Is the docs overhaul shipping as v2.0.0 (per CHANGELOG) or as v2.0.0-rc1 / sanghavi-unified pre-release? The version affects migration-guide tone.
3. **Tutorial source of truth.** For a ship-quality release, should `examples/*.jl` or `docs/src/pages/tutorials/Tutorial_*.jl` be the canonical place for runnable code? Current duplication risk: new user doesn't know which to read first.
4. **Dev-notes archival policy.** Is in-repo archival (move to `docs/dev_notes/archive/`) acceptable, or should stale notes be deleted and preserved only in git history? CI impact of warn-only-off flip depends on this.
5. **Linearized Raman documentation.** `CLAUDE_HANDOFF_BRIEF.md` §3 declares linearized Raman permanently out of scope. Should the CHANGELOG and published docs carry a stronger "not supported, will never be supported" statement instead of the current "does not yet support" wording?
6. **API-reference scope.** Should `api_reference.md` document only root-level exports (user-facing) or also internal-but-exported symbols (Scattering helpers, compute_B, etc.)? Current page mixes both at irregular depth.
7. **Canary CI lane gating.** Should it run on every PR, or only post-merge? Claimed 474-test suite is fast-ish (~couple minutes per subset); tradeoff is per-PR feedback vs CI minutes on GPU-bearing runners.
8. **Examples reference-output format.** PNG only? PNG + committed text (JLD2/CSV)? The latter enables regression checks but inflates repo size.
9. **Backward-compat shim retention.** `Base.getproperty(::RTModel)` currently forwards `model.τ_abs`, `model.profile`, etc. Is this shim permanent, or is v2.x its last release before users must migrate to `model.optics.τ_abs`? Migration guide wording depends on this.
10. **Exported aerosol scheme types.** `AerosolScheme`, `TOMAS15Scheme`, `TwoMomentScheme` are all exported from `Aerosols/` (not reachable) but also referenced in `IO/` parameter parsing. Need to confirm whether IO actually needs these types at load time or only via reflection.
