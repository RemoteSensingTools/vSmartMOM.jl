# Release readiness audit — 2026-09-06

**Verdict: do not release or register this commit yet.** The integration is a useful release base, with substantial automated coverage, but reproducible numerical inconsistencies, a downstream compatibility regression, and a failing documentation build remain. Registration metadata and release notes also need reconciliation.

Candidate: `integration/surface-split-multisensor`, commit `b524a0d85ae36ab0eecf8b23a939f4d6ce62c456`. Review branch: `review/release-readiness`; worktree: `/tmp/vsmartmom-release-review`. The original checkout's unfinished merge was left intact. No solver fixes, tags, remote pushes, registration requests, or deployments were made during this audit.

The assessment covers installation/registry metadata, CI, executed CPU/CUDA checks, targeted numerical invariants, the public API and example workflows, documentation generation and consistency, source organization, and Jacobian performance. It is not a proof of every scientific configuration or an independent rederivation of the full radiative-transfer theory. The conventions page was read before assessing comparison guidance.

## Validation record

| Check | Result and scope |
|---|---|
| Exact candidate GitHub CI | All nine test jobs passed: Julia 1.10/1.11/1.12 on Linux/macOS/Windows. Taplo also passed. [Run 33457316128](https://github.com/RemoteSensingTools/vSmartMOM.jl/actions/runs/33457316128). No Documentation check ran for this commit. |
| Manifest-free candidate resolution and load | Passed on Julia 1.12.6 using registry dependencies. `pathof(vSmartMOM)` confirms this worktree's source. Shared package/artifact cache was used; this was not an empty-depot network test. |
| Dependency provenance | AtmosphericAbsorption 0.1.2 and CanopyOptics 0.2.0 resolved to General's registered tree hashes. CUDA resolved to 5.11.3. No dependency URL/path override was needed. |
| Local complete CPU suite | **Passed: 3,686 assertions, 14 entries reported as Broken/skipped, across 57 top-level testsets; process exit 0.** CUDA devices hidden; `PHASE1B_CPU=1` enabled the optional Raman reference comparison. VLIDORT baseline: 42/42 passed. JET warned about exceeding its advisory baseline (see below). |
| Shipped GPU runner | Failed immediately: five include-path errors before any test assertions. |
| GPU tests through audit harness | **20,291/20,291 passed**, in 6m50.7s on an NVIDIA A100. Harness changes only test include paths, not package code. Includes Mie, fused Raman, multisensor, forward Raman and GPU Jacobians. |
| Quickstart | Passed: `R = 0.02848592913380708`, `T = 0.1557755644589617`, both `(1,1,1)`. Explicit external-solar TOA gives the same R in this scene. |
| Documentation | **Failed** in Documenter before rendering: missing `IntensityConvergence` documentation and unresolved `IntensityConvergence`/`StokesConvergence` references. Julia 1.12.6, local deployment disabled. |
| Static relative Markdown targets | No missing relative file targets found in README or `docs/src` by the path scan. This does not validate every external URL or heading anchor. |
| npm dependency audit | `npm ci` succeeded; `npm audit` reports 6 findings: 5 moderate, 1 high (Vite, transitive). Scope is docs tooling, not Julia RT execution. |
| Machine portability scan | No `/home/`, `/Users/`, `/net/`, or `/mnt/` paths found in tracked source, public docs, or shipped configs. |

The machine's default `julia` command points to a broken juliaup backup launcher. Tests used `/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia` directly. That launcher problem is environmental, not a candidate defect. Metal hardware was not available; macOS CPU CI does not establish Metal execution coverage.

## Findings requiring action before a release

### R1 — Forward and full-Jacobian aerosol models implement different reference normalization (P1)

Locations: `src/CoreRT/tools/model_from_parameters.jl:480`, `src/CoreRT/tools/lin_model_from_parameters.jl:391`, `src/IO/Parameters.jl:1209`.

The forward builder evaluates reference extinction using the common configured `scattering.n_ref`. The full linearized builder instead uses each aerosol's own refractive index. This changes the forward state before any derivative is evaluated. The parser defaults the common reference index to the **first** aerosol's index, so the mismatch is not confined to an unusual explicit override: it occurs in ordinary multispecies scenes with differing indices.

Reproduction with two aerosols, indices 1.3 and 1.5, each `τ_ref=0.04`, and no explicit `n_ref`:

- First species' AOD agrees.
- Second species' first-wavelength AOD: forward `0.10971730066954853`; full linearized model `0.040000000000000036`.
- Resulting maximum radiance difference, normalized by maximum absolute forward radiance: **13.27%**.
- Running the forward and derivative solvers on the *same linearized-built model* agrees to `2.4e-16` relative, isolating the discrepancy to model construction.

An explicit `n_ref=1.5` case with aerosol index 1.3 also reproduces the problem: AOD `0.01458` versus `0.04`, with **4.07%** radiance difference. Setting `compute_aerosol_microphysics_jacobians=false` restores the forward builder's normalization in this case. A computation-saving switch should not change the simulated forward physics.

Required remedy: define reference-index semantics once and share the forward calculation and its derivative across both builders and batch updates. Test ordinary multispecies defaults, explicit distinct `n_ref`, and finite differences that respect whether the reference index is fixed or linked. Merely adjusting a tolerance will not resolve this mismatch.

Evidence: [probes.jl](release_audit_2026_09_06/probes.jl), [extra_probes.jl](release_audit_2026_09_06/extra_probes.jl).

### R2 — Cox–Munk adds glint with no illumination and violates solar scaling (P1)

Locations: `src/CoreRT/rt_run.jl:809`, `src/CoreRT/rt_run_split.jl:435`, `src/CoreRT/Surfaces/coxmunk_surface.jl:509`.

`apply_ss_correction!` adds an exact-minus-Fourier glint correction without an incident-source factor. It runs whenever the surface is Cox–Munk and SFI is enabled, including a `NoSource()` scene. The source comment acknowledges the unit-beam limitation, but the public API neither rejects nor corrects nonunit illumination.

Small allowed configuration: 3 streams, Stokes IQU, two layers, wind 5 m/s, SZA/VZA 30°, relative azimuth 0°:

| Source | TOA I | TOA Q |
|---|---:|---:|
| Unit solar beam | 0.0166148143 | 0.0243883684 |
| Double solar beam | 0.0288913891 | 0.0212812163 |
| `NoSource()` | 0.0043382394 | 0.0274955205 |

The dark scene should return zero. The vector-norm scaling defect `maximum(abs,R2-2R1)/maximum(abs,2R1)` is **56.37%** in this probe. The deliberately low stream count exposes the correction; this is not a claim about all ocean-scene errors.

Required remedy: make the correction use the same source normalization as the surface operator, including the zero-source case, and apply the same contract to cache replay. Add dark-scene and scalar/spectral source-linearity regressions.

Evidence: [probes.jl](release_audit_2026_09_06/probes.jl).

### R3 — Cox–Munk linearized radiance omits the forward glint correction (P1)

Locations: `src/CoreRT/rt_run.jl:809`, `src/CoreRT/rt_run_lin.jl:515`.

The forward driver applies `apply_ss_correction!` after the Fourier sum; the linearized driver does not. This is separate from R2's source scaling and persists for the default unit solar beam.

On the same linearized-built model used by both solvers:

- Forward: I `0.0166148143`, Q `0.0243883684`.
- Linearized forward output: I `0.0122765749`, Q `-0.0031071521`.
- The discrepancy is the omitted correction, not different Mie optics or model construction.

The `TMSCorrection` rejection in the linearized driver does not protect this case: Cox–Munk's built-in correction is separate from the configured aerosol `TMSCorrection` strategy. The branch's own `forward_lin_parity_policy.md` requires matching behavior or explicit rejection.

Required remedy: propagate the correction and relevant tangents consistently, or explicitly reject the unsupported ocean Jacobian mode until it is implemented. Test full forward/linearized radiance equality and wind/optical-depth finite differences near glint.

Evidence: [coxmunk_parity.jl](release_audit_2026_09_06/coxmunk_parity.jl).

### R4 — Importing vSmartMOM breaks valid NNlib calls in other packages (P1)

Locations: `src/CoreRT/CoreRT.jl:42`, `src/CoreRT/tools/cpu_batched.jl:65`.

The package adds a method to the external `NNlib.batched_mul` function for ordinary Julia arrays of BLAS element types. Its batch-dimension equality assertion removes NNlib's supported singleton-batch broadcasting behavior throughout the Julia process.

Reproduction:

```julia
using NNlib
A = ones(2, 2, 3); B = ones(2, 2, 1)
size(NNlib.batched_mul(A, B))  # (2, 2, 3)
using vSmartMOM
NNlib.batched_mul(A, B)        # AssertionError: batch dim mismatch: 3 vs 1
```

Required remedy: put the optimized implementation behind a package-owned wrapper/backend operation, or at minimum preserve the full NNlib input contract. Include a downstream coexistence regression. Passing RT tests with equal-sized batches does not test this behavior.

Evidence: [nnlib_compat.jl](release_audit_2026_09_06/nnlib_compat.jl).

### R5 — Documentation cannot build at the candidate commit (P1 release gate)

Locations: `src/CoreRT/types.jl:1045`, `src/CoreRT/types.jl:1079`, `docs/src/pages/api/core_rt.md:35`, `.github/workflows/Documentation.yml:4`.

The local strict build terminates with `makedocs encountered errors [:docs_block, :cross_references]`. The bare `raw` string blocks preceding the convergence types are not registered as the docstrings Documenter expects; `StokesConvergence` also needs an API documentation entry. Attach the docstrings explicitly and include both public strategies in the manual. Keep strict checking enabled.

The docs workflow's push branch list excludes the integration branch. That explains why green candidate CI did not expose this failure. Require a successful Documentation check on the actual release candidate, through a PR or an appropriate branch trigger.

### R6 — Release identity and registry history disagree (release gate)

Locations: `Project.toml:3`, `CHANGELOG.md:3`, `CHANGELOG.md:61`.

The candidate still declares `2.1.0`. Tag `v2.1.0` already identifies commit `1680085d` from May 8, while the changelog describes an unreleased externalAbsorption branch targeting `2.2.0`. General currently contains releases only through **1.1.0**. The repository's latest GitHub release is also 1.1.0; a Git tag alone is not a registry release. [General version history](https://raw.githubusercontent.com/JuliaRegistries/General/master/V/vSmartMOM/Versions.toml).

Required remedy: choose a new, noncolliding version, synchronize metadata and release notes, and explain migration from the last registered version. If retaining the planned `2.2.0` number, address General's sequential-version AutoMerge rule explicitly: no 2.0/2.1 versions are currently registered. This is an AutoMerge/manual-review issue, not a blanket ban on registration. Do not move or reuse the existing 2.1.0 tag. [RegistryCI guidelines](https://juliaregistries.github.io/RegistryCI.jl/stable/guidelines/).

The release is breaking relative to registry 1.1.0: model/result APIs, configuration semantics, minimum Julia version, and scientific behavior have changed. The migration narrative should start from that baseline, not assume that users received an unregistered 2.0 release. Also assess the compatibility promise made by the existing public 2.1.0 tag when choosing a version.

### R7 — The shipped GPU runner cannot execute its tests (P2)

Location: `test/local/gpu/runtests.jl:33` and the four subsequent includes.

`include("local/gpu/test_mie_gpu.jl")` inside the GPU runner resolves relative to that source file, producing `test/local/gpu/local/gpu/test_mie_gpu.jl`. Changing the working directory does not change Julia's source-relative include behavior. The documented command fails with five errors. Use paths based on `@__DIR__` or bare sibling filenames.

The audit harness demonstrated that all **20,291** assertions pass once these paths are corrected. That is encouraging evidence for the tested CUDA kernels, but the shipped command must work. Hardware-enabled CI or an explicit release GPU job should fail on kernel compilation/scalar-indexing regressions; `test_jacobians_GPU.jl:91` currently turns several such failures into skips even on functional hardware.

Evidence: [gpu_runner.jl](release_audit_2026_09_06/gpu_runner.jl).

### R8 — Spectral Greek arrays cause quadratic cache growth in the forward solver (P1)

Locations: `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:203`, `src/Scattering/compute_Z_matrices.jl:53`.

The forward cache chooses the angular expansion order with `length(greek_coefs.β)`. For matrix-valued coefficients this counts both angular and spectral axes. `ZMomentTables` then constructs six arrays quadratic in that inflated order. This affects forward solves on models with wavelength-dependent Greek arrays, including the linearized-built model used for forward/full parity in this audit.

Reproduced shape: `β` is **5×512**, but the cache requests angular length **2,560**. Its six base tables alone require **1,887,436,800 bytes**, versus **7,200 bytes** at the actual five-term angular order. A forward timer sample assigns 3.14 s and 2.35 GiB of cumulative allocations to cache initialization. Increasing the number of wavelengths should not increase scatterer-independent angular tables this way.

Using the existing `_Z_TABLES_ENABLED` diagnostic switch to disable this cache reduces the warmed 512-wavelength CPU forward solve from 3.53 s to **0.460 s**, with **bit-identical R and T**. This is a diagnostic bypass, not an implemented release fix. Correct the dimension calculation to use the angular axis and cover vector and matrix Greek arrays in cache-size and radiance-parity tests. Audit equivalent length calculations for the same shape assumption.

Evidence: [cache_dimensions.jl](release_audit_2026_09_06/cache_dimensions.jl), [cache probe output](release_audit_2026_09_06/evidence/cache-dimensions.log), and the [profiling report](jacobian_bottlenecks_2026-09-06.md).

## Documentation, onboarding and API consistency

The manual's task-first entry points, Concepts arc, schema pages, named observer results, and theory map are valuable. The tiny CPU quickstart is reproducible and independent of local science datasets. These are strong foundations; the main weakness is conflicting versions of the same contract.

| Issue | Evidence | Required correction |
|---|---|---|
| README Jacobian example throws | `README.md:95` dereferences `params.scattering_params.rt_aerosols`, but the chosen ocean config has `scattering_params === nothing`. Reproduced `FieldError`. | Guard the absent scattering block and choose a validated Jacobian scene; also address R3 before advertising this ocean workflow. |
| README install-to-example path is fragile | `README.md:84` uses repository-relative `config/quickstart.yaml` immediately after `Pkg.add` instructions. | Use `joinpath(pkgdir(vSmartMOM), "config", "quickstart.yaml")`, as the manual already does. |
| Azimuth instructions contradict the conventions page | `config/quickstart.yaml:29`, `config/ocean_coxmunk.yaml:54`, and many test configs describe 0°/180° with the opposite same-side/opposite-side wording to `docs/src/pages/conventions.md:88`. | Establish one authoritative diagram/text, remove copied contradictory banners, and validate any cross-code conversion instructions. This audit did not independently resolve the physical convention. |
| Retired truncation angle described as active | `CHANGELOG.md:43` and the “Correctness Fixes” section of release notes say a nonzero `Δ_angle` is applied; `src/Scattering/types.jl:292` warns and forces it to zero. | Update current release notes and migration guidance. Do not imply existing nonzero-angle workflows preserve their numerical behavior. |
| CUDA described as a weak dependency | `docs/src/pages/concepts/07_architecture.md:116` versus `Project.toml:8`; CanopyOptics also depends on CUDA. | Explain the actual install/load contract, or complete the weak-dependency work across the dependency chain. CPU use works, but installation still includes CUDA. |
| External-solar default is stale in agent/developer onboarding | `AGENTS.md:44,222`, `CLAUDE.md:60` versus constructor defaults at `model_from_parameters.jl:257` and `lin_model_from_parameters.jl:100`. | Both constructors default to `false`. `rt_run_toa(model_from_parameters(params))` throws unless explicitly opted in. Bring onboarding into agreement with the already-correct newer manual pages. |
| Schema rejects new numerical strategies | `schemas/vsmartmom-parameters.schema.json:196` sets `additionalProperties=false` but lists only four numerics fields. Parser accepts Fourier convergence/tolerance/guard/consecutive settings and `ss_correction` at `src/IO/Parameters.jl:1408`. | Add these fields and accepted values to the JSON Schema and schema docs; validate actual configurations, not only key presence. |
| Performance promise is too broad | `docs/src/pages/concepts/06_linearization.md:112,141`, `AGENTS.md:132`, `doubling_lin.jl:48` promise full Jacobians at `<2×` forward cost. | Replace with scoped, reproducible benchmark results. Shared inverses save work, but per-parameter products, allocations and transfers remain. See the [profiling supplement](jacobian_bottlenecks_2026-09-06.md). |
| Version/WIP narrative is inconsistent | Changelog refers to 2.0/2.1 and a different branch; release notes say thermal emission is not implemented despite a production source and tests. | Write one candidate-specific release/migration summary; separate supported, experimental and intentionally rejected combinations. |

The tutorial generator uses ordinary `julia` fences (`docs/make.jl:591`). A successful docs build alone will not execute every displayed example. Promote a small set of installation, forward, Jacobian, batch-update and source examples to executable smoke tests or doctests.

## Jacobian performance

The [dedicated profiling report](jacobian_bottlenecks_2026-09-06.md) contains warmed forward/full/selected measurements, CPU sampling, allocation attribution, backend comparisons, reproducible scripts and a prioritized optimization plan. For a 64-wavelength, five-layer, 6×6 operator scene, fourteen-column Jacobians take 1.568 s on one CPU thread versus 0.0888 s forward, allocating 1.55 GB cumulatively. Four-column selection reduces the solve to 0.450 s and 0.440 GB, with exact agreement against the selected full outputs. Doubling and interaction account for about 90% of measured CPU section time. Prioritize their allocating per-parameter matrix-product chains and derivative-slice copies; broad claims about shared inverse cost do not describe this workload.

## Registry and release automation

- Package name/UUID remain consistent. The root Apache-2.0 license is present. Runtime dependencies have bounded compatibility entries, and registry-only resolution succeeded.
- Current direct dependency constraints allow CUDA 6 and ForwardDiff 1, but CanopyOptics 0.2.0's registered constraints require CUDA 5 and ForwardDiff 0.10. Those advertised wider combinations are not reachable with the present dependency graph. Coordinate dependency compatibility rather than treating a broad root entry as tested support. [CanopyOptics compatibility](https://raw.githubusercontent.com/JuliaRegistries/General/master/C/CanopyOptics/Compat.toml).
- `.github/workflows/TagBot.yml:11` uses only `GITHUB_TOKEN`. Tags created with that configuration do not trigger the downstream tag-based docs workflow. Configure the supported deploy-key mechanism if automatic versioned docs are part of the release, and modernize the schedule-only workflow using the current TagBot template. Actual repository secret/permission configuration was not inspected. [TagBot setup](https://github.com/JuliaRegistries/TagBot#setup).
- Registration requires a deliberate Registrator action; TagBot reacts to registration and does not perform it. Prepare release notes explaining breaking changes before that action. [Registrator](https://juliaregistries.github.io/Registrator.jl/stable/).
- `npm audit` findings need triage against the docs toolchain's actual usage. Do not equate the Vite development-server advisory with a vulnerability in the Julia solver or a confirmed compromise of the published static site. Prefer reproducible `npm ci` in docs CI once the lockfile is accepted.

## Code readability and maintainability

The architectural direction is sound: backend extensions, Stokes types, source and convergence strategies, optical-property composition, and explicit selective Jacobian layouts make the solver's intent visible. Reused forward inverses, finite-element exponentials using stable forms, and source-boundary separation are useful design choices. Existing tests exercise many meaningful invariants rather than only implementation details.

The most valuable cleanup is to reduce duplicate scientific definitions. Forward/linearized/reference normalization already diverged (R1), and surface handling is mirrored across full, split and linearized drivers (R2/R3). Shared helpers for these calculations would improve correctness and readability together. Preserve the mathematical identifiers and float-type generic interfaces while doing so.

The tracked package contains 154 Julia source/extension files and about 54,473 lines; seven files exceed 1,000 lines. Large files alone are not blockers, but `model_from_parameters`, `update_model`, `rt_run` and the type collection mix several independent responsibilities. Extract cohesive profile/absorption/aerosol construction and result/source assembly helpers, with parity tests, instead of undertaking a wholesale rewrite during release stabilization.

There is also substantial historical prose embedded in production files: phase/PR labels, past review discussions, obsolete line anchors, duplicate formulas and superseded performance claims. Keep the current contract and numerical rationale close to the implementation; move historical discussion to development notes. `rt_run_lin.jl`'s opening bibliographic entries should be checked against `theory_references.md` rather than copied as authoritative citations.

Test-maintenance gaps:

- `test_fourier_convergence.jl` is included twice by `test/runtests.jl`.
- Aqua is run in both `test_quality.jl` and `test_aqua.jl`; ambiguity and persistent-task checks are disabled in both.
- `test_aqua.jl` says ambiguities are tracked by `test_quality.jl`, but that file also disables them.
- JET is advisory unless `VSMARTMOM_JET_STRICT=1`; the local run reports **280 findings against a baseline of 232**, so its default successful test does not mean the static-analysis baseline is met. Pin a dedicated analysis environment before making it a strict regression gate.
- Selected GPU tests catch compiler and scalar-indexing errors as skips. Distinguish absent hardware from failure on hardware that is available.
- CI job durations for this candidate ranged from roughly 33 to 81 minutes. Separate fast invariants, full numerical regression and hardware jobs to improve feedback without dropping release coverage.

Tracked content totals about 58 MiB across 1,027 files: approximately 32 MiB of tests, 11 MiB of docs and 8 MiB of sandbox workflows. No individual tracked file exceeded 4 MiB in this snapshot. Most payload is explainable by scientific fixtures; size is a secondary packaging concern. Keep external truth products outside the package and make required fixture provenance explicit.

## Scope of later Sanghavi work

This audit intentionally did not merge `suniti_multi_sensor` into the candidate. Its September work includes SIF truth normalization/restart changes, tapered CO2 campaigns and packaging. For example, `a2968e26` changes RRS_XCO2 workflow scripts and tests rather than the core `src/` tree. If the new release claims that the bundled OCO/SIF workflow is production-ready, review and port the relevant corrections with their validation. Do not infer that the integration branch already contains them, and do not replace the integration base with the less-complete development branch.

## Recommended release sequence

1. Resolve R1–R4 and R8 with targeted invariants, cache-dimension checks and finite-difference/downstream coexistence regressions.
2. Repair docs generation and GPU test invocation, then require both checks on the release commit.
3. Reconcile source/default/azimuth/truncation documentation, fix the installed-package examples and complete the numerical-strategy schema.
4. Use the Jacobian profile to prioritize measured bottlenecks; publish performance claims only for named workloads, parameter counts, precision and hardware. Avoid a broad performance rewrite that delays correctness fixes.
5. Select the version and migration baseline, reconcile General's history/AutoMerge requirements and repair versioned-docs automation.
6. Run final package tests, the applicable CPU/CUDA/Metal release matrix, strict docs, and a fresh-environment install on the fixed release commit. Prepare the registration request and changelog for final review, then register/tag through the agreed release process.

Selected logs and profiles are preserved in [the evidence bundle](release_audit_2026_09_06/evidence/); the complete raw evidence directory is `/tmp/vsmartmom-release-evidence/`. Reproduction scripts are retained beside this report in [release_audit_2026_09_06](release_audit_2026_09_06/). They are standalone audit tools, not additions to the default package test suite. Run them with the candidate environment and run package tests from `test/` as required by the repository conventions.
