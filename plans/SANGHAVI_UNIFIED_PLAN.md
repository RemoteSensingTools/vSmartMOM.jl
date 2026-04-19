# `sanghavi-unified` merge plan

*Draft 2026-04-19. For review by Codex. Supersedes `docs/dev_notes/sanghavi_unified_merge_plan.md` (Christian's 2026-03-21 plan), which was written before the sanghavi branch owner clarified scope. Phases 2–4 below are faithful to Christian's technical content; Phase 0 and Phase 5 are new; Phase 5 of the original (linearized Raman) is removed from scope entirely.*

Repository: https://github.com/RemoteSensingTools/vSmartMOM.jl

## 1. Scope and goals

Create a new branch `sanghavi-unified` that combines the engineering quality of `unified-vsmartmom` (RTModel hierarchy, linearized RT, GPU infrastructure, 474-test suite, Cox-Munk/canopy/RossLi/RPV surfaces, Aerosols module, IO module, Float32 support, HITRAN version tracking) with the newer Raman + SIF + application physics from `sanghavi` (EMIT, Balsamic, OCO forward simulations, O2 A/B band Raman+SIF grids, CarbonI prototype, hydrolight comparisons). The merged branch eventually replaces `main`.

The merge has three tightly coupled outcomes:

1. **Physics preservation.** The Raman + SIF + application workloads that run on `sanghavi` today continue to run on `sanghavi-unified`, within a tolerance bar defined per quantity at execution time. Acceptance is script-based, not test-unit-based: specific user-supplied scripts must produce equivalent results.
2. **Forward Raman performance.** The GPU Raman path is rewritten to eliminate ~19,000 per-run GPU allocations and replace sequential-over-Raman-index kernels with gathered batched ones. Targets are the *forward* Raman path only — linearization remains elastic-only permanently.
3. **Clean release-quality documentation and examples.** Every module has usable docstrings, every supported use case has a runnable example in `examples/` with a narrative `.md` companion, the top-level README is a fast on-ramp, and a dedicated migration guide helps the sanghavi-branch author (and future contributors arriving from older forks) orient in the new code.

## 2. Non-goals (explicit)

- **Linearized Raman is out of scope.** The unified linearized path remains `RS_type = noRS()` only, as it is today. No `AddedLayerRS`-style linearized constructors will be added; no `rt_kernel_lin.jl` Raman branch. This is a deliberate scope decision, not a deferral — it is not planned to land in a later PR either.
- **Elastic-path performance retuning.** The existing `⊠` notation in the elastic path works well with Julia's JIT and is valuable for physics readability. Phase 2 and Phase 3 touch the Raman path only.
- **New physics capabilities.** The merge preserves and consolidates existing physics; it does not add new scattering regimes, new surface types, or new retrieval operators.
- **Format/lint churn.** No mass reformat or lint sweep. Code changes are tied to concrete merge or performance objectives.

## 3. Phase structure (summary)

| Phase | Scope                                                 | Risk   | Key outcome                                    |
|-------|-------------------------------------------------------|--------|-----------------------------------------------|
| 0     | Pre-merge reconciliation with sanghavi's recent work  | Low    | No duplicated work; `InteractionWorkspace` reuse resolved |
| 1     | Branch creation + sanghavi port + migration guide     | Low    | `sanghavi-unified` branch exists and tests pass |
| 2     | Eliminate allocations in inelastic path               | Low    | Bit-exact; ~0 GPU allocs inside Raman loops    |
| 3     | Flatten 4D → 3D with `RamanIndexMap`                  | Medium | ~nRaman× GPU speedup; same numerics (tolerance)|
| 4     | Reduce GPU synchronization                            | Low    | <10 syncs/run; further speedup                 |
| 5     | Documentation and examples overhaul                   | Low    | Release-quality docs; runnable examples; API ref |

Phases 0–4 are serial. Phase 5 runs continuously from mid-Phase 1 onward, folding each completed capability into the docs as it lands, with a consolidation pass after Phase 4.

## 4. Phase 0 — Pre-merge reconciliation

**Problem:** The sanghavi branch has recent commits that overlap with Christian's Phase 2–3 design. If we jump straight to Phase 1 without reconciling, Phase 2 will reinvent work that already exists.

**Overlapping commits on `sanghavi`:**

- `854b44c` *Add Float32 support, InteractionWorkspace, and optimize get_n₀_n₁*
- `d75dacb` *Add per-direction CPU-staging for interaction workspace*
- `3e9926f` *Add Raman optimization test infrastructure*
- `e89ec1c` *Update baseline script and add batched ops benchmark*
- `083353b` *Update lin+Raman* (physics content still relevant even though we're not porting linearized Raman)

**Tasks:**

- 0a. Read each of the above commits on `sanghavi`. Produce a short note (committed as `docs/dev_notes/sanghavi_optimization_inventory.md` on the new branch) listing, for each commit, what it changed, what abstraction it introduced (e.g., sanghavi's `InteractionWorkspace`), and how it relates to Christian's proposed `InelasticWorkspace` in Phase 2a.
- 0b. Decide per-item whether to (i) lift the sanghavi implementation wholesale, (ii) use it as a design input for a rewritten version on unified's API, or (iii) discard in favor of Christian's approach. Record the decision and rationale in the same note.
- 0c. If `InteractionWorkspace` on sanghavi is structurally compatible with Christian's `InelasticWorkspace` design, rename/adapt rather than duplicate. If incompatible, pick one — whichever survives becomes the name used throughout Phase 2.
- 0d. Feed the reconciliation note to Christian for sign-off before starting Phase 1. This is a coordination step, not a code step; target one working day.

**Exit criterion:** A single-page note exists in the new branch's `docs/dev_notes/` that explains, per sanghavi overlap-commit, what happens to it in the merge. Christian has seen and commented.

## 5. Phase 1 — Branch creation and sanghavi compatibility

Drawn from Christian's Phase 1, with deliverables sharpened.

- 1a. **Branch.** Create `sanghavi-unified` from `unified-vsmartmom` HEAD.
- 1b. **Port remaining sanghavi physics content** to unified's API:
  - Single-scatter approximation kernel → adapt to `RTModel` conventions, land in `src/CoreRT/CoreKernel/rt_kernel_ss.jl` and `interaction_ss.jl`.
  - SIF infrastructure: confirm that the `RS_type.SIF₀` plumbing used by `creategrid_O2Aband_RamanSIF.jl` is present or port it. Commit `src/SIF_emission/sif-spectra.csv` (and any other real data files — only `.dat~` backups are in the repo today). If SIF is intended to stay outside the `RS_type` system in the long term, decide that now and document it.
  - EMIT / Balsamic application scripts → `examples/EMIT/`, `examples/Balsamic/` (adapted to the new API; see Phase 5 for the narrative companions).
  - `src/Testing/perturb_parameters.jl` (sanghavi's finite-difference Jacobian utility) → commit to unified under `test/utils/perturb_parameters.jl` so the EMIT script and any other consumer can include it portably.
  - Verify remaining physics diffs in `elemental_inelastic.jl`, `raman_atmo_prop.jl`; fold in anything sanghavi has that unified lacks.
- 1c. **Migration guide.** Write `docs/dev_notes/migration_from_sanghavi.md` documenting the old→new API, covering:
  - Struct changes: `vSmartMOM_Model` → `RTModel{ARCH,FT}` with sub-structs.
  - Field access: `model.τ_abs` via `Base.getproperty` → `model.optics.τ_abs`; `model.obs_geom` → `model.geometry`.
  - Function signature changes: `model_from_parameters(LinMode(), params)` returns `(model, lin_model)`.
  - Naming: `@unpack` → `(; field) = struct`, `J₀⁺/J₀⁻` → `j₀⁺/j₀⁻` (lowercase for `AddedLayer`), `FT<:Union{AbstractFloat,Dual}` → `FT<:Real`.
  - New accessors: `CoreRT.polarization_type(m)`, `float_type(m)`, `n_aerosols(m)`, `get_surface(m, i)`, `get_spec_bands(m)`.
  - Removed: `vSmartMOM_Model`, `rt_run_bck.jl`, `types_inelastic.jl`.
  - Include **every** converted sanghavi example script as a worked old-vs-new pair in the guide.
- 1d. **Sanghavi script conversion.** Convert the script list the sanghavi-branch owner provides (see Phase 5 open-questions and Claude Code-at-execution instructions in §10). For each converted script, capture baseline outputs (R, T, ieR, ieT, dR/dT where applicable) into `test/regression/baselines/` and wire a regression test.
- 1e. **Cleanup in new branch:**
  - Drop sanghavi-side hygiene detritus from the merged branch: top-level PNGs, `Manifest copy*.toml`, `shettle_*.dat~`, `notes.txt`, `tmpfile.jld`, `SSA.png`, `SolarSpectrum_kurusz.png`, backup scripts like `rt_run_bck.jl`, `inelastic_helper_old.jl`, `src/Inelastic/src/prototype.jl`, `src/Inelastic/src/plots/`.
  - Replace absolute paths (`/home/sanghavi/code/github/...`) in ported scripts with project-relative references.
  - Remove `@show` debug statements from interaction kernels.
  - Remove commented-out dead code blocks (e.g., `#bla`, the `repeat(...) ⊠ reshape(...)` attempt in `doubling_inelastic.jl`).
  - Convert lingering `@unpack` → `(; field) = struct`.
- 1f. **Forward-Raman smoke tests.** Add minimal-signature tests ensuring `rt_run_test(RS_type, model, iBand)` returns finite, nonzero R, T, ieR, ieT for each `RS_type ∈ {RRS, VS_0to1, VS_1to0, RRS_plus, VS_0to1_plus, VS_1to0_plus}` and at least one polarization and one quadrature configuration. These are regression guards for Phases 2–4.

**Exit criterion:** On `sanghavi-unified`, (i) full 474-test suite passes, (ii) all ported sanghavi scripts run and produce outputs matching their sanghavi-branch baselines to the tolerances set at execution time (see §8 and §10), (iii) the migration guide renders and the worked examples it references actually run.

## 6. Phase 2 — Eliminate allocations in the forward Raman path

*Objective:* drop per-run GPU allocations inside the forward Raman path to ~0, without changing numerical outputs.

*Scope:* inelastic path only. Elastic `⊠` notation in `doubling.jl`, `interaction.jl`, `elemental.jl` is untouched.

*Grounding:* `docs/dev_notes/raman_gpu_optimization.md` identifies ~19,000 GPU allocations per Raman run. Breakdown:

- `interaction_inelastic.jl:293-305` allocates 13 arrays via `similar()` and zeroes them per call; ~640 calls/run → ~8,320 allocs.
- `doubling_inelastic.jl:62-81, 97-110` creates `⊠`-expression temporaries inside `for Δn ∈ 1:nRaman` nested in `for n ∈ 1:ndoubl`; ~360–3,600 depending on nRaman.
- `ext/gpu_batched_cuda.jl:34, 63-64` has `batch_inv!` allocating pivot arrays per call; ~2,000 allocs.
- Miscellaneous `⊠` temporaries in nested expressions; ~5,000+.

*Approach:*

- 2a. **Reconcile workspace abstraction with Phase 0 decision.** Either (i) adopt and extend the `InteractionWorkspace` already added on sanghavi (commit `854b44c`), or (ii) introduce Christian's proposed `InelasticWorkspace`. The two designs differ in scope — sanghavi's focuses on interaction temporaries; Christian's also covers doubling. Proposed reconciliation: keep sanghavi's name (`InteractionWorkspace`) and widen it to cover doubling buffers too. Fields (drawn from Christian's Phase 2a proposal):

  ```julia
  mutable struct InteractionWorkspace{FT, AT3, AT4}
      # Doubling temporaries (3D elastic-sized)
      gp_refl::AT3
      tt_gp::AT3
      J₁⁺::AT3; J₁⁻::AT3                # (nQ, 1, nSpec)
      ieJ₁⁺::AT4; ieJ₁⁻::AT4            # (nQ, 1, nSpec, nRaman)
      # Interaction temporaries (ScatteringInterface_11)
      tmp_inv::AT3
      T_inv::AT3                         # T01_inv or T21_inv
      tmpJ₀⁺::AT3; tmpJ₀⁻::AT3
      tmpR⁻⁺::AT3; tmpR⁺⁻::AT3
      tmpT⁺⁺::AT3; tmpT⁻⁻::AT3
      tmpieJ₀⁺::AT4; tmpieJ₀⁻::AT4
      tmpieR⁻⁺::AT4; tmpieR⁺⁻::AT4
      tmpieT⁺⁺::AT4; tmpieT⁻⁻::AT4
      # Batched_mul output buffers (reused across ⊠ calls)
      buf3d_a::AT3; buf3d_b::AT3; buf3d_c::AT3
      # batch_inv! pivot/info (shared with elastic RTWorkspace if present)
      pivot::AbstractMatrix{Cint}
      info::AbstractVector{Cint}
  end
  ```
- 2b. **Allocate in `model_from_parameters()`.** In `src/CoreRT/tools/model_from_parameters.jl`, construct a single `InteractionWorkspace` alongside the elastic `RTWorkspace` when `RS_type` is not `noRS`. Thread it through the RT pipeline (stored on a per-run context struct, not on the immutable model).
- 2c. **Refactor `doubling_inelastic.jl`** to read from the workspace instead of calling `similar()`. Replace `⊠` expressions inside the hot loop with `batched_mul!` writing into workspace buffers. Zero workspace buffers only where logically required; skip if the write-pattern guarantees overwrite.
- 2d. **Refactor `interaction_inelastic.jl`** analogously. Replace the 13 `similar() + .=0` at the entry of `ScatteringInterface_11` with workspace fields plus explicit zeroing.
- 2e. **Pass the pivot/info workspace to `batch_inv!`** in the inelastic path. The workspace-aware `batch_inv!(X, A, ws)` already exists in `ext/gpu_batched_cuda.jl`; the inelastic path just never calls it.

*Files touched:* `src/CoreRT/types.jl`, `src/CoreRT/tools/model_from_parameters.jl`, `src/CoreRT/CoreKernel/doubling_inelastic.jl`, `src/CoreRT/CoreKernel/interaction_inelastic.jl`, `src/CoreRT/rt_run.jl`.

*Acceptance:* (i) all existing tests pass with outputs bit-exact to the pre-Phase-2 baseline — no physics change; (ii) `@allocated` / `CUDA.@time` inside the doubling and interaction loops reports ~0 allocations; (iii) wall-clock speedup is measurable even without Phase 3.

*Known risk:* an "overwrite-only" invariant on workspace buffers can be subtly violated by a future contributor reading a buffer before writing. Mitigation: add assertions in debug builds and document the invariant on `InteractionWorkspace`.

## 7. Phase 3 — Flatten 4D → 3D with `RamanIndexMap`

*Objective:* replace the `(nQ, nQ, nSpec, nRaman)` inelastic arrays with `(nQ, nQ, nTotal)` where `nTotal = Σ_{valid (n₁, Δn)} 1`. Replace per-Δn sequential BLAS calls with a single gathered batched operation. Target ~nRaman× GPU speedup in doubling/interaction.

*Approach:*

- 3a. **Introduce `RamanIndexMap`** in `src/Inelastic/types.jl`:

  ```julia
  struct RamanIndexMap{IT<:AbstractVector{Int}}
      nTotal::Int                     # Σ valid (n₁, Δn) pairs
      el_n₁::IT                      # elastic index at scattered wavelength
      el_n₀::IT                      # elastic index at incident wavelength
      ranges::Vector{UnitRange{Int}}  # slice of nTotal per Δn
  end
  ```

  Built from `i_λ₁λ₀` + `get_n₀_n₁()` during `getRamanSSProp!()` in `src/Inelastic/inelastic_helper.jl`. Reuses sanghavi's `get_n₀_n₁` optimization from commit `854b44c`.
- 3b. **Flatten `AddedLayerRS` / `CompositeLayerRS`** fields: `ier⁻⁺`, `ieR⁻⁺`, `iej₀⁺`, `ieJ₀⁺`, etc. become 3D (`nQ, nQ, nTotal`).
- 3c. **Replace per-Δn loops** in `doubling_inelastic.jl`, `interaction_inelastic.jl`, `elemental_inelastic.jl` with:
  - A `gather!(dst, src, idx)` helper for pulling `r⁻⁺[:, :, idx.el_n₁]` and `r⁻⁺[:, :, idx.el_n₀]` slabs.
  - One batched `⊠` (or `batched_mul!`) per logical multiplication, operating on the gathered flat 3D arrays.
- 3d. **Update `default_matrix_ie()`** and helper constructors in `src/CoreRT/tools/rt_helper_functions.jl` to produce flat 3D allocations.
- 3e. **Memory budget.** The gather buffers add roughly `nQ² × nTotal × sizeof(FT)` bytes. For `NquadN=30, nSpec=30000, nRaman=100`, that's ~0.72 GB per buffer at Float64; a half-dozen buffers ≈ 4–5 GB. Modern 16+ GB GPUs handle this, but the numbers deserve a concrete check before the implementation lands. If OOM risk is material, add a chunked fallback that processes Δn-ranges in groups.

*Files touched:* `src/Inelastic/types.jl`, `src/Inelastic/inelastic_helper.jl`, `src/CoreRT/types.jl`, `src/CoreRT/tools/rt_helper_functions.jl`, `src/CoreRT/CoreKernel/doubling_inelastic.jl`, `src/CoreRT/CoreKernel/interaction_inelastic.jl`, `src/CoreRT/CoreKernel/elemental_inelastic.jl`.

*Acceptance:* (i) all tests pass within tolerance (see §8 on why this changes from bit-exact); (ii) Rayleigh = Cabannes + RRS energy conservation check continues to hold; (iii) CUDA profiler shows single `gemm_strided_batched` calls replacing per-Δn sequential loops; (iv) `@btime` on a representative forward-Raman run shows ~nRaman× speedup in doubling/interaction kernel time.

*Known risk:* reordering BLAS operations from sequential `for Δn` to a batched gather can change rounding. Mathematically the result is identical, but in floating-point the per-element accumulation order can differ. This is a real departure from the "bit-for-bit" claim in Christian's plan; we honestly declare a tolerance bar here (see §8).

## 8. Phase 4 — Reduce GPU synchronization

Drawn from Christian's Phase 4, minimal changes.

- 4a. Remove redundant `synchronize_if_gpu()` calls inside `batch_inv!` (`ext/gpu_batched_cuda.jl`). Keep one sync at the end of each doubling step, not after every LU call.
- 4b. Consolidate sync points in `doubling_inelastic.jl` so multiple layer operations can queue before the next barrier.
- 4c. Profile to identify any remaining sync bottlenecks; target <10 syncs per full RT run (down from ~2,000).

*Files touched:* `ext/gpu_batched_cuda.jl`, `src/CoreRT/CoreKernel/doubling_inelastic.jl`.

*Acceptance:* (i) all tests pass; (ii) profile confirms <10 syncs per RT run on GPU; (iii) measurable further wall-clock improvement on top of Phase 3.

## 9. Phase 5 — Documentation and examples overhaul

This phase is distinct from Christian's original and is expanded from the single migration-guide deliverable he planned. The goal is a release-quality documentation surface that a new user can onboard to in an evening, and that a returning contributor can navigate without re-reading the source.

Phase 5 runs *continuously from mid-Phase 1*. Each time a script is ported (1d), a new surface is integrated, or an API shift lands, the corresponding docs section is updated in the same PR. A final consolidation pass happens after Phase 4. Don't leave all the documentation to the end.

### 9a. Top-level README refresh

Replace the current `README.md` with a fast on-ramp:

- One-paragraph description: what vSmartMOM is, who uses it (Earth, planetary, astrophysics RT researchers), what it computes (polarized radiances + analytic Jacobians, elastic + Raman, CPU/GPU).
- Install block (`Pkg.add`, or `Pkg.develop` for the merged branch pre-release).
- 20-line runnable quick-start producing a plot. Something like: load a bundled YAML, build the model, run forward RT, plot.
- Feature matrix (table: forward ✅ / linearized ✅ elastic only / Raman ✅ forward only / GPU ✅ / Float32 ✅ / polarization ✅ / surfaces: Lambertian, Cox-Munk, canopy, RossLi, RPV, Rahman, fresnel, water refraction / aerosols: TOMAS15, two-moment, GEOS-Chem).
- Pointers to the user guide, example index, API reference, and migration guide.

### 9b. User guide (Documenter.jl pages in `docs/src/`)

Structure:

- **Getting started** — install, quickstart, first modifications.
- **Concepts** — what the Matrix Operator Method computes, what layers and the adding-doubling scheme mean, what polarization / quadrature settings control, what Raman physics is covered (RRS, VS, forward only), what analytic Jacobians are (elastic only).
- **Setting up a run** — walk through a YAML parameter file end to end, naming every field, pointing to the types it populates.
- **Choosing architecture** — CPU vs GPU, when to prefer each, Float32 vs Float64 trade-offs.
- **Surfaces** — one page per BRDF type with the math, references, and a working snippet: Lambertian, Cox-Munk, canopy, RossLi, RPV, Rahman, fresnel, water refraction.
- **Aerosols** — TOMAS15, two-moment, GEOS-Chem integration; refractive-index conventions.
- **Gas absorption** — HITRAN-driven cross-sections, line shapes (Voigt, Doppler, Lorentz), edition management.
- **Inelastic scattering (Raman)** — RRS vs VS, A-band, Cabannes decomposition; what's supported on GPU; **explicit note that linearized Raman is not supported**.
- **Solar model + SIF** — `F₀` plumbing, `SIF₀` plumbing, data files.
- **Linearized RT** — Sanghavi-Stephens / Sanghavi-Davis-Eldering methodology, parameter layout (`ParameterLayout`: aerosol sub-parameters, gas VMR, surface), output shape of `dR`/`dT`, elastic-only limitation stated up front.
- **Canopy and multisensor** — when to use them, what each pipeline computes.
- **Performance notes** — how to size batches, how to check allocations, CPU vs GPU tuning knobs.
- **Troubleshooting / FAQ.**

### 9c. `examples/` directory

Every example has three artifacts: the runnable `.jl` script, a companion `.md` walking through the science and the code, and a committed reference output (figure or small JSON of scalar outputs). Examples must be runnable end-to-end from a fresh clone — no hardcoded absolute paths, no data dependencies that aren't either bundled or documented.

Proposed starter set (each drawn from a sanghavi acceptance script or a unified feature):

- `examples/EMIT/emit_aerosol_forward.jl` — forward RT over the EMIT spectral range with aerosols + gas absorption (Raman-free variant). Adapted from `prototype_EMIT_aer_ht.jl` but with `RS_type = noRS()` since linearized Raman is not supported. A linearized-elastic variant lives next to it.
- `examples/O2_A_band/raman_sif_grid.jl` — forward Raman + SIF over a surface albedo × surface pressure × SZA grid. Adapted from `creategrid_O2Aband_RamanSIF.jl`, parameterized so truncated spectral ranges are easy to use for development.
- `examples/O2_B_band/raman_sif.jl` — analogous for O2 B band.
- `examples/O3_Huggins/polarized_raman.jl` — O3 Huggins band polarized Raman. From `O3Huggins_polRaman.jl`.
- `examples/OCO/oco_rayleigh.jl` — OCO-style Rayleigh-only forward simulation. From `creategrid_O2Aband_OCORayl.jl`.
- `examples/CarbonI/carbonI_prototype.jl` — CarbonI prototype. From `prototype_CarbonI.jl`.
- `examples/canopy/rami_canopy.jl` — RAMI benchmark canopy RT.
- `examples/ocean/coxmunk_ocean.jl` — Cox-Munk polarized ocean.
- `examples/aerosols/tomas_aerosols.jl` — TOMAS15 aerosol setup.
- `examples/linearization/simple_jacobian.jl` — linearized elastic RT with a minimal aerosol+albedo Jacobian, exercising the `ParameterLayout` ordering.

The sanghavi branch owner selects the authoritative list at execution time (see §10); this is a starting proposal.

### 9d. API reference

`docs/src/api/` autogenerated from docstrings via Documenter.jl. Every public export in `src/vSmartMOM.jl` must have a docstring covering: signature, purpose, argument types, return shape, a `# Examples` block that doctest-runs, and a physics reference where applicable. A pre-merge audit produces a checklist; each missing docstring is an explicit item in Phase 5 backlog.

### 9e. Developer docs

- `docs/dev_notes/migration_from_sanghavi.md` — from Phase 1c.
- `docs/dev_notes/sanghavi_optimization_inventory.md` — from Phase 0.
- `docs/dev_notes/architecture.md` — lift the `CLAUDE.md` content into human-readable form; keep `CLAUDE.md` as the AI-oriented quick index.
- `docs/dev_notes/extending_surfaces.md` — how to add a new BRDF.
- `docs/dev_notes/extending_scattering.md` — how to add a new scattering regime.
- `docs/dev_notes/gpu_notes.md` — the Raman GPU optimization audit and rationale for Phases 2–4 (lift from `raman_gpu_optimization.md` and `batched_kernel_benchmarks.md`).

### 9f. Inline comments in ported physics

For the sanghavi-branch author specifically, every block ported from sanghavi carries a comment of the form:

```julia
# SANGHAVI-PHYSICS: <short name of the physics operation>
# Source: sanghavi branch, <filename>:<linerange>, commit <abbrev SHA>
# See docs/dev_notes/migration_from_sanghavi.md#<anchor> for the pre-/post-refactor comparison.
<code>
```

These are greppable (`# SANGHAVI-PHYSICS:`) and let the sanghavi-branch author jump from a physics concept to its place in the merged code without reading the whole file.

### 9g. Tutorial pass

A single narrative tutorial at `docs/src/tutorials/first_retrieval.md` that walks from raw inputs to a plotted forward spectrum and a one-parameter analytic Jacobian, tying together the concepts pages. Should run in under 30 seconds on CPU to keep the tutorial loop fast.

**Acceptance for Phase 5:** (i) `make docs` (Documenter build) succeeds without warnings; (ii) every script in `examples/` runs from a fresh clone on a clean environment; (iii) every public export has a docstring; (iv) the top-level README's quickstart copy-pastes into a Julia REPL and produces the plot it claims.

## 10. Testing and verification strategy

### 10a. Tolerance bar

**Deliberately not pre-specified here.** Per-quantity tolerances (radiances I, Q, U, V; Jacobians dR/dp for each parameter class; Raman-scattered intensities ieR, ieT; SIF contributions) are decided by the sanghavi-branch author at execution time. Claude Code on the lab server must ask her, per acceptance script and per quantity, what rtol/atol bar constitutes equivalence. Baselines captured on the `sanghavi` branch are the reference.

Phases 2 is claimed as *bit-exact* in Christian's plan; keep that claim only if the Phase 0 reconciliation with sanghavi's `InteractionWorkspace` doesn't introduce floating-point-order changes. Phase 3 (4D → 3D with gather/batched) *will* change accumulation order in floating-point and should honestly declare a tolerance bar, not bit-exactness. Phase 4 should be bit-exact again since it only removes sync barriers.

### 10b. Regression harness

- `test/regression/runner.jl` — loads baselines, runs each acceptance script (full or spectrally-truncated variant), compares outputs via `isapprox` with per-quantity rtol/atol, records runtime via `@timed` or BenchmarkTools, records GPU allocation counts via `CUDA.@time`.
- `test/regression/baselines/` — reference outputs committed to the repo (JLD2 or NPZ, small enough to live in git; Git LFS only if total footprint exceeds ~50 MB).
- `test/regression/truncated_configs/` — YAML overrides that narrow spectral ranges for fast iteration.

### 10c. Long-running test exclusion

Introduce a `ENV["VSMARTMOM_SKIP_SLOW"]` flag (default `true` for interactive runs, `false` on CI). Tag slow tests with `@testset "slow begin ... end` and honor the flag. Wigner 3j symbol validation test (~60 s) is an existing instance.

### 10d. Gate per phase

| Phase | Test gate                                                                                   |
|-------|--------------------------------------------------------------------------------------------|
| 0     | Reconciliation note exists; no code change yet                                             |
| 1     | Full 474-test suite + ported sanghavi scripts + forward-Raman smoke tests all pass         |
| 2     | Phase 1 gate + bit-exact baseline match (unless Phase 0 reconciliation introduces changes) |
| 3     | Phase 1 gate + within-tolerance baseline match + Rayleigh = Cabannes + RRS energy check    |
| 4     | Phase 3 gate + bit-exact relative to end of Phase 3                                        |
| 5     | Phase 4 gate + `make docs` green + every `examples/*.jl` runs from clean clone             |

### 10e. Claude Code's role at execution

The plan will be handed to Claude Code on the lab server (where GitHub access and Julia environments are available). Claude Code **must** ask the sanghavi-branch author the following before starting work:

- **Acceptance script list.** Which of the scripts I've proposed in §9c (and any others from the sanghavi benchmarks folder) should be in the regression suite. Default proposal: `prototype_EMIT_aer_ht.jl` (with `RS_type = noRS()`), `creategrid_O2Aband_RamanSIF.jl`, plus a subset she picks from the full sanghavi benchmarks list.
- **Truncation strategy per script.** For each acceptance script, what spectral range / bands to narrow to for fast iteration.
- **Per-quantity tolerance bars.** rtol and atol for I, Q, U, V, ieR, ieT, dR, dT, SIF as relevant.
- **Coordination touchpoints with Christian.** Which phases should go through Christian for sign-off before merging into `unified-vsmartmom` (plausibly Phase 0 and Phase 4 at minimum).

## 11. Risks and open questions

- **Phase 0 reconciliation may take longer than one day** if sanghavi's `InteractionWorkspace` is architecturally incompatible with the rest of unified. Mitigation: box the worst-case outcome — if reconciliation stalls, ship the merge with the unified-side design and defer sanghavi's optimization to a follow-up; do not let Phase 0 block the branch creation indefinitely.
- **Phase 3 memory budget is not verified at the top-end problem size.** The `NquadN=30, nSpec=30000, nRaman=100` figures in `raman_large_scale.jl` may stress 16 GB GPUs once gather buffers are added. Plan the chunked-fallback path in Phase 3 even if not implemented initially.
- **`SIF₀` plumbing status is unclear.** Need to confirm on a live branch whether `RS_type.SIF₀` from sanghavi's `prototype_O2ABand_RRS_SIF.jl` has survived into unified's `RS_type` hierarchy. If not, Phase 1 has additional scope.
- **Documenter build may reveal unexpected docstring gaps.** Many public exports on unified currently lack docstrings. Phase 5 budget may grow once a concrete audit runs.
- **`examples/` script selection** is ultimately a sanghavi-branch-author call; the list in §9c is a proposal. Claude Code will ask at execution time.

## 12. Differences from Christian's 2026-03-21 plan

Summary for Codex and for Christian's review:

- **Out of scope: linearized Raman.** Christian's Phase 5 (deferred linearized Raman) is removed. The sanghavi-branch author confirms the EMIT / SIF / Raman workflows she cares about are all forward-only. The linearized path stays elastic-only permanently; this is a product decision, not a timing decision.
- **New: Phase 0 reconciliation with recent sanghavi optimization commits** (`854b44c`, `d75dacb`, `3e9926f`, `e89ec1c`, `083353b`). Prevents Phase 2 from re-inventing work already done.
- **New: Phase 5 is a full documentation and examples overhaul** (not just the migration guide from Christian's Phase 1c). Includes README rewrite, Documenter user guide, curated runnable examples with `.md` companions, API reference audit, developer docs, tutorial, and greppable `# SANGHAVI-PHYSICS:` inline markers.
- **Tolerance bars are deliberately deferred** to execution time, to be decided by the sanghavi-branch author per acceptance script and per quantity. Christian's "bit-for-bit identical" claim is retained only where defensible (Phase 2 and Phase 4, not Phase 3).
- **Hygiene cleanup** of sanghavi-side detritus is made explicit in Phase 1e. Christian's plan mentions `@show` removal but not the PNGs, `Manifest copy*.toml`, `.dat~` backups, and absolute paths.
- **Smoke tests for forward Raman** are made an explicit Phase 1 deliverable so Phases 2–4 have regression guards from day one.

## 13. Review questions for Codex

1. Does dropping linearized Raman as a permanent product decision raise any red flags you'd want articulated to the sanghavi-branch author before she commits? The reasoning is that the EMIT / SIF / Raman workflows she cares about are all forward-only, and the linearized path's elastic-only scope is already stable. Is there a science use case we should flag that this closes off?
2. Is the Phase 0 reconciliation scope right, or should it be broader (e.g., audit *all* 50 commits unique to sanghavi for overlap with unified, not just the five recent optimization commits)?
3. Phase 3 declares a tolerance bar rather than bit-exactness because reordering BLAS changes floating-point accumulation. Is this the right call, or is there a tighter ordering that preserves bit-exactness?
4. Phase 5 is ambitious — README + user guide + examples + API ref + developer docs + tutorial. Is it realistically sizeable alongside the technical phases, or should it be split into a Phase 5 (docs scaffolding) and a post-merge Phase 6 (examples and tutorials)?
5. Are there phase-ordering alternatives that would reduce overall risk — e.g., running Phase 5 documentation ahead of Phase 2–4 so the API surface freezes before performance work touches it?
6. Any missing phases, missing deliverables, or missing non-goals that a reviewer with vSmartMOM / GPU RT / Julia-scientific-package experience would flag?
