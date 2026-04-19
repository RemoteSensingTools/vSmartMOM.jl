# Inventory F — Christian plan delta

**As-of date:** 2026-04-19
**Author of plan under review:** Christian Frankenberg's team, 2026-03-21
**Plan under review:** `/home/sanghavi/code/github/uni_vSmartMOM/docs/dev_notes/sanghavi_unified_merge_plan.md`
**Companion document:** `/home/sanghavi/code/github/uni_vSmartMOM/docs/dev_notes/raman_gpu_optimization.md`

**Authoritative scope document (supersedes Christian's plan):** `/home/sanghavi/code/github/uni_vSmartMOM/plans/CLAUDE_HANDOFF_BRIEF.md`

## Branch tips audited

| Worktree | Branch | Tip |
|---|---|---|
| `/home/sanghavi/code/github/uni_vSmartMOM/` | `unified-vsmartmom` | `a4e4187` ("Add batched-kernel and Raman scaling benchmarks with writeup") |
| `/home/sanghavi/code/github/vSmartMOM.jl/` | `sanghavi` | `8419745` visible in `git log` (handoff brief §4 confirms same tip) |

Both tips carry commits dated 2026-04-19, four weeks after Christian's plan was written. The staleness of the plan is therefore expected and the purpose of this inventory is to spell out exactly which references can and cannot be trusted, phase by phase.

No Julia was executed in producing this inventory; all checks are static.

---

## 1. Staleness audit of concrete references

Ten+ concrete references lifted from Christian's two documents, checked against the current worktrees. "Location in Christian's plan" uses shorthand `MP:§X` for `sanghavi_unified_merge_plan.md` and `RGO:§X` for `raman_gpu_optimization.md`.

| # | Reference | Claim in Christian's plan | Current-code status | Impact on plan |
|---|---|---|---|---|
| 1 | `interaction_inelastic.jl:293-305` — 13 `similar() + .= 0` at entry of `ScatteringInterface_11` | RGO:§1 / MP:GPU audit table: "13× `similar()` + `.=0` per call … ~8,300 GPU allocs/run" | **Holds on unified** (uni_vSmartMOM line 293–305 in `ScatteringInterface_11` RRS branch has the 13 `similar()` at those exact line numbers). **Does NOT hold on sanghavi** — sanghavi's `interaction_inelastic.jl` has `InteractionWorkspace` and the RRS `ScatteringInterface_11` handler (line 319+) takes `workspace` and only falls back to `similar()` when `workspace === nothing` (line 357–369). | Phase 2's framing ("eliminate 13 `similar()` from scratch") is already obsolete for sanghavi. Recast as "land sanghavi's InteractionWorkspace pattern in unified's architecture." |
| 2 | `interaction_inelastic.jl:70-86` — 2D `*` inside `for n₁ … for Δn …` | RGO:§4 / MP:GPU audit: "nSpec × nRaman sequential BLAS calls instead of one batched call" | **Holds on unified** at the cited line range (ScatteringInterface_01 branch, lines 70–86, still uses 2D `*` inside nested loops). Also holds on sanghavi (the `ScatteringInterface_01` family wasn't re-worked by sanghavi's InteractionWorkspace PR, which focused on `ScatteringInterface_11`). | Phase 3 batching target is still valid, but scoping note: sanghavi's workspace fix only touched `ScatteringInterface_11` RRS; `_01`/`_10` remain candidates on both branches. |
| 3 | `doubling_inelastic.jl:62-110` — `tmp3/tmp4/tmp5/tmp6` temporaries from `⊠` in nested `Δn` loops | RGO:§2 / MP:GPU audit table: "~3,600 allocs/run" | **Holds on both branches.** unified line 55–113 and sanghavi `doubling_inelastic.jl` line 55+ still allocate `tmp3/tmp4/tmp5/tmp6` via `⊠` in the `for Δn = 1:nRaman` loop. `similar()` count is 12 per file on both branches. | Doubling-side allocation work remains open on BOTH branches. Sanghavi's InteractionWorkspace does NOT cover doubling. New plan should explicitly note that Phase 2 on sanghavi-unified still has doubling headroom. |
| 4 | `ext/gpu_batched_cuda.jl:34,63-64` — `batch_inv!` allocates pivot/info | RGO:§3 / MP:GPU audit table: "~2,000 allocs/run" | **Partially holds on unified.** Line 32–39 of `uni_vSmartMOM/ext/gpu_batched_cuda.jl` still allocates pivot/info (via `CUDA.CUBLAS.getrf_strided_batched!` which returns new pivot) in the base `batch_inv!(X, A)` overload. A workspace-aware overload `batch_inv!(X, A, ws::RTWorkspace)` exists at line 47–55 and shares the same behaviour of re-allocating pivot/info (it just takes the workspace but does not actually read `ws.pivot`/`ws.info`). The Float32/Float64 pointer-based overloads at line 58–89 explicitly allocate `CuArray{Cint}` + `CUDA.zeros(Cint,...)` on every call. **On sanghavi, `ext/` does not exist at all** — the GPU extension is unified-only. | Phase 2e ("pass RTWorkspace to `batch_inv!` in inelastic path") still stands as a pending optimization; the overload plumbing is in place on unified but the workspace arrays are not actually consumed yet. Sanghavi does not carry this work. |
| 5 | `batch_inv!` has "2 sync barriers per call" (MP:GPU audit table, ~2000 stalls) | MP:Phase 4 / RGO:§6 | **Holds on unified.** Every `batch_inv!` overload in `uni_vSmartMOM/ext/gpu_batched_cuda.jl` calls `synchronize_if_gpu()` twice (after `getrf`, after `getri`); 12 `synchronize_if_gpu` references at lines 11, 21, 28, 35, 38, 51, 54, 67, 71, 84, 88. | Phase 4 sync-reduction target is still valid on unified. Sanghavi does not have this extension at all, so for the merged branch only the unified copy matters. |
| 6 | Proposed `InelasticWorkspace` struct layout with 23 fields (`gp_refl`, `tt_gp`, `J₁⁺/⁻`, `ieJ₁⁺/⁻`, `tmp_inv`, `T_inv`, `tmpJ₀⁺/⁻`, `tmpR⁻⁺/⁺⁻`, `tmpT⁺⁺/⁻⁻`, `tmpieJ₀⁺/⁻`, `tmpieR⁻⁺/⁺⁻`, `tmpieT⁺⁺/⁻⁻`, `buf3d_a/b/c`, `pivot`, `info`) | MP:§Phase 2a | **Stale.** Sanghavi already has a competing design in production: `InteractionWorkspace{A3,A4_mat,A4_src,A3_src}` at `vSmartMOM.jl/src/CoreRT/CoreKernel/interaction_inelastic.jl:16-36` with 15 fields covering the same responsibility PLUS a `staged::Bool` flag for per-direction CPU staging and three CPU-backed output arrays. It deliberately does NOT include doubling fields (`gp_refl`, `J₁⁺/⁻`, `ieJ₁⁺/⁻`) — those stay in `RTWorkspace` or remain as local `similar()`s in doubling. Christian's union of "doubling + interaction" temporaries in one struct is not what sanghavi adopted. | Phase 2a (invent a new struct) is replaced by Phase 2a (adopt sanghavi's struct, extend if doubling coverage is wanted). See Inventory A for the commit-by-commit rationale. |
| 7 | "Allocate once in `model_from_parameters()`" and "store it in the RT run context" (MP:§Phase 2b) | Intended allocation site | **Stale.** Sanghavi already allocates the workspace inside `rt_run.jl` — not `model_from_parameters.jl` — at `vSmartMOM.jl/src/CoreRT/rt_run.jl:122` (main SFI path) and `:379` (second code path). Allocation is conditional on `!(RS_type <: noRS)`. Two workspaces (`_interaction_ws`, `_interaction_ws2`) are allocated per run because the two code paths have different composite-layer shapes. | Phase 2b's allocation-site proposal is overridden by the existing pattern. New plan should adopt rt_run-side allocation; whether to migrate it into `model_from_parameters` is a separate call (potentially harmful because the composite-layer shape isn't known until rt_run walks the band loop). |
| 8 | "Files to modify" for Phase 1 includes `src/CoreRT/CoreKernel/rt_kernel_ss.jl` and `src/CoreRT/CoreKernel/interaction_ss.jl` as **new** files (ports from sanghavi) | MP:§Phase 1 Files to create/modify | **Stale — already landed.** Both files already exist on unified-vsmartmom: `uni_vSmartMOM/src/CoreRT/CoreKernel/rt_kernel_ss.jl` and `interaction_ss.jl`. They also exist on sanghavi. The "port single-scatter approximation kernel" task in Phase 1b is done (or at least file-scaffolded). | Phase 1b item on SS kernel port can be downgraded to "verify unified's SS kernel matches sanghavi's physics"; the file creation is not needed. |
| 9 | "Files to modify" for Phase 2 includes `src/CoreRT/tools/model_from_parameters.jl` | MP:§Phase 2 Files to modify | Path exists on unified (`uni_vSmartMOM/src/CoreRT/tools/model_from_parameters.jl`). On sanghavi the layout is flat (`vSmartMOM.jl/src/CoreRT/model_from_parameters.jl`) — there is no `tools/` subdirectory. File path is therefore unified-specific. | Low-severity correction; when merging, path on sanghavi-unified will be unified's. |
| 10 | "~19,000 GPU allocations per Raman run" | RGO:§"Estimated allocation counts per run" (total row) | **Stale and unverifiable.** The number was estimated (not measured) against unified's snapshot at an earlier date. Its sub-components (8,300 for interaction, 3,600 for doubling, 2,000 for batch_inv, 5,000+ for `⊠` temporaries) are each individually plausible from static inspection but the total was not captured by instrumentation. Sanghavi's InteractionWorkspace should eliminate the ~8,300 interaction count on the sanghavi side — sanghavi's current number is almost certainly much lower. | Treat "19k" as folkloric. Both branches need re-measurement. Cross-reference Inventory B; this is its domain. |
| 11 | "All tests pass — **bit-for-bit identical** results" after Phase 2 (MP:§Testing & Verification Strategy, item 7) | Bit-exactness claim | **Will hold** *only* if Phase 2 is restricted to workspace threading that preserves write order and accumulation order. Sanghavi's `staged=true` mode (which routes pass-1 outputs through CPU-backed arrays and reuses GPU 4D buffers across passes) introduces additional `Array(gpu) -> gpu_buf` copies but preserves the math. So the claim is defensible for pure workspace threading. | Keep bit-exactness as target for Phase 2 only in the non-staged path. Note that the CPU-staging mode needs its own regression check. |
| 12 | "All tests pass — **bit-for-bit identical** (same math, different memory layout)" after Phase 3 (MP:§Testing & Verification Strategy, item 10) | Bit-exactness claim for 4D→3D + gathered batched_mul | **False by construction.** `gemm_strided_batched` on a gathered 3D array performs its matmuls in a different accumulation order than nSpec × nRaman sequential 2D `mul!` calls. FP64 sum-reduction order matters; bit-exactness is not preserved. Polarized Stokes V in particular is sensitive to FP rounding drift (cf. commit `a5e0de5` which fixed a `mod`→`mod1` bug for exactly this reason). | Flag for Christian: Phase 3 needs tolerance-based verification, not bit-exact. The actual tolerance is out of scope here (user-driven per handoff brief §6, §7). |
| 13 | `RamanIndexMap` struct proposal with `nTotal::Int`, `el_n₁::IT`, `el_n₀::IT`, `ranges::Vector{UnitRange{Int}}` in `src/Inelastic/types.jl` | MP:§Phase 3a | Not yet present on either branch. `i_λ₁λ₀` and `get_n₀_n₁` exist on both (see grep results: `get_n₀_n₁` in `src/Inelastic/inelastic_helper.jl`, `InelasticScattering.jl`, `stellar_inelastic_helper.jl`, `src/CoreRT/CoreKernel/interaction_inelastic.jl`, `doubling_inelastic.jl`, `interaction_multisensor.jl`, `interlayer_flux.jl` on both worktrees). sanghavi commit `854b44c` optimized `get_n₀_n₁` but did not introduce a flattened map. | Phase 3's index-map design is still unbuilt on both branches; the proposal remains valid as a forward design. Note that sanghavi's `get_n₀_n₁` optimization should be inspected to ensure the flat-map design subsumes or improves on it (see Inventory A). |
| 14 | "`AddedLayerRS` and `CompositeLayerRS`: `ier⁻⁺: (nQ, nQ, nSpec, nRaman) → (nQ, nQ, nTotal)`" | MP:§Phase 3b | Structs currently 4D on both branches — `uni_vSmartMOM/src/CoreRT/types.jl:215` (`CompositeLayerRS`) and `:245` (`AddedLayerRS`). Flattening not begun. | Forward design still valid; no stale references. |
| 15 | Phase 5: "Linearized Raman (Deferred — separate PR)" | MP:§Phase 5 | **Scope-removed, not deferred.** The handoff brief §3 (Non-goals) is explicit: "Linearized Raman. Raman + analytic Jacobians together is permanently out of scope. The linearized RT path stays `RS_type = noRS()` only. No `AddedLayerRS`-style linearized constructors, no Raman branch in `rt_kernel_lin.jl`. This is a product decision by the sanghavi-branch author, not a deferral — it will not land in a follow-up PR either." | Phase 5 must be removed entirely — not re-scheduled. See §3 below. |
| 16 | "remove `@show` debug statements from interaction kernels" | MP:§Phase 1e Cleanup | Still open on unified — `@show` appears in `interaction_inelastic.jl` (e.g., `@show "interaction 00"` at line 80 of sanghavi copy, similar statements on unified's copy). `#bla` and the commented `repeat(r⁻⁺,1,1,1,nRaman) ⊠ reshape(ieJ₁⁻,...)` survive in `uni_vSmartMOM/src/CoreRT/CoreKernel/doubling_inelastic.jl` at lines 55 and 86. | Cleanup item still valid, pointer targets still accurate. |
| 17 | Removed types/files list: "`vSmartMOM_Model`, `rt_run_bck.jl`, `types_inelastic.jl`" (MP:§Phase 1c migration guide bullet) | Files removed on unified | `rt_run_bck.jl`: **does not exist on unified**, **still exists on sanghavi** (`vSmartMOM.jl/src/CoreRT/rt_run_bck.jl`). `types_inelastic.jl`: **does not exist on unified**, **still exists on sanghavi** (`vSmartMOM.jl/src/CoreRT/types_inelastic.jl`). `vSmartMOM_Model`: removed on unified in favour of `RTModel{ARCH,FT}` per `uni_vSmartMOM/src/CoreRT/types.jl:843`. | Migration-guide bullet is correct. Cleanup note: when merging sanghavi → sanghavi-unified, `rt_run_bck.jl` and `types_inelastic.jl` should be dropped (they are dead code even on sanghavi — `rt_run_bck.jl` is a backup of rt_run.jl and unreferenced). |

---

## 2. Phase-by-phase reconciliation

### Phase 1 — Branch Creation & Sanghavi Compatibility

- **Christian's intent:** Create `sanghavi-unified` branch, port single-scatter (SS) kernel and EMIT/Balsamic scripts from sanghavi, write migration guide, clean up debug residue.
- **What has since landed / shifted:**
  - SS kernel files already exist on both branches (`rt_kernel_ss.jl`, `interaction_ss.jl`). No port needed — physics parity audit only.
  - EMIT/Balsamic/RamanSIF scripts list is much richer than Christian's plan assumed; see `test/benchmarks/` listings in Inventory B §1a — there are ~25 Raman-aware scripts on sanghavi, not two.
  - The sanghavi top level has research detritus (PNGs, `Manifest copy*.toml`, `.dat~` files, `nohup.out`, etc. — see the top-level listing) that needs to be excluded from the merge rather than carried.
  - `RTModel` migration docs exist at `uni_vSmartMOM/CLAUDE.md`, `uni_vSmartMOM/docs/dev_notes/RAMAN_CODE_HANDOFF.md`; Christian's "write a migration guide" task is partially satisfied.
- **Verdict:** Phase re-scopes but does not go away. Rename to "Phase 1 — Branch creation, port remaining sanghavi content (benchmarks, SIF plumbing, single scatter parity), migration guide."

### Phase 2 — Eliminate allocations in inelastic path

- **Christian's intent:** Invent a new `InelasticWorkspace` (doubling + interaction temporaries merged), allocate in `model_from_parameters`, refactor `doubling_inelastic.jl` and `interaction_inelastic.jl` to consume it. Claim bit-exact results.
- **What has since landed / shifted:**
  - Sanghavi has `InteractionWorkspace` (interaction side only) already threaded through `rt_run.jl` (both SFI code paths). Covers `ScatteringInterface_11` RRS only.
  - Sanghavi's struct design differs from Christian's proposal: 15 fields + per-pass CPU staging, NOT 23 fields with doubling merged in.
  - Doubling-side allocations (`tmp3`…`tmp6`, `gp_refl`, `J₁⁺/⁻`, `ieJ₁⁺/⁻`) are **still present on sanghavi** — sanghavi's optimization did not extend to doubling.
  - `batch_inv!` pivot/info allocations are unified-only (no `ext/` on sanghavi); the workspace-aware overload exists but is not actually wired (pivot/info from the workspace aren't consumed).
- **Verdict:** Phase 2 is **recast** per the handoff brief: "land sanghavi's InteractionWorkspace in unified's architecture as the starting point; extend to cover doubling temporaries and `batch_inv!` pivot/info where clear headroom exists." This is NOT "eliminate from scratch." The handoff brief §7 is explicit that reverting sanghavi's work is unacceptable.

### Phase 3 — Flatten 4D → 3D with index map

- **Christian's intent:** Build `RamanIndexMap`, reshape `AddedLayerRS` / `CompositeLayerRS` arrays, replace nRaman-sequential kernels with single gathered `batched_mul`.
- **What has since landed / shifted:**
  - Neither branch has attempted the flattening. sanghavi's `get_n₀_n₁` optimization (commit `854b44c`) is the closest relative.
  - `test/benchmarks/raman_batched_ops_benchmark.jl` on sanghavi explicitly benchmarks the "should we batch this?" question at NQ=9, NSPEC=5500, NRAMAN=963, NL=5000, comparing `⊠` broadcast vs `batched_mul!` vs 2D `mul!` vs `batch_inv!` and a full doubling-like simulation. That benchmark is the data Phase 3 should reference.
- **Verdict:** Phase still makes sense. **But the bit-exactness claim is false by construction** (see audit item #12). Recommend tolerance-based verification; handoff brief §6 leaves the tolerance numbers for the user to specify at execution time, so this inventory does not propose numbers.

### Phase 4 — Reduce GPU synchronization

- **Christian's intent:** Remove redundant `synchronize_if_gpu()` calls from `batch_inv!` and consolidate doubling syncs.
- **What has since landed / shifted:**
  - `ext/gpu_batched_cuda.jl` still has 12 `synchronize_if_gpu` calls across the various `batch_inv!` overloads (two per overload, as Christian noted). No consolidation has been done.
  - `doubling_inelastic.jl` still has the pre-existing sync pattern on both branches.
- **Verdict:** Phase still makes sense as-scoped. Targets and file pointers are current.

### Phase 5 — Linearized Raman

- **Christian's intent:** Deferred to a later PR; design notes say "use flat 3D layout from the start."
- **What the handoff brief mandates:** **Removed from scope.** Permanent. Not a deferral. Handoff brief §3 (Non-goals) first bullet:
  > "Linearized Raman. Raman + analytic Jacobians together is permanently out of scope. The linearized RT path stays `RS_type = noRS()` only. No `AddedLayerRS`-style linearized constructors, no Raman branch in `rt_kernel_lin.jl`. This is a product decision by the sanghavi-branch author, not a deferral — it will not land in a follow-up PR either."
- **Verdict:** **Delete Phase 5.** Do not renumber it to a future phase. Do not label as "deferred." The user has explicitly directed this be treated as dead-lettered.

---

## 3. "Files to modify" verification

| Christian's file path (in plan) | Phase | Exists on unified? | Exists on sanghavi? | Notes |
|---|---|:---:|:---:|---|
| `src/CoreRT/CoreKernel/rt_kernel_ss.jl` | 1 (new/port) | ✅ | ✅ | Already present; treat as physics-parity audit target, not port target. |
| `src/CoreRT/CoreKernel/interaction_ss.jl` | 1 (new/port) | ✅ | ✅ | Already present; same as above. |
| `examples/EMIT/`, `examples/Balsamic/` | 1 (new) | ❌ | ❌ (content lives under `test/benchmarks/` on sanghavi) | Create under `examples/` on sanghavi-unified; port from `test/benchmarks/` on sanghavi. |
| `docs/dev_notes/migration_from_sanghavi.md` | 1 (new) | ❌ | ❌ | Related material scattered across `uni_vSmartMOM/CLAUDE.md` and `uni_vSmartMOM/docs/dev_notes/RAMAN_CODE_HANDOFF.md`; new guide should consolidate. |
| `src/CoreRT/types.jl` | 2 (add `InelasticWorkspace`) | ✅ | ✅ | On unified: `RTModel` at line 843, `CompositeLayerRS` at 215, `AddedLayerRS` at 245. No InelasticWorkspace yet. |
| `src/CoreRT/tools/model_from_parameters.jl` | 2 (allocate workspace) | ✅ | ❌ (sanghavi has flat `src/CoreRT/model_from_parameters.jl`) | Path mismatch; on sanghavi-unified the unified path wins. Allocation site is a design call (see §2 Phase 2 verdict). |
| `src/CoreRT/CoreKernel/doubling_inelastic.jl` | 2 (use workspace) | ✅ | ✅ | Both files: 12 `similar()` calls; `tmp3/tmp4/tmp5/tmp6` allocations in nRaman loops intact. Doubling is the biggest remaining headroom post sanghavi's InteractionWorkspace. |
| `src/CoreRT/CoreKernel/interaction_inelastic.jl` | 2 (use workspace) | ✅ | ✅ | On unified: legacy 13 `similar()` at lines 293–305 and 440–452 (`ScatteringInterface_11` RRS + VS variants). On sanghavi: `InteractionWorkspace` struct at 16–36, used at 319–369 (RRS) and 536–548 (VS) with an explicit `workspace === nothing` fallback. |
| `src/CoreRT/rt_run.jl` | 2 (thread workspace) | ✅ | ✅ | On sanghavi: `InteractionWorkspace` created at lines 122 (`_interaction_ws`) and 379 (`_interaction_ws2`). On unified: no workspace plumbing yet. |
| `src/Inelastic/types.jl` | 3 (add `RamanIndexMap`) | ✅ | ✅ | No RamanIndexMap on either branch. Path is consistent. |
| `src/Inelastic/inelastic_helper.jl` | 3 (build index map) | ✅ | ✅ | Both carry `get_n₀_n₁`. |
| `src/CoreRT/tools/rt_helper_functions.jl` | 3 (add `gather!`, update `default_matrix_ie()`) | ✅ | ❌ (sanghavi: `src/CoreRT/rt_helper_functions.jl`) | Same path mismatch as `model_from_parameters`. |
| `src/CoreRT/CoreKernel/elemental_inelastic.jl` | 3 (write to flat layout) | ✅ | ✅ | Present on both. |
| `ext/gpu_batched_cuda.jl` | 4 (remove inner syncs) | ✅ | ❌ (no `ext/` on sanghavi; the CUDA extension is unified-only) | Sanghavi has `src/CoreRT/gpu_batched.jl` instead — different layout, different dispatch. The extension-driven path is unified-only. |

No files referenced by Christian's plan have been renamed or deleted on unified. Path-level staleness is limited to (a) sanghavi's flat `src/CoreRT/` layout vs. unified's `src/CoreRT/tools/` split, which the merge resolves by adopting unified's layout, and (b) the `ext/` directory being unified-only. Christian's references are usable once rebased against unified's file layout.

---

## 4. Allocation-count re-measurement note

The "~19,000 GPU allocations per Raman run" figure in `raman_gpu_optimization.md` §"Estimated allocation counts per run" is an **estimate**, not a measurement. Its sub-components add to 18,900:

- `interaction_helper!` `similar()`: ~8,300
- `doubling tmp3/4/5/6` in Δn loop: ~3,600
- `batch_inv!` pivot/info: ~2,000
- `⊠` temporaries in expressions: ~5,000+

On the current branches this figure is neither verified nor representative. In particular:

- On **sanghavi**, the `InteractionWorkspace` PR (commit `854b44c` + follow-ups) should have eliminated most of the 8,300 interaction-side allocations for the RRS `ScatteringInterface_11` path. It does NOT cover `_01`/`_10` interfaces nor doubling, so those lines in the budget remain.
- On **unified**, the 8,300 number may still be approximately current because InteractionWorkspace was never merged.
- `batch_inv!` pivot/info allocations are **unified-only** (sanghavi has no `ext/`). The 2,000 line item does not apply to sanghavi's CPU-first stack.
- The `⊠`-temporary line of 5,000+ is mostly inherited from the same code paths that the above addresses; independent measurement is the only way to separate them.

**Recommendation (for Inventory B to own the execution):**

1. Re-measure both branches under matched workload. Canonical target is the `O2_parameters2_1band_opttest.yaml` workload used by `raman_optimization_baseline.jl` on sanghavi (see Inventory B §1a). The `timing_baseline.txt` in `raman_opttest_output/` is stale (it records `RRS = 1671.851 s`, `noRS = 2.836 s` from a previous run on a different tip of sanghavi) and explicitly flagged as non-authoritative in Inventory B.
2. Use `CUDA.@time` or `NVTX.@range` for allocation counts; static inspection is not sufficient.
3. Publish the measured floor into Inventory B and reference it from the implementation plan; do NOT carry "19k" forward in Christian's form.

---

## 5. Final delta summary — what Christian needs to know

### Scope changes made by the sanghavi-branch owner (binding)

- **Linearized Raman is permanently out of scope.** Not deferred. Phase 5 of Christian's plan is struck, not re-scheduled. Reason: user product decision documented in the handoff brief §3. No `AddedLayerRS`-side linearized constructors, no Raman branch in `rt_kernel_lin.jl`.
- **The elastic path is not to be touched.** `⊠` in `doubling.jl`, `interaction.jl`, `elemental.jl` stays as-is. Inelastic path only. Handoff brief §7.
- **Baselines come from sanghavi, not from unified.** The regression harness compares against numbers captured on the sanghavi branch for each acceptance script. Handoff brief §7. Unified's pre-port numbers are not the target.
- **Sanghavi's existing GPU Raman optimizations must be preserved or improved, never reverted.** Specifically: `InteractionWorkspace` struct and its rt_run-level threading, `get_n₀_n₁` optimization, per-direction CPU-staging, the `raman_batched_ops_benchmark.jl` benchmark infrastructure. Any plan that reimplements against Christian's older `InelasticWorkspace` proposal without a concrete improvement over sanghavi's design is a regression.
- **Phase 3's bit-exactness claim is dropped.** 4D→3D + gathered batched_mul changes BLAS accumulation order. Tolerance-based verification only, per a user-specified tolerance at execution time. Tolerances are not pre-specified in any inventory or plan.

### Staleness findings by severity

- **High severity** (material to the plan's content or conclusions):
  - Audit items **#1, #6, #7, #11, #12, #15** — `InelasticWorkspace` design is superseded by sanghavi's `InteractionWorkspace`; allocation site is `rt_run.jl` not `model_from_parameters.jl`; Phase 5 is not a deferral; Phase 3 bit-exactness is false.
- **Medium severity** (targets still valid, ancillary claims need correction):
  - Audit items **#3, #4, #5, #10** — doubling allocations untouched on both branches, `batch_inv!` workspace overload exists but is not plumbed, sync-reduction targets still accurate, "19k" is estimate-not-measurement.
- **Low severity** (editorial / path corrections):
  - Audit items **#8, #9, #13, #14, #16, #17** — SS kernel already ported, unified path differs from sanghavi path, flattening work not yet started, debug residue still present, removed-files list accurate.

### Phase-by-phase recommendation

- **Phase 1** — Proceed as scoped, but scope the SS-kernel work down to a physics-parity audit (files already exist) and expand the benchmarks/SIF port to match the real content set (Inventory B §1a).
- **Phase 2** — **Recast**: "land sanghavi's `InteractionWorkspace` in unified's architecture; extend coverage where clear headroom remains (doubling temporaries, `batch_inv!` pivot/info, `ScatteringInterface_01`/`_10` 2D-mul loops)." Drop Christian's proposed struct layout in favour of sanghavi's existing one. Bit-exactness target applies only to the non-staged path.
- **Phase 3** — Proceed as scoped. **Replace bit-exactness claim with tolerance-based verification** (tolerance numbers to come from the user at execution time, per handoff brief §6 and §7).
- **Phase 4** — Proceed as scoped. Targets in `ext/gpu_batched_cuda.jl` and `doubling_inelastic.jl` are current.
- **Phase 5** — **Remove.** Not "defer." The linearized Raman path is not to be built.

---

## 6. Open questions

1. **Christian-side sign-off checkpoints.** Handoff brief §6 item 4 asks the user which phases go through Christian for review. At minimum this inventory (Inventory F) is a review artifact; the full implementation plan is another. Does Christian also want to sign off on the Phase 2 struct design before code lands? Open for the user to decide.
2. **Scope of doubling-side work.** Sanghavi's `InteractionWorkspace` does not cover doubling. The handoff brief does not prohibit extending it. Does Christian's team have a preference between (a) a separate `DoublingInelasticWorkspace`, (b) a unified `InelasticWorkspace` that subsumes interaction + doubling, or (c) keeping doubling temporaries local pending Phase 3 (because Phase 3's 4D→3D reshape will invalidate doubling's temporary shapes anyway)? Not decided here.
3. **Per-pass CPU staging.** Sanghavi's `staged::Bool` mode trades GPU memory for CPU↔GPU bandwidth. Is that trade acceptable to Christian's target hardware (noting that the 16+ GB GPUs assumed by Christian's plan Phase 3 memory trade-off may or may not match sanghavi's deployment targets)? Open.
4. **Does Phase 4 (sync reduction) go in before or after Phase 3?** Christian ordered 2→3→4 strictly. Phase 4 is low-risk and independent of Phase 3. It could run in parallel with Phase 3 development. Not pre-decided here.
5. **Disposition of dead code on sanghavi-unified.** `rt_run_bck.jl`, `types_inelastic.jl`, `inelastic_helper_old.jl` exist only on sanghavi and are unreferenced. Drop on merge (default) or preserve in a history tag for archaeology? Open question for the user / Christian's team.
