# Plan amendments — sanghavi-unified merge

**Date:** 2026-04-19
**Status:** Authoritative. Supersedes the current `plans/IMPLEMENTATION_PLAN.md` on points of conflict.
**Inputs this document consolidates:**
- `plans/CLAUDE_HANDOFF_BRIEF.md` (unchanged; canonical mission document)
- `plans/inventories/A_sanghavi_optimization_inventory.md`
- `plans/inventories/B_sanghavi_performance_baselines.md`
- `plans/inventories/C_sif_plumbing_status.md`
- `plans/inventories/D_physics_delta_summary.md`
- `plans/inventories/E_docs_state_audit.md`
- `plans/inventories/F_christian_plan_delta.md`
- User review session of 2026-04-19 (seven decisions below)

**What to do with this document:** Read it before revising `IMPLEMENTATION_PLAN.md`. Apply every directive under §1–§3. Use §4 as the amended phase-by-phase scaffold. Use §5 as the regression-check invariants. §6 lists the editorial passes needed on the inventories before they go to Christian.

---

## 1. Authority rule (binding, applies to every phase)

**Sanghavi is the authority for the inelastic path. Full stop.**

"Inelastic path" means, concretely:

- `src/Inelastic/` in full. All files, including nested sub-module structure. `InelasticScattering` is a sub-module (not a separate registered package), so this is a file copy, not a dependency bump.
- All `*_inelastic.jl` files in `src/CoreRT/CoreKernel/`: `interaction_inelastic.jl`, `elemental_inelastic.jl`, `elemental_inelastic_plus.jl`, `doubling_inelastic.jl`, `raman_atmo_prop.jl`.
- The inelastic-specific portions of `src/CoreRT/atmo_prof.jl`, `src/CoreRT/model_from_parameters.jl`, `src/CoreRT/types.jl`.
- The `RS_type <: noRS` branches and `InteractionWorkspace` allocation inside `src/CoreRT/rt_run.jl`. `rt_run.jl` itself is a merge, not a copy: unified's `RTModel` scaffolding plus sanghavi's inelastic branches.

This rule settles, by fiat, every per-item decision in Inventories A §5 and D §6 on inelastic files. The planning agent should not re-evaluate those decisions file-by-file; the answer for inelastic content is "sanghavi." Per-item rationale is retained in the inventories for historical reference only.

**Exceptions to the rule** (explicit carve-outs):

1. `reduce_profile` default flips — see §2.1.
2. `apply_D_matrix_elemental!` scalar shortcut from commit 9a26002 — dropped per §2.2.
3. `hem_R` / `hem_T` hemispheric API change from commit 9a26002 — deferred post-merge per §2.3.
4. `inelastic_helper_old.jl` on sanghavi is a 924-line frozen reference file; does not port.

---

## 2. Resolved open questions

Each decision below closes an open question that the inventories flagged. The planning agent should treat these as settled; do not re-raise them.

### 2.1. `reduce_profile` default — interpolation wins

**Decision:** Linear interpolation on uniform pressure half-levels (sanghavi's `reduce_profile`) becomes the default on `sanghavi-unified`. Bin-averaging (unified's current default) is preserved as `reduce_profile_binavg` (or equivalent) available via an opt-in keyword.

**Consequence:** Phase 1's physics-homogenization sub-phase broadens. It now covers two default switches, not one — Bodhaine Rayleigh and profile interpolation. The 474-test re-baselining workstream gets bigger. See §4 Phase 1a for the updated scope.

**TODO comment** to land at the `reduce_profile` call site during Phase 1:

```julia
# TODO: reduce_profile and getRayleighLayerOptProp both deserve refactoring
# beyond this merge. The current defaults are physics-forward (Bodhaine
# 1999 Eq. 30, interpolated profile) but the implementation layout mixes
# legacy and modern code paths via keyword args and name suffixes. Cleaner
# architecture: a single configurable ProfileReduction strategy (dispatch on
# type) and a single RayleighFormula strategy, both YAML-configurable.
# Tracked as a future PR; not a merge blocker.
```

### 2.2. `apply_D_matrix_elemental!` scalar shortcut — dropped

**Decision:** Do not port commit 9a26002's scalar shortcut that skips the D-matrix kernel launch when `n_stokes == 1`. The dispatch-correctness risk in a merge is not worth the small performance gain. If scalar-mode performance becomes a concern post-merge, the shortcut is a one-line follow-up with a dedicated test.

### 2.3. `hem_R` / `hem_T` API change — deferred

**Decision:** Commit 9a26002 bundles several unrelated changes: EMIT scripts, a Float32 YAML config, the D-matrix shortcut (dropped per §2.2), and the `hem_R`/`hem_T` hemispheric-integrated radiance output (breaking `rt_run` API change, 4-tuple → 6-tuple return). Port everything except `hem_R`/`hem_T`. The hemispheric API design is work-in-progress and will land as a separate post-merge PR. The merge does not need to decide between breaking change, keyword flag, or accessor approach.

Explicitly: `hem_R` / `hem_T` is a non-goal of this merge.

### 2.4. `α̅` `2π` factor — verified, adopt

**Decision:** Sanghavi's removal of the `2π` factor in the `α̅` frequency correction (commit 083353b) is **verified correct** by the user. Port via §1's authority rule. At the port site in `src/Inelastic/src/inelastic_cross_section.jl`, add a one-line comment citing the verification:

```julia
# Frequency correction (1 - (c·ν_eff/ω₀)²). ν_eff is wavenumber (cm⁻¹),
# ω₀ is stored in wavenumber units; no 2π factor. Verified correct against
# Buldakov et al. 1996 Eqs. 36a-39b by [author/source].
```

The planning agent should insert the actual verification source (notebook, paper section, commit message, email thread) during Phase 1.

Inventories A §8 Q1/Q2 and D §7 Q1/Q2 on this question are now closed.

### 2.5. SIF `0.5π / max` rescaling — intentional hack, document

**Decision:** The `0.5π / maximum(J_SIF)` rescaling in the SIF benchmark scripts is an intentional normalization that makes SIF magnitude data-independent for grid generation. It normalizes shape but discards absolute physical magnitude. Leave as-is for this merge.

**Comment to add** above the rescaling line in `creategrid_O2Aband_RamanSIF.jl` (and sibling scripts that use the same pattern):

```julia
# TODO: The (0.5π / maximum(J_SIF)) rescaling is an intentional hack to make
# SIF magnitude data-independent for grid generation. This normalizes shape
# but discards absolute physical magnitude. Revisit: confirm the downstream
# physics depends only on SIF shape (not absolute flux), or replace with
# physical units (mW/m²/cm⁻¹). Not a merge blocker.
```

**Consequence for the regression harness:** SIF correctness is a *pattern* check (ieR changes when SIF is nonzero), not an absolute-magnitude check. Revisit when the rescaling is addressed.

Inventory C §6 Q1 is closed.

### 2.6. Canopy + Raman — not coupled, document

**Decision:** Canopy + Raman coupling is a future project, not in scope for this merge. Unified's canopy code is noRS-only; sanghavi has no canopy. The combination `CanopySurface + RRS` is untested and has undefined behavior. Do not add a code guard; do not add a smoke test. Add a comment at the canopy dispatch site.

**Comment to land** during Phase 1 at the canopy dispatch site (find the right anchor during implementation — likely `src/CoreRT/CoreKernel/rt_kernel.jl` or `src/CoreRT/Surfaces/canopy_surface.jl`):

```julia
# TODO: Canopy + Raman (RRS/VS) coupling is not currently implemented. The
# canopy BRDF code is noRS-only; calling rt_run with CanopySurface + RRS has
# undefined behavior. Future work: couple canopy to the inelastic path.
```

Inventory D §7 Q5 is closed.

### 2.7. Aerosols module — user-facing WIP, wire in now

**Decision:** `src/Aerosols/Aerosols.jl` is user-facing but work in progress. It gets `include`d in `src/vSmartMOM.jl` so `using vSmartMOM.Aerosols` works. Header comment on the module documents WIP status. The orphan test `test_Aerosols.jl` wires into `runtests.jl`.

This moves out of Phase 7 (where the current plan has it) into Phase 1.

**Header comment** for `src/Aerosols/Aerosols.jl`:

```julia
# WIP: This module is user-facing but the API may evolve in follow-up PRs.
# Landed as part of the sanghavi-unified merge to close the
# `using vSmartMOM.Aerosols` export gap; further cleanup is a known
# follow-up workstream.
```

Inventory E §1.7 / §8 Q1 is closed.

---

## 3. Christian coordination — final shape

Per user's review: Christian agrees with the structure below. The planning agent should not add new Christian checkpoints.

- **Phase 0 checkpoint:** Christian signs off on the Inventory A + Inventory F delta and on this amendments document before branch `sanghavi-unified` is created. This is the main scope-alignment checkpoint.
- **Phase 4 checkpoint:** Christian reviews the `InteractionWorkspace`-as-adapted-to-`RTModel` design *before* kernel-wiring commits land. Framing for his review: "here's sanghavi's `InteractionWorkspace` ported to unified's architecture; these are the additional fields from your `InelasticWorkspace` proposal we've layered on; any more?" It's a review of additions, not a design review — the base design is sanghavi's by authority.
- **Phase 5 checkpoint (conditional):** If Phase 5 runs (i.e., Phase 4 measurements show headroom), Christian reviews the specific optimization targets and measurement protocol per optimization.
- **Final merge PR:** Does not require Christian's sign-off per user's §6 Q4 answer in the original brief. User drives.

---

## 4. Amended phase scaffold

The planning agent should rewrite `IMPLEMENTATION_PLAN.md` using this scaffold. Content below is directive for each phase; full phase-level wording (regression checks, exit criteria, etc.) comes from the existing plan plus the adjustments noted.

### Phase 0 — Reconciliation (no code)

**Additions beyond current plan:**

- Editorial pass on Inventories A, D, F, and 0_SUMMARY to reframe per-item inelastic decisions under §1's authority rule. Do this before Christian sees the inventory set. Details in §6.
- `σ_Rayl_coeff` storage convention confirmed via grep over `src/Inelastic/src/` for future-maintainer clarity. The authority rule resolves which code lands; this confirms which convention the code assumes, for documentation purposes.

**Exit criterion update:** "User approves this plan; Christian has acked the A/F delta, this amendments document, and the editorial inventory pass in writing."

### Phase 1 — Physics truth sync + Inelastic port

The current plan's Phase 1 is restructured into five sub-phases, in order.

#### Phase 1a — Physics default homogenization

**Scope:** Two default switches land together in one logical commit.

1. Migrate `getRayleighLayerOptProp` (or wherever the old `0.00864 × …` fit lives on unified) to Bodhaine 1999 Eq. 30.
2. Migrate `reduce_profile` default to linear interpolation on uniform pressure half-levels (sanghavi's current default). Preserve bin-averaging as `reduce_profile_binavg` available via opt-in keyword.
3. Grep the 474-test suite for reference values that depend on Rayleigh formula or layer averaging (most forward-RT tests, natraj, 6SV1, any multi-layer profile config).
4. Re-baseline affected tests. Review diff as a unit. Spot-check ~10 tests; verify shifts are consistent with combined Bodhaine + interpolation expectations; commit the batch.
5. Land the `reduce_profile` refactor TODO comment from §2.1 at the call site.

**Commit as one logical unit:** "Migrate Rayleigh formula and profile reduction defaults; re-baseline regression tests." Splitting the commits creates intermediate states where some tests pass on mixed-old-and-new defaults, making the diagnostic harder.

**Regression check:** 474 tests pass after re-baselining. Spot-check documentation records which tests were re-baselined and the reviewer's initials.

#### Phase 1b — Inelastic port (sanghavi authority)

**Scope:**

- Bulk-copy `src/Inelastic/` from sanghavi. Confirm the `include` path in `src/vSmartMOM.jl` resolves (nested `src/Inelastic/src/` structure is retained).
- Copy all `*_inelastic.jl` files in `src/CoreRT/CoreKernel/` from sanghavi.
- Adopt sanghavi's inelastic-specific portions of `src/CoreRT/atmo_prof.jl`, `src/CoreRT/model_from_parameters.jl`, `src/CoreRT/types.jl`.
- `src/CoreRT/rt_run.jl` is a merge, not a copy: unified's `RTModel` scaffolding, sanghavi's `RS_type` branches, sanghavi's `InteractionWorkspace` allocation site. Note: the actual `InteractionWorkspace` *landing and threading* is Phase 4, not Phase 1. Phase 1 leaves `_interaction_ws = nothing` at the allocation site and the kernel fallbacks to `similar()` per sanghavi's existing `workspace === nothing` branch.
- Port the rest of commit 9a26002's non-`_lin` content *minus* the D-matrix scalar shortcut (§2.2) and *minus* `hem_R`/`hem_T` (§2.3). EMIT scripts and the Float32 YAML config port normally.
- Land the `α̅` verification comment from §2.4.
- Land the canopy + Raman TODO comment from §2.6.
- Do not port `inelastic_helper_old.jl`.

**Regression check:** `test_forward_raman.jl` (currently orphan on unified — wire in during this sub-phase) matches sanghavi reference within user tolerance on a toy 1-band RRS config. Reference JLD2 committed to `test/reference/phase1b_inelastic_port.jld2`.

#### Phase 1c — `rt_run_ss` driver port

**Scope:**

- Copy the `rt_run_ss` driver body from sanghavi's `src/CoreRT/rt_run.jl` (~lines 258 and 450) into unified's `src/CoreRT/rt_run.jl`.
- Adapt `vSmartMOM_Model` field access to `RTModel` sub-struct access.
- Add smoke test `test_forward_ss.jl` and wire into `runtests.jl`.

**Regression check:** `test_forward_ss.jl` passes. The export leak (`rt_run_ss` exported on unified but undefined, per Inventory D §5) is closed.

#### Phase 1d — Aerosols module wire-in

**Scope:**

- Add `include("Aerosols/Aerosols.jl")` to `src/vSmartMOM.jl`.
- Add the WIP header comment from §2.7 to `src/Aerosols/Aerosols.jl`.
- Wire `test_Aerosols.jl` into `runtests.jl`.

**Regression check:** `using vSmartMOM.Aerosols; <one exported symbol>` resolves. `test_Aerosols.jl` runs and passes.

#### Phase 1e — `perturb_parameters.jl` port

**Scope:**

- Port `src/Testing/perturb_parameters.jl` from sanghavi.
- One smoke test.

This is carried forward from the current plan's Phase 3; it fits better here adjacent to the physics port.

**Regression check:** Smoke test passes.

**Phase 1 exit criterion:** All sub-phases complete; 474 tests pass post-re-baselining; orphan tests wired where scoped land; references committed.

### Phase 2 — Baseline capture + regression harness

**Amendments to current plan:**

- Sample count: N ≥ 5 for wall-clock measurements (was N = 3). Report median + MAD; trim outliers.
- `CUDA.math_mode!(:strict)` or equivalent for bit-exactness checks to avoid cuBLAS algorithm nondeterminism.
- Budget ~30 min per full-range RRS warm run in capture planning; reduced YAMLs target < 60 s per warm run.
- Primary gate metric: GPU used-memory delta. Secondary: CPU allocation count. Others (kernel launches, peak memory) diagnostic.

**Phase 2 splits into two steps:**

#### Phase 2a — Diagnostic baseline

Run the harness on unmodified `sanghavi-unified` (= post-Phase-1) vs sanghavi. Record the delta. This verifies Inventory D's physics-fix list was complete: if the post-Phase-1 delta exceeds what D predicted, D was incomplete and must be revisited before Phase 3 proceeds. This is a diagnostic step, not a gate on a number; it's a falsification test for Inventory D.

#### Phase 2b — Baseline capture

Once Phase 2a passes review, freeze the baseline JSONs under `test/benchmarks/baseline_output/sanghavi_9ee9a75/` and `test/benchmarks/baseline_output/sanghavi-unified_<sha>/`. These are the non-regression floor for Phase 4 and Phase 8.

**Regression check:** Physics correctness (sanghavi-unified matches sanghavi within user tolerance on I/Q/U/V/ieR/ieT for both acceptance scripts) passes.

### Phase 3 — SIF + benchmark script port

**Amendments:**

- `perturb_parameters.jl` port moves to Phase 1e (above). Remove from Phase 3.
- Before porting `creategrid_O2Aband_RamanSIF.jl`, land the `0.5π / max` TODO comment from §2.5.
- SIF regression test is a *pattern* check (ieR responds when SIF is nonzero), not an absolute-magnitude check, until the rescaling is revisited.

Rest of Phase 3 stands as currently written.

### Phase 4 — Workspace landing

**Amendments:**

- Design is sanghavi's `InteractionWorkspace` adapted to unified's `RTModel` architecture. Christian's `InelasticWorkspace` proposal contributes candidate additions (doubling-tmp fields, `batch_inv!` pivot/info sharing, shared 3D buffers with elastic `RTWorkspace`). Additions layer on top of sanghavi's design *if compatible*. If not compatible, sanghavi wins. This is the authority rule applied to workspace design.
- `staged::Bool = true` default. Phase 5 may eventually replace the memory win with batching, but Phase 5 is conditional and may exit as a no-op. Don't trade a certain 3.5 GB win for an uncertain one.
- Christian checkpoint (per §3) is "here's what's landing; any additions you want on top?" Not a design review.
- Phase 4 can run in parallel with a subset of Phase 3 — workspace-side work doesn't touch the SIF Lambertian injection or the benchmark-script ports. If resources permit, parallelize.

Rest of Phase 4 stands.

### Phase 5 — Further GPU headroom (conditional)

No structural change. Stands as currently written. Explicit exit-as-no-op path preserved.

### Phase 6 — Remaining scripts port + detritus cleanup

**Addition:** Before Phase 6 starts, produce `plans/phase6_script_port_list.md` — a concrete table drawn from Inventory B §1a listing each remaining sanghavi benchmark script with its runtime, data-file dependencies, and keep/port/drop decision. User signs off on the list before porting starts. Otherwise Phase 6 slides — "we thought 8, it's actually 23" is a real failure mode for research-code ports.

### Phase 7 — Docs overhaul

**Amendments:**

- Aerosols module `include` fix is no longer a Phase 7 item — it moved to Phase 1d (§2.7).
- Pre-flight step at the start of Phase 7: verify each `Tutorial_*.jl` under `docs/src/pages/tutorials/` runs against current `sanghavi-unified` tip before starting docstring backfill. Literate.jl tutorials that error through the `Base.getproperty` shim need fixing before `warnonly=false` flips.
- `warnonly=true` removal: dry-run `julia --project=docs docs/make.jl` *without* `warnonly=true` once at the start of Phase 7 and count actual errors. Scope the backfill work from the actual count rather than the current plan's ~30–50 estimate.
- CHANGELOG deepening covers HITRAN-2024, Float32, GPU Mie, Canopy, GEOSChem, batched kernels in addition to the `RTModel` break.

Rest of Phase 7 stands.

### Phase 8 — Final merge PR

**Additions:**

- Version bump in `Project.toml`.
- Git tag at merge commit.
- Release notes summarizing the merge.

Rest of Phase 8 stands.

---

## 5. Cross-phase regression invariants

The current plan's five invariants stand. Add one:

**6. Tolerance claims are honest.** Bit-exactness applies only to phases that preserve execution path verbatim (workspace threading with identical write order, sync removal). Phases that reorder accumulation (staging mode, 4D→3D batching) use tolerance-based verification with user-specified numbers.

---

## 6. Editorial passes on inventories before sending to Christian

The inventories were written before the user's review session locked the authority rule in §1. They still frame inelastic decisions as per-item evaluations. Before Christian reads them, apply these edits:

- **Inventory A §5**: Reframe "Adopt wholesale for items 1–10" as "these items come with the wholesale Inelastic port per the authority rule." §8 Q1 (staging) is resolved: `staged::Bool = true` default. §8 Q4 (Bodhaine default): answered in Phase 1a. §8 Q5 (reduce_profile default): answered in Phase 1a. §8 Q6 (delete `inelastic_helper_old.jl`): yes, by not porting.

- **Inventory D §6 (port list)**: Reframe P0/P1/P2 items as "these land via the Inelastic port per the authority rule." §7 Q1/Q2 on `α̅` unit convention: closed per §2.4. §7 Q5 on canopy+Raman: closed per §2.6.

- **Inventory F §5 ("Final delta summary")**: Update the Phase 2 recast wording. Current wording says "land sanghavi's InteractionWorkspace... extend to cover doubling temporaries and `batch_inv!` pivot/info where clear headroom exists." That's correct but softens under the authority rule — sanghavi's workspace design is the base by authority, not by reconciliation. Christian's doubling/pivot additions are candidates to layer on, not a synthesis.

- **Inventory 0_SUMMARY finding #3**: Currently reads "Build the merged workspace to Christian's layout + add `staged` as an opt-in flag." Under the authority rule, rewrite as "Build the merged workspace to sanghavi's `InteractionWorkspace` layout adapted to `RTModel`. Christian's proposed additional fields layer on if compatible. `staged::Bool = true` default."

- **Inventory 0_SUMMARY finding #4**: Phase 3 bit-exactness claim is already correctly flagged as false. No edit needed.

These passes are editorial, not factual. The evidence in each inventory stands.

---

## 7. Restatement of hard constraints (unchanged from handoff brief)

Repeated here for the planning agent's convenience; no changes from brief §3/§7:

- No linearized Raman, ever. Not a phase, not a deferred item, not a stretch goal.
- No touching the elastic path. `⊠` in `doubling.jl`, `interaction.jl`, `elemental.jl` stays as-is.
- No pre-specifying tolerances. Per-script, per-quantity, at execution time.
- Baselines come from sanghavi, not from unified.
- Preserve or improve sanghavi's existing GPU Raman optimizations; never roll back.
- Produce the plan for human review before any code change.

---

## 8. What to produce next

Write `plans/IMPLEMENTATION_PLAN_v2.md` incorporating §1–§7 above. Structure it as phases with explicit gates, files touched, acceptance criteria, and regression checks, matching the current plan's format. Cite this amendments document in the v2 plan's References section. Submit v2 to user for approval before any code lands on `sanghavi-unified`.
