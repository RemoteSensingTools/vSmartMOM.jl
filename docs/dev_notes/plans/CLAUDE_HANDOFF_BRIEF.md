# Claude handoff brief — `sanghavi-unified` merge

*For a Claude Opus Plan agent running on the lab server. This is a framing document, not an implementation plan. Your job is to generate the implementation plan after inventorying the two branches, asking the user the open questions below, and respecting the guardrails at the bottom. Do not make code changes until the plan you produce has been reviewed by the user.*

---

## 1. Mission

vSmartMOM.jl has two active branches that need to be merged into a new branch `sanghavi-unified` that eventually replaces `main`. The merge must (a) preserve all forward-path physics that runs on the sanghavi branch today, within user-defined numerical tolerance, (b) preserve or improve the GPU Raman compute and memory performance sanghavi has already achieved, (c) live inside the unified-vsmartmom architecture (`RTModel{ARCH,FT}` hierarchy, Documenter-driven docs, 474-test suite, Cox-Munk/canopy/RossLi/RPV surfaces, Aerosols module, IO module, Float32 support, HITRAN version tracking), and (d) ship with release-quality documentation and runnable examples so first-time users can onboard without asking the maintainers.

Success is measured script-by-script: a list of sanghavi acceptance scripts (EMIT, O2 A/B band Raman+SIF grids, OCO, CarbonI, O3 Huggins, and others the user selects) runs on the merged branch and reproduces the sanghavi-branch baselines within tolerances the user specifies at execution time.

## 2. Scope (in)

- Branch creation (`sanghavi-unified` from unified-vsmartmom) and physics/infra port from sanghavi.
- Preservation or improvement of sanghavi's existing GPU Raman optimizations (see §4 for what already exists on sanghavi).
- Further GPU Raman optimization beyond what sanghavi has today, if an inventory shows clear headroom (allocation elimination, 4D→3D flattening with index-map-based gather, sync reduction).
- Migration guide (old-sanghavi API → new-unified API) with every converted script as a worked pair.
- Full documentation overhaul: README rewrite, Documenter user guide, curated `examples/` with runnable `.jl` + narrative `.md` + committed reference output, API reference audit, developer docs, a narrative tutorial.
- Regression harness capturing numerical outputs, runtime, and allocation counts for each acceptance script. Baselines captured on sanghavi before any refactoring touches the new branch.
- Hygiene cleanup of sanghavi-side research detritus (top-level PNGs, `Manifest copy*.toml`, `.dat~` backups, hardcoded absolute paths in scripts, `rt_run_bck.jl`, `inelastic_helper_old.jl`, in-source prototype/plot folders).

## 3. Non-goals (out of scope — do not plan these)

- **Linearized Raman.** Raman + analytic Jacobians together is permanently out of scope. The linearized RT path stays `RS_type = noRS()` only. No `AddedLayerRS`-style linearized constructors, no Raman branch in `rt_kernel_lin.jl`. This is a product decision by the sanghavi-branch author, not a deferral — it will not land in a follow-up PR either.
- **Regressing sanghavi's existing GPU Raman optimizations.** Repeat: any plan that strips out sanghavi's `InteractionWorkspace`, `get_n₀_n₁` optimization, per-direction CPU-staging, or Raman benchmark infrastructure in favor of reimplementing against unified's older design is unacceptable. Preserve or improve; do not roll back. See §4 for specifics.
- **Elastic-path performance retuning.** The elastic `⊠` notation in `doubling.jl`, `interaction.jl`, `elemental.jl` works well with Julia's JIT and is valued for readability. Touch only the inelastic path.
- **New physics.** No new scattering regimes, no new surface types, no new retrieval operators. This merge consolidates; it does not extend.
- **Mass reformat / lint sweep.** Code changes must be tied to concrete merge or performance objectives.
- **Pre-specifying tolerance bars.** Tolerances are user-facing choices per-script and per-quantity. Ask at execution time (see §6).

## 4. Verified repo facts

The lab server has both branches checked out as linked worktrees:

- **sanghavi worktree:** `/home/sanghavi/code/github/vSmartMOM.jl/` — branch `sanghavi`, tip `8419745` (2026-04-19, "Ignore .jld2/.jld files to keep large data out of the repo"). Owned by the sanghavi-branch author.
- **unified-vsmartmom worktree:** `/home/sanghavi/code/github/uni_vSmartMOM/` — branch `unified-vsmartmom`, tip `a4e4187` (2026-04-19, "Add batched-kernel and Raman scaling benchmarks with writeup"). Maintained by Christian Frankenberg's team.

Both branches have commits dated today. Both are actively developed. You must read both when inventorying; assumptions drawn from only one worktree are likely wrong.

### 4.1 Already on unified-vsmartmom (do not re-do)

- Hierarchical `RTModel{ARCH,FT}` (Oceananigans-style) with `solver`, `geometry`, `quad_points`, `atmosphere`, `optics`, `surfaces` sub-structs. Flat `vSmartMOM_Model` is removed with a `Base.getproperty` back-compat shim.
- Linearized RT (elastic only): `rt_run(model, lin_model, NAer, NGas, NSurf)` returns `(R, T, dR, dT)`. `ParameterLayout` struct for Jacobian index arithmetic. See `uni_vSmartMOM/src/CoreRT/rt_run_lin.jl`.
- Forward Raman physics ported from sanghavi in commit `59f8de8` ("Port sanghavi branch physics fixes to forward kernels").
- Cox-Munk polarized ocean, canopy, RossLi, RPV, Rahman, fresnel, water refraction surfaces.
- Aerosols module with TOMAS15, two-moment, GEOS-Chem integration, refractive index handling.
- GPU-accelerated NAI2 Mie scattering via KernelAbstractions (commit `d6bd45c`).
- HITRAN direct download with version tracking (commit `df19b06`).
- Float32/Float64 flexibility.
- 474-test suite at `uni_vSmartMOM/test/runtests.jl`.
- Developer docs under `uni_vSmartMOM/docs/dev_notes/`: `sanghavi_unified_merge_plan.md` (Christian's 2026-03-21 plan — *historical, superseded on scope AND stale on code references; see §4.5*), `raman_gpu_optimization.md` (GPU allocation audit — ~19k allocs/run on unified's Raman path, same staleness caveat), `RAMAN_CODE_HANDOFF.md` (addressed to the sanghavi-branch author), `LINEARIZATION_BUGS.md` (bugs found porting sanghavi's Jacobian code), `batched_kernel_benchmarks.md` (2026-04-19 writeup; current).
- GPU batched extension at `uni_vSmartMOM/ext/gpu_batched_cuda.jl`; elastic path uses pre-allocated `RTWorkspace`; inelastic path does not currently use the workspace-aware `batch_inv!`.
- `CLAUDE.md` at `uni_vSmartMOM/CLAUDE.md` — architecture and terminology reference for AI assistants.

### 4.2 Already on sanghavi — MUST PRESERVE OR IMPROVE

The sanghavi branch has the following optimization work already landed. Any plan must account for this when scoping Phase-2-style work; any plan that ignores or regresses it is wrong.

- **`InteractionWorkspace` struct and threading.** Sanghavi has already introduced a Raman workspace abstraction and threads it through the forward Raman path. Files to inventory:
  - `vSmartMOM.jl/src/CoreRT/CoreKernel/interaction_inelastic.jl` — sanghavi-side `InteractionWorkspace` implementation and staged Raman memory optimization.
  - `vSmartMOM.jl/src/CoreRT/rt_run.jl` — sanghavi forward Raman path already threads the workspace.
- **`get_n₀_n₁` optimization** (commit `854b44c`, "Add Float32 support, InteractionWorkspace, and optimize get_n₀_n₁").
- **Per-direction CPU-staging for interaction workspace** (commit `d75dacb`).
- **Raman optimization test infrastructure** (commit `3e9926f`).
- **Batched-ops benchmark** (commit `e89ec1c`; related file `vSmartMOM.jl/test/benchmarks/raman_batched_ops_benchmark.jl`).
- **Baseline script updates** for what the Raman optimization work is measured against.

Any plan must start by *inventorying what exists on sanghavi here*, what benchmark numbers sanghavi already achieves, and designing Phase 2 (or its equivalent) as "land sanghavi's optimizations in the new branch and go further where there's clear headroom." Not as "eliminate allocations from scratch using unified's design" — that framing assumes the work hasn't been done and risks reverting it.

### 4.3 Also on sanghavi (to port to the merged branch)

- Single-scatter approximation kernel.
- EMIT / Balsamic application scripts.
- Rich `test/benchmarks/` folder: EMIT, O2 A/B band Raman+SIF grids, OCO, CarbonI, O3 Huggins, Aerosol+Raman evaluations, RamanSIF spectra/maps/worldmaps, hydrolight comparison, Lelli comparisons.
- SIF infrastructure: the `RS_type.SIF₀` plumbing used by `creategrid_O2Aband_RamanSIF.jl`. The SIF data file (`src/SIF_emission/sif-spectra.csv`) is referenced but only `.dat~` backups are committed — needs to be committed as part of the port.
- `src/Testing/perturb_parameters.jl` — FD-Jacobian verification utility; `prototype_EMIT_aer_ht.jl` references it by absolute path but it is not committed on unified.

### 4.4 Christian's 2026-03-21 plan is stale — do not trust its code references

Christian's `uni_vSmartMOM/docs/dev_notes/sanghavi_unified_merge_plan.md` and the co-written `raman_gpu_optimization.md` were authored against a snapshot of the unified-vsmartmom branch that has since shifted. The user has confirmed the codebase was modified after those documents were written. Treat any file path, line number range, struct field listing, or specific code quotation from those documents as a *starting hint* that you must verify against the current worktree before using. In particular:

- Line-number references like `interaction_inelastic.jl:293-305` in `raman_gpu_optimization.md` may now point to different code.
- The proposed `InelasticWorkspace` struct layout in §2a of Christian's plan may differ from what's sensible today, particularly given sanghavi already has `InteractionWorkspace` threaded in (see §4.2).
- File lists like "Files to modify" in each phase section may reference files that have since been renamed, moved, or consolidated.
- Allocation-count estimates (~19k/run) were measured on the old unified snapshot; the current counts on both branches should be re-measured as part of Inventory B.

Any claim in Claude's final plan that derives from Christian's plan must be re-grounded in the current tip of the relevant worktree. Do not propagate stale references.

### 4.5 Branch relationship

- Merge base: `702cbc3d` (2022-04-07, "Added multisensor simulations with Raman") — ~4 years of divergence.
- `main` is entirely contained in `unified-vsmartmom` (zero commits unique to main). Treat `unified-vsmartmom` as the next main.
- Commits unique to `unified-vsmartmom` relative to sanghavi: 250.
- Commits unique to `sanghavi` relative to unified-vsmartmom: 50.

## 5. Required inventories (do these before planning)

Do not produce an implementation plan until the following inventories are written out and shared with the user for validation. Each inventory becomes a short markdown file on the new branch under `docs/dev_notes/`. Write the inventory files before the plan.

### Inventory A — `sanghavi_optimization_inventory.md`

For each of these commits on sanghavi, document what changed, what abstraction was introduced, where it lives in the source tree, and how it relates to any analogous design on unified-vsmartmom. Explicitly compare sanghavi's `InteractionWorkspace` against Christian's `InelasticWorkspace` proposal in `uni_vSmartMOM/docs/dev_notes/sanghavi_unified_merge_plan.md` §Phase 2a. Decide per-item: adopt sanghavi's implementation wholesale, use it as a design input for a rewrite on unified's API, or discard in favor of the unified-side design. Record the rationale.

Commits to inventory:
- `854b44c` Add Float32 support, InteractionWorkspace, and optimize get_n₀_n₁
- `d75dacb` Add per-direction CPU-staging for interaction workspace
- `3e9926f` Add Raman optimization test infrastructure
- `e89ec1c` Update baseline script and add batched ops benchmark
- `083353b` Update lin+Raman (physics content — even though we're not porting linearized Raman, the underlying forward-Raman updates are in scope)

Plus the broader question: are there other commits in the 50 unique to sanghavi that introduce optimizations or physics improvements not covered by the recent five? Scan the log and flag anything that looks non-trivial.

### Inventory B — `sanghavi_performance_baselines.md`

Capture, from the sanghavi worktree, the current measured GPU / CPU performance numbers and allocation counts for forward Raman on the representative problem sizes in `raman_batched_ops_benchmark.jl` and `raman_large_scale.jl` (the latter lives on unified but the physics it exercises is sanghavi's). These numbers become the non-regression floor for the merged branch. Any plan must commit to "`sanghavi-unified` matches or beats these" — not "`sanghavi-unified` matches unified's pre-port numbers."

### Inventory C — `sif_plumbing_status.md`

Read `vSmartMOM.jl/test/benchmarks/creategrid_O2Aband_RamanSIF.jl` and the Inelastic `RS_type` definitions on both branches. Document whether `RS_type.SIF₀` plumbing exists on unified (grep for it in `uni_vSmartMOM/src/Inelastic/`), whether sanghavi-side changes are needed to reach parity, and whether the SIF data file (`sif-spectra.csv`) is anywhere in the repo or only on disk. Report back what additional port work SIF needs.

### Inventory D — `physics_delta_summary.md`

Short audit of the non-optimization physics differences between branches. Check in particular:
- `src/CoreRT/CoreKernel/elemental_inelastic.jl`, `raman_atmo_prop.jl` — any numerical differences between sanghavi and unified.
- `src/CoreRT/CoreKernel/rt_kernel.jl` on both branches — has sanghavi added kernels or dispatches unified lacks?
- `src/Inelastic/` — types, `computeRamanZλ!`, `getRamanSSProp!`, cross-section helpers.

Anything sanghavi has that unified lacks, or any subtle bug-fix on one side that didn't propagate, goes here.

### Inventory E — `docs_state_audit.md`

For unified-vsmartmom:
- Public exports in `uni_vSmartMOM/src/vSmartMOM.jl`. Which have docstrings, which don't.
- Existing `uni_vSmartMOM/docs/` structure. What is already Documenter-ready, what is legacy dev notes.
- The 474-test suite: which tests are "slow" (candidates for a `VSMARTMOM_SKIP_SLOW` flag), which exercise Raman, which exercise linearization, which are canary / quickstart candidates.

This inventory drives the docs-overhaul phase; without it, the docs plan is guesswork.

### Inventory F — `christian_plan_delta.md`

Compare Christian's `sanghavi_unified_merge_plan.md` on unified-vsmartmom (2026-03-21) against the scope above AND against the current code. Specifically:
- **Staleness audit.** The plan and its companion `raman_gpu_optimization.md` were written against an older snapshot of unified-vsmartmom. Sample ten or so concrete references from the plan (line numbers, specific file lists, struct field counts, allocation numbers) and check each one against the current worktree. Report which still hold, which don't, and how that shifts the plan's recommendations.
- Christian's plan assumes Phase 2 eliminates allocations from scratch. Sanghavi already has `InteractionWorkspace`. Reconcile.
- Christian's plan includes Phase 5 = Linearized Raman (deferred). That phase is *removed from scope here*, not deferred. Record the reason (user's product decision).
- Christian's plan's testing strategy claims "bit-for-bit identical" after Phase 2 and Phase 3. Phase 3's 4D→3D + gathered batched_mul changes BLAS accumulation order, so bit-exactness is unlikely. Flag this.
- What else in Christian's plan needs adjustment given the corrections above.

This inventory is the diff the user will use to coordinate with Christian before any branch creation.

## 6. Open questions (ask the user before planning — use AskUserQuestion or equivalent)

Do not guess these; do not pre-specify. Ask, and block on the answers.

1. **Acceptance script list.** Which scripts from `vSmartMOM.jl/test/benchmarks/` form the regression suite? Proposed starting set: `prototype_EMIT_aer_ht.jl` (with `RS_type = noRS()`), `creategrid_O2Aband_RamanSIF.jl`, `creategrid_O2Bband_RamanSIF.jl`, `O3Huggins_polRaman.jl`, `prototype_CarbonI.jl`, plus whatever else the user selects. Ask her to confirm or amend.
2. **Spectral truncation per script.** For each acceptance script, what narrower spectral range or band subset makes the script fast enough for dev iteration? Sanghavi scripts typically read YAML configs that define bands — ask per script.
3. **Per-quantity tolerance bars.** rtol and atol for each quantity relevant to each acceptance script: I, Q, U, V radiances; ieR, ieT Raman-scattered intensities; dR, dT Jacobians (elastic only — linearized Raman is out of scope); SIF contributions.
4. **Coordination with Christian Frankenberg.** Which phases of the eventual plan go through Christian for sign-off before merging into `unified-vsmartmom`? At minimum, the Inventory A reconciliation note (before any code change) and the final merge PR. Ask whether additional checkpoints are wanted.
5. **Truncation / cleanup confirmation.** The sanghavi branch has research detritus at top level (PNGs, Manifest copies, `.dat~` backups, hardcoded `/home/sanghavi/...` paths in scripts). Confirm that these should be dropped from `sanghavi-unified`, and that the canonical versions of `src/SIF_emission/sif-spectra.csv` and `src/Testing/perturb_parameters.jl` should be committed to the merged branch.

## 7. Planning guardrails (hard rules for the plan you produce)

- **Preserve or improve sanghavi GPU Raman optimizations.** Any plan whose net effect reverts sanghavi's `InteractionWorkspace`, `get_n₀_n₁`, per-direction CPU-staging, or Raman benchmark infrastructure is wrong. Phase wording like "eliminate allocations" must be recast as "preserve sanghavi's allocation elimination where present, extend where there's headroom, land it in unified's architecture."
- **Do not include linearized Raman work.** Not a phase, not a deferred item, not a stretch goal. Treat as a dead-lettered product decision.
- **Do not pre-specify tolerances.** Even if you have sensible defaults in mind. Ask.
- **Read both worktrees.** Every claim in the plan must be citable to a file path rooted at `/home/sanghavi/code/github/vSmartMOM.jl/` or `/home/sanghavi/code/github/uni_vSmartMOM/`. No assumptions drawn from only one worktree.
- **Do not touch the elastic path.** `⊠` in `doubling.jl`, `interaction.jl`, `elemental.jl` stays as-is. Inelastic path only.
- **Produce the plan for human review before any code change.** The plan you emit must be a document the user reads, comments on, and approves. Do not start editing Julia files as part of producing the plan.
- **Bit-exactness claims must be honest.** Phases that reorder BLAS operations, change accumulation order, or change data layout cannot be bit-exact by construction. Use tolerance-based claims for those; reserve bit-exactness for phases that genuinely preserve the execution path (workspace threading with identical write order, sync removal, etc.).
- **Verify every code reference lifted from Christian's 2026-03-21 plan.** That plan and its companion `raman_gpu_optimization.md` were written against an older snapshot of unified-vsmartmom. Line numbers, file lists, struct layouts, and allocation counts in those documents may not match the current code. Re-ground every such reference against the current worktree before putting it in the plan you produce.
- **Baselines come from sanghavi, not from unified.** The regression harness compares against numbers captured *on sanghavi* for each acceptance script. Sanghavi is the physics truth source.
- **Phase 0 is a reconciliation phase, not a code phase.** No branch creation, no commits, until Inventories A–F exist and are approved by the user.

## 8. Reference documents for Claude to read before planning

All already checked in, either on the sanghavi branch or on unified-vsmartmom:

- `uni_vSmartMOM/CLAUDE.md` — architecture and terminology.
- `uni_vSmartMOM/CHANGELOG.md` — v2.0.0 breaking changes and "Known Limitations".
- `uni_vSmartMOM/docs/dev_notes/sanghavi_unified_merge_plan.md` — Christian's 2026-03-21 plan (historical; scope now superseded).
- `uni_vSmartMOM/docs/dev_notes/raman_gpu_optimization.md` — GPU allocation audit (done against unified's Raman path; may already be partially obsoleted by sanghavi's InteractionWorkspace).
- `uni_vSmartMOM/docs/dev_notes/RAMAN_CODE_HANDOFF.md` — how Raman is wired in unified; addressed to the sanghavi-branch author.
- `uni_vSmartMOM/docs/dev_notes/LINEARIZATION_BUGS.md` — bugs caught during the sanghavi → unified Jacobian port.
- `uni_vSmartMOM/docs/dev_notes/batched_kernel_benchmarks.md` — 2026-04-19 benchmark writeup.
- `/sessions/awesome-friendly-mayer/mnt/vSmartMOM/MERGE_INVESTIGATION.md` — structural investigation of both branches (evidence base, written 2026-04-19 before Codex review).
- `/sessions/awesome-friendly-mayer/mnt/vSmartMOM/SANGHAVI_UNIFIED_PLAN.md` — earlier draft plan. **Read only as background, not as the canonical plan**; its Phase 2 framing pre-dates the sanghavi InteractionWorkspace correction and is the artifact this handoff brief supersedes.

---

## 9. Paste-ready prompt for the Claude Opus Plan agent

Paste this directly into `claude --model opus /plan` (or the `/plan` slash command inside a Claude Code Opus session) on the lab server:

```
You are the planning agent for a merge between two branches of vSmartMOM.jl,
a Julia radiative transfer package. Read the handoff brief at
/home/sanghavi/code/github/vSmartMOM.jl/plans/CLAUDE_HANDOFF_BRIEF.md
(or wherever the user has placed it — confirm the path before proceeding)
and treat it as authoritative over any older plan documents.

Both branches are checked out as linked worktrees:
- sanghavi worktree: /home/sanghavi/code/github/vSmartMOM.jl/
- unified-vsmartmom worktree: /home/sanghavi/code/github/uni_vSmartMOM/

Work in three stages:

Stage 1 — Inventory. Before writing any implementation plan, produce
Inventories A–F described in §5 of the handoff brief. Write each
inventory as a short markdown file and place them in
docs/dev_notes/ on the sanghavi-unified branch (which you will
create in Stage 2, not yet). For Stage 1, you may either stage the
inventory files in a holding directory or present them inline and
commit later. Validate each inventory with the user before moving on.

Stage 2 — Ask the user the open questions in §6 of the brief and block
on the answers. Do not guess tolerances, do not pick the acceptance
script list, do not decide Christian-coordination touchpoints unilaterally.

Stage 3 — Produce the implementation plan. Structure it as phases with
explicit gates, files touched, acceptance criteria, and regression checks.
Obey the guardrails in §7 — in particular:
- Preserve or improve sanghavi's existing GPU Raman optimizations
  (InteractionWorkspace, get_n₀_n₁, per-direction CPU-staging, Raman
  benchmark infrastructure). Do not reimplement them from Christian's
  older design.
- No linearized Raman, ever.
- No tolerance pre-specification.
- No touching the elastic path.
- Plan only; do not begin code changes until the user has reviewed and
  approved the plan.

Baselines for the regression harness are captured on the sanghavi branch,
not on unified-vsmartmom. The sanghavi branch is the physics truth source.

Write the plan as a single markdown file at
/home/sanghavi/code/github/vSmartMOM.jl/plans/IMPLEMENTATION_PLAN.md,
numbered phases, cross-referencing each inventory, with an explicit
per-phase "regression check" section that confirms no sanghavi
optimization is being rolled back.

Before you start Stage 1, restate the mission back to me in two sentences
so I can confirm you have it right.
```

*End of prompt.*

---

## 10. What this brief is not

This brief deliberately does *not* contain:

- Phase-level implementation detail (no `mutable struct InteractionWorkspace` layout, no specific batched-mul call shapes, no file-by-file refactor sketch). Claude generates those after inventorying.
- Specific tolerance numbers.
- A committed acceptance script list.
- A claim about what Phase 2's exit criterion is in terms of absolute allocation counts. Those numbers become concrete after Inventory B captures sanghavi's current baselines.
- A timeline. Phase lengths depend on what Inventories A–F reveal.

The previous `SANGHAVI_UNIFIED_PLAN.md` contained most of those details. It is useful evidence and background, but its Phase 2 framing pre-dates the corrected understanding that sanghavi already has `InteractionWorkspace`. Treat that document as a draft that informs the plan Claude produces, not as a plan to hand to Claude.
