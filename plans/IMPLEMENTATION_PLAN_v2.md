# sanghavi-unified merge — implementation plan v2

**Authored 2026-04-19. Status: draft, awaiting user approval before any code lands.**

**Supersedes** `plans/IMPLEMENTATION_PLAN.md` dated 2026-04-19. Consolidates `plans/PLAN_AMENDMENTS_2026-04-19.md` §1–§7. On any conflict between this document and v1, v2 governs. On any conflict between this document and the inventories, the amendments doc governs (via v2, since v2 mirrors it).

Worktree tips as of this draft:
- sanghavi: `/home/sanghavi/code/github/vSmartMOM.jl/` @ `9ee9a75` (branch `sanghavi`)
- unified-vsmartmom: `/home/sanghavi/code/github/uni_vSmartMOM/` @ `a4e4187` (branch `unified-vsmartmom`)

---

## Context

vSmartMOM.jl has two active branches with ~4 years of divergence since merge base `702cbc3d`:

- `unified-vsmartmom` (250 commits unique) owns the v2.0.0 architecture: `RTModel{ARCH,FT}` hierarchy, `SolverConfig`/`Atmosphere`/`Optics`/`ObsGeometry`/`QuadPoints` sub-structs, Aerosols/IO modules, 474-test suite, Cox-Munk/canopy/RossLi/RPV/Rahman/fresnel surfaces, HITRAN direct download + version tracking, Float32/Float64 flexibility, GPU-accelerated NAI2 Mie, GPU batched-kernel extension, forward Raman physics partially ported from sanghavi.
- `sanghavi` (50 commits unique) owns the physics truth source for Raman + several GPU optimizations landed on top: `InteractionWorkspace` struct with per-direction CPU-staging, `get_n₀_n₁` arithmetic rewrite, Float32 support, Raman benchmark infrastructure, and forward-physics corrections that did NOT propagate to unified via commit `59f8de8` (see Inventory D §1 & §3).

The goal is to branch `sanghavi-unified` off `unified-vsmartmom` HEAD, port sanghavi's remaining contributions on top (physics fixes + optimizations + SIF + scripts + `perturb_parameters.jl` + single-scatter driver), produce a release-quality doc set, ship a regression harness that pins sanghavi's numbers as the non-regression floor, and eventually promote `sanghavi-unified` to `main`.

**Hard constraints** (from handoff brief §3 and §7, and amendments §7):
- Linearized Raman is permanently out of scope — no `AddedLayerRS`-style constructors, no Raman branch in `rt_kernel_lin.jl`.
- Elastic path (`⊠` notation in `doubling.jl`, `interaction.jl`, `elemental.jl`) stays as-is.
- Sanghavi's existing GPU Raman optimizations must be preserved or improved — never rolled back.
- Baselines come from sanghavi at current tip, not from unified.
- No pre-specifying tolerances. Per-script, per-quantity, at execution time.
- Produce the plan for human review before any code change.

---

## Authority rule (binding, applies to every phase)

**Sanghavi is the authority for the inelastic path. Full stop.** (amendments §1)

"Inelastic path" means, concretely:

- `src/Inelastic/` in full. All files, including nested sub-module structure. `InelasticScattering` is a sub-module (not a separate registered package), so this is a file copy, not a dependency bump.
- All `*_inelastic.jl` files in `src/CoreRT/CoreKernel/`: `interaction_inelastic.jl`, `elemental_inelastic.jl`, `elemental_inelastic_plus.jl`, `doubling_inelastic.jl`, `raman_atmo_prop.jl`.
- The inelastic-specific portions of `src/CoreRT/atmo_prof.jl`, `src/CoreRT/model_from_parameters.jl`, `src/CoreRT/types.jl`.
- The `RS_type <: noRS` branches and `InteractionWorkspace` allocation inside `src/CoreRT/rt_run.jl`. `rt_run.jl` itself is a merge, not a copy: unified's `RTModel` scaffolding plus sanghavi's inelastic branches.

The authority rule settles, by fiat, every per-item decision that the inventories framed as "adopt wholesale / rewrite / discard" on inelastic files. Per-item rationale is retained in the inventories for historical reference only; the answer for any inelastic file is "sanghavi."

**Exceptions to the rule** (explicit carve-outs from amendments §2):

1. `reduce_profile` default flips — see §Phase 1a.
2. `apply_D_matrix_elemental!` scalar shortcut from commit `9a26002` — dropped per §Phase 1b (amendments §2.2).
3. `hem_R` / `hem_T` hemispheric API change from commit `9a26002` — deferred post-merge per §Phase 1b (amendments §2.3).
4. `inelastic_helper_old.jl` on sanghavi is a 924-line frozen reference file; does not port.

---

## Branch strategy

```
unified-vsmartmom @ a4e4187
        │
        └── sanghavi-unified      (created at Phase 0 exit, after inventories + this plan are user + Christian approved)
```

No rebasing. No merging sanghavi wholesale. Port specific commits/files/patterns as described in each phase. Commits on `sanghavi-unified` use their own short messages; co-authorship cited to sanghavi where physics or optimization code is imported.

---

## Acceptance-script configuration

Both scripts need reduced spectral ranges and possibly reduced vertical-layer counts for dev iteration. Concrete numbers are phase-gated: the YAML config edits land in **Phase 3 (script port)** and are filled in at that time based on what runs in <60 s on the reference GPU. Placeholder table:

| Script | YAML baseline | Reduced range | Reduced n_layers | Notes |
|---|---|---|---|---|
| `prototype_EMIT_aer_ht.jl` | inherits from its YAML | TBD @ Phase 3 | TBD @ Phase 3 | noRS; aerosol profile retained |
| `creategrid_O2Aband_RamanSIF.jl` | O2 A-band ~755–775 nm, ~72 layers | TBD @ Phase 3 | TBD @ Phase 3 | RRS + SIF grid; coarsest grid point only during dev |

The full-range versions remain runnable (separate YAML or a `--full` flag) for final validation at Phase 8.

---

## Tolerance table (TBD — user fills during Phase 2 draft)

Per handoff brief §7 and amendments §7, tolerances are a user-facing choice per-script, per-quantity. Placeholders to be populated before Phase 2b baselines are committed:

| Quantity | `prototype_EMIT_aer_ht.jl` (noRS) | `creategrid_O2Aband_RamanSIF.jl` (RRS+SIF) |
|---|---|---|
| I | rtol=?, atol=? | rtol=?, atol=? |
| Q | rtol=?, atol=? | rtol=?, atol=? |
| U | rtol=?, atol=? | rtol=?, atol=? |
| V | rtol=?, atol=? | rtol=?, atol=? |
| ieR (Raman reflected) | n/a | rtol=?, atol=? |
| ieT (Raman transmitted) | n/a | rtol=?, atol=? |
| dR (Jacobian wrt R) | rtol=?, atol=? | n/a (linearized Raman out of scope) |
| dT (Jacobian wrt T) | rtol=?, atol=? | n/a |
| SIF contribution (pattern; see Phase 3) | n/a | rtol=?, atol=? |
| Wall-clock (warm) | regression ≤ +X% of sanghavi | ≤ +X% |
| GPU alloc count | ≤ sanghavi count | ≤ sanghavi count |
| GPU peak memory | ≤ sanghavi peak + Y MB | ≤ sanghavi peak + Y MB |

Per invariant 6 (below), tolerance claims must be honest: bit-exactness applies only to phases that preserve execution path verbatim; phases that reorder accumulation (staging, 4D→3D batching) use the tolerances above.

---

## Phase 0 — Reconciliation (no code)

**Purpose:** Ensure all parties agree on scope before `sanghavi-unified` exists. No branch, no commits.

**Preconditions / Gate:** None. This is the starting phase.

**Actions:**
1. User reads inventories A–F at `plans/inventories/` and this plan v2.
2. **Editorial pass on inventories** (per amendments §6) before Christian sees them:
   - Inventory A §5: reframe "Adopt wholesale for items 1–10" as "these items come with the wholesale Inelastic port per the authority rule." §8 Q1 (staging) resolved: `staged::Bool = true` default. §8 Q4 (Bodhaine default): answered in Phase 1a. §8 Q5 (reduce_profile default): answered in Phase 1a. §8 Q6 (delete `inelastic_helper_old.jl`): yes, by not porting.
   - Inventory D §6 (port list): reframe P0/P1/P2 items as "these land via the Inelastic port per the authority rule." §7 Q1/Q2 on `α̅` unit convention: closed per amendments §2.4. §7 Q5 on canopy+Raman: closed per amendments §2.6.
   - Inventory F §5 ("Final delta summary"): rewrite the Phase 2 recast wording so sanghavi's workspace design is the base *by authority*, not by reconciliation. Christian's doubling/pivot additions are candidates to layer on, not a synthesis.
   - Inventory 0_SUMMARY finding #3: rewrite as "Build the merged workspace to sanghavi's `InteractionWorkspace` layout adapted to `RTModel`. Christian's proposed additional fields layer on if compatible. `staged::Bool = true` default."
   - Inventory 0_SUMMARY finding #4: Phase 3 bit-exactness claim is already correctly flagged as false. No edit needed.
   - These passes are editorial, not factual. The evidence in each inventory stands.
3. `σ_Rayl_coeff` storage convention confirmed via grep over `src/Inelastic/src/` for future-maintainer clarity. The authority rule resolves which code lands; this confirms which convention the code assumes, for documentation purposes.
4. User flags any remaining factual corrections; we patch the relevant inventory in place.
5. User forwards the edited Inventory A + Inventory F + this plan v2 + `plans/PLAN_AMENDMENTS_2026-04-19.md` to Christian Frankenberg and collects explicit ack.

**Regression check:** n/a (no code).

**Christian sign-off:** YES — blocking. Phase 0 does not exit without Christian's written ack of the editorial pass and the amendments doc.

**Exit criterion:** User approves this plan; Christian has acked the Inventory A + F delta, the amendments document, and the editorial inventory pass in writing.

**Estimated scope:** 1–3 days of human calendar; all coordination.

---

## Phase 1 — Physics truth sync + Inelastic port

**Purpose:** Land the forward-physics defaults that sanghavi carries, bulk-port the inelastic path from sanghavi per the authority rule, port the `rt_run_ss` driver, wire in the Aerosols module, and port `perturb_parameters.jl`. Do all of this BEFORE baseline capture so baselines reflect correct physics.

Phase 1 is structured as five sub-phases, in order. Each sub-phase is a logical commit unit (possibly multiple commits internally).

**Preconditions / Gate:** Phase 0 exited. Branch `sanghavi-unified` is created off `unified-vsmartmom @ a4e4187`.

### Phase 1a — Physics default homogenization

**Scope:** Two default switches land together in one logical commit.

1. Migrate `getRayleighLayerOptProp` (or wherever the old `0.00864 × …` fit lives on unified) to Bodhaine 1999 Eq. 30.
2. Migrate `reduce_profile` default to linear interpolation on uniform pressure half-levels (sanghavi's current default). Preserve bin-averaging as `reduce_profile_binavg` (or equivalent) available via an opt-in keyword (amendments §2.1).
3. Grep the 474-test suite for reference values that depend on Rayleigh formula or layer averaging (most forward-RT tests, natraj, 6SV1, any multi-layer profile config).
4. Re-baseline affected tests. Review diff as a unit. Spot-check ~10 tests; verify shifts are consistent with combined Bodhaine + interpolation expectations; commit the batch.
5. Land this `reduce_profile` refactor TODO comment at the call site:

   ```julia
   # TODO: reduce_profile and getRayleighLayerOptProp both deserve refactoring
   # beyond this merge. The current defaults are physics-forward (Bodhaine
   # 1999 Eq. 30, interpolated profile) but the implementation layout mixes
   # legacy and modern code paths via keyword args and name suffixes. Cleaner
   # architecture: a single configurable ProfileReduction strategy (dispatch on
   # type) and a single RayleighFormula strategy, both YAML-configurable.
   # Tracked as a future PR; not a merge blocker.
   ```

**Commit as one logical unit:** "Migrate Rayleigh formula and profile reduction defaults; re-baseline regression tests." Splitting the commits creates intermediate states where some tests pass on mixed-old-and-new defaults, making the diagnostic harder.

**Regression check:** 474 tests pass after re-baselining. Spot-check documentation records which tests were re-baselined and the reviewer's initials. No sanghavi optimization present yet; no rollback risk.

### Phase 1b — Inelastic port (sanghavi authority)

**Scope (per authority rule §1 + amendments §§2.2–2.6):**

- Bulk-copy `src/Inelastic/` from sanghavi. Confirm the `include` path in `src/vSmartMOM.jl` resolves (nested `src/Inelastic/src/` structure is retained).
- Copy all `*_inelastic.jl` files in `src/CoreRT/CoreKernel/` from sanghavi: `interaction_inelastic.jl`, `elemental_inelastic.jl`, `elemental_inelastic_plus.jl`, `doubling_inelastic.jl`, `raman_atmo_prop.jl`.
- Adopt sanghavi's inelastic-specific portions of `src/CoreRT/atmo_prof.jl`, `src/CoreRT/model_from_parameters.jl`, `src/CoreRT/types.jl`.
- `src/CoreRT/rt_run.jl` is a **merge**, not a copy: unified's `RTModel` scaffolding, sanghavi's `RS_type` branches, sanghavi's `InteractionWorkspace` allocation site. The actual `InteractionWorkspace` *landing and threading* is Phase 4, not Phase 1. Phase 1 leaves `_interaction_ws = nothing` at the allocation site and the kernel falls back to `similar()` per sanghavi's existing `workspace === nothing` branch.
- Port the rest of commit `9a26002`'s non-`_lin` content **minus** the D-matrix scalar shortcut (amendments §2.2) and **minus** `hem_R`/`hem_T` (amendments §2.3). EMIT scripts and the Float32 YAML config port normally.
- Land the `α̅` verification comment from amendments §2.4 at the port site in `src/Inelastic/src/inelastic_cross_section.jl` (the planning agent should insert the actual verification source — notebook, paper section, commit message, email thread — during this sub-phase):

  ```julia
  # Frequency correction (1 - (c·ν_eff/ω₀)²). ν_eff is wavenumber (cm⁻¹),
  # ω₀ is stored in wavenumber units; no 2π factor. Verified correct against
  # Buldakov et al. 1996 Eqs. 36a-39b by [author/source].
  ```

- Land the canopy + Raman TODO comment from amendments §2.6 (find the right anchor during implementation — likely `src/CoreRT/CoreKernel/rt_kernel.jl` or `src/CoreRT/Surfaces/canopy_surface.jl`):

  ```julia
  # TODO: Canopy + Raman (RRS/VS) coupling is not currently implemented. The
  # canopy BRDF code is noRS-only; calling rt_run with CanopySurface + RRS has
  # undefined behavior. Future work: couple canopy to the inelastic path.
  ```

- Do NOT port `inelastic_helper_old.jl` (924-line frozen reference file on sanghavi).

**Regression check:** `test_forward_raman.jl` (currently orphan on unified — wire in during this sub-phase) matches sanghavi reference within user tolerance on a toy 1-band RRS config. Reference JLD2 committed to `test/reference/phase1b_inelastic_port.jld2`. No sanghavi optimization rolled back — `InteractionWorkspace` threading is Phase 4; this sub-phase preserves sanghavi's `workspace === nothing` fallback path.

### Phase 1c — `rt_run_ss` driver port

**Scope:**

- Copy the `rt_run_ss` driver body from sanghavi's `src/CoreRT/rt_run.jl` (~lines 258 and 450) into unified's `src/CoreRT/rt_run.jl`.
- Adapt `vSmartMOM_Model` field access to `RTModel` sub-struct access.
- Add smoke test `test_forward_ss.jl` and wire into `runtests.jl`.

**Regression check:** `test_forward_ss.jl` passes. The export leak (`rt_run_ss` exported on unified but undefined, per Inventory D §5) is closed.

### Phase 1d — Aerosols module wire-in (per amendments §2.7)

**Scope:**

- Add `include("Aerosols/Aerosols.jl")` to `src/vSmartMOM.jl`.
- Add this WIP header comment to `src/Aerosols/Aerosols.jl`:

  ```julia
  # WIP: This module is user-facing but the API may evolve in follow-up PRs.
  # Landed as part of the sanghavi-unified merge to close the
  # `using vSmartMOM.Aerosols` export gap; further cleanup is a known
  # follow-up workstream.
  ```

- Wire `test_Aerosols.jl` into `runtests.jl`.

**Regression check:** `using vSmartMOM.Aerosols; <one exported symbol>` resolves. `test_Aerosols.jl` runs and passes. The live bug (Aerosols exported without `include`, per Inventory E) is closed.

### Phase 1e — `perturb_parameters.jl` port

**Scope:**

- Port `src/Testing/perturb_parameters.jl` from sanghavi (handoff brief §4.3).
- One smoke test.

This is carried forward from v1's Phase 3; it fits better here adjacent to the physics port.

**Regression check:** Smoke test passes.

### Phase 1 exit

**Christian sign-off:** Optional (already covered by Phase 0 ack). Physics fixes are uncontroversial by the authority rule.

**Exit criterion:** All five sub-phases complete; 474 tests pass post-re-baselining; orphan tests wired where in-scope land; reference JLD2 committed; Aerosols module callable end-to-end.

**Estimated scope:** ~1.5–2 weeks. The 474-test re-baseline in 1a is the biggest single task.

---

## Phase 2 — Baseline capture + regression harness

**Purpose:** Capture the non-regression floor on `sanghavi-unified` (post-Phase-1) AND on the sanghavi branch itself. Before Phase 4 can gate on "≤ sanghavi's allocation floor," we need numbers on both.

Phase 2 is structured as two sub-phases.

**Preconditions / Gate:** Phase 1 complete. Tolerance table in this doc populated by user.

### Phase 2a — Diagnostic baseline (falsification test for Inventory D)

**Scope:**

- Run the harness on unmodified `sanghavi-unified` (= post-Phase-1) vs sanghavi.
- Record the delta.

This verifies Inventory D's physics-fix list was complete: if the post-Phase-1 delta exceeds what D predicted, D was incomplete and must be revisited before Phase 3 proceeds. This is a diagnostic step, not a gate on a number; it's a falsification test for Inventory D.

**Files touched:**
- `test/benchmarks/harness/run_benchmarks.jl` (new) — wraps each acceptance script, captures wall-clock, CPU allocs, GPU used-delta via `CUDA.Mem.info()` before/after, GPU peak via `CUDA.@profile trace=true`, kernel-launch count, I/Q/U/V and ieR/ieT arrays.
- `test/benchmarks/harness/metrics.jl` (new) — metric computation + dumping.
- `test/benchmarks/harness/report.jl` (new) — loads two baseline JSONs and produces a side-by-side tolerance check report.
- `test/benchmarks/harness/scenarios.toml` (new) — declares which scripts + which reduced YAML configs to run.

**Harness requirements (per amendments §4 Phase 2):**
- Sample count N ≥ 5 for wall-clock measurements. Report median + MAD; trim outliers.
- `CUDA.math_mode!(:strict)` (or equivalent) for bit-exactness checks to avoid cuBLAS algorithm nondeterminism.
- Budget ~30 min per full-range RRS warm run in capture planning; reduced YAMLs target <60 s per warm run.
- Primary gate metric: **GPU used-memory delta.** Secondary: CPU allocation count. Others (kernel launches, peak memory) diagnostic.
- Harness records GPU/CPU identifier in output JSON; re-capture if hardware changes.

**Regression check:** Harness invocation is idempotent and deterministic (warmed, N ≥ 5 runs, median reported). Physics correctness: sanghavi-unified matches sanghavi within user tolerance on I/Q/U/V/ieR/ieT for both acceptance scripts. **This gate must pass before Phase 2b.**

### Phase 2b — Baseline capture

**Scope:**

Once Phase 2a passes review, freeze the baselines. These are the non-regression floor for Phase 4 and Phase 8.

**Files touched:**
- `test/benchmarks/baseline_output/sanghavi_9ee9a75/` (new) — committed JSON baselines from the sanghavi worktree.
- `test/benchmarks/baseline_output/sanghavi-unified_<sha>/` (new) — committed JSON baselines on the merged branch.

**Actions:**
1. Run harness on sanghavi-worktree. Commit outputs under `baseline_output/sanghavi_9ee9a75/`.
2. Run harness on sanghavi-unified (still no workspace port yet; Raman code runs via unified's existing path with Phase 1 physics fixes applied). Commit outputs under `baseline_output/sanghavi-unified_<commit>/`.
3. Record expected delta narrative: at this point sanghavi-unified's wall-clock and alloc counts will be **worse** than sanghavi's — that's the gap Phase 4 closes.

**Regression check:** `report.jl` green on physics correctness. No sanghavi optimization rolled back (sanghavi's optimizations aren't on the branch yet).

### Phase 2 exit

**Christian sign-off:** None needed.

**Exit criterion:** Two baseline directories committed; `report.jl` green on physics; user reviews the wall-clock / alloc gap numbers and confirms the Phase 4 target.

**Estimated scope:** 3–5 days to build harness, plus re-run time on both branches.

---

## Phase 3 — SIF + benchmark script port

**Purpose:** Land the 2-line SIF injection, commit SIF data files, port the two acceptance scripts with reduced YAML configs, document the `0.5π / max` rescaling.

**Preconditions / Gate:** Phase 2 harness green.

**Files touched:**

| File | Change | Source |
|---|---|---|
| `src/CoreRT/Surfaces/lambertian_surface.jl` | Add `(1/π) * SIF₀` injection into `J₀⁻` (sanghavi lines 67–68 and 157–158; `_lin` variant stays as-is per out-of-scope policy). | Inv. C §2 |
| `src/SIF_emission/sif-spectra.csv` | Commit (14 KB). | Inv. C §3 |
| `src/SIF_emission/PC1_SIFSpectra_allSpecies.csv` | Commit (8 KB). | Inv. C §3 |
| `src/SIF_emission/ficus_refl.dat` | Commit (35 KB). Drop `.dat~` backup. | Inv. C §3 |
| `src/SIF_emission/ficus_refl_600to800nm.dat` | Commit (3 KB). Drop `.dat~` backup. | Inv. C §3 |
| `src/SIF_emission/sif_loader.jl` | New helper — `readdlm`+unit-conversion boilerplate into a single function. Export from `vSmartMOM`. | Inv. C §6 (optional recommendation) |
| `test/benchmarks/prototype_EMIT_aer_ht.jl` | Port from sanghavi; rewrite hardcoded paths to `joinpath(pkgdir(vSmartMOM), ...)`; add reduced-range YAML variant. | Inv. C §4 |
| `test/benchmarks/creategrid_O2Aband_RamanSIF.jl` | Port; rewrite paths; add reduced-range YAML variant. Keep `SIF₀ .= 0.0` override togglable. | Inv. C §4 |
| Associated YAML configs in `test/test_parameters/` | Ported + reduced variants (populate spectral range + layers during this phase). | |

**Actions:**
1. Port the 2-line Lambertian SIF injection. Small commit.
2. Commit the four real SIF data files. Separate commit. Drop `.dat~` backups in the same commit.
3. Add `sif_loader.jl` helper and route both benchmark scripts through it.
4. Before porting `creategrid_O2Aband_RamanSIF.jl`, land this comment above the rescaling line (amendments §2.5):

   ```julia
   # TODO: The (0.5π / maximum(J_SIF)) rescaling is an intentional hack to make
   # SIF magnitude data-independent for grid generation. This normalizes shape
   # but discards absolute physical magnitude. Revisit: confirm the downstream
   # physics depends only on SIF shape (not absolute flux), or replace with
   # physical units (mW/m²/cm⁻¹). Not a merge blocker.
   ```

5. Port the two acceptance scripts. Determine reduced spectral ranges and layer counts empirically: pick ranges such that **one full warm run completes in <60 s on the reference GPU** (user-specified hardware; record in harness README). Commit reduced YAMLs alongside full-range YAMLs.
6. Re-run harness on sanghavi-unified. Update `baseline_output/sanghavi-unified_<commit>/`.

**Regression check:**
- Acceptance scripts now complete on sanghavi-unified.
- **SIF correctness is a *pattern* check** (ieR responds when SIF is nonzero), not an absolute-magnitude check, until the rescaling is revisited (amendments §2.5).
- `rt_run_ss` smoke test still passes.
- Wall-clock / alloc delta vs sanghavi captured but not yet closed.
- No sanghavi optimization rolled back (still nothing to roll back).

**Christian sign-off:** None needed.

**Exit criterion:** Both acceptance scripts run on sanghavi-unified to completion in <60 s each (reduced YAML) and produce outputs within user tolerance vs sanghavi reference.

**Estimated scope:** ~1 week.

---

## Phase 4 — Workspace landing

**Purpose:** Land sanghavi's `InteractionWorkspace`, adapted to unified's `RTModel` architecture, threaded through `rt_run.jl`, `elemental_inelastic.jl`, `doubling_inelastic.jl`, `interaction_inelastic.jl`. Christian's `InelasticWorkspace` proposal contributes candidate additional fields that layer on top *if compatible*; sanghavi's design is the base **by authority**.

**Preconditions / Gate:**
- Phase 3 complete and green.
- **Christian has reviewed** the sanghavi-workspace-adapted-to-`RTModel` design and signed off on any additional fields to layer on. Framing: "here's sanghavi's `InteractionWorkspace` ported to unified's architecture; these are the additional fields from your `InelasticWorkspace` proposal we've layered on; any more?" Review of additions, not a design review — the base design is sanghavi's by authority.

**Design (per authority rule + amendments §4 Phase 4):**
1. Base = sanghavi's `InteractionWorkspace` 15-field layout. Additions from Christian's `InelasticWorkspace` proposal (doubling-tmp fields, `batch_inv!` pivot/info sharing, shared 3D buffers with elastic `RTWorkspace`) layer on top **if compatible**. If not compatible, sanghavi wins.
2. CPU-staging from sanghavi's `d75dacb`: `staged::Bool = true` **default-on** (per amendments §2 / §4 Phase 4). Phase 5 may eventually replace the memory win with batching, but Phase 5 is conditional and may exit as a no-op. Don't trade a certain 3.5 GB win for an uncertain one.
3. Allocation site = `rt_run` (not `model_from_parameters`), consistent with sanghavi. Rationale: workspace lifetime matches a single forward RT solve; `rt_run` knows the problem size.
4. Scattering-interface coverage = extend to `ScatteringInterface_00`, `ScatteringInterface_11`, `ScatteringInterface_10`, `ScatteringInterface_01` if any of these fall within sanghavi's authority surface. Extend beyond sanghavi's `_11`-only coverage if and only if doing so is additive with sanghavi's design, not a redesign.
5. Thread through all Raman kernel call sites. Elastic path untouched (per hard constraint).

**Files touched:**
- `src/CoreRT/types.jl` — `InteractionWorkspace` struct definition (ported from sanghavi; adapted field types to `RTModel`'s `FT` / array-type parameters).
- `src/CoreRT/rt_run.jl` — allocate + pass through (replaces the `_interaction_ws = nothing` placeholder landed in Phase 1b).
- `src/CoreRT/CoreKernel/interaction_inelastic.jl` — consume workspace in place of `similar()` allocations (matching sanghavi's design).
- `src/CoreRT/CoreKernel/elemental_inelastic.jl` — same.
- `src/CoreRT/CoreKernel/doubling_inelastic.jl` — same; adds the doubling-tmp fields from Christian's design that sanghavi doesn't have yet, if they layer on compatibly.

**Actions:**
1. Commit workspace struct + allocator (no call-site changes yet — should compile + tests should still pass).
2. Commit per-kernel wiring incrementally (one kernel file per commit). After each, re-run harness and confirm physics output matches Phase 3 tolerance.
3. Commit the `staged::Bool` default-on path + dedicated test exercising both `staged=true` and `staged=false` (`test_workspace_staged.jl`).
4. Re-run harness on sanghavi-unified, update `baseline_output/sanghavi-unified_<commit>/`.

**Phase 4 can run in parallel with a subset of Phase 3** — workspace-side work doesn't touch the SIF Lambertian injection or the benchmark-script ports. If resources permit, parallelize.

**Regression check:**
- **Sanghavi's `InteractionWorkspace` capabilities are preserved.** `ScatteringInterface_11` RRS path shows ≤ sanghavi's allocation count. If it shows MORE allocations than sanghavi, Phase 4 is not done — root-cause and fix before moving on. This is the single most important regression gate in the entire plan.
- CPU-staging (`staged=true` default) reproduces sanghavi's ~3.5 GB savings on a large FP64 config. Record measured savings.
- Doubling-tmp fields not present on sanghavi — new ground. Confirm no numerical regression vs Phase 3 output (tolerance-based; Phase 3 physics is reference).
- Wall-clock ≤ sanghavi's non-staged wall-clock within user tolerance. Staged mode: ≤ sanghavi's staged wall-clock within tolerance.
- Bit-exactness claim per invariant 6: only for the workspace threading commits that preserve execution path verbatim. Not for the staging mode, which serializes mid-kernel.

**Christian sign-off:** YES — blocking. Review the adapted struct design before Action 1.

**Exit criterion:** Harness green on physics + allocation counts ≤ sanghavi's floor + CPU-staging default-on works.

**Estimated scope:** 2–3 weeks. This is the workhorse phase.

---

## Phase 5 — Further GPU headroom (conditional)

**Purpose:** Extend beyond sanghavi's optimization set IF Phase 4 numbers show clear headroom. Candidates identified by inventories:

- Doubling-side allocations (`tmp3–tmp6`, `gp_refl`, `ieJ₁⁺/⁻`) not yet absorbed into workspace.
- `batched_mul!` vs `⊠` micro-benchmark re-run (sanghavi's `e89ec1c` showed 5.5× slowdown for 15×15 GEMM — **do NOT go the batched_mul direction on that shape**). But there may be shapes or call sites where `batched_mul!` wins.
- Sync-point reduction on GPU (non-staged path).
- 4D→3D flattening + gather-based batched_mul — **only if re-benchmarking proves it's faster for the actual problem shapes**.

**Preconditions / Gate:**
- Phase 4 complete and green.
- **Christian signs off on the specific optimization(s) to attempt** and the measurement protocol, per-optimization.
- Phase 4 measured numbers show actual headroom above sanghavi's floor.

**Actions:**
1. Re-benchmark `batched_mul!` vs `⊠` on CURRENT Phase-4 shapes (not sanghavi's 15×15 specifically; the merged branch's `creategrid_O2Aband_RamanSIF.jl` may hit different shapes). Commit the benchmark + its output.
2. For each candidate that shows empirical win: land it as a separate commit with its own micro-benchmark entry and harness re-run.
3. If no candidate shows a win, Phase 5 exits as a no-op with a committed `dev_notes/phase5_headroom.md` explaining why further optimization isn't justified. This is a valid Phase 5 outcome.

**Regression check:**
- Each optimization commit: allocation count ≤ pre-commit count; wall-clock ≤ pre-commit wall-clock; physics output within user tolerance.
- **Sanghavi-floor check:** running allocation total must remain ≤ sanghavi's floor, regardless of what Phase 5 attempts.
- Bit-exactness claim per invariant 6: NOT for reorder-based phases. Tolerance-based only.

**Christian sign-off:** YES — blocking, per-optimization.

**Exit criterion:** Either (a) one or more commits landed with measured wins and no regressions, or (b) a no-op with a dev note explaining why.

**Estimated scope:** Variable — could be zero (no headroom) or 1–2 weeks (multiple wins).

---

## Phase 6 — Remaining scripts port + detritus cleanup

**Purpose:** Port the rest of sanghavi's benchmark scripts (OCO, CarbonI, O3 Huggins, RamanSIF variants, retrieval comparisons, EMIT/Balsamic applications), drop all research detritus.

**Preconditions / Gate:**
- Phase 5 exited.
- **Pre-flight deliverable:** `plans/phase6_script_port_list.md` — a concrete table drawn from Inventory B §1a listing each remaining sanghavi benchmark script with its runtime, data-file dependencies, and keep/port/drop decision. **User signs off on the list before porting starts.** "We thought 8, it's actually 23" is a real failure mode for research-code ports (amendments §4 Phase 6).

**Files touched:**
- `test/benchmarks/` — port remaining Raman + SIF scripts from sanghavi per the user-approved `phase6_script_port_list.md`. **Paths rewritten** to `joinpath(pkgdir(vSmartMOM), ...)`. Scripts NOT on the acceptance list still must compile + run on reduced YAMLs.
- `src/` — remove `rt_run_bck.jl`, `inelastic_helper_old.jl`, any `*_old.jl` / `*_bck.jl` files.
- Top-level `.png` files on sanghavi — not ported.
- `Manifest copy*.toml` files — not ported.
- `.dat~` backup files — not ported (already dropped in Phase 3 for SIF; cover remainder here).
- In-source prototype/plot folders — not ported.
- `.gitignore` — confirm excludes `*.jld`, `*.jld2`, top-level `*.png`, `Manifest copy*.toml`.

**Actions:**
1. Produce `plans/phase6_script_port_list.md`. User signs off.
2. Port each approved script with absolute-path rewrites. Commit scripts in logical groupings (SIF variants together, retrieval comparisons together, etc.).
3. Ensure every ported script has a reduced-YAML entry usable from the Phase 2 harness.
4. One cleanup commit removes detritus files not ported.
5. `.gitignore` tightened.

**Regression check:**
- Every ported script runs to completion on reduced YAML without errors.
- Acceptance scripts (Phase 3) still green.
- No sanghavi optimization rolled back (Phase 4's workspace remains intact).

**Christian sign-off:** None needed.

**Exit criterion:** `test/benchmarks/` contains the full set of ported scripts per the approved list; detritus is gone; `.gitignore` updated.

**Estimated scope:** ~1 week.

---

## Phase 7 — Docs overhaul

**Purpose:** Release-quality documentation: README rewrite, Documenter user guide, curated examples, API reference audit, narrative tutorial, migration guide (old-sanghavi API → new-unified API).

**Preconditions / Gate:** Phase 6 exited.

### Pre-flight (do this at the start of Phase 7, per amendments §4 Phase 7)

1. **Tutorial compatibility check.** Verify each `Tutorial_*.jl` under `docs/src/pages/tutorials/` runs against current `sanghavi-unified` tip. Literate.jl tutorials that error through the `Base.getproperty` shim need fixing before `warnonly=false` flips.
2. **`warnonly=true` removal dry-run.** Run `julia --project=docs docs/make.jl` *without* `warnonly=true` once and count actual errors. Scope the backfill work from the actual count — do not use v1's "30–50" estimate.

**Files touched:**
- `README.md` — rewrite. Cover install, quick start, `RTModel` + forward/linearized RT, surface types (one subsection each), Raman + SIF, GPU selection, Float32, HITRAN-2024 edition switching, link to examples + docs.
- `CHANGELOG.md` — **deepen** the v2.0.0 (or v2.1.0) entry beyond the `vSmartMOM_Model → RTModel` paragraph. Cover HITRAN-2024, Float32, GPU Mie, Canopy, GEOSChem, batched kernels, merged Raman from sanghavi (per amendments §4 Phase 7).
- `docs/make.jl` — remove `warnonly=true`; fix every missing `@docs` reference the dry-run surfaced. Add `@docs` entries for `RTModelLin`, `SolverConfig`, `Optics{,Lin}`, `ParameterLayout` accessors, `FwdMode`/`LinMode`, SolarModel exports.
- `docs/src/` — new pages: (a) migration guide (old-API → new-API pair for every converted script — at minimum the 2 acceptance scripts, with diffs), (b) Raman + SIF narrative tutorial, (c) per-surface quickstart, (d) GPU quickstart.
- `examples/` — 5 new examples covering the golden path:
  1. `lambertian_forward.jl` (fastest warm-up; committed JSON reference output).
  2. `coxmunk_polarized.jl` (ocean polarization).
  3. `forward_jacobians.jl` (linearized RT elastic).
  4. `rrs_sif.jl` (Raman + SIF; flagship showcase).
  5. `gpu_cuda.jl` (GPU quick-start).
  Each example: `.jl` + narrative `.md` + committed JSON/CSV reference output. Tests wire each into a quick regression gate.
- `src/Aerosols/aerosol_integration_example.jl` — revive now that the module is loaded (Phase 1d closed the `include` gap).
- **Docstring gaps** (from dry-run count): fill `AbstractPolarizationType`, `NAI2`, `PCW`, Absorption error-function structs, `FwdMode`, `LinMode`, `RTModelLin`, `SolverConfig`, `Optics`, etc.
- `docs/dev_notes/CUDA_SETUP.md` → migrate to `docs/src/advanced/gpu.md`.
- `docs/dev_notes/JACOBIAN_TEST_WORKFLOW.md` → migrate to `docs/src/advanced/jacobians.md`.
- `docs/dev_notes/*_stale.md` (see Inv. E §3) — archive to `docs/dev_notes/archive/` or delete, per user preference at Phase 7 kick-off.
- `test/runtests.jl` — wire in remaining orphan test files (those not wired in Phase 1): `test_forward_lin.jl`, `test_hybrid_ad.jl`, `test_jacobians_GPU.jl`, `test_mie_gpu.jl`, `test_performance.jl`. Triage each first; if any are stale, fix before wiring. (`test_Aerosols.jl`, `test_forward_raman.jl`, `test_forward_ss.jl` were wired in Phase 1.)

**Actions:**
1. Run the pre-flight tutorial check + `warnonly=true` dry-run. Scope remaining actions from the actual error count.
2. Turn off `warnonly=true`; fix the cascade of missing-docstring / bad-`@ref` errors.
3. Write the 5 examples with reference outputs.
4. Rewrite README; update CHANGELOG with the deepened v2.x.0 entry.
5. Author migration guide + narrative tutorial.
6. Triage + wire remaining orphan test files.
7. Archive stale dev notes.

**Regression check:**
- `julia --project=docs docs/make.jl` exits clean with no warnings.
- Example scripts run + match committed reference outputs.
- Test suite count records the new total.

**Christian sign-off:** None required (but invite review of README + migration guide).

**Exit criterion:** Clean Documenter build; 5 examples green; test suite expanded; README reviewed by user.

**Estimated scope:** ~2 weeks, scope-gated on the dry-run error count.

---

## Phase 8 — Final merge PR

**Purpose:** Promote `sanghavi-unified` to `main`.

**Preconditions / Gate:** Phases 0–7 exited green.

**Actions:**
1. Full-range re-run of both acceptance scripts on sanghavi-unified; commit outputs alongside reduced-range outputs in `baseline_output/`.
2. User validates full-range outputs match sanghavi reference within tolerance.
3. **Version bump** in `Project.toml` (e.g., 2.0.0 → 2.1.0).
4. Update CHANGELOG with the release entry summarizing the merge.
5. Open `sanghavi-unified → main` PR. CI runs full test suite + Documenter build.
6. User merges when green.
7. **Git tag** at the merge commit.
8. **Release notes** summarizing the merge (can reuse the CHANGELOG entry).
9. Delete `sanghavi-unified` branch post-merge (or keep tagged for history).

**Regression check:**
- Full-range acceptance scripts pass user tolerance vs sanghavi reference.
- Test suite 100% green.
- Documenter 100% clean.

**Christian sign-off:** **NOT required** (per amendments §3 / original §6 Q4). User drives the merge.

**Exit criterion:** `main` pointer on `sanghavi-unified` tip; git tag landed; release notes published.

**Estimated scope:** 2–3 days including CI cycles.

---

## Cross-phase regression invariants

These hold across every phase as invariants:

1. **Sanghavi's GPU Raman optimizations are never rolled back.** After Phase 4, allocation counts never exceed sanghavi's floor.
2. **Elastic path `⊠` notation is not touched.** `doubling.jl`, `interaction.jl`, `elemental.jl` (non-inelastic variants) see only docstring edits in Phase 7, if anything.
3. **Linearized Raman does not appear anywhere.** No new `_lin` files in `Inelastic/`, no Raman branch in `rt_kernel_lin.jl`, no `AddedLayerRS`-style lin constructors.
4. **Baselines are from sanghavi, not unified.** Every tolerance check compares against `baseline_output/sanghavi_9ee9a75/` (or whatever sanghavi tip we freeze).
5. **Every phase exit is user-validated.** No auto-promotion between phases — the assistant pauses after each phase's regression check and surfaces the numbers for user review.
6. **Tolerance claims are honest.** Bit-exactness applies only to phases that preserve execution path verbatim (workspace threading with identical write order, sync removal). Phases that reorder accumulation (staging mode, 4D→3D batching) use tolerance-based verification with user-specified numbers.

---

## Christian coordination — final shape

Per user's 2026-04-19 review session (amendments §3):

- **Phase 0 checkpoint:** Christian signs off on the (edited) Inventory A + Inventory F delta and on `plans/PLAN_AMENDMENTS_2026-04-19.md` before branch `sanghavi-unified` is created. This is the main scope-alignment checkpoint.
- **Phase 4 checkpoint:** Christian reviews the `InteractionWorkspace`-as-adapted-to-`RTModel` design *before* kernel-wiring commits land. Framing: review of additions, not a design review — the base design is sanghavi's by authority.
- **Phase 5 checkpoint (conditional):** If Phase 5 runs (i.e., Phase 4 measurements show headroom), Christian reviews the specific optimization targets and measurement protocol per optimization.
- **Final merge PR (Phase 8):** Does not require Christian's sign-off. User drives.

No new Christian checkpoints beyond these.

---

## Risk log

| Risk | Phase | Mitigation |
|---|---|---|
| Phase 1a re-baselining reveals a test whose numerical shift exceeds what Bodhaine + interpolation should produce, indicating a bug. | 1a | Spot-check ~10 tests with expected-vs-actual analytic comparison before batch-committing; if a shift is unexplained, isolate that test before the batch lands. |
| Christian disagrees with the authority-rule-adapted workspace design at the Phase 4 checkpoint | 4 | Phase 0 Christian-ack covers the authority rule itself; Phase 4 review is scoped to "which additional fields layer on," so disagreement should be bounded. |
| Phase 4 can't match sanghavi's allocation floor | 4 | Root-cause in place; escalate to user before moving to Phase 5. Don't paper over. |
| Phase 5 no-headroom outcome | 5 | Explicitly allowed; commit a dev note instead of a code change. |
| Reduced YAML configurations hide a bug that full-range would reveal | 3 | Phase 8 mandatory full-range re-run. |
| Phase 7 `warnonly=true` dry-run surfaces many more errors than expected | 7 | Pre-flight dry-run scopes the work; if count is prohibitive, scope-trim with user before starting. |
| Orphan tests contain stale code that fails | 1d / 1e / 7 | Triage first; gate wiring on fix. |
| User's hardware differs from baseline-capture hardware over time | 2 | Harness records GPU/CPU identifier in output JSON; re-capture if hardware changes. |
| SIF `1/π` factor sign/unit mismatch | 3 | Port verbatim from sanghavi (two well-tested lines); verify against sanghavi reference in the harness with a pattern check (amendments §2.5). |
| Phase 6 script list balloons beyond expectations | 6 | `phase6_script_port_list.md` pre-flight deliverable; user signs off on list before porting starts. |

---

## Open items — to be resolved during execution

(Items the plan defers to specific phases rather than pre-specifying. Items resolved by amendments §2 are **not** listed here.)

1. **Tolerance numbers** — user fills table above before Phase 2b baseline capture.
2. **Reduced spectral ranges / layer counts** for the two acceptance scripts — determined empirically during Phase 3 (target: <60 s per warm run on reference GPU).
3. **`α̅` verification source citation** — fill in the specific source (notebook, paper section, commit message, email thread) at the Phase 1b port site.
4. **Phase 5 target list** — derived empirically from Phase 4's measured numbers; no pre-specified list.
5. **Stale dev notes archive vs delete** — user preference at Phase 7 kick-off.
6. **Exact `phase6_script_port_list.md` contents** — produced at Phase 6 start.

The following v1 open items are **resolved** by amendments §2 and no longer appear:
- `hem_R` / `hem_T` exposure approach — **deferred post-merge** (amendments §2.3).
- `α̅` unit-convention check — **closed, verified correct** (amendments §2.4).
- `reduce_profile` default — **linear interpolation wins** (amendments §2.1).
- `apply_D_matrix_elemental!` shortcut — **dropped** (amendments §2.2).
- Canopy + Raman coupling — **out of scope, documented** (amendments §2.6).
- SIF `0.5π / max` rescaling — **intentional hack, documented** (amendments §2.5).
- Aerosols module `include` — **moved to Phase 1d** (amendments §2.7).

---

## References

- [plans/PLAN_AMENDMENTS_2026-04-19.md](PLAN_AMENDMENTS_2026-04-19.md) — **authoritative**. Consolidates the 2026-04-19 user review session decisions. This v2 plan mirrors its §1–§7.
- [plans/CLAUDE_HANDOFF_BRIEF.md](CLAUDE_HANDOFF_BRIEF.md) — canonical mission / scope / constraints document. Unchanged.
- [plans/IMPLEMENTATION_PLAN.md](IMPLEMENTATION_PLAN.md) — v1, superseded by this document on points of conflict.
- [plans/inventories/A_sanghavi_optimization_inventory.md](inventories/A_sanghavi_optimization_inventory.md) — optimization commits inventory + workspace comparison. Editorial pass pending per Phase 0.
- [plans/inventories/B_sanghavi_performance_baselines.md](inventories/B_sanghavi_performance_baselines.md) — benchmark infrastructure + harness recommendation.
- [plans/inventories/C_sif_plumbing_status.md](inventories/C_sif_plumbing_status.md) — SIF port gaps.
- [plans/inventories/D_physics_delta_summary.md](inventories/D_physics_delta_summary.md) — non-optimization physics differences. Editorial pass pending per Phase 0.
- [plans/inventories/E_docs_state_audit.md](inventories/E_docs_state_audit.md) — docs/exports/tests audit (drives Phase 7).
- [plans/inventories/F_christian_plan_delta.md](inventories/F_christian_plan_delta.md) — staleness audit of Christian's 2026-03-21 plan. Editorial pass pending per Phase 0.
- [plans/inventories/0_SUMMARY.md](inventories/0_SUMMARY.md) — inventories landing page. Editorial pass pending per Phase 0.
- `docs/dev_notes/sanghavi_unified_merge_plan.md` (Christian, 2026-03-21) — historical; superseded.
- `docs/dev_notes/raman_gpu_optimization.md` (Christian) — historical; allocation figure non-authoritative.
- [CLAUDE.md](../CLAUDE.md) — architecture & terminology reference.

---

## Next actions after user approval

1. Proceed to Phase 0: editorial pass on Inventories A, D, F, 0_SUMMARY per §6 of amendments; forward edited inventories + amendments doc + this v2 plan to Christian.
2. Block on Christian's written ack before creating branch `sanghavi-unified`.
3. On ack, create `sanghavi-unified` off `unified-vsmartmom @ a4e4187` and begin Phase 1a.
