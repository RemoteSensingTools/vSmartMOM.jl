# sanghavi-unified merge — implementation plan

**Authored 2026-04-19. Status: draft, awaiting user approval before any code lands.**

Worktree tips as of this draft:
- sanghavi: `/home/sanghavi/code/github/vSmartMOM.jl/` @ `9ee9a75` (branch `sanghavi`)
- unified-vsmartmom: `/home/sanghavi/code/github/uni_vSmartMOM/` @ `a4e4187` (branch `unified-vsmartmom`)

---

## Context

vSmartMOM.jl has two active branches with ~4 years of divergence since merge base `702cbc3d`:

- `unified-vsmartmom` (250 commits unique) owns the v2.0.0 architecture: `RTModel{ARCH,FT}` hierarchy, `SolverConfig`/`Atmosphere`/`Optics`/`ObsGeometry`/`QuadPoints` sub-structs, Aerosols/IO modules, 474-test suite, Cox-Munk/canopy/RossLi/RPV/Rahman/fresnel surfaces, HITRAN direct download + version tracking, Float32/Float64 flexibility, GPU-accelerated NAI2 Mie, GPU batched-kernel extension, forward Raman physics partially ported from sanghavi.
- `sanghavi` (50 commits unique) owns the physics truth source for Raman + several GPU optimizations landed on top: `InteractionWorkspace` struct with per-direction CPU-staging, `get_n₀_n₁` arithmetic rewrite, Float32 support, Raman benchmark infrastructure, and a set of **forward-physics corrections that did NOT propagate to unified via commit `59f8de8`** (see Inventory D §1 & §3).

The goal is branch `sanghavi-unified` off `unified-vsmartmom` HEAD, port sanghavi's remaining contributions on top (physics fixes + optimizations + SIF + scripts + `perturb_parameters.jl` + single-scatter driver), produce a release-quality doc set, ship a regression harness that pins sanghavi's numbers as the non-regression floor, and eventually promote `sanghavi-unified` to `main`.

**Hard constraints** (from handoff brief §3 and §7):
- Linearized Raman is permanently out of scope — no `AddedLayerRS`-style constructors, no Raman branch in `rt_kernel_lin.jl`.
- Elastic path (`⊠` notation in `doubling.jl`, `interaction.jl`, `elemental.jl`) stays as-is.
- Sanghavi's existing GPU Raman optimizations (`InteractionWorkspace`, `get_n₀_n₁`, per-direction CPU-staging, Raman benchmark infra) must be preserved or improved — never rolled back.
- Baselines come from sanghavi at current tip, not from unified.

**User answers to §6 questions** (captured 2026-04-19):
- Acceptance scripts: `prototype_EMIT_aer_ht.jl` (with `RS_type = noRS()`) and `creategrid_O2Aband_RamanSIF.jl` — with **reduced spectral ranges and (if needed) fewer vertical layers** for dev-iteration speed.
- Physics-fix ordering: land the 083353b fixes on `sanghavi-unified` BEFORE baseline capture.
- Christian coordination: sign-off required on (a) Inventory A+F delta before branch creation, (b) Phase 4 workspace design, (c) Phase 5 further GPU-optimization decisions. NOT required on the final merge PR.
- Research detritus: drop all of it (PNGs, Manifest copies, `.dat~`, hardcoded abs paths, `rt_run_bck.jl`, `inelastic_helper_old.jl`, in-source prototype/plot folders).

---

## Branch strategy

```
unified-vsmartmom @ a4e4187
        │
        └── sanghavi-unified      (created in Phase 0 exit, after inventories approved)
```

No rebasing. No merging sanghavi wholesale. Port specific commits/files/patterns as described in each phase. Commits on `sanghavi-unified` use their own short messages; co-authorship cited to sanghavi where we import physics or optimization code.

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

Per handoff brief §7, tolerances are a user-facing choice per-script, per-quantity. Placeholders to be populated before Phase 2 baselines are committed:

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
| SIF contribution | n/a | rtol=?, atol=? |
| Wall-clock (warm) | regression ≤ +X% of sanghavi | ≤ +X% |
| GPU alloc count | ≤ sanghavi count | ≤ sanghavi count |
| GPU peak memory | ≤ sanghavi peak + Y MB | ≤ sanghavi peak + Y MB |

---

## Phase 0 — Reconciliation (no code)

**Purpose:** Ensure all parties agree on scope before `sanghavi-unified` exists. No branch, no commits.

**Preconditions / Gate:** None. This is the starting phase.

**Actions:**
1. User reads inventories A–F (already produced) at `uni_vSmartMOM/plans/inventories/`.
2. User flags any factual corrections needed; we patch the relevant inventory in place.
3. User forwards Inventory A + Inventory F + this plan to Christian Frankenberg and collects explicit ack on: (a) sanghavi's `InteractionWorkspace` supersedes his `InelasticWorkspace` proposal modulo CPU-staging opt-in; (b) Phase 3 4D→3D batching direction is conditional on Phase 5 re-measurement (not a given); (c) Phase 5 (linearized Raman) permanently dropped; (d) Christian's "~19k allocations" figure is non-authoritative and will be re-measured; (e) unported forward-physics fixes on sanghavi's `083353b` will land on `sanghavi-unified` in Phase 1.
4. Revise this plan if Christian flags anything unexpected.

**Regression check:** n/a (no code).

**Christian sign-off:** YES — blocking. Phase 0 does not exit without Christian's explicit ack on the five items above.

**Exit criterion:** User approves this plan; Christian has acked the A/F delta in writing.

**Estimated scope:** 1–2 days of human calendar; all coordination.

---

## Phase 1 — Physics truth sync on sanghavi-unified

**Purpose:** Port the forward-physics fixes from sanghavi (primarily in commit `083353b`) that did NOT land via `59f8de8`, plus the `get_n₀_n₁` micro-optimization, plus the `rt_run_ss` driver. Do this BEFORE baseline capture so baselines reflect correct physics.

**Preconditions / Gate:** Phase 0 exited. Branch `sanghavi-unified` is created off `unified-vsmartmom @ a4e4187`.

**Files touched** (paths on the sanghavi-unified worktree — adapted from sanghavi-worktree equivalents):

| File | Change | Source on sanghavi |
|---|---|---|
| `src/Inelastic/inelastic_helper.jl` | Fix `γ_mol_Cabannes` / `γ_C_Rayl` inversion. Rename accessor to make it unambiguous. Update `compute_ϖ_Cabannes`. | `vSmartMOM.jl/src/Inelastic/inelastic_helper.jl` via Inv. D §1 |
| `src/Inelastic/inelastic_cross_section.jl` | Drop stray `2π` in `α̅` denominator (pending `ω₀` unit-convention confirmation — **explicit sanity check before merging**). | Inv. D §3 |
| `src/CoreRT/CoreKernel/raman_atmo_prop.jl` | Switch Rayleigh formula to Bodhaine et al. 1999 Eq. 30. Update band-specific γ constants via new `compute_γ_air_Cabannes!` / `compute_γ_air_Rayleigh!` exports. | Inv. A §083353b, Inv. D |
| `src/Inelastic/inelastic_helper.jl` (rename `compute_γ_mol_Rayleigh!` → `compute_γ_mol_Cabannes!`) | Rename + re-wire callers. Ensure no back-compat stub remains after callers migrated. | Inv. A |
| `src/Inelastic/constants` (or wherever `nm_per_m` lives) | Rename `nm_per_m` → `nm_per_cm` (unit-correctness). Grep all call sites. | Inv. A |
| `src/Absorption/*` | LUT range-clamping fix. Locate sanghavi's change, port minimally. | Inv. A |
| `src/Inelastic/inelastic_helper.jl` (or wherever `get_n₀_n₁` lives) | Adopt sanghavi's arithmetic rewrite verbatim (`max(1, 1-Δ):min(nSpec, nSpec-Δ)`). | Inv. A `854b44c` |
| `src/CoreRT/rt_run.jl` | Add `rt_run_ss(model)` driver (currently undefined on unified despite `export`). Adapt from sanghavi `rt_run.jl` lines ~258 & ~450 to `RTModel` field access. | Inv. D §4 |
| `src/CoreRT/CoreKernel/elemental.jl` | `n_stokes == 1` scalar shortcut in `apply_D_matrix_elemental!` to skip GPU kernel dispatch. Inexpensive correctness + micro-perf win. | Inv. A `9a26002`, Inv. D |
| `src/CoreRT/rt_run.jl` | **Decision point: `hem_R` / `hem_T` hemispheric outputs.** Either (a) land the 6-tuple return and update all callers + doc, or (b) expose them via a separate `rt_run_hemispheric` helper keeping `rt_run` at 4-tuple. Default recommendation: (b) — additive helper, no breaking API change. Tag this as **user decision required in Phase 1 kick-off**. | Inv. A `9a26002` |

**Actions:**
1. Create branch: `git checkout -b sanghavi-unified unified-vsmartmom` (NB: the existing VSCode checkout is on `unified-vsmartmom` — branching adds no new worktree).
2. Commit-by-commit port (one logical fix per commit), each with a smoke-test run of the existing Raman `test_forward_noRS.jl` and `test_forward_raman.jl` (the latter must first be wired into `runtests.jl` — see Phase 7) to confirm no regressions beyond the expected physics change.
3. Land the `hem_R`/`hem_T` decision per the user's response in the Phase 1 kick-off.
4. Confirm `α̅` unit-convention before removing `2π`. Check `ω₀` definition in both branches; if unified uses angular frequency and sanghavi uses wavenumber (or vice versa) the `2π` factor may be correct on one side. **Do not blind-port** until this is checked.
5. Commit the `get_n₀_n₁` micro-optimization separately.
6. Commit `rt_run_ss` with an accompanying smoke test (new `test_forward_ss.jl`).

**Regression check:**
- All existing 474 tests in `test/runtests.jl` still pass.
- `test_forward_raman.jl` output now matches sanghavi-branch reference output (up to float precision) for a toy 1-band RRS run. Capture this reference before the branch and commit it as `test/reference/phase1_cabannes_fix.jld2`.
- Sanghavi's `InteractionWorkspace` is NOT yet landed — that's Phase 4. But the physics changes land NOT on sanghavi's optimized code path, but on unified's current Raman path. Confirm the physics fix is computing-strategy-independent (same formula regardless of workspace layout).
- **No sanghavi optimization rolled back:** only unified code is edited in Phase 1. Sanghavi's optimization code isn't present yet; no rollback risk.

**Christian sign-off:** Optional (already covered by Phase 0 ack). Physics fixes are uncontroversial.

**Exit criterion:** All tests pass; new `test_forward_ss.jl` passes; unit-convention sanity check recorded in a code comment above the `α̅` formula; reference JLD2 committed.

**Estimated scope:** Medium — ~6–10 commits, ~1 week human work.

---

## Phase 2 — Baseline capture + regression harness

**Purpose:** Capture the non-regression floor on `sanghavi-unified` (post-Phase-1-physics) AND on the sanghavi branch itself. The comparison target is "`sanghavi-unified` matches sanghavi within user tolerances," but before that's meaningful we need numbers on both.

**Preconditions / Gate:** Phase 1 complete. Tolerance table in this doc populated by user.

**Files touched:**
- `test/benchmarks/harness/run_benchmarks.jl` (new) — wraps each acceptance script, captures wall-clock, CPU allocs, GPU used-delta via `CUDA.Mem.info()` before/after, GPU peak via `CUDA.@profile trace=true`, kernel-launch count, I/Q/U/V and ieR/ieT arrays.
- `test/benchmarks/harness/metrics.jl` (new) — metric computation + dumping.
- `test/benchmarks/harness/report.jl` (new) — loads two baseline JSONs and produces a side-by-side tolerance check report.
- `test/benchmarks/harness/scenarios.toml` (new) — declares which scripts + which reduced YAML configs to run.
- `test/benchmarks/baseline_output/sanghavi_9ee9a75/` (new) — committed JSON baselines from the sanghavi worktree.
- `test/benchmarks/baseline_output/sanghavi-unified_<commit>/` (new) — committed JSON baselines on the merged branch.

**Actions:**
1. User fills the tolerance table above.
2. Author the harness skeleton on `sanghavi-unified`. Keep it minimal — a ~200 line `run_benchmarks.jl` that invokes each acceptance script via `include` under a captured-output wrapper.
3. Run harness on sanghavi-worktree (`cd /home/sanghavi/code/github/vSmartMOM.jl && julia --project=... test/benchmarks/harness/run_benchmarks.jl`). Commit outputs under `baseline_output/sanghavi_9ee9a75/`.
4. Run harness on sanghavi-unified (still no workspace port yet; Raman code runs via unified's existing path with Phase 1 physics fixes applied). Commit outputs under `baseline_output/sanghavi-unified_<commit>/`.
5. Author `report.jl` that compares two baseline directories and emits a tolerance-based pass/fail.
6. Record expected delta narrative: at this point sanghavi-unified's wall-clock and alloc counts will be **worse** than sanghavi's — that's the gap Phase 4 closes.

**Regression check:**
- Harness invocation is idempotent and deterministic (warmed, ≥3 runs, median reported).
- Physics correctness: `report.jl` confirms sanghavi-unified's output matches sanghavi's output within user tolerance on I/Q/U/V/ieR/ieT for both acceptance scripts. **This gate must pass before Phase 3.**
- **No sanghavi optimization rolled back:** nothing yet — sanghavi's optimizations aren't on the branch.

**Christian sign-off:** None needed.

**Exit criterion:** Two baseline directories committed; `report.jl` green on physics correctness; user reviews the wall-clock / alloc gap numbers and confirms the Phase 4 target.

**Estimated scope:** Medium — 1–2 days to build harness, plus re-run time on both branches (runtime dominated by the full-range RRS+SIF script; reduced-range YAML keeps it sane).

---

## Phase 3 — SIF port + single-scatter driver + acceptance-script port

**Purpose:** Land the 2-line SIF injection, commit SIF data files, port the two acceptance scripts with reduced YAML configs, finalize `rt_run_ss` wiring.

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
| `src/Testing/perturb_parameters.jl` | Port from sanghavi for FD-Jacobian verification. | Brief §4.3 |

**Actions:**
1. Port the 2-line Lambertian SIF injection. Small commit.
2. Commit the four real SIF data files. Separate commit. Drop `.dat~` backups in the same commit.
3. Add `sif_loader.jl` helper and route both benchmark scripts through it.
4. Port the two acceptance scripts. Determine reduced spectral ranges and layer counts empirically: pick ranges such that **one full warm run completes in <60 s on the reference GPU** (user-specified hardware; record in harness README). Commit reduced YAMLs alongside full-range YAMLs.
5. Port `perturb_parameters.jl` and add one smoke test.
6. Re-run harness on sanghavi-unified. Update `baseline_output/sanghavi-unified_<commit>/`.

**Regression check:**
- Acceptance scripts now complete on sanghavi-unified.
- SIF injection matches sanghavi's output when `SIF₀` is non-zero (add a dedicated `test_sif_injection.jl` with a minimal configuration).
- `rt_run_ss` smoke test still passes.
- Wall-clock / alloc delta vs sanghavi captured but not yet closed.
- **No sanghavi optimization rolled back:** still nothing to roll back.

**Christian sign-off:** None needed.

**Exit criterion:** Both acceptance scripts run on sanghavi-unified to completion in <60 s each (reduced YAML) and produce outputs within user tolerance vs sanghavi reference.

**Estimated scope:** Medium — ~1 week.

---

## Phase 4 — Workspace landing (`InelasticWorkspace`)

**Purpose:** Land sanghavi's `InteractionWorkspace` (plus doubling/pivot extensions from Christian's `InelasticWorkspace` design) threaded through `rt_run.jl`, `elemental_inelastic.jl`, `doubling_inelastic.jl`, `interaction_inelastic.jl`, in a way that slots into unified's `RTModel` architecture.

**Preconditions / Gate:**
- Phase 3 complete and green.
- **Christian signs off on the workspace design** (struct field list, `staged::Bool` opt-in, doubling-side fields added per his proposal, `ScatteringInterface_XX` coverage decisions).

**Design (to be ratified with Christian before implementation):**
1. Field set = (sanghavi's `InteractionWorkspace` 15 fields) ∪ (Christian's `InelasticWorkspace` doubling/pivot fields) — exhaustive superset.
2. CPU-staging from sanghavi's `d75dacb`: add `staged::Bool` opt-in flag, 3-buffer GPU alias set, 6 CPU staging arrays. Default `false` (non-staged) on GPU; `true` is a memory-saving mode the user opts into when running large `creategrid_*RamanSIF.jl` configurations and is willing to accept ~4% wall-clock cost for ~3.5 GB FP64 memory savings.
3. Allocation site = `rt_run` (not `model_from_parameters`), consistent with sanghavi. Rationale: workspace lifetime matches a single forward RT solve; `rt_run` knows the problem size.
4. Scattering-interface coverage = `ScatteringInterface_00`, `ScatteringInterface_11`, `ScatteringInterface_10`, `ScatteringInterface_01` — extend beyond sanghavi's `_11`-only coverage in this phase to avoid a follow-up.
5. Thread through all Raman kernel call sites. Elastic path untouched (per hard constraint).

**Files touched:**
- `src/CoreRT/types.jl` — `InelasticWorkspace` struct definition.
- `src/CoreRT/rt_run.jl` — allocate + pass through.
- `src/CoreRT/CoreKernel/interaction_inelastic.jl` — consume workspace in place of `similar()` allocations at lines ~293–305 (matching sanghavi's `interaction_inelastic.jl:16-36` design).
- `src/CoreRT/CoreKernel/elemental_inelastic.jl` — same.
- `src/CoreRT/CoreKernel/doubling_inelastic.jl` — same; adds the doubling-tmp fields from Christian's design that sanghavi doesn't have yet.

**Actions:**
1. Commit workspace struct + allocator (no call-site changes yet — should compile + tests should still pass).
2. Commit per-kernel wiring incrementally (one kernel file per commit). After each, re-run harness and confirm physics output matches Phase 3 tolerance.
3. Commit the `staged::Bool` opt-in path + dedicated test exercising it (`test_workspace_staged.jl`).
4. Re-run harness on sanghavi-unified, update `baseline_output/sanghavi-unified_<commit>/`.

**Regression check:**
- **Sanghavi's `InteractionWorkspace` capabilities are preserved:** `ScatteringInterface_11` RRS path shows ≤ sanghavi's allocation count. If it shows MORE allocations than sanghavi, Phase 4 is not done — root-cause and fix before moving on. This is the single most important regression gate in the entire plan.
- CPU-staging opt-in reproduces sanghavi's ~3.5 GB savings on a large FP64 config. Record measured savings.
- Doubling-tmp fields not present on sanghavi — new ground. Confirm no numerical regression vs Phase 3 output (tolerance-based; Phase 3 physics is reference).
- Wall-clock ≤ sanghavi's non-staged wall-clock within user tolerance. Staged mode: ≤ sanghavi's staged wall-clock within tolerance.

**Christian sign-off:** YES — blocking. Review the struct design before Action 1.

**Exit criterion:** Harness green on physics + allocation counts ≤ sanghavi's floor + CPU-staging opt-in works.

**Estimated scope:** Large — 2–3 weeks. This is the workhorse phase.

---

## Phase 5 — Further GPU headroom (conditional)

**Purpose:** Extend beyond sanghavi's optimization set IF Phase 4 numbers show clear headroom. Candidates identified by inventories:
- Doubling-side allocations (`tmp3–tmp6`, `gp_refl`, `ieJ₁⁺/⁻`) not yet absorbed into workspace.
- `batched_mul!` vs `⊠` micro-benchmark re-run (sanghavi's `e89ec1c` showed 5.5× slowdown for 15×15 GEMM — **do NOT go the batched_mul direction on that shape**). But there may be shapes or call sites where `batched_mul!` wins.
- Sync-point reduction on GPU (non-staged path).
- 4D→3D flattening + gather-based batched_mul — **only if re-benchmarking proves it's faster for the actual problem shapes**.

**Preconditions / Gate:**
- Phase 4 complete and green.
- **Christian signs off on the specific optimization(s) to attempt** and the measurement protocol.
- Phase 4 measured numbers show there's actual headroom above sanghavi's floor.

**Actions:**
1. Re-benchmark `batched_mul!` vs `⊠` on CURRENT Phase-4 shapes (not sanghavi's 15×15 specifically; the merged branch's `creategrid_O2Aband_RamanSIF.jl` may hit different shapes). Commit the benchmark + its output.
2. For each candidate that shows empirical win: land it as a separate commit with its own micro-benchmark entry and harness re-run.
3. If no candidate shows a win, Phase 5 exits as a no-op with a committed `dev_notes/phase5_headroom.md` explaining why further optimization isn't justified. This is a valid Phase 5 outcome.

**Regression check:**
- Each optimization commit: allocation count ≤ pre-commit count; wall-clock ≤ pre-commit wall-clock; physics output within user tolerance.
- **Sanghavi-floor check:** running allocation total must remain ≤ sanghavi's floor, regardless of what Phase 5 attempts.

**Christian sign-off:** YES — blocking, per-optimization.

**Exit criterion:** Either (a) one or more commits landed with measured wins and no regressions, or (b) a no-op with a dev note explaining why.

**Estimated scope:** Variable — could be zero (no headroom) or 1–2 weeks (multiple wins).

---

## Phase 6 — Remaining scripts port + detritus cleanup

**Purpose:** Port the rest of sanghavi's benchmark scripts (OCO, CarbonI, O3 Huggins, etc.), port EMIT/Balsamic application scripts, drop all research detritus.

**Preconditions / Gate:** Phase 5 exited.

**Files touched:**
- `test/benchmarks/` — port remaining Raman + SIF scripts from sanghavi (from Inventory B's full inventory: RamanSIFspectra.jl, RamanSIFmaps.jl, RamanSIFworldmaps.jl, O3Huggins_polRaman.jl, prototype_CarbonI.jl, compLelli.jl, compLelli2.jl, RT_vs_hydrolight.jl, 6SV1_R_trues.jl, etc.). **Paths rewritten** to `joinpath(pkgdir(vSmartMOM), ...)`. Scripts NOT on the acceptance list still must compile + run on reduced YAMLs.
- `src/` — remove `rt_run_bck.jl`, `inelastic_helper_old.jl`, any `*_old.jl` / `*_bck.jl` files.
- Top-level `.png` files on sanghavi — not ported.
- `Manifest copy*.toml` files — not ported.
- `.dat~` backup files — not ported (already dropped in Phase 3 for SIF; cover remainder here).
- In-source prototype/plot folders — not ported.
- `.gitignore` — confirm excludes `*.jld`, `*.jld2`, top-level `*.png`, `Manifest copy*.toml`.

**Actions:**
1. Enumerate remaining sanghavi benchmark scripts to port (derive from Inventory B).
2. Port each with absolute-path rewrites. Commit scripts in logical groupings (SIF variants together, retrieval comparisons together, etc.).
3. Ensure every ported script has a reduced-YAML entry usable from harness.
4. One cleanup commit removes detritus files not ported.
5. `.gitignore` tightened.

**Regression check:**
- Every ported script runs to completion on reduced YAML without errors.
- Acceptance scripts (Phase 3) still green.
- No sanghavi optimization rolled back (Phase 4's workspace remains intact).

**Christian sign-off:** None needed.

**Exit criterion:** `test/benchmarks/` contains the full set of ported scripts; detritus is gone; `.gitignore` updated.

**Estimated scope:** Medium — 1 week.

---

## Phase 7 — Docs overhaul

**Purpose:** Release-quality documentation: README rewrite, Documenter user guide, curated examples, API reference audit, narrative tutorial, migration guide (old-sanghavi API → new-unified API).

**Preconditions / Gate:** Phase 6 exited (scripts stable enough to document).

**Files touched:**
- `README.md` — rewrite. Cover install, quick start, `RTModel` + forward/linearized RT, surface types (one subsection each), Raman + SIF, GPU selection, Float32, HITRAN-2024 edition switching, link to examples + docs.
- `CHANGELOG.md` — expand v2.0.0 entry with the full break list (not just the `vSmartMOM_Model → RTModel` paragraph): HITRAN-2024, Float32, GPU Mie, Canopy, GEOSChem, batched-kernel, merged Raman from sanghavi, etc.
- `docs/make.jl` — remove `warnonly=true`; fix every missing `@docs` reference that surfaces. Add `@docs` entries for `RTModelLin`, `SolverConfig`, `Optics{,Lin}`, `ParameterLayout` accessors, `FwdMode`/`LinMode`, SolarModel exports.
- `docs/src/` — new pages: (a) migration guide (old-API → new-API pair for every converted script — at minimum the 2 acceptance scripts, with diffs), (b) Raman + SIF narrative tutorial, (c) per-surface quickstart, (d) GPU quickstart.
- `examples/` — 5 new examples covering the golden path:
  1. `lambertian_forward.jl` (fastest warm-up example; committed JSON reference output).
  2. `coxmunk_polarized.jl` (ocean polarization).
  3. `forward_jacobians.jl` (linearized RT elastic).
  4. `rrs_sif.jl` (Raman + SIF; the flagship showcase).
  5. `gpu_cuda.jl` (GPU quick-start).
  Each example: `.jl` + narrative `.md` + committed JSON/CSV reference output. Tests wire each into a quick regression gate.
- `src/vSmartMOM.jl` — fix the `src/Aerosols/Aerosols.jl` inclusion bug (Inv. E finding — currently exports 11 symbols from a module that isn't `include`d).
- `src/Aerosols/aerosol_integration_example.jl` — revive now that the module is loaded.
- Docstring gaps: fill `AbstractPolarizationType`, `NAI2`, `PCW`, Absorption error-function structs, `FwdMode`, `LinMode`, `RTModelLin`, `SolverConfig`, `Optics`, etc.
- `docs/dev_notes/CUDA_SETUP.md` → migrate to `docs/src/advanced/gpu.md`.
- `docs/dev_notes/JACOBIAN_TEST_WORKFLOW.md` → migrate to `docs/src/advanced/jacobians.md`.
- `docs/dev_notes/*_stale.md` (see Inv. E §3) — archive to `docs/dev_notes/archive/` or delete, per user preference.
- `test/runtests.jl` — wire in the 7 orphan test files: `test_Aerosols.jl`, `test_forward_lin.jl`, `test_forward_raman.jl`, `test_hybrid_ad.jl`, `test_jacobians_GPU.jl`, `test_mie_gpu.jl`, `test_performance.jl`. Triage each first; if any are stale, fix before wiring.

**Actions:**
1. Fix Aerosols module `include` bug first — it's a live export-without-definition issue.
2. Turn off `warnonly=true`; fix the cascade of missing-docstring / bad-`@ref` errors (expect ~30–50 to fix).
3. Write the 5 examples with reference outputs.
4. Rewrite README; update CHANGELOG.
5. Author migration guide + narrative tutorial.
6. Triage + wire 7 orphan test files.
7. Archive stale dev notes.

**Regression check:**
- `julia --project=docs docs/make.jl` exits clean with no warnings.
- Example scripts run + match committed reference outputs.
- Test suite count is now 474 + examples' wired tests + orphans' wired tests. Record new total.

**Christian sign-off:** None required (but invite review of README + migration guide).

**Exit criterion:** Clean Documenter build; 5 examples green; test suite expanded; READMEs reviewed by user.

**Estimated scope:** Large — 2 weeks.

---

## Phase 8 — Final merge PR

**Purpose:** Promote `sanghavi-unified` to `main`.

**Preconditions / Gate:** Phases 0–7 exited green.

**Actions:**
1. Full-range re-run of both acceptance scripts on sanghavi-unified; commit outputs alongside reduced-range outputs in `baseline_output/`.
2. User validates full-range outputs match sanghavi reference within tolerance.
3. Update CHANGELOG with the v2.1.0 (or whatever version) entry summarizing the merge.
4. Open `sanghavi-unified → main` PR. CI runs full test suite + Documenter build.
5. User merges when green.
6. Delete `sanghavi-unified` branch post-merge (or keep tagged for history).

**Regression check:**
- Full-range acceptance scripts pass user tolerance vs sanghavi reference.
- Test suite 100% green.
- Documenter 100% clean.

**Christian sign-off:** Per user's §6 answer, NOT required for the final merge PR. User drives the merge.

**Exit criterion:** `main` pointer on `sanghavi-unified` tip.

**Estimated scope:** Small — 2–3 days including CI cycles.

---

## Cross-phase regression guarantees

These hold across every phase as invariants:

1. **Sanghavi's GPU Raman optimizations are never rolled back.** After Phase 4, allocation counts never exceed sanghavi's floor.
2. **Elastic path `⊠` notation is not touched.** `doubling.jl`, `interaction.jl`, `elemental.jl` (non-inelastic variants) see only docstring edits in Phase 7, if anything.
3. **Linearized Raman does not appear anywhere.** No new `_lin` files in `Inelastic/`, no Raman branch in `rt_kernel_lin.jl`, no `AddedLayerRS`-style lin constructors.
4. **Baselines are from sanghavi, not unified.** Every tolerance check compares against `baseline_output/sanghavi_9ee9a75/` (or whatever sanghavi tip we freeze).
5. **Every phase exit is user-validated.** No auto-promotion between phases — the assistant pauses after each phase's regression check and surfaces the numbers for user review.

---

## Risk log

| Risk | Phase | Mitigation |
|---|---|---|
| `α̅` `2π` fix is wrong on unified-side unit convention | 1 | Explicit unit-convention sanity check before landing; code comment documenting the check. |
| Christian disagrees with workspace design during Phase 4 gate | 4 | Phase 0 delta review surfaces this early; Phase 4 kick-off is gated on explicit ack. |
| Phase 4 can't match sanghavi's allocation floor | 4 | Root-cause in place; escalate to user before moving to Phase 5. Don't paper over. |
| Phase 5 no-headroom outcome | 5 | Explicitly allowed; commit a dev note instead of a code change. |
| Reduced YAML configurations hide a bug that full-range would reveal | 3 | Phase 8 mandatory full-range re-run. |
| 7 orphan tests contain stale code that fails | 7 | Triage first; gate wiring on fix. |
| User's hardware differs from baseline-capture hardware over time | 2 | Harness records GPU/CPU identifier in output JSON; re-capture if hardware changes. |
| SIF `1/π` factor sign/unit mismatch | 3 | Port verbatim from sanghavi (two well-tested lines); verify against sanghavi reference in the harness. |

---

## Open items — to be resolved during execution

(These are items the plan defers to specific phases rather than pre-specifying.)

1. **Tolerance numbers** — user fills table above before Phase 2 baseline capture.
2. **Reduced spectral ranges / layer counts** for the two acceptance scripts — determined empirically during Phase 3 (target: <60 s per warm run on reference GPU).
3. **`hem_R` / `hem_T` exposure approach** — default recommendation is (b) additive helper; confirm with user at Phase 1 kick-off.
4. **`α̅` unit-convention check** — resolved during Phase 1 with an explicit code comment.
5. **Phase 5 target list** — derived empirically from Phase 4's measured numbers; no pre-specified list.
6. **Stale dev notes archive vs delete** — user preference at Phase 7.

---

## References

- `plans/CLAUDE_HANDOFF_BRIEF.md` — canonical mission / scope / constraints document (supersedes Christian's 2026-03-21 plan).
- `plans/inventories/A_sanghavi_optimization_inventory.md` — optimization commits inventory + workspace comparison.
- `plans/inventories/B_sanghavi_performance_baselines.md` — benchmark infrastructure + harness recommendation.
- `plans/inventories/C_sif_plumbing_status.md` — SIF port gaps.
- `plans/inventories/D_physics_delta_summary.md` — non-optimization physics differences (drives Phase 1).
- `plans/inventories/E_docs_state_audit.md` — docs/exports/tests audit (drives Phase 7).
- `plans/inventories/F_christian_plan_delta.md` — staleness audit of Christian's 2026-03-21 plan.
- `plans/inventories/0_SUMMARY.md` — inventories landing page.
- `docs/dev_notes/sanghavi_unified_merge_plan.md` (Christian, 2026-03-21) — historical; superseded.
- `docs/dev_notes/raman_gpu_optimization.md` (Christian) — historical; allocation figure non-authoritative.
- `CLAUDE.md` — architecture & terminology reference.
