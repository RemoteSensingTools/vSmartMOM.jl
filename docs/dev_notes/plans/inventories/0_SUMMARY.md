# Inventories summary — sanghavi-unified merge (Stage 1 output)

As-of 2026-04-19. Worktree tips: sanghavi `9ee9a75`, unified-vsmartmom `a4e4187`.

**Editorial-pass note (2026-04-19):** This summary was written before the 2026-04-19 user review session established the **authority rule** — *sanghavi is the authority for the inelastic path*. The recommendations in Findings #3 and in "What this means for the implementation plan shape" have been rewritten to match the authority rule. The *evidence* and *factual observations* stand unchanged. Per-item rationale in the body of each inventory is retained for historical reference; for binding decisions see `plans/PLAN_AMENDMENTS_2026-04-19.md` and `plans/IMPLEMENTATION_PLAN_v2.md`.

Six inventories produced per §5 of `plans/CLAUDE_HANDOFF_BRIEF.md`. Each is a standalone markdown file you should read and validate before the implementation plan is written.

- [A — sanghavi_optimization_inventory.md](A_sanghavi_optimization_inventory.md)
- [B — sanghavi_performance_baselines.md](B_sanghavi_performance_baselines.md)
- [C — sif_plumbing_status.md](C_sif_plumbing_status.md)
- [D — physics_delta_summary.md](D_physics_delta_summary.md)
- [E — docs_state_audit.md](E_docs_state_audit.md)
- [F — christian_plan_delta.md](F_christian_plan_delta.md)

---

## Top findings that reshape the merge plan

### 1. Forward Raman physics is NOT fully ported to unified

The handoff brief (§4.1) claims "Forward Raman physics ported from sanghavi in commit `59f8de8`." **Inventories A and D both flag this as outdated.** Sanghavi's later commit `083353b` ("Update lin+Raman") contains forward-physics corrections that are missing from unified:

- **`ϖ_Cabannes` / `γ_mol_Cabannes` bug fix** (D §1, P0). Unified still treats `mol.effCoeff.γ_C_Rayl` as γ_mol_Cabannes when it actually stores γ_mol_Rayleigh — a long-standing inversion that corrupts all RRS and VS output. Sanghavi renamed the accessor to make the fix unambiguous.
- **Rayleigh formula switched to Bodhaine et al. 1999 Eq. 30** (A).
- **`α̅` frequency correction**: sanghavi dropped a stray `2π` factor in the `(1-(c·ν_eff/ω₀)²)` denominator (D §3, P0 pending unit-convention sanity check on `ω₀`).
- **Band-specific effective-T N₂/O₂ constants** via new `compute_γ_air_Cabannes!` / `compute_γ_air_Rayleigh!` exports (A).
- **Constant rename** `nm_per_m` → `nm_per_cm` — unit-correctness fix (A).
- **Absorption LUT range-clamping** (A).

**Implication:** Baselines captured on sanghavi would disagree with unified-as-is even before any optimization work. These physics fixes must land on `sanghavi-unified` before regression comparisons are meaningful.

### 2. `rt_run_ss` is advertised on unified but undefined

Inventory D §4: unified exports `rt_run_ss` and ships `rt_kernel_ss.jl` + `interaction_ss.jl`, but the driver function `rt_run_ss(model)` has no definition — calling it raises `UndefVarError`. The single-scatter approximation kernel mentioned in the handoff brief §4.3 must be ported from sanghavi's `rt_run.jl` (~lines 258 & 450), adapted from `vSmartMOM_Model` field access to `RTModel` sub-struct access.

### 3. Sanghavi's `InteractionWorkspace` is narrower than Christian's proposal

Inventory A: Christian's proposed `InelasticWorkspace` is a strict superset of sanghavi's `InteractionWorkspace` plus extra doubling/pivot fields, EXCEPT sanghavi's `d75dacb` adds a CPU-staging mode (3-buffer GPU aliasing + 6 CPU staging arrays + `staged::Bool` flag) that saves ~3.5 GB FP64 at ~4.2% runtime cost. Sanghavi covers `ScatteringInterface_11` RRS only.

**Recommendation (per authority rule, 2026-04-19):** Build the merged workspace to **sanghavi's `InteractionWorkspace` layout adapted to `RTModel`**. Christian's proposed additional fields (doubling-tmp, `batch_inv!` pivot/info, shared 3D buffers with elastic `RTWorkspace`) layer on **if compatible**; if not compatible, sanghavi wins. `staged::Bool = true` **default-on** (not an opt-in). This is a review-of-additions, not a design-review, at the Phase 4 Christian checkpoint.

### 4. Christian's "bit-exact after Phase 2/3" claim is false by construction

Inventory F: Phase 3's 4D→3D + gathered `batched_mul` reorders FP accumulation relative to nSpec×nRaman sequential 2D `mul!`. Tolerance-based verification is required, not bit-exact. Also: sanghavi's own micro-benchmark (`e89ec1c`) shows **`NNlib.batched_mul!` is 5.5× slower than `⊠` for 15×15 GEMM** — so the naive "replace loops with batched_mul" direction in Christian's plan is actively anti-performance; Phase 3 needs revalidation before being put in the plan.

### 5. Missed 4th-major commit on sanghavi

Inventory A flags commit `9a26002` (not in the original 5-commit list): adds `hem_R` / `hem_T` hemispheric-integrated radiance outputs, **changing `rt_run`'s return signature from 4-tuple to 6-tuple**. This is a breaking API change that must get its own line item in the implementation plan.

### 6. SIF plumbing at parity by types, one 2-line gap in kernel

Inventory C: `SIF₀` field exists on all RS_type variants on both branches; the only missing code is a 2-line `(1/π) * SIF₀` injection into `J₀⁻` in `lambertian_surface.jl`. Real SIF data files (14 KB `sif-spectra.csv`, 8 KB `PC1_SIFSpectra_allSpecies.csv`, 35 KB `ficus_refl.dat`, 3 KB `ficus_refl_600to800nm.dat`) exist untracked on sanghavi; .dat~ backups are what's tracked. Commit the four real files, drop the backups, port the two `creategrid_*RamanSIF.jl` benchmarks with absolute paths rewritten to package-relative.

### 7. Docs state has concrete blocking gaps

Inventory E: `docs/make.jl` builds with `warnonly=true`, silencing missing docstrings. API reference is missing `@docs` entries for `RTModelLin`, `SolverConfig`, `Optics{,Lin}`, `ParameterLayout` accessors, `FwdMode`/`LinMode`, most SolarModel exports. 7 test files exist on disk but aren't wired into `runtests.jl`. README doesn't teach `RTModel`, Jacobians, surfaces, Float32, GPU selection, or HITRAN-2024 switching. **Live bug:** `src/Aerosols/Aerosols.jl` is not `include`d in `src/vSmartMOM.jl`, so `using vSmartMOM.Aerosols` and `aerosol_integration_example.jl` are dead today.

### 8. Benchmark infrastructure partially exists on sanghavi

Inventory B: sanghavi already has `raman_optimization_baseline.jl` + `raman_optimization_compare.jl` + `raman_reference.jld2` as a thin regression harness, but it covers only a 1-band opttest, captures wall-clock only (no GPU memory / no allocation counts), and `timing_baseline.txt` (RRS=1671.9 s / noRS=2.8 s on `NquadN=15, nSpec=5424, nRaman=172`) predates the `9ee9a75` tip and must be regenerated. Unified has additional Raman benchmarks sanghavi doesn't have (`raman_large_scale.jl`, `raman_memory_benchmark.jl`) that should land on the merged branch too. Christian's "~19k allocations" figure is stale and non-authoritative.

---

## What this means for the implementation plan shape

1. **Phase 0 (reconciliation, no code).** User validates inventories; Christian sign-off on the A/F delta before any branch creation.
2. **Phase 1 (physics truth sync + inelastic port).** Under the authority rule: bulk-port `src/Inelastic/` + inelastic `*_inelastic.jl` kernels from sanghavi. Sub-phases per v2 plan: 1a (Bodhaine + `reduce_profile` default switches, re-baseline 474 tests), 1b (inelastic port — drops `hem_R`/`hem_T` and D-matrix scalar shortcut per amendments §2.2/§2.3; adds α̅ verification comment per §2.4; adds canopy+Raman TODO per §2.6), 1c (`rt_run_ss` driver port), 1d (Aerosols module wire-in per §2.7), 1e (`perturb_parameters.jl` port). This must come BEFORE baseline capture or baselines compare against physics sanghavi has since corrected.
3. **Phase 2 (baseline capture).** Regenerate Inventory B's timing numbers and allocation counts on **both branches at current tips**, on a thin extended harness that also captures GPU memory and allocation deltas. Commit JSON/CSV baselines.
4. **Phase 3 (SIF + single-scatter port).** 2-line Lambertian SIF injection + 4 data files committed + two `creategrid_*RamanSIF.jl` benchmarks ported with path rewrites. Single-scatter driver adapted from sanghavi to RTModel access.
5. **Phase 4 (workspace landing).** Build the merged `InteractionWorkspace` by adapting sanghavi's layout to `RTModel` per the authority rule; Christian's additions (doubling-tmp, `batch_inv!` pivot/info) layer on if compatible. `staged::Bool = true` default-on. Thread through elemental / doubling / interaction. Validate against Phase 2 baselines.
6. **Phase 5 (further GPU headroom).** ONLY if Phase 4 shows headroom and Phase 2 empirically proves `batched_mul` isn't the right tool. Doubling-side allocations (`tmp3-tmp6`, `gp_refl`, `ieJ₁⁺/⁻`) remain open even post-Phase-4.
7. **Phase 6 (scripts + detritus cleanup).** Port EMIT/Balsamic + perturb_parameters.jl + remaining benchmarks; drop PNGs, Manifest copies, `.dat~`, hardcoded paths.
8. **Phase 7 (docs overhaul).** Fix `warnonly=true`, fill docstring gaps, wire 7 orphan tests, migrate CUDA_SETUP + JACOBIAN_TEST_WORKFLOW to Documenter, rewrite README, add committed-output examples covering the golden path, `include` the Aerosols module.
9. **Phase 8 (final merge PR).** Christian reviews, tag `sanghavi-unified` as new `main`.

Each phase gets explicit regression checks confirming no sanghavi optimization is rolled back and a gate requiring Phase N-1 to pass before Phase N starts.

---

## What still needs user input (§6 open questions)

**Resolution status (2026-04-19 user review session + amendments):**

1. ~~Acceptance script list~~ — **Resolved.** `prototype_EMIT_aer_ht.jl` (noRS) and `creategrid_O2Aband_RamanSIF.jl` (RRS + SIF) are the acceptance pair for this merge. Full ported-script list for Phase 6 is a separate pre-flight deliverable (`plans/phase6_script_port_list.md`, user-signed) per v2 plan.
2. **Per-script spectral truncation** — **Still open.** Determined empirically during Phase 3 (target: <60 s per warm run on reference GPU).
3. **Per-quantity tolerance bars** — **Still open.** User fills the tolerance table in `plans/IMPLEMENTATION_PLAN_v2.md` before Phase 2b baseline capture.
4. ~~Christian coordination~~ — **Resolved (amendments §3):** Phase 0 sign-off on Inventory A/F delta + amendments doc; Phase 4 review of workspace-design-adapted-to-`RTModel` *additions* (not a design review — base is sanghavi's by authority); Phase 5 per-optimization reviews if conditional phase runs; NOT required on the final merge PR.
5. ~~Cleanup confirmation~~ — **Resolved.** Drop all research detritus (PNGs, Manifest copies, `.dat~`, hardcoded absolute paths, `rt_run_bck.jl`, `inelastic_helper_old.jl`, in-source prototype/plot folders).
