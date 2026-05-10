# VLIDORT golden-standard validation suite — design

**Document version:** v0.4
**Status:** focused implementation plan; standalone, can be built before the bigger refactor.

**Document set (delivered together):**
- `vsmartmom_dispatch_design_v0_6.md` — architecture this baseline must protect.
- `standalone_ss_solver_plan.md` — Piece A; the kernel-based exact-SS solver.
- `vlidort_baseline_suite.md` — *this document* (Piece B).
- `dev_notes/exact_ss_reference/` — standalone four-path reference.

**Pre-existing committed artifacts referenced (in PyVLIDORT-main distribution):**
- `vlidort_v_test/V2p8p3_Siewert2000_validation.f90` + `.cfg` — Siewert 2000 Problem IIA setup.
- `vlidort_v_test/saved_results/gfortran/results_Siewert2000_validation.all` — peer-reviewed reference (Tables 2-4 from Siewert 2000).
- `vlidort_s_test/2p8p3_solar_tester.f90` + `.cfg` — multi-task scalar setup.
- `vlidort_s_test/saved_results/gfortran/results_solar_tester.all` — multi-task scalar reference.
- `vlidort_v_test/V2p8p3_solar_tester.f90` — multi-task vector setup.
- `vlidort_v_test/saved_results/gfortran/nstokes3/results_solar_tester_IQU0.all` — multi-task vector reference.

### Changes from v0.3
- **§3 (validation scopes): explicit alignment with Piece A's scope-bounded back-correction.** Tasks 1↔3 in `solar_tester` validate **FO-equivalent back-correction** (Piece A §7); Tasks 1↔4 (with full sphericity) require comparison via `ExactSFIPhase` (architecture doc §6). Two distinct validation tracks per case.
- **§4 (Stage 1 cases): refined task selection.** The most useful comparison for our SS work is Task 1 vs Task 3 (FO Regular PS difference) for the FO-equivalent comparison; Task 1 vs Task 4 (Enhanced PS) requires sphericity treatment we don't have in plane-parallel vSmartMOM.

### Changes carried forward from v0.3
- **§3** uses VLIDORT-shipped gold-standard fixtures (Siewert 2000, solar_tester). No PyVLIDORT installation needed initially.
- **§4** Stage 1 leads with the fastest-time-to-data approach.

---

## 0. What this is

A reproducible test suite that compares vSmartMOM RT outputs against gold-standard reference data shipped with the VLIDORT distribution. Records per-case agreement and flags divergence.

Built **before** the bigger refactor. Acts as both:

1. **Regression baseline** through the refactor.
2. **Existing-bug discovery** in current `sanghavi-unified`.

**Key insight**: VLIDORT distribution ships gold-standard fixtures already. We don't have to generate them via PyVLIDORT — the committed `saved_results/` tree contains both peer-reviewed benchmark reproductions and VLIDORT-internal multi-configuration comparisons. We parse these and run vSmartMOM with matching configurations.

---

## 1. The shipped gold standards

The PyVLIDORT distribution at `src/Components/rtms/RTSI/VLIDORT2p8p3/` includes two categories of pre-computed validation data:

### 1.1 Category A — Peer-reviewed benchmark reproductions

**Siewert (2000) Problem IIA** — vector RT validation problem from a peer-reviewed paper:
- Single-layer aerosol slab, τ_total = 1.0
- Lambertian albedo = 0.0 (black surface)
- ω = 0.973527, plane-parallel
- μ₀ = 0.6 (cos SZA), 11 viewing angles, 3 azimuths = 33 geometries
- 7 τ levels for vertical output
- **Polarized: I, Q, U components in Tables 2, 3, 4**
- Greek matrix coefficients embedded in Fortran source as DATA statements (`PROBLEM_IIA(6,0:11)`)

This is **peer-reviewed truth**. Siewert's paper publishes the numbers; VLIDORT reproduces them; if vSmartMOM also reproduces them, vSmartMOM is correct on this case (modulo polarization sign conventions).

The Greek coefficients are committed *in the Fortran source*, so we extract them into Julia and reproduce the configuration in vSmartMOM without running PyVLIDORT.

### 1.2 Category B — Multi-task self-consistency benchmarks

`solar_tester` runs the same atmospheric setup through 6 solver configurations:

| Task | FO Correction | δ-M Scaling | Sphericity |
|---|---|---|---|
| 1 | none | none | plane-parallel |
| 2 | none | yes | plane-parallel |
| 3 | regular PS | yes | sun curved |
| 4 | enhanced PS | yes | sun + LOS curved |
| 5 | as 4 + Solution Saving | yes | sun + LOS curved |
| 6 | as 5 + BVP Telescoping | yes | sun + LOS curved |

The output file records I (and Q, U for vector) for each task at multiple geometries × levels.

**Why this is so valuable for us:**

- **Task 1**: closest to "vSmartMOM-no-correction-no-DM, plane-parallel" — direct apples-to-apples
- **Task 2**: "vSmartMOM with δ-M only, no SS correction"
- **Task 3**: Task 2 + Regular PS FO correction (sphericity from sun side; LOS still planar)
- **Task 4**: Task 3 + Enhanced PS (also LOS sphericity)

The **Task 1 → Task 3 difference** is closest to "what FO correction adds in plane-parallel-like conditions" — and this is what the standalone solver's back-correction (Piece A §7) is designed to reproduce. (Task 3 still has sun-side sphericity, so the comparison isn't perfect; see §3.2 below.)

### 1.3 Category C — Other shipped artifacts

- Stokes-3 lattice geometries (`saved_results/.../nstokes3/`)
- Doublet-geometry post-processing
- Observation-geometry post-processing
- BRDF self-tests and BRDF+ Jacobian tests
- `ProblemIII.Moms` — Mie F-matrix expansion for an additional published problem

These extend our reference set as the suite expands.

---

## 2. Why this comes first

Three reasons:

1. **Bugs before changes.** Q3 (missing upward attenuation in Cox-Munk SS correction) is a candidate; surfacing existing bugs before the refactor is cleaner than after.

2. **Regression detection during refactor.** Each major step gets validated; regressions caught immediately.

3. **Confidence at the end.** Same baseline produces same numbers post-refactor → strongest "we didn't break anything" claim.

The suite is **independent of the standalone solver** (Piece A); both proceed in parallel.

---

## 3. Validation scopes

This is the most subtle part of the suite design. Different validation comparisons have different scopes; conflating them produces nonsensical conclusions.

### 3.1 Scope mapping

| What is being validated | Reference data | vSmartMOM configuration | Notes |
|---|---|---|---|
| Plane-parallel truncated MS | Task 1 in `solar_tester` (no FO, no DM, plane-parallel) | `FullMOM(sfi_phase=TruncatedSFIPhase, ss_paths=:none)` with no δ-M | Cleanest comparison; both codes do the same thing |
| Plane-parallel + δ-M | Task 2 | `FullMOM(sfi_phase=TruncatedSFIPhase, ss_paths=:none)` with δ-M | Tests δ-M handling |
| Standalone solver paths 1+2 (FO-equivalent) | Task 1 → Task 3 difference | `run_exact_ss(...; paths=:paths_1_2)` standalone | The first-scatter contribution; sphericity differences caveat below |
| Back-corrected vSmartMOM (FO-equivalent) | Task 3 | `FullMOM(...) + apply_back_correction!` | Validates the back-correction adapter (Piece A §7) |
| `ExactSFIPhase` full-MOM correction | Custom run via PyVLIDORT (Stage 2) | `FullMOM(sfi_phase=ExactSFIPhase)` | Higher-order paths corrected; not bit-equivalent to FO; Task 4 has sphericity |
| Peer-reviewed polarized reference | `results_Siewert2000_validation.all` Tables 2-4 | `FullMOM(sfi_phase=ExactSFIPhase)` for Siewert config | Plane-parallel + black surface + single layer; clean test |

### 3.2 Sphericity caveat

VLIDORT's `solar_tester` Tasks 3-6 use sphericity (Regular PS or Enhanced PS); vSmartMOM is plane-parallel. So:

- **Task 1 vs Task 2 difference** is δ-M effect — directly comparable to vSmartMOM δ-M-on minus δ-M-off in plane-parallel mode.
- **Task 1 vs Task 3 difference** is *δ-M + sun-side sphericity + FO-correction*. To isolate the FO part for direct comparison with our standalone solver, we'd need either (a) a custom plane-parallel-no-sphericity-with-FO VLIDORT run via PyVLIDORT, or (b) accept that Task 1 vs Task 3 is approximate when sphericity effects are small (they are, for low SZA and moderate τ).

For Stage 1, we use Task 1 vs Task 3 as the FO-equivalent reference, with the sphericity caveat documented. Tasks 5, 6 don't add new physics relevant to our SS work.

For Stage 2, we add custom plane-parallel runs via PyVLIDORT to get exact comparisons.

### 3.3 Polarization caveat

VLIDORT and vSmartMOM have known convention differences for U/V components. Initial validation focuses on scalar I; vector I in cases where conventions are confirmed to match (the Siewert 2000 case is documented well enough that we can verify); per-case tolerance notes flag suspected convention mismatches.

---

## 4. Stage 1 — Parse shipped fixtures, validate vSmartMOM

### 4.1 Stage 1 cases

Cases come from the shipped `saved_results/` data:

**Case A — Siewert 2000 Problem IIA (polarized vector RT):**
- Source: `V2p8p3_Siewert2000_validation.f90` + `results_Siewert2000_validation.all`
- vSmartMOM configuration: extract Greek coefficients from Fortran DATA statements; build single-layer model
- Reference: peer-reviewed; published in Siewert 2000; reproduced in VLIDORT
- Compare: I, Q, U at 33 geometries × 7 τ levels

**Cases B1-B6 — Solar tester multi-task (scalar):**
- Source: `2p8p3_solar_tester.f90` + `results_solar_tester.all`
- vSmartMOM configuration: 23-layer aerosol atmosphere from `input_atmos.dat`, Lambertian albedo=0.05
- For each Task 1, 2, 3 (skipping 4-6 because sphericity comparison is harder), compare vSmartMOM-with-matching-config to VLIDORT output

**Cases C1-C3 — Solar tester multi-task (vector, Stokes-3):**
- Source: `V2p8p3_solar_tester.f90` + `saved_results/.../nstokes3/results_solar_tester_IQU0.all`
- Same atmosphere as B but with full I, Q, U output

**Case D — FO-equivalent back-correction validation:**
- Compute vSmartMOM-Task1-config + back-correction
- Compare to Task 3 (with sphericity caveat from §3.2)
- Validates Piece A §7's back-correction adapter

Total Stage 1: 1 + 3 + 3 + 1 = **8 reference comparison points** ready from shipped data.

### 4.2 Stage 1 implementation phases

**Step 1 (~1 day): Parsers for VLIDORT result files**
- `results_Siewert2000_validation.all`: tabular fixed-width with `Table 2/3/4` headers (I, Q, U)
- `results_solar_tester.all`: tabular fixed-width with `(geometry, direction, level, tasks 1-6)`
- ~100 lines of Julia parsing code total

**Step 2 (~1 day): Configuration extractors**
- Read Fortran DATA statements from `V2p8p3_Siewert2000_validation.f90` for Greek coefficients
- Read `2p8p3_solar_tester.f90` plus `2p8p3_VLIDORT_ReadInput.cfg` for atmospheric profile
- Parse `input_atmos.dat` for layer-by-layer profile data
- One-time effort; afterward configurations are committed as Julia files

**Step 3 (~1 day): vSmartMOM model construction + comparison harness**
- Build `RTModel` matching each VLIDORT case
- Map Greek matrix indices to vSmartMOM's `Z⁺⁺/Z⁻⁺` representation (or Greek-coefficient pathway, depending on contributor type)
- Generic comparison utility

**Step 4 (~1 day): Run, report, triage**
- Execute all 8 cases
- Produce report with per-case verdicts
- Categorize failures: pre-refactor fix / tolerated / deferred
- Q3 verdict from cases where Cox-Munk-like geometry is exercised (Stage 2 will have direct Cox-Munk cases; Stage 1 gives early signal)

**Total: ~3-4 days** focused work.

### 4.3 Stage 1 deliverables

- Parsers for `results_*.all` files
- Configuration extractors (committed as Julia files)
- 8 case configurations in `test/vlidort_baseline/stage1/`
- Comparison harness
- Stage 1 results report:
  - Siewert 2000: vSmartMOM vs published I, Q, U
  - Solar tester Tasks 1, 2, 3 (and Stokes-3 versions)
  - Back-correction validation (Case D)
- Triage list

### 4.4 Stage 1 success criteria

- All 8 cases run end-to-end
- Siewert 2000 polarized I passes at `rtol = 1e-4`
- Solar tester Task 1 passes at `rtol = 1e-5` (cleanest comparison)
- Solar tester Task 3 (with sphericity caveat) passes at `rtol = 1e-3` (sphericity differences contribute)
- Back-correction validation (Case D): back-corrected vSmartMOM vs Task 3 — quantifies FO-equivalent reproduction quality

If Task 1 passes at 1e-5, then any Task 3 disagreement quantifies "FO + sphericity + δ-M". If we have an estimate of sphericity contribution (small for the test atmosphere), the residual quantifies the FO-equivalent error and tells us how the back-correction is doing.

---

## 5. Stage 2 — Custom PyVLIDORT runs (deferred)

After Stage 1 demonstrates the workflow.

### 5.1 What Stage 2 might add

**Cases not in shipped saved_results:**
- Cox-Munk surface validation (shipped suite uses Lambertian)
- Plane-parallel runs of Tasks 3-6 (without sphericity) for clean FO comparison
- Glint-angle geometries
- Wind-speed sweep (1, 3, 5, 10 m/s)
- Multi-aerosol combinations
- Higher Nstreams sensitivity
- **Validation of `ExactSFIPhase` (architecture doc §6) against custom VLIDORT runs:** plane-parallel + DO_FOCORR=True + DO_FOCORR_OUTGOING=True; compare against `FullMOM(sfi_phase=ExactSFIPhase())`. The expected result: closer agreement than back-correction (FO-equivalent) achieves, because `ExactSFIPhase` corrects more paths.

**Wider parameter sweeps:**
- Albedo sensitivity
- Aerosol asymmetry parameter sweep
- τ_total sweep

### 5.2 PyVLIDORT setup (Stage 2 only)

By Stage 2, Stage 1 has established the comparison harness, parsing logic, and case-data conventions. Adding new cases is straightforward.

---

## 6. Test infrastructure

[Section 6 unchanged from v0.3 — fixture format, comparison harness, test runner. Per-case Julia files in `test/vlidort_baseline/stage{1,2}/case_*.jl`.]

---

## 7. Tolerances

[Section 7 unchanged from v0.3.]

- Scalar I: `rtol = 1e-5` (clean cases), `1e-4` (harder)
- Vector Q: `rtol = 1e-4`
- Vector U: `rtol = 1e-3` (sign conventions can differ)
- Absolute floor: `atol = 1e-9`

Per-case overrides for known divergences. Tighten as understanding improves.

---

## 8. CI integration

[Section 8 unchanged from v0.3.]

After Stage 1 demonstrates the workflow, suite runs on PRs touching `src/CoreRT/`. Manual trigger available.

---

## 9. Cross-references

**Delivered together with this document:**
- `vsmartmom_dispatch_design_v0_6.md` — architecture; especially §6 (`ExactSFIPhase`) which is the production correction validated separately from this Stage 1's FO-equivalent comparisons.
- `standalone_ss_solver_plan.md` — Piece A; uses the same Stage 1 fixtures for cross-validation. The back-correction adapter (Piece A §7) is validated against Stage 1 Case D.
- `dev_notes/exact_ss_reference/` — alternate validation source for SS-only cases.

**Pre-existing committed in PyVLIDORT-main distribution:**
- `vlidort_v_test/V2p8p3_Siewert2000_validation.f90` — Siewert 2000 setup
- `vlidort_v_test/saved_results/gfortran/results_Siewert2000_validation.all` — Siewert 2000 reference
- `vlidort_s_test/2p8p3_solar_tester.f90` — multi-task scalar setup
- `vlidort_s_test/saved_results/gfortran/results_solar_tester.all` — multi-task scalar reference
- `vlidort_v_test/V2p8p3_solar_tester.f90` — multi-task vector setup
- `vlidort_v_test/saved_results/gfortran/nstokes3/results_solar_tester_IQU0.all` — multi-task vector reference
- `vlidort_main/regular/vlidort_inputs.f90:4138-4141` — `NMOMENTS = min(2·NSTREAMS-1, NGREEK_MOMENTS_INPUT)`
- `vlidort_focode/FO_VectorSS_RTCalcs_I.f90` — VLIDORT FO source

---

## 10. Open questions

**BQ1**: PyVLIDORT version pinning. The shipped saved_results were generated with VLIDORT 2.8.3; document this and only update when intentionally regenerating.

**BQ2**: Number of Stage 2 cases — diminishing returns past ~80; bias toward fewer, cleaner cases.

**BQ3**: Suite as separate package vs `test/` in vSmartMOM. Current proposal: `test/`.

**BQ4**: CI policy — Stage 1 on every PR; full Stage 2 nightly or weekly.

**BQ5**: Polarization tolerance — Initial scalar; vector with per-case tolerance.

**BQ6**: Vector RT scope — `n_stokes=3` initial; `n_stokes=4` (with V) deferred indefinitely.

**BQ7**: Should we run vSmartMOM through the back-correction (Piece A §7) on Stage 1 cases, in addition to vanilla vSmartMOM, to validate the adapter? Yes — that's Case D in §4.1. Adds a third comparison column to results.

**BQ8 (new)**: When `ExactSFIPhase` is implemented, validate against custom plane-parallel-with-FO PyVLIDORT runs as Stage 2 expansion. The expected validation outcome: `ExactSFIPhase` matches VLIDORT-with-FO better than the back-correction does, because it corrects higher-order paths the back-correction can't reach. Quantifying this difference is the strongest validation of the `ExactSFIPhase` architecture choice.

**BQ9 (new)**: Sphericity in `solar_tester` Tasks 3-6. Stage 1 accepts the caveat; Stage 2 generates custom plane-parallel comparison runs for clean validation.
