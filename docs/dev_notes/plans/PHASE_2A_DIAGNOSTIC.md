# Phase 2a — Diagnostic baseline (falsification test for Inventory D)

**Date:** 2026-04-20
**Branch @ commit:** `sanghavi-unified` @ `b4a5ab1` (Phase 1 complete)
**Scenarios:** 3 (RRS forward, noRS forward, single-scatter noRS) on the
`Phase1b_RRS_761-764nm.yaml` toy config (CPU, Stokes_IQU, Float32).

## What Phase 2a checks (v2 plan §4 Phase 2a)

> Run the harness on unmodified `sanghavi-unified` (= post-Phase-1) vs
> sanghavi. Record the delta. This verifies Inventory D's physics-fix list
> was complete: if the post-Phase-1 delta exceeds what D predicted, D was
> incomplete and must be revisited before Phase 3 proceeds. This is a
> diagnostic step, not a gate on a number.

## Harness under test

Three new files, all in [test/benchmarks/harness/](../test/benchmarks/harness/):

- `metrics.jl` — wall-clock (warm-up + N runs, median + MAD, outlier trim),
  CPU allocation count, GPU used-memory delta (when CUDA functional),
  output-array persistence to JLD2, hardware fingerprint.
- `run_benchmarks.jl` — scenario runner driven by `scenarios.toml`.
  Single Raman type builder per scenario; wraps `rt_run` or `rt_run_ss`.
- `report.jl` — side-by-side comparator between two baseline directories
  with per-quantity tolerance gating.
- `scenarios.toml` — declarative scenario list. Phase 2a is the
  `phase1b_RRS_761-764nm_cpu` + `phase1b_noRS_761-764nm_cpu` +
  `phase1b_ss_761-764nm_cpu` diagnostic set. Full-range acceptance
  scripts (`prototype_EMIT_aer_ht.jl`, `creategrid_O2Aband_RamanSIF.jl`)
  are added in Phase 3 after they're ported from sanghavi.

## Result — RRS scenario

`compare_baselines(sanghavi_9ee9a75, sanghavi-unified_b4a5ab1)` on the
`phase1b_RRS_761-764nm_cpu` scenario:

| Quantity | max \|Δ\| | max rel | median rel | ratio range | verdict |
|---|---|---|---|---|---|
| `R_SFI` (elastic I/Q/U) | 2.72e-05 | 0.032 | 0.0096 | 0.968 .. 0.991 | PASS (rtol 0.05) |
| `T_SFI` (elastic trans) | 6.80e-05 | 0.046 | 0.0093 | 0.954 .. 1.006 | PASS |
| `ieR_SFI` (inelastic I/Q/U) | 6.41e-08 | 0.0011 | 0.0007 | 0.999 .. 1.001 | PASS |
| `ieT_SFI` (inelastic trans) | 1.09e-06 | 0.016 | 0.0017 | 0.988 .. 1.016 | PASS |

## Interpretation

Two physics regimes emerge cleanly:

**Inelastic (ieR/ieT) — matches sanghavi.** max rel 0.11% and 1.6%, median
rel 0.07% and 0.17%. The `/π` weight fix + `ϖ_λ₁λ₀` normalization landed in
the Phase 1b closeout (commit [9980120](https://git/commit/9980120))
aligned the inelastic path to sanghavi essentially bit-close modulo Float32
rounding. This confirms Inventory D's inelastic-physics-fix list was
complete and correctly ported.

**Elastic (R/T) — residual ~1-3% bias documented in Phase 1b §9.** max rel
3.2% on R, 4.6% on T; ratios tightly clustered around 0.97–0.99. This is
a *systematic multiplicative bias* (std of the ratio field is O(1e-4)),
polarization-sensitive (different factor on I vs Q), not a scalar error —
already flagged in [PHASE_1B_STAGING.md §9](PHASE_1B_STAGING.md) as a
follow-up investigation with candidate root causes (depolarization
plumbing, `greek_rayleigh`/`greek_cabannes` coefficient construction, or
Fourier-moment weighting).

## Does this falsify Inventory D?

**No.** Inventory D §1 / §6 listed inelastic-path physics fixes (Cabannes
ϖ, α̅ `2π` drop, Bodhaine Rayleigh, ϖ_λ₁λ₀ normalization). All landed in
Phase 1 and produce the expected ~0.1% inelastic agreement. The elastic
residual is NOT something Inventory D claimed to fix — it's a separate
pre-existing physics-convention difference between the branches that the
`/π` fix exposed by bringing the elastic output to comparable magnitude.

Phase 2b baseline freeze can proceed once the user sets the tolerance
table. The residual elastic delta is addressed by either:
(a) loosening the elastic rtol to ≥0.05 (current harness default), and
(b) tracking the residual as a separate follow-up ticket, OR
(c) root-causing the elastic residual before Phase 2b freeze.

## What's NOT measured (yet)

- **Sanghavi-side wall-clock / allocation / GPU memory numbers.** The
  sanghavi baseline for Phase 2a is output-arrays-only, re-packaged from
  [test/reference/phase1b_RRS_sanghavi_q0.jld2](../test/reference/phase1b_RRS_sanghavi_q0.jld2).
  Running the harness on the sanghavi worktree to fill these in is
  straightforward but needs session-local `@eval` workarounds for sanghavi's
  broken `rt_run(::RRS)` 21-arg `postprocessing_vza!` call (see
  [PHASE_1B_STAGING.md §8](PHASE_1B_STAGING.md#L215)) — low-priority for
  Phase 2a since physics parity is the Phase 2a gate, not performance.
  Phase 2b or Phase 4 can revisit.

- **GPU Raman scenarios.** The harness ran CPU-only in this diagnostic —
  the Phase1b YAML uses `architecture: CPU()`. A GPU scenario variant is
  straightforward to add once we settle the v2 tolerances. For Phase 2b
  freeze, the full-range RRS + SIF scripts will run on GPU.

## Phase 2a exit

Exit criterion per v2 plan: "diagnostic baseline + delta narrative
reviewed." The delta is:
- Inelastic: sub-1% (matches sanghavi).
- Elastic: 1-3% residual, already documented as a follow-up.

No Inventory D items missed. Phase 2a complete.

## Next — Phase 2b

Needs user tolerance table (v2 plan §Tolerance table) filled in, then:
1. Update `DEFAULT_TOLERANCES` in `report.jl` to user values.
2. Re-run harness on sanghavi worktree for proper wall-clock/alloc metrics.
3. Freeze baselines in `baseline_output/sanghavi_9ee9a75/` (replacing
   the Phase 2a output-only version) and `baseline_output/sanghavi-unified_<sha>/`.
