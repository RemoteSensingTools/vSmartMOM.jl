# Phase 2b — Baseline freeze (sanghavi + sanghavi-unified)

**Date:** 2026-04-21
**Sanghavi branch @:** `9ee9a75` (read-only worktree)
**sanghavi-unified @:** `3f876c9` (Phase 1 complete + Phase 2a diagnostic)
**Scenarios:** 3 — `phase1b_RRS_761-764nm_cpu`, `phase1b_noRS_761-764nm_cpu`,
`phase1b_ss_761-764nm_cpu`. All on the `Phase1b_RRS_761-764nm.yaml` toy config
(CPU, Stokes_IQU, Float32, 103 spec pts, 15 layers).

## User-specified tolerances (2026-04-21)

- **I (intensity):** `rtol = 1e-4` (better than 0.01% relative).
- **Q, U, V:** `atol = 1e-8` (absolute).
- **Performance:** sanghavi is the ceiling; any wall / alloc / GPU metric
  above sanghavi is **flagged**, not auto-failed. Iterate.

## Sanghavi-side capture

[test/benchmarks/harness/capture_sanghavi.jl](../test/benchmarks/harness/capture_sanghavi.jl)
is a one-off runner launched from the sanghavi worktree with `--project=.`.
It loads vSmartMOM from sanghavi's tree, applies a session-local monkey-patch
that forwards the 21-arg `postprocessing_vza!` call (which sanghavi's rt_run
does for RRS) to sanghavi's 18-arg RRS implementation, and writes metrics +
output JLD2 into this repo's [test/benchmarks/baseline_output/sanghavi_9ee9a75/](../test/benchmarks/baseline_output/sanghavi_9ee9a75/).

## Results — per-Stokes tolerance gate

Full report: `include("test/benchmarks/harness/report.jl"); compare_baselines(...)`.
Condensed:

### `phase1b_RRS_761-764nm_cpu`

| Quantity | max\|Δ\| | max rel | Gate | Verdict |
|---|---|---|---|---|
| R_SFI[I] | 2.72e-5 | 0.97% | rtol≤1e-4 | **FAIL** |
| R_SFI[Q] | 4.89e-6 | 3.19% | atol≤1e-8 | **FAIL** |
| R_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |
| T_SFI[I] | 6.80e-5 | 2.42% | rtol≤1e-4 | **FAIL** |
| T_SFI[Q] | 7.15e-6 | 4.63% | atol≤1e-8 | **FAIL** |
| T_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |
| ieR_SFI[I] | 6.41e-8 | 0.108% | rtol≤1e-4 | **FAIL** |
| ieR_SFI[Q] | 5.35e-10 | 0.108% | atol≤1e-8 | PASS |
| ieR_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |
| ieT_SFI[I] | 1.09e-6 | 1.56% | rtol≤1e-4 | **FAIL** |
| ieT_SFI[Q] | 7.44e-9 | 1.45% | atol≤1e-8 | PASS |
| ieT_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |

Perf:

| Metric | sanghavi | sanghavi-unified | Ratio | Flag |
|---|---|---|---|---|
| wall_median_s | 5.858 | 6.364 | 1.09× | **FLAG** (>sanghavi) |
| cpu_alloc_bytes | 9.38 GB | 8.98 GB | 0.96× | ok |

### `phase1b_noRS_761-764nm_cpu`

| Quantity | max\|Δ\| | max rel | Gate | Verdict |
|---|---|---|---|---|
| R_SFI[I] | 1.24e-6 | 0.043% | rtol≤1e-4 | **FAIL** |
| R_SFI[Q] | 4.22e-8 | 0.028% | atol≤1e-8 | **FAIL** |
| R_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |
| T_SFI[I] | 4.35e-5 | 1.54% | rtol≤1e-4 | **FAIL** |
| T_SFI[Q] | 2.38e-6 | 1.58% | atol≤1e-8 | **FAIL** |
| T_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |
| ieR/ieT · all | 0 | 0 | — | PASS (noRS) |

Perf:

| Metric | sanghavi | sanghavi-unified | Ratio | Flag |
|---|---|---|---|---|
| wall_median_s | 0.3137 | 0.2987 | 0.95× | ok |
| cpu_alloc_bytes | 461 MB | 400 MB | 0.87× | ok |

### `phase1b_ss_761-764nm_cpu` (single-scatter)

| Quantity | max\|Δ\| | max rel | Gate | Verdict |
|---|---|---|---|---|
| R_SFI[I] | 2.54e-3 | **91.5%** | rtol≤1e-4 | **FAIL** |
| R_SFI[Q] | 1.38e-4 | **91.5%** | atol≤1e-8 | **FAIL** |
| R_SFI[U] | 0 | 0 | atol≤1e-8 | PASS |
| T_SFI[I] | 2.54e-3 | **91.5%** | rtol≤1e-4 | **FAIL** |
| T_SFI[Q] | 1.38e-4 | **91.5%** | atol≤1e-8 | **FAIL** |
| hem_R[I] | 1.10e-2 | 91.3% | rtol≤1e-4 | **FAIL** |
| hem_T[I] | 1.10e-2 | 91.3% | rtol≤1e-4 | **FAIL** |
| ieR/ieT · all | 0 | 0 | — | PASS (noRS) |

Perf:

| Metric | sanghavi | sanghavi-unified | Ratio | Flag |
|---|---|---|---|---|
| wall_median_s | 0.0430 | 0.0466 | 1.08× | **FLAG** |
| cpu_alloc_bytes | 123 MB | 57 MB | 0.47× | ok |

## Flagged exceedences — iteration backlog

Physics (all **blocking** at user tolerance; iterate):

1. **Elastic R[I]/R[Q]/T[I]/T[Q] residual on RRS + noRS** — already flagged in
   [PHASE_1B_STAGING.md §9](PHASE_1B_STAGING.md#L229). Systematic
   polarization-sensitive multiplicative bias, ratios 0.968-0.991. The
   Phase 2b data shows it is WORSE on RRS (1-3% rel) than on noRS (<1.6% rel),
   suggesting a Raman-Cabannes coupling term in the elastic path. Candidates
   remain: depol plumbing, greek_rayleigh / greek_cabannes construction,
   Fourier moment weighting.
2. **Inelastic ieR/ieT[I] at 0.1-1.6% rel** — new information. The Phase 1b
   ~0.1% agreement was under a looser rtol=0.05 gate. At the user's 1e-4
   target, inelastic I also exceeds — by 10× on ieR and 100× on ieT. Likely
   same root cause as (1) since ieR/ieQ agree to the same 0.1% factor.
3. **Single-scatter 91% error** — **this is new and large.** The ported
   `rt_run_ss` produces numbers that differ from sanghavi by a factor
   close to 1 - 0.085 ≈ 0.915 on both R and T and hem. Hypothesis: the
   Phase 1c port's 6-tuple `hem_R`/`hem_T` return picked up a factor of
   π or similar that sanghavi doesn't have in the SS path; or an
   inconsistency in the weight = `0.5/π` vs `0.5` convention survived the
   port in one of the two SS driver bodies. Needs a dedicated root-cause
   before Phase 3.

Performance (flags only — not blocking):

4. **RRS wall-clock +9%** — sanghavi-unified is slower on RRS at CPU/Float32.
   Sanghavi's `InteractionWorkspace` landing is Phase 4 and will likely
   close this. No action required at Phase 2b.
5. **SS wall-clock +8%** — similar. The SS driver path is narrow enough
   that this can be revisited if (3) root-causing touches the driver body.
6. **Alloc counts all ≤ sanghavi** — good. RRS 0.96×, noRS 0.87×, SS 0.47×.
   The Phase 4 workspace port should preserve or improve these.

## What Phase 2b freezes

- [test/benchmarks/baseline_output/sanghavi_9ee9a75/](../test/benchmarks/baseline_output/sanghavi_9ee9a75/) —
  3 scenarios with wall/alloc/output. Supersedes the Phase 2a output-only
  version.
- [test/benchmarks/baseline_output/sanghavi-unified_3f876c9/](../test/benchmarks/baseline_output/sanghavi-unified_3f876c9/) —
  3 scenarios from post-Phase-1 sanghavi-unified tip.
- [test/benchmarks/harness/report.jl](../test/benchmarks/harness/report.jl) —
  DEFAULT_TOLERANCES frozen with per-Stokes-component gates (I rtol 1e-4;
  Q/U/V atol 1e-8).

## Phase 2b exit criterion

Per v2 plan §Phase 2 exit: "Two baseline directories committed;
`report.jl` green on physics; user reviews the wall-clock / alloc gap
numbers and confirms the Phase 4 target."

**Not exited yet.** The user said "work iteratively" — Phase 2b baselines
are frozen and committed, but the physics gate is RED. Next-session
actions:

1. User picks which FAIL to tackle first. Recommended order:
   - SS 91% error (biggest, newest, likely cleanest root cause).
   - Elastic residual on RRS (ripples into ieR/ieT at the 1e-4 gate).
   - Inelastic I residual at 1e-4 (depends on whether (2) closes (1)).
2. Root-cause + fix each item; re-capture baselines; report delta.
