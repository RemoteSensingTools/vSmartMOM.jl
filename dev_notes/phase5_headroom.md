# Phase 5 GPU headroom — baseline findings (2026-04-22)

Captured against `sanghavi-unified @ 16d074e` (post-Phase 4 InteractionWorkspace + post-greek_cabannes fix). 2× A100-PCIE-40GB, AMD EPYC 7H12, Julia 1.12.5.

## Baseline numbers (Phase1b_RRS_761-764nm, 103 spectral points)

| Scenario | Wall median | CPU alloc | GC time | GPU peak delta |
|---|---|---|---|---|
| RRS rt_run CPU  | 6.88 s | 9.02 GB | 1.76 s | – |
| RRS rt_run GPU  | 5.12 s | 9.12 GB | **2.45 s** | **2.98 GB** |
| noRS rt_run CPU | 0.30 s | 0.38 GB | 0.05 s | – |
| noRS rt_run GPU | 0.19 s | 0.12 GB | 0.00 s | 0 |
| noRS rt_run_ss CPU | 0.042 s | 0.07 GB | 0.00 s | – |
| noRS rt_run_ss GPU | 0.041 s | 0.09 GB | 0.00 s | 0 |

Stored metrics: [baseline_output/sanghavi-unified_16d074e_gpu/](../test/benchmarks/baseline_output/sanghavi-unified_16d074e_gpu/) — RRS, noRS, ss.

## Is there headroom?

**Yes, for RRS only.** noRS and single-scatter are already tight; the elastic path fully amortizes the adding-doubling CUBLAS calls.

RRS on GPU still allocates 9.12 GB of host memory per run and spends **2.45 s of 5.12 s total wall in GC** — more than the CPU path's GC (1.76 s). The 2.98 GB GPU peak delta shows significant CuArray turnover that isn't absorbed into the InteractionWorkspace. At 47% GC-fraction on the GPU run, the kernel is starved.

The 103-point RRS config is small enough that per-band fixed overheads (kernel launches, reductions) dominate. On larger RRS configs (OCO O₂-A ~6800 points, EMIT-Balsamic) the allocation cost scales as O(nλ²) for the coupling matrices, so the absolute gap gets worse. This baseline is the *lower bound* on Phase 5 benefit.

## Optimization targets (ordered by expected yield)

From the v2 plan amendments §3 — each gated on Christian sign-off.

1. **Absorb doubling-side intermediates into InteractionWorkspace.** `tmp3-tmp6`, `gp_refl`, `ieJ₁⁺/⁻` are currently per-call `similar(...)` allocations inside `doubling_inelastic!`. Candidate for the same workspace-hoisting that Phase 4 applied to `interaction_inelastic!`. Expected: drop a large chunk of the 2.98 GB GPU-peak delta and cut GC time.

2. **`batched_mul!` vs `⊠` on current shapes.** Sanghavi reported a 5.5× slowdown from a naive `@kernel` rewrite on 15×15 (memory confirms: "Kernel approach failed — custom @kernel 8x slower than CUBLAS for 15×15"). Do NOT replace `⊠`. The specific measurement target: whether `batched_mul!` with pre-reshaped contiguous 3-D buffers beats `⊠` on current shapes (15×15×nSpec). Benchmark-only phase first.

3. **Sync-point reduction on non-staged GPU path.** `staged=true` paged 4-D ie buffers through CPU between passes. The non-staged variant has more host-device syncs than necessary. Audit `@sync`/`CUDA.synchronize` calls.

4. **4D → 3D flattening + gather-based `batched_mul`.** The 4-D ieR/ieT buffers could be flattened to 3-D via a gather kernel, enabling CUBLAS `batched_mul!` end-to-end without per-slice dispatch. Larger refactor — only if (1) and (2) don't close the gap.

## Suggested measurement protocol (to propose to Christian)

For each candidate:

1. **Before-after wall + alloc + GC** on the three Phase 1b GPU scenarios (this baseline).
2. **Scale test** on a larger RRS config (creategrid_O2Aband_RamanSIF or similar) — the small config can hide gains that would be visible at O(nλ²) scale.
3. **Correctness gate**: Phase 1b RRS regression (already `rtol=0.02` post-fix, tight enough to catch regressions without being noise-limited).
4. **Merge if**: ≥20% wall reduction OR ≥30% GPU peak reduction OR ≥50% GC reduction on the scale test, with no Phase 1b regression. Otherwise: revert and document in this note.

## Status

**Paused pending Christian sign-off** per `plans/PLAN_AMENDMENTS_2026-04-19.md` §3. This file is the artifact for that decision.

## Files

- Metrics: `test/benchmarks/baseline_output/sanghavi-unified_16d074e_gpu/`
- CPU baselines: `test/benchmarks/baseline_output/sanghavi-unified_d84b789/`
- Scenarios: `test/benchmarks/harness/scenarios.toml` (CPU) and `scenarios_gpu.toml` (GPU)
