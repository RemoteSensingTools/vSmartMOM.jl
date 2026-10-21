# Fused adding-doubling kernels — plan (2026-08-20)

## Problem, measured

The CoreRT adding-doubling solver is **kernel-launch-bound on GPU**. Ungated
Nsight trace of two warm 801-pt nstr=8 IQU F32 solves (A100, driven from the
Google-RT satellite campaign):

- ~69 000 kernel launches, total GPU busy **≈ 0.45 s inside ~50 s wall** —
  the GPU is idle > 99 % of the solve.
- Kernel inventory: ~53 k `gpu_broadcast_kernel_linear` (4–6 µs each),
  ~15 k Boolean `partial_mapreduce_grid` (sync-forcing checks), a few
  hundred real physics kernels (cross-section, Mie, Greek, fused GP).
- `cuLaunchKernel` API is only 5.7 µs/launch; the cost is the ~0.5 ms of
  host-side Julia work around each launch (dispatch, `⊠` output
  allocation, GC). This is the measured **npts-independent ≈ 18–20 s per
  solve** — spectral points are nearly free (25× points → 1.38× time),
  launches are not.

Line-level attribution (isolated, production shapes 24×24×801 F32,
`test/benchmarks/launch_attribution_benchmark.jl`):

| expression | kernels/call |
|---|---|
| pure dotted chain (broadcast fusion works) | 1 |
| single `⊠` (NNlib batched_mul, allocates output) | 4 |
| `r⁻⁺ .= r⁻⁺ .+ tt_gp ⊠ r⁻⁺ ⊠ t⁺⁺` | 7 |
| `j₀⁻ .= j₀⁻ .+ tt_gp ⊠ (j₁⁻ .+ r⁻⁺ ⊠ j₀⁺)` | 8 |
| `compute_geometric_progression!` — **already fused** (`ka_fused_gp_solve!`) | **2** |
| `any(isnan, ·)` Bool reduce | 3 + sync |
| **`doubling_source_update!` total** | **18** |
| **`doubling_rt_update!` total** | **12** |

Per doubling step ≈ 30–50 launches (more with per-source `j₀_by_src`
slots), × ndoubl × scattering layers × (m_max+1) ⇒ the observed totals.

`compute_geometric_progression!` is the existence proof: its hand-written
fused KA kernel does the mathematically hardest step (build E−R·R, LU
solve, right-multiply) in 2 launches where the fallback needs ~10 — and it
is **backend-agnostic** (KernelAbstractions), which is a requirement:
NVIDIA-only solutions (CUDA graphs) are a side experiment, not the plan.

## Plan

Phase 0 — **preallocation pass** (prerequisite for everything, incl. the
CUDA-graph experiment on the sibling branch):
- in-place `batched_mul!` with preallocated outputs everywhere `⊠` appears
  in the hot path (`rt_helpers.jl`, `interaction.jl` — 32 `⊠` sites);
- hoist/batch the Boolean reductions (`any`/`isnan`-style checks) out of
  the per-step path (accumulate a device-side flag, check once per layer
  or per m).

Phase 1 — **fused doubling-step kernels** (the payoff):
- `ka_fused_doubling_source!`: one KA kernel per spectral point doing
  j₁± scaling, `r⁻⁺·j₀⁺`, `j₁⁻ + …`, `tt_gp·(…)`, and the accumulate —
  the 24×24 tiles fit in local memory exactly like `ka_fused_gp_solve!`.
  Replaces 18 launches with 1–2.
- `ka_fused_doubling_rt!`: same treatment for the R/T update (12 → 1–2).
- Optional: fold both + the fused GP solve into ONE per-doubling-step
  kernel (the whole step is ~10 small matmuls on the same tiles).
- Same fits-gating + fallback pattern as `_use_fused_gp`.

Phase 2 — **interaction fusion**: `interaction.jl`'s ⊠ chains (32 sites),
same recipe, after Phase 1 proves the approach.

Validation contract (repo convention): bitwise or tight-tolerance
equivalence vs the unfused path per backend, KA-CPU + CUDA, all precision
policies — mirroring the batched-Mie-seam validation discipline.

Isolated speed tests live in `test/benchmarks/launch_attribution_benchmark.jl`
(this branch); extend it with fused-vs-unfused A/B rows as kernels land.

## Status (2026-08-20)

Phase 0 ✔ and Phase 1 (doubling) ✔ landed on this branch:

| per call, 24×24×801 F32 | kernels | GPU time | host API |
|---|---|---|---|
| source_update, Phase 0 (_bmm!) | 18 | 0.090 ms | 0.220 ms |
| **source_update, Phase 1 fused** | **1** | **0.009 ms** | 0.066 ms |
| rt_update, Phase 0 (_bmm!) | 12 | 0.115 ms | 0.151 ms |
| **rt_update, Phase 1 fused** | **1** | **0.058 ms** | 0.076 ms |

A doubling step is now 3 kernels (fused GP + fused source + fused rt) vs
~32, at ~0.19 ms/step eager — and with cuBLAS out of the captured region
the full 12-step ladder capture/instantiate/RELAUNCHES as a CUDA graph
(all five bisect stages OK; see perf/cuda-graph-capture).

Validation: tolerance-equivalence tests on KA-CPU + real CUDA, F64@1e-13 /
F32@2e-5 (test_fused_doubling.jl); END-TO-END GPU rt_run A/B fused vs
kill-switch on DefaultParameters: max rel diff 1.1e-15. Kernels written
math-visible (adding equations inline at each compute block, Sanghavi
et al. 2014 Eq. references) per collaborator requirement.

Remaining: Phase 2 (interaction ⊠ chains, 32 sites — same recipe).

## Expected effect

Launches per doubling step: ~30–50 → ~3–5. End-to-end solve ceiling
(GPU-busy bound): ~30–100× at 801–4001 pts. For the satellite campaign:
nstr=8 full-band per-cell RT ~2.5 min → ~5–15 s; a full 3-band C90 day
~90 → ~5–10 GPU-days.

## Context: how JAX avoids this class of problem

XLA traces the whole jitted program and fuses elementwise work across
statement/function boundaries (no `⊠` fusion barrier), and `lax.scan`
keeps loops on-device inside one compiled executable — there is no per-op
host dispatch to pay. The Julia-native analogue is Reactant.jl (Julia →
XLA/MLIR); a Reactant port of the adding core is a heavier, longer-term
alternative to hand-fused KA kernels, noted here for the record.
