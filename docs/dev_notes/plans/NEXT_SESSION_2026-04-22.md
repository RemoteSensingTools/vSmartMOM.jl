# Pickup plan for 2026-04-22

Written 2026-04-21 end of session. Branch `sanghavi-unified` @ `1e08840`. All commits pushed to the local branch; `origin/sanghavi-unified` is behind by 13 (unpushed).

## State of the tree

Phases 0 → 4 complete. Phase 5 conditional on GPU measurement + Christian sign-off (per amendments §3).

| Area | Status |
|---|---|
| SS fix (91% error) | Closed in `776ff00`. R_SFI matches sanghavi to FP32 noise. |
| SIF injection + acceptance scripts (Phase 3) | `0178a7d..ef1db27`. test_sif.jl 21/21 pass; both EMIT bands + RRS+SIF grid run to completion. |
| Bodhaine + profile alignment | `55b3001`. g₀=9.8032465, per-band greek_rayleigh from molecular constants, dropped q/1000 divisor. |
| linRT bounds fix | `89a8fb4`. `surface_index(layout, 1)` instead of `surface_index(layout, iBand)`. |
| Phase 4 workspace port | `d84b789..1e08840`. `InteractionWorkspace` wired; 12/12 phase1b pass; CPU allocs 26.5GB → 9.78GB (-63%); fresh harness baselines @ `baseline_output/sanghavi-unified_d84b789/`. |

## Known open items

### (a) Elastic R/T residual — CLOSED 2026-04-22

Root cause: unified's `compEffectiveLayerProperties.jl:22` unconditionally used `greek_rayleigh` (depol ≈ 0.028) for the Rayleigh Z moments, while sanghavi branches on `RS_type <: noRS` and uses `greek_cabannes` (depol ≈ 0.007) for RRS/VS. The Cabannes phase matrix has ~3% larger polarization-sensitive greek coefs (β, δ) — exactly the observed residual.

Fix: [compEffectiveLayerProperties.jl](../src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl) now picks `greek_cabannes` when `RS_type` is Raman-active. Phase 1b RRS ratios all ≈ 1.000 within FP32 precision (mean <1e-4 deviation). Test tolerance tightened from `rtol=0.05` → `rtol=0.02`. See [PHASE_1B_STAGING.md §10](PHASE_1B_STAGING.md).

### (b) Phase 5 — GPU headroom (conditional)

Christian sign-off blocking per-optimization per amendments §3. Candidates from the plan:
- `doubling-side` allocations (`tmp3-tmp6`, `gp_refl`, `ieJ₁⁺/⁻`) not yet absorbed into workspace
- `batched_mul!` vs `⊠` on actual problem shapes (sanghavi's 5.5× slowdown was on 15×15)
- sync-point reduction on non-staged GPU path
- 4D→3D flattening + gather-based batched_mul

**Pre-work**: run the harness on GPU (instead of CPU) to establish the GPU baseline post-Phase 4. If GPU shows significant headroom remaining beyond sanghavi, Phase 5 is real; otherwise, Phase 5 exits as a `dev_notes/phase5_headroom.md` no-op.

### (c) Phase 6 — Remaining script port

Pre-flight: produce `plans/phase6_script_port_list.md` — a triage of sanghavi's `test/benchmarks/*.jl` scripts (OCO, CarbonI, O3 Huggins, RamanSIF variants, EMIT/Balsamic applications, retrieval comparisons). For each: runtime, data-file dependencies, keep/port/drop decision. User signs off before porting starts ("we thought 8, it's actually 23" is a real failure mode for research-code ports).

### (d) Phase 7 — Docs overhaul

Pre-flight tasks:
1. Verify each `docs/src/pages/tutorials/Tutorial_*.jl` runs against current `sanghavi-unified` tip. Literate.jl tutorials that error through `Base.getproperty` shim need fixing before `warnonly=false` flips.
2. Run `julia --project=docs docs/make.jl` WITHOUT `warnonly=true` once and count actual errors. Scope the backfill work from the count — don't use v1's "30-50" estimate.

Then: README rewrite, CHANGELOG deepening, 5 examples + reference outputs, migration guide, narrative tutorial. Test-suite orphan wire-in (`test_forward_lin.jl`, `test_hybrid_ad.jl`, `test_jacobians_GPU.jl`, `test_mie_gpu.jl`, `test_performance.jl`).

### (e) Phase 8 — Merge PR

Bump version in `Project.toml` (2.0.0 → 2.1.0), update CHANGELOG, open PR `sanghavi-unified → main`. No Christian sign-off required.

## Recommended ordering for 2026-04-22

1. ~~Push today's commits~~ — done (30 commits pushed to origin/sanghavi-unified).
2. ~~Bisect the elastic R/T residual (item a)~~ — **done**, root-caused and fixed (greek_cabannes dispatch). Phase 1b RRS 12/12 pass at `rtol=0.02`.
3. **Run harness on GPU** for Phase 5 measurement. If headroom exists, brief Christian on specific optimization targets and measurement protocol. If not, write `dev_notes/phase5_headroom.md` no-op and proceed to Phase 6.
4. **Phase 6 pre-flight**: `plans/phase6_script_port_list.md` triage.

## Files / commits to reference

- `plans/IMPLEMENTATION_PLAN_v2.md` — master phase scaffold.
- `plans/PHASE_1B_STAGING.md §9` — elastic residual hypotheses.
- `plans/SS_DIAGNOSIS_WIP.md` — SS fix resolution (completed 2026-04-21).
- `plans/PLAN_AMENDMENTS_2026-04-19.md` — authority rule.
- `test/benchmarks/baseline_output/sanghavi-unified_d84b789/` — Phase 4 baselines (supersedes `_3f876c9` and `_b4a5ab1`).
- `test/reference/phase1b_RRS_sanghavi_q0.jld2` — sanghavi-authoritative physics gate.

## Environment notes

- `test/Project.toml` has JSON, TOML, StatsBase. vSmartMOM is dev'd to `..`. If resolver reverts to registered v1.1.0 (lacks altitude-form z₀/σ₀), run `julia --project=test -e 'using Pkg; Pkg.develop(path=".")'`.
- MCP julia session auto-dev's vSmartMOM via project path. Standalone `julia --project=test ...` needs the explicit dev.
- 2× A100 GPUs available; CUDA functional in sessions with `using CUDA` first.
- Sanghavi worktree at `/home/sanghavi/code/github/vSmartMOM.jl/` is READ-ONLY. Any re-runs require `@eval` monkey-patching on the session side.
