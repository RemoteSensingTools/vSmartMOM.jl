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

### (a) Elastic R/T residual — not closed by Bodhaine alignment

Phase 1b RRS regression vs sanghavi reference:
- Stokes I ratio: **0.9907** (≈ 1% delta)
- Stokes Q ratio: **0.9683** (≈ 3% delta)
- ieR/ieT ratios: within 0.1% (effectively matched)
- Both ratios are constant across VZA / spectrum (std < 5e-5)

The **depol delta** between user-supplied 0.028 and computed
`depol_air_Rayleigh = 0.0281` is 0.5% → propagates to <0.01% on phase
coefs. Not the culprit.

**Next hypotheses to bisect** (from [PHASE_1B_STAGING.md §9](PHASE_1B_STAGING.md)):
1. **Elastic Cabannes phase routing**: unified `compEffectiveLayerProperties.jl` uses `greek_rayleigh` (now per-band, depol≈0.0281) for the Rayleigh Z-moments. Sanghavi does the same at `rt_run.jl:84`. Verify bit-exact agreement on `Rayl𝐙⁺⁺` / `Rayl𝐙⁻⁺` between branches for the same layer.
2. **`construct_atm_layer` elastic weighting**: sanghavi `atmo_prof.jl:410-414` weights τ_rayl × ϖ_Cabannes × Rayl𝐙⁺⁺. Unified `compEffectiveLayerProperties.jl:24` does `CoreScatteringOpticalProperties(τ_rayl, ϖ_Cabannes, Rayl𝐙⁺⁺, Rayl𝐙⁻⁺)`. Verify the downstream consumer applies the same weighting.
3. **Q-specific terms (δ, γ)**: the ratio-of-ratios Q/I = 0.977 is suggestive of a fixed 2.3% multiplicative bias on Q. Check `elemental!`'s `δ * dpl_r * 1.5` term and `apply_D_matrix_elemental!` Stokes sign convention.
4. **Fourier-moment weight `wct02`**: m=0 → 0.5, m>0 → 0.25 in both branches; verify identical.

**Suggested attack**: write a minimal diagnostic that dumps `layer_opt[iz].Z⁺⁺` and `computed_layer_properties.Z⁺⁺` (first layer, m=0) from both branches on identical YAML, diff elementwise. If they're bit-identical, the residual is in elemental / doubling / interaction. If they differ, localize to `compEffectiveLayerProperties` or `compute_Z_moments`.

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

1. **Push today's commits** (`git push origin sanghavi-unified`) — 13 commits currently unpushed.
2. **Bisect the elastic R/T residual** (item a). Two hours tops before switching. If not root-caused by then, write up findings in `PHASE_1B_STAGING.md §10` and move on.
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
