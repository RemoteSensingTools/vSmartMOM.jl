# Inventory A — Sanghavi Branch Optimization & Physics Delta

**As of:** 2026-04-19
**Purpose:** Reconcile what is on the `sanghavi` branch today vs what Christian Frankenberg's merge plan
(`docs/dev_notes/sanghavi_unified_merge_plan.md`, 2026-03-21) assumed. The merge plan predates
the `InteractionWorkspace` work (Apr 2026) and several physics updates; this inventory documents
the actual state so the `sanghavi-unified` merge does not regress sanghavi's optimizations or
physics.

**Editorial-pass note (2026-04-19):** The 2026-04-19 user review session established the **authority rule** — *sanghavi is the authority for the inelastic path*. The per-item "adopt wholesale / rewrite / discard" framing in §5 below and several of the open questions in §8 are **superseded** by that rule and by `plans/PLAN_AMENDMENTS_2026-04-19.md`. Specific redirections are marked inline. The evidence and factual observations stand unchanged.

## Branch tips

| Branch | Tip commit | Tip message | Path |
|--------|-----------|-------------|------|
| `sanghavi` | `9ee9a75` | "add plan" | `/home/sanghavi/code/github/vSmartMOM.jl/` |
| `unified-vsmartmom` | `a4e4187` | "Add batched-kernel and Raman scaling benchmarks with writeup" | `/home/sanghavi/code/github/uni_vSmartMOM/` |
| Divergence base | `702cbc3` (Apr 2022) | "Added multisensor simulations with Raman" | — |
| Unique sanghavi commits (`unified-vsmartmom..sanghavi`) | **51** | — | — |

---

## 1. Commit `854b44c` — Float32 support, InteractionWorkspace, optimize `get_n₀_n₁`

### Summary

Landed three coupled changes: (a) relaxed several `where FT` constraints so Float32 and mixed
precision work through the inelastic pipeline, (b) introduced an `InteractionWorkspace` struct to
pre-allocate GPU buffers for `ScatteringInterface_11` (the Raman interaction path), and (c)
replaced the `findall`-based `get_n₀_n₁` with pure index arithmetic. Verified Float64 bit-exact
against its own baseline; Float32 achieves ~50% memory reduction. Threads workspace through
`rt_kernel! → interaction! → interaction_helper!`.

### Abstractions introduced

- **`mutable struct InteractionWorkspace{A3, A4_mat, A4_src, A3_src}`** (first form, before the
  staging split): 3D buffers `tmp_inv, tmpR⁻⁺, tmpR⁺⁻, tmpT⁻⁻, tmpT⁺⁺, tmpJ₀⁻, tmpJ₀⁺` and 4D
  buffers `tmpieR⁻⁺, tmpieR⁺⁻, tmpieT⁻⁻, tmpieT⁺⁺, tmpieJ₀⁻, tmpieJ₀⁺` (shapes
  `nQuad × nQuad × nSpec [× nRaman]` and `nQuad × 1 × nSpec [× nRaman]`).
- **`InteractionWorkspace(composite_layer, added_layer)` constructor** — allocates via `similar()`
  against a CompositeLayerRS.
- **`reset!(ws::InteractionWorkspace)`** — zeroes all buffers.
- **Optional `workspace` keyword** added to `interaction!`, `interaction_helper!`, and every
  `rt_kernel!` method. Falls back to old `similar()` path when `workspace === nothing`.
- **Converting constructor** `CoreScatteringOpticalProperties(τ, ϖ, Z⁺⁺, Z⁻⁺)` that
  `promote_type`s mixed float inputs.
- **Relaxed operators** `Base.:+(x::Core…{FT1}, y::Core…{FT2})` and `Base.:*(x, y::Core…)` to
  permit mixed-precision scalar multiplication.
- **`get_n₀_n₁(ieJ₁⁺, Δ)`** rewritten as `n₁_start = max(1, 1 - Δ); n₁_end = min(nSpec, nSpec - Δ)`
  (no `findall`, no heap allocation — GPU-safe arithmetic).
- **Pol_type / Z-matrix FT coercion** inside `rt_run` and
  `constructCoreOpticalProperties`: if `FT ≠ eltype(pol_type.D)` the pol_type is re-wrapped, and
  `Rayl𝐙⁺⁺ / Rayl𝐙⁻⁺ / AerZ⁺⁺ / AerZ⁻⁺ / τ_abs` are `FT.(…)` on the fly.
- **`q::Union{Array{FT,1}, Nothing}`** on `AtmosphericProfile` and `vSmartMOM_Parameters` — makes
  specific humidity optional.

### Files & line ranges (at `854b44c`)

- `src/CoreRT/CoreKernel/interaction_inelastic.jl` — +55 lines at head (workspace struct + ctor + reset!),
  and workspace branches at lines ≈296–320, 447–475 (two `ScatteringInterface_11` methods: `RRS`,
  `VS_0to1_plus`/`VS_1to0_plus`).
- `src/CoreRT/CoreKernel/rt_kernel.jl` — added `workspace=nothing` kwarg to all four `rt_kernel!`
  methods (noRS legacy, noRS `CoreScattering…`, inelastic legacy, inelastic `CoreScattering…`)
  and plumbed it to `interaction!`.
- `src/CoreRT/rt_run.jl` — allocates `_interaction_ws2 = (typeof(RS_type) <: noRS) ? nothing :
  InteractionWorkspace(composite_layer, added_layer)` once before the `for iz = 1:Nz` layer loop
  and passes to both `rt_kernel!` (per layer) and the final `interaction!` with the surface
  layer. Same pattern added to `rt_run_bck.jl` and `rt_run_lin.jl`.
- `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl` — FT coercion for pol_type,
  Z matrices, τ_abs, ϖ_Cabannes scalar (`FT(ϖ_Cabannes[iB])`).
- `src/CoreRT/types.jl` — `q::Union{Array{FT,1}, Nothing}`, converting ctor, relaxed `+`/`*`
  operators (around lines 22, 416, 623, 716).
- `src/Inelastic/inelastic_helper.jl` — new `get_n₀_n₁` (lines 13–19); relaxed
  `getRamanAtmoConstants(ν̃::AbstractFloat, T::AbstractFloat)` with `promote`.
- `src/CoreRT/atmo_prof.jl`, `src/CoreRT/model_from_parameters.jl`,
  `src/CoreRT/lin_model_from_parameters.jl`, `src/CoreRT/parameters_from_yaml.jl`,
  `src/CoreRT/rt_run_bck.jl`, `src/CoreRT/rt_run_lin.jl`,
  `src/CoreRT/show_utils.jl`, `src/CoreRT/types_inelastic.jl`,
  `src/Inelastic/stellar_inelastic_helper.jl`, `src/Testing/perturb_parameters.jl` —
  FT-relaxation edits.
- `test/test_parameters/O2_parameters2_1band_opttest_f32.yaml` — new Float32 YAML (71 lines).
- `test/benchmarks/benchmark_natraj.jl`, `test/benchmarks/natraj.yaml` — minor test updates.

### Unified counterpart

- Unified already has Float32 support at its own API level (commits `c581224` "Add Float32/type
  stability tests" and `37bce51` "Ensure Float32 type stability") and a parametric `FT` stored
  on the `RTModel`. But unified still carries the OLD `get_n₀_n₁` with `findall` at
  `src/Inelastic/inelastic_helper.jl:3-11` and still has 12× `similar() + .=0` at the entry of
  both `ScatteringInterface_11` methods in `src/CoreRT/CoreKernel/interaction_inelastic.jl:293-305`
  and `:440-452`.
- Unified already has a pre-allocated workspace — but only for the **elastic** path:
  `RTWorkspace{FT, AT3}` in `src/CoreRT/types.jl:1147-1170` with fields `gp_refl, tt_gp, J₁⁺, J₁⁻,
  pivot, info, tmp_inv, T_inv, tmp3d_a, tmp3d_b`. The constructor is `make_rt_workspace(FT,
  arr_type, NquadN, nSpec)` at line 1177. There is NO analog for the inelastic side; Christian's
  merge plan flagged this gap and proposed `InelasticWorkspace` (see §Cross-comparison).

### Recommendation

**Adopt the concept wholesale; rewrite against unified's API.** Specifically:

1. **`get_n₀_n₁` arithmetic rewrite** — adopt verbatim. It is a pure perf/correctness improvement
   with no API surface: replace lines 3–11 of `src/Inelastic/inelastic_helper.jl` on unified with
   sanghavi's arithmetic version. This also removes a hidden CPU allocation from every
   `ScatteringInterface_11` call on GPU runs.
2. **`InteractionWorkspace` struct** — do **not** port verbatim. Unified already has `RTWorkspace`
   covering the 3D (`tmp_inv`, `tmpR*`, `tmpT*`, `tmpJ*`) portion that is shared between elastic
   and inelastic. Extend the unified-side design with a sibling `InelasticWorkspace` (or fields
   added to `RTWorkspace`) covering only the 4D `ie*` buffers plus the CPU staging arrays (see
   commit `d75dacb` below). Use sanghavi's struct layout as **design input**, not verbatim code,
   because the unified conventions are lowercase added-layer / uppercase composite-layer and use
   `(; field) = struct` unpacking rather than `@unpack`.
3. **Float32 relaxations to `CoreScatteringOpticalProperties` + operators + converting ctor** —
   re-implement on unified (its `CoreScatteringOpticalProperties{FT,FT2,FT3}` at
   `src/CoreRT/types.jl:1017-1118` is 3-parametric, so the signature change differs). The FT
   coercion in `constructCoreOpticalProperties` should be trivial to port because it's localized
   to Z-matrix creation.
4. **`q::Union{..., Nothing}` relaxation** — port to unified; no reason unified should require a
   non-nothing `q`. Low risk.
5. **Pol_type FT coercion in rt_run** — port as a safety check; after merge, unified's model
   carries FT parametrically so this may become a no-op, but keep the defensive cast.

---

## 2. Commit `d75dacb` — Per-direction CPU-staging for interaction workspace

### Summary

Evolves `InteractionWorkspace` to reduce peak GPU memory by ~3.5 GB (FP64). Splits the six 4D
output buffers into two passes: Pass 1 (downwelling — computes `ieR⁻⁺, ieT⁻⁻, ieJ₀⁻`) uses 3 GPU
buffers, results are synchronously `copyto!`'d to CPU arrays, GPU buffers are zeroed and reused
for Pass 2 (upwelling — `ieT⁺⁺, ieR⁺⁻, ieJ₀⁺`). At end, all 6 CPU-staged results are copied back
to the composite_layer GPU fields. Verified bit-exact; reports ~4.2% runtime overhead from extra
D2H/H2D transfers.

### Abstractions introduced

- **Revised `mutable struct InteractionWorkspace`** with fields:
  `tmp_inv, tmpR⁻⁺, tmpR⁺⁻, tmpT⁻⁻, tmpT⁺⁺, tmpJ₀⁻, tmpJ₀⁺` (3D GPU, always present);
  `gpu_ie_mat_A, gpu_ie_mat_B` (2 × 4D matrix GPU, reused across passes);
  `gpu_ie_src` (1 × 4D source GPU, reused);
  `cpu_ieR⁻⁺, cpu_ieT⁻⁻, cpu_ieJ₀⁻, cpu_ieR⁺⁻, cpu_ieT⁺⁺, cpu_ieJ₀⁺` (6 × `Array` CPU staging);
  `staged::Bool` flag.
- **Constructor keyword** `InteractionWorkspace(composite_layer, added_layer; staged::Bool=false)`.
- **Two code paths** inside `ScatteringInterface_11::RRS` method: when `staged=true`, buffer
  aliasing (`tmpieR⁻⁺ = ws.gpu_ie_mat_A`, etc.) is used, with a `synchronize_if_gpu(); copyto!()`
  barrier between Pass 1 and Pass 2.
- Callers in `rt_run.jl`, `rt_run_bck.jl`, `rt_run_lin.jl` default to `staged=true`.

### Files & line ranges (at `d75dacb`)

- `src/CoreRT/CoreKernel/interaction_inelastic.jl`:
  - Lines ≈10–60 (struct + ctor + reset!) — ~123 insertions/33 deletions total.
  - Lines ≈328–400 (Pass 1 block in `ScatteringInterface_11::RRS`).
  - Lines ≈415–430 (inter-pass staging copy).
  - Lines ≈488–510 (final CPU→GPU copy of all 6 fields).
  - Note: the `VS_0to1_plus/VS_1to0_plus` variant of `ScatteringInterface_11` was **not** updated
    in this commit — only the pure-`RRS` path is staged. Open question whether VS+ needs the same
    treatment for large-nRaman runs.
- `src/CoreRT/rt_run.jl`, `rt_run_bck.jl`, `rt_run_lin.jl` — changed to `staged=true`.

### Unified counterpart

None. Unified's `RTWorkspace` is elastic-only and has no CPU staging story. Christian's
`InelasticWorkspace` proposal in the merge plan §Phase 2a does not include CPU staging; it lists
6 GPU 4D buffers (`tmpieR⁻⁺, tmpieR⁺⁻, tmpieT⁺⁺, tmpieT⁻⁻, tmpieJ₀⁺, tmpieJ₀⁻`) — matching
sanghavi's *earlier* 854b44c design, not the staged 3-buffer design.

### Recommendation

**Adopt as design input; re-examine whether staging is the right solution on unified.** The 4.2%
wall-clock overhead is small but the memory win (~3.5 GB FP64 on typical Raman runs) matters on
16 GB GPUs. Two reasons to treat it as design input rather than port verbatim:

- The staging mutex between the two halves of `ScatteringInterface_11` forces a hard sync that
  fights Phase 3 (flat 3D batching) of Christian's plan. If Phase 3 lands, nRaman collapses into
  a single `batched_mul!` and the 4D layout may shrink enough that GPU-only fits again.
- The sanghavi staged code has two copies of the logic (staged vs non-staged) inside one 150-line
  function. Unified conventions would prefer a smaller, clearer boundary — e.g. an
  `InelasticWorkspace` with `on_gpu::Bool` backing for `ie*` fields, and a `finalize!` that
  copies CPU→GPU at once.

Decision for the user: **do we merge staging as-is, or skip staging and rely on Phase 3 batching
to fix the memory pressure?** If the user chooses to skip staging, we lose a bit-exact 3.5 GB
headroom *today* but save ~150 lines of branchy kernel code.

**Recommended path:** include `staged` as an opt-in flag on the new unified `InelasticWorkspace`
(default `false` on machines with sufficient GPU memory), keep the hot path single-branched, and
let Phase 3 address memory on smaller GPUs.

---

## 3. Commit `3e9926f` — Raman optimization test infrastructure

### Summary

Purely additive: a single-band O2 A-band YAML config, a baseline-capture JLD2 script, a
regression-comparison script, and `.gitignore`'d output directory. No source-tree changes. This
is the harness that will let the merge be verified bit-exact against sanghavi's reference.

### Abstractions introduced

- **YAML config** `test/test_parameters/O2_parameters2_1band_opttest.yaml` — 1-band,
  GaussLegQuad, Stokes_IQU, `l_trunc=5`, `max_m=1`, `Architectures.GPU()`, 34-layer T/p/q
  profile, one LUT (`/net/fluo/data1/ABSCO_CS_Database/v5.2_final/o2_v52_v2.jld2`). Fast enough
  for iterative testing but exercises the full RRS pipeline.
- **`test/benchmarks/raman_optimization_baseline.jl`** — runs one RRS and one noRS computation,
  saves `(R, T, ieR_SFI, ieT_SFI)` + timing/allocation stats to
  `raman_opttest_output/baseline.jld2`. Hard-coded to project root
  `/home/sanghavi/code/github/vSmartMOM.jl`.
- **`test/benchmarks/raman_optimization_compare.jl`** — loads the baseline, reruns, and checks
  bit-exact equality of R, T, ieR, ieT with reporting on max abs diff.
- **`test/benchmarks/raman_opttest_output/.gitignore`** — ignores `*.jld2`, `*.log`.

### Files & line ranges

- `test/benchmarks/raman_optimization_baseline.jl` — new, 161 lines.
- `test/benchmarks/raman_optimization_compare.jl` — new, 182 lines.
- `test/test_parameters/O2_parameters2_1band_opttest.yaml` — new, 73 lines.
- `test/benchmarks/raman_opttest_output/.gitignore` — new, 3 lines.

### Unified counterpart

None. Unified has its own 474-test suite via `test/runtests.jl` with helpers
`test/test_helpers.jl`, but no Raman-specific GPU baseline capture.

### Recommendation

**Adopt wholesale, but rewrite the scripts against unified's RTModel API.** The YAML can be
dropped in as-is to `test/test_parameters/` on unified. The two driver scripts use the old
sanghavi API (`vSmartMOM_Model`, `model.τ_abs`, etc.) and need to be updated — which is exactly
the kind of "example scripts from the branch owner" that Christian's Phase 1c of the merge plan
asks for. This ticket doubles as both a migration example and a regression gate.

The hard-coded project root (`/home/sanghavi/code/github/vSmartMOM.jl`) must be rewritten to
accept a command-line argument or use `dirname(@__DIR__)`.

---

## 4. Commit `e89ec1c` — Update baseline script and add batched ops benchmark

### Summary

Simplified the baseline capture (single RRS + single noRS run), and added a 365-line
micro-benchmark that proves `NNlib.batched_mul!` is **5.5× slower** than the allocating `⊠` form
at Raman's characteristic shapes (NquadN=15 × NquadN=15 × nSpec≈5500 matrices). This empirical
result **invalidates** the naive "just make it in-place" rewrite: for 15×15 batched GEMM on GPU,
CUDA's memory pool handles small-allocation reuse better than `batched_mul!` handles fixed-stride
work decomposition. The commit message calls out that the optimization path must instead be
"custom GPU kernels to parallelize the sequential Δn loop."

Baseline numbers reported at this commit: n_Raman=172, nSpec=5424, NquadN=15, RRS wall time 1672 s,
noRS wall time 2.8 s — **589× slower** with Raman.

### Abstractions introduced

None in `src/`. Benchmark helpers only:

- Constants in `test/benchmarks/raman_batched_ops_benchmark.jl`: `NQ=9, NSPEC=5500, NRAMAN=963,
  NL=5000`.
- Test functions `test_allocating_expr`, similar variants for in-place and view-based ops.

### Files & line ranges

- `test/benchmarks/raman_batched_ops_benchmark.jl` — new, 365 lines.
- `test/benchmarks/raman_optimization_baseline.jl` — rewrite, −137/+181.

### Unified counterpart

Possibly relevant: `a4e4187` on unified ("Add batched-kernel and Raman scaling benchmarks with
writeup"). Its existence suggests unified has already begun a similar empirical probe; worth
reading its contents and reconciling.

### Recommendation

**Preserve the benchmark findings** as guidance for Phase 3 of Christian's merge plan. Port the
micro-benchmark script (`raman_batched_ops_benchmark.jl`) to the unified tree at
`test/benchmarks/` — its result reshapes the whole GPU optimization strategy. Specifically:

- **Christian's Phase 3 proposal ("gather + single batched_mul")** should be validated against
  this result before implementation. Sanghavi's data says the naive `batched_mul!` swap loses;
  the right move is likely either custom `@kernel` ops (see user's memory: *"Custom @kernel 8x
  slower than CUBLAS for 15x15; don't replace ⊠"*) OR collapsing over Δn so CUBLAS batches over a
  dimension that's far bigger than 15.
- The "589× Raman slowdown vs noRS" number is a useful top-line regression target.

Flag for the user: **unified's own `a4e4187` benchmark needs to be cross-checked against
sanghavi's numbers before the Phase 3 direction is locked in.** (That cross-check is Inventory
B's remit.)

---

## 5. Commit `083353b` — Update lin+Raman (forward-Raman physics updates in scope)

### Summary

Large commit (+1484/−522, 22 files) mixing forward-path physics fixes, inelastic library
cleanups, and linearized Raman. **Only the forward-path physics and helper changes are in scope
for this inventory** (per user: "linearized Raman is OUT OF SCOPE, but forward updates are in
scope"). The in-scope changes are:

1. **Rayleigh cross-section formula switched to Bodhaine et al. 1999 Eq. 30**, replacing the
   older `0.00864 × (psurf/1013.25) × λ^(…)` fit. The Bodhaine form is more accurate,
   particularly in the UV/blue, and has an implicit depolarization factor that is recomputed
   explicitly here for flexibility.
2. **New `reduce_profile` using linear interpolation on uniform pressure half-levels**; the old
   bin-averaging method is preserved as `reduce_profile_old`.
3. **`compute_ϖ_Cabannes` rewritten** to use the unambiguous formula
   `ϖ_Cabannes = 1 - σ_RRS / σ_Rayl` (or `1 - (σ_VRS+σ_RVRS)/σ_Rayl` for VS). The commit comment
   notes that the previous `σ_elastic` was *incorrectly* assumed to be the purely-elastic
   (Cabannes) cross-section; it was actually the total Rayleigh cross-section. This is a physics
   bug fix.
4. **`model_from_parameters` now computes band-specific effective T-averaged N₂/O₂ constants**
   via `InelasticScattering.getRamanAtmoConstants(νₘ, effT)` where
   `effT = (profile.vcd_dry' * profile.T) / sum(profile.vcd_dry)`, and calls new
   `compute_γ_air_Cabannes!` + `compute_γ_air_Rayleigh!` accessors separately. Replaces the old
   single-call `compute_γ_air_Rayleigh!(λₘ)` that returned both ϖ_Cabannes and two γ values from
   the same invocation.
5. **Constant rename `nm_per_m → nm_per_cm`** in `InelasticScattering` module (and related files).
   This is a **units correction**: the numeric value `1.0e7` converts cm⁻¹↔nm, not m⁻¹↔nm, so the
   old name was misleading. Must be paired with the same rename on unified or the two codebases
   won't link.
6. **`apply_gridlines!` allocation cleanups** in `stellar_inelastic_helper.jl`: uses
   `grid_in_collected = collect(grid_in)` once at the top, computes `nz_mask = abs.(xin) .> 0`
   once, and uses `findall(>(0), σ_tmp)` instead of `findall(x -> x in σ_tmp[σ_tmp.>0], σ_tmp)`.
7. **`compute_energy_levels!` performance rewrite** in `inelastic_cross_section.jl`: precomputes
   `E₁_pow` and `E₂_pow` powers instead of recomputing `E₁^(l-1)` and `(v+0.5)^(k-1)` inside the
   hot nested loop. Wraps loop in `@inbounds`. No physics change; 2×–3× speedup expected.
8. **`compute_effective_coefficents!` physics fix**: the formula `α̅(2πcν, T) = α̅₀₀(1 + …)/(1 −
   (2πcν/ω₀)²)` is changed to `… /(1 − (cν/ω₀)²)`. The old `2π` factor was wrong (cross-multiplies
   incorrectly with the way ν is stored as cm⁻¹). Physics fix.
9. **`AbstractArray → AbstractVector/AbstractMatrix` type-firming** across `raman_constants.jl`,
   and `inelastic_helper_old.jl` (+924 lines) added as a frozen copy of the old helper for
   reference. This is a readability / dispatch-stability improvement.
10. **Absorption LUT out-of-range clamping + range-masking**: `compute_absorption_cross_section`
    now detects grid points outside the LUT ν-range and returns zero there instead of
    extrapolating / erroring.
11. **Renamed `q::AbstractArray{FT}` → `q::Union{AbstractArray{FT}, Nothing}`** — partial overlap
    with 854b44c.
12. **`atmo_prof.jl`** `cat(..., dims=(1))` replaced with `vcat(...)` for `i_λ₁λ₀_all`
    construction — clarity only.

Out of scope: `src/CoreRT/atmo_prof_lin.jl`, `src/CoreRT/CoreKernel/elemental_lin.jl`,
`src/CoreRT/lin_model_from_parameters.jl`, `src/CoreRT/CoreKernel/doubling_lin.jl`,
`src/CoreRT/CoreKernel/interaction_lin.jl`, `src/CoreRT/CoreKernel/rt_kernel_lin.jl` and any
`_lin` test configs — these are linearized Raman and must not be carried over.

### Abstractions introduced (in scope)

- Function `reduce_profile_old(n, profile)` — kept as escape hatch.
- Function `reduce_profile(n, profile)` — new linear-interpolation implementation.
- New `InelasticScattering` exports: `compute_γ_air_Cabannes!`, `compute_γ_air_Rayleigh!`.
- Local helper `_interp(data)` inside `reduce_profile`.

### Files & line ranges (in-scope subset)

- `src/CoreRT/atmo_prof.jl` — lines ≈94–220 (`reduce_profile_old` + new `reduce_profile` + new
  `getRayleighLayerOptProp` Bodhaine formula, ~+90 lines).
- `src/CoreRT/model_from_parameters.jl` — lines ≈140–180 and ≈840–890 (in-scope) — replace
  single-call γ_air with split Cabannes/Rayleigh.
- `src/Inelastic/inelastic_helper.jl` — ~500-line refactor; in-scope portions: rename constant,
  rewrite `compute_ϖ_Cabannes` (Lines ≈70–180), bug fix at `σ_RRS` assembly (uses `o2.effCoeff`
  not `n2.effCoeff` for O₂ lines).
- `src/Inelastic/raman_atmo_prop.jl` — `end:-1:1` → `reverse(…)` and matching `cat → vcat`.
- `src/Inelastic/InelasticScattering.jl` — add exports (line ≈29).
- `src/Inelastic/src/inelastic_cross_section.jl` — `α̅` formula 2π fix, `compute_energy_levels!`
  loop reorder.
- `src/Inelastic/src/raman_constants.jl` — type-firming of struct fields.
- `src/Inelastic/stellar_inelastic_helper.jl` — `apply_gridlines!` allocation cleanup (in-scope
  for stellar forward).
- `src/Inelastic/inelastic_helper_old.jl` — new file, 924 lines, frozen reference copy. Consider
  deleting after merge.
- `src/Absorption/compute_absorption_cross_section.jl` — range-masking for out-of-LUT ν
  (≈line 160–170).
- `src/Absorption/constants/constants.jl` — minor constant rename.
- `src/Absorption/make_model_helpers.jl` — minor.

### Unified counterpart

- `src/CoreRT/tools/atmo_prof.jl` on unified still uses the **OLD** Rayleigh formula (`tau_scat =
  0.00864 * …`) at line 175, and the old bin-averaging `reduce_profile` at line 96.
- `src/Inelastic/inelastic_helper.jl` on unified still has `nm_per_m` (not `nm_per_cm`) and the
  old `compute_ϖ_Cabannes` formula. Christian's merge plan claims unified "already ported
  sanghavi forward physics (commit 59f8de8)" — that claim is **outdated**; these specific physics
  fixes are newer than 59f8de8.
- `src/Absorption/compute_absorption_cross_section.jl` on unified likely does not have the
  range-masking (needs a diff check).
- The `α̅` 2π bug fix and `compute_energy_levels!` power rewrite are almost certainly not on
  unified (they are in `src/Inelastic/src/inelastic_cross_section.jl` which did not appear in
  59f8de8's port list).

### Recommendation

**Per authority rule (2026-04-19):** Items 1, 3, 4, 5, 6, 7, 8, 9, 10 (and item 2 per amendments §2.1 with bin-averaging kept as `reduce_profile_binavg`) **land via the wholesale Inelastic port** in Phase 1a / 1b per `plans/IMPLEMENTATION_PLAN_v2.md`. No per-item evaluation is required.

Carve-outs from the authority rule:

- **Item 2 (new `reduce_profile`)** — default flips to linear interpolation per amendments §2.1 / Phase 1a. Bin-averaging preserved as `reduce_profile_binavg` via opt-in keyword. Lands in Phase 1a alongside the Bodhaine Rayleigh switch as a single commit with 474-test re-baselining.
- **Item 8 (α̅ `2π` drop)** — verified correct by the user per amendments §2.4. Port with the verification comment specified there.
- **`inelastic_helper_old.jl` (924-line frozen reference)** — does not port per amendments §1 exception; drop is achieved by not-porting, not by a delete commit.
- **_lin files** — do **not** port (linearized Raman permanently out of scope per amendments §7). Flag any cross-contamination at Phase 1b port site.

---

## 6. Cross-comparison: sanghavi `InteractionWorkspace` vs Christian's proposed `InelasticWorkspace`

### Field-by-field

| Purpose | Sanghavi has (post-`d75dacb`) | Christian's plan proposed (§2a) | Recommendation | Why |
|---|---|---|---|---|
| `(I − R⁺⁻r⁻⁺)⁻¹` | `tmp_inv :: AT3` | `tmp_inv :: AT3` | **Unify with `RTWorkspace.tmp_inv`** (unified already has it) | No reason to allocate twice |
| `T_inv` (T01_inv/T21_inv staging) | Not a separate field; computed inline | `T_inv :: AT3` | **Adopt Christian's `T_inv`** | Sanghavi code computes into temporaries that become workspace fields anyway; giving it a name is cleaner |
| 3D R/T/J scratch | `tmpR⁻⁺, tmpR⁺⁻, tmpT⁻⁻, tmpT⁺⁺, tmpJ₀⁺, tmpJ₀⁻` (6 × AT3) | `tmpR⁻⁺, tmpR⁺⁻, tmpT⁺⁺, tmpT⁻⁻, tmpJ₀⁺, tmpJ₀⁻` (6 × AT3) | **Identical; adopt either naming** (but use unified's lowercase/uppercase convention) | Same struct, different spellings |
| 4D ie* matrix buffers | `gpu_ie_mat_A, gpu_ie_mat_B` (2 × AT4 reused) | `tmpieR⁻⁺, tmpieR⁺⁻, tmpieT⁺⁺, tmpieT⁻⁻` (4 × AT4 distinct) | **Keep distinct fields** (Christian's layout) unless memory pressure forces staging | Aliasing 4 named buffers onto 2 physical buffers is what forces sanghavi's hard sync between passes; on ≥24 GB GPUs this is unnecessary |
| 4D ie* source buffers | `gpu_ie_src` (1 × AT4 reused for J₀⁻ then J₀⁺) | `tmpieJ₀⁺, tmpieJ₀⁻` (2 × AT4) | **Keep distinct fields** (Christian's layout) | Same reasoning as above |
| CPU staging buffers | `cpu_ieR⁻⁺, cpu_ieT⁻⁻, cpu_ieJ₀⁻, cpu_ieR⁺⁻, cpu_ieT⁺⁺, cpu_ieJ₀⁺` + `staged::Bool` | Not proposed | **Add as optional feature** (`staged::Bool = false` flag; lazy-alloc the 6 CPU Arrays only when `staged=true`) | Needed to hit 16 GB GPUs with large nRaman; but Phase 3 batching may obviate the need |
| `buf3d_a, buf3d_b, buf3d_c` reusable 3D bufs | Not present | Yes, proposed | **Adopt Christian's proposal**; share with elastic `RTWorkspace.tmp3d_a/tmp3d_b` if possible | These are the `⊠`-result staging bufs inside `doubling_inelastic.jl`; sanghavi's code allocates them fresh every call |
| `batch_inv!` pivot/info | Not present in `InteractionWorkspace` (allocated inside `batch_inv!`) | Proposed: share with elastic workspace | **Adopt proposal; share with `RTWorkspace.pivot/info`** | Sanghavi's runs pay ~2000 pivot allocations per run (see merge plan audit table); massive win |
| Doubling-specific bufs `gp_refl, tt_gp, J₁⁺, J₁⁻, ieJ₁⁺, ieJ₁⁻` | Not in `InteractionWorkspace` (`doubling_inelastic.jl` still `similar()`s them per layer) | Yes, proposed | **Adopt proposal** | Doubling is called once per layer; currently allocates 6 fresh arrays per layer call, totalling ~3600 GPU allocs/run per merge plan audit |

### Summary of the reconciliation

Christian's proposed `InelasticWorkspace` covers a **superset** of sanghavi's `InteractionWorkspace`:

- ✓ Both cover the 3D interaction scratch (`tmp_inv`, `tmpR*`, `tmpT*`, `tmpJ*`).
- ✓ Both cover the 4D interaction scratch (`tmpie*`) — but sanghavi aliases to save 3.5 GB.
- ✗ Only Christian covers the doubling scratch (`gp_refl, tt_gp, J₁⁺/⁻, ieJ₁⁺/⁻`).
- ✗ Only Christian covers the `batch_inv!` pivot/info.
- ✗ Only Sanghavi has CPU staging.

**Recommendation:** Build `InelasticWorkspace` on `sanghavi-unified` to Christian's layout plus
two additions from sanghavi:
1. An optional `staged::Bool` flag with lazy-allocated CPU buffers (default `false`).
2. Reuse sanghavi's `reset!` pattern with `@inline`-friendly `.=0` fills.

Name and naming conventions: use unified's uppercase-for-composite / lowercase-for-added-layer
convention, and unified's `(; field) = struct` unpacking (not `@unpack`).

---

## 7. Other notable commits (scan of ~50 unique sanghavi commits)

The inventory above covers the 5 commits called out by the user. Beyond those, from
`git log unified-vsmartmom..sanghavi --oneline` (51 commits), the following are load-bearing for
the merge:

| Commit | Title | Why it matters | Recommendation |
|---|---|---|---|
| `9a26002` | update Raman, EMIT, Misc. | **Largest recent commit** (+4035/−97, 36 files). Contains (a) `hem_R, hem_T` hemispheric-integrated radiance outputs added to `rt_run` return signature and `postprocessing_vza!`; (b) `apply_D_matrix_elemental!` Stokes_I (scalar) shortcut avoiding GPU kernel dispatch when `n_stokes==1`; (c) EMIT benchmark scripts (`emit_modtran_noRS_scenarios.jl` 804 lines, `compare_rt_EMIT.jl` 350 lines, `create_HITRAN_LUTs.jl`); (d) Float32 test config `O2_parameters2_SIF_grid_float32.yaml`; (e) many `*.log` artifacts in `raman_opttest_output/`. | **Needs own line-item ticket.** Treat as a 4th optimization commit. The API change `rt_run` → returns 6-tuple not 4-tuple (`R_SFI, T_SFI, ieR_SFI, ieT_SFI, hem_R, hem_T`) is a **breaking change** to any caller. Must coordinate with unified's `rt_run.jl` signature. The `n_stokes==1` scalar shortcut is a clear perf win; port wholesale. EMIT scripts → `examples/` per Christian's Phase 1b. Delete `*.log` artifacts. |
| `ad16041` | add hemispherically integrated radiance | 1-line fix in `emit_modtran_noRS_scenarios.jl`; depends on `9a26002`. | Follows `9a26002`; port together. |
| `22b425d` | adapt T-p boundaries for absorption XS | 6-line `clamp(pressure, …); clamp(temperature, …)` in `compute_absorption_cross_section`. Physics-correctness for high-altitude layers where profile p/T fall outside LUT. | **Adopt wholesale**; low-risk numerical safety. |
| `8419745` | Ignore .jld2/.jld files | `.gitignore` update — trivial but important for the merge (large binary JLD2s should not be committed). | Adopt. |
| `9ee9a75` | add plan | Just a plans file. | Ignore for code merge; keep for archival. |

Remaining commits (dcea2a6, da7c321, 44f881c, 080d5e9, 50cb910, 0c638b9) are all linearized RT
work — **out of scope** per user direction. None of these should be ported to `sanghavi-unified`.

Older application/research commits (b6d8fba, 70c0c4d, 87467a1, f40eb59, cb3ad9b, 41e01a5,
118eb65, 3a77b99, dc0bf7c, bbb577a, aae70e0, df569bf, e745c35, b0de496, e366837, 5616518,
a549a13, e7eef97, cd2de96, 2b9a0a4, 7376422, 298c63c, 2cda7ed, e61b318, 2fc2ae2, 89e34c5,
190967c, 77497d5, 1e154b9, d6afa16, a04d2a2, 295cf47, cc7b7b7, 3dd423c) are either:
- Application-code for specific papers (OCO, EMIT, Balsamic, SIF, Raman paper figures) → belong
  in `examples/` per Christian's Phase 1b, not in the library merge.
- Single-scatter approximation (295cf47, a04d2a2, d6afa16) → Christian's plan already calls this
  out in Phase 1b; non-trivial library code to port, but **already on Christian's radar**.
- Incremental Raman debug/improve commits that are likely already in unified's port
  (59f8de8 merge).

---

## 8. Open questions for user validation

**Resolution status as of 2026-04-19 user review session:** Q1, Q4, Q5, Q6 are closed by `plans/PLAN_AMENDMENTS_2026-04-19.md`. Remaining questions (Q2, Q3, Q7, Q8, Q9, Q10) are still open at the noted severity.

1. ~~**Staging (commit `d75dacb`): Keep or drop?**~~ **Resolved (amendments §4 Phase 4):** Keep. `staged::Bool = true` **default-on**. Phase 5 may eventually replace the memory win with batching if re-measurement shows benefit; until then, don't trade a certain 3.5 GB win for an uncertain one.

2. **Scope of `9a26002`:** The `hem_R, hem_T` additions change `rt_run`'s return signature from
   4-tuple to 6-tuple (SFI=true) and 2-tuple to 4-tuple (SFI=false). On unified, `rt_run(model)`
   returns `(R, T)`. Do we:
   (a) also add `hem_R, hem_T` to unified's return (breaking current callers), or
   (b) add a keyword `return_hemispheric::Bool=false` to opt in, or
   (c) expose `hem_R, hem_T` via a new accessor that reads from the composite layer?
   Decision affects every downstream caller and the unified test suite.

3. **`n_stokes == 1` scalar shortcut in `apply_D_matrix_elemental!`** (`9a26002`): The shortcut
   does `r⁺⁻[:] = r⁻⁺; t⁻⁻[:] = t⁺⁺; return` without launching the KernelAbstractions kernel.
   This assumes the scalar case has no sign flips — verify against unified's canopy/RAMI test
   cases, which use scalar Stokes_I. Does it still match `apply_D_matrix!` in `doubling.jl`
   mod1 logic (commit `a5e0de5`)?

4. ~~**Rayleigh formula switch** (Bodhaine, item 1 of commit `083353b`)~~ **Resolved (amendments §4 Phase 1a):** Bodhaine 1999 Eq. 30 becomes the default on `sanghavi-unified`. 474-test suite is re-baselined against the new Bodhaine values as one logical commit, paired with the `reduce_profile` default switch (Q5).

5. ~~**`reduce_profile` semantics** (item 2 of `083353b`)~~ **Resolved (amendments §2.1 / §4 Phase 1a):** Linear interpolation on uniform pressure half-levels becomes the default. Bin-averaging preserved as `reduce_profile_binavg` (or equivalent) available via an opt-in keyword. Re-baselining lands in the same commit as the Bodhaine switch.

6. ~~**Delete `src/Inelastic/inelastic_helper_old.jl`**~~ **Resolved (amendments §1 exception):** Not ported. Drop is achieved by not-porting, not by a delete commit; the file stays in git history on the sanghavi branch.

7. **Sanghavi's baseline numbers** (e89ec1c commit message): `RRS: 1672s, noRS: 2.8s` on n_Raman=172,
   nSpec=5424, NquadN=15. Inventory B is supposed to compare unified's current numbers against
   these. Is that captured elsewhere, or should Inventory A flag this as a regression target?

8. **Delete `raman_opttest_output/*.log` artifacts** (added in `9a26002`, ~14 log files, ~2000
   lines total)? They should be gitignored, not committed.

9. **`q::Union{…, Nothing}` relaxation**: The `AtmosphericProfile.q` field becomes optional.
   Does unified's IO layer (YAML parsing, NetCDF readers) handle `q === nothing` gracefully, or
   does this need defensive code in `AtmosProfile.jl`?

10. **Const rename `nm_per_m → nm_per_cm`**: Will break any user script that directly referenced
    `InelasticScattering.nm_per_m`. Breaking change; needs release note.

---

*Inventory A ends here. Inventory B (tolerance and benchmark reconciliation) is out of scope for
this document.*
