# Inventory D — Physics / Semantic Delta Summary

**As-of:** 2026-04-19
**Sanghavi worktree:** `/home/sanghavi/code/github/vSmartMOM.jl`
**Sanghavi branch tip:** `9ee9a75` (branch `sanghavi`) — "add plan" (2026-04-18)
**Unified worktree:** `/home/sanghavi/code/github/uni_vSmartMOM`
**Unified branch tip:** `a4e4187` (branch `unified-vsmartmom`) — "Add batched-kernel and Raman scaling benchmarks with writeup"
**Reference merge commit (ported portion):** `59f8de8` "Port sanghavi branch physics fixes to forward kernels" (elemental.jl, elemental_lin.jl, rt_kernel.jl, rt_run.jl — `mod`→`mod1`, off-diagonal L'Hôpital limit, wavelength-dependent F₀).

Scope: non-optimization physics differences — bug fixes, new dispatches, missing functions, new numerical formulas — between `sanghavi` and `unified-vsmartmom`. Out of scope: Inventory A (performance/workspace), Inventory C (SIF plumbing), linearized Raman (dropped).

**Editorial-pass note (2026-04-19):** The 2026-04-19 user review session established the **authority rule** — *sanghavi is the authority for the inelastic path* — which settles by fiat the per-item P0/P1/P2 framing in §6 below. §7 open questions Q1, Q2, Q5 are closed by `plans/PLAN_AMENDMENTS_2026-04-19.md`. Specific redirections are marked inline. Evidence and factual observations stand.

---

## 1. Inelastic CoreKernel files

### 1.1 `src/CoreRT/CoreKernel/elemental_inelastic.jl`

| File | First-on | What differs | Assessment |
|------|----------|--------------|------------|
| `elemental_inelastic.jl` | sanghavi | `apply_D_matrix_elemental!` for `Union{RRS,RRS_plus}` **and** `Union{VS_0to1_plus,VS_1to0_plus}` has an **early-return scalar branch** `if n_stokes == 1 … ier⁺⁻[:] = ier⁻⁺; iet⁻⁻[:] = iet⁺⁺; return nothing; end`. Unified lacks this scalar shortcut. | **Correctness-preserving optimization with physics-safety guard.** The D-matrix for scalar (nStokes=1) is identity, so `ier⁺⁻ = ier⁻⁺` and `iet⁻⁻ = iet⁺⁺`. Without the early return, unified launches a kernel whose body depends on the D-matrix; if the kernel code path for `n_stokes==1` is right, results match. Low-risk-to-port; medium-value. |
| `elemental_inelastic.jl` | equal (semantically) | Sanghavi uses inline `wct02 = m == 0 ? FT(0.50) : FT(0.25)`; unified factored into `fourier_weight`/`scaled_weights` helpers in `rt_helpers.jl`. | Equivalent. No action. |
| `elemental_inelastic.jl` | unified | Unified adds comprehensive docstrings citing Sanghavi & Frankenberg 2023, JQSRT 311, 108791, Eq. 14. Sanghavi has only internal `#Suniti:` comments. | Doc-only. |
| `elemental_inelastic.jl` | sanghavi | Many `@show`/debug comments on sanghavi; unified cleaned. | No action. |

No formula differences in `get_elem_rt_RRS!`, `get_elem_rt_SFI_RRS!`, `get_elem_rt_VS!`, `get_elem_rt_SFI_VS!`.

### 1.2 `src/CoreRT/CoreKernel/elemental_inelastic_plus.jl`

Same pattern as 1.1: sanghavi has leftover `@kernel function get_elem_rt_RRS!` commented-out duplicate; unified cleaned. Numerics identical.

### 1.3 `src/CoreRT/CoreKernel/doubling_inelastic.jl`

| File | First-on | What differs | Assessment |
|------|----------|--------------|------------|
| `doubling_inelastic.jl` | equal | Unified uses unified `j₀⁺`/`j₀⁻` (lowercase) for AddedLayer source terms per the branch-wide AddedLayer rename; sanghavi still uses `J₀⁺`/`J₀⁻` on AddedLayer. | Naming only. Unified has the rename consistent with `AddedLayer` struct. |
| `doubling_inelastic.jl` | equal | Unified docstrings reference Part II equations; sanghavi references "Raman paper draft" equation numbers. | Doc-only. |

No formula differences in `doubling_inelastic!`, `apply_D_matrix_elemental!`, `apply_D_matrix_SFI!`.

### 1.4 `src/CoreRT/CoreKernel/interaction_inelastic.jl`

| Item | First-on | What differs | Assessment |
|------|----------|--------------|------------|
| `composite_layer.R⁺⁻ .= tmpR⁺⁻` (was `tmpR⁻⁺`) | unified (`e0ab098`) and sanghavi (both have it) | R⁺⁻ bug-fix (Part II Eq. 19) in both `interaction_helper!(RRS,::ScatteringInterface_11,…)` and `interaction_helper!(VS_*,::ScatteringInterface_11,…)`. | **Present on both.** No action. |
| `InteractionWorkspace` + `staged` path | sanghavi | Per-direction CPU-staging workspace optimization. | **Inventory A scope** (optimization). Exclude here. |
| `j₀⁺`/`j₀⁻` rename | unified | Unified uses lowercase AddedLayer source fields. | Naming consistent with unified AddedLayer rename. |
| Docstrings citing Sanghavi & Frankenberg 2023, JQSRT 311, 108791 | unified | Documentation add. | Doc-only. |

No formula differences.

### 1.5 Missing `*_inelastic.jl` files

Both branches have `elemental_inelastic.jl`, `elemental_inelastic_plus.jl`, `doubling_inelastic.jl`, `interaction_inelastic.jl`. No missing files.

However:
- Unified has `elemental_canopy.jl`, `interaction_hdrf.jl`, `rt_helpers.jl` that sanghavi does not. These are canopy/hotspot BRDF additions on unified (noRS-only), not Raman physics.

---

## 2. `rt_kernel.jl` and `rt_kernel_lin.jl`

### 2.1 `rt_kernel.jl` — active dispatches

| Dispatch | sanghavi | unified | Status |
|---|---|---|---|
| `rt_kernel!(::noRS, …, computed_layer_properties)` (legacy precomputed) | yes | yes | equal (cosmetic) |
| `rt_kernel_canopy!(::noRS, …)` | no | **yes** (unified-only) | canopy addition in unified, no Raman coupling |
| `rt_kernel!(::Union{RRS,VS_0to1,VS_1to0}, …, computed_layer_properties)` (legacy precomputed) | yes | yes | equal |
| `rt_kernel!(::noRS, …, ::CoreScatteringOpticalProperties, scattering_interface, τ_sum, …)` (new hybrid path) | yes | yes | equal (unified refactored via `init_layer`, `get_dtau_ndoubl`) |
| `rt_kernel!(::Union{RRS,VS_0to1,VS_1to0,RRS_plus,VS_0to1_plus,VS_1to0_plus}, …, ::CoreScatteringOpticalProperties, …)` (hybrid, Raman+_plus union) | **yes** (single merged method) | **no** | **Unified has split into two separate dispatches** (RRS/VS and RRS_plus/VS_plus) |
| `rt_kernel!(::Union{RRS_plus,VS_0to1_plus,VS_1to0_plus}, …)` (dedicated _plus method) | commented out in `#= … =#` | **yes** (active) | Unified's dedicated _plus method is active; sanghavi merged it |

Both approaches are **numerically equivalent** but organized differently. Unified's split is cleaner; sanghavi's merged union is shorter. Neither is a physics bug. Pick one and retain.

`F₀` is now passed through all call sites on both branches (per `59f8de8`). No remaining forward-only leaks.

### 2.2 `rt_kernel_lin.jl`

Out of scope (linearized Raman). Not audited.

### 2.3 `rt_kernel_multisensor.jl` and `rt_kernel_ss.jl`

Both present on both branches. Only cosmetic/naming diffs (`J₀⁺`→`j₀⁺`, `Array()`→`collect()`, `@unpack`→`(; …)`). No physics differences.

### 2.4 `interaction_ss.jl`

Both present. Unified adds `atype()` wrapping of `composite_layer.ieJ₀⁺/ieJ₀⁻` before dispatch, and uses the new AddedLayer `j₀⁺`/`j₀⁻` naming. Cosmetic.

---

## 3. `src/Inelastic/`

### 3.1 Types (`types.jl`)

| Field / change | First-on | What differs | Assessment |
|---|---|---|---|
| `abstract type AbstractRamanType` (no FT param) | sanghavi | Unified: `AbstractRamanType{FT}` (parameterized) | Unified is better Julia; all concrete types subtype `AbstractRamanType{FT}`. Not a physics difference but interface-breaking. |
| `bandSpecLim = []` default | sanghavi | Unified: `bandSpecLim::Vector{UnitRange{Int}} = UnitRange{Int}[]` | Type-stable default. Not a physics difference. |
| `iBand::Array{Int,1} = [1]` | sanghavi | Unified: `iBand::Vector{Int} = Int[1]` | Cosmetic. |
| `noRS` `ϖ_Cabannes` default | sanghavi: `[1.0]` (length 1) | unified: `[1.0, 1.0, 1.0]` (length 3) | **Semantic**: unified assumes 3 bands by default, sanghavi 1. Impacts single-band runs where index >1 is accessed. Sanghavi's single-element default was a source of bugs for multi-band. |
| `F₀ = []` on sanghavi `noRS_plus` | sanghavi | unified: `F₀::Array{FT,2} = zeros(FT, 1, 1)` | Type-stable default. Unified is safer. |
| `grid_in::Array{StepRangeLen{FT},1}` | sanghavi | unified: `grid_in::Vector{AbstractRange{Float64}} = AbstractRange{Float64}[]` (for `VS_0to1_plus`, `VS_1to0_plus`, `sol_VS_0to1_plus`, `sol_VS_1to0_plus`) | Cosmetic. |

**All docstrings on unified types** (cite Sanghavi & Frankenberg 2023) — doc-only on unified.

No new physics fields on sanghavi. Both have `F₀`, `SIF₀`, `ϖ_Cabannes`, `ϖ_λ₁λ₀`, `i_λ₁λ₀`, `Z⁺⁺_λ₁λ₀`, `Z⁻⁺_λ₁λ₀`, etc.

### 3.2 `inelastic_helper.jl` — IMPORTANT physics-adjacent differences

| Function | First-on | What differs | Assessment |
|---|---|---|---|
| `compute_ϖ_Cabannes(RS_type::Union{RRS, RRS_plus}, λ₀)` | sanghavi (renamed/re-derived) | Sanghavi: `σ_Rayl = Σ(vmr_i * σ_Rayl_coeff_i) * ν₀⁴; ϖ_Cabannes = 1.0 - σ_RRS/σ_Rayl`. Unified: `σ_elastic = Σ(vmr_i * σ_Rayl_coeff_i) * ν₀⁴; ϖ_Cabannes = σ_elastic/(σ_RRS+σ_elastic)`. | **Numerically different if `σ_Rayl_coeff` is the _Rayleigh_ (Cabannes+RRS) cross-section (as sanghavi comment claims) rather than pure-elastic Cabannes.** Sanghavi docstring: `σ_Rayl_coeff` stores `128π⁵α̅² (1 + 2γ_C)/(3 − 4γ_C)` which is the FULL Rayleigh. Then `σ_elastic = σ_Rayl − σ_RRS`, and `σ_elastic/σ_Rayl = 1 − σ_RRS/σ_Rayl` (sanghavi) ≠ `σ_elastic/(σ_elastic+σ_RRS)` (unified, with `σ_elastic` labelled but computed as `σ_Rayl`). **Sanghavi's is the corrected physics**; unified's is the pre-fix code with wrong labelling. This pattern repeats in the VS and _plus variants. |
| `compute_ϖ_Cabannes(RS_type::Union{VS_0to1_plus, VS_1to0_plus}, λ₀)` | sanghavi | Same issue: sanghavi uses `ϖ_Cabannes_VS = 1.0 − (σ_RVRS+σ_VRS)/σ_Rayl`; unified uses `σ_elastic/(σ_RVRS+σ_VRS+σ_elastic)` | Same physics fix needed. |
| `compute_ϖ_Cabannes(RS_type::Union{RRS,VS,…}, depol, λ₀)` (3-arg depolarization variant) | sanghavi has it commented out with a debug version that returns `(4-3*depol)/(4+2*depol)`; unified has an active version returning the depolarization formula | — | Both retain the 3-arg variant as dead/debug; active callers use the 2-arg `(RS_type, λ₀)` signature. |
| `compute_γ_air_Rayleigh!(λ₀, RS_type)` and `compute_γ_air_Rayleigh!(λ₀, n2, o2)` | **sanghavi** | Sanghavi rewrote this: directly computes `γ_air_Rayleigh` from N₂/O₂ polarizabilities using a molecule-weighted combination formula with Rayleigh cross-sections (returns `(γ_air_Rayleigh, σ_air_Rayleigh)`). **Unified**: `compute_γ_air_Rayleigh!(λ₀)` internally calls `getRamanAtmoConstants(ν̃,300.0)` (hard-coded T=300K) and derives γ from `γ_air_Cabannes` and `ϖ_Cabannes`. | **Sanghavi's version is more accurate and returns the air-Rayleigh cross-section as well.** Unified uses a hard-coded 300 K effective T. |
| `compute_γ_air_Cabannes!(λ₀, RS_type)` / `compute_γ_air_Cabannes!(λ₀, n2, o2)` | **sanghavi** | Sanghavi: computes from molecule-specific Cabannes γ's using `compute_γ_mol_Cabannes!`, returns `(γ_air_Cabannes, ϖ_air_Cabannes)`. Unified `compute_γ_air_Cabannes!(RS_type)` and `compute_γ_air_Cabannes!(n2,o2)` return only γ_air_Cabannes, and use `γ_C_Rayl` as if it were already the Cabannes γ. | **Sanghavi's corrects an identified bug:** sanghavi comment line 411 reads: *"In the old version, `mol.effCoeff.γ_C_Rayl` was assumed to be γ_mol_Cabannes, but it actually was indeed γ_mol_Rayleigh. This assumption has now been corrected, and the correct γ_mol_Cabannes is computed in the function below."* Unified still has the bug. |
| `compute_γ_mol_Cabannes!(λ₀, mol)` (sanghavi) vs `compute_γ_mol_Rayleigh!(λ₀, mol)` (unified) | **sanghavi renames + rewrites** | Sanghavi: treats `mol.effCoeff.γ_C_Rayl` as γ_mol_Rayleigh, and derives γ_mol_Cabannes via `tmpN/tmpD` from `ϖ_Cabannes`. Unified treats `γ_C_Rayl` as γ_mol_Cabannes and derives γ_mol_Rayleigh from it. | **Inverted interpretation of what `γ_C_Rayl` stores.** Sanghavi's interpretation matches its `compute_effective_coefficents!` which defines `γ_C_Rayl = 3γ²/(45α̅² + 4γ²)` (the Rayleigh depolarization). This is a physics correctness issue. |
| `get_n₀_n₁(ieJ₁⁺, Δ)` | — | Unified uses the optimized branch-free form; sanghavi uses `findall` | Both equivalent; **unified faster** per commit `854b44c`. No action. |

### 3.3 `raman_atmo_prop.jl`

| Item | First-on | What differs | Assessment |
|---|---|---|---|
| Formulas in `getRamanSSProp!` for `RRS` and `VS` | equal | `RS_type.ϖ_λ₁λ₀ = reverse(atmo_σ_RRS)/atmo_σ_Rayl` (sanghavi) vs `atmo_σ_RRS[end:-1:1]/atmo_σ_Rayl` (unified) | Semantically identical. |
| `RS_type.i_λ₁λ₀ = reverse(index_raman_grid)` (sanghavi RRS) vs `index_raman_grid[end:-1:1]` (unified) | equal | — | Equivalent. |
| `nm_per_cm` (sanghavi) vs `nm_per_m` (unified) constant naming | — | Both = `1e7`. Unified moved these to `raman_constants.jl`. | No physics difference. |
| Type-stable init `grid_in = AbstractRange{Float64}[]` | unified | Sanghavi: `grid_in = []` (Vector{Any}) | Cosmetic/perf. |
| `i_λ₁λ₀_all = unique(vcat(…))` (sanghavi) vs `unique(cat(…, dims=(1)))` (unified) | — | Equivalent. | — |

No physics differences.

### 3.4 `raman_stellar_prop.jl` and `stellar_inelastic_helper.jl`

Type-stable init, `nm_per_m` vs `nm_per_cm` naming, docstring additions. No physics differences.

### 3.5 `src/Inelastic/src/inelastic_cross_section.jl` — CRITICAL FORMULA DIFFERENCE

| Line | sanghavi | unified | Assessment |
|---|---|---|---|
| Mean polarizability α̅ frequency correction | `α̅ = α̅₀₀*(1 + α_b*T + α_c*T²)/(1-(c*ν_eff/ω₀)²)` | `α̅ = α̅₀₀*(1 + α_b*T + α_c*T²)/(1-(2π*c*ν_eff/ω₀)²)` | **Physics change:** sanghavi removed the `2π` factor in commit `083353b` ("Update lin+Raman", 2026-04-01, sanghavi). Buldakov et al. 1996 Eqs. 36a-39b use `(ω/ω₀)²` where `ω = 2π·c·ν̃`. If `ν_eff` is in angular frequency units, sanghavi is wrong; if `ν_eff` is in wavenumber units (cm⁻¹) and ω₀ is also wavenumber, sanghavi is correct. **Need to confirm unit convention of `ν_eff` (called with `ν̃` in wavenumber in callers — so sanghavi's form is likely correct).** |
| `compute_energy_levels!` loop | swapped `v/J` order + `E₁_pow/E₂_pow` accumulators | standard nested loop with `^` | Optimization; numerically identical. |

### 3.6 `src/Inelastic/src/apply_lineshape.jl`

Physical constants moved to `raman_constants.jl` on unified. No formula differences.

### 3.7 `inelastic_helper_old.jl` (sanghavi only)

Legacy file kept by sanghavi; not imported. Ignore.

---

## 4. Bug fixes scan (sanghavi commits with "fix"/"correct"/"bug"/"sign" keywords)

| Sanghavi commit | Files touched | Status on unified |
|---|---|---|
| `e0ab098` "Fix R⁺⁻ bug in inelastic ScatteringInterface_11" | `interaction_inelastic.jl` | **on unified** (native to unified, sanghavi also has it) |
| `59f8de8` "Port sanghavi branch physics fixes to forward kernels" | `elemental.jl`, `elemental_lin.jl`, `rt_kernel.jl`, `rt_run.jl` | **on unified** (is the port commit) |
| `b8227c2` "Fix refractive index sign convention in Mie scattering" | Scattering | **on unified** |
| `e00e49e` "Regenerate PCW reference with corrected sign convention" | Scattering test references | **on unified** |
| `a5e0de5` "Fix mod→mod1 in doubling kernels" | `doubling.jl` | **on unified** |
| `7376422` "Fix bugs for Rayleigh = Cabannes+RRS match" | `doubling_inelastic.jl`, `elemental.jl`, `elemental_inelastic.jl`, `interaction_inelastic.jl`, `rt_kernel.jl`, `compEffectiveLayerProperties.jl`, `postprocessing_vza.jl`, `rt_run.jl`, `inelastic_helper.jl`, `compute_Z_matrices.jl` | **partially on unified.** Kernel-level changes: unified. Inelastic_helper `compute_ϖ_Cabannes` / `compute_γ_air_Rayleigh` / `compute_γ_air_Cabannes` / `compute_γ_mol_Cabannes` changes: **MISSING on unified** (see §3.2). |
| `083353b` "Update lin+Raman" | `inelastic_cross_section.jl` (drop `2π` in α̅ formula), etc. | **MISSING on unified** (α̅ formula still has `2π`). |
| `e745c35` "correct δBGE truncation" | `Scattering/truncate_phase.jl` | Check Inventory A/E; likely Scattering-side. |
| `d6afa16`, `a04d2a2`, `295cf47` "Add single scatter" | `rt_kernel_ss.jl`, `interaction_ss.jl`, `rt_run_ss` in `rt_run.jl`, `rt_run_lin.jl` | `rt_kernel_ss.jl` + `interaction_ss.jl` **on unified** (files included). `rt_run_ss` entry-point driver **MISSING on unified** — see §5. |
| `c472d20` "Adjust scattering test for epsilon sign convention" | Tests | not audited (test scope). |

---

## 5. Single-scatter approximation kernel

### What it does

`rt_run_ss` computes **only the single-scattering component** of the signal by using the elemental layer directly without doubling (`ndoubl_ss = 0`), then calls `interaction_ss!` / `interaction_inelastic_ss!` to accumulate source terms top-down. It supports `noRS`, `RRS`, and `VS_0to1`/`VS_1to0` types and is intended as a reference/debug path and for fast approximate computations.

### Files

**Kernel code — present on BOTH:**
- `src/CoreRT/CoreKernel/rt_kernel_ss.jl` — `rt_kernel_ss!(::noRS, …)`, `rt_kernel_ss!(::Union{RRS,VS_0to1,VS_1to0}, …)` (legacy + hybrid dispatches)
- `src/CoreRT/CoreKernel/interaction_ss.jl` — `interaction_ss!`, `interaction_inelastic_ss!(::RRS, …)`, `interaction_inelastic_ss!(::Union{VS_0to1_plus,VS_1to0_plus}, …)`

**Driver function — only on sanghavi (MISSING on unified):**
- `rt_run_ss(model::vSmartMOM_Model; i_band=1)` and `rt_run_ss(RS_type::AbstractRamanType, model, iBand)` in `src/CoreRT/rt_run.jl` (sanghavi lines 258–266 and 450–600+). Also a `rt_run_test_ss` helper.
- A linearized variant `rt_run_ss` exists in sanghavi `src/CoreRT/rt_run_lin.jl`. (Out of scope per linearized-Raman policy.)

**Export leak on unified:** `src/CoreRT/CoreRT.jl` line 129 exports `rt_run_ss`, but the symbol is **not defined** anywhere in unified's source tree. A `using vSmartMOM; rt_run_ss(…)` call will throw `UndefVarError`.

### Dependencies

`rt_run_ss` depends on:
1. Existing `rt_kernel_ss!` dispatches — **on unified.**
2. `interaction_ss!` + `interaction_inelastic_ss!` — **on unified.**
3. `make_added_layer`, `make_composite_layer`, `expandOpticalProperties`, `constructCoreOpticalProperties`, `extractEffectiveProps` — all present on unified (though with different signatures; verify at port time).
4. `postprocessing_vza!` with a single-scatter-aware call path — sanghavi's `rt_run_ss` calls `postprocessing_vza!` the same way as `rt_run`. No ss-specific postprocessing code is needed.
5. References `model.ϖ_Cabannes[iBand[1]]` and `RS_type.ϖ_λ₁λ₀ .*= (1 - model.ϖ_Cabannes[iBand[1]])/sum(RS_type.ϖ_λ₁λ₀)` — present on both.

### Port effort

Low-to-medium. Copy the `rt_run_ss` body from sanghavi `src/CoreRT/rt_run.jl` into unified's `src/CoreRT/rt_run.jl`, adapt `vSmartMOM_Model` → `RTModel` field-access (`model.obs_geom` → `model.geometry`, `model.quad_points` → `model.quad_points`, `model.params` → accessors, `model.τ_abs` → `model.optics.τ_abs`, `model.ϖ_Cabannes` → `model.optics.rayleigh.ϖ_Cabannes` or similar). Then verify `iBand` handling and `rt_kernel_ss!` dispatch signatures match unified. Out-of-scope for linearized variant.

---

## 6. Port list for sanghavi-unified

**Per authority rule (2026-04-19):** All items touching the inelastic path (`src/Inelastic/`, `*_inelastic.jl`, inelastic-specific portions of `atmo_prof.jl`, `model_from_parameters.jl`, `types.jl`, and the `RS_type <: noRS` branches of `rt_run.jl`) **land via the wholesale Inelastic port in Phase 1b** of `plans/IMPLEMENTATION_PLAN_v2.md`. No per-item P0/P1/P2 evaluation is required. The priority-ordered framing below is retained for historical reference.

**Items 1 and 2 below** (Cabannes/depolarization physics fixes; α̅ frequency correction) land via the wholesale port per the authority rule. Item 2's unit-convention question is **closed — sanghavi is correct** per amendments §2.4 (user verification).

**Item 3** (`rt_run_ss` driver port) is **not inelastic per se** but an export-leak fix; it lands in **Phase 1c** of the v2 plan.

**Item 4** (`apply_D_matrix_elemental!` scalar shortcut) is **dropped** per amendments §2.2 / v2 Phase 1b. Exception from the authority rule: the dispatch-correctness risk in a merge is not worth the small performance gain; revisit as a one-line follow-up PR if scalar-mode performance becomes a concern.

**Items 5 and 6** (organizational decisions about `rt_kernel!` dispatch structure and `noRS.ϖ_Cabannes` default length) remain organizational; resolve at port time in Phase 1b.

Original priority framing retained below for historical reference.

---

**P0 — correctness bugs affecting Raman numerics, must port:**

1. **`inelastic_helper.jl` Cabannes/depolarization physics fixes** (sanghavi commit `7376422`). Rewrite `compute_ϖ_Cabannes(RRS, λ₀)`, `compute_ϖ_Cabannes(VS_plus, λ₀)`, `compute_γ_air_Cabannes!`, `compute_γ_air_Rayleigh!`, and rename/rewrite `compute_γ_mol_Rayleigh!` → `compute_γ_mol_Cabannes!` (inverting the γ_C_Rayl interpretation). See §3.2 for exact formula diffs. **Without this port, `γ_C_Rayl` is used as if it were γ_mol_Cabannes but the stored value is γ_mol_Rayleigh, corrupting the Cabannes phase-matrix depolarization for RRS and VS.**

2. **`inelastic_cross_section.jl` α̅ frequency correction** (sanghavi commit `083353b`). Change `1-(2π*c*ν_eff/ω₀)^2` → `1-(c*ν_eff/ω₀)^2`. ~~**Verify unit convention first**~~ — **Closed per amendments §2.4: sanghavi is correct.** `ν_eff` is wavenumber (cm⁻¹), `ω₀` is stored in wavenumber units; no `2π` factor. Verified against Buldakov et al. 1996 Eqs. 36a–39b by the user. Port with the one-line verification comment specified in amendments §2.4.

**P1 — missing feature / export leak:**

3. **Port `rt_run_ss` driver function** from sanghavi `src/CoreRT/rt_run.jl` to unified. Adapt field accesses to unified's `RTModel` structure. Also port the `rt_run_test_ss` wrapper. Currently `rt_run_ss` is exported from unified but undefined → import-time silent, runtime `UndefVarError`. Lands in Phase 1c of v2 plan.

**P2 — correctness nice-to-have:**

4. ~~**Port `apply_D_matrix_elemental!` scalar early-return**~~ — **Dropped per amendments §2.2.** The dispatch-correctness risk in a merge is not worth the small performance gain. One-line follow-up post-merge if scalar-mode performance becomes a concern.

**P3 — organizational (pick one style; no numerical effect):**

5. Decide: keep unified's split `rt_kernel!` dispatches (separate RRS/VS and RRS_plus/VS_plus methods) vs. sanghavi's merged union. Unified's is cleaner; keep unified.

6. Decide: unified has `noRS.ϖ_Cabannes = [1.0, 1.0, 1.0]` (3-element default) vs sanghavi `[1.0]` (1-element). Unified's is more tolerant; keep unified — but ensure `noRS_plus` and the multi-band expansion still work with single-band configs.

**Out of scope per policy:**

- Linearized Raman (`rt_run_lin` path for RRS/VS) — sanghavi has a `rt_run_ss` linearized variant; do not port.
- Inventory A items: `InteractionWorkspace` staging, `get_n₀_n₁` optimization, FP32 path.
- Inventory C items: `SIF₀` plumbing through `solar_spec` and elemental_inelastic_plus.

---

## 7. Open questions

**Resolution status as of 2026-04-19 user review session:** Q1, Q2, Q5 are closed by `plans/PLAN_AMENDMENTS_2026-04-19.md`. Q3, Q4, Q6 remain open and resolve in-flight during the Phase 1b port.

1. ~~**Unit convention for `ν_eff` / `ω₀` in `compute_effective_coefficents!`.**~~ **Resolved (amendments §2.4):** `ν_eff` is wavenumber (cm⁻¹); `ω₀` is stored in wavenumber units; no `2π` factor. Verified correct against Buldakov et al. 1996 Eqs. 36a–39b by the user. Port sanghavi's form with the verification comment from amendments §2.4.

2. ~~**Which storage convention does `σ_Rayl_coeff` follow?**~~ **Resolved by authority rule:** sanghavi's `inelastic_helper.jl` is the code that lands; the convention is whatever sanghavi's code assumes. Separately, amendments §4 Phase 0 adds a grep over `src/Inelastic/src/` to confirm the convention for future-maintainer documentation — that's a documentation task, not a decision gate.

3. **Should `compute_ϖ_Cabannes` take `depol` argument path be live or dead?** Sanghavi commented it out; unified kept it live but returns a simple depolarization-derived value. Who calls the 3-argument variant? A grep on each branch would clarify the call-site contract. Resolve in-flight during Phase 1b port.

4. **`rt_run_ss` for inelastic_plus_plus (RRS_plus/VS_plus)?** Sanghavi's `rt_run_ss` uses `RS_type<:Union{RRS,RRS_plus}` to dispatch ϖ_λ₁λ₀ normalization but the `rt_kernel_ss!` has no explicit RRS_plus/VS_plus dispatch (only RRS, VS_0to1, VS_1to0). Is single-scatter needed for `_plus` variants or can it be noRS/RRS-only? Resolve in Phase 1c when porting the driver.

5. ~~**Does unified's `interaction_hdrf!` / `elemental_canopy!` / canopy-surface code need to interact with Raman at all?**~~ **Resolved (amendments §2.6):** Canopy + Raman coupling is a future project, not in scope for this merge. Unified's canopy code is noRS-only; sanghavi has no canopy. The combination `CanopySurface + RRS` is untested and has undefined behavior. No code guard, no smoke test; land the TODO comment specified in amendments §2.6 at the canopy dispatch site during Phase 1b.

6. **Noted missing: `compute_γ_air_Cabannes!` / `compute_γ_air_Rayleigh!` callers on each branch.** Confirm they are not on some cold path that hides the bug before declaring P0 priority. A grep sweep over CoreRT/atmo_prof/compEffectiveLayerProperties will find the call sites and tell us whether these enter forward RT for RRS/VS or only SIF/diagnostic paths. Resolve in-flight during Phase 1b port.
