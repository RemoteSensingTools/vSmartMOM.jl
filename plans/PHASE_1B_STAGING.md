# Phase 1b staging plan — Inelastic port (sanghavi authority)

**Authored 2026-04-19.** Phase 1a committed as `8eee415`. This document stages Phase 1b so the next session has a concrete file-by-file path forward.

Authority rule (from `PLAN_AMENDMENTS_2026-04-19.md` §1): **sanghavi is the authority for the inelastic path.** Per-file decisions are already settled; this doc enumerates the concrete work, the compatibility surfaces, and the commit boundaries.

---

## 1. What lives where (inventory)

Worktree paths:
- unified (sanghavi-unified @ `8eee415`): `/home/sanghavi/code/github/uni_vSmartMOM/`
- sanghavi (@ `9ee9a75`): `/home/sanghavi/code/github/vSmartMOM.jl/`

### 1.1. Files that fully replace (sanghavi → unified, wholesale copy)

| File | uni lines | san lines | Notes |
|---|---|---|---|
| `src/Inelastic/InelasticScattering.jl` | 51 | 52 | Module entry; confirm `include` paths match unified's `src/Inelastic/src/` layout. |
| `src/Inelastic/inelastic_helper.jl` | 911 | 914 | Shared helpers. |
| `src/Inelastic/raman_atmo_prop.jl` | 480 | 581 | Atmospheric Raman properties. |
| `src/Inelastic/raman_stellar_prop.jl` | 483 | 495 | Stellar Raman properties. |
| `src/Inelastic/stellar_inelastic_helper.jl` | 755 | 714 | Stellar Raman helpers. |
| `src/Inelastic/stellar_types.jl` | 245 | 367 | Stellar Raman types. |
| `src/Inelastic/types.jl` | 581 | 521 | **Structural API change**, see §2.1. |
| `src/Inelastic/src/apply_lineshape.jl` | 61 | 71 | Lineshape application. |
| `src/Inelastic/src/inelastic_cross_section.jl` | 349 | 297 | **Add `α̅` verification comment**, see §3.2. |
| `src/Inelastic/src/molecular_constructors.jl` | 212 | 212 | 3-line delta. |
| `src/Inelastic/src/raman_constants.jl` | 153 | 153 | Constants. |
| `src/CoreRT/CoreKernel/doubling_inelastic.jl` | *(uni has)* | | Kernel doubling for RRS/VS. |
| `src/CoreRT/CoreKernel/elemental_inelastic.jl` | *(uni has)* | | Elemental kernel. |
| `src/CoreRT/CoreKernel/elemental_inelastic_plus.jl` | *(uni has)* | | Elemental kernel (plus variant). |
| `src/CoreRT/CoreKernel/interaction_inelastic.jl` | *(uni has)* | | Interaction kernel. |
| `src/CoreRT/CoreKernel/raman_kernel_test.jl` | *(uni has)* | | Test kernel (verify still relevant; drop if orphan). |

### 1.2. Files not ported from sanghavi

- `src/Inelastic/inelastic_helper_old.jl` — 924-line frozen reference on sanghavi. **Do not port.** (Already locked in carve-outs.)
- `src/Inelastic/src/plots/` — dev detritus. **Do not port.**
- `src/Inelastic/src/prototype.jl` — dev detritus. **Do not port.**

### 1.3. Files that merge (not copy)

These have substantial unified-side scaffolding that must be preserved. Each one needs a targeted diff, not a wholesale replace.

| File | uni lines | san lines | Risk | Merge focus |
|---|---|---|---|---|
| `src/CoreRT/rt_run.jl` | 288 | 625 | **High** | Keep unified's `RTModel` scaffolding; graft sanghavi's `RS_type <: noRS` branches; land `_interaction_ws = nothing` placeholder (workspace is Phase 4). See §3.1. |
| `src/CoreRT/tools/model_from_parameters.jl` | ~1200 | 1225 (flat) | **High** | Adopt sanghavi's inelastic-specific `model_from_parameters(RS_type, …)` methods (VS variants, `γ_air_Cabannes`, `compute_γ_air_Rayleigh!`). The `noRS` path stays unified. |
| `src/CoreRT/tools/atmo_prof.jl` | 527 | 457 | Medium | Only inelastic-specific helpers to pull over (e.g. `construct_atm_layer` RRS branch if sanghavi has extra). Rayleigh + reduce_profile already Phase 1a. |
| `src/CoreRT/types.jl` | ~1200 | ~1100 | Low–Medium | ieR/ieT field layouts match (both branches). Confirm no sanghavi-only fields missing (e.g. `fscattRayl` flags). |

---

## 2. Compatibility surfaces (cross-file mechanical work)

### 2.1. `AbstractRamanType{FT}` → `AbstractRamanType` (strip type parameter)

**Unified declares** `abstract type AbstractRamanType{FT} end` and subtypes with `<: AbstractRamanType{FT}`.
**Sanghavi declares** `abstract type AbstractRamanType end` with unparameterized subtypes.

Per authority rule: sanghavi wins. Stripping the `{FT}` touches **44 occurrences across 7 files**:

| File | Count | Touchpoint kind |
|---|---|---|
| `src/Inelastic/types.jl` | 30 | Replaced wholesale (§1.1). |
| `src/Inelastic/stellar_types.jl` | 7 | Replaced wholesale. |
| `src/Inelastic/inelastic_helper.jl` | 1 | Replaced wholesale. |
| `src/CoreRT/rt_run.jl` | 2 | Manual merge (§3.1). |
| `src/CoreRT/rt_run_lin.jl` | 2 | Strip `{FT}` from constraint; lin-Raman branches stay excluded per scope. |
| `src/CoreRT/rt_run_multisensor.jl` | 1 | Strip `{FT}`; keep unified logic. |
| `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl` | 1 | Strip `{FT}`. |

**How to apply:** Do the strip as part of the commit that replaces `src/Inelastic/types.jl`, so no intermediate state compiles with mixed declarations. Mechanical: `sed -i 's/AbstractRamanType{FT}/AbstractRamanType/g'` for the non-Inelastic callsites, but eyeball each one — a few (e.g. `rt_run_lin.jl`) also have function signatures that need review for lin-Raman scope.

---

## 3. File-specific notes

### 3.1. `src/CoreRT/rt_run.jl` merge

Unified's 288 lines are the `RTModel`-aware forward driver (Phase 0 baseline). Sanghavi's 625 lines add:
- `rt_run_bck` (older driver body, reuses v1 `vSmartMOM_Model`). Not ported.
- `rt_run_ss` single-scatter driver. Ported in **Phase 1c**, not 1b.
- `InteractionWorkspace` allocation site around line 122: `_interaction_ws = (typeof(RS_type) <: noRS) ? nothing : InteractionWorkspace(composite_layer, added_layer; staged=true)`.
- `RS_type.bandSpecLim`, `RS_type.F₀`, `RS_type.Z⁺⁺_λ₁λ₀`/`Z⁻⁺_λ₁λ₀` set-up, `InelasticScattering.compute_Z_moments`, `RS_type.fscattRayl` updates, `interaction_inelastic!` branches.

**Phase 1b scope for rt_run.jl:**
1. Take unified's `rt_run(RS_type::AbstractRamanType, model, iBand)` body (already present, `noRS()` path works).
2. Graft sanghavi's `RS_type`-branch pre-setup (bandSpecLim, F₀, Z⁺⁺_λ₁λ₀, fscattRayl).
3. Leave `_interaction_ws = nothing` at the allocation site (per amendments §4 Phase 1b; workspace landing is Phase 4).
4. Pass `workspace=_interaction_ws` kwarg through `rt_kernel!` / `interaction!` / surface — kernels already accept `workspace === nothing` (verify).
5. Strip the `{FT}` type parameter from the two `AbstractRamanType{FT}` occurrences.

**Do NOT** port `rt_run_bck` or `rt_run_ss` in Phase 1b. `rt_run_ss` is Phase 1c.

### 3.2. `src/Inelastic/src/inelastic_cross_section.jl` α̅ verification comment

Per amendments §2.4, land this comment at the `α̅` frequency-correction site (the `(1 - (c·ν_eff/ω₀)²)` expression):

```julia
# Frequency correction (1 - (c·ν_eff/ω₀)²). ν_eff is wavenumber (cm⁻¹),
# ω₀ is stored in wavenumber units; no 2π factor. Verified correct against
# Buldakov et al. 1996 Eqs. 36a-39b by [AUTHOR/SOURCE].
```

**Blocking input required from user** before commit: the `[AUTHOR/SOURCE]` citation (notebook, paper section, commit message, email thread). The plan explicitly flags this as an open item to fill during implementation.

### 3.3. Canopy + Raman TODO comment

Per amendments §2.6, land this at the canopy dispatch site (likely `src/CoreRT/Surfaces/canopy_surface.jl` — candidate anchor; search for `create_surface_layer!(::CanopySurface, …)`):

```julia
# TODO: Canopy + Raman (RRS/VS) coupling is not currently implemented. The
# canopy BRDF code is noRS-only; calling rt_run with CanopySurface + RRS has
# undefined behavior. Future work: couple canopy to the inelastic path.
```

### 3.4. σ_Rayl_coeff storage convention

Per amendments §Phase 0, grep `src/Inelastic/src/` for `σ_Rayl_coeff` (or synonyms) and document in a header comment the unit convention assumed. This is a documentation step, not a code change.

---

## 4. Proposed commit boundaries for Phase 1b

Each numbered bullet is one logical commit.

1. **"Port Inelastic types + strip AbstractRamanType{FT}"** — replaces `src/Inelastic/types.jl` and `stellar_types.jl` wholesale; strips `{FT}` from the 5 non-Inelastic callsites. Package should compile after this.
2. **"Bulk-port Inelastic helpers (sanghavi authority)"** — replaces `InelasticScattering.jl`, `inelastic_helper.jl`, `raman_atmo_prop.jl`, `raman_stellar_prop.jl`, `stellar_inelastic_helper.jl`, and all `src/Inelastic/src/*.jl`. Include the α̅ verification comment.
3. **"Port inelastic CoreKernel files"** — replaces the 5 `*_inelastic*.jl` files under `src/CoreRT/CoreKernel/`.
4. **"Merge rt_run.jl: sanghavi inelastic branches on RTModel scaffolding"** — the rt_run.jl merge per §3.1. Leave workspace placeholder.
5. **"Port inelastic-specific model_from_parameters methods"** — add sanghavi's VS/RRS-specific methods; keep unified's noRS method unchanged.
6. **"Land canopy+Raman TODO; σ_Rayl_coeff convention note"** — tiny doc-only commit.
7. **"Wire test_forward_raman.jl into runtests.jl; commit phase1b_inelastic_port.jld2"** — closes the Phase 1b gate.

Each commit should leave the package compiling and `test_forward_noRS` passing. The full regression gate (`test_forward_raman.jl` vs sanghavi reference) lands at step 7.

---

## 5. Phase 1b exit criterion

From `IMPLEMENTATION_PLAN_v2.md` §Phase 1b:

> `test_forward_raman.jl` (currently orphan on unified — wire in during this sub-phase) matches sanghavi reference within user tolerance on a toy 1-band RRS config. Reference JLD2 committed to `test/reference/phase1b_inelastic_port.jld2`.

Concrete checks before flipping Phase 1b complete:
- [ ] `test_forward_noRS.jl` still passes (no elastic-path regression).
- [ ] `test_forward_raman.jl` passes on CPU against `phase1b_inelastic_port.jld2` reference; user confirms tolerance per-quantity at execution time.
- [ ] CUDA test re-runs: if Raman GPU test was omitted during Phase 1a per `feedback_phase1a_test_subset.md`, re-enable and verify it matches the same JLD2 reference.
- [ ] No `InteractionWorkspace` fields populated yet — `_interaction_ws = nothing` is the expected allocation state.

---

## 6. Open items for user on Phase 1b kickoff

1. **α̅ verification citation** (§3.2) — needed to fill `[AUTHOR/SOURCE]` placeholder before commit 2 lands.
2. **Reference JLD2 generation** — need user to approve the toy 1-band RRS config (spectral range, layers) that generates `phase1b_inelastic_port.jld2`. Target <60 s warm run per amendments §4 Phase 2, but for Phase 1b it's just a regression artifact, so smaller is fine.
3. **Per-quantity tolerance on I/Q/U/V/ieR/ieT** for the phase1b regression check. Per the plan, user supplies at execution time — this is the moment.

---

## 7. Session handoff

- Phase 1a: committed `8eee415` on `sanghavi-unified`, `runtests.jl` green (442 pass + 1 expected skip).
- Phase 1b: this document is the next-session launchpad.
- Tooling: `julia-mcp` available via `.mcp.json`; persistent session at `env_path=/home/sanghavi/code/github/uni_vSmartMOM/test/` (TestEnv auto-activates).
- Memory updated: `project_sanghavi_unified_phase_status.md` now points at Phase 1b.
