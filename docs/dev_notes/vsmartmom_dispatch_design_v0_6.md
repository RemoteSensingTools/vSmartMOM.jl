# vSmartMOM.jl — architecture design and Fourier work

**Document version:** v0.6
**Status:** working draft for Sanghavi review.
**Branch target:** `sanghavi-unified` (with `unified-vsmartmom` as the architectural merge base, physics layered in from `sanghavi`).

**Document set (delivered together; cross-references resolve within this set):**
- `vsmartmom_dispatch_design_v0_6.md` — *this document* (architecture).
- `standalone_ss_solver_plan.md` — kernel-based exact-SS solver: pol-type-generic, dual-form optical properties aligned with v0.5/v0.6, AD-vs-handcoded design experiment, FO-equivalent diagnostic for SS-paths-to-sensor validation.
- `vlidort_baseline_suite.md` — VLIDORT golden-standard validation using the shipped Siewert (2000) and solar_tester fixtures.
- `source_terms_architecture_v0_6.md` — source terms as first-class affine RHS contributors: solar SFI, thermal emission, SIF, lidar/lunar sources, and the AD/hand-written differentiation boundary.
- `dev_notes/exact_ss_reference/` — standalone four-path reference (Julia + Python validation companion + derivations).

**Pre-existing committed artifacts referenced:**
- `docs/dev_notes/LINEARIZATION_BUGS.md` — linearization bug catalogue and redesign track (committed in `vSmartMOM.jl`).
- `src/CoreRT/types.jl` — current type definitions; especially `CoreScatteringOpticalProperties{FT,FT2,FT3}` at line 1052 with `τ::FT, ϖ::FT2, Z⁺⁺::FT3, Z⁻⁺::FT3`.
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl` — current linearization kernel; the "Bug 19" note at line 32 about the Z chain rule needing to be applied *before* doubling-adding mixes indices is the precedent for our chain-rule design in §6 of `standalone_ss_solver_plan.md`.

---

## Preamble — for Sanghavi

This is the sixth iteration. The architectural skeleton has held; v0.6 sharpens the SS-correction story substantially through three major shifts:

1. **Solver moves into the model.** `model.solver::AbstractRTSolver` is now part of the model's complete description. `run_rt` takes only the model, the per-call mode (forward vs lin), and the output spec.

2. **Precheck-driven defaults** for which corrections are applied per band, with thresholds derived empirically from the required-Nquad framework (`standalone_ss_solver_plan.md` §5) rather than picked by hand.

3. **`ExactFirstOrderOnly` as a first-class RT solver.** SS-only computation (skipping the MS solver entirely) becomes a peer of `FullMOM` under `AbstractRTSolver`.

Plus:
- **`ExactSFIPhase` as the full-MOM correction architecture, retained from v0.5.** The correction `(I - M_trunc)^-1 (J_exact - J_trunc)` propagates the corrected first-scatter source through the MS operator. This corrects strictly more paths than VLIDORT FO and is the architectural reason the standalone solver's post-hoc back-correction (`standalone_ss_solver_plan.md` §7) is treated as a *diagnostic against VLIDORT FO* — equivalent to FO scope — and not as a substitute for `ExactSFIPhase`.

### What v0.6 changes from v0.5

1. **Solver in the model** (§1.3, §7). `AbstractRTSolver` is a model field, not a `run_rt` kwarg. `run_rt(model; mode, outputs)` is a two-kwarg entry point.

2. **Precheck for default correction selection** (§3.6). `precheck_truncation(model, iBand)` walks contributors and surface, returns which corrections are needed. Thresholds derived from the required-Nquad framework in `standalone_ss_solver_plan.md` §5, not picked by hand.

3. **`AbstractRTSolver` as user-facing dispatch axis** (§7.3): `FullMOM` (default), `ExactFirstOrderOnly` (SS-only).

4. **Standalone exact-SS solver as separate Piece A** (`standalone_ss_solver_plan.md`): pol-type-generic from day one, dual-form optical properties aligned with this document's §5.5 / v0.5 §6.

5. **VLIDORT baseline suite as separate Piece B** (`vlidort_baseline_suite.md`): leverages shipped fixtures (Siewert 2000, solar_tester multi-task), no PyVLIDORT install needed initially.

6. **Possible bug found in existing Cox-Munk correction** (§3.4): missing `exp(-τ_total/μ_v)` factor. **See Q3.**

### What v0.6 retains from v0.5 unchanged

- **`AbstractSFIPhase` as the full-MOM correction mechanism** (§6) — the first solar scattering event source `J0±` can be built from either truncated or untruncated phase. `ExactSFIPhase` builds it from untruncated `(τ, ϖ, Z)` (all three; not a mix); `TruncatedSFIPhase` builds it from primed `(τ', ϖ', Z')`. The dual-form `TruncatedAndExactScatteringOpticalProperties` carries both. The MS operator stays truncated; the *source* changes.
- **Mathematical framing** (§6.2): `ΔI = (I - M_trunc)^-1 (J_exact - J_trunc)` propagating the corrected source through the MS operator. Strictly more paths than VLIDORT FO; not bit-equivalent to FO.
- **Energy-budget normalization** (§6.3): `J_exact` requires *all three* unprimed; mixing forms double-counts.

### Six theory/philosophy questions for you

**Q1 — Solver in the model.** v0.6 puts `AbstractRTSolver` into the model. Clean enough?

**Q2 — Precheck-driven correction defaults.** Framework right; thresholds derived from required-Nquad framework (Piece A §5)? Defensible?

**Q3 — Possible bug in existing Cox-Munk SS correction.** Missing upward attenuation `exp(-τ_total/μ_v)`. Standalone reference includes both attenuations with BRDF reciprocity passing to ~10⁻⁶. Was the omission intentional? *Detail: §3.4.*

**Q4 — Mixed truncation: TMS for FO-equivalent paths 1+2 (when used standalone), untruncated for paths 3+4.** Defensible? Note: this is the *standalone solver* question (Piece A); for full-MOM `ExactSFIPhase` everything is unprimed throughout. *Detail: §6.3, also Piece A §3.*

**Q5 — Canopy hotspot is silently missing today.** Same as v0.5 Q4. v0.6's `CanopyHotspotPath` under the standalone solver's `AbstractSSPath` is the prototype seam; full integration via `ExactSFIPhase` source modification when canopy is the surface.

**Q6 — `NoScattering` layer kind, band-level.** Same as v0.5 Q5. Layer kind classified at construction.

### Three narrower questions

**Q7 — `n_fourier_moments(::CoxMunk) = 2·Nstreams`.** Matches VLIDORT.

**Q8 — `GaussQuadFullSphere` retirement.**

**Q9 — Threshold values for the precheck.** Derived from Piece A §5's required-Nquad framework (empirically calibrated), or hand-set?

### What v0.6 is *not*

- Not a fait accompli. Open questions in §13.
- Not yet implemented; current `sanghavi-unified` has none of these changes.
- Not the linearization redesign — separate track in committed `docs/dev_notes/LINEARIZATION_BUGS.md`. The standalone SS solver's §6 hybrid AD-handcoded experiment (Piece A) is the small-scale prototype that informs the broader redesign.

---

## 1. The dispatch rule (governing principle)

### 1.1 The rule

> **Use multiple dispatch at physical and workflow boundaries; keep `if/else` only for scalar numerical cases.**

Physical and workflow boundaries become method dispatch on small concrete strategy structs. `if/else` is fine for scalar numerical cases (Fourier normalization weights, NaN guards, runtime user predicates).

### 1.2 The strategy struct shape

Strategy objects are small concrete structs held in tuples, dispatched via single-method functions. Concrete tuples (not abstract vectors), single-method functions (not switch tables), naming consistency across families.

### 1.3 The model carries solver state

```julia
struct RTModel{ATM, SUR, GEO, QP, SOL, ARCH}
    atmosphere::ATM
    surfaces::SUR
    geometry::GEO
    quad_points::QP
    solver::SOL              # AbstractRTSolver; carries SFI-phase choice, ss_paths spec, l_trunc
    architecture::ARCH
    # ... bands, aerosol_optics, n_layers, etc.
end
```

`run_rt(model; mode=ForwardMode(), outputs=(DirectionalTOA(),))`. Two kwargs. Everything else is in `model`.

### 1.4 The `RTContext` object

Built per `run_rt` call from `model + mode + outputs`. Contains the resolved `ss_paths` tuple (when running `ExactFirstOrderOnly` or building diagnostic SS outputs) and the resolved `sfi_phase` choice (when running `FullMOM` with `ExactSFIPhase`).

### 1.5 What the rule unlocks

- **Step 1A–C** — small contained pieces (§§2–4).
- **Step 1.5** — optical-properties algebra and layer-kind dispatch (§5).
- **Step 2** — `AbstractSFIPhase` source dispatch for full MOM (§6); standalone SS solver kernels in Piece A.
- **Step 3** — `run_rt` consolidation (§7).

### 1.6 Branches that stay (the "fine" list)

- `m == 0` Fourier weight normalization
- Cache initialization
- `maximum(τ·ϖ) > 2eps(FT)` check at layer construction
- NaN/denorm guards
- The `ss_paths === :auto` resolution branch in context construction

---

## 2. Step 1A — per-component `n_fourier_moments`

[Section 2 unchanged from v0.5 — Fourier moment trait per component, `Nstreams ≠ Nquad`, `max_m::Union{Nothing,Int}`, `GaussQuadFullSphere` retirement, source-type forward-compatibility invariant.]

Summary:
- `n_fourier_moments(component, ctx)` trait dispatches Fourier-moment count per component
- VLIDORT's `NMOMENTS = min(2·NSTREAMS - 1, NGREEK_MOMENTS_INPUT)` (verified in `vlidort_main/regular/vlidort_inputs.f90:4138-4141`) translates to `2·Nstreams` (count, not order)
- `Nstreams` (base) separated from `Nquad` (augmented with VZA/SZA)
- `GaussQuadFullSphere` deprecated in v2.0, removed in v2.1

---

## 3. Step 1B — `AbstractSSPath` and the precheck

### 3.1 The current state

`rt_run.jl:310-315`:

```julia
if brdf isa CoxMunkSurface && SFI
    @timeit "SS Correction" apply_ss_correction!(
        R_SFI, brdf, pol_type, vza, vaz, μ₀,
        Array(τ_sum_all[:,end]), max_m, nSpec)
end
```

Three problems: `isa` gate (extending requires more `if`s); inconsistent across entry points (silent Jacobian bug at glint geometries); hardcoded to one correction.

### 3.2 What VLIDORT FO does and where vSmartMOM should diverge

Read directly from `vlidort_focode/FO_VectorSS_RTCalcs_I.f90`. Per the file's own comment at line 82: *"For the Atmospheric Solar Single-scatter and Surface Direct-Beam solar sources."* VLIDORT FO covers exactly two paths:
- **Path 1 — atmospheric SS (sun→atm→sensor)**
- **Path 2 — surface direct-beam (sun→surface→sensor)**

VLIDORT FO does not cover paths 3 or 4. These are absorbed into the truncated MS solve.

**vSmartMOM's full-MOM architecture is `ExactSFIPhase` (§6), not a post-hoc FO correction.** `ExactSFIPhase` corrects the first-scatter *source* and propagates through the truncated MS operator — strictly more paths than FO covers, but the MS operator stays truncated.

The four-path framework (`AbstractSSPath`, with concrete subtypes `AtmosphericSSPath`, `SurfaceDirectBeamPath`, `AtmosphereToSurfacePath`, `SurfaceToAtmospherePath`, `CanopyHotspotPath`) lives in the *standalone solver* (Piece A) for two purposes:
1. Powering `ExactFirstOrderOnly` (SS-only mode, no MS solver).
2. Providing FO-equivalent diagnostic capability for validating the truncated-MS solve against VLIDORT FO outputs (paths 1+2 only).

### 3.3 Path/correction type hierarchy

Lives in Piece A's standalone module:

```julia
abstract type AbstractSSPath end

struct AtmosphericSSPath <: AbstractSSPath end
struct SurfaceDirectBeamPath{S<:AbstractSurfaceType} <: AbstractSSPath
    surface::S
end
struct AtmosphereToSurfacePath <: AbstractSSPath
    inner_quad::Int
end
struct SurfaceToAtmospherePath{S<:AbstractSurfaceType} <: AbstractSSPath
    surface::S
    inner_quad::Int
end
struct CanopyHotspotPath <: AbstractSSPath end
```

`AbstractSSPath` is the dispatch surface for `apply_correction!`, `exact_ss`, and `truncated_reconstruction`. See `standalone_ss_solver_plan.md` §3 for full pol-type-generic interface.

### 3.4 Sanghavi's existing Cox-Munk correction: prototype and a possible bug

Reading `coxmunk_surface.jl:505-546`, the existing `apply_ss_correction!`:

```julia
for s in 1:nSpec
    atten = μ₀ * exp(-τ_total[s] / μ₀)        # ← downward only
    for si in 1:n
        correction = atten * (M_exact[si, 1] - M_fourier[si, 1])
        R_SFI[iv, si, s] += correction
    end
end
```

⚠️ **Q3 — possible bug**: `atten = μ₀ × exp(-τ_total/μ₀)` is the downward direct-beam attenuation only. There appears to be no upward `exp(-τ_total/μ_v)` factor for surface→TOA propagation.

The standalone reference (`dev_notes/exact_ss_reference/exact_ss_reference.jl`) includes both attenuations:

```julia
function path2_surface_direct(τ_total, μ₀, μ_v, albedo, I₀)
    return (μ₀ * I₀ * albedo / π) * exp(-τ_total/μ₀) * exp(-τ_total/μ_v)
end
```

with BRDF reciprocity tests passing to ~10⁻⁶.

### 3.5 The unified call site

`rt_run.jl:310-315` becomes (after both Step 1B refactor and Step 3 entry-point consolidation):

```julia
@timeit "SS path corrections (diagnostic / ExactFirstOrderOnly)" foreach(
    p -> apply_correction!(R_SFI, p, model.solver, model, ctx),
    ctx.ss_paths)
```

`ctx.ss_paths` is populated at context construction from `model.solver.ss_paths`. `model.solver` dispatches `apply_correction!` between `FullMOM` (when running diagnostically — see §6.5 for the relationship to `ExactSFIPhase`) and `ExactFirstOrderOnly` (direct add).

For `FullMOM` production runs, the primary correction is `ExactSFIPhase` (§6); the standalone-solver path corrections are diagnostic/validation tools.

### 3.6 The precheck

`precheck_truncation(model, iBand)` walks contributors and surface, returns a NamedTuple of which corrections are needed. Thresholds derived from `standalone_ss_solver_plan.md` §5 required-Nquad framework rather than hand-set:

```julia
function precheck_truncation(model, iBand)
    contributors = model.atmosphere.contributors[iBand]
    surface = model.surfaces[iBand]
    Nstreams = model.quad_points.Nstreams
    
    # Required-Nquad framework (Piece A §5) inverts the question
    L_required = maximum(determine_required_l_aerosol(c) for c in contributors
                         if c isa AbstractScatteringContributor; init=2)
    Nstreams_required = ceil(Int, (L_required + 1) / 2)
    needs_atm = Nstreams < Nstreams_required
    
    Nbrdf_required = determine_required_nbrdf(surface)
    needs_surf = needs_specular_correction_via_required_n(surface, Nstreams, Nbrdf_required)
    
    surface_reflects = !is_black(surface)
    needs_atm_to_surf = needs_atm && surface_reflects
    needs_surf_to_atm = needs_atm_to_surf
    needs_canopy_hs = surface isa CanopySurface && !(surface.hotspot isa NoHotSpot)

    return (; needs_atm, needs_surf, needs_atm_to_surf, needs_surf_to_atm, needs_canopy_hs,
              L_required, Nstreams_required, Nbrdf_required)  # diagnostic
end
```

The diagnostic record (last three fields) tells the user *why* the precheck flagged each correction.

### 3.7 SS-paths resolution at context construction

`model.solver` carries the SS-paths spec:

```julia
struct FullMOM{S, P} <: AbstractRTSolver
    sfi_phase::P             # AbstractSFIPhase: ExactSFIPhase() or TruncatedSFIPhase()
    ss_paths::S              # :auto, :none, :all, or NTuple{N, AbstractSSPath} — DIAGNOSTIC ONLY
    l_trunc::Int
end
```

For `FullMOM`, `ss_paths` is *diagnostic-only*: it produces post-hoc FO-equivalent corrections that match VLIDORT FO's scope (paths 1+2 to sensor, not propagated through MS). The primary correction in `FullMOM` is `sfi_phase=ExactSFIPhase()`. See §6 and Piece A §7 for the relationship.

For `ExactFirstOrderOnly`, `ss_paths` *is* the computation — there's no MS solver to propagate through, so post-hoc summing of paths is exactly the right thing.

```julia
struct ExactFirstOrderOnly{S} <: AbstractRTSolver
    ss_paths::S              # :all is the sensible default
end
```

### 3.8 What the precheck delivers

- **Defaults are right.** Sanghavi's typical case (Cox-Munk + Rayleigh + modest aerosol) returns `needs_surf=true, needs_atm=false`. `ExactSFIPhase` becomes the primary `FullMOM` correction; the diagnostic ss_paths defaults to `(SurfaceDirectBeamPath(coxmunk),)` for FO-comparison.
- **Per-band selection automatic.** Different bands have different aerosols; precheck runs per band.
- **Override via `solver.ss_paths` and `solver.sfi_phase` at construction.**
- **Logging hook.** `precheck_truncation(model, iBand)` returns NamedTuple including `L_required, Nstreams_required, Nbrdf_required` — exposed via structured output field for transparency.

---

## 4. Step 1C — `AbstractHotSpot` field on `CanopySurface`

[Section 4 unchanged from v0.5 — `hotspot::HS` parameter on `CanopySurface`, plumbed through `model_from_parameters` and YAML, default `NoHotSpot{FT}()`.]

For `FullMOM + ExactSFIPhase`: when surface is canopy with non-NoHotSpot, the SFI source for the *direct-sun→direct-view* component uses joint gap probability (Kuusk hotspot) instead of independent Beer-Lambert. Diffuse contributions untouched.

For `ExactFirstOrderOnly`: `CanopyHotspotPath` from §3.3 is the direct dispatch.

Q5: are there current canopy benchmarks tuned to no-hotspot output? If so, enabling will change those numbers.

---

## 5. Step 1.5 — optical-properties algebra and layer-kind dispatch

[Section 5 mostly unchanged from v0.5; §5.5 sharpened to align with §6 dual-form requirement.]

### 5.5 Dual-form optical properties — for `ExactSFIPhase` and for the standalone solver

The existing `CoreScatteringOpticalProperties{FT,FT2,FT3}` at `types.jl:1052` carries `τ::FT, ϖ::FT2, Z⁺⁺::FT3, Z⁻⁺::FT3` (forward and backward Z matrices). Already pol-type-generic — `Z⁺⁺` and `Z⁻⁺` are matrices whose dimensions reflect the polarization mode.

For `ExactSFIPhase`, both forms are required as full triplets:

```julia
struct TruncatedAndExactScatteringOpticalProperties{T, E, FT}
    truncated::T            # CoreScatteringOpticalProperties — primed (τ', ϖ', Z'); for MS solver
    exact::E                # CoreScatteringOpticalProperties — unprimed (τ, ϖ, Z); for J_exact
    truncfac::Vector{FT}    # f per layer; metadata for derivation, not a separate input
end
```

**Critical**: the truncated and exact forms are *full* `CoreScatteringOpticalProperties` triplets. This aligns with the v0.5 §6.3 requirement: `J_trunc` uses *all primed* `(τ', ϖ', Z')`, `J_exact` uses *all unprimed* `(τ, ϖ, Z)`. Mixing forms (e.g., primed τ with unprimed ϖ) double-counts the scattering integral and breaks the energy budget.

The MS solver consumes `truncated` directly (current behavior).

`ExactSFIPhase` (§6) reads `exact` to build `J_exact`; reads `truncated` to build `J_trunc`; computes `(J_exact - J_trunc)`; propagates through `(I - M_trunc)^-1`.

The standalone solver (Piece A) reads `exact` for paths 3+4 and uses TMS reconstruction `(τ', ϖ, Z)` with `truncfac` for FO-equivalent paths 1+2 (when running diagnostically against VLIDORT FO). Detail in `standalone_ss_solver_plan.md` §3.

§5.2 (documentation + bug fix at `types.jl:1137`), §5.3 (dimensioned zero), §5.4 (contributor pattern with `exact_phase_function` trait), §5.6 (storage form as type parameter), §5.7 (`NoScattering` layer kind, band-level) unchanged from v0.5.

---

## 6. Step 2 — `AbstractSFIPhase` source dispatch (for `FullMOM`)

[Section retained from v0.5 §6, updated to clarify relationship to standalone solver.]

### 6.1 Motivation

vSmartMOM's MS operators `r/t` use truncated/δ-M phase moments — that's correct for energy conservation in the MS solve. But the SFI source vectors `J0±` are also currently built from truncated phase, which means the *first solar scattering event* is approximated even though we have the exact phase function available.

This is the same architectural gap that VLIDORT FO addresses, but the right vSmartMOM-native abstraction is different from VLIDORT's. vSmartMOM's `J0±` are already in `elemental.jl` and propagated by the MOM doubling-adding algebra. Fixing them at construction (rather than adding a separate FO solver) corrects more paths than VLIDORT FO does.

### 6.2 The mathematical framing

The standard truncated SFI solve is:
```
I_trunc = (I - M_trunc)^-1 J_trunc
```

Using exact first-scatter source:
```
I_exact_first = (I - M_trunc)^-1 J_exact
```

The correction (added to the truncated solve):
```
ΔI = (I - M_trunc)^-1 (J_exact - J_trunc)
```

This corrects the contribution of the first solar scattering event for *all paths beginning with that scatter*:
```
sun → atm scatter → sensor                            (FO-equivalent first-order)
sun → atm scatter → surface → sensor                  (FO-equivalent)
sun → atm scatter → atm scatter → sensor              (higher-order; FO doesn't fix)
sun → atm scatter → surface → atm scatter → sensor    (higher-order; FO doesn't fix)
...
```

**This is hybrid:** the *first solar scattering event* is exact (untruncated phase, untruncated optical depth, untruncated ϖ); *later propagation* uses the truncated MS operator. Strictly more paths corrected than VLIDORT FO; not the same as exact full RT.

### 6.3 The energy-budget normalization

When δ-M truncates the phase function:
- Phase becomes `P_trunc(Θ) = (P(Θ) - f δ(Θ))/(1-f)`, normalized to `4π`
- Optical depth becomes `τ' = (1 - fϖ) τ`
- Single-scattering albedo becomes `ϖ' = (1-f) ϖ / (1 - fϖ)`

For **`ExactSFIPhase` to work correctly**, `J_exact` must be built from *all three unprimed quantities* — `τ, ϖ, Z`. This is exactly the dual-form §5.5 type carries. Mixing forms double-counts.

### 6.4 The `AbstractSFIPhase` dispatch

```julia
abstract type AbstractSFIPhase end
struct ExactSFIPhase <: AbstractSFIPhase end
struct TruncatedSFIPhase <: AbstractSFIPhase end       # current behavior
```

Set on `model.solver`; resolved at context construction; dispatched in `elemental.jl` source-vector construction:

```julia
build_J_source!(J0_plus, J0_minus, ::TruncatedSFIPhase, optical_props, ...) =
    # Use optical_props.truncated (current behavior)
    
build_J_source!(J0_plus, J0_minus, ::ExactSFIPhase, optical_props, ...) =
    # Use optical_props.exact (new)
```

### 6.5 Relationship to standalone solver's post-hoc back-correction (Piece A §7)

The standalone solver's back-correction (`standalone_ss_solver_plan.md` §7) adds `(exact_p1+p2 - truncated_p1+p2)` to `R_SFI` *after* `rt_run`. This corrects the **directly-escaped first-scatter only** — equivalent to VLIDORT FO scope. It does *not* correct first-scatter-then-anything paths because those require propagation through `(I - M_trunc)^-1`, which only happens inside the MS solver.

**Therefore the back-correction is a diagnostic for FO-equivalent comparison**, not a substitute for `ExactSFIPhase`. Specifically:

- Back-correction validates against **VLIDORT FO outputs** (Tasks 1→3 difference in solar_tester) — useful and well-scoped.
- `ExactSFIPhase` validates against **VLIDORT-with-DO_FOCORR-on outputs** that include the higher-order coupling — more comprehensive and the production target.

Both have value. Piece A uses back-correction as a fast first-pass test; production `FullMOM` uses `ExactSFIPhase` for the substantive correction.

---

## 7. Step 3 — `AbstractRTSolver`, `AbstractRTMode`, and the unified `run_rt`

### 7.1 Unified entry point

```julia
run_rt(model::RTModel; mode::AbstractRTMode = ForwardMode(), outputs = (DirectionalTOA(),))
```

Two kwargs.

### 7.2 `AbstractRTSolver`: `FullMOM` and `ExactFirstOrderOnly`

```julia
abstract type AbstractRTSolver end

struct FullMOM{S, P, FT} <: AbstractRTSolver
    sfi_phase::P             # AbstractSFIPhase
    ss_paths::S              # :auto / :none / :all / tuple — DIAGNOSTIC for FullMOM
    l_trunc::Int
end
FullMOM(; sfi_phase=ExactSFIPhase(), ss_paths=:none, l_trunc=15) =
    FullMOM(sfi_phase, ss_paths, l_trunc)

struct ExactFirstOrderOnly{S} <: AbstractRTSolver
    ss_paths::S              # :all is the sensible default
end
ExactFirstOrderOnly(; ss_paths=:all) = ExactFirstOrderOnly(ss_paths)
```

Dispatch in `run_rt`:

```julia
function run_rt(model; mode=ForwardMode(), outputs=(DirectionalTOA(),))
    results = []
    for iBand in 1:n_bands(model)
        ctx = build_rt_context(model, mode, outputs, iBand)
        R_SFI = init_R_SFI(model, ctx)

        if model.solver isa FullMOM
            # MS solver; SFI source built per model.solver.sfi_phase dispatch
            run_ms_solver!(R_SFI, model, ctx, mode)
        end
        # ExactFirstOrderOnly: R_SFI stays zero; only postprocessing contributes

        # Standalone-solver corrections (diagnostic for FullMOM; primary for ExactFirstOrderOnly)
        foreach(p -> apply_correction!(R_SFI, p, model.solver, model, ctx), ctx.ss_paths)

        push!(results, extract_outputs(R_SFI, outputs, model, ctx))
    end
    return results
end
```

### 7.3 The four solver×mode combinations

| Solver | Mode | What it does | Use case |
|---|---|---|---|
| `FullMOM(sfi_phase=ExactSFIPhase())` | `ForwardMode` | Full RT with exact first-scatter source | Standard production |
| `FullMOM(sfi_phase=ExactSFIPhase())` | `LinMode` | Full RT + Jacobians | Inversion / retrieval |
| `FullMOM(sfi_phase=TruncatedSFIPhase())` | `ForwardMode/LinMode` | Current behavior — for regression testing |
| `ExactFirstOrderOnly` | `ForwardMode` | SS-only, all four paths | Fast forward where MS is negligible |
| `ExactFirstOrderOnly` | `LinMode` | SS-only + Jacobians | Fast Jacobians for retrieval inner loops |

### 7.4 Backwards compat

`rt_run`, `rt_run_lin`, `rt_run_multisensor` get unified as `run_rt`. Wrappers:

```julia
function rt_run(model; kwargs...)
    @warn "rt_run is deprecated; use run_rt(model)" maxlog=1
    return run_rt(model; mode=ForwardMode(), kwargs...)
end
```

---

## 8. Extension cookbook

[Same as v0.5 — adding new BRDFs, contributors, etc. via dispatch.]

---

## 9. Remaining type-gates inventory

[Same as v0.5 — list of remaining `isa` gates and migration path.]

---

## 10. Future architecture (deferred)

- `AbstractSourceType` for thermal/MW source dispatch
- Pseudo-spherical sphericity (`DO_FOCORR_OUTGOING` analog)
- Monte Carlo backend
- Adjoint-mode lin

---

## 11. Validation infrastructure

### 11.1 Reference set

[Same as v0.5: full-MS reference tables. Coulson-Dave-Sekera 1960, Garcia-Siewert 1989, Spurr 2002, Wauben-Hovenier 1992.]

### 11.2 VLIDORT comparison

Detailed in `vlidort_baseline_suite.md` (Piece B). Leverages shipped fixtures:
- `vlidort_v_test/V2p8p3_Siewert2000_validation.f90` + `saved_results/.../results_Siewert2000_validation.all` — peer-reviewed Siewert (2000) Problem IIA.
- `vlidort_s_test/2p8p3_solar_tester.f90` + `saved_results/.../results_solar_tester.all` — multi-task solver-configuration sweep enabling Task1↔Task3/4 = SS-correction validation.

### 11.3 Internal cross-checks

- `ExactFirstOrderOnly` vs `FullMOM(sfi_phase=ExactSFIPhase())` thin-limit
- `FullMOM(sfi_phase=TruncatedSFIPhase())` vs current `rt_run` (regression)
- BRDF reciprocity for Lambertian + plane-parallel + ϖ=1

### 11.4 Standalone reference

`dev_notes/exact_ss_reference/`: pol-scalar Julia + Python reference for paths 1+2+3+4 with Lambertian. Test fixture for Piece A's implementation.

### 11.5 Sanghavi regression

After Step 1B refactor, existing Cox-Munk test cases produce output equal to (with Q3 fix) pre-migration output.

---

## 12. Execution plan

### Step 1 (small refactors)
- 1A: `n_fourier_moments` dispatch
- 1B: `AbstractSSPath` family
- 1C: `AbstractHotSpot` field

### Step 1.5 (optical-properties redesign)
§5.2–5.7

### Step 2
- §6 `AbstractSFIPhase` source dispatch in `elemental.jl`
- Standalone solver kernels (Piece A)

### Step 3 (`run_rt` consolidation)

---

## 13. Open questions

**Q1**: Solver in the model — clean enough?

**Q2**: Precheck-driven defaults — framework right; thresholds derived from required-Nquad framework defensible?

**Q3**: Bug in existing Cox-Munk SS correction — intentional or oversight?

**Q4**: Mixed truncation in standalone solver — TMS for FO-equivalent paths 1+2, untruncated for paths 3+4 — defensible? Note: `ExactSFIPhase` for full-MOM uses unprimed throughout; the mixed-truncation question is for the standalone solver's FO-equivalent diagnostic mode only.

**Q5**: Canopy hotspot regression — current canopy benchmarks tuned to no-hotspot output?

**Q6**: `NoScattering` band-level dispatch acceptable?

**Q7**: `n_fourier_moments(::CoxMunk) = 2·Nstreams` empirical sanity-check.

**Q8**: `GaussQuadFullSphere` retirement.

**Q9**: Threshold values — derived from required-Nquad framework empirically calibrated, or hand-set?

**Q10**: Lin compatibility for the precheck — function-barrier pattern (precheck on primal, tuple fixed for the run).

**Q11**: Validation against published references for `ExactFirstOrderOnly` — Siewert 2000 is the natural one; for paths 3+4, the standalone reference is the only validation source.

---

## 14. References

### 14.1 Primary RT references
- Chandrasekhar 1960; Liou 2002; Hansen & Travis 1974

### 14.2 δ-M truncation
- Wiscombe 1977; Nakajima & Tanaka 1988

### 14.3 VLIDORT and FO
- Spurr 2002, 2008
- VLIDORT 2.8.3 source: `vlidort_focode/FO_VectorSS_RTCalcs_I.f90`

### 14.4 Validation references
- Coulson, Dave, Sekera 1960 (Rayleigh tables)
- Garcia & Siewert 1989; Wauben & Hovenier 1992
- Siewert 2000 (Problem IIA — peer-reviewed; in shipped VLIDORT distribution)

### 14.5 Cox-Munk
- Cox & Munk 1954

### 14.6 Canopy hotspot
- Kuusk 1991

### 14.7 vSmartMOM-internal companion documents

**Pre-existing committed in `vSmartMOM.jl`:**
- `docs/dev_notes/LINEARIZATION_BUGS.md` — linearization bug catalogue and redesign track.
- `src/CoreRT/types.jl` — current type definitions (especially `CoreScatteringOpticalProperties` at line 1052).
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl` — current linearization kernel; "Bug 19" note at line 32.

**Delivered with this design (cross-references resolve within set):**
- `standalone_ss_solver_plan.md` — Piece A implementation plan.
- `vlidort_baseline_suite.md` — Piece B validation suite.
- `dev_notes/exact_ss_reference/` — standalone reference Julia + Python + derivations.

**Existing prototype:**
- `coxmunk_surface.jl:505-546` — Sanghavi's `apply_ss_correction!` (the prototype that becomes `apply_correction!(::SurfaceDirectBeamPath{<:CoxMunkSurface}, ...)`).
- `rt_run.jl:310-315` — the existing call site.
