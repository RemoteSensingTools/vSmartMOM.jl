# vSmartMOM.jl — architecture design and Fourier work

**Document version:** v0.5
**Status:** working draft for Sanghavi review.
**Branch target:** `sanghavi-unified` (with `unified-vsmartmom` as the architectural merge base, physics layered in from `sanghavi`).
**Companion:** `LINEARIZATION_BUGS.md`, the linearization redesign track.

---

## Preamble — for Sanghavi

This is the fifth iteration of this document. The architectural skeleton (the dispatch rule, the broad strategy) has held across review rounds and v0.5 is the first one that's ready for your input on the **theory and philosophy**. The earlier rounds (Codex against the codebase, GPT against the design) caught a lot of implementation mistakes; you'll likely catch the physics mistakes I haven't yet noticed, and several of the choices below are not technically forced — they reflect judgment about how the *code* should express the *physics*. Your view of where the science is going should shape that.

**The most important change in v0.5:** Step 2 has been reframed. v0.4 described it as "subtract analytic SS from the MS solve" — which would have been wrong, since the MS solve needs the direct beam as a source for higher-order scattering paths. v0.5 reframes Step 2 as **dispatch on SFI source phase**: the multiple-scatter `r/t` operators continue to use truncated/δ-M phase functions, but the SFI source vectors `J0±` can be built from either the truncated phase (current behavior) or the untruncated phase (new). The correction `(I - M_trunc)^-1 (J_exact - J_trunc)` then propagates the exact first-scatter event through the full MOM solve. **This is structurally a different correction than VLIDORT's FO** (FO replaces escaped first-order at the sensor; this corrects the source for *all* paths beginning with a first scatter, including ones that propagate through higher-order scattering).

**The second most important change:** Step 1 has been split into Step 1 (small contained refactors) and Step 1.5 (the optical-properties algebra redesign). The split is a response to GPT's correct observation that v0.4's Step 1 conflated mechanical refactors with type-system restructuring; the latter deserves its own milestone with its own validation boundary.

### Six theory/philosophy questions for you

These are the choices that need your physics judgment before code lands. Each is one paragraph; pointers to the technical detail are references, not required reading for the philosophical pass.

**Q1 — The architectural commitment to multiple dispatch.** The proposal makes "dispatch on types at physical/workflow boundaries; `if/else` only for scalar numerical cases" the *rule* the codebase commits to. It applies to surface kernels, optical contributors, layer kinds, storage forms, output products, RT modes, RT solvers, and SFI source phase strategies — not just outer dispatch. The cost is that the type system grows: vSmartMOM v2.0 will have ~10 abstract type families dispatching at run time. The benefit is that physics-extension PRs become "define one struct, write one method." Are you on board with the rule as a hard commitment, or is there a class of physics where the freedom to branch via `if` is worth keeping? *Detail: §1.*

**Q2 — `ExactSFIPhase` corrects more paths than VLIDORT FO does, and is a different correction.** The proposal builds the SFI source vectors `J0±` from the untruncated aerosol phase (and the untruncated optical depth and ϖ — see §6.3 for why the τ/ϖ also need to be untruncated to keep the energy budget consistent). The corrected source then propagates through the *truncated* MS operator. This corrects single-scatter-then-anything paths (sun→atm→sensor, sun→atm→surface→sensor, sun→atm→atm→sensor, sun→atm→surface→atm→sensor, …) — strictly more than what VLIDORT FO corrects, but also not bit-for-bit equivalent to FO. The correction is *not* exact RT (the MS operator is still truncated), but it's a strict improvement on the truncated MS solve. The validation framing has to match: convergence to the truncated MS solve when truncation is negligible (small `f`), and improved agreement with full-MS reference tables at glint/forward-peak geometries. Does this match your physical intuition for what the correction does, and what's the right validation regime? *Detail: §6.*

**Q3 — `AerosolForwardPeakCorrection` from v0.4 is obsoleted by `ExactSFIPhase`, and v0.5 drops it.** v0.4 had `AerosolForwardPeakCorrection` as a deferred type intended to land later (a post-processing correction analog to the Cox-Munk SS correction). With `ExactSFIPhase` properly built (carrying both truncated and untruncated optical-property forms), the forward-peak correction for first-scatter events is structurally automatic — it falls out of the SFI-phase dispatch. A separate post-processing correction would only be relevant for *higher-order* path corrections (e.g., correcting truncation in second-scattering events too), which is a more sophisticated and probably unnecessary extension. **v0.5 removes `AerosolForwardPeakCorrection` from the planned `AbstractSSCorrection` types entirely.** Is there a residual post-processing component you intended this correction to do that `ExactSFIPhase` won't capture? *Detail: §3.3.*

**Q4 — `NoScattering` layer kind is a *band-level* property, not a per-layer-mixed-spectrum property.** GPT correctly caught that v0.4's `layer_kind(props) = props.ϖ ≈ 0 ? NoScattering : Scattering` was wrong because `ϖ` is per-wavelength. The existing code classifies via `maximum(τ .* ϖ) > 2eps(FT)` — *any* spectral bin scattering forces the layer onto the scattering path. v0.5 keeps this exact predicate but gives the result a type-level handle. `NoScattering` therefore fires only when the *whole band* is non-scattering across all wavelengths — which is the genuine thermal/MW regime, not a mixed-spectrum optimization. This is a smaller win than v0.4 implied, but it's still useful: it makes thermal-only RT a natural fall-out of the dispatch architecture, and it eliminates the runtime predicate from inner kernels. Does the band-level framing match how you think about the regimes? *Detail: §5.6.*

**Q5 — Canopy hotspot is silently missing today; the fix has cost.** vSmartMOM's `elemental_canopy!` uses independent Beer-Lambert (`exp(-(k_s + k_o) L)`) where it should use joint gap probability (`exp(-(k_s + k_o) L) · C_hs(α, L)`). CanopyOptics ships `KuuskHotSpot` with ForwardDiff hooks but vSmartMOM doesn't call `joint_gap_probability` anywhere. *Any current canopy benchmark of yours that compares against RAMI/PROSAIL hotspot-aware references is either wrong or was tuned to match the no-hotspot result.* The proposal defers full hotspot wiring to Step 2 (it requires the SFI-source dispatch to be in place — the hotspot is a multiplicative correction on the direct-beam-to-direct-view contribution, which lives in `J0±`), but Step 1C plumbs the field through `CanopySurface`. We'd like to know if any hotspot-tuned benchmarks exist before Step 2 lands the fix — those would need to be regenerated, not silently corrected. *Detail: §4.*

**Q6 — Validation against full-MS reference tables only; 2OS references explicitly excluded.** vSmartMOM is a full multiple-scattering code (matrix operator method, doubling-adding to convergence). Reference tables for validation must be from full-MS codes; 2OS approximations like Natraj-Spurr 2007 are *not appropriate* — comparing full-MS output against 2OS would penalize vSmartMOM for capturing higher-order scattering correctly. Proposed reference set: Coulson-Dave-Sekera 1960 (Rayleigh canonical), Garcia-Siewert 1989 (benchmark suite), Spurr 2002 VLIDORT validation appendix, Wauben-Hovenier 1992 (vector with aerosols), and the VLIDORT distribution `results/` directory. The metadata-plus-data-file harness (§11) treats validation regime as a type-dispatchable property. Does this match the reference set you've been using, or are there additional sources we should add? *Detail: §11.*

### Two narrower technical questions

These are smaller but block specific commits.

**Q7 — `n_fourier_moments(::CoxMunk) = 2 · Nstreams`.** Translating VLIDORT's `NMOMENTS = 2·NSTREAMS - 1` (Fortran inclusive loop, `2·NSTREAMS` iterations) into vSmartMOM's count convention gives `2 · Nstreams` where `Nstreams` is the *base* polar discrete-ordinate count (not `model.quad_points.Nquad`, which is augmented with VZA/SZA nodes). For `l_trunc = 60, Nstreams = 30`, the count is 60. Does this match your physical intuition? Empirical sanity: run a Cox-Munk + Rayleigh test with `n_moments` set to 60 vs. 30 vs. 100 and check convergence. *Detail: §2.*

**Q8 — `GaussQuadFullSphere` retirement, including the YAML default change.** It constructs `gausslegendre(2·Nquad)` over `[-1, 1]` and discards the lower hemisphere; the retained nodes/weights are *not* an optimal Gauss-Legendre quadrature for `[0, 1]`. Currently the *default* in `DefaultParameters.yaml:12`. Proposal: deprecate in v2.0, make `GaussLegQuad` the default, remove in v2.1. The default change is itself a behavior shift for users inheriting from the default. Do you concur, or is there a use case to preserve? *Detail: §2.8.*

### What v0.5 is *not*

- Not a fait accompli. The open-questions list (§13) flags everything we know is open.
- Not yet implemented. Current `sanghavi-unified` has none of these changes.
- Not the linearization redesign. That's a separate track. This proposal is upstream of it.
- Not exhaustive. `AbstractSourceType` for thermal/MW source dispatch is flagged as future work and explicitly *not* refactored in this round, but Step 1's design is structured so resurrecting it later is mechanical.

The rest of the document treats the technical detail. Sections 1–4 are Step 1. Section 5 is Step 1.5 (optical-properties redesign). Section 6 is Step 2 (SFI source phase). Section 7 is Step 3 (output products and modes). Sections 8–10 are extension recipes, migration debt, deferred future work. Section 11 is validation. Section 12 is the execution plan. Section 13 is open questions.

---

## Iteration log

- **v0.5** — major revision. Step 1 split into Step 1 + Step 1.5. Step 2 reframed as `AbstractSFIPhase` dispatch with dual-form optical properties (`TruncatedAndExactScatteringOpticalProperties`); the v0.4 "subtract analytic SS from MS solve" framing was physically wrong. `AbstractRTSolver` axis added (`FullMOM`, `ExactFirstOrderOnly`). `AerosolForwardPeakCorrection` from v0.4 dropped — obsoleted by `ExactSFIPhase`. GPT v0.4 review folded throughout: `layer_kind` predicate corrected to band-level max-based, `ScatteringInterface` cumulative state spelled out, `_expand_layer_rayleigh!` extraction marked as explicit migration task, `RTContext` sketch marked schematic, GPL framing softened to conservative engineering policy. `+(abs, abs)` overload and `y.G` bug fix at `types.jl:1137` confirmed in §5.2.
- **v0.4** — `AbstractLayerKind` (`Scattering`/`NoScattering`) added; storage form promoted to type parameter; three-step structure with direct-beam subsystem as Step 2. GPT v0.3 fixes folded.
- **v0.3** — Sanghavi-facing preamble; Codex v0.2 review folded; Step 1D (optical-properties algebra) added.
- **v0.2** — promoted dispatch rule to governing principle; Codex v0.1 review folded.
- **v0.1** — initial draft.

---

## Table of contents

1. [The dispatch rule (governing principle)](#1-dispatch-rule)
2. [Step 1A — per-component `n_fourier_moments`](#2-step-1a)
3. [Step 1B — `AbstractSSCorrection` strategies](#3-step-1b)
4. [Step 1C — `AbstractHotSpot` field on `CanopySurface`](#4-step-1c)
5. [Step 1.5 — optical-properties algebra and layer-kind dispatch](#5-step-1-5)
6. [Step 2 — `AbstractSFIPhase` source dispatch](#6-step-2)
7. [Step 3 — output products and RT modes/solvers as dispatched strategies](#7-step-3)
8. [Extension cookbook (aspirational)](#8-extension-cookbook)
9. [Remaining type-gates inventory (migration debt)](#9-type-gates)
10. [Future architecture (deferred)](#10-future-architecture)
11. [Validation infrastructure](#11-validation)
12. [Execution plan](#12-execution-plan)
13. [Open questions](#13-open-questions)
14. [References](#14-references)

---

<a id="1-dispatch-rule"></a>
## 1. The dispatch rule (governing principle)

### 1.1 The rule

> **Use multiple dispatch at physical and workflow boundaries; keep `if/else` only for scalar numerical cases.**

Physical and workflow boundaries become method dispatch on small concrete strategy structs:

- Surface physics (`if brdf isa CoxMunkSurface`)
- Fourier bounds per component
- Corrections — SS, hotspot, future thermal/Raman
- **Layer kinds** — `Scattering` vs. `NoScattering`, dispatched at the kernel
- **Storage forms** — `Compact` vs. `Expanded` Z, dispatched at the algebra level
- **SFI source phase** — `Truncated` vs. `Exact`, dispatched at source construction
- **RT solvers** — `FullMOM` vs. `ExactFirstOrderOnly`, dispatched at the top-level coordinator
- Output products
- RT modes (forward, lin)
- Output topology (today: separate `rt_run*.jl` files)
- Atmospheric optical contributors

`if/else` is fine — and clearer than dispatch — for scalar numerical cases:

- Fourier normalization weight (`m == 0 ? 0.5 : 1.0`)
- Numerical guards (NaN, denorm, divide-by-zero)
- Cache state checks (`if isnothing(brdf._cache)`)
- Optional benchmark/diagnostic instrumentation
- **Runtime user-data predicates** (`if any(a -> a.fᵗ > 0.01, aerosols)` — data check on values, not type gate)
- **Layer-kind classification at construction** (`maximum(τ .* ϖ) > 2eps(FT) ? Scattering() : NoScattering()` — one comparison per layer at construction time, then dispatch from there)

The rule eliminates *type and physics gates*, not all branches. Litmus test: if the branch's truth is determined at model-construction time by which types are in play, dispatch on types. If determined at evaluation time by numerical values or user toggles, branch.

**Applied at every level.** The rule is not just for outer dispatch (surface kernels, output products). It applies *all the way down* — to layer kinds, storage forms, SFI source phases, kernel specialization. The compounding payoff: the compiler can specialize through the entire stack when types are concrete, and the human reader can find each physics case by following a single dispatch chain.

### 1.2 The strategy struct shape

Strategy objects are **small concrete structs held in tuples, dispatched via single-method functions.**

**(a) Concrete tuples, not abstract vectors.** A `Vector{AbstractSSCorrection}` forces runtime dispatch in inner loops. A `Tuple{BRDFSpecularCorrection{CoxMunkSurface{Float64}}}` lets the compiler specialize.

**Honest about type stability:** when a strategy set is built conditionally and stored in a long-lived struct, end-to-end compile-time dispatch requires either parameterizing the storage struct on the strategy-tuple type (invasive) or accepting runtime dispatch through a function-barrier when unpacked. **Pragmatic answer: accept runtime dispatch at outer levels (per-band, per-run) where the strategy work itself is much more expensive than dispatch overhead. Compile-time dispatch is achieved *inside* each strategy method via concrete type parameters, not at the outer driver level.** Layer-kind, storage-form, and SFI-phase dispatch *do* propagate to compile-time because they're parameterized on `Layer`, `CoreScatteringOpticalProperties`, and the run-level invocation respectively.

**(b) Single-method functions, not switch tables.** Each strategy owns its behavior:

```julia
apply_correction!(c::BRDFSpecularCorrection, ctx)
apply_correction!(::NoSSCorrection, ctx) = nothing
build_sfi_source!(::TruncatedSFIPhase, layer, ctx)
build_sfi_source!(::ExactSFIPhase, layer, ctx)
```

**(c) Naming consistency across families.**

| Family | Trait | Apply |
|---|---|---|
| Fourier bounds | `n_fourier_moments(component, ctx)` | (count only) |
| SS corrections | `needs_correction(c, ctx)` | `apply_correction!(c, ctx)` |
| Optical contributors | `optical_properties(c, ctx)` | composes via `+` algebra |
| Layer kinds | (read off via existing predicate) | `elemental!(::Scattering, ...)` / `elemental!(::NoScattering, ...)` |
| Storage forms | (type parameter) | `expandOpticalProperties(::Compact)` returns `Expanded` |
| SFI source phase | (run kwarg) | `build_sfi_source!(::TruncatedSFIPhase, ...)` / `build_sfi_source!(::ExactSFIPhase, ...)` |
| RT solver | (run kwarg) | `solve_rt!(::FullMOM, ctx)` / `solve_rt!(::ExactFirstOrderOnly, ctx)` |
| Output products | (presence in tuple is signal) | `accumulate_output!(state, output, ctx)` |
| Hotspot models | (parameter on `CanopySurface`) | `joint_gap_probability(model, …)` (in CanopyOptics) |

### 1.3 The `RTContext` object

*(Schematic — implementation parameterizes concrete types throughout to maintain type stability. Per GPT review.)*

Strategy methods need a consistent argument convention. `RTContext` is constructed once before the Fourier loop, mutated as `m` advances.

**Lifecycle:**

- Constructed once before the Fourier loop. Carries everything that doesn't depend on `m`.
- Field `m` is mutated as the loop advances. Strategy methods read it, don't mutate it.
- `composite_layer` accumulates **vertically within one Fourier moment** (across `iz` layers), reset at `iz == 1` on the next `m`.
- **Output arrays accumulate across `m`** — these belong to `AbstractRTOutput` strategies (Step 3).
- `τ_sum_all` is m-independent.

```julia
mutable struct RTContext{FT, S, P, A, ML, AS, WS, ...}   # parameterized concretely in implementation
    model::ML                  # concrete RTModel parameterization
    iBand::Int
    m::Int
    quad_points::QuadPoints{FT}
    composite_layer::CL        # concrete composite-layer storage
    added_layer_surface::AS
    τ_sum_all::Matrix{FT}
    geometry::Geometry         # observation: vza, vaz
    source::S                  # solar today; thermal/lunar later via §2.9 invariant
    pol_type::P                # concrete subtype of AbstractPolarizationType
    arch::A                    # concrete subtype of AbstractArchitecture
    workspace::WS              # InteractionWorkspace, allocated once
    Nstreams::Int              # base polar streams (NOT augmented Nquad)
    Nquad::Int                 # augmented quadrature including VZA/SZA
    nSpec::Int                 # spectral grid size
    n_layers::Int
end
```

Exact field set is open — see §13 Q11.

### 1.4 What the rule unlocks

- **Step 1 (§§2–4)** sketches the idiom on small contained pieces. Three work items, mostly mechanical refactors.
- **Step 1.5 (§5)** applies the rule at the optical-properties level: contributor pattern, dimensioned zero, storage-form type parameter, layer-kind dispatch, dual-form properties for SFI-phase support. Substantial type-system restructure; its own validation boundary.
- **Step 2 (§6)** applies the rule to SFI source phase: dispatch `build_sfi_source!` on truncated vs. exact phase. Real solver work; validates against full-MS references at glint/forward-peak.
- **Step 3 (§7)** applies the rule to output products, RT modes, and RT solvers. Collapses 1004 lines of forked `rt_run*.jl` into one coordinator.
- **The linearization redesign** uses the same rule. Lin Jacobians are an `AbstractRTMode`, orthogonal to outputs, orthogonal to solver, orthogonal to phase.
- **Future corrections** (thermal SS, Raman SS, multiple-scatter δ-M, polarization SS for thermal) all fit the same `Abstract*Correction` shape.

The rule has compounding returns. Step 1 is the first application; everything downstream gets cheaper.

### 1.5 Branches that stay (the "fine" list)

- `m == 0` Fourier weight normalization. Numerical, stays.
- Cache initialization (`canopy._cache === nothing`). State, not type/physics dispatch.
- `maximum(τ .* ϖ) > 2eps(FT)` check at layer construction to classify `Scattering` vs. `NoScattering`. One comparison per layer; dispatch from there.
- Numerical NaN/denorm guards in kernels.

The rule eliminates type and physics gates. Not all branches.

---

<a id="2-step-1a"></a>
## 2. Step 1A — per-component `n_fourier_moments`

### 2.1 Convention: count, not order

Existing code uses `max_m` as a count: loops are `for m = 0:max_m - 1` in `rt_run.jl:208`, etc. The trait must match.

**Decision:** `n_fourier_moments` returns a count. Loop convention stays `for m = 0:n - 1`.

| Component | `n_fourier_moments` | Why |
|---|---|---|
| Lambertian (any flavor) | 1 | Only `m = 0`; no azimuthal structure |
| `RayleighScattering` | 3 | `m = 0, 1, 2` cover Rayleigh greek coefs |
| `AerosolOptics` (greek_coefs.β of length L+1) | min(L+1, 2·Nstreams) | Greek-coef tail capped at solver resolution |
| `CoxMunkSurface`, `rpvSurfaceScalar`, `RossLiSurfaceScalar` | `2 · Nstreams` | VLIDORT cap (§2.3) |
| `CanopySurface` | `2 · Nstreams` (provisional) | Proper LAD-driven derivation deferred |

### 2.2 `Nstreams` is *not* `Nquad`

`model.quad_points.Nquad` is *augmented* — `rt_set_streams.jl:40` redefines `Nquad = length(qp_μ)` after appending VZA/SZA nodes. If `n_fourier_moments(::CoxMunk, ctx) = 2 * ctx.Nquad`, the count silently varies with the user's requested viewing angles.

**Fix:** introduce a separate `Nstreams` capturing the polar discrete-ordinate count *before* VZA/SZA augmentation:

```julia
mutable struct QuadPoints{FT}
    # ... existing fields ...
    Nstreams::Int   # base solver streams, never augmented
    Nquad::Int      # may equal Nstreams (Gauss hemisphere) or include VZA/SZA (Radau)
end
```

Computed at quadrature setup time as `Nstreams = (l_trunc + 1) ÷ 2` *before* any VZA/SZA augmentation.

### 2.3 The VLIDORT count formula

VLIDORT 2.8.3, `vbrdf_sup_masters.f90:2305`:

```fortran
NMOMENTS = 2 * NSTREAMS - 1
```

Followed by Fortran-style **inclusive** loop `DO M = 0, NMOMENTS` — runs `0..2·NSTREAMS - 1`, **`2·NSTREAMS` iterations**. Translated to vSmartMOM's count convention:

```
n_fourier_moments_BRDF = 2 · Nstreams
```

For `Nstreams = (l_trunc + 1) ÷ 2`:
- Even `l_trunc = 60`: `Nstreams = 30`, count = 60.
- Odd `l_trunc = 59`: `Nstreams = 30`, count = 60.

Empirical sanity check planned (Q7).

**Rationale:** the BRDF Fourier series projects onto the discrete-ordinate basis. A polar grid of `Nstreams` cannot resolve azimuthal modes higher than `2·Nstreams - 1`. Past that aliases. Cap is set by *solver resolving power*, not surface complexity.

### 2.4 The grids (kept separate)

| Quantity | Role | Source | Notes |
|---|---|---|---|
| `Nstreams` | Solver polar discrete-ordinate count, base | `rt_set_streams.jl` (new field) | Used for VLIDORT cap |
| `Nquad` | Augmented quadrature including VZA/SZA | `rt_set_streams.jl:40` | Used for kernel sizing, Z-matrix dimensions |
| `nQuad_BRDF` | Azimuth integration grid for each Fourier coefficient | Hard-coded `100` in surface files | Quadrature accuracy, deferred to per-surface field |
| `n_fourier_moments` | Count of orders the loop runs over | New trait | The trait being added |

### 2.5 Consumers reading the wrong value, and `max_m` semantics

**Two-part fix:**

**Part 1 — the consumer read.** `rt_run.jl:97` calls `get_max_m(model)` returning the scalar `model.solver.max_m`, not the per-band vector. Replace with `n_fourier_moments(model, iBand)` in `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`.

**Part 2 — the `max_m` override semantics.** YAML schema currently requires `radiative_transfer.max_m` as an integer; no way to distinguish "user explicitly capped" from "default value." Old configs with `max_m: 30` would still silently cap Cox-Munk to 30 even though the trait says 60.

Fix: change `max_m` to `Union{Nothing, Int}`:

```julia
struct SolverConfig{...}
    max_m::Union{Nothing, Int}   # nothing = "use trait"; Int = "user override cap"
end

function n_fourier_moments(model::RTModel, iBand::Int)
    ctx = (; Nstreams = model.quad_points.Nstreams)
    n = n_fourier_moments(model.optics.rayleigh, ctx)
    for a in model.optics.aerosols.aerosol_optics[iBand]
        n = max(n, n_fourier_moments(a, ctx))
    end
    n = max(n, n_fourier_moments(model.surfaces[iBand], ctx))
    n = min(n, 2 * ctx.Nstreams)
    return isnothing(model.solver.max_m) ? n : min(n, model.solver.max_m)
end
```

YAML default becomes `max_m: null` (or omits the field; loader treats absence as `nothing`). **Old configs with the historical default integer must be regenerated** — flag in changelog.

### 2.6 Proposed trait methods

```julia
n_fourier_moments(::RayleighScattering,        ctx) = 3

n_fourier_moments(::LambertianSurfaceScalar,   ctx) = 1
n_fourier_moments(::LambertianSurfaceSpectrum, ctx) = 1
n_fourier_moments(::LambertianSurfaceLegendre, ctx) = 1
n_fourier_moments(::LambertianSurfaceSpline,   ctx) = 1

# VLIDORT vbrdf_sup_masters.f90:2305 invariant
n_fourier_moments(::CoxMunkSurface,           ctx) = 2 * ctx.Nstreams
n_fourier_moments(::rpvSurfaceScalar,         ctx) = 2 * ctx.Nstreams
n_fourier_moments(::RossLiSurfaceScalar,      ctx) = 2 * ctx.Nstreams
n_fourier_moments(::CanopySurface,            ctx) = 2 * ctx.Nstreams

n_fourier_moments(a::AerosolOptics,           ctx) =
    min(length(a.greek_coefs.β), 2 * ctx.Nstreams)
```

### 2.7 Knock-on simplifications

- **`apply_correction!` reads `n_fourier_moments` from `ctx`**, never recomputes — guarantees the correction subtracts the same Fourier sum that was added.
- **Lambertian-only bands skip `m > 0`**. Loop runs once.
- **`vSmartMOM_Parameters.max_m` becomes a user-override knob.**

### 2.8 `GaussQuadFullSphere` retirement (expanded scope)

Reading `rt_set_streams.jl:63-86`: `GaussQuadFullSphere` constructs `gausslegendre(2·Nquad)` over `[-1, 1]` and discards the lower hemisphere. Retained nodes/weights are *not* an optimal Gauss-Legendre quadrature for `[0, 1]`. Strictly worse than `GaussLegQuad`.

**Files touched by retirement:**
- `src/CoreRT/types.jl:117` — type definition
- `src/CoreRT/CoreRT.jl:154` — export
- `src/CoreRT/tools/rt_set_streams.jl:53,63` — dispatched method
- `src/CoreRT/DefaultParameters.yaml:12` — **currently the default**
- `test/test_parameters/ThreeBandsParameters.yaml:73`
- `docs/src/pages/concepts/02_rt_theory.md:76`
- `docs/src/pages/IO/Schema.md:32,84`
- `docs/src/pages/IO/Overview.md:33,55`
- `docs/src/pages/geoschem_integration.md:79`
- `src/IO/Parameters.jl:94,99` — YAML loader

**Retirement plan (v2.0):**

1. `Base.depwarn` in the `rt_set_streams(::GaussQuadFullSphere, ...)` method.
2. `DefaultParameters.yaml` default changes to `GaussLegQuad()`. **Behavior change for users inheriting from default.** Document in changelog.
3. `ThreeBandsParameters.yaml` test fixture updates.
4. Docs updates (4 files).
5. YAML loader keeps the string mapping with deprecation warning.
6. Audit reference benchmarks; regenerate with `GaussLegQuad`.

**v2.1:** remove the type, the `rt_set_streams` method, and the YAML mapping.

Asks Sanghavi for confirmation (Q8) given the YAML default change.

### 2.9 Source-type forward-compatibility (deferred refactor, no rename)

`DNI` and `SFI` types exist in `CoreRT/types.jl:120-126` as `AbstractSourceType` subtypes but are essentially vestigial — the `SFI::Bool` parameter does the real work, and `"DNI"`/`"SFI"` strings aren't in the YAML schema. v0.4 had proposed renaming; v0.5 keeps the existing names because the rename was churn without removing the boolean.

The forward-compat invariant remains:

> *Any code touching solar geometry or source-term assembly in Step 1 must funnel through `ctx.source`, not through hardcoded `μ₀`/`I₀` references or positional arguments.*

When thermal RT eventually lands (post-v2.0), the path is: add `ThermalSource` type, add `source_term(::ThermalSource, ctx)` methods, retire `SFI::Bool` and the existing `DNI`/`SFI` types together with a proper migration. No re-architecture needed — Step 1's strategy methods already go through `ctx.source`.

The `AbstractLayerKind` work in §5.6 (`NoScattering` dispatch) is the *complementary* future-compat piece.

---

<a id="3-step-1b"></a>
## 3. Step 1B — `AbstractSSCorrection` strategies

### 3.1 The current state

`rt_run.jl:310-315`:

```julia
if brdf isa CoxMunkSurface && SFI
    @timeit "SS Correction" apply_ss_correction!(
        R_SFI, brdf, pol_type, vza, vaz, μ₀,
        Array(τ_sum_all[:,end]), max_m, nSpec)
end
```

Three problems:

1. **`isa` gate.** Adding any new BRDF needing the same correction means another `if`.
2. **Inconsistent across entry points.** Present in `rt_run.jl`, missing from `rt_run_lin.jl` and `rt_run_multisensor.jl`. Missing from `rt_run_lin.jl` is a silent Jacobian bug at glint geometries.
3. **`if brdf isa CoxMunkSurface` hides physics.** Currently the only BRDF post-processing correction; making it surface-agnostic enables future BRDFs (Cox-Munk-Bréon, etc.) that need similar treatment.

### 3.2 What VLIDORT does

VLIDORT factors single-scattering into **FO — First-Order Code**. Peer of the multiple-scatter solver. Control surface in `vlidort_inputs_def.f90:444-481`:

- `DO_FOCORR` — master on/off
- `DO_FOCORR_EXTERNAL` — external pre-computed exact-SS hook
- `DO_FOCORR_NADIR` / `DO_FOCORR_OUTGOING` — sphericity treatment
- `DO_SSCORR_TRUNCATION` — whether SS correction also accounts for δ-M truncation
- `DO_SSCORR_USEFMAT` — exact phase function from F-matrix vs. Greek coefficients

Step 1B adopts the **post-processing surface correction** dimension of FO. The **source-construction phase correction** dimension (the `ExactSFIPhase` work) lives in Step 2 (§6) — these are different architectural seams attacking different truncation effects.

### 3.3 Strategy hierarchy

```julia
abstract type AbstractSSCorrection end

struct NoSSCorrection <: AbstractSSCorrection end

struct BRDFSpecularCorrection{S<:AbstractSurfaceType} <: AbstractSSCorrection
    surface::S
    nQuad_BRDF::Int
end

struct ExternalSSCorrection{F} <: AbstractSSCorrection
    callback::F
end
```

**`AerosolForwardPeakCorrection` is removed from v0.5.** v0.4 had this as a deferred `AbstractSSCorrection` type intended to land later. With `ExactSFIPhase` (§6) properly built — carrying both truncated and untruncated optical-property forms via the dual-form contributor pattern (§5) — the forward-peak correction for first-scatter events is structurally automatic at source construction. A separate post-processing correction would only be relevant for higher-order corrections, which is a more sophisticated and probably unnecessary extension. **Q3** asks Sanghavi whether there's a residual post-processing component the v0.4 framing intended that this misses.

Per-surface trait controlling defaults:

```julia
has_specular_peak(::AbstractSurfaceType) = false
has_specular_peak(::CoxMunkSurface)      = true
```

### 3.4 Storage and the call site

**Storage decision:** `ss_corrections` does **not** live as a field on `SolverConfig`. Instead, computed via `default_ss_corrections(model)` at `run_rt` time, with user override via kwarg:

```julia
run_rt(model;
       solver=FullMOM(),
       mode=ForwardMode(),
       sfi_phase=TruncatedSFIPhase(),       # see §6
       outputs=(DirectionalTOA(), DirectionalBOA()),
       ss_corrections=default_ss_corrections(model),
       iBand=1)
```

Reasons:
- Avoids type-stability dance of storing concrete tuples in a long-lived struct
- Cleaner separation: model = atmospheric state; `run_rt` kwargs = run options
- User override is one kwarg, not a model mutation

Default selection:

```julia
function default_ss_corrections(model)
    map(eachindex(model.aerosol_optics)) do b
        cs = ()
        if has_specular_peak(model.surfaces[b])
            cs = (cs..., BRDFSpecularCorrection(model.surfaces[b], 100))
        end
        isempty(cs) ? (NoSSCorrection(),) : cs
    end
end
```

Dispatch:

```julia
apply_correction!(::NoSSCorrection, ctx) = nothing

function apply_correction!(c::BRDFSpecularCorrection{<:CoxMunkSurface}, ctx)
    # Existing coxmunk_surface.jl:505-546 body, parameterized over ctx.
    # Reads ctx.n_moments_band — same value the loop used.
end

apply_correction!(c::ExternalSSCorrection, ctx) = c.callback(ctx)
```

`rt_run.jl:310-315` becomes:

```julia
SFI && @timeit "SS corrections" foreach(c -> apply_correction!(c, ctx),
                                        ctx.ss_corrections[iBand])
```

For lin: `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::LinContext)` is **a separate method** updating both `R_SFI` and `Ṙ_SFI`. Real implementation work — the call site is one line; the dispatched method has its own derivative math. Multisensor follows the same pattern.

### 3.5 What this gives the linearization redesign

Each correction is its own self-contained unit of differentiation. Wrap `apply_correction!` in a ForwardDiff seed → get `∂(correction)/∂(params)` → add to chain rule via `lin_added_layer_all_params.jl`. With strategy structs, the hack becomes a typed unit the linearization framework can see.

---

<a id="4-step-1c"></a>
## 4. Step 1C — `AbstractHotSpot` field on `CanopySurface`

### 4.1 What CanopyOptics ships

`CanopyOptics.jl/src/canopy_structure/hotspot.jl`:

- `abstract type AbstractHotSpot{FT<:Real}`
- `struct NoHotSpot{FT}` — independent gaps, `C_hs = 1`
- `struct KuuskHotSpot{FT}` with size parameter `h`
- `joint_gap_probability(model, k_s, k_o, μ_s, μ_o, dϕ, L)` — `exp(-(k_s+k_o)L) · C_hs(L)`
- ForwardDiff hooks (`_hotspot_value(x::ForwardDiff.Dual)`)

**Abstraction shipped. Implementation correct. Linearization anticipated. Missing: vSmartMOM-side wiring.**

### 4.2 Why this is *not* a subtype of `AbstractSSCorrection`

The two corrections share a *design pattern* but live on different seams:

| | `AbstractSSCorrection` | `AbstractHotSpot` |
|---|---|---|
| **When** | After Fourier loop completes | Inside the canopy SFI source assembly |
| **How** | Additive | Multiplicative on joint gap probability |
| **What it patches** | Fourier truncation losing specular-peak resolution | Beer-Lambert assuming independent gaps |

Conceptual similarity captured by **shared dispatch idiom**, not shared type.

### 4.3 The wiring: physical reasoning and seam choice

**Physical reasoning:** the hotspot is a correlated-path-pair phenomenon. `P_so(L)` describes the joint probability that the *same gap* is open along both sun and view paths. When sun and view coincide (`α → 0`), every gap that lets the photon in also lets it out — paths perfectly correlated.

**This correlation only exists for single-scattering events observed directly without further interaction.** Once a photon scatters once into the diffuse field, its outgoing direction decorrelates. So the hotspot correctly applies to:
- Direct-solar → direct-view contribution

It does **not** apply to:
- Direct-solar → diffuse field → view contribution
- Diffuse-source → diffuse field → view contribution

**The natural seam is the SFI source dispatch (Step 2).** With `build_sfi_source!` dispatched on phase strategy, the canopy hotspot becomes another dispatch dimension on the same hook: `build_sfi_source!(canopy_layer, sfi_phase, hotspot, ctx)` reads the hotspot from the canopy surface's parameter and applies the multiplicative correction at source construction.

This is a clean unification — three previously-separate concerns (phase truncation correction, surface BRDF specular correction, canopy hotspot correction) all live on the same `build_sfi_source!` dispatch point, each as its own dispatch dimension.

### 4.4 Step 1C does the field, Step 2 does the wiring

**In Step 1C (this round):**

1. Add `hotspot::HS` parameter to `CanopySurface`, where `HS<:AbstractHotSpot{FT}`. Default `NoHotSpot{FT}()`. **Parameterized on hotspot type** so dispatch is compile-time when methods reach `surface.hotspot`.
2. Plumb through `model_from_parameters` and YAML so `CanopySurface(..., hotspot=KuuskHotSpot(0.1))` is constructable.
3. Tests verifying field is constructable but runtime behavior is bit-for-bit identical to today (`NoHotSpot` is the default and no code path acts on it yet).
4. Document explicitly that wiring is deferred to Step 2.

~30 lines.

**In Step 2 (next round, §6):**

1. The `build_sfi_source!` dispatch automatically picks up the hotspot from `surface.hotspot` for canopy layers.
2. Validate against RAMI HOM26/HOM27.

⚠️ **For Sanghavi (Q5):** any current canopy benchmark tuned to match the no-hotspot result needs flagging before Step 2 lands.

---

<a id="5-step-1-5"></a>
## 5. Step 1.5 — optical-properties algebra and layer-kind dispatch

This is its own milestone, peer in scope and risk to Step 2 — separated from Step 1 per GPT review because it restructures the core type system and kernel dispatch surface. Six pieces: documentation + bug fix (5.2), dimensioned zero (5.3), contributor pattern (5.4), dual-form properties for SFI-phase support (5.5), storage form as type parameter (5.6), `NoScattering` layer kind (5.7).

### 5.1 Motivation

`constructCoreOpticalProperties` in `compEffectiveLayerProperties.jl:11-65` hardcodes Rayleigh + aerosols + absorption as the contributors. Two consequences:

1. New contributors (clouds, ice, sea-salt) require editing the file. Rayleigh-vs-Cabannes branch is also a special case here.
2. The thermal/MW case has no Rayleigh, but the current code uses Rayleigh as the loop initializer.

Plus a third reason discovered in the SFI-phase analysis (§6): **the optical-properties pipeline needs to carry both truncated and untruncated forms for `ExactSFIPhase` to work correctly.** The current pipeline discards untruncated information after δ-M truncation is applied. Step 1.5 keeps both forms.

The operator algebra (`+` and `*`) on `CoreScatteringOpticalProperties` already does most of the right work — contributors *combine via `+`*. What's missing: making the contributor *list* the dispatched axis, formalizing the algebra, fixing existing bugs, supporting dual forms, and adding layer-kind dispatch.

### 5.2 The existing algebra (what the docs need to say, plus the bug)

Five overloads in `types.jl:1083-1137`:

- `+(scattering, scattering)`: τ-and-ϖ-weighted average of phase functions, with `wx == 0` short-circuit at line 1099
- `*(scattering, scattering)`: vertical/band concatenation via `cat`
- `+(scattering, absorption)` and `+(absorption, scattering)`: thicken layer, leave Z untouched (commutative)
- `*(scalar, scattering)`: scale τ only, leave ϖ and Z untouched

**Bug at types.jl:1137** (GPT caught this): the scalar `*` overload references `y.G` but `y` is `CoreScatteringOpticalProperties` which has no `G` field — only `CoreDirectionalScatteringOpticalProperties` does. Real bug; calling `scalar * CoreScatteringOpticalProperties(...)` would error. Never exercised because nothing actually calls scalar multiplication on the non-directional type. **Fix as part of Step 1.5.1** — remove the `y.G` reference.

**Storage-form duality** (the part that's not formally documented):

- `CoreScatteringOpticalProperties` may carry Z as `(n_pq, n_pq, 1)` (compact, broadcast across the spectral grid) or `(n_pq, n_pq, nSpec)` (expanded, per-wavelength).
- All operators accept both via broadcasting; `expandOpticalProperties` converts compact → expanded explicitly when kernels need it.

Algebraic properties:
- `+` is commutative and associative on scattering+scattering, commutative on scattering+absorption (no identity element today)
- `*` is associative but **not** commutative — vertical/band stacking, order matters

**Missing from existing algebra (added in Step 1.5.1):**

- `+(absorption, absorption)` — needed once per-species absorption contributors exist. Trivial: `CoreAbsorptionOpticalProperties(x.τ .+ y.τ)`.

### 5.3 Dimensioned zero

To start the contributor loop with empty `props` (so thermal/MW with no Rayleigh works), the algebra needs a zero element.

**Dimensioned zero — a `CoreScatteringOpticalProperties` value with τ = 0.** Reasons:

1. **The existing `+` short-circuit at line 1099 already implements zero-element behavior.** When `x.τ = 0`, `wx = 0`, the short-circuit returns y's data.
2. **No new sentinel type** — would require ~10 new method signatures.
3. **Compact form keeps cost low.** Z is `(n_pq, n_pq, 1)` for the zero, ~256 KB per band per Fourier moment.
4. **Zero participates only in `+`, never in `*`.**

**Correct dimensions:**

```julia
function zero_optical_properties(ctx)
    n_pq = ctx.pol_type.n * ctx.Nquad   # NOT Nstreams — Z uses augmented quadrature
    CoreScatteringOpticalProperties(
        τ   = zeros(ctx.FT, ctx.nSpec),  # NOT n_layers — τ is per-wavelength
        ϖ   = zero(ctx.FT),               # placeholder; never read because τ=0
        Z⁺⁺ = zeros(ctx.FT, n_pq, n_pq, 1),
        Z⁻⁺ = zeros(ctx.FT, n_pq, n_pq, 1),
    )
end
```

**Pure-absorption case** (thermal IR with only absorption, no scattering anywhere): the algebra produces a degenerate `CoreScatteringOpticalProperties` with τ from absorption, ϖ=0, Z=0 from the zero. **This is not the same as `CoreAbsorptionOpticalProperties`.** §5.7's `NoScattering` layer-kind dispatch correctly classifies this case via the existing max-based predicate, and the kernel dispatches to a closed-form Beer-Lambert path that ignores Z. So the degenerate scattering type is *numerically equivalent* to absorption-only behavior, achieved via dispatch. **One pipeline, dispatched kernels handle the regime difference.**

### 5.4 The contributor pattern

```julia
abstract type AbstractOpticalContributor end

# Single trait:
optical_properties(c::AbstractOpticalContributor, ctx) -> AbstractOpticalProperties

struct RayleighContributor{FT, S} <: AbstractOpticalContributor
    greek_source::S         # carries Rayleigh-vs-Cabannes choice as data
    τ_rayl::Vector{FT}
    ϖ::FT
end

struct AerosolContributor{FT, M} <: AbstractOpticalContributor
    aerosol_optics::M
    τ_aer::Vector{FT}
    apply_delta_m::Bool
end

struct AbsorptionContributor{FT} <: AbstractOpticalContributor
    τ_abs::Vector{FT}
end
```

**Scope:** the contributor pattern is **for atmospheric layers only.** `CoreDirectionalScatteringOpticalProperties` (used by canopy in `canopy_surface.jl:408`) has no `+`/`*` algebra and is its own pipeline. Step 1.5 doesn't refactor canopy's optical-properties construction.

### 5.5 Dual-form optical properties for SFI-phase support

This is the structural change required to support `ExactSFIPhase` (§6). The MS solver needs truncated quantities for energy-budget consistency; the exact-SFI source needs untruncated quantities. The current pipeline discards untruncated information after δ-M is applied.

**Solution:** each contributor produces a `TruncatedAndExactScatteringOpticalProperties` value carrying both forms:

```julia
struct TruncatedAndExactScatteringOpticalProperties{T, E}
    truncated::T   # CoreScatteringOpticalProperties or CoreAbsorptionOpticalProperties
    exact::E       # same — equal to truncated when no truncation applies
end
```

The aerosol contributor:

```julia
function optical_properties(c::AerosolContributor, ctx)
    Z⁺⁺_full, Z⁻⁺_full = compute_Z_moments(ctx.pol_type, ctx.μ,
                                            c.aerosol_optics.greek_coefs,
                                            ctx.m, arr_type=ctx.arr_type)
    if c.apply_delta_m
        (; fᵗ, ω̃) = c.aerosol_optics
        τ_trunc = (1 - fᵗ * ω̃) .* c.τ_aer
        ϖ_trunc = (1 - fᵗ) * ω̃ / (1 - fᵗ * ω̃)
        Z⁺⁺_trunc, Z⁻⁺_trunc = compute_Z_moments(ctx.pol_type, ctx.μ,
                                                  truncate_greek(c.aerosol_optics.greek_coefs, fᵗ),
                                                  ctx.m, arr_type=ctx.arr_type)
        TruncatedAndExactScatteringOpticalProperties(
            truncated = CoreScatteringOpticalProperties(τ_trunc, ϖ_trunc, Z⁺⁺_trunc, Z⁻⁺_trunc),
            exact     = CoreScatteringOpticalProperties(c.τ_aer, c.aerosol_optics.ω̃, Z⁺⁺_full, Z⁻⁺_full),
        )
    else
        # No truncation — both forms are identical
        same = CoreScatteringOpticalProperties(c.τ_aer, c.aerosol_optics.ω̃, Z⁺⁺_full, Z⁻⁺_full)
        TruncatedAndExactScatteringOpticalProperties(truncated=same, exact=same)
    end
end
```

Rayleigh and absorption contributors produce dual forms with `truncated == exact` (no truncation issue).

The `+` algebra extends to dual forms:

```julia
function +(x::TruncatedAndExactScatteringOpticalProperties,
           y::TruncatedAndExactScatteringOpticalProperties)
    TruncatedAndExactScatteringOpticalProperties(
        truncated = x.truncated + y.truncated,
        exact     = x.exact     + y.exact,
    )
end
```

Each form combined independently with the existing `CoreScatteringOpticalProperties` algebra. Storage-form duality (5.6) and layer-kind dispatch (5.7) apply to each side independently.

**Memory cost:** ~2× for aerosol Z matrices in atmospheres with δ-M truncation. Compact-form Z is small (~256 KB per band per moment), so not a serious burden, but not free either. When `apply_delta_m=false` for all contributors, both forms share storage.

**Consumers:**

- MS solver and `TruncatedSFIPhase` source construction read `props.truncated`
- `ExactSFIPhase` source construction reads both `props.truncated` and `props.exact`

This data structure is the precondition for Step 2.

### 5.6 Storage form promoted to type parameter

`expandOpticalProperties` converts compact `(n_pq, n_pq, 1)` to expanded `(n_pq, n_pq, nSpec)`. Today the *type* doesn't track which form is in play.

**Promote to type parameter:**

```julia
abstract type AbstractStorageForm end
struct Compact <: AbstractStorageForm end       # Z is (n_pq, n_pq, 1)
struct Expanded <: AbstractStorageForm end      # Z is (n_pq, n_pq, nSpec)

struct CoreScatteringOpticalProperties{SF<:AbstractStorageForm, FT, FT2, FT3} <: AbstractOpticalProperties
    τ::FT
    ϖ::FT2
    Z⁺⁺::FT3
    Z⁻⁺::FT3
end
```

Operators dispatch on storage form (4 `+` overloads instead of 1 with broadcast):

```julia
+(x::CoreScatteringOpticalProperties{Compact},  y::CoreScatteringOpticalProperties{Compact})  # stays compact
+(x::CoreScatteringOpticalProperties{Compact},  y::CoreScatteringOpticalProperties{Expanded}) # promotes
+(x::CoreScatteringOpticalProperties{Expanded}, y::CoreScatteringOpticalProperties{Compact})  # mirror
+(x::CoreScatteringOpticalProperties{Expanded}, y::CoreScatteringOpticalProperties{Expanded}) # per-wavelength

expandOpticalProperties(x::CoreScatteringOpticalProperties{Compact}, nSpec) :: CoreScatteringOpticalProperties{Expanded}
expandOpticalProperties(x::CoreScatteringOpticalProperties{Expanded}, _) = x
```

What this gives: compiler specializes kernel code; kernels requiring expanded form can require it at the type level; `expandOpticalProperties` becomes a type-level transformation.

### 5.7 `NoScattering` layer kind (band-level)

The current code already uses `maximum(τ .* ϖ) > 2eps(FT)` as the predicate for whether a layer scatters (`compEffectiveLayerProperties.jl:87`). The classification is correct; what's missing is exposing it as type-level dispatch.

```julia
abstract type AbstractLayerKind end
struct Scattering   <: AbstractLayerKind end
struct NoScattering <: AbstractLayerKind end

# Layer carries kind as type parameter (zero runtime cost)
struct Layer{LK<:AbstractLayerKind, ...}
    optical_properties::AbstractOpticalProperties
    # ... other fields ...
end
```

Layer kind determined at construction using the **existing band-level max-based predicate**:

```julia
function layer_kind(props::CoreScatteringOpticalProperties{<:Any, FT}) where FT
    maximum(props.τ .* props.ϖ) > 2eps(FT) ? Scattering() : NoScattering()
end

layer_kind(props::CoreAbsorptionOpticalProperties) = NoScattering()
```

**`NoScattering` fires only when the *whole band* is non-scattering across all wavelengths.** Per Q4 — this is a smaller win than v0.4 implied (it's not a per-spectrum optimization), but it's still useful for:

- The genuine thermal/MW regime (no Rayleigh, no scattering aerosols, no clouds across the band)
- Eliminating the runtime `maximum(τ .* ϖ) > 2eps(FT)` predicate from inner kernels
- Making thermal-only RT a natural fall-out of the dispatch architecture

Kernels dispatch:

```julia
function elemental!(::Scattering, layer, ctx)
    # full doubling/elemental as today
end

function elemental!(::NoScattering, layer, ctx)
    # closed-form: r = 0, t = exp(-τ/μ) * I
    # if thermal source: j = B(T) * (1 - exp(-τ/μ))
end

function doubling!(::Scattering, layer, ctx)
    # full doubling
end

function doubling!(::NoScattering, layer, ctx)
    # trivial: t' = t·t, r' = 0
end
```

**Interaction with existing `ScatteringInterface` cumulative state** (per GPT): the existing code uses `scattering_interface = get_scattering_interface(prev_state, current_scatter, iz)` to track whether the composite has scattered up to this point. This is a state machine, not a per-layer property.

`LayerKind` feeds `get_scattering_interface` via a predicate:

```julia
scatter_predicate(::Scattering)   = true
scatter_predicate(::NoScattering) = false

scattering_interface = get_scattering_interface(prev_interface,
                                                 scatter_predicate(layer_kind(lods[iz])),
                                                 iz)
```

The lin/SS/multisensor paths that currently have hardcoded `scatter = true` cases need auditing as part of Step 1.5.7 — they should query `LayerKind` instead.

### 5.8 Open design choice: `fScattRayleigh` extraction

The current code computes `fScattRayleigh = [Array(rayl[i].τ ./ combo[i].τ) ...]` as a Rayleigh-specific quantity used downstream by `_expand_layer_rayleigh!` for inelastic Raman.

**This needs a real migration plan, not a one-line answer.** Two options:

- **(A) Trait-based.** `f_scattering(c::AbstractOpticalContributor, props, ctx) -> Vector{FT}` defined per contributor. Generic.
- **(B) Push to use site.** Compute inside `_expand_layer_rayleigh!` directly. Requires plumbing changes — the use site doesn't currently carry both the Rayleigh contributor and the assembled props.

Per GPT review, Option B is *not* the low-scaffolding answer it appears to be. v0.5 marks this as an explicit migration task in Step 1.5 with its own substep (1.5.8) and defers the choice between A and B to that substep's implementation.

### 5.9 Operator style preserved

The existing `+` and `*` operators stay (modulo the `y.G` bug fix at line 1137). The contributor pattern adds *one new abstract type* (`AbstractOpticalContributor`) and *one new trait* (`optical_properties`). Storage-form dispatch refines the existing operators with explicit form-promotion methods. Layer-kind dispatch is at the kernel level. Dual-form properties wrap the existing single-form types.

No semantic change to `+` or `*` for standard cases; the algebra is *richer and more typed*, not *different*.

---

<a id="6-step-2"></a>
## 6. Step 2 — `AbstractSFIPhase` source dispatch

### 6.1 Motivation

vSmartMOM's MS operators `r/t` use truncated/δ-M phase moments — that's correct for energy conservation in the MS solve. But the SFI source vectors `J0±` are also currently built from truncated phase, which means the *first solar scattering event* is approximated even though we have the exact phase function available.

This is the same architectural gap that VLIDORT FO addresses, but the right vSmartMOM-native abstraction is different from VLIDORT's. vSmartMOM's `J0±` are already in `elemental.jl` and are propagated by the MOM doubling/adding algebra. Fixing them at construction (rather than adding a separate FO solver) corrects more paths than VLIDORT FO does.

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

This corrects the contribution of the first solar scattering event for **all paths beginning with that scatter**:

```
sun → atm scatter → sensor                               (FO-equivalent first-order)
sun → atm scatter → surface → sensor                     (FO-equivalent first-order)
sun → atm scatter → atm scatter → sensor                 (higher-order; FO doesn't fix)
sun → atm scatter → surface → atm scatter → sensor       (higher-order; FO doesn't fix)
...
```

**This is hybrid:** the *first solar scattering event* is exact (untruncated phase, untruncated optical depth, untruncated ϖ); *later propagation* uses the truncated MS operator. Strictly more paths corrected than VLIDORT FO; not the same as exact full RT (the MS operator stays truncated).

### 6.3 The energy-budget normalization (Q2 detail)

This is the technically delicate part.

When δ-M truncates the phase function:
- Phase becomes `P_trunc(Θ) = (P(Θ) - f δ(Θ))/(1-f)`, normalized to `4π`
- Optical depth becomes `τ' = (1 - fϖ) τ` (less attenuation, because forward-truncated part is treated as unscattered)
- Single-scattering albedo becomes `ϖ' = (1-f) ϖ / (1 - fϖ)`

The truncated SFI source for layer `iz` uses *all three* primed quantities to keep energy consistent.

For **`ExactSFIPhase` to work correctly**, `J_exact` must be built from *all three unprimed quantities* — `P`, `τ`, `ϖ`. The correction `J_exact - J_trunc` then has the right energy budget: it's the difference between "what the radiance would have been with no truncation anywhere" and "what truncated single-scatter gives." Mixing forms (e.g., untruncated `P` with primed `ϖ'`) double-counts the scattering integral and breaks the energy budget.

**This is why §5.5's dual-form contributor pattern is required.** Each contributor produces `TruncatedAndExactScatteringOpticalProperties` carrying both forms together; `ExactSFIPhase` reads both, builds `J_trunc` from `.truncated` and `J_exact` from `.exact`, and stores the difference for propagation.

### 6.4 Implementation

```julia
abstract type AbstractSFIPhase end

struct TruncatedSFIPhase <: AbstractSFIPhase end   # current behavior
struct ExactSFIPhase <: AbstractSFIPhase end

# Dispatch in elemental.jl source construction:
function build_sfi_source!(::TruncatedSFIPhase, layer, ctx)
    # Use layer.optical_properties.truncated
    # Build J0± as before
end

function build_sfi_source!(::ExactSFIPhase, layer, ctx)
    # Read both .truncated and .exact
    # Build J_trunc and J_exact at the actual viewing geometry
    # Store (J_exact - J_trunc) in the source vector for propagation
end
```

The MS operator is unchanged — it always uses `.truncated`. Only source construction varies.

### 6.5 `AbstractRTSolver` and `ExactFirstOrderOnly`

Step 2 introduces a third orthogonal axis: the RT solver itself.

```julia
abstract type AbstractRTSolver end
struct FullMOM <: AbstractRTSolver end
struct ExactFirstOrderOnly <: AbstractRTSolver end
```

Combinations:

| Solver + Phase | Result |
|---|---|
| `FullMOM` + `TruncatedSFIPhase` | Current standard solve |
| `FullMOM` + `ExactSFIPhase` | Exact first solar event propagated through MOM |
| `ExactFirstOrderOnly` + `ExactSFIPhase` | Fast exact first-order/direct-beam solution only |
| `ExactFirstOrderOnly` + `TruncatedSFIPhase` | (no sensible interpretation; runtime error or restrict types) |

`ExactFirstOrderOnly` is genuinely useful for: quick-look simulations, debugging SFI/source-term physics, validating exact phase evaluation, FO-only comparisons against VLIDORT, high-speed approximations where higher-order scattering is negligible.

**v0.5 implements `FullMOM` + both phases.** `ExactFirstOrderOnly` is deferred to a follow-on step (Step 4 or similar — could be post-v2.0). The dispatch infrastructure (the `AbstractRTSolver` type and its single-method functions) lands in Step 2; the `ExactFirstOrderOnly` body lands when there's demand.

### 6.6 Canopy hotspot wiring (deferred from Step 1C)

With `build_sfi_source!` dispatched and the canopy surface carrying `hotspot::HS`, the hotspot multiplication is another dispatch dimension on the same hook:

```julia
function build_sfi_source!(sfi_phase, canopy_layer, ctx)
    base_source = compute_base_sfi_source(sfi_phase, canopy_layer, ctx)
    apply_hotspot!(base_source, canopy_layer.surface.hotspot, ctx)
    return base_source
end

apply_hotspot!(source, ::NoHotSpot, ctx) = nothing  # default no-op
function apply_hotspot!(source, h::KuuskHotSpot, ctx)
    # Multiply direct-beam-to-direct-view contribution by C_hs
end
```

The hotspot is independent of `sfi_phase` choice — it's a gap-probability correction, separate from the phase-function-truncation correction. Both apply.

### 6.7 Validation

**Per Q2 framing:** the validation is *not* bit-for-bit VLIDORT FO parity in `FullMOM` mode (different correction).

**Validation regime for `FullMOM + ExactSFIPhase`:**

1. **Convergence to truncated.** When `f → 0` (smooth phase function, no δ-M needed), `J_exact → J_trunc` and the result must converge to `FullMOM + TruncatedSFIPhase`. Bit-for-bit when `apply_delta_m = false` for all contributors.
2. **Energy conservation.** For an atmosphere with a single forward-peaked aerosol layer (e.g., `f = 0.6, ϖ = 0.95, τ = 1.0`): compute total reflected energy two ways and compare to a full-MS-no-truncation reference (expensive but doable for one test case). The exact-phase result must lie between truncated-phase and reference, closer to reference. `R + T + A = 1` must hold to numerical precision.
3. **Improved agreement with full-MS reference tables** at glint and forward-peak geometries (Coulson-Dave-Sekera, Wauben-Hovenier).

**Validation regime for `ExactFirstOrderOnly + ExactSFIPhase`** (when implemented):
- Bit-for-bit comparison against VLIDORT FO output for matching geometries.

### 6.8 Lin-mode interaction

`ExactSFIPhase` changes how `J0±` is built; if those vectors feed into `Ṙ`/`Ṫ` derivatives via the lin path, the lin path needs an exact-phase source-construction method. **Real implementation task, flagged as Step 2.5.**

---

<a id="7-step-3"></a>
## 7. Step 3 — output products and RT modes/solvers as dispatched strategies

API consolidation. Builds on Steps 1, 1.5, 2.

### 7.1 The leak

`rt_run.jl:136-140` unconditionally allocates `R`, `T`, `R_SFI`, `T_SFI`, `ieR_SFI`, `ieT_SFI`, `hdr`, `bhr_uw`, `bhr_dw`. Return tuple at line 325 is conditional on `SFI`, but allocations aren't. 1004 lines across `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` for what is structurally one loop.

### 7.2 Three orthogonal axes

```julia
abstract type AbstractRTSolver end
struct FullMOM <: AbstractRTSolver end
struct ExactFirstOrderOnly <: AbstractRTSolver end

abstract type AbstractRTMode end
struct ForwardMode <: AbstractRTMode end
struct LinMode{PL} <: AbstractRTMode
    parameter_layout::PL
end

abstract type AbstractRTOutput end
struct DirectionalTOA <: AbstractRTOutput end
struct DirectionalBOA <: AbstractRTOutput end
struct HDRFOutput     <: AbstractRTOutput end
struct BHROutput      <: AbstractRTOutput end
struct MultiSensorOutput{...} <: AbstractRTOutput
    sensor_levels::Vector{...}
end
struct JacobianTOA <: AbstractRTOutput end
```

Plus `AbstractSFIPhase` from Step 2. **Four orthogonal dispatch axes** — solver, mode, sfi_phase, outputs. Each with its own family of strategy types.

### 7.3 Public API

```julia
run_rt(model;
       solver=FullMOM(),
       mode=ForwardMode(),
       sfi_phase=TruncatedSFIPhase(),
       outputs=(DirectionalTOA(), DirectionalBOA()),
       ss_corrections=default_ss_corrections(model),
       iBand=1)
```

**Defaults vs. legacy compatibility:**

- `run_rt`'s defaults are the *new clean API* — minimal, opt-in for additional outputs. Return type is `NamedTuple` with named output fields.
- `rt_run(model)` becomes a compat wrapper calling `run_rt` with the full legacy output set, assembling the legacy 7-tuple shape on return. Existing user code keeps working.
- `rt_run_lin` and `rt_run_multisensor` similarly become compat wrappers.

The two APIs have different defaults and different return shapes by design. Wrappers bridge.

### 7.4 What this collapses

| Current entry point | Becomes |
|---|---|
| `rt_run(model)` | `run_rt(model, outputs=(DirectionalTOA(), ...))` (with legacy outputs) |
| `rt_run` with `SFI=true` and RAMI | `run_rt(model, outputs=(DirectionalTOA(), HDRFOutput(), BHROutput()))` |
| `rt_run_multisensor.jl` (177 lines) | `run_rt(model, outputs=(MultiSensorOutput(levels),))` |
| `rt_run_lin.jl` (303 lines) | `run_rt(model, mode=LinMode(layout), outputs=(DirectionalTOA(), JacobianTOA()))` |

1004 lines collapse into one coordinator (~150 lines) plus per-output strategy methods.

### 7.5 Migration substeps

1. Define `AbstractRTOutput`, `AbstractRTMode`, `AbstractRTSolver`. No call sites changed.
2. Implement `DirectionalTOA`, `DirectionalBOA`. Both dispatch to existing internals. No behavior change.
3. Refactor current `rt_run` to take strategy kwargs. Introduce `run_rt` as the new entry point.
4. Implement `HDRFOutput`, `BHROutput`. Migrate RAMI workflows.
5. Implement `MultiSensorOutput`. `rt_run_multisensor` becomes wrapper.
6. Implement `LinMode` + `JacobianTOA`. `rt_run_lin` becomes wrapper. **Riskiest substep.**

Each substep is independently shippable; benchmarks at each boundary.

---

<a id="8-extension-cookbook"></a>
## 8. Extension cookbook (aspirational)

Recipes describing what the API will look like after Steps 1, 1.5, 2, 3 land.

### Recipe: adding a new BRDF surface kernel

```julia
struct MyBRDF{FT} <: AbstractSurfaceType
    # ... your parameters ...
end

function reflectance(brdf::MyBRDF, n_stokes, μᵢ, μⱼ, ϕ; ...)
    # ... your physics ...
end

n_fourier_moments(::MyBRDF, ctx) = 2 * ctx.Nstreams
has_specular_peak(::MyBRDF) = true
```

### Recipe: adding a new optical contributor (e.g. cloud water)

```julia
struct CloudWaterContributor{FT, M} <: AbstractOpticalContributor
    cloud_optics::M
    τ_cloud::Vector{FT}
    apply_delta_m::Bool
end

function optical_properties(c::CloudWaterContributor, ctx)
    # Build both truncated and exact forms (§5.5 pattern)
    # Return TruncatedAndExactScatteringOpticalProperties
end
```

### Recipe: adding a new SS correction type

```julia
struct MyCorrection <: AbstractSSCorrection
    # ...
end

function apply_correction!(c::MyCorrection, ctx)
    # ... mutate ctx.R_SFI in place ...
end
```

### Recipe: adding a new output product

```julia
struct PhotolysisRate <: AbstractRTOutput end
struct PhotolysisState
    # ...
end

allocate_output(::PhotolysisRate, model, nSpec) = PhotolysisState(...)
accumulate_output!(s::PhotolysisState, ::PhotolysisRate, ctx) = ...
finalize_output(s::PhotolysisState, ::PhotolysisRate, ctx) = ...
```

User code:

```julia
results = run_rt(model, outputs=(DirectionalTOA(), PhotolysisRate()))
```

### Recipe: adding a new RT solver

```julia
struct MyRTSolver <: AbstractRTSolver
    # ... configuration ...
end

function solve_rt!(::MyRTSolver, ctx)
    # ... your solver implementation ...
end
```

### Recipe: adding a new validation benchmark

(See §11.)

```toml
[source]
type = "published_full_ms"
citation = "Author et al. 2024, JQSRT N(M), pp-pp, Table X"
doi = "10.xxxx/yyyy"
regime = ["thin", "thick"]

[atmosphere]
# ...

[expected]
data_file = "my_new_case.h5"
```

---

<a id="9-type-gates"></a>
## 9. Remaining type-gates inventory (migration debt)

Beyond Step 1, 1.5, 2, 3, these places need future migration:

| Smell | Files | Replace with | Status |
|---|---|---|---|
| `if brdf isa CanopySurface` (4 occurrences) | `rt_run.jl:170, 175, 182, 245` | `requires_canopy_setup` trait | Step 3 |
| `if RS_type isa AbstractRamanType` style branches | `interaction_inelastic.jl`, `elemental_inelastic.jl` | Already partially dispatched — finish it | Future |
| `arch isa GPU` checks | scattered | Already mostly dispatched on `AbstractArchitecture`; audit holdouts | Future |
| `if SFI` boolean threading | every `rt_run*.jl`, `postprocessing_*.jl` | Goes away under Step 3 | Step 3 |
| Hardcoded `scatter = true` in lin/SS/multisensor | scattered | Query `LayerKind` via `scatter_predicate` | Step 1.5.7 |
| Polarization-type branching via `pol_type.n` | many places | Numerical property — likely fine, audit only | Audit only |
| `Stokes_IQU` implementation completeness | scattered | IQU code paths most likely to have rot | Future |
| `AbstractSourceType` resurrection | `types.jl:120-126` | `SolarSource`, `ThermalSource`, etc. | Future (thermal RT trigger) |

The forward-compat invariant from §2.9 keeps the `AbstractSourceType` resurrection mechanical when it happens.

---

<a id="10-future-architecture"></a>
## 10. Future architecture (deferred)

**C1 — Optical-properties algebra extension** (now Step 1.5, see §5). *In scope this round.*

**C2 — Layer-kind dispatch** (now §5.7). *In scope this round.*

**C3 — Validation-as-data** (now §11). *In scope this round as enabling work for Step 2.*

**C4 — Layer-kernel interface stabilization.** Doubling/elemental/interaction kernels coupled to specific layer-storage types. Future-proof: `AbstractLayerStorage` interface. Belongs in linearization redesign track.

**C5 — Physics audit traits.** A `physics_summary(model)` returning structured description of every active assumption. Useful for reproducibility, sanity checks. Cheap given dispatch architecture.

**C6 — Mie greek-coef interpolation.** Aerosol-specific. Compute greek coefs at sparse "anchor" wavelengths, interpolate in size-parameter space. The contributor pattern handles this naturally.

**C7 — Versioned interface contracts.** External BRDF/aerosol/surface implementers need stable contracts across vSmartMOM minor versions. Documented interface methods with `Required`/`Optional` annotations and a `check_interface(MyBRDF)` test helper.

**C8 — Thermal RT.** `ThermalSource` adds Planck per-layer source dispatching through `source_term(::ThermalSource, ctx)`. With `NoScattering` layer-kind already in place from Step 1.5 and the `ctx.source` invariant from §2.9, this slots in mechanically.

**C9 — `ExactFirstOrderOnly` solver implementation.** The dispatch infrastructure (`AbstractRTSolver`) lands in Step 2; the actual `ExactFirstOrderOnly` body lands when there's demand. Useful for fast first-order computations, debugging, FO-comparisons against VLIDORT.

---

<a id="11-validation"></a>
## 11. Validation infrastructure

### 11.1 The framing distinction

vSmartMOM is a full multiple-scattering code. **Reference tables for validation must be from full-MS codes.**

2OS approximations (Natraj-Spurr 2007) are *not appropriate* — comparing full-MS output against 2OS would penalize vSmartMOM for capturing higher-order scattering correctly.

| Reference type | What it validates | Where it applies |
|---|---|---|
| Full-MS reference | Total radiance | Always |
| Single-scatter-only reference (e.g. VLIDORT FO output) | `ExactFirstOrderOnly` solver mode | Future (when ExactFirstOrderOnly is implemented) |
| 2OS reference | (nothing in vSmartMOM directly) | Don't use |
| Analytic Lambertian SS | Direct-beam in trivial case | Step 2 sanity |

### 11.2 The metadata-plus-data-file pattern

Each reference benchmark is a TOML metadata file plus an HDF5/CSV data file. One harness iterates over all of them.

```toml
# test/reference_data/published/coulson_dave_sekera_table_iv.metadata.toml
[source]
type = "published_full_ms"
citation = "Coulson, Dave & Sekera 1960, Tables Related to Radiation Emerging from a Planetary Atmosphere with Rayleigh Scattering, U. California Press, Table IV"
regime = ["thin", "thick"]

[atmosphere]
profile = "rayleigh_only"
τ_rayleigh = 0.5
n_layers = 1
contributors = ["rayleigh"]

[surface]
type = "lambertian"
albedo = 0.0

[geometry]
sza_deg = 60.0
vza_deg = [0, 30, 60]
vaz_deg = [0, 90, 180]

[expected]
data_file = "coulson_dave_sekera_table_iv.h5"
fields = ["I", "Q", "U"]
shape = [3, 3, 3]
tolerance_relative = 1e-4
```

### 11.3 Type-dispatchable validation

```julia
abstract type AbstractValidationRegime end
struct FullMS <: AbstractValidationRegime end
struct SingleScatterOnly <: AbstractValidationRegime end
struct TwoOrdersOnly <: AbstractValidationRegime end

validate_against_reference(::FullMS, vsmartmom_output, expected) = ...
validate_against_reference(::SingleScatterOnly, vsmartmom_output, expected) = ...
function validate_against_reference(::TwoOrdersOnly, vsmartmom_output, expected)
    @warn "Skipping 2OS reference; vSmartMOM is full-MS"
end
```

Cross-regime comparisons are caught by the type system.

### 11.4 Reference set for v2.0

**Published, committed to the repo:**

- **Coulson, Dave & Sekera 1960** — Rayleigh-only canonical. Full-MS. τ ∈ [0.05, 1.0]. Stokes I, Q.
- **Garcia & Siewert 1989** — Benchmark suite with varying complexity. Full-MS.
- **Spurr 2002** VLIDORT validation appendix — Full-MS cross-validation against multiple codes.
- **Wauben & Hovenier 1992** — Vector RT with aerosols. Full-MS.

**Locally-generated (gitignored, dev-machine only):**

- VLIDORT distribution `results/` directory output.
- VLIDORT FO output for Step 2 direct-beam validation (when ExactFirstOrderOnly lands).
- Cox-Munk specific configurations not in the published set.
- RAMI HOM26/HOM27 canopy hotspot cases.

```
test/reference_data/
  published/                                ← committed to repo
    coulson_dave_sekera_table_iv.{metadata.toml, h5}
    garcia_siewert_1989_case1.{metadata.toml, h5}
    spurr_2002_appendix_a.{metadata.toml, h5}
    wauben_hovenier_1992_table2.{metadata.toml, h5}
  gitignored/                               ← .gitignored
    cfranken_local_cox_munk_w15_fo.{metadata.toml, h5}
    cfranken_local_canopy_rami_hom26.{metadata.toml, h5}
```

CI runs only `published/`. Dev validation runs both.

### 11.5 Conservative engineering policy on VLIDORT-derived data

- Don't commit VLIDORT source code (clear).
- Published numerical tables from peer-reviewed papers are the safest committed baseline. Citation handles provenance.
- Locally-generated VLIDORT outputs go in the gitignored directory regardless of any "facts vs. source" reasoning, until licensing for *those specific outputs* is confirmed.
- The VLIDORT distribution `results/` directory contains numerical outputs but they're shipped under the package's GPL terms; treat as gitignored until provenance is sorted.

This is more conservative than v0.4's framing. The win is small (most validation references are published tables anyway) and the risk is real.

### 11.6 What enables this in Step 1

- Step 1A.0.6 — write the metadata format spec, the harness, and one test case (Coulson-Dave-Sekera) as proof-of-concept.
- Step 1B.4 onward — every Step 1 substep that changes numerics gets a validation against the published set.
- Step 2 onward — validates `FullMOM + ExactSFIPhase` against published full-MS references at glint/forward-peak, plus the energy-conservation sanity check.

---

<a id="12-execution-plan"></a>
## 12. Execution plan

### Step 1 — Fourier work, SS-correction plumbing, validation harness, small algebra fixes

| # | What | Risk | Validation |
|---|---|---|---|
| 1A.0 | Rename `AbstractFourierDecompositionType` → `AbstractMieDecompositionAlgorithm`. Add deprecation alias. | ~50–100 lines including tests, docs, IO imports. Touches exported API. | All existing tests pass. |
| 1A.0.5 | Deprecate `GaussQuadFullSphere`. `Base.depwarn`. Update `DefaultParameters.yaml` default to `GaussLegQuad`. Update fixtures. Update 4 docs files. | **YAML default change** — behavior shift. | Audit benchmarks for `FullSphere` use. |
| 1A.0.6 | Validation harness foundation: metadata-plus-data-file format, type-dispatchable validation regime, one published reference (Coulson-Dave-Sekera) as proof-of-concept. | New infrastructure. | Harness round-trips on the proof-of-concept reference. |
| 1A.1 | Introduce `FourierDecomposition` module skeleton. | Plumbing. | Module loads. |
| 1A.2 | Add `RTContext` struct (§1.3, parameterized concretely). | New abstraction. | Trivial unit tests. |
| 1A.3 | Add `Nstreams` field to `QuadPoints`. | One new field, two value-producers. | All quadrature setups produce both values. |
| 1A.4 | Change `SolverConfig.max_m` to `Union{Nothing, Int}`. YAML accepts `null` or omission. | **Behavior change for old configs** with explicit integers. | Existing tests pass. Document in changelog. |
| 1A.5 | Add `n_fourier_moments` trait + methods. | New code. | Unit tests for each surface type. |
| 1A.6 | Wire `n_fourier_moments(model, iBand)` into consumers. | **Behavior change risk.** Cox-Munk bands may compute more moments. | 6SV1 / Natraj reference benchmarks; published harness. |
| 1A.7 | (Deferred) Promote `nQuad_BRDF` to per-surface field. | Trivial. | Convergence sweeps. |
| 1B.1 | Refactor BRDF Fourier integration into shared helper. | Pure refactor. | Bit-for-bit equivalence. |
| 1B.2 | Add `AbstractSSCorrection` hierarchy. Migrate Cox-Munk apply. **Read `n_moments` from ctx.** | Pure relocation. | Bit-for-bit equivalence. |
| 1B.3 | Wire `apply_correction!` into `rt_run.jl` via `foreach`. Add `default_ss_corrections`. | Pure call-site swap. | Existing benchmarks unchanged. |
| 1B.4 | Add lin-context apply method updating `R_SFI` and `Ṙ_SFI` consistently. | New physics in lin path. | Validate Jacobians vs. finite-difference at glint. |
| 1B.5 | Add multisensor-context apply method. | New physics in multisensor path. | Validate against multisensor reference cases. |
| 1C.1 | Add `hotspot::HS` parameter to `CanopySurface`, parameterized, default `NoHotSpot`. Plumb through `model_from_parameters` and YAML. | Pure plumbing; no-op runtime. | All existing tests pass bit-for-bit. |
| 1D.1 | Document optical-properties algebra (§5.2). Add `+(absorption, absorption)` overload. **Fix `y.G` bug at types.jl:1137.** | Mostly docs + 1 trivial overload + 1 bug fix. | All existing tests pass. |

### Step 1.5 — Optical-properties algebra and layer-kind dispatch

Peer in scope and risk to Step 2. Validates with its own equivalence test boundary.

| # | What | Risk | Validation |
|---|---|---|---|
| 1.5.1 | Add dimensioned `zero_optical_properties` constructor (§5.3). Verify dimensions (`τ` is `nSpec`, `n_pq` uses augmented `Nquad`). | New constructor. | Unit tests for zero-element behavior. |
| 1.5.2 | Add `AbstractOpticalContributor` + `optical_properties` trait. | New abstraction. | Unit tests for each contributor type. |
| 1.5.3 | Migrate `constructCoreOpticalProperties` to contributor pattern. **Scope: non-directional algebra only.** Canopy stays separate. | **Touches the core pipeline.** | Bit-for-bit equivalence. |
| 1.5.4 | Add `TruncatedAndExactScatteringOpticalProperties` dual-form type (§5.5). Migrate contributors to produce dual forms. | Type-system change. | Bit-for-bit equivalence (consumers read `.truncated`). |
| 1.5.5 | Promote storage form to type parameter (`Compact` / `Expanded`). Refactor `+`/`*` overloads. | Algebra restructure (4 overloads instead of 1). | Bit-for-bit equivalence. |
| 1.5.6 | Add `AbstractLayerKind` with `Scattering`/`NoScattering`. Use existing band-level max-based predicate. Implement `elemental!(::NoScattering, ...)` and `doubling!(::NoScattering, ...)`. | New kernel methods. | Bit-for-bit equivalence on default scattering atmospheres; new pure-absorption test case. |
| 1.5.7 | Audit lin/SS/multisensor paths for hardcoded `scatter = true`; migrate to `scatter_predicate(layer_kind(...))`. | Touches multiple paths. | Bit-for-bit equivalence. |
| 1.5.8 | `_expand_layer_rayleigh!` migration: choose Option A (trait) or Option B (push to use site) and implement. | Plumbing change. | Bit-for-bit equivalence on Raman runs. |

### Step 2 — `AbstractSFIPhase` source dispatch

| # | What | Risk | Validation |
|---|---|---|---|
| 2.1 | Define `AbstractSFIPhase` types. Dispatch `build_sfi_source!` on phase. | New abstraction. | Unit tests; `TruncatedSFIPhase` is bit-for-bit current. |
| 2.2 | Implement `ExactSFIPhase` source construction reading dual-form properties (§5.5). | Real physics work. | Convergence to truncated when `f → 0`; bit-for-bit when `apply_delta_m=false`. |
| 2.3 | Add `AbstractRTSolver` types. Dispatch infrastructure for `FullMOM`. | New abstraction. | Unit tests. |
| 2.4 | Wire canopy hotspot via `apply_hotspot!` dispatch on `surface.hotspot`. | New physics for canopy. | RAMI HOM26/HOM27. |
| 2.5 | Lin-context `build_sfi_source!(::ExactSFIPhase, ::LinContext)`. | New physics in lin path. | Validate Jacobians against finite-difference. |
| 2.6 | Energy-budget validation: forward-peak aerosol case, R+T+A = 1, agreement with full-MS reference. | Physics validation. | Energy conservation; full-MS table comparison. |
| 2.7 | Glint/forward-peak validation against published full-MS tables. | Physics validation. | Coulson-Dave-Sekera, Wauben-Hovenier. |

### Step 3 — Output products and RT modes/solvers

Per §7.5. Six substeps. Lin migration last.

### Future steps (post-v2.0)

- Step 4: `ExactFirstOrderOnly` solver implementation.
- Step 5: Thermal RT (C8).
- Other items from §10 as triggered.

---

<a id="13-open-questions"></a>
## 13. Open questions

In rough order of "blocks the next commit" → "nice to settle eventually":

1. **Q1 (preamble)** — dispatch rule as governing principle, applied at every level.
2. **Q2 (preamble)** — `ExactSFIPhase` corrects more paths than VLIDORT FO; validation regime. **Blocks Step 2.6 / 2.7 validation plan.**
3. **Q3 (preamble)** — `AerosolForwardPeakCorrection` obsoleted by `ExactSFIPhase`; v0.5 drops it. Residual post-processing component intended? **Blocks Step 1B.2 type definitions.**
4. **Q4 (preamble)** — `NoScattering` layer kind is band-level (max-based predicate matches existing code). **Blocks Step 1.5.6 design.**
5. **Q5 (preamble)** — canopy hotspot regression risk. Are there benchmarks tuned to no-hotspot result? **Blocks Step 1C.1 changelog and Step 2.4.**
6. **Q6 (preamble)** — full-MS reference set for v2.0. Additional sources?
7. **Q7 (preamble)** — `n_fourier_moments(::CoxMunk) = 2 · Nstreams` derivation. Empirical sanity check ordering.
8. **Q8 (preamble)** — `GaussQuadFullSphere` retirement, including YAML default change.
9. **`l_trunc` semantics** — empirical check via Q7 sanity test.
10. **`fScattRayleigh` extraction** — Option A (trait) vs. Option B (push to use site). Step 1.5.8.
11. **`ExternalSSCorrection` callback signature** — mutation-via-context vs. return-value API.
12. **`RTContext` lifecycle** — exact field set and m-mutation contract.
13. **Naming harmonization** — `AbstractTruncationType`/`NoTruncation` vs. `Abstract*Correction`/`No*Correction`.
14. **`m_max(::CanopySurface)` provisional vs. proper** — placeholder is `2·Nstreams`. LAD-driven derivation in CanopyOptics.
15. **`run_rt` API surface area** — four orthogonal axes (solver, mode, sfi_phase, outputs) is a lot of knobs. Are all combinations worth exposing, or do some stay internal?

---

<a id="14-references"></a>
## 14. References

- **VLIDORT 2.8.3 source:**
  - `vsup/vbrdf/vbrdf_sup_routines.f90` — `VBRDF_FOURIER`
  - `vsup/vbrdf/vbrdf_sup_masters.f90:2305` — `NMOMENTS = 2·NSTREAMS - 1` invariant
  - `vlidort_focode/VFO_Master.f90` — first-order (exact SS) subsystem
  - `vlidort_def/vlidort_inputs_def.f90:444-481` — `DO_FOCORR*` flag hierarchy
- **Nakajima & Tanaka 1988**, *Algorithms for radiative intensity calculations in moderately thick atmospheres using a truncation approximation*, JQSRT 40(1), 51-69. Original δ-M / TMS scheme.
- **Spurr 2002**, *VLIDORT: A linearized pseudo-spherical vector discrete ordinate radiative transfer code*. Original VLIDORT linearization.
- **Coulson, Dave & Sekera 1960**, *Tables Related to Radiation Emerging from a Planetary Atmosphere with Rayleigh Scattering*, U. California Press.
- **Garcia & Siewert 1989**, *Benchmark results in radiative transfer*, Transp. Theor. Stat. Phys.
- **Wauben & Hovenier 1992**, *Polarized radiation of an atmosphere containing randomly-oriented spheroids*, A&A.
- **Kuusk 1991**, *The hot spot effect in plant canopy reflectance*. Original Kuusk hotspot model.
- **CanopyOptics.jl** `src/canopy_structure/hotspot.jl`.
- **RAMI** (Radiation transfer Model Intercomparison) test scenes — HOM26 / HOM27.
- **Codex review of v0.1, v0.2, v0.3** and **GPT review of v0.3, v0.4** and **SFI-phase note**, integrated throughout this document.
