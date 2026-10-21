# vSmartMOM.jl — architecture design and Fourier work

**Document version:** v0.4
**Status:** working draft for Sanghavi review.
**Branch target:** `sanghavi-unified` (with `unified-vsmartmom` as the architectural merge base, physics layered in from `sanghavi`).
**Companion:** `LINEARIZATION_BUGS.md`, the linearization redesign track.

---

## Preamble — for Sanghavi

This document is a design proposal for vSmartMOM v2.0.0. It has now been through three rounds of review (Codex against current `sanghavi-unified`; GPT against v0.3) and reaches you with the architectural skeleton settled. Your sign-off on the **direction** matters more than implementation detail, because:

- You own the theory. The architectural choices are decisions about how the *code* expresses the *physics* — if the physics framing is wrong, the architecture is too.
- Several choices are not technically forced. They reflect judgment about what kind of codebase vSmartMOM should be in two years. Your view of where the science is going should shape that.

### What this document is trying to do

The proposal has two intertwined parts:

1. **Establish a single architectural rule** — multiple dispatch at physical and workflow boundaries; `if/else` only for scalar numerical cases — and apply it consistently *at every level of the pipeline*. Surface dispatch (BRDF), atmospheric dispatch (optical contributors), layer dispatch (`Scattering` vs. `NoScattering`), storage dispatch (`Compact` vs. `Expanded` Z), output dispatch (which products to record), mode dispatch (forward vs. linearized). Every architectural seam exposed as a type-dispatchable boundary.

2. **Apply this rule to fix five real, currently-broken or currently-missing pieces of physics**:
   - Per-component Fourier-moment count (currently surface-blind — Cox-Munk and Lambertian both get the same `max_m` from aerosol Greek coefficients)
   - Single-scattering correction (currently a hack in `rt_run.jl:311` missing from linearized and multisensor entry points)
   - Canopy hotspot (CanopyOptics ships it; vSmartMOM doesn't wire it — every current canopy result silently uses independent Beer-Lambert)
   - Optical-properties algebra (the `+`/`*` operators on `CoreScatteringOpticalProperties` are elegant but undocumented and contain a real bug at `types.jl:1137`; thermal/MW work needs a proper layer-kind dispatch)
   - Output-products and RT-mode topology (1004 lines of forked code across `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` for what is structurally one loop)

The proposal lands in **three steps**:

- **Step 1** (this round): Fourier work, SS corrections, optical-properties algebra, layer-kind dispatch, `GaussQuadFullSphere` retirement, validation harness. Mostly refactor and small physics fixes.
- **Step 2** (next round): Direct-beam separation subsystem (analog of VLIDORT's FO). Real solver work; validates against VLIDORT FO. Enables canopy hotspot wiring.
- **Step 3** (final round): Output products and RT modes as dispatched strategies. The `run_rt(model; mode, outputs)` API. Deduplicates the three forked `rt_run*.jl` files.

### Six questions for you up front

These are the choices that affect everything else. Each is one paragraph; pointers to the technical detail are references, not required reading for the philosophical pass.

**Q1 — Dispatch architecture as governing principle, applied at every level.** The proposal makes "dispatch on types at physical/workflow boundaries; `if/else` only for numerical cases" the *rule* the codebase commits to. It applies to surface kernels, optical contributors, layer kinds, storage forms, output products, and run modes — not just the outer levels. New BRDFs, new corrections, new outputs, new aerosol contributors, new layer kinds (e.g. `NoScattering`) all follow the same pattern: small concrete strategy struct, single-method dispatch. Are you on board with the rule as a hard commitment going forward, or is there a class of physics where you think the rule is wrong-headed and the freedom to branch is worth keeping? *Detail: §1.*

**Q2 — Canopy hotspot is silently missing today; the fix has cost.** vSmartMOM's `elemental_canopy!` uses independent Beer-Lambert (`exp(-(k_s + k_o) L)`) where it should use joint gap probability (`exp(-(k_s + k_o) L) · C_hs(α, L)`). CanopyOptics ships `KuuskHotSpot` with ForwardDiff hooks but vSmartMOM doesn't call `joint_gap_probability` anywhere. *Any current canopy benchmark of yours that compares against RAMI/PROSAIL hotspot-aware references is either wrong or was tuned to match the no-hotspot result.* We'd like to know if any such benchmarks exist before Step 2 lands the fix — those would need to be regenerated, not silently corrected. *Detail: §4.*

**Q3 — `GaussQuadFullSphere` retirement.** It constructs `gausslegendre(2·Nquad)` over `[-1, 1]` and discards the lower hemisphere; the retained nodes/weights are *not* an optimal Gauss-Legendre quadrature for `[0, 1]`. For RT problems where upper-hemisphere integration is what matters (essentially all plane-parallel RT), it's strictly worse than `GaussLegQuad`. Currently it's even the *default* in `DefaultParameters.yaml:12`. Proposal: deprecate in v2.0 (warning, still works), make `GaussLegQuad` the default, remove in v2.1. Do you concur, or is there a use case to preserve? *Detail: §2.8.*

**Q4 — Optical-properties algebra: `NoScattering` layer kind unifies thermal/MW.** The `+`/`*` operators on `CoreScatteringOpticalProperties` are mathematically clean but undocumented, contain a bug at `types.jl:1137` (the scalar `*` overload references `y.G` on a type that doesn't have that field), and don't compose with thermal/MW physics (no Rayleigh in IR/MW). Proposal: (a) document the algebra formally; (b) introduce `AbstractLayerKind` with `Scattering` and `NoScattering` subtypes, dispatched at the kernel level — `NoScattering` layers go through closed-form Beer-Lambert + Planck source, no doubling needed; (c) promote storage form (compact vs. expanded Z) to a type parameter so kernels specialize; (d) introduce `AbstractOpticalContributor` so atmospheric components are responsible for their own optical-properties contribution via the existing `+` algebra. **One pipeline, dispatched kernels handle thermal and scattering atmospheres.** Does this match how you think about RT, or is there a structural reason to prefer a separate thermal solver? *Detail: §5.*

**Q5 — `n_fourier_moments(::CoxMunk) = 2 · Nstreams` derivation.** Translating VLIDORT's `NMOMENTS = 2·NSTREAMS - 1` (Fortran inclusive loop, `2·NSTREAMS` iterations) into vSmartMOM's count convention gives `2 · Nstreams` where `Nstreams` is the *base* polar discrete-ordinate count (not `model.quad_points.Nquad`, which is augmented with VZA/SZA nodes). For `l_trunc = 60, Nstreams = 30`, the count is 60. Does this match your physical intuition? The empirical sanity check is to run a Cox-Munk + Rayleigh test with `n_moments` set to 60 vs. 30 vs. 100 and check convergence; we'd value your view on the expected reference numerics before running that test. *Detail: §2.*

**Q6 — Validation against full-MS reference tables only.** vSmartMOM is a full multiple-scattering code (matrix operator method, doubling-adding to convergence). Reference tables for validation must be from full-MS codes; 2OS approximations like Natraj-Spurr 2007 are *not appropriate* even though they share VLIDORT lineage — comparing full-MS output against 2OS would penalize vSmartMOM for capturing higher-order scattering correctly. Proposed reference set: Coulson-Dave-Sekera 1960 (Rayleigh canonical), Garcia-Siewert 1989 (benchmark suite), Spurr 2002 VLIDORT validation appendix, Wauben-Hovenier 1992 (vector with aerosols), and the VLIDORT distribution `results/` directory directly. The metadata-plus-data-file harness (§11) treats validation regime as a type-dispatchable property so cross-regime comparisons are caught. Does this match the reference set you've been using, or do you have additional sources we should add? *Detail: §11.*

### What this document is *not*

- Not a fait accompli. The open-questions list (§12) flags everything we know is open.
- Not yet implemented. Current `sanghavi-unified` has none of these changes.
- Not the linearization redesign. That's a separate track. This proposal is upstream of it — landing the dispatch architecture cleanly makes the linearization redesign tractable.
- Not exhaustive of every smell. Some pieces (e.g. `AbstractSourceType` for thermal/MW source dispatch) are flagged as future work and explicitly *not* refactored in this round, but Step 1's design is structured so resurrecting them later is mechanical, not architectural.

The rest of the document treats the technical detail. Sections 1–5 are the architecture and Step 1 work items. Section 6 is Step 2 (direct-beam subsystem). Section 7 is Step 3 (output products and modes). Sections 8–10 are extension recipes, migration debt, and deferred future work. Section 11 is validation infrastructure. Section 12 is the execution plan. Section 13 is open questions.

---

## Iteration log

- **v0.4** — major revision after GPT review of v0.3. Folded throughout: `max_m` `Union{Nothing,Int}` semantics, correct dimensions in `zero_optical_properties` (`τ` is `nSpec`, Z uses augmented `Nquad`), Step 1D scoped to non-directional algebra (with bug fix at `types.jl:1137`), source-type rename dropped from Step 1 (forward-compat invariant kept), `run_rt` defaults vs. legacy tuple shape clarified, `ss_corrections` as kwarg not `SolverConfig` field, `GaussQuadFullSphere` audit scope expanded. Added: `AbstractLayerKind` with `Scattering`/`NoScattering` dispatch (one solver, dispatched kernels — replaces v0.3's "thermal/MW gets its own solver entry point" framing). Storage form (`Compact`/`Expanded`) promoted to type parameter. Three-step plan: Step 2 (direct-beam subsystem) split out from Step 3 (output products/modes) since direct-beam is real solver work with its own validation boundary. Validation infrastructure section added (§11) with full-MS reference set; 2OS references explicitly excluded. Mie greek-coef interpolation as new C6 future work.
- **v0.3** — Sanghavi-facing preamble; Codex v0.2 review folded; Step 1D (optical-properties algebra) added.
- **v0.2** — promoted dispatch rule to governing principle; Codex v0.1 review folded.
- **v0.1** — initial draft.

---

## Table of contents

1. [The dispatch rule (governing principle)](#1-dispatch-rule)
2. [Step 1A — per-component `n_fourier_moments`](#2-step-1a)
3. [Step 1B — `AbstractSSCorrection` strategies](#3-step-1b)
4. [Step 1C — `AbstractHotSpot` field on `CanopySurface`](#4-step-1c)
5. [Step 1D — optical-properties algebra and layer-kind dispatch](#5-step-1d)
6. [Step 2 — direct-beam separation subsystem](#6-step-2)
7. [Step 3 — output products and RT modes as dispatched strategies](#7-step-3)
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
- **Layer-kind classification at construction** (`props.ϖ ≈ 0 ? NoScattering() : Scattering()` — one comparison per layer at construction time, then dispatch from there)

The rule eliminates *type and physics gates*, not all branches. Litmus test: if the branch's truth is determined at model-construction time by which types are in play, dispatch on types. If determined at evaluation time by numerical values or user toggles, branch.

**Applied at every level.** The rule is not just for outer dispatch (surface kernels, output products). It applies *all the way down* — to layer kinds, storage forms, kernel specialization. The compounding payoff: the compiler can specialize through the entire stack when types are concrete, and the human reader can find each physics case by following a single dispatch chain.

### 1.2 The strategy struct shape

Strategy objects are **small concrete structs held in tuples, dispatched via single-method functions.**

**(a) Concrete tuples, not abstract vectors.** A `Vector{AbstractSSCorrection}` forces runtime dispatch in inner loops. A `Tuple{BRDFSpecularCorrection{CoxMunkSurface{Float64}}}` lets the compiler specialize.

**Honest about type stability:** when a strategy set is built conditionally and stored in a long-lived struct like `SolverConfig`, end-to-end compile-time dispatch requires either parameterizing the storage struct on the strategy-tuple type (invasive) or accepting runtime dispatch through a function-barrier when unpacked. **Pragmatic answer: accept runtime dispatch at outer levels (per-band, per-run) where the strategy work itself is much more expensive than dispatch overhead. Compile-time dispatch is achieved *inside* each strategy method via concrete type parameters, not at the outer driver level.** Layer-kind and storage-form dispatch *do* propagate to compile-time because they're parameterized on `Layer` and `CoreScatteringOpticalProperties` respectively.

**(b) Single-method functions, not switch tables.** Each strategy owns its behavior:

```julia
apply_correction!(c::BRDFSpecularCorrection, ctx)
apply_correction!(c::AerosolForwardPeakCorrection, ctx)
apply_correction!(::NoSSCorrection, ctx) = nothing
```

**(c) Naming consistency across families.**

| Family | Trait | Apply |
|---|---|---|
| Fourier bounds | `n_fourier_moments(component, ctx)` | (count only) |
| SS corrections | `needs_correction(c, ctx)` | `apply_correction!(c, ctx)` |
| Optical contributors | `optical_properties(c, ctx)` | composes via `+` algebra |
| Layer kinds | (read off from `ϖ`) | `elemental!(::Scattering, ...)` / `elemental!(::NoScattering, ...)` |
| Storage forms | (type parameter) | `expandOpticalProperties(::Compact)` returns `Expanded` |
| Output products | (presence in tuple is signal) | `accumulate_output!(state, output, ctx)` |
| Hotspot models | (parameter on `CanopySurface`) | `joint_gap_probability(model, …)` (in CanopyOptics) |

### 1.3 The `RTContext` object

Strategy methods need a consistent argument convention. `RTContext` is constructed once before the Fourier loop, mutated as `m` advances.

**Lifecycle:**

- Constructed once before the Fourier loop. Carries everything that doesn't depend on `m`.
- Field `m` is mutated as the loop advances. Strategy methods read it, don't mutate it.
- `composite_layer` accumulates **vertically within one Fourier moment** (across `iz` layers), reset at `iz == 1` on the next `m`.
- **Output arrays accumulate across `m`** — these belong to `AbstractRTOutput` strategies (Step 3).
- `τ_sum_all` is m-independent.

```julia
mutable struct RTContext{FT, S<:AbstractSourceType, ...}
    model::RTModel
    iBand::Int
    m::Int
    quad_points::QuadPoints
    composite_layer            # vertically accumulated within m, reset per-m
    added_layer_surface
    τ_sum_all::Matrix{FT}      # cumulative optical depth, m-independent
    geometry::Geometry         # observation: vza, vaz
    source::S                  # solar today; thermal/lunar later via §2.9 invariant
    pol_type::AbstractPolarizationType
    arch::AbstractArchitecture
    workspace                  # InteractionWorkspace, allocated once
    Nstreams::Int              # base polar streams (NOT augmented Nquad)
    Nquad::Int                 # augmented quadrature including VZA/SZA
    nSpec::Int                 # spectral grid size — needed for τ dimensions
    n_layers::Int
    FT::Type{FT}
end
```

Exact field set is open — see §13 Q11.

### 1.4 What the rule unlocks

- **Step 1 (§§2–5)** sketches the idiom on small contained pieces. Five work items, each <250 lines.
- **Step 2 (§6)** applies the rule to the direct-beam separation subsystem (analog of VLIDORT's FO).
- **Step 3 (§7)** applies the rule to output products and RT modes. Collapses 1004 lines of forked `rt_run*.jl` into one coordinator.
- **The linearization redesign** uses the same rule. Lin Jacobians are an `AbstractRTMode`, orthogonal to `AbstractRTOutput`.
- **Future corrections** (thermal SS, Raman SS, multiple-scatter δ-M, polarization SS for thermal) all fit the same `Abstract*Correction` shape.

The rule has compounding returns. Step 1 is the first application; everything downstream gets cheaper.

### 1.5 Branches that stay (the "fine" list)

- `m == 0` Fourier weight normalization. Numerical, stays.
- Cache initialization (`canopy._cache === nothing`). State, not type/physics dispatch.
- `props.ϖ ≈ 0` check at layer construction to classify `Scattering` vs. `NoScattering`. One comparison per layer at construction; dispatch from there.
- Numerical NaN/denorm guards in kernels.

The rule eliminates type and physics gates. Not all branches.

---

<a id="2-step-1a"></a>
## 2. Step 1A — per-component `n_fourier_moments`

### 2.1 Convention: count, not order

Existing code uses `max_m` as a count: loops are `for m = 0:max_m - 1` in `rt_run.jl:208`, etc. The trait must match.

**Decision:** the trait is `n_fourier_moments` and returns a count. Loop convention stays `for m = 0:n - 1`.

| Component | `n_fourier_moments` | Why |
|---|---|---|
| Lambertian (any flavor) | 1 | Only `m = 0`; no azimuthal structure |
| `RayleighScattering` | 3 | `m = 0, 1, 2` cover Rayleigh greek coefs |
| `AerosolOptics` (greek_coefs.β of length L+1) | min(L+1, 2·Nstreams) | Greek-coef tail capped at solver resolution |
| `CoxMunkSurface`, `rpvSurfaceScalar`, `RossLiSurfaceScalar` | `2 · Nstreams` | VLIDORT cap (§2.3) |
| `CanopySurface` | `2 · Nstreams` (provisional) | Proper LAD-driven derivation deferred |

### 2.2 `Nstreams` is *not* `Nquad`

`model.quad_points.Nquad` is *augmented* — `rt_set_streams.jl:40` redefines `Nquad = length(qp_μ)` after appending VZA/SZA nodes. If `n_fourier_moments(::CoxMunk, ctx) = 2 * ctx.Nquad`, the count silently varies with the user's requested viewing angles. Different VZAs give different numbers of Fourier modes — not what VLIDORT's invariant says, not physically meaningful.

**Fix:** introduce a separate `Nstreams` capturing the polar discrete-ordinate count *before* VZA/SZA augmentation:

```julia
mutable struct QuadPoints{FT}
    # ... existing fields ...
    Nstreams::Int   # base solver streams, never augmented
    Nquad::Int      # may equal Nstreams (Gauss hemisphere) or include VZA/SZA (Radau)
end
```

Computed at quadrature setup time as `Nstreams = (l_trunc + 1) ÷ 2` *before* any VZA/SZA augmentation. Existing `Nquad` field stays for kernel sizing.

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

Empirical sanity check planned (Q5).

### 2.4 The grids (kept separate)

| Quantity | Role | Source | Notes |
|---|---|---|---|
| `Nstreams` | Solver polar discrete-ordinate count, base | `rt_set_streams.jl` (new field) | Used for VLIDORT cap |
| `Nquad` | Augmented quadrature including VZA/SZA | `rt_set_streams.jl:40` | Used for kernel sizing, Z-matrix dimensions |
| `nQuad_BRDF` | Azimuth integration grid for each Fourier coefficient | Hard-coded `100` in surface files | Quadrature accuracy, deferred to per-surface field |
| `n_fourier_moments` | Count of orders the loop runs over | New trait | The trait being added |

### 2.5 Consumers reading the wrong value, and `max_m` semantics

**Two-part fix** (Codex Issue 3 + GPT Issue 1):

**Part 1 — the consumer read.** Even with `max_m_bands` populated correctly in `model_from_parameters`, the loop driver in `rt_run.jl:97` calls `get_max_m(model)` which returns `model.solver.max_m` (the *scalar*), not the per-band vector.

Wiring change: replace `max_m = get_max_m(model)` with `n_moments = n_fourier_moments(model, iBand)` in `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`.

**Part 2 — the `max_m` override semantics.** Currently the YAML schema requires `radiative_transfer.max_m` as an integer (`Parameters.jl:493`). There's no way to distinguish "user explicitly capped" from "default value." Old configs with `max_m: 30` would still silently cap Cox-Munk to 30 even though the trait says 60.

Fix: change `max_m` to `Union{Nothing, Int}`:

```julia
struct SolverConfig{...}
    # ...
    max_m::Union{Nothing, Int}   # nothing = "use trait"; Int = "user override cap"
end

function n_fourier_moments(model::RTModel, iBand::Int)
    ctx = (; Nstreams = model.quad_points.Nstreams)
    n = n_fourier_moments(model.optics.rayleigh, ctx)
    for a in model.optics.aerosols.aerosol_optics[iBand]
        n = max(n, n_fourier_moments(a, ctx))
    end
    n = max(n, n_fourier_moments(model.surfaces[iBand], ctx))
    n = min(n, 2 * ctx.Nstreams)        # solver resolution cap
    return isnothing(model.solver.max_m) ? n : min(n, model.solver.max_m)
end
```

YAML default becomes `max_m: null` (or omits the field entirely; loader treats absence as `nothing`). Existing configs with explicit integers continue to work as user overrides. **Old configs with the historical default integer must be regenerated** — flag this in the changelog.

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

n_fourier_moments(::CanopySurface,            ctx) = 2 * ctx.Nstreams   # provisional

n_fourier_moments(a::AerosolOptics,           ctx) =
    min(length(a.greek_coefs.β), 2 * ctx.Nstreams)
```

### 2.7 Knock-on simplifications

- **`apply_correction!` reads `n_fourier_moments` from `ctx`**, never recomputes — guarantees the correction subtracts the same Fourier sum that was added.
- **Lambertian-only bands skip `m > 0`**. Loop runs once.
- **`vSmartMOM_Parameters.max_m` becomes a user-override knob**, not the primary path.

### 2.8 `GaussQuadFullSphere` retirement (expanded scope)

Reading `rt_set_streams.jl:63-86`: `GaussQuadFullSphere` constructs `gausslegendre(2·Nquad)` over `[-1, 1]` and discards the lower hemisphere. The retained nodes/weights are *not* an optimal Gauss-Legendre quadrature for `[0, 1]` — they're half of an optimal full-sphere rule. For RT problems where upper-hemisphere integration is what matters, it's strictly worse than `GaussLegQuad`.

**Audit scope** (GPT Issue 8 corrected v0.3 undersizing):

Files touched by the retirement:
- `src/CoreRT/types.jl:117` — type definition
- `src/CoreRT/CoreRT.jl:154` — export
- `src/CoreRT/tools/rt_set_streams.jl:53,63` — dispatched method
- `src/CoreRT/DefaultParameters.yaml:12` — **currently the default!**
- `test/test_parameters/ThreeBandsParameters.yaml:73`
- `docs/src/pages/concepts/02_rt_theory.md:76`
- `docs/src/pages/IO/Schema.md:32,84`
- `docs/src/pages/IO/Overview.md:33,55`
- `docs/src/pages/geoschem_integration.md:79`
- `src/IO/Parameters.jl:94,99` — YAML loader

**Retirement plan:**

1. **v2.0 (Step 1A.0.5):**
   - `Base.depwarn` in the `rt_set_streams(::GaussQuadFullSphere, ...)` method (fires when the type is used regardless of source).
   - `DefaultParameters.yaml` default changes to `GaussLegQuad()`. **This is itself a behavior change for any user inheriting from the default.** Document in changelog.
   - `ThreeBandsParameters.yaml` test fixture updates.
   - Docs updates (4 files).
   - YAML loader keeps the string mapping with deprecation warning, still works.

2. **Audit reference benchmarks.** Any test currently using `GaussQuadFullSphere` gets regenerated with `GaussLegQuad`. Diff in changelog.

3. **v2.1:** remove the type, the `rt_set_streams` method, and the YAML mapping.

Asks Sanghavi for confirmation (Q3) given the YAML default change.

### 2.9 Source-type forward-compatibility (deferred refactor, no rename)

Today `DNI` and `SFI` types exist in `CoreRT/types.jl:120-126` as `AbstractSourceType` subtypes but are essentially vestigial — the `SFI::Bool` parameter does the real work, and `"DNI"`/`"SFI"` strings aren't in the YAML schema (GPT Issue 4 confirmed).

Eventually the source axis needs proper dispatch (`SolarSource`, `ThermalSource`, `LunarSource`, `CompositeSource`) for thermal/MW work — but **not in this round.** v0.4 keeps the existing type names; the rename was dropped from v0.3 because it would create churn without removing the real boolean path.

The forward-compat invariant remains:

> *Any code touching solar geometry or source-term assembly in Step 1 must funnel through `ctx.source`, not through hardcoded `μ₀`/`I₀` references or positional arguments.*

Concretely:

1. `RTContext` carries a `source::AbstractSourceType` field, parameterized.
2. Today populated with an internal `DefaultSolarSource{FT}` value carrying current `μ₀, I₀` — *this is an internal implementation detail, not a user-visible type rename*.
3. Strategy methods (SS corrections, etc.) read solar geometry through `ctx.source` only.
4. `SFI::Bool` parameter retained unchanged.

When thermal RT eventually lands (post-v2.0), the path is: add `ThermalSource` type, add `source_term(::ThermalSource, ctx)` methods, retire `SFI::Bool` and the existing `DNI`/`SFI` types together with a proper migration. No re-architecture needed because Step 1's strategy methods already go through `ctx.source`.

The `AbstractLayerKind` work in §5 (`NoScattering` dispatch) is the *complementary* future-compat piece — once both the source axis and the layer-kind axis dispatch correctly, thermal RT slots in by writing one new source type and one new layer-kind method.

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

1. **`isa` gate.** Adding any new BRDF needing the same correction (Cox-Munk-Bréon, glitter) means another `if`.
2. **Inconsistent across entry points.** Present in `rt_run.jl`, missing from `rt_run_lin.jl` and `rt_run_multisensor.jl`. Missing from `rt_run_lin.jl` is a silent Jacobian bug at glint geometries.
3. **`if brdf isa CoxMunkSurface` hides physics.** The same correction is needed for δ-M-truncated aerosol forward peaks. vSmartMOM has the first half (`delta_m_truncation.jl`) but not the second.

### 3.2 What VLIDORT does

VLIDORT factors single-scattering into **FO — First-Order Code**. Peer of the multiple-scatter solver. Control surface in `vlidort_inputs_def.f90:444-481`:

- `DO_FOCORR` — master on/off
- `DO_FOCORR_EXTERNAL` — external pre-computed exact-SS hook
- `DO_FOCORR_NADIR` / `DO_FOCORR_OUTGOING` — sphericity treatment
- `DO_SSCORR_TRUNCATION` — whether SS correction also accounts for δ-M truncation
- `DO_SSCORR_USEFMAT` — exact phase function from F-matrix vs. Greek coefficients

Lessons:
1. SS correction is a **subsystem with its own master flag**, not a per-surface afterthought.
2. **Not surface-specific.** Same pass corrects truncated aerosol forward-peaks and surface specularities.
3. **External** hook for user-supplied corrections.

### 3.3 Strategy hierarchy

```julia
abstract type AbstractSSCorrection end

struct NoSSCorrection <: AbstractSSCorrection end

struct BRDFSpecularCorrection{S<:AbstractSurfaceType} <: AbstractSSCorrection
    surface::S
    nQuad_BRDF::Int
end

struct AerosolForwardPeakCorrection{FT} <: AbstractSSCorrection
    fᵗ_threshold::FT
end

struct ExternalSSCorrection{F} <: AbstractSSCorrection
    callback::F
end
```

Per-surface trait controlling defaults:

```julia
has_specular_peak(::AbstractSurfaceType) = false
has_specular_peak(::CoxMunkSurface)      = true
```

### 3.4 Storage and the call site

**Storage decision** (GPT Issue 6): `ss_corrections` does **not** live as a field on `SolverConfig`. Instead, it's computed via `default_ss_corrections(model)` at `run_rt` time, with user override via kwarg:

```julia
run_rt(model; mode=ForwardMode(),
              outputs=(DirectionalTOA(), DirectionalBOA()),
              ss_corrections=default_ss_corrections(model),
              iBand=1)
```

Reasons:
- Avoids the type-stability dance of storing concrete tuples in a long-lived struct
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
        # AerosolForwardPeakCorrection NOT included until Step 1B.7 implements
        # the apply method. Type defined ≠ ready to enable.
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

For lin: `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::LinContext)` is **a separate method** updating both `R_SFI` and `Ṙ_SFI` consistently. Real implementation work, not a wiring one — the call site is one line; the dispatched method has its own derivative math. Multisensor follows the same pattern.

### 3.5 What this gives the linearization redesign

Each correction is its own self-contained unit of differentiation. Wrap `apply_correction!` in a ForwardDiff seed → get `∂(correction)/∂(params)` → add to chain rule via `lin_added_layer_all_params.jl`. Current implementation makes that hard because the correction lives *outside* the lin seam. With strategy structs, the hack becomes a typed unit the linearization framework can see.

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
| **When** | After Fourier loop completes | Inside the canopy direct-beam pipeline |
| **How** | Additive | Multiplicative on joint gap probability |
| **What it patches** | Fourier truncation losing specular-peak resolution | Beer-Lambert assuming independent gaps |

Conceptual similarity captured by **shared dispatch idiom**, not shared type.

### 4.3 The wiring: physical reasoning and seam choice

**Physical reasoning:** the hotspot is a correlated-path-pair phenomenon. `P_so(L)` describes the joint probability that the *same gap* is open along both sun and view paths. When sun and view coincide (`α → 0`), every gap that lets the photon in also lets it out — paths perfectly correlated.

**This correlation only exists for single-scattering events observed directly without further interaction.** Once a photon scatters once into the diffuse field, its outgoing direction decorrelates. So the hotspot correctly applies to:
- Direct-solar → direct-view contribution (single scattering observed straight back)

It does **not** apply to:
- Direct-solar → diffuse field → view contribution
- Diffuse-source → diffuse field → view contribution

The wiring requires direct-beam vs. diffuse separation at postprocessing — which is **Step 2's direct-beam subsystem** (§6). Not Step 1.

### 4.4 Step 1C does the field, Step 2 does the wiring

**In Step 1C (this round):**

1. Add `hotspot::HS` parameter to `CanopySurface`, where `HS<:AbstractHotSpot{FT}`. Default `NoHotSpot{FT}()`. **Parameterized on hotspot type** (GPT Issue 7) so dispatch is compile-time when methods reach `surface.hotspot`.
2. Plumb through `model_from_parameters` and YAML so `CanopySurface(..., hotspot=KuuskHotSpot(0.1))` is constructable.
3. Tests verifying field is constructable but runtime behavior is bit-for-bit identical to today (`NoHotSpot` is the default and no code path acts on it yet).
4. Document explicitly that wiring is deferred to Step 2.

~30 lines.

**In Step 2 (next round, §6):**

1. Implement direct-beam vs. diffuse separation in postprocessing as a real subsystem (analog of VLIDORT's FO).
2. Apply hotspot multiplicatively in postprocessing on the directional direct-beam term only.
3. Validate against VLIDORT FO tables and RAMI HOM26/HOM27.

⚠️ **For Sanghavi (Q2):** any current canopy benchmark tuned to match the no-hotspot result needs flagging before Step 2 lands.

---

<a id="5-step-1d"></a>
## 5. Step 1D — optical-properties algebra and layer-kind dispatch

This is the consolidated algebra chapter. Five pieces, all related: documentation + bug fix (5.2), dimensioned zero (5.3), contributor pattern (5.4), storage form as type parameter (5.5), `NoScattering` layer kind (5.6).

### 5.1 Motivation

`constructCoreOpticalProperties` in `compEffectiveLayerProperties.jl:11-65` hardcodes Rayleigh + aerosols + absorption as the contributors. Two consequences:

1. New contributors (clouds, ice, sea-salt with custom phase function) require editing the file. Rayleigh-vs-Cabannes branch (`_rayleigh_greek_source`) is also a special case here.
2. The thermal/MW case has no Rayleigh, but the current code uses Rayleigh as the loop initializer — there's no way to express "no Rayleigh" cleanly.

The operator algebra (`+` and `*`) on `CoreScatteringOpticalProperties` already does most of the right work — contributors *combine via `+`*. What's missing: making the contributor *list* the dispatched axis, formalizing the algebra for new contributors and reviewers, fixing existing bugs, and adding layer-kind dispatch so thermal/MW falls out naturally.

### 5.2 The existing algebra (what the docs need to say, plus the bug)

Five overloads in `types.jl:1083-1137`:

- `+(scattering, scattering)`: τ-and-ϖ-weighted average of phase functions, with `wx == 0` short-circuit (line 1099)
- `*(scattering, scattering)`: vertical/band concatenation via `cat` along dim 3
- `+(scattering, absorption)` and `+(absorption, scattering)`: thicken layer, leave Z untouched (commutative)
- `*(scalar, scattering)`: scale τ only, leave ϖ and Z untouched

**Bug at types.jl:1137** (GPT Issue 3 caught this): the scalar `*` overload references `y.G` but `y` is `CoreScatteringOpticalProperties` which has no `G` field — only `CoreDirectionalScatteringOpticalProperties` does. This is a real bug in existing code; calling `scalar * CoreScatteringOpticalProperties(...)` would error. It's never exercised because nothing actually calls scalar multiplication on the non-directional type. **Fix as part of Step 1D.1** — remove the `y.G` reference. Worth flagging that this kind of latent bug is exactly what the documentation work surfaces.

**Storage-form duality** (the part that's not formally documented):

- `CoreScatteringOpticalProperties` may carry Z as `(n_pq, n_pq, 1)` (compact, broadcast across the spectral grid) or `(n_pq, n_pq, nSpec)` (expanded, per-wavelength).
- All operators accept both via broadcasting; `expandOpticalProperties` converts compact → expanded explicitly when kernels need it.
- `+` of compact and expanded broadcasts correctly. `*` calls `expandOpticalProperties` on both first.

Algebraic properties:
- `+` is commutative and associative on scattering+scattering, commutative on scattering+absorption (no identity element today)
- `*` is associative but **not** commutative — vertical/band stacking, order matters

**Missing from the existing algebra (added in Step 1D.1):**

- `+(absorption, absorption)` — needed the moment per-species absorption contributors exist (separate O₂, H₂O, CO₂ contributors). Trivial: `CoreAbsorptionOpticalProperties(x.τ .+ y.τ)`.

### 5.3 Dimensioned zero

To start the contributor loop with empty `props` (so thermal/MW with no Rayleigh works), the algebra needs a zero element.

**Decision (per Q4):** dimensioned zero — a `CoreScatteringOpticalProperties` value with τ = 0. Reasons:

1. **The existing `+` short-circuit `all(wx .== 0.0)` at line 1099 already implements zero-element behavior.** When `x.τ = 0`, `wx = 0`, the short-circuit returns y's data. The current code *already implements* zero-element semantics; we just give it a name.
2. **No new sentinel type.** A sentinel `EmptyOpticalProperties` would need ~10 new method signatures for all `(Empty, Scattering, Absorption)` combinations.
3. **Compact form keeps the cost low.** Z is `(n_pq, n_pq, 1)` for the zero, ~256 KB per band per Fourier moment. Expanded only on demand.
4. **Zero participates only in `+`, never in `*`.** Contributor loop uses `+`; band concatenation via `*` happens between always-real values.

**Correct dimensions** (GPT Issue 2 — v0.3 was wrong about both):

```julia
function zero_optical_properties(ctx)
    n_pq = ctx.pol_type.n * ctx.Nquad   # NOT Nstreams — Z uses augmented quadrature
    CoreScatteringOpticalProperties(
        τ   = zeros(ctx.FT, ctx.nSpec),  # NOT n_layers — τ is per-wavelength
        ϖ   = zero(ctx.FT),               # placeholder; never read because τ=0
        Z⁺⁺ = zeros(ctx.FT, n_pq, n_pq, 1),  # compact form
        Z⁻⁺ = zeros(ctx.FT, n_pq, n_pq, 1),
    )
end
```

The construction loop iterates *over layers* in the outer dimension; each layer's `CoreScatteringOpticalProperties` has `τ` of length `nSpec` (the spectral grid).

**Pure-absorption case** (no scattering contributors anywhere — plausible for thermal IR with only absorption): the algebra produces a degenerate `CoreScatteringOpticalProperties` with τ from absorption, ϖ=0, Z=0 from the zero element. **This is not the same as `CoreAbsorptionOpticalProperties`** (GPT correctly caught v0.3's overclaim).

The handling: §5.6's `NoScattering` layer-kind dispatch *correctly classifies* this case at layer construction (`ϖ ≈ 0 → NoScattering`), and the kernel dispatches to a closed-form Beer-Lambert path that ignores Z entirely. So the degenerate scattering type is *numerically equivalent* to absorption-only behavior, achieved via dispatch rather than via a separate type. **One pipeline; the layer-kind dispatch handles the regime difference.**

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

Implementations:

```julia
function optical_properties(c::RayleighContributor, ctx)
    Z⁺⁺, Z⁻⁺ = compute_Z_moments(ctx.pol_type, ctx.μ, c.greek_source[ctx.iBand],
                                  ctx.m, arr_type=ctx.arr_type)
    CoreScatteringOpticalProperties(c.τ_rayl, c.ϖ, Z⁺⁺, Z⁻⁺)
end

function optical_properties(c::AerosolContributor, ctx)
    Z⁺⁺, Z⁻⁺ = compute_Z_moments(ctx.pol_type, ctx.μ, c.aerosol_optics.greek_coefs,
                                  ctx.m, arr_type=ctx.arr_type)
    if c.apply_delta_m
        (; fᵗ, ω̃) = c.aerosol_optics
        τ_mod = (1 - fᵗ * ω̃) .* c.τ_aer
        ϖ_mod = (1 - fᵗ) * ω̃ / (1 - fᵗ * ω̃)
        CoreScatteringOpticalProperties(τ_mod, ϖ_mod, Z⁺⁺, Z⁻⁺)
    else
        CoreScatteringOpticalProperties(c.τ_aer, c.aerosol_optics.ω̃, Z⁺⁺, Z⁻⁺)
    end
end

optical_properties(c::AbsorptionContributor, ctx) =
    CoreAbsorptionOpticalProperties(c.τ_abs)
```

Construction:

```julia
function band_optical_properties(atm::Atmosphere, iBand, m, ctx)
    n_layers = length(atm.profile)
    layer_props = Vector{CoreScatteringOpticalProperties}(undef, n_layers)
    for iz in 1:n_layers
        local_ctx = (; ctx..., iBand, m, iz)
        props = zero_optical_properties(local_ctx)
        for c in atm.contributors
            props += optical_properties(c, local_ctx)
        end
        layer_props[iz] = props
    end
    return layer_props
end
```

Per-band concatenation via `*` happens at the outer level, unchanged.

**Scope** (GPT Issue 3 corrected): **the contributor pattern is for atmospheric layers only.** `CoreDirectionalScatteringOpticalProperties` (used by canopy in `canopy_surface.jl:408`) has no `+`/`*` algebra and is its own pipeline. Step 1D doesn't refactor canopy's optical-properties construction. Canopy stays as a surface boundary condition with bespoke handling.

### 5.5 Storage form promoted to type parameter

`expandOpticalProperties` converts compact `(n_pq, n_pq, 1)` to expanded `(n_pq, n_pq, nSpec)`. Today the *type* doesn't track which form is in play — so kernels assume the worst and call `expandOpticalProperties` defensively, the compiler can't specialize, and there's no type-level enforcement.

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

Operators dispatch on storage form:

```julia
+(x::CoreScatteringOpticalProperties{Compact}, y::CoreScatteringOpticalProperties{Compact}) = ...    # stays compact
+(x::CoreScatteringOpticalProperties{Compact}, y::CoreScatteringOpticalProperties{Expanded}) = ...   # promotes to expanded
+(x::CoreScatteringOpticalProperties{Expanded}, y::CoreScatteringOpticalProperties{Compact}) = ...   # mirror
+(x::CoreScatteringOpticalProperties{Expanded}, y::CoreScatteringOpticalProperties{Expanded}) = ... # per-wavelength

expandOpticalProperties(x::CoreScatteringOpticalProperties{Compact}, nSpec) :: CoreScatteringOpticalProperties{Expanded}
expandOpticalProperties(x::CoreScatteringOpticalProperties{Expanded}, _) = x   # no-op
```

What this gives:
- Compiler specializes kernel code for both paths
- Kernels requiring expanded form can require it at the type level
- `expandOpticalProperties` becomes a type-level transformation, not a runtime check

Cost: 4 `+` overloads instead of 1 (with broadcast). Same total LOC because the broadcast handling now lives in named methods rather than implicit.

### 5.6 `NoScattering` layer kind

The current code has implicit "is there scattering?" branching in `elemental.jl`, `doubling.jl`, `interaction.jl`. There's already `AbstractScatteringInterface` with `ScatteringInterface_00/01/10/11` for layer-pair dispatch — the *infrastructure exists* but the *single-layer kind* doesn't.

**Promote to dispatch:**

```julia
abstract type AbstractLayerKind end
struct Scattering   <: AbstractLayerKind end
struct NoScattering <: AbstractLayerKind end
```

Layer carries kind as a type parameter (zero runtime cost):

```julia
struct Layer{LK<:AbstractLayerKind, ...}
    optical_properties::AbstractOpticalProperties
    # ... other fields ...
end
```

Layer kind determined at construction:

```julia
layer_kind(props::CoreScatteringOpticalProperties) = props.ϖ ≈ 0 ? NoScattering() : Scattering()
layer_kind(props::CoreAbsorptionOpticalProperties) = NoScattering()
```

The runtime check happens **once at layer construction**, then dispatch propagates. Per §1.5 ("branches that stay"), this is the right kind of branch — numerical predicate at construction, type dispatch from there.

Kernels dispatch:

```julia
function elemental!(::Scattering, layer, ctx)
    # full doubling/elemental as today
end

function elemental!(::NoScattering, layer, ctx)
    # closed-form: r = 0, t = exp(-τ/μ) * I
    # if thermal source: j = B(T) * (1 - exp(-τ/μ))
    # if no thermal source: j = 0
end

function doubling!(::Scattering, layer, ctx)
    # full doubling
end

function doubling!(::NoScattering, layer, ctx)
    # trivial: t' = t·t, r' = 0
end
```

The existing `ScatteringInterface_00/01/10/11` dispatch at the layer-pair level continues to work — when both layers are `NoScattering`, the resulting interface is `_00`, which already has the simplified composition method.

**What this unlocks:**

- **Eliminates `if ϖ == 0` numerical guards in elemental** — the no-scattering case is its own dispatched method
- **Thermal-only RT becomes a natural fall-out** — `NoScattering + ThermalSource` (when thermal lands) dispatches to a closed-form Beer-Lambert + Planck source kernel
- **Pure-absorption atmospheres are correctly handled** — the contributor algebra produces a degenerate `CoreScatteringOpticalProperties`, but `layer_kind` correctly classifies it as `NoScattering`, and the kernel dispatches to the simplified path. **One solver, dispatched kernels** — no separate "thermal entry point" needed. (This corrects v0.3's framing.)
- **Composes with the existing `ScatteringInterface_*` dispatch at no cost** — those types handle the layer-pair combinations already

Implementation cost: ~80–120 lines total. Worth doing in Step 1 because it removes a class of runtime branches, makes the thermal RT future trivial, and validates the dispatch rule at the kernel level.

### 5.7 Open design choice: `fScattRayleigh` extraction

The current code computes `fScattRayleigh = [Array(rayl[i].τ ./ combo[i].τ) ...]` as a Rayleigh-specific quantity used downstream by `_expand_layer_rayleigh!` for inelastic Raman. Two options:

- **(A) Trait-based.** `f_scattering(c::AbstractOpticalContributor, props, ctx) -> Vector{FT}` defined per contributor.
- **(B) Push to use site.** Compute `fScattRayleigh` inside `_expand_layer_rayleigh!` directly from the Rayleigh contributor and assembled props.

Recommended **(B)** — only consumer; no scaffolding for hypothetical future cases. Open for review (Q9).

### 5.8 Operator style preserved

The existing `+` and `*` operators stay exactly as written today (modulo the bug fix at line 1137). The contributor pattern adds *one new abstract type* (`AbstractOpticalContributor`) and *one new trait* (`optical_properties`). Storage-form dispatch refines the existing operators with explicit form-promotion methods. Layer-kind dispatch is at the kernel level, downstream of the algebra.

No semantic change to `+` or `*` for the standard cases; the algebra is *richer and more typed*, not *different*.

---

<a id="6-step-2"></a>
## 6. Step 2 — direct-beam separation subsystem

This is real solver work, not API consolidation. Sized between Step 1 and Step 3 because:

1. It's new physics with its own validation surface (VLIDORT FO).
2. It enables the canopy hotspot wiring deferred from Step 1C.
3. It's cleaner reviewed independently — Sanghavi can validate the physics; Step 3 can build on top.

### 6.1 What it does

VLIDORT factors single-scattering into a first-order code (FO) that runs separately from the multiple-scatter solver. The FO computes the exact direct-beam contribution: the unscattered solar beam attenuated through the atmosphere to the surface, reflected by the surface (or its single-scatter equivalent), back to the sensor — all in closed form, no doubling needed.

vSmartMOM today does *not* have a separate direct-beam computation. The direct-beam contribution is mixed into the multiple-scatter source terms `j₀±` in `elemental!` and `elemental_canopy!`. This is functionally correct but loses information — once mixed in, you can't recover "what was direct-beam vs. what was diffuse" without rerunning.

**Step 2 introduces a direct-beam subsystem that:**

1. Computes the direct-beam contribution analytically (closed form: `(μ₀ · I₀ / π) · exp(-τ_total/μ₀) × surface_term × exp(-τ_total/μ_v)` for each viewing geometry, with appropriate single-scatter phase function evaluation).
2. Subtracts the analytic SS contribution from the standard MS solve, so the MS solver computes only the diffuse contribution.
3. Adds the analytic SS back at postprocessing — exact at the actual viewing geometry (no Fourier truncation), like VLIDORT's FO.
4. Exposes both the direct-beam and the MS-diffuse as accessible quantities (eventually as Step 3 outputs).

### 6.2 Why this is solver work, not just an output

GPT Issue 7 was right: separating direct-beam is not "add an output accumulator." It requires *changing the kernel's source-term assembly*. The current `j₀±` assembly mixes contributions; the new design has the kernel either (a) skip the direct-beam contribution to `j₀±` entirely (because it'll be added analytically later) or (b) record it separately. Either way, the kernel is changed. The `DirectionalDirectBeam` *output* is the visible artifact; the *subsystem* is the real work.

Validation is against VLIDORT FO, not against existing benchmarks. Real physics derivation, real validation boundary.

### 6.3 What this enables

**Canopy hotspot wiring.** With direct-beam available separately, the hotspot factor multiplies the direct-beam term only at postprocessing — clean, local, matches VLIDORT/Verstraete-Pinty-Walthall convention.

**Improved SS correction.** The current Cox-Munk SS correction reconstructs `M_fourier` by re-evaluating Fourier coefficients (wasteful). With the direct-beam subsystem, the "correction" becomes "the direct-beam computation already gives you exact SS — just use it." Replaces the existing Cox-Munk hack with the proper physics.

**FO-style sphericity treatment.** VLIDORT's FO has separate methods for plane-parallel vs. pseudo-spherical attenuation. The direct-beam subsystem provides a natural place to add this without disturbing the MS solver.

### 6.4 Scope and validation

- **Plane-parallel Lambertian** — analytic single-scatter has closed form; reference is a one-line formula.
- **Plane-parallel Cox-Munk** — VLIDORT FO at varying wind speeds and glint geometries.
- **Plane-parallel canopy with `KuuskHotSpot`** — RAMI HOM26/HOM27.
- **Plane-parallel mixed atmosphere** — Rayleigh + aerosol + Lambertian, validated against VLIDORT FO direct output.

Step 2 is sized as ~5-7 substeps including the validation harness work (which uses the §11 infrastructure). It's the largest step in the plan.

---

<a id="7-step-3"></a>
## 7. Step 3 — output products and RT modes as dispatched strategies

This is the API consolidation. Builds on Steps 1 and 2.

### 7.1 The leak

`rt_run.jl:136-140` unconditionally allocates `R`, `T`, `R_SFI`, `T_SFI`, `ieR_SFI`, `ieT_SFI`, `hdr`, `bhr_uw`, `bhr_dw`. Return tuple at line 325 is conditional on `SFI`, but the allocations aren't. 1004 lines across `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` for what is structurally one loop.

### 7.2 Two orthogonal axes

The linearization redesign isn't just "another output." It changes model inputs, layer construction, composite storage, kernel calls, the Fourier loop body itself. That's a **mode**, not an output.

```julia
abstract type AbstractRTMode end
struct ForwardMode <: AbstractRTMode end
struct LinMode{PL} <: AbstractRTMode
    parameter_layout::PL
end

abstract type AbstractRTOutput end
struct DirectionalTOA <: AbstractRTOutput end
struct DirectionalBOA <: AbstractRTOutput end
struct DirectionalDirectBeam <: AbstractRTOutput end   # Step 2 enables; Step 3 exposes
struct HDRFOutput     <: AbstractRTOutput end
struct BHROutput      <: AbstractRTOutput end
struct MultiSensorOutput{...} <: AbstractRTOutput
    sensor_levels::Vector{...}
end
struct JacobianTOA <: AbstractRTOutput end       # paired with LinMode

allocate_output(::DirectionalTOA, model, nSpec) = ...
accumulate_output!(state, ::DirectionalTOA, ctx) = ...
finalize_output(state, ::DirectionalTOA, ctx) = ...
```

### 7.3 Public API

```julia
run_rt(model;
       mode=ForwardMode(),
       outputs=(DirectionalTOA(), DirectionalBOA()),
       ss_corrections=default_ss_corrections(model),
       iBand=1)
```

**Defaults vs. legacy compatibility** (GPT Issue 5):

- `run_rt`'s defaults are the *new clean API* — `(DirectionalTOA(), DirectionalBOA())`, return type is `NamedTuple` with named output fields.
- `rt_run(model)` becomes a compat wrapper that internally calls `run_rt` with the full legacy output set and assembles the legacy 7-tuple shape `(R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw, bhr_dw)` on return. Existing user code keeps working.
- `rt_run_lin` and `rt_run_multisensor` similarly become compat wrappers.

The two APIs have different defaults and different return shapes by design. Wrappers bridge.

### 7.4 What this collapses

| Current entry point | Becomes |
|---|---|
| `rt_run(model)` | `run_rt(model, outputs=(DirectionalTOA(), ...))` (with all legacy outputs requested) |
| `rt_run` with `SFI=true` and RAMI | `run_rt(model, outputs=(DirectionalTOA(), HDRFOutput(), BHROutput()))` |
| `rt_run_multisensor.jl` (177 lines) | `run_rt(model, outputs=(MultiSensorOutput(levels),))` |
| `rt_run_lin.jl` (303 lines) | `run_rt(model, mode=LinMode(layout), outputs=(DirectionalTOA(), JacobianTOA()))` |

1004 lines collapse into one coordinator (~150 lines) plus per-output strategy methods.

### 7.5 Migration substeps

1. Define `AbstractRTOutput`, `AbstractRTMode`. No call sites changed.
2. Implement `DirectionalTOA`, `DirectionalBOA`. Both dispatch to existing `postprocessing_vza!` internals. No behavior change.
3. Refactor current `rt_run` to take `outputs::Tuple` argument. Introduce `run_rt` as the new entry point.
4. Implement `HDRFOutput`, `BHROutput`. Migrate RAMI workflows.
5. Implement `MultiSensorOutput`. `rt_run_multisensor` becomes wrapper.
6. Implement `DirectionalDirectBeam` (uses Step 2's subsystem). Wire canopy hotspot.
7. Implement `LinMode` + `JacobianTOA`. `rt_run_lin` becomes wrapper. **Riskiest substep** — most lin methods to move.

Each substep is independently shippable; benchmarks at each boundary.

---

<a id="8-extension-cookbook"></a>
## 8. Extension cookbook (aspirational)

Recipes describing what the API will look like after Steps 1, 2, 3 land. Aspirational — they're a target for the architecture and a contributor reference.

### Recipe: adding a new BRDF surface kernel

```julia
struct MyBRDF{FT} <: AbstractSurfaceType
    # ... your parameters ...
end

function reflectance(brdf::MyBRDF, n_stokes, μᵢ, μⱼ, ϕ; ...)
    # ... your physics ...
end

n_fourier_moments(::MyBRDF, ctx) = 2 * ctx.Nstreams   # or whatever
has_specular_peak(::MyBRDF) = true   # or false
```

That's it. No edits to `rt_run`, `model_from_parameters`, or any *_lin variants.

### Recipe: adding a new optical contributor (e.g. cloud water)

```julia
struct CloudWaterContributor{FT, M} <: AbstractOpticalContributor
    cloud_optics::M
    τ_cloud::Vector{FT}
    apply_delta_m::Bool
end

function optical_properties(c::CloudWaterContributor, ctx)
    # ... your conversion ...
end
```

User code:

```julia
atm = Atmosphere(profile, spec_bands, contributors = (
    RayleighContributor(...),
    AerosolContributor(...),
    CloudWaterContributor(...),
    AbsorptionContributor(...),
))
```

### Recipe: adding a new SS correction type

```julia
struct MyCorrection <: AbstractSSCorrection
    # ... parameters ...
end

function apply_correction!(c::MyCorrection, ctx)
    # ... mutate ctx.R_SFI in place ...
end

needs_correction(c::MyCorrection, ctx) = ...
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

### Recipe: adding a new validation benchmark

(See §11 for the harness pattern.)

```toml
# test/reference_data/published/my_new_case.metadata.toml
[source]
type = "published_full_ms"
citation = "Author et al. 2024, JQSRT N(M), pp-pp, Table X"
doi = "10.xxxx/yyyy"
regime = ["thin", "thick"]

[atmosphere]
# ...

[expected]
data_file = "my_new_case.h5"
# ...
```

One TOML file per benchmark. Run via `julia test/run_reference_validation.jl`. New benchmark = new file.

### Recipe: adding a new layer kind

(Genuinely speculative — `Scattering` and `NoScattering` are likely the only two for a long time. But the pattern is documented.)

```julia
struct ApproximatePartialScattering <: AbstractLayerKind end

# Define how this kind is recognized:
function layer_kind(props::CoreScatteringOpticalProperties)
    if props.ϖ ≈ 0
        NoScattering()
    elseif props.ϖ < 0.1
        ApproximatePartialScattering()  # if some clever approximation is worth it here
    else
        Scattering()
    end
end

# Define the kernel methods:
function elemental!(::ApproximatePartialScattering, layer, ctx)
    # ... your specialized kernel ...
end
function doubling!(::ApproximatePartialScattering, layer, ctx)
    # ...
end
```

---

<a id="9-type-gates"></a>
## 9. Remaining type-gates inventory (migration debt)

Beyond Step 1, 2, 3, these places need future migration:

| Smell | Files | Replace with | Status |
|---|---|---|---|
| `if brdf isa CanopySurface` (4 occurrences) | `rt_run.jl:170, 175, 182, 245` | `requires_canopy_setup(::AbstractSurfaceType)` trait | Step 3 |
| `if RS_type isa AbstractRamanType` style branches | `interaction_inelastic.jl`, `elemental_inelastic.jl` | Already partially dispatched — finish it | Future |
| `arch isa GPU` checks | scattered (`array_type`, `synchronize_if_gpu`) | Already mostly dispatched on `AbstractArchitecture`; audit holdouts | Future |
| `if SFI` boolean threading | every `rt_run*.jl`, `postprocessing_*.jl` | Goes away under Step 3 | Step 3 |
| Polarization-type branching via `pol_type.n` | many places | Numerical property — likely fine, audit only | Audit only |
| `Stokes_IQU` implementation completeness | scattered | IQU code paths most likely to have rot | Future |
| `AbstractSourceType` resurrection | `types.jl:120-126` | `SolarSource`, `ThermalSource`, etc. | Future (thermal RT trigger) |

The forward-compat invariant from §2.9 (all solar geometry through `ctx.source`) keeps the `AbstractSourceType` resurrection mechanical when it happens.

---

<a id="10-future-architecture"></a>
## 10. Future architecture (deferred)

Flagged for post-v2.0 roadmap. Each is the natural next application of the dispatch rule but out of scope for this round.

**C1 — Optical-properties algebra extension** (now Step 1D, see §5). *In scope this round.*

**C2 — Layer-kind dispatch** (now §5.6). *In scope this round.*

**C3 — Validation-as-data** (now §11). *In scope this round as enabling work for Step 2.*

**C4 — Layer-kernel interface stabilization.** Doubling/elemental/interaction kernels coupled to specific layer-storage types. Each new layer flavor (RRS, VRS, future inelastic) requires a new kernel variant. Future-proof: kernels operate on `AbstractLayerStorage` interface. Belongs in linearization redesign track.

**C5 — Physics audit traits.** A `physics_summary(model)` returning structured description of every active assumption, dispatched through every component. Useful for reproducibility, sanity checks, downstream tooling. Cheap to add given dispatch architecture.

**C6 — Mie greek-coef interpolation.** Aerosol-specific. Today, every spectrum point requires its own Mie computation. Smarter: compute greek coefs at sparse "anchor" wavelengths, interpolate in size-parameter space (greek coefs are smooth in size parameter). The contributor pattern (§5.4) handles this naturally — `AerosolContributor` carries either pre-computed greek coefs *or* a `GreekCoefInterpolation` object. Optional optimization, no architectural change.

**C7 — Versioned interface contracts.** External BRDF/aerosol/surface implementers need stable contracts across vSmartMOM minor versions. Documented interface methods with `Required`/`Optional` annotations and a `check_interface(MyBRDF)` test helper. Compare to `AbstractArrays`. Matters when vSmartMOM grows an ecosystem.

**C8 — Thermal RT.** `ThermalSource` adds Planck per-layer source dispatching through `source_term(::ThermalSource, ctx)`. With `NoScattering` layer-kind already in place from Step 1D and the `ctx.source` invariant from §2.9, this slots in mechanically: write `source_term` methods for thermal, retire `SFI::Bool`. No new architecture.

---

<a id="11-validation"></a>
## 11. Validation infrastructure

### 11.1 The framing distinction

vSmartMOM is a full multiple-scattering code (matrix operator method, doubling-adding to convergence). **Reference tables for validation must be from full-MS codes.**

2OS approximations (Natraj-Spurr 2007) are *not appropriate* for validating vSmartMOM, even though they share VLIDORT lineage. 2OS computes `I_1 + I_2` exactly and neglects orders 3+; vSmartMOM computes `I_total = I_1 + I_2 + I_3 + ...`. For optically thick cases, the difference is *signal*, not error. Validating against 2OS would penalize vSmartMOM for capturing higher-order scattering correctly.

Per-regime appropriateness:

| Reference type | What it validates | Where it applies |
|---|---|---|
| Full-MS reference | Total radiance | Always |
| Single-scatter-only reference (e.g. VLIDORT FO output) | Direct-beam component, SS correction | Step 2 only |
| 2OS reference | (nothing in vSmartMOM directly) | Don't use |
| Analytic Lambertian SS | Direct-beam in trivial case | Step 2 sanity |

### 11.2 The metadata-plus-data-file pattern

Each reference benchmark is a TOML metadata file plus an HDF5 (or CSV) data file with the numerical reference output. One harness iterates over all of them.

```toml
# test/reference_data/published/coulson_dave_sekera_table_iv.metadata.toml
[source]
type = "published_full_ms"
citation = "Coulson, Dave & Sekera 1960, Tables Related to Radiation Emerging from a Planetary Atmosphere with Rayleigh Scattering, U. California Press, Table IV"
regime = ["thin", "thick"]
notes = "Rayleigh-only canonical reference; covers τ from 0.05 to 1.0"

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

The data file holds just numerical content extracted from the published table. Citation in metadata makes provenance transparent.

### 11.3 Type-dispatchable validation

Apply the dispatch rule to the harness itself:

```julia
abstract type AbstractValidationRegime end
struct FullMS <: AbstractValidationRegime end
struct SingleScatterOnly <: AbstractValidationRegime end
struct TwoOrdersOnly <: AbstractValidationRegime end

validate_against_reference(::FullMS, vsmartmom_output, expected) = ...   # direct comparison
validate_against_reference(::SingleScatterOnly, vsmartmom_output, expected) = ...   # compare SS component only
function validate_against_reference(::TwoOrdersOnly, vsmartmom_output, expected)
    @warn "Skipping 2OS reference; vSmartMOM is full-MS"
end
```

Cross-regime comparisons are caught by the type system. Adding a new validation regime is one new method.

### 11.4 Reference set for v2.0

**Published, committed to the repo:**

- **Coulson, Dave & Sekera 1960** — Rayleigh-only canonical reference. Full-MS. τ ∈ [0.05, 1.0]. Stokes I, Q.
- **Garcia & Siewert 1989** — Benchmark suite with varying complexity. Full-MS.
- **Spurr 2002** VLIDORT validation appendix — Full-MS cross-validation against multiple codes.
- **Wauben & Hovenier 1992** — Vector RT with aerosols. Full-MS.

**Locally-generated (gitignored, dev-machine only):**

- VLIDORT distribution `results/` directory output — apples-to-apples with direct VLIDORT comparison. Full-MS.
- VLIDORT FO output for Step 2 direct-beam validation. Single-scatter-only.
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

### 11.5 What enables this in Step 1

- Step 1A.0.6 — write the metadata format spec, the harness, and one test case (Coulson-Dave-Sekera) as proof-of-concept.
- Step 1B.4 onward — every Step 1 substep that changes numerics gets a validation against the published set.
- Step 2 onward — the VLIDORT FO validation uses the harness from day one.

### 11.6 No GPL contagion

VLIDORT is GPL. Published reference tables (Coulson-Dave-Sekera, Garcia-Siewert, Wauben-Hovenier, Spurr 2002 appendix) are *facts in published papers*, not VLIDORT source code, and can be committed freely. The VLIDORT distribution `results/` directory contains numerical outputs (also facts). Locally-running VLIDORT to generate additional cases is a private dev-machine activity — those outputs go in `gitignored/`. No GPL touches the vSmartMOM repo.

---

<a id="12-execution-plan"></a>
## 12. Execution plan

### Step 1 — Fourier work, optical-properties algebra, validation harness (this round)

| # | What | Risk | Validation |
|---|---|---|---|
| 1A.0 | Rename `AbstractFourierDecompositionType` → `AbstractMieDecompositionAlgorithm`. Add deprecation alias. | ~50–100 lines including tests, docs, IO imports. Touches exported API. | All existing tests pass. |
| 1A.0.5 | Deprecate `GaussQuadFullSphere`. `Base.depwarn` in `rt_set_streams` method. Update `DefaultParameters.yaml` default to `GaussLegQuad`. Update test fixture. Update 4 docs files. | **YAML default change** — behavior shift for users inheriting defaults. | Audit benchmarks for `FullSphere` use; document any diffs. |
| 1A.0.6 | Validation harness foundation: metadata-plus-data-file format, type-dispatchable validation regime, one published reference (Coulson-Dave-Sekera) as proof-of-concept. | New infrastructure. | Harness round-trips on the proof-of-concept reference. |
| 1A.1 | Introduce `FourierDecomposition` module skeleton. | Plumbing. | Module loads. |
| 1A.2 | Add `RTContext` struct (§1.3), including `source::AbstractSourceType` field (§2.9 forward-compat invariant). | New abstraction. | Trivial unit tests. |
| 1A.3 | Add `Nstreams` field to `QuadPoints` separate from augmented `Nquad`. | One new field, two value-producers. | All quadrature setups produce both values. |
| 1A.4 | Change `SolverConfig.max_m` to `Union{Nothing, Int}`. YAML schema accepts `null` or omission. | **Behavior change for existing configs** with explicit `max_m: 30`-style values — those continue as user overrides. | Existing tests pass. Document in changelog. |
| 1A.5 | Add `n_fourier_moments` trait + methods (§2.6). | New code. | Unit tests for each surface type. |
| 1A.6 | **Wire `n_fourier_moments(model, iBand)` into the consumers.** Replace `get_max_m(model)` with the trait call in `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`. | **Behavior change risk.** Cox-Munk bands may now compute more moments than before. | 6SV1 / Natraj reference benchmarks before and after; published set from harness. |
| 1A.7 | (Deferred) Promote `nQuad_BRDF` to per-surface field. | Trivial. | Convergence sweeps. |
| 1B.1 | Refactor BRDF Fourier integration into shared helper. | Pure refactor. | Bit-for-bit equivalence. |
| 1B.2 | Add `AbstractSSCorrection` hierarchy. Migrate Cox-Munk apply method. **Read `n_moments` from ctx.** | Pure relocation. | Bit-for-bit equivalence. |
| 1B.3 | Wire `apply_correction!` into `rt_run.jl` via `foreach`. Add `default_ss_corrections` (no `AerosolForwardPeakCorrection`). | Pure call-site swap. | Existing benchmarks unchanged. |
| 1B.4 | Add lin-context apply method updating `R_SFI` and `Ṙ_SFI` consistently. | New physics in lin path. | Validate Jacobians against finite-difference at glint. |
| 1B.5 | Add multisensor-context apply method, per-sensor-level attenuation. | New physics in multisensor path. | Validate against multisensor reference cases. |
| 1B.6 | (Deferred to later round) Implement `AerosolForwardPeakCorrection`. | New physics. | PyVLIDORT high-AOD cases. |
| 1B.7 | (Deferred) Cache Fourier coefficients for SS correction reuse. | Optimization only. | No behavior change. |
| 1C.1 | Add `hotspot::HS` parameter to `CanopySurface` (parameterized, not field with abstract type), default `NoHotSpot`. Plumb through `model_from_parameters` and YAML. | Pure plumbing; field constructable but no-op. | All existing tests pass bit-for-bit. |
| 1D.1 | Document optical-properties algebra (§5.2). Add `+(absorption, absorption)` overload. **Fix `y.G` bug at types.jl:1137.** | Mostly docs + 1 trivial overload + 1 bug fix. | All existing tests pass. |
| 1D.2 | Add dimensioned `zero_optical_properties` constructor (§5.3). Verify dimensions against actual usage (`τ` is `nSpec`, `n_pq` uses augmented `Nquad`). | New constructor. | Unit tests for zero-element behavior. |
| 1D.3 | Add `AbstractOpticalContributor` + `optical_properties` trait (§5.4). | New abstraction. | Unit tests for each contributor type. |
| 1D.4 | Migrate `constructCoreOpticalProperties` to use the contributor pattern. **Scope: non-directional algebra only.** Canopy stays separate. | **Touches the optical-properties pipeline.** | Bit-for-bit equivalence with current default contributors. |
| 1D.5 | Promote storage form to type parameter (`Compact` / `Expanded`). Refactor `+`/`*` overloads to dispatch on form. | Algebra restructure (4 overloads instead of 1). | Bit-for-bit equivalence. |
| 1D.6 | Add `AbstractLayerKind` with `Scattering`/`NoScattering`. Implement `layer_kind(props)`. Implement `elemental!(::NoScattering, ...)` and `doubling!(::NoScattering, ...)`. | New kernel methods. | Bit-for-bit equivalence on default scattering atmospheres; new pure-absorption test case. |

**Step 1 boundaries for benchmarks:**

- After 1A.0–1A.4: all existing tests pass bit-for-bit. Pure plumbing + YAML default change (1A.0.5).
- After 1A.6: 6SV1 / Natraj benchmarks for any band where `n_moments` could change.
- After 1B.3: existing Cox-Munk TOA reference cases.
- After 1B.4: lin Jacobian regression at glint geometries.
- After 1B.5: multisensor + Cox-Munk reference cases.
- After 1D.4: end-to-end equivalence test on a representative atmosphere.
- After 1D.6: bit-for-bit equivalence on scattering atmospheres; new test for pure-absorption layer.

### Step 2 — Direct-beam separation subsystem (next round)

Per §6.4. Substeps:

1. Define direct-beam computation interface (analytic SS evaluation per geometry).
2. Implement plane-parallel direct-beam for Lambertian surfaces. Validate against analytic formula.
3. Implement plane-parallel direct-beam for Cox-Munk surfaces. Validate against VLIDORT FO.
4. Hook direct-beam into kernel source-term assembly (subtract analytic SS from MS source, add back at postprocessing).
5. Implement canopy direct-beam with `KuuskHotSpot` dispatch. Validate against RAMI HOM26/HOM27.
6. Update SS correction methods to use direct-beam subsystem (replaces existing Cox-Munk hack).

### Step 3 — Output products and RT modes (final round)

Per §7.5. Seven substeps. Lin migration last and riskiest.

---

<a id="13-open-questions"></a>
## 13. Open questions

In rough order of "blocks the next commit" → "nice to settle eventually":

1. **Q1 (preamble)** — dispatch rule as governing principle, applied at every level.
2. **Q2 (preamble)** — canopy hotspot regression risk. Are there reference benchmarks tuned to the no-hotspot result? **Blocks step 1C.1 changelog work and Step 2.5.**
3. **Q3 (preamble)** — `GaussQuadFullSphere` retirement, including the YAML default change. **Blocks step 1A.0.5.**
4. **Q4 (preamble)** — `NoScattering` layer-kind dispatch as the unifying mechanism for thermal/MW. Match physics intuition? **Blocks step 1D.6 design.**
5. **Q5 (preamble)** — `n_fourier_moments(::CoxMunk) = 2 · Nstreams` derivation. Empirical sanity check ordering. **Blocks step 1A.5 validation plan.**
6. **Q6 (preamble)** — full-MS reference set for v2.0. Are there additional sources Sanghavi has been using? **Blocks step 1A.0.6 reference selection.**
7. **`l_trunc` semantics** — does vSmartMOM's `l_trunc` correspond to VLIDORT's "highest Legendre order in the truncated phase function expansion" or to VLIDORT's `NSTREAMS`? Empirical check via the Q5 sanity test.
8. **`fScattRayleigh` extraction** — trait-based (Option A) or push-to-use-site (Option B, recommended). §5.7.
9. **`ExternalSSCorrection` callback signature** — mutation-via-context vs. return-value API.
10. **`AerosolForwardPeakCorrection` threshold** — default `fᵗ > 0.01`. Reasonable starting point?
11. **`RTContext` lifecycle** — exact field set and m-mutation contract. §1.3.
12. **Naming harmonization** — existing `AbstractTruncationType`/`NoTruncation` vs. proposed `Abstract*Correction`/`No*Correction`.
13. **`m_max(::CanopySurface)` provisional vs. proper** — placeholder is `2·Nstreams`. Proper LAD-driven derivation lives in CanopyOptics.

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
- **Codex review of v0.1, v0.2, v0.3** and **GPT review of v0.3**, integrated throughout this document.
