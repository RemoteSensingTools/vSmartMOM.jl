# vSmartMOM.jl — architecture design and Fourier work

**Document version:** v0.3
**Status:** working draft — first version with substantive content for Sanghavi.
**Branch target:** `sanghavi-unified` (with `unified-vsmartmom` as the architectural merge base, physics layered in from `sanghavi`).
**Companion:** `LINEARIZATION_BUGS.md`, the linearization redesign track.

---

## Preamble — for Sanghavi

This document is a design proposal for vSmartMOM v2.0.0. It has been iterated through several rounds of code review (against the current `sanghavi-unified` branch, against VLIDORT 2.8.3 as the established reference, and against CanopyOptics.jl) but it has not yet had your input. Your sign-off on the **direction** matters more than the technical detail, because:

- You own the theory. The architectural choices below are decisions about how the *code* expresses the *physics* — but if the physics framing is wrong, the architecture is too. We need to know that early, not after a hundred commits.
- Several of the choices below are not technically forced. They reflect judgment calls about what kind of codebase vSmartMOM should be in two years. Your view of where the science is going should shape that.

### What this document is trying to do

The proposal has two parts that go together:

1. **Establish a single architectural rule** — multiple dispatch at physical and workflow boundaries; `if/else` only for scalar numerical cases — and apply it consistently across the codebase. This eliminates a class of bugs (the `if brdf isa CoxMunkSurface` smell that has caused silent inconsistencies between `rt_run.jl`, `rt_run_lin.jl`, and `rt_run_multisensor.jl`) and makes future extensions cheap (new BRDFs, new aerosol models, new outputs each become "define one struct, write one method").

2. **Apply this rule to four real, currently-broken or currently-missing pieces of physics**:
   - Per-component Fourier-moment count (currently surface-blind — Cox-Munk and Lambertian both get the same `max_m` from the aerosol Greek coefficients)
   - Single-scattering correction (currently a hack in `rt_run.jl:311` that's missing from the linearized and multisensor entry points)
   - Canopy hotspot (CanopyOptics has it; vSmartMOM doesn't wire it — every current canopy result is using independent Beer-Lambert silently)
   - Optical-properties algebra extension (the `+`/`*` operators on `CoreScatteringOpticalProperties` are elegant but undocumented; the algebra needs a true zero element for thermal/MW work where there's no Rayleigh)

### Five questions for you up front

These are the choices that affect everything else. We'd appreciate your view before we draft v0.4 or commit code. Each is one paragraph. Detailed treatment is later in the doc; the references are pointers, not required reading for the philosophical pass.

**Q1 — Dispatch architecture as governing principle.** The proposal makes "dispatch on types at physical/workflow boundaries; `if/else` only for numerical cases" the *rule* the codebase commits to. New BRDFs, new corrections, new outputs, new aerosol contributors all follow the same pattern (small concrete strategy struct, single-method dispatch). The alternative is to keep the current mixed style (some dispatch, some `isa` branches, some boolean parameters) and accept that bug-classes will keep recurring. Are you on board with the rule as a hard commitment going forward, or is there a class of physics where you think the rule is wrong-headed and the freedom to branch is worth keeping? *Detail: §1.*

**Q2 — Canopy hotspot is silently missing today; the fix has cost.** vSmartMOM's `elemental_canopy!` uses independent Beer-Lambert (`exp(-(k_s + k_o) L)`) where it should use joint gap probability (`exp(-(k_s + k_o) L) · C_hs(α, L)`). CanopyOptics ships `KuuskHotSpot` with ForwardDiff hooks but vSmartMOM doesn't call `joint_gap_probability` anywhere. Any current canopy benchmark of yours that compares against RAMI/PROSAIL hotspot-aware references is either wrong or was tuned to match the no-hotspot result. *We'd like to know if any such benchmarks exist before we quietly fix the bug — those would need to be regenerated, not silently corrected.* The proposal defers full hotspot wiring to Step 2 (it requires direct-beam vs. diffuse separation in postprocessing, which is naturally an output-strategy concern), but the fix is not free physics-wise. *Detail: §4.*

**Q3 — `GaussQuadFullSphere` retirement.** Christian flagged this as "a wild idea, not really correct." Reading the implementation in `rt_set_streams.jl:63-86`: it constructs `gausslegendre(2·Nquad)` over `[-1, 1]` and discards the lower hemisphere, but the retained nodes/weights are *not* an optimal Gauss-Legendre quadrature for `[0, 1]` — they're half of an optimal full-sphere rule, which is a different polynomial fit. For RT problems where upper-hemisphere integration is what matters (essentially all plane-parallel RT), it's strictly worse than `GaussLegQuad`. Proposal: deprecate in v2.0 (warning, still works), remove in v2.1. Do you concur, or is there a use case for `FullSphere` we should preserve? *Detail: §2.8.*

**Q4 — Optical-properties algebra: zero element and storage form.** The `+` and `*` overloads on `CoreScatteringOpticalProperties` are mathematically clean but undocumented. Two issues we want to address: (a) thermal/MW work won't have Rayleigh, so the algebra needs a defined zero — proposal is a *dimensioned zero* (a `CoreScatteringOpticalProperties` with τ = 0, exploiting the existing `wx == 0` short-circuit at `types.jl:1099`); (b) the duality between "compact" Z storage `(n_pq, n_pq, 1)` and "expanded" Z storage `(n_pq, n_pq, n_layers)` is real and exploited by the existing operators but never explained — proposal is to add a docstring section formalizing it. Does the dimensioned-zero approach match how you think about the algebra, or is there a structural reason to prefer a sentinel `EmptyOpticalProperties` type? *Detail: §5.*

**Q5 — `n_fourier_moments(::CoxMunk) = 2 · Nstreams` derivation.** Translating VLIDORT's `NMOMENTS = 2·NSTREAMS - 1` (Fortran inclusive loop, `2·NSTREAMS` iterations) into vSmartMOM's count convention gives `2 · Nstreams` where `Nstreams` is the *base* polar discrete-ordinate count (not `model.quad_points.Nquad`, which is augmented with VZA/SZA nodes). For `l_trunc = 60, Nstreams = 30`, the count is 60. Does this match your physical intuition for what Cox-Munk needs? The empirical sanity check is to run a Cox-Munk + Rayleigh test with the count set to 60 vs. 30 vs. 100 and check convergence; we'd value your view on what the reference numerics should look like before we run that test. *Detail: §2.*

### What this document is *not*

- Not a fait accompli. Every choice is explicitly revisable and the open-questions list at the end (§11) flags everything we know is open.
- Not yet implemented. The current `sanghavi-unified` branch has none of these changes; it has only the pre-existing fragmented architecture.
- Not the linearization redesign. That's a separate track in `LINEARIZATION_BUGS.md`. This proposal is upstream of it — landing the dispatch architecture cleanly makes the linearization redesign tractable, but doesn't replace it.
- Not exhaustive of every smell. Some pieces (e.g. `AbstractSourceType` for thermal/MW source dispatch) are flagged as future work and explicitly *not* refactored in this round, but the design is structured so resurrecting them later is mechanical, not architectural.

The rest of the document treats the technical detail. Sections 1–5 are the architecture and the four work items. Section 6 is the bigger Step 2 deduplication. Sections 7–9 are extension recipes, migration debt, and deferred future work. Section 10 is the concrete execution plan.

---

## Iteration log

- **v0.3** — substantial rewrite for Sanghavi's first read. Codex review of v0.2 folded throughout (real type names, `Nstreams` separated from augmented `Nquad`, `RTContext` lifecycle corrected, `AbstractRTMode` orthogonal to `AbstractRTOutput`, hotspot seam reconsidered). Added §5 (optical-properties algebra extension), §7 (extension cookbook, aspirational), §8 (type-gates inventory), §9 (deferred future work). `GaussQuadFullSphere` retirement included (§2.8). `AbstractSourceType` forward-compat invariant added (§2.9) without doing the refactor.
- **v0.2** — promoted dispatch rule to governing principle; two-step plan; Codex v0.1 review folded.
- **v0.1** — initial draft after VLIDORT/CanopyOptics/vSmartMOM code review.

---

## Table of contents

1. [The dispatch rule (governing principle)](#1-dispatch-rule)
2. [Step 1A — per-component `n_fourier_moments`](#2-step-1a)
3. [Step 1B — `AbstractSSCorrection` strategies](#3-step-1b)
4. [Step 1C — `AbstractHotSpot` field on `CanopySurface` (wiring deferred to Step 2)](#4-step-1c)
5. [Step 1D — optical-properties algebra: docs, dimensioned zero, contributor pattern](#5-step-1d)
6. [Step 2 — output products and RT modes as dispatched strategies](#6-step-2)
7. [Extension cookbook (aspirational)](#7-extension-cookbook)
8. [Remaining type-gates inventory (migration debt)](#8-type-gates)
9. [Future architecture (deferred)](#9-future-architecture)
10. [Execution plan](#10-execution-plan)
11. [Open questions](#11-open-questions)
12. [References](#12-references)

---

<a id="1-dispatch-rule"></a>
## 1. The dispatch rule (governing principle)

### 1.1 The rule

> **Use multiple dispatch at physical and workflow boundaries; keep `if/else` only for scalar numerical cases.**

Physical and workflow boundaries — places where the codebase currently uses `isa` checks or symbol-comparison conditionals — become method dispatch on small concrete strategy structs. These are:

- Surface physics (`if brdf isa CoxMunkSurface`, `if brdf isa CanopySurface`)
- Fourier bounds per component (currently surface-blind)
- Corrections — SS, hotspot, future thermal/Raman/multiple-scatter δ-M
- Output products (`if SFI`, `if RAMI` branches scattered through postprocessing)
- Output topology (separate `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` files for what is structurally one loop)
- Atmospheric optical contributors (currently hardcoded Rayleigh + aerosols + absorption in `constructCoreOpticalProperties`)
- *Source types* (deferred — see §2.9)

`if/else` is fine — and clearer than dispatch — for scalar numerical cases:

- Fourier normalization weight (`m == 0 ? 0.5 : 1.0`)
- Numerical guards (NaN, denorm, divide-by-zero)
- Cache state checks (`if isnothing(brdf._cache)`)
- Optional benchmark/diagnostic instrumentation
- **Runtime user-data predicates** (`if any(a -> a.fᵗ > 0.01, aerosols)` — data check on values, not type gate)

The rule eliminates *type and physics gates*, not all branches. Litmus test: if the branch's truth is determined at model-construction time by which types are in play, dispatch on those types. If determined at evaluation time by numerical values or user toggles, branch.

### 1.2 The strategy struct shape

Strategy objects are **small concrete structs held in tuples, dispatched via single-method functions.** Three implementation rules:

**(a) Concrete tuples, not abstract vectors.** A `Vector{AbstractSSCorrection}` forces runtime dispatch in inner loops. A `Tuple{BRDFSpecularCorrection{CoxMunkSurface{Float64}}}` lets the compiler specialize.

```julia
# Wrong — abstract eltype, runtime dispatch:
corrections = AbstractSSCorrection[BRDFSpecularCorrection(surf)]

# Right — concrete tuple, compile-time dispatch:
corrections = (BRDFSpecularCorrection(surf),)
```

**Honest about type stability:** when a strategy set is built conditionally and stored in a long-lived struct like `SolverConfig`, end-to-end compile-time dispatch requires either parameterizing the storage struct on the strategy-tuple type (invasive) or accepting runtime dispatch through a function-barrier when unpacked. **Pragmatic answer: runtime dispatch at band granularity is acceptable** — the strategy `apply_correction!` body is much more expensive (milliseconds) than the dispatch overhead (nanoseconds). Compile-time dispatch is achieved *inside* each strategy method via concrete type parameters, not at the outer driver level.

**(b) Single-method functions, not switch tables.** Each strategy owns its behavior:

```julia
apply_correction!(c::BRDFSpecularCorrection, ctx)
apply_correction!(c::AerosolForwardPeakCorrection, ctx)
apply_correction!(::NoSSCorrection, ctx) = nothing
```

Driver:

```julia
foreach(c -> apply_correction!(c, ctx), corrections)
```

**(c) Naming consistency across families.**

| Family | Trait | Apply |
|---|---|---|
| Fourier bounds | `n_fourier_moments(component, ctx)` | (count only) |
| SS corrections | `needs_correction(c, ctx)` | `apply_correction!(c, ctx)` |
| Optical contributors | `optical_properties(c, ctx)` | (composes via `+` algebra) |
| Output products | (presence in tuple is signal) | `accumulate_output!(state, output, ctx)` |
| Hotspot models | (carried as field) | `joint_gap_probability(model, …)` (in CanopyOptics) |

### 1.3 The `RTContext` object

Strategy methods need a consistent argument convention. `RTContext` is constructed once before the Fourier loop, mutated as `m` advances, passed to every dispatched method.

**Lifecycle:**

- Constructed once before the Fourier loop. Carries everything that doesn't depend on `m`.
- Field `m` is mutated as the loop advances. Strategy methods read it, don't mutate it.
- `composite_layer` accumulates **vertically within one Fourier moment** (across `iz` layers), reset at `iz == 1` on the next `m`.
- **Output arrays accumulate across `m`** — these belong to `AbstractRTOutput` strategies (Step 2).
- `τ_sum_all` is m-independent (the existing comment at `rt_run.jl:204`).

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
    source::S                  # see §2.9 — solar today, thermal/lunar/etc later
    pol_type::AbstractPolarizationType
    arch::AbstractArchitecture
    workspace                  # InteractionWorkspace, allocated once
    Nstreams::Int              # base polar streams (NOT augmented Nquad)
    n_layers::Int
    FT::Type{FT}
end
```

Exact field set is open — see Q10 in §11.

### 1.4 What the rule unlocks

- **Step 1 (§§2–5)** sketches the idiom on small contained pieces. Four work items, each <200 lines of new code.
- **Step 2 (§6)** applies the rule to output products, collapsing 1004 lines of forked `rt_run*.jl` into one coordinator + per-output strategy methods.
- **The linearization redesign** uses the same rule. Lin Jacobians are an `AbstractRTMode`, orthogonal to `AbstractRTOutput`. Not "another output."
- **Future corrections** (thermal SS, Raman SS, multiple-scatter δ-M, polarization SS for thermal) all fit the same `Abstract*Correction` shape without architectural work.

The rule has compounding returns. Step 1 is the first application; everything downstream gets cheaper.

### 1.5 Branches that stay (the "fine" list)

To preempt over-application:

- `m == 0` Fourier weight normalization (`weight = m == 0 ? 0.5/π : 1.0/π`). Numerical, stays.
- BRDF integration prefactor (`m == 0 ? FT(1.0) : FT(2.0)`). Numerical, stays.
- **Lambertian-surface kernel zeroing higher moments must stay** (`if m == 0` guard inside `lambertian_surface.jl:37`). The Fourier loop driver still runs `m > 0` for Rayleigh + aerosol contributions; the Lambertian *surface* must zero its own contribution because it has no higher moments. Numerical truth, not a physics gate.
- Cache initialization (`canopy._cache === nothing`). State, not type/physics dispatch.
- `if SFI` for source-function integration: prime target for Step 2 — becomes presence/absence of `DirectionalTOA()` in the output tuple, *and* a source-type concern (§2.9).

The rule eliminates type and physics gates. Not all branches.

---

<a id="2-step-1a"></a>
## 2. Step 1A — per-component `n_fourier_moments`

### 2.1 Convention: count, not order

Existing code uses `max_m` as a **count**: loops are `for m = 0:max_m - 1` in `rt_run.jl:208`, `rt_run_lin.jl:197`, `coxmunk_surface.jl:522`. The trait must match.

**Decision:** the trait is `n_fourier_moments` and returns a count. Loop convention stays `for m = 0:n - 1`.

| Component | `n_fourier_moments` | Why |
|---|---|---|
| Lambertian (any flavor) | 1 | Only `m = 0` contributes; no azimuthal structure |
| `RayleighScattering` | 3 | `m = 0, 1, 2` cover the Rayleigh phase matrix Greek coefs |
| `AerosolOptics` (greek_coefs.β of length L+1) | min(L+1, 2·Nstreams) | From Greek-coef tail, capped at solver resolution |
| `CoxMunkSurface`, `rpvSurfaceScalar`, `RossLiSurfaceScalar` | `2 · Nstreams` | VLIDORT cap (see §2.3) |
| `CanopySurface` | `2 · Nstreams` (provisional) | Proper LAD-driven derivation deferred |

### 2.2 `Nstreams` is *not* `Nquad`

`model.quad_points.Nquad` is *augmented* — `rt_set_streams.jl:40` redefines `Nquad = length(qp_μ)` after appending VZA/SZA nodes:

```julia
qp_μ = unique(FT[qp_μ[Nquad + 1:end]; cosd.(vza); μ₀])
n_eff = length(qp_μ) - length(wt_μ[Nquad + 1:end])
wt_μ = FT[wt_μ[Nquad + 1:end]; zeros(FT, n_eff)]
Nquad = length(qp_μ)
```

If `n_fourier_moments(::CoxMunk, ctx) = 2 * ctx.Nquad`, the count silently varies with the user's requested viewing angles. Different VZAs give different numbers of Fourier modes — not what VLIDORT's invariant says, not physically meaningful.

**Fix:** introduce a separate `Nstreams` capturing the polar discrete-ordinate count *before* VZA/SZA augmentation. Use it for the VLIDORT cap.

```julia
mutable struct QuadPoints{FT}
    # ... existing fields ...
    Nstreams::Int   # base solver streams, never augmented
    Nquad::Int      # may equal Nstreams (Gauss hemisphere) or include VZA/SZA (Radau)
end
```

Computed at quadrature setup time as `Nstreams = (l_trunc + 1) ÷ 2` *before* any VZA/SZA augmentation. Existing `Nquad` field stays for code needing the augmented value (kernel calls, Z-matrix dimensions).

### 2.3 The VLIDORT count formula

VLIDORT 2.8.3, `vbrdf_sup_masters.f90:2305`:

```fortran
NMOMENTS = 2 * NSTREAMS - 1
```

Followed by Fortran-style **inclusive** loop `DO M = 0, NMOMENTS` — runs `0..2·NSTREAMS - 1`, **`2·NSTREAMS` iterations**. Translated to vSmartMOM's count convention:

```
n_fourier_moments_BRDF = 2 · Nstreams
```

For standard `Nstreams = (l_trunc + 1) ÷ 2`:
- Even `l_trunc = 60`: `Nstreams = 30`, `2·Nstreams = 60`.
- Odd `l_trunc = 59`: `Nstreams = 30`, `2·Nstreams = 60`.

Empirical sanity check: run a Cox-Munk + Rayleigh test with `n_moments` set to 60 vs. 30 vs. 100 and confirm convergence to reference. *We'd value Sanghavi's view on the expected reference numerics before running this test* — see Q5.

**Rationale:** the BRDF Fourier series projects onto the discrete-ordinate basis. A polar grid of `Nstreams` cannot resolve azimuthal modes higher than `2·Nstreams - 1` regardless of the kernel. Past that aliases. Cap is set by *the solver's resolving power*, not surface complexity. Same formula for all VLIDORT BRDF kernels.

### 2.4 The two grids (kept separate)

| Quantity | Role | Source | Notes |
|---|---|---|---|
| `Nstreams` | Solver polar discrete-ordinate count, base | `rt_set_streams.jl` (new field) | Used for VLIDORT cap |
| `Nquad` | Augmented quadrature including VZA/SZA | `rt_set_streams.jl:40` | Used for kernel sizing |
| `nQuad_BRDF` | Azimuth integration grid for each Fourier coefficient | Hard-coded `100` in surface files | Quadrature accuracy, not physical limit |
| `n_fourier_moments` | Count of orders the loop runs over | New trait | The trait being added |

`nQuad_BRDF` recommended promotion to per-surface field (Step 1A.7, deferred — doesn't block).

### 2.5 Consumers reading the wrong value

Even with `max_m_bands` populated correctly in `model_from_parameters`, the loop driver in `rt_run.jl:97` calls `get_max_m(model)` which returns `model.solver.max_m` (the *scalar*), not the per-band vector. The `model.max_m` property forwarder at `types.jl:922` does give the per-band vector but it's not what the loop reads.

**Wiring change:** consumers in `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` must read the per-band value. Replace `max_m = get_max_m(model)` with `n_moments = n_fourier_moments(model, iBand)` (a new method on the trait, dispatching on model + band).

The scalar `model.solver.max_m` stays as a **user override** (hard cap, only consulted if explicitly set). Trait is the default; override clamps it down.

### 2.6 Proposed trait methods

```julia
# CoreRT/FourierDecomposition/n_fourier_moments.jl

"""
    n_fourier_moments(component, ctx) -> Int

Number of Fourier moments needed for an accurate RT solve given the context.
Loop convention is `for m = 0 : n_fourier_moments(...) - 1`.

The framework caps the per-band aggregate at `2·Nstreams` (VLIDORT's NMOMENTS+1
identity). Above this, modes alias on the polar discrete-ordinate basis.

`ctx.Nstreams` is the *base* polar stream count, NOT `model.quad_points.Nquad`
(which is augmented with VZA/SZA nodes after stream setup).
"""
n_fourier_moments(::RayleighScattering,        ctx) = 3

n_fourier_moments(::LambertianSurfaceScalar,   ctx) = 1
n_fourier_moments(::LambertianSurfaceSpectrum, ctx) = 1
n_fourier_moments(::LambertianSurfaceLegendre, ctx) = 1
n_fourier_moments(::LambertianSurfaceSpline,   ctx) = 1

# VLIDORT vbrdf_sup_masters.f90:2305 invariant
n_fourier_moments(::CoxMunkSurface,           ctx) = 2 * ctx.Nstreams
n_fourier_moments(::rpvSurfaceScalar,         ctx) = 2 * ctx.Nstreams
n_fourier_moments(::RossLiSurfaceScalar,      ctx) = 2 * ctx.Nstreams

# Provisional. Proper LAD-driven derivation lives in CanopyOptics.
n_fourier_moments(::CanopySurface,            ctx) = 2 * ctx.Nstreams

# From Greek-coef tail length, capped at 2·Nstreams.
n_fourier_moments(a::AerosolOptics,           ctx) =
    min(length(a.greek_coefs.β), 2 * ctx.Nstreams)
```

Aggregator:

```julia
function n_fourier_moments(model::RTModel, iBand::Int)
    ctx = (; Nstreams = model.quad_points.Nstreams,
             l_trunc = model.solver.l_trunc)
    n = n_fourier_moments(model.optics.rayleigh, ctx)
    for a in model.optics.aerosols.aerosol_optics[iBand]
        n = max(n, n_fourier_moments(a, ctx))
    end
    n = max(n, n_fourier_moments(model.surfaces[iBand], ctx))
    n = min(n, 2 * ctx.Nstreams)        # solver resolution cap
    return min(n, model.solver.max_m)    # user override
end
```

### 2.7 Knock-on simplifications

After the trait lands and consumers are wired:

- **`apply_correction!` argument list shortens.** Reads `n_fourier_moments(model, iBand)` from `ctx`, no positional `max_m`. **Reads the same value the loop used** — no drift between what was added and what gets subtracted in corrections.
- **Lambertian-only bands skip `m > 0`.** With `n_fourier_moments = 1`, the loop runs once. Lambertian surface kernel's `m == 0` guard inside `lambertian_surface.jl` stays intact (right guard for the right reason).
- **`vSmartMOM_Parameters.max_m` becomes a user-override knob**, not the primary path.

### 2.8 `GaussQuadFullSphere` retirement

Reading `rt_set_streams.jl:63-86`: `GaussQuadFullSphere` constructs `gausslegendre(2·Nquad)` over `[-1, 1]` and discards the lower hemisphere. The retained nodes/weights are *not* an optimal Gauss-Legendre quadrature for `[0, 1]` — they're half of an optimal full-sphere rule, which is a different polynomial fit. For RT problems where upper-hemisphere integration is what matters (essentially all plane-parallel RT), it's strictly worse than `GaussLegQuad`: same number of streams, lower accuracy, no compensating benefit.

**Retirement plan:**

1. **v2.0:** mark deprecated. `Base.depwarn` in the constructor: *"GaussQuadFullSphere is deprecated and will be removed in v2.1. Use GaussLegQuad — the mathematically correct hemispherical Gauss-Legendre quadrature."* YAML loader (`IO/Parameters.jl:99`) still accepts the string but emits the same warning.
2. **Audit reference benchmarks.** Any test currently using `GaussQuadFullSphere` gets regenerated with `GaussLegQuad`. Diff documented in changelog.
3. **v2.1:** remove the type, the `rt_set_streams` method, and the YAML mapping.

This fits naturally into the dispatch architecture — `AbstractQuadratureType` is already a strategy family; removing one concrete subtype is mechanical. *Q3 confirms with Sanghavi.*

### 2.9 Source-type forward-compatibility (deferred refactor)

Today `DNI` and `SFI` types exist in `CoreRT/types.jl:123,126` as `AbstractSourceType` subtypes but are essentially vestigial — the `SFI::Bool` parameter does the real work. Eventually the source axis needs proper dispatch:

- **`SolarSource`** — direct solar beam (current behavior)
- **`ThermalSource`** — Planck emission per layer (when thermal RT lands)
- **`LunarSource`** — for nighttime atmospheric measurements
- **`CompositeSource`** — solar + thermal mixed (real for SWIR bands)

Each has fundamentally different RTE physics (thermal has `B(T)` source distributed throughout; solar has `δ(μ - μ₀)` boundary condition). This is a future-proofing concern for thermal/MW work.

**We don't refactor source types in this round.** But Step 1's design must not foreclose the future refactor. The forward-compat invariant:

> *Any code touching solar geometry or source-term assembly in Step 1 must funnel through `ctx.source`, not through hardcoded `μ₀`/`I₀` references or positional arguments.*

Concretely:

1. `RTContext` carries a `source::AbstractSourceType` field, parameterized.
2. Today populated with `DefaultSolarSource{FT}` carrying current `μ₀, I₀`.
3. Strategy methods (SS corrections, etc.) read solar geometry through `ctx.source` only.
4. **Rename** `DNI` → `DefaultSolarSource` and `SFI` → `LegacySFISource` (or similar — exact names TBD with Sanghavi). Old YAML strings deprecated with warnings.
5. `SFI::Bool` parameter retained unchanged for now — still does the existing work.

When thermal lands later, the path is: add `ThermalSource` type, add `source_term(::ThermalSource, ctx)` methods, retire `SFI::Bool`. No re-architecture needed because Step 1's strategy methods already go through `ctx.source`.

The renames in (4) are the only Step 1 change touching user-visible YAML keywords. Old strings (`"DNI"`, `"SFI"`) emit deprecation warnings but still work; v2.1 removes the aliases.

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
3. There is an **"external" hook**.

### 3.3 Strategy hierarchy

```julia
# CoreRT/FourierDecomposition/ss_correction.jl

abstract type AbstractSSCorrection end

"""Default no-op."""
struct NoSSCorrection <: AbstractSSCorrection end

"""
    BRDFSpecularCorrection{S<:AbstractSurfaceType}

Replaces the truncated Fourier reconstruction of the surface BRDF SS
contribution with its exact evaluation. Nakajima-Tanaka TMS for surfaces.
"""
struct BRDFSpecularCorrection{S<:AbstractSurfaceType} <: AbstractSSCorrection
    surface::S
    nQuad_BRDF::Int
end

"""
    AerosolForwardPeakCorrection{FT}

The δ-M companion: replaces the truncated SS aerosol contribution with the
exact phase function evaluation. Type defined in Step 1B but **not enabled in
defaults until Step 1B.7 implements the apply method**.
"""
struct AerosolForwardPeakCorrection{FT} <: AbstractSSCorrection
    fᵗ_threshold::FT
end

"""
    ExternalSSCorrection{F}

User-supplied correction. Equivalent to VLIDORT's `DO_FOCORR_EXTERNAL`.
"""
struct ExternalSSCorrection{F} <: AbstractSSCorrection
    callback::F
end
```

Per-surface trait controlling defaults:

```julia
has_specular_peak(::AbstractSurfaceType) = false
has_specular_peak(::CoxMunkSurface)      = true
```

Default selection:

```julia
function default_ss_corrections(surfaces, aerosol_optics)
    map(eachindex(aerosol_optics)) do b
        cs = ()
        if has_specular_peak(surfaces[b])
            cs = (cs..., BRDFSpecularCorrection(surfaces[b], 100))
        end
        # AerosolForwardPeakCorrection NOT included until Step 1B.7 implements
        # the apply method. Type defined ≠ ready to enable.
        isempty(cs) ? (NoSSCorrection(),) : cs
    end
end
```

### 3.4 The unified call site

`rt_run.jl:310-315` becomes:

```julia
SFI && @timeit "SS corrections" foreach(c -> apply_correction!(c, ctx),
                                        model.solver.ss_corrections[iBand])
```

For the lin run: same call, but `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::LinContext)` is **a separate method** that updates both `R_SFI` and `Ṙ_SFI` consistently. Lin SS correction is a real implementation work item, not a wiring one — only the *call site* is one line; the dispatched method has its own derivative math.

For multisensor: similar, with per-sensor-level attenuation. Also a real method.

### 3.5 What this gives the linearization redesign

Each correction is its own self-contained unit of differentiation. Wrap `apply_correction!` in a ForwardDiff seed → get `∂(correction)/∂(Mie params, surface params)` → add to chain rule via `lin_added_layer_all_params.jl`. Current implementation makes that hard because the correction lives *outside* the lin seam. With strategy structs, the hack becomes a typed unit the linearization framework can see.

### 3.6 Caching deferred

Current Cox-Munk `apply_ss_correction!` re-evaluates Fourier coefficients in `_fourier_coeff_element` rather than caching from the main loop. Wasteful but safe. Recommendation: keep recomputing for the first landing.

---

<a id="4-step-1c"></a>
## 4. Step 1C — `AbstractHotSpot` field on `CanopySurface` (wiring deferred to Step 2)

### 4.1 What CanopyOptics already ships

`CanopyOptics.jl/src/canopy_structure/hotspot.jl`:

- `abstract type AbstractHotSpot{FT<:Real}`
- `struct NoHotSpot{FT}` — independent gaps, `C_hs = 1`
- `struct KuuskHotSpot{FT}` with size parameter `h`
- `joint_gap_probability(model, k_s, k_o, μ_s, μ_o, dϕ, L)` — `exp(-(k_s+k_o)L) · C_hs(L)`
- ForwardDiff hooks (`_hotspot_value(x::ForwardDiff.Dual)`) for future linearization

The author was thinking ahead about linearization. Kuusk math is closed-form; Jacobians w.r.t. `h` (and through `k_s, k_o` to LAI/G) drop out under ForwardDiff.

**Abstraction shipped. Implementation correct. Linearization story anticipated. Missing: vSmartMOM-side wiring.**

### 4.2 Why this is *not* a subtype of `AbstractSSCorrection`

The temptation: add `KuuskHotSpotCorrection <: AbstractSSCorrection`. It would be wrong.

| | `AbstractSSCorrection` | `AbstractHotSpot` |
|---|---|---|
| **When** | After Fourier loop completes | Inside the canopy direct-beam pipeline |
| **How** | Additive: `(M_exact − M_fourier)` per geometry | Multiplicative: `C_hs · P_so/P_o` on joint gap probability |
| **What it patches** | Fourier truncation losing specular-peak resolution | Beer-Lambert assuming independent gaps |

If `AbstractHotSpot` were forced into `AbstractSSCorrection`, the apply hook would either widen to handle in-loop calls (abstraction leaks) or operate on `R_SFI` after the fact and reverse-engineer source assembly (hack).

Conceptual similarity captured by **shared dispatch idiom**, not shared type.

### 4.3 The wiring is harder than first thought

v0.2 sketched "one multiplicative correction in `get_canopy_elem_rt_SFI!`" — that turned out to be wrong. The SFI kernel at `elemental_canopy.jl:129` sees only quadrature streams and Fourier mode, not view geometry (`vza, vaz, dϕ`). Geometry isn't materialized until `postprocessing_vza!` after the Fourier loop.

**Physical reasoning** (the real determinant): the hotspot is a correlated-path-pair phenomenon. `P_so(L)` describes the joint probability that the *same gap* is open along both sun and view paths. When sun and view coincide (`α → 0`), every gap that lets the photon in also lets it out — paths are perfectly correlated, gaps don't multiply, and `P_so → exp(-min(k_s, k_o)L)` rather than the independent-paths product.

**This correlation only exists for single-scattering events observed directly without further interaction.** Once a photon scatters once into the diffuse field, its outgoing direction decorrelates from the incoming direction; subsequent gap probabilities are again independent.

So the hotspot correctly applies to:
- Direct-solar → direct-view contribution (single scattering observed straight back)

It does **not** apply to:
- Direct-solar → diffuse field → view contribution (multiple scattering)
- Diffuse-source → diffuse field → view contribution

A per-Fourier-mode hotspot expansion would over-correct: it would impose path-correlation on contributions where photons have already decorrelated. This is what VLIDORT's FO subsystem effectively does for canopies (single-scatter pass, separate from MS solve), and what the canonical canopy RT codes (Verstraete-Pinty-Walthall, Kuusk's papers, SAIL-with-hotspot) do.

### 4.4 Step 1 does the field, Step 2 does the wiring

**In Step 1C (this round):**

1. Add `hotspot::AbstractHotSpot{FT}` field to `CanopySurface` in `CoreRT/types.jl`. Default `NoHotSpot{FT}()`.
2. Plumb through `model_from_parameters` and YAML config so `CanopySurface(..., hotspot=KuuskHotSpot(0.1))` is constructable.
3. Tests verifying the field is constructable but the runtime behavior is bit-for-bit identical to today (since `NoHotSpot` is the default and no code path acts on it yet).
4. Document explicitly that wiring is deferred to Step 2.

That's ~30 lines. Honest scope.

**In Step 2 (next round):**

1. Implement direct-beam vs. diffuse separation in postprocessing. The current `j₀±` source terms in `elemental_canopy!` mix both contributions. Cleanest factoring: compute direct-beam analytically as a closed-form `exp(-k_o · L_cumulative) × Z(μ_s → μ_v)` that bypasses the Fourier loop entirely, matching FO's structure.
2. Apply hotspot multiplicatively in postprocessing on the directional direct-beam term only.
3. Validate against RAMI HOM26/HOM27.

⚠️ **For Sanghavi (Q2):** any current canopy benchmark tuned to match the no-hotspot result needs to be flagged before Step 2 lands — the fix is silent unless we explicitly call it out in the changelog. Specifically asking: are there canopy reference numbers in your test suite that would change?

---

<a id="5-step-1d"></a>
## 5. Step 1D — optical-properties algebra: docs, dimensioned zero, contributor pattern

### 5.1 The motivation

`constructCoreOpticalProperties` in `compEffectiveLayerProperties.jl:11-65` hardcodes Rayleigh + aerosols + absorption as the atmospheric optical contributors:

```julia
combo = rayl
for i = 1:nAero
    aer = createAero(...)
    combo = combo .+ aer
end
combo = combo .+ [CoreAbsorptionOpticalProperties(...) for ...]
```

Two consequences:

1. **Adding a new contributor (clouds, ice crystals, sea-salt with custom phase function) requires editing this file.** The Rayleigh-vs-Cabannes branch (`_rayleigh_greek_source`) is also a special case here. Every new physics adds a special case.
2. **The thermal/MW case has no Rayleigh.** At 10 μm, Rayleigh optical depth is ~10⁻⁶. The current code assumes Rayleigh is always present and uses it as the loop initializer. There's no way to express "no Rayleigh" cleanly.

The operator algebra (`+` and `*`) on `CoreScatteringOpticalProperties` is already doing most of the right work — the contributors *combine via `+`*. What's missing is making the contributor *list* the dispatched axis rather than hardcoded sequence.

### 5.2 The existing algebra (what the docs need to say)

Five overloads in `types.jl:1083-1137`:

- `+(scattering, scattering)`: τ-and-ϖ-weighted average of phase functions, with the `wx == 0` short-circuit corner case (line 1099)
- `*(scattering, scattering)`: vertical/band concatenation via `cat` along dim 3
- `+(scattering, absorption)` and `+(absorption, scattering)`: thicken layer, leave Z untouched (commutative)
- `*(scalar, scattering)`: scale τ only, leave ϖ and Z untouched

What's *not* explicitly documented but is real:

- **Storage-form duality.** `CoreScatteringOpticalProperties` may carry Z as `(n_pq, n_pq, 1)` (compact, broadcast across layers) or `(n_pq, n_pq, n_layers)` (expanded, per-layer). All operators accept both via broadcasting; `expandOpticalProperties` converts compact → expanded explicitly when kernels need it.
- **Closure under combination.** `+` of compact and expanded forms broadcasts correctly. `*` calls `expandOpticalProperties` on both operands first.
- **Algebraic properties.** `+` is commutative and associative on scattering+scattering, commutative on scattering+absorption (no identity element today). `*` is associative but **not** commutative (vertical/band stacking — order matters).
- **Missing absorption + absorption.** Two `CoreAbsorptionOpticalProperties` cannot currently be added. Comes up the moment per-species absorption contributors exist (e.g. separate O₂, H₂O, CO₂).

**Layer-1 work for Step 1D: write the docs.** A module-level docstring formalizing the algebra. Maybe ~150 lines of docstrings, no code change beyond the absorption+absorption overload. This pays off immediately for current contributors and is the precondition for the contributor pattern.

### 5.3 Dimensioned zero (no new sentinel type)

To start the contributor loop with an empty `props` (so thermal/MW with no Rayleigh works), the algebra needs a zero element.

**Decision (per Q4):** use a dimensioned zero — a `CoreScatteringOpticalProperties` value with τ = 0. Reasons:

1. **The existing `+` short-circuit `all(wx .== 0.0)` (line 1099) already implements zero-element behavior correctly.** When `x.τ = 0`, `wx = 0`, the short-circuit returns y's data. The current code *already implements* the zero-element semantics — it just doesn't have a name for "the value that triggers this short-circuit."
2. **No new type means no new overloads.** A sentinel `EmptyOpticalProperties` would require ~8-10 new method signatures (`Empty + Scattering`, `Scattering + Empty`, ...). The dimensioned zero needs zero new overloads.
3. **The cost is small.** Z is stored *compact* (shape `(n_pq, n_pq, 1)`) for the zero, not expanded. For `n_pq = 128`, that's ~256 KB per band per Fourier moment, not 256 MB. The existing code already produces compact Z from contributors and only expands when needed.
4. **Restrict zero to `+` only.** The contributor loop uses `+`; band concatenation via `*` happens outside the loop, between always-real operands. The zero never reaches `*`, avoiding the "phantom zero band" semantic question.

```julia
function zero_optical_properties(ctx)
    n_pq = ctx.pol_type.n * ctx.Nstreams
    CoreScatteringOpticalProperties(
        τ   = zeros(ctx.FT, ctx.n_layers),
        ϖ   = zero(ctx.FT),                       # placeholder; never read because τ=0
        Z⁺⁺ = zeros(ctx.FT, n_pq, n_pq, 1),       # compact form
        Z⁻⁺ = zeros(ctx.FT, n_pq, n_pq, 1),
    )
end
```

**Pure-absorption case** (no scattering contributors at all — plausible for thermal IR with only absorption, no clouds): the algebra produces `CoreAbsorptionOpticalProperties` (since `zero + abs + abs = abs` via the existing absorption+absorption overload, once added). That's a valid output type for an absorption-only solver, just not for the scattering RT. **Step 1D adds an explicit guard:** if the final `props isa CoreAbsorptionOpticalProperties` after the contributor loop, error before entering the scattering RT solver with a clear message. Future work routes that case to a different solver entry point.

### 5.4 The contributor pattern

Each contributor is a small concrete struct that knows how to convert itself to optical properties for a given context:

```julia
abstract type AbstractOpticalContributor end

# Single trait: each contributor implements optical_properties(c, ctx)
# returning an AbstractOpticalProperties value. The construction loop
# uses the existing + algebra.

struct RayleighContributor{FT, S} <: AbstractOpticalContributor
    greek_source::S         # carries Rayleigh-vs-Cabannes choice as data
    τ_rayl::Vector{FT}
    ϖ::FT
end

struct AerosolContributor{FT, O} <: AbstractOpticalContributor
    aerosol_optics::O
    τ_aer::Vector{FT}
    apply_delta_m::Bool
end

struct AbsorptionContributor{FT} <: AbstractOpticalContributor
    τ_abs::Vector{FT}
end

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

The construction loop:

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

Per-band concatenation via `*` happens at the outer level, unchanged:

```julia
layer_opt = [prod([band_props[i][iz] for i in iBand]) for iz in 1:n_layers]
```

### 5.5 What this gives

- **Rayleigh-vs-Cabannes branch eliminated.** `RayleighContributor` *carries* its `greek_source` as a field, populated correctly at model construction time. The construction loop never sees the branch.
- **δ-M decision is per-contributor data.** A future contributor (clouds, ice) that wants different truncation can decline by setting `apply_delta_m = false` or use a different scheme.
- **Thermal/MW works for free.** No Rayleigh contributor in the list → loop starts from zero → produces an absorption-only or clouds-only result. The pure-absorption guard catches the "no scattering anywhere" edge case.
- **Per-species absorption becomes natural.** Today's aggregate `AbsorptionContributor` can be split into `O2AbsorptionContributor`, `H2OAbsorptionContributor`, etc. — each contributes via `+`, and the new `+(absorption, absorption)` overload composes them. Enables future per-species Jacobians without architectural change.
- **Adding cloud water optics** is 2 steps: define `CloudWaterContributor`, define `optical_properties(::CloudWaterContributor, ctx)`. Done.

### 5.6 Open design choice: `fScattRayleigh` extraction

The current code computes `fScattRayleigh = [Array(rayl[i].τ ./ combo[i].τ) ...]` as a side effect — a Rayleigh-specific quantity used downstream by `_expand_layer_rayleigh!` for inelastic Raman. Two options:

- **(A) Trait-based.** `f_scattering(c::AbstractOpticalContributor, props, ctx) -> Vector{FT}` defined per contributor. Generic; lets future contributors define their own scattering-fraction diagnostics.
- **(B) Push to use site.** Compute `fScattRayleigh` inside `_expand_layer_rayleigh!` (the only consumer) directly from the Rayleigh contributor and the assembled props.

We recommend **(B)** for simplicity — the only consumer that needs it computes it locally; no scaffolding for hypothetical future cases. *Open for review (Q5).*

### 5.7 Operator style preserved

The existing `+` and `*` operators stay exactly as written today. The contributor pattern adds *one new abstract type* (`AbstractOpticalContributor`) and *one new trait* (`optical_properties`). The combination semantics are the existing operator algebra. New contributors just need to produce `CoreScatteringOpticalProperties` or `CoreAbsorptionOpticalProperties` values; everything else falls out of `+`.

---

<a id="6-step-2"></a>
## 6. Step 2 — output products and RT modes as dispatched strategies

This is the larger refactor — 1004 lines across `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`. Lands after Step 1 has established the dispatch idiom on small contained pieces.

### 6.1 The leak

`rt_run.jl:136-140` unconditionally allocates `R`, `T`, `R_SFI`, `T_SFI`, `ieR_SFI`, `ieT_SFI`, `hdr`, `bhr_uw`, `bhr_dw`. Return tuple at line 325 is conditional on `SFI`, but the allocations aren't:

```julia
return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw[1,:], bhr_dw[1,:]) : (R, T)
```

Three problems:

1. Every run pays allocation cost for every output, used or not.
2. Adding a new output requires touching every `rt_run*.jl` file (1004 lines for what is structurally one loop).
3. The choice is implicit in `if SFI`, `if RAMI` scattered through postprocessing — not a workflow boundary the user controls explicitly.

### 6.2 Two orthogonal axes: mode and outputs

The linearization redesign isn't just "another output." It changes:
- Model inputs (lin parameters)
- Layer construction (`AddedLayerLin` with tangent fields)
- Composite storage (derivative state)
- Kernel calls (`elemental_lin!`, `doubling_lin!`, `interaction_lin!`)
- The Fourier loop body itself

That's a **mode**, not an output. Two abstract types:

```julia
abstract type AbstractRTMode end
struct ForwardMode <: AbstractRTMode end
struct LinMode{PL} <: AbstractRTMode
    parameter_layout::PL
end

abstract type AbstractRTOutput end
struct DirectionalTOA <: AbstractRTOutput end
struct DirectionalBOA <: AbstractRTOutput end
struct DirectionalDirectBeam <: AbstractRTOutput end   # enables canopy hotspot wiring
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

A `LinMode` run can record `DirectionalTOA` *and* `JacobianTOA`. A `ForwardMode` run can only record `DirectionalTOA`. The two axes are independent dimensions of "what kind of run is this."

### 6.3 The new public API

```julia
run_rt(model; mode=ForwardMode(), outputs=(DirectionalTOA(), DirectionalBOA()), iBand=1)
```

- Defaults to forward mode + standard directional outputs (matches current `rt_run(model)` behavior).
- `rt_run`, `rt_run_lin`, `rt_run_multisensor` become thin compatibility wrappers calling `run_rt` with the appropriate mode/outputs.
- Deprecation cycle: v2.0 introduces `run_rt`; existing entry points keep working as wrappers; v2.x removes the wrappers (date TBD with downstream users).

### 6.4 What this collapses

| Current entry point | Becomes |
|---|---|
| `rt_run(model)` | `run_rt(model, outputs=(DirectionalTOA(), DirectionalBOA()))` |
| `rt_run` with `SFI=true` and RAMI | `run_rt(model, outputs=(DirectionalTOA(), HDRFOutput(), BHROutput()))` |
| `rt_run_multisensor.jl` (177 lines) | `run_rt(model, outputs=(MultiSensorOutput(levels),))` |
| `rt_run_lin.jl` (303 lines) | `run_rt(model, mode=LinMode(layout), outputs=(DirectionalTOA(), JacobianTOA()))` |

The 1004 lines collapse into one coordinator (~150 lines) plus per-output strategy methods. **The lin variant disappearing as a separate file is the deduplication that unlocks the linearization redesign.**

### 6.5 `DirectionalDirectBeam` enables canopy hotspot

Step 2 introduces `DirectionalDirectBeam` as an output strategy that records the direct-solar-to-direct-view single-scatter contribution separately from diffuse. Once it exists, the canopy hotspot wiring (deferred from Step 1C) is straightforward: the hotspot multiplicatively corrects the `DirectionalDirectBeam` accumulator only, leaving diffuse contributions untouched.

This is a real architectural reason to do hotspot in Step 2 rather than Step 1: the direct-beam separation is generally useful (it's how FO works in VLIDORT), not specific to canopy.

### 6.6 Migration substeps

Cannot land in one commit. Order:

1. Define `AbstractRTOutput`, `AbstractRTMode`, `RTContext`. No call sites changed.
2. Implement `DirectionalTOA`, `DirectionalBOA`. Both dispatch to existing `postprocessing_vza!` internals. No behavior change.
3. Refactor current `rt_run` to take `outputs::Tuple` argument with defaults matching today. Existing callers unaffected. Introduce `run_rt` as the new entry point.
4. Implement `HDRFOutput`, `BHROutput`. Migrate RAMI workflows. Validate against existing RAMI benchmarks.
5. Implement `MultiSensorOutput`. Wrapper makes `rt_run_multisensor` call `run_rt`.
6. Implement `DirectionalDirectBeam`. Wire canopy hotspot multiplicatively. Validate against RAMI HOM26/HOM27.
7. Implement `LinMode` + `JacobianTOA`. Wrapper makes `rt_run_lin` call `run_rt`. **Riskiest substep** — most lin methods to move. Validate against `LINEARIZATION_BUGS.md` test set.

Each substep is independently shippable; benchmarks at each boundary.

---

<a id="7-extension-cookbook"></a>
## 7. Extension cookbook (aspirational)

Recipes for common future changes, **describing what the API will look like after Step 1 + Step 2 land**, not what's currently required. Aspirational by design — they're a target for the architecture and a contributor reference.

### Recipe: adding a new BRDF surface kernel

```julia
# 1. Define the type
struct MyBRDF{FT} <: AbstractSurfaceType
    # ... your parameters ...
end

# 2. Implement the existing reflectance interface (already in vSmartMOM)
function reflectance(brdf::MyBRDF, n_stokes, μᵢ, μⱼ, ϕ; ...)
    # ... your physics ...
end

# 3. Declare the Fourier-moment count
n_fourier_moments(::MyBRDF, ctx) = 2 * ctx.Nstreams  # or whatever your kernel needs

# 4. Declare specular peak presence (true if you want SS correction by default)
has_specular_peak(::MyBRDF) = true   # or false
```

That's it. No edits to `rt_run.jl`, `model_from_parameters.jl`, or any *_lin variants.

### Recipe: adding a new optical contributor (e.g. cloud water)

```julia
# 1. Define the type
struct CloudWaterContributor{FT, M} <: AbstractOpticalContributor
    cloud_optics::M
    τ_cloud::Vector{FT}
    apply_delta_m::Bool
end

# 2. Implement optical_properties (returns either CoreScatteringOpticalProperties
#    or CoreAbsorptionOpticalProperties)
function optical_properties(c::CloudWaterContributor, ctx)
    # ... your conversion ...
end
```

Then user code:

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
# 1. Define the type
struct MyCorrection <: AbstractSSCorrection
    # ... parameters ...
end

# 2. Implement apply_correction!
function apply_correction!(c::MyCorrection, ctx)
    # ... mutate ctx.R_SFI in place ...
end

# 3. Optionally implement default selection
needs_correction(c::MyCorrection, ctx) = ...
```

### Recipe: adding a new output product

```julia
# 1. Define the type and what it accumulates
struct PhotolysisRate <: AbstractRTOutput end

struct PhotolysisState
    # ... whatever state you need ...
end

# 2. Implement the (allocate, accumulate!, finalize) triple
allocate_output(::PhotolysisRate, model, nSpec) = PhotolysisState(...)
accumulate_output!(s::PhotolysisState, ::PhotolysisRate, ctx) = ...
finalize_output(s::PhotolysisState, ::PhotolysisRate, ctx) = ...
```

User code:

```julia
results = run_rt(model, outputs=(DirectionalTOA(), PhotolysisRate()))
```

### Recipe: adding a new validation benchmark

(Aspirational target — see §9 / C3 for the validation-as-data future direction.)

```toml
# benchmarks/cox_munk_glint_high_wind.toml
[atmosphere]
profile = "us_standard"
contributors = ["rayleigh", "aerosol_dust"]

[surface]
type = "cox_munk"
wind_speed = 15.0

[geometry]
sza = 30.0
vza = [0, 30, 60]
vaz = [0, 90, 180]

[expected]
source = "vlidort_2.8.3"
file = "reference/cox_munk_15ms.h5"
tolerance = 1e-5
```

One TOML file per benchmark. Run via `julia --project test/run_benchmarks.jl benchmarks/`. New benchmark = new file.

---

<a id="8-type-gates"></a>
## 8. Remaining type-gates inventory (migration debt)

Beyond what's addressed in Step 1 and Step 2, these are the remaining places where `isa` checks or boolean-as-physics-gate patterns need future migration. Listed here so they don't accumulate further during the migration.

| Smell | Files | Replace with | Status |
|---|---|---|---|
| `if brdf isa CanopySurface` (4 occurrences) | `rt_run.jl:170, 175, 182, 245` | `requires_canopy_setup(::AbstractSurfaceType)` trait + dispatched setup methods | Step 2 |
| `if RS_type isa AbstractRamanType` style branches | `interaction_inelastic.jl`, `elemental_inelastic.jl` | Already partially dispatched — finish it | Future |
| `arch isa GPU` checks | scattered (`array_type`, `synchronize_if_gpu`) | Already mostly dispatched on `AbstractArchitecture`; audit holdouts | Future |
| `if SFI` boolean threading | every `rt_run*.jl`, `postprocessing_*.jl` | Goes away under Step 2's output strategies + source-type axis | Step 2 + future |
| Polarization-type branching via `pol_type.n` | many places | Numerical property — likely fine, but worth auditing | Audit only |
| `Stokes_IQU` implementation completeness | scattered | Implementation audit (IQU code paths most likely to have rot) | Future |
| `AbstractSourceType` resurrection | `types.jl:120-126` | New `SolarSource`, `ThermalSource`, etc. | Future (thermal RT trigger) |

The forward-compat invariant from §2.9 (all solar geometry through `ctx.source`) keeps the `AbstractSourceType` resurrection mechanical when it happens.

---

<a id="9-future-architecture"></a>
## 9. Future architecture (deferred)

These are flagged for the post-v2.0 roadmap. Each is the natural next application of the dispatch rule but is out of scope for this round.

**C1 — Optical-properties algebra extension** (now Step 1D, see §5). *In scope this round.* Listed here for completeness.

**C2 — Layer-kernel interface stabilization.** Doubling/elemental/interaction kernels are coupled to specific layer-storage types (`AddedLayer`, `AddedLayerRS`, `AddedLayerLin`). Each new layer flavor (RRS, VRS, future inelastic variants) requires a new kernel variant. A more future-proof structure: kernels operate on an `AbstractLayerStorage` interface (get/set primitives, no concrete-type assumptions); storage type is a parameter. Belongs in the linearization redesign track.

**C3 — Validation-as-data.** Reference benchmarks (6SV1, Natraj, RAMI) currently hand-written in Julia files. Future-proof setup: each benchmark is a TOML/YAML data file; one generic test driver runs all of them. New benchmark = new file. Right granularity for collaborative scientific software. The 6SV1/Natraj benchmark restoration (mentioned in the project memory as deleted from `unified-vsmartmom`) is the moment to fix this, not just restore.

**C4 — Physics audit traits.** A class of future-proofing specific to RT: traits letting runtime introspection answer *"does my model handle polarization correctly?"*, *"is δ-M truncation active?"*, *"what's the assumed BRDF treatment?"*. A `physics_summary(model)` returning a structured description of every active assumption, dispatched through every component. Useful for reproducibility (saved with results), sanity checks, downstream tooling. Cheap to add given dispatch architecture, expensive to retrofit.

**C5 — Versioned interface contracts.** When BRDF / aerosol / surface implementers extend the package externally — which the Mie-decomposition rename and BRDF dispatch arguably enable — they need to know which methods are required vs. optional, with a stable contract across vSmartMOM minor versions. Julia idiom: documented interface methods with `Required`/`Optional` annotations and a `check_interface(MyBRDF)` test helper. Compare to `AbstractArrays` documentation. Matters most if vSmartMOM grows an ecosystem.

---

<a id="10-execution-plan"></a>
## 10. Execution plan

### Step 1 — Fourier work and contained refactors (this round)

| # | What | Risk | Validation |
|---|---|---|---|
| 1A.0 | Rename `AbstractFourierDecompositionType` → `AbstractMieDecompositionAlgorithm`. Add deprecation alias. | ~50–100 lines including tests, docs, IO imports. | All existing tests pass. |
| 1A.0.5 | Deprecate `GaussQuadFullSphere`. `Base.depwarn` in constructor; YAML loader warns + still works. | Touches user-visible YAML keyword. | Audit benchmarks for `FullSphere` use. |
| 1A.1 | Introduce `FourierDecomposition` module skeleton. | Plumbing. | Module loads. |
| 1A.2 | Add `RTContext` struct (§1.3). Include `source::AbstractSourceType` field (§2.9). | New abstraction. | Trivial unit tests. |
| 1A.3 | Rename `DNI` → `DefaultSolarSource`, `SFI` → `LegacySFISource` (§2.9). YAML strings deprecated, still work. | Touches user-visible names. | Existing tests pass; deprecation warnings appear. |
| 1A.4 | Add `Nstreams` field to `QuadPoints` separate from augmented `Nquad`. | One new field, two value-producers. | All quadrature setups produce both values. |
| 1A.5 | Add `n_fourier_moments` trait + methods (§2.6). | New code. | Unit tests for each surface type, count convention asserted. |
| 1A.6 | **Wire `n_fourier_moments(model, iBand)` into the consumers.** Replace `get_max_m(model)` with the trait call in `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`. | **Behavior change risk.** Cox-Munk bands may now compute more moments than before. | 6SV1 / Natraj reference benchmarks before and after. |
| 1A.7 | (Deferred) Promote `nQuad_BRDF` to per-surface field. | Trivial. | Convergence sweeps. |
| 1B.1 | Refactor BRDF Fourier integration into shared helper (`_brdf_fourier_coefficient` in `FourierDecomposition`). | Pure refactor. | Existing tests pass bit-for-bit. |
| 1B.2 | Add `AbstractSSCorrection` hierarchy (§3.3). Migrate Cox-Munk `apply_ss_correction!` to `BRDFSpecularCorrection{<:CoxMunkSurface}`. **Read `n_moments` from ctx**, not via recomputation. | Pure relocation. | Existing Cox-Munk tests pass bit-for-bit. |
| 1B.3 | Wire `apply_correction!` into `rt_run.jl`. Add `default_ss_corrections` (without `AerosolForwardPeakCorrection`). | Pure call-site swap. | Existing benchmarks unchanged. |
| 1B.4 | Add lin-context `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::LinContext)` updating `R_SFI` and `Ṙ_SFI` consistently. | New physics in lin path. | Validate Jacobians against finite-difference at glint geometry. |
| 1B.5 | Add multisensor-context method, per-sensor-level attenuation. | New physics in multisensor path. | Validate against multisensor reference cases if any. |
| 1B.6 | (Deferred to later round) Implement `AerosolForwardPeakCorrection`. | New physics. | PyVLIDORT high-AOD strongly-forward-peaked aerosol cases. |
| 1B.7 | (Deferred) Cache Fourier coefficients for SS correction reuse. | Optimization only. | No behavior change. |
| 1C.1 | Add `hotspot::AbstractHotSpot{FT}` field to `CanopySurface`, default `NoHotSpot{FT}()`. Plumb through `model_from_parameters` and YAML. | Pure plumbing; field is constructable but no-op. | All existing tests pass bit-for-bit. |
| 1D.1 | Document optical-properties algebra (§5.2). Add absorption + absorption overload. | Mostly docs + 1 trivial overload. | All existing tests pass. |
| 1D.2 | Add dimensioned `zero_optical_properties` constructor (§5.3). | New constructor. | Unit tests for zero-element behavior. |
| 1D.3 | Add `AbstractOpticalContributor` + `optical_properties` trait (§5.4). | New abstraction. | Unit tests for each contributor type. |
| 1D.4 | Migrate `constructCoreOpticalProperties` to use the contributor pattern. | **Touches the optical-properties pipeline.** | Bit-for-bit numerical equivalence with current default contributors (Rayleigh + aerosols + absorption). |

**Step 1 boundaries for benchmarks:**
- After 1A.0–1A.5: all existing tests pass bit-for-bit. Pure plumbing.
- After 1A.6: 6SV1 / Natraj benchmarks for any band where `n_moments` could change.
- After 1B.3: existing Cox-Munk TOA reference cases.
- After 1B.4: lin Jacobian regression at glint geometries.
- After 1B.5: multisensor + Cox-Munk reference cases.
- After 1D.4: end-to-end equivalence test on a representative atmosphere.

### Step 2 — Output products as strategies (next round)

Per §6.6. Seven-substep migration, each independently shippable, benchmark boundary at each. Order: `Directional*` → `HDRF/BHR` → `MultiSensor` → `DirectionalDirectBeam` (enables canopy hotspot wiring) → `LinMode/Linearized`. Lin migration is the riskiest substep, lands last.

Step 2 starts after Step 1 lands. Step 1B.6 (`AerosolForwardPeakCorrection`) and Step 2 can proceed in parallel — they touch different code paths.

---

<a id="11-open-questions"></a>
## 11. Open questions

In rough order of "blocks the next commit" → "nice to settle eventually":

1. **Q1 (preamble)** — dispatch rule as governing principle. Sanghavi's view on the philosophical commitment.
2. **Q2 (preamble)** — canopy hotspot regression risk. Are there reference benchmarks tuned to the no-hotspot result? **Blocks step 1C.1 changelog work and step 2.6.**
3. **Q3 (preamble)** — `GaussQuadFullSphere` retirement. Does Sanghavi concur, or is there a use case to preserve? **Blocks step 1A.0.5.**
4. **Q4 (preamble)** — dimensioned zero vs. sentinel for the optical-properties algebra. Blocks step 1D.2.
5. **Q5 (preamble)** — `n_fourier_moments(::CoxMunk) = 2 · Nstreams` derivation. Empirical sanity check ordering. **Blocks step 1A.5 validation plan.**
6. **`l_trunc` semantics** — does vSmartMOM's `l_trunc` correspond to VLIDORT's "highest Legendre order in the truncated phase function expansion" or to VLIDORT's `NSTREAMS`? Easiest answer: derive both formulas from `Nstreams = (l_trunc+1)÷2` and check which matches current numerics. Blocks step 1A.5 if Q5 surfaces ambiguity.
7. **`fScattRayleigh` extraction** — trait-based (Option A) or push-to-use-site (Option B, recommended). §5.6.
8. **`ExternalSSCorrection` callback signature** — mutation-via-context vs. return-value API. Affects downstream tooling.
9. **`AerosolForwardPeakCorrection` threshold** — default `fᵗ > 0.01`. Reasonable starting point? Worth tuning during validation.
10. **`RTContext` lifecycle** — exact field set and m-mutation contract. §1.3.
11. **Naming harmonization** — existing `AbstractTruncationType`/`NoTruncation` vs. proposed `Abstract*Correction`/`No*Correction`. Harmonize during this work or keep as historical exceptions?
12. **`m_max(::CanopySurface)` provisional vs. proper** — placeholder is `2·Nstreams`. Proper LAD-driven derivation lives in CanopyOptics. Defer or do properly now?
13. **Source-type rename strings** — proposed `DefaultSolarSource` and `LegacySFISource`. Better names? Worth bikeshedding once Sanghavi sees the YAML examples.

---

<a id="12-references"></a>
## 12. References

- **VLIDORT 2.8.3 source:**
  - `vsup/vbrdf/vbrdf_sup_routines.f90` — `VBRDF_FOURIER` (numerical Fourier integration)
  - `vsup/vbrdf/vbrdf_sup_masters.f90:2305` — `NMOMENTS = 2·NSTREAMS - 1` invariant
  - `vlidort_focode/VFO_Master.f90` — first-order (exact SS) subsystem
  - `vlidort_def/vlidort_inputs_def.f90:444-481` — `DO_FOCORR*` flag hierarchy
- **Nakajima & Tanaka 1988**, *Algorithms for radiative intensity calculations in moderately thick atmospheres using a truncation approximation*, JQSRT 40(1), 51-69. Original δ-M / TMS scheme.
- **Spurr 2002**, *VLIDORT: A linearized pseudo-spherical vector discrete ordinate radiative transfer code* — VLIDORT linearization, relevant for the `f₁` ForwardDiff redesign.
- **Kuusk 1991**, *The hot spot effect in plant canopy reflectance*. Original Kuusk hotspot model.
- **CanopyOptics.jl** `src/canopy_structure/hotspot.jl` — `AbstractHotSpot`, `KuuskHotSpot` with ForwardDiff hooks.
- **RAMI** (Radiation transfer Model Intercomparison) test scenes for canopy hotspot validation. HOM26 / HOM27 specifically.
- **Codex review of v0.1 and v0.2**, integrated throughout v0.3:
  - v0.1 issues 1–10 → resolved in v0.2/v0.3
  - v0.2 Issue 1 (`Nstreams` separation) → §2.2
  - v0.2 Issue 2 (canopy hotspot seam) → §4.3
  - v0.2 Issue 3 (`RTContext` lifecycle) → §1.3
  - v0.2 Issue 4 (`AbstractRTMode` orthogonal to outputs) → §6.2
  - v0.2 Issue 5 (type stability) → §1.2
  - v0.2 Issue 6 (real type names) → §2.6
