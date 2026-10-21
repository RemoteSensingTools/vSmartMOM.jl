# vSmartMOM.jl — Fourier decomposition, corrections, and dispatch architecture

**Document version:** v0.2
**Status:** working draft for review with Sanghavi and GPT.
**Scope:** part of the v2.0.0-rc1 release plan. Companion to `LINEARIZATION_BUGS.md` and the linearization redesign.
**Branch target:** `sanghavi-unified` (with `unified-vsmartmom` as the merge base for the architecture, physics layered in from `sanghavi`).
**Reference codes consulted:** VLIDORT 2.8.3 (Spurr et al., RT Solutions Inc.) for established BRDF Fourier-decomposition and SS-correction architecture; CanopyOptics.jl for the existing hotspot abstraction; the current `sanghavi-unified` branch of vSmartMOM.jl for present-state code.

**Iteration log:**
- **v0.2** — promoted the multiple-dispatch rule to governing principle (was an emergent observation in v0.1). Reorganized into two work steps. Folded in Codex review feedback: count vs. order convention fixed (`n_fourier_moments` is a count, matching current loops), VLIDORT translation off-by-one corrected, consumers-vs-construction wiring distinguished, Lambertian guard preserved, SS correction must read actual loop bound, lin/multisensor are not one-liners and become future output-strategy work, rename scope estimate corrected.
- **v0.1** — initial draft after VLIDORT/CanopyOptics/vSmartMOM code review.

---

## TL;DR

There is **one architectural rule** and **two steps** built on it.

**The rule (Section 1):** use multiple dispatch at physical and workflow boundaries; keep `if/else` only for scalar numerical cases. Surface physics, Fourier bounds, corrections, hotspot strategies, output products are *strategy objects with methods*. `if brdf isa CoxMunkSurface` in `rt_run.jl` is the smell this rule eliminates; `if m == 0` for Fourier weight normalization is fine and stays.

**Step 1 — Fourier decomposition and corrections (Sections 2–4).** Apply the rule to three currently-broken or missing pieces:
- **Per-component `n_fourier_moments` trait** (Section 2). Replace the surface-blind aerosol-only `max_m_bands` derivation with dispatched traits asking each component for its Fourier-moment count. Following VLIDORT, BRDF kernels saturate at `2·Nquad` moments (= `l_trunc + 1` for the standard `Nquad = (l_trunc+1)÷2`).
- **`AbstractSSCorrection` strategies** (Section 3). Promote the `apply_ss_correction!` hack inside `coxmunk_surface.jl` to first-class strategy structs in the new `FourierDecomposition` module. Default `NoSSCorrection`; concrete `BRDFSpecularCorrection`, `AerosolForwardPeakCorrection`, `ExternalSSCorrection`. One unified call site replaces the current scattered conditionals.
- **Wire `AbstractHotSpot` into `CanopySurface`** (Section 4). CanopyOptics already ships `AbstractHotSpot`/`NoHotSpot`/`KuuskHotSpot` with ForwardDiff hooks; vSmartMOM's canopy path uses none of it. Wire it on the SFI direct-beam term where the joint gap probability is physically active.

**Step 2 — Output products as dispatched strategies (Section 5).** Apply the same rule to a much bigger leak: `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` are 1004 lines of forked variations of the same loop, with all output arrays unconditionally pre-allocated regardless of what the user requested. Recast outputs as `AbstractRTOutput` strategy structs each owning `(allocate, accumulate!, finalize)` methods. Standard forward RT requests `(DirectionalTOA(),)`; RAMI requests `(HDRFOutput(), BHROutput())`; multisensor becomes `MultiSensorOutput(...)` passed into the same `rt_run` rather than its own file. The lin variant becomes another output strategy rather than a forked file. This is the deduplication that *unlocks* the linearization redesign.

Step 1 lands first because it's smaller, has a clean correctness story (VLIDORT validation), and proves out the dispatch idiom on contained pieces. Step 2 follows because it's larger, riskier (touches every output path), and is best done after the dispatch vocabulary is established by Step 1.

---

## 1. Dispatch architecture (the governing rule)

### 1.1 The rule

> **Use multiple dispatch at physical and workflow boundaries; keep `if/else` only for scalar numerical cases.**

Physical and workflow boundaries are exactly the places where the codebase currently uses `isa` checks or symbol-comparison conditionals. Concretely, they are:

- Surface physics (`if brdf isa CoxMunkSurface`, `if brdf isa CanopySurface`)
- Fourier bounds per component (currently surface-blind, ignoring the surface entirely)
- Corrections — SS, hotspot, future thermal/Raman (`if brdf isa CoxMunkSurface && SFI` in `rt_run.jl:311`)
- Output products (`if SFI` decisions threaded everywhere; `if RAMI` branches in postprocessing)
- Output topology (separate `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl` files for what is structurally the same loop)

These all become method dispatch on small concrete strategy structs.

`if/else` is fine — and clearer than dispatch — for scalar numerical cases:

- Fourier normalization weight (`m == 0 ? 0.5 : 1.0`)
- Numerical guards (NaN, denorm, divide-by-zero protection)
- Cache-initialized vs. not
- Optional benchmark/diagnostic instrumentation
- **Runtime user-data predicates** (`if any(a -> a.fᵗ > 0.01, aerosols)` for default-correction selection — this is a data check on values, not a type gate)

The rule eliminates *type and physics gates*, not all branches. The litmus test: if the branch's truth is determined at model-construction time by which types are in play, dispatch on those types. If it's determined at evaluation time by numerical values or user toggles, branch.

### 1.2 The strategy struct shape

Strategy objects in this codebase are **small concrete structs held in tuples, dispatched via single-method functions.** Three implementation rules:

**Rule (a) — concrete tuples, not abstract vectors.** A `Vector{AbstractSSCorrection}` forces runtime dispatch in the inner loop. A `Tuple{BRDFSpecularCorrection{CoxMunkSurface{Float64},Float64}}` lets the compiler specialize. Construction must produce concrete tuple types:

```julia
# Wrong — abstract eltype, runtime dispatch:
corrections = AbstractSSCorrection[BRDFSpecularCorrection(surf)]

# Right — concrete tuple, compile-time dispatch:
corrections = (BRDFSpecularCorrection(surf),)
```

When the set is built conditionally, build a tuple with `(items...,)` splat from a typed accumulator, or use `Tuple` construction with explicit element types. (Codex's Issue 8 in the v0.1 review correctly flagged the original sketch as type-unstable; this rule resolves it.)

**Rule (b) — single-method functions, not switch tables.** Each strategy owns its own behavior via dispatch:

```julia
apply_correction!(c::BRDFSpecularCorrection, ctx)   # one method per concrete strategy
apply_correction!(c::AerosolForwardPeakCorrection, ctx)
apply_correction!(c::NoSSCorrection, ctx) = nothing  # explicit no-op
```

The driver is one line:

```julia
foreach(c -> apply_correction!(c, ctx), corrections)
```

No `if c isa ...`, no `Symbol`-based switch. New correction types are added by writing one method.

**Rule (c) — naming consistency across strategy families.** Use parallel verb-form for each family:

| Family | Trait | Apply |
|---|---|---|
| Fourier bounds | `n_fourier_moments(component, ctx)` | (no apply, just a count) |
| SS corrections | `needs_correction(c, ctx)` | `apply_correction!(c, ctx)` |
| Output products | (no trait, presence in tuple is signal) | `accumulate_output!(state, output, ctx)` |
| Hotspot models | (no trait, presence in struct is signal) | `joint_gap_probability(model, …)` (already exists in CanopyOptics) |

Same grammar across families means contributors can guess the right name without consulting docs.

### 1.3 The `RTContext` object

Strategy methods need a consistent argument convention or the unification is theoretical. Introduce a `RTContext` struct constructed once before the Fourier loop, mutated as `m` advances, passed to every dispatched method:

```julia
mutable struct RTContext{FT, …}
    model::RTModel              # solver config, atmosphere, optics, surfaces
    iBand::Int                  # current band index
    m::Int                      # current Fourier order
    quad_points::QuadPoints
    composite_layer             # accumulated RT state
    added_layer_surface
    τ_sum_all::Matrix{FT}       # cumulative optical depth, m-independent
    geometry::Geometry          # μ₀, vza, vaz
    pol_type::AbstractPolarizationType
    arch::AbstractArchitecture
    workspace                   # InteractionWorkspace, allocated once
end
```

**Lifecycle contract (worth pinning down with Sanghavi):**

- Constructed once before the Fourier loop. Carries everything that doesn't depend on `m`.
- Field `m` is mutated as the loop advances. Strategy methods read it but don't mutate it.
- `composite_layer` accumulates across `m` iterations — that is the existing pattern, not a change.
- `τ_sum_all` is m-independent (the existing comment at `rt_run.jl:204` already establishes this).
- Constructed once, passed everywhere. No "thread `μ₀, vza, vaz, max_m, nSpec` through 7 positional arguments" — that is the current `apply_ss_correction!` smell.

Without this, every strategy method invents its own threading convention and the unification breaks down. With it, strategy struct definitions stay small (just the parameters specific to the strategy) because everything else is in `ctx`.

### 1.4 What the rule unlocks

This is the forward-looking part — naming the rule explicitly is what makes Step 2 a natural follow-on instead of a rewrite.

- **Step 1 sketches the idiom on small contained pieces** (~3 work items, each <100 lines moved).
- **Step 2 applies the same rule to output products**, collapsing 1004 lines of forked `rt_run*.jl` into a single coordinator + per-output strategy methods.
- **The linearization redesign uses the same rule**: lin Jacobians become an output strategy `LinearizedOutput(parameter_layout)` rather than a forked file, slotting into the same `rt_run` as forward outputs.
- **Future corrections (thermal SS, Raman SS, multiple-scatter δ-M)** all fit the same `Abstract*Correction` shape without further architectural work.

The rule has compounding returns. The Fourier/SS-correction work is the first application; everything downstream gets cheaper.

### 1.5 Branches that stay (the "fine" list)

To preempt over-application of the rule:

- `m == 0` Fourier weight normalization (`weight = m == 0 ? 0.5/π : 1.0/π`). Numerical, stays.
- Loop initialization (`m == 0 ? FT(1.0) : FT(2.0)` for the BRDF integration prefactor). Numerical, stays.
- The Lambertian-surface kernel zeroing higher moments **must stay** (`if m == 0` guard inside `lambertian_surface.jl:37`). The Fourier loop driver still runs `m > 0` for the Rayleigh and aerosol contributions; the Lambertian *surface* must zero its own contribution for those iterations because it has no higher moments. This is a numerical truth about the surface kernel, not a physics gate. (Codex's Issue 4 in the v0.1 review correctly flagged that deleting this guard would be wrong; this rule's framing makes clear why it's not the kind of branch the rule targets.)
- Cache initialization (`canopy._cache === nothing`). State, not type/physics dispatch.
- `if SFI` for whether to compute source-function integration *at all*. **This is a workflow boundary and is the prime target for Step 2** — `SFI` becomes presence/absence of `DirectionalTOA()` in the output tuple.

The rule eliminates type and physics gates. It does not eliminate all branches.

---

## 2. Per-component Fourier-moment count

### 2.1 The convention question (must settle first)

vSmartMOM's existing code uses `max_m` as a **count**: loops are `for m = 0:max_m - 1` in `rt_run.jl:208`, `rt_run_lin.jl:197`, the Cox-Munk SS reconstruction at `coxmunk_surface.jl:522`, and elsewhere. The trait being introduced **must match this convention** to avoid breaking the existing loops.

**Decision:** the trait is named `n_fourier_moments` and returns a count. Loop convention stays `for m = 0:n - 1`.

| Component | `n_fourier_moments` (count) |
|---|---|
| Lambertian (any flavor) | 1 (only `m = 0` contributes) |
| Rayleigh | 3 (`m = 0, 1, 2`) |
| AerosolOptics (with greek_coefs.β of length L+1) | min(L, l_trunc) + 1 |
| Cox-Munk, RPV, Ross-Li | `2·Nquad` (= `l_trunc + 1` for the standard `Nquad = (l_trunc+1)÷2`) |
| Canopy | `l_trunc + 1` (provisional — proper LAD-driven derivation deferred) |

This corrects the v0.1 draft, which used `m_max` as a highest-order quantity (Lambertian = 0). Under the count convention, Lambertian = 0 would mean "loop zero times" — silently disabling Lambertian surfaces. (Codex's Issue 1.)

### 2.2 Why `2·Nquad` for the BRDFs (the user's physics question)

VLIDORT 2.8.3, `vbrdf_sup_masters.f90:2305`:

```fortran
NMOMENTS = 2 * NSTREAMS - 1
```

Followed by a Fortran-style **inclusive** loop `DO M = 0, NMOMENTS`. That runs `0..2·NSTREAMS - 1`, which is **`2·NSTREAMS` iterations**. Translated to vSmartMOM's count convention, with `NSTREAMS ↔ Nquad`:

```
n_fourier_moments_BRDF = 2 · Nquad
```

For the standard `Nquad = (l_trunc + 1) ÷ 2`, this evaluates to `l_trunc + 1` (for even `l_trunc`; integer division shaves one for odd). 

**Rationale:** the BRDF Fourier series is projected onto the discrete-ordinate basis during the RT solve. A polar grid with `Nquad` points cannot resolve azimuthal modes higher than `2·Nquad - 1` regardless of how rich the kernel itself is. Anything past that aliases. The cap is set by *the solver's resolving power*, not by the surface's intrinsic angular complexity. This is why VLIDORT uses the same formula for all 17 of its BRDF kernels.

(Codex's Issue 2 in the v0.1 review caught the off-by-one in v0.1's `m_max = l_trunc` formula. Under the count convention, the correct value is `l_trunc + 1`, which works out to `2·Nquad`. I'm derivation-first here so we don't carry the same bug into v0.3.)

### 2.3 Cox-Munk, RPV, Ross-Li — same formula, different reasons

For Cox-Munk specifically, you cannot get a closed-form Fourier expansion: the slope distribution couples non-trivially with the Fresnel matrix and an azimuth-dependent shadowing term. There's no clean recursion the way Rayleigh has.

RPV and Ross-Li are slightly better behaved analytically — the Rahman-Pinty-Verstraete kernel and the Ross-Li volumetric/geometric kernels have known analytic forms in `cos ϕ`, `cos 2ϕ` pieces — but every production code (VLIDORT, 6SV1, MYSTIC) integrates them numerically anyway because:

- The cost is negligible compared to the RT solve.
- It unifies the code path across surfaces.
- It lets you swap kernels without rederiving Fourier coefficients.

So the same `2·Nquad` count applies to all three, computed via the same numerical-azimuth-quadrature recipe.

### 2.4 The two grids

Worth keeping these separate in the module documentation so future contributors don't merge them:

| Quantity | Role | Where set | Notes |
|---|---|---|---|
| `Nquad` | Polar discrete-ordinate count | `rt_set_streams.jl:32` | Already correct; user's invariant |
| `nQuad_BRDF` | Azimuth integration grid for each Fourier coefficient | Hard-coded `100` in `coxmunk_surface.jl`, `rpv_surface.jl` | Quadrature accuracy, not physical limit. VLIDORT default is also 100. |
| `n_fourier_moments` | Count of Fourier orders the loop runs over | `model_from_parameters.jl:263` (currently aerosol-blind) | The trait being added |

`nQuad_BRDF` is a property of the BRDF's roughness, not the solver. CoxMunk at wind speed 0.5 m/s might genuinely need 200 points; CoxMunk at 15 m/s is fine with 50. **Recommended:** promote it to a per-surface field with default 100 (Step 1.5, deferred). Doesn't block anything.

### 2.5 What's already half-built

| Piece | Status | Location |
|---|---|---|
| Polar stream rule `Nquad = (l_trunc+1)÷2` | Correct | `rt_set_streams.jl:32` |
| Per-mode azimuth integration for Cox-Munk | Correct, `nQuad = 100` hard-coded | `coxmunk_surface.jl:_fourier_coeff_element` |
| Per-mode azimuth integration for RPV | Correct, `nQuad = 100` hard-coded | `rpv_surface.jl:reflectance` |
| Canopy `compute_Z_matrices(canopy, μ, 0:m_max)` | Cleanly parameterized | CanopyOptics `z_matrices.jl` |
| Per-band `max_m_bands` from aerosol Greek coefs | Computed but **surface-blind** | `model_from_parameters.jl:254-263` |
| `n_fourier_moments(::Surface)` dispatch trait | Doesn't exist | — |
| **Consumers reading `max_m_bands`** | **Broken — return scalar instead** | `types.jl:906` `get_max_m` returns `solver.max_m`, not the per-band vector |

The last row is critical and was missed in v0.1. (Codex's Issue 3.) Even with `max_m_bands` populated correctly in `model_from_parameters`, the Fourier loop driver in `rt_run.jl:97` calls `get_max_m(model)` which returns the *scalar*, not the per-band vector. The `model.max_m` property forwarder at `types.jl:922` does give the per-band vector, but it's not what the loop reads.

**Wiring change needed:** `rt_run.jl`, `rt_run_lin.jl`, and `rt_run_multisensor.jl` consumers must read the per-band value, not the scalar. Easiest way: replace `max_m = get_max_m(model)` with `n_moments = n_fourier_moments(model, iBand)` (a new method on the trait). This dispatches on the model + band, returns the per-band count, and goes through the same trait used for default-construction.

### 2.6 Proposed trait methods

```julia
# CoreRT/FourierDecomposition/n_fourier_moments.jl

"""
    n_fourier_moments(component, ctx) -> Int

Number of Fourier moments (count) needed for an accurate RT solve given the
context. Loop convention is `for m = 0 : n_fourier_moments(...) - 1`.

The framework caps the per-band aggregate at `2·Nquad`, the highest count
the polar discrete-ordinate basis can resolve (VLIDORT's NMOMENTS+1
identity). Anything above this aliases.
"""
n_fourier_moments(::Rayleigh,                ctx) = 3

n_fourier_moments(::LambertianSurfaceScalar,   ctx) = 1
n_fourier_moments(::LambertianSurfaceSpectrum, ctx) = 1
n_fourier_moments(::LambertianSurfaceLegendre, ctx) = 1
n_fourier_moments(::LambertianSurfaceSpline,   ctx) = 1

# VLIDORT vbrdf_sup_masters.f90:2305 invariant: NMOMENTS = 2·NSTREAMS - 1,
# inclusive loop → 2·NSTREAMS iterations. Translated to vSmartMOM's count
# convention with NSTREAMS ↔ Nquad = (l_trunc+1)÷2.
n_fourier_moments(::CoxMunkSurface,           ctx) = 2 * Nquad(ctx)
n_fourier_moments(::rpvSurfaceScalar,         ctx) = 2 * Nquad(ctx)
n_fourier_moments(::RossLiSurfaceScalar,      ctx) = 2 * Nquad(ctx)

# Provisional placeholder. Proper LAD-driven derivation lives in CanopyOptics
# and depends on the leaf angle distribution; for now match the BRDF default.
n_fourier_moments(::CanopySurface,            ctx) = 2 * Nquad(ctx)

# From Greek-coef tail length, capped at 2·Nquad.
n_fourier_moments(a::AerosolOptics,           ctx) =
    min(length(a.greek_coefs.β), 2 * Nquad(ctx))
```

Aggregator (per-band):

```julia
"""
    n_fourier_moments(model, iBand) -> Int

Per-band Fourier-loop count, derived from all components present in the band.
"""
function n_fourier_moments(model::RTModel, iBand::Int)
    ctx = (; l_trunc = model.solver.l_trunc, Nquad = model.quad_points.Nquad)
    n = n_fourier_moments(model.optics.rayleigh, ctx)
    for a in model.optics.aerosols.aerosol_optics[iBand]
        n = max(n, n_fourier_moments(a, ctx))
    end
    n = max(n, n_fourier_moments(model.surfaces[iBand], ctx))
    return min(n, 2 * ctx.Nquad)   # cap
end
```

User override stays as `model.solver.max_m`: if explicitly set, take `min(n, max_m)`.

### 2.7 Knock-on simplifications

After the trait lands and the consumers are wired:

- **`apply_correction!` argument list shortens.** It reads `n_fourier_moments(model, iBand)` from `ctx`, no positional `max_m` argument. **And critically, it reads the same value the loop used** — eliminating the drift Codex's Issue 5 warned about (where my v0.1 plan had the correction recompute `m_max` and risked subtracting a different sum than was added).
- **Lambertian-only bands skip `m > 0` work.** The loop driver can short-circuit:
  ```julia
  for m = 0 : n_fourier_moments(model, iBand) - 1
      ...
  end
  ```
  With `n_fourier_moments = 1` for Lambertian-only bands, the loop runs once, and the Lambertian surface kernel's existing `m == 0` guard inside `lambertian_surface.jl` stays intact (it's the right guard for the right reason).
- **`vSmartMOM_Parameters.max_m`** can be deprecated to a "user override only" knob, with the trait as the default.

---

## 3. Single-scattering corrections as dispatched strategies

### 3.1 The current state

`rt_run.jl:310-315`:

```julia
if brdf isa CoxMunkSurface && SFI
    @timeit "SS Correction" apply_ss_correction!(
        R_SFI, brdf, pol_type, vza, vaz, μ₀,
        Array(τ_sum_all[:,end]), max_m, nSpec)
end
```

Three problems, all of which the dispatch rule eliminates:

1. **`isa` gate.** Adding any new BRDF that needs a similar correction (Cox-Munk-Bréon, glitter, etc.) means another `if` branch.
2. **Inconsistent across entry points.** Present in `rt_run.jl`, missing from `rt_run_lin.jl` and `rt_run_multisensor.jl`. Missing from `rt_run_lin.jl` is a silent Jacobian bug at glint geometries — the forward radiance gets the correction but its tangent doesn't.
3. **`if brdf isa CoxMunkSurface` hides physics.** The same kind of correction is *also* needed for δ-M-truncated aerosol forward peaks. vSmartMOM has the first half (`delta_m_truncation.jl`) but not the second (the exact-SS replacement). The current SS correction was solving Cox-Munk specifically without anyone noticing the same pattern was missing for aerosols.

### 3.2 What VLIDORT does (the reference)

VLIDORT factors single-scattering into a first-class subsystem: **FO — First-Order Code, exact single-scattering and direct-thermal**. Peer of the multiple-scatter solver. Control surface in `vlidort_inputs_def.f90:444-481`:

- `DO_FOCORR` — master on/off
- `DO_FOCORR_EXTERNAL` — external pre-computed exact-SS hook
- `DO_FOCORR_NADIR` / `DO_FOCORR_OUTGOING` — sphericity treatment
- `DO_SSCORR_TRUNCATION` — whether SS correction also accounts for δ-M truncation
- `DO_SSCORR_USEFMAT` — exact phase function from F-matrix vs. Greek coefficients

Lessons:
1. SS correction is a **subsystem with its own master flag**, not a per-surface afterthought.
2. **Not surface-specific**. Same pass corrects truncated aerosol forward-peaks and truncated surface specularities.
3. There is an **"external" hook** for users to override the built-in correction.

### 3.3 Proposed strategy hierarchy

```julia
# CoreRT/FourierDecomposition/ss_correction.jl

"""
    AbstractSSCorrection

Marker for a single-scattering correction strategy applied after the Fourier
loop. The Fourier expansion truncates contributions with sharp angular features
(BRDF specular peaks, aerosol forward peaks); an SS correction restores them by
replacing the truncated single-scatter contribution with its exact evaluation.
"""
abstract type AbstractSSCorrection end

"""Default no-op."""
struct NoSSCorrection <: AbstractSSCorrection end

"""
    BRDFSpecularCorrection{S<:AbstractSurfaceType}

Replaces the truncated Fourier reconstruction of the surface BRDF SS
contribution with its exact evaluation. Nakajima-Tanaka TMS applied to
surfaces — generalizes the existing Cox-Munk hack.
"""
struct BRDFSpecularCorrection{S<:AbstractSurfaceType} <: AbstractSSCorrection
    surface::S
    nQuad_BRDF::Int   # default 100
end

"""
    AerosolForwardPeakCorrection{FT}

The δ-M companion: replaces the truncated SS aerosol contribution with the
exact phase function evaluation. **Not in defaults until Step 1.7 lands.**
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

Single-method dispatch:

```julia
apply_correction!(::NoSSCorrection, ctx) = nothing

function apply_correction!(c::BRDFSpecularCorrection{<:CoxMunkSurface}, ctx)
    # Existing coxmunk_surface.jl:505-546 body, parameterized over ctx.
    # Reads ctx.n_moments_band so it subtracts the SAME Fourier sum the loop added.
    # ...
end

# Future:
function apply_correction!(c::AerosolForwardPeakCorrection, ctx)
    # ...
end

function apply_correction!(c::ExternalSSCorrection, ctx)
    c.callback(ctx)
end
```

Per-surface trait controlling defaults:

```julia
has_specular_peak(::AbstractSurfaceType) = false
has_specular_peak(::CoxMunkSurface)      = true
```

Default selection (note the explicit deferral of aerosol correction):

```julia
function default_ss_corrections(surfaces, aerosol_optics)
    map(eachindex(aerosol_optics)) do b
        cs = ()
        if has_specular_peak(surfaces[b])
            cs = (cs..., BRDFSpecularCorrection(surfaces[b], 100))
        end
        # AerosolForwardPeakCorrection NOT included until Step 1.7 implements it.
        # See Codex Issue 6 — type defined ≠ ready to enable.
        isempty(cs) ? (NoSSCorrection(),) : cs
    end
end
```

The construction returns a `Vector` of concretely-typed tuples. Each tuple's element type is determined by which strategies are actually selected for that band, so dispatch in the apply loop is compile-time per band. (Codex Issue 8 resolution.)

### 3.4 The unified call site

`rt_run.jl:310-315` becomes:

```julia
SFI && @timeit "SS corrections" foreach(c -> apply_correction!(c, ctx),
                                        model.solver.ss_corrections[iBand])
```

For the lin run, same call, but `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::LinContext)` is a separate method that updates both `R_SFI` and `Ṙ_SFI` consistently. **Lin SS correction is a real implementation work item, not a wiring one.** (Codex Issue 7 resolution. The v0.1 doc misleadingly called the lin wiring "the same one-liner" — that's only the *call site*; the dispatched method has its own derivative math.)

For multisensor, similar: `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::MultiSensorContext)` reaches the per-sensor outputs, with the attenuation path computed at each sensor level rather than only at TOA. Also a real implementation work item.

These extra methods are why the lin and multisensor SS-correction wiring is not free, but they are *additive* methods on the existing strategy types — no architectural rework, just more methods.

### 3.5 What this gives the linearization redesign

Each correction is its own self-contained unit of differentiation. Wrap `apply_correction!` in a ForwardDiff seed → get `∂(correction)/∂(Mie params, surface params)` for free → add to the chain rule via `lin_added_layer_all_params.jl`. The current implementation makes that hard because the SS correction hack happens *outside* the lin seam. With the strategy structs, the hack becomes a typed unit the linearization framework can see.

### 3.6 Caching deferred

Current Cox-Munk `apply_ss_correction!` re-evaluates Fourier coefficients in `_fourier_coeff_element` rather than caching from the main loop. Wasteful but safe.

**Recommendation:** keep recomputing for the first landing. Abstraction first, optimization second. Cost is negligible compared to the RT solve. Add caching later once the structure is in place.

---

## 4. Canopy hotspot — same idiom, different seam

### 4.1 Why this is *not* a subtype of `AbstractSSCorrection`

The temptation is to add `KuuskHotSpotCorrection <: AbstractSSCorrection` and call it done. **It would be wrong.** The same temptation will recur for other corrections (Raman SS, thermal SS, multiple-scatter δ-M); the answer is the same.

The two corrections share a *design pattern* (the dispatch rule from Section 1) but live on **completely different seams** in the RT pipeline:

| | `AbstractSSCorrection` | `AbstractHotSpot` |
|---|---|---|
| **When** | After Fourier loop completes | Inside `elemental_canopy!` SFI source assembly |
| **Where** | Post-processing on `R_SFI` | Per-layer · per-Fourier-mode |
| **How** | Additive: `(M_exact − M_fourier)` per geometry | Multiplicative: `C_hs · P_so/P_o` on joint gap probability |
| **What it patches** | Fourier truncation losing specular-peak resolution | Beer-Lambert assuming independent gaps |
| **Where the math lives** | `vSmartMOM/CoreRT/FourierDecomposition` | `CanopyOptics/canopy_structure/hotspot.jl` |

If `AbstractHotSpot` were forced into `AbstractSSCorrection`, the `apply_correction!` hook would either widen to handle in-loop calls (abstraction leaks) or operate on `R_SFI` after the fact and reverse-engineer what should have happened inside source assembly (hack).

The conceptual similarity is captured by **shared dispatch idiom**, not shared type.

### 4.2 What CanopyOptics already ships

`CanopyOptics.jl/src/canopy_structure/hotspot.jl` defines:

- `abstract type AbstractHotSpot{FT<:Real}`
- `struct NoHotSpot{FT}` — independent gaps, `C_hs = 1`
- `struct KuuskHotSpot{FT}` with size parameter `h`
- `joint_gap_probability(model, k_s, k_o, μ_s, μ_o, dϕ, L)` — `exp(-(k_s+k_o)L) · C_hs(L)`
- ForwardDiff hooks: `_hotspot_value(x::ForwardDiff.Dual)` already in place

The author was thinking ahead about linearization. Kuusk math is a closed-form 4-line expression in `(h, k_s, k_o, μ_s, μ_o, dϕ, L)`; Jacobians w.r.t. `h` (and through `k_s, k_o` to LAI/G) drop out cleanly under ForwardDiff.

**The abstraction is shipped. The implementation is correct. The linearization story is anticipated. What's missing is wiring on the vSmartMOM side.**

### 4.3 What's missing in vSmartMOM (and what's harder than v0.1 implied)

Greps confirm zero references to `hotspot`, `joint_gap`, `kuusk` anywhere in `canopy_surface.jl`, `elemental_canopy.jl`, `interaction_hdrf.jl`. `elemental_canopy!` uses independent Beer-Lambert via `arr_type(τ_sum)` multiplication.

**Wiring is harder than the v0.1 doc implied.** (Codex Issue 9.) `elemental_canopy!` works in **Fourier/quadrature space**, while `joint_gap_probability` is geometry-dependent (view direction, relative azimuth). The plain `P_so/P_o` sketch in v0.1 is not enough.

Two paths to consider, one acceptable, one not:

**(a) Per-Fourier-mode hotspot expansion (heavy).** The hotspot has its own Fourier expansion in `dϕ` — Kuusk's `α(μ_s, μ_o, dϕ, h)` is a function of relative azimuth, so it expands in `cos(m·dϕ)`. Doable but it's its own decomposition pass with its own quadrature, derivative chain, and validation surface. Substantial new physics work.

**(b) Hotspot only on the SFI direct-beam term, not on the diffuse multiple-scatter term (recommended).** This is what VLIDORT-style codes typically do. Physical justification: the hotspot is a *direct-beam-only* effect — it's the joint gap probability for *correlated* sun/view paths, and only the SFI term sees both paths simultaneously. The diffuse multiple-scatter sees decorrelated paths, so independent Beer-Lambert is correct there.

Under (b), the hotspot multiplies the direct-beam source `j₀⁺` (the SFI branch in `elemental_canopy!`'s `get_canopy_elem_rt_SFI!` kernel), and the multiple-scatter `r⁻⁺, t⁺⁺, ...` matrices stay untouched. Geometry `(μ_v, dϕ)` is available at the SFI assembly point because the SFI kernel already produces directional radiance for the requested `vza × vaz`.

**Recommendation: path (b).** Reasons:
- Matches the dominant convention in established RT codes (VLIDORT, 6SV1).
- Math is local — one multiplicative correction in `get_canopy_elem_rt_SFI!`, no new Fourier pass.
- Linearization through `joint_gap_probability` is already supported in CanopyOptics via the ForwardDiff hooks.
- Validates against published RAMI test cases (HOM26/HOM27) which use this convention.

**Confirm with Sanghavi:** path (b) is the right choice. If RAMI validation requires the per-mode expansion, fall back to path (a) as a Step-1.6.5 deferral.

### 4.4 Proposed wiring (path b)

```julia
# CoreRT/types.jl — add hotspot field to CanopySurface
mutable struct CanopySurface{FT, ...} <: AbstractSurfaceType
    # ... existing fields ...
    hotspot::AbstractHotSpot{FT}   # default NoHotSpot{FT}() — bit-for-bit current behavior
end
```

In `elemental_canopy.jl`'s `get_canopy_elem_rt_SFI!` kernel, where the SFI source `j₀⁺, j₀⁻` is currently assembled with implicit Beer-Lambert via `τ_sum`:

```julia
# Today (independent gaps via τ_sum multiplication):
#   j₀⁺[...] = ... * exp_τ_sum_term * ...

# Proposed — multiplicative correction on the SFI branch only:
P_so = joint_gap_probability(canopy.hotspot, k_s, k_o, μ_s, μ_v, dϕ, τ_sum)
P_o  = exp(-k_o * τ_sum)
direct_source_scale = P_so / P_o
# multiply into the direct-beam source contribution
```

With `NoHotSpot()` (the default), `direct_source_scale = 1.0` exactly, and behavior is bit-for-bit identical. With `KuuskHotSpot()`, the user gets the corrected joint gap.

Same wiring needed in `interaction_hdrf_canopy!` for HDRF/BHR.

### 4.5 Validation and regression risk

**This is a correctness fix, not a refactor.** Current vSmartMOM canopy results compared against literature with hotspot effects (RAMI test cases, NDVI/red-edge over crop canopies, MODIS BRDF retrievals at small phase angle) are wrong by whatever `C_hs` would have been.

⚠️ **Grep the vSmartMOM test suite before quietly fixing this.** If any reference benchmark was tuned to match the current independent-gap result, that benchmark needs to be regenerated with the corrected physics, and the discrepancy needs to be documented in the changelog.

Validation references:
- **RAMI** test scenes HOM26 / HOM27 (homogeneous canopy with explicit hotspot tests)
- **PROSAIL-RT** turbid-medium reference

---

## 5. Step 2 — Output products as dispatched strategies

This is the larger refactor and lands after Step 1 has established the dispatch idiom on small contained pieces.

### 5.1 The leak

`rt_run.jl` (lines 136-140 and elsewhere) unconditionally allocates `R`, `T`, `R_SFI`, `T_SFI`, `ieR_SFI`, `ieT_SFI`, `hdr`, `bhr_uw`, `bhr_dw`. The return tuple at line 325 is conditional on `SFI`, but the allocations aren't:

```julia
return SFI ? (R_SFI, T_SFI, ieR_SFI, ieT_SFI, hdr, bhr_uw[1,:], bhr_dw[1,:]) : (R, T)
```

Three problems:

1. **Every run pays the allocation cost for every output product**, whether it's used or not.
2. **Adding a new output product requires touching every `rt_run*.jl` file** (currently 1004 lines across three files for what is structurally the same loop).
3. **The choice is implicit** in `if SFI`, `if RAMI`, scattered throughout postprocessing — not a workflow boundary the user controls explicitly.

### 5.2 The strategy

Recast each output as an `AbstractRTOutput` strategy struct owning three methods:

```julia
abstract type AbstractRTOutput end

struct DirectionalTOA <: AbstractRTOutput end
struct DirectionalBOA <: AbstractRTOutput end
struct HDRFOutput     <: AbstractRTOutput end
struct BHROutput      <: AbstractRTOutput end
struct MultiSensorOutput{...} <: AbstractRTOutput
    sensor_levels::Vector{...}
end
struct LinearizedOutput{PL} <: AbstractRTOutput
    parameter_layout::PL
end

allocate_output(::DirectionalTOA, model, nSpec) = ...
accumulate_output!(state, ::DirectionalTOA, ctx) = ...
finalize_output(state, ::DirectionalTOA, ctx) = ...
```

`rt_run` becomes the coordinator:

```julia
function rt_run(model, outputs::Tuple; iBand = 1)
    states = map(o -> allocate_output(o, model, nSpec), outputs)
    ctx = RTContext(model, iBand, ...)

    for m = 0 : n_fourier_moments(model, iBand) - 1
        ctx.m = m
        # ... existing layer + surface solve ...
        foreach((s, o) -> accumulate_output!(s, o, ctx), states, outputs)
    end

    foreach(c -> apply_correction!(c, ctx), model.solver.ss_corrections[iBand])

    return map((s, o) -> finalize_output(s, o, ctx), states, outputs)
end
```

### 5.3 What this collapses

| Current entry point | What replaces it |
|---|---|
| `rt_run(model)` | `rt_run(model, (DirectionalTOA(), DirectionalBOA()))` |
| `rt_run` with `SFI=true` and RAMI | `rt_run(model, (DirectionalTOA(), HDRFOutput(), BHROutput()))` |
| `rt_run_multisensor.jl` (177 lines, separate file) | `rt_run(model, (MultiSensorOutput(levels),))` |
| `rt_run_lin.jl` (303 lines, separate file) | `rt_run(model, (DirectionalTOA(), LinearizedOutput(layout)))` |

The 1004 lines across three files collapse into one coordinator (~150 lines) plus per-output strategy methods. **The lin variant disappearing as a separate file is the substantial deduplication that unlocks the linearization redesign.**

### 5.4 Migration strategy

Step 2 cannot land in one commit. Order:

1. Define `AbstractRTOutput` and `RTContext`. No call sites changed.
2. Implement `DirectionalTOA`, `DirectionalBOA`. Both dispatch to existing `postprocessing_vza!` internals. No behavior change.
3. Refactor current `rt_run` to take `outputs::Tuple` argument with a default of `(DirectionalTOA(), DirectionalBOA())` for backward compat. Existing callers unaffected.
4. Implement `HDRFOutput`, `BHROutput`. Migrate RAMI workflows. Validate against existing RAMI benchmarks.
5. Implement `MultiSensorOutput`. Delete `rt_run_multisensor.jl`. Validate against existing multisensor cases.
6. Implement `LinearizedOutput`. Delete `rt_run_lin.jl`. **This is the riskiest step — the lin path has the most methods to move.** Validate against `LINEARIZATION_BUGS.md` test set.

Each step is independently shippable; benchmarks at each step boundary.

### 5.5 Why this is Step 2, not Step 1

Three reasons:

1. **Risk surface.** Step 1's three work items each touch <100 lines. Step 2 touches 1004 lines across three files. Doing them simultaneously means mixing low-risk and high-risk changes in the same commit boundary.
2. **Validation surface.** Step 1's benchmarks are local (Cox-Munk SS, RAMI canopy, 6SV1/Natraj for Fourier moments). Step 2's benchmarks span every output path that exists. Doing them after Step 1 means the dispatch idiom is already validated in three small applications before being applied to the big one.
3. **Vocabulary.** Step 2's strategy structs (`DirectionalTOA`, `MultiSensorOutput`, etc.) follow the same naming conventions as Step 1's (`BRDFSpecularCorrection`, etc.). Doing Step 1 first lets contributors see the pattern in a small example before extending it to a big one.

---

## 6. Execution plan

### Step 1 — Fourier work (this round)

| # | What | Risk | Validation |
|---|---|---|---|
| 1.1 | Rename `AbstractFourierDecompositionType` → `AbstractMieDecompositionAlgorithm`. Add deprecation alias. | Touches exported API. ~50–100 lines including tests, docs, IO imports. (Codex Issue 10 — bigger than v0.1 claimed.) | All existing tests pass. |
| 1.2 | Introduce `FourierDecomposition` module skeleton. | Plumbing. | Module loads. |
| 1.3 | Add `RTContext` struct (Section 1.3). | New abstraction. | Trivial unit tests. |
| 1.4 | Add `n_fourier_moments` trait + methods (Section 2.6). | New code. | Unit tests for each surface type, count convention asserted. |
| 1.5 | **Wire `n_fourier_moments(model, iBand)` into the consumers** — `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`. (Codex Issue 3.) Replace `get_max_m(model)` with the trait call. Run 6SV1/Natraj reference benchmarks before and after — Cox-Munk bands may now compute more moments than before. | **Behavior change risk.** | 6SV1, Natraj benchmarks. |
| 1.6 | Move BRDF Fourier integration to shared helper (`_brdf_fourier_coefficient` in `FourierDecomposition`). | Pure refactor. | Existing tests pass bit-for-bit. |
| 1.7 | Add `AbstractSSCorrection` hierarchy (Section 3.3). Migrate Cox-Munk `apply_ss_correction!` to `BRDFSpecularCorrection{<:CoxMunkSurface}`. **Read `n_moments` from ctx**, not via recomputation. (Codex Issue 5.) | Pure relocation. | Existing Cox-Munk tests pass bit-for-bit. |
| 1.8 | Wire `apply_correction!` into `rt_run.jl` via `foreach`. Add `default_ss_corrections` (without `AerosolForwardPeakCorrection` — Codex Issue 6). | Pure call-site swap. | Existing benchmarks unchanged. |
| 1.9 | Add lin-context method `apply_correction!(::BRDFSpecularCorrection{<:CoxMunkSurface}, ::LinContext)` updating `R_SFI` and `Ṙ_SFI` consistently. (Codex Issue 7 — real implementation, not one-liner.) | New physics in lin path. | Validate Jacobians against finite-difference at glint geometry. |
| 1.10 | Add multisensor-context method, similar pattern. Per-sensor-level attenuation, not TOA. | New physics in multisensor path. | Validate against multisensor reference cases if any exist. |
| 1.11 | Wire `AbstractHotSpot` into `CanopySurface` via path (b) (Section 4.4). Multiplicative on SFI direct-beam. Default `NoHotSpot()` preserves current behavior bit-for-bit. | **Correctness fix in disguise.** Grep test suite first for hotspot-dependent benchmarks. | RAMI HOM26/HOM27. |
| 1.12 | Implement `AerosolForwardPeakCorrection`. Genuine new physics (the missing δ-M companion). | New physics. | PyVLIDORT high-AOD strongly-forward-peaked aerosol cases. |
| 1.13 | (Deferred) Promote `nQuad_BRDF` to per-surface field. | Trivial. | Convergence sweeps. |
| 1.14 | (Deferred) Cache Fourier coefficients for SS correction reuse. | Optimization only. | No behavior change. |

**Step 1 boundaries for benchmarks:**
- After 1.1–1.4: all existing tests must pass bit-for-bit. Pure plumbing.
- After 1.5: 6SV1 / Natraj benchmarks for any band where `n_moments` could change.
- After 1.8: any existing Cox-Munk TOA reference cases.
- After 1.9: lin Jacobian regression at glint geometries.
- After 1.10: multisensor + Cox-Munk reference cases.
- After 1.11: RAMI HOM26/HOM27.
- After 1.12: PyVLIDORT high-AOD cases.

### Step 2 — Output products as strategies (next round)

Per Section 5.4. Six-substep migration, each independently shippable, each with its own benchmark boundary. Order: `Directional*` → `HDRF/BHR` → `MultiSensor` → `Linearized`. The lin migration is the riskiest substep and lands last.

**Step 2 starts after Step 1.11 lands** (the canopy hotspot wiring is the last Step 1 piece that can stay forward-only). Step 1.12 (aerosol forward peak) and Step 2 can proceed in parallel — they touch different code paths.

---

## 7. Open questions for Sanghavi

In rough order of "blocks the next commit" → "nice to settle eventually":

1. **Count convention confirmation.** v0.2 uses `n_fourier_moments` as a count to match existing loops (`for m = 0:n-1`). This means Lambertian = 1, Rayleigh = 3, Cox-Munk = `2·Nquad`. Is this consistent with how you think about Fourier moments in your physics work, or should the trait expose the highest order with `n_moments = highest_order + 1` derived where needed? Either works; the codebase loops decide. **Blocks step 1.4.**

2. **VLIDORT translation sanity check.** `n_fourier_moments(::CoxMunk, ctx) = 2 · Nquad(ctx)`. For `l_trunc = 60, Nquad = 30`, this gives 60 moments (highest order 59). Does this match any internal reference numerics? Easiest verification: run a Cox-Munk + Rayleigh test with `n_moments` set to 60 vs. 30 vs. 100 and check which converges to the reference radiance.

3. **Hotspot wiring path.** Path (b) — multiplicative on SFI direct-beam only — is the recommended approach (Section 4.3). Path (a) — full per-Fourier-mode hotspot expansion — is heavier but more physically complete. RAMI HOM26/HOM27 should disambiguate. **Blocks step 1.11.**

4. **`ExternalSSCorrection` callback signature.** Current sketch: `c.callback(ctx)` returns nothing, mutates `R_SFI` in `ctx`. Is the mutation-via-context pattern OK with you, or do you want a return-value API? Affects how downstream tooling integrates.

5. **`AerosolForwardPeakCorrection` threshold.** Default `fᵗ > 0.01` to skip when truncation is negligible. Reasonable starting point? Worth tuning during validation.

6. **Naming harmonization.** Existing `AbstractTruncationType`/`NoTruncation` vs. proposed `Abstract*Correction`/`No*Correction`. Harmonize during this work, or keep them as historical exceptions? Affects how `δBGE` migrates from `Scattering` to `FourierDecomposition`.

7. **Canopy hotspot regression risk.** Any vSmartMOM reference benchmark whose published numbers depend on the missing-hotspot bug? Need to identify before quietly fixing. **Blocks step 1.11.**

8. **`m_max(::CanopySurface)` provisional vs. proper.** Placeholder is `2·Nquad`. Proper LAD-driven value lives in CanopyOptics. Worth doing properly now, or defer? I'd defer — placeholder doesn't change current numerics.

9. **`δBGE` and `NoTruncation` migration.** Prior draft flagged these as "probable second migration" from `Scattering/` to `FourierDecomposition/`. Confirm during this work or keep as separate followup?

10. **`RTContext` lifecycle.** Section 1.3 specifies fields and lifecycle. Worth confirming this matches your intuition about what's m-independent vs. m-mutated, and whether `composite_layer` belongs in `ctx` or stays as a separate threaded argument.

---

## 8. References

- **VLIDORT 2.8.3 source:**
  - `vsup/vbrdf/vbrdf_sup_routines.f90` — `VBRDF_FOURIER` routine (numerical Fourier integration)
  - `vsup/vbrdf/vbrdf_sup_masters.f90:2305` — `NMOMENTS = 2·NSTREAMS - 1` invariant
  - `vlidort_focode/VFO_Master.f90` — first-order (exact SS) subsystem
  - `vlidort_def/vlidort_inputs_def.f90:444-481` — `DO_FOCORR*` flag hierarchy
- **Nakajima & Tanaka 1988**, *Algorithms for radiative intensity calculations in moderately thick atmospheres using a truncation approximation*, JQSRT 40(1), 51-69. Original δ-M / TMS scheme.
- **CanopyOptics.jl** `src/canopy_structure/hotspot.jl` — `AbstractHotSpot`, `KuuskHotSpot`.
- **Kuusk 1991**, *The hot spot effect in plant canopy reflectance*. Original Kuusk hotspot model.
- **RAMI** (Radiation transfer Model Intercomparison) test scenes for canopy hotspot validation.
- **Spurr 2002**, *VLIDORT: A linearized pseudo-spherical vector discrete ordinate radiative transfer code* — VLIDORT linearization, relevant for the `f₁` ForwardDiff redesign.
- **Codex review of v0.1**, integrated throughout v0.2. Specific references:
  - Issue 1 — count vs. order convention → Section 2.1
  - Issue 2 — VLIDORT off-by-one → Section 2.2
  - Issue 3 — consumers reading scalar `max_m` → Section 2.5, Step 1.5
  - Issue 4 — Lambertian guard preservation → Section 1.5
  - Issue 5 — SS correction must read actual loop bound → Section 3.3, Step 1.7
  - Issue 6 — `AerosolForwardPeakCorrection` not in defaults until implemented → Section 3.3, Step 1.8
  - Issue 7 — lin/multisensor are not one-liners → Section 3.4, Steps 1.9–1.10
  - Issue 8 — type stability → Section 1.2
  - Issue 9 — canopy hotspot wiring is harder than v0.1 implied → Section 4.3
  - Issue 10 — rename scope → Step 1.1
