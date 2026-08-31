# Fourier decomposition, SS correction, and canopy hotspot — design notes

**Status:** working draft for review with Sanghavi and GPT.
**Scope:** part of the `FourierDecomposition` module proposed in the v2.0.0-rc1 release plan (Step 0.5 in the prior refactor scoping doc). Companion to `LINEARIZATION_BUGS.md` and the linearization redesign.
**Branch target:** `sanghavi-unified` (with `unified-vsmartmom` as the merge base for the architecture, physics layered in from `sanghavi`).
**Reference codes consulted:** VLIDORT 2.8.3 (Spurr et al., RT Solutions Inc.) for established BRDF Fourier-decomposition and SS-correction architecture; CanopyOptics.jl for the existing hotspot abstraction; the current `sanghavi-unified` branch of vSmartMOM.jl for present-state code.
**Open for discussion:** every formula, threshold, and naming choice below — particularly the Cox-Munk `m_max` formula, which differs from the prior draft.

---

## TL;DR

Three related work items, one design vocabulary:

1. **Per-component `m_max` trait** — replace the surface-blind aerosol-only `max_m_bands` derivation with a dispatched trait that asks each component for its own Fourier-order ceiling. Following VLIDORT, the BRDF kernels (Cox-Munk, RPV, Ross-Li) all saturate at `l_trunc`. This corrects the prior draft's `2·l_trunc - 1`.
2. **`AbstractSSCorrection` hierarchy** — promote the `apply_ss_correction!` hack inside `coxmunk_surface.jl` to a first-class abstract type in the new `FourierDecomposition` module. Default `NoSSCorrection`; concrete `BRDFSpecularCorrection`, `AerosolForwardPeakCorrection`, `ExternalSSCorrection`. One unified call site replaces three currently-inconsistent ones across `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`.
3. **Wire `AbstractHotSpot` into `CanopySurface`** — CanopyOptics already ships `AbstractHotSpot`/`NoHotSpot`/`KuuskHotSpot`, but vSmartMOM's canopy path uses none of it. Parallel work item (different seam in the pipeline, same dispatch idiom).

The unifying design pattern:

> Identify a physical effect silently truncated or approximated by the solver's main path. Lift the correction to an abstract type with a `No···` default. Let `model_from_parameters` pick a sensible default per band/per surface. Let the user override. Have one apply hook called from one place.

This pattern is going to keep showing up — δ-M aerosol forward-peak correction, multiple-scatter δ-M correction, polarization SS approximation for thermal sources, Raman SS corrections. Each instance gets its own `Abstract*Correction` type living in the framework where its physics belongs.

---

## 1. Number of Fourier moments for complex surfaces

### 1.1 Three knobs, often conflated

The prior refactor draft talked about a single `m_max`. There are actually three independent quantities. Conflating them is what produced the `2·l_trunc - 1` formula in the draft, which I think is wrong by a factor of two.

| Quantity | What it is | Where it's set today | Notes |
|---|---|---|---|
| `Nquad` | Polar discrete-ordinate count (μᵢ on `[0,1]`) | `rt_set_streams.jl:32`: `Nquad = (l_trunc+1) ÷ 2` | **Already correct.** This is the user's quoted invariant. |
| `nQuad_BRDF` | Azimuth integration grid for computing each Fourier coefficient | Hard-coded `nQuad = 100` in `coxmunk_surface.jl` and `rpv_surface.jl` | Quadrature accuracy, not a physical limit. |
| `m_max` | Number of Fourier orders the solver loops over (`m ∈ 0:m_max`) | `model_from_parameters.jl:263`: `max_m_bands[i] = ⌈(l_max+1)/2⌉` from aerosol Greek coefs **only** | Surface is ignored. This is the gap. |

Three separate concerns. The trait being added in Step 1 of the prior draft is `m_max`, not the others. Worth making this explicit in the module documentation so future contributors don't merge them.

### 1.2 What VLIDORT does, and why I think it's the right reference

I checked the VLIDORT 2.8.3 vector BRDF supplement (`vbrdf_sup_masters.f90`, `vbrdf_sup_routines.f90`) because that codebase has 20+ years of working answers to exactly this question. Two findings:

**On how many moments to compute.** From `vbrdf_sup_masters.f90:2305`:

```fortran
NMOMENTS = 2 * NSTREAMS - 1
```

That's it. No surface-specific formula. No per-kernel exception. *All 17 non-Lambertian kernels* (Cox-Munk, RPV, Ross-Li, Hapke, Ross-Thick-Li-Sparse, Mod-Fresnel, etc.) get the same upper bound, set by the polar stream count.

**Rationale:** the BRDF Fourier series gets projected onto the discrete-ordinate basis during the RT solve. A polar grid with `NSTREAMS` points cannot resolve azimuthal modes higher than `2·NSTREAMS - 1` regardless of how rich the kernel itself is. Anything past that aliases. So the cap is set by *the solver's resolving power*, not by the surface's intrinsic angular complexity.

Translated to vSmartMOM's vocabulary, since `Nquad = (l_trunc+1)÷2`:

```
m_max_BRDF = 2·Nquad - 1 = 2·(l_trunc+1)/2 - 1 = l_trunc
```

For even `l_trunc` this is exact; for odd `l_trunc` the integer division shaves one mode (which is fine — `m_max` should round down, not up, to stay within the resolvable basis).

**This differs from the prior draft's `m_max(::CoxMunk) = 2·l_trunc - 1`.** I think the draft conflated `l_trunc` with `Nquad` — the formula `2·Nquad - 1` is correct, but the substitution into `l_trunc` was missed. **Worth confirming with Sanghavi** that vSmartMOM's `l_trunc` semantically corresponds to VLIDORT's "highest Legendre order in the truncated phase function expansion" rather than to `NSTREAMS`. If it does, then `m_max(::CoxMunk) = l_trunc` is right; if there's an internal convention I'm missing where `l_trunc` is defined as `Nquad`, then the prior draft was right. Either way, the framework should pin this down with a comment that derives the formula from VLIDORT's invariant rather than asserting it.

**On how each moment is computed.** VLIDORT precomputes the BRDF on a 3-D grid `BRDFUNC[stokes, μᵢ, μⱼ, kϕ]` over a Gauss-Legendre azimuth quadrature `X_BRDF` of size `NSTREAMS_BRDF`, and for each `m` builds weighted sums:

```fortran
BRDF_COSAZMFAC(I) = A_BRDF(I) * COS(M * X_BRDF(I))
L_BRDF_F(OM,I,J) = HELP * DOT_PRODUCT(BRDFUNC(OM,I,J,1:NK), BRDF_COSAZMFAC(1:NK))
```

`NSTREAMS_BRDF` defaults to **100** in their test suite — far higher than `NSTREAMS` itself. The polar grid is for the RT solve; the azimuth grid is for resolving glint specularities in the BRDF integration. They are decoupled by design.

vSmartMOM is doing the same thing already — `nQuad = 100` Gauss-Legendre on `[0, π]` in `_fourier_coeff_element`. The infrastructure is correct in kind. What's missing is making it explicit and unified across surfaces, instead of being hard-coded inside each surface file.

### 1.3 Why this is also the right answer for Cox-Munk specifically (the user's physics question)

For Cox-Munk you cannot get a closed-form Fourier expansion: the slope distribution couples non-trivially with the Fresnel matrix and an azimuth-dependent shadowing term, and there's no clean recursion the way Rayleigh has (`m_max(Rayleigh) = 2`).

RPV and Ross-Li are slightly better behaved analytically — the Rahman-Pinty-Verstraete kernel and the Ross-Li volumetric/geometric kernels have known analytic forms in `cos ϕ`, `cos 2ϕ` pieces — but in practice every production code (VLIDORT, 6SV1, MYSTIC) integrates them numerically anyway. Reasons:

- The cost is negligible compared to the RT solve.
- It unifies the code path across surfaces.
- It lets you swap kernels without rederiving Fourier coefficients.

So vSmartMOM's existing approach of "same numerical-quadrature recipe per surface" is correct. The fix is just to expose it cleanly.

### 1.4 What's already half-built in vSmartMOM

| Piece | Status | Location |
|---|---|---|
| Polar stream rule `Nquad = (l_trunc+1)÷2` | Correct | `rt_set_streams.jl:32` |
| Per-mode azimuth integration for Cox-Munk | Correct, `nQuad = 100` hard-coded | `coxmunk_surface.jl:_fourier_coeff_element` |
| Per-mode azimuth integration for RPV | Correct, `nQuad = 100` hard-coded | `rpv_surface.jl:reflectance` |
| Canopy: `compute_Z_matrices(canopy, μ, 0:m_max)` | Cleanly parameterized | CanopyOptics `z_matrices.jl` |
| Per-band `max_m_bands` from aerosol Greek coefs | Computed but **surface-blind** | `model_from_parameters.jl:254-263` |
| `m_max(::Surface)` dispatch trait | Doesn't exist yet | — |
| Canopy hotspot wired into elemental | **Missing** (uses independent Beer-Lambert) | `elemental_canopy.jl` (no hotspot calls) |

The Fourier-moment gap is small: just (6) — `model_from_parameters` doesn't ask the surface what *it* needs, so a Lambertian surface and a Cox-Munk surface look identical to the Fourier loop driver. The proposed `FourierDecomposition` module fixes this in roughly 12 lines of trait definitions.

### 1.5 Proposed `m_max` trait methods

```julia
# CoreRT/FourierDecomposition/m_max_traits.jl

"""
    m_max(component, l_trunc) -> Int

Per-component upper bound on the Fourier order m needed for an accurate
RT solve given a Legendre truncation `l_trunc`. The framework caps the
aggregate at `l_trunc` (no point computing more modes than the polar
discrete-ordinate basis can resolve, per VLIDORT's NMOMENTS = 2·NSTREAMS - 1
identity, which evaluates to l_trunc when Nquad = (l_trunc+1)÷2).
"""
m_max(::Rayleigh,                    l_trunc::Int) = 2

m_max(::LambertianSurfaceScalar,     l_trunc::Int) = 0
m_max(::LambertianSurfaceSpectrum,   l_trunc::Int) = 0
m_max(::LambertianSurfaceLegendre,   l_trunc::Int) = 0
m_max(::LambertianSurfaceSpline,     l_trunc::Int) = 0

# Following VLIDORT vbrdf_sup_masters.f90:2305 — NMOMENTS = 2·NSTREAMS - 1.
# With Nquad = (l_trunc+1)÷2, this evaluates to l_trunc.
m_max(::CoxMunkSurface,              l_trunc::Int) = l_trunc
m_max(::rpvSurfaceScalar,            l_trunc::Int) = l_trunc
m_max(::RossLiSurfaceScalar,         l_trunc::Int) = l_trunc

# Provisional placeholder. Proper LAD-driven derivation lives in CanopyOptics
# and depends on the leaf angle distribution; for now match the BRDF default.
m_max(::CanopySurface,               l_trunc::Int) = l_trunc

# From Greek-coef tail length, capped at l_trunc.
m_max(a::AerosolOptics,              l_trunc::Int) =
    min(length(a.greek_coefs.β) - 1, l_trunc)
```

Aggregator:

```julia
"""
    m_max_bands(rayleigh, aerosol_optics, surfaces, l_trunc) -> Vector{Int}

Per-band Fourier loop bound. Each band picks the maximum m_max across all
components present in that band, capped at l_trunc.
"""
function m_max_bands(rayleigh, aerosol_optics, surfaces, l_trunc::Int)
    n_bands = length(aerosol_optics)
    out = Vector{Int}(undef, n_bands)
    cap = cap_from_l_trunc(l_trunc)
    for b in 1:n_bands
        m = m_max(rayleigh, l_trunc)
        for a in aerosol_optics[b]
            m = max(m, m_max(a, l_trunc))
        end
        m = max(m, m_max(surfaces[b], l_trunc))
        out[b] = min(m, cap)
    end
    out
end

cap_from_l_trunc(l_trunc::Int) = l_trunc   # was (l_trunc+1)÷2 in the prior draft
```

**Two changes vs. the prior draft worth flagging:**

- **Cap formula.** Draft had `cap_from_l_trunc(l_trunc) = (l_trunc + 1) ÷ 2`. That's `Nquad`, not `m_max`. The correct cap is `l_trunc` itself per VLIDORT's invariant. With the draft's cap, `m_max` saturates at half the value VLIDORT uses, which under-resolves the BRDF.
- **Cox-Munk formula.** Draft had `m_max(::CoxMunk) = 2·l_trunc - 1`. Same conflation — should be `l_trunc`.

If I'm reading the draft's intent wrong and `l_trunc` *is* meant to be `Nquad` in vSmartMOM's internal convention, then both formulas are correct as written and *my* derivation is wrong. This is the single most important question to settle with Sanghavi before any code lands.

### 1.6 Open question on `nQuad_BRDF`

This is the *integration* grid for `∫₀^π ρ(μᵢ, μⱼ, ϕ) cos(mϕ) dϕ`, not the loop bound. Two options:

1. **Leave as-is** (pragmatic). 100-point Gauss-Legendre on `[0, π]` has spectral convergence; for any physically reasonable BRDF — including Cox-Munk at low wind speed (the worst case, sharp glint specular peak) — the integration error is below numerical noise. VLIDORT also defaults to 100. This works.
2. **Promote it to an `Int` field on `AbstractSurfaceType`** with default 100. Per-surface, not per-solver — it's a property of the BRDF's roughness, not of the solver. CoxMunk at wind speed 0.5 m/s might genuinely need 200 points; CoxMunk at 15 m/s is fine with 50.

Recommendation: option 2, but as Step 1.5 deferred. Doesn't block anything, low cost (~3 lines per surface), gives users a knob, and makes convergence testable.

### 1.7 Knock-on simplifications

Once the trait lands:

- **`apply_ss_correction!` argument list shortens** — currently takes `max_m::Int` positionally; can call `m_max(surf, l_trunc)` itself.
- **Lambertian `if m == 0` short-circuit (Step 2 of prior draft) becomes structurally unreachable.** With `m_max(::LambertianSurface*) = 0`, the loop simply doesn't enter `m > 0` for Lambertian-only bands. Net `−40 / +0` per the draft's estimate.
- **`vSmartMOM_Parameters.max_m`** can be deprecated to a "user override" only consulted if explicitly set, with `m_max_bands` as the default. Or dropped entirely as the draft proposes — the per-band `Vector{Int}` is strictly more informative.

---

## 2. Single-scattering correction as a first-class abstraction

### 2.1 The smell

Current state in `rt_run.jl:310-315`:

```julia
# Single-scattering correction for Cox-Munk specular hotspot (TMS)
if brdf isa CoxMunkSurface && SFI
    @timeit "SS Correction" apply_ss_correction!(
        R_SFI, brdf, pol_type, vza, vaz, μ₀,
        Array(τ_sum_all[:,end]), max_m, nSpec)
end
```

Three reasons this is wrong architecturally:

1. **Surface-specific dispatch via `isa`.** Adding a second BRDF that needs the same correction (e.g. fancier Cox-Munk-Bréon) means another `if` branch.
2. **Inconsistent across entry points.** `rt_run.jl` has it. `rt_run_lin.jl` does not — no linearized SS correction at all (silent Jacobian bug at glint geometries waiting to happen). `rt_run_multisensor.jl` does not — multisensor + Cox-Munk silently skips the correction.
3. **The `if brdf isa CoxMunkSurface` gate hides physics.** The same kind of correction is *also* needed for δ-M-truncated aerosol forward peaks (the Nakajima-Tanaka second-half scheme). vSmartMOM has the first half (`delta_m_truncation.jl`) but not the second. The current "SS correction" code was solving Cox-Munk specifically without anyone noticing the same pattern was missing for aerosols.

### 2.2 What VLIDORT does, again as the reference

VLIDORT factors single-scattering into a first-class subsystem called **FO — First-Order Code, exact single-scattering and direct-thermal**. It's a peer of the multiple-scatter solver, not a postscript to Cox-Munk. The control surface is a hierarchy of flags in `vlidort_inputs_def.f90:444-481`:

- **`DO_FOCORR`** — master on/off. Without it, "VLIDORT will perform a truncated pseudo-spherical SS calculation" — i.e. you accept the Fourier-truncation residual as your error budget.
- **`DO_FOCORR_EXTERNAL`** — external pre-computed exact-SS hook. Someone else does the FO calc; VLIDORT just adds it.
- **`DO_FOCORR_NADIR`** / **`DO_FOCORR_OUTGOING`** — sphericity treatment (mutually exclusive).
- **`DO_SSCORR_TRUNCATION`** — whether the SS correction also accounts for δ-M phase function truncation.
- **`DO_SSCORR_USEFMAT`** — evaluate the exact phase function from the F-matrix vs. from Greek coefficients. This is the choice between "compute exact SS from truncated Legendre series" vs. "from full phase matrix"; the latter is what gives you the glint hotspot back.

Lessons:

1. SS correction is a **subsystem with its own master flag**, not a per-surface afterthought.
2. It is **not surface-specific**. The same pass corrects truncated aerosol forward-peaks and truncated surface specularities.
3. There needs to be an **"external" hook** so users (or test harnesses, or downstream tools) can override the built-in correction with a pre-computed one.

### 2.3 What's already in place, and what's missing

- ✅ `apply_ss_correction!` in `coxmunk_surface.jl:505-546` does the Nakajima-Tanaka math correctly: `(M_exact − M_fourier) × beam_attenuation`. Numerically right, just architecturally trapped.
- ✅ `delta_m_forward` in `LayerOpticalProperties/delta_m_truncation.jl` does the **first half** of Nakajima-Tanaka for aerosols (modify τ and ϖ̃ to compensate for the truncated forward peak), with linearization (`delta_m_truncation_lin`).
- ❌ The **second half for aerosols** — replacing the truncated SS aerosol contribution with the exact phase function evaluation — is missing entirely. Aerosols currently get the truncation but not the correction. Quiet ~1% bias at high optical depth, much larger near the forward peak.
- ❌ **Linearized SS correction** missing throughout. Will matter for Jacobians once the `f₁` ForwardDiff redesign lands.
- ❌ **Multisensor SS correction** missing throughout.
- ⚠️ **Name collision**: `Scattering/types.jl:63` defines an `AbstractFourierDecompositionType` whose subtypes are `NAI2` and `PCW`. That's the *Mie phase function decomposition algorithm strategy* (Siewert vs Domke), not the per-mode RT loop the new module owns. Two different concepts sharing a name. Recommend renaming the older one to `AbstractMieDecompositionAlgorithm` before the new module lands.

### 2.4 Note on `interaction_ss.jl`

There's a file called `interaction_ss.jl` in `CoreKernel/`. **It is not a post-Fourier correction.** It's a per-layer kernel for the unscattered solar→sensor direct-beam contribution within the multiple-scatter solve. Different concept. Worth being explicit about the namespacing so the new abstraction doesn't get conflated with it during code review.

### 2.5 Proposed type hierarchy

```julia
# CoreRT/FourierDecomposition/ss_correction.jl

"""
    AbstractSSCorrection

Marker for a single-scattering correction strategy applied after the Fourier
loop completes. The Fourier expansion truncates contributions with sharp
azimuthal or angular features (BRDF specular peaks, aerosol forward peaks);
an SS correction restores them by replacing the truncated single-scatter
contribution with its exact evaluation at the actual viewing geometry.

Concrete subtypes implement `apply_ss_correction!` and may override
`needs_correction`.
"""
abstract type AbstractSSCorrection end

"""Default. Apply nothing — used when the truncation residual is acceptable."""
struct NoSSCorrection <: AbstractSSCorrection end

"""
    BRDFSpecularCorrection{S<:AbstractSurfaceType, FT}

Replaces the truncated Fourier reconstruction of the surface BRDF
single-scatter contribution with its exact evaluation. Nakajima-Tanaka TMS
applied to surfaces — what the current `apply_ss_correction!` in
coxmunk_surface.jl already does, but generic over surface type.

Activated by default whenever `has_specular_peak(surf) == true`.
"""
struct BRDFSpecularCorrection{S<:AbstractSurfaceType, FT} <: AbstractSSCorrection
    surface::S
    nQuad_BRDF::Int   # default 100, tunable per-surface
end

"""
    AerosolForwardPeakCorrection{FT}

The δ-M companion: replaces the truncated SS aerosol contribution with the
exact phase function evaluation. Currently missing from vSmartMOM. With this
in place, `delta_m_truncation.jl` becomes a complete Nakajima-Tanaka scheme
rather than just its first half.

Activated by default whenever any aerosol in the band has `fᵗ > threshold`.
"""
struct AerosolForwardPeakCorrection{FT} <: AbstractSSCorrection
    fᵗ_threshold::FT   # below this, don't bother
end

"""
    ExternalSSCorrection{F}

User-supplied correction. Stores a callable returning the additive
contribution. Equivalent to VLIDORT's `DO_FOCORR_EXTERNAL`. Useful for
testing, coupling, and validation.
"""
struct ExternalSSCorrection{F} <: AbstractSSCorrection
    callback::F
end
```

Two traits each subtype answers:

```julia
"""
    needs_correction(c::AbstractSSCorrection, model, iBand) -> Bool

Cheap predicate. Lets the loop driver skip composition+allocation when
nothing fires. Default: true; NoSSCorrection overrides to false.
"""
needs_correction(::NoSSCorrection,     model, iBand) = false
needs_correction(::AbstractSSCorrection, model, iBand) = true

"""
    apply_ss_correction!(R_SFI, c, model, iBand) -> nothing

In-place additive update of the post-Fourier reflectance. Each subtype owns
the math for its own correction; nothing is gated on `if surf isa ...`.
"""
function apply_ss_correction! end   # methods per concrete subtype
```

Per-surface trait controlling the default:

```julia
has_specular_peak(::AbstractSurfaceType) = false
has_specular_peak(::CoxMunkSurface)      = true
# Future: has_specular_peak(::ModFresnelSurface) = true, etc.
```

### 2.6 Composition over inheritance

A band may need *both* aerosol forward-peak and BRDF specular corrections. They're additive and physically independent. So `SolverConfig` carries a tuple per band:

```julia
ss_corrections::Vector{Tuple{Vararg{AbstractSSCorrection}}}
```

with the loop:

```julia
function apply_ss_corrections!(R_SFI, model, iBand)
    for c in model.solver.ss_corrections[iBand]
        needs_correction(c, model, iBand) || continue
        apply_ss_correction!(R_SFI, c, model, iBand)
    end
    return R_SFI
end
```

Called *once* from each `rt_run*.jl`, after the Fourier loop. No more `if brdf isa CoxMunkSurface`. Multisensor gets the call for free; linearized version gets it (with lin methods dispatched on the same types) for free.

**Why tuple, not Vector{AbstractSSCorrection}?** Type stability. The set of corrections is fixed at model-construction time; using a `Tuple` lets the compiler specialize the loop. `Vector{AbstractSSCorrection}` would force runtime dispatch in the inner loop, which matters for the Cox-Munk apply specifically (it's already inside an `iv × max_m × n × n` quadruple loop).

### 2.7 Default selection

```julia
function default_ss_corrections(surfaces, aerosol_optics, l_trunc)
    n_bands = length(aerosol_optics)
    out = Vector{Tuple{Vararg{AbstractSSCorrection}}}(undef, n_bands)
    for b in 1:n_bands
        cs = AbstractSSCorrection[]
        if has_specular_peak(surfaces[b])
            push!(cs, BRDFSpecularCorrection(surfaces[b], 100))
        end
        if any(a -> a.fᵗ > 0.01, aerosol_optics[b])
            push!(cs, AerosolForwardPeakCorrection(0.01))
        end
        out[b] = isempty(cs) ? (NoSSCorrection(),) : Tuple(cs)
    end
    out
end
```

User can override the whole vector — including supplying `ExternalSSCorrection(my_callback)` for testing or for coupling to an external SS solver.

### 2.8 What this gives you for the linearization redesign

This is the second-order benefit and worth flagging now.

The SS correction is one of the awkward pieces for the `f₁` ForwardDiff redesign because:
1. It depends on the **exact** BRDF / phase function, not the truncated Greek-coefficient form.
2. It depends on **viewing geometry** `(μ₀, μ_v, Δϕ)`, not just on quadrature points.
3. It produces an **additive contribution** to the radiance that needs its own derivative chain.

If `BRDFSpecularCorrection` and `AerosolForwardPeakCorrection` are first-class types with a clean `apply_ss_correction!` interface, the linearization story becomes:

> Wrap each correction's apply method in a ForwardDiff seed → get `∂(correction)/∂(Mie params, surface params)` for free → add to the existing chain rule via `lin_added_layer_all_params.jl`.

The current implementation makes that hard because the SS correction hack happens *outside* the `lin_added_layer_all_params.jl` seam. There's no clean way to ForwardDiff through `if brdf isa CoxMunkSurface`. With the abstraction, each correction is its own self-contained unit of differentiation. Likely earns its own item in `LINEARIZATION_BUGS.md` if not addressed.

### 2.9 Optimization deferred: cache vs. recompute

Current Cox-Munk `apply_ss_correction!` reconstructs `M_fourier` by *re-evaluating* the Fourier coefficients on the fly inside `_fourier_coeff_element`. Wasteful — those coefficients were already computed inside the Fourier loop. Caching would speed the correction up substantially, but couples the loop and the correction.

Recommendation: **keep recomputing** for the first landing. Abstraction first, optimization second. Cost is negligible compared to the RT solve. Add caching as a follow-up once the structure is in place. **Worth confirming with Sanghavi.**

---

## 3. Canopy hotspot — same idiom, different seam

### 3.1 The temptation, and why it's wrong

When I first thought about the canopy hotspot fitting into `AbstractSSCorrection`, the temptation was to add `KuuskHotSpotCorrection <: AbstractSSCorrection` to the tuple and call it done. **It would be wrong.** Worth being explicit about why, because the same temptation will recur for other corrections (Raman SS, thermal SS, etc.) and the answer is the same.

The two corrections live on **completely different seams** in the RT pipeline:

| | `AbstractSSCorrection` | `AbstractHotSpot` |
|---|---|---|
| **When** | After Fourier loop completes | Inside `elemental_canopy!` source assembly |
| **Where** | Post-processing on `R_SFI` | Per-layer · per-Fourier-mode |
| **How** | Additive: `(M_exact − M_fourier)` per geometry | Multiplicative: `C_hs · P_so/P_o` on joint gap probability |
| **What it patches** | Fourier truncation losing specular-peak resolution | Beer-Lambert assuming independent gaps |
| **Where the math lives** | `vSmartMOM/CoreRT/FourierDecomposition` | `CanopyOptics/canopy_structure/hotspot.jl` |

If I tried to force the hotspot into `AbstractSSCorrection`, the apply hook signature would have to widen to either:
- (a) handle the in-loop case (which means it's no longer a post-Fourier correction — abstraction is leaking), or
- (b) operate on `R_SFI` after the fact and reverse-engineer what should have happened inside the source assembly (which is a hack of the kind we're trying to escape).

The conceptual similarity ("both correct truncations the RT solver makes for tractability") is real, and that similarity is captured by the **shared dispatch idiom** — `Abstract*Correction → No··· · Concrete··· · External···` — not by a shared type.

### 3.2 What CanopyOptics already ships

`CanopyOptics.jl/src/canopy_structure/hotspot.jl` defines:

- `abstract type AbstractHotSpot{FT<:Real}`
- `struct NoHotSpot{FT}` — independent gaps, `C_hs = 1`
- `struct KuuskHotSpot{FT}` with size parameter `h` (≈ leaf size / canopy height)
- `canopy_extinction(G_value, μ) = G_value / μ`
- `hotspot_separation(μ_s, μ_o, dϕ, h)` — the angular separation parameter
- `hotspot_correction(model, k_s, k_o, μ_s, μ_o, dϕ, L)` — the multiplicative correction
- `joint_gap_probability(model, k_s, k_o, μ_s, μ_o, dϕ, L)` — `exp(-(k_s+k_o)L) · C_hs(L)`

The author was thinking ahead about the linearization redesign too — there's a `_hotspot_value(x::ForwardDiff.Dual)` plumbing helper already in place. The Kuusk math is a closed-form 4-line expression in `(h, k_s, k_o, μ_s, μ_o, dϕ, L)`; Jacobians w.r.t. `h` (and through `k_s, k_o` to LAI/G) drop out cleanly under ForwardDiff.

The docstring on `joint_gap_probability` even explains how the consumer is supposed to wire it up. Quoting:

> Hotspot corrections belong in the direct-beam source assembly of the canopy RT consumer, not in compute_Z_matrices.

So the abstraction is shipped, the implementation is correct, the linearization story is anticipated. **What's missing is the wiring on the vSmartMOM side.**

### 3.3 What's missing in vSmartMOM

I grepped `src/CoreRT/Surfaces/canopy_surface.jl`, `src/CoreRT/CoreKernel/elemental_canopy.jl`, and `src/CoreRT/CoreKernel/interaction_hdrf.jl` for any reference to `hotspot`, `Hotspot`, `joint_gap`, `kuusk`, `nohot`. **Zero matches.**

`elemental_canopy!` uses `τ_sum` (cumulative LAI) for direct-beam attenuation via `arr_type(τ_sum)` multiplication — that's the **independent-gap** case (`NoHotSpot` implicit). There's no path for the user to opt in to `KuuskHotSpot` because `CanopySurface` doesn't carry a hotspot field, and `joint_gap_probability` is never called.

**Implication:** any vSmartMOM canopy result currently being compared against literature with hotspot effects (RAMI test cases, NDVI/red-edge over crop canopies, MODIS BRDF retrievals at small phase angle) is wrong by whatever `C_hs` would have been. **This is a correctness fix, not a refactor.**

⚠️ Worth grepping the vSmartMOM test suite before quietly fixing this — if any reference benchmark was tuned to match the current independent-gap result, that benchmark needs to be regenerated with the corrected physics, and the discrepancy needs to be documented.

### 3.4 Proposed wiring

```julia
# CoreRT/types.jl — add hotspot field to CanopySurface

mutable struct CanopySurface{FT, ...} <: AbstractSurfaceType
    # ... existing fields ...
    hotspot::AbstractHotSpot{FT}   # default NoHotSpot{FT}() to preserve current behavior
end
```

```julia
# CoreRT/CoreKernel/elemental_canopy.jl — replace implicit Beer-Lambert

# Today (implicit independent gaps):
#   direct_source uses τ_sum directly via arr_type(τ_sum)

# Proposed:
P_so = joint_gap_probability(canopy.hotspot, k_s, k_o, μ_s, μ_o, dϕ, τ_sum)
P_o  = exp(-k_o * τ_sum)
direct_source_scale = P_so / P_o
# multiply into the direct-beam source term before the canopy scattering kernel
```

With `NoHotSpot()` the result is bit-for-bit identical to today. With `KuuskHotSpot()` the user gets the corrected joint gap.

Same wiring needed in `interaction_hdrf_canopy!` for the HDRF/BHR computation — the hotspot affects directional-hemispherical reflectance too, since it changes how much direct beam couples to the view path.

### 3.5 Validation

Reference benchmark against published canopy hotspot test cases:

- **RAMI** (Radiation transfer Model Intercomparison) — has published results for both NoHotSpot and KuuskHotSpot configurations across canopy depth and sun-view angle separation.
- **PROSAIL-RT** — turbid-medium reference.
- Specific test scene: HOM26 / HOM27 from RAMI-IV (homogeneous canopy with explicit hotspot tests).

This is the step where reference benchmarks are non-negotiable, similar to the other "physics actually changes" steps elsewhere in the v2.0.0 plan.

---

## 4. The unifying design pattern

The three work items above share a single architectural pattern that I think should keep being applied:

> **Identify a physical effect silently truncated or approximated by the solver's main path. Lift the correction to an abstract type with a `No···` default. Let `model_from_parameters` pick a sensible default per band/per surface. Let the user override. Have one apply hook called from one place.**

Instances of this pattern that exist or are coming:

| Correction | Abstract type | Status | Lives in |
|---|---|---|---|
| BRDF specular peak (Cox-Munk glint) | `AbstractSSCorrection`/`BRDFSpecularCorrection` | Migrating from hack | `FourierDecomposition` |
| Aerosol forward peak (δ-M companion) | `AbstractSSCorrection`/`AerosolForwardPeakCorrection` | **Missing entirely** | `FourierDecomposition` |
| Canopy hotspot | `AbstractHotSpot`/`KuuskHotSpot` | Shipped in CanopyOptics, **not wired** | `CanopyOptics` |
| User-supplied SS | `ExternalSSCorrection` | Proposed | `FourierDecomposition` |
| δ-M phase truncation | `AbstractTruncationType`/`δBGE`/`NoTruncation` | Already exists | `Scattering` (probable migration to `FourierDecomposition`) |
| Polarization SS approx for thermal | TBD | Future | TBD |
| Raman SS correction | TBD | Future | `Inelastic` |

Pattern reuse means consistent vocabulary for users, consistent dispatch idiom for contributors, consistent place to attach linearization hooks. Worth being explicit about it in the module-level docstrings as the framework grows.

**Naming consistency note for review:** the existing `AbstractTruncationType` (with `δBGE` and `NoTruncation`) follows the same idiom but with different naming (`No···Type` vs `No···Correction`). Worth deciding if these get harmonized during the FourierDecomposition migration that's already on the table for that type.

---

## 5. Combined plan for Claude Code

Step ordering picks pure refactors first, then physics changes, isolating the points where reference benchmarks are mandatory.

### Phase 1: Plumbing (no behavior change)

1. **Rename collision.** `Scattering/types.jl::AbstractFourierDecompositionType` → `AbstractMieDecompositionAlgorithm` (with `NAI2`, `PCW` as subtypes). ~10 lines, no behavior change. Frees up `FourierDecomposition` as an unambiguous module name. All tests pass.

2. **Step 0.5 (per prior draft):** introduce `FourierDecomposition` module skeleton at `src/CoreRT/FourierDecomposition/FourierDecomposition.jl`. Includes, exports, module declaration. ~30 lines.

3. **Step 1 (per prior draft):** add `m_max(component, l_trunc)` trait and methods inside the module. ~12 lines of trait definitions. Add unit tests:
   - `m_max(LambertianSurfaceScalar(0.3), 60) == 0`
   - `m_max(CoxMunkSurface(...), 60) == 60`
   - `m_max_bands` aggregates correctly across components

4. **Refactor BRDF Fourier integration into shared helper.** Move the `nQuad = 100` Gauss-Legendre integration in `coxmunk_surface.jl` and `rpv_surface.jl` into a single helper `_brdf_fourier_coefficient(surf, n_stokes, μ_grid, m; nQuad)` inside `FourierDecomposition`. Pure refactor, exact numerical equivalence verifiable by existing tests.

5. **Add `AbstractSSCorrection` hierarchy.** `NoSSCorrection`, `BRDFSpecularCorrection`, `AerosolForwardPeakCorrection`, `ExternalSSCorrection`, plus `apply_ss_correction!` dispatch + `needs_correction` + `has_specular_peak`. No call sites change yet.

6. **Migrate Cox-Munk `apply_ss_correction!`** from `coxmunk_surface.jl` to `FourierDecomposition/corrections/brdf_specular.jl` as the `BRDFSpecularCorrection{<:CoxMunkSurface}` method. Pure relocation; numerics identical. Existing tests pass.

### Phase 2: Wiring (potential behavior change at known points)

7. **Wire `m_max_bands` into `model_from_parameters.jl` and `lin_model_from_parameters.jl`.** This is the moment a Cox-Munk band's behavior could change because it'll now compute more Fourier orders than before. **Run 6SV1 / Natraj reference benchmarks before and after.** If previous default `params.max_m` was set high enough by users, nothing changes. If too low, results get *better* against reference.

8. **Wire `apply_ss_corrections!` into `rt_run.jl`.** Replace the `if brdf isa CoxMunkSurface` block with `apply_ss_corrections!(R_SFI, model, iBand)`. Ensure `model_from_parameters` populates `ss_corrections` via `default_ss_corrections`. Existing benchmarks unchanged for CoxMunk.

9. **Wire into `rt_run_multisensor.jl`.** One new call site. **Net behavior change** for multisensor + Cox-Munk (gains the correction it was silently missing). Document in changelog; verify against any existing multisensor reference cases.

10. **Wire into `rt_run_lin.jl`.** Same one-liner. Lin methods come later; for the forward run it just calls the same apply hook.

11. **Delete Lambertian `if m == 0` short-circuits** per Step 2 of prior draft. Now structurally unreachable. Net `−40 / +0`.

### Phase 3: New physics

12. **Add `AerosolForwardPeakCorrection` implementation.** Genuinely new physics in vSmartMOM — the missing companion to δ-M truncation. **Validate against PyVLIDORT or 6SV1** with high-AOD strongly-forward-peaked aerosols (large dust particles are the obvious test case). Reference benchmarks are non-negotiable.

13. **Wire `AbstractHotSpot` into `CanopySurface`.** Add `hotspot` field defaulting to `NoHotSpot()`. Replace implicit Beer-Lambert in `elemental_canopy!` and `interaction_hdrf_canopy!` with `joint_gap_probability` calls. With `NoHotSpot()` (default), bit-for-bit identical to current behavior.

14. **Validate canopy hotspot against RAMI.** Use HOM26/HOM27 or similar published canopy hotspot test cases. Document any vSmartMOM test scenes whose previous numbers depended on the missing-hotspot bug.

### Phase 4: Optional / deferred

15. **Promote `nQuad_BRDF` to per-surface field** with default 100. ~3 lines per surface. Gives users a knob for very-rough surfaces.

16. **Cache Fourier coefficients** for SS correction reuse (vs. recomputing in `_fourier_coeff_element`). Optimization-only; structure must already be in place.

17. **Linearized SS correction methods** (`apply_ss_correction!` for `Lin*` types). Falls out of the abstraction once the `f₁` ForwardDiff redesign lands.

### Phase boundaries for benchmarks

- **End of Phase 1:** All existing tests must pass bit-for-bit. Pure plumbing.
- **After step 7:** 6SV1 / Natraj benchmarks for any band where `max_m` could change.
- **After step 9:** any existing multisensor + Cox-Munk reference cases.
- **After step 12:** PyVLIDORT or 6SV1 high-AOD strongly-forward-peaked aerosol cases.
- **After step 13:** RAMI HOM26/HOM27 canopy hotspot cases.

---

## 6. Open questions for Sanghavi

In rough order of "blocks the next commit" → "nice to settle eventually":

1. **`l_trunc` semantics.** Does vSmartMOM's `l_trunc` correspond to VLIDORT's "highest Legendre order in the truncated phase function expansion" (in which case `m_max(::CoxMunk) = l_trunc`) or to VLIDORT's `NSTREAMS` (in which case the prior draft's `2·l_trunc - 1` was right)? **Blocks step 3.** Easiest answer: derive both formulas from `Nquad = (l_trunc+1)÷2` and the VLIDORT invariant `m_max = 2·Nquad - 1`, see which one matches the current numerics on a single Cox-Munk + Rayleigh test.

2. **Cap formula.** Related to (1). The prior draft's `cap_from_l_trunc(l_trunc) = (l_trunc+1)÷2` looks like `Nquad`, not the cap. Is this a typo or a deliberate convention difference?

3. **Cox-Munk SS correction caching.** Keep recomputing Fourier coefficients in `_fourier_coeff_element` (simple, safe, negligible cost) or cache from the main Fourier loop (faster, couples loop and correction)? Recommend deferring caching but want confirmation.

4. **`AerosolForwardPeakCorrection` threshold.** Default `fᵗ > 0.01` to skip the correction when truncation is negligible. Reasonable? Worth tuning during validation.

5. **Naming harmonization.** Should existing `AbstractTruncationType`/`NoTruncation` be renamed to fit the proposed `Abstract*Correction`/`No*Correction` idiom? Or keep them as historical exceptions? Affects how `δBGE` migrates from `Scattering` to `FourierDecomposition`.

6. **Canopy hotspot regression risk.** Is there any vSmartMOM test or reference benchmark whose published numbers would shift when `KuuskHotSpot` becomes wireable? If so, those need regenerating with documentation, not silent fixing.

7. **`m_max(::CanopySurface)` provisional vs. proper.** Placeholder is `l_trunc`. The "proper" LAD-driven value lives in CanopyOptics. Is it worth doing the proper derivation now or deferring? I'd defer — the placeholder doesn't change current numerics.

8. **Module location for `δBGE` and `NoTruncation`.** Prior draft flagged these as "probable second migration" from `Scattering/` to `FourierDecomposition/`. Worth confirming the move happens during this work or stays a separate followup.

---

## 7. References

- **VLIDORT 2.8.3 source**, particularly:
  - `vsup/vbrdf/vbrdf_sup_routines.f90` — `VBRDF_FOURIER` routine (numerical Fourier integration of BRDF)
  - `vsup/vbrdf/vbrdf_sup_masters.f90:2305` — `NMOMENTS = 2 * NSTREAMS - 1` invariant
  - `vlidort_focode/VFO_Master.f90` — first-order (exact SS) subsystem
  - `vlidort_def/vlidort_inputs_def.f90:444-481` — `DO_FOCORR*` flag hierarchy
- **Nakajima & Tanaka 1988**, *Algorithms for radiative intensity calculations in moderately thick atmospheres using a truncation approximation*, JQSRT 40(1), 51-69. Original δ-M / TMS scheme.
- **CanopyOptics.jl** `src/canopy_structure/hotspot.jl` — `AbstractHotSpot`, `KuuskHotSpot`.
- **Kuusk 1991**, *The hot spot effect in plant canopy reflectance* (referenced in CanopyOptics docstrings). Original Kuusk hotspot model.
- **RAMI** (Radiation transfer Model Intercomparison) test scenes for canopy hotspot validation.
- **Spurr 2002**, *VLIDORT: A linearized pseudo-spherical vector discrete ordinate radiative transfer code* — original VLIDORT linearization approach, relevant for the `f₁` ForwardDiff redesign.

---

## 8. Document status / iteration log

- **v0.1** — initial draft after VLIDORT/CanopyOptics/vSmartMOM code review. Open questions (Section 6) need Sanghavi's input before Phase 1 step 3 commits.
- *(future entries here as the design firms up)*
