# Standalone exact single-scattering kernel solver — design

**Document version:** v0.4
**Status:** focused implementation plan; standalone, can be built before the bigger refactor.

**Document set (delivered together):**
- `vsmartmom_dispatch_design_v0_6.md` — architecture this eventually integrates with.
- `standalone_ss_solver_plan.md` — *this document* (Piece A).
- `vlidort_baseline_suite.md` — VLIDORT golden-standard validation (Piece B; independent track).
- `dev_notes/exact_ss_reference/` — standalone four-path reference (Julia + Python + derivations).

**Pre-existing committed artifacts referenced:**
- `docs/dev_notes/LINEARIZATION_BUGS.md` — linearization bug catalogue.
- `src/CoreRT/types.jl` line 1052 — `CoreScatteringOpticalProperties{FT,FT2,FT3}` with `τ::FT, ϖ::FT2, Z⁺⁺::FT3, Z⁻⁺::FT3`. *This is the existing pol-type-generic core abstraction.*
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl` line 32 — "Bug 19" note: Z chain rule must be applied before doubling-adding mixes indices. *This is the design precedent for our chain-rule structure in §6.*

### Changes from v0.3

This version addresses GPT review findings on v0.3:

- **§3 (data structures): pol-type-generic from day one**, aligned with existing `CoreScatteringOpticalProperties`. Z⁺⁺/Z⁻⁺ matrix representation, n_stokes parametric. Phase 1 implementation specializes for n_stokes=1 but the framework is generic.
- **§5 (dual-form interface): align with v0.6 §5.5 / v0.5 §6.3** — `TruncatedAndExactScatteringOpticalProperties` carries *both forms as full triplets* `(τ, ϖ, Z⁺⁺, Z⁻⁺)`. No mixing of primed and unprimed within a single form.
- **§6 (linearization): chain rule contracts properly**, no element-wise broadcasting across optical-property axes. Cites `rt_kernel_lin.jl:32` Bug 19 precedent.
- **§6 (path 1 derivative): sign error fixed** — same-layer derivative is `+exp(-τ_bot·a)·(ω·P·I₀)/(4π·μ_v)`, positive after the `1/a` prefactor cancels with the `+a` from differentiation. Lower-layer attenuation contributions are negative (separate term).
- **§7 (back-correction): scope narrowed.** Back-correction is a **diagnostic for FO-equivalent comparison only** — adds at the sensor; does not propagate the corrected source through `(I - M_trunc)^-1`. The full-MOM correction is `ExactSFIPhase` (architecture doc §6), not this back-correction. Validation against VLIDORT FO (Tasks 1↔3 in solar_tester) uses the back-correction; validation against full FullMOM-with-FO uses `ExactSFIPhase`.
- **Cross-references**: distinguish committed artifacts (in `vSmartMOM.jl`) from delivered-with-this-document and to-be-built.

### Changes carried forward from v0.2/v0.3
- **§5 (required-Nquad framework)** — empirically derive Nquad from contributors and surface, both for standalone solver inner integrals and for full MOM stream count.
- **§6 (AD-handcoded linearization design experiment)** — standalone solver as testbed for hybrid AD-handcoded architecture.
- **§7 (back-correction)** — exists, now properly scoped.

---

## 0. What this is

A standalone Julia module implementing exact single-scattering through plane-parallel atmosphere. Self-contained, fast, ForwardDiff-friendly, with a hand-coded f₂ linearization at the core-optical-properties seam designed as a *testbed* for the broader vSmartMOM linearization redesign.

Designed so that:

1. **It runs by itself** — no dependency on existing `rt_run` machinery.
2. **It's fast** — kernel-based dispatch over `(viewing_geometry, spectral_point)`; CPU+GPU friendly.
3. **It auto-determines required Nquad** for both its own paths-3+4 inner quadrature and (as byproduct) the full MOM solver's stream count.
4. **It supports hybrid AD-handcoded linearization** at the core-optical-properties seam — testbed for broader lin redesign.
5. **It validates two ways**: against `dev_notes/exact_ss_reference/` (correctness oracle) and against VLIDORT shipped fixtures (Piece B).
6. **It's pol-type generic from day one** — Phase 1 specializes for n_stokes=1 but the framework matches the existing `CoreScatteringOpticalProperties{FT,FT2,FT3}` abstraction.
7. **It enables FO-equivalent diagnostic back-correction** of current vSmartMOM. Scope-bounded: validates against VLIDORT FO outputs (Task1↔Task3 difference); does not substitute for `ExactSFIPhase` (which is the production full-MOM correction, defined in v0.6 architecture doc §6).

---

## 1. Why this comes first

Three reasons:

1. **Concrete is better than abstract.** Working kernel makes the architecture "wrap this in dispatch."
2. **Independently useful** — fast tool for SS-only computations, retrieval inner loops.
3. **Three high-leverage experiments** before the refactor:
   - **Back-correction (§7)** — FO-equivalent diagnostic, no plumbing needed.
   - **Required-Nquad framework (§5)** — derives precheck thresholds empirically.
   - **Hybrid AD-handcoded linearization (§6)** — small-scale prototype for broader redesign in `LINEARIZATION_BUGS.md`.

---

## 2. Module layout

```
src/StandaloneSS/
├── StandaloneSS.jl                 # Module entry; exports public API
├── types.jl                        # ExactSSConfig; aligns with vSmartMOM.CoreRT types
├── kernels.jl                      # KernelAbstractions kernels (pol-type generic)
├── path1.jl                        # exact_ss for AtmosphericSSPath
├── path2.jl                        # exact_ss for SurfaceDirectBeamPath{<:LambertianSurface},
│                                   #          SurfaceDirectBeamPath{<:CoxMunkSurface}
├── path3.jl                        # exact_ss for AtmosphereToSurfacePath
├── path4.jl                        # exact_ss for SurfaceToAtmospherePath{<:Surface}
├── canopy_hotspot.jl               # exact_ss for CanopyHotspotPath
├── quadrature_required.jl          # Required-Nquad framework (§5)
├── linearized_f2.jl                # f₂ handcoded derivatives (§6); contraction-correct
├── chain_rule.jl                   # Combiner: handcoded f₂ ∘ AD f₁ — contraction (§6)
├── solver.jl                       # Top-level orchestration
├── back_correction.jl              # FO-equivalent diagnostic (§7) — scope-bounded
└── public_api.jl                   # Public-facing run functions

test/StandaloneSS/
├── runtests.jl
├── test_vs_reference.jl            # vs dev_notes/exact_ss_reference/
├── test_vs_vlidort.jl              # vs Piece B fixtures (Siewert 2000, solar_tester)
├── test_brdf_reciprocity.jl
├── test_thin_limit.jl
├── test_forwarddiff.jl
├── test_handcoded_vs_AD.jl         # Three-way validation (§6)
├── test_chain_rule_contraction.jl  # Contraction-correctness explicit test (§6)
├── test_quadrature_required.jl
├── test_back_correction.jl         # Diagnostic against VLIDORT FO scope (§7)
└── fixtures/
```

---

## 3. Data structures — pol-type-generic from day one

### 3.1 Aligning with existing `CoreScatteringOpticalProperties`

The existing type at `src/CoreRT/types.jl:1052`:

```julia
Base.@kwdef struct CoreScatteringOpticalProperties{FT,FT2,FT3} <: AbstractOpticalProperties
    τ::FT       # layer optical depth
    ϖ::FT2      # single-scattering albedo
    Z⁺⁺::FT3    # forward-scattering Z matrix (n_stokes × n_stokes × n_geometries × ...)
    Z⁻⁺::FT3    # backward-scattering Z matrix
end
```

This is already pol-type-generic — `Z⁺⁺` and `Z⁻⁺` are matrix arrays whose dimensions reflect `n_stokes`. The standalone solver consumes this type directly; no parallel "Greek-coefficient" representation.

The Z matrices represent the phase function evaluated at the actual geometry (forward and backward semi-spheres). For the SS computation at the actual viewing geometry, what we need is `Z⁺⁻(cosΘ, geometry)` — the Z matrix at the SS scattering angle. This is a property of the contributor (aerosol, Rayleigh, etc.) and computed via:

```julia
"""
    exact_phase_function(contributor, cosΘ, n_stokes) -> SMatrix{n_stokes, n_stokes}

Returns the F-matrix (4×4 for full polarization, 1×1 for scalar) at scattering angle cosΘ.
For aerosols, evaluated from Greek coefficients via Wigner-d functions.
For Rayleigh, the analytic 4×4 Rayleigh F-matrix.
"""
exact_phase_function(::AbstractContributor, cosΘ, n_stokes) end
```

The standalone solver's path-1 kernel uses `exact_phase_function` rather than reconstructing from Z⁺⁺/Z⁻⁺ at the discrete-ordinate angles, because the SS scattering angle isn't generally one of the quadrature angles.

### 3.2 Dual-form `TruncatedAndExactScatteringOpticalProperties`

Aligning with v0.6 architecture doc §5.5 and v0.5 §6.3:

```julia
"""
    TruncatedAndExactScatteringOpticalProperties{T, E, FT}

Both truncated and exact forms as full triplets. Each is a CoreScatteringOpticalProperties
with consistent (τ, ϖ, Z⁺⁺, Z⁻⁺) — no mixing of primed and unprimed within a form.

The MS solver consumes `truncated`. ExactSFIPhase (full-MOM correction in arch doc §6)
reads `exact` to build J_exact and `truncated` to build J_trunc. The standalone solver's
FO-equivalent paths use `truncated` for paths 1+2 in TMS reconstruction mode and `exact`
for paths 3+4 (which don't have a clean TMS interpretation).
"""
struct TruncatedAndExactScatteringOpticalProperties{T, E, FT}
    truncated::T   # CoreScatteringOpticalProperties — primed (τ', ϖ', Z'); for MS solver / J_trunc
    exact::E       # CoreScatteringOpticalProperties — unprimed (τ, ϖ, Z); for paths 3+4 / J_exact
    truncfac::Vector{FT}    # f per layer; metadata
end
```

**Critical**: each form is a full triplet. The truncated triplet has `τ' = (1-fϖ)τ`, `ϖ' = (1-f)ϖ/(1-fϖ)`, `Z'` = δ-M-scaled. The exact triplet has the unprimed `(τ, ϖ, Z)`. Mixing forms (e.g., primed τ with unprimed ϖ) breaks the energy budget — see v0.6 §6.3 for derivation.

### 3.3 The standalone solver's input context

```julia
struct ExactSSConfig{FT, OPT, GEO, SUR, CTRB}
    optical_properties::OPT       # TruncatedAndExactScatteringOpticalProperties
    contributors::CTRB            # Tuple for `exact_phase_function` evaluation per layer
    surface::SUR
    geometry::GEO
    n_stokes::Int                 # 1, 3, or 4 — pol-type from day one
    I0::Vector{FT}
    paths::Symbol                 # :auto, :all, :paths_1_2, :all_four
    quadrature::Symbol            # :auto, :required, or explicit Int
    quadrature_tolerance::FT
    jacobian_method::Symbol       # :forwarddiff, :handcoded_f2_ad_f1, :both
end
```

The kernels dispatch on `n_stokes` via the type of `Z⁺⁺/Z⁻⁺` matrices. A scalar implementation uses 1×1 matrices; a vector implementation uses 4×4. The kernel code is the same; the matrix size dispatches.

---

## 4. Kernel-based dispatch

[Section 4 unchanged from v0.3 — KernelAbstractions, `(iv, s)` natural batching unit, CPU/GPU agnostic.]

Per-path kernels launched as data-parallel kernels over geometry × spectral grid. Inner loop over layers (and over inner quadrature points for paths 3+4) is sequential per `(iv, s)` thread.

---

## 5. Required-Nquad framework

[Section 5 unchanged from v0.3 — empirically-derived Nquad from contributors and surface, both for standalone inner quadrature and for full MOM stream count.]

### 5.1 Why this matters

VLIDORT enforces `NMOMENTS = min(2·NSTREAMS - 1, NGREEK_MOMENTS_INPUT)` (`vlidort_main/regular/vlidort_inputs.f90:4138-4141`) but does not auto-derive NSTREAMS from the aerosols or surface. We can — modern Julia ecosystem makes empirical convergence-driven Nquad selection tractable.

### 5.2 Aerosol Greek-coefficient series convergence

```julia
function determine_required_l_aerosol(aerosol, target_relative_error=1e-4)
    β = aerosol.greek_coefficients     # β_l
    L_max = length(β) - 1
    tail_cumulative = reverse(cumsum(reverse(β.^2)))
    threshold = tail_cumulative[1] * target_relative_error^2
    L = findfirst(i -> tail_cumulative[i+1] < threshold, 1:L_max-1)
    return isnothing(L) ? L_max : L
end
```

### 5.3 Surface BRDF Fourier convergence

```julia
function determine_required_nbrdf_coxmunk(wind_speed_ms, target_relative_error=1e-4)
    σ² = 0.003 + 0.00512 * wind_speed_ms     # Cox & Munk 1954
    m_max = ceil(Int, log(1/target_relative_error) / sqrt(σ²))
    return clamp(m_max, 16, 200)
end
```

### 5.4 Combined for the MOM solver

```julia
function determine_required_nstreams(model, iBand; target_relative_error=1e-4)
    contributors = model.atmosphere.contributors[iBand]
    surface = model.surfaces[iBand]
    L_aerosol = maximum(determine_required_l_aerosol(c, target_relative_error)
                        for c in contributors if c isa AbstractScatteringContributor;
                        init=2)
    Nstreams_aero = ceil(Int, (L_aerosol + 1) / 2)
    Nbrdf = determine_required_nbrdf(surface, target_relative_error)
    return (; Nstreams=Nstreams_aero, NSTREAMS_BRDF=Nbrdf,
              L_aerosol, surface_kind=typeof(surface).name.name)
end
```

### 5.5 Inner-quadrature Nquad for paths 3+4

```julia
function determine_required_nquad_inner(τ_total, contributors, target_relative_error=1e-4)
    g_residual = maximum(residual_asymmetry(c) for c in contributors;
                         init=zero(typeof(target_relative_error)))
    τ_factor = log(1 + τ_total)
    N = ceil(Int, 8 + 4 * g_residual + 6 * τ_factor)
    return clamp(N, 8, 64)
end
```

### 5.6 Feeds the v0.6 precheck

The framework directly answers v0.6 Q9: `forward_peak_threshold` and `cox_munk_wind_threshold` are derived from this rather than hand-set. Diagnostic record (last fields of the NamedTuple) tells users *why* the precheck flagged each correction.

---

## 6. Linearization as design experiment for hybrid AD-handcoded

### 6.1 The design question

The chain of computations:

```
parameters  ──f₁──→  core optical properties (τ, ϖ, Z⁺⁺, Z⁻⁺, surface params)  ──f₂──→  radiance
```

Three candidate architectures:
- **(A) Pure AD.** ForwardDiff seeds parameters; chain rule propagates through both f₁ and f₂.
- **(B) Pure handcoded.** Manual derivatives everywhere.
- **(C) Hybrid.** Handcoded f₂ (closed-form; speed-critical). AD f₁ (messy upstream physics).

The standalone SS solver tests architecture C as a small-scale prototype for the broader vSmartMOM lin redesign in `docs/dev_notes/LINEARIZATION_BUGS.md`. The seam is `CoreScatteringOpticalProperties` (the existing type — see §3.1).

### 6.2 Path 1 derivative w.r.t. layer τ — derivation done correctly

The closed form for path 1, layer iz contribution:

$$
\Delta I_{\text{layer}} = \frac{\varpi P}{4\pi \mu_v a} I_0 \left[\exp(-\tau_{\text{top}} a) - \exp(-\tau_{\text{bot}} a)\right]
$$

where $a = 1/\mu_0 + 1/\mu_v > 0$ and $\tau_{\text{top}} < \tau_{\text{bot}}$.

**Same-layer derivative.** Differentiating w.r.t. $\Delta\tau_{iz} = \tau_{\text{bot}} - \tau_{\text{top}}$ holding $\tau_{\text{top}}$ fixed (since $\tau_{\text{top}}$ for layer iz depends only on layers above it):

$$
\frac{\partial}{\partial \Delta\tau_{iz}} \left[\exp(-\tau_{\text{top}} a) - \exp(-\tau_{\text{bot}} a)\right]
= 0 - \exp(-\tau_{\text{bot}} a) \cdot (-a)
= +a \exp(-\tau_{\text{bot}} a)
$$

After the prefactor:

$$
\frac{\partial \Delta I_{\text{layer}}}{\partial \Delta\tau_{iz}} = \frac{\varpi P I_0}{4\pi \mu_v} \exp(-\tau_{\text{bot}} a)
$$

**Positive.** The `1/a` from the prefactor cancels with the `+a` from differentiation.

**Lower-layer derivative.** For a layer iz' *below* iz (so the layer-above-iz' has its $\tau_{\text{top}}$ and $\tau_{\text{bot}}$ shifted by $\Delta\tau_{iz}$ when $\Delta\tau_{iz}$ changes):

$$
\frac{\partial}{\partial \Delta\tau_{iz}} \left[\exp(-\tau_{\text{top}}^{(iz')} a)\right] = -a \exp(-\tau_{\text{top}}^{(iz')} a)
$$

After the prefactor for layer iz':

$$
\frac{\partial \Delta I_{\text{layer iz'}}}{\partial \Delta\tau_{iz}} = -\frac{\varpi^{(iz')} P^{(iz')} I_0}{4\pi \mu_v} \left[\exp(-\tau_{\text{top}}^{(iz')} a) - \exp(-\tau_{\text{bot}}^{(iz')} a)\right]
$$

**Negative** (since the integrand on the right is positive). This is the cumulative-attenuation contribution: increasing $\Delta\tau_{iz}$ shifts all layers below it deeper in optical depth, attenuating their contributions more.

The total path-1 derivative w.r.t. $\Delta\tau_{iz}$ is the sum: same-layer (positive) + sum over lower layers (negative). They typically partially cancel; the sign of the total depends on geometry and atmosphere.

### 6.3 Implementation pattern

```julia
function path1_derivative_dtau!(∂L_∂τ::Array{FT,4}, config)
    nGeom, nStokes, nSpec = output_dims(config)
    nLayers = config.n_layers
    
    for iv in 1:nGeom
        μv = cos(config.geometry.vza[iv])
        μ0 = config.geometry.μ₀
        a = 1/μ0 + 1/μv
        cosΘ = scattering_angle_cosine(μ0, μv, config.geometry.Δϕ[iv])
        
        for s in 1:nSpec, iz in 1:nLayers
            ω = config.optical_properties.exact.ϖ[iz, s]
            P = exact_phase_function(config.contributors, iz, cosΘ, config.n_stokes)
            τ_top = cumulative_τ(config, iz, s)        # primed (truncated)
            τ_bot = cumulative_τ(config, iz+1, s)
            f = config.optical_properties.truncfac[iz, s]
            tms = 1 / (1 - f * ω)
            prefactor = ω * P[1,1] * tms / (4π * μv) * config.I0[s]
            
            # Same-layer (positive)
            ∂L_∂τ[iv, 1, s, iz] += prefactor * exp(-τ_bot * a)
            
            # Lower-layer attenuation (negative): for each iz' > iz, this layer's
            # contribution attenuates more when our layer's τ increases
            for iz_below in (iz+1):nLayers
                ω_b = config.optical_properties.exact.ϖ[iz_below, s]
                P_b = exact_phase_function(config.contributors, iz_below, cosΘ, config.n_stokes)
                τ_top_b = cumulative_τ(config, iz_below, s)
                τ_bot_b = cumulative_τ(config, iz_below+1, s)
                f_b = config.optical_properties.truncfac[iz_below, s]
                tms_b = 1 / (1 - f_b * ω_b)
                pref_b = ω_b * P_b[1,1] * tms_b / (4π * μv) * config.I0[s]
                ∂L_∂τ[iv, 1, s, iz] -= pref_b * (exp(-τ_top_b * a) - exp(-τ_bot_b * a))
            end
        end
    end
end
```

The output `∂L_∂τ` has shape `(nGeom, nStokes, nSpec, nLayers)`. This is the f₂ Jacobian element for τ — one element of the dictionary returned by `exact_ss_path1_with_jacobians`.

Sister functions for ω (simpler, no dependency between layers) and Z (per-moment dispatch) similarly populate their respective Jacobian arrays.

### 6.4 The chain-rule combiner — contraction-correct

The f₂ Jacobian has shape `(nGeom, nStokes, nSpec, nOpticalAxis)` where the optical axis is layer-indexed for τ/ϖ, layer-and-moment-indexed for Z, and surface-parameter-indexed for surface. The f₁ Jacobian has shape `(nOpticalAxis, nParameter)` (or per-spec, depending on whether the parameter is wavelength-resolved). The combined Jacobian must contract over `nOpticalAxis` to produce `(nGeom, nStokes, nSpec, nParameter)`.

**Element-wise broadcasting is wrong** (GPT Finding 4 — confirmed by the existing lin code's "Bug 19" warning at `rt_kernel_lin.jl:32`: matrix products mix indices, so `.*` chain rule becomes incorrect after operator mixing). For SS paths 1+2 (closed-form, no operator mixing) the issue is less severe but the implementation must still contract correctly:

```julia
"""
    chain_rule_combine_dτ(∂L_∂τ, ∂τ_∂p) -> ∂L_∂p

Contracts ∂L_∂τ[iv, is, isp, iLayer] with ∂τ_∂p[iLayer, isp, iParam] over iLayer
(and broadcasts over isp where appropriate) to produce ∂L_∂p[iv, is, isp, iParam].

Both inputs assume per-spec optical-property dependence; if a parameter affects
all layers uniformly across spec, that's encoded in ∂τ_∂p[:, :, iParam] having
the appropriate shape.
"""
function chain_rule_combine_dτ(∂L_∂τ::AbstractArray{FT,4},
                                ∂τ_∂p::AbstractArray{FT,3}) where FT
    nGeom, nStokes, nSpec, nLayers = size(∂L_∂τ)
    nLayers2, nSpec2, nParam = size(∂τ_∂p)
    @assert nLayers == nLayers2 && nSpec == nSpec2
    
    ∂L_∂p = zeros(FT, nGeom, nStokes, nSpec, nParam)
    for ig in 1:nGeom, is in 1:nStokes, isp in 1:nSpec, ip in 1:nParam
        s = zero(FT)
        for il in 1:nLayers
            s += ∂L_∂τ[ig, is, isp, il] * ∂τ_∂p[il, isp, ip]
        end
        ∂L_∂p[ig, is, isp, ip] = s
    end
    return ∂L_∂p
end
```

(Equivalent tensor-product expressions via `tensor` macros or `TensorOperations.jl` are fine; the explicit loop is shown for clarity. Performance-tuned versions can use `mul!` or batched matmul where shapes allow.)

For Z derivatives, the Jacobian shape depends on the Z storage convention — Greek-coefficient per moment, F-matrix per (geometry, layer, spec), or other. Each has its own contraction structure. The combiner has one method per optical-property kind (`τ, ϖ, Z, surface_param`); they assemble their Jacobian contributions, sum into the combined output Jacobian.

**Test for contraction correctness explicitly**:

```julia
@testset "Chain rule contraction is dimensionally correct" begin
    # Setup: ∂L/∂τ from f₂; ∂τ/∂AOD from f₁
    config = standard_test_config()
    layout = ParameterLayout([AODParameter()])  # 1 parameter
    
    L, ∂L_∂τ = path1_with_jacobian_dτ(config)            # (nGeom, nStokes, nSpec, nLayers)
    ∂τ_∂AOD = aod_to_layer_τ_jacobian(config)            # (nLayers, nSpec, 1)
    
    # Combined
    ∂L_∂AOD_chain = chain_rule_combine_dτ(∂L_∂τ, ∂τ_∂AOD)  # (nGeom, nStokes, nSpec, 1)
    
    # Compare to pure ForwardDiff over the whole pipeline
    f(p) = run_path1(set_aod(config, p[1]))[:total]
    ∂L_∂AOD_AD = ForwardDiff.jacobian(f, [config.aod])
    
    @test isapprox(reshape(∂L_∂AOD_chain, :), reshape(∂L_∂AOD_AD, :), rtol=1e-12)
end
```

If the contraction is wrong (e.g., element-wise broadcasting), this test fails immediately — different Jacobian shape or different values.

### 6.5 The benchmark experiment

Three implementations of the same final answer (∂L/∂physics_params for representative cases):

- **(A) Pure ForwardDiff** — seed parameters, propagate through f₁ ∘ f₂.
- **(B) Pure handcoded** — handcoded f₁ + handcoded f₂ + chain rule.
- **(C) Hybrid** — handcoded f₂ + AD f₁ + chain rule.

Run on representative cases: small parameters (1–5), medium (10–20), large (50+).

Measure: forward speed; lin speed; lin speed normalized by parameter count; memory; correctness (all three agree to ~1e-12).

The findings inform the broader `LINEARIZATION_BUGS.md` redesign:

- **f₁/f₂ seam structure** — does `CoreScatteringOpticalProperties` work as the architectural boundary?
- **Chain-rule combiner** — pure Julia, reusable in broader system.
- **Handcoded-f₂ effort and payoff** — derivation time vs speedup; informs whether to extend to the f₂ of the full RT (doubling-adding, etc.).
- **AD-f₁ overhead** — ForwardDiff cost for upstream chain; informs broader AD optimization needs.

### 6.6 Implementation phases for §6

- **Phase 5a**: handcoded f₂ for paths 1+2 (closed forms; ~3-4 days). Path 1 derivation per §6.2.
- **Phase 5b**: chain-rule combiner with contraction-correct implementation. Test for dimensional correctness (§6.4 test). ~1 day.
- **Phase 5c**: methods A, B, C; benchmark experiment; report. ~1 week.
- **Phase 5d (conditional)**: extend handcoded f₂ to paths 3+4 if speedup justifies.

---

## 7. Back-correction of current vSmartMOM — scope-bounded

### 7.1 What this is and what it is NOT

**Is**: post-hoc correction adding `(exact_p1+p2 - truncated_p1+p2)` to `R_SFI` *after* `rt_run` completes. Validates against **VLIDORT FO outputs** (which have the same "first-scatter to sensor only" scope).

**Is NOT**: a substitute for `ExactSFIPhase` (architecture doc §6). The post-hoc back-correction adds at the sensor only — it does **not** propagate the corrected first-scatter source through `(I - M_trunc)^-1`, so first-scatter-then-anything paths (sun→atm→atm→sensor, sun→atm→surface→atm→sensor, etc.) are **not** corrected. These higher-order paths are corrected by `ExactSFIPhase`'s internal source modification, which lives inside `rt_run`'s solver loop.

Per v0.6 architecture doc §6.2: the production full-MOM correction is `ΔI = (I - M_trunc)^-1 (J_exact - J_trunc)`, structurally distinct from the post-hoc subtraction. Mixing the two architectures is what GPT Finding 1 flagged as physically wrong.

### 7.2 Why the back-correction is still useful

Three scoped uses, all valid:

1. **Validation against VLIDORT FO**: VLIDORT FO has the exact same scope as the post-hoc back-correction (first-scatter to sensor). The Task1↔Task3 difference in `solar_tester` *is* VLIDORT's FO contribution. If our back-correction matches that difference, the standalone solver's path 1+2 implementation is correct.

2. **Validation of the standalone solver against real vSmartMOM optical setups**: feeds real model contributors and surfaces through the standalone solver, compares against reference computations. Catches bugs that synthetic fixtures might miss.

3. **Fast first-pass tool**: even before `ExactSFIPhase` is implemented in `rt_run`, users can call `apply_back_correction!` on existing `rt_run` outputs to get a *partial* correction. Documented as partial (FO-equivalent, not full-MOM-correction) so users don't mistake it for production.

### 7.3 The truncated reconstruction

To compute `(exact_p1+p2 - truncated_p1+p2)`, we need both forms.

`exact_p1+p2` comes from `run_exact_ss(config; paths=:paths_1_2).total` — the standalone solver in its primary mode.

`truncated_p1+p2` is reconstructed: same closed-form structure as `exact_p1+p2`, but using the truncated phase representation that vSmartMOM's MS solver internally produces. Specifically, for path 1, evaluate the truncated phase at the actual scattering angle from primed Greek coefficients up to `l_trunc`. For path 2, reconstruct the BRDF Fourier sum at the actual viewing geometry from m=0 to m=max_m-1 — Sanghavi's existing `_fourier_coeff_element` is the prototype (`coxmunk_surface.jl:505-546`).

```julia
"""
    truncated_ss_path1(config, l_trunc) -> Array

Compute path 1 SS with phase function truncated at l_trunc moments.
This is what vSmartMOM's MS solver implicitly produces for the SS-into-view
contribution at viewing geometry — when ExactSFIPhase is OFF.
"""
function truncated_ss_path1(config, l_trunc) end

"""
    truncated_ss_path2(config, max_m) -> Array

Path 2 with BRDF Fourier-reconstructed at viewing geometry from m = 0..max_m-1.
"""
function truncated_ss_path2(config, max_m) end
```

### 7.4 Validation that the truncated reconstruction matches vSmartMOM

```julia
# Sanity: in a case where vSmartMOM agrees with VLIDORT-no-FO (both compute truncated solve),
# our truncated reconstruction should match the SS portion of vSmartMOM's R_SFI.
let
    # Run vSmartMOM with sfi_phase=TruncatedSFIPhase (current behavior)
    result_trunc = run_rt(model_with_truncated_sfi)
    
    # Run VLIDORT with DO_FOCORR=False
    result_vlidort_no_fo = vlidort_task1_outputs   # From shipped solar_tester results
    
    @test result_trunc.R_SFI ≈ result_vlidort_no_fo  rtol=1e-5
    
    # Now reconstruct
    truncated_p12 = truncated_ss_path1(config, l_trunc) + truncated_ss_path2(config, max_m)
    
    # The truncated reconstruction should match what's "in" R_SFI as the SS contribution
    # (this requires extracting the SS portion specifically; see §7.5)
end
```

This is the test for whether we understand vSmartMOM's discretization. If we pass, the back-correction is well-defined.

### 7.5 The adapter

```julia
"""
    apply_back_correction!(R_SFI, model, vsmartmom_metadata)

FO-EQUIVALENT DIAGNOSTIC. Adds (exact - truncated) for paths 1+2 to R_SFI.
This corrects the directly-escaped first-scatter only (FO scope); does NOT
correct first-scatter-then-anything paths. For full-MOM correction including
higher-order paths, use FullMOM(sfi_phase=ExactSFIPhase()) which lives inside
rt_run rather than this post-hoc adapter.
"""
function apply_back_correction!(R_SFI, model, vsmartmom_metadata)
    config = config_from_vsmartmom_model(model, vsmartmom_metadata.iBand)
    exact_p12 = run_exact_ss(config; paths=:paths_1_2).total
    truncated_p12 = (truncated_ss_path1(config, vsmartmom_metadata.l_trunc) +
                     truncated_ss_path2(config, vsmartmom_metadata.max_m))
    R_SFI .+= exact_p12 .- truncated_p12
end
```

Used as a diagnostic, not a production substitute. Docstring is explicit about scope.

### 7.6 What we learn from the back-correction validation

- **Back-correction matches VLIDORT FO Task1↔Task3 difference**: standalone solver paths 1+2 are correct. Refactor proceeds with confidence at the FO layer.

- **Doesn't match**: standalone solver or truncated reconstruction has a bug. Investigate which before proceeding with `ExactSFIPhase` implementation.

- **Matches FO but vSmartMOM-with-back-correction still doesn't match VLIDORT-with-FO-on**: differences are higher-order paths that FO doesn't fix and our back-correction doesn't either. This is *expected*; the residual should match what VLIDORT-with-FO vs VLIDORT-no-FO-no-FO doesn't capture either, and quantifies the value of full `ExactSFIPhase` over plain FO.

In all three outcomes we have actionable evidence.

---

## 8. Validation against existing references

### 8.1 Against `dev_notes/exact_ss_reference/`

The reference (delivered with this document set) is the correctness oracle. Each kernel passes the corresponding reference fixture to ~1e-8 relative.

### 8.2 Against VLIDORT shipped fixtures (Piece B)

Detailed in `vlidort_baseline_suite.md`. Standalone solver in `paths=:paths_1_2` mode matches:
- Siewert (2000) Problem IIA (peer-reviewed; in `vlidort_v_test/saved_results/.../results_Siewert2000_validation.all`).
- Solar tester Task1↔Task3 difference (in `vlidort_s_test/saved_results/.../results_solar_tester.all`).

### 8.3 Internal consistency

- BRDF reciprocity for total: `total(μ₀, μ_v)/μ₀ = total(μ_v, μ₀)/μ_v` for Lambertian + plane-parallel + ϖ=1.
- Thin-limit: SS contributions linear in τ at small τ.
- Energy conservation: hemispherical reflectance bounded.
- Path-by-path: black surface → only path 1; conservative absorption → only path 2.

---

## 9. Development phases

### Phase 1 — Forward solver, Lambertian, paths 1+2 (pol-type generic from day one)
Module skeleton; `ExactSSConfig` aligning with `CoreScatteringOpticalProperties`; kernels for path 1 and path 2 with n_stokes-parametric Z matrix multiplication; orchestration. Tests vs reference. **Phase 1 specializes for n_stokes=1** for first concrete implementation, but the data structures and kernel signatures are pol-type generic.
**Deliverable**: validated Lambertian paths 1+2; framework ready for vector extension.

### Phase 2 — Required-Nquad framework
`determine_required_l_aerosol`, `determine_required_nbrdf_*`, `determine_required_nstreams`. Empirical convergence calibration.
**Deliverable**: framework calibrated against test aerosols and surfaces.

### Phase 3 — Paths 3+4 (Lambertian) + inner-quadrature auto-Nquad
Kernels using framework from Phase 2.
**Deliverable**: all four paths working for Lambertian.

### Phase 4 — Cox-Munk path 2 + ForwardDiff compatibility audit
Cox-Munk path 2 kernel; ForwardDiff audit; per-path Jacobian via FD comparison.
**Deliverable**: Cox-Munk path 2; pure-AD method (Architecture A) works end-to-end.

### Phase 5 — Hybrid linearization design experiment (§6)
- 5a: handcoded f₂ for paths 1+2 (sign-correct per §6.2)
- 5b: chain-rule combiner with contraction-correctness test (§6.4)
- 5c: methods A, B, C benchmark; report
- 5d (conditional): extend handcoded f₂ to paths 3+4 if 5c justifies

**Deliverable**: design experiment complete; report informs `LINEARIZATION_BUGS.md` redesign.

### Phase 6 — FO-equivalent back-correction (§7)
`truncated_ss_path1`, `truncated_ss_path2`, `apply_back_correction!` with scope-bounded documentation. Validation:
- Truncated reconstruction matches vSmartMOM-no-FO output (§7.4)
- Back-corrected vSmartMOM matches VLIDORT FO Task1↔Task3 (§7.6 outcome 1)

**Deliverable**: validated back-correction with bounded scope; ready for use as diagnostic alongside the (separate, in-progress) `ExactSFIPhase` implementation.

### Phase 7 — Vector polarization (n_stokes=3) — when integration is ready
The framework is generic from Phase 1. Concrete vector implementations of paths 1, 2 specialize the kernel matrix arithmetic for 4×4 Z matrices (or 3×3 for n_stokes=3 omitting V). Validates against Siewert (2000) Tables 3, 4 (Q, U). Coordinated with the broader vSmartMOM polarization story.

Each phase small, testable. Estimated total: 8–10 weeks focused work.

---

## 10. What this enables

After all phases:

- Fast standalone exact-SS solver. Pol-type generic.
- Validated against `dev_notes/exact_ss_reference/` AND against VLIDORT shipped fixtures (Siewert 2000 + solar_tester) AND via back-correction of real vSmartMOM cases.
- ForwardDiff-friendly. Hand-coded f₂ derivatives for paths 1+2 (sign-correct, contraction-correct chain rule). Hybrid AD-handcoded architecture proven at small scale.
- Required-Nquad framework empirically calibrated; informs v0.6 precheck thresholds.
- FO-equivalent back-correction available as diagnostic; full-MOM `ExactSFIPhase` implementation proceeds independently.
- Empirical answers for `LINEARIZATION_BUGS.md` broader redesign: should it use hybrid AD-handcoded? How should the f₁/f₂ seam (`CoreScatteringOpticalProperties`) be structured?

When the bigger refactor (`vsmartmom_dispatch_design_v0_6.md`) integrates the SS subsystem, the kernels are already built, validated, and proven against the existing pol-type-generic abstraction.

---

## 11. Open questions

**SQ1**: Where exactly should the standalone solver live? Current proposal: `src/StandaloneSS/` in vSmartMOM.

**SQ2**: Z storage in path-1 implementation. Current proposal: `exact_phase_function` trait dispatching on contributor; returns SMatrix evaluated at the actual scattering angle. Compatible with both Greek-coefficient and direct-evaluation contributor representations.

**SQ3**: How many physics-parameter cases for the §6 benchmark experiment? Sweep small (5), medium (20), large (100+); each scales differently.

**SQ4**: Auto-quadrature heuristic constants — calibration set. Pick 10 representative aerosols + surfaces; tune until empirical convergence is matched.

**SQ5**: Multi-band handling. Per-band orchestration; batched-across-bands deferred.

**SQ6**: f₁ source for §6 benchmark experiment. Current proposal: parametrized scaling of precomputed properties (AOD scales τ, HG g param, etc.) — avoids Mie complexity, isolates the linearization architecture question.

**SQ7**: Phase 7 (vector polarization) timing — coordinate with broader vSmartMOM polarization work? Or implement standalone first then propagate?

**SQ8**: Should the `truncated_ss_path1` reconstruction match vSmartMOM's discretization to numerical precision, or just to "good enough for diagnostic"? Numerical precision is harder but cleaner; "good enough" is faster to implement.

---

## 12. Cross-references

**Delivered together with this document:**
- `vsmartmom_dispatch_design_v0_6.md` — broader architecture; especially §5.5 (dual-form optical properties) and §6 (`AbstractSFIPhase`).
- `vlidort_baseline_suite.md` — golden-standard validation; uses shipped Siewert 2000 + solar_tester fixtures.
- `dev_notes/exact_ss_reference/` — correctness oracle.

**Pre-existing committed in `vSmartMOM.jl`:**
- `docs/dev_notes/LINEARIZATION_BUGS.md` — broader redesign track that §6 informs.
- `src/CoreRT/types.jl` line 1052 — `CoreScatteringOpticalProperties{FT,FT2,FT3}` with `τ::FT, ϖ::FT2, Z⁺⁺::FT3, Z⁻⁺::FT3` — the pol-type-generic core abstraction.
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl` line 32 — "Bug 19" precedent for chain rule applied before index-mixing.
- `coxmunk_surface.jl:505-546` — Sanghavi's existing `apply_ss_correction!`.
- `vlidort_main/regular/vlidort_inputs.f90:4138-4141` — `NMOMENTS = min(2·NSTREAMS-1, NGREEK_MOMENTS_INPUT)` formula.
- `vlidort_focode/FO_VectorSS_RTCalcs_I.f90` — VLIDORT FO source.
