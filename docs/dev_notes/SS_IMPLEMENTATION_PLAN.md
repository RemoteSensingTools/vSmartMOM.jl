# SS implementation plan — kernel-based, non-conflicting development

**Document version:** v0.2
**Status:** working plan; reviewable independently of the architecture doc.

**Document set (delivered together):**
- `vsmartmom_dispatch_design_v0_6.md` — architecture doc.
- `standalone_ss_solver_plan.md` — Piece A: the standalone kernel solver this plan eventually integrates.
- `vlidort_baseline_suite.md` — Piece B: validation suite.
- `SS_IMPLEMENTATION_PLAN.md` — *this document*; integration plan from standalone to production.

**Pre-existing committed in `vSmartMOM.jl`:**
- `docs/dev_notes/LINEARIZATION_BUGS.md` — linearization redesign.
- `src/CoreRT/types.jl` line 1052 — `CoreScatteringOpticalProperties{FT,FT2,FT3}`.
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl` line 32 — Bug 19 chain-rule precedent.

### Changes from v0.1
- **§2 (module structure)**: `src/CoreRT/ExactSS/` builds *on top of* the standalone solver from Piece A, rather than reimplementing it. The standalone solver is the kernel library; the ExactSS module is the integration layer.
- **§4 (truncated reconstruction)**: aligned with corrected scope — back-correction is for FO-equivalent diagnostic, `ExactSFIPhase` (architecture doc §6) is for full-MOM correction. Two integration paths, separately tested.
- **§7 (development phases)**: integration phase (Phase 5) is the swap-in PR after Piece A and Piece B validation are complete.

---

## 0. Goals and constraints

This document specifies *how to build* the integrated SS subsystem in vSmartMOM, after the standalone solver (Piece A) and validation suite (Piece B) are complete. Two goals:

1. **Develop in parallel with current code** — no merge conflicts. SS code lives in new files; no existing files modified during development.

2. **Build on the validated standalone solver** — kernel implementations from Piece A are reused; the integration layer wires them into vSmartMOM's dispatch.

Constraints:
- Match `dev_notes/exact_ss_reference/` to ~10⁻⁸.
- Match Piece B Stage 1 fixtures (Siewert 2000, solar_tester) within established tolerances.
- ForwardDiff-friendly throughout.
- CPU/GPU agnostic via the standalone solver's KernelAbstractions infrastructure.
- Final integration with `rt_run.jl` is a single small PR.

---

## 1. Relationship to other documents

This plan covers the **integration layer** between Piece A's standalone solver and vSmartMOM's `rt_run`. It does *not* cover:

- The kernel implementations themselves (those are in Piece A).
- The validation infrastructure (that's in Piece B).
- The architectural decisions (those are in the v0.6 dispatch design).

What it covers:
- The `src/CoreRT/ExactSS/` module that depends on Piece A's standalone solver and provides the dispatch surface for `apply_correction!` and `ExactSFIPhase` source modification.
- The integration test suite ensuring the integration layer doesn't break existing physics.
- The phased integration plan with Piece B validation gates between phases.

---

## 2. Module structure

```
src/CoreRT/ExactSS/                            # Integration layer (new)
├── ExactSS.jl                                 # Module entry; imports StandaloneSS
├── apply_correction.jl                        # apply_correction! dispatch on AbstractSSPath
├── exact_sfi_phase.jl                         # ExactSFIPhase source modification (architecture §6)
├── precheck_integration.jl                    # Wires precheck_truncation into rt_run flow
└── back_correction_adapter.jl                 # FO-equivalent diagnostic (Piece A §7) wired up

src/StandaloneSS/                              # Kernel library (Piece A)
├── ... (per Piece A §2)

test/CoreRT/ExactSS/                           # Integration tests
├── runtests.jl
├── test_apply_correction_dispatch.jl
├── test_exact_sfi_phase_vs_truncated.jl       # Regression: ExactSFIPhase off = current behavior
├── test_full_pipeline_vs_vlidort.jl           # Uses Piece B fixtures
└── test_back_correction_diagnostic.jl
```

**Critical**: `src/CoreRT/ExactSS/` *imports* from `src/StandaloneSS/`. The standalone solver is the kernel library; the ExactSS module is the integration glue. This keeps the kernels self-contained (Piece A independently testable) while the integration layer is the only piece that touches `rt_run.jl`.

---

## 3. The integration surface

Three integration points, all adding new code without modifying existing files until Phase 5:

### 3.1 `apply_correction!` dispatch (for `ExactFirstOrderOnly` and FO-equivalent diagnostic)

In `src/CoreRT/ExactSS/apply_correction.jl`:

```julia
"""
    apply_correction!(R_SFI, path::AbstractSSPath, solver, model, ctx)

Dispatched on path type and solver type. For ExactFirstOrderOnly: adds the
exact path contribution. For FullMOM: adds (exact - truncated) — FO-equivalent
diagnostic; does NOT propagate through (I - M_trunc)^-1; for that use ExactSFIPhase.
"""
function apply_correction!(R_SFI, p::AtmosphericSSPath, ::ExactFirstOrderOnly, model, ctx)
    # Adds the exact path-1 contribution
    L = StandaloneSS.exact_ss_path1(make_ss_config(model, ctx))
    R_SFI .+= L
end

function apply_correction!(R_SFI, p::AtmosphericSSPath, ::FullMOM, model, ctx)
    # Diagnostic FO-equivalent (exact - truncated); Piece A §7
    config = make_ss_config(model, ctx)
    exact = StandaloneSS.exact_ss_path1(config)
    truncated = StandaloneSS.truncated_ss_path1(config, model.solver.l_trunc)
    R_SFI .+= exact .- truncated
end

# Similar for SurfaceDirectBeamPath, AtmosphereToSurfacePath, etc.
```

### 3.2 `ExactSFIPhase` source modification (in `elemental.jl`)

This is the *production* full-MOM correction (architecture doc §6). It's not a separate adapter; it's a modification to how the `J0±` source vectors are constructed in `elemental.jl`.

Currently `elemental.jl` uses truncated phase moments to build `J0±`. The change:

```julia
# Current behavior (becomes TruncatedSFIPhase)
function build_J_source!(J0_plus, J0_minus, ::TruncatedSFIPhase,
                         optical_props::TruncatedAndExactScatteringOpticalProperties, ...)
    # Use optical_props.truncated (current behavior preserved)
    _build_J_source_internal!(J0_plus, J0_minus, optical_props.truncated, ...)
end

# New behavior
function build_J_source!(J0_plus, J0_minus, ::ExactSFIPhase,
                         optical_props::TruncatedAndExactScatteringOpticalProperties, ...)
    # Use optical_props.exact (untruncated τ, ϖ, Z); the MS operator stays truncated
    _build_J_source_internal!(J0_plus, J0_minus, optical_props.exact, ...)
end
```

`_build_J_source_internal!` is the existing source-construction machinery; it operates on a `CoreScatteringOpticalProperties` triplet. The dispatch on `AbstractSFIPhase` selects which triplet (truncated or exact) to feed it.

The MS solver `(I - M_trunc)` continues to use `optical_props.truncated`. Only the source `J_*` changes.

### 3.3 `precheck_truncation` wired into `rt_run`

Calls `precheck_truncation` (architecture doc §3.6) at context construction; resolves `:auto` to a concrete `ss_paths` tuple plus diagnostic record. Stored in `ctx.precheck_results` for transparency.

---

## 4. Truncated reconstruction — for the back-correction diagnostic

The back-correction adapter (Piece A §7) needs a *truncated* SS reconstruction at viewing geometry. This is the standalone-solver's `truncated_ss_path1(config, l_trunc)` and `truncated_ss_path2(config, max_m)` — both implemented in Piece A, called from the integration layer's `apply_correction!(::AbstractSSPath, ::FullMOM, ...)`.

**Important scoping**: this back-correction is FO-equivalent only. It validates against VLIDORT FO outputs (Piece B Case D). It does **not** substitute for `ExactSFIPhase`. The integration tests must verify this distinction:

```julia
@testset "Back-correction is FO-equivalent only" begin
    model = standard_test_model()
    
    # FullMOM with ExactSFIPhase (production correction)
    model_exact_sfi = with_solver(model, FullMOM(sfi_phase=ExactSFIPhase()))
    result_exact_sfi = run_rt(model_exact_sfi)
    
    # FullMOM with TruncatedSFIPhase + back-correction (FO-equivalent diagnostic)
    model_trunc_with_bc = with_solver(model,
        FullMOM(sfi_phase=TruncatedSFIPhase(), ss_paths=:auto))
    result_trunc_with_bc = run_rt(model_trunc_with_bc)
    
    # These should NOT be equal — they correct different sets of paths
    diff = maximum(abs.(result_exact_sfi.R_SFI - result_trunc_with_bc.R_SFI))
    @test diff > 1e-7  # Non-zero — they do different things
    @test diff < 1e-2  # But not wildly different — both are corrections of similar magnitude
    
    # The difference is "first-scatter-then-anything" paths that the back-correction misses
end
```

This test is critical for catching architectural confusion in the future.

---

## 5. The dual-form `TruncatedAndExactScatteringOpticalProperties`

[See architecture doc §5.5 for full type definition.]

This type carries both forms as full triplets. Layer construction populates both:

```julia
function build_layer_optical_properties(layer_data, contributors, l_trunc, ...)
    # Compute truncated form: δ-M-scaled τ, ϖ, Z
    τ′, ϖ′, Z⁺⁺′, Z⁻⁺′, f = compute_truncated_form(layer_data, l_trunc)
    truncated = CoreScatteringOpticalProperties(τ=τ′, ϖ=ϖ′, Z⁺⁺=Z⁺⁺′, Z⁻⁺=Z⁻⁺′)
    
    # Compute exact form: untruncated τ, ϖ, Z
    τ, ϖ, Z⁺⁺, Z⁻⁺ = compute_exact_form(layer_data)
    exact = CoreScatteringOpticalProperties(τ=τ, ϖ=ϖ, Z⁺⁺=Z⁺⁺, Z⁻⁺=Z⁻⁺)
    
    return TruncatedAndExactScatteringOpticalProperties(truncated, exact, f)
end
```

When the user disables `ExactSFIPhase` (running with `TruncatedSFIPhase`), only the truncated form is consumed. When `ExactSFIPhase` is enabled, both forms are consumed: truncated for the MS operator, exact for the source.

The construction cost of computing both forms is small (mostly the δ-M scaling factors); both are cached on the layer.

---

## 6. Integration testing strategy

Integration tests live in `test/CoreRT/ExactSS/`. They build on Piece A's standalone tests (which validate the kernels themselves) and Piece B's VLIDORT comparison (which validates physics agreement).

### 6.1 Regression: TruncatedSFIPhase preserves current behavior

```julia
@testset "TruncatedSFIPhase preserves existing rt_run output" begin
    # Build identical models, one with current rt_run, one with new run_rt + TruncatedSFIPhase
    model_current = build_current_model()
    model_new = build_new_model_with_truncated_sfi()
    
    result_current = rt_run(model_current)
    result_new = run_rt(model_new)
    
    @test result_current.R_SFI ≈ result_new.R_SFI  rtol=1e-12
end
```

This is the *most important* integration test. If `TruncatedSFIPhase` does anything different from current `rt_run`, the integration broke something.

### 6.2 Forward direction: ExactSFIPhase improves agreement with full-physics references

```julia
@testset "ExactSFIPhase reduces error vs Siewert 2000" begin
    siewert_config = load_siewert_2000_fixture()
    siewert_reference_I = parse_results_file("results_Siewert2000_validation.all", :I)
    
    # With TruncatedSFIPhase (current behavior)
    result_trunc = run_rt(with_solver(siewert_config, FullMOM(sfi_phase=TruncatedSFIPhase())))
    error_trunc = relative_error(result_trunc.I, siewert_reference_I)
    
    # With ExactSFIPhase (production correction)
    result_exact = run_rt(with_solver(siewert_config, FullMOM(sfi_phase=ExactSFIPhase())))
    error_exact = relative_error(result_exact.I, siewert_reference_I)
    
    @test error_exact < error_trunc  # ExactSFIPhase should be strictly better
    @test error_exact < 1e-4  # Should reach published-reference quality
end
```

This is the substantive validation that `ExactSFIPhase` is the right architecture.

### 6.3 Lin-direction: regression and ExactSFIPhase Jacobians

```julia
@testset "Jacobians via run_rt match existing rt_run_lin output (TruncatedSFIPhase)" begin
    # Same regression as §6.1 but for Jacobians
end

@testset "Jacobians with ExactSFIPhase agree with finite differences" begin
    # ExactSFIPhase Jacobians validated against perturbation FD
end
```

---

## 7. Development phases

Phases 1-4 happen entirely in `src/CoreRT/ExactSS/` and `src/StandaloneSS/` — no modification of existing files. Phase 5 is the integration PR.

### Phase 1 (after Piece A Phases 1-2 complete) — Basic integration
- `src/CoreRT/ExactSS/` module skeleton
- `apply_correction!(::AtmosphericSSPath, ::ExactFirstOrderOnly, ...)` calling `StandaloneSS.exact_ss_path1`
- `apply_correction!(::SurfaceDirectBeamPath{<:LambertianSurface}, ::ExactFirstOrderOnly, ...)`
- Configuration adapters from RTModel to `ExactSSConfig`

**Deliverable**: integration layer can call standalone kernels for paths 1+2 Lambertian. Tests pass.

### Phase 2 — `TruncatedAndExactScatteringOpticalProperties` type + dual-form construction
- `TruncatedAndExactScatteringOpticalProperties` type
- Layer construction populates both forms
- Backward compatibility: existing single-form constructors still work for non-`ExactSFIPhase` paths

**Deliverable**: dual-form optical properties; existing rt_run unchanged.

### Phase 3 — `AbstractSFIPhase` source modification in `elemental.jl`
- `build_J_source!` dispatch on `AbstractSFIPhase`
- `TruncatedSFIPhase` (current behavior) and `ExactSFIPhase` implementations
- *Note: this phase modifies `elemental.jl`* — it's the smallest change to existing code that's structurally needed for `ExactSFIPhase`.

**Deliverable**: ExactSFIPhase source modification working.

### Phase 4 — Back-correction adapter integration + Piece B validation
- `apply_correction!(::AbstractSSPath, ::FullMOM, ...)` for FO-equivalent diagnostic
- Integration tests against Piece B Case D
- Documentation: scope-bounded back-correction; not substitute for ExactSFIPhase

**Deliverable**: back-correction adapter integrated with explicit scope.

### Phase 5 — `run_rt` consolidation (integration PR)
This is the only phase touching `rt_run.jl`, `rt_run_lin.jl`, `rt_run_multisensor.jl`:
- Replace `apply_ss_correction!` call with `apply_correction!` dispatch loop
- Wire `model.solver` field through `RTModel`
- Add backward-compat wrappers for `rt_run`, `rt_run_lin`, `rt_run_multisensor`
- Run full Piece B Stage 1 + 2 validation suite

**Deliverable**: unified `run_rt(model; mode, outputs)` entry point. All Piece B tests pass. No regressions on existing test suite.

### Phase 6 — Cox-Munk path 2 + canopy hotspot integration
- After Piece A Phase 3 (Cox-Munk path 2) and Phase 6 (back-correction validation) complete
- Cox-Munk path 2 specialization in integration layer
- Canopy hotspot integration via dispatch

**Deliverable**: Full surface coverage in dispatch.

Phases 1-4 can proceed in parallel with Piece A's Phases 3-6. Phase 5 sequencing depends on completion of both Piece A and Piece B Stage 1.

---

## 8. Non-conflicting development

Phases 1-4 modify only files under `src/CoreRT/ExactSS/`, `src/StandaloneSS/`, and `test/`. Phase 3 also modifies `elemental.jl` for `AbstractSFIPhase` source dispatch — this is the *only* existing-file modification before Phase 5. Phase 3 is intentionally minimal: just the source-construction dispatch, nothing else.

Concurrent work on `sanghavi-unified` can proceed without merge conflicts on Phases 1, 2, 4. Phase 3's `elemental.jl` change should be a small, focused PR that's easy to review and unlikely to conflict.

The `RTModel.solver` field addition is part of Phase 5; until then `ExactSS` works with a substituted solver field for testing.

---

## 9. Summary

Six phases from "standalone kernels work" to "production integration":

- **Phase 1** (~3 weeks after Piece A Phase 1-2): basic integration layer
- **Phase 2** (~3 weeks): dual-form `TruncatedAndExactScatteringOpticalProperties`
- **Phase 3** (~3 weeks): `AbstractSFIPhase` source dispatch in `elemental.jl`
- **Phase 4** (~3 weeks): back-correction adapter integrated
- **Phase 5** (~4 weeks): `run_rt` consolidation; full validation
- **Phase 6** (~4 weeks): Cox-Munk + canopy hotspot

Total: ~20 weeks integration work, after Piece A and Piece B Stage 1 are complete.

Phases 1-4 happen in parallel-development territory: only Phase 3 touches existing code minimally. Phase 5 is the integration PR.

The result: a kernel-based, ForwardDiff-friendly SS subsystem that uses Piece A's validated standalone kernels, integrates cleanly with the v0.6 architecture, supports both `FullMOM(sfi_phase=ExactSFIPhase())` (production) and `FullMOM(sfi_phase=TruncatedSFIPhase())` (regression baseline), and provides FO-equivalent back-correction as a scope-bounded diagnostic alongside the production correction.
