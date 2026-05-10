# Plan: Merge sanghavi + unified-vsmartmom → sanghavi-unified

*Created: 2026-03-21. See also: `raman_gpu_optimization.md` for the GPU performance audit.*

## Context

Two branches diverged from a common ancestor (702cbc3d, Apr 2022):
- **unified-vsmartmom** (244 commits): RTModel hierarchy, GPU/Float32 support, linearized RT, 474 tests, docs. Already ported sanghavi forward physics (commit 59f8de8).
- **sanghavi** (20 commits, plus unpushed local changes): Raman physics improvements (RRS, VS, A-band), single-scatter approx, application code (EMIT, Balsamic). Structures/notations have changed and existing sanghavi scripts need updating.

**Goal**: Create `sanghavi-unified` branch combining both, with Raman GPU optimization, eventually becoming the new `main`.

**Critical finding**: The Raman path has severe GPU performance issues — not just the batching problem, but pervasive memory allocation and deallocation inside hot loops that will dominate GPU runtime. See `raman_gpu_optimization.md` for the full audit (~19,000 GPU allocations per Raman run).

---

## GPU Performance Audit Summary

See `docs/dev_notes/raman_gpu_optimization.md` for full details. Key issues:

| Issue | Location | Impact |
|-------|----------|--------|
| 13× `similar()` + `.=0` per call | `interaction_inelastic.jl:293-305` | ~8,300 GPU allocs/run |
| `⊠` temps in nested Δn loops | `doubling_inelastic.jl:62-110` | ~3,600 allocs/run |
| `batch_inv!` allocates pivots | `gpu_batched_cuda.jl:34,63` | ~2,000 allocs/run |
| 2D `*` instead of batched `⊠` | `interaction_inelastic.jl:70-86` | nSpec×nRaman sequential BLAS |
| No `InelasticWorkspace` | (missing) | All above are avoidable |
| Excessive `synchronize_if_gpu()` | `batch_inv!`, doubling | ~2,000 GPU stalls/run |

**Elastic (noRS) path is fine** — uses pre-allocated `RTWorkspace`, `⊠` notation works well with JIT. The `⊠` readability is a priority; don't convert elastic path to in-place ops.

---

## Phase 1: Branch Creation & Sanghavi Compatibility

### 1a. Create branch
- Create `sanghavi-unified` from `unified-vsmartmom`

### 1b. Port remaining sanghavi content
- Single-scatter approximation kernel → adapt to RTModel conventions
- EMIT/Balsamic scripts → `examples/` directory, adapted to RTModel API
- Verify any remaining physics diffs in `elemental_inelastic.jl`, `raman_atmo_prop.jl`

### 1c. Migration Guide for Sanghavi Branch Owner

Write a detailed **migration guide** (`docs/dev_notes/migration_from_sanghavi.md`) documenting:
- **Struct changes**: `vSmartMOM_Model` → `RTModel{ARCH,FT}` with sub-structs (`SolverConfig`, `Atmosphere`, `Optics`, `RayleighScattering`, `AerosolState`)
- **Field access changes**: e.g., `model.τ_abs` now routes through `Base.getproperty` to `model.optics.τ_abs`; `model.obs_geom` → `model.geometry`
- **Function signature changes**: `model_from_parameters()` returns `RTModel` directly; `LinMode()` variant returns `(model, lin_model)`
- **Naming convention changes**: `@unpack` → `(; field) = struct`, `J₀⁺/J₀⁻` → `j₀⁺/j₀⁻` (lowercase for AddedLayer), `FT<:Union{AbstractFloat,Dual}` → `FT<:Real`
- **New accessor functions**: `CoreRT.polarization_type(m)`, `CoreRT.float_type(m)`, `CoreRT.n_aerosols(m)`, `get_surface(m,i)`, `get_spec_bands(m)`
- **Removed types/files**: `vSmartMOM_Model`, `rt_run_bck.jl`, `types_inelastic.jl`

**TODO**: Sanghavi branch owner will provide example scripts that need conversion. These will serve as concrete migration examples showing old → new API usage.

### 1d. Sanghavi script conversion
- Convert provided sanghavi scripts to work with unified RTModel API
- Each converted script serves as a living example in the migration guide
- **Before any refactoring**: run converted scripts and capture reference Raman outputs (R, T, ieR, ieT) for regression testing

### 1e. Cleanup
- Remove `@show` debug statements from interaction kernels
- Remove commented-out code (`#bla`, `repeat(...) ⊠ reshape(...)` attempt)
- Convert any remaining `@unpack` → `(; field) = struct`

### Files to create/modify
- `docs/dev_notes/migration_from_sanghavi.md` (new — migration guide)
- `src/CoreRT/CoreKernel/rt_kernel_ss.jl` (port from sanghavi)
- `src/CoreRT/CoreKernel/interaction_ss.jl` (port from sanghavi)
- `examples/EMIT/`, `examples/Balsamic/` (new, port from sanghavi)
- Sanghavi example scripts (TBD — will be provided by branch owner)

---

## Phase 2: Eliminate Allocations in Inelastic Path (Critical Performance Fix)

**Scope: Inelastic (Raman) path only.** The elastic (noRS) path's `⊠` expressions don't appear to cause real allocations in practice (Julia's JIT optimizes them), and the `⊠` notation is valuable for physics readability. We can revisit the elastic path later if profiling shows issues.

This is the **highest-impact optimization** — eliminate `similar()` and temporary allocations from hot loops before touching the batching structure.

### 2a. Create `InelasticWorkspace` struct

Add to `src/CoreRT/types.jl`:
```julia
mutable struct InelasticWorkspace{FT, AT3, AT4}
    # Doubling temporaries (3D elastic-sized)
    gp_refl::AT3
    tt_gp::AT3
    J₁⁺::AT3      # (nQ, 1, nSpec)
    J₁⁻::AT3
    ieJ₁⁺::AT4    # (nQ, 1, nSpec, nRaman)
    ieJ₁⁻::AT4
    # Interaction temporaries (for ScatteringInterface_11)
    tmp_inv::AT3
    T_inv::AT3     # T01_inv or T21_inv
    tmpJ₀⁺::AT3;  tmpJ₀⁻::AT3
    tmpR⁻⁺::AT3;  tmpR⁺⁻::AT3
    tmpT⁺⁺::AT3;  tmpT⁻⁻::AT3
    tmpieJ₀⁺::AT4; tmpieJ₀⁻::AT4
    tmpieR⁻⁺::AT4; tmpieR⁺⁻::AT4
    tmpieT⁺⁺::AT4; tmpieT⁻⁻::AT4
    # Batched_mul output temporaries (reuse across ⊠ calls)
    buf3d_a::AT3;  buf3d_b::AT3;  buf3d_c::AT3
    # batch_inv! pivot/info (shared with elastic workspace)
    pivot::AbstractMatrix{Cint}
    info::AbstractVector{Cint}
end
```

### 2b. Allocate once in `model_from_parameters()`

In `src/CoreRT/tools/model_from_parameters.jl`, create `InelasticWorkspace` alongside `RTWorkspace` when RS_type is not `noRS`. Store it in the RT run context or pass through the pipeline.

### 2c. Refactor `doubling_inelastic.jl`

**Before:**
```julia
gp_refl = similar(t⁺⁺)           # ALLOC
ieJ₁⁺ = similar(ieJ₀⁺); ieJ₁⁺.=0  # ALLOC + ZERO
...
tmp3 = ieJ₁⁺[:,:,n₁,Δn] + (tt⁺⁺_gp_refl[:,:,n₁] ⊠ (...))  # ALLOC from ⊠
```

**After:**
```julia
(; gp_refl, ieJ₁⁺, ieJ₁⁻, buf3d_a) = ws_ie
ieJ₁⁺ .= 0  # zero pre-allocated buffer (or skip if guaranteed clean)
...
batched_mul!(buf3d_a, tt⁺⁺_gp_refl_slice, rhs_slice)  # in-place ⊠
ieJ₀⁺[:,:,n₁,Δn] .= ieJ₁⁺[:,:,n₁,Δn] .+ buf3d_a    # no alloc
```

### 2d. Refactor `interaction_inelastic.jl`

Replace all 13 `similar() + .=0` at entry of `ScatteringInterface_11` with workspace fields:
```julia
(; tmp_inv, tmpieJ₀⁺, tmpieJ₀⁻, tmpieR⁻⁺, ...) = ws_ie
tmpieJ₀⁺ .= 0; tmpieR⁻⁺ .= 0; ...  # still need zeroing, but no allocation
```

### 2e. Pass RTWorkspace to `batch_inv!` in inelastic path

Change `batch_inv!(gp_refl, I_static .- r⁻⁺ ⊠ r⁻⁺)` to use the workspace-aware version with pre-allocated pivots.

### Files to modify
- `src/CoreRT/types.jl` — add `InelasticWorkspace`
- `src/CoreRT/tools/model_from_parameters.jl` — allocate workspace
- `src/CoreRT/CoreKernel/doubling_inelastic.jl` — use workspace, eliminate `similar()`
- `src/CoreRT/CoreKernel/interaction_inelastic.jl` — use workspace, eliminate `similar()`
- `src/CoreRT/rt_run.jl` — thread workspace through RT pipeline

---

## Phase 3: Flatten 4D → 3D with Index Map (Batching Optimization)

### 3a. Create `RamanIndexMap`

In `src/Inelastic/types.jl`:
```julia
struct RamanIndexMap{IT<:AbstractVector{Int}}
    nTotal::Int                    # sum of valid (n₁, Δn) pairs
    el_n₁::IT                     # elastic index at scattered wavelength
    el_n₀::IT                     # elastic index at incident wavelength
    ranges::Vector{UnitRange{Int}} # which slice of nTotal → each Δn
end
```

Built from existing `i_λ₁λ₀` + `get_n₀_n₁()` during `getRamanSSProp!()`.

### 3b. Replace 4D inelastic arrays with 3D

In `AddedLayerRS` and `CompositeLayerRS`:
- `ier⁻⁺`: `(nQ, nQ, nSpec, nRaman)` → `(nQ, nQ, nTotal)`
- Same for all ie* fields

### 3c. Replace loops with gathered batched_mul

**Before** (nRaman sequential launches):
```julia
for Δn = 1:nRaman
    n₀, n₁ = get_n₀_n₁(...)
    tmp = tt⁺⁺_gp_refl[:,:,n₁] ⊠ (ier⁻⁺[:,:,n₁,Δn] ⊠ r⁻⁺[:,:,n₀] + ...)
end
```

**After** (1 batched launch):
```julia
gather!(ws.buf_r_n₁, r⁻⁺, idx.el_n₁)    # r⁻⁺[:,:,n₁] for all pairs
gather!(ws.buf_r_n₀, r⁻⁺, idx.el_n₀)    # r⁻⁺[:,:,n₀] for all pairs
iet⁺⁺_flat .= tt_gp_n₁ ⊠ (ier⁻⁺_flat ⊠ ws.buf_r_n₀ .+ ...) .+ ...
```

### 3d. Convert interaction 2D `*` → batched `⊠`

All `for n₁ ... for Δn ...` loops with 2D `*` become single `⊠` on gathered 3D arrays.

### Memory trade-off
Full batch only (no chunked fallback). Modern GPUs (16+ GB) can handle the extra ~1 GB for gather buffers. Add chunking later if needed.

### Files to modify
- `src/Inelastic/types.jl` — add `RamanIndexMap`
- `src/Inelastic/inelastic_helper.jl` — build index map
- `src/CoreRT/types.jl` — flatten AddedLayerRS/CompositeLayerRS
- `src/CoreRT/tools/rt_helper_functions.jl` — add `gather!`, update `default_matrix_ie()`
- `src/CoreRT/CoreKernel/doubling_inelastic.jl` — batched ops
- `src/CoreRT/CoreKernel/interaction_inelastic.jl` — batched ops
- `src/CoreRT/CoreKernel/elemental_inelastic.jl` — write to flat layout

---

## Phase 4: Reduce GPU Synchronization

- Remove redundant `synchronize_if_gpu()` inside `batch_inv!` (keep one at end of doubling step, not after every LU call)
- Batch multiple layer operations before syncing where possible
- Profile to identify remaining sync bottlenecks

### Files to modify
- `ext/gpu_batched_cuda.jl` — remove inner syncs from batch_inv!
- `src/CoreRT/CoreKernel/doubling_inelastic.jl` — consolidate sync points

---

## Phase 5: Linearized Raman (Deferred — separate PR)

No linearized Raman path exists yet. After Phases 2-4 stabilize, create `_lin` variants using the flat 3D layout from the start.

---

## Testing & Verification Strategy

### Pre-merge testing (before ANY refactoring)
1. Run all 474 tests on unified-vsmartmom: `cd test && julia --project=. runtests.jl`
2. Run sanghavi's Raman scripts on unified codebase — capture reference R, T, ieR, ieT outputs
3. Document any sanghavi scripts that fail due to struct/notation changes

### After Phase 1 (branch + sanghavi delta)
4. All tests still pass
5. Sanghavi Raman scripts produce correct output (or have been updated)
6. EMIT/Balsamic examples run under new API

### After Phase 2 (allocation elimination)
7. All tests pass — **bit-for-bit identical** results (no physics change)
8. Profile with `@allocated` or `CUDA.@time`: verify ~0 allocations inside doubling/interaction loops
9. Wall-clock speedup measurable even without batching change

### After Phase 3 (flat 3D batching)
10. All tests pass — **bit-for-bit identical** (same math, different memory layout)
11. Rayleigh = Cabannes + RRS energy conservation check
12. `@btime` comparison: expect ~nRaman× speedup in doubling/interaction on GPU
13. CUDA profiler: verify single `gemm_strided_batched` calls replacing sequential loops

### After Phase 4 (sync reduction)
14. All tests pass
15. Profile sync stalls: should be <10 per RT run instead of ~2000

---

## Implementation Priority

**Phase 2 (allocation elimination) should come FIRST** — it's the highest-impact, lowest-risk change. Even without restructuring the 4D→3D layout, eliminating `similar()` inside loops will dramatically reduce GPU overhead. This is also fully testable with bit-exact regression.

**Phase 3 (flat batching) is the algorithmic change** — higher risk, requires careful validation, but provides the ~nRaman× kernel-launch reduction.

```
Phase 1: Branch + sanghavi compat     → low risk, enables testing
Phase 2: Kill allocations             → HIGH IMPACT, low risk, bit-exact
Phase 3: Flatten 4D→3D + batch        → HIGH IMPACT, medium risk
Phase 4: Reduce sync barriers         → medium impact, low risk
Phase 5: Linearized Raman             → deferred
```

---

## Branch Comparison Summary

**What unified-vsmartmom already has from sanghavi** (commit 59f8de8):
- Core Raman physics (RRS, VS, A-band improvements)
- Rayleigh = Cabannes + RRS matching
- R⁺⁻ bug fix in ScatteringInterface_11

**What sanghavi has that unified does NOT**:
- Single-scatter approximation kernel
- EMIT/Balsamic application scripts
- Some unpushed local changes (TBD)

**What unified has that sanghavi does NOT**:
- RTModel{ARCH,FT} hierarchy with sub-structs
- Linearized RT (Jacobians) with ParameterLayout
- RTWorkspace pre-allocation (elastic path)
- Cox-Munk ocean surface, canopy surface
- GPU weak dependency via KernelAbstractions
- Float32/Float64 flexibility
- Comprehensive test suite (474 tests)
- Docstrings, tutorials, tree display
