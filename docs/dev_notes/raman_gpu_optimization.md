# Raman GPU Optimization Analysis

## Date: 2026-03-21

## Problem Summary

The Raman (inelastic scattering) path has severe GPU performance issues compared to the elastic (noRS) path. The elastic path uses pre-allocated workspace and single batched operations; the Raman path allocates memory inside hot loops and uses sequential operations.

## Specific Issues Found

### 1. `similar()` + `.= 0` inside hot functions

**`interaction_inelastic.jl:293-305`** (ScatteringInterface_11) allocates **13 arrays** via `similar()` and zeros them on every call (~640 calls/run):

```julia
tmp_inv   = similar(t⁺⁺); tmp_inv.=0;
tmpieJ₀⁻ = similar(ieJ₀⁻); tmpieJ₀⁻.=0;
tmpieR⁻⁺ = similar(ieR⁻⁺); tmpieR⁻⁺.=0;
# ... 10 more
```

**Impact**: ~8,320 GPU `cudaMalloc`/`cudaFree` cycles + 8,320 zero-fill kernel launches per run.

**Elastic comparison**: `doubling.jl:21-22` uses pre-allocated `temp1`, `temp2`, `dbl_gp_refl` from AddedLayer — zero allocations.

### 2. Temporary arrays from `⊠` in nested loops

**`doubling_inelastic.jl:62-81, 97-110`** — Inside `for Δn = 1:nRaman` nested inside `for n = 1:ndoubl`:

```julia
tmp3 = ieJ₁⁺[:,:,n₁,Δn] + (tt⁺⁺_gp_refl[:,:,n₁] ⊠ (...))
```

Each `⊠` returns a new array. Expression `A ⊠ (B + C ⊠ D)` creates 3 temporaries. With nRaman=10, ndoubl=3: ~360 allocations just in doubling SFI.

### 3. `batch_inv!` allocates pivot arrays every call

**`ext/gpu_batched_cuda.jl:34,63-64`**:
```julia
pivot, info = CUDA.CUBLAS.getrf_strided_batched!(A, true)  # allocates each call
```

The workspace-aware version exists (`batch_inv!(X, A, ws::RTWorkspace)`) but the inelastic path never uses it.

### 4. 2D `*` instead of batched `⊠` in interaction

**`interaction_inelastic.jl:70-86`** uses scalar matrix multiply `*` on 2D slices inside nested `for n₁ ... for Δn` loops — nSpec × nRaman sequential BLAS calls instead of one batched call.

### 5. 4D loop-and-slice prevents GPU batching

Inelastic arrays are 4D `(nQ, nQ, nSpec, nRaman)`. All operations loop over `Δn = 1:nRaman` and operate on 3D slices — this is nRaman sequential GPU kernel launches instead of 1.

### 6. Excessive GPU synchronization

`batch_inv!` has 2 sync barriers per call. Doubling syncs after every step. ~2000+ stalls/run.

## Estimated allocation counts per run

| Source | Allocs/run |
|--------|-----------|
| interaction_helper! `similar()` | ~8,300 |
| doubling tmp3/4/5/6 in Δn loop | ~3,600 |
| batch_inv! pivot/info | ~2,000 |
| ⊠ temporaries in expressions | ~5,000+ |
| **Total** | **~19,000** |

## Proposed fix priority

1. **Kill allocations** — Create `InelasticWorkspace` with pre-allocated buffers (like elastic `RTWorkspace`). Highest impact, lowest risk, bit-exact results.
2. **Flatten 4D→3D** — Precompute `RamanIndexMap`, use gather+batched_mul. ~nRaman× fewer kernel launches.
3. **Reduce syncs** — Consolidate `synchronize_if_gpu()` calls.
