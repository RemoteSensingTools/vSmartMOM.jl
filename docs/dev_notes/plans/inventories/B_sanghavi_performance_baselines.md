# Inventory B — sanghavi performance baselines

**As-of date:** 2026-04-19
**Purpose:** Catalogue benchmark infrastructure on both branches so that numeric
non-regression floors can be captured from the `sanghavi` branch and re-measured
on `sanghavi-unified` after merge. This inventory does **not** contain measured
numbers; it documents what must be measured and how.

## Branch tips audited

| Branch | Worktree | Tip |
|---|---|---|
| `sanghavi` (performance truth source for Raman) | `/home/sanghavi/code/github/vSmartMOM.jl/` | `9ee9a75` |
| `unified-vsmartmom` (current merge target) | `/home/sanghavi/code/github/uni_vSmartMOM/` | `a4e4187` |

Both worktrees are at the commits named in the handoff brief. No Julia runs
were performed in assembling this inventory.

---

## 1. Benchmark file inventory

### 1a. `sanghavi` worktree — `test/benchmarks/`

| File | Purpose (1-line) | Raman? |
|---|---|:---:|
| `raman_batched_ops_benchmark.jl` | Micro-benchmark: 6 tests of batched GPU ops at `NQ=9, NSPEC=5500, NRAMAN=963, NL=5000` comparing broadcasting `⊠` vs in-place `NNlib.batched_mul!`, vs 2D `mul!` loops, vs `batch_inv!` (CUBLAS), and a full doubling-like Raman loop simulation. Float64. GPU-only (`device!(1)`). | ✅ |
| `raman_optimization_baseline.jl` | End-to-end reference generator: runs RRS + noRS on 1-band `O2_parameters2_1band_opttest.yaml`, saves `R,T,ieR,ieT` to `raman_opttest_output/raman_reference.jld2` for bit-exact comparison. Records `t_rrs` and `t_nors`. Float64. GPU (`device!(1)`). | ✅ |
| `raman_optimization_compare.jl` | Regression check: reloads the reference JLD2, re-runs RRS + noRS, reports bit-exact match on all six arrays, and prints wall-clock speedup relative to `timing_baseline.txt`. Float64, GPU. | ✅ |
| `raman_opttest_output/timing_baseline.txt` | Pre-existing snapshot (`RRS = 1671.851 s`, `noRS = 2.836 s`) from a previous sanghavi run. **Not authoritative** — must be re-captured on current `9ee9a75` to anchor the merge floor. | ✅ |
| `O3Huggins_polRaman.jl` | Polarized O3 Huggins-band RRS retrieval prototype (286 Raman references). | ✅ |
| `RamanSIFmaps.jl`, `RamanSIFspectra.jl`, `RamanSIFworldmaps.jl`, `RamanSIFPolyCoeffMaps.jl`, `RamanSIF_oco_retr.jl` | SIF retrieval scripts leveraging Raman filling-in of Fraunhofer lines. | ✅ |
| `creategrid_O2Aband_RamanSIF.jl`, `creategrid_O2Bband_RamanSIF.jl`, `creategrid_betwnAB_RamanSIF.jl`, `test_creategrid_O2Aband_RamanSIF.jl`, `test_creategrid_noabs_RamanSIF.jl` | LUT-grid builders for the OCO O2 A/B and inter-band retrievals with Raman. | ✅ |
| `evalAer_O2Aband_Raman.jl`, `evalAer_O2ABbands_Raman.jl`, `testAer_O2Aband_Raman.jl` | Aerosol sensitivity evaluations with Raman on. | ✅ |
| `prototype_O2ABand_RRS.jl`, `prototype_O2ABand_RRS_SIF.jl`, `prototype_O2BBand_RRS_SIF.jl` | RRS prototypes for the OCO-2 A/B bands. | ✅ |
| `prototype_VS_O2Aband.jl`, `prototype_fraunhofer_VS_spectrum.jl`, `prototype_inelastic.jl`, `prototype_inelastic_OCO2.jl`, `prototype_inelastic_VS2.jl`, `prototype_inelastic_ms.jl`, `prototype_inelastic_sunspot_VS.jl` | Vibrational-scattering (VS) and inelastic prototypes. | ✅ |
| `sif_raman.jl`, `testRayCabRaman.jl`, `plot_ramanXS.jl` | Raman + Cabannes helpers and plotting. | ✅ |
| `compLelli.jl`, `compLelli2.jl` | Lelli (Coulson) slab reference intercomparison. | ❌ |
| `compare_rt_EMIT.jl`, `emit_modtran_noRS_scenarios.jl`, `prototype_EMIT_aer_ht.jl` | EMIT noRS runs + MODTRAN comparison. `emit_modtran_noRS_scenarios.jl` runs on GPU (`device!(1)`). | ❌ |
| `prototype_CarbonI.jl` | CarbonI instrument prototype (no Raman grep match). | ❌ |
| `creategrid_CO2Wband_OCORayl.jl`, `creategrid_O2Aband_OCORayl.jl` | OCO Rayleigh-only LUT builders. | ❌ |
| `benchmark_6SV1.jl`, `6SV1_1.yaml`, `6SV1_R_trues.jl` | 6SV1 intercomparison. | ❌ |
| `benchmark_natraj.jl`, `natraj.yaml`, `natraj_I.yaml`, `natraj_v2.yaml`, `natraj_trues.jl`, `plot_MODvsMOM.jl`, `TestSZA.jl` | Natraj/plane-parallel standards. | ❌ |
| `testAerosol.jl`, `testCPU.jl` | Developer sanity scripts. | ❌ |
| `create_HITRAN_LUTs.jl` | HITRAN LUT builder. | ❌ |
| `RT_vs_hydrolight.jl` | Hydrolight water-leaving radiance comparison. | ❌ |
| `Manifest.toml`, `Project.toml`, `test.ipynb`, `plot_spec_RamanSIF.py` | Env/notebook/plotting artifacts. | — |

### 1b. `unified-vsmartmom` worktree — `test/benchmarks/` + `test/`

| File | Purpose (1-line) | Raman? |
|---|---|:---:|
| `test/benchmarks/raman_large_scale.jl` | Large-scale Raman memory strategy bakeoff at `NquadN=30, nSpec=30000, nRaman=100, ndoubl=8, FT=Float64`. Compares five data-placement strategies for the 4D ie_r/ie_t tensors: (A) full 4D on GPU, (B) `Vec{CuArray3D}` all on GPU, (C) CPU-pinned 4D + per-Δn H2D/D2H, (D) CPU-pinned Vec3D + per-Δn copies, (E) tiled (tile=4,16,32); also runs (F) Vec{CuArray3D} at Float32 for the FP32/FP64 speedup. Skips Strategy A automatically if it would exceed 85% of free VRAM. | ✅ |
| `test/benchmarks/raman_memory_benchmark.jl` | Toy mimic of the inelastic doubling/interaction allocation + compute pattern. `PRESETS = {small, medium, large, production}` parameterize `NquadN=15, nSpec∈{68,500,2000,5424}, nRaman∈{20,80,172}, ndoubl=8`. Individual sub-benchmarks: (1) 4D allocation cost, (2) batched matmul throughput (incl. 3D⊠4D-view vs copy-then-mul), (3) CPU↔GPU transfer vs on-GPU 4D access. | ✅ |
| `test/benchmarks/batched_fused_benchmark.jl` | Custom KA kernel for fused `X = B * inv(I − A₁A₂)` (LU + parallel triangular solve) vs cuBLAS 3-step. Tests doubling + adding patterns. Not Raman-specific but exercises the same 15×15 shape. | ⚪ (shape-related) |
| `test/benchmarks/batched_fused_v2_benchmark.jl` | v2 of above: direct-solve `X*M=B` (never forms M⁻¹) vs v1 inv+mul vs cuBLAS. | ⚪ |
| `test/benchmarks/batched_inv_benchmark.jl` | cuBLAS `getrf+getri` vs custom KA Gauss-Jordan vs custom KA LU+solve, single-kernel batched inversion. | ⚪ |
| `test/benchmarks/batched_matmul_benchmark.jl` | cuBLAS `gemm_strided_batched` vs custom KA batched matmul (element-per-thread and shared-memory variants). | ⚪ |
| `test/benchmarks/interaction11_bundled_benchmark.jl` | End-to-end `ScatteringInterface_11` forward update: explicit-inverse path vs algebraically reduced cuBLAS path vs bundled KA single-kernel reduction. Elastic-only. | ⚪ |
| `test/benchmarks/prototype_EMIT_aer_ht.jl` | EMIT aerosol height prototype. | ❌ |
| `test/benchmarks/6SV1_1.yaml`, `6SV1_1_simple.yaml`, `6SV1_R_trues.jl`, `natraj*.yaml`, `natraj_trues.jl` | Standard intercomparison inputs. | ❌ |
| `test/benchmark_rt.jl` | Top-level benchmark: forward noRS, forward RRS (O2 A-band), linearized RT; CPU and optional GPU. Uses `test_parameters/PureRayleighParameters.yaml` + others. | ✅ (RRS) |
| `test/benchmark_fwd_vs_lin.jl` | Forward vs linearized RT full-pipeline, CPU + GPU, O2Parameters.yaml, 5 timed runs. | ❌ (elastic only in this file) |

### 1c. Cross-cutting observations

- The sanghavi branch has a **published regression harness**
  (`raman_optimization_baseline.jl` + `raman_optimization_compare.jl`) with a
  JLD2 reference and a bit-exact-match checker. This is the nearest thing to a
  non-regression harness today. It is **single-benchmark** (1-band O2 opttest)
  and **wall-clock only** — it does not record allocations or peak GPU memory.
- The unified branch has **no end-to-end Raman regression harness** in the same
  sense; `raman_large_scale.jl` and `raman_memory_benchmark.jl` are isolated
  microbenchmarks with synthetic data, and `benchmark_rt.jl` times the noRS/RRS
  forward pipeline without persistent baselines.
- Neither branch currently records **GPU allocation counts** or **peak GPU
  memory** persistently. The scripts print to stdout only; there is no JSON/CSV
  output directory aside from `raman_opttest_output/`.
- The unified branch has Raman-relevant micro-benchmarks that the sanghavi
  branch does not (and vice versa). For cross-branch regression capture, both
  sets must be able to run on `sanghavi-unified` post-merge.

---

## 2. Quantities to measure

For each acceptance benchmark, capture the following on **both** `sanghavi` and
`sanghavi-unified` (after merge), on the same GPU, same Julia version, same
CUDA.jl version:

### 2a. Per-run metrics

| Metric | Source | Notes |
|---|---|---|
| Wall-clock (warm) | `@elapsed` after 1 warmup + `CUDA.synchronize()` | Report min / median / max of N≥3 runs. |
| CPU allocation count & bytes | `@allocations`, `@allocated` | Must wrap the same expression as the timed run. |
| GPU memory pool — used delta | `CUDA.Mem.info()` before/after | `(total − free)` difference across the measured region. |
| GPU memory pool — peak | `CUDA.pool_status()` **or** `CUDA.@profile` | Prefer `CUDA.@profile trace=true` for kernel-level peak. |
| GPU kernel-launch count | `CUDA.@profile` | Sanity check that optimizations reduce launches, not just bytes. |
| Correctness hash | bit-exact vs saved JLD2 reference | Already implemented in `raman_optimization_compare.jl`; extend to other benchmarks. |

### 2b. Problem-size axes (explicit for each script)

| Axis | `raman_batched_ops` | `raman_large_scale` | `raman_memory` (production) | `raman_optimization` | `benchmark_rt` (RRS) |
|---|:---:|:---:|:---:|:---:|:---:|
| `n_quad` (NquadN) | 9 | 30 | 15 | derived from yaml | derived from yaml |
| `n_stokes` | 3 (IQU, implicit in NQ=9) | — | — | from yaml | from yaml |
| `n_spec` | 5500 | 30000 | 5424 | ~1600 (opttest) | ~6800 (O2 A-band) |
| `n_layers` | — | — | 12 | profile-driven | profile-driven |
| `ndoubl` | implicit | 8 | 8 | solver-default | solver-default |
| `n_Raman` | 963 (sim: 100) | 100 | 172 | from `getRamanSSProp!` | from `getRamanSSProp!` |
| `RS_type` | synthetic (no RS_type object) | synthetic | synthetic | `RRS` + `noRS` | `RRS` + `noRS` |
| `float_type` | Float64 | Float64 **and** Float32 | Float64 (preset-driven) | from yaml | from yaml |
| Architecture | GPU only | GPU only | GPU only | GPU only | CPU + GPU |

### 2c. Warm-up / measure / record pattern (recommended)

```julia
# 1. Load config / build model / build all GPU buffers
setup(...)

# 2. Warmup: one full run, discard result
run(...)
CUDA.synchronize()
CUDA.reclaim()

# 3. Measure N times (N ≥ 3), capturing each metric:
for i in 1:N
    GC.gc(); CUDA.reclaim()
    m0 = CUDA.Mem.info()
    CUDA.synchronize()
    t = @elapsed begin
        result = run(...)
        CUDA.synchronize()
    end
    m1 = CUDA.Mem.info()
    alloc_cpu = @allocated run(...)  # separate invocation to avoid
                                     #   perturbing timed one
    # store: (i, t, alloc_cpu, m0, m1)
end

# 4. Separate `CUDA.@profile` invocation for kernel counts/peak GPU mem
CUDA.@profile trace=true begin
    run(...)
end
```

This matches the structure already used in `raman_optimization_compare.jl` and
`raman_batched_ops_benchmark.jl` (both perform 1 warmup + timed run with
`CUDA.reclaim()` between). It does **not** currently capture GPU memory or
allocation counts systematically — extending it is the harness work in §3.

---

## 3. Recommended harness (thin)

### 3a. Location

Create a harness package rooted at:

```
test/benchmarks/harness/
├── run_benchmarks.jl        # entry point, CLI arg → which benchmark(s) to run
├── metrics.jl               # capture helpers (wall, alloc, GPU mem, peak)
├── report.jl                # JSON/CSV writer
└── scenarios.toml           # declarative list of acceptance benchmarks
```

This keeps the existing benchmark files untouched (they remain directly
runnable) and lets the harness `include()` them while capturing metrics.

### 3b. Output directory

```
test/benchmarks/baseline_output/
├── sanghavi_9ee9a75/
│   ├── raman_optimization.json
│   ├── raman_batched_ops.json
│   ├── raman_large_scale.json
│   ├── benchmark_rt_rrs.json
│   └── summary.csv
└── unified_merged_<sha>/
    └── ... (same layout)
```

Committing the JSON/CSVs under `baseline_output/` gives a reviewable diff when
`sanghavi-unified` is benchmarked post-merge.

### 3c. Pseudo-code sketch

```julia
# test/benchmarks/harness/metrics.jl
struct RunMetrics
    name::String
    problem_size::Dict{Symbol,Any}
    wall_s::Vector{Float64}          # per-run
    alloc_cpu_bytes::Vector{Int}
    gpu_used_delta_gb::Vector{Float64}
    gpu_peak_mib::Float64            # from CUDA.@profile
    n_kernels::Int                   # from CUDA.@profile
    passes_correctness::Bool         # vs reference JLD2
    branch::String
    commit::String
    gpu_name::String
    julia_version::String
    ts::DateTime
end

function capture(name, setup_fn, run_fn; nruns=3, ref_path=nothing)
    setup_fn()                       # model build, GPU buffers, etc.
    run_fn()                         # warmup
    CUDA.synchronize(); CUDA.reclaim()

    walls = Float64[]; allocs = Int[]; mem_deltas = Float64[]
    for i in 1:nruns
        GC.gc(); CUDA.reclaim()
        free0, _ = CUDA.Mem.info()
        CUDA.synchronize()
        t = @elapsed (run_fn(); CUDA.synchronize())
        free1, _ = CUDA.Mem.info()
        push!(walls, t)
        push!(mem_deltas, (free0 - free1) / 1e9)
        push!(allocs, @allocated run_fn())
    end

    # Profile invocation for kernel + peak
    prof = CUDA.@profile trace=true run_fn()
    peak_mib = extract_peak(prof)    # parse @profile output
    nk      = count_kernels(prof)

    # Optional bit-exact check
    passes = ref_path === nothing ? true : check_reference(run_fn(), ref_path)

    return RunMetrics(name, ..., walls, allocs, mem_deltas,
                      peak_mib, nk, passes, ...)
end

# test/benchmarks/harness/run_benchmarks.jl
scenarios = TOML.parsefile("scenarios.toml")["scenarios"]
for sc in scenarios
    include(sc["script"])             # exposes setup / run / ref_path
    m = capture(sc["name"], setup, run; nruns=sc["nruns"], ref_path=sc["ref"])
    write_json(joinpath(OUT_DIR, sc["name"] * ".json"), m)
end
```

### 3d. Running it

```bash
# Baseline (on sanghavi worktree):
cd /home/sanghavi/code/github/vSmartMOM.jl
julia --project=. test/benchmarks/harness/run_benchmarks.jl \
      --out test/benchmarks/baseline_output/sanghavi_9ee9a75/ \
      --branch sanghavi --commit 9ee9a75

# Regression check (on sanghavi-unified worktree, post-merge):
cd /home/sanghavi/code/github/uni_vSmartMOM
julia --project=. test/benchmarks/harness/run_benchmarks.jl \
      --out test/benchmarks/baseline_output/unified_merged_<sha>/ \
      --branch sanghavi-unified --commit <sha>
julia --project=. test/benchmarks/harness/compare.jl \
      --baseline test/benchmarks/baseline_output/sanghavi_9ee9a75/summary.csv \
      --candidate test/benchmarks/baseline_output/unified_merged_<sha>/summary.csv
```

`compare.jl` should emit pass/fail per scenario; tolerance selection is Stage 2
and explicitly out of scope for this inventory.

---

## 4. Baseline capture checklist

Run these (in this order) on the sanghavi worktree at `9ee9a75`. All Julia
commands assume `cd /home/sanghavi/code/github/vSmartMOM.jl`.

- [ ] **Record environment**: Julia version, CUDA.jl version, NVIDIA driver,
      GPU name (`CUDA.name(CUDA.device())`), free/total VRAM. Write to
      `baseline_output/sanghavi_9ee9a75/env.txt`.
- [ ] **Confirm GPU 1 is free** (matches `device!(1)` in the scripts). If a
      different device index is required on the target machine, update the
      scripts or parameterize via ENV before the run.
- [ ] **Regenerate the 1-band reference** (since `raman_reference.jld2` may be
      stale for `9ee9a75`):
      `julia --project=. test/benchmarks/raman_optimization_baseline.jl`
- [ ] **Measure end-to-end 1-band RRS + noRS** (3 warm runs):
      `raman_optimization_compare.jl` with the harness wrapper, capture wall +
      alloc + GPU peak. Record `RRS (s)`, `noRS (s)`, `ratio RRS/noRS`.
- [ ] **Measure batched-ops microbench** (6 tests × sizes as coded):
      `raman_batched_ops_benchmark.jl` via the harness. Record every `t*` and
      `alloc*` printed as a numeric row per test.
- [ ] **Measure full O2 A-band RRS forward** from `benchmark_rt.jl` (CPU and
      GPU), captured via the harness. This is the "real" per-band runtime
      claim, distinct from the 1-band opttest.
- [ ] **Dump `raman_large_scale.jl` results** for Strategies A–F on the same
      GPU to set a memory-strategy baseline for the large-grid case.
- [ ] **Persist** all JSON + the human-readable `timing_baseline.txt` under
      `baseline_output/sanghavi_9ee9a75/` and commit to the `sanghavi-unified`
      branch once the merge begins (so regression diffs live in-repo).

Once the checklist is filled, the numeric "non-regression floor" is frozen.
Re-run the identical harness on `sanghavi-unified` after each merge milestone.

---

## 5. Open questions

1. **Authoritative allocation count**: Christian's `raman_gpu_optimization.md`
   quoted ~19k GPU allocations/run, measured on the older unified snapshot.
   **Do not quote that as an authoritative floor.** The sanghavi branch with
   `InteractionWorkspace` should measurably reduce allocations vs both the
   old unified and current `unified-vsmartmom` tips. Both branches need fresh
   `CUDA.@profile` numbers before merge; until those land, treat allocation
   counts as "to-be-measured" and not a pass/fail gate.
2. **GPU device selection**: multiple scripts hardcode `device!(1)`. Decide
   whether the harness overrides via ENV (recommended) or leaves the
   hardcoded device — this matters if the target machine has a different
   topology.
3. **`benchmark_rt.jl` RRS path**: it exists on the unified worktree and
   advertises an RRS case, but the corresponding sanghavi worktree
   `test/benchmarks/` does not have the same file. We should either port
   `benchmark_rt.jl` to the sanghavi worktree for fair baseline capture, or
   accept that its baseline will only exist on the unified side.
4. **Float32 scope**: `raman_large_scale.jl` has an explicit FP32 arm; the
   sanghavi regression harness is FP64-only. Decide whether the non-regression
   floor needs separate FP32/FP64 tracks (the float32 investigation memory
   thread suggests yes, but that may be a Stage-2 decision).
5. **Reference freshness**: `raman_opttest_output/raman_reference.jld2` and
   `timing_baseline.txt` predate `9ee9a75`. They must be regenerated before
   capture; otherwise the bit-exact check in `raman_optimization_compare.jl`
   will fail for non-performance reasons (code churn since the last snapshot).
6. **CPU vs GPU coverage**: `benchmark_rt.jl` and `benchmark_fwd_vs_lin.jl`
   run on both; the Raman-specific benchmarks run GPU-only. Decide whether
   the non-regression floor includes CPU for noRS (probably yes) and CPU for
   RRS (probably no — too slow).
7. **Peak GPU mem capture**: `CUDA.pool_status()` prints to stdout and is
   per-invocation; `CUDA.@profile trace=true` yields a richer trace. Preferred
   source needs to be fixed before harness implementation; suggest
   `CUDA.@profile` as the authoritative source with `pool_status()` as a
   sanity check.
