# AI Handoff Memo — vSmartMOM.jl (Jacobians & Tests)

**Purpose:** Handoff for another AI tool to continue work on analytic Jacobians, tests, AD comparison, and GPU support.

---

## Repo & environment

- **Repo:** https://github.com/RemoteSensingTools/vSmartMOM.jl  
- **Branch:** `unified-vsmartmom`  
- **Language:** Julia (1.9+). Use project: `julia --project=.`  
- **Setup:** `julia --project=. -e 'using Pkg; Pkg.instantiate()'`

---

## What this branch does

vSmartMOM.jl is a radiative transfer (RT) code. This branch adds **analytic Jacobians** (dR/d params) for TOA radiance R. Parameters include: τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp (per aerosol), gas VMR, surface albedo.

- **Forward RT:** `model_from_parameters(params)`, then `rt_run(model, ...)` (no `LinMode()`).
- **Linearized RT (analytic dR):** `model_from_parameters(LinMode(), params)` → `rt_run(model, lin_model, ...)` returns `R, T, dR, dT`.

---

## Test suite (run in this order)

| Command | What it does | Status |
|--------|----------------|--------|
| `julia --project=. test/test_jacobians_unit.jl` | Analytic vs FD for albedo, τ_ref, p₀, nᵣ, σ_dist; finite checks | ✅ 29 tests pass |
| `julia --project=. test/test_jacobians_AD_compare.jl` | Table: analytic vs central FD for τ_ref, albedo, nᵣ; timing | ✅ Runs (no AD column) |
| `julia --project=. test/test_jacobians_GPU.jl` | Linearized RT on GPU; GPU vs CPU | ⚠️ Skips: scalar indexing in doubling |

- **Config:** All use `test/test_parameters/JacobianTestFast.yaml` (4 wavelengths, 1 aerosol).
- **More:** `docs/JACOBIAN_TEST_WORKFLOW.md`, `docs/SESSION_HANDOFF.md`.

---

## Key paths (for edits / search)

- **Tests:** `test/test_jacobians_unit.jl`, `test/test_jacobians_AD_compare.jl`, `test/test_jacobians_GPU.jl`
- **Linearized RT entry:** `model_from_parameters(LinMode(), params)` and `rt_run(..., lin_model, ...)` — see `src/CoreRT/rt_run_lin.jl`, `src/CoreRT/model_lin.jl` (or similar).
- **RT kernel (linearized):** `src/CoreRT/CoreKernel/rt_kernel_lin.jl` — elemental → chain_rule → doubling_allparams! → interaction.
- **Elemental (linearized):** `src/CoreRT/CoreKernel/elemental_lin.jl` — kernel `get_elem_rt!` uses `FT = eltype(r⁻⁺)` (device-safe).
- **Doubling (linearized):** `src/CoreRT/CoreKernel/doubling_lin.jl` — `doubling_allparams_helper!`; **GPU scalar indexing** happens here (batched_mul / views of 4D CuArrays).
- **CUBLAS (GPU):** `CoreRT.CUBLAS_ref` in `src/CoreRT/constants.jl`; set in `ext/vSmartMOMCUDAExt.jl`. Used in `src/CoreRT/tools/rt_helper_functions_lin.jl` and `rt_helper_functions.jl`.
- **Bugs / theory:** `docs/LINEARIZATION_BUGS.md`, `docs/SESSION_HANDOFF.md`.

---

## Current limitations (for next AI)

1. **AD (ForwardDiff) in AD compare test**  
   - AD is **not** run in `test_jacobians_AD_compare.jl` because the **linearized** path (`model_from_parameters(LinMode(), params)`) builds Float64-only arrays; Dual-valued params cause `MethodError`.  
   - **Desired:** Run AD on the **forward-only** path (no `LinMode()`): wrap a function that takes a scalar param (e.g. τ_ref), builds `params` with that param set to the scalar (or to a `Dual`), calls `model_from_parameters(params)` (no LinMode) and `rt_run(model, ...)`, returns `R[1,1,1]`. Then `ForwardDiff.derivative(f, x0)` will work. Compare AD column to analytic (and FD) in the same test script; document that AD uses forward model, analytic/FD use linearized.

2. **GPU scalar indexing**  
   - Linearized RT on GPU fails with “Scalar indexing is disallowed” / “scalar indexing of a GPU array” in `doubling_allparams_helper!` (in `doubling_lin.jl`), around batched_mul / SubArray of 4D CuArray.  
   - **Task:** Find the exact indexing that triggers scalar access (e.g. a loop over a dimension that indexes into a GPU array, or a view that forces host-side iteration) and replace with a GPU-safe pattern (e.g. batched ops, or a small kernel that avoids host indexing).

3. **Optional: FT = Float64 hardcoding**  
   - User noted that “in the past we just used FT types that could either be floats or duals.” If you add AD on the forward path, check that the **forward** model uses generic `FT` (or `AbstractFloat`/Dual-friendly types) and is not hardcoding `Float64` so that ForwardDiff can propagate Duals.

---

## Suggested first steps for next AI

1. **Re-enable AD in the AD compare test (forward model):**  
   - In `test/test_jacobians_AD_compare.jl`, add a wrapper that runs the **forward** RT only (no `LinMode()`): e.g. `run_fwd_R_forward(params)` that builds `model = model_from_parameters(params)`, runs `rt_run(model, ...)` (signature without `lin_model`), returns `R`.  
   - For each param (τ_ref, albedo, nᵣ), define `f(x)` that sets that param to `x`, calls `run_fwd_R_forward(params)` and returns `R[1,1,1]`.  
   - Call `ForwardDiff.derivative(f, x0)` and add an “AD” column to the table. Keep analytic and FD as they are (linearized path). Note in comments/output that AD is via forward model, analytic/FD via linearized.

2. **Locate GPU scalar indexing in doubling_lin.jl:**  
   - Search for `getindex`, `view`, `@views`, or loops over indices that slice 4D arrays (e.g. `ap_ṫ⁺⁺[iparam,:,:,:]`) in `doubling_allparams_helper!` and the batched_mul call site.  
   - The stack trace points to `doubling_lin.jl:246` and batched_mul with `SubArray{..., CuArray{Float64, 4}, ...}`. Replace the pattern that causes host-side indexing (e.g. loop over params that does scalar/row indexing) with a single batched op or a kernel.

---

## One-line summary

**Repo:** vSmartMOM.jl, branch `unified-vsmartmom`. Jacobian unit tests pass; AD compare runs analytic vs FD only; GPU test skips due to scalar indexing in `doubling_lin.jl`. Next: add AD column using forward-model wrapper; fix GPU scalar indexing in `doubling_allparams_helper!`.
