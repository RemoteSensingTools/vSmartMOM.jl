# Jacobian Test Workflow

This document describes how to run the full Jacobian test suite (unit tests, AD comparison, and GPU).

## Prerequisites

- Julia 1.9+ with project instantiated: `julia --project=. -e 'using Pkg; Pkg.instantiate()'`
- For GPU tests: CUDA-capable device and CUDA.jl (optional; tests skip if unavailable)

## Test Scripts

| Script | Purpose | Typical time |
|--------|---------|----------------|
| `test/test_jacobians_unit.jl` | Unit tests: analytic vs FD for surface, τ_ref, p₀, nᵣ, σ_dist | ~5–8 min |
| `test/test_jacobians_AD_compare.jl` | Compare analytic vs FD (central); report timing | ~5–10 min |
| `test/test_jacobians_GPU.jl` | Same linearized RT on GPU; finite check + GPU vs CPU consistency | ~5–8 min |

## Recommended order

```bash
# 1. Unit tests (analytic vs finite difference)
julia --project=. test/test_jacobians_unit.jl

# 2. AD comparison (analytic vs FD only; timing)
julia --project=. test/test_jacobians_AD_compare.jl

# 3. GPU tests (if you have CUDA)
julia --project=. test/test_jacobians_GPU.jl
```

## What each test does

### test_jacobians_unit.jl

- Runs one linearized RT with `JacobianTestFast.yaml` (4 wavelengths, CPU).
- Compares analytic `dR` to one-sided or central FD for:
  - Surface albedo (strict tolerance)
  - τ_ref, nᵣ, σ_dist (relaxed; p₀ has known catastrophic cancellation, see Bug 23 in `LINEARIZATION_BUGS.md`).
- Ensures forward R and dR are finite.

### test_jacobians_AD_compare.jl

- Runs linearized RT once, then for selected parameters (τ_ref, albedo, nᵣ):
  - **Central FD**: same code path as analytic (linearized); two runs per parameter.
- Prints a table: analytic vs FD and relative error (analytic vs FD).
- Prints timing for FD total. Analytic gives all columns in one run; FD needs two runs per parameter.
- **Note:** ForwardDiff (AD) through the linearized model is not supported: `model_from_parameters(LinMode(), params)` expects Float64; Dual-valued params cause MethodError in array assignment.

### test_jacobians_GPU.jl

- Loads same YAML, sets `params.architecture = GPU()`.
- Runs linearized RT on GPU; checks R and dR are finite.
- Runs same config on CPU and compares R and dR with `rtol=1e-4`.
- Skips with `@test_skip` if CUDA is not available or not functional.

## Known limitations

- **p₀ (and σp)**: Large relative errors and wrong sign with surface albedo > 0 due to catastrophic cancellation (Bug 23). Unit test uses relaxed tolerance; AD comparison can include p₀ with `run_comparison(; skip_p₀=false)` (slower).
- **τ_ref with surface**: ~1–2% error from quadrature discretization; τ_ref at albedo=0 matches FD very well.
- **GPU**: CUBLAS is wired via `CoreRT.CUBLAS_ref`. The elemental kernel `get_elem_rt!` uses `FT = eltype(r⁻⁺)` so it compiles on device. Full linearized RT on GPU still fails with "Scalar indexing is disallowed" in `doubling_allparams_helper!` (batched_mul / views of 4D CuArrays). The GPU test skips with a message when it hits this. CPU linearized RT is the reference.

## References

- Full bug list and fixes: `docs/LINEARIZATION_BUGS.md`
- Session handoff and file map: `docs/SESSION_HANDOFF.md`
