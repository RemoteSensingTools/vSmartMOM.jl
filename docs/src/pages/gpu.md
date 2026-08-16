# Run on GPU

**For:** users who want to run the solver on GPU backends.

**Next:** [Quick Start](quickstart.md), [IO API](IO/Overview.md), [Architecture-Agnostic Code (Concepts)](concepts/07_architecture.md).

Every RT array in vSmartMOM has shape `(NquadN, NquadN, nSpec)` — the
spectral axis is the batch axis, and the elemental/doubling/interaction
kernels run as one batched matrix operation per layer per Fourier moment.
GPU support attaches to that same kernel source via
[KernelAbstractions.jl](https://juliagpu.github.io/KernelAbstractions.jl)
and Julia package extensions, with no second code path. See
[Architecture-Agnostic Code (Concepts)](concepts/07_architecture.md) for the
full story with file:line evidence; this page is the practical how-to.

## Selecting an architecture

Set `radiative_transfer.architecture` in the scene configuration:

```yaml
radiative_transfer:
  architecture: Architectures.CPU()       # default; always available
  # architecture: Architectures.GPU()       # CUDA
  # architecture: Architectures.MetalGPU()  # Apple Silicon
```

or override an already-loaded `vSmartMOM_Parameters` in a script:

```julia
params.architecture = vSmartMOM.Architectures.GPU()
```

`default_architecture()` picks `GPU()` if a functional CUDA device was
detected when the CUDA extension loaded, `MetalGPU()` if Metal was
detected, and `CPU()` otherwise; it is used when `architecture` is left
unset.

## CUDA vs Metal

Both backends attach through Julia package extensions
(`ext/vSmartMOMCUDAExt.jl`, `ext/vSmartMOMMetalExt.jl`): the dispatch
methods `Architectures.devi`, `Architectures.array_type`, and
`Architectures.architecture` are only defined for `GPU()` / `MetalGPU()`
once the matching backend package has been loaded in the same session
(`using CUDA` or `using Metal`). Calling `GPU()`/`MetalGPU()` without the
matching package loaded throws an actionable error (`Architectures.ka_backend`,
`default_mie_precision_policy`) rather than failing inside a kernel launch.

- **CUDA (`GPU()`)** is the mature, most-tested path: CUBLAS-backed batched
  matmul/inverse, a GPU Mie pipeline (NAI2), and ForwardDiff `Dual` support
  through the batched-matmul kernels for linearized (Jacobian) runs.
- **Metal (`MetalGPU()`)** is a first-pass Apple Silicon path: a portable
  KernelAbstractions matmul + LU-with-partial-pivoting kernel, no vendor
  BLAS. The batched-inverse kernel uses `@localmem` threadgroup memory
  with a conservative 32 KiB guard, so `Float32` matrices with
  `N = Nquad * n_stokes ≥ 64` are rejected with a clear local-memory error
  instead of a driver crash.

## Float32 vs Float64

- `MetalGPU()` requires `Float32` scene parameters (`float_type: Float32`) —
  Metal has no Float64 hardware support at all. `default_mie_precision_policy`
  auto-selects `NativeFloat32` (pure Float32, never Float64-widened) for
  Float32 inputs on Metal, and throws a clear `ArgumentError` for Float64
  inputs rather than silently failing inside a device-array allocation.
- `GPU()` (CUDA) supports both. The GPU Mie precision policy auto-selects
  `NativeFloat64` for `Float64` models and a Float32-native double-single
  path (`DSEmulated`) for `Float32` models (`default_mie_precision_policy`).
  Consumer/F32-throughput-limited CUDA GPUs can opt into `NativeFloat32`
  explicitly (`precision_policy = Scattering.NativeFloat32()`) for the
  caller-node Mie API (`compute_aerosol_optical_properties_nodes`/
  `compute_aerosol_extinction_nodes`) — end-to-end Float32, zero Float64
  device arrays, at a further accuracy cost vs `DSEmulated` (see
  `vSmartMOM.Scattering.MiePrecisionPolicy` and `test/test_mie_nodes.jl` for
  measured error bands).
- `CPU()` supports both, with no precision-policy distinction.

### `NativeFloat32` measured accuracy (caller-node Mie API, vs Float64 CPU reference)

Measured on this host (A100, `test/test_mie_nodes.jl`), KA-CPU backend, standard
log-normal set (nquad=300) and a wide-range (0.005-6 μm) TOMAS-like set:

| Metric | standard set | wide TOMAS set |
|---|---|---|
| `k` (extinction) relerr | 8.1e-8 | 7.3e-9 |
| `ω̃` (SSA) relerr | 7.8e-8 | 1.9e-8 |
| Greek coefficients, max abs (α/β/δ/ζ worst) | 5.0e-3 | 4.3e-3 |
| P11 (reconstructed phase function), max relerr | 0.89% | 0.19% |
| P12/P11, max abs | 1.6e-4 | 2.5e-4 |
| P33/P11, max abs | **1.4%** | **3.8%** |
| P34/P11, max abs | 3.1e-4 | 2.8e-4 |

`k`/`ω̃` land close to `DSEmulated`'s established floor (both are dominated by
Neumaier-compensated summation, not the Dₙ recursion). Greek coefficients and
the P33/P11 polarized ratio are the clearest signal of `NativeFloat32`'s
reduced accuracy vs `DSEmulated` (no double-single Dₙ emulation, and — unlike
the log-normal `MieModel` GPU path below — the reduction itself is never
Float64-widened, by design, to keep zero Float64 device arrays). P12/P11 and
P34/P11 stay small throughout. See `test/test_mie_nodes.jl`'s "NativeFloat32:
accuracy..." and "NativeFloat32: polarized phase-matrix validation gate"
testsets for the exact reproduction recipe.

On the log-normal `MieModel` GPU path (`compute_aerosol_optical_properties_gpu`),
`NativeFloat32` lands at `k`/`ω̃` relerr ~1-2e-7 and Greek coefficients ~7e-5
abs — much closer to `NativeFloat64` than to the caller-node path's numbers
above — because that path's host-side size-distribution reduction is ALWAYS
Float64-widened regardless of precision policy (see `test/local/gpu/test_mie_gpu.jl`'s
"NativeFloat32 (log-normal MieModel GPU path)" testsets).

## What is GPU-safe today

- **Forward elastic RT** (elemental → doubling → interaction, `noRS`) —
  fully GPU-friendly on CUDA; this is the hot path the batched design
  targets.
- **Linearized RT (Jacobians)** — GPU-friendly on CUDA
  (`test/local/gpu/test_jacobians_GPU.jl`); Metal Jacobian workflows have
  not been validated yet, so use `CPU()` or `GPU()` for linearized runs.
- **Raman (RRS/VS)** — runs on GPU (`test/local/gpu/test_forward_raman_gpu.jl`)
  but carries materially higher memory pressure than the elastic path — the
  4-D cross-wavelength coupling arrays scale with `nSpec * n_Raman` — so
  watch VRAM on large spectral batches.
- **Atmosphere/surface split caches** (`rt_run_atmosphere` / `rt_run_surface`,
  see [Fast Re-runs & Batch Processing](batch_processing.md)) keep their
  `AtmosphereRTCache` snapshot device-resident: the cached R/T/J arrays and
  `τ_sum_surf` stay on the architecture the cache was built on, so a
  GPU-built cache replays on GPU without a host round-trip.
- **Mie scattering** — automatic GPU path on CUDA
  (`make_mie_model(...; architecture=GPU())`, NAI2 only; PCW and the
  ForwardDiff AD path fall back to CPU). `MetalGPU()` now also has a
  registered GPU Mie pipeline (`has_gpu_mie(::MetalGPU) = true`, Float32-only,
  `NativeFloat32`); `has_gpu_mie` is a single architecture-level trait shared
  by the log-normal `MieModel` GPU path AND the caller-node Mie API, so both
  route through Metal — and BOTH are `NativeFloat32`-aware end to end (shared
  `_mie_kernel1` Kernel-1 dispatch, shared Metal-only-Float64-device-array
  guard). The caller-node Mie API (`compute_aerosol_optical_properties_nodes`/
  `compute_aerosol_extinction_nodes`) is the primary target of the Metal Mie
  route (device-resident reduction, so its `RA` reduction type must never be
  Float64 under `NativeFloat32`); the log-normal path's reduction is
  host-side regardless of policy, so it was never at risk of an illegal
  Float64 device array in the first place.
- **Gas absorption** — the production pipeline (`model_from_parameters` →
  `rt_run`) computes line-by-line absorption via
  [AtmosphericAbsorption.jl](https://github.com/RemoteSensingTools/AtmosphericAbsorption.jl),
  the sole supported absorption dependency going forward, with its own backend
  dispatch. The internal standalone `Absorption` module is legacy-only and will
  be removed; do not use its GPU Voigt kernel for new work.

## CPU fallback

`CPU()` is always available and is the safe default for debugging scenes
and for any workflow not yet validated on a GPU backend (Metal Jacobians).
Pass `architecture=CPU()` to `make_mie_model` to force CPU Mie regardless of
the RT solver's own architecture. The local GPU test suite
(`test/local/gpu/`, not part of CI) self-skips cleanly when no CUDA device
is present.

## Performance Notes

The main architecture switch is the `radiative_transfer.architecture` field in
the scene configuration, or `params.architecture = GPU()` after loading a
configuration. GPU runs are most useful when enough spectral points, layers, or
viewing geometries are batched to amortize kernel-launch and transfer costs.
Use `CPU()` for small debugging scenes and for workflows that are not yet
GPU-safe.

For Apple Silicon experiments, load Metal.jl and use `MetalGPU()` with
`Float32` scene parameters. The current Metal path focuses on the core batched
matrix multiply and inverse operations; fused kernels can be added after Mac
validation. The portable inverse kernel uses Metal threadgroup memory and is
intended for modest stream/Stokes dimensions; larger matrices fail early with a
clear local-memory error instead of a driver launch failure. With the current
32 KiB guard, Float32 matrices with `N = Nquad * nStokes >= 64` are rejected.
Metal Jacobian workflows have not been validated yet, so use CPU or CUDA for
linearized runs.
