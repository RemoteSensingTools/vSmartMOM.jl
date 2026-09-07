# Overall speed change through e5a4d38e

The original 10,000-point absorption-free fixture was rerun on the current
implementation to establish a direct historical timing comparison. The YAML
is unchanged from integration commit `b524a0d8`; the benchmark applies the same
radius cap (3 μm), median radius (0.15 μm), geometric width (1.4), disabled Greek
cutoff, and Float64 scalar/IQU choices as the original driver.

## Original fixture: first recorded reference versus current implementation

A100 PCIe 40 GB, Julia 1.12.6, four Julia threads, one BLAS thread, three warmed
synchronized samples. Five layers, one aerosol, Lambertian surface, embedded
solar, principal-plane views, 14 native columns (including five zero gas slots),
and all 10,000 wavelengths in one spectral batch. Mie/model construction is
outside timing; optical mixing, adding–doubling, surface and endpoints are in.
Current Jacobians use `jacobian_basis=:local, jacobian_adding=:source`.

| Complete forward + Jacobian solve | Initial reference | Current | Speedup |
|---|---:|---:|---:|
| Scalar I | 16.716 s | 0.434 s | **38.50×** |
| Polarized IQU | 38.268 s | 2.343 s | **16.33×** |

The matching forward-only historical/current measurements are 0.15486 / 0.10747 s
for scalar (1.44× faster), and 0.80079 / 0.28012 s for IQU (2.86× faster).
These historical forward medians come from the same initial log files as the
Jacobian references; subsequent reports also contain other forward samples.
Current Jacobian/forward ratios on this fixture are 4.04× scalar and 8.37× IQU.
The <3× end-to-end result on the larger real-gas benchmark therefore must not be
presented as a universal ratio for every batch, optical scene and timing boundary.

Current forward/Jacobian endpoint radiances agree within 5.90e-17. This compares
the two current paths, not numerical parity with historical code: scientific
normalization and tangent corrections were also made during development.
Historical and current timings were taken on a shared machine at different
times; they are workload-specific medians, not dedicated-machine guarantees.

## Later realistic benchmark and construction improvements

The 20-layer IQU case has 69 columns, including 60 nonzero CO₂/CH₄/H₂O profile
columns, and ten independently prepared 1,000-point chunks over 6150–6250 cm⁻¹.
Its earliest recorded prepared-model physical propagation took 61.044 s;
current source adding takes 4.892 s: **12.48× faster**. See
[local basis](local_basis.md) and [source adding](source_adding.md).

For full repeated construction plus solving with direct HITRAN, the recorded
32.456 s combined runtime fell to 6.473 s: **5.01× faster** during the subsequent
construction optimization. Forward-only fell from 28.143 to 3.174 s. Those
32.456/28.143 s measurements already had optimized source adding; they are not
measurements of the original integration implementation. See
[construction report](construction_cost.md), including the cached-LUT result
(5.263 s combined / 1.834 s forward) and Mie derivative speedup (14.4×).

There is no first-integration full-construction benchmark for the realistic
scene. Consequently, no single first-to-current end-to-end multiplier is
established, and the solver and construction speedups must not be multiplied.

## Evidence and reproduction

Initial timings: [scalar log](evidence/benchmark-cuda-n10000.log),
[IQU log](evidence/benchmark-cuda-iqu-n10000.log).
Current timings: [measurements](evidence/source_adding/original-fixture-current.toml),
[run log](evidence/source_adding/original-fixture-current.log).

From `test/`:

```sh
CUDA_VISIBLE_DEVICES=0 AUDIT_BACKEND=cuda AUDIT_NSPEC=10000 \
  JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  julia --project=. ../docs/dev_notes/jacobian_batched/historical_fixture_benchmark.jl
```
