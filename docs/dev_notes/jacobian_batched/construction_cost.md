# Model construction: spectroscopy reuse and Mie derivatives

Before spectroscopy caching, the full-rebuild benchmark paid mostly for
repeated HITRAN text parsing. Downloaded files were already present, but each constructor called `load_lines`
again for CO₂, CH₄ and H₂O. Neither constructor supplies a wavelength-window
filter to that loader. In AtmosphericAbsorption 0.1.2, `load_lines` calls
`parse_par`, converts the text columns, sorts the line centers and constructs
new line-database arrays. That dependency loader does not cache its results;
the constructors now supply the reuse boundary described below.

The constructors are in `src/CoreRT/tools/model_from_parameters.jl` and
`lin_model_from_parameters.jl`. The loaded dependency source and sampled call
stacks are recorded in [the profile log](evidence/source_adding/construction-profile.log).

## Baseline warmed profile (before caching)

Same A100/Float64 IQU scene as the source-adding report: 20 layers, three
weighted streams, one aerosol, CO₂/CH₄/H₂O, and the first 1,000-point chunk
of the 10,000-point 6150–6250 cm⁻¹ grid. Three warmed constructions per mode;
compilation and downloads are excluded. Total timings synchronize CUDA.
Section times below use the existing constructor timers and are approximate.

| Operation | Forward constructor | Linearized constructor |
|---|---:|---:|
| Complete construction, median | 2.687 s | 2.610 s |
| Read/parse all three HITRAN files | ≈2.46 s | ≈2.30 s |
| Evaluate gas absorption across the profile | ≈0.21 s | ≈0.23 s |
| Two timed Mie-node evaluations | ≈0.002 s | ≈0.062 s |

Thus line loading accounts for roughly 90% of construction. It allocates about
864 MiB on the host per build, out of approximately 886 MiB forward and 955 MiB
linearized total. Garbage collection also costs about 0.48 s forward / 0.28 s
linearized in the median samples; that time is already included above, not an
additional cost. The untimed remainder includes reference-extinction and other
setup work; the two Mie-node timers should not be interpreted as every operation
related to aerosols.

The sampled stacks confirm that much of the time is in `parse_par` and numeric
text parsing. This is repeated preparation on already downloaded files, not
network-download latency or an inherent multi-second Mie cost for this scene.

## Why the 10,000-point total was so large

The benchmark deliberately rebuilt each of ten 1,000-point chunks independently.
It therefore paid the same file-parsing cost ten times per mode. The previously
reported roughly 25–26 seconds of construction is that sum, not a measurement
of one efficient 10,000-point construction.

This also explains the small 1.15× full-rebuild Jacobian/forward ratio: both
timers contain a large shared setup cost. The prepared-model RT ratios remain
the relevant evidence for the propagation improvements. Removing repeated
setup should reduce absolute runtimes while potentially increasing this
particular full-rebuild ratio.

## Implemented spectroscopy reuse

`src/CoreRT/tools/hitran_cache.jl` owns a process-local parsed line cache keyed
by resolved file path and floating-point type. Forward, linearized, Raman
construction and `BatchContext` initialization all use it. The lock covers the
initial load, preventing duplicate parsing by concurrent constructors. Each
scene still creates its own cheap absorption-model wrapper with its own
backend, line shape and broadening settings. No spectral clipping is added.

There is no automatic eviction or reload. Each database stays resident until
`clear_spectroscopy_cache!()` or process exit. This follows the requested
persistent reuse policy; memory grows with the distinct file/type pairs used.
Changing editions selects a new path. Replacing a downloaded file in place
requires explicitly clearing the cache before constructing new models.
Existing models retain their original data; clearing does not mutate them.

LUTs remain caller-owned: load them once into parameters and reuse those objects
across constructions. `AUDIT_LUT_DIR` in the benchmark exercises this ownership
with existing legacy HITRAN `.jld2` tables and reports one-time load costs
separately. The LUTs and current line-by-line setup may differ in isotopes,
resolution and broadening; the benchmark tests forward/Jacobian parity within
each setup, not numerical equivalence between the two spectroscopy sources.

## Mie derivatives

The same aerosol node, evaluated explicitly on CPU for both modes (Float64,
λ = 1.626016 μm, 100 radius nodes, five warmed samples), initially cost
0.892 ms forward versus 28.943 ms including four native Mie derivatives.
A sampled profile pointed to scalar phase-amplitude derivative evaluations.
`Aerosol.size_distribution` has an abstract distribution field; unresolved
quadrature types propagated into the numerical loops and boxed arithmetic.

The linearized NAI-2 entry point now prepares the radius quadrature and
normalized weights, then calls `_nai2_bulk_optics_lin` with concrete arrays.
This function boundary lets Julia specialize the existing numerical body.
No equations, quadrature rules or summation order were changed.

| CPU Mie | Before | After |
|---|---:|---:|
| Forward | 0.892 ms | 0.860 ms |
| Forward + native derivatives | 28.943 ms | 2.009 ms |
| Derivative-path host allocation | 27.74 MB | 2.97 MB |

That is a 14.4× derivative-path speedup, with a remaining 2.34× cost relative
to CPU forward Mie for this node. These are Mie evaluations, not complete
construction or RT timings. The CUDA constructor still uses its existing
forward Mie backend; these CPU numbers isolate the derivative implementation
from backend differences. Larger particles and other distributions have not
been assigned this speedup.

Reproduce from `test/`:

```sh
CUDA_VISIBLE_DEVICES=0 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  AUDIT_BACKEND=cuda AUDIT_NSPEC=1000 AUDIT_CHUNKS=10 AUDIT_LAYERS=20 \
  AUDIT_GASES=CO2,CH4,H2O AUDIT_NU_MIN=6150 AUDIT_NU_MAX=6250 \
  julia --project=. ../docs/dev_notes/jacobian_batched/construction_profile.jl
```

## Cached-LUT end-to-end result

A100, Float64 IQU, 20 layers, three weighted streams, one aerosol, 69 Jacobian
columns (60 nonzero gas columns), embedded solar and a Lambertian surface.
The 10,000-point 6150–6250 cm⁻¹ grid is processed as ten 1,000-point chunks.
Each timed chunk independently rebuilds parameters and the model; three warmed
trials are summed across chunks and their median is reported. This follows the
same endpoint/diagnostic boundary as the earlier embedded-solar benchmark:
forward uses public `rt_run` including hemispheric diagnostics; Jacobians return
TOA/BOA endpoints with `jacobian_basis=:local, jacobian_adding=:source`.

| Cached legacy LUT setup | Forward | Forward + 69 Jacobians |
|---|---:|---:|
| Full construction + RT | 1.834 s | 5.263 s |
| Construction subtotal | 0.169 s | 0.266 s |
| RT subtotal | 1.661 s | 4.989 s |

The full-rebuild ratio is **2.870×**. Subtotals are separately medianed, so they
need not add exactly to the full median. The maximum forward-radiance mismatch
between independent forward and linearized constructors/solves is 3.82e-17.
This is evidence for the stated scene, not a universal ratio or the separate
<2× target for a solve from supplied core optical properties.

The caller-owned tables were `/home/sanghavi/data/HITRAN_LUTs/{CO2,CH4,H2O}.jld2`.
They have 3,285,715 spectral nodes spaced by 0.01 cm⁻¹, covering approximately
2857–35714 cm⁻¹, pressure 0.01–1080.01 hPa and temperature 180–360 K. The input
profile is inside these ranges. Stored isotope IDs are CO₂: 1, CH₄/H₂O: -1
(all). The legacy LUT evaluator has no varying H₂O broadener axis; the current
line-by-line setup can differentiate moist self-broadening. A spectroscopy
accuracy comparison would need matched table-generation settings.

One-time table loading took 10.38 s in this run (including initial loader
compilation) and allocated about 5.98 GB on the host. The three loaded tables
stay resident across every chunk, forward and Jacobian construction. These
startup costs are excluded from the warmed timings above; they must be
amortized over repeated scenes.

Evidence: [LUT measurements](evidence/source_adding/lut-end-to-end.toml),
[load metadata and parity log](evidence/source_adding/lut-end-to-end.log),
[Mie before](evidence/source_adding/mie-before.log),
[Mie after](evidence/source_adding/mie-after.log),
[bitwise baseline parity](evidence/source_adding/mie-parity.log), and
[targeted CPU regressions](evidence/source_adding/regressions.log).

Reproduce the LUT benchmark from `test/`:

```sh
CUDA_VISIBLE_DEVICES=0 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 \
  AUDIT_BACKEND=cuda AUDIT_NSPEC=1000 AUDIT_CHUNKS=10 AUDIT_LAYERS=20 \
  AUDIT_GASES=CO2,CH4,H2O AUDIT_NU_MIN=6150 AUDIT_NU_MAX=6250 \
  AUDIT_END_TO_END=true AUDIT_LUT_DIR=/home/sanghavi/data/HITRAN_LUTs \
  julia --project=. ../docs/dev_notes/jacobian_batched/local_basis_benchmark.jl
```

For the CPU Mie profile, run `mie_derivative_profile.jl` with
`CUDA_VISIBLE_DEVICES=-1 AUDIT_NSPEC=1000 AUDIT_CHUNKS=10 AUDIT_NU_MIN=6150
AUDIT_NU_MAX=6250 JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1`.

## Cached direct-HITRAN end-to-end result

The same scene, timing protocol and source-adding path using current direct
HITRAN evaluation, with parsed line data retained across every chunk:

| Cached direct HITRAN | Forward | Forward + 69 Jacobians |
|---|---:|---:|
| Full construction + RT | 3.174 s | 6.473 s |
| Construction subtotal | 1.564 s | 1.646 s |
| RT subtotal | 1.606 s | 4.832 s |

The complete ratio is **2.039×**, and independent forward-radiance parity is
4.16e-17. All 60 gas columns are nonzero in each chunk. Loading/parsing the
three line databases happens in initial fixture preparation; it is excluded
from these warmed repeated-use timings, just as one-time LUT loading is.
The original direct-HITRAN reconstruction took 28.143 s forward / 32.456 s
including Jacobians. The new totals are about 8.9× / 5.0× faster, respectively,
for this scene. Construction subtotals fall from 26.156 / 24.829 s to
1.564 / 1.646 s. The larger new Jacobian/forward ratio reflects removal of the
shared parser cost; both absolute runtimes are much shorter.

For reproduction, use the LUT command above without `AUDIT_LUT_DIR`. Evidence:
[cached-HITRAN measurements](evidence/source_adding/hitran-cached-end-to-end.toml)
and [parity log](evidence/source_adding/hitran-cached-end-to-end.log).

## Validation

The numerical Mie body is bitwise identical to the `8bc60941` baseline in eight
cases spanning Float32/Float64, two wavelengths and absorbing/nonabsorbing
particles. The targeted CPU regression run passed 161 checks: initial parsed
cache tests, aerosol reference finite differences, source-adding tests, and
H₂O self-broadening. The extended cache test passed 13/13, including concurrent
first-load reuse, type separation, independent fresh-data parity, BatchContext
sharing and clearing ownership without mutating existing users. The latter
includes the 11 initial cache checks rather than 13 additional distinct tests.
Both 10,000-point CUDA benchmarks passed independent forward/LinMode radiance
parity. The strict local Documenter/Vitepress build passed. Logs:
[extended cache](evidence/source_adding/construction-cache-final.log),
[docs](evidence/source_adding/construction-docs.log).
