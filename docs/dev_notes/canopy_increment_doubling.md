# Canopy transmission increments

The directional canopy elemental operator approaches the identity as its
thickness tends to zero. Rounding `T = I + O(δ)` in Float32 and then repeatedly
squaring can amplify lost diagonal increments into apparent absorption or gain.
Promoting only the doubling stage cannot recover information lost in the seed.
Coarser seeds avoid that cancellation but increase finite single-scattering
initialization error. The canopy path therefore carries the increment directly.

## Algebra

All operators below use the existing D-transformed basis during doubling.
Let `Δ = T-I`, `S = (I-RR)⁻¹`, `H = S-I = S RR` and
`C = T S-I = Δ+H+ΔH`. For identical slabs, standard adding-doubling gives

```
Q = R + RΔ
R_new = R + Q + C Q
Δ_new = C + Δ + C Δ
```

No commutativity is assumed. `H` is obtained by solving `(I-RR)H=RR`.
For the elemental diagonal, form `expm1(-x) + exp(-x)*scattering`
directly, where `x = δτ G/μ`; never subtract one from an already rounded T.
Off-diagonal transmission is already an increment.

The solar attenuation increment is `β=expm1(-δτ G₀/μ₀)`. Its doubling
update is `β_new=β(2+β)`. If `j_d` and `j_u` are the current slab sources,
the sources of the lower copy are `j_d1=j_d+β j_d`, `j_u1=j_u+β j_u`.
Using old R and old sources throughout each right-hand side:

```
v_u = R j_d + j_u1
j_u_new = j_u + v_u + C v_u
v_d = R j_u1 + j_d
j_d_new = j_d1 + v_d + C v_d
```

Only after all doublings is `T=I+Δ` reconstructed. Existing D-matrix
restoration then supplies the reverse operators and upward solar source.
There is no phase, flux, absorption or operator normalization correction.

## Scope and tests

`test_canopy_increment.jl` compares this algebra with ordinary doubling at
matched Float64 seed thickness and compares Float32 with Float64 at thin
seeds. It spans Fourier orders 0–2, Stokes I/IQU/IQUV, diffuse and combined
diffuse/specular optics, one and three spectral points, and multiple LAIs.
The optional CUDA run requires `VSMART_CANOPY_INCREMENT_GPU=true` and disables
scalar indexing. `test_canopy_precision_surface.jl` checks the actual canopy
wrapper against a conservative photon budget, including multiple physical
layers over a perfectly reflecting soil. Settings and tolerances are in TOML.

Local validation on Julia 1.12.6 and an NVIDIA L40S: 3,744 operator/phase checks
on CPU, 16 vertical-solar-beam checks, and 3,552 CUDA checks passed, including
zero-thickness layers and finite entries throughout the cached phase matrices. The production
surface tests passed 48 CPU and 8 CUDA checks; the existing canopy suite passed
135 checks, including its soil-albedo derivative regression. For the dense
conservative fixture, direct-budget error at 64 streams was approximately
1.2e-6 in Float32 and 3.2e-6 in Float64. The 16-stream angular discretization
retained approximately 1.8e-4 error in both precisions; the strict 1e-4 budget
target is checked at 32 and 64 streams. No result is rescaled to hit that target.

The combined specular Float32 cases used the separate CanopyOptics fix in
[PR #64](https://github.com/RemoteSensingTools/CanopyOptics.jl/pull/64).
The complete CPU package suite also passed 25,512 checks, with 16 conditional
skips reported by Julia in the `Broken` column (CUDA was masked in that run).
That full run used the azimuth fix from PR #64; after the additional grazing
limit fix, the focused counts above were rerun against the updated dependency.
The full CanopyOptics suite passed 86,099 checks with both fixes. Upstream CI
and release dependencies remain separate gates.

These tests establish numerical agreement over their stated cases, not
physical correctness of every specular parameterization or infinite-resolution
conservation. Angular resolution, leaf-angle quadrature and elemental thickness
must be converged separately. Metal and arbitrary dual-number differentiation
have not been validated for this path. The existing canopy soil-albedo
linearization uses finite differences of the forward wrapper; this change does
not add leaf or canopy analytic Jacobians.

This implementation carries the solar source only and explicitly rejects
additional source slots. A distributed fluorescence source cannot use the
solar attenuation recurrence: its two sublayers emit according to their own
source distribution. It needs a separate source law and elemental integral.
The reflection/transmission algebra is reusable when those sources are added.

## CPU solve and measured costs

The square CPU geometric-series system now uses pivoted LU, with its factorization
in the existing scratch matrix and its right-hand side copied into the output.
This specialization is local to the canopy increment path. Other backends and
non-BLAS element types retain `batch_solve!`; the general atmospheric solve API
is unchanged. After this change, all 3,744 operator/phase, 16 vertical-beam and
48 production-surface CPU checks passed. The complete package-suite result above
predates this CPU specialization; it is not a full-suite claim for the new solve.

`test/bench_canopy_increment.jl` measures warmed doubling calls, including a
completion synchronization on CUDA. Phase-cache construction, elemental seeding,
host transfers, compilation and explicit pre-sample garbage collection are
outside the timed region. The fixture sets two warm-ups and five measured calls;
the CSV records minimum and median times and minimum host allocations. Host
allocations are not GPU device-memory measurements. Run from `test/`:

```sh
julia --project=.. bench_canopy_increment.jl
VSMART_CANOPY_BENCH_GPU=true julia --project=.. bench_canopy_increment.jl
```

On Julia 1.12.6, the Sapphire Rapids host, one Julia thread and one BLAS thread,
the 64-spectral-point CPU batches gave these median doubling times:

| Precision | Streams | Stokes | Previous QR (ms) | New LU (ms) | LU host allocation (MB) |
|---|---:|---:|---:|---:|---:|
| Float32 | 8 | 1 | 12.77 | 4.63 | 0.248 |
| Float32 | 8 | 4 | 94.90 | 60.19 | 0.642 |
| Float32 | 32 | 1 | 86.92 | 47.80 | 0.633 |
| Float32 | 32 | 4 | 1112.08 | 1099.58 | 2.183 |
| Float64 | 8 | 1 | 14.08 | 4.49 | 0.249 |
| Float64 | 8 | 4 | 115.52 | 72.49 | 0.643 |
| Float64 | 32 | 1 | 102.57 | 67.67 | 0.634 |
| Float64 | 32 | 4 | 2004.48 | 1955.21 | 2.186 |

For the largest Float64 case, host allocations fell from 406.04 MB to 2.19 MB,
while elapsed time changed little: matrix products remain substantial work.
The ordinary transmission formulation remains faster in most of these cases,
but loses thin-seed Float32 accuracy as documented above. Its timings are a
cost baseline, not evidence that its inaccurate solution is interchangeable.

Raw [QR](data/canopy_increment_cpu_qr.csv) and
[LU](data/canopy_increment_cpu_lu.csv) CSVs include singleton batches and the
ordinary baseline. Values are shared-host snapshots, not ESM throughput or a
claim about ideal parallel scaling. End-to-end and device-allocation measurements
remain necessary before selecting an ESM backend. GPU timing also needs to
account for the ordinary path's existing fused solves and updates; the initial
increment implementation does not yet reuse those optimizations.

The [initial L40S GPU timings](data/canopy_increment_gpu_initial.csv) cover the
same 32 method/precision/batch cases with scalar indexing disabled. For example,
the Float32 8-stream, four-Stokes, 64-point case took 5.54 ms in increment form
versus 1.93 ms in ordinary form. The precision fix is not presently a speed
optimization. A next implementation step is to reuse the existing fused
geometric-series solve for `H = R² (I−R²)⁻¹`; this is algebraically equivalent
to `(I−R²)⁻¹R²` because the two factors commute. It must retain the small `R²`
right-hand side, rather than subtracting identity from a rounded inverse, and
must be checked with the same precision and conservation tests before adoption.
