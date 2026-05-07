# vSmartMOM nSpec=1 multi-layer Stokes_I numerical issue

**Status**: Surfaced 2026-05-05 by `test/vlidort_baseline/cases/case_B_solar_tester.jl`.
**PARTIALLY** patched 2026-05-06 in `src/CoreRT/tools/cpu_batched.jl` per Codex
investigation — singleton-batch defers to `NNlib.batched_mul`. That fix
unblocks **Case A** (single-layer Stokes_IQUV, nSpec=1, GaussQuadHemisphere)
and **Case C** (multi-layer Stokes_IQU, nSpec=1, then tested under RadauQuad), with a dramatic
speedup (Case A 12 min → 2 min). But **Case B** (multi-layer **Stokes_I**,
nSpec=1) still produces ~50× wrong intensity. Bug is therefore narrower than
"all singleton-batch CPU multiplies": specific to scalar mode + multi-layer
+ nSpec=1. Branch `SS-exact`.

**Definition-audit update (2026-05-06):** the solar-tester VLIDORT reference
uses `NSTREAMS=8` Gauss-Legendre half-space streams, `NMOMENTS=15`
(`m=0..15`), Task 1 has `DO_DELTAM_SCALING = .false.`, and the saved
reference geometries have relative azimuths `[0, 90, 180]` degrees. The
baseline suite now has a preflight check for these invariants and no longer
runs Radau variants as VLIDORT-equivalent tests. Full Case B RT still needs a
fresh timed run after the definition fixes; the last MCP attempt exceeded the
120 s tool ceiling.

## Symptom

When all of the following hold simultaneously:

- `polarization_type = Stokes_I()`
- 23 atmospheric layers (multi-layer)
- Per-layer optical properties injected post-`model_from_parameters` via
  `model.τ_rayl[1][:, n] .= …`, `model.τ_abs[1][:, n] .= …`,
  `model.τ_aer[1][1, n] = …`, plus `model.aerosol_optics[1][1] = …`
- `LambertianSurfaceScalar(0.05)`
- `radiative_transfer.spec_bands` resolves to `nSpec = 1`
  (e.g. `"[18867.92]"`)

`rt_run(model, i_band=1)[1]` returns intensity ~50× smaller than the VLIDORT
2.8.3 reference for the same atmosphere + geometry.

Switching only to `nSpec = 2` (e.g. `"[18867.92 18867.93]"`) — same physics,
just a duplicated spectral point — restores agreement to ~5e-4.

Concrete numbers at `sza=35°`, `raz=0°`, `vza=10°` (VLIDORT
`results_solar_tester.all` Task 1, geom 1, level 1, dir 1):

| nSpec | modeled `R[1,1,1]` | truth | rel err |
| ----- | ------------------ | ----- | ------- |
| 1     | 1.309e-3           | 6.460e-2 | 0.98 |
| 2     | 6.457e-2           | 6.460e-2 | 4.6e-4 |

Same per-layer τ_rayl / τ_abs / τ_aer values are confirmed in both runs (the
injection prints them — see `inject_solar_tester_optics!` in the case file).

## Why it doesn't trip Case A (Siewert)

Case A is single-layer + `Stokes_IQUV()` and passes at ~1e-6 with `nSpec=1`.
The bug is specific to the conjunction above; either the multi-layer adding
loop, the Stokes_I scalar code path, or the per-layer injection produces a
different result for `nSpec=1` than for `nSpec=2`.

## Workaround in current Case B

`test/vlidort_baseline/configs/solar_tester.yaml` carries
`spec_bands: ["[18867.92 18867.93]"]` (two-point spectral band). The case
file compares `R[:, 1, 1]` (first spec point); both points produce identical
intensity by construction. This is documented in the case file header and
this dev note.

## Triage starting points

- Compare per-layer added-layer matrices `(r⁻⁺, t⁺⁺, j₀⁺)` between nSpec=1
  and nSpec=2 builds — they should differ only by trivial broadcasting.
- Inspect whether the `Stokes_I` scalar code path collapses a `(1,1,nSpec)`
  array to a scalar incorrectly when `nSpec=1`.
- Check whether `compEffectiveLayerProperties.jl` or the doubling kernel
  uses `nSpec` to size temporaries that misbehave at `nSpec=1`.
- Try the same Case B with `Stokes_IQU` to confirm/reject the
  Stokes_I-specific hypothesis.

## Root cause (Codex investigation, 2026-05-06)

The CPU override of `NNlib.batched_mul` for dense `Array{T,3}` BLAS-floats in
`src/CoreRT/tools/cpu_batched.jl:64` parallelizes over the batch (third) axis
with `Threads.@threads for k in 1:size(C,3)` and a per-slice `mul!`. At
`size(C,3) == 1` (singleton spectral batch), this path produces wrong values
— effectively a ~1/50 magnitude error on TOA-up I — while the same kernel
with `size(C,3) >= 2` is correct.

Codex traced dispatch: in the multi-layer + elastic + Stokes_I + CPU path,
all the per-Fourier-moment kernel composition (`r⁻⁺ ⊠ t⁺⁺` etc. in
`interaction.jl:176-211` and `rt_helpers.jl:102-163`) funnels through
`batched_mul`. Single-layer Case A bypasses this because there's no
multi-layer accumulation; `nSpec=2` bypasses by avoiding the singleton-batch
specialization.

## Patch

Defer singleton batches to NNlib's generic path:

```julia
function batched_mul(A::Array{T,3}, B::Array{T,3}) where {T<:LinearAlgebra.BLAS.BlasFloat}
    @assert size(A,3) == size(B,3) ...
    @assert size(A,2) == size(B,1) ...
    if size(A, 3) == 1
        return invoke(NNlib.batched_mul,
                      Tuple{AbstractArray{T,3}, AbstractArray{T,3}},
                      A, B)
    end
    C = Array{T,3}(undef, size(A,1), size(B,2), size(A,3))
    Threads.@threads for k in 1:size(C,3)
        @views mul!(C[:,:,k], A[:,:,k], B[:,:,k])
    end
    return C
end
```

The `invoke(NNlib.batched_mul, Tuple{AbstractArray{T,3}, AbstractArray{T,3}}, A, B)`
form bypasses our specialized method and dispatches to NNlib's own
implementation, which has been independently tested and handles singleton
batches correctly.

Historical verification before the definition-audit correction (combined run
with quadrature-swapped YAMLs):
- Case A (Siewert single-layer IQUV, GaussQuadHemisphere, nSpec=1): 8/8
  pass in 1m53s (down from 12m09s). ✓
- Case B (solar_tester scalar Stokes_I, RadauQuad, nSpec=1): 0/2 pass —
  same 0.98 rel-err pattern as the original bug. The cpu_batched fix did
  NOT resolve Case B's path. Bug is elsewhere.
- Case C (solar_tester vector Stokes_IQU, RadauQuad, nSpec=1): 2/2 pass.
  Case C did not previously have a documented nSpec=1 issue; combined
  with Case A's improvement, suggests the cpu_batched fix is real
  hygiene that just doesn't reach Case B's specific scalar code path.

## Remaining triage for Case B

The Stokes_I (n=1) path uses `pol_type.n*Nquad`-sized matrices that are
**smaller** than IQU/IQUV (12×12 vs 36×36 vs 48×48 here). The bug must
live in code that is exercised by multi-layer doubling/interaction
specifically when Stokes_I is in use:

- `compute_Z_moments` / `construct_B_matrix` returns a *scalar* for
  `Stokes_I` (`construct_B_matrix(::Stokes_I, ...) = β[l]`), unlike the
  3×3/4×4 `SMatrix` returned for the polarized cases. Could a downstream
  reshape / batched op mishandle that scalar-vs-matrix distinction at
  nSpec=1 specifically?
- `batch_inv!` / `batch_solve!` in `cpu_batched.jl` — same `Threads.@threads
  for i = 1:size(A,3)` pattern as the (now-fixed) `batched_mul`. May
  warrant the same singleton-batch guard for completeness.
- The Stokes_I doubling/interaction kernel may have a scalar-vs-matrix
  type mismatch only at nSpec=1.
