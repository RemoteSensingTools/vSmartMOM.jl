# Local optical directions and the remaining adding cost

Development branch: `perf/jacobian-batched-propagation`. Local optical factoring
is committed as `fe6b865e`, with the medium CUDA gate in `ce3d38a4` and
the regular-limit correction in `4df62396`. CUDA product selection and
solve-owned pointer batches are committed as `e06a1a54`.
This extends the [IQU investigation](iqu_followup.md) and the
[paper review](paper_review.md). The target is the **complete forward plus
Jacobian solve from supplied core optics**, below twice forward cost. The
measurements below do **not** meet that target.

## Production ordering

The default `jacobian_basis=:auto` uses a local basis when its dimension is
smaller than the requested atmospheric layout. `:physical` retains the direct
retrieval-column path for comparison; `:local` forces factoring. Factoring
currently assembles one band per solve. A multi-band combined call retains the
physical path, while ordinary separate band calls can each use the local basis.
Selected layouts keep their existing column order and coordinate transforms.

For aerosol mode i, let Sᵢ=τᵢϖᵢ after truncation, S the total scattering depth
including Rayleigh, and αᵢ=Sᵢ/S. Using the fixed Rayleigh phase matrix Zᵣ,

```math
Z = Z_r + \sum_i \alpha_i (Z_i-Z_r),\qquad
\dot Z = \sum_i [\dot\alpha_i(Z_i-Z_r)+\alpha_i\dot Z_i],
```

```math
\dot\alpha_i=(\dot S_i-\alpha_i\dot S)/S,\qquad
\dot\varpi=(\dot S-\varpi\dot\tau)/\tau.
```

These are derivatives of S2014 (C.22)–(C.24), with the chain rule of
(C.25)–(C.26); see the [verified references](../theory_references.md).
The choice of a fixed Rayleigh reference is an implementation derivation,
not a claim about the historical Fortran ordering.

There are at most **2 + 5 N_aerosol** directions: τ, ϖ, and one phase difference
plus four truncated microphysical phase tangents per aerosol. An active
retrieval layout retains only its selected microphysical directions, giving
**2 + N_aerosol + N_selected_microphysics**. Three species with fixed
microphysics therefore need five directions instead of seventeen. Selection
is compiled before workspace allocation and applies to both matrix and source
adding; `:auto` compares this reduced count with the atmospheric layout.
When a species has no selected microphysics, only its forward phase is evaluated.
All species' mixture directions remain: changing one scattering contribution
also changes the normalized weights of the others. Selection is structural,
so a zero aerosol loading does not discard its potentially nonzero derivative.
An unselected native layout still carries all four microphysical directions
per aerosol, independent of the number of gas/profile columns.
An independently retrieved phase shape or truncation parameter would require
its complete phase direction and coefficient map; it is not a free extra
column. The generic `LocalOpticalJacobian` boundary accepts supplied bases.
The phase basis is shared across all layers for a Fourier order. The scalar
coefficient tensor C[wavelength,basis,parameter] is built once before the
Fourier loop. Gas columns have exactly zero phase coefficients.

The exact finite-thickness elemental map and every doubling operate on this
basis. Complete doubled matrix tangents are then contracted with C. This is
valid because the tangent map is linear in its input perturbation; an
**elementwise** phase partial contracted after doubling would lose angular
coupling and is not equivalent. Adding currently receives retrieval columns.

During local doubling the above-layer beam attenuation E=exp(-τ_above/μ₀)
is fixed. After contraction, the physical source tangent receives
`-J * dτ_above / μ₀` exactly once. Local optical-depth attenuation is already
included in the doubled basis tangent. Both the scalar truncation-factor
terms and the derivative of the normalized truncated Greek coefficients are
retained; no fᵗ chain term is omitted.

## Supporting changes

- Phase-node interpolation now runs on the target backend, with the same
  two-node linear / three-node natural-cubic map as the host reference.
  Forward matrices and microphysical tangents use the same linear map.
- Elastic matrix/source buffers are allocated and zeroed on their backend.
- The retrieval-sized added-layer carrier omits doubling scratch and unused
  elemental core partials when doubling is performed in the local workspace.
- `elemental_lin.jl` now contains the driver and symmetry operations;
  production fused evaluations and reference partial evaluations have their
  own files. Optical caches, phase evaluation, and local-basis construction
  likewise have separate, documented responsibilities.

## Measured complete solves

A100 PCIe 40 GB, Float64, four Julia threads, one BLAS thread. Timings are
three warmed samples with explicit CUDA synchronization. Upstream Mie and
spectroscopy are prepared before timing; optical mixing, all Fourier orders,
adding–doubling, the Lambertian surface and endpoint output are included.
The physical and local modes use the same forward implementation. Forward
optical mixtures retain the same successive weighted-average order as optical
property `+`, including the intermediate `τ·ϖ` products. The local derivative
basis still uses phase differences `Zᵢ-Zᵣ`, but the forward mixture is not
reconstructed by subtracting and adding the Rayleigh reference. In Float32,
these mathematically equivalent reorderings can change an optical input by an
ulp and produce a larger radiance difference after repeated doubling. Cached
spectral mixing weights preserve the forward order with one phase output array
per layer/block; no retrieval-sized phase tangent is introduced. Regression
checks require identical forward τ, ϖ and phase arrays across the two modes.
Rayleigh's conservative albedo is `one(FT)` in the physical Jacobian path;
the former literal `1.0` promoted its Float32 mixtures to Float64. That hidden
promotion was exposed by the O₂/SIF study replay and violated the common
precision contract before elemental propagation.
Adding pure absorption preserves both `Z` and its existing scattering
derivatives. The physical path now appends exact zero gas-phase columns,
matching the local basis, instead of expanding and cancelling a quotient. This
implements `dZ_after = [dZ_before, 0_gas]` directly and avoids both Float32
cancellation residues and large temporary phase tensors.
The five-layer fixtures are absorption-free: their fourteen-column layout
contains five zero gas columns. They exercise aerosol/surface derivatives,
but the larger gas case below is needed to assess active profile columns.

With the expanded small-operator gate, the scalar 10,000-point, five-layer,
fourteen-column check measures 0.767 s for physical propagation, 0.498 s for
the local basis and 0.108 s for forward. The same-run physical reference
without shared-memory products takes 0.848 s. Maximum physical/local
disagreement is 4.45e-16; local/forward is 4.60.

| IQU scene | Physical columns through doubling | Local basis through doubling | Forward |
|---|---:|---:|---:|
| 10,000 points, 5 layers, 14 columns (embedded solar) | 5.728 s | 3.229 s | 0.283 s |
| 10,000 points, 20 layers, 69 columns, real gas absorption | 61.044 s | 24.589 s | 1.739 s |

The first row improves the original integration Jacobian timing of 38.268 s
by about twelvefold, but is still 11.40 times its now faster forward baseline.
The historical scene used principal-plane views; this follow-up uses a 37°
relative azimuth for one view to exercise nonzero U. All fresh physical/local
comparisons use exactly the same scene.
Its physical/local maximum absolute disagreement is 3.93e-16.

The second scene spans 6150–6250 cm⁻¹ with CO₂, CH₄ and humidity derivatives;
all 60 gas-profile columns are nonzero. It uses ten independent 1,000-point
chunks. **Each chunk prepares its own model and Mie interpolation anchors**;
this tests an identical chunking procedure for both paths, not slicing a
single globally prepared 10,000-point model. Total times are medians of the
three summed chunk timings. Maximum disagreement is 2.77e-14. Cumulative
GPU allocation per chunk falls from 191.20 GB to 11.36 GB; these counters
measure allocation traffic, not peak resident memory.

A separate controlled experiment supplies three core directions (τ, ϖ, β₂)
per layer, a black lower boundary, 20 layers and m=0:2. All 10,000 points are
batched together; optical assembly and microphysics are outside timing.

| Retrieval columns | Physical solve | Local solve | Forward | Local / forward |
|---|---:|---:|---:|---:|
| 8 | 8.885 s | 4.958 s | 0.801 s | 6.19 |
| 20 | 20.663 s | 6.978 s | 0.790 s | 8.84 |
| 60 | 59.874 s | 13.480 s | 0.899 s | 15.00 |

This isolates the remaining state-size scaling in adding. It is not an
independent aerosol retrieval benchmark. The CPU companion checks every
column against central finite differences at both endpoints.

## Larger IQU operators

For 512 wavelengths, five layers and fourteen columns, the medium-operator
CUDA path combines wavelength and direction in one launch and stages 16×16
tiles in shared memory. It avoids the per-column fallback used above 32.

| Streams / diffuse size | Previous physical path | Blocked physical path | Blocked local basis | Forward |
|---|---:|---:|---:|---:|
| 8 / 33×33 | 19.414 s | 8.396 s | 5.339 s | 0.428 s |
| 16 / 57×57 | 94.378 s | 49.031 s | 30.514 s | 5.908 s |

All four endpoint fields agree within 3.89e-16 across the tangent paths.
At 16 streams, cumulative allocation is 8.394 TB for the previous path,
352.29 GB for blocked physical propagation and 53.61 GB for the local basis.
These operator sizes include appended viewing/solar nodes.

A subsequent cuBLAS pointer-batch path improves the same complete solves
further. Its device pointer vectors fold wavelength × direction into one
batch, repeating forward spectral pointers without copying matrix values.
One workspace owns the cache and its backing arrays; equivalent active-prefix
views reuse entries, and deep copies rebuild metadata. Source vectors retain
the fused portable kernel. Cache construction is included in these timings.

| Streams / diffuse size | cuBLAS physical path | cuBLAS local basis | Forward | Local / forward |
|---|---:|---:|---:|---:|
| 8 / 33×33 | 5.572 s | 3.580 s | 0.436 s | 8.21 |
| 16 / 57×57 | 26.451 s | 16.639 s | 5.865 s | 2.84 |

The physical path agrees exactly with the blocked reference in these runs;
local-basis disagreement is at most 3.89e-16. Relative to blocked products,
the sixteen-stream local solve is 1.83 times faster. The greater relative
importance of forward matrix work improves amortization, but still does not
establish the <2 target.

The default CUDA gate admits 33–64 operators at 512 or more spectral
points, using cached cuBLAS products when its array-layout requirements hold
and blocked kernels otherwise. CPU behavior and the Metal limit are unchanged. Tests cover ragged
tile edges at sizes 33, 48, 57 and 64, Float32/Float64, mixed 3D/4D operands,
accumulation and the fused product rule. This gate describes tested coverage,
not a measured crossover for every configuration.

## Small external-solar IQU operators

Three weighted streams plus two viewing nodes produce a 15×15 diffuse IQU
operator in external-solar mode. The previous shared-memory gate began at 16,
so this geometry missed operand reuse. Isolated product rules show a 2.47×
improvement at N=15 in Float64 (10,000 wavelengths, three directions).
The gate now begins at six; tests cover N=6, 9, 12, 15, 16, 18 and 32 in
both precisions, including accumulation and active-prefix views.

For the complete absorption-free, five-layer external-solar scene at 10,000
wavelengths, physical propagation falls from **5.877 s to 4.052 s** with
identical outputs. The local basis takes **2.289 s**, versus **0.174 s**
forward. Local/physical disagreement is at most 3.93e-16. This path computes
TOA only; these timings must not be treated as a direct comparison with the
embedded-solar two-endpoint scene. The faster forward path leaves the
local/forward ratio at **13.14**, despite the absolute Jacobian improvement.

## Where the remaining time goes

A separate CUDA trace of the 69-column, 1,000-point scene finds adding to be
the dominant host timer section. The trace includes approximately 779 ms
in tiled product rules, 438 ms in local-to-physical matrix contraction,
277 ms in full tangent additions and 229 ms in full tangent copies. The
fused inverse consumes about 19 ms. These are instrumented device times,
not a partition of the uninstrumented wall-time median.

An exact-zero coefficient guard subsequently removes unnecessary phase-basis
reads during contraction. It does not threshold small physical derivatives.
The larger opportunity is to avoid forming retrieval-sized matrix tangents
for adding at all.

## Layer-response experiment

`layer_response_check.jl` and `layer_response_operators.jl` test a different
association of the same adding equations. Cache forward prefixes Pᵢ and
suffixes Sᵢ, then differentiate `(Pᵢ ⊕ Lᵢ) ⊕ Sᵢ` with local directions in
Lᵢ. All lower **solar** source vectors share exp(-τᵢ/μ₀), so the suffix
source tangent is `-J_suffix * dτᵢ / μ₀`; its operator tangents are zero.
Contract each completed local radiance response with C at the output and
sum over layers.

`source_response_operators.jl` additionally rewrites each local operator
perturbation as equivalent up/down source vectors evaluated on the fixed
incident fields. Two ordinary adding resolvents then propagate these source
vectors; adding itself carries only vector tangents. Its interface-balance
derivation is recorded beside the evaluations. The original matrix-response
association remains available with `AUDIT_RESPONSE_MODE=matrix`. Both
formulations pass the CPU parity and endpoint finite-difference checks.

This makes the number of matrix directions in these adding operations
independent of retrieval size, without an adjoint. Coefficient storage and
output contraction still scale with the state vector. It trades cached
forward subcolumns and more forward compositions for fewer matrix tangents.
It does not claim that Sanghavi's Fortran used this association.

The initial CPU experiment agrees with direct propagation to 1.67e-16,
including non-principal-plane U and parameters that perturb multiple layers;
all 1/18-column endpoint finite differences pass. This remains a development
experiment for solar-only, black-boundary columns. The source-vector formulation also passes all 60-column scalar endpoint
finite differences across 20 layers. On CUDA with IQU, 20 layers and 512
wavelengths, it takes 0.583 s for one retrieval column and 0.581 s for 60;
maximum endpoint-Jacobian disagreement with direct propagation is 3.34e-16.
Its complete cost is still about five times forward in this small batch.
At 10,000 points, the response solve takes 3.934 s for one retrieval column
and 4.017 s for 60, versus forward medians of 0.789 and 0.833 s. This large
run omits the dense comparison workspace after the separate parity/FD checks;
its forward endpoints are still compared in-process. The controlled core
fixture repeats its optical values across wavelength to isolate batch and
state-size costs, unlike the spectrally varying gas benchmark above.

At sixteen streams (57×57 operators), the same 20-layer, 512-point source
response takes 4.516 s for one column and 4.514 s for 60, versus 1.977 s
and 1.973 s forward: approximately 2.29×. The separate 1/18-column CUDA
comparison agrees with direct propagation within 5.56e-16; the sixteen-stream
CPU checks pass all endpoint finite differences. A single physical tangent
already costs only 1.54× forward in this controlled scene. That does not
establish the target for a complete local basis or a large retrieval layout.

### Equivalent-source sweeps

`source_sweep_operators.jl` tests a leaner association with
`AUDIT_RESPONSE_MODE=sweep`. A backward pass recovers each layer's fixed
incident fields using the forward prefixes and one cached resolvent per layer:

```math
G=(I-R_L^\uparrow R_P^\downarrow)^{-1},\quad
u=G(R_L^\uparrow J_P^\downarrow+T_L^\uparrow U+j_L^\uparrow),\quad
D=R_P^\downarrow u+J_P^\downarrow.
```

The same local operator forcing `dR D + dT U + dj` is contracted to retrieval
**vectors**. A forward source-only adding pass then uses the cached G:

```math
v=G(R_L^\uparrow\delta J_P^\downarrow+f^\uparrow),\quad
\delta J_{\rm new}^\uparrow=\delta J_P^\uparrow+T_P^\uparrow v,
```

```math
\delta J_{\rm new}^\downarrow=f^\downarrow+
T_L^\downarrow(\delta J_P^\downarrow+R_P^\downarrow v).
```

These follow from the same affine interface system and S2014 (C.6); this is
an implementation derivation. Above-layer solar attenuation is restored as
`-j_L dτ_above/μ₀` during each forcing contraction. No suffix cache is needed.
Unlike the independent-layer response, this adding pass scales with retrieval
size, but its tangents are vectors rather than matrices. The initial CPU
1/18-column parity and all endpoint finite-difference checks pass, including U.
The sixteen-stream checks also pass for all three response formulations after
the shared driver was reorganized. At 512 points and 20 layers on CUDA, the
sweep takes 4.191 s for one column and 4.286 s for 60, compared with 1.997 s
and 1.975 s forward (2.10× and 2.17×). The separate 1/18-column tangent
comparison agrees within 5.56e-16. The 60-column sweep saves about 5% relative
to the independent-source response in this regime; it still misses <2×.
At 10,000 points with the smaller 18×18 operator, the sweep takes 3.708 s
for one column and 4.116 s for 60, versus 0.794 s and 0.815 s forward.
The independent-source response previously took 4.017 s for 60 columns;
the sweep therefore has no demonstrated advantage for that large-state,
small-operator scene. This argues for retaining both measured formulations
as experiments, rather than selecting a universal replacement prematurely.

Before integrating either source-response formulation into production, four
contracts need to be established:

- Treat the surface as a differentiable boundary, with its own response and
  retrieval-column map; validate Lambertian albedo before general BRDFs.
- Keep solar attenuation separate from thermal and emission sources. The
  `-J_suffix dτ/μ₀` shortcut is valid only for the solar component.
- Project responses onto the requested endpoints/interior levels and retain
  the solver's Fourier convergence policy and external-solar representation.
- Budget cached subcolumns and layer tangents by operator size, layers and
  wavelengths. State-size independence does not imply small memory use; a
  16-stream, 10,000-point cache needs chunking or recomputation on this GPU.

The controlled core timings allocate reusable workspaces outside the timed
solve for both forward and Jacobian paths. They include constructing the
forward prefixes/suffixes and all layer responses on each call. These timings
therefore have a different preparation boundary from `rt_run` measurements.

## Inverse-cache experiment omitted from production

Reusing cuBLAS inverse pointers and LU pivot/info buffers passed 132 CUDA
checks and the 416-assertion CPU propagation regression. In a same-run
sixteen-stream comparison, however, physical propagation took 26.451 s with
the cache versus 26.321 s without it. Cumulative GPU allocation fell only
from 352.289 GB to 352.025 GB. The local solve took 16.591 s, comparable to
the earlier 16.639 s measurement. This is no demonstrated complete-solve gain.
The extra production fields, dispatch and flag were removed. The isolated
patch, reproduction driver and measurements remain in the evidence directory.
The subsequent cached-sweep timing was stopped once this experiment was
rejected; no result is claimed for that unfinished follow-up.

## Elemental derivative limits

Review found a separate correctness defect in the fused and reference
elemental partials: divisions `j/ϖ` and `j/(ZF₀)` returned zero when the
denominator vanished, even when the derivative was nonzero. The equal-μ
thickness derivative also used `j*(1/τ-1/μ)`, which was undefined at τ=0.
Direct analytic products now retain all three regular limits. A shared
source-factor helper documents SF2023-II (11) and serves both embedded and
external solar tangents. The diffuse albedo partials follow the same rule.

Independent forward-kernel finite differences exposed 28 failures before the
fix; all 160 Float32/Float64 checks pass afterward on each of CPU and CUDA. The physical-Jacobian,
external-solar, interior-sensor, source-routing and local-basis regression run
then passed 1,674 assertions; the 48 selected-layout checks also passed.
Full-solve checks also pass at zero aerosol loading, using second-order
one-sided finite differences. They cover physical/local tangents, both
endpoints with embedded solar, and TOA with external solar. This is significant for an
aerosol-free initial guess: its higher forward phase moments can be zero
while their aerosol optical-depth tangents are nonzero.
This correction is separate from the performance factorization.

## Validation and reproduction

The final CPU suite passes **8,360 assertions**, with 15 tests reported as
skipped/broken (including unavailable GPU checks). The CUDA follow-up passes
108 pointer-plan assertions, 84 small-operator product assertions and 12
geometric-inverse assertions. Pointer tests exercise active-prefix reuse,
inactive-column preservation, mixed 3D/4D operands and copied-cache lifetime.
Both development core drivers now share `core_column_fixture.jl`; after this
extraction, the response experiment again passes all 1/18-column endpoint
finite differences and direct-propagation parity at nonzero U geometry.
Raw samples and stage-specific coverage are retained in
[the evidence directory](evidence/local_basis/README.md).

**Coverage correction:** the initial new local-basis fixtures placed
`external_solar` in YAML, which this integration base does not parse. Those
runs exercised embedded solar in both nominal modes. The existing 73-check
external-solar suite used the constructor keyword and did pass. The new
fixtures now also use `model_from_parameters(...; external_solar=...)` and
assert the actual quadrature flag; external-solar comparisons check TOA and
explicitly unavailable BOA fields. The initial counts did not establish
aerosol external-solar coverage. The corrected
fixtures pass **1,634 CPU and 311 CUDA assertions**, including two aerosol
modes and explicit external-solar phase-column reconstruction. The zero-aerosol
case now covers both solar representations, with TOA-only output
for external solar.

Phase interpolation passed 81 checks on each backend; medium-operator products
passed 192 on each backend. Backend allocation checks cover Float32, Float64
and ForwardDiff.Dual. Metal has not been tested on hardware.

The strict documentation build passes with deployment disabled, including
the backend-workspace additions, after attaching both convergence-strategy
docstrings and listing both in the API manual.
The other numerical and registration findings in the
[release audit](../release_readiness_2026-09-06.md) remain separate release gates.

From `test/`, with a suitable Julia installation and an available GPU:

```bash
VSMARTMOM_JACOBIAN_GPU_TEST=true julia --project=. test_local_jacobian.jl
VSMARTMOM_JACOBIAN_GPU_TEST=true julia --project=. test_jacobian_vendor.jl
julia --project=. test_jacobian_tiled.jl
julia --project=. test_phase_interpolation.jl
AUDIT_BACKEND=cuda AUDIT_NSPEC=10000 julia --project=. \
  ../docs/dev_notes/jacobian_batched/local_basis_benchmark.jl
AUDIT_BACKEND=cuda AUDIT_NSPEC=10000 AUDIT_EXTERNAL_SOLAR=true \
  AUDIT_COMPARE_TILES=true julia --project=. \
  ../docs/dev_notes/jacobian_batched/local_basis_benchmark.jl
AUDIT_BACKEND=cuda AUDIT_NSPEC=512 AUDIT_STREAMS=16 \
  AUDIT_COMPARE_VENDOR=true julia --project=. \
  ../docs/dev_notes/jacobian_batched/local_basis_benchmark.jl
AUDIT_BACKEND=cuda AUDIT_NSPEC=1000 AUDIT_CHUNKS=10 AUDIT_LAYERS=20 \
  AUDIT_GASES=CO2,CH4,H2O AUDIT_NU_MIN=6150 AUDIT_NU_MAX=6250 \
  julia --project=. ../docs/dev_notes/jacobian_batched/local_basis_benchmark.jl
AUDIT_BACKEND=cuda AUDIT_NSPEC=10000 AUDIT_AZIMUTH=37 julia --project=. \
  ../docs/dev_notes/jacobian_batched/local_core_benchmark.jl
AUDIT_NSPEC=5 AUDIT_AZIMUTH=37 julia --project=. \
  ../docs/dev_notes/jacobian_batched/layer_response_check.jl
```
