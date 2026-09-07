# IQU Jacobian follow-up

The smaller IQU improvement had four concrete causes: repeated polarized phase
construction on the host, large zero tangent tensors staged through host memory,
redundant operand loads in the GPU matrix products, and an older inversion
path that the forward solver had already replaced. The full physical-state
benchmark and the full solve from supplied core optics measure different work;
neither should be described using the other's forward/Jacobian ratio.

## Measured physical-Jacobian result

On the same 10,000-wavelength Float64 IQU fixture (five layers, one aerosol,
Lambertian surface, embedded solar, fourteen full/four selected columns):

| Implementation | Full Jacobian | Selected Jacobian |
|---|---:|---:|
| Original reference, historical median | 38.268 s | 8.164 s |
| Previous branch result, after surface batching | 16.441 s | 4.119 s |
| This follow-up, including fused inverse | **7.131 s** | **1.873 s** |

The full result improves 2.31× over the previous branch and 5.37× over the
original reference. Selected columns improve 2.20× and 4.36× respectively.
Full host allocation falls from 25.88 to 3.64 GB. Full cumulative device
allocation is 121.00 GB (not peak memory), down from 1,911.63 GB in the original
reference. The fresh reference/optimized comparison agrees to 3.82e-17 in
radiance and 6.67e-16 in Jacobians; selected dR matches the corresponding full
columns exactly.

Timings are medians of three warmed synchronized samples on an A100 PCIe 40 GB,
Julia 1.12.6, four Julia threads and one BLAS thread. The machine was shared;
forward medians varied across runs. Use the raw samples and allocation/kernel
measurements alongside these historical end-to-end comparisons. The latest
forward-only median for this fixture was 0.573 s. The separate core-optics
fixture below has a different phase model and Fourier count; its times must
not be subtracted from this table to estimate upstream cost.

The matching scalar fixture now takes 1.082 s full and 0.414 s selected,
with a 0.141 s forward median. Full host allocation is 0.587 GB, cumulative
device allocation 14.20 GB, and maximum absolute Jacobian disagreement with
the reference is 4.45e-16. These remain physical-state timings.

## What the profile found

For the earlier 10,000-wavelength IQU implementation (16.44 s warmed median),
the separate instrumented run spent 6.27 s in optical properties, 5.93 s in
doubling, 1.86 s in interaction, and 2.17 s creating layer arrays. These are
instrumented section timings, not an additive decomposition of the warmed median.
Host allocation was 25.88 GB per full physical-Jacobian solve.

The polarized phase routine rebuilt spherical-function tables for every
wavelength and copied small matrix slices for each angular pair, harmonic degree
and microphysical tangent. Scalar phase evaluation did not incur those matrix
slice allocations. The replacement builds the angular tables once and accumulates
whole static Stokes blocks. It shares the exact same angular evaluator between
Greek coefficients and their tangents, exploiting linearity in the coefficients.

Several phase-lifting, gas-mixing and layer constructors also allocated CPU zeros
before uploading the full tensors. They now allocate and zero directly on the
owning backend and preserve the requested floating type. Together, these changes
reduced the full IQU median to 11.34 s and host allocation to 3.64 GB.

The original GPU tangent products assigned one output element per thread, but
reloaded the same operands for every dot product. A CUDA workgroup now stages
one square matrix for a wavelength/parameter pair in shared memory. The isolated
Float64 product-rule benchmark improved from 7.48 to 3.44 ms at N=18 and from
27.43 to 7.82 ms at N=32, with 10,000 wavelengths and fourteen tangents. Vendor
BLAS was only modestly faster than the old N=18 kernel, and was slower for
smaller and source-vector products, so a blanket BLAS replacement was unsuitable.

The shared-memory dispatch is restricted to CUDA, square N=16–32 and at least
512 wavelengths. Source-vector products and smaller operators keep the simpler
kernel. Float32 and Float64, inactive prefixes, mixed forward/tangent operands,
and accumulation into an existing output are covered by tests. Metal has not
been hardware-tested.

## Reusing the forward solver's fused inverse machinery

GPU event timings exposed a cost that was understated in the host profile.
At N=18 and 10,000 wavelengths, cuBLAS inversion took 4.80 ms, the portable
pivoted LU inverse took 0.43 ms, and the existing fused product/right-solve
kernel took 0.56 ms. The latter timing also includes forming the product.
Asynchronous GPU work can be charged to a subsequent synchronization point in
host section timers; the small host times under getrf/getri were misleading.

The Jacobian path now reuses the forward kernel: solving
G(I−R₁R₂)=I gives the geometric inverse G required by all tangent directions.
This combines product construction and pivoted solve into one kernel after
preparing the identity RHS. It avoids the old cuBLAS getrf/getri path for CUDA
operators up to 32×32, while honoring the forward fused-solve enable flag.
The CPU and other backend paths retain the existing inverse implementation.

## The algebraic saving

At a fixed wavelength and Fourier order, write Q=R₁R₂,
G=(I−Q)⁻¹ and H=TG. The inverse rule and product rule give

```math
 dG = G(dQ)G, \qquad dH = (dT)G + T(dG)
     = [dT + H(dQ)]G.
```

The last expression saves two tangent matrix products and the dG buffer. Matrix
order is preserved; only associativity and distributivity are used. Doubling
sets R₁=R₂ to the D-transformed reflection. General adding evaluates the two
reflection orderings separately. Every source and operator tangent still
propagates through the full multiple-scattering calculation.

## Scientific comments and equation references

The comments adjacent to the evaluations now show the finite-thickness elemental
partials, upstream/core chain rule, inverse/product rules, both adding directions,
solar attenuation derivative, D-matrix transformation, phase expansion and
normalized scatterer mixing. Two older explanations were corrected: the upward
and downward inverses have different factor orderings, and the current elemental
kernel differentiates finite-thickness solar formulas rather than literally
implementing the infinitesimal/thermal expressions in the 2014 appendix.

Equation numbers were checked against the local PDFs:

- [Sanghavi, Davis & Eldering (2014), JQSRT 133, 412–433](https://doi.org/10.1016/j.jqsrt.2013.09.004):
  (23)–(28) adding; (C.6) product rule; (C.7) inverse rule;
  (C.11)–(C.16) adding tangents; (C.17)–(C.20) D symmetry;
  (C.22)–(C.26) mixing/chain rule; (C.40) phase tangents.
- [Sanghavi (2014), JQSRT 136, 16–27](https://doi.org/10.1016/j.jqsrt.2013.12.015):
  (14)–(16) phase expansion, spherical-function matrices and Greek matrices.
  This is the distinct Fourier-expansion paper, not the JQSRT 133 paper.
- [Sanghavi & Frankenberg (2023), Part II, JQSRT 311, 108791](https://doi.org/10.1016/j.jqsrt.2023.108791):
  (10)–(11) finite-thickness elastic operators and solar sources;
  (12) elastic adding relations.

The factored dH expression above is an algebraic consequence of the cited rules,
not a separately numbered equation attributed to a paper. The maintained code
map is [theory_references.md](../theory_references.md).

## Core-optics benchmark and the <2× target

`core_properties_benchmark.jl` measures a complete five-layer atmospheric solve
from supplied core optical properties, including elemental construction, every
doubling/adding step, solar SFI, all Fourier orders m=0:2, and TOA/BOA diffuse
radiance reconstruction. The lower boundary is black. It excludes model/Mie
construction, optical mixing, parameter transformation, input uploads, cumulative
optical-depth preparation, and persistent workspace allocation from both timings.
The timed inputs and scratch are already resident on the selected backend.

The Float64 scene has total optical depth 0.305, single-scatter albedo 0.9,
Rayleigh-form baseline phase matrices, three weighted streams, embedded solar
at SZA=45°, and views at 0° and 30°. The operator is 6×6 for I and 18×18 for IQU.
It supplies these explicit derivative bases:

- One common optical-depth perturbation across all five layers.
- Three common perturbations: optical depth, albedo, and degree-2 Greek β.
- Fifteen independent perturbations: those three quantities in each layer.

The β₂ perturbation preserves β₀ and tests a phase-shape direction through the
full solve. It is **not** the dense derivative with respect to every Z entry.
The inputs are post-truncation quantities; no independent raw fᵗ perturbation
is supplied. This measures the existing directional engine, not an independent
reusable core-Jacobian assembler.
A full matrix-valued Z Jacobian has more than three scalar columns. The legacy
three-slot doubling helper cannot justify a full core-optics speed claim: the
entry-local elemental Z partial becomes a matrix-coupled derivative after
multiple scattering. Production contracts supplied directions before doubling.

The <2× development target concerns this full core-optics solve, not just the
elemental partials and not the physical-state benchmark including upstream work.
It remains an unmet target for all the directional cases measured
here. Inverse reuse alone is insufficient: propagation still performs matrix
products per direction, with additional source, memory and launch costs.

Final CUDA IQU measurements at 10,000 wavelengths, after the fused-inverse change:

| Supplied directions | Forward | Forward + tangents | Ratio |
|---|---:|---:|---:|
| 1 common τ | 0.189 s | 0.491 s | 2.60× |
| 3 common τ, ϖ, phase | 0.187 s | 0.924 s | 4.93× |
| 15 independent layer directions | 0.185 s | 3.529 s | 19.10× |

These ratios are not lower bounds for a different core-Jacobian representation.
The [paper and code review](paper_review.md) explains the distinction and verifies
that fᵗ already enters the upstream chain rule. Factored phase/core bases and
layer-local sparsity should be evaluated alongside an adjoint, rather than
concluding that only an adjoint could meet the target.

For a dense core-optics Jacobian with few requested radiances, a radiance-adjoint
implementation is a candidate for subsequent work: propagate output sensitivities
back through adding/doubling, then contract them with elemental core partials.
Its cost would scale with the requested output seeds rather than all input phase
entries. This would require its own implementation, performance measurement and
finite-difference validation; no <2× guarantee follows from that proposal.

## Validation and reproduction

- The default suite passed 6,749 tests with three existing broken/skipped cases
  after phase tables, backend allocation and tiled products were added. That
  full run preceded the subsequent dH factorization and fused-inverse changes.
- After dH factorization, 158 physical-Jacobian finite-difference, selective,
  external-solar CPU and Float32 regressions passed. The phase/operator run
  passed 2,188 CPU/CUDA assertions, including bitwise phase-reference parity.
- After the final fused-inverse change, 888 operator, allocation, tiled-product
  and inverse-residual assertions passed on CPU/CUDA. The external-solar CUDA
  suite passed all 73 tests with scalar indexing disabled.
- The supplied-core CPU test checks all 19 seeded directions at both TOA and
  BOA against central differences, including azimuth 37° to exercise nonzero U.
  Full physical 10,000-wavelength scalar/IQU runs also check reference R/T/dR/dT
  parity and selected-column parity. These do not resolve the separate
  scientific issues recorded in the release audit.
- The subsequent paper review adds 96 checks of delayed basis contraction
  through doubling, including solar attenuation, and 22 independent finite
  differences through the scalar and matrix parts of truncation. See
  [paper_review.md](paper_review.md) for their scope.

Run from `test/`, using the test environment and a working Julia binary:

```sh
AUDIT_BACKEND=cuda AUDIT_NSPEC=10000 AUDIT_POL=IQU JULIA_NUM_THREADS=4 julia --project=. ../docs/dev_notes/jacobian_batched/core_properties_benchmark.jl
AUDIT_BACKEND=cpu AUDIT_NSPEC=16 AUDIT_POL=IQU AUDIT_AZIMUTH=37 JULIA_NUM_THREADS=4 julia --project=. ../docs/dev_notes/jacobian_batched/core_properties_benchmark.jl
VSMARTMOM_JACOBIAN_GPU_TEST=true JULIA_NUM_THREADS=4 julia --project=. -e 'include("test_jacobian_batched.jl"); include("test_z_jacobian_tables.jl")'
EXTERNAL_SOLAR_TEST_GPU=1 JULIA_NUM_THREADS=4 julia --project=. -e 'using CUDA; CUDA.allowscalar(false); include("test_external_solar_sfi.jl")'
```

The physical benchmark command and fixture are documented in [README.md](README.md).
Final raw measurements are `evidence/iqu-fused-inverse-final.log`,
`evidence/scalar-fused-inverse-final.log`, and
`evidence/core-iqu-fused-inverse.log`; older stage logs are retained and named
separately. Structured physical/core results are appended to `results.json`.
The server was shared and no Metal hardware validation was available.
