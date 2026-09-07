# Equivalent-source adding in the Lambertian retrieval path

The opt-in `jacobian_adding=:source` path runs through the ordinary linearized
`rt_run`, with `jacobian_basis=:local`. It retains the established elemental and
doubling kernels, optical-property assembly, Fourier convergence, and endpoint
postprocessing. Matrix adding remains the default.

```julia
result = rt_run(model, lin_model, NAer, NGas, NSurf;
                jacobian_basis=:local, jacobian_adding=:source)
```

Supported configurations are elastic solar illumination, scalar/Legendre
Lambertian surfaces, and endpoint observers. Embedded solar provides TOA and BOA;
external solar retains the existing TOA-only contract. Other source and observer
configurations are rejected explicitly. The keyword is also forwarded by
`rt_run_lin` and the planned-Jacobian entry point.

## Mathematical order

1. Build each layer's local optical basis and coefficient tensor `C[s,b,p]` as in
   [the local-basis report](local_basis.md). This includes both the scalar
   truncation-factor chain and derivatives of normalized truncated phase moments.
2. Double the complete local matrix/source directions, holding the *value* of
   above-layer solar attenuation fixed and its tangent zero.
3. Retain forward layer operators and forward prefixes. A backward interface
   solve recovers diffuse illumination `D,U` at both faces of each layer.
4. At these fixed incident fields, form equivalent-source perturbations:

   ```text
   f⁻ = dr⁻⁺ D + dt⁻⁻ U + dj⁻
   f⁺ = dt⁺⁺ D + dr⁺⁻ U + dj⁺.
   ```

   This is S2014 (C.6), applied to the affine layer balance. Contract these
   **vectors** with `C`, then append `−j dτ_above/μ₀` once. No elementwise
   post-doubling phase-matrix chain rule is used.
5. Propagate the forcing vectors through the fixed forward column. For prefix P
   above layer L, with `G=(I−r_L⁻⁺ R_P⁺⁻)⁻¹`, the source recurrence is

   ```text
   v = G (r_L⁻⁺ δJ_P⁺ + f⁻)
   δJ_new⁻ = δJ_P⁻ + T_P⁻⁻ v
   δJ_new⁺ = f⁺ + t_L⁺⁺ (δJ_P⁺ + R_P⁺⁻ v).
   ```

   These are S2014 (27)–(28) / SF2023-II (12) applied to equivalent sources.
   The backward/forward ordering is our derivation, not a claim about the
   implementation in the papers or historical Fortran.
6. The surface is the last affine operator. Only its albedo directions are
   constructed as matrix tangents. Their forcing vectors are scattered into
   retrieval slots, and the total-column solar-attenuation tangent is appended
   once. The surface convention `t⁺⁺=I` preserves BOA downwelling.

The implementation is in
[`source_adding_lin.jl`](../../../src/CoreRT/CoreKernel/source_adding_lin.jl).
Workspace allocation, layer construction, incident-field recovery, forcing
formation, and source propagation are separate functions with equations above
the key evaluations. The equation references follow the
[verified reference map](../theory_references.md).

## Cost and limits

Atmospheric matrix directions scale with the local optical basis. Surface matrix
directions scale with the number of albedo coefficients. Coefficients, source
vectors and radiance outputs still scale with retrieval length: this is tangent
mode, not an adjoint. The implementation caches local layer tangents for the
backward pass, so peak memory grows with layers × wavelengths × local directions.
A spectral batch must fit memory; wavelength chunking remains useful.

The <2× full-solve target is not a guarantee. Complete-model benchmarks include
optical mixing and allocations inside `rt_run`; Mie and spectroscopy are prepared
before timing. Use `AUDIT_COMPARE_SOURCE=true` with `local_basis_benchmark.jl` to
compare local matrix adding, source adding, and forward solves on the same scene.
The embedded forward baseline is the public `rt_run`, including its usual
HDR/BHR diagnostics; the linearized call returns endpoint radiances and their
Jacobians. External solar uses the lean `rt_run_toa` forward entry point.

## Real absorption benchmark

A100 PCIe 40 GB, Float64 IQU, three weighted streams, two views (second azimuth
37°), 20 layers and 69 columns. CO₂/CH₄/H₂O contribute 60 nonzero layer-gas
columns. The 6150–6250 cm⁻¹ grid has 10,000 points processed as ten 1,000-point
chunks, each preparing its own Mie interpolation anchors before timing. All
modes use the same chunk inputs. Times are medians of three summed warmed
runs, including allocations and complete `rt_run` work.

| Solar representation | Local matrix adding | Source adding | Forward | Source/forward |
|---|---:|---:|---:|---:|
| Embedded, TOA + BOA | 22.6374 s | 4.8919 s | 1.6669 s | 2.935× |
| External, TOA only | 15.6473 s | 4.3402 s | 1.6165 s | 2.685× |

The embedded source path is 4.63× faster than local matrix adding. Maximum
Jacobian disagreement is 1.78×10⁻¹⁴. Per-chunk CUDA pool used-memory high-water
marks are approximately 11.356 GB (matrix), 7.125 GB (source), and 0.726 GB
(forward). These include allocated objects awaiting Julia garbage collection;
they are not the minimum live workspace sizes or total device/context memory.
They are measured by resetting `MEMPOOL_ATTR_USED_MEM_HIGH` before each timed
call and reading it afterward. Cumulative allocated bytes are recorded separately.
The <2× target remains unmet; this benchmark also includes optical assembly
rather than starting from already mixed core properties.

External solar gives a 3.61× speedup over local matrix adding, with a maximum
Jacobian disagreement of 2.13×10⁻¹⁴. Its per-chunk pool high-water marks are
9.130 GB (matrix), 5.693 GB (source) and 0.603 GB (forward), with the same
measurement caveats. These two rows have different endpoint output contracts;
they are not an isolated comparison of solar representations.

### Larger angular operator

With 16 weighted streams (57×57 IQU operators), 20 layers, 69 columns and
512 wavelengths spanning the same real-absorption band, source adding takes
29.3663 s versus 11.0652 s forward: **2.654×**. These are warmed complete
`rt_run` calls on prepared optics. Forward radiances agree within
5.21×10⁻¹⁸. The 69-column matrix-adding reference was not run at this size;
the smaller CUDA regressions provide Jacobian parity and finite differences.

Pool high-water marks are 21.553 GB combined and 18.727 GB forward, with the
allocator caveats above. Cumulative device allocations are much larger:
157.876 GB and 18.727 GB respectively. A separate instrumented source trace
shows 9.69 s in the leading GEMM kernel, 4.02 s in LU factorization, 3.40 s in
matrix inversion, and 4.31 s in the two leading 4D broadcast kernels. These
trace times are not a partition of the uninstrumented median. The remaining
work is concentrated in matrix propagation and allocations; reducing phase
interpolation alone would not address those costs.

## Baseline model construction plus the complete solve (before caching)

The `AUDIT_END_TO_END=true` run independently rebuilds both models from the
same parameters. It includes parameter parsing, Mie, absorption, upstream
optical-property derivatives for LinMode, and the full surface-coupled solve.
JIT compilation and artifact downloads are warmed. The scene is the embedded
20-layer, 69-column IQU case above, again totaling 10,000 wavelengths in ten
independently prepared 1,000-point chunks.

| Timing boundary | Forward only | Forward + Jacobians | Ratio |
|---|---:|---:|---:|
| Independent construction + complete solve | 28.1427 s | 32.4561 s | **1.153×** |

The maximum forward-radiance difference between independently constructed
models is 4.16×10⁻¹⁷. This meets the <3× end-to-end milestone for this scene
with the opt-in source-adding path. It does not establish a universal ratio
across stream counts, aerosol populations, or output/source configurations.

These timings precede spectroscopy caching. The [construction optimization
report](construction_cost.md) identifies repeated HITRAN parsing as roughly 90%
of this setup cost and replaces it with shared parsed data. Updated totals are
3.174 s forward / 6.473 s with Jacobians (2.039×) for cached direct HITRAN, and
1.834 / 5.263 s (2.870×) with caller-owned LUTs. The latter spectroscopy setup
has different table-generation assumptions; see that report for scope.

Construction dominates this particular full-rebuild workload: synchronized
construction subtimers are 26.1563 s (forward) and 24.8286 s (LinMode); solve
subtimers are 2.0001 s and 7.7210 s. These are medians of summed samples, so
their medians need not add to the median total. These solve subtimers include
any garbage collection immediately following construction; the cause of their
higher solve times has not been isolated. They should not replace the separately
warmed prepared-model comparison above. Retrieval workflows that reuse prepared models have that
different cost balance; the original <2× core-solve target remains open.

Raw timings and reproduction commands are in
[the evidence directory](evidence/source_adding/README.md).

## Spectral nodes and shared phase matrices

An aerosol's phase matrix at a given wavelength is shared across layers. Each
layer changes its Rayleigh/aerosol scattering fractions, not the aerosol phase
matrix itself. The local Jacobian basis is shared across layers too. Spectral
nodes reduce angular evaluations; their phase matrices and tangents are expanded
once onto the wavelength grid for each Fourier order. Full mixed layer matrices
are then materialized for the elemental kernels.

`phase_storage_benchmark.jl` isolates these costs on the A100 with 10,000
wavelengths **in one batch**, 20 layers, Float64 IQU, three weighted streams,
one aerosol, two phase nodes, embedded solar and Fourier order m=1. Absorption
is synthetic; this is a phase-storage profile, not a spectroscopy benchmark.
The entries are separate warmed measurements (medians of five), not additive
subtimers of one call.

| Operation | Time |
|---|---:|
| Evaluate angular phase matrices and tangents at nodes | 0.494 ms |
| Evaluate nodes and expand across wavelengths | 1.391 ms |
| Construct shared local phase basis | 2.640 ms |
| Materialize mixed phase matrices for all layers | 9.478 ms |
| Complete phase assembly | 13.670 ms |

The shared phase basis contains 362.88 MB; the mixed layer phase arrays contain
1,036.80 MB. The profile asserts that all layer Jacobians reference the same
basis object. Thus node interpolation adds some work, but materializing layer
mixtures costs more in this scene. A useful later optimization is to consume
the shared basis and mixing fractions directly in the elemental kernel. That
would remove mixed-array storage while preserving wavelength-dependent aerosol
physics. It is not implemented by this change.

## Correctness repairs accompanying integration

The full linearized constructor now uses the same fixed configured `n_ref` as
the forward and selective constructors. For `q=k/k_ref`, the retrieved aerosol
index changes only the numerator; native lognormal size coordinates change
both. No artificial unity anchor is inserted when the reference wavelength is
inside the band. Such an anchor contradicted the common-reference convention
for species with different refractive indices and could produce negative
Float32 optical depths near a band endpoint. Nearly coincident phase knots
are also suppressed at the spectral grid's precision.

Normalized Mie quadrature weights now retain the raw normalization `S` when
differentiated: `dw=(du-w*sum(du))/S`. The size directions are derivatives with
respect to the native LogNormal parameters `μ=log(median radius)` and
`σ=log(geometric width)`, not the exponentiated YAML inputs.

The Legendre linearized surface now follows the scalar surface's diffuse-source
convention. Its BOA collimated carrier and attenuation derivative are appended
once outside the Fourier sum. The common forward surface scaffold retains its
carrier for hemispheric diagnostics. Tests include a view on the embedded solar
ordinate so an extra or missing carrier cannot hide in off-solar comparisons.

Validation completed with 8,501 CPU passes and 15 existing skips/broken tests
across contiguous suite segments, plus 92/92 targeted CUDA checks and a strict
Documenter/Vitepress build. The source regression covers scalar/IQU,
Float32/Float64, embedded/external solar, two aerosols, both supported surface
types, independent forward parity, gas finite differences and a zero-albedo
finite difference. Metal was not hardware-tested. Historical failures, their
repairs and the final successful logs are identified in the evidence README.
