# Jacobian architecture and extension review

Reviewed 2026-09-07 on `perf/jacobian-batched-propagation`. This is a review
of the current implementation and a staged design proposal, not a claim that
all proposed extension interfaces already exist.

## Assessment

Keep the analytic MOM core and its supplied-tangent interface. The expensive
solver already differentiates optical operators independently of retrieval
names. Local optical factorization and equivalent-source adding exploit that
separation successfully. The remaining extensibility work belongs mainly at
the model-construction, parameter-mapping, source, and instrument boundaries.
A larger generic rewrite of elemental/doubling algebra is not justified by
this review.

The current implementation is extensible for new selections of its existing
physical parameters. It is not yet a general derivative provider for arbitrary
new atmospheric physics, source families, geometry, or retrieval coordinates.
In particular, an abstract type or AD-mode declaration is not evidence that
all propagation paths implement its derivatives.

## The actual chain rule

Let `x` be retrieval coordinates, `p` physical model inputs, `q` layer optics
and source values, `R` the RT output, and `y` the instrument measurement:

```text
x → p → q = (τ, ϖ, Z++, Z-+, direct-beam phase columns, source terms)
      → elemental → doubling → adding → views/Fourier reconstruction → y
```

For a fixed instrument operator H, the atmospheric contribution is
`Kx = H · DR(q) · Dq(p) · Dp(x)`. If instrument parameters are retrieved,
add `(dH/dx)R` and any offset/source terms. Reordering columns is not a
coordinate transform, and a coordinate transform need not be diagonal.

| Boundary | Current implementation | Extension responsibility |
|---|---|---|
| Retrieval identity | `ParameterKey`, `ActiveParameterLayout`, `JacobianPlan` | Stable names and band maps; units and coordinates remain external |
| Optical derivatives | `lin_model_from_parameters.jl`, scattering AD, pressure/profile formulas, absorption tangents | Supply every physical dependency, including truncation and interpolation |
| Optical compression | `local_jacobian{,_cache}.jl` | Factor physical optical tangents into local seeds and coefficients |
| MOM propagation | `CoreKernel/{elemental,doubling,interaction}_lin.jl` | Propagate supplied directions; no retrieval-specific branches |
| Source adding | `CoreKernel/source_adding_lin.jl` | Supported solar/SIF Lambertian endpoint configurations only |
| Measurement | Study `VSmartMOMForward.jl` | Stokes/instrument mapping, unit conversion, log coordinates, surface basis transform |

The local representation is `Dq = B C`: double the small seed basis `B`,
then contract with `C`. Matrix adding contracts complete doubled operator
tangents before adding. Equivalent-source adding instead contracts the
forcing vectors formed with the solved illumination. Both preserve angular
coupling; multiplying a phase derivative by an elemental scalar after the
complete solve would not be equivalent.

## Concrete gaps and changes in this review

### Disabled upstream derivatives could be selected

Flavor traits can disable Mie microphysics or the q-driven H2O tangent block.
The constructor previously accepted a plan selecting these disabled columns;
structural zero placeholders could then masquerade as valid sensitivities.
Construction now checks native selections against the **effective** upstream
options, including keyword overrides. It rejects incompatible requests and
invalid band/native-column counts. This checks declared availability, never
whether derivatives happen to be numerically zero at one state.

The low-level `PlannedRTModelLin(base, plan)` wrapper still accepts supplied
optical tangents: callers constructing it manually own their provenance.
Do not infer availability from zeros in these arrays; tests deliberately
supply synthetic gas tangents through this interface.

### Forward source support did not establish linearized support

`source_ad_mode` defaults to `AnalyticSourceJacobian`. Thermal emission has a
forward source-slot implementation, but its full linearized propagation is
unfinished. The linearized driver currently constructs legacy solar layer
sources and handles SIF at the surface; it does not pass arbitrary prepared
sources through every tangent stage. It now rejects unsupported source types,
including compositions containing them, before allocating layer workspaces.

The internal `_linearized_source_supported` gate is separate from the
narrower source-adding gate. Enabling a new type requires implementing and
validating its complete path, not merely overriding this predicate. Reserved
`ForwardDiffSourceJacobian` vocabulary does not itself implement AD.

## Priorities for future development

### 1. Freeze retrieval meaning, validate model compatibility

Compile active keys once for an inversion. The OCO flavor currently selects
CO2 layers from trial-model geometric centers below 10 km and recompiles on
every model construction. Near that threshold a changing pressure/profile
can change the active set. A length check catches some changes but cannot
prove column identity. Store and compare the ordered keys and per-band maps,
or pass a frozen selection into construction. A deliberate change of vertical
grid must require recompilation and explicit state/prior remapping.

Keep plan metadata separate from trial buffers. A reusable compiled plan
should record the component ordering, vertical grid identity, spectral bands,
surface/source parameterization, and backend/precision requirements relevant
to each cache. It must not reuse state-dependent optical values by accident.
Current public plan vectors are mutable; avoid treating object identity as a
sufficient cache key.

### 2. Replace native positional assumptions at the optical boundary

The native aerosol block has seven fields, with microphysics in positions
2:5. The local phase basis assumes four LogNormal/Mie microphysics directions,
a fixed Rayleigh phase reference, and layer-independent species phase shapes.
The gas block reserves its first species for q-driven H2O. These assumptions
appear in `parameter_layout.jl`, `optical_jacobian_cache.jl`, and both local
basis files, so adding a field to `ParameterKey` alone cannot extend physics.

Introduce component-provided descriptors for native physical directions and
an optical-tangent provider that returns values, tangent identities, seeds,
and coefficients. Let each provider use analytic formulas, upstream AD, or
precomputed derivatives. Keep dense supplied optical tangents as the general
reference/fallback for providers that cannot factor into a shared phase basis.
Compile the upstream work request from these descriptors; the current
all-or-none microphysics trait cannot avoid work separately for each species.

Examples requiring explicit extension are layer-dependent particle size,
nonspherical/cloud phase providers, Rayleigh depolarization derivatives, and
temperature changes coupling opacity, density, geometry, and thermal emission.
A shared fixed phase basis must not silently omit their extra directions.

### 3. Make source tangents explicit and composable

A future source handoff should carry prepared forward source values and
physical/source-coordinate tangents with named column maps. Propagation must
retain each source's attenuation law and surface behavior. For thermal
emission, temperature changes both the Planck source and optical properties;
for SIF, surface emission must not acquire the direct-solar overburden factor.
The existing SIF regression protects that distinction.

Generalize source-adding capabilities only after the matrix reference handles
the new source, surface, and observer combination. A new BRDF also needs its
reflection/transmission derivatives and coupling to source illumination.
Do not infer Jacobian support from forward surface parameter count alone.

### 4. Keep retrieval transforms and instrument derivatives outside MOM

Use an explicit physical-to-retrieval map `Dp/Dx`, with units and reference
wavelengths attached at the adapter boundary. The study's fixed-microphysics
identity `dF/dlog(AOD760) = τ_ref dF/dτ_ref` needs additional cross terms if
microphysics changes the extinction ratio at 760 nm. Profile EOFs, shared
parameters across species/bands, and constrained states require non-diagonal
maps. Apply their covariance transformations consistently.

Instrument shifts, line-shape width, calibration gains, and polarization
mixing should contribute at measurement formation. Computing RT columns for
parameters that affect only H would waste work and blur ownership.

### 5. Support arbitrary seed directions before adding a second solver

The supplied-tangent core can underpin Jacobian-vector products by accepting
seed matrices at the optical/retrieval boundary. Start with exact equivalence
to multiplying the assembled Jacobian by a seed matrix. Parameter chunking
can bound memory for large states; spectral chunking must preserve full-band
Fourier stopping decisions and instrument convolution overlap.

A reverse/adjoint path may eventually help scalar objectives with many state
parameters, but it needs its own transpose propagation, storage/recomputation
policy, and dot-product tests. Equivalent-source adding is a forward tangent
method, not an implemented adjoint. Dense OE still needs K or equivalent
information-matrix products; an adjoint is not automatically a replacement.

## Validation contract for an extension

1. Check physical values and upstream tangents independently, including the
   normalization, interpolation, and truncation chain rules.
2. Check RT forward/linearized forward parity with identical numerical
   settings. Adaptive Fourier, truncation, doubling counts, and grid changes
   introduce discrete boundaries; diagnose finite differences on both sides
   and state what is held fixed.
3. Compare compact/local/source paths against direct physical propagation.
   Use directional finite differences in addition to column checks, across
   more than one step size. A shared reference can share the same upstream bug.
4. Verify named maps, units, inactive exact zeros, and all requested source
   contributions. Include zero-loading and mixed fixed/active species cases.
5. Test the backend/precision actually supported. CPU/CUDA tests do not prove
   Metal parity or memory behavior. Benchmark dimensions and optical basis
   explicitly; do not promise a universal '<2× forward' cost.
6. For retrieval-facing changes, compare full convergence, rejected steps,
   final cost/state, posterior products, and trial-template isolation.

Read-only LUT sharing is opt-in and already tested across A→B→A trial
sequences. Future caches must declare their dependencies and lifetimes:
geometry, spectroscopy, and fixed microphysics may be reusable under explicit
conditions; pressure/temperature/VMR-dependent absorption and mixed layer
optics must be recomputed when their inputs change. Shared LUT mutation is
still a caller violation, not prevented by an immutable outer Julia struct.

## Evidence

The paired retrieval driver and results are in
[`evidence/convergence/`](evidence/convergence/README.md). The source study is
read only. The new guards have an independent extension-contract test in
`test/test_jacobian_upstream_contract.jl`; existing selective Jacobian and
source tests provide supported-path regression coverage.
