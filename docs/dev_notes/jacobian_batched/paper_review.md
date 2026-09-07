# Sanghavi core derivatives and the Claude review

Follow-up: the user supplied an older Fortran/C++ implementation. Its
[verified call ordering](fortran_ordering.md) shows delayed molecular-optics
contraction after doubling, alongside early aerosol-parameter expansion.
It provides architectural evidence, not a replacement numerical baseline.

Checked 2026-09-06 against the local paper PDFs and the active implementation on
`origin/sanghavi`, commit `fc74e3f0685b6cdaac9db963278fa3ecfe10cb8d`.

The review's claim that an adjoint is the only route to combined forward and
Jacobian cost below 2× forward is not justified. The papers support factoring
the chain rule through core optical properties, including truncation. However,
neither those equations nor the current Sanghavi Julia implementation establish
that the complete RT Jacobian requires only four scalar tangents independent of
the retrieval state size. The distinction is the representation of the phase
derivative and when the chain rule is contracted.

## What the papers actually report

| Paper | Location | Linearized / forward runtime | Parameters |
|---|---|---:|---:|
| Sanghavi, Martonchik, Davis & Diner (2013), scalar smartMOM | §4.6, p.14, rough ocean | 1.98× | 15 |
| Same, Lambertian | §4.6, p.14 | 1.64× | 14 |
| Same, modified RPV | §4.6, p.14 | 1.29× | 16 |
| Sanghavi, Davis & Eldering (2014), vector vSmartMOM | §5.3, p.423 | approximately 5× | 10 |

Sources: [2013 scalar Jacobian paper](https://doi.org/10.1016/j.jqsrt.2012.10.021)
and [2014 vector Jacobian paper](https://doi.org/10.1016/j.jqsrt.2013.09.004).
These are tangent-linear results. The scalar paper attributes substantial savings
to reusing the forward matrix inverse and explicitly notes that the ratios may
increase with further forward-code optimization. This is particularly relevant
after the recent forward kernel fusion in this repository. Published sub-2×
examples refute an adjoint-only necessity claim; they do not guarantee that ratio
for today's IQU GPU solver, every parameter count, or every requested output.

## Core factorization, including the truncation factor

The scalar paper §3.2.1, (68)–(69), and vector paper Appendix C,
(C.25)–(C.26), factor a microphysical derivative through the single-scattering
quantities. In the vector notation these are the blocks

```math
 x_{SS}=(\tau_i,\omega_i,Z_i,\beta_i),\qquad
 \frac{\partial X}{\partial x_\mu}
 =\frac{\partial X}{\partial x_{SS}}
  \frac{\partial x_{SS}}{\partial x_\mu}.
```

The truncation fraction called β_i there is `fᵗ` in this code. It is distinct
from the Greek coefficients β_l that describe the angular phase expansion.
The paper gives the mixture derivatives in (C.22)–(C.24), truncation/mixing
partials in (C.27)–(C.39), and phase and truncated-coefficient derivatives in
(C.40)–(C.42). The specific δ-M truncation prescription in (C.41) must not be
confused with the separate δ-BGE fitting implementation.

For a raw aerosol with τ, ω, and f, let a=1−fω. The transformed quantities and
their differentials are

```math
 \tau^*=a\tau,\qquad \omega^*=\frac{(1-f)\omega}{a},
```

```math
 d\tau^*=a\,d\tau-\tau(f\,d\omega+\omega\,df),\qquad
 d\omega^*=\frac{(1-f)d\omega-\omega(1-\omega)df}{a^2}.
```

These are already implemented in
`src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties_lin.jl`,
`_createAero_invariant`, including the `ḟᵗ_block` terms. The truncated phase
derivative, including its normalization derivative, comes from
`src/Scattering/truncate_phase_lin.jl`. Downstream RT receives the modified
τ*, ω*, Z* and their tangents; a separate fourth f slot is unnecessary at that
boundary because its effect has already been included. New comments now make
this explicit beside the optical-depth and albedo evaluations.

The complete matrix-side handoff was also traced: `lin_model_from_parameters`
calls `truncate_phase` on the value/tangent pair before storing single-node
optics; for spectral optics it truncates the endpoints and optional reference
node first. `_spectralize_truncated_endpoints` preserves those normalized
tangents both in `lin_greek_coefs` and `phase_lin_greek`. Finally,
`_compute_aerosol_phase_blocks_lin` passes those fields into `compute_Z_moments`
and, when needed, `compute_Z_source_moments`. It does not substitute the raw
Mie coefficients at this boundary. NoTruncation deliberately passes them
through with f=df=0.

`truncation_chain_check.jl` independently perturbs raw Greek coefficients while
holding their zeroth β coefficient normalized, with a nonzero df. It compares
the production truncation tangents against finite differences through the
forward fit: all six normalized Greek families, scalar/IQU Z matrices at
m=0,1,2, and the actual scalar cache with simultaneous τ/ω changes. All 22
assertions pass; the largest absolute Greek-tangent discrepancy is 5.45e-11.
This isolates the truncation chain and is not a validation of the separate
upstream Mie/AOD reference-normalization issues in the release audit.

Computing a reusable core Jacobian and contracting it with a retrieval mapping
later is mathematically valid. At fixed optical state and core basis, its RT
construction can be independent of how many retrieval parameters are subsequently
mapped onto that basis. Producing the final dense state Jacobian still requires
the mapping and output storage. The cost of the core Jacobian itself depends on
the number of layers, angular/phase basis functions, and requested outputs.

## Why four blocks do not imply four scalar derivative arrays

Z is a matrix-valued phase function. Its elemental derivative may be local to
an angular entry, but adding and doubling couple those entries. For example,
if G(Z)=(I−Z)⁻¹, then

```math
 D G(Z)[\Delta Z]=G\,\Delta Z\,G.
```

This is a linear map on phase perturbations, not an entrywise product with one
matrix called ∂G/∂Z. A single perturbed entry generally affects many entries of
G. One can represent that map explicitly, in a chosen phase basis, as a sequence
of linear operations, or via output adjoints. An unqualified three/four-slot
array does not retain it after multiple scattering.

The vector paper itself states after (C.26) that elemental derivatives with
respect to optical thickness and microphysical parameters are propagated through
doubling and adding. Thus its chain-rule statement does not require contraction
to occur only after the full solve.

## What the Sanghavi Julia branch does

The checked remote head was `fc74e3f0685b6cdaac9db963278fa3ecfe10cb8d`;
`git ls-remote origin refs/heads/sanghavi` agreed with the local remote-tracking
commit. The following links are pinned to that snapshot:

- [rt_kernel_lin.jl, lines 197–212](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/fc74e3f0685b6cdaac9db963278fa3ecfe10cb8d/src/CoreRT/CoreKernel/rt_kernel_lin.jl#L197)
  expands derivatives with `lin_added_layer_all_params!` **before** calling
  `doubling!`.
- [doubling_lin.jl, lines 32–64](https://github.com/RemoteSensingTools/vSmartMOM.jl/blob/fc74e3f0685b6cdaac9db963278fa3ecfe10cb8d/src/CoreRT/CoreKernel/doubling_lin.jl#L32)
  obtains `Nparams`, allocates parameter-sized derivative buffers, and loops
  `iparam = 1:Nparams` for inverse and transmission product derivatives.
- `interaction_lin.jl` likewise propagates parameter columns.

These are active code paths, not the commented-out earlier implementation near
the start of `rt_kernel_lin.jl`. The current performance branch batches those
directions; it did not introduce state-vector scaling into an otherwise
state-independent Sanghavi solver. This observation concerns the checked Julia
branch, not an assertion about every historical implementation used for the papers.

## What our measurements do and do not establish

The supplied-core benchmark currently seeds 1, 3 or 15 directions into that
directional engine. Its three common directions are optical depth, albedo, and
a degree-2 Greek phase perturbation; fifteen makes these independent per layer.
It starts after truncation and has no independent raw f direction. It is **not**
a benchmark of an independently assembled, reusable full core Jacobian.
Consequently, its ratios cannot rule out a faster factored core method or be
presented as proof that the requested <2× target requires an adjoint.

The physical-state benchmark is also deliberately limited: five layers, no gas
absorption, and fourteen native columns including five unused gas slots.
It measures the improvements recorded in [iqu_followup.md](iqu_followup.md),
but does not establish the bottleneck ranking for a realistic multilayer gas
retrieval. The review is correct about that limitation.

## How this changes the proposed work

Claude's subsequent refinement proposes the useful intermediate contraction
point: **after local doubling, before atmospheric adding**. At fixed forward
state, write the elemental derivative as a combination of local basis tangents,
`dL/dx_p = Σ_b (dL/dq_b) C_bp`, with coefficients allowed to vary by wavelength
and layer. The derivative of doubling is linear, so

```math
 D\mathcal{D}(L)\left[\sum_b L_b C_{bp}\right]
 =\sum_b D\mathcal{D}(L)[L_b] C_{bp}.
```

Here each phase basis tangent is a complete matrix (and, for external solar,
its matching source column), not an elementwise scalar partial. For fixed
Rayleigh properties and four Mie parameters per aerosol mode, a conservative
basis consists of two scalar core directions, Rayleigh phase, each aerosol's
phase, and four truncated aerosol phase tangents per mode: `3 + 5 NAer`.
Normalization and linear dependence can reduce that count. Additional variable
phase physics or geometry may enlarge it; eight is not a universal dimension.
Gas VMRs and vertical-profile parameters then change the local coefficients,
without adding independent phase directions under these assumptions.

For the solar source, with `E = exp(-τ_above/μ₀)` and fixed geometry,

```math
 dJ_p=E\sum_b (d\widehat J/dq_b)C_{bp}
       -\frac{d\tau_{above,p}}{\mu_0}J.
```

The second term can be appended after doubling because the source is linear in
the incident beam. Local beam attenuation still has to be differentiated inside
doubling with the core optical-depth direction. Sources depending on additional
parameters require their own basis terms; the solar result is not automatically
a complete rule for all source types.

The executable `core_basis_check.jl` verifies this contraction order using the
actual doubling implementation: eight independent tangents expanded to sixty
state columns early versus after doubling, with wavelength-dependent coefficients
and independent above-layer attenuation derivatives. All 96 CPU assertions pass
for scalar/IQU, Float32/Float64 and three/six doubling steps. This is an algebra
experiment, not yet a production basis builder or an end-to-end speed benchmark.

Under this design the expensive angular phase assembly and local doubling can
avoid a retrieval-state-sized matrix axis. Scalar coefficient construction,
contraction, dense adding and output still depend on state size. Thus even the
assembly is not literally independent of Nparams in every operation.
The suggested ~1.7×/~4× savings are estimates of propagation work, not measured
full-solve gains: doubling and adding have different costs, and contraction,
workspace and source work must also be counted.

The revised assertion that *only* an adjoint can make adding independent of
Nparams is still too categorical. A global core basis with later contraction
is another mathematical representation, although it can be much larger and more
expensive. Dense adding remains state-sized in the specific local-basis design
above; an adjoint is an attractive alternative when there are few output seeds.

One historical detail also needs qualification. The scalar paper §4.1 reports
301 viewing directions, but opposite azimuths share zenith cosines. Its earlier
discussion allows either interpolation or adding zero-weight angular nodes.
Those passages do not establish that the reported run inserted 301 distinct
nodes into one operator. The inverse-dominated timing is explicitly reported;
that particular explanation of the exact matrix dimension remains an inference.
Likewise, a larger-stream modern IQU speed ratio needs measurement.

The next architectural comparison should include a factored optical-property
boundary and a compact layer-local core/phase basis, with retrieval contraction
placed where it is cheapest. Gas-profile sparsity can reduce local doubling work;
overlying direct-beam attenuation and the accumulated atmosphere still need their
correct source/adding derivatives. Factoring phase assembly alone reduces tensor
allocation but does not automatically eliminate state-column propagation.

An adjoint remains a useful alternative when few output radiances are requested,
but is not a prerequisite established by either paper. Compare both cost models
after profiling an absorbing, multilayer retrieval and accounting for the phase
basis dimension, output count and memory use. No adjoint or state-independent
core-Jacobian assembler has been implemented by this follow-up.

The review's interpolation ranking also needs measurement on that workload.
`_interpolate_phase_nodes` does allocate per wavelength, but it is conditional on
phase-node optics; the three-node route is natural-cubic interpolation, so a
replacement must preserve that rule. The completed phase-table and backend-zero
changes overlap some allocation concerns, but do not implement the proposed
factored phase boundary. Dead reductions, excess workspaces and repeated uploads
remain reasonable code-level targets without implying that each dominates runtime.
