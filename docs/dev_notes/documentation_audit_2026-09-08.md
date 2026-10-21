# Documentation and v2.2 follow-up audit — 2026-09-08

This audit follows the local release integration at `326c1e91`. It changes
manual/tutorial text, docstrings, and documentation checks. It does **not**
change radiative-transfer algorithms or complete the missing derivatives below.
Scientific formulation and linearization credit remains with Suniti Sanghavi;
the audit and execution tooling are supporting maintenance.

## Coverage and limits

The manual contains 63 Markdown pages, approximately 12,000 lines, and over
200 ordinary Julia fences. The whole manual was scanned for API/support/units
contradictions and repository-link targets. The detailed implementation review
focused on Jacobians, source and surface contracts, the ten Literate tutorials,
GPU boundaries, configuration examples, and related concepts/conventions.

This is **not** a claim that every equation, historical developer note, or all
999 registered docstrings have been independently verified. The earlier export
inventory and strict Documenter build were structural checks; they did not run
ordinary Julia fences. Some snippets are intentionally schematic or need
external datasets. A passing tutorial script also does not execute examples
that exist only in its prose/comments. External-paper equations, every numeric
line anchor, browser interactions, and Metal hardware remain separate audits.

## Findings corrected in this audit

| Finding | Consequence | Change |
|---|---|---|
| Hybrid AD and scattering tutorials advertised retired Mie `autodiff=true`; the public docstring did too | Documented calls throw | Use the hand-linearized `LinMode()` Mie path, state its native lognormal coordinates, correct the public docstring |
| Hybrid tutorial used undefined `model_cpu`; surface tutorial omitted imports; canopy tutorial accessed removed `model.params`; absorption tutorial collided with Makie's `Absorption` binding | Tutorials could render successfully but fail during execution | Fix references/imports and add an execution runner for all ten tutorials |
| `docs/test_examples.jl` still expected layout offsets from before surface-pressure column 1 | Five documentation assertions failed, outside the old docs build gate | Correct offsets; execute the contracts and tutorials in documentation CI |
| Source schema used an invalid keyword-only blackbody constructor, wrong spectral argument, and a thermal stub description | Broken examples and wrong physical interpretation | Runnable source example; distinguish a Planck-shaped collimated beam from atmospheric thermal emission; match polarization dimension |
| Source table multiplied the blackbody factor by an extra π; main-driver legacy SIF behavior was stale | Wrong source normalization/migration guidance | Match `F₀ = factor*scale*B`, document ignored legacy `RS_type.SIF₀` on `rt_run` |
| Main tutorial/plot labels called output radiance reflectance; Doppler 1/e half-width was called standard deviation | Unit/normalization errors in interpretation | Correct radiance labels and Gaussian width relations; distinguish measurement noise sigma from line width |
| Three Plotly panels used hand-built curves without clear synthetic labeling | Figures could be mistaken for computed HITRAN, BRDF, or canopy RT results | Label synthetic temperature, angular-surface, and red-edge illustrations in captions and plot titles |
| Surface Fourier docstring reversed the m=0/higher-m weights | Wrong extension recipe | Correct weights; add absolute `μ₀ F₀ BRDF` normalization and complete tangent checks to extension guide |
| MOM illustrative equal-μ transmission formula ignored different Stokes components sharing a direction | Misleading polarized-kernel example | State the Kronecker-delta limit and distinguish equal direction from equal matrix index |
| GEOS-Chem schema called nonexistent `IO.read_atmosphere` with an invalid source constructor | External-data example could never run | Use indexed `GeosChemSource`, `geoschem_to_dict`, and `read_atmos_profile_dict` |
| Architecture prose conflated optional CUDA Dual algebra with production analytic RT; AGENTS described spectral concatenation as vertical stacking | Misleading extension model | Separate batched Dual support from the plain-array RT tangent path; correct spectral concatenation |
| Mie deep-dive recommended the retired exclusion-angle knob | Tuning advice had no effect | Document its removal and use convergence/size-distribution-tail checks |
| One repository URL pointed to nonexistent `delta_m_truncation_lin.jl` | Broken equation-to-code link | Point to `delta_m_truncation.jl` |

The public Jacobian guide now describes the instrument chain rule, including
`(∂C/∂x)R` for fitted shifts/line shapes; matched-grid versus native-grid
precision tests; covariance whitening; and the local state-bias estimate.
It also makes the existing pressure-derivative limitation prominent rather
than leaving it only in the absorption schema.

## Highest-value remaining work for v2.2

### 1. Complete pressure-dependent absorption tangents

`src/CoreRT/tools/lin_model_from_parameters.jl` explicitly holds ordinary
line cross sections fixed when forming `τ̇_abs_psurf`. It differentiates the
bottom molecular column, but omits the pressure dependence of line shape/LUT
interpolation that a full forward rebuild includes. CIA and MT_CKD have
separate analytic pressure-scaling terms; their abundance response is also
incomplete.

The RT operator tangent can be correct while the upstream pressure derivative
is incomplete. This is a scientific limitation, not a documentation-only bug.
Prioritize the upstream cross-section derivative and centered finite
differences with real line absorption/ABSCO. Until then, describe the limited
coordinate explicitly and do not promise a complete pressure-retrieval
Jacobian. No new retrieval bias in ppm has been measured for this omission in
this audit.

### 2. Reject silently ignored source requests

A small CPU quickstart reproduction established:

- `SolarBeam(sza=75)` gives exactly the model-geometry result: source `sza` is ignored.
- `SolarBeam()+SolarBeam()` gives exactly one beam, and **half** the radiance
  obtained with the correctly summed irradiance. `extract_solar_F₀` uses the
  first beam only.
- Adding nonzero prescribed `SurfaceSIF` over an RPV surface produces exactly
  the same result as no SIF; the generic surface contribution is a no-op.

These behaviors are now documented, but not repaired. A useful v2.2 safeguard
is explicit rejection of incompatible SZA, multiple beam sources, and SIF on
unsupported surfaces, with tests across forward, linearized, split, and TOA
entry points. Supporting multiple geometries is a separate scientific feature;
for same-geometry beams, callers can sum irradiance into one beam today.

### 3. Make backend capability claims precise

`CoreRT._to_aa_arch` has CPU and CUDA mappings only, and the Metal extension
adds none. Production absorption paths using that adapter therefore lack
Metal support even though portable RT and Float32 Metal Mie kernels exist.
The GPU guide now states this limitation. A complete Metal absorption adapter
and real Apple-hardware forward/linearized validation are required before
claiming end-to-end gas-retrieval support there.

### 4. Keep instrument/state metadata at the retrieval boundary

A stable future extension should carry physical parameter identity, units,
native/retrieval coordinates, and band mapping together. Preserve the current
`ParameterKey`/`JacobianPlan` seam; avoid new hard-coded column offsets.
For any reusable instrument operator, make convolution, resampling, spectral
units, and shift/width derivatives part of one tested contract. Gate precision
on convolved radiance/Jacobians and converged retrievals, not only raw spectra
or a per-sample sigma threshold.

### 5. Separate computed figures and scientific references from illustrations

The immediate misleading labels are fixed. Replacing the three synthetic
panels with figures generated from the actual tutorial outputs would be more
valuable than adding more illustrative plots. Store fixture/configuration,
precision, source units, and result provenance with exported figures.

Numeric `file:line` anchors still need a systematic symbol-aware refresh:
links that point within a file can be stale even when a link checker passes.
Prefer symbol references for maintenance navigation and commit-pinned anchors
for archived numerical evidence. The full equation-to-paper review remains a
separate task; this audit did not revalidate every cited equation.

### 6. Keep larger integration work explicitly scoped

The GCHP sectional-AOD branch contains valuable reader, refractive-index,
and optical-property work, but its full sectional-to-RT bridge is incomplete
and it carries older RT changes. Port it as a bounded adapter with column-AOD,
phase, and RT validation rather than merging obsolete solver code wholesale.
Thermal and Raman Jacobians also require complete scientific propagation and
validation; adding an AD trait or exported name does not implement them.

## Validation record

Local logs: `/home/cfranken/vsmartmom-release-evidence/docs-audit-2026-09-08/`.
The original failed attempts are retained separately from corrected runs.
Final checks on Julia 1.12.6, CPU (CUDA hidden), with four Julia threads and
one BLAS thread:

| Check | Result |
|---|---|
| All ten Literate tutorials, including absorption/scattering/surface/canopy plots | **26/26 assertions, exit 0**, 3m50.8s; optional GIFs and CUDA sections excluded |
| Documentation contracts, including the actual source-schema code fence | **20/20 assertions, exit 0** |
| Strict Documenter/VitePress build, including generated plot assets | **Exit 0**, local deployment disabled |
| Export/docstring inventory | **475 exports, 999 registered docstrings, zero structural failures**; CUDA extension loaded, Metal unloaded |
| Repository `blob/main` link targets across all 63 manual pages | **No missing paths** after correction; heading/line semantics not certified |
| Documentation workflow YAML and `git diff --check` | Pass |

The initial tutorial run failed on an error in the new runner plus missing
imports/removed fields; corrected runs are retained with distinct filenames.
The final all-ten run is a single successful run. The first source probe also
needed an explicit RPV import; `source-probe-final.log` contains the complete
four observations.

[Inventory](documentation_audit_2026-09-08/inventory.json) and the
[source reproduction](documentation_audit_2026-09-08/source_contract_probe.jl)
are tracked. The reproduction runs from `test/` with the package/test environment.
Documentation CI now runs both execution checks before the strict build; hosted
CI has not run on these local commits yet.

No release was tagged, published, or registered. The prior numerical CPU/CUDA
record is in [release validation](release_validation_2026-09-07.md), including
its non-green full CPU run followed by passing affected-group corrections.
Those results are not relabeled as a new complete numerical run here.
