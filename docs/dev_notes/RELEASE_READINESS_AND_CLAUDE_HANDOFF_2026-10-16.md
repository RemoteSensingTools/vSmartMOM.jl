# vSmartMOM release audit and implementation handoff — 2026-10-16

**To:** Suniti Sanghavi and her implementation agent, using Claude Opus 5.5 or Sol 6.1.

**Purpose:** consolidate the stable solver and supported workflows into a clean public release, with independently reviewed scientific corrections, reproducible tests, accurate documentation, and Julia General registration.

**Status:** audit and implementation plan, not an implementation record. The audit was read-only. This handoff and its evidence archive were subsequently saved at the user's request; no source fixes, branch merges, commits, pushes, GitHub issue changes, tags, registration, or publication were performed.

**Read this without the chat:** this document records the release decisions, findings, qualifications, source locations, numerical evidence, branch dispositions, proposed work, and acceptance criteria. Treat all repository/GitHub observations as a snapshot from October 16, 2026, and recheck state before implementation.

**Follow-on ocean work:** see the [atmosphere–ocean implementation plan and OceanOptics preparation](ATMOSPHERE_OCEAN_IMPLEMENTATION_PLAN_2026-10-16.md). It reviews the six supplied PDFs and subsequent Chowdhary benchmark ZIP in `~/oceanPapers/`, maps Fell's MOM interface equations, records the OceanOptics source inspection and reference-data inconsistencies, and defines commits, tests and review gates. The benchmark supplement is now acquired; its source-entry issues require a reviewed import. OceanOptics preparation can proceed independently; vSmartMOM integration waits for the consolidated release update. The accepted first coupled milestone is elastic polarized RT, with water Raman and fluorescence in later stages. This follow-on effort does not add ocean coupling to the current atmospheric release requirements.

## 1. Release verdict and instructions for resuming

### Verdict

**Do not release the current candidate unchanged.** Its existing CPU, CUDA, and executable-documentation checks pass, but independent checks found scientific correctness defects in Cox–Munk and incorrect HITRAN cache reuse. Legacy LUT construction, documentation, and test selection/reporting also need corrections.

**The core RT solver already has substantial independent scientific benchmark coverage.** The passing default suite includes 6SV1, Natraj polarized Rayleigh, Siewert IIA, and scalar/vector multilayer VLIDORT references. The targeted release findings do not invalidate that evidence. The remaining validation task is to map coverage to the advertised physics, derivatives, and execution paths, and extend it where the new findings expose gaps. It is not to establish external validation from scratch.

Preserve the integrated solver, selectively consolidate subsequent workflow changes, correct the confirmed defects, make supported derivative contracts explicit, and release the reviewed result as v2.2.0.

“Consolidated” means every relevant contribution has a documented disposition: incorporated, superseded, deliberately deferred, or retained as a study. It does not require merging every historical branch.

### Audited state

| Item | Audited value |
| --- | --- |
| Working repository | `/home/cfranken/code/gitHub/vSmartMOM-release` |
| Branch | `integration/release-candidate` |
| Commit | `d441061ba29cec333ec5c57228c5e0f408c86906` |
| Version in source | `2.2.0`, unpublished candidate |
| Recommended consolidation baseline | `integration/round5-jacobian-speedup`, `7acab57000dae259207a6760faae156cfa1734f6` |
| Latest General-registered version | `1.1.0` |
| Existing unregistered tag | `v2.1.0` |
| External comparison checkout | `/home/cfranken/code/gitHub/PyVLIDORT` |
| PyVLIDORT commit | `20d97e948f10d191edc69a8a9b1f3eeb91287227` |
| Fresh local runtime | Julia 1.12.6 |
| CUDA hardware used | NVIDIA L40S |
| Audit date | October 16, 2026 |

The Round 5 commit is **one commit after** the audited candidate; its parent is `d441061b`. It changes 14 workflow files, with 2,136 additions and four deletions. `src/`, `ext/`, `test/`, `config/`, and `.github/` are identical between these two commits.

The only preexisting untracked file in the release checkout was:

```text
docs/dev_notes/SUNITI_V2_2_JACOBIAN_INTEGRATION_MEMO_2026-09-09.md
```

Preserve it. Other worktrees contain unrelated dirty work and must not be reset or cleaned.

### Decisions already accepted by the user

- **HITRAN:** deprecate the legacy standalone line implementation with migration guidance; retain compatibility during 2.x and plan removal for 3.0.
- **Scope:** stable solver plus workflow compatibility; defer unvalidated research features.
- **Backends:** CPU and CUDA supported; Metal experimental.
- **Review:** independent review after each implementation batch, followed by review of corrections; Suniti approves scientific changes.
- **Derivatives:** support an explicit subset. Preserve documented native frozen-cross-section partials, but reject requests implying unsupported complete humidity/self-broadening/CIA/continuum abundance derivatives.
- **Compatibility:** preserve embedded-solar defaults and documented result/tuple compatibility. Keep opt-in capabilities opt-in.
- **Attribution:** scientific formulation, linearization, Raman physics, and retrieval methodology belong to Suniti Sanghavi. Preserve original authorship and distinguish supporting engineering. No co-author trailers.

Do not reopen these decisions merely because the implementation agent has a different default preference.

### Reading order and authority

1. [AGENTS.md](../../AGENTS.md)
2. [CLAUDE.md](../../CLAUDE.md)
3. [Conventions](../src/pages/conventions.md), before any cross-code comparison
4. This report
5. [September session handoff](SESSION_HANDOFF.md), [review corrections](review_corrections_2026-09-08.md), and the relevant theory/extension notes
6. The preserved [Round 5 integration memo](SUNITI_V2_2_JACOBIAN_INTEGRATION_MEMO_2026-09-09.md)

Older audit documents contain findings subsequently fixed. Do not treat every historical finding as outstanding.

This document is a handoff, not authorization to execute its full plan. When implementation is separately requested, use an isolated branch/worktree, preserve other users' work, and follow the review gates below. Merge, registration, and publication remain distinct release operations.

## 2. Bugs, scientific evidence, and limitations

### Finding summary

| ID | Priority | Finding | Evidence |
| --- | --- | --- | --- |
| B1 | P1 | Cox–Munk treats total mean-square slope as per-axis variance, doubling integrated total slope variance | Independent integration and VLIDORT source |
| B2 | P1 | Gaussian Smith shadowing has inconsistent normalization | Independent projected-Gaussian integral |
| B3 | P1 | Cox–Munk polarization geometry/reference planes and caller argument order are inconsistent | Compiled VLIDORT, independent meridian projection, vacuum RT |
| B4 | P1 | HITRAN cache reuse ignores isotopologues; artifact resolution bypasses coverage checks | Isolated offline reproductions |
| B5 | P2 | Legacy LUT builder ignores `vmr`; wavelength construction fails and has inconsistent conversion | Direct numerical reproduction |
| B6 | P2 | Documentation examples and docstrings describe incorrect calls, formulas, or types | Inspection and runtime reproduction |
| B7 | P2 | Documentation lockfile has four high-severity npm audit entries | Fresh npm audit |
| B8 | P2 | Combined-band driver applies the first selected surface to all selected bands | Inspected limitation; no new full numerical reproduction |

P1 denotes a release-blocking correctness issue. B8 is explicitly an inspected limitation rather than a newly reproduced end-to-end failure.

### B1 — Cox–Munk slope variance

**Source:** [coxmunk_surface.jl](../../src/CoreRT/Surfaces/coxmunk_surface.jl), especially `wind_to_sigma2` near line 22 and `cox_munk_pdf`.

The implementation combines:

```julia
σ² = 0.003 + 0.00512U
P = exp(-(zx^2 + zy^2) / (2σ²)) / (2πσ²)
```

The PDF treats `σ²` as the variance of **each horizontal slope component**, so:

```text
E[zx² + zy²] = 2σ²
```

The empirical wind relationship supplies total mean-square slope. With this PDF convention, use:

```text
σ²_axis = 0.5 × (0.003 + 0.00512U)
```

At `U = 5 m/s`:

| Quantity | Value |
| --- | ---: |
| Empirical total mean-square slope | 0.0286 |
| Numerically integrated total from current PDF | 0.05719999999999995 |
| Required per-axis variance | 0.0143 |

Independent quadrature reproduced the factor of two. The local VLIDORT GISS Cox–Munk source explicitly documents `0.5*(0.003 + 0.00512*WindSpeed)` as its per-axis input near line 1232.

**Impact:** glint width and amplitude change. This is not a universal radiance scaling error: the PDF changes differently across geometry. At the tested specular geometry, current intensity was half the standard-variance reference.

**Correction and acceptance:**

- Correct the per-axis variance and the wind derivative together: `dσ²_axis/dU = 0.00256`.
- Correct tutorial equations and derivative explanations.
- Test PDF normalization and integrated second moments independently of the implementation.
- Compare direct BRDF and vacuum RT intensity with approved external cases.
- Test multiple winds and both precisions; do not certify the correction from a derivative-vs-FD test alone.

### B2 — Smith shadowing normalization

**Source:** `smith_Λ` and `smith_Λ_dσ²` in [coxmunk_surface.jl](../../src/CoreRT/Surfaces/coxmunk_surface.jl).

The code defines:

```text
ν = cot(θ)/(√2 σ)
```

but evaluates:

```text
Λ = ½ [exp(-ν²)/(√(2π) ν) − erfc(ν)]
```

For that definition of `ν`, the Gaussian result uses `√π`:

```text
Λ = ½ [exp(-ν²)/(√π ν) − erfc(ν)]
```

The independent check integrated the projected Gaussian slope distribution:

```text
Λ = ∫ from c to ∞ of (s/c − 1) fσ(s) ds,  c = cot(θ)
```

Using the current `σ² = 0.0286` to isolate B2 from B1:

| μ | Current Λ | Independent integral |
| ---: | ---: | ---: |
| 0.1 | 0.1216765043 | 0.2864652175 |
| 0.2 | 0.0 | 0.0458159841 |
| 0.3 | 0.0 | 0.0065980654 |

The incorrect normalization makes the raw expression negative where shadowing is still nonzero. Clamping hides that discrepancy.

**Correction and acceptance:**

- Correct forward expression and analytic derivative together.
- Verify against quadrature over multiple variances and zenith angles.
- Check centered finite differences away from boundaries.
- Test bistatic shadowing and grazing behavior.
- Keep B1 and B2 acceptance separate so one correction cannot conceal the other.

### B3 — Polarized Cox–Munk geometry and argument order

**Source:** `coxmunk_geometry`, Mueller construction, and callers in the [surface implementation](../../src/CoreRT/Surfaces/coxmunk_surface.jl) and related surface builders.

Inspection found inconsistent outgoing-direction signs in plane calculations and rotations described in terms of scattering/facet planes where the solver needs the appropriate meridian references. Some production callers reverse incident/reflected inputs relative to the documented function contract.

Evidence came from:

1. The compiled, unmodified VLIDORT GISS surface kernel.
2. A separate geometric outgoing-meridian projection.
3. A complete vacuum-atmosphere TOA reproduction.

The tests deliberately matched variance to isolate B3 from B1, disabled whitecaps and shadowing, and used a real refractive index of `1.5`. The exactly representable index avoided an unrelated precision issue in the reference complex-index routine.

For SZA 30°, VZA 45°, relative azimuth 20°:

| Calculation | Q/I | U/I |
| --- | ---: | ---: |
| Current direct kernel, documented input order | 0.3987386 | 0.4364532 |
| VLIDORT after convention conversion | 0.5536500 | 0.2072573 |
| Independent meridian projection | 0.5536500 | 0.2072573 |

The production reversed input order produced another differing polarization result, approximately `Q/I=-0.42587018`, `U/I=-0.41002253`.

The vacuum forward solve returned:

```text
Current I,Q,U:
[0.05051108583996471, -0.021511165115848024, -0.02071068298021525]

Matched-variance reference:
[0.050511085839964734, 0.027965462046390544, 0.01046878916972477]

norm(Δ[Q,U])/I = 1.1577970833496574
```

Intensity agrees to numerical precision while polarization differs substantially. An intensity-only external comparison therefore cannot certify this surface.

**Convention handling used:**

- Follow [the repository conventions](../src/pages/conventions.md).
- Retain the raw relative-azimuth mapping used by the established kernel fixtures.
- Treat mission-wrapper conversion of solar bearing by 180° separately. Some Python wrappers form relative azimuth from view bearing minus solar propagation bearing.
- Do not impose a blanket 180° shift on every raw VLIDORT kernel call.
- Convert Stokes bases consistently. The full Mueller basis conversion is `D*M*D`, with `D=diag(1,-1,-1,-1)` for the convention conversion used here.
- The demonstrated unpolarized-source comparison exercises the first matrix column; it does not establish parity for every Mueller element.
- Divide the relevant reference π-scaled BRDF output by π.
- Keep BRDF normalization distinct from the `μ₀ F₀` factor in vacuum radiance.

**Correction and acceptance:**

- Define one incident/outgoing contract and apply it to all callers.
- Construct reference-plane rotations explicitly.
- Test off-principal-plane geometries, since principal-plane agreement can hide rotation errors.
- Validate full Mueller matrices, unpolarized first columns, vacuum RT, and wind tangents.
- Include normal incidence, specular, nonspecular, and grazing cases.
- Expand refractive-index and whitecap/shadow cases beyond the deliberately isolated reproduction.
- Require Suniti's geometry/convention approval.

### B4 — HITRAN cache identity and coverage

**Sources:** [hitran_api.jl](../../src/Artifacts/hitran_api.jl), especially the reuse predicate near lines 108–127, and [artifact_helper.jl](../../src/Artifacts/artifact_helper.jl), near lines 48–54.

Two wrong selections were reproduced offline:

1. Metadata describing isotopologue IDs `[7]` was accepted for a request for `[7,8]` over the same interval.
2. `artifact("CO2")` accepted an existing finite-band scratch file without checking that its coverage satisfied the intended dataset.

The fetch predicate checks numeric spectral bounds but ignores `iso_ids`. The artifact path largely checks existence. A narrower file can therefore enter the production model path as if it were the required dataset.

The isolated fixtures were synthetic. They establish cache-selection behavior, not scientific cross-section accuracy.

**Additional inspected weaknesses, not new interrupted-transfer reproductions:**

- Downloads open/truncate the final path before success is established.
- Nonempty content is accepted without fixed-width record validation.
- Some transport-error responses can retain a successful HTTP status.
- The requested “edition” string labels the cache; the URL does not itself select a historical HITRAN edition.

**Correction and acceptance:**

- Apply one coverage predicate to fetch and artifact resolution.
- Validate spectral coverage, requested isotopologue inclusion, metadata, and content identity.
- Treat missing/inconsistent provenance as insufficient for verified reuse.
- Download into temporary storage and validate before replacing a valid cache.
- Preserve previous valid data on failed replacement.
- Make metadata/content agreement checkable with hashes, including crash/failure cases.
- Invalidate package-owned parsed spectroscopy after managed replacement.
- Document that a cache label does not prove a historical database identity.
- Test all identity/failure paths offline using controlled responses.

An initial experimental probe had a Julia world-age problem and reported unsuccessful reproduction. It is superseded by `independent-final.log`, where `Base.invokelatest` reached the intended runtime fixture override. Do not cite the earlier false result as evidence of safety.

The parsed cache's documented requirement to call `clear_spectroscopy_cache!` after an external file replacement is a separate existing contract; it is not independently classified here as a new defect.

### B5 — Legacy interpolation-model construction

**Source:** [make_model_helpers.jl](../../src/Absorption/make_model_helpers.jl), `make_interpolation_model`.

**Ignored VMR:** the builder accepts `vmr` but omits it from `make_hitran_model`.

Reproduction used CO₂ isotopologue 1, lines from 6199–6201 cm⁻¹, evaluation grid `6199.5:0.25:6200.5`, pressure grid `600:100:800` hPa, and temperature grid `260:10:280` K. LUTs for `vmr=0` and `vmr=0.5` were identical; direct models changed by a relative norm of `0.017280872228360014`.

**Broken wavelength construction:** with `wavelength_flag=true`, the builder converts the input to wavenumber, passes the wavelength flag again during cross-section evaluation, and eventually stores a nonuniform vector in an `AbstractRange` field.

Observed exception:

```text
MethodError: Cannot convert Vector{Float64} to AbstractRange{<:Real}
```

**Correction and acceptance:**

- Pass `vmr` into the line model.
- Define exactly one spectral-coordinate conversion boundary.
- Preserve loading of existing serialized LUTs.
- For this release, reject unsupported wavelength-grid LUT construction early with a uniform-wavenumber migration recipe rather than casually changing serialization types.
- Preserve supported direct wavelength cross-section evaluation.
- Compare LUT values with direct evaluation at nodes; demonstrate nonzero VMR sensitivity where physically expected.

### B6 — Incorrect documentation and docstrings

Confirmed examples:

- [Absorption schema](../src/pages/IO/Schema/absorption.md), near line 130, calls `set_hitran_edition!(2020)`; the method requires a string and throws `MethodError`.
- The same paragraph calls 2020 the default; the actual default is `"artifact"`, resolving legacy HITRAN 2016 data.
- [Surface tutorial source](../src/pages/tutorials/Tutorial_Surfaces.jl) repeats the B1 variance and B2 Smith formulas and the old wind derivative.
- `_derive_m_max_bands_via_traits` in [model construction](../../src/CoreRT/tools/model_from_parameters.jl) documents `stream_l_cap=2*l_max-1`, whereas the corrected code uses `2*nstreams-1` and additional caps.
- The main `rt_run` docstring describes a legacy global Fourier bound instead of current derived per-band behavior.
- The Smith docstring names `smith_shadowing` and describes a shadow factor while attached to `smith_Λ`, which returns the shadowing integral.
- Three [Raman type docstrings](../../src/Inelastic/types.jl) describe the wrong types.

Runtime documentation lookup confirmed:

```text
RRS     → documentation for VS_0to1
VS_0to1 → documentation for VS_1to0
VS_1to0 → documentation for RRS
```

The onboarding guide also instructs callers to enter `test/` and run `julia --project=. -e 'using Pkg; Pkg.test()'`. This was reproduced failing because `test/Project.toml` is unnamed:

```text
The Project.toml of the package being tested must have a name and a UUID entry
```

A corrected invocation while remaining in `test/` is:

```bash
julia --project=.. -e 'using Pkg; Pkg.test()'
```

Other remaining drift includes legacy absorption examples despite the production AtmosphericAbsorption route, incomplete CIA/continuum descriptions, stale numeric source anchors, and accuracy statements whose scope exceeds the cited finite-resolution benchmark evidence.

**Acceptance:** correct the concrete examples/formulas/attachments, execute affected examples, and distinguish supported capabilities from experimental and historical material. Structural docstring presence is not semantic verification.

### B7 — Documentation dependencies

Fresh `npm audit --package-lock-only --omit=dev --json` reported four high-severity package entries:

| Entry | Relevant advisory information observed |
| --- | --- |
| `@vue/server-renderer` | [GHSA-g2v6-rqmx-r4w6](https://github.com/advisories/GHSA-g2v6-rqmx-r4w6); fixed range begins at 3.5.42 |
| `vue` | Shares the server-renderer advisory chain |
| `brace-expansion` | Multiple recursion/expansion advisories; the 5.x fixes require at least 5.0.12 for all observed entries |
| `source-map-js` | [GHSA-68fv-2mgg-jv7q](https://github.com/advisories/GHSA-68fv-2mgg-jv7q); fixed range begins at 1.2.2 |

These are four package entries, not proof of four independent exploitable flaws in the deployed static site. Recheck advisories when updating.

The pending brace-expansion PR to 5.0.9 does not resolve the later advisory ranges observed here.

**Acceptance:** update documentation dependencies and lockfile coherently, rerun the audit, rebuild and inspect representative pages/plots, and reconcile superseded bot PRs. Keep unrelated Julia major upgrades out of the release cleanup.

### B8 — Combined-band surface selection

**Source:** `_rt_run_column` in [rt_run.jl](../../src/CoreRT/rt_run.jl), near line 374.

The implementation selects:

```julia
brdf = get_surface(model, iBand[1])
```

and logs that multiple BRDFs are not implemented when multiple bands are selected. The configuration/model interfaces nevertheless describe surfaces per band.

This is an inspected correctness risk: a combined-band call with different surfaces can receive an unintended boundary condition. A log message is not a sufficient user contract.

**Release treatment:** reject unsupported combined-band surface combinations explicitly. Run independent bands separately; coupled inelastic calculations need a documented supported common-boundary contract. No new complete numerical reproduction of this limitation was performed.

### Scientific limitations: distinguish them from confirmed new bugs

| Area | Current meaning | Release treatment |
| --- | --- | --- |
| Gas abundance Jacobians | Native partials hold cross sections and other documented quantities fixed | Preserve and identify the partial derivative |
| H₂O self broadening | Forward absorption responds to moist mole fraction; native abundance tangent omits the full cross-section response | Reject requests claiming complete physical humidity derivatives |
| CIA/MT-CKD abundance | Forward terms exist; complete abundance tangents are absent | Guard unsupported retrieval requests |
| Surface pressure | Includes line/LUT response on the final-grid bottom-interface coordinate | Preserve; regridding is a separate chain rule |
| Raman/thermal Jacobians | Forward support does not establish complete tangent support | Explicitly reject unsupported combinations |
| Metal absorption | CPU/CUDA AA mappings exist; complete Metal production absorption does not | Experimental designation |
| Ross–Li grazing angles | Kernel zenith angles are capped at 80° | Document package regularization |
| BLAS thread setting | Per-model option changes persistent process-wide BLAS state | State application-level ownership |
| Optional CUDA inverse cache | Mutable scratch is not designed for concurrent identical-size calls | Keep opt-in; defer per-solve workspace redesign |
| Legacy profile bin averaging | Some gas VMR means are unweighted | Document/deprecate legacy path; any correction needs column-conservation tests |

A targeted H₂O experiment, with fixed dry column and broadener-changing abundance perturbations, found approximately **5.86% and 11.22%** norm discrepancies between native frozen-cross-section tangents and complete forward perturbations. The line-center ratios were approximately 1.0655 and 1.1320. This demonstrates different derivative definitions, not a measured retrieval bias.

Relevant code is in [atmo_prof_lin.jl](../../src/CoreRT/tools/atmo_prof_lin.jl), near lines 308–338, and [lin_model_from_parameters.jl](../../src/CoreRT/tools/lin_model_from_parameters.jl), near lines 339–379. The [Jacobian guide](../src/pages/jacobians.md) already states the frozen-cross-section coordinate. The missing step is an enforceable capability contract for requests that imply a different physical derivative.

The user chose the explicit subset rather than requiring every missing physical derivative before release.

### Existing independent core benchmarks and positive external controls

The initial overall assessment understated the existing independent benchmark coverage. After the user's correction, the actual default includes, reference suites, and complete CPU log were rechecked. The following are external scientific references, distinct from CPU/GPU parity or internal forward/tangent consistency:

| Reference | Coverage visible in source | Audit execution |
| --- | --- | --- |
| 6SV1 | Six wavelength/optical-depth/albedo configurations, three solar angles, three azimuths, sixteen viewing angles; 864 intensity comparisons summarized by six gates | Included in the passing complete CPU suite |
| Natraj 2009 Rayleigh | Optical depth 0.5, seven azimuths, sixteen viewing cosines, I/Q/U | Included in the passing complete CPU suite |
| Siewert 2000 Problem IIA | Aerosol slab, eleven viewing angles, three azimuths; available I/Q/U/V reference tables and convention conversions | Included for Float64 and Float32 CPU |
| VLIDORT solar_tester scalar | 23-layer molecular/aerosol atmosphere; selected Task 1 geometry with TOA-up/BOA-down reference values | Included for Float64 and Float32 CPU |
| VLIDORT solar_tester vector | 23-layer atmosphere with Problem III aerosol moments; selected Task 1 geometry, I/Q/U, TOA-up/BOA-down | Included for Float64 and Float32 CPU |
| Broader Natraj CDS tables | Seven optical depths, three albedos, seven solar cosines, sixteen viewing cosines, seven azimuths, I/Q/U in both directions | Manual benchmark/data available; full sweep not freshly run in this audit |

See [test_CoreRT.jl](../../test/test_CoreRT.jl), [the VLIDORT cases](../../test/vlidort_baseline/runtests.jl), [Natraj CDS inventory](../../test/benchmarks/natraj_cds/README.md), and [the full CDS driver](../../test/benchmarks/benchmark_natraj_cds.jl).

The default runner explicitly includes the VLIDORT suite near its end. A comment inside the standalone VLIDORT runner still says it is not included; that comment is stale and must not be used to infer execution status. The harness can additionally select CUDA when available, but the fresh complete CPU audit deliberately hid CUDA and therefore establishes CPU execution of these references.

The fresh Float64 log reports Siewert I relative errors around 5e-5–8e-5 and polarized maxima through approximately 3.3e-4 under that case's regularized metric. The multilayer vector case reports TOA I/Q/U maxima around 3.5e-4–4.9e-4 and BOA maxima through approximately 1.21e-3. These are results at the configured resolution, not universal accuracy limits. Keep precision, geometry selection, near-zero treatment, and numerical settings attached to each benchmark claim.

**Assessment correction:** describe vSmartMOM as having an independently benchmarked numerical core and remaining gaps in particular surface/spectroscopy paths, derivative contracts, and release integration. Do not describe it as generally lacking external scientific validation, or infer core immaturity from the Cox–Munk findings.

The newly executed external surface comparison additionally found agreement:

- RPV matched the corresponding three-parameter VLIDORT configuration with `ρc=ρ0`.
- RossThick and LiSparse matched at tested angles through 80°, approximately to 10⁻¹³ or better.
- The land comparison passed **18 assertions** across six geometries and three kernels.

At SZA 85°, VZA 70°, azimuth 60°, RossThick and LiSparse differed by approximately **25.9% and 26.4%** because of the package's angle cap. Local VLIDORT source does not support describing that cap as universal “MODIS/VLIDORT practice.”

These results validate particular kernel configurations. They do not establish universal surface parity or a newly run complete VLIDORT atmospheric campaign.

## 3. Streamlining tests without losing scientific coverage

### Inventory and interpretation

The inspected tree contains approximately:

- **164 Julia scripts and 44,013 lines under `test/`**, including benchmarks and utilities.
- **83 Julia scripts and 15,780 lines directly under `test/`.**
- Separate CPU, local hardware/data, workflow, tutorial, documentation-contract, and README-example runners.

These counts do not prove scientific redundancy. The accumulation mixes tests, benchmarks, fixtures, and studies. No byte-identical scene/profile fixtures were found in the checked `config/`, `test/test_parameters/`, `test/local/test_parameters/`, and `test/test_profiles/` files. Do not delete fixtures on a line-count argument.

Preserve the 6SV1, Natraj/CDS, Siewert, and LIDORT/VLIDORT golden/reference coverage as named scientific integration profiles. Simplification should remove duplicated setup/execution and clarify ownership, not discard independent references because many values are checked by a small number of assertions.

### Concrete findings and changes

| ID | Observation | Recommended change |
| --- | --- | --- |
| T1 | Tiled Jacobian file included twice when the CUDA flag is enabled | One runner owner |
| T2 | Shared helper names are redefined with different return contracts | One explicit support module |
| T3 | Some “forward” helpers invoke the linearized solver | Explicit names plus independent forward checks |
| T4 | Local GPU runner omits flag-gated CUDA Jacobian paths | One complete hardware profile |
| T5 | Shared imports, globals, and current directory leak between files | Isolated files/modules and local fixture paths |
| T6 | Legacy scripts allow NaN summaries or broad historical residuals | Repair/classify/archive before treating as accuracy gates |
| T7 | Local data runner catches broad errors and continues | Skip known missing resources; fail unexpected errors |
| T8 | Unit/integration/static/optional checks share one orchestration path | Explicit profiles with preserved default coverage |
| T9 | Test file presence is confused with execution | Active/optional/benchmark/historical inventory |

**T1: proved duplicate execution.** [test_jacobian_batched.jl](../../test/test_jacobian_batched.jl), near line 100, includes `test_jacobian_tiled.jl` when `VSMARTMOM_JACOBIAN_GPU_TEST=true`. [runtests.jl](../../test/runtests.jl), near lines 57–58, includes it again under the same flag.

**T2–T3: helper contracts.** [test_helpers.jl](../../test/test_helpers.jl), [test_jacobians_unit.jl](../../test/test_jacobians_unit.jl), and [test_hybrid_ad.jl](../../test/test_hybrid_ad.jl) define overlapping `run_lin_rt`, `run_fwd_only`, or `rel_errors` helpers. Relative-error helpers return different tuple/named-tuple shapes and different near-zero fallback behavior.

`run_fwd_only` obtains radiance through `LinMode()`. That can validate tangent-path consistency, but it is not an independent ordinary-forward check. `run_fwd_noRS` assigns an entire `ObserverRTResult` to a variable documented as a reflectance array. These are test-tooling contracts, not evidence that ordinary production RT is broken.

**T4: missing hardware profile coverage.** [The local GPU runner](../../test/local/gpu/runtests.jl) does not enable several root-level CUDA Jacobian arms. The additional audit run exercised them successfully; the suite structure should make their execution routine and explicit.

**T5: state and setup.** Root tests share imports and top-level helper definitions. Many paths are relative to `test/`. Some algorithm switches are correctly restored with `try/finally`; make that uniform. [Float32 checks](../../test/test_float32.jl) rebuild the same basic Float32/Float64 scenes multiple times; reuse immutable fixture preparation where safe, and copy mutable inputs.

**T6: historical tests are not automatically current gates.** [test_forward_lin.jl](../../test/test_forward_lin.jl) is referenced by the tutorial but excluded from the default runner. It contains broad 15–50% residual thresholds and `isnan` allowances with historical bug comments. Those comments are not a fresh measurement of current defects. Review the scientific purpose and repair or archive the cases before promoting them to required gates.

Matching-NaN options in comparisons should be paired with explicit finiteness checks; agreeing NaNs alone must never certify an accuracy test.

**T7: skipped setup can hide bugs.** [test/local/runtests.jl](../../test/local/runtests.jl) catches arbitrary parser/build errors as local setup limitations. Distinguish known absent data from unexpected `MethodError`, invalid configuration, or solver failures. An explicitly required external profile must not succeed after executing no relevant cases.

**T9: unowned or non-default scripts.** These root `test_*.jl` files are not directly in the default include list:

- [test_absorption_batched.jl](../../test/test_absorption_batched.jl): externally dependent script, activates another environment, and can exit successfully without testing; it also requires an AA API not guaranteed by the registered pin.
- [test_forward_lin.jl](../../test/test_forward_lin.jl): legacy derivative checks, above.
- [test_forward_raman.jl](../../test/test_forward_raman.jl): older configuration paths and data-dependent workflow; classify explicitly.
- [test_helpers.jl](../../test/test_helpers.jl): support, not a standalone suite.
- [test_performance.jl](../../test/test_performance.jl): identifies itself as a diagnostic benchmark, not a correctness test.

Do not move everything in `test/benchmarks/` indiscriminately: [test_CoreRT.jl](../../test/test_CoreRT.jl) consumes committed 6SV1/Natraj reference fixtures from that directory.

### Proposed profiles and runner contract

Use the existing numerical test project and a small explicit runner/catalog. Do not introduce a new testing framework during release cleanup.

| Profile | Contents | Required behavior |
| --- | --- | --- |
| `unit` | Small algebra, parsing, layouts, contracts, synthetic optics, portable CPU kernels | Offline after installation; focused files run independently |
| `integration` | Full model builds, Jacobians, observers, batch updates, external numerical fixtures, frozen regressions | Explicit reference provenance and tolerances |
| `cuda` | All CUDA paths, including flagged root tests and CPU comparison oracles | Fail if explicitly requested hardware is unavailable |
| `external` | ABSCO, covariance/SIF products, campaigns | Declared prerequisites; required execution cannot silently skip all |
| `quality` | Aqua and structured JET reporting | Distinguish package health from scientific accuracy |
| `docs` | Examples, contracts, tutorials, doctests | Execute published behavior |
| `benchmarks` | Timing/allocation studies | Report performance; separate from correctness acceptance |

Keep `Pkg.test()` as the canonical complete CPU entry point. Initially preserve default-suite coverage and the current nine-cell CPU matrix.

For each case, catalog its identity, profile, data requirements, backend/precision coverage, and reference type. Detect duplicate execution by case identity plus backend/precision. One file may provide distinct CPU/CUDA cases; that is different from executing the same tiled cases twice.

Use isolated modules/functions and one support module. Restore global switches in `finally`. Avoid concurrent in-process cases that mutate process-global algorithm flags, BLAS settings, or CUDA policy.

Only simplify the hosted matrix after demonstrating coverage equivalence. Expensive advisory analysis can then run once in a pinned environment while required scientific/platform coverage remains enforced. A warm-cache unit profile under about five minutes is a useful target, not a measured guarantee.

### Environment and dependency discipline

- Resolve independently for Julia 1.10, 1.11, and 1.12.
- Do not reuse a Julia 1.12 manifest as a Julia 1.10 lockfile. Earlier evidence already demonstrated incompatible resolved precompile dependencies.
- Bootstrap the checkout explicitly for focused tests; do not depend on the user's global environment.
- Keep site/plot dependencies in the docs environment.
- Remove test dependencies only after checking imports and profile ownership.
- Retain `AtmosphericAbsorption = "=0.1.2"` until its private prepared-line interface is replaced by an approved public seam.
- Do not advertise optional CUDA installation while CUDA remains a direct package dependency.
- Preserve Julia 1.10-compatible setup; do not adopt newer-only workspace behavior as an undocumented requirement.

The supported Pkg workflow merges package and test dependencies into a test environment: [Pkg test documentation](https://pkgdocs.julialang.org/v1.11/creating-packages/).

### JET and inference

Fresh JET reported **294 findings**, versus a historical baseline of 232. Current JET was 0.11.6; the old baseline used a different version. The difference is not 62 confirmed new runtime bugs.

- Save structured findings grouped by symbol/diagnostic, rather than printing only the first 20.
- Separate actionable undefined/error paths from macro/inference noise.
- Keep broad cross-version analysis advisory.
- Gate a small pinned set of representative calls with actionable diagnostics.
- Do not increase a count limit just to obtain green status.
- Separate float preservation, inference, allocations, and scientific accuracy in test names and claims.

Avoid a broad type-parameter rewrite. First measure one boundary with inference/allocation/runtime/compile-latency evidence. Abstract fields alone do not prove a practically important performance regression.

## 4. Documentation, docstrings, and source readability

### Current documentation execution

The manual contains **63 Markdown pages, 202 ordinary Julia fences, and zero jldoctest blocks**.

Existing execution is useful but incomplete: ten Literate tutorials, selected documentation contracts, selected README examples, and an export/docstring inventory. A strict rendered build does not execute all ordinary fences, and structural docstring attachment does not establish correct scientific content.

Before release:

1. Correct B6's calls, equations, defaults, and attachments.
2. Convert small deterministic public examples to doctests/executable examples.
3. Keep expensive demonstrations in explicit tutorial/contract runners.
4. Mark schematic snippets and external-data prerequisites.
5. Assert outputs, shapes, units, and parameter meaning rather than completion alone.
6. Generate tutorial Markdown from one maintained source.
7. Build generated docs in staging so verification does not rewrite tracked source pages.
8. Use symbol-aware navigation; preserve commit-pinned numeric anchors for archived evidence.

See [Documenter doctests](https://documenter.juliadocs.org/stable/man/doctests/).

### A clearer beginner-facing API

Organize the primary reference around:

```text
read_parameters
    → model_from_parameters
    → rt_run
    → named radiance results

JacobianPlan / LinMode
    → linearized construction
    → rt_run_lin
    → layout-aware Jacobians
```

The inventory's 478 exports span package modules; they are not 478 beginner operations. Preserve public compatibility while presenting stable, advanced, experimental, and deprecated reference tiers. Put implementation kernels in advanced documentation.

Important public docstrings should state signature/defaults, physical quantity/units, axes/ordering, mutation/ownership, return/unavailable-field behavior, supported physics/backends, derivative coordinate, held-fixed quantities, and one executable minimal example.

Clarify these boundaries:

- Absorption `(spectrum,layer)` versus contributor inputs `(layer,spectrum)`.
- Diffuse operators versus rectangular external-solar columns.
- Radiance, normalized radiance, and reflectance.
- Named `ObserverRTResult` fields and historical tuple iteration.
- Scalar `obs_alt=0` versus vector `[0]`.
- Fourier orders versus legacy counts.
- Physical, local optical, and retrieval-coordinate Jacobian layouts.

Preserve Unicode notation that maps to papers. Add a symbol glossary and local shape comments instead of romanizing the mathematics.

### Source organization

The model-construction file is approximately 1,261 lines, batch updating 1,136, the forward driver 1,252, and central types 1,830. Size alone is not a bug. The concern is mixed responsibilities and duplicated preparation across forward, tangent, and update builders.

Refactor at physical boundaries, one reviewed extraction at a time:

- Absorption preparation/evaluation.
- Molecular/Rayleigh preparation.
- One Fourier-moment solve.
- Surface interaction.
- Observer result packaging.

Preserve exact finite-thickness elemental expressions, product ordering, D-matrix symmetry, scattering-controlled doubling, and intentional wider intermediates. Require rebuild/update parity and finite differences for upstream extractions. Keep structural refactors separate from scientific corrections.

Move completed Phase A/B/C narratives to linked notes, retaining concise comments about present equations, stability, ownership, update order, and provenance. Keep the existing Concepts arc as the scientific narrative instead of duplicating it across onboarding files.

Larger Raman fixed/prepared/scratch-state redesign, carrier specialization, executor-level threading, and mutable CUDA workspace redesign are later maintenance work, not prerequisites for solving B1–B6.

### Already corrected historical findings

Do not redo or revert the September fixes:

- Strict standalone HITRAN records, including malformed fields and special isotopologue notation.
- `RTModel.propertynames` aliases.
- Standard tolerance keywords on `GreekCoefs.isapprox`.
- Broken legacy `noRS` Cabannes assignment.
- Shared analytic-surface scaffolding.
- Explicit float-type profile-reader support.
- Transactional batch updates and owned inputs.
- Batched algebra mutation/aliasing contracts.
- Pressure-dependent line/LUT surface-pressure response.
- Shared guards rejecting ignored source requests.

The newly confirmed wrong Raman type attachments still require correction despite earlier documentation improvements.

## 5. Branch consolidation and GitHub disposition

### Contribution ledger

| Source | Observation | Disposition |
| --- | --- | --- |
| Candidate `d441061b` | Audited integrated solver | Preserve corrected core |
| Round 5 `7acab570` | Direct child; 14 workflow changes | Start consolidated work here |
| `origin/suniti_multi_sensor`, `c26c85c32cc6654ba981f7f6bf710b9d2c168c19` | October study changes on older divergent core | Port selected deltas, not whole branch |
| Pure-gas predecessor `8814f305ce7ec8051a6b336a44fef4af32d88e52` | Bounded helper contribution | Review separately |
| Sanghavi active `uni_vSmartMOM`, local `f14371a8` | Dirty workflow changes and untracked campaigns | Preserve explicit snapshots |
| Sanghavi Round 5 worktree | Same integrated Round 5 lineage | Already represented |
| Sanghavi HWO worktree | Dirty Rayleigh/API/docs/parser work | Review bounded contribution; preserve newer parser |
| Cox–Munk scaling worktree | Relevant fixes mixed with older work | Use evidence; transplant reviewed corrections only |
| Jacobian performance `0540d37a` | Performance work already integrated with ancestry/author changes | Avoid duplicate ports |
| GCHP `c2c88a83` | Diagnostics/LUT changes; full bridge incomplete | Defer bounded adapter |
| Earlier absorption/Mie/unified lineages | Major work already incorporated | Preserve provenance; do not remerge |

The October branch has **135 candidate-only and 15 branch-only commits**, with **113 differing source/extension paths**. Raw divergence includes relocated or integrated work and does not equal missing functionality.

Its latest campaign commit changes 114 files with about 21,529 additions, without changing core source in that commit. Its hosted matrix is red due to an older finite-difference step issue already corrected in the candidate: do not restore the old `h=1e-5` setting from that lineage in place of the reviewed candidate setting.

Older externalAbsorption, caller-node Mie, Raman sanity, unified, and exoplanet contributions were already represented in candidate ancestry. The GCHP adapter and unvalidated study work need explicit deferred dispositions, not wholesale integration.

### Worktrees and dirty changes

Relevant Sanghavi paths:

```text
/home/sanghavi/code/github/uni_vSmartMOM
/home/sanghavi/code/github/uni_vSmartMOM_round5_jacobian
/home/sanghavi/code/github/run_worktrees/uni_vSmartMOM_suniti_HWO
```

The active checkout had five tracked dirty files, **79 additions and 10 deletions**, involving flexible state/prior validation, optional active-to-core output mapping, component-path arguments, and SIF-template configuration:

```text
RRS_XCO2/inversion/RetrievalOutput.jl
RRS_XCO2/inversion/RetrievalState.jl
RRS_XCO2/inversion/VSmartMOMForward.jl
RRS_XCO2/inversion/test_optimal_estimation.jl
RRS_XCO2/scripts/common.jl
```

There were also untracked campaign/source/product directories; some result paths were unreadable. The inventory is not a certification of every private/untracked product.

The HWO checkout had nine tracked changed files, **497 additions and 48 deletions**. It introduces `AbstractMolecularRayleighMethod` / `RefractiveIndexRayleigh`, expands molecular provenance, and changes parser/docs/tests. An isolated helper check passed **14 assertions**, not the full HWO workflow. New Ar/CO₂ scientific coverage still needs review of supported ranges and provenance; preserve the existing Earth defaults.

One proposed citation correction is DOI `10.1016/j.jqsrt.2019.05.031`, replacing an unrelated earlier reference. Verify the cited equation and applicability before accepting new scientific claims.

Preserve the candidate's strict HITRAN parser; it already handles `0/A/B` isotopologue notation. Do not overwrite it with a weaker older parser.

Other Christian worktrees:

- `/home/cfranken/code/gitHub/vSmartMOM.jl`: preexisting unfinished `.gitignore` conflict and staged hooks.
- `/home/cfranken/code/gitHub/vSmartMOM-coxmunk-scaling`: dirty Cox–Munk/source-scaling work.
- `/home/cfranken/code/gitHub/vSmartMOM-jacobian-perf`: performance provenance already incorporated.

Do not reset, clean, switch, or resolve these unrelated worktrees. Some Sanghavi index permissions prevent normal status operations; read-only logs/diffs were inspected using per-command safe-directory configuration, without global configuration changes.

### Workflow compatibility acceptance

Port to the relocated `sandbox/workflows/RRS_XCO2/` paths. Preserve:

- Ordered parameter identities and local-to-global maps, not just column counts.
- Physical/local/retrieval coordinates and their chain rules.
- Prior means, covariance transforms, and active-state mappings.
- Instrument convolution/resampling and fitted-instrument derivative terms.
- Surface basis transforms when extending the solve grid.
- Independent trial parameters and intentional read-only LUT sharing.
- Existing Round 5 fixed-SIF semantics and source-adding opt-in choices.

Do not turn the 14 skipped portable workflow cases into presumed campaign validation. Require the external-data profile for whichever workflow capabilities are claimed as release-supported.

### GitHub issues and PRs

Two open issues were observed:

- [#208: read_parameters export](https://github.com/RemoteSensingTools/vSmartMOM.jl/issues/208): candidate export coverage exists; close when users can obtain the fix publicly.
- [#176: individual gas absorption storage](https://github.com/RemoteSensingTools/vSmartMOM.jl/issues/176): enhancement, not fulfilled merely by gas Jacobian arrays; defer explicitly unless implemented.

Eleven open PRs:

| PRs | Treatment |
| --- | --- |
| #223 | Docstrings largely superseded by candidate; reconcile after integration |
| #224, #225 | Requested LogExpFunctions/Parameters compat already allowed |
| #228, #238 | DataInterpolations major upgrades; defer pending validation |
| #233 | TimerOutputs major upgrade; defer pending validation |
| #226, #232, #234, #235, #236 | Reconcile against a current coherent docs dependency update |

Current main protection requires one approval and selected strict checks, but does not dismiss stale approvals or require last-push approval. Code-owner review is enabled without a CODEOWNERS file. Required checks currently include Linux 1.11/1.12, Taplo, and docs build; API/workflow and full supported-platform gates are not all required.

Make the final release PR's review requirements effective: fresh independent approval after substantive corrections, explicit Suniti scientific approval, and required release checks. Do not mistake the current code-owner flag for an enforced science-review assignment.

## 6. Implementation sequence, compatibility, and review contract

### Ordered batches

| Batch | Work | Acceptance |
| --- | --- | --- |
| 0 | Preserve handoff/evidence; record branch/source disposition and dirty-work boundaries | Another agent can resume without the chat |
| 1 | Establish a Round 5-based consolidation branch and contribution ledger | No unexplained solver regressions or duplicate ports |
| 2 | Correct B1–B3, tangents, and convention docs | Independent BRDF/Mueller/vacuum references and focused CPU/CUDA checks |
| 3 | Correct B4–B5; introduce standalone line-engine deprecation/migration | Offline cache/failure tests, direct/LUT comparisons, old LUT loading |
| 4 | Enforce supported derivative and combined-band contracts | Unsupported physical requests reject; supported fixed-humidity workflows remain valid |
| 5 | Port stable October/dirty workflow changes with provenance | Layout identities, priors, instrument chain rules, output compatibility |
| 6 | Remove proven duplicate execution, centralize helpers, repair B6/onboarding | Focused files run independently; documented examples/commands execute |
| 7 | Update docs dependencies and finish packaging/release notes | Current audit, strict docs build, migration review |
| 8 | Independent final review and complete final-revision validation | Blockers resolved; evidence tied to the approved release commit |

Batch 0's handoff/evidence preservation is now done. The implementation batches remain proposed.

Keep HWO research additions and the GCHP bridge deferred unless a bounded contribution receives scientific review and its own acceptance evidence.

### HITRAN migration policy

The production RT pipeline already uses AtmosphericAbsorption.jl for ordinary line-by-line calculations. Keeping a second standalone line engine indefinitely invites differences in partition functions, broadening conventions, defects, and docs.

For v2.2:

- Deprecate legacy standalone line-engine entry points with actionable migration examples.
- Retain compatibility throughout 2.x, including existing serialized LUT loading.
- Keep CIA/MT-CKD/LUT responsibilities distinct from retirement of duplicated line-by-line code.
- State partition-function/data-provenance differences. The internal and production engines do not necessarily use identical spectroscopy conventions.
- Preserve current spectroscopy defaults unless a change is explicitly documented and scientifically approved.
- Plan line-engine removal for 3.0.

The exact AA pin protects use of its internal prepared-line interface for pressure tangents. Obtain a reviewed public upstream seam before relaxing the pin. A local development checkout still labelled 0.1.2 is not equivalent to the registered tree.

Do not introduce new default numerical behavior as part of test reorganization or docstring cleanup. Cox–Munk bug corrections should be identified clearly in release notes because users may observe changed results.

### Per-batch commit and review discipline

1. State finding IDs, intended behavior, and unchanged compatibility contracts.
2. Add meaningful reproductions for correctness defects before fixing them.
3. Separate scientific corrections, structural refactors, dependency changes, and documentation rewrites into coherent commits.
4. Run focused tests and affected integration/backend checks.
5. Obtain an independent review of the actual diff and evidence.
6. Address findings in a follow-up commit.
7. Have those corrections reviewed again.
8. Record validation against the resulting commit.

An author rereading their own patch does not satisfy independent review. The reviewer must inspect equations/contracts and whether the tests can detect the claimed defect.

Suniti must approve surface physics, polarization conventions, derivative coordinates/missing terms, Raman/molecular formulations, reference values/tolerances, and spectroscopy-default interpretation.

Commit rules:

- One coherent purpose per commit; avoid formatting unrelated scientific code.
- Preserve contribution authorship and record source commits/paths for ports.
- Do not redo historical attribution corrections or fabricate attribution.
- No co-author trailers.
- Exclude unrelated dirty files, local manifests, and generated campaign products.
- Explain every skip/tolerance change; never regenerate truth simply to match new output.
- Demonstrate finite-difference step convergence when changing perturbations.
- Restore mutated global state even when a test throws.

Run tests appropriate to each change, then the complete final-revision gates. Do not repeat full suites after report-only edits, or substitute many cheap mirrored assertions for an independent scientific check.

Each validation record must state commit, Julia/dependency versions, backend/hardware, thread settings, data provenance, command, exit status, counts, and skip reasons. A passing defect-reproduction assertion is evidence of the defect, not acceptance of current behavior.

### Final release acceptance

Before the main merge:

- Resolve B1–B5 and B6's concrete documentation defects.
- Resolve docs dependency advisories or record a specific reviewed disposition.
- Enforce the B8 and derivative capability boundaries.
- Account for every selected contribution and every deferred branch.
- Pass complete CPU validation for Julia 1.10, 1.11, and 1.12 on Linux, macOS, and Windows.
- Pass real CUDA validation including the separate flagged Jacobian paths.
- Run the frozen RRS gate and explicitly required scientific profiles.
- Pass Aqua, schemas, examples, tutorials, and strict documentation rendering.
- Test installation in a clean environment and explicitly check cold data/artifact acquisition.
- Publish support/migration tables with defaults and limitations.
- Obtain final independent review of the exact release revision.

Metal remains experimental; absence of Apple hardware does not block the selected CPU/CUDA release.

### General registration and publication

The package is already registered; General currently ends at 1.1.0: [registered versions](https://github.com/JuliaRegistries/General/blob/master/V/vSmartMOM/Versions.toml). This is a new-version registration, not first-time registration. The GitHub release list includes v1.1.0 from December 2024; GitHub's “latest” API designation is not itself a complete version inventory.

Keep intended version 2.2.0 unless maintainers explicitly change release identity. Do not reuse the existing v2.1.0 tag.

Registering 2.2.0 from General's current history can trigger the sequential-version check. Prepare migration/breaking-change notes and use RegistryCI's documented package-author approval mechanism if requested; do not fabricate intermediate releases. The `[merge approved]` mechanism addresses that specific check, not scientific approval or every registry requirement: [RegistryCI guidelines](https://juliaregistries.github.io/RegistryCI.jl/stable/guidelines/).

After an approved merge and separately authorized release operation:

1. Register the exact reviewed revision.
2. Ensure tag/release refer to the registered tree.
3. Publish matching versioned docs.
4. Close #208 after the fix is publicly obtainable.
5. Reconcile superseded PRs; keep #176 and deferred science explicitly tracked.

## 7. Evidence, reproduction, and audit limits

### Fresh validation results

| Check | Result | Qualification |
| --- | --- | --- |
| Complete default CPU suite | **12,793 passed; 16 broken/skipped; zero failures/errors** | Julia 1.12.6; 20m25.3s |
| Existing external core references | **6SV1, Natraj, Siewert, scalar/vector VLIDORT cases passed within that suite** | Broader manual CDS sweep was not freshly run |
| Local GPU runner | **20,323/20,323 passed** | Functional L40S; includes host checks |
| Additional explicit CUDA Jacobian profile | **13,256/13,256 passed** | Eight files; 7m48.2s; overlaps CPU coverage |
| Documentation contracts | **20/20 passed** | Selected executable contracts |
| Ten Literate tutorials | **26/26 passed** | 5m46.9s; optional sections excluded |
| Portable workflows | **120 passed; 14 data-related skips** | External covariance/SIF coverage incomplete |
| API inventory | **478 exports, 1,006 registered docstrings, zero structural failures** | Structural, not semantic certification |
| Schema checks | **74 scenes and 11 contract cases passed** | Dedicated Python requirements |
| Taplo and diff checks | Passed | Lint, format check, git diff check |
| Fresh temporary package environment | Resolution, loading, finite quickstart passed | Existing shared depot/artifact cache |
| Independent scalar land kernels | **18/18 passed** | Matched configurations through 80° |
| Cox–Munk independent probes | Confirmed B1–B3 | Reproduction checks demonstrate defects |
| Offline HITRAN probes | Confirmed B4 | Synthetic identity/coverage fixtures |
| Legacy LUT probe | Confirmed B5 | Sensitivity and exception reproduction |
| Isolated HWO molecular helper | **14 assertions passed** | Not complete workflow/scientific certification |

Do not sum assertion counts into a unique-test total; coverage overlaps. The legacy suite title “GPU / Metal” does not imply real Metal validation.

### Evidence archive

Original temporary artifacts were copied to:

```text
/home/cfranken/vsmartmom-release-evidence/release-audit-2026-10-16/
```

The archive has **18 log/inventory files**, totaling **479,052 bytes**, plus [EVIDENCE_MANIFEST.json](/home/cfranken/vsmartmom-release-evidence/release-audit-2026-10-16/EVIDENCE_MANIFEST.json) recording their SHA-256 checksums, source commit, runtime, and interpretation notes.

This is a same-host evidence directory, not a public repository URL. Transfer the bundle with the handoff when moving to another machine. All essential findings and acceptance criteria are also described here.

| File suffix after `vsmartmom-audit-20261016-` | Purpose |
| --- | --- |
| `cpu.log` | Complete default CPU run |
| `gpu.log` | Local GPU suite |
| `cuda-optin.log` | Additional flagged CUDA Jacobian coverage |
| `doc-contracts.log` | Executable documentation contracts |
| `tutorials.log` | Ten tutorials |
| `workflows.log` | Portable workflow contracts and skips |
| `api.log`, `api.json` | Structural API/docstring inventory |
| `fresh-install.log` | Fresh temporary environment and dependency identities |
| `independent-final.log` | Final variance/Smith and offline cache reproductions |
| `legacy-lut.log` | VMR and wavelength LUT failures |
| `vlidort-surface.log` | Surface intensity controls and polarization discrepancies |
| `vlidort-vacuum-final.log` | Vacuum polarization and independent projection |
| `vlidort-land-final.log` | Scalar land controls and angle-cap differences |
| `probes.log` | H₂O/documentation/Metal observations plus superseded cache experiment |
| Other non-final logs | Preserved earlier attempts; use final counterparts |

The archive includes earlier failed/superseded attempts for transparency. Original ad hoc probe scripts were not saved as standalone files. Their methods, inputs, expected discrepancies, and external-kernel interface are documented here; the future regression tests must recreate them.

Fresh package resolution selected:

```text
AtmosphericAbsorption 0.1.2
  tree: 3f9b48e459c504dda42304551fe4743200f5ad26
CUDA 5.11.3
KernelAbstractions 0.9.45
DataInterpolations 8.10.0
CanopyOptics 0.2.0
JET 0.11.6
Aqua 0.8.16
```

### Current commands

Run from `test/`, following the repository convention. The audited Julia binary was:

```text
/home/cfranken/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/bin/julia
```

CPU baseline:

```bash
env JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES='' \
  julia --startup-file=no --project=. runtests.jl
```

Explicit GPU suite:

```bash
env JULIA_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=0 \
  julia --startup-file=no --project=. \
  -e 'using CUDA; @assert CUDA.functional(); include("local/gpu/runtests.jl")'
```

Documentation execution:

```bash
julia --project=../docs ../docs/test_examples.jl
julia --project=../docs ../tools/check_tutorials.jl
```

Portable workflow contracts:

```bash
env CO2_COVARIANCE_FILE='' julia --project=. local/workflows/runtests.jl
```

The additional CUDA run used `JULIA_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1`, `CUDA_VISIBLE_DEVICES=0`, and `VSMARTMOM_JACOBIAN_GPU_TEST=true`, with:

```julia
using Test, CUDA, LinearAlgebra, Statistics, Distributions, Random
using vSmartMOM, vSmartMOM.CoreRT, vSmartMOM.Architectures, vSmartMOM.Scattering
@assert CUDA.functional()
@testset "Explicit CUDA Jacobian coverage" begin
    for file in (
        "test_jacobian_batched.jl",
        "test_elemental_boundary_jacobians.jl",
        "test_jacobian_vendor.jl",
        "test_jacobian_blocked.jl",
        "test_phase_interpolation.jl",
        "test_local_jacobian.jl",
        "test_z_jacobian_tables.jl",
        "test_lambertian_jacobian_batched.jl",
    )
        @testset "$file" begin
            include(file)
        end
    end
end
```

The nested batched file includes the tiled tests once; do not also include them separately.

For frozen RRS release validation, explicitly enable `PHASE1B_CPU=1`. The Wigner stress checks use `VSMARTMOM_FULL_TESTS=true`; assess these independently of the default CPU count.

For schema validation use Python 3.12 with [tools/requirements-checks.txt](../../tools/requirements-checks.txt); system Python 3.9 did not supply the needed environment. The older `checks-venv` executable was broken on this host, so it was not evidence of a package defect. Taplo used 0.10.0.

### Rebuilding the external surface oracle

The Python wrapper was not ready to run as a full installed model. The audit compiled unmodified local Fortran sources into temporary build artifacts instead.

Use the active `vsup/vbrdf/` sources, not `vbrdf_pre_29July21/`. Compile in this order with GNU Fortran; the audit compiler was GCC 11.5.0:

```text
vlidort_def/vlidort_pars.f90
vsup/vbrdf/vbrdf_sup_outputs_def.f90
vsup/vbrdf/vbrdf_sup_aux.f90
vsup/vbrdf/vbrdf_sup_kernels.f90
```

Flags were `-O2 -fPIC -shared -ffree-line-length-none`, with `-J` and `-I` pointing to a temporary module/build directory.

The existing temporary library is:

```text
/tmp/vsmartmom-vlidort-audit-20261016-1i7niqlp/libvbrdf_audit.so
```

Relevant exported symbol:

```text
__vbrdf_sup_kernels_m_MOD_gisscoxmunk_vfunction
```

The call takes maximum parameter count, parameter count, a Float64 parameter vector, Stokes-square count, incident cosine/sine, outgoing cosine/sine, azimuth/cosine/sine, and a 16-element Float64 output. Integer scalars use 32-bit Fortran integers. The GISS parameter vector uses per-axis variance, real refractive index, and shadow flag. Preserve the declared argument order in the source.

The 16 outputs are row-major relative to Julia reshape: use `Matrix(reshape(out,4,4)')`, then the appropriate π and Stokes-basis conversions. For the isolated reproduction use wind 5 m/s, matched variance, refractive index 1.5, shadow false, whitecaps false, and the geometries in B3.

Preserve source/version/copyright provenance. Do not copy the Fortran implementation into the Julia package as part of the fix. Normal package tests should use reviewed, provenance-labelled reference values; an optional oracle rebuild can remain external.

### Hosted evidence and limits

Historical success on the exact candidate:

- [CPU/quality run 34280992287](https://github.com/RemoteSensingTools/vSmartMOM.jl/actions/runs/34280992287), September 8, exact `d441061b`.
- [Documentation run 34280992189](https://github.com/RemoteSensingTools/vSmartMOM.jl/actions/runs/34280992189), September 8, same commit.
- Round 5 test/docs runs 34414483580 and 34414483625 succeeded on `7acab570`.
- [October study failure 37966232746](https://github.com/RemoteSensingTools/vSmartMOM.jl/actions/runs/37966232746), exact `c26c85c3`, includes the older FD issue.

Limits:

- Fresh complete local CPU testing used Julia 1.12.6; hosted supported-platform evidence is historical, not a fresh October matrix.
- No fresh local strict rendered docs build was run because the current build rewrites tracked tutorial Markdown. Tutorials/contracts were executed.
- Not every ordinary Julia fence, equation, numeric source anchor, or docstring was verified.
- No actual Metal hardware validation occurred.
- Fresh package installation reused depot/artifact caches; the default CO₂ artifact URL responded successfully, but that is not a cold-depot download certification.
- Optional CPU scientific gates, including frozen RRS, can skip by default.
- External covariance/SIF/ABSCO-dependent retrieval campaigns remain unvalidated by the portable subset.
- The compiled oracle covered surface kernels and a vacuum check, not a fresh full atmospheric VLIDORT campaign.
- Opt-in CUDA scratch concurrency was inspected, not a reproduced race.
- The broad historical source-quality audits are not fresh proof that every deferred refactor is necessary.

**Next implementation:** verify current state, preserve other worktrees, start from the reviewed Round 5 lineage, and take the Cox–Munk corrections through independent review. Keep the support boundaries and evidence qualifications intact throughout consolidation.
