# v2.2.0 candidate validation — 2026-09-07

**Status: integration and local validation complete; final hosted CI and publication gates remain. Nothing has been published.**

Worktree: `vSmartMOM-release`; branch: `integration/release-candidate`.
This record supersedes the defect status in the September 6 audit of `b524a0d8`.
The [consolidation record](release_consolidation_2026-09-07.md) documents branch
scope and scientific attribution. Suniti Sanghavi is credited for the scientific
linearization/Raman work; integration and performance engineering retain helper
credit. Eight later Suniti workflow commits were ported with author metadata intact.

## Executed checks

| Check | Result and practical scope |
|---|---|
| Fresh Julia 1.10.12 resolution/install | Passed; no copied Manifest. Shipped quickstart forward RT finite. Shared depot/artifact cache used. |
| Fresh Julia 1.12.6 resolution/install | Passed; no copied Manifest. Shipped quickstart analytic Jacobians finite. Shared depot/artifact cache used. |
| Registered dependencies | Both resolve AtmosphericAbsorption 0.1.2, CanopyOptics 0.2.0, CUDA 5.11.3; only vSmartMOM is developed from the local candidate. |
| Complete CPU suite, optional Raman enabled | Full run: **12,429 passed, 1 failed, 1 errored, 14 skipped** in 25m10.9s. Both issues were test-harness defects (ambiguous CPU import and obsolete Taplo scope assertion). After test-only fixes, the affected groups passed **96/96**, reproducing the conflicting imports. No complete post-fix CPU rerun was performed; hosted CI must confirm the final tip. Combining the unaffected groups and corrected groups covers 12,441 passing assertions and 14 skips. |
| Complete CUDA runner | Final shipped runner passed **20,301/20,301** in 7m13.2s, including new no-pointer inverse dispatch coverage (`gpu-release.log`). |
| Cox–Munk exact/source/Jacobian contracts | 10/10 CPU and 10/10 CUDA. Checks exact pure-absorption radiance, whitecap/Lambertian normalization, incident-source scaling, dark source, cache replay, forward/linearized parity, and wind/optical-depth finite differences. |
| Scene schema | 74 YAML/TOML scenes plus 11 positive/negative contract cases passed. Includes public and local/data-dependent scene structure; no external spectroscopy is loaded by this check. |
| Julia numerics parsing | 34 additional Float32/Float64 assertions passed; invalid keys, aliases, convergence thresholds and choices reject explicitly. |
| Aerosol ingestion guard | Existing reader/RI unit tests: 59/59; new unsupported-optics guards: 10/10. Placeholder optics are removed. Heavy NetCDF ingestion is not tested without its fixture. |
| Portable retrieval workflows | 120 assertions passed; 14 explicit skips without external study inputs. |
| Extended retrieval workflows | 1,139/1,139 passed using original covariance/SIF inputs read-only, with temporary synthetic truth/measurement/prior/provenance products. No actual study inversion or HPC submission was executed. |
| API/docstring inventory | 475 exported bindings, 999 registered docstrings, zero mechanical failures, including loaded CUDA extension. Metal extension unavailable on this host. |
| Dispatch ambiguity inventory | 14 initial ambiguities reduced to zero on Julia 1.10 and 1.12; direct CUDA inverse dispatch 10/10 passed. Aqua ambiguity checking enabled. |
| README forward/Jacobian examples | 7/7 passed by executing the code blocks from README. The harness world-age/dimension-assumption errors were corrected separately. |
| Strict Documenter/VitePress | Final build passed with strict export/reference checks (Node 22.23.2, npm 12.0.2). |
| Docs dependencies | `npm audit`: zero vulnerabilities. Vite 6.4.3 override lies outside VitePress 1.6.4's declared Vite range; strict build plus development/preview HTTP smoke tests passed. These are not comprehensive browser interaction tests. |
| Relative documentation links | 243 README/manual relative file targets resolve. External URLs and heading anchors are outside this static check. |
| Taplo | Pinned 0.10.0 binary SHA-256 verified; lint and format checks pass. Parsed TOML equality confirmed for formatting-only changes. |

## Defects corrected during integration

- Cox–Munk direct glint had inconsistent absolute normalization, source scaling,
  attenuation and linearized treatment. The shared correction now includes
  exact-minus-Fourier BRDF, both atmospheric transits, incident irradiance,
  and wind/optical-depth tangents. Ocean baselines must be revalidated.
- CPU optimization previously extended NNlib's public batched multiplication.
  RT now owns its wrapper and preserves NNlib broadcasting/coexistence.
- The shipped GPU runner had broken includes; functional-hardware Jacobian
  compilation/scalar-index errors could be swallowed. Both are corrected.
- A broad DoubleSingle conversion and CUDA/CPU inverse overload intersection
  caused 14 dispatch ambiguities. Narrow conversion dispatch and an explicit
  CUDA no-pointer method resolve them.
- The experimental aerosol ingestion adapter returned fabricated Mie/number
  densities/asymmetry and mislabeled units. It now raises an explanatory error;
  production Mie/RT code is separate and remains available.
- Scene schema/parser drift, an invalid Draft-4 exclusiveMinimum declaration,
  four ambiguous YAML scientific-number spellings, stale azimuth banners,
  CUDA dependency/default claims, and native aerosol-size derivative labels
  were corrected. Legacy nonzero delta-BGE exclusion angles are retired.
- Ported campaign checks assumed old source paths and private test helpers.
  Public fixture tests and explicit producer-input paths retain checksum and
  canonical-checkout isolation requirements after relocation.

## Accuracy and performance evidence retained

The existing Jacobian precision, convolution, and complete-OE convergence
records remain under `jacobian_batched/`. Historical archived-versus-replayed
Suniti inversion timing shows approximately 12–20× speedups in those workloads,
not a fresh matched-worker old/new A/B. The strict 0.01-noise-sigma diagnostic
still failed 3/6 cases (worst 0.02076 sigma); it has not been relabeled as passing.
The separate full retrieval experiment reported approximately -0.0184 ppm
Float32/Float64 displacement with native grids and +0.00249 ppm with matched
grids, versus 0.314 ppm posterior sigma. These are retained scoped experiments,
not new end-to-end measurements from this release-validation pass.

## Remaining publication/platform gates

- The exact final candidate still needs the hosted Julia 1.10/1.11/1.12 ×
  Linux/macOS/Windows CI matrix. Earlier base-commit CI is historical evidence,
  not validation of this tip. Metal device execution was unavailable locally.
- JET reports **285 findings versus its historical baseline of 232** in this
  environment. It is advisory by default; the strict baseline is not met.
  Counts vary with Julia/JET versions, so the difference is not automatically
  53 new runtime defects. Do not describe a passing default suite as clean
  static analysis.
- General's latest registered version was 1.1.0 when checked. Existing public
  `v2.1.0` is preserved; 2.2.0 skips the registry's sequential-version
  guideline. RegistryCI supports a package-author approval override; otherwise
  the registration needs manual review. The GitHub latest-release endpoint returned v1.0.5,
  a separate publication channel from tags and the registry.
- TagBot is updated to issue-comment/manual triggers. No repository-level
  secrets were listed by the read-only audit; deployment-key availability is
  not established. Confirm a write-enabled DOCUMENTER_KEY and the documentation
  deployment path before registration/tagging. TagBot needs a deploy key to
  trigger docs on tag creation; commits changing workflows may require manual
  GitHub-release creation or an appropriately scoped publication credential.
- `gchp-io`, unsuccessful CUDA graph experiments, old divergent development
  histories, and private research products were deliberately not blindly merged.
  The GCHP sectional-AOD-to-full-RT bridge remains a separate workstream.
- Full OCO evaluator state-mapping tests need real spectroscopy/solar/prior
  inputs and are separate from the 1,139 workflow-contract assertions.

Primary references checked: [General versions](https://github.com/JuliaRegistries/General/blob/master/V/vSmartMOM/Versions.toml),
[RegistryCI version/author-approval guidance](https://juliaregistries.github.io/RegistryCI.jl/stable/guidelines/),
[TagBot setup](https://github.com/JuliaRegistries/TagBot#setup),
[TagBot deploy keys](https://github.com/JuliaRegistries/TagBot#ssh-deploy-keys).

## Reproduction and evidence

See [tools/README.md](../../tools/README.md) for portable commands and external
input switches. Full local logs are retained in
`/home/cfranken/vsmartmom-release-evidence/`. Compact [API, dependency-audit, and source-fingerprint inventories](release_evidence_2026-09-07/)
are checked in; no private input tables or campaign products are copied into
the release. The source manifest identifies the reviewed implementation tree;
the later evidence-only commit changes no numerical code.

The first CPU attempt exposed an obsolete exact-error-message assertion. A
subsequent run was superseded when dispatch fixes enabled Aqua ambiguity checks;
`cpu-release.log` records the completed full run;
`cpu-harness-corrections.log` records the subsequent 96/96 affected-group rerun.
The two failing/erroring groups are superseded by this focused rerun, not
represented as an originally successful full-suite process. An initial combined CUDA no-pointer overload passed Julia 1.12 static detection
but failed its direct GPU regression and Julia 1.10 static detection; separate
Float32/Float64 overloads fixed it. Earlier workflow harness
attempts failed on missing private inputs and relocated paths; only the final
completed run is reported as passing. These failures were corrected, not waived.
