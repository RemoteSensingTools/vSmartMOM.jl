# Internal API Cleanup — Triage Brief

Created: 2026-04-28
Status: **non-urgent** — for a follow-up release after v2.0.0 registration. Not blocking docs work or Vitepress cutover.

## Background

Commit 42d9366 (warnonly retirement, missing-docs phase) chose to *document
all 48 missing-docstring exports* rather than de-export the internal-only
ones. That was the right risk call — zero API breakage, conservative for
registration. The 48 now live in
[`docs/src/pages/internal_api_coverage.md`](../src/pages/internal_api_coverage.md)
under a "Developer API Coverage" header that distinguishes them from the
canonical public API.

But every symbol in that file is now an implicit semver promise. Several
are genuinely internal — numerical accumulator helpers, internal precision-
policy machinery, dispatch traits — and shouldn't be public.

This brief organizes the 48 into four action buckets so a triage with the
source authors can be done in one sitting rather than per-symbol over time.

## Triage rules

Each symbol is assigned one of:

- **🟥 De-export** — internal implementation detail. Should be removed from
  the module's `export` line and demoted from `internal_api_coverage.md`.
  Stays accessible via fully-qualified name for anyone who really wants it.
- **🟧 Probably de-export, confirm with author** — looks internal but a
  contributor or external collaborator may already depend on it. Needs a
  one-question check before removing.
- **🟨 Keep public, move to main API reference** — genuinely user-facing
  but currently buried under the "developer coverage" header. Promote.
- **🟩 Keep public in developer coverage** — extension-author territory
  (traits, abstract types, linearized internals power users may want).

## The 48 by bucket

### 🟥 De-export (clear internals — 11 symbols)

| Module | Symbol | One-line justification |
|---|---|---|
| Scattering | `neumaier_add` | Kahan-Neumaier compensated-sum primitive. Numerical detail, not API. |
| Scattering | `neumaier_sum` | Same as above. |
| Scattering | `NeumaierAccum` | The accumulator type backing the above. |
| Scattering | `compute_Sl` | PCW-internal series sum. Implementation detail. |
| Scattering | `compute_mie_π_τ` | Mie integration kernel. Used inside Mie routines, not externally. |
| Scattering | `MiePrecisionPolicy` | Internal type controlling Float64 / DoubleSingle precision selection. |
| Scattering | `NativeFloat64` | Concrete `MiePrecisionPolicy` instance. |
| Scattering | `DSEmulated` | Concrete `MiePrecisionPolicy` instance. |
| Scattering | `DoubleSingle` | DoubleSingle arithmetic helper used inside PCW. |
| Scattering | `ComplexDS` | Complex DoubleSingle. Same scope as above. |
| Architectures | `@hascuda` | Internal CUDA-detection macro. |

Estimated commit: one per-module `export` cleanup, ~15 line edit.

### 🟧 Probably de-export, confirm with author (5 symbols)

| Module | Symbol | Why uncertain |
|---|---|---|
| InelasticScattering | `has_inelastic` | Dispatch trait. Useful for extension authors adding new `RS_type` subtypes — but only if there *are* extension authors. Default: keep until contributor lands a third-party Raman mode. |
| InelasticScattering | `uses_cabannes_phase` | Same logic as above. |
| InelasticScattering | `needs_interaction_workspace` | Same logic as above. |
| InelasticScattering | `needs_rayleigh_expansion` | Same logic as above. |
| InelasticScattering | `normalize_raman_weights!` | Helper for adding new Raman modes. Same gating question. |

The audit specifically recommended de-exporting these. The argument
*for* keeping them is "future extension authors need them"; the
argument *against* is "the package has no third-party Raman modes
today and the traits can be re-exported when needed." Recommendation:
de-export now; re-export when a contributor needs them.

### 🟨 Keep public, move to main API reference (12 symbols)

These are user-facing but currently sit under "Developer API Coverage."
They should be promoted into `pages/api_reference.md` (or the relevant
module Overview page) under proper sections.

| Module | Symbol | Where it should live |
|---|---|---|
| vSmartMOM | (module docstring) | Top of `api_reference.md` |
| CoreRT | (module docstring) | Top of CoreRT section in `api_reference.md` |
| IO | (module docstring) | Top of `IO/Schema.md` or IO section in `api_reference.md` |
| InelasticScattering | (module docstring) | `pages/Inelastic/Overview.md` |
| Aerosols | (module docstring) | `pages/Aerosols/Overview.md` |
| Absorption | `AbstractCrossSectionModel` | Absorption section, types subsection |
| Absorption | `AbstractComplexErrorFunction` | Absorption section, types subsection |
| Absorption | `save_interpolation_model` | Absorption section, methods |
| Absorption | `load_interpolation_model` | Absorption section, methods |
| Aerosols | `compute_optical_properties` | Aerosols overview — main user-facing entry point for both schemes |
| Architectures | `AbstractArchitecture` | New small `pages/Architectures.md` page or Surfaces section |
| Scattering | `compute_aerosol_optical_properties_gpu` | Scattering section — GPU variant of the documented public function |

Estimated commit: per-module promotion, ~30 line edits across
`api_reference.md` and module overview pages.

### 🟩 Keep public in developer coverage (20 symbols)

These belong in the "Developer API Coverage" section as currently scoped.
Extension authors and power users genuinely want them.

**CoreRT linearized internals:**
- `RTModelLin`, `OpticalPropertyJacobian`, `RawAerosolJacobian`
- `lin_added_layer_all_params!`, `delta_m_truncation_lin`, `delta_m_forward`
- `rt_run_ss` (the single-scatter solver — note: this is *also* documented in jacobians.md per the SF2023-II Eq. 32 discussion; could be promoted to main API ref instead)

**InelasticScattering computation internals:**
- `get_greek_raman_VS`, `get_n₀_n₁`, `getRamanSSProp!`
- `compute_effective_coefficents!`, `compute_σ_RoVibRaman_coeff!`
- `sol_RRS`, `sol_VS_0to1`, `sol_VS_1to0`, `sol_VS_0to1_plus`, `sol_VS_1to0_plus`

**IO extension point:**
- `IO.Formats.load_config` — the function third-party format adapters
  implement to register new file types.

These don't need to move; they're correctly scoped under "developer
coverage." If anything, the IO `load_config` could get a callout as the
extension hook for new formats — already covered in the
`extending/raman.md` or a new `extending/io_formats.md` page if you ever
want one.

## Suggested commit plan

Six small commits, in priority order:

1. **`api: trim Scattering numerical helpers from public exports`**
   Removes the 5 Scattering numerical helpers + 5 precision-policy types
   from `Scattering.jl`'s `export` line. Removes them from
   `internal_api_coverage.md`. ~15 lines.

2. **`api: drop Architectures.@hascuda from public exports`**
   One-line export removal + coverage page edit.

3. **`api: confirm Inelastic trait exports with source author`**
   Per the audit's recommendation. After confirmation, remove the 5
   trait exports.

4. **`api: promote module docstrings to main API reference`**
   Move `vSmartMOM`, `CoreRT`, `IO`, `InelasticScattering`, `Aerosols`
   module docstring blocks from `internal_api_coverage.md` to the
   appropriate places in the public docs.

5. **`api: promote Absorption + Aerosols + Architectures public symbols`**
   Move the 7 symbols listed in the 🟨 bucket to the main API
   reference / module overviews. Drop the duplicate entries from
   `internal_api_coverage.md`.

6. **`api: clarify Developer API Coverage header`**
   Update the page header to make explicit that the remaining 20
   symbols are extension-author / power-user territory and may change
   between minor versions.

After commit 6, `internal_api_coverage.md` should have ~20 entries,
all genuinely extension-or-power-user.

## When to do this

Not before:
- Vitepress cutover (would conflict with `api_reference.md` moves)
- v2.0.0 registration (let the registry tag a stable surface first)

Best window: as a v2.0.1 or v2.1.0 follow-up commit set, alongside
whatever other cleanup batches with it. Or earlier if a contributor
wants to land a third-party Raman mode and needs the trait exports
clarified.

## Hand-off

When ready to execute:
1. Run buckets 🟧 by source-author confirmation (one Slack / email).
2. Land commits 1–6 above. Each is small enough to revert independently.
3. Update [theory_references.md](theory_references.md) if any moved
   symbol's location changed.
4. Tag a minor release (v2.0.1 or v2.1.0) so the API change is
   semver-visible.

## Out of scope

- Renaming `vSmartMOM_Lin` for naming consistency with `RTModel` — that
  is a separate refactor flagged in CLAUDE.md memory.
- Changing the actual implementation of any of these symbols — this is
  a pure surface-area decision.
- Documenting symbols that aren't exported. `checkdocs = :exports` only
  cares about exports; non-exported symbols can have or lack docstrings
  freely.
