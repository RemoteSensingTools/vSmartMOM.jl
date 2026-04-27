# vSmartMOM.jl Pre-Registration Hardening Plan

Created: 2026-04-24  
Branch context: `sanghavi-unified`  
Audience: maintainers, Codex, Claude, and other AI-assisted development sessions

## Purpose

This file is a working instruction document for turning the current
`sanghavi-unified` branch into a release-quality `main` candidate for a
registered Julia package.

The goal is not to rewrite vSmartMOM.jl. The goal is to make the current package
reliable, testable, easier to navigate, and honest in its documentation before it
is presented as the new mainline version.

The most important theme is to make physical modes explicit in the type system:

- source mode: direct solar / source function integration
- Raman mode: no Raman, RRS, RRS plus, VS variants
- architecture: CPU / GPU
- surface model
- aerosol scheme
- linearization mode

Where possible, replace scattered `if`, `isa`, and `typeof` checks with small
multiple-dispatch methods. This should reduce code size and make future
extensions less copy-paste heavy.

## Ground Rules

- Do not commit any `Manifest.toml` file. Manifests are intentionally local and
  ignored. Registration should rely on `Project.toml` compat bounds.
- Prefer small PR-sized changes with tests over broad mechanical rewrites.
- Keep numerical kernel rewrites conservative unless regression tests are strong.
- Do not remove reference data or examples unless the replacement path is clear.
- When touching user-facing input validation, use stable exceptions such as
  `ArgumentError` or `DomainError`, not `@assert`.
- Every example shown in README or docs should either run in CI or be marked as
  schematic.
- Every bug fix should include a regression test or a documented reason why a
  regression test is not practical.
- Avoid adding new abstractions unless they remove meaningful duplication or
  encode a real physical distinction.

## Current Audit Snapshot

These are the known issues from the current audit. Recheck them before editing,
because line numbers and symptoms may shift.

### Test Failures

`Pkg.test()` currently fails in `test/test_CoreRT.jl` before the later refactor
tests run.

Known failures:

- six 6SV1 checks around `test/test_CoreRT.jl`
- three Natraj checks around `test/test_CoreRT.jl`
- Natraj deltas observed during audit:
  - `I`: about `0.0213` versus tolerance `0.002`
  - `Q`: about `0.056` versus tolerance `0.008`
  - `U`: about `0.064` versus tolerance `0.008`

Required action:

- Determine whether the reference data or the refactored solver is wrong.
- If references are stale, regenerate them with a documented script.
- If solver behavior changed unintentionally, fix the solver and keep the old
  references.

### Manifests

Local manifests exist but are ignored:

- `Manifest.toml`
- `test/Manifest.toml`

They must stay untracked.

Useful checks:

```bash
git ls-files '*Manifest.toml'
git check-ignore -v Manifest.toml test/Manifest.toml docs/Manifest.toml
```

Acceptance criterion:

- `git ls-files '*Manifest.toml'` prints nothing.

### VS-Plus Workspace Bug

The current non-`noRS` workspace logic can allocate an `InteractionWorkspace`
that does not match what the VS-plus path reads.

Relevant files:

- `src/CoreRT/rt_run.jl`
- `src/CoreRT/CoreKernel/interaction_inelastic.jl`

Observed problem:

- VS-plus code reads fields such as `tmpieJ...` / `tmpieR...`.
- The workspace currently has a different field set.

Required action:

- Either gate the existing workspace to only the modes it actually supports, or
  port VS-plus to the staged workspace correctly.
- Add a minimal VS-plus runtime smoke test. Parse-only tests are not enough.

Recommended implementation direction:

```julia
interaction_workspace(::noRS, args...) = nothing
interaction_workspace(::noRS_plus, args...) = nothing
interaction_workspace(::RRS, args...) = InteractionWorkspace(args...)
interaction_workspace(::RRS_plus, args...) = InteractionWorkspace(args...)

# Until the VS-plus path is fully ported:
interaction_workspace(::VS_0to1_plus, args...) = nothing
interaction_workspace(::VS_1to0_plus, args...) = nothing
```

Then update the call sites to dispatch on Raman mode instead of relying on
`!(RS_type isa noRS)`.

### Export Typo

Check `src/Inelastic/InelasticScattering.jl`.

Known issue:

- `sol_VS_0to1_plus` appears to be exported twice.
- `sol_VS_1to0_plus` is likely missing from exports.

Acceptance criterion:

- Public exports match the implemented public functions.
- Add a tiny import/export test if possible.

### SolarModel Path Bug

Calling:

```bash
julia --project=. -e 'using vSmartMOM.SolarModel; default_solar_transmission([13000.0])'
```

has produced an `UndefVarError` involving `RadiativeTransfer`.

Required action:

- Replace stale module path assumptions with `pathof(vSmartMOM)`, `pkgdir`, an
  artifact, or a scratch/data path that is valid in this package.

Acceptance criterion:

- The command above runs from a clean project environment.

### RAMI Test Path Bug

`test/test_rami_smoke.jl` uses a cwd-sensitive path like:

```julia
"../sandbox/rami/RamiNoGas.yaml"
```

Required action:

- Replace with a package-root-relative path, for example:

```julia
joinpath(pkgdir(vSmartMOM), "sandbox", "rami", "RamiNoGas.yaml")
```

Acceptance criterion:

- The RAMI smoke test passes from the repository root and from `test/`.

### CI And Compat

Current package compat allows Julia `1.9`, `1.10`, `1.11`, and `1.12`.

Project decision:

- It is acceptable to require Julia `1.12`.

Required action:

- Update `Project.toml` compat to Julia `1.12`.
- Update CI to test Julia `1.12`.
- Process pending dependency compat updates.
- Keep compat bounds accurate and not unnecessarily broad.

### Aqua And JET

Project decision:

- Add Aqua.jl and JET.jl as test tools.

Recommended approach:

- Add both to `test/Project.toml`.
- Add a focused `test/test_quality.jl`.
- Make Aqua blocking after obvious cleanup.
- Make JET targeted at first. Do not require whole-package JET cleanliness until
  dynamic and optional-code paths are better controlled.

Possible structure:

```julia
using Test
using Aqua
using JET
using vSmartMOM

@testset "Aqua" begin
    Aqua.test_all(vSmartMOM)
end

@testset "JET smoke" begin
    # Start with stable, concrete entry points.
    # Add more report_call/report_package checks as false positives are fixed.
end
```

### Docs Truth Gap

Known issues:

- `docs/make.jl` currently generates tutorials with Literate but does not appear
  to guarantee example execution.
- `makedocs(... warnonly=true)` lets documentation warnings pass.
- `docs/src/index.md` contains placeholder examples such as
  `parameters_from_yaml("path/to/your/params.yaml")`.
- Docs still refer to Julia `1.9+`.

Required action:

- Make docs examples executable or clearly schematic.
- Fail docs CI on broken examples once the examples are cleaned up.
- Update docs to Julia `1.12`.
- Prefer a small number of reliable tutorials over many stale ones.

## Time-Boxed Plan

## Week 1: Stabilize And Rebaseline

### 1. Establish Clean Baseline

Run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Also run syntax parsing over source, tests, docs examples, and examples:

```bash
rg --files -g '*.jl' src test docs/src/pages/tutorials examples |
  julia --project=. -e '
      failed = String[]
      for file in eachline(stdin)
          try
              Meta.parseall(read(file, String))
          catch err
              push!(failed, file)
              @warn "Parse failure" file err
          end
      end
      isempty(failed) || error("Parse failures: " * join(failed, ", "))
  '
```

Acceptance criterion:

- Parse check passes.
- `Pkg.test()` failures are understood and categorized as either real bugs or
  stale references.

### 2. Fix Release Blockers

Fix in this order:

1. CoreRT reference failures.
2. VS-plus workspace bug.
3. duplicate/missing VS-plus export.
4. SolarModel path bug.
5. RAMI cwd-sensitive test.
6. stale manifest confusion in local workflows.

Acceptance criterion:

- `Pkg.test()` reaches all normal test files.
- Any remaining failures are documented with specific owner/next action.

### 3. Move To Julia 1.12

Files to inspect:

- `Project.toml`
- `test/Project.toml`
- `.github/workflows/*.yml`
- `docs/Project.toml`
- docs and README setup instructions

Acceptance criterion:

- CI uses Julia `1.12`.
- package compat says Julia `1.12`.
- docs say Julia `1.12`.
- dependency compat PRs have been addressed or explicitly deferred.

### 4. Add Aqua And JET

Recommended staging:

1. Add dependencies to `test/Project.toml`.
2. Add `test/test_quality.jl`.
3. Include it from `test/runtests.jl`.
4. Initially use targeted JET checks, not broad whole-package analysis.

Acceptance criterion:

- Aqua and selected JET checks run in CI.
- Any suppressed Aqua/JET checks have comments explaining why.

## Week 2: Dispatch Refactors With High Payoff

The week 2 goal is to remove branch-heavy glue code without changing numerical
intent.

### 5. Replace `SFI::Bool` With Source-Mode Dispatch

Current smell:

- `SFI::Bool` is passed through many RT paths.
- The code already has source-mode marker types such as `DNI` and `SFI`.

Preferred direction:

```julia
abstract type AbstractSourceMode end

source_mode(::Val{true}) = SFI()
source_mode(::Val{false}) = DNI()

apply_source_terms!(::SFI, args...) = ...
apply_source_terms!(::DNI, args...) = nothing

postprocess_source!(::SFI, args...) = ...
postprocess_source!(::DNI, args...) = nothing
```

If `DNI` and `SFI` already subtype a local abstract source type, reuse that
instead of introducing another abstraction.

Migration plan:

1. Introduce source-mode helper methods.
2. Convert one high-level RT path.
3. Keep compatibility wrappers for old Boolean call sites.
4. Remove the Boolean branches after tests are green.

Acceptance criterion:

- Repeated `if SFI` branches are reduced.
- Existing public behavior is unchanged.
- New tests cover both source modes if both are supported.

### 6. Dispatch Raman/No-Raman Behavior

Current smell:

- Scattered checks such as `RS_type isa noRS`, `RS_type isa noRS_plus`, and
  `typeof(RS_type) <: noRS`.

Preferred direction:

```julia
has_inelastic(::noRS) = Val(false)
has_inelastic(::noRS_plus) = Val(false)
has_inelastic(::RRS) = Val(true)
has_inelastic(::RRS_plus) = Val(true)

rayleigh_greek_source(::noRS, model) = model.greek_rayleigh
rayleigh_greek_source(::noRS_plus, model) = model.greek_rayleigh
rayleigh_greek_source(::RRS, model) = model.greek_cabannes
rayleigh_greek_source(::RRS_plus, model) = model.greek_cabannes

inelastic_workspace(::noRS, args...) = nothing
inelastic_workspace(::noRS_plus, args...) = nothing
inelastic_workspace(::RRS, args...) = InteractionWorkspace(args...)
inelastic_workspace(::RRS_plus, args...) = InteractionWorkspace(args...)
```

Use explicit VS methods. Do not allow VS-plus to accidentally reuse an RRS-only
workspace shape.

Acceptance criterion:

- No-Raman, RRS, RRS-plus, and VS-plus code paths choose workspace behavior by
  dispatch.
- At least one smoke test exercises the VS-plus runtime path.

### 7. Clean IO Dispatch And Validation

Current smell:

- `load_config(src::IOSource)` branches on concrete source type.
- YAML validation uses many `@assert`s.

Preferred direction:

```julia
load_config(src::FileSource) = ...
load_config(src::DictSource) = ...
load_config(src::IOSource) =
    throw(ArgumentError("Unsupported IO source: $(typeof(src))"))
```

Validation helper direction:

```julia
require_key(config, key, context)
require_type(value, ::Type{T}, context) where {T}
require_option(value, allowed, context)
```

Acceptance criterion:

- Invalid YAML produces clear, stable errors.
- IO source handling uses dispatch.
- README/docs examples that involve YAML match the actual schema.

### 8. Surface Parsing Cleanup

Current smell:

- YAML surface parsing uses string maps, closures, and special-case branches.
- Surface implementations repeat scalar/vector normalization logic.

Preferred direction:

Keep strings at the YAML boundary, then dispatch on surface constructor type:

```julia
surface_type_from_name("Lambertian") = LambertianSurfaceScalar
surface_type_from_name("CoxMunk") = CoxMunkSurface

parse_surface(::Type{LambertianSurfaceScalar}, args, ::Type{FT}) where {FT} = ...
parse_surface(::Type{CoxMunkSurface}, args, ::Type{FT}) where {FT} = ...
```

Normalize repeated scalar/vector logic:

```julia
grid_values(x::Number, grid, ::Type{FT}) where {FT} = fill(FT(x), length(grid))
grid_values(x::AbstractVector, grid, ::Type{FT}) where {FT} = FT.(x)
```

Acceptance criterion:

- Surface YAML errors identify the bad surface and bad field.
- Surface constructors become easier to document.
- No behavior change for valid existing configs.

## Week 3: Structural Cleanup And Documentation

Week 3 is for improvements that are valuable before the package becomes more
visible, but only after the core tests are green.

### 9. Consolidate Inelastic Elemental Drivers

Current smell:

- `elemental_inelastic.jl` and `elemental_inelastic_plus.jl` duplicate driver
  structure.

Do not start with a physics rewrite. Start by extracting transition metadata:

```julia
inelastic_transitions(::RRS)
inelastic_transitions(::RRS_plus)
inelastic_transitions(::VS_0to1)
inelastic_transitions(::VS_1to0)
inelastic_transitions(::VS_0to1_plus)
inelastic_transitions(::VS_1to0_plus)
```

Then move generic driver loops to shared helpers.

Acceptance criterion:

- Duplicate driver code is reduced.
- Mode-specific physics remains readable.
- Tests cover every supported Raman/VS mode that the public API claims to
  support.

### 10. Clean Dead Or Development-Only Files

Candidates to inspect:

- empty surface files
- unincluded stellar type files
- kernel scratch/test files under `src`
- download scripts that belong under `dev/`, `sandbox/`, or artifact tooling

Required action:

- For each candidate, decide one of:
  - keep as runtime package code
  - move to `dev/`
  - move to `sandbox/`
  - remove

Acceptance criterion:

- `src/` contains runtime package code, not scratch experiments.
- Any moved file has a discoverable new location.

### 11. Repository Data And Artifact Hygiene

Known concern:

- The repository contains large sandbox/reference data files and generated
  outputs.

Required action:

- Classify tracked data as one of:
  - required package runtime data
  - small test fixture
  - large reference artifact
  - generated output that should not be tracked
  - sandbox-only input

Acceptance criterion:

- Large data needed for tests/docs is either justified or moved toward Julia
  Artifacts/Scratch/external fixtures.
- Generated logs and local outputs are untracked.

### 12. Public API Review

Required action:

- Review all exported names.
- Remove accidental exports where possible.
- Add missing exports where intended.
- Document the public workflow rather than every internal helper.

Suggested docs structure:

- Quick start from YAML.
- Forward RT.
- Polarized RT.
- Raman / VS mode.
- Linearized RT / Jacobians.
- RAMI-style scenario.
- Output structure and interpretation.

Acceptance criterion:

- A first-time user can run one documented example without knowing internal
  source layout.
- Public exports match the documented API.

### 13. Performance And Allocation Pass

Scope:

- Focus on obvious allocation hotspots and workspace construction.
- Avoid broad numerical kernel rewrites until regression tests are stronger.

Recommended targets:

- inelastic workspace allocation
- repeated temporary arrays in Raman interaction/doubling
- GPU staging and synchronization points
- avoid accidental CPU scalar indexing in GPU paths

Acceptance criterion:

- Add benchmark scripts for representative small and medium cases.
- Benchmarks are not part of normal tests unless they are lightweight.
- Performance changes have numerical regression coverage.

## README "In Progress" Suggestions

Add a concise section to the main README so users know what is actively being
worked on. Suggested bullets:

- Julia 1.12 modernization and dependency compat cleanup.
- Aqua.jl and JET.jl quality checks in CI.
- Runnable documentation examples and refreshed tutorials.
- Raman, rotational Raman, and vibrational Raman regression validation.
- GPU memory/workspace cleanup for inelastic RT paths.
- Surface model cleanup and clearer YAML schema.
- Aerosol and GEOS-Chem integration hardening.
- Artifact-based handling of large reference data and benchmark fixtures.

Keep this section honest. If an item is not actively being worked on, remove it
or move it to a roadmap document.

## Recommended AI Working Protocol

When an AI assistant picks up this plan:

1. Start with `git status --short`.
2. Read this file and any linked dev note relevant to the task.
3. Re-run only the checks needed for the task unless establishing a new baseline.
4. Make the smallest coherent change.
5. Add or update tests.
6. Run focused tests.
7. Update this file if the plan changes materially.
8. Summarize:
   - files changed
   - tests run
   - remaining risks
   - next recommended task

Do not silently rewrite unrelated code. If a task reveals a larger architectural
problem, document it here and keep the immediate patch focused.

## Commands To Keep Handy

Clean status:

```bash
git status --short
```

Tracked manifests:

```bash
git ls-files '*Manifest.toml'
```

Ignored manifest check:

```bash
git check-ignore -v Manifest.toml test/Manifest.toml docs/Manifest.toml
```

Full tests:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

Docs:

```bash
julia --project=docs docs/make.jl
```

One targeted test file:

```bash
julia --project=. -e 'using Test; include("test/test_file_name.jl")'
```

If a test depends on `test/Project.toml`, prefer:

```bash
julia --project=test -e 'using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.test()'
```

Search for branch-heavy dispatch candidates:

```bash
rg '\bisa\b|typeof\(|::Bool|if .*SFI|if .*RS_type|elseif .*RS_type' src test
```

Find largest Julia files:

```bash
rg --files -g '*.jl' src | xargs wc -l | sort -nr | head -30
```

Find potentially accidental tracked generated files:

```bash
git ls-files '*.log' 'nohup.out' '*.jld2' '*.jld' '*.nc' '*.nc4' '*.png'
```

## What Not To Do Before First Registered Main

- Do not attempt a full unification of forward, linearized, multisensor, and
  single-scatter kernels unless the regression suite is much stronger.
- Do not make whole-package JET strict on day one if it creates excessive noise.
- Do not remove large data files without preserving reproducibility.
- Do not broaden public API just to make internal refactors convenient.
- Do not replace readable physics notation with low-level mutation everywhere
  unless profiling shows a real cost.

## Definition Of Ready For New Main

The branch is ready to become the new main candidate when:

- Julia `1.12` is the tested and documented version.
- `Pkg.test()` passes in CI.
- Aqua passes, or every exception is explicit and justified.
- JET has at least targeted checks on stable public entry points.
- README examples run or are explicitly schematic.
- Docs build without hiding important warnings.
- Manifests are untracked.
- VS-plus behavior is either tested and working or clearly marked unsupported.
- Public exports match documented public API.
- Large data and generated files have a documented policy.
- The README has an honest "In Progress" section.
