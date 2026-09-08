# v2.2 code-quality follow-up — 2026-09-08

This pass applies the bounded, high-confidence items from the
[source-quality audit](source_quality_audit_2026-09-08.md). Its goal is a code
path a scientist can read in physical order, an engineer can maintain without
hidden ownership rules, and a future agent can extend without inferring
unstated contracts.

The pass deliberately preserves the numerical kernels, equation ordering,
Unicode scientific notation and public compatibility behavior. Idiomatic Julia
here means using multiple dispatch at genuine capability boundaries, standard
interfaces for package-owned types, concrete typed construction where useful,
and explicit mutation. It does not mean turning every conditional into a type,
parameterizing every carrier without measurement, or applying a cosmetic
formatter across sensitive scientific code.

## Completed changes

### Transactional batch updates

`BatchContext` owns reusable trial optical-depth arrays. `update_model!`
prepares the candidate profile, geometry, sources, Rayleigh state, aerosols and
every absorption band before copying any candidate into the live model. Source
shape validation uses dispatch for ordinary sources, thermal emission and
`SourceSet`; the same methods protect both the trial and commit paths.

Input temperature, pressure, humidity and VMR arrays are copied into owned
state. Failed updates therefore cannot mutate caller inputs through aliasing.
The successful commit preserves the model object and its leaf-array identities,
so callers may safely retain references to preallocated solver storage.

The regression intentionally fails in the second of two absorption bands after
the first band has accepted the new temperature. It verifies preservation of
the previous scene and radiance, then retries a valid update and compares it
with a fresh model for both the unreduced and one-layer cases (**44/44**).

### Strict standalone HITRAN records

`Absorption.read_hitran` now builds concrete named columns directly. It rejects
short/non-ASCII records and malformed required fields with filename, line and
field diagnostics. Optional blank statistical weights remain zero by explicit
format rule. The isotopologue field handles digits plus the traditional
`0`/`A`/`B` notation. Existing and new `test_Absorption.jl` test sets all pass.

### Julia interfaces and multiple-dispatch boundaries

- `propertynames(::RTModel)` advertises the same aliases as `getproperty`, so
  `hasproperty` and interactive discovery agree with direct access.
- `isapprox(::GreekCoefs, ::GreekCoefs; kwargs...)` accepts standard tolerances
  and avoids a temporary Boolean array.
- The legacy `noRS` Cabannes helper now fills its vector with an element-type
  appropriate value instead of assigning a scalar to vector storage.
- Generic analytic BRDF Fourier integration and surface-layer construction live
  in `Surfaces/analytic_surface.jl`. Surface-specific files provide reflectance
  physics through dispatch; they no longer host the universal scaffold.
- Standalone profile readers accept `FT=Float32` or `Float64`, with `Float64`
  retained as the default.

Interface regressions pass **6/6**, profile I/O passes **42/42**, and RPV/Ross-Li
smoke checks pass.

### Batched algebra ownership

CPU and CUDA docstrings now define output, scratch, read-only and aliasing roles
for `batch_solve!` and `batch_inv!`, plus their return value and backend-owned
singular handling. The CUDA workspace overload no longer claims to reuse
pivot/info storage that its current CUBLAS call does not consume. A new
mutation-contract set passes **7/7**; the complete batched-kernel file passes
**46/46** on CPU. GPU parity remains covered by the existing hardware suite,
but workspace allocation improvement requires a new measured CUDA change.

### Process-wide numerical setting

`numerics.blas_threads` remains compatible, but its type and public docs now
state the real ownership: a run changes persistent process-wide BLAS state.
Applications should select one value rather than vary it among concurrent
models. A future executor-level API may improve this design without pretending
a per-model field is thread-local.

## Deliberately deferred

- Additional type parameters on `RTModel`, layer carriers and `GreekCoefs` need
  `@code_warntype`/JET, allocation, runtime and compile-latency evidence at one
  hot boundary before expanding specialization.
- Sharing more absorption/Rayleigh preparation between forward, linearized and
  update builders is valuable, but each extraction needs full-build/update and
  finite-difference parity. It should not become a boolean-controlled mega-path.
- Raman's mixed fixed/prepared/scratch state and historical formulas need the
  scientific author's review. Only the demonstrably broken compatibility
  overload and incorrect type labels were changed.
- Per-band surface compatibility, dry-column-conserving legacy VMR reduction,
  end-to-end Metal absorption and CUDA scratch ownership each require a focused
  contract, hardware or scientific acceptance test.
- Broad formatting, comment removal and variable renaming are excluded. The
  current Unicode symbols and explicit operator products make the papers easier
  to map to code; non-obvious stability and mutation comments are retained.

## Validation

Focused Julia 1.12.6 CPU checks use four Julia threads, one BLAS thread and run
from `test/`:

- transactional update regression: **44/44**;
- existing batch-update tests: all test sets pass;
- Julia interfaces: **6/6**;
- standalone absorption/HITRAN: every test set passes;
- batched kernels and ownership contract: **46/46**;
- profile/IO validation: **42/42**;
- complete CPU suite: **12,743 passed, 16 broken/skipped, zero failures/errors**
  in 25m07.6s;
- Aqua ambiguity gate: pass; broad JET remains advisory at the unchanged
  **285 findings**;
- executable documentation contracts: **20/20**; strict
  Documenter/VitePress build: pass;
- A100 CUDA batched ownership/accuracy check: **32/32**;
- clean Julia 1.10.12 resolution/precompile plus interface checks: **6/6**;
- package precompile and `git diff --check`: pass.

The complete CPU rerun covers the combined pressure/source and code-quality
implementation tree. One final non-ASCII parser assertion was then added and
passed in its focused malformed-record set (**13/13**); it changes no source
behavior. The focused A100 result covers the changed batched-algebra contract;
the preceding full **20,301/20,301** GPU runner in `SESSION_HANDOFF.md` remains
the latest whole-GPU-suite result.

The local, untracked manifest was generated by Julia 1.12 and selects
`PrecompileTools` 1.3.4, which itself requires Julia 1.12. Reusing that manifest
directly with Julia 1.10 fails before loading vSmartMOM. A fresh Julia 1.10
environment correctly resolved `PrecompileTools` 1.2.1, precompiled the package
and passed the new interface checks. Version-support CI must resolve each Julia
minor independently rather than treating a 1.12 manifest as a 1.10 lockfile.
