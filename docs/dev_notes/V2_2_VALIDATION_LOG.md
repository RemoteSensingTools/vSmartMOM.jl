# vSmartMOM v2.2 consolidation — validation log

One entry per commit on `release/v2.2-consolidation`. A commit cannot carry
its own SHA, so each entry names the change; the SHA appears in `git log`.
Environment unless stated: host curry.gps.caltech.edu, Julia 1.12.5,
AtmosphericAbsorption 0.1.2, CPU (`CUDA_VISIBLE_DEVICES=''`),
`JULIA_NUM_THREADS=4`, `OPENBLAS_NUM_THREADS=1`, tests run from `test/`.

Review gate for every commit: focused tests → independent Claude reviewer →
Codex (`gpt-6.1-sol`, reasoning effort xhigh, read-only) → fixes → re-review
by both.

## 1. docs: decision log, ExoOptics audit, handoff copy

Docs only (`docs/dev_notes/`). Handoff copied verbatim from Christian
Frankenberg's checkout (SHA-256 identical to the source file).
- Claude reviewer: factual errors in §1.1, §1.3, §1.5, §2.2 → fixed; re-review
  confirmed, final wording nits applied as suggested.
- Codex: APPROVE-WITH-FIXES (four documentation findings) → fixed; re-review
  APPROVE-WITH-FIXES (two wording nits) → applied as suggested.

## 2. fix(coxmunk): per-axis slope variance (B1) and Smith shadowing normalization (B2)

- `test_coxmunk.jl`: 256/256 pass; `test_coxmunk_rt_contract.jl`: 10/10;
  `StandaloneSS/runtests.jl`: 242 pass, 1 pre-existing `@test_broken`.
  Reviewer runs also: `test_surface_split.jl` 64/64,
  `test_component_m_max.jl` 27/27, `test_schema_docs.jl` 67/67,
  external-solar parity 73/73.
- Defect detection: with HEAD's helpers restored, the new testsets fail
  (Claude reviewer: 55 failures + 4 errors; Codex: 12 failures from the old
  variance, 45 from the old Smith functions; specular closed form 0.5× in all
  12 cases).
- Independent references: PDF normalization and second moment by direct
  quadrature (Float64, Float32); Smith Λ by quadrature of the projected-slope
  integral and the handoff's values; complete specular BRDF in closed form
  with Fresnel reflectance from Snell's law.
- Claude reviewer and Codex: no code defects; documentation, release-note
  and test additions requested and made; both re-reviews passed.
- Not run in this commit: CUDA, the full suite, documentation build.
