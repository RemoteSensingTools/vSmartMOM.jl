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

## 3. fix(coxmunk): polarized geometry, argument order, Fourier storage sign and V terms (B3)

- `test_coxmunk.jl` + `test_coxmunk_polarization.jl` + `test_coxmunk_rt_contract.jl`:
  992/992 pass (Claude reviewer run); polarization file alone 726/726.
  `test_surface_split.jl`, `test_component_m_max.jl` (with the Cox–Munk files):
  357/357. `StandaloneSS/runtests.jl`: 242 pass, 1 pre-existing `@test_broken`.
- VLIDORT baseline Case D (new): 48/48 — Rayleigh layer (τ = 0.1, 0.3,
  ρ = 0.03) over GISS Cox–Munk (U = 2, 7 m/s, shadowing on/off), Lambertian
  and Ross–Li (80° cap), NSTREAMS 16, TOA and BOA (BOA excludes the μ = μ₀
  view, which carries the direct beam in vSmartMOM): I within 7e-5, Q and U
  within 1e-5 of max I. References from an external driver calling
  unmodified VLIDORT 2.8.3; driver, inputs and full outputs archived in
  `/home/sanghavi/vsmartmom-release-evidence/2026-10-17/` with SHA256SUMS.
- Defect detection: HEAD fails 32/40 stored VLIDORT kernel matrices (Codex)
  and 274/521 of the original polarization tests in all six testsets (Claude
  reviewer); HEAD violates full-matrix reciprocity at order 1. With only the
  old Fourier sign restored, Case D fails 22/48 (I off by up to 3.2e-3,
  Q/U by 1.2e-3).
- Independent references: Jones-vector field calculation whose Stokes basis
  is anchored to the code's own Rayleigh phase matrix (with and without the
  Fourier kernel); stored VLIDORT GISS and GCMCRI (complex index) kernels;
  Hovenier reciprocity; the V sign was also cross-checked by the reviewer
  against vSmartMOM's Mie code for a large absorbing sphere.
- Reviews: Claude reviewer APPROVE-WITH-FIXES (docs, latent StandaloneSS
  kernel, extra tests) → fixed → APPROVE. Codex APPROVE-WITH-FIXES (Float32
  cancellation in cos θ_local, conventions page, comment) → fixed → APPROVE.
- Not run: CUDA hardware, Documenter build, full suite.
