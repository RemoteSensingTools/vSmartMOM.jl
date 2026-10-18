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

## 4. fix(absorption): legacy LUT builder vmr and grids (B5); deprecate standalone line model

- `test_Absorption.jl` + `test_quality.jl` + `test_schema_docs.jl`: 228/228.
  New testset "make_interpolation_model: vmr and spectral grid": 21/21; on
  HEAD's builder 11–12 fail (identical vmr=0/0.5 tables; `MethodError`
  instead of an early `ArgumentError` for wavelength grids).
- Migration recipe run against AtmosphericAbsorption 0.1.2 (CO₂, HITRAN2016
  artifact): agrees with the legacy engine to 2e-6 at vmr = 0; at vmr > 0
  the engines differ through AA's mixed air/self pressure shift (reviewers:
  0.05 % at vmr 0.01 to 2–3.5 % at vmr 0.5, ~11 % at vmr 1).
- The deprecation notice fires once per session from `make_hitran_model`
  only; forward and linearized `model_from_parameters` emit none.
  `make_interpolation_model` stays silent: it remains the builder for YAML
  `LUTfiles` (which read OCO-2 ABSCO HDF, AA `.absco`, and vSmartMOM JLD2).
- GPU check (A100, GPU1): `test_coxmunk_rt_contract.jl` with
  `VSMARTMOM_COXMUNK_ARCH="GPU()"` 10/10 on the B3 commit.
- Reviews: Claude reviewer APPROVE-WITH-FIXES (engine-difference cause, LUT
  route, install hint, vmr docstring) → fixed → APPROVE. Codex
  APPROVE-WITH-FIXES ×3 (grid recipe endpoint, difference wording, LUT
  guidance across README/docs) → fixed.

## 5. docs(citation): HITRAN and related spectroscopic references

- CITATION.bib entries generated from Crossref BibTeX (HITRAN2016/2020/2024,
  TIPS-2021, TIPS-2017, HAPI, HITRAN CIA, MT_CKD 2012 and 2023); DOIs,
  volumes, pages verified against Crossref; `pdflatex`+`bibtex` with
  `\nocite{*}` renders all 15 entries (Claude reviewer).
- Edition guidance checked against `src/Artifacts` (default "artifact" =
  HITRAN2016; downloads = current hitran.org edition; edition string is a
  cache label). Relative documentation links resolve.
- Reviews: included in the reviews of entry 4 (Claude APPROVE, Codex
  citation metadata confirmed).

## 6. merge: molecular Rayleigh API (`release/molecular-rayleigh`)

- Branch commits `c4c14770` (verbatim HWO port), `eaa340d3` (UsingRefInd /
  UsingPolTen), `37f18dd6` (review follow-ups); merged as `a405d87c`.
- `test_molecular_rayleigh.jl` 445/445; with Raman Rayleigh/Cabannes, Aqua and
  IO-export tests 521 assertions (implementer). RI outputs identical to the
  HWO source across all species, ~300 wavelengths, Float64/Float32
  (Claude reviewer); every ExoOptics call form runs, including the A–D RI
  audit loop (36,802 rows).
- UsingPolTen anchored to a closed form from the raw tensor constants
  (rtol 1e-12) and pinned values at 400 nm; RI vs PolTen for N₂/O₂ at
  350–800 nm: σ within −1.5…+0.5 %, ρ within −15…+11 % (anisotropy has no
  dispersion in the Raman-module constants).
- Reviews: Codex APPROVE (both commits, after fixes); Claude reviewer
  APPROVE-WITH-FIXES → fixed (37f18dd6) → APPROVE. Follow-up CHANGELOG
  bullet and docstring note: Codex APPROVE; Claude reviewer
  APPROVE-WITH-FIXES (provenance scope, "two-argument" wording) → fixed.
- Open scientific points for Suniti are in the decision log.

## 7. fix(rt): combined-band solves require one supported common surface (B8)

- `test_lambertian_surfaces.jl` 63/63 (new testset 10/10, message-matched).
  Raman, multisensor, canopy and SS suites (`test_multisensor_heights`,
  `test_ray_cab_raman`, `test_inelastic_vs_smoke`,
  `test_forward_raman_phase1b`, `test_external_solar_raman`, `test_canopy`,
  `test_forward_ss`): 379 pass, 1 pre-existing `@test_broken`.
- Combined elastic solve with a common scalar surface, and with one spectral
  albedo spanning both bands, equals the separate band solves (rtol 1e-12).
- Reviews: Codex REJECT (generic value equality) → rewritten to identity +
  per-type grid checks → APPROVE-WITH-FIXES (scalar-soil canopy) → final
  contract → APPROVE. Claude reviewer APPROVE-WITH-FIXES (spanning spectral
  surfaces, Legendre, early check, message tests) → fixed → APPROVE.
- Downstream: RRS_XCO2, OCORaman and ExoOptics call one band at a time.

## 8. merge: HITRAN cache identity and coverage (`release/b4-hitran-cache`, B4)

- Branch commits `505dffcf`, `20a99b59`, `b6128b2c`, `8e534660`; merged as
  `08347c62`. `test_hitran_cache.jl` 202/202 on the merged tree; Aqua and IO
  exports 38/38.
- Live checks against hitran.org (throwaway depots, 2 + 3 requests): the first
  exposed a blocker (Content-Length counts one extra byte per record, so every
  download raised a transport error); after the fix, O2 13000–13010 (23
  records) and the full O2 set (30,762 records) download, validate, verify by
  SHA-256 and are reused without re-download.
- Reviews: Codex REJECT ×2 → APPROVE (round 5); Claude reviewer REJECT (live
  blocker) → APPROVE.

## 9. cleanup: Raman type docstrings, Fourier-bound docstrings, linearized band selection

- Docstrings of `RRS`, `VS_0to1`, `VS_1to0` were rotated (B6), and those of
  `sol_VS_0to1`, `sol_VS_1to0`, `VS_0to1_plus`, `VS_1to0_plus`,
  `sol_VS_0to1_plus` and `sol_VS_1to0_plus` described `RRS`; `Base.Docs.doc`
  now returns the right text for all nine. `_derive_m_max_bands_via_traits` and the main `rt_run`
  docstring now describe the derived per-band Fourier bound (B6).
- `test_lin_band_selection.jl` 3/3: `[1]` equals band 1 bit for bit; `[1, 1]`
  is rejected with a clear message (previously an unclear downstream error).
