# AGENTS.md — Quick onboarding for AI coding agents

This file is the **first read** for any agent (Claude Code, Codex, etc.) joining
work on vSmartMOM.jl. The goal is a 3-minute orientation: what the package is,
the narrative thread it's organized around, where to find each concept, and the
conventions that keep the codebase coherent.

**Sister files:**
- [CLAUDE.md](CLAUDE.md) — project conventions, build/test commands, file structure.
- [docs/dev_notes/theory_references.md](docs/dev_notes/theory_references.md) — verified equation↔code crib sheet (the source of truth for paper citations and `file.jl:LINE` anchors).
- [docs/src/pages/concepts/](docs/src/pages/concepts/) — the user-facing narrative arc; same thread, polished prose.

## What vSmartMOM.jl is

A Julia package for **vectorized polarized atmospheric radiative transfer** using
the **Matrix Operator Method (MOM)**. Computes Stokes-vector radiances (I, Q, U,
V) leaving an inhomogeneous, scattering, absorbing, polarized atmosphere
illuminated by the sun, plus **analytic Jacobians**, on **CPU / CUDA / Metal**
via `KernelAbstractions`. Used for retrieval of trace gases, aerosols, surface
properties from hyperspectral satellite instruments (OCO-2/3, GOSAT, EMIT, …).

## The narrative thread (memorize this — everything hangs off it)

```
PROBLEM      Polarized radiance from a layered scattering+absorbing atmosphere.
                │
                ▼
WHY MOM?     Per-layer matrix operator → linear-algebra composition →
             full multiple scattering free, large τ cheap (logarithmic via
             doubling), exact analytic Jacobians.
                │
                ▼
DISCRETIZE   Fourier in azimuth (m=0..max_m), Gauss/Radau quadrature in μ.
             Per Fourier moment m, each layer reduces to FOUR ARRAYS of
             shape (NquadN, NquadN, nSpec):
               τ      — optical depth          ϖ      — single-scatter albedo
               Z⁺⁺   — forward phase matrix    Z⁻⁺   — backscatter phase matrix
                │
                ▼
BUILD LAYER  ┌──────────────┐    ┌──────────────────────────────┐
OPTICS       │ Gas opacity  │    │ Scattering (Rayleigh+aerosol)│
             │  τ_abs(ν,z)  │    │  Mie/NAI2/PCW → GreekCoefs   │
             │  HITRAN LBL  │    │  compute_Z_moments → Z⁺⁺,Z⁻⁺ │
             └──────┬───────┘    └──────┬───────────────────────┘
                    └──────┬────────────┘
                           ▼
                constructCoreOpticalProperties:
                • τϖ-weighted mixing of scatterers   (operator + on
                  CoreScatteringOpticalProperties)
                • δ-M truncation (forward-peak removal)
                • Add gas absorption to total τ
                • Vertical concatenation of layers   (operator *)
                → Vector{CoreScatteringOpticalProperties}, one per layer
                │
                ▼
SOLVE (MOM)  For each Fourier moment m, for each layer iz from TOA to BOA:
             1. Elemental:    thin-layer (r, t, j) — EXACT finite-δ
                              single-scatter formulas (Fell 1997 /
                              SF2023-II Eqs 10-11), NOT the linear
                              S2014 Eqs 19-20 limit.
             2. Doubling:     thicken via geometric series (I − R·R)⁻¹
                              + D-matrix symmetry → halves doubling cost
             3. Adding:       stack layer onto composite via
                              ScatteringInterface_{00,01,10,11} dispatch
                │
                ▼
SURFACE      BRDF as the bottom-most AddedLayer. Cox-Munk (polarized),
             RPV, Ross-Li, Lambertian, canopy. Same matrix-operator language.
                │
                ▼
LINEARIZE    Each operator (elemental/doubling/interaction) has an EXACT
             tangent-linear partner. ForwardDiff is used UPSTREAM at the
             optical-property boundary; the RT kernel is hand-differentiated.
             ParameterLayout names the Jacobian columns.
                │
                ▼
EVERYWHERE   ONE source tree, three GPU backends. Kernels written once with
ON ANY       KernelAbstractions @kernel. CUDA / Metal injected via Julia 1.9+
HARDWARE     package extensions. ALL WAVELENGTHS PROCESSED IN PARALLEL: the
             third array dim is spectral; batched_mul broadcasts over it.
             LBL on GPU is the design point. ForwardDiff Duals flow through
             GPU batched paths. No Fortran. No two-language problem.
```

That's the spine. Every page in `docs/src/pages/concepts/` is one segment of
that thread.

## Two non-obvious tricks that compound (read this twice)

These two design choices make line-by-line RT on the GPU practical. They're
load-bearing in the prose; if you're explaining the package to anyone, lead
with these:

1. **Elemental kernel uses the EXACT finite-δ single-scatter formulas, not
   the linear approximation.** S2014 Eqs (19)–(20) are written in the
   `δ → 0` limit (linear in δ; equivalent to `1−exp(−x) ≈ x`). Many MOM codes
   stop there and need very thin elemental layers (large `N_doubl`) for that
   approximation to hold. vSmartMOM's `elemental.jl:207–252` instead uses
   `1 − exp(−δτ(1/μᵢ+1/μⱼ))` and `exp(−δτ/μᵢ) − exp(−δτ/μⱼ)` written via
   `-expm1(-x)` and `expdiff_neg(a,b)` for numerical stability. Result:
   thicker elemental layer at the same single-scatter accuracy → smaller
   `N_doubl` → less doubling round-off, stable in `Float32`.

2. **The elemental layer is sized by SCATTERING optical depth, not total
   optical depth.** `get_dtau_ndoubl` (`rt_kernel.jl:245–253`) picks
   `N_doubl` from `τ_scat·ϖ`; absorption is layered in via
   `τ_λ = τ_abs + τ_scat`, `ϖ_λ = τ_scat·ϖ / τ_λ`. The transmission factors
   inside the elemental kernel use the per-wavelength `dτ_λ`. Consequence:
   a layer with `τ_abs ≈ 50` (deep gas absorption) and `τ_scat ≈ 0.05`
   doesn't blow up `N_doubl`. SF2023-II Eqs (8)–(9) call this the
   "constant-`N_doubl` trick"; it's why a hyperspectral grid runs as one
   batched call per layer per Fourier moment.

These two compose: small `N_doubl` (from #1) *and* `N_doubl` constant across
the spectral axis (from #2) → all wavelengths run together on GPU at full
bandwidth.

## What sets vSmartMOM apart (anchored)

Not marketing — each claim has a `file:line`. Mirror this list in user-facing
material when explaining the package:

1. **Operator-level analytic linearization.** RT kernel itself is hand-differentiated; AD is upstream-only. **Combined forward + linearized run costs less than 2× a forward-only run** — the closed-form chain rule on adding-doubling reuses the already-computed batched matrix inverses, so the dominant LU work is paid once. ForwardDiff would pay `(1+N_params)×`; finite differences `(1+N_state)×`. — `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl`, `{elemental,doubling,interaction}_lin.jl`. See `docs/src/pages/concepts/06_linearization.md` § "Why this is fast". *Status: production codebase implements the kernel-level fast path; a cleaner AD-upstream → analytic-downstream boundary (more idiomatic handoff struct, less manual chain-rule expansion in upstream code) is in active development.*
2. **One `@kernel` source compiles for CPU, CUDA, and Metal.** — `src/Architectures.jl:33–96`, `ext/vSmartMOMCUDAExt.jl:21–27`, `ext/vSmartMOMMetalExt.jl:19–22`.
3. **Hybrid AD across the GPU boundary.** `ForwardDiff.Dual` flows through `NNlib.batched_mul` on `CuArray`. — `ext/gpu_batched_cuda.jl:141–177`.
4. **Polarization is a type, not a runtime branch.** `Stokes_I/IQ/IQU/IQUV` specialize the kernels at compile time. — `src/Scattering/types.jl:92–143`.
5. **Optical properties as algebra.** `+` mixes scatterers; `*` stacks layers. — `src/CoreRT/types.jl:1063–1101`.
6. **Weak GPU dependency** (Julia 1.9+ package extensions). CPU-only installs first-class. — `ext/vSmartMOMCUDAExt.jl`.
7. **Three `(τ, ϖ, Z)` core variables per layer.** RT kernel differentiates against these directly. — `src/CoreRT/types_lin.jl:119–149`.
8. **Exact finite-δ elemental** (point #1 above). — `src/CoreRT/CoreKernel/elemental.jl:207–252`.

## Where everything is

### The Concepts arc (the docs version of the spine above)

| # | Page | What it covers |
|---|---|---|
| 01 | [`docs/src/pages/concepts/01_overview.md`](docs/src/pages/concepts/01_overview.md) | Problem statement; MOM thesis; vSmartMOM differentiators |
| 02 | [`docs/src/pages/concepts/02_rt_theory.md`](docs/src/pages/concepts/02_rt_theory.md) | Vector RTE; Fourier+quadrature discretization; polarization-as-type |
| 03 | [`docs/src/pages/concepts/03_layer_optics.md`](docs/src/pages/concepts/03_layer_optics.md) | The `(τ, ϖ, Z⁺⁺, Z⁻⁺)` per-layer abstraction (the bridge) |
| 03a | [`docs/src/pages/concepts/03a_absorption.md`](docs/src/pages/concepts/03a_absorption.md) | How `τ_abs` is computed (HITRAN LBL, line shapes) |
| 03b | [`docs/src/pages/concepts/03b_scattering.md`](docs/src/pages/concepts/03b_scattering.md) | How `Z⁺⁺/Z⁻⁺` and `ϖ`, `k` are computed (Mie, Rayleigh, GreekCoefs, NAI-2 vs PCW) |
| 03c | [`docs/src/pages/concepts/03c_mixing.md`](docs/src/pages/concepts/03c_mixing.md) | Mixing scatterers, scattering-vs-total τ split, δ-M truncation |
| 04 | [`docs/src/pages/concepts/04_mom_solver.md`](docs/src/pages/concepts/04_mom_solver.md) | Elemental → Doubling → Adding (the kernel sequence) |
| 05 | [`docs/src/pages/concepts/05_surfaces.md`](docs/src/pages/concepts/05_surfaces.md) | Surface BRDFs as the bottom AddedLayer |
| 06 | [`docs/src/pages/concepts/06_linearization.md`](docs/src/pages/concepts/06_linearization.md) | Operator-level chain rule, ParameterLayout |
| 07 | [`docs/src/pages/concepts/07_architecture.md`](docs/src/pages/concepts/07_architecture.md) | Architecture-agnostic kernels, GPU, Julia advantages |
| 08 | [`docs/src/pages/concepts/08_inelastic.md`](docs/src/pages/concepts/08_inelastic.md) | Raman / Cabannes (brief) |

### Verified equation ↔ code map

[`docs/dev_notes/theory_references.md`](docs/dev_notes/theory_references.md) is the source of truth: every paper equation that the docs quote has a `file.jl:LINE` reference. When in doubt about a citation or symbol, check there first.

### Code map (key anchors)

| Concept | Source |
|---|---|
| Top-level RT loop (bands × m × layers) | `src/CoreRT/rt_run.jl:53–329` |
| Per-layer kernel dispatch | `src/CoreRT/CoreKernel/rt_kernel.jl:48–229` |
| Build per-layer (τ, ϖ, Z) | `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties.jl:11–65` |
| τϖ-weighted scatterer mixing (`+`) | `src/CoreRT/types.jl:1063–1093` |
| Vertical layer stacking (`*`) | `src/CoreRT/types.jl:1096+` |
| δ-M truncation | `src/CoreRT/LayerOpticalProperties/delta_m_truncation.jl:44–48` |
| Mie / Greek (NAI-2 vs PCW) | `src/Scattering/compute_NAI2.jl:44`, `compute_PCW.jl:28` |
| Phase matrix Z from Greek + μ | `src/Scattering/compute_Z_matrices.jl::compute_Z_moments` |
| HITRAN absorption | `src/Absorption/compute_absorption_cross_section.jl:32–280` |
| Elemental (r, t, j) | `src/CoreRT/CoreKernel/elemental.jl:207–252` |
| Doubling (geom series + D-symmetry) | `src/CoreRT/CoreKernel/doubling.jl + rt_helpers.jl:88–122` |
| Interaction (4 cases) | `src/CoreRT/CoreKernel/interaction.jl:14–136` |
| Surface BRDFs | `src/CoreRT/Surfaces/` (10 files) |
| Surface coupling | `src/CoreRT/CoreKernel/interaction_hdrf.jl:1–42` |
| Postprocessing (azimuth, VZA) | `src/CoreRT/tools/postprocessing_vza*.jl:23–99` |
| Quadrature (Gauss / Radau) | `src/CoreRT/tools/rt_set_streams.jl:24–110` |
| Linearization kernels | `src/CoreRT/CoreKernel/{elemental,doubling,interaction}_lin.jl` |
| Chain-rule expansion | `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl:1–100` |
| Jacobian column layout | `src/CoreRT/parameter_layout.jl:1–67` |
| Architecture types | `src/Architectures.jl:1–98` |
| `@kernel` apply_D! (D-matrix) | `src/CoreRT/CoreKernel/doubling.jl:85–110` |
| `@kernel` line shape | `src/Absorption/compute_absorption_cross_section.jl:229–280` |
| `@kernel` Mie coefficients | `src/Scattering/gpu_mie_kernels.jl:30–100` |
| `@kernel` portable LU (Metal) | `src/CoreRT/tools/ka_batched_kernels.jl:100–178` |
| Batched matmul CPU | `src/CoreRT/tools/cpu_batched.jl:24–73` |
| Batched matmul CUDA | `ext/gpu_batched_cuda.jl:122–139` |
| Batched matmul Metal | `ext/vSmartMOMMetalExt.jl:29–54` |
| ForwardDiff.Dual through GPU | `ext/gpu_batched_cuda.jl:141–177` |
| CUDA extension init | `ext/vSmartMOMCUDAExt.jl:21–66` |
| Metal extension init | `ext/vSmartMOMMetalExt.jl:19–101` |

## Conventions you must follow

(Detailed in CLAUDE.md; condensed here.)

- **Float-type generic.** Most types/functions are parameterized by `FT`. Don't hard-code `Float64`.
- **`_lin` suffix** = linearized (Jacobian) variant of a function or type.
- **Unicode in identifiers** is used directly: `τ`, `ϖ`, `μ`, `μ₀`, `Z⁺⁺`, `Z⁻⁺`, `R⁻⁺`, `T⁺⁺`, `r⁻⁺`, `t⁺⁺`. Don't romanize.
- **Sign convention:** `+` = incoming/downward, `−` = outgoing/upward.
- **Layer naming:** `CompositeLayer` fields are uppercase (`R, T, J`); `AddedLayer` fields are lowercase (`r, t, j`).
- **3D RT array layout:** `(NquadN, NquadN, nSpec)` where `NquadN = Nquad * n_stokes`. The third dim is **spectral** — that's what enables batched-matmul over wavelengths.
- **Spectral units:** wavenumber (cm⁻¹) internally; wavelength (μm) for Mie.
- **Profile direction:** TOA-to-BOA internally. GEOSChem data is BOA-to-TOA and gets flipped on read.
- **No co-author attribution in commits.**
- **Test discipline:** must `cd test/` before running tests — test files use relative paths (`test_profiles/`) that only resolve from `test/`.

## How to run

```bash
# Tests
cd test && julia --project=. runtests.jl

# Forward RT
julia --project=. -e '
  using vSmartMOM
  params = parameters_from_yaml("config/quickstart.yaml")
  model  = model_from_parameters(params)
  R, T   = rt_run(model)
'

# Linearized RT (Jacobians)
julia --project=. -e '
  using vSmartMOM
  params = parameters_from_yaml("config/ocean_coxmunk.yaml")
  model, lin_model = model_from_parameters(LinMode(), params)
  NAer = length(params.scattering_params.rt_aerosols)
  NGas = size(lin_model.tau_dot_abs[1], 1)
  R, T, dR, dT = rt_run(model, lin_model, NAer, NGas, 1)
'

# Build docs
cd docs && julia --project=. make.jl
```

## When in doubt

1. **What is the package about?** → Read the spine above (top of this file).
2. **What does this code do?** → Find it in the code-anchor table above; read the relevant Concepts page; cross-reference `theory_references.md` for the paper equation.
3. **Why is the code shaped this way?** → Check the **design-choice passages** at the bottom of `theory_references.md`. There are seven; they explain non-obvious trade-offs (no solar SFI in J, constant-`N_doubl`, Cabannes vs Rayleigh greek, linear-in-inelastic, `rt_run_ss` raison d'être, δBGE-fit vs δ-m, exact finite-δ elemental).
4. **Where does this paper equation map in the code?** → `theory_references.md` Sections A–J.
5. **How should I write a new feature?** → Check `CLAUDE.md` for the "Common Workflows" patterns (adding a surface BRDF, adding a test, modifying YAML parsing).

## Maintenance contract for this file

When the docs evolve, this file evolves with them. Specifically:

- If a new Concepts page is added or the arc reorders, update the Concepts arc table.
- If a kernel's source location moves, update the code-anchor table.
- If a new differentiator is added or a current claim no longer holds (e.g., Mie becomes default-GPU), update the differentiator list.
- If a new design-choice passage is added to `theory_references.md`, add it to "When in doubt" §3.

This file is for agents who don't have time to read the whole codebase. If it's stale, it's actively misleading. Keep it honest.
