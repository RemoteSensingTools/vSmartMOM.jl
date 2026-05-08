# Refactoring roadmap (post-v2.1.0)

These items were surfaced during the v2.1.0 release-readiness audit but
deferred — they're real opportunities, but each is a project, not a quick
win, and bundling them with the source-term + Fourier-stream-resolution
work would have made the release brittle. Each entry names the *what* and
the *why-deferred* so a future PR can pick it up cleanly.

The audit lens was: code cleanliness, conciseness, code duplication,
Julian style, and ease of understanding for physical scientists.

---

## 1. Forward / linearized kernel duplication (`*_lin` variants)

**Files**

| forward                            | linearized                            | overlap |
|------------------------------------|---------------------------------------|---------|
| `src/CoreRT/CoreKernel/elemental.jl`   | `…/elemental_lin.jl`     | ~70% structural |
| `src/CoreRT/CoreKernel/doubling.jl`    | `…/doubling_lin.jl`      | ~70% structural |
| `src/CoreRT/CoreKernel/interaction.jl` | `…/interaction_lin.jl`   | ~70% structural |
| every `Surfaces/*_surface.jl`          | `…/*_surface_lin.jl`     | ~80% structural |

The `_lin` variants re-implement the same math and add product-rule
derivative propagation alongside the forward computation.  Aggregate
duplication is ~1500-2000 lines.

**What a future PR would do.**  Introduce a shared `*_core!(::Val{:fwd}
or ::Val{:lin}, …)` kernel.  The `:lin` path branches into derivative
assembly only at the few product-rule sites; everything else lives once.
Surfaces would benefit from an `@autolinearize` macro that wraps the
forward routine and derives `_lin` mechanically (the analytic surface
BRDFs all follow the same chain-rule pattern).

**Why deferred.**  High blast radius (every Jacobian regression test must
re-validate to FD tolerance), and the structure has to remain bit-equal
across forward and lin so users in the field don't see numerical drift.
Best done in its own PR with an `--exact` Jacobian regression gate.

---

## 2. NAI2 / NAI2_lin / NAI2_gpu trio

**Files**: `src/Scattering/compute_NAI2.jl` (~580 lines),
`compute_NAI2_lin.jl` (~645 lines), `compute_NAI2_gpu.jl` (~225 lines).
Combined ~1450 lines of largely the same Mie + Greek-coefficient logic
implemented three times.

**What a future PR would do.**  Extract the core
`compute_mie_and_greeks(::MieModel, ::Val{FT})` kernel that is generic in
the float type — ForwardDiff `Dual` makes the linearized variant
disappear automatically; KernelAbstractions makes the GPU variant
disappear automatically.  Should leave one base implementation plus thin
dispatchers.

**Why deferred.**  Larger and more invasive than item 1, and the GPU
variant has hand-tuned launch shapes that don't map 1:1 onto the CPU
loop structure — needs careful benchmarking before merging.

---

## 3. `compute_Z_moments` overload pair

**Location**: `src/Scattering/compute_Z_matrices.jl`, two near-duplicate
overloads (lines 1–89 and 95–218).  The second adds μ₀ direct-beam
handling on top of the first; ~80% of the A-matrix assembly loop is
copied verbatim.

**What a future PR would do.**  Factor out a private
`assemble_A_matrices!(A⁺⁺, A⁻⁺, μ_set, 𝐁_all, l_range)` helper; call it
once with the quadrature `μ` set, optionally a second time with `[μ₀]`
gated by an `include_direct=false` keyword.

**Why deferred.**  Small, low-risk pilot — would make a good
"refactor practice" PR, but didn't fit the v2.1 scope freeze.  Pick this
one first when starting on item 1.

---

## 4. Mutable parameter / config structs

**Locations**: `src/CoreRT/types.jl` lines 502 (`CanopySurface`), 601
(`AbsorptionParameters`), 629 (`ScatteringParameters`), 702
(`vSmartMOM_Parameters`); `src/Scattering/types.jl:79` (`Aerosol`);
`src/CoreRT/types.jl:83` (`RT_Aerosol`).

**What a future PR would do.**  Audit each: if the only mutation point
is during construction, change `Base.@kwdef mutable struct` →
`Base.@kwdef struct`.  `RTWorkspace` (line 1278) is genuinely mutable
(pre-allocated buffers) — leave it alone.

**Why deferred.**  Needs a careful read of every code path that touches
each struct to confirm no in-place updates are made post-construction;
some downstream tools may rely on `setfield!` access for config
overrides.  Worth doing in a focused PR with a CI gate that exercises
`Pkg.test` after each individual struct conversion.

---

## 5. `Aerosol.nᵣ`, `Aerosol.nᵢ` untyped fields

**Location**: `src/Scattering/types.jl:79–86`, `mutable struct Aerosol{}`.

Fields hold real and imaginary refractive indices but lack type
constraints, defeating specialization.

**What a future PR would do.**  Parameterize: `struct
Aerosol{FT<:Real}` with `nᵣ::FT`, `nᵢ::FT`.  Combines naturally with
item 4.

**Why deferred.**  Tiny diff but leaks into every Mie call site that
constructs an `Aerosol`; do as part of item 4 to keep ripple-fixes in
one PR.

---

## 6. Pre-existing diagnostics in CoreKernel

The Julia language-server flags ~50 `UnusedBinding` /
`UnusedFunctionArgument` / `IncorrectCallArgs` lints across
`elemental.jl`, `elemental_lin.jl`, `doubling.jl`, `doubling_lin.jl`,
`interaction.jl`, and the surface BRDFs.  Most are stub arguments
preserved across the forward/lin signature contract.

**What a future PR would do.**  Sweep with `_` prefixes
(`_unused_arg`) for arguments that are intentionally part of a method
contract but not used by this dispatch, and remove genuinely dead
local bindings.  Aqua doesn't catch these — the IDE does.

**Why deferred.**  Touches every dispatch in CoreKernel; better to fold
it into item 1 (the forward/lin unification) so we don't have to do the
churn twice.

---

## 7. `sandbox/` 79 MB of tracked content

**Files**: `sandbox/rami/RAMI4ATM_experiments_v1.0.json` (46 MB),
`sandbox/rami/ref_kurucz_2006.nc` (41 MB), plus `*.f90`, `.png`s, and
prototype `.jl` scripts (≈ 35 MB total).  Every `git clone` pulls all of
this.

**What a future PR would do.**  Move the two large data files to the
same artifact / scratch-cache pattern used for HITRAN
(`src/Artifacts/artifact_helper.jl`), or to git-lfs.  Keep the
prototype `.jl` and `.f90` files in tree (they're small and useful for
reference).

**Why deferred.**  User explicitly chose to leave `sandbox/` alone for
v2.1.  Likely fits a v2.2 hygiene PR.

---

## 8. Physicist-readability sweep beyond CoreKernel + 2 BRDFs

In v2.1 we added module-header docstring blocks to `elemental.jl`,
`doubling.jl`, `interaction.jl`, and physics-named-step blocks to
`lambertian_surface.jl` and `rpv_surface.jl`.  The same treatment
would benefit:

- `src/CoreRT/Surfaces/coxmunk_surface.jl` (Cox–Munk wave-slope ocean BRDF)
- `src/CoreRT/Surfaces/rossli_surface.jl` (kernel-driven sparse vegetation)
- `src/CoreRT/Surfaces/canopy_surface.jl` (multi-layer canopy)
- `src/Scattering/compute_Z_matrices.jl` (phase-matrix Fourier moments)
- `src/Scattering/compute_NAI2.jl` (Greek-coefficient assembly)

**Why deferred.**  Documentation-only, no urgency, low risk.  Do
opportunistically (a good "while I'm in here" task for any future PR
that touches these files).

---

## 9. Linearized Mie Jacobian column-count bug

**Location**: `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties_lin.jl`
around lines 322-330.

```julia
ω̃̇ = arr_type(zeros(n,7))
ω̃̇[:,2:5] .= arr_type(collect(tmpω̃̇'))   # hard-coded col 2:5 (4 wide)

ḟᵗ_vec = collect(ḟᵗ)  # comment says "keep as 4-element vector from Mie"
ḟᵗ = arr_type(zeros(n,7))
for k in eachindex(ḟᵗ_vec)
    ḟᵗ[:,1+k] .= ḟᵗ_vec[k]              # writes col 1+k; out-of-bounds when len(ḟᵗ_vec) > 6
end
```

The 7-column matrix is sized for the 7-parameter aerosol layout
(`[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]`) but the loop assumes `ḟᵗ_vec`
length 4 (Mie-only).  Under the v0.7 migration of `JacobianTestFast.yaml`
(legacy `max_m=3, l_trunc=5, Δ_angle=2.0` → `nstreams=3, truncation: auto`)
the docs build's `_jacobian_aod_plot` triggers a `BoundsError` here.

**Tactical workaround in v2.1**: `test/test_parameters/JacobianTestFast.yaml`
was held back on the legacy schema — the only one of 48 fixtures not
migrated.  All other Jacobian regression tests pass on the migrated
fixtures.

**What a future PR would do.**  Replace the hard-coded `2:5` and `1+k`
with the actual Mie-derivative column indices (probably `2:size(ḟᵗ_vec,1)+1`
clamped to the matrix width, plus an `@assert` that the input length
matches `ParameterLayout`'s expectation).  Then re-migrate
`JacobianTestFast.yaml` to v0.7.

**Why deferred.**  The fix touches the analytic-Jacobian construction
path; needs Jacobian regression validation across every YAML in
`test/test_parameters/` to confirm no numerical drift.  Out of scope
for v2.1's "quick wins" release.

---

## Quick wins delivered in v2.1.0 *(for reference)*

So the next person doesn't redo these:

- CHANGELOG header rename, `Project.toml` 2.0.0 → 2.1.0
- CI matrix Julia 1.10/1.11/1.12 + `actions/checkout@v4` + cache step
- Dropped dead `using UnicodePlots` from `src/vSmartMOM.jl`,
  `src/CoreRT/CoreRT.jl`; removed `UnicodePlots`, `OrderedCollections`
  from `Project.toml` (Aqua was about to flag them)
- Migrated `config/quickstart.yaml` to v0.7 schema (`nstreams: 3`,
  `truncation: NoTruncation()`)
- README: added forward + Jacobian copy-paste snippets
- CoreKernel + Lambertian/RPV physics-tour comment headers
