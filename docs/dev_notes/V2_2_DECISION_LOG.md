# vSmartMOM v2.2 consolidation — decision log and open items

Maintained during the v2.2 consolidation (branch `release/v2.2-consolidation`,
started 2026-10-17 from `7acab570`, the Round 5 integration). Each entry says
*what* was decided, *why*, and *what to look at* if the decision needs to be
revisited. Decisions are Suniti Sanghavi's unless stated otherwise.

The full release plan is Christian Frankenberg's handoff
`RELEASE_READINESS_AND_CLAUDE_HANDOFF_2026-10-16.md`; this file records only
decisions taken on top of it.

## 1. Decisions

### 1.1 Molecular Rayleigh: two methods, selected by dispatch

**Planned (not yet implemented in this branch).** The HWO refractive-index
Rayleigh code (`src/Scattering/molecular_rayleigh.jl` in the HWO worktree)
will be ported with two method types:

| Type | Physics | Species available now |
|---|---|---|
| `UsingRefInd()` | refractive index + (effective) King factor | He, Ar, N₂, O₂, CO₂ |
| `UsingPolTen()` | polarizability tensor (mean polarizability ᾱ, anisotropy γ) | N₂, O₂, H₂ (constants from the Raman module) |

`RefractiveIndexRayleigh` will stay as an alias of `UsingRefInd` so existing
ExoOptics code keeps working. The two-argument call
`molecular_rayleigh_properties(gas, λ)` will remain the refractive-index method
(ExoOptics calls it for He, which has no polarizability tensor yet).

Note on what the RT model uses today: the default *elastic* Rayleigh optical
depth is the Bodhaine et al. (1999) Eq. 30 air formula scaled by
p_surf/1013.25 and rescaled by the ratio of King factors that takes
Bodhaine's built-in depolarization (0.0279) to the configured one
(`getRayleighLayerOptProp`).
With the default `depol: -1`, that depolarization — which also sets the
Rayleigh Greek coefficients — comes from the N₂/O₂ polarizability constants
of the Raman module (`model_from_parameters.jl`, `getRamanAtmoConstants`).
The Cabannes/RRS/VRS paths use the polarizability tensors directly.

### 1.2 Common names between vSmartMOM and ExoOptics

Methods live in their most general form in vSmartMOM; ExoOptics derives from
them. Where ExoOptics duplicates a vSmartMOM method under another name, the
ExoOptics name changes to the vSmartMOM one (or delegates to it). Where no
exact parallel exists, a vSmartMOM method is added via multiple dispatch.
An audit of ExoOptics is in progress (see §3).

### 1.3 Solar direction: embedded μ₀ is the default everywhere

`external_solar` controls the *direction* of the solar beam in the
quadrature, not the solar spectrum:

- `external_solar=false` (**default**, embedded): μ₀ is an extra quadrature
  node — zero-weight with Gauss quadrature, a weighted node with Radau.
  Supports every output (TOA, BOA, interior sensors, full R/T, Radau,
  atmosphere/surface split).
- `external_solar=true` (external): μ₀ stays outside the stream set and the
  exact direct-beam source is evaluated at the requested μ₀ through
  rectangular source columns; TOA-only, Gauss-only, through `rt_run_toa`.
  It saves one operator node per solve (unless μ₀ already coincides with a
  node). Supported only for forward `noRS` and `RRS` through `rt_run_toa`,
  linearized `noRS` through `rt_run(model, lin_model, …)`, Lambertian
  surfaces, no `rt_run_ss`, no interior sensors. Retained for quick-and-dirty use only.

**Not to be confused with stream interpolation.** Reusing one solve for many
solar geometries by interpolating between quadrature nodes is a *separate*
algorithm, `rt_run_streams` (per-Fourier-moment results at all streams). It
runs with the embedded representation (`external_solar=false`). ExoOptics'
production stellar path (`PerPixelRT._per_pixel_from_streams`) uses it; its
documented interpolation error (≈9 % DOLP at α = 90° with 5 streams) belongs
to that path, not to `external_solar`.

The solar *spectrum* is a separate setting: `SolarBeam(F₀=...)` (default
F₀ = 1) or the Raman `RS_type.F₀`. Raman runs need the real solar spectrum;
ExoOptics uses F₀ = 1 so that incident starlight stays generic. Both rules
are independent of `external_solar`.

Embedded is the default for the release *and* for ExoOptics. The HWO work
should also have used embedded; the ExoOptics production paths that force
the external mode are being inventoried (§3) before switching.
The Horak pure-gas diagnostic keeps `external_solar=true` explicitly (one
keyword added to its `model_from_parameters` call).

### 1.4 Workflow layout

All workflow source lives in `sandbox/workflows/<NAME>/` with its own
documentation. Generated data stays outside git and is located through
environment variables (`RRS_XCO2_DATA_ROOT`, `VSMARTMOM_SIF_DATA_DIR`).
Scripts that read workflow configs/data are updated to these locations.

### 1.5 `Δ_angle` retired permanently

The δBGE forward-peak exclusion angle is redundant: the fit is effectively
bounded by the smallest angle of the quadrature that represents the phase
function, so a separate exclusion cone is at best redundant and at worst
error-prone (the linearized path never honoured it, so non-zero values made
Jacobians inconsistent with the forward model). Its numerical effect is
already retired in this branch (the δBGE constructor forces it to 0 with a
one-time warning). Still to do: remove remaining user-facing documentation
and tests that treat it as meaningful, while keeping compatibility metadata
(the YAML key, the `δBGE(l, Δ)` signature, the stored field) accepted and
ignored, so existing scripts keep running.

### 1.6 `Aerosols.compute_optical_properties` placeholder

The old TOMAS-15→optics placeholder returned made-up numbers (approximate
Mie efficiencies, fixed g = 0.7, a rough number density). The release raises
an explicit "not implemented" error; production Mie optics live in
`Scattering`. Remnant tests/examples that relied on the placeholder are
removed or rewritten.

### 1.7 Other

- Release worktree: `/home/sanghavi/code/github/vSmartMOM-v2.2-release`,
  local commits only; push/PR/tag/registration need explicit approval.
- HITRAN: citation only for now (no change to data hosting or the HITRAN2016
  default).
- TMS (Nakajima–Tanaka) GPU efficiency: deferred to Christian.
- Compute: free A100 (GPU 1 on curry) may be used.
- Review gate before every commit: focused tests → independent Claude review
  → Codex review (`gpt-6.1-sol`, xhigh) → fixes → re-review.

## 2. Explanations kept for later decisions

### 2.1 Why the Round 6/7 runners are pinned to commit `7acab570`

RRS_XCO2 Round 7 is complete (80 scenes × 11 members = 880 retrievals,
finished 2026-09-28 on curry and wurst, no failures; deployment tarball in
`bottom_layer_XCO2_retrievals/round7_fixed_sif_imperfect_correction/deployment/`).

The runners (`Round7RetrievalCampaign.jl`, `round6_environment.sh`) refuse to
start unless the loaded vSmartMOM checkout is exactly `7acab570`, its tree is
clean, and its file-tree hash matches. This is a reproducibility lock: the
campaign identity hash written into every output file then proves which code
produced it, and nobody can accidentally extend Round 7 with different code.

**When this matters to you:** only if Round 7 is re-run or extended. Keep the
clean checkout `/home/sanghavi/code/github/uni_vSmartMOM_round5_jacobian`
(at `7acab570`) as the archive environment. New campaigns get their own pin
(e.g. Round 7_altgrd pins its own commit; later rounds would pin v2.2.0).

### 2.2 Surface-pressure Jacobian: what the code actually does

Both the older core lineage (`suniti_multi_sensor`, `f14371a8` — dated
2026-10-09, but branched from older core code) and the Round 5+
core (`7acab570`) use the same vertical convention: on the *final* model
grid, all interface pressures except the bottom one are held fixed, together
with T, humidity and VMR on that grid, and only `p_half[end] = p_surf`
moves. (The analytic tangent does not include how profile reduction or
observer-interface insertion would move other interfaces if the model were
rebuilt from perturbed *input* parameters; a finite-difference check must
perturb the final grid, not the input file.) Upper layers are not entirely
inert: the Rayleigh optical depth is redistributed by column fraction, and
aerosol vertical allocations are renormalized over all layers (for a profile
defined in altitude all interfaces above BOA move), so upper layers carry
small nonzero tangents. Everything above the bottom layer keeps its pressure levels and
thickness; only the bottom layer gains or loses air (its column, and its
thickness Δz). Neither version holds levels fixed in altitude.

What changed in commit `36e7d1e0` (2026-09-08, "harden v2.2 pressure, source
and batch contracts"): the bottom layer's gas *cross-section* σ(p, T) now
also responds to its shifted mid-layer pressure (dσ/dp via the LUT/line
model, `absorption_pressure.jl`). Without it, the Jacobian is not the
derivative of the model's own forward calculation. The CIA and MT_CKD
density factors were already present in the older core.

Timeline: Round 4 outputs are dated 2026-09-07 and the Round 3-style
`retrievals_acos_*` outputs 2026-09-05, i.e. before this term existed;
Rounds 5–7 ran on `7acab570`, which includes it. (To be confirmed against
each round's provenance — task in progress.)

**Fixed-altitude grid (new):** for a fixed-height flux tower or a
fixed-altitude airborne observer the layering must not move around the
observer when p_surf changes. A `FixedAltitudeGrid()` convention (levels
fixed in altitude; a p_surf change propagates hydrostatically to every
level) is being implemented alongside the current `FixedPressureGrid()`
default, and Round 7 is being repeated as a pilot ("Round 7_altgrd": one
noiseless and one noisy member per scene). Whether Rounds 3–4 (and possibly
all rounds) need re-running is decided after comparing the pilot with
Round 7.

## 3. Open items / TODO

- [ ] **Polarizability tensors for more species.** Populate ᾱ and γ (and
  their wavelength dependence) for as many species as possible — at least
  Ar, CO₂, He, H₂O, CH₄ — for use in Rayleigh (`UsingPolTen`) and Raman
  computations, and add intercomparisons between `UsingPolTen` and
  `UsingRefInd` for every species that has both.
- [ ] ExoOptics/vSmartMOM name unification (audit in progress).
- [ ] Switch ExoOptics production paths to the embedded solar direction;
  quantify embedded vs external differences on the A–D case first.
- [ ] Fix ExoOptics provenance: `ExoEarthRunner.jl:946` records and requires
  a clean `../uni_vSmartMOM`, which is not the vSmartMOM actually loaded.
  Use the loaded package's own path.
- [ ] Move remaining workflows to `sandbox/workflows/` and update the scripts
  that reference old paths (e.g. `generate_truth_map.jl:42`).
- [ ] Fixed-altitude grid + Round 7_altgrd pilot (in progress); decide on
  re-running Rounds 3–4 / all rounds afterwards.
- [ ] TMS GPU efficiency — Christian.
- [x] **For Suniti (convention, documented 2026-10-17):** `conventions.md` §1
  said Q > 0 means polarization parallel to the meridian plane. The code's
  convention (calibrated against its own Rayleigh phase matrix and VLIDORT)
  is Q = |E_φ|² − |E_θ|², U = 2Re(E_φE_θ*), V = 2Im(E_φE_θ*), i.e. Q > 0
  *perpendicular* to the meridian plane. The page now states this; please
  confirm the wording.
- [ ] **For Suniti (surface definition):** vSmartMOM's Ross–Li caps both zenith
  angles at 80° inside the kernels; VLIDORT does not (its Li-sparse kernel
  turns strongly negative at grazing quadrature nodes). The cap changes TOA
  intensity by 0.35–0.84 % in the VLIDORT comparison (BOA by up to 1.1 %).
  Keep the cap (physical, non-negative BRDF) or follow the MODIS definition?
- [ ] ExoOptics audit item 16 (`_read_hitran2024`): delegate to
  AtmosphericAbsorption's `load_lines(HitranPort(path))`, not to vSmartMOM's
  legacy `read_hitran` (the standalone line engine is deprecated in 2.2).
- [ ] Generic analytic surfaces (`analytic_surface.jl`, RPV and Ross–Li) pass
  the reflected (row) direction as `μᵢ`; harmless for these reciprocal scalar
  BRDFs, but fix the order and docstring in a separate commit.
- [ ] Extremely grazing Float32 Cox–Munk inputs can overflow the
  1/(4μᵢμᵣcos⁴β) weight (pre-existing; reviewer finding) — add a tail guard.
- [ ] Ocean/glint products made before the Cox–Munk corrections (B1 slope
  variance, B2 shadowing, B3 polarization) are stale, e.g. the outputs of
  `sandbox/workflows/RRS_XCO2/scripts/generate_truth_map_land_glint.jl`;
  regenerate before use.
