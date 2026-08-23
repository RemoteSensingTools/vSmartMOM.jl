# TMS Discovery Audit (memo §1)

Audit of `feat/surface-split` (post multisensor merge) against the approved
"In-line TMS Single-Scattering Correction" memo. Paths/lines at commit time
of the TMS implementation.

1. **SFI beam source assembly.** `elemental!`
   (`src/CoreRT/CoreKernel/elemental.jl`): documented separable beam term
   `j₀± = (ϖ·Z(μ,μ₀)/4μ)·F₀·e^{−τ̄/μ₀}` per Fourier mode; the MS term is
   driven separately by the diffuse field through doubling/interaction.
   The per-layer beam term at a view row, taken right after `elemental!`
   and before `doubling!`, is exactly the truncated single-scattering
   source — the masking point used here (`rt_kernel.jl`, noRS method).

2. **Zero-weight view nodes.** `rt_set_streams`
   (`src/CoreRT/tools/rt_set_streams.jl:132,155`) appends VZA cosines with
   weight 0; robust identification is `wt_μN .== 0` (Stokes-expanded). In
   the legacy embedded representation the SZA node is also zero-weight —
   masking it too is correct (its beam source is only read back at its own
   output, i.e. only when it *is* a view). Radau keeps μ₀ weighted →
   unmasked. External-solar mode has no μ₀ node at all; the same weight
   test selects exactly the VZAs.

3. **δ-M truncation.** `f` per scatterer: `AerosolOptics.fᵗ` from
   `truncate_phase` (δBGE, full angular domain — Δ_angle retired). τ/ω
   scaling: `createAero` / `build_m_invariant_cache` pre-apply
   `(1−fᵗω̃)`. **No per-layer effective f exists — and none was added:**
   since `Δτ′ = Δτ(1−ωf)`, the N–T factor `ω/(1−ωf)` equals
   `β_scat/Δτ′`, so the accumulation uses the solver's own scaled
   `τ_sum_all` with UNSCALED scattering coefficients
   (`single_scattering.jl`, module header). One less derived quantity
   that could drift from the optics.

4. **Untruncated scattering matrix.** Retained as
   `AerosolOptics.phase_greek_raw` (node-resolved raw Greek sets, attached
   by `truncate_phase` and threaded through
   `_spectralize_truncated_endpoints`). Meridian rotation lives in
   `Scattering.phase_matrix_first_column` /
   `_rotation_from_scattering_plane`
   (`src/Scattering/analytic_phase_functions.jl`) with the memo's
   upwelling convention `cosΘ = −μ_v μ₀ + √(1−μ_v²)√(1−μ₀²)cosΔφ` —
   reused, not re-derived. For an unpolarized beam the incident-side
   rotation is the identity on `[1,0,0,0]`, so the first Z column fully
   determines the SS source (per Sanghavi: no Fourier moments are needed
   once zenith angles are replaced by the scattering angle).

5. **Beam transmittance convention.** Plane-parallel everywhere (no
   Chapman factors in the tree). `τ_sum_all` (δ-scaled cumulative τ per
   layer boundary) is computed once per run in `extractEffectiveProps` and
   was already saved "for TMS correction" (`rt_run.jl`). The accumulation
   consumes it directly.

6. **`ExactSFIPhase` (v0.5).** Does not exist on this branch. Prior art:
   `rt_run_ss` — a Fourier-space solver over *truncated* per-mode optics
   (elemental-only, `interaction_ss!`), i.e. the truncated SS diagnostic,
   NOT an exact reference (a peaked kernel's azimuth expansion converges
   too slowly for that to ever work); and `apply_ss_correction!` — the
   Cox-Munk glint subtract/add surface correction (separate mechanism,
   untouched). The exact reference role is now filled by
   `rt_run_ss_exact` (same engine as TMS, `f = 0` semantics, no Fourier
   loop).

7. **Fourier summation / post-sum assembly.** `postprocessing_vza!`
   accumulates azimuth-weighted moments into `R_SFI (nvza, nStokes,
   nSpec)`; the post-sum insertion point (shared with the glint
   correction) is where `tms_correction!` adds the exact SS.

8. **Linearization touchpoints.** `rt_run_lin` **rejects**
   `TMSCorrection` explicitly (parity policy rule 2) — the correction is
   forward-only in this deliverable. The f₁ additions
   (`phase_greek_raw` retention, `_exact_phase_columns`,
   `_ss_path_factor` with `expm1`) are `FT`-generic and Dual-safe.

9. **Branch state.** Memo deviations, approved: (a) no `ExactSFIPhase`
   (see 6); (b) two beam representations exist (embedded-μ₀ and
   external-solar Z₀ columns) — the weight-test mask covers both;
   (c) spectral-node optics — `Z_exact` is evaluated at the stored
   `phase_ν` nodes and linearly interpolated in ν, mirroring
   `phase_greek` consumption.

## Known scope limits (all rejected loudly, per policy)

- Forward-elastic (`noRS`) + SFI only; unpolarized solar source only
  (first-column evaluation).
- Upwelling TOA only: `j₀⁺` (downwelling) keeps the truncated SS, so BOA
  `T` is unchanged from the historical solver.
- Atmosphere/surface split cache builds reject TMS (the post-sum addition
  must move to the replay — follow-up; it will fold into the Lambertian
  closure's `R_black`).
- Interior sensor heights (multisensor) reject TMS.
- IMS (second-order forward-peak term) out of scope; the small-Θ, high-τ
  residual of TMS is the known IMS term (Nakajima & Tanaka 1988).
