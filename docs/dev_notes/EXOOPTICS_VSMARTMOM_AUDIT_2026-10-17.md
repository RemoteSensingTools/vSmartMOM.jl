# ExoOptics ↔ vSmartMOM audit (2026-10-17)

Read-only audit for the v2.2 consolidation. Nothing was modified and nothing
was executed; findings come from reading the code. "(uncertain)" or "not
verified" marks inference.

**Rule being applied (Suniti, 2026-10-17):** methods live in their most
general form in vSmartMOM, and ExoOptics derives from them.
- Where ExoOptics duplicates a vSmartMOM method under another name, rename it
  to the vSmartMOM name or delegate to it.
- Where no exact parallel exists, add a vSmartMOM method via multiple
  dispatch.
- Embedded solar direction (`external_solar=false`) is the default
  everywhere. The external mode is kept for quick-and-dirty analyses only.

Trees audited:
- **ExoOptics:** `/home/sanghavi/code/github/ExoOptics.jl`, branch
  `suniti_HWO`. The working tree is dirty: `Manifest.toml`,
  `generate_uvrayleigh_horak_dataset.jl` and `SixPhaseInversion.jl` are
  modified, and `exo_degeneracy/nuv_revision/` is untracked.
- **REL:** `vSmartMOM-v2.2-release/src`, the release candidate.
- **HWO:** `run_worktrees/uni_vSmartMOM_suniti_HWO/src`, the vSmartMOM that
  ExoOptics currently loads (`c26c85c3`, dirty).

ExoOptics paths are relative to its repository root. Abbreviations:
- `1D` = `exo_Earth/src/AToD/OneDimensionalDiagnostic.jl`
- `Disk` = `exo_Earth/src/AToD/DiskIntegratedDiagnostic.jl`
- `AFGL` = `exo_Earth/src/AToD/AFGLTropicalProfiles.jl`
- `Runner` = `exo_Earth/src/ExoEarthRunner.jl`
- `PPRT` = `src/ForwardModel/PerPixelRT.jl`
- `UVX` = `exo_degeneracy/exposure/src/UVExposureTimes.jl`
- `DR` = `examples/diagnostics/disk_resolved/scripts`

## Headline findings

1. **The two vSmartMOM trees have opposite solar-direction defaults.**
   - HWO defaults to `external_solar=true` (HWO
     `CoreRT/tools/model_from_parameters.jl:253`). In external mode, any
     non-TOA call throws: `rt_run`, `rt_run_ss`, and the `rt_run_streams`
     callback (HWO `rt_run.jl:333-345`).
   - REL defaults to `false` (REL `model_from_parameters.jl:257`,
     `rt_set_streams.jl:45/201`).
   - On switching to REL, every ExoOptics path that relies on the default
     becomes embedded with no edit.
2. **REL already provides `remake_geometry(model, params; sza, vza, vaz)`**
   (REL `model_from_parameters.jl:703`, exported). It reuses the optics,
   rebuilds only the quadrature and geometry, and builds the embedded form.
   - It replaces `_with_external_solar_sfi` / `_model_at_geometry`
     (Runner:703/768, Disk:366/393).
   - It also replaces the per-node `model_from_parameters` rebuilds in about
     12 DR scripts.
3. **`run_forward`'s production stellar path interpolates across μ₀.**
   - It uses `PPRT._per_pixel_from_streams` (`rt_run_streams` plus bilinear
     (μᵥ, μ₀) interpolation). This is the only production ExoOptics path that
     does so.
   - Measured against exact per-pixel SFI (`SFI_VS_DNI.md:40-57,84-94`):
     about 8.6–9% low DOLP at α=90° with nstreams=5, and +34% to +83%
     intensity-shape error at α ≥ 135°. REVIEW_NEEDED C6 reports DOLP
     converged by nstreams≈13 at its test geometry.
4. **Provenance points at the wrong checkout.**
   - Runner:946 records and requires a clean `../uni_vSmartMOM`
     (`suniti_multi_sensor` @ `f14371a8`, dirty). The loaded tree is HWO.
   - Both are worktrees of one repository, so the recorded hash is a real
     commit from the wrong branch.
   - The committed `Manifest.toml` still names `…/uni_vSmartMOM`; only the
     uncommitted copy points at HWO.
5. **The Stokes rotation matrix is duplicated.**
   `hovenier_rotation(c2β, s2β)` (`src/Geometry/StokesRotation.jl:227`) is the
   same matrix, with the same sign convention, as REL
   `CoreRT.stokes_rotation_matrix(φ, n)` (`CoreRT/Surfaces/fresnel.jl:97`).
6. **Physical constants disagree.**
   - amu: 1.66053907e-27 (PPRT:90), 1.66053906660e-27 (AFGL:20),
     1.66053873e-27 (REL `Absorption/constants/constants.jl:8`).
   - k_B: 1.380649e-23 (1D:53, UVX:29, REL `cia.jl:17`) vs 1.3806503e-23
     (REL `constants.jl:14`).
   - c₂: 1.4387752 vs 1.4387769.
   - Dry-air molar mass: 28.96546 (`src/Inputs/PhysicalInputs.jl:113`) vs
     28.9644 (REL `atmo_prof.jl:48`).
7. **Side bug.** `CrossSectionGrid.sigma_at` (`src/Opacity/CrossSectionGrid.jl:230`)
   interpolates linearly in P, although its docstring says log10 P.

## Task 1 — Duplicates

37 items: 7 rename/delegate (a), 22 new vSmartMOM dispatch methods (b),
8 keep in ExoOptics (c).

Closeness codes:
- **I** = identical
- **C** = same physics, different convention or units
- **S+** = ExoOptics is a superset
- **N** = no vSmartMOM counterpart

### 1A. Rayleigh, depolarization, Greek coefficients, phase and rotation

| # | ExoOptics | vSmartMOM counterpart | Close | Rec. | Call sites |
|---|---|---|---|---|---|
| 1 | `hovenier_rotation` `src/Geometry/StokesRotation.jl:227` | REL `CoreRT.stokes_rotation_matrix(φ, n)` `fresnel.jl:97` | I | (b) Move it to Scattering, export it, and add `stokes_rotation_matrix(cos2φ, sin2φ, n)`. Rename in ExoOptics and keep a deprecated alias. | DiskIntegrate.jl:115, HorakDiskIntegrate.jl:82, Runner:292, Disk:355, DR/pure_gas_horak_observables.jl:129, DR/generate_uvrayleigh_horak_dataset.jl:293, exports, test_stokes_rotation.jl |
| 2 | `stokes_rotation_angles`, `cos_view_angle`, `cos_solar_zenith`, `relative_azimuth` | — | N | (c) Disk-pixel geometry | — |
| 3 | `RAYLEIGH_METHOD = RefractiveIndexRayleigh()`, `rayleigh_tables` 1D:28,57,222-238; audit_ri script | HWO `molecular_rayleigh.jl` | I | (a) Use `UsingRefInd()` after the port | 1D:28-29,57,228,293 |
| 4 | `_rayleigh_extrapolation_masks` 1D:240-250 | HWO `molecular_rayleigh_provenance` | I | (a) Derive the masks from the provenance ranges | 1D ~943 |
| 5 | `rayleigh_optics` 1D:253-283 | REL `getRayleighLayerOptProp` (Bodhaine, air only) | S+ | (b) `rayleigh_layer_optical_depth(method, λ, columns)` | 1D:949 |
| 6 | `_mixed_rayleigh_greek` 1D:286-316 | REL `get_greek_rayleigh(depol)` | S+ | (b) `get_greek_rayleigh(method, composition, λ)` and `mix_greek` | 1D:790, Disk:640 |
| 7 | `_copy_greek_coefficients!` and in-place τ copies 1D:700-767 | — | N | (b) Model hook such as `set_rayleigh!(model, i_band; greek, τ_rayl)` (check whether `update_model!` already serves) | 1D:744,766-767 |
| 8 | `_tau400_diagnostic` Runner:776-789; audit_nuv `air_properties` | Inlined in REL `model_from_parameters.jl:333-341` | I | (b) `air_rayleigh_depolarization(λ, T)` and `depolarization_from_γ` | Runner:1049, audit_nuv:21 |
| 9 | Inverse of `get_greek_rayleigh`, Runner:805-810; audit_nuv `coefficient_depol` | — | I | (b) `rayleigh_depolarization(g::GreekCoefs)` | Runner:808, audit_nuv:31,87 |
| 10 | `phase_matrix(gas, λ)` in DR/export_pure_gas_rayleigh_phase_matrices.jl:22 | `molecular_rayleigh_properties` → `get_greek_rayleigh` → `reconstruct_phase` | I | (a) Call it explicitly with `UsingRefInd()` | :23; DR/pure_gas_horak_observables.jl:115,118,232,233 |

### 1B. Absorption, cross sections, HITRAN, CIA

| # | ExoOptics | vSmartMOM counterpart | Close | Rec. | Call sites |
|---|---|---|---|---|---|
| 11 | `CrossSectionGrid`, `sigma_at(_logP)`, `sigma_mix`, `interp_weights` (`src/Opacity/CrossSectionGrid.jl`); `_interp_sigma_to_ν` | `InterpolationModel`, `_layer_absorption_cross_section` | C/S+ | (b) `MultiPerturberXSTable <: AbstractCrossSectionModel`; ExoMol I/O stays in ExoOptics. Fix the log-P bug. | ExoMolBackend.jl:121,134; exports; test_phase_f1/f3 |
| 12 | `fill_tau_abs!` PPRT:76-129; `_fill_monomer_opacity!` 1D:486-593 | REL `compute_absorption_profile!` | C | (b) A method taking an explicit column vector | PPRT:305,431,517; 1D:966; tests |
| 13 | Hydrostatic column Δp/(g·m̄): PPRT:90-101, AFGL:403-410 | REL `compute_atmos_profile_fields` (fixed g₀ and dry mass) | S+ | (b) Thread `g₀` and `dry_molar_mass` through the parameters | PPRT:83-101,258; AFGL:405-410; build_profiles.jl |
| 14 | `molecular_mass`, `GasMixture`, `mean_molar_mass` | REL `mol_weight` (per isotopologue) | S+ | (b) Move into vSmartMOM and reconcile the air value | PPRT:161; AFGL; 1D; build_profiles.jl |
| 15 | `specific_humidity_from_vmr` Runner:196; 1D `_specific_humidity` | Inverse of REL `atmo_prof.jl:66` | I/C | (b) `CoreRT.specific_humidity_from_vmr` | Runner:748; 1D:730 |
| 16 | `_read_hitran2024`, `_decode_hitran_local_isotope` 1D:448-484 | REL `read_hitran` already decodes 0/A..Z | I | (a) Delete; call `read_hitran` | 1D:532 |
| 17 | Electronic cross sections (SG14 O₃, HITRAN xsc; 1D:318-394 and build_uv_opacity_cache.jl, near-copies) | — | N | (b) `read_hitran_xsc` and `ElectronicCrossSectionModel` | 1D, build_uv_opacity_cache.jl, Runner:895,1013 |
| 18 | `cia_density_column_product` 1D:596-614 | REL `compute_τ_cia!` (Earth columns) | C | (b) Method with explicit densities and columns | 1D:670 |
| 19 | `_fill_cia_opacity!` (family stitching) 1D:622-691 | REL `load_cia_table` | S+ | (b) A `reference_families` keyword | 1D:979 |
| 20 | `validate_opacity_bounds` Runner:231-246 | REL returns zero out of range | S+ | (b) A `bounds=:zero \| :error` keyword | tests |

### 1C. Stellar spectra, Planck, constants

| # | ExoOptics | vSmartMOM counterpart | Close | Rec. |
|---|---|---|---|---|
| 21 | `flux_scaling` (GravitationalDarkening.jl:282) | REL `planck_spectrum_wn` | C | (b) Scalar expm1 method plus exported c₂ |
| 22 | `planck_lambda` UVX:159 | REL `planck_spectrum_wl` | C | (b) Owner's call: UVX is documented as standalone |
| 23 | Photon conversion UVX:171 | REL `watts_to_photons` (not exported) | S+ | (b) Generalise and export |
| 24 | `SolarSpectrum`, `load_solar_spectrum` (E490) | REL default solar spectrum | N/C | (b, low priority) |
| 25 | Scattered constants (amu, k_B, h, c, c₂) | Scattered and inconsistent in REL | I/C | (b) An exported `vSmartMOM.PhysicalConstants` (CODATA 2018) |

### 1D. Quadrature, geometry, model construction, profiles, ILS

| # | ExoOptics | vSmartMOM counterpart | Close | Rec. |
|---|---|---|---|---|
| 26 | `_with_external_solar_sfi`, `_model_at_geometry` (Runner:703-774; Disk:366-400) | REL `remake_geometry` | I | (a) Use `remake_geometry` |
| 27 | Per-node `deepcopy(params)` plus full rebuild in about 12 DR scripts, test_phase_v3, PPRT:430 | REL `remake_geometry` | I | (a) Build once, then `remake_geometry` |
| 28 | `ExoParams.T_layers/p_levels/ν_grid` duplicating the model's profile and bands | `model.profile`, `get_spec_bands` | I | (a) Read from the built model |
| 29 | `remap_pressure_column_mean` AFGL:284 | REL `reduce_profile(_binavg)` | C | (b) `reduce_profile(p_half_target, profile, ::ColumnMeanRemap)` |
| 30 | `convolve_constant_R`, `convolved_wavelengths` (ExoEarthLUT.jl) | — | N | (b), borderline (c) |
| 31 | Horak node loops (Runner, Disk, DR scripts) | — | intra-ExoOptics | (c) Consolidate onto `horak_disk_integrate` |
| 32 | `_stream_interp_weights` PPRT:663 and four script copies | — | intra-ExoOptics | (c) One shared helper, quick mode only |
| 33 | Disk, Fibonacci and Horak quadratures, spheroid, GD models | — | N | (c) |
| 34 | ExoMol I/O, Dummy backends | — | N | (c) |
| 35 | Campaign profiles and grids (AFGL, Runner, ExoEarthLUT, 1D) | — | N | (c) SHA-pinned study inputs |
| 36 | Exposure, statistics, retrieval utilities (UVX, SixPhaseInversion, two solar-twin copies) | — | N | (c) Merge the solar-twin copies |
| 37 | Git provenance helpers (1D, Disk, Runner, scripts) | — | intra-ExoOptics | (c) Consolidate (see Task 3) |

Already delegating correctly: `build_thermal_emission` (`ThermalEmission`),
`build_solar_beam` (`SolarBeam` / `BlackbodySource`), and UV CIA
(`load_cia_table`, `compute_τ_cia!`).

## Task 2 — Solar-direction inventory

Baseline facts:
- **Same TOA output in both modes.** Embedded `first(rt_run(m))` has the same
  [nVZA × nStokes × nSpec] shape as `rt_run_toa(m)`.
- **Parity in Float64.** REL `test_external_solar_sfi.jl:118-120`: the two
  agree at rtol 2e-13 (Float64, CPU). Float32/GPU is not verified.
- **Cost of the embedded operator.** 5 streams plus 1 VZA in IQU goes from
  18×18 to 21×21, about 1.6× in the dense steps (estimate).
- **Extra work in embedded `rt_run`.** It also post-processes BOA, T, HDR and
  BHR.
- **Optional vSmartMOM change:** allow the TOA-only path on embedded models.
  This needs a parity test.

| # | Location | Purpose | Embedded replacement | Numerical change |
|---|---|---|---|---|
| T1 | Runner `_with_external_solar_sfi` :703-731, `_model_at_geometry` :768-774, `rt_run_toa` :1030; production_base.yaml; tests | **Production** Exo-Earth | `remake_geometry` per Horak node, then `first(rt_run(...))`; shapes unchanged; update the Nquad/operator asserts and attributes | Roundoff |
| T2 | PPRT `_run_thermal_only_bulk` :304-312 | **Production** (thermal-only) | Drop the branch; embedded automatically on REL | Roundoff |
| T3 | PPRT `_per_pixel_stokes_cube` :420-444; test_phase_v2/v3 | Validation | `remake_geometry` plus `rt_run` | Roundoff |
| T4 | PPRT `_per_pixel_from_streams` :488-653, from `run_forward` | **Production** stellar | Exact per-pixel, SZA buckets, or Horak nodes; keep streams behind an explicit `method=:streams_interp` | **Yes** (see headline 3) |
| T5 | 1D :739 `external_solar=true`, :803 `rt_run_toa` | Diagnostic (A–D 1-D) | Drop the keyword; asserts become 7/21 | Roundoff (Float32 not verified) |
| T6 | Disk :366-400, :657 | Diagnostic (A–D disk) | `remake_geometry` plus `rt_run`; validator 7/21 | Roundoff |
| T7 | DR/pure_gas_horak_observables.jl:139-156 | Diagnostic | Decision 2026-10-17: keep external, adding `external_solar=true` (quick analysis) | — |
| T8 | About 9 DR scripts using the default model with `rt_run` | Diagnostics | Become embedded on REL with no edit (they throw under HWO) | None vs embedded |
| T9 | DR build_GL/build_uniform/scan_* (`rt_run_streams` plus μ₀ interpolation) | Quick-mode diagnostics | Keep; pass `external_solar=false` explicitly | Interpolation error |
| T10 | generate_uvrayleigh_horak_dataset.jl, audit_nuv_depolarization.jl | Manuscript dataset and audit (already embedded) | None | None |

## Task 3 — Provenance path

- **Wrong:** Runner:945-946 uses `../uni_vSmartMOM`. Minimal fix:
  `repo_root = pkgdir(ExoOptics)` and `vsmart_root = pkgdir(VM)`. Also
  record `pkgversion(VM)` and the path, and fail clearly when there is no
  `.git` (a registry install).
- **Stale after the switch:**
  - `exo_degeneracy/nuv_revision/launch_thin_horak.py:180` hard-codes the
    HWO path.
  - `audit_ri_rayleigh_cross_sections.jl:49-51` hard-codes the
    `Scattering/molecular_rayleigh.jl` file name; keep that path in the port.
- **Correct:** 1D:1245, Disk:841 and audit_nuv:83 already use
  `pathof(vSmartMOM)`.
- **Uncommitted edit:** `Manifest.toml` needs committing or retargeting.
- **Documentation only:** STATUS.md, SESSION_2026-05-19.md and
  exo_Earth/README.md have stale `../uni_vSmartMOM` links.
