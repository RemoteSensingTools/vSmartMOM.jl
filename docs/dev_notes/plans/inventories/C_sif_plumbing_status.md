# Inventory C — SIF Plumbing Status

**As of:** 2026-04-19

**Branch tips:**

| Branch | Worktree | Tip | Subject |
|---|---|---|---|
| `sanghavi` | `/home/sanghavi/code/github/vSmartMOM.jl/` | `9ee9a75` | add plan |
| `unified-vsmartmom` | `/home/sanghavi/code/github/uni_vSmartMOM/` | `a4e4187` | Add batched-kernel and Raman scaling benchmarks with writeup |

**Scope of this inventory:** catalog the state of Solar-Induced Fluorescence (SIF) plumbing on both branches to scope the port into `sanghavi-unified`. This covers `RS_type.SIF₀` field definitions, kernel-side usage, the `src/SIF_emission/` data files, and the benchmark scripts that consume them. Linearized-SIF work is out of scope.

---

## 1. SIF₀ plumbing on `sanghavi`

### 1a. Struct field definitions (`vSmartMOM.jl/src/Inelastic/types.jl`)

`SIF₀::Array{FT,2}` (shape: `(pol_type.n, nSpec)`, Stokes × spectral) appears as a **live, uncommented** field on the following Raman-type structs on `sanghavi`:

| Line | Struct | Status |
|---|---|---|
| 34  | `RRS{FT}` (mutable) | live |
| 63  | `VS_0to1{FT}` | live |
| 91  | `VS_1to0{FT}` | live |
| 128 | `noRS{FT}` (mutable) | live |
| 219 | `VS_0to1_plus{FT}` (mutable) | live, with `= zeros(FT,1,1)` default |
| 272 | `VS_1to0_plus{FT}` (mutable) | live |
| 282 | `noRS_plus{FT}` (mutable) | live |

Commented-out `#SIF₀` declarations (not active) appear on:
- line 309 (`sol_RRS`), 336 (`sol_VS_0to1`), 362 (`sol_VS_1to0`), 433 (second `RVRS` definition block), 475 (`sol_VS_0to1_plus`), 520 (`sol_VS_1to0_plus`).

### 1b. Field also echoed in `stellar_types.jl`

`vSmartMOM.jl/src/Inelastic/stellar_types.jl` lines 32, 59, 85, 156, 198, 243, 265, 282, 299, 318, 342, 366 all show `#SIF₀::Array{FT,2}` — all commented out. Stellar-side plumbing is placeholder-only.

### 1c. Kernel-side use (Lambertian surface only)

`RS_type.SIF₀` is read (not written) only in the Lambertian surface kernel on the sanghavi side:

- `vSmartMOM.jl/src/CoreRT/lambertian_surface.jl`
  - L64, L67–68 (first `create_surface_layer!` for the basic Lambertian case): `added_layer.J₀⁻[:,1,:] .+= (1/π)*repeat(arr_type(RS_type.SIF₀),Nquad) * unweight` — **live**.
  - L157–158 (second variant, `LambertianSurfaceLegendre` or spectrum branch): same injection, **live**.
- `vSmartMOM.jl/src/CoreRT/lambertian_surface_lin.jl`
  - L93–95 (first variant in linearized file): **commented out**.
  - L189–190 (second variant in linearized file): **live**.

No other surface kernel (RPV, RossLi, Cox-Munk, Canopy) references `SIF₀` on sanghavi. Elemental/doubling/interaction kernels never touch `SIF₀` directly — the SIF contribution enters purely through the surface source.

### 1d. Population paths (benchmark scripts only)

`RS_type.SIF₀` is populated exclusively from the benchmark scripts (no src/ loader). Active writers:

- `test/benchmarks/creategrid_O2Aband_RamanSIF.jl` L116–117, 132 (struct construction), L150–151 (resize to `(n, length(ν))`), L165–166 (assignment; currently forced to zero at L162 with live assignment commented out).
- `test/benchmarks/creategrid_O2Bband_RamanSIF.jl` L113–114, 129 (construction), L146–147 (resize), L161–162 (assignment; forced to zero at L158).
- `test/benchmarks/emit_modtran_noRS_scenarios.jl` L709: `RS_type.SIF₀ .= SIF₀_vec` — live.
- `test/benchmarks/prototype_O2ABand_RRS_SIF.jl`, `prototype_O2BBand_RRS_SIF.jl`, `test_creategrid_O2Aband_RamanSIF.jl`, `test_creategrid_noabs_RamanSIF.jl`, `creategrid_betwnAB_RamanSIF.jl`, `creategrid_O2Aband_OCORayl.jl`: all read `sif-spectra.csv` via `readdlm(...)` and follow the same pattern (see §3, §4).

No `src/` function exists that computes or loads `SIF₀` for an RS-type; the plumbing is entirely user-driven.

---

## 2. SIF₀ plumbing on `unified-vsmartmom`

### 2a. Struct field definitions (`uni_vSmartMOM/src/Inelastic/types.jl`)

**Good news:** the `SIF₀` field already exists on unified in the core Raman types, cleaned up to use `AbstractRamanType{FT}`:

| Line | Struct | Default |
|---|---|---|
| 39  | `RRS{FT}` (mutable) | `Array{FT,2}` required |
| 73  | `VS_0to1{FT}` | required |
| 106 | `VS_1to0{FT}` | required |
| 150 | `noRS{FT}` (mutable) | `= zeros(Float64, 1, 1)` |
| 248 | `VS_0to1_plus{FT}` | required |
| 305 | `VS_1to0_plus{FT}` | required |
| 321 | `noRS_plus{FT}` | required |

Commented-out `#SIF₀::Array{FT,2}` on lines 357, 387, 416, 487, 532, 580 mirrors sanghavi (sol_* / second `RVRS` block).

### 2b. `stellar_types.jl` — same as sanghavi

All `SIF₀` declarations are commented out (lines 32, 59, 85, 156, 198, 243). No behavioral change from sanghavi.

### 2c. Kernel-side use — MISSING

Grep for `SIF` across `uni_vSmartMOM/src/CoreRT/` returns only **one match**, and it is commented out:

- `uni_vSmartMOM/src/CoreRT/Surfaces/lambertian_surface_lin.jl` L121–123:
  ```
  # for SIF
  if ...
      #added_layer.j₀⁻[:,1,:] .+= (1/π)*repeat(arr_type(RS_type.SIF₀),Nquad) * unweight
  ```

The non-linearized `uni_vSmartMOM/src/CoreRT/Surfaces/lambertian_surface.jl` has **no SIF reference whatsoever** — the surface kernel writes `j₀⁻[:,1,:] .= μ₀*(R_surf*I₀_NquadN) .* exp.(-τ_sum/μ₀)'` at L55 and never adds the `(1/π) * SIF₀` isotropic emission term.

No other surface file (`rpv_surface.jl`, `rossli_surface.jl`, `coxmunk_surface.jl`, `canopy_surface.jl`, `rpv_rahman.jl`) references SIF on either branch.

### 2d. Population — MISSING from CSV

`RS_type.SIF₀` is initialized to zeros and never populated from any data file on unified:

- `uni_vSmartMOM/test/test_forward_raman.jl` L60 / L77: `SIF₀ = zeros(FT, nPol, nSpec)` → passed into `RRS(...)` — stays zero.
- `uni_vSmartMOM/test/test_forward_raman_gpu.jl`: same pattern.
- No grep hit for `sif-spectra`, `SIF_emission`, or a SIF loader function in `uni_vSmartMOM/src/`.

### 2e. Summary: unified has the **field, not the behavior**

| Piece | sanghavi | unified | Parity gap |
|---|---|---|---|
| `SIF₀` on `RRS`, `VS_0to1`, `VS_1to0`, `noRS` | yes (live) | yes (live) | — |
| `SIF₀` on `*_plus` variants | yes (live) | yes (live) | — |
| Lambertian surface injection (non-lin) | yes, L67–68 & L157–158 | **missing** (no SIF code at all) | need to port |
| Lambertian surface injection (lin) | partial (L93 commented, L189 live) | commented out at L123 | needs review + re-enable |
| Other BRDF surfaces: SIF emission term | none | none | ask: should RPV/RossLi/Canopy also emit SIF? |
| CSV loader helper in `src/` | none | none | candidate for new `src/SIF_emission/loader.jl` |
| `sif-spectra.csv` committed | no (only `.dat~` backups) | no (no `SIF_emission/` dir at all) | must be committed |

---

## 3. SIF data files — commit plan

### 3a. Files on disk (sanghavi worktree)

Under `vSmartMOM.jl/src/SIF_emission/`:

| File | Size | Purpose | Git status |
|---|---|---|---|
| `sif-spectra.csv` | 13 695 B | SIF emission spectrum (λ, J_SIF), the primary driver used by every `*RamanSIF*` benchmark | **UNTRACKED** |
| `PC1_SIFSpectra_allSpecies.csv` | 8 288 B | Alternate SIF PCA basis, all benchmarks have it as a commented-out alternative to `sif-spectra.csv` | **UNTRACKED** |
| `ficus_refl.dat` | 35 161 B | Ficus leaf reflectance (full range) for the Canopy BRDF | **UNTRACKED** (only `.dat~` backup is tracked) |
| `ficus_refl_600to800nm.dat` | 3 215 B | Ficus leaf reflectance trimmed to O2A band | **UNTRACKED** (only `.dat~` backup is tracked) |
| `ficus_refl.dat~` | 0 B | backup | tracked (empty) |
| `ficus_refl_600to800nm.dat~` | 35 161 B | backup | tracked |

Note the oddity: on sanghavi the `.dat~` backups are tracked and they are actually wrong (`ficus_refl.dat~` is empty; `ficus_refl_600to800nm.dat~` holds what should be in `ficus_refl.dat`).

### 3b. Files on disk (unified worktree)

`uni_vSmartMOM/src/SIF_emission/` does not exist. Neither the directory nor the data files are anywhere on the branch.

### 3c. Recommended canonical home for the merged branch

All four active files are small (≤ 35 KB). **None need Git LFS**; all can be committed directly to git.

Recommendation for `sanghavi-unified`:

1. `src/SIF_emission/sif-spectra.csv` — commit directly. Required by all `*RamanSIF*` benchmarks.
2. `src/SIF_emission/PC1_SIFSpectra_allSpecies.csv` — commit directly. Alternate basis, consistently referenced (commented-out) in benchmarks; users may swap.
3. `src/SIF_emission/ficus_refl.dat` — commit directly (Canopy BRDF depends on this). Note: the tracked `.dat~` backup is empty; the real data is in the untracked file.
4. `src/SIF_emission/ficus_refl_600to800nm.dat` — commit directly (O2A-band Canopy use).
5. Drop both `.dat~` backup files from the tree (they are misnamed backups of other files).
6. Add a one-line `src/SIF_emission/README` or header comment describing provenance, units (`J_SIF` in mW/m²/nm per nm column), and the columns `λ (nm), J_SIF`.

### 3d. Path-resolution cleanup

Every benchmark references these files via **hardcoded absolute paths** like `/home/sanghavi/code/github/vSmartMOM.jl/src/SIF_emission/sif-spectra.csv`. No `@__DIR__` / `joinpath` / `pkgdir` resolution anywhere in `src/` or in the benchmarks. As part of the port, these should be rewritten to `joinpath(pkgdir(vSmartMOM), "src", "SIF_emission", "sif-spectra.csv")` or similar, so the scripts work on any machine.

---

## 4. SIF benchmark scripts

### 4a. `test/benchmarks/creategrid_O2Aband_RamanSIF.jl` (927 lines)

- **RS types used:** `InelasticScattering.noRS` (elastic reference) and `InelasticScattering.RRS` (full rotational Raman). Both constructed with `SIF₀ = zeros(FT,1,1)` at L116–117 and L132; later resized at L150–151 to `zeros(pol_type.n, length(ν))`.
- **SIF₀ source:** `readdlm(".../src/SIF_emission/sif-spectra.csv", ',')` at L43; wavelengths reversed to wavenumber grid, normalized by `(0.5π / max)`, unit-converted to mW/m²/cm⁻¹, wrapped in a `LinearInterpolation`. At L165–166 the interpolated SIF is assigned to `RS_type0.SIF₀[1,:]` and `RS_type1.SIF₀[1,:]`.
- **Currently:** L162 forces `SIF₀ .= 0.0` so the active code path is a no-SIF run (the real `SIF_interp.(ν)` call is commented out — sanghavi runs it with SIF disabled while keeping the plumbing intact).
- **Solar model:** `solar_transmission_from_file(".../src/SolarModel/solar.out")` → `LinearInterpolation` → multiplied with Planck function at 5777 K to form `F₀`.
- **Output:** `writedlm` → `/home/sanghavi/data/RamanSIFgrid_new/rayl_sza{sza}_vza{vza}_vaz{vaz}_alb{alb}_psurf{psurf}hpa_{nors|rrs}_ABO2.dat`. Sweeps 14 SZAs × 21 albedos × 3 surface pressures × N_vaz × 2 (RRS / no-RRS). Rows: `[ν R_I R_Q R_U ieR_I ieR_Q ieR_U F₀]`.
- **YAML config consumed:** `test/test_parameters/O2_parameters2_SIF_grid_float32.yaml` (L22).
- **Hardcoded absolute paths requiring rewrite:**
  - L22 (YAML config)
  - L43–44 (SIF CSVs)
  - L49 (solar.out)
  - L187–191 (`/home/sanghavi/data/RamanSIFgrid_new/...` output dir — this is user-scratch, should become a CLI arg or env var).

### 4b. `test/benchmarks/creategrid_O2Bband_RamanSIF.jl` (907 lines)

Structurally identical to the A-band script:
- Same noRS/RRS struct construction (L108–129), same SIF₀ resize / assignment pattern (L144–147, 161–162).
- Same `sif-spectra.csv` loader (L39), same solar-model load (L45).
- **YAML config:** `test/test_parameters/O2_parameters1_SIF_grid.yaml` (L19) — B-band variant.
- **Output:** same `/home/sanghavi/data/RamanSIFgrid_new/...` naming but for the B-band.
- Same hardcoded-path cleanup needed (L19, L39–40, L45, and the output dir string).

### 4c. Related SIF/Raman benchmark scripts on sanghavi (for completeness)

`creategrid_betwnAB_RamanSIF.jl`, `test_creategrid_O2Aband_RamanSIF.jl`, `test_creategrid_noabs_RamanSIF.jl`, `prototype_O2ABand_RRS_SIF.jl`, `prototype_O2BBand_RRS_SIF.jl`, `RamanSIFspectra.jl`, `RamanSIFmaps.jl`, `RamanSIFPolyCoeffMaps.jl`, `RamanSIFworldmaps.jl`, `RamanSIF_oco_retr.jl`, `sif_raman.jl`, `testAer_O2Aband_Raman.jl`, `emit_modtran_noRS_scenarios.jl`, `creategrid_O2Aband_OCORayl.jl` all read `sif-spectra.csv` via the same hardcoded absolute path. `emit_modtran_noRS_scenarios.jl` is the **only** script that actually assigns a nonzero SIF (L708–709).

### 4d. YAML configs used by the two scripts

Neither YAML is committed on `unified-vsmartmom`:
- `test/test_parameters/O2_parameters2_SIF_grid_float32.yaml` (A-band) — sanghavi only
- `test/test_parameters/O2_parameters1_SIF_grid.yaml` (B-band) — sanghavi only

The full set of 12 `*SIF*.yaml` configs under `test/test_parameters/` on sanghavi (listed in §1d / git ls-files) has no counterpart on unified and must come over as part of the port.

---

## 5. Gap list for `sanghavi-unified` (port order)

Priority 1 — must have, low risk:

1. **Commit SIF data files** into `src/SIF_emission/`: `sif-spectra.csv`, `PC1_SIFSpectra_allSpecies.csv`, `ficus_refl.dat`, `ficus_refl_600to800nm.dat`. Drop the two tracked `.dat~` backups. Sizes are all ≤ 35 KB; no LFS needed.
2. **Port the Lambertian SIF injection** from `vSmartMOM.jl/src/CoreRT/lambertian_surface.jl` L67–68 and L157–158 into the corresponding spots in `uni_vSmartMOM/src/CoreRT/Surfaces/lambertian_surface.jl` (add `(1/π) * repeat(arr_type(RS_type.SIF₀), Nquad) * unweight` to `j₀⁻[:,1,:]`). This is the single behavioral change needed to make the existing unified `RS_type.SIF₀` field actually do something.
3. **Port the two benchmark scripts** `creategrid_O2Aband_RamanSIF.jl` and `creategrid_O2Bband_RamanSIF.jl` into `uni_vSmartMOM/test/benchmarks/`, **replacing hardcoded absolute paths with `joinpath(pkgdir(vSmartMOM), ...)`** for the YAML config, `sif-spectra.csv`, and `solar.out` inputs. Leave the scratch output directory as a user-supplied CLI/env var.
4. **Port the two YAML configs** `O2_parameters2_SIF_grid_float32.yaml` and `O2_parameters1_SIF_grid.yaml` into `uni_vSmartMOM/test/test_parameters/`.

Priority 2 — useful, moderate scope:

5. **Add a tiny `src/SIF_emission/sif_loader.jl`** exposing something like `load_sif_spectrum(path = default_sif_csv_path()) → (ν, jSIF_interp)` so scripts and tests can share the CSV reader instead of duplicating the `readdlm` / `reverse` / `jSIF .*= 1e7./ν.^2` chain. Export from `vSmartMOM.SolarModel` or a new `vSmartMOM.SIF` submodule.
6. **Port the remaining `*SIF*.yaml` configs** (`O2Parameters1_SIF.yaml`, `O2Parameters2_SIF.yaml`, `O2_parameters1p5_SIF_grid.yaml`, `O2_parameters2_SIF_grid.yaml`, `O2_parameters2_SIF_grid_all_layers.yaml`, `O2_parameters2_SIF_grid_debug.yaml`, `WCO2_parameters_SIF_grid.yaml`, `abs_parameters2_SIF_grid.yaml`, `aerO2_parameters2_SIF_grid.yaml`, `noabs_parameters2_SIF_grid.yaml`) if/when the corresponding benchmark scripts are ported.

Priority 3 — housekeeping:

7. **Revisit `lambertian_surface_lin.jl`** L121–123 on unified: decide whether the commented-out linearized SIF injection should be enabled (out of scope per task, flag for a follow-up).
8. **Document** in `src/SIF_emission/README` the units (`J_SIF` in mW/m²/nm input; converted to mW/m²/cm⁻¹ downstream), the normalization factor `0.5π / maximum(J_SIF)` currently hardcoded in benchmarks, and provenance of the CSV.

Explicitly NOT in scope (per task brief): porting linearized SIF work; implementing the actual physics changes; porting RPV/RossLi/Canopy SIF emission.

---

## 6. Open questions

1. **Is the `0.5π / maximum(J_SIF)` rescaling in the benchmarks intentional normalization or a hack?** It makes the SIF magnitude independent of the CSV's absolute scale, which is unusual for a physics benchmark. Clarify before porting so the merged script does the right thing.
2. **Should non-Lambertian surfaces also emit SIF?** Real canopy SIF is not isotropic and should plug into `canopy_surface.jl`. Neither branch has that; confirm whether the merged branch gets a "SIF only via Lambertian" simplification or a full emission-through-BRDF port.
3. **`noRS` vs `RRS` SIF₀ default:** unified's `noRS` has `SIF₀ = zeros(Float64, 1, 1)` default but `RRS` requires `SIF₀` at construction. Should they be made consistent (both default-initialized) so users can skip SIF arg entirely when it is zero?
4. **Should `SIF₀` move off `RS_type` entirely?** The CLAUDE handoff plan mentions "If SIF is intended to stay outside the `RS_type` system in the long term, decide that now". SIF is a surface emission source, not a Raman property — architecturally it sits awkwardly on `RS_type`. This inventory documents the current layout; a refactor to move `SIF₀` onto the surface struct (or a new `emission` field on `Optics`) is a separate design call.
5. **Tracked `.dat~` files:** are these intentional (some workflow that reads `.dat~`) or leftover mistakes? Confirm before dropping them in the merged branch.
6. **`emit_modtran_noRS_scenarios.jl`:** this is the only script that actually runs with nonzero SIF (`RS_type.SIF₀ .= SIF₀_vec`). Should it be the canonical SIF integration test on `sanghavi-unified` (vs. the currently-SIF-disabled `creategrid_*` scripts)?
