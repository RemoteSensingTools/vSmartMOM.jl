# GCHP → vSmartMOM Aerosol Coupling — Student Brief & Framework

**Audience:** summer student coupling GEOS-Chem High Performance (GCHP / GOCART–TOMAS)
aerosol output into vSmartMOM for radiative-transfer (RT) forward modeling, ultimately to
generate a large ML training set of (atmospheric state → reflected radiances).

**Status of this document:** framework / outlay only. No production code changes are
proposed yet. It maps what already exists, names the design forks you must decide, and
gives a staged roadmap and a standalone-package skeleton you can start from.

**Author context:** assembled from a survey of the `gchp-io` branch, the core RT path on
`RamanLuT`/`main`, and the real file `GEOSChem.Custom.20190702_0000z.nc4` (C24 cubed
sphere, 72 levels, 6 faces, TOMAS-15).

---

## 0. TL;DR — where things stand

- **Reading a GCHP column and turning it into a per-layer, per-bin, speciated aerosol
  state already works.** It lives on the **`gchp-io`** branch (not your current
  `RamanLuT` checkout). Base your work there.
- **Computing column AOD from that state already works** (real Mie extinction, effective
  refractive index via a selectable mixing rule, per-bin integration). This is the
  "use the bins as they are" path.
- **What does NOT exist yet** is the bridge from the sectional aerosol state to the
  *full per-layer scattering optics* (τ, single-scattering albedo ϖ, and the Greek /
  Z-matrix phase coefficients) that the vSmartMOM RT solver actually consumes. The AOD
  code says so itself at the top of the file: *"This does not produce phase matrices or
  wire aerosols into RT."*
- Two **core vSmartMOM enhancements** are implied by the goal: (a) let the **refractive
  index vary across a spectral window** (today only the size parameter varies); (b)
  decide whether the **size-distribution shape** varies with height or only the loading.
- The end goal (≈1M scenes → ML training package) needs an **input/output schema with
  dimensionality reduction** (EOF/SVD on the profiles), designed up front.

---

## 1. Branch reality check (do this first)

| Branch | What's there for aerosols |
|---|---|
| `RamanLuT` (your current checkout) | An **older, mostly-stubbed** `src/Aerosols/`: placeholder Mie, hardcoded asymmetry `g=0.7`, no real mixing rule. **Do not build on this copy.** |
| `gchp-io` | The **real** GCHP coupling work: `GCHPFile`/`GCHPScene` reader, `SectionalAerosolData`, mixing rules, per-bin Mie extinction, column-AOD writers, TOMAS species table. |
| `main` / `RamanLuT` core | The RT solver, `Scattering` (Mie → `AerosolOptics`), `CoreRT/model_from_parameters.jl`. This is shared and is what the new package depends on. |

The two branches have **diverged substantially** in both directions (RamanLuT has RRS/GPU
perf work that gchp-io lacks; gchp-io has the aerosol work RamanLuT lacks). At some point
someone has to decide how to reconcile them — see §7. For the student's work, **branch off
`gchp-io`.**

---

## 2. What the GCHP file gives you (verified against the real `.nc4`)

`GEOSChem.Custom.20190702_0000z.nc4`: cubed-sphere **C24**, dims `(time=1, lev=72, nf=6,
Ydim=24, Xdim=24)`. All aerosol fields are `mol mol⁻¹ dry` mixing ratios, named
`SpeciesConcVV_<SPECIES><bin:02d>` for bins `01..15`:

| Component | Variable prefix | Meaning |
|---|---|---|
| Sulfate | `SF01..SF15` | hydrophilic |
| Sea salt | `SS01..SS15` | hydrophilic |
| Dust | `DUST01..DUST15` | hydrophobic |
| Organic carbon | `OCOB`/`OCIL` `01..15` | hydrophobic / hydrophilic (in-cloud processed) |
| Black carbon | `ECOB`/`ECIL` `01..15` | hydrophobic / hydrophilic |
| Aerosol water | `AW01..AW15` | hygroscopic water per bin → **RH / refractive-index effect** |
| **Number** | `NK01..NK15` | particles per bin (`1000 × particles/mol_air`) |

Plus meteorology you need to convert mixing ratio → number/mass per m³:
`Met_AD` (dry air mass, kg), `Met_AIRVOL` (box volume, m³), `Met_DELP`/`Met_DELPDRY`,
`Met_T`, `Met_SPHU`, `Met_PS*`.

**Three independent things vary down a column** (your framing, confirmed):
1. **Loading profile** — total aerosol concentration vs height (72 levels).
2. **Size distribution** — the 15-bin shape per layer (can change shape, not just scale).
3. **Speciation** — per-bin mass fractions of the components above, each with its own
   complex refractive index. This is what mostly drives the **effective refractive index**.

---

## 3. What already works on `gchp-io` (use it, don't rebuild it)

- **Column reader** — `src/IO/NetCDF/GCHPScene.jl`: `GCHPFile` (open-once) +
  `scene_at(f, ix, iy, iface)` → a `GCHPScene` carrying `p_half, T, q, air_mass,
  air_vol, vmr`, and an optional `aerosols::SectionalAerosolData`.
- **Sectional aerosol state** — `src/Aerosols/abstract_types.jl`: `SectionalAerosolData`,
  `AerosolSizeGrid`, `SpeciesComposition`, plus the policy types
  `{DirectBinSum, ConstantIntegrationPerBin, LinearIntegrationPerBin, LogNormalFit}` and
  `{ExternalMixing, VolumeWeightedMixing, MaxwellGarnettMixing, BruggemanMixing}`.
- **Speciation → effective refractive index** — `effective_ri(comp, db, λ, mixing_rule,
  scheme)` in `src/Aerosols/aod_diagnostics.jl`.
- **Per-bin Mie extinction → column AOD** — `compute_column_aod(...)` (real Mie via
  `Scattering.compute_mie_ab!`), optional Mie-extinction LUT
  (`src/Aerosols/LUT/mie_extinction_lut.jl`).
- **Bulk grid writers** — `write_gchp_aod_diagnostic(...)` and a scene-loop driver
  `generate_benchmark(...)` in `src/IO/Benchmark/`.
- **TOMAS metadata single source of truth** — `data/aerosols/tomas_species.yaml`
  (verified against GCHP `tomas_mod.F90`), loaded by `TOMASScheme(:tomas15)`.
- **Refractive-index database** — `data/refractive_indices_database.yaml` + linear
  interpolation to arbitrary λ (`get_refractive_index`).

**Net:** the entire **"bins-as-they-are → AOD"** path is built. That is the natural
**Stage 1** validation target.

---

## 4. The central design fork: how to represent the aerosol optical properties

This is the part to lay out explicitly for the student. There are **two families**, and the
recommendation is to support both behind one interface, starting with the first.

### Mode A — Sectional (use the bins as they are)
Keep the 15 TOMAS bins; integrate optics within/across them.
- **A1 — piecewise within bin** (already implemented: `DirectBinSum`,
  `ConstantIntegrationPerBin`, `LinearIntegrationPerBin`). Bin shape used as-is.
- **A2 — spline across bin centers** (your "don't discretize the way they are, more like a
  spline"): fit a smooth `dN/dlogD` through the 15 bin values, then quadrature on the
  smooth curve. Reduces the staircase artifact of raw bins; still non-parametric. *New.*

**Pros:** faithful to GCHP microphysics; no fitting assumption; speciation is naturally
per-bin. **Cons:** more Mie evaluations; produces optics only via summation (need the
bridge in §5).

### Mode B — Parametric fit (continuous distribution)
Fit a **lognormal** (or a small sum of lognormals / modal distribution) to the binned
`dN/dlogD` per layer → continuous distribution with a handful of parameters
(median radius, geometric width) — exactly the classic vSmartMOM "parameter file" shape
(`Aerosol(LogNormal(μ, σ), nᵣ, nᵢ)`, see `src/Scattering/types.jl:88`). The
`LogNormalFit()` policy type exists but is a **stub** — this is the implementable piece.

**Pros:** maps directly onto the existing core Mie path; cheap; compact for ML inputs.
**Cons:** an approximation of the true distribution; **speciation is the hard part** — a
single fitted mode has to carry a *mixed* refractive index, so you either (i) fit one mode
per species (multi-modal), or (ii) fit one mode for the size and derive a single effective
RI from the bulk composition. Your instinct — *"maybe the speciation is better done within
the bin"* — is sound: a clean hybrid is **B for the size shape, A for the speciation**
(fit the size continuum but carry composition/RI from the sectional data).

> **Decision to record:** start Mode A1 (validate), add A2 and B behind a common
> `AerosolRepresentation` interface so RT can consume any of them identically.

---

## 5. The missing bridge: sectional state → per-layer RT optics

The RT solver does **not** want AOD. Per band and per layer it wants an `AerosolOptics`
object: extinction `k`, single-scattering albedo `ϖ`, and the **Greek coefficients**
(`α, β, γ, δ, ϵ, ζ`) of the scattering matrix (see how
`compute_aerosol_optical_properties(mie_model)` and `GreekCoefs` are used in
`src/CoreRT/tools/model_from_parameters.jl:447-520`). The per-layer optical depth array
`τ_aer[iBand][iAer, nSpec, iZ]` is already 3-D, i.e. the core already supports
**wavelength-dependent k(λ)** within a band.

**Two ways to build the bridge:**
1. **Per-bin Mie → summed Greek** (pairs with Mode A): Mie each bin at its effective RI,
   weight by number density, sum the phase matrices into one per-layer `AerosolOptics`.
   Most faithful; reuses `Scattering` Mie machinery; new code is the summation/normalization.
2. **Lognormal-fit per layer → core Mie** (pairs with Mode B): emit a per-layer
   `Aerosol(LogNormal, n_eff)` and run the existing `make_mie_model` /
   `compute_aerosol_optical_properties`. Less new code; inherits the fit approximation.

### Per-layer Z-matrices are unavoidable anyway
Even if the aerosol *type* were constant with height, the **effective layer optics still
change per layer** because aerosol optics get mixed with Rayleigh, and the Rayleigh
fraction varies with pressure/height and with wavelength. So "compute the Z-matrix once and
reuse" is not on the table. The real efficiency question is narrower:

> **Does the size-distribution _shape_ vary with height, or only the loading?**
> - *Cheap:* one AOD-weighted mean size distribution for the whole column; scale by the
>   per-layer concentration. Fewer Mie calls.
> - *Faithful:* per-layer size distribution shape (GCHP gives you this directly). More Mie
>   calls, physically correct.

Recommend: build faithful, expose the cheap mode as an option, and **measure** the radiance
difference on a few columns to decide what the ML set needs.

---

## 6. Core vSmartMOM enhancement: refractive index across the window (the "Sanghavi" extension)

Today the band-endpoint scheme (`model_from_parameters.jl:447-520`) evaluates Mie at
`λ_lo` and `λ_hi`, averages Greek coefs, and interpolates `k(λ)` — **with a fixed
refractive index**. For speciated/mixed aerosols the **effective RI itself changes across
the window** (and across RH via `AW`). The enhancement:

- Thread an **`n_eff(λ)`** (from the effective-RI mixing + the RI database) into the
  endpoint Mie calls so `_mie_fwd` uses `n_eff(λ_lo)` and `n_eff(λ_hi)`, mirroring what
  Sanghavi already does for AOD but now for the refractive index.
- Note there's already a separate **reference index** `n_ref` (used for `k_ref`
  normalization) — the new path generalizes the in-band index from "constant" to
  "endpoint-evaluated", consistent with the existing 3-D `τ_aer` k(λ) scaling.

This is a small, well-scoped core change; do it after Stage 1–3 so you have a validated
optics path to test it against.

---

## 7. Standalone package: why and how

**Why:** the GCHP I/O + sectional→optics machinery is *not* core RT. Keeping it in a
separate package keeps vSmartMOM focused on radiative computation, lets the student own and
iterate fast, and avoids entangling the diverged `gchp-io` branch with core releases.

**The one rule that keeps it clean:**
> **`GCHPAerosolOptics.jl` depends on `vSmartMOM`. vSmartMOM never depends back.**
The new package consumes `vSmartMOM.Scattering` (`make_mie_model`,
`compute_aerosol_optical_properties`, `compute_mie_ab!`, `Aerosol`, `GreekCoefs`,
`AerosolOptics`) and `NCDatasets`.

**What moves out** (all currently on `gchp-io`): `GCHPScene.jl`, `abstract_types.jl`,
`aod_diagnostics.jl`, `schemes/tomas15.jl`, the Mie-extinction LUT, the Benchmark writers,
and `data/aerosols/tomas_species.yaml`. **What stays in vSmartMOM:** Scattering / CoreRT /
Absorption and the `Aerosol`/`AerosolOptics` types.

**Mechanics for the student (with Claude's help):**
1. `git init` the new package (a skeleton is provided alongside this brief — see
   `package_skeleton/`).
2. In its `Project.toml`, add `vSmartMOM` and `NCDatasets` (+ `Distributions`,
   `Interpolations`, `YAML`).
3. Develop against your local vSmartMOM checkout: `pkg> dev /path/to/vSmartMOM.jl`.
4. Move the `gchp-io` files in, swap `..Scattering`/`..Aerosols` internal references for
   `using vSmartMOM`.
5. Port the existing tests (`test/test_aerosol_generic_tomas_aod.jl`,
   `test/test_layer_resolved_aerosol_optics.jl`) as the regression baseline.

> The **core enhancement in §6** is the exception — that one lands *in vSmartMOM*, not the
> new package.

---

## 8. The downstream goal: ML training package (design it up front)

Target: ~10⁶ scenes across many GCHP columns/dates → ML training set mapping
**input state → output reflected radiances**. Two design problems to settle before
generating data:

**(a) Input/output schema.**
- *Output:* TOA reflectances/radiances over the spectral window × viewing/solar geometry
  sweep (+ polarization if used). Store the wavelength and geometry grids once, per-scene
  radiances as the payload.
- *Input:* the atmospheric + aerosol state. This is where the dimensionality explodes.

**(b) Dimensionality reduction (the key issue you flagged).**
Raw per-layer fields (size-dist parameters, composition, loading over 50–72 layers) make a
hyper-dimensional input. Use an **EOF / SVD / PCA** basis per profile field:
- Assemble a training matrix of profiles across many scenes; SVD → keep the leading
  eigenvectors (shared basis) + per-scene scores (a handful of numbers).
- Store the **basis** with the dataset so inputs are reversible/interpretable.
- Candidate fields to compress: AOD/extinction profile, effective-radius profile,
  per-species fraction profiles, RH/`AW` profile.

> You have prior art for exactly this pattern in the AtmosTransport AI-training package
> (binary-only metric + EOF-style reduction). Reuse that thinking: define
> *this-is-input / this-is-output*, freeze a basis, and version it.

---

## 9. Staged roadmap (suggested order of work)

| Stage | Goal | Mostly exists? |
|---|---|---|
| 0 | Orient on `gchp-io`; run the AOD example on the real `.nc4`; reproduce a column AOD | ✅ run it |
| 1 | Sectional (A1) optics → **column AOD validated** against GCHP's own AOD diagnostic | ✅ exists, validate |
| 2 | Add **A2 spline** + **B lognormal-fit** behind one `AerosolRepresentation` interface | 🟡 fit is a stub |
| 3 | **Bridge:** sectional/parametric → per-layer `AerosolOptics` (τ, ϖ, Greek) | 🔴 new — core piece |
| 4 | Core enhancement: **RI(λ) across the window** at band endpoints (§6) | 🔴 new, in vSmartMOM |
| 5 | **End-to-end RT** from one GCHP column (the "model setup file") | 🔴 after 3–4 |
| 6 | Scale to many columns → **ML dataset** + EOF/SVD input reduction (§8) | 🔴 design now |

Decision checkpoints to bring back to the group: (i) Mode A vs B as the *primary*
representation; (ii) mixing rule for internal mixtures (volume vs Maxwell-Garnett vs
Bruggeman; core/shell for BC coating?); (iii) height-varying size shape vs loading-only;
(iv) ML input schema + basis.

---

## 10. Pointers (file:line, all paths in `vSmartMOM.jl`)

- Reader: `gchp-io:src/IO/NetCDF/GCHPScene.jl`
- Sectional types + policies: `gchp-io:src/Aerosols/abstract_types.jl`
- AOD + effective RI: `gchp-io:src/Aerosols/aod_diagnostics.jl` (line 2 = scope note)
- TOMAS scheme + species table: `gchp-io:src/Aerosols/schemes/tomas15.jl`,
  `gchp-io:data/aerosols/tomas_species.yaml`
- Mie-extinction LUT: `gchp-io:src/Aerosols/LUT/mie_extinction_lut.jl`
- Bulk writers: `gchp-io:src/IO/Benchmark/{aod_diagnostic,netcdf_writer}.jl`
- Core Mie + optics: `src/Scattering/{types.jl,make_mie_model.jl}`,
  `compute_aerosol_optical_properties`
- Band-endpoint Mie + k(λ) interp (the RI-fixed Sanghavi path): 
  `src/CoreRT/tools/model_from_parameters.jl:447-520`
- Aerosol parameter file shape: `src/CoreRT/DefaultParameters.yaml` (`aerosols:` block)
- Example to run: `gchp-io:examples/gchp_aod_diagnostic.jl`
