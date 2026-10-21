# OCO-like retrieval linearization handoff

- **Prepared:** 2026-07-15
- **Resume on:** `suniti_multi_sensor`
- **Last code commit:** `ee21b857` (`fix: align aerosol profile tangents with forward model`)

**Status:** the analytic RT kernel is suitable, but the requested 33-state,
three-band retrieval is **not yet correct end to end**. Start with the
multi-aerosol mixing failure below.

## 1. Agreed retrieval state

Use one shared state vector for the O₂ A band, weak CO₂ band, and strong CO₂
band. Keep the existing vSmartMOM convention of aerosol columns first, then
gas columns, then surface columns:

| Global columns | State variables |
|---|---|
| `1:7` | aerosol 1: `τ_ref, nᵣ, nᵢ, rₘ, σ_g, z₀, σ₀` |
| `8:14` | aerosol 2: same seven variables |
| `15:21` | aerosol 3: same seven variables |
| `22:28` | aerosol 4: same seven variables |
| `29` | CO₂ profile scale `s_CO2` |
| `30` | H₂O profile scale `s_H2O` |
| `31` | O₂ A-band Lambertian reflectance `ρ_O2A` |
| `32` | weak-CO₂-band Lambertian reflectance `ρ_WCO2` |
| `33` | strong-CO₂-band Lambertian reflectance `ρ_SCO2` |

Recommended state semantics:

- `τ_ref` is each aerosol mode's actual column AOD at `λ_ref`.
- `(rₘ, σ_g)` are the user-facing median radius and geometric standard
  deviation, not the internal logarithmic parameters stored by
  `Distributions.LogNormal`.
- `(z₀, σ₀)` describe a LogNormal density in geometric altitude (km),
  normalized over the modeled column.
- `s_CO2` and `s_H2O` multiply fixed baseline VMR profiles while holding
  pressure, temperature, and dry-air column density fixed. This is simpler and
  better defined than calling the H₂O state a derivative with respect to `q`.
- Each Lambertian parameter affects only its own band.

If different semantics are wanted—especially a true `q` retrieval—settle that
before changing the gas tangent builder. A `q` derivative also changes the
H₂O conversion, dry-air column, Rayleigh optical depth, and other-gas columns.

## 2. Branch and workspace safety

The feature branch is pushed and synchronized:

```text
suniti_multi_sensor
  4c0758e7  feat: restore height-aware observer radiances
  0cae2db1  feat: linearize height-resolved observer radiances
  ee21b857  fix: align aerosol profile tangents with forward model
```

Do not stage or overwrite the unrelated OCORaman work currently present in the
shared worktree:

```text
M  workflows/OCORaman/scripts/chunked_createRamanLUT_O2A.jl
?? workflows/OCORaman/scripts/compile_RamanLUT_psurf_sza_shards.jl
?? workflows/OCORaman/scripts/compile_RamanLUT_unified.jl
?? workflows/OCORaman/scripts/run_all_o2a_psurf500_wurst.sh
?? workflows/OCORaman/scripts/run_o2a_psurf1000_vza_block.sh
?? workflows/OCORaman/scripts/run_o2a_psurf1000_vza_blocks.sh
?? workflows/OCORaman/scripts/run_o2a_psurf_vza_block.sh
?? workflows/OCORaman/scripts/run_o2a_psurf_vza_blocks.sh
```

Stage files explicitly. Do not use `git add .`.

## 3. What is already working

### Height-aware analytic RT

- Scalar/vector observer-height conventions, height-anchored profile framing,
  reduced-layer preservation, and near-TOA/near-BOA behavior are implemented.
- Elastic multisensor radiances and their analytic Jacobians are implemented.
- Solar-aligned output includes a separate unscattered direct-beam radiance and
  analytic tangent.
- The interlayer solve, diffuse fields, surface derivative, absorption
  derivative, and direct-beam derivative have finite-difference coverage.

### Single-aerosol derivative chain

For one Mie aerosol, the code constructs seven tangent columns and propagates
them through δ-M scaling and the MOM kernel:

```text
τ_ref, nᵣ, nᵢ, size location, size width,
vertical-profile location, vertical-profile width
```

`ParameterLayout` already concatenates seven columns per aerosol, so its
arithmetic produces the desired four blocks `1:7`, `8:14`, `15:21`, and
`22:28`.

### Profile tangent parity fixed on this branch

Commit `ee21b857` made LinMode use the same normalized aerosol profile as the
forward builder. Pressure `Normal(p₀, σp)` and stored
`LogNormal(log(z₀), σ₀)` tangents now reproduce the *current* forward model.
The H=5 km independent `τ_ref` finite-difference comparison closed to:

| Output | Maximum relative error | Maximum absolute error |
|---|---:|---:|
| upwelling | `3.64e-9` | `2.81e-10` |
| diffuse downwelling | `5.18e-9` | `9.93e-12` |
| direct beam | `2.29e-7` | `2.19e-13` |
| total downwelling | `5.24e-9` | `1.00e-11` |

The direct derivative is only about `9.56e-7`, so its absolute closure is the
more useful diagnostic.

### Baseline regression state

Before this handoff, the branch passed:

- multisensor suite: `225/225`
- Jacobian unit suite: `29/29`
- forward/LinMode builder parity: `34/34`
- perturbation utility: `12/12`
- Aqua package health: `9/9`

## 4. Blocker 1 — multispectral mixing fails at aerosol 2

This is the first task tomorrow.

### Reproduction

Run from `test/`:

```julia
using vSmartMOM
using vSmartMOM.CoreRT

p = parameters_from_yaml("test_parameters/JacobianTestFast.yaml")
base = only(p.scattering_params.rt_aerosols)
aerosols = [deepcopy(base) for _ in 1:4]

for (i, aer) in enumerate(aerosols)
    aer.τ_ref = 0.01i
    aer.aerosol.nᵣ = 1.25 + 0.02i
end
p.scattering_params.rt_aerosols = aerosols

model, lin_model = model_from_parameters(LinMode(), p)
@show size(model.τ_aer[1])       # (4, 4, 5)
@show size(lin_model.τ̇_aer[1]) # (4, 7, 4, 5)

rt_run(model, lin_model, 4, size(lin_model.τ̇_abs[1], 1), 1)
```

The builder succeeds. During `m=0`, aerosol 2 fails in
`src/CoreRT/types_lin.jl` with:

```text
DimensionMismatch: new dimensions (6, 6, 1, 1) must be consistent
with array size 144
```

The stack enters the scattering-plus-scattering tangent mixer around
`types_lin.jl:263-349`, specifically the reshape near line 338, from
`compEffectiveLayerProperties_lin.jl:111`.

### Root cause

After Rayleigh is mixed with aerosol 1, the combined forward phase matrices
are spectrally expanded:

```text
Z       : (nμ, nμ, nSpec)
Z tangent: (nμ, nμ, nSpec, 7)
```

The `ẋ !== nothing` branch then treats the already-mixed `x.Z` as a 2-D,
spectrally invariant matrix and forces it to `(nμ, nμ, 1, 1)`. This happens to
work for `nSpec=1` and fails for `nSpec>1`.

### Intended fix

Normalize both operands before applying the quotient rule:

- forward `Z`: accept `(nμ,nμ)` or `(nμ,nμ,nSpec)` and present it as
  `(nμ,nμ,nSpec,1)`, using a singleton spectral axis only when broadcasting is
  actually required;
- tangent `Ż`: accept `(nμ,nμ,nParam)` or
  `(nμ,nμ,nSpec,nParam)` and present it as
  `(nμ,nμ,nSpec,nParam)`;
- do this symmetrically for the existing mixture `x` and newly added aerosol
  `y`, for both `Z⁺⁺` and `Z⁻⁺`.

Keep the helper array-backend generic. Do not force GPU arrays through host
`Array` conversion.

### Tests required for this fix

1. Operator-level test with two aerosols and `nSpec=1`.
2. Operator-level test with two aerosols and `nSpec=4`.
3. Full RT smoke with four aerosols and `nSpec=4`.
4. Verify the tangent dimension is `28` before gas/surface columns.
5. Perturb one AOD in each aerosol independently and verify only its seven-
   column block carries the expected seed.
6. Central finite differences for at least one parameter in aerosol 1 and
   aerosol 4, then expand to all 28 columns after the semantic fixes below.

## 5. Blocker 2 — aerosol state semantics are not yet retrieval-correct

Fixing the reshape makes the code run, but it does not by itself make the
four-mode retrieval correct.

### 5.1 Reference-AOD normalization differs between forward and LinMode

`ScatteringParameters.n_ref` is one shared complex value. When omitted, the
parser takes it from aerosol 1 (`src/IO/Parameters.jl:1046-1049`). The forward
builder uses that shared `n_ref` to compute `k_ref` for **every** aerosol
(`model_from_parameters.jl:419-435`). LinMode instead computes `k_ref` using
each aerosol's own microphysics (`lin_model_from_parameters.jl:280-288`).

For modes 2–4 with distinct refractive indices, forward and LinMode therefore
do not represent the same base optical state.

Recommended resolution: define `τ_ref` as each mode's AOD at `λ_ref`, compute
that mode's `k_ref` using its own size distribution and refractive index in
both builders, and retain the quotient-rule derivative of `k(λ)/k_ref`. This
keeps `τ(λ_ref) = τ_ref` when refractive index or size changes. If the legacy
shared-`n_ref` convention must remain, it needs an explicit opt-in and a
matching LinMode derivative definition.

### 5.2 Size-distribution columns are in stored log coordinates

YAML `μ, σ` are user-facing median radius and geometric standard deviation,
but the parser stores them as:

```julia
LogNormal(log(μ), log(σ))
```

(`src/IO/Parameters.jl:550-552`). The current Mie weight derivatives in
`src/Scattering/mie_helper_functions_lin.jl:235-249` differentiate with
respect to the stored `LogNormal.μ` and `.σ`. Therefore columns 4–5 currently
mean derivatives with respect to `log(rₘ)` and `log(σ_g)`, while public docs
label them `rₘ` and `σ_g`.

For the agreed physical state, apply the chain rule:

```math
∂/∂rₘ = (1/rₘ) ∂/∂log(rₘ),
∂/∂σ_g = (1/σ_g) ∂/∂log(σ_g).
```

Alternatively, deliberately retrieve log-parameters and rename the state
columns. Do not leave the current name/derivative mismatch implicit.

There is also a normalization issue to fix while this code is open.
`compute_wₓ` normalizes `wₓ` and then uses `sum(wₓ)`—now equal to one—in the
quotient rule for `ẇₓ`. Preserve the raw normalization `S = sum(w_raw)` and
use

```math
∂(w/S) = (∂w - (w/S) ∂S)/S.
```

Validate this helper directly against finite differences before relying on
the full Mie/RT comparison.

### 5.3 The altitude LogNormal still runs on pressure

The parser stores `(z₀, σ₀)` as `LogNormal(log(z₀), σ₀)`, but its own note at
`src/IO/Parameters.jl:537-542` records that the forward helper evaluates this
distribution on `profile.p_full`. Commit `ee21b857` correctly linearizes that
current behavior; it does not turn it into a physical altitude profile.

Recommended altitude implementation:

1. Obtain geometric half-level altitudes with `half_level_altitudes(profile)`.
2. Integrate the LogNormal between each pair of height interfaces, preferably
   using CDF differences rather than center sampling.
3. Normalize over the finite model column so the layer AODs sum to `τ_ref`.
4. Differentiate the normalized layer weights with respect to user-facing
   `(z₀, σ₀)`.
5. Confirm that an inserted observer height is an exact CDF integration
   boundary and remains present under layer reduction.

Pressure-form `Normal(p₀,σp)` should remain available as its own explicit
profile-coordinate convention.

### 5.4 Mie implementation limits

- LinMode Mie dispatch is implemented for `NAI2`, not `PCW`
  (`src/Scattering/compute_NAI2_lin.jl`). Use `NAI2()` for this work.
- Analytic phase-function aerosols are explicitly rejected by LinMode.
- Linearized Mie optical properties are built on CPU even when the RT arrays
  run on CUDA (`lin_model_from_parameters.jl:272-278`). This is a performance
  limitation, not an RT correctness blocker.

## 6. Blocker 3 — gas columns need named, profile-scale seeds

### Current behavior

`τ̇_abs[band]` is allocated as:

```text
(1 + number_of_unique_variable_species, nSpec, nLayer)
```

Row 1 is reserved for q-derived H₂O. A variable molecule at local position
`j` in a band is written to row `j+1`. Assignment is by local list position,
not by molecule name.

The narrow OCO-like convention

```text
variable_molecules = [[], ["CO2"], ["CO2"]]
```

therefore gives local row 1 = H₂O and row 2 = CO₂, with a zero CO₂ row in the
O₂ A band. This can be made to work, but reordered or disjoint species lists
can silently swap/collapse columns. There is no species metadata in
`ParameterLayout`.

The current seed

```math
τ̇ = σ(p,T,ν) VCD_dry
```

is a derivative with respect to an additive absolute VMR applied coherently at
each layer. It is not the desired derivative with respect to a multiplier of a
nonuniform baseline profile.

### Required behavior

Introduce an explicit global gas layout, at minimum:

```text
CO2 -> global column 29
H2O -> global column 30
```

For a profile multiplier at base scale 1, seed each layer with:

```math
∂τ_species(ν,z)/∂s_species = τ_species,base(ν,z).
```

That requires retaining or constructing species-separated absorption before
it is summed into total `τ_abs`.

For this first retrieval implementation, define the scale as changing gas
absorption while holding `T`, `p`, dry VCD, and Rayleigh fixed. Do not call it
a derivative with respect to `q`.

### H₂O contributions to decide explicitly

- line absorption: must be included;
- MT_CKD continuum: currently contributes only to the forward optical depth,
  with no Jacobian;
- CIA: currently fixed, with no Jacobian;
- pressure/temperature/line-shape feedback: outside the two-scalar gas state.

Record whether continuum and any H₂O-dependent CIA are included in `s_H2O`.
If omitted initially, make the omission explicit in the result metadata/docs.

### Gas tests

1. Assert name-to-column mapping is invariant to per-band molecule-list order.
2. O₂ A band: CO₂ column is exactly zero; H₂O column behaves as configured.
3. Weak/strong CO₂ bands: both share global CO₂ column 29.
4. Central FD of a nonuniform CO₂ profile multiplier.
5. Central FD of a nonuniform H₂O profile multiplier.
6. If continuum is enabled, FD the complete line+continuum result.

## 7. Blocker 4 — assemble one global three-band Jacobian

The public linearized entry point accepts one integer `i_band` per call. For
the current model each band would return:

```text
28 aerosol + 2 gas + 1 local Lambertian = 31 columns
```

The surface code deliberately seeds `surface_index(layout, 1)` for every
per-band call (`src/CoreRT/rt_run_lin.jl:332-338`). There is no result with
three independent surface columns or a shared 33-column layout.

### Recommended implementation

Keep one RT solve per band, but make every solve use the same global layout:

```text
28 aerosol + 2 named gas + 3 surface = 33 columns
```

Add an explicit active surface slot rather than inferring it from `i_band`:

```julia
rt_run(model, lin_model, 4, 2, 3;
       i_band=iband, surface_slot=iband)
```

Preserve `surface_slot=1` as the backward-compatible default. The O₂ A call
then seeds column 31, weak CO₂ column 32, and strong CO₂ column 33; the two
inactive surface columns remain exactly zero.

Add a thin multiband wrapper that:

1. validates a common 33-column layout;
2. runs the three bands;
3. returns band-labelled results (prefer this over hiding unequal wavelength
   grids in one array);
4. optionally concatenates radiances/Jacobians along the spectral dimension
   when VZA and Stokes dimensions match;
5. retains band spectral grids and column names.

Do not merely concatenate three 31-column Jacobians—the local surface column
and current local gas order must first be mapped into the global state.

## 8. Recommended implementation sequence

Use separate, reviewable commits:

1. **Fix multispectral multi-aerosol mixing.** Add operator tests and the
   four-aerosol RT smoke reproduction.
2. **Make aerosol semantics consistent.** Resolve per-mode `k_ref`, size
   parameter chain factors, and true altitude LogNormal profiles; verify
   forward/LinMode base parity.
3. **Introduce named gas profile-scale tangents.** Keep the initial state at
   fixed atmospheric structure and document continuum scope.
4. **Add global multiband/surface layout.** Implement `surface_slot` and a
   three-band result wrapper with the agreed 33 columns.
5. **Add the retrieval benchmark.** Validate every column against independent
   central finite differences.
6. **Run GPU parity after CPU correctness.** Linearized Mie remains CPU-built,
   but the 33-column RT propagation must still work on CUDA/Float32.

## 9. Acceptance criteria

Do not call the retrieval supported until all of these hold:

- [ ] Four distinct aerosols run with `nSpec > 1` without shape errors.
- [ ] Forward-only and LinMode base optical properties are equal for all four
      modes, including distinct refractive indices and size distributions.
- [ ] AOD at `λ_ref` equals each mode's `τ_ref` independently.
- [ ] Size columns have explicitly documented physical or log-space semantics
      and match those semantics by FD.
- [ ] Altitude `(z₀,σ₀)` profiles are evaluated in km, normalized, and split
      exactly at inserted sensor heights.
- [ ] CO₂ and H₂O are mapped by name to global columns 29 and 30.
- [ ] Gas columns are derivatives of the agreed profile multipliers.
- [ ] Each band has a 33-column Jacobian with only its own Lambertian column
      active among columns 31–33.
- [ ] All 28 aerosol, 2 gas, and 3 surface columns pass independent central FD
      checks at representative wavelengths/geometries.
- [ ] Near-TOA/near-BOA and solar-aligned direct-beam tests remain green.
- [ ] CPU Float64, CPU Float32, and CUDA RT agree within established tolerances.
- [ ] Existing benchmark and package-health suites remain green.

## 10. Fast validation commands

Run tests from `test/` because profiles and YAML paths are relative:

```bash
cd test

# Existing height/multisensor and analytic-Jacobian regression gates
env -u LD_LIBRARY_PATH julia --project=. --startup-file=no test_multisensor_heights.jl
env -u LD_LIBRARY_PATH julia --project=. --startup-file=no test_jacobians_unit.jl

# All CI-buildable forward/LinMode configurations
env -u LD_LIBRARY_PATH julia --project=. --startup-file=no test_lin_forward_loop_parity.jl

# FD parameter perturbation semantics and package health
env -u LD_LIBRARY_PATH julia --project=. --startup-file=no test_perturb_parameters.jl
env -u LD_LIBRARY_PATH julia --project=. --startup-file=no test_aqua.jl
```

For the new work, add a dedicated test such as
`test_oco_retrieval_linearization.jl` instead of overloading the height test
with expensive four-aerosol/three-band setup.

Use a cheap first fixture (`nquad_radius≈100`, 2–4 wavelengths per band,
`nstreams=3`, CPU Float64). Once every column closes by FD, add the realistic
OCO grids and CUDA/Float32 checks.

## 11. First actions tomorrow

1. Confirm the branch and dirty-file exclusions with `git status`.
2. Add the four-aerosol reproduction as a failing regression test.
3. Fix the spectral `Z/Ż` normalization in the scattering tangent mixer.
4. Re-run the new smoke test plus `test_jacobians_unit.jl`.
5. Then settle the three state-semantic decisions before expanding the API:
   - physical versus log size parameters;
   - true altitude LogNormal implementation;
   - exact H₂O scale/continuum definition.

The immediate runtime bug is localized. The larger task is making the
retrieval state names, forward model, tangent seeds, and multiband column map
describe the same 33 physical quantities.
