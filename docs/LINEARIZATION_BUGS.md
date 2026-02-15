# Linearized RT Bug Log

This document catalogs all bugs found and fixed during the integration of the
analytic Jacobian (linearized RT) code from the `sanghavi` branch into the
unified `TOMAS-aerosols` branch.

---

## Bug 1: Loop range in `compEffectiveLayerProperties_lin.jl`

**File:** `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties_lin.jl`  
**Line:** ~447  
**Symptom:** Only the last parameter's `τ̇_sum` was accumulated; all others were skipped.  
**Root cause:** `for ip=nParams` iterates only over the scalar value `nParams` (a single
iteration at `ip = nParams`), instead of iterating over the range `1:nParams`.  
**Fix:** Changed `for ip=nParams` → `for ip=1:nParams`.  
**Impact:** All per-parameter cumulative optical depth derivatives were incorrect
except the last one.

---

## Bug 2: Loop variable and broadcasting in `rt_kernel_lin.jl` (old version)

**File:** `src/CoreRT/CoreKernel/rt_kernel_lin.jl` (commented-out old code block)  
**Lines:** ~223, ~225  
**Symptom:** `temp_lin` multiplication used `*` instead of `.*` (element-wise),
and the loop iterated over `length(τ_λ)` instead of `length(τ)`.  
**Root cause:** Copy-paste from the forward code without adjusting for the
linearized variable names (`τ_λ` was renamed to `τ` in the new code path).  
**Fix:** Changed `*` → `.*` for broadcasting; changed `length(τ_λ)` → `length(τ)`.  
**Impact:** Incorrect non-scattering layer transmission derivatives; dimension
mismatch errors at runtime.

---

## Bug 3: Wrong denominator in Z⁻⁺ source derivative in `elemental_lin.jl`

**File:** `src/CoreRT/CoreKernel/elemental_lin.jl`  
**Line:** ~426  
**Symptom:** The Z-derivative of the downward source function `J̇₀⁻[3,i,1,n]`
was computed as `J₀⁻[i,1,n] / Z⁺⁺_I₀` instead of `J₀⁻[i,1,n] / Z⁻⁺_I₀`.  
**Root cause:** Typo — `Z⁺⁺_I₀` was used where `Z⁻⁺_I₀` is required.  
The downward source `J₀⁻` depends on Z⁻⁺ (not Z⁺⁺), so the derivative
`∂J₀⁻/∂Z` must divide by `Z⁻⁺_I₀`.  
**Fix:** Changed `J₀⁻[i, 1, n] / Z⁺⁺_I₀` → `J₀⁻[i, 1, n] / Z⁻⁺_I₀`.  
**Impact:** Incorrect Z-derivative of downward source in every elemental layer,
propagating errors through doubling and interaction to the final Jacobians.

---

## Bug 4: (Verified correct) Return values in `rt_run_lin.jl`

**File:** `src/CoreRT/rt_run_lin.jl`  
**Symptom:** None — verified that `(R, T, Ṙ, Ṫ)` return values and `iBand`
handling are correct.

---

## Bug 5: `F₀` (solar irradiance) initialization in `rt_run_lin.jl`

**File:** `src/CoreRT/rt_run_lin.jl`  
**Symptom:** `BoundsError: attempt to access 1×1 Matrix{Float64} at index [2, 1]`
in `elemental_lin.jl:390`.  
**Root cause:** `RS_type.F₀` was initialized as a 1×1 matrix but the elemental
code indexes it as `F₀[pol_index, spec_index]` where `pol_index` can be > 1
(for Stokes) and `spec_index` can be > 1 (for multiple spectral points).  
**Fix:** Added logic to check `if size(F₀) != (pol_type.n, nSpec)` and reinitialize
`F₀` as `ones(pol_type.n, nSpec)` when dimensions don't match.  
**Impact:** Runtime crash when running with multiple spectral points or
polarization components.

---

## Bug 6: Active `@show` debug statements

**Files:** Multiple `_lin` files:
- `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties_lin.jl`
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl`
- `src/CoreRT/rt_run_lin.jl`
- `src/CoreRT/tools/lin_model_from_parameters.jl`
- `src/CoreRT/tools/model_from_parameters.jl`
- `src/CoreRT/Surfaces/lambertian_surface_lin.jl`

**Symptom:** Excessive console output during runs, some `@show` statements
printing large arrays causing significant slowdown.  
**Fix:** Commented out all active `@show` statements.  
**Impact:** Performance degradation and noisy output; no numerical effect.

---

## Bug 7: Case mismatch — `J₀⁺`/`J₀⁻` vs `j₀⁺`/`j₀⁻` field names

**Files:**
- `src/CoreRT/CoreKernel/elemental_lin.jl`
- `src/CoreRT/CoreKernel/doubling_lin.jl`
- `src/CoreRT/CoreKernel/interaction_lin.jl`
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl`
- `src/CoreRT/Surfaces/lambertian_surface_lin.jl`

**Symptom:** `FieldError: type vSmartMOM.CoreRT.AddedLayer has no field J₀⁺`  
**Root cause:** The `AddedLayer` struct uses lowercase field names `j₀⁺` and `j₀⁻`,
but the linearized code (from the `sanghavi` branch) referenced uppercase `J₀⁺`
and `J₀⁻`. This is a convention difference between branches.  
**Fix:** Replaced all `added_layer.J₀⁺` → `added_layer.j₀⁺` and
`added_layer.J₀⁻` → `added_layer.j₀⁻` in the affected files.  
**Impact:** Runtime crash (FieldError) whenever the linearized code accessed the
added layer's source terms.

---

## Bug 8: `AddedLayer` constructor argument count

**File:** `src/CoreRT/tools/rt_helper_functions_lin.jl`  
**Symptom:** `MethodError: no method matching vSmartMOM.CoreRT.AddedLayer(...)`
(6 arguments provided, 10 expected).  
**Root cause:** The `sanghavi` branch's `AddedLayer` had fewer fields than the
unified branch's version (which includes `j₀⁺`, `j₀⁻`, `r⁺⁻`, `r⁻⁺`,
`t⁺⁺`, `t⁻⁻` and additional scattering-related fields).  
**Fix:** Modified `make_added_layer` to pass all 10 required arguments.  
**Impact:** Runtime crash during layer construction.

---

## Bug 9: Type parameter mismatch in linearized optical property structs

**File:** `src/CoreRT/types_lin.jl`  
**Symptom:** `MethodError: no method matching
vSmartMOM.CoreRT.UmbrellaCoreScatteringOpticalProperties(...)`  
**Root cause:** The struct had overly restrictive type parameters that didn't
match the types produced by the linearized model construction.  
**Fix:** Relaxed type parameters in `UmbrellaCoreScatteringOpticalProperties`
and `UmbrellaCoreAbsorptionOpticalProperties`.  
**Impact:** Runtime crash during model construction.

---

## Bug 10: `τ_rayl` dimension mismatch (forward model)

**File:** `src/CoreRT/tools/model_from_parameters.jl`  
**Line:** ~42  
**Symptom:** `DimensionMismatch: array could not be broadcast to match destination`
when computing Rayleigh optical depth.  
**Root cause:** `τ_rayl` was initialized with `length(params.T)` (the raw profile
size before reduction), but `getRayleighLayerOptProp` returns data sized to
`length(profile.p_full)` (the reduced profile size, after `reduce_profile`).  
**Fix:** Changed `length(params.T)` → `length(profile.p_full)`.  
**Impact:** Runtime crash when `profile_reduction_n` is set in the YAML config.

---

## Bug 11: `τ_rayl` dimension mismatch (linearized model) — 4 occurrences

**File:** `src/CoreRT/tools/lin_model_from_parameters.jl`  
**Lines:** ~164, ~635, ~991, ~1286 (one per RS_type code path: noRS, VRS, RRS, RVRS)  
**Symptom:** Same as Bug 10 — `DimensionMismatch` during Rayleigh optical depth.  
**Root cause:** Same as Bug 10 — used `length(profile.T)` instead of
`length(profile.p_full)`.  
**Fix:** Changed all 4 occurrences to use `length(profile.p_full)`.  
**Impact:** Runtime crash when `profile_reduction_n` is set.

---

## Bug 12: Division by zero in elemental linearized derivatives

**File:** `src/CoreRT/CoreKernel/elemental_lin.jl`  
**Lines:** ~284, ~321, ~341, ~401, ~403, ~412, ~414, ~424, ~426  
**Symptom:** `NaN` values appearing in the Jacobian output (`Ṙ`) for aerosol
and gas parameters.  
**Root cause:** The Z-derivatives and ω-derivatives in the elemental layer are
computed as ratios, e.g.:
```julia
J̇₀⁺[3, i, 1, n] = J₀⁺[i, 1, n] / Z⁺⁺_I₀
J̇₀⁺[2, i, 1, n] = J₀⁺[i, 1, n] / ϖ_λ[n]
```
When `Z⁺⁺_I₀ == 0` or `ϖ_λ[n] == 0` (which happens for non-scattering layers
or spectral points with zero single-scattering albedo), these produce `0/0 = NaN`
or `X/0 = Inf`.  
**Fix:** Added safe division guards:
```julia
J̇₀⁺[3, i, 1, n] = Z⁺⁺_I₀ == 0 ? FT(0) : J₀⁺[i, 1, n] / Z⁺⁺_I₀
J̇₀⁺[2, i, 1, n] = ϖ_λ[n] == 0 ? FT(0) : J₀⁺[i, 1, n] / ϖ_λ[n]
```
Applied to all 9 affected locations.  
**Impact:** `NaN` propagation through doubling and interaction, corrupting all
downstream Jacobian values.

---

## Bug 13: Depolarization constant mismatch between forward and linearized paths

**File:** `src/CoreRT/tools/lin_model_from_parameters.jl`  
**Line:** ~205 (noRS code path)  
**Symptom:** Forward reflectance `R` computed via the linearized model
(`model_from_parameters(LinMode(), params)`) differed from the standalone
forward model (`model_from_parameters(params)`) by up to 0.035 (relative).
Specifically, `τ_rayl` sums differed:
```
τ_rayl fwd sum: 4.38797649085687
τ_rayl lin sum: 4.486852374654676
```
**Root cause:** The forward `model_from_parameters` computes Rayleigh optical
depth using `params.depol` (the depolarization factor from the YAML config),
while the linearized `model_from_parameters(LinMode(), ...)` computed it using
`depol_air_Rayleigh` — a value derived from molecular constants via
`InelasticScattering.compute_γ_air_Rayleigh!()`. These give slightly different
values (~0.0279 from params vs ~0.0286 from molecular constants).  
**Fix:** Changed the noRS code path in `lin_model_from_parameters.jl` to use
`params.depol` for consistency with the forward model:
```julia
τ_rayl[i_band] .= getRayleighLayerOptProp(profile.p_half[end],
                      curr_band_λ, params.depol, profile.vcd_dry)
```
**Note:** The VRS/RRS/RVRS paths intentionally use the molecular-derived
`depol_air_Rayleigh` for their specific physics; this fix only applies to noRS.  
**Impact:** Forward R mismatch between forward-only and linearized paths,
making it impossible to validate Jacobians via finite-difference comparison.

---

## Bug 14: Chain rule uses full `τ̇` instead of elemental `dτ̇` (CRITICAL)

**File:** `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl`  
**Lines:** 51, 54, 58, 61, 65, 68 (every occurrence of `τ̇[iparam,:]`)  
**Symptom:** Jacobian values for all parameters are incorrect (too large by a
factor of 2^ndoubl) for non-TOA layers. Combined with already-large Mie
parameter derivatives (O(10⁸–10¹¹)), this can cause numerical overflow and `NaN`.  

**Root cause:** After the elemental + doubling steps, the core derivatives
(`ṫ⁺⁺[1,:,:,:]`, `ṙ⁻⁺[1,:,:,:]`, etc.) are derivatives with respect to the
**elemental** optical depth `dτ = τ / 2^ndoubl`. The chain rule to physical
parameters therefore requires:

```
∂R_layer/∂param = (∂R/∂dτ) × (∂dτ/∂param) + (∂R/∂ϖ) × (∂ϖ/∂param) + (∂R/∂Z) × (∂Z/∂param)
                = ṙ[1] × dτ̇ + ṙ[2] × ϖ̇ + ṙ[3] × Ż
```

where `dτ̇ = τ̇ / 2^ndoubl`.

The **TOA code** in `rt_kernel_lin.jl` (line 230) correctly uses `dτ̇`:
```julia
dτ̇ = τ̇ ./ 2^ndoubl   # line 155
...
composite_layer_lin.Ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ṫ⁺⁺[1,:,:,:] .* reshape(dτ̇[iparam,:], ...)
```

But `lin_added_layer_all_params!` (used for **all non-TOA layers**) incorrectly
used the full layer `τ̇`:
```julia
ap_ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ṫ⁺⁺[1,:,:,:] .* reshape(τ̇[iparam,:], ...)  # WRONG
```

This means the τ-derivative contribution for non-TOA layers was multiplied by
2^ndoubl too much. Since `ndoubl` is typically 3–8 (giving factors of 8–256),
and the Mie parameter derivatives are already O(10⁸–10¹¹), the products can
exceed Float64 range, producing `Inf` → `NaN` in subsequent arithmetic.

**Fix:** Added `ndoubl` parameter to `lin_added_layer_all_params!` and compute
`dτ̇ = τ̇ / 2^ndoubl` inside:
```julia
function lin_added_layer_all_params_helper!(..., ndoubl::Int) where ...
    ...
    dτ̇ = τ̇ ./ FT(2^ndoubl)
    ...
    ap_ṫ⁺⁺[iparam,:,:,:] .= ṫ⁺⁺[1,:,:,:] .* reshape(dτ̇[iparam,:], ...) .+ ...
```
Updated the call site in `rt_kernel_lin.jl` to pass `ndoubl`.

**Impact:** All Jacobians for non-TOA layers were incorrect by a factor of
2^ndoubl in the optical-depth contribution. This is the most likely root cause
of the persistent `NaN` in `Ṙ` for aerosol/gas parameters, since these have
large `τ̇` values that, when not divided by 2^ndoubl, overflow during matrix
operations in the interaction step.

---

## Bug 15: `parse_surface_str` regex didn't handle type parameters

**File:** `src/IO/Parameters.jl`  
**Symptom:** `AssertionError: Invalid surface specification:
'LambertianSurfaceScalar{Float64}(0.0)'`  
**Root cause:** The regex for parsing surface strings from YAML didn't account
for Julia's type parameter syntax `{Float64}`.  
**Fix:** Updated regex to optionally match `{...}` after the type name.  
**Impact:** Runtime crash when loading YAML files that specify surfaces with
explicit type parameters.

---

## Bug 16: Division by zero in Z-derivative of r and t matrices (CRITICAL)

**File:** `src/CoreRT/CoreKernel/elemental_lin.jl`  
**Lines:** ~286, ~323, ~343  
**Symptom:** NaN in `ṙ⁻⁺[3,:,:,:]` and `ṫ⁺⁺[3,:,:,:]` (core derivative 3 = Z-derivative)
at every layer (iz=1 through iz=5), propagating through doubling and interaction
to corrupt all final Jacobians.  

Step-by-step debug tracing revealed:
```
[DEBUG iz=1] NaN AFTER ELEMENTAL: ṙ⁻⁺=true, ṫ⁺⁺=true, J̇₀⁺=false, J̇₀⁻=false
    core[1]: ṙ⁻⁺=false, ṫ⁺⁺=false   ← τ-derivative OK
    core[2]: ṙ⁻⁺=false, ṫ⁺⁺=false   ← ϖ-derivative OK
    core[3]: ṙ⁻⁺=true, ṫ⁺⁺=true    ← Z-derivative NaN!
```

**Root cause:** The Z-derivative of reflection and transmission was computed as
a ratio `r/Z` or `t/Z`:
```julia
ṙ⁻⁺[3,i,j,n] = r⁻⁺[i,j,n] / Z⁻⁺[i,j,n2]   # line 286
ṫ⁺⁺[3,i,j,n] = t⁺⁺[i,j,n] / Z⁺⁺[i,j,n2]   # lines 323, 343
```
When `Z[i,j] = 0` (which occurs for many angular entries of the scattering
matrix), this produces `0/0 = NaN` (since r and t also contain Z as a factor).

Note: Bug 12 fixed the same issue for the source term derivatives (`J̇₀⁺`, `J̇₀⁻`)
using `Z == 0 ? 0 : val/Z` guards. But the fix was NOT applied to `ṙ` and `ṫ`,
which have the same problem.

**Fix:** Replaced all ratio formulas with mathematically equivalent direct formulas
that avoid division entirely. Since `r = ϖ·Z·f(τ,μ)`, the derivative
`∂r/∂Z = ϖ·f(τ,μ)` is computed directly:

```julia
# Line 286: was r⁻⁺[i,j,n] / Z⁻⁺[i,j,n2]
ṙ⁻⁺[3,i,j,n] = ϖ_λ[n] * (qp_μN[j]/(qp_μN[i]+qp_μN[j])) * wct[j] *
                 (1 - exp(-dτ_λ[n] * (1/qp_μN[i] + 1/qp_μN[j])))

# Line 323 (i≠j, same μ): was t⁺⁺[i,j,n] / Z⁺⁺[i,j,n2]
ṫ⁺⁺[3,i,j,n] = exp(-dτ_λ[n]/qp_μN[j]) * ϖ_λ[n] * (dτ_λ[n]/qp_μN[i]) * wct[j]

# Line 343 (i≠j, different μ): was t⁺⁺[i,j,n] / Z⁺⁺[i,j,n2]
ṫ⁺⁺[3,i,j,n] = ϖ_λ[n] * (qp_μN[j]/(qp_μN[i]-qp_μN[j])) * wct[j] *
                 (exp(-dτ_λ[n]/qp_μN[i]) - exp(-dτ_λ[n]/qp_μN[j]))
```

These direct formulas:
1. Never produce NaN (no division by potentially-zero values)
2. Are correct even when Z=0 (give the true derivative ϖ·f(τ,μ))
3. Are slightly more efficient (avoid division)

**Impact:** This was the root cause of all NaN values in the final Jacobians.
After this fix, all 9 Jacobian parameters across all layers and Fourier orders
are NaN-free and finite.

---

## Bug 17: Wrong variable references in ScatteringInterface_11 interaction (CRITICAL)

**File:** `src/CoreRT/CoreKernel/interaction_lin.jl`  
**Lines:** 189-191 (J̇₀⁻), 225 (J̇₀⁺)  
**Symptom:** Surface albedo Jacobian was exactly 2x the correct value (finite-difference).  

**Root cause:** In the `ScatteringInterface_11` case (both composite and added layers
have scattering), the source function derivative used the wrong variable:

Line 189: Used `ap_J̇₀⁻` (added layer derivative) instead of `J̇₀⁻` (composite layer
derivative) as the leading term. For the surface albedo parameter, `J̇₀⁻ = 0`
(atmosphere doesn't depend on surface), but `ap_J̇₀⁻ ≠ 0` (surface source depends
on albedo), causing the derivative to be counted twice.

Line 191: Used `ap_J̇₀⁺` (added layer) instead of `J̇₀⁺` (composite) for the
upwelling source derivative inside the reflection-correction term.

Line 225: Same bug pattern for the upwelling equation (J̇₀⁺).

**Fix:** Replaced three instances:
```julia
# Line 189: was ap_J̇₀⁻, now J̇₀⁻
tmpap_J̇₀⁻[iparam,...] .= J̇₀⁻[iparam,...] .+ ...

# Line 191: was ap_J̇₀⁺, now J̇₀⁺  
T01_inv ⊠ (... .+ r⁻⁺ ⊠ J̇₀⁺[iparam,...] .+ ...)

# Line 225: was ap_J̇₀⁺, now J̇₀⁺
T21_inv ⊠ (J̇₀⁺[iparam,...] .+ ...)
```

**Impact:** This affected ALL Jacobian parameters through the surface interaction
step. The surface albedo Jacobian was 2x too large; other parameters were also
affected when the surface interaction uses ScatteringInterface_11.

**Verification:** After fix, surface albedo Jacobian matches finite-difference to
5.9e-6 relative error across all spectral points.

---

## Bug 18: Index Mapping Error in `createAero` δ-M Derivative Assembly

**File:** `src/CoreRT/LayerOpticalProperties/compEffectiveLayerProperties_lin.jl`  
**Severity:** CRITICAL — all aerosol Jacobians were wrong  

**Symptom:** Aerosol τ_ref Jacobian was ~5.5x larger than finite-difference, and all
other aerosol sub-parameter Jacobians were cross-contaminated with neighboring parameters.

**Root cause:** The `createAero` function computes δ-M scaled aerosol optical property
derivatives for 7 sub-parameters `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]`. Three bugs:

1. **`τ̇_mod[1,:]`** was set to `(1-fᵗω̃)` — the bare scaling coefficient — instead
   of `(1-fᵗω̃) * τ̇Aer[1,:]` (the actual τ_ref derivative). This made parameter 1
   represent "d/d(τ_aer)" instead of "d/d(τ_ref)".

2. **`τ̇_mod[2:5,:]`** used `τ̇Aer[1:4,:]` (mapping τ_ref, nᵣ, nᵢ, rₘ) but combined
   with Mie derivatives `ω̃̇[2:5,:], ḟᵗ[2:5]` (mapping nᵣ, nᵢ, rₘ, σᵣ). This **off-by-one**
   index shift mixed τ_ref's τ_aer chain with nᵣ's Mie chain, nᵣ's τ_aer chain
   with nᵢ's Mie chain, etc.

3. **`τ̇_mod[6:7,:]`** used `τ̇Aer[5:6,:]` (σᵣ, p₀) instead of `τ̇Aer[6:7,:]` (p₀, σp).
   This meant the σp derivative was completely **lost** and the p₀ derivative was
   contaminated with σᵣ.

**Fix:**
```julia
# Before (WRONG):
τ̇_mod[1,:]   = (1 .- fᵗ * ω̃)                              # Missing τ̇Aer[1,:]
τ̇_mod[2:5,:] = (1-fᵗω̃)' .* τ̇Aer[1:4,:] .- tmp .* τAer'  # Off by one
τ̇_mod[6:7,:] = (1-fᵗω̃)' .* τ̇Aer[5:6,:]                  # Off by one

# After (CORRECT):
τ̇_mod[1,:]   .= (1 .- fᵗ * ω̃) .* τ̇Aer[1,:]               # τ_ref chain
τ̇_mod[2:5,:] .= (1-fᵗω̃)' .* τ̇Aer[2:5,:] .- tmp .* τAer'  # Mie params
τ̇_mod[6:7,:] .= (1-fᵗω̃)' .* τ̇Aer[6:7,:]                  # Profile params
```

**Impact:** ALL aerosol Jacobians (τ_ref, Mie parameters, profile parameters) were
incorrect. Each parameter's derivative was contaminated by or swapped with its neighbor's.

---

### Bug 19: Z Chain Rule Applied After Doubling (Element-wise is Wrong)

**Status:** CONFIRMED — causes ~12% max rel error for τ_ref, worse for p₀/σp  
**File:** `src/CoreRT/CoreKernel/rt_kernel_lin.jl`, `lin_added_layer_all_params.jl`  
**Reference:** Sanghavi, Davis & Eldering (2014); Sanghavi & Stephens (2015) Eq. (8)

**Root Cause:**

The RT kernel applies the chain rule from 3 core parameters (τ, ϖ, Z) to N physical
parameters **AFTER** the doubling step:

```
elemental! → doubling! → lin_added_layer_all_params! (chain rule) → interaction!
```

The chain rule for the Z term uses element-wise multiplication:
```julia
ap_ṫ⁺⁺[iparam,:,:,:] .= ṫ⁺⁺[1,:,:,:] .* dτ̇[iparam,:] .+
                          ṫ⁺⁺[2,:,:,:] .* ϖ̇[iparam,:] .+
                          ṫ⁺⁺[3,:,:,:] .* Ż⁺⁺[iparam,:,:,:]   ← ELEMENT-WISE
```

At the **elemental** level, `ṫ⁺⁺[3,i,j]` = ∂t_ij/∂Z_ij is correct because each t_ij
depends linearly on Z_ij only (the 4th-rank derivative tensor ∂t_ij/∂Z_kl is diagonal).

After **doubling**, the matrix products `T_doubled = T · G · T` mix all (i,j) indices:
```
∂T_doubled_ij/∂Z_kl ≠ 0  for (k,l) ≠ (i,j)
```

So `ṫ_doubled[3,i,j]` is NOT `∂T_ij/∂Z_ij` anymore — it contains mixed contributions
from all Z elements. The element-wise multiplication with Ż[iparam,i,j] then gives the
wrong result, ignoring all off-diagonal terms.

**Evidence:**
- Surface albedo (Ż = 0): error ≈ 1e-6 (essentially perfect) ✓
- τ_ref (Ż ≠ 0 from Rayleigh-aerosol Z mixing): error ≈ 12% max, 5.7% mean ✗
- p₀ (Ż ≠ 0 from aerosol profile redistribution): error ≈ 177% max ✗
- Gas VMR (Ż = 0 after gas addition): should be correct ✓

The pattern is clear: errors appear ONLY for parameters with nonzero Ż.

**The τ and ϖ chain rule terms are CORRECT** after doubling because they are scalars.
The Z term is the only matrix quantity, and its element-wise application is wrong after
the index-mixing matrix products in the doubling step.

**Fix:** Apply the chain rule BEFORE doubling:
```
elemental! → lin_added_layer_all_params! (chain rule) → doubling (N params) → interaction!
```

At the elemental level, the Z chain rule IS element-wise and correct. After applying
the chain rule, we have N physical parameter derivative matrices, which can then be
propagated through doubling using standard matrix products.

**Cost:** Doubling iterates over N parameters instead of 3. For typical N ≈ 10-15,
this is ~4× slower for the doubling step, but correctness is essential.

**Alternative (approximate):** For parameters where Ż is small relative to τ̇ and ϖ̇
(e.g., weak aerosol, Rayleigh-dominated), the error is small. But for optically thick
aerosol layers or near-backscatter geometries, the error can be large.

---

## Bug 20: Wrong product/quotient rule in Mie derivative interpolation

**File:** `src/CoreRT/tools/lin_model_from_parameters.jl`, lines 501-511

**Symptom:** Microphysical parameter Jacobians (nᵣ, nᵢ, rₘ, σₚ) are wrong.
Surface albedo, τ_ref, and profile parameters (p₀, σp) are unaffected because they
bypass the Mie derivative interpolation entirely.

**Root cause — two interrelated errors:**

### Bug 20a: k̇sca uses product of derivatives instead of product rule

The forward k_sca is correctly computed as `k_sca = k_ext * ω̃`. But the linearized
version used `k̇_sca = k̇_ext * ω̃̇` — the product of two derivatives — instead of the
product rule:

```
d(k_ext * ω̃)/dp = dk_ext/dp * ω̃ + k_ext * dω̃/dp = k̇_ext * ω̃ + k_ext * ω̃̇
```

### Bug 20b: ω̃̇ recovery divides by k̇ instead of k

The forward ω̃ is correctly recovered as `ω̃ = k_sca / k_ext`. But the linearized
version used `ω̃̇ = k̇sca_interp / k̇ext_interp` — dividing by the DERIVATIVE k̇
instead of the VALUE k. The correct quotient rule is:

```
dω̃/dp = d(k_sca/k_ext)/dp = (dk_sca/dp - ω̃ * dk_ext/dp) / k_ext
                            = (k̇_sca - ω̃ * k̇_ext) / k_ext
```

### Self-cancellation at grid points

Remarkably, at the two Mie grid points (no interpolation), the bugs cancel:
`(k̇ * ω̃̇) / k̇ = ω̃̇`, recovering the correct Mie output. With interpolation between
grid points, the cancellation is imperfect — the code computes a k̇-weighted average
of ω̃̇ instead of the proper quotient rule. Additionally, when `k̇ ≈ 0` (e.g., near
an extinction minimum w.r.t. a parameter), the division produces Inf/NaN.

### Diagnostic evidence

This bug only affects parameters that go through the Mie derivative chain (nᵣ, nᵢ,
rₘ, σₚ). The following parameters bypass this chain and thus are unaffected:
- Surface albedo: only affects boundary → perfect Jacobian
- τ_ref: only scales aerosol amount, Mie properties unchanged → ~5-12% error (from Z chain rule)
- p₀, σp: only redistribute aerosol vertically → peaks match ~10%

**Fix:** Applied correct product rule for k̇sca and quotient rule for ω̃̇.

---

## Bug 21: Missing k_band factor in τ̇_aer for microphysical parameters

**File:** `src/CoreRT/tools/lin_model_from_parameters.jl`, lines 562-568

**Symptom:** τ̇_aer for microphysical parameters (nᵣ, nᵢ, rₘ, σₚ) is spectrally wrong.

The aerosol optical depth per layer is:
```
τ_aer(λ, z) = (τ_ref / k_ref) * k(λ) * τₚ(z)
```

The derivative w.r.t. a microphysical parameter p (which affects both k(λ) and k_ref):
```
dτ_aer/dp = τ_ref * d(k(λ)/k_ref)/dp * τₚ
          = (τ_ref/k_ref * dk(λ)/dp  −  τ_ref * k(λ)/k_ref² * dk_ref/dp) * τₚ
```

The code had:
```julia
(τ_ref/k_ref * k̇_band  −  τ_ref/k_ref² * k̇_ref) * τₚ  # WRONG: missing k(λ)
```

The second term was missing the `k(λ)` factor (`aerosol_optics[i_band][i_aer].k`).
When k(λ) ≈ k_ref (narrow band near reference wavelength), the error is small
(factor ≈ 1). For wider bands or wavelengths far from the reference, the error
scales as `|1 − k(λ)/k_ref|`.

**Fix:** Added `.* aerosol_optics[i_band][i_aer].k` to the k_ref normalization term.

---

## Bug 22: Beam Attenuation Derivative Used Wrong Parameter Index (CRITICAL)

**Status:** FIXED  
**File:** `src/CoreRT/CoreKernel/elemental_lin.jl` (removal), `src/CoreRT/CoreKernel/rt_kernel_lin.jl` (addition)

**Symptom:** p₀ convergence test showed relative error *increasing* with decreasing FD step size,
stabilizing at ~425%. This is a hallmark of a wrong analytic derivative (not just truncation error).

**Root cause:** In `elemental_lin.jl` (lines 489-494), the beam attenuation derivative
`exp(-τ_sum/μ₀) * (-τ̇_sum/μ₀)` was baked into the 3-core source derivatives J̇₀⁺[1], J̇₀⁻[1],
J̇₀⁺[2], J̇₀⁻[2], J̇₀⁺[3], J̇₀⁻[3] using only `τ̇_sum[1,n]` — the FIRST parameter's cumulative
optical depth derivative. This is wrong because:

1. The beam attenuation depends on the FULL cumulative optical depth from TOA, not just one parameter's contribution
2. Each physical parameter has its OWN τ̇_sum, which must be used per-parameter
3. Baking this into the 3-core framework then applying the chain rule gives:
   `τ̇_sum[1] * dτ̇[iparam]` instead of the correct `τ̇_sum[iparam]`

**Fix (two parts):**

1. **Removed** the incorrect `τ̇_sum[1,n]` terms from all 6 core derivatives in `elemental_lin.jl`
2. **Added** per-parameter τ̇_sum beam attenuation in `rt_kernel_lin.jl` AFTER the chain rule:
```julia
for iparam = 1:nparams_τ_sum
    ap_J̇₀⁺[iparam,:,1,:] .+= j₀⁺[:,1,:] .* reshape(-τ̇_sum[iparam,:] ./ μ₀, 1, nspec)
    ap_J̇₀⁻[iparam,:,1,:] .+= j₀⁻[:,1,:] .* reshape(-τ̇_sum[iparam,:] ./ μ₀, 1, nspec)
end
```

**Impact:** Reduced p₀ error from ~425% (stabilized) to ~120% at small FD steps, and from
8.5% to 3.0% mean relative error in the unit test.

---

## Bug 23: Catastrophic Cancellation for Profile Redistribution Parameters (FUNDAMENTAL LIMITATION)

**Status:** DIAGNOSED — not a code bug but an inherent numerical precision limitation  
**Parameters affected:** p₀, σp (and any parameter that redistributes optical depth while conserving total)

**Symptom:** After fixing all code bugs (1-22), profile parameters (p₀, σp) still show
large Jacobian errors:
- With surface albedo > 0: ~3000% error, WRONG SIGN
- With surface albedo = 0: 28-76% error for most wavelengths, but 1.7% for wavelength 4

**Diagnosis (comprehensive):**

The following were all verified to be correct:
1. ✅ Layer-level τ̇, ϖ̇, Ż for p₀ — match FD to <1%
2. ✅ Cumulative τ̇_sum — matches FD
3. ✅ Combined layer optical properties (τ̇, ϖ̇, Ż after Rayleigh+aerosol+gas) — match FD to <1%
4. ✅ Chain rule algebra — verified correct
5. ✅ Doubling algebra — verified correct
6. ✅ Interaction (Adding) algebra — verified correct
7. ✅ Surface interaction algebra — verified correct
8. ✅ τ_ref Jacobian at albedo=0 — **perfect** (0.0% error)
9. ✅ τ_ref Jacobian at albedo=0.05 — **1.6% error** (from quadrature discretization)

**Root cause — catastrophic cancellation in the Adding method:**

Profile parameters like p₀ redistribute aerosol optical depth across layers while conserving
the total column optical depth (Σ τ̇[iz] ≈ 0). This means:

- Some layers have **negative** τ̇ (less aerosol when p₀ shifts)
- Other layers have **positive** τ̇ (more aerosol)
- The corresponding J̇₀⁻ contributions from these layers are large and nearly opposite

As layers are combined through the Adding method:
```
After layer 3: J̇₀⁻[6, quad4] ≈ -7.1e-7  (large negative)
After layer 4: J̇₀⁻[6, quad4] ≈ -3.6e-7  (partial cancellation)
After layer 5: J̇₀⁻[6, quad4] ≈ -5.5e-10 (near-complete cancellation — 3 orders of magnitude!)
```

This ~1000× cancellation means any absolute error at the O(1e-7) level produces
an O(1e-10) residual that dominates the tiny true derivative O(1e-10).

**Evidence — constant absolute error across wavelengths:**

| Wavelength | True dR/dp₀ | Analytic dR/dp₀ | Absolute Error | Relative Error |
|:---:|:---:|:---:|:---:|:---:|
| s=1 | -5.0e-11 | -8.8e-11 | 3.8e-11 | 76% |
| s=2 | -1.3e-10 | -1.7e-10 | 3.8e-11 | 29% |
| s=3 | -5.2e-11 | -9.1e-11 | 3.8e-11 | 73% |
| s=4 | -2.0e-9 | -2.0e-9 | 3.5e-11 | **1.7%** |

The absolute error is remarkably consistent (~3.8e-11), confirming it arises from the
cancellation of O(1e-7) quantities, not from a systematic algorithmic error.

**Surface amplification:**

The surface interaction amplifies this error because:
1. J̇₀⁺[6] (downwelling derivative) is numerically imprecise — O(1e-6) magnitude but
   contaminated by O(1e-11) absolute error from incomplete cancellation
2. Surface reflection converts this O(1e-6) J̇₀⁺ into an O(1e-9) contribution to J̇₀⁻
3. This spurious O(1e-9) term overwhelms the true dR of O(1e-11)
4. Result: wrong sign and ~3000% error when surface albedo > 0

**Why τ_ref doesn't have this problem:**

τ_ref scales ALL layers equally (positive τ̇ everywhere). There is no cancellation
between layers, so the derivative accumulates monotonically through the Adding method.
The 1.6% error for τ_ref with surface is from quadrature discretization (only 3 Gaussian
points per hemisphere), not from cancellation.

**Recommended solutions (in order of effort):**

1. **Accept limitation for profile params** — use central FD for p₀, σp; use analytic for
   τ_ref, surface albedo, gas VMR (which all work well)
2. **Increase quadrature** — more quadrature points reduce the discretization error for
   τ_ref (~1.6%) but do NOT fix the cancellation issue for profile params
3. **Adjoint (reverse-mode) method** — the long-term solution. The adjoint computes
   ∂R/∂p₀ as a scalar sum of (adjoint × layer_derivative) products, where the cancellation
   happens in a simple scalar sum that can use compensated arithmetic (e.g., Kahan summation).
   This avoids the matrix-level cancellation in the Adding method.
4. **Extended precision** — use BigFloat or compensated arithmetic for the Adding steps.
   Feasible but slow (~100× slowdown for BigFloat).

---

## Bug 24: Incorrect τ derivative in diagonal elemental T⁺⁺ branch (CRITICAL)

**Status:** FIXED  
**File:** `src/CoreRT/CoreKernel/elemental_lin.jl` (line ~351)

**Symptom:** Core-only directional finite-difference checks showed order-1 relative
error in the τ derivative for `T⁺⁺`, while `R⁻⁺`, `J₀⁺`, and `J₀⁻` were near machine precision.
This propagated to large τ mismatches after chain-rule + doubling.

**Root cause:** In the `i == j` branch of `get_elem_rt!`, the derivative formula
multiplied the entire bracket by `wct[i]`:
```julia
... * (-1 + ϖ*Z*(1 - dτ/μ)) * wct[i]
```
This incorrectly scaled the pure extinction term `-1` by `wct[i]`.
From the forward formula
`t = exp(-dτ/μ) * (1 + ϖ*Z*(dτ/μ)*wct)`,
the correct derivative is:
```julia
exp(-dτ/μ) * (1/μ) * (-1 + ϖ*Z*wct*(1 - dτ/μ))
```

**Fix:** Updated the expression to:
```julia
(-1 + ϖ_λ[n] * Z⁺⁺[i,i,n2] * wct[i] * (1 - dτ_λ[n] / qp_μN[i]))
```

**Verification:**
- `test/debug_tau_stagewise.jl` (before):
  - elemental `T++` τ derivative: max rel ~`1.0`
  - after chain+doubling `T++`: max rel ~`1.06`
- `test/debug_tau_stagewise.jl` (after):
  - elemental `T++` τ derivative: max rel `5.13e-9`
  - after chain+doubling `T++`: max rel `1.47e-6`
- `test/debug_core_kernel_fd.jl` (after):
  - τ case: `T++` max rel `1.47e-6`, `R-+` `1.09e-8`, `J0+` `1.11e-7`, `J0-` `1.01e-8`

**Impact:** Removed the dominant core τ analytic-derivative error.

---

## Bug 25: GPU scalar indexing from 4D CuArray views in linearized doubling (HIGH)

**Status:** FIXED  
**Files:**
- `ext/gpu_batched_cuda.jl` (new `batched_mul` overloads for CuArray views)
- Trigger site: `src/CoreRT/CoreKernel/doubling_lin.jl` (view slices like `ap_ṙ⁻⁺[iparam,:,:,:]`)

**Symptom:** GPU Jacobian run failed with:
- "Scalar indexing is disallowed"
- "scalar indexing of a GPU array"

**Root cause:** `doubling_allparams_helper!` passes 3D `SubArray` views from 4D CuArrays
to `batched_mul`. There was no CuArray-view overload, so fallback behavior triggered
host-side scalar access.

**Fix:** Added CuArray-view overloads:
- `batched_mul(SubArray{...,<:CuArray}, CuArray)`
- `batched_mul(CuArray, SubArray{...,<:CuArray})`
- `batched_mul(SubArray{...,<:CuArray}, SubArray{...,<:CuArray})`

Views are materialized to contiguous 3D CuArrays (`copy(view)`) before CUBLAS
strided batched GEMM.

**Verification:** `julia --project=. test/test_jacobians_GPU.jl` now passes (`4/4`)
with CPU vs GPU Jacobians matching within configured tolerance.

---

## Bug 23 Status Update (2026-02-15)

The earlier "WRONG SIGN / fundamentally broken" interpretation for profile derivatives
was partially confounded by Bug 24. After the τ-derivative fix:

- Core derivatives are now consistent in stagewise and directional checks.
- `p₀` comparisons are strongly step-size/method dependent (finite-difference numerics).
- Coarse one-sided FD and relative-only metrics near tiny true derivatives can still
  report large relative error.

Representative central-FD results from `test/test_p0_convergence.jl`:
- `δ/p₀ = 1e-3`: mean rel `6.94e-4`, max rel `1.35e-3`
- Raw wavelength check (first VZA/Stokes): rel errors `1.42e-4`, `1.07e-4`, `3.20e-4`, `6.78e-6`

So the current evidence supports: remaining large `%` values are largely FD validation
sensitivity/cancellation artifacts, not an obvious analytic core-kernel bug.

---

## Summary of Bug Categories (Updated 2026-02-15)

| Category | Bugs | Severity |
|----------|------|----------|
| Loop/range errors | #1, #2 | High — silent wrong results |
| Math errors | #3, #14, #18, #19, #20, #21, #22, #24 | Critical — wrong Jacobians |
| Index mapping | #17, #18, #22 | Critical — wrong variable/param references |
| Division by zero | #12, #16 | Critical — NaN propagation |
| Naming/casing | #7, #8 | High — runtime crash |
| Dimension mismatch | #5, #10, #11 | High — runtime crash |
| Type system | #9, #15 | Medium — runtime crash |
| Consistency | #13 | Medium — forward R mismatch |
| Structural (chain rule order) | #19 | Critical — all aerosol Jacobians affected |
| Beam attenuation derivative | #22 | Critical — profile params affected |
| Mie derivative interpolation | #20, #21 | Critical — microphysical Jacobians only |
| GPU derivative infrastructure | #25 | High — GPU Jacobian run failure |
| Numerical precision / FD sensitivity | #23 (revised) | Validation caveat |
| Debug artifacts | #6 | Low — performance only |
| Verified OK | #4 | N/A |

## Current Jacobian Status (Updated 2026-02-15)

| Check | Result | Evidence |
|------|--------|----------|
| Core τ/ϖ/Z/tau_sum directional derivatives | Good (`~1e-6` to `~1e-8` rel) | `test/debug_core_kernel_fd.jl` |
| Stagewise τ isolation (elemental → chain → doubling) | Good after Bug 24 fix | `test/debug_tau_stagewise.jl` |
| Unit Jacobian suite | Pass (`29/29`) | `test/test_jacobians_unit.jl` |
| GPU Jacobian suite | Pass (`4/4`) | `test/test_jacobians_GPU.jl` |
| `p₀` Jacobian exactness claim | Sensitive to FD setup; central FD can be near-exact for representative channels | `test/test_p0_convergence.jl` |
