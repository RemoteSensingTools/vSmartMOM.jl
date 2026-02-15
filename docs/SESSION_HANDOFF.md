# Session Handoff: Analytic Jacobian Debugging for vSmartMOM.jl

## Repository & Branch
- **Repo:** `https://github.com/RemoteSensingTools/vSmartMOM.jl`
- **Branch:** `unified-vsmartmom`
- **Latest commit:** `f269ddc` — "Fix Bugs 1-21 in analytic Jacobians + comprehensive test suite"

## Project Context

vSmartMOM.jl is a Julia radiative transfer (RT) code using the Matrix Operator Method. The `unified-vsmartmom` branch merges features from the `TOMAS-aerosols` and `sanghavi` branches, specifically integrating **analytic Jacobians** of the forward model (from the sanghavi branch) to replace slow automatic differentiation.

The analytic Jacobians compute derivatives of TOA radiance with respect to physical parameters:
- **7 aerosol sub-params** per aerosol: `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]`
- **Gas VMR** (e.g., O₂)
- **Surface albedo**

The linearization follows Sanghavi, Davis & Eldering (2014) — derivatives propagate through:
1. **Elemental layer** (single-scattering, Sanghavi & Stephens 2013 Eqs. 19-34)
2. **Doubling** (de Haan et al. 1987) with 3 core derivatives (τ, ϖ, Z)
3. **Chain rule** mapping 3 core → N physical parameters
4. **Adding/interaction** combining layers from TOA down

## What Has Been Done

### 21 Bugs Fixed
All documented in `docs/LINEARIZATION_BUGS.md` with root causes and fixes.

**Categories:**
- **Bugs 1-17:** Compilation/runtime (loop ranges, dimensions, types, division by zero, naming)
- **Bug 18:** Index mapping in `createAero` — off-by-one mixing τ_ref with Mie param derivatives
- **Bug 19:** Z chain rule structural issue — **IDENTIFIED but fix is UNTESTED** (see below)
- **Bug 20a/b:** Mie derivative interpolation in `lin_model_from_parameters.jl:498-515`:
  - `k̇sca` used `k̇ × ω̃̇` (product of derivatives) instead of product rule `k̇×ω̃ + k×ω̃̇`
  - `ω̃̇` recovered via `k̇sca/k̇` (divides by derivative!) instead of quotient rule `(k̇sca - ω̃×k̇)/k`
- **Bug 21:** `τ̇_aer` for microphysical params missing `k_band` factor in `k_ref` normalization

### Docstrings
LaTeX-equation docstrings added to all key linearized RT functions referencing the relevant papers.

### Test Results (4-wavelength fast config)

| Parameter | Pathway | Max Rel Err | Mean Rel Err | Status |
|-----------|---------|-------------|--------------|--------|
| Surface albedo | Boundary only | 6e-6 | 3e-6 | ✅ Perfect |
| τ_ref | createAero → RT kernel | ~12% | ~5.8% | ⚠️ Residual |
| p₀ | Profile → createAero → RT | peaks ~10% | inflated by zeros | ⚠️ OK-ish |
| nᵣ | **Mie → interp → createAero → RT** | ~22% | ~10% | ⚠️ Residual |

### Key Diagnostic Insight
The error pattern reveals which code paths are problematic:
- **Surface albedo is perfect** → adding/interaction code is correct
- **τ_ref ~10% error** → RT kernel linearization mostly works; residual is from Z chain rule
- **nᵣ ~10% error after Mie fixes** → Mie interpolation bugs (20/21) are fixed; same residual
- **All parameters with nonzero Ż show ~10% error** → Bug 19 (Z chain rule after doubling) is the remaining issue

## What Needs To Be Done Next

### Priority 1: Run Tests on Server
```bash
git clone https://github.com/RemoteSensingTools/vSmartMOM.jl
cd vSmartMOM.jl
git checkout unified-vsmartmom
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. test/test_jacobians_unit.jl    # ~5 min, @testset-based
julia --project=. test/test_mie_jacobians.jl     # ~15 min, FD convergence for nᵣ, nᵢ, μ, σ
```

### Priority 2: Fix Bug 19 — Z Chain Rule After Doubling
This is the **root cause of the remaining ~10% error** for all aerosol parameters.

**The problem:** After the elemental step, the 3-core derivative `ṙ⁻⁺[3,:,:,:]` represents `∂r_ij/∂Z_ij` — an element-wise relationship. The chain rule (`lin_added_layer_all_params!`) uses element-wise multiplication `ṙ[3] .* Ż` to map to physical params. But currently this chain rule is applied **after doubling**, where matrix products have mixed indices and the relationship is no longer element-wise (it's a 4th-rank tensor `∂R_ij/∂Z_kl`).

**The fix strategy** (prototype code exists in `doubling_lin.jl` as `doubling_allparams_helper!`):
1. Apply chain rule **before** doubling (at elemental level where element-wise is correct)
2. Propagate N physical-parameter derivatives through a new `doubling_allparams!` function
3. Skip the 3-core doubling (no longer needed)

**Flow change:**
```
CURRENT (buggy):  elemental → doubling(3-core) → chain_rule → interaction
NEW (correct):    elemental → chain_rule → doubling(N-params) → interaction
```

**Key files:**
- `src/CoreRT/CoreKernel/rt_kernel_lin.jl` — orchestrates the flow per layer
- `src/CoreRT/CoreKernel/doubling_lin.jl` — contains both old `doubling_helper!` and new `doubling_allparams_helper!`
- `src/CoreRT/CoreKernel/lin_added_layer_all_params.jl` — the chain rule function
- `src/CoreRT/CoreKernel/elemental_lin.jl` — computes 3-core derivatives

**Important subtleties:**
- The D matrix (`Diagonal{1,1,-1,-1}` for Stokes) transforms r⁻⁺→r⁺⁻ and t⁺⁺→t⁻⁻. After elemental, D is already applied. After doubling, D must be reapplied. The new code handles this.
- For SFI (Source Function Integration), per-parameter beam attenuation derivatives `d(exp(-τ/μ₀))/dp_j` must be computed and propagated through doubling.
- The TOA layer (iz==1) currently applies chain rule inline; it should use the same approach.

### Priority 3: Test Gas VMR Jacobians
Not yet tested. These go through `τ̇_abs` (absorption optical depth derivatives) and should be relatively simple since gas doesn't scatter (Ż=0 for gas-only parameters).

### Priority 4: CUDA/GPU Testing
The code uses `KernelAbstractions.jl` for GPU portability. Key GPU files:
- `src/CoreRT/CoreKernel/elemental_lin.jl` — KA kernels for elemental step
- The `⊠` operator (batched matrix multiply) is the main GPU-accelerated operation

## File Map

### Core Linearized RT (where the bugs are)
```
src/CoreRT/CoreKernel/
  elemental_lin.jl          — Single-scattering R,T,J and 3-core derivatives
  doubling_lin.jl           — Doubling method for R,T,J + derivatives
  lin_added_layer_all_params.jl — Chain rule: 3 core → N physical params
  rt_kernel_lin.jl          — Orchestrates elemental → doubling → chain → interaction per layer
  interaction_lin.jl        — Adding method: combine current layer with composite above
```

### Optical Property Construction
```
src/CoreRT/LayerOpticalProperties/
  compEffectiveLayerProperties_lin.jl — createAero (δ-M scaling), constructCoreOpticalProperties
src/CoreRT/tools/
  lin_model_from_parameters.jl       — Mie derivatives, interpolation, τ̇_aer (Bugs 20, 21)
```

### Mie Scattering
```
src/Scattering/
  compute_NAI2_lin.jl    — Linearized Mie calculation (ForwardDiff for an, bn)
  mie_helper_functions_lin.jl — Size distribution weight derivatives (ẇₓ)
```

### Types
```
src/CoreRT/types_lin.jl — AddedLayerLin, CompositeLayerLin, CoreScatteringOpticalPropertiesLin
                          Also: Base.:+ operators for combining optical properties
```

### Tests
```
test/
  test_jacobians_unit.jl     — @testset unit tests (surface, τ_ref, p₀, nᵣ, σ) — RUN THIS FIRST
  test_mie_jacobians.jl      — FD convergence for microphysical params
  test_fd_convergence.jl     — FD convergence for τ_ref, p₀, σp, albedo
  test_jacobian_diagnostic.jl — Multi-level diagnostic
  test_parameters/
    JacobianTestFast.yaml    — 4 wavelengths, fast (~5 min total)
    JacobianTest.yaml        — 169 wavelengths, full
```

### Documentation
```
docs/
  LINEARIZATION_BUGS.md  — All 21 bugs with detailed root cause analysis
  SESSION_HANDOFF.md     — This file
```

## Key Data Structures

```julia
# Forward RT matrices per layer
AddedLayer: r⁻⁺, r⁺⁻, t⁺⁺, t⁻⁻, j₀⁺, j₀⁻  # [nμ × nμ × nSpec]

# Linearized: 3-core + N all-params
AddedLayerLin:
  ṙ⁻⁺, ṫ⁺⁺, ...        # [3 × nμ × nμ × nSpec] — core (τ, ϖ, Z)
  ap_ṙ⁻⁺, ap_ṫ⁺⁺, ...  # [Nparams × nμ × nμ × nSpec] — physical params

# Per-layer optical properties and their derivatives
CoreScatteringOpticalProperties: τ, ϖ, Z⁺⁺, Z⁻⁺
CoreScatteringOpticalPropertiesLin: τ̇, ϖ̇, Ż⁺⁺, Ż⁻⁺  # [Nparams × ...]
```

## Parameter Indexing in dR Array
`dR[param_idx, vza, stokes, wavelength]` where:
- `1` = τ_ref
- `2` = nᵣ (real refractive index)
- `3` = nᵢ (imaginary refractive index)
- `4` = μ (log-mean of LogNormal size distribution)
- `5` = σ (log-std of LogNormal size distribution)
- `6` = p₀ (pressure center of aerosol profile)
- `7` = σp (pressure width of aerosol profile)
- `8..8+NGas-1` = gas VMR parameters
- `last` = surface albedo

## Papers Referenced
- Sanghavi & Stephens (2013): Elemental RT formulas (Eqs. 19-34)
- Sanghavi, Davis & Eldering (2014): Full linearization framework
- Sanghavi & Stephens (2015): δ-M truncation scaling
- de Haan, Bosma & Hovenier (1987): Doubling method (Eq. 25)
