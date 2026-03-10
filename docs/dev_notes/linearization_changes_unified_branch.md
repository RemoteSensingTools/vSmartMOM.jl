# CoreRT Linearization Changes: `sanghavi` → `unified-vsmartmom`

**Document Version:** 1.0  
**Date:** February 19, 2026  
**Authors:** vSmartMOM Development Team

---

## Overview

This document describes all changes to the linearized radiative transfer (RT) code in the `CoreRT` module between the `sanghavi` and `unified-vsmartmom` branches. The changes include critical bug fixes, architectural improvements, and code modernization.

The linearized RT computes both forward radiances and their analytic Jacobians with respect to physical parameters using the Matrix Operator Method (MOM) with the linearization approach of Sanghavi & Stephens (2013) and Sanghavi, Davis & Eldering (2014).

---

## Table of Contents

1. [Summary Statistics](#1-summary-statistics)
2. [Architecture & Type System Overhaul](#2-architecture--type-system-overhaul)
3. [Critical Bug Fixes](#3-critical-bug-fixes)
4. [New Doubling Implementation for N Parameters](#4-new-doubling-implementation-for-n-parameters)
5. [Elemental Layer Updates](#5-elemental-layer-updates)
6. [Interaction Function Updates](#6-interaction-function-updates)
7. [Layer Optical Properties](#7-layer-optical-properties)
8. [Code Organization & File Restructuring](#8-code-organization--file-restructuring)
9. [Syntax Modernization](#9-syntax-modernization)
10. [Theory Summary](#10-theory-summary)
11. [References](#11-references)

---

## 1. Summary Statistics

| Metric | Value |
|--------|-------|
| Files changed | 64 |
| Lines added | ~6,541 |
| Lines removed | ~9,269 |
| Net change | -2,728 lines |

### Key Files Modified

| File | Changes |
|------|---------|
| `types_lin.jl` | -496 lines (cleanup) |
| `rt_run_lin.jl` | -450 lines (cleanup) |
| `rt_kernel_lin.jl` | Major rewrite |
| `doubling_lin.jl` | +200 lines (new N-param doubling) |
| `lin_added_layer_all_params.jl` | Bug fixes |
| `elemental_lin.jl` | Simplified interface |
| `interaction_lin.jl` | Documentation + fixes |

---

## 2. Architecture & Type System Overhaul

### 2.1 Streamlined `vSmartMOM_Lin` Struct

**Before (sanghavi):**
```julia
mutable struct vSmartMOM_Lin
    τ̇_abs::AbstractArray{AbstractArray}
    τ̇_aer::AbstractArray{AbstractArray}
    lin_aerosol_optics::AbstractArray{AbstractArray{linAerosolOptics}}
end
```

**After (unified-vsmartmom):**
```julia
"""
Holds linearized (Jacobian) model parameters: derivatives of optical depths
and aerosol properties w.r.t. physical state-vector elements.
"""
mutable struct vSmartMOM_Lin{A,B,C}
    "∂τ_abs/∂x per band: Vector of arrays [NGas × nSpec × nLayers]"
    τ̇_abs::A
    "∂τ_aer/∂x per band: Vector of arrays [NAer × 7 × nSpec × nLayers]"
    τ̇_aer::B
    "Linearized aerosol optics per band per aerosol"
    lin_aerosol_optics::C
end
```

### 2.2 New `CoreScatteringOpticalPropertiesLin` Struct

This struct serves as the **AD boundary** — the clean interface between physical parameter derivatives and RT kernel propagation:

```julia
"""
Per-layer Jacobian of the four core optical properties (τ, ϖ, Z⁺⁺, Z⁻⁺)
with respect to the retrieval state vector x.

This struct is the **AD boundary**: everything upstream of it (Mie code,
absorption cross-sections, atmospheric profiles) may use ForwardDiff.Dual
numbers, but by the time values reach the RT kernels they are extracted
into plain Float64/Float32 arrays stored here.
"""
Base.@kwdef struct CoreScatteringOpticalPropertiesLin{T1,T2,T3} <: AbstractOpticalPropertiesLin
    "∂τ/∂x — [Nparams] or [Nparams × nSpec]"
    τ̇::T1
    "∂ϖ/∂x — [Nparams] or [Nparams × nSpec]"
    ϖ̇::T2
    "∂Z⁺⁺/∂x — [nμ × nμ × nSpec] or [Nparams × nμ × nμ × nSpec]"
    Ż⁺⁺::T3
    "∂Z⁻⁺/∂x — [nμ × nμ × nSpec] or [Nparams × nμ × nμ × nSpec]"
    Ż⁻⁺::T3
end
```

### 2.3 Removed Legacy Types

The following unused/commented-out types were removed:
- `RT_Aerosol_Lin`
- `ScatteringParametersLin`
- `vSmartMOM_Lin_Parameters`
- `ComputedAtmospherePropertiesLin`
- `CompositeLayerMS` (multisensor linearized)

---

## 3. Critical Bug Fixes

### 3.1 Bug 14: Elemental Optical Depth Scaling

**File:** `lin_added_layer_all_params.jl`

**Problem:** The chain rule was using full-layer τ derivatives instead of elemental τ derivatives.

**Root Cause:** Core derivatives from `elemental!` are computed w.r.t. the **elemental** optical depth `dτ = τ/2^n_d`, not the full layer τ.

**Fix:**
```julia
# Bug 14 fix: scale τ̇ by 2^ndoubl for elemental optical depth
dτ̇ = τ̇ ./ FT(2^ndoubl)

# Chain rule now uses dτ̇ instead of τ̇
@views ap_ṫ⁺⁺[iparam,:,:,:] .= added_layer_lin.ṫ⁺⁺[1,:,:,:].*reshape(dτ̇[iparam,:],1,1,nspec) .+ ...
```

### 3.2 Bug 17: Source Function Variable References

**File:** `interaction_lin.jl`

**Problem:** In `ScatteringInterface_11` (full scattering case), the source function interaction was using `ap_J̇₀⁻, ap_J̇₀⁺` from `added_layer_lin` instead of `J̇₀⁻, J̇₀⁺` from `composite_layer_lin`.

**Impact:** Double-counting of surface derivatives in the source function propagation.

**Fix:**
```julia
# Before (incorrect):
tmpap_J̇₀⁻[iparam,:,:,:] .= ap_J̇₀⁻[iparam,:,:,:] .+ ...
T01_inv ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ J₀⁺ .+ r⁻⁺ ⊠ ap_J̇₀⁺[iparam,:,:,:] .+ ...)

# After (correct):
tmpap_J̇₀⁻[iparam,:,:,:] .= J̇₀⁻[iparam,:,:,:] .+ ...
T01_inv ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ J₀⁺ .+ r⁻⁺ ⊠ J̇₀⁺[iparam,:,:,:] .+ ...)
```

### 3.3 Bug 18: Index Mapping in `createAero`

**File:** `compEffectiveLayerProperties_lin.jl`

**Problem:** Off-by-one error in aerosol parameter derivative indexing, causing τ_ref derivatives to be mixed with Mie parameter derivatives.

**Fix:**
```julia
# Before (incorrect):
τ̇_mod[1+iparam,:] = (1 .- fᵗ .* ω̃) .* τ̇Aer[1+iparam,:]

# After (correct):
τ̇_mod[1+iparam,:] = (1 .- fᵗ .* ω̃) .* τ̇Aer[iparam,:]
```

### 3.4 Bug 19: Z Chain Rule Must Be Applied Before Doubling (Critical)

**Files:** `rt_kernel_lin.jl`, `doubling_lin.jl`

**Problem:** The Z-matrix derivative `∂r/∂Z` is a 4th-rank tensor that is **diagonal** at the elemental level (each `r_{ij}` depends only on `Z_{ij}`). After doubling, matrix products mix all `(i,j)` indices, making the element-wise chain rule `ṙ[3] .* Ż` mathematically incorrect.

**Mathematical Explanation:**

At elemental level:
```
r_{ij} = f(Z_{ij})  →  ∂r_{ij}/∂Z_{kl} = δ_{ik}δ_{jl} × (∂r_{ij}/∂Z_{ij})
```

After one doubling step:
```
R_{ij} = r_{ij} + Σ_k t_{ik} × G_{kk'} × r_{k'j} × t_{ji}
```
Now `R_{ij}` depends on **all** `Z_{kl}` entries, not just `Z_{ij}`.

**Old Approach (incorrect):**
1. `elemental!` → 3-derivative arrays `(ṙ[1], ṙ[2], ṙ[3])`
2. `doubling!` → propagate 3-derivative arrays through doubling
3. `lin_added_layer_all_params!` → apply chain rule: `ap_ṙ = ṙ[1]×τ̇ + ṙ[2]×ϖ̇ + ṙ[3]×Ż`

**New Approach (correct):**
1. `elemental!` → 3-derivative arrays `(ṙ[1], ṙ[2], ṙ[3])`
2. `lin_added_layer_all_params!` → apply chain rule **at elemental level**
3. `doubling_allparams!` → propagate N-parameter arrays through doubling with proper matrix products

**Implementation:**
```julia
# rt_kernel_lin.jl - Bug 19 fix
@timeit "chain_rule" lin_added_layer_all_params!(
    RS_type, pol_type, SFI, quad_points, 
    computed_layer_properties_lin, 
    added_layer_lin, architecture, ndoubl)  # Apply BEFORE doubling

@timeit "doubling" doubling_allparams!(pol_type, SFI, expk,
    ndoubl, added_layer, added_layer_lin, 
    I_static, architecture, dτ̇_allparams, μ₀)  # New N-param doubling
```

### 3.5 Bug 22: Per-Parameter τ̇_sum Beam Attenuation

**File:** `rt_kernel_lin.jl`

**Problem:** The beam attenuation factor `exp(-τ_sum/μ₀)` derivative was only computed for parameter index 1, not all parameters.

**Fix:**
```julia
# Bug 22 fix: Add per-parameter τ̇_sum beam attenuation contribution to SFI
if SFI
    nparams_τ_sum = size(τ̇_sum, 1)
    for iparam = 1:nparams_τ_sum
        @views added_layer_lin.ap_J̇₀⁺[iparam,:,1,:] .+= 
            added_layer.j₀⁺[:,1,:] .* reshape(-τ̇_sum[iparam,:] ./ μ₀, 1, nspec_here)
        @views added_layer_lin.ap_J̇₀⁻[iparam,:,1,:] .+= 
            added_layer.j₀⁻[:,1,:] .* reshape(-τ̇_sum[iparam,:] ./ μ₀, 1, nspec_here)
    end
end
```

---

## 4. New Doubling Implementation for N Parameters

### 4.1 New Function: `doubling_allparams_helper!`

**File:** `doubling_lin.jl`

This function propagates **N physical-parameter** derivatives through the doubling method, unlike the original `doubling_helper!` which only propagated the 3 core derivatives (τ, ϖ, Z).

**Key Features:**
- Per-parameter beam attenuation derivatives for SFI
- Proper matrix product propagation for N parameters
- D-matrix transformation for polarization

**Signature:**
```julia
function doubling_allparams_helper!(pol_type, SFI, expk, ndoubl::Int, 
                                    added_layer::AddedLayer,
                                    added_layer_lin::AddedLayerLin,
                                    I_static::AbstractArray{FT}, 
                                    architecture, dτ̇::AbstractArray, μ₀::FT)
```

**Doubling Formulas (de Haan et al. 1987):**

For a homogeneous layer with reflection **R** and transmission **T**, the doubled layer has:

$$\mathbf{G} = (\mathbf{I} - \mathbf{R} \mathbf{R})^{-1}$$

$$\mathbf{R}_{2\tau} = \mathbf{R} + \mathbf{T} \, \mathbf{G} \, \mathbf{R} \, \mathbf{T}$$

$$\mathbf{T}_{2\tau} = \mathbf{T} \, \mathbf{G} \, \mathbf{T}$$

**Linearized Doubling:**

For each parameter $p_j$:

$$\dot{\mathbf{G}}_{p_j} = \mathbf{G} (\dot{\mathbf{R}}_{p_j} \mathbf{R} + \mathbf{R} \dot{\mathbf{R}}_{p_j}) \mathbf{G}$$

$$\dot{\mathbf{R}}_{2\tau,p_j} = \dot{\mathbf{R}}_{p_j} + \dot{\mathbf{T}}_{p_j} \mathbf{G} \mathbf{R} \mathbf{T} + \mathbf{T} \dot{\mathbf{G}}_{p_j} \mathbf{R} \mathbf{T} + \ldots$$

---

## 5. Elemental Layer Updates

### 5.1 Simplified Interface

**File:** `elemental_lin.jl`

The elemental function now uses the combined `computed_layer_properties` struct instead of individual arrays:

```julia
# Before:
function elemental!(pol_type, SFI, τ_sum, dτ_λ, dτ, ϖ_λ, ϖ, Z⁺⁺, Z⁻⁺,
                   τ̇_sum, dτ̇_λ, dτ̇, ϖ̇_λ, ϖ̇, Ż⁺⁺, Ż⁻⁺, F₀, m, ndoubl, ...)

# After:
function elemental!(pol_type, SFI, τ_sum, τ̇_sum, dτ, F₀,
                   computed_layer_properties, m, ndoubl, scatter, quad_points,
                   added_layer, added_layer_lin, architecture)
```

### 5.2 Safe Derivative Formulas

Guards against division by zero when ϖ=0 or Z=0:

```julia
# derivative wrt ϖ: safe guard against ϖ=0
ṙ⁻⁺[2, i, j, n] = ϖ_λ[n] == 0 ? FT(0) : r⁻⁺[i, j, n] / ϖ_λ[n]

# derivative wrt Z: direct formula avoids 0/0 when Z=0
ṙ⁻⁺[3,i,j,n] = ϖ_λ[n] * (qp_μN[j] / (qp_μN[i] + qp_μN[j])) * wct[j] * 
               (1 - exp(-dτ_λ[n] * ((1/qp_μN[i]) + (1/qp_μN[j]))))
```

### 5.3 Field Naming Convention

Source function fields use lowercase for elemental layer, uppercase for composite:
- `added_layer.j₀⁺, j₀⁻` (elemental)
- `composite_layer.J₀⁺, J₀⁻` (composite)

---

## 6. Interaction Function Updates

### 6.1 Enhanced Documentation

**File:** `interaction_lin.jl`

All `ScatteringInterface` dispatch functions now have comprehensive docstrings:

```julia
"""
    interaction_helper!(::ScatteringInterface_11, ...)

**Full scattering interaction (general case).** Both composite and added layers scatter.
This requires the full Adding Method with geometric series inversion:

```math
\\mathbf{G} = (\\mathbf{I} - \\mathbf{R}^{-+}_\\text{add} \\, \\mathbf{R}^{+-}_\\text{comp})^{-1}
```
"""
```

### 6.2 Type Signatures

Simplified from `FT<:Union{AbstractFloat, ForwardDiff.Dual}` to `FT<:Real`:

```julia
# Before:
function interaction_helper!(...) where {FT<:Union{AbstractFloat, ForwardDiff.Dual},FT2}

# After:
function interaction_helper!(...) where {FT<:Real, FT2}
```

### 6.3 Performance Annotations

Added `@inbounds` for performance-critical loops:

```julia
@inbounds for iparam=1:Nparams
    tmp_inv_lin[iparam,:,:,:] .= tmp_inv ⊠ (ap_ṙ⁻⁺[iparam,:,:,:] ⊠ R⁺⁻ .+ ...)
end
```

---

## 7. Layer Optical Properties

### 7.1 AD Architecture Boundary

**File:** `compEffectiveLayerProperties_lin.jl`

This file defines the interface between physical parameter differentiation and RT kernel propagation:

```julia
#=
AD Architecture Boundary
========================
Physical parameters (VMR, aerosol τ_ref/n_r/n_i/size/profile, surface albedo) produce
derivatives ∂(τ,ϖ,Z)/∂x_phys via Mie theory (ForwardDiff), Beer-Lambert law, and
mixing rules (quotient rule). These derivatives are computed HERE.

The RT kernel (elemental!, doubling!, interaction!) then only needs to propagate
w.r.t. the 3 core variables (τ, ϖ, Z) using the chain rule:
    ∂R/∂x = ∂R/∂τ · ∂τ/∂x + ∂R/∂ϖ · ∂ϖ/∂x + ∂R/∂Z · ∂Z/∂x
=#
```

### 7.2 δ-M Truncation Chain Rule

The `createAero` function computes δ-M scaled aerosol properties and their derivatives:

**δ-M Scaling Formulas:**

$$\tau_{\text{mod}} = (1 - f^t \tilde{\omega}) \cdot \tau_{\text{aer}}$$

$$\varpi_{\text{mod}} = \frac{(1 - f^t) \tilde{\omega}}{1 - f^t \tilde{\omega}}$$

**Derivative Chain Rule:**

$$\frac{\partial \tau_{\text{mod}}}{\partial p_j} = (1 - f^t \tilde{\omega}) \frac{\partial \tau_{\text{aer}}}{\partial p_j} - \tau_{\text{aer}} \left(f^t \frac{\partial \tilde{\omega}}{\partial p_j} + \tilde{\omega} \frac{\partial f^t}{\partial p_j}\right)$$

---

## 8. Code Organization & File Restructuring

### 8.1 File Relocations

| Old Location | New Location |
|--------------|--------------|
| `src/CoreRT/atmo_prof_lin.jl` | `src/CoreRT/tools/atmo_prof_lin.jl` |
| `src/CoreRT/lin_model_from_parameters.jl` | `src/CoreRT/tools/lin_model_from_parameters.jl` |
| `src/CoreRT/lambertian_surface_lin.jl` | `src/CoreRT/Surfaces/lambertian_surface_lin.jl` |
| `src/CoreRT/postprocessing_vza_lin.jl` | `src/CoreRT/tools/postprocessing_vza_lin.jl` |
| `src/CoreRT/rt_helper_functions_lin.jl` | `src/CoreRT/tools/rt_helper_functions_lin.jl` |
| `src/CoreRT/gpu_batched.jl` | `ext/gpu_batched_cuda.jl` (extension) |

### 8.2 New Files

- `src/CoreRT/CoreKernel/rt_helpers.jl` — Common RT helper functions
- `src/CoreRT/tools/cpu_batched.jl` — CPU batched operations
- `src/CoreRT/Surfaces/canopy_surface.jl` — Canopy RT surface model

### 8.3 Removed Files (Merged/Deleted)

- `src/CoreRT/model_from_parameters.jl` — Merged into tools/
- `src/CoreRT/parameters_from_yaml.jl` — Moved to main module
- `src/CoreRT/postprocessing_vza_ms.jl` — Consolidated

---

## 9. Syntax Modernization

### 9.1 Named Tuple Destructuring

Julia 1.7+ style:

```julia
# Before:
@unpack obs_alt, sza, vza, vaz = model.obs_geom
@unpack τ, ϖ, Z⁺⁺, Z⁻⁺ = computed_layer_properties

# After:
(; obs_alt, sza, vza, vaz) = model.obs_geom
(; τ, ϖ, Z⁺⁺, Z⁻⁺) = computed_layer_properties
```

### 9.2 Array Type Handling

```julia
# Before:
Array(...)

# After:
collect(...)      # For CPU
arr_type(...)     # For GPU compatibility
```

### 9.3 Type-Stable Zeros

```julia
# Before:
r⁻⁺ .= 0.0

# After:
r⁻⁺ .= zero(FT)
```

### 9.4 Helper Functions

New helper functions for common operations:

```julia
fourier_weight(m, FT) = m == 0 ? FT(0.50) : FT(0.25)
scaled_weights(m, wt_μN) = m == 0 ? wt_μN/2 : wt_μN/4
```

---

## 10. Theory Summary

### 10.1 Linearized RT Equation

The linearized RT computes both forward radiances and Jacobians via:

$$\frac{\partial \mathbf{R}}{\partial p_j} = \frac{\partial \mathbf{R}}{\partial \tau} \frac{\partial \tau}{\partial p_j} + \frac{\partial \mathbf{R}}{\partial \varpi} \frac{\partial \varpi}{\partial p_j} + \frac{\partial \mathbf{R}}{\partial \mathbf{Z}} \frac{\partial \mathbf{Z}}{\partial p_j}$$

where $p_j$ is any physical parameter in the state vector.

### 10.2 Key Architectural Change

The chain rule is now applied **at the elemental level** (before doubling), which is mathematically correct because the Z-derivative tensor is diagonal only at the single-scattering level.

### 10.3 Parameter Layout in Jacobians

The `Nparams` derivative dimension in `dR` is ordered as:

1. **Gas parameters** (NGas): `[q, VMR_1, VMR_2, ...]`
2. **Aerosol sub-parameters** (7 per aerosol): `[τ_ref, nᵣ, nᵢ, rₘ, σᵣ, p₀, σp]` for each aerosol
3. **Surface parameters** (NSurf): `[albedo_band1, albedo_band2, ...]`

---

## 11. References

1. **Sanghavi, S. & Stephens, G.** (2013). "Adaptation of the delta-M and Successive Order of Scattering methods to the matrix operator method." *JQSRT*, 117, 1–12.

2. **Sanghavi, S., Davis, A. & Eldering, A.** (2014). "vSmartMOM: A linearized discrete ordinate radiative transfer implementation." *JQSRT*, 146, 182–207.

3. **de Haan, J.F., Bosma, P.B. & Hovenier, J.W.** (1987). "The adding method for multiple scattering calculations of polarized light." *A&A*, 183, 371–391.

4. **Nakajima, T. & Tanaka, M.** (1988). "Algorithms for radiative intensity calculations in moderately thick atmospheres using a truncation approximation." *JQSRT*, 40, 51–69.

---

## Appendix A: Bug Summary Table

| Bug # | File(s) | Description | Impact |
|-------|---------|-------------|--------|
| 14 | `lin_added_layer_all_params.jl` | Elemental τ scaling | ~10% error in τ derivatives |
| 17 | `interaction_lin.jl` | Source function variable refs | Surface derivative double-counting |
| 18 | `compEffectiveLayerProperties_lin.jl` | Index mapping in createAero | τ_ref/Mie derivative mixing |
| 19 | `rt_kernel_lin.jl`, `doubling_lin.jl` | Z chain rule timing | Z derivatives incorrect after doubling |
| 22 | `rt_kernel_lin.jl` | Per-param τ̇_sum beam attenuation | Missing beam attenuation derivatives |

---

## Appendix B: Testing

### B.1 Validation Test File

See `test/test_forward_lin.jl` for comprehensive FD vs analytic Jacobian validation at two levels:

1. **Model-level**: Tests `τ̇_aer`, `k̇`, `ω̃̇` from `model_from_parameters`
2. **RT-level**: Tests `dR` (TOA reflectance Jacobians) from `rt_run`

### B.2 Expected Accuracy

| Parameter Type | Expected Mean Error |
|----------------|---------------------|
| Surface albedo | < 1% |
| τ_ref | < 15% (Bug 19 residual) |
| Mie params (nᵣ, nᵢ, μ, σ) | < 20% |
| Profile params (p₀, σp) | < 15% |

---

*Document generated from git diff analysis between `sanghavi` and `unified-vsmartmom` branches.*
