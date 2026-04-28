# Inelastic Scattering

**For:** users choosing Raman/Cabannes modes and developers extending inelastic scattering.

**Next:** [Core RT Theory](../vSmartMOM/CoreRTTheory.md), [Add a Raman Mode](../extending/raman.md), [API Reference](../api_reference.md).

The `InelasticScattering` module defines the Raman mode types and molecular precomputations used by the CoreRT solver. The forward model can run without Raman redistribution, with rotational Raman scattering, or with vibrational Raman scattering modes.

This page documents the current public mode layer. Solar-induced fluorescence is intentionally not documented here until the product/data policy for shipped fixtures is settled.

## Mode Types

| Type | Meaning | Typical use |
| --- | --- | --- |
| [`noRS`](@ref) | Pure elastic scattering mode with no Raman redistribution | Default mode for ordinary `rt_run(model)` calls |
| [`RRS`](@ref) | Rotational Raman scattering for an N2/O2 atmosphere | Filling-in calculations over one spectral grid |
| [`VS_0to1`](@ref) | Vibrational Stokes-style mode | Single vibrational transition family |
| [`VS_1to0`](@ref) | Vibrational anti-Stokes-style mode | Single vibrational transition family |
| [`noRS_plus`](@ref) | Elastic mode for concatenated-grid code paths | Plus-mode workflows that share Raman interfaces |
| [`RRS_plus`](@ref) | Rotational Raman for concatenated spectral bands | Multi-band Raman calculations |
| [`VS_0to1_plus`](@ref) | Vibrational Stokes-style plus mode | Concatenated-grid vibrational workflows |
| [`VS_1to0_plus`](@ref) | Vibrational anti-Stokes-style plus mode | Concatenated-grid vibrational workflows |

The mode type is the dispatch key. CoreRT asks trait methods such as `has_inelastic`, `uses_cabannes_phase`, and `needs_interaction_workspace` what workspaces and optical-property transforms are needed for the selected mode.

## What Gets Precomputed

Raman modes carry:

- molecular constants for the active species;
- Raman and Cabannes/Rayleigh optical-property weights;
- wavelength-redistribution indices;
- Raman phase-function Greek coefficients;
- solar or stellar irradiation arrays used by the inelastic source terms.

Molecular constants are stored as `Float64` even when the RT work arrays use `Float32`. Raman cross-section prefactors can sit near the `Float32` subnormal range, so this mixed precision is intentional: molecular state remains stable while the solver arrays keep their configured element type.

## Current Solver Integration

The inelastic path is not a separate solver. It is an alternate dispatch path through the same adding-doubling machinery:

1. the selected Raman mode builds or normalizes the required optical properties;
2. CoreRT allocates only the workspaces required by that mode;
3. the elemental, doubling, and interaction kernels propagate the source terms through the layer stack.

For the current code map, see the inelastic section of [Core RT Theory](../vSmartMOM/CoreRTTheory.md).

## Useful APIs

- [`AbstractRamanType`](@ref)
- [`noRS`](@ref), [`RRS`](@ref), [`VS_0to1`](@ref), [`VS_1to0`](@ref)
- [`noRS_plus`](@ref), [`RRS_plus`](@ref), [`VS_0to1_plus`](@ref), [`VS_1to0_plus`](@ref)
- [`compute_γ_air_Cabannes!`](@ref), [`compute_γ_air_Rayleigh!`](@ref), [`apply_lineshape!`](@ref)
